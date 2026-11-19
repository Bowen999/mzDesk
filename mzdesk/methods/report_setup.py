"""RepoRT 准备，一次性运行：python -m mzdesk.methods.report_setup

RepoRT（michaelwitting/RepoRT，数据许可 CC BY-SA 4.0）是小分子保留时间库，Known methods 页用它查"同一化合物 / 结构最相近的化合物用过什么色谱条件"。
这里用固定提交号做稀疏克隆（git --filter=blob:none，只取每个数据集的 rtdata_canonical_success / info / metadata / gradient 四个 tsv，约 100 MB）；没有 git 就下载该提交的 GitHub 压缩包（约 190 MB，只留同样的 tsv），
再建一个紧凑的本地索引（report_index.npz + report_datasets.json，几十 MB），并把 LICENSE 与引用说明放在索引旁边。
文件都放在 ~/.mzdesk/data/methods/report/（环境变量 MZDESK_REPORT_DIR 或 MZDESK_METHODS_DIR 可改目录）。数据不入库。
电荷类别用 MolGpKa（权重已安装时，多进程，十分钟以内）；没装权重就用官能团规则，装好后用 --force 重建。
"""
import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
import time
from datetime import date
from pathlib import Path

from mzdesk.local import archive, paths
from mzdesk.methods import compound_index as ci
from mzdesk.methods import report

COMMIT = "1115f14202b3ff29a0fb3dd0d410a03969fe0d89"      # michaelwitting/RepoRT master, 2026-10-21 ("missing metadata yaml for 0392")
REPO_GIT = report.REPO_URL + ".git"
SPARSE = ["/LICENSE", "/README.md",
          "/processed_data/*/*_rtdata_canonical_success.tsv", "/processed_data/*/*_info.tsv",
          "/processed_data/*/*_metadata.tsv", "/processed_data/*/*_gradient.tsv"]
SETUP_COMMAND = report.SETUP_COMMAND
MIN_DATASETS = 300                                       # sanity check on the checkout (the pinned commit has 421)


# ---------------------------------------------------------------- download
def raw_dir():
    return report.report_dir() / "repo"


def _git(repo, *args, capture=False):
    """Run git in ``repo``; ``capture`` returns stdout and keeps stderr quiet (for queries that may fail)."""
    env = dict(os.environ, GIT_TERMINAL_PROMPT="0")
    cmd = [paths.working_git() or "git", "-C", str(repo), *args]
    return subprocess.run(cmd, check=True, env=env, text=True, stdout=subprocess.PIPE if capture else None,
                          stderr=subprocess.DEVNULL if capture else None).stdout


def fetch(repo, commit=COMMIT):
    """Sparse, blob-less checkout of ``commit`` in ``repo`` (about 100 MB); skipped when that commit is already checked out.
    Without a usable git: the commit's archive from codeload.github.com instead (about 190 MB to download, the same tsv files are kept)."""
    if paths.working_git() is None:
        return fetch_archive(repo, commit)
    repo.mkdir(parents=True, exist_ok=True)
    if not (repo / ".git").exists():
        _git(repo, "init", "-q")
    remotes = _git(repo, "remote", capture=True).split()
    _git(repo, "remote", "set-url" if "origin" in remotes else "add", "origin", REPO_GIT)
    _git(repo, "config", "core.sparseCheckout", "true")
    _git(repo, "config", "core.sparseCheckoutCone", "false")
    (repo / ".git" / "info").mkdir(exist_ok=True)
    (repo / ".git" / "info" / "sparse-checkout").write_text("\n".join(SPARSE) + "\n")
    try:
        head = _git(repo, "rev-parse", "HEAD", capture=True).strip()
    except subprocess.CalledProcessError:
        head = ""
    if head == commit:
        _git(repo, "-c", "advice.detachedHead=false", "checkout", "-q", "--force", commit)      # re-apply the sparse patterns
        print(f"RepoRT {commit[:10]} is already downloaded in {repo}")
        return
    print(f"Downloading RepoRT {commit[:10]} (sparse clone, only the needed tsv files) into {repo} ...", flush=True)
    try:
        _git(repo, "fetch", "--depth", "1", "--filter=blob:none", "origin", commit)
        _git(repo, "-c", "advice.detachedHead=false", "checkout", "-q", "--force", "FETCH_HEAD")
    except subprocess.CalledProcessError as exc:
        raise SystemExit(f"git failed ({exc}). Check the network connection and run again.") from None
    head = _git(repo, "rev-parse", "HEAD", capture=True).strip()
    if head != commit:
        raise SystemExit(f"Checked out {head}, expected {commit}.")


def fetch_archive(repo, commit=COMMIT):
    """No git: stream the commit's tar.gz and keep only the files SPARSE names."""
    if archive.has_commit(repo, commit):
        print(f"RepoRT {commit[:10]} is already downloaded in {repo}")
        return
    print(f"Downloading RepoRT {commit[:10]} (archive, about 190 MB; only the needed tsv files are kept) into {repo} ...", flush=True)
    try:
        n = archive.fetch(report.REPO_URL, commit, repo, SPARSE, log=lambda text: print(text, flush=True), timeout=120)
    except archive.ArchiveError as exc:
        raise SystemExit(str(exc)) from None
    print(f"{n} files kept")


def folder_size(path):
    return sum(f.stat().st_size for f in Path(path).rglob("*") if f.is_file())


# ---------------------------------------------------------------- parse the downloaded files
def _rows(path):
    """Rows of a RepoRT tsv as dicts (no quote handling: names contain quote characters)."""
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        yield from csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE)


def _num(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def _eluents(meta):
    out = {}
    for letter in "ABCD":
        solvents = {k: v for k, _ in report.SOLVENTS if (v := _num(meta.get(f"eluent.{letter}.{k}"))) and v > 0}
        additives = []
        for k, _ in report.ADDITIVES:
            v = _num(meta.get(f"eluent.{letter}.{k}"))
            if v and v > 0:
                additives.append([k, v, (meta.get(f"eluent.{letter}.{k}.unit") or "").strip()])
        ph = _num(meta.get(f"eluent.{letter}.pH"))
        if solvents:
            out[letter] = {"solvents": solvents, "additives": additives, "pH": ph if ph and ph > 0 else None}
    return out


def _gradient(path):
    t, comp = [], {L: [] for L in "ABCD"}
    if not path.is_file():
        return {}
    for r in _rows(path):
        row = {k.split(" ")[0]: _num(v) for k, v in r.items() if k}
        if row.get("t") is None:
            continue
        t.append(row["t"])
        for L in "ABCD":
            comp[L].append(row.get(L) if row.get(L) is not None else 0.0)
    if not t:
        return {}
    return {"t": t, **{L: v for L, v in comp.items() if any(v)}}


def parse_dataset(ddir):
    """processed_data/<id>/ -> (dataset dict, rt rows [(smiles, rt, name)]); None when the info or retention times are missing."""
    ddir = Path(ddir)
    ds_id = ddir.name
    info_path, rt_path, meta_path = (ddir / f"{ds_id}_{x}.tsv" for x in ("info", "rtdata_canonical_success", "metadata"))
    if not (info_path.is_file() and rt_path.is_file()):
        return None
    info = next(iter(_rows(info_path)), {}) or {}
    meta = next(iter(_rows(meta_path)), {}) if meta_path.is_file() else {}
    column = {"name": (meta.get("column.name") or "").strip(), "usp": (meta.get("column.usp.code") or "").strip(),
              "length": _num(meta.get("column.length")), "id": _num(meta.get("column.id")), "particle": _num(meta.get("column.particle.size")),
              "temp": _num(meta.get("column.temperature")), "flow": _num(meta.get("column.flowrate"))}
    ds = {"id": ds_id, "name": info.get("name", ""), "type": info.get("method.type", ""), "source": info.get("source", ""),
          "url": info.get("url", ""), "pmid": info.get("pmid", ""), "authors": info.get("authors", ""), "column": column,
          "eluents": _eluents(meta), "gradient": _gradient(ddir / f"{ds_id}_gradient.tsv"),
          "start": {L: v for L in "ABCD" if (v := _num(meta.get(f"gradient.start.{L}"))) is not None},
          "end": {L: v for L in "ABCD" if (v := _num(meta.get(f"gradient.end.{L}"))) is not None}}
    rows = [((r.get("smiles.std") or "").strip(), _num(r.get("rt")), (r.get("name") or "").strip()) for r in _rows(rt_path)]
    rows = [r for r in rows if r[0]]
    if not rows:
        return None
    ds["n_rows"] = len(rows)
    return ds, rows


# ---------------------------------------------------------------- build the index
def build_index(raw_root, out_dir, *, commit=COMMIT, method=None, workers=1, log=print, cfg=ci.CONFIG):
    """processed_data folder -> report_index.npz + report_datasets.json in ``out_dir``. Returns a stats dict.

    method: ci.METHOD_MOLGPKA / ci.METHOD_FG (default: MolGpKa when its weights are installed).
    """
    raw_root, out_dir = Path(raw_root), Path(out_dir)
    method = method or ci.best_charge_method()
    if method == ci.METHOD_MOLGPKA and not ci.molgpka_available()[0]:
        raise SystemExit(f"MolGpKa weights are not installed: {ci.molgpka_available()[1]}")
    started = time.time()
    datasets, rows = [], []
    for ddir in sorted(p for p in (raw_root / "processed_data").iterdir() if p.is_dir()):
        parsed = parse_dataset(ddir)
        if parsed is None:
            continue
        ds, rts = parsed
        ref = len(datasets)
        datasets.append(ds)
        rows.extend((smiles, ref, rt if rt is not None else float("nan"), name) for smiles, rt, name in rts)
    log(f"Read {len(datasets)} datasets, {len(rows)} retention-time entries ({time.time() - started:.0f} s)")
    started_idx = time.time()

    last = [0]

    def progress(done, total):
        if done - last[0] >= max(1, total // 10) or done == total:          # about ten lines
            last[0] = done
            log(f"  structures {done}/{total} ({time.time() - started_idx:.0f} s)")
    meta = {"source": "RepoRT", "repo": report.REPO_URL, "commit": commit, "licence": report.LICENCE, "built": date.today().isoformat(),
            "n_datasets": len(datasets)}
    index, skipped, fallbacks = ci.CompoundIndex.build(rows, meta, method, cfg, workers=workers, progress=progress)
    index.meta["skipped_entries"] = skipped
    out_dir.mkdir(parents=True, exist_ok=True)
    index.save(out_dir / report.INDEX_FILE)
    tmp = out_dir / (report.DATASETS_FILE + ".part")
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump({"commit": commit, "datasets": datasets}, f, ensure_ascii=False, separators=(",", ":"))
    os.replace(tmp, out_dir / report.DATASETS_FILE)
    return {"datasets": len(datasets), "entries": len(rows), "compounds": len(index), "skipped_entries": skipped,
            "charge_method": method, "charge_fallbacks": fallbacks, "index_seconds": time.time() - started_idx,
            "seconds": time.time() - started, "index_bytes": (out_dir / report.INDEX_FILE).stat().st_size + (out_dir / report.DATASETS_FILE).stat().st_size}


def write_attribution(out_dir, raw_root, stats, commit=COMMIT):
    """Licence and attribution next to the index."""
    out_dir = Path(out_dir)
    licence = Path(raw_root) / "LICENSE"
    if licence.is_file():
        shutil.copyfile(licence, out_dir / "LICENSE")
    (out_dir / "ATTRIBUTION.md").write_text(
        "# RepoRT data used by mzDesk (Known methods page)\n\n"
        f"* Source: {report.REPO_URL}, commit {commit}\n"
        f"* Licence: {report.LICENCE} ({report.LICENCE_URL}); the full licence text is in LICENSE next to this file.\n"
        f"* Citation: {report.CITATION}\n"
        "* Each dataset keeps its own origin (source, DOI) in its `info.tsv`; mzDesk shows it in the Reference column.\n"
        f"* Downloaded files: processed_data/*/*_rtdata_canonical_success.tsv, *_info.tsv, *_metadata.tsv, *_gradient.tsv ({stats['datasets']} datasets).\n"
        "* report_index.npz and report_datasets.json are a derived work (normalised structures, fingerprints, descriptors, "
        "dataset metadata); they are shared under the same licence (CC BY-SA 4.0: attribution, share alike).\n"
        f"* Built {date.today().isoformat()} with the {stats['charge_method']} charge classes.\n", encoding="utf-8")


def setup(force=False, workers=None, method=None, clean=False):
    out = report.report_dir()
    if report.installed(out) and not force:
        print(f"RepoRT index already exists in {out} (use --force to rebuild).")
        return out
    started = time.time()
    repo = raw_dir()
    fetch(repo)
    git_bytes = folder_size(repo / ".git")
    print(f"Download / checkout took {time.time() - started:.0f} s; {repo} holds {folder_size(repo) / 1e6:.0f} MB "
          f"({git_bytes / 1e6:.0f} MB of git data, {(folder_size(repo) - git_bytes) / 1e6:.0f} MB of tsv files)")
    n_datasets = sum(1 for p in (repo / "processed_data").iterdir() if p.is_dir()) if (repo / "processed_data").is_dir() else 0
    if n_datasets < MIN_DATASETS:
        raise SystemExit(f"Only {n_datasets} datasets found in {repo}/processed_data (expected about 421). Delete {repo} and run again.")
    workers = workers or max(1, (os.cpu_count() or 2) - 1)
    stats = build_index(repo, out, workers=workers, method=method)
    write_attribution(out, repo, stats)
    print(f"Index: {stats['compounds']} compounds, {stats['entries']} entries from {stats['datasets']} datasets "
          f"({stats['index_bytes'] / 1e6:.0f} MB) in {stats['seconds']:.0f} s with {workers} workers; charge classes: {stats['charge_method']}"
          + (f" ({stats['charge_fallbacks']} structures fell back to functional groups)" if stats["charge_fallbacks"] else "")
          + f"; {stats['skipped_entries']} entries skipped (structure not usable).")
    started_jp = time.time()
    from mzdesk.methods import jp18
    jp18.load()                                          # builds / refreshes the JP18 index cache too, so the first page run is fast
    print(f"JP18 index ready ({time.time() - started_jp:.0f} s).")
    if stats["charge_method"] != ci.METHOD_MOLGPKA:
        print("MolGpKa weights are not installed, so the functional-group fallback was used for the charge classes. "
              "Run `python -m mzdesk.platform_fit.pka_setup` and then this command with --force for pKa-based classes.")
    if clean:
        shutil.rmtree(repo, ignore_errors=True)
        print(f"Removed {repo}.")
    print(f"RepoRT is {report.LICENCE} (attribution in {out / 'ATTRIBUTION.md'}). Citation: {report.CITATION}")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="Download RepoRT (pinned commit, needed files only) and build the local index for the Known methods page.")
    ap.add_argument("--force", action="store_true", help="rebuild the index even if it exists (reuses the downloaded files)")
    ap.add_argument("--workers", type=int, default=None, help="processes for the structure descriptors (default: CPU count - 1)")
    ap.add_argument("--functional-groups", action="store_true", help="use the functional-group charge classes even if MolGpKa is installed")
    ap.add_argument("--clean", action="store_true", help="delete the downloaded tsv files after the index is built")
    args = ap.parse_args(argv)
    setup(args.force, args.workers, ci.METHOD_FG if args.functional_groups else None, args.clean)


if __name__ == "__main__":
    sys.exit(main())
