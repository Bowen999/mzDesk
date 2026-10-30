"""CCS 数据准备，一次性运行：python -m mzdesk.orthogonal.ccs_setup

CCSbase 的源数据（dylanhross/c3sdb，MIT 许可，Ross, Cho, Xu. Anal. Chem. 2020, 92, 4548）：26 个已发表的数据集，全部在 N2 中测得，
DT / TW / TIMS 三种平台。CCS 页用它查"同一化合物、同一加合物"的参考 CCS。
这里用固定提交号做稀疏克隆（git --filter=blob:none，只取 LICENSE、README、26 个 json 和 SMILES 缓存，约 5 MB），
读 json（只读数据，克隆里的任何代码都不导入、不运行），用 RDKit 规范化结构（最大片段、中性化、InChIKey 前 14 位 = 连接键），
建一个紧凑的本地索引（ccs_index.json.gz），再训练一个 CCS 预测模型（没有参考值时的后备；joblib，本机生成）。
文件都放在 ~/.mzdesk/data/ccs/（环境变量 MZDESK_CCS_DIR 可改目录），LICENSE.txt 放在索引旁边。数据不入库。
没有 SMILES 的条目跳过（不联网查）。--force 重建（重用已下载的文件）。
"""
import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import time
from datetime import date
from pathlib import Path

import numpy as np

from mzdesk.methods import compound_index as ci
from mzdesk.orthogonal import ccs
from mzdesk.orthogonal.ccs import CONFIG, Row

COMMIT = "5f1393eec6f0acea9af9c592b911d5602a9236ea"      # dylanhross/c3sdb master, 2024-12-03
REPO_GIT = ccs.REPO_URL + ".git"
SPARSE = ["/LICENSE.txt", "/README.md", "/c3sdb/_include/smiles_search_cache.json", "/c3sdb/_include/src_data/*.json"]
MIN_DATASETS = 20                                        # 检出后的自检（固定提交有 26 个）
CCS_TYPES = ("DT", "TW", "TIMS")
SKIP_REASONS = ("bad_value", "no_structure", "duplicate", "bad_structure")


# ---------------------------------------------------------------- 下载
def raw_dir():
    return ccs.ccs_dir() / "repo"


def _git(repo, *args, capture=False):
    """在 repo 里运行 git；capture 返回 stdout 并静音 stderr（用于可能失败的查询）。"""
    env = dict(os.environ, GIT_TERMINAL_PROMPT="0")
    cmd = ["git", "-C", str(repo), *args]
    return subprocess.run(cmd, check=True, env=env, text=True, stdout=subprocess.PIPE if capture else None,
                          stderr=subprocess.DEVNULL if capture else None).stdout


def fetch(repo, commit=COMMIT):
    """稀疏、无 blob 的检出 commit；已经检出同一提交就跳过。"""
    if shutil.which("git") is None:
        raise SystemExit("git is needed to download the CCS data (sparse clone of github.com/dylanhross/c3sdb). Install git and run again.")
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
        _git(repo, "-c", "advice.detachedHead=false", "checkout", "-q", "--force", commit)
        print(f"c3sdb {commit[:10]} is already downloaded in {repo}")
        return
    print(f"Downloading c3sdb {commit[:10]} (sparse clone, only the data files) into {repo} ...", flush=True)
    try:
        _git(repo, "fetch", "--depth", "1", "--filter=blob:none", "origin", commit)
        _git(repo, "-c", "advice.detachedHead=false", "checkout", "-q", "--force", "FETCH_HEAD")
    except subprocess.CalledProcessError as exc:
        raise SystemExit(f"git failed ({exc}). Check the network connection and run again.") from None
    head = _git(repo, "rev-parse", "HEAD", capture=True).strip()
    if head != commit:
        raise SystemExit(f"Checked out {head}, expected {commit}.")


# ---------------------------------------------------------------- 读下载的文件（只读 json）
def _read_json(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise SystemExit(f"{path} is missing or is not a regular file.")
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def read_source(raw_root):
    """克隆目录 -> (数据集列表 [(metadata, rows)], SMILES 缓存 {名字: SMILES})。"""
    inc = Path(raw_root) / "c3sdb" / "_include"
    cache = _read_json(inc / "smiles_search_cache.json")
    files = sorted((inc / "src_data").glob("*.json"))
    datasets = []
    for f in files:
        d = _read_json(f)
        meta = d["metadata"]
        if meta.get("ccs_type") not in CCS_TYPES:
            raise SystemExit(f"{f.name}: unknown ccs_type {meta.get('ccs_type')!r}.")
        datasets.append((meta, d["data"]))
    return datasets, cache


def _positive(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) and x > 0 else None


def build_rows(datasets, cache, log=print):
    """数据集 -> (Row 列表, 数据集信息 {src_tag: …}, 统计)。没有 SMILES / 数值不对 / 重复 / 结构用不了的行跳过并计数。"""
    norms, rows, sources = {}, [], {}
    skipped = dict.fromkeys(SKIP_REASONS, 0)
    n_total = 0
    for meta, data in datasets:
        tag = meta["src_tag"]
        seen = set()
        sources[tag] = {"reference": meta.get("reference", ""), "url": meta.get("url", ""), "ccs_type": meta["ccs_type"],
                        "ccs_method": meta.get("ccs_method", ""), "n": 0}
        for r in data:
            n_total += 1
            name = str(r.get("name") or "").strip()
            adduct = ccs.map_adduct(r.get("adduct") or "")
            mz, value = _positive(r.get("mz")), _positive(r.get("ccs"))
            if not name or not adduct or mz is None or value is None:
                skipped["bad_value"] += 1
                continue
            smiles = (r.get("smi") or cache.get(name) or "").strip()
            if not smiles:
                skipped["no_structure"] += 1
                continue
            sig = (name, adduct, value)
            if sig in seen:                                    # c3sdb 同样按 名字+加合物+CCS+类型+数据集 去重
                skipped["duplicate"] += 1
                continue
            seen.add(sig)
            if smiles not in norms:
                try:
                    norms[smiles] = ci.normalize(smiles)
                except ValueError:
                    norms[smiles] = None
            norm = norms[smiles]
            if norm is None:
                skipped["bad_structure"] += 1
                continue
            rows.append(Row(key=norm.key14, name=name, smiles=norm.smiles, adduct=adduct, mz=mz, z=abs(ccs.adduct_charge(adduct)) or 1,
                            ccs=value, ccs_type=meta["ccs_type"], ccs_method=meta.get("ccs_method", ""), src_tag=tag))
            sources[tag]["n"] += 1
    by_type = {t: sum(1 for r in rows if r.ccs_type == t) for t in CCS_TYPES}
    compounds_by_type = {t: len({r.key for r in rows if r.ccs_type == t}) for t in CCS_TYPES}
    stats = {"rows_in_source": n_total, "rows_kept": len(rows), "rows_skipped": skipped, "rows_by_type": by_type,
             "compounds": len({r.key for r in rows}), "compounds_by_type": compounds_by_type, "datasets": len(datasets)}
    log(f"Read {len(datasets)} datasets, {n_total} rows; kept {len(rows)} ({', '.join(f'{t} {n}' for t, n in by_type.items())}), "
        f"{stats['compounds']} compounds; skipped " + ", ".join(f"{k.replace('_', ' ')} {v}" for k, v in skipped.items()))
    return rows, sources, stats


# ---------------------------------------------------------------- 写文件
def write_attribution(out_dir, raw_root, stats, commit=COMMIT):
    out_dir = Path(out_dir)
    licence = Path(raw_root) / "LICENSE.txt"
    if licence.is_file():
        shutil.copyfile(licence, out_dir / "LICENSE.txt")
    (out_dir / "ATTRIBUTION.md").write_text(
        "# CCSbase source data used by mzDesk (CCS page)\n\n"
        f"* Source: {ccs.REPO_URL}, commit {commit}\n"
        f"* Licence: {ccs.LICENCE} ({ccs.LICENCE_URL}); the full licence text is in LICENSE.txt next to this file.\n"
        f"* Citation: {ccs.CITATION}\n"
        "* Each dataset keeps its own origin (reference, URL); mzDesk shows it in the row's ? on the CCS page.\n"
        f"* Downloaded files: LICENSE.txt, README.md, c3sdb/_include/smiles_search_cache.json, c3sdb/_include/src_data/*.json ({stats['datasets']} datasets).\n"
        "* ccs_index.json.gz, ccs_model.joblib and ccs_domain.npz are derived from that data (normalised structures, a prediction model "
        "trained on it, fingerprints of its training compounds); they are shared under the same licence (MIT).\n"
        f"* Built {date.today().isoformat()}.\n", encoding="utf-8")


def write_meta(out_dir, meta):
    out_dir = Path(out_dir)
    tmp = out_dir / (ccs.META_FILE + ".part")
    tmp.write_text(json.dumps(meta, ensure_ascii=False, indent=1), encoding="utf-8")
    os.replace(tmp, out_dir / ccs.META_FILE)


# ---------------------------------------------------------------- 训练
def _abs_rel(pred, true):
    return np.abs(pred - true) / true


def error_summary(err):
    """绝对相对误差 -> {n, median_pct, q95_pct}。"""
    err = np.asarray(err, dtype=float)
    if not len(err):
        return {"n": 0, "median_pct": None, "q95_pct": None}
    return {"n": int(len(err)), "median_pct": round(float(np.median(err)) * 100, 3), "q95_pct": round(float(np.quantile(err, 0.95)) * 100, 3)}


def round_up(value, step):
    """向上取整到 step 的整数倍（容差：测试集 95 分位 -> 0.1 %）。"""
    return round(math.ceil(round(value / step, 6)) * step, 6)


# 候选模型（在验证集上按中位绝对相对误差选）。SVR 的 gamma 是对标准化后的特征而言（约 50 维），epsilon 是对标准化后的 ln(CCS) 而言
GRID = ([("SVR", {"C": C, "gamma": g, "epsilon": 0.02}) for C in (3.0, 10.0, 30.0) for g in (0.002, 0.005, 0.01, 0.02)]
        + [("HistGradientBoosting", {"learning_rate": lr, "max_iter": it, "max_leaf_nodes": 31, "l2_regularization": 1.0}) for lr in (0.05, 0.1) for it in (300, 800)])


def _pct(v):
    return "—" if v is None else f"{v:.2f}"


def train(rows, out_dir, cfg=CONFIG, log=print, grid=GRID):
    """训练预测模型（HistGradientBoosting 与 RBF-SVR 在验证集上选），在测试集上评估，写 ccs_model.joblib / ccs_domain.npz。返回 meta["model"]。

    划分按连接键分组：训练 / 验证 / 测试里没有同一个化合物。验证集选模型，选好后用训练 + 验证集重训（这就是装上的模型），测试集只评估一次。
    """
    import joblib
    import sklearn
    from sklearn.compose import TransformedTargetRegressor
    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.model_selection import GroupShuffleSplit
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import FunctionTransformer, StandardScaler
    from sklearn.svm import SVR

    started = time.time()
    descs = {}
    for r in rows:
        if r.smiles not in descs:
            descs[r.smiles] = ccs.descriptor_vector(ci.normalize(r.smiles).mol)
    keys = np.array([r.key for r in rows])
    y = np.array([r.ccs for r in rows])
    types = np.array([r.ccs_type for r in rows])
    adduct_of = np.array([r.adduct for r in rows])

    split = GroupShuffleSplit(n_splits=1, test_size=cfg.test_fraction, random_state=cfg.seed)
    trainval, test = next(split.split(np.zeros(len(rows)), groups=keys))
    inner = GroupShuffleSplit(n_splits=1, test_size=cfg.val_fraction / (1 - cfg.test_fraction), random_state=cfg.seed)
    tr, va = next(inner.split(np.zeros(len(trainval)), groups=keys[trainval]))
    train_idx, val_idx = trainval[tr], trainval[va]
    assert not (set(keys[trainval]) & set(keys[test])) and not (set(keys[train_idx]) & set(keys[val_idx]))

    counts = {a: int(n) for a, n in sorted(zip(*np.unique(adduct_of[trainval], return_counts=True)), key=lambda x: -x[1])}
    adducts = {a: n for a, n in counts.items() if n >= cfg.min_adduct_rows}
    log(f"Split by compound: train {len(train_idx)} / validation {len(val_idx)} / test {len(test)} rows; adducts with >= {cfg.min_adduct_rows} "
        f"training rows: {', '.join(adducts)}")
    adduct_list = list(adducts)

    def matrix(idx):
        return np.array([ccs.feature_vector(descs[rows[i].smiles], rows[i].mz, rows[i].adduct, adduct_list) for i in idx])

    keep = lambda idx: np.array([i for i in idx if adduct_of[i] in adducts], dtype=int)      # noqa: E731
    train_idx, val_idx, trainval_k, test_k = keep(train_idx), keep(val_idx), keep(trainval), keep(test)
    X = {name: matrix(idx) for name, idx in (("train", train_idx), ("val", val_idx), ("trainval", trainval_k), ("test", test_k))}

    def make(kind, params):
        target = Pipeline([("log", FunctionTransformer(np.log, np.exp, validate=True)), ("scale", StandardScaler())])
        if kind == "SVR":
            reg = Pipeline([("scale", StandardScaler()), ("svr", SVR(kernel="rbf", **params))])
        else:
            reg = HistGradientBoostingRegressor(random_state=cfg.seed, **params)
        return TransformedTargetRegressor(regressor=reg, transformer=target)

    results = []
    for kind, params in grid:
        est = make(kind, params).fit(X["train"], y[train_idx])
        err = error_summary(_abs_rel(est.predict(X["val"]), y[val_idx]))
        results.append({"name": kind, "params": params, "val_median_pct": err["median_pct"], "val_q95_pct": err["q95_pct"]})
        log(f"  validation {kind} {params}: median {_pct(err['median_pct'])} %, q95 {_pct(err['q95_pct'])} %")
    best = min(results, key=lambda r: (r["val_median_pct"], r["val_q95_pct"]))
    log(f"Chosen: {best['name']} {best['params']} (validation median {_pct(best['val_median_pct'])} %)")

    est = make(best["name"], best["params"]).fit(X["trainval"], y[trainval_k])         # 选好以后用训练 + 验证集重训，测试集只评估一次
    err = _abs_rel(est.predict(X["test"]), y[test_k])
    test_summary = error_summary(err)
    test_summary["by_type"] = {t: error_summary(err[types[test_k] == t]) for t in CCS_TYPES}
    test_summary["rows_without_model_adduct"] = int(len(test) - len(test_k))
    tolerance_pct = round_up(float(np.quantile(err, 0.95)) * 100, cfg.tol_round)         # 用没取整的 95 分位，再向上取整到 0.1 %

    # 训练集指纹 + 测试集上的领域分析
    train_keys = sorted({keys[i] for i in trainval_k})
    smiles_of = {}
    for i in trainval_k:
        smiles_of.setdefault(keys[i], rows[i].smiles)
    fps = np.stack([np.frombuffer(ci.fingerprint(ci.normalize(smiles_of[k]).mol), dtype=np.uint8) for k in train_keys])
    big = fps.view(np.uint64)
    pops = ci.popcount(big).sum(axis=1).astype(np.uint16)
    sim_of = {}
    for i in test_k:
        if keys[i] not in sim_of:
            q = np.frombuffer(ci.fingerprint(ci.normalize(rows[i].smiles).mol), dtype=np.uint8).view(np.uint64)
            sim_of[keys[i]] = float(ci.tanimoto(big, pops, q, int(ci.popcount(q).sum())).max())
    sims = np.array([sim_of[keys[i]] for i in test_k])
    inside = sims >= cfg.domain_min_tanimoto
    edges = [0.0, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0001]
    bins = [{"from": lo, "to": min(hi, 1.0), **error_summary(err[(sims >= lo) & (sims < hi)])} for lo, hi in zip(edges, edges[1:])]
    domain = {"threshold": cfg.domain_min_tanimoto, "train_compounds": len(train_keys),
              "test_inside": error_summary(err[inside]), "test_outside": error_summary(err[~inside]), "test_bins": bins,
              "test_max_tanimoto_median": round(float(np.median(sims)), 3)}

    out_dir = Path(out_dir)
    joblib.dump(est, out_dir / (ccs.MODEL_FILE + ".part"))
    os.replace(out_dir / (ccs.MODEL_FILE + ".part"), out_dir / ccs.MODEL_FILE)
    with open(out_dir / (ccs.DOMAIN_FILE + ".part"), "wb") as f:
        np.savez_compressed(f, fps=fps)
    os.replace(out_dir / (ccs.DOMAIN_FILE + ".part"), out_dir / ccs.DOMAIN_FILE)
    info = {"name": best["name"], "params": best["params"], "features": "m/z, adduct one-hot, 42 MQNs", "target": "ln(CCS), standardised",
            "adducts": adducts, "min_adduct_rows": cfg.min_adduct_rows, "sklearn": sklearn.__version__, "seed": cfg.seed,
            "split": {"by": "connectivity key", "train_rows": int(len(train_idx)), "val_rows": int(len(val_idx)), "test_rows": int(len(test_k)),
                      "trainval_compounds": len(train_keys), "test_compounds": len(sim_of)},
            "validation": results, "test": test_summary, "tolerance_pct": tolerance_pct, "domain": domain,
            "seconds": round(time.time() - started, 1)}
    log(f"Test set ({test_summary['n']} rows, held-out compounds): median {_pct(test_summary['median_pct'])} %, q95 {_pct(test_summary['q95_pct'])} % "
        f"-> prediction tolerance ±{tolerance_pct:g} %; " + "; ".join(f"{t} n={s['n']} median {_pct(s['median_pct'])} % q95 {_pct(s['q95_pct'])} %"
                                                                         for t, s in test_summary["by_type"].items()))
    log(f"Domain (max Tanimoto to the training compounds, threshold {cfg.domain_min_tanimoto:g}): inside n={domain['test_inside']['n']} "
        f"median {_pct(domain['test_inside']['median_pct'])} % q95 {_pct(domain['test_inside']['q95_pct'])} %; outside n={domain['test_outside']['n']} "
        f"median {_pct(domain['test_outside']['median_pct'])} % q95 {_pct(domain['test_outside']['q95_pct'])} %")
    return info


# ---------------------------------------------------------------- 主流程
def setup(force=False, clean=False, skip_model=False):
    out = ccs.ccs_dir()
    if ccs.installed(out) and not force:
        print(f"The CCS index already exists in {out} (use --force to rebuild).")
        return out
    started = time.time()
    repo = raw_dir()
    fetch(repo)
    src = repo / "c3sdb" / "_include" / "src_data"
    n_datasets = sum(1 for _ in src.glob("*.json")) if src.is_dir() else 0
    if n_datasets < MIN_DATASETS:
        raise SystemExit(f"Only {n_datasets} datasets found in {src} (expected 26). Delete {repo} and run again.")
    datasets, cache = read_source(repo)
    rows, sources, stats = build_rows(datasets, cache)
    out.mkdir(parents=True, exist_ok=True)
    ccs.save_index(out / ccs.INDEX_FILE, rows, sources)
    meta = {"source": "CCSbase source data (c3sdb)", "repo": ccs.REPO_URL, "commit": COMMIT, "licence": ccs.LICENCE, "licence_url": ccs.LICENCE_URL,
            "citation": ccs.CITATION, "built": date.today().isoformat(), "schema": ccs.SCHEMA, **stats}
    write_attribution(out, repo, stats)
    write_meta(out, meta)                                       # 先写不带模型的 meta：索引已经可用
    if skip_model:
        print("Model training skipped; the page then gives no predictions.")
    else:
        try:
            meta["model"] = train(rows, out)
        except ImportError as exc:
            print(f"scikit-learn is needed to train the prediction model ({exc}). The index works without it; install scikit-learn and run with --force.")
        write_meta(out, meta)
    if clean:
        shutil.rmtree(repo, ignore_errors=True)
        print(f"Removed {repo}.")
    print(f"CCS index: {stats['rows_kept']} rows, {stats['compounds']} compounds, {stats['datasets']} datasets in {time.time() - started:.0f} s ({out}).")
    print(f"The data are {ccs.LICENCE} (attribution in {out / 'ATTRIBUTION.md'}). Citation: {ccs.CITATION}")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="Download the CCSbase source data (c3sdb, pinned commit), build the local CCS index and train the prediction model.")
    ap.add_argument("--force", action="store_true", help="rebuild the index and the model even if they exist (reuses the downloaded files)")
    ap.add_argument("--clean", action="store_true", help="delete the downloaded files after the index is built")
    ap.add_argument("--skip-model", action="store_true", help="build the index only (no scikit-learn needed; no predictions)")
    args = ap.parse_args(argv)
    setup(args.force, args.clean, args.skip_model)


if __name__ == "__main__":
    sys.exit(main())
