"""谱库下载与建索引，一次性运行：python -m mzdesk.ms2.library_setup [MassBank] [GNPS]

MassBank：GitHub 最新 release 的 MassBank_NISTformat.msp（约 137 MB；NIST 格式带 Spectrum_type / 仪器 / 碰撞能量）。
GNPS：ALL_GNPS_NO_PROPOGATED.mgf（约 2.2 GB，只含原始提交的谱，不含传播的）。
原始文件放 ~/.mzdesk/data/libraries/raw/（环境变量 MZDESK_LIBRARY_DIR 改根目录，支持断点续传），
索引放 <根目录>/<库名>/（格式见 library_index.py），meta.json 记录来源 URL、版本、日期。数据不入库。
不带库名 = 两个都装。--rebuild 用已下载的原始文件重建索引；--force 重新下载再建；--delete-raw 建好后删原始文件。
"""
import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

from mzdesk.ms2 import library_index as lx
from mzdesk.ms2 import library_parse as lp

MASSBANK_API = "https://api.github.com/repos/MassBank/MassBank-data/releases/latest"
MASSBANK_ASSET = "MassBank_NISTformat.msp"
GNPS_URL = "https://external.gnps2.org/gnpslibrary/ALL_GNPS_NO_PROPOGATED.mgf"
GNPS_FILE = "ALL_GNPS_NO_PROPOGATED.mgf"
NAMES = ("MassBank", "GNPS")
_UA = {"User-Agent": "mzdesk-library-setup"}


def raw_dir():
    return lx.library_dir() / "raw"


def _request(url, headers=None, method="GET"):
    return urllib.request.Request(url, headers={**_UA, **(headers or {})}, method=method)


def source_info(name):
    """{"url", "file", "version", "date", "size"}（size 可能是 None）；要联网。"""
    if name == "MassBank":
        with urllib.request.urlopen(_request(MASSBANK_API), timeout=30) as resp:
            release = json.load(resp)
        asset = next((a for a in release.get("assets", []) if a["name"] == MASSBANK_ASSET), None)
        if asset is None:
            raise SystemExit(f"{MASSBANK_ASSET} is not in the latest MassBank release ({release.get('tag_name')}).")
        return {"url": asset["browser_download_url"], "file": MASSBANK_ASSET, "version": release["tag_name"],
                "date": str(release.get("published_at", ""))[:10], "size": asset.get("size")}
    with urllib.request.urlopen(_request(GNPS_URL, method="HEAD"), timeout=30) as resp:
        modified = resp.headers.get("Last-Modified", "")
        size = resp.headers.get("Content-Length")
    try:
        date = time.strftime("%Y-%m-%d", time.strptime(modified, "%a, %d %b %Y %H:%M:%S %Z"))
    except ValueError:
        date = modified
    return {"url": GNPS_URL, "file": GNPS_FILE, "version": date, "date": date, "size": int(size) if size else None}


def _download(url, dest, size=None):
    """断点续传下载：先写 dest.part，下完再改名。"""
    part = dest.with_name(dest.name + ".part")
    have = part.stat().st_size if part.exists() else 0
    if size and have > size:
        part.unlink()
        have = 0
    headers = {"Range": f"bytes={have}-"} if have else {}
    print(f"Downloading {url}\n  -> {dest}" + (f" ({size / 1e6:.0f} MB)" if size else "") + (f", resuming at {have / 1e6:.0f} MB" if have else ""), flush=True)
    try:
        resp = urllib.request.urlopen(_request(url, headers), timeout=60)
    except urllib.error.HTTPError as exc:
        if exc.code == 416 and have:                           # 已经下完，只差改名
            part.rename(dest)
            return
        raise
    mode = "ab" if have and resp.status == 206 else "wb"
    done = have if mode == "ab" else 0
    report = done + (100 << 20)
    with resp, open(part, mode) as out:
        while block := resp.read(1 << 20):
            out.write(block)
            done += len(block)
            if done >= report:
                print(f"  {done >> 20} MB" + (f" / {size >> 20} MB" if size else ""), flush=True)
                report += 100 << 20
    if size and done != size:
        raise SystemExit(f"Download stopped at {done} of {size} bytes. Run the command again to resume.")
    part.rename(dest)


def _sidecar(path):
    return path.with_name(path.name + ".json")


def get_raw(name, force=False):
    """(原始文件路径, 来源信息)。已下载过就直接用（不联网，除非旁边的 .json 不在）。"""
    raw = raw_dir()
    raw.mkdir(parents=True, exist_ok=True)
    file = MASSBANK_ASSET if name == "MassBank" else GNPS_FILE
    path = raw / file
    side = _sidecar(path)
    if path.is_file() and not force:
        if side.is_file():
            return path, json.loads(side.read_text())
        try:                                                    # 手动下的文件：补上来源信息
            info = source_info(name)
            if info["size"] in (None, path.stat().st_size):
                side.write_text(json.dumps(info, indent=1))
                return path, info
            print(f"{path} has {path.stat().st_size} bytes, the current release has {info['size']}: downloading again.")
            path.unlink()
        except OSError as exc:
            print(f"Could not read the release information ({exc}); using the file as it is.")
            return path, {"url": "", "file": file, "version": "unknown", "date": "", "size": path.stat().st_size}
    info = source_info(name)
    if force and path.exists():
        path.unlink()
    _download(info["url"], path, info["size"])
    side.write_text(json.dumps(info, indent=1))
    return path, info


def build_from(name, path, info, progress=print):
    """原始文件 -> 索引；返回 (meta, 跳过原因计数)。"""
    skipped = {}

    def kept():
        for rec, why in lp.iter_records(name, path):
            if rec is None:
                skipped[why] = skipped.get(why, 0) + 1
            else:
                yield rec

    t0 = time.time()
    # skipped 是同一个 dict：build 读完所有记录后才写 meta.json，那时计数已齐
    meta = lx.build(name, kept(), {"source_url": info.get("url", ""), "version": info.get("version", ""),
                                   "date": info.get("date", ""), "source_file": Path(path).name,
                                   "source_size": Path(path).stat().st_size, "skipped": skipped}, progress=progress)
    meta["build_seconds"] = round(time.time() - t0, 1)
    return meta


def setup(name, force=False, rebuild=False, delete_raw=False):
    if name not in NAMES:
        raise SystemExit(f"Unknown library {name!r}. Choose from: {', '.join(NAMES)}.")
    if lx.is_installed(name) and not (force or rebuild):
        print(f"{name} is already installed ({lx.index_dir(name)}). Use --rebuild to rebuild the index or --force to download it again.")
        return lx.load(name).meta
    try:
        path, info = get_raw(name, force=force)
    except OSError as exc:
        raise SystemExit(f"Could not download {name} ({exc}). Check the network connection and run the command again.") from None
    print(f"Building the {name} index from {path.name} ({info.get('version')}) ...", flush=True)
    meta = build_from(name, path, info)
    size = lx.folder_size(lx.index_dir(name))
    c = meta["counts"]
    print(f"{name} {meta.get('version')}: {meta['n_spectra']:,} spectra kept "
          f"(MS/MS positive {c['msms_positive']:,}, MS/MS negative {c['msms_negative']:,}, EI {c['ei']:,}); "
          f"index {size / 1e6:.0f} MB in {lx.index_dir(name)} ({meta['build_seconds']:.0f} s).")
    for why, n in sorted(meta.get("skipped", {}).items()):
        print(f"  skipped {n:,}: {why}")
    if delete_raw:
        path.unlink(missing_ok=True)
        _sidecar(path).unlink(missing_ok=True)
        print(f"  deleted {path.name}")
    return meta


def main(argv=None):
    ap = argparse.ArgumentParser(description="Download MassBank / GNPS and build the local spectrum-library index.")
    ap.add_argument("libraries", nargs="*", metavar="LIBRARY", help="MassBank and/or GNPS (default: both)")
    ap.add_argument("--rebuild", action="store_true", help="rebuild the index from the downloaded file")
    ap.add_argument("--force", action="store_true", help="download again and rebuild")
    ap.add_argument("--delete-raw", action="store_true", help="delete the downloaded source file after the index is built")
    args = ap.parse_args(argv)
    for name in args.libraries or NAMES:
        setup(name, force=args.force, rebuild=args.rebuild, delete_raw=args.delete_raw)


if __name__ == "__main__":
    sys.exit(main())
