"""谱库本地索引：建索引、按需打开（mmap）、取候选谱。

目录 ~/.mzdesk/data/libraries/<库名>/（环境变量 MZDESK_LIBRARY_DIR 改根目录），文件：
  meta.json      来源 URL、版本、日期、条数、加合物词表……
  prec_mz.npy    前体 m/z（float64）。记录按 (kind, prec_mz) 排序：先全部 MS/MS（前体升序），再全部 EI
  kind.npy       0 = MS/MS，1 = EI（int8）
  ion_mode.npy   +1 / -1 / 0 = 未知（int8）
  adduct_id.npy  规范加合物在 meta.json 词表里的序号，-1 = 未知（int16）
  copy_flag.npy  1 = 这是 GNPS 里转存的 MassBank 记录（副本；原件的 ID 在 meta.sqlite 的 copy_of 列），否则 0（int8）
  offsets.npy    每条谱的峰起点，长度 N+1（int64）
  peak_mz.npy    所有峰的 m/z（float64，每条谱内升序）
  peak_int.npy   所有峰的相对强度（float32，每条谱基峰 = 1）
  meta.sqlite    名称、登录号、InChIKey、SMILES、加合物、仪器、碰撞能量；id = 上面数组里的行号
数组用 np.load(mmap_mode="r") 打开，常驻内存只有元数据；一次检索只读到候选谱的那几页。
"""
import json
import os
import shutil
import sqlite3
import time
from pathlib import Path

import numpy as np

from mzdesk.ms2.library_parse import EI, MSMS

FORMAT = 2
SETUP_COMMAND = "python -m mzdesk.ms2.library_setup"
COPY_SOURCE = "MassBank"            # copy_flag 标的副本是哪个库的记录
ARRAYS = ("prec_mz", "kind", "ion_mode", "adduct_id", "copy_flag", "offsets", "peak_mz", "peak_int")
_COLUMNS = ("accession", "name", "inchikey", "smiles", "adduct", "instrument", "collision_energy", "copy_of")


def library_dir():
    return Path(os.environ.get("MZDESK_LIBRARY_DIR") or Path.home() / ".mzdesk" / "data" / "libraries")


def index_dir(name):
    return library_dir() / name


def short_id(accession):
    """MassBank 登录号 MSBNK-LCSB-LU119903 -> LU119903（GNPS 副本的名字里写的就是这一段）。"""
    parts = str(accession).split("-", 2)
    return parts[2] if len(parts) == 3 else str(accession)


class LibraryMissing(RuntimeError):
    """索引没装（或格式过旧）。消息直接给用户看，里面带安装命令。"""


# ---------------------------------------------------------------- 建索引
def _write_sorted_peaks(tmp_dir, out_dir, order, npk):
    """按 order（新顺序 -> 原序号）把临时峰文件重排成最终的 peak_*.npy；分块做，内存只占一块。"""
    n = len(order)
    off_src = np.zeros(len(npk) + 1, dtype=np.int64)
    np.cumsum(npk, out=off_src[1:])
    total = int(off_src[-1])
    mz_src = np.memmap(tmp_dir / "peak_mz.bin", dtype=np.float64, mode="r", shape=(total,)) if total else np.zeros(0)
    in_src = np.memmap(tmp_dir / "peak_int.bin", dtype=np.float32, mode="r", shape=(total,)) if total else np.zeros(0, np.float32)
    lens = npk[order]
    offsets = np.zeros(n + 1, dtype=np.int64)
    np.cumsum(lens, out=offsets[1:])
    mz_out = np.lib.format.open_memmap(out_dir / "peak_mz.npy", mode="w+", dtype=np.float64, shape=(total,))
    in_out = np.lib.format.open_memmap(out_dir / "peak_int.npy", mode="w+", dtype=np.float32, shape=(total,))
    chunk = 20000
    for a in range(0, n, chunk):
        b = min(n, a + chunk)
        ids = order[a:b]
        starts, ln = off_src[ids], lens[a:b]
        flat_len = int(ln.sum())
        rel = np.arange(flat_len) - np.repeat(np.cumsum(ln) - ln, ln)
        idx = np.repeat(starts, ln) + rel
        lo, hi = int(offsets[a]), int(offsets[b])
        mz_out[lo:hi] = mz_src[idx]
        in_out[lo:hi] = in_src[idx]
    mz_out.flush()
    in_out.flush()
    del mz_out, in_out, mz_src, in_src
    return offsets, total


def build(name, records, info, progress=None):
    """records（library_parse 产出的记录，已规范化）-> 索引目录。info 写进 meta.json。返回 meta。

    先写到同级的临时目录，成功后整体替换，所以建到一半失败不会弄坏已有的索引。
    """
    out = index_dir(name)
    out.parent.mkdir(parents=True, exist_ok=True)
    work = out.parent / f".{name}.building"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir()
    try:
        meta = _build(name, records, info, work, progress)
    except BaseException:
        shutil.rmtree(work, ignore_errors=True)
        raise
    old = out.parent / f".{name}.old"
    shutil.rmtree(old, ignore_errors=True)
    if out.exists():
        out.rename(old)
    work.rename(out)
    shutil.rmtree(old, ignore_errors=True)
    clear_cache()
    return meta


def _build(name, records, info, work, progress):
    tmp = work / "tmp"
    tmp.mkdir()
    db = sqlite3.connect(tmp / "build.sqlite")
    db.execute("CREATE TABLE t (src INTEGER PRIMARY KEY, accession, name, inchikey, smiles, adduct, instrument, collision_energy, copy_of)")
    prec, kind, mode, aid, npk, copy = [], [], [], [], [], []
    vocab = {}
    batch = []
    n = 0
    t0 = time.time()
    with open(tmp / "peak_mz.bin", "wb") as f_mz, open(tmp / "peak_int.bin", "wb") as f_in:
        for rec in records:
            rec["mz"].tofile(f_mz)
            rec["intensity"].tofile(f_in)
            prec.append(rec["precursor_mz"] or 0.0)
            kind.append(rec["kind"])
            mode.append(rec["ion_mode"])
            canon = rec["adduct"]
            aid.append(-1 if canon is None else vocab.setdefault(canon, len(vocab)))
            npk.append(len(rec["mz"]))
            copy.append(0 if rec.get("copy_of") is None else 1)
            batch.append((n, rec["accession"], rec["name"], rec["inchikey"], rec["smiles"], canon or rec["adduct_raw"] or "",
                          rec["instrument"], rec["collision_energy"], rec.get("copy_of") or ""))
            n += 1
            if len(batch) >= 5000:
                db.executemany("INSERT INTO t VALUES (?,?,?,?,?,?,?,?,?)", batch)
                batch = []
            if progress and n % 50000 == 0:
                progress(f"  {n:,} spectra kept ({time.time() - t0:.0f} s)")
    if batch:
        db.executemany("INSERT INTO t VALUES (?,?,?,?,?,?,?,?,?)", batch)
    db.commit()
    if n == 0:
        raise ValueError(f"{name}: no usable spectra in the source file.")

    prec_a, kind_a = np.asarray(prec, np.float64), np.asarray(kind, np.int8)
    mode_a, aid_a, npk_a = np.asarray(mode, np.int8), np.asarray(aid, np.int16), np.asarray(npk, np.int64)
    copy_a = np.asarray(copy, np.int8)
    # 排序：MS/MS 在前，按前体升序；EI 在后（前体都是 0，按原顺序）
    order = np.lexsort((np.arange(n), prec_a, kind_a))
    np.save(work / "prec_mz.npy", prec_a[order])
    np.save(work / "kind.npy", kind_a[order])
    np.save(work / "ion_mode.npy", mode_a[order])
    np.save(work / "adduct_id.npy", aid_a[order])
    np.save(work / "copy_flag.npy", copy_a[order])
    offsets, n_peaks = _write_sorted_peaks(tmp, work, order, npk_a)
    np.save(work / "offsets.npy", offsets)

    # 元数据：按新行号重排进 meta.sqlite
    db.execute("CREATE TABLE perm (new INTEGER PRIMARY KEY, src INTEGER)")
    db.executemany("INSERT INTO perm VALUES (?,?)", enumerate(order.tolist()))
    db.execute(f"ATTACH DATABASE ? AS final", (str(work / "meta.sqlite"),))
    db.execute("CREATE TABLE final.meta (id INTEGER PRIMARY KEY, accession, name, inchikey, smiles, adduct, instrument, collision_energy, copy_of)")
    db.execute("INSERT INTO final.meta SELECT perm.new, t.accession, t.name, t.inchikey, t.smiles, t.adduct, t.instrument, t.collision_energy, t.copy_of "
               "FROM perm JOIN t ON t.src = perm.src ORDER BY perm.new")
    db.commit()
    db.close()
    shutil.rmtree(tmp)

    sk = kind_a[order]
    ms = mode_a[order]
    n_msms = int((sk == MSMS).sum())
    counts = {"msms_positive": int(((sk == MSMS) & (ms > 0)).sum()), "msms_negative": int(((sk == MSMS) & (ms < 0)).sum()),
              "msms_unknown_mode": int(((sk == MSMS) & (ms == 0)).sum()), "ei": int((sk == EI).sum()),
              "msms_without_adduct": int(((sk == MSMS) & (aid_a[order] < 0)).sum()), "massbank_copies": int(copy_a.sum())}
    meta = {"format": FORMAT, "library": name, **info, "built": time.strftime("%Y-%m-%d %H:%M:%S"), "n_spectra": n,
            "n_msms": n_msms, "n_peaks": n_peaks, "counts": counts, "adducts": list(vocab)}
    (work / "meta.json").write_text(json.dumps(meta, indent=1, ensure_ascii=False))
    return meta


def folder_size(path):
    return sum(f.stat().st_size for f in Path(path).rglob("*") if f.is_file())


# ---------------------------------------------------------------- 打开 / 检索
class Library:
    """一个已建好的索引（数组用 mmap 打开）。"""

    def __init__(self, name, path, meta):
        self.name, self.path, self.meta = name, path, meta
        a = {k: np.load(path / f"{k}.npy", mmap_mode="r") for k in ARRAYS}
        self.prec_mz, self.kind, self.ion_mode, self.adduct_id = a["prec_mz"], a["kind"], a["ion_mode"], a["adduct_id"]
        self.copy_flag = a["copy_flag"]
        self._short_ids = None
        self.offsets, self.peak_mz, self.peak_int = a["offsets"], a["peak_mz"], a["peak_int"]
        self.n = int(meta["n_spectra"])
        self.n_msms = int(meta["n_msms"])
        self._adduct_ids = {a: i for i, a in enumerate(meta["adducts"])}

    @property
    def version(self):
        return str(self.meta.get("version") or "")

    def label(self):
        return f"{self.name} {self.version}".strip()

    def adduct_id_of(self, canonical):
        """规范加合物 -> 词表序号；库里一条也没有返回 -1。"""
        return self._adduct_ids.get(canonical, -1)

    def msms_candidates(self, precursor_mz, tol_ppm, positive, adduct=None):
        """前体窗口内、同离子模式（给了加合物就同加合物）的 MS/MS 行号，和 Δppm = (查询 − 库) / 库。

        窗口以库里的质量为基准：|查询 − 库| / 库 <= tol_ppm。
        """
        t = tol_ppm * 1e-6
        seg = self.prec_mz[:self.n_msms]
        lo = int(np.searchsorted(seg, precursor_mz / (1 + t), side="left"))
        hi = int(np.searchsorted(seg, precursor_mz / (1 - t), side="right"))
        ids = np.arange(lo, hi)
        if ids.size:
            keep = np.asarray(self.ion_mode[lo:hi]) == (1 if positive else -1)
            if adduct is not None:
                wanted = self.adduct_id_of(adduct)
                keep &= (np.asarray(self.adduct_id[lo:hi]) == wanted) if wanted >= 0 else np.zeros(hi - lo, dtype=bool)
            ids = ids[keep]
        lib_mz = np.asarray(self.prec_mz[ids])
        return ids, (precursor_mz - lib_mz) / lib_mz * 1e6

    def short_id_rows(self, short_ids):
        """MassBank 的短 ID（登录号去掉 "MSBNK-<机构>-" 前缀，如 LU119903）-> 行号列表。第一次用时从 meta.sqlite 读一遍登录号并缓存。"""
        if self._short_ids is None:
            con = sqlite3.connect(f"{(self.path / 'meta.sqlite').resolve().as_uri()}?mode=ro", uri=True)
            try:
                table = {}
                for row, acc in con.execute("SELECT id, accession FROM meta"):
                    table.setdefault(short_id(acc), []).append(row)
            finally:
                con.close()
            self._short_ids = table
        return {s: self._short_ids.get(s, []) for s in short_ids}

    def ei_ids(self):
        return np.arange(self.n_msms, self.n)

    def gather(self, ids):
        """行号 -> (mz 平铺, 相对强度平铺, 每条谱的峰数)。"""
        ids = np.asarray(ids, dtype=np.int64)
        starts = np.asarray(self.offsets[ids])
        lens = np.asarray(self.offsets[ids + 1]) - starts
        total = int(lens.sum())
        if total == 0:
            return np.zeros(0), np.zeros(0, np.float32), lens
        if ids.size and np.array_equal(ids, np.arange(ids[0], ids[0] + ids.size)):      # 连续行（EI 全库）：直接切片
            lo = int(starts[0])
            return np.asarray(self.peak_mz[lo:lo + total]), np.asarray(self.peak_int[lo:lo + total]), lens
        idx = np.repeat(starts, lens) + (np.arange(total) - np.repeat(np.cumsum(lens) - lens, lens))
        return np.asarray(self.peak_mz[idx]), np.asarray(self.peak_int[idx]), lens

    def peaks(self, i):
        a, b = int(self.offsets[i]), int(self.offsets[i + 1])
        return np.asarray(self.peak_mz[a:b]), np.asarray(self.peak_int[a:b])

    def info(self, ids):
        """行号 -> 元数据字典列表（顺序同 ids）。"""
        ids = [int(i) for i in ids]
        if not ids:
            return []
        con = sqlite3.connect(f"{(self.path / 'meta.sqlite').resolve().as_uri()}?mode=ro", uri=True)       # as_uri 转义路径里的空格、# 等
        rows = {}
        try:
            for a in range(0, len(ids), 500):                      # sqlite 对 ? 的个数有上限，分块查
                part = ids[a:a + 500]
                marks = ",".join("?" * len(part))
                rows.update({r[0]: dict(zip(_COLUMNS, r[1:])) for r in
                             con.execute(f"SELECT id, {', '.join(_COLUMNS)} FROM meta WHERE id IN ({marks})", part)})
        finally:
            con.close()
        return [rows.get(i, dict.fromkeys(_COLUMNS, "")) for i in ids]


_CACHE = {}


def clear_cache():
    _CACHE.clear()


def is_installed(name):
    return (index_dir(name) / "meta.json").is_file()


def load(name):
    """已建好的索引，进程内缓存（meta.json 变了会重新打开）。没装 -> LibraryMissing。"""
    path = index_dir(name)
    meta_file = path / "meta.json"
    if not meta_file.is_file():
        raise LibraryMissing(f"{name} is not installed. Run {SETUP_COMMAND} {name}")
    key = (str(path), meta_file.stat().st_mtime_ns)
    lib = _CACHE.get(key)
    if lib is None:
        meta = json.loads(meta_file.read_text())
        if meta.get("format") != FORMAT:
            raise LibraryMissing(f"The {name} index has an old format. Run {SETUP_COMMAND} {name} --rebuild")
        try:
            lib = Library(name, path, meta)
        except (OSError, ValueError, KeyError) as exc:
            raise LibraryMissing(f"The {name} index is incomplete ({exc}). Run {SETUP_COMMAND} {name} --force") from None
        for k in [k for k in _CACHE if k[0] == key[0]]:
            del _CACHE[k]
        _CACHE[key] = lib
    return lib
