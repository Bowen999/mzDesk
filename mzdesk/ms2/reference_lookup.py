"""参考谱查询：一个结构 -> 本地 MassBank / GNPS 索引里同一化合物的全部实测谱（不要查询谱）。

同一化合物 = InChIKey 前 14 位（2D 结构块）相同：立体异构体、同位素以外的差别都算不同化合物，立体构型不同的记录也列出来，
"Stereo" 列说明完整 InChIKey 是否也相同。盐 / 混合物取最大的组分（与 MRM 页、各引擎同一规则：mzdesk.ms2.mrm.structure_key）。
库里没有 InChIKey 的记录查不到（GNPS 约 20 % 的记录没有）。

同时选了 MassBank 和 GNPS 时，GNPS 里转存的 MassBank 记录只要原件也在结果里就隐藏（规则同 library_search._hide_massbank_copies）。
给了查询谱时，每条参考谱都算加权谱熵相似度和 CosineGreedy（MS/MS 用所给碎片容差，EI 固定 0.5 Da），结果按谱熵降序。
"""
from __future__ import annotations

import math

import numpy as np

from mzdesk import adducts
from mzdesk.ms2 import library_index as lx
from mzdesk.ms2 import library_search as ls
from mzdesk.ms2 import metric_names as mn
from mzdesk.ms2 import ms2_similarity as sim

KINDS = ("MS/MS", "EI")
ION_MODES = ("Any", "Positive", "Negative")
METRICS = (mn.ENTROPY_WEIGHTED, mn.COSINE_GREEDY)
_MODE_NAME = {1: "Positive", -1: "Negative", 0: ""}


def lookup(smiles, libraries, *, kinds=KINDS, ion_mode="Any", query=None, tolerance=0.02):
    """-> {"key": 完整 InChIKey, "key14", "records": [...], "libraries": [{"name", "version", "found"}], "hidden_copies", "warnings"}。

    records 的每项：library, row, accession, name, kind ("MS/MS" / "EI"), ion_mode, adduct, precursor_mz（EI 为 None）, instrument,
    collision_energy, inchikey, same_stereo, smiles, copy_of, n_peaks, scores（{指标: 分数或 None}，没有查询谱时为空）, matches（同上，None）。
    query：[(mz, intensity)] 或 None。库全没装 -> lx.LibraryMissing；ValueError = 输入不对。
    """
    from rdkit import Chem
    from rdkit.Chem import inchi as rd_inchi

    from mzdesk.ms2.mrm import structure_key

    kinds = tuple(k for k in KINDS if k in set(kinds or ()))
    if not kinds:
        raise ValueError("Select at least one spectrum type.")
    if ion_mode not in ION_MODES:
        raise ValueError(f"Ion mode must be one of {', '.join(ION_MODES)}.")
    key14 = structure_key(smiles)
    mol = Chem.MolFromSmiles(smiles)
    parts = Chem.GetMolFrags(mol, asMols=True)
    if len(parts) > 1:
        from rdkit.Chem import Descriptors
        mol = max(parts, key=lambda part: (part.GetNumHeavyAtoms(), Descriptors.ExactMolWt(part)))
    full_key = rd_inchi.MolToInchiKey(mol) or ""

    warnings, loaded, missing = [], [], []
    for name in dict.fromkeys(libraries):
        try:
            loaded.append(lx.load(name))
        except lx.LibraryMissing as exc:
            missing.append(str(exc))
    if not loaded:
        raise lx.LibraryMissing(" ".join(missing) if missing else "Select at least one library.")
    warnings += [f"{msg} Searching {', '.join(l.name for l in loaded)} only." for msg in missing]

    want_kind = {0 if k == "MS/MS" else 1 for k in kinds}
    want_mode = {"Positive": 1, "Negative": -1}.get(ion_mode)
    found = []
    for lib in loaded:
        ids = lib.structure_rows(key14)
        if ids.size:
            keep = np.isin(np.asarray(lib.kind[ids]), list(want_kind))
            if want_mode is not None:          # 离子模式只筛 MS/MS；EI 谱都是正离子 [M]+•
                keep &= (np.asarray(lib.kind[ids]) == 1) | (np.asarray(lib.ion_mode[ids]) == want_mode)
            ids = ids[keep]
        found.append((lib, ids, np.zeros(ids.size)))
    found, hidden = ls._hide_massbank_copies(found)

    q = _scorers(query, tolerance)

    records, summary = [], []
    for lib, ids, _ in found:
        summary.append({"name": lib.name, "version": lib.version, "found": int(ids.size)})
        if not ids.size:
            continue
        mz_flat, rel_flat, lens = lib.gather(ids)
        off = np.concatenate([[0], np.cumsum(lens)])
        for k, (row, info) in enumerate(zip(ids.tolist(), lib.info(ids))):
            rec = _record(lib, row, info, int(lens[k]), full_key)
            if q is not None and lens[k]:
                _score(rec, q, mz_flat[off[k]:off[k + 1]], rel_flat[off[k]:off[k + 1]])
            records.append(rec)

    lib_order = {l.name: i for i, l in enumerate(loaded)}
    mode_order = {"Positive": 0, "Negative": 1, "": 2}

    def plain(r):
        # within a library and type: ion mode, records with an adduct before those without, then more peaks first
        # (the first row is what the page plots, so it should be an informative spectrum)
        return (lib_order[r["library"]], r["kind"] != "MS/MS", mode_order.get(r["ion_mode"], 2), r["adduct"] == "", r["adduct"],
                -r["n_peaks"], r["accession"])

    if q is not None:
        records.sort(key=lambda r: (-(r["scores"].get(mn.ENTROPY_WEIGHTED) if r["scores"].get(mn.ENTROPY_WEIGHTED) is not None else -1.0),
                                    -(r["matches"] or 0), plain(r)))
    else:
        records.sort(key=plain)
    return {"key": full_key, "key14": key14, "records": records, "libraries": summary, "hidden_copies": hidden, "warnings": warnings}


def _scorers(query, tolerance):
    """查询谱 -> MS/MS 与 EI 两个打分器（同 lookup）；没有查询谱 -> None。"""
    if not query:
        return None
    q_mz, q_int = sim.clean_peaks(*sim.normalize_peaks_input(query))
    tol = float(tolerance)
    if not (math.isfinite(tol) and tol > 0):
        raise ValueError("Fragment tolerance must be a positive number.")
    return {"msms": ls._Scorer(q_mz, q_int, None, tol, list(METRICS), False),
            "ei": ls._Scorer(q_mz, q_int, None, ls.EI_TOLERANCE, list(METRICS), True)}


def _record(lib, row, info, n_peaks, full_key):
    """库里一行 -> lookup 的一条记录（scores 空，matches None）。"""
    ei = int(lib.kind[row]) == 1
    return {"library": lib.name, "row": int(row), "accession": info["accession"], "name": info["name"],
            "kind": "EI" if ei else "MS/MS", "ion_mode": _MODE_NAME.get(int(lib.ion_mode[row]), ""),
            "adduct": "" if ei else (adducts.canonical(info["adduct"]) or info["adduct"] or ""),
            "precursor_mz": None if ei or not lib.prec_mz[row] > 0 else float(lib.prec_mz[row]),
            "instrument": info["instrument"] or "", "collision_energy": info["collision_energy"] or "",
            "inchikey": info["inchikey"] or "", "same_stereo": bool(full_key) and info["inchikey"] == full_key,
            "smiles": info.get("smiles") or "", "copy_of": info.get("copy_of") or "",
            "n_peaks": int(n_peaks), "scores": {}, "matches": None}


def _score(rec, q, mz, rel):
    """给一条记录打分（就地写 scores 和 matches）。rel：相对强度，基峰 = 1。"""
    got, rec["matches"] = q["ei" if rec["kind"] == "EI" else "msms"].score(mz, np.asarray(rel).astype(np.float64) * 100.0, None)
    rec["scores"] = {m: (v if v is not None and math.isfinite(v) else None) for m, v in got.items() if m in METRICS}


def record(library, row, accession, *, inchikey="", query=None, tolerance=0.02):
    """表里点开的一条参考谱：库名 + 行号 -> 与 lookup 同样的记录（给了查询谱就打分）。

    accession 必须与该行一致，否则 LookupError（库在查询和点开之间换了版本）。inchikey：查询结构的完整 InChIKey，决定 same_stereo。
    库没装 -> lx.LibraryMissing；ValueError = 输入不对。
    """
    lib = lx.load(library)
    row = int(row)
    if not 0 <= row < lib.n:
        raise LookupError(f"{library} has no record number {row}. Run the page again.")
    info = lib.info([row])[0]
    if info["accession"] != accession:
        raise LookupError(f"{library} record {accession} is no longer at this place; the library was updated. Run the page again.")
    mz, rel = lib.peaks(row)
    rec = _record(lib, row, info, len(mz), inchikey or "")
    q = _scorers(query, tolerance)
    if q is not None and len(mz):
        _score(rec, q, mz, rel)
    return rec


def peaks_of(record):
    """record -> [(mz, 相对强度, 基峰 = 100)]。"""
    lib = lx.load(record["library"])
    mz, rel = lib.peaks(record["row"])
    return [(float(m), round(float(r) * 100.0, 2)) for m, r in zip(mz, rel)]
