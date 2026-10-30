"""谱库解析：MassBank（NIST 格式 MSP）和 GNPS（MGF）流式读取，统一成同一种记录。

文件很大（GNPS 2.2 GB），所以一次只读一条记录，不整体载入内存。只用标准库 + numpy。
每条记录 -> normalise_*() -> (record, None) 或 (None, 跳过原因)。record 是 dict：
  accession name inchikey smiles adduct(规范写法或 None) adduct_raw precursor_mz ion_mode(+1/-1/0) kind(MSMS / EI)
  instrument collision_energy mz(float64，升序) intensity(float32，基峰 = 1)
  copy_of  None = 不是副本；字符串 = GNPS 里转存的 MassBank 记录（ID 例如 "LU119903"，认不出 ID 时为 ""）
"""
import re
from pathlib import Path

import numpy as np

MSMS, EI = 0, 1
KIND_NAMES = {MSMS: "MS/MS", EI: "EI"}


# ---------------------------------------------------------------- 加合物写法统一
# 同一种加合物在不同库里写法不同：M+H / [M+H]+ / [M+H]1+；[M+2H]++ / [M+2H]2+；[M+H-H2O]+ / [M-H2O+H]+；
# 甲酸根：[M+FA-H]- / [M+HCOO]- / [M+FA]-（GNPS）/ M+formate（GNPS）；不带方括号时电荷也可能写在末尾：M+Cl- / M+FA-
_TERM = re.compile(r"([+-])(\d*)([A-Za-z][A-Za-z0-9()]*)")
_ALIAS = {"FA": "HCOOH", "HFA": "HCOOH", "HAc": "CH3COOH", "Hac": "CH3COOH", "HOAc": "CH3COOH", "AcOH": "CH3COOH",
          "Ac": "CH3COO", "HCO2": "HCOO", "HCO2H": "HCOOH", "C2H3O2": "CH3COO", "CHO2": "HCOO", "C2H4O2": "CH3COOH",
          "formate": "HCOO"}
# 加上这些就是正电荷（减去则相反）；Cl / Br / HCOO 等是负电荷
_POS = {"H": 1, "Na": 1, "K": 1, "Li": 1, "NH4": 1, "Cs": 1}
_NEG = {"Cl": 1, "Br": 1, "I": 1, "F": 1, "HCOO": 1, "CH3COO": 1, "OH": 1, "CF3COO": 1}


def _charge_of_terms(terms):
    """没写电荷时，由各项推出电荷：每个 +H / +Na 等 +1，每个 −H +(−1)，每个 +Cl 等 −1。"""
    charge = 0
    for sign, count, name in terms:
        step = 1 if sign == "+" else -1
        if name in _POS:
            charge += step * count * _POS[name]
        elif name in _NEG:
            charge -= step * count * _NEG[name]
    return charge


def normalize_adduct(text):
    """规范化的加合物写法（如 "[M+H]+"），认不出来返回 None。"""
    s = re.sub(r"\s+", "", str(text or ""))
    if not s or s.upper() in ("N/A", "NA", "NONE", "NULL", "UNKNOWN"):
        return None
    s = s.rstrip("*•.")                                       # 自由基标记 [M]+*
    m = re.fullmatch(r"\[(.+?)\](\d*)([+-]*)", s)
    tail = None
    if m:
        body, digits, signs = m.group(1), m.group(2), m.group(3)
        charge = None
        if signs:
            n = int(digits) if digits else len(signs)
            charge = n * (1 if signs[0] == "+" else -1)
        elif digits:
            return None
    elif re.fullmatch(r"\d*M.*", s):
        body, charge = s, None
        tail = re.fullmatch(r"(\d*M[+-].*[A-Za-z])([+-]+)", s)        # 不带方括号时末尾的符号是电荷：M+Cl- / M+FA- / M-H-
        if tail:
            body, signs = tail.group(1), tail.group(2)
            charge = len(signs) * (1 if signs[0] == "+" else -1)
    else:
        return None
    head = re.match(r"(\d*)M", body)
    if not head or "/" in body:
        return None
    n_m = int(head.group(1)) if head.group(1) else 1
    rest = body[head.end():]
    terms = []
    pos = 0
    for t in _TERM.finditer(rest):
        if t.start() != pos:
            return None
        pos = t.end()
        count = int(t.group(2)) if t.group(2) else 1
        terms.append([t.group(1), count, _ALIAS.get(t.group(3), t.group(3))])
    if pos != len(rest):
        return None
    # 加酸再减 H = 加酸根：+HCOOH −H -> +HCOO，+CH3COOH −H -> +CH3COO
    # 负离子只写了加酸（[M+FA]- / [M+HCOOH]-）也是酸根加合物：照字面是 M + 酸 + 电子的自由基负离子，ESI 里不出现；
    # GNPS 的 [M+FA]- 记录（2026-03 版 34,611 条）前体全是 M + 44.9982（甲酸根），不是 M + 46.0060
    for acid, anion in (("HCOOH", "HCOO"), ("CH3COOH", "CH3COO")):
        plus = [t for t in terms if t[0] == "+" and t[2] == acid]
        minus = [t for t in terms if t[0] == "-" and t[2] == "H"]
        if plus and minus and plus[0][1] == 1 and minus[0][1] == 1:
            terms.remove(plus[0])
            terms.remove(minus[0])
            terms.append(["+", 1, anion])
        elif plus and not minus and plus[0][1] == 1 and charge == -1 and _charge_of_terms([t for t in terms if t is not plus[0]]) == 0:
            terms.remove(plus[0])
            terms.append(["+", 1, anion])
    # 同名项合并，加项在前、减项在后，各自按名字排序
    merged = {}
    for sign, count, name in terms:
        merged[(sign, name)] = merged.get((sign, name), 0) + count
    ordered = sorted(merged.items(), key=lambda kv: (kv[0][0] == "-", kv[0][1]))
    implied = _charge_of_terms([(sg, c, nm) for (sg, nm), c in ordered])
    if charge is None:
        charge = implied
        if charge == 0 and ordered:
            return None                                       # 没写电荷又推不出来
    elif tail and implied != charge:
        return None                                           # 末尾符号和组成对不上（GNPS 的 M-2H- 其实是 [M-2H]2-，M-Ac-H- 是丢乙酸）：不猜
    text_terms = "".join(f"{sg}{c if c > 1 else ''}{nm}" for (sg, nm), c in ordered)
    mag = abs(charge)
    return f"[{n_m if n_m > 1 else ''}M{text_terms}]{mag if mag > 1 else ''}{'+' if charge > 0 else '-'}"


def adduct_charge_sign(adduct):
    """规范加合物的电荷符号：+1 / -1 / 0（不知道）。"""
    if not adduct:
        return 0
    return 1 if adduct.endswith("+") else -1 if adduct.endswith("-") else 0


def _ion_mode(text, adduct=None):
    t = str(text or "").strip().lower()
    if t.startswith("pos"):
        return 1
    if t.startswith("neg"):
        return -1
    return adduct_charge_sign(adduct)


def _num(text):
    """正的有限数；"449.2 1234"（PEPMASS 后面跟强度）只取第一个。"""
    try:
        v = float(str(text).split()[0])
    except (TypeError, ValueError, IndexError):
        return None
    return v if np.isfinite(v) and v > 0 else None


def _clean_peaks(mz, intensity):
    """去掉非法 / 零强度峰，按 m/z 升序，强度除以基峰。没有有效峰返回 None。"""
    mz = np.asarray(mz, dtype=np.float64)
    intensity = np.asarray(intensity, dtype=np.float64)
    ok = np.isfinite(mz) & np.isfinite(intensity) & (mz > 0) & (intensity > 0)
    if not ok.any():
        return None
    mz, intensity = mz[ok], intensity[ok]
    order = np.argsort(mz, kind="stable")
    mz, intensity = mz[order], intensity[order]
    return mz, (intensity / intensity.max()).astype(np.float32)


def _peaks_from_lines(lines):
    """峰行（"m/z 强度" 以空白分隔）-> (mz, intensity)。"""
    if not lines:
        return None
    try:
        flat = np.fromstring(" ".join(lines), sep=" ")
        if flat.size % 2 or flat.size != 2 * len(lines):
            raise ValueError
        arr = flat.reshape(-1, 2)
    except (ValueError, DeprecationWarning):
        rows = []
        for ln in lines:
            tok = ln.replace(",", " ").split()
            try:
                rows.append((float(tok[0]), float(tok[1])))
            except (IndexError, ValueError):
                continue
        if not rows:
            return None
        arr = np.array(rows, dtype=np.float64)
    return _clean_peaks(arr[:, 0], arr[:, 1])


def _open(path):
    return open(path, "r", encoding="utf-8", errors="replace", newline=None)


# ---------------------------------------------------------------- MassBank（NIST 格式 MSP）
def iter_massbank_msp(path):
    """MassBank_NISTformat.msp -> {字段名: 值, "_peaks": [行, ...]}。空行结束一条记录；Synon 行不要。"""
    rec, peaks, n_left = {}, [], 0
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if n_left > 0:
                if line.strip():
                    peaks.append(line)
                    n_left -= 1
                    continue
                n_left = 0
            if not line.strip():
                if rec:
                    rec["_peaks"] = peaks
                    yield rec
                rec, peaks = {}, []
                continue
            key, sep, value = line.partition(":")
            if not sep:
                continue
            key = key.strip()
            if key == "Num Peaks":
                try:
                    n_left = int(value.strip())
                except ValueError:
                    n_left = 0
            elif key != "Synon" and key not in rec:
                rec[key] = value.strip()
    if rec:
        rec["_peaks"] = peaks
        yield rec


def _is_ei(instrument_type):
    """GC-EI-TOF、EI-B 这类仪器类型（按 "-" 分段，含 EI 一段）。"""
    return "EI" in re.split(r"[-_/ ]", str(instrument_type or "").upper())


def normalise_massbank(rec):
    """MassBank 记录 -> (record, None) | (None, 跳过原因)。"""
    level = (rec.get("Spectrum_type") or "").strip().upper()
    itype = rec.get("Instrument_type", "")
    adduct_raw = rec.get("Precursor_type", "")
    adduct = normalize_adduct(adduct_raw)
    if level == "MS2":
        precursor = _num(rec.get("PrecursorMZ"))
        if precursor is None:
            return None, "MS/MS without a precursor m/z"
        kind = MSMS
    elif level == "MS" and _is_ei(itype):
        precursor, kind = None, EI
    else:
        return None, "not MS/MS or EI (MS1 of other ion sources, MS3 and above)"
    peaks = _peaks_from_lines(rec.get("_peaks"))
    if peaks is None:
        return None, "no valid peaks"
    ce = rec.get("Collision_energy", "")
    return {"accession": rec.get("DB#", ""), "name": rec.get("Name", ""), "inchikey": rec.get("InChIKey", ""),
            "smiles": rec.get("SMILES", ""), "adduct": adduct, "adduct_raw": adduct_raw,
            "precursor_mz": precursor, "ion_mode": _ion_mode(rec.get("Ion_mode"), adduct), "kind": kind,
            "instrument": " · ".join(x for x in (itype, rec.get("Instrument", "")) if x),
            "collision_energy": ce, "copy_of": None, "mz": peaks[0], "intensity": peaks[1]}, None


# ---------------------------------------------------------------- GNPS（MGF）
def iter_gnps_mgf(path):
    """ALL_GNPS_NO_PROPOGATED.mgf -> {KEY: value, "_peaks": [行, ...]}。"""
    rec, peaks, inside = None, [], False
    with _open(path) as fh:
        for line in fh:
            line = line.strip()
            if not inside:
                if line == "BEGIN IONS":
                    rec, peaks, inside = {}, [], True
                continue
            if line == "END IONS":
                rec["_peaks"] = peaks
                yield rec
                inside = False
                continue
            if not line:
                continue
            c = line[0]
            if c.isdigit() or c == "." or c == "-":
                peaks.append(line)
            else:
                key, sep, value = line.partition("=")
                if sep and key not in rec:
                    rec[key.strip().upper()] = value.strip()


def _gnps_na(value):
    v = str(value or "").strip()
    return "" if v.upper() in ("N/A", "NA", "NONE", "NULL") else v


_INCHIKEY = re.compile(r"[A-Z]{14}-[A-Z]{10}-[A-Z]")
# GNPS 里转存的 MassBank 记录：NAME 以 "Massbank:LU119903 ..." / "MassbankEU:SM866601 ..." 开头，ORGANISM 是 MASSBANK / MASSBANKEU
_MB_COPY_NAME = re.compile(r"^\s*MassBank(?:EU)?\s*:\s*([^\s|]+)", re.I)
MASSBANK_ORGANISMS = ("MASSBANK", "MASSBANKEU")
# 名字末尾的加合物："Caffeine M+H"、"... [M-H]-"
_NAME_ADDUCT = re.compile(r"\s(\[?\d*M(?:[+-]\d*[A-Za-z][A-Za-z0-9()]*)*\]?\d*[+-]*)\s*$")


def normalise_gnps(rec):
    """GNPS 记录 -> (record, None) | (None, 跳过原因)。"""
    level = str(rec.get("MSLEVEL", "")).strip()
    instrument = _gnps_na(rec.get("SOURCE_INSTRUMENT"))
    adduct_raw = _gnps_na(rec.get("ADDUCT")) or _gnps_na(rec.get("PRECURSOR_TYPE"))
    name = _gnps_na(rec.get("NAME"))
    if not adduct_raw:
        m = _NAME_ADDUCT.search(name)
        if m and normalize_adduct(m.group(1)):
            adduct_raw = m.group(1)
    adduct = normalize_adduct(adduct_raw)
    if level == "2":
        precursor = _num(rec.get("PEPMASS"))
        if precursor is None:
            return None, "MS/MS without a precursor m/z"
        kind = MSMS
    elif level == "1" and _is_ei(instrument):
        precursor, kind = None, EI
    else:
        return None, "not MS/MS or EI"
    peaks = _peaks_from_lines(rec.get("_peaks"))
    if peaks is None:
        return None, "no valid peaks"
    # GNPS 的 InChIKey 有时放在 INCHI / INCHIAUX 字段里
    inchikey = next((v for v in (_gnps_na(rec.get(k)) for k in ("INCHIKEY", "INCHI", "INCHIAUX")) if _INCHIKEY.fullmatch(v)), "")
    m = _MB_COPY_NAME.match(name)
    is_copy = bool(m) or str(rec.get("ORGANISM", "")).strip().upper() in MASSBANK_ORGANISMS
    return {"accession": rec.get("SPECTRUMID", ""), "name": name, "inchikey": inchikey,
            "smiles": _gnps_na(rec.get("SMILES")), "adduct": adduct, "adduct_raw": adduct_raw,
            "precursor_mz": precursor, "ion_mode": _ion_mode(rec.get("IONMODE"), adduct), "kind": kind,
            "instrument": instrument, "collision_energy": _gnps_na(rec.get("COLLISION_ENERGY")),
            "copy_of": (m.group(1) if m else "") if is_copy else None, "mz": peaks[0], "intensity": peaks[1]}, None


SOURCES = {
    "MassBank": (iter_massbank_msp, normalise_massbank),
    "GNPS": (iter_gnps_mgf, normalise_gnps),
}


def iter_records(library, path):
    """流式产出 (record, skip_reason)。"""
    reader, normaliser = SOURCES[library]
    for rec in reader(Path(path)):
        yield normaliser(rec)
