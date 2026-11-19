"""同位素分布引擎（Theoretical 页和 Mass search 页共用）：中性分子式 + 加合物 -> 离子式 -> IsoSpecPy 精细结构 -> 按分辨率合并 -> 峰表；
以及观测同位素峰与理论峰的余弦相似度。

约定：
  - "M" = 每个元素取丰度最高的同位素（11B、79Br、35Cl …）的离子，与 msbuddy / PubChem 的单同位素质量一致；比 M 轻的峰标 "M−1"（如 10B）。
  - m/z = (质量 − z·电子) / |z|，电子质量与 msbuddy_engine 相同；加合物定义直接用 msbuddy 的（[M]+ / [M]- 单独处理）。
  - 元素不限于 msbuddy 的 12 个：IsoSpecPy 认识的都行；带同位素标记（D、[13C]）或未知元素 -> ValueError。
  - IsoSpecPy 只在用到时才导入。
"""
import re
from functools import lru_cache

RESOLUTIONS = {"Orbitrap 60k": 60000.0, "Q-TOF 30k": 30000.0, "Unit mass": None}      # None = 按名义质量合并
DEFAULT_RESOLUTION = "Orbitrap 60k"
UNIT_TOL_DA = 0.5               # Unit mass 时匹配峰的容差（Da）
MIN_REL = 0.1                   # 保留 >= 最强峰的 0.1 %
PROB_COVER = 0.99999            # IsoSpecPy 覆盖的总概率
MAX_ATOMS = 3000                # 再大的分子式不算（精细结构会爆）
_ELECTRON = 0.00054857990946    # 与 msbuddy_engine._ELECTRON 一致
ISOTOPE_LABELS = {"D", "T"}      # IsoSpecPy 的表里有 D / T，这里不当元素
_MINUS = "−"


def _tables():
    from IsoSpecPy import PeriodicTbl
    return PeriodicTbl.symbol_to_masses, PeriodicTbl.symbol_to_probs


@lru_cache(maxsize=None)
def _mono_mass(element):
    """丰度最高的同位素的质量（B -> 11B）。"""
    masses, probs = _tables()
    return masses[element][max(range(len(probs[element])), key=probs[element].__getitem__)]


def parse_formula(formula):
    """'C8H10N4O2' / 'C5H14NO+' -> ({元素: 个数}, 公式自带的电荷符号 0/+1/-1)；不认识 -> ValueError。"""
    s = (formula or "").strip()
    charge = 0
    m = re.search(r"([+-])(\d*)$", s)
    if m:
        charge = (1 if m.group(1) == "+" else -1) * int(m.group(2) or 1)
        s = s[:m.start()]
    if not s:
        raise ValueError("Empty formula.")
    if not re.fullmatch(r"(?:[A-Z][a-z]?\d*)+", s):
        raise ValueError(f"Cannot read the formula {formula!r}. Use plain element symbols and counts, e.g. C8H10N4O2 (no isotope labels).")
    masses, _ = _tables()
    counts = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", s):
        if el not in masses or el in ISOTOPE_LABELS:
            raise ValueError(f"Unknown element or isotope label {el!r} in {formula!r}. Isotope labels such as D are not supported.")
        counts[el] = counts.get(el, 0) + int(n or 1)
    if sum(counts.values()) > MAX_ATOMS:
        raise ValueError(f"{formula!r} has too many atoms for an isotope pattern (limit {MAX_ATOMS}).")
    return counts, charge


def ion_composition(formula, adduct):
    """(离子的 {元素: 个数}, 电荷 z)；加合物不适用（要去掉分子式里没有的原子）或认不出 -> ValueError。"""
    from mzdesk.ms1 import msbuddy_engine as bd
    counts, _ = parse_formula(formula)
    if adduct in ("[M]+", "[M]-"):
        return counts, (1 if adduct == "[M]+" else -1)
    ion = bd._parse_adduct(adduct)
    if ion is None:
        raise ValueError(f"Unknown adduct {adduct}.")
    ion_counts = {el: n * int(ion.m) for el, n in counts.items()}
    for el, n in zip(bd._ALPHABET, ion.net_formula.array):
        if n:
            ion_counts[el] = ion_counts.get(el, 0) + int(n)
    if any(n < 0 for n in ion_counts.values()):
        raise ValueError(f"{formula} cannot form {adduct} (the adduct removes atoms the formula does not have).")
    return {el: n for el, n in ion_counts.items() if n}, int(ion.charge)


def _formula_string(counts):
    """Hill 顺序（C、H 在前，其余按字母），带电荷符号放在外面。"""
    order = [e for e in ("C", "H") if e in counts] + sorted(e for e in counts if e not in ("C", "H"))
    return "".join(f"{e}{counts[e] if counts[e] != 1 else ''}" for e in order)


def ion_mz(formula, adduct):
    """单同位素（M）离子的 (m/z, 带电荷符号的离子式)；认不出 -> ValueError。"""
    counts, z = ion_composition(formula, adduct)
    mass = sum(_mono_mass(el) * n for el, n in counts.items())
    return (mass - z * _ELECTRON) / abs(z), _formula_string(counts) + ("+" if z > 0 else "-")


def _label(offset):
    return "M" if offset == 0 else f"M+{offset}" if offset > 0 else f"M{_MINUS}{-offset}"


def _merge(fine, resolution):
    """fine: [(m/z, 概率)] 按 m/z 升序 -> 合并后的 [(m/z, 概率)]。R：与簇质心相差不到 FWHM = m/z / R 的并入该簇（贪心）；Unit mass 在调用处按名义质量分簇。"""
    fwhm_over = 1.0 / resolution
    clusters = []       # [sum(p), sum(p*mz)]
    for mz, p in fine:
        if clusters:
            c = clusters[-1]
            if mz - c[1] / c[0] < mz * fwhm_over:
                c[0] += p
                c[1] += p * mz
                continue
        clusters.append([p, p * mz])
    return [(c[1] / c[0], c[0]) for c in clusters]


@lru_cache(maxsize=256)
def _pattern(formula, adduct, resolution):
    import IsoSpecPy
    if resolution not in RESOLUTIONS:
        raise ValueError(f"Unknown resolution {resolution!r}.")
    counts, z = ion_composition(formula, adduct)
    mono_mz, ion_formula = ion_mz(formula, adduct)
    iso = IsoSpecPy.IsoTotalProb(formula=_formula_string(counts), prob_to_cover=PROB_COVER, get_confs=False)
    zabs = abs(z)
    fine = sorted(((m - z * _ELECTRON) / zabs, p) for m, p in zip(iso.masses, iso.probs))
    R = RESOLUTIONS[resolution]
    if R is None:
        by_offset = {}
        for mz, p in fine:
            c = by_offset.setdefault(round((mz - mono_mz) * zabs), [0.0, 0.0])
            c[0] += p
            c[1] += p * mz
        merged = sorted((c[1] / c[0], c[0]) for c in by_offset.values())
    else:
        merged = _merge(fine, R)
    top = max(p for _, p in merged)
    merged = [(mz, p / top * 100.0) for mz, p in merged if p / top * 100.0 >= MIN_REL]
    offsets = [round((mz - mono_mz) * zabs) for mz, _ in merged]
    seen = {}
    for o in offsets:
        seen[o] = seen.get(o, 0) + 1
    used, peaks = {}, []
    for (mz, ab), o in zip(merged, offsets):
        label = _label(o)
        if seen[o] > 1:
            used[o] = used.get(o, 0) + 1
            label += f" ({chr(96 + used[o])})"
        peaks.append({"label": label, "offset": o, "mz": mz, "abundance": ab})
    return {"mono_mz": mono_mz, "ion_formula": ion_formula, "charge": z, "resolution": resolution, "peaks": peaks}


def pattern(formula, adduct, resolution=DEFAULT_RESOLUTION):
    """{"mono_mz", "ion_formula", "charge", "resolution", "peaks": [{"label", "offset", "mz", "abundance"(最强 = 100)}]}，按 m/z 升序。

    同一名义质量有几个分开的峰时，标签带 (a) (b)。ValueError：分子式 / 加合物 / 分辨率不能算。
    """
    res = _pattern((formula or "").strip(), adduct, resolution)
    return {**res, "peaks": [dict(p) for p in res["peaks"]]}


def match_tolerance(resolution, tol_ppm):
    """匹配用的容差：(ppm 或 Da, 是否 Da)。"""
    return (UNIT_TOL_DA, True) if RESOLUTIONS[resolution] is None else (tol_ppm, False)


def similarity(peaks, observed, resolution=DEFAULT_RESOLUTION, tol_ppm=5.0):
    """理论峰 vs 观测 [(mz, 强度)] 的余弦相似度（0-1）。

    每个理论峰（丰度高的先挑）匹配容差内最强且未被用过的观测峰；没匹配上的理论峰取 0；观测里没被匹配的峰忽略（可能是别的离子）。
    返回 {"score", "rows": [{"label", "mz", "abundance", "obs_mz", "obs_rel", "ppm"}]}；obs_rel 是匹配峰相对最强匹配峰的 %。
    """
    tol, in_da = match_tolerance(resolution, tol_ppm)
    used, hit = set(), {}
    for i in sorted(range(len(peaks)), key=lambda k: -peaks[k]["abundance"]):
        t = peaks[i]["mz"]
        window = tol if in_da else t * tol * 1e-6
        cand = [j for j, (mz, inten) in enumerate(observed) if j not in used and inten > 0 and abs(mz - t) <= window]
        if cand:
            j = max(cand, key=lambda k: observed[k][1])
            used.add(j)
            hit[i] = j
    obs = [observed[hit[i]][1] if i in hit else 0.0 for i in range(len(peaks))]
    theo = [p["abundance"] for p in peaks]
    dot = sum(a * b for a, b in zip(theo, obs))
    norm = (sum(a * a for a in theo) * sum(b * b for b in obs)) ** 0.5
    score = dot / norm if norm else 0.0
    top = max(obs) if obs and max(obs) > 0 else 1.0
    rows = []
    for i, p in enumerate(peaks):
        j = hit.get(i)
        rows.append({"label": p["label"], "mz": p["mz"], "abundance": p["abundance"],
                     "obs_mz": observed[j][0] if j is not None else None,
                     "obs_rel": obs[i] / top * 100.0 if j is not None else None,
                     "ppm": (observed[j][0] - p["mz"]) / p["mz"] * 1e6 if j is not None else None})
    return {"score": min(1.0, score), "rows": rows}
