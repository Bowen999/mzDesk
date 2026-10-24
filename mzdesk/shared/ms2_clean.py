"""MS2 峰表清洗：噪声阈值、去前体离子、归一化、排序。每个输入峰都有“保留 / 去掉及原因”的记录。"""
import math

NORMALIZE = ["None", "Base peak = 100", "Base peak = 999", "Base peak = 1", "Sum = 1"]
SORT_BY = ["m/z", "Intensity"]
_BASE = {"Base peak = 100": 100.0, "Base peak = 999": 999.0, "Base peak = 1": 1.0}


def clean_spectrum(peaks, *, noise_pct=0.0, normalize="None", precursor_mz=None, remove_precursor=False,
                   precursor_tol=0.02, sort_by="m/z"):
    """Clean ``peaks`` [(mz, intensity)]. Returns ``{"peaks", "rows", "n_in", "n_out"}``.

    Invalid peaks (non-finite or <= 0) are dropped and duplicate m/z (5 dp) collapsed. A peak is
    removed as "precursor" when within ± ``precursor_tol`` Da of ``precursor_mz`` (only if
    ``remove_precursor``), else as "noise" when below ``noise_pct`` % of the input base peak.
    ``normalize`` is applied to the kept peaks; ``rows`` (m/z order) hold every input peak with
    ``mz``, ``intensity`` (same scale as the kept peaks), ``rel`` (% of the input base peak),
    ``kept`` and ``reason`` ("noise", "precursor" or ""). ``peaks`` is ordered by ``sort_by``.
    """
    if normalize not in NORMALIZE:
        raise ValueError(f"normalize must be one of {', '.join(NORMALIZE)}")
    if sort_by not in SORT_BY:
        raise ValueError(f"sort_by must be one of {', '.join(SORT_BY)}")
    agg = {}
    for mz, inten in peaks:
        if all(isinstance(v, (int, float)) and math.isfinite(v) for v in (mz, inten)) and mz > 0 and inten > 0:
            key = round(float(mz), 5)
            agg[key] = max(agg.get(key, 0.0), float(inten))
    base = max(agg.values(), default=1.0)
    cut = noise_pct / 100.0 * base
    rows = []
    for mz, inten in sorted(agg.items()):
        if remove_precursor and precursor_mz and abs(mz - precursor_mz) <= precursor_tol:
            reason = "precursor"
        else:
            reason = "noise" if inten < cut else ""
        rows.append({"mz": mz, "intensity": inten, "rel": inten / base * 100.0, "kept": not reason, "reason": reason})
    kept = [r["intensity"] for r in rows if r["kept"]]
    if normalize in _BASE:
        scale = _BASE[normalize] / max(kept, default=1.0)
    elif normalize == "Sum = 1":
        scale = 1.0 / (sum(kept) or 1.0)
    else:
        scale = 1.0
    for r in rows:
        r["intensity"] *= scale
    out = [(r["mz"], r["intensity"]) for r in rows if r["kept"]]
    out.sort(key=(lambda p: p[0]) if sort_by == "m/z" else (lambda p: -p[1]))
    return {"peaks": out, "rows": rows, "n_in": len(rows), "n_out": len(out)}
