"""MS2 镜像图（Plotly）：上半为 query，下半为 reference；深色 = 按 CosineGreedy 规则匹配上的峰。峰标签和悬停卡片由前端 Figure.tsx 画（和其它谱图一样）。"""
from __future__ import annotations

import numpy as np
import plotly.graph_objects as go

from mzdesk.shared import colors as theme
from mzdesk.shared import ms2_plot as base

Y_TOP = 118.0
COMPACT_BELOW = 440          # px: a mirror plot lower than this gets the short y-axis titles ("Experimental · rel. %")


def _as_peaks(peaks):
    array = np.asarray(peaks, dtype=float)
    if array.size == 0:
        return np.array([]), np.array([])
    if array.ndim != 2 or array.shape[1] < 2:
        raise ValueError("peaks must be a sequence of (mz, intensity) pairs")

    mz, intensity = array[:, 0], array[:, 1]
    valid = np.isfinite(mz) & np.isfinite(intensity) & (mz > 0) & (intensity > 0)
    mz, intensity = mz[valid], intensity[valid]
    if mz.size == 0:
        return np.array([]), np.array([])

    order = np.argsort(mz)
    mz, intensity = mz[order], intensity[order]
    intensity = intensity / intensity.max() * 100.0
    return mz, intensity


def cosine_greedy_matches(mz1, rel1, mz2, rel2, tolerance):
    """Match peaks 1:1 like matchms CosineGreedy with default parameters."""
    if tolerance is None or mz1.size == 0 or mz2.size == 0:
        return []

    candidate_i, candidate_j = [], []
    lowest = 0
    for i, mz in enumerate(mz1):
        low, high = mz - tolerance, mz + tolerance
        for j in range(lowest, mz2.size):
            if mz2[j] > high:
                break
            if mz2[j] < low:
                lowest = j + 1
            else:
                candidate_i.append(i)
                candidate_j.append(j)
    if not candidate_i:
        return []

    candidate_i = np.asarray(candidate_i)
    candidate_j = np.asarray(candidate_j)
    products = (rel1[candidate_i] / 100.0) * (rel2[candidate_j] / 100.0)
    order = np.argsort(products, kind="mergesort")[::-1]

    used1, used2, matches = set(), set(), []
    for index in order:
        i, j = int(candidate_i[index]), int(candidate_j[index])
        if i in used1 or j in used2:
            continue
        used1.add(i)
        used2.add(j)
        matches.append((i, j))
    return matches


def match_summary(peaks_a, peaks_b, tolerance=0.02):
    mz1, rel1 = _as_peaks(peaks_a)
    mz2, rel2 = _as_peaks(peaks_b)
    matches = cosine_greedy_matches(mz1, rel1, mz2, rel2, tolerance)

    query_total = len(mz1)
    ref_total = len(mz2)
    query_matched = len({i for i, _ in matches})
    ref_matched = len({j for _, j in matches})
    return {
        "matched_pairs": len(matches),
        "query_total": query_total,
        "query_matched": query_matched,
        "query_matched_pct": query_matched / query_total * 100.0 if query_total else 0.0,
        "ref_total": ref_total,
        "ref_matched": ref_matched,
        "ref_matched_pct": ref_matched / ref_total * 100.0 if ref_total else 0.0,
    }


def _valid_precursor(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) and value > 0 else None


def _card(label, mass, intensity, matched, precursor_mz, other_label, other_mass, sign):
    """Hover card of one peak (Figure.tsx PeakCard); ``rel`` is the signed height, so the reference half is labelled below its tips."""
    rows = [["Rel. intensity", f"{base.fmt_rel(intensity)} %"]]
    if other_mass is not None:
        rows += [[f"{other_label} m/z", base.fmt_mz(other_mass)], ["Δ m/z", f"{other_mass - mass:+.4g} Da"]]
    if precursor_mz is not None and precursor_mz - mass > 0:
        rows.append(["Precursor − peak", f"{precursor_mz - mass:.6g}"])
    if float(f"{mass:.4f}") != mass:
        rows.append(["m/z (more precision)", base.fmt_mz(mass)])
    return {"mz": float(mass), "rel": float(sign * intensity), "on": bool(matched),
            "kind": f"{label} · {'matched' if matched else 'unmatched'}", "rows": rows}


def _precursor_text(value, other):
    """Four decimals, unless that would make two different precursors read the same."""
    short = f"{value:.4f}"
    return base.fmt_mz(value) if other is not None and other != value and f"{other:.4f}" == short else short


def build_mirror_figure(peaks_a, peaks_b, precursor_mz1=None, precursor_mz2=None,
                        tolerance=0.02, labels=("Query", "Reference"), height=520):
    """Mirror plot -> ``(fig, cards)``; pass cards to ``common.plotly(fig, peaks=cards)``: Figure.tsx labels the peaks
    and shows the hover card. Query up, reference down; dark = matched (CosineGreedy), light = unmatched."""
    mz1, rel1 = _as_peaks(peaks_a)
    mz2, rel2 = _as_peaks(peaks_b)
    label_a, label_b = labels
    precursor_a = _valid_precursor(precursor_mz1)
    precursor_b = _valid_precursor(precursor_mz2)

    if mz1.size == 0 and mz2.size == 0:
        return go.Figure(), []

    matches = cosine_greedy_matches(mz1, rel1, mz2, rel2, tolerance)
    partner_a = {i: j for i, j in matches}
    partner_b = {j: i for i, j in matches}

    x_range = base.mz_range(np.concatenate([mz1, mz2, [p for p in (precursor_a, precursor_b) if p is not None]]))
    width = base.stick_px(np.concatenate([mz1, mz2]), x_range)

    figure = go.Figure()
    annotations, cards = [], []
    halves = (
        (mz1, rel1, +1, label_a, precursor_a, partner_a, mz2, label_b),
        (mz2, rel2, -1, label_b, precursor_b, partner_b, mz1, label_a),
    )
    for mz, rel, sign, label, precursor, partner, other_mz, other_label in halves:
        if mz.size == 0:
            continue
        heights = sign * rel
        is_matched = np.array([index in partner for index in range(mz.size)], dtype=bool)
        for mask, color, suffix in ((~is_matched, base.LIGHT, "unmatched"), (is_matched, base.DARK, "matched")):
            if mask.any():
                figure.add_trace(base.stick_trace(mz[mask], heights[mask], color, f"{label} · {suffix}", width=width))
        # one hover trace per half, both named "peaks": Figure.tsx picks the half the pointer is in
        figure.add_trace(base.hover_markers(mz, heights, customdata=list(range(len(cards), len(cards) + mz.size))))
        cards += [_card(label, mz[index], rel[index], is_matched[index], precursor, other_label,
                        other_mz[partner[index]] if index in partner else None, sign) for index in range(mz.size)]

    shapes = []
    same_precursor = precursor_a is not None and precursor_a == precursor_b
    for precursor, other, sign, text in (
            (precursor_a, precursor_b, +1, "precursor" if same_precursor else f"{label_a} precursor"),
            (None if same_precursor else precursor_b, precursor_a, -1, f"{label_b} precursor")):
        if precursor is None:
            continue
        y0, y1 = (-100.0, 100.0) if same_precursor else (0.0, sign * 100.0)
        shapes.append(base.precursor_line(precursor, y0, y1))
        annotations.append(base.precursor_label(
            precursor, f"{text} {_precursor_text(precursor, other)}", x_range, y=sign * Y_TOP, yref="y",
            yanchor="top" if sign > 0 else "bottom", yshift=0))

    half_label = dict(
        xref="paper", yref="paper", x=0, xanchor="right", xshift=-44,
        yanchor="middle", textangle=-90, showarrow=False,
        font=dict(family=base.FONT, color=theme.MUTED, size=12),
    )
    unit = "rel. intensity (%)" if height >= COMPACT_BELOW else "rel. %"      # each half of a low plot is too short for the long title
    annotations += [
        dict(half_label, text=f"{label_a} · {unit}", y=0.75),
        dict(half_label, text=f"{label_b} · {unit}", y=0.25),
    ]

    ticks = [-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100]
    figure.update_layout(**base.base_layout(
        height, x_range=x_range, y_range=[-Y_TOP, Y_TOP], y_title=None, shapes=shapes, annotations=annotations,
        margin=dict(l=80, r=22, t=16, b=50),
        yaxis=dict(tickvals=ticks, ticktext=[str(abs(tick)) for tick in ticks],
                   zeroline=True, zerolinecolor=theme.AXIS, zerolinewidth=1)))
    return figure, cards


if __name__ == "__main__":
    demo_a = [(69.0335, 12), (85.0289, 100), (110.0340, 18), (141.0650, 8)]
    demo_b = [(85.0289, 95), (110.0340, 20), (128.0572, 6), (200.1000, 12)]
    figure, _cards = build_mirror_figure(
        demo_a, demo_b, precursor_mz1=250.12, precursor_mz2=250.12, tolerance=0.02,
    )
    figure.write_html("/tmp/mirror_demo.html", include_plotlyjs="cdn")
    print("wrote /tmp/mirror_demo.html")
