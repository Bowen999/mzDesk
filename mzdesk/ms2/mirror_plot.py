"""MS2 镜像图（Plotly）：上半为 query，下半为 reference；深色 = 按 CosineGreedy 规则匹配上的峰。"""
from __future__ import annotations

import numpy as np
import plotly.graph_objects as go

from mzdesk.shared import colors as theme

MATCH_COLOR = theme.PLOT_DARK
UNMATCHED_COLOR = theme.PLOT_LIGHT
GRID = theme.SURFACE_2
FONT = theme.PLOT_FONT
MONO = theme.PLOT_MONO

STICK_PX = 2
MIN_STICK = 1.6
Y_TOP = 118.0
LABELS_PER_HALF = 8
LABEL_MIN_REL = 5.0
LABEL_SEP = 0.085


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


def _sticks(mz, heights, color, name):
    xs, ys = [], []
    for mass, height in zip(mz, heights):
        xs += [float(mass), float(mass), None]
        ys += [0.0, float(height), None]
    return go.Scatter(
        x=xs, y=ys, mode="lines", line=dict(color=color, width=STICK_PX),
        hoverinfo="skip", name=name, showlegend=False,
    )


def _tip_markers(mz, heights, rel, label, is_matched, precursor_mz, partners):
    customdata = []
    for mass, intensity, matched, partner in zip(mz, rel, is_matched, partners):
        loss = ""
        if precursor_mz is not None and precursor_mz - mass > 0:
            loss = f"<br>Neutral loss {precursor_mz - mass:.4f} Da"
        title = f"{label} · {'matched' if matched else 'unmatched'}"
        customdata.append([title, float(mass), round(float(intensity), 1), loss, partner])
    return go.Scatter(
        x=mz,
        y=heights,
        mode="markers",
        marker=dict(size=10, color="rgba(0,0,0,0)"),
        customdata=customdata,
        hovertemplate=(
            "<b>%{customdata[0]}</b><br>"
            "m/z %{customdata[1]:.4f}<br>"
            "Rel. intensity %{customdata[2]:.1f}%"
            "%{customdata[3]}%{customdata[4]}"
            "<extra></extra>"
        ),
        name=f"{label} peaks",
        showlegend=False,
    )


def _label_indices(mz, rel, x_span):
    min_separation = LABEL_SEP * x_span
    chosen = []
    for index in np.argsort(-rel, kind="stable"):
        if rel[index] < LABEL_MIN_REL or len(chosen) >= LABELS_PER_HALF:
            break
        if all(abs(mz[index] - mz[selected]) >= min_separation for selected in chosen):
            chosen.append(int(index))
    return chosen


def _anchor(x, x_range):
    fraction = (x - x_range[0]) / (x_range[1] - x_range[0])
    return "right" if fraction > 0.8 else ("left" if fraction < 0.2 else "center")


def build_mirror_figure(peaks_a, peaks_b, precursor_mz1=None, precursor_mz2=None,
                        tolerance=0.02, labels=("Query", "Reference"), height=520):
    mz1, rel1 = _as_peaks(peaks_a)
    mz2, rel2 = _as_peaks(peaks_b)
    label_a, label_b = labels
    precursor_a = _valid_precursor(precursor_mz1)
    precursor_b = _valid_precursor(precursor_mz2)

    if mz1.size == 0 and mz2.size == 0:
        return go.Figure()

    matches = cosine_greedy_matches(mz1, rel1, mz2, rel2, tolerance)
    partner_a = {i: j for i, j in matches}
    partner_b = {j: i for i, j in matches}

    xs = np.concatenate([mz1, mz2, [p for p in (precursor_a, precursor_b) if p is not None]])
    lo, hi = float(xs.min()), float(xs.max())
    padding = max((hi - lo) * 0.05, 5.0)
    x_range = [max(0.0, lo - padding), hi + padding]
    x_span = x_range[1] - x_range[0]

    figure = go.Figure()
    annotations = []
    halves = (
        (mz1, rel1, +1, label_a, precursor_a, partner_a, mz2, label_b),
        (mz2, rel2, -1, label_b, precursor_b, partner_b, mz1, label_a),
    )
    for mz, rel, sign, label, precursor, partner, other_mz, other_label in halves:
        if mz.size == 0:
            continue
        heights = sign * np.maximum(rel, MIN_STICK)
        is_matched = np.array([index in partner for index in range(mz.size)], dtype=bool)
        partners = [
            f"<br>Matched ↔ {other_label} m/z {other_mz[partner[index]]:.4f} "
            f"(Δ {other_mz[partner[index]] - mz[index]:+.4f} Da)" if index in partner else ""
            for index in range(mz.size)
        ]
        for mask, color, suffix in ((~is_matched, UNMATCHED_COLOR, "unmatched"),
                                    (is_matched, MATCH_COLOR, "matched")):
            if mask.any():
                figure.add_trace(_sticks(mz[mask], heights[mask], color, f"{label} · {suffix}"))
        figure.add_trace(_tip_markers(mz, heights, rel, label, is_matched, precursor, partners))

        for index in _label_indices(mz, rel, x_span):
            annotations.append(dict(
                x=float(mz[index]), y=float(heights[index]), xref="x", yref="y",
                yanchor="bottom" if sign > 0 else "top", yshift=3 * sign, showarrow=False,
                text=f"{mz[index]:.4f}", font=dict(family=MONO, size=10, color=theme.INK),
            ))

    shapes = []
    same_precursor = (precursor_a is not None and precursor_b is not None and
                      abs(precursor_a - precursor_b) < 1e-6)
    for precursor, sign, text in (
            (precursor_a, +1, "precursor" if same_precursor else f"{label_a} precursor"),
            (None if same_precursor else precursor_b, -1, f"{label_b} precursor")):
        if precursor is None:
            continue
        y0, y1 = (-100.0, 100.0) if same_precursor else (0.0, sign * 100.0)
        shapes.append(dict(
            type="line", xref="x", yref="y", x0=precursor, x1=precursor, y0=y0, y1=y1,
            layer="below", line=dict(color=theme.MUTED, width=1, dash="dot"),
        ))
        annotations.append(dict(
            x=precursor, y=sign * Y_TOP, xref="x", yref="y", showarrow=False,
            yanchor="top" if sign > 0 else "bottom", xanchor=_anchor(precursor, x_range),
            text=f"{text} {precursor:.4f}", font=dict(family=MONO, size=10, color=theme.MUTED),
        ))

    half_label = dict(
        xref="paper", yref="paper", x=0, xanchor="right", xshift=-44,
        yanchor="middle", textangle=-90, showarrow=False,
        font=dict(family=FONT, color=theme.MUTED, size=12),
    )
    annotations += [
        dict(half_label, text=f"{label_a} · rel. intensity (%)", y=0.75),
        dict(half_label, text=f"{label_b} · rel. intensity (%)", y=0.25),
    ]

    ticks = [-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100]
    axis = dict(
        showline=True, linecolor=theme.INK, linewidth=1, ticks="outside",
        tickcolor=theme.INK, ticklen=4, tickfont=dict(color=theme.MUTED, size=11),
    )
    figure.update_layout(
        template="simple_white",
        font=dict(family=FONT, color=theme.INK, size=12),
        plot_bgcolor=theme.SURFACE,
        paper_bgcolor=theme.SURFACE,
        height=height,
        margin=dict(l=80, r=22, t=16, b=50),
        showlegend=False,
        hovermode="x",
        hoverdistance=16,
        dragmode="zoom",
        hoverlabel=dict(
            bgcolor=theme.SURFACE, bordercolor=theme.LINE,
            font=dict(family=FONT, color=theme.INK, size=12), align="left",
        ),
        shapes=shapes,
        annotations=annotations,
        xaxis=dict(
            axis, range=x_range, hoverformat=".4f",
            title=dict(text="m/z", font=dict(color=theme.MUTED, size=12)),
            showgrid=False, zeroline=False,
        ),
        yaxis=dict(
            axis, range=[-Y_TOP, Y_TOP], tickvals=ticks, ticktext=[str(abs(tick)) for tick in ticks],
            showgrid=True, gridcolor=GRID,
            zeroline=True, zerolinecolor=theme.INK, zerolinewidth=1,
        ),
    )
    return figure


if __name__ == "__main__":
    demo_a = [(69.0335, 12), (85.0289, 100), (110.0340, 18), (141.0650, 8)]
    demo_b = [(85.0289, 95), (110.0340, 20), (128.0572, 6), (200.1000, 12)]
    figure = build_mirror_figure(
        demo_a, demo_b, precursor_mz1=250.12, precursor_mz2=250.12, tolerance=0.02,
    )
    figure.write_html("/tmp/mirror_demo.html", include_plotlyjs="cdn")
    print("wrote /tmp/mirror_demo.html")
