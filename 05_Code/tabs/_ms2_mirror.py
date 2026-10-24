"""Interactive mirror MS² plot for the MS² Similarity tab.

``build_mirror_figure`` returns a Plotly figure; ``render_mirror`` draws it with
``st.plotly_chart``. The query spectrum points up and the reference spectrum
points down.

Matched peaks come from ``cosine_greedy_matches``, which reproduces the 1:1 peak
assignment of matchms ``CosineGreedy`` (default parameters) exactly, so the
highlighted peaks and the "Matched pairs" count agree with the CosineGreedy
score card.

Visual language shared with the MS² Interpretation tab: peaks are 2-px sticks
(constant pixel width, so neighbouring peaks never merge when zooming and every
peak stays visible when zoomed out); matched peaks use ``theme.PLOT_DARK`` and
unmatched peaks ``theme.PLOT_LIGHT``; the strongest peaks carry m/z labels;
dotted lines mark the precursor m/z; hover snaps to the nearest peak along m/z.
"""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
import streamlit as st

try:
    import theme
except ImportError:  # pragma: no cover - standalone run from tabs/ (demo below)
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import theme


MATCH_COLOR = theme.PLOT_DARK       # matched within tolerance
UNMATCHED_COLOR = theme.PLOT_LIGHT  # no partner in the other spectrum
GRID = theme.SURFACE_2
FONT = theme.PLOT_FONT
MONO = theme.PLOT_MONO

STICK_PX = 2          # stick line width (px)
MIN_STICK = 1.6       # % - tiny peaks are drawn at least this tall so they stay visible
Y_TOP = 118.0         # y-range half height: headroom for m/z and precursor labels
LABELS_PER_HALF = 8   # at most this many m/z labels per spectrum
LABEL_MIN_REL = 5.0   # only peaks >= 5 % relative intensity get a label
LABEL_SEP = 0.085     # labels in one half stay >= 8.5 % of the m/z span apart

_PLOT_CONFIG = {
    "displaylogo": False,
    # Double-click always returns to the padded initial view (never Plotly's autorange).
    "doubleClick": "reset",
    "modeBarButtonsToRemove": ["select2d", "lasso2d", "autoScale2d"],
    "toImageButtonOptions": {"format": "png", "filename": "mirror_ms2", "scale": 2},
}


def _as_peaks(peaks):
    """Coerce ``[(mz, intensity), ...]`` into sorted m/z and relative intensity (0-100)."""
    arr = np.asarray(peaks, dtype=float)
    if arr.size == 0:
        return np.array([]), np.array([])
    if arr.ndim != 2 or arr.shape[1] < 2:
        raise ValueError("peaks must be a sequence of (mz, intensity) pairs")

    mz, intensity = arr[:, 0], arr[:, 1]
    valid = np.isfinite(mz) & np.isfinite(intensity) & (mz > 0) & (intensity > 0)
    mz, intensity = mz[valid], intensity[valid]
    if mz.size == 0:
        return np.array([]), np.array([])

    order = np.argsort(mz)
    mz, intensity = mz[order], intensity[order]
    intensity = intensity / intensity.max() * 100.0
    return mz, intensity


def cosine_greedy_matches(mz1, rel1, mz2, rel2, tolerance):
    """1:1 peak assignment identical to matchms ``CosineGreedy`` (mz_power=0, intensity_power=1).

    Every query/reference pair with ``|Δm/z| <= tolerance`` is a candidate. Candidates
    are ranked by intensity product (stable sort, highest first - ties are resolved
    exactly like matchms) and accepted greedily so each peak is used at most once.
    ``rel1``/``rel2`` are relative intensities (0-100) as returned by ``_as_peaks``;
    they are rescaled to matchms' base-peak = 1 normalisation so products compare
    bit-for-bit like matchms. Returns ``[(i, j), ...]`` index pairs.
    """
    if tolerance is None or mz1.size == 0 or mz2.size == 0:
        return []

    cand_i, cand_j = [], []
    lowest = 0
    for i, mz in enumerate(mz1):
        low, high = mz - tolerance, mz + tolerance
        for j in range(lowest, mz2.size):
            if mz2[j] > high:
                break
            if mz2[j] < low:
                lowest = j + 1
            else:
                cand_i.append(i)
                cand_j.append(j)
    if not cand_i:
        return []

    cand_i = np.asarray(cand_i)
    cand_j = np.asarray(cand_j)
    products = (rel1[cand_i] / 100.0) * (rel2[cand_j] / 100.0)
    order = np.argsort(products, kind="mergesort")[::-1]

    used1, used2, matches = set(), set(), []
    for k in order:
        i, j = int(cand_i[k]), int(cand_j[k])
        if i in used1 or j in used2:
            continue
        used1.add(i)
        used2.add(j)
        matches.append((i, j))
    return matches


def match_summary(peaks_a, peaks_b, tolerance=0.02):
    """Matched-peak counts and percentages (CosineGreedy assignment) for both spectra."""
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
    """One trace of vertical sticks from 0 to ``heights`` (hover handled by tip markers)."""
    xs, ys = [], []
    for m, h in zip(mz, heights):
        xs += [float(m), float(m), None]
        ys += [0.0, float(h), None]
    return go.Scatter(x=xs, y=ys, mode="lines", line=dict(color=color, width=STICK_PX),
                      hoverinfo="skip", name=name, showlegend=False)


def _tip_markers(mz, heights, rel, label, is_matched, precursor_mz, partners):
    """Invisible markers at the stick tips; hovermode "x" snaps to the nearest one."""
    customdata = []
    for m, r, matched, partner in zip(mz, rel, is_matched, partners):
        loss = ""
        if precursor_mz is not None and precursor_mz - m > 0:
            loss = f"<br>Neutral loss {precursor_mz - m:.4f} Da"
        title = f"{label} · {'matched' if matched else 'unmatched'}"
        customdata.append([title, float(m), round(float(r), 1), loss, partner])
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
    """Most intense peaks (>= LABEL_MIN_REL %) whose m/z labels cannot collide."""
    min_sep = LABEL_SEP * x_span
    chosen = []
    for k in np.argsort(-rel, kind="stable"):
        if rel[k] < LABEL_MIN_REL or len(chosen) >= LABELS_PER_HALF:
            break
        if all(abs(mz[k] - mz[c]) >= min_sep for c in chosen):
            chosen.append(int(k))
    return chosen


def _anchor(x, x_range):
    frac = (x - x_range[0]) / (x_range[1] - x_range[0])
    return "right" if frac > 0.8 else ("left" if frac < 0.2 else "center")


def build_mirror_figure(
    peaks_a,
    peaks_b,
    precursor_mz1=None,
    precursor_mz2=None,
    tolerance=0.02,
    labels=("Query", "Reference"),
    height=520,
):
    """Build the interactive mirror-plot figure.

    Parameters
    ----------
    peaks_a, peaks_b : sequences of ``(mz, intensity)`` pairs (query up, reference down).
    precursor_mz1, precursor_mz2 : optional floats; drawn as dotted lines and used for
        the neutral-loss hover text.
    tolerance : Da window used to identify matched peaks (CosineGreedy assignment).
    labels : display names of the query and reference spectra.
    height : figure height in pixels.
    """
    mz1, rel1 = _as_peaks(peaks_a)
    mz2, rel2 = _as_peaks(peaks_b)
    label_a, label_b = labels
    prec1 = _valid_precursor(precursor_mz1)
    prec2 = _valid_precursor(precursor_mz2)

    if mz1.size == 0 and mz2.size == 0:
        return go.Figure()

    matches = cosine_greedy_matches(mz1, rel1, mz2, rel2, tolerance)
    partner_a = {i: j for i, j in matches}
    partner_b = {j: i for i, j in matches}

    # x-range: all peaks and precursors, padded (5 %, at least 5 Da) so edge peaks and
    # their labels are never clipped.
    xs = np.concatenate([mz1, mz2, [p for p in (prec1, prec2) if p is not None]])
    lo, hi = float(xs.min()), float(xs.max())
    pad = max((hi - lo) * 0.05, 5.0)
    x_range = [max(0.0, lo - pad), hi + pad]
    x_span = x_range[1] - x_range[0]

    fig = go.Figure()
    annotations = []
    halves = (
        (mz1, rel1, +1, label_a, prec1, partner_a, mz2, label_b),
        (mz2, rel2, -1, label_b, prec2, partner_b, mz1, label_a),
    )
    for mz, rel, sign, label, prec, partner, other_mz, other_label in halves:
        if mz.size == 0:
            continue
        heights = sign * np.maximum(rel, MIN_STICK)
        is_matched = np.array([k in partner for k in range(mz.size)], dtype=bool)
        partners = [
            f"<br>Matched ↔ {other_label} m/z {other_mz[partner[k]]:.4f} "
            f"(Δ {other_mz[partner[k]] - mz[k]:+.4f} Da)" if k in partner else ""
            for k in range(mz.size)
        ]
        for mask, color, suffix in ((~is_matched, UNMATCHED_COLOR, "unmatched"),
                                    (is_matched, MATCH_COLOR, "matched")):
            if mask.any():
                fig.add_trace(_sticks(mz[mask], heights[mask], color, f"{label} · {suffix}"))
        fig.add_trace(_tip_markers(mz, heights, rel, label, is_matched, prec, partners))

        for k in _label_indices(mz, rel, x_span):
            annotations.append(dict(
                x=float(mz[k]), y=float(heights[k]), xref="x", yref="y",
                yanchor="bottom" if sign > 0 else "top", yshift=3 * sign, showarrow=False,
                text=f"{mz[k]:.4f}", font=dict(family=MONO, size=10, color=theme.INK),
            ))

    # Precursor m/z: dotted guide(s) below the sticks, labelled in the headroom band.
    shapes = []
    same_precursor = prec1 is not None and prec2 is not None and abs(prec1 - prec2) < 1e-6
    for prec, sign, text in ((prec1, +1, "precursor" if same_precursor else f"{label_a} precursor"),
                             (None if same_precursor else prec2, -1, f"{label_b} precursor")):
        if prec is None:
            continue
        y0, y1 = (-100.0, 100.0) if same_precursor else (0.0, sign * 100.0)
        shapes.append(dict(type="line", xref="x", yref="y", x0=prec, x1=prec, y0=y0, y1=y1,
                           layer="below", line=dict(color=theme.MUTED, width=1, dash="dot")))
        annotations.append(dict(
            x=prec, y=sign * Y_TOP, xref="x", yref="y", showarrow=False,
            yanchor="top" if sign > 0 else "bottom", xanchor=_anchor(prec, x_range),
            text=f"{text} {prec:.4f}", font=dict(family=MONO, size=10, color=theme.MUTED),
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
    axis = dict(showline=True, linecolor=theme.INK, linewidth=1.5, ticks="outside",
                tickcolor=theme.INK, ticklen=5, tickfont=dict(color=theme.MUTED, size=11))
    fig.update_layout(
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
            bgcolor=theme.SURFACE,
            bordercolor=theme.INK,
            font=dict(family=FONT, color=theme.INK, size=12),
            align="left",
        ),
        shapes=shapes,
        annotations=annotations,
        xaxis=dict(
            axis, range=x_range, hoverformat=".4f",
            title=dict(text="m/z", font=dict(color=theme.MUTED, size=12)),
            showgrid=False, zeroline=False,
        ),
        yaxis=dict(
            axis, range=[-Y_TOP, Y_TOP], tickvals=ticks, ticktext=[str(abs(t)) for t in ticks],
            showgrid=True, gridcolor=GRID,
            zeroline=True, zerolinecolor=theme.INK, zerolinewidth=1,
        ),
    )
    return fig


def render_mirror(
    peaks_a,
    peaks_b,
    precursor_mz1=None,
    precursor_mz2=None,
    tolerance=0.02,
    height=520,
    labels=("Query", "Reference"),
):
    """Render the mirror plot with Streamlit's interactive Plotly chart."""
    try:
        fig = build_mirror_figure(
            peaks_a,
            peaks_b,
            precursor_mz1=precursor_mz1,
            precursor_mz2=precursor_mz2,
            tolerance=tolerance,
            labels=labels,
            height=height,
        )
    except ValueError as exc:
        st.warning(str(exc))
        return

    if not fig.data:
        st.info("Add at least one valid peak to either spectrum to display the mirror plot.")
        return

    # theme=None: the figure is fully styled from theme.py tokens above; Streamlit's
    # chart theme would otherwise swap the hover-label font and colours.
    st.plotly_chart(fig, width="stretch", height=height, theme=None, config=_PLOT_CONFIG)


if __name__ == "__main__":
    demo_a = [(69.0335, 12), (85.0289, 100), (110.0340, 18), (141.0650, 8)]
    demo_b = [(85.0289, 95), (110.0340, 20), (128.0572, 6), (200.1000, 12)]
    fig = build_mirror_figure(
        demo_a, demo_b, precursor_mz1=250.12, precursor_mz2=250.12, tolerance=0.02
    )
    fig.write_html("/tmp/mirror_demo.html", include_plotlyjs="cdn")
    print("wrote /tmp/mirror_demo.html")
