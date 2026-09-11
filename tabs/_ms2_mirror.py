"""Interactive mirror MS² plot for the MS2 Similarity tab.

The MS² Similarity tab imports this module as ``_ms2_mirror`` and calls
``render_mirror(...)``. ``build_mirror_figure`` returns a Plotly figure;
``render_mirror`` draws it with Streamlit's native interactive chart.

Colour logic:
  * Query Spec (up)  — gray        ``#B8BEC7``
  * Ref Spec (down)  — light blue  ``#AFC7EA``
  * Matched peaks    — orange      ``#E28A3C``
"""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
import streamlit as st


INK       = "#0A1633"
MUTED     = "#5C6B8A"
GRID      = "#E6EAF2"
QUERY_UNMATCHED = "#B8BEC7"   # gray
REF_UNMATCHED   = "#AFC7EA"   # light blue
MATCH           = "#E28A3C"   # orange
FONT      = "Inter, 'Helvetica Neue', Arial, sans-serif"


def _as_peaks(peaks):
    """Coerce ``[(mz, intensity), ...]`` into two cleaned numpy arrays."""
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


def _greedy_matches(mz1, int1, mz2, int2, tolerance):
    """Return a 1:1 greedy list of ``(i, j)`` matched peak-index pairs."""
    if tolerance is None or mz1.size == 0 or mz2.size == 0:
        return []

    used2 = set()
    matches = []
    for i in np.argsort(-int1):  # strongest query peaks first
        lo = np.searchsorted(mz2, mz1[i] - tolerance)
        hi = np.searchsorted(mz2, mz1[i] + tolerance, side="right")
        best = None
        for j in range(lo, hi):
            if j in used2:
                continue
            if best is None or int2[j] > int2[best]:
                best = j
        if best is not None:
            used2.add(int(best))
            matches.append((int(i), int(best)))
    return matches


def match_summary(peaks_a, peaks_b, tolerance=0.02):
    """Return matched-peak counts and percentages for query and reference spectra."""
    mz1, int1 = _as_peaks(peaks_a)
    mz2, int2 = _as_peaks(peaks_b)
    matches = _greedy_matches(mz1, int1, mz2, int2, tolerance)

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


def _split_by_mask(mz, intensity, matched_idx):
    """Split a spectrum into (unmatched, matched) arrays."""
    mask = np.zeros(len(mz), dtype=bool)
    mask[list(matched_idx)] = True
    return mz[~mask], intensity[~mask], mz[mask], intensity[mask]


def _bar_trace(mz, intensity, sign, color, name, precursor_mz, width, showlegend):
    """Mass-spectrum bars with neutral-loss hover information."""
    customdata = []
    for m, i in zip(mz, intensity):
        if precursor_mz is not None:
            loss = float(precursor_mz) - float(m)
            nl = f"{loss:.4f}" if loss > 0 else "—"
        else:
            nl = "—"
        customdata.append([round(float(i), 1), nl])

    hovertemplate = (
        "<b>%{fullData.name}</b><br>"
        "m/z %{x:.4f}<br>"
        "Rel. intensity %{customdata[0]:.1f}%<br>"
        "Neutral loss %{customdata[1]}"
        "<extra></extra>"
    )
    return go.Bar(
        x=mz,
        y=sign * intensity,
        marker_color=color,
        width=width,
        customdata=customdata,
        hovertemplate=hovertemplate,
        name=name,
        showlegend=showlegend,
    )


def _matched_connector(mz1, int1, mz2, int2, matches):
    """Dotted purple lines connecting matched query/reference peaks."""
    if not matches:
        return None
    xs, ys = [], []
    for i, j in matches:
        xs.extend([mz1[i], mz2[j], None])
        ys.extend([float(int1[i]), -float(int2[j]), None])
    return go.Scatter(
        x=xs,
        y=ys,
        mode="lines",
        line=dict(color=MATCH, width=1, dash="dot"),
        hoverinfo="skip",
        name="Matched",
        showlegend=True,
    )


def build_mirror_figure(
    peaks_a,
    peaks_b,
    precursor_mz1=None,
    precursor_mz2=None,
    tolerance=0.02,
    labels=("Query Spec", "Ref Spec"),
    height=520,
):
    """Build the interactive mirror-plot figure.

    Parameters
    ----------
    peaks_a, peaks_b : sequences of ``(mz, intensity)`` pairs.
    precursor_mz1, precursor_mz2 : optional floats for neutral-loss hover.
    tolerance : Da window used to identify matched peaks.
    labels : legend names for the query and reference spectra.
    height : figure height in pixels.
    """
    mz1, int1 = _as_peaks(peaks_a)
    mz2, int2 = _as_peaks(peaks_b)
    label_a, label_b = labels

    if mz1.size == 0 and mz2.size == 0:
        return go.Figure()

    matches = _greedy_matches(mz1, int1, mz2, int2, tolerance)
    matched_a = {i for i, _ in matches}
    matched_b = {j for _, j in matches}

    mz1_u, int1_u, mz1_m, int1_m = _split_by_mask(mz1, int1, matched_a)
    mz2_u, int2_u, mz2_m, int2_m = _split_by_mask(mz2, int2, matched_b)

    mz_min = min(mz1.min(), mz2.min()) if mz1.size and mz2.size else (
        mz1.min() if mz1.size else mz2.min()
    )
    mz_max = max(mz1.max(), mz2.max()) if mz1.size and mz2.size else (
        mz1.max() if mz1.size else mz2.max()
    )
    bar_width = max((mz_max - mz_min) * 0.002, 0.008)

    fig = go.Figure()

    # Query Spec: unmatched blue, matched purple.
    fig.add_trace(_bar_trace(
        mz1_u, int1_u, +1, QUERY_UNMATCHED, label_a, precursor_mz1, bar_width, True,
    ))
    fig.add_trace(_bar_trace(
        mz1_m, int1_m, +1, MATCH, "Matched (Query)", precursor_mz1, bar_width, False,
    ))

    # Ref Spec: unmatched red, matched purple.
    fig.add_trace(_bar_trace(
        mz2_u, int2_u, -1, REF_UNMATCHED, label_b, precursor_mz2, bar_width, True,
    ))
    fig.add_trace(_bar_trace(
        mz2_m, int2_m, -1, MATCH, "Matched (Ref)", precursor_mz2, bar_width, False,
    ))

    connector = _matched_connector(mz1, int1, mz2, int2, matches)
    if connector is not None:
        fig.add_trace(connector)

    fig.update_layout(
        title=dict(text="Mirror MS² plot", font=dict(color=INK, size=18)),
        font=dict(family=FONT, color=INK),
        plot_bgcolor="#FFFFFF",
        paper_bgcolor="#FFFFFF",
        height=height,
        margin=dict(l=52, r=24, t=58, b=52),
        bargap=0,
        hovermode="closest",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1,
        ),
        xaxis=dict(
            title=dict(text="m/z", font=dict(color=MUTED, size=13)),
            showline=True,
            linecolor=INK,
            linewidth=1,
            gridcolor=GRID,
            zeroline=False,
            tickfont=dict(color=MUTED),
        ),
        yaxis=dict(
            title=dict(text="Relative intensity (%)", font=dict(color=MUTED, size=13)),
            tickvals=[-100, -75, -50, -25, 0, 25, 50, 75, 100],
            ticktext=["100", "75", "50", "25", "0", "25", "50", "75", "100"],
            tickfont=dict(color=MUTED),
            showline=True,
            linecolor=INK,
            linewidth=1,
            gridcolor=GRID,
            zeroline=True,
            zerolinecolor=INK,
            zerolinewidth=1,
            range=[-108, 108],
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
    labels=("Query Spec", "Ref Spec"),
):
    """Render the mirror plot with Streamlit's interactive chart."""
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

    st.plotly_chart(fig, width="stretch", height=height)


if __name__ == "__main__":
    demo_a = [(69.0335, 12), (85.0289, 100), (110.0340, 18), (141.0650, 8)]
    demo_b = [(85.0289, 95), (110.0340, 20), (128.0572, 6), (200.1000, 12)]
    fig = build_mirror_figure(
        demo_a, demo_b, precursor_mz1=250.12, precursor_mz2=250.12, tolerance=0.02
    )
    fig.write_html("/tmp/mirror_demo.html", include_plotlyjs="cdn")
    print("wrote /tmp/mirror_demo.html")
