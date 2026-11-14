"""MS2 谱图基础可视化（x = m/z，y = 强度）。清洗页、镜像图、MetFrag 图都从这里取样式、棒图、悬停点和前体离子标记。"""
import math

import numpy as np
import plotly.graph_objects as go

from mzdesk.shared import colors as theme

DARK = theme.PLOT_DARK        # kept / explained / matched
LIGHT = theme.PLOT_LIGHT      # removed / unexplained / unmatched
FONT, MONO = theme.PLOT_FONT, theme.PLOT_MONO
HEIGHT = 460
STICK_PX = 2                  # stick width; dense spectra get thinner sticks (stick_px)
HEADROOM = 1.15               # y axis top = tallest peak x this, room for the m/z labels
PRECURSOR_MATCH = 0.05        # default window (Da) for indicating a peak near the precursor m/z


def axis_style():
    """x axis: a quiet baseline with short outside ticks, muted tick font (the sticks stand on the baseline)."""
    return dict(showline=True, linecolor=theme.AXIS, linewidth=1, ticks="outside", tickcolor=theme.AXIS,
                ticklen=4, tickfont=dict(color=theme.MUTED, size=11), zeroline=False)


def y_axis_style():
    """y axis: no line, only the light gridlines; the invisible ticks keep the tick labels off the plot."""
    return dict(showline=False, ticks="outside", tickcolor="rgba(0,0,0,0)", ticklen=6,
                tickfont=dict(color=theme.MUTED, size=11), zeroline=False, showgrid=True, gridcolor=theme.SURFACE_2)


def axis_title(text):
    return dict(text=text, font=dict(color=theme.MUTED, size=12))


def mz_range(mz, pad_min=5.0):
    """x range around the m/z values: 5 % padding (at least ``pad_min``), never below 0."""
    lo, hi = float(min(mz)), float(max(mz))
    pad = max((hi - lo) * 0.05, pad_min)
    return [max(0.0, lo - pad), hi + pad]


def percent_ticks(top=100):
    """y ticks every 20 % up to 100: the headroom above the base peak is for labels, not for values."""
    return dict(tickvals=list(range(0, int(min(top, 100)) + 1, 20)))


def base_layout(height=HEIGHT, *, x_range, y_range, y_title, yaxis=None, margin=None, shapes=(), annotations=()):
    """Layout dict for ``fig.update_layout``: transparent, light y grid, no legend, hover by x, drag to zoom."""
    y = dict(y_axis_style(), title=axis_title(y_title) if y_title else None, range=y_range)
    y.update(yaxis or {})
    return dict(
        template="none", height=height, margin=margin or dict(l=62, r=22, t=30, b=50),
        paper_bgcolor=theme.SURFACE, plot_bgcolor=theme.SURFACE, font=dict(family=FONT, color=theme.INK, size=12),
        showlegend=False, hovermode="x", hoverdistance=16, dragmode="zoom",
        hoverlabel=dict(bgcolor=theme.PAPER, bordercolor=theme.LINE, font=dict(family=FONT, color=theme.INK, size=12), align="left"),
        shapes=list(shapes), annotations=list(annotations),
        xaxis=dict(axis_style(), title=axis_title("m/z"), range=x_range, showgrid=False, hoverformat=".4f"), yaxis=y)


def stick_px(mz, x_range):
    """Stick width for this many peaks across ``x_range``: 2 px, thinner when sticks would run together (unit-mass EI)."""
    per_100 = len(mz) / max(x_range[1] - x_range[0], 1e-9) * 100
    return STICK_PX if per_100 <= 20 else 1.5 if per_100 <= 45 else 1.0


def stick_trace(mz, y, color, name, y0=0.0, width=STICK_PX):
    """All sticks as ONE lines trace (None breaks): constant pixel width, never merges when zoomed."""
    xs, ys = [], []
    for m, h in zip(mz, y):
        xs += [float(m), float(m), None]
        ys += [y0, float(h), None]
    return go.Scatter(x=xs, y=ys, mode="lines", line=dict(color=color, width=width), hoverinfo="skip",
                      name=name, showlegend=False)


def hover_markers(mz, y, *, name="peaks", customdata=None, hovertemplate=None):
    """Invisible marker per peak at its tip. Trace 'peaks' + customdata = peak index drives the hover card
    and auto-labels in web/src/components/Figure.tsx; give ``hovertemplate`` for a native Plotly hover instead."""
    hover = dict(hovertemplate=hovertemplate) if hovertemplate else dict(hoverinfo="none")
    return go.Scatter(x=mz, y=y, mode="markers", name=name, showlegend=False, marker=dict(size=10, color="rgba(0,0,0,0)"),
                      customdata=list(range(len(mz))) if customdata is None else customdata, **hover)


def precursor_line(x, y0, y1):
    """Dotted vertical line from y0 to y1 (Figure.tsx finds the precursor by this dash and keeps labels off it)."""
    return dict(type="line", xref="x", yref="y", x0=x, x1=x, y0=y0, y1=y1, layer="below",
                line=dict(color=theme.MUTED, width=1, dash="dot"))


def precursor_label(x, text, x_range, **pos):
    """Small mono label; default position just above the plot, anchored away from the edges."""
    frac = (x - x_range[0]) / (x_range[1] - x_range[0])
    ann = dict(x=x, xref="x", y=1, yref="paper", yanchor="bottom", yshift=4, showarrow=False, text=text,
               xanchor="right" if frac > 0.8 else ("left" if frac < 0.2 else "center"),
               font=dict(family=MONO, size=10, color=theme.MUTED))
    return dict(ann, **pos)


def fmt_intensity(y):
    return f"{y:,.0f}" if abs(y) >= 1e4 else f"{y:.4g}"


def fmt_mz(value):
    """Keep distinct float values distinct, with at least four decimal places when practical."""
    text = repr(float(value))
    if "e" in text:
        return text  # Avoid an unwieldy run of zeros for very small or large values.
    whole, fraction = text.split(".")
    return whole + "." + fraction.ljust(4, "0")


def fmt_rel(value):
    """A small positive percentage must not be rounded to zero."""
    return f"{value:.3g}" if 0 < abs(value) < 1 else f"{value:.1f}"


def spectrum_figure(peaks, *, precursor_mz=None, precursor_tol=PRECURSOR_MATCH, y_title="Intensity", dim=None, height=HEIGHT):
    """Stick spectrum of ``peaks`` [(mz, y)] -> ``(fig, cards)``; pass cards to ``common.plotly(fig, peaks=cards)``.

    ``dim`` marks removed peaks (drawn light, card kind "Removed"): a list of indices, a bool mask,
    or {index: reason} (kind "Removed: reason"). The y axis starts at 0 with headroom for labels
    and uses SI ticks for large raw intensities. A dotted line + label marks ``precursor_mz``;
    ``precursor_tol`` controls the window used to indicate a nearby peak, not peak processing.
    Zero-height peaks keep their positions and cards; their relative intensity is undefined when
    every peak is zero. Invalid peak values are rejected rather than silently changing indices.
    """
    peaks = [(float(m), float(y)) for m, y in peaks]
    if not peaks:
        return go.Figure(), []
    if any(not math.isfinite(m) or not math.isfinite(y) or m <= 0 or y < 0 for m, y in peaks):
        raise ValueError("peaks must have positive finite m/z and non-negative finite intensity")
    n = len(peaks)
    if isinstance(dim, dict):
        off = {int(i): str(r) for i, r in dim.items()}
    elif dim is None:
        off = {}
    else:
        arr = np.asarray(dim)
        off = {int(i): "" for i in (np.flatnonzero(arr) if arr.dtype == bool else arr)}
    mz, ys = [p[0] for p in peaks], [p[1] for p in peaks]
    ymax = max(ys)
    top = ymax * HEADROOM if ymax > 0 else 1.0
    xr = mz_range(mz + ([precursor_mz] if precursor_mz else []))

    fig = go.Figure()
    width = stick_px(mz, xr)
    for name, color, idx in (("dimmed", LIGHT, sorted(off)), ("kept", DARK, [i for i in range(n) if i not in off])):
        if idx:
            fig.add_trace(stick_trace([mz[i] for i in idx], [ys[i] for i in idx], color, name, width=width))
    fig.add_trace(hover_markers(mz, ys))

    shapes, annotations, prec = [], [], None
    if precursor_mz:
        near = [i for i in range(n) if abs(mz[i] - precursor_mz) <= precursor_tol]
        prec = min(near, key=lambda i: abs(mz[i] - precursor_mz)) if near else None
        y0 = 0.0 if prec is None else ys[prec] + 0.08 * top      # start above the precursor peak's label
        if top - y0 >= 0.1 * top:                                 # no stub between a tall precursor peak and the top
            shapes.append(precursor_line(precursor_mz, y0, top))
        annotations.append(precursor_label(precursor_mz, f"precursor {precursor_mz:.4f}", xr))
    yaxis = dict(tickformat="~s") if ymax >= 1e4 else percent_ticks() if abs(ymax - 100) < 0.5 else None
    fig.update_layout(**base_layout(
        height, x_range=xr, y_range=[0, top], y_title=y_title, shapes=shapes, annotations=annotations, yaxis=yaxis))

    cards = []
    for i, (m, y) in enumerate(peaks):
        kind = (("Removed: " + off[i] if off[i] else "Removed") if i in off else
                ("Precursor ion" if m == precursor_mz else "Near precursor m/z") if i == prec else "Peak")
        rows = [["Intensity", fmt_intensity(y)],
                ["Rel. intensity", f"{fmt_rel(y / ymax * 100)} %" if ymax > 0 else "—"]]
        # Figure.tsx displays four decimals in its heading. Keep the stored m/z and expose
        # extra digits in the card without changing the shared frontend's contract.
        if float(f"{m:.4f}") != m:
            rows.append(["m/z (more precision)", fmt_mz(m)])
        cards.append({"mz": m, "rel": y, "on": i not in off, "kind": kind,
                      "rows": rows})
    return fig, cards
