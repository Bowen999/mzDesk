"""Separation 页的两张图：logD 与净电荷随 pH 的曲线（Plotly），以及标出电离位点的结构图（RDKit SVG）。"""
from typing import Optional

import numpy as np
import plotly.graph_objects as go
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from mzdesk.platform_fit.separation_rules import MINUS, charge_text, logd_text
from mzdesk.shared import colors as theme, depict

LOGD_COLOR = theme.DATA_BLUE
CHARGE_COLOR = theme.DATA_AMBER
ACID = theme.HL_ACID                           # coral tint
BASE = theme.HL_BASE                           # teal tint
BOTH = theme.HL_BOTH                           # violet tint: the same atom is an acid and a base site
LEGEND = "coral: acid site, teal: base site, violet: both"
HEIGHT = 400


def ph_figure(a, ph, state) -> go.Figure:
    """logD(pH) (left axis) and expected net charge (right axis) over pH 0-14, with the chosen pH marked on both."""
    grid = a.ph                                    # python lists below: plotly 6 would otherwise write numpy arrays as base64 "bdata"
    has_logd = bool(np.isfinite(a.logd).any())
    fig = go.Figure()
    axis = dict(showline=True, linecolor=theme.INK, linewidth=1, ticks="outside", tickcolor=theme.INK, ticklen=4,
                tickfont=dict(color=theme.MUTED, size=11), zeroline=False)

    if has_logd:
        fig.add_trace(go.Scatter(x=grid.tolist(), y=a.logd.tolist(), mode="lines", name="logD", line=dict(color=LOGD_COLOR, width=2.5),
                                 hovertemplate="pH %{x:.1f}<br>logD %{y:.2f}<extra></extra>"))
        lo = max(float(np.nanmin(a.logd)), a.logp - 8)
        y_range = [min(lo, -1.0) - 0.3, a.logp + 0.6]
    fig.add_trace(go.Scatter(x=grid.tolist(), y=a.mean_charge.tolist(), mode="lines", name="Net charge", yaxis="y2" if has_logd else "y",
                             line=dict(color=CHARGE_COLOR, width=2.5, dash="dash"),
                             hovertemplate="pH %{x:.1f}<br>net charge %{y:.2f}<extra></extra>"))
    # the chosen pH
    marks = [(ph, state["mean_charge"], "y2" if has_logd else "y", CHARGE_COLOR, f"net charge {state['mean_charge']:+.2f}".replace("-", MINUS))]
    if has_logd and state["logd"] is not None:
        marks.append((ph, state["logd"], "y", LOGD_COLOR, f"logD {logd_text(state['logd'])}"))
    for x, y, ax, color, text in marks:
        fig.add_trace(go.Scatter(x=[x], y=[y], mode="markers", yaxis=ax, showlegend=False, marker=dict(size=10, color=color, line=dict(color="white", width=1.5)),
                                 hovertemplate=f"pH {ph:g}<br>{text}<extra></extra>"))
    q_lo, q_hi = min(float(a.charges.min()), -1.0) - 0.4, max(float(a.charges.max()), 1.0) + 0.4
    ticks = list(range(int(np.ceil(q_lo)), int(np.floor(q_hi)) + 1))
    ticks = ticks if len(ticks) <= 9 else ticks[::2]
    charge_axis = dict(tickvals=ticks, ticktext=[charge_text(t) for t in ticks])

    layout = dict(
        template="none", height=HEIGHT, margin=dict(l=62, r=66 if has_logd else 22, t=34, b=50), showlegend=False,
        paper_bgcolor=theme.SURFACE, plot_bgcolor=theme.SURFACE, font=dict(family=theme.PLOT_FONT, color=theme.INK, size=12),
        hovermode="closest", dragmode="zoom",
        hoverlabel=dict(bgcolor=theme.PAPER, bordercolor=theme.LINE, font=dict(family=theme.PLOT_FONT, color=theme.INK, size=12), align="left"),
        xaxis=dict(axis, title=dict(text="pH", font=dict(color=theme.MUTED, size=12)), range=[0, 14], dtick=2, showgrid=False),
        shapes=[dict(type="line", xref="x", yref="paper", x0=ph, x1=ph, y0=0, y1=1, layer="below", line=dict(color=theme.MUTED, width=1, dash="dot"))],
        annotations=[dict(x=ph, xref="x", y=1, yref="paper", yanchor="bottom", yshift=4, showarrow=False, text=f"pH {ph:g}",
                          xanchor="right" if ph > 11 else ("left" if ph < 3 else "center"), font=dict(family=theme.PLOT_MONO, size=10, color=theme.MUTED))])
    if has_logd:
        layout["yaxis"] = dict(axis, title=dict(text="logD", font=dict(color=LOGD_COLOR, size=12)), range=y_range, showgrid=True, gridcolor=theme.SURFACE_2)
        layout["yaxis2"] = dict(axis, title=dict(text="Net charge", font=dict(color=CHARGE_COLOR, size=12)), range=[q_lo, q_hi], overlaying="y", side="right",
                                showgrid=False, **charge_axis)
    else:
        layout["yaxis"] = dict(axis, title=dict(text="Net charge", font=dict(color=CHARGE_COLOR, size=12)), range=[q_lo, q_hi], showgrid=True, gridcolor=theme.SURFACE_2,
                               **charge_axis)
        layout["annotations"].append(dict(x=0.5, xref="paper", y=0.97, yref="paper", showarrow=False, yanchor="top",
                                          text="logD is not defined: permanent charge, no neutral form", font=dict(size=11, color=theme.MUTED)))
    fig.update_layout(**layout)
    return fig


def sites_svg(a, size=(560, 420)) -> Optional[str]:
    """The neutral molecule with every ionizable site highlighted (acid / base / both) and labelled with atom number and pKa."""
    mol = Chem.Mol(a.mol)
    by_atom = {}
    for s in a.sites:
        by_atom.setdefault(s.atom, []).append(s)
    if not by_atom:
        return None
    colors = {}
    for idx, group in by_atom.items():
        kinds = {s.kind for s in group}
        colors[idx] = BOTH if len(kinds) == 2 else (ACID if kinds == {"acid"} else BASE)
        mol.GetAtomWithIdx(idx).SetProp("atomNote", f"{group[0].label} " + "/".join(f"{s.pka:.1f}" for s in group))
    drawer = depict.svg_drawer(*size)
    opts = drawer.drawOptions()
    opts.annotationFontScale = 0.75
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol, highlightAtoms=list(colors), highlightAtomColors=colors)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()
