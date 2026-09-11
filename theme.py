"""theme.py — unified design system for the Metabolomics Workbench app.

Design language: NEO-BRUTALISM x MINIMALISM, cool-white + indigo ("high-end blue").
  - thick navy-ink borders, hard offset shadows (no blur), near-square corners
  - uppercase mono kickers / section labels, restrained indigo-blue palette
  - generous whitespace, hairline dividers, few decorative elements

PUBLIC API (import as `import theme`):

  Colours (hex strings):
    theme.INK, BG, SURFACE, SURFACE_2, MUTED, LINE,
    theme.BLUE, BLUE_BRIGHT, BLUE_SOFT, BLUE_GHOST
    theme.CATEGORY_ACCENT  -> {"Lipidome","Metabolome","Both","Neither","Invalid"}
    theme.PLOT_DARK, PLOT_LIGHT     (MS² bar colours)
    theme.FRAG_HL                   (RDKit substructure-highlight RGB 0..1 tuples)

  Render helpers (call st.* directly — use inside a tab / page):
    theme.inject_css()                         # once, right after set_page_config
    theme.hero(title, subtitle=None, kicker=None)
    theme.tab_header(title, subtitle=None, kicker=None)
    theme.section(label)                       # small uppercase section label
    theme.metric_cards(items)                  # items: [(label, value), ...]
    theme.stat_grid(rows)                      # rows:  [(label, value), ...]
    theme.verdict(category, subtitle, accent)  # big result banner
    theme.chips(items)                         # items: [(label, css_color), ...]
    theme.similarity_score_cards(groups)       # premium score cards for grouped metrics
    theme.copy_fields(fields, heading=None)    # identity rows with copy buttons
    theme.under_construction(title, note=None, eta=None)
    theme.divider()
    theme.note(text)

  HTML-string builders (return str; caller wraps in st.markdown(..., unsafe_allow_html=True)):
    theme.card_html, section_label_html, stat_grid_html, metric_cards_html,
    theme.chips_html, similarity_score_cards_html, verdict_html
"""
from __future__ import annotations

import base64
from html import escape
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

import streamlit as st
import streamlit.components.v1 as components

# --------------------------------------------------------------------------
# Colour tokens
# --------------------------------------------------------------------------
INK        = "#0A1633"   # near-black navy — borders & primary text
BG         = "#F3F5FA"   # cool near-white background
SURFACE    = "#FFFFFF"   # card surface
SURFACE_2  = "#E9EDF7"   # subtle panel / input fill
MUTED      = "#5C6B8A"   # secondary text (cool grey-blue)
LINE       = "#C7D0E4"   # hairline divider

BLUE       = "#1B3A8F"   # primary indigo
BLUE_BRIGHT= "#3B5BDB"   # hover / bright accent
BLUE_SOFT  = "#A9BEEA"   # soft fill / unannotated
BLUE_GHOST = "#EEF2FC"   # very light blue tint

CATEGORY_ACCENT = {
    "Lipidome":   "#1B3A8F",  # deep indigo
    "Metabolome": "#0E7C9B",  # teal-cyan
    "Both":       "#5566CC",  # periwinkle
    "Neither":    "#8A94A6",  # cool grey
    "Invalid":    "#C0392B",  # functional error red
}

# MS² spectrum bar colours
PLOT_DARK  = BLUE          # annotated peaks
PLOT_LIGHT = BLUE_SOFT     # unannotated peaks

# RDKit substructure highlight colours (RGB 0..1), harmonised to the cool palette
FRAG_HL = {
    "acyl":       (0.106, 0.227, 0.561),  # indigo   >=C12 acyl chain
    "isoprenoid": (0.055, 0.486, 0.608),  # teal     sterol / isoprenoid
    "head":       (0.333, 0.400, 0.800),  # periwinkle polar head group
    "glycosyl":   (0.902, 0.624, 0.208),  # amber    glycosyl ring (single warm anchor)
}

# --------------------------------------------------------------------------
# Global CSS
# --------------------------------------------------------------------------
GLOBAL_CSS = f"""
<style>
:root{{
  --ink:{INK}; --bg:{BG}; --surface:{SURFACE}; --surface-2:{SURFACE_2};
  --muted:{MUTED}; --line:{LINE};
  --blue:{BLUE}; --blue-bright:{BLUE_BRIGHT}; --blue-soft:{BLUE_SOFT}; --blue-ghost:{BLUE_GHOST};
}}

/* ---- base ---- */
html, body, [data-testid="stAppViewContainer"], [data-testid="stApp"]{{ background:var(--bg); }}
[data-testid="stHeader"]{{ background:transparent; }}
#MainMenu, footer, [data-testid="stToolbar"], [data-testid="stDecoration"]{{ visibility:hidden; }}
.block-container{{ padding-top:1.4rem; padding-bottom:3rem; max-width:1240px; }}
html, body, [class*="css"]{{
  font-family:"Inter","Helvetica Neue","Segoe UI",Arial,sans-serif;
  color:var(--ink);
}}
.mono{{ font-family:ui-monospace,"SF Mono",Menlo,Consolas,monospace; }}
a{{ color:var(--blue); }}

/* ---- hero / headers ---- */
.nb-kicker{{
  font-family:ui-monospace,"SF Mono",Menlo,monospace; font-size:.72rem;
  text-transform:uppercase; letter-spacing:.20em; color:var(--blue);
  font-weight:700; margin:0 0 .5rem 0;
}}
.nb-hero{{
  text-transform:uppercase; letter-spacing:.005em; font-weight:900;
  font-size:clamp(2.1rem,4vw,3rem); line-height:1.0; margin:0; color:var(--ink);
}}
.nb-hero .accent{{ color:var(--blue); }}
.nb-sub{{ font-size:.95rem; color:var(--muted); margin:.55rem 0 0 0; max-width:70ch; }}
.nb-rule{{ height:3px; background:var(--ink); border:0; margin:1.1rem 0 1.4rem 0; }}
.nb-hairline{{ height:1px; background:var(--line); border:0; margin:1.3rem 0; }}

.nb-tabhead{{ margin:.2rem 0 1.1rem 0; }}
.nb-tabhead .t{{ font-size:1.55rem; font-weight:900; text-transform:uppercase;
  letter-spacing:.01em; line-height:1.05; margin:0; }}
.nb-tabhead .s{{ font-size:.9rem; color:var(--muted); margin:.35rem 0 0 0; }}

/* ---- section label ---- */
.nb-sec{{
  display:flex; align-items:center; gap:.55rem;
  font-size:.72rem; text-transform:uppercase; letter-spacing:.14em;
  font-weight:800; color:var(--muted); margin:.2rem 0 .7rem 0;
}}
.nb-sec::before{{ content:""; width:11px; height:11px; background:var(--blue);
  border:2px solid var(--ink); display:inline-block; }}

/* ---- cards ---- */
.nb-card{{
  background:var(--surface); border:3px solid var(--ink);
  box-shadow:5px 5px 0 var(--ink); border-radius:2px;
  padding:16px 18px; margin-bottom:16px;
}}
.nb-card.ghost{{ background:var(--blue-ghost); }}

/* ---- verdict banner ---- */
.nb-verdict{{
  border:4px solid var(--ink); box-shadow:8px 8px 0 var(--ink);
  border-radius:2px; padding:20px 24px; margin-bottom:18px; color:#fff;
}}
.nb-verdict .cat{{ font-size:2.5rem; font-weight:900; text-transform:uppercase;
  line-height:1; letter-spacing:.01em; }}
.nb-verdict .sub{{ font-size:1rem; margin-top:8px; opacity:.92; font-weight:500; }}

/* ---- stat grid ---- */
.nb-statgrid{{ display:grid; grid-template-columns:repeat(2,1fr); gap:12px; }}
.nb-statgrid.wide{{ grid-template-columns:repeat(auto-fit,minmax(160px,1fr)); }}
.nb-stat{{ background:var(--surface); border:3px solid var(--ink);
  box-shadow:4px 4px 0 var(--ink); border-radius:2px; padding:10px 13px; }}
.nb-stat .lab{{ font-size:.66rem; text-transform:uppercase; letter-spacing:.08em;
  font-weight:800; color:var(--muted); }}
.nb-stat .val{{ font-size:1.4rem; font-weight:800; color:var(--ink);
  font-family:ui-monospace,Menlo,monospace; margin-top:3px; word-break:break-word; }}

/* ---- metric cards (KPI row) ---- */
.nb-metricrow{{ display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr));
  gap:12px; margin-bottom:6px; }}
.nb-metric{{ background:var(--ink); color:#fff; border:3px solid var(--ink);
  box-shadow:5px 5px 0 var(--blue); border-radius:2px; padding:13px 15px; }}
.nb-metric .lab{{ font-size:.64rem; text-transform:uppercase; letter-spacing:.1em;
  font-weight:800; color:var(--blue-soft); }}
.nb-metric .val{{ font-size:1.55rem; font-weight:900;
  font-family:ui-monospace,Menlo,monospace; margin-top:4px; word-break:break-word; }}

/* ---- chips ---- */
.nb-chip{{ display:inline-flex; align-items:center; gap:7px; border:2px solid var(--ink);
  border-radius:2px; padding:4px 9px; margin:4px 7px 4px 0; font-size:.71rem;
  font-weight:800; text-transform:uppercase; letter-spacing:.03em; background:var(--surface); }}
.nb-chip .sw{{ width:12px; height:12px; border:1.5px solid var(--ink); display:inline-block; }}

/* ---- identity panel with copy buttons ---- */
.nb-identity{{
  background:var(--surface); border:3px solid var(--ink);
  box-shadow:5px 5px 0 var(--ink); border-radius:2px;
  padding:14px 16px 10px; margin-bottom:16px;
}}
.nb-identity .nb-sec{{ margin:.1rem 0 .55rem 0; }}
.nb-identity .row{{
  display:flex; align-items:flex-start; gap:12px;
  border-top:1px solid var(--line); padding:8px 0;
}}
.nb-identity .row.first{{ border-top:0; padding-top:2px; }}
.nb-identity .lab{{
  flex:0 0 150px; font-size:.64rem; text-transform:uppercase;
  letter-spacing:.08em; font-weight:800; color:var(--muted); padding-top:3px;
}}
.nb-identity .val{{
  flex:1; font-family:ui-monospace,"SF Mono",Menlo,Consolas,monospace;
  font-size:.76rem; line-height:1.55; color:var(--ink);
  word-break:break-all; min-width:0;
}}
.nb-identity .copy-btn{{
  flex:none; width:20px; height:20px; margin-top:-1px; padding:0;
  display:inline-flex; align-items:center; justify-content:center;
  background:var(--blue-ghost); border:1.5px solid var(--ink); border-radius:2px;
  box-shadow:2px 2px 0 var(--ink); cursor:pointer; color:var(--blue);
}}
.nb-identity .copy-btn:hover{{ background:var(--blue-soft); }}
.nb-identity .copy-btn:active{{ transform:translate(2px,2px); box-shadow:none; }}
.nb-identity .copy-btn svg{{ width:11px; height:11px; display:block; }}
.nb-identity .copy-btn .icon-ok{{ display:none; }}
.nb-identity .copy-btn.copied{{ color:#1C7A3D; }}
.nb-identity .copy-btn.copied .icon-copy{{ display:none; }}
.nb-identity .copy-btn.copied .icon-ok{{ display:block; }}
.nb-identity .copy-spacer{{ flex:none; width:20px; }}

/* ---- under construction ---- */
.nb-uc{{ border:3px dashed var(--ink); border-radius:2px; background:var(--blue-ghost);
  padding:44px 32px; text-align:center; margin-top:8px; }}
.nb-uc .tag{{ font-family:ui-monospace,Menlo,monospace; font-size:.72rem; font-weight:700;
  text-transform:uppercase; letter-spacing:.2em; color:var(--blue);
  border:2px solid var(--ink); background:var(--surface); display:inline-block;
  padding:4px 12px; box-shadow:3px 3px 0 var(--ink); }}
.nb-uc .t{{ font-size:1.9rem; font-weight:900; text-transform:uppercase; margin:18px 0 6px 0; }}
.nb-uc .n{{ font-size:.92rem; color:var(--muted); max-width:52ch; margin:0 auto; }}
.nb-uc .eta{{ margin-top:16px; font-size:.75rem; font-weight:800; text-transform:uppercase;
  letter-spacing:.08em; color:var(--ink); }}

.nb-note{{ font-size:.82rem; color:var(--muted); margin:.4rem 0; }}

/* ---- Streamlit widgets ---- */
.stButton>button, [data-testid="stBaseButton-secondary"]{{
  background:var(--blue); color:#fff; border:3px solid var(--ink);
  box-shadow:4px 4px 0 var(--ink); border-radius:2px; font-weight:800;
  text-transform:uppercase; letter-spacing:.04em; padding:.5rem 1.3rem; transition:none;
}}
.stButton>button:hover, [data-testid="stBaseButton-secondary"]:hover{{
  background:var(--blue-bright); color:#fff; border-color:var(--ink); }}
.stButton>button:active{{ transform:translate(4px,4px); box-shadow:none; }}
[data-testid="stBaseButton-primary"]{{
  background:var(--ink); color:#fff; border:3px solid var(--ink);
  box-shadow:4px 4px 0 var(--blue); border-radius:2px; font-weight:800;
  text-transform:uppercase; letter-spacing:.04em; }}
[data-testid="stBaseButton-primary"]:hover{{ background:#111f42; color:#fff; }}
[data-testid="stBaseButton-primary"]:active{{ transform:translate(4px,4px); box-shadow:none; }}

[data-testid="stTextInput"] input, [data-testid="stNumberInput"] input,
[data-testid="stTextArea"] textarea{{
  border:3px solid var(--ink)!important; border-radius:2px!important;
  box-shadow:3px 3px 0 var(--ink); background:var(--surface)!important;
  font-family:ui-monospace,Menlo,monospace; color:var(--ink)!important; }}
[data-testid="stTextInput"] input:focus, [data-testid="stNumberInput"] input:focus,
[data-testid="stTextArea"] textarea:focus{{ border-color:var(--blue)!important; }}
div[data-baseweb="select"]>div{{
  border:3px solid var(--ink)!important; border-radius:2px!important;
  box-shadow:3px 3px 0 var(--ink); background:var(--surface)!important; }}
label, [data-testid="stWidgetLabel"] p{{ font-weight:700!important; color:var(--ink)!important;
  font-size:.82rem!important; }}

/* number-input steppers */
[data-testid="stNumberInputStepUp"], [data-testid="stNumberInputStepDown"]{{
  border-left:2px solid var(--ink)!important; }}

/* slider */
[data-testid="stSlider"] [data-baseweb="slider"] div[role="slider"]{{
  border:3px solid var(--ink); background:var(--blue); box-shadow:2px 2px 0 var(--ink); }}

/* sidebar */
section[data-testid="stSidebar"]{{ background:var(--surface-2); border-right:3px solid var(--ink); }}
section[data-testid="stSidebar"] .block-container{{ padding-top:2rem; }}

/* images */
[data-testid="stImage"] img{{ border:3px solid var(--ink); box-shadow:5px 5px 0 var(--ink);
  border-radius:2px; background:#fff; }}

/* tabs */
[data-baseweb="tab-list"]{{ gap:8px; border-bottom:3px solid var(--ink);
  padding-bottom:0; margin-bottom:1.2rem; }}
[data-baseweb="tab-list"] button[data-baseweb="tab"]{{
  background:transparent; border:3px solid transparent; border-radius:2px 2px 0 0;
  padding:9px 16px; margin-bottom:-3px; font-weight:800; text-transform:uppercase;
  letter-spacing:.04em; font-size:.82rem; color:var(--muted); }}
[data-baseweb="tab-list"] button[data-baseweb="tab"]:hover{{ color:var(--ink); background:var(--blue-ghost); }}
[data-baseweb="tab-list"] button[aria-selected="true"]{{
  background:var(--blue); color:#fff!important; border:3px solid var(--ink);
  border-bottom-color:var(--blue); box-shadow:4px 4px 0 var(--ink); }}
[data-baseweb="tab-list"] button[aria-selected="true"] p{{ color:#fff!important; font-weight:800; }}
[data-baseweb="tab-highlight"], [data-baseweb="tab-border"]{{ background:transparent!important; height:0!important; }}
[data-baseweb="tab-list"] button[data-baseweb="tab"] p{{ font-size:.82rem; font-weight:800; }}

/* expander */
[data-testid="stExpander"] details{{ border:3px solid var(--ink)!important; border-radius:2px!important;
  box-shadow:4px 4px 0 var(--ink); background:var(--surface); }}
[data-testid="stExpander"] summary{{ font-weight:800; text-transform:uppercase;
  letter-spacing:.04em; font-size:.8rem; }}

/* dataframe */
[data-testid="stDataFrame"]{{ border:3px solid var(--ink); border-radius:2px;
  box-shadow:4px 4px 0 var(--ink); }}

/* bordered containers (premium input panels) */
[data-testid="stVerticalBlockBorderWrapper"], [data-testid="stHorizontalBlockBorderWrapper"]{{
  border:3px solid var(--ink)!important; border-radius:2px!important;
  box-shadow:5px 5px 0 var(--ink)!important; background:var(--surface)!important;
}}

/* ---- similarity score tables ---- */
.nb-score-grid{{
  display:grid; grid-template-columns:repeat(auto-fit,minmax(400px,1fr));
  gap:16px; align-items:start; margin-top:4px;
}}
.nb-scorecard{{
  background:var(--surface); border:3px solid var(--ink);
  box-shadow:5px 5px 0 var(--ink); border-radius:2px;
  padding:14px 16px 7px; break-inside:avoid;
}}
.nb-scorecard-hd{{
  display:flex; justify-content:space-between; align-items:baseline;
  gap:12px; margin-bottom:2px;
}}
.nb-scorecard-title{{
  font-size:.74rem; font-weight:900; text-transform:uppercase;
  letter-spacing:.11em; color:var(--blue);
}}
.nb-scorecard-count{{
  font-size:.64rem; font-weight:800; text-transform:uppercase;
  letter-spacing:.08em; color:var(--muted); white-space:nowrap;
}}
.nb-score-row{{ border-top:1px solid var(--line); padding:9px 0 8px; }}
.nb-score-row:first-child{{ border-top:0; padding-top:7px; }}
.nb-score-top{{ display:flex; align-items:baseline; gap:12px; }}
.nb-score-name{{ flex:1; font-size:.82rem; font-weight:800; color:var(--ink); min-width:0; }}
.nb-score-val{{
  font-family:ui-monospace,Menlo,monospace; font-size:1.05rem;
  font-weight:900; color:var(--ink); white-space:nowrap;
}}
.nb-score-row.pass .nb-score-name,
.nb-score-row.pass .nb-score-val{{
  color:var(--ink); font-weight:900;
}}
.nb-score-row.fail .nb-score-name,
.nb-score-row.fail .nb-score-val,
.nb-score-row.fail .nb-score-badge{{
  color:var(--muted);
}}
.nb-score-row.fail .nb-score-name{{ font-weight:700; }}
.nb-score-badge{{
  font-size:.6rem; font-weight:800; text-transform:uppercase;
  letter-spacing:.06em; color:var(--muted); border:1.5px solid var(--line);
  border-radius:2px; padding:2px 6px; white-space:nowrap;
}}
.nb-score-track{{
  height:7px; background:var(--surface-2); border:1.5px solid var(--ink);
  margin-top:6px; position:relative; overflow:hidden;
}}
.nb-score-fill{{ height:100%; }}
.nb-score-note{{ margin-top:5px; font-size:.68rem; color:var(--muted); line-height:1.35; }}

/* file uploader */
[data-testid="stFileUploaderDropzone"]{{ border:3px dashed var(--ink)!important;
  border-radius:2px!important; background:var(--blue-ghost)!important; }}

/* alerts keep a brutalist edge */
[data-testid="stAlert"]{{ border:3px solid var(--ink); border-radius:2px;
  box-shadow:4px 4px 0 var(--ink); }}
</style>
"""


# --------------------------------------------------------------------------
# HTML-string builders
# --------------------------------------------------------------------------
def section_label_html(label: str) -> str:
    return f'<div class="nb-sec">{escape(str(label))}</div>'


def card_html(inner_html: str, ghost: bool = False) -> str:
    cls = "nb-card ghost" if ghost else "nb-card"
    return f'<div class="{cls}">{inner_html}</div>'


def stat_grid_html(rows: Sequence[Tuple[str, str]], wide: bool = False) -> str:
    cells = "".join(
        f'<div class="nb-stat"><div class="lab">{escape(str(l))}</div>'
        f'<div class="val">{escape(str(v))}</div></div>'
        for l, v in rows
    )
    cls = "nb-statgrid wide" if wide else "nb-statgrid"
    return f'<div class="{cls}">{cells}</div>'


def metric_cards_html(items: Sequence[Tuple[str, str]]) -> str:
    cells = "".join(
        f'<div class="nb-metric"><div class="lab">{escape(str(l))}</div>'
        f'<div class="val">{escape(str(v))}</div></div>'
        for l, v in items
    )
    return f'<div class="nb-metricrow">{cells}</div>'


def chips_html(items: Sequence[Tuple[str, str]]) -> str:
    parts = []
    for label, color in items:
        parts.append(
            f'<span class="nb-chip"><span class="sw" style="background:{escape(str(color))};"></span>'
            f'{escape(str(label))}</span>'
        )
    return f'<div style="margin-top:6px;">{"".join(parts)}</div>'


def _similarity_score_color(score: Optional[float]) -> str:
    """Low-saturation indigo scale for the score bar, with a muted red for negative scores."""
    if score is None:
        return BLUE_GHOST
    if score < 0:
        return "#C29393"
    if score >= 0.8:
        return "#5B78B8"
    if score >= 0.5:
        return "#7F99CC"
    if score >= 0.2:
        return "#A9BEEA"
    return "#C9D6F2"


def similarity_score_card_html(
    category: str,
    rows: Sequence[dict],
    threshold: float = 0.7,
) -> str:
    """Build one premium score card for a metric category.

    ``rows`` are dictionaries with keys ``method``, ``score``, ``matches`` and ``note``.
    Scores are shown as a filled horizontal bar; ``score=None`` is rendered as N/A.
    """
    row_html: List[str] = []
    for r in rows:
        method = escape(str(r.get("method", "")))
        score = r.get("score")
        matches = r.get("matches")
        note = escape(str(r.get("note") or ""))

        if score is None:
            value_html = '<span class="nb-score-na">N/A</span>'
            width = 0.0
        else:
            value_html = f"<span>{float(score):.4f}</span>"
            width = max(0.0, min(100.0, abs(float(score)) * 100.0))

        badge = "—"
        if matches is not None:
            badge = f"{int(matches)} matched"

        passed = score is not None and float(score) >= threshold
        row_class = "nb-score-row pass" if passed else "nb-score-row fail"
        fill_color = _similarity_score_color(score)
        note_html = f'<div class="nb-score-note">{note}</div>' if note else ""
        row_html.append(
            f'<div class="{row_class}">'
            f'<div class="nb-score-top">'
            f'<div class="nb-score-name">{method}</div>'
            f'<div class="nb-score-badge">{badge}</div>'
            f'<div class="nb-score-val">{value_html}</div>'
            '</div>'
            f'<div class="nb-score-track"><div class="nb-score-fill" '
            f'style="width:{width:.1f}%;background:{fill_color};"></div></div>'
            f'{note_html}'
            '</div>'
        )

    return (
        '<div class="nb-scorecard">'
        '<div class="nb-scorecard-hd">'
        f'<div class="nb-scorecard-title">{escape(str(category))}</div>'
        f'<div class="nb-scorecard-count">{len(rows)} methods</div>'
        '</div>'
        f'{"".join(row_html)}'
        '</div>'
    )


def similarity_score_cards_html(
    groups: Sequence[Tuple[str, Sequence[dict]]],
    threshold: float = 0.7,
) -> str:
    cards = "".join(
        similarity_score_card_html(category, rows, threshold=threshold)
        for category, rows in groups
    )
    return f'<div class="nb-score-grid">{cards}</div>'


def property_card_html(label: str, value: str, value_style: str = "") -> str:
    """A bordered card for long identifiers such as InChI / InChI Key."""
    style = "word-break:break-all;line-height:1.55;" + value_style
    return card_html(
        section_label_html(label)
        + f'<div class="mono" style="{style}">{escape(str(value))}</div>'
    )


_COPY_ICON = (
    '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" '
    'stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">'
    '<rect x="9" y="9" width="11" height="11" rx="2"></rect>'
    '<path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"></path></svg>'
)

_CHECK_ICON = (
    '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" '
    'stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">'
    '<polyline points="20 6 9 17 4 12"></polyline></svg>'
)


def copy_fields_html(
    fields: Sequence[Tuple[str, str, bool]],
    heading: Optional[str] = None,
) -> str:
    """Build an identity card whose copyable rows have small copy buttons.

    Each field is ``(label, value, copyable)``. Displayed values are
    HTML-escaped and clipboard payloads are base64-encoded, so user-supplied
    SMILES cannot inject markup or script into the rendered block.
    """
    rows: List[str] = []
    for idx, (label, value, copyable) in enumerate(fields):
        lab = escape(str(label))
        val = escape(str(value))
        if copyable:
            payload = base64.b64encode(str(value).encode("utf-8")).decode("ascii")
            control = (
                f'<button class="copy-btn" type="button" data-copy="{payload}" '
                f'title="Copy {lab}" aria-label="Copy {lab}">'
                f'<span class="icon-copy">{_COPY_ICON}</span>'
                f'<span class="icon-ok">{_CHECK_ICON}</span></button>'
            )
        else:
            control = '<span class="copy-spacer"></span>'
        row_class = "row first" if idx == 0 else "row"
        rows.append(
            f'<div class="{row_class}"><div class="lab">{lab}</div>'
            f'<div class="val">{val}</div>{control}</div>'
        )

    heading_html = (
        f'<div class="nb-sec">{escape(str(heading))}</div>' if heading else ""
    )
    script = (
        "<script>(function(){"
        "function b64ToText(s){try{"
        "var bin=atob(s);var bytes=new Uint8Array(bin.length);"
        "for(var i=0;i<bin.length;i++){bytes[i]=bin.charCodeAt(i);}"
        "if(typeof TextDecoder==='undefined'){return bin;}"
        "return new TextDecoder().decode(bytes);}catch(e){return s;}}"
        "function fallbackCopy(text){var ta=document.createElement('textarea');"
        "ta.value=text;ta.style.position='fixed';ta.style.opacity='0';"
        "document.body.appendChild(ta);ta.select();var ok=false;"
        "try{ok=document.execCommand('copy');}catch(e){}"
        "document.body.removeChild(ta);return ok;}"
        "function flash(btn){if(!btn)return;btn.classList.add('copied');"
        "window.setTimeout(function(){btn.classList.remove('copied');},1200);}"
        "function copy(text,btn){"
        "if(navigator.clipboard&&window.isSecureContext){"
        "navigator.clipboard.writeText(text).then(function(){flash(btn);},"
        "function(){if(fallbackCopy(text)){flash(btn);}});}"
        "else if(fallbackCopy(text)){flash(btn);}}"
        "document.querySelectorAll('.nb-identity .copy-btn').forEach(function(btn){"
        "btn.addEventListener('click',function(){"
        "copy(b64ToText(btn.getAttribute('data-copy')),btn);});"
        "});})();</script>"
    )
    return f'<div class="nb-identity">{heading_html}{"".join(rows)}</div>{script}'


def verdict_html(category: str, subtitle: str, accent: str) -> str:
    return (
        f'<div class="nb-verdict" style="background:{escape(str(accent))}">'
        f'<div class="cat">{escape(str(category))}</div>'
        f'<div class="sub">{escape(str(subtitle))}</div></div>'
    )


# --------------------------------------------------------------------------
# Render helpers (call st.* directly)
# --------------------------------------------------------------------------
def inject_css() -> None:
    st.markdown(GLOBAL_CSS, unsafe_allow_html=True)


def hero(title: str, subtitle: Optional[str] = None, kicker: Optional[str] = None) -> None:
    html = ""
    if kicker:
        html += f'<div class="nb-kicker">{escape(kicker)}</div>'
    html += f'<h1 class="nb-hero">{title}</h1>'  # title may contain <span class="accent">
    if subtitle:
        html += f'<div class="nb-sub">{escape(subtitle)}</div>'
    html += '<hr class="nb-rule"/>'
    st.markdown(html, unsafe_allow_html=True)


def tab_header(title: str, subtitle: Optional[str] = None, kicker: Optional[str] = None) -> None:
    html = '<div class="nb-tabhead">'
    if kicker:
        html += f'<div class="nb-kicker">{escape(kicker)}</div>'
    html += f'<div class="t">{escape(title)}</div>'
    if subtitle:
        html += f'<div class="s">{escape(subtitle)}</div>'
    html += '</div>'
    st.markdown(html, unsafe_allow_html=True)


def section(label: str) -> None:
    st.markdown(section_label_html(label), unsafe_allow_html=True)


def metric_cards(items: Sequence[Tuple[str, str]]) -> None:
    st.markdown(metric_cards_html(items), unsafe_allow_html=True)


def stat_grid(rows: Sequence[Tuple[str, str]], wide: bool = False) -> None:
    st.markdown(stat_grid_html(rows, wide=wide), unsafe_allow_html=True)


def verdict(category: str, subtitle: str, accent: str) -> None:
    st.markdown(verdict_html(category, subtitle, accent), unsafe_allow_html=True)


def chips(items: Sequence[Tuple[str, str]]) -> None:
    st.markdown(chips_html(items), unsafe_allow_html=True)


def similarity_score_cards(
    groups: Sequence[Tuple[str, Sequence[dict]]],
    threshold: float = 0.7,
) -> None:
    st.markdown(similarity_score_cards_html(groups, threshold=threshold), unsafe_allow_html=True)


def property_card(label: str, value: str, value_style: str = "") -> None:
    st.markdown(property_card_html(label, value, value_style), unsafe_allow_html=True)


def copy_fields(
    fields: Sequence[Tuple[str, str, bool]],
    heading: Optional[str] = None,
) -> None:
    st.html(copy_fields_html(fields, heading), unsafe_allow_javascript=True)


_3DMOL_JS_PATH = Path(__file__).resolve().parent / "lib" / "3Dmol-min.js"


def _load_3dmol_js() -> str:
    try:
        return _3DMOL_JS_PATH.read_text(encoding="utf-8")
    except OSError:
        return ""


def molecule_3d_html(molblock: str, height: int = 420) -> str:
    """Self-contained 3Dmol.js viewer for an MDL molblock (with toolbar)."""
    payload = base64.b64encode(molblock.encode("utf-8")).decode("ascii")
    lib_js = _load_3dmol_js()
    lib_tag = (
        f"<script>{lib_js}</script>"
        if lib_js
        else '<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>'
    )

    return f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
{lib_tag}
<style>
html, body {{ margin:0; padding:0; height:100%; background:{SURFACE}; }}
.wrap {{ display:flex; flex-direction:column; height:100%; box-sizing:border-box;
  border:3px solid {INK}; border-radius:2px; box-shadow:5px 5px 0 {INK};
  background:{SURFACE}; }}
.stage {{ flex:1; min-height:0; position:relative; background:{SURFACE}; }}
.stage canvas {{ display:block; }}
.bar {{ display:flex; flex-wrap:wrap; gap:6px; align-items:center;
  padding:8px 10px; border-top:3px solid {INK}; background:{BLUE_GHOST}; }}
.vbtn {{ font-family:ui-monospace,"SF Mono",Menlo,monospace; font-size:11px;
  font-weight:800; background:{SURFACE}; color:{INK}; border:2px solid {INK};
  border-radius:2px; padding:4px 9px; cursor:pointer; box-shadow:2px 2px 0 {INK}; }}
.vbtn:hover {{ background:{BLUE_SOFT}; }}
.vbtn.active {{ background:{INK}; color:#FFFFFF; }}
.spacer {{ flex:1; }}
.hint {{ font-size:10px; color:{MUTED}; text-transform:uppercase; letter-spacing:.08em; }}
</style>
</head>
<body>
<div class="wrap">
  <div id="viewer" class="stage"></div>
  <div class="bar">
    <button class="vbtn active" data-style="stick" type="button">Stick</button>
    <button class="vbtn" data-style="ballstick" type="button">Ball &amp; Stick</button>
    <button class="vbtn" data-style="sphere" type="button">Spacefill</button>
    <button class="vbtn" data-style="line" type="button">Line</button>
    <span class="spacer"></span>
    <span class="hint">drag · scroll</span>
    <button class="vbtn" id="spin" type="button">Spin</button>
    <button class="vbtn" id="reset" type="button">Reset</button>
  </div>
</div>
<script>
(function(){{
  var molblock = atob("{payload}");
  if (typeof $3Dmol === "undefined") {{
    document.getElementById("viewer").innerHTML = "3D viewer failed to load.";
    return;
  }}
  var viewer = $3Dmol.createViewer(document.getElementById("viewer"), {{
    backgroundColor: "{SURFACE}"
  }});
  viewer.addModel(molblock, "sdf");
  var style = {{ stick: {{ radius: 0.15 }} }};
  function applyStyle(name) {{
    if (name === "ballstick") {{
      style = {{ stick: {{ radius: 0.12 }}, sphere: {{ scale: 0.25 }} }};
    }} else if (name === "sphere") {{
      style = {{ sphere: {{ scale: 0.42 }} }};
    }} else if (name === "line") {{
      style = {{ line: {{ linewidth: 2 }} }};
    }} else {{
      style = {{ stick: {{ radius: 0.15 }} }};
    }}
    viewer.setStyle({{}}, style);
    viewer.render();
  }}
  viewer.setStyle({{}}, style);
  viewer.zoomTo();
  viewer.render();
  window.addEventListener("resize", function() {{
    viewer.resize();
    viewer.render();
  }});
  var spinning = false;
  document.querySelectorAll("[data-style]").forEach(function(btn) {{
    btn.addEventListener("click", function() {{
      document.querySelectorAll("[data-style]").forEach(function(x) {{
        x.classList.remove("active");
      }});
      btn.classList.add("active");
      applyStyle(btn.dataset.style);
    }});
  }});
  document.getElementById("spin").addEventListener("click", function() {{
    spinning = !spinning;
    this.classList.toggle("active", spinning);
    viewer.spin(spinning ? "y" : false, 1);
  }});
  document.getElementById("reset").addEventListener("click", function() {{
    viewer.zoomTo();
    viewer.render();
  }});
}})();
</script>
</body>
</html>"""


def molecule_3d(molblock: str, height: int = 420) -> None:
    components.html(
        molecule_3d_html(molblock, height),
        height=height,
        scrolling=False,
    )


def divider(hairline: bool = True) -> None:
    st.markdown('<hr class="nb-hairline"/>' if hairline else '<hr class="nb-rule"/>',
                unsafe_allow_html=True)


def note(text: str) -> None:
    st.markdown(f'<div class="nb-note">{escape(text)}</div>', unsafe_allow_html=True)


def under_construction(title: str, note: Optional[str] = None, eta: Optional[str] = None) -> None:
    html = ('<div class="nb-uc"><span class="tag">Under Construction</span>'
            f'<div class="t">{escape(title)}</div>')
    if note:
        html += f'<div class="n">{escape(note)}</div>'
    if eta:
        html += f'<div class="eta">{escape(eta)}</div>'
    html += '</div>'
    st.markdown(html, unsafe_allow_html=True)
