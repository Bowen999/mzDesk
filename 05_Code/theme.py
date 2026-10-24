"""theme.py — unified design system for the Metabolomics Workbench app.

Design language: NEO-BRUTALISM x MINIMALISM, cool-white + indigo ("high-end blue").
  - thick navy-ink borders, hard offset shadows (no blur), near-square corners
  - uppercase mono kickers / section labels, restrained indigo-blue palette
  - generous whitespace, hairline dividers, few decorative elements

Two things sit on top of that look and are deliberately not decoration:

  TWO HEADING LEVELS. ``section()`` is a page section (indigo square + label);
  ``subsection()`` is a block inside one (no square, hairline rule, lighter).
  Using one level for both is what made a card title read with the same weight
  as a page section.

  MOTION. Switching tabs, and a result arriving after a run, are the two moments
  a reader can lose their place, so both are animated. Nothing else moves — the
  button press stays instant — and everything honours ``prefers-reduced-motion``.

PUBLIC API (import as `import theme`):

  Colours (hex strings):
    theme.INK, BG, SURFACE, SURFACE_2, MUTED, LINE,
    theme.BLUE, BLUE_BRIGHT, BLUE_SOFT, BLUE_GHOST
    theme.CATEGORY_ACCENT  -> {"Lipidome","Metabolome","Both","Neither","Invalid"}
    theme.PLOT_DARK, PLOT_LIGHT     (spectrum stick colours: annotated / unannotated)
    theme.PLOT_FONT, PLOT_MONO      (font stacks shared by the spectrum charts)
    theme.FRAG_HL                   (RDKit substructure-highlight RGB 0..1 tuples)

  Render helpers (call st.* directly — use inside a tab / page):
    theme.inject_css()                         # once, right after set_page_config
    theme.hero(title, subtitle=None, kicker=None)
    theme.tab_header(title, subtitle=None, kicker=None)
    theme.section(label, aux=None)             # level 1: square marker + label
    theme.subsection(label, aux=None)          # level 2: block title inside a section
    with theme.panel(key): ...                 # bordered input panel (keyed st.container)
    theme.metric_cards(items)                  # items: [(label, value), ...]
    theme.stat_grid(rows, wide=False, max_cols=None, min_px=160)  # rows: [(label, value), ...]
    theme.verdict(category, subtitle, accent)  # big result banner
    theme.chips(items)                         # items: [(label, css_color), ...]
    theme.chip_line(label, items)              # one-line status row: mono label + chips
    theme.reasons(items)                       # bulleted conclusions, balanced columns
    theme.notes(items)                         # muted paragraphs in balanced columns
    theme.similarity_score_cards(groups)       # premium score cards for grouped metrics
    theme.copy_fields(fields, heading=None, variant=None)  # identity / adduct rows
    theme.molecule_3d(molblock, height=420, title=None, caption=None)  # 3Dmol.js viewer
    theme.formula_verdict(formula, subtitle, accent)  # verdict that keeps formula case
    theme.under_construction(title, note=None, eta=None)
    theme.divider()
    theme.note(text)

  HTML-string builders (return str; caller wraps in st.markdown(..., unsafe_allow_html=True)):
    theme.card_html, section_label_html, subsection_label_html, stat_grid_html,
    theme.metric_cards_html, chips_html, chip_line_html, reasons_html, notes_html,
    similarity_score_cards_html, verdict_html
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
    "Neither":    "#6B7588",  # cool grey (dark enough for white verdict text, 4.6:1)
    "Invalid":    "#C0392B",  # functional error red
}

# Spectrum stick colours. PLOT_LIGHT is a shade darker than BLUE_SOFT so the 2-px sticks of
# unannotated / unmatched peaks stay visible on white (~2.6:1) while ranking below PLOT_DARK.
PLOT_DARK  = BLUE          # annotated / explained / matched peaks
PLOT_LIGHT = "#8A9FD0"     # unannotated / unexplained / unmatched peaks
# Font stacks shared by the three spectrum charts. Helvetica first: the MS² chart renders in an
# iframe where Streamlit's bundled Source Sans is not loaded, and the charts should all match.
PLOT_FONT  = '"Helvetica Neue",Helvetica,Arial,sans-serif'
PLOT_MONO  = 'ui-monospace,"SF Mono",Menlo,Consolas,monospace'

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
/* The header is transparent and its toolbar hidden: let clicks reach the content scrolled under it. */
[data-testid="stHeader"]{{ background:transparent; pointer-events:none; }}
#MainMenu, footer, [data-testid="stToolbar"], [data-testid="stDecoration"]{{ visibility:hidden; }}
/* ...but keep Streamlit's run indicator (running icon + Stop) visible as the global busy signal. */
[data-testid="stStatusWidget"]{{
  visibility:visible; pointer-events:auto; background:var(--surface);
  border:2px solid var(--ink); border-radius:2px; box-shadow:3px 3px 0 var(--ink);
  padding:0 .45rem;
}}
.block-container{{ padding-top:1.4rem; padding-bottom:3rem; max-width:1240px; }}
/* Streamlit's bundled Source Sans (config font = sans-serif): inherited text matches markdown text. */
html, body{{
  font-family:"Source Sans","Source Sans Pro",-apple-system,BlinkMacSystemFont,"Segoe UI",
    "Helvetica Neue",Arial,sans-serif;
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
/* scoped under the markdown container so it beats Streamlit's own h1 size/weight/padding */
.nb-hero, [data-testid="stMarkdownContainer"] h1.nb-hero{{
  text-transform:uppercase; letter-spacing:.005em; font-weight:900;
  font-size:clamp(2.1rem,4vw,3rem); line-height:1.0; margin:0; padding:0; color:var(--ink);
}}
.nb-hero [data-testid="stHeaderActionElements"]{{ display:none; }}
.nb-hero .accent{{ color:var(--blue); }}
.nb-sub{{ font-size:.95rem; color:var(--muted); margin:.55rem 0 0 0; max-width:70ch; }}
/* scoped likewise: Streamlit's markdown <hr> style (1px line, 2em margins) otherwise wins */
hr.nb-rule, [data-testid="stMarkdownContainer"] hr.nb-rule{{
  height:3px; background:var(--ink); border:0; margin:1.1rem 0 1.4rem 0; }}
hr.nb-hairline, [data-testid="stMarkdownContainer"] hr.nb-hairline{{
  height:1px; background:var(--line); border:0; margin:1.3rem 0; }}

.nb-tabhead{{ margin:.2rem 0 1.1rem 0; }}
.nb-tabhead .t{{ font-size:1.55rem; font-weight:900; text-transform:uppercase;
  letter-spacing:.01em; line-height:1.05; margin:0; }}
.nb-tabhead .s{{ font-size:.9rem; color:var(--muted); margin:.35rem 0 0 0; }}

/* ---- section label: LEVEL 1 (a page section) ---- */
.nb-sec{{
  display:flex; align-items:center; gap:.55rem;
  font-size:.72rem; text-transform:uppercase; letter-spacing:.14em;
  font-weight:800; color:var(--muted); margin:.2rem 0 .7rem 0;
}}
.nb-sec::before{{ content:""; width:11px; height:11px; background:var(--blue);
  border:2px solid var(--ink); display:inline-block; flex:none; }}
.nb-sec .aux{{
  margin-left:auto; font-family:ui-monospace,"SF Mono",Menlo,monospace;
  font-size:.64rem; letter-spacing:.08em; color:var(--muted);
  font-weight:700; white-space:nowrap;
}}

/* ---- section label: LEVEL 2 (a block inside a section) ----
   No square marker, and a hairline under it instead: it has to read as subordinate
   to .nb-sec, or a card title carries the same weight as a page section. */
.nb-subsec{{
  display:flex; align-items:baseline; gap:.5rem;
  font-size:.65rem; text-transform:uppercase; letter-spacing:.12em;
  font-weight:800; color:var(--muted);
  border-bottom:1px solid var(--line); padding-bottom:5px; margin:.85rem 0 .6rem 0;
}}
.nb-subsec .aux{{
  margin-left:auto; font-family:ui-monospace,"SF Mono",Menlo,monospace;
  font-size:.6rem; letter-spacing:.08em; color:var(--muted);
  font-weight:700; white-space:nowrap;
}}

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
.nb-statgrid{{ display:grid; grid-template-columns:repeat(2,1fr); gap:12px;
  /* every row takes the height of the tallest one: a card whose label wraps
     ("Isoprenoid / sterol") otherwise leaves its neighbour short, and a grid of
     ragged rows reads as misaligned even though every card is the same width. */
  grid-auto-rows:1fr; }}
.nb-statgrid.wide{{ grid-template-columns:repeat(auto-fit,minmax(160px,1fr)); }}
.nb-stat{{ background:var(--surface); border:3px solid var(--ink);
  box-shadow:4px 4px 0 var(--ink); border-radius:2px; padding:10px 13px; }}
.nb-stat .lab{{ font-size:.7rem; text-transform:uppercase; letter-spacing:.08em;
  font-weight:800; color:var(--muted); line-height:1.25; }}
.nb-stat .val{{ font-size:1.4rem; font-weight:800; color:var(--ink);
  font-family:ui-monospace,Menlo,monospace; margin-top:3px; word-break:break-word; }}

/* ---- metric cards (KPI row) ---- */
.nb-metricrow{{ display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr));
  gap:12px; margin-bottom:6px; }}
.nb-metric{{ background:var(--ink); color:#fff; border:3px solid var(--ink);
  box-shadow:5px 5px 0 var(--blue); border-radius:2px; padding:13px 15px; }}
.nb-metric .lab{{ font-size:.7rem; text-transform:uppercase; letter-spacing:.1em;
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
  padding:14px 16px 10px; margin-bottom:16px; container-type:inline-size;
}}
.nb-identity .nb-subsec{{ margin:.1rem 0 .55rem 0; }}
.nb-identity .row{{
  display:flex; align-items:flex-start; gap:12px;
  border-top:1px solid var(--line); padding:8px 0;
}}
.nb-identity .row.first{{ border-top:0; padding-top:2px; }}
/* scoped to .row: the block heading is a label too, and must not take the row width */
.nb-identity .row .lab{{
  flex:0 0 124px; font-size:.64rem; text-transform:uppercase;
  letter-spacing:.08em; font-weight:800; color:var(--muted); padding-top:3px;
}}
.nb-identity .row .val{{
  flex:1; font-family:ui-monospace,"SF Mono",Menlo,Consolas,monospace;
  font-size:.76rem; line-height:1.55; color:var(--ink);
  word-break:break-all; min-width:0;
}}
.nb-identity .row .aux{{
  flex:none; font-family:ui-monospace,"SF Mono",Menlo,Consolas,monospace;
  font-size:.68rem; color:var(--muted); white-space:nowrap;
}}
/* copy button with a CSS-drawn icon (st.html's sanitizer strips inline <svg>) */
.nb-identity .copy-btn{{
  position:relative; flex:none; width:20px; height:20px; margin-top:-1px; padding:0;
  background:var(--blue-ghost); border:1.5px solid var(--ink); border-radius:2px;
  box-shadow:2px 2px 0 var(--ink); cursor:pointer; color:var(--blue);
}}
.nb-identity .copy-btn::before, .nb-identity .copy-btn::after{{
  content:""; position:absolute; box-sizing:border-box; width:8px; height:9px;
  border:1.5px solid currentColor; border-radius:1px;
}}
.nb-identity .copy-btn::before{{ left:3px; top:2px; }}
.nb-identity .copy-btn::after{{ left:6px; top:5px; background:var(--blue-ghost); }}
.nb-identity .copy-btn:hover, .nb-identity .copy-btn:hover::after{{ background:var(--blue-soft); }}
.nb-identity .copy-btn:active{{ transform:translate(2px,2px); box-shadow:none; }}
.nb-identity .copy-btn.copied{{ color:#1C7A3D; }}
.nb-identity .copy-btn.copied::before{{ display:none; }}
.nb-identity .copy-btn.copied::after{{
  left:6px; top:2px; width:5px; height:9px; background:none;
  border-width:0 2px 2px 0; border-radius:0; transform:rotate(45deg);
}}
.nb-identity .copy-spacer{{ flex:none; width:20px; }}
/* Adduct tables: ion notation is data, so mono and not upper-cased. These rows are a grid
   rather than flex — an m/z is one unbreakable token, so in a narrow card a flex row lets it
   overflow on top of the delta note instead of shrinking. Grid tracks cannot collide. */
.nb-identity.adducts .row{{
  display:grid; grid-template-columns:minmax(0,auto) 1fr auto auto;
  align-items:center; gap:0 10px; padding:6px 0;
}}
.nb-identity.adducts .row.first{{ padding-top:2px; }}
.nb-identity.adducts .row .lab{{
  flex:none; padding-top:0; font-family:ui-monospace,"SF Mono",Menlo,Consolas,monospace;
  font-size:.74rem; text-transform:none; letter-spacing:0; color:var(--ink);
  white-space:nowrap; overflow:hidden; text-overflow:ellipsis;
}}
.nb-identity.adducts .row .val{{
  font-size:.82rem; font-weight:800; text-align:right;
  word-break:normal; white-space:nowrap;
}}
.nb-identity.adducts .row .aux{{ min-width:0; overflow:hidden; }}
/* a narrow adduct card drops the delta note rather than crushing the m/z */
@container (max-width:250px){{
  .nb-identity.adducts .row .aux{{ display:none; }}
}}
/* narrow card: identifier labels go above their values */
@container (max-width:340px){{
  .nb-identity:not(.adducts) .row{{ flex-wrap:wrap; row-gap:2px; }}
  .nb-identity:not(.adducts) .row .lab{{ flex-basis:100%; }}
}}

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

/* ---- running text in columns ----
   A full-width section whose body is one narrow column of text leaves more than
   half the canvas blank. CSS multi-column lets the browser balance the flow
   instead of the caller guessing where to split it; both collapse to one column
   before the text gets too narrow to read. */
.nb-reasons{{ margin:.15rem 0 .3rem 0; padding-left:1.15rem; column-gap:38px; }}
.nb-reasons li{{ font-size:.87rem; line-height:1.5; color:var(--ink);
  margin:0 0 .42rem 0; break-inside:avoid; }}
.nb-reasons li::marker{{ color:var(--blue); }}
.nb-reasons li.last{{ font-weight:800; }}
.nb-notegrid{{ column-gap:38px; }}
.nb-notegrid .nb-note{{ break-inside:avoid; margin:0 0 .7rem 0; }}
@media (max-width:820px){{
  .nb-reasons, .nb-notegrid{{ columns:1!important; }}
}}

/* a chips row with a leading mono label: sits on one line with the action button
   beside it, without the rule a subsection label would draw across the row */
.nb-chipline{{ display:flex; align-items:center; gap:.65rem; flex-wrap:wrap; }}
.nb-chipline .lb{{ font-size:.65rem; text-transform:uppercase; letter-spacing:.12em;
  font-weight:800; color:var(--muted); white-space:nowrap; }}
.nb-chipline .chips{{ display:flex; align-items:center; flex-wrap:wrap; gap:8px; }}
.nb-chipline .nb-chip{{ margin:0; }}

/* ---- Streamlit widgets ---- */
/* buttons: primary = filled indigo CTA, secondary = white outline; both keep the ink frame and hard
   shadow and press into it on :active. data-testid selectors only: the old `.stButton>button` rule
   out-ranked the primary rule, so every button rendered with the same style. */
[data-testid="stBaseButton-primary"], [data-testid="stBaseButton-secondary"]{{
  border:3px solid var(--ink); border-radius:2px; box-shadow:4px 4px 0 var(--ink);
  padding:.5rem 1.3rem; transition:background-color .12s ease, color .12s ease; }}
[data-testid="stBaseButton-primary"] p, [data-testid="stBaseButton-secondary"] p{{
  font-weight:800; text-transform:uppercase; letter-spacing:.05em; font-size:.88rem; }}
[data-testid="stBaseButton-primary"]{{ background:var(--blue); color:#fff; }}
[data-testid="stBaseButton-primary"]:hover{{ background:var(--blue-bright); color:#fff; border-color:var(--ink); }}
[data-testid="stBaseButton-secondary"]{{ background:var(--surface); color:var(--ink); }}
[data-testid="stBaseButton-secondary"]:hover{{ background:var(--blue-ghost); color:var(--blue); border-color:var(--ink); }}
/* the press stays instant on purpose: the block slams into its own shadow */
[data-testid="stBaseButton-primary"]:active, [data-testid="stBaseButton-secondary"]:active{{
  transform:translate(4px,4px); box-shadow:none; transition:none; }}
[data-testid="stBaseButton-primary"]:focus-visible, [data-testid="stBaseButton-secondary"]:focus-visible{{
  outline:3px solid var(--blue-bright); outline-offset:3px; }}
[data-testid="stBaseButton-primary"]:disabled, [data-testid="stBaseButton-secondary"]:disabled{{
  background:var(--surface-2); color:var(--muted); border-color:var(--muted);
  box-shadow:none; transform:none; cursor:not-allowed; }}

/* text / number inputs + textarea: the ink frame sits on the BaseWeb root, so the hard shadow is
   not clipped and the number-input steppers sit inside the frame */
[data-testid="stTextInputRootElement"], [data-testid="stTextAreaRootElement"],
[data-testid="stNumberInputContainer"]{{
  border:3px solid var(--ink)!important; border-radius:2px!important;
  box-shadow:3px 3px 0 var(--ink); background:var(--surface)!important;
  transition:border-color .12s ease; }}
[data-testid="stTextInputRootElement"]:focus-within, [data-testid="stTextAreaRootElement"]:focus-within,
[data-testid="stNumberInputContainer"]:focus-within{{ border-color:var(--blue-bright)!important; }}
[data-testid="stTextInputRootElement"] [data-baseweb="base-input"],
[data-testid="stTextAreaRootElement"] [data-baseweb="base-input"],
[data-testid="stNumberInputContainer"] [data-baseweb="input"],
[data-testid="stNumberInputContainer"] [data-baseweb="base-input"]{{
  border:0!important; background:var(--surface)!important; }}
[data-testid="stTextInput"] input, [data-testid="stNumberInput"] input,
[data-testid="stTextArea"] textarea{{
  background:var(--surface)!important; font-family:ui-monospace,Menlo,monospace;
  color:var(--ink)!important; }}
/* Whole-pixel line box so a fixed height never clips a peak-list row through the middle.
   st.text_area(height=H) gives the textarea H minus Streamlit's label strip; with this line
   box the peak-list heights in the tabs are whole numbers of rows. See README. */
[data-testid="stTextArea"] textarea{{
  line-height:20px; padding-top:9px; padding-bottom:9px; }}
div[data-baseweb="select"]>div{{
  border:3px solid var(--ink)!important; border-radius:2px!important;
  box-shadow:3px 3px 0 var(--ink); background:var(--surface)!important;
  transition:border-color .12s ease; }}
div[data-baseweb="select"]>div:focus-within{{ border-color:var(--blue-bright)!important; }}
/* selectbox dropdown menu */
[data-testid="stSelectboxVirtualDropdown"]{{ border:3px solid var(--ink)!important;
  background:var(--surface)!important; }}
[data-baseweb="popover"]:has([data-testid="stSelectboxVirtualDropdown"]){{
  box-shadow:4px 4px 0 var(--ink)!important; }}
label, [data-testid="stWidgetLabel"] p{{ font-weight:700!important; color:var(--ink)!important;
  font-size:.82rem!important; }}

/* number-input steppers */
[data-testid="stNumberInputStepUp"], [data-testid="stNumberInputStepDown"]{{
  border-left:2px solid var(--ink)!important; }}

/* checkbox box */
[data-testid="stCheckbox"] label[data-baseweb="checkbox"] > span:first-child{{
  border:2px solid var(--ink)!important; border-radius:2px!important; }}

/* segmented control: square ink segments, active = indigo fill */
[data-testid="stBaseButton-segmented_control"], [data-testid="stBaseButton-segmented_controlActive"]{{
  border:2px solid var(--ink)!important; border-radius:0!important;
  transition:background-color .12s ease, color .12s ease; }}
[data-testid="stButtonGroup"] [data-baseweb="button-group"] > button + button{{ margin-left:-2px; }}
[data-testid="stBaseButton-segmented_control"]{{ background:var(--surface)!important; color:var(--ink)!important; }}
[data-testid="stBaseButton-segmented_control"]:hover{{ background:var(--blue-ghost)!important;
  color:var(--blue)!important; }}
[data-testid="stBaseButton-segmented_controlActive"]{{ background:var(--blue)!important; color:#fff!important; }}
[data-testid="stButtonGroup"] [data-baseweb="button-group"] button p{{ font-weight:700; }}

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
  letter-spacing:.04em; font-size:.82rem; color:var(--muted);
  transition:background-color .14s ease, color .14s ease; }}
/* hover restyles inactive tabs only: the active tab's white label must never sit on a pale fill */
[data-baseweb="tab-list"] button[data-baseweb="tab"]:not([aria-selected="true"]):hover{{
  color:var(--ink); background:var(--blue-ghost); }}
[data-baseweb="tab-list"] button[aria-selected="true"]{{
  background:var(--blue); color:#fff!important; border:3px solid var(--ink);
  border-bottom-color:var(--blue); box-shadow:4px 4px 0 var(--ink); }}
[data-baseweb="tab-list"] button[aria-selected="true"] p{{ color:#fff!important; font-weight:800; }}
[data-baseweb="tab-list"] button[data-baseweb="tab"]:focus-visible{{
  outline:2px dashed currentColor; outline-offset:-6px; }}
[data-baseweb="tab-highlight"], [data-baseweb="tab-border"]{{ background:transparent!important; height:0!important; }}
[data-baseweb="tab-list"] button[data-baseweb="tab"] p{{ font-size:.82rem; font-weight:800; }}

/* expander (the label is a markdown <p> that resets size/weight, so style it too) */
[data-testid="stExpander"] details{{ border:3px solid var(--ink)!important; border-radius:2px!important;
  box-shadow:4px 4px 0 var(--ink); background:var(--surface); }}
[data-testid="stExpander"] summary, [data-testid="stExpander"] summary p{{ font-weight:800;
  text-transform:uppercase; letter-spacing:.04em; font-size:.8rem; }}
[data-testid="stExpander"] summary{{ transition:background-color .12s ease; }}
[data-testid="stExpander"] summary:hover{{ background:var(--blue-ghost); }}

/* dataframe */
[data-testid="stDataFrame"]{{ border:3px solid var(--ink); border-radius:2px;
  box-shadow:4px 4px 0 var(--ink); }}

/* bordered panels: theme.panel(key) == st.container(border=True, key="nb-panel-<key>").
   Streamlit 1.5x draws container borders without a DOM hook of their own (the old
   stVerticalBlockBorderWrapper test id no longer exists), so the key class is the hook. */
[class*="st-key-nb-panel-"]{{
  border:3px solid var(--ink)!important; border-radius:2px!important;
  box-shadow:5px 5px 0 var(--ink)!important; background:var(--surface)!important;
}}

/* ---- similarity score tables ---- */
.nb-score-grid{{
  /* masonry-style columns: short cards stack beside tall ones; 1 column when narrow */
  columns:2 360px; column-gap:22px; margin-top:4px;
}}
.nb-scorecard{{
  background:var(--surface); border:3px solid var(--ink);
  box-shadow:5px 5px 0 var(--ink); border-radius:2px;
  padding:14px 16px 7px; break-inside:avoid; margin:0 0 20px 0;
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
  font-size:.64rem; font-weight:800; text-transform:uppercase;
  letter-spacing:.06em; color:var(--muted); border:1.5px solid var(--line);
  border-radius:2px; padding:2px 6px; white-space:nowrap;
}}
.nb-score-track{{
  height:7px; background:var(--surface-2); border:1.5px solid var(--ink);
  margin-top:6px; position:relative; overflow:hidden;
}}
.nb-score-fill{{ height:100%; }}
/* tick at the score threshold: on a near-full bar 0.99 and 1.00 are otherwise identical */
.nb-score-tick{{
  position:absolute; top:0; bottom:0; width:1.5px; background:var(--ink); opacity:.5;
}}
.nb-score-note{{ margin-top:5px; font-size:.68rem; color:var(--muted); line-height:1.35; }}
.nb-score-row.na .nb-score-name{{ color:var(--muted); font-weight:700; }}
.nb-score-na{{ font-size:.82rem; letter-spacing:.06em; color:var(--muted); }}
.nb-score-row.na .nb-score-track{{ border:1.5px dashed var(--line); background:transparent; }}
.nb-score-tag{{
  display:inline-block; margin-left:7px; padding:1px 5px; vertical-align:1px;
  font-size:.6rem; font-weight:900; text-transform:uppercase; letter-spacing:.08em;
  color:#fff; background:var(--blue); border:1.5px solid var(--ink); border-radius:2px;
}}

/* file uploader */
[data-testid="stFileUploaderDropzone"]{{ border:3px dashed var(--ink)!important;
  border-radius:2px!important; background:var(--blue-ghost)!important; }}

/* alerts keep a brutalist edge; the tinted inner box is squared off so no corner notches show */
[data-testid="stAlert"]{{ border:3px solid var(--ink); border-radius:2px;
  box-shadow:4px 4px 0 var(--ink); background:var(--surface); }}
[data-testid="stAlertContainer"]{{ border-radius:0; }}

/* ============================================================
   MOTION
   Two moments lose the reader: switching tabs (the whole canvas changes at
   once) and a result arriving after a run (new blocks appear, often below the
   fold). Both are animated. Nothing else moves.
   ============================================================ */
@keyframes nb-rise{{
  from{{ opacity:0; transform:translateY(10px); }}
  to{{ opacity:1; transform:none; }}
}}
@keyframes nb-fade{{ from{{ opacity:0; }} to{{ opacity:1; }} }}
/* one pass of the bright accent through a warning's shadow, so a stale-result
   banner is noticed without it blinking forever */
@keyframes nb-flag{{
  0%{{ box-shadow:4px 4px 0 var(--ink); }}
  45%{{ box-shadow:4px 4px 0 var(--blue-bright); }}
  100%{{ box-shadow:4px 4px 0 var(--ink); }}
}}

/* A tab panel that has just been revealed. st.tabs switches on the client, so the
   panel is unhidden rather than rebuilt and a plain CSS animation would never
   replay — MOTION_JS re-adds this class on each switch to retrigger it. */
[data-baseweb="tab-panel"].nb-enter [data-testid="stVerticalBlock"]:first-of-type{{
  animation:nb-rise .30s cubic-bezier(.22,.7,.3,1) both;
}}

/* Result blocks. Streamlit reuses a DOM node whose output is unchanged, so these
   fire when a result is genuinely new and stay still on an unrelated rerun. */
.nb-verdict{{ animation:nb-rise .34s cubic-bezier(.22,.7,.3,1) both; }}
.nb-metricrow{{ animation:nb-rise .34s cubic-bezier(.22,.7,.3,1) .04s both; }}
.nb-uc{{ animation:nb-fade .3s ease both; }}
/* the structure figure crossfades when the 2D / 3D view is switched */
[data-testid="stImage"] img{{ animation:nb-fade .22s ease both; }}
/* an alert arriving (stale result, engine skipped, parse error) has to be noticed */
[data-testid="stAlert"]{{ animation:nb-rise .26s cubic-bezier(.22,.7,.3,1) both; }}
[data-testid="stAlert"]:has([data-testid="stAlertContentWarning"]){{
  animation:nb-rise .26s cubic-bezier(.22,.7,.3,1) both, nb-flag .9s ease .28s 1;
}}

@media (prefers-reduced-motion: reduce){{
  *, *::before, *::after{{
    animation-duration:.001ms!important; animation-iteration-count:1!important;
    transition-duration:.001ms!important;
  }}
}}

/* ---- narrow screens: keep all five tabs on one row ---- */
@media (max-width: 900px){{
  [data-baseweb="tab-list"]{{ gap:4px; }}
  [data-baseweb="tab-list"] button[data-baseweb="tab"]{{ padding:8px 10px; letter-spacing:.02em; }}
  [data-baseweb="tab-list"] button[data-baseweb="tab"] p{{ font-size:.74rem; }}
}}
@media (max-width: 640px){{
  [data-baseweb="tab-list"] button[data-baseweb="tab"]{{ padding:7px 7px; letter-spacing:0; }}
  [data-baseweb="tab-list"] button[data-baseweb="tab"] p{{ font-size:.68rem; }}
}}
</style>
"""

# Retriggers the tab-panel entry animation (see the MOTION block above).
# Keep "<" out of the script body: st.html's DOMPurify pass silently drops any
# <script> whose text contains "<" followed by a letter, "/" or "!".
MOTION_JS = """<script>(function(){
  if (window.__nbMotion) { return; }
  window.__nbMotion = 1;
  function markActive(){
    var found = false;
    document.querySelectorAll('[data-baseweb="tab-panel"]').forEach(function(p){
      p.classList.remove('nb-enter');
      if (!found && p.getBoundingClientRect().height > 0) {
        found = true;
        void p.offsetWidth;
        p.classList.add('nb-enter');
      }
    });
  }
  document.addEventListener('click', function(e){
    var t = e.target;
    if (t && t.closest && t.closest('[data-baseweb="tab"]')) {
      window.setTimeout(markActive, 0);
    }
  }, true);
})();</script>"""


# --------------------------------------------------------------------------
# HTML-string builders
# --------------------------------------------------------------------------
def section_label_html(label: str, aux: Optional[str] = None) -> str:
    """Level 1: a page section — indigo square marker, uppercase label."""
    aux_html = f'<span class="aux">{escape(str(aux))}</span>' if aux else ""
    return f'<div class="nb-sec">{escape(str(label))}{aux_html}</div>'


def subsection_label_html(label: str, aux: Optional[str] = None) -> str:
    """Level 2: a block inside a section — no marker, hairline rule, lighter."""
    aux_html = f'<span class="aux">{escape(str(aux))}</span>' if aux else ""
    return f'<div class="nb-subsec">{escape(str(label))}{aux_html}</div>'


def card_html(inner_html: str, ghost: bool = False) -> str:
    cls = "nb-card ghost" if ghost else "nb-card"
    return f'<div class="{cls}">{inner_html}</div>'


def stat_grid_html(rows: Sequence[Tuple[str, str]], wide: bool = False,
                   max_cols: Optional[int] = None, min_px: int = 160) -> str:
    """Stat cards. ``wide`` auto-fits columns; ``max_cols`` caps a wide row (ten cards
    lay out 5 + 5 instead of 6 + 4) while still wrapping to fewer columns when narrow.

    ``min_px`` is the narrowest a card may get before the grid drops a column. Lower it
    to fit more cards on one row: the Basic results canvas is ~758 px, which fits five
    160 px cards only if the floor comes down to ~132 px.
    """
    cells = "".join(
        f'<div class="nb-stat"><div class="lab">{escape(str(l))}</div>'
        f'<div class="val">{escape(str(v))}</div></div>'
        for l, v in rows
    )
    cls = "nb-statgrid wide" if wide else "nb-statgrid"
    style = ""
    if wide:
        floor = max(96, int(min_px))
        track = (f"max({floor}px,calc({100 / max_cols:g}% - 16px))" if max_cols
                 else f"{floor}px")
        style = f' style="grid-template-columns:repeat(auto-fit,minmax({track},1fr));"'
    return f'<div class="{cls}"{style}>{cells}</div>'


def metric_cards_html(items: Sequence[Tuple[str, str]]) -> str:
    cells = "".join(
        f'<div class="nb-metric"><div class="lab">{escape(str(l))}</div>'
        f'<div class="val">{escape(str(v))}</div></div>'
        for l, v in items
    )
    return f'<div class="nb-metricrow">{cells}</div>'


def chips_html(items: Sequence[Tuple[str, str]], inline: bool = False) -> str:
    parts = []
    for label, color in items:
        parts.append(
            f'<span class="nb-chip"><span class="sw" style="background:{escape(str(color))};"></span>'
            f'{escape(str(label))}</span>'
        )
    body = "".join(parts)
    if inline:
        return f'<div class="chips">{body}</div>'
    return f'<div style="margin-top:6px;">{body}</div>'


def chip_line_html(label: str, items: Sequence[Tuple[str, str]]) -> str:
    """A chips row with a leading mono label, for a single-line status row."""
    return (f'<div class="nb-chipline"><span class="lb">{escape(str(label))}</span>'
            f'{chips_html(items, inline=True)}</div>')


def reasons_html(items: Sequence[str], columns: int = 2) -> str:
    """A conclusion list that fills the width; the last entry is the verdict.

    Rendered as one CSS multi-column list so the browser balances the columns.
    """
    lis = "".join(
        f'<li class="last">{escape(str(t))}</li>' if i == len(items) - 1
        else f'<li>{escape(str(t))}</li>'
        for i, t in enumerate(items)
    )
    style = f' style="columns:{int(columns)};"' if columns > 1 else ""
    return f'<ul class="nb-reasons"{style}>{lis}</ul>'


def notes_html(items: Sequence[str], columns: int = 2) -> str:
    """Muted paragraphs in balanced columns (a full-width 'how it works' block)."""
    style = f' style="columns:{int(columns)};"' if columns > 1 else ""
    body = "".join(f'<p class="nb-note">{escape(str(t))}</p>' for t in items)
    return f'<div class="nb-notegrid"{style}>{body}</div>'


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

    ``rows`` are dictionaries with keys ``method``, ``score``, ``matches``, ``note`` and
    optionally ``primary`` (adds a "Primary" tag). Scores are shown as a filled horizontal
    bar carrying a tick at ``threshold`` — without it a 0.99 and a 1.00 bar look identical.
    ``score=None`` (or a non-finite score) is rendered as N/A.
    """
    from math import isfinite

    tick = max(0.0, min(100.0, float(threshold) * 100.0))
    row_html: List[str] = []
    for r in rows:
        method = escape(str(r.get("method", "")))
        score = r.get("score")
        matches = r.get("matches")
        note = escape(str(r.get("note") or ""))
        if score is not None and not isfinite(float(score)):
            score = None

        if score is None:
            value_html = '<span class="nb-score-na">N/A</span>'
            width = 0.0
            row_class = "nb-score-row na"
        else:
            value_html = f"<span>{float(score):.4f}</span>"
            width = max(0.0, min(100.0, abs(float(score)) * 100.0))
            row_class = "nb-score-row pass" if float(score) >= threshold else "nb-score-row fail"

        # Only metrics that report a peak-match count get a badge.
        badge_html = (f'<div class="nb-score-badge">{int(matches)} matched</div>'
                      if matches is not None else "")
        tag_html = '<span class="nb-score-tag">Primary</span>' if r.get("primary") else ""
        fill_color = _similarity_score_color(score)
        note_html = f'<div class="nb-score-note">{note}</div>' if note else ""
        tick_html = ("" if score is None
                     else f'<div class="nb-score-tick" style="left:{tick:.1f}%;"></div>')
        row_html.append(
            f'<div class="{row_class}">'
            f'<div class="nb-score-top">'
            f'<div class="nb-score-name">{method}{tag_html}</div>'
            f'{badge_html}'
            f'<div class="nb-score-val">{value_html}</div>'
            '</div>'
            f'<div class="nb-score-track"><div class="nb-score-fill" '
            f'style="width:{width:.1f}%;background:{fill_color};"></div>{tick_html}</div>'
            f'{note_html}'
            '</div>'
        )

    n_rows = len(rows)
    return (
        '<div class="nb-scorecard">'
        '<div class="nb-scorecard-hd">'
        f'<div class="nb-scorecard-title">{escape(str(category))}</div>'
        f'<div class="nb-scorecard-count">{n_rows} method{"s" if n_rows != 1 else ""}</div>'
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
        subsection_label_html(label)
        + f'<div class="mono" style="{style}">{escape(str(value))}</div>'
    )


def copy_fields_html(
    fields: Sequence[tuple],
    heading: Optional[str] = None,
    variant: Optional[str] = None,
) -> str:
    """Build an identity card whose copyable rows have small copy buttons.

    Each field is ``(label, value, copyable)`` or ``(label, value, copyable, aux)``
    where ``aux`` is a short muted note shown after the value. ``variant`` adds a
    modifier class (``"adducts"``: mono, non-uppercase labels). Displayed values
    are HTML-escaped and clipboard payloads are base64-encoded, so user-supplied
    SMILES cannot inject markup or script into the rendered block.

    The copy icon is drawn in CSS: st.html sanitises with DOMPurify's HTML
    profile, which strips inline <svg>.
    """
    import hashlib

    rows: List[str] = []
    for idx, field in enumerate(fields):
        label, value, copyable = field[:3]
        aux = field[3] if len(field) > 3 else None
        lab = escape(str(label))
        val = escape(str(value))
        aux_html = f'<div class="aux">{escape(str(aux))}</div>' if aux else ""
        if copyable:
            payload = base64.b64encode(str(value).encode("utf-8")).decode("ascii")
            control = (
                f'<button class="copy-btn" type="button" data-copy="{payload}" '
                f'title="Copy {lab}" aria-label="Copy {lab}"></button>'
            )
        else:
            control = '<span class="copy-spacer"></span>'
        row_class = "row first" if idx == 0 else "row"
        rows.append(
            f'<div class="{row_class}"><div class="lab">{lab}</div>'
            f'<div class="val">{val}</div>{aux_html}{control}</div>'
        )

    heading_html = subsection_label_html(heading) if heading else ""
    body = heading_html + "".join(rows)
    # Content-derived id: the script binds only this block's buttons, and an
    # unchanged block keeps its DOM (and listeners) across Streamlit reruns.
    block_id = "nbid-" + hashlib.sha1(body.encode("utf-8")).hexdigest()[:12]
    cls = "nb-identity" + (f" {escape(variant)}" if variant else "")
    # Keep "<" out of the script body: st.html's DOMPurify pass silently drops
    # any <script> whose text contains "<" followed by a letter, "/" or "!".
    script = (
        "<script>(function(){"
        f"var root=document.getElementById('{block_id}');if(!root){{return;}}"
        "function b64ToText(s){try{"
        "var bin=atob(s);"
        "var bytes=Uint8Array.from(bin,function(c){return c.charCodeAt(0);});"
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
        "root.querySelectorAll('.copy-btn').forEach(function(btn){"
        "if(btn.dataset.bound){return;}btn.dataset.bound='1';"
        "btn.addEventListener('click',function(){"
        "copy(b64ToText(btn.getAttribute('data-copy')),btn);});"
        "});})();</script>"
    )
    return f'<div class="{cls}" id="{block_id}">{body}</div>{script}'


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
    # st.markdown strips <script>, so the motion hook goes in through st.html.
    st.html(MOTION_JS, unsafe_allow_javascript=True)


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


def section(label: str, aux: Optional[str] = None) -> None:
    """Level 1 section label — indigo square marker. Use for page sections."""
    st.markdown(section_label_html(label, aux), unsafe_allow_html=True)


def subsection(label: str, aux: Optional[str] = None) -> None:
    """Level 2 label — no marker, hairline rule. Use for blocks inside a section."""
    st.markdown(subsection_label_html(label, aux), unsafe_allow_html=True)


def panel(key: str):
    """Bordered brutalist panel for grouped inputs: ``with theme.panel("setup"): ...``.

    A thin wrapper around ``st.container(border=True)``; the key adds the
    ``st-key-nb-panel-<key>`` class that GLOBAL_CSS styles (ink frame + hard shadow).
    ``key`` must be unique on the page.
    """
    return st.container(border=True, key=f"nb-panel-{key}")


def metric_cards(items: Sequence[Tuple[str, str]]) -> None:
    st.markdown(metric_cards_html(items), unsafe_allow_html=True)


def stat_grid(rows: Sequence[Tuple[str, str]], wide: bool = False,
              max_cols: Optional[int] = None, min_px: int = 160) -> None:
    st.markdown(stat_grid_html(rows, wide=wide, max_cols=max_cols, min_px=min_px),
                unsafe_allow_html=True)


def verdict(category: str, subtitle: str, accent: str) -> None:
    st.markdown(verdict_html(category, subtitle, accent), unsafe_allow_html=True)


def chips(items: Sequence[Tuple[str, str]]) -> None:
    st.markdown(chips_html(items), unsafe_allow_html=True)


def chip_line(label: str, items: Sequence[Tuple[str, str]]) -> None:
    st.markdown(chip_line_html(label, items), unsafe_allow_html=True)


def reasons(items: Sequence[str], columns: int = 2) -> None:
    """Bulleted conclusions, flowed into balanced columns. Last entry = the verdict."""
    st.markdown(reasons_html(items, columns=columns), unsafe_allow_html=True)


def notes(items: Sequence[str], columns: int = 2) -> None:
    st.markdown(notes_html(items, columns=columns), unsafe_allow_html=True)


def similarity_score_cards(
    groups: Sequence[Tuple[str, Sequence[dict]]],
    threshold: float = 0.7,
) -> None:
    st.markdown(similarity_score_cards_html(groups, threshold=threshold), unsafe_allow_html=True)


def property_card(label: str, value: str, value_style: str = "") -> None:
    st.markdown(property_card_html(label, value, value_style), unsafe_allow_html=True)


def copy_fields(
    fields: Sequence[tuple],
    heading: Optional[str] = None,
    variant: Optional[str] = None,
) -> None:
    st.html(copy_fields_html(fields, heading, variant), unsafe_allow_javascript=True)


_3DMOL_JS_PATH = Path(__file__).resolve().parent / "lib" / "3Dmol-min.js"


@st.cache_data(show_spinner=False, max_entries=1)
def _load_3dmol_js() -> str:
    try:
        return _3DMOL_JS_PATH.read_text(encoding="utf-8")
    except OSError:
        return ""


def molecule_3d_html(
    molblock: str,
    height: int = 420,
    title: Optional[str] = None,
    caption: Optional[str] = None,
) -> str:
    """Self-contained 3Dmol.js viewer card for an MDL molblock.

    Optional header (``title`` left, muted ``caption`` right) and a style
    toolbar. The card leaves 5 px at the right/bottom of the iframe so its hard
    shadow is not clipped; atoms without bonds (ions, lone atoms) are drawn as
    spheres in every style; Reset restores the initial view. The output is
    deterministic, so an unchanged viewer keeps its state across reruns.
    """
    payload = base64.b64encode(molblock.encode("utf-8")).decode("ascii")
    lib_js = _load_3dmol_js()
    lib_tag = (
        f"<script>{lib_js}</script>"
        if lib_js
        else '<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>'
    )
    head_html = ""
    if title or caption:
        head_html = (
            '<div class="head">'
            f'<span class="ttl">{escape(title or "")}</span><span class="spacer"></span>'
            f'<span class="cap">{escape(caption or "")}</span></div>'
        )

    return f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
{lib_tag}
<style>
*, *::before, *::after {{ box-sizing:border-box; }}
html, body {{ margin:0; padding:0; height:100%; background:transparent; overflow:hidden;
  font-family:"Source Sans","Source Sans Pro",-apple-system,BlinkMacSystemFont,"Segoe UI",
    "Helvetica Neue",Arial,sans-serif; color:{INK}; }}
/* leave 5px inside the iframe for the hard shadow */
.wrap {{ display:flex; flex-direction:column; height:calc(100% - 5px); margin-right:5px;
  border:3px solid {INK}; border-radius:2px; box-shadow:5px 5px 0 {INK};
  background:{SURFACE}; overflow:hidden; animation:vfade .24s ease both; }}
/* matches the crossfade the 2D depiction gets, so switching views reads as one swap */
@keyframes vfade {{ from {{ opacity:0; }} to {{ opacity:1; }} }}
@media (prefers-reduced-motion: reduce) {{ .wrap {{ animation-duration:.001ms; }} }}
.head {{ display:flex; align-items:center; gap:8px; padding:13px 16px 4px; }}
.ttl {{ display:flex; align-items:center; gap:.55rem; font-size:11.5px; font-weight:800;
  text-transform:uppercase; letter-spacing:.14em; color:{MUTED}; }}
.ttl::before {{ content:""; width:11px; height:11px; background:{BLUE}; border:2px solid {INK}; }}
.cap {{ font-size:10px; font-weight:700; text-transform:uppercase; letter-spacing:.08em;
  color:{MUTED}; }}
.stage {{ flex:1; min-height:0; position:relative; background:{SURFACE}; }}
.stage canvas {{ display:block; }}
.msg {{ padding:24px 16px; font-size:13px; color:{MUTED}; }}
.bar {{ display:flex; flex-wrap:wrap; gap:6px; align-items:center;
  padding:8px 10px; border-top:3px solid {INK}; background:{BLUE_GHOST}; }}
.vbtn {{ font-family:ui-monospace,"SF Mono",Menlo,monospace; font-size:11px;
  font-weight:800; background:{SURFACE}; color:{INK}; border:2px solid {INK};
  border-radius:2px; padding:4px 9px; cursor:pointer; box-shadow:2px 2px 0 {INK};
  transition:background-color .12s ease; }}
.vbtn:hover {{ background:{BLUE_SOFT}; }}
.vbtn:active {{ transform:translate(2px,2px); box-shadow:none; transition:none; }}
.vbtn.active {{ background:{INK}; color:#FFFFFF; }}
.spacer {{ flex:1; }}
/* style buttons left, view actions right; when narrow the actions wrap as a group */
.grp {{ display:flex; flex-wrap:wrap; gap:6px; align-items:center; }}
.grp.end {{ margin-left:auto; }}
.hint {{ font-size:10px; font-weight:700; color:{MUTED}; text-transform:uppercase;
  letter-spacing:.08em; }}
@media (max-width:540px) {{ .hint {{ display:none; }} }}
@media (max-width:400px) {{ .cap {{ display:none; }} }}
</style>
</head>
<body>
<div class="wrap">
  {head_html}
  <div id="viewer" class="stage"></div>
  <div class="bar">
    <div class="grp">
      <button class="vbtn active" data-style="stick" type="button">Stick</button>
      <button class="vbtn" data-style="ballstick" type="button">Ball &amp; Stick</button>
      <button class="vbtn" data-style="sphere" type="button">Spacefill</button>
      <button class="vbtn" data-style="line" type="button">Line</button>
    </div>
    <div class="grp end">
      <span class="hint">drag · scroll</span>
      <button class="vbtn" id="spin" type="button">Spin</button>
      <button class="vbtn" id="reset" type="button">Reset</button>
    </div>
  </div>
</div>
<script>
(function(){{
  var molblock = atob("{payload}");
  if (typeof $3Dmol === "undefined") {{
    document.getElementById("viewer").innerHTML =
      '<div class="msg">3D viewer failed to load.</div>';
    return;
  }}
  var viewer = $3Dmol.createViewer(document.getElementById("viewer"), {{
    backgroundColor: "{SURFACE}"
  }});
  viewer.addModel(molblock, "sdf");
  // atoms without bonds (ions, lone atoms) are invisible as sticks/lines
  var LONE = {{ sphere: {{ scale: 0.35 }} }};
  function applyStyle(name) {{
    var style;
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
    if (name !== "sphere") {{
      viewer.setStyle({{ bonds: 0 }}, LONE);
    }}
    viewer.render();
  }}
  applyStyle("stick");
  viewer.zoomTo();
  viewer.render();
  var home = viewer.getView();
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
    viewer.setView(home);
    viewer.render();
  }});
}})();
</script>
</body>
</html>"""


def molecule_3d(
    molblock: str,
    height: int = 420,
    title: Optional[str] = None,
    caption: Optional[str] = None,
) -> None:
    components.html(
        molecule_3d_html(molblock, height, title=title, caption=caption),
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


# --------------------------------------------------------------------------
# Structure Annotation: result banner for a chemical formula (appended by the
# Structure tab). Same look as verdict(), but keeps element case - the uppercase
# verdict style would print Cl / Na / Br as "CL" / "NA" / "BR".
# --------------------------------------------------------------------------
def formula_verdict_html(formula: str, subtitle: str, accent: str) -> str:
    return (
        f'<div class="nb-verdict" style="background:{escape(str(accent))}">'
        f'<div class="cat" style="text-transform:none">{escape(str(formula))}</div>'
        f'<div class="sub">{escape(str(subtitle))}</div></div>'
    )


def formula_verdict(formula: str, subtitle: str, accent: str) -> None:
    st.markdown(formula_verdict_html(formula, subtitle, accent), unsafe_allow_html=True)
