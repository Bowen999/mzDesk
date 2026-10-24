from __future__ import annotations

import base64
import hashlib
from html import escape
from math import isfinite
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import streamlit as st
import streamlit.components.v1 as components

INK = "#1F2937"
BG = "#FFFFFF"
SURFACE = "#FFFFFF"
SURFACE_2 = "#F6F7F9"
MUTED = "#667085"
LINE = "#E4E7EC"
BLUE = "#2563EB"
BLUE_BRIGHT = "#1D4ED8"
BLUE_SOFT = "#BFDBFE"
BLUE_GHOST = "#EFF6FF"

CATEGORY_ACCENT = {
    "Lipidome": "#2563EB",
    "Metabolome": "#0F766E",
    "Both": "#7C3AED",
    "Neither": "#64748B",
    "Invalid": "#B42318",
}

PLOT_DARK = BLUE
PLOT_LIGHT = "#94A3B8"
PLOT_FONT = '"Helvetica Neue",Helvetica,Arial,sans-serif'
PLOT_MONO = 'ui-monospace,"SF Mono",Menlo,Consolas,monospace'

FRAG_HL = {
    "acyl": (0.106, 0.227, 0.561),
    "isoprenoid": (0.055, 0.486, 0.608),
    "head": (0.333, 0.400, 0.800),
    "glycosyl": (0.902, 0.624, 0.208),
}

GLOBAL_CSS = f"""
<style>
:root{{
  --ink:{INK}; --bg:{BG}; --surface:{SURFACE}; --surface-2:{SURFACE_2};
  --muted:{MUTED}; --line:{LINE}; --blue:{BLUE}; --blue-bright:{BLUE_BRIGHT};
  --blue-soft:{BLUE_SOFT}; --blue-ghost:{BLUE_GHOST};
}}
html, body, [data-testid="stAppViewContainer"], [data-testid="stApp"]{{background:var(--bg);}}
html, body{{
  font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","Helvetica Neue",Arial,sans-serif;
  color:var(--ink);
}}
[data-testid="stHeader"]{{background:transparent;pointer-events:none;}}
#MainMenu, footer, [data-testid="stToolbar"], [data-testid="stDecoration"]{{visibility:hidden;}}
[data-testid="stStatusWidget"]{{visibility:visible;pointer-events:auto;}}
.block-container{{max-width:1200px;padding-top:1.4rem;padding-bottom:2.5rem;}}
h1{{font-size:1.8rem!important;font-weight:650!important;letter-spacing:-.02em;}}
a{{color:var(--blue);}}
.mono{{font-family:ui-monospace,"SF Mono",Menlo,Consolas,monospace;}}
.mz-section{{
  display:flex;align-items:baseline;gap:.75rem;margin:.9rem 0 .55rem;
  font-size:1rem;font-weight:650;color:var(--ink);
}}
.mz-subsection{{
  display:flex;align-items:baseline;gap:.75rem;margin:.6rem 0 .45rem;
  font-size:.84rem;font-weight:600;color:var(--ink);
}}
.mz-section .aux, .mz-subsection .aux{{
  margin-left:auto;font-size:.72rem;font-weight:400;color:var(--muted);white-space:nowrap;
}}
.mz-note{{font-size:.82rem;line-height:1.45;color:var(--muted);margin:.35rem 0 .6rem;}}
.mz-reasons{{margin:.2rem 0 .6rem;padding-left:1.2rem;}}
.mz-reasons li{{font-size:.88rem;line-height:1.5;margin:.25rem 0;}}
.mz-reasons li:last-child{{font-weight:600;}}
.mz-verdict{{
  border:1px solid var(--line);border-left:4px solid var(--blue);border-radius:6px;
  background:var(--surface);padding:14px 16px;margin:.7rem 0 1rem;
}}
.mz-verdict .cat{{font-size:1.55rem;font-weight:650;line-height:1.2;}}
.mz-verdict .sub{{font-size:.86rem;color:var(--muted);margin-top:4px;}}
.mz-statgrid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:8px;margin:.3rem 0 .8rem;}}
.mz-statgrid.wide{{grid-template-columns:repeat(auto-fit,minmax(150px,1fr));}}
.mz-stat{{min-width:0;border:1px solid var(--line);border-radius:5px;background:var(--surface);padding:9px 11px;}}
.mz-stat .lab{{font-size:.7rem;line-height:1.3;color:var(--muted);}}
.mz-stat .val{{
  font-family:ui-monospace,"SF Mono",Menlo,monospace;font-size:1rem;font-weight:600;
  line-height:1.35;color:var(--ink);margin-top:3px;overflow-wrap:anywhere;
}}
.mz-chips{{display:flex;align-items:center;flex-wrap:wrap;gap:6px;margin:.25rem 0 .5rem;}}
.mz-chip{{
  display:inline-flex;align-items:center;gap:5px;border:1px solid var(--line);
  border-radius:4px;background:var(--surface);padding:3px 7px;font-size:.72rem;color:var(--muted);
}}
.mz-chip .sw{{width:9px;height:9px;border-radius:2px;display:inline-block;}}
.mz-chipline{{display:flex;align-items:center;flex-wrap:wrap;gap:.6rem;}}
.mz-chipline .lb{{font-size:.78rem;color:var(--muted);}}
.mz-chipline .mz-chips{{margin:0;}}
.mz-fields{{
  container-type:inline-size;border:1px solid var(--line);border-radius:6px;
  background:var(--surface);padding:10px 12px;margin-bottom:10px;
}}
.mz-fields .mz-subsection{{margin-top:0;}}
.mz-fields .row{{
  display:flex;align-items:flex-start;gap:10px;border-top:1px solid var(--line);padding:7px 0;
}}
.mz-fields .row.first{{border-top:0;}}
.mz-fields .row .lab{{flex:0 0 112px;font-size:.72rem;color:var(--muted);padding-top:2px;}}
.mz-fields .row .val{{
  flex:1;min-width:0;font-family:ui-monospace,"SF Mono",Menlo,monospace;
  font-size:.75rem;line-height:1.5;overflow-wrap:anywhere;
}}
.mz-fields .row .aux{{flex:none;font-size:.68rem;color:var(--muted);white-space:nowrap;}}
.mz-fields.adducts .row{{display:grid;grid-template-columns:minmax(0,auto) 1fr auto 20px;align-items:center;}}
.mz-fields.adducts .row .lab{{padding-top:0;font-family:ui-monospace,"SF Mono",Menlo,monospace;color:var(--ink);}}
.mz-fields.adducts .row .val{{text-align:right;white-space:nowrap;}}
.mz-copy{{
  position:relative;flex:none;width:18px;height:18px;padding:0;border:0;border-radius:3px;
  background:transparent;color:var(--muted);cursor:pointer;
}}
.mz-copy:hover{{background:var(--blue-ghost);color:var(--blue);}}
.mz-copy::before, .mz-copy::after{{
  content:"";position:absolute;box-sizing:border-box;width:8px;height:9px;
  border:1.5px solid currentColor;border-radius:1px;
}}
.mz-copy::before{{left:2px;top:2px;}}
.mz-copy::after{{left:6px;top:6px;background:var(--surface);}}
.mz-copy.copied::before{{display:none;}}
.mz-copy.copied::after{{
  left:6px;top:2px;width:5px;height:10px;background:none;
  border-width:0 2px 2px 0;transform:rotate(45deg);
}}
.mz-copy-spacer{{flex:none;width:18px;}}
@container (max-width:340px){{
  .mz-fields:not(.adducts) .row{{flex-wrap:wrap;row-gap:2px;}}
  .mz-fields:not(.adducts) .row .lab{{flex-basis:100%;}}
  .mz-fields.adducts .row .aux{{display:none;}}
}}
[data-testid="stBaseButton-primary"], [data-testid="stBaseButton-secondary"]{{
  border:1px solid var(--line);border-radius:5px;box-shadow:none;padding:.45rem 1rem;
}}
[data-testid="stBaseButton-primary"]{{background:var(--blue);border-color:var(--blue);color:#fff;}}
[data-testid="stBaseButton-primary"]:hover{{background:var(--blue-bright);border-color:var(--blue-bright);color:#fff;}}
[data-testid="stBaseButton-secondary"]{{background:var(--surface);color:var(--ink);}}
[data-testid="stBaseButton-secondary"]:hover{{background:var(--surface-2);color:var(--ink);border-color:var(--line);}}
[data-testid="stBaseButton-primary"]:active, [data-testid="stBaseButton-secondary"]:active{{transform:none;box-shadow:none;}}
[data-testid="stBaseButton-primary"]:focus-visible, [data-testid="stBaseButton-secondary"]:focus-visible{{
  outline:2px solid var(--blue);outline-offset:2px;
}}
[data-testid="stTextInputRootElement"], [data-testid="stTextAreaRootElement"],
[data-testid="stNumberInputContainer"]{{
  border:1px solid var(--line)!important;border-radius:5px!important;box-shadow:none;background:var(--surface)!important;
}}
[data-testid="stTextInputRootElement"]:focus-within, [data-testid="stTextAreaRootElement"]:focus-within,
[data-testid="stNumberInputContainer"]:focus-within{{border-color:var(--blue)!important;}}
[data-testid="stTextInputRootElement"] [data-baseweb="base-input"],
[data-testid="stTextAreaRootElement"] [data-baseweb="base-input"],
[data-testid="stNumberInputContainer"] [data-baseweb="input"],
[data-testid="stNumberInputContainer"] [data-baseweb="base-input"]{{border:0!important;background:var(--surface)!important;}}
[data-testid="stTextInput"] input, [data-testid="stNumberInput"] input,
[data-testid="stTextArea"] textarea{{background:var(--surface)!important;color:var(--ink)!important;}}
[data-testid="stTextArea"] textarea{{line-height:20px;padding-top:9px;padding-bottom:9px;}}
div[data-baseweb="select"]>div{{
  border:1px solid var(--line)!important;border-radius:5px!important;box-shadow:none;background:var(--surface)!important;
}}
div[data-baseweb="select"]>div:focus-within{{border-color:var(--blue)!important;}}
label, [data-testid="stWidgetLabel"] p{{font-weight:500!important;color:var(--ink)!important;font-size:.82rem!important;}}
[data-testid="stNumberInputStepUp"], [data-testid="stNumberInputStepDown"]{{border-left:1px solid var(--line)!important;}}
[data-testid="stCheckbox"] label[data-baseweb="checkbox"]>span:first-child{{border:1px solid var(--line)!important;}}
[data-testid="stBaseButton-segmented_control"], [data-testid="stBaseButton-segmented_controlActive"]{{
  border:1px solid var(--line)!important;border-radius:4px!important;box-shadow:none;
}}
[data-testid="stBaseButton-segmented_control"]{{background:var(--surface)!important;color:var(--ink)!important;}}
[data-testid="stBaseButton-segmented_controlActive"]{{background:var(--blue)!important;color:#fff!important;}}
[data-testid="stSlider"] [data-baseweb="slider"] div[role="slider"]{{border:1px solid var(--blue);background:var(--blue);box-shadow:none;}}
section[data-testid="stSidebar"]{{background:var(--surface-2);border-right:1px solid var(--line);}}
[data-testid="stImage"] img{{border:1px solid var(--line);border-radius:5px;background:#fff;}}
[data-baseweb="tab-list"]{{gap:18px;border-bottom:1px solid var(--line);margin-bottom:.8rem;}}
[data-baseweb="tab-list"] button[data-baseweb="tab"]{{
  background:transparent;border:0;border-bottom:2px solid transparent;border-radius:0;
  padding:8px 0;margin-bottom:-1px;color:var(--muted);
}}
[data-baseweb="tab-list"] button[data-baseweb="tab"]:hover{{background:transparent;color:var(--ink);}}
[data-baseweb="tab-list"] button[aria-selected="true"]{{background:transparent;color:var(--blue)!important;border-bottom-color:var(--blue);box-shadow:none;}}
[data-baseweb="tab-list"] button[aria-selected="true"] p{{color:var(--blue)!important;}}
[data-baseweb="tab-list"] button[data-baseweb="tab"] p{{font-size:.84rem;font-weight:500;}}
[data-baseweb="tab-highlight"], [data-baseweb="tab-border"]{{background:transparent!important;height:0!important;}}
[data-testid="stExpander"] details{{border:1px solid var(--line)!important;border-radius:5px!important;box-shadow:none;background:var(--surface);}}
[data-testid="stExpander"] summary, [data-testid="stExpander"] summary p{{font-weight:500;font-size:.82rem;}}
[data-testid="stDataFrame"]{{border:1px solid var(--line);border-radius:5px;box-shadow:none;}}
[class*="st-key-mz-panel-"]{{border:1px solid var(--line)!important;border-radius:6px!important;box-shadow:none!important;background:var(--surface)!important;}}
[data-testid="stFileUploaderDropzone"]{{border:1px dashed var(--line)!important;border-radius:5px!important;background:var(--surface-2)!important;}}
[data-testid="stAlert"]{{border:1px solid var(--line);border-radius:5px;box-shadow:none;}}
[data-testid="stAlertContainer"]{{border-radius:4px;}}
@media (max-width:900px){{
  [data-baseweb="tab-list"]{{gap:12px;}}
  [data-baseweb="tab-list"] button[data-baseweb="tab"] p{{font-size:.76rem;}}
}}
@media (max-width:640px){{
  .mz-statgrid{{grid-template-columns:1fr;}}
  [data-baseweb="tab-list"]{{gap:9px;}}
  [data-baseweb="tab-list"] button[data-baseweb="tab"] p{{font-size:.7rem;}}
}}
</style>
"""


def inject_css() -> None:
    st.markdown(GLOBAL_CSS, unsafe_allow_html=True)


def _label_html(label: str, aux: Optional[str], class_name: str) -> str:
    aux_html = f'<span class="aux">{escape(str(aux))}</span>' if aux else ""
    return f'<div class="{class_name}">{escape(str(label))}{aux_html}</div>'


def section(label: str, aux: Optional[str] = None) -> None:
    st.markdown(_label_html(label, aux, "mz-section"), unsafe_allow_html=True)


def subsection(label: str, aux: Optional[str] = None) -> None:
    st.markdown(_label_html(label, aux, "mz-subsection"), unsafe_allow_html=True)


def panel(key: str):
    return st.container(border=True, key=f"mz-panel-{key}")


def stat_grid(rows: Sequence[Tuple[str, str]], wide: bool = False,
              max_cols: Optional[int] = None, min_px: int = 150) -> None:
    cells = "".join(
        f'<div class="mz-stat"><div class="lab">{escape(str(label))}</div>'
        f'<div class="val">{escape(str(value))}</div></div>'
        for label, value in rows
    )
    class_name = "mz-statgrid wide" if wide else "mz-statgrid"
    style = ""
    if wide:
        floor = max(96, int(min_px))
        track = (f"max({floor}px,calc({100 / max_cols:g}% - 10px))" if max_cols
                 else f"{floor}px")
        style = f' style="grid-template-columns:repeat(auto-fit,minmax({track},1fr));"'
    st.markdown(f'<div class="{class_name}"{style}>{cells}</div>', unsafe_allow_html=True)


def metric_cards(items: Sequence[Tuple[str, str]]) -> None:
    stat_grid(items, wide=True)


def _verdict(value: str, subtitle: str, accent: str, preserve_case: bool = False) -> None:
    style = "text-transform:none" if preserve_case else ""
    st.markdown(
        f'<div class="mz-verdict" style="border-left-color:{escape(str(accent))}">'
        f'<div class="cat" style="{style}">{escape(str(value))}</div>'
        f'<div class="sub">{escape(str(subtitle))}</div></div>',
        unsafe_allow_html=True,
    )


def verdict(category: str, subtitle: str, accent: str) -> None:
    _verdict(category, subtitle, accent)


def formula_verdict(formula: str, subtitle: str, accent: str) -> None:
    _verdict(formula, subtitle, accent, preserve_case=True)


def chips(items: Sequence[Tuple[str, str]]) -> None:
    body = "".join(
        f'<span class="mz-chip"><span class="sw" style="background:{escape(str(color))};"></span>'
        f'{escape(str(label))}</span>'
        for label, color in items
    )
    st.markdown(f'<div class="mz-chips">{body}</div>', unsafe_allow_html=True)


def chip_line(label: str, items: Sequence[Tuple[str, str]]) -> None:
    body = "".join(
        f'<span class="mz-chip"><span class="sw" style="background:{escape(str(color))};"></span>'
        f'{escape(str(item))}</span>'
        for item, color in items
    )
    st.markdown(
        f'<div class="mz-chipline"><span class="lb">{escape(str(label))}</span>'
        f'<div class="mz-chips">{body}</div></div>',
        unsafe_allow_html=True,
    )


def reasons(items: Sequence[str]) -> None:
    body = "".join(f'<li>{escape(str(item))}</li>' for item in items)
    st.markdown(f'<ul class="mz-reasons">{body}</ul>', unsafe_allow_html=True)


def note(text: str) -> None:
    st.markdown(f'<div class="mz-note">{escape(str(text))}</div>', unsafe_allow_html=True)


def divider() -> None:
    st.divider()


def similarity_score_cards(groups: Sequence[Tuple[str, Sequence[dict]]],
                           threshold: float = 0.7) -> None:
    rows = []
    for category, results in groups:
        for result in results:
            score = result.get("score")
            if score is not None and not isfinite(float(score)):
                score = None
            method = str(result.get("method", ""))
            if result.get("primary"):
                method += " (primary)"
            matches = result.get("matches")
            rows.append({
                "Category": category,
                "Method": method,
                "Score": "N/A" if score is None else f"{float(score):.4f}",
                "Pass": "—" if score is None else ("Yes" if float(score) >= threshold else "No"),
                "Matched": "—" if matches is None else str(int(matches)),
                "Note": result.get("note") or "",
            })
    st.dataframe(rows, width="stretch", hide_index=True)


def copy_fields(fields: Sequence[tuple], heading: Optional[str] = None,
                variant: Optional[str] = None) -> None:
    rows: List[str] = []
    for index, field in enumerate(fields):
        label, value, copyable = field[:3]
        aux = field[3] if len(field) > 3 else None
        lab = escape(str(label))
        val = escape(str(value))
        aux_html = f'<div class="aux">{escape(str(aux))}</div>' if aux else ""
        if copyable:
            payload = base64.b64encode(str(value).encode("utf-8")).decode("ascii")
            control = (
                f'<button class="mz-copy" type="button" data-copy="{payload}" '
                f'title="Copy {lab}" aria-label="Copy {lab}"></button>'
            )
        else:
            control = '<span class="mz-copy-spacer"></span>'
        row_class = "row first" if index == 0 else "row"
        rows.append(
            f'<div class="{row_class}"><div class="lab">{lab}</div>'
            f'<div class="val">{val}</div>{aux_html}{control}</div>'
        )

    heading_html = _label_html(heading, None, "mz-subsection") if heading else ""
    body = heading_html + "".join(rows)
    block_id = "mz-fields-" + hashlib.sha1(body.encode("utf-8")).hexdigest()[:12]
    class_name = "mz-fields" + (f" {escape(variant)}" if variant else "")
    script = (
        "<script>(function(){"
        f"var root=document.getElementById('{block_id}');if(!root){{return;}}"
        "function decode(s){try{var bin=atob(s);"
        "var bytes=Uint8Array.from(bin,function(c){return c.charCodeAt(0);});"
        "if(typeof TextDecoder==='undefined'){return bin;}"
        "return new TextDecoder().decode(bytes);}catch(e){return s;}}"
        "function fallback(text){var ta=document.createElement('textarea');"
        "ta.value=text;ta.style.position='fixed';ta.style.opacity='0';"
        "document.body.appendChild(ta);ta.select();var ok=false;"
        "try{ok=document.execCommand('copy');}catch(e){}"
        "document.body.removeChild(ta);return ok;}"
        "function flash(btn){btn.classList.add('copied');"
        "window.setTimeout(function(){btn.classList.remove('copied');},1200);}"
        "function copy(text,btn){if(navigator.clipboard&&window.isSecureContext){"
        "navigator.clipboard.writeText(text).then(function(){flash(btn);},"
        "function(){if(fallback(text)){flash(btn);}});}else if(fallback(text)){flash(btn);}}"
        "root.querySelectorAll('.mz-copy').forEach(function(btn){"
        "if(btn.dataset.bound){return;}btn.dataset.bound='1';"
        "btn.addEventListener('click',function(){copy(decode(btn.getAttribute('data-copy')),btn);});"
        "});})();</script>"
    )
    st.html(f'<div class="{class_name}" id="{block_id}">{body}</div>{script}',
            unsafe_allow_javascript=True)


_3DMOL_JS_PATH = Path(__file__).resolve().parents[1] / "lib" / "3Dmol-min.js"


@st.cache_data(show_spinner=False, max_entries=1)
def _load_3dmol_js() -> str:
    try:
        return _3DMOL_JS_PATH.read_text(encoding="utf-8")
    except OSError:
        return ""


def molecule_3d_html(molblock: str, height: int = 420, title: Optional[str] = None,
                     caption: Optional[str] = None) -> str:
    payload = base64.b64encode(molblock.encode("utf-8")).decode("ascii")
    lib_js = _load_3dmol_js()
    lib_tag = (f"<script>{lib_js}</script>" if lib_js else
               '<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>')
    head_html = ""
    if title or caption:
        head_html = (
            '<div class="head">'
            f'<span>{escape(title or "")}</span><span class="spacer"></span>'
            f'<span class="cap">{escape(caption or "")}</span></div>'
        )

    return f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
{lib_tag}
<style>
*, *::before, *::after {{box-sizing:border-box;}}
html, body{{margin:0;padding:0;height:100%;background:transparent;overflow:hidden;
  font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Arial,sans-serif;color:{INK};}}
.wrap{{display:flex;flex-direction:column;height:100%;border:1px solid {LINE};
  border-radius:6px;background:{SURFACE};overflow:hidden;}}
.head{{display:flex;gap:8px;padding:10px 12px 2px;font-size:12px;}}
.cap{{font-size:11px;color:{MUTED};}}
.stage{{flex:1;min-height:0;position:relative;background:{SURFACE};}}
.stage canvas{{display:block;}}
.msg{{padding:20px 14px;font-size:13px;color:{MUTED};}}
.bar{{display:flex;flex-wrap:wrap;gap:6px;align-items:center;padding:7px 9px;
  border-top:1px solid {LINE};background:{SURFACE_2};}}
.vbtn{{font-size:11px;background:{SURFACE};color:{INK};border:1px solid {LINE};
  border-radius:4px;padding:4px 8px;cursor:pointer;}}
.vbtn:hover{{background:{BLUE_GHOST};}}
.vbtn.active{{background:{BLUE};border-color:{BLUE};color:#fff;}}
.spacer{{flex:1;}}
.grp{{display:flex;flex-wrap:wrap;gap:6px;align-items:center;}}
.grp.end{{margin-left:auto;}}
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
      <button class="vbtn" id="spin" type="button">Spin</button>
      <button class="vbtn" id="reset" type="button">Reset</button>
    </div>
  </div>
</div>
<script>
(function(){{
  var molblock = atob("{payload}");
  if (typeof $3Dmol === "undefined") {{
    document.getElementById("viewer").innerHTML = '<div class="msg">3D viewer failed to load.</div>';
    return;
  }}
  var viewer = $3Dmol.createViewer(document.getElementById("viewer"), {{backgroundColor:"{SURFACE}"}});
  viewer.addModel(molblock, "sdf");
  var lone = {{sphere:{{scale:0.35}}}};
  function applyStyle(name) {{
    var style;
    if (name === "ballstick") {{
      style = {{stick:{{radius:0.12}},sphere:{{scale:0.25}}}};
    }} else if (name === "sphere") {{
      style = {{sphere:{{scale:0.42}}}};
    }} else if (name === "line") {{
      style = {{line:{{linewidth:2}}}};
    }} else {{
      style = {{stick:{{radius:0.15}}}};
    }}
    viewer.setStyle({{}},style);
    if (name !== "sphere") {{viewer.setStyle({{bonds:0}},lone);}}
    viewer.render();
  }}
  applyStyle("stick");
  viewer.zoomTo();
  viewer.render();
  var home = viewer.getView();
  window.addEventListener("resize",function(){{viewer.resize();viewer.render();}});
  var spinning = false;
  document.querySelectorAll("[data-style]").forEach(function(btn){{
    btn.addEventListener("click",function(){{
      document.querySelectorAll("[data-style]").forEach(function(x){{x.classList.remove("active");}});
      btn.classList.add("active");
      applyStyle(btn.dataset.style);
    }});
  }});
  document.getElementById("spin").addEventListener("click",function(){{
    spinning=!spinning;
    this.classList.toggle("active",spinning);
    viewer.spin(spinning ? "y" : false,1);
  }});
  document.getElementById("reset").addEventListener("click",function(){{viewer.setView(home);viewer.render();}});
}})();
</script>
</body>
</html>"""


def molecule_3d(molblock: str, height: int = 420, title: Optional[str] = None,
                caption: Optional[str] = None) -> None:
    components.html(molecule_3d_html(molblock, height, title, caption),
                    height=height, scrolling=False)
