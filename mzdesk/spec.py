"""页面定义（唯一来源）。前端按这里的 pipelines / modules 渲染输入和输出格式。

字段类型：structure, precursor, spectrum, ms1, number, text, textarea, select, multiselect, radio,
checkbox, slider, file, table, heading。
输出块类型：table, cards, figure, code, formats, notes, button, heading, portal, chips, reaction, grouptable。
  任何块都可带 help（"?" 弹窗，字符串或字符串列表；运行时数据里的 data.help 会并进去）。数据形状见各 helper 的 docstring，
  运行时数据可用 mzdesk/api/blocks.py 里的 builder 生成。
  formats：同一份内容的多种格式，数据 {"formats": {格式名: 文本}}，前端画成一排标签 + 一个代码框。
  id：后端结果按 id 填进这个块；没有 id 的块只画格式。
  key：块的标签，前端 Outputs 可按 key 排除块；pipeline 的 stage 不再按 key 过滤，步骤的块由后端运行时按 id 给出。
模块：group 是侧栏大类（modules / pipelines）；subgroup 是 Modules 下的小类（molecular / detectability / ms1 / ms2 / orthogonal / gnps，以及 2.0 才开放的 dmpk / inorganic；ms2 的显示名是 Spectrum），
  None = 直接列在 Modules 下（Calculator）。DEFERRED_* 里的模块已做好但本版不显示。live=True 表示已接引擎（后端 mzdesk/api/runners/ 包里有对应 runner）。
  example：一组可直接运行的示例输入，前端 "Load example" 按钮填入。
  one_liner：一行 "<已有的输入> → <得到的结果>"，纯文本，不写引擎名（引擎名放 based_on）；没有哪个页面再标 "(Local version only)"：Docker / Mono / FRIGID 的页面在 mzdesk.com 的服务器上也装好了（deploy/install-engines.sh）。
  based_on：页脚 "Based on …" 的名单（module() 的参数；PIPELINES 的每个字典里也有）。页面真正依赖的外部引擎 / 模型 / 数据在前，纯 mzDesk 自己的规则写 "mzDesk"，最多 4 个。
  help / column_help / caption 只写读结果所需的一两句短话；删掉的长说明原样存在 mzdesk/docs/<id>.md（应用内 Docs 页显示）。
"""
import base64

from mzdesk.shared.ms2_clean import NORMALIZE, SORT_BY
# 注意：resolve.py 导入时不联网，所以这里可以直接导入 KINDS；不要在 spec 里导入任何联网代码
from mzdesk.identity.resolve import KINDS as IDENTITY_KINDS
# cfmid_engine 顶层只用标准库，导入它不会拉进 rdkit / Docker
from mzdesk.ms2.cfmid_engine import ADDUCTS as ESI_ADDUCTS
# msfinder_engine 顶层只用标准库，这里只读它的加合物名单
from mzdesk.ms2.msfinder_engine import ADDUCTS as MSFINDER_ADDUCT_NAMES
# frigid_engine 顶层也只用标准库（memguard / heavy_lock 都是纯标准库），这里只读加合物名单和 Quick / Full 的显示名
from mzdesk.ms2.frigid_engine import ADDUCTS as FRIGID_ADDUCT_NAMES, MODE_LABELS as FRIGID_MODES, workflow_default as frigid_workflow_default
# 指标名字只用标准库（和 ms2_similarity.py 共用一份），所以这里可以直接导入
from mzdesk.ms2 import metric_names as METRIC_NAMES
# isotopes 顶层只用标准库（IsoSpecPy 在用到时才导入）
from mzdesk.ms1.isotopes import RESOLUTIONS
# Inorganic MS 的选项是纯静态数据；同位素表在 runner 内才读取
from mzdesk.inorganic.data import ELEMENT_SYMBOLS, MATRIX_ELEMENTS, PARTICLE_PRESETS, SENSITIVITY_DWELL, SENSITIVITY_RATE

ADDUCTS_POS = ["Auto", "[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M]+", "[M+H-H2O]+", "[2M+H]+"]
ADDUCTS_NEG = ["Auto", "[M-H]-", "[M+Cl]-", "[M+HCOO]-", "[M+CH3COO]-", "[M-H-H2O]-", "[2M-H]-"]
# 引擎各自支持的加合物（Auto = [M+H]+ / [M-H]-）
METFRAG_ADDUCTS = {"Positive": ["Auto", "[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M]+"],
                   "Negative": ["Auto", "[M-H]-", "[M+Cl]-", "[M+HCOO]-", "[M+CH3COO]-"]}
BUDDY_ADDUCTS = {"Positive": ["Auto", "[M+H]+", "[M+Na]+", "[M+NH4]+", "[M+K]+", "[M+H-H2O]+"],
                 "Negative": ["Auto", "[M-H]-", "[M+Cl]-", "[M+HCOO]-", "[M+CH3COO]-", "[M+Br]-"]}
# MS-FINDER 用自己的加合物写法（FA = formate, Hac = acetate）；Auto = [M+H]+ / [M-H]-
MSFINDER_ADDUCTS = {mode: ["Auto"] + names for mode, names in MSFINDER_ADDUCT_NAMES.items()}
# FRIGID 只做正离子（负离子时 runner 返回说明），所以 Negative 下只有 Auto；Auto = [M+H]+
FRIGID_ADDUCTS = {"Positive": ["Auto"] + FRIGID_ADDUCT_NAMES, "Negative": ["Auto"]}
LIBRARIES = ["MassBank", "GNPS"]            # 本地索引见 mzdesk/ms2/library_setup.py
LIBRARY_METRICS = METRIC_NAMES.LIBRARY_ORDER    # Library search 页的指标多选 = Similarity 页得分表里的全部 16 项（同一份名字，谱熵排第一）
LIBRARY_DEFAULT_METRICS = [METRIC_NAMES.ENTROPY_WEIGHTED, METRIC_NAMES.COSINE_GREEDY]
# CCS: platform / gas / TWIMS calibration choices, shared by the CCS page and the ESI MS/MS Annotation workflow's Ion Mobility section
# (mzdesk/orthogonal/ccs.py PLATFORMS, GASES and TW_CALIBRATIONS; tests/test_annotation_ccs.py checks they agree)
CCS_PLATFORMS = ["DTIMS (stepped field)", "DTIMS (single field)", "TIMS", "TWIMS", "FAIMS / DMS"]
CCS_GASES = ["N₂", "He"]
CCS_TW_CALIBRATIONS = ["N₂ reference values", "He reference values (N₂→He)"]
DATABASES = ["HMDB", "PubChem", "LIPID MAPS", "KEGG", "ChEBI", "DrugBank", "NORMAN SusDat"]
# The extra structure-field mode of Lipid Fragment Explanation: a lipid shorthand, read by mzdesk.lipid_name_map (local, no lookup)
LIPID_NAME_MODE = "Lipid name"


# ---- 字段 ----
def structure(id="structure", label="Compound", optional=False):
    """network_modes: input modes the backend resolves online (SMILES / InChI and formulas are parsed locally).
    The frontend asks before using the internet when the switch is off; the wire value is components/Fields.tsx NAME_MODE.
    label_field=True adds an optional label input inside the field, whose text names the compound on the result cards
    (read from the field value's "label" key by mzdesk/api/runners/structure.py)."""
    return {"type": "structure", "id": id, "label": label, "optional": optional, "network_modes": ["Name / ID / Formula"]}
def precursor(id="precursor", adduct=True, label="Precursor m/z", adducts=None, mz=True):
    """mz=False: 不显示 m/z 输入框（只选离子模式和加合物）；该字段的加合物列表里没有 Auto 时，默认取列表第一项。"""
    f = {"type": "precursor", "id": id, "adduct": adduct, "label": label,
         "adducts": adducts or {"Positive": ADDUCTS_POS, "Negative": ADDUCTS_NEG}}
    if not mz:
        f["mz"] = False
    return f
def spectrum(id="ms2", label="MS2 spectrum", multi=False):
    """network_sources: input sources fetched online (see mzdesk/api/common.py spectrum())."""
    return {"type": "spectrum", "id": id, "label": label, "multi": multi, "network_sources": ["USI"]}
def ms1(id="ms1"):
    return {"type": "ms1", "id": id, "label": "MS1 isotope pattern (optional)", "placeholder": "m/z intensity per line"}
def number(id, label, value=None, step=None, placeholder=None):
    return {"type": "number", "id": id, "label": label, "value": value, "step": step, "placeholder": placeholder}
def text(id, label, value="", placeholder=None):
    return {"type": "text", "id": id, "label": label, "value": value, "placeholder": placeholder}
def textarea(id, label, value="", placeholder=None, rows=5):
    return {"type": "textarea", "id": id, "label": label, "value": value, "placeholder": placeholder, "rows": rows}
def select(id, label, options, value=None):
    return {"type": "select", "id": id, "label": label, "options": options, "value": value or options[0]}
def multiselect(id, label, options, value=None, network_options=None):
    """network_options: options that use the internet when selected (e.g. the deep-learning scores)."""
    f = {"type": "multiselect", "id": id, "label": label, "options": options, "value": value if value is not None else options}
    if network_options:
        f["network_options"] = list(network_options)
    return f
def radio(id, label, options, value=None):
    return {"type": "radio", "id": id, "label": label, "options": options, "value": value or options[0]}
def checkbox(id, label, value=False, network=False):
    """network: ticking it uses the internet."""
    f = {"type": "checkbox", "id": id, "label": label, "value": value}
    if network:
        f["network"] = True
    return f
def slider(id, label, lo, hi, step, value):
    return {"type": "slider", "id": id, "label": label, "min": lo, "max": hi, "step": step, "value": value}
def file(id, label, accept):
    return {"type": "file", "id": id, "label": label, "accept": accept}
def grid(id, label, columns, caption=None):
    return {"type": "table", "id": id, "label": label, "columns": columns, "caption": caption}
def heading(label):
    return {"type": "heading", "label": label}


# ---- 输出块 ----
def _extras(block, **extra):
    """Add the optional keys that are set; a block that uses none of them serialises exactly as before."""
    block.update({k: v for k, v in extra.items() if v})
    return block


# help (every block): string or list of strings, shown in a "?" popover after the block title (after the label for figure / code / formats).
#   A list renders as bullets. The page may add runtime text in data["help"] (string or list); static and runtime text are merged.
#   Text may use <sub> <sup> <i> <b>. Items inside cards / chips / portal tiles / table rows / groups take their own "help" in the data.
def table(columns, title=None, caption=None, key=None, id=None, help=None, row_help=False, page_size=None, compact=False, column_help=None):
    """Table. data: {"rows": [row, ...], "count"?: int, "columns"?: [str, ...]}   count -> small pill next to the title.
    data["columns"] replaces the spec columns for that run (the spec columns are what the page shows before a run); each row then has one cell per data column.
    cell: "text" | {"text", "url"} link (new tab) | {"badge": text, "tone": jp|report|good|warn|bad|neutral} | {"text", "help"} inline "?"
          | {"text", "hover_svg": svg} the 2D structure is drawn on hover over the text | {"img": data-URL}; badge / text cells may also carry
          "url" and "help".   (the old {"text", "href"} link still works)
    row:  [cell, ...] | {"cells": [cell, ...], "help": str|list}.  row_help=True adds a trailing "?" column that shows each row's help.
    page_size paginates and filters returned preview rows; compact lets long text expand without changing the data.
    column_help: {column name: short text} puts a small "?" after that column's name."""
    return _extras({"type": "table", "title": title, "columns": columns, "caption": caption, "key": key, "id": id},
                   help=help, row_help=row_help, page_size=page_size, compact=compact, column_help=column_help)
def cards(items, title=None, key=None, id=None, help=None, style=None):
    """Cards. items: (label, unit) or (label, unit, help). data: {"values": {label: "value" | {"value", "note"?, "help"?}}}
    style="metric": compact metric cards (muted label, 20px value, small note / unit line). The "?" sits after each label."""
    its = []
    for label, unit, *rest in items:
        it = {"label": label, "unit": unit}
        if rest and rest[0]:
            it["help"] = rest[0]
        its.append(it)
    return _extras({"type": "cards", "title": title, "items": its, "key": key, "id": id}, help=help, style=style)
def figure(label, key=None, id=None, help=None, span=None):
    """span="half": a structure drawing or other compact figure; the page puts it next to another block on wide screens."""
    return _extras({"type": "figure", "label": label, "key": key, "id": id}, help=help, span=span)
def compare(title=None, id=None, help=None):
    """One card per compound, for side-by-side comparison (Structure & formats).
    data: {"compounds": [{"label", "formula"?, "svg"?: square 2D SVG, "mol3d"?: molblock}, ...]} (1-3 items).
    Each card shows the 2D depiction and the 3D viewer in equal square boxes; one compound puts them side by side."""
    return _extras({"type": "compare", "title": title, "id": id}, help=help)
def code(label, title=None, key=None, id=None, help=None):
    return _extras({"type": "code", "title": title, "label": label, "key": key, "id": id}, help=help)
def formats(label, title=None, key=None, id=None, help=None):
    return _extras({"type": "formats", "title": title, "label": label, "key": key, "id": id}, help=help)
def notes(title, key=None, id=None, help=None, collapsed=False):
    """Notes. data: {"items": [str, ...]}.  collapsed=True: no visible list; the items sit in a "?" next to the title
    (a small "Details" disclosure when title is None)."""
    return _extras({"type": "notes", "title": title, "key": key, "id": id}, help=help, collapsed=collapsed)
def downloads(title, id=None, help=None):
    """Downloadable text files, each with a name, content and MIME type."""
    return _extras({"type": "downloads", "title": title, "id": id}, help=help)
def button(label, key=None):
    return {"type": "button", "label": label, "key": key}
def head(title, key=None, help=None):
    return _extras({"type": "heading", "title": title, "key": key}, help=help)
def portal(title=None, id=None, help=None, key=None):
    """Grid of status tiles (2-4 per row, responsive).
    data: {"tiles": [{"label", "status": "yes"|"partial"|"no", "value"?: short text, "chips"?: [str], "help"?: str|list}]}
    yes = teal tint, partial = amber tint, no = quiet grey (no fill)."""
    return _extras({"type": "portal", "title": title, "key": key, "id": id}, help=help)
def chips(title=None, id=None, help=None, key=None):
    """Row of coloured status pills.
    data: {"chips": [{"label", "level": good|possible|poor, "value"?: str, "help"?: str|list}]}
    level aliases: strong / yes = good, moderate / partial = possible, weak / no = poor (green / amber / grey).
    A chip with help is clickable: hover or focus shows its popover, click pins it."""
    return _extras({"type": "chips", "title": title, "key": key, "id": id}, help=help)
def reaction(title=None, id=None, help=None, key=None):
    """Reaction scheme next to the parent structure; all four CIL channel options show the complete substrate and product instead.
    data: {"svg"?: fallback parent structure SVG with the reactive groups highlighted,
           "options": [{"label": "CIL · dansyl", "steps": [{"from": "R–NH₂", "reagent": "Dns-Cl", "to": "R–NH–Dns", "count"?: 2}],
                        "formula"?, "mass"?, "heavy"?, "help"?: str|list, "svg"?: parent structure with only this option's sites highlighted,
                        "product_svg"?: complete product structure; when present, render substrate —reagent→ product instead of functional-group steps}]}
    Options render as selectable pills (the first is selected); the selected option selects its own SVG when present and draws either full structures or each step as from —reagent→ to ×count.
    A small line shows "formula · mass · heavy mass" + the option's "?". mass / heavy: bare value (a number is shown with 4 decimals),
    the renderer adds the "mass" / "heavy mass" labels. Strings may use unicode subscripts or <sub> / <sup>."""
    return _extras({"type": "reaction", "title": title, "key": key, "id": id}, help=help)
def groupmap(title=None, id=None, help=None, key=None):
    """Functional groups drawn on the input structure.
    data: {"svg": input structure, "found": [{"family", "name", "count", "smarts", "atoms": "3, 2; 5, 4", "svg": the group alone,
           "all_svg"?: every site tinted on the input (only with several sites),
           "sites": [{"atoms": "3, 2", "svg": input with this site tinted and its atoms numbered}], "more"?: sites without a picture}]}
    A list of groups (grouped by family) next to one large structure; choosing a group tints its sites, choosing a site shows that one with atom numbers."""
    return _extras({"type": "groupmap", "title": title, "key": key, "id": id}, help=help)
def grouptable(columns, title=None, id=None, help=None, caption=None, key=None, row_help=False):
    """Collapsible groups of rows that share the same columns (collapsed by default; chevron or header click toggles; keyboard accessible).
    data: {"groups": [{"title", "meta"?: short text on the right e.g. "0.66 · 3 methods", "help"?: str|list, "rows": [row, ...]}], "count"?: int}
    rows use the same cell / row formats as table (row_help=True adds the "?" column)."""
    return _extras({"type": "grouptable", "title": title, "columns": columns, "caption": caption, "key": key, "id": id}, help=help, row_help=row_help)


# ---- 模块 ----
def module(id, title, cls, one_liner, inputs, outputs, group="modules", subgroup=None, live=False, example=None, network=None, based_on=None):
    """network: a sentence saying what the page sends online; the frontend asks first when Internet search is off.
    based_on: short names of what the page relies on, the most important first, at most 4 (external engines / models / data, then "mzDesk" for its own rules);
    the frontend writes "Based on …" at the foot of the page."""
    m = {"id": id, "title": title, "cls": cls, "group": group, "subgroup": subgroup, "live": live,
         "one_liner": one_liner, "inputs": inputs, "outputs": outputs, "example": example}
    if network:
        m["network"] = network
    if based_on:
        m["based_on"] = list(based_on)
    return m


# 示例输入（与前端 initialValues 的形状一致）
def _smiles(s):
    return {"mode": "SMILES / InChI", "value": s}
def _name(q):
    """Name / ID / Formula mode (resolved online via PubChem, then cached). Examples use names: the page's first input mode."""
    return {"mode": "Name / ID / Formula", "value": q}
_formula = _name
def _prec(mz="", ion_mode="Positive", adduct="Auto"):
    return {"precursor_mz": mz, "ion_mode": ion_mode, "adduct": adduct}
def _paste(peaks):
    return {"source": "Paste", "value": peaks}


CAFFEINE = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
CAFFEINE_MS2 = "42.0338 5\n56.0495 4\n69.0447 6\n83.0604 9\n110.0713 18\n138.0662 100\n195.0877 32"
CAFFEINE_MS1 = "195.0877 100\n196.0902 10.3\n197.0922 0.9"
# Structures for the other page examples (PubChem canonical / isomeric SMILES; tests resolve the names from EXAMPLE_STRUCTURES offline).
ASPIRIN = "CC(=O)OC1=CC=CC=C1C(=O)O"
TESTOSTERONE = "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C"
KAEMPFEROL = "C1=CC(=CC=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O"
QUININE = "COC1=CC2=C(C=CN=C2C=C1)[C@H]([C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)O"
METFORMIN = "CN(C)C(=N)N=C(N)N"
SULFAMETHOXAZOLE = "CC1=CC(=NO1)NS(=O)(=O)C2=CC=C(C=C2)N"
ACETAMINOPHEN = "CC(=O)NC1=CC=C(C=C1)O"
PROPRANOLOL = "CC(C)NCC(COC1=CC=CC2=CC=CC=C21)O"
CHOLESTEROL = "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C"
# Metformin [M+H]+ isotope pattern at unit-mass resolution (IsoSpecPy; M = 130.1087).
METFORMIN_MS1 = "130.1087 100\n131.1103 6.3\n132.1116 0.2"
# MassBank MSBNK-ACES_SU-AS000877: carbamazepine [M+H]+, LC-ESI-QFT, ramp 20-70 % (CC BY 4.0).
# m/z as recorded, intensity as relative intensity (base peak = 999); the MS1 pattern is IsoSpecPy at Orbitrap 60k.
CARBAMAZEPINE_MASSBANK_MS2 = ("62.3173 3\n75.9558 3\n94.357 3\n101.0235 4\n104.4793 3\n107.373 3\n143.0337 7\n165.0704 5\n"
                               "167.0846 9\n175.0603 6\n177.0682 5\n179.0733 15\n191.073 4\n192.0805 201\n193.0891 45\n"
                               "194.0967 809\n195.1005 63\n220.0767 8\n221.0656 4\n237.1025 999\n238.1051 80")
CARBAMAZEPINE_MS1 = "237.1022 100\n238.0993 0.7\n238.1056 16.6\n239.1082 1.6"
# MassBank MSBNK-ACES_SU-AS000069: sulfamethoxazole [M-H]-, LC-APCI-QFT, ramp 20-70 % (CC BY 4.0).
# m/z as recorded, intensity as relative intensity (base peak = 999).
SULFAMETHOXAZOLE_MASSBANK_MS2 = ("63.9626 319\n92.0507 63\n156.0129 999\n157.0136 15\n157.0166 54\n158.0089 20\n"
                                  "213.8455 9\n213.8606 11\n252.0448 362\n253.0486 37")
# MassBank 里的真实咖啡因 [M+H]+ MS/MS（MSBNK-LCSB-LU119903，Q Exactive，CE 45；NIST 格式的 0-999 强度）：Library search 页的示例
CAFFEINE_MASSBANK_MS2 = ("56.0496 2\n69.0447 18\n83.0603 10\n97.0076 2\n97.0395 1\n109.0397 3\n110.0713 101\n111.0554 5\n112.0873 1\n122.0711 1\n"
                         "123.0428 3\n135.0668 1\n136.0747 1\n138.0663 999\n151.0978 5\n156.0768 2\n180.0642 2\n181.0719 3\n195.0878 847")
# MS-FINDER 页的示例：MassBank 记录 MSBNK-LCSB-LU119905 的真实咖啡因 [M+H]+ MS/MS（Q Exactive，CE 75，CC BY 4.0；作者 Elapavalore, Kondić, Singh, Schymanski）；
# m/z 取记录里的 4 位小数，强度取记录的 rel.int.（0-999）
CAFFEINE_MSFINDER_MS2 = ("54.0339 7\n56.0495 88\n57.0448 9\n67.029 9\n68.0369 7\n69.0446 320\n70.0286 4\n80.0129 2\n81.0447 38\n82.0525 16\n"
                         "83.0603 316\n93.0447 7\n97.0075 27\n97.0397 21\n108.0557 14\n109.0396 90\n110.0712 590\n111.0552 39\n113.0346 4\n"
                         "122.0711 3\n123.0427 144\n135.0663 5\n136.0744 3\n138.0662 999\n154.0607 4\n156.0762 2\n195.0876 58")
# Example names and the structure PubChem gives for each (checked by InChIKey). Tests use this table to run the examples offline.
EXAMPLE_STRUCTURES = {
    "caffeine": CAFFEINE,
    "aspirin": ASPIRIN,
    "testosterone": TESTOSTERONE,
    "kaempferol": KAEMPFEROL,
    "quinine": QUININE,
    "metformin": METFORMIN,
    "sulfamethoxazole": SULFAMETHOXAZOLE,
    "acetaminophen": ACETAMINOPHEN,
    "propranolol": PROPRANOLOL,
    "cholesterol": CHOLESTEROL,
    "L-tyrosine": "N[C@@H](Cc1ccc(O)cc1)C(=O)O",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "lidocaine": "CCN(CC)CC(=O)Nc1c(C)cccc1C",
    "POPC": r"CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\CCCCCCCC",
    "D-glucose": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
    "crisaborole": "B1(C2=C(CO1)C=C(C=C2)OC3=CC=C(C=C3)C#N)O",     # PubChem CID 44591583，C14H10BNO3
    "captopril": "C[C@H](CS)C(=O)N1CCC[C@H]1C(=O)O",
    "DOPG": r"CCCCCCCC/C=C\CCCCCCCC(=O)OC[C@H](COP(=O)(O)OCC(O)CO)OC(=O)CCCCCCC/C=C\CCCCCCCC",   # PG 36:2
}
PG_36_2 = r"CCCCCCCC\C=C/CCCCCCCC(=O)OCC(COP(O)(=O)OCC(O)CO)OC(=O)CCCCCCC\C=C/CCCCCCCC"
PG_36_2_MS2 = ("139.00977 3\n152.99622 24\n168.04965 3\n279.22757 48\n280.2258 27\n281.24707 1101\n"
               "282.24652 24\n283.26199 14\n491.25638 10\n509.2948 17\n599.81622 3\n713.46429 10\n773.53589 416")
# GNPS / MSnLib L-leucine [M+H]+, CCMSLIB00012361594. Retrieved with its public USI.
LEUCINE_GNPS_MS2 = ("41.03854 39.840\n41.04063 1.207\n41.04158 0.942\n42.03378 2.887\n43.05418 39.781\n43.05638 1.101\n"
                    "43.05739 0.818\n44.04942 100\n44.05167 2.780\n44.0528 2.151\n44.05395 0.824\n45.03341 1.259\n"
                    "55.01777 1.415\n55.05413 2.384\n56.04937 1.796\n57.0572 0.609\n57.93464 0.994\n58.06499 2.649\n"
                    "59.93002 0.700\n61.92777 0.641\n62.98968 1.084\n65.03836 0.513\n68.04931 0.550\n69.06963 1.236\n"
                    "73.06475 0.644\n86.05969 0.677\n86.09624 5.773\n97.00741 1.578\n115.96404 0.856\n116.97169 0.973\n117.95972 0.536")
_MASSQL_EXAMPLE = base64.b64encode(b"BEGIN IONS\nTITLE=example\nPEPMASS=132.102\nSCANS=42\n44.04942 100\n86.09624 5.773\nEND IONS\n").decode("ascii")


MODULES = [
    module("gnps_fastmasst", "Spectrum in Which Set", "A", "MS2 spectrum → matching public datasets and reference library hits",
        [dict(spectrum("spectrum", "MS2 spectrum"), required=True, section="Query"),
         dict(number("precursor_mz", "Precursor m/z", step=0.0001), header_fallback=True, section="Query"),
         dict(text("name", "Query name"), section="Query"),
         dict(select("library", "Search", ["both", "repository", "library"], "both"), section="Search"),
         dict(number("cosine", "Cosine threshold", 0.7, 0.01), section="Search"),
         dict(number("min_peaks", "Minimum matched peaks", 3, 1), section="Search"),
         dict(number("pm_tol", "Precursor tolerance (Da)", 0.05, 0.001), section="Search"),
         dict(number("frag_tol", "Fragment tolerance (Da)", 0.05, 0.001), section="Search"),
         dict(number("query_min_rel", "Query noise threshold (% base peak)", 0, 0.1), section="Search"),
         dict(text("diagnostic", "Diagnostic ions (m/z)", placeholder="86.0964, 44.0495"), section="Optional verification"),
         dict(text("diag_required", "Required ions (m/z)", placeholder="86.0964"), section="Optional verification"),
         dict(text("diag_rule", "Rule: all, any, k, or k-of-n", "all"), section="Optional verification"),
         dict(number("diag_tol", "Diagnostic tolerance (Da)", 0.01, 0.001), section="Optional verification"),
         dict(number("diag_min_rel", "Minimum intensity (% base peak)", 0, 0.1), section="Optional verification"),
         dict(number("verify_max", "Max spectra to verify per search", 20, 1), section="Optional verification"),
         dict(checkbox("redu_join", "Join Pan-ReDU metadata"), section="Pan-ReDU"),
         dict(select("redu_bg", "Background", ["pos_lc_or_unknown", "pos_lc", "neg_lc_or_unknown", "neg_lc", "all"]),
              section="Pan-ReDU", visible_when={"field": "redu_join", "equals": True},
              help="Choose the background for your ion mode.")],
        [table(["Search", "Passing hits", "Before filter", "Best cosine", "Task ID"], "Searches", id="searches"),
         downloads("Full results", id="downloads"),
         table(["Search", "Dataset", "Hits"], "Top datasets", id="datasets", page_size=10, compact=True),
         table(["Name", "Cosine", "Matched peaks", "Adduct", "Accession", "Entries"], "Top reference-library names", id="reference_hits", compact=True),
         table(["Name", "Cosine", "Matched peaks", "Adduct", "Accession", "Entries"], "Top suspect names", id="suspect_hits",
               help="Suspect names are propagated annotations.", compact=True),
         table(["Search", "Dataset", "USI", "Cosine", "Matched peaks", "Name", "Adduct", "Verification", "Diagnostic pass", "ReDU match"], "Passing hits", id="hits", page_size=10, compact=True),
         table(["Search", "Candidates", "Verified", "Diagnostic pass", "Not selected", "Fetch failed"], "Diagnostic verification", id="verification"),
         table(["Search", "Hit set", "Column", "Value", "Hit files", "Hit datasets", "File enrichment", "Dataset enrichment", "Dataset-weighted enrichment"], "Pan-ReDU enrichment", id="redu", page_size=10, compact=True)], subgroup="gnps", live=True,
        example={"spectrum": _paste(LEUCINE_GNPS_MS2), "precursor_mz": "132.102", "name": "L-leucine"},
        network="The spectrum is sent to GNPS2 fastMASST, and hits are verified through the public GNPS USI service.",
        based_on=["fastMASST", "GNPS", "Pan-ReDU"]),

    module("gnps_redu", "Dataset Finder", "B", "Sample metadata + polarity → matching public datasets",
        [text("organism", "Organism name or NCBI taxid", placeholder="Mus musculus or 10090"),
         text("body_part", "Body part", placeholder="brain"),
         dict(textarea("where", "Other filters", placeholder="SampleType=animal", rows=3),
              help="One column=value per line; commas separate alternatives."),
         number("min_files", "Minimum matching files per dataset", 1, 1),
         select("polarity", "Required polarity", ["any", "positive", "negative", "both"])],
        [cards([("Cache date", ""), ("Cached rows", ""), ("Matched files", ""), ("Passed datasets", ""), ("Undetermined", "")], "Summary", id="summary", style="metric"),
         downloads("Full results", id="downloads"),
         table(["Status", "Dataset", "Source", "Files", "Positive", "Negative", "Alternating", "Unknown", "Other"], "Datasets", id="datasets", page_size=10, compact=True),
         table(["Column", "Value", "Files"], "Matched metadata values", id="values", page_size=10, compact=True)],
        subgroup="gnps", live=True, example={"organism": "Mus musculus", "body_part": "brain", "min_files": 13, "polarity": "both"},
        based_on=["Pan-ReDU", "GNPS"]),

    module("gnps_massql", "Fragment Pattern", "A", "MS file + fragment or neutral-loss query → matching MS2 spectra",
        [dict(textarea("query", "MassQL query", "QUERY scaninfo(MS2DATA) WHERE MS2PROD=44.04942:TOLERANCEMZ=0.01", rows=4), required=True),
         dict(file("input_file", "mzML / mzXML / MGF file", ".mzml,.mzxml,.mgf"), required=True, binary=True, max_bytes=20_000_000,
              help="Processed on the mzDesk server, not sent to GNPS2.")],
        [cards([("Input file", ""), ("MS2 spectra", ""), ("Matches", "")], "Summary", id="summary", style="metric"),
         downloads("Full results", id="downloads"),
         table(["Scan", "Precursor m/z", "RT (min)", "Polarity", "Charge", "Peaks", "Title", "Matched conditions"], "Matching spectra", id="matches", page_size=10, compact=True)], subgroup="gnps", live=True,
        example={"query": "QUERY scaninfo(MS2DATA) WHERE MS2PROD=44.04942:TOLERANCEMZ=0.01",
                 "input_file": {"name": "leucine-example.mgf", "base64": _MASSQL_EXAMPLE}},
        based_on=["MassQL", "pyOpenMS"]),

    module("identity", "Identity", "B", "Compound → names, IDs, classification",
        [structure(), select("kind", "Input type", IDENTITY_KINDS)],
        [cards([("Name", ""), ("Formula", ""), ("PubChem CID", "")], "Resolved", id="resolved", style="metric"),
         figure("2D depiction", id="depiction", span="half"),
         table(["Name type", "Value"], "Names", id="names"),
         table(["Database", "ID"], "Database IDs", id="database_ids"),
         table(["PubChem CID", "Title"], "Other matches for this name", id="candidates"),
         table(["Synonym"], "Synonyms", id="synonyms"),
         table(["System", "Level", "Term"], "Classification", id="classification", help="Each system classifies on its own."),
         table(["Origin", "Evidence", "Source"], "Origin", id="origin"),
         table(["Pathway", "Database", "Link"], "Pathways", id="pathways"),
         table(["Biofluid", "Normal concentration", "Unit", "Source"], "Biofluid concentrations", id="biofluids"),
         table(["Species", "Source"], "Species", id="species"),
         table(["Vendor", "CAS", "Purity"], "Standards & availability", id="standards"),
         table(["Source", "Signal", "Hazard statements"], "GHS safety", id="ghs")],
        subgroup="molecular", live=True, example={"structure": _name("aspirin"), "kind": "Auto"},
        network="Identity looks the compound up in PubChem and UniChem, then loads classification and background from public sources.",
        based_on=["PubChem", "UniChem", "LOTUS", "ClassyFire"]),

    module("structure_formats", "Structure & Formats", "B", "Up to 3 structures → 2D / 3D views, SMILES, InChI, InChIKey, MOL",
        [dict(structure(label="Compound 1"), required=True, label_field=True),
         {**structure("structure_2", "Compound 2", optional=True), "label_field": True},
         {**structure("structure_3", "Compound 3", optional=True), "label_field": True}],
        [compare("Compounds", id="compare", help="3D: one low-energy conformer."),
         formats("Structure in every format", "Structure formats", id="formats")],
        subgroup="molecular", live=True, example={"structure": _name("testosterone")}, based_on=["RDKit"]),

    module("properties", "Properties & Groups", "B", "Compound → formula, masses, descriptors, functional groups",
        [structure()],
        [cards([("Formula", ""), ("Monoisotopic mass", "Da"), ("Average MW", "Da"), ("Net charge", "formal")], "Composition", id="summary", style="metric"),
         table(["Group", "Descriptor", "Value"], "Descriptors", id="descriptors"),
         groupmap("Functional groups", id="functional_groups",
                help=["Click a group to see its sites on the structure.", "Overlapping groups are all listed."])],
        subgroup="molecular", live=True, example={"structure": _name("kaempferol")}, based_on=["RDKit", "mzDesk"]),

    module("lipid_short_name", "Lipid Short Name", "B", "Lipid shorthand → standardized names and an illustrative structure",
        [dict(text("name", "Lipid name", placeholder="PC 16:0/18:1"), required=True)],
        [cards([("L1", "Sum composition"), ("L2", "Chains"), ("L3", "Chains with groups")], "LIPID+", id="lipid_plus", style="metric",
               help="No sn positions. — = the name does not give that level."),
         cards([("Species", "Sum composition"), ("Molecular species", "Chains")], "LIPID MAPS", id="lipid_maps", style="metric",
               help="— = cannot be written from this input; the ? says why."),
         figure("Illustrative structure", id="structure"),
         notes(None, id="no_structure"),
         notes("Assumptions", id="assumptions", collapsed=True),
         notes("Flags", id="flags", collapsed=True)],
        subgroup="molecular", live=True, example={"name": "PC 16:0/18:1"}, based_on=["LIPID+", "RDKit"]),

    # ---- Detectability：7 页。只有 Platform 的脂质组 / 代谢组部分接了引擎，其余都是占位 ----
    module("platform_fit", "Platform", "A", "Compound → suitable platforms: lipidomics, metabolomics, CIL Metabolomics, GC-MS",
        [structure()],
        [portal("Platform", id="platform",
                help=["Yes: expected to work. Partial: after derivatization. No: not expected.", "The ? on a tile gives the reason."]),
         reaction("Why and how", id="reactions",
                  help=["Pick a reagent to highlight its reactive groups.", "Masses are monoisotopic."])],
        subgroup="detectability", live=True, example={"structure": _name("L-tyrosine")}, based_on=["RDKit", "mzDesk"]),

    module("separation", "Separation", "B", "Compound + pH → charge, logD, suitable separation modes",
        [structure(), slider("ph", "Mobile phase pH", 1, 13, 0.1, 2.7)],
        [cards([("logD at pH", "", "logP corrected for the neutral fraction at this pH."),
                ("Dominant charge state", "", "Most probable net charge at this pH."),
                ("Fraction neutral", "", "Share of molecules with no charge.")],
               "pH behaviour", id="ph_summary", style="metric"),
         figure("logD and net charge vs pH", id="ph_curves", help="Solid: logD (left axis). Dashed: net charge (right axis). Dotted: chosen pH."),
         figure("Ionizable sites on the structure", id="sites_structure", span="half", help="Coral: acid. Teal: base. Violet: both. Label: atom and pKa."),
         chips("Separation mode", id="modes",
               help="Green: good. Amber: possible. Grey: poor or not needed. Hover for the reason."),
         table(["Site", "pKa", "Type"], "Ionizable sites", id="sites",
               help=["Site: element and atom number (from 1).", "Acid loses H; base takes H."]),
         notes("Notes", id="rule_hits", collapsed=True)],
        subgroup="detectability", live=True, example={"structure": _name("ibuprofen"), "ph": 2.7}, based_on=["MolGpKa", "RDKit", "mzDesk"]),

    module("ionization", "Ionization", "B", "Compound + pH → ESI polarity and EI usability",
        [structure(), slider("ph", "Mobile phase pH", 1, 13, 0.1, 2.7)],
        [chips("Ionization", id="chips",
               help=["Polarity: the better ESI mode (+, −, both); poor = weak in both.", "Hover a chip for the reason."]),
         figure("Ionizable sites", id="sites_structure", span="half", help="Coral: acid. Teal: base. Violet: both. Label: atom and pKa."),
         notes("Rule hits", id="rule_hits", collapsed=True)],
        subgroup="detectability", live=True, example={"structure": _name("lidocaine"), "ph": 2.7}, based_on=["MolGpKa", "NEIMS", "mzDesk"]),

    module("adducts", "Adducts", "B", "Compound + ion mode + additive → expected ions, strongest first",
        [structure(), radio("ion_mode", "Ion mode", ["Positive", "Negative"]),
         select("additive", "Mobile phase additive", ["Formic acid", "Ammonium formate", "Ammonium acetate", "Acetic acid", "None"]),
         checkbox("alkali", "Sodium / potassium present", True)],
        [cards([("Dominant ion", ""), ("m/z", "")], "Expected base ion", id="base", style="metric"),
         table(["Ion", "Formula", "m/z", "Likelihood"], "Ions", id="ions", row_help=True,
               help="Ranked by expected intensity (high / medium / low)."),
         notes("Rule hits", id="rule_hits", collapsed=True)],
        subgroup="detectability", live=True,
        example={"structure": _name("POPC"),
                 "ion_mode": "Negative", "additive": "Ammonium formate", "alkali": True},
        based_on=["MolGpKa", "mzDesk"]),

    module("non_ms", "Non-MS Detection", "B", "Compound → suitability for UV, fluorescence, CAD / ELSD / RI",
        [structure(), select("mobile_phase", "Mobile phase", ["Water/acetonitrile + phosphate", "Water/methanol + formic acid"])],
        [portal("Detectors", id="detectors",
                help=["yes = good, partial = possible, no = poor.", "UV/PDA: chromophore strength and λmax. FLD: native fluorescence."]),
         figure("Chromophores", id="chromophore_structure", span="half", help="Blue: strong. Amber: weak. Numbers are atom numbers (from 1)."),
         notes("Rule hits", id="rule_hits", collapsed=True)],
        subgroup="detectability", live=True, example={"structure": _name("quinine"), "mobile_phase": "Water/methanol + formic acid"},
        based_on=["RDKit", "mzDesk"]),

    module("derivatization", "Derivatization", "B", "Compound + reagents → reactive groups, derivative formula and mass",
        [structure(),
         multiselect("reagents", "Reagents",
                     ["GC · TMS (MSTFA / BSTFA)", "GC · MeOx + TMS", "GC · TBDMS (MTBSTFA)", "GC · Acylation (TFAA)",
                      "CIL · Dansyl chloride (amine / phenol)", "CIL · Dansyl chloride, base-activated (hydroxyl)",
                      "CIL · DmPA bromide (carboxyl)", "CIL · Dansylhydrazine (carbonyl)",
                      "UV/FLD · OPA + 2-mercaptoethanol (primary amine)", "UV/FLD · FMOC-Cl (amine)", "UV/FLD · DNPH (carbonyl)"],
                     ["GC · TMS (MSTFA / BSTFA)", "GC · MeOx + TMS", "CIL · Dansyl chloride (amine / phenol)",
                      "CIL · Dansyl chloride, base-activated (hydroxyl)", "CIL · DmPA bromide (carboxyl)", "CIL · Dansylhydrazine (carbonyl)"])],
        [table(["Reagent", "Labels", "Formula", "Mass", "Heavy mass"], "Derivatives", id="deriv_rows", row_help=True,
               help=["Labels = reagent groups added.", "Heavy mass: ¹³C₂-labelled reagent (CIL only)."]),
         reaction("Reactions", id="deriv_reactions",
                  help=["×n = number of groups.", "Pick a reagent to highlight its reactive groups."])],
        subgroup="detectability", live=True, example={"structure": _name("D-glucose")}, based_on=["RDKit", "mzDesk"]),

    module("known_methods", "Known Methods", "B", "Compound → published methods, or those of similar compounds",
        [structure(), multiselect("sources", "Sources", ["Japanese Pharmacopoeia", "RepoRT"])],
        [cards([("Match", "", "Exact: the compound itself. Similar: only similar compounds. None: neither."),
                ("Methods found", "", "Exact matches plus the methods of similar compounds."),
                ("Sources", "")], "Summary", id="summary", style="metric"),
         table(["Source", "Test / dataset", "Column", "Mobile phase", "Detector", "RT (min)"], "Exact matches", id="exact", row_help=True,
               help=["Same compound; stereo is noted per row.", "? on a row: reference and conditions. RT: RepoRT only."]),
         grouptable(["Source", "Test / dataset", "Column", "Mobile phase", "Detector", "RT (min)"], "Similar compounds", id="similar", row_help=True,
                    help="Click a compound to see its methods.")],
        subgroup="detectability", live=True, example={"structure": _name("captopril"), "sources": ["Japanese Pharmacopoeia", "RepoRT"]},
        based_on=["Japanese Pharmacopoeia", "RepoRT", "RDKit"]),

    # ---- MS1 ----
    module("ms1_search", "Mass Search", "A", "Precursor m/z → matching database compounds",
        [precursor(), number("tol_ppm", "Mass tolerance (ppm)", 5.0, 0.1), multiselect("databases", "Exact-mass databases", DATABASES, ["HMDB", "PubChem"]),
         ms1(), select("resolution", "Resolution", list(RESOLUTIONS)), number("isotope_min", "Isotope score threshold", 0.8, 0.05)],
        [cards([("Compounds", "", "Counted per PubChem CID."),
                ("Formulas", ""), ("Closest hit", ""), ("Databases searched", "")],
               "Summary", id="summary", style="metric"),
         grouptable(["Name", "IDs", "PubChem CID"], "Exact-mass hits", id="hits",
                    help=["One group per formula and adduct, closest first.", "Groups below the isotope threshold are hidden."])],
        subgroup="ms1", live=True,
        example={"precursor": _prec("130.1087"), "tol_ppm": "5", "databases": ["HMDB", "PubChem"],
                 "ms1": METFORMIN_MS1, "resolution": "Unit mass", "isotope_min": "0.8"},
        network="Mass Search queries PubChem through NCBI E-utilities and PUG REST; every database is searched in its PubChem copy.",
        based_on=["PubChem", "IsoSpecPy"]),

    module("ms1_theory", "Theoretical", "A", "Compound or formula + adduct → precursor m/z, isotope pattern",
        [structure(label="Compound or formula"),
         precursor(label="Adduct", adducts={"Positive": ADDUCTS_POS[1:], "Negative": ADDUCTS_NEG[1:]}, mz=False),
         select("resolution", "Resolution", list(RESOLUTIONS))],
        [cards([("Precursor m/z", "theoretical", "Ion with the most abundant isotope of each element (M)."),
                ("Ion formula", ""), ("Most abundant peak", "")], "Theory", id="theory", style="metric"),
         table(["Peak", "m/z", "Relative abundance (%)"], "Isotope pattern", id="pattern",
               help=["M uses the most abundant isotope of every element.", "(a), (b): peaks resolved at this resolution."]),
         figure("Theoretical isotope pattern", id="pattern_plot")],
        subgroup="ms1", live=True, example={"structure": _name("crisaborole"), "precursor": _prec("", adduct="[M+H]+"), "resolution": "Orbitrap 60k"},
        based_on=["IsoSpecPy", "mzDesk"]),

    # ---- Spectrum（subgroup id 仍是 ms2）----
    module("ms2_clean", "Spectrum Cleanup", "B", "MS2 spectrum → cleaned spectrum in every common format",
        [spectrum("spectrum", "MS2 spectrum (any common format)"), number("precursor_mz", "Precursor m/z (optional)", step=0.0001),
         radio("ion_mode", "Ion mode (for exported headers)", ["Positive", "Negative"]),
         number("noise_pct", "Noise threshold (% of base peak)", 0, 0.1), select("normalize", "Normalize", NORMALIZE),
         number("precursor_tol", "Precursor tolerance (Da)", 0.02, 0.001), select("sort_by", "Sort by", SORT_BY),
         checkbox("remove_precursor", "Remove precursor")],
        [cards([("Detected format", ""), ("Peaks in", ""), ("Peaks out", ""), ("Removed", "noise · precursor")], id="summary", style="metric"),
         figure("Cleaned spectrum", id="spectrum", help="Grey: removed peaks."),
         table(["m/z", "Intensity", "Rel. intensity (%)", "Kept", "Reason"], "Peaks", "Rel. intensity: % of the input base peak. Reason: noise, precursor", id="peaks"),
         formats("Cleaned spectrum in every format", "Formats", id="formats")],
        subgroup="ms2", live=True,
        example={"spectrum": _paste(SULFAMETHOXAZOLE_MASSBANK_MS2), "precursor_mz": "252.0448", "ion_mode": "Negative",
                 "noise_pct": "1", "normalize": "Base peak = 100", "remove_precursor": True},
        based_on=["mzDesk"]),

    module("ms2_similarity", "Similarity", "A", "Two MS2 spectra → similarity scores and a mirror plot",
        [spectrum("query", "Query spectrum"), number("query_precursor", "Query precursor m/z (optional)", step=0.0001),
         spectrum("reference", "Reference spectrum"), number("reference_precursor", "Reference precursor m/z (optional)", step=0.0001),
         number("tol_da", "Fragment tolerance (Da)", 0.02, 0.001), number("threshold", "Score threshold", 0.7, 0.05),
         radio("ion_mode", "Ion mode", ["Positive", "Negative"]),
         checkbox("deep", "Deep-learning scores", network=True)],
        [cards([("Primary score", "spectral entropy, weighted"), ("Metrics ≥ threshold", "n / N"), ("Matched peaks", "pairs"),
                ("Query matched", "%"), ("Reference matched", "%")], id="summary", style="metric"),
         figure("Mirror plot", id="mirror", help="Up: query. Down: reference. Dark: matched peaks."),
         table(["Category", "Metric", "Score", "Matched peaks", "≥ threshold", "Note"], "Scores", id="scores")],
        subgroup="ms2", live=True,
        example={"query": _paste("197.04643 10\n255.22775 92\n256.23184 31\n277.60822 10\n283.24814 20"),
                 "reference": _paste("197.04643 8\n255.22775 100\n256.23184 28\n277.60822 7\n283.24814 22")},
        based_on=["matchms", "ms_entropy", "MS2DeepScore", "Spec2Vec"]),

    module("ms2_library", "Library Search", "A", "MS/MS or EI spectrum → ranked library matches",
        [radio("spectrum_type", "Spectrum type", ["MS/MS", "EI"]), spectrum(label="Query spectrum"),
         precursor(label="Precursor m/z (MS/MS only)"),
         number("tol_ppm", "Precursor tolerance (ppm, MS/MS only)", 10.0, 0.1),
         number("tol_da", "Fragment tolerance (Da, MS/MS only)", 0.02, 0.001),
         multiselect("libraries", "Libraries", LIBRARIES),
         multiselect("metrics", "Metrics", LIBRARY_METRICS, LIBRARY_DEFAULT_METRICS, network_options=METRIC_NAMES.DEEP),
         number("threshold", "Score threshold", 0.7, 0.05)],
        [cards([("Best hit", "name"), ("Rank score", "ranking metric"), ("Hits ≥ threshold", "n"),
                ("Candidates scored", "n"), ("Libraries", "version")], id="summary", style="metric"),
         table(["Rank", "Library", "Accession", "Name", "Adduct", "Precursor m/z", "Δ ppm", "Matched peaks"] + LIBRARY_DEFAULT_METRICS + ["Instrument / CE"],
               "Library hits", id="hits",
               help=["Hit: any selected metric ≥ threshold.", "Ranked by the first selected metric; ties: more matched peaks."]),
         figure("Mirror: query vs top hit", id="mirror", help="Up: query. Down: top hit. Dark: matched peaks.")],
        subgroup="ms2", live=True,
        example={"spectrum_type": "MS/MS", "ms2": _paste(CAFFEINE_MASSBANK_MS2), "precursor": _prec("195.0877")},
        based_on=["MassBank", "GNPS", "matchms", "ms_entropy"]),

    module("ms1_formula", "Formula Annotation", "A", "Precursor m/z + isotopes + MS2 → ranked molecular formulas",
        [precursor(adducts=BUDDY_ADDUCTS), ms1(), spectrum(label="MS2 spectrum (optional)"),
         select("instrument", "Instrument", ["Orbitrap", "Q-TOF", "FT-ICR"]),
         text("formula", "Formula (optional)", placeholder="C8H10N4O2"),
         checkbox("halogen", "Include halogens (F, Cl, Br, I)")],
        [cards([("Top formula", ""), ("Adduct used", ""), ("MS1 isotope peaks used", "")], id="summary", style="metric"),
         table(["Rank", "Formula", "Theoretical m/z", "Δ ppm", "Probability", "FDR", "MS1 isotope similarity"], "Formula candidates", id="candidates"),
         figure("MS2 fragment formulas · top candidate", id="fragments_plot"),
         table(["m/z", "Rel. intensity (%)", "Fragment formula", "Δ ppm", "Neutral loss", "Other fits"], id="fragments")],
        subgroup="ms2", live=True,
        example={"precursor": _prec("237.1022"), "ms1": CARBAMAZEPINE_MS1, "ms2": _paste(CARBAMAZEPINE_MASSBANK_MS2)},
        based_on=["BUDDY"]),

    module("ms2_structure", "Structure Prediction", "A", "Precursor + MS2 spectrum → ranked candidate structures",
        [precursor(adducts=MSFINDER_ADDUCTS), spectrum(label="MS2 spectrum"), ms1(),
         checkbox("halogen", "Include halogens (F, Cl, Br, I)"), checkbox("silicon", "Include Si"),
         checkbox("extended", "Extended element range")],
        [cards([("Top candidate", ""), ("Top formula", ""), ("Adduct used", ""), ("Elements searched", ""), ("Candidates", "structures / formulas")], id="summary", style="metric"),
         table(["Rank", "Structure", "Name", "SMILES", "Formula", "Score", "Ontology", "InChIKey", "Databases"], "Structure candidates", id="structures",
               help=["Score: fragmentation fit to the MS2; higher is better.", "Databases: where it is listed (? shows IDs)."]),
         table(["Rank", "Formula", "Neutral mass", "Mass error (mDa)", "Score"], "Formula candidates", id="formulas",
               help=["Score: mass accuracy, isotope pattern and MS2 fit.", "Mass error (mDa) = observed − theoretical."])],
        subgroup="ms2", live=True,
        example={"precursor": _prec("195.0877"), "ms2": _paste(CAFFEINE_MSFINDER_MS2)},
        based_on=["MS-FINDER"]),

    module("ms2_denovo", "De Novo Structure", "A", "Precursor + MS2 spectrum → generated candidate structures",
        [precursor(label="Precursor m/z (positive mode only)", adducts=FRIGID_ADDUCTS), spectrum(), ms1(),
         text("formula", "Formula (optional)", placeholder="C8H10N4O2"),
         select("mode", "Mode", FRIGID_MODES), checkbox("halogen", "Include halogens (F, Cl, Br, I)")],
        [cards([("Top candidate", ""), ("Formula used", ""), ("Device", ""), ("Samples", "generated / shown")], id="summary", style="metric"),
         table(["Rank", "Structure", "SMILES", "Formula", "Score", "InChIKey (2D)"], "Generated candidates", id="candidates",
               help="Score: similarity to the fingerprint predicted from the spectrum.")],
        subgroup="ms2", live=True,
        example={"precursor": _prec("195.0877"), "ms2": _paste(CAFFEINE_MSFINDER_MS2), "formula": "C8H10N4O2"},
        based_on=["FRIGID", "BUDDY"]),

    module("ms2_ei", "EI-MS Prediction", "A", "Structure → predicted 70 eV EI spectrum",
        [structure()],
        [cards([("Ion", ""), ("Formula", ""), ("Molecular ion m/z", ""), ("Base peak m/z", "")], id="summary", style="metric"),
         figure("Predicted EI spectrum", id="spectrum", help="Unit mass resolution."),
         table(["m/z", "Rel. intensity (%)"], "Peaks", id="peaks"),
         formats("Predicted spectrum in every format", "Formats", id="formats")],
        subgroup="ms2", live=True, example={"structure": _name("cholesterol")}, based_on=["NEIMS"]),

    module("ms2_esi", "MS/MS Prediction", "A", "Structure + adduct → predicted MS/MS at 10 / 20 / 40 eV",
        [structure(), select("adduct", "Adduct", ESI_ADDUCTS, "[M+H]+")],
        [cards([("Engine", ""), ("Precursor m/z", ""), ("Ion formula", ""), ("Peaks", "")], id="summary", style="metric"),
         figure("Low energy (10 eV)", id="spectrum_low"), figure("Medium energy (20 eV)", id="spectrum_medium"),
         figure("High energy (40 eV)", id="spectrum_high"),
         table(["Energy", "m/z", "Rel. intensity (%)", "Fragment"], "Peaks", "Fragment: formula · SMILES", id="peaks"),
         formats("Predicted spectra (one record per energy)", "Formats", id="formats")],
        subgroup="ms2", live=True, example={"structure": _name("caffeine"), "adduct": "[M+H]+"}, based_on=["CFM-ID", "MSRB"]),

    module("ms2_fragments", "Fragment Explanation", "A", "Structure + MS2 spectrum → explained peaks and fragments",
        [structure(label="Candidate structure"), precursor(label="Precursor m/z (optional)", adducts=METFRAG_ADDUCTS), spectrum(),
         slider("tree_depth", "Tree depth (bond cleavages)", 1, 4, 1, 2),
         number("tol_abs", "Absolute tolerance (Da)", 0.002, 0.001), number("tol_ppm", "Relative tolerance (ppm)", 10.0, 1)],
        [cards([("Explained peaks", "n / N"), ("Explained intensity", "%"), ("Theoretical precursor m/z", ""), ("Ion formula", "")], id="summary", style="metric"),
         figure("Observed spectrum", id="spectrum"),
         table(["Observed m/z", "Rel. intensity (%)", "Kind", "Ion formula", "Theoretical m/z", "Δ ppm", "Neutral loss", "Fragment structure"], id="peaks")],
        subgroup="ms2", live=True,
        example={"structure": _name("DOPG"), "precursor": _prec("", "Negative", "[M-H]-"), "ms2": _paste(PG_36_2_MS2)},
        based_on=["MetFrag"]),

    module("exclusivity", "Alternative Explanations", "A", "Compound or formula + m/z → same-formula compounds, isobars, isomers",
        [structure(label="Compound or formula", optional=True), precursor(), number("tol_ppm", "Tolerance (ppm)", 5.0, 0.1),
         multiselect("databases", "Databases", DATABASES, ["HMDB", "PubChem"]),
         spectrum(label="MS2 spectrum (optional)"), number("tol_da", "Fragment tolerance (Da)", 0.01, 0.001)],
        [table(["InChIKey (first block)", "Name", "Databases", "Stereoisomers / isotopes collapsed"], "Same-formula candidates",
               id="same_formula",
               help=["One row per connectivity (stereoisomers folded).", "Databases: which selected databases list it."]),
         table(["Formula", "Name", "Database", "Exact m/z", "Δ ppm"], "Isobaric interferences", "Different formula, same m/z within tolerance",
               id="isobaric",
               help="Exact m/z and Δ ppm use the closest adduct."),
         table(["Isomer", "Fragments", "Source", "SMILES", "Explained peaks", "Distinguishable by MS2", "Distinguishing fragments", "Atoms covered"], "Isomers",
               id="isomers",
               column_help={
                   "Explained peaks": "Peaks explained, and their share of the intensity.",
                   "Distinguishable by MS2": "yes: a peak is explained for your compound but not for this isomer.",
                   "Distinguishing fragments": "Those peaks and the fragments that explain them.",
                   "Atoms covered": "Share of heavy atoms inside the coloured fragments.",
                   "Fragments": "View: explained peaks drawn on the structure."})],
        subgroup="ms2", live=True,
        example={"structure": _name("caffeine"), "precursor": _prec("195.0877"), "ms2": _paste(CAFFEINE_MSFINDER_MS2)},
        network="Name / ID inputs are resolved online. The same-formula and isobaric tables and isomer candidate search query PubChem with the formula or mass; MetFrag runs locally.",
        based_on=["PubChem", "MetFrag"]),

    # ---- Orthogonal ----
    module("ccs", "CCS", "A", "Compound + adduct + observed CCS → reference or predicted CCS, match check",
        [structure(), select("adduct", "Adduct", ADDUCTS_POS[1:] + ADDUCTS_NEG[1:]),
         text("observed", "Observed CCS (Å², optional)", placeholder="152.3 or 152.3, 160.1"),
         select("platform", "Platform", CCS_PLATFORMS),
         select("gas", "Drift gas", CCS_GASES),
         select("tw_calibration", "TWIMS calibration (TWIMS only)", CCS_TW_CALIBRATIONS)],
        [cards([("Reference CCS", "Å²", "The measured value closest to yours, else the median; a model prediction if none."),
                ("Source", "", "Measured: published N₂ value. Predicted: model estimate."),
                ("Observed", "Å²", "With several values, the closest to each reference is used."),
                ("Δ", "%", "(observed − reference) / reference, in %."),
                ("Verdict", "", "Consistent: within tolerance. Inconsistent: outside it. Not comparable: He, FAIMS / DMS.")],
               "Summary", id="summary", style="metric", help="Observed vs. the N₂ reference CCS of the same adduct."),
         table(["Source dataset", "CCS type", "Method", "CCS (Å²)", "Δ %", "Tolerance %", "Verdict"], "Reference values", id="reference", row_help=True,
               help=["Measured N₂ values for this compound and adduct.", "The ? on a row gives the publication."]),
         table(["Adduct", "CCS (Å²)", "CCS type", "Source dataset"], "Other ions of this compound", id="other_ions",
               help="Other adducts of this compound; not compared."),
         table(["Predicted CCS (Å²)", "Tolerance %", "Test median error %", "Test 95th percentile %", "Training domain"], "Prediction", id="prediction", row_help=True,
               help="Shown only when no measured value exists; the ? gives the model and its error."),
         notes("Assumptions", id="notes", help="Comparison rules, tolerances, model and data source.")],
        subgroup="orthogonal", live=True, example={"structure": _name("caffeine"), "adduct": "[M+H]+", "observed": "140.4", "platform": "DTIMS (single field)", "gas": "N₂"},
        based_on=["CCSbase", "mzDesk"]),

    module("mrm", "MRM Transitions", "B", "Pick <b>precursor–product pairs</b> and <b>collision energies</b> for a targeted method.",
        [radio("start_from", "Start from", ["Structure (predicted fragments)", "Library spectrum", "My spectrum"]), structure(optional=True), spectrum(label="MS2 spectrum (if starting from my spectrum)"),
         radio("ion_mode", "Ion mode", ["Positive", "Negative"]), select("instrument", "Instrument", ["Triple quadrupole (generic)", "QTRAP", "Q-TOF (PRM)", "Orbitrap (PRM)"]),
         number("n_transitions", "Transitions", 3, 1), select("ce_unit", "CE unit", ["eV", "NCE"])],
        [cards([("Polarity", ""), ("Transitions", "n")], id="summary", style="metric"),
         table(["Precursor m/z", "Product m/z", "CE", "Role", "Fragment formula", "Source", "Relative intensity"],
               caption="Role: quantifier / qualifier. Source: library spectrum or predicted", id="transitions"),
         button("Export method table (CSV)"),
         formats("Method table (CSV)", id="method_table"),
         notes("Assumptions", id="notes")], subgroup=None, live=True,
        example={"start_from": "My spectrum", "ion_mode": "Positive", "n_transitions": 3, "ce_unit": "eV",
                 "ms2": _paste("BEGIN IONS\nPEPMASS=195.0877\n" + CAFFEINE_MASSBANK_MS2 + "\nEND IONS")},
        based_on=["CFM-ID", "MassBank", "MetFrag"]),

    module("calculator", "Calculator", "B", "Masses, concentrations, retention times → ppm error, dilution, CE, Kovats RI",
        [heading("ppm error"), number("ppm_theoretical", "Theoretical m/z", step=0.0001), number("ppm_observed", "Observed m/z", step=0.0001),
         heading("Isotope labelling"), text("lab_smiles", "SMILES"), text("lab_label", "Label", placeholder="13C6, d4, 15N2"), select("lab_adduct", "Adduct", ADDUCTS_POS[1:]),
         heading("Concentration"), number("conc_value", "Value"), select("conc_from", "From", ["µM", "mM", "ng/mL", "µg/mL", "mg/L"]), select("conc_to", "To", ["µM", "mM", "ng/mL", "µg/mL", "mg/L"]),
         number("conc_mw", "Molar mass (g/mol)", step=0.01, placeholder="for mass ↔ molar"),
         heading("Dilution"), number("target_conc", "Target concentration"), number("final_volume", "Final volume (mL)", step=0.1),
         heading("Collision energy"), number("ce_nce", "NCE"), number("ce_ev", "eV"), number("ce_mz", "Precursor m/z", step=0.0001), number("ce_z", "Charge", 1, 1),
         heading("GC retention index (Kovats)"), grid("alkanes", "Alkane series", ["Carbon number", "RT"]), number("ri_rt", "Analyte RT", step=0.01)],
        [cards([("ppm error", "ppm")], "ppm error", id="ppm_error", style="metric"),
         table(["Isotopologue", "m/z", "Abundance"], "Isotope labelling", id="isotope_labels"), table(["Internal standard", "Analyte M+n", "Crosstalk %"], id="crosstalk"),
         cards([("Converted", "")], "Concentration", id="concentration", style="metric"), table(["Target conc.", "Stock conc.", "Final volume", "Stock volume", "Diluent volume"], id="dilution", caption="Dilution and stock preparation"),
         cards([("NCE ↔ eV", "")], "Collision energy", id="collision_energy", style="metric"), cards([("Kovats RI", "")], "GC retention index", id="kovats", style="metric"),
         notes("Assumptions", id="notes")], subgroup=None, live=True,
        example={"ppm_theoretical": 260.1645, "ppm_observed": 260.1651, "lab_smiles": PROPRANOLOL, "lab_label": "13C3", "lab_adduct": "[M+H]+",
                 "conc_value": 25.935, "conc_from": "mg/L", "conc_to": "µM", "conc_mw": 259.35, "target_conc": 10, "final_volume": 10,
                 "ce_nce": 30, "ce_mz": 260.1645, "ce_z": 1, "alkanes": "C8 8.00\nC9 9.00\nC10 10.00\nC11 11.00", "ri_rt": 9.50},
        based_on=["IsoSpecPy", "mzDesk"]),

    # ---- DMPK（非核心）----
    module("met_shift", "Mass-Shift Search", "B", "Search a peak or feature list for <b>phase I / II metabolites</b> of the parent.",
        [structure(label="Parent compound or formula"), precursor(label="Measured precursor m/z (optional; blank = theoretical)"),
         grid("features", "Peak / feature list", ["m/z", "RT", "Intensity"], "Paste or upload: one feature per line"), file("features_file", "Feature table (CSV, optional)", ".csv"),
         multiselect("shifts", "Mass shifts", ["+O (hydroxylation)", "−2H (dehydrogenation)", "+2H (reduction)", "−CH2 (demethylation)", "+O −CH2", "+C6H8O6 glucuronide", "+SO3 sulfate", "+C2H2O acetyl", "+C10H15N3O6S GSH", "+C2H3NO glycine"]),
         number("tol_ppm", "Tolerance (ppm)", 10.0, 0.1)],
        [table(["Shift", "Δ mass", "Expected m/z", "Found m/z", "Δ ppm", "RT", "Intensity", "Candidate reaction"], id="hits"),
         figure("Features with metabolite hits marked (m/z vs RT)", id="features_plot"), notes("Assumptions", id="notes")], subgroup="dmpk", live=True,
        example={"structure": _name("acetaminophen"), "precursor": _prec("152.0706", "Positive", "[M+H]+"),
                 "features": "m/z RT Intensity\n152.0706 3.2 100000\n168.0655 2.1 18000\n328.1027 2.8 9000\n232.0274 3.5 4000",
                 "shifts": ["+O (hydroxylation)", "−2H (dehydrogenation)", "+C6H8O6 glucuronide", "+SO3 sulfate"], "tol_ppm": 10},
        based_on=["IsoSpecPy", "mzDesk"]),

    module("som", "Site of Metabolism", "B", "Which <b>atoms</b> are metabolized first.",
        [structure(), multiselect("enzymes", "Enzymes", ["CYP3A4", "CYP2D6", "CYP2C9", "CYP2C19", "CYP1A2", "UGT", "All"], ["All"])],
        [table(["Atom", "Rank", "Score", "Reaction type", "Enzyme", "Tool"], caption="SMARTCyp / FAME"), figure("Structure with predicted sites highlighted")], subgroup="dmpk"),

    module("stability", "Metabolic Stability", "B", "<b>t½</b> and <b>intrinsic clearance</b> from an in vitro incubation.",
        [grid("timecourse", "Time course", ["Time (min)", "% remaining"]), select("system", "System", ["Liver microsomes", "Hepatocytes", "S9", "Plasma"]),
         number("protein", "Protein (mg/mL) or cells (10⁶/mL)", 0.5, 0.1), select("species", "Species", ["Human", "Rat", "Mouse", "Dog", "Monkey"])],
        [cards([("t½", "min"), ("CLint", "µL/min/mg or µL/min/10⁶ cells"), ("Scaled CLint", "mL/min/kg"), ("Predicted hepatic CL", "mL/min/kg, well-stirred"), ("Extraction ratio", "")],
               id="summary", style="metric"),
         figure("ln(% remaining) vs time with fit", id="fit_plot"),
         table(["Parameter", "Value", "Unit", "Assumption"], caption="Scaling factors used", id="parameters"),
         notes("Assumptions", id="notes")], subgroup="dmpk", live=True,
        example={"timecourse": "Time (min)\t% remaining\n0\t100\n5\t82\n15\t61\n30\t40\n45\t26",
                 "system": "Liver microsomes", "protein": 0.5, "species": "Human"},
        based_on=["mzDesk"]),

    module("nca", "NCA Pharmacokinetics", "B", "<b>Non-compartmental parameters</b> from a concentration–time profile.",
        [grid("profile", "Concentration–time", ["Time (h)", "Concentration (ng/mL)"]), number("dose", "Dose (mg)", step=0.01),
         select("route", "Route", ["IV bolus", "IV infusion", "Oral"]), number("lambda_points", "Points for λz", 3, 1)],
        [cards([("Cmax", "ng/mL"), ("Tmax", "h"), ("AUC0-t", "ng·h/mL"), ("AUC0-∞", "ng·h/mL"), ("t½", "h"), ("CL/F", "L/h"), ("Vz/F", "L"), ("MRT", "h")],
               id="summary", style="metric"),
         figure("Concentration–time (semi-log) with λz fit", id="profile_plot"), table(["Parameter", "Value", "Unit", "Method"], id="parameters"),
         notes("Assumptions", id="notes")], subgroup="dmpk", live=True,
        example={"profile": "Time (h) Concentration (ng/mL)\n0.5 42\n1 88\n2 131\n4 152\n6 121\n8 84\n12 41\n24 9",
                 "dose": 100, "route": "Oral", "lambda_points": 4},
        based_on=["mzDesk"]),

    module("sil_is", "Isotope-Labelled Internal Standard", "B", "Which <b>labelling</b> to buy: mass shift, stable positions, crosstalk.",
        [structure(label="Analyte"), spectrum(label="Analyte MS2 (optional)"), number("min_shift", "Minimum mass shift (Da)", 3, 1), checkbox("check_vendors", "Check vendor availability", True)],
        [table(["Labelling", "Mass shift", "Positions", "Exchangeable-H risk", "Metabolic-loss risk", "Label kept in quantifier fragment", "Crosstalk with analyte M+n (%)", "Vendor"]),
         figure("Structure with labelled positions highlighted")], subgroup="dmpk"),

    # ---- Inorganic MS（非核心）----
    module("icp_isotopes", "Isotopes & Internal Standard", "B", "Choose an ICP-MS isotope from natural abundance, atomic isobars and matrix polyatomics; match an internal standard by mass and first ionization energy.",
        [select("element", "Element", list(ELEMENT_SYMBOLS), "As"),
         multiselect("matrix", "Plasma / solvent / matrix elements", list(ELEMENT_SYMBOLS), list(MATRIX_ELEMENTS))],
        [cards([("Analyte", ""), ("Recommended isotope", ""), ("Top internal standard", ""), ("Analyte first IE", "eV")], "Recommendation", id="summary", style="metric",
                help=["Ranking is a transparent rule, not a measurement: isotopes with an atomic isobar from another element rank below all isotopes without one. Within each group, score = natural abundance fraction ÷ (1 + sum of polyatomic likelihood scores).",
                      "A likelihood score is the product of natural isotope abundances, not an intensity. Real isotope choice depends on the actual sample matrix, instrument, resolution, cell mode and sensitivity."]),
         table(["Isotope", "Exact mass (Da)", "Natural abundance (%)", "Atomic isobars", "Polyatomics (n)", "Strongest polyatomic"], "Isotopes", id="isotope_table",
               help="Atomic isobars are isotopes of other elements with the same nominal mass and natural abundance > 0. Polyatomic counts and the strongest row come from the Interferences page engine with 2–3 atoms and the selected matrix."),
         table(["Rank", "Isotope", "Rule score", "Reason"], "Recommended isotope(s)", id="recommendations",
               help="Higher rule score is better within the no-atomic-isobar and atomic-isobar groups. The reason lists the natural abundance, total polyatomic likelihood penalty and strongest candidate."),
         table(["Rank", "Internal standard", "|Δ mass| (Da)", "First IE (eV)", "|Δ IE| (eV)", "Combined distance"], "Internal standard suggestions", id="internal_standards",
               help=["Candidates: ⁷Li / ⁶Li, ⁹Be, ⁴⁵Sc, ⁷²Ge / ⁷⁴Ge, ⁸⁹Y, ¹⁰³Rh, ¹¹⁵In, ¹⁵⁹Tb, ¹⁶⁵Ho, ¹⁷⁵Lu, ¹⁸⁵Re, ²⁰⁹Bi; the analyte element itself is excluded.",
                     "Combined distance = 0.5 × (|candidate isotope mass − analyte average mass| / analyte average mass + |candidate IE − analyte IE| / analyte IE); lower is better. First ionization energies are from NIST ASD for all offered elements.",
                     "This ranks only mass and first-ionization-energy closeness. It does not model recovery, ionization suppression, contamination, stability or instrument-specific response."])],
        subgroup="inorganic", live=True, example={"element": "As", "matrix": list(MATRIX_ELEMENTS)}, based_on=["IsoSpecPy", "NIST ASD", "mzDesk"]),

    module("icp_interferences", "Interferences", "B", "Enumerate isotope-specific atomic isobars, M²⁺ ions and 2–3-atom polyatomic ions at an ICP-MS target nominal mass.",
        [text("target", "Target isotope", "75As", "75As or As-75"),
         multiselect("matrix", "Plasma / solvent / matrix elements", list(ELEMENT_SYMBOLS), list(MATRIX_ELEMENTS)),
         select("max_atoms", "Maximum atoms per polyatomic ion", ["2", "3"], "3"),
         multiselect("charge_states", "Charge states considered", ["1+", "2+"], ["1+", "2+"])],
        [cards([("Interferences found", "n"), ("Strongest", "likelihood score"), ("Highest resolving power", "m/|Δm|")], "Summary", id="summary", style="metric"),
         table(["Ion", "Type", "Exact m/z", "Δm (mDa)", "Required R", "Likelihood score"], "Interferences", id="interferences",
               help=["Candidates match the target nominal mass: other-element atomic isobars (1+), isotopes at twice the target nominal mass (2+), and singly charged 2–3-atom combinations of the selected matrix elements. This includes ArX⁺, Ar₂⁺, oxides, hydroxides and hydrides.",
                     "Exact m/z subtracts one electron mass per charge. Δm = candidate − target; required R = target m/z ÷ |Δm|. An exact mass coincidence has infinite required R.",
                     "Likelihood score = product of constituent natural isotope abundances. It is not an intensity or ion-formation probability. Rows are score-sorted and capped at 50. No collision/reaction-cell recommendation is made."]),
         notes("Scope & assumptions", id="notes", collapsed=True)],
        subgroup="inorganic", live=True,
        example={"target": "75As", "matrix": list(MATRIX_ELEMENTS), "max_atoms": "3", "charge_states": ["1+", "2+"]}, based_on=["IsoSpecPy", "mzDesk"]),

    module("icp_idms", "Isotope Dilution (IDMS)", "B", "Single-spike, mole-based isotope dilution with composition-specific molar masses, ratio sensitivity and an optimal-blend hint.",
        [select("element", "Element", list(ELEMENT_SYMBOLS), "Ag"),
         text("isotope_a", "Reference isotope a", "107Ag"), text("isotope_b", "Spike isotope b", "109Ag"),
         number("a_x", "Sample A_x (blank = natural)", step="any"), number("b_x", "Sample B_x (blank = natural)", step="any"),
         number("a_y", "Spike A_y", 0.01, "any"), number("b_y", "Spike B_y", 0.99, "any"),
         number("spike_concentration", "Spike concentration c_y", 10.0), select("spike_unit", "Spike concentration unit", ["µg/g", "µmol/g"]),
         number("m_x", "Sample mass m_x (g)", 1.0), number("m_y", "Spike mass m_y (g)", 0.5),
         number("r_b", "Measured blend ratio R_B = a/b (mass-bias corrected)", 1.0, "any")],
        [cards([("Sample concentration", "µg/g"), ("Sample amount content", "µmol/g"),
                ("Error magnification", "|∂ln c_x / ∂ln R_B|"), ("Optimal blend R_opt", "√(R_x·R_y)")], "IDMS result", id="summary", style="metric",
                help=["c_x = c_y·(m_y/m_x)·(A_y − R_B·B_y)/(R_B·B_x − A_x). c_x and c_y are amount contents (µmol/g); m_x and m_y must use the same mass unit.",
                      "µg/g spike input is converted with the spike molar mass from its isotopic composition; sample µg/g output uses the sample composition. Unlisted isotopes retain natural relative proportions within the remaining atom fraction.",
                      "Error magnification = |∂c_x/∂R_B · R_B/c_x|. R_opt = √(R_x·R_y), where R_x = A_x/B_x and R_y = A_y/B_y. R_B must already be mass-bias corrected; this is not a full uncertainty budget."]),
         table(["Quantity", "Value", "Unit / meaning"], "Intermediate values", id="intermediates"),
         notes("Equation & limits", id="notes", collapsed=True)],
        subgroup="inorganic", live=True,
        example={"element": "Ag", "isotope_a": "107Ag", "isotope_b": "109Ag", "a_x": "", "b_x": "", "a_y": "0.01", "b_y": "0.99",
                 "spike_concentration": "10", "spike_unit": "µg/g", "m_x": "1", "m_y": "0.5", "r_b": "1"},
        based_on=["IsoSpecPy", "mzDesk"]),

    module("icp_massbias", "Mass Bias & δ Values", "B", "Correct isotope ratios with three mass-bias laws or standard–sample bracketing, and convert between ratios and δ values.",
        [heading("Mass-bias factor from a standard"),
         text("mb_a", "Numerator isotope a", "107Ag"), text("mb_b", "Denominator isotope b", "109Ag"),
         number("mb_measured_std", "Measured standard ratio R_meas,std", 0.99, "any"),
         number("mb_true_std", "Certified standard ratio R_true,std", 1.0, "any"),
         number("mb_measured_sample", "Measured sample ratio", 0.995, "any"),
         heading("Standard–sample bracketing"),
         number("bracket_sample", "Measured sample ratio", 0.995, "any"),
         number("bracket_before", "Measured standard before", 0.990, "any"),
         number("bracket_after", "Measured standard after", 0.992, "any"),
         number("bracket_true", "Certified standard ratio", 1.0, "any"),
         heading("δ value and inverse"),
         number("delta_sample", "Sample ratio R_sample", 1.002, "any"),
         number("delta_standard", "Standard ratio R_standard", 1.0, "any"),
         number("delta_per_mil", "δ value for inverse (‰)", 2.0, "any")],
        [head("Mass-bias correction"),
         cards([("K = R_true/R_meas", ""), ("Linear ε", "per Da"), ("Power ε", "per Da"), ("Russell β", "")], "Mass-bias factors", id="mass_bias", style="metric",
               help=["K = R_true,std/R_meas,std. Linear: R_corr = R_meas·[1 + ε(m_a − m_b)]. Power: R_corr = R_meas·(1 + ε)^(m_a − m_b).",
                     "Exponential (Russell): R_corr = R_meas·(m_a/m_b)^β, with β = ln(K)/ln(m_a/m_b). Isotope masses come from IsoSpecPy by isotope label.",
                     "All laws require positive ratios and different isotope masses. They are alternative empirical mass-bias models, not interchangeable physical measurements; use the law validated for the instrument and element."]),
         table(["Law", "Parameter", "Value", "Correction factor", "Corrected sample ratio"], "Corrected sample ratios", id="corrected_ratios"),
         head("Standard–sample bracketing"),
         cards([("Corrected sample ratio", ""), ("Mean measured standard", "")], "Bracketing result", id="bracketing", style="metric",
               help="Corrected R_sample = R_meas,sample·R_true,std/mean(R_meas,std before, after). Bracketing assumes the mean of the two standards represents drift across the sample; it does not estimate uncertainty or nonlinear drift."),
         head("δ value and inverse"),
         cards([("δ from ratios", "‰"), ("R_sample/R_standard from δ", ""), ("R_sample from δ", "")], "δ result", id="delta", style="metric",
               help="δ(‰) = (R_sample/R_standard − 1) × 1000. Inverse: R_sample/R_standard = 1 + δ/1000; multiply by R_standard for R_sample. The inverse requires δ > −1000‰ for a positive ratio.")],
        subgroup="inorganic", live=True,
        example={"mb_a": "107Ag", "mb_b": "109Ag", "mb_measured_std": "0.99", "mb_true_std": "1", "mb_measured_sample": "0.995",
                 "bracket_sample": "0.995", "bracket_before": "0.990", "bracket_after": "0.992", "bracket_true": "1",
                 "delta_sample": "1.002", "delta_standard": "1", "delta_per_mil": "2"},
        based_on=["IsoSpecPy", "mzDesk"]),

    module("icp_sp", "Single-Particle ICP-MS", "B", "Convert detected particle-event signals to spherical diameter and number concentration using dissolved-standard calibration and transport efficiency.",
        [select("material", "Particle material", list(PARTICLE_PRESETS) + ["Custom"], "Au"),
         text("custom_formula", "Custom composition formula (Custom only)", "", "TiO2"), number("density", "Custom density (g/cm³; Custom only)", step=0.01),
         select("element", "Measured element (must occur in the formula)", list(ELEMENT_SYMBOLS), "Au"),
         number("sensitivity", "Dissolved calibration sensitivity S", 1000.0),
         select("sensitivity_unit", "Sensitivity unit", [SENSITIVITY_DWELL, SENSITIVITY_RATE]),
         number("dwell_ms", "Dwell time t_dw (ms)", 1.0), number("flow", "Sample flow q (mL/min)", 0.3),
         number("transport_pct", "Transport efficiency η (%)", 5.0), number("background", "Background counts per dwell I_bg", 2.0),
         textarea("intensities", "Particle event intensities (counts; one per line or a single value)", "5054.728\n6200\n4500", "5054.7\n6200\n4500", 6),
         number("n_events", "Number of events N (blank = listed events)", 3, 1), number("acquisition_s", "Acquisition time t_acq (s)", 60.0),
         number("threshold", "Threshold intensity for size LOD (counts; optional)", 3000.0)],
        [cards([("Mean diameter", "nm"), ("Median diameter", "nm"), ("Number concentration", "particles/mL"), ("Size LOD", "nm")], "Single-particle results", id="summary", style="metric",
                help=["Pace et al., Anal. Chem. 2011, 83, 9361. m_el[µg] = (I − I_bg)·η·q[L/s]·t_dw[s]/S, where S is counts per dwell per µg/L; counts/s sensitivity is first multiplied by t_dw in seconds. m_p = m_el/X and d = (6m_p/(πρ))^(1/3).",
                      "C_N = N/(η·q·t_acq), with q in mL/s for particles/mL. X is the measured-element mass fraction from the formula and IsoSpecPy average atomic masses.",
                      "Only existing event signals are converted: no raw time-series processing, threshold detection, dissolved-background fitting, coating/dissolution or spectral-interference correction. Size LOD is the optional threshold converted to diameter, not a statistical detection limit."]),
         table(["Event", "Intensity (counts)", "Net counts", "Element mass (fg)", "Particle mass (fg)", "Diameter (nm)"], "Particle events", id="events",
               help="Element and particle masses are shown in fg (1 fg = 10⁻¹⁵ g). The table is capped at 200 events; means, medians and the histogram use every supplied event."),
         figure("Diameter histogram", id="histogram", help="One histogram of all supplied particle-event diameters; summary statistics are not based on the capped table."),
         notes("Method & limits", id="notes", collapsed=True)],
        subgroup="inorganic", live=True,
        example={"material": "Au", "custom_formula": "", "density": "", "element": "Au", "sensitivity": "1000",
                 "sensitivity_unit": SENSITIVITY_DWELL, "dwell_ms": "1", "flow": "0.3", "transport_pct": "5", "background": "2",
                 "intensities": "5054.728\n6200\n4500", "n_events": "3", "acquisition_s": "60", "threshold": "3000"},
        based_on=["IsoSpecPy", "mzDesk"]),
]


# ---- Form metadata (Spectrum pages) ----
# Keys are optional and only added when used, so the field serialisation of other pages is unchanged.
#   section / form_order : the frontend groups fields into sections (in this order)
#   required             : a star, aria-required and a pre-flight check before POST
#   required_when        : conditional required ({"field": ..., "equals"/"in"/"empty": ...})
#   header_fallback      : a blank m/z/number is accepted when the spectrum header supplies it
#   visible_when         : hide controls that the chosen mode ignores
#   auto_label / ion_modes / allow_formula / help : display and input capability metadata
def _configure_forms(module_id, sections, *, required=(), required_when=None, header_fallback=(),
                     visible=None, auto_labels=None, ion_modes=None, allow_formula=(), help=None):
    module = next(m for m in MODULES if m["id"] == module_id)
    fields = {f.get("id"): f for f in module["inputs"] if f.get("id")}
    for order, (title, ids) in enumerate(sections):
        for fid in ids:
            if fid in fields:
                fields[fid]["section"] = title
                fields[fid]["form_order"] = order
    for fid in required:
        if fid in fields:
            fields[fid]["required"] = True
    for fid, condition in (required_when or {}).items():
        if fid in fields:
            fields[fid]["required_when"] = condition
    for fid in header_fallback:
        if fid in fields:
            fields[fid]["header_fallback"] = True
    for fid, condition in (visible or {}).items():
        if fid in fields:
            fields[fid]["visible_when"] = condition
    for fid, label in (auto_labels or {}).items():
        if fid in fields:
            fields[fid]["auto_label"] = label
    for fid, modes in (ion_modes or {}).items():
        if fid in fields:
            fields[fid]["ion_modes"] = modes
    for fid in allow_formula:
        if fid in fields:
            fields[fid]["allow_formula"] = True
    for fid, text in (help or {}).items():
        if fid in fields:
            fields[fid]["help"] = text


_configure_forms("ms2_clean",
    [("Spectrum", ["spectrum", "precursor_mz", "ion_mode"]),
     ("Cleaning", ["noise_pct", "normalize", "remove_precursor", "precursor_tol", "sort_by"])],
    required=("spectrum",), header_fallback=("precursor_mz",),
    help={"precursor_mz": "Blank = read from the spectrum header."})
_configure_forms("ms2_similarity",
    [("Query", ["query", "query_precursor"]), ("Reference", ["reference", "reference_precursor"]),
     ("Matching & scores", ["tol_da", "threshold", "ion_mode", "deep"])],
    required=("query", "reference"), header_fallback=("query_precursor", "reference_precursor"),
    visible={"ion_mode": {"field": "deep", "equals": True}},
    help={"query": "A header precursor m/z is used when that field is blank."})
_configure_forms("ms2_library",
    [("Query", ["spectrum_type", "ms2", "precursor", "tol_ppm", "tol_da"]),
     ("Search", ["libraries", "metrics", "threshold"])],
    required=("ms2", "libraries", "metrics"), required_when={"precursor": {"field": "spectrum_type", "equals": "MS/MS"}},
    header_fallback=("precursor",),
    visible={"precursor": {"field": "spectrum_type", "equals": "MS/MS"},
             "tol_ppm": {"field": "spectrum_type", "equals": "MS/MS"},
             "tol_da": {"field": "spectrum_type", "equals": "MS/MS"}},
    auto_labels={"precursor": "Auto (any adduct)"},
    help={"precursor": "Required unless the spectrum header has it. Auto = any adduct.",
          "metrics": "Hits are ranked by the first selected metric in list order."})
_configure_forms("ms1_formula",
    [("Query", ["precursor", "ms1", "ms2"]), ("Formula search", ["instrument", "formula", "halogen"])],
    required=("precursor",), header_fallback=("precursor",),
    help={"precursor": "Required unless the MS2 header has it.",
          "formula": "Skips ranking; the MS2 is explained with this neutral formula.",
          "halogen": "Ignored when a formula is given."})
_configure_forms("ms2_structure",
    [("Query", ["precursor", "ms2", "ms1"]), ("Search scope", ["halogen", "silicon", "extended"])],
    required=("precursor", "ms2"), header_fallback=("precursor",),
    auto_labels={"precursor": "Auto ([M+H]+ / [M-H]-)"},
    help={"extended": "Allows more F per carbon (PFAS) and other unusual compositions."})
_configure_forms("ms2_denovo",
    [("Query", ["precursor", "ms2", "ms1"]), ("Generation", ["formula", "mode", "halogen"])],
    required=("precursor", "ms2"), header_fallback=("precursor",),
    ion_modes={"precursor": ["Positive"]},
    auto_labels={"precursor": "Auto ([M+H]+)"},
    help={"ms1": "Used only when Formula is blank.",
          "formula": "Blank = BUDDY's top formula. Enter one for [M]+ or [M+H-2H2O]+.",
          "halogen": "Used only when Formula is blank."})
_configure_forms("ms2_ei",
    [("Structure", ["structure"])],
    required=("structure",))
_configure_forms("ms2_esi",
    [("Structure", ["structure"]), ("Ion", ["adduct"])],
    required=("structure",),
    help={"adduct": "Lipids: seven adducts. Other compounds: [M+H]+ and [M-H]- only."})
_configure_forms("ms2_fragments",
    [("Candidate", ["structure", "precursor"]), ("Spectrum", ["ms2"]), ("Fragmentation", ["tree_depth", "tol_abs", "tol_ppm"])],
    required=("structure", "ms2"), header_fallback=("precursor",),
    auto_labels={"precursor": "Auto ([M+H]+ / [M-H]-)"},
    help={"precursor": "Blank = theoretical m/z."})
_configure_forms("exclusivity",
    [("Input", ["structure", "precursor"]), ("Search", ["tol_ppm", "databases"]), ("Isomer check", ["ms2", "tol_da"])],
    allow_formula=("structure",),
    help={"structure": "Structure or formula. With only a precursor m/z, only isobars are searched.",
          "precursor": "Blank = theoretical m/z of the compound.",
          "databases": "Used by the same-formula and isobaric tables.",
          "ms2": "Needed for the isomer check.",
          "tol_da": "Isomer check only."})

# Also expose the formula capability on the pages whose runners pass allow_formula=True.
for _mid, _fid in (("ms1_theory", "structure"), ("met_shift", "structure")):
    _m = next((m for m in MODULES if m["id"] == _mid), None)
    if _m:
        _f = next((f for f in _m["inputs"] if f.get("id") == _fid), None)
        if _f:
            _f["allow_formula"] = True


# ---- Pipelines ----
def stage(n, title, expected_s=1.0):
    """一个步骤：后端（mzdesk/annotation/pipeline.py）按步骤 id 报告状态、用时、说明和输出块。modules 留空。
    expected_s：该步骤的典型用时（秒），只给进度条用——每步在进度条上占的份额与它成正比，运行中的步骤按它推算填充速度；不是时限。"""
    return {"id": n, "title": title, "modules": [], "expected_s": float(expected_s)}


def _section(title, *fields, collapsed=False):
    """Put the fields under one form heading (the pipeline page draws one group per heading).
    collapsed: the section starts closed, behind its heading (a <details>, like "Advanced"); its fields carry section_collapsed."""
    return [{**f, "section": title, **({"section_collapsed": True} if collapsed else {})} for f in fields]


# adduct_rules.ADDITIVES (tests/test_pipeline.py checks they agree; spec does not import the rule engine)
PIPELINE_ADDITIVES = ["Formic acid", "Ammonium formate", "Ammonium acetate", "Acetic acid", "None"]
# library_step.METRIC_CHOICES (tests/test_pipeline.py checks they agree)
PIPELINE_LIBRARY_METRICS = [METRIC_NAMES.ENTROPY_WEIGHTED, METRIC_NAMES.ENTROPY_UNWEIGHTED, METRIC_NAMES.COSINE_GREEDY, METRIC_NAMES.COSINE_HUNGARIAN,
                            METRIC_NAMES.MOD_COSINE_GREEDY, METRIC_NAMES.MOD_COSINE_HUNGARIAN, METRIC_NAMES.NEUTRAL_LOSSES, METRIC_NAMES.NIST_DOT,
                            METRIC_NAMES.REVERSE_DOT, METRIC_NAMES.BINNED_COSINE]
_PIPELINE_ADVANCED = [
    {**select("additive", "Mobile phase additive", ["Unknown"] + PIPELINE_ADDITIVES), "help": "Unknown = the best additive for each adduct."},
    {**select("alkali", "Sodium / potassium present", ["Unknown", "Present", "Absent"]), "help": "Unknown counts as present."},
    {**number("ph", "Mobile phase pH", 2.7, 0.1), "help": "Default 2.7 (0.1 % formic acid)."},
    {**multiselect("lib_metrics", "Library similarity algorithms", PIPELINE_LIBRARY_METRICS, [METRIC_NAMES.ENTROPY_WEIGHTED]), "checkboxes": True, "required": True,
     "help": "A hit passes when any ticked algorithm reaches the threshold; the first ticked one orders the hits."},
    {**number("lib_threshold", "Similarity threshold", 0.80, 0.01), "min": 0, "max": 1,
     "help": "Minimum similarity for a Level-2 library match."},
    {**number("lib_min_peaks", "Minimum matched peaks", 4, 1), "min": 1,
     "help": "A library hit needs at least this many matched peaks (all of them if the spectrum has fewer)."},
    # ticked once FRIGID is installed (local version), never on the hosted server (MZDESK_ANNOTATION_FRIGID); spec() reads it again on every call.
    # needs: ticking it with FRIGID missing offers the download (POST /api/local/needs, mzdesk/local/targets.py)
    {**checkbox("run_frigid", "Run FRIGID", value=frigid_workflow_default()), "needs": ["frigid"],
     "help": "De novo structures, when the library and MS-FINDER do not settle it. Positive mode only; slow."}]
# MS/MS Annotation, optional Ion Mobility fields (the end of the Advanced section): one observed CCS checks every structure candidate in step 6 (mzdesk/annotation/ccs_step.py)
_PIPELINE_CCS = [
    {**number("ccs_observed", "Observed CCS (Å²)", None, None, "e.g. 152.3"), "min": 0, "new_row": True,
     "help": "One N₂ value. For TIMS, use the CCS your software reports, not 1/K0."},
    {**select("ccs_platform", "Platform", [""] + CCS_PLATFORMS, value=""), "placeholder": "Select the platform",
     "help": "Required with a CCS."},
    select("ccs_gas", "Drift gas", CCS_GASES),
    {**select("ccs_tw_calibration", "TWIMS calibration", CCS_TW_CALIBRATIONS), "visible_when": {"field": "ccs_platform", "equals": "TWIMS"}}]
_PIPELINE_SPECTRUM = _section("Spectrum", {**precursor(), "required": True, "header_fallback": True, "auto_label": "Auto (rank the adducts)"},
                              {**spectrum(label="MS2 spectrum"), "required": True}, ms1())

# EI MS Annotation / Targeted Method (mzdesk/workflows/). The option lists are read by the workflows' parse_inputs; spec stays free of heavy imports.
TARGETED_INSTRUMENTS = ["Any", "LC-MS", "GC-MS"]
EI_DERIVATIZATIONS = ["Auto (from marker ions)", "None", "TMS", "MeOx + TMS", "TBDMS"]
# Real EI spectra for the EI MS Annotation examples, from the local MassBank index. Licences checked in the MassBank-data records: PR010014 CC BY-SA
# (Kusano M, Fukushima A, RIKEN; Glycine, 3 TMS, GC-EI-TOF), NL0043 CC BY (MILK, NILU; Q Exactive GC, high-resolution EI). Both records are in the local
# library, so each example matches itself; the replicate spectra of other labs show up next to it.
EI_GLYCINE_3TMS = "60 17\n61 11\n66 10\n67 1\n69 1\n70 11\n71 9\n72 25\n73 999\n74 86\n75 50\n76 3\n77 10\n83 1\n84 8\n85 7\n86 399\n87 41\n88 16\n89 2\n95 2\n98 1\n99 4\n100 221\n101 44\n102 30\n103 14\n104 2\n105 4\n113 6\n114 6\n115 9\n116 13\n117 38\n118 7\n119 12\n120 1\n121 1\n128 1\n129 3\n130 39\n131 42\n132 12\n133 116\n134 17\n135 9\n136 1\n142 1\n143 1\n144 10\n145 3\n146 5\n147 319\n148 49\n149 24\n150 2\n158 12\n159 4\n160 7\n161 2\n172 7\n173 2\n174 750\n175 148\n176 73\n177 15\n178 2\n188 7\n189 2\n190 2\n202 3\n203 1\n204 4\n205 1\n246 11\n247 4\n248 146\n249 42\n250 20\n251 4\n276 49\n277 13\n278 6\n279 1"
EI_DIPENTYL_PHTHALATE = "104.0257 8\n105.03365 10\n121.02855 21\n122.03618 5\n123.04401 5\n149.02344 999\n150.02673 88\n167.03389 3\n176.04677 3\n219.10176 12\n237.11231 24\n238.11569 3"


def _caffeine_example(libraries, ms2):
    return {"precursor": _prec("195.0877"), "ms2": _paste(ms2), "ms1": CAFFEINE_MS1, "libraries": libraries, "tol_ppm": 10.0, "tol_da": 0.02,
            "omics": "Either / unknown", "halogen": False, "local_isomers": False, "additive": "Unknown", "alkali": "Unknown", "ph": 2.7,
            "lib_metrics": [METRIC_NAMES.ENTROPY_WEIGHTED], "lib_threshold": 0.80, "lib_min_peaks": 4, "run_frigid": frigid_workflow_default()}


PIPELINES = [
    {"id": "esi_ms2_annotation", "title": "MS/MS Annotation",
     "one_liner": "ESI MS/MS spectrum → adduct, formula and ranked structure candidates",
     "based_on": ["BUDDY", "MS-FINDER", "FRIGID", "MassBank/GNPS"],
     "inputs": (_PIPELINE_SPECTRUM
                + _section("Search", multiselect("libraries", "Libraries", LIBRARIES),
                           number("tol_ppm", "Precursor tolerance (ppm)", 10.0, 0.1), number("tol_da", "Fragment tolerance (Da)", 0.02, 0.001),
                           select("omics", "Omics", ["Either / unknown", "Metabolomics", "Lipidomics"]))
                + _section("Options", checkbox("halogen", "Allow Cl/Br in formulas"),
                           checkbox("local_isomers", "Compare with database isomers"))
                + _section("Advanced", *_PIPELINE_ADVANCED, *_PIPELINE_CCS)),       # closed by default; the observed CCS (Ion Mobility) is its last four fields
     # expected_s: Adduct includes the formula search (step 2 only reads it back): CASMI 2022 500 challenges median 5.0 s, caffeine cold 9.6-10.3 s;
     # Merge and re-rank: 0.03-0.2 s on the examples; an observed CCS adds about 0.01 s per structure, plus 1.0 s once per server process for loading the CCS prediction model
     # (not worth a larger share of the bar). Library search: first query of a server loads the index (7.5 s cold, 0.03 s warm); MS-FINDER / FRIGID run only when the library does not decide
     # (CASMI medians 44 s / 95 s when run, caffeine 6-11 s): conservative middles; Isomer check re-runs MS-FINDER (CASMI median 7.8 s, caffeine 4.5-9 s)
     "stages": [stage(1, "Adduct", 6.0), stage(2, "Formula", 0.5), stage(3, "Library search", 4.0), stage(4, "MS-FINDER structures", 20.0),
                stage(5, "FRIGID structures", 40.0), stage(6, "Merge and re-rank", 0.5), stage(7, "Isomer check", 8.0), stage(8, "Answer", 0.1)],
     "examples": [{"label": "Caffeine, library match", "values": _caffeine_example(LIBRARIES, CAFFEINE_MASSBANK_MS2)},
                  {"label": "Caffeine, no library (MS-FINDER)", "values": _caffeine_example([], CAFFEINE_MSFINDER_MS2)}]},
    {"id": "ei_ms_annotation", "title": "EI MS Annotation",
     "one_liner": "EI spectrum (GC-MS) → candidate metabolites with evidence",
     "based_on": ["NEIMS", "MassBank", "RDKit", "mzDesk"],
     "inputs": (_section("Spectrum",
                         {**spectrum("spectrum", "EI spectrum"), "required": True},
                         {**select("derivatization", "Derivatization", EI_DERIVATIZATIONS), "help": "Auto reads the marker ions; MeOx and TBDMS need the molecular ion."},
                         {**number("molecular_ion", "Molecular ion m/z", None, 0.0001), "help": "M+• of the derivative, if derivatized."},
                         {**textarea("candidates", "Candidate structures", placeholder="One SMILES or name per line", rows=4), "network_names": True})
                + _section("Retention index",
                           # no required_when: Fields.tsx marks any field that has one with a permanent star, and both of these are optional; the pairing is
                           # enforced by ei_ms_annotation.parse_inputs (HTTP 400)
                           grid("alkanes", "Alkane series", ["Carbon number", "RT"]),
                           {**number("sample_rt", "Sample RT", None, 0.01), "help": "Same unit as the alkane RTs. Needs the alkane series."})),
     # expected_s: both examples run in-process: EI library search is 6.2 s on the first query of a process (index load / JIT; up to ~15 s on a cold disk
     # cache) and 0.25-0.3 s afterwards, every other step is 0.01-0.13 s (NEIMS 0.13 s cold); stats DB: 3 runs, median 7.1 s in total
     "stages": [stage(1, "Spectrum cleanup", 0.1), stage(2, "EI library search", 7.0), stage(3, "Derivatization markers", 0.05), stage(4, "Candidates", 0.2),
                stage(5, "Derivatives", 0.2), stage(6, "NEIMS prediction", 0.5), stage(7, "Similarity and rerank", 0.1), stage(8, "Kovats RI", 0.05), stage(9, "Answer", 0.1)],
     "examples": [{"label": "Glycine, 3 TMS (MassBank MSBNK-RIKEN-PR010014)",
                   "values": {"spectrum": _paste(EI_GLYCINE_3TMS), "derivatization": EI_DERIVATIZATIONS[0], "molecular_ion": "291"}},
                  {"label": "Dipentyl phthalate, high resolution (MassBank MSBNK-NILU-NL0043)",
                   "values": {"spectrum": _paste(EI_DIPENTYL_PHTHALATE), "derivatization": EI_DERIVATIZATIONS[0]}}]},
    {"id": "targeted_method", "title": "Targeted Method",
     "one_liner": "Compound → method sheet: platform, column, pH, polarity, main ion, methods",
     "based_on": ["MolGpKa", "PubChem", "RepoRT", "mzDesk"],
     # step 1 always asks PubChem / UniChem, also for a SMILES (workflows/targeted_method.py _lookup): the page sentence makes the "Internet search: off" dialog ask
     "network": "Targeted Method looks the compound up in PubChem and UniChem by its InChIKey or SMILES; when that fails, the sheet uses values computed locally.",
     "inputs": ([{**structure(label="Compound"), "required": True}]
                + _section("Conditions",
                           {**number("ph", "pH", None, 0.1, "2.7"), "min": 0, "max": 14, "help": "Blank = 2.7 (0.1 % formic acid)."},
                           select("additive", "Mobile-phase additive", PIPELINE_ADDITIVES),
                           select("instrument", "Instrument", TARGETED_INSTRUMENTS))),
     # expected_s: both examples run in-process (cold / warm): 0.0-0.4 s per step except Adducts (2.4 s cold) and an uncached name lookup in Identity (1-2 s);
     # stats DB: 5 runs, median 4.7 s in total
     "stages": [stage(1, "Identity", 1.5), stage(2, "Platform", 0.2), stage(3, "Separation", 0.4), stage(4, "Ionization", 0.2), stage(5, "Adducts", 2.0),
                stage(6, "Non-MS detection", 0.1), stage(7, "Derivatization", 0.1), stage(8, "Known methods", 0.3), stage(9, "Theoretical m/z", 0.1)],
     "examples": [{"label": "Caffeine, any platform", "values": {"structure": _name("caffeine"), "additive": "Formic acid", "instrument": "Any"}},
                  {"label": "D-Glucose, GC-MS", "values": {"structure": _name("D-glucose"), "additive": "Formic acid", "instrument": "GC-MS"}}]},
    {"id": "lipid_fragments", "title": "Lipid Fragment Explanation",
     "one_liner": "Lipid name or structure + MS2 spectrum → names, illustrative structure, explained peaks",
     "based_on": ["LIPID+", "MetFrag", "RDKit"],
     # modes: the structure field offers a lipid shorthand first (no Name / ID lookup, so nothing goes online); a SMILES or a drawing gives no names
     "inputs": (_section("Lipid",
                         {**structure("lipid", "Lipid"), "modes": [LIPID_NAME_MODE, "Draw", "SMILES / InChI"], "network_modes": [], "required": True,
                          "placeholder": "PC 16:0/18:1",
                          "help": "A short name gives the names and an illustrative structure. SMILES or a drawing gives no names."})
                + _section("Spectrum",
                           {**precursor(label="Precursor m/z (optional)", adducts=METFRAG_ADDUCTS), "header_fallback": True,
                            "auto_label": "Auto (match the precursor m/z)",
                            "help": "Blank = theoretical m/z. Auto picks the adduct that fits the precursor m/z; without one, [M+H]+ / [M-H]-."},
                           {**spectrum(), "required": True})
                + _section("Fragmentation",
                           {**slider("tree_depth", "Tree depth (bond cleavages)", 1, 4, 1, 3), "help": "Default 3. Higher is slower and can over-explain."},
                           number("tol_abs", "Absolute tolerance (Da)", 0.002, 0.001), number("tol_ppm", "Relative tolerance (ppm)", 10.0, 1))),
     # expected_s: names and the drawing are local (0.05-0.3 s); one MetFrag run at depth 3 on a lipid is 1.5-3 s (JVM start included)
     "stages": [stage(1, "Lipid name", 0.1), stage(2, "Structure", 0.3), stage(3, "Fragment explanation", 2.5)],
     "examples": [{"label": "PG 18:1/18:1, negative",
                   "values": {"lipid": {"mode": LIPID_NAME_MODE, "value": "PG 18:1/18:1"}, "precursor": _prec("773.5359", "Negative", "Auto"),
                              "ms2": _paste(PG_36_2_MS2), "tree_depth": 3}},
                  {"label": "PG 36:2, chains assumed",
                   "values": {"lipid": {"mode": LIPID_NAME_MODE, "value": "PG 36:2"}, "precursor": _prec("773.5359", "Negative", "Auto"),
                              "ms2": _paste(PG_36_2_MS2), "tree_depth": 3}}]},
]


GROUPS = [{"id": "modules", "title": "Modules"}, {"id": "pipelines", "title": "Workflows"}]
# Modules 下的小类，侧栏按这里的顺序列出（组内已接引擎的模块排前面）；Workflows 排在 Modules 之后
SUBGROUPS = [{"id": "molecular", "title": "Compound"}, {"id": "detectability", "title": "Detectability"},
             {"id": "ms1", "title": "MS1"}, {"id": "ms2", "title": "Spectrum"}, {"id": "orthogonal", "title": "Orthogonal"}, {"id": "dmpk", "title": "DMPK"},
             {"id": "inorganic", "title": "ICP-MS"}, {"id": "gnps", "title": "Dataset Search"}]


# Built but not offered in this release (left for 2.0): DMPK, ICP-MS and MRM Transitions. They stay defined here (with their engines and tests) so
# they can be switched on again by emptying these two sets; the page list (spec()) and the run endpoint (api/app.py) skip them.
DEFERRED_SUBGROUPS = frozenset({"dmpk", "inorganic"})
DEFERRED_MODULES = frozenset({"mrm"})


def is_deferred(module):
    return module["subgroup"] in DEFERRED_SUBGROUPS or module["id"] in DEFERRED_MODULES


def _pipelines_now():
    """PIPELINES with the run_frigid default read now: FRIGID can be downloaded while the server runs, and the checkbox (and the examples that
    set it) follow. Files only, no subprocess. Everything else is shared with PIPELINES, not copied."""
    ticked = frigid_workflow_default()
    out = []
    for pipe in PIPELINES:
        if any(f["id"] == "run_frigid" for f in pipe["inputs"]):
            pipe = {**pipe, "inputs": [{**f, "value": ticked} if f["id"] == "run_frigid" else f for f in pipe["inputs"]],
                    "examples": [{**e, "values": {**e["values"], "run_frigid": ticked}} if "run_frigid" in e["values"] else e for e in pipe.get("examples", [])]}
        out.append(pipe)
    return out


def spec():
    return {"groups": GROUPS, "pipelines": _pipelines_now(),
            "subgroups": [g for g in SUBGROUPS if g["id"] not in DEFERRED_SUBGROUPS],
            "modules": [m for m in MODULES if not is_deferred(m)]}
