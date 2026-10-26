"""页面定义（唯一来源）。前端按这里的 pipelines / modules 渲染输入和输出格式。

字段类型：structure, precursor, spectrum, ms1, number, text, textarea, select, multiselect, radio,
checkbox, slider, file, table, heading。
输出块类型：table, cards, figure, code, formats, notes, button, heading, portal, chips, reaction, grouptable。
  任何块都可带 help（"?" 弹窗，字符串或字符串列表；运行时数据里的 data.help 会并进去）。数据形状见各 helper 的 docstring，
  运行时数据可用 mzdesk/api/blocks.py 里的 builder 生成。
  formats：同一份内容的多种格式，数据 {"formats": {格式名: 文本}}，前端画成一排标签 + 一个代码框。
  id：后端结果按 id 填进这个块；没有 id 的块只画格式。
  key：pipeline 的 stage 用 exclude 按 key 过滤块。
模块：group 是侧栏大类（modules / pipelines）；subgroup 是 Modules 下的小类（molecular / detectability / ms1 / ms2 / orthogonal / dmpk / inorganic / tools；ms2 的显示名是 Spectrum），
  None = 直接列在 Modules 下。live=True 表示已接引擎（后端 mzdesk/api/runners.py 里有对应 runner）。
  example：一组可直接运行的示例输入，前端 "Load example" 按钮填入。
"""
from mzdesk.shared.ms2_clean import NORMALIZE, SORT_BY
# 注意：resolve.py 导入时不联网，所以这里可以直接导入 KINDS；不要在 spec 里导入任何联网代码
from mzdesk.identity.resolve import KINDS as IDENTITY_KINDS
# cfmid_engine 顶层只用标准库，导入它不会拉进 rdkit / Docker
from mzdesk.ms2.cfmid_engine import ADDUCTS as ESI_ADDUCTS
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
LIBRARIES = ["MassBank", "GNPS"]            # 本地索引见 mzdesk/ms2/library_setup.py
LIBRARY_METRICS = METRIC_NAMES.LIBRARY_ORDER    # Library search 页的指标多选 = Similarity 页得分表里的全部 16 项（同一份名字，谱熵排第一）
LIBRARY_DEFAULT_METRICS = [METRIC_NAMES.ENTROPY_WEIGHTED, METRIC_NAMES.COSINE_GREEDY]
DATABASES = ["HMDB", "PubChem", "LIPID MAPS", "KEGG", "ChEBI", "DrugBank", "NORMAN SusDat"]


# ---- 字段 ----
def structure(id="structure", label="Compound", optional=False):
    return {"type": "structure", "id": id, "label": label, "optional": optional}
def precursor(id="precursor", adduct=True, label="Precursor m/z", adducts=None, mz=True):
    """mz=False: 不显示 m/z 输入框（只选离子模式和加合物）；该字段的加合物列表里没有 Auto 时，默认取列表第一项。"""
    f = {"type": "precursor", "id": id, "adduct": adduct, "label": label,
         "adducts": adducts or {"Positive": ADDUCTS_POS, "Negative": ADDUCTS_NEG}}
    if not mz:
        f["mz"] = False
    return f
def spectrum(id="ms2", label="MS2 spectrum", multi=False):
    return {"type": "spectrum", "id": id, "label": label, "multi": multi}
def ms1(id="ms1"):
    return {"type": "ms1", "id": id, "label": "MS1 isotope pattern (m/z intensity per line, optional)"}
def number(id, label, value=None, step=None, placeholder=None):
    return {"type": "number", "id": id, "label": label, "value": value, "step": step, "placeholder": placeholder}
def text(id, label, value="", placeholder=None):
    return {"type": "text", "id": id, "label": label, "value": value, "placeholder": placeholder}
def textarea(id, label, value="", placeholder=None, rows=5):
    return {"type": "textarea", "id": id, "label": label, "value": value, "placeholder": placeholder, "rows": rows}
def select(id, label, options, value=None):
    return {"type": "select", "id": id, "label": label, "options": options, "value": value or options[0]}
def multiselect(id, label, options, value=None):
    return {"type": "multiselect", "id": id, "label": label, "options": options, "value": value if value is not None else options}
def radio(id, label, options, value=None):
    return {"type": "radio", "id": id, "label": label, "options": options, "value": value or options[0]}
def checkbox(id, label, value=False):
    return {"type": "checkbox", "id": id, "label": label, "value": value}
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
def table(columns, title=None, caption=None, key=None, id=None, help=None, row_help=False):
    """Table. data: {"rows": [row, ...], "count"?: int, "columns"?: [str, ...]}   count -> small pill next to the title.
    data["columns"] replaces the spec columns for that run (the spec columns are what the page shows before a run); each row then has one cell per data column.
    cell: "text" | {"text", "url"} link (new tab) | {"badge": text, "tone": jp|report|good|warn|bad|neutral} | {"text", "help"} inline "?"
          | {"img": data-URL}; badge / text cells may also carry "url" and "help".   (the old {"text", "href"} link still works)
    row:  [cell, ...] | {"cells": [cell, ...], "help": str|list}.  row_help=True adds a trailing "?" column that shows each row's help."""
    return _extras({"type": "table", "title": title, "columns": columns, "caption": caption, "key": key, "id": id}, help=help, row_help=row_help)
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
def code(label, title=None, key=None, id=None, help=None):
    return _extras({"type": "code", "title": title, "label": label, "key": key, "id": id}, help=help)
def formats(label, title=None, key=None, id=None, help=None):
    return _extras({"type": "formats", "title": title, "label": label, "key": key, "id": id}, help=help)
def notes(title, key=None, id=None, help=None, collapsed=False):
    """Notes. data: {"items": [str, ...]}.  collapsed=True: no visible list; the items sit in a "?" next to the title
    (a small "Details" disclosure when title is None)."""
    return _extras({"type": "notes", "title": title, "key": key, "id": id}, help=help, collapsed=collapsed)
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
def grouptable(columns, title=None, id=None, help=None, caption=None, key=None, row_help=False):
    """Collapsible groups of rows that share the same columns (collapsed by default; chevron or header click toggles; keyboard accessible).
    data: {"groups": [{"title", "meta"?: short text on the right e.g. "0.66 · 3 methods", "help"?: str|list, "rows": [row, ...]}], "count"?: int}
    rows use the same cell / row formats as table (row_help=True adds the "?" column)."""
    return _extras({"type": "grouptable", "title": title, "columns": columns, "caption": caption, "key": key, "id": id}, help=help, row_help=row_help)


# ---- 模块 ----
def module(id, title, cls, one_liner, inputs, outputs, group="modules", subgroup=None, live=False, example=None):
    return {"id": id, "title": title, "cls": cls, "group": group, "subgroup": subgroup, "live": live,
            "one_liner": one_liner, "inputs": inputs, "outputs": outputs, "example": example}


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
# MassBank 里的真实咖啡因 [M+H]+ MS/MS（MSBNK-LCSB-LU119903，Q Exactive，CE 45；NIST 格式的 0-999 强度）：Library search 页的示例
CAFFEINE_MASSBANK_MS2 = ("56.0496 2\n69.0447 18\n83.0603 10\n97.0076 2\n97.0395 1\n109.0397 3\n110.0713 101\n111.0554 5\n112.0873 1\n122.0711 1\n"
                         "123.0428 3\n135.0668 1\n136.0747 1\n138.0663 999\n151.0978 5\n156.0768 2\n180.0642 2\n181.0719 3\n195.0878 847")
# Example names and the structure PubChem gives for each (checked by InChIKey). Tests use this table to run the examples offline.
EXAMPLE_STRUCTURES = {
    "caffeine": CAFFEINE,
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


MODULES = [
    module("identity", "Identity", "B", "Who is this compound: <b>names</b>, <b>database IDs</b>, synonyms, classification and background (online lookup).",
        [structure(), select("kind", "Input type", IDENTITY_KINDS)],
        [cards([("Name", ""), ("Formula", ""), ("PubChem CID", "")], "Resolved", id="resolved"),
         figure("2D depiction", id="depiction", span="half"),
         table(["Name type", "Value"], "Names", id="names"),
         table(["Database", "ID"], "Database IDs", id="database_ids"),
         table(["PubChem CID", "Title"], "Other matches for this name", id="candidates"),
         table(["Synonym"], "Synonyms", id="synonyms"),
         table(["System", "Level", "Term"], "Classification", id="classification",
               help=["ClassyFire and LOTUS records are source annotations; NPClassifier is a natural-product-oriented structure classifier.",
                     "LIPID MAPS rows appear for LMSD records. mzDesk rows come from the local lipid / metabolite rules. Different systems need not agree."]),
         table(["Origin", "Evidence", "Source"], "Origin", id="origin",
               help="Evidence of association or occurrence, not proof of exclusive origin. Human-associated records do not by themselves prove endogenous biosynthesis."),
         table(["Pathway", "Database", "Link"], "Pathways", id="pathways",
               help="KEGG, Reactome and HMDB / SMPDB annotations. Reactome uses the numeric ChEBI identifier and human pathway projection."),
         table(["Biofluid", "Normal concentration", "Unit", "Source"], "Biofluid concentrations", id="biofluids",
               help=["Read from a local HMDB XML only; values are never inferred. Set MZDESK_HMDB_XML to the unzipped HMDB metabolites XML path.",
                     "Normal ranges depend on cohort, matrix, method and units. Check the linked HMDB record before reuse."]),
         table(["Species", "Source"], "Species", id="species",
               help="Reported species-level LOTUS occurrences, deduplicated across taxonomic providers. Absence of a row is not evidence of absence."),
         table(["Vendor", "CAS", "Purity"], "Standards & availability", id="standards",
               help="PubChem vendor availability is a link-out, not a certified-standard claim. Product grade and purity must be checked with each vendor."),
         table(["Source", "Signal", "Hazard statements"], "GHS safety", id="ghs",
               help="Reporter-specific GHS classifications from PubChem. Reports can differ by jurisdiction, impurities and additives; this is not a safety data sheet.")],
        subgroup="molecular", live=True, example={"structure": _name("caffeine"), "kind": "Auto"}),

    module("structure_formats", "Structure & formats", "B", "Draw the structure in <b>2D / 3D</b> and export it as <b>SMILES, InChI, InChIKey or MOL</b>.",
        [structure()],
        [figure("2D depiction", id="depiction", span="half"), figure("3D viewer", id="viewer3d", span="half"),
         formats("Structure in every format", "Structure formats", id="formats")],
        subgroup="molecular", live=True, example={"structure": _name("caffeine")}),

    module("properties", "Properties & groups", "B", "<b>Formula</b>, <b>masses</b>, physicochemical descriptors and <b>functional groups</b>.",
        [structure()],
        [cards([("Formula", ""), ("Monoisotopic mass", "Da"), ("Average MW", "Da"), ("Net charge", "formal")], "Composition", id="summary"),
         table(["Group", "Descriptor", "Value"], "Descriptors", id="descriptors"),
         table(["Group", "Count", "Atoms", "SMARTS"], "Functional groups", id="functional_groups",
               help=["Transparent local SMARTS search; all matching rules are shown, including overlapping generic and specific groups.",
                     "Atoms are one-based indices in the parsed input; semicolons separate matches. At most eight atom sets are shown per row.",
                     "The curated rule library is broad but not a universally exhaustive functional-group ontology."])],
        subgroup="molecular", live=True, example={"structure": _name("caffeine")}),

    # ---- Detectability：7 页。只有 Platform 的脂质组 / 代谢组部分接了引擎，其余都是占位 ----
    module("platform_fit", "Platform", "A", "Which <b>workflow</b> can measure it: lipidomics, metabolomics, CIL LC-MS or GC-MS (decided mainly by <b>extraction</b> and <b>derivatization</b>).",
        [structure()],
        [portal("Platform", id="platform",
                help=["Rule-based screen of which workflow should work, not a measurement.",
                      "Yes: expected to work. Partial: only after derivatization. No: not expected. The ? on a tile gives its rule hits."]),
         reaction("Why and how", id="reactions",
                  help=["Shown when CIL or GC-MS applies. Choosing a reagent updates the structure to highlight only its reactive groups.",
                        "All four CIL channels show complete substrate and product structures, including every label; other reagents use functional-group schemes. Masses are monoisotopic."])],
        subgroup="detectability", live=True, example={"structure": _name("L-tyrosine")}),

    module("separation", "Separation", "B", "Which <b>column</b>, and what <b>charge</b> the molecule carries at each pH (MolGpKa pKa prediction; needs a one-time weights download).",
        [structure(), slider("ph", "Mobile phase pH", 1, 13, 0.1, 2.7)],
        [cards([("logD at pH", "", "logD = logP (Crippen) + log10(fraction of fully uncharged molecules). Ions are assumed not to partition and zwitterions count as charged, "
                                   "so logD is a lower bound at low or high pH."),
                ("Dominant charge state", "", "Most probable net charge at this pH, summed over the predicted sites (a zwitterion is net 0)."),
                ("Fraction neutral", "", "Share of the molecules with no charge at all; a zwitterion counts as charged.")],
               "pH behaviour", id="ph_summary", style="metric"),
         figure("logD and net charge vs pH", id="ph_curves", help=["Solid line, left axis: logD. Dashed line, right axis: net charge.", "Dotted line: the chosen pH."]),
         figure("Ionizable sites on the structure", id="sites_structure", span="half", help=["Coral: acid site. Teal: base site. Violet: both.", "Label: atom and pKa."]),
         chips("Separation mode", id="modes",
               help=["Simple rules on logD, net charge, TPSA, MW and stereocentres (thresholds: mzdesk/platform_fit/separation_rules.py).",
                     "Green: good. Amber: possible. Grey: poor or not needed. Hover or click a mode for its reason."]),
         table(["Site", "pKa", "Type"], "Ionizable sites", id="sites",
               help=["Site: element and atom number, counted from 1 as in the SMILES.", "Acid: loses H. Base: takes H.",
                     "Base sites on amide, sulfonamide and pyrrole-type N are ignored (see Notes)."]),
         notes("Notes", id="rule_hits", collapsed=True)],
        subgroup="detectability", live=True, example={"structure": _name("ibuprofen"), "ph": 2.7}),

    module("ionization", "Ionization", "B", "<b>ESI polarity</b> for LC-MS and <b>EI usability</b> for GC-MS.",
        [structure(), slider("ph", "Mobile phase pH", 1, 13, 0.1, 2.7)],
        [chips("Ionization", id="chips",
               help=["Polarity: the better ESI mode (+, −, or both); poor = weak in both",
                     "Simple rules at the chosen pH: pKa sites (MolGpKa; textbook group values without its weights), O/N counts, GC amenability and the NEIMS molecular ion",
                     "Thresholds: mzdesk/platform_fit/ionization_rules.py"]),
         figure("Ionizable sites", id="sites_structure", span="half", help=["Coral: acid site. Teal: base site. Violet: both.", "Label: atom and pKa.", "acids above pKa 10 are left out"]),
         notes("Rule hits", id="rule_hits", collapsed=True)],
        subgroup="detectability", live=True, example={"structure": _name("lidocaine"), "ph": 2.7}),

    module("adducts", "Adducts", "B", "Which <b>ions</b> appear (adducts, in-source losses, dimers) for the chosen polarity and mobile phase additive, and which is <b>strongest</b>. Rule-based: pKa sites (MolGpKa) and lipid class.",
        [structure(), radio("ion_mode", "Ion mode", ["Positive", "Negative"]),
         select("additive", "Mobile phase additive", ["Formic acid", "Ammonium formate", "Ammonium acetate", "Acetic acid", "None"]),
         checkbox("alkali", "Sodium / potassium present (glassware, salts, biological matrix)", True)],
        [cards([("Dominant ion", ""), ("m/z", "")], "Expected base ion", id="base", style="metric"),
         table(["Ion", "Formula", "m/z", "Likelihood"], "Ions", id="ions", row_help=True,
               help=["Ranked by expected intensity: a rule score (high / medium / low), not a measured intensity.",
                     "m/z computed as on the Theoretical m/z page."]),
         notes("Rule hits", id="rule_hits", collapsed=True)],
        subgroup="detectability", live=True,
        example={"structure": _name("POPC"),
                 "ion_mode": "Negative", "additive": "Ammonium formate", "alkali": True}),

    module("non_ms", "Non-MS detection", "B", "Can it be measured <b>without MS</b>: <b>UV</b> chromophore, <b>fluorescence</b>, or <b>universal detectors</b> (CAD / ELSD / RI). Rule-based.",
        [structure(), select("mobile_phase", "Mobile phase (changes the low-UV advice)", ["Water/acetonitrile + phosphate", "Water/methanol + formic acid"])],
        [portal("Detectors", id="detectors",
                help=["yes = good, partial = possible, no = poor",
                      "UV/PDA value: chromophore strength and estimated λmax. strong: aromatic or two or more conjugated π units; weak: one isolated C=O, C=C or amide; none: no π system",
                      "FLD value: native fluorescence (likely / possible / unlikely)",
                      "Rules and thresholds: mzdesk/platform_fit/non_ms_rules.py"]),
         figure("Chromophores", id="chromophore_structure", span="half", help=["blue: strong, amber: weak", "numbers = atoms (as in the SMILES, from 1)"]),
         notes("Rule hits", id="rule_hits", collapsed=True)],
        subgroup="detectability", live=True, example={"structure": _name("caffeine"), "mobile_phase": "Water/acetonitrile + phosphate"}),

    module("derivatization", "Derivatization", "B", "<b>Reactive groups</b> and the <b>derivatized molecule</b> for GC-MS, CIL LC-MS or LC-UV / FLD.",
        [structure(),
         multiselect("reagents", "Reagents",
                     ["GC · TMS (MSTFA / BSTFA)", "GC · MeOx + TMS", "GC · TBDMS (MTBSTFA)", "GC · Acylation (TFAA)",
                      "CIL · Dansyl chloride (amine / phenol)", "CIL · Dansyl chloride, base-activated (hydroxyl)",
                      "CIL · DmPA bromide (carboxyl)", "CIL · Dansylhydrazine (carbonyl)",
                      "UV/FLD · OPA + 2-mercaptoethanol (primary amine)", "UV/FLD · FMOC-Cl (amine)", "UV/FLD · DNPH (carbonyl)"],
                     ["GC · TMS (MSTFA / BSTFA)", "GC · MeOx + TMS", "CIL · Dansyl chloride (amine / phenol)",
                      "CIL · Dansyl chloride, base-activated (hydroxyl)", "CIL · DmPA bromide (carboxyl)", "CIL · Dansylhydrazine (carbonyl)"])],
        [table(["Reagent", "Labels", "Formula", "Mass", "Heavy mass"], "Derivatives", id="deriv_rows", row_help=True,
               help=["Full derivatization, one row per selected reagent that has a reactive group. Labels: reagent groups added.",
                     "Mass is monoisotopic; heavy mass (CIL only) is with a ¹³C₂-labelled reagent. The ? on a row lists the groups labelled, EI key ions and caveats."]),
         reaction("Reactions", id="deriv_reactions",
                  help=["All four CIL channels show complete substrate and product structures, including every label. Other reagents use functional-group schemes; ×n = number of groups.",
                        "Choose a reagent above the scheme; the substrate highlights only its reactive groups."])],
        subgroup="detectability", live=True, example={"structure": _name("D-glucose")}),

    module("known_methods", "Known methods", "B", "Is there a <b>published method</b> (pharmacopoeia, RepoRT)? If not, what do the <b>most similar compounds</b> use?",
        [structure(), multiselect("sources", "Sources", ["Japanese Pharmacopoeia", "RepoRT"])],
        [cards([("Match", "", "Exact: the compound itself has a published method. Similar: only similar compounds do. None: neither."),
                ("Methods found", "", "Exact matches plus the methods listed under the similar compounds."),
                ("Sources", "", "Which sources were searched; licence and citation.")], "Summary", id="summary", style="metric"),
         table(["Source", "Test / dataset", "Column", "Mobile phase", "Detector", "RT (min)"], "Exact matches", id="exact", row_help=True,
               help=["The compound itself: same connectivity key (first 14 characters of the InChIKey after largest fragment and neutralisation); stereo is noted in the row's ?.",
                     "The ? on a row gives match details, the reference and the printed JP conditions or the RepoRT gradient. RepoRT is LC-MS: RT is the compound's retention time in that dataset (min). JP18 has no RT."]),
         grouptable(["Source", "Test / dataset", "Column", "Mobile phase", "Detector", "RT (min)"], "Similar compounds", id="similar", row_help=True,
                    help="Click a compound to see its methods. RT is that compound's retention time.")],
        subgroup="detectability", live=True, example={"structure": _name("captopril"), "sources": ["Japanese Pharmacopoeia", "RepoRT"]}),

    # ---- MS1 ----
    module("ms1_search", "Mass search", "A", "Which <b>database compounds</b> match the precursor m/z, optionally filtered by the <b>isotope pattern</b>.",
        [precursor(), number("tol_ppm", "Mass tolerance (ppm)", 5.0, 0.1), multiselect("databases", "Exact-mass databases", DATABASES, ["HMDB", "PubChem"]),
         ms1(), select("resolution", "Resolution", list(RESOLUTIONS)), number("isotope_min", "Isotope score threshold", 0.8, 0.05)],
        [cards([("Compounds", "", "One row per PubChem compound (CID). Stereoisomers and salts have their own CIDs, so they count separately."),
                ("Formulas", "", "Different molecular formulas among the hits."),
                ("Closest hit", "", "The compound with the smallest |Δ ppm|; among equals, one listed in a chosen database comes first."),
                ("Databases searched", "", "Compounds PubChem finds in the mass window, per database (before the 100-per-database limit).")],
               "Summary", id="summary", style="metric"),
         grouptable(["Name", "IDs", "PubChem CID"], "Exact-mass hits", id="hits",
                    help=["Every database is searched online through PubChem, so hits come from PubChem's copy of each database "
                          "(HMDB, DrugBank, ChEBI, KEGG, LIPID MAPS; NORMAN is its Suspect List Exchange, which includes SusDat).",
                          "One group per formula and adduct, closest first. Δ ppm = (observed m/z − theoretical m/z) / theoretical m/z × 10<sup>6</sup>; "
                          "the theoretical m/z comes from the compound's monoisotopic mass and the adduct, as on the Theoretical page.",
                          "Auto tries every adduct of the ion mode. A compound is listed once, under its best-fitting adduct. "
                          "Permanently charged compounds (e.g. choline) only match [M]+.",
                          "At most 100 compounds are listed per database. PubChem alone can have thousands: compounds deposited by more "
                          "PubChem sources come first (depositor count), then lower CIDs (older, usually better known).",
                          "IDs are the databases' own IDs, taken from PubChem cross-references; PubChem has no ID for NORMAN.",
                          "With an isotope MS, each group's theoretical pattern (IsoSpecPy, merged to the chosen resolution) is compared with it: every theoretical peak "
                          "takes the most intense observed peak within the mass tolerance (Unit mass: 0.5 Da), each observed peak once; the score is the cosine "
                          "between theoretical abundances and matched intensities, and a theoretical peak with no match counts as 0.",
                          "Observed peaks that match no theoretical peak are ignored (they may be other ions). Groups below the threshold are hidden; "
                          "groups whose formula cannot be scored (e.g. isotope labels) are kept and marked N/A."])],
        subgroup="ms1", live=True,
        example={"precursor": _prec("195.0877"), "tol_ppm": "5", "databases": ["HMDB", "PubChem"],
                 "ms1": CAFFEINE_MS1, "resolution": "Q-TOF 30k", "isotope_min": "0.8"}),

    module("ms1_theory", "Theoretical", "A", "<b>Precursor m/z</b> and <b>isotope pattern</b> of a structure or formula for an adduct.",
        [structure(label="Compound or formula"),
         precursor(label="Adduct", adducts={"Positive": ADDUCTS_POS[1:], "Negative": ADDUCTS_NEG[1:]}, mz=False),
         select("resolution", "Resolution", list(RESOLUTIONS))],
        [cards([("Precursor m/z", "theoretical", "m/z of the ion made of the most abundant isotope of each element (M), as the Mass search page uses it."),
                ("Ion formula", ""), ("Most abundant peak", "")], "Theory", id="theory", style="metric"),
         table(["Peak", "m/z", "Relative abundance (%)"], "Isotope pattern", id="pattern",
               help=["Computed with IsoSpecPy from the ion formula (isotope fine structure), then merged to the chosen resolution: "
                     "peaks closer than the FWHM (m/z ÷ resolution) become one peak (intensity-weighted m/z, summed abundance); Unit mass merges by nominal mass.",
                     "M is the ion with the most abundant isotope of every element (<sup>11</sup>B, <sup>79</sup>Br, <sup>35</sup>Cl …). Lighter peaks are M−1 (e.g. <sup>10</sup>B).",
                     "Peaks that stay apart at the chosen resolution carry (a), (b) after the label. Peaks below 0.1 % of the largest are dropped."]),
         figure("Theoretical isotope pattern", id="pattern_plot", help="Hover a peak for its m/z and relative abundance.")],
        subgroup="ms1", live=True, example={"structure": _name("crisaborole"), "precursor": _prec("", adduct="[M+H]+"), "resolution": "Orbitrap 60k"}),

    # ---- Spectrum（subgroup id 仍是 ms2）----
    module("ms2_clean", "Clean & Vis", "B", "<b>Clean</b> an MS2 spectrum, <b>see</b> it, and get it in <b>every common format</b> (MS-DIAL, GNPS, MassBank, HMDB, MSP, MGF …).",
        [spectrum("spectrum", "MS2 spectrum (any common format)"), number("precursor_mz", "Precursor m/z (optional)", step=0.0001),
         radio("ion_mode", "Ion mode (for MSP / MGF headers)", ["Positive", "Negative"]),
         number("noise_pct", "Noise threshold (% of base peak)", 0, 0.1), select("normalize", "Normalize", NORMALIZE),
         checkbox("remove_precursor", "Remove precursor"), number("precursor_tol", "Precursor tolerance (Da)", 0.02, 0.001),
         select("sort_by", "Sort by", SORT_BY)],
        [cards([("Detected format", ""), ("Peaks in", ""), ("Peaks out", ""), ("Removed", "noise · precursor")], id="summary"),
         figure("Cleaned spectrum", id="spectrum", help=["Grey: removed peaks.", "Hover a peak for details."]),
         table(["m/z", "Intensity", "Rel. intensity (%)", "Kept", "Reason"], "Peaks", "Rel. intensity: % of the input base peak. Reason: noise, precursor", id="peaks"),
         formats("Cleaned spectrum in every format", "Formats", id="formats")],
        subgroup="ms2", live=True,
        example={"spectrum": _paste("42.0338:5 56.0495:4 61.0012:0.3 69.0447:6 83.0604:9 97.5521:0.6 110.0713:18 138.0662:100 195.0877:32"),
                 "precursor_mz": "195.0877", "noise_pct": "1", "normalize": "Base peak = 100", "remove_precursor": True}),

    module("ms2_ei", "EI MS", "A", "Predict the <b>70 eV electron-ionization spectrum</b> of a structure (NEIMS). One ion type: [M]+•.",
        [structure()],
        [cards([("Ion", ""), ("Formula", ""), ("Molecular ion m/z", ""), ("Base peak m/z", "")], id="summary"),
         figure("Predicted EI spectrum", id="spectrum", help=["Unit mass resolution.", "Hover a peak for details."]),
         table(["m/z", "Rel. intensity (%)"], "Peaks", id="peaks"),
         formats("Predicted spectrum in every format", "Formats", id="formats")],
        subgroup="ms2", live=True, example={"structure": _name("caffeine")}),

    module("ms2_esi", "ESI MS/MS", "A", "Predict the <b>ESI MS/MS spectrum</b> of a structure at <b>10 / 20 / 40 eV</b> (CFM-ID 4; MSRB rules for lipids). Needs Docker.",
        [structure(), select("adduct", "Adduct", ESI_ADDUCTS, "[M+H]+")],
        [cards([("Engine", ""), ("Precursor m/z", ""), ("Ion formula", ""), ("Peaks", "")], id="summary"),
         figure("Low energy (10 eV)", id="spectrum_low"), figure("Medium energy (20 eV)", id="spectrum_medium"),
         figure("High energy (40 eV)", id="spectrum_high"),
         table(["Energy", "m/z", "Rel. intensity (%)", "Fragment"], "Peaks", "Fragment: formula · SMILES from the engine", id="peaks"),
         formats("Predicted spectra (one record per energy)", "Formats", id="formats")],
        subgroup="ms2", live=True, example={"structure": _name("caffeine"), "adduct": "[M+H]+"}),

    module("ms2_library", "Library search", "A", "Search a spectrum against <b>MassBank</b> and <b>GNPS</b> with the similarity metrics you choose.",
        [radio("spectrum_type", "Spectrum type", ["MS/MS", "EI"]), spectrum(label="Query spectrum"),
         precursor(label="Precursor m/z (MS/MS only)"),
         number("tol_ppm", "Precursor tolerance (ppm, MS/MS only)", 10.0, 0.1),
         number("tol_da", "Fragment tolerance (Da; EI: 0.5)", 0.02, 0.001),
         multiselect("libraries", "Libraries", LIBRARIES), multiselect("metrics", "Metrics", LIBRARY_METRICS, LIBRARY_DEFAULT_METRICS),
         number("threshold", "Score threshold", 0.7, 0.05)],
        [cards([("Best hit", "name", "The first row of the table: the hit with the highest score on the ranking metric."),
                ("Rank score", "ranking metric"), ("Hits ≥ threshold", "n"),
                ("Candidates scored", "n"), ("Libraries", "version")], id="summary", style="metric"),
         table(["Rank", "Library", "Accession", "Name", "Adduct", "Precursor m/z", "Δ ppm", "Matched peaks"] + LIBRARY_DEFAULT_METRICS + ["Instrument / CE"],
               "Library hits", id="hits",
               help=["A spectrum is a hit when ANY selected metric is at or above the threshold (N/A scores are ignored). The table shows the best 100.",
                     "Hits are ranked by the first selected metric that applies, in the order of the Metrics list (not the order you clicked them); a metric that is N/A for every hit is skipped. Ties: more matched peaks first.",
                     "Matched peaks: peak pairs matched one-to-one within the fragment tolerance (greedy, as in CosineGreedy). Δ ppm: query minus library precursor.",
                     "MS/MS: the precursor m/z is required (a header in an MGF / MSP spectrum counts). Candidates are the library MS/MS spectra within the precursor tolerance, in the same ion mode "
                     "(and with the chosen adduct, unless Auto); every candidate is scored with every selected metric. With both libraries selected, GNPS copies of MassBank spectra are hidden when the MassBank original is also a candidate.",
                     "EI: the precursor is not used. Fragments are matched at unit mass (0.5 Da). Candidates = top 5000 by entropy similarity; a spectrum with low entropy but high other scores can be missed.",
                     "Deep-learning scores need torch / gensim and a one-time model download; they are N/A where those are missing."]),
         figure("Mirror: query vs top hit", id="mirror", help=["Up: query. Down: the top-ranked hit.", "Dark: matched peaks."])],
        subgroup="ms2", live=True,
        example={"spectrum_type": "MS/MS", "ms2": _paste(CAFFEINE_MASSBANK_MS2), "precursor": _prec("195.0877")}),

    module("ms2_similarity", "Similarity & mirror", "A", "How alike are <b>two spectra</b>: 14 classic scores, optional deep-learning scores, <b>mirror plot</b>.",
        [spectrum("query", "Query spectrum"), number("query_precursor", "Query precursor m/z (optional)", step=0.0001),
         spectrum("reference", "Reference spectrum"), number("reference_precursor", "Reference precursor m/z (optional)", step=0.0001),
         number("tol_da", "Fragment tolerance (Da)", 0.02, 0.001), number("threshold", "Score threshold", 0.7, 0.05),
         radio("ion_mode", "Ion mode", ["Positive", "Negative"]),
         checkbox("deep", "Deep-learning scores (MS2DeepScore, Spec2Vec; first use downloads the models)")],
        [cards([("Primary score", "spectral entropy, weighted"), ("Metrics ≥ threshold", "n / N"), ("Matched peaks", "pairs"),
                ("Query matched", "%"), ("Reference matched", "%")], id="summary"),
         figure("Mirror plot", id="mirror", help=["Up: query. Down: reference.", "Dark: matched peaks."]),
         table(["Category", "Metric", "Score", "Matched peaks", "≥ threshold", "Note"], "Scores", id="scores")],
        subgroup="ms2", live=True,
        example={"query": _paste("197.04643 10\n255.22775 92\n256.23184 31\n277.60822 10\n283.24814 20"),
                 "reference": _paste("197.04643 8\n255.22775 100\n256.23184 28\n277.60822 7\n283.24814 22")}),

    module("ms2_fragments", "Fragments (MetFrag)", "A", "Which <b>peaks</b> the candidate structure explains, and with which <b>fragment</b> (MetFrag).",
        [structure(label="Candidate structure"), precursor(label="Precursor m/z (optional; blank = theoretical)", adducts=METFRAG_ADDUCTS), spectrum(),
         slider("tree_depth", "Tree depth (bond cleavages)", 1, 4, 1, 2),
         number("tol_abs", "Absolute tolerance (Da)", 0.002, 0.001), number("tol_ppm", "Relative tolerance (ppm)", 10.0, 1)],
        [cards([("Explained peaks", "n / N"), ("Explained intensity", "%"), ("Theoretical precursor m/z", ""), ("Ion formula", "")], id="summary"),
         figure("Observed spectrum", id="spectrum", help="Hover a peak: fragment formula, Δ ppm, neutral loss and fragment structure."),
         table(["Observed m/z", "Rel. intensity (%)", "Kind", "Ion formula", "Theoretical m/z", "Δ ppm", "Neutral loss", "Fragment structure"], id="peaks")],
        subgroup="ms2", live=True,
        example={"structure": _name("DOPG"), "precursor": _prec("", "Negative", "[M-H]-"), "ms2": _paste(PG_36_2_MS2)}),

    module("ms1_formula", "Formula (BUDDY)", "A", "<b>Rank molecular formulas</b> for the precursor from its MS1 isotopes and MS2 (BUDDY).",
        [precursor(adducts=BUDDY_ADDUCTS), ms1(), spectrum(label="MS2 spectrum (optional)"),
         select("instrument", "Instrument", ["Orbitrap", "Q-TOF", "FT-ICR"]), checkbox("halogen", "Include halogens (F, Cl, Br, I)"),
         text("formula", "Formula (optional: skip ranking, explain the MS2 with it)", placeholder="C8H10N4O2")],
        [cards([("Top formula", ""), ("Adduct used", ""), ("MS1 isotope peaks used", "")], id="summary"),
         table(["Rank", "Formula", "Theoretical m/z", "Δ ppm", "Probability", "FDR", "MS1 isotope similarity"], "Formula candidates", id="candidates"),
         figure("MS2 fragment formulas · top candidate", id="fragments_plot", help="Hover a peak for its formula, Δ ppm and neutral loss."),
         table(["m/z", "Rel. intensity (%)", "Fragment formula", "Δ ppm", "Neutral loss", "Other fits"], id="fragments")],
        subgroup="ms2", live=True,
        example={"precursor": _prec("195.0877"), "ms1": CAFFEINE_MS1, "ms2": _paste(CAFFEINE_MS2)}),

    module("ms2_denovo", "De novo", "A", "<b>Structure candidates</b> when no structure is given.",
        [precursor(), spectrum(), text("formula", "Formula (optional)")],
        [table(["Rank", "SMILES", "Formula", "Score", "Tool"], "De novo candidates")],
        subgroup="ms2"),

    module("exclusivity", "Exclusivity", "A", "Is there <b>another explanation</b>: same formula, isobaric interference, ion relations, isomers, metabolites, lipid acyl combinations.",
        [structure(label="Compound or formula", optional=True), precursor(), number("tol_ppm", "Tolerance (ppm)", 5.0, 0.1),
         number("mz2", "Second m/z (relation check, optional)", step=0.0001), multiselect("databases", "Databases", DATABASES, ["HMDB", "PubChem"]),
         spectrum(label="MS2 spectrum (optional, for isomer check)")],
        [table(["InChIKey (first block)", "Name", "Databases", "Stereoisomers collapsed"], "Same-formula candidates"),
         table(["Formula", "Name", "Database", "Exact m/z", "Δ ppm"], "Isobaric interferences", "Different formula, same m/z within tolerance"),
         table(["Relation", "Δ mass", "Explanation", "Confidence"], "Relation between two m/z", "Adduct / in-source fragment / dimer / neutral loss"),
         table(["Isomer SMILES", "Distinguishable by MS2", "Distinguishing fragments"], "Isomers"),
         table(["Direction", "Parent", "Product", "Reaction", "Tool"], "Metabolite / transformation products", "BioTransformer, SyGMa"),
         table(["Combination", "Diagnostic ions", "Observed"], "Lipids: sum composition → acyl combinations")], subgroup="tools"),

    # ---- Orthogonal ----
    module("rt", "RT / RI", "A", "Does the <b>retention</b> agree: RT / RI prediction and projection across systems.",
        [structure(), select("system", "Chromatographic system", ["RP C18 generic", "HILIC generic", "GC (Kovats)", "Custom (upload calibration)"]),
         number("rt_min", "Observed RT (min, optional)", step=0.01), file("calibration", "Calibration set for projection (optional, CSV: compound, RT)", ".csv")],
        [cards([("Predicted RT / RI", ""), ("Observed", ""), ("Δ", "")], "Retention"), table(["System", "Predicted", "Projected", "Observed", "Δ"])],
        subgroup="orthogonal"),

    module("ccs", "CCS", "A", "Does the <b>ion mobility</b> agree: CCS prediction and lookup.",
        [structure(), select("adduct", "Adduct", ADDUCTS_POS[1:] + ADDUCTS_NEG[1:]), number("ccs", "Observed CCS (Å², optional)", step=0.1)],
        [cards([("Predicted CCS", "Å²"), ("Database CCS", "Å²"), ("Observed", "Å²"), ("Δ", "%")], "CCS"), table(["Adduct", "Predicted", "Database", "Observed", "Δ %"])],
        subgroup="orthogonal"),

    module("mrm", "MRM / targeted transitions", "B", "Pick <b>precursor–product pairs</b> and <b>collision energies</b> for a targeted method.",
        [radio("start_from", "Start from", ["Structure (predicted fragments)", "Library spectrum", "My spectrum"]), structure(optional=True), spectrum(label="MS2 spectrum (if starting from my spectrum)"),
         radio("ion_mode", "Ion mode", ["Positive", "Negative"]), select("instrument", "Instrument", ["Triple quadrupole (generic)", "QTRAP", "Q-TOF (PRM)", "Orbitrap (PRM)"]),
         number("n_transitions", "Transitions", 3, 1), select("ce_unit", "CE unit", ["eV", "NCE"])],
        [cards([("Predicted RT", "min"), ("Polarity", ""), ("Transitions", "n")]),
         table(["Precursor m/z", "Product m/z", "CE", "Role", "Fragment formula", "Source", "Relative intensity"], caption="Role: quantifier / qualifier. Source: library spectrum or predicted"),
         button("Export method table (CSV)")], subgroup="tools"),

    module("calculator", "Calculator", "B", "Everyday calculations. Own inputs, nothing shared.",
        [heading("ppm error"), number("ppm_theoretical", "Theoretical m/z", step=0.0001), number("ppm_observed", "Observed m/z", step=0.0001),
         heading("Isotope labelling"), text("lab_smiles", "SMILES"), text("lab_label", "Label", placeholder="13C6, d4, 15N2"), select("lab_adduct", "Adduct", ADDUCTS_POS[1:]),
         heading("Concentration"), number("conc_value", "Value"), select("conc_from", "From", ["µM", "mM", "ng/mL", "µg/mL", "mg/L"]), select("conc_to", "To", ["µM", "mM", "ng/mL", "µg/mL", "mg/L"]),
         heading("Collision energy"), number("ce_nce", "NCE"), number("ce_ev", "eV"), number("ce_mz", "Precursor m/z", step=0.0001), number("ce_z", "Charge", 1, 1),
         heading("GC retention index (Kovats)"), grid("alkanes", "Alkane series", ["Carbon number", "RT"]), number("ri_rt", "Analyte RT", step=0.01)],
        [cards([("ppm error", "ppm")], "ppm error"),
         table(["Isotopologue", "m/z", "Abundance"], "Isotope labelling"), table(["Internal standard", "Analyte M+n", "Crosstalk %"]),
         cards([("Converted", "")], "Concentration"), table(["Target conc.", "Stock conc.", "Final volume", "Stock volume", "Diluent volume"], caption="Dilution and stock preparation"),
         cards([("NCE ↔ eV", "")], "Collision energy"), cards([("Kovats RI", "")], "GC retention index")], subgroup="tools"),

    module("report", "Report", "A/B", "Collect the evidence, give a <b>confidence level</b>, export. Plus name normalization, figures, MassBank record.",
        [select("source", "Evidence from", ["Last pipeline run", "Upload report JSON"]), file("report_json", "Report JSON (optional)", ".json")],
        [cards([("Schymanski level", "1–5"), ("Evidence items", "n"), ("Conflicts", "n")], "Annotation report (A)"),
         table(["Module", "Evidence", "Result", "Supports", "Weight"], caption="One row per piece of evidence, in pipeline order"),
         table(["Parameter", "Value", "Module"], caption="Provenance"), table(["Tool", "Version", "Citation"], caption="What to cite"),
         button("Export JSON"), button("Export CSV"), button("Export PDF"),
         table(["Convention", "Name"], "Name normalization (B)", "Goslin lipid shorthand, LIPID MAPS, RefMet"),
         button("Spectrum with structures (SVG)"), button("Mirror plot (SVG)"),
         code("MassBank record", "MassBank record generator (B)")], subgroup="tools"),

    # ---- DMPK（非核心）----
    module("met_shift", "Metabolite mass-shift search", "B", "Search a peak or feature list for <b>phase I / II metabolites</b> of the parent.",
        [structure(label="Parent compound or formula"), precursor(),
         grid("features", "Peak / feature list", ["m/z", "RT", "Intensity"], "Paste or upload: one feature per line"), file("features_file", "Feature table (CSV, optional)", ".csv"),
         multiselect("shifts", "Mass shifts", ["+O (hydroxylation)", "−2H (dehydrogenation)", "+2H (reduction)", "−CH2 (demethylation)", "+O −CH2", "+C6H8O6 glucuronide", "+SO3 sulfate", "+C2H2O acetyl", "+C10H15N3O6S GSH", "+C2H3NO2 glycine"]),
         number("tol_ppm", "Tolerance (ppm)", 10.0, 0.1)],
        [table(["Shift", "Δ mass", "Expected m/z", "Found m/z", "Δ ppm", "RT", "Intensity", "Candidate reaction"]),
         figure("Features with metabolite hits marked (m/z vs RT)")], subgroup="dmpk"),

    module("som", "Site of metabolism", "B", "Which <b>atoms</b> are metabolized first.",
        [structure(), multiselect("enzymes", "Enzymes", ["CYP3A4", "CYP2D6", "CYP2C9", "CYP2C19", "CYP1A2", "UGT", "All"], ["All"])],
        [table(["Atom", "Rank", "Score", "Reaction type", "Enzyme", "Tool"], caption="SMARTCyp / FAME"), figure("Structure with predicted sites highlighted")], subgroup="dmpk"),

    module("stability", "Metabolic stability", "B", "<b>t½</b> and <b>intrinsic clearance</b> from an in vitro incubation.",
        [grid("timecourse", "Time course", ["Time (min)", "% remaining"]), select("system", "System", ["Liver microsomes", "Hepatocytes", "S9", "Plasma"]),
         number("protein", "Protein (mg/mL) or cells (10⁶/mL)", 0.5, 0.1), select("species", "Species", ["Human", "Rat", "Mouse", "Dog", "Monkey"])],
        [cards([("t½", "min"), ("CLint", "µL/min/mg or µL/min/10⁶ cells"), ("Scaled CLint", "mL/min/kg"), ("Predicted hepatic CL", "mL/min/kg, well-stirred"), ("Extraction ratio", "")]),
         figure("ln(% remaining) vs time with fit"), table(["Parameter", "Value", "Unit", "Assumption"], caption="Scaling factors used")], subgroup="dmpk"),

    module("nca", "NCA pharmacokinetics", "B", "<b>Non-compartmental parameters</b> from a concentration–time profile.",
        [grid("profile", "Concentration–time", ["Time (h)", "Concentration (ng/mL)"]), number("dose", "Dose (mg)", step=0.01),
         select("route", "Route", ["IV bolus", "IV infusion", "Oral"]), number("lambda_points", "Points for λz", 3, 1)],
        [cards([("Cmax", "ng/mL"), ("Tmax", "h"), ("AUC0-t", "ng·h/mL"), ("AUC0-∞", "ng·h/mL"), ("t½", "h"), ("CL/F", "L/h"), ("Vz/F", "L"), ("MRT", "h")]),
         figure("Concentration–time (semi-log) with λz fit"), table(["Parameter", "Value", "Unit", "Method"])], subgroup="dmpk"),

    module("sil_is", "Isotope-labelled internal standard", "B", "Which <b>labelling</b> to buy: mass shift, stable positions, crosstalk.",
        [structure(label="Analyte"), spectrum(label="Analyte MS2 (optional)"), number("min_shift", "Minimum mass shift (Da)", 3, 1), checkbox("check_vendors", "Check vendor availability", True)],
        [table(["Labelling", "Mass shift", "Positions", "Exchangeable-H risk", "Metabolic-loss risk", "Label kept in quantifier fragment", "Crosstalk with analyte M+n (%)", "Vendor"]),
         figure("Structure with labelled positions highlighted")], subgroup="dmpk"),

    # ---- Inorganic MS（非核心）----
    module("icp_isotopes", "Isotopes & internal standard", "B", "Choose an ICP-MS isotope from natural abundance, atomic isobars and matrix polyatomics; match an internal standard by mass and first ionization energy.",
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
        subgroup="inorganic", live=True, example={"element": "As", "matrix": list(MATRIX_ELEMENTS)}),

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
        example={"target": "75As", "matrix": list(MATRIX_ELEMENTS), "max_atoms": "3", "charge_states": ["1+", "2+"]}),

    module("icp_idms", "Isotope dilution (IDMS)", "B", "Single-spike, mole-based isotope dilution with composition-specific molar masses, ratio sensitivity and an optimal-blend hint.",
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
                 "spike_concentration": "10", "spike_unit": "µg/g", "m_x": "1", "m_y": "0.5", "r_b": "1"}),

    module("icp_massbias", "Mass bias & δ values", "B", "Correct isotope ratios with three mass-bias laws or standard–sample bracketing, and convert between ratios and δ values.",
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
                 "delta_sample": "1.002", "delta_standard": "1", "delta_per_mil": "2"}),

    module("icp_sp", "Single-particle ICP-MS", "B", "Convert detected particle-event signals to spherical diameter and number concentration using dissolved-standard calibration and transport efficiency.",
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
                 "intensities": "5054.728\n6200\n4500", "n_events": "3", "acquisition_s": "60", "threshold": "3000"}),
]


# ---- Pipelines ----
def stage(title, modules, exclude=None):
    """modules: 一个或多个模块 id，按顺序显示它们的输出块；exclude 按块的 key 过滤。"""
    return {"title": title, "modules": [modules] if isinstance(modules, str) else list(modules), "exclude": exclude or []}


PIPELINES = [
    {"id": "library_search", "title": "Library search", "one_liner": "I have a <b>spectrum</b> and no structure. Find it in the <b>libraries</b> and check the hit.",
     "inputs": [spectrum(), precursor(), multiselect("libraries", "Libraries", LIBRARIES), number("tol_ppm", "Tolerance (ppm)", 10.0, 0.1)],
     "stages": [stage("Spectrum cleanup", "ms2_clean"), stage("MS2 library search", ["ms2_library", "ms2_similarity"]),
                stage("MS1 check of the hits", "ms1_theory"), stage("Exclusivity", "exclusivity"), stage("Report", "report")]},
    {"id": "de_novo", "title": "De novo annotation", "one_liner": "I have a spectrum, the libraries do not have it. Infer <b>formula</b>, then <b>structure candidates</b>.",
     "inputs": [precursor(), ms1(), spectrum(), number("rt_min", "RT (min, optional)", step=0.01), number("ccs", "CCS (Å², optional)", step=0.1),
                text("elements", "Element constraints", "C0-50 H0-100 N0-10 O0-20 P0-3 S0-3"),
                multiselect("databases", "Candidate databases", DATABASES, ["HMDB", "PubChem"])],
     "stages": [stage("Spectrum cleanup", "ms2_clean"), stage("MS1: formula", "ms1_formula"), stage("MS2: fragment explanation / de novo", ["ms2_fragments", "ms2_denovo"]),
                stage("Exclusivity", "exclusivity"), stage("Orthogonal", ["rt", "ccs"]), stage("Report", "report")]},
    {"id": "verify", "title": "Result verification", "one_liner": "I have a <b>candidate structure</b>. Does <b>every piece of evidence</b> support it?",
     "inputs": [structure(label="Candidate structure"), precursor(), ms1(), spectrum(),
                number("rt_min", "RT (min, optional)", step=0.01), number("ccs", "CCS (Å², optional)", step=0.1)],
     "stages": [stage("Molecular", ["structure_formats", "properties"]), stage("Platform fit", "platform_fit"), stage("MS1", "ms1_theory"),
                stage("MS2", ["ms2_esi", "ms2_library", "ms2_similarity", "ms2_fragments"]),
                stage("Exclusivity", "exclusivity"), stage("Orthogonal", ["rt", "ccs"]), stage("Report", "report")]},
    {"id": "feature_prediction", "title": "MS feature prediction", "one_liner": "I have a <b>structure</b> and no spectrum. Predict <b>what the instrument will see</b>: adducts, precursor m/z, isotope pattern, MS2.",
     "inputs": [structure(), radio("ion_mode", "Ion mode", ["Positive", "Negative", "Both"]), multiselect("adducts", "Adducts", ["Auto"] + ADDUCTS_POS[1:] + ADDUCTS_NEG[1:], ["Auto"]),
                select("instrument", "Instrument / CE", ["Q-TOF 20 eV", "Q-TOF 40 eV", "Orbitrap HCD NCE 30", "Triple quadrupole"]), select("out_format", "Export format", ["MSP", "MGF", "MassBank record", "JSON"])],
     "stages": [stage("Molecular", ["structure_formats", "properties"]), stage("Adduct prediction", "adducts"),
                stage("MS1 theory: precursor m/z and isotope pattern", "ms1_theory"),
                stage("Spectrum prediction (ESI MS/MS, EI MS)", ["ms2_esi", "ms2_ei"]),
                stage("Export", "ms2_clean")]},
]

GROUPS = [{"id": "modules", "title": "Modules"}, {"id": "pipelines", "title": "Pipelines"}]
# Modules 下的小类，侧栏按这里的顺序列出（组内已接引擎的模块排前面）；Pipelines 排在 Modules 之后
SUBGROUPS = [{"id": "molecular", "title": "Molecular"}, {"id": "detectability", "title": "Detectability"},
             {"id": "ms1", "title": "MS1"}, {"id": "ms2", "title": "Spectrum"}, {"id": "orthogonal", "title": "Orthogonal"}, {"id": "dmpk", "title": "DMPK"},
             {"id": "inorganic", "title": "Inorganic MS"}, {"id": "tools", "title": "Tools"}]


def spec():
    return {"groups": GROUPS, "subgroups": SUBGROUPS, "pipelines": PIPELINES, "modules": MODULES}
