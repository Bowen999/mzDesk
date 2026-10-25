"""页面定义（唯一来源）。前端按这里的 pipelines / modules 渲染输入和输出格式。

字段类型：structure, precursor, spectrum, ms1, number, text, select, multiselect, radio,
checkbox, slider, file, table, heading。
输出块类型：table, cards, figure, code, formats, notes, button, heading。
  formats：同一份内容的多种格式，数据 {"formats": {格式名: 文本}}，前端画成一排标签 + 一个代码框。
  id：后端结果按 id 填进这个块；没有 id 的块只画格式。
  key：pipeline 的 stage 用 exclude 按 key 过滤块。
模块：group 是侧栏大类（pipelines / modules）；subgroup 是 Modules 下的小类（molecular / detectability / ms1 / ms2 / orthogonal / dmpk；ms2 的显示名是 Spectrum），
  None = 直接列在 Modules 下。live=True 表示已接引擎（后端 mzdesk/api/runners.py 里有对应 runner）。
  example：一组可直接运行的示例输入，前端 "Load example" 按钮填入。
"""
from mzdesk.shared.ms2_clean import NORMALIZE, SORT_BY
# 注意：resolve.py 导入时不联网，所以这里可以直接导入 KINDS；不要在 spec 里导入任何联网代码
from mzdesk.identity.resolve import KINDS as IDENTITY_KINDS
# cfmid_engine 顶层只用标准库，导入它不会拉进 rdkit / Docker
from mzdesk.ms2.cfmid_engine import ADDUCTS as ESI_ADDUCTS

ADDUCTS_POS = ["Auto", "[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M]+", "[M+H-H2O]+", "[2M+H]+"]
ADDUCTS_NEG = ["Auto", "[M-H]-", "[M+Cl]-", "[M+HCOO]-", "[M+CH3COO]-", "[M-H-H2O]-", "[2M-H]-"]
# 引擎各自支持的加合物（Auto = [M+H]+ / [M-H]-）
METFRAG_ADDUCTS = {"Positive": ["Auto", "[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M]+"],
                   "Negative": ["Auto", "[M-H]-", "[M+Cl]-", "[M+HCOO]-", "[M+CH3COO]-"]}
BUDDY_ADDUCTS = {"Positive": ["Auto", "[M+H]+", "[M+Na]+", "[M+NH4]+", "[M+K]+", "[M+H-H2O]+"],
                 "Negative": ["Auto", "[M-H]-", "[M+Cl]-", "[M+HCOO]-", "[M+CH3COO]-", "[M+Br]-"]}
LIBRARIES = ["MassBank", "MoNA", "GNPS"]
DATABASES = ["HMDB", "PubChem", "LIPID MAPS", "KEGG", "ChEBI", "DrugBank", "NORMAN SusDat"]


# ---- 字段 ----
def structure(id="structure", label="Compound", optional=False):
    return {"type": "structure", "id": id, "label": label, "optional": optional}
def precursor(id="precursor", adduct=True, label="Precursor m/z", adducts=None):
    return {"type": "precursor", "id": id, "adduct": adduct, "label": label,
            "adducts": adducts or {"Positive": ADDUCTS_POS, "Negative": ADDUCTS_NEG}}
def spectrum(id="ms2", label="MS2 spectrum", multi=False):
    return {"type": "spectrum", "id": id, "label": label, "multi": multi}
def ms1(id="ms1"):
    return {"type": "ms1", "id": id, "label": "MS1 isotope pattern (m/z intensity per line, optional)"}
def number(id, label, value=None, step=None, placeholder=None):
    return {"type": "number", "id": id, "label": label, "value": value, "step": step, "placeholder": placeholder}
def text(id, label, value="", placeholder=None):
    return {"type": "text", "id": id, "label": label, "value": value, "placeholder": placeholder}
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
def table(columns, title=None, caption=None, key=None, id=None):
    return {"type": "table", "title": title, "columns": columns, "caption": caption, "key": key, "id": id}
def cards(items, title=None, key=None, id=None):
    return {"type": "cards", "title": title, "items": [{"label": l, "unit": u} for l, u in items], "key": key, "id": id}
def figure(label, key=None, id=None):
    return {"type": "figure", "label": label, "key": key, "id": id}
def code(label, title=None, key=None, id=None):
    return {"type": "code", "title": title, "label": label, "key": key, "id": id}
def formats(label, title=None, key=None, id=None):
    return {"type": "formats", "title": title, "label": label, "key": key, "id": id}
def notes(title, key=None, id=None):
    return {"type": "notes", "title": title, "key": key, "id": id}
def button(label, key=None):
    return {"type": "button", "label": label, "key": key}
def head(title, key=None):
    return {"type": "heading", "title": title, "key": key}


# ---- 模块 ----
def module(id, title, cls, one_liner, inputs, outputs, group="modules", subgroup=None, live=False, example=None):
    return {"id": id, "title": title, "cls": cls, "group": group, "subgroup": subgroup, "live": live,
            "one_liner": one_liner, "inputs": inputs, "outputs": outputs, "example": example}


# 示例输入（与前端 initialValues 的形状一致）
def _smiles(s):
    return {"mode": "SMILES / InChI", "value": s}
def _formula(f):
    return {"mode": "Name / ID / Formula", "value": f}
def _prec(mz="", ion_mode="Positive", adduct="Auto"):
    return {"precursor_mz": mz, "ion_mode": ion_mode, "adduct": adduct}
def _paste(peaks):
    return {"source": "Paste", "value": peaks}


CAFFEINE = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
CAFFEINE_MS2 = "42.0338 5\n56.0495 4\n69.0447 6\n83.0604 9\n110.0713 18\n138.0662 100\n195.0877 32"
CAFFEINE_MS1 = "195.0877 100\n196.0902 10.3\n197.0922 0.9"
PG_36_2 = r"CCCCCCCC\C=C/CCCCCCCC(=O)OCC(COP(O)(=O)OCC(O)CO)OC(=O)CCCCCCC\C=C/CCCCCCCC"
PG_36_2_MS2 = ("139.00977 3\n152.99622 24\n168.04965 3\n279.22757 48\n280.2258 27\n281.24707 1101\n"
               "282.24652 24\n283.26199 14\n491.25638 10\n509.2948 17\n599.81622 3\n713.46429 10\n773.53589 416")


MODULES = [
    module("identity", "Identity", "B", "Who is this compound: names, database IDs, synonyms, classification and background (online lookup).",
        [structure(), select("kind", "Input type", IDENTITY_KINDS)],
        [cards([("Name", ""), ("Formula", ""), ("PubChem CID", "")], "Resolved", id="resolved"),
         figure("2D depiction", id="depiction"),
         table(["Name type", "Value"], "Names", id="names"),
         table(["Database", "ID"], "Database IDs", id="database_ids"),
         table(["PubChem CID", "Title"], "Other matches for this name", id="candidates"),
         table(["Synonym"], "Synonyms", id="synonyms"),
         table(["System", "Level", "Term"], "Classification", "ClassyFire, NPClassifier, LIPID MAPS"),
         table(["Origin", "Evidence", "Source"], "Origin", "Endogenous / drug / food / pollutant / microbial / plant"),
         table(["Pathway", "Database", "Link"], "Pathways"),
         table(["Biofluid", "Normal concentration", "Unit", "Source"], "Biofluid concentrations"),
         table(["Species", "Source"], "Species", "LOTUS"),
         table(["Vendor", "CAS", "Purity", "GHS hazard"], "Standards & safety")],
        subgroup="molecular", live=True, example={"structure": {"mode": "Name / ID / Formula", "value": "caffeine"}, "kind": "Auto"}),

    module("structure_formats", "Structure & formats", "B", "Draw the structure in 2D / 3D and export it as SMILES, InChI, InChIKey or MOL.",
        [structure()],
        [figure("2D depiction", id="depiction"), figure("3D viewer", id="viewer3d"),
         formats("Structure in every format", "Structure formats", id="formats")],
        subgroup="molecular", live=True, example={"structure": _smiles(CAFFEINE)}),

    module("properties", "Properties & groups", "B", "Formula, masses, physicochemical descriptors and functional groups.",
        [structure()],
        [cards([("Formula", ""), ("Monoisotopic mass", "Da"), ("Average MW", "Da"), ("Net charge", "formal")], "Composition", id="summary"),
         table(["Group", "Descriptor", "Value"], "Descriptors", id="descriptors"),
         table(["Group", "Count", "Atoms", "SMARTS"], "Functional groups", "Substructure search")],
        subgroup="molecular", live=True, example={"structure": _smiles(CAFFEINE)}),

    # ---- Detectability：7 页。只有 Platform 的脂质组 / 代谢组部分接了引擎，其余都是占位 ----
    module("platform_fit", "Platform", "A", "Which workflow can measure it: lipidomics, metabolomics, CIL LC-MS or GC-MS (decided mainly by extraction and derivatization).",
        [structure()],
        [cards([("Best platform", "Lipidome / Metabolome / Both / Neither")], "Lipidomics vs metabolomics", key="detect", id="detect"),
         figure("Structure with lipid features highlighted (indigo: acyl chain, teal: isoprenoid / sterol, periwinkle: polar head, amber: sugar)", key="detect", id="detect_structure"),
         table(["Feature", "Value"], key="detect", id="detect_features"),
         notes("Rule hits", key="detect", id="detect_reasons"),
         cards([("CIL LC-MS", "yes / no"), ("Labelable groups", "amine / phenol / carboxyl / carbonyl / hydroxyl")], "Chemical isotope labelling (CIL)", key="cil", id="cil"),
         cards([("GC-MS", "yes / after derivatization / no"), ("Volatility", ""), ("Thermal stability", "")], "GC-MS", key="gc", id="gc"),
         notes("CIL / GC-MS rule hits", key="gc", id="gc_cil_reasons")],
        subgroup="detectability", live=True, example={"structure": _smiles("CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCCCC")}),

    module("separation", "Separation", "B", "Which column, and what charge the molecule carries at each pH (MolGpKa pKa prediction; needs a one-time weights download).",
        [structure(), slider("ph", "Mobile phase pH", 1, 13, 0.1, 2.7)],
        [cards([("logD at pH", "logP + log10(fraction neutral)"), ("Dominant charge state", "net charge"), ("Fraction neutral", "fully uncharged")], "pH behaviour", key="ph", id="ph_summary"),
         table(["Site", "pKa", "Type"], "Ionizable sites", "Site: element and atom number (as in the SMILES, from 1). Acid: loses H, base: takes H. Base sites on amide, sulfonamide and pyrrole-type N are ignored (listed in Rule hits)", key="ph", id="sites"),
         figure("logD (left) and expected net charge (right) vs pH; dotted line = chosen pH", key="ph", id="ph_curves"),
         figure("Structure with ionizable sites highlighted (red: acid, teal: base, violet: both; label = atom and pKa)", key="ph", id="sites_structure"),
         table(["Mode", "Recommendation", "Reason"], "Separation mode", "Simple rules on logD, net charge, TPSA, MW and stereocentres; thresholds in mzdesk/platform_fit/separation_rules.py", key="column", id="modes"),
         notes("Rule hits", key="column", id="rule_hits")],
        subgroup="detectability", live=True, example={"structure": _smiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O"), "ph": 2.7}),

    module("ionization", "Ionization", "B", "Which ion source and polarity: ESI + / −, APCI, APPI for LC; EI or CI for GC.",
        [structure(), radio("platform", "Platform", ["LC-MS", "GC-MS"])],
        [cards([("Best source", "ESI / APCI / APPI / EI / CI"), ("Polarity", "+ / − / both")], "Recommendation"),
         table(["Source", "Polarity", "Suitability", "Reason"], "All sources"),
         notes("Rule hits")],
        subgroup="detectability"),

    module("adducts", "Adducts", "B", "Which adducts appear for the chosen polarity and mobile phase additive, and which is strongest.",
        [structure(), radio("ion_mode", "Ion mode", ["Positive", "Negative"]),
         select("additive", "Mobile phase additive", ["Formic acid", "Ammonium formate", "Ammonium acetate", "Acetic acid", "None"])],
        [cards([("Dominant adduct", ""), ("m/z", "")], "Expected base ion"),
         table(["Adduct", "Ion formula", "m/z", "Likelihood", "Reason"], "Adducts", "Ranked by expected intensity"),
         notes("Rule hits")],
        subgroup="detectability"),

    module("non_ms", "Non-MS detection", "B", "Can it be measured without MS: UV chromophore, fluorescence, or universal detectors (CAD / ELSD / RI).",
        [structure()],
        [cards([("UV chromophore", "yes / no"), ("Estimated λmax", "nm"), ("Native fluorescence", "likely / unlikely")], "Optical detection"),
         table(["Detector", "Suitability", "Reason"], "Detectors", "UV / PDA / FLD / CAD / ELSD / RI"),
         figure("Structure with chromophore highlighted")],
        subgroup="detectability"),

    module("derivatization", "Derivatization", "B", "Reactive groups and the derivatized molecule for GC-MS, CIL LC-MS or LC-UV / FLD.",
        [structure(),
         multiselect("reagents", "Reagents",
                     ["GC · TMS (MSTFA / BSTFA)", "GC · MeOx + TMS", "GC · TBDMS (MTBSTFA)", "GC · Acylation (TFAA)",
                      "CIL · Dansyl chloride (amine / phenol)", "CIL · Dansyl chloride, base-activated (hydroxyl)",
                      "CIL · DmPA bromide (carboxyl)", "CIL · Dansylhydrazine (carbonyl)",
                      "UV/FLD · OPA + 2-mercaptoethanol (primary amine)", "UV/FLD · FMOC-Cl (amine)", "UV/FLD · DNPH (carbonyl)"],
                     ["GC · TMS (MSTFA / BSTFA)", "GC · MeOx + TMS", "CIL · Dansyl chloride (amine / phenol)",
                      "CIL · Dansyl chloride, base-activated (hydroxyl)", "CIL · DmPA bromide (carboxyl)", "CIL · Dansylhydrazine (carbonyl)"])],
        [table(["Group", "Count", "Atoms"], "Reactive groups", id="deriv_groups"),
         table(["Reagent", "Derivatized groups", "Labels", "Formula", "Monoisotopic mass", "Heavy (13C) mass", "Notes"], "Derivatives",
               "Full derivatization; one row per reagent", id="deriv_rows"),
         figure("Structure with reactive groups highlighted", id="deriv_structure")],
        subgroup="detectability", live=True, example={"structure": _smiles("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O")}),

    module("known_methods", "Known methods", "B", "Is there a published method (pharmacopoeia, RepoRT)? If not, what do the most similar compounds use?",
        [structure(), multiselect("sources", "Sources", ["Japanese Pharmacopoeia", "RepoRT"])],
        [cards([("Match", "exact / similar / none"), ("Methods found", ""), ("Sources searched", "")], "Summary", id="summary"),
         table(["Source", "Compound", "Match", "Similarity", "Test / dataset", "Column", "Mobile phase", "Detector", "RT (min)", "Reference"], "Methods",
               "Exact hits first, then the methods of similar compounds (similarity score: see Why these compounds are similar). RepoRT is LC-MS; its RT is that compound's retention time in the dataset", id="methods"),
         table(["Block", "Test", "Condition text"], "JP condition text", "Printed JP18 conditions of the JP methods above (cross-referenced blocks show the referenced text too)", id="jp_conditions"),
         notes("Why these compounds are similar", id="similar_notes"),
         notes("Sources & licence", id="sources_notes")],
        subgroup="detectability", live=True, example={"structure": _smiles("C[C@H](CS)C(=O)N1CCC[C@H]1C(=O)O"), "sources": ["Japanese Pharmacopoeia", "RepoRT"]}),

    # ---- MS1 ----
    module("ms1_theory", "Theoretical m/z", "A", "Precursor m/z of a structure or formula for an adduct, and its isotope pattern.",
        [structure(label="Compound or formula"), precursor(label="Observed precursor m/z (optional)")],
        [cards([("Precursor m/z", "theoretical"), ("Ion formula", ""), ("Error vs observed", "ppm")], id="theory"),
         table(["Adduct", "Ion formula", "Theoretical m/z", "Δ ppm vs observed"], "All adducts", "Δ ppm shown when the observed m/z is within 100 ppm", id="adducts"),
         table(["Isotope", "Theoretical m/z", "Theoretical abundance"], "Isotope pattern"), figure("Theoretical isotope pattern")],
        subgroup="ms1", live=True, example={"structure": _smiles(CAFFEINE), "precursor": _prec("195.0877")}),

    module("ms1_match", "Isotope match", "A", "Score the observed isotope pattern against theory; flag Cl / Br / S.",
        [structure(label="Compound or formula"), precursor(), ms1(), number("tol_ppm", "Mass tolerance (ppm)", 5.0, 0.1)],
        [cards([("Isotope score", "0–1"), ("Cl / Br", "present?"), ("S", "present?")]),
         table(["Isotope", "Theoretical m/z", "Observed m/z", "Δ ppm", "Theoretical %", "Observed %"]),
         figure("Theoretical vs observed isotope pattern")],
        subgroup="ms1"),

    module("ms1_search", "Exact-mass search", "A", "Which database compounds match the precursor m/z.",
        [precursor(), number("tol_ppm", "Mass tolerance (ppm)", 5.0, 0.1), multiselect("databases", "Exact-mass databases", DATABASES, ["HMDB", "PubChem"])],
        [table(["Formula", "Database", "ID", "Name", "Δ ppm"], "Exact-mass hits")],
        subgroup="ms1"),

    module("ms1_formula", "Formula (BUDDY)", "A", "Rank molecular formulas for the precursor from its MS1 isotopes and MS2 (BUDDY).",
        [precursor(adducts=BUDDY_ADDUCTS), ms1(), spectrum(label="MS2 spectrum (optional)"),
         select("instrument", "Instrument", ["Orbitrap", "Q-TOF", "FT-ICR"]), checkbox("halogen", "Include halogens (F, Cl, Br, I)"),
         text("formula", "Formula (optional: skip ranking, explain the MS2 with it)", placeholder="C8H10N4O2")],
        [cards([("Top formula", ""), ("Adduct used", ""), ("MS1 isotope peaks used", "")], id="summary"),
         table(["Rank", "Formula", "Theoretical m/z", "Δ ppm", "Probability", "FDR", "MS1 isotope similarity"], "Formula candidates", id="candidates"),
         figure("MS2 fragment formulas of the top candidate (hover for details)", id="fragments_plot"),
         table(["m/z", "Rel. intensity (%)", "Fragment formula", "Δ ppm", "Neutral loss", "Other fits"], id="fragments")],
        subgroup="ms1", live=True,
        example={"precursor": _prec("195.0877"), "ms1": CAFFEINE_MS1, "ms2": _paste(CAFFEINE_MS2)}),

    # ---- Spectrum（subgroup id 仍是 ms2）----
    module("ms2_clean", "Clean & Vis", "B", "Clean an MS2 spectrum, see it, and get it in every common format (MS-DIAL, GNPS, MassBank, HMDB, MSP, MGF …).",
        [spectrum("spectrum", "MS2 spectrum (any common format)"), number("precursor_mz", "Precursor m/z (optional)", step=0.0001),
         radio("ion_mode", "Ion mode (for MSP / MGF headers)", ["Positive", "Negative"]),
         number("noise_pct", "Noise threshold (% of base peak)", 0, 0.1), select("normalize", "Normalize", NORMALIZE),
         checkbox("remove_precursor", "Remove precursor"), number("precursor_tol", "Precursor tolerance (Da)", 0.02, 0.001),
         select("sort_by", "Sort by", SORT_BY)],
        [cards([("Detected format", ""), ("Peaks in", ""), ("Peaks out", ""), ("Removed", "noise · precursor")], id="summary"),
         figure("Cleaned spectrum; light grey = removed (hover a peak for details)", id="spectrum"),
         table(["m/z", "Intensity", "Rel. intensity (%)", "Kept", "Reason"], "Peaks", "Rel. intensity: % of the input base peak. Reason: noise, precursor", id="peaks"),
         formats("Cleaned spectrum in every format", "Formats", id="formats")],
        subgroup="ms2", live=True,
        example={"spectrum": _paste("42.0338:5 56.0495:4 61.0012:0.3 69.0447:6 83.0604:9 97.5521:0.6 110.0713:18 138.0662:100 195.0877:32"),
                 "precursor_mz": "195.0877", "noise_pct": "1", "normalize": "Base peak = 100", "remove_precursor": True}),

    module("ms2_ei", "EI MS", "A", "Predict the 70 eV electron-ionization spectrum of a structure (NEIMS). One ion type: [M]+•.",
        [structure()],
        [cards([("Ion", ""), ("Formula", ""), ("Molecular ion m/z", ""), ("Base peak m/z", "")], id="summary"),
         figure("Predicted EI spectrum, unit mass (hover a peak for details)", id="spectrum"),
         table(["m/z", "Rel. intensity (%)"], "Peaks", id="peaks"),
         formats("Predicted spectrum in every format", "Formats", id="formats")],
        subgroup="ms2", live=True, example={"structure": _smiles(CAFFEINE)}),

    module("ms2_esi", "ESI MS/MS", "A", "Predict the ESI MS/MS spectrum of a structure at 10 / 20 / 40 eV (CFM-ID 4; MSRB rules for lipids). Needs Docker.",
        [structure(), select("adduct", "Adduct", ESI_ADDUCTS, "[M+H]+")],
        [cards([("Engine", ""), ("Precursor m/z", ""), ("Ion formula", ""), ("Peaks", "")], id="summary"),
         figure("Low energy (10 eV)", id="spectrum_low"), figure("Medium energy (20 eV)", id="spectrum_medium"),
         figure("High energy (40 eV)", id="spectrum_high"),
         table(["Energy", "m/z", "Rel. intensity (%)", "Fragment"], "Peaks", "Fragment: formula · SMILES from the engine", id="peaks"),
         formats("Predicted spectra (one record per energy)", "Formats", id="formats")],
        subgroup="ms2", live=True, example={"structure": _smiles(CAFFEINE), "adduct": "[M+H]+"}),

    module("ms2_library", "Library search", "A", "Find the spectrum in MassBank / MoNA / GNPS and check the best hit.",
        [spectrum(), precursor(), number("tol_da", "Fragment tolerance (Da)", 0.02, 0.001), multiselect("libraries", "Libraries", LIBRARIES)],
        [table(["Library", "Accession", "Name", "Score", "Matched peaks", "Precursor Δ", "Instrument / CE"], "Library hits"),
         figure("Mirror: observed vs best hit")],
        subgroup="ms2"),

    module("ms2_similarity", "Similarity & mirror", "A", "How alike are two spectra: 12 classic scores, optional deep-learning scores, mirror plot.",
        [spectrum("query", "Query spectrum"), number("query_precursor", "Query precursor m/z (optional)", step=0.0001),
         spectrum("reference", "Reference spectrum"), number("reference_precursor", "Reference precursor m/z (optional)", step=0.0001),
         number("tol_da", "Fragment tolerance (Da)", 0.02, 0.001), number("threshold", "Score threshold", 0.7, 0.05),
         radio("ion_mode", "Ion mode", ["Positive", "Negative"]),
         checkbox("deep", "Deep-learning scores (MS2DeepScore, Spec2Vec; first use downloads the models)")],
        [cards([("Primary score", "spectral entropy, weighted"), ("Metrics ≥ threshold", "n / N"), ("Matched peaks", "pairs"),
                ("Query matched", "%"), ("Reference matched", "%")], id="summary"),
         figure("Mirror: query (up) vs reference (down); dark = matched", id="mirror"),
         table(["Category", "Metric", "Score", "Matched peaks", "≥ threshold", "Note"], "Scores", id="scores")],
        subgroup="ms2", live=True,
        example={"query": _paste("197.04643 10\n255.22775 92\n256.23184 31\n277.60822 10\n283.24814 20"),
                 "reference": _paste("197.04643 8\n255.22775 100\n256.23184 28\n277.60822 7\n283.24814 22")}),

    module("ms2_fragments", "Fragments (MetFrag)", "A", "Which peaks the candidate structure explains, and with which fragment (MetFrag).",
        [structure(label="Candidate structure"), precursor(label="Precursor m/z (optional; blank = theoretical)", adducts=METFRAG_ADDUCTS), spectrum(),
         slider("tree_depth", "Tree depth (bond cleavages)", 1, 4, 1, 2),
         number("tol_abs", "Absolute tolerance (Da)", 0.002, 0.001), number("tol_ppm", "Relative tolerance (ppm)", 10.0, 1)],
        [cards([("Explained peaks", "n / N"), ("Explained intensity", "%"), ("Theoretical precursor m/z", ""), ("Ion formula", "")], id="summary"),
         figure("Observed spectrum (hover: fragment formula, Δ ppm, neutral loss, structure)", id="spectrum"),
         table(["Observed m/z", "Rel. intensity (%)", "Kind", "Ion formula", "Theoretical m/z", "Δ ppm", "Neutral loss", "Fragment structure"], id="peaks")],
        subgroup="ms2", live=True,
        example={"structure": _smiles(PG_36_2), "precursor": _prec("", "Negative", "[M-H]-"), "ms2": _paste(PG_36_2_MS2)}),

    module("ms2_denovo", "De novo", "A", "Structure candidates when no structure is given.",
        [precursor(), spectrum(), text("formula", "Formula (optional)")],
        [table(["Rank", "SMILES", "Formula", "Score", "Tool"], "De novo candidates")],
        subgroup="ms2"),

    module("exclusivity", "Exclusivity", "A", "Is there another explanation: same formula, isobaric interference, ion relations, isomers, metabolites, lipid acyl combinations.",
        [structure(label="Compound or formula", optional=True), precursor(), number("tol_ppm", "Tolerance (ppm)", 5.0, 0.1),
         number("mz2", "Second m/z (relation check, optional)", step=0.0001), multiselect("databases", "Databases", DATABASES, ["HMDB", "PubChem"]),
         spectrum(label="MS2 spectrum (optional, for isomer check)")],
        [table(["InChIKey (first block)", "Name", "Databases", "Stereoisomers collapsed"], "Same-formula candidates"),
         table(["Formula", "Name", "Database", "Exact m/z", "Δ ppm"], "Isobaric interferences", "Different formula, same m/z within tolerance"),
         table(["Relation", "Δ mass", "Explanation", "Confidence"], "Relation between two m/z", "Adduct / in-source fragment / dimer / neutral loss"),
         table(["Isomer SMILES", "Distinguishable by MS2", "Distinguishing fragments"], "Isomers"),
         table(["Direction", "Parent", "Product", "Reaction", "Tool"], "Metabolite / transformation products", "BioTransformer, SyGMa"),
         table(["Combination", "Diagnostic ions", "Observed"], "Lipids: sum composition → acyl combinations")]),

    # ---- Orthogonal ----
    module("rt", "RT / RI", "A", "Does the retention agree: RT / RI prediction and projection across systems.",
        [structure(), select("system", "Chromatographic system", ["RP C18 generic", "HILIC generic", "GC (Kovats)", "Custom (upload calibration)"]),
         number("rt_min", "Observed RT (min, optional)", step=0.01), file("calibration", "Calibration set for projection (optional, CSV: compound, RT)", ".csv")],
        [cards([("Predicted RT / RI", ""), ("Observed", ""), ("Δ", "")], "Retention"), table(["System", "Predicted", "Projected", "Observed", "Δ"])],
        subgroup="orthogonal"),

    module("ccs", "CCS", "A", "Does the ion mobility agree: CCS prediction and lookup.",
        [structure(), select("adduct", "Adduct", ADDUCTS_POS[1:] + ADDUCTS_NEG[1:]), number("ccs", "Observed CCS (Å², optional)", step=0.1)],
        [cards([("Predicted CCS", "Å²"), ("Database CCS", "Å²"), ("Observed", "Å²"), ("Δ", "%")], "CCS"), table(["Adduct", "Predicted", "Database", "Observed", "Δ %"])],
        subgroup="orthogonal"),

    module("mrm", "MRM / targeted transitions", "B", "Pick precursor–product pairs and collision energies for a targeted method.",
        [radio("start_from", "Start from", ["Structure (predicted fragments)", "Library spectrum", "My spectrum"]), structure(optional=True), spectrum(label="MS2 spectrum (if starting from my spectrum)"),
         radio("ion_mode", "Ion mode", ["Positive", "Negative"]), select("instrument", "Instrument", ["Triple quadrupole (generic)", "QTRAP", "Q-TOF (PRM)", "Orbitrap (PRM)"]),
         number("n_transitions", "Transitions", 3, 1), select("ce_unit", "CE unit", ["eV", "NCE"])],
        [cards([("Predicted RT", "min"), ("Polarity", ""), ("Transitions", "n")]),
         table(["Precursor m/z", "Product m/z", "CE", "Role", "Fragment formula", "Source", "Relative intensity"], caption="Role: quantifier / qualifier. Source: library spectrum or predicted"),
         button("Export method table (CSV)")]),

    module("calculator", "Calculator", "B", "Everyday calculations. Own inputs, nothing shared.",
        [heading("ppm error"), number("ppm_theoretical", "Theoretical m/z", step=0.0001), number("ppm_observed", "Observed m/z", step=0.0001),
         heading("Isotope labelling"), text("lab_smiles", "SMILES"), text("lab_label", "Label", placeholder="13C6, d4, 15N2"), select("lab_adduct", "Adduct", ADDUCTS_POS[1:]),
         heading("Concentration"), number("conc_value", "Value"), select("conc_from", "From", ["µM", "mM", "ng/mL", "µg/mL", "mg/L"]), select("conc_to", "To", ["µM", "mM", "ng/mL", "µg/mL", "mg/L"]),
         heading("Collision energy"), number("ce_nce", "NCE"), number("ce_ev", "eV"), number("ce_mz", "Precursor m/z", step=0.0001), number("ce_z", "Charge", 1, 1),
         heading("GC retention index (Kovats)"), grid("alkanes", "Alkane series", ["Carbon number", "RT"]), number("ri_rt", "Analyte RT", step=0.01)],
        [cards([("ppm error", "ppm")], "ppm error"),
         table(["Isotopologue", "m/z", "Abundance"], "Isotope labelling"), table(["Internal standard", "Analyte M+n", "Crosstalk %"]),
         cards([("Converted", "")], "Concentration"), table(["Target conc.", "Stock conc.", "Final volume", "Stock volume", "Diluent volume"], caption="Dilution and stock preparation"),
         cards([("NCE ↔ eV", "")], "Collision energy"), cards([("Kovats RI", "")], "GC retention index")]),

    module("report", "Report", "A/B", "Collect the evidence, give a confidence level, export. Plus name normalization, figures, MassBank record.",
        [select("source", "Evidence from", ["Last pipeline run", "Upload report JSON"]), file("report_json", "Report JSON (optional)", ".json")],
        [cards([("Schymanski level", "1–5"), ("Evidence items", "n"), ("Conflicts", "n")], "Annotation report (A)"),
         table(["Module", "Evidence", "Result", "Supports", "Weight"], caption="One row per piece of evidence, in pipeline order"),
         table(["Parameter", "Value", "Module"], caption="Provenance"), table(["Tool", "Version", "Citation"], caption="What to cite"),
         button("Export JSON"), button("Export CSV"), button("Export PDF"),
         table(["Convention", "Name"], "Name normalization (B)", "Goslin lipid shorthand, LIPID MAPS, RefMet"),
         button("Spectrum with structures (SVG)"), button("Mirror plot (SVG)"),
         code("MassBank record", "MassBank record generator (B)")]),

    # ---- DMPK（非核心）----
    module("met_shift", "Metabolite mass-shift search", "B", "Search a peak or feature list for phase I / II metabolites of the parent.",
        [structure(label="Parent compound or formula"), precursor(),
         grid("features", "Peak / feature list", ["m/z", "RT", "Intensity"], "Paste or upload: one feature per line"), file("features_file", "Feature table (CSV, optional)", ".csv"),
         multiselect("shifts", "Mass shifts", ["+O (hydroxylation)", "−2H (dehydrogenation)", "+2H (reduction)", "−CH2 (demethylation)", "+O −CH2", "+C6H8O6 glucuronide", "+SO3 sulfate", "+C2H2O acetyl", "+C10H15N3O6S GSH", "+C2H3NO2 glycine"]),
         number("tol_ppm", "Tolerance (ppm)", 10.0, 0.1)],
        [table(["Shift", "Δ mass", "Expected m/z", "Found m/z", "Δ ppm", "RT", "Intensity", "Candidate reaction"]),
         figure("Features with metabolite hits marked (m/z vs RT)")], subgroup="dmpk"),

    module("som", "Site of metabolism", "B", "Which atoms are metabolized first.",
        [structure(), multiselect("enzymes", "Enzymes", ["CYP3A4", "CYP2D6", "CYP2C9", "CYP2C19", "CYP1A2", "UGT", "All"], ["All"])],
        [table(["Atom", "Rank", "Score", "Reaction type", "Enzyme", "Tool"], caption="SMARTCyp / FAME"), figure("Structure with predicted sites highlighted")], subgroup="dmpk"),

    module("stability", "Metabolic stability", "B", "t½ and intrinsic clearance from an in vitro incubation.",
        [grid("timecourse", "Time course", ["Time (min)", "% remaining"]), select("system", "System", ["Liver microsomes", "Hepatocytes", "S9", "Plasma"]),
         number("protein", "Protein (mg/mL) or cells (10⁶/mL)", 0.5, 0.1), select("species", "Species", ["Human", "Rat", "Mouse", "Dog", "Monkey"])],
        [cards([("t½", "min"), ("CLint", "µL/min/mg or µL/min/10⁶ cells"), ("Scaled CLint", "mL/min/kg"), ("Predicted hepatic CL", "mL/min/kg, well-stirred"), ("Extraction ratio", "")]),
         figure("ln(% remaining) vs time with fit"), table(["Parameter", "Value", "Unit", "Assumption"], caption="Scaling factors used")], subgroup="dmpk"),

    module("nca", "NCA pharmacokinetics", "B", "Non-compartmental parameters from a concentration–time profile.",
        [grid("profile", "Concentration–time", ["Time (h)", "Concentration (ng/mL)"]), number("dose", "Dose (mg)", step=0.01),
         select("route", "Route", ["IV bolus", "IV infusion", "Oral"]), number("lambda_points", "Points for λz", 3, 1)],
        [cards([("Cmax", "ng/mL"), ("Tmax", "h"), ("AUC0-t", "ng·h/mL"), ("AUC0-∞", "ng·h/mL"), ("t½", "h"), ("CL/F", "L/h"), ("Vz/F", "L"), ("MRT", "h")]),
         figure("Concentration–time (semi-log) with λz fit"), table(["Parameter", "Value", "Unit", "Method"])], subgroup="dmpk"),

    module("sil_is", "Isotope-labelled internal standard", "B", "Which labelling to buy: mass shift, stable positions, crosstalk.",
        [structure(label="Analyte"), spectrum(label="Analyte MS2 (optional)"), number("min_shift", "Minimum mass shift (Da)", 3, 1), checkbox("check_vendors", "Check vendor availability", True)],
        [table(["Labelling", "Mass shift", "Positions", "Exchangeable-H risk", "Metabolic-loss risk", "Label kept in quantifier fragment", "Crosstalk with analyte M+n (%)", "Vendor"]),
         figure("Structure with labelled positions highlighted")], subgroup="dmpk"),
]


# ---- Pipelines ----
def stage(title, modules, exclude=None):
    """modules: 一个或多个模块 id，按顺序显示它们的输出块；exclude 按块的 key 过滤。"""
    return {"title": title, "modules": [modules] if isinstance(modules, str) else list(modules), "exclude": exclude or []}


PIPELINES = [
    {"id": "library_search", "title": "Library search", "one_liner": "I have a spectrum and no structure. Find it in the libraries and check the hit.",
     "inputs": [spectrum(), precursor(), multiselect("libraries", "Libraries", LIBRARIES), number("tol_ppm", "Tolerance (ppm)", 10.0, 0.1)],
     "stages": [stage("Spectrum cleanup", "ms2_clean"), stage("MS2 library search", ["ms2_library", "ms2_similarity"]),
                stage("MS1 check of the hits", ["ms1_theory", "ms1_match"]), stage("Exclusivity", "exclusivity"), stage("Report", "report")]},
    {"id": "de_novo", "title": "De novo annotation", "one_liner": "I have a spectrum, the libraries do not have it. Infer formula, then structure candidates.",
     "inputs": [precursor(), ms1(), spectrum(), number("rt_min", "RT (min, optional)", step=0.01), number("ccs", "CCS (Å², optional)", step=0.1),
                text("elements", "Element constraints", "C0-50 H0-100 N0-10 O0-20 P0-3 S0-3"),
                multiselect("databases", "Candidate databases", DATABASES, ["HMDB", "PubChem"])],
     "stages": [stage("Spectrum cleanup", "ms2_clean"), stage("MS1: formula", "ms1_formula"), stage("MS2: fragment explanation / de novo", ["ms2_fragments", "ms2_denovo"]),
                stage("Exclusivity", "exclusivity"), stage("Orthogonal", ["rt", "ccs"]), stage("Report", "report")]},
    {"id": "verify", "title": "Result verification", "one_liner": "I have a candidate structure. Does every piece of evidence support it?",
     "inputs": [structure(label="Candidate structure"), precursor(), ms1(), spectrum(),
                number("rt_min", "RT (min, optional)", step=0.01), number("ccs", "CCS (Å², optional)", step=0.1)],
     "stages": [stage("Molecular", ["structure_formats", "properties"]), stage("Platform fit", "platform_fit"), stage("MS1", ["ms1_theory", "ms1_match"]),
                stage("MS2", ["ms2_esi", "ms2_library", "ms2_similarity", "ms2_fragments"]),
                stage("Exclusivity", "exclusivity"), stage("Orthogonal", ["rt", "ccs"]), stage("Report", "report")]},
    {"id": "feature_prediction", "title": "MS feature prediction", "one_liner": "I have a structure and no spectrum. Predict what the instrument will see: adducts, precursor m/z, isotope pattern, MS2.",
     "inputs": [structure(), radio("ion_mode", "Ion mode", ["Positive", "Negative", "Both"]), multiselect("adducts", "Adducts", ["Auto"] + ADDUCTS_POS[1:] + ADDUCTS_NEG[1:], ["Auto"]),
                select("instrument", "Instrument / CE", ["Q-TOF 20 eV", "Q-TOF 40 eV", "Orbitrap HCD NCE 30", "Triple quadrupole"]), select("out_format", "Export format", ["MSP", "MGF", "MassBank record", "JSON"])],
     "stages": [stage("Molecular", ["structure_formats", "properties"]), stage("Adduct prediction", "adducts"),
                stage("MS1 theory: precursor m/z and isotope pattern", "ms1_theory"),
                stage("Spectrum prediction (ESI MS/MS, EI MS)", ["ms2_esi", "ms2_ei"]),
                stage("Export", "ms2_clean")]},
]

GROUPS = [{"id": "pipelines", "title": "Pipelines"}, {"id": "modules", "title": "Modules"}]
# Modules 下的小类，侧栏在该小类第一个模块的位置展开
SUBGROUPS = [{"id": "molecular", "title": "Molecular"}, {"id": "detectability", "title": "Detectability"},
             {"id": "ms1", "title": "MS1"}, {"id": "ms2", "title": "Spectrum"}, {"id": "orthogonal", "title": "Orthogonal"}, {"id": "dmpk", "title": "DMPK"}]


def spec():
    return {"groups": GROUPS, "subgroups": SUBGROUPS, "pipelines": PIPELINES, "modules": MODULES}
