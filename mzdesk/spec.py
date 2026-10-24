"""页面定义（唯一来源）。前端按这里的 pipelines / modules 渲染输入和输出格式。

字段类型：structure, precursor, spectrum, ms1, number, text, select, multiselect, radio,
checkbox, slider, file, table, heading。
输出块类型：table, cards, figure, code, button, heading。块可带 key，pipeline 的 stage 用 exclude 过滤。
"""
ADDUCTS_POS = ["Auto", "[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M]+", "[M+H-H2O]+", "[2M+H]+"]
ADDUCTS_NEG = ["Auto", "[M-H]-", "[M+Cl]-", "[M+HCOO]-", "[M+CH3COO]-", "[M-H-H2O]-", "[2M-H]-"]
LIBRARIES = ["MassBank", "MoNA", "GNPS"]
DATABASES = ["HMDB", "PubChem", "LIPID MAPS", "KEGG", "ChEBI", "DrugBank", "NORMAN SusDat"]


# ---- 字段 ----
def structure(id="structure", label="Compound", optional=False):
    return {"type": "structure", "id": id, "label": label, "optional": optional}
def precursor(id="precursor", adduct=True):
    return {"type": "precursor", "id": id, "adduct": adduct}
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
def rt_ccs():
    return [number("rt_min", "RT (min, optional)", step=0.01), number("ccs", "CCS (Å², optional)", step=0.1)]


# ---- 输出块 ----
def table(columns, title=None, caption=None, key=None):
    return {"type": "table", "title": title, "columns": columns, "caption": caption, "key": key}
def cards(items, title=None, key=None):
    return {"type": "cards", "title": title, "items": [{"label": l, "unit": u} for l, u in items], "key": key}
def figure(label, key=None):
    return {"type": "figure", "label": label, "key": key}
def code(label, title=None, key=None):
    return {"type": "code", "title": title, "label": label, "key": key}
def button(label, key=None):
    return {"type": "button", "label": label, "key": key}
def head(title, key=None):
    return {"type": "heading", "title": title, "key": key}


# ---- 模块 ----
def module(id, title, cls, one_liner, inputs, outputs, group="modules"):
    return {"id": id, "title": title, "cls": cls, "group": group, "one_liner": one_liner, "inputs": inputs, "outputs": outputs}


MODULES = [
    module("identity", "Identity & background", "B", "Who is this compound, where does it come from, what is known about it.",
        [structure()],
        [head("Identifiers"),
         table(["Name", "IUPAC", "CAS", "SMILES", "InChI", "InChIKey", "PubChem CID", "HMDB", "KEGG", "ChEBI", "LIPID MAPS", "DrugBank"]),
         table(["Database", "ID", "Link", "Synonyms"], "Cross-links & synonyms"),
         table(["System", "Level", "Term"], "Classification", "ClassyFire, NPClassifier, LIPID MAPS"),
         table(["Origin", "Evidence", "Source"], "Origin", "Endogenous / drug / food / pollutant / microbial / plant"),
         table(["Pathway", "Database", "Link"], "Biology"), table(["Biofluid", "Normal concentration", "Unit", "Source"]), table(["Species", "Source (LOTUS)"]),
         table(["Vendor", "CAS", "Purity", "GHS hazard"], "Standard availability")]),

    module("properties", "Molecular properties", "B", "Representation, descriptors, functional groups, structure formats.",
        [structure()],
        [figure("2D depiction"), figure("3D viewer"),
         cards([("Formula", ""), ("Exact mass", "Da"), ("logP", ""), ("Charge", "")], "Properties"),
         table(["Descriptor", "Value"], caption="MW, HBD, HBA, TPSA, rotatable bonds, rings"),
         table(["Group", "Count", "Atoms", "SMARTS"], "Functional groups (substructure search)"),
         table(["Format", "Value"], "Structure format conversion", "SMILES, MOL, InChI, InChIKey — each copyable")]),

    module("platform_fit", "Platform fit", "A", "Can my platform see it: pH behaviour, column choice, ionization.",
        [structure(), slider("ph", "pH", 1, 13, 0.1, 7.4), radio("ion_mode", "Ion mode", ["Positive", "Negative"])],
        [cards([("logD at pH", ""), ("Dominant charge state", "")], "pH behaviour", key="ph"),
         table(["Site", "pKa", "Type"], key="ph"), figure("logD vs pH", key="ph"),
         table(["Column", "Recommendation", "Reason"], "Column", "RP / HILIC", key="column"),
         cards([("ESI preference", "+ / − / both"), ("GC volatile", "yes / no / after derivatization"), ("Detectability", "Lipidomics / Metabolomics / Both / Neither")], "Ionization", key="ion"),
         table(["Adduct", "m/z", "Likelihood"], caption="Dominant adducts for the chosen ion mode", key="ion")]),

    module("gc_deriv", "GC derivatization", "B", "Mass and formula of the derivatized molecule.",
        [structure(), multiselect("reagents", "Reagents", ["TMS (MSTFA / BSTFA)", "MeOx + TMS", "TBDMS (MTBSTFA)", "Acylation (TFAA / PFPA)"], ["TMS (MSTFA / BSTFA)", "MeOx + TMS"])],
        [table(["Reagent", "Derivative", "Formula", "MW", "Derivatized groups", "Monoisotopic mass", "Key EI ions"], caption="One row per reagent")]),

    module("cleanup", "Spectrum cleanup & conversion", "B", "Clean one or more spectra, merge collision energies, convert format. First step of every pipeline.",
        [spectrum("spectra", "Spectrum", multi=True), precursor(adduct=False),
         number("noise_pct", "Noise threshold (% base peak)", 1.0, 0.1), checkbox("deisotope", "Deisotope", True),
         checkbox("remove_precursor", "Remove precursor", True), checkbox("merge_ce", "Merge collision energies"),
         checkbox("normalize", "Normalize to 100", True), select("out_format", "Output format", ["MSP", "MGF", "MassBank record", "JSON"])],
        [cards([("Peaks in", ""), ("Peaks out", ""), ("Removed", "noise / isotope / precursor")]),
         table(["m/z", "Intensity", "Kept", "Reason"], caption="Reason: noise, isotope peak, precursor, merged from CE"),
         figure("Before vs after (mirror)"), code("converted spectrum", "Converted record"), button("Download converted file")]),

    module("ms1", "MS1", "A", "Does the precursor and its isotope pattern fit: theory, match, exact-mass search, formula inference.",
        [structure(label="Compound or formula", optional=True), precursor(), ms1(),
         number("tol_ppm", "Mass tolerance (ppm)", 5.0, 0.1), text("elements", "Element constraints", "C0-50 H0-100 N0-10 O0-20 P0-3 S0-3"),
         multiselect("databases", "Exact-mass databases", DATABASES, ["HMDB", "PubChem"])],
        [cards([("Precursor m/z", "from structure + adduct"), ("Formula", ""), ("Error vs observed", "ppm")], "Theoretical", key="theory"),
         table(["Isotope", "Theoretical m/z", "Theoretical abundance"], key="theory"), figure("Theoretical vs observed isotope pattern", key="theory"),
         cards([("Isotope score", "0–1"), ("Cl / Br", "present?"), ("S", "present?")], "Isotope match", key="match"),
         table(["Isotope", "Theoretical m/z", "Observed m/z", "Δ ppm", "Theoretical %", "Observed %"], key="match"),
         table(["Formula", "Database", "ID", "Name", "Δ ppm"], "Exact-mass hits", key="hits"),
         table(["Rank", "Formula", "Score", "FDR", "RDBE", "Element check", "Δ ppm"], "Formula candidates (BUDDY)", key="formula")]),

    module("ms2", "MS2", "A", "Does the fragment spectrum support the structure: predicted spectrum, library search, similarity, fragment explanation, de novo.",
        [structure(optional=True), precursor(), spectrum(),
         number("tol_da", "Fragment tolerance (Da)", 0.02, 0.001), text("ce", "Collision energy", placeholder="20 eV or NCE 30"),
         multiselect("libraries", "Libraries", LIBRARIES),
         multiselect("metrics", "Similarity metrics", ["Cosine", "Modified cosine", "Spectral entropy", "MS2DeepScore", "Spec2Vec"], ["Cosine", "Spectral entropy"])],
        [figure("Observed spectrum (hover: fragment formula, Δ ppm, neutral loss, structure)", key="observed"),
         table(["m/z", "Fragment formula", "Relative intensity", "Structure"], "Predicted spectrum (from structure)", key="predicted"), figure("Mirror: observed vs predicted", key="predicted"),
         table(["Library", "Accession", "Name", "Score", "Matched peaks", "Precursor Δ", "Instrument / CE"], "Library search", key="library"), figure("Mirror: observed vs best hit", key="library"),
         table(["Metric", "Score", "Threshold", "Pass"], "Similarity", key="similarity"),
         cards([("Explained peaks", "n / N"), ("Explained intensity", "%")], "Fragment explanation (MetFrag)", key="metfrag"),
         table(["Observed m/z", "Fragment formula", "Theoretical m/z", "Δ ppm", "Neutral loss", "Fragment structure"], key="metfrag"),
         table(["Rank", "SMILES", "Formula", "Score", "Tool"], "De novo (no structure given)", key="denovo")]),

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

    module("orthogonal", "Orthogonal", "A", "Do retention and mobility agree: RT / RI prediction and projection, CCS prediction and lookup.",
        [structure(), select("adduct", "Adduct (for CCS)", ADDUCTS_POS[1:] + ADDUCTS_NEG[1:]),
         select("system", "Chromatographic system", ["RP C18 generic", "HILIC generic", "GC (Kovats)", "Custom (upload calibration)"]),
         *rt_ccs(), file("calibration", "Calibration set for projection (optional, CSV: compound, RT)", ".csv")],
        [cards([("Predicted RT / RI", ""), ("Observed", ""), ("Δ", "")], "Retention"), table(["System", "Predicted", "Projected", "Observed", "Δ"]),
         cards([("Predicted CCS", "Å²"), ("Database CCS", "Å²"), ("Observed", "Å²"), ("Δ", "%")], "CCS"), table(["Adduct", "Predicted", "Database", "Observed", "Δ %"])]),

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
         figure("Features with metabolite hits marked (m/z vs RT)")], group="dmpk"),

    module("som", "Site of metabolism", "B", "Which atoms are metabolized first.",
        [structure(), multiselect("enzymes", "Enzymes", ["CYP3A4", "CYP2D6", "CYP2C9", "CYP2C19", "CYP1A2", "UGT", "All"], ["All"])],
        [table(["Atom", "Rank", "Score", "Reaction type", "Enzyme", "Tool"], caption="SMARTCyp / FAME"), figure("Structure with predicted sites highlighted")], group="dmpk"),

    module("stability", "Metabolic stability", "B", "t½ and intrinsic clearance from an in vitro incubation.",
        [grid("timecourse", "Time course", ["Time (min)", "% remaining"]), select("system", "System", ["Liver microsomes", "Hepatocytes", "S9", "Plasma"]),
         number("protein", "Protein (mg/mL) or cells (10⁶/mL)", 0.5, 0.1), select("species", "Species", ["Human", "Rat", "Mouse", "Dog", "Monkey"])],
        [cards([("t½", "min"), ("CLint", "µL/min/mg or µL/min/10⁶ cells"), ("Scaled CLint", "mL/min/kg"), ("Predicted hepatic CL", "mL/min/kg, well-stirred"), ("Extraction ratio", "")]),
         figure("ln(% remaining) vs time with fit"), table(["Parameter", "Value", "Unit", "Assumption"], caption="Scaling factors used")], group="dmpk"),

    module("nca", "NCA pharmacokinetics", "B", "Non-compartmental parameters from a concentration–time profile.",
        [grid("profile", "Concentration–time", ["Time (h)", "Concentration (ng/mL)"]), number("dose", "Dose (mg)", step=0.01),
         select("route", "Route", ["IV bolus", "IV infusion", "Oral"]), number("lambda_points", "Points for λz", 3, 1)],
        [cards([("Cmax", "ng/mL"), ("Tmax", "h"), ("AUC0-t", "ng·h/mL"), ("AUC0-∞", "ng·h/mL"), ("t½", "h"), ("CL/F", "L/h"), ("Vz/F", "L"), ("MRT", "h")]),
         figure("Concentration–time (semi-log) with λz fit"), table(["Parameter", "Value", "Unit", "Method"])], group="dmpk"),

    module("sil_is", "Isotope-labelled internal standard", "B", "Which labelling to buy: mass shift, stable positions, crosstalk.",
        [structure(label="Analyte"), spectrum(label="Analyte MS2 (optional)"), number("min_shift", "Minimum mass shift (Da)", 3, 1), checkbox("check_vendors", "Check vendor availability", True)],
        [table(["Labelling", "Mass shift", "Positions", "Exchangeable-H risk", "Metabolic-loss risk", "Label kept in quantifier fragment", "Crosstalk with analyte M+n (%)", "Vendor"]),
         figure("Structure with labelled positions highlighted")], group="dmpk"),
]


# ---- Pipelines ----
def stage(title, module_id, exclude=None):
    return {"title": title, "module": module_id, "exclude": exclude or []}


PIPELINES = [
    {"id": "library_search", "title": "Library search", "one_liner": "I have a spectrum and no structure. Find it in the libraries and check the hit.",
     "inputs": [spectrum(), precursor(), multiselect("libraries", "Libraries", LIBRARIES), number("tol_ppm", "Tolerance (ppm)", 10.0, 0.1)],
     "stages": [stage("Spectrum cleanup", "cleanup"), stage("MS2 library search", "ms2", ["predicted", "metfrag", "denovo"]),
                stage("MS1 check of the hits", "ms1", ["formula"]), stage("Exclusivity", "exclusivity"), stage("Report", "report")]},
    {"id": "de_novo", "title": "De novo annotation", "one_liner": "I have a spectrum, the libraries do not have it. Infer formula, then structure candidates.",
     "inputs": [precursor(), ms1(), spectrum(), *rt_ccs(), text("elements", "Element constraints", "C0-50 H0-100 N0-10 O0-20 P0-3 S0-3"),
                multiselect("databases", "Candidate databases", DATABASES, ["HMDB", "PubChem"])],
     "stages": [stage("Spectrum cleanup", "cleanup"), stage("MS1: formula", "ms1", ["theory"]), stage("MS2: fragment explanation / de novo", "ms2", ["predicted", "library"]),
                stage("Exclusivity", "exclusivity"), stage("Orthogonal", "orthogonal"), stage("Report", "report")]},
    {"id": "verify", "title": "Result verification", "one_liner": "I have a candidate structure. Does every piece of evidence support it?",
     "inputs": [structure(label="Candidate structure"), precursor(), ms1(), spectrum(), *rt_ccs()],
     "stages": [stage("Molecular properties", "properties"), stage("Platform fit", "platform_fit"), stage("MS1", "ms1", ["formula"]), stage("MS2", "ms2", ["denovo"]),
                stage("Exclusivity", "exclusivity"), stage("Orthogonal", "orthogonal"), stage("Report", "report")]},
    {"id": "feature_prediction", "title": "MS feature prediction", "one_liner": "I have a structure and no spectrum. Predict what the instrument will see: adducts, precursor m/z, isotope pattern, MS2.",
     "inputs": [structure(), radio("ion_mode", "Ion mode", ["Positive", "Negative", "Both"]), multiselect("adducts", "Adducts", ["Auto"] + ADDUCTS_POS[1:] + ADDUCTS_NEG[1:], ["Auto"]),
                select("instrument", "Instrument / CE", ["Q-TOF 20 eV", "Q-TOF 40 eV", "Orbitrap HCD NCE 30", "Triple quadrupole"]), select("out_format", "Export format", ["MSP", "MGF", "MassBank record", "JSON"])],
     "stages": [stage("Molecular properties", "properties"), stage("Adduct prediction", "platform_fit", ["ph", "column"]),
                stage("MS1 theory: precursor m/z and isotope pattern", "ms1", ["match", "hits", "formula"]),
                stage("MS2 prediction", "ms2", ["observed", "library", "similarity", "metfrag", "denovo"]),
                stage("Export", "cleanup", [])]},
]

GROUPS = [{"id": "pipelines", "title": "Pipelines"}, {"id": "modules", "title": "Modules"}, {"id": "dmpk", "title": "DMPK"}]


def spec():
    return {"groups": GROUPS, "pipelines": PIPELINES, "modules": MODULES}
