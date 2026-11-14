"""Unit tests for the JP18 known-methods build (mzdesk/methods/jp_build): small text snippets only, no PDFs, no network."""
from __future__ import annotations

from mzdesk.methods.jp_build import analytes, header, merge, parse, resolve, segment
from mzdesk.methods.jp_build.extract import span_text


def _span(font, text):
    return {"font": font, "text": text}


# ---------------------------------------------------------------- glyph mapping
def test_glyph_mapping_special_fonts():
    assert span_text(_span("Gj202Pif", "z")) == "%"                     # the percent / degree font
    assert span_text(_span("Gj202Pif", "9")) == "°"
    assert span_text(_span("Oh107Pi", "?")) == "′"                      # prime
    assert span_text(_span("Go068Pif", "m")) == "μ"                     # Symbol-encoded Greek: 5 μm
    assert span_text(_span("Go068Pif", "d")) == "δ"
    assert span_text(_span("Go087Pif", "Identiˆcation")) == "Identification"   # fi ligature glyph
    assert span_text(_span("IwaGoTPT-Md", "＝＋－")) == "=+-"             # full-width operators


# ---------------------------------------------------------------- monograph header
def test_header_cas_formula_chemical_name():
    rows = ["アセブトロール塩酸塩", "C18H28N2O4.HCl: 372.89",
            "N-{3-Acetyl-4-[(2RS)-2-hydroxy-", "3-(1-methylethyl)aminopropyloxy]phenyl}butanamide",
            "monohydrochloride", "[34381-68-5]",
            "Acebutolol Hydrochloride, when dried, contains", "not less than 98.0% and not more than 102.0% of",
            "acebutolol hydrochloride (C18H28N2O4.HCl)."]
    h = header.parse_header(rows, "Acebutolol Hydrochloride")
    assert h["jp_name"] == "アセブトロール塩酸塩"
    assert h["formula"] == "C18H28N2O4.HCl" and h["mw"] == 372.89
    assert h["cas"] == ["34381-68-5"]
    assert h["chemical_name"].startswith("N-{3-Acetyl-4-[(2RS)-2-hydroxy-3-(1-methylethyl)aminopropyloxy]phenyl}butanamide")
    assert h["chemical_name"].endswith("monohydrochloride")
    assert h["substance"] == "acebutolol hydrochloride"


def test_header_cas_note_and_preparation_statement():
    h = header.parse_header(["ピペミド酸水和物", "C14H17N5O3.3H2O: 357.36", "8-Ethyl-5-oxo-2-(piperazin-1-yl)-",
                             "5,8-dihydropyrido[2,3-d]pyrimidine-", "6-carboxylic acid trihydrate", "[51940-44-4, anhydride]"], "Pipemidic Acid Hydrate")
    assert h["cas"] == ["51940-44-4"] and h["cas_note"] == "anhydride"
    assert "5,8-dihydropyrido" in h["chemical_name"] and "pyrimidine-6-carboxylic" in h["chemical_name"]
    t = header.parse_header(["ペミロラストカリウム錠", "Pemirolast Potassium Tablets contain not less than 95.0% and not more than 105.0% of the labeled",
                             "amount of pemirolast potassium (C10H7KN6O:", "266.30)."], "Pemirolast Potassium Tablets")
    assert t["substance"] == "pemirolast potassium" and t["formula"] is None
    assert t["formula_stmt"] == "C10H7KN6O" and t["mw_stmt"] == 266.30


def test_header_combination_preparation():
    h = header.parse_header(["x", "Foo Tablets contain not less than 95.0% of the labeled amount of candesartan cilexetil (C33H34N6O6: 610.66) and "
                             "amlodipine besylate (C20H25ClN2O5.C6H6O3S: 567.05)."], "Candesartan Cilexetil and Amlodipine Besylate Tablets")
    assert [s["name"] for s in h["substances"]] == ["candesartan cilexetil", "amlodipine besylate"]


# ---------------------------------------------------------------- gradient table
def test_gradient_table_parse():
    text = ("Operating conditions—\nFlowing of mobile phase: Control the gradient as directed in the following table.\n"
            "Time after injection  Mobile phase A  Mobile phase B\nof sample (min)  (vol%)  (vol%)\n"
            "0 – 40  60 →0  40 →100\n40 – 50  0  100\nFlow rate: 1.0 mL per minute.")
    assert parse.gradient_table(text) == [[0, 40, 60, 0, 40, 100], [40, 50, 0, 0, 100, 100]]
    assert parse.gradient_table("Detector: UV\nColumn: x") is None


def test_parse_block_full_lc():
    text = ("Operating conditions—\nDetector: An ultraviolet absorption photometer (wavelength: 254 nm).\n"
            "Column: A stainless steel column 4.6 mm in inside diameter and 15 cm in length, packed with octadecylsilanized\n"
            "silica gel for liquid chromatography (5 μm in particle diameter).\nColumn temperature: A constant temperature of about\n40 °C.\n"
            "Mobile phase: Adjust the pH of 0.05 mol/L potassium dihydrogen phosphate TS to 7.0 with triethylamine. To 350 mL of this solution add\n"
            "350 mL of methanol and 300 mL of acetonitrile.\nFlow rate: Adjust so that the retention time of nicergoline\nis about 25 minutes.\n"
            "Time span of measurement: About 2 times as long as the retention time of nicergoline.")
    p = parse.parse_block({"text": text, "technique_hint": "LC"})
    assert (p["technique"], p["block_type"], p["detector_type"], p["wavelength_nm"]) == ("LC", "full", "UV", 254)
    assert (p["column_id_mm"], p["column_length_mm"], p["packing_class"], p["particle_or_film_um"], p["column_temp_c"]) == (4.6, 150, "C18", 5, 40)
    assert p["mobile_phase_ph"] == 7.0                                   # not 0.05 (a concentration)
    assert p["organic_solvents"] == ["acetonitrile", "methanol"] and p["elution_mode"] == "isocratic"
    assert (p["flow_target_compound"], p["flow_target_rt_min"], p["run_time_rt_multiple"]) == ("nicergoline", 25, 2)
    assert p["flow_rate_value"] is None


def test_parse_block_partial_xref_and_tolerance():
    p = parse.parse_block({"text": "Operating conditions—\nDetector, column, column temperature, mobile phase, and\nflow rate: Proceed as directed in the operating "
                                   "conditions in the Assay (3) under Goshajinkigan Extract.\nTime span of measurement: About 3 times as long as the retention time of X.",
                           "technique_hint": "LC"})
    assert p["block_type"] == "partial_xref" and p["xref_target"] == "Assay (3) under Goshajinkigan Extract"
    assert p["xref_fields"] == ["column", "column_temperature", "detector", "flow_rate", "mobile_phase"] and p["run_time_rt_multiple"] == 3
    q = parse.parse_block({"text": "Operating conditions—\nDetector: A differential refractometer.\nColumn: A stainless steel column 7.8 mm in inside diameter and "
                                   "30 cm in length, packed with strongly acidic ion-exchange resin (9 μm in particle diameter).\n"
                                   "Column temperature: 85 ± 2 °C.\nMobile phase: Water.\nFlow rate: 0.5 mL per minute.", "technique_hint": "LC"})
    assert q["column_temp_c"] == 85 and q["packing_class"] == "cation-exchange" and q["flow_rate_value"] == 0.5


# ---------------------------------------------------------------- instructions / supplement override
def test_parse_instruction_scopes():
    a = segment.parse_instruction("Change the Purity (3) as follows:")
    assert (a["action"], a["kind"], a["scope"]) == ("change", "partial", [("Purity", 3)])
    b = segment.parse_instruction("Change the Purity (2) and Assay (1) as follows:")
    assert b["scope"] == [("Purity", 2), ("Assay", 1)]
    c = segment.parse_instruction("Change the Assay (1), (3) as follows:")
    assert c["scope"] == [("Assay", 1), ("Assay", 3)]
    assert segment.parse_instruction("Delete the following Monographs:")["action"] == "delete"
    assert segment.parse_instruction("Add the following:")["kind"] == "whole"
    assert segment.parse_instruction("Delete the arsenic item in the Purity.")["scope"] == []


def _rec(bid, mono, section, item=None, source="JP18", event=None, entry=None):
    return {"block_id": bid, "monograph": mono, "entry": entry, "section": section, "item_no": item, "source": source,
            "event_id": event, "block_type": "full", "technique": "LC"}


def _ev(eid, action, kind, mono, scope=(), source="JP18 Supplement I", entry=None):
    return {"event_id": eid, "source": source, "action": action, "kind": kind, "monograph": mono, "entry": entry, "scope": list(scope),
            "delete_items": [], "page": 1, "raw": "x", "ambiguous": False}


def test_supplement_override_rule():
    main = [_rec("M1", "Drug A", "Purity", 3), _rec("M2", "Drug A", "Purity", 5), _rec("M3", "Drug A", "Assay"),
            _rec("M4", "Drug B", "Assay"), _rec("M5", "Drug C", "Assay")]
    s1 = [_rec("S1", "Drug A", "Purity", 3, "JP18 Supplement I", "e1"), _rec("S2", "Drug B", "Assay", None, "JP18 Supplement I", "e2")]
    ev1 = [_ev("e1", "change", "partial", "Drug A", [("Purity", 3)]), _ev("e2", "change", "whole", "Drug B"),
           _ev("e3", "delete", "whole", "Drug C")]
    s2 = [_rec("T1", "Drug A", "Purity", 3, "JP18 Supplement II", "f1"), _rec("T2", "Drug D", "Assay", None, "JP18 Supplement II", "f2")]
    ev2 = [_ev("f1", "change", "partial", "Drug A", [("Purity", 3)], "JP18 Supplement II"),
           _ev("f2", "add", "whole", "Drug D", (), "JP18 Supplement II")]
    log: list = []
    cur, sup = merge.apply_supplements(main, [(s1, ev1), (s2, ev2)], log)
    ids = {r["block_id"] for r in cur}
    assert ids == {"M2", "M3", "S2", "T1", "T2"}                     # Purity (3) replaced twice; Purity (5), Assay kept; B replaced; C deleted; D added
    assert {r["block_id"] for r in sup} == {"M1", "S1", "M4", "M5"}
    t1 = next(r for r in cur if r["block_id"] == "T1")
    assert t1["revised_in"] == "JP18 Supplement II" and t1["replaces"] == ["S1"]
    assert next(r for r in sup if r["block_id"] == "M1")["superseded_by"] == "S1"
    assert {m["decision"] for m in log} == {"replaced", "deleted", "added"}


def test_ambiguous_instruction_keeps_both():
    main = [_rec("M1", "Drug A", "Purity", 3)]
    s1 = [_rec("S1", "Drug A", "Purity", 3, "JP18 Supplement I", "e1")]
    ev = _ev("e1", "other", "other", "Drug A")
    ev["ambiguous"] = True
    cur, sup = merge.apply_supplements(main, [(s1, [ev])], [])
    assert {r["block_id"] for r in cur} == {"M1", "S1"} and not sup


# ---------------------------------------------------------------- cross-reference resolution
def _parsed(bid, mono, section, item, text, entry=None):
    p = parse.parse_block({"text": text, "technique_hint": "LC"})
    return {"block_id": bid, "monograph": mono, "entry": entry, "section": section, "item_no": item, "source": "JP18",
            "cond_label": None, "sub_label": None, **p}


ASSAY = ("Operating conditions—\nDetector: An ultraviolet absorption photometer (wavelength: 230 nm).\nColumn: A stainless steel column 4.6 mm in inside "
         "diameter and 15 cm in length, packed with octadecylsilanized silica gel for liquid chromatography (5 μm in particle diameter).\n"
         "Column temperature: A constant temperature of about 40 °C.\nMobile phase: A mixture of water and methanol (1:1).\nFlow rate: 1.0 mL per minute.\n"
         "Time span of measurement: About 3 times as long as the retention time of X.")


def test_xref_resolution_whole_partial_and_named_monograph():
    recs = [_parsed("A1", "Drug X", "Assay", None, ASSAY),
            _parsed("P1", "Drug X", "Purity", 2, "Operating conditions—\nProceed as directed in the operating conditions in the Assay."),
            _parsed("P2", "Drug X Tablets", "Purity", None,
                    "Operating conditions—\nDetector, column, column temperature, mobile phase, and flow rate: Proceed as directed in the "
                    "operating conditions in the Assay under Drug X.\nTime span of measurement: About 5 times as long as the retention time of X."),
            _parsed("P3", "Drug Y", "Assay", None, "Operating conditions—\nProceed as directed in the operating conditions in the Assay (2).")]
    merge.resolve_xrefs(recs, [])
    p1, p2, p3 = recs[1], recs[2], recs[3]
    assert p1["resolved"] and p1["resolved_from"] == "A1" and p1["column_length_mm"] == 150 and p1["flow_rate_value"] == 1.0
    assert p1["run_time_rt_multiple"] == 3                           # whole-block reference also covers the time span
    assert p2["resolved"] and p2["resolved_from"] == "A1" and p2["detector_type"] == "UV"
    assert p2["run_time_rt_multiple"] == 5                           # stated in the block itself, not overwritten
    assert p3["resolved"] is False and p3["resolve_note"]            # no such monograph / test: unresolved, not guessed


def test_xref_ambiguous_target_is_not_guessed():
    other = ASSAY.replace("230 nm", "280 nm")
    recs = [_parsed("A1", "Drug X", "Assay", 1, ASSAY), _parsed("A2", "Drug X", "Assay", 2, other),
            _parsed("P1", "Drug X", "Purity", None, "Operating conditions—\nProceed as directed in the operating conditions in the Assay."),
            _parsed("P2", "Drug X", "Purity", 2, "Operating conditions—\nProceed as directed in the operating conditions in the Assay (2).")]
    merge.resolve_xrefs(recs, [])
    assert recs[2]["resolved"] is False and "ambiguous" in recs[2]["resolve_note"]
    assert recs[3]["resolved"] and recs[3]["resolved_from"] == "A2" and recs[3]["wavelength_nm"] == 280


def test_parse_target_forms():
    assert merge.parse_target("Assay (3) under Goshajinkigan Extract") == {"sec": "Assay", "item": 3, "sub": None, "paren": None,
                                                                          "mono": "Goshajinkigan Extract"}
    assert merge.parse_target("operating conditions 1") == {"cond": "1"}
    assert merge.parse_target("4.3. under Crude Drugs Test")["chapter_item"] == "4.3"
    assert merge.parse_target("(1)") == {"sub_only": "(1)"}


# ---------------------------------------------------------------- analytes / mapping helpers
def test_analytes_crude_and_reagent_and_drug():
    crude = {"part": "crude", "monograph": "Peony Root", "entry": None, "section": "Assay", "item_title": "Paeoniflorin", "ctx_text": ""}
    assert analytes.analytes_for(crude, None)[0] == ["paeoniflorin"]
    two = dict(crude, item_title="Ephedrine and pseudoephedrine")
    assert analytes.analytes_for(two, None)[0] == ["ephedrine", "pseudoephedrine"]
    reagent = {"part": "gen", "monograph": "9.41 Reagents, Test Solutions", "entry": "Benzoic acid for assay", "section": "Purity", "item_title": None}
    assert analytes.analytes_for(reagent, None)[0] == ["benzoic acid"]
    drug = {"part": "mono", "monograph": "Acemetacin Tablets", "entry": None, "section": "Purity", "item_title": "Related substances"}
    assert analytes.analytes_for(drug, {"substance": "acemetacin"})[0] == ["acemetacin", "related substances"]


def test_resolve_helpers_offline():
    assert resolve.cas_ok("34381-68-5") and not resolve.cas_ok("34381-68-4")
    v = resolve.name_variants("Acebutolol Hydrochloride")
    assert v[0] == ("Acebutolol Hydrochloride", "name") and ("Acebutolol", "name-stripped") in v
    cache = resolve.Cache(None, offline=True)
    m = resolve.map_one(cache, "insulin human (genetical recombination)", [])
    assert m["map_method"] == "none" and "biologic" in m["reason"]
    m2 = resolve.map_one(cache, "related substances", [])
    assert m2["cid"] is None and m2["map_method"] == "none"
    s = resolve.standardize("CC(=O)Oc1ccccc1C(=O)[O-].[Na+]")                # sodium salt -> neutral parent
    assert s["parent_smiles"] == "CC(=O)Oc1ccccc1C(=O)O"


# ---------------------------------------------------------------- GC / multiple columns / segmentation
def test_parse_block_gc_and_guard_column():
    gc = parse.parse_block({"text": "Operating conditions—\nDetector: A hydrogen flame-ionization detector.\nColumn: A glass column 3 mm in inside diameter and 2 m in "
                                    "length, packed with siliceous earth for gas chromatography (177–250 μm in particle diameter) coated in 20% with 50% "
                                    "phenyl-methyl silicone polymer for gas chromatography.\nColumn temperature: A constant temperature of about 120 °C.\n"
                                    "Carrier gas: Helium.\nFlow rate: about 50 mL per minute (the retention time of X is 3 to 4 minutes.)\n"
                                    "Injection port temperature: 140 °C.\nSplit ratio: 1:5.", "technique_hint": "GC"})
    assert gc["technique"] == "GC" and gc["detector_type"] == "FID" and gc["packing_class"] == "GC-other"
    assert (gc["column_id_mm"], gc["column_length_mm"], gc["particle_or_film_um"], gc["column_temp_c"]) == (3, 2000, 177, 120)
    assert (gc["carrier_gas"], gc["flow_rate_value"], gc["flow_target_rt_min"], gc["injection_temp_c"], gc["split_ratio"]) == ("helium", 50, 3.5, 140, "1:5")
    pov = parse.parse_block({"text": "Operating conditions—\nDetector: An ultraviolet absorption photometer (wavelength: 205 nm).\nColumn: Stainless steel columns 4.0 mm in "
                                     "inside diameter and 10 mm in length, and 4.6 mm in inside diameter and 150 mm in length, packed with octadecylsilanized "
                                     "silica gel for liquid chromatography (5 μm in particle diameter), and use them as a guard column and a separation column, "
                                     "respectively.\nMobile phase: A mixture of water and methanol (19:1).\nFlow rate: 0.8 mL per min.", "technique_hint": "LC"})
    assert (pov["column_id_mm"], pov["column_length_mm"], pov["column_count"], pov["flow_rate_value"]) == (4.6, 150, 2, 0.8)
    assert "guard" in pov["notes"] or "separation" in pov["notes"]


SEG_TEXT = """@@PAGE 5 403
### Foo Tablets
フォー錠
Foo Tablets contain not less than 95.0% and not more than 105.0% of the labeled amount of foo (C1H2: 3.00).
**Assay** Weigh accurately foo, and determine as directed under Liquid Chromatography <2.01>.
Operating conditions—
Detector: An ultraviolet absorption photometer (wavelength: 254 nm).
Column: A stainless steel column 4.6 mm in inside diameter and 15 cm in length, packed with octadecylsilanized silica gel (5 μm in particle diameter).
Column temperature: A constant temperature of about 40 °C.
Mobile phase: A mixture of water and methanol (1:1).
Flow rate: 1.0 mL per minute.
System suitability—
System performance: resolution not less than 1.5.
@@PAGE 6 404
**Purity** (1) Heavy metals <1.07>—Proceed with 1.0 g of foo.
**(2)** Related substances—Dissolve 10 mg of foo in 10 mL of methanol. Perform the test as directed under Liquid Chromatography <2.01>.
Operating conditions—
Detector, column, column temperature, mobile phase, and flow rate: Proceed as directed in the operating conditions in the Assay.
Time span of measurement: About 3 times as long as the retention time of foo.
System suitability—
@@I Change the Purity (2) as follows:
**(2)** Related substances—Dissolve 5 mg of foo.
Operating conditions—
Proceed as directed in the operating conditions in the Assay.
"""


def test_segment_context_and_blocks():
    rows = segment.load_rows(SEG_TEXT)
    blocks, events, monos = segment.segment(rows, source="JP18 Supplement I", file="x", default_part="mono")
    assert len(blocks) == 3 and monos[0]["monograph"] == "Foo Tablets"
    a, p, q = blocks
    assert (a["section"], a["item_no"], a["page"], a["page_printed"]) == ("Assay", None, 5, "403")
    assert a["text"].count("\n") == 5 and "System" not in a["text"]                 # block stops before "System suitability"
    assert (p["section"], p["item_no"], p["item_title"], p["page"]) == ("Purity", 2, "Related substances", 6)
    assert q["instr"]["action"] == "change" and q["instr"]["kind"] == "partial"        # governed by "Change the Purity (2) as follows:"
    assert [e.scope for e in events if e.kind == "partial"] == [[("Purity", 2)]]
    h = header.parse_header(monos[0]["header_rows"], "Foo Tablets")
    assert h["jp_name"] == "フォー錠" and h["substance"] == "foo" and h["formula_stmt"] == "C1H2"


def test_time_span_and_film_forms():
    base = "Operating conditions—\nDetector: A hydrogen flame-ionization detector.\nColumn: A fused silica column 0.25 mm in inside diameter and 30 m in length, " \
           "coated the inside surface with a layer about 0.25 μm thick of polymethylsiloxane.\nCarrier gas: Helium.\nFlow rate: About 35 cm/second.\n"
    for tail, mult, minutes in (("About twice as long as the retention time of X, beginning after the solvent peak.", 2, None),
                                ("As long as about 4 times of the retention time of X.", 4, None),
                                ("40 minutes, beginning after the solvent peak.", None, 40),
                                ("For about 30 minutes after the sample is injected.", None, 30),
                                ("20 minutes after elution of the peak of Y.", None, None)):
        p = parse.parse_block({"text": base + "Time span of measurement: " + tail, "technique_hint": "GC"})
        assert (p["run_time_rt_multiple"], p["run_time_min"]) == (mult, minutes), tail
        assert p["particle_or_film_um"] == 0.25 and p["flow_rate_value"] == 35 and p["flow_rate_unit"] == "cm/s"


# ---------------------------------------------------------------- audit fixes: packing, analytes, sub-items, program / stepwise, multi-wavelength
def _col(text, gc=False):
    return parse.column_info(text, text.lower(), gc)["packing_class"]


def test_packing_class_never_falls_back_to_silica_for_modified_gels():
    base = "A stainless steel column 4.6 mm in inside diameter and 15 cm in length, packed with %s for liquid chromatography (5 μm in particle diameter)."
    assert _col(base % "silica gel") == "silica" and _col(base % "porous silica gel") == "silica"
    for gel, want in (("sulfonamide group bound to hexadecylsilanized silica gel", "other"),      # "sulfo" is not a cation exchanger
                      ("human albumin chemically bonded silica gel", "chiral"),
                      ("α1-acid glycoprotein binding silica gel", "chiral"),
                      ("carbamoyl group bound silica gel", "other"), ("trimethylsilanized silica gel", "other"),
                      ("octadesilcylanized silica gel", "C18"),                                    # the source's misspelling
                      ("stylene-divinylbenzene copolymer", "polymer-RP"),
                      ("slightly acidic ion-exchange silica gel", "cation-exchange"),
                      ("strongly acidic ion-exchange resin composed with a sulfonated styrene-divinylbenzene copolymer", "cation-exchange")):
        assert _col(base % gel) == want, gel
    gc = "A fused silica column 0.32 mm in inside diameter and 30 m in length, coated with %s for gas chromatography in 1.8 μm thickness."
    assert _col(gc % "6% cyanopropyl phenyl- 94% dimethyl silicone polymer", True) == "GC-6%cyanopropylphenyl"
    assert _col(gc % "14% cyanopropylphenyl-86% dimethyl silicone polymer", True) == "GC-other"
    assert _col(gc % "5% phenyl-95% dimethylpolysiloxane", True) == "GC-5%phenyl"


def test_program_stepwise_and_mobile_phase_table():
    p = parse.parse_block({"text": "Operating conditions—\nDetector: A visible spectrophotometer (wavelength: 440 nm and 570 nm).\nColumn: A stainless steel column 4 mm in "
                                   "inside diameter and 8 cm in length, packed with strongly acidic ion-exchange resin (3 μm in particle diameter).\n"
                                   "Column temperature: Varied between 50 °C and 65 °C.\nMobile phase: Buffer solutions A, B, C and D. The ion concentration of the "
                                   "mobile phase is changed stepwise from 0.10 mol/L to 1.26 mol/L.\nComponents of buffer solutions\nBuffer solution:  A  B  C  D\n"
                                   "Citric acid  8.85 g  7.72 g  6.10 g  —\nEthanol (95)  60.0 mL  —  —  60.0 mL\nReaction reagent: Mix 407 g of lithium acetate.\n"
                                   "Flow rate of mobile phase: 0.25 mL per minute.", "technique_hint": "LC"})
    assert p["column_temp_c"] == "program" and p["elution_mode"] == "gradient" and p["n_mobile_phases"] == 4
    assert p["organic_solvents"] == ["ethanol"]


def test_wavelength_options_follow_the_analyte():
    opts = parse.wavelength_options("An ultraviolet absorption photometer (wavelength: 231 nm for aconitine, hypaconitine and mesaconitine; 254 nm for jesaconitine)")
    assert opts == [(231, "aconitine, hypaconitine and mesaconitine"), (254, "jesaconitine")]
    assert parse.wavelength_options("A visible spectrophotometer (wavelength: 440 nm and 570 nm)") == []


def _drug(mono, section, title, ctx="", item_no=None, **kw):
    return {"part": "mono", "monograph": mono, "entry": None, "section": section, "item_title": title, "item_no": item_no, "ctx_text": ctx, **kw}


def test_analytes_are_what_the_test_measures():
    bulk = {"monograph_is_bulk": True, "bulk_name": "Povidone", "formula_source": None}
    # named impurity: the title compound wins over the monograph substance
    assert analytes.analytes_for(_drug("Povidone", "Purity", "1-Vinyl-2-pyrrolidone"), bulk)[0] == ["1-vinyl-2-pyrrolidone"]
    # impurity measured against the substance's own peak: both
    paro = {"monograph_is_bulk": True, "bulk_name": "Paroxetine Hydrochloride Hydrate"}
    ctx = "the area of the peak, having the relative retention time of about 0.8 to paroxetine, obtained from the sample solution is not larger"
    assert analytes.analytes_for(_drug("Paroxetine Hydrochloride Hydrate", "Purity", "4-(4-Fluorophenyl)-1-methyl-1,2,3,6-tetrahydropyridine", ctx), paro)[0] \
        == ["4-(4-fluorophenyl)-1-methyl-1,2,3,6-tetrahydropyridine", "paroxetine hydrochloride hydrate"]
    # title lists compounds and "related substances"
    pg = {"monograph_is_bulk": True, "bulk_name": "Propylene Glycol"}
    assert analytes.analytes_for(_drug("Propylene Glycol", "Purity", "Ethylene glycol, diethylene glycol and related substances"), pg)[0] \
        == ["ethylene glycol", "diethylene glycol", "related substances"]
    assert analytes.analytes_for(_drug("Sodium L-Lactate Solution", "Purity", "Methanol"), {"monograph_is_bulk": True, "bulk_name": "Sodium L-Lactate Solution"})[0] == ["methanol"]
    # "Residual solvents": the compounds named in the test's own quantity formulae
    ctx = "Amount (%) of acetone (CH3COCH3)\n=1/MT × QTa/QSa × 400 × 0.79\nAmount (%) of methanol (CH3OH)\n=1/MT × QTb/QSb × 100 × 0.79\nMT: Amount (mg) taken"
    assert analytes.analytes_for(_drug("Cilastatin Sodium", "Purity", "Residual solvents", ctx),
                                 {"monograph_is_bulk": True, "bulk_name": "Cilastatin Sodium"})[0] == ["acetone", "methanol"]
    # an unnamed impurity class keeps the substance
    assert analytes.analytes_for(_drug("Calcium Levofolinate Hydrate", "Purity", "Diastereomer"),
                                 {"monograph_is_bulk": True, "bulk_name": "Calcium Levofolinate Hydrate"})[0] == ["calcium levofolinate hydrate"]
    # crude drug whose title names a material: use the compound the test measures
    crude = {"part": "crude", "monograph": "Olive Oil", "entry": None, "section": "Purity", "item_title": "Peanut oil",
             "ctx_text": "Measure the peak heights, HT and HS, of methyl behenate of respective solutions: HT is not higher than HS."}
    assert analytes.analytes_for(crude, None)[0] == ["methyl behenate"]
    # reagent entries: related substances too
    reag = {"part": "gen", "monograph": "9.41 Reagents, Test Solutions", "entry": "Saikosaponin d for assay", "section": "Purity", "item_title": "Related substances"}
    assert analytes.analytes_for(reag, None)[0] == ["saikosaponin d", "related substances"]


def test_analytes_combination_preparation_and_title_forms():
    combo = {"prep_substances": [("cefoperazone sodium", []), ("sulbactam sodium", [])], "formula_source": "combination preparation (see analytes)"}
    ctx = "Identification (1) The retention times of cefoperazone in the chromatogram obtained from the sample solution and the standard solution are the same"
    assert analytes.analytes_for(_drug("Cefoperazone Sodium and Sulbactam Sodium for Injection", "Identification", None, ctx), combo)[0] == ["cefoperazone sodium"]
    both = "determine the peak areas of imipenem, and those of cilastatin in each solution"
    c2 = {"prep_substances": [("imipenem hydrate", []), ("cilastatin sodium", [])], "formula_source": None}
    assert analytes.analytes_for(_drug("Imipenem and Cilastatin Sodium for Injection", "Assay", None, both), c2)[0] == ["imipenem hydrate", "cilastatin sodium"]
    assert header.bulk_name("Isosorbide Mononitrate 70% /Lactose 30%") == "Isosorbide Mononitrate"
    assert header.bulk_name("Biphasic Isophane Insulin Human (Genetical Recombination)") == "Insulin Human (Genetical Recombination)"
    assert analytes.split_markers("Aconitum diester alkaloids (aconitine, jesaconitine and mesaconitine)") == ["aconitine", "jesaconitine", "mesaconitine"]
    assert analytes.title_compounds("Relative content of stearic acid and palmitic acid")[0] == ["stearic acid", "palmitic acid"]


SUBITEM_TEXT = """@@PAGE 3 401
### Foo Sodium
フー
**Molecular mass** Calculate the molecular mass of Foo Sodium.
(i) Creation of calibration curve—Weigh 20 mg of the standard.
Operating conditions—
Detector: A differential refractometer.
Column: A stainless steel column 7.5 mm in inside diameter and 30 cm in length, packed with silica gel for liquid chromatography (5 μm in particle diameter).
Flow rate: 0.5 mL per minute.
(ii) Determination of molecular mass—Dissolve 20 mg of Foo Sodium.
Operating conditions—
Detector, column, column temperature, mobile phase, and flow rate: Proceed as directed in the operating conditions in (i)
Creation of calibration curve.
**Purity (1)** Heavy metals—Proceed with 1.0 g.
**(3)** Related substances (ii)—Dissolve 10 mg.
Operating conditions—
Detector: A hydrogen flame-ionization detector.
**(4)** Total alkaloids (Alpha and beta, or gamma and
delta)—Weigh accurately about 1 g.
Operating conditions—
Detector: An ultraviolet absorption photometer (wavelength: 254 nm).
"""


def test_sub_item_titles_are_kept():
    rows = segment.load_rows(SUBITEM_TEXT)
    blocks, _, _ = segment.segment(rows, source="JP18", file="x", default_part="mono")
    assert [b["sub_disp"] for b in blocks][:2] == ["(i) Creation of calibration curve", "(ii) Determination of molecular mass"]
    from mzdesk.methods.jp_build import build
    recs = []
    for b in blocks:
        r = {"part": b["part"], "monograph": b["monograph"], "entry": b["entry"], "section": b["section"], "item_no": b["item_no"],
             "item_title": b["item_title"], "sub_disp": b["sub_disp"]}
        recs.append(build.item_label(r))
    assert recs[1] == ("(ii) Determination of molecular mass", "Molecular mass (ii) Determination of molecular mass")
    assert recs[2][0] == "(3) Related substances (ii)"
    assert blocks[3]["item_title"] == "Total alkaloids (Alpha and beta, or gamma and delta)"       # title wrapped over two lines
