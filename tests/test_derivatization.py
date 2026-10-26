"""衍生化引擎、官能团反应式，以及 Platform 页（四个 tile + Why and how 反应块）和 Derivatization 页（表 + 反应块）的测试（质量均为单同位素，容差 0.002）。"""
import pytest
from fastapi.testclient import TestClient

from mzdesk.api.app import app
from mzdesk.platform_fit import derivatization as d

client = TestClient(app)
GLUCOSE = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"
GLUCOSE_OPEN = "OCC(O)C(O)C(O)C(O)C=O"
CAFFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
PC = "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCCCC"


def _row(smiles, option):
    return d.derive(smiles, [option])["results"][0]


def _platform(smiles):
    r = client.post("/api/modules/platform_fit/run", json={"structure": {"mode": "SMILES / InChI", "value": smiles}})
    assert r.status_code == 200, r.text
    return r.json()["outputs"]


def _tiles(smiles):
    """{tile label: tile} of the Platform page."""
    return {t["label"]: t for t in _platform(smiles)["platform"]["tiles"]}


def _options(smiles):
    """{reaction option label: option} of the Platform page's Why and how block ({} when there is none)."""
    return {o["label"]: o for o in _platform(smiles).get("reactions", {}).get("options", [])}


def _help(tile):
    h = tile.get("help", [])
    return " ".join([h] if isinstance(h, str) else h)


def test_glycine_tms():
    r = _row("NCC(=O)O", "GC · TMS (MSTFA / BSTFA)")
    assert r["n"] == 3 and r["formula"] == "C11H29NO2Si3" and r["mass"] == pytest.approx(291.151, abs=0.002)


def test_glucose_tms_and_meox():
    r = _row(GLUCOSE, "GC · TMS (MSTFA / BSTFA)")
    assert r["n"] == 5 and r["formula"] == "C21H52O6Si5" and r["mass"] == pytest.approx(540.261, abs=0.002)
    r = _row(GLUCOSE_OPEN, "GC · MeOx + TMS")
    assert r["n_meox"] == 1 and r["n_tms"] == 5
    assert r["formula"] == "C22H55NO6Si5" and r["mass"] == pytest.approx(569.288, abs=0.002)


def test_dansyl_on_glycine_with_heavy_label():
    r = _row("NCC(=O)O", "CIL · Dansyl chloride (amine / phenol)")
    assert r["n"] == 1 and r["formula"] == "C14H16N2O4S" and r["mass"] == pytest.approx(308.083, abs=0.002)
    assert r["heavy"] - r["mass"] == pytest.approx(2.007, abs=0.001)


def test_dmpa_on_acetic_acid():
    r = _row("CC(=O)O", "CIL · DmPA bromide (carboxyl)")
    assert r["formula"] == "C12H15NO3" and r["mass"] == pytest.approx(221.105, abs=0.002)


def test_opa_delta_and_no_group_row():
    assert d.formula_string(d.OPA) == "C10H8OS"
    res = d.derive(CAFFEINE, ["GC · TMS (MSTFA / BSTFA)"])
    r = res["results"][0]
    assert res["n"] == [0] and r["formula"] is None and r["mass"] is None and r["steps"] == [] and r["text"] == ""


def test_platform_tiles_for_caffeine_glucose_pc_and_glycine():
    caffeine = _tiles(CAFFEINE)
    assert list(caffeine) == ["Lipidomics", "Metabolomics", "CIL LC-MS", "GC-MS"]
    assert [caffeine[k]["status"] for k in caffeine] == ["no", "yes", "no", "yes"] and caffeine["GC-MS"]["value"] == "Volatile as is"
    glucose = _tiles(GLUCOSE)
    assert glucose["GC-MS"]["status"] == "partial" and glucose["GC-MS"]["value"] == "After derivatization" and glucose["GC-MS"]["chips"] == ["MeOx + TMS"]
    assert glucose["CIL LC-MS"]["status"] == "yes" and glucose["CIL LC-MS"]["chips"] == ["hydroxyl"]
    pc = _tiles(PC)
    assert [pc[k]["status"] for k in pc] == ["yes", "no", "no", "no"] and pc["Lipidomics"]["value"] == "Glycerophospholipid"
    assert pc["GC-MS"]["value"] == "Not volatile" and "quaternary ammonium" in _help(pc["GC-MS"])
    gly = _tiles("NCC(=O)O")
    assert gly["CIL LC-MS"]["status"] == "yes" and gly["CIL LC-MS"]["chips"] == ["amine / phenol", "carboxyl"]
    assert "amine / phenol (dansyl)" in _help(gly["CIL LC-MS"]) and "carboxyl (DmPA)" in _help(gly["CIL LC-MS"])


def test_platform_tile_help_has_the_rule_hits_and_the_lipid_features():
    pc = _tiles(PC)
    lipid, metab = _help(pc["Lipidomics"]), _help(pc["Metabolomics"])
    assert "Glycerophospholipid backbone → Lipidome" in lipid and "≥C12 tail: C16 of C16; acyl-linked" in lipid and "phosphocholine" in lipid
    assert "logP 10.61" in metab and "MW within" in metab and "≥C12" not in metab
    both = _tiles("C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2O")          # testosterone: both platforms report it
    assert both["Lipidomics"]["status"] == both["Metabolomics"]["status"] == "yes" and both["Lipidomics"]["value"] == "Sterol lipid"
    neither = _tiles("[Na+].[Cl-]")
    assert neither["Lipidomics"]["status"] == neither["Metabolomics"]["status"] == "no" and "Neither" in _help(neither["Metabolomics"])


def test_platform_reaction_block_for_tyrosine_glucose_caffeine_and_pc():
    tyr = _platform("NC(Cc1ccc(O)cc1)C(=O)O")
    assert tyr["reactions"]["svg"].lstrip().startswith(("<?xml", "<svg"))
    options = {o["label"]: o for o in tyr["reactions"]["options"]}
    assert list(options) == ["CIL · Dansyl chloride", "CIL · DmPA bromide", "GC · TMS"]
    dns, gc = options["CIL · Dansyl chloride"], options["GC · TMS"]
    assert all(o.get("svg") for o in options.values()) and dns["svg"] != gc["svg"]
    assert [(s["from"], s["to"]) for s in dns["steps"]] == [("R–NH₂", "R–NH–Dns"), ("Ar–OH", "Ar–O–Dns")]
    assert dns["formula"] == "C₃₃H₃₃N₃O₇S₂" and dns["mass"] == pytest.approx(647.176, abs=0.002)
    assert dns["heavy"] - dns["mass"] == pytest.approx(4.0134, abs=0.001)
    assert [(s["from"], s["to"]) for s in gc["steps"]] == [("R–COOH", "R–COO–TMS"), ("Ar–OH", "Ar–O–TMS"), ("R–NH₂", "R–N(TMS)₂")]
    assert "heavy" not in gc and gc["mass"] == pytest.approx(469.232, abs=0.002)
    assert any("EI key ions: m/z 73" in h for h in gc["help"]) and any("only 1 silyl group" in h for h in gc["help"])
    glu = _options(GLUCOSE)
    assert list(glu) == ["CIL · Dansyl chloride (base)", "GC · MeOx + TMS"]
    assert glu["CIL · Dansyl chloride (base)"]["steps"] == [{"from": "R–OH", "reagent": "Dns-Cl + base", "to": "R–O–Dns", "count": 5}]
    assert glu["GC · MeOx + TMS"]["steps"][0] == {"from": "ring hemiacetal", "reagent": "MeOx", "to": "R–CH=N–OCH₃"}
    assert glu["GC · MeOx + TMS"]["steps"][1]["count"] == 5 and glu["GC · MeOx + TMS"]["mass"] == pytest.approx(569.288, abs=0.002)
    for smiles in (CAFFEINE, PC):                                       # nothing reacts: the block is returned with an empty options list, no placeholder text
        assert _platform(smiles)["reactions"] == {"options": []}


def test_dansyl_option_highlights_only_its_own_reactive_groups(monkeypatch):
    from mzdesk.api.runners import structure as runner

    res = d.derive("NC(Cc1ccc(O)cc1)C(=O)O",
                   ["CIL · Dansyl chloride (amine / phenol)", "CIL · DmPA bromide (carboxyl)"])
    seen = []

    def fake_highlight(mol, groups):
        keys = sorted(groups)
        seen.append(set(keys))
        return "<svg>" + ",".join(keys) + "</svg>"

    monkeypatch.setattr(d, "highlight_svg", fake_highlight)
    options = runner._reaction_options(res["results"], detail=True, mol=res["mol"], groups=res["groups"])
    assert seen == [{"amine1", "phenol"}, {"acid"}]
    assert options[0]["svg"] == "<svg>amine1,phenol</svg>"
    assert options[1]["svg"] == "<svg>acid</svg>"


def test_platform_lipid_with_an_amine_keeps_cil_yes_and_says_lipidomics():
    pe = "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OCCN)OC(=O)CCCCCCCCCCCCCCC"
    tiles = _tiles(pe)
    assert tiles["Lipidomics"]["status"] == "yes" and tiles["CIL LC-MS"]["status"] == "yes"
    assert "usually measured by lipidomics" in _help(tiles["CIL LC-MS"])
    assert "CIL · Dansyl chloride" in _options(pe)


def test_gc_option_is_optional_for_a_small_volatile_molecule():
    opts = _options("CCCCO")                                      # butanol: one active H, MW 74 -> GC yes, derivatization optional
    assert _tiles("CCCCO")["GC-MS"]["status"] == "yes" and list(opts)[-1] == "GC · TMS"
    assert opts["GC · TMS"]["help"][0].startswith("Optional")


def test_derivatization_page_outputs():
    r = client.post("/api/modules/derivatization/run", json={"structure": {"mode": "SMILES / InChI", "value": GLUCOSE},
                                                           "reagents": ["GC · TMS (MSTFA / BSTFA)", "CIL · DmPA bromide (carboxyl)", "CIL · Dansyl chloride, base-activated (hydroxyl)"]})
    assert r.status_code == 200, r.text
    out = r.json()["outputs"]
    assert set(out) == {"deriv_rows", "deriv_reactions"}
    rows = out["deriv_rows"]["rows"]
    assert out["deriv_rows"]["count"] == 2 and len(rows) == 2                                       # DmPA has no carboxyl to react with
    cells = rows[0]["cells"]
    assert cells[0] == {"badge": "GC · TMS", "tone": "jp"} and cells[1] == "5" and cells[2] == "C₂₁H₅₂O₆Si₅"
    assert cells[3] == "540.2610" and cells[4] == "—"
    assert any("alcohol OH ×5" in h for h in rows[0]["help"]) and any("EI key ions" in h for h in rows[0]["help"])
    base = rows[1]["cells"]
    assert base[0] == {"badge": "CIL · Dansyl chloride (base)", "tone": "report"} and base[1] == "5" and float(base[4]) - float(base[3]) == pytest.approx(10.0335, abs=0.001)
    assert "DmPA" in out["deriv_rows"]["help"] and "No reactive group" in out["deriv_rows"]["help"]
    rx = out["deriv_reactions"]
    assert rx["svg"] and [o["label"] for o in rx["options"]] == ["GC · TMS", "CIL · Dansyl chloride (base)"]
    assert all(o.get("svg") for o in rx["options"])
    assert all("formula" not in o and "mass" not in o and "heavy" not in o for o in rx["options"])                  # the table has them
    assert rx["options"][0]["steps"] == [{"from": "R–OH", "reagent": "MSTFA / BSTFA", "to": "R–O–TMS", "count": 5}]
    assert rx["help"] == "Reactive groups found: alcohol OH ×5."


def test_derivatization_page_without_a_reactive_group_or_reagent():
    body = client.post("/api/modules/derivatization/run", json={"structure": {"mode": "SMILES / InChI", "value": CAFFEINE},
                                                              "reagents": ["GC · TMS (MSTFA / BSTFA)"]}).json()
    assert body["outputs"]["deriv_rows"]["rows"] == [] and body["outputs"]["deriv_rows"]["count"] == 0
    assert body["outputs"]["deriv_reactions"] == {"options": []}
    assert any("No reactive group" in w for w in body["warnings"])
    none = client.post("/api/modules/derivatization/run", json={"structure": {"mode": "SMILES / InChI", "value": GLUCOSE}, "reagents": []}).json()
    assert none["outputs"]["deriv_rows"]["rows"] == [] and any("No reagent selected" in w for w in none["warnings"])
    assert none["outputs"]["deriv_reactions"] == {"options": []}
    default = client.post("/api/modules/derivatization/run", json={"structure": {"mode": "SMILES / InChI", "value": GLUCOSE}}).json()          # no reagents key: the defaults
    assert default["outputs"]["deriv_rows"]["count"] == 3 and default["outputs"]["deriv_reactions"]["options"]
    no_acid = client.post("/api/modules/derivatization/run", json={"structure": {"mode": "SMILES / InChI", "value": GLUCOSE},
                                                                    "reagents": ["CIL · DmPA bromide (carboxyl)"]}).json()
    assert any("None of the selected reagents" in w for w in no_acid["warnings"])


# ---------------------------------------------------------------- functional-group schemes (shared by both pages)
def test_every_reagent_has_a_scheme_and_every_step_names_a_real_group():
    assert {v[0] for v in d.REAGENTS.values()} == set(d.SCHEMES)
    for name, scheme in d.SCHEMES.items():
        assert scheme["abbr"] and scheme["steps"]
        for key, reagent, frm, to in scheme["steps"]:
            assert key in d.GROUPS or key in ("carbonyl", "hemiacetal"), (name, key)
            assert reagent and frm and to and frm != to, (name, key)
    # every group a reagent labels has a step, so the scheme explains every label of the mass
    for option, (name, rule, _delta, _cil) in d.REAGENTS.items():
        keys = {k for k, *_ in d.SCHEMES[name]["steps"]}
        assert {("carbonyl" if k in d.CARBONYL else k) for k in rule} <= keys, option


def test_scheme_counts_are_numbers_of_groups_and_agree_with_the_labels():
    res = d.derive("NC(CO)C(=O)O", list(d.REAGENTS))["results"]                  # serine: amine, alcohol, acid
    by = {r["reagent"]: r for r in res}
    tms = by["GC · TMS (MSTFA / BSTFA)"]
    assert [(f, t, n) for f, _r, t, n in tms["steps"]] == [("R–COOH", "R–COO–TMS", 1), ("R–OH", "R–O–TMS", 1), ("R–NH₂", "R–N(TMS)₂", 1)] and tms["n"] == 4
    base = by["CIL · Dansyl chloride, base-activated (hydroxyl)"]
    assert [(f, n) for f, _r, _t, n in base["steps"]] == [("R–OH", 1), ("R–NH₂", 1)] and base["n"] == 2
    assert by["CIL · Dansyl chloride (amine / phenol)"]["steps"] == [("R–NH₂", "Dns-Cl", "R–NH–Dns", 1)]
    assert by["CIL · Dansylhydrazine (carbonyl)"]["steps"] == [] and by["CIL · Dansylhydrazine (carbonyl)"]["n"] == 0
    ketone = d.derive("CC(=O)CC(=O)C", ["GC · MeOx + TMS", "CIL · Dansylhydrazine (carbonyl)", "UV/FLD · DNPH (carbonyl)"])["results"]    # two C=O: counted once per group
    assert ketone[0]["steps"][0] == ("R₂C=O", "MeOx", "R₂C=N–OCH₃", 2) and ketone[1]["steps"] == [("R₂C=O", "Dns-NHNH₂", "R₂C=N–NH–Dns", 2)]
    assert ketone[2]["steps"] == [("R₂C=O", "DNPH", "R₂C=N–NH–DNP", 2)]


def test_reacting_groups_and_group_text():
    res = d.derive(GLUCOSE, ["CIL · DmPA bromide (carboxyl)", "GC · MeOx + TMS"])
    assert d.group_text(res["groups"]) == ["alcohol OH ×5"]
    used = d.reacting_groups(res["results"], res["groups"])
    assert list(used) == ["alcohol"]                                              # DmPA has nothing to react with, MeOx + TMS labels the OH
    assert d.reacting_groups(d.derive(CAFFEINE, [])["results"], {}) == {}


# ---------------------------------------------------------------- CIL 通道规则（Zhao 2019）与中和
DANSYL_AMINE = "CIL · Dansyl chloride (amine / phenol)"


def test_histidine_dansyl_labels_amine_and_imidazole_nh():
    r = _row("NC(Cc1cnc[nH]1)C(=O)O", DANSYL_AMINE)
    assert r["n"] == 2 and r["labels"] == {"amine1": 1, "arom_nh_cil": 1}
    assert r["formula"] == "C30H31N5O6S2" and r["mass"] == pytest.approx(621.172, abs=0.002)
    # TMS 里咪唑 NH 只算一次：COOH 1 + NH2 2 + [nH] 1
    assert _row("NC(Cc1cnc[nH]1)C(=O)O", "GC · TMS (MSTFA / BSTFA)")["n"] == 4


def test_amide_like_aromatic_nh_is_not_a_cil_amine():
    assert _row("O=c1cc[nH]c(=O)[nH]1", DANSYL_AMINE)["n"] == 0


def test_cysteine_thiol_goes_to_hydroxyl_channel():
    r = _row("NC(CS)C(=O)O", "CIL · Dansyl chloride, base-activated (hydroxyl)")
    assert r["labels"] == {"thiol": 1, "amine1": 1} and r["n"] == 2
    res = d.derive("NC(CS)C(=O)O", [])
    assert "hydroxyl (base-activated dansyl)" in d.cil_assess(res["groups"])["channels"]


def test_enol_oh_counted_once():
    res = d.derive("CC(O)=CC(C)=O", ["GC · TMS (MSTFA / BSTFA)", "CIL · Dansyl chloride, base-activated (hydroxyl)"])
    assert [r["n"] for r in res["results"]] == [1, 1] and "alcohol" not in res["groups"]
    assert "hydroxyl (base-activated dansyl)" in d.cil_assess(res["groups"])["channels"]


def test_aromatic_ketone_is_a_carbonyl():
    flavone = "O=c1cc(-c2ccccc2)oc2ccccc12"
    assert _row(flavone, "CIL · Dansylhydrazine (carbonyl)")["labels"] == {"ketone_ar": 1}
    r = _row(flavone, "GC · MeOx + TMS")
    assert r["n_meox"] == 1 and r["formula"] == "C16H13NO2"        # C15H10O2 + CH3N
    assert _row(flavone, "UV/FLD · DNPH (carbonyl)")["formula"] == "C21H14N4O5"


def test_zwitterion_is_neutralized():
    res = d.derive("[NH3+]CC(=O)[O-]", ["GC · TMS (MSTFA / BSTFA)"])
    r = res["results"][0]
    assert res["neutralized"] and r["n"] == 3 and r["formula"] == "C11H29NO2Si3" and r["mass"] == pytest.approx(291.151, abs=0.002)
    r = client.post("/api/modules/platform_fit/run", json={"structure": {"mode": "SMILES / InChI", "value": "[NH3+]CC(=O)[O-]"}}).json()
    assert {t["label"]: t["status"] for t in r["outputs"]["platform"]["tiles"]}["GC-MS"] == "partial"
    assert any("neutralized" in w for w in r["warnings"])


def test_quaternary_ammonium_stays_charged_and_lipid_cil_note():
    tiles = _tiles(PC)
    assert tiles["GC-MS"]["status"] == "no" and tiles["CIL LC-MS"]["status"] == "no"
    assert "Lipid: usually measured by lipidomics, not CIL." in tiles["CIL LC-MS"]["help"]
    assert "quaternary ammonium" in _help(tiles["GC-MS"]) and "formal charge" in _help(tiles["GC-MS"])


# ---------------------------------------------------------------- MeOx 打开还原糖的环
def test_meox_opens_reducing_sugar_ring():
    r = _row(GLUCOSE, "GC · MeOx + TMS")
    assert r["n_meox"] == 1 and r["n_tms"] == 5 and r["n"] == 6
    assert r["formula"] == "C22H55NO6Si5" and r["mass"] == pytest.approx(569.288, abs=0.002)
    assert "hemiacetal ×1 (ring opened, MeOx)" in d.derive(GLUCOSE, ["GC · MeOx + TMS"])["results"][0]["text"]
    # 只加 TMS、以及 dansylhydrazine / DNPH 不受影响
    assert _row(GLUCOSE, "GC · TMS (MSTFA / BSTFA)")["mass"] == pytest.approx(540.261, abs=0.002)
    assert _row(GLUCOSE, "CIL · Dansylhydrazine (carbonyl)")["n"] == 0
    assert _row(GLUCOSE, "UV/FLD · DNPH (carbonyl)")["n"] == 0
    assert "569.29" in _help(_tiles(GLUCOSE)["GC-MS"])


def test_non_reducing_ring_gets_no_meox():
    methyl_glucoside = "CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"
    r = _row(methyl_glucoside, "GC · MeOx + TMS")
    assert r["n_meox"] == 0 and r["n"] == 4
    assert r["mass"] == pytest.approx(_row(methyl_glucoside, "GC · TMS (MSTFA / BSTFA)")["mass"], abs=1e-6)
