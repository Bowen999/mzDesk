"""衍生化引擎、官能团反应式，以及 Platform 页（四个 tile + Why and how 反应块）和 Derivatization 页（表 + 反应块）的测试（质量均为单同位素，容差 0.002）。"""
import pytest
from fastapi.testclient import TestClient
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

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
    assert list(caffeine) == ["Lipidomics", "Metabolomics", "CIL Metabolomics", "GC-MS"]
    assert [caffeine[k]["status"] for k in caffeine] == ["no", "yes", "no", "yes"] and caffeine["GC-MS"]["value"] == "Volatile as is"
    glucose = _tiles(GLUCOSE)
    assert glucose["GC-MS"]["status"] == "partial" and glucose["GC-MS"]["value"] == "After derivatization" and glucose["GC-MS"]["chips"] == ["MeOx + TMS"]
    assert glucose["CIL Metabolomics"]["status"] == "yes" and glucose["CIL Metabolomics"]["chips"] == ["hydroxyl"]
    pc = _tiles(PC)
    assert [pc[k]["status"] for k in pc] == ["yes", "no", "no", "no"] and pc["Lipidomics"]["value"] == "Glycerophospholipid"
    assert pc["GC-MS"]["value"] == "Not volatile" and "quaternary ammonium" in _help(pc["GC-MS"])
    gly = _tiles("NCC(=O)O")
    assert gly["CIL Metabolomics"]["status"] == "yes" and gly["CIL Metabolomics"]["chips"] == ["amine / phenol", "carboxyl"]
    assert "amine / phenol (dansyl)" in _help(gly["CIL Metabolomics"]) and "carboxyl (DmPA)" in _help(gly["CIL Metabolomics"])


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


@pytest.mark.parametrize(("smiles", "reagent", "n_tags"), [
    ("NC(Cc1ccc(O)cc1)C(=O)O", "CIL · Dansyl chloride (amine / phenol)", 2),
    ("NC(Cc1cnc[nH]1)C(=O)O", "CIL · Dansyl chloride (amine / phenol)", 2),
    (GLUCOSE, "CIL · Dansyl chloride, base-activated (hydroxyl)", 5),
    ("NC(CS)C(=O)O", "CIL · Dansyl chloride, base-activated (hydroxyl)", 2),
])
def test_dansyl_full_product_has_every_tag_and_the_same_formula(smiles, reagent, n_tags):
    res = d.derive(smiles, [reagent])
    r = res["results"][0]
    product = d.reaction_product(res["mol"], r, res["groups"])
    assert r["n"] == n_tags and product is not None
    assert d.formula_string(d.parent_formula(product)) == r["formula"]
    assert sum(a.GetAtomicNum() == 16 for a in product.GetAtoms()) == sum(a.GetAtomicNum() == 16 for a in res["mol"].GetAtoms()) + n_tags
    targets = [i for key in r["labels"] for i in res["groups"][key]["atoms"]]
    assert all(any(n.GetAtomicNum() == 16 for n in product.GetAtomWithIdx(i).GetNeighbors()) for i in targets)
    svg = d.reaction_product_svg(res["mol"], r, res["groups"])
    assert svg.lstrip().startswith(("<?xml", "<svg"))


@pytest.mark.parametrize(("smiles", "n_tags"), [("CC(=O)O", 1), ("OC(=O)CCC(=O)O", 2)])
def test_dmpa_full_product_has_every_ester_and_the_same_formula(smiles, n_tags):
    res = d.derive(smiles, ["CIL · DmPA bromide (carboxyl)"])
    r = res["results"][0]
    product = d.reaction_product(res["mol"], r, res["groups"])
    assert r["n"] == n_tags and product is not None
    assert d.formula_string(d.parent_formula(product)) == r["formula"]
    original_atoms = res["mol"].GetNumAtoms()
    targets = [i for key in r["labels"] for i in res["groups"][key]["atoms"]]
    assert all(any(n.GetIdx() >= original_atoms and n.GetAtomicNum() == 6
                   for n in product.GetAtomWithIdx(i).GetNeighbors()) for i in targets)
    assert d.reaction_product_svg(res["mol"], r, res["groups"]).lstrip().startswith(("<?xml", "<svg"))


@pytest.mark.parametrize("smiles", ["CC(=O)CC(=O)C", "O=Cc1ccccc1", "O=c1cc(-c2ccccc2)oc2ccccc12"])
def test_dansylhydrazine_full_product_has_every_hydrazone_and_the_same_formula(smiles):
    res = d.derive(smiles, ["CIL · Dansylhydrazine (carbonyl)"])
    r = res["results"][0]
    product = d.reaction_product(res["mol"], r, res["groups"])
    assert product is not None and d.formula_string(d.parent_formula(product)) == r["formula"]
    assert product.GetNumAtoms() == res["mol"].GetNumAtoms() + 17 * r["n"]
    assert len(product.GetSubstructMatches(Chem.MolFromSmarts("[#6]=[N][N][S]"))) == r["n"]
    assert d.reaction_product_svg(res["mol"], r, res["groups"]).lstrip().startswith(("<?xml", "<svg"))


_SVG = ("<?xml", "<svg")
_PRODUCT_CASES = ["NCC(=O)O", "NCCc1ccc(O)c(O)c1", "NCCCCC(N)C(=O)O", GLUCOSE, "NC(CS)C(=O)O", "NC(Cc1c[nH]c2ccccc12)C(=O)O",
                  "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C", "O=Cc1ccccc1", "OC(=O)C1CCCN1", "C[N+](C)(C)CC(O)CC(=O)[O-]", "O=c1cc[nH]c(=O)[nH]1"]


@pytest.mark.parametrize("smiles", _PRODUCT_CASES)
def test_every_reagent_full_product_has_the_reported_formula(smiles):
    """All eleven reagents (GC, CIL, UV / FLD): the drawn product is the molecule whose formula derive() reports."""
    res = d.derive(smiles, list(d.REAGENTS))
    for r in res["results"]:
        if not r["n"]:
            assert d.reaction_product(res["mol"], r, res["groups"]) is None
            continue
        product = d.reaction_product(res["mol"], r, res["groups"])
        assert d.formula_string(d.parent_formula(product)) == r["formula"], (smiles, r["name"])
        assert d.reaction_product_svg(res["mol"], r, res["groups"]).lstrip().startswith(_SVG), (smiles, r["name"])


@pytest.mark.parametrize(("smiles", "reagent", "product"), [
    ("NCCc1ccc(O)c(O)c1", "UV/FLD · OPA + 2-mercaptoethanol (primary amine)", "OCCSc1c2ccccc2cn1CCc1ccc(O)c(O)c1"),
    ("NCCc1ccc(O)c(O)c1", "GC · Acylation (TFAA)", "O=C(NCCc1ccc(OC(=O)C(F)(F)F)c(OC(=O)C(F)(F)F)c1)C(F)(F)F"),
    ("NCC(=O)O", "UV/FLD · FMOC-Cl (amine)", "O=C(O)CNC(=O)OCC1c2ccccc2-c2ccccc21"),
    ("CC(=O)C", "UV/FLD · DNPH (carbonyl)", "CC(C)=NNc1ccc([N+](=O)[O-])cc1[N+](=O)[O-]"),
    ("NCC(=O)O", "GC · TMS (MSTFA / BSTFA)", "C[Si](C)(C)OC(=O)CN([Si](C)(C)C)[Si](C)(C)C"),
])
def test_gc_and_uv_full_products(smiles, reagent, product):
    res = d.derive(smiles, [reagent])
    got = d.reaction_product(res["mol"], res["results"][0], res["groups"])
    assert Chem.MolToSmiles(got) == Chem.MolToSmiles(Chem.MolFromSmiles(product))


def test_full_product_svgs_on_the_platform_and_derivatization_pages():
    platform = _options("NC(Cc1ccc(O)cc1)C(=O)O")
    for label in ("CIL · Dansyl chloride", "CIL · DmPA bromide", "GC · TMS"):
        assert platform[label]["product_svg"].lstrip().startswith(_SVG)

    body = client.post("/api/modules/derivatization/run", json={
        "structure": {"mode": "SMILES / InChI", "value": GLUCOSE},
        "reagents": ["GC · TMS (MSTFA / BSTFA)", "GC · MeOx + TMS", "CIL · Dansyl chloride, base-activated (hydroxyl)"],
    })
    assert body.status_code == 200, body.text
    options = {o["label"]: o for o in body.json()["outputs"]["deriv_reactions"]["options"]}
    for label in ("CIL · Dansyl chloride (base)", "GC · TMS", "GC · MeOx + TMS"):
        assert options[label]["product_svg"].lstrip().startswith(_SVG)

    body = client.post("/api/modules/derivatization/run", json={
        "structure": {"mode": "SMILES / InChI", "value": "CC(=O)CC(=O)C"},
        "reagents": ["CIL · Dansylhydrazine (carbonyl)", "UV/FLD · DNPH (carbonyl)"],
    })
    assert body.status_code == 200, body.text
    options = {o["label"]: o for o in body.json()["outputs"]["deriv_reactions"]["options"]}
    for label in ("CIL · Dansylhydrazine", "UV/FLD · DNPH"):
        assert options[label]["product_svg"].lstrip().startswith(_SVG)


def test_platform_lipid_with_an_amine_keeps_cil_yes_and_says_lipidomics():
    pe = "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OCCN)OC(=O)CCCCCCCCCCCCCCC"
    tiles = _tiles(pe)
    assert tiles["Lipidomics"]["status"] == "yes" and tiles["CIL Metabolomics"]["status"] == "yes"
    assert "usually measured by lipidomics" in _help(tiles["CIL Metabolomics"])
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
    assert tiles["GC-MS"]["status"] == "no" and tiles["CIL Metabolomics"]["status"] == "no"
    assert "Lipid: usually measured by lipidomics, not CIL." in tiles["CIL Metabolomics"]["help"]
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


# ---------------------------------------------------------------- GC derivative structures (derivatize / underivatize), used by the EI MS Annotation workflow
TMS_OPTION, MEOX_OPTION, TBDMS_OPTION = d.EI_MODES["TMS"], d.EI_MODES["MeOx + TMS"], d.EI_MODES["TBDMS"]
METABOLITES = {"glycine": "NCC(=O)O", "L-alanine": "C[C@H](N)C(=O)O", "lactic acid": "CC(O)C(=O)O", "L-serine": "N[C@@H](CO)C(=O)O", "succinic acid": "OC(=O)CCC(O)=O",
               "citric acid": "OC(=O)CC(O)(CC(O)=O)C(O)=O", "pyruvic acid": "CC(=O)C(=O)O", "urea": "NC(N)=O", "indole": "c1ccc2[nH]ccc2c1",
               "L-tryptophan": "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O", "benzoic acid": "OC(=O)c1ccccc1", "cholesterol": "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
               "acetoacetic acid": "CC(=O)CC(O)=O", "glucose (open chain)": GLUCOSE_OPEN, "glucose (ring)": GLUCOSE, "fructose": "OC[C@H]1OC(O)(CO)[C@@H](O)[C@@H]1O",
               "caffeine": CAFFEINE, "O-phosphoserine": "N[C@@H](COP(O)(O)=O)C(O)=O", "4-pyridone": "O=c1cc[nH]cc1"}


def _key(smiles):
    return Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))


@pytest.mark.parametrize("name", list(METABOLITES))
@pytest.mark.parametrize("option", [TMS_OPTION, MEOX_OPTION, TBDMS_OPTION])
def test_the_derivative_structure_has_the_formula_and_label_count_of_derive(name, option):
    smiles = METABOLITES[name]
    built, expected = d.derivatize(smiles, option), _row(smiles, option)
    assert built["n_labels"] == expected["n"]
    assert built["formula"] == (expected["formula"] or d.formula_string(d.parent_formula(Chem.MolFromSmiles(smiles))))
    assert built["mass"] == pytest.approx(expected["mass"] if expected["mass"] else d.formula_mass(d.parent_formula(Chem.MolFromSmiles(smiles))), abs=1e-6)
    assert d.derivative_smiles(smiles, option) == (built["smiles"], built["n_labels"])


def test_known_derivative_formulas():
    assert d.derivatize("NCC(=O)O", TMS_OPTION)["formula"] == "C11H29NO2Si3"                                   # glycine, 3 TMS (M 291)
    assert d.derivatize("CC(O)C(=O)O", TMS_OPTION)["formula"] == "C9H22O3Si2"                                  # lactic acid, 2 TMS (M 234)
    assert d.derivatize(GLUCOSE, TMS_OPTION)["formula"] == "C21H52O6Si5"                                       # glucose, ring form, 5 TMS
    glucose_meox = d.derivatize(GLUCOSE, MEOX_OPTION)
    assert glucose_meox["formula"] == "C22H55NO6Si5" and (glucose_meox["n_silyl"], glucose_meox["n_meox"], glucose_meox["n_labels"]) == (5, 1, 6)
    assert d.derivatize(GLUCOSE_OPEN, MEOX_OPTION)["formula"] == "C22H55NO6Si5"
    assert d.derivatize("C[C@H](N)C(=O)O", TBDMS_OPTION)["formula"] == "C21H49NO2Si3"                          # alanine, 3 TBDMS
    assert d.nominal_mass("C11H29NO2Si3") == 291 and d.nominal_mass("C22H55NO6Si5") == 569 and d.nominal_mass(d.parse_formula("C6H12O6")) == 180


def test_a_primary_amine_takes_two_silyl_groups_and_the_structure_shows_it():
    smiles = d.derivatize("NCC(=O)O", TMS_OPTION)["smiles"]
    assert Chem.MolFromSmiles(smiles).HasSubstructMatch(Chem.MolFromSmarts("[NX3]([Si](C)(C)C)([Si](C)(C)C)C"))
    assert Chem.MolFromSmiles(smiles).HasSubstructMatch(Chem.MolFromSmarts("C(=O)O[Si](C)(C)C"))
    tbdms = d.derivatize("NCC(=O)O", TBDMS_OPTION)["smiles"]
    assert Chem.MolFromSmiles(tbdms).HasSubstructMatch(Chem.MolFromSmarts("N([Si](C)(C)C(C)(C)C)[Si](C)(C)C(C)(C)C"))


def test_meox_turns_a_carbonyl_into_an_oxime_ether_and_opens_a_sugar_ring():
    pyruvic = d.derivatize("CC(=O)C(=O)O", MEOX_OPTION)
    assert pyruvic["smiles"] == "CON=C(C)C(=O)O[Si](C)(C)C" and (pyruvic["n_silyl"], pyruvic["n_meox"]) == (1, 1)
    assert d.derivatize("CC(=O)C(=O)O", TMS_OPTION)["smiles"] == "CC(=O)C(=O)O[Si](C)(C)C"                       # TMS alone leaves the ketone
    ring = Chem.MolFromSmiles(d.derivatize(GLUCOSE, MEOX_OPTION)["smiles"])
    assert not ring.GetRingInfo().NumRings() and ring.HasSubstructMatch(Chem.MolFromSmarts("[CH]=NO[CH3]"))
    assert d.carbonyl_count(GLUCOSE) == 1 and d.carbonyl_count(GLUCOSE_OPEN) == 1 and d.carbonyl_count("NCC(=O)O") == 0 and d.carbonyl_count("CC(=O)C(=O)O") == 1
    assert d.carbonyl_count("not a smiles") == 0


def test_a_zwitterion_is_neutralised_and_a_compound_without_a_site_comes_back_unchanged():
    assert d.derivatize("[NH3+]CC(=O)[O-]", TMS_OPTION)["formula"] == "C11H29NO2Si3" and d.derivatize("[NH3+]CC(=O)[O-]", TMS_OPTION)["neutralized"]
    caffeine = d.derivatize(CAFFEINE, TMS_OPTION)
    assert (caffeine["n_labels"], caffeine["n_silyl"], caffeine["formula"]) == (0, 0, "C8H10N4O2") and _key(caffeine["smiles"]) == _key(CAFFEINE)


def test_derivatize_rejects_a_bad_structure_or_reagent():
    with pytest.raises(ValueError, match="could not parse"):
        d.derivatize("not a smiles", TMS_OPTION)
    with pytest.raises(ValueError, match="could not parse"):
        d.derivatize("", TMS_OPTION)
    with pytest.raises(ValueError, match="not a GC silylation reagent"):
        d.derivatize("NCC(=O)O", "GC · Acylation (TFAA)")


@pytest.mark.parametrize("name", [n for n in METABOLITES if n not in ("glucose (ring)", "fructose", "4-pyridone")])     # MeOx opens a sugar ring / changes a pyridone's tautomer
def test_underivatize_gives_back_the_same_compound(name):
    smiles = METABOLITES[name]
    for option in (TMS_OPTION, MEOX_OPTION, TBDMS_OPTION):
        derivative = d.derivatize(smiles, option)
        assert _key(d.underivatize(derivative["smiles"])) == _key(smiles), (name, option)


def test_underivatize_info_counts_the_groups_and_names_the_derivative():
    info = d.underivatize_info(d.derivatize("NCC(=O)O", TMS_OPTION)["smiles"])
    assert (info["n_tms"], info["n_tbdms"], info["n_meox"], info["mode"], info["changed"]) == (3, 0, 0, "TMS", True)
    assert _key(info["parent"]) == _key("NCC(=O)O")
    assert d.underivatize_info(d.derivatize("NCC(=O)O", TBDMS_OPTION)["smiles"])["mode"] == "TBDMS"
    both = d.underivatize_info(d.derivatize(GLUCOSE_OPEN, MEOX_OPTION)["smiles"])
    assert (both["mode"], both["n_tms"], both["n_meox"]) == ("MeOx + TMS", 5, 1) and _key(both["parent"]) == _key(GLUCOSE_OPEN)
    assert d.underivatize_info("CC(=O)C(=O)O")["mode"] is None and not d.underivatize_info("CC(=O)C(=O)O")["changed"]
    ring = d.underivatize(d.derivatize(GLUCOSE, MEOX_OPTION)["smiles"])
    assert _key(ring)[:14] == _key(GLUCOSE_OPEN)[:14] and _key(ring)[:14] != _key(GLUCOSE)[:14]                   # MeOx opened the ring: the open chain comes back


def test_a_library_style_derivative_smiles_maps_back_to_its_compound():
    # MassBank stores LACTIC ACID-DI-TMS and an N-silylated indole this way (explicit [H] and aromatic nitrogen included)
    assert _key(d.underivatize("O=C(O[Si](C)(C)C)C([H])(C)O[Si](C)(C)C")) == _key("CC(O)C(=O)O")
    assert _key(d.underivatize("O=C(O[Si](C)(C)C)Cc(c1)c(c2)c(ccc2)n1[Si](C)(C)C")) == _key("OC(=O)Cc1c[nH]c2ccccc12")
    assert _key(d.underivatize("C[Si](C)(C)Nc(c1)cc(cc1)O[Si](C)(C)C")) == _key("Nc1cccc(O)c1")
    assert _key(d.underivatize("CCCCCC(=O)O[Si](C)(C)C")) == _key("CCCCCC(=O)O")


@pytest.mark.parametrize("smiles", ["CCCCCCCC[Si](OCC)(OCC)OCC", "C[Si](C)(C)C(Cl)I", "c1ccccc1[SiH3]", "C[Si](C)(C)[Si](C)(C)C"])
def test_an_organosilane_is_not_a_derivative(smiles):
    info = d.underivatize_info(smiles)
    assert not info["changed"] and info["mode"] is None and _key(info["parent"]) == _key(smiles)


@pytest.mark.parametrize("smiles", ["CC1=CC(=O)c2ccccc2O1", "O=c1ccocc1", "O=c1cccccc1O", "O=c1cc(-c2ccccc2)oc2ccccc12", "O=c1cc[nH]cc1"])
def test_meox_on_a_ring_carbonyl_that_rdkit_writes_as_aromatic_is_undone(smiles):
    """chromones, pyranones, tropolones and pyridones come out of derivatize() as CON=c1...: the C=N-OCH3 test must accept an aromatic carbon."""
    derivative = d.derivatize(smiles, MEOX_OPTION)["smiles"]
    assert "CON=c" in derivative or "CON=C" in derivative
    info = d.underivatize_info(derivative)
    assert info["n_meox"] == 1 and _key(info["parent"]) == _key(smiles)


@pytest.mark.parametrize("smiles", ["C[Si](C)(C)O[Si](C)(C)C", "C[Si](C)(C)N[Si](C)(C)C", "C[Si](C)(C)O", "CC(C)(C)[Si](C)(C)O", "CC(C)(C)[Si](C)(C)O[Si](C)(C)C(C)(C)C"])
def test_a_siloxane_silazane_or_silanol_is_not_a_derivative_of_water_or_ammonia(smiles):
    """Hexamethyldisiloxane, hexamethyldisilazane and the silanols are reagent / column-bleed compounds of their own; they used to come back as "O" / "N"."""
    info = d.underivatize_info(smiles)
    assert not info["changed"] and info["mode"] is None and info["parent"] == Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    assert d.underivatize("CCO[Si](C)(C)C") == "CCO" and d.underivatize_info("C[Si](C)(C)n1ccnc1")["n_tms"] == 1       # a silyl group on a carbon-bearing O / N still is one


def test_underivatize_rejects_a_bad_structure():
    with pytest.raises(ValueError, match="could not parse"):
        d.underivatize("not a smiles")


# ---------------------------------------------------------------- phosphorus: only a free P-OH makes a compound non-volatile
TRIESTERS = [("dichlorvos", "COP(=O)(OC)OC=C(Cl)Cl"), ("triethyl phosphate", "CCOP(=O)(OCC)OCC"), ("tributyl phosphate", "CCCCOP(=O)(OCCCC)OCCCC"),
             ("triphenyl phosphate", "O=P(Oc1ccccc1)(Oc1ccccc1)Oc1ccccc1"), ("paraoxon", "CCOP(=O)(OCC)Oc1ccc([N+](=O)[O-])cc1"),
             ("DMMP", "COP(C)(=O)OC"), ("parathion", "CCOP(=S)(OCC)Oc1ccc([N+](=O)[O-])cc1")]
P_ACIDS = [("glyphosate", "OC(=O)CNCP(=O)(O)O"), ("glucose 6-phosphate", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1OP(=O)(O)O"),
           ("diethyl phosphate", "CCOP(=O)(O)OCC"), ("methyl phosphate dianion", "COP(=O)([O-])[O-]"), ("methylphosphonic acid", "CP(=O)(O)O"),
           ("diethyl thiophosphate", "CCOP(=S)(O)OCC")]


def _gc(smiles):
    mol, _ = d.neutralize(Chem.MolFromSmiles(smiles))
    return d.gc_assess(mol, d.find_groups(mol))


@pytest.mark.parametrize("name, smiles", TRIESTERS)
def test_phosphate_triesters_and_phosphonate_esters_run_in_gc(name, smiles):
    """The labile-group rule was P(=O)(O)O, which every phosphorus ester matches: dichlorvos, TEP, TBP, TPP, paraoxon and DMMP were 'not volatile'."""
    gc = _gc(smiles)
    assert gc["gc"] == "yes" and gc["volatility"] == "volatile as is", name
    tile = _tiles(smiles)["GC-MS"]
    assert tile["status"] == "yes" and tile["value"] == "Volatile as is", name


@pytest.mark.parametrize("name, smiles", P_ACIDS)
def test_phosphorus_with_a_free_p_oh_stays_not_volatile(name, smiles):
    gc = _gc(smiles)
    assert gc["gc"] == "no" and "phosphate / phosphonate" in " ".join(gc["reasons"]), name
    assert _tiles(smiles)["GC-MS"]["status"] == "no", name


# ---------------------------------------------------------------- salts and mixtures: the largest fragment, as on the other pages
PARTS_NOTE = "2 disconnected parts (salt or mixture): only the largest is analysed."
DOPAMINE = "NCCc1ccc(O)c(O)c1"


def _page(module, smiles, **extra):
    r = client.post(f"/api/modules/{module}/run", json={"structure": {"mode": "SMILES / InChI", "value": smiles}, **extra})
    assert r.status_code == 200, r.text
    return r.json()


def _table(body):
    return [[c if isinstance(c, str) else c.get("badge", "") for c in row["cells"]] for row in body["outputs"]["deriv_rows"]["rows"]]


def test_a_salt_is_derivatized_by_its_largest_fragment():
    """Dopamine hydrochloride used to give a TMS product with the chloride in the formula (C20H44ClNO2Si4, 477.2138)."""
    res = d.derive(DOPAMINE + ".Cl", [TMS_OPTION])
    r = res["results"][0]
    assert res["n_parts"] == 2 and r["n"] == 4 and r["formula"] == "C20H43NO2Si4" and r["mass"] == pytest.approx(441.2371, abs=0.002)
    product = d.derivatize(DOPAMINE + ".Cl", TMS_OPTION)
    assert product["formula"] == "C20H43NO2Si4" and product["mass"] == pytest.approx(441.2371, abs=0.002) and "Cl" not in product["smiles"]
    assert d.carbonyl_count("O=CCc1ccc(O)c(O)c1.Cl") == 1


def test_derivatization_page_gives_a_salt_the_result_and_the_warning_of_the_free_base():
    salt, free = _page("derivatization", DOPAMINE + ".Cl"), _page("derivatization", DOPAMINE)
    assert _table(salt) == _table(free) and any(row[2] == "C₂₀H₄₃NO₂Si₄" and row[3] == "441.2371" for row in _table(salt))
    assert salt["warnings"].count(PARTS_NOTE) == 1 and PARTS_NOTE not in free["warnings"]
    assert "No reactive group found in this structure." not in salt["warnings"]


def test_the_salt_warning_is_the_one_of_the_other_pages():
    from mzdesk.platform_fit import pka_engine as pe
    assert pe.prepare(DOPAMINE + ".Cl").notes[0] == PARTS_NOTE
    for module in ("platform_fit", "derivatization"):
        assert _page(module, DOPAMINE + ".[Cl-]")["warnings"].count(PARTS_NOTE) == 1, module


def test_sodium_benzoate_has_a_carboxyl_channel_and_a_gc_verdict():
    """It used to report 'No reactive group found' (Derivatization) and CIL no / GC no (Platform): the carboxylate was left as an anion."""
    salt = "[Na+].[O-]C(=O)c1ccccc1"
    body = _page("derivatization", salt)
    assert "No reactive group found in this structure." not in body["warnings"] and PARTS_NOTE in body["warnings"]
    rows = {row[0]: row for row in _table(body)}
    assert rows["CIL · DmPA bromide"][2] == "C₁₇H₁₇NO₃" and rows["GC · TMS"][2] == "C₁₀H₁₄O₂Si"        # benzoic acid + DmPA / + TMS
    tiles = _tiles(salt)
    assert tiles["CIL Metabolomics"]["status"] == "yes" and tiles["GC-MS"]["status"] != "no"
    platform = _page("platform_fit", salt)
    assert platform["warnings"].count(PARTS_NOTE) == 1
    assert {o["label"] for o in platform["outputs"]["reactions"]["options"]} == {"CIL · DmPA bromide", "GC · TMS"}


# ---------------------------------------------------------------- a quaternary N+ next to a carboxylate: protonate the carboxylate, keep the N+
QUATERNARY_ACIDS = {"betaine": "C[N+](C)(C)CC(=O)[O-]", "L-carnitine": "C[N+](C)(C)C[C@H](O)CC(=O)[O-]", "acetylcarnitine": "CC(=O)OC(CC(=O)[O-])C[N+](C)(C)C",
                    "trigonelline": "C[n+]1cccc(c1)C(=O)[O-]"}


def test_neutralize_protonates_a_carboxylate_next_to_a_quaternary_n():
    mol, changed = d.neutralize(Chem.MolFromSmiles(QUATERNARY_ACIDS["betaine"]))
    assert changed and Chem.MolToSmiles(mol) == "C[N+](C)(C)CC(=O)O"
    assert Chem.MolToSmiles(d.neutralize(Chem.MolFromSmiles("[NH3+]CC(=O)[O-]"))[0]) == "NCC(=O)O"            # an ordinary zwitterion still ends up neutral
    assert not d.neutralize(Chem.MolFromSmiles("C[N+](=O)[O-]"))[1]                                           # nitro: not a charge to remove
    pc, _ = d.neutralize(Chem.MolFromSmiles(PC))
    assert d.find_groups(pc) == {}


@pytest.mark.parametrize("name", list(QUATERNARY_ACIDS))
def test_quaternary_ammonium_acids_get_the_carboxyl_channel_and_stay_not_gc_amenable(name):
    smiles = QUATERNARY_ACIDS[name]
    tiles = _tiles(smiles)
    assert tiles["CIL Metabolomics"]["status"] == "yes" and "carboxyl" in tiles["CIL Metabolomics"]["chips"], name
    assert tiles["GC-MS"]["status"] == "no" and "Rule:" in _help(tiles["GC-MS"]), name         # [NX4+] / ionic rule untouched
    r = _row(smiles, "CIL · DmPA bromide (carboxyl)")
    assert r["n"] == 1 and r["charge"] == 1, name
    assert "No reactive group found in this structure." not in _page("derivatization", smiles)["warnings"], name


def test_a_cationic_derivative_reports_the_ions_formula_and_mz():
    """Carnitine + DmPA: the ester is a cation, C17H27N2O4+, m/z 323.1965 (not a neutral formula with an H too many or too few; RDKit as the check)."""
    carnitine = QUATERNARY_ACIDS["L-carnitine"]
    r = _row(carnitine, "CIL · DmPA bromide (carboxyl)")
    ester = Chem.MolFromSmiles("CN(C)c1ccc(C(=O)COC(=O)CC(O)C[N+](C)(C)C)cc1")
    assert r["formula"] == "C17H27N2O4" and rdMolDescriptors.CalcMolFormula(ester) == "C17H27N2O4+"
    assert r["mass"] == pytest.approx(Descriptors.ExactMolWt(ester), abs=0.0005) == pytest.approx(323.1965, abs=0.0005)
    assert r["heavy"] == pytest.approx(r["mass"] + 2 * d.C13_SHIFT, abs=1e-6)
    tms = d.derivatize(carnitine, TMS_OPTION)                       # two silyl groups: COOH and OH
    product = Chem.MolFromSmiles(tms["smiles"])
    assert tms["n_silyl"] == 2 and tms["charge"] == 1 and tms["formula"] == "C13H32NO3Si2"
    assert rdMolDescriptors.CalcMolFormula(product) == "C13H32NO3Si2+" and tms["mass"] == pytest.approx(Descriptors.ExactMolWt(product), abs=0.0005)
    body = _page("derivatization", carnitine)
    rows = {row[0]: row for row in _table(body)}
    assert rows["CIL · DmPA bromide"][2:4] == ["C₁₇H₂₇N₂O₄⁺", "323.1965"] and rows["GC · TMS"][2] == "C₁₃H₃₂NO₃Si₂⁺"
    assert any("Net charge +1" in w for w in body["warnings"])
    option = next(o for o in _page("platform_fit", carnitine)["outputs"]["reactions"]["options"] if o["label"] == "CIL · DmPA bromide")
    assert option["formula"] == "C₁₇H₂₇N₂O₄⁺" and option["mass"] == pytest.approx(323.1965, abs=0.0005)


def test_a_neutral_derivative_has_no_charge_and_an_unchanged_mass():
    r = _row("NCC(=O)O", TMS_OPTION)
    assert r["charge"] == 0 and r["formula"] == "C11H29NO2Si3" and r["mass"] == pytest.approx(291.151, abs=0.002)
    assert {row[0]: row for row in _table(_page("derivatization", "NCC(=O)O"))}["GC · TMS"][2] == "C₁₁H₂₉NO₂Si₃"
