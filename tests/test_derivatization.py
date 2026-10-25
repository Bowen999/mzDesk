"""衍生化引擎与 Platform 页 CIL / GC-MS 两节的测试（质量均为单同位素，容差 0.002）。"""
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
    assert res["n"] == [0] and res["rows"][0][6] == "no reactive group" and res["rows"][0][3] == "—"


def test_platform_gc_and_cil():
    out = _platform(CAFFEINE)
    assert out["gc"]["values"]["GC-MS"]["value"] == "yes"
    assert out["cil"]["values"]["CIL LC-MS"]["value"] == "no"
    assert _platform(GLUCOSE)["gc"]["values"]["GC-MS"]["value"] == "after derivatization"
    pc = _platform(PC)
    assert pc["gc"]["values"]["GC-MS"]["value"] == "no" and pc["cil"]["values"]["CIL LC-MS"]["value"] == "no"
    gly = _platform("NCC(=O)O")["cil"]["values"]
    assert gly["CIL LC-MS"]["value"] == "yes"
    assert "amine" in gly["Labelable groups"]["value"] and "carboxyl" in gly["Labelable groups"]["value"]


def test_derivatization_page_outputs():
    r = client.post("/api/modules/derivatization/run", json={"structure": {"mode": "SMILES / InChI", "value": GLUCOSE},
                                                           "reagents": ["GC · TMS (MSTFA / BSTFA)"]})
    out = r.json()["outputs"]
    assert out["deriv_groups"]["rows"][0][:2] == ["alcohol OH", "5"] and out["deriv_rows"]["rows"][0][2] == "5" and out["deriv_structure"]["svg"]


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
    assert r["outputs"]["gc"]["values"]["GC-MS"]["value"] == "after derivatization"
    assert any("neutralized" in w for w in r["warnings"])


def test_quaternary_ammonium_stays_charged_and_lipid_cil_note():
    out = _platform(PC)
    assert out["gc"]["values"]["GC-MS"]["value"] == "no"
    note = out["cil"]["values"]["CIL LC-MS"]["note"]
    assert note == "lipid: usually measured by lipidomics, not CIL"
    assert any("lipid" in x for x in out["gc_cil_reasons"]["items"])


# ---------------------------------------------------------------- MeOx 打开还原糖的环
def test_meox_opens_reducing_sugar_ring():
    r = _row(GLUCOSE, "GC · MeOx + TMS")
    assert r["n_meox"] == 1 and r["n_tms"] == 5 and r["n"] == 6
    assert r["formula"] == "C22H55NO6Si5" and r["mass"] == pytest.approx(569.288, abs=0.002)
    row = d.derive(GLUCOSE, ["GC · MeOx + TMS"])["rows"][0]
    assert "hemiacetal ×1 (ring opened, MeOx)" in row[1]
    # 只加 TMS、以及 dansylhydrazine / DNPH 不受影响
    assert _row(GLUCOSE, "GC · TMS (MSTFA / BSTFA)")["mass"] == pytest.approx(540.261, abs=0.002)
    assert _row(GLUCOSE, "CIL · Dansylhydrazine (carbonyl)")["n"] == 0
    assert _row(GLUCOSE, "UV/FLD · DNPH (carbonyl)")["n"] == 0
    assert "569.29" in _platform(GLUCOSE)["gc"]["values"]["Volatility"]["value"]


def test_non_reducing_ring_gets_no_meox():
    methyl_glucoside = "CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"
    r = _row(methyl_glucoside, "GC · MeOx + TMS")
    assert r["n_meox"] == 0 and r["n"] == 4
    assert r["mass"] == pytest.approx(_row(methyl_glucoside, "GC · TMS (MSTFA / BSTFA)")["mass"], abs=1e-6)
