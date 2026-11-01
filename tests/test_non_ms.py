"""Non-MS detection 页的测试：规则表、Woodward-Fieser / Fieser-Kuhn / Scott、荧光、检测器规则、API，以及 JP18 小型离线对照。

离线对照（FIXTURE）：JP18 的 27 个分析物（SMILES、JP 检测器与波长取自 mzdesk/data/known_methods/jp18_methods.jsonl.gz）加经典化合物。
JP 常用比 λmax 更低的波长换灵敏度，所以 JP 波长只作参考，不作断言。
JP18 全库对照（UV 分析物有无色基、RI / 柱后衍生化分析物是否 none / weak、|JP − λmax| 分布）在 validation/non_ms/jp18_validate.py；
这里的最后一个测试用同一批数据做粗阈值回归，数据文件不在时跳过。
"""
import gzip
import json
import statistics
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from mzdesk.api.app import app
from mzdesk.platform_fit import non_ms_rules as R
from mzdesk.spec import MODULES

client = TestClient(app)
ROOT = Path(__file__).parent.parent
JP_FILE = ROOT / "mzdesk" / "data" / "known_methods" / "jp18_methods.jsonl.gz"
ACN, MEOH = "Water/acetonitrile + phosphate", "Water/methanol + formic acid"
GLUCOSE_RING = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"
GLUCOSE_OPEN = "OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C=O"
BETA_CAROTENE = "CC1=C(C(CCC1)(C)C)C=CC(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C1=C(CCCC1(C)C)C)C)C)C)C"
CAFFEINE = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"


def _a(smiles, phase=None):
    return R.analyze(smiles, phase)


def _run(smiles, **extra):
    return client.post("/api/modules/non_ms/run", json={"structure": {"mode": "SMILES / InChI", "value": smiles}, **extra})


# (name, SMILES, expected UV strength, lambda window (nm) or None, JP detector, JP wavelength)
# first block: JP18 analytes; second block: classic compounds (windows from literature values)
FIXTURE = [
    ("nicotinamide", "C1=CC(=CN=C1)C(=O)N", "strong", (245, 275), "UV", 254),
    ("prednisolone", "C[C@]12C[C@@H]([C@H]3[C@H]([C@@H]1CC[C@@]2(C(=O)CO)O)CCC4=CC(=O)C=C[C@]34C)O", "strong", (238, 250), "UV", 247),
    ("testosterone propionate", "CCC(=O)O[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4=CC(=O)CC[C@]34C)C", "strong", (238, 250), "UV", 241),
    ("hydrocortisone", "C[C@]12CCC(=O)C=C1CC[C@@H]3[C@@H]2[C@H](C[C@]4([C@H]3CC[C@@]4(C(=O)CO)O)C)O", "strong", (238, 250), "UV", 254),
    ("cefalexin", "CC1=C(N2[C@@H]([C@@H](C2=O)NC(=O)[C@@H](C3=CC=CC=C3)N)SC1)C(=O)O", "strong", (250, 270), "UV", 254),
    ("ciprofloxacin hydrochloride hydrate", "C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O.O.Cl", "strong", (270, 300), "UV", 278),
    ("chlorpromazine hydrochloride", "CN(C)CCCN1C2=CC=CC=C2SC3=C1C=C(C=C3)Cl.Cl", "strong", (245, 270), "UV", 256),
    ("metformin hydrochloride", "CN(C)C(=N)N=C(N)N.Cl", "strong", (225, 240), "UV", 235),
    ("digoxin", "C[C@@H]1[C@H]([C@H](C[C@@H](O1)O[C@@H]2[C@H](O[C@H](C[C@@H]2O)O[C@@H]3[C@H](O[C@H](C[C@@H]3O)O[C@H]4CC[C@]5([C@@H](C4)CC[C@@H]6[C@@H]5C[C@H]([C@]7([C@@]6(CC[C@@H]7C8=CC(=O)OC8)O)C)O)C)C)C)O)O",
     "strong", (212, 225), "UV", 220),
    ("acetaminophen", "CC(=O)NC1=CC=C(C=C1)O", "strong", (235, 255), "UV", 225),
    ("amlodipine besilate", "CCOC(=O)C1=C(NC(=C(C1C2=CC=CC=C2Cl)C(=O)OC)C)COCCN.C1=CC=C(C=C1)S(=O)(=O)O", "strong", (230, 260), "UV", 237),   # DHP 238 and the chlorophenyl ring 255: the largest lambda is reported
    ("warfarin potassium", "CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)[O-].[K+]", "strong", (295, 325), "UV", 260),
    ("propranolol hydrochloride", "CC(C)NCC(COC1=CC=CC2=CC=CC=C21)O.Cl", "strong", (265, 300), "UV", 292),
    ("diazepam", "CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3", "strong", (240, 300), "UV", 254),
    ("menatetrenone", "CC1=C(C(=O)C2=CC=CC=C2C1=O)C/C=C(\\C)/CC/C=C(\\C)/CC/C=C(\\C)/CCC=C(C)C", "strong", (240, 275), "UV", 270),
    ("ergocalciferol", "C[C@H](/C=C/[C@H](C)C(C)C)[C@H]1CC[C@@H]\\2[C@@]1(CCC/C2=C\\C=C/3\\C[C@H](CCC3=C)O)C", "strong", (255, 285), "UV", 254),
    ("tocopherol", "CC1=C(C2=C(CCC(O2)(C)CCCC(C)CCCC(C)CCCC(C)C)C(=C1O)C)C", "strong", (255, 300), "UV", 292),
    ("ibuprofen piconol", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)OCC2=CC=CC=N2", "strong", (245, 275), "UV", 254),
    ("erythromycin", "CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O",
     "weak", None, "UV", 215),
    ("isosorbide mononitrate", "C1[C@@H]([C@@H]2[C@H](O1)[C@@H](CO2)O[N+](=O)[O-])O", "weak", None, "UV", 214),
    ("cyclophosphamide hydrate (JP uses 205 nm: end absorption)", "C1CNP(=O)(OC1)N(CCCl)CCCl.O", "none", None, "UV", 205),
    ("purified glucose", "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O", "none", None, "RI", None),
    ("d-mannitol", "C([C@H]([C@H]([C@@H]([C@@H](CO)O)O)O)O)O", "none", None, "RI", None),
    ("maltose hydrate", "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@@H]2[C@H](OC([C@@H]([C@H]2O)O)O)CO)O)O)O)O.O", "none", None, "RI", None),
    ("voglibose (post-column fluorescence)", "C1[C@@H]([C@@H]([C@H]([C@@H]([C@]1(CO)O)O)O)O)NC(CO)CO", "none", None, "fluorescence", 350),
    ("arbekacin sulfate (post-column fluorescence)",
     "C1C[C@H]([C@H](O[C@@H]1CN)O[C@@H]2[C@H](C[C@H]([C@@H]([C@H]2O)O[C@@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)N)O)NC(=O)[C@H](CCN)O)N)N", "weak", None, "fluorescence", 340),
    ("l-alanine (post-column ninhydrin)", "C[C@@H](C(=O)O)N", "weak", None, "visible", 570),
    # classic compounds
    ("benzene", "c1ccccc1", "strong", (245, 265), None, None),
    ("naphthalene", "c1ccc2ccccc2c1", "strong", (265, 285), None, None),
    ("anthracene", "c1ccc2cc3ccccc3cc2c1", "strong", (345, 365), None, None),
    ("caffeine", CAFFEINE, "strong", (263, 283), None, None),
    ("glucose (open-chain drawing)", GLUCOSE_OPEN, "none", None, None, None),
    ("cholesterol", "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C", "weak", None, None, None),
    ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "strong", (245, 275), None, None),
    ("quinine", "COc1ccc2nccc([C@H](O)[C@@H]3C[C@@H]4CC[N@]3C[C@@H]4C=C)c2c1", "strong", (300, 335), None, None),
    ("tryptophan", "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O", "strong", (270, 290), None, None),
    ("beta-carotene", BETA_CAROTENE, "strong", (450, 456), None, None),
    ("mesityl oxide (Woodward-Fieser 239; observed ~237)", "CC(C)=CC(C)=O", "strong", (237, 241), None, None),
]
FIXTURE_IDS = [f[0] for f in FIXTURE]


@pytest.mark.parametrize("name, smiles, strength, window, jp_detector, jp_wl", FIXTURE, ids=FIXTURE_IDS)
def test_fixture_strength_and_lambda(name, smiles, strength, window, jp_detector, jp_wl):
    a = _a(smiles)
    assert a.strength == strength, (a.strength, [(c.name, c.lam) for c in a.chromophores])
    if window:
        assert window[0] <= a.lam <= window[1], (a.lam, a.lam_method)
    elif strength == "none":
        assert a.lam is None


def test_fixture_size_and_jp_context():
    assert len(FIXTURE) == 38 and sum(f[4] is not None for f in FIXTURE) == 27
    d = [abs(f[5] - _a(f[1]).lam) for f in FIXTURE if f[4] == "UV" and f[5] and f[5] >= 230 and f[3]]
    assert statistics.median(d) <= 25      # JP wavelengths are the ones chosen for sensitivity: context only, loose bound


# ---------------------------------------------------------------- class table
def test_class_table_is_wellformed():
    names = [c[0] for c in R.CHROMOPHORES]
    assert len(names) == len(set(names))
    assert all(c[3] in (R.STRONG, R.WEAK) and c[2] > 0 for c in R.CHROMOPHORES)
    assert len(R._PATTERNS) == len(R.CHROMOPHORES) and all(p is not None for ps in R._PATTERNS for p in ps)
    assert R.WHOLE_RING <= set(names)


@pytest.mark.parametrize("smiles, name", [
    ("c1ccccc1", "Benzene ring"), ("Oc1ccccc1", "Phenol / aryl ether"), ("Nc1ccccc1", "Aniline"), ("CC(=O)Nc1ccccc1", "Anilide"),
    ("[O-][N+](=O)c1ccccc1", "Nitroarene"), ("O=C(c1ccccc1)c1ccccc1", "Diaryl ketone"), ("NC(=O)c1ccccc1", "Benzamide"), ("N#Cc1ccccc1", "Benzonitrile"),
    ("C=Cc1ccccc1", "Styrene (aryl alkene)"), ("c1ccccc1/C=C/c1ccccc1", "Stilbene"), ("O=C(/C=C/c1ccccc1)c1ccccc1", "Chalcone"),
    ("c1ccc2ccccc2c1", "Naphthalene"), ("c1ccc2cc3ccccc3cc2c1", "Anthracene"), ("c1ccc2nccc(c2c1)", "Quinoline / isoquinoline / diazanaphthalene"),
    ("c1ccc2[nH]ccc2c1", "Indole"), ("c1ccncc1", "Pyridine"), ("O=c1cc[nH]c(=O)[nH]1", "Pyrimidine (uracil, cytosine, thymine)"),
    ("Nc1ncnc2[nH]cnc12", "Purine (adenine, guanine)"), ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "Xanthine (caffeine, theophylline)"),
    ("c1ncc[nH]1", "Imidazole"), ("c1ccsc1", "Thiophene"), ("O=c1ccc2ccccc2o1", "Coumarin"), ("c1ccc(cc1)N=Nc1ccccc1", "Azobenzene (aryl azo)"),
    ("O=C1C=CC(=O)C=C1", "Quinone (benzo- / naphtho- / anthraquinone)"), ("CC(C)=O", "Ketone (isolated C=O)"), ("CC=O", "Aldehyde (isolated C=O)"),
    ("CC(=O)OC", "Acid / ester / lactone (isolated C=O)"), ("CC(=O)N", "Amide / urea / carbamate (isolated C=O)"), ("CCC=C", "Isolated C=C"),
])
def test_class_hits(smiles, name):
    names = [c.name for c in _a(smiles).chromophores]
    assert name in names, names


def test_weak_azole_rings_are_weak_alone_but_strong_when_fused_or_conjugated():
    assert _a("c1ncc[nH]1").strength == "weak" and _a("c1ccoc1").strength == "weak"
    assert _a("c1ccc2[nH]cnc2c1").strength == "strong"                 # benzimidazole
    assert _a("c1ccccc1-c1ncc[nH]1").strength == "strong"              # conjugated with benzene
    assert _a("c1cscn1").strength == "strong"                          # thiazole


def test_strength_follows_the_conjugation_rule():
    assert _a("CCCCCC").strength == "none" and _a("OCC(O)CO").strength == "none"
    assert _a("CCCCC=C").strength == "weak"                            # isolated C=C
    assert _a("CC=CC=CC").strength == "strong"                         # 2 conjugated C=C
    assert _a("CC(=O)C(C)=O").strength == "weak"                       # 1,2-dicarbonyl is not read as conjugated
    assert _a("CC#N").strength == "none"                               # nitrile absorbs below 200 nm only
    assert _a("COC(=O)C=CC").strength == "strong"                      # unsaturated ester


def test_open_chain_sugar_is_not_a_carbonyl():
    assert _a(GLUCOSE_OPEN).strength == "none" == _a(GLUCOSE_RING).strength
    assert _a("CCC=O").strength == "weak" and _a("OCC(=O)CO").strength == "weak"
    note = _a(GLUCOSE_OPEN).notes[0]
    assert "hemiacetal" in note


# ---------------------------------------------------------------- Woodward-Fieser / Fieser-Kuhn / Scott
@pytest.mark.parametrize("smiles, expected", [
    ("C=CC=C", 214), ("C=CC(C)=C", 219), ("C1=CC=CCC1", 263), ("CC=CC=CC", 224),
    ("O=C1CCCC=C1", 227), ("O=C1CCCC(C)=C1", 239), ("O=C1CCC=C1", 214), ("O=CC=C", 207), ("CC=CC=O", 219),
    ("CC(C)=CC(C)=O", 239), ("COC(=O)C=CC", 207), ("COC(=O)C(C)=C", 205),
])
def test_woodward_fieser(smiles, expected):
    a = _a(smiles)
    assert a.lam_method == "Woodward–Fieser" and a.lam == pytest.approx(expected, abs=0.01)


def test_woodward_fieser_exocyclic_and_ring_residues_steroid_enones():
    t = _a("C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C")       # testosterone: 215 + 2 x 12 (beta ring residues) + 5 (exocyclic)
    assert t.lam == pytest.approx(244) and t.lam_method == "Woodward–Fieser"
    p = _a("C[C@]12C=CC(=O)C=C1CC[C@@H]3[C@@H]2[C@H](C[C@]4([C@H]3CC[C@@]4(C(=O)CO)O)C)O")   # 1,4-dien-3-one
    assert p.lam == pytest.approx(244)


def test_woodward_fieser_does_not_run_on_macrocycle_rings_or_beta_amino_esters():
    macro = _a("C[C@@H]1C/C=C/C=C/[C@@H](O)CC(=O)O1")       # diene in a macrolactone: no homoannular term
    assert macro.strength == "strong"
    dhp = _a("CCOC(=O)C1=C(NC(=C(C1c1ccccc1)C(=O)OC)C)C")    # beta-amino acrylate: class table (1,4-dihydropyridine), not Woodward-Fieser
    assert any(c.name.startswith("1,4-Dihydropyridine") and c.method == "class table" for c in dhp.chromophores)
    assert not any(c.method == "Woodward–Fieser" for c in dhp.chromophores)


def test_fieser_kuhn_for_long_polyenes():
    bc = _a(BETA_CAROTENE)
    assert bc.lam_method == "Fieser–Kuhn" and bc.lam == pytest.approx(453.3, abs=0.5)
    retinol = _a("CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/CO)/C)/C")
    assert retinol.lam_method == "Fieser–Kuhn" and 320 <= retinol.lam <= 330        # observed 325


@pytest.mark.parametrize("smiles, expected", [
    ("CC(=O)c1ccccc1", 246), ("O=Cc1ccccc1", 250), ("OC(=O)c1ccccc1", 230), ("Nc1ccc(cc1)C(=O)O", 288),
    ("CC(=O)c1ccc(O)cc1", 271), ("COc1cc(C=O)ccc1O", 282), ("CC(=O)c1ccc(Cl)cc1", 256), ("CC(=O)c1ccccc1C", 249),
    ("O=C1CCCc2ccccc12", 246 + 3), ("COC(=O)c1ccc(O)cc1", 255), ("OC(=O)c1ccccc1O", 237),
])
def test_scott(smiles, expected):
    a = _a(smiles)
    assert a.lam_method == "Scott" and a.lam == pytest.approx(expected, abs=0.01)


def test_scott_skips_amides_and_diaryl_ketones():
    assert _a("NC(=O)c1ccccc1").lam_method == "class table"
    assert _a("O=C(c1ccccc1)c1ccccc1").lam_method == "class table"


def test_headline_is_the_largest_lambda_of_the_strongest_tier():
    a = _a("CC(=O)CCc1ccc(O)cc1")        # phenol (strong, 270) beats the isolated ketone (weak, 280)
    assert a.strength == "strong" and a.lam == 270
    a = _a("CC(=O)CCCC(=O)OC")           # only weak groups: ketone 280 is the largest
    assert a.strength == "weak" and a.lam == 280


def test_weak_chromophore_below_the_usable_window_counts_as_none():
    a = _a("CCCC#N")
    assert a.strength == "none" and a.lam is None
    assert "not usable" in a.detectors[0][2]


def test_generic_fallback_for_unlisted_systems_and_extra_conjugation():
    a = _a("c1ccc2ncc3ccccc3c2c1")                       # phenanthridine: no own entry beyond quinoline
    assert a.strength == "strong"
    styrene_like = _a("C=CC=Cc1ccccc1")                  # one conjugated C=C is not covered by a hit: +25
    assert styrene_like.lam >= 255


# ---------------------------------------------------------------- fluorescence
@pytest.mark.parametrize("smiles, level", [
    ("c1ccc2cc3ccccc3cc2c1", "likely"), ("c1ccc2c(c1)ccc1ccccc12", "likely"), ("c1ccc2nc3ccccc3cc2c1", "likely"), ("c1ccc2c(c1)[nH]c1ccccc12", "likely"),
    ("C1c2ccccc2-c2ccccc12", "likely"), ("O=c1ccc2ccccc2o1", "likely"), ("O=c1ccc2ccc(O)cc2o1", "likely"), ("c1ccc2[nH]ccc2c1", "likely"),
    ("COc1ccc2nccc(C)c2c1", "likely"), ("Oc1cccc2ccccc12", "likely"), ("O=c1cc[nH]c2ccccc12", "likely"), ("O=C1OC2(c3ccccc13)c1ccc(O)cc1Oc1cc(O)ccc12", "likely"),
    ("c1ccc2ccccc2c1", "possible"), ("c1ccc2nccc(C)c2c1", "possible"), ("Oc1ccccc1", "possible"), ("Nc1ccccc1", "possible"), ("O=c1cc(-c2ccccc2)oc2ccccc12", "possible"),
    ("c1ccc2[nH]cnc2c1", "possible"), ("c1ccccc1-c1ccccc1", "possible"),
    ("c1ccccc1", "unlikely"), ("Cc1ccccc1", "unlikely"), (CAFFEINE, "unlikely"), ("c1ccncc1", "unlikely"), ("CCCCCC", "unlikely"), (BETA_CAROTENE, "unlikely"),
    ("[O-][N+](=O)c1cccc2ccccc12", "unlikely"), ("Brc1cccc2ccccc12", "unlikely"), ("Ic1ccc2cc3ccccc3cc2c1", "unlikely"), ("O=c1ccc2ccc([N+](=O)[O-])cc2o1", "unlikely"),
    ("c1ccc(cc1)N=Nc1ccc2ccccc2c1", "unlikely"),
])
def test_fluorescence_rules(smiles, level):
    a = _a(smiles)
    assert a.fluor == level, (a.fluor, a.fluor_reason)


def test_conjugated_ketone_lowers_fluorescence_one_step():
    anthracene = _a("c1ccc2cc3ccccc3cc2c1")
    anthracenone = _a("O=C(C)c1ccc2cc3ccccc3cc2c1")          # 2-acetylanthracene
    assert anthracene.fluor == "likely" and anthracenone.fluor == "possible"
    assert "ketone" in anthracenone.fluor_reason


def test_quencher_reason_is_reported():
    a = _a("[O-][N+](=O)c1cccc2ccccc12")
    assert a.fluor == "unlikely" and "nitro" in a.fluor_reason


# ---------------------------------------------------------------- detectors
def _rows(a):
    return {r[0]: r for r in a.detectors}


def test_detector_table_shape_and_values():
    rows = _a(CAFFEINE).detectors
    assert [r[0] for r in rows] == ["UV/PDA", "FLD", "CAD", "ELSD", "RI"]
    assert all(len(r) == 3 and r[1] in ("good", "possible", "poor") and r[2] for r in rows)


def test_uv_and_fld_follow_strength_and_fluorescence():
    caffeine = _rows(_a(CAFFEINE))
    assert caffeine["UV/PDA"][1] == "good" and caffeine["FLD"][1] == "poor"
    anthracene = _rows(_a("c1ccc2cc3ccccc3cc2c1"))
    assert anthracene["FLD"][1] == "good"
    naphthalene = _rows(_a("c1ccc2ccccc2c1"))
    assert naphthalene["FLD"][1] == "possible"
    erythromycin = _rows(_a(next(f[1] for f in FIXTURE if f[0] == "erythromycin")))
    assert erythromycin["UV/PDA"][1] == "possible"
    assert _rows(_a(GLUCOSE_RING))["UV/PDA"][1] == "poor"


def test_cad_and_elsd_use_the_gc_verdict_and_mw():
    caffeine = _a(CAFFEINE)
    assert caffeine.gc == "yes"
    assert _rows(caffeine)["CAD"][1] == "poor" and _rows(caffeine)["ELSD"][1] == "poor"
    glucose = _a(GLUCOSE_RING)
    assert glucose.gc != "yes" and glucose.mw >= 150
    assert _rows(glucose)["CAD"][1] == "good" and _rows(glucose)["ELSD"][1] == "good"
    assert "volatile" in _rows(glucose)["CAD"][2]                      # the mobile phase must be volatile
    glycine = _a("NCC(=O)O")
    assert glycine.mw < 150 and _rows(glycine)["CAD"][1] == "possible"
    ionic = _a("CCCCCCCCCCCC[N+](C)(C)Cc1ccccc1")
    assert _rows(ionic)["CAD"][1] == "good"


def test_cad_and_elsd_follow_the_corrected_gc_verdict_for_phosphate_esters():
    """Triethyl phosphate / triphenyl phosphate used to be 'GC no' (every P(=O)(O)O ester matched), so CAD / ELSD were 'good'; as GC-volatile compounds they are
    treated as volatile like any other (the CAD / ELSD rule itself is unchanged). A phosphate with a free P-OH is still not volatile."""
    for smiles in ("CCOP(=O)(OCC)OCC", "O=P(Oc1ccccc1)(Oc1ccccc1)Oc1ccccc1"):
        a = _a(smiles)
        assert a.gc == "yes" and _rows(a)["CAD"][1] == "poor" and _rows(a)["ELSD"][1] == "poor" and "GC verdict yes" in _rows(a)["CAD"][2], smiles
    g6p = _a("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1OP(=O)(O)O")
    assert g6p.gc == "no" and _rows(g6p)["CAD"][1] == "good"


def test_ri_is_always_possible_and_never_good():
    for s in (CAFFEINE, GLUCOSE_RING, "CCCCCC", BETA_CAROTENE):
        assert _rows(_a(s))["RI"][1] == "possible"
    assert "isocratic" in _rows(_a(GLUCOSE_RING))["RI"][2]


def test_mobile_phase_changes_the_low_uv_advice():
    glycine_acn = _a("NCC(=O)O", ACN)
    glycine_meoh = _a("NCC(=O)O", MEOH)
    assert _rows(glycine_acn)["UV/PDA"][1] == "possible" and _rows(glycine_meoh)["UV/PDA"][1] == "poor"
    assert any("Low-UV rule" in n for n in glycine_meoh.notes) and any("phosphate" in n for n in glycine_acn.notes)
    # a strong chromophore near the cutoff: fine with acetonitrile, only possible with methanol / formic acid
    ester = "COC(=O)C(C)=CC"
    assert _rows(_a(ester, ACN))["UV/PDA"][1] != _rows(_a(ester, MEOH))["UV/PDA"][1] or _a(ester, ACN).lam < 210
    # a strong benzene chromophore is unaffected
    assert _rows(_a("c1ccccc1", ACN))["UV/PDA"][1] == _rows(_a("c1ccccc1", MEOH))["UV/PDA"][1] == "good"


def test_derivatization_pointers_only_without_a_strong_chromophore():
    notes = " ".join(_a("NCC(=O)O").notes)
    assert "OPA" in notes and "FMOC" in notes and "Derivatization page" in notes
    assert "DNPH" in " ".join(_a("CC(C)CC(=O)C").notes)
    assert "no group for OPA / FMOC / DNPH" in " ".join(_a(GLUCOSE_RING).notes)
    assert "Derivatization option" not in " ".join(_a("c1ccccc1N").notes)          # strong chromophore: no hint
    assert "OPA" in _rows(_a("NCC(=O)O"))["FLD"][2]


def test_notes_carry_the_rules_and_no_general_caveats():
    notes = " ".join(_a(CAFFEINE).notes)
    assert "±10–20 nm" in notes and "UV rule" in notes and "λmax rule" in notes and "CAD / ELSD rule" in notes
    assert "Caveat" not in notes and "solvent and pH" not in notes


# ---------------------------------------------------------------- structures and errors
def test_salt_is_reduced_to_its_largest_part_and_keeps_input_atom_numbers():
    a = _a("[Na+].[O-]C(=O)c1ccccc1")
    assert any("disconnected" in w for w in a.warnings)
    assert a.strength == "strong" and a.lam == 230
    atoms = {i for c in a.chromophores for i in c.atoms}
    assert 0 not in atoms and atoms >= {1, 2, 3, 4, 5, 6, 7, 8, 9}       # Na is atom 1 of the input; the benzoate carries on from atom 2
    b = _a("c1ccccc1C(=O)[O-].[Na+]")
    assert max(i for c in b.chromophores for i in c.atoms) <= 8


def test_charge_neutralization_warning():
    assert any("neutralized" in w for w in _a("CC(=O)[O-]").warnings)


def test_inorganic_and_tiny_inputs_do_not_fail():
    for s in ("O", "[Na+].[Cl-]", "[H][H]", "C", "CS(C)=O"):
        a = _a(s)
        assert a.strength in ("none", "weak") and len(a.detectors) == 5


def test_errors():
    with pytest.raises(R.NonMsError):
        _a("not a smiles")
    with pytest.raises(R.NonMsError):
        _a("c1ccccc1", "Water/ethanol")


def test_figure_has_the_chromophore_atoms_highlighted():
    a = _a("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
    assert a.svg.startswith("<?xml") or "<svg" in a.svg
    hexes = ["#%02X%02X%02X" % tuple(int(c * 255) for c in rgb) for rgb in (R.STRONG_COLOR, R.WEAK_COLOR)]
    assert all(h in a.svg for h in hexes)                       # benzene ring (strong) and the acid C=O (weak) are both highlighted
    assert _a("CCCCCC").svg and not any(h in _a("CCCCCC").svg for h in hexes)


# ---------------------------------------------------------------- API
def test_spec_matches_engine_and_runner():
    mod = next(m for m in MODULES if m["id"] == "non_ms")
    assert mod["live"] and mod["example"]["structure"]["value"] == "quinine"
    phase = next(f for f in mod["inputs"] if f["id"] == "mobile_phase")
    assert phase["options"] == list(R.MOBILE_PHASES) and phase["value"] == R.DEFAULT_MOBILE_PHASE
    ids = {b["id"]: b for b in mod["outputs"]}
    assert [b["id"] for b in mod["outputs"]] == ["detectors", "chromophore_structure", "rule_hits"]          # no separate optical cards
    assert ids["detectors"]["type"] == "portal" and ids["detectors"]["title"] == "Detectors" and ids["detectors"]["help"]
    assert ids["chromophore_structure"]["type"] == "figure" and ids["chromophore_structure"]["help"]
    assert ids["rule_hits"]["type"] == "notes" and ids["rule_hits"]["collapsed"] is True


def test_api_fills_every_block():
    mod = next(m for m in MODULES if m["id"] == "non_ms")
    body = client.post("/api/modules/non_ms/run", json=mod["example"]).json()
    assert body["status"] == "ok"
    out = body["outputs"]
    assert set(out) == {b["id"] for b in mod["outputs"]}
    assert set(out) == {"detectors", "chromophore_structure", "rule_hits"}
    tiles = out["detectors"]["tiles"]
    assert [t["label"] for t in tiles] == ["UV/PDA", "FLD", "CAD", "ELSD", "RI"]
    assert [t["status"] for t in tiles] == ["yes", "yes", "yes", "yes", "partial"]        # quinine: UV and native fluorescence; non-volatile for CAD / ELSD; RI possible
    assert tiles[0]["value"] == "strong · ~315 nm" and tiles[1]["value"] == "likely" and all(t["help"] for t in tiles)
    assert "Quinoline" in tiles[0]["help"][0] and tiles[0]["help"][1] == "λmax by class table: ±10–20 nm"      # class and method in the UV help
    assert len(tiles[1]["help"]) == 1 and "quinoline-type" in tiles[1]["help"][0]                                                        # the FLD help is the reason only
    fig = out["chromophore_structure"]
    assert fig["svg"] and fig["help"][0].startswith("Quinoline") and "atoms 3–12" in fig["help"][0] and "~315 nm" in fig["help"][0]
    assert out["rule_hits"]["items"]


def test_api_woodward_fieser_method_in_the_uv_tile_help():
    uv = _run("CC(C)=CC(C)=O").json()["outputs"]["detectors"]["tiles"][0]
    assert uv["value"] == "strong · ~239 nm" and any(h.startswith("λmax by Woodward–Fieser") for h in uv["help"])


def test_api_no_chromophore_gives_a_none_tile_and_a_plain_figure():
    body = _run(GLUCOSE_RING).json()
    out = body["outputs"]
    assert out["chromophore_structure"]["svg"] and out["chromophore_structure"]["help"] == "no chromophore found"
    uv = out["detectors"]["tiles"][0]
    assert (uv["status"], uv["value"]) == ("no", "none") and not any(h.startswith("λmax") for h in uv["help"])


def test_api_mobile_phase_input():
    glycine = "NCC(=O)O"
    acn = _run(glycine, mobile_phase=ACN).json()["outputs"]["detectors"]["tiles"][0]
    meoh = _run(glycine, mobile_phase=MEOH).json()["outputs"]["detectors"]["tiles"][0]
    assert (acn["status"], meoh["status"]) == ("partial", "no")
    assert _run(glycine).json()["outputs"]["detectors"]["tiles"][0]["status"] == "partial"       # default = acetonitrile / phosphate
    assert _run(glycine, mobile_phase="Water/ethanol").status_code == 400


def _tiles(smiles):
    return {t["label"]: t for t in _run(smiles).json()["outputs"]["detectors"]["tiles"]}


def test_api_tiles_follow_the_engine_rows_and_values():
    for smiles in (CAFFEINE, GLUCOSE_RING, "NCC(=O)O", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"):
        a = _a(smiles)
        tiles = _tiles(smiles)
        assert [tiles[r[0]]["status"] for r in a.detectors] == [{"good": "yes", "possible": "partial", "poor": "no"}[r[1]] for r in a.detectors]
        assert [tiles[r[0]]["help"][0] for r in a.detectors] == [r[2] for r in a.detectors]       # the reason is the tile's help
    assert _tiles(CAFFEINE)["CAD"]["value"] == "volatile" and _tiles(GLUCOSE_RING)["CAD"]["value"] == "non-volatile"
    assert _tiles("NCC(=O)O")["ELSD"]["value"] == "small, non-volatile" and _tiles(CAFFEINE)["RI"]["value"] == "isocratic only"


def test_api_derivatization_pointers_sit_in_the_uv_and_fld_tile_help():
    glycine = _tiles("NCC(=O)O")                                                          # weak chromophore, primary amine: OPA / FMOC
    assert any("OPA" in h and "FMOC" in h for h in glycine["UV/PDA"]["help"]) and any("OPA" in h for h in glycine["FLD"]["help"])
    assert not any("Derivatization adds" in h for t in _tiles(CAFFEINE).values() for h in t["help"])     # strong chromophore: no pointer
    assert not any("Derivatization adds" in h for h in _tiles(GLUCOSE_RING)["UV/PDA"]["help"])           # no amine / carbonyl: nothing to point to


def test_api_uv_and_fld_tile_values_for_other_cases():
    t = _tiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
    assert t["UV/PDA"]["value"] == "strong · ~255 nm" and "Benzene ring" in t["UV/PDA"]["help"][0]
    assert t["FLD"]["value"] == "unlikely" and "benzene ring without a donor or a fused fluorophore" in t["FLD"]["help"][0]
    assert _tiles("NCC(=O)O")["UV/PDA"]["value"] == "weak · end absorption"                  # glycine: only the C=O band, at the end absorption
    assert _tiles("CC(C)=CC(C)=O")["FLD"]["value"] == "unlikely"
    lido = _tiles("CCN(CC)CC(=O)Nc1c(C)cccc1C")
    assert lido["FLD"]["value"] == "possible" and lido["FLD"]["status"] == "partial"


def test_atom_list_is_compressed_to_ranges():
    from mzdesk.api.runners.non_ms import _atoms
    assert _atoms([1, 2, 3, 4, 6, 7, 9]) == "2–5, 7, 8, 10" and _atoms([0]) == "1" and _atoms([5, 3, 4]) == "4–6"


def test_api_errors_and_warnings():
    assert _run("not a smiles").status_code == 400
    assert client.post("/api/modules/non_ms/run", json={}).status_code == 400
    body = _run("Cl.CN(C)CCCN1c2ccccc2Sc2ccc(Cl)cc21").json()
    assert any("disconnected" in w for w in body["warnings"])


def test_api_name_input_uses_the_resolver(monkeypatch):
    from mzdesk.identity import resolve
    monkeypatch.setattr(resolve, "to_smiles", lambda s: CAFFEINE)
    body = client.post("/api/modules/non_ms/run", json={"structure": {"mode": "Name / ID / Formula", "value": "caffeine"}}).json()
    assert body["outputs"]["detectors"]["tiles"][0]["value"] == "strong · ~273 nm"


# ---------------------------------------------------------------- JP18 regression (coarse thresholds, same data as validation/non_ms/jp18_validate.py)
@pytest.mark.skipif(not JP_FILE.is_file(), reason="JP18 dataset not present")
def test_jp18_dataset_regression():
    import re

    recs = [json.loads(line) for line in gzip.open(JP_FILE, "rt", encoding="utf-8")]
    cache = {}

    def run(smiles):
        if smiles not in cache:
            cache[smiles] = R.analyze(smiles)
        return cache[smiles]

    uv, universal = {}, {}
    for r in recs:
        if r["technique"] != "LC":
            continue
        text = r["condition_text"] or ""
        derivatized = bool(re.search(r"reaction coil|post-column|derivati[sz]|ninhydrin|o-phthal", text, re.I))
        for s in r["analyte_structures"] or []:
            smi = s.get("smiles")
            if not smi or not re.search(r"[Cc]", re.sub(r"Cl|Ca|Cu|Co|Cr|Cd|Cs|Ce", "", smi)):
                continue
            if r["detector_type"] in ("UV", "PDA") and not derivatized:
                uv[(smi, r["wavelength_nm"])] = r["wavelength_nm"]
            elif r["detector_type"] == "RI" or derivatized or (r["detector_type"] == "visible" and r["wavelength_nm"] == 570):
                universal[smi] = r["detector_type"]
    uv_smiles = {k[0] for k in uv}
    with_chromophore = sum(run(s).strength != "none" for s in uv_smiles) / len(uv_smiles)
    assert len(uv_smiles) > 500 and with_chromophore >= 0.95                           # measured: 98.7 %
    none_or_weak = sum(run(s).strength in ("none", "weak") for s in universal) / len(universal)
    assert len(universal) >= 25 and none_or_weak >= 0.9                                # measured: 100 %
    d = [abs(wl - run(smi).lam) for (smi, wl) in uv if wl >= 230 and run(smi).lam]
    assert len(d) > 500 and statistics.median(d) <= 20                                 # measured: 13 nm (unique analyte x wavelength)
