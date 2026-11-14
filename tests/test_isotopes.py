"""同位素引擎（mzdesk/ms1/isotopes.py）：与加合物算术一致、分辨率合并、M−1（硼）、余弦相似度、Theoretical / Mass search 页。"""
import pytest
from fastapi.testclient import TestClient

from mzdesk.api import verdicts
from mzdesk.api.app import app
from mzdesk.api.runners.ms1 import run_search
from mzdesk.ms1 import isotopes as iso
from mzdesk.ms1 import msbuddy_engine as bd
from mzdesk.spec import ADDUCTS_NEG, ADDUCTS_POS, CAFFEINE_MS1, METFORMIN_MS1, MODULES, PIPELINES

client = TestClient(app)
ALL = ADDUCTS_POS[1:] + ADDUCTS_NEG[1:]


def by_offset(p):
    out = {}
    for x in p["peaks"]:
        out[x["offset"]] = out.get(x["offset"], 0.0) + x["abundance"]
    return out


@pytest.mark.parametrize("adduct", ALL)
def test_monoisotopic_mz_is_the_adduct_arithmetic(adduct):
    for formula in ("C8H10N4O2", "C42H82NO8P"):
        mz, ion = bd.ion_mz_and_formula(formula, adduct)
        p = iso.pattern(formula, adduct, "Unit mass")
        assert p["mono_mz"] == pytest.approx(mz, rel=1e-7)          # 0.1 ppm
        assert p["ion_formula"] == ion
        top = next(x for x in p["peaks"] if x["offset"] == 0)
        assert top["mz"] == pytest.approx(mz, abs=0.02)


def test_caffeine_m_plus_1():
    # 手算：8*0.0107/0.9893 + 4*0.00364/0.99636 + 11*0.000115 + 2*0.00038/0.99762 ≈ 0.103
    for res in ("Unit mass", "Q-TOF 30k", "Orbitrap 60k"):
        p = iso.pattern("C8H10N4O2", "[M+H]+", res)
        assert 9.0 <= by_offset(p)[1] <= 11.0 and by_offset(p)[1] == pytest.approx(10.3, abs=0.5)
    import IsoSpecPy
    d = IsoSpecPy.IsoTotalProb(formula="C8H11N4O2", prob_to_cover=0.99999, get_confs=False)
    m1 = sum(pr for m, pr in zip(d.masses, d.probs) if 0.5 < m - min(d.masses) < 1.5)
    assert by_offset(iso.pattern("C8H10N4O2", "[M+H]+", "Unit mass"))[1] == pytest.approx(m1 / max(d.probs) * 100, rel=1e-3)


def test_resolution_merging():
    q, o, u = (iso.pattern("C8H10N4O2", "[M+H]+", r) for r in ("Q-TOF 30k", "Orbitrap 60k", "Unit mass"))
    assert [x["label"] for x in u["peaks"]] == [x["label"] for x in q["peaks"]] == ["M", "M+1", "M+2"]
    assert len(o["peaks"]) > len(q["peaks"])                              # 60k 分开 13C 和 15N
    assert [x["label"] for x in o["peaks"] if x["offset"] == 1] == ["M+1 (a)", "M+1 (b)"]
    assert max(x["abundance"] for x in u["peaks"]) == 100 and min(x["abundance"] for x in u["peaks"]) >= 0.1


@pytest.mark.parametrize("formula,adduct,m2", [("C6H5Cl", "[M]+", (30.0, 34.0)), ("C6H5Br", "[M]+", (93.0, 103.0))])
def test_halogen_m_plus_2(formula, adduct, m2):
    p = by_offset(iso.pattern(formula, adduct, "Unit mass"))
    assert m2[0] <= p[2] / p[0] * 100 <= m2[1]        # 35Cl/37Cl ≈ 32 %，79Br/81Br ≈ 97 %


def test_boron_has_an_m_minus_1_peak():
    # 氯马司特 C14H10BNO3：10B（19.9 %）在 M−1，约为 11B 峰的 25 %；M = 11B（与 PubChem 单同位素质量一致）
    p = iso.pattern("C14H10BNO3", "[M+H]+", "Orbitrap 60k")
    labels = {x["label"]: x for x in p["peaks"]}
    assert p["ion_formula"] == "C14H11BNO3+" and p["mono_mz"] == pytest.approx(251.0753733 + 1.00727646 + 0.0000, abs=2e-3)
    assert labels["M−1"]["abundance"] == pytest.approx(24.9, abs=1.0) and labels["M−1"]["mz"] < p["mono_mz"]
    assert max(p["peaks"], key=lambda x: x["abundance"])["label"].startswith("M ")
    assert iso.pattern("C14H10BNO3", "[M+H]+", "Unit mass")["peaks"][0]["label"] == "M−1"


def test_charged_formula_and_errors():
    p = iso.pattern("C5H14NO+", "[M]+", "Unit mass")
    assert p["ion_formula"] == "C5H14NO+" and p["mono_mz"] == pytest.approx(104.1070, abs=1e-4)
    for bad, adduct in [("C6H5D", "[M+H]+"), ("C8H10N4O2[13C]", "[M+H]+"), ("Xx2", "[M+H]+"), ("", "[M+H]+"),
                        ("CH4", "[M-2H2O+H]+"), ("C8H10N4O2", "[M+Foo]+")]:
        with pytest.raises(ValueError):
            iso.pattern(bad, adduct)
    with pytest.raises(ValueError):
        iso.pattern("C8H10N4O2", "[M+H]+", "R=1")


def test_similarity_is_a_cosine_with_unmatched_as_zero():
    p = iso.pattern("C8H10N4O2", "[M+H]+", "Q-TOF 30k")["peaks"]
    perfect = [(x["mz"], x["abundance"]) for x in p]
    assert iso.similarity(p, perfect, "Q-TOF 30k")["score"] == pytest.approx(1.0)
    assert iso.similarity(p, [(x["mz"] + 0.0001, 5 * x["abundance"]) for x in p], "Q-TOF 30k")["score"] == pytest.approx(1.0)   # 与强度尺度无关
    only_m = iso.similarity(p, perfect[:1], "Q-TOF 30k")
    assert only_m["score"] == pytest.approx(100 / (100 ** 2 + 10.39 ** 2 + 0.9 ** 2) ** 0.5, abs=1e-3)
    assert [r["obs_mz"] is None for r in only_m["rows"]] == [False, True, True]
    junk = perfect + [(300.0, 1e6), (195.0877 + 0.0004, 1)]            # 无匹配的峰被忽略；每个观测峰只用一次
    assert iso.similarity(p, junk, "Q-TOF 30k")["score"] == pytest.approx(1.0)
    assert iso.similarity(p, [(500.0, 100)], "Q-TOF 30k")["score"] == 0.0
    wrong = iso.similarity(p, [(x["mz"], v) for x, v in zip(p, (100, 80, 60))], "Q-TOF 30k")["score"]
    assert wrong < 0.9
    assert iso.similarity(p, [(195.30, 100)], "Unit mass")["score"] > 0.9 and iso.similarity(p, [(195.30, 100)], "Q-TOF 30k")["score"] == 0.0


# ---------------------------------------------------------------- pages
def _theory(**kw):
    body = {"structure": {"mode": "SMILES / InChI", "value": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"},
            "precursor": {"precursor_mz": "", "ion_mode": "Positive", "adduct": "[M+H]+"}, "resolution": "Q-TOF 30k", **kw}
    return client.post("/api/modules/ms1_theory/run", json=body)


def test_theory_page_outputs():
    r = _theory()
    assert r.status_code == 200, r.text
    body = r.json()
    o = body["outputs"]
    assert o["theory"]["values"]["Precursor m/z"] == {"value": "195.0877", "note": "[M+H]+"}
    assert o["theory"]["values"]["Ion formula"] == "C8H11N4O2+"
    assert o["theory"]["values"]["Most abundant peak"]["value"] == "M"
    assert [r[0] for r in o["pattern"]["rows"]] == ["M", "M+1", "M+2"]
    fig = o["pattern_plot"]
    assert fig["plotly"]["data"] and [c["kind"] for c in fig["peaks"]] == ["M", "M+1", "M+2"]
    assert "adducts" not in o and body["verdict"].startswith("<b>[M+H]+</b> at <b>m/z 195.0877</b> (C8H11N4O2+); 3 isotope peaks")


def test_theory_example_is_crisaborole_and_offline():
    m = next(m for m in MODULES if m["id"] == "ms1_theory")
    assert m["example"]["structure"]["value"] == "crisaborole" and m["example"]["precursor"]["adduct"] == "[M+H]+"
    pf = next(f for f in m["inputs"] if f["type"] == "precursor")
    assert pf["mz"] is False and all("Auto" not in v for v in pf["adducts"].values())
    r = client.post("/api/modules/ms1_theory/run", json=m["example"])
    assert r.status_code == 200, r.text
    o = r.json()["outputs"]
    assert o["theory"]["values"]["Ion formula"] == "C14H11BNO3+" and o["theory"]["values"]["Precursor m/z"]["value"] == "252.0826"
    assert o["pattern"]["rows"][0][0] == "M−1"


def test_theory_accepts_any_element_and_rejects_bad_input():
    r = _theory(structure={"mode": "Name / ID / Formula", "value": "C6H5Cl"}, precursor={"precursor_mz": "", "ion_mode": "Positive", "adduct": "[M+H]+"})
    assert r.status_code == 200 and r.json()["outputs"]["theory"]["values"]["Ion formula"] == "C6H6Cl+"
    assert _theory(resolution="nope").status_code == 400
    assert _theory(precursor={"precursor_mz": "", "ion_mode": "Positive", "adduct": "[M-H]-"}).status_code == 400
    assert _theory(structure={"mode": "Name / ID / Formula", "value": "C2H6Xx"}).status_code in (400, 502)


def test_theory_auto_from_old_inputs_means_the_default_adduct():
    r = _theory(precursor={"precursor_mz": "", "ion_mode": "Negative", "adduct": "Auto"})
    assert r.status_code == 200 and r.json()["outputs"]["theory"]["values"]["Ion formula"] == "C8H9N4O2-"


def test_isotope_match_page_is_gone_and_formula_moved():
    ids = [m["id"] for m in MODULES]
    assert "ms1_match" not in ids
    assert [m["id"] for m in MODULES if m["subgroup"] == "ms1"] == ["ms1_search", "ms1_theory"]
    assert [m["id"] for m in MODULES if m["subgroup"] == "ms2"].index("ms1_formula") == [m["id"] for m in MODULES if m["subgroup"] == "ms2"].index("ms2_library") + 1
    assert next(m["title"] for m in MODULES if m["id"] == "ms1_search") == "Mass Search"
    assert next(m["title"] for m in MODULES if m["id"] == "ms1_theory") == "Theoretical"
    assert all("ms1_match" not in s["modules"] for p in PIPELINES for s in p["stages"])


# ---------------------------------------------------------------- Mass search with isotope MS
def _search(**kw):
    body = {"precursor": {"precursor_mz": "195.0877", "ion_mode": "Positive", "adduct": "Auto"}, "tol_ppm": "5",
            "databases": ["HMDB", "PubChem"], **kw}
    return client.post("/api/modules/ms1_search/run", json=body)


def _groups(body):
    return {g["title"]: g for g in body["outputs"]["hits"]["groups"]}


def test_search_without_isotope_ms_has_no_scores():
    body = _search().json()
    g = _groups(body)["C8H10N4O2 · [M+H]+"]
    assert "isotope" not in g["meta"] and "help" not in g and body["warnings"] == []


def test_search_scores_and_keeps_good_groups():
    r = _search(ms1=CAFFEINE_MS1, resolution="Q-TOF 30k", isotope_min="0.8")
    body = r.json()
    g = _groups(body)["C8H10N4O2 · [M+H]+"]
    assert "isotope 1.00" in g["meta"] and g["help"][0].startswith("Isotope score 1.00")
    assert body["warnings"] == [] and body["outputs"]["summary"]["values"]["Compounds"] == "3"


def test_search_hides_groups_below_the_threshold():
    bad = "195.0877 100\n196.0902 80\n197.0922 60"
    body = _search(ms1=bad, resolution="Q-TOF 30k", isotope_min="0.99").json()
    assert body["outputs"]["hits"] == {"groups": [], "count": 0}
    v = body["outputs"]["summary"]["values"]
    assert v["Compounds"]["value"] == "0" and "3 compounds in 1 formula group hidden: isotope score below 0.99" in v["Compounds"]["note"]
    assert v["Closest hit"] == "—"
    assert body["verdict"] == "All 3 compounds within 5 ppm fail the isotope score (≥ 0.99)."
    shown = _search(ms1=bad, resolution="Q-TOF 30k", isotope_min="0.5").json()
    assert shown["outputs"]["summary"]["values"]["Compounds"] == "3"


def test_search_threshold_note_in_verdict_and_unit_mass_tolerance():
    body = _search(ms1=CAFFEINE_MS1, resolution="Unit mass").json()
    assert "isotope 1.00" in _groups(body)["C8H10N4O2 · [M+H]+"]["meta"]
    far = _search(ms1="195.3 100\n196.3 10", resolution="Unit mass", isotope_min="0.5").json()       # 0.5 Da 内
    assert "isotope" in _groups(far)["C8H10N4O2 · [M+H]+"]["meta"]


def test_search_warns_without_precursor_peak_and_bad_ms1():
    body = _search(ms1="300.1 100\n301.1 10", resolution="Q-TOF 30k", isotope_min="0").json()
    assert any("no peak within 5 ppm of the precursor" in w for w in body["warnings"])
    assert "isotope 0.00" in _groups(body)["C8H10N4O2 · [M+H]+"]["meta"]
    body = _search(ms1="not peaks").json()
    assert any("No valid isotope MS peaks" in w for w in body["warnings"]) and body["outputs"]["summary"]["values"]["Compounds"] == "3"
    assert _search(ms1=CAFFEINE_MS1, resolution="nope").status_code == 400
    assert _search(ms1=CAFFEINE_MS1, isotope_min="2").status_code == 400


def test_search_keeps_unscorable_groups_as_na(monkeypatch):
    from mzdesk.ms1 import isotopes
    real = isotopes.pattern

    def pattern(formula, adduct, resolution="Orbitrap 60k"):
        if formula == "C8H10N4O2":
            raise ValueError("isotope labels are not supported")
        return real(formula, adduct, resolution)
    monkeypatch.setattr(isotopes, "pattern", pattern)
    body = _search(ms1=CAFFEINE_MS1, resolution="Q-TOF 30k", isotope_min="0.99").json()
    g = _groups(body)["C8H10N4O2 · [M+H]+"]
    assert g["meta"].endswith("isotope N/A") and "isotope labels" in g["help"]
    assert any("No isotope score for C8H10N4O2 [M+H]+" in w for w in body["warnings"])


def test_search_permanent_cation_scores():
    body = _search(precursor={"precursor_mz": "104.1070", "ion_mode": "Positive", "adduct": "[M]+"},
                   ms1="104.1070 100\n105.1100 6.0", resolution="Unit mass").json()
    g = _groups(body)["C5H14NO+ · [M]+"]
    assert "isotope 1.00" in g["meta"] or "isotope 0.9" in g["meta"]


def test_example_runs_with_isotopes():
    m = next(m for m in MODULES if m["id"] == "ms1_search")
    assert m["example"]["ms1"] == METFORMIN_MS1
    body = client.post("/api/modules/ms1_search/run", json=m["example"]).json()
    assert "isotope 1.00" in _groups(body)["C4H11N5 · [M+H]+"]["meta"]
