"""API 与页面定义的测试：spec 与 runner 对得上、示例能跑、封装结果与直接调引擎一致、错误输入返回 400。

运行：python -m pytest tests
MetFrag 需要 Java；没有 Java 时相关测试跳过。
"""
import shutil

import pytest
from fastapi.testclient import TestClient

from mzdesk.api.app import app
from mzdesk.api.runners import RUNNERS
from mzdesk.spec import MODULES, PIPELINES, SUBGROUPS

client = TestClient(app)
LIVE = [m for m in MODULES if m["live"]]
HAS_JAVA = shutil.which("java") is not None
OPTIONAL_BLOCKS = {"viewer3d"}          # 3D embedding can fail for some structures


def _run(module_id, inputs):
    return client.post(f"/api/modules/{module_id}/run", json=inputs)


def _paste(peaks):
    return {"source": "Paste", "value": peaks}


# ---------------------------------------------------------------- spec
def test_live_modules_and_runners_match():
    assert {m["id"] for m in LIVE} == set(RUNNERS)


def test_module_ids_unique_and_subgroups_known():
    ids = [m["id"] for m in MODULES]
    assert len(ids) == len(set(ids))
    known = {s["id"] for s in SUBGROUPS}
    assert all(m["subgroup"] in known for m in MODULES if m["subgroup"])


def test_ms2_is_split_into_its_own_subgroup():
    ms2 = [m["id"] for m in MODULES if m["subgroup"] == "ms2"]
    assert ms2 == ["ms2_predicted", "ms2_library", "ms2_similarity", "ms2_fragments", "ms2_denovo"]
    assert "ms2" not in {m["id"] for m in MODULES}


def test_pipeline_stages_point_to_existing_modules():
    ids = {m["id"] for m in MODULES}
    for p in PIPELINES:
        for s in p["stages"]:
            assert s["modules"] and set(s["modules"]) <= ids, (p["id"], s["title"])


def test_block_ids_unique_and_examples_match_inputs():
    for m in MODULES:
        block_ids = [b["id"] for b in m["outputs"] if b.get("id")]
        assert len(block_ids) == len(set(block_ids)), m["id"]
        if m["example"]:
            assert set(m["example"]) <= {f["id"] for f in m["inputs"] if f.get("id")}, m["id"]


# ---------------------------------------------------------------- live modules
@pytest.mark.parametrize("mod", LIVE, ids=[m["id"] for m in LIVE])
def test_example_runs_and_fills_every_block(mod):
    if mod["id"] == "ms2_fragments" and not HAS_JAVA:
        pytest.skip("MetFrag needs Java")
    assert mod["example"], f"{mod['id']} has no example"
    r = _run(mod["id"], mod["example"])
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["status"] == "ok"
    expected = {b["id"] for b in mod["outputs"] if b.get("id")} - OPTIONAL_BLOCKS
    assert expected <= set(body["outputs"]), expected - set(body["outputs"])


def test_placeholder_module_answers_placeholder():
    body = _run("ms2_predicted", {}).json()
    assert body["status"] == "placeholder" and body["outputs"] == {}


def test_similarity_matches_engine():
    from mzdesk.ms2.ms2_similarity import compute_all_similarities
    from mzdesk.shared.spectrum import parse_spectrum

    q, ref = "100 1\n150 3\n200 5", "100 2\n151 3\n200 4"
    body = _run("ms2_similarity", {"query": _paste(q), "reference": _paste(ref), "query_precursor": "250",
                                   "reference_precursor": "251", "tol_da": "0.02"}).json()
    direct = compute_all_similarities(parse_spectrum(q), parse_spectrum(ref), 250, 251, tolerance=0.02, use_deep=False)
    api = {row[1]: row[2] for row in body["outputs"]["scores"]["rows"]}
    for r in direct:
        expected = "—" if r["score"] is None else f"{r['score']:.4f}"
        assert api[r["method"]] == expected, r["method"]


@pytest.mark.skipif(not HAS_JAVA, reason="MetFrag needs Java")
def test_metfrag_matches_engine():
    from mzdesk.ms2 import metfrag_engine as mf
    from mzdesk.shared.spectrum import parse_spectrum
    from mzdesk.spec import PG_36_2, PG_36_2_MS2

    body = _run("ms2_fragments", {"structure": {"mode": "SMILES / InChI", "value": PG_36_2},
                                  "precursor": {"precursor_mz": "", "ion_mode": "Negative", "adduct": "[M-H]-"},
                                  "ms2": _paste(PG_36_2_MS2)}).json()
    direct = mf.run_metfrag(PG_36_2, parse_spectrum(PG_36_2_MS2), "[M-H]-", None, 2, 0.002, 10.0)
    api = [(row[0], row[3]) for row in body["outputs"]["peaks"]["rows"]]
    assert api == [(f"{a['mz']:.4f}", a["ion_formula"] or "—") for a in direct["annotation"]]


def test_buddy_matches_engine():
    from mzdesk.ms1 import msbuddy_engine as bd
    from mzdesk.spec import CAFFEINE_MS1, CAFFEINE_MS2
    from mzdesk.shared.spectrum import parse_spectrum

    body = _run("ms1_formula", {"precursor": {"precursor_mz": "195.0877", "ion_mode": "Positive", "adduct": "Auto"},
                                "ms1": CAFFEINE_MS1, "ms2": _paste(CAFFEINE_MS2)}).json()
    direct = bd.annotate_formula(195.0877, True, None, None, "", tuple(parse_spectrum(CAFFEINE_MS2)),
                                 tuple(parse_spectrum(CAFFEINE_MS1)), "orbitrap", False)
    assert [row[1] for row in body["outputs"]["candidates"]["rows"]] == [c["formula"] for c in direct["candidates"]]


def test_theory_handles_permanent_cation():
    # msbuddy reads [M]+ as [M+H]+; the page computes it itself
    body = _run("ms1_theory", {"structure": {"mode": "SMILES / InChI", "value": "C[N+](C)(C)CCO"},
                               "precursor": {"precursor_mz": "104.1070", "ion_mode": "Positive", "adduct": "[M]+"}}).json()
    assert body["outputs"]["theory"]["values"]["Precursor m/z"]["value"] == "104.1070"


# ---------------------------------------------------------------- bad input
@pytest.mark.parametrize("module_id,inputs", [
    ("properties", {"structure": {"mode": "SMILES / InChI", "value": "C1CC"}}),
    ("properties", {"structure": {"mode": "Name / ID / Formula", "value": "caffeine"}}),
    ("ms1_formula", {"precursor": {"precursor_mz": "2000", "ion_mode": "Positive", "adduct": "Auto"}}),
    ("ms2_similarity", {"query": _paste("100 1"), "reference": _paste("not a spectrum")}),
    ("ms2_fragments", {"structure": {"mode": "SMILES / InChI", "value": "CCO"},
                       "precursor": {"precursor_mz": "", "ion_mode": "Positive", "adduct": "[2M+H]+"}, "ms2": _paste("47.05 100")}),
    ("ms2_fragments", {"structure": {"mode": "SMILES / InChI", "value": "CCO"},
                       "precursor": {"precursor_mz": "", "ion_mode": "Positive", "adduct": "Auto"},
                       "ms2": {"source": "USI", "value": "mzspec:GNPS:x"}}),
])
def test_bad_input_is_400(module_id, inputs):
    r = _run(module_id, inputs)
    assert r.status_code == 400
    assert r.json()["detail"]


def test_unknown_module_is_404():
    assert _run("nope", {}).status_code == 404
