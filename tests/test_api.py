"""API 与页面定义的测试：spec 与 runner 对得上、示例能跑、封装结果与直接调引擎一致、错误输入返回 400。

运行：python -m pytest tests
MetFrag 需要 Java；没有 Java 时相关测试跳过。
"""
import shutil
import types

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
    assert ms2 == ["ms2_clean", "ms2_predicted", "ms2_library", "ms2_similarity", "ms2_fragments", "ms2_denovo"]
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


def test_clean_page_reads_any_format_and_removes_peaks():
    body = _run("ms2_clean", {"spectrum": _paste("[[100.1, 5], [200.2, 100], [250.0, 40], [300.3, 0.5]]"), "precursor_mz": "250",
                              "noise_pct": "1", "remove_precursor": True, "normalize": "Base peak = 999", "sort_by": "Intensity"}).json()
    out = body["outputs"]
    assert out["summary"]["values"]["Detected format"] == "GNPS peaks_json"
    assert out["summary"]["values"]["Removed"]["note"] == "1 noise · 1 precursor"
    assert [r[3] for r in out["peaks"]["rows"]] == ["yes", "yes", "no", "no"]
    assert out["formats"]["formats"]["Peak list"].splitlines() == ["200.2\t999", "100.1\t49.95"]


def test_clean_page_reads_precursor_from_pasted_record():
    body = _run("ms2_clean", {"spectrum": _paste("BEGIN IONS\nPEPMASS=250.5\n100 1\n250.5 9\nEND IONS"), "remove_precursor": True}).json()
    assert body["outputs"]["summary"]["values"]["Peaks out"] == "1"
    assert any("250.5000" in w for w in body["warnings"])


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
    ("properties", {"structure": {"mode": "Name / ID / Formula", "value": ""}}),
    ("ms1_formula", {"precursor": {"precursor_mz": "2000", "ion_mode": "Positive", "adduct": "Auto"}}),
    ("ms2_similarity", {"query": _paste("100 1"), "reference": _paste("not a spectrum")}),
    ("ms2_clean", {"spectrum": _paste("not a spectrum")}),
    ("ms2_clean", {"spectrum": _paste("100 1"), "noise_pct": "150"}),
    ("ms2_clean", {"spectrum": _paste("100 1"), "normalize": "Everything = 7"}),
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


# ---------------------------------------------------------------- identity / name resolver
def _fake_resolved(**kw):
    base = dict(query="caffeine", kind="Name", smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                inchi="InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
                inchikey="RYYVLZVUVIJVGH-UHFFFAOYSA-N", cid=2519, name="caffeine", iupac="1,3,7-trimethylpurine-2,6-dione",
                ids={"CAS": ["58-08-2"], "HMDB": ["HMDB0001847"], "KEGG": ["C07481"]},
                synonyms=["caffeine", "1,3,7-trimethylxanthine"], candidates=[(2519, "Caffeine"), (1, "Other")],
                stereo_fallback=False, warnings=[],
                links={"PubChem": "https://pubchem.ncbi.nlm.nih.gov/compound/2519",
                       "HMDB": "https://hmdb.ca/metabolites/HMDB0001847", "CAS": "https://commonchemistry.cas.org/detail?cas_rn=58-08-2"})
    base.update(kw)
    return types.SimpleNamespace(**base)


def test_identity_runner_offline(monkeypatch):
    from mzdesk.identity import resolve
    monkeypatch.setattr(resolve, "resolve", lambda q, kind="Auto", **kw: _fake_resolved(query=q))
    body = _run("identity", {"structure": {"mode": "Name / ID / Formula", "value": "caffeine"}, "kind": "Auto"}).json()
    assert body["status"] == "ok"
    out = body["outputs"]
    rows = out["identifiers"]["rows"]
    assert len(rows) == 14
    hmdb = next(r for r in rows if r[0] == "HMDB")
    assert hmdb[1] == "HMDB0001847" and hmdb[2] == {"text": "HMDB0001847", "href": "https://hmdb.ca/metabolites/HMDB0001847"}
    pc = next(r for r in rows if r[0] == "PubChem CID")
    assert pc[2] == {"text": "2519", "href": "https://pubchem.ncbi.nlm.nih.gov/compound/2519"}
    assert next(r for r in rows if r[0] == "ChEMBL")[1:] == ["—", "—"]
    assert out["resolved"]["values"] == {"Input read as": "Name", "Name": "caffeine", "PubChem CID": "2519", "Formula": "C8H10N4O2"}
    assert out["candidates"]["rows"] == [["1", "Other"]]
    assert out["depiction"]["svg"]


def test_identity_runner_empty_and_resolver_error_are_400(monkeypatch):
    from mzdesk.identity import resolve
    assert _run("identity", {"structure": {"mode": "Name / ID / Formula", "value": " "}}).status_code == 400

    def boom(*a, **k):
        raise resolve.ResolveError("No compound found for 'zzz'.")
    monkeypatch.setattr(resolve, "resolve", boom)
    r = _run("identity", {"structure": {"mode": "Name / ID / Formula", "value": "zzz"}})
    assert r.status_code == 400 and "zzz" in r.json()["detail"]


def test_structure_name_mode_uses_resolver(monkeypatch):
    from mzdesk.api import common
    from mzdesk.identity import resolve
    monkeypatch.setattr(resolve, "to_smiles", lambda q, kind="Auto": "CCO")
    got = common.structure({"structure": {"mode": "Name / ID / Formula", "value": "ethanol"}})
    assert got == {"kind": "smiles", "smiles": "CCO"}


def test_structure_kegg_id_not_treated_as_formula(monkeypatch):
    from mzdesk.api import common
    from mzdesk.identity import resolve
    calls = []
    monkeypatch.setattr(resolve, "to_smiles", lambda q, kind="Auto": calls.append(q) or "OC(=O)CC(O)=O")
    f = lambda v: common.structure({"structure": {"mode": "Name / ID / Formula", "value": v}}, allow_formula=True)
    assert f("C00001")["kind"] == "smiles" and calls == ["C00001"]
    assert f("C8H10N4O2") == {"kind": "formula", "formula": "C8H10N4O2"} and calls == ["C00001"]


def test_resolver_input_error_is_400(monkeypatch):
    from mzdesk.identity import resolve

    def boom(q, kind="Auto"):
        raise resolve.ResolveError("not found")
    monkeypatch.setattr(resolve, "to_smiles", boom)
    r = _run("properties", {"structure": {"mode": "Name / ID / Formula", "value": "zzz"}})
    assert r.status_code == 400 and r.json()["detail"] == "Structure: not found"
