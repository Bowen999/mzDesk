"""API 与页面定义的测试：spec 与 runner 对得上、示例能跑、封装结果与直接调引擎一致、错误输入返回 400。

运行：python -m pytest tests
MetFrag 需要 Java；ESI MS/MS 需要 Docker 和 CFM-ID 镜像；EI MS 需要 NEIMS 权重（python -m mzdesk.ms2.neims_setup）；
Separation 需要 MolGpKa 权重（python -m mzdesk.platform_fit.pka_setup）。没有时相关测试跳过。
"""
import shutil
import types
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from mzdesk.api.app import app
from mzdesk.api.runners import RUNNERS
from mzdesk.ms2 import cfmid_engine, neims_engine
from mzdesk.platform_fit import pka_engine
from mzdesk.spec import CAFFEINE, MODULES, PIPELINES, SUBGROUPS


def _docker_unavailable():
    """Why the CFM-ID image cannot run here (same check as tests/test_cfmid.py), or None."""
    try:
        cfmid_engine.check_docker()
    except cfmid_engine.CFMIDError as exc:
        return str(exc)
    return None


client = TestClient(app)
LIVE = [m for m in MODULES if m["live"]]
HAS_JAVA = shutil.which("java") is not None
NO_DOCKER = _docker_unavailable()
NO_NEIMS = not Path(neims_engine.weights_path()).is_file()
NO_MOLGPKA = not Path(pka_engine.weights_path()).is_file()
needs_docker = [pytest.mark.docker, pytest.mark.skipif(NO_DOCKER is not None, reason=f"Docker unavailable: {NO_DOCKER}")]
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
    assert ms2 == ["ms2_clean", "ms2_ei", "ms2_esi", "ms2_library", "ms2_similarity", "ms2_fragments", "ms2_denovo"]
    assert "ms2" not in {m["id"] for m in MODULES}
    assert next(s["title"] for s in SUBGROUPS if s["id"] == "ms2") == "Spectrum"


def test_predicted_spectrum_page_is_gone():
    assert "ms2_predicted" not in {m["id"] for m in MODULES} | set(RUNNERS)
    assert _run("ms2_predicted", {}).status_code == 404
    assert "ms2_predicted" not in {mid for p in PIPELINES for s in p["stages"] for mid in s["modules"]}


def test_pipelines_use_the_spectrum_prediction_pages():
    stages = {(p["id"], s["title"]): s["modules"] for p in PIPELINES for s in p["stages"]}
    assert stages[("verify", "MS2")][0] == "ms2_esi"
    assert stages[("feature_prediction", "Spectrum prediction (ESI MS/MS, EI MS)")] == ["ms2_esi", "ms2_ei"]


def test_molecular_and_detectability_subgroups():
    assert [m["id"] for m in MODULES if m["subgroup"] == "molecular"] == ["identity", "structure_formats", "properties"]
    assert [m["id"] for m in MODULES if m["subgroup"] == "detectability"] == ["platform_fit", "separation", "ionization", "adducts", "non_ms", "derivatization", "known_methods"]
    assert [s["id"] for s in SUBGROUPS][:2] == ["molecular", "detectability"]


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
@pytest.mark.parametrize("mod", [pytest.param(m, id=m["id"], marks=needs_docker if m["id"] == "ms2_esi" else []) for m in LIVE])
def test_example_runs_and_fills_every_block(mod):
    if mod["id"] == "ms2_fragments" and not HAS_JAVA:
        pytest.skip("MetFrag needs Java")
    if mod["id"] == "ms2_esi" and NO_DOCKER:
        pytest.skip(f"CFM-ID needs Docker: {NO_DOCKER}")
    if mod["id"] == "ms2_ei" and NO_NEIMS:
        pytest.skip("NEIMS weights not installed (python -m mzdesk.ms2.neims_setup)")
    if mod["id"] == "separation" and NO_MOLGPKA:
        pytest.skip("MolGpKa weights not installed (python -m mzdesk.platform_fit.pka_setup)")
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
    body = _run("ms2_library", {}).json()
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


def test_ion_mz_and_formula_is_the_shared_ion_calculation():
    from mzdesk.ms1 import msbuddy_engine as bd
    mz, ion = bd.ion_mz_and_formula("C8H10N4O2", "[M+H]+")
    assert mz == pytest.approx(195.0877, abs=1e-4) and ion == "C8H11N4O2+"
    assert bd.ion_mz_and_formula("C8H10N4O2", "[M+H]+")[0] == pytest.approx(bd.ion_mz("C8H10N4O2", "[M+H]+"))
    mz, ion = bd.ion_mz_and_formula("C5H14NO", "[M]+")                    # choline: permanent cation, no H added
    assert mz == pytest.approx(104.1070, abs=1e-4) and ion == "C5H14NO+"
    assert bd.ion_mz_and_formula("CH4", "[M-2H2O+H]+") == (None, None)      # the adduct removes atoms the formula lacks
    assert bd.ion_mz_and_formula("not a formula", "[M+H]+") == (None, None)


def test_record_name_is_the_typed_name_else_the_formula():
    from mzdesk.api import common
    assert common.record_name({"structure": {"mode": "SMILES / InChI", "value": CAFFEINE}}, CAFFEINE) == "C8H10N4O2"
    assert common.record_name({"structure": {"mode": "Name / ID / Formula", "value": " caffeine "}}, CAFFEINE) == "caffeine"
    assert common.record_name({"structure": {"mode": "Name / ID / Formula", "value": ""}}, CAFFEINE) == "C8H10N4O2"
    assert common.record_name({}, "CCO") == "C2H6O"


# ---------------------------------------------------------------- bad input
@pytest.mark.parametrize("module_id,inputs", [
    ("properties", {"structure": {"mode": "SMILES / InChI", "value": "C1CC"}}),
    ("properties", {"structure": {"mode": "Name / ID / Formula", "value": ""}}),
    ("ms1_formula", {"precursor": {"precursor_mz": "2000", "ion_mode": "Positive", "adduct": "Auto"}}),
    ("ms2_similarity", {"query": _paste("100 1"), "reference": _paste("not a spectrum")}),
    ("ms2_ei", {"structure": {"mode": "SMILES / InChI", "value": "C1CC"}}),
    ("ms2_esi", {"structure": {"mode": "SMILES / InChI", "value": "C1CC"}, "adduct": "[M+H]+"}),
    ("ms2_esi", {"structure": {"mode": "SMILES / InChI", "value": CAFFEINE}, "adduct": "[M+K]+"}),
    ("ms2_esi", {"structure": {"mode": "SMILES / InChI", "value": "C[Si](C)C"}, "adduct": "[M+H]+"}),
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


@pytest.mark.docker
@pytest.mark.skipif(NO_DOCKER is not None, reason=f"Docker unavailable: {NO_DOCKER}")
def test_esi_msrb_only_adduct_on_caffeine_is_400_with_a_clear_message():
    r = _run("ms2_esi", {"structure": {"mode": "SMILES / InChI", "value": CAFFEINE}, "adduct": "[M+Na]+"})
    assert r.status_code == 400
    detail = r.json()["detail"]
    assert "[M+Na]+" in detail and "MSRB rules do not cover this structure" in detail and "[M+H]+ and [M-H]-" in detail


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
    names, ids = out["names"]["rows"], out["database_ids"]["rows"]
    assert [r[0] for r in names] == ["Name", "IUPAC", "CAS"]
    assert names[2][1] == {"text": "58-08-2", "href": "https://commonchemistry.cas.org/detail?cas_rn=58-08-2"}
    assert names[0][1] == "caffeine"
    assert len(ids) == 9 and all(len(r) == 2 for r in ids)          # no separate Link column
    assert ids[2] == ["HMDB", {"text": "HMDB0001847", "href": "https://hmdb.ca/metabolites/HMDB0001847"}]
    assert ids[1] == ["PubChem CID", {"text": "2519", "href": "https://pubchem.ncbi.nlm.nih.gov/compound/2519"}]
    assert next(r for r in ids if r[0] == "ChEMBL")[1] == "—"
    # SMILES / InChI live on the Structure page; InChIKey is the first Database IDs row
    assert not {"SMILES", "InChI"} & {r[0] for r in names + ids}
    assert ids[0] == ["InChIKey", "RYYVLZVUVIJVGH-UHFFFAOYSA-N"]
    assert out["resolved"]["values"] == {"Name": {"value": "caffeine", "note": "read as Name"}, "Formula": "C8H10N4O2",
                                         "PubChem CID": "2519"}
    assert out["candidates"]["rows"] == [["1", "Other"]]
    assert out["depiction"]["svg"]


def test_structure_formats_has_four_tabs():
    body = _run("structure_formats", {"structure": {"mode": "SMILES / InChI", "value": CAFFEINE}}).json()
    fm = body["outputs"]["formats"]["formats"]
    assert list(fm) == ["SMILES", "InChI", "InChIKey", "MOL block"]
    assert fm["InChIKey"] == "RYYVLZVUVIJVGH-UHFFFAOYSA-N" and fm["InChI"].startswith("InChI=1S/C8H10N4O2")
    assert body["outputs"]["depiction"]["svg"] and "viewer3d" in body["outputs"]


def test_properties_page_has_no_duplicates():
    body = _run("properties", {"structure": {"mode": "SMILES / InChI", "value": CAFFEINE}}).json()
    out = body["outputs"]
    assert set(out) == {"summary", "descriptors"}
    assert out["summary"]["values"]["Monoisotopic mass"] == "194.0804"
    rows = out["descriptors"]["rows"]
    assert all(len(r) == 3 for r in rows) and not any("mass" in r[1].lower() or "MW" in r[1] for r in rows)
    assert [r[0] for r in rows] == ["Polarity", "", "", "", "Size & topology", "", "", ""]


def test_platform_fit_section_is_renamed():
    mod = next(m for m in MODULES if m["id"] == "platform_fit")
    detect = next(b for b in mod["outputs"] if b.get("id") == "detect")
    assert detect["title"] == "Lipidomics vs metabolomics" and detect["items"][0]["label"] == "Best platform"
    assert detect["items"][0]["unit"] == "Lipidome / Metabolome / Both / Neither"


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
