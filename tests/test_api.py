"""API 与页面定义的测试：spec 与 runner 对得上、示例能跑、封装结果与直接调引擎一致、错误输入返回 400。

运行：python -m pytest tests
MetFrag 需要 Java；ESI MS/MS 需要 Docker 和 CFM-ID 镜像；EI MS 需要 NEIMS 权重（python -m mzdesk.ms2.neims_setup）；
Separation 需要 MolGpKa 权重（python -m mzdesk.platform_fit.pka_setup）；Library search 的示例需要 MassBank 索引（python -m mzdesk.ms2.library_setup MassBank）；
Structure prediction 需要 Mono + MS-FINDER（python -m mzdesk.ms2.msfinder_setup，只在本机）。没有时相关测试跳过。
De novo structure 需要 FRIGID（python -m mzdesk.ms2.frigid_setup，只在本机），而且真实运行默认不跑：设 MZDESK_RUN_FRIGID=1 才跑（约 20–40 s、4–5 GB 内存）。
"""
import json
import os
import shutil
import types
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from mzdesk.api.app import app
from mzdesk.api.runners import RUNNERS
from mzdesk.ms2 import cfmid_engine, frigid_engine, library_index, msfinder_engine, neims_engine
from mzdesk.platform_fit import pka_engine
from mzdesk import spec
from mzdesk.spec import CAFFEINE, MODULES, PIPELINES, SUBGROUPS


def _docker_unavailable():
    """Why the CFM-ID image cannot run here (same check as tests/test_cfmid.py), or None."""
    try:
        cfmid_engine.check_docker()
    except cfmid_engine.CFMIDError as exc:
        return str(exc)
    return None


def _msfinder_unavailable():
    """Why MS-FINDER cannot run here (Mono or the install is missing), or None."""
    try:
        msfinder_engine.check_msfinder()
    except msfinder_engine.MsFinderError as exc:
        return str(exc)
    return None


def _frigid_unavailable():
    """Why FRIGID cannot run here (not installed), or None."""
    try:
        frigid_engine.check_frigid()
    except frigid_engine.FrigidError as exc:
        return str(exc)
    return None


client = TestClient(app)
LIVE = [m for m in MODULES if m["live"]]
HAS_JAVA = shutil.which("java") is not None
NO_DOCKER = _docker_unavailable()
NO_NEIMS = not Path(neims_engine.weights_path()).is_file()
NO_MOLGPKA = not Path(pka_engine.weights_path()).is_file()
NO_MASSBANK = not library_index.is_installed("MassBank")
NO_MSFINDER = _msfinder_unavailable()
NO_FRIGID = _frigid_unavailable()
RUN_FRIGID = os.environ.get("MZDESK_RUN_FRIGID") == "1"        # a real FRIGID run takes 20-40 s and 4-5 GB: opt-in
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
    assert ms2 == ["ms2_clean", "ms2_similarity", "ms2_library", "ms1_formula", "ms2_structure", "ms2_denovo", "ms2_ei", "ms2_esi", "ms2_fragments", "exclusivity"]
    assert "ms2" not in {m["id"] for m in MODULES}
    assert next(s["title"] for s in SUBGROUPS if s["id"] == "ms2") == "Spectrum"


def test_predicted_spectrum_page_is_gone():
    assert "ms2_predicted" not in {m["id"] for m in MODULES} | set(RUNNERS)
    assert _run("ms2_predicted", {}).status_code == 404
    assert "ms2_predicted" not in {mid for p in PIPELINES for s in p["stages"] for mid in s["modules"]}


def test_the_workflows():
    assert [p["id"] for p in PIPELINES] == ["esi_ms2_annotation", "ei_ms_annotation", "targeted_method"]
    esi, ei_ms, targeted = PIPELINES
    assert [p["title"] for p in PIPELINES] == ["MS/MS Annotation", "EI MS Annotation", "Targeted Method"]
    assert [(s["id"], s["title"]) for s in esi["stages"]] == [
        (1, "Adduct"), (2, "Formula"), (3, "Library search"), (4, "MS-FINDER structures"), (5, "FRIGID structures"),
        (6, "Merge and re-rank"), (7, "Isomer check"), (8, "Answer")]
    assert "structure" not in {f["id"] for f in esi["inputs"] if f.get("id")}       # the old "Candidate structure" field is gone
    assert len(esi["examples"]) == 2
    assert [s["title"] for s in targeted["stages"]] == ["Identity", "Platform", "Separation", "Ionization", "Adducts", "Non-MS detection",
                                                        "Derivatization", "Known methods", "Theoretical m/z"]
    assert [s["title"] for s in ei_ms["stages"]] == ["Spectrum cleanup", "EI library search", "Derivatization markers", "Candidates", "Derivatives",
                                                  "NEIMS prediction", "Similarity and rerank", "Kovats RI", "Answer"]
    assert len({p["id"] for p in PIPELINES}) == len(PIPELINES)
    assert all([s["id"] for s in p["stages"]] == list(range(1, len(p["stages"]) + 1)) for p in PIPELINES)


def test_molecular_and_detectability_subgroups():
    assert [m["id"] for m in MODULES if m["subgroup"] == "molecular"] == ["identity", "structure_formats", "properties", "lipid_short_name"]
    assert [m["id"] for m in MODULES if m["subgroup"] == "detectability"] == ["platform_fit", "separation", "ionization", "adducts", "non_ms", "derivatization", "known_methods"]
    assert [s["id"] for s in SUBGROUPS][:2] == ["molecular", "detectability"]


def test_identity_offers_comptox_and_chembl_inputs():
    identity = next(m for m in MODULES if m["id"] == "identity")
    kinds = next(f for f in identity["inputs"] if f["id"] == "kind")["options"]
    assert kinds[-2:] == ["CompTox", "ChEMBL"]


def test_pipeline_examples_match_the_form():
    for p in PIPELINES:
        ids = {f["id"] for f in p["inputs"] if f.get("id")}
        for example in p["examples"]:
            assert set(example["values"]) <= ids, (p["id"], example["label"])


def test_block_ids_unique_and_examples_match_inputs():
    for m in MODULES:
        block_ids = [b["id"] for b in m["outputs"] if b.get("id")]
        assert len(block_ids) == len(set(block_ids)), m["id"]
        if m["example"]:
            assert set(m["example"]) <= {f["id"] for f in m["inputs"] if f.get("id")}, m["id"]
            for v in m["example"].values():   # structure examples are names (the first input mode) with a known structure
                if isinstance(v, dict) and "mode" in v:
                    assert v["mode"] == "Name / ID / Formula" and v["value"] in spec.EXAMPLE_STRUCTURES, m["id"]


# ---------------------------------------------------------------- live modules
@pytest.mark.parametrize("mod", [pytest.param(m, id=m["id"], marks=needs_docker if m["id"] == "ms2_esi" else []) for m in LIVE])
def test_example_runs_and_fills_every_block(mod, monkeypatch):
    if mod["id"] == "identity":
        from mzdesk.identity import enrich, resolve
        monkeypatch.setattr(resolve, "resolve", lambda q, kind="Auto", **kw: _fake_resolved(query=q))
        monkeypatch.setattr(enrich, "enrich", lambda r: {"outputs": enrich.empty_outputs(), "warnings": []})
    if mod["id"] == "exclusivity":          # recorded MS-FINDER answer; real runs: MZDESK_RUN_ISOMERS=1 tests/test_isomers.py
        run = json.loads((Path(__file__).parent / "fixtures" / "isomers" / "caffeine_local.json").read_text())["explain_candidates"]
        monkeypatch.setattr(msfinder_engine, "explain_candidates", lambda *a, **kw: run)
    if mod["id"] == "ms2_fragments" and not HAS_JAVA:
        pytest.skip("MetFrag needs Java")
    if mod["id"] == "ms2_esi" and NO_DOCKER:
        pytest.skip(f"CFM-ID needs Docker: {NO_DOCKER}")
    if mod["id"] == "ms2_structure" and NO_MSFINDER:
        pytest.skip(f"MS-FINDER not available: {NO_MSFINDER}")
    if mod["id"] == "ms2_denovo" and not RUN_FRIGID:
        pytest.skip("real FRIGID runs are opt-in: set MZDESK_RUN_FRIGID=1 (about 20-40 s and 4-5 GB of memory)")
    if mod["id"] == "ms2_denovo" and NO_FRIGID:
        pytest.skip(f"FRIGID not available: {NO_FRIGID}")
    if mod["id"] == "ms2_ei" and NO_NEIMS:
        pytest.skip("NEIMS weights not installed (python -m mzdesk.ms2.neims_setup)")
    if mod["id"] == "ms2_library" and NO_MASSBANK:
        pytest.skip("MassBank index not installed (python -m mzdesk.ms2.library_setup MassBank)")
    if mod["id"] == "separation" and NO_MOLGPKA:
        pytest.skip("MolGpKa weights not installed (python -m mzdesk.platform_fit.pka_setup)")
    assert mod["example"], f"{mod['id']} has no example"
    r = _run(mod["id"], mod["example"])
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["status"] == "ok"
    expected = {b["id"] for b in mod["outputs"] if b.get("id")} - OPTIONAL_BLOCKS
    assert expected <= set(body["outputs"]), expected - set(body["outputs"])
    from mzdesk.api import verdicts
    assert body["verdict"] or mod["id"] in verdicts.NO_VERDICT, "no one-line answer above the results"


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
    body = _run("som", {}).json()
    assert body["status"] == "placeholder" and body["outputs"] == {}
    assert not next(m for m in MODULES if m["id"] == "som")["live"]


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


def test_draw_mode_takes_the_editor_smiles():
    from mzdesk.api import common
    assert common.structure({"structure": {"mode": "Draw", "value": "C1C=CC=CC=1", "molfile": "x"}}) == {"kind": "smiles", "smiles": "C1C=CC=CC=1"}
    assert common.structure({"structure": {"mode": "Draw", "value": ""}}, required=False) is None
    r = _run("structure_formats", {"structure": {"mode": "Draw", "value": "C1C=CC=CC=1"}})
    assert r.status_code == 200


# ---------------------------------------------------------------- bad input
@pytest.mark.parametrize("module_id,inputs", [
    ("properties", {"structure": {"mode": "SMILES / InChI", "value": "C1CC"}}),
    ("properties", {"structure": {"mode": "Name / ID / Formula", "value": ""}}),
    ("properties", {"structure": {"mode": "Draw", "value": ""}}),
    ("ms1_formula", {"precursor": {"precursor_mz": "2000", "ion_mode": "Positive", "adduct": "Auto"}}),
    ("ms2_similarity", {"query": _paste("100 1"), "reference": _paste("not a spectrum")}),
    ("ms2_ei", {"structure": {"mode": "SMILES / InChI", "value": "C1CC"}}),
    ("ms2_esi", {"structure": {"mode": "SMILES / InChI", "value": "C1CC"}, "adduct": "[M+H]+"}),
    ("ms2_esi", {"structure": {"mode": "SMILES / InChI", "value": CAFFEINE}, "adduct": "[M+K]+"}),
    ("ms2_esi", {"structure": {"mode": "SMILES / InChI", "value": "C[Si](C)C"}, "adduct": "[M+H]+"}),
    ("ms2_structure", {"ms2": _paste("138.07 100")}),
    ("ms2_structure", {"precursor": {"precursor_mz": "195.0877", "ion_mode": "Positive", "adduct": "[M-H]-"}, "ms2": _paste("138.07 100")}),
    ("ms2_structure", {"precursor": {"precursor_mz": "195.0877", "ion_mode": "Positive", "adduct": "Auto"}, "ms2": _paste("138.07 100"), "halogen": "yes"}),
    ("ms2_structure", {"precursor": {"precursor_mz": "195.0877", "ion_mode": "Positive", "adduct": "Auto"}, "ms2": _paste("138.07 100"), "silicon": 1}),
    ("ms2_denovo", {"ms2": _paste("138.07 100")}),
    ("ms2_denovo", {"precursor": {"precursor_mz": "195.0877", "ion_mode": "Negative", "adduct": "Auto"}, "ms2": _paste("138.07 100"), "formula": "C8H10N4O2"}),
    ("ms2_denovo", {"precursor": {"precursor_mz": "195.0877", "ion_mode": "Positive", "adduct": "Auto"}, "ms2": _paste("138.07 100"), "formula": "C8H10Zn"}),
    ("ms2_denovo", {"precursor": {"precursor_mz": "195.0877", "ion_mode": "Positive", "adduct": "Auto"}, "ms2": _paste("138.07 100"), "formula": "C8H10N4O2", "mode": "Slow"}),
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
    from mzdesk.identity import enrich, resolve
    monkeypatch.setattr(resolve, "resolve", lambda q, kind="Auto", **kw: _fake_resolved(query=q))
    monkeypatch.setattr(enrich, "enrich", lambda r: {"outputs": enrich.empty_outputs(), "warnings": []})
    body = _run("identity", {"structure": {"mode": "Name / ID / Formula", "value": "caffeine"}, "kind": "Auto"}).json()
    assert body["status"] == "ok"
    out = body["outputs"]
    assert set(enrich.OUTPUT_KEYS) <= set(out)
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


def test_identity_runner_forwards_the_comptox_kind(monkeypatch):
    from mzdesk.identity import enrich, resolve
    seen = {}

    def fake_resolve(q, kind="Auto", **kw):
        seen["kind"] = kind
        return _fake_resolved(query=q, kind=kind, ids={}, links={"PubChem": "https://pubchem.ncbi.nlm.nih.gov/compound/2519"})

    monkeypatch.setattr(resolve, "resolve", fake_resolve)
    monkeypatch.setattr(enrich, "enrich", lambda r: {"outputs": enrich.empty_outputs(), "warnings": []})
    body = _run("identity", {"structure": {"mode": "Name / ID / Formula", "value": "DTXSID7020182"}, "kind": "CompTox"}).json()
    assert body["status"] == "ok" and seen["kind"] == "CompTox"
    assert body["outputs"]["resolved"]["values"]["Name"]["note"] == "read as CompTox"


def test_structure_formats_has_four_tabs():
    body = _run("structure_formats", {"structure": {"mode": "SMILES / InChI", "value": CAFFEINE}}).json()
    compounds = body["outputs"]["compare"]["compounds"]
    assert len(compounds) == 1 and compounds[0]["formula"] == "C8H10N4O2"
    assert "width='420px'" in compounds[0]["svg"] and "height='420px'" in compounds[0]["svg"]   # square 2D canvas
    assert compounds[0]["mol3d"]                                                                  # square 3D box in the page
    fm = body["outputs"]["formats"]["sets"][0]["formats"]
    assert list(fm) == ["SMILES", "InChI", "InChIKey", "MOL block"]
    assert fm["InChIKey"] == "RYYVLZVUVIJVGH-UHFFFAOYSA-N" and fm["InChI"].startswith("InChI=1S/C8H10N4O2")


def test_structure_formats_labels_name_the_compounds():
    body = _run("structure_formats", {
        "structure": {"mode": "SMILES / InChI", "value": CAFFEINE, "label": "  Caffeine\n  citrate "},
        "structure_2": {"mode": "Name / ID / Formula", "value": "CCO"},        # no label: the searched name is shown
        "structure_3": {"mode": "SMILES / InChI", "value": spec.TESTOSTERONE, "label": "   "},
    }).json()
    labels = [s["label"] for s in body["outputs"]["formats"]["sets"]]
    assert labels == ["Caffeine citrate", "CCO", "Compound 3"]                                     # label > searched name > Compound n
    assert [c["label"] for c in body["outputs"]["compare"]["compounds"]] == labels


def test_structure_formats_has_a_label_input_on_every_compound():
    fields = next(m for m in spec.MODULES if m["id"] == "structure_formats")["inputs"]
    assert [f["id"] for f in fields if f.get("label_field")] == ["structure", "structure_2", "structure_3"]
    assert all(f["type"] == "structure" for f in fields if f.get("label_field"))


def test_structure_formats_compares_up_to_three_compounds():
    body = _run("structure_formats", {
        "structure": {"mode": "SMILES / InChI", "value": CAFFEINE},
        "structure_2": {"mode": "SMILES / InChI", "value": spec.ASPIRIN},
        "structure_3": {"mode": "SMILES / InChI", "value": spec.TESTOSTERONE},
    }).json()
    compounds = body["outputs"]["compare"]["compounds"]
    assert [c["formula"] for c in compounds] == ["C8H10N4O2", "C9H8O4", "C19H28O2"]
    assert [s["label"] for s in body["outputs"]["formats"]["sets"]] == ["Compound 1", "Compound 2", "Compound 3"]
    assert all(c["svg"] and c["mol3d"] for c in compounds)
    # holes are allowed: only the first and the third compound
    body = _run("structure_formats", {
        "structure": {"mode": "SMILES / InChI", "value": CAFFEINE},
        "structure_3": {"mode": "SMILES / InChI", "value": spec.TESTOSTERONE},
    }).json()
    assert [c["formula"] for c in body["outputs"]["compare"]["compounds"]] == ["C8H10N4O2", "C19H28O2"]
    assert _run("structure_formats", {}).status_code == 400


def test_properties_page_has_no_duplicates():
    body = _run("properties", {"structure": {"mode": "SMILES / InChI", "value": CAFFEINE}}).json()
    out = body["outputs"]
    assert set(out) == {"summary", "descriptors", "functional_groups"}
    assert out["summary"]["values"]["Monoisotopic mass"] == "194.0804"
    rows = out["descriptors"]["rows"]
    assert all(len(r) == 3 for r in rows) and not any("mass" in r[1].lower() or "MW" in r[1] for r in rows)
    assert [r[0] for r in rows] == ["Polarity", "", "", "", "Size & topology", "", "", ""]
    fg = out["functional_groups"]
    assert fg["svg"].startswith("<?xml") and fg["found"]
    for g in fg["found"]:
        assert g["svg"] and g["smarts"] and g["family"] and g["sites"]
        assert len(g["sites"]) == min(g["count"], 6) and all(s["svg"] and s["atoms"] for s in g["sites"])
        assert ("all_svg" in g) == (g["count"] > 1)


def test_platform_page_is_a_portal_and_a_reaction_block():
    mod = next(m for m in MODULES if m["id"] == "platform_fit")
    portal, reaction = mod["outputs"]
    assert (portal["type"], portal["id"], portal["title"]) == ("portal", "platform", "Platform") and portal["help"]
    assert (reaction["type"], reaction["id"], reaction["title"]) == ("reaction", "reactions", "Why and how") and reaction["help"]
    out = _run("platform_fit", mod["example"]).json()["outputs"]
    tiles = {t["label"]: t for t in out["platform"]["tiles"]}
    assert list(tiles) == ["Lipidomics", "Metabolomics", "CIL LC-MS", "GC-MS"]
    # the example (tyrosine) shows everything: CIL chips, GC partial and the reaction block
    assert mod["example"]["structure"] == {"mode": "Name / ID / Formula", "value": "L-tyrosine"}
    assert tiles["CIL LC-MS"]["status"] == "yes" and tiles["CIL LC-MS"]["chips"] and tiles["GC-MS"]["status"] == "partial"
    assert [o["label"] for o in out["reactions"]["options"]] == ["CIL · Dansyl chloride", "CIL · DmPA bromide", "GC · TMS"]


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


def test_verdict_never_breaks_a_run_and_reads_the_outputs():
    from mzdesk.api import verdicts
    assert verdicts.build("separation", {}, {"modes": "not a dict"}) is None            # a malformed block gives no verdict, no error
    assert verdicts.build("unknown_module", {}, {}) is None
    out = {"modes": {"chips": [{"label": "RP", "level": "good"}, {"label": "HILIC", "level": "poor"}]},
           "ph_summary": {"values": {"Fraction neutral": {"value": "98.0 %"}, "logD at pH": {"value": "3.1"}}}}
    assert verdicts.build("separation", {"ph": 2.7}, out) == "<b>RP</b> fits best at pH 2.7: mostly neutral (98.0 %), logD 3.1."


def test_spec_carries_the_version():
    from mzdesk import __version__
    assert client.get("/api/spec").json()["version"] == __version__


def test_stats_record_runs_and_views_and_need_the_password(monkeypatch):
    from mzdesk.api import stats
    monkeypatch.setattr(stats, "check_password", lambda pw: pw == "right")
    client.post("/api/track", json={"path": "/modules/separation", "ref": "https://www.google.com/search"})
    client.post("/api/track", json={"path": "/stats"})                                    # the stats page itself is not counted
    client.post("/api/modules/properties/run", json={"structure": {"mode": "SMILES / InChI", "value": "CCO"}})
    client.post("/api/modules/properties/run", json={"structure": {"mode": "SMILES / InChI", "value": "not a smiles"}})
    assert client.post("/api/stats", json={"password": "wrong"}).status_code == 401
    s = client.post("/api/stats", json={"password": "right", "days": 7}).json()
    assert len(s["daily"]) == 7 and s["daily"][-1]["runs"] == 2 and s["daily"][-1]["views"] == 1
    assert s["modules"] == [{"module": "properties", "runs": 2, "errors": 1, "avg_ms": s["modules"][0]["avg_ms"]}]
    assert s["pages"] == [{"path": "/modules/separation", "views": 1}] and s["referrers"] == [{"host": "www.google.com", "views": 1}]
    assert s["totals"]["today"]["visitors"] == 1          # same IP + browser on one day = one visitor
    assert s["errors"][0]["status"] == 400


def test_stats_password_matches_the_page_lock():
    import re
    from pathlib import Path
    from mzdesk.api import stats
    lock_ts = (Path(__file__).parents[1] / "web" / "src" / "lock.ts").read_text()
    assert stats.SALT_HEX in lock_ts and stats.HASH_HEX in lock_ts
    assert re.search(r"ITERATIONS = 200_000", lock_ts) and stats.ITERATIONS == 200_000
    assert not stats.check_password("") and not stats.check_password("wrong")


# ---------------------------------------------------------------- Spectrum input-regression cases
def test_ms1_formula_formula_adduct_and_charge_validation():
    from mzdesk.spec import CAFFEINE
    base = {"precursor": {"precursor_mz": "153.0", "ion_mode": "Negative", "adduct": "[M-H]-"}, "formula": "CCl4"}
    r = _run("ms1_formula", base)
    assert r.status_code == 400 and "cannot form [M-H]-" in r.json()["detail"]
    for formula in ("C8H10N4O2+", "C8H10N4O2•"):
        r = _run("ms1_formula", {"precursor": {"precursor_mz": "195.0877", "ion_mode": "Positive", "adduct": "Auto"},
                                 "formula": formula})
        assert r.status_code == 400 and "neutral formula" in r.json()["detail"]
    r = _run("ms1_formula", {"precursor": {"precursor_mz": "195.0877", "ion_mode": "Positive", "adduct": "[M-H]-"},
                             "formula": CAFFEINE})
    assert r.status_code == 400 and "not a positive-mode adduct" in r.json()["detail"]


def test_ms1_formula_reads_precursor_from_the_ms2_header_and_warns_about_unused_ms1():
    mgf = "BEGIN IONS\nPEPMASS=195.0877\n138.0662 100\nEND IONS"
    r = _run("ms1_formula", {"precursor": {"precursor_mz": "", "ion_mode": "Positive", "adduct": "Auto"},
                             "ms2": _paste(mgf), "formula": "C8H10N4O2", "ms1": "195.0877 100\n197.09 90"})
    assert r.status_code == 200
    body = r.json()
    assert any("read from the MS2 spectrum header" in w for w in body["warnings"])
    assert any("does not use the MS1 isotope pattern" in w for w in body["warnings"])


def test_ms2_clean_ion_mode_is_case_insensitive_and_writes_the_charge():
    body = _run("ms2_clean", {"spectrum": _paste("100 1\n200 10"), "ion_mode": "negative"}).json()
    assert "CHARGE=1-" in body["outputs"]["formats"]["formats"]["MGF"]
    assert _run("ms2_clean", {"spectrum": _paste("100 1"), "ion_mode": "sideways"}).status_code == 400


def test_malformed_composite_and_spectrum_values_are_400_not_500():
    for module, inputs in (
        ("ms2_clean", {"spectrum": "100 1\n200 2"}),
        ("ms2_similarity", {"query": _paste("100 1"), "reference": _paste("100 1"), "query_precursor": "x"}),
        ("ms2_structure", {"precursor": "195.1", "ms2": _paste("100 1")}),
    ):
        r = _run(module, inputs)
        assert r.status_code == 400, (module, r.text)


def test_subnormal_intensities_do_not_500():
    r = _run("ms2_clean", {"spectrum": _paste("100 1e-308\n200 2e-308"), "normalize": "Base peak = 100"})
    assert r.status_code == 200
    assert all(row[1] == "0" or row[1] == "100" or row[1] for row in r.json()["outputs"]["peaks"]["rows"])


def test_similarity_tolerance_is_capped():
    r = _run("ms2_similarity", {"query": _paste("100 100"), "reference": _paste("107 100"), "tol_da": "10"})
    assert r.status_code == 400 and "at most 5" in r.json()["detail"]


def test_similarity_reads_headers_and_says_so():
    q = "BEGIN IONS\nPEPMASS=250\n100 1\n150 3\nEND IONS"
    r = "BEGIN IONS\nPEPMASS=251\n100 2\n151 3\nEND IONS"
    body = _run("ms2_similarity", {"query": _paste(q), "reference": _paste(r)}).json()
    assert any("Query precursor m/z 250.0000 is read" in w for w in body["warnings"])

