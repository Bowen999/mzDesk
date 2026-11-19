"""Lipid Fragment Explanation workflow (mzdesk/workflows/lipid_fragments.py): the spec entry, the form, the Auto adduct, the three steps
(names, illustrative structure, MetFrag) and a run through the job API. The MetFrag runs need Java and are skipped without it."""
import shutil
import time

import pytest
from fastapi.testclient import TestClient

from mzdesk import spec
from mzdesk.api import app as app_module
from mzdesk.api.common import InputError
from mzdesk.api.runners.lipid import run_lipid_short_name
from mzdesk.workflows import lipid_fragments as lf

client = TestClient(app_module.app)
HAS_JAVA = shutil.which("java") is not None
needs_java = pytest.mark.skipif(not HAS_JAVA, reason="MetFrag needs Java")
PIPE = next(p for p in spec.PIPELINES if p["id"] == "lipid_fragments")
DOPG = r"CCCCCCCC/C=C\CCCCCCCC(=O)OC[C@H](COP(=O)(O)OCC(O)CO)OC(=O)CCCCCCC/C=C\CCCCCCCC"     # PG 18:1(9Z)/18:1(9Z)
MS2 = {"source": "Paste", "value": spec.PG_36_2_MS2}


@pytest.fixture(autouse=True)
def clean_jobs():
    app_module._jobs.clear()
    yield
    app_module._jobs.clear()


def form(lipid="PG 18:1/18:1", mode="Lipid name", mz="773.5359", ion="Negative", adduct="Auto", **over):
    return {"lipid": {"mode": mode, "value": lipid}, "precursor": {"precursor_mz": mz, "ion_mode": ion, "adduct": adduct}, "ms2": MS2, **over}


def wait_done(job_id, timeout=120):
    t0 = time.time()
    while time.time() - t0 < timeout:
        body = client.get(f"/api/pipelines/jobs/{job_id}").json()
        if body["status"] not in ("running", "queued"):
            return body
        time.sleep(0.2)
    raise AssertionError("the job did not finish")


# ---------------------------------------------------------------- spec
def test_the_spec_entry():
    assert PIPE["title"] == "Lipid Fragment Explanation" and PIPE["one_liner"].startswith("Lipid name or structure + MS2 spectrum →")
    assert [(s["id"], s["title"]) for s in PIPE["stages"]] == [(1, "Lipid name"), (2, "Structure"), (3, "Fragment explanation")]
    assert lf.STEP_TITLES == {1: "Lipid name", 2: "Structure", 3: "Fragment explanation"}
    fields = {f["id"]: f for f in PIPE["inputs"]}
    lipid = fields["lipid"]
    assert lipid["type"] == "structure" and lipid["modes"] == ["Lipid name", "Draw", "SMILES / InChI"] and lipid["required"]
    assert lipid["network_modes"] == [], "a lipid name is read locally: nothing goes online"
    assert fields["tree_depth"]["value"] == 3 and (fields["tree_depth"]["min"], fields["tree_depth"]["max"]) == (1, 4)
    assert fields["precursor"]["adducts"] == spec.METFRAG_ADDUCTS and fields["precursor"]["header_fallback"]
    assert fields["ms2"]["required"] and not fields["precursor"].get("required")
    assert "network" not in PIPE and PIPE["based_on"] == ["LIPID+", "MetFrag", "RDKit"]
    assert [e["values"]["lipid"]["value"] for e in PIPE["examples"]] == ["PG 18:1/18:1", "PG 36:2"]


# ---------------------------------------------------------------- form
def test_a_name_is_read_and_the_numbers_default():
    q = lf.parse_inputs(form())
    assert (q.name, q.smiles, q.mode, q.lipid["L2"]) == ("PG 18:1/18:1", "", "Lipid name", "PG 18:1_18:1")
    assert (q.precursor_mz, q.positive, q.adduct, q.depth, q.tol_abs, q.tol_ppm) == (773.5359, False, None, 3, 0.002, 10.0)
    assert len(q.peaks) == 13 and q.warnings == []


def test_a_smiles_or_a_drawing_gives_a_structure_and_no_name():
    for mode in ("SMILES / InChI", "Draw"):
        q = lf.parse_inputs(form(DOPG, mode))
        assert (q.smiles, q.lipid, q.name, q.mode) == (DOPG, None, "", mode)


def test_the_header_precursor_is_used_when_the_field_is_blank():
    q = lf.parse_inputs(form(mz="", ms2={"source": "Paste", "value": "BEGIN IONS\nPEPMASS=773.5359\n281.2471 100\n152.9962 20\nEND IONS"}))
    assert q.precursor_mz == 773.5359 and q.warnings == ["Precursor m/z 773.5359 is read from the spectrum."]
    assert "PEPMASS" not in q.ms2["value"] and q.ms2["source"] == "Paste"


@pytest.mark.parametrize("bad, message", [
    (form(""), "Lipid is required"),
    (form("POPC"), "Common names are not supported"),
    (form("PC 16:0/18:1", adduct="[M+H]+"), "does not match the ion mode"),
    (form("PC 16:0/18:1", ion="Positive", adduct="[M+2H]2+"), "MetFrag does not support"),
    (form("C1CC", "SMILES / InChI"), "invalid SMILES"),
    (form(lipid=None, mode="Name / ID / Formula"), "unknown input mode"),
    (form(tree_depth=7), "at most 4"),
    ({**form(), "ms2": {"source": "Paste", "value": ""}}, "no valid peaks")])
def test_a_bad_form_is_an_input_error(bad, message):
    with pytest.raises(InputError, match=message):
        lf.parse_inputs(bad)


def test_a_bad_form_is_a_400_and_starts_no_job():
    r = client.post("/api/pipelines/lipid_fragments/run", json=form("POPC"))
    assert r.status_code == 400 and "Use a shorthand" in r.json()["detail"] and not app_module._jobs


# ---------------------------------------------------------------- Auto adduct
def test_auto_picks_the_adduct_that_fits_the_precursor():
    pc = "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\\CCCCCCCC"           # PC 16:0/18:1
    assert lf.pick_adduct(pc, True, 760.5851) == ("[M+H]+", None)
    assert lf.pick_adduct(pc, True, 782.5670) == ("[M+Na]+", None)
    assert lf.pick_adduct(pc, False, 804.5760) == ("[M+HCOO]-", None)
    assert lf.pick_adduct(pc, False, 818.5917) == ("[M+CH3COO]-", None)
    assert lf.pick_adduct(pc, True, None) == ("[M+H]+", None) and lf.pick_adduct(pc, False, None) == ("[M-H]-", None)
    adduct, note = lf.pick_adduct(pc, True, 700.1)
    assert adduct == "[M+H]+" and "fits none of the positive adducts" in note and "760.5851" in note


# ---------------------------------------------------------------- steps
def test_names_come_from_the_lipid_short_name_page():
    q = lf.parse_inputs(form())
    status, note, blocks, got = lf.step_names(q)
    page = run_lipid_short_name({"name": "PG 18:1/18:1"})["outputs"]
    assert status == "done" and note == "LIPID+ PG 36:2 · PG 18:1_18:1"
    assert got["cards"]["lipid_plus"] == page["lipid_plus"] and got["cards"]["lipid_maps"] == page["lipid_maps"]
    assert [b["id"] for b in blocks] == ["lipid_plus", "lipid_maps"] and blocks[0]["title"] == "LIPID+"


def test_the_structure_of_a_name_is_the_illustration_and_says_what_is_assumed():
    status, note, blocks, got = lf.step_structure(lf.parse_inputs(form("PG 36:2")))
    assert status == "done" and got["illustrative"] and "<svg" in got["svg"]
    assert got["assumptions"][:2] == ["Chains split by the default rule", "C=C positions assumed"]
    assert note.startswith("C42H79O10P, illustrative: chains split by the default rule; C=C positions assumed")
    assert blocks[0]["label"] == "Illustrative structure" and blocks[0]["id"] == "structure"


def test_a_structure_is_drawn_as_given():
    status, note, blocks, got = lf.step_structure(lf.parse_inputs(form(DOPG, "SMILES / InChI")))
    assert (status, note, got["smiles"], got["illustrative"], got["assumptions"]) == ("done", "Structure as given.", DOPG, False, [])
    assert blocks[0]["label"] == "Structure" and "help" not in blocks[0]


def test_a_class_without_a_template_keeps_its_names_and_skips_the_rest():
    q = lf.parse_inputs(form("ST 27:1;O", ion="Positive", mz=""))
    r = lf.run_query(q)
    assert [s["status"] for s in r["steps"]] == ["done", "skipped", "skipped"]
    assert "No illustrative structure for this lipid class" in r["steps"][1]["note"] and "SMILES" in r["steps"][1]["note"]
    assert r["answer"]["statement"].startswith("ST 27:1;O has no illustrative structure") and "draw it" in r["answer"]["statement"]
    assert [b["id"] for b in r["blocks"]][:2] == ["lipid_plus", "lipid_maps"] and "structure" not in {b["id"] for b in r["blocks"]}
    assert r["warnings"] == []


# ---------------------------------------------------------------- whole runs (MetFrag)
@needs_java
def test_a_name_runs_end_to_end():
    r = lf.run_query(lf.parse_inputs(form()))
    assert [s["status"] for s in r["steps"]] == ["done", "done", "done"]
    assert r["answer"]["statement"] == "The illustrative structure of PG 18:1/18:1 explains 4 of 13 peaks (91.6 % of the intensity)."
    assert r["answer"]["level"] == "Assumed parts"
    assert dict(r["answer"]["facts"]) == {"Lipid": "PG 18:1/18:1", "Formula": "C42H79O10P", "Adduct": "[M-H]- (Auto)", "Tree depth": "3"}
    assert [b["id"] for b in r["blocks"]] == ["lipid_plus", "lipid_maps", "structure", "summary", "spectrum", "peaks", "assumptions"]
    blocks = {b["id"]: b for b in r["blocks"]}
    assert blocks["summary"]["data"]["values"]["Explained peaks"] == "4 / 13" and blocks["summary"]["title"] == "Fragment explanation"
    assert len(blocks["peaks"]["data"]["rows"]) == 13 and blocks["spectrum"]["data"]["plotly"]
    assert blocks["lipid_plus"]["data"]["values"] == {"L1": "PG 36:2", "L2": "PG 18:1_18:1", "L3": "PG 18:1_18:1"}
    assert r["warnings"] == []


@needs_java
def test_a_smiles_runs_with_dash_names_and_the_same_explanation():
    r = lf.run_query(lf.parse_inputs(form(DOPG, "SMILES / InChI")))
    assert [s["status"] for s in r["steps"]] == ["skipped", "done", "done"] and r["steps"][0]["note"] == lf.NO_NAME_HELP
    assert r["answer"]["statement"] == "The structure explains 4 of 13 peaks (91.6 % of the intensity)." and r["answer"]["level"] is None
    plus = {b["id"]: b for b in r["blocks"]}["lipid_plus"]["data"]["values"]
    assert set(plus) == {"L1", "L2", "L3"} and all(v == {"value": "—", "help": lf.NO_NAME_HELP} for v in plus.values())


@needs_java
def test_a_far_precursor_falls_back_with_one_warning():
    r = lf.run_query(lf.parse_inputs(form(mz="700.1")))
    assert len(r["warnings"]) == 1 and "fits none of the negative adducts" in r["warnings"][0]


@needs_java
def test_the_job_api_runs_the_workflow():
    r = client.post("/api/pipelines/lipid_fragments/run", json=form(tree_depth=2))
    assert r.status_code == 200
    body = wait_done(r.json()["job_id"])
    assert body["status"] == "done" and body["error"] is None
    result = body["result"]
    assert result["pipeline"] == "lipid_fragments" and [s["status"] for s in result["steps"]] == ["done", "done", "done"]
    assert ("Tree depth", "2") in [tuple(f) for f in result["answer"]["facts"]]
    assert body["steps"][2]["blocks"][0]["id"] == "summary"                              # a step's blocks reach the job while it runs
