"""Scaffold of the EI MS Annotation and Targeted Method workflows (mzdesk/workflows/): the spec entries, the form parsing, the step runner
and the generic result, and EI MS Annotation run end to end through the job API on a bare spectrum. The real step logic has its own tests.
"""
import logging
import subprocess
import sys
import time
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from mzdesk import spec
from mzdesk.api import app as app_module
from mzdesk.api import blocks as blk
from mzdesk.api.common import InputError
from mzdesk.platform_fit import adduct_rules
from mzdesk.workflows import common as wc
from mzdesk.workflows import ei_ms_annotation as ems
from mzdesk.workflows import targeted_method as tm

client = TestClient(app_module.app)
REPO = Path(__file__).resolve().parents[1]
TARGETED = next(p for p in spec.PIPELINES if p["id"] == "targeted_method")
EI = next(p for p in spec.PIPELINES if p["id"] == "ei_ms_annotation")
NEW_WORKFLOWS = [("ei_ms_annotation", ems), ("targeted_method", tm)]


@pytest.fixture(autouse=True)
def clean_jobs():
    app_module._jobs.clear()
    yield
    app_module._jobs.clear()


def targeted_form(**over):
    base = {"structure": {"mode": "SMILES / InChI", "value": "CCO"}, "ph": "", "additive": "Formic acid", "instrument": "Any"}
    return {**base, **over}


def ei_form(**over):
    base = {"spectrum": {"source": "Paste", "value": "73 999\n147 250\n204 120\n217 90"}, "derivatization": "Auto (from marker ions)",
            "molecular_ion": "", "candidates": "", "alkanes": "", "sample_rt": ""}
    return {**base, **over}


def fields(pipeline):
    return {f["id"]: f for f in pipeline["inputs"] if f.get("id")}


# ---------------------------------------------------------------- spec
def test_the_pipelines_are_registered_in_order():
    assert [p["id"] for p in spec.PIPELINES] == ["esi_ms2_annotation", "ei_ms_annotation", "targeted_method"]
    assert [p["id"] for p in spec.spec()["pipelines"]] == ["esi_ms2_annotation", "ei_ms_annotation", "targeted_method"]
    assert (EI["title"], TARGETED["title"]) == ("EI MS Annotation", "Targeted Method")
    assert TARGETED["one_liner"].startswith("Compound →") and EI["one_liner"].startswith("EI spectrum")
    assert isinstance(TARGETED["examples"], list) and isinstance(EI["examples"], list)


def test_both_workflows_have_nine_stages_with_the_agreed_titles():
    assert [(s["id"], s["title"]) for s in TARGETED["stages"]] == list(enumerate(
        ["Identity", "Platform", "Separation", "Ionization", "Adducts", "Non-MS detection", "Derivatization", "Known methods", "Theoretical m/z"], start=1))
    assert [(s["id"], s["title"]) for s in EI["stages"]] == list(enumerate(
        ["Spectrum cleanup", "EI library search", "Derivatization markers", "Candidates", "Derivatives", "NEIMS prediction",
         "Similarity and rerank", "Kovats RI", "Answer"], start=1))
    assert all(s["modules"] == [] for p in (TARGETED, EI) for s in p["stages"])


def test_the_targeted_method_form():
    f = fields(TARGETED)
    assert list(f) == ["structure", "ph", "additive", "instrument"]
    assert f["structure"]["type"] == "structure" and f["structure"]["required"] and f["structure"]["label"] == "Compound"
    assert f["structure"]["network_modes"] == ["Name / ID / Formula"]                       # a name or ID asks before using the internet
    assert [f[k].get("section") for k in ("ph", "additive", "instrument")] == ["Conditions"] * 3 and "section" not in f["structure"]
    assert f["ph"]["type"] == "number" and f["ph"]["value"] is None and f["ph"]["placeholder"] == "2.7" and not f["ph"].get("required")
    assert f["additive"]["options"] == list(adduct_rules.ADDITIVES) == spec.PIPELINE_ADDITIVES and f["additive"]["value"] == "Formic acid"
    assert f["instrument"]["options"] == ["Any", "LC-MS", "GC-MS"] == spec.TARGETED_INSTRUMENTS and f["instrument"]["value"] == "Any"
    assert tm.DEFAULT_PH == 2.7 and tm.DEFAULT_ADDITIVE == f["additive"]["value"]


def test_the_ei_ms_form():
    f = fields(EI)
    assert list(f) == ["spectrum", "derivatization", "molecular_ion", "candidates", "alkanes", "sample_rt"]
    assert f["spectrum"]["type"] == "spectrum" and f["spectrum"]["required"] and f["spectrum"]["network_sources"] == ["USI"]
    assert f["derivatization"]["options"] == ["Auto (from marker ions)", "None", "TMS", "MeOx + TMS", "TBDMS"] == spec.EI_DERIVATIZATIONS
    assert f["derivatization"]["value"] == "Auto (from marker ions)"
    assert f["molecular_ion"]["type"] == "number" and f["candidates"]["type"] == "textarea"
    assert f["alkanes"]["type"] == "table" and f["alkanes"]["columns"] == ["Carbon number", "RT"] and f["sample_rt"]["type"] == "number"
    assert {f[k]["section"] for k in ("alkanes", "sample_rt")} == {"Retention index"}
    assert not any(f[k].get("required") or f[k].get("required_when") for k in ("alkanes", "sample_rt"))     # the pair is enforced by parse_inputs (the UI would star both)


def test_the_form_ids_are_unique_and_every_field_has_a_label():
    for p in (TARGETED, EI):
        ids = [f["id"] for f in p["inputs"] if f.get("id")]
        assert len(ids) == len(set(ids)), p["id"]
        assert all(f.get("label") for f in p["inputs"]), p["id"]


@pytest.mark.parametrize("pipeline_id, module", NEW_WORKFLOWS)
def test_the_step_titles_come_from_the_spec(pipeline_id, module):
    stages = next(p for p in spec.PIPELINES if p["id"] == pipeline_id)["stages"]
    assert module.PIPELINE_ID == pipeline_id and module.STEP_TITLES == {s["id"]: s["title"] for s in stages} and len(module.STEP_TITLES) == 9


def test_importing_the_workflows_loads_no_engine():
    code = ("import sys; import mzdesk.spec, mzdesk.workflows.targeted_method, mzdesk.workflows.ei_ms_annotation; "
            "bad = [m for m in sys.modules if m.split('.')[0] == 'torch' or m.startswith(('mzdesk.annotation', 'mzdesk.platform_fit', 'mzdesk.calculator', "
            "'mzdesk.ms2.library_index', 'mzdesk.ms2.neims_engine'))]; print(bad); sys.exit(1 if bad else 0)")
    run = subprocess.run([sys.executable, "-c", code], cwd=REPO, capture_output=True, text=True, timeout=120)
    assert run.returncode == 0, run.stdout + run.stderr


# ---------------------------------------------------------------- Targeted Method: parse_inputs
def test_targeted_parse_inputs_accepts_a_valid_form():
    q = tm.parse_inputs(targeted_form())
    assert isinstance(q, tm.Query) and q.smiles == "CCO" and q.name == ""
    assert (q.ph, q.ph_given, q.additive, q.instrument, q.warnings) == (2.7, False, "Formic acid", "Any", [])
    q = tm.parse_inputs(targeted_form(ph="6.8", additive="Ammonium acetate", instrument="LC-MS"))
    assert (q.ph, q.ph_given, q.additive, q.instrument) == (6.8, True, "Ammonium acetate", "LC-MS")
    q = tm.parse_inputs({"structure": {"mode": "SMILES / InChI", "value": "CCO"}})                 # the conditions are all optional
    assert (q.ph, q.additive, q.instrument) == (2.7, "Formic acid", "Any")


def test_targeted_parse_inputs_resolves_a_name_and_keeps_what_was_typed(monkeypatch):
    from mzdesk.identity import resolve
    monkeypatch.setattr(resolve, "to_smiles", lambda text: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    q = tm.parse_inputs(targeted_form(structure={"mode": "Name / ID / Formula", "value": "caffeine"}))
    assert q.smiles == "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" and q.name == "caffeine"


@pytest.mark.parametrize("over, message", [
    ({"structure": {"mode": "SMILES / InChI", "value": ""}}, "Compound is required"),
    ({"structure": None}, "Compound is required"),
    ({"structure": {"mode": "Draw", "value": ""}}, "draw a structure first"),
    ({"structure": "CCO"}, "invalid input"),
    ({"structure": {"mode": "SMILES / InChI", "value": "C1CC"}}, "invalid SMILES"),
    ({"structure": {"mode": "Name / ID / Formula", "value": "C8H10N4O2"}}, "molecular formula"),
    ({"ph": "acidic"}, "must be a number"), ({"ph": 15}, "at most 14"), ({"ph": -1}, "at least 0"),
    ({"additive": "Water"}, "additive must be"), ({"instrument": "CE-MS"}, "Instrument must be")])
def test_targeted_parse_inputs_rejects_bad_forms(over, message):
    with pytest.raises(InputError, match=message):
        tm.parse_inputs(targeted_form(**over))


def test_targeted_parse_inputs_rejects_an_empty_form():
    with pytest.raises(InputError, match="Compound is required"):
        tm.parse_inputs({})


# ---------------------------------------------------------------- EI MS Annotation: parse_inputs
def test_ei_parse_inputs_accepts_a_valid_form():
    q = ems.parse_inputs(ei_form())
    assert isinstance(q, ems.Query) and q.ei == [(73.0, 999.0), (147.0, 250.0), (204.0, 120.0), (217.0, 90.0)]
    assert (q.derivatization, q.molecular_ion, q.candidates, q.alkanes, q.sample_rt, q.warnings) == (None, None, [], [], None, [])
    q = ems.parse_inputs({"spectrum": {"source": "Paste", "value": "73 100"}})                      # everything but the spectrum is optional
    assert len(q.ei) == 1 and q.derivatization is None


def test_ei_parse_inputs_reads_the_optional_fields():
    q = ems.parse_inputs(ei_form(derivatization="MeOx + TMS", molecular_ion="319.1", candidates="CCO\n\n# a comment\n  alanine \nCCO\nC(C(=O)O)N",
                                alkanes="C10 10.0\nC8 8.0\nC9 9.1", sample_rt="9.5"))
    assert q.derivatization == "MeOx + TMS" and q.molecular_ion == 319.1
    assert q.candidates == ["CCO", "alanine", "C(C(=O)O)N"]                                         # blank and comment lines dropped, repeats merged
    assert q.alkanes == [(8, 8.0), (9, 9.1), (10, 10.0)] and q.sample_rt == 9.5
    assert ems.parse_inputs(ei_form(derivatization="None")).derivatization == "None"
    assert ems.parse_inputs(ei_form(derivatization="TBDMS")).derivatization == "TBDMS"


def test_ei_parse_inputs_reads_an_msp_record_and_warns_about_extra_spectra():
    msp = "Name: A\nPrecursorMZ: 319\nNum Peaks: 2\n73 999\n147 250\n\nName: B\nNum Peaks: 1\n45 100\n"
    q = ems.parse_inputs(ei_form(spectrum={"source": "Paste", "value": msp}))
    assert q.ei == [(73.0, 999.0), (147.0, 250.0)] and len(q.warnings) == 1 and "2 spectra" in q.warnings[0]


@pytest.mark.parametrize("over, message", [
    ({"spectrum": {"source": "Paste", "value": ""}}, "no valid peaks"),
    ({"spectrum": {"source": "Paste", "value": "no peaks here"}}, "no valid peaks"),
    ({"spectrum": {"source": "USI", "value": ""}}, "enter a USI"),
    ({"spectrum": {"source": "Elsewhere", "value": "73 100"}}, "unknown source"),
    ({"spectrum": "73 100"}, "invalid input"),
    ({"derivatization": "BSTFA"}, "Derivatization must be"),
    ({"molecular_ion": "heavy"}, "must be a number"), ({"molecular_ion": 0}, "must be positive"), ({"molecular_ion": -5}, "must be positive"),
    ({"alkanes": "C8 8.0\nC10 10.0"}, "Sample RT is required"),
    ({"sample_rt": "9.5"}, "alkane series is required"),
    ({"alkanes": "C8 8.0\nC10 10.0", "sample_rt": "late"}, "must be a number"),
    ({"alkanes": "C8 8.0\nC10 10.0", "sample_rt": 0}, "must be positive"),
    ({"alkanes": "C8 8.0", "sample_rt": "9"}, "at least two"),
    ({"alkanes": "C8 8.0\nC10 7.0", "sample_rt": "9"}, "must increase"),
    ({"alkanes": "C8 8.0\nC8 9.0", "sample_rt": "9"}, "once"),
    ({"alkanes": "C8 8.0\nCx nine", "sample_rt": "9"}, "Alkane series: Alkane line 2")])
def test_ei_parse_inputs_rejects_bad_forms(over, message):
    with pytest.raises(InputError, match=message):
        ems.parse_inputs(ei_form(**over))


def test_ei_parse_inputs_rejects_an_empty_form():
    with pytest.raises(InputError, match="EI spectrum"):
        ems.parse_inputs({})


def test_ei_parse_inputs_limits_the_size_of_the_input():
    big = "\n".join(f"{50 + i * 0.01:.2f} 10" for i in range(ems.MAX_INPUT_PEAKS + 1))
    with pytest.raises(InputError, match="at most 5000"):
        ems.parse_inputs(ei_form(spectrum={"source": "Paste", "value": big}))
    with pytest.raises(InputError, match="At most 50 candidate"):
        ems.parse_inputs(ei_form(candidates="\n".join("C" * (i + 1) for i in range(ems.MAX_CANDIDATES + 1))))


# ---------------------------------------------------------------- the step runner and the generic result
STEPS = {1: "One", 2: "Two", 3: "Three", 4: "Four"}


def test_step_runner_records_every_step_and_reports_progress():
    events = []
    steps = wc.StepRunner(STEPS, lambda *a: events.append(a))
    assert [r["status"] for r in steps.records.values()] == ["queued"] * 4 and steps.records[1] == {"id": 1, "title": "One", "status": "queued", "seconds": 0.0, "note": ""}
    block = wc.notes_block("n", "Notes", ["x"])
    assert steps.guarded(1, lambda: (wc.STATUS_DONE, "ok", [block], 42)) == 42
    steps.skip(2, "not needed")
    steps.start(3)
    steps.emit(3, "running", "scored 1/2")
    steps.finish(3, wc.STATUS_DONE, "done", [], seconds=1.25)
    assert [(e[0], e[1]) for e in events] == [(1, "running"), (1, "done"), (2, "skipped"), (3, "running"), (3, "running"), (3, "done")]
    assert events[1][2] == "ok" and events[1][3] == [block] and events[2][2] == "not needed" and events[4][2] == "scored 1/2"
    listed = steps.as_list()
    assert [r["id"] for r in listed] == [1, 2, 3, 4] and [r["status"] for r in listed] == ["done", "skipped", "done", "queued"]
    assert listed[0]["note"] == "ok" and listed[0]["seconds"] >= 0 and listed[2]["seconds"] == 1.25 and listed[3]["note"] == ""


def test_a_failing_step_does_not_stop_the_run(caplog):
    steps = wc.StepRunner(STEPS)

    def boom():
        raise RuntimeError("engine exploded")
    with caplog.at_level(logging.CRITICAL):
        assert steps.guarded(1, boom) is None
    assert steps.records[1]["status"] == "failed" and steps.records[1]["note"] == "One failed: RuntimeError: engine exploded"
    assert steps.guarded(2, lambda: (wc.STATUS_DONE, "fine", [], "next")) == "next"             # the next step still runs
    assert [r["id"] for r in steps.failed()] == [1]


def test_a_broken_progress_listener_does_not_stop_the_run(caplog):
    def listener(*args):
        raise ValueError("listener bug")
    steps = wc.StepRunner(STEPS, listener)
    with caplog.at_level(logging.CRITICAL):
        assert steps.guarded(1, lambda: (wc.STATUS_DONE, "ok", [], 1)) == 1
    assert steps.records[1]["status"] == "done"


def test_make_result_has_the_generic_shape():
    steps = wc.StepRunner(STEPS)
    steps.skip(1, "n/a")
    table = wc.table_block("t", "A table", ["Name", "Fit"], [["Formic acid", blk.badge("good", "good")], ["Acetic acid", "possible"]])
    assert table == {"type": "table", "id": "t", "title": "A table", "columns": ["Name", "Fit"],
                     "data": blk.table_data([["Formic acid", blk.badge("good", "good")], ["Acetic acid", "possible"]], count=2)}
    result = wc.make_result("demo", steps, wc.make_answer("Done.", "A", "Level A: strong", [["Platform", "LC-MS"], ["n", 3]]), [table], ["careful"])
    assert list(result) == ["pipeline", "steps", "answer", "blocks", "warnings"] and result["pipeline"] == "demo"
    assert result["steps"] == steps.as_list() and len(result["steps"]) == 4
    assert result["answer"] == {"statement": "Done.", "level": "A", "level_label": "Level A: strong", "facts": [["Platform", "LC-MS"], ["n", "3"]]}
    assert result["blocks"] == [table] and result["warnings"] == ["careful"]
    plain = wc.make_result("demo", steps.as_list(), wc.make_answer("Just text."))                # a list of records works too, no blocks or warnings
    assert plain["answer"] == {"statement": "Just text.", "level": None, "level_label": None, "facts": []} and plain["blocks"] == [] and plain["warnings"] == []


def test_make_result_drops_a_malformed_or_repeated_block(caplog):
    good = wc.notes_block("n", "Notes", ["x"])
    bad = [{"type": "notes", "title": "no id", "data": {"items": []}}, {"id": "x", "data": {}}, {"type": "table", "id": "y"}, "text", good]
    with caplog.at_level(logging.CRITICAL):
        result = wc.make_result("demo", [], wc.make_answer("s"), bad + [dict(good)], [])
    assert result["blocks"] == [good] and len(result["warnings"]) == 5


# ---------------------------------------------------------------- the stubs through the job API
def wait_done(job_id, timeout=15):
    end = time.time() + timeout
    while time.time() < end:
        body = client.get(f"/api/pipelines/jobs/{job_id}").json()
        if body["status"] not in ("running", "queued"):
            return body
        time.sleep(0.02)
    raise AssertionError("job did not finish")


def test_every_new_workflow_is_registered_with_the_api():
    for pipeline_id, module in NEW_WORKFLOWS:
        assert app_module._PIPELINE_RUNNERS[pipeline_id] is module


def test_ei_ms_annotation_runs_end_to_end_through_the_job_api(monkeypatch, tmp_path):
    """The real workflow on a bare spectrum (no library installed, no candidates): the generic result, every step reports, nothing is a stub any more.
    The workflow itself is tested in tests/test_workflow_ei_ms.py."""
    from mzdesk.ms2 import library_index as lx
    monkeypatch.setenv("MZDESK_LIBRARY_DIR", str(tmp_path / "no_libs"))
    lx.clear_cache()
    pipeline_id = "ei_ms_annotation"
    r = client.post(f"/api/pipelines/{pipeline_id}/run", json=ei_form())
    assert r.status_code == 200 and r.json()["job_id"] and r.json()["status"] in ("running", "queued")
    body = wait_done(r.json()["job_id"])
    assert body["status"] == "done" and body["error"] is None
    result = body["result"]
    assert result["pipeline"] == pipeline_id and set(result) == {"pipeline", "steps", "answer", "blocks", "warnings"}
    assert len(result["steps"]) == 9 and [s["id"] for s in result["steps"]] == list(range(1, 10))
    assert [s["title"] for s in result["steps"]] == [s["title"] for s in next(p for p in spec.PIPELINES if p["id"] == pipeline_id)["stages"]]
    assert all(s["status"] in ("done", "skipped") and "not implemented" not in s["note"] and s["seconds"] >= 0 for s in result["steps"])
    assert set(result["answer"]) == {"statement", "level", "level_label", "facts"} and result["answer"]["statement"].startswith("No candidate.")
    assert result["answer"]["level"] is None and all(isinstance(f, list) and len(f) == 2 for f in result["answer"]["facts"])
    assert all(b["id"] and b["type"] and isinstance(b["data"], dict) for b in result["blocks"])
    assert [s["title"] for s in body["steps"]] == [s["title"] for s in result["steps"]]
    assert all(set(s) == {"id", "title", "status", "seconds", "note", "blocks", "expected_s"} for s in body["steps"])
    lx.clear_cache()


@pytest.mark.parametrize("pipeline_id, form, message", [
    ("targeted_method", {}, "Compound is required"),
    ("targeted_method", targeted_form(ph=20), "at most 14"),
    ("ei_ms_annotation", {}, "EI spectrum"),
    ("ei_ms_annotation", ei_form(alkanes="C8 8\nC9 9"), "Sample RT is required"),
    ("ei_ms_annotation", ei_form(sample_rt="9"), "alkane series is required")])
def test_a_bad_form_is_a_400_and_starts_no_job(pipeline_id, form, message):
    r = client.post(f"/api/pipelines/{pipeline_id}/run", json=form)
    assert r.status_code == 400 and message in r.json()["detail"] and not app_module._jobs


def test_a_failing_step_of_a_new_workflow_is_reported_and_the_run_goes_on(monkeypatch):
    def run_query(q, progress=None):
        steps = wc.StepRunner(tm.STEP_TITLES, progress)

        def boom():
            raise RuntimeError("lookup failed")
        steps.guarded(1, boom)
        for n in range(2, 10):
            steps.guarded(n, lambda: (wc.STATUS_DONE, "ok", [wc.notes_block("n", "Notes", ["x"])], None))
        return wc.make_result(tm.PIPELINE_ID, steps, wc.make_answer("Partial."), [], [f"Step {r['id']} failed: {r['note']}" for r in steps.failed()])
    monkeypatch.setattr(tm, "run_query", run_query)
    body = wait_done(client.post("/api/pipelines/targeted_method/run", json=targeted_form()).json()["job_id"])
    assert body["status"] == "done" and [s["status"] for s in body["steps"]] == ["failed"] + ["done"] * 8
    assert body["steps"][0]["note"] == "Identity failed: RuntimeError: lookup failed"
    assert body["result"]["warnings"] == ["Step 1 failed: Identity failed: RuntimeError: lookup failed"]
    assert body["steps"][1]["blocks"][0]["id"] == "n"                                               # a step's blocks reach the job while it runs
