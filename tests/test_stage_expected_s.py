"""expected_s: the typical seconds of every pipeline stage (mzdesk/spec.py stage(...)). The progress bar of the workflow page gives each stage a share
proportional to it and fills a running stage at a pace derived from it, so every stage needs a positive number and the job API has to hand it on.
The numbers themselves come from measurements (see the comments next to the stages in spec.py); only their shape is tested here.
"""
import math
import threading
import time
import types

import pytest
from fastapi.testclient import TestClient

from mzdesk import spec
from mzdesk.api import app as app_module

client = TestClient(app_module.app)
PIPELINE_IDS = [p["id"] for p in spec.PIPELINES]


@pytest.fixture(autouse=True)
def clean_jobs():
    app_module._jobs.clear()
    yield
    app_module._jobs.clear()


def test_stage_has_a_positive_default():
    s = spec.stage(1, "Anything")
    assert s["expected_s"] > 0 and s["id"] == 1 and s["title"] == "Anything" and s["modules"] == []
    assert spec.stage(2, "Slow", 30)["expected_s"] == 30.0 and isinstance(spec.stage(2, "Slow", 30)["expected_s"], float)


def test_three_pipelines_are_registered():
    assert PIPELINE_IDS == ["esi_ms2_annotation", "ei_ms_annotation", "targeted_method"]


@pytest.mark.parametrize("stages", [p["stages"] for p in spec.PIPELINES] + [p["stages"] for p in spec.spec()["pipelines"]],
                         ids=PIPELINE_IDS * 2)
def test_every_stage_of_every_pipeline_has_a_positive_expected_s(stages):
    assert stages
    for s in stages:
        assert isinstance(s["expected_s"], float) and math.isfinite(s["expected_s"]) and 0 < s["expected_s"] <= 600, s


def test_the_spec_numbers_are_not_all_the_default():
    """A pipeline whose stages all keep the default 1.0 would draw equal slices: each one needs at least one real estimate."""
    for p in spec.PIPELINES:
        values = [s["expected_s"] for s in p["stages"]]
        assert len(set(values)) > 1, f"{p['id']}: every stage has the same expected_s"


@pytest.mark.parametrize("pipeline_id", PIPELINE_IDS)
def test_the_job_api_returns_expected_s_for_every_step(pipeline_id, monkeypatch):
    stages = next(p for p in spec.PIPELINES if p["id"] == pipeline_id)["stages"]
    midway, go = threading.Event(), threading.Event()

    def run_query(query, progress):
        progress(1, "running", None, None)
        midway.set()
        go.wait(5)
        return {"steps": [{"id": s["id"], "title": s["title"], "status": "done", "seconds": 0.5, "note": ""} for s in stages]}

    stub = types.SimpleNamespace(PIPELINE_ID=pipeline_id, parse_inputs=lambda inputs: object(), run_query=run_query)
    monkeypatch.setitem(app_module._PIPELINE_RUNNERS, pipeline_id, stub)
    r = client.post(f"/api/pipelines/{pipeline_id}/run", json={})
    assert r.status_code == 200
    job_id = r.json()["job_id"]
    assert midway.wait(5)
    running = client.get(f"/api/pipelines/jobs/{job_id}").json()
    assert [(s["id"], s["expected_s"]) for s in running["steps"]] == [(s["id"], s["expected_s"]) for s in stages]
    assert running["steps"][0]["status"] == "running" and all(s["expected_s"] > 0 for s in running["steps"])
    go.set()
    end = time.time() + 10
    while time.time() < end:
        body = client.get(f"/api/pipelines/jobs/{job_id}").json()
        if body["status"] not in ("running", "queued"):
            break
        time.sleep(0.02)
    assert body["status"] == "done"
    assert [(s["id"], s["expected_s"], s["status"], s["seconds"]) for s in body["steps"]] == [(s["id"], s["expected_s"], "done", 0.5) for s in stages]
