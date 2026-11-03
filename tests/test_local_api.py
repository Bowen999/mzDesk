"""Local version: /api/local/* and /api/update. Response shapes (the frontend depends on them), the Host / Origin guard, unknown ids, mzdesk.com
(hosted) answers, and an install job through the API with a fake component. Nothing is installed; HOME is a temporary folder."""
import dataclasses
import os
import sys
import time

import pytest
from fastapi.testclient import TestClient

from mzdesk.api.app import app
from mzdesk.local import components, jobs, prereq


@pytest.fixture(autouse=True)
def home(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("USERPROFILE", str(tmp_path))
    for key in list(os.environ):
        if key.startswith("MZDESK_"):
            monkeypatch.delenv(key)
    monkeypatch.setenv("MZDESK_NO_UPDATE_CHECK", "1")
    monkeypatch.setattr(prereq, "docker_state", lambda image, use_cache=True: "ready")
    monkeypatch.setattr(jobs, "STORE", jobs.JobStore())
    return tmp_path


@pytest.fixture
def client():
    return TestClient(app, base_url="http://localhost:8731")


def fake_registry(monkeypatch, run, **flags):
    state = {}
    reg = {cid: dataclasses.replace(c, check=lambda v=flags.get(cid, False): v) for cid, c in components.REGISTRY.items()}

    def runner(ctx):
        run(ctx)
        state["done"] = True
    reg["molgpka"] = dataclasses.replace(reg["molgpka"], check=lambda: state.get("done", False), run=runner)
    monkeypatch.setattr(components, "REGISTRY", reg)
    return state


def wait_done(client, job_id, seconds=30):
    deadline = time.monotonic() + seconds
    while time.monotonic() < deadline:
        body = client.get(f"/api/local/jobs/{job_id}").json()
        if body["status"] != "running":
            return body
        time.sleep(0.05)
    raise AssertionError("job did not finish")


# ---------------------------------------------------------------- status
def test_status_shape(client):
    body = client.get("/api/local/status").json()
    assert body["local"] is True and body["platform"] in ("darwin", "linux", "win32")
    assert [c["id"] for c in body["components"]] == ["libraries", "msfinder", "frigid", "cfmid", "neims", "molgpka", "ccs", "repo"]
    for c in body["components"]:
        assert set(c) == {"id", "title", "size", "installed", "supported", "ask", "note"}
        assert isinstance(c["installed"], bool) and isinstance(c["supported"], bool) and isinstance(c["ask"], bool) and isinstance(c["note"], str)
    assert next(c for c in body["components"] if c["id"] == "libraries")["size"] == "2.3 GB"


def test_hosted_answers_local_false_and_404(client, monkeypatch):
    monkeypatch.setenv("MZDESK_HOSTED", "1")
    assert client.get("/api/local/status").json() == {"local": False}
    assert client.post("/api/local/needs", json={"target": "module:ccs"}).status_code == 404
    assert client.post("/api/local/install", json={"ids": ["ccs"]}).status_code == 404
    assert client.get("/api/local/jobs/abc").status_code == 404
    assert client.post("/api/local/jobs/abc/cancel").status_code == 404
    body = client.get("/api/update").json()
    assert body["newer"] is False and body["latest"] is None                        # no check on mzdesk.com


# ---------------------------------------------------------------- needs
def test_needs_shape_and_ask(client, monkeypatch):
    body = client.post("/api/local/needs", json={"target": "workflow:esi_ms2_annotation",
                                                 "inputs": {"libraries": ["MassBank"], "run_frigid": True}}).json()
    assert set(body) == {"missing", "blocked", "unsupported"}
    for item in body["missing"]:
        assert set(item) == {"id", "title", "size", "ask", "why"} and isinstance(item["ask"], bool)
    ask = {m["id"]: m["ask"] for m in body["missing"]}
    assert ask["frigid"] is True and ask.get("libraries", False) is False and ask["msfinder"] is False
    monkeypatch.setattr(prereq, "docker_state", lambda image, use_cache=True: "no-docker")
    blocked = client.post("/api/local/needs", json={"target": "module:ms2_esi", "inputs": {}}).json()["blocked"]
    assert blocked == [{"id": "cfmid", "title": components.REGISTRY["cfmid"].title, "message": blocked[0]["message"], "link": "https://docs.docker.com/get-docker/"}]


def test_needs_rejects_bad_requests(client):
    assert client.post("/api/local/needs", json={"target": "module:nope"}).status_code == 404
    assert client.post("/api/local/needs", json={"target": 5}).status_code == 400
    assert client.post("/api/local/needs", json={}).status_code == 400
    assert client.post("/api/local/needs", json={"target": "module:ccs", "inputs": "text"}).status_code == 200      # odd inputs are ignored


# ---------------------------------------------------------------- guard
@pytest.mark.parametrize("path,body", [("/api/local/needs", {"target": "module:ccs"}), ("/api/local/install", {"ids": ["molgpka"]}),
                                       ("/api/local/jobs/x/cancel", {})])
def test_posts_refuse_a_foreign_host(client, path, body):
    for host in ("evil.example", "evil.example:8731", "localhost.evil.example", "192.168.1.5:8731", "127.0.0.1.evil.example"):
        assert client.post(path, json=body, headers={"Host": host}).status_code == 403, host
    for host in ("localhost", "localhost:8731", "127.0.0.1:8599", "[::1]:8731", "LOCALHOST:9"):
        assert client.post(path, json=body, headers={"Host": host}).status_code != 403, host


@pytest.mark.parametrize("origin", ["http://evil.example", "http://localhost:3000", "https://localhost:8731", "null", "http://127.0.0.1:8731",
                                    "http://localhost.evil.example:8731"])
def test_posts_refuse_a_cross_site_origin(client, origin):
    for path, body in (("/api/local/needs", {"target": "module:ccs"}), ("/api/local/install", {"ids": ["molgpka"]}), ("/api/local/jobs/x/cancel", {})):
        assert client.post(path, json=body, headers={"Origin": origin}).status_code == 403, (path, origin)
    assert jobs.STORE.running() is None                                              # nothing started


def test_same_origin_posts_pass(client):
    ok = {"Origin": "http://localhost:8731"}
    assert client.post("/api/local/needs", json={"target": "module:ccs"}, headers=ok).status_code == 200
    ok = {"Origin": "http://127.0.0.1:5175", "Host": "127.0.0.1:5175"}                 # the vite dev server keeps the browser's Host
    assert client.post("/api/local/needs", json={"target": "module:ccs"}, headers=ok).status_code == 200
    assert client.post("/api/local/needs", json={"target": "module:ccs"}).status_code == 200      # curl: no Origin


# ---------------------------------------------------------------- install
def test_install_rejects_bad_ids(client):
    assert client.post("/api/local/install", json={"ids": ["nope"]}).status_code == 404
    assert client.post("/api/local/install", json={"ids": ["libraries", "../../etc/passwd"]}).status_code == 404
    for bad in ({}, {"ids": []}, {"ids": "ccs"}, {"ids": [1]}, {"ids": [["ccs"]]}):
        assert client.post("/api/local/install", json=bad).status_code == 400, bad
    assert jobs.STORE.running() is None


def test_install_job_through_the_api(client, monkeypatch):
    def run(ctx):
        ctx.run([sys.executable, "-u", "-c", "print('hello from the setup')"], what="fake setup")
    state = fake_registry(monkeypatch, run)
    first = client.post("/api/local/install", json={"ids": ["molgpka"]})
    assert first.status_code == 200 and set(first.json()) == {"job_id"}
    body = wait_done(client, first.json()["job_id"])
    assert body["status"] == "done" and body["current"] is None and body["error"] == "" and "hello from the setup" in body["log"]
    assert body["components"] == [{"id": "molgpka", "title": "MolGpKa pKa model", "status": "done", "detail": "installed", "step": "", "progress": None,
                                   "got_mb": None, "total_mb": None}]
    status = {c["id"]: c for c in client.get("/api/local/status").json()["components"]}
    assert status["molgpka"]["installed"] is True and state["done"]
    assert client.get("/api/local/jobs/unknown").status_code == 404


def test_a_second_install_call_returns_the_running_job(client, monkeypatch):
    import threading
    gate = threading.Event()
    fake_registry(monkeypatch, lambda ctx: gate.wait(20))
    a = client.post("/api/local/install", json={"ids": ["molgpka"]}).json()["job_id"]
    b = client.post("/api/local/install", json={"ids": ["ccs"]}).json()["job_id"]
    assert a == b
    running = client.get(f"/api/local/jobs/{a}").json()
    assert running["status"] == "running" and running["current"] == "molgpka" and running["components"][0]["status"] == "running"
    gate.set()
    assert wait_done(client, a)["status"] == "done"


def test_failed_install_reports_the_error(client, monkeypatch):
    fake_registry(monkeypatch, lambda ctx: ctx.run([sys.executable, "-c", "print('network is down'); raise SystemExit(1)"], what="fake setup"))
    body = wait_done(client, client.post("/api/local/install", json={"ids": ["molgpka"]}).json()["job_id"])
    assert body["status"] == "error" and "fake setup failed" in body["error"] and "network is down" in body["error"]
    assert body["components"][0]["status"] == "error"


def test_cancel_through_the_api(client, monkeypatch):
    fake_registry(monkeypatch, lambda ctx: ctx.run([sys.executable, "-u", "-c", "print('go', flush=True); import time; time.sleep(60)"], what="slow"))
    job_id = client.post("/api/local/install", json={"ids": ["molgpka"]}).json()["job_id"]
    deadline = time.monotonic() + 20
    while "go" not in client.get(f"/api/local/jobs/{job_id}").json()["log"] and time.monotonic() < deadline:
        time.sleep(0.05)
    assert client.post(f"/api/local/jobs/{job_id}/cancel").json() == {"ok": True}
    body = wait_done(client, job_id)
    assert body["status"] == "cancelled"
    assert client.post("/api/local/jobs/unknown/cancel").status_code == 404


def test_unsupported_component_is_refused_with_its_reason(client, monkeypatch):
    import dataclasses
    monkeypatch.setattr(components, "REGISTRY", {k: dataclasses.replace(c, platforms=components.POSIX_ONLY) if k == "msfinder" else c
                                                 for k, c in components.REGISTRY.items()})
    monkeypatch.setattr(components, "current_platform", lambda: "win32")
    reply = client.post("/api/local/install", json={"ids": ["msfinder"]})
    assert reply.status_code == 400 and "macOS / Linux only for now" in reply.json()["detail"]


# ---------------------------------------------------------------- /api/update
def test_update_endpoint_shape_and_no_check_switch(client):
    body = client.get("/api/update").json()
    assert set(body) == {"current", "latest", "newer", "released", "notes", "command", "checked_at"}
    from mzdesk import __version__
    assert body["current"] == __version__ and body["latest"] is None and body["newer"] is False
    assert body["notes"] == "https://mzdesk.com/docs/changelog" and body["command"] == "mzdesk update"
