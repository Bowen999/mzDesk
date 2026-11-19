"""Crash reports and the maintainer's reports page: what is recorded, what the user sees, who may read it."""
import json

import pytest
from fastapi.testclient import TestClient

from mzdesk.api import admin, app as app_module, feedback, reports, stats
from mzdesk.api.app import app

client = TestClient(app, raise_server_exceptions=False)
PASSWORD = "correct horse battery staple"


@pytest.fixture(autouse=True)
def admin_password(monkeypatch):
    salt, digest = admin.hash_password(PASSWORD)
    monkeypatch.setenv("MZDESK_ADMIN_SALT", salt)
    monkeypatch.setenv("MZDESK_ADMIN_HASH", digest)
    monkeypatch.setattr(admin.time, "sleep", lambda s: None)
    monkeypatch.setattr(reports, "forward_async", lambda ref: None)
    feedback._CRASH_HITS.clear()
    feedback._HITS.clear()
    yield
    feedback._CRASH_HITS.clear()
    feedback._HITS.clear()


def _boom(inputs):
    raise RuntimeError("engine exploded at /var/lib/mzdesk/models/neims/weights.h5")


def test_a_crashing_module_gives_the_user_a_reference_and_no_server_folders(monkeypatch):
    monkeypatch.setitem(app_module.RUNNERS, "properties", _boom)
    r = client.post("/api/modules/properties/run", json={"structure": {"mode": "SMILES / InChI", "value": "CCO"}})
    assert r.status_code == 500
    detail = r.json()["detail"]
    assert detail.startswith("properties failed: RuntimeError: engine exploded at weights.h5 (ref E-") and "/var/lib" not in detail
    ref = detail.split("ref ")[1].rstrip(")")
    [row] = reports.listing()
    assert row["id"] == ref and row["kind"] == "crash" and row["title"].startswith("module properties: RuntimeError")
    full = reports.detail(ref)
    assert "Traceback" in full["body"] and "_boom" in full["body"]
    assert "structure (" in full["detail"] and "CCO" not in json.dumps(full)          # field names and sizes, never the values


def test_the_same_crash_again_counts_up_the_first_report_and_is_forwarded_once_an_hour(monkeypatch):
    monkeypatch.setitem(app_module.RUNNERS, "properties", _boom)
    sent = []
    monkeypatch.setattr(reports, "forward_async", sent.append)
    refs = {client.post("/api/modules/properties/run", json={"structure": {"value": "CCO"}}).json()["detail"].split("ref ")[1].rstrip(")") for _ in range(3)}
    assert len(refs) == 1
    [row] = reports.listing()
    assert row["count"] == 3 and sent == [row["id"]]


def test_a_different_crash_is_a_different_report(monkeypatch):
    def other(inputs):
        raise ValueError("something else")
    monkeypatch.setitem(app_module.RUNNERS, "properties", _boom)
    monkeypatch.setitem(app_module.RUNNERS, "identity", other)
    client.post("/api/modules/properties/run", json={})
    client.post("/api/modules/identity/run", json={})
    assert len(reports.listing()) == 2


def test_a_route_that_raises_answers_with_a_reference(monkeypatch):
    monkeypatch.setattr(app_module, "_depict", lambda smiles: (_ for _ in ()).throw(KeyError("x")))
    r = client.get("/api/depict", params={"smiles": "CCO"})
    assert r.status_code == 500 and r.json()["detail"].startswith("Something went wrong (ref E-")
    assert reports.listing()[0]["title"].startswith("GET /api/depict")


def _failing_step(number):
    import logging
    try:
        raise OSError("disk gone")
    except OSError:
        logging.getLogger("mzdesk.annotation.pipeline").exception("step %s failed", number)


def test_log_exception_anywhere_in_mzdesk_files_a_crash_report():
    _failing_step(3)
    [row] = reports.listing()
    assert row["kind"] == "crash" and row["title"] == "mzdesk.annotation.pipeline: step %s failed: OSError: disk gone"
    _failing_step(4)
    assert len(reports.listing()) == 1 and reports.listing()[0]["count"] == 2          # the step number does not make it a new crash


def test_a_page_that_crashed_in_the_browser_is_reported():
    stack = "TypeError: Cannot read properties of undefined (reading 'map')\n    at Table (assets/index-abc.js:1:2345)\n    at div"
    r = client.post("/api/crash", json={"message": "Cannot read properties of undefined (reading 'map')", "stack": stack, "page": "/modules/ms1_search"},
                    headers={"user-agent": "Mozilla/5.0 Test"})
    assert r.status_code == 200 and r.json()["ref"].startswith("E-")
    [row] = reports.listing()
    assert row["kind"] == "crash" and row["page"] == "/modules/ms1_search"
    assert reports.detail(row["id"])["browser"] == "Mozilla/5.0 Test"
    client.post("/api/crash", json={"message": "Cannot read properties of undefined (reading 'map')", "stack": stack, "page": "/modules/ms1_search"})
    assert reports.listing()[0]["count"] == 2


def test_browser_crash_reports_are_capped():
    assert client.post("/api/crash", json={"message": ""}).status_code == 400
    assert client.post("/api/crash", content=b"[1]", headers={"content-type": "application/json"}).status_code == 400
    big = client.post("/api/crash", json={"message": "m" * 5000, "stack": "s" * 90_000, "page": "p" * 900}).json()
    detail = reports.detail(big["ref"])
    assert len(detail["title"]) <= 140 and len(detail["body"]) <= 6000 and len(detail["page"]) <= 200
    for i in range(feedback.CRASH_HOURLY_LIMIT):
        client.post("/api/crash", json={"message": f"distinct {i}", "stack": f"at {i}"})
    assert client.post("/api/crash", json={"message": "one more", "stack": "x"}).status_code == 429


def test_the_reports_pages_need_the_admin_password():
    ref = reports.add_feedback("Hello there, something broke.", "/modules/x", "UA", [("a.png", "image/png", b"\x89PNG")])
    assert client.post("/api/reports", json={}).status_code == 401
    assert client.post("/api/reports", json={"password": "wrong"}).status_code == 401
    assert client.post(f"/api/reports/{ref}", json={"password": "wrong"}).status_code == 401
    assert client.post(f"/api/reports/{ref}/delete", json={}).status_code == 401
    assert reports.detail(ref) is not None

    listing = client.post("/api/reports", json={"password": PASSWORD}).json()
    assert [r["id"] for r in listing["reports"]] == [ref] and listing["channels"] == []
    one = client.post(f"/api/reports/{ref}", json={"password": PASSWORD}).json()
    assert one["body"] == "Hello there, something broke." and one["images"][0]["src"].startswith("data:image/png;base64,")
    assert client.post("/api/reports/F-NOPE00", json={"password": PASSWORD}).status_code == 404
    assert client.post(f"/api/reports/{ref}/delete", json={"password": PASSWORD}).json() == {"ok": True}
    assert reports.listing() == []


def test_wrong_passwords_are_throttled_per_address():
    for _ in range(admin.MAX_WRONG):
        assert client.post("/api/reports", json={"password": "nope"}).status_code == 401
    assert client.post("/api/reports", json={"password": "nope"}).status_code == 429
    assert client.post("/api/reports", json={"password": PASSWORD}).status_code == 429          # even the right one waits
    assert client.post("/api/stats", json={"password": PASSWORD}).status_code == 429


def test_the_server_password_replaces_the_one_in_the_page_lock(monkeypatch):
    assert stats.check_password(PASSWORD) and not stats.check_password("")
    monkeypatch.delenv("MZDESK_ADMIN_SALT")
    monkeypatch.delenv("MZDESK_ADMIN_HASH")
    assert not stats.check_password(PASSWORD)                    # without the server's pair the page-lock pair applies (tests/test_api.py)
    assert admin.hash_password("x", "00" * 16) == admin.hash_password("x", "00" * 16) != admin.hash_password("x")


def test_old_reports_go_first(monkeypatch):
    monkeypatch.setattr(reports, "KEEP", 3)
    refs = [reports.add_feedback(f"message number {i}", "/", "", [("a.png", "image/png", b"x")]) for i in range(5)]
    assert {r["id"] for r in reports.listing()} == set(refs[2:])
    assert reports.detail(refs[0]) is None
    assert reports._db().execute("SELECT COUNT(*) FROM images").fetchone()[0] == 3


def test_screenshots_have_a_total_budget_and_the_oldest_reports_lose_theirs_first(monkeypatch):
    monkeypatch.setattr(reports, "MAX_IMAGE_BYTES", 2500)
    refs = [reports.add_feedback(f"message number {i}", "/", "", [("a.png", "image/png", b"x" * 1000)]) for i in range(5)]
    assert [reports.detail(r)["images"] != [] for r in refs] == [False, False, False, True, True]       # 2 x 1000 B fit; the text of all five stays
    assert all(reports.detail(r)["body"] for r in refs)


def test_crash_emails_have_an_overall_hourly_budget(monkeypatch):
    monkeypatch.setattr(reports, "CRASH_FORWARDS_PER_HOUR", 3)
    sent = []
    monkeypatch.setattr(reports, "forward_async", sent.append)
    for i in range(6):
        client.post("/api/crash", json={"message": f"a different crash {i}", "stack": f"at place {i}"})
    assert len(reports.listing()) == 6 and len(sent) == 3                                              # all stored, only three forwarded


def test_the_password_script_hashes_what_it_reads_from_stdin(monkeypatch, capsys):
    import io
    monkeypatch.setattr("sys.stdin", io.StringIO("a typed password\n"))
    admin.main(["--stdin"])
    out = dict(line.split("=", 1) for line in capsys.readouterr().out.split())
    monkeypatch.setenv("MZDESK_ADMIN_SALT", out["MZDESK_ADMIN_SALT"])
    monkeypatch.setenv("MZDESK_ADMIN_HASH", out["MZDESK_ADMIN_HASH"])
    assert stats.check_password("a typed password") and not stats.check_password("a typed password ")        # the line end is not part of it
    monkeypatch.setattr("sys.stdin", io.StringIO("abc\n"))                        # short is allowed, with a note
    admin.main(["--stdin"])
    assert "MZDESK_ADMIN_HASH=" in capsys.readouterr().out
    monkeypatch.setattr("sys.stdin", io.StringIO("\n"))
    with pytest.raises(SystemExit):                                              # but not nothing
        admin.main(["--stdin"])
