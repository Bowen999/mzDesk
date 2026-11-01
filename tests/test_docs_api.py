"""The Docs page of a module or workflow: GET /api/docs/<id> returns mzdesk/docs/<id>.md."""
from fastapi.testclient import TestClient

from mzdesk.api.app import DOCS, app
from mzdesk.spec import DEFERRED_MODULES, DEFERRED_SUBGROUPS, MODULES, PIPELINES

client = TestClient(app)
# read at import time: conftest.py switches the held-back pages on for the other tests
HELD_BACK = {m["id"] for m in MODULES if m["subgroup"] in DEFERRED_SUBGROUPS} | set(DEFERRED_MODULES)


def test_every_visible_page_has_notes():
    ids = [m["id"] for m in MODULES if m["id"] not in HELD_BACK] + [p["id"] for p in PIPELINES]
    missing = [i for i in ids if not (DOCS / f"{i}.md").is_file()]
    assert not missing, f"no docs file for {missing}"


def test_notes_are_served_as_markdown():
    r = client.get("/api/docs/identity")
    assert r.status_code == 200
    body = r.json()
    assert body["id"] == "identity" and body["markdown"].startswith("# Identity")
    assert client.get("/api/docs/esi_ms2_annotation").status_code == 200       # workflows too


def test_unknown_and_deferred_pages_answer_404():
    assert client.get("/api/docs/nope").status_code == 404
    assert client.get("/api/docs/icp_isotopes").status_code == 404             # deferred module: not offered in this release
    assert "markdown" not in client.get("/api/docs/..%2Fspec").text            # no path tricks
