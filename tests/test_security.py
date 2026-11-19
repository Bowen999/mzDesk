"""Security headers, request size limits, hidden API docs, error text without server folders, odd page addresses."""
import asyncio

import pytest
from fastapi.testclient import TestClient

from mzdesk.api import request_limits, security
from mzdesk.api.app import app

client = TestClient(app)


def test_every_answer_carries_the_security_headers():
    for r in (client.get("/api/health"), client.get("/api/spec"), client.post("/api/track", json={"path": "/x"}), client.get("/api/docs/nope")):
        h = r.headers
        assert h["x-content-type-options"] == "nosniff" and h["x-frame-options"] == "SAMEORIGIN"
        assert h["referrer-policy"] == "strict-origin-when-cross-origin" and h["cross-origin-opener-policy"] == "same-origin"
        assert "camera=()" in h["permissions-policy"]
        csp = h["content-security-policy"]
        assert "object-src 'none'" in csp and "frame-ancestors 'self'" in csp and "base-uri 'self'" in csp and "connect-src 'self'" in csp
        assert "script-src 'self' 'sha256-" in csp and "'unsafe-inline'" not in csp.split("style-src")[0]       # eval only for the structure editor; no inline script without a hash
        assert h["cache-control"] == "no-store"                                  # API answers are not cached


def test_the_inline_scripts_of_the_page_are_allowed_by_hash_not_by_unsafe_inline():
    html = "<html><script>window.a = 1</script><script src='/x.js'></script><script type='module'>import './m.js'</script></html>"
    hashes = security.inline_script_hashes(html)
    assert len(hashes) == 2 and all(h.startswith("'sha256-") for h in hashes)
    csp = security.content_security_policy(html)
    script_src = next(part for part in csp.split("; ") if part.startswith("script-src"))
    assert all(h in script_src for h in hashes) and "'unsafe-inline'" not in script_src


def test_the_depiction_svg_cannot_run_script():
    r = client.get("/api/depict", params={"smiles": "CCO"})
    assert r.status_code == 200 and r.headers["content-type"].startswith("image/svg+xml")
    assert r.headers["content-security-policy"].startswith("default-src 'none'")


def test_the_api_browser_and_schema_are_off_unless_asked_for():
    for path in ("/redoc", "/openapi.json"):                  # (/docs/<id> is a page of the app: the Docs pages)
        r = client.get(path)
        assert '"paths"' not in r.text and "swagger" not in r.text.lower()
    assert app.openapi_url is None and app.docs_url is None and app.redoc_url is None


def test_a_body_over_the_limit_is_refused_with_the_headers_still_on(monkeypatch):
    monkeypatch.setattr(request_limits, "DEFAULT_REQUEST_LIMIT", 100)
    r = client.post("/api/track", content=b'{"path": "' + b"a" * 200 + b'"}', headers={"content-type": "application/json"})
    assert r.status_code == 413 and r.json() == {"detail": "The request is too large."}
    assert r.headers["x-content-type-options"] == "nosniff"
    assert client.post("/api/track", json={"path": "/ok"}).status_code == 200
    assert client.get("/api/health").status_code == 200                      # GET bodies are not counted


def test_the_default_limit_applies_without_content_length():
    messages = iter([{"type": "http.request", "body": b"x" * 5_000_000, "more_body": True},
                     {"type": "http.request", "body": b"x" * 5_000_000, "more_body": False}])
    sent = []

    async def receive():
        return next(messages)

    async def send(message):
        sent.append(message)

    async def should_not_run(*args):
        raise AssertionError("oversized request reached the app")

    scope = {"type": "http", "method": "POST", "path": "/api/modules/properties/run", "headers": []}
    asyncio.run(request_limits.RequestLimit(should_not_run)(scope, receive, send))
    assert sent[0]["status"] == 413


def test_massql_keeps_its_larger_allowance():
    assert request_limits.limit_for("/api/modules/gnps_massql/run") == 27_000_000 > request_limits.limit_for("/api/modules/identity/run") == 8_000_000


def test_error_text_loses_the_folders_of_the_server():
    assert security.scrub("No such file: '/var/lib/mzdesk/models/neims/w.h5'") == "No such file: 'w.h5'"
    assert security.scrub("File \"/opt/mzdesk/.venv/lib/python3.12/site-packages/rdkit/x.py\", line 3") == "File \"x.py\", line 3"
    assert security.scrub("C:\\Users\\me\\data\\x.mgf not found") == "x.mgf not found"
    assert security.scrub("a\n\n  b") == "a b" and len(security.scrub("x" * 1000)) == 240
    assert security.strip_paths("InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H") == "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"        # chemistry with slashes is untouched
    assert security.strip_paths("m/z 100/200 at /dev/null") == "m/z 100/200 at /dev/null"


@pytest.mark.parametrize("path", ["/%00", "/..%2F..%2Fetc%2Fpasswd", "/%2e%2e/%2e%2e/etc/passwd", "/assets/..%2f..%2findex.html", "/" + "a" * 5000])
def test_odd_page_addresses_never_leak_files_or_crash(path):
    r = client.get(path)
    assert r.status_code in (200, 404, 503) and "root:" not in r.text
