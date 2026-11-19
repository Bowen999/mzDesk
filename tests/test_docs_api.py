"""The Docs: GET /api/docs/<id> returns mzdesk/docs/<id>.md (general pages, modules, workflows) and GET /api/docs-img/<name> the figures in mzdesk/docs/img."""
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from mzdesk.api import app as app_module
from mzdesk.api.app import DOCS, GENERAL_PAGES, app
from mzdesk.spec import DEFERRED_MODULES, DEFERRED_SUBGROUPS, MODULES, PIPELINES

client = TestClient(app)
# read at import time: conftest.py switches the held-back pages on for the other tests
HELD_BACK = {m["id"] for m in MODULES if m["subgroup"] in DEFERRED_SUBGROUPS} | set(DEFERRED_MODULES)
DOCS_CSS = Path(__file__).resolve().parents[1] / "web" / "src" / "docs" / "docs.css"


@pytest.fixture
def docs_dir(tmp_path, monkeypatch):
    """A throw-away docs folder in place of mzdesk/docs, so a test does not depend on which pages are written."""
    (tmp_path / "img").mkdir()
    monkeypatch.setattr(app_module, "DOCS", tmp_path)
    return tmp_path


def test_every_visible_page_has_notes():
    ids = [m["id"] for m in MODULES if m["id"] not in HELD_BACK] + [p["id"] for p in PIPELINES]
    missing = [i for i in ids if not (DOCS / f"{i}.md").is_file()]
    assert not missing, f"no docs file for {missing}"


def test_every_general_page_has_a_file():
    assert len(GENERAL_PAGES) == len(set(GENERAL_PAGES)) == 14
    missing = [i for i in GENERAL_PAGES if not (DOCS / f"{i}.md").is_file()]
    assert not missing, f"no docs file for the general pages {missing}"


def test_notes_are_served_as_markdown():
    r = client.get("/api/docs/identity")
    assert r.status_code == 200
    body = r.json()
    assert body["id"] == "identity" and body["markdown"].startswith("# Identity")
    assert client.get("/api/docs/esi_ms2_annotation").status_code == 200       # workflows too


def test_general_pages_are_served_once_their_file_exists(docs_dir):
    assert client.get("/api/docs/quick-start").status_code == 404              # on the list, no file yet
    (docs_dir / "quick-start.md").write_text("# Quick Start\n\nRun `mzdesk`.\n", encoding="utf-8")
    r = client.get("/api/docs/quick-start")
    assert r.status_code == 200
    assert r.json() == {"id": "quick-start", "markdown": "# Quick Start\n\nRun `mzdesk`.\n"}


def test_only_listed_pages_are_served(docs_dir):
    (docs_dir / "secret.md").write_text("# Not a page", encoding="utf-8")
    (docs_dir / "overview.md").write_text("# Overview", encoding="utf-8")
    assert client.get("/api/docs/secret").status_code == 404                   # a file that is not on any list
    assert client.get("/api/docs/overview").status_code == 200


def test_unknown_and_deferred_pages_answer_404():
    assert client.get("/api/docs/nope").status_code == 404
    assert client.get("/api/docs/icp_isotopes").status_code == 404             # deferred module: not offered in this release
    assert "markdown" not in client.get("/api/docs/..%2Fspec").text            # no path tricks
    assert "markdown" not in client.get("/api/docs/%2e%2e%2fspec").text


def test_figures_are_served_as_svg(docs_dir):
    (docs_dir / "img" / "flow-1.svg").write_text('<svg viewBox="0 0 720 10"><title>t</title></svg>', encoding="utf-8")
    r = client.get("/api/docs-img/flow-1")
    assert r.status_code == 200
    assert r.headers["content-type"].startswith("image/svg+xml")
    assert r.text.startswith("<svg")
    csp = r.headers["content-security-policy"]
    assert "default-src 'none'" in csp and "sandbox" in csp                    # opened on its own, a figure cannot run anything


def test_figure_names_are_checked(docs_dir):
    for path in [docs_dir / "img" / "ok.svg", docs_dir / "secret.svg", docs_dir / "img" / "Upper.svg", docs_dir / "img" / "a_b.svg"]:
        path.write_text("<svg>SECRET</svg>", encoding="utf-8")
    assert client.get("/api/docs-img/ok").status_code == 200
    for name in ["nope", "Upper", "a_b", "ok.svg", "ok%00", "a%20b"]:
        assert client.get(f"/api/docs-img/{name}").status_code == 404, name
    for name in ["..%2Fsecret", "%2e%2e%2fsecret", "..", "img%2Fok"]:          # not even routed here: whatever answers, it is not a figure
        r = client.get(f"/api/docs-img/{name}")
        assert "SECRET" not in r.text and not r.headers["content-type"].startswith("image/svg"), name


def figures():
    return sorted((DOCS / "img").glob("*.svg"))


@pytest.mark.parametrize("path", figures(), ids=lambda p: p.stem)
def test_published_figures_follow_the_figure_rules(path):
    """Hand-written, theme-driven SVG: 720 wide, a title and a description, classes only from docs.css, no scripts, styles, colours or fonts of its own."""
    src = path.read_text(encoding="utf-8")
    assert src.lstrip().startswith("<svg") and re.fullmatch(r"[a-z0-9-]+", path.stem)
    root = ET.fromstring(src)                                                  # well-formed XML
    ns = "{http://www.w3.org/2000/svg}"
    assert root.tag == ns + "svg" and re.fullmatch(r"0 0 720 \d+", root.get("viewBox", "")), "viewBox must be 0 0 720 H"
    assert root.find(ns + "title") is not None and (root.find(ns + "title").text or "").strip()
    assert root.find(ns + "desc") is not None and (root.find(ns + "desc").text or "").strip()
    known = set(re.findall(r"\.(fg-[a-z0-9-]+)", DOCS_CSS.read_text(encoding="utf-8")))
    for el in root.iter():
        tag = el.tag.removeprefix(ns)
        assert tag not in {"script", "style", "foreignObject", "image", "a", "animate", "set"}, f"<{tag}> is not allowed"
        for name, value in el.attrib.items():
            assert not name.startswith("on"), f"{name} on <{tag}>"
            assert name not in {"style", "font-family"}, f"{name} on <{tag}>"
            assert not re.search(r"#[0-9a-fA-F]{3,8}\b|rgb|hsl", value) or name in {"d", "points"}, f"colour in {name}={value!r}"
            if name == "class":
                assert set(value.split()) <= known, f"unknown class in {value!r}"
            if name in {"fill", "stroke"}:
                assert value in {"none", "currentColor"} or value.startswith("url(#"), f"{name}={value!r}: colours come from the fg-* classes"
    assert "http://" not in src.replace("http://www.w3.org/2000/svg", "") and "https://" not in src
