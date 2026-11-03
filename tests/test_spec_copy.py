"""Page copy (mzdesk/spec.py): every shown page says in one line what goes in and what comes out, names what it is based on, and has none of the
"online lookup" / "one-time download" remarks. The long text that left the pages is kept in mzdesk/docs/<id>.md, one file per shown page.
The text a run adds (warnings, notes, help) is short and plain as well: the last test runs the page examples.
(conftest.py switches the held-back pages on inside tests, so the shown pages are worked out once, when this file is collected.)"""
import json
from pathlib import Path

import pytest

from mzdesk import spec

DOCS = Path(spec.__file__).resolve().parent / "docs"
SHOWN_MODULES = [m for m in spec.MODULES if not spec.is_deferred(m)]
SHOWN = SHOWN_MODULES + list(spec.PIPELINES)
LOCAL_ONLY = set()          # Mono / FRIGID / Docker pages are installed on the hosted site too (deploy/install-engines.sh): no page is marked local-only
BANNED = ("online lookup", "(online)", "one-time", "one time", "first use")


def _ids(pages):
    return [p["id"] for p in pages]


def _strings(node, key=""):
    """Every string in a page definition except the network consent sentences and flags (they drive the confirmation dialog and may say "online")."""
    if key.startswith("network"):
        return
    if isinstance(node, str):
        yield node
    elif isinstance(node, dict):
        for k, v in node.items():
            yield from _strings(v, str(k))
    elif isinstance(node, (list, tuple)):
        for v in node:
            yield from _strings(v, key)


def test_the_shown_pages_leave_out_what_waits_for_2_0():
    assert _ids(spec.PIPELINES) == ["esi_ms2_annotation", "ei_ms_annotation", "targeted_method", "lipid_fragments"]
    assert not {"mrm", "met_shift", "stability", "nca", "icp_isotopes", "icp_sp"} & set(_ids(SHOWN_MODULES)) and len(SHOWN_MODULES) > 20


@pytest.mark.parametrize("page", SHOWN, ids=_ids(SHOWN))
def test_every_shown_page_names_what_it_is_based_on(page):
    based_on = page.get("based_on")
    assert isinstance(based_on, list) and 1 <= len(based_on) <= 4, page["id"]
    assert all(isinstance(x, str) and x.strip() == x and x for x in based_on) and len(set(based_on)) == len(based_on), based_on


@pytest.mark.parametrize("page", SHOWN, ids=_ids(SHOWN))
def test_every_one_liner_is_a_plain_have_to_get_line(page):
    line = page["one_liner"]
    have, arrow, get = line.partition("→")
    assert arrow and have.strip() and get.strip() and "→" not in get, line      # "<what I have> → <what I get>", one arrow
    assert "<" not in line and not line.endswith("."), line                       # plain text, no trailing period
    assert len(line) <= 90, line
    assert "Local version only" not in line, line      # every page runs on the hosted site as well


def test_no_page_text_has_an_online_lookup_or_one_time_download_remark():
    for page in [*spec.MODULES, *spec.PIPELINES]:
        for text in _strings(page):
            assert not any(phrase in text.lower() for phrase in BANNED), (page["id"], text)


@pytest.mark.parametrize("page", SHOWN, ids=_ids(SHOWN))
def test_every_shown_page_has_its_docs_file(page):
    path = DOCS / f"{page['id']}.md"
    assert path.is_file(), path
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    assert lines[0] == f"# {page['title']}" and lines[1] == "" and lines[2].strip(), path        # title, blank line, the original one-liner
    assert "<b>" not in text and "<i>" not in text and text.endswith("\n"), path
    assert {line for line in lines if line.startswith("## ")} <= {"## Inputs", "## Results"}, path


# ---------------------------------------------------------------- the text of a run
RUN_TEXT_MAX = 200
RUN_BANNED = ("online lookup", "one-time", "first use")
NEEDS_NETWORK = {"gnps_fastmasst"}              # the example sends the spectrum to GNPS2 (tests/test_api.py runs it live; this test does not need to)
DATA_KEYS = {"rows", "groups", "fragments", "formats", "files", "peaks", "plotly", "svg", "sets", "compounds"}       # data and quoted sources (JP18 printed conditions), not sentences
PROSE_KEYS = {"help", "note", "items", "caption"}


def _texts(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from _texts(item)


def _prose(node):
    """The sentences a reader gets from the outputs of a run: the items of the notes blocks and every help / note / caption."""
    if isinstance(node, dict):
        for key, value in node.items():
            if key in DATA_KEYS:
                continue
            if key in PROSE_KEYS:
                yield from _texts(value)
            elif isinstance(value, (dict, list)):
                yield from _prose(value)
    elif isinstance(node, list):
        for item in node:
            yield from _prose(item)


@pytest.mark.parametrize("module", SHOWN_MODULES, ids=_ids(SHOWN_MODULES))
def test_the_example_of_every_shown_page_runs_with_short_plain_text(module, monkeypatch):
    """Warnings, notes and help of a run stay short and carry no "online lookup" / "one-time" / "first use" remark. The example needs the engine the
    page needs: what is not installed here (Docker, Mono, FRIGID, NEIMS, MassBank, MolGpKa, Java) or needs the network is skipped, as in tests/test_api.py."""
    import test_api as api

    if module["id"] in NEEDS_NETWORK:
        pytest.skip("the example needs the network")
    skip = {"ms2_esi": api.NO_DOCKER, "ms2_structure": api.NO_MSFINDER, "ms2_denovo": None if api.RUN_FRIGID and not api.NO_FRIGID else "FRIGID is opt-in (MZDESK_RUN_FRIGID=1)",
            "ms2_ei": None if not api.NO_NEIMS else "NEIMS weights not installed", "ms2_library": None if not api.NO_MASSBANK else "MassBank index not installed",
            "separation": None if not api.NO_MOLGPKA else "MolGpKa weights not installed", "ms2_fragments": None if api.HAS_JAVA else "MetFrag needs Java"}.get(module["id"])
    if skip:
        pytest.skip(skip)
    if module["id"] == "identity":                      # the same offline stand-ins as tests/test_api.py
        from mzdesk.identity import enrich, resolve
        monkeypatch.setattr(resolve, "resolve", lambda q, kind="Auto", **kw: api._fake_resolved(query=q))
        monkeypatch.setattr(enrich, "enrich", lambda r: {"outputs": enrich.empty_outputs(), "warnings": []})
    if module["id"] == "exclusivity":
        from mzdesk.ms2 import msfinder_engine
        run = json.loads((Path(api.__file__).parent / "fixtures" / "isomers" / "caffeine_local.json").read_text())["explain_candidates"]
        monkeypatch.setattr(msfinder_engine, "explain_candidates", lambda *a, **kw: run)
    response = api.client.post(f"/api/modules/{module['id']}/run", json=module["example"])
    assert response.status_code == 200, response.text
    body = response.json()
    sentences = [*body["warnings"], *_prose(body["outputs"])]
    assert all(isinstance(s, str) for s in sentences)
    too_long = [s for s in sentences if len(s) > RUN_TEXT_MAX]
    assert not too_long, (module["id"], too_long)
    said = [s for s in [*sentences, body.get("verdict") or ""] if any(phrase in s.lower() for phrase in RUN_BANNED)]
    assert not said, (module["id"], said)
