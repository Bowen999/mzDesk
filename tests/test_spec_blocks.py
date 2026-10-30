"""Output-block helpers in mzdesk/spec.py and the data builders in mzdesk/api/blocks.py: shapes and backward compatibility."""
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest

from mzdesk import spec as S
from mzdesk.api import blocks as B
from mzdesk.api.common import clean

KNOWN_TYPES = {"table", "cards", "figure", "compare", "code", "formats", "notes", "button", "heading", "portal", "chips", "reaction", "grouptable", "downloads"}


# ---------------------------------------------------------------- backward compatibility
def test_old_calls_serialise_exactly_as_before():
    assert S.table(["a", "b"], "T", "cap", "k", "i") == {"type": "table", "title": "T", "columns": ["a", "b"], "caption": "cap", "key": "k", "id": "i"}
    assert S.cards([("L", "u")], "T", "k", "i") == {"type": "cards", "title": "T", "items": [{"label": "L", "unit": "u"}], "key": "k", "id": "i"}
    assert S.figure("F", "k", "i") == {"type": "figure", "label": "F", "key": "k", "id": "i"}
    assert S.code("C", "T", "k", "i") == {"type": "code", "title": "T", "label": "C", "key": "k", "id": "i"}
    assert S.formats("F", "T", "k", "i") == {"type": "formats", "title": "T", "label": "F", "key": "k", "id": "i"}
    assert S.notes("N", "k", "i") == {"type": "notes", "title": "N", "key": "k", "id": "i"}
    assert S.head("H", "k") == {"type": "heading", "title": "H", "key": "k"}
    assert S.button("B", "k") == {"type": "button", "label": "B", "key": "k"}


def test_optional_keys_appear_only_when_set():
    for block in (S.table(["a"]), S.cards([("L", "")]), S.notes("N"), S.figure("F"), S.portal(), S.chips(), S.reaction(), S.grouptable(["a"])):
        assert not {"help", "style", "row_help", "collapsed"} & set(block), block["type"]


# ---------------------------------------------------------------- spec helpers
def test_help_on_every_block_type():
    helps = [S.table(["a"], help="h"), S.cards([("L", "")], help="h"), S.figure("F", help="h"), S.code("C", help="h"), S.formats("F", help="h"),
             S.notes("N", help="h"), S.head("H", help="h"), S.portal(help="h"), S.chips(help="h"), S.reaction(help="h"), S.grouptable(["a"], help="h")]
    assert all(b["help"] == "h" for b in helps)
    assert S.table(["a"], help=["x", "y"])["help"] == ["x", "y"]


def test_table_row_help_and_cards_style_and_item_help():
    assert S.table(["a"], row_help=True)["row_help"] is True
    assert S.table(["a"], page_size=10, compact=True)["page_size"] == 10
    assert S.table(["a"], page_size=10, compact=True)["compact"] is True
    assert not {"page_size", "compact"} & set(S.table(["a"]))
    c = S.cards([("A", "u", "why"), ("B", "")], style="metric")
    assert c["style"] == "metric" and c["items"] == [{"label": "A", "unit": "u", "help": "why"}, {"label": "B", "unit": ""}]


def test_notes_collapsed():
    assert S.notes("Rule hits", collapsed=True)["collapsed"] is True
    assert S.notes(None, id="x", collapsed=True)["title"] is None


@pytest.mark.parametrize("helper, extra", [(S.portal, {}), (S.chips, {}), (S.reaction, {})])
def test_portal_chips_reaction_signature(helper, extra):
    # documented positional order: (title, id, help)
    b = helper("Title", "bid", "hlp", key="kk")
    assert b == {"type": helper.__name__, "title": "Title", "key": "kk", "id": "bid", "help": "hlp"}


def test_grouptable_signature():
    b = S.grouptable(["A", "B"], "Title", "gid", "hlp", row_help=True)
    assert b == {"type": "grouptable", "title": "Title", "columns": ["A", "B"], "caption": None, "key": None, "id": "gid", "help": "hlp", "row_help": True}


def test_every_module_block_is_well_formed():
    """Catches typos in other pages' blocks: known type, str / list[str] help, table-like blocks have columns."""
    for m in S.MODULES:
        for b in m["outputs"]:
            where = (m["id"], b.get("id"), b["type"])
            assert b["type"] in KNOWN_TYPES, where
            h = b.get("help")
            assert h is None or isinstance(h, str) or (isinstance(h, list) and all(isinstance(x, str) for x in h)), where
            if b["type"] in ("table", "grouptable"):
                assert b["columns"], where
            if b["type"] == "cards":
                assert all("label" in it for it in b["items"]), where


def test_restyled_pages_have_unique_icons_and_metric_cards():
    groups = {"molecular", "ms1", "ms2", "gnps", "tools"}
    modules = [m for m in S.MODULES if m["subgroup"] in groups]
    path = Path(__file__).resolve().parents[1] / "web/src/assets/icons.svg"
    ids = [node.get("id") for node in ET.parse(path).getroot()]
    assert len(ids) == len(set(ids))
    assert {m["id"] for m in modules} | {"esi_ms2_annotation"} <= set(ids)
    for m in modules:
        assert all(b.get("style") == "metric" for b in m["outputs"] if b["type"] == "cards"), m["id"]


# ---------------------------------------------------------------- data builders
def test_cells_and_rows():
    assert B.link("JP18", "https://x.org") == {"text": "JP18", "url": "https://x.org"}
    assert B.badge("JP", "jp") == {"badge": "JP", "tone": "jp"}
    assert B.badge("x") == {"badge": "x", "tone": "neutral"}
    assert B.hinted("ODS", "column") == {"text": "ODS", "help": "column"}
    assert B.row(["a", "b"]) == ["a", "b"]
    assert B.row(["a"], help=["why"]) == {"cells": ["a"], "help": ["why"]}
    assert B.table_data([["a"]], count=0) == {"rows": [["a"]], "count": 0}
    assert B.table_data([], help="h") == {"rows": [], "help": "h"}


def test_grouptable_data():
    g = B.group("Caffeine", [["a"]], meta="0.66 · 3 methods", help="h")
    assert g == {"title": "Caffeine", "meta": "0.66 · 3 methods", "help": "h", "rows": [["a"]]}
    assert B.grouptable_data([g, B.group("T", [])], count=2) == {"groups": [g, {"title": "T", "rows": []}], "count": 2}


def test_portal_and_chips_data():
    t = B.tile("GC-MS", "yes", "ok", chips=["TMS"], help="h")
    assert t == {"label": "GC-MS", "status": "yes", "value": "ok", "chips": ["TMS"], "help": "h"}
    assert B.tile("X", "no") == {"label": "X", "status": "no"}
    assert B.portal_data([t]) == {"tiles": [t]}
    c = B.chip("UV", "good", "273 nm", help=["a", "b"])
    assert c == {"label": "UV", "level": "good", "value": "273 nm", "help": ["a", "b"]}
    assert B.chips_data([c, B.chip("RI", "weak")]) == {"chips": [c, {"label": "RI", "level": "weak"}]}


@pytest.mark.parametrize("bad", ["maybe", "", "Yes"])
def test_tile_status_vocabulary(bad):
    with pytest.raises(ValueError):
        B.tile("X", bad)


@pytest.mark.parametrize("bad", ["great", "partial ", None])
def test_chip_level_vocabulary(bad):
    with pytest.raises(ValueError):
        B.chip("X", bad)


def test_reaction_data():
    s = B.step("R–NH₂", "Dns-Cl", "R–NH–Dns", 2)
    assert s == {"from": "R–NH₂", "reagent": "Dns-Cl", "to": "R–NH–Dns", "count": 2}
    assert B.step("a", None, "b") == {"from": "a", "to": "b"}          # count 1 and no reagent are left out
    o = B.reaction_option("CIL · dansyl", [s], formula="C₁₂H₁₂N₂O₂S", mass=248.0619, heavy=250.07, help="h")
    assert o == {"label": "CIL · dansyl", "steps": [s], "formula": "C₁₂H₁₂N₂O₂S", "mass": 248.0619, "heavy": 250.07, "help": "h"}
    assert B.reaction_data([o], svg="<svg/>") == {"svg": "<svg/>", "options": [o]}
    assert B.reaction_data([o]) == {"options": [o]}
    assert B.reaction_option("with SVG", [s], svg="<svg id='sites'/>")["svg"] == "<svg id='sites'/>"
    assert B.reaction_option("with product", [s], product_svg="<svg id='product'/>")["product_svg"] == "<svg id='product'/>"
    with pytest.raises(ValueError):
        B.reaction_option("x", [])
    with pytest.raises(ValueError):
        B.reaction_data([])


def test_builder_output_survives_the_api_cleaner_and_json():
    data = B.grouptable_data([B.group("G", [B.row([B.link("a", "https://x.org"), B.badge("JP", "jp")], help="h")], meta="m")], count=1)
    assert json.loads(json.dumps(clean(data))) == data
