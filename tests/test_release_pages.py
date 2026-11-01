"""v0.3.3: what the page list offers. DMPK, ICP-MS and MRM Transitions are built but left for 2.0; Report is gone; Calculator stands alone;
Alternative Explanations sits in Spectrum. (conftest.py switches the held-back pages on for every other test file, not this one.)"""
import pytest
from fastapi.testclient import TestClient

from mzdesk.api.app import app
from mzdesk.spec import DEFERRED_MODULES, DEFERRED_SUBGROUPS, MODULES, PIPELINES, spec

client = TestClient(app)


def test_page_list_leaves_out_what_waits_for_2_0():
    public = spec()
    ids = {m["id"] for m in public["modules"]}
    held_back = {m["id"] for m in MODULES if m["subgroup"] in DEFERRED_SUBGROUPS} | DEFERRED_MODULES
    assert DEFERRED_MODULES == {"mrm"} and DEFERRED_SUBGROUPS == {"dmpk", "inorganic"} and len(held_back) > 3
    assert not held_back & ids and "report" not in {m["id"] for m in MODULES}
    assert [s["title"] for s in public["subgroups"]] == ["Compound", "Detectability", "MS1", "Spectrum", "Orthogonal", "Dataset Search"]
    assert all(m["subgroup"] is None or m["subgroup"] in {s["id"] for s in public["subgroups"]} for m in public["modules"])


def test_calculator_stands_alone_and_alternative_explanations_sit_in_spectrum():
    mods = {m["id"]: m for m in spec()["modules"]}
    assert mods["calculator"]["subgroup"] is None and mods["exclusivity"]["subgroup"] == "ms2"
    assert [m["id"] for m in spec()["modules"] if m["subgroup"] == "ms2"][-1] == "exclusivity"


@pytest.mark.parametrize("module_id", ["mrm", "met_shift", "stability", "nca", "icp_isotopes", "icp_sp", "report"])
def test_held_back_modules_are_not_served(module_id):
    assert client.post(f"/api/modules/{module_id}/run", json={}).status_code == 404


def test_the_workflow_is_called_ms_ms_annotation():
    assert PIPELINES[0]["id"] == "esi_ms2_annotation" and PIPELINES[0]["title"] == "MS/MS Annotation"
