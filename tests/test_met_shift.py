"""Mass-shift search API: parent ion, feature parsing, shift matching and error handling."""
from fastapi.testclient import TestClient

from mzdesk import spec
from mzdesk.api.app import app
from mzdesk.calculator import mass_shift


client = TestClient(app)


def _module(module_id):
    return next(module for module in spec.MODULES if module["id"] == module_id)


def _run(inputs):
    return client.post("/api/modules/met_shift/run", json=inputs)


def test_every_spec_shift_option_is_known_to_the_engine():
    options = next(field for field in _module("met_shift")["inputs"] if field["id"] == "shifts")["options"]
    assert len(options) == len(set(options))
    result = mass_shift.search([], 200.0, 1, options, 10.0)
    assert [row["label"] for row in result["rows"]] == options
    glycine = next(row for row in result["rows"] if "glycine" in row["label"])
    assert abs(glycine["delta"] - 57.02146) < 1e-4
    response = _run({
        "structure": {"mode": "Name / ID / Formula", "value": "C8H10N4O2"},
        "precursor": {"precursor_mz": "", "ion_mode": "Positive", "adduct": "Auto"},
        "shifts": options,
        "tol_ppm": 10,
    })
    assert response.status_code == 200, response.text
    assert len(response.json()["outputs"]["hits"]["rows"]) == len(options)


def test_example_matches_three_shifts_and_three_features():
    response = _run(_module("met_shift")["example"])
    assert response.status_code == 200, response.text
    body = response.json()
    assert body["status"] == "ok"
    rows = body["outputs"]["hits"]["rows"]
    assert len(rows) == 4
    assert rows[0][0] == "+O (hydroxylation)"
    assert rows[0][3] == "168.06550"
    assert rows[1][0] == "−2H (dehydrogenation)"
    assert rows[1][3] == "—"
    assert body["outputs"]["features_plot"]["plotly"]["data"]
    assert "3/4" in body["verdict"]
    assert body["outputs"]["notes"]["items"]


def test_formula_parent_and_no_features():
    response = _run({
        "structure": {"mode": "Name / ID / Formula", "value": "C8H10N4O2"},
        "precursor": {"precursor_mz": "", "ion_mode": "Positive", "adduct": "Auto"},
        "shifts": ["+O (hydroxylation)"],
        "tol_ppm": 10,
    })
    assert response.status_code == 200, response.text
    body = response.json()
    rows = body["outputs"]["hits"]["rows"]
    assert len(rows) == 1
    assert rows[0][2] == "211.08257"
    assert rows[0][3] == "—"
    assert "load a feature table" in body["verdict"]


def test_two_column_features_warn_and_still_match():
    response = _run({
        "structure": {"mode": "SMILES / InChI", "value": spec.CAFFEINE},
        "features": "195.0877 100000\n211.0826 12000",
        "shifts": ["+O (hydroxylation)"],
    })
    assert response.status_code == 200, response.text
    body = response.json()
    assert any("RT is left blank" in warning for warning in body["warnings"])
    assert body["outputs"]["hits"]["rows"][0][3] == "211.08260"


def test_precursor_mismatch_warns():
    response = _run({
        "structure": {"mode": "SMILES / InChI", "value": spec.CAFFEINE},
        "precursor": {"precursor_mz": "196.0877", "ion_mode": "Positive", "adduct": "[M+H]+"},
        "features": "211.0826 6.0 1000",
        "shifts": ["+O (hydroxylation)"],
    })
    assert response.status_code == 200, response.text
    assert any("ppm from the theoretical" in warning for warning in response.json()["warnings"])


def test_invalid_inputs_return_400():
    base = {"structure": {"mode": "SMILES / InChI", "value": spec.CAFFEINE}}
    assert _run({**base, "shifts": []}).status_code == 400
    assert _run({**base, "features": "not a feature", "shifts": ["+O (hydroxylation)"]}).status_code == 400
    assert _run({**base, "features": "195.0877 6 100", "shifts": ["unknown shift"]}).status_code == 400
