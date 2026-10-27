"""Calculator engine and API: ppm, isotope labels, concentration, CE, Kovats, feature parsing."""
import pytest
from fastapi.testclient import TestClient

from mzdesk import spec
from mzdesk.api.app import app
from mzdesk.calculator import calculations, mass_shift


client = TestClient(app)


def _module(module_id):
    return next(module for module in spec.MODULES if module["id"] == module_id)


def test_ppm_error_sign_and_zero():
    assert calculations.ppm_error(195.0881, 195.08765) == pytest.approx(2.306655, rel=1e-5)
    assert calculations.ppm_error(100.0, 100.0) == 0.0
    with pytest.raises(ValueError):
        calculations.ppm_error(100.0, 0.0)


def test_concentration_units_and_dilution():
    assert calculations.concentration(19.419, "mg/L", "µM", 194.19) == pytest.approx(100.0)
    assert calculations.concentration(1.0, "µg/mL", "mg/L") == pytest.approx(1.0)
    assert calculations.concentration(1.0, "mM", "µM") == pytest.approx(1000.0)
    with pytest.raises(ValueError):
        calculations.concentration(1.0, "mg/L", "µM")
    mix = calculations.dilution(100.0, "µM", 10.0, "µM", 10.0)
    assert mix["stock_in_target"] == pytest.approx(100.0)
    assert mix["stock_volume"] == pytest.approx(1.0)
    assert mix["diluent_volume"] == pytest.approx(9.0)
    with pytest.raises(ValueError):
        calculations.dilution(10.0, "µM", 20.0, "µM", 10.0)
    with pytest.raises(ValueError):
        calculations.dilution(100.0, "µM", 10.0, "µM", 0.0)


def test_collision_energy_round_trip():
    ev = calculations.nce_to_ev(30.0, 195.0877, 1)
    assert ev == pytest.approx(11.705262, rel=1e-6)
    assert calculations.ev_to_nce(ev, 195.0877, 1) == pytest.approx(30.0)
    with pytest.raises(ValueError):
        calculations.nce_to_ev(30.0, 0.0, 1)


def test_kovats_interpolation_and_extrapolation():
    points = [(8, 8.0), (9, 9.0), (10, 10.0)]
    assert calculations.kovats_ri(points, 9.5)["ri"] == pytest.approx(950.0)
    assert calculations.kovats_ri(points, 9.0)["ri"] == pytest.approx(900.0)
    below = calculations.kovats_ri(points, 7.0)
    assert below["ri"] == pytest.approx(700.0)
    assert below["extrapolated"] is True
    with pytest.raises(ValueError):
        calculations.kovats_ri([(8, 8.0), (9, 7.0)], 8.5)
    assert calculations.parse_alkanes("C8 8\nC9 9\nC10 10") == [(8, 8.0), (9, 9.0), (10, 10.0)]
    assert calculations.parse_alkanes("Carbon number, RT\nC8, 8\nC9, 9") == [(8, 8.0), (9, 9.0)]
    with pytest.raises(ValueError):
        calculations.parse_alkanes("C8 8")


def test_isotope_labels_distribution_and_crosstalk():
    labels = calculations.parse_labels("13C3, d4, 15N2")
    assert [(label.element, label.count) for label in labels] == [("C", 3), ("H", 4), ("N", 2)]
    assert labels[0].nominal_shift == 3
    assert labels[0].exact_shift == pytest.approx(3.0100645, rel=1e-6)
    formula, mass = calculations.smiles_formula_mass(spec.CAFFEINE)
    assert formula == "C8H10N4O2"
    assert mass == pytest.approx(194.08037556, rel=1e-8)
    result = calculations.labelled_distribution(formula, "[M+H]+", calculations.parse_labels("13C3"))
    assert result["mono_mz"] == pytest.approx(198.0977165, rel=1e-8)
    assert result["peaks"][0]["abundance"] == pytest.approx(100.0)
    assert result["peaks"][1]["abundance"] == pytest.approx(7.118655, rel=1e-4)
    assert calculations.crosstalk(formula, "[M+H]+", 3) == pytest.approx(0.056366, rel=1e-3)
    with pytest.raises(ValueError):
        calculations.validate_labels(calculations.parse_labels("13C9"), formula)
    with pytest.raises(ValueError):
        calculations.parse_labels("13C0")
    with pytest.raises(ValueError):
        calculations.parse_labels("foo")


def test_adduct_mass_and_charge():
    assert calculations.adduct_charge("[M+H]+") == 1
    assert calculations.adduct_charge("[M-H]-") == -1
    assert calculations.adduct_mz("C8H10N4O2", "[M+H]+") == pytest.approx(195.087652, abs=1e-6)
    assert calculations.adduct_mz("C8H10N4O2", "[M-H]-") == pytest.approx(193.073099, abs=1e-6)


def test_feature_parser_header_two_column_and_errors():
    features, warnings = mass_shift.parse_features("m/z (Da)\tRT (min)\tIntensity\n195.0877\t6.0\t1000\n")
    assert features == [{"mz": 195.0877, "rt": 6.0, "intensity": 1000.0, "line": 2}]
    assert warnings == []
    features, warnings = mass_shift.parse_features("195.0877 1000\n211.0826 20")
    assert features[0]["rt"] is None
    assert warnings
    with pytest.raises(ValueError):
        mass_shift.parse_features("foo bar")
    with pytest.raises(ValueError):
        mass_shift.parse_features("m/z RT\n")


def test_mass_shift_search_matches_and_unmatched_rows():
    features = [
        {"mz": 195.0877, "rt": 6.0, "intensity": 100000.0, "line": 2},
        {"mz": 211.0826, "rt": 4.8, "intensity": 12000.0, "line": 3},
        {"mz": 181.0721, "rt": 5.5, "intensity": 8000.0, "line": 4},
    ]
    result = mass_shift.search(features, 195.087652, 1, ["+O (hydroxylation)", "−CH2 (demethylation)"], 10)
    assert result["matched_shifts"] == 2
    assert len(result["matched_features"]) == 2
    assert [row["label"] for row in result["rows"]] == ["+O (hydroxylation)", "−CH2 (demethylation)"]
    assert all(row["feature"] is not None for row in result["rows"])
    unmatched = mass_shift.search(features, 195.087652, 1, ["−2H (dehydrogenation)"], 10)
    assert unmatched["matched_shifts"] == 0
    assert unmatched["rows"][0]["feature"] is None
    assert mass_shift.search(features, 195.087652, 1, ["-2H (dehydrogenation)"], 10)["rows"][0]["label"] == "−2H (dehydrogenation)"
    with pytest.raises(ValueError):
        mass_shift.search(features, 195.087652, 1, [], 10)
    with pytest.raises(ValueError):
        mass_shift.search(features, 195.087652, 1, ["nope"], 10)


def test_calculator_api_example():
    response = client.post("/api/modules/calculator/run", json=_module("calculator")["example"])
    assert response.status_code == 200, response.text
    body = response.json()
    assert body["status"] == "ok"
    assert body["outputs"]["ppm_error"]["values"]["ppm error"]["value"] == "+2.31"
    assert body["outputs"]["isotope_labels"]["rows"][0][0] == "M (13C3)"
    assert body["outputs"]["crosstalk"]["rows"][0][2] == "0.05637"
    assert body["outputs"]["concentration"]["values"]["Converted"]["value"] == "100"
    assert body["outputs"]["dilution"]["rows"][0][3] == "1"
    assert body["outputs"]["collision_energy"]["values"]["NCE ↔ eV"]["value"] == "30 NCE · 11.71 eV"
    assert body["outputs"]["kovats"]["values"]["Kovats RI"]["value"] == "950.0"
    assert body["outputs"]["notes"]["items"]


def test_calculator_api_errors_and_missing_molar_mass():
    response = client.post("/api/modules/calculator/run", json={"lab_smiles": spec.CAFFEINE, "lab_label": "13C99"})
    assert response.status_code == 400
    response = client.post("/api/modules/calculator/run",
                           json={"conc_value": 10, "conc_from": "mg/L", "conc_to": "µM"})
    assert response.status_code == 200
    assert response.json()["warnings"]
