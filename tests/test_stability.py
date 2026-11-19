"""Metabolic stability: time-course parsing, hand-checked fit / CLint / scaling arithmetic, and the runner."""
import pytest
from fastapi.testclient import TestClient

from mzdesk import spec
from mzdesk.api.app import app
from mzdesk.dmpk import stability


client = TestClient(app)
EXAMPLE = "Time (min)\t% remaining\n0\t100\n5\t82\n15\t61\n30\t40\n45\t26"


def _run(module_id, inputs):
    return client.post(f"/api/modules/{module_id}/run", json=inputs)


def _module(module_id):
    return next(module for module in spec.MODULES if module["id"] == module_id)


# ---------------------------------------------------------------- parsing
def test_parse_accepts_header_comments_blank_lines_and_separators():
    text = "# HLM, 0.5 mg/mL\n\nTime (min), % remaining\n0, 100\n5, 82  # second point\n15, 61\n"
    assert stability.parse_timecourse(text) == [(0.0, 100.0), (5.0, 82.0), (15.0, 61.0)]
    assert stability.parse_timecourse("Time (min)\t% remaining\n0\t100\n5\t82\n15\t61") == [(0.0, 100.0), (5.0, 82.0), (15.0, 61.0)]
    assert stability.parse_timecourse("\ufefftime_min %_remaining\n0 100\n5 82\n15 61") == [(0.0, 100.0), (5.0, 82.0), (15.0, 61.0)]
    assert stability.parse_timecourse("0 100\n5 82\n15 61") == [(0.0, 100.0), (5.0, 82.0), (15.0, 61.0)]


@pytest.mark.parametrize("text,message", [
    ("0 100\n5 82", "At least three time points"),
    ("0 100\n5 82\n10 nope", "Line 3"),
    ("0 100\n5 82\n5 60", "Line 3"),
    ("0 100\n5 82\n10 101", "Line 3"),
    ("0 100\n5 82\n10 0", "Line 3"),
    ("0 100\n5 82 7", "Line 2"),
    ("0 100\n-5 82\n10 60", "Line 2"),
])
def test_parse_rejects_bad_rows_and_names_the_line(text, message):
    with pytest.raises(ValueError, match=message):
        stability.parse_timecourse(text)


# ---------------------------------------------------------------- fit and clearance arithmetic
def test_fit_and_half_life_are_hand_checked():
    # x = 0, 10, 20 min; % = 100, 70, 60.  y = ln(%/100) = 0, -0.3566749, -0.5108256
    # x̄ = 10, ȳ = -0.2891669; Sxx = 200, Sxy = -5.1082562 -> slope = -0.02554128 min⁻¹, k = 0.02554128
    # t½ = ln 2 / k = 27.1383 min; R² = Sxy²/(Sxx·Syy) = 0.950214
    result = stability.fit_curve([(0.0, 100.0), (10.0, 70.0), (20.0, 60.0)])
    assert result["k"] == pytest.approx(0.02554128, rel=1e-6)
    assert result["t_half"] == pytest.approx(27.13831, rel=1e-6)
    assert result["r2"] == pytest.approx(0.950214, abs=1e-6)
    assert result["n"] == 3


def test_a_slope_that_is_not_negative_is_rejected():
    with pytest.raises(ValueError, match="not negative"):
        stability.fit_curve([(0.0, 50.0), (10.0, 60.0), (20.0, 100.0)])


def test_clint_in_vitro_uses_the_amount_per_ml():
    # k × 1000 / amount: 0.02 min⁻¹ at 0.5 mg/mL -> 40 µL/min/mg; at 0.25 × 10⁶ cells/mL -> 80 µL/min/10⁶ cells
    assert stability.clint_invitro(0.02, "Liver microsomes", 0.5) == pytest.approx(40.0)
    assert stability.clint_invitro(0.02, "S9", 0.5) == pytest.approx(40.0)
    assert stability.clint_invitro(0.02, "Hepatocytes", 0.25) == pytest.approx(80.0)
    assert stability.clint_invitro(0.02, "Plasma", 0.5) is None
    assert stability.clint_invitro(0.02, "Liver microsomes", 0.0) is None
    assert stability.clint_unit("Liver microsomes") == "µL/min/mg protein"
    assert stability.clint_unit("Hepatocytes") == "µL/min/10⁶ cells"
    assert stability.clint_unit("Plasma") is None


def test_scaling_to_ml_min_kg_is_hand_checked():
    # microsomes: 40 µL/min/mg × 45 mg/g × 26 g/kg / 1000 µL/mL = 46.8 mL/min/kg (human)
    assert stability.scale_clint(40.0, "Liver microsomes", "Human") == pytest.approx(46.8)
    assert stability.scale_clint(40.0, "S9", "Human") == pytest.approx(46.8)          # S9 reuses the microsomal MPPGL
    # hepatocytes: 80 µL/min/10⁶ cells × 99 × 10⁶ cells/g × 40 g/kg / 1000 = 316.8 mL/min/kg (rat)
    assert stability.scale_clint(80.0, "Hepatocytes", "Rat") == pytest.approx(316.8)
    assert stability.scale_clint(None, "Liver microsomes", "Human") is None
    assert stability.scale_clint(40.0, "Plasma", "Human") is None


def test_well_stirred_clearance_is_hand_checked():
    # CLh = 20.7 × 46.8 / (20.7 + 46.8) = 14.352 mL/min/kg; ER = 14.352 / 20.7 = 0.69333 (human Qh)
    clh, er = stability.hepatic_clearance(46.8, "Human")
    assert clh == pytest.approx(14.352)
    assert er == pytest.approx(0.693333, rel=1e-6)
    assert stability.hepatic_clearance(None, "Human") == (None, None)
    flow_limited, ratio = stability.hepatic_clearance(1e9, "Human")
    assert flow_limited == pytest.approx(20.7) and ratio == pytest.approx(1.0)


def test_compute_covers_fit_clint_scaling_and_clearance():
    result = stability.compute(stability.parse_timecourse(EXAMPLE), "Liver microsomes", 0.5, "Human")
    assert result["k"] == pytest.approx(0.02945081, rel=1e-6)
    assert result["t_half"] == pytest.approx(23.53576, rel=1e-6)
    assert result["clint"] == pytest.approx(58.90163, rel=1e-6)
    assert result["clint_scaled"] == pytest.approx(68.91490, rel=1e-6)
    assert result["clh"] == pytest.approx(15.91854, rel=1e-6)
    assert result["er"] == pytest.approx(0.769012, rel=1e-6)
    assert result["liver_fraction"] == pytest.approx(0.026) and result["qh"] == pytest.approx(20.7)


def test_plasma_yields_no_clearance_and_a_non_positive_amount_is_not_an_error():
    points = stability.parse_timecourse(EXAMPLE)
    plasma = stability.compute(points, "Plasma", 0.5, "Human")
    assert (plasma["clint"], plasma["clint_scaled"], plasma["clh"], plasma["er"]) == (None, None, None, None)
    invalid = stability.compute(points, "Liver microsomes", 0.0, "Human")
    assert invalid["clint"] is None and invalid["clh"] is None


# ---------------------------------------------------------------- runner
def test_example_runs_and_fills_every_block():
    module = _module("stability")
    response = _run("stability", module["example"])
    assert response.status_code == 200, response.text
    body = response.json()
    assert body["status"] == "ok"
    expected = {block["id"] for block in module["outputs"] if block.get("id")}
    assert expected == {"summary", "fit_plot", "parameters", "notes"}
    assert expected <= set(body["outputs"])
    assert body["warnings"] == []


def test_runner_card_values_and_verdict_come_from_the_engine():
    body = _run("stability", _module("stability")["example"]).json()
    values = body["outputs"]["summary"]["values"]
    assert values["t½"]["value"] == "23.5"
    assert values["CLint"]["value"] == "58.9"
    assert values["Scaled CLint"]["value"] == "68.9"
    assert values["Predicted hepatic CL"]["value"] == "15.9"
    assert values["Extraction ratio"]["value"] == "0.769"
    assert body["verdict"] == "t½ <b>23.5 min</b> · CLint <b>58.9 µL/min/mg</b> · predicted hepatic CL <b>15.9 mL/min/kg</b>."
    rows = {row[0]: row for row in body["outputs"]["parameters"]["rows"]}
    assert list(rows) == ["k", "t½", "R²", "CLint (in vitro)", "MPPGL", "Liver weight fraction", "CLint scaled", "Qh", "fu",
                          "Predicted hepatic CL", "Extraction ratio"]
    assert rows["MPPGL"][1] == "45" and "Barter et al. 2007" in rows["MPPGL"][3]
    assert rows["Liver weight fraction"][1] == "0.026" and "Davies & Morris 1993" in rows["Liver weight fraction"][3]
    assert rows["fu"][1] == "1" and "assumed 1" in rows["fu"][3]

    traces = body["outputs"]["fit_plot"]["plotly"]["data"]
    assert [trace["mode"] for trace in traces] == ["lines", "markers"]
    # the figure plots ln(% remaining), so 100 % is ln 100 and the fit line has slope -k
    assert len(traces[1]["x"]) == 5 and traces[1]["y"][0] == pytest.approx(4.6051702)
    assert traces[0]["x"] == [0.0, 45.0]
    slope = (traces[0]["y"][1] - traces[0]["y"][0]) / 45.0
    assert slope == pytest.approx(-0.02945081, rel=1e-6)


def test_plasma_reports_half_life_only():
    body = _run("stability", {**_module("stability")["example"], "system": "Plasma"}).json()
    values = body["outputs"]["summary"]["values"]
    assert values["t½"]["value"] == "23.5"
    assert values["CLint"]["value"] == "—" and values["Predicted hepatic CL"]["value"] == "—"
    assert "plasma stability: t½ only" in values["CLint"]["note"]
    assert body["verdict"].endswith("plasma stability: t½ only.")
    rows = {row[0]: row for row in body["outputs"]["parameters"]["rows"]}
    assert rows["CLint (in vitro)"][1] == "—" and "plasma stability" in rows["CLint (in vitro)"][3]


def test_hepatocytes_and_s9_switch_the_scaling_factor_row():
    hepatocytes = _run("stability", {**_module("stability")["example"], "system": "Hepatocytes", "species": "Rat"}).json()
    rows = {row[0]: row for row in hepatocytes["outputs"]["parameters"]["rows"]}
    assert "HPGL" in rows and rows["HPGL"][1] == "99"
    assert rows["CLint (in vitro)"][2] == "µL/min/10⁶ cells"
    assert hepatocytes["outputs"]["summary"]["values"]["Scaled CLint"]["value"] == "233"
    s9 = _run("stability", {**_module("stability")["example"], "system": "S9"}).json()
    s9_rows = {row[0]: row for row in s9["outputs"]["parameters"]["rows"]}
    assert "MPPGL" in s9_rows and "simplification" in s9_rows["MPPGL"][3]


@pytest.mark.parametrize("timecourse,warning", [
    ("0 100\n10 40\n20 60", "Poor ln-linear fit"),
    ("0 100\n30 95\n60 91", "The last point is 91% remaining"),
    ("0 100\n10 85\n20 80", "longer than twice the incubation"),
])
def test_runner_warns_about_unreliable_estimates(timecourse, warning):
    body = _run("stability", {"timecourse": timecourse, "system": "Liver microsomes", "protein": 0.5, "species": "Human"}).json()
    assert any(warning in message for message in body["warnings"]), body["warnings"]


def test_runner_warns_and_reports_no_clearance_for_a_non_positive_amount():
    for protein in (0, -1):
        body = _run("stability", {**_module("stability")["example"], "protein": protein}).json()
        assert any("must be positive" in message for message in body["warnings"])
        values = body["outputs"]["summary"]["values"]
        assert values["CLint"]["value"] == "—" and values["Extraction ratio"]["value"] == "—"


@pytest.mark.parametrize("inputs,message", [
    ({"timecourse": ""}, "paste at least three"),
    ({"timecourse": "0 100\n5 82"}, "At least three"),
    ({"timecourse": "0 100\n5 82\n10 x"}, "Line 3"),
    ({"timecourse": "0 100\n5 82\n5 60"}, "must increase"),
    ({"timecourse": "0 50\n10 60\n20 100"}, "not negative"),
    ({"timecourse": EXAMPLE, "system": "Kidney"}, "Unknown system"),
    ({"timecourse": EXAMPLE, "species": "Zebrafish"}, "Unknown species"),
])
def test_bad_input_is_a_400_with_the_engine_message(inputs, message):
    response = _run("stability", inputs)
    assert response.status_code == 400
    assert message in response.json()["detail"]
