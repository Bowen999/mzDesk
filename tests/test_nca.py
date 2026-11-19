"""NCA engine and API: profile parsing, every formula by hand, route conventions and user-facing errors."""
import json
import math

import pytest
from fastapi.testclient import TestClient

from mzdesk import spec
from mzdesk.api.app import app
from mzdesk.dmpk import nca


client = TestClient(app)
# C(t) = 100 * 2^(-t): every hand calculation below is exact for this profile.
PERFECT = [(0.0, 100.0), (1.0, 50.0), (2.0, 25.0), (4.0, 6.25)]


def _module():
    return next(module for module in spec.MODULES if module["id"] == "nca")


def _run(inputs):
    return client.post("/api/modules/nca/run", json=inputs)


def test_parse_profile_reads_header_comments_and_separators():
    text = ("\ufeff# oral profile\nTime (h) Concentration (ng/mL)\n\n0.5, 42\n1\t88\n2 131 # second sample\n4, 152\n")
    assert nca.parse_profile(text) == [(0.5, 42.0), (1.0, 88.0), (2.0, 131.0), (4.0, 152.0)]
    assert nca.parse_profile("0 100\n1 50\n2 25") == [(0.0, 100.0), (1.0, 50.0), (2.0, 25.0)]
    assert nca.parse_profile("Time (h),Concentration (ng/mL)\n0,100\n1,80\n2,40") == [(0.0, 100.0), (1.0, 80.0), (2.0, 40.0)]


@pytest.mark.parametrize("text, message", [
    ("", "empty"),
    ("Time (h) Concentration (ng/mL)", "At least 3"),
    ("0 100\n1 50", "At least 3"),
    ("0 100\n1 abc\n2 25", r"Row 2.*not a number"),
    ("0 100\n1 50 7\n2 25", r"Row 2.*expected 2 values"),
    ("0 100\n2 50\n1 25", r"Row 3.*strictly increase"),
    ("0 100\n1 0\n2 25", r"Row 2.*> 0"),
    ("0 nan\n1 50\n2 25", r"Row 1.*finite"),
])
def test_parse_profile_errors_name_the_row(text, message):
    with pytest.raises(ValueError, match=message):
        nca.parse_profile(text)


def test_linear_trapezoid_matches_hand_calculation():
    # 0.5*(100+50) + 1*(50+25) + 2*(25+6.25) = 75 + 37.5 + 31.25
    assert nca.linear_trapezoid([0, 1, 2, 4], [100, 50, 25, 6.25]) == pytest.approx(143.75)
    assert nca.linear_trapezoid([0, 1, 2], [0, 100, 100]) == pytest.approx(150.0)
    with pytest.raises(ValueError):
        nca.linear_trapezoid([0], [1])
    with pytest.raises(ValueError):
        nca.linear_trapezoid([0, 1], [1])


def test_lambda_z_of_a_perfect_exponential():
    fit = nca.fit_lambda_z([0, 1, 2, 4], [100, 50, 25, 6.25], 3)
    assert fit["lambda_z"] == pytest.approx(math.log(2))
    assert fit["slope"] == pytest.approx(-math.log(2))
    assert fit["intercept"] == pytest.approx(math.log(100))
    assert fit["r2"] == pytest.approx(1.0)
    assert fit["n_points"] == 3
    assert fit["times"] == [1.0, 2.0, 4.0]
    assert fit["predicted"] == pytest.approx([50.0, 25.0, 6.25])


def test_lambda_z_r2_matches_hand_regression():
    # ln C = 2.302585, 2.079442, 1.386294 at t = 0, 1, 2; slope -0.458145..., R² 0.919369...
    fit = nca.fit_lambda_z([0, 1, 2], [10, 8, 4], 3)
    assert fit["lambda_z"] == pytest.approx(0.45814536593707766)
    assert fit["intercept"] == pytest.approx(2.3809190312016684)
    assert fit["r2"] == pytest.approx(0.9193685919816033)
    assert fit["predicted"] == pytest.approx([10.81483747120199, 6.839903786706788, 4.325934988480795])


def test_lambda_z_rejects_a_non_declining_or_invalid_window():
    with pytest.raises(ValueError, match="does not decline"):
        nca.fit_lambda_z([0, 1, 2], [10, 20, 30], 3)
    with pytest.raises(ValueError, match="at least 3"):
        nca.fit_lambda_z([0, 1, 2], [10, 5, 2], 2)
    with pytest.raises(ValueError, match="cannot exceed"):
        nca.fit_lambda_z([0, 1, 2], [10, 5, 2], 4)
    with pytest.raises(ValueError, match="integer"):
        nca.fit_lambda_z([0, 1, 2], [10, 5, 2], 3.0)
    with pytest.raises(ValueError, match="> 0"):
        nca.fit_lambda_z([0, 1, 2], [10, 0, 2], 3)
    with pytest.raises(ValueError, match="strictly increase"):
        nca.fit_lambda_z([0, 2, 1], [10, 5, 2], 3)


def test_analyse_perfect_exponential_every_formula():
    result = nca.analyse(PERFECT, 10.0, "Oral", 3)
    assert result["cmax"] == 100.0 and result["tmax"] == 0.0
    assert result["auc_t"] == pytest.approx(143.75)
    assert result["aumc_t"] == pytest.approx(150.0)             # 25 + 50 + 75
    assert result["lambda_z"] == pytest.approx(math.log(2))
    assert result["t_half"] == pytest.approx(1.0)
    # AUC0-∞ = 143.75 + 6.25/ln2; extrapolated = 9.0168/152.7668 = 5.902 %
    assert result["auc_inf"] == pytest.approx(152.766844005556)
    assert result["auc_extrapolated_pct"] == pytest.approx(5.902356669244333)
    # AUMC0-∞ = 150 + 4*6.25/ln2 + 6.25/ln2²
    assert result["aumc_inf"] == pytest.approx(199.07593215350914)
    assert result["mrt"] == pytest.approx(1.3031357258795562)
    # CL/F = 10 mg * 1e6 / 152.7668 ng·h/mL / 1000 = 65.459 L/h; Vz/F = CL/λz
    assert result["cl"] == pytest.approx(65.45923014313438)
    assert result["vz"] == pytest.approx(94.43770670790931)
    assert result["r2"] == pytest.approx(1.0)
    assert result["n_points"] == 4 and result["n_lambda"] == 3


def test_cmax_takes_the_first_observation_on_ties():
    result = nca.analyse([(0, 10), (1, 30), (2, 30), (4, 5)], 1.0, "IV bolus", 3)
    assert result["cmax"] == 30.0 and result["tmax"] == 1.0


def test_analyse_validates_dose_route_and_points():
    with pytest.raises(ValueError, match="positive"):
        nca.analyse(PERFECT, 0, "Oral", 3)
    with pytest.raises(ValueError, match="route"):
        nca.analyse(PERFECT, 10, "Subcutaneous", 3)
    with pytest.raises(ValueError, match="At least 3"):
        nca.analyse(PERFECT[:2], 10, "Oral", 3)
    with pytest.raises(ValueError, match="cannot exceed"):
        nca.analyse(PERFECT, 10, "Oral", 5)
    with pytest.raises(ValueError, match="strictly increase"):
        nca.analyse([(0, 10), (2, 5), (1, 2)], 10, "Oral", 3)


# ---------------------------------------------------------------- API
def test_example_runs_and_fills_every_block():
    module = _module()
    response = _run(module["example"])
    assert response.status_code == 200, response.text
    body = response.json()
    assert body["status"] == "ok"
    block_ids = {block["id"] for block in module["outputs"] if block.get("id")}
    assert block_ids == {"summary", "profile_plot", "parameters", "notes"}
    assert block_ids <= set(body["outputs"])

    summary = body["outputs"]["summary"]["values"]
    assert summary["Cmax"] == "152"
    assert summary["Tmax"] == "4"
    assert summary["AUC0-t"] == "1,450"
    assert summary["AUC0-∞"] == "1,520"
    assert summary["t½"] == "4.88"
    assert summary["CL/F"]["value"] == "65.9"
    assert summary["Vz/F"]["value"] == "464"
    assert summary["MRT"]["value"] == "8.14"
    assert body["warnings"] == []

    for needle in ("4.88 h", "1,520 ng·h/mL", "65.9 L/h"):
        assert needle in body["verdict"]


def test_example_parameters_table_and_figure():
    body = _run(_module()["example"]).json()
    rows = body["outputs"]["parameters"]["rows"]
    assert [row[0] for row in rows] == [
        "Cmax", "Tmax", "AUC0-t", "AUC0-∞", "AUC extrapolated (%)", "t½", "λz", "CL/F", "Vz/F", "MRT",
        "Points used", "R²",
    ]
    by_name = {row[0]: row for row in rows}
    assert by_name["t½"][1] == "4.88" and by_name["t½"][2] == "h"
    assert by_name["λz"][1] == "0.1421" and by_name["λz"][2] == "1/h"
    assert by_name["AUC extrapolated (%)"][1] == "4.18"
    assert by_name["Points used"][1] == "4"
    assert by_name["R²"][1] == "0.9924"
    assert "OLS" in by_name["λz"][3] and "last 4 points" in by_name["λz"][3]

    plot = body["outputs"]["profile_plot"]["plotly"]
    assert plot["layout"]["yaxis"]["type"] == "log"
    assert len(plot["data"]) == 2
    observed, fit = plot["data"]
    assert observed["mode"] == "lines+markers"
    assert observed["y"] == [42, 88, 131, 152, 121, 84, 41, 9]
    assert fit["line"]["dash"] == "dash"
    assert len(fit["x"]) == 4
    assert "t½" in plot["layout"]["annotations"][0]["text"]

    notes = " ".join(body["outputs"]["notes"]["items"])
    assert "BLQ" in notes and "C0" in notes and "linear trapezoidal rule" in notes and "steady-state" in notes


def test_route_notes_are_shown_in_the_cards():
    example = _module()["example"]
    for route, needles in [
        ("IV bolus", ("IV: this is CL", "IV: this is Vz", "AUMC0-∞ / AUC0-∞")),
        ("IV infusion", ("IV: this is CL", "T_inf/2", "not corrected")),
        ("Oral", ("apparent CL/F", "apparent Vz/F", "MAT + MRT")),
    ]:
        response = _run({**example, "route": route})
        assert response.status_code == 200, response.text
        summary = response.json()["outputs"]["summary"]["values"]
        shown = json.dumps(summary, ensure_ascii=False)
        for needle in needles:
            assert needle in shown, (route, needle)


def test_warnings_for_a_short_poor_terminal_window():
    response = _run({"profile": "0 100\n1 20\n2 30", "dose": 10, "route": "Oral", "lambda_points": 3})
    assert response.status_code == 200, response.text
    warnings = " ".join(response.json()["warnings"])
    assert "extrapolates" in warnings
    assert "R²" in warnings
    assert "only 3 points" in warnings


def test_bad_input_is_400_and_names_the_row():
    response = _run({"profile": "0 100\n1 nope\n2 25", "dose": 10})
    assert response.status_code == 400
    assert "Row 2" in response.json()["detail"]

    assert _run({"profile": "0 100\n1 50\n2 25"}).status_code == 400          # dose missing
    assert _run({"profile": "0 100\n1 50\n2 25", "dose": 10, "lambda_points": 2}).status_code == 400
    assert _run({"profile": "0 100\n1 50\n2 25", "dose": 10, "lambda_points": 2.5}).status_code == 400
    assert _run({"profile": "0 100\n1 50\n2 25", "dose": 10, "route": "Subcutaneous"}).status_code == 400
