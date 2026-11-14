"""Inorganic MS：同位素 / 干扰、IDMS、质量歧视和 spICP-MS 的数值与 API 合约。"""
import math
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from mzdesk.api.app import app
from mzdesk.inorganic import idms, massbias, single_particle
from mzdesk.inorganic.core import all_isotopes, isotope_label, parse_isotope
from mzdesk.inorganic.data import (ELEMENT_SYMBOLS, INTERNAL_STANDARDS, IONIZATION_EV, MATRIX_ELEMENTS,
                                   PARTICLE_PRESETS, SENSITIVITY_DWELL, SENSITIVITY_RATE)
from mzdesk.inorganic.interferences import find_interferences
from mzdesk.inorganic.isotopes import element_report, internal_standards
from mzdesk.spec import MODULES, SUBGROUPS

client = TestClient(app)
INORGANIC = [m for m in MODULES if m["subgroup"] == "inorganic"]
IDS = ["icp_isotopes", "icp_interferences", "icp_idms", "icp_massbias", "icp_sp"]


def _run(module_id, inputs):
    return client.post(f"/api/modules/{module_id}/run", json=inputs)


def _module(module_id):
    return next(m for m in MODULES if m["id"] == module_id)


# ---------------------------------------------------------------- page definitions and API shapes
def test_inorganic_subgroup_and_live_modules():
    ids = [s["id"] for s in SUBGROUPS]
    assert ids[ids.index("dmpk") + 1:] == ["inorganic", "gnps"]
    assert [m["id"] for m in INORGANIC] == IDS
    assert next(s["title"] for s in SUBGROUPS if s["id"] == "inorganic") == "ICP-MS"
    for module in INORGANIC:
        assert module["cls"] == "B" and module["live"] is True and module["example"]
        assert any(block.get("help") for block in module["outputs"]), module["id"]
    intensity = next(f for f in _module("icp_sp")["inputs"] if f.get("id") == "intensities")
    assert intensity["type"] == "textarea"
    icons = (Path(__file__).parents[1] / "web" / "src" / "assets" / "icons.svg").read_text()
    assert 'symbol id="inorganic"' in icons


@pytest.mark.parametrize("module", INORGANIC, ids=lambda m: m["id"])
def test_examples_run_and_fill_render_shaped_outputs(module):
    response = _run(module["id"], module["example"])
    assert response.status_code == 200, response.text
    body = response.json()
    assert body["status"] == "ok" and body["verdict"]
    expected = {block["id"] for block in module["outputs"] if block.get("id")}
    assert set(body["outputs"]) == expected
    for block in module["outputs"]:
        if not block.get("id"):
            continue
        data = body["outputs"][block["id"]]
        if block["type"] == "table":
            assert isinstance(data["rows"], list) and data.get("count", len(data["rows"])) >= len(data["rows"])
            for row in data["rows"]:
                cells = row["cells"] if isinstance(row, dict) else row
                assert len(cells) == len(block["columns"])
        elif block["type"] == "cards":
            assert {item["label"] for item in block["items"]} == set(data["values"])
        elif block["type"] == "figure":
            assert data["plotly"]["data"] and data["plotly"]["layout"]
        elif block["type"] == "notes":
            assert data["items"] and all(isinstance(item, str) for item in data["items"])


@pytest.mark.parametrize("module_id,inputs", [
    ("icp_isotopes", {"element": "Xx"}),
    ("icp_interferences", {"target": "not an isotope"}),
    ("icp_idms", {"element": "Ag", "isotope_a": "107Ag", "isotope_b": "109Ag", "a_y": "0.01", "b_y": "0.99",
                  "spike_concentration": "10", "spike_unit": "µg/g", "m_x": "1", "m_y": "1", "r_b": "-1"}),
    ("icp_massbias", {"mb_a": "107Ag"}),
    ("icp_sp", {"material": "Au", "element": "Au", "sensitivity": "1000", "sensitivity_unit": SENSITIVITY_DWELL,
                "dwell_ms": "1", "flow": "0.3", "transport_pct": "5", "background": "2", "intensities": "1",
                "n_events": "1", "acquisition_s": "60"}),
])
def test_bad_input_is_400(module_id, inputs):
    response = _run(module_id, inputs)
    assert response.status_code == 400 and response.json()["detail"]


# ---------------------------------------------------------------- isotope data and selection
def test_reference_data_cover_every_offered_element():
    assert len(ELEMENT_SYMBOLS) == 84 and set(IONIZATION_EV) == set(ELEMENT_SYMBOLS)
    assert all(value > 0 for value in IONIZATION_EV.values())
    assert IONIZATION_EV["As"] == 9.78855 and IONIZATION_EV["Lu"] == 5.425871 and IONIZATION_EV["Pa"] == 5.89
    assert {isotope.element for isotope in all_isotopes()} == set(ELEMENT_SYMBOLS)
    assert not {"D", "T", "E", "Me", "Pn"} & set(ELEMENT_SYMBOLS)
    assert [parse_isotope(label).label for label in INTERNAL_STANDARDS] == list(INTERNAL_STANDARDS)


def test_isotope_labels_and_element_report():
    assert parse_isotope("75As") == parse_isotope("As-75") == parse_isotope("⁷⁵As")
    assert isotope_label(parse_isotope("75As")) == "⁷⁵As"
    with pytest.raises(ValueError):
        parse_isotope("75Ge", "As")
    report = element_report("As")
    assert len(report["isotopes"]) == 1
    row = report["isotopes"][0]
    assert row["isotope"].label == "75As" and row["isotope"].abundance == 1
    assert not row["isobars"] and row["poly_count"] > 0 and row["strongest_poly"]
    assert report["recommendations"][0]["isotope"].label == "75As"
    assert "No atomic isobar" in report["recommendations"][0]["reason"]


def test_internal_standard_ranking_uses_mass_and_ionization_energy():
    result = internal_standards("As")
    rows = result["standards"]
    assert rows[0]["isotope"].label == "74Ge"
    assert [row["rank"] for row in rows] == list(range(1, len(rows) + 1))
    assert [row["score"] for row in rows] == sorted(row["score"] for row in rows)
    assert all(row["mass_difference"] >= 0 and row["ie_difference"] >= 0 for row in rows)
    assert all(row["isotope"].element != "Ge" for row in internal_standards("Ge")["standards"])


# ---------------------------------------------------------------- interference engine
def test_as75_arcl_interference_mass_and_resolution():
    rows = find_interferences("75As", MATRIX_ELEMENTS)["interferences"]
    row = next(row for row in rows if row["formula"] == "⁴⁰Ar³⁵Cl⁺")
    ar40, cl35 = parse_isotope("40Ar"), parse_isotope("35Cl")
    neutral_mass_sum = ar40.mass + cl35.mass
    assert neutral_mass_sum == pytest.approx(74.9312, abs=0.0001)
    assert row["mz"] == pytest.approx(74.9312, abs=0.0006)  # 再扣一个电子质量
    assert row["mz"] == pytest.approx(neutral_mass_sum - 0.00054857990946, rel=1e-12)
    assert 7.7e3 <= row["required_r"] <= 7.8e3
    assert row["delta_mda"] == pytest.approx(9.640152, abs=1e-6)
    assert row["score"] == pytest.approx(ar40.abundance * cl35.abundance, rel=1e-12)
    assert row["type"] == "polyatomic"


def test_atomic_isobar_doubly_charged_and_charge_selection():
    rows = find_interferences("64Zn", [], charge_states=(1, 2))["interferences"]
    assert any(row["formula"] == "⁶⁴Ni⁺" and row["type"] == "isobar" for row in rows)
    assert any(row["type"] == "doubly charged" and row["formula"].endswith("²⁺") for row in rows)
    singly = find_interferences("64Zn", [], charge_states=("1+",))["interferences"]
    doubly = find_interferences("64Zn", [], charge_states=("2+",))["interferences"]
    assert singly and all(row["type"] == "isobar" for row in singly)
    assert doubly and all(row["type"] == "doubly charged" for row in doubly)


def test_argon_dimer_trimer_oxides_and_max_atoms():
    rows = find_interferences("80Se", ("Ar", "H", "O"), 3, (1,))["interferences"]
    assert "⁴⁰Ar⁴⁰Ar⁺" in {row["formula"] for row in rows}
    oxide = find_interferences("56Fe", ("Ar", "Ca", "O"), 2, (1,))["interferences"]
    assert {"⁴⁰Ar¹⁶O⁺", "⁴⁰Ca¹⁶O⁺"} <= {row["formula"] for row in oxide}
    hydroxide = find_interferences("57Fe", ("Ar", "O", "H"), 3, (1,))["interferences"]
    assert "⁴⁰Ar¹⁶O¹H⁺" in {row["formula"] for row in hydroxide}
    hydride = find_interferences("41K", ("Ar", "H"), 2, (1,))["interferences"]
    assert "⁴⁰Ar¹H⁺" in {row["formula"] for row in hydride}
    rows79 = find_interferences("79Br", ("Ar", "H"), 3, (1,))["interferences"]
    trimer = next(row for row in rows79 if row["formula"] == "⁴⁰Ar³⁸Ar¹H⁺")
    expected_score = math.prod(parse_isotope(label).abundance for label in ("40Ar", "38Ar", "1H"))
    assert trimer["score"] == pytest.approx(expected_score, rel=1e-12)
    pairs = find_interferences("75As", ("Ar", "Cl"), 2, (1,))["interferences"]
    assert any(row["formula"] == "⁴⁰Ar³⁵Cl⁺" for row in pairs)
    assert len(find_interferences("79Br", ("Ar", "H"), 2, (1,))["interferences"]) < len(rows79)


def test_interference_sorting_and_api_cap():
    rows = find_interferences("80Se")["interferences"]
    assert [row["score"] for row in rows] == sorted((row["score"] for row in rows), reverse=True)
    assert len(rows) > 50
    body = _run("icp_interferences", {"target": "80Se", "matrix": list(MATRIX_ELEMENTS), "max_atoms": "3",
                                      "charge_states": ["1+", "2+"]}).json()
    assert len(body["outputs"]["interferences"]["rows"]) == 50
    assert body["outputs"]["interferences"]["count"] == len(rows)
    assert body["outputs"]["summary"]["values"]["Interferences found"]["value"] == str(len(rows))


def test_interference_notes_do_not_recommend_cell_conditions():
    body = _run("icp_interferences", _module("icp_interferences")["example"]).json()
    notes = " ".join(body["outputs"]["notes"]["items"])
    assert "not an intensity" in notes and "No collision/reaction-cell recommendation" in notes
    arcl = next(row for row in body["outputs"]["interferences"]["rows"] if row[0] == "⁴⁰Ar³⁵Cl⁺")
    assert arcl[2] == "74.9306873" and 7.7e3 < float(arcl[4]) < 7.8e3


# ---------------------------------------------------------------- IDMS
def test_idms_blend_round_trip_and_mass_concentration_unit():
    element, a, b = "Ag", parse_isotope("107Ag"), parse_isotope("109Ag")
    ax, bx, ay, by = a.abundance, b.abundance, 0.01, 0.99
    expected_sample_mass = ax * a.mass + bx * b.mass
    expected_spike_mass = ay * a.mass + by * b.mass
    expected_cx, cy, mx, my = 2.3456789, 8.7654321, 0.7, 0.2
    nx, ny = expected_cx * mx, cy * my
    ratio_b = (ax * nx + ay * ny) / (bx * nx + by * ny)
    result = idms.calculate(element, a.label, b.label, ax, bx, ay, by, cy, "µmol/g", mx, my, ratio_b)
    assert result["concentration_umol_g"] == pytest.approx(expected_cx, rel=1e-9)
    assert result["sample_molar_mass"] == pytest.approx(expected_sample_mass, rel=1e-12)
    assert result["spike_molar_mass"] == pytest.approx(expected_spike_mass, rel=1e-12)
    assert result["concentration_ug_g"] == pytest.approx(expected_cx * expected_sample_mass, rel=1e-12)
    assert result["ratio_opt"] == pytest.approx(math.sqrt((ax / bx) * (ay / by)), rel=1e-12)
    assert result["sensitivity"] > 0

    spike_ug_g = cy * expected_spike_mass
    converted = idms.calculate(element, a.label, b.label, None, None, ay, by, spike_ug_g, "µg/g", mx, my, ratio_b)
    assert converted["spike_amount_umol_g"] == pytest.approx(cy, rel=1e-12)
    assert converted["concentration_umol_g"] == pytest.approx(expected_cx, rel=1e-9)
    assert converted["fraction_ax"] == a.abundance and converted["fraction_bx"] == b.abundance


def test_idms_rejects_undefined_or_unphysical_blends():
    a, b = parse_isotope("107Ag"), parse_isotope("109Ag")
    with pytest.raises(ValueError, match="undefined"):
        idms.calculate("Ag", a.label, b.label, None, None, 0.01, 0.99, 10, "µmol/g", 1, 1, a.abundance / b.abundance)
    with pytest.raises(ValueError, match="negative"):
        idms.calculate("Ag", a.label, b.label, None, None, 0.01, 0.99, 10, "µmol/g", 1, 1, 2)
    with pytest.raises(ValueError, match="sum"):
        idms.calculate("Ag", a.label, b.label, None, None, 0.5, 0.6, 10, "µmol/g", 1, 1, 1)


# ---------------------------------------------------------------- mass bias and delta
@pytest.mark.parametrize("law,parameter", [("Linear", 0.001), ("Power", 0.002), ("Exponential (Russell)", -1.5)])
def test_mass_bias_laws_round_trip(law, parameter):
    a, b = parse_isotope("107Ag"), parse_isotope("109Ag")
    if law == "Linear":
        factor = 1 + parameter * (a.mass - b.mass)
    elif law == "Power":
        factor = (1 + parameter) ** (a.mass - b.mass)
    else:
        factor = (a.mass / b.mass) ** parameter
    true_standard, true_sample = 1.08, 0.91
    result = massbias.calibrate(a.label, b.label, true_standard / factor, true_standard, true_sample / factor)
    row = next(row for row in result["laws"] if row["law"] == law)
    assert result["k"] == pytest.approx(factor, rel=1e-12)
    assert row["parameter"] == pytest.approx(parameter, rel=1e-10, abs=1e-12)
    assert row["corrected_sample"] == pytest.approx(true_sample, rel=1e-12)


def test_bracketing_delta_and_inverse():
    corrected = massbias.bracketing(0.995, 0.990, 0.992, 1.0)
    assert corrected == pytest.approx(0.995 / 0.991, rel=1e-12)
    delta = massbias.delta_value(1.002, 1.0)
    assert delta == pytest.approx(2.0, abs=1e-12)
    assert massbias.ratio_from_delta(delta) == pytest.approx(1.002, rel=1e-12)
    assert massbias.ratio_from_delta(-500) == 0.5
    with pytest.raises(ValueError):
        massbias.ratio_from_delta(-1000)


def test_mass_bias_api_runs_all_three_calculators():
    body = _run("icp_massbias", _module("icp_massbias")["example"]).json()
    outputs = body["outputs"]
    assert [row[0] for row in outputs["corrected_ratios"]["rows"]] == list(massbias.LAWS)
    assert outputs["bracketing"]["values"]["Corrected sample ratio"] != "—"
    assert outputs["delta"]["values"]["δ from ratios"]["value"].startswith("2")
    assert outputs["delta"]["values"]["R_sample from δ"]["value"].startswith("1.002")
    blank = _run("icp_massbias", {})
    assert blank.status_code == 200 and blank.json()["outputs"]["corrected_ratios"]["rows"] == []


# ---------------------------------------------------------------- single-particle ICP-MS
def test_particle_presets_and_tio2_mass_fraction():
    expected = {"Au": 19.30, "Ag": 10.49, "Pt": 21.45, "TiO₂ rutile": 4.23, "TiO₂ anatase": 3.89,
                "CeO₂": 7.22, "ZnO": 5.61, "SiO₂ amorphous": 2.20, "Fe₃O₄": 5.17, "CuO": 6.31}
    assert {name: density for name, (_, density) in PARTICLE_PRESETS.items()} == expected
    assert single_particle.mass_fraction("TiO2", "Ti") == pytest.approx(0.5993, abs=0.0001)
    rutile = single_particle.resolve_material("TiO₂ rutile", measured_element="Ti")
    assert rutile["density"] == 4.23 and rutile["mass_fraction"] == pytest.approx(0.5993, abs=0.0001)
    custom = single_particle.resolve_material("Custom", "Fe3O4", 5.17, "Fe")
    assert custom["mass_fraction"] == pytest.approx(0.7236, abs=0.0002)
    with pytest.raises(ValueError, match="not present"):
        single_particle.mass_fraction("TiO2", "Au")


@pytest.mark.parametrize("unit,sensitivity", [(SENSITIVITY_DWELL, 1000.0), (SENSITIVITY_RATE, 100000.0)])
def test_au_50nm_forward_signal_recovers_size_and_number_concentration(unit, sensitivity):
    material = single_particle.resolve_material("Au", measured_element="Au")
    dwell_ms, flow, eta_pct, background = 1.0, 0.3, 5.0, 2.0
    particle_mass_g = math.pi / 6 * (50.0e-7) ** 3 * 19.30
    assert particle_mass_g == pytest.approx(1.263e-15, rel=0.001)
    sensitivity_dwell = sensitivity if unit == SENSITIVITY_DWELL else sensitivity * dwell_ms / 1000
    eta = eta_pct / 100
    intensity = particle_mass_g * 1e6 * sensitivity_dwell / (eta * (flow / 60000) * (dwell_ms / 1000)) + background
    threshold = background + (intensity - background) * 0.25
    result = single_particle.calculate(material, sensitivity, unit, dwell_ms, flow, eta_pct, background,
                                       [intensity], 1, 60, threshold)
    event = result["events"][0]
    assert event["particle_mass_g"] == pytest.approx(particle_mass_g, rel=1e-12)
    assert event["diameter_nm"] == pytest.approx(50.0, rel=1e-10)
    assert result["mean_diameter_nm"] == result["median_diameter_nm"] == pytest.approx(50.0, rel=1e-10)
    assert result["number_concentration_ml"] == pytest.approx(1 / (eta * (flow / 60) * 60), rel=1e-12)
    assert result["number_concentration_l"] == pytest.approx(result["number_concentration_ml"] * 1000, rel=1e-12)
    assert result["size_lod_nm"] == pytest.approx(50.0 * 0.25 ** (1 / 3), rel=1e-10)


def test_sp_api_histogram_events_lod_and_conversion_only_scope():
    body = _run("icp_sp", _module("icp_sp")["example"]).json()
    assert body["status"] == "ok"
    outputs = body["outputs"]
    assert len(outputs["events"]["rows"]) == 3 and outputs["events"]["count"] == 3
    assert outputs["summary"]["values"]["Number concentration"]["value"] == "200"
    assert outputs["summary"]["values"]["Number concentration"]["note"] == "200000 particles/L"
    assert outputs["summary"]["values"]["Size LOD"] != "—"
    assert outputs["histogram"]["plotly"]["data"][0]["type"] == "histogram"
    assert len(outputs["histogram"]["plotly"]["data"][0]["x"]) == 3
    notes = " ".join(outputs["notes"]["items"])
    assert "No raw time-resolved processing" in notes and "Pace et al." in notes


def test_sp_rejects_bad_event_lists_and_counts():
    material = single_particle.resolve_material("Au", measured_element="Au")
    with pytest.raises(ValueError, match="Event 1"):
        single_particle.calculate(material, 1000, SENSITIVITY_DWELL, 1, 0.3, 5, 2, [2], 1, 60)
    with pytest.raises(ValueError, match="N must be an integer"):
        single_particle.calculate(material, 1000, SENSITIVITY_DWELL, 1, 0.3, 5, 2, [5000], 0.5, 60)
    with pytest.raises(ValueError, match="Threshold"):
        single_particle.calculate(material, 1000, SENSITIVITY_DWELL, 1, 0.3, 5, 2, [5000], 1, 60, 1)
    with pytest.raises(ValueError, match="positive atom counts"):
        single_particle.resolve_material("Custom", "Ti0O2", 4.23, "Ti")
