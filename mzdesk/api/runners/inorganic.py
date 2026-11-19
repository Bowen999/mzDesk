"""Inorganic MS 五页：同位素 / 干扰 / IDMS / 质量歧视 / 单颗粒 ICP-MS。"""
from functools import wraps
import math
import re

from mzdesk.api import blocks, common
from mzdesk.api.common import InputError


def _input_errors(function):
    @wraps(function)
    def wrapped(inputs):
        try:
            return function(inputs)
        except InputError:
            raise
        except ValueError as exc:
            raise InputError(str(exc)) from None
    return wrapped


def _g(value, digits=8):
    return common.MISSING if value is None or not math.isfinite(value) else f"{value:.{digits}g}"


def _sci(value, digits=4):
    return common.MISSING if value is None or not math.isfinite(value) else f"{value:.{digits}e}"


def _selection(inputs, key, default):
    value = inputs.get(key)
    if value is None:
        return list(default)
    return [value] if isinstance(value, str) else list(value)


def _integer(inputs, key, label, default):
    value = common.number(inputs, key, label, default=default)
    if value != int(value):
        raise InputError(f"{label} must be an integer.")
    return int(value)


def _result(outputs, verdict, warnings=None):
    return {"outputs": outputs, "warnings": warnings or [], "verdict": verdict}


@_input_errors
def run_isotopes(inputs):
    from mzdesk.inorganic import isotopes as engine
    from mzdesk.inorganic.data import MATRIX_ELEMENTS

    element = common.text(inputs, "element") or "As"
    matrix = _selection(inputs, "matrix", MATRIX_ELEMENTS)
    report = engine.element_report(element, matrix)
    standards = engine.internal_standards(element)
    recommended = report["recommendations"][0]
    top_standard = standards["standards"][0]
    isotope_rows = []
    for row in report["isotopes"]:
        isotope = row["isotope"]
        strongest = row["strongest_poly"]
        isotope_rows.append([
            isotope.display, common.fmt(isotope.mass, 7), common.fmt(isotope.abundance * 100, 6),
            ", ".join(i.display for i in row["isobars"]) or common.MISSING,
            str(row["poly_count"]),
            f"{strongest['formula']} (score {_g(strongest['score'], 4)})" if strongest else common.MISSING,
        ])
    recommendation_rows = [[str(row["rank"]), row["isotope"].display, _g(row["score"], 6), row["reason"]]
                           for row in report["recommendations"]]
    standard_rows = [[str(row["rank"]), row["isotope"].display, common.fmt(row["mass_difference"], 4),
                      common.fmt(row["ionization_energy"], 6), common.fmt(row["ie_difference"], 4), _g(row["score"], 5)]
                     for row in standards["standards"]]
    outputs = {
        "summary": {"values": {
            "Analyte": {"value": element, "note": f"average mass {standards['analyte_mass']:.5f} Da"},
            "Recommended isotope": {"value": recommended["isotope"].display,
                                    "note": "no atomic isobar" if not recommended["isobars"] else "best scored despite atomic isobar"},
            "Top internal standard": {"value": top_standard["isotope"].display,
                                      "note": f"|Δmass| {top_standard['mass_difference']:.3f} Da · |ΔIE| {top_standard['ie_difference']:.3f} eV"},
            "Analyte first IE": common.fmt(standards["analyte_ie"], 6),
        }},
        "isotope_table": blocks.table_data(isotope_rows, count=len(isotope_rows)),
        "recommendations": blocks.table_data(recommendation_rows, count=len(recommendation_rows)),
        "internal_standards": blocks.table_data(standard_rows, count=len(standard_rows)),
    }
    verdict = (f"<b>{recommended['isotope'].display}</b> is the top-ranked analyte isotope; "
               f"<b>{top_standard['isotope'].display}</b> is the closest internal-standard match.")
    return _result(outputs, verdict)


@_input_errors
def run_interferences(inputs):
    from mzdesk.inorganic.data import MATRIX_ELEMENTS
    from mzdesk.inorganic.interferences import find_interferences

    target_text = common.text(inputs, "target") or "75As"
    matrix = _selection(inputs, "matrix", MATRIX_ELEMENTS)
    max_atoms = _integer(inputs, "max_atoms", "Maximum atoms per polyatomic ion", 3)
    charges = _selection(inputs, "charge_states", ("1+", "2+"))
    result = find_interferences(target_text, matrix, max_atoms, charges)
    target, rows = result["target"], result["interferences"]
    strongest = rows[0] if rows else None
    exact = next((row for row in rows if row["required_r"] is None), None)
    highest = exact or (max(rows, key=lambda row: row["required_r"]) if rows else None)
    highest_value = "∞" if exact else _sci(highest["required_r"], 4) if highest else common.MISSING
    kind_names = {"isobar": "Atomic isobar", "doubly charged": "M²⁺", "polyatomic": "Polyatomic"}
    outputs = {
        "summary": {"values": {
            "Interferences found": {"value": str(len(rows)), "note": f"showing {min(50, len(rows))} in the table"},
            "Strongest": {"value": strongest["formula"], "note": f"likelihood score {_g(strongest['score'], 5)}"} if strongest else common.MISSING,
            "Highest resolving power": {"value": highest_value, "note": highest["formula"] + (" · exact mass coincidence" if exact else "")} if highest else common.MISSING,
        }},
        "interferences": blocks.table_data([
            [row["formula"], kind_names[row["type"]], common.fmt(row["mz"], 7), common.signed(row["delta_mda"], 4),
             "∞" if row["required_r"] is None else _sci(row["required_r"], 4), _sci(row["score"], 5)]
            for row in rows[:50]
        ], count=len(rows)),
        "notes": {"items": [
            f"Target {target.display}⁺: exact m/z {target.mass - 0.00054857990946:.7f}; all candidates match its nominal mass.",
            "Likelihood score = product of the constituent isotopes' natural abundances. It is not an intensity, ion-formation probability, or cell-response prediction.",
            "Exact m/z subtracts one electron mass per positive charge; required R = target m/z ÷ |candidate m/z − target m/z|.",
            "No collision/reaction-cell recommendation is made in this version.",
        ] + ([f"Table capped at the top 50 of {len(rows)} candidates, sorted by likelihood score."] if len(rows) > 50 else [])},
    }
    verdict = (f"{len(rows)} interference{'s' if len(rows) != 1 else ''} for <b>{target.display}⁺</b>"
               + (f"; strongest likelihood score: <b>{strongest['formula']}</b>." if strongest else "."))
    return _result(outputs, verdict)


@_input_errors
def run_idms(inputs):
    from mzdesk.inorganic import idms

    element = common.text(inputs, "element") or "Ag"
    result = idms.calculate(
        element, common.text(inputs, "isotope_a"), common.text(inputs, "isotope_b"),
        common.number(inputs, "a_x", "Sample atom fraction A_x", minimum=0, maximum=1),
        common.number(inputs, "b_x", "Sample atom fraction B_x", minimum=0, maximum=1),
        common.number(inputs, "a_y", "Spike atom fraction A_y", required=True, minimum=0, maximum=1),
        common.number(inputs, "b_y", "Spike atom fraction B_y", required=True, minimum=0, maximum=1),
        common.number(inputs, "spike_concentration", "Spike concentration", required=True, positive=True),
        common.text(inputs, "spike_unit") or "µg/g",
        common.number(inputs, "m_x", "Sample mass m_x", required=True, positive=True),
        common.number(inputs, "m_y", "Spike mass m_y", required=True, positive=True),
        common.number(inputs, "r_b", "Measured blend ratio R_B", required=True, positive=True),
    )
    intermediates = [
        ("Reference isotope a", result["a"].display, "A/B ratio numerator"),
        ("Spike isotope b", result["b"].display, "A/B ratio denominator"),
        ("A_x", _g(result["fraction_ax"], 10), "sample atom fraction"),
        ("B_x", _g(result["fraction_bx"], 10), "sample atom fraction"),
        ("A_y", _g(result["fraction_ay"], 10), "spike atom fraction"),
        ("B_y", _g(result["fraction_by"], 10), "spike atom fraction"),
        ("R_x = A_x/B_x", _g(result["ratio_x"], 10), "sample isotope ratio"),
        ("R_y = A_y/B_y", _g(result["ratio_y"], 10), "spike isotope ratio"),
        ("c_y", _g(result["spike_amount_umol_g"], 10), "µmol/g"),
        ("M_x", _g(result["sample_molar_mass"], 10), "g/mol, from sample composition"),
        ("M_y", _g(result["spike_molar_mass"], 10), "g/mol, from spike composition"),
        ("m_y/m_x", _g(result["mass_ratio"], 10), "masses in the same unit"),
        ("A_y − R_B·B_y", _g(result["numerator"], 10), "numerator"),
        ("R_B·B_x − A_x", _g(result["denominator"], 10), "denominator"),
    ]
    outputs = {
        "summary": {"values": {
            "Sample concentration": {"value": _g(result["concentration_ug_g"], 10), "note": _g(result["concentration_umol_g"], 10) + " µmol/g"},
            "Sample amount content": _g(result["concentration_umol_g"], 10),
            "Error magnification": {"value": _g(result["sensitivity"], 6), "note": "|∂ln(c_x)/∂ln(R_B)|"},
            "Optimal blend R_opt": {"value": _g(result["ratio_opt"], 8), "note": "√(R_x·R_y)"},
        }},
        "intermediates": blocks.table_data(intermediates),
        "notes": {"items": [
            "c_x = c_y·(m_y/m_x)·(A_y − R_B·B_y)/(R_B·B_x − A_x); c_x and c_y are amount contents in µmol/g.",
            "R_B must already be mass-bias corrected. Sample and spike molar masses are calculated from their isotopic compositions.",
            "Isotopes other than a and b retain their natural relative proportions, scaled to the atom fraction left over after A + B.",
            "Error magnification is the dimensionless relative sensitivity |∂c_x/∂R_B · R_B/c_x|; it is not a complete uncertainty budget.",
        ]},
    }
    verdict = (f"Sample concentration: <b>{_g(result['concentration_ug_g'], 8)} µg/g</b> "
               f"({_g(result['concentration_umol_g'], 8)} µmol/g).")
    return _result(outputs, verdict)


def _mass_bias_section(inputs):
    labels = (("mb_a", "Isotope a"), ("mb_b", "Isotope b"))
    texts = [common.text(inputs, key) for key, _ in labels]
    numbers = [common.number(inputs, key, label, positive=True) for key, label in (
        ("mb_measured_std", "Measured standard ratio"), ("mb_true_std", "Certified standard ratio"),
        ("mb_measured_sample", "Measured sample ratio"))]
    if not any(texts + numbers):
        return None
    if not all(texts) or any(value is None for value in numbers):
        raise InputError("Complete isotopes a/b, measured and certified standard ratios, and measured sample ratio for mass-bias correction.")
    return texts + numbers


def _bracketing_section(inputs):
    values = [common.number(inputs, key, label, positive=True) for key, label in (
        ("bracket_sample", "Bracketing measured sample ratio"), ("bracket_before", "Standard ratio before sample"),
        ("bracket_after", "Standard ratio after sample"), ("bracket_true", "Bracketing certified standard ratio"))]
    if not any(value is not None for value in values):
        return None
    if any(value is None for value in values):
        raise InputError("Complete all four standard–sample bracketing ratios.")
    return values


@_input_errors
def run_massbias(inputs):
    from mzdesk.inorganic import massbias

    outputs = {}
    verdict = None
    section = _mass_bias_section(inputs)
    if section:
        result = massbias.calibrate(*section)
        parameters = {row["law"]: row["parameter"] for row in result["laws"]}
        outputs["mass_bias"] = {"values": {
            "K = R_true/R_meas": _g(result["k"], 10),
            "Linear ε": _g(parameters["Linear"], 8),
            "Power ε": _g(parameters["Power"], 8),
            "Russell β": _g(parameters["Exponential (Russell)"], 8),
        }}
        outputs["corrected_ratios"] = blocks.table_data([
            [row["law"], row["parameter_name"], _g(row["parameter"], 10), _g(row["factor"], 10), _g(row["corrected_sample"], 10)]
            for row in result["laws"]
        ])
        russell = next(row for row in result["laws"] if row["law"] == "Exponential (Russell)")
        verdict = f"Russell-corrected sample ratio: <b>{_g(russell['corrected_sample'], 10)}</b>."
    else:
        outputs["mass_bias"] = {"values": {key: common.MISSING for key in ("K = R_true/R_meas", "Linear ε", "Power ε", "Russell β")}}
        outputs["corrected_ratios"] = blocks.table_data([])

    bracket = _bracketing_section(inputs)
    if bracket:
        corrected = massbias.bracketing(*bracket)
        mean_standard = (bracket[1] + bracket[2]) / 2
        outputs["bracketing"] = {"values": {
            "Corrected sample ratio": _g(corrected, 10),
            "Mean measured standard": _g(mean_standard, 10),
        }}
        verdict = verdict or f"Bracketing-corrected sample ratio: <b>{_g(corrected, 10)}</b>."
    else:
        outputs["bracketing"] = {"values": {"Corrected sample ratio": common.MISSING, "Mean measured standard": common.MISSING}}

    delta_sample = common.number(inputs, "delta_sample", "δ sample ratio", positive=True)
    delta_standard = common.number(inputs, "delta_standard", "δ standard ratio", positive=True)
    delta_input = common.number(inputs, "delta_per_mil", "δ value")
    if (delta_sample is None) != (delta_standard is None):
        raise InputError("Provide both sample and standard ratios to calculate δ.")
    delta = massbias.delta_value(delta_sample, delta_standard) if delta_sample is not None else None
    relative = massbias.ratio_from_delta(delta_input) if delta_input is not None else None
    absolute = relative * delta_standard if relative is not None and delta_standard is not None else None
    outputs["delta"] = {"values": {
        "δ from ratios": {"value": _g(delta, 8), "note": "(R_sample/R_standard − 1) × 1000"} if delta is not None else common.MISSING,
        "R_sample/R_standard from δ": _g(relative, 10),
        "R_sample from δ": {"value": _g(absolute, 10), "note": "uses the δ standard ratio above"} if absolute is not None else common.MISSING,
    }}
    if verdict is None and delta is not None:
        verdict = f"δ = <b>{_g(delta, 8)}‰</b>."
    elif verdict is None and relative is not None:
        verdict = f"R_sample/R_standard = <b>{_g(relative, 10)}</b> from the supplied δ value."
    return _result(outputs, verdict or "Enter at least one complete calculator section.")


def _intensities(inputs):
    raw = common.text(inputs, "intensities")
    if not raw:
        raise InputError("Particle event intensities are required; paste one count value per line or a single value.")
    tokens = [token for token in re.split(r"[\s,;]+", raw) if token]
    try:
        return [float(token) for token in tokens]
    except ValueError:
        raise InputError("Particle event intensities must be numeric count values, one per line.") from None


def _histogram(diameters):
    import plotly.graph_objects as go
    from mzdesk.shared import colors
    from mzdesk.shared.ms2_plot import base_layout

    fig = go.Figure(go.Histogram(
        x=diameters, autobinx=True, marker=dict(color=colors.DATA_BLUE, line=dict(color=colors.PAPER, width=0.5)),
        hovertemplate="Diameter %{x:.2f} nm<br>Events %{y}<extra></extra>",
    ))
    fig.update_layout(**base_layout(height=380, x_range=None, y_range=None, y_title="Events"))
    fig.update_xaxes(title_text="Diameter (nm)", hoverformat=".2f")
    fig.update_yaxes(rangemode="tozero")
    return common.plotly(fig)


@_input_errors
def run_sp(inputs):
    from mzdesk.inorganic import single_particle
    from mzdesk.inorganic.data import SENSITIVITY_DWELL

    preset = common.text(inputs, "material") or "Au"
    element = common.text(inputs, "element") or "Au"
    material = single_particle.resolve_material(
        preset, common.text(inputs, "custom_formula"),
        common.number(inputs, "density", "Custom density", positive=True), element)
    intensities = _intensities(inputs)
    n_events = common.number(inputs, "n_events", "Number of events N", default=len(intensities), minimum=0)
    result = single_particle.calculate(
        material,
        common.number(inputs, "sensitivity", "Dissolved calibration sensitivity S", required=True, positive=True),
        common.text(inputs, "sensitivity_unit") or SENSITIVITY_DWELL,
        common.number(inputs, "dwell_ms", "Dwell time", required=True, positive=True),
        common.number(inputs, "flow", "Sample flow q", required=True, positive=True),
        common.number(inputs, "transport_pct", "Transport efficiency η", required=True, positive=True, maximum=100),
        common.number(inputs, "background", "Background counts per dwell", default=0, minimum=0),
        intensities, n_events,
        common.number(inputs, "acquisition_s", "Acquisition time", required=True, positive=True),
        common.number(inputs, "threshold", "Threshold intensity", positive=True),
    )
    events = result["events"]
    diameters = [event["diameter_nm"] for event in events]
    outputs = {
        "summary": {"values": {
            "Mean diameter": _g(result["mean_diameter_nm"], 7),
            "Median diameter": _g(result["median_diameter_nm"], 7),
            "Number concentration": {"value": _g(result["number_concentration_ml"], 7),
                                     "note": _g(result["number_concentration_l"], 7) + " particles/L"},
            "Size LOD": {"value": _g(result["size_lod_nm"], 7), "note": "from the supplied threshold"} if result["size_lod_nm"] is not None else common.MISSING,
        }},
        "events": blocks.table_data([
            [str(event["event"]), _g(event["intensity"], 8), _g(event["net_intensity"], 8),
             _sci(event["element_mass_ug"] * 1e9, 5), _sci(event["particle_mass_g"] * 1e15, 5), _g(event["diameter_nm"], 7)]
            for event in events[:200]
        ], count=len(events)),
        "histogram": _histogram(diameters),
        "notes": {"items": [
            f"{material['name']}: composition {material['formula']}, density {material['density']:g} g/cm³, measured-element mass fraction X = {material['mass_fraction']:.7f}.",
            "Pace et al., Anal. Chem. 2011, 83, 9361: dissolved-standard calibration plus transport efficiency. S is first converted to counts per dwell per µg/L when supplied as counts/s per µg/L.",
            "m_el[µg] = (I − I_bg)·η·q[L/s]·t_dw[s]/S; m_p = m_el/X; d = (6m_p/(πρ))^(1/3), assuming compact spheres and bulk density.",
            "C_N = N/(η·q·t_acq), with q converted from mL/min to mL/s. Number concentration is reported per mL and per L.",
            "Only supplied event signals are converted. No raw time-resolved processing, event detection, threshold optimization, dissolved-background fitting, coating, or dissolution correction is performed.",
            "Size LOD is only the diameter corresponding to the optional threshold intensity; it is not estimated from the event list.",
        ] + ([f"Event table capped at 200 of {len(events)} events; summary statistics and the histogram use all events."] if len(events) > 200 else [])},
    }
    verdict = (f"Mean diameter <b>{_g(result['mean_diameter_nm'], 6)} nm</b>; "
               f"number concentration <b>{_g(result['number_concentration_ml'], 6)} particles/mL</b>.")
    return _result(outputs, verdict)


RUNNERS = {
    "icp_isotopes": run_isotopes,
    "icp_interferences": run_interferences,
    "icp_idms": run_idms,
    "icp_massbias": run_massbias,
    "icp_sp": run_sp,
}
