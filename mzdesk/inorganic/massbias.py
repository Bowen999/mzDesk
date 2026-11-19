"""同位素比质量歧视校正、标准–样品 bracketing 与 δ 值。"""
import math

from mzdesk.inorganic.core import parse_isotope, require_finite

LAWS = ("Linear", "Power", "Exponential (Russell)")


def _positive(value, label):
    value = require_finite(float(value), label)
    if value <= 0:
        raise ValueError(f"{label} must be positive.")
    return value


def correction_factor(law, mass_a, mass_b, parameter):
    """R_true / R_meas；linear / power 参数为每 Da 的 ε，Russell 参数为 β。"""
    mass_a = _positive(mass_a, "Isotope mass m_a")
    mass_b = _positive(mass_b, "Isotope mass m_b")
    parameter = require_finite(float(parameter), "Mass-bias parameter")
    difference = mass_a - mass_b
    try:
        if law == "Linear":
            factor = 1 + parameter * difference
        elif law == "Power":
            factor = (1 + parameter) ** difference
        elif law == "Exponential (Russell)":
            factor = (mass_a / mass_b) ** parameter
        else:
            raise ValueError(f"Unknown mass-bias law {law!r}.")
    except (OverflowError, ZeroDivisionError):
        raise ValueError("The mass-bias parameter gives an undefined correction.") from None
    if not math.isfinite(factor) or factor <= 0:
        raise ValueError("The mass-bias parameter gives a non-positive or non-finite correction factor.")
    return factor


def calibrate(a_label, b_label, measured_standard, true_standard, measured_sample):
    a = parse_isotope(a_label)
    b = parse_isotope(b_label)
    if a == b:
        raise ValueError("The numerator and denominator isotopes must be different.")
    measured_standard = _positive(measured_standard, "Measured standard ratio")
    true_standard = _positive(true_standard, "Certified standard ratio")
    measured_sample = _positive(measured_sample, "Measured sample ratio")
    k = true_standard / measured_standard
    difference = a.mass - b.mass
    if difference == 0 or a.mass == b.mass:
        raise ValueError("The isotope masses must be different.")
    parameters = {
        "Linear": (k - 1) / difference,
        "Power": k ** (1 / difference) - 1,
        "Exponential (Russell)": math.log(k) / math.log(a.mass / b.mass),
    }
    rows = []
    for law in LAWS:
        parameter = parameters[law]
        factor = correction_factor(law, a.mass, b.mass, parameter)
        rows.append({
            "law": law,
            "parameter_name": "β" if law == "Exponential (Russell)" else "ε",
            "parameter": parameter,
            "factor": factor,
            "corrected_sample": measured_sample * factor,
        })
    return {"a": a, "b": b, "k": k, "laws": rows}


def bracketing(measured_sample, measured_before, measured_after, true_standard):
    measured_sample = _positive(measured_sample, "Measured sample ratio")
    before = _positive(measured_before, "Standard ratio before sample")
    after = _positive(measured_after, "Standard ratio after sample")
    true_standard = _positive(true_standard, "Certified standard ratio")
    return measured_sample * true_standard / ((before + after) / 2)


def delta_value(sample_ratio, standard_ratio):
    sample_ratio = _positive(sample_ratio, "Sample ratio")
    standard_ratio = _positive(standard_ratio, "Standard ratio")
    return (sample_ratio / standard_ratio - 1) * 1000


def ratio_from_delta(delta):
    """由 δ 求 R_sample / R_standard；再乘标准比即得 R_sample。"""
    delta = require_finite(float(delta), "δ value")
    if delta <= -1000:
        raise ValueError("δ must be greater than −1000‰ to give a positive ratio.")
    return 1 + delta / 1000
