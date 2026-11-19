"""单尖峰同位素稀释质谱（IDMS），全部中间量使用摩尔和原子分数。"""
import math

from mzdesk.inorganic.core import isotopic_molar_mass, parse_isotope, require_finite


def _fraction(value, natural, label):
    value = natural if value is None else require_finite(float(value), label)
    if not 0 <= value <= 1:
        raise ValueError(f"{label} must be between 0 and 1.")
    return value


def calculate(element, a_label, b_label, fraction_ax, fraction_bx, fraction_ay, fraction_by,
              spike_concentration, spike_unit, mass_x, mass_y, ratio_b):
    a = parse_isotope(a_label, element)
    b = parse_isotope(b_label, element)
    if a == b:
        raise ValueError("Reference isotope a and spike isotope b must be different.")
    ax = _fraction(fraction_ax, a.abundance, "Sample atom fraction A_x")
    bx = _fraction(fraction_bx, b.abundance, "Sample atom fraction B_x")
    ay = _fraction(fraction_ay, 0.0, "Spike atom fraction A_y")
    by = _fraction(fraction_by, 0.0, "Spike atom fraction B_y")
    if ax + bx > 1 + 1e-10 or ay + by > 1 + 1e-10:
        raise ValueError("The a and b atom fractions in each material cannot sum to more than 1.")
    if bx <= 0 or by <= 0:
        raise ValueError("B_x and B_y must be positive to define the isotope ratios and optimal blend.")
    spike_concentration = require_finite(float(spike_concentration), "Spike concentration")
    mass_x = require_finite(float(mass_x), "Sample mass m_x")
    mass_y = require_finite(float(mass_y), "Spike mass m_y")
    ratio_b = require_finite(float(ratio_b), "Measured blend ratio R_B")
    if spike_concentration <= 0 or mass_x <= 0 or mass_y <= 0 or ratio_b <= 0:
        raise ValueError("Spike concentration, m_x, m_y and R_B must be positive.")

    molar_x = isotopic_molar_mass(element, a, b, ax, bx)
    molar_y = isotopic_molar_mass(element, a, b, ay, by)
    if spike_unit == "µg/g":
        spike_umol_g = spike_concentration / molar_y
    elif spike_unit == "µmol/g":
        spike_umol_g = spike_concentration
    else:
        raise ValueError("Spike concentration unit must be µg/g or µmol/g.")
    numerator = ay - ratio_b * by
    denominator = ratio_b * bx - ax
    if denominator == 0:
        raise ValueError("R_B·B_x − A_x is zero; the IDMS equation is undefined.")
    amount_x = spike_umol_g * (mass_y / mass_x) * numerator / denominator
    if amount_x < 0:
        raise ValueError("The supplied fractions and R_B give a negative sample concentration; check isotope assignments and inputs.")
    ratio_x = ax / bx
    ratio_y = ay / by
    sensitivity = None
    if numerator != 0:
        sensitivity = abs(ratio_b * ((-by * denominator - numerator * bx) / (numerator * denominator)))
    return {
        "concentration_umol_g": amount_x,
        "concentration_ug_g": amount_x * molar_x,
        "sample_molar_mass": molar_x,
        "spike_molar_mass": molar_y,
        "spike_amount_umol_g": spike_umol_g,
        "mass_ratio": mass_y / mass_x,
        "numerator": numerator,
        "denominator": denominator,
        "ratio_x": ratio_x,
        "ratio_y": ratio_y,
        "ratio_opt": math.sqrt(ratio_x * ratio_y),
        "sensitivity": sensitivity,
        "fraction_ax": ax,
        "fraction_bx": bx,
        "fraction_ay": ay,
        "fraction_by": by,
        "a": a,
        "b": b,
    }
