"""Pace et al. (2011) 的 spICP-MS 信号 → 粒径 / 数量浓度换算，不处理原始时间序列。"""
import math
import statistics

from mzdesk.inorganic.core import ELEMENT_SYMBOLS, formula_average_mass, formula_counts, require_finite
from mzdesk.inorganic.data import PARTICLE_PRESETS, SENSITIVITY_DWELL, SENSITIVITY_RATE


def mass_fraction(formula, measured_element):
    counts = formula_counts(formula)
    if measured_element not in ELEMENT_SYMBOLS:
        raise ValueError(f"Unknown measured element {measured_element!r}.")
    if measured_element not in counts:
        raise ValueError(f"Measured element {measured_element} is not present in {formula}.")
    averages_mass = formula_average_mass(counts)
    element_mass = formula_average_mass({measured_element: counts[measured_element]})
    return element_mass / averages_mass


def resolve_material(preset, custom_formula=None, custom_density=None, measured_element=None):
    if preset == "Custom":
        if not custom_formula:
            raise ValueError("Custom particle composition is required.")
        if custom_density is None or float(custom_density) <= 0:
            raise ValueError("Custom particle density must be positive.")
        formula, density, name = custom_formula, float(custom_density), "Custom"
    elif preset in PARTICLE_PRESETS:
        formula, density = PARTICLE_PRESETS[preset]
        name = preset
    else:
        raise ValueError(f"Unknown particle material {preset!r}.")
    fraction = mass_fraction(formula, measured_element)
    return {"name": name, "formula": formula, "density": density, "mass_fraction": fraction}


def _positive(value, label):
    value = require_finite(float(value), label)
    if value <= 0:
        raise ValueError(f"{label} must be positive.")
    return value


def _event(intensity, background, sensitivity_dwell, eta, flow_ml_min, dwell_ms, density, fraction):
    net = intensity - background
    if net <= 0:
        raise ValueError("Every event intensity and the optional size-LOD threshold must exceed the background.")
    dwell_s = dwell_ms / 1000
    flow_l_s = flow_ml_min / 60000
    element_mass_ug = net * eta * flow_l_s * dwell_s / sensitivity_dwell
    particle_mass_g = element_mass_ug * 1e-6 / fraction
    diameter_nm = (6 * particle_mass_g / (math.pi * density)) ** (1 / 3) * 1e7
    return {"net_intensity": net, "element_mass_ug": element_mass_ug,
            "particle_mass_g": particle_mass_g, "diameter_nm": diameter_nm}


def calculate(material, sensitivity, sensitivity_unit, dwell_ms, flow_ml_min, transport_pct,
              background, intensities, n_events, acquisition_s, threshold=None):
    sensitivity = _positive(sensitivity, "Dissolved calibration sensitivity S")
    dwell_ms = _positive(dwell_ms, "Dwell time")
    flow_ml_min = _positive(flow_ml_min, "Sample flow q")
    acquisition_s = _positive(acquisition_s, "Acquisition time")
    transport_pct = _positive(transport_pct, "Transport efficiency")
    if transport_pct > 100:
        raise ValueError("Transport efficiency cannot exceed 100%.")
    background = require_finite(float(background), "Background counts per dwell")
    if background < 0:
        raise ValueError("Background counts per dwell cannot be negative.")
    if sensitivity_unit == SENSITIVITY_DWELL:
        sensitivity_dwell = sensitivity
    elif sensitivity_unit == SENSITIVITY_RATE:
        sensitivity_dwell = sensitivity * dwell_ms / 1000
    else:
        raise ValueError("Sensitivity unit must be counts per dwell or counts/s per µg/L.")
    if not intensities:
        raise ValueError("Provide at least one particle event intensity.")
    events = []
    for index, intensity in enumerate(intensities, 1):
        intensity = require_finite(float(intensity), f"Event intensity {index}")
        try:
            event = _event(intensity, background, sensitivity_dwell, transport_pct / 100,
                           flow_ml_min, dwell_ms, material["density"], material["mass_fraction"])
        except ValueError:
            raise ValueError(f"Event {index} intensity must exceed the background.") from None
        events.append({"event": index, "intensity": intensity, **event})
    n_events = require_finite(float(n_events), "Number of events N")
    if n_events < len(events) or n_events != int(n_events):
        raise ValueError("N must be an integer greater than or equal to the number of listed event intensities.")
    n_events = int(n_events)
    eta = transport_pct / 100
    number_ml = n_events / (eta * (flow_ml_min / 60) * acquisition_s)
    diameters = [event["diameter_nm"] for event in events]
    size_lod = None
    if threshold is not None:
        threshold = require_finite(float(threshold), "Threshold intensity")
        if threshold <= background:
            raise ValueError("Threshold intensity must exceed the background.")
        size_lod = _event(threshold, background, sensitivity_dwell, eta, flow_ml_min, dwell_ms,
                          material["density"], material["mass_fraction"])["diameter_nm"]
    return {
        "events": events,
        "mean_diameter_nm": statistics.fmean(diameters),
        "median_diameter_nm": statistics.median(diameters),
        "number_concentration_ml": number_ml,
        "number_concentration_l": number_ml * 1000,
        "size_lod_nm": size_lod,
        "sensitivity_dwell": sensitivity_dwell,
    }
