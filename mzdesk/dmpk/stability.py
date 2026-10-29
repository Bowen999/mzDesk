"""Metabolic stability: ln-linear fit of an in vitro time course, CLint and scaling to hepatic clearance.

Pure calculations only: the runner parses the pasted table and formats the outputs, this module owns the
numbers and the assumptions behind them.

Constants (one source per column; every number also appears in the page's parameter table):

* MPPGL = 45 mg microsomal protein per g liver and HPGL = 99 x 10^6 hepatocytes per g liver are the human
  consensus values of Barter et al. 2007 (Curr Drug Metab 8:33-45; the hepatocyte data reviewed there go
  back to Wilson et al. 2003).  They are applied to every species here: species-specific microsomal
  recovery / hepatocyte yield is not modelled, so CLint_scaled for a non-human species differs only
  through liver weight and hepatic blood flow.
* Liver weight (g liver per kg body weight) and hepatic blood flow Qh (mL/min/kg body weight) are from
  Davies & Morris 1993 (Pharm Res 10:1093-1095).  The human values are 26 g/kg (2.6 % of body weight) and
  20.7 mL/min/kg.

S9 deliberately reuses the microsomal MPPGL: S9 is not microsomal protein, so this is a simplification and
the page shows it as an explicit assumption.
"""
from __future__ import annotations

import math
import re

SYSTEMS = ("Liver microsomes", "Hepatocytes", "S9", "Plasma")

# Species constants.  mppgl / hpgl: Barter et al. 2007 (human consensus, applied to all species, see the
# module docstring).  liver_g_per_kg / qh: Davies & Morris 1993.
SPECIES = {
    "Human": {"mppgl": 45.0, "hpgl": 99.0, "liver_g_per_kg": 26.0, "qh": 20.7},
    "Rat": {"mppgl": 45.0, "hpgl": 99.0, "liver_g_per_kg": 40.0, "qh": 55.2},
    "Mouse": {"mppgl": 45.0, "hpgl": 99.0, "liver_g_per_kg": 50.0, "qh": 90.0},
    "Dog": {"mppgl": 45.0, "hpgl": 99.0, "liver_g_per_kg": 32.0, "qh": 30.9},
    "Monkey": {"mppgl": 45.0, "hpgl": 99.0, "liver_g_per_kg": 26.0, "qh": 43.6},
}

FU = 1.0   # no fu input on the page: plasma protein binding is assumed away

_LN2 = math.log(2.0)

# Header detection: "Time (min) / % remaining" in its common spellings.
_TIME_KEYS = {"time", "time(min)", "timemin", "t", "t(min)", "min", "minutes"}
_PCT_KEYS = {"%remaining", "%remaining(%)", "percentremaining", "percentremaining(%)", "remaining",
             "remaining(%)", "%", "percent", "pct"}


def _tokens(line: str) -> list[str]:
    """Split on tab / comma / semicolon, else on whitespace; quotes and empty cells are dropped."""
    if "\t" in line:
        parts = line.split("\t")
    elif "," in line:
        parts = line.split(",")
    elif ";" in line:
        parts = line.split(";")
    else:
        parts = re.split(r"\s+", line)
    return [part.strip().strip('"').strip("'") for part in parts if part.strip()]


def _number(token: str) -> float | None:
    try:
        value = float(token)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _looks_like_header(tokens: list[str]) -> bool:
    """True for a text header.  A line with a numeric token is never a header: it is a data row to check."""
    if not tokens:
        return False
    keys = {re.sub(r"[\s_\[\]]+", "", token.lower()) for token in tokens}
    if keys & _TIME_KEYS and keys & _PCT_KEYS:
        return True
    return all(_number(token) is None for token in tokens)


def parse_timecourse(text: str) -> list[tuple[float, float]]:
    """Pasted table -> [(time_min, % remaining), ...].  Raises ValueError naming the offending line.

    One optional header line, comma / tab / semicolon / whitespace separators, blank lines and # comments
    are accepted.  Exactly two numeric columns, at least three rows, strictly increasing time and
    % remaining in (0, 100].
    """
    rows: list[tuple[float, float]] = []
    previous_time = None
    first = True
    for lineno, raw in enumerate((text or "").lstrip("\ufeff").splitlines(), 1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        tokens = _tokens(line)
        if first:
            first = False
            if _looks_like_header(tokens):
                continue
        if len(tokens) != 2:
            raise ValueError(f"Line {lineno}: expected two columns (time and % remaining), found {len(tokens)}.")
        time, percent = _number(tokens[0]), _number(tokens[1])
        if time is None or percent is None:
            raise ValueError(f"Line {lineno}: time and % remaining must both be numbers, got '{' '.join(tokens)}'.")
        if time < 0:
            raise ValueError(f"Line {lineno}: time must not be negative (got {time:g} min).")
        if previous_time is not None and time <= previous_time:
            raise ValueError(f"Line {lineno}: time must increase; {time:g} min is not after the previous {previous_time:g} min.")
        if not 0 < percent <= 100:
            raise ValueError(f"Line {lineno}: % remaining must be in (0, 100], got {percent:g}.")
        rows.append((time, percent))
        previous_time = time
    if len(rows) < 3:
        raise ValueError(f"At least three time points are needed; got {len(rows)}.")
    return rows


def fit_curve(points: list[tuple[float, float]]) -> dict:
    """Unweighted OLS of ln(% remaining / 100) vs time.  Raises ValueError if the slope is not negative."""
    if len(points) < 3:
        raise ValueError(f"At least three time points are needed for the fit; got {len(points)}.")
    xs = [float(time) for time, _ in points]
    ys = [math.log(percent / 100.0) for _, percent in points]
    n = len(xs)
    x_mean = sum(xs) / n
    y_mean = sum(ys) / n
    sxx = sum((x - x_mean) ** 2 for x in xs)
    sxy = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys))
    if sxx <= 0:
        raise ValueError("Time must increase so that a slope can be fitted.")
    slope = sxy / sxx
    if slope >= 0:
        raise ValueError(f"The fitted slope is {slope:.4g} min⁻¹, not negative: % remaining does not fall with "
                         "time, so no half-life can be derived.")
    intercept = y_mean - slope * x_mean
    ss_res = sum((y - (intercept + slope * x)) ** 2 for x, y in zip(xs, ys))
    ss_tot = sum((y - y_mean) ** 2 for y in ys)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    k = -slope
    return {"k": k, "slope": slope, "intercept": intercept, "r2": r2, "n": n, "t_half": _LN2 / k}


def clint_unit(system: str) -> str | None:
    """Unit of the in vitro CLint for a system (None for plasma, which has no CLint)."""
    if system == "Hepatocytes":
        return "µL/min/10⁶ cells"
    if system == "Plasma":
        return None
    return "µL/min/mg protein"


def clint_invitro(k: float, system: str, protein_or_cells: float | None) -> float | None:
    """k (min⁻¹) -> CLint per mg protein or per 10⁶ cells; None for plasma or a non-positive amount."""
    if system not in SYSTEMS:
        raise ValueError(f"Unknown system '{system}'. Choose one of: {', '.join(SYSTEMS)}.")
    if system == "Plasma":
        return None
    if protein_or_cells is None or not math.isfinite(protein_or_cells) or protein_or_cells <= 0:
        return None
    return k * 1000.0 / protein_or_cells


def scale_clint(clint: float | None, system: str, species: str) -> float | None:
    """In vitro CLint -> mL/min/kg body weight.

    µL/min/mg × mg/g liver × g liver/kg body weight / 1000 µL/mL = mL/min/kg; the hepatocyte path uses
    10⁶ cells/g liver instead of mg/g, which leaves the unit arithmetic unchanged.
    """
    if clint is None or system == "Plasma":
        return None
    constants = SPECIES[species]
    per_gram = constants["hpgl"] if system == "Hepatocytes" else constants["mppgl"]
    return clint * per_gram * constants["liver_g_per_kg"] / 1000.0


def hepatic_clearance(clint_scaled: float | None, species: str, fu: float = FU) -> tuple[float | None, float | None]:
    """Well-stirred hepatic clearance and extraction ratio; (None, None) when there is no scaled CLint."""
    if clint_scaled is None:
        return None, None
    qh = SPECIES[species]["qh"]
    clh = qh * fu * clint_scaled / (qh + fu * clint_scaled)
    return clh, clh / qh


def compute(points: list[tuple[float, float]], system: str, protein_or_cells: float | None, species: str) -> dict:
    """Everything the page needs: fit, in vitro CLint, scaled CLint, well-stirred CLh and extraction ratio."""
    if system not in SYSTEMS:
        raise ValueError(f"Unknown system '{system}'. Choose one of: {', '.join(SYSTEMS)}.")
    if species not in SPECIES:
        raise ValueError(f"Unknown species '{species}'. Choose one of: {', '.join(SPECIES)}.")
    result = fit_curve(points)
    clint = clint_invitro(result["k"], system, protein_or_cells)
    scaled = scale_clint(clint, system, species)
    clh, er = hepatic_clearance(scaled, species)
    constants = SPECIES[species]
    result.update({
        "system": system,
        "species": species,
        "protein_or_cells": protein_or_cells,
        "clint": clint,
        "clint_unit": clint_unit(system),
        "clint_scaled": scaled,
        "clh": clh,
        "er": er,
        "fu": FU,
        "mppgl": constants["mppgl"],
        "hpgl": constants["hpgl"],
        "liver_g_per_kg": constants["liver_g_per_kg"],
        "liver_fraction": constants["liver_g_per_kg"] / 1000.0,
        "qh": constants["qh"],
        "last_time": points[-1][0],
        "last_percent": points[-1][1],
    })
    return result
