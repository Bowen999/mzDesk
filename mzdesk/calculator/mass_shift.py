"""Mass-shift search: parent ion, feature-table parsing and fixed biotransformation shifts."""
from __future__ import annotations

import math
import re

from mzdesk.calculator import calculations
from mzdesk.ms1 import isotopes


_H = isotopes._mono_mass("H")
_C = isotopes._mono_mass("C")
_N = isotopes._mono_mass("N")
_O = isotopes._mono_mass("O")
_S = isotopes._mono_mass("S")

SHIFT_SPECS = (
    ("+O (hydroxylation)", _O, "Hydroxylation"),
    ("−2H (dehydrogenation)", -2 * _H, "Dehydrogenation"),
    ("+2H (reduction)", 2 * _H, "Reduction"),
    ("−CH2 (demethylation)", -(_C + 2 * _H), "Demethylation"),
    ("+O −CH2", _O - (_C + 2 * _H), "Oxidative demethylation"),
    ("+C6H8O6 glucuronide", 6 * _C + 8 * _H + 6 * _O, "Glucuronidation"),
    ("+SO3 sulfate", _S + 3 * _O, "Sulfation"),
    ("+C2H2O acetyl", 2 * _C + 2 * _H + _O, "Acetylation"),
    ("+C10H15N3O6S GSH", 10 * _C + 15 * _H + 3 * _N + 6 * _O + _S, "Glutathione conjugation"),
    ("+C2H3NO2 glycine", 2 * _C + 3 * _H + _N + 2 * _O, "Glycine conjugation"),
)
SHIFTS = {label: delta for label, delta, _ in SHIFT_SPECS}
REACTIONS = {label: reaction for label, _, reaction in SHIFT_SPECS}

_MZ_KEYS = {"m/z", "mz", "mass", "observedm/z", "featurem/z", "precursormz", "precursorm/z", "precursor"}
_RT_KEYS = {"rt", "retentiontime", "retention", "time", "tr"}
_INT_KEYS = {"intensity", "int", "area", "height", "abundance", "response", "signal"}
MAX_FEATURES = 20000


def _key(token: str) -> str:
    text = token.strip().strip('"').strip("'").lower()
    text = re.sub(r"[\[(][^)\]]*[)\]]", "", text)
    return re.sub(r"[\s_]+", "", text)


def _tokens(line: str) -> list[str]:
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


def _columns(tokens: list[str]) -> tuple[int, int | None, int | None] | None:
    keys = [_key(token) for token in tokens]
    mz = next((i for i, key in enumerate(keys) if key in _MZ_KEYS), None)
    if mz is None:
        return None
    rt = next((i for i, key in enumerate(keys) if key in _RT_KEYS), None)
    intensity = next((i for i, key in enumerate(keys) if key in _INT_KEYS), None)
    return mz, rt, intensity


def parse_features(text: str) -> tuple[list[dict], list[str]]:
    lines = [line.strip() for line in (text or "").lstrip("\ufeff").splitlines()]
    lines = [line for line in lines if line and not line.startswith("#")]
    if not lines:
        return [], []
    warnings = []
    columns = _columns(_tokens(lines[0]))
    start = 1 if columns else 0
    two_column = False
    features = []
    skipped = 0
    for number, line in enumerate(lines[start:], start=start + 1):
        tokens = _tokens(line)
        if not tokens:
            continue
        if columns:
            mz_i, rt_i, int_i = columns
            if mz_i >= len(tokens):
                skipped += 1
                continue
            mz = _number(tokens[mz_i])
            rt = _number(tokens[rt_i]) if rt_i is not None and rt_i < len(tokens) else None
            intensity = _number(tokens[int_i]) if int_i is not None and int_i < len(tokens) else None
        else:
            mz = _number(tokens[0])
            if len(tokens) == 1:
                rt, intensity = None, 1.0
            elif len(tokens) == 2:
                rt, intensity = None, _number(tokens[1])
                two_column = True
            else:
                rt, intensity = _number(tokens[1]), _number(tokens[2])
        if mz is None or mz <= 0 or (rt is not None and not math.isfinite(rt)) or (intensity is not None and intensity < 0):
            skipped += 1
            continue
        features.append({"mz": mz, "rt": rt, "intensity": intensity, "line": number})
    if len(features) > MAX_FEATURES:
        raise ValueError(f"The feature list has {len(features)} rows; the limit is {MAX_FEATURES}.")
    if skipped:
        warnings.append(f"{skipped} feature line(s) were skipped because they were not numeric or had m/z ≤ 0.")
    if two_column:
        warnings.append("Two-column features were read as m/z and intensity; RT is left blank.")
    if not features:
        raise ValueError("No valid features found. Use one row per line: m/z [RT] [intensity].")
    return features, warnings


def parent_identity(struct: dict) -> dict:
    if struct["kind"] == "formula":
        formula = struct["formula"]
        return {"formula": formula, "mass": calculations.formula_mass(formula), "smiles": None}
    formula, mass = calculations.smiles_formula_mass(struct["smiles"])
    return {"formula": formula, "mass": mass, "smiles": struct["smiles"]}


def _canonical(label: str) -> str:
    if label in SHIFTS:
        return label
    alternate = label.replace("-", "−")
    if alternate in SHIFTS:
        return alternate
    raise ValueError(f"Unknown mass shift: {label}.")


def search(features: list[dict], base_mz: float, charge: int, selected: list[str], tol_ppm: float) -> dict:
    if not selected:
        raise ValueError("Choose at least one mass shift.")
    selected = list(dict.fromkeys(_canonical(label) for label in selected))
    if charge == 0:
        raise ValueError("The adduct charge cannot be zero.")
    if tol_ppm < 0:
        raise ValueError("The match tolerance cannot be negative.")
    rows, matched, matched_features = [], 0, set()
    for label in selected:
        delta = SHIFTS[label]
        expected = base_mz + delta / abs(charge)
        if expected <= 0:
            raise ValueError(f"The expected m/z for {label} is not positive.")
        matches = []
        for index, feature in enumerate(features):
            ppm = (feature["mz"] - expected) / expected * 1e6
            if abs(ppm) <= tol_ppm:
                matches.append((abs(ppm), ppm, index, feature))
        matches.sort(key=lambda item: item[0])
        if matches:
            matched += 1
        for _, ppm, index, feature in matches:
            matched_features.add(index)
            rows.append({"label": label, "delta": delta, "expected": expected, "feature": feature, "ppm": ppm,
                         "reaction": REACTIONS[label]})
        if not matches:
            rows.append({"label": label, "delta": delta, "expected": expected, "feature": None, "ppm": None,
                         "reaction": REACTIONS[label]})
    return {"rows": rows, "matched_shifts": matched, "matched_features": matched_features}
