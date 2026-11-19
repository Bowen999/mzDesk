"""Calculator engines: ppm error, isotope labels, concentration, collision energy and Kovats RI."""
from __future__ import annotations

from dataclasses import dataclass
import math
import re

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from mzdesk.ms1 import isotopes
from mzdesk.platform_fit.derivatization import neutralize
from mzdesk.platform_fit.structure_prep import largest_fragment, parts_note


ELECTRON = isotopes._ELECTRON
MOLAR_UNITS = {"µM": 1e-6, "mM": 1e-3}
MASS_UNITS = {"ng/mL": 1e-6, "µg/mL": 1e-3, "mg/L": 1e-3}


def ppm_error(observed: float, theoretical: float) -> float:
    if theoretical <= 0:
        raise ValueError("Theoretical m/z must be positive.")
    return (observed - theoretical) / theoretical * 1e6


def concentration(value: float, from_unit: str, to_unit: str, molar_mass: float | None = None) -> float:
    if from_unit not in MOLAR_UNITS | MASS_UNITS or to_unit not in MOLAR_UNITS | MASS_UNITS:
        raise ValueError("Unknown concentration unit.")
    if from_unit in MOLAR_UNITS and to_unit in MOLAR_UNITS:
        return value * MOLAR_UNITS[from_unit] / MOLAR_UNITS[to_unit]
    if from_unit in MASS_UNITS and to_unit in MASS_UNITS:
        return value * MASS_UNITS[from_unit] / MASS_UNITS[to_unit]
    if molar_mass is None or molar_mass <= 0:
        raise ValueError("Molar mass is required to convert between mass and molar concentration.")
    if from_unit in MOLAR_UNITS:
        grams_per_l = value * MOLAR_UNITS[from_unit] * molar_mass
        return grams_per_l / MASS_UNITS[to_unit]
    moles_per_l = value * MASS_UNITS[from_unit] / molar_mass
    return moles_per_l / MOLAR_UNITS[to_unit]


def dilution(stock: float, stock_unit: str, target: float, target_unit: str, final_volume_ml: float,
             molar_mass: float | None = None) -> dict:
    if final_volume_ml <= 0:
        raise ValueError("Final volume must be positive.")
    stock_in_target = concentration(stock, stock_unit, target_unit, molar_mass)
    if target > stock_in_target:
        raise ValueError("Target concentration cannot exceed the stock concentration.")
    stock_volume = target * final_volume_ml / stock_in_target
    return {
        "stock_in_target": stock_in_target,
        "target": target,
        "final_volume": final_volume_ml,
        "stock_volume": stock_volume,
        "diluent_volume": final_volume_ml - stock_volume,
    }


def nce_to_ev(nce: float, precursor_mz: float, charge: int) -> float:
    if precursor_mz <= 0 or charge <= 0:
        raise ValueError("Precursor m/z and charge must be positive.")
    return nce * precursor_mz * charge / 500.0


def ev_to_nce(ev: float, precursor_mz: float, charge: int) -> float:
    if precursor_mz <= 0 or charge <= 0:
        raise ValueError("Precursor m/z and charge must be positive.")
    return ev * 500.0 / (precursor_mz * charge)


def parse_alkanes(text: str) -> list[tuple[int, float]]:
    points = []
    for number, line in enumerate((text or "").splitlines(), start=1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        tokens = [token for token in re.split(r"[\s,;]+", line) if token]
        if len(tokens) < 2:
            continue
        carbon = tokens[0].strip().upper().removeprefix("C")
        try:
            n = int(carbon)
            rt = float(tokens[1])
        except ValueError:
            if number == 1:
                continue
            raise ValueError(f"Alkane line {number}: use a carbon number and a retention time.") from None
        if n <= 0 or not math.isfinite(rt):
            raise ValueError(f"Alkane line {number}: carbon number and RT must be positive finite numbers.")
        points.append((n, rt))
    if len(points) < 2:
        raise ValueError("Enter at least two n-alkanes (carbon number, RT).")
    return points


def kovats_ri(alkanes: list[tuple[int, float]], analyte_rt: float) -> dict:
    if len(alkanes) < 2:
        raise ValueError("At least two n-alkanes are required.")
    points = sorted(alkanes)
    if len({n for n, _ in points}) != len(points):
        raise ValueError("Each carbon number must appear once.")
    for (n1, rt1), (n2, rt2) in zip(points, points[1:]):
        if n2 <= n1:
            raise ValueError("Carbon numbers must increase.")
        if rt2 <= rt1:
            raise ValueError("Alkane retention times must increase with carbon number.")
    pair = None
    extrapolated = False
    for left, right in zip(points, points[1:]):
        if left[1] <= analyte_rt <= right[1]:
            pair = (left, right)
            break
    if pair is None:
        extrapolated = True
        pair = (points[0], points[1]) if analyte_rt < points[0][1] else (points[-2], points[-1])
    (n1, rt1), (n2, rt2) = pair
    ri = 100.0 * n1 + 100.0 * (analyte_rt - rt1) / (rt2 - rt1) * (n2 - n1)
    return {"ri": ri, "lower": (n1, rt1), "upper": (n2, rt2), "extrapolated": extrapolated}


@dataclass(frozen=True)
class IsotopeLabel:
    text: str
    element: str
    mass_number: int
    count: int
    isotope_mass: float
    mono_mass: float

    @property
    def exact_shift(self) -> float:
        return (self.isotope_mass - self.mono_mass) * self.count

    @property
    def nominal_shift(self) -> int:
        return (self.mass_number - round(self.mono_mass)) * self.count


def _isotope_mass(element: str, mass_number: int) -> float:
    from IsoSpecPy import PeriodicTbl
    masses = PeriodicTbl.symbol_to_masses.get(element)
    probs = PeriodicTbl.symbol_to_probs.get(element)
    if not masses:
        raise ValueError(f"Unknown element {element!r} in isotope label.")
    candidates = [mass for mass in masses if round(mass) == mass_number]
    if not candidates:
        raise ValueError(f"{element}-{mass_number} is not a known isotope.")
    return max(candidates, key=lambda mass: probs[masses.index(mass)])


def parse_labels(text: str) -> list[IsotopeLabel]:
    tokens = [token for token in re.split(r"[\s,;]+", (text or "").strip()) if token]
    if not tokens:
        raise ValueError("Enter an isotope label such as 13C6, 15N2 or d4.")
    labels = []
    for token in tokens:
        match = re.fullmatch(r"[dD](\d*)", token)
        if match:
            element, mass_number, count = "H", 2, int(match.group(1) or 1)
        else:
            match = re.fullmatch(r"(\d+)([A-Z][a-z]?)(\d*)", token)
            if not match:
                raise ValueError(f"Cannot read isotope label {token!r}; use 13C6, 15N2 or d4.")
            mass_number, element, count = int(match.group(1)), match.group(2), int(match.group(3) or 1)
        if count <= 0:
            raise ValueError(f"Isotope label {token!r} must have a positive count.")
        mono_mass = isotopes._mono_mass(element)
        labels.append(IsotopeLabel(token, element, mass_number, count, _isotope_mass(element, mass_number), mono_mass))
    return labels


def _neutral_formula_counts(formula: str) -> dict[str, int]:
    counts, _ = isotopes.parse_formula(re.sub(r"[+-]\d*$", "", formula))
    return counts


def validate_labels(labels: list[IsotopeLabel], formula: str) -> None:
    counts = _neutral_formula_counts(formula)
    wanted: dict[str, int] = {}
    for label in labels:
        wanted[label.element] = wanted.get(label.element, 0) + label.count
    for element, count in wanted.items():
        if count > counts.get(element, 0):
            raise ValueError(f"The label needs {count} {element} atom(s), but the molecule has {counts.get(element, 0)}.")


def smiles_formula_mass(smiles: str, warnings: list | None = None) -> tuple[str, float]:
    """Neutral formula and monoisotopic mass of a SMILES. A salt or mixture is taken by its largest fragment (the other pages' choice; the note goes
    to ``warnings`` when given), and a charged fragment (benzoate, protonated amine) by its neutral form. A net-neutral zwitterion (betaine) is M as drawn."""
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        raise ValueError("The SMILES could not be parsed.")
    mol, n_parts = largest_fragment(mol)
    if n_parts > 1 and warnings is not None:
        warnings.append(parts_note(n_parts))
    if Chem.GetFormalCharge(mol):
        mol, _ = neutralize(mol)
    formula = re.sub(r"[+-]\d*$", "", rdMolDescriptors.CalcMolFormula(mol))
    return formula, float(Descriptors.ExactMolWt(mol))


def formula_mass(formula: str) -> float:
    counts, _ = isotopes.parse_formula(re.sub(r"[+-]\d*$", "", formula))
    return sum(isotopes._mono_mass(element) * count for element, count in counts.items())


def adduct_charge(adduct: str) -> int:
    match = re.search(r"(\d*)([+-])$", (adduct or "").strip())
    if not match:
        raise ValueError(f"Unknown adduct {adduct}.")
    magnitude = int(match.group(1) or 1)
    return magnitude if match.group(2) == "+" else -magnitude


def adduct_mz(formula: str, adduct: str) -> float:
    counts, charge = isotopes.ion_composition(formula, adduct)
    mass = sum(isotopes._mono_mass(element) * count for element, count in counts.items())
    return (mass - charge * ELECTRON) / abs(charge)


def _distribution(formula: str, adduct: str, labels: list[IsotopeLabel] | None = None) -> dict:
    counts, charge = isotopes.ion_composition(formula, adduct)
    extra_mass = 0.0
    remaining = dict(counts)
    for label in labels or ():
        remaining[label.element] = remaining.get(label.element, 0) - label.count
        if remaining[label.element] < 0:
            raise ValueError(f"The label needs more {label.element} atoms than the ion contains.")
        extra_mass += label.isotope_mass * label.count
    remaining = {element: count for element, count in remaining.items() if count}
    base_mass = sum(isotopes._mono_mass(element) * count for element, count in remaining.items())
    if remaining:
        import IsoSpecPy
        iso = IsoSpecPy.IsoTotalProb(formula=isotopes._formula_string(remaining), prob_to_cover=0.999999, get_confs=False)
        masses, probs = iso.masses, iso.probs
    else:
        masses, probs = [base_mass], [1.0]
    groups: dict[int, list[float]] = {}
    for mass, probability in zip(masses, probs):
        offset = round(mass - base_mass)
        group = groups.setdefault(offset, [0.0, 0.0])
        group[0] += probability
        group[1] += probability * mass
    total = sum(group[0] for group in groups.values())
    peaks = []
    for offset, (probability, mass_sum) in sorted(groups.items()):
        mass = mass_sum / probability
        peaks.append({
            "offset": offset,
            "mz": (mass + extra_mass - charge * ELECTRON) / abs(charge),
            "fraction": probability / total,
        })
    top = max(peak["fraction"] for peak in peaks)
    for peak in peaks:
        peak["abundance"] = peak["fraction"] / top * 100.0
    return {"peaks": peaks, "charge": charge, "base_mass": base_mass + extra_mass}


def labelled_distribution(formula: str, adduct: str, labels: list[IsotopeLabel]) -> dict:
    validate_labels(labels, formula)
    result = _distribution(formula, adduct, labels)
    result["label_text"] = ", ".join(label.text for label in labels)
    result["exact_shift"] = sum(label.exact_shift for label in labels)
    result["nominal_shift"] = sum(label.nominal_shift for label in labels)
    result["mono_mz"] = (result["base_mass"] - result["charge"] * ELECTRON) / abs(result["charge"])
    return result


def crosstalk(formula: str, adduct: str, nominal_shift: int) -> float:
    result = _distribution(formula, adduct)
    by_offset = {peak["offset"]: peak["fraction"] for peak in result["peaks"]}
    base = by_offset.get(0)
    if not base:
        raise ValueError("The natural M peak is not available.")
    return by_offset.get(nominal_shift, 0.0) / base * 100.0
