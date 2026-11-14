"""Minimal, self-contained MIST-style subformula assignment (numpy only).

Re-implementation of the logic used by

    samgoldman97/mist, branch main_v2
    src/mist/subformulae/assign_subformulae.py   (MIT licence, (c) Samuel Goldman et al.)
      -> src/mist/utils/spectra_utils.py :: process_spec_file, max_inten_spec, assign_subforms, get_output_dict
      -> src/mist/utils/chem_utils.py    :: get_all_subsets, rdbe_filter, clipped_ppm, ion masses

Reference copies of those three files are kept (read-only, for audit) in ~/tools/frigid/reference/.

It produces the same JSON that FRIGID's PeakFormula featurizer reads:

    {"cand_form": "C8H10N4O2", "cand_ion": "[M+H]+",
     "output_tbl": {"mz": [...], "ms2_inten": [...], "mono_mass": [...], "abs_mass_diff": [...],
                    "mass_diff": [...], "formula": [...], "ions": [...]}}      # output_tbl is null if no peak matched

Parameters follow the task brief (and the script's defaults): 20 ppm (relative to max(m/z, 200), as in MIST's
`clipped_ppm`), keep the top 50 peaks by (sqrt-normalised) intensity, intensity threshold 0.001.

Faithfulness notes (deliberately kept so features match what the encoder was trained on):
  * Peaks are merged at 4-decimal m/z (max intensity wins), peaks above precursor+1 are dropped, intensities are
    divided by the maximum and then square-rooted; the 0.001 threshold is applied to the SQRT intensity.
  * Candidate sub-formulae are all element-count sub-vectors of the neutral formula with RDBE >= 0.
  * The ion mass is added to each sub-formula mass; the nearest sub-formula within the ppm window is assigned.
  * Several peaks assigned to the same formula are merged exactly as MIST does (including its chained-sum quirk).
  * Upstream quirk: the original script exposes `--max-formulae 50` but its .ms-directory branch calls
    `utils.process_spec_file(meta, tuples)` without passing it, so the upstream effective cap is 60. Here the
    cap is a parameter (default 50, as requested). FRIGID then keeps only the 10 most intense assigned peaks
    (max_peaks=10), so this has no practical influence.
"""
from __future__ import annotations

import re
from functools import reduce
from typing import Iterable, Sequence

import numpy as np

# --- element table (same 18 elements / order as MIST's VALID_ELEMENTS) -------------------------------------------
VALID_ELEMENTS = ["C", "H", "As", "B", "Br", "Cl", "Co", "F", "Fe", "I", "K", "N", "Na", "O", "P", "S", "Se", "Si"]
# most-common-isotope masses, copied from RDKit's periodic table (what MIST's P_TBL.GetMostCommonIsotopeMass returns),
# so sub-formula masses are bit-identical to those used when the encoder features were computed.
_MONO = {"C": 12.0, "H": 1.007825032, "As": 74.9215965, "B": 11.0093054, "Br": 78.9183371, "Cl": 34.96885268,
         "Co": 58.933195, "F": 18.99840322, "Fe": 55.9349375, "I": 126.904473, "K": 38.96370668, "N": 14.003074,
         "Na": 22.98976928, "O": 15.99491462, "P": 30.97376163, "S": 31.972071, "Se": 79.9165213, "Si": 27.97692653}
VALID_MONO_MASSES = np.array([_MONO[e] for e in VALID_ELEMENTS])
_EL_IDX = {e: i for i, e in enumerate(VALID_ELEMENTS)}
ELEMENT_VECTORS = np.eye(len(VALID_ELEMENTS))
ELECTRON_MASS = 0.00054858
_FORMULA_RE = re.compile(r"([A-Z][a-z]*)([0-9]*)")

# RDBE multipliers: C +2, N +1, P +1, H/Cl/Br/I/F -1
_RDBE = np.zeros(len(VALID_ELEMENTS))
for _e, _w in zip(["C", "N", "P", "H", "Cl", "Br", "I", "F"], [2, 1, 1, -1, -1, -1, -1, -1]):
    _RDBE[_EL_IDX[_e]] = _w

ION_LST = ["[M+H]+", "[M+Na]+", "[M+K]+", "[M-H2O+H]+", "[M+H3N+H]+", "[M]+", "[M-H4O2+H]+"]
ION_REMAP = {i: i for i in ION_LST}
ION_REMAP.update({
    "[M+NH4]+": "[M+H3N+H]+", "M+H": "[M+H]+", "M+Na": "[M+Na]+", "M+H-H2O": "[M-H2O+H]+",
    "M-H2O+H": "[M-H2O+H]+", "M+NH4": "[M+H3N+H]+", "M-2H2O+H": "[M-H4O2+H]+", "[M-2H2O+H]+": "[M-H4O2+H]+",
})
_m = _MONO
ION_TO_MASS = {
    "[M+H]+": _m["H"] - ELECTRON_MASS,
    "[M+Na]+": _m["Na"] - ELECTRON_MASS,
    "[M+K]+": _m["K"] - ELECTRON_MASS,
    "[M-H2O+H]+": -_m["O"] - _m["H"] - ELECTRON_MASS,
    "[M+H3N+H]+": _m["N"] + _m["H"] * 4 - ELECTRON_MASS,
    "[M]+": 0 - ELECTRON_MASS,
    "[M-H4O2+H]+": -_m["O"] * 2 - _m["H"] * 3 - ELECTRON_MASS,
}
MAX_SUBSETS = 5_000_000  # safety cap on the sub-formula lattice (memory guard)


def standardize_adduct(adduct: str) -> str:
    a = adduct.replace(" ", "")
    a = ION_REMAP.get(a, a)
    if a not in ION_LST:
        raise ValueError(f"adduct {adduct!r} is not one of the positive adducts FRIGID/MIST was trained on: {ION_LST}")
    return a


def formula_to_dense(formula: str) -> np.ndarray:
    """'C8H10N4O2' -> count vector over VALID_ELEMENTS. Raises ValueError for unsupported elements / charges."""
    formula = formula.strip()
    consumed = "".join(s + n for s, n in _FORMULA_RE.findall(formula))
    if not formula or consumed != formula:
        raise ValueError(f"cannot parse formula {formula!r} (expected a neutral molecular formula such as C8H10N4O2)")
    vec = np.zeros(len(VALID_ELEMENTS))
    for sym, num in _FORMULA_RE.findall(formula):
        if sym not in _EL_IDX:
            raise ValueError(f"element {sym!r} is not supported by the MIST/FRIGID featurizer {VALID_ELEMENTS}")
        vec[_EL_IDX[sym]] += int(num) if num else 1
    return vec


def vec_to_formula(vec: np.ndarray) -> str:
    out = ""
    for i in np.argwhere(vec > 0).flatten():
        ct = int(vec[i])
        out += f"{VALID_ELEMENTS[i]}{ct}" if ct > 1 else VALID_ELEMENTS[i]
    return out


def _cross_sum(x, y):
    return (np.expand_dims(x, 0) + np.expand_dims(y, 1)).reshape(-1, y.shape[-1])


def get_all_subsets(formula: str):
    """All sub-formula count vectors (RDBE >= 0) of `formula` and their monoisotopic masses (MIST get_all_subsets)."""
    dense = formula_to_dense(formula)
    n_subsets = int(np.prod(dense + 1))
    if n_subsets > MAX_SUBSETS:
        raise ValueError(f"formula {formula} has {n_subsets:,} sub-formulae (> {MAX_SUBSETS:,}); too large for assignment")
    parts = []
    for ind in np.argwhere(dense > 0).flatten():
        parts.append(ELEMENT_VECTORS[ind] * np.arange(0, dense[ind] + 1).reshape(-1, 1))
    cross = reduce(_cross_sum, parts, np.zeros((1, ELEMENT_VECTORS.shape[-1])))
    keep = np.argwhere(1 + 0.5 * cross.dot(_RDBE) >= 0).flatten()  # rdbe_filter
    cross = cross[keep]
    return cross, cross.dot(VALID_MONO_MASSES)


def clipped_ppm(mass_diff: np.ndarray, mz: np.ndarray) -> np.ndarray:
    d = mz * 1.0
    d[mz < 200] = 200
    return mass_diff / d * 1e6


def max_inten_spec(spec: np.ndarray, max_num_inten: int = 60, inten_thresh: float = 0.0) -> np.ndarray:
    masses, intens = spec[:, 0], spec[:, 1]
    order = np.argsort(intens)[::-1]
    if max_num_inten is not None:
        order = order[:max_num_inten]
    masses, intens = masses[order], intens[order]
    mask = intens > inten_thresh
    return np.vstack([masses[mask], intens[mask]]).transpose(1, 0)


def process_spectrum(peaks: Iterable[Sequence[float]], precursor_mz: float, precision: int = 4,
                     max_peaks: int = 50, inten_thresh: float = 0.001):
    """MIST process_spec_file: merge at 4 decimals, drop > precursor+1, /max, sqrt, top-N, threshold."""
    merged: dict[float, list[float]] = {}
    order: list[float] = []
    for mz, inten in peaks:
        key = float(np.round(mz, precision))
        if key not in merged:
            merged[key] = [float(mz), float(inten)]
            order.append(key)
        elif inten > merged[key][1]:
            merged[key][1] = float(inten)
    if not merged:
        return None
    spec = np.array([merged[k] for k in order], dtype=float)
    spec = spec[spec[:, 0] <= precursor_mz + 1]
    if spec.shape[0] == 0 or spec[:, 1].max() <= 0:
        return None
    spec[:, 1] = np.sqrt(spec[:, 1] / spec[:, 1].max())
    return max_inten_spec(spec, max_num_inten=max_peaks, inten_thresh=inten_thresh)


def assign_subformulae(formula: str, adduct: str, peaks: Iterable[Sequence[float]], precursor_mz: float,
                       ppm: float = 20.0, max_peaks: int = 50, inten_thresh: float = 0.001) -> dict:
    """Return the MIST/FRIGID subformula dict for one spectrum."""
    ion = standardize_adduct(adduct)
    formula_to_dense(formula)  # validates
    out = {"cand_form": formula, "cand_ion": ion, "output_tbl": None}
    spec = process_spectrum(peaks, precursor_mz, max_peaks=max_peaks, inten_thresh=inten_thresh)
    if spec is None or spec.shape[0] == 0:
        return out

    cross, masses = get_all_subsets(formula)
    spec_masses, spec_intens = spec[:, 0].copy(), spec[:, 1].copy()
    masses_with_ion = masses + ION_TO_MASS[ion]
    mass_diffs = np.abs(spec_masses[:, None] - masses_with_ion[None, :])
    formula_inds = mass_diffs.argmin(-1)
    min_diff = mass_diffs[np.arange(len(mass_diffs)), formula_inds]
    rel = clipped_ppm(min_diff, spec_masses)

    valid = rel < ppm
    spec_masses, spec_intens = spec_masses[valid], spec_intens[valid]
    min_diff, rel, formula_inds = min_diff[valid], rel[valid], formula_inds[valid]
    forms = np.array([vec_to_formula(v) for v in cross[formula_inds]])
    form_masses = masses_with_ion[formula_inds]

    # merge peaks sharing a formula (MIST's exact, chained, procedure)
    idx_of: dict[str, int] = {}
    uniq = []
    for idx, f in enumerate(forms):
        uniq.append(f not in idx_of)
        gi = idx_of.get(f)
        if gi is not None:
            spec_intens[gi] += spec_intens[idx]
        idx_of[f] = idx
    uniq = np.array(uniq, dtype=bool)
    if uniq.sum() == 0:
        return out
    out["output_tbl"] = {
        "mz": [float(x) for x in spec_masses[uniq]],
        "ms2_inten": [float(x) for x in spec_intens[uniq]],
        "mono_mass": [float(x) for x in form_masses[uniq]],
        "abs_mass_diff": [float(x) for x in min_diff[uniq]],
        "mass_diff": [float(x) for x in rel[uniq]],
        "formula": [str(x) for x in forms[uniq]],
        "ions": [ion] * int(uniq.sum()),
    }
    return out
