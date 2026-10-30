"""Step 5 of the ESI MS/MS Annotation workflow: FRIGID structure hypotheses for the top-1 formula.

run_frigid_layer() asks FRIGID (mzdesk/ms2/frigid_engine.py: a MIST encoder predicts a fingerprint from the spectrum, a diffusion language
model writes SMILES for the formula and that fingerprint) for structures of the top-1 formula of step 2. Positive ions only, so a negative-mode
adduct or one outside fe.ADDUCTS ends the layer as "skipped" with the reason. Every candidate is a hypothesis (Candidate.hypothesis, Level H):
the formula is the only evidence, the spectrum is used through a predicted fingerprint, not matched to the structure.

The score is FRIGID's own: the Tanimoto similarity of the generated structure's fingerprint to the predicted one, 0-1 and independent of the
other candidates. layer_score_norm = score_norm(score): linear, 0 at 0.5 and 1 at 1.0, clipped. Calibration: 18 real Quick-mode runs with the
true formula on 18 spectra of a seeded random sample of positive-mode MassBank 2026.03 records (14 [M+H]+, 3 [M+Na]+, 1 [M+NH4]+; no CASMI
data) and the caffeine run of tests/fixtures/frigid. Two runs gave no structure of the formula at all (a 220-peak spectrum, a sodium citrate); the
other 17, by the score of their top candidate:
      raw score      norm          runs    top-1 is the true structure
      < 0.50         0.00             3        0
      0.50 - 0.70    0.00 - 0.40      5        1   (aldicarb [M+Na]+, 0.68)
      0.70 - 0.85    0.40 - 0.70      2        1   (0.74 right, 0.79 wrong)
      >= 0.85        >= 0.70          7        7   (0.88 - 1.00)
In 9 of the 17 the true structure was generated and ranked first; in the other 8 it was not generated (best scores 0.36 - 0.79). A small sample of
known, mostly drug-like compounds: read the table as "a score of 0.85 or more has always been the right structure so far", not as a probability.
The re-rank multiplies a pure FRIGID hypothesis by 0.7 on top of this.

Time and memory: 15-35 s per run on an idle machine (17 s on the Apple GPU), up to 80 s on a busy one (CPU), about 5 GB; one heavy job at a time
(heavy lock, shared with MS-FINDER). All real engine rules (formula check, fingerprint, device choice) are the engine's, not repeated here.
"""
from __future__ import annotations

import math
import re
import time

from mzdesk import adducts
from mzdesk.annotation.msfinder_step import clean_peaks
from mzdesk.annotation.schema import LAYER_FRIGID, STATUS_DONE, STATUS_FAILED, STATUS_SKIPPED, Candidate, LayerResult
from mzdesk.ms2 import frigid_engine as fe

# ---------------------------------------------------------------- constants
TOP_N = 5                           # owner's design: top-5 structures
SCORE_LOW, SCORE_HIGH = 0.5, 1.0    # the FRIGID scores where layer_score_norm is 0 and 1 (calibration table above)
MODE_LABELS = {"quick": "Quick", "full": "Full"}


def score_norm(score):
    """FRIGID score (Tanimoto to the predicted fingerprint) -> 0-1: linear from 0 at SCORE_LOW to 1 at SCORE_HIGH, clipped (absolute)."""
    if score is None or not math.isfinite(score):
        return 0.0
    return min(1.0, max(0.0, (score - SCORE_LOW) / (SCORE_HIGH - SCORE_LOW)))


def _ppm(formula, adduct, precursor_mz):
    """(observed - theoretical) / theoretical x 1e6 of the precursor against `formula` as `adduct`, with the engine's own arithmetic."""
    from rdkit import Chem

    table = Chem.GetPeriodicTable()
    theo = fe.ion_mz(sum(n * table.GetMostCommonIsotopeMass(el) for el, n in fe.formula_counts(formula).items()), adduct)
    return (precursor_mz - theo) / theo * 1e6


# ---------------------------------------------------------------- the step
def run_frigid_layer(precursor_mz, adduct, ms2_peaks, formula, *, top_n=TOP_N, mode="quick"):
    """Step 5. Returns a LayerResult; never raises (a missing or failing engine gives status "failed" and a plain note).

    precursor_mz  measured precursor m/z
    adduct        any spelling mzdesk.adducts knows; sent to FRIGID as adducts.engine_name(adduct, "frigid") (positive ions only)
    ms2_peaks     [(m/z, intensity)]
    formula       neutral molecular formula M, the top-1 formula of step 2 (None / empty: skipped)
    top_n         structures kept (one per 2D structure, best score first)
    mode          "quick" (16 samples, 15-50 s) or "full" (up to 100 samples)

    status "skipped": no formula, no MS2 peaks, a formula FRIGID cannot use, or an adduct it does not support (negative ions; anything outside
    fe.ADDUCTS); "failed": FRIGID is missing, busy or fails; else "done" (also with no candidates: FRIGID found no structure of the formula).
    """
    started = time.perf_counter()

    def result(status, note, candidates=()):
        return LayerResult(LAYER_FRIGID, status, note, time.perf_counter() - started, list(candidates))

    try:
        return _run(result, precursor_mz, adduct, ms2_peaks, formula, top_n, mode)
    except fe.FrigidError as exc:
        return result(STATUS_FAILED, str(exc))
    except Exception as exc:                # noqa: BLE001 - a layer must never take the workflow down
        return result(STATUS_FAILED, f"FRIGID step failed: {type(exc).__name__}: {exc}")


def _run(result, precursor_mz, adduct, ms2_peaks, formula, top_n, mode):
    canonical, frigid_adduct = adducts.canonical(adduct), adducts.engine_name(adduct, "frigid")
    if not frigid_adduct:
        if adducts.polarity(adduct) == "Negative":
            return result(STATUS_SKIPPED, f"FRIGID supports positive ions only; {canonical} is a negative-ion adduct.")
        return result(STATUS_SKIPPED, f"FRIGID does not support the adduct {canonical or adduct or '(none given)'}; it takes {', '.join(fe.ADDUCTS)}.")
    formula = (formula or "").strip()
    if not formula:
        return result(STATUS_SKIPPED, "There is no molecular formula to start from (step 2 gave none), so FRIGID was not run.")
    try:
        formula = fe.check_formula(formula)
    except ValueError as exc:
        return result(STATUS_SKIPPED, f"FRIGID cannot use this formula: {exc}")
    peaks = clean_peaks(ms2_peaks, fe.MAX_PEAKS)
    if not peaks:
        return result(STATUS_SKIPPED, "There are no MS2 peaks for FRIGID to work from.")
    try:
        res = fe.predict(formula, precursor_mz, frigid_adduct, peaks, mode)
    except ValueError as exc:               # input FRIGID rejects after all
        return result(STATUS_SKIPPED, str(exc))

    ppm = _ppm(formula, canonical, precursor_mz)
    candidates = []
    for c in res["candidates"][:max(1, int(top_n))]:
        key = c["inchikey_2d"] if re.fullmatch(r"[A-Z]{14}", c.get("inchikey_2d") or "") else None
        candidates.append(Candidate(
            id=f"frigid_{len(candidates) + 1}", layer=LAYER_FRIGID, smiles=c["smiles"], ik14=key, name="", formula=c.get("formula") or formula,
            adduct=canonical, mz_error_ppm=ppm, layer_score=c["score"], layer_score_norm=score_norm(c["score"]),
            evidence={LAYER_FRIGID: {"score": c["score"], "rank": c["rank"], "mode": res["mode"], "n_samples": res["n_samples"],
                                     "n_formula_matched": res["n_formula_matched"], "n_unique_matched": res["n_unique_matched"], "device": res["device"]}},
            hypothesis=True))
    label = f"{MODE_LABELS.get(res['mode'], res['mode'])} mode, {res['n_samples']} samples"
    if not candidates:
        return result(STATUS_DONE, f"FRIGID generated no structure with the formula {formula} ({label}).")
    best, total = candidates[0], len(res["candidates"])
    count = f"{len(candidates)} structure{'' if total == 1 else 's'}" if len(candidates) == total else f"{len(candidates)} of {total} structures"
    return result(STATUS_DONE, f"{count} for {formula} ({label}); best fingerprint similarity {best.layer_score:.2f}, "
                               f"normalised {best.layer_score_norm:.2f}: hypotheses, not identifications.", candidates)
