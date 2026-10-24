"""MS1 小类：理论 m/z（分子式 + 加合物）、BUDDY 分子式推断。"""
import re

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from mzdesk.api import common
from mzdesk.api.common import InputError, fmt, signed
from mzdesk.shared import colors
from mzdesk.shared.spectrum import parse_spectrum
from mzdesk.spec import ADDUCTS_NEG, ADDUCTS_POS

INSTRUMENTS = {"Orbitrap": "orbitrap", "Q-TOF": "qtof", "FT-ICR": "fticr"}
ELEMENTS_NOTE = "Supported elements: C H N O P S F Cl Br I Na K."
ADDUCT_TABLE_PPM = 100      # the all-adducts table shows Δ ppm only this close to the observed m/z


def _neutral_formula(struct, warnings):
    """Formula of the structure / formula input, without a charge sign."""
    if struct["kind"] == "formula":
        return struct["formula"]
    mol = Chem.MolFromSmiles(struct["smiles"])
    if len(Chem.GetMolFrags(mol)) > 1:
        warnings.append("The structure has several disconnected parts; the formula covers all of them.")
    charge = Chem.GetFormalCharge(mol)
    if charge:
        warnings.append(f"The structure carries a net charge ({charge:+d}). Its formula is used as is: "
                        "pick [M]+ if the structure is the observed ion.")
    return re.sub(r"[+-]\d*$", "", rdMolDescriptors.CalcMolFormula(mol))


def _ion(bd, formula, adduct):
    """(theoretical m/z, ion formula with charge sign) or (None, None) when the adduct does not apply."""
    if adduct in ("[M]+", "[M]-"):      # msbuddy reads [M]+ as [M+H]+, so the radical / permanent ion is done here
        arr = bd._safe_read_formula(formula)
        if arr is None:
            return None, None
        charge = 1 if adduct == "[M]+" else -1
        return bd._mz_from_counts(arr, charge), bd.form_arr_to_str(np.asarray(arr, dtype=np.int16)) + adduct[-1]
    counts, charge = bd.ion_composition(formula, adduct)
    if counts is None:
        return None, None
    return bd._mz_from_counts(counts, charge), bd.form_arr_to_str(counts.astype(np.int16)) + ("+" if charge > 0 else "-")


def run_theory(inputs):
    from mzdesk.ms1 import msbuddy_engine as bd

    warnings = []
    struct = common.structure(inputs, label="Compound or formula", allow_formula=True)
    observed, positive, adduct = common.precursor(inputs, label="Observed precursor m/z")
    formula = _neutral_formula(struct, warnings)
    if not bd.formula_is_valid(formula):
        raise InputError(f"{formula}: cannot compute m/z. {ELEMENTS_NOTE}")

    rows, ions = [], {}
    for a in (ADDUCTS_POS if positive else ADDUCTS_NEG)[1:]:
        mz, ion_formula = _ion(bd, formula, a)
        ions[a] = (mz, ion_formula)
        ppm = bd.ppm_error(observed, mz) if (observed and mz) else None
        near = ppm is not None and abs(ppm) <= ADDUCT_TABLE_PPM
        rows.append([a, ion_formula or common.MISSING, fmt(mz, 4), signed(ppm, 1) if near else common.MISSING])

    if adduct is None:       # Auto: the adduct closest to the observed m/z, else [M+H]+ / [M-H]-
        default = "[M+H]+" if positive else "[M-H]-"
        if observed:
            fits = [(abs(bd.ppm_error(observed, mz)), a) for a, (mz, _) in ions.items() if mz]
            adduct = min(fits)[1] if fits else default
        else:
            adduct = default
        note = f"{adduct} (Auto)"
    else:
        note = adduct
    mz, ion_formula = ions.get(adduct) or _ion(bd, formula, adduct)
    if mz is None:
        raise InputError(f"{formula} cannot form {adduct} (the adduct removes atoms the formula does not have).")
    ppm = bd.ppm_error(observed, mz) if observed else None
    if ppm is not None and abs(ppm) > 10:
        warnings.append(f"Observed m/z {observed:.4f} is {ppm:+.1f} ppm from {formula} {adduct} ({mz:.4f}). "
                        "Check the structure, adduct and ion mode.")
    return {"outputs": {
        "theory": {"values": {"Precursor m/z": {"value": fmt(mz, 4), "note": note}, "Ion formula": ion_formula,
                              "Error vs observed": signed(ppm, 1) if ppm is not None else common.MISSING}},
        "adducts": {"rows": rows},
    }, "warnings": warnings}


def _formula_mismatch(bd, formula, precursor_mz, positive, adduct, tol_ppm):
    """Warning text when a user-given formula does not fit the precursor (from the legacy MS1 tab)."""
    theoretical = bd.ion_mz(formula, adduct)
    if theoretical is None:
        return None
    error = bd.ppm_error(precursor_mz, theoretical)
    if abs(error) <= tol_ppm:
        return None
    others = bd.ADDUCTS_POSITIVE if positive else bd.ADDUCTS_NEGATIVE
    fits = [a for a in others if a != adduct and (m := bd.ion_mz(formula, a)) is not None
            and abs(bd.ppm_error(precursor_mz, m)) <= tol_ppm]
    msg = (f"{formula} as {adduct}: theoretical m/z {theoretical:.4f}, measured {precursor_mz:.4f} "
           f"({error:+,.1f} ppm; tolerance {tol_ppm:g} ppm).")
    return msg + (f" {fits[0]} fits." if fits else " Check formula, adduct and precursor.")


def run_formula(inputs):
    from mzdesk.ms1 import msbuddy_engine as bd

    warnings = []
    precursor_mz, positive, adduct = common.precursor(inputs, required=True)
    ms1_text = common.text(inputs, "ms1")
    ms1 = parse_spectrum(ms1_text)
    if ms1_text and not ms1:
        warnings.append("No valid MS1 peaks; the isotope pattern is ignored.")
    ms2, file_precursor = common.spectrum(inputs, "ms2", "MS2 spectrum", warnings, required=False)
    instrument = INSTRUMENTS.get(inputs.get("instrument") or "Orbitrap", "orbitrap")
    ms1_tol, ms2_tol = bd.INSTRUMENT_TOLS[instrument]
    halogen = bool(inputs.get("halogen"))
    formula = common.text(inputs, "formula")
    if formula and not bd.formula_is_valid(formula):
        raise InputError(f"Invalid formula {formula}. {ELEMENTS_NOTE}")
    lo, hi = bd.PRECURSOR_MZ_RANGE
    if not formula and not lo < precursor_mz < hi:
        raise InputError(f"BUDDY supports precursor m/z between {lo:.0f} and {hi:.0f}.")
    if file_precursor and abs(file_precursor - precursor_mz) > 0.01:
        warnings.append(f"The MS2 file says precursor m/z {file_precursor:.4f}; the form value {precursor_mz:.4f} is used.")
    if formula:
        mismatch = _formula_mismatch(bd, formula, precursor_mz, positive, adduct or bd.default_adduct(positive), ms1_tol)
        if mismatch:
            warnings.append(mismatch)
    if ms1:
        m0_mz, m0_int = min(ms1, key=lambda p: abs(p[0] - precursor_mz))
        err = bd.ppm_error(m0_mz, precursor_mz)
        if abs(err) > ms1_tol:
            warnings.append(f"MS1 has no peak within {ms1_tol:g} ppm of the precursor (closest {m0_mz:.4f}, "
                            f"{err:+,.1f} ppm); BUDDY ignores the isotope pattern.")
        elif not halogen and not formula:
            m2 = [i for mz, i in ms1 if 1.98 <= mz - m0_mz <= 2.02]
            if m2 and max(m2) >= 0.2 * m0_int:
                warnings.append(f"M+2 is {max(m2) / m0_int:.0%} of M: tick Include halogens if Cl / Br is possible.")

    try:
        res = bd.annotate_formula(precursor_mz, positive, adduct, None, formula, tuple(ms2), tuple(ms1), instrument, halogen)
    except ValueError as exc:
        raise InputError(str(exc)) from None
    cands, used = res["candidates"], res["adduct_used"]
    if not cands:
        warnings.append("BUDDY returned no candidate. Check the precursor m/z, adduct and instrument.")
        return {"outputs": {"summary": {"values": {"Top formula": common.MISSING, "Adduct used": used,
                                                   "MS1 isotope peaks used": common.MISSING}}}, "warnings": warnings}
    top = cands[0]
    ms1_used = res.get("ms1_used")
    outputs = {
        "summary": {"values": {
            "Top formula": {"value": top["formula"], "note": "input formula, ranking skipped" if res["from_input"] else f"rank 1 of {len(cands)}"},
            "Adduct used": used,
            "MS1 isotope peaks used": common.MISSING if ms1_used is None else str(ms1_used)}},
        "candidates": {"rows": [[str(c["rank"]), c["formula"], fmt(c["mz_theo"], 4), signed(c["ppm"], 2),
                                 fmt(c["prob"], 4), fmt(c["fdr"], 4), fmt(c["isotope_sim"], 4)] for c in cands]},
    }
    if ms2:
        explained = bd.explain_fragments(tuple(mz for mz, _ in ms2), top["formula"], used, ms2_tol)
        rows = bd.build_annotation_rows(ms2, explained, precursor_mz)
        fig = bd.build_figure(rows, precursor_mz, adduct=used, dark=colors.PLOT_DARK, light=colors.PLOT_LIGHT,
                              ink=colors.INK, muted=colors.MUTED, grid=colors.SURFACE_2, highlight=colors.BLUE_BRIGHT)
        outputs["fragments_plot"] = common.plotly(fig)
        outputs["fragments"] = {"rows": [[fmt(r["mz"], 4), fmt(r["rel"], 1),
                                          (r["formula"] + (" (precursor)" if r["precursor"] else "")) or common.MISSING,
                                          signed(r["ppm"], 1) if r["annotated"] else common.MISSING,
                                          fmt(r["nl"], 4), ", ".join(r["alternatives"]) or common.MISSING] for r in rows]}
    return {"outputs": outputs, "warnings": warnings}


RUNNERS = {"ms1_theory": run_theory, "ms1_formula": run_formula}
