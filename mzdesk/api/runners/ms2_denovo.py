"""De novo structure 页：FRIGID-base（本机 PyTorch）从一张正离子 MS/MS 谱 + 分子式生成候选结构（假设，不是鉴定）。

分子式没填就先用 BUDDY（和 Formula annotation 页同一个引擎）取第一名，结果里写明用了哪个分子式、从哪来。
引擎见 mzdesk/ms2/frigid_engine.py；没装 FRIGID 时返回一条说明（502）。
"""
from mzdesk.api import common
from mzdesk.api.common import EngineError, InputError, fmt

IMG_SIZE = (200, 140)
SMILES_IN_CARD = 44         # longer SMILES are cut in the card; the whole string is in the card's "?" and in the table
MISMATCH_PPM = 20.0         # FRIGID's own sub-formula assignment window
DEVICE_NAMES = {"mps": "MPS (Apple GPU)", "cpu": "CPU"}


def _cut(text, n=SMILES_IN_CARD):
    return text if len(text) <= n else text[:n - 1] + "…"


def _formula_mismatch(fe, formula, precursor_mz, adduct):
    """Warning text when the formula as the chosen adduct does not land on the precursor m/z (within FRIGID's 20 ppm window)."""
    from rdkit import Chem
    table = Chem.GetPeriodicTable()
    mass = sum(n * table.GetMostCommonIsotopeMass(el) for el, n in fe.formula_counts(formula).items())
    theoretical = fe.ion_mz(mass, adduct)
    ppm = (precursor_mz - theoretical) / theoretical * 1e6
    if abs(ppm) <= MISMATCH_PPM:
        return None
    return (f"{formula} as {adduct} has m/z {theoretical:.4f}, but the precursor is {precursor_mz:.4f} ({ppm:+,.1f} ppm; FRIGID assigns fragments within "
            f"{MISMATCH_PPM:g} ppm). Check the formula and the adduct.")


def _formula_from_buddy(fe, precursor_mz, adduct, ms2, warnings):
    """(formula, source note) from BUDDY's top-1 formula; InputError when BUDDY cannot give one."""
    from mzdesk.ms1 import msbuddy_engine as bd

    if adduct not in bd.ADDUCTS_POSITIVE:
        raise InputError(f"No formula was entered, and BUDDY cannot infer one for {adduct} (it handles {', '.join(bd.ADDUCTS_POSITIVE)}). "
                         "Enter the molecular formula.")
    try:
        res = bd.annotate_formula(precursor_mz, True, adduct, None, "", tuple(ms2), (), "orbitrap", False)
    except ValueError as exc:
        raise InputError(f"No formula was entered, and BUDDY could not be used: {exc} Enter the molecular formula.") from None
    cands = res["candidates"]
    if not cands:
        raise InputError(f"No formula was entered, and BUDDY found none for m/z {precursor_mz:.4f} as {adduct}. "
                         "Check the precursor m/z and the adduct, or enter the molecular formula.")
    top = cands[0]
    try:
        formula = fe.check_formula(top["formula"])
    except ValueError as exc:
        raise InputError(f"BUDDY's top formula {top['formula']} cannot be used by FRIGID: {exc} Enter the molecular formula.") from None
    prob = top.get("prob")
    note = f"BUDDY top-1 of {len(cands)}" + (f" · probability {prob:.2f}" if isinstance(prob, (int, float)) else "")
    others = [c["formula"] for c in cands[1:3]]
    warnings.append(f"No formula was entered, so BUDDY's top-1 formula {formula}" + (f" (probability {prob:.2f})" if isinstance(prob, (int, float)) else "")
                    + " was used (Orbitrap settings, no halogens). If it is wrong, every structure below is wrong too"
                    + (f"; BUDDY's next candidates: {', '.join(others)}" if others else "") + ". Enter the formula to override it.")
    return formula, note


def _rows(candidates):
    from mzdesk.ms2.metfrag_engine import render_frag_png

    return [[str(c["rank"]), common.png(render_frag_png(c["smiles"], IMG_SIZE)), c["smiles"], c["formula"] or common.MISSING, fmt(c["score"], 3),
             (c["inchikey_2d"] or common.MISSING)] for c in candidates]


def _dropped_text(dropped):
    parts = [(dropped["isotope"], "isotope-labelled"), (dropped["radical"], "radical"), (dropped["charged"], "wrong net charge"),
             (dropped["duplicate"], "duplicate")]
    return ", ".join(f"{n} {label}" for n, label in parts if n)


def run_denovo(inputs):
    from mzdesk.ms2 import frigid_engine as fe

    warnings = []
    field_mz, positive, adduct = common.precursor(inputs)
    if not positive:
        raise InputError(fe.NEGATIVE_MESSAGE)
    ms2, file_mz = common.spectrum(inputs, "ms2", "MS2 spectrum", warnings)
    precursor_mz = field_mz or file_mz
    if precursor_mz is None:
        raise InputError("Precursor m/z is required. Enter it, or use a spectrum whose MGF / MSP header has it.")
    if field_mz is None:
        warnings.append(f"Precursor m/z {file_mz:.4f} is read from the spectrum.")
    elif file_mz and abs(file_mz - field_mz) > 0.01:
        warnings.append(f"The MS2 file says precursor m/z {file_mz:.4f}; the form value {field_mz:.4f} is used.")
    try:
        adduct, _ = fe.resolve_adduct(adduct)
        mode = fe.resolve_mode(common.text(inputs, "mode"))
    except ValueError as exc:
        raise InputError(str(exc)) from None

    typed = common.text(inputs, "formula")
    if typed:
        if not common.is_formula(typed):
            raise InputError(f"Invalid formula {typed}. Use a plain neutral formula such as C8H10N4O2.")
        try:
            formula = fe.check_formula(typed)
        except ValueError as exc:
            raise InputError(str(exc)) from None
        source = "entered by you"
        mismatch = _formula_mismatch(fe, formula, precursor_mz, adduct)
        if mismatch:
            warnings.append(mismatch)
    else:
        formula, source = _formula_from_buddy(fe, precursor_mz, adduct, ms2, warnings)

    try:
        res = fe.predict(formula, precursor_mz, adduct, ms2, mode)
    except fe.FrigidError as exc:
        raise EngineError(str(exc), exc.log) from None
    except ValueError as exc:
        raise InputError(str(exc)) from None

    cands, dropped = res["candidates"], res["dropped"]
    # the driver's "no generated molecule matches the formula" is said better below; its other notes (few assigned peaks, ...) pass through
    warnings.extend(w[:1].upper() + w[1:] for w in res["warnings"] if not w.startswith("no generated molecule matches"))
    if res["n_peaks_assigned"] == 0:
        warnings.append("No peak could be assigned a fragment formula of this molecular formula: the formula or the adduct is probably wrong.")
    hidden = _dropped_text(dropped)
    if not cands:
        warnings.append(f"FRIGID generated no usable structure for {formula} in {res['n_samples']} samples" + (f" ({hidden} hidden)" if hidden else "") + ". "
                        + ("Try Full mode or check the formula; " if res["mode"] == "quick" else "Check the formula; ")
                        + "the compound may also be outside what FRIGID can generate.")
    elif hidden:
        warnings.append(f"Hidden: {hidden} (not valid candidates for this formula and adduct).")
    if cands and res["mode"] == "quick" and len(cands) < 3:
        warnings.append(f"Only {len(cands)} formula-matched structure{'s' if len(cands) != 1 else ''} in {res['n_samples']} samples; Full mode samples up to 100.")

    top = cands[0] if cands else None
    gb = (res["peak_footprint_mb"] or 0) / 1024
    device = DEVICE_NAMES.get(res["device"], str(res["device"] or common.MISSING).upper())
    device_help = ["Device auto picks the Apple GPU (MPS) only when the Mac has at least 12 GB of memory available and swap to spare; otherwise the CPU. "
                   "MZDESK_FRIGID_DEVICE = auto | mps | cpu overrides it."]
    if res.get("device_reason"):
        device_help.insert(0, "This run: " + res["device_reason"] + ".")
    sample_help = [f"{res['n_samples']} molecules sampled, {res['n_valid']} valid SMILES, {res['n_formula_matched']} with the right formula "
                   f"({res['n_unique_matched']} different structures)." + (f" Hidden: {hidden}." if hidden else "")]
    summary = {
        "Top candidate": ({"value": _cut(top["smiles"]), "note": f"score {fmt(top['score'], 3)} · {top['formula']}", "help": ["SMILES: " + top["smiles"]]}
                          if top else common.MISSING),
        "Formula used": {"value": formula, "note": source},
        "Device": {"value": device, "note": f"{res['seconds']:.0f} s" + (f" · peak memory {gb:.1f} GB" if gb else "") + f" · {res['device_requested']}",
                   "help": device_help},
        "Samples": {"value": f"{res['n_samples']} / {len(cands)}", "note": f"generated / shown · {res['mode']} mode", "help": sample_help}}
    return {"outputs": {
        "summary": {"values": summary},
        "candidates": {"rows": _rows(cands), "count": len(cands)}}, "warnings": warnings}


RUNNERS = {"ms2_denovo": run_denovo}
