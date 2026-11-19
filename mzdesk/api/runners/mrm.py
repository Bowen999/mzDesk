"""MRM transitions 页：从预测碎片 / 库谱 / 自己的谱选离子对，给出角色、CE、碎片式和来源。

三条路径：
  结构（预测碎片）  CFM-ID 三个能量合并（m/z ±0.01 Da，取各碎片最强强度和出现能量），按强度取前 N 个不同 m/z；
                    前体 m/z 用 MetFrag 的理论值。CE 单位 eV 时用能量，NCE 时用通用换算。
  库谱              RDKit InChIKey 前 14 位 + 同一加合物 + 前体 ±10 ppm，取峰最多的一条 MassBank 记录（并列取 |Δppm| 最小）；
                    CE 按记录的 collision_energy 原样显示。
  自己的谱          从粘贴谱头部读前体 m/z，取最强的前 N 个峰；CE 未知，来源 My spectrum。
加合物不由页面输入：正离子固定 [M+H]+，负离子固定 [M-H]-。
"""
from mzdesk.api import blocks, common
from mzdesk.api.common import EngineError, InputError, MISSING, fmt
from mzdesk.ms2 import mrm as engine

START_FROM = ("Structure (predicted fragments)", "Library spectrum", "My spectrum")
CE_UNITS = ("eV", "NCE")
ADDUCT = {"Positive": "[M+H]+", "Negative": "[M-H]-"}
COLUMNS = ["Precursor m/z", "Product m/z", "CE", "Role", "Fragment formula", "Source", "Relative intensity"]
LIBRARY = "MassBank"
LIBRARY_TOL_PPM = 10.0
DEFAULT_TRANSITIONS = 3


# ---------------------------------------------------------------- inputs
def _start_from(inputs):
    value = common.text(inputs, "start_from") or START_FROM[0]
    if value not in START_FROM:
        raise InputError(f"Unknown starting point {value!r}. Choose one of: {', '.join(START_FROM)}.")
    return value


def _ion_mode(inputs):
    mode = common.text(inputs, "ion_mode") or "Positive"
    if mode not in ADDUCT:
        raise InputError("Ion mode must be Positive or Negative.")
    return mode, ADDUCT[mode], mode == "Positive"


def _transitions(inputs):
    value = common.number(inputs, "n_transitions", "Transitions", default=DEFAULT_TRANSITIONS, minimum=1)
    if value != int(value):
        raise InputError("Transitions must be a whole number.")
    return int(value)


def _ce_unit(inputs):
    unit = common.text(inputs, "ce_unit") or "eV"
    if unit not in CE_UNITS:
        raise InputError(f"CE unit must be {' or '.join(CE_UNITS)}.")
    return unit


def _structure_smiles(inputs, start_from):
    """Structure of a structure-requiring path; a path-specific message when it is missing."""
    raw = str((inputs.get("structure") or {}).get("value") or "").strip()
    if not raw:
        if start_from == "Library spectrum":
            raise InputError("The library path needs a structure: it identifies the compound by its InChIKey.")
        raise InputError("The predicted-fragments path needs a structure. Draw or enter the compound, or start from My spectrum.")
    return common.structure(inputs, label="Structure")["smiles"]


def _theoretical_precursor(smiles, adduct, warnings):
    """Theoretical precursor m/z of the neutral structure for the fixed adduct (MetFrag, no run)."""
    from mzdesk.ms2 import metfrag_engine as mf

    info = mf.precursor_info(smiles, adduct)
    if info.get("error"):
        raise InputError(info["error"])
    dropped = info.get("dropped") or []
    if dropped:
        warnings.append("The structure has several disconnected parts; only the largest is used. Ignored: "
                        + ", ".join(dropped) + ".")
    return info["theo_mz"]


# ---------------------------------------------------------------- shared result shaping
def _select(products, precursor_mz, n_transitions, warnings):
    """Engine ranking with the warnings the page promises (precursor peaks, too few products)."""
    try:
        selection = engine.select_transitions(products, precursor_mz, n_transitions)
    except ValueError as exc:
        raise InputError(str(exc)) from None
    skipped = selection["n_precursor_skipped"]
    if skipped:
        warnings.append(f"{skipped} peak{'s' if skipped != 1 else ''} at or above the precursor m/z left out: "
                        "the precursor ion is not a product.")
    found = len(selection["transitions"])
    if found < n_transitions:
        warnings.append(f"Fewer product ions below the precursor m/z than requested: {found} of "
                        f"{n_transitions} transitions available.")
    return selection["transitions"]


def _row(precursor_mz, transition, ce, formula, source):
    return [fmt(precursor_mz, 4), fmt(transition["mz"], 4), ce or MISSING, transition["role"].capitalize(),
            formula or MISSING, source, f"{transition['relative']:.1f} %"]


# ---------------------------------------------------------------- path 1: predicted fragments
def _predicted(warnings, smiles, precursor_mz, adduct, n_transitions, ce_unit):
    from mzdesk.ms2 import cfmid_engine as cf

    try:
        result = cf.run_cfmid(smiles, adduct)
    except cf.CFMIDError as exc:
        raise EngineError(str(exc), exc.log) from None
    except ValueError as exc:
        raise InputError(str(exc)) from None
    warnings += result["warnings"]
    try:
        products = engine.merge_energies(result["spectra"])
    except ValueError as exc:
        raise InputError(str(exc)) from None
    transitions = _select(products, precursor_mz, n_transitions, warnings)
    rows = [_row(precursor_mz, transition, engine.ce_label(transition["ev"], precursor_mz, ce_unit),
                 engine.fragment_formula(transition.get("annotation")), "CFM-ID") for transition in transitions]
    note = ("Rows come from CFM-ID predicted fragments: the 10 / 20 / 40 eV spectra are merged by m/z, "
            "and CE is the energy where the product ion is most intense.")
    return rows, note


# ---------------------------------------------------------------- path 2: library spectrum
def _library(warnings, smiles, precursor_mz, adduct, positive, n_transitions, ce_unit):
    from mzdesk.ms2 import library_index as lx

    if not lx.is_installed(LIBRARY):
        raise InputError("MassBank index is not installed. Run: python -m mzdesk.ms2.library_setup MassBank")
    try:
        key = engine.structure_key(smiles)
    except ValueError as exc:
        raise InputError(str(exc)) from None
    try:
        library = lx.load(LIBRARY)
    except lx.LibraryMissing as exc:
        raise InputError(str(exc)) from None
    ids, deltas = library.msms_candidates(precursor_mz, LIBRARY_TOL_PPM, positive, adduct=adduct)
    infos = library.info(ids)
    matches = [(int(row), float(delta), info) for row, delta, info in zip(ids, deltas, infos)
               if str(info.get("inchikey") or "")[:14] == key]
    others = len(infos) - len(matches)
    if others:
        warnings.append(f"{others} MassBank MS/MS spectra within {LIBRARY_TOL_PPM:g} ppm as {adduct} have a "
                        "different InChIKey and were ignored.")
    if not matches:
        raise InputError(f"No MassBank MS/MS spectrum of InChIKey {key} as {adduct} within "
                         f"{LIBRARY_TOL_PPM:g} ppm of m/z {precursor_mz:.4f}.")
    ranked = []
    for row, delta, info in matches:
        mz, _ = library.peaks(row)
        ranked.append((len(mz), -abs(delta), row, info))
    _, _, row, info = max(ranked, key=lambda item: (item[0], item[1]))
    if len(matches) > 1:
        warnings.append(f"{len(matches)} MassBank spectra match the compound and adduct; the one with the most "
                        f"peaks was used ({info['accession']}).")
    peaks_mz, peaks_intensity = library.peaks(row)
    products = [{"mz": float(mz), "intensity": float(intensity)} for mz, intensity in zip(peaks_mz, peaks_intensity)]
    transitions = _select(products, precursor_mz, n_transitions, warnings)
    ce = str(info.get("collision_energy") or "").strip() or None
    if ce is None:
        warnings.append(f"MassBank {info['accession']} has no recorded collision energy; CE is left blank.")
    if ce_unit == "NCE":
        warnings.append("MassBank collision energies are shown as recorded; the CE unit selection applies to "
                        "the predicted-fragment path only.")
    source = f"MassBank {info['accession']}"
    rows = [_row(precursor_mz, transition, ce, None, source) for transition in transitions]
    note = (f"Rows come from MassBank {info['accession']}, a recorded spectrum of the same compound and adduct. "
            f"Only MassBank spectra recorded as {adduct} are considered; CE is shown as recorded and not converted.")
    return rows, note


# ---------------------------------------------------------------- path 3: my spectrum
def _my_spectrum(inputs, warnings, n_transitions):
    if str((inputs.get("structure") or {}).get("value") or "").strip():
        warnings.append("The structure is ignored when starting from my spectrum; the rows come from the pasted spectrum.")
    raw = str((inputs.get("ms2") or {}).get("value") or "")
    if not raw.strip():
        raise InputError("Starting from my spectrum needs an MS2 spectrum: paste the peak list, with a PEPMASS "
                         "or BEGIN IONS header for the precursor m/z.")
    peaks, file_precursor = common.spectrum(inputs, "ms2", "MS2 spectrum", warnings)
    if file_precursor is None:
        raise InputError("The pasted spectrum has no precursor m/z. Add a PEPMASS or BEGIN IONS header line, "
                         "or start from the structure path instead.")
    warnings.append(f"Precursor m/z {file_precursor:.4f} is read from the spectrum header.")
    products = [{"mz": mz, "intensity": intensity} for mz, intensity in peaks]
    transitions = _select(products, file_precursor, n_transitions, warnings)
    rows = [_row(file_precursor, transition, None, None, "My spectrum") for transition in transitions]
    note = ("Rows come from the pasted MS2 spectrum; the precursor m/z is read from its header and the "
            "collision energy is unknown (CE is blank).")
    return file_precursor, rows, note


# ---------------------------------------------------------------- blocks
def _notes(path_note, predicted, instrument):
    ce_note = ("Collision energies are instrument-specific: transfer them to the instrument and optimize them there. "
               "This page does not adapt CE to the instrument.")
    items = [
        path_note,
        "Adduct: fixed to [M+H]+ for Positive and [M-H]- for Negative mode. This page has no adduct input; "
        "use the Adducts page for [M+Na]+, [M+Cl]- and other ions.",
        f"Instrument: {instrument}. {ce_note}" if instrument else ce_note,
        "The quantifier (most intense product) / qualifier choice is an intensity heuristic, not a validated "
        "choice; confirm the transitions with a standard and the sample matrix.",
        "Product ions are limited to peaks below the precursor m/z, so the precursor ion is never a product. "
        "Relative intensity is normalised to the base peak of the chosen spectrum (100).",
    ]
    if predicted:
        items.append("The predicted fragments are a CFM-ID model output, not a measured spectrum.")
        items.append("eV to NCE uses the generic conversion NCE = CE × 500 / (precursor m/z × z) with z = 1. "
                     "The Calculator page uses the same conversion; it is not instrument-validated.")
    return items


def run_mrm(inputs):
    warnings = []
    start_from = _start_from(inputs)
    mode, adduct, positive = _ion_mode(inputs)
    n_transitions = _transitions(inputs)
    ce_unit = _ce_unit(inputs)
    if start_from == "My spectrum":
        precursor_mz, rows, path_note = _my_spectrum(inputs, warnings, n_transitions)
    else:
        smiles = _structure_smiles(inputs, start_from)
        precursor_mz = _theoretical_precursor(smiles, adduct, warnings)
        if start_from == "Structure (predicted fragments)":
            rows, path_note = _predicted(warnings, smiles, precursor_mz, adduct, n_transitions, ce_unit)
        else:
            rows, path_note = _library(warnings, smiles, precursor_mz, adduct, positive, n_transitions, ce_unit)
    outputs = {
        "summary": {"values": {
            "Polarity": mode,
            "Transitions": str(len(rows))}},
        "transitions": blocks.table_data(rows, count=len(rows)),
        "method_table": {"formats": {"Method table": engine.method_table_csv(COLUMNS, rows)}},
        "notes": {"items": _notes(path_note, start_from == "Structure (predicted fragments)",
                                  common.text(inputs, "instrument"))},
    }
    return {"outputs": outputs, "warnings": warnings}


RUNNERS = {"mrm": run_mrm}
