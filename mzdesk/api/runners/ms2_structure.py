"""Structure prediction 页：MS-FINDER（本机 Mono）对一张 MS/MS 谱给出分子式候选和结构候选。

不联网：结构只来自 MS-FINDER 自带的本地结构库。引擎见 mzdesk/ms2/msfinder_engine.py；服务器没有 Mono 时返回一条说明（502）。
"""
from mzdesk.api import common
from mzdesk.api.common import EngineError, InputError, fmt, signed
from mzdesk.shared.spectrum import parse_spectrum

STRUCTURE_ROWS = 50         # rows shown (MS-FINDER reports up to 100 per formula)
FORMULA_ROWS = 50
DB_NAMES_SHOWN = 4
HELP_IDS = 40
IMG_SIZE = (200, 140)
BASE_ELEMENTS = "C H N O P S"


def _databases_cell(databases):
    """'HMDB, ChEBI, DrugBank, KNApSAcK +9' with the IDs in the cell's "?"."""
    if not databases:
        return common.MISSING
    names = [name for name, _ in databases]
    text = ", ".join(names[:DB_NAMES_SHOWN]) + (f" +{len(names) - DB_NAMES_SHOWN}" if len(names) > DB_NAMES_SHOWN else "")
    lines = [f"{name}: {ident}" if ident else name for name, ident in databases[:HELP_IDS]]
    if len(databases) > HELP_IDS:
        lines.append(f"… {len(databases) - HELP_IDS} more")
    return {"text": text, "help": lines}


def _structure_rows(structures):
    from mzdesk.ms2.metfrag_engine import render_frag_png

    rows = []
    for rank, s in enumerate(structures[:STRUCTURE_ROWS], 1):
        rows.append([str(rank), common.png(render_frag_png(s["smiles"], IMG_SIZE)), s["name"] or common.MISSING, s["smiles"] or common.MISSING,
                     s["formula"] or common.MISSING, fmt(s["score"], 4), s["ontology"] or common.MISSING, s["inchikey"] or common.MISSING,
                     _databases_cell(s["databases"])])
    return rows


def _formula_rows(formulas):
    return [[str(rank), f["formula"], fmt(f["neutral_mass"], 4),
             signed(None if f["mass_error_da"] is None else f["mass_error_da"] * 1000, 2), fmt(f["score"], 3)]
            for rank, f in enumerate(formulas[:FORMULA_ROWS], 1)]


def _elements(halogen, silicon):
    """Elements the formula search used, for the summary card."""
    return " ".join([BASE_ELEMENTS] + (["F Cl Br I"] if halogen else []) + (["Si"] if silicon else []))


def _element_hint(halogen, silicon):
    """Sentence that tells the user which element switch to try; empty when both are already on."""
    if halogen and silicon:
        return ""
    if not halogen and not silicon:
        return "Include halogens (F, Cl, Br, I) and Include Si are off; tick them if the compound may contain those elements."
    if not halogen:
        return "Include halogens is off; tick it if the compound may contain F, Cl, Br or I."
    return "Include Si is off; tick it if the compound may contain silicon."


def run_structure(inputs):
    from mzdesk.ms2 import msfinder_engine as mf

    warnings = []
    halogen, silicon = common.flag(inputs, "halogen", "Include halogens"), common.flag(inputs, "silicon", "Include Si")
    field_mz, positive, adduct = common.precursor(inputs)
    ms2, file_mz = common.spectrum(inputs, "ms2", "MS2 spectrum", warnings)
    precursor_mz = field_mz or file_mz
    if precursor_mz is None:
        raise InputError("Precursor m/z is required. Enter it, or use a spectrum whose MGF / MSP header has it.")
    if field_mz is None:
        warnings.append(f"Precursor m/z {file_mz:.4f} is read from the spectrum.")
    elif file_mz and abs(file_mz - field_mz) > 0.01:
        warnings.append(f"The MS2 file says precursor m/z {file_mz:.4f}; the form value {field_mz:.4f} is used.")

    ms1_text = common.text(inputs, "ms1")
    ms1 = parse_spectrum(ms1_text)
    if ms1_text and not ms1:
        warnings.append("No valid MS1 peaks; the isotope pattern is ignored.")
    if ms1:
        m0_mz, _ = min(ms1, key=lambda p: abs(p[0] - precursor_mz))
        if abs(m0_mz - precursor_mz) > mf.MS1_MATCH_DA:
            warnings.append(f"MS1 has no peak within {mf.MS1_MATCH_DA * 1000:g} mDa of the precursor m/z (closest {m0_mz:.4f}, "
                            f"{(m0_mz - precursor_mz) * 1000:+.1f} mDa); MS-FINDER cannot match the isotope pattern. Use the same m/z in both fields.")

    try:
        res = mf.predict(precursor_mz, positive, adduct, ms2, ms1, halogen=halogen, silicon=silicon)
    except mf.MsFinderError as exc:
        raise EngineError(str(exc), exc.log) from None
    except ValueError as exc:
        raise InputError(str(exc)) from None

    structures, formulas = res["structures"], res["formulas"]
    elements, hint = _elements(halogen, silicon), _element_hint(halogen, silicon)
    if not formulas:
        warnings.append(f"MS-FINDER found no molecular formula for m/z {precursor_mz:.4f} as {res['adduct']} within 1 mDa (elements {elements}). "
                        "Check the precursor m/z and the adduct." + (" " + hint if hint else ""))
    elif not structures:
        warnings.append(f"MS-FINDER's local structure databases have no candidate for {formulas[0]['formula']} and the other top formulas. "
                        "PubChem is not used, so the compound may simply not be in them." + (" " + hint if hint else ""))
    if len(structures) > STRUCTURE_ROWS:
        warnings.append(f"MS-FINDER returned {len(structures)} structures; the table shows the best {STRUCTURE_ROWS}.")

    top = structures[0] if structures else None
    summary = {
        "Top candidate": ({"value": top["name"] or top["inchikey"] or common.MISSING, "note": f"{top['formula']} · score {fmt(top['score'], 2)}"}
                          if top else common.MISSING),
        "Top formula": ({"value": formulas[0]["formula"], "note": f"score {fmt(formulas[0]['score'], 3)}"} if formulas else common.MISSING),
        "Adduct used": ({"value": res["adduct"], "note": "Auto: " + ("positive" if positive else "negative") + " mode"} if adduct is None else res["adduct"]),
        "Elements searched": {"value": elements, "note": "halogens " + ("on" if halogen else "off") + " · Si " + ("on" if silicon else "off")},
        "Candidates": {"value": f"{len(structures)} / {len(formulas)}",
                       "note": f"MS2 peaks {res['ms2_count']}" + (f" · MS1 peaks {res['ms1_count']}" if res["ms1_count"] else "")}}
    return {"outputs": {
        "summary": {"values": summary},
        "structures": {"rows": _structure_rows(structures), "count": len(structures)},
        "formulas": {"rows": _formula_rows(formulas), "count": len(formulas)}}, "warnings": warnings}


RUNNERS = {"ms2_structure": run_structure}
