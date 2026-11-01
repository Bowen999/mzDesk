"""ESI MS/MS 页：CFM-ID 4 + MSRB 预测谱（低 / 中 / 高三个碰撞能量），输入结构 + 加合物。"""
import re

from mzdesk.api import common
from mzdesk.api.common import EngineError, InputError, fmt
from mzdesk.shared import ms2_formats, ms2_plot

RECORD_FORMATS = ["MSP", "MGF", "MassBank", "SIRIUS .ms"]       # formats that can hold several spectra


def _label(spectrum):
    return f"{spectrum['energy']} ({spectrum['ev']} eV)" if spectrum["ev"] else spectrum["energy"]


def _charge(adduct):
    """Absolute charge of an adduct notation: [M-2H]2- -> 2, [M+H]+ -> 1."""
    m = re.search(r"(\d*)([+-])$", adduct or "")
    return max(1, int(m.group(1) or 1)) if m else 1


def _formats(res, name):
    """{format: text}; one record per energy, separated by a blank line."""
    mode = "Positive" if res["adduct"].endswith("+") else "Negative"
    spectra = [s for s in res["spectra"] if s["peaks"]]
    return {f: "\n\n".join(ms2_formats.to_format(s["peaks"], f, precursor_mz=res["precursor_mz"], ion_mode=mode,
                                                 charge=_charge(res["adduct"]),
                                                 name=f"{name} {res['adduct']} | {_label(s)}") for s in spectra)
            for f in RECORD_FORMATS}


def run_esi(inputs):
    from mzdesk.ms2 import cfmid_engine as cf

    warnings = []
    smiles = common.structure(inputs, label="Structure", warnings=warnings)["smiles"]
    common.note_isotope_labels(warnings, smiles, "The precursor m/z and ion formula are")
    adduct = common.text(inputs, "adduct") or "[M+H]+"
    try:
        res = cf.run_cfmid(smiles, adduct)
    except cf.CFMIDError as exc:
        raise EngineError(str(exc), exc.log) from None
    except ValueError as exc:
        raise InputError(str(exc)) from None

    outputs = {"summary": {"values": {
        "Engine": res["engine"],
        "Precursor m/z": {"value": fmt(res["precursor_mz"], 4), "note": adduct},
        "Ion formula": res["ion_formula"],
        "Peaks": {"value": " · ".join(str(len(s["peaks"])) for s in res["spectra"]),
                  "note": " · ".join(f"{s['ev']} eV" for s in res["spectra"])}}}}
    rows = []
    for s in res["spectra"]:
        fig, cards = ms2_plot.spectrum_figure(s["peaks"], precursor_mz=res["precursor_mz"], y_title="Relative intensity (%)")
        for card, fragment in zip(cards, s["fragments"]):
            if fragment:
                card["rows"].append(["Fragment", fragment])
        outputs[f"spectrum_{s['energy']}"] = common.plotly(fig, peaks=cards)
        rows += [[_label(s), fmt(mz, 4), fmt(rel, 1), fragment or common.MISSING]
                 for (mz, rel), fragment in zip(s["peaks"], s["fragments"])]
    outputs["peaks"] = {"rows": rows}
    outputs["formats"] = {"formats": _formats(res, common.record_name(inputs, res["smiles"]))}
    return {"outputs": outputs, "warnings": warnings + res["warnings"]}


RUNNERS = {"ms2_esi": run_esi}
