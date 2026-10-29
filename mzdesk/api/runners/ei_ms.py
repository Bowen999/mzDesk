"""EI MS 页面：结构 → NEIMS 预测的 70 eV 单位质量 EI 谱（[M]+•，无需选加合物）。"""
from mzdesk.api import common
from mzdesk.api.common import EngineError, InputError
from mzdesk.shared import ms2_formats, ms2_plot


def run_ei(inputs):
    from mzdesk.ms2 import neims_engine as ne

    warnings = []
    smiles = common.structure(inputs, warnings=warnings)["smiles"]
    try:
        res = ne.predict(smiles)
    except ne.NEIMSError as exc:
        raise EngineError(str(exc)) from None
    except ValueError as exc:
        raise InputError(str(exc)) from None

    peaks = [(mz, round(rel, 1)) for mz, rel in res["peaks"]]
    mol_ion = res["molecular_ion_mz"]
    base = max(peaks, key=lambda p: p[1])[0] if peaks else None
    name = common.record_name(inputs, res["smiles"])

    fig, cards = ms2_plot.spectrum_figure(peaks, y_title="Relative intensity (%)")
    for card in cards:
        card["rows"] = [["Rel. intensity", f"{card['rel']:.1f} %"]]
        if card["mz"] == mol_ion:
            card["kind"] = "Molecular ion [M]+•"
    return {"outputs": {
        "summary": {"values": {
            "Ion": {"value": "[M]+•", "note": "70 eV EI"},
            "Formula": res["formula"],
            "Molecular ion m/z": {"value": str(mol_ion), "note": "nominal"},
            "Base peak m/z": common.MISSING if base is None else str(base)}},
        "spectrum": common.plotly(fig, peaks=cards),
        "peaks": {"rows": [[str(mz), f"{rel:.1f}"] for mz, rel in peaks]},
        "formats": {"formats": ms2_formats.all_formats(peaks, name=name, ion_mode="Positive", precursor_mz=mol_ion)},
    }, "warnings": warnings + res["warnings"]}


RUNNERS = {"ms2_ei": run_ei}
