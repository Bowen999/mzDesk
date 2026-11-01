"""Calculator page: independent everyday calculations, all optional."""
import math

from mzdesk.api import blocks, common
from mzdesk.api.common import InputError


def _num(value, nd=6):
    if value is None or not math.isfinite(value):
        return common.MISSING
    return f"{value:.{nd}g}"


def _pct(value):
    if value is None or not math.isfinite(value):
        return common.MISSING
    if value == 0:
        return "0"
    if abs(value) < 0.001:
        return f"{value:.2e}"
    return f"{value:.4g}"


def _card(value, note=None):
    return {"value": value, "note": note} if note else value


def run_calculator(inputs):
    from mzdesk.calculator import calculations

    warnings = []
    verdicts = []
    outputs = {}

    theoretical = common.number(inputs, "ppm_theoretical", "Theoretical m/z", positive=True)
    observed = common.number(inputs, "ppm_observed", "Observed m/z", positive=True)
    ppm = None
    if (theoretical is None) != (observed is None):
        warnings.append("Enter both theoretical and observed m/z for the ppm error.")
    if theoretical is not None and observed is not None:
        ppm = calculations.ppm_error(observed, theoretical)
        verdicts.append(f"ppm error <b>{common.signed(ppm, 2)} ppm</b>.")
        ppm_card = _card(common.signed(ppm, 2), f"observed {observed:.4f} vs theoretical {theoretical:.4f}")
    else:
        ppm_card = common.MISSING
    outputs["ppm_error"] = {"values": {"ppm error": ppm_card}}

    lab_smiles = common.text(inputs, "lab_smiles")
    lab_label = common.text(inputs, "lab_label")
    lab_adduct = common.text(inputs, "lab_adduct") or "[M+H]+"
    isotope_rows, crosstalk_rows = [], []
    if lab_smiles or lab_label:
        if not lab_smiles or not lab_label:
            warnings.append("Isotope labelling needs both a SMILES and a label.")
        else:
            try:
                formula, _ = calculations.smiles_formula_mass(lab_smiles, warnings)
                labels = calculations.parse_labels(lab_label)
                result = calculations.labelled_distribution(formula, lab_adduct, labels)
                crosstalk = calculations.crosstalk(formula, lab_adduct, result["nominal_shift"])
            except ValueError as exc:
                raise InputError(str(exc)) from None
            for peak in result["peaks"]:
                if peak["abundance"] < 0.1:
                    continue
                offset = peak["offset"]
                name = f"M ({result['label_text']})" if offset == 0 else f"M+{offset}" if offset > 0 else f"M−{-offset}"
                isotope_rows.append([name, common.fmt(peak["mz"], 5), common.fmt(peak["abundance"], 2)])
            isotope_rows = isotope_rows[:20]
            crosstalk_rows = [[result["label_text"], f"M+{result['nominal_shift']}", _pct(crosstalk)]]
            verdicts.append(f"Labelled ion <b>m/z {result['mono_mz']:.4f}</b>; natural M+{result['nominal_shift']} crosstalk {_pct(crosstalk)}%.")
            if crosstalk >= 0.1:
                warnings.append(f"The unlabelled analyte's natural M+{result['nominal_shift']} peak is {_pct(crosstalk)}% of M; check chromatographic separation.")
    outputs["isotope_labels"] = blocks.table_data(isotope_rows, count=len(isotope_rows))
    outputs["crosstalk"] = blocks.table_data(crosstalk_rows, count=len(crosstalk_rows))

    conc_value = common.number(inputs, "conc_value", "Concentration value", positive=True)
    conc_from = inputs.get("conc_from") or "µM"
    conc_to = inputs.get("conc_to") or "µM"
    conc_mw = common.number(inputs, "conc_mw", "Molar mass", positive=True)
    converted = None
    if conc_value is not None:
        try:
            converted = calculations.concentration(conc_value, conc_from, conc_to, conc_mw)
        except ValueError as exc:
            warnings.append(str(exc))
        else:
            verdicts.append(f"<b>{converted:.6g} {conc_to}</b>.")
    outputs["concentration"] = {"values": {"Converted": _card(_num(converted), f"{conc_to} · from {conc_value:g} {conc_from}") if converted is not None else common.MISSING}}

    target = common.number(inputs, "target_conc", "Target concentration", positive=True)
    final_volume = common.number(inputs, "final_volume", "Final volume", positive=True)
    dilution_rows = []
    if target is not None or final_volume is not None:
        if converted is None or target is None or final_volume is None:
            warnings.append("Dilution needs a stock concentration, a target concentration and a final volume.")
        else:
            try:
                mix = calculations.dilution(conc_value, conc_from, target, conc_to, final_volume, conc_mw)
            except ValueError as exc:
                warnings.append(str(exc))
            else:
                dilution_rows = [[f"{target:g} {conc_to}", f"{mix['stock_in_target']:.6g} {conc_to}",
                                  f"{final_volume:g}", f"{mix['stock_volume']:.4g}", f"{mix['diluent_volume']:.4g}"]]
    outputs["dilution"] = blocks.table_data(dilution_rows, count=len(dilution_rows))

    nce = common.number(inputs, "ce_nce", "NCE", minimum=0)
    ev = common.number(inputs, "ce_ev", "eV", positive=True)
    ce_mz = common.number(inputs, "ce_mz", "Precursor m/z", positive=True)
    charge_value = common.number(inputs, "ce_z", "Charge", default=1, positive=True)
    if charge_value != int(charge_value):
        raise InputError("Charge must be an integer.")
    charge = int(charge_value)
    ce_text = common.MISSING
    ce_note = None
    if nce is not None or ev is not None:
        if ce_mz is None:
            warnings.append("Precursor m/z is required for NCE ↔ eV.")
            parts = []
            if nce is not None:
                parts.append(f"{nce:g} NCE")
            if ev is not None:
                parts.append(f"{ev:g} eV")
            ce_text = " · ".join(parts)
        else:
            if nce is not None:
                ev_calc = calculations.nce_to_ev(nce, ce_mz, charge)
                ce_text = f"{nce:g} NCE · {ev_calc:.2f} eV"
                ce_note = f"m/z {ce_mz:g}, z {charge}"
                if ev is not None and abs(ev - ev_calc) > max(0.1, 0.02 * ev_calc):
                    warnings.append(f"The entered {ev:g} eV does not match {ev_calc:.2f} eV from NCE {nce:g}.")
                verdicts.append(f"<b>{nce:g} NCE ↔ {ev_calc:.2f} eV</b>.")
            else:
                nce_calc = calculations.ev_to_nce(ev, ce_mz, charge)
                ce_text = f"{nce_calc:.2f} NCE · {ev:g} eV"
                ce_note = f"m/z {ce_mz:g}, z {charge}"
                verdicts.append(f"<b>{nce_calc:.2f} NCE ↔ {ev:g} eV</b>.")
    outputs["collision_energy"] = {"values": {"NCE ↔ eV": _card(ce_text, ce_note) if ce_text != common.MISSING else common.MISSING}}

    alkanes_text = common.text(inputs, "alkanes")
    analyte_rt = common.number(inputs, "ri_rt", "Analyte RT", minimum=0)
    kovats_text = common.MISSING
    kovats_note = None
    if alkanes_text or analyte_rt is not None:
        if not alkanes_text or analyte_rt is None:
            warnings.append("Kovats RI needs the alkane series and the analyte RT.")
        else:
            try:
                points = calculations.parse_alkanes(alkanes_text)
                result = calculations.kovats_ri(points, analyte_rt)
            except ValueError as exc:
                raise InputError(str(exc)) from None
            kovats_text = f"{result['ri']:.1f}"
            n1, rt1 = result["lower"]
            n2, rt2 = result["upper"]
            kovats_note = f"C{n1} {rt1:g} → C{n2} {rt2:g}" + (" · extrapolated" if result["extrapolated"] else "")
            verdicts.append(f"Kovats RI <b>{result['ri']:.1f}</b>.")
            if result["extrapolated"]:
                warnings.append("The analyte RT is outside the alkane series; the Kovats RI is extrapolated.")
    outputs["kovats"] = {"values": {"Kovats RI": _card(kovats_text, kovats_note) if kovats_text != common.MISSING else common.MISSING}}

    outputs["notes"] = {"items": [
        "ppm error = (observed − theoretical) / theoretical × 10⁶.",
        "Isotope labelling assumes complete incorporation of the label; the table shows the natural isotope distribution of the remaining atoms, normalised to the most abundant labelled isotopologue.",
        "Crosstalk is the unlabelled analyte's natural M+n abundance relative to its M peak, where n is the label's nominal mass shift.",
        "Mass ↔ molar concentration conversion needs a molar mass; 1 µg/mL = 1 mg/L. Dilution uses the converted stock concentration and C1V1 = C2V2.",
        "NCE ↔ eV uses CE = NCE × precursor m/z × charge / 500; instrument vendors define NCE differently, so treat this as a generic conversion.",
        "Kovats RI is linearly interpolated between adjacent n-alkanes; values outside the measured range are extrapolated and flagged.",
    ]}
    return {"outputs": outputs, "warnings": warnings, "verdict": " ".join(verdicts) or "Enter at least one calculator section."}


RUNNERS = {"calculator": run_calculator}
