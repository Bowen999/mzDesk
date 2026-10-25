"""MS1 小类：理论 m/z（分子式 + 加合物）、精确质量检索（PubChem）、BUDDY 分子式推断。"""
import re

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from mzdesk.api import blocks, common
from mzdesk.api.common import EngineError, InputError, fmt, signed
from mzdesk.shared import colors
from mzdesk.shared.spectrum import parse_spectrum
from mzdesk.spec import ADDUCTS_NEG, ADDUCTS_POS, DATABASES

INSTRUMENTS = {"Orbitrap": "orbitrap", "Q-TOF": "qtof", "FT-ICR": "fticr"}
ELEMENTS_NOTE = "Supported elements: C H N O P S F Cl Br I Na K."
ADDUCT_TABLE_PPM = 100      # the all-adducts table shows Δ ppm only this close to the observed m/z
SEARCH_MAX_PPM = 100        # exact-mass search: widest tolerance


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
        mz, ion_formula = bd.ion_mz_and_formula(formula, a)
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
    mz, ion_formula = ions.get(adduct) or bd.ion_mz_and_formula(formula, adduct)
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


def _id_cell(hit, em):
    """IDs of the databases that list the compound (the database name where PubChem gives no ID); links to the first database page."""
    tokens, url = [], None
    for db in em.SOURCES:
        if db not in hit.dbs:
            continue
        native = hit.ids.get(db)
        tokens.append(f"KEGG {native}" if db == "KEGG" and native else native or db)
        if native and url is None:
            url = em.DB_URL[db].format(native)
    return blocks.link(" · ".join(tokens), url or em.PUBCHEM_URL.format(hit.cid)) if tokens else common.MISSING


def run_search(inputs):
    from mzdesk.ms1 import exact_mass as em

    mz, positive, adduct = common.precursor(inputs, required=True)
    tol = common.number(inputs, "tol_ppm", "Mass tolerance", default=5.0, positive=True, maximum=SEARCH_MAX_PPM)
    dbs = inputs.get("databases")
    dbs = ["HMDB", "PubChem"] if dbs is None else [dbs] if isinstance(dbs, str) else list(dbs)
    if not dbs:
        raise InputError("Choose at least one database.")
    if any(d not in DATABASES for d in dbs):
        raise InputError(f"Unknown database: {', '.join(d for d in dbs if d not in DATABASES)}.")
    pool = (ADDUCTS_POS if positive else ADDUCTS_NEG)[1:]
    if adduct is not None and adduct not in pool:
        raise InputError(f"{adduct} is not a {'positive' if positive else 'negative'}-mode adduct.")
    adducts = pool if adduct is None else [adduct]
    if not em.mass_windows(mz, adducts, tol):
        raise InputError("The precursor m/z is too small for the chosen adduct.")
    try:
        res = em.search(mz, adducts, tol, dbs)
    except em.ExactMassError as exc:
        raise EngineError(str(exc)) from None

    groups = {}         # (formula, adduct) -> hits; the engine's order is closest first
    for h in res.hits:
        groups.setdefault((h.formula, h.adduct), []).append(h)
    table = []
    for (formula, ad), hs in groups.items():
        hs.sort(key=lambda h: h.order)
        n = len(hs)
        table.append(blocks.group(
            f"{formula} · {ad}",
            [[h.name, _id_cell(h, em), blocks.link(str(h.cid), em.PUBCHEM_URL.format(h.cid))] for h in hs],
            meta=f"m/z {hs[0].mz:.4f} · Δ {signed(hs[0].ppm, 1)} ppm · {n} compound{'s' if n != 1 else ''}"))

    chosen = [d for d in em.SOURCES if d in dbs]
    per_db = " · ".join(f"{d} {res.totals[d]:,}" + (f" ({res.shown[d]} listed)" if res.totals[d] > res.shown[d] else "") for d in chosen)
    capped = any(res.totals[d] > res.shown[d] for d in chosen)
    top = res.hits[0] if res.hits else None
    summary = {
        "Compounds": {"value": str(len(res.hits)), "note": f"capped at {em.MAX_PER_DB} per database"} if capped else str(len(res.hits)),
        "Formulas": str(len({h.formula for h in res.hits})),
        "Closest hit": {"value": top.name, "note": f"{top.formula} {top.adduct} · {signed(top.ppm, 1)} ppm"} if top else common.MISSING,
        "Databases searched": {"value": ", ".join(chosen), "note": per_db},
    }
    return {"outputs": {"summary": {"values": summary}, "hits": blocks.grouptable_data(table, count=len(res.hits))},
            "warnings": res.warnings}


RUNNERS = {"ms1_theory": run_theory, "ms1_search": run_search, "ms1_formula": run_formula}
