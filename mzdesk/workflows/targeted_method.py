"""Targeted method workflow: one compound in, a method sheet out.

Steps (spec.py stages): 1 Identity  2 Platform  3 Separation  4 Ionization  5 Adducts  6 Non-MS detection  7 Derivatization
8 Known methods  9 Theoretical m/z.   Result contract and the step runner: mzdesk/workflows/common.py.

The steps only chain existing engines (Identity, Platform, Separation, Ionization, Adducts, Non-MS detection, Derivatization, Known methods and the
Theoretical m/z code); the decisions that combine their answers (can it be measured, column, pH, main ion, derivatization) are the pure rules of
targeted_rules.py. A step that fails is marked failed and the sheet says "not available" for the rows that need it. No MRM transitions and no
collision energies anywhere.

Form field ids: structure (the structure field: {"mode", "value"}), ph (blank = DEFAULT_PH), additive, instrument.
Instrument: Any = the best platform overall; LC-MS or GC-MS = judge only that platform (the other platform's steps are skipped when they are not needed).
"""
from __future__ import annotations

import copy
import logging
import re
from dataclasses import dataclass, field

from mzdesk.api import blocks as blk
from mzdesk.api import common
from mzdesk.api.common import InputError
from mzdesk.spec import MODULES, PIPELINES, TARGETED_INSTRUMENTS
from mzdesk.workflows import targeted_rules as rules
from mzdesk.workflows.common import STATUS_DONE, STATUS_SKIPPED, StepRunner, make_answer, make_result, table_block

log = logging.getLogger("mzdesk.workflows")

PIPELINE_ID = "targeted_method"
STEP_TITLES = {s["id"]: s["title"] for s in next(p for p in PIPELINES if p["id"] == PIPELINE_ID)["stages"]}
DEFAULT_PH = 2.7                # 0.1 % formic acid, the Ionization page's default
PH_MIN, PH_MAX = 0.0, 14.0
DEFAULT_ADDITIVE = "Formic acid"
NAME_MODE = "Name / ID / Formula"

ID_TIMEOUT = 8.0                # seconds per PubChem / UniChem request in the Identity step
ISOTOPE_TOLERANCE = 0.0005      # Da: the adduct code and the isotope code must give the same ion m/z within this
EXACT_ROWS_SHOWN = 5            # rows of the exact-match table in step 8 (the Known methods page lists all)
SHEET_ITEMS = ["Can it be measured", "Platform", "Column", "Mobile phase / pH", "Polarity", "Main ion", "Theoretical m/z", "Derivatization",
               "Non-MS alternative", "Published methods"]
VERDICT_TONE = {rules.VERDICT_YES: "good", rules.VERDICT_POSSIBLE: "warn", rules.VERDICT_UNLIKELY: "bad"}
NA = rules.NOT_AVAILABLE


@dataclass
class Query:
    """The form, parsed. smiles is the resolved structure (a name or ID is looked up in parse_inputs, online); name is what the user typed in
    Name / ID mode, else "". ph is DEFAULT_PH when blank (ph_given False). additive is one of adduct_rules.ADDITIVES, instrument one of
    spec.TARGETED_INSTRUMENTS."""
    smiles: str
    name: str = ""
    ph: float = DEFAULT_PH
    ph_given: bool = False
    additive: str = DEFAULT_ADDITIVE
    instrument: str = "Any"
    warnings: list = field(default_factory=list)


def parse_inputs(inputs):
    """Form JSON -> Query. Raises InputError (HTTP 400) with a message for the user."""
    from mzdesk.platform_fit import adduct_rules

    warnings = []
    found = common.structure(inputs, "structure", "Compound", warnings=warnings)          # a name or ID is resolved online here
    typed = inputs.get("structure") or {}
    name = str(typed.get("value") or "").strip() if typed.get("mode") == NAME_MODE else ""
    ph = common.number(inputs, "ph", "pH", minimum=PH_MIN, maximum=PH_MAX)
    additive = common.text(inputs, "additive") or DEFAULT_ADDITIVE
    if additive not in adduct_rules.ADDITIVES:
        raise InputError(f"Mobile-phase additive must be one of: {', '.join(adduct_rules.ADDITIVES)}.")
    instrument = common.text(inputs, "instrument") or TARGETED_INSTRUMENTS[0]
    if instrument not in TARGETED_INSTRUMENTS:
        raise InputError(f"Instrument must be one of: {', '.join(TARGETED_INSTRUMENTS)}.")
    return Query(smiles=found["smiles"], name=name, ph=DEFAULT_PH if ph is None else ph, ph_given=ph is not None, additive=additive,
                 instrument=instrument, warnings=warnings)


# ---------------------------------------------------------------- helpers
def _inputs(smiles, **extra):
    """The form of a module runner for this structure."""
    return {"structure": {"mode": "SMILES / InChI", "value": smiles}, **extra}


def _module_block(module_id, block_id, data, new_id=None):
    """The spec block `block_id` of a module page with the data its runner returned; the id is made unique within the workflow."""
    module = next(m for m in MODULES if m["id"] == module_id)
    spec_block = next(b for b in module["outputs"] if b.get("id") == block_id)
    return {**copy.deepcopy(spec_block), "id": new_id or f"{module_id}_{block_id}", "data": data}


def _clip(text, n=110):
    text = " ".join(str(text).split())
    return text if len(text) <= n else text[:n - 1].rstrip() + "…"


# ---------------------------------------------------------------- step 1: Identity
def _lookup(inchikey, smiles):
    """PubChem record of the structure. The InChIKey finds the parent compound (metformin -> Metformin); a SMILES search can answer with an isotopologue
    ('[14C]metformin'), so it is only the fallback for a structure whose key PubChem does not know (or that has none)."""
    from mzdesk.identity import resolve

    if inchikey:
        try:
            return resolve.resolve(inchikey, "InChIKey", timeout=ID_TIMEOUT)
        except resolve.ResolveError as exc:
            if "not found" not in str(exc):         # unreachable, timeout ...: the SMILES search would fail the same way
                raise
    return resolve.resolve(smiles, "SMILES", timeout=ID_TIMEOUT)


def step_identity(q):
    """Name, formula, monoisotopic mass, InChIKey and database IDs. PubChem / UniChem through resolve.resolve (not the slow enrich step);
    offline or a failed lookup falls back to RDKit values with a note."""
    from rdkit import Chem
    from rdkit.Chem import inchi, rdMolDescriptors

    mol = Chem.MolFromSmiles(q.smiles)
    if mol is None:
        raise ValueError("RDKit could not parse the structure.")
    try:
        key = inchi.MolToInchiKey(mol) or ""
    except Exception:           # noqa: BLE001 - InChI generation can fail for odd structures; the key is optional
        key = ""
    ident = {"name": q.name, "formula": rdMolDescriptors.CalcMolFormula(mol), "mass": rdMolDescriptors.CalcExactMolWt(mol), "inchikey": key,
             "cid": None, "ids": {}, "links": {}, "online": False}
    problem = ""
    try:
        r = _lookup(key, q.smiles)
    except Exception as exc:    # noqa: BLE001 - offline or any lookup failure must not stop the workflow
        log.info("identity lookup failed: %s", exc)
        r, problem = None, " ".join(str(exc).split())
    if r is not None:
        ident.update(name=r.name or q.name, inchikey=r.inchikey or key, cid=r.cid, ids=dict(r.ids), links=dict(r.links), online=r.cid is not None)
        if r.cid is None:
            problem = "the structure is not in PubChem"
    ident["name"] = ident["name"] or ident["formula"]

    def first(label):
        values = ident["ids"].get(label)
        return (values[0] if isinstance(values, (list, tuple)) else values) if values else None

    rows = [["Name", ident["name"]], ["Formula", ident["formula"]], ["Monoisotopic mass", f"{ident['mass']:.4f}"], ["InChIKey", ident["inchikey"] or NA]]
    if ident["cid"]:
        rows.append(["PubChem CID", blk.link(str(ident["cid"]), ident["links"].get("PubChem"))])
    for label in ("CAS", "HMDB", "KEGG", "ChEBI"):
        value = first(label)
        if value:
            rows.append([label, blk.link(str(value), ident["links"].get(label)) if ident["links"].get(label) else str(value)])
    why = "PubChem unreachable" if "unreachable" in problem else _clip(problem, 60)
    note = (ident["name"] if ident["name"] == ident["formula"] else f"{ident['name']} · {ident['formula']}") + ("" if ident["online"] else f" · RDKit values only ({why})")
    return STATUS_DONE, note, [table_block("identity", "Identity", ["Item", "Value"], rows)], ident


# ---------------------------------------------------------------- step 2: Platform
def step_platform(q):
    """The four Platform tiles (Lipidomics, Metabolomics, CIL Metabolomics, GC-MS) and the verdicts behind the decisions. The verdicts are those of the
    neutral form of the largest fragment (a salt or mixture is analysed by its main component)."""
    from rdkit import Chem

    from mzdesk.api.runners import structure
    from mzdesk.platform_fit import derivatization as deriv
    from mzdesk.platform_fit import ionization_rules as ir

    gc, mol = ir.gc_verdict(q.smiles)
    main = Chem.MolToSmiles(mol)
    if Chem.MolFromSmiles(main) is None:             # an odd structure (organometallic ...) may not survive the round trip: keep the input
        main = q.smiles
    out = structure.run_platform_fit(_inputs(main))["outputs"]
    tiles = {t["label"]: t for t in out["platform"]["tiles"]}
    groups = deriv.find_groups(mol)
    cil = deriv.cil_assess(groups)
    option = deriv.gc_option(mol, groups)
    lipid, metabolite = tiles["Lipidomics"]["status"] == "yes", tiles["Metabolomics"]["status"] == "yes"
    channels = [(short, deriv.reagent_label(opt).split(" · ", 1)[1]) for short, opt in zip(cil["short"], cil["options"])]
    value = {"gc": gc, "cil": cil, "channels": channels, "gc_option": option, "gc_reagent": deriv.REAGENTS[option][0], "lipid": lipid,
             "metabolite": metabolite, "omics": rules.omics_label(lipid, metabolite), "main_smiles": main, "reactions": out.get("reactions")}
    note = (value["omics"] or "outside lipidomics and metabolomics") + f"; GC-MS {gc['gc']}; CIL " + ("yes" if cil["cil"] else "no")
    return STATUS_DONE, note, [_module_block("platform_fit", "platform", out["platform"], "platform")], value


# ---------------------------------------------------------------- step 3: Separation
def step_separation(q):
    """Separation modes at the pH (MolGpKa pKa -> logD, charge) and the robust pH ranges. Without the MolGpKa weights the step is skipped."""
    from mzdesk.api.runners import structure
    from mzdesk.platform_fit import pka_engine as pe
    from mzdesk.platform_fit import separation_rules as sr

    try:
        a = pe.analyze(q.smiles)
    except pe.PkaError as exc:
        return STATUS_SKIPPED, f"MolGpKa weights not installed: separation mode and pH range not available ({_clip(exc, 90)})", [], None
    st = a.at(q.ph)
    rows = sr.mode_rows(ph=q.ph, logd=st["logd"], tpsa=a.tpsa, mw=a.mw, charge=st["charge"], f_neutral=st["f_neutral"], permanent=a.permanent,
                        fixed_charge=a.fixed_charge, n_stereo=a.n_stereo, n_unassigned=a.n_unassigned, meso=a.meso)
    levels = rules.mode_levels(rows)
    mode, level = rules.best_mode(levels)
    why = {r[0]: r[2] for r in rows}
    ranges = sr.robust_ph_ranges(a.sites, a.ph)
    out = structure.run_separation(_inputs(q.smiles, ph=q.ph))["outputs"]
    blocks = [_module_block("separation", "modes", out["modes"], "separation_modes"), _module_block("separation", "ph_summary", out["ph_summary"], "separation_ph"),
              _module_block("separation", "rule_hits", out["rule_hits"], "separation_notes")]
    value = {"rows": rows, "levels": levels, "mode": mode, "level": level, "why": why, "ranges": ranges, "n_sites": len(a.sites),
             "n_stereo": 0 if a.meso else a.n_stereo, "warnings": list(a.notes)}
    return STATUS_DONE, f"best mode {mode} ({level}); " + ("no ionizable site" if not a.sites else f"robust pH {rules.range_text(ranges)}" if ranges else "no robust pH"), blocks, value


# ---------------------------------------------------------------- step 4: Ionization
def step_ionization(q):
    """ESI+ / ESI- response at the pH, polarity, GC amenability and the EI molecular ion (ionization_rules.assess)."""
    from mzdesk.api.runners import ionization
    from mzdesk.platform_fit import ionization_rules as ir

    a = ir.assess(q.smiles, q.ph)
    out = ionization.run_ionization(_inputs(q.smiles, ph=q.ph))["outputs"]
    blocks = [_module_block("ionization", "chips", out["chips"], "ionization_chips"), _module_block("ionization", "rule_hits", out["rule_hits"], "ionization_notes")]
    value = {"pos": a.pos.level, "neg": a.neg.level, "polarity": rules.polarity_char(a.polarity), "polarity_note": a.polarity_note, "gc": a.gc,
             "ion": a.ion, "warnings": list(a.warnings)}
    return STATUS_DONE, f"ESI+ {a.pos.level}, ESI{ir.MINUS} {a.neg.level}; polarity {a.polarity}", blocks, value


# ---------------------------------------------------------------- step 5: Adducts
def step_adducts(q, polarity):
    """Ranked ions for both polarities; the main ion is the best one of the recommended polarity (the higher score for both / poor)."""
    from mzdesk.api.runners import adducts
    from mzdesk.platform_fit import adduct_rules as ar

    res = ar.predict(q.smiles, "Both", q.additive)
    sign, ion, why = rules.pick_main_ion(polarity or "both", res.best["+"], res.best["-"])
    out = adducts.run_adducts(_inputs(q.smiles, ion_mode="Both", additive=q.additive))["outputs"]
    blocks = [_module_block("adducts", "ions", out["ions"], "adducts_ions")]
    value = {"sign": sign, "ion": ion, "why": why, "formula": res.formula, "warnings": list(res.warnings)}
    note = f"{ion.adduct} in {rules.esi_label(sign)} ({ion.likelihood})" if ion else why
    return STATUS_DONE, note, blocks, value


# ---------------------------------------------------------------- step 6: Non-MS detection
def _non_ms_texts(best, a):
    """(recommendation, basis) of the best non-MS detector: short sentences from the analysis."""
    name, suit, reason = best
    if suit == "poor":
        return "None suitable", _clip(reason.split(":")[0].split(". ")[0], 90)
    if name == "UV/PDA":
        return f"UV/PDA at ~{a.lam:.0f} nm", f"{a.strength} chromophore, λmax ~{a.lam:.0f} nm" + ("" if suit == "good" else " (low sensitivity)")
    if name == "FLD":
        return "Fluorescence (FLD)", f"native fluorescence {a.fluor}: {_clip(a.fluor_reason, 60)}"
    if name in ("CAD", "ELSD"):
        return "CAD or ELSD", f"non-volatile (MW {a.mw:.0f}), no chromophore needed"
    return "RI (isocratic only)", "universal but insensitive; every other detector is poor"


def step_non_ms(q):
    """UV / FLD / CAD / ELSD / RI suitability (non_ms_rules); the best one is the non-MS alternative."""
    from mzdesk.api.runners import non_ms
    from mzdesk.platform_fit import non_ms_rules as nm

    try:
        a = nm.analyze(q.smiles)
    except nm.NonMsError as exc:
        raise ValueError(str(exc)) from None
    best = rules.best_non_ms(a.detectors)
    out = non_ms.run_non_ms(_inputs(q.smiles))["outputs"]
    text, basis = _non_ms_texts(best, a) if best else (NA, "")
    value = {"best": best, "text": text, "basis": basis, "warnings": list(a.warnings)}
    return STATUS_DONE, f"best: {best[0]} ({best[1]})" if best else "no detector rated", [_module_block("non_ms", "detectors", out["detectors"], "non_ms_detectors")], value


# ---------------------------------------------------------------- step 7: Derivatization
def step_derivatization(q, platform, ionization, decision):
    """LC-MS: a CIL label only when the ESI response is weak; GC-MS: from the GC verdict, with the derivative's formula and mass."""
    from mzdesk.platform_fit import derivatization as deriv

    gc = (platform or {}).get("gc") or (ionization or {}).get("gc")
    sign = _chosen_sign(ionization)
    level = rules.level_in(sign, ionization["pos"], ionization["neg"]) if sign else None
    lc = rules.lc_derivatization(level, platform["channels"] if platform else [])
    gcd = rules.gc_derivatization(gc["gc"] if gc else None, platform["gc_reagent"] if platform else "TMS", gc["volatility"] if gc else "")
    derivative = None
    if gcd and gcd.kind == "needed" and platform:
        res = deriv.derive(platform["main_smiles"], [platform["gc_option"]])["results"]
        if res and res[0]["n"] and res[0]["formula"]:
            derivative = {"formula": res[0]["formula"], "mass": res[0]["mass"], "name": res[0]["name"], "ions": res[0]["ions"], "warns": res[0]["warns"]}
    rows = [["LC-MS", lc.text if lc else NA, lc.basis if lc else ""], ["GC-MS", gcd.text if gcd else NA, _clip(gcd.basis, 110) if gcd else ""]]
    blocks = [table_block("derivatization", "Derivatization", ["Platform", "Verdict", "Reason"], rows)]
    wanted = (decision.ms_route == rules.ROUTE_GC and gcd and gcd.kind == "needed") or (decision.ms_route == rules.ROUTE_LC and lc and lc.kind == "suggested")
    reactions = (platform or {}).get("reactions")
    if wanted and reactions and reactions.get("options"):
        blocks.append(_module_block("platform_fit", "reactions", reactions, "derivatization_reactions"))
    value = {"lc": lc, "gc": gcd, "derivative": derivative}
    chosen = gcd if decision.ms_route == rules.ROUTE_GC else lc
    note = (chosen.text if chosen else "not available") + (f" ({derivative['formula']})" if derivative and chosen is gcd else "")
    return STATUS_DONE, note, blocks, value


def _chosen_sign(ionization):
    """'+' or '-': the polarity the LC-MS route uses (the better ESI response; '+' when both are equal or the step is missing)."""
    if not ionization:
        return None
    pol = ionization["polarity"]
    if pol in ("+", "-"):
        return pol
    return "-" if rules.ESI_RANK[ionization["neg"]] > rules.ESI_RANK[ionization["pos"]] else "+"


# ---------------------------------------------------------------- step 8: Known methods
def method_counts(res):
    """-> (JP18 exact methods, RepoRT exact methods): the true totals of the search (Result.exact_totals). Result.exact lists at most
    cfg.max_exact_per_source methods per compound."""
    return res.exact_totals.get("JP18", 0), res.exact_totals.get("RepoRT", 0)


def step_known_methods(q):
    """Published exact-match methods (JP18, RepoRT when installed) and, without an exact match, the methods of the most similar compounds."""
    from mzdesk.methods import known_methods as km

    res = km.search(q.smiles)
    out = km.page_outputs(res)
    n_jp, n_report = method_counts(res)
    listed, total = len(res.exact), n_jp + n_report                  # the search lists at most cfg.max_exact_per_source methods per compound
    summary = copy.deepcopy(out["summary"])
    if total != listed:                                              # say the true totals in the cards and the table, as the sheet does
        n_similar = len(res.methods) - listed
        summary["values"]["Methods found"] = {"value": str(total + n_similar), "note": f"{total} exact · {n_similar} for {len(res.similar)} similar compounds"}
    exact = dict(out["exact"], rows=out["exact"]["rows"][:EXACT_ROWS_SHOWN], count=total)
    if total > EXACT_ROWS_SHOWN:
        exact["help"] = f"The first {EXACT_ROWS_SHOWN} of {total}; the Known Methods page lists up to {res.cfg.max_exact_per_source} per compound and source."
    blocks = [_module_block("known_methods", "summary", summary, "methods_summary"), _module_block("known_methods", "exact", exact, "methods_exact")]
    if not listed and res.similar:
        blocks.append(_module_block("known_methods", "similar", dict(out["similar"], groups=out["similar"]["groups"][:3]), "methods_similar"))
    value = {"res": res, "exact": list(res.exact), "n_jp": n_jp, "n_report": n_report, "report_installed": res.report_installed}
    note = f"{total} exact, {len(res.similar)} similar" + (" (JP18 only: RepoRT not installed)" if res.report_installed is False else "")
    return STATUS_DONE, note, blocks, value


# ---------------------------------------------------------------- step 9: Theoretical m/z
ELECTRON = rules.ELECTRON       # u, as in adduct_rules / msbuddy_engine


def _plain_formula(formula):
    """C8H10N4O2 from 'C8H10N4O2+' (the charge suffix removed)."""
    return re.sub(r"[+-]\d*$", "", formula or "")


def step_theoretical_mz(q, decision, platform, adducts, derivatization):
    """m/z and ion formula of the main ion, checked against the isotope code (mzdesk/ms1/isotopes.py ion_mz): both must agree within 0.0005.
    LC route: the ion of step 5 (adduct code). GC route: M+• of what is injected (the derivative of step 7, else the compound): its mass minus one electron."""
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    from mzdesk.ms1 import isotopes

    gc_route = decision.ms_route == rules.ROUTE_GC
    if not gc_route:
        ion = (adducts or {}).get("ion")
        if ion is None or ion.mz is None:
            return STATUS_SKIPPED, "no main ion to calculate", [], None
        neutral, adduct, mz, ion_formula, label = adducts["formula"], ion.adduct, ion.mz, ion.formula, ion.adduct
    else:
        derivative = (derivatization or {}).get("derivative")
        if derivative:
            neutral, mass = _plain_formula(derivative["formula"]), derivative["mass"]
        else:
            mol = Chem.MolFromSmiles((platform or {}).get("main_smiles") or q.smiles)
            neutral, mass = _plain_formula(rdMolDescriptors.CalcMolFormula(mol)), rdMolDescriptors.CalcExactMolWt(mol)
        adduct, label, mz, ion_formula = "[M]+", "M+•", mass - ELECTRON, None
    try:
        iso_mz, iso_formula = isotopes.ion_mz(neutral, adduct)
    except ValueError as exc:
        check, delta = f"not checked ({_clip(exc, 60)})", None
        iso_formula = None
    else:
        delta = iso_mz - mz
        check = "agrees" if abs(delta) <= ISOTOPE_TOLERANCE else f"differs by {delta:+.4f}"
    if gc_route:
        ion_formula = (iso_formula or neutral + "+").replace("+", "+•")
    value = {"mz": mz, "ion_formula": ion_formula, "ion": label, "check": check, "delta": delta}
    if delta is not None and abs(delta) > ISOTOPE_TOLERANCE:
        value["warning"] = f"The adduct code and the isotope code disagree on the m/z of {label} by {delta:+.4f}."
    rows = [[label, ion_formula or NA, f"{mz:.4f}", check]]
    return STATUS_DONE, f"{mz:.4f} {label}; {check}", [table_block("main_ion", "Main ion", ["Ion", "Formula", "m/z", "Isotope check"], rows)], value


# ---------------------------------------------------------------- the method sheet
def _gc_verdict(platform, ionization):
    """The GC verdict (yes / after derivatization / no): the Platform step's, else the Ionization step's; None when neither ran."""
    gc = (platform or {}).get("gc") or (ionization or {}).get("gc")
    return gc["gc"] if gc else None


def _decide(q, got, with_non_ms=True):
    """rules.decide on what the steps found (a missing step is None)."""
    platform, ionization, non_ms = got.get(2), got.get(4), got.get(6)
    esi = rules.esi_best(ionization["pos"], ionization["neg"]) if ionization else None
    best = non_ms["best"] if with_non_ms and non_ms and non_ms.get("best") else None
    return rules.decide(q.instrument, esi, _gc_verdict(platform, ionization), (best[0], best[1]) if best else None, (platform or {}).get("omics", ""))


def _row(item, recommendation, basis=""):
    return [item, recommendation, basis]


def _basis(text, *numbers):
    """'reason (step 4)' / 'reason (steps 3, 8)'."""
    where = ("step " if len(numbers) == 1 else "steps ") + ", ".join(str(n) for n in numbers)
    return f"{text} ({where})" if text else f"({where})"


def _column_row(q, decision, got):
    """-> (ColumnChoice, row)."""
    platform, separation, ionization, methods = got.get(2), got.get(3), got.get(4), got.get(8)
    gc_route = decision.ms_route == rules.ROUTE_GC
    levels = separation["levels"] if separation else None
    choice = rules.choose_column(rules.ROUTE_GC if gc_route else rules.ROUTE_LC, levels, methods["exact"] if methods else [],
                                 _gc_verdict(platform, ionization) == "yes", separation["n_stereo"] if separation else 0)
    if choice.text is None:
        return choice, _row("Column", NA, _basis("separation modes not available and no published column", 3, 8))
    note = f"; {choice.note}" if choice.note else ""
    if choice.source == "published":
        ref = choice.method.reference[0]
        if gc_route:
            return choice, _row("Column", choice.text, _basis(f"published GC method: {ref}", 8))
        if levels is None:
            return choice, _row("Column", choice.text, _basis(f"published: {ref}; separation rules not available", 8))
        mode = rules.CLASS_MODE[rules.column_class(choice.text)]
        rt = f"RT {choice.method.rt} min; " if getattr(choice.method, "rt", None) else ""
        return choice, _row("Column", choice.text, _basis(f"published: {ref}; {rt}{mode} is rated {levels[mode]}{note}", 3, 8))
    if gc_route:
        if _gc_verdict(platform, ionization) == "yes":
            return choice, _row("Column", choice.text, _basis("no published GC column; default phase", 2, 8))
        return choice, _row("Column", choice.text, _basis("default phase; derivatives run on the same column", 7, 8))
    return choice, _row("Column", choice.text, _basis(f"rule: {separation['mode']} is rated {separation['level']} ({_clip(separation['why'][separation['mode']], 70)}){note}", 3))


def _ion_rows(q, decision, got):
    """Polarity and Main ion rows."""
    ionization, adducts, deriv = got.get(4), got.get(5), got.get(7)
    if decision.ms_route == rules.ROUTE_GC:
        if deriv and deriv.get("derivative"):
            d = deriv["derivative"]
            fragment = rules.silyl_fragment(d["name"], d["mass"])
            stronger = f"; [{fragment[0]}]+ {fragment[1]:.4f} is usually stronger" if fragment else ""
            main = _row("Main ion", f"M+• ({d['name']} derivative)", _basis("EI molecular ion of the derivative; it can be weak" + stronger, 7))
        else:
            ion = ionization["ion"] if ionization else None
            main = _row("Main ion", "M+•", _basis(f"NEIMS predicts M+• {ion.value}" if ion and ion.value else "EI molecular ion", 4))
        return [_row("Polarity", "EI (positive ions)", "GC-MS platform"), main]
    sign = adducts["sign"] if adducts else None
    if ionization is None:
        polarity = _row("Polarity", NA, _basis("ionization step not available", 4))
    else:
        levels = f"ESI+ {ionization['pos']}, ESI{rules.MINUS} {ionization['neg']}"
        pol = ionization["polarity"]
        if pol == "poor":
            polarity = _row("Polarity", f"{rules.esi_label(sign)} (weak)" if sign else "ESI response weak", _basis(f"{levels}; ESI response weak, consider APCI or GC", 4))
        elif pol == "both":
            polarity = _row("Polarity", rules.esi_label(sign) if sign else "ESI+ or ESI\u2212",
                            _basis(f"{levels}; both respond" + (f", {rules.esi_label(sign)} scores higher" if sign else ""), 4, 5))
        else:
            polarity = _row("Polarity", rules.esi_label(pol), _basis(levels, 4))
    if adducts is None:
        main = _row("Main ion", NA, _basis("adducts step not available", 5))
    elif adducts["ion"] is None:
        main = _row("Main ion", "none expected", _basis(adducts["why"], 5))
    else:
        i = adducts["ion"]
        weak = ionization is not None and ionization["polarity"] == "poor"
        main = _row("Main ion", i.adduct, _basis(f"{rules.esi_label(adducts['sign'])} top ion, {i.likelihood} likelihood, score {i.score:g}"
                                                  + (f"; {adducts['why']}" if weak else ""), 5))
    return [polarity, main]


def build_sheet(q, got, decision):
    """-> (rows, parts): the ten rows [Item, Recommendation, Basis] and the pieces the statement needs. got = {step number: value, None when the
    step failed or was skipped}."""
    separation, adducts, non_ms, deriv, methods, mz = got.get(3), got.get(5), got.get(6), got.get(7), got.get(8), got.get(9)
    platform = got.get(2)
    gc_route = decision.ms_route == rules.ROUTE_GC
    parts = {"sign": adducts["sign"] if adducts else None}

    cited = (2, 4, 6) if decision.route == rules.ROUTE_NONMS else (2, 4)
    verdict_cell = blk.badge(decision.verdict, VERDICT_TONE[decision.verdict]) if decision.verdict else NA
    how = {"Any": "best platform overall", "LC-MS": "instrument set to LC-MS", "GC-MS": "instrument set to GC-MS"}[q.instrument]
    omics = f"; Platform tiles: {platform['omics']}" if platform and platform["omics"] and not gc_route else ""
    rows = [_row("Can it be measured", verdict_cell, _basis(decision.reason, *cited)), _row("Platform", decision.platform, _basis(how + omics, 2, 4))]

    choice, column = _column_row(q, decision, got)
    rows.append(column)
    parts["column_kind"], parts["column_source"] = choice.kind, choice.source

    if gc_route:
        rows.append(_row("Mobile phase / pH", "Not applicable (GC)", "GC-MS platform"))
    else:
        phc = rules.ph_choice(q.ph, q.ph_given, separation["ranges"] if separation else None, separation["n_sites"] if separation else 1)
        parts["ph"] = q.ph
        rows.append(_row("Mobile phase / pH", f"{rules.ADDITIVE_TEXT[q.additive]}, pH {q.ph:g}", _basis(phc.text, 3)))
    rows += _ion_rows(q, decision, got)

    if mz:
        parts["mz"], parts["ion_label"] = mz["mz"], mz["ion"]
        rows.append(_row("Theoretical m/z", f"{mz['mz']:.4f}", _basis(f"{mz['ion_formula'] or mz['ion']}; isotope check: {mz['check']}", 9)))
    else:
        rows.append(_row("Theoretical m/z", NA, _basis("no main ion", 9)))

    chosen = (deriv["gc"] if gc_route else deriv["lc"]) if deriv else None
    if chosen is None:
        rows.append(_row("Derivatization", NA, _basis("verdict not available", 7)))
    else:
        d = deriv.get("derivative") if gc_route else None
        extra = ""
        if not gc_route and deriv["gc"] and deriv["gc"].kind == "needed":
            extra = f"; GC-MS would need {deriv['gc'].text.removeprefix('Needed: ')}"
        rows.append(_row("Derivatization", chosen.text + (f" ({d['formula']}, {d['mass']:.2f} Da)" if d else ""), _basis(_clip(chosen.basis, 90) + extra, 7)))
        parts["deriv_kind"], parts["deriv_text"] = chosen.kind, chosen.text

    if non_ms and non_ms.get("best"):
        rows.append(_row("Non-MS alternative", non_ms["text"], _basis(non_ms["basis"], 6)))
    else:
        rows.append(_row("Non-MS alternative", NA, _basis("detector step not available", 6)))

    if methods:
        res, exact = methods["res"], methods["exact"]
        total = methods["n_jp"] + methods["n_report"]
        gc_exact = [m for m in exact if rules.column_class(m.column) == "GC"]
        if gc_route:                    # a GC-MS user needs a GC method: the LC methods (RepoRT) only show that the compound is known
            top = choice.method if choice.source == "published" else gc_exact[0] if gc_exact else None
            basis = f"top: {top.reference[0]}" if top else "the exact matches are LC methods" if exact else "no exact match"
        else:
            top = choice.method if choice.source == "published" else exact[0] if exact else None
            basis = f"top: {top.reference[0]}" if top else None
        if basis is None:
            basis = f"best similar: {res.similar[0].name} (score {res.similar[0].score:.2f})" if res.similar else "no similar compound found either"
        if methods["report_installed"] is False:
            basis += "; JP18 only, RepoRT not installed"
        if gc_route:
            text = rules.gc_methods_text(len(gc_exact), total - len(gc_exact))
            parts["methods"], parts["methods_gc"] = len(gc_exact), True
        else:
            text = rules.methods_count_text(methods["n_jp"], methods["n_report"], len(res.similar))
            parts["methods"], parts["methods_gc"] = total, False
        rows.append(_row("Published methods", text, _basis(basis, 8)))
    else:
        rows.append(_row("Published methods", NA, _basis("known-methods search not available", 8)))
    assert [r[0] for r in rows] == SHEET_ITEMS
    return rows, parts


DERIV_CLAUSE = {"not needed": "no derivatization", "none": "no derivatization available", "not feasible": "derivatization does not make it GC-amenable"}


def statement(name, decision, parts):
    """One short sentence for the answer: verdict and platform, ion and m/z, column and pH, derivatization, published methods. It says what the
    sheet says: no ion for a route without MS (HPLC with a non-MS detector, nothing suitable), no derivatization for a non-MS detector, the published GC
    column when the sheet recommends it."""
    if decision.verdict is None:
        return f"{name}: the method sheet is incomplete; see the steps that failed."
    verb = {rules.VERDICT_YES: "measurable by", rules.VERDICT_POSSIBLE: "possible by", rules.VERDICT_UNLIKELY: "unlikely to be measured by"}[decision.verdict]
    gc_route = decision.ms_route == rules.ROUTE_GC
    no_route = decision.route == rules.ROUTE_NONE
    out = f"{name}: {verb} {'LC-MS or GC-MS' if no_route else decision.platform}"
    if parts.get("mz") is not None and (gc_route or parts.get("sign")) and decision.route not in (rules.ROUTE_NONMS, rules.ROUTE_NONE):
        out += f", {'EI' if gc_route else rules.esi_label(parts['sign'])} {parts['ion_label']} m/z {parts['mz']:.4f}"
    if parts.get("column_kind") and not no_route:
        if gc_route:
            out += " on the published GC column" if parts.get("column_source") == "published" else " on a DB-5ms-type column"
        else:
            out += f" on a {parts['column_kind']} column" + (f" at pH {parts['ph']:g}" if "ph" in parts else "")
    kind = None if decision.route == rules.ROUTE_NONMS else parts.get("deriv_kind")           # a label raises the ESI response: it says nothing about a non-MS detector
    if kind in DERIV_CLAUSE:
        out += "; " + DERIV_CLAUSE[kind]
    elif kind == "suggested":
        out += f"; {parts['deriv_text']} would raise the response"
    elif kind == "needed":
        out += f"; derivatization: {parts['deriv_text'].removeprefix('Needed: ')}"
    if "methods" in parts:
        n, what = parts["methods"], "published GC method" if parts.get("methods_gc") else "published method"
        out += "; " + (f"{n} {what}{'' if n == 1 else 's'}" if n else f"no {what}")
    return out + "."


def run_query(q, progress=None):
    """The nine steps in order; each is wrapped in steps.guarded, so a failing engine marks its step failed and the sheet says "not available"
    for the rows that need it. Steps that the instrument choice makes unnecessary are skipped with a note."""
    steps = StepRunner(STEP_TITLES, progress)
    got = {}
    got[1] = steps.guarded(1, lambda: step_identity(q))
    got[2] = steps.guarded(2, lambda: step_platform(q))
    if q.instrument == "GC-MS":
        steps.skip(3, "GC-MS selected: no LC separation needed.")
    else:
        got[3] = steps.guarded(3, lambda: step_separation(q))
    got[4] = steps.guarded(4, lambda: step_ionization(q))

    decision = _decide(q, got, with_non_ms=False)          # which MS route the ion steps follow (the non-MS detector only matters for the verdict)
    if decision.ms_route == rules.ROUTE_GC:
        steps.skip(5, "GC-MS route: the ions are EI, not ESI adducts.")
    else:
        got[5] = steps.guarded(5, lambda: step_adducts(q, (got[4] or {}).get("polarity")))
    got[6] = steps.guarded(6, lambda: step_non_ms(q))
    decision = _decide(q, got)                              # final verdict, now with the non-MS fallback
    got[7] = steps.guarded(7, lambda: step_derivatization(q, got[2], got[4], decision))
    got[8] = steps.guarded(8, lambda: step_known_methods(q))
    got[9] = steps.guarded(9, lambda: step_theoretical_mz(q, decision, got[2], got.get(5), got[7]))

    rows, parts = build_sheet(q, got, decision)
    ident = got[1]
    name = (ident or {}).get("name") or q.name or "The compound"
    facts = [["Name", name]] + ([["Formula", ident["formula"]], ["Monoisotopic mass", f"{ident['mass']:.4f}"], ["InChIKey", ident["inchikey"] or NA]] if ident else [])
    sheet = table_block("method_sheet", "Method sheet", ["Item", "Recommendation", "Basis"], rows, span="full",
                        help="Rule-based starting point built from the nine steps.")
    return make_result(PIPELINE_ID, steps, make_answer(statement(name, decision, parts), facts=facts), [sheet], _warnings(q, steps, got, decision))


NOISE = ("Input neutralized", "Permanent charge")        # engine notes about the preparation of the structure: normal, not worth a warning here


def _warnings(q, steps, got, decision):
    """The result's warnings: the form's, the failed steps, the engines' own (the same remark from two engines once; the NEIMS ones only on the GC
    route, where EI matters), a failed m/z check, a structure the rules do not fit (no carbon, isotope labels) and a pH that does not match the additive."""
    lc_route = decision.ms_route != rules.ROUTE_GC
    out = list(q.warnings)
    out += [f"Step {r['id']} ({r['title']}) failed: {r['note']}" for r in steps.failed()]
    for n in (3, 4, 5, 6):
        for w in (got.get(n) or {}).get("warnings", []):
            if not w.startswith(NOISE) and not (lc_route and "NEIMS" in w) and not any(w[:30] == seen[:30] for seen in out):
                out.append(w)
    if (got.get(9) or {}).get("warning"):
        out.append(got[9]["warning"])
    out += _structure_warnings(q, got, lc_route)
    if lc_route:
        out += _ph_warnings(q)
    return out


def _structure_warnings(q, got, lc_route):
    """Remarks on structures the rules are not written for: the analysed part has no carbon (NaCl is analysed as HCl, which GC-MS can see, but the
    salt is not measured that way), or the structure carries isotope labels (the ESI ion calculation ignores them)."""
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    out = []
    analysed = Chem.MolFromSmiles(((got.get(2) or {}).get("main_smiles")) or q.smiles)
    if analysed is not None and not any(a.GetSymbol() == "C" for a in analysed.GetAtoms()):
        out.append(f"The analysed part ({rdMolDescriptors.CalcMolFormula(analysed)}) contains no carbon. The rules are written for organic compounds; "
                   "inorganic ions are usually measured by ion chromatography or ICP-MS.")
    given = Chem.MolFromSmiles(q.smiles)
    if lc_route and given is not None and any(a.GetIsotope() for a in given.GetAtoms()):
        out.append("The structure has isotope labels. The ESI m/z is calculated for the unlabelled compound; add the label mass to it.")
    return out


def _ph_warnings(q):
    """The pH used and the additive chosen disagree (more than PH_TOLERANCE): the sheet lists both as given, so say so."""
    from mzdesk.platform_fit import adduct_rules

    assumed = adduct_rules.additive(q.additive).ph
    if abs(assumed - q.ph) <= rules.PH_TOLERANCE:
        return []
    this = "no additive (unbuffered)" if q.additive == "None" else q.additive.lower()
    if q.ph_given:
        return [f"pH {q.ph:g} does not match {this}, which gives about pH {assumed:g}. Choose the additive that gives pH {q.ph:g}, or change the pH."]
    return [f"No pH was given, so pH {q.ph:g} was used for the separation and ionization rules, but {this} gives about pH {assumed:g}. Enter the pH to match."]
