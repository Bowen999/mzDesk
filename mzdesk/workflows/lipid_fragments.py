"""Lipid Fragment Explanation workflow: a lipid (shorthand name, SMILES or drawing) and its MS2 spectrum in; names, a structure and the explained
peaks out.

Steps (spec.py stages): 1 Lipid name  2 Structure  3 Fragment explanation.   Result contract and the step runner: mzdesk/workflows/common.py.

The steps only chain two existing pages:
  1 + 2  Lipid Short Name (mzdesk/api/runners/lipid.py): LIPID+ L1 / L2 / L3, LIPID MAPS species / molecular species and the illustrative structure.
         A name only at the sum level (PC 34:1) is split into chains by the page's default rule; the drawing marks every assumed part.
  3      Fragment Explanation (mzdesk/api/runners/ms2.py run_fragments, MetFrag) on that structure, tree depth 3 by default.
A SMILES or a drawing gives no names (step 1 is skipped and the name cards show "—"); its structure is used as given.

Auto adduct: with a precursor m/z (typed or read from the spectrum header), the MetFrag adduct of the chosen polarity whose theoretical m/z is
closest within the precursor tolerance; without one, or when none is close, [M+H]+ / [M-H]-.

Form field ids: lipid ({"mode": LIPID_MODE | "Draw" | "SMILES / InChI", "value"}), precursor, ms2, tree_depth, tol_abs, tol_ppm.
"""
from __future__ import annotations

import copy
import logging
from dataclasses import dataclass, field

from mzdesk.api import common
from mzdesk.api.common import InputError
from mzdesk.spec import LIPID_NAME_MODE, MODULES, PIPELINES
from mzdesk.workflows.common import STATUS_DONE, STATUS_SKIPPED, StepRunner, block, make_answer, make_result, notes_block

log = logging.getLogger("mzdesk.workflows")

PIPELINE_ID = "lipid_fragments"
STEP_TITLES = {s["id"]: s["title"] for s in next(p for p in PIPELINES if p["id"] == PIPELINE_ID)["stages"]}
LIPID_MODE = LIPID_NAME_MODE
DEFAULT_DEPTH = 3
DEPTH_MIN, DEPTH_MAX = 1, 4
DEFAULT_TOL_ABS, DEFAULT_TOL_PPM = 0.002, 10.0
PRECURSOR_TOL_DA, PRECURSOR_TOL_PPM = 0.01, 10.0     # Auto adduct: the theoretical m/z must be this close to the precursor (the larger of the two)
NO_NAME_HELP = "Names come only from a lipid short name."
NO_STRUCTURE_HINT = "Enter the structure as SMILES or draw it to explain the fragments."


@dataclass
class Query:
    """The form, parsed. name: the lipid shorthand ("" for a structure); lipid: its normalize() result (None for a structure); smiles: the structure
    given as SMILES or drawn ("" for a name). precursor_mz: typed, else read from the spectrum header, else None (theoretical). adduct: chosen, or
    None for Auto. ms2: the spectrum field as sent, handed on to Fragment Explanation."""
    name: str
    lipid: dict | None
    smiles: str
    precursor_mz: float | None
    positive: bool
    adduct: str | None
    ms2: dict
    peaks: list
    depth: int = DEFAULT_DEPTH
    tol_abs: float = DEFAULT_TOL_ABS
    tol_ppm: float = DEFAULT_TOL_PPM
    mode: str = LIPID_MODE
    warnings: list = field(default_factory=list)


def parse_inputs(inputs):
    """Form JSON -> Query. Raises InputError (HTTP 400): no lipid, an unreadable name (the Lipid Short Name page's sentence), a bad SMILES,
    no peaks, an adduct MetFrag does not have or that does not fit the ion mode, numbers out of range."""
    from mzdesk.api.runners import lipid as lipid_page
    from mzdesk.ms2 import metfrag_engine as mf

    warnings = []
    value = inputs.get("lipid")
    if value is not None and not isinstance(value, dict):
        raise InputError("Lipid: invalid input; expected an object with mode and value.")
    value = value or {}
    mode = value.get("mode") or LIPID_MODE
    name, lipid, smiles = "", None, ""
    if mode == LIPID_MODE:
        name = str(value.get("value") or "").strip()
        if not name:
            raise InputError("Lipid is required.")
        lipid = lipid_page.read_name(name)
    elif mode in ("SMILES / InChI", "Draw"):
        smiles = common.structure({"lipid": value}, "lipid", "Lipid", warnings=warnings)["smiles"]
    else:
        raise InputError(f"Lipid: unknown input mode {mode!r}; use {LIPID_MODE}, Draw or SMILES / InChI.")

    precursor_mz, positive, adduct = common.precursor(inputs)
    peaks, file_mz = common.spectrum(inputs, "ms2", "MS2 spectrum", warnings)
    if precursor_mz is None and file_mz:
        precursor_mz = file_mz
        warnings.append(f"Precursor m/z {file_mz:.4f} is read from the spectrum.")
    elif precursor_mz is not None and file_mz and abs(precursor_mz - file_mz) > 0.01:
        warnings.append(f"The MS2 file says precursor m/z {file_mz:.4f}; the form value {precursor_mz:.4f} is used.")
    if adduct is not None:
        if adduct not in mf.ADDUCTS:
            raise InputError(f"MetFrag does not support {adduct}. Supported: {', '.join(mf.ADDUCTS)}.")
        if mf.ADDUCTS[adduct]["positive"] != positive:
            raise InputError(f"{adduct} does not match the ion mode.")
    depth = int(common.number(inputs, "tree_depth", "Tree depth", default=DEFAULT_DEPTH, minimum=DEPTH_MIN, maximum=DEPTH_MAX))
    tol_abs = common.number(inputs, "tol_abs", "Absolute tolerance", default=DEFAULT_TOL_ABS, minimum=0)
    tol_ppm = common.number(inputs, "tol_ppm", "Relative tolerance", default=DEFAULT_TOL_PPM, minimum=0)
    # the header precursor is already in precursor_mz: hand the peaks on as plain text so Fragment Explanation does not say it twice
    ms2 = {"source": "Paste", "value": "\n".join(f"{mz:.6f} {inten:g}" for mz, inten in peaks)}
    return Query(name=name, lipid=lipid, smiles=smiles, precursor_mz=precursor_mz, positive=positive, adduct=adduct, ms2=ms2, peaks=peaks,
                 depth=depth, tol_abs=tol_abs, tol_ppm=tol_ppm, mode=mode, warnings=warnings)


# ---------------------------------------------------------------- helpers
def _module_block(module_id, block_id, data, title=None, **extra):
    """The spec block `block_id` of a module page with the data its runner returned; the id stays (one workflow, one block per id)."""
    module = next(m for m in MODULES if m["id"] == module_id)
    spec_block = copy.deepcopy(next(b for b in module["outputs"] if b.get("id") == block_id))
    if title is not None:
        spec_block["title" if spec_block["type"] != "figure" else "label"] = title
    return {**spec_block, **extra, "data": data}


def pick_adduct(smiles, positive, precursor_mz):
    """Auto adduct -> (adduct, note or None). The MetFrag adduct of this polarity whose theoretical m/z is closest to the precursor, when it is
    within PRECURSOR_TOL_DA or PRECURSOR_TOL_PPM; else [M+H]+ / [M-H]- with a note. [M]+ is left out: it is for permanently charged structures."""
    from mzdesk.ms2 import metfrag_engine as mf

    default = "[M+H]+" if positive else "[M-H]-"
    if precursor_mz is None:
        return default, None
    tol = max(PRECURSOR_TOL_DA, precursor_mz * PRECURSOR_TOL_PPM * 1e-6)
    best = None
    for name, add in mf.ADDUCTS.items():
        if add["positive"] != positive or name == "[M]+":
            continue
        info = mf.precursor_info(smiles, name)
        if info.get("error"):
            return default, None            # Fragment Explanation reports the structure problem itself
        diff = abs(info["theo_mz"] - precursor_mz)
        if best is None or diff < best[1]:
            best = (name, diff)
    if best and best[1] <= tol:
        return best[0], None
    theo = mf.precursor_info(smiles, default)["theo_mz"]
    return default, (f"Precursor m/z {precursor_mz:.4f} fits none of the {'positive' if positive else 'negative'} adducts of this structure "
                     f"({default} is {theo:.4f}); {default} was used. Choose the adduct, or check the precursor m/z and the lipid.")


# ---------------------------------------------------------------- step 1: Lipid name
def step_names(q):
    from mzdesk.api.runners import lipid as lipid_page

    cards = lipid_page.name_cards(q.lipid)
    flags = lipid_page.flag_lines(q.lipid, q.name)
    blocks = [_module_block("lipid_short_name", "lipid_plus", cards["lipid_plus"]),
              _module_block("lipid_short_name", "lipid_maps", cards["lipid_maps"])]
    plus = " · ".join(dict.fromkeys(v for v in (q.lipid.get("L1"), q.lipid.get("L2"), q.lipid.get("L3")) if v))
    return STATUS_DONE, f"LIPID+ {plus}", blocks, {"cards": cards, "flags": flags}


def no_name_cards():
    """The name cards for a structure input: every value "—" with the reason in its "?"."""
    dash = {"value": common.MISSING, "help": NO_NAME_HELP}
    return {"lipid_plus": {"values": {"L1": dash, "L2": dash, "L3": dash}},
            "lipid_maps": {"values": {"Species": dash, "Molecular species": dash}}}


# ---------------------------------------------------------------- step 2: Structure
def step_structure(q):
    """-> {"smiles", "svg", "assumptions", "illustrative": bool}; a name the templates cannot draw skips the step with the page's sentence."""
    from mzdesk.api.runners import lipid as lipid_page

    if q.lipid is None:
        svg = common.depiction_svg(q.smiles, size=(720, 300))
        got = {"smiles": q.smiles, "svg": svg, "assumptions": [], "illustrative": False}
        blocks = [structure_block(got)] if svg else []
        return STATUS_DONE, "Structure as given.", blocks, got
    ill, svg, assumptions, note = lipid_page.illustrate(q.lipid, q.name)
    if ill is None:
        return STATUS_SKIPPED, f"{note} {NO_STRUCTURE_HINT}", [], None
    # the shaded parts of the drawing, in words (the drawing's legend says the same)
    shaded = [text for atoms, text in ((ill.split_atoms, "Chains split by the default rule"), (ill.db_atoms, "C=C positions assumed"),
                                       (ill.fg_atoms, "Group positions assumed")) if atoms]
    assumptions = shaded + [a for a in assumptions if a not in shaded]
    got = {"smiles": ill.smiles, "svg": svg, "assumptions": assumptions, "illustrative": True}
    said = "; ".join(a[0].lower() + a[1:] if a[1:2].islower() else a for a in assumptions)      # "Chains …" -> "chains …", "C=C …" stays
    return STATUS_DONE, f"{ill.formula}, illustrative{': ' + said if said else ''}.", [structure_block(got)], got


def structure_block(got):
    title = "Illustrative structure" if got["illustrative"] else "Structure"
    help = (["One structure that fits the name, built from templates.", "Shaded parts are assumed: the name does not give them."]
            if got["illustrative"] else None)
    return _module_block("lipid_short_name", "structure", {"svg": got["svg"]}, title=title, **({"help": help} if help else {}))


# ---------------------------------------------------------------- step 3: Fragment explanation
def step_fragments(q, smiles):
    from mzdesk.api.runners import ms2

    adduct, adduct_note = (q.adduct, None) if q.adduct else pick_adduct(smiles, q.positive, q.precursor_mz)
    form = {"structure": {"mode": "SMILES / InChI", "value": smiles},
            "precursor": {"precursor_mz": q.precursor_mz, "ion_mode": "Positive" if q.positive else "Negative", "adduct": adduct},
            "ms2": q.ms2, "tree_depth": q.depth, "tol_abs": q.tol_abs, "tol_ppm": q.tol_ppm}
    out = ms2.run_fragments(form)
    o = out["outputs"]
    blocks = [_module_block("ms2_fragments", "summary", o["summary"], title="Fragment explanation"),
              _module_block("ms2_fragments", "spectrum", o["spectrum"]),
              _module_block("ms2_fragments", "peaks", o["peaks"], title="Peaks")]
    values = o["summary"]["values"]
    note = f"{values['Explained peaks']} peaks explained ({values['Explained intensity']} % of the intensity) as {adduct}, tree depth {q.depth}."
    warnings = [adduct_note] if adduct_note else []
    warnings += [w for w in out["warnings"] if not (adduct_note and w.startswith("Precursor m/z") and "does not match" in w)]
    return STATUS_DONE, note, blocks, {"adduct": adduct, "auto": q.adduct is None, "values": values, "warnings": warnings, "blocks": blocks}


# ---------------------------------------------------------------- run
def run_query(q, progress=None):
    """The three steps in order. A name the templates cannot draw skips steps 2 and 3 with the reason; the names still show."""
    steps = StepRunner(STEP_TITLES, progress)
    got = {}
    if q.lipid is not None:
        got[1] = steps.guarded(1, lambda: step_names(q))
    else:
        steps.skip(1, NO_NAME_HELP)
    got[2] = steps.guarded(2, lambda: step_structure(q))
    if got[2]:
        got[3] = steps.guarded(3, lambda: step_fragments(q, got[2]["smiles"]))
    else:
        steps.skip(3, "No structure to explain the spectrum with.")

    cards = (got.get(1) or {}).get("cards") or no_name_cards()
    blocks = [_module_block("lipid_short_name", "lipid_plus", cards["lipid_plus"]),
              _module_block("lipid_short_name", "lipid_maps", cards["lipid_maps"])]
    if got[2] and got[2]["svg"]:
        blocks.append(structure_block(got[2]))
    if got.get(3):
        blocks += got[3]["blocks"]
    if got[2] and got[2]["assumptions"]:
        blocks.append(notes_block("assumptions", "Assumptions", got[2]["assumptions"], collapsed=True))
    flags = (got.get(1) or {}).get("flags") or []
    if flags:
        blocks.append(notes_block("flags", "Flags", flags, collapsed=True))
    return make_result(PIPELINE_ID, steps, answer(q, got), blocks, _warnings(q, steps, got))


def answer(q, got):
    label = q.name
    frag, struct = got.get(3), got.get(2)
    facts = [["Lipid", q.name or ("Drawn structure" if q.mode == "Draw" else "SMILES")]]
    if struct and struct.get("smiles"):
        from mzdesk.ms2 import metfrag_engine as mf
        info = mf.precursor_info(struct["smiles"], (frag or {}).get("adduct") or ("[M+H]+" if q.positive else "[M-H]-"))
        if not info.get("error"):
            facts.append(["Formula", info["formula"]])
    if frag:
        facts += [["Adduct", frag["adduct"] + (" (Auto)" if frag["auto"] else "")], ["Tree depth", str(q.depth)]]
        v = frag["values"]
        n, total = (s.strip() for s in v["Explained peaks"].split("/"))
        who = f"The illustrative structure of {label}" if struct and struct["illustrative"] else "The structure"
        statement = f"{who} explains {n} of {total} peaks ({v['Explained intensity']} % of the intensity)."
        level, level_label = None, None
        if struct and struct["illustrative"] and struct["assumptions"]:
            level, level_label = "Assumed parts", "The name does not give every detail; the shaded parts of the structure are assumed."
        return make_answer(statement, level, level_label, facts)
    if struct is None and q.lipid is not None:
        return make_answer(f"{label} has no illustrative structure, so the fragments were not explained. {NO_STRUCTURE_HINT}", facts=facts)
    return make_answer("The fragments could not be explained; see the steps below.", facts=facts)


def _warnings(q, steps, got):
    out = list(q.warnings)
    out += [f"Step {r['id']} ({r['title']}) failed: {r['note']}" for r in steps.failed()]
    for w in (got.get(3) or {}).get("warnings", []):
        if w not in out:
            out.append(w)
    return out
