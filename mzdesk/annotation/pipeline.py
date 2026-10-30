"""The ESI MS/MS Annotation workflow: runs the 8 steps in order and builds the Result JSON.

  1 adduct  2 formula   run_adduct_and_formula (BUDDY, one call for both)
  3 library             run_library; the Level-2 gate
  4 MS-FINDER           only if the gate failed
  5 FRIGID              only if the gate failed and MS-FINDER was not good enough (positive mode)
  6 merge + re-rank     every layer's candidates, the sub-threshold library hits included; the layers' scores decide the order,
                        detectability only breaks near-ties (rerank.py)
  7 isomers             one MS-FINDER run for the top-1 and its same-formula candidates; information only, never changes a level
  8 answer              answer.build_answer

A library hit that passes the gate skips steps 4-5; steps 6-7 still run on the hits and the result stays Level 2. Every step may fail
without stopping the run (status "failed" + a note); the result always has an answer (formula_only / mass_only / none at worst).
The engines run one after another, each under the heavy-job lock of its own engine module.

run(inputs, progress=None) -> Result:
  {"answer": answer.build_answer(), "candidates": [Candidate.to_dict()] (ranked, with isomer annotations; formula rows when no structure),
   "steps": [{"id", "title", "status": done | skipped | failed, "seconds", "note"}], "warnings": [str],
   "context": {"precursor_mz", "polarity", "adduct", "retained_adducts", "omics", "libraries", "gate_passed"},
   "isomer": the dict of isomers.discriminate() without its engine log, or None when step 7 did not run}
progress(step_id, status, note, blocks) is called when a step starts (status "running"), when its note changes ("scored 7/25") and when it
ends. status is running | done | skipped | failed; note and blocks are None when unchanged. A block is
{"type": "table" | "notes", "id", "title", "columns"?, "data": {...}}: the spec block and its data (mzdesk/api/blocks.py shapes) together.
parse_inputs(inputs) -> Query and run_query(query, progress) are the two halves of run(): the API parses first, so a bad form is a 400.
"""
from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field

from mzdesk import adducts
from mzdesk.annotation import answer as ans
from mzdesk.annotation import formula_step, frigid_step, library_step, msfinder_step
from mzdesk.annotation.levels import level_name
from mzdesk.annotation.merge import formula_only_candidates, merge_candidates
from mzdesk.annotation.rerank import DEFAULT_PH, rerank
from mzdesk.annotation.schema import LAYER_LABELS, STATUS_DONE, STATUS_FAILED, STATUS_SKIPPED
from mzdesk.api import blocks as blk
from mzdesk.api.common import InputError, flag, fmt, number, precursor, spectrum, text
from mzdesk.exclusivity import isomers
from mzdesk.shared.spectrum import parse_spectrum
from mzdesk.spec import LIBRARIES, PIPELINES

log = logging.getLogger("mzdesk.annotation")

PIPELINE_ID = "esi_ms2_annotation"
STEP_TITLES = {s["id"]: s["title"] for s in next(p for p in PIPELINES if p["id"] == PIPELINE_ID)["stages"]}
OMICS_CHOICES = ("Either / unknown", "Metabolomics", "Lipidomics")
ALKALI_CHOICES = {"unknown": None, "present": True, "absent": False}
MAX_ENGINE_PEAKS = 100          # most intense MS2 peaks sent to every engine step (1-2, 4, 5, 7): MS-FINDER passed 4 GB on 141 peaks. The library search sees all peaks
MAX_INPUT_PEAKS = 5000          # a longer MS2 spectrum is rejected (400): it is not a cleaned MS2 spectrum
MSFINDER_TOP_N = 10
FRIGID_TOP_N = 5
LEVEL_TONES = {"2a": "good", "3": "neutral", "H": "warn", "4": "neutral", "5": "neutral"}


# ---------------------------------------------------------------- inputs
@dataclass
class Query:
    """The form, parsed. adduct None = Auto; ms1 may be empty; alkali None = unknown."""
    precursor_mz: float
    positive: bool
    adduct: str | None
    ms2: list
    ms1: list
    libraries: list
    tol_ppm: float
    tol_da: float
    omics: str
    additive: str
    alkali: bool | None
    ph: float
    halogen: bool
    local_isomers: bool
    warnings: list = field(default_factory=list)

    @property
    def polarity(self):
        return "Positive" if self.positive else "Negative"


def parse_inputs(inputs):
    """Form JSON -> Query. Raises InputError (HTTP 400) with a message for the user."""
    from mzdesk.platform_fit import adduct_rules

    warnings = []
    mz, positive, adduct = precursor(inputs)
    ms2, header_mz = spectrum(inputs, "ms2", "MS2 spectrum", warnings)
    if len(ms2) > MAX_INPUT_PEAKS:
        raise InputError(f"The MS2 spectrum has {len(ms2)} peaks; at most {MAX_INPUT_PEAKS} are accepted. Clean it first (Spectrum cleanup page) or keep the strongest peaks.")
    mz = mz if mz is not None else header_mz
    if mz is None:
        raise InputError("Precursor m/z is required (or put it in the header of the MS2 spectrum).")
    if adduct is not None:
        if not adducts.canonical(adduct):
            raise InputError(f"Unknown adduct {adduct!r}.")
        if adducts.polarity(adduct) != ("Positive" if positive else "Negative"):
            raise InputError(f"{adduct} is not a {'positive' if positive else 'negative'}-ion adduct.")
    libraries = inputs.get("libraries") or []
    if not isinstance(libraries, list) or any(name not in LIBRARIES for name in libraries):
        raise InputError(f"Libraries must be a list of: {', '.join(LIBRARIES)}.")
    omics = text(inputs, "omics") or OMICS_CHOICES[0]
    if omics not in OMICS_CHOICES:
        raise InputError(f"Omics must be one of: {', '.join(OMICS_CHOICES)}.")
    additive = text(inputs, "additive") or "Unknown"
    if additive != "Unknown" and additive not in adduct_rules.ADDITIVES:
        raise InputError(f"Mobile phase additive must be Unknown or one of: {', '.join(adduct_rules.ADDITIVES)}.")
    alkali = text(inputs, "alkali").lower() or "unknown"
    if alkali not in ALKALI_CHOICES:
        raise InputError("Sodium / potassium must be Unknown, Present or Absent.")
    return Query(
        precursor_mz=mz, positive=positive, adduct=adduct, ms2=ms2, ms1=parse_spectrum(text(inputs, "ms1")), libraries=list(dict.fromkeys(libraries)),
        tol_ppm=number(inputs, "tol_ppm", "Precursor tolerance (ppm)", default=10.0, positive=True),
        tol_da=number(inputs, "tol_da", "Fragment tolerance (Da)", default=0.02, positive=True),
        omics=omics, additive=additive, alkali=ALKALI_CHOICES[alkali],
        ph=number(inputs, "ph", "Mobile phase pH", default=DEFAULT_PH, minimum=0, maximum=14),
        halogen=flag(inputs, "halogen", "Allow Cl/Br in formulas"),
        local_isomers=flag(inputs, "local_isomers", "Also compare with database isomers"), warnings=warnings)


# ---------------------------------------------------------------- blocks (shown while the job runs)
def _table(block_id, title, columns, rows):
    return {"type": "table", "id": block_id, "title": title, "columns": columns, "data": blk.table_data(rows, count=len(rows))}


def _notes(block_id, title, items):
    return {"type": "notes", "id": block_id, "title": title, "data": {"items": list(items)}}


def _level_cell(cand):
    return blk.badge(level_name(cand.level), LEVEL_TONES.get(cand.level, "neutral")) if cand.level else "—"


def _adduct_blocks(res):
    rows = []
    for h in res.hypotheses:
        ev = h.evidence or {}
        found = [s["label"] for s in (ev.get("ms2") or {}).get("signals", []) if s.get("found")]
        kept = "no" if not h.retained else "yes (ambiguous)" if h.ambiguous else "yes"
        rows.append([h.rank, h.adduct, fmt(h.score, 3), kept, ", ".join(found) or "—",
                     (ev.get("formula") or {}).get("formula") or "—", fmt((ev.get("isotope") or {}).get("similarity"), 2)])
    return [_table("adduct_table", "Adduct hypotheses", ["Rank", "Adduct", "Score", "Kept", "MS2 signals", "Best formula", "Isotope fit"], rows)]


def _formula_blocks(res):
    rows = [[f.adduct, f.formula, f.rank, fmt(f.prob, 3), fmt(f.ppm, 2), fmt(f.isotope_sim, 2)] for f in res.formulas]
    return [_table("formula_table", "Formulas", ["Adduct", "Formula", "Rank", "BUDDY probability", "Error (ppm)", "Isotope fit"], rows)] if rows else []


def _library_blocks(layer, gate):
    out, rows = [], []
    for c in layer.candidates:
        ev = c.evidence["library"]
        rows.append([c.name or "—", c.smiles or "—", c.adduct or "—", fmt(ev["similarity"], 3), ev["matched_peaks"], fmt(ev.get("precursor_error_ppm"), 2),
                     blk.badge("passed", "good") if ev["gate_passed"] else blk.badge("failed", "neutral"), "; ".join(ev["failed_checks"]) or "—"])
    if rows:
        out.append(_table("library_table", "Library hits",
                          ["Name", "SMILES", "Adduct", "Similarity", "Matched peaks", "Precursor error (ppm)", "Gate", "Failed checks"], rows))
    if gate is not None and gate.checks:
        out.append(_table("library_gate", "Level-2 gate (best hit)", ["Check", "Result", "Detail"], [[c["label"], c["status"], c["detail"]] for c in gate.checks]))
    return out


def _engine_blocks(prefix, title, layer):
    rows = [[i, c.name or "—", c.smiles or "—", c.formula or "—", c.adduct or "—", fmt(c.layer_score, 3), fmt(c.layer_score_norm, 2)]
            for i, c in enumerate(layer.candidates, start=1)]
    return [_table(f"{prefix}_table", title, ["Rank", "Name", "SMILES", "Formula", "Adduct", "Score", "Score (0-1)"], rows)] if rows else []


def _ranked_blocks(ranked):
    rows = []
    for c in ranked:
        d = c.detect or {}
        rows.append([c.rank, c.name or "—", c.formula or "—", c.adduct or "—", LAYER_LABELS[c.layer], _level_cell(c), fmt(c.base_score, 3),
                     fmt(c.final_score, 3), d.get("omics", "—"), d.get("ionizable", "—"), d.get("adduct", "—")])
    columns = ["Rank", "Name", "Formula", "Adduct", "Layer", "Level", "Base score", "Final score", "Omics", "Ionizable", "Adduct fit"]
    if not rows:
        return []
    note = ("Candidates within 10% of the best base score are ordered by base score x detectability; the others follow in base-score order, "
            "so the Final score column does not always go down the list.")
    return [_table("ranked_table", "Ranked candidates", columns, rows), _notes("ranked_note", "How the rank is set", [note])]


def _isomer_blocks(res):
    by_id = {i["id"]: i for i in res["isomers"]}
    challengers = {ch["id"] for ch in res["challengers"]}
    rows = []
    for c in res["candidates"]:
        iso = by_id.get(c["id"])
        if c["id"] == res["top_id"]:
            role = "top-1"
        elif c["id"] in res["group"]:
            role = "not separated"
        else:
            role = "separated" if iso and iso["distinguished"] else "not scored"
        rows.append([c["name"] or c["smiles"], c["smiles"], c["source"], role + (" · challenger" if c["id"] in challengers else ""),
                     fmt(c["explained_fraction"], 2), len(iso["distinguishing_peaks"]) if iso else "—"])
    items = list(res["notes"])
    for ch in reversed(res["challengers"]):
        mzs = ", ".join(f"{p['mz']:.4f}" for p in ch["peaks"])
        items.insert(0, f"Another isomer explains peaks the top-1 cannot: {ch['name'] or ch['id']} (m/z {mzs}).")
    blocks = [_table("isomer_table", "Isomers compared", ["Candidate", "SMILES", "Source", "Role", "Explained intensity", "Peaks only the top-1 explains"], rows)]
    return blocks + [_notes("isomer_notes", "Notes", items)]


# ---------------------------------------------------------------- the run
class _Steps:
    """Per-step bookkeeping and the live progress calls."""

    def __init__(self, progress):
        self.progress = progress
        self.records = {n: {"id": n, "title": title, "status": "queued", "seconds": 0.0, "note": ""} for n, title in STEP_TITLES.items()}
        self.started = {}

    def emit(self, n, status, note=None, blocks=None):
        if self.progress:
            try:
                self.progress(n, status, note, blocks)
            except Exception:           # noqa: BLE001 - a broken listener must not stop the run
                log.exception("progress callback failed")

    def start(self, n):
        self.started[n] = time.perf_counter()
        self.records[n]["status"] = "running"
        self.emit(n, "running")

    def finish(self, n, status, note, blocks=(), seconds=None):
        if seconds is None:
            seconds = time.perf_counter() - self.started.get(n, time.perf_counter())
        self.records[n].update(status=status, note=note, seconds=round(seconds, 3))
        self.emit(n, status, note, list(blocks))

    def skip(self, n, note):
        self.records[n].update(status=STATUS_SKIPPED, note=note)
        self.emit(n, STATUS_SKIPPED, note, [])

    def guarded(self, n, fn):
        """Run fn() -> (status, note, blocks, value) as step n; an exception makes the step "failed". Returns the value (None on failure)."""
        self.start(n)
        try:
            status, note, blocks, value = fn()
        except Exception as exc:        # noqa: BLE001 - a layer must never take the workflow down
            log.exception("step %s failed", n)
            status, note, blocks, value = STATUS_FAILED, f"{STEP_TITLES[n]} failed: {type(exc).__name__}: {exc}", [], None
        self.finish(n, status, note, blocks)
        return value


def run(inputs, progress=None):
    """Form JSON -> Result (module docstring). Raises InputError for a bad form; nothing else escapes."""
    return run_query(parse_inputs(inputs), progress)


def run_query(q, progress=None):
    steps = _Steps(progress)
    warnings = list(q.warnings)
    ms2_engine = msfinder_step.clean_peaks(q.ms2, MAX_ENGINE_PEAKS)         # what the engine steps see; q.ms2 stays whole for the library search
    if len(ms2_engine) < len(q.ms2):
        warnings.append(f"The MS2 spectrum has {len(q.ms2)} peaks; the engine steps (formula, MS-FINDER, FRIGID, isomers) used the {len(ms2_engine)} most intense.")

    # ---- steps 1 + 2 (one BUDDY pass serves both)
    steps.start(1)
    steps.start(2)
    adduct_res = formula_res = None
    try:
        adduct_res, formula_res = formula_step.run_adduct_and_formula(
            q.precursor_mz, q.positive, ms2_engine, q.ms1 or None, adduct=q.adduct, tol_ppm=q.tol_ppm, tol_da=q.tol_da, halogen=q.halogen)
    except Exception as exc:            # noqa: BLE001
        log.exception("steps 1-2 failed")
        steps.finish(1, STATUS_FAILED, f"Adduct step failed: {type(exc).__name__}: {exc}")
        steps.finish(2, STATUS_FAILED, "No formulas: step 1 failed.")
    if adduct_res is not None:
        steps.finish(1, adduct_res.status, adduct_res.note, _adduct_blocks(adduct_res), adduct_res.seconds)
        steps.finish(2, formula_res.status, formula_res.note, _formula_blocks(formula_res), formula_res.seconds)
    retained = (adduct_res.retained_adducts if adduct_res is not None else []) or [q.adduct or ("[M+H]+" if q.positive else "[M-H]-")]
    top_adduct = retained[0]
    formulas = formula_res.formulas if formula_res is not None else []
    formula_dicts = [f.to_dict() for f in formulas]

    # ---- step 3
    layers, gate = [], None
    if not q.libraries:
        steps.skip(3, "No library was selected, so no library spectrum was searched.")
    else:
        def library():
            layer, summary = library_step.run_library(q.ms2, q.precursor_mz, q.positive, q.libraries, adducts=retained, formulas=formulas,
                                                      ms1_peaks=q.ms1 or None, tol_ppm=q.tol_ppm, tol_da=q.tol_da, halogen=q.halogen)
            return layer.status, layer.note, _library_blocks(layer, summary), (layer, summary)
        got = steps.guarded(3, library)
        if got:
            layers.append(got[0])
            gate = got[1]
            warnings += gate.warnings
    gate_passed = bool(gate and gate.passed)

    # ---- steps 4 + 5
    if gate_passed:
        steps.skip(4, "Library gate passed.")
        steps.skip(5, "Library gate passed.")
    else:
        def msfinder():
            layer, summary = msfinder_step.run_msfinder_layer(q.precursor_mz, top_adduct, ms2_engine, formulas, tol_ppm=q.tol_ppm, tol_da=q.tol_da,
                                                              top_n=MSFINDER_TOP_N, ms1_peaks=q.ms1 or None)
            return layer.status, layer.note, _engine_blocks("msfinder", "MS-FINDER structures", layer), (layer, summary)
        got, good_enough = steps.guarded(4, msfinder), False
        if got:
            layers.append(got[0])
            good_enough = bool(got[1].get("good_enough"))
            warnings += got[1].get("warnings", [])
        if good_enough:
            steps.skip(5, "MS-FINDER result good enough.")
        else:
            formula = formulas[0].formula if formulas else None
            adduct5 = formulas[0].adduct if formulas else top_adduct

            def frigid():
                layer = frigid_step.run_frigid_layer(q.precursor_mz, adduct5, ms2_engine, formula, top_n=FRIGID_TOP_N)
                return layer.status, layer.note, _engine_blocks("frigid", "FRIGID structures (hypotheses)", layer), layer
            layer = steps.guarded(5, frigid)
            if layer:
                layers.append(layer)
    for n in (3, 4, 5):
        if steps.records[n]["status"] == STATUS_FAILED:
            warnings.append(f"Step {n} ({STEP_TITLES[n]}) failed: {steps.records[n]['note']}")

    # ---- step 6
    def merge_rerank():
        merged = merge_candidates(layers)
        if not merged:
            rows = formula_only_candidates(formula_dicts)
            return STATUS_DONE, ("No structure was proposed; the formulas are listed instead." if rows else "No structure and no formula."), [], ([], rows)

        def scored(done, total):
            steps.emit(6, "running", f"scored {done}/{total}")
        try:
            ranked = rerank(merged, polarity=q.polarity, adduct=top_adduct, omics=q.omics, ph=q.ph, additive=q.additive, alkali=q.alkali, progress=scored)
            n_layers = len({s["layer"] for c in merged for s in c.sources})
            note = f"{len(ranked)} structure{'' if len(ranked) == 1 else 's'} from {n_layers} layer{'' if n_layers == 1 else 's'}, in the layers' order; detectability only orders near-ties."
            status = STATUS_DONE
        except Exception as exc:        # noqa: BLE001 - keep the merged list, ordered by the layers' own scores
            log.exception("re-rank failed")
            ranked = sorted(merged, key=lambda c: -c.layer_score_norm)
            for i, c in enumerate(ranked, start=1):
                c.rank = i
            status, note = STATUS_FAILED, f"Re-rank failed ({type(exc).__name__}: {exc}); candidates are ordered by the layer scores."
            warnings.append(note)
        return status, note, _ranked_blocks(ranked), (ranked, [])
    got = steps.guarded(6, merge_rerank)
    ranked, formula_rows = got if got else ([], formula_only_candidates(formula_dicts))
    if gate_passed and ranked and ranked[0].level != "2a":
        warnings.append("The library gate passed, but re-ranking put a candidate that did not pass the gate first; the answer follows the ranking.")

    # ---- step 7
    iso = None
    if not ranked:
        steps.skip(7, "No structure candidate to compare.")
    else:
        def isomer_check():
            top = ranked[0]
            if not top.smiles:
                return STATUS_SKIPPED, "The top-ranked candidate has no structure, so isomers cannot be compared.", [], None
            pool = [{"id": c.id, "smiles": c.smiles, "name": c.name} for c in ranked
                    if c.smiles and (c is top or library_step.same_formula(c.formula, top.formula))]
            res = isomers.discriminate(pool, ms2_engine, q.precursor_mz, top.adduct or top_adduct, tol_da=q.tol_da, tol_ppm=q.tol_ppm,
                                       engine="msfinder", local_isomers=q.local_isomers)
            if res["verdict"] == isomers.INSUFFICIENT:
                status = STATUS_FAILED if res.get("error") else STATUS_SKIPPED
                return status, (res["notes"] or ["The isomer check could not be done."])[0], [_notes("isomer_notes", "Notes", res["notes"])], res
            note = ("Fragment peaks tell the top-1 apart from every isomer compared." if res["verdict"] == isomers.SINGLE
                    else f"{len(res['group'])} structures cannot be separated by the spectrum.")
            if res["challengers"]:
                note += " Another isomer explains peaks the top-1 cannot."
            return STATUS_DONE, note, _isomer_blocks(res), res
        iso = steps.guarded(7, isomer_check)
        if iso is None:
            iso = {"verdict": isomers.INSUFFICIENT, "top_id": ranked[0].id, "group": [], "isomers": [], "challengers": [], "candidates": [],
                   "notes": [steps.records[7]["note"]]}
        ans.attach_isomer(ranked, iso)

    # ---- step 8
    def build():
        answer = ans.build_answer(ranked, iso, formula_list=formula_dicts, precursor_mz=q.precursor_mz)
        note = f"{answer['level_label'] or 'No level'} ({answer['type'].replace('_', ' ')})."
        return STATUS_DONE, note, [_notes("answer_notes", "Answer", [answer["statement"]])], answer
    answer = steps.guarded(8, build)
    if answer is None:                  # build_answer itself failed: the result still carries an answer
        answer = {"type": ans.ANSWER_NONE, "ids": [], "layer": "", "layer_label": "", "level": "", "level_label": "", "top_id": None,
                  "isomer_verdict": None, "statement": "The answer could not be built: " + steps.records[8]["note"]}

    # Display only: these figures use the parsed query and the best judged library hit, never the editable form.
    try:
        from mzdesk.annotation.visualization import spectrum_figures

        figures = spectrum_figures(q.ms2, q.precursor_mz, q.tol_da, gate.reference if gate else None)
    except Exception:  # noqa: BLE001 - a display failure must not affect the annotation
        log.exception("workflow figures failed")
        figures = None

    return {
        "answer": answer,
        "candidates": [c.to_dict() for c in (ranked or formula_rows)],
        "steps": [steps.records[n] for n in sorted(steps.records)],
        "warnings": warnings,
        "context": {"precursor_mz": q.precursor_mz, "polarity": q.polarity, "adduct": top_adduct, "retained_adducts": retained,
                    "omics": q.omics, "libraries": q.libraries, "gate_passed": gate_passed},
        "isomer": {k: v for k, v in iso.items() if k != "log"} if iso else None,
        "visualization": figures,
    }
