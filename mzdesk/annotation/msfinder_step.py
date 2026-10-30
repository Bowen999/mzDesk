"""Step 4 of the ESI MS/MS Annotation workflow: structure candidates from MS-FINDER's local databases, fragment-scored, for the shared formulas.

run_msfinder_layer() searches MS-FINDER's local structure databases (PubChem is never queried) for every formula of the shared list of step 2
and scores each structure by in-silico fragmentation against the MS/MS spectrum (MS-FINDER's "Total score"). Candidates are Level 3 (one
per compound, best score first); "good enough" (best layer_score_norm >= GOOD_ENOUGH) means FRIGID (step 5) is not needed.

How the formula list reaches MS-FINDER (decided on real runs: MassBank spectra, BUDDY's top-5 lists, quiet machine, 2026-10-28)
  MS-FINDER cannot take a formula list: it finds its own formulas around the precursor and runs the structure search for the best few.
  (a)  one query per formula, the precursor set to that formula's theoretical m/z with a 0.2 mDa formula window, all queries in ONE process
       (the console app loads its 2.3 GB of databases once and takes every query file of the folder); results filtered to the formula;
  (a2) the same, one process per formula;   (b1) one run at the measured precursor, shipped settings (1 mDa window, 5 formulas);
  (b2) one run with a 5 ppm window and 20 formulas.
                            seconds for the run / BUDDY formulas that got structures
                            (a)         (a2)        (b1)       (b2)
      colchicine, 5 f.      11.8 / 2    27.8 / 2    9.9 / 1    10.3 / 2     (b1 missed C23H21N5O2: not among MS-FINDER's own top 5)
      apigenin-Me, 5 f.      6.4 / 2    33.7 / 2    6.1 / 2     6.2 / 2
      2-OH-estrone, 2 f.     6.9 / 1    10.5 / 1    7.0 / 1     7.1 / 1
      tyrosine, 1 f.         4.5 / 1     4.2 / 1    5.5 / 1     4.9 / 1
  (a) is used: it searches exactly the shared formulas (a run around the precursor depends on MS-FINDER's own formula ranking, and its wider
  window also searches formulas that are not in the list), it covers all of them, and it costs about one run (a2 is 2-5 times slower).
  The target formula ranked first of MS-FINDER's formulas in 70 of 70 queries, so 2 formulas per query are searched. All six CASMI adducts and
  [M+Na]+, [M+NH4]+ came out within 0.01 mDa of the computed m/z, so the 0.2 mDa window is safe. A formula's search time is read from its .sfd
  file (created when the search starts, written when it ends; None where the file system keeps no creation time).
  No formula list (step 2 failed): one run at the measured precursor with MS-FINDER's own formula search (window = tol_ppm), within tol_ppm.

Other settings
  * ms1_peaks is not sent to MS-FINDER when a formula list exists: MS-FINDER folds the isotope fit into the structure score (a wrong pattern
    lowered tyrosine's 9.14 to 8.34 and kept the formula), which would count BUDDY's isotope evidence twice and move the scores off the
    calibrated scale. It is used only in the run without a list.
  * At most MAX_PEAKS (100) most intense MS2 peaks and MAX_FORMULAS (6) formulas. A flavonoid glycoside with 141 peaks passed 4 GB and was
    stopped; with its 100 most intense peaks it finished (4.0 GB, 48 s). Scores moved by 0.0-0.3 for the five calibration spectra above 100 peaks.
  * A structure search starts no new structure after SEARCH_MINUTES and keeps the partial result (reported as "time limit"); structures already
    in flight finish, so a search ends up to ~20 s later. Heavy formulas (C20H30O2: slow, 100 structures) can reach it. The run is stopped
    after ENGINE_TIMEOUT_S or RSS_CAP_MB and keeps the searches that had finished.
  * The shipped parameter file has TreeDepth=2 (owner's design). For large glycosides (CASMI 497, C48H76O20, m/z 971; also 236, 328) the fragment
    enumeration of ONE structure runs for minutes: the time limit above is checked only after a structure is scored, and a formula's .sfd is
    written only when its search ends, so the ENGINE_TIMEOUT_S kill leaves nothing for the unfinished formulas. A run stopped by the TIME limit
    is therefore repeated ONCE for the unfinished formulas only (the whole run without a formula list) with TreeDepth=1 and RETRY_TIMEOUT_S;
    the same 497 search took ~50 s at TreeDepth=1. Memory stops are not repeated: the retry's marker says "timed out", and a memory stop (a
    wide spectrum, not the tree) was not shown to depend on the depth. Retry candidates are marked (formula info "tree_depth"/"note",
    candidate evidence "tree_depth"), because TreeDepth=1 scores are not on the TreeDepth=2 scale.
  * Fragment tolerance = tol_da; the halogen switch follows the formulas (BUDDY has no Si); everything else is the shipped parameter file.
  * Typical run: 5-12 s for 3-5 formulas on a quiet machine (about 2.3 GB); a 600 Da glycoside takes ~50 s and up to 4 GB.

layer_score_norm = score_norm(Total score): linear, 0 at 5.0 and 1 at 9.0, clipped; GOOD_ENOUGH 0.6 is raw 7.4.
  Calibration: 70 random high-resolution MassBank 2026.03 records (seed 20261028; one per compound, C H N O P S F Cl Br I, 100-800 Da, at
  least 6 peaks; 30 [M+H]+, 22 [M-H]-, 5 [M+Na]+, 5 [M+NH4]+, 5 [M+HCOO]-, 3 [M+Cl]-; no CASMI data) run with the true formula.
  69 gave structures; the true compound was among them for 62 and top-1 for 48 (70 %). The best candidate's score:
      raw score      norm          spectra   top-1 right   compound not in the database
      < 5.5          < 0.12           5        3 (60 %)       1
      5.5 - 6.0      0.12 - 0.25      7        4 (57 %)       2
      6.0 - 6.5      0.25 - 0.38      7        5 (71 %)       1
      6.5 - 7.0      0.38 - 0.50     17       13 (76 %)       2
      7.0 - 7.5      0.50 - 0.62     18       11 (61 %)       1
      7.5 - 8.0      0.62 - 0.75      9        7 (78 %)       0
      >= 8.0         >= 0.75          6        5 (83 %)       0
  Correct top-1 scores: min 4.6, p25 6.5, median 7.0, p75 7.5, p90 8.1, max 9.1; the 7 spectra whose compound was not in the database reached
  at most 7.2. From raw 7.4 (norm 0.6) up, 14 of 18 top-1 were right (78 %; 26 % of the spectra) and none of those 7 was among them.
  Read it as weak evidence: the Total score separates right from wrong only at the top (wrong top-1 are mostly isomers that explain the same
  fragments). The gap to the second candidate says more (top-1 right by gap: < 0.1 9/23, 0.1-0.3 11/15, >= 0.3 15/17, a single structure
  13/14), so summary["margin"] reports it; step 7 (isomers) checks what the spectrum can separate.
"""
from __future__ import annotations

import math
import re
import time

from mzdesk import adducts
from mzdesk.annotation.schema import LAYER_MSFINDER, STATUS_DONE, STATUS_FAILED, STATUS_SKIPPED, Candidate, LayerResult
from mzdesk.ms2 import msfinder_engine as mf

# ---------------------------------------------------------------- constants
GOOD_ENOUGH = 0.6                   # owner's design: the best normalised score at or above this means FRIGID (step 5) is not needed
SCORE_LOW, SCORE_HIGH = 5.0, 9.0    # the Total scores where layer_score_norm is 0 and 1 (calibration table above)
MAX_FORMULAS = 6                    # formulas searched in one run (the shared list holds up to 5 per target adduct, interleaved by rank)
MAX_PEAKS = 100                     # most intense MS2 peaks sent to MS-FINDER
FORMULA_WINDOW_DA = 0.0002          # Ms1Tolerance of a per-formula query: MS-FINDER sees that formula and (almost) nothing else
DEFAULT_WINDOW_DA = 0.001           # Ms1Tolerance of the shipped parameter file; the floor for the run without a formula list
FORMULAS_PER_QUERY = 2              # TryTopNMolecularFormulaSearch of a per-formula query
OWN_FORMULAS = 5                    # the same for the run without a formula list (the shipped value)
SEARCH_MINUTES = 0.5                # StructurePredictionTimeOut: no new structure is started after this long
ENGINE_TIMEOUT_S = 90               # the whole run
RETRY_TIMEOUT_S = 120               # owner's design 2026-10-30: the rerun of the unfinished formulas with TreeDepth=1
RETRY_TREE_DEPTH = 1                # owner's design 2026-10-30: the shipped parameter file has 2
RETRY_NOTE = "scored with TreeDepth=1 after the TreeDepth=2 search timed out"
RSS_CAP_MB = 6144                   # memory cap of the run (the engine's 4096 stops 600 Da glycosides)
DEFAULT_TOL_DA = 0.01               # MS-FINDER's own fragment tolerance
HALOGENS = {"F", "Cl", "Br", "I"}
NAME_MAX = 120


def score_norm(score):
    """MS-FINDER Total score -> 0-1: linear from 0 at SCORE_LOW to 1 at SCORE_HIGH, clipped (absolute: no other candidate is involved)."""
    if score is None or not math.isfinite(score):
        return 0.0
    return min(1.0, max(0.0, (score - SCORE_LOW) / (SCORE_HIGH - SCORE_LOW)))


# ---------------------------------------------------------------- helpers
def clean_peaks(peaks, limit=MAX_PEAKS):
    """[(m/z, intensity)] with finite positive values only, as floats, the `limit` most intense ones, sorted by m/z."""
    keep = [(float(mz), float(i)) for mz, i in peaks or [] if math.isfinite(mz) and math.isfinite(i) and mz > 0 and i > 0]
    return sorted(sorted(keep, key=lambda p: -p[1])[:limit])


def _counts(formula):
    out = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula or ""):
        out[el] = out.get(el, 0) + int(n or 1)
    return out


def _same(a, b):
    """Two formula strings with the same element counts (MS-FINDER's element order is not BUDDY's)."""
    return bool(a and b) and _counts(a) == _counts(b)


def _ion_mz(formula, adduct):
    """Theoretical m/z of neutral `formula` as `adduct` (msbuddy_engine arithmetic, the same as the BUDDY step's ppm); None if it cannot be computed."""
    try:
        from mzdesk.ms1 import msbuddy_engine as bd

        return bd.ion_mz_and_formula(formula, adduct)[0]
    except ImportError:
        return None


def _formula_items(formulas, default_adduct):
    """[(formula, adduct)] from the shared formula list: strings, dicts with "formula" (and "adduct"), or objects with those attributes.
    A formula comes once per adduct, in list order; items without a readable formula are dropped."""
    out = []
    for f in formulas or []:
        formula, adduct = (f, None) if isinstance(f, str) else (f.get("formula"), f.get("adduct")) if isinstance(f, dict) \
            else (getattr(f, "formula", None), getattr(f, "adduct", None))
        formula = (formula or "").strip()
        adduct = adducts.canonical(adduct) or default_adduct
        if formula and re.fullmatch(r"(?:[A-Z][a-z]?\d*)+", formula) and not any(_same(formula, g) and adduct == a for g, a in out):
            out.append((formula, adduct))
    return out


def _first_name(name):
    """MS-FINDER joins synonyms with ';': keep the first."""
    return (name or "").split(";")[0].strip()[:NAME_MAX]


def _plural(n, word):
    return f"{n} {word}" + ("" if n == 1 else "s")


# ---------------------------------------------------------------- the step
def run_msfinder_layer(precursor_mz, adduct, ms2_peaks, formulas, *, tol_ppm, tol_da, top_n=10, ms1_peaks=None):
    """Step 4. Returns (LayerResult, summary); never raises (a missing or failing engine gives status "failed" and a plain note).

    precursor_mz  measured precursor m/z
    adduct        any spelling mzdesk.adducts knows; sent to MS-FINDER as adducts.engine_name(adduct, "msfinder")
    ms2_peaks     [(m/z, intensity)]
    formulas      the shared formula list of step 2, best first: strings, dicts with "formula" (and "adduct"), or objects with those
                  attributes. A formula whose own adduct differs from `adduct` (the second adduct of step 1) is searched as that adduct.
                  Empty: MS-FINDER chooses its own formulas (one run at the measured precursor, within tol_ppm).
    tol_ppm       precursor tolerance: a formula whose theoretical m/z is further from precursor_mz is not searched
    tol_da        fragment tolerance (MS-FINDER's Ms2Tolerance)
    top_n         candidates kept (one per compound, best score first)
    ms1_peaks     isotope peaks; used only in the run without a formula list (with a list, step 2 has used them, see the module docstring)

    status: "skipped" when MS-FINDER does not apply (an adduct it cannot use, no MS2 peaks, no formula within tol_ppm), "failed" when the engine
    is missing, busy or fails, else "done" (also with no candidates). summary:
      {"good_enough": best_norm >= GOOD_ENOUGH, "best_norm", "best_score", "margin": best score minus the second (None: one candidate),
       "n_candidates", "mode": "formula_list" | "own_formulas",
       "formulas": [{"formula", "adduct", "mz_theoretical", "ppm", "accepted": MS-FINDER's formula rules took it (None: run stopped),
                     "n_structures", "best_score", "search_seconds", "status": searched | no structures | time limit | not accepted | not finished}],
       "skipped": [{"formula", "reason"}], "warnings": [plain sentences for the result], "seconds", "engine_seconds", "peak_rss_mb",
       "retried": formulas rerun with TreeDepth=1 after a time-limit stop (None: no rerun; [] with "own_formulas": the whole run)}
    Each entry of "formulas" also has "tree_depth" (1 for a rerun formula, else None: the shipped value) and "note" (RETRY_NOTE for a formula
    that finished in the rerun).
    """
    started = time.perf_counter()
    summary = {"good_enough": False, "best_norm": None, "best_score": None, "margin": None, "n_candidates": 0, "mode": None, "formulas": [],
               "skipped": [], "warnings": [], "seconds": 0.0, "engine_seconds": None, "peak_rss_mb": None, "retried": None}
    try:
        status, note, candidates = _search(summary, precursor_mz, adduct, ms2_peaks, formulas, tol_ppm, tol_da, top_n, ms1_peaks)
    except mf.MsFinderError as exc:
        status, note, candidates = STATUS_FAILED, str(exc), []
    except Exception as exc:                # noqa: BLE001 - a layer must never take the workflow down
        status, note, candidates = STATUS_FAILED, f"MS-FINDER step failed: {type(exc).__name__}: {exc}", []
    if candidates:
        summary["best_score"] = candidates[0].layer_score
        summary["best_norm"] = max(c.layer_score_norm for c in candidates)
        summary["good_enough"] = summary["best_norm"] >= GOOD_ENOUGH
        summary["margin"] = candidates[0].layer_score - candidates[1].layer_score if len(candidates) > 1 else None
    summary["n_candidates"] = len(candidates)
    summary["seconds"] = time.perf_counter() - started
    return LayerResult(LAYER_MSFINDER, status, note, summary["seconds"], candidates), summary


def _search(summary, precursor_mz, adduct, ms2_peaks, formulas, tol_ppm, tol_da, top_n, ms1_peaks):
    """-> (status, note, candidates) and fills `summary`. Raises MsFinderError when the engine is missing or fails."""
    canonical, msf_adduct = adducts.canonical(adduct), adducts.engine_name(adduct, "msfinder")
    if not msf_adduct:
        return STATUS_SKIPPED, f"MS-FINDER cannot use the adduct {adduct or '(none given)'}.", []
    peaks = clean_peaks(ms2_peaks)
    if not peaks:
        return STATUS_SKIPPED, "There are no MS2 peaks for MS-FINDER to fragment against.", []
    items = _formula_items(formulas, canonical)
    summary["mode"] = "formula_list" if items else "own_formulas"

    plan = []                                   # one search per formula of the list
    for formula, a in items:
        theo, name = _ion_mz(formula, a), adducts.engine_name(a, "msfinder")
        ppm = None if not theo else (precursor_mz - theo) / theo * 1e6
        if not name or not theo:
            summary["skipped"].append({"formula": formula, "reason": f"{formula} cannot be searched as {a}"})
        elif tol_ppm and abs(ppm) > tol_ppm:
            summary["skipped"].append({"formula": formula, "reason": f"{formula} as {a} is {ppm:+.1f} ppm from the precursor (limit {tol_ppm:g} ppm)"})
        elif len(plan) >= MAX_FORMULAS:
            summary["skipped"].append({"formula": formula, "reason": f"only the first {MAX_FORMULAS} formulas are searched"})
        else:
            plan.append({"formula": formula, "adduct": a, "name": name, "theo": theo, "ppm": ppm})
    if items and not plan:
        return STATUS_SKIPPED, "No formula of the list lies within the precursor tolerance, so MS-FINDER was not run.", []

    values = {"Ms2Tolerance": tol_da if tol_da and tol_da > 0 else DEFAULT_TOL_DA, "StructurePredictionTimeOut": SEARCH_MINUTES}
    if plan:                                    # MS-FINDER cannot take a formula list: each query's precursor is that formula's own m/z
        queries = [{"key": f"f{n}", "precursor_mz": round(p["theo"], 5), "positive": adducts.polarity(p["adduct"]) == "Positive", "adduct": p["name"],
                    "ms2": peaks} for n, p in enumerate(plan)]
        values.update(Ms1Tolerance=FORMULA_WINDOW_DA, TryTopNMolecularFormulaSearch=FORMULAS_PER_QUERY)
        elements = set(re.findall(r"[A-Z][a-z]?", "".join(p["formula"] for p in plan)))
    else:                                       # no list: MS-FINDER's own formula search at the measured precursor
        queries = [{"key": "own", "precursor_mz": precursor_mz, "positive": adducts.polarity(canonical) == "Positive", "adduct": msf_adduct,
                    "ms2": peaks, "ms1": clean_peaks(ms1_peaks, mf.MAX_PEAKS) or None}]
        values.update(Ms1Tolerance=max(DEFAULT_WINDOW_DA, (tol_ppm or 0) * 1e-6 * precursor_mz), TryTopNMolecularFormulaSearch=OWN_FORMULAS)
        elements = set()
    halogen = bool(elements & HALOGENS)
    res = _engine(queries, halogen, values, ENGINE_TIMEOUT_S)
    stopped, runs = res["stopped"], [res]
    groups = _groups(res, queries, plan, canonical, precursor_mz, tol_ppm)

    if _timed_out(stopped):                     # a time-limit stop: rerun what did not finish, ONCE, with the coarser fragment tree
        todo = [n for n, g in enumerate(groups) if g["info"]["status"] == "not finished"] if plan else ([] if groups else [0])
        if todo:
            rq = [queries[n] for n in todo]
            retry = _engine(rq, halogen, {**values, "TreeDepth": RETRY_TREE_DEPTH}, RETRY_TIMEOUT_S, any_error=True)
            if res["failed"] and retry["failed"]:           # nothing finished in either run: the step fails with both reasons
                raise mf.MsFinderError(f"{stopped} The rerun with TreeDepth={RETRY_TREE_DEPTH} failed too: {retry['stopped']}")
            runs.append(retry)
            again = _groups(retry, rq, [plan[n] for n in todo] if plan else [], canonical, precursor_mz, tol_ppm, RETRY_TREE_DEPTH)
            if plan:
                for n, g in zip(todo, again):
                    groups[n] = g
            else:
                groups = again
            summary["retried"] = [plan[n]["formula"] for n in todo] if plan else []
            stopped = retry["stopped"]
    summary["engine_seconds"] = sum(r["seconds"] for r in runs if r["seconds"] is not None) or None
    summary["peak_rss_mb"] = max((r["peak_rss_mb"] for r in runs if r["peak_rss_mb"] is not None), default=None)
    summary["formulas"] = [g["info"] for g in groups]
    if res["stopped"]:
        summary["warnings"].append(_stopped_warning(groups, plan, summary["retried"], res["stopped"], stopped))
    summary["warnings"] += [g["warning"] for g in groups if g["warning"]]
    summary["warnings"] += _skipped_warnings(summary["skipped"])

    candidates, seen = [], set()
    for c in sorted((c for g in groups for c in g["candidates"]), key=lambda c: -c.layer_score):         # stable: formula order breaks ties
        key = c.ik14 or c.smiles or id(c)
        if key not in seen:
            seen.add(key)
            candidates.append(c)
    candidates = candidates[:max(1, int(top_n))]
    for n, c in enumerate(candidates, start=1):
        c.id = f"msfinder_{n}"
    return STATUS_DONE, _note(candidates, groups, summary, stopped), candidates


def _timed_out(reason):
    """True for the engine's time-limit stop (run_watched: "... did not finish within N s and was stopped."); a memory stop, busy or missing is not."""
    return bool(reason) and "did not finish within" in reason


def _engine(queries, halogen, values, timeout, any_error=False):
    """mf.predict_queries; a run stopped by the time limit with nothing finished (the engine raises then) comes back as an empty answer
    with "failed": True. Any other MsFinderError (busy, missing, memory limit) is raised; with any_error (the rerun) it is returned the same way."""
    try:
        res = mf.predict_queries(queries, halogen=halogen, values=values, timeout=timeout, rss_cap_mb=RSS_CAP_MB)
        return {**res, "failed": False}
    except mf.MsFinderError as exc:
        if not (any_error or _timed_out(str(exc))):
            raise
        empty = {"formulas": [], "structures": [], "searches": {}}
        return {"results": {q["key"]: empty for q in queries}, "seconds": None, "peak_rss_mb": None, "log": "", "stopped": str(exc), "failed": True}


def _groups(res, queries, plan, canonical, precursor_mz, tol_ppm, depth=None):
    """The per-formula groups of one engine answer (depth: TreeDepth of a rerun, else None)."""
    if plan:
        return [_group(res["results"][q["key"]], p, res["stopped"], depth) for q, p in zip(queries, plan)]
    groups = []                                 # own formulas: one group per formula that has structures
    got = res["results"]["own"]
    for formula in dict.fromkeys(s["formula"] for s in got["structures"]):
        theo = _ion_mz(formula, canonical)
        ppm = None if not theo else (precursor_mz - theo) / theo * 1e6
        if theo and not (tol_ppm and abs(ppm) > tol_ppm):
            groups.append(_group(got, {"formula": formula, "adduct": canonical, "theo": theo, "ppm": ppm}, res["stopped"], depth))
    return groups


def _stopped_warning(groups, plan, retried, first, last):
    """The warning of a run that hit a limit: which formulas finished (and at which TreeDepth) and which did not, after the rerun if there was one.
    A formula finished when its search wrote a result file (status is not "not finished"); the unfinished ones have lost their candidates.
    retried: the formulas rerun (None: no rerun); first, last: the stop reasons of the first run and of the last one (None: it finished)."""
    done = [g["info"] for g in groups if g["info"]["status"] != "not finished"]
    lost = [g["info"]["formula"] for g in groups if g["info"]["status"] == "not finished"]
    names = lambda infos: ", ".join(i["formula"] for i in infos) or "none"          # noqa: E731
    text = "MS-FINDER was stopped before it finished (time or memory limit)"
    if retried is None:
        text += "; only the searches that had finished are used."
        text += f" Finished: {names(done)}."
        if lost:
            text += f" Not finished: {', '.join(lost)} (MS-FINDER writes a formula's structures only when its search ends, so they are lost)."
        elif not plan:
            text += " Any formula still being searched is lost, as MS-FINDER writes a formula's structures only when its search ends."
        return text
    shallow = [i for i in done if i["tree_depth"] == RETRY_TREE_DEPTH]
    text += f"; the unfinished {'formulas' if plan else 'run'} {'were' if plan else 'was'} repeated once with TreeDepth={RETRY_TREE_DEPTH} (the shipped file has 2)."
    text += f" Finished at the shipped TreeDepth: {names([i for i in done if i['tree_depth'] != RETRY_TREE_DEPTH])}."
    text += f" Finished at TreeDepth={RETRY_TREE_DEPTH}: {names(shallow)}" + (" (their scores come from a coarser fragment search and are marked)." if shallow else ".")
    if lost or (not plan and not done):
        text += f" Still not finished: {', '.join(lost) or 'the whole run'}" + (f" ({last.rstrip('.')})." if last else ".")
    return text


def _skipped_warnings(skipped):
    """One warning per reason: formulas skipped for the same reason (e.g. the formula limit) are listed together."""
    by_reason = {}
    for s in skipped:
        by_reason.setdefault(s["reason"], []).append(s["formula"])
    out = []
    for reason, formulas in by_reason.items():
        formulas = list(dict.fromkeys(formulas))
        head = reason[0].upper() + reason[1:]
        out.append(f"{head}; it was not searched." if len(formulas) == 1 else f"{head}; {len(formulas)} were not searched: {', '.join(formulas)}.")
    return out


def _group(got, p, stopped, depth=None):
    """One formula's structure search: its candidates, a summary entry and, when something went wrong, a warning.
    got: the engine's result for the query; p: {"formula", "adduct", "theo", "ppm"}; stopped: the engine's "stopped" reason or None;
    depth: the TreeDepth of a rerun (the group is then marked with RETRY_NOTE), None for the shipped value."""
    formula = p["formula"]
    rows = [s for s in got["structures"] if _same(s["formula"], formula)]
    fscore = next((f["score"] for f in got["formulas"] if _same(f["formula"], formula)), None)
    accepted = None if stopped else fscore is not None
    seconds = next((v for k, v in got["searches"].items() if _same(k, formula)), None)
    how, warning = "searched", None
    finished = any(_same(k, formula) for k in got["searches"])      # a result file exists: MS-FINDER writes it when the search ends
    if not rows and stopped and not finished:
        how = "not finished"
    elif accepted is False:
        how, warning = "not accepted", f"MS-FINDER's element rules do not accept {formula}, so no structure was searched for it."
    elif not rows:
        how = "no structures"
    elif seconds is not None and seconds >= 0.9 * SEARCH_MINUTES * 60:
        how, warning = "time limit", (f"The structure search for {formula} reached the {SEARCH_MINUTES * 60:g} s limit, "
                                      "so only part of its candidate structures were scored.")
    candidates = []
    for rank, s in enumerate(rows, start=1):
        inchikey = s["inchikey"] if re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", s["inchikey"] or "") else None
        candidates.append(Candidate(
            id="", layer=LAYER_MSFINDER, smiles=s["smiles"] or None, inchikey=inchikey, ik14=inchikey[:14] if inchikey else None,
            name=_first_name(s["name"]), formula=formula, adduct=p["adduct"], mz_error_ppm=p["ppm"], layer_score=s["score"],
            layer_score_norm=score_norm(s["score"]),
            evidence={LAYER_MSFINDER: {"total_score": s["score"], "subscores": dict(s["subscores"]), "databases": [db for db, _ in s["databases"]],
                                       "database_ids": {db: ident.split(";")[0] for db, ident in s["databases"]}, "ontology": s["ontology"] or None,
                                       "formula_score": fscore, "rank_in_formula": rank, "formula_structures": len(rows),
                                       **({"tree_depth": depth, "note": RETRY_NOTE} if depth else {})}},
            hypothesis=False))
    info = {"formula": formula, "adduct": p["adduct"], "mz_theoretical": p["theo"], "ppm": p["ppm"], "accepted": accepted, "n_structures": len(rows),
            "best_score": rows[0]["score"] if rows else None, "search_seconds": seconds, "status": how, "tree_depth": None, "note": None}
    if depth:                                   # a rerun: never silently mixed with the shipped TreeDepth
        info["tree_depth"] = depth
        info["note"] = RETRY_NOTE if how != "not finished" else f"not finished at TreeDepth={depth} either"
    return {"candidates": candidates, "info": info, "warning": warning}


def _note(candidates, groups, summary, stopped):
    """The one-sentence note of the step table."""
    caveats = []
    lost = [g["info"]["formula"] for g in groups if g["info"]["status"] == "not finished"]
    if summary["retried"] is not None:          # what finally happened: depth 2, depth 1 or nothing for each formula
        shallow = [g["info"]["formula"] for g in groups if g["info"]["tree_depth"] and g["info"]["status"] != "not finished"]
        if shallow:
            caveats.append(f"TreeDepth={RETRY_TREE_DEPTH} for {', '.join(shallow)} after the TreeDepth=2 search timed out")
        if lost:
            caveats.append(f"still not finished at TreeDepth={RETRY_TREE_DEPTH}: {', '.join(lost)}")
        elif stopped:
            caveats.append("MS-FINDER was stopped early")
    elif stopped:
        caveats.append("MS-FINDER was stopped early")
    limit = [g["info"]["formula"] for g in groups if g["info"]["status"] == "time limit"]
    if limit:
        caveats.append(f"time limit reached for {', '.join(limit)}")
    if summary["skipped"]:
        caveats.append(_plural(len(summary["skipped"]), "formula") + " not searched")
    tail = f" ({'; '.join(caveats)})" if caveats else ""
    if candidates:
        best = candidates[0]
        verdict = "good enough" if best.layer_score_norm >= GOOD_ENOUGH else f"not good enough (below {GOOD_ENOUGH:g})"
        return (f"{_plural(len(candidates), 'structure')} from {_plural(len(groups), 'formula')}; best {best.name or best.formula} ({best.formula}), "
                f"score {best.layer_score:.2f}, normalised {best.layer_score_norm:.2f}: {verdict}{tail}.")
    if summary["mode"] == "own_formulas":
        return f"MS-FINDER found no structure for this spectrum with its own formula search{tail}."
    names = ", ".join(g["info"]["formula"] for g in groups)
    if all(g["info"]["status"] == "no structures" for g in groups):
        return f"MS-FINDER's local structure databases hold no structure for {names}{tail}."
    return f"MS-FINDER gave no structure for {names}: " + "; ".join(f"{g['info']['formula']} {g['info']['status']}" for g in groups) + f"{tail}."
