"""Step 3 of the ESI MS/MS Annotation workflow: library search and the Level-2 gate.

run_library() searches the selected local libraries (MassBank, GNPS; mzdesk.ms2.library_search) with the query spectrum and judges every
hit that scores at least LIST_THRESHOLD against the gate. A hit passes only if ALL six checks hold:

  similarity  >= the similarity threshold (default GATE_SIMILARITY) on AT LEAST ONE of the chosen algorithms (default: weighted spectral
              entropy, the Similarity page's primary score; METRIC_CHOICES lists the ones offered). An algorithm that cannot be computed is
              left out (all of them: CosineGreedy replaces them). The best of the chosen scores is the hit's similarity (evidence["library"]
              ["metric"] names that algorithm, ["scores"] holds every chosen one); the first chosen algorithm orders the hits. Hits join the pool from the threshold - LIST_MARGIN
              (default LIST_THRESHOLD), again when any chosen algorithm reaches it.
  peaks       >= GATE_MIN_MATCHED matched peaks (greedy 1:1 pairing at the fragment tolerance); all of them if the query has fewer
  precursor   |observed - theoretical m/z of the hit's formula + adduct| <= tol_ppm; for a hit without a structure the library's own
              recorded precursor m/z stands in
  adduct      the hit's canonical adduct is one of the retained adduct hypotheses of step 1 (evidence["library"]["adduct_source"] = "step1"),
              or step 1 missed it but the library match confirms it: a known adduct of the polarity whose own formula + adduct fits the
              precursor within tol_ppm while every other check passes (adduct_source = "library"). The candidate keeps its own adduct.
  formula     the hit's formula is in the shared formula list of step 2 (the formulas BUDDY gave for the hit's adduct); when BUDDY gave
              none for that adduct (failed, not asked: only the best adduct, a close second and any retained adduct whose best BUDDY
              formula has prob >= 0.9 get formulas), the formula is outside BUDDY's search space (formula_step.buddy_covers: halogens,
              metals, huge molecules), or it is not in BUDDY's short list but its own mass fits the precursor, the check is "not
              evaluated" and the precursor check stands in
  isotope     cosine similarity (isotopes.py) of the hit's theoretical isotope pattern against the MS1 peaks >= GATE_ISOTOPE.
              "not evaluated" without a usable MS1 pattern; that never blocks.

A hit that passes is Level 2a; one that fails still joins the pool as a Level 3 candidate with the failed checks listed. A hit without
a structure (no usable SMILES or InChIKey) can pass but is capped at Level 3 (levels.py). The library
index has no formula column, so a hit's formula is computed from its SMILES (largest fragment, neutralised); a hit without a SMILES is
kept with smiles=None and no formula, and then cannot pass the formula check when BUDDY gave formulas for its adduct.

Candidates are one per compound (its best spectrum), gate-passing ones first, then by similarity, matched peaks, a known structure.
layer_score = layer_score_norm = the hit's rank score: its score on the first chosen algorithm (already 0-1 and absolute); the hits are ordered by it
too. The best of the chosen scores (evidence["library"]["similarity"]) only decides the gate.
Peak lists are [(mz, intensity), ...] (parse_spectrum output).
"""
from __future__ import annotations

import re
import time
from dataclasses import asdict, dataclass, field
from functools import lru_cache

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.rdBase import BlockLogs

from mzdesk.adducts import canonical, polarity
from mzdesk.annotation.adduct_step import clean_peaks
from mzdesk.annotation.formula_step import buddy_covers, ms1_usable
from mzdesk.annotation.merge import has_structure, structure_keys
from mzdesk.annotation.schema import STATUS_DONE, STATUS_FAILED, STATUS_SKIPPED, Candidate, LayerResult
from mzdesk.ms2 import library_index as lx
from mzdesk.ms2 import library_search as ls
from mzdesk.ms2 import metric_names as mn

# ---------------------------------------------------------------- constants (owner's design unless a comment says otherwise)
PRIMARY_METRIC = mn.ENTROPY_WEIGHTED        # the Similarity page's primary score (api/runners/ms2.py PRIMARY_METHOD)
FALLBACK_METRIC = mn.COSINE_GREEDY          # used when entropy cannot be computed (ms_entropy missing)
LIST_MARGIN = 0.10                          # hits scoring at least (gate similarity - this) join the pool, gate passed or not
GATE_SIMILARITY = 0.80
LIST_THRESHOLD = round(GATE_SIMILARITY - LIST_MARGIN, 4)      # the default list threshold, 0.70
# the algorithms the form offers (spec.py _PIPELINE_ADVANCED), in the Library search page's order; the deep-learning scores are left out (slow, model downloads)
METRIC_CHOICES = (mn.ENTROPY_WEIGHTED, mn.ENTROPY_UNWEIGHTED, mn.COSINE_GREEDY, mn.COSINE_HUNGARIAN, mn.MOD_COSINE_GREEDY, mn.MOD_COSINE_HUNGARIAN,
                  mn.NEUTRAL_LOSSES, mn.NIST_DOT, mn.REVERSE_DOT, mn.BINNED_COSINE)
assert all(m in mn.LIBRARY_ORDER for m in METRIC_CHOICES) and list(METRIC_CHOICES) == [m for m in mn.LIBRARY_ORDER if m in METRIC_CHOICES]


def as_metrics(metric):
    """One algorithm name or several -> a tuple in the choices' order, without repeats. Raises ValueError for an empty or unknown choice."""
    names = (metric,) if isinstance(metric, str) else tuple(metric or ())
    if not names:
        raise ValueError("Choose at least one library similarity algorithm.")
    bad = [m for m in names if m not in METRIC_CHOICES]
    if bad:
        raise ValueError(f"The library similarity algorithm must be one of: {', '.join(METRIC_CHOICES)} (got {bad[0]!r}).")
    return tuple(m for m in METRIC_CHOICES if m in names)
GATE_MIN_MATCHED = 4
GATE_ISOTOPE = 0.80
MAX_CANDIDATES = 20                         # compounds returned; a crowded precursor window can hold dozens of look-alike spectra
SEARCH_TOP_N = 200                          # spectra taken from the search before they are reduced to one per compound
ISOTOPE_RESOLUTION_PPM = 5.0                # tol_ppm up to this: Orbitrap / FT data (R 60k), above it Q-TOF (R 30k) for the isotope fine structure

PASS, FAIL, NOT_EVALUATED = "pass", "fail", "not evaluated"
CHECKS = (("similarity", "Similarity"), ("peaks", "Matched peaks"), ("precursor", "Precursor m/z"), ("adduct", "Adduct"),
          ("formula", "Formula"), ("isotope", "Isotope pattern"))


@dataclass
class GateSummary:
    """The Level-2 gate in one object.

    passed       at least one hit passed all checks (the orchestrator then skips steps 4-5)
    best_id      id of the best passing candidate (the first candidate), None if none passed
    top_id       id of the first candidate: the best passing one, else the highest similarity; `checks` describe this hit
    thresholds   the numbers used: {"metric" (the first chosen algorithm), "metrics" (all chosen), "similarity", "min_matched_peaks" (after the "all of them" rule), "precursor_ppm",
                 "isotope", "list"}
    checks       [{"id", "label", "status": pass | fail | not evaluated, "detail"}] of the top hit, in CHECKS order
    warnings     warnings of the library search (a library that is not installed, capped candidates, ...)
    """
    passed: bool = False
    best_id: str | None = None
    top_id: str | None = None
    thresholds: dict = field(default_factory=dict)
    checks: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    reference: dict | None = None              # best judged hit, for the spectrum display only

    def to_dict(self):
        return asdict(self)


# ---------------------------------------------------------------- helpers
@lru_cache(maxsize=4096)
def hit_formula(smiles):
    """Molecular formula of the neutral molecule M behind a library SMILES, or None.

    Largest fragment, neutralised (a stored [M+H]+ or [M-H]- structure becomes M), charge sign dropped; a permanent cation
    (quaternary ammonium) keeps its formula, which is right for its [M]+ adduct. Isotope labels stay in square brackets
    (C2[13C]H6O2): no m/z is computed from such a formula and BUDDY cannot list it.
    """
    if not smiles:
        return None
    try:
        with BlockLogs():
            mol = Chem.MolFromSmiles(smiles)
            if mol is None or not mol.GetNumHeavyAtoms():
                return None
            mol = max(Chem.GetMolFrags(mol, asMols=True), key=lambda p: p.GetNumHeavyAtoms())
            return re.sub(r"[+-]\d*$", "", rdMolDescriptors.CalcMolFormula(rdMolStandardize.Uncharger().uncharge(mol), True, False))
    except Exception:           # an odd structure must not stop the search
        return None


_INCHIKEY = re.compile(r"[A-Z]{14}-[A-Z]{10}-[A-Z]")


@lru_cache(maxsize=4096)
def _key_of_smiles(smiles):
    return structure_keys(smiles)[0]


def _inchikey(hit):
    """The hit's InChIKey if it is one (the libraries also store "N/A" and the like), else the one derived from its SMILES
    (GNPS rows carry a SMILES but no key), else None."""
    given = hit.get("inchikey") or ""
    if _INCHIKEY.fullmatch(given):
        return given
    return _key_of_smiles(hit["smiles"]) if hit.get("smiles") else None


def _has_structure(hit):
    """True when the hit carries a usable SMILES or an InChIKey; without either it cannot be verified (Level 3 at most, see levels.py)."""
    return has_structure(hit.get("smiles")) or bool(_inchikey(hit))


def _element_counts(formula):
    counts = {}
    for element, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula or ""):
        counts[element] = counts.get(element, 0) + int(n or 1)
    return counts


def same_formula(a, b):
    """True when two neutral formula strings have the same element counts (the order of the elements does not matter)."""
    return bool(a and b) and _element_counts(a) == _element_counts(b)


def _formula_lists(formulas):
    """{adduct: [formula, ...]} from FormulaCandidates or dicts with "formula" and "adduct"; plain strings (no adduct) go under None."""
    lists = {}
    for item in formulas or []:
        if isinstance(item, str):
            formula, adduct = item, None
        elif isinstance(item, dict):
            formula, adduct = item["formula"], item.get("adduct")
        else:
            formula, adduct = item.formula, item.adduct
        lists.setdefault(canonical(adduct), []).append(formula)
    return lists


def _ion_mz(formula, adduct):
    """Theoretical m/z of neutral `formula` as `adduct` (msbuddy_engine arithmetic, as for BUDDY's ppm); None if it cannot be computed."""
    if not formula or not adduct or "[" in formula:
        return None
    try:
        from mzdesk.ms1 import msbuddy_engine as bd

        return bd.ion_mz_and_formula(formula, adduct)[0]
    except ImportError:
        return None


def _isotope_similarity(formula, adduct, ms1, resolution, tol_ppm):
    """Cosine similarity of the isotope pattern of formula + adduct against the MS1 peaks; None if the pattern cannot be computed."""
    from mzdesk.ms1 import isotopes

    try:
        pattern = isotopes.pattern(formula, adduct, resolution)
        return isotopes.similarity(pattern["peaks"], ms1, resolution, tol_ppm)["score"]
    except (ValueError, ImportError):
        return None


def list_threshold(gate_similarity):
    """Hits from this similarity on join the pool: the gate threshold minus LIST_MARGIN, not below 0."""
    return max(0.0, round(gate_similarity - LIST_MARGIN, 4))


def check_settings(metric, gate_similarity):
    """Raise ValueError (a clear message) for an algorithm outside METRIC_CHOICES (none, or an unknown one) or a threshold outside 0-1."""
    as_metrics(metric)
    if isinstance(gate_similarity, bool) or not isinstance(gate_similarity, (int, float)) or not 0 <= gate_similarity <= 1:
        raise ValueError(f"The similarity threshold must be a number from 0 to 1 (got {gate_similarity!r}).")


# ---------------------------------------------------------------- the gate
class _Gate:
    """The gate of one query. judge(hit) runs the six checks on one library hit and returns them with what was derived on the way."""

    def __init__(self, *, precursor_mz, tol_ppm, n_query, adducts, formulas, ms1, positive=True, halogen=False, metric=PRIMARY_METRIC,
                 gate_similarity=GATE_SIMILARITY, min_matched_peaks=GATE_MIN_MATCHED):
        self.precursor_mz, self.tol_ppm, self.halogen = precursor_mz, tol_ppm, halogen
        self.metrics = self.requested_metrics = as_metrics(metric)      # the chosen algorithms; self.metrics becomes what the search really used
        self.gate_similarity = gate_similarity
        self.list_threshold = list_threshold(gate_similarity)
        self.polarity = "Positive" if positive else "Negative"
        self.min_matched = min(min_matched_peaks, n_query)                      # "all of them" when the query has fewer peaks
        self.adducts = [a for a in (canonical(getattr(x, "adduct", x)) for x in adducts or []) if a]
        self.formulas = _formula_lists(formulas)
        self.ms1 = ms1                                                          # usable MS1 peaks, or None
        self.resolution = "Orbitrap 60k" if tol_ppm <= ISOTOPE_RESOLUTION_PPM else "Q-TOF 30k"

    @property
    def metric(self):
        """The first chosen algorithm (the one the list is ordered by)."""
        return self.metrics[0]

    @property
    def requested_metric(self):
        return self.requested_metrics[0]

    @property
    def metric_text(self):
        """"A", or "A or B or C" when several algorithms count."""
        return " or ".join(self.metrics)

    def thresholds(self):
        return {"metric": self.metric, "metrics": list(self.metrics), "similarity": self.gate_similarity, "min_matched_peaks": self.min_matched,
                "precursor_ppm": self.tol_ppm, "isotope": GATE_ISOTOPE, "list": self.list_threshold}

    def judge(self, hit):
        scores = {m: hit["scores"].get(m) for m in self.metrics}
        got = {m: v for m, v in scores.items() if v is not None}
        best = max(got, key=lambda m: (got[m], -self.metrics.index(m))) if got else None       # the best score; a tie goes to the earlier algorithm
        similarity = got[best] if best else None
        rank_score = next((scores[m] for m in self.metrics if scores[m] is not None), None)       # the first chosen algorithm that scored it: orders the hits
        formula = hit_formula(hit.get("smiles"))
        adduct = canonical(hit.get("adduct"))
        theoretical = _ion_mz(formula, adduct)
        if theoretical:
            error, basis = (self.precursor_mz - theoretical) / theoretical * 1e6, f"{formula} {adduct}"
        else:
            error, basis = hit["ppm"], "the library's recorded precursor m/z"
        iso = _isotope_similarity(formula, adduct, self.ms1, self.resolution, self.tol_ppm) if self.ms1 and formula and adduct else None

        res = {}
        ok = similarity is not None and similarity >= self.gate_similarity
        if len(self.metrics) == 1:
            res["similarity"] = (PASS if ok else FAIL, f"{self.metric} {'n/a' if similarity is None else format(similarity, '.2f')} "
                                                       f"{'>=' if ok else '<'} {self.gate_similarity:.2f}")
        else:
            listed = ", ".join(f"{m} {'n/a' if v is None else format(v, '.2f')}" for m, v in scores.items())
            res["similarity"] = (PASS if ok else FAIL, f"{listed}; {'at least one' if ok else 'none'} {'>=' if ok else 'reaches'} {self.gate_similarity:.2f}")
        ok = hit["matches"] >= self.min_matched
        res["peaks"] = (PASS if ok else FAIL, f"{hit['matches']} matched, {self.min_matched} needed")
        ok = abs(error) <= self.tol_ppm
        res["precursor"] = (PASS if ok else FAIL, f"{error:+.1f} ppm against {basis} (limit {self.tol_ppm:g} ppm)")
        listed = list(dict.fromkeys(self.formulas.get(adduct, []) + self.formulas.get(None, [])))       # BUDDY's formulas for this adduct
        if not listed:
            res["formula"] = (NOT_EVALUATED, f"BUDDY gave no formula for {adduct or 'this adduct'}; the precursor check stands in for it")
        elif formula is None:
            res["formula"] = (FAIL, "the library entry has no usable structure, so its formula is unknown")
        elif not buddy_covers(formula, self.halogen):
            res["formula"] = (NOT_EVALUATED, f"{formula} is outside BUDDY's search space; the precursor check stands in for it")
        elif any(same_formula(formula, f) for f in listed):
            res["formula"] = (PASS, f"{formula} is among the BUDDY formulas for {adduct}")
        elif theoretical and abs(error) <= self.tol_ppm:
            # BUDDY lists at most 5 formulas, so a formula missing from that short list is weak evidence when its own mass fits the precursor
            res["formula"] = (NOT_EVALUATED, f"{formula} is not among BUDDY's top {len(listed)} formulas for {adduct} ({', '.join(listed)}), "
                                             f"but it matches the precursor at {error:+.1f} ppm; the precursor check stands in for it")
        else:
            res["formula"] = (FAIL, f"{formula} is not among the BUDDY formulas for {adduct} ({', '.join(listed)})")
        if self.ms1 is None:
            res["isotope"] = (NOT_EVALUATED, "no usable MS1 isotope pattern was given")
        elif iso is None:
            res["isotope"] = (NOT_EVALUATED, "the isotope pattern needs the hit's formula and adduct")
        else:
            ok = iso >= GATE_ISOTOPE
            res["isotope"] = (PASS if ok else FAIL, f"isotope similarity {iso:.2f} {'>=' if ok else '<'} {GATE_ISOTOPE:.2f}")

        # The adduct check goes last: step 1 may have missed the adduct, and a good library match then confirms the hit's own adduct.
        retained, source = ", ".join(self.adducts) or "none", "step1"
        if adduct in self.adducts:
            res["adduct"] = (PASS, f"{adduct} is one of the retained adducts ({retained})")
        elif adduct is None:
            raw = (hit.get("adduct") or "").strip()
            res["adduct"] = (FAIL, f"the library entry's adduct {raw!r} is not a known adduct" if raw else "the library entry records no adduct")
        elif polarity(adduct) == self.polarity and theoretical and all(status != FAIL for status, _ in res.values()):
            source = "library"
            res["adduct"] = (PASS, f"{adduct} is not one of the retained adducts ({retained}), but it is a known {self.polarity.lower()}-mode adduct "
                                   f"and the library match confirms it ({error:+.1f} ppm against {basis})")
        else:
            res["adduct"] = (FAIL, f"{adduct} is not one of the retained adducts ({retained})")

        checks = [{"id": key, "label": label, "status": res[key][0], "detail": res[key][1]} for key, label in CHECKS]
        return {"checks": checks, "passed": all(c["status"] != FAIL for c in checks), "formula": formula, "adduct": adduct,
                "similarity": similarity, "rank_score": rank_score, "metric": best, "scores": scores, "error": error, "error_basis": basis, "error_from_formula": bool(theoretical),
                "isotope": res["isotope"][0], "isotope_similarity": iso, "adduct_source": source}


# ---------------------------------------------------------------- the step
def _search(peaks, libraries, precursor_mz, positive, tol_ppm, tol_da, metric=PRIMARY_METRIC, threshold=LIST_THRESHOLD):
    """library_search.search with the algorithm(s) `metric` (a name or several; it gets the precursor m/z, which ModifiedCosine needs).
    An algorithm that cannot be computed is left out: the entropy ones when ms_entropy is missing, then everything when nothing could be
    scored, and CosineGreedy replaces the lot. -> (result, the algorithm(s) used: a name when one was given, else a tuple)"""
    one = isinstance(metric, str)
    names = as_metrics(metric)

    def run(chosen):
        return ls.search(peaks, "MS/MS", libraries, list(chosen), threshold, precursor_mz=precursor_mz, positive=positive, adduct=None,
                         precursor_tol_ppm=tol_ppm, tolerance=tol_da, top_n=SEARCH_TOP_N)

    def used(chosen):
        return chosen[0] if one else tuple(chosen)

    if names == (FALLBACK_METRIC,):
        return run(names), used(names)
    try:
        res = run(names)
        if not res["n_candidates"] or res["n_metrics_scored"]:
            return res, used(names)
    except ImportError:                          # ms_entropy is not installed: keep the other chosen algorithms
        rest = tuple(m for m in names if m not in mn.ENTROPY)
        if rest:
            res = run(rest)
            if not res["n_candidates"] or res["n_metrics_scored"]:
                return res, used(rest)
    return run((FALLBACK_METRIC,)), used((FALLBACK_METRIC,))


def _reference_spectrum(hit):
    """Read the selected hit's indexed spectrum for display; never affect the search result."""
    try:
        mz, intensity = lx.load(hit["library"]).peaks(hit["row"])
        peaks = [(float(m), float(i)) for m, i in zip(mz, intensity)]
        if not peaks:
            return None
        return {"peaks": peaks, "library": hit["library"], "accession": hit["accession"],
                "name": hit.get("name") or "", "precursor_mz": hit.get("precursor_mz"),
                "smiles": hit["smiles"] if has_structure(hit.get("smiles")) else None}
    except Exception:  # noqa: BLE001 - a missing display spectrum must not change the annotation
        return None


def run_library(ms2_peaks, precursor_mz, positive, libraries, *, adducts, formulas, ms1_peaks=None, tol_ppm=10.0, tol_da=0.02, halogen=False,
                metric=PRIMARY_METRIC, gate_similarity=GATE_SIMILARITY, min_matched_peaks=GATE_MIN_MATCHED):
    """Step 3 -> (LayerResult, GateSummary).

    ms2_peaks     [(mz, intensity)] query spectrum
    precursor_mz  observed precursor m/z
    positive      True = positive ion mode
    libraries     names of the libraries to search (spec.LIBRARIES: "MassBank", "GNPS"); empty -> status "skipped"
    adducts       the retained adduct hypotheses of step 1: names, or AdductHypothesis objects (AdductResult.retained_adducts). A hit with
                  another known adduct of the polarity still passes when the library match confirms that adduct (see the module docstring)
    formulas      the shared formula list of step 2: FormulaCandidates, dicts with "formula" and "adduct", or strings (valid for every adduct);
                  a hit whose adduct has no formula in it is judged by the precursor check only
    ms1_peaks     [(mz, intensity)] or None; the isotope check runs when it holds a pattern at the precursor (formula_step.ms1_usable)
    tol_ppm       precursor tolerance (ppm): the library window and the gate's precursor check
    tol_da        fragment tolerance (Da) for the spectrum comparison
    halogen       the halogen setting step 2 used (run_adduct_and_formula(halogen=...)): without it BUDDY cannot list a Cl / Br / F / I formula
    metric        the similarity algorithm of the gate and the search: one name or several (a hit passes when at least one reaches the threshold),
                  all from METRIC_CHOICES; one that cannot be computed is left out and CosineGreedy replaces all of them if none can
    min_matched_peaks  the fewest matched peaks a hit may have (default GATE_MIN_MATCHED; all of them when the query has fewer)
    gate_similarity  the similarity threshold, 0-1; hits from this minus LIST_MARGIN (not below 0) on join the pool
    Raises ValueError for no algorithm, one outside METRIC_CHOICES, or a threshold outside 0-1.

    LayerResult: layer "library"; candidates best first, each with evidence["library"] = {gate_passed, failed_checks, similarity, metric,
    matched_peaks, has_structure (False caps the level at 3), precursor_error_ppm, isotope ("pass" | "fail" | "not evaluated"), isotope_similarity, adduct_source ("step1" | "library"),
    checks, library, accession, ...}.
    status: skipped (no library selected or installed, no peaks), failed (the search crashed), done.
    """
    check_settings(metric, gate_similarity)
    started = time.perf_counter()
    peaks = clean_peaks(ms2_peaks)
    ms1_ok, _ = ms1_usable(ms1_peaks, precursor_mz, tol_ppm)
    gate = _Gate(precursor_mz=precursor_mz, tol_ppm=tol_ppm, n_query=len(peaks), adducts=adducts, formulas=formulas,
                 ms1=clean_peaks(ms1_peaks) if ms1_ok else None, positive=positive, halogen=halogen, metric=metric,
                 gate_similarity=gate_similarity, min_matched_peaks=min_matched_peaks)

    def finish(status, note, candidates=(), summary=None):
        layer = LayerResult("library", status, note, time.perf_counter() - started, list(candidates))
        return layer, summary or GateSummary(thresholds=gate.thresholds())

    libraries = list(dict.fromkeys(libraries or []))
    if not libraries:
        return finish(STATUS_SKIPPED, "No library was selected, so no library spectrum was searched.")
    if not peaks:
        return finish(STATUS_SKIPPED, "The MS2 spectrum has no peaks, so no library spectrum was searched.")
    try:
        res, used = _search(peaks, libraries, precursor_mz, positive, tol_ppm, tol_da, gate.metrics, gate.list_threshold)
        gate.metrics = as_metrics(used)                  # what the search really used (a stand-in in tests may answer with one name)
    except lx.LibraryMissing as exc:
        return finish(STATUS_SKIPPED, str(exc))
    except Exception as exc:                                                    # a broken library must not stop the workflow
        return finish(STATUS_FAILED, f"The library search failed ({type(exc).__name__}: {exc}).")

    best = {}                                    # compound -> [sort key, hit, verdict, number of spectra]: the best spectrum of each
    for hit in res["hits"]:
        verdict = gate.judge(hit)
        key = (structure_keys(hit.get("smiles"))[1] or (_inchikey(hit) or "")[:14] or (hit.get("name") or "").strip().lower()
               or (hit["library"], hit["accession"]))
        order = (not verdict["passed"], verdict["passed"] and not _has_structure(hit), -round(verdict.get("rank_score", verdict["similarity"]), 4), -hit["matches"], verdict["formula"] is None, abs(verdict["error"]))
        if key not in best:
            best[key] = [order, hit, verdict, 1]
            continue
        best[key][3] += 1
        if order < best[key][0]:
            best[key][:3] = order, hit, verdict
    chosen = sorted(best.values(), key=lambda b: b[0])[:MAX_CANDIDATES]
    candidates = [_candidate(f"lib{n}", hit, verdict, count, verdict.get("metric") or gate.metric) for n, (_, hit, verdict, count) in enumerate(chosen, start=1)]

    passing = [c for c in candidates if c.evidence["library"]["gate_passed"]]
    warnings = list(res["warnings"]) + [f"{name}: {text}" for name, text in res["notes"].items()]
    if gate.metrics != gate.requested_metrics:
        left_out = [m for m in gate.requested_metrics if m not in gate.metrics]
        warnings.append(f"{', '.join(left_out)} could not be computed, so {', '.join(gate.metrics)} {'was' if len(gate.metrics) == 1 else 'were'} "
                        f"used for the library search and the gate.")
    if res["n_kept"] > SEARCH_TOP_N:
        warnings.append(f"{res['n_kept']:,} library spectra scored {gate.list_threshold:.2f} or higher; only the best {SEARCH_TOP_N} were judged by the gate.")
    summary = GateSummary(passed=bool(passing), best_id=passing[0].id if passing else None, top_id=candidates[0].id if candidates else None,
                          thresholds=gate.thresholds(), checks=chosen[0][2]["checks"] if chosen else [], warnings=warnings,
                          reference=_reference_spectrum(chosen[0][1]) if chosen else None)
    return finish(STATUS_DONE, _note(res, candidates, passing, tol_ppm, precursor_mz, gate), candidates, summary)


def _candidate(cid, hit, verdict, n_spectra, metric):
    """The Candidate of one compound, from its best spectrum."""
    inchikey = _inchikey(hit)
    smiles = hit["smiles"] if has_structure(hit.get("smiles")) else None                  # the libraries also store "N/A" and the like
    failed = [f"{c['label']}: {c['detail']}" for c in verdict["checks"] if c["status"] == FAIL]
    evidence = {"gate_passed": verdict["passed"], "failed_checks": failed, "similarity": verdict["similarity"], "metric": metric, "scores": verdict.get("scores"),
                "matched_peaks": hit["matches"], "precursor_error_ppm": verdict["error"], "precursor_error_basis": verdict["error_basis"],
                "isotope": verdict["isotope"], "isotope_similarity": verdict["isotope_similarity"], "adduct_source": verdict["adduct_source"],
                "checks": verdict["checks"], "has_structure": _has_structure(hit),
                "library": hit["library"], "accession": hit["accession"], "instrument": hit["instrument"],
                "collision_energy": hit["collision_energy"], "library_adduct": hit.get("adduct") or None,
                "library_precursor_mz": hit["precursor_mz"], "n_spectra": n_spectra}
    return Candidate(id=cid, layer="library", smiles=smiles, inchikey=inchikey, ik14=inchikey[:14] if inchikey else None, name=hit.get("name") or "",
                     formula=verdict["formula"], adduct=verdict["adduct"],
                     mz_error_ppm=verdict["error"] if verdict["error_from_formula"] else None, layer_score=verdict.get("rank_score", verdict["similarity"]),
                     layer_score_norm=min(1.0, max(0.0, verdict.get("rank_score", verdict["similarity"]))), evidence={"library": evidence})


def _note(res, candidates, passing, tol_ppm, precursor_mz, gate):
    if not res["n_candidates"]:
        return f"No library spectrum within {tol_ppm:g} ppm of m/z {precursor_mz:.4f}."
    if not candidates:
        return f"{res['n_candidates']:,} library spectra compared; none scored {gate.list_threshold:.2f} or higher on {gate.metric_text}."
    text = (f"{res['n_kept']:,} library spectra scored {gate.list_threshold:.2f} or higher on {gate.metric_text} ({len(candidates)} compounds listed); "
            f"the gate needs {gate.gate_similarity:.2f}.")
    if passing:
        ev = passing[0].evidence["library"]
        return text + f" Gate passed by {passing[0].name or passing[0].id} (similarity {ev['similarity']:.2f}, {ev['matched_peaks']} matched peaks)."
    top = candidates[0]
    return text + f" Gate not passed: the best hit {top.name or top.id} failed " + "; ".join(top.evidence["library"]["failed_checks"]) + "."
