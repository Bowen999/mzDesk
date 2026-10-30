"""Step 3 of the ESI MS/MS Annotation workflow: library search and the Level-2 gate.

run_library() searches the selected local libraries (MassBank, GNPS; mzdesk.ms2.library_search) with the query spectrum and judges every
hit that scores at least LIST_THRESHOLD against the gate. A hit passes only if ALL six checks hold:

  similarity  >= GATE_SIMILARITY on the primary metric: weighted spectral entropy (the Similarity page's primary score), CosineGreedy
              if entropy cannot be computed. evidence["library"]["metric"] says which one was used.
  peaks       >= GATE_MIN_MATCHED matched peaks (greedy 1:1 pairing at the fragment tolerance); all of them if the query has fewer
  precursor   |observed - theoretical m/z of the hit's formula + adduct| <= tol_ppm; for a hit without a structure the library's own
              recorded precursor m/z stands in
  adduct      the hit's canonical adduct is one of the retained adduct hypotheses of step 1 (evidence["library"]["adduct_source"] = "step1"),
              or step 1 missed it but the library match confirms it: a known adduct of the polarity whose own formula + adduct fits the
              precursor within tol_ppm while every other check passes (adduct_source = "library"). The candidate keeps its own adduct.
  formula     the hit's formula is in the shared formula list of step 2 (the formulas BUDDY gave for the hit's adduct); when BUDDY gave
              none for that adduct (failed, not asked: only the best adduct and a close second get formulas), or the formula is outside
              BUDDY's search space (formula_step.buddy_covers: halogens, metals, huge molecules), the precursor check stands in
  isotope     cosine similarity (isotopes.py) of the hit's theoretical isotope pattern against the MS1 peaks >= GATE_ISOTOPE.
              "not evaluated" without a usable MS1 pattern; that never blocks.

A hit that passes is Level 2a; one that fails still joins the pool as a Level 3 candidate with the failed checks listed. The library
index has no formula column, so a hit's formula is computed from its SMILES (largest fragment, neutralised); a hit without a SMILES is
kept with smiles=None and no formula, and then cannot pass the formula check when BUDDY gave formulas for its adduct.

Candidates are one per compound (its best spectrum), gate-passing ones first, then by similarity, matched peaks, a known structure.
layer_score = layer_score_norm = the similarity of the primary metric (already 0-1 and absolute).
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
LIST_THRESHOLD = 0.70                       # hits scoring at least this join the pool, gate passed or not
GATE_SIMILARITY = 0.80
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
    thresholds   the numbers used: {"metric", "similarity", "min_matched_peaks" (after the "all of them" rule), "precursor_ppm",
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


def _inchikey(hit):
    """The hit's InChIKey if it is one (the libraries also store "N/A" and the like), else None."""
    return hit["inchikey"] if _INCHIKEY.fullmatch(hit.get("inchikey") or "") else None


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


# ---------------------------------------------------------------- the gate
class _Gate:
    """The gate of one query. judge(hit) runs the six checks on one library hit and returns them with what was derived on the way."""

    def __init__(self, *, precursor_mz, tol_ppm, n_query, adducts, formulas, ms1, positive=True, halogen=False):
        self.precursor_mz, self.tol_ppm, self.metric, self.halogen = precursor_mz, tol_ppm, PRIMARY_METRIC, halogen
        self.polarity = "Positive" if positive else "Negative"
        self.min_matched = min(GATE_MIN_MATCHED, n_query)                      # "all of them" when the query has fewer peaks
        self.adducts = [a for a in (canonical(getattr(x, "adduct", x)) for x in adducts or []) if a]
        self.formulas = _formula_lists(formulas)
        self.ms1 = ms1                                                          # usable MS1 peaks, or None
        self.resolution = "Orbitrap 60k" if tol_ppm <= ISOTOPE_RESOLUTION_PPM else "Q-TOF 30k"

    def thresholds(self):
        return {"metric": self.metric, "similarity": GATE_SIMILARITY, "min_matched_peaks": self.min_matched,
                "precursor_ppm": self.tol_ppm, "isotope": GATE_ISOTOPE, "list": LIST_THRESHOLD}

    def judge(self, hit):
        similarity = hit["scores"].get(self.metric)
        formula = hit_formula(hit.get("smiles"))
        adduct = canonical(hit.get("adduct"))
        theoretical = _ion_mz(formula, adduct)
        if theoretical:
            error, basis = (self.precursor_mz - theoretical) / theoretical * 1e6, f"{formula} {adduct}"
        else:
            error, basis = hit["ppm"], "the library's recorded precursor m/z"
        iso = _isotope_similarity(formula, adduct, self.ms1, self.resolution, self.tol_ppm) if self.ms1 and formula and adduct else None

        res = {}
        ok = similarity is not None and similarity >= GATE_SIMILARITY
        res["similarity"] = (PASS if ok else FAIL, f"{self.metric} {similarity:.2f} {'>=' if ok else '<'} {GATE_SIMILARITY:.2f}")
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
                "similarity": similarity, "error": error, "error_basis": basis, "error_from_formula": bool(theoretical),
                "isotope": res["isotope"][0], "isotope_similarity": iso, "adduct_source": source}


# ---------------------------------------------------------------- the step
def _search(peaks, libraries, precursor_mz, positive, tol_ppm, tol_da):
    """library_search.search with the primary metric; with the fallback metric when the primary one cannot be computed. -> (result, metric)"""
    def run(metric):
        return ls.search(peaks, "MS/MS", libraries, [metric], LIST_THRESHOLD, precursor_mz=precursor_mz, positive=positive, adduct=None,
                         precursor_tol_ppm=tol_ppm, tolerance=tol_da, top_n=SEARCH_TOP_N)

    try:
        res = run(PRIMARY_METRIC)
        if not res["n_candidates"] or res["n_metrics_scored"]:
            return res, PRIMARY_METRIC
    except ImportError:                          # ms_entropy is not installed
        pass
    return run(FALLBACK_METRIC), FALLBACK_METRIC


def _reference_spectrum(hit):
    """Read the selected hit's indexed spectrum for display; never affect the search result."""
    try:
        mz, intensity = lx.load(hit["library"]).peaks(hit["row"])
        peaks = [(float(m), float(i)) for m, i in zip(mz, intensity)]
        if not peaks:
            return None
        return {"peaks": peaks, "library": hit["library"], "accession": hit["accession"],
                "name": hit.get("name") or "", "precursor_mz": hit.get("precursor_mz")}
    except Exception:  # noqa: BLE001 - a missing display spectrum must not change the annotation
        return None


def run_library(ms2_peaks, precursor_mz, positive, libraries, *, adducts, formulas, ms1_peaks=None, tol_ppm=10.0, tol_da=0.02, halogen=False):
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

    LayerResult: layer "library"; candidates best first, each with evidence["library"] = {gate_passed, failed_checks, similarity, metric,
    matched_peaks, precursor_error_ppm, isotope ("pass" | "fail" | "not evaluated"), isotope_similarity, adduct_source ("step1" | "library"),
    checks, library, accession, ...}.
    status: skipped (no library selected or installed, no peaks), failed (the search crashed), done.
    """
    started = time.perf_counter()
    peaks = clean_peaks(ms2_peaks)
    ms1_ok, _ = ms1_usable(ms1_peaks, precursor_mz, tol_ppm)
    gate = _Gate(precursor_mz=precursor_mz, tol_ppm=tol_ppm, n_query=len(peaks), adducts=adducts, formulas=formulas,
                 ms1=clean_peaks(ms1_peaks) if ms1_ok else None, positive=positive, halogen=halogen)

    def finish(status, note, candidates=(), summary=None):
        layer = LayerResult("library", status, note, time.perf_counter() - started, list(candidates))
        return layer, summary or GateSummary(thresholds=gate.thresholds())

    libraries = list(dict.fromkeys(libraries or []))
    if not libraries:
        return finish(STATUS_SKIPPED, "No library was selected, so no library spectrum was searched.")
    if not peaks:
        return finish(STATUS_SKIPPED, "The MS2 spectrum has no peaks, so no library spectrum was searched.")
    try:
        res, gate.metric = _search(peaks, libraries, precursor_mz, positive, tol_ppm, tol_da)
    except lx.LibraryMissing as exc:
        return finish(STATUS_SKIPPED, str(exc))
    except Exception as exc:                                                    # a broken library must not stop the workflow
        return finish(STATUS_FAILED, f"The library search failed ({type(exc).__name__}: {exc}).")

    best = {}                                    # compound -> [sort key, hit, verdict, number of spectra]: the best spectrum of each
    for hit in res["hits"]:
        verdict = gate.judge(hit)
        key = (structure_keys(hit.get("smiles"))[1] or (_inchikey(hit) or "")[:14] or (hit.get("name") or "").strip().lower()
               or (hit["library"], hit["accession"]))
        order = (not verdict["passed"], -round(verdict["similarity"], 4), -hit["matches"], verdict["formula"] is None, abs(verdict["error"]))
        if key not in best:
            best[key] = [order, hit, verdict, 1]
            continue
        best[key][3] += 1
        if order < best[key][0]:
            best[key][:3] = order, hit, verdict
    chosen = sorted(best.values(), key=lambda b: b[0])[:MAX_CANDIDATES]
    candidates = [_candidate(f"lib{n}", hit, verdict, count, gate.metric) for n, (_, hit, verdict, count) in enumerate(chosen, start=1)]

    passing = [c for c in candidates if c.evidence["library"]["gate_passed"]]
    warnings = list(res["warnings"]) + [f"{name}: {text}" for name, text in res["notes"].items()]
    if res["n_kept"] > SEARCH_TOP_N:
        warnings.append(f"{res['n_kept']:,} library spectra scored {LIST_THRESHOLD:.2f} or higher; only the best {SEARCH_TOP_N} were judged by the gate.")
    summary = GateSummary(passed=bool(passing), best_id=passing[0].id if passing else None, top_id=candidates[0].id if candidates else None,
                          thresholds=gate.thresholds(), checks=chosen[0][2]["checks"] if chosen else [], warnings=warnings,
                          reference=_reference_spectrum(chosen[0][1]) if chosen else None)
    return finish(STATUS_DONE, _note(res, candidates, passing, tol_ppm, precursor_mz), candidates, summary)


def _candidate(cid, hit, verdict, n_spectra, metric):
    """The Candidate of one compound, from its best spectrum."""
    inchikey = _inchikey(hit)
    smiles = hit["smiles"] if has_structure(hit.get("smiles")) else None                  # the libraries also store "N/A" and the like
    failed = [f"{c['label']}: {c['detail']}" for c in verdict["checks"] if c["status"] == FAIL]
    evidence = {"gate_passed": verdict["passed"], "failed_checks": failed, "similarity": verdict["similarity"], "metric": metric,
                "matched_peaks": hit["matches"], "precursor_error_ppm": verdict["error"], "precursor_error_basis": verdict["error_basis"],
                "isotope": verdict["isotope"], "isotope_similarity": verdict["isotope_similarity"], "adduct_source": verdict["adduct_source"],
                "checks": verdict["checks"],
                "library": hit["library"], "accession": hit["accession"], "instrument": hit["instrument"],
                "collision_energy": hit["collision_energy"], "library_adduct": hit.get("adduct") or None,
                "library_precursor_mz": hit["precursor_mz"], "n_spectra": n_spectra}
    return Candidate(id=cid, layer="library", smiles=smiles, inchikey=inchikey, ik14=inchikey[:14] if inchikey else None, name=hit.get("name") or "",
                     formula=verdict["formula"], adduct=verdict["adduct"],
                     mz_error_ppm=verdict["error"] if verdict["error_from_formula"] else None, layer_score=verdict["similarity"],
                     layer_score_norm=min(1.0, max(0.0, verdict["similarity"])), evidence={"library": evidence})


def _note(res, candidates, passing, tol_ppm, precursor_mz):
    if not res["n_candidates"]:
        return f"No library spectrum within {tol_ppm:g} ppm of m/z {precursor_mz:.4f}."
    if not candidates:
        return f"{res['n_candidates']:,} library spectra compared; none scored {LIST_THRESHOLD:.2f} or higher."
    text = f"{res['n_kept']:,} library spectra scored {LIST_THRESHOLD:.2f} or higher ({len(candidates)} compounds listed)."
    if passing:
        ev = passing[0].evidence["library"]
        return text + f" Gate passed by {passing[0].name or passing[0].id} (similarity {ev['similarity']:.2f}, {ev['matched_peaks']} matched peaks)."
    top = candidates[0]
    return text + f" Gate not passed: the best hit {top.name or top.id} failed " + "; ".join(top.evidence["library"]["failed_checks"]) + "."
