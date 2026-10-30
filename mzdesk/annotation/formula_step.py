"""Step 2 of the ESI MS/MS Annotation workflow: molecular formulas from BUDDY (msbuddy, in process), coupled with step 1.

Why one entry function: the isotope evidence of step 1 is BUDDY's isotope fit of each adduct's best formula, and step 2 needs the
adduct ranking to know which adducts get a formula list. run_adduct_and_formula() therefore does both:

  1. pre-score the adducts (adduct_step.score_adducts: prior + MS2 signals);
  2. run BUDDY for every retained adduct and feed back what it found for the best formula (its probability, and its isotope fit when
     a usable MS1 pattern is given), until every retained adduct has been evaluated (a hypothesis that is not yet evaluated must not
     win over an evaluated one just because nothing lowered it);
  3. the shared formula list for steps 3-5: BUDDY's top-5 formulas for the best adduct, for the second one if its score is at least
     SECOND_ADDUCT_RATIO x the best, and for any other retained adduct whose own best formula has a BUDDY probability >= EXTRA_ADDUCT_PROB
     (auto mode; formulas that miss the user's precursor tolerance are dropped). Within an adduct the formulas are ordered by BUDDY's
     probability x the mass factor (adduct_step.mass_error_factor: a formula near the tolerance limit is demoted), then the list is
     interleaved across these adducts (A#1, B#1, C#1, A#2, ...), so the consumers that take only the first formulas (MS-FINDER: 6,
     FRIGID: 1) see every target adduct. The "rank" of a formula stays BUDDY's own rank.

Mind that [M+NH4]+ and [M+H-H2O]+ give the same ion composition as [M+H]+ of another neutral (M+NH3, M-H2O): BUDDY's score and the
isotope fit are then identical, and only the MS2 loss (NH3 17.0265) and the prior can tell these hypotheses apart.

BUDDY gets adducts.engine_name(adduct, "buddy"), never a raw label. msbuddy_engine.annotate_formula takes its own lock and caches.
BUDDY missing or failing is not an error of the workflow: the formula list is empty, the step says why, and the library formula check
falls back to the precursor mass.

Peak lists are [(mz, intensity), ...] (parse_spectrum / api.common.spectrum output). ms1_peaks is the pasted MS1 isotope pattern (or a
wider MS1 list), None / [] when not given; it is used only if it has at least MIN_MS1_PEAKS peaks and one of them is at the precursor.
"""
from __future__ import annotations

import re
import time
from dataclasses import asdict, dataclass, field

from mzdesk import adducts
from mzdesk.annotation import adduct_step
from mzdesk.annotation.adduct_step import clean_peaks

SECOND_ADDUCT_RATIO = 0.5       # owner's design: the second adduct also gets formulas if its score >= 0.5 x the best
EXTRA_ADDUCT_PROB = 0.9         # owner's design 2026-10-30: any retained adduct whose best BUDDY formula has prob >= 0.9 also gets formulas
MIN_MS1_PEAKS = 2               # an isotope pattern is the monoisotopic peak and at least one more (BUDDY needs two as well)
# What BUDDY can produce: the element ranges of msbuddy 0.3.9 MsbuddyConfig; the halogens only with halogen=True. A formula outside this
# space cannot be in BUDDY's list, so a list says nothing about it (library_step then judges it by the precursor mass).
BUDDY_LIMITS = {"C": 80, "H": 150, "N": 20, "O": 30, "P": 10, "S": 15}
BUDDY_HALOGEN_LIMITS = {"F": 20, "Cl": 15, "Br": 10, "I": 10}


# ---------------------------------------------------------------- result types
@dataclass
class FormulaCandidate:
    """One BUDDY formula (one table row). The field names match what merge.formula_only_candidates reads from a formula dict.

    formula      neutral molecular formula M, as BUDDY writes it
    adduct       canonical adduct the formula was computed for
    rank         BUDDY's rank within that adduct (1 = best)
    prob         BUDDY's estimated probability that the formula is the right one (its score), 0-1
    ppm          (observed - theoretical) / theoretical x 1e6 of the precursor m/z
    mz_theo      theoretical m/z of formula + adduct
    isotope_sim  BUDDY's MS1 isotope similarity of the formula (1 - L1 distance of the sum-normalised first isotope peaks); None when
                 no MS1 pattern was used
    mass_factor  adduct_step.mass_error_factor(ppm, tolerance): 1 for a perfect mass, falling to MASS_FLOOR at the tolerance limit. Within
                 an adduct the list is ordered by prob x mass_factor; prob stays BUDDY's own.
    """
    formula: str
    adduct: str
    rank: int
    prob: float | None = None
    ppm: float | None = None
    mz_theo: float | None = None
    isotope_sim: float | None = None
    mass_factor: float | None = None

    def to_dict(self):
        return asdict(self)


@dataclass
class FormulaResult:
    """Step 2 result: the shared formula list, interleaved across the adducts by position (best adduct first); within an adduct ordered by prob x mass_factor.

    adducts   the adducts a list was asked for, in adduct-score order (the best one, the second if it scored high enough, and any other
              retained adduct with a confident BUDDY formula)
    status    done | skipped (BUDDY does not apply: unsupported adduct, precursor outside its range) | failed (BUDDY missing or crashed)
    ms1_used  the MS1 list had a usable isotope pattern and was given to BUDDY
    An empty list (any status) means: the library formula check uses the precursor mass instead.
    """
    formulas: list = field(default_factory=list)
    adducts: list = field(default_factory=list)
    status: str = "done"
    note: str = ""
    seconds: float = 0.0
    ms1_used: bool = False

    @property
    def formula_strings(self):
        """Distinct neutral formulas of the list, best first."""
        return list(dict.fromkeys(f.formula for f in self.formulas))

    def to_dict(self):
        return {"status": self.status, "note": self.note, "seconds": self.seconds, "ms1_used": self.ms1_used,
                "adducts": self.adducts, "formulas": [f.to_dict() for f in self.formulas]}


# ---------------------------------------------------------------- MS1 and BUDDY
def ms1_usable(ms1_peaks, precursor_mz, tol_ppm):
    """(True, "") when the MS1 list can serve as an isotope pattern; else (False, reason). (False, "") when no MS1 was given."""
    peaks = clean_peaks(ms1_peaks)
    if not peaks:
        return False, ""
    if len(peaks) < MIN_MS1_PEAKS:
        return False, "The MS1 list has a single peak; an isotope pattern needs at least two, so MS1 is not used."
    window = precursor_mz * tol_ppm * 1e-6
    if not any(abs(mz - precursor_mz) <= window for mz, _ in peaks):
        return False, f"The MS1 list has no peak within {tol_ppm:g} ppm of the precursor m/z, so MS1 is not used."
    return True, ""


class _Buddy:
    """BUDDY runs of one query: one per adduct, each run once. run(adduct) -> (status, payload):
    ("ok", annotate_formula result) | ("unsupported", text) | ("range", text) | ("error", text)."""

    def __init__(self, precursor_mz, positive, ms2, ms1, tol_ppm, halogen):
        self.args = (precursor_mz, positive, tuple(ms2), tuple(ms1), tol_ppm, halogen)
        self._runs = {}
        self._broken = None

    def run(self, adduct):
        if adduct not in self._runs:
            self._runs[adduct] = self._run(adduct)
        return self._runs[adduct]

    def _run(self, adduct):
        if self._broken:
            return "error", self._broken
        try:
            from mzdesk.ms1 import msbuddy_engine as bd

            name = adducts.engine_name(adduct, "buddy")
            if name is None:
                return "unsupported", f"BUDDY does not support the adduct {adduct}."
            precursor_mz, positive, ms2, ms1, tol_ppm, halogen = self.args
            # the tightest instrument preset whose MS1 window covers the user's precursor tolerance
            fits = [(tols[0], key) for key, tols in bd.INSTRUMENT_TOLS.items() if tols[0] >= tol_ppm]
            preset = min(fits)[1] if fits else max((tols[0], key) for key, tols in bd.INSTRUMENT_TOLS.items())[1]
            return "ok", bd.annotate_formula(precursor_mz, positive, name, None, "", ms2, ms1, preset, halogen)
        except ValueError as exc:                      # BUDDY only annotates 1 < m/z < 1500
            return "range", str(exc)
        except Exception as exc:                       # missing package, crash: do not try again for the other adducts
            self._broken = f"{type(exc).__name__}: {exc}"
            return "error", self._broken


def buddy_covers(formula, halogen=False):
    """True when `formula` lies in BUDDY's search space (C H N O P S, plus F Cl Br I with halogen=True, within msbuddy's element ranges).

    An isotope-labelled formula (square brackets) is outside it.
    """
    if "[" in (formula or ""):
        return False
    limits = {**BUDDY_LIMITS, **(BUDDY_HALOGEN_LIMITS if halogen else {})}
    counts = {}
    for element, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula or ""):
        counts[element] = counts.get(element, 0) + int(n or 1)
    return bool(counts) and all(limits.get(element, -1) >= n for element, n in counts.items())


def _buddy_fit(run, tol_ppm):
    """What adduct_step needs from one adduct run: {"formula", "prob", "ppm", "isotope_sim", "mass_factor", "note"} of the best formula within tol_ppm.

    "Best" = the highest BUDDY probability x mass factor (adduct_step.mass_error_factor of its ppm); "prob" stays BUDDY's own.

    A run that gave no formula at all within the tolerance is evidence in itself (probability 0: nothing explains the precursor under
    this adduct); a run that could not be made (unsupported adduct, BUDDY failing) is not (all None).
    """
    status, payload = run
    if status != "ok":
        return {"formula": None, "prob": None, "ppm": None, "isotope_sim": None, "note": f"not used: no BUDDY formula ({payload})"}
    within = [c for c in payload["candidates"] if c["ppm"] is not None and abs(c["ppm"]) <= tol_ppm]
    if not within:
        return {"formula": None, "prob": 0.0, "ppm": None, "isotope_sim": None, "mass_factor": None,
                "note": f"BUDDY found no formula within {tol_ppm:g} ppm"}
    # the best formula maximises probability x mass factor (the first of equals: BUDDY's order)
    best = max(within, key=lambda c: (c["prob"] or 0.0) * adduct_step.mass_error_factor(c["ppm"], tol_ppm))
    return {"formula": best["formula"], "prob": best["prob"], "ppm": best["ppm"], "isotope_sim": best["isotope_sim"],
            "mass_factor": adduct_step.mass_error_factor(best["ppm"], tol_ppm), "note": ""}


def _candidates(adduct, payload, tol_ppm):
    """FormulaCandidates of one BUDDY run that fit the user's precursor tolerance, best first by prob x mass factor (BUDDY's order on ties)."""
    found = [FormulaCandidate(formula=c["formula"], adduct=adduct, rank=c["rank"], prob=c["prob"], ppm=c["ppm"], mz_theo=c["mz_theo"],
                              isotope_sim=c["isotope_sim"], mass_factor=adduct_step.mass_error_factor(c["ppm"], tol_ppm))
             for c in payload["candidates"] if c["ppm"] is not None and abs(c["ppm"]) <= tol_ppm]
    return sorted(found, key=lambda f: -(f.prob or 0.0) * f.mass_factor)          # a stable sort


# ---------------------------------------------------------------- steps 1 + 2
def run_adduct_and_formula(precursor_mz, positive, ms2_peaks, ms1_peaks=None, *, adduct=None, tol_ppm=10.0, tol_da=0.02, halogen=False):
    """Steps 1 and 2 -> (AdductResult, FormulaResult).

    precursor_mz  observed precursor m/z
    positive      True = positive ion mode
    ms2_peaks     [(mz, intensity)]
    ms1_peaks     [(mz, intensity)] isotope pattern, or None
    adduct        None / "Auto" = rank the adducts of the polarity; a name (any spelling) = use it as is. ValueError if it is unknown
                  or of the other polarity.
    tol_ppm       precursor tolerance (ppm): picks BUDDY's instrument preset and drops formulas that miss it
    tol_da        fragment tolerance (Da) for the MS2 diagnostic signals
    halogen       let BUDDY consider Cl / Br / F / I (as the Formula page's "Include halogens")
    """
    started = time.perf_counter()
    ms2, ms1 = clean_peaks(ms2_peaks), clean_peaks(ms1_peaks)
    use_ms1, ms1_note = ms1_usable(ms1, precursor_mz, tol_ppm)
    buddy = _Buddy(precursor_mz, positive, ms2, ms1 if use_ms1 else [], tol_ppm, halogen)

    def score(fits=None):
        return adduct_step.score_adducts(precursor_mz, positive, ms2, adduct=adduct, tol_da=tol_da, buddy_fits=fits)

    adduct_result = score()
    if adduct_result.mode == "auto":
        fits = {}
        while True:
            todo = [a for a in adduct_result.retained_adducts if a not in fits]
            if not todo:
                break
            for name in todo:
                fits[name] = _buddy_fit(buddy.run(name), tol_ppm)
            adduct_result = score(fits)
    if ms1_note:
        adduct_result.note += " " + ms1_note
    adduct_result.seconds = time.perf_counter() - started

    started = time.perf_counter()
    kept = adduct_result.retained
    wanted = {kept[0].adduct}
    if len(kept) > 1 and kept[1].score >= SECOND_ADDUCT_RATIO * kept[0].score:
        wanted.add(kept[1].adduct)
    if adduct_result.mode == "auto":                    # these runs are cached from the evaluation above: no extra BUDDY cost
        wanted |= {h.adduct for h in kept if (_buddy_fit(buddy.run(h.adduct), tol_ppm)["prob"] or 0.0) >= EXTRA_ADDUCT_PROB}
    targets = [h.adduct for h in kept if h.adduct in wanted]
    lists, problems = [], []
    for name in targets:
        status, payload = buddy.run(name)
        if status == "ok":
            lists.append(_candidates(name, payload, tol_ppm))
        else:
            problems.append((name, status, payload))
    # owner's design 2026-10-30: round-robin over the adducts (A#1, B#1, C#1, A#2, ...), not adduct after adduct; each adduct's list is
    # already ordered by probability x mass factor (_candidates)
    formulas = [lst[i] for i in range(max(map(len, lists), default=0)) for lst in lists if i < len(lst)]
    result = FormulaResult(formulas=formulas, adducts=targets, ms1_used=use_ms1)
    result.status, result.note = _summary(result, targets, problems, use_ms1, tol_ppm)
    result.seconds = time.perf_counter() - started
    return adduct_result, result


def _join(items, word):
    """'a', 'a and b', 'a, b and c'."""
    return items[0] if len(items) == 1 else f"{', '.join(items[:-1])} {word} {items[-1]}"


def _summary(result, targets, problems, use_ms1, tol_ppm):
    """(status, note) of the formula step."""
    fallback = "there is no formula list, so the library formula check uses the precursor mass"
    if len(problems) == len(targets):
        status, text = problems[0][1:]
        if status == "error":
            return "failed", f"BUDDY failed ({text}); {fallback}."
        return "skipped", f"{text.rstrip('.')}; {fallback}."
    asked = [a for a in targets if a not in {p[0] for p in problems}]
    how = "BUDDY" + (" with the MS1 isotope pattern" if use_ms1 else "")
    if result.formulas:
        parts = [f"{n} formula{'' if n == 1 else 's'} for {a}" for a in asked for n in [sum(f.adduct == a for f in result.formulas)]]
        note = f"{how}: {_join(parts, 'and')}, interleaved by rank." if len(parts) > 1 else f"{how}: {parts[0]}."
    else:
        note = f"{how} found no formula within {tol_ppm:g} ppm for {_join(asked, 'or')}; {fallback}."
    return "done", note + "".join(f" {name}: {text}" for name, _, text in problems)
