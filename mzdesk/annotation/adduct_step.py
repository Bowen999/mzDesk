"""Step 1 of the ESI MS/MS Annotation workflow: which adduct is the precursor ion?

Adduct "Auto": every adduct of the precursor's polarity in the app's master list (spec.ADDUCTS_POS / ADDUCTS_NEG) is a hypothesis.
Its score is a relative likelihood, the product of five pieces of evidence, scaled so that all hypotheses sum to 1:

    score  ~  prior  x  MS2 factor  x  formula factor  x  isotope factor  x  mass factor

  prior           how often MassBank + GNPS annotate the adduct (LIBRARY_ADDUCT_COUNTS), flattened (PRIOR_EXPONENT) so that it only
                  breaks ties.
  MS2 factor      the adduct's diagnostic signals in the MS/MS spectrum (DIAGNOSTICS): one seen -> boost, none seen -> small penalty.
  formula factor  BUDDY's probability that its best formula under the adduct is the right one (formula_factor): if no formula explains
                  precursor and fragments once the adduct mass is taken off, the adduct is unlikely. formula_step.py runs BUDDY and
                  passes the result in (buddy_fits).
  isotope factor  fit of the MS1 isotope pattern to that best formula (BUDDY's isotope similarity). It partly overlaps the formula
                  factor, because BUDDY's own score already uses the isotope fit when MS1 is given.
  mass factor     how far the precursor is from that best formula's theoretical m/z, relative to the precursor tolerance
                  (mass_error_factor): a formula that only just fits the tolerance is a weaker explanation than one that fits well.
                  The "best formula" is the one that maximises BUDDY's probability x mass factor (formula_step._buddy_fit).
No MS1, no BUDDY or an adduct BUDDY does not support: that factor is 1 (no evidence).

An adduct chosen by the user is the only hypothesis: score 1, evidence "user". The ranking keeps the top RETAIN_MAX hypotheses whose
score is at least RETAIN_RATIO x the best, and is "ambiguous" when the best is less than AMBIGUOUS_RATIO x the second.

Peak lists are [(mz, intensity), ...], as mzdesk.shared.spectrum.parse_spectrum / api.common.spectrum return them (any intensity scale).
"""
from __future__ import annotations

import math
import time
from dataclasses import asdict, dataclass, field

from mzdesk.adducts import canonical, polarity

# ---------------------------------------------------------------- constants (owner's design unless a comment says otherwise)
RETAIN_MAX = 3              # keep at most the top 3 hypotheses ...
RETAIN_RATIO = 0.2          # ... whose score is at least 0.2 x the best
AMBIGUOUS_RATIO = 3.0       # best / second below 3: the evidence does not decide between the two

# Weak prior: how many MS/MS spectra carry each adduct in the local libraries, counted 2026-10-28 from the indexes
# (MassBank 2026.03: 113,030 MS/MS; GNPS 2026-03-29: 955,289 MS/MS without its 65,371 copies of MassBank records; canonical names).
LIBRARY_ADDUCT_COUNTS = {
    "[M+H]+": 529425, "[M+Na]+": 143712, "[M+NH4]+": 27674, "[M+H-H2O]+": 16239, "[M+K]+": 7923, "[M]+": 5969, "[2M+H]+": 1899,
    "[M-H]-": 192573, "[M+HCOO]-": 36630, "[M+Cl]-": 16332, "[M+CH3COO]-": 1806, "[2M-H]-": 2627, "[M-H-H2O]-": 97,
}
PRIOR_EXPONENT = 0.5        # square root of the library count: [M+H]+ vs [M+Na]+ stays about 2:1 (raw 3.7:1), so a signal can outvote it

# MS2 diagnostic signals. A signal is (kind, value, label): "loss" = peak at precursor - value (neutral loss), "ion" = peak at m/z value.
# An adduct's factor is `seen` if any of its signals is found, else `none_seen`. Library check of the numbers: the share of
# spectra in which ms2_evidence sees a signal, among 3,000 random MassBank + GNPS spectra of the adduct vs of [M+H]+ / [M-H]-
# (fragment tolerance 0.02 Da, peak >= 2 % of the base peak):
#   [M+NH4]+ 11 % vs 6 %     [M+HCOO]- 77 % vs 4 %     [M+Cl]- 60 % vs 2 %     [M+CH3COO]- 13 % vs 3 %     K+ ion 3 % vs 0 %
# so the measured ratios are 2 / 18 / 28 / 5; `seen` is capped at 10 and `none_seen` is weaker than measured (0.24 / 0.41 for formate /
# chloride) because user data differ from the libraries. NH3 is also lost by plain [M+H]+ of amines, hence only 2. The Na+ ion is
# never recorded (m/z < 50): it can only add, not subtract.
LOSS, ION = "loss", "ion"
DIAGNOSTICS = {                  # adduct: (signals, seen, none_seen)
    "[M+NH4]+":    (((LOSS, 17.0265, "loss of NH3"),), 2.0, 1.0),
    "[M+Na]+":     (((ION, 22.9892, "Na+"),), 5.0, 1.0),
    "[M+K]+":      (((ION, 38.9632, "K+"),), 5.0, 1.0),
    "[M+HCOO]-":   (((LOSS, 46.0055, "loss of HCOOH"), (ION, 44.9982, "HCOO-")), 10.0, 0.5),
    "[M+Cl]-":     (((LOSS, 35.9767, "loss of HCl"), (ION, 34.9689, "Cl-")), 10.0, 0.5),
    "[M+CH3COO]-": (((LOSS, 60.0211, "loss of CH3COOH"),), 5.0, 0.9),
}
# "Little neutral-loss fragmentation" (alkali adducts hold on to their cation): at most SPARSE_MAX_PEAKS fragment peaks. In the same
# library sample 62 % of [M+Na]+ spectra are sparse against 39 % of [M+H]+ (ratios 1.6 sparse, 0.6 not sparse); the [M+K]+ share
# (11 %) is skewed by one data set that records many peaks, other [M+K]+ data sets are sparse like [M+Na]+.
SPARSE_ADDUCTS = ("[M+Na]+", "[M+K]+")
SPARSE_FACTORS = (1.5, 0.7)      # (sparse, not sparse)
SPARSE_MAX_PEAKS = 3
SPARSE_MIN_REL = 5.0             # % of the base peak: a fragment peak for the sparseness count
SIGNAL_MIN_REL = 2.0             # % of the base peak: a diagnostic peak must be at least this strong
PRECURSOR_CLUSTER_DA = 1.5       # closer to the precursor than this = the precursor and its isotopes, not a fragment

# Isotope fit (BUDDY's MS1 isotope similarity of the adduct's best formula) -> factor in [ISOTOPE_FLOOR, 1], linear between the bounds.
ISOTOPE_GOOD = 0.95              # fit at or above this: no penalty
ISOTOPE_POOR = 0.80              # fit at or below this: the floor (the library gate's 0.80, here on BUDDY's scale; a poor fit lowers, never vetoes)
ISOTOPE_FLOOR = 0.25             # a hypothesis whose best formula fits the pattern poorly keeps a quarter of its score

# BUDDY's probability for the best formula (within the precursor tolerance) under the adduct -> factor FORMULA_FLOOR + (1 - FORMULA_FLOOR) x p.
# Weight: a weak, bounded evidence source. Library check (166 MassBank + GNPS spectra (40 each of [M+H]+, [M-H]-, [M+NH4]+,
# [M+HCOO]-, 6 of [M+Na]+) x every BUDDY-supported adduct of the polarity, qtof preset, 10 ppm): the best-formula probability is >= 0.05 for 90-100 % of the true adducts but also for 74-95 % of the
# wrong ones (median 0.85-0.98 vs 0.98-1.0), because neutral formulas of wrong adducts are often chemically fine too. It only separates
# when no formula fits (probability 0: about 1 in 4-5 wrong [M+H]+ / [M-H]- / [M+Na]+ hypotheses, up to 1 in 10 true ones). So the factor runs
# from 0.1 (probability 0) to 1 (probability 1): at most a 10-fold penalty, usually 1-2-fold.
FORMULA_FLOOR = 0.1

# Precursor mass error of the best formula -> factor in [MASS_FLOOR, 1], owner's design 2026-10-30: a Gaussian in ppm with sigma = tol_ppm / 2,
# lifted to MASS_FLOOR so that a bad mass is penalised, never vetoed. Relative to the user's tolerance, so a narrower window penalises
# the same fraction of it (10 ppm: x0.97 at 1.3 ppm, x0.85 at 3 ppm, x0.65 at 5 ppm, x0.27 at 9.2 ppm, x0.22 at the limit). Why: BUDDY weighs
# the mass error only mildly inside its window, so a wrong formula near the edge of the window (e.g. +9 ppm) can keep a high probability.
# The penalty is mild up to ~3 ppm, where the precursors of a well-calibrated high-resolution instrument usually fall. No instrument offset is
# modelled: that would be instrument-specific.
MASS_SIGMA_DIVISOR = 2.0        # sigma = tol_ppm / 2
MASS_FLOOR = 0.1                # at most a 10-fold penalty, like FORMULA_FLOOR


# ---------------------------------------------------------------- result types
@dataclass
class AdductHypothesis:
    """One adduct hypothesis (one table row).

    score      relative likelihood, 0-1, summing to 1 over the scored hypotheses (1.0 for an adduct the user chose)
    rank       1 = best
    retained   in the top RETAIN_MAX and score >= RETAIN_RATIO x the best: these go on to the formula and library steps
    ambiguous  retained and less than AMBIGUOUS_RATIO times weaker than the best (the best one itself is flagged when it has such a rival)
    prior / ms2_factor / formula_factor / isotope_factor / mass_factor   the multipliers (None for a user-chosen adduct)
    evidence   per source: "prior" {share, library_spectra, note}, "ms2" {factor, signals [{label, mz, found, rel_int}], note},
               "formula" {factor, probability, formula, ppm, note}, "isotope" {factor, similarity, formula, note},
               "mass" {factor, ppm, formula, note}; or "user" {note}
    """
    adduct: str
    score: float
    rank: int = 1
    retained: bool = True
    ambiguous: bool = False
    prior: float | None = None
    ms2_factor: float | None = None
    formula_factor: float | None = None
    isotope_factor: float | None = None
    mass_factor: float | None = None
    evidence: dict = field(default_factory=dict)

    def to_dict(self):
        return asdict(self)


@dataclass
class AdductResult:
    """Step 1 result. mode: "auto" | "user". hypotheses: every scored adduct, best first (retained ones come first)."""
    mode: str
    hypotheses: list
    ambiguous: bool = False
    status: str = "done"
    note: str = ""
    seconds: float = 0.0

    @property
    def retained(self):
        return [h for h in self.hypotheses if h.retained]

    @property
    def retained_adducts(self):
        """Canonical names of the retained hypotheses, best first (the adducts a library hit may have, and BUDDY runs for)."""
        return [h.adduct for h in self.retained]

    def to_dict(self):
        return {"mode": self.mode, "ambiguous": self.ambiguous, "status": self.status, "note": self.note, "seconds": self.seconds,
                "retained_adducts": self.retained_adducts, "hypotheses": [h.to_dict() for h in self.hypotheses]}


# ---------------------------------------------------------------- evidence
def clean_peaks(peaks):
    """[(mz, intensity)] with finite positive values only, as floats, sorted by m/z; None / empty -> []. (All three steps use this.)"""
    return sorted((float(mz), float(inten)) for mz, inten in peaks or [] if math.isfinite(mz) and math.isfinite(inten) and mz > 0 and inten > 0)


def master_adducts(positive):
    """The adducts of the polarity in the app's master list (spec.ADDUCTS_POS / ADDUCTS_NEG), canonical names, without "Auto"."""
    from mzdesk import spec          # imported here: spec pulls in the engine tables

    return [canonical(a) for a in (spec.ADDUCTS_POS if positive else spec.ADDUCTS_NEG) if canonical(a)]


def prior_shares(universe):
    """{adduct: share} of the library frequency^PRIOR_EXPONENT over `universe`, summing to 1 (an adduct the table lacks counts as 1)."""
    weights = {a: LIBRARY_ADDUCT_COUNTS.get(a, 1) ** PRIOR_EXPONENT for a in universe}
    total = sum(weights.values())
    return {a: w / total for a, w in weights.items()}


def _nearest(peaks, target, tol):
    """The most intense peak within `tol` of `target`: (mz, % of base peak) or None. `peaks` = [(mz, % of base peak)]."""
    near = [p for p in peaks if abs(p[0] - target) <= tol and p[1] >= SIGNAL_MIN_REL]
    return max(near, key=lambda p: p[1]) if near else None


def ms2_evidence(adduct, precursor_mz, ms2_peaks, tol_da=0.02):
    """MS2 diagnostic evidence for one adduct -> {"factor", "signals": [{"label", "mz", "found", "rel_int"}], "note"}.

    Needs at least one fragment peak (a peak more than PRECURSOR_CLUSTER_DA below the precursor); a spectrum without one cannot show
    the absence of a signal, so the factor is 1.
    """
    signals, seen_f, none_f = DIAGNOSTICS.get(adduct, ((), 1.0, 1.0))
    sparse_rule = adduct in SPARSE_ADDUCTS
    if not signals and not sparse_rule:
        return {"factor": 1.0, "signals": [], "note": "no diagnostic signal is defined for this adduct"}
    base = max((i for _, i in ms2_peaks), default=0.0)
    fragments = [(mz, i / base * 100.0) for mz, i in ms2_peaks if mz < precursor_mz - PRECURSOR_CLUSTER_DA] if base > 0 else []
    if not fragments:
        return {"factor": 1.0, "signals": [], "note": "the MS2 spectrum has no fragment peak, so its signals cannot be checked"}

    found = []
    for kind, value, label in signals:
        target = precursor_mz - value if kind == LOSS else value
        hit = _nearest(fragments, target, tol_da)
        found.append({"label": f"{label} ({value})", "mz": round(target, 4), "found": hit is not None,
                      "rel_int": None if hit is None else round(hit[1], 1)})
    factor = (seen_f if any(s["found"] for s in found) else none_f) if signals else 1.0
    notes = []
    if signals:
        seen = [f"{s['label']} at {s['rel_int']:g} %" for s in found if s["found"]]
        notes.append("seen: " + "; ".join(seen) if seen else "none of " + ", ".join(s["label"] for s in found) + " seen")
    if sparse_rule:
        n = sum(1 for _, rel in fragments if rel >= SPARSE_MIN_REL)
        sparse = n <= SPARSE_MAX_PEAKS
        factor *= SPARSE_FACTORS[0 if sparse else 1]
        notes.append(f"{n} fragment peaks >= {SPARSE_MIN_REL:g} % ({'few' if sparse else 'many'}: alkali adducts fragment little)")
    return {"factor": factor, "signals": found, "note": "; ".join(notes)}


def isotope_factor(similarity):
    """Factor in [ISOTOPE_FLOOR, 1] from BUDDY's isotope similarity of an adduct's best formula; None (no evidence) -> 1."""
    if similarity is None:
        return 1.0
    t = (similarity - ISOTOPE_POOR) / (ISOTOPE_GOOD - ISOTOPE_POOR)
    return ISOTOPE_FLOOR + (1.0 - ISOTOPE_FLOOR) * min(1.0, max(0.0, t))


def formula_factor(probability):
    """Factor in [FORMULA_FLOOR, 1] from BUDDY's probability (0-1) for the best formula under an adduct; None (no evidence) -> 1."""
    if probability is None:
        return 1.0
    return FORMULA_FLOOR + (1.0 - FORMULA_FLOOR) * min(1.0, max(0.0, probability))


def mass_error_factor(ppm, tol_ppm):
    """Factor in [MASS_FLOOR, 1] from the precursor mass error `ppm` of a formula, relative to the tolerance `tol_ppm`; None (no formula) -> 1.

    1 at 0 ppm, falling with |ppm| (Gaussian, sigma = tol_ppm / MASS_SIGMA_DIVISOR) towards MASS_FLOOR, which it reaches well beyond the limit.
    """
    if ppm is None or not tol_ppm or tol_ppm <= 0:
        return 1.0
    z = abs(ppm) / (tol_ppm / MASS_SIGMA_DIVISOR)
    return MASS_FLOOR + (1.0 - MASS_FLOOR) * math.exp(-0.5 * z * z)


# ---------------------------------------------------------------- scoring
def _user_adduct(adduct, positive):
    """Canonical name of the adduct the user chose; None for Auto. ValueError for a name that is unknown or of the wrong polarity."""
    if adduct is None or str(adduct).strip().lower() in ("", "auto"):
        return None
    name = canonical(adduct)
    if name is None:
        raise ValueError(f"Unknown adduct notation: {adduct}.")
    if polarity(name) != ("Positive" if positive else "Negative"):
        raise ValueError(f"{adduct} does not match the {'positive' if positive else 'negative'} ion mode.")
    return name


def score_adducts(precursor_mz, positive, ms2_peaks, *, adduct=None, tol_da=0.02, buddy_fits=None, candidates=None):
    """Score the adduct hypotheses -> AdductResult.

    precursor_mz  observed precursor m/z
    positive      True = positive ion mode
    ms2_peaks     [(mz, intensity)]
    adduct        None / "Auto" = score the hypotheses; a name (any spelling) = use that adduct as is
    tol_da        fragment tolerance for the MS2 signals (Da)
    buddy_fits    {adduct: {"formula": str | None, "prob": float | None, "ppm": float | None, "isotope_sim": float | None,
                  "mass_factor": float | None, "note": str}}:
                  what BUDDY found for the best formula under each adduct, filled by formula_step for the adducts it ran for;
                  an adduct that is missing, or whose prob / isotope_sim is None, counts as "no evidence" for that factor
    candidates    the hypotheses to score (default: master_adducts(positive))
    """
    started = time.perf_counter()
    chosen = _user_adduct(adduct, positive)
    if chosen:
        hyp = AdductHypothesis(adduct=chosen, score=1.0, evidence={"user": {"note": "Chosen in the form; used as is."}})
        return AdductResult("user", [hyp], False, "done", f"Adduct {chosen} was chosen in the form; used as is.", time.perf_counter() - started)

    universe = [c for c in (canonical(a) for a in (candidates or master_adducts(positive))) if c]
    shares = prior_shares(universe)
    rows = []
    for name in universe:
        ms2 = ms2_evidence(name, precursor_mz, ms2_peaks, tol_da)
        fit = (buddy_fits or {}).get(name) or {}
        prob, sim, found = fit.get("prob"), fit.get("isotope_sim"), fit.get("formula")
        iso, form = isotope_factor(sim), formula_factor(prob)
        mass = fit.get("mass_factor")
        mass = 1.0 if mass is None else mass
        counts = LIBRARY_ADDUCT_COUNTS.get(name)
        if not fit:
            form_note = iso_note = mass_note = "not used: BUDDY was not run for this adduct"
        elif found is None:
            form_note = iso_note = mass_note = fit.get("note") or "not used: no BUDDY formula"
        else:
            form_note = (f"BUDDY's best formula {found}: probability {prob:.3f}" + (f", {fit['ppm']:+.1f} ppm" if fit.get("ppm") is not None else "")
                         if prob is not None else "not used: BUDDY gave no probability")
            iso_note = f"best BUDDY formula {found}: isotope similarity {sim:.3f}" if sim is not None else "not used: no usable MS1 isotope pattern"
            mass_note = (f"{found}: mass error {fit['ppm']:+.1f} ppm \u2192 \u00d7{mass:.2f}" if fit.get("ppm") is not None
                         else "not used: no mass error")
        evidence = {
            "prior": {"share": round(shares[name], 4), "library_spectra": counts,
                      "note": f"{counts:,} library spectra" if counts else "not in the library table"},
            "ms2": ms2,
            "formula": {"factor": round(form, 4), "probability": fit.get("prob"), "formula": fit.get("formula"), "ppm": fit.get("ppm"),
                        "note": form_note},
            "isotope": {"factor": round(iso, 4), "similarity": fit.get("isotope_sim"), "formula": fit.get("formula"), "note": iso_note},
            "mass": {"factor": round(mass, 4), "ppm": fit.get("ppm"), "formula": fit.get("formula"), "note": mass_note},
        }
        rows.append((shares[name] * ms2["factor"] * form * iso * mass, name, shares[name], ms2["factor"], form, iso, mass, evidence))

    total = sum(r[0] for r in rows)
    ranked = sorted(rows, key=lambda r: -r[0])               # a stable sort: ties keep the master-list order
    best = ranked[0][0] / total
    hypotheses = []
    for rank, (raw, name, share, ms2_f, form_f, iso_f, mass_f, evidence) in enumerate(ranked, start=1):
        score = raw / total
        retained = rank <= RETAIN_MAX and score >= RETAIN_RATIO * best
        # Keep the full score for downstream formula selection; formatting is display-only.
        hypotheses.append(AdductHypothesis(adduct=name, score=score, rank=rank, retained=retained,
                                           ambiguous=retained and (rank == 1 or best < AMBIGUOUS_RATIO * score),
                                           prior=round(share, 4), ms2_factor=round(ms2_f, 4), formula_factor=round(form_f, 4),
                                           isotope_factor=round(iso_f, 4), mass_factor=round(mass_f, 4), evidence=evidence))
    ambiguous = len(hypotheses) > 1 and ranked[0][0] < AMBIGUOUS_RATIO * ranked[1][0]
    if not ambiguous:
        hypotheses[0].ambiguous = False                      # the best one is flagged only when it has a rival within the ratio
    kept = [h for h in hypotheses if h.retained]
    note = f"Scored {len(hypotheses)} adducts; kept " + ", ".join(f"{h.adduct} ({h.score:.2f})" for h in kept) + "."
    if ambiguous:
        note += f" Ambiguous: {hypotheses[0].adduct} is less than {AMBIGUOUS_RATIO:g} times more likely than {hypotheses[1].adduct}."
    return AdductResult("auto", hypotheses, ambiguous, "done", note, time.perf_counter() - started)
