"""Step 6b: re-rank the merged candidates: the layers' evidence decides, detectability breaks near-ties.

    final_score = base_score x f_omics x f_ionizable x f_adduct

base_score   the evidence the layers gave: the best layer_score_norm over the contributing layers, plus a bonus for every extra
             layer that found the same structure, capped; a pure FRIGID hypothesis is multiplied by a factor below 1.
f_omics      could the chosen omics platform detect this molecule?        lipid_metabolite_classifier.analyze
f_ionizable  does it respond in ESI of the precursor polarity?            ionization_rules (the ESI part of assess)
f_adduct     is the observed adduct plausible for this structure?         adduct_rules.predict

The factors are a tie-breaker, not a second opinion (owner decision from a 70-spectrum MS-FINDER evaluation: layer order 48 top-1 correct,
full re-rank by final_score 43, tie-breaker 49). Candidates whose base_score is within RerankConfig.tiebreak_window (10 %) of the best base_score
form the window and are ordered by final_score; the candidates outside the window follow in base_score order, so a candidate far behind on the
layers' evidence is never promoted above the window, however detectable it is. Outside the window the detect details are still computed and
shown, but final_score is set to base_score (the factors did not take part in the order) and detect["reasons"] gets a fourth line saying so.
The order is therefore given by `rank`, not by final_score alone.

Re-ranking only sorts. It never removes a candidate and never changes its layer or level. A dimension that cannot be evaluated
(no valid SMILES, no adduct, or a rule that fails on the structure) counts RerankConfig.not_evaluated and says so in detect["reasons"].

The adduct rules need a mobile-phase context the workflow form does not ask for. `additive` and `alkali` default to "unknown", and
unknown never lowers a score: the adduct is then scored under every additive that can offer it and the best result counts (an NH4+
adduct needs an ammonium additive, a formate or acetate adduct that anion, any other adduct can appear under all of them; the
search stops at the first "high", which cannot be beaten), and sodium / potassium count as present. Adduct names are made canonical
first ([M+FA-H]- is [M+HCOO]-, mzdesk.adducts); a name the adduct rules do not list reads as "not expected".

An observed CCS (optional, ccs_step.py) can move candidates after this re-rank: apply_ccs() below. It is supporting evidence only, so it never
changes a score or a level; see its docstring for the rule.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass
from functools import lru_cache

from mzdesk.adducts import canonical
from mzdesk.annotation.merge import has_structure
from mzdesk.platform_fit import adduct_rules as ar
from mzdesk.platform_fit import display_en as en
from mzdesk.platform_fit import ionization_rules as ir
from mzdesk.platform_fit import lipid_metabolite_classifier as clf


# ---------------------------------------------------------------- constants (the owner's start values; change them here only)
@dataclass(frozen=True)
class RerankConfig:
    """Start values of the re-rank. They are the owner's design, not fitted to any data set (never to CASMI, which leaks).

    f_omics: could the platform the user chose detect the molecule? (lipid_metabolite_classifier)
                          classifier says  lipid-type    metabolite-type   both        neither
        Lipidomics                         match         mismatch          match       neither
        Metabolomics                       mismatch      match             match       neither
        Either / unknown                   match         match             match       neither
      omics_match      detectable by the chosen platform ("both" counts), or, with the platform unknown, by any
      omics_mismatch   detectable only by the other platform
      omics_neither    detectable by neither (outside the 50-2000 Da LC-MS window, or no carbon)
    f_ionizable (ESI response in the precursor polarity)
      ion_strong / ion_moderate / ion_weak
    f_adduct (likelihood of the observed adduct)
      adduct_high / adduct_medium / adduct_low, and adduct_not_expected when the rules do not list the adduct at all
    order
      tiebreak_window   a candidate is in the window when base_score >= (1 - tiebreak_window) x the best base_score; only the window is
                        ordered by final_score, the rest stays in base_score order behind it (0 = only ties of the best base_score)
    base score
      layer_bonus      added for every layer beyond the first that found the same structure (the layers are independent evidence)
      base_cap          the base score never exceeds this
      hypothesis_factor multiplier for a structure that only FRIGID proposed (formula-level evidence, Candidate.hypothesis)
    not_evaluated       multiplier for each dimension that could not be evaluated (no valid SMILES, no adduct, or the rule failed), so
                        a candidate without a structure gets this three times
    """
    omics_match: float = 1.0
    omics_mismatch: float = 0.5
    omics_neither: float = 0.2
    ion_strong: float = 1.0
    ion_moderate: float = 0.8
    ion_weak: float = 0.5
    adduct_high: float = 1.0
    adduct_medium: float = 0.75
    adduct_low: float = 0.4
    adduct_not_expected: float = 0.2
    layer_bonus: float = 0.05
    base_cap: float = 1.0
    hypothesis_factor: float = 0.7
    not_evaluated: float = 0.8
    tiebreak_window: float = 0.10


DEFAULT_PH = 2.7                    # 0.1 % formic acid, the Ionization page's default
UNKNOWN = "unknown"

# omics the user can choose (the form's "Either / unknown" | "Metabolomics" | "Lipidomics")
OMICS_EITHER, OMICS_METABOLOMICS, OMICS_LIPIDOMICS = "either", "metabolomics", "lipidomics"
OMICS_LABELS = {OMICS_EITHER: "Either / unknown", OMICS_METABOLOMICS: "Metabolomics", OMICS_LIPIDOMICS: "Lipidomics"}
# what the classifier says about a molecule (detect["omics"])
LIPIDOME, METABOLOME, BOTH, NEITHER = "lipidome", "metabolome", "both", "neither"
CLASS_OF = {"Lipidome": LIPIDOME, "Metabolome": METABOLOME, "Both": BOTH, "Neither": NEITHER}      # display_en.CATEGORY_EN values
NOT_EVALUATED, NOT_EXPECTED = "not evaluated", "not expected"      # "not expected": the adduct rules do not list this adduct
ADDUCT_RANK = {NOT_EXPECTED: 0, "low": 1, "medium": 2, "high": 3}


# ---------------------------------------------------------------- the three platform_fit checks (cached: a structure is scored once)
@lru_cache(maxsize=2048)
def _omics_check(smiles):
    """(lipidome | metabolome | both | neither, one-line verdict) from the lipid / metabolite rules; (None, "") if unparsable."""
    res = clf.analyze(smiles)
    cls = CLASS_OF.get(en.category_en(res.category))
    return cls, (en.tile_reasons_en(res)[0][0] if cls else "")


@lru_cache(maxsize=2048)
def _ion_check(smiles, sign, ph):
    """(strong | moderate | weak, the rule's reason) of the ESI response in the polarity `sign`; (None, error text) if it cannot run.

    The same verdict as ionization_rules.assess(smiles, ph).pos / .neg, without the EI half of assess (NEIMS and the GC rules).
    """
    try:
        profile = ir.profile(smiles)
        r = ir.esi_positive(profile, ph) if sign == "+" else ir.esi_negative(profile, ph)
    except Exception as exc:        # an odd structure must not stop the ranking; the reason says what happened
        return None, f"{type(exc).__name__}: {exc}"
    return r.level, r.reason


@lru_cache(maxsize=2048)
def _predicted_ions(smiles, additive, alkali):
    """((adduct, likelihood, score, reason), ...) of both polarities, strongest first, and an error text ("" when it ran)."""
    try:
        pred = ar.predict(smiles, "Both", additive, alkali)
    except Exception as exc:        # same boundary as _ion_check
        return (), f"{type(exc).__name__}: {exc}"
    return tuple((i.adduct, i.likelihood, i.score, i.reason) for i in pred.all_ions), ""


def _adduct_check(smiles, adduct, additive, alkali):
    """(high | medium | low | not expected, rule score, reason) of the observed adduct under one additive; (None, 0, error text) if the rules cannot run."""
    ions, error = _predicted_ions(smiles, additive, alkali)
    if error:
        return None, 0.0, error
    for name, likelihood, score, reason in ions:
        if name == adduct:
            return likelihood, score, reason
    return NOT_EXPECTED, 0.0, "the adduct rules do not list it for this structure"


def _additives_for(adduct):
    """Additives under which the adduct rules can list `adduct`, in adduct_rules.ADDITIVES order (formic acid first)."""
    adds = [ar.additive(name) for name in ar.ADDITIVES]
    if "NH4" in adduct:                                      # NH4+ comes from an ammonium additive
        adds = [a for a in adds if a.ammonium]
    elif adduct in ("[M+HCOO]-", "[M+CH3COO]-"):             # formate / acetate comes from the additive's own anion
        adds = [a for a in adds if f"[M+{a.anion}]-" == adduct]
    return [a.name for a in adds]


def _best_adduct(smiles, adduct, additive, alkali):
    """(high | medium | low | not expected, reason) of the observed adduct: under the given additive, or for "unknown" under the best one.

    "Best" is the highest likelihood, then the highest rule score; the search stops at the first "high". (None, error text) if the rules cannot run.
    """
    best, tried = None, []
    for context in (_additives_for(adduct) if additive == UNKNOWN else [additive]):
        likelihood, score, reason = _adduct_check(smiles, adduct, context, alkali)
        if likelihood is None:
            return None, reason
        tried.append(context)
        if best is None or (ADDUCT_RANK[likelihood], score) > best[0]:
            best = ((ADDUCT_RANK[likelihood], score), likelihood, context, score, reason)
        if likelihood == "high":
            break
    _, likelihood, context, score, reason = best
    if likelihood == NOT_EXPECTED:
        return likelihood, f"{reason} (checked with {', '.join(tried)})"
    return likelihood, f"with {context}; rule score {score:g}; {reason}"


def _in_window(base, best, cfg):
    """True when a base score is within the tie-break window of the best one (the test rerank() and apply_ccs() share)."""
    return base >= (1 - cfg.tiebreak_window) * best - 1e-9


# ---------------------------------------------------------------- the factors
def _omics_factor(cls, chosen, cfg):
    if cls == NEITHER:
        return cfg.omics_neither
    detectable = chosen == OMICS_EITHER or cls == BOTH or (cls == LIPIDOME) == (chosen == OMICS_LIPIDOMICS)
    return cfg.omics_match if detectable else cfg.omics_mismatch


def _base_score(cand, cfg):
    """Best layer_score_norm + the bonus for every extra layer (capped), x the hypothesis factor."""
    sources = cand.sources or [{"layer": cand.layer, "layer_score_norm": cand.layer_score_norm}]
    best = max(float(s.get("layer_score_norm") or 0.0) for s in sources)
    base = min(cfg.base_cap, best + cfg.layer_bonus * (len({s["layer"] for s in sources}) - 1))
    return max(0.0, base * cfg.hypothesis_factor if cand.hypothesis else base)


def _detect(cand, *, sign, chosen, ph, adduct, additive, alkali, cfg):
    """The detect block of one candidate."""
    omics = ion = likelihood = None
    reasons = []
    used = cand.adduct or adduct
    if not has_structure(cand.smiles):
        why = "no valid structure, so detectability cannot be scored"
        reasons = [f"omics: {why}", f"ionization: {why}", f"adduct: {why}"]
    else:
        omics, verdict = _omics_check(cand.smiles)
        reasons.append(f"omics: {verdict} Chosen: {OMICS_LABELS[chosen]}." if omics else "omics: not evaluated")
        ion, ion_why = _ion_check(cand.smiles, sign, ph)
        reasons.append(f"ionization: ESI{sign} response {ion}; {ion_why}" if ion else f"ionization: not evaluated ({ion_why})")
        if used:
            likelihood, add_why = _best_adduct(cand.smiles, used, additive, True if alkali is None else bool(alkali))
            reasons.append(f"adduct: {used} {likelihood}; {add_why}" if likelihood else f"adduct: not evaluated ({add_why})")
        else:
            reasons.append("adduct: not evaluated (no adduct given)")
    factors = {
        "omics": _omics_factor(omics, chosen, cfg) if omics else cfg.not_evaluated,
        "ionizable": {ir.STRONG: cfg.ion_strong, ir.MODERATE: cfg.ion_moderate, ir.WEAK: cfg.ion_weak}.get(ion, cfg.not_evaluated),
        "adduct": {"high": cfg.adduct_high, "medium": cfg.adduct_medium, "low": cfg.adduct_low,
                   NOT_EXPECTED: cfg.adduct_not_expected}.get(likelihood, cfg.not_evaluated),
    }
    return {"omics": omics or NOT_EVALUATED, "ionizable": ion or NOT_EVALUATED, "adduct": likelihood or NOT_EVALUATED,
            "factors": factors, "factor": round(factors["omics"] * factors["ionizable"] * factors["adduct"], 6), "reasons": reasons}


# ---------------------------------------------------------------- inputs
def _sign(polarity):
    p = str(polarity or "").strip().lower()
    if p in ("positive", "pos", "+"):
        return "+"
    if p in ("negative", "neg", "-", ir.MINUS):
        return "-"
    raise ValueError(f"Polarity must be Positive or Negative, not {polarity!r}.")


def _chosen_omics(omics):
    o = str(omics or "").strip().lower()
    if o in ("", "either", "unknown", "either / unknown", "either/unknown"):
        return OMICS_EITHER
    if o.startswith("metab"):
        return OMICS_METABOLOMICS
    if o.startswith("lipid"):
        return OMICS_LIPIDOMICS
    raise ValueError(f"Omics must be Either / unknown, Metabolomics or Lipidomics, not {omics!r}.")


# ---------------------------------------------------------------- the re-rank
def rerank(candidates, *, polarity, adduct, omics, ph=DEFAULT_PH, additive=UNKNOWN, alkali=None, config=None, progress=None):
    """Score and order the merged candidates -> new Candidates (the input is not changed), best first, rank = 1..n.

    polarity   Positive | Negative: which ESI response counts for f_ionizable
    adduct     adduct of the run, used for a candidate that has no adduct of its own (a candidate's own adduct wins);
               None or Auto: such a candidate is not scored for the adduct
    omics      Either / unknown | Metabolomics | Lipidomics
    ph         mobile-phase pH for the ionization rules
    additive   "unknown" (the best additive counts) or one of adduct_rules.ADDITIVES; alkali: True / False, None = unknown = present
    config     RerankConfig (default values: the owner's start values)
    progress   optional callback progress(done, total), called after each structure is scored (the step takes seconds per structure)

    Fills base_score, final_score, rank and detect (see Candidate). Order: the window (base_score >= (1 - tiebreak_window) x the best) first,
    sorted by final_score, then base_score, then the input order; then everything else by base_score, then the input order. Outside the
    window final_score = base_score (module docstring). A window of one candidate (or a single candidate) changes nothing.
    """
    cfg = config or RerankConfig()
    sign, chosen, adduct = _sign(polarity), _chosen_omics(omics), canonical(adduct)
    additive = UNKNOWN if not additive or str(additive).strip().lower() == UNKNOWN else additive
    if not ir.PH_MIN <= ph <= ir.PH_MAX:
        raise ValueError(f"pH must be between {ir.PH_MIN:g} and {ir.PH_MAX:g}.")
    if additive != UNKNOWN and additive not in ar.ADDITIVES:
        raise ValueError(f"Additive must be unknown or one of {', '.join(ar.ADDITIVES)}, not {additive!r}.")
    if not 0 <= cfg.tiebreak_window < 1:
        raise ValueError("The tie-break window must be at least 0 and below 1.")

    candidates = list(candidates)
    scored = []                                             # [base, final, input index, candidate]
    for index, cand in enumerate(candidates):
        new = copy.deepcopy(cand)
        base = _base_score(new, cfg)
        new.detect = _detect(new, sign=sign, chosen=chosen, ph=ph, adduct=adduct, additive=additive, alkali=alkali, cfg=cfg)
        scored.append([base, base * new.detect["factor"], index, new])
        if progress:
            progress(index + 1, len(candidates))
    best = max((row[0] for row in scored), default=0.0)
    in_window = [_in_window(row[0], best, cfg) for row in scored]
    for row, inside in zip(scored, in_window):              # not a field, never serialized: apply_ccs() must see this decision, not one redone on the rounded base_score
        row[3]._in_tiebreak_window = inside
    window = sorted((row for row, inside in zip(scored, in_window) if inside), key=lambda r: (-round(r[1], 6), -round(r[0], 6), r[2]))
    rest = sorted((row for row, inside in zip(scored, in_window) if not inside), key=lambda r: (-round(r[0], 6), r[2]))
    for row in rest:                                        # the factors took no part in this position
        row[1] = row[0]
        row[3].detect["reasons"].append(f"re-rank: base score {row[0]:.2f} is more than {cfg.tiebreak_window:.0%} below the best "
                                        f"({best:.2f}), so detectability did not change its position")
    ordered = window + rest
    for rank, (base, final, _, new) in enumerate(ordered, start=1):
        new.base_score, new.final_score, new.rank = round(base, 4), round(final, 4), rank
    return [row[3] for row in ordered]


# ---------------------------------------------------------------- the CCS rule (step 6c)
CCS_INCONSISTENT = "Inconsistent"        # = mzdesk.orthogonal.ccs.INCONSISTENT (tests/test_annotation_ccs.py checks they agree)
CCS_MEASURED, CCS_PREDICTED = "measured", "predicted"
FLAG_CCS_MEASURED, FLAG_CCS_PREDICTED = "CCS ≠ measured", "CCS ≠ predicted"


def apply_ccs(ranked, config=None):
    """Move candidates whose CCS disagrees with their reference -> new Candidates (the input is not changed), rank renumbered 1..n.

    `ranked` is the order rerank() returned (the tie-break window first) and every candidate carries `ccs` (ccs_step.check). The rule:

      measured, strong   A candidate whose CCS is Inconsistent with a measured reference goes behind every candidate that is not
                         measured-inconsistent, wherever it stood (also outside the window). Those candidates keep their relative order.
      predicted, weak    A candidate whose CCS is Inconsistent with a predicted value only moves inside the tie-break window: behind the
                         window's other members. A predicted-inconsistent candidate outside the window stays where it is, and one that
                         is also behind a measured-inconsistent one is simply behind it.
      everything else    Consistent, No reference, No structure, Not checked (and a candidate without a `ccs`): not moved by the rule.
                         Candidates only change place by the others leaving theirs; relative order is otherwise stable.

    The window is the re-rank's own (base_score >= (1 - tiebreak_window) x the best base_score), measured on the list as given, so a demoted
    candidate still counts for it. rerank() leaves its decision on each candidate (it was made on the exact base score, which base_score only holds rounded
    to 4 decimals); a candidate that did not come from rerank() is judged by its base_score. Scores, levels and `detect` are not touched: CCS is supporting evidence, not a score. A candidate the rule
    applies to gets ccs["flag"] ("CCS ≠ measured" / "CCS ≠ predicted"); one that a candidate behind it overtook also gets ccs["moved_from"], its rank before
    (a flagged candidate can still gain places when a measured-inconsistent one ahead of it leaves; it then has a flag and no moved_from).
    """
    cfg = config or RerankConfig()
    cands = [copy.deepcopy(c) for c in ranked]
    best = max((c.base_score or 0.0 for c in cands), default=0.0)
    window = {i for i, c in enumerate(cands)
              if (c._in_tiebreak_window if hasattr(c, "_in_tiebreak_window") else _in_window(c.base_score or 0.0, best, cfg))}     # rerank()'s own decision; recomputed only for hand-made candidates

    def inconsistent_with(i, source):
        x = cands[i].ccs or {}
        return x.get("verdict") == CCS_INCONSISTENT and x.get("source") == source

    strong = [i for i in range(len(cands)) if inconsistent_with(i, CCS_MEASURED)]
    weak = [i for i in sorted(window) if inconsistent_with(i, CCS_PREDICTED) and i not in strong]
    order = [i for i in range(len(cands)) if i not in strong and i not in weak]
    if weak:
        at = max((k for k, i in enumerate(order) if i in window), default=-1) + 1          # right behind the window's other members
        order[at:at] = weak
    order += strong

    place = {i: new_pos for new_pos, i in enumerate(order)}
    for i, c in enumerate(cands):
        c.rank = place[i] + 1
        if i in strong or i in weak:
            passed = any(place[j] < place[i] for j in range(i + 1, len(cands)))          # a candidate that stood behind it is now ahead of it
            c.ccs = {**c.ccs, "flag": FLAG_CCS_MEASURED if i in strong else FLAG_CCS_PREDICTED, "moved_from": i + 1 if passed else None}
    return [cands[i] for i in order]
