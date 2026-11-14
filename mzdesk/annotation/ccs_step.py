"""Step 6c: check the ranked candidates against an observed CCS (optional; supporting evidence only).

The form may carry one observed CCS (Å², one value) with the platform and drift gas it was measured on (the "Ion Mobility" section). Then each
ranked candidate that has a structure is looked up in the local CCS index with the CCS page's own engine (mzdesk/orthogonal/ccs.py, whose
comparison rules are the owner's and are not repeated here): the same compound (InChIKey-14) and the same adduct, N2 scale only, a measured
reference first and a model prediction only when there is none. The candidate's own adduct is used (the run's adduct when it has none).

What comes out is one `ccs` dict per candidate (schema.Candidate.ccs). The ranking effect is rerank.apply_ccs(); this module only measures:

  - Not comparable (He, FAIMS / DMS, TWIMS calibrated with He): no candidate is looked up, nothing moves, one note says why.
  - The CCS index is not installed (or has an old format): same, and the note carries the setup command.
  - A candidate without a structure: "No structure". One the engine cannot read: "Not checked" with the engine's message in `note`.
  - Everything else: Consistent / Inconsistent / No reference, with the numbers behind it.

CCS never changes a candidate's level.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field

from mzdesk.annotation.merge import has_structure
from mzdesk.api.common import InputError, number, text
from mzdesk.orthogonal import ccs as engine

log = logging.getLogger("mzdesk.annotation")

NO_STRUCTURE, NOT_CHECKED = "No structure", "Not checked"
FORM_KEYS = ("ccs_observed", "ccs_platform", "ccs_gas", "ccs_tw_calibration")        # the "Ion Mobility" fields of the spec


@dataclass(frozen=True)
class CcsQuery:
    """The Ion Mobility section of the form, parsed: one observed CCS and how it was measured."""
    observed: float
    platform: str
    gas: str = engine.GASES[0]
    tw_calibration: str = engine.TW_CALIBRATIONS[0]

    @property
    def label(self):
        """"152.3 Å² (TIMS, N₂)": what the notes say about the value."""
        return f"{self.observed:g} Å² ({self.platform}, {self.gas})"


@dataclass
class Outcome:
    """What check() found. checks is None when the CCS could not be used at all; message then says why (full sentence) and short why in a few words."""
    checks: dict | None = None            # candidate id -> the ccs dict
    message: str = ""
    short: str = ""
    counts: dict = field(default_factory=dict)


# ---------------------------------------------------------------- the form
def parse_form(inputs, warnings):
    """Form JSON -> CcsQuery, or None when no observed CCS was given (then every other Ion Mobility field is ignored). Raises InputError."""
    observed = number(inputs, "ccs_observed", "Observed CCS (Å²)")
    if observed is None:
        return None
    try:
        engine.parse_observed(observed)                       # the CCS page's limits: 1/K0-like values (too small) and non-molecular values are refused
    except ValueError as exc:
        raise InputError(str(exc)) from None
    platform = text(inputs, "ccs_platform")
    if not platform:
        raise InputError("Platform is required when an observed CCS is given (Ion Mobility section): the comparison depends on how the CCS was measured.")
    gas = text(inputs, "ccs_gas") or engine.GASES[0]
    tw = text(inputs, "ccs_tw_calibration") or engine.TW_CALIBRATIONS[0]
    for value, options, label in ((platform, tuple(engine.PLATFORMS), "Platform"), (gas, engine.GASES, "Drift gas"),
                                  (tw, engine.TW_CALIBRATIONS, "TWIMS calibration")):
        if value not in options:
            raise InputError(f"{label} must be one of: {', '.join(options)}.")
    if platform != "TWIMS" and tw != engine.TW_CALIBRATIONS[0]:
        warnings.append("The TWIMS calibration choice only applies to TWIMS and was ignored.")
        tw = engine.TW_CALIBRATIONS[0]
    return CcsQuery(observed=observed, platform=platform, gas=gas, tw_calibration=tw)


# ---------------------------------------------------------------- the check
def _from_result(res):
    """engine.Result -> the ccs dict of a candidate (verdict and the numbers behind it)."""
    d = {"verdict": res.verdict, "source": res.source or None, "adduct": res.adduct, "observed": res.observed[0] if res.observed else None,
         "reference": None, "delta_pct": None, "tolerance_pct": None, "ref_label": None, "dataset": None, "n_refs": len(res.rows), "in_domain": None,
         "note": "", "flag": None, "moved_from": None}
    used = res.used
    if used is not None and res.verdict in (engine.CONSISTENT, engine.INCONSISTENT):
        d.update(reference=round(used.ccs, 2), delta_pct=round(used.delta_pct, 3), tolerance_pct=used.tolerance)
        if res.source == "measured":
            d.update(ref_label=engine.TYPE_LABEL[used.ref_type], dataset=used.row.src_tag)
        elif res.prediction is not None:
            d.update(ref_label=res.prediction.model, in_domain=res.prediction.in_domain)
    elif res.verdict == engine.NO_REFERENCE:
        d["note"] = res.prediction_note or "none in the index"
    return d


def _check_one(cand, q, run_adduct):
    if not has_structure(cand.smiles):
        return {"verdict": NO_STRUCTURE, "note": "a CCS can only be looked up for a structure"}
    adduct = cand.adduct or run_adduct
    if not adduct:
        return {"verdict": NOT_CHECKED, "note": "no adduct to look up"}
    try:
        res = engine.analyse(cand.smiles, adduct, [q.observed], q.platform, q.gas, q.tw_calibration)
    except engine.QueryError as exc:
        return {"verdict": NOT_CHECKED, "adduct": adduct, "note": str(exc)}
    except engine.IndexProblem:
        raise
    except Exception as exc:        # noqa: BLE001 - one odd structure must not take the other candidates' CCS checks down
        log.exception("CCS check of %s failed", cand.id)
        return {"verdict": NOT_CHECKED, "adduct": adduct, "note": f"{type(exc).__name__}: {exc}"}
    return _from_result(res)


def check(candidates, q, run_adduct, progress=None):
    """Outcome: the ccs dict of every candidate (by id), or no checks at all with the reason (not comparable, index missing).

    candidates  the ranked Candidates (not changed); run_adduct  adduct of a candidate that has none; progress(done, total) after each candidate.
    """
    comparable, short, reason = engine.comparability(q.platform, q.gas, q.tw_calibration)
    if not comparable:
        return Outcome(message=f"The observed CCS ({q.label}) was not used: {reason}", short=short)
    if not engine.installed():
        return Outcome(message=f"The observed CCS ({q.label}) was not used: the CCS reference data is not installed. Run `{engine.SETUP_COMMAND}`, "
                               "then run the workflow again.", short="CCS data not installed")
    try:
        engine.load()
        checks = {}
        for done, cand in enumerate(candidates, start=1):
            checks[cand.id] = _check_one(cand, q, run_adduct)
            if progress:
                progress(done, len(candidates))
    except engine.IndexProblem as exc:
        return Outcome(message=f"The observed CCS ({q.label}) was not used: {exc}", short="CCS index needs a rebuild")
    counts = {}
    for c in checks.values():
        key = c["verdict"] if c["verdict"] != engine.INCONSISTENT else f"{c['verdict']} {c['source']}"
        counts[key] = counts.get(key, 0) + 1
    return Outcome(checks=checks, counts=counts)


# ---------------------------------------------------------------- words
def summary(q, outcome, moved):
    """One sentence for the step note: the value, how many candidates agree / disagree / have no reference, and how many CCS moved."""
    n = outcome.counts
    parts = []
    for key, text_ in ((engine.CONSISTENT, "consistent"), (f"{engine.INCONSISTENT} measured", "inconsistent with a measured value"),
                       (f"{engine.INCONSISTENT} predicted", "inconsistent with a predicted value"), (engine.NO_REFERENCE, "without a reference")):
        if n.get(key):
            parts.append(f"{n[key]} {text_}")
    other = n.get(NO_STRUCTURE, 0) + n.get(NOT_CHECKED, 0)
    if other:
        parts.append(f"{other} not checked")
    sentence = f"Observed CCS {q.label}: {', '.join(parts) or 'no candidate checked'}."
    return sentence + (f" {moved} moved by CCS." if moved else "")


def describe(c):
    """The disagreement of one candidate in a phrase: "+9.8 % from the measured 140.5 Å² (DTCCS N2, ±1.0 %)"."""
    x = c.ccs or {}
    if x.get("delta_pct") is None:
        return ""
    kind = "measured" if x.get("source") == "measured" else "predicted"
    sign = "+" if x["delta_pct"] >= 0 else "−"
    return f"{sign}{abs(x['delta_pct']):.1f} % from the {kind} {x['reference']:.1f} Å² (±{x['tolerance_pct']:.1f} %)"
