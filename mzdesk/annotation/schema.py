"""Candidate and LayerResult: the data contract of the ESI MS/MS Annotation workflow.

Who fills what
  steps 3-5  every layer (library, MS-FINDER, FRIGID) returns a LayerResult holding its own Candidates. A layer fills the identity
             fields, mz_error_ppm, layer_score, layer_score_norm, evidence and hypothesis, and nothing else.
  step 6a    merge.merge_candidates() makes one Candidate per structure (ik14) and fills id, sources, flags, level and level_reason.
  step 6b    rerank.rerank() fills detect, base_score, final_score and rank.
  step 6c    ccs_step.check() fills ccs when the form has an observed CCS; rerank.apply_ccs() may then move the candidate and renumber rank.
  step 7     answer.attach_isomer() fills isomer.
  step 8     Candidate.to_dict() is one row of the ranked list; answer.build_answer() writes the answer block.

How a layer fills layer_score_norm
  1. Absolute: it depends only on the candidate's own raw score and the engine's documented score range, never on the other
     candidates of the run. (Dividing by the best score of the list would make the winner of every list 1.0 and break both the
     re-rank base and the "good enough" tests.)
  2. Monotone increasing in the raw score and clipped to [0, 1]; 1 is the strongest evidence this layer can give, 0 is none.
  3. Library: the similarity of the primary metric (already 0-1). MS-FINDER and FRIGID: a documented monotone map of the engine
     score onto [0, 1], calibrated on real score ranges and written next to the step's code.
  4. Norms of different layers are not comparable one to one: merge keeps the best norm per layer, rerank takes the max over layers.

How a layer fills evidence
  A dict keyed by the layer name, so merged candidates keep every layer's evidence side by side:
    evidence["library"]  = {"gate_passed": bool (required: it decides Level 2a against 3), "failed_checks": [str],
                            "similarity": float, "metric": str, "matched_peaks": int, "precursor_error_ppm": float, ...}
    evidence["msfinder"] = {"total_score": float, ...}
    evidence["frigid"]   = {...}
  Everything except gate_passed is for display only.
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field

from mzdesk.adducts import canonical as canonical_adduct

LAYER_LIBRARY, LAYER_MSFINDER, LAYER_FRIGID = "library", "msfinder", "frigid"
LAYER_FORMULA = "formula_only"          # a row that only carries a molecular formula (no structure was proposed)
# The structural layers, strongest evidence first (an experimental spectrum, then in-silico fragmentation, then a hypothesis).
LAYERS = (LAYER_LIBRARY, LAYER_MSFINDER, LAYER_FRIGID)
LAYER_LABELS = {LAYER_LIBRARY: "library", LAYER_MSFINDER: "MS-FINDER", LAYER_FRIGID: "FRIGID", LAYER_FORMULA: "formula-only"}

STATUS_DONE, STATUS_SKIPPED, STATUS_FAILED = "done", "skipped", "failed"
STATUSES = (STATUS_DONE, STATUS_SKIPPED, STATUS_FAILED)

FLAG_NO_STRUCTURE = "no_structure"      # no valid SMILES: kept in the list, but detectability cannot be scored


@dataclass
class Candidate:
    """One row of the ranked list. JSON-able through to_dict().

    Filled by the layer that proposes it
      id               unique within one result (a layer may use any string; merge renumbers: "c1", "c2", ...)
      layer            library | msfinder | frigid, or formula_only for a formula row. After merge: the contributing layer with
                       the strongest evidence (library > MS-FINDER > FRIGID); the others are in `sources`
      smiles           structure of the neutral molecule M (None when the layer has none); the first valid one wins in a merge
      inchikey         full InChIKey (merge computes it from a valid SMILES when the layer gives none)
      ik14             InChIKey first block = connectivity only (no stereo, no isotopes): the key candidates are merged by
      name             compound name or "" (a FRIGID structure has none)
      formula          molecular formula of M
      adduct           canonical adduct name as in the app's adduct lists: [M+FA-H]- becomes [M+HCOO]- when the Candidate is
                       made (mzdesk.adducts.canonical)
      mz_error_ppm     (observed - theoretical) / theoretical x 1e6: precursor m/z against this formula and adduct
      layer_score      the layer's raw score (library similarity, MS-FINDER total score, ...); None when it has none
      layer_score_norm the same score as 0-1, see "How a layer fills layer_score_norm" above
      evidence         dict keyed by layer, see "How a layer fills evidence" above
      hypothesis       True for a structure that rests on formula-level evidence only (FRIGID); after a merge: True only when
                       every contributing layer is such a hypothesis

    Filled by merge (step 6a)
      sources          one entry per contributing layer, strongest first:
                       {"layer", "layer_score", "layer_score_norm", "rank_in_layer"} (rank_in_layer: 1-based position in
                       the layer's own list, so a layer must return its candidates best first)
      flags            short codes; FLAG_NO_STRUCTURE = no valid SMILES (kept, never scored by detectability)
      level            "2a" | "3" | "H" | "4" | "5" (levels.py)
      level_reason     one sentence why

    Filled by rerank (step 6b)
      detect           {"omics": lipidome | metabolome | both | neither | not evaluated,
                        "ionizable": strong | moderate | weak | not evaluated   (ESI response in the precursor polarity),
                        "adduct": high | medium | low | not expected | not evaluated   (likelihood of the observed adduct),
                        "factors": {"omics", "ionizable", "adduct"} the three multipliers,
                        "factor": their product,
                        "reasons": [short English strings, one per dimension]}
      base_score       max layer_score_norm over the layers + a bonus per extra layer, x the hypothesis factor
      final_score      base_score x detect["factor"]; detectability only orders candidates whose base_score is within 10% of the best, the rest keep layer-score order
      rank             1-based position in the ranked list

    Filled by the CCS check (step 6c), only when the form has an observed CCS (otherwise None, and to_dict() leaves the key out)
      ccs              {"verdict": Consistent | Inconsistent | No reference | No structure | Not checked,
                        "source": "measured" | "predicted" | None   (what the verdict is against),
                        "adduct": the adduct that was looked up, "observed": the observed CCS (Å²),
                        "reference": the reference or predicted CCS (Å²), "delta_pct": (observed - reference) / reference x 100,
                        "tolerance_pct": ±%, "ref_label": CCS type of the measured row ("DTCCS N2") or the model name,
                        "dataset": source dataset of the measured row, "n_refs": measured values in the index,
                        "in_domain": prediction only, False = outside the model's training domain, "note": why there is no verdict,
                        "flag": "CCS ≠ measured" | "CCS ≠ predicted" | None   (the rule of rerank.apply_ccs applied to this candidate),
                        "moved_from": the rank before the rule moved it, or None}
                       CCS is supporting evidence only: it never changes `level`.

    Filled by the isomer check (step 7), answer.attach_isomer()
      isomer           {"verdict": single | indistinguishable | insufficient | not run (the same on every row),
                        "role": "top" (the top-1), "separated" (told apart from the top-1), "not separated" (cannot be told
                                from it), or None (not compared),
                        "challenger": True for an isomer that explains observed peaks the top-1 cannot (information only),
                        "explained_fraction": share of the considered intensity this structure explains, None when unknown}
    """
    id: str
    layer: str
    smiles: str | None = None
    inchikey: str | None = None
    ik14: str | None = None
    name: str = ""
    formula: str | None = None
    adduct: str | None = None
    mz_error_ppm: float | None = None
    layer_score: float | None = None
    layer_score_norm: float = 0.0
    evidence: dict = field(default_factory=dict)
    hypothesis: bool = False
    sources: list = field(default_factory=list)
    flags: list = field(default_factory=list)
    detect: dict | None = None
    base_score: float | None = None
    final_score: float | None = None
    rank: int | None = None
    level: str | None = None
    level_reason: str = ""
    isomer: dict | None = None
    ccs: dict | None = None

    def __post_init__(self):
        if self.layer not in LAYERS + (LAYER_FORMULA,):
            raise ValueError(f"Unknown layer {self.layer!r}; choose one of {', '.join(LAYERS + (LAYER_FORMULA,))}.")
        self.adduct = canonical_adduct(self.adduct) or self.adduct
        if self.layer_score_norm is None:
            self.layer_score_norm = 0.0

    def to_dict(self):
        """JSON-able dict: every field (ccs only when the CCS check ran), plus the display names layer_label and level_label. NaN / inf become None."""
        from mzdesk.annotation.levels import level_label      # levels imports this module, so import it here

        d = _jsonable(asdict(self))
        if d.get("ccs") is None:                    # no observed CCS in the form: the row is exactly what it was before the CCS check existed
            d.pop("ccs", None)
        d["layer_label"] = LAYER_LABELS[self.layer]
        d["level_label"] = level_label(self.level)
        return d


@dataclass
class LayerResult:
    """What a layer returns: its status, a one-line note for the step table and its candidates (best first).

    layer       library | msfinder | frigid
    status      done | skipped | failed ("skipped": the layer did not apply or was not needed; the note says why)
    note        one plain sentence for the step table (why skipped or failed, how many candidates, ...)
    seconds     wall time of the layer
    candidates  at most one per compound (merge keeps the best one if a layer repeats a compound), best first
    """
    layer: str
    status: str = STATUS_DONE
    note: str = ""
    seconds: float = 0.0
    candidates: list = field(default_factory=list)

    def __post_init__(self):
        if self.layer not in LAYERS:
            raise ValueError(f"Unknown layer {self.layer!r}; choose one of {', '.join(LAYERS)}.")
        if self.status not in STATUSES:
            raise ValueError(f"Unknown status {self.status!r}; choose one of {', '.join(STATUSES)}.")

    def to_dict(self):
        return {"layer": self.layer, "status": self.status, "note": self.note, "seconds": _jsonable(self.seconds),
                "candidates": [c.to_dict() for c in self.candidates]}


def evidence_of(candidate, layer):
    """The evidence dict of one layer: evidence[layer] when the layer keyed it, else the whole dict (a layer that forgot to key it)."""
    ev = candidate.evidence or {}
    inner = ev.get(layer)
    return inner if isinstance(inner, dict) else ev


def _jsonable(value):
    """Plain JSON types: NaN / inf -> None, numpy scalars -> Python numbers, tuples -> lists."""
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if hasattr(value, "item") and not isinstance(value, (str, bytes)):         # numpy scalar
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value
