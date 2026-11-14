"""Ranking and levels of the EI MS Annotation workflow (step 7 and the answer). Pure functions on ei_candidates.Candidate; no engine at import.

Order (rerank): evidence tier first, then consistency with the typed molecular ion, then the score.
  tier        a candidate with an experimental EI library spectrum match ("library") ranks above one with only a NEIMS-predicted spectrum ("neims")
  M           with a molecular ion given, a candidate whose derivative nominal mass is not M (and is not a partly silylated form of it) goes after the
              consistent ones of its tier and carries the flag "M mismatch"; one without a structure (nothing to compare) sits between the two
  score       library tier: the library similarity, the NEIMS similarity breaks ties; NEIMS tier: the NEIMS similarity
Both similarities are the primary metric of the library search (weighted spectral entropy, unit mass, 0.5 Da).

Levels (ei_level) use the ESI workflow's scheme and its Level-2 gate numbers (annotation/levels.py, annotation/library_step.py); there is no
Level 1 (no reference standard, the library has no retention index):
  Level 2a  library tier, similarity >= the gate (0.80), at least min(4, query peaks) matched peaks, the entry has a structure, the entry's
            derivatization state is that of the query, and M does not contradict the derivative
  Level 3   any other library-tier candidate, and a NEIMS-tier candidate whose predicted spectrum scores >= NEIMS_MATCH_MIN with an M that does
            not contradict it
  none      everything else (no identification)
NEIMS_MATCH_MIN is a working value (not a published cut-off): with the 0.80 gate for measured spectra a predicted spectrum is held to a lower bar.
"""
from __future__ import annotations

NEIMS_MATCH_MIN = 0.40
TIER_LIBRARY, TIER_NEIMS = "library", "neims"
M_OK, M_PARTIAL, M_MISMATCH, M_UNKNOWN = "ok", "partial", "mismatch", "unknown"       # unknown: M is given but the candidate has no derivative structure to compare
SILYL_NOMINAL = {"TMS": 72, "MeOx + TMS": 72, "TBDMS": 114}
M_TOLERANCE = 0.5               # Da between a typed high-resolution M and the derivative's monoisotopic mass


def gate():
    """-> (similarity gate, minimum matched peaks, list threshold): the ESI workflow's Level-2 numbers."""
    from mzdesk.annotation import library_step as ls
    return ls.GATE_SIMILARITY, ls.GATE_MIN_MATCHED, ls.LIST_THRESHOLD


def m_status(variant, molecular_ion):
    """Does the derivative ``variant`` (dict: mode, nominal, mass, n_silyl) fit the typed molecular ion? -> (status, k): M_OK | M_PARTIAL (M is the
    derivative with only k of its n_silyl silyl groups) | M_MISMATCH, or (None, None) without M or without a derivative structure."""
    if molecular_ion is None or not variant or variant.get("nominal") is None:
        return None, None
    nominal, n = variant["nominal"], variant.get("n_silyl") or 0
    target = int(molecular_ion + 0.5)
    if nominal == target or abs(variant.get("mass", nominal) - molecular_ion) <= M_TOLERANCE:
        return M_OK, n
    step = SILYL_NOMINAL.get(variant.get("mode"), 0)
    for k in range(n - 1, 0, -1):
        if step and nominal - (n - k) * step == target:
            return M_PARTIAL, k
    return M_MISMATCH, None


def sort_key(c):
    """Smaller sorts first: (tier, M fit, -score, -NEIMS score, name). A candidate without a score sorts last in its group."""
    library = c.tier == TIER_LIBRARY
    primary = c.lib_score if library else c.neims_score
    secondary = (c.neims_score or 0.0) if library else 0.0
    m_rank = {M_MISMATCH: 2, M_UNKNOWN: 1}.get(c.m_status, 0)
    return (0 if library else 1, m_rank, -(primary if primary is not None else -1.0), -secondary, c.name.lower())


def rerank(candidates):
    """Candidates in rank order, ``rank`` set (1 = best). The input list is not changed."""
    ordered = sorted(candidates, key=sort_key)
    for i, c in enumerate(ordered, 1):
        c.rank = i
    return ordered


def state_matches(selected_mode, reference_state):
    """True when a reference spectrum shows the derivatization family of the query (selected_mode: the effective mode, reference_state: "TMS" |
    "TBDMS" | "None" or None when unknown, which counts as matching)."""
    from mzdesk.ms2 import ei_markers as em
    return reference_state is None or em.family(selected_mode) == reference_state


def library_problems(c, n_query_peaks):
    """Short phrases for why a library-tier candidate is not Level 2a ([] when it is)."""
    sim_gate, min_matched, _ = gate()
    lib, problems = c.library, []
    if lib["score"] < sim_gate:
        problems.append(f"score below {sim_gate:.2f}")
    needed = min(min_matched, n_query_peaks)
    if lib["matches"] < needed:
        problems.append(f"only {lib['matches']} matched peaks (needs {needed})")
    if c.parent_smiles is None:
        problems.append("the library entry has no structure")
    if not c.state_ok:
        problems.append("the reference is in another derivatization state")
    if c.m_status == M_MISMATCH:
        problems.append("the derivative mass does not fit the molecular ion")
    return problems


def ei_level(c, n_query_peaks):
    """-> (level "2a" | "3" | None, one-sentence reason) for a candidate whose tier, scores, M status and state flag are set (see the module docstring)."""
    from mzdesk.annotation.levels import LEVEL_2A, LEVEL_3
    if c.tier == TIER_LIBRARY and c.library:
        lib = c.library
        problems = library_problems(c, n_query_peaks)
        detail = f"similarity {lib['score']:.2f}, {lib['matches']} matched peaks"
        if not problems:
            return LEVEL_2A, f"EI library spectrum match ({detail})."
        return LEVEL_3, f"EI library spectrum match ({detail}), but {'; '.join(problems)}."
    if c.neims_score is not None and c.neims_score >= NEIMS_MATCH_MIN and c.m_status != M_MISMATCH:
        return LEVEL_3, f"NEIMS-predicted spectrum match only (similarity {c.neims_score:.2f}); no library spectrum."
    if c.neims_score is None:
        return None, "No library spectrum and no predicted spectrum to compare."
    why = "the derivative mass does not fit the molecular ion" if c.m_status == M_MISMATCH else f"similarity {c.neims_score:.2f} is below {NEIMS_MATCH_MIN:.2f}"
    return None, f"The NEIMS-predicted spectrum does not match: {why}."
