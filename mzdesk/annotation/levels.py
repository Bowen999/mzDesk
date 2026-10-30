"""Confidence levels of an annotation (Schymanski et al., Environ. Sci. Technol. 2014) as the workflow uses them.

  2a  a library spectrum passed the Level-2 gate (step 3)
  3   an MS-FINDER candidate, or a library hit that failed the gate
  H   a FRIGID structure hypothesis: formula-level evidence only, so Level 4 at most
  4   nothing structural, but a molecular formula
  5   only an exact mass

A merged candidate gets the best level of its contributing layers. The isomer verdict is reported next to the level and never
changes it, and neither does re-ranking.
"""
from mzdesk.annotation.schema import (LAYER_FORMULA, LAYER_FRIGID, LAYER_LABELS, LAYER_LIBRARY, LAYER_MSFINDER, evidence_of)

LEVEL_2A, LEVEL_3, LEVEL_H, LEVEL_4, LEVEL_5 = "2a", "3", "H", "4", "5"
LEVELS = (LEVEL_2A, LEVEL_3, LEVEL_H, LEVEL_4, LEVEL_5)       # best evidence first
LEVEL_NAMES = {LEVEL_2A: "Level 2a", LEVEL_3: "Level 3", LEVEL_H: "Hypothesis", LEVEL_4: "Level 4", LEVEL_5: "Level 5"}
LEVEL_MEANINGS = {
    LEVEL_2A: "probable structure (library spectrum match)",
    LEVEL_3: "tentative candidate",
    LEVEL_H: "formula-level evidence only (Level 4 at most)",
    LEVEL_4: "molecular formula only",
    LEVEL_5: "exact mass only",
}


def level_name(level):
    """Short display name of a level ("" for None): Level 2a, Level 3, Hypothesis, Level 4, Level 5."""
    return LEVEL_NAMES.get(level, "")


def level_label(level):
    """Display name with its meaning ("" for None), e.g. "Level 3: tentative candidate"."""
    return f"{LEVEL_NAMES[level]}: {LEVEL_MEANINGS[level]}" if level in LEVEL_NAMES else ""


def contribution_level(layer, evidence):
    """Level that one layer's finding earns. `evidence` is that layer's evidence dict (library: gate_passed decides 2a against 3)."""
    if layer == LAYER_LIBRARY:
        return LEVEL_2A if evidence.get("gate_passed") else LEVEL_3
    if layer == LAYER_MSFINDER:
        return LEVEL_3
    if layer == LAYER_FRIGID:
        return LEVEL_H
    if layer == LAYER_FORMULA:
        return LEVEL_4
    raise ValueError(f"Unknown layer {layer!r}.")


def fallback_level(has_formulas, has_mass):
    """Level when no structure was proposed: a formula list gives 4, only an exact mass gives 5, nothing gives None."""
    return LEVEL_4 if has_formulas else LEVEL_5 if has_mass else None


def _reason(level, layer, ev, others):
    """One sentence for the winning contribution; `others` are the other contributing layers."""
    if layer == LAYER_LIBRARY and level == LEVEL_2A:
        detail = [f"similarity {ev['similarity']:.2f}"] if isinstance(ev.get("similarity"), (int, float)) else []
        if ev.get("matched_peaks") is not None:
            detail.append(f"{ev['matched_peaks']} matched peaks")
        text = "Library spectrum passed the Level-2 gate" + (f" ({', '.join(detail)})" if detail else "")
    elif layer == LAYER_LIBRARY:
        failed = ev.get("failed_checks") or []
        text = "Library hit that did not pass the Level-2 gate" + (f" (failed: {'; '.join(failed)})" if failed else "")
    elif layer == LAYER_MSFINDER:
        text = "MS-FINDER structure candidate: in-silico fragment fit, not confirmed by a library spectrum"
    elif layer == LAYER_FRIGID:
        text = "FRIGID hypothesis built from the formula, with no spectral evidence for this structure"
    else:
        text = "Molecular formula only, no structure was proposed"
    if others:
        text += "; also proposed by " + ", ".join(LAYER_LABELS[o] for o in others)
    return text + "."


def assign_level(candidate):
    """Set candidate.level and candidate.level_reason from the layers that contribute to it (candidate.sources, else candidate.layer).

    Returns the level. The best level wins; on a tie the layer with the stronger evidence (library, MS-FINDER, FRIGID) is named.
    """
    layers = [s["layer"] for s in candidate.sources] or [candidate.layer]
    earned = [(contribution_level(layer, evidence_of(candidate, layer)), layer) for layer in layers]
    level, layer = min(earned, key=lambda e: LEVELS.index(e[0]))
    candidate.level = level
    candidate.level_reason = _reason(level, layer, evidence_of(candidate, layer), [x for x in layers if x != layer])
    return level
