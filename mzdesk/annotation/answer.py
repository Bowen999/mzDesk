"""Steps 7-8: the answer block of the result, and the isomer annotation of the ranked rows.

The isomer result is the dict of mzdesk/exclusivity/isomers.discriminate(); only these keys are used:
  verdict      single | indistinguishable | insufficient
  top_id       id of the top-1 candidate the check was run for
  group        ids that cannot be separated from the top-1, the top-1 included (verdict indistinguishable)
  isomers      [{"id", "distinguished": True | False | None}]: the top-1 against each isomer (None: the engine did not score it)
  challengers  isomers that explain observed peaks the top-1 cannot, information only (dicts with an "id", or ids)
  candidates   [{"id", "explained_fraction"}]: the share of the considered intensity each structure explains
  notes        plain sentences; the first one is quoted when the check could not be done

The answer is always relative to the candidates that were considered; every statement about isomers says so.
"""
from __future__ import annotations

from mzdesk.annotation.levels import fallback_level, level_label, level_name
from mzdesk.annotation.merge import formula_only_candidates
from mzdesk.annotation.schema import LAYER_FORMULA, LAYER_LABELS

ANSWER_SINGLE, ANSWER_GROUP, ANSWER_FORMULA, ANSWER_MASS, ANSWER_NONE = "single", "group", "formula_only", "mass_only", "none"
VERDICT_SINGLE, VERDICT_GROUP, VERDICT_NOT_RUN = "single", "indistinguishable", "not run"      # insufficient: anything else
ROLE_TOP, ROLE_SEPARATED, ROLE_NOT_SEPARATED = "top", "separated", "not separated"


def _ids(items):
    """Ids from a list of ids or of dicts with an "id"."""
    return [x if isinstance(x, str) else x["id"] for x in items or []]


def _label(cand):
    """How the statement names a candidate: its name, else its SMILES, else its id."""
    return cand.name or cand.smiles or cand.id


def _join(labels):
    return labels[0] if len(labels) == 1 else ", ".join(labels[:-1]) + " and " + labels[-1]


def _plural(n, word):
    return f"{n} {word}" + ("" if n == 1 else "s")


# ---------------------------------------------------------------- the isomer annotation of the rows
def attach_isomer(ranked, isomer_result):
    """Fill Candidate.isomer on every row, in place, from the isomer check; returns `ranked`. Without a result the rows stay as they are.

    isomer = {"verdict": the run's verdict (the same on every row),
              "role": "top" for the top-1; "separated" for an isomer the top-1 is told apart from; "not separated" for one it is not
                      (these are the group); None for a row that was not compared (another formula, or no check),
              "challenger": True for an isomer that explains peaks the top-1 cannot,
              "explained_fraction": share of the considered intensity this structure explains, None when unknown}
    """
    if not isomer_result:
        return ranked
    top_id = isomer_result.get("top_id")
    group = set(_ids(isomer_result.get("group")))
    separated = {i["id"] for i in isomer_result.get("isomers") or [] if i.get("distinguished") is True}
    challengers = set(_ids(isomer_result.get("challengers")))
    fractions = {c["id"]: c.get("explained_fraction") for c in isomer_result.get("candidates") or []}
    for cand in ranked:
        role = ROLE_TOP if cand.id == top_id else ROLE_NOT_SEPARATED if cand.id in group else ROLE_SEPARATED if cand.id in separated else None
        cand.isomer = {"verdict": isomer_result.get("verdict") or VERDICT_NOT_RUN, "role": role, "challenger": cand.id in challengers,
                       "explained_fraction": fractions.get(cand.id)}
    return ranked


# ---------------------------------------------------------------- the answer block
def build_answer(ranked, isomer_result=None, *, formula_list=None, precursor_mz=None):
    """The `answer` block of the result.

    ranked         the re-ranked structural candidates (rerank.rerank); formula-only rows in it are ignored
    isomer_result  the dict of the isomer check, None when it was not run
    formula_list   step 2's shared formula list (strings, or dicts with "formula", "adduct", ...): the answer when no structure exists
    precursor_mz   the answer (Level 5) when there is no structure and no formula

    Returns {"type": single | group | formula_only | mass_only | none,
             "ids": row ids the answer points at (single: the top-1; group: the group, top-1 first; formula_only: the ids of
                    formula_only_candidates(formula_list)),
             "layer", "layer_label", "level", "level_label": those of the top-1 (the isomer verdict never changes the level),
                    or Level 4 / 5 for a formula / mass answer, "" for none,
             "statement": plain English, with the isomer caveat whenever isomers are mentioned,
             "top_id": id of the top-1 or None, "isomer_verdict": single | indistinguishable | insufficient | not run, or None}
    """
    structural = [c for c in ranked if c.layer != LAYER_FORMULA]
    if structural:
        return _structure_answer(structural, isomer_result)
    rows = formula_only_candidates(formula_list)
    level = fallback_level(bool(rows), precursor_mz is not None)
    if rows:
        shown = [f"{r.formula} ({r.adduct})" if r.adduct else r.formula for r in rows]
        statement = f"No structure could be proposed. The data support {_plural(len(rows), 'molecular formula')}: {_join(shown)}."
        return _answer(ANSWER_FORMULA, [r.id for r in rows], LAYER_FORMULA, level, statement)
    if level:
        return _answer(ANSWER_MASS, [], "", level, f"Only the exact mass is available: precursor m/z {float(precursor_mz):.4f}. "
                                                    "No formula or structure could be proposed.")
    return _answer(ANSWER_NONE, [], "", None, "No annotation: no structure, formula or precursor mass is available.")


def _answer(kind, ids, layer, level, statement, top_id=None, verdict=None):
    return {"type": kind, "ids": ids, "layer": layer, "layer_label": LAYER_LABELS.get(layer, ""), "level": level or "",
            "level_label": level_label(level), "statement": statement, "top_id": top_id, "isomer_verdict": verdict}


def _structure_answer(structural, isomer_result):
    result = isomer_result or {}
    verdict = result.get("verdict") or VERDICT_NOT_RUN
    by_id = {c.id: c for c in structural}
    top = by_id.get(result.get("top_id")) or structural[0]
    where = f"{LAYER_LABELS[top.layer]}, {level_name(top.level)}"
    group = [top.id] + [i for i in _ids(result.get("group")) if i != top.id and i in by_id]

    if verdict == VERDICT_GROUP and len(group) > 1:
        names = [_label(by_id[i]) for i in group]
        statement = f"The spectrum cannot separate {_plural(len(group), 'isomer')}: {_join(names)}. The top-ranked one is {names[0]} ({where})."
        kind, ids = ANSWER_GROUP, group
    else:
        statement = f"{_label(top)} is the top-ranked structure ({where}). "
        if verdict == VERDICT_SINGLE:
            statement += "Fragment peaks tell it apart from every other isomer that was considered."
        elif not result:
            statement += "The isomer check was not run."
        else:
            notes = result.get("notes") or ["no isomer comparison was possible"]
            why = notes if isinstance(notes, str) else notes[0]
            statement += f"The isomer check could not be done ({why.rstrip('.')})."
        kind, ids = ANSWER_SINGLE, [top.id]

    challengers = [_label(by_id[i]) for i in _ids(result.get("challengers")) if i in by_id and i != top.id]
    if challengers:
        statement += (f" For information only: {_join(challengers)} {'explains' if len(challengers) == 1 else 'explain'} "
                      f"observed peaks that {_label(top)} does not.")
    return _answer(kind, ids, top.layer, top.level, statement, top.id, verdict)
