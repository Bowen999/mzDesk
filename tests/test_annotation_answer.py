"""ESI MS/MS Annotation core, part 3: the answer block and the isomer annotation of the rows (mzdesk/annotation/answer.py).

The isomer result is the dict of mzdesk/exclusivity/isomers.discriminate() (verdict, top_id, group, challengers, notes), built by hand here.
"""
import json

import pytest

from mzdesk.annotation.answer import ISOMER_CAVEAT, attach_isomer, build_answer
from mzdesk.annotation.levels import assign_level
from mzdesk.annotation.merge import formula_only_candidates, merge_candidates
from mzdesk.annotation.schema import Candidate, LayerResult

CAFFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
PARAXANTHINE_LIKE = "Cn1cnc2c1c(=O)[nH]c(=O)n2C"          # a C7H8N4O2 isomer family stand-in (other formula; the verdict is hand-made anyway)
GLUCOSE = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"
POPC = r"CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\CCCCCCCC"


def row(cid, name, layer="library", formula="C6H12O6", smiles=GLUCOSE, gate=True, rank=1):
    c = Candidate(id=cid, layer=layer, smiles=smiles, name=name, formula=formula, adduct="[M+H]+", layer_score=0.9, layer_score_norm=0.9,
                  evidence={layer: {"gate_passed": gate} if layer == "library" else {}}, hypothesis=(layer == "frigid"), rank=rank)
    assign_level(c)
    return c


@pytest.fixture
def ranked():
    """Three isomers (one formula) and a different compound, best first."""
    return [row("c1", "alpha-glucose", rank=1), row("c2", "beta-glucose", "msfinder", rank=2),
            row("c3", "fructose", "msfinder", rank=3), row("c4", "caffeine", formula="C8H10N4O2", smiles=CAFFEINE, rank=4)]


def isomer(verdict, **kw):
    return {"verdict": verdict, "top_id": "c1", "group": [], "challengers": [], "notes": [], "engine": "stub", **kw}


# ---------------------------------------------------------------- single
def test_single_when_the_isomer_check_separates_the_top_candidate(ranked):
    a = build_answer(ranked, isomer("single", group=["c1"]))
    assert (a["type"], a["ids"], a["top_id"], a["isomer_verdict"]) == ("single", ["c1"], "c1", "single")
    assert (a["layer"], a["layer_label"], a["level"]) == ("library", "library", "2a")
    assert a["level_label"] == "Level 2a: probable structure (library spectrum match)"
    assert a["statement"].startswith("alpha-glucose is the top-ranked structure (library, Level 2a).")
    assert "tell it apart from every other isomer" in a["statement"] and ISOMER_CAVEAT in a["statement"]
    assert "only holds for the candidates that were considered" in a["statement"]


def test_the_isomer_verdict_does_not_change_the_level(ranked):
    levels = {v: build_answer(ranked, isomer(v, group=["c1", "c2"] if v == "indistinguishable" else ["c1"]))["level"]
              for v in ("single", "indistinguishable", "insufficient")}
    assert set(levels.values()) == {"2a"} and build_answer(ranked)["level"] == "2a"


@pytest.mark.parametrize("result, expected", [
    (isomer("insufficient", notes=["No MS2 peaks above 1 % of the base peak."]),
     "The isomer check could not be done (No MS2 peaks above 1 % of the base peak), so other isomers of C6H12O6 are not ruled out."),
    (isomer("insufficient", notes="the engine failed."), "The isomer check could not be done (the engine failed), so other isomers of C6H12O6 are not ruled out."),
    (isomer("insufficient"), "The isomer check could not be done (no isomer comparison was possible), so other isomers of C6H12O6 are not ruled out."),
    (None, "The isomer check was not run, so other isomers of C6H12O6 are not ruled out."),
])
def test_an_isomer_check_that_could_not_be_done_still_gives_the_top_candidate_and_says_isomers_are_not_ruled_out(ranked, result, expected):
    a = build_answer(ranked, result)
    assert a["type"] == "single" and a["ids"] == ["c1"] and a["isomer_verdict"] == ("not run" if result is None else "insufficient")
    assert a["statement"] == f"alpha-glucose is the top-ranked structure (library, Level 2a). {expected}"


def test_the_top_candidate_is_the_one_the_isomer_check_was_run_for(ranked):
    a = build_answer(ranked, isomer("single", top_id="c3"))
    assert a["top_id"] == "c3" and a["ids"] == ["c3"] and a["layer"] == "msfinder" and a["level"] == "3"
    unknown = build_answer(ranked, isomer("single", top_id="not-in-the-list"))
    assert unknown["top_id"] == "c1"


def test_a_candidate_without_a_name_is_called_by_its_smiles(ranked):
    ranked[0].name = ""
    assert build_answer(ranked)["statement"].startswith(f"{GLUCOSE} is the top-ranked structure")


def test_a_frigid_only_top_candidate_is_a_hypothesis(ranked):
    hypothesis = [row("c1", "", "frigid", formula="C6H12O6")]
    a = build_answer(hypothesis, isomer("insufficient"))
    assert (a["type"], a["layer"], a["layer_label"], a["level"]) == ("single", "frigid", "FRIGID", "H")
    assert "(FRIGID, Hypothesis)" in a["statement"]


# ---------------------------------------------------------------- group
def test_group_when_the_top_candidate_cannot_be_separated_from_an_isomer(ranked):
    a = build_answer(ranked, isomer("indistinguishable", group=["c3", "c1"]))
    assert (a["type"], a["ids"], a["top_id"], a["isomer_verdict"]) == ("group", ["c1", "c3"], "c1", "indistinguishable")
    assert a["layer"] == "library" and a["level"] == "2a"
    assert a["statement"].startswith("The spectrum cannot separate 2 isomers: alpha-glucose and fructose.")
    assert "The top-ranked one is alpha-glucose (library, Level 2a)." in a["statement"] and a["statement"].endswith(ISOMER_CAVEAT)
    three = build_answer(ranked, isomer("indistinguishable", group=["c1", "c2", "c3"]))
    assert three["ids"] == ["c1", "c2", "c3"] and "alpha-glucose, beta-glucose and fructose" in three["statement"]


def test_a_group_ignores_ids_that_are_not_in_the_ranked_list(ranked):
    a = build_answer(ranked, isomer("indistinguishable", group=["c1", "c9", "c2"]))
    assert a["ids"] == ["c1", "c2"]


def test_challengers_are_named_for_information_only_whatever_their_form(ranked):
    for challengers in (["c2", "c3"], [{"id": "c2", "peaks": [(60.04, 12.0)]}, {"id": "c3"}]):
        a = build_answer(ranked, isomer("single", challengers=challengers, group=["c1"]))
        assert a["type"] == "single" and a["ids"] == ["c1"]
        assert a["statement"].endswith("For information only: beta-glucose and fructose explain observed peaks that alpha-glucose does not.")
    one = build_answer(ranked, isomer("indistinguishable", group=["c1", "c2"], challengers=["c2"]))
    assert "For information only: beta-glucose explains observed peaks" in one["statement"] and one["ids"] == ["c1", "c2"]


# ---------------------------------------------------------------- no structure
def test_formula_only_when_there_is_no_structure_but_a_formula_list():
    a = build_answer([], None, formula_list=[{"formula": "C8H10N4O2", "adduct": "[M+H]+", "prob": 0.8}, "C7H6N6O2"], precursor_mz=195.0877)
    assert (a["type"], a["ids"], a["layer"], a["layer_label"], a["level"], a["top_id"]) == ("formula_only", ["f1", "f2"], "formula_only", "formula-only", "4", None)
    assert a["statement"] == ("No structure could be proposed. The data support 2 molecular formulas: C8H10N4O2 ([M+H]+) and C7H6N6O2. "
                              "This is formula-level evidence only.")
    assert a["ids"] == [r.id for r in formula_only_candidates(["C8H10N4O2", "C7H6N6O2"])]
    one = build_answer([], None, formula_list=["C8H10N4O2"])
    assert "1 molecular formula: C8H10N4O2." in one["statement"]


def test_formula_rows_in_the_ranked_list_do_not_count_as_a_structure():
    rows = formula_only_candidates(["C8H10N4O2"])
    a = build_answer(rows, None, formula_list=["C8H10N4O2"])
    assert a["type"] == "formula_only" and a["ids"] == ["f1"]


def test_mass_only_when_there_is_neither_structure_nor_formula():
    for formulas in (None, []):
        a = build_answer([], None, formula_list=formulas, precursor_mz=195.0877)
        assert (a["type"], a["ids"], a["layer"], a["level"], a["level_label"]) == ("mass_only", [], "", "5", "Level 5: exact mass only")
        assert a["statement"] == "Only the exact mass is available: precursor m/z 195.0877. No formula or structure could be proposed."


def test_none_when_there_is_nothing_at_all():
    a = build_answer([], None)
    assert (a["type"], a["ids"], a["layer"], a["level"], a["level_label"], a["top_id"], a["isomer_verdict"]) == ("none", [], "", "", "", None, None)
    assert a["statement"] == "No annotation: no structure, formula or precursor mass is available."


def test_every_answer_is_strict_json(ranked):
    answers = [build_answer(ranked, isomer("single", group=["c1"])), build_answer(ranked, isomer("indistinguishable", group=["c1", "c2"])),
               build_answer(ranked), build_answer([], None, formula_list=["C8H10N4O2"]), build_answer([], None, precursor_mz=195.0877), build_answer([], None)]
    for a in answers:
        json.dumps(a, allow_nan=False)
        assert set(a) == {"type", "ids", "layer", "layer_label", "level", "level_label", "statement", "top_id", "isomer_verdict"}
    assert [a["type"] for a in answers] == ["single", "group", "single", "formula_only", "mass_only", "none"]


# ---------------------------------------------------------------- the isomer annotation of the rows
def test_attach_isomer_marks_the_top_the_group_the_separated_isomers_and_the_challengers(ranked):
    result = isomer("indistinguishable", group=["c1", "c2"], challengers=[{"id": "c3", "name": "fructose", "peaks": []}, {"id": "c2"}],
                    isomers=[{"id": "c2", "distinguished": False}, {"id": "c3", "distinguished": True}],
                    candidates=[{"id": "c1", "explained_fraction": 0.91}, {"id": "c2", "explained_fraction": 0.88}, {"id": "c3", "explained_fraction": None}])
    out = attach_isomer(ranked, result)
    assert out is ranked
    assert [c.isomer for c in ranked] == [
        {"verdict": "indistinguishable", "role": "top", "challenger": False, "explained_fraction": 0.91},
        {"verdict": "indistinguishable", "role": "not separated", "challenger": True, "explained_fraction": 0.88},
        {"verdict": "indistinguishable", "role": "separated", "challenger": True, "explained_fraction": None},
        {"verdict": "indistinguishable", "role": None, "challenger": False, "explained_fraction": None}]
    assert ranked[0].to_dict()["isomer"]["role"] == "top"


def test_attach_isomer_keeps_an_isomer_the_engine_did_not_score_in_the_group(ranked):
    result = isomer("indistinguishable", group=["c1", "c2"], isomers=[{"id": "c2", "distinguished": None}])
    attach_isomer(ranked, result)
    assert [c.isomer["role"] for c in ranked] == ["top", "not separated", None, None]
    attach_isomer(ranked, isomer("insufficient"))                                # a failed check: only the top-1 and the verdict are known
    assert [c.isomer["role"] for c in ranked] == ["top", None, None, None] and ranked[1].isomer["verdict"] == "insufficient"


def test_attach_isomer_without_a_result_leaves_the_rows_alone(ranked):
    attach_isomer(ranked, None)
    assert all(c.isomer is None for c in ranked)


def test_the_real_isomer_check_without_peaks_gives_a_result_the_answer_can_use(ranked):
    """discriminate() stops before any engine when there is no MS2 peak: its real output has the keys this module reads."""
    from mzdesk.exclusivity.isomers import discriminate
    result = discriminate([{"id": c.id, "smiles": c.smiles, "name": c.name} for c in ranked], [], 181.0707, "[M+H]+")
    assert result["verdict"] == "insufficient" and result["group"] == [] and result["challengers"] == [] and result["notes"][0] == "No MS2 peaks."
    attach_isomer(ranked, result)
    a = build_answer(ranked, result)
    assert a["type"] == "single" and a["ids"] == ["c1"] and a["isomer_verdict"] == "insufficient"
    assert "The isomer check could not be done (No MS2 peaks), so other isomers of C6H12O6 are not ruled out." in a["statement"]
    assert all(c.isomer == {"verdict": "insufficient", "role": None, "challenger": False, "explained_fraction": None} for c in ranked)


# ---------------------------------------------------------------- the whole flow
def test_from_the_layers_to_the_answer(monkeypatch):
    """merge -> (rules stubbed) rerank -> isomer annotation -> answer."""
    from mzdesk.annotation import rerank as rr
    monkeypatch.setattr(rr, "_omics_check", lambda smiles: ("metabolome", "Stub verdict."))
    monkeypatch.setattr(rr, "_ion_check", lambda smiles, sign, ph: ("strong", "stub"))
    monkeypatch.setattr(rr, "_adduct_check", lambda smiles, adduct, additive, alkali: ("high", 80.0, "stub"))

    def layer_cand(layer, smiles, name, norm, gate=True):
        return Candidate(id=f"{layer}-{name}", layer=layer, smiles=smiles, name=name, formula="C8H10N4O2" if smiles == CAFFEINE else "C6H12O6",
                         adduct="[M+H]+", layer_score=norm, layer_score_norm=norm, hypothesis=(layer == "frigid"),
                         evidence={layer: {"gate_passed": gate} if layer == "library" else {}})
    merged = merge_candidates([LayerResult("library", candidates=[layer_cand("library", CAFFEINE, "caffeine", 0.93)]),
                               LayerResult("frigid", candidates=[layer_cand("frigid", POPC, "", 0.5), layer_cand("frigid", CAFFEINE, "", 0.9)])])
    ranked = rr.rerank(merged, polarity="Positive", adduct="[M+H]+", omics="Metabolomics")
    result = isomer("single", top_id=ranked[0].id, group=[ranked[0].id])
    attach_isomer(ranked, result)
    a = build_answer(ranked, result, precursor_mz=195.0877)
    assert [c.layer for c in ranked] == ["library", "frigid"] and ranked[0].base_score == 0.98
    assert (a["type"], a["ids"], a["layer"], a["level"]) == ("single", [ranked[0].id], "library", "2a")
    assert "caffeine is the top-ranked structure" in a["statement"] and ranked[0].isomer["role"] == "top"
