"""ESI MS/MS Annotation core, part 1: the Candidate / LayerResult contract (schema.py), the confidence levels (levels.py) and the
merge of the layers' candidates (merge.py). Pure Python: synthetic candidates, no engine, no network.
"""
import json
import math

import pytest

from mzdesk.annotation import levels
from mzdesk.annotation.merge import formula_only_candidates, has_structure, merge_candidates, structure_keys
from mzdesk.annotation.schema import (FLAG_NO_STRUCTURE, LAYER_LABELS, Candidate, LayerResult, evidence_of)

CAFFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
CAFFEINE_KEKULE = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"          # the same molecule written another way
CAFFEINE_IK14 = "RYYVLZVUVIJVGH"
GLUCOSE = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"
L_LEUCINE = "CC(C)C[C@H](N)C(=O)O"
D_LEUCINE = "CC(C)C[C@@H](N)C(=O)O"
ISOLEUCINE = "CC[C@H](C)[C@H](N)C(=O)O"
POPC = r"CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\CCCCCCCC"
SALT = "[Na+].CC(=O)[O-]"
INVALID = "this is not a smiles"


def cand(layer, smiles, norm, *, cid=None, name="", gate=True, **kw):
    """A candidate as a layer would return it (evidence keyed by layer; the library one carries the gate verdict)."""
    evidence = {layer: {"gate_passed": gate, "failed_checks": [] if gate else ["adduct"]}} if layer == "library" else {layer: {}}
    return Candidate(id=cid or f"{layer}-{name or (smiles or 'none')[:6]}", layer=layer, smiles=smiles, name=name, formula=kw.pop("formula", None),
                     adduct=kw.pop("adduct", "[M+H]+"), layer_score=norm, layer_score_norm=norm, evidence=evidence,
                     hypothesis=(layer == "frigid"), **kw)


def lib(*cands, **kw):
    return LayerResult("library", candidates=list(cands), **kw)


def msf(*cands, **kw):
    return LayerResult("msfinder", candidates=list(cands), **kw)


def frg(*cands, **kw):
    return LayerResult("frigid", candidates=list(cands), **kw)


# ---------------------------------------------------------------- dedupe
def test_the_same_structure_from_three_layers_is_one_candidate_that_keeps_every_layer():
    layers = [lib(cand("library", GLUCOSE, 0.74, name="glucose", gate=False), cand("library", CAFFEINE, 0.93, name="caffeine")),
              msf(cand("msfinder", CAFFEINE_KEKULE, 0.80, name="")),
              frg(cand("frigid", ISOLEUCINE, 0.9), cand("frigid", CAFFEINE, 0.6))]
    merged = merge_candidates(layers)
    assert len(merged) == 3                                           # caffeine, glucose, isoleucine
    caf = next(c for c in merged if c.ik14 == CAFFEINE_IK14)
    assert caf.layer == "library" and caf.name == "caffeine"          # the strongest layer supplies the fields
    assert [(s["layer"], s["layer_score_norm"], s["rank_in_layer"]) for s in caf.sources] == [
        ("library", 0.93, 2), ("msfinder", 0.80, 1), ("frigid", 0.6, 2)]
    assert caf.layer_score_norm == 0.93 and caf.inchikey == "RYYVLZVUVIJVGH-UHFFFAOYSA-N"
    assert set(caf.evidence) == {"library", "msfinder", "frigid"}
    assert caf.hypothesis is False                                    # a library hit supports it: no longer a bare hypothesis
    assert caf.level == "2a"
    assert caf.flags == []


def test_ids_are_new_unique_and_the_input_is_not_changed():
    first = cand("library", CAFFEINE, 0.9, cid="same")
    second = cand("msfinder", GLUCOSE, 0.7, cid="same")
    merged = merge_candidates([lib(first), msf(second)])
    assert [c.id for c in merged] == ["c1", "c2"]
    assert first.id == second.id == "same" and first.sources == [] and first.level is None and first.ik14 is None


def test_stereoisomers_merge_but_constitutional_isomers_do_not():
    merged = merge_candidates([lib(cand("library", L_LEUCINE, 0.9, name="L-leucine")), msf(cand("msfinder", D_LEUCINE, 0.7), cand("msfinder", ISOLEUCINE, 0.6))])
    assert len(merged) == 2
    leu = next(c for c in merged if c.ik14 == structure_keys(L_LEUCINE)[1])
    assert [s["layer"] for s in leu.sources] == ["library", "msfinder"] and leu.smiles == L_LEUCINE
    assert structure_keys(ISOLEUCINE)[1] != structure_keys(L_LEUCINE)[1]


def test_a_salt_is_keyed_by_its_largest_fragment():
    assert structure_keys(SALT)[1] == structure_keys("CC(=O)O")[1]
    assert len(merge_candidates([lib(cand("library", SALT, 0.8)), msf(cand("msfinder", "CC(=O)O", 0.7))])) == 1


def test_a_layer_that_repeats_a_compound_counts_once_with_its_best_score():
    layer = lib(cand("library", CAFFEINE, 0.95, cid="a"), cand("library", GLUCOSE, 0.8, cid="b"), cand("library", CAFFEINE_KEKULE, 0.85, cid="c"))
    merged = merge_candidates([layer])
    assert len(merged) == 2
    assert merged[0].sources == [{"layer": "library", "layer_score": 0.95, "layer_score_norm": 0.95, "rank_in_layer": 1}]


def test_the_inchikey_given_by_a_layer_is_used_when_there_is_no_smiles():
    sibling = cand("msfinder", None, 0.7, inchikey="RYYVLZVUVIJVGH-UHFFFAOYSA-N")
    only_block = cand("frigid", None, 0.5, ik14=CAFFEINE_IK14.lower())
    merged = merge_candidates([lib(cand("library", CAFFEINE, 0.9)), msf(sibling), frg(only_block)])
    assert len(merged) == 1 and [s["layer"] for s in merged[0].sources] == ["library", "msfinder", "frigid"]
    assert merged[0].smiles == CAFFEINE and merged[0].flags == []


def test_a_candidate_without_a_valid_smiles_is_kept_flagged_and_never_merged_blindly():
    merged = merge_candidates([lib(cand("library", INVALID, 0.9, name="broken one")),
                               msf(cand("msfinder", None, 0.8, name="no structure"), cand("msfinder", CAFFEINE, 0.7))])
    assert len(merged) == 3
    broken, none_at_all, good = merged
    assert broken.flags == [FLAG_NO_STRUCTURE] and none_at_all.flags == [FLAG_NO_STRUCTURE] and good.flags == []
    assert broken.smiles == INVALID and broken.ik14 is None          # what the layer gave is kept for display
    assert not has_structure(INVALID) and not has_structure("") and not has_structure(None) and has_structure(CAFFEINE)


def test_an_invalid_smiles_takes_the_structure_of_its_sibling_from_another_layer():
    broken = cand("library", INVALID, 0.9, name="caffeine", inchikey="RYYVLZVUVIJVGH-UHFFFAOYSA-N")
    merged = merge_candidates([lib(broken), msf(cand("msfinder", CAFFEINE, 0.7))])
    assert len(merged) == 1
    assert merged[0].smiles == CAFFEINE and merged[0].flags == [] and merged[0].layer == "library" and merged[0].name == "caffeine"


def test_missing_fields_are_filled_from_the_next_layer():
    top = cand("library", CAFFEINE, 0.9, name="", formula=None, adduct=None)
    below = cand("msfinder", CAFFEINE, 0.7, name="1,3,7-trimethylxanthine", formula="C8H10N4O2", adduct="[M+Na]+", mz_error_ppm=0.0)
    merged = merge_candidates([lib(top), msf(below)])[0]
    assert (merged.name, merged.formula, merged.adduct, merged.mz_error_ppm) == ("1,3,7-trimethylxanthine", "C8H10N4O2", "[M+Na]+", 0.0)


def test_hypothesis_stays_only_when_every_contributor_is_one():
    both = merge_candidates([frg(cand("frigid", CAFFEINE, 0.8))])[0]
    mixed = merge_candidates([msf(cand("msfinder", CAFFEINE, 0.5)), frg(cand("frigid", CAFFEINE, 0.8))])[0]
    assert both.hypothesis is True and mixed.hypothesis is False


def test_skipped_and_failed_layers_add_nothing():
    layers = [lib(cand("library", CAFFEINE, 0.9)), msf(status="skipped", note="gate passed"), frg(status="failed", note="not available")]
    assert [c.id for c in merge_candidates(layers)] == ["c1"]
    assert merge_candidates([]) == [] and merge_candidates([lib()]) == []


def test_a_formula_row_cannot_sit_in_a_layer_result():
    row = Candidate(id="f1", layer="formula_only", formula="C8H10N4O2")
    with pytest.raises(ValueError, match="formula_only"):
        merge_candidates([lib(row)])


# ---------------------------------------------------------------- levels
@pytest.mark.parametrize("layers, level", [
    ([("library", True)], "2a"),
    ([("library", False)], "3"),
    ([("msfinder", None)], "3"),
    ([("frigid", None)], "H"),
    ([("library", False), ("frigid", None)], "3"),
    ([("msfinder", None), ("frigid", None)], "3"),
    ([("library", True), ("msfinder", None), ("frigid", None)], "2a"),
])
def test_the_merged_candidate_gets_the_best_level_of_its_layers(layers, level):
    results = []
    for layer, gate in layers:
        c = cand(layer, CAFFEINE, 0.8, gate=bool(gate))
        results.append(LayerResult(layer, candidates=[c]))
    merged = merge_candidates(results)
    assert len(merged) == 1 and merged[0].level == level and merged[0].level_reason.endswith(".")


def test_the_level_reason_names_the_gate_checks_that_failed_and_the_other_layers():
    failed = cand("library", CAFFEINE, 0.74)
    failed.evidence["library"].update(gate_passed=False, failed_checks=["adduct not retained", "isotope similarity 0.62"])
    merged = merge_candidates([lib(failed), msf(cand("msfinder", CAFFEINE, 0.8))])[0]
    assert "adduct not retained; isotope similarity 0.62" in merged.level_reason and "also proposed by MS-FINDER" in merged.level_reason
    passed = cand("library", CAFFEINE, 0.93)
    passed.evidence["library"].update(similarity=0.93, matched_peaks=9)
    assert "similarity 0.93, 9 matched peaks" in merge_candidates([lib(passed)])[0].level_reason


def test_assign_level_works_on_a_candidate_a_layer_has_not_merged_yet():
    c = cand("library", CAFFEINE, 0.9)
    assert levels.assign_level(c) == "2a" and c.level == "2a" and c.level_reason
    flat = Candidate(id="x", layer="library", evidence={"gate_passed": True})        # evidence a layer forgot to key by layer
    assert evidence_of(flat, "library") == {"gate_passed": True} and levels.assign_level(flat) == "2a"


def test_fallback_levels_and_labels():
    assert levels.fallback_level(True, True) == "4" and levels.fallback_level(True, False) == "4"
    assert levels.fallback_level(False, True) == "5" and levels.fallback_level(False, False) is None
    assert [levels.level_name(x) for x in levels.LEVELS] == ["Level 2a", "Level 3", "Hypothesis", "Level 4", "Level 5"]
    assert levels.level_label("3") == "Level 3: tentative candidate" and levels.level_label(None) == ""
    assert levels.LEVELS == ("2a", "3", "H", "4", "5")


# ---------------------------------------------------------------- formula-only rows
def test_formula_only_rows_are_level_4_ranked_in_order_and_not_scored():
    rows = formula_only_candidates(["C8H10N4O2", {"formula": "C7H6N6O2", "adduct": "[M+H]+", "prob": 0.31, "ppm": -1.2, "isotope_sim": None}])
    assert [(r.id, r.layer, r.level, r.rank, r.formula) for r in rows] == [
        ("f1", "formula_only", "4", 1, "C8H10N4O2"), ("f2", "formula_only", "4", 2, "C7H6N6O2")]
    assert rows[1].layer_score == 0.31 and rows[1].mz_error_ppm == -1.2 and rows[1].evidence == {"formula_only": {"prob": 0.31, "ppm": -1.2}}
    assert all(r.flags == [FLAG_NO_STRUCTURE] and r.final_score is None and r.smiles is None for r in rows)
    assert formula_only_candidates(None) == [] and formula_only_candidates([]) == []


# ---------------------------------------------------------------- schema
def test_to_dict_is_strict_json_with_display_names():
    class NumpyLike:                      # np.float32 and friends are not float subclasses; they have .item()
        def item(self):
            return 0.5

    c = cand("library", CAFFEINE, 0.9, mz_error_ppm=float("nan"))
    c.evidence["library"].update(similarity=float("inf"), score=NumpyLike(), peaks=((1.0, 2.0),))
    merged = merge_candidates([lib(c)])[0]
    d = merged.to_dict()
    text = json.dumps(d, allow_nan=False)
    assert json.loads(text)["id"] == "c1"
    assert d["mz_error_ppm"] is None and d["evidence"]["library"]["similarity"] is None
    assert d["evidence"]["library"]["score"] == 0.5 and d["evidence"]["library"]["peaks"] == [[1.0, 2.0]]
    assert (d["layer_label"], d["level_label"]) == ("library", "Level 2a: probable structure (library spectrum match)")
    for field in ("id", "layer", "smiles", "inchikey", "ik14", "name", "formula", "adduct", "mz_error_ppm", "layer_score", "layer_score_norm",
                  "evidence", "hypothesis", "sources", "detect", "final_score", "rank", "level", "level_reason", "isomer"):
        assert field in d
    assert math.isfinite(d["layer_score_norm"])


def test_adduct_aliases_and_a_missing_score_are_normalised_when_the_candidate_is_made():
    c = Candidate(id="x", layer="msfinder", adduct="[M+FA-H]-", layer_score_norm=None)
    assert c.adduct == "[M+HCOO]-" and c.layer_score_norm == 0.0
    assert Candidate(id="y", layer="library", adduct="[M+H]+").adduct == "[M+H]+" and Candidate(id="z", layer="library").adduct is None
    merged = merge_candidates([msf(Candidate(id="m", layer="msfinder", smiles=CAFFEINE, adduct="[M+Hac-H]-"))])
    assert merged[0].adduct == "[M+CH3COO]-"


def test_layer_result_to_dict_and_validation():
    result = frg(cand("frigid", CAFFEINE, 0.5), status="done", note="1 hypothesis", seconds=31.2)
    d = result.to_dict()
    assert (d["layer"], d["status"], d["note"], d["seconds"]) == ("frigid", "done", "1 hypothesis", 31.2) and d["candidates"][0]["layer_label"] == "FRIGID"
    json.dumps(d, allow_nan=False)
    with pytest.raises(ValueError, match="status"):
        LayerResult("library", status="running")
    with pytest.raises(ValueError, match="layer"):
        LayerResult("formula_only")
    with pytest.raises(ValueError, match="layer"):
        Candidate(id="x", layer="metfrag")
    assert LAYER_LABELS == {"library": "library", "msfinder": "MS-FINDER", "frigid": "FRIGID", "formula_only": "formula-only"}
