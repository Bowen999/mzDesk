"""Owner rules of the Level-2 gate: a formula missing from BUDDY's short list but matching the precursor is not a failure,
and a library hit with no structure is never better than Level 3."""
import pytest

from mzdesk.annotation import library_step as L
from mzdesk.annotation.levels import assign_level
from mzdesk.annotation.merge import merge_candidates
from mzdesk.annotation.schema import LayerResult
from mzdesk.ms1 import msbuddy_engine as bd

CAFFEINE_SMILES = "CN1C=NC2=C1C(=O)N(C)C(=O)N2C"
CAFFEINE_KEY = "RYYVLZVUVIJVGH-UHFFFAOYSA-N"
THEOPHYLLINE = "Cn1c2nc[nH]c2c(=O)n(C)c1=O"                  # same formula as caffeine minus CH2: C7H8N4O2
BUDDY_TOP = ["C7H14O6", "C6H10N5O2", "C9H14N2O3"]


def hit(**kw):
    base = {"rank": 1, "library": "MassBank", "row": 0, "accession": "ACC1", "name": "Caffeine", "adduct": "[M+H]+", "precursor_mz": 195.0877,
            "ppm": 0.0, "matches": 12, "scores": {L.PRIMARY_METRIC: 0.93}, "instrument": "QTOF", "collision_energy": "30",
            "inchikey": CAFFEINE_KEY, "smiles": CAFFEINE_SMILES}
    base.update(kw)
    return base


def gate(**kw):
    args = {"precursor_mz": bd.ion_mz("C8H10N4O2", "[M+H]+"), "tol_ppm": 10.0, "n_query": 15, "adducts": ["[M+H]+"], "ms1": None,
            "formulas": [{"formula": f, "adduct": "[M+H]+"} for f in BUDDY_TOP]}
    args.update(kw)
    return L._Gate(**args)


def formula_check(verdict):
    return next(c for c in verdict["checks"] if c["id"] == "formula")


# ---------------------------------------------------------------- rule 1: formula outside BUDDY's list
def test_formula_missing_from_the_buddy_list_but_matching_the_precursor_is_not_evaluated():
    v = gate().judge(hit())
    check = formula_check(v)
    assert check["status"] == "not evaluated" and v["passed"]
    assert "C8H10N4O2 is not among BUDDY's top 3 formulas for [M+H]+" in check["detail"] and "matches the precursor at +0.0 ppm" in check["detail"]


def test_formula_missing_from_the_buddy_list_and_not_matching_the_precursor_fails():
    v = gate(precursor_mz=bd.ion_mz("C8H10N4O2", "[M+H]+") + 0.01).judge(hit())         # about 51 ppm off
    assert formula_check(v)["status"] == "fail" and not v["passed"]
    assert "is not among the BUDDY formulas for [M+H]+" in formula_check(v)["detail"]


def test_the_mass_match_uses_the_formula_not_the_library_recorded_mz():
    # the library records an m/z that fits the query, but the hit's own formula does not (caffeine mass against a theophylline-sized precursor)
    v = gate(precursor_mz=181.0720).judge(hit(ppm=0.0, precursor_mz=181.0720))
    assert v["error_from_formula"] and formula_check(v)["status"] == "fail" and not v["passed"]
    # without a computable adduct the error is the library's own, which must not rescue the formula
    no_adduct = gate(formulas=BUDDY_TOP).judge(hit(adduct="not an adduct"))             # plain strings count for every adduct
    assert not no_adduct["error_from_formula"] and formula_check(no_adduct)["status"] == "fail"


def test_formula_in_the_buddy_list_still_passes():
    v = gate(formulas=[{"formula": "C8H10N4O2", "adduct": "[M+H]+"}]).judge(hit())
    assert formula_check(v)["status"] == "pass"


# ---------------------------------------------------------------- rule 2: no structure, at most Level 3
def test_a_hit_with_no_smiles_and_no_inchikey_is_level_3_even_when_the_gate_passes():
    g = gate(formulas=[])
    bare = hit(smiles=None, inchikey=None, adduct="[M+H]+")
    v = g.judge(bare)
    assert v["passed"]
    cand = L._candidate("lib1", bare, v, 1, L.PRIMARY_METRIC)
    assert cand.evidence["library"]["gate_passed"] and cand.evidence["library"]["has_structure"] is False
    assert assign_level(cand) == "3"
    assert "no structure" in cand.level_reason and "stays at Level 3" in cand.level_reason
    merged = merge_candidates([LayerResult("library", candidates=[cand])])[0]
    assert merged.level == "3" and "no structure" in merged.level_reason


def test_a_hit_with_a_smiles_or_an_inchikey_stays_level_2a():
    g = gate(formulas=[])
    for given in (hit(), hit(smiles=None), hit(inchikey=None), hit(smiles="N/A", inchikey=CAFFEINE_KEY)):
        cand = L._candidate("lib1", given, g.judge(given), 1, L.PRIMARY_METRIC)
        assert cand.evidence["library"]["has_structure"] is True and assign_level(cand) == "2a", given


def test_a_structureless_hit_ranks_after_a_passing_structured_one(monkeypatch):
    bare = hit(smiles=None, inchikey=None, name="No structure", accession="A", scores={L.PRIMARY_METRIC: 0.99})
    full = hit(name="Caffeine", accession="B", scores={L.PRIMARY_METRIC: 0.90})
    res = {"hits": [bare, full], "n_candidates": 2, "n_kept": 2, "warnings": [], "notes": {}, "n_metrics_scored": 1}
    monkeypatch.setattr(L, "_search", lambda *a, **k: (res, L.PRIMARY_METRIC))
    monkeypatch.setattr(L, "_reference_spectrum", lambda h: None)
    layer, summary = L.run_library([(100.0, 10.0)], bd.ion_mz("C8H10N4O2", "[M+H]+"), True, ["MassBank"], adducts=["[M+H]+"], formulas=[])
    assert [c.name for c in layer.candidates] == ["Caffeine", "No structure"] and summary.passed and summary.best_id == "lib1"
    assert [assign_level(c) for c in layer.candidates] == ["2a", "3"]


# ---------------------------------------------------------------- the metric and the gate threshold are settings of the step
def test_the_gate_uses_the_chosen_threshold_and_metric_and_lists_from_threshold_minus_0_10():
    g = gate(metric="CosineGreedy", gate_similarity=0.90)
    assert g.thresholds()["metric"] == "CosineGreedy" and g.thresholds()["similarity"] == 0.9 and g.thresholds()["list"] == 0.8
    below = g.judge(hit(scores={"CosineGreedy": 0.85, L.PRIMARY_METRIC: 0.99}))                 # only the chosen metric counts
    assert [c["status"] for c in below["checks"] if c["id"] == "similarity"] == ["fail"]
    assert "CosineGreedy 0.85 < 0.90" in below["checks"][0]["detail"]
    assert gate(metric="CosineGreedy", gate_similarity=0.90).judge(hit(scores={"CosineGreedy": 0.92}))["checks"][0]["status"] == "pass"
    assert gate().thresholds()["similarity"] == L.GATE_SIMILARITY and gate().thresholds()["list"] == L.LIST_THRESHOLD == 0.7


def test_the_list_threshold_is_the_gate_minus_0_10_and_never_negative():
    assert L.list_threshold(0.8) == 0.7 and L.list_threshold(1.0) == 0.9 and L.list_threshold(0.05) == 0.0 and L.list_threshold(0) == 0.0


def test_an_unavailable_score_fails_the_similarity_check_instead_of_crashing():
    v = gate(metric="ModifiedCosineGreedy").judge(hit(scores={"ModifiedCosineGreedy": None}))
    assert v["checks"][0]["status"] == "fail" and "n/a" in v["checks"][0]["detail"]


def test_run_library_hands_the_metric_and_the_list_threshold_to_the_search(monkeypatch):
    seen = {}
    modified = "ModifiedCosineGreedy"

    def fake_search(*args):
        seen["args"] = args
        return {"hits": [hit(scores={modified: 0.86})], "n_candidates": 1, "n_kept": 1, "warnings": [], "notes": {}, "n_metrics_scored": 1}, modified
    monkeypatch.setattr(L, "_search", fake_search)
    monkeypatch.setattr(L, "_reference_spectrum", lambda h: None)
    layer, summary = L.run_library([(100.0, 10.0)], bd.ion_mz("C8H10N4O2", "[M+H]+"), True, ["MassBank"], adducts=["[M+H]+"], formulas=[],
                                   metric=modified, gate_similarity=0.85)
    assert seen["args"][-2:] == (modified, 0.75)
    assert summary.passed and summary.thresholds["metric"] == modified and summary.thresholds["similarity"] == 0.85 and summary.thresholds["list"] == 0.75
    assert "scored 0.75 or higher on ModifiedCosineGreedy" in layer.note and "the gate needs 0.85" in layer.note
    assert layer.candidates[0].evidence["library"]["metric"] == modified and not summary.warnings


def test_a_metric_that_could_not_be_computed_is_reported_as_the_fallback(monkeypatch):
    res = {"hits": [hit(scores={"CosineGreedy": 0.93})], "n_candidates": 1, "n_kept": 1, "warnings": [], "notes": {}, "n_metrics_scored": 1}
    monkeypatch.setattr(L, "_search", lambda *a: (res, "CosineGreedy"))
    monkeypatch.setattr(L, "_reference_spectrum", lambda h: None)
    _, summary = L.run_library([(100.0, 10.0)], bd.ion_mz("C8H10N4O2", "[M+H]+"), True, ["MassBank"], adducts=["[M+H]+"], formulas=[],
                               metric="SpectralEntropy (unweighted)")
    assert summary.thresholds["metric"] == "CosineGreedy" and any("SpectralEntropy (unweighted) could not be computed" in w for w in summary.warnings)


@pytest.mark.parametrize("kw", [{"metric": "Spec2Vec"}, {"metric": ""}, {"gate_similarity": 1.1}, {"gate_similarity": -0.01}, {"gate_similarity": "0.8"},
                                {"gate_similarity": float("nan")}])
def test_run_library_rejects_an_unknown_metric_and_a_threshold_outside_0_to_1(kw):
    with pytest.raises(ValueError, match="similarity"):
        L.run_library([(100.0, 10.0)], 195.0877, True, ["MassBank"], adducts=["[M+H]+"], formulas=[], **kw)
