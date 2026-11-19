"""Step 5 over several adducts (owner's design 2026-10-30): FRIGID runs on the #1 formula of each of the best 2 adducts and the runs are one layer.

The FRIGID engine (fe.predict) is replaced by hand-made answers; nothing starts. pick_targets() chooses the (formula, adduct) pairs,
run_frigid_layers() runs them one after the other and merges them.
"""
from types import SimpleNamespace

import pytest

from mzdesk.annotation import frigid_step
from mzdesk.annotation.merge import merge_candidates
from mzdesk.ms2 import frigid_engine as fe

PEAKS = [(100.0, 1000.0), (150.0, 500.0)]
CAFFEINE, ETHER, THIRD = ("Cn1c(=O)c2c(ncn2C)n(C)c1=O", "RYYVLZVUVIJVGH"), ("CCOCC", "RYPKRALMXUUSPQ"), ("CCO", "LFQSCWFLJHTTHZ")


def cand(smiles_ik, score, rank, formula):
    return {"smiles": smiles_ik[0], "inchikey_2d": smiles_ik[1], "score": score, "rank": rank, "formula": formula}


def result(formula, adduct, candidates):
    return {"formula": formula, "adduct": adduct, "driver_adduct": adduct, "mode": "quick", "n_samples": 16, "n_formula_matched": len(candidates),
            "n_unique_matched": len(candidates), "device": "cpu", "candidates": candidates}


def f(formula, adduct, rank=1):
    return SimpleNamespace(formula=formula, adduct=adduct, rank=rank)


@pytest.fixture
def engine(monkeypatch):
    """install({(formula, adduct): result | Exception}) -> the list of (formula, adduct) calls, in order."""
    calls = []

    def install(answers):
        def fake(formula, precursor_mz, adduct, peaks, mode="quick", **kw):
            calls.append((formula, adduct))
            answer = answers[(formula, adduct)]
            if isinstance(answer, Exception):
                raise answer
            return answer
        monkeypatch.setattr(fe, "predict", fake)
        return calls
    return install


# ---------------------------------------------------------------- pick_targets
def test_one_formula_per_adduct_in_adduct_score_order_at_most_two():
    # round-robin list: the order between adducts must not matter, only "best first within one adduct"
    formulas = [f("C8H10N4O2", "[M+H]+"), f("C7H6N4O3", "[M+Na]+"), f("C8H10N4O2X", "[M+H]+", 2), f("C9H9N5O", "[M+K]+"), f("C7H7N5O2", "[M+Na]+", 2)]
    assert frigid_step.pick_targets(formulas, ["[M+Na]+", "[M+H]+", "[M+K]+"], "[M+H]+") == [("C7H6N4O3", "[M+Na]+"), ("C8H10N4O2", "[M+H]+")]
    assert frigid_step.pick_targets(formulas, ["[M+H]+", "[M+Na]+", "[M+K]+"], "[M+H]+") == [("C8H10N4O2", "[M+H]+"), ("C7H6N4O3", "[M+Na]+")]
    assert frigid_step.MAX_ADDUCTS == 2


def test_adducts_without_formulas_are_left_out_and_no_formula_gives_the_default():
    formulas = [f("C8H10N4O2", "[M+Na]+")]
    assert frigid_step.pick_targets(formulas, ["[M+H]+", "[M+Na]+"], "[M+H]+") == [("C8H10N4O2", "[M+Na]+")]
    assert frigid_step.pick_targets([], ["[M+H]+"], "[M+H]+") == [(None, "[M+H]+")]


# ---------------------------------------------------------------- run_frigid_layers
def test_two_adducts_make_two_calls_one_merged_deduped_layer(engine):
    calls = engine({("C8H10N4O2", "[M+H]+"): result("C8H10N4O2", "[M+H]+", [cand(CAFFEINE, 0.60, 1, "C8H10N4O2"), cand(ETHER, 0.55, 2, "C8H10N4O2")]),
                    ("C4H10O", "[M+Na]+"): result("C4H10O", "[M+Na]+", [cand(ETHER, 0.90, 1, "C4H10O"), cand(THIRD, 0.40, 2, "C4H10O")])})
    layer, warns = frigid_step.run_frigid_layers(100.0, [("C8H10N4O2", "[M+H]+"), ("C4H10O", "[M+Na]+")], PEAKS)
    assert calls == [("C8H10N4O2", "[M+H]+"), ("C4H10O", "[M+Na]+")] and warns == []
    assert layer.layer == "frigid" and layer.status == "done"
    assert [c.id for c in layer.candidates] == ["frigid_1", "frigid_2", "frigid_3"]
    assert [(c.ik14, c.layer_score, c.formula, c.adduct) for c in layer.candidates] == [
        (ETHER[1], 0.90, "C4H10O", "[M+Na]+"), (CAFFEINE[1], 0.60, "C8H10N4O2", "[M+H]+"), (THIRD[1], 0.40, "C4H10O", "[M+Na]+")]
    assert all(c.hypothesis and c.layer == "frigid" for c in layer.candidates)
    assert "C8H10N4O2 as [M+H]+" in layer.note and "C4H10O as [M+Na]+" in layer.note and "3 structures in all" in layer.note
    assert [c.adduct for c in merge_candidates([layer])] == ["[M+Na]+", "[M+H]+", "[M+Na]+"]          # the per-candidate adduct survives the merge


def test_one_failed_run_still_gives_a_done_layer_and_a_warning(engine):
    engine({("C8H10N4O2", "[M+H]+"): fe.FrigidError("FRIGID is busy"),
            ("C4H10O", "[M+Na]+"): result("C4H10O", "[M+Na]+", [cand(ETHER, 0.90, 1, "C4H10O")])})
    layer, warns = frigid_step.run_frigid_layers(100.0, [("C8H10N4O2", "[M+H]+"), ("C4H10O", "[M+Na]+")], PEAKS)
    assert layer.status == "done" and [c.adduct for c in layer.candidates] == ["[M+Na]+"]
    assert "failed" in layer.note and "FRIGID is busy" in layer.note
    assert len(warns) == 1 and "C8H10N4O2 as [M+H]+" in warns[0] and "FRIGID is busy" in warns[0]


def test_all_runs_failed_fails_the_layer(engine):
    engine({("C8H10N4O2", "[M+H]+"): fe.FrigidError("one"), ("C4H10O", "[M+Na]+"): fe.FrigidError("two")})
    layer, warns = frigid_step.run_frigid_layers(100.0, [("C8H10N4O2", "[M+H]+"), ("C4H10O", "[M+Na]+")], PEAKS)
    assert layer.status == "failed" and not layer.candidates and warns == [] and "one" in layer.note and "two" in layer.note


def test_unsupported_adduct_is_skipped_with_a_note_without_failing(engine):
    calls = engine({("C8H10N4O2", "[M+H]+"): result("C8H10N4O2", "[M+H]+", [cand(CAFFEINE, 0.9, 1, "C8H10N4O2")])})
    layer, warns = frigid_step.run_frigid_layers(100.0, [("C8H10N4O2", "[M+H]+"), ("C4H10O", "[M+2H]2+")], PEAKS)
    assert calls == [("C8H10N4O2", "[M+H]+")] and warns == []                  # the engine was not called for the unsupported adduct
    assert layer.status == "done" and len(layer.candidates) == 1
    assert "skipped" in layer.note and "does not support" in layer.note


def test_negative_mode_is_skipped_without_any_engine_call(engine):
    calls = engine({})
    layer, warns = frigid_step.run_frigid_layers(100.0, [("C8H10N4O2", "[M-H]-"), ("C4H10O", "[M+HCOO]-")], PEAKS)
    assert calls == [] and warns == []
    assert layer.status == "skipped" and not layer.candidates and "positive ions only" in layer.note


def test_a_single_target_is_exactly_the_single_run(engine):
    engine({("C8H10N4O2", "[M+H]+"): result("C8H10N4O2", "[M+H]+", [cand(CAFFEINE, 0.9, 1, "C8H10N4O2")])})
    layer, warns = frigid_step.run_frigid_layers(100.0, [("C8H10N4O2", "[M+H]+")], PEAKS)
    single = frigid_step.run_frigid_layer(100.0, "[M+H]+", PEAKS, "C8H10N4O2")
    assert (layer.status, layer.note, [c.id for c in layer.candidates]) == (single.status, single.note, [c.id for c in single.candidates]) and warns == []
