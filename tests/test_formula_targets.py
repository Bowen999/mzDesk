"""Step 2 target adducts and list order: the confident-adduct rule (EXTRA_ADDUCT_PROB) and the round-robin interleaving by BUDDY rank.

BUDDY is replaced by a fake with a per-adduct table; adduct_step is replaced by a fake that returns given hypotheses, so the tests do not
depend on the real adduct priors.
"""
import pytest

from mzdesk.annotation import adduct_step as A
from mzdesk.annotation import formula_step as F
from mzdesk.ms1 import msbuddy_engine as bd

NAMES = ["[M+H]+", "[M+Na]+", "[M+NH4]+", "[M+K]+"]


class FakeBuddy:
    def __init__(self, table):
        self.table, self.calls = table, []

    def __call__(self, precursor_mz, positive, adduct, rt, formula_input, ms2_peaks, ms1_peaks, ms_instr="orbitrap", halogen=False):
        self.calls.append(adduct)
        rows = self.table.get(adduct, [])
        return {"candidates": [{"rank": i, "formula": f, "prob": p, "normed_prob": p, "fdr": 0.0, "isotope_sim": None, "mz_theo": 300.0, "ppm": ppm}
                               for i, (f, p, ppm) in enumerate(rows, start=1)], "adduct_used": adduct, "from_input": False, "ms1_used": 0}


def setup(monkeypatch, table, scores, mode="auto"):
    """scores: {adduct: score}, retained in the given order. Returns the fake BUDDY."""
    fake = FakeBuddy(table)
    monkeypatch.setattr(bd, "annotate_formula", fake)

    def score_adducts(precursor_mz, positive, ms2, adduct=None, tol_da=0.02, buddy_fits=None):
        hyps = [A.AdductHypothesis(a, s, rank=i, retained=True) for i, (a, s) in enumerate(scores.items(), start=1)]
        return A.AdductResult(mode=mode, hypotheses=hyps)
    monkeypatch.setattr(A, "score_adducts", score_adducts)
    return fake


def rows(*probs_ppm):
    return [(f"C{10 + i}H20O5", p, ppm) for i, (p, ppm) in enumerate(probs_ppm)]


def run(**kw):
    return F.run_adduct_and_formula(300.0, True, [(100.0, 1.0)], **kw)


def test_a_far_retained_adduct_with_a_confident_best_formula_gets_formulas(monkeypatch):
    table = {"[M+H]+": rows((0.8, 0.5), (0.1, 1.0)), "[M+Na]+": rows((0.95, 0.5), (0.02, 1.0)), "[M+NH4]+": rows((0.5, 0.5))}
    fake = setup(monkeypatch, table, {"[M+H]+": 1.0, "[M+NH4]+": 0.3, "[M+Na]+": 0.22})         # both far below 0.5 x best
    _, result = run()
    assert result.adducts == ["[M+H]+", "[M+Na]+"]                                           # adduct score order, [M+NH4]+ (prob 0.5) stays out
    assert max(fake.calls.count(a) for a in fake.calls) == 1                                  # BUDDY ran once per adduct (cached)
    assert "2 formulas for [M+H]+ and 2 formulas for [M+Na]+" in result.note


def test_the_prob_threshold_is_inclusive_and_uses_the_best_formula_within_the_tolerance(monkeypatch):
    # [M+Na]+: its rank-1 formula misses the 10 ppm tolerance, the best one within it has prob 0.9 -> counts; [M+K]+ 0.8999 does not
    table = {"[M+H]+": rows((0.8, 0.5)), "[M+Na]+": rows((0.99, 50.0), (0.9, 1.0)), "[M+K]+": rows((0.8999, 1.0)), "[M+NH4]+": rows((0.95, 20.0))}
    setup(monkeypatch, table, {"[M+H]+": 1.0, "[M+Na]+": 0.2, "[M+K]+": 0.15, "[M+NH4]+": 0.1})
    _, result = run()
    assert result.adducts == ["[M+H]+", "[M+Na]+"] and [f.formula for f in result.formulas if f.adduct == "[M+Na]+"] == ["C11H20O5"]


def test_interleaving_is_round_robin_by_buddy_rank_with_the_best_adduct_first(monkeypatch):
    table = {"[M+H]+": rows((0.9, 0.1), (0.05, 0.1), (0.01, 0.1)), "[M+Na]+": rows((0.95, 0.1)), "[M+NH4]+": rows((0.92, 0.1), (0.03, 0.1))}
    setup(monkeypatch, table, {"[M+H]+": 1.0, "[M+Na]+": 0.2, "[M+NH4]+": 0.1})
    _, result = run()
    assert result.adducts == ["[M+H]+", "[M+Na]+", "[M+NH4]+"]
    assert [(f.adduct, f.rank) for f in result.formulas] == [("[M+H]+", 1), ("[M+Na]+", 1), ("[M+NH4]+", 1), ("[M+H]+", 2), ("[M+NH4]+", 2),
                                                             ("[M+H]+", 3)]
    assert result.formulas[0].adduct == "[M+H]+" and result.formulas[0].rank == 1
    assert "3 formulas for [M+H]+, 1 formula for [M+Na]+ and 2 formulas for [M+NH4]+" in result.note


def test_ranks_keep_gaps_from_the_tolerance_filter_and_follow_list_position(monkeypatch):
    table = {"[M+H]+": rows((0.9, 0.1), (0.5, 40.0), (0.1, 0.2)), "[M+Na]+": rows((0.95, 0.1), (0.04, 0.3))}
    setup(monkeypatch, table, {"[M+H]+": 1.0, "[M+Na]+": 0.6})
    _, result = run()
    assert [(f.adduct, f.rank) for f in result.formulas] == [("[M+H]+", 1), ("[M+Na]+", 1), ("[M+H]+", 3), ("[M+Na]+", 2)]


def test_non_auto_mode_is_unchanged(monkeypatch):
    table = {"[M+H]+": rows((0.9, 0.1), (0.05, 0.1)), "[M+Na]+": rows((0.99, 0.1))}
    fake = setup(monkeypatch, table, {"[M+H]+": 1.0, "[M+Na]+": 0.1}, mode="user")
    _, result = run(adduct="[M+H]+")
    assert result.adducts == ["[M+H]+"] and [f.rank for f in result.formulas] == [1, 2] and fake.calls == ["[M+H]+"]


def test_second_adduct_rule_is_kept_and_a_failing_target_is_reported(monkeypatch):
    table = {"[M+H]+": rows((0.3, 0.1)), "[M+Na]+": rows((0.2, 0.1)), "[M+K]+": rows((0.2, 0.1))}
    setup(monkeypatch, table, {"[M+H]+": 1.0, "[M+Na]+": 0.5, "[M+K]+": 0.4})                  # 0.5 x best is included, 0.4 is not
    _, result = run()
    assert result.adducts == ["[M+H]+", "[M+Na]+"] and [f.adduct for f in result.formulas] == ["[M+H]+", "[M+Na]+"]


def test_summary_text_with_three_adducts_and_no_formulas():
    r = F.FormulaResult(formulas=[], adducts=["A", "B", "C"])
    status, note = F._summary(r, ["A", "B", "C"], [], False, 10.0)
    assert status == "done" and "for A, B or C;" in note
