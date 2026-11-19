"""The mass-error factor of steps 1-2: a formula that only just fits the precursor tolerance is weaker evidence than one that fits well.

adduct_step.mass_error_factor (shape, scaling with the tolerance), its use in the adduct score (the best formula maximises BUDDY's
probability x factor; the factor is part of that adduct's evidence) and in the formula list (each adduct's formulas ordered by
probability x factor). BUDDY is replaced by a fake with a per-adduct table, as in test_formula_targets.py.
"""
import json
import math

import pytest

from mzdesk.annotation import adduct_step as A
from mzdesk.annotation import formula_step as F
from mzdesk.ms1 import msbuddy_engine as bd

RICH = [(80.0, 20), (100.0, 40), (120.0, 100), (140.0, 60), (160.0, 30), (180.0, 15)]       # 6 fragment peaks >= 5 %: not sparse
NA_ION = (22.9892, 100.0)                                                                    # the Na+ ion: [M+Na]+ evidence


class FakeBuddy:
    """table = {adduct: [(formula, prob, ppm), ...]} in BUDDY's order."""

    def __init__(self, table):
        self.table = table

    def __call__(self, precursor_mz, positive, adduct, rt, formula_input, ms2_peaks, ms1_peaks, ms_instr="orbitrap", halogen=False):
        rows = self.table.get(adduct, [])
        return {"candidates": [{"rank": i, "formula": f, "prob": p, "normed_prob": p, "fdr": 0.0, "isotope_sim": None, "mz_theo": 233.0, "ppm": ppm}
                               for i, (f, p, ppm) in enumerate(rows, start=1)], "adduct_used": adduct, "from_input": False, "ms1_used": 0}


def run(monkeypatch, table, extra=(), **kw):
    monkeypatch.setattr(bd, "annotate_formula", FakeBuddy(table))
    return F.run_adduct_and_formula(233.117, True, RICH + [(233.117, 50.0)] + list(extra), **kw)


# ---------------------------------------------------------------- the function
def test_the_factor_is_one_at_zero_error_and_close_to_the_floor_at_the_limit():
    assert A.mass_error_factor(0.0, 10.0) == 1.0
    at_limit = A.mass_error_factor(10.0, 10.0)
    assert A.MASS_FLOOR < at_limit < 0.25                                                      # a quarter of the range at most: never a veto
    assert A.mass_error_factor(50.0, 10.0) == pytest.approx(A.MASS_FLOOR) and A.mass_error_factor(50.0, 10.0) >= A.MASS_FLOOR


def test_the_factor_falls_with_the_size_of_the_error_and_ignores_its_sign():
    values = [A.mass_error_factor(x / 10, 10.0) for x in range(0, 101)]
    assert values == sorted(values, reverse=True) and values[0] > values[-1]
    assert all(A.MASS_FLOOR <= v <= 1.0 for v in values)
    assert A.mass_error_factor(-3.3, 10.0) == A.mass_error_factor(3.3, 10.0)


def test_typical_instrument_errors_are_hardly_penalised():
    assert A.mass_error_factor(1.3, 10.0) > 0.95 and A.mass_error_factor(3.2, 10.0) > 0.8        # CASMI 2022: median 1.3 ppm, 95 % under 3.2


def test_the_factor_scales_with_the_tolerance():
    for ppm in (0.5, 2.0, 4.0, 7.5):
        assert A.mass_error_factor(ppm, 10.0) == pytest.approx(A.mass_error_factor(ppm / 2, 5.0))
        assert A.mass_error_factor(ppm, 10.0) == pytest.approx(A.mass_error_factor(ppm * 2, 20.0))
    assert A.mass_error_factor(9.2, 20.0) > 2.5 * A.mass_error_factor(9.2, 10.0)                 # the same error is less serious with a wider window


def test_no_formula_or_no_tolerance_means_no_evidence():
    assert A.mass_error_factor(None, 10.0) == 1.0 and A.mass_error_factor(3.0, 0.0) == 1.0 and A.mass_error_factor(3.0, None) == 1.0


def test_the_gaussian_has_the_documented_width():
    assert A.mass_error_factor(10.0 / A.MASS_SIGMA_DIVISOR, 10.0) == pytest.approx(A.MASS_FLOOR + (1 - A.MASS_FLOOR) * math.exp(-0.5))


# ---------------------------------------------------------------- step 1: the adduct score
def fits(ppm_h, ppm_na):
    return {"[M+H]+": {"formula": "C14H16O3", "prob": 0.998, "ppm": ppm_h, "isotope_sim": None, "mass_factor": A.mass_error_factor(ppm_h, 10.0)},
            "[M+Na]+": {"formula": "C12H18O3", "prob": 0.736, "ppm": ppm_na, "isotope_sim": None, "mass_factor": A.mass_error_factor(ppm_na, 10.0)}}


def test_a_formula_near_the_tolerance_limit_loses_to_one_that_fits_well():
    peaks = RICH + [(233.117, 50.0), NA_ION]
    no_penalty = {a: {**f, "mass_factor": None} for a, f in fits(-1.1, 9.2).items()}
    before = A.score_adducts(233.117, True, peaks, buddy_fits=no_penalty)
    after = A.score_adducts(233.117, True, peaks, buddy_fits=fits(-1.1, 9.2))
    assert before.hypotheses[0].adduct == "[M+Na]+"                                              # the Na+ ion outvotes the prior without a penalty
    na, h = after.hypotheses[1], after.hypotheses[0]
    assert h.adduct == "[M+H]+" and h.score > 2 * na.score                                       # protonation leads clearly
    assert na.adduct == "[M+Na]+" and na.mass_factor < 0.3 and h.mass_factor > 0.95
    assert na.evidence["mass"]["ppm"] == 9.2 and na.evidence["mass"]["formula"] == "C12H18O3"
    assert na.evidence["mass"]["note"] == f"C12H18O3: mass error +9.2 ppm → ×{na.mass_factor:.2f}"


def test_the_same_mass_error_for_every_adduct_changes_nothing_between_them():
    peaks = RICH + [(233.117, 50.0)]
    same = A.score_adducts(233.117, True, peaks, buddy_fits=fits(4.0, 4.0))
    none = A.score_adducts(233.117, True, peaks, buddy_fits={a: {**f, "mass_factor": None} for a, f in fits(4.0, 4.0).items()})
    assert score_ratio(same) == pytest.approx(score_ratio(none), rel=0.01)


def score_ratio(result):
    by = {h.adduct: h.score for h in result.hypotheses}
    return by["[M+H]+"] / by["[M+Na]+"]


def test_an_adduct_without_a_mass_error_gets_factor_one_and_a_note():
    res = A.score_adducts(233.117, True, RICH, buddy_fits={"[M+H]+": {"formula": None, "prob": 0.0, "ppm": None, "isotope_sim": None, "note": "BUDDY found no formula within 10 ppm"}})
    h = next(h for h in res.hypotheses if h.adduct == "[M+H]+")
    assert h.mass_factor == 1.0 and "no formula" in h.evidence["mass"]["note"]
    assert next(h for h in res.hypotheses if h.adduct == "[M+K]+").evidence["mass"]["note"].startswith("not used")
    json.dumps(res.to_dict())


def test_the_best_formula_of_an_adduct_maximises_probability_times_mass_factor(monkeypatch):
    # rank 1 (prob 0.60) sits at +9.0 ppm; rank 2 (prob 0.30) at +0.5 ppm wins: 0.30 x 0.99 > 0.60 x 0.28
    table = {"[M+H]+": [("C10H20O5", 0.60, 9.0), ("C9H16N2O4", 0.30, 0.5), ("C11H24O4", 0.05, 0.2)]}
    adduct, result = run(monkeypatch, table, adduct=None)
    h = next(h for h in adduct.hypotheses if h.adduct == "[M+H]+")
    assert h.evidence["formula"]["formula"] == "C9H16N2O4" and h.evidence["formula"]["probability"] == 0.30       # BUDDY's own probability is kept
    assert h.evidence["mass"]["ppm"] == 0.5 and h.mass_factor == pytest.approx(A.mass_error_factor(0.5, 10.0), abs=1e-4)
    assert h.evidence["isotope"]["formula"] == "C9H16N2O4"


def test_a_confident_formula_at_the_limit_still_beats_an_improbable_good_one(monkeypatch):
    table = {"[M+H]+": [("C10H20O5", 0.99, 9.0), ("C9H16N2O4", 0.001, 0.0)]}
    adduct, _ = run(monkeypatch, table)
    assert next(h for h in adduct.hypotheses if h.adduct == "[M+H]+").evidence["formula"]["formula"] == "C10H20O5"


# ---------------------------------------------------------------- the whole of steps 1-2 on a challenge-16-like case
def test_a_challenge_16_like_case_goes_to_protonation(monkeypatch):
    table = {"[M+H]+": [("C14H16O3", 0.998, -1.1), ("C12H15F3O", 0.001, 9.4)], "[M+Na]+": [("C12H18O3", 0.736, 9.2), ("C10H15FN4", 0.033, -1.4)],
             "[M+H-H2O]+": [("C14H18O4", 0.998, -1.1)]}
    adduct, result = run(monkeypatch, table, extra=[NA_ION])
    assert adduct.retained_adducts[0] == "[M+H]+" and result.formulas[0].adduct == "[M+H]+" and result.formulas[0].formula == "C14H16O3"
    na = next(h for h in adduct.hypotheses if h.adduct == "[M+Na]+")
    assert na.mass_factor < 0.3 and na.score < 0.5 * adduct.hypotheses[0].score


def test_a_wider_tolerance_forgives_the_same_error(monkeypatch):
    table = {"[M+H]+": [("C14H16O3", 0.998, -1.1)], "[M+Na]+": [("C12H18O3", 0.736, 9.2)]}
    narrow, _ = run(monkeypatch, table, extra=[NA_ION], tol_ppm=10.0)
    wide, _ = run(monkeypatch, table, extra=[NA_ION], tol_ppm=30.0)
    f = lambda r: next(h for h in r.hypotheses if h.adduct == "[M+Na]+").mass_factor
    assert f(wide) > 2.5 * f(narrow)


# ---------------------------------------------------------------- step 2: the formula list
def test_each_adducts_formulas_are_ordered_by_probability_times_mass_factor_and_keep_their_buddy_rank(monkeypatch):
    table = {"[M+H]+": [("C10H20O5", 0.60, 9.0), ("C9H16N2O4", 0.30, 0.5), ("C11H24O4", 0.05, 0.2)]}
    _, result = run(monkeypatch, table)
    h = [f for f in result.formulas if f.adduct == "[M+H]+"]
    assert [f.formula for f in h] == ["C9H16N2O4", "C10H20O5", "C11H24O4"]                          # 0.30 x 1.0 > 0.60 x 0.28 > 0.05 x 1.0
    assert [f.rank for f in h] == [2, 1, 3] and [f.prob for f in h] == [0.30, 0.60, 0.05]           # BUDDY's rank and probability are untouched
    assert h[0].mass_factor == pytest.approx(A.mass_error_factor(0.5, 10.0)) and h[1].mass_factor < 0.3
    assert result.to_dict()["formulas"][0]["mass_factor"] == h[0].mass_factor


def test_equal_mass_factors_keep_buddys_order(monkeypatch):
    table = {"[M+H]+": [("C10H20O5", 0.5, 2.0), ("C9H16N2O4", 0.4, -2.0), ("C11H24O4", 0.3, 2.0)]}
    _, result = run(monkeypatch, table)
    assert [f.rank for f in result.formulas if f.adduct == "[M+H]+"] == [1, 2, 3]


def test_the_interleaving_takes_the_reordered_lists(monkeypatch):
    table = {"[M+H]+": [("C10H20O5", 0.60, 9.5), ("C9H16N2O4", 0.30, 0.5)], "[M+Na]+": [("C12H18O3", 0.95, 0.5), ("C8H18N4O2", 0.02, 0.5)]}
    monkeypatch.setattr(bd, "annotate_formula", FakeBuddy(table))

    def score_adducts(precursor_mz, positive, ms2, adduct=None, tol_da=0.02, buddy_fits=None):
        hyps = [A.AdductHypothesis("[M+H]+", 1.0, rank=1), A.AdductHypothesis("[M+Na]+", 0.8, rank=2)]
        return A.AdductResult(mode="auto", hypotheses=hyps)
    monkeypatch.setattr(A, "score_adducts", score_adducts)
    _, result = F.run_adduct_and_formula(233.117, True, RICH)
    assert [(f.adduct, f.formula) for f in result.formulas] == [("[M+H]+", "C9H16N2O4"), ("[M+Na]+", "C12H18O3"), ("[M+H]+", "C10H20O5"),
                                                                ("[M+Na]+", "C8H18N4O2")]          # the mass-fitting formula leads its adduct's list


def test_non_auto_mode_orders_by_the_factor_too_without_an_adduct_factor(monkeypatch):
    table = {"[M+H]+": [("C10H20O5", 0.60, 9.0), ("C9H16N2O4", 0.30, 0.5)]}
    adduct, result = run(monkeypatch, table, adduct="[M+H]+")
    assert adduct.mode == "user" and adduct.hypotheses[0].mass_factor is None
    assert [f.formula for f in result.formulas] == ["C9H16N2O4", "C10H20O5"]
