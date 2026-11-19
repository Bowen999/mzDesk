"""ESI MS/MS Annotation, steps 1-3: adduct hypotheses (adduct_step), BUDDY formulas (formula_step), library search and the Level-2 gate
(library_step).

Pure tests use hand-made spectra. BUDDY is replaced by a recorder (msbuddy_engine.annotate_formula is monkeypatched); the library tests
index a few hand-written records into a temporary MZDESK_LIBRARY_DIR and run the real library_search on them. Nothing here starts
MS-FINDER, FRIGID, MetFrag or CFM-ID.
Real runs (real BUDDY and the local MassBank / GNPS indexes; the first BUDDY call of a process takes 15-20 s) are opt-in:
MZDESK_RUN_ANNOTATION_STEPS=1. They skip when BUDDY or a library is not installed.
"""
import json
import os

import pytest

from mzdesk import adducts, spec
from mzdesk.annotation import adduct_step as A
from mzdesk.annotation import formula_step as F
from mzdesk.annotation import library_step as L
from mzdesk.annotation.levels import assign_level
from mzdesk.annotation.merge import formula_only_candidates, merge_candidates
from mzdesk.ms1 import msbuddy_engine as bd
from mzdesk.ms2 import library_index as lx
from mzdesk.ms2 import library_search as ls
from mzdesk.ms2 import library_setup as setup_mod
from mzdesk.shared.spectrum import parse_spectrum

RUN_REAL = os.environ.get("MZDESK_RUN_ANNOTATION_STEPS") == "1"
CAFFEINE_MS2 = parse_spectrum(spec.CAFFEINE_MASSBANK_MS2)          # MassBank MSBNK-LCSB-LU119903, 19 peaks, base peak 138.0663
CAFFEINE_MS1 = parse_spectrum(spec.CAFFEINE_MS1)                   # M, M+1, M+2 of C8H10N4O2 [M+H]+
CAFFEINE_MZ = 195.0877
CAFFEINE_SMILES = "CN1C=NC2=C1C(=O)N(C)C(=O)N2C"
CAFFEINE_KEY = "RYYVLZVUVIJVGH-UHFFFAOYSA-N"
RICH = [(80.0, 20), (100.0, 40), (120.0, 100), (140.0, 60), (160.0, 30), (180.0, 15)]       # 6 fragment peaks >= 5 %: not sparse
SPARSE = [(150.0, 100), (200.0, 30)]                                                      # 2 fragment peaks: sparse


def peaks_with(precursor, *extra, base=RICH):
    """A fragment spectrum for `precursor` (m/z) plus the precursor itself and the given extra (mz, intensity) peaks."""
    return list(base) + [(precursor, 50.0)] + list(extra)


def score_of(result, adduct):
    return next(h for h in result.hypotheses if h.adduct == adduct)


# ================================================================ step 1: adduct hypotheses
def test_auto_scores_every_adduct_of_the_polarity_in_the_master_list_and_the_scores_sum_to_one():
    for positive, listed in ((True, spec.ADDUCTS_POS), (False, spec.ADDUCTS_NEG)):
        res = A.score_adducts(300.0, positive, peaks_with(300.0))
        assert res.mode == "auto" and res.status == "done"
        assert sorted(h.adduct for h in res.hypotheses) == sorted(adducts.canonical(a) for a in listed if a != "Auto")
        assert sum(h.score for h in res.hypotheses) == pytest.approx(1.0, abs=1e-3)
        assert [h.rank for h in res.hypotheses] == list(range(1, len(res.hypotheses) + 1))
    scored = {h.adduct for positive in (True, False) for h in A.score_adducts(300.0, positive, RICH).hypotheses}
    assert set(adducts.CASMI_ADDUCTS) <= scored                                                 # the six CASMI adducts are all hypotheses


def test_every_master_list_adduct_has_a_library_count():
    assert all(a in A.LIBRARY_ADDUCT_COUNTS for a in A.master_adducts(True) + A.master_adducts(False))


def test_prior_is_the_flattened_library_frequency():
    shares = A.prior_shares(["[M+H]+", "[M+Na]+", "[M+NH4]+"])
    assert sum(shares.values()) == pytest.approx(1.0)
    raw = A.LIBRARY_ADDUCT_COUNTS["[M+H]+"] / A.LIBRARY_ADDUCT_COUNTS["[M+Na]+"]
    assert shares["[M+H]+"] / shares["[M+Na]+"] == pytest.approx(raw ** A.PRIOR_EXPONENT)       # weaker than the raw ratio
    assert 1 < shares["[M+H]+"] / shares["[M+Na]+"] < raw
    assert A.prior_shares(["[M+H]+", "[M+Br]-"])["[M+Br]-"] > 0                               # an adduct the table lacks counts as 1


def test_without_evidence_the_prior_puts_protonation_first_and_leaves_the_result_ambiguous():
    pos = A.score_adducts(300.0, True, peaks_with(300.0))
    assert [h.adduct for h in pos.retained] == ["[M+H]+", "[M+Na]+", "[M+NH4]+"] and pos.ambiguous
    assert pos.hypotheses[0].ambiguous and score_of(pos, "[M+Na]+").ambiguous and not score_of(pos, "[M+NH4]+").ambiguous
    assert score_of(pos, "[M+K]+").retained is False and "Ambiguous" in pos.note
    neg = A.score_adducts(300.0, False, peaks_with(300.0))
    assert neg.hypotheses[0].adduct == "[M-H]-" and "[M+HCOO]-" in neg.retained_adducts
    assert all(h.evidence["isotope"]["factor"] == 1.0 for h in pos.hypotheses)                 # no MS1: no isotope evidence


def test_loss_of_ammonia_supports_the_ammonium_adduct_but_does_not_decide():
    without = score_of(A.score_adducts(300.0, True, peaks_with(300.0)), "[M+NH4]+")
    res = A.score_adducts(300.0, True, peaks_with(300.0, (300.0 - 17.0265, 40.0)))
    nh4 = score_of(res, "[M+NH4]+")
    assert nh4.ms2_factor == 2.0 and nh4.score > without.score and res.hypotheses[0].adduct == "[M+H]+"
    found = nh4.evidence["ms2"]["signals"][0]
    assert found["found"] and found["rel_int"] == pytest.approx(40.0, abs=0.1) and "seen" in nh4.evidence["ms2"]["note"]       # 40 % of the base peak


def test_formate_loss_makes_the_formate_adduct_first_and_its_absence_counts_against_it():
    seen = A.score_adducts(400.0, False, peaks_with(400.0, (400.0 - 46.0055, 80.0)))
    assert seen.hypotheses[0].adduct == "[M+HCOO]-" and seen.hypotheses[0].ms2_factor == 10.0 and not seen.ambiguous
    ion = A.score_adducts(400.0, False, peaks_with(400.0, (44.9982, 10.0)))              # the HCOO- ion alone is enough
    assert score_of(ion, "[M+HCOO]-").ms2_factor == 10.0
    none = A.score_adducts(400.0, False, peaks_with(400.0))
    assert score_of(none, "[M+HCOO]-").ms2_factor == 0.5 and score_of(none, "[M+Cl]-").ms2_factor == 0.5
    assert none.hypotheses[0].adduct == "[M-H]-"


def test_chloride_and_acetate_signals():
    cl = A.score_adducts(400.0, False, peaks_with(400.0, (400.0 - 35.9767, 60.0)))
    assert cl.hypotheses[0].adduct == "[M+Cl]-" and score_of(cl, "[M+Cl]-").ms2_factor == 10.0
    assert score_of(A.score_adducts(400.0, False, peaks_with(400.0, (34.9689, 5.0))), "[M+Cl]-").ms2_factor == 10.0
    ac = A.score_adducts(400.0, False, peaks_with(400.0, (400.0 - 60.0211, 60.0)))
    assert score_of(ac, "[M+CH3COO]-").ms2_factor == 5.0
    assert score_of(A.score_adducts(400.0, False, peaks_with(400.0)), "[M+CH3COO]-").ms2_factor == 0.9


def test_alkali_adducts_get_their_ion_and_a_sparse_spectrum_helps_them():
    sparse = A.score_adducts(300.0, True, peaks_with(300.0, base=SPARSE))
    na = score_of(sparse, "[M+Na]+")
    assert na.ms2_factor == 1.5 and "few" in na.evidence["ms2"]["note"]                      # 2 fragment peaks
    rich = score_of(A.score_adducts(300.0, True, peaks_with(300.0)), "[M+Na]+")
    assert rich.ms2_factor == 0.7 and "many" in rich.evidence["ms2"]["note"]                  # 6 fragment peaks
    ion = score_of(A.score_adducts(300.0, True, peaks_with(300.0, (22.9892, 6.0))), "[M+Na]+")
    assert ion.ms2_factor == pytest.approx(5.0 * 0.7)
    k = score_of(A.score_adducts(300.0, True, peaks_with(300.0, (38.9632, 6.0))), "[M+K]+")
    assert k.ms2_factor == pytest.approx(5.0 * 0.7) and k.evidence["ms2"]["signals"][0]["label"].startswith("K+")


def test_signals_respect_the_fragment_tolerance_and_a_minimum_intensity():
    off = peaks_with(300.0, (300.0 - 17.0265 + 0.05, 40.0))                                   # 0.05 Da away
    assert score_of(A.score_adducts(300.0, True, off, tol_da=0.02), "[M+NH4]+").ms2_factor == 1.0
    assert score_of(A.score_adducts(300.0, True, off, tol_da=0.1), "[M+NH4]+").ms2_factor == 2.0
    weak = peaks_with(300.0, (300.0 - 17.0265, 1.0))                                          # 1 % of the base peak: below the 2 % rule
    assert score_of(A.score_adducts(300.0, True, weak), "[M+NH4]+").ms2_factor == 1.0
    assert score_of(A.score_adducts(300.0, True, weak[:-1] + [(300.0 - 17.0265, 2.5)]), "[M+NH4]+").ms2_factor == 2.0


def test_a_spectrum_without_fragments_cannot_show_the_absence_of_a_signal():
    only_precursor = [(400.0, 100.0), (401.003, 10.0)]
    res = A.score_adducts(400.0, False, only_precursor)
    assert all(h.ms2_factor == 1.0 for h in res.hypotheses)
    assert "no fragment peak" in score_of(res, "[M+HCOO]-").evidence["ms2"]["note"]
    assert all(h.ms2_factor == 1.0 for h in A.score_adducts(400.0, False, []).hypotheses)


def test_an_adduct_chosen_by_the_user_is_used_as_is():
    res = A.score_adducts(240.0, False, peaks_with(240.0), adduct="[M+FA-H]-")
    assert res.mode == "user" and res.retained_adducts == ["[M+HCOO]-"] and not res.ambiguous
    top = res.hypotheses[0]
    assert top.score == 1.0 and top.rank == 1 and top.retained and set(top.evidence) == {"user"}
    assert A.score_adducts(240.0, False, peaks_with(240.0), adduct="Auto").mode == "auto"
    assert A.score_adducts(240.0, False, peaks_with(240.0), adduct=None).mode == "auto"
    with pytest.raises(ValueError, match="does not match the negative ion mode"):
        A.score_adducts(240.0, False, peaks_with(240.0), adduct="[M+H]+")
    with pytest.raises(ValueError, match="Unknown adduct"):
        A.score_adducts(240.0, False, peaks_with(240.0), adduct="banana")


def test_isotope_factor_is_one_for_a_good_fit_and_falls_to_the_floor_for_a_poor_one():
    assert A.isotope_factor(None) == 1.0
    assert A.isotope_factor(1.0) == A.isotope_factor(A.ISOTOPE_GOOD) == 1.0
    assert A.isotope_factor(A.ISOTOPE_POOR) == A.isotope_factor(0.1) == A.ISOTOPE_FLOOR
    values = [A.isotope_factor(s / 100) for s in range(70, 101)]
    assert values == sorted(values) and 0.25 < A.isotope_factor(0.875) < 1.0


def test_the_isotope_fit_of_the_best_formula_reorders_the_hypotheses():
    base = A.score_adducts(300.0, True, peaks_with(300.0))
    fits = {"[M+H]+": {"formula": "C1", "prob": None, "isotope_sim": 0.80}, "[M+Na]+": {"formula": "C2", "prob": None, "isotope_sim": 0.99}}
    res = A.score_adducts(300.0, True, peaks_with(300.0), buddy_fits=fits)
    assert base.hypotheses[0].adduct == "[M+H]+" and res.hypotheses[0].adduct == "[M+Na]+"
    assert score_of(res, "[M+H]+").isotope_factor == A.ISOTOPE_FLOOR and score_of(res, "[M+Na]+").isotope_factor == 1.0
    iso = score_of(res, "[M+H]+").evidence["isotope"]
    assert iso["similarity"] == 0.80 and iso["formula"] == "C1" and "isotope similarity 0.800" in iso["note"]
    assert score_of(res, "[M+NH4]+").isotope_factor == 1.0                                 # not evaluated: no evidence, no penalty


def test_formula_factor_runs_from_the_floor_at_probability_zero_to_one():
    assert A.formula_factor(None) == 1.0 and A.formula_factor(1.0) == 1.0 and A.formula_factor(0.0) == A.FORMULA_FLOOR
    assert A.formula_factor(0.5) == pytest.approx(A.FORMULA_FLOOR + (1 - A.FORMULA_FLOOR) * 0.5) and A.formula_factor(7) == 1.0


def test_buddys_probability_of_the_best_formula_is_evidence_and_can_split_an_ambiguous_pair():
    base = A.score_adducts(300.0, True, peaks_with(300.0, base=SPARSE))                       # sparse: [M+Na]+ is 0.78 x [M+H]+ -> ambiguous
    assert base.ambiguous
    fits = {"[M+H]+": {"formula": "C10H20N2O5", "prob": 0.95, "ppm": 0.5}, "[M+Na]+": {"formula": None, "prob": 0.0, "ppm": None}}
    res = A.score_adducts(300.0, True, peaks_with(300.0, base=SPARSE), buddy_fits=fits)
    assert res.hypotheses[0].adduct == "[M+H]+" and not res.ambiguous and not res.hypotheses[0].ambiguous
    na = score_of(res, "[M+Na]+")
    assert na.formula_factor == A.FORMULA_FLOOR and na.evidence["formula"]["probability"] == 0.0
    assert score_of(res, "[M+H]+").evidence["formula"]["note"].startswith("BUDDY's best formula C10H20N2O5: probability 0.950")
    assert "not used" in score_of(res, "[M+NH4]+").evidence["formula"]["note"] and score_of(res, "[M+NH4]+").formula_factor == 1.0


def test_retention_keeps_at_most_three_hypotheses_that_reach_a_fifth_of_the_best(monkeypatch):
    monkeypatch.setattr(A, "PRIOR_EXPONENT", 1.0)
    monkeypatch.setitem(A.LIBRARY_ADDUCT_COUNTS, "[M+H]+", 1000)
    for name, n in (("[M+Na]+", 500), ("[M+K]+", 300), ("[M+NH4]+", 250), ("[M]+", 200), ("[M+H-H2O]+", 199)):
        monkeypatch.setitem(A.LIBRARY_ADDUCT_COUNTS, name, n)
    res = A.score_adducts(300.0, True, [], candidates=["[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M]+", "[M+H-H2O]+"])
    assert res.retained_adducts == ["[M+H]+", "[M+Na]+", "[M+K]+"]                           # [M+NH4]+ is 0.25 of the best but the 4th
    res = A.score_adducts(300.0, True, [], candidates=["[M+H]+", "[M]+", "[M+H-H2O]+"])
    assert res.retained_adducts == ["[M+H]+", "[M]+"]                                       # 0.20 of the best is kept, 0.199 is not


@pytest.mark.parametrize("count, ambiguous", [(299, True), (300, False), (301, False)])
def test_the_result_is_ambiguous_when_the_best_is_less_than_three_times_the_second(monkeypatch, count, ambiguous):
    monkeypatch.setattr(A, "PRIOR_EXPONENT", 1.0)
    monkeypatch.setitem(A.LIBRARY_ADDUCT_COUNTS, "[M+H]+", count)
    monkeypatch.setitem(A.LIBRARY_ADDUCT_COUNTS, "[M+Na]+", 100)
    res = A.score_adducts(300.0, True, [], candidates=["[M+H]+", "[M+Na]+"])
    assert res.ambiguous is ambiguous and res.hypotheses[0].ambiguous is ambiguous and res.hypotheses[1].ambiguous is ambiguous


def test_master_adducts_follow_the_master_list(monkeypatch):
    monkeypatch.setattr(spec, "ADDUCTS_POS", ["Auto", "[M+H]+", "[M+Na]+"])
    monkeypatch.setattr(spec, "ADDUCTS_NEG", ["Auto", "[M-H]-"])
    assert A.master_adducts(True) == ["[M+H]+", "[M+Na]+"] and A.master_adducts(False) == ["[M-H]-"]
    assert [h.adduct for h in A.score_adducts(300.0, True, RICH).hypotheses] == ["[M+H]+", "[M+Na]+"]


def test_adduct_result_is_json():
    res = A.score_adducts(300.0, True, peaks_with(300.0, (283.0, 30.0)))
    again = json.loads(json.dumps(res.to_dict()))
    assert again["retained_adducts"] == res.retained_adducts and again["hypotheses"][0]["evidence"]["prior"]["library_spectra"]


def test_clean_peaks_drops_invalid_values():
    assert A.clean_peaks([(2.0, 1.0), (1.0, 5.0), (3.0, 0.0), (float("nan"), 1.0), (4.0, -1.0)]) == [(1.0, 5.0), (2.0, 1.0)]
    assert A.clean_peaks(None) == [] and A.clean_peaks([]) == []


# ================================================================ step 2: BUDDY formulas (BUDDY replaced by a recorder)
class FakeBuddy:
    """Stands in for msbuddy_engine.annotate_formula. table = {adduct: [(formula, prob, ppm, isotope_sim), ...]}; calls records each run."""

    def __init__(self, table=None, error=None):
        self.table, self.error, self.calls = table or {}, error, []

    def __call__(self, precursor_mz, positive, adduct, rt, formula_input, ms2_peaks, ms1_peaks, ms_instr="orbitrap", halogen=False):
        self.calls.append({"adduct": adduct, "positive": positive, "ms2": len(ms2_peaks), "ms1": len(ms1_peaks), "preset": ms_instr, "halogen": halogen})
        if self.error:
            raise self.error
        rows = self.table.get(adduct, [])
        return {"candidates": [{"rank": i, "formula": f, "prob": p, "normed_prob": p, "fdr": 0.0, "isotope_sim": iso,
                                "mz_theo": precursor_mz - ppm * 1e-6 * precursor_mz, "ppm": ppm} for i, (f, p, ppm, iso) in enumerate(rows, start=1)],
                "adduct_used": adduct, "from_input": False, "ms1_used": len(ms1_peaks)}


@pytest.fixture
def buddy(monkeypatch):
    def install(table=None, error=None):
        fake = FakeBuddy(table, error)
        monkeypatch.setattr(bd, "annotate_formula", fake)
        return fake
    return install


MS1_PATTERN = [(300.0, 100.0), (301.003, 12.0), (302.006, 1.0)]
TABLE = {"[M+H]+": [("C10H20N2O5", 0.9, 0.5, 0.99), ("C9H16N4O4", 0.05, -3.0, 0.97)], "[M+Na]+": [("C11H16N3O3", 0.01, 4.0, 0.85)],
         "[M+NH4]+": [("C10H17NO5", 0.89, 0.5, 0.99)], "[M+K]+": [("C9H12N2O6", 0.01, 2.0, 0.99)]}


def test_buddy_runs_for_every_retained_adduct_but_only_the_best_gets_formulas_when_the_second_is_far_behind(buddy):
    fake = buddy(TABLE)
    adduct, formula = F.run_adduct_and_formula(300.0, True, peaks_with(300.0))             # rich spectrum: [M+Na]+ is far behind (prob 0.01)
    # BUDDY ran for [M+H]+, [M+Na]+ and [M+NH4]+ (retained before its evidence); [M+Na]+ then falls out of the cut (probability 0.01)
    assert [c["adduct"] for c in fake.calls] == ["[M+H]+", "[M+Na]+", "[M+NH4]+"] and adduct.retained_adducts == ["[M+H]+", "[M+NH4]+"]
    assert formula.adducts == ["[M+H]+"]
    assert [(f.formula, f.adduct, f.rank) for f in formula.formulas] == [("C10H20N2O5", "[M+H]+", 1), ("C9H16N4O4", "[M+H]+", 2)]
    assert formula.status == "done" and not formula.ms1_used
    assert score_of(adduct, "[M+Na]+").evidence["formula"]["probability"] == 0.01 and score_of(adduct, "[M+Na]+").formula_factor < 0.2
    assert "2 formulas for [M+H]+" in formula.note


def test_a_close_second_adduct_gets_formulas_too(buddy):
    table = {**TABLE, "[M+Na]+": [("C11H16N3O3", 0.9, 0.5, None)]}                             # BUDDY likes [M+Na]+ as much as [M+H]+ (same mass error: no mass penalty)
    fake = buddy(table)
    adduct, formula = F.run_adduct_and_formula(300.0, True, peaks_with(300.0, base=SPARSE))     # sparse: [M+Na]+ is 0.78 x [M+H]+
    assert {c["adduct"] for c in fake.calls} == set(adduct.retained_adducts) and formula.adducts == ["[M+H]+", "[M+Na]+"]
    assert [f.adduct for f in formula.formulas] == ["[M+H]+", "[M+Na]+", "[M+H]+"]            # interleaved by BUDDY rank, best adduct first
    assert formula.formula_strings == ["C10H20N2O5", "C11H16N3O3", "C9H16N4O4"]
    low = F.SECOND_ADDUCT_RATIO * adduct.hypotheses[0].score
    assert adduct.hypotheses[1].score >= low


def test_formula_order_preserves_adduct_priority_and_unrounded_buddy_ranks_after_filtering(buddy):
    buddy({"[M+H]+": [("C10H20N2O5", 0.900004, 0.5, None), ("C8H20O5", 0.900003, 12.0, None),
                       ("C9H16N4O4", 0.900002, -0.5, None)],
           "[M+Na]+": [("C11H16N3O3", 0.99, 0.2, None)]})
    adduct, result = F.run_adduct_and_formula(300.0, True, peaks_with(300.0, base=SPARSE), tol_ppm=10.0)
    assert result.adducts == ["[M+H]+", "[M+Na]+"] == adduct.retained_adducts[:2]
    assert [(f.adduct, f.formula, f.rank) for f in result.formulas] == [
        ("[M+H]+", "C10H20N2O5", 1), ("[M+Na]+", "C11H16N3O3", 1), ("[M+H]+", "C9H16N4O4", 3)]
    assert [r.formula for r in formula_only_candidates([f.to_dict() for f in result.formulas])] == [f.formula for f in result.formulas]


@pytest.mark.parametrize("second, include", [(0.24996, False), (0.25, True), (0.25004, True)])
def test_second_adduct_formula_cut_uses_unrounded_scores(buddy, monkeypatch, second, include):
    names = ["[M+H]+", "[M+Na]+", "[M+NH4]+", "[M+K]+"]
    shares = dict(zip(names, [0.5, second, 0.125, 0.375 - second]))
    monkeypatch.setattr(A, "master_adducts", lambda positive: names)
    monkeypatch.setattr(A, "prior_shares", lambda universe: shares)
    monkeypatch.setattr(A, "ms2_evidence", lambda *a, **k: {"factor": 1.0, "signals": []})
    buddy({name: [("C10H20N2O5", 0.89, 0.1, None)] for name in names})            # prob < EXTRA_ADDUCT_PROB: only the ratio rule decides
    adduct, result = F.run_adduct_and_formula(300.0, True, peaks_with(300.0))
    assert ("[M+Na]+" in result.adducts) is include
    assert adduct.hypotheses[1].score == pytest.approx(second, rel=0, abs=1e-10)


def test_adduct_rank_uses_full_precision_even_when_displayed_scores_tie(monkeypatch):
    names = ["[M+Na]+", "[M+H]+"]  # the higher score is second in input order
    monkeypatch.setattr(A, "prior_shares", lambda universe: {"[M+Na]+": 0.499996, "[M+H]+": 0.500004})
    monkeypatch.setattr(A, "ms2_evidence", lambda *a, **k: {"factor": 1.0, "signals": []})
    result = A.score_adducts(300.0, True, peaks_with(300.0), candidates=names)
    assert [h.adduct for h in result.hypotheses] == list(reversed(names))
    assert result.hypotheses[0].score > result.hypotheses[1].score
    assert [f"{h.score:.3f}" for h in result.hypotheses] == ["0.500", "0.500"]


def test_with_ms1_buddy_runs_for_every_retained_adduct_and_the_isotope_fit_feeds_the_ranking(buddy):
    fake = buddy(TABLE)
    adduct, formula = F.run_adduct_and_formula(300.0, True, peaks_with(300.0), MS1_PATTERN)
    assert [c["adduct"] for c in fake.calls] == ["[M+H]+", "[M+Na]+", "[M+NH4]+"] and all(c["ms1"] == 3 for c in fake.calls)
    na = score_of(adduct, "[M+Na]+")
    assert na.evidence["isotope"]["similarity"] == 0.85 and na.evidence["isotope"]["formula"] == "C11H16N3O3"
    assert na.isotope_factor == pytest.approx(A.isotope_factor(0.85)) and score_of(adduct, "[M+H]+").isotope_factor == 1.0
    assert formula.ms1_used and "MS1 isotope pattern" in formula.note
    assert formula.formulas[0].isotope_sim == 0.99 and formula.formulas[0].ppm == 0.5


def test_an_adduct_that_enters_the_cut_after_the_isotope_penalty_is_evaluated_too(buddy):
    poor = [("C10H20N2O5", 0.9, 0.5, 0.80)]                                                 # the isotope floor: x 0.25
    fake = buddy({"[M+H]+": poor, "[M+Na]+": poor, "[M+NH4]+": poor, "[M+K]+": poor, "[M+H-H2O]+": poor})
    adduct, _ = F.run_adduct_and_formula(300.0, True, peaks_with(300.0), MS1_PATTERN)
    called = [c["adduct"] for c in fake.calls]
    # [M+H]+, [M+Na]+ and [M+NH4]+ are retained first and all fit badly; the next one (not evaluated, x 1) climbs into the cut and is evaluated too
    assert called[:3] == ["[M+H]+", "[M+Na]+", "[M+NH4]+"] and len(called) > 3 and len(called) == len(set(called))
    supported = [a for a in adduct.retained_adducts if adducts.engine_name(a, "buddy")]
    assert set(supported) <= set(called) and all(score_of(adduct, a).evidence["isotope"]["similarity"] == 0.80 for a in supported)
    # [M]+ cannot be run in BUDDY, so nothing lowers it: with every other adduct penalised it enters the cut, unevaluated and never sent to BUDDY
    assert "[M]+" in adduct.retained_adducts and "[M]+" not in called and score_of(adduct, "[M]+").isotope_factor == 1.0


def test_buddy_gets_the_engine_spelling_of_the_adduct_and_the_result_the_canonical_one(buddy):
    fake = buddy({"[M+HCOO]-": [("C6H12O6", 0.9, 0.1, None)]})
    adduct, formula = F.run_adduct_and_formula(225.0612, False, peaks_with(225.0612), adduct="[M+FA-H]-")
    assert [c["adduct"] for c in fake.calls] == [adducts.engine_name("[M+FA-H]-", "buddy")] == ["[M+HCOO]-"]
    assert fake.calls[0]["positive"] is False
    assert formula.formulas[0].adduct == "[M+HCOO]-" and formula.adducts == ["[M+HCOO]-"] and adduct.mode == "user"


def test_formulas_that_miss_the_precursor_tolerance_are_dropped(buddy):
    buddy({"[M+H]+": [("C10H20N2O5", 0.9, 3.0, None), ("C9H16N4O4", 0.05, -12.0, None)]})
    _, formula = F.run_adduct_and_formula(300.0, True, peaks_with(300.0), tol_ppm=10.0)
    assert formula.formula_strings == ["C10H20N2O5"]
    _, strict = F.run_adduct_and_formula(300.0, True, peaks_with(300.0), tol_ppm=2.0)
    assert strict.formulas == [] and strict.status == "done" and "no formula within 2 ppm" in strict.note
    _, wide = F.run_adduct_and_formula(300.0, True, peaks_with(300.0), tol_ppm=20.0)
    assert wide.formula_strings == ["C10H20N2O5", "C9H16N4O4"]


def test_the_no_formula_note_names_the_tolerance_buddy_really_searched(buddy):
    fake = buddy({"[M+H]+": []})
    adduct, formula = F.run_adduct_and_formula(300.0, True, peaks_with(300.0), adduct="[M+H]+", tol_ppm=20.0)
    assert "no formula within 20 ppm" in formula.note and bd._tolerances(fake.calls[0]["preset"])[0] == 20.0


# Up to 10 ppm the tightest preset that covers the tolerance. Above 10 ppm no preset covers it (qtof is the widest, 10 ppm): BUDDY must then get the
# user's tolerance itself as (ms1_tol, ms2_tol) -- not fall back to qtof, which searched only +-10 ppm and silently lost formulas at 10-50 ppm.
# The MS2 tolerance stays at the widest preset's 20 ppm.
@pytest.mark.parametrize("tol_ppm, preset", [(1.0, "fticr"), (2.0, "fticr"), (3.0, "orbitrap"), (5.0, "orbitrap"),
                                             (8.0, "qtof"), (10.0, "qtof"), (10.5, (10.5, 20.0)), (20.0, (20.0, 20.0)), (50.0, (50.0, 20.0))])
def test_the_instrument_preset_follows_the_precursor_tolerance(buddy, tol_ppm, preset):
    fake = buddy(TABLE)
    F.run_adduct_and_formula(300.0, True, peaks_with(300.0), adduct="[M+H]+", tol_ppm=tol_ppm)
    assert fake.calls[0]["preset"] == preset


@pytest.mark.parametrize("tol_ppm", [0.5, 2.0, 3.0, 10.0, 10.5, 15.0, 30.0, 100.0])
def test_buddy_never_searches_a_narrower_window_than_the_precursor_tolerance(buddy, tol_ppm):
    """The tolerance in "BUDDY found no formula within X ppm" and in the table is the user's; BUDDY must have searched at least that."""
    fake = buddy(TABLE)
    F.run_adduct_and_formula(300.0, True, peaks_with(300.0), adduct="[M+H]+", tol_ppm=tol_ppm)
    assert bd._tolerances(fake.calls[0]["preset"])[0] >= tol_ppm


def test_halogen_option_is_passed_on(buddy):
    fake = buddy(TABLE)
    F.run_adduct_and_formula(300.0, True, peaks_with(300.0), adduct="[M+H]+", halogen=True)
    assert fake.calls[0]["halogen"] is True


def test_a_failing_buddy_gives_an_empty_list_and_a_note_and_is_not_tried_again(buddy):
    fake = buddy(error=RuntimeError("msbuddy exploded"))
    adduct, formula = F.run_adduct_and_formula(300.0, True, peaks_with(300.0, base=SPARSE), MS1_PATTERN)
    assert formula.status == "failed" and formula.formulas == [] and "msbuddy exploded" in formula.note
    assert "precursor mass" in formula.note and len(fake.calls) == 1                           # the other adducts were not tried
    assert adduct.status == "done" and adduct.retained_adducts                                 # the adduct step still ranks them
    assert all(h.isotope_factor == 1.0 for h in adduct.hypotheses)                              # and has no isotope evidence


def test_negative_mode_goes_to_buddy_as_negative_with_the_formate_name(buddy):
    fake = buddy({"[M+HCOO]-": [("C20H24O6", 0.99, 0.0, None)], "[M-H]-": [("C21H26O7", 0.2, 3.0, None)]})
    adduct, formula = F.run_adduct_and_formula(405.1555, False, peaks_with(405.1555, (405.1555 - 46.0055, 80.0)))
    assert adduct.hypotheses[0].adduct == "[M+HCOO]-" and fake.calls[0]["adduct"] == "[M+HCOO]-" and all(c["positive"] is False for c in fake.calls)
    assert [c["adduct"] for c in fake.calls] == ["[M+HCOO]-", "[M-H]-"] and set(adduct.retained_adducts) <= {c["adduct"] for c in fake.calls}
    assert [(f.formula, f.adduct) for f in formula.formulas] == [("C20H24O6", "[M+HCOO]-")]


def test_without_buddy_the_library_gate_falls_back_to_the_precursor_mass(buddy, tiny):
    buddy(error=RuntimeError("msbuddy exploded"))
    adduct, formula = F.run_adduct_and_formula(CAFFEINE_MZ, True, CAFFEINE_MS2)
    layer, summary = L.run_library(CAFFEINE_MS2, CAFFEINE_MZ, True, ["MassBank"], adducts=adduct.retained_adducts, formulas=formula.formulas)
    assert formula.status == "failed" and formula.formulas == [] and summary.passed
    assert summary.checks[4]["status"] == "not evaluated" and "precursor check stands in" in summary.checks[4]["detail"]
    assert summary.checks[2]["status"] == "pass" and layer.candidates[0].evidence["library"]["gate_passed"]


def test_a_missing_buddy_is_a_failure_with_a_note_too(buddy):
    buddy(error=ModuleNotFoundError("No module named 'msbuddy'"))
    _, formula = F.run_adduct_and_formula(300.0, True, peaks_with(300.0))
    assert formula.status == "failed" and "msbuddy" in formula.note and formula.formulas == []


def test_a_precursor_outside_the_buddy_range_is_skipped(buddy):
    buddy(error=ValueError("BUDDY supports precursor m/z between 1 and 1500."))
    _, formula = F.run_adduct_and_formula(1600.0, True, peaks_with(1600.0))
    assert formula.status == "skipped" and "between 1 and 1500" in formula.note and formula.formulas == []


def test_an_adduct_buddy_does_not_support_is_skipped_without_a_run(buddy):
    fake = buddy(TABLE)
    adduct, formula = F.run_adduct_and_formula(300.0, True, peaks_with(300.0), adduct="[M]+")
    assert fake.calls == [] and formula.status == "skipped" and "[M]+" in formula.note and adduct.retained_adducts == ["[M]+"]


def test_ms1_that_cannot_serve_as_a_pattern_is_ignored_with_a_note(buddy):
    fake = buddy(TABLE)
    single = F.run_adduct_and_formula(300.0, True, peaks_with(300.0), [(300.0, 100.0)])
    assert fake.calls[0]["ms1"] == 0 and not single[1].ms1_used and "single peak" in single[0].note
    far = F.run_adduct_and_formula(300.0, True, peaks_with(300.0), [(250.0, 100.0), (251.0, 10.0)])
    assert not far[1].ms1_used and "no peak within 10 ppm of the precursor" in far[0].note
    assert all(c["ms1"] == 0 for c in fake.calls) and len(fake.calls) % 2 == 0     # MS1 is not given to BUDDY
    assert F.ms1_usable(None, 300.0, 10.0) == (False, "") and F.ms1_usable(MS1_PATTERN, 300.0, 10.0) == (True, "")
    assert F.ms1_usable([(300.0004, 100.0), (301.0, 5.0)], 300.0, 1.0)[0] is False             # 1.3 ppm away with a 1 ppm tolerance


def test_the_explicit_adduct_needs_one_run_only_and_ignores_the_ms1_loop(buddy):
    fake = buddy(TABLE)
    F.run_adduct_and_formula(300.0, True, peaks_with(300.0), MS1_PATTERN, adduct="[M+Na]+")
    assert [c["adduct"] for c in fake.calls] == ["[M+Na]+"] and fake.calls[0]["ms1"] == 3


def test_the_formula_list_feeds_the_formula_only_rows_of_the_merge(buddy):
    buddy(TABLE)
    _, formula = F.run_adduct_and_formula(300.0, True, peaks_with(300.0), MS1_PATTERN, adduct="[M+H]+")
    rows = formula_only_candidates([f.to_dict() for f in formula.formulas])
    assert [r.formula for r in rows] == ["C10H20N2O5", "C9H16N4O4"] and rows[0].adduct == "[M+H]+"
    assert rows[0].layer_score == 0.9 and rows[0].mz_error_ppm == 0.5 and rows[0].evidence["formula_only"]["isotope_sim"] == 0.99


def test_step_results_are_json(buddy):
    buddy(TABLE)
    adduct, formula = F.run_adduct_and_formula(300.0, True, peaks_with(300.0), MS1_PATTERN)
    json.dumps(adduct.to_dict())
    again = json.loads(json.dumps(formula.to_dict()))
    assert again["formulas"][0]["formula"] == "C10H20N2O5" and again["adducts"] == ["[M+H]+"]


# ================================================================ step 3: library search and the gate
def hit(**kw):
    """A hit as library_search returns it (a caffeine [M+H]+ hit unless changed)."""
    base = {"rank": 1, "library": "MassBank", "row": 0, "accession": "ACC1", "name": "Caffeine", "adduct": "[M+H]+", "precursor_mz": 195.0877, "ppm": 0.0,
            "matches": 12, "scores": {L.PRIMARY_METRIC: 0.93}, "instrument": "QTOF", "collision_energy": "30", "inchikey": CAFFEINE_KEY,
            "smiles": CAFFEINE_SMILES}
    base.update(kw)
    return base


def gate(**kw):
    args = {"precursor_mz": CAFFEINE_MZ, "tol_ppm": 10.0, "n_query": 15, "adducts": ["[M+H]+"], "ms1": None,
            "formulas": [{"formula": "C8H10N4O2", "adduct": "[M+H]+"}]}
    args.update(kw)
    return L._Gate(**args)


def status(verdict):
    return {c["id"]: c["status"] for c in verdict["checks"]}


def test_a_hit_that_meets_every_requirement_passes_the_gate():
    v = gate().judge(hit())
    assert v["passed"] and status(v) == {"similarity": "pass", "peaks": "pass", "precursor": "pass", "adduct": "pass", "formula": "pass",
                                         "isotope": "not evaluated"}
    assert v["formula"] == "C8H10N4O2" and v["adduct"] == "[M+H]+" and v["error_from_formula"] and abs(v["error"]) < 1.0
    assert [c["label"] for c in v["checks"]] == ["Similarity", "Matched peaks", "Precursor m/z", "Adduct", "Formula", "Isotope pattern"]


@pytest.mark.parametrize("similarity, ok", [(0.80, True), (0.7999, False), (0.95, True), (0.72, False)])
def test_gate_similarity_boundary(similarity, ok):
    v = gate().judge(hit(scores={L.PRIMARY_METRIC: similarity}))
    assert (status(v)["similarity"] == "pass") is ok and v["passed"] is ok


@pytest.mark.parametrize("matches, n_query, ok", [(4, 15, True), (3, 15, False), (19, 19, True), (3, 3, True), (2, 3, False), (1, 1, True),
                                                  (4, 4, True), (3, 4, False)])
def test_gate_needs_four_matched_peaks_or_all_of_them_when_the_query_has_fewer(matches, n_query, ok):
    g = gate(n_query=n_query)
    assert g.min_matched == min(4, n_query)
    assert (status(g.judge(hit(matches=matches)))["peaks"] == "pass") is ok


def test_gate_precursor_error_is_measured_against_the_hit_formula_and_adduct():
    theo = bd.ion_mz("C8H10N4O2", "[M+H]+")                          # the query is 12 ppm off; the library's own recorded m/z (ppm=0) agrees with it
    v = gate(precursor_mz=theo * (1 + 12e-6)).judge(hit(ppm=0.0))
    assert status(v)["precursor"] == "fail" and v["error"] == pytest.approx(12.0, abs=1e-6) and not v["passed"]
    assert "against C8H10N4O2 [M+H]+" in v["checks"][2]["detail"]
    assert status(gate(precursor_mz=theo * (1 + 9.9e-6)).judge(hit()))["precursor"] == "pass"
    assert status(gate(precursor_mz=theo * (1 + 10.5e-6), tol_ppm=11.0).judge(hit()))["precursor"] == "pass"
    assert status(gate(precursor_mz=theo * (1 - 10.5e-6)).judge(hit()))["precursor"] == "fail"          # both signs count


def test_gate_precursor_error_of_a_hit_without_structure_is_the_library_one():
    v = gate(formulas=[]).judge(hit(smiles=None, ppm=-4.2))
    assert v["error"] == -4.2 and not v["error_from_formula"] and status(v)["precursor"] == "pass"
    assert "library's recorded precursor" in v["checks"][2]["detail"]
    assert status(gate(formulas=[]).judge(hit(smiles=None, ppm=11.0)))["precursor"] == "fail"


def test_gate_adduct_is_a_retained_hypothesis_in_any_spelling():
    g = gate(adducts=["[M+H]+", "[M+Na]+"])
    assert status(g.judge(hit(adduct="[M+Na]+")))["adduct"] == "pass"
    v = g.judge(hit(adduct="[M+NH4]+"))
    assert status(v)["adduct"] == "fail" and "[M+NH4]+ is not one of the retained adducts ([M+H]+, [M+Na]+)" in v["checks"][3]["detail"]
    assert v["adduct_source"] == "step1" and g.judge(hit())["adduct_source"] == "step1"
    assert status(g.judge(hit(adduct="")))["adduct"] == "fail" and g.judge(hit(adduct=None))["checks"][3]["detail"] == "the library entry records no adduct"
    assert "'M-' is not a known adduct" in g.judge(hit(adduct="M-"))["checks"][3]["detail"]
    f = gate(adducts=["[M+HCOO]-"], formulas=[])                                                 # any spelling of formate
    assert status(f.judge(hit(adduct="[M+FA-H]-", smiles="OCC(O)C(O)C(O)C(O)CO")))["adduct"] == "pass"
    assert gate(adducts=[A.AdductHypothesis("[M+H]+", 1.0)]).adducts == ["[M+H]+"]                  # hypothesis objects are accepted


def test_a_hit_with_an_adduct_step_1_missed_passes_when_the_library_match_confirms_it():
    theo = bd.ion_mz("C8H10N4O2", "[M+Na]+")                                         # the query is caffeine [M+Na]+, step 1 kept [M+H]+ only
    v = gate(precursor_mz=theo, adducts=["[M+H]+"], formulas=[]).judge(hit(adduct="[M+Na]+"))
    assert v["passed"] and v["adduct_source"] == "library" and v["adduct"] == "[M+Na]+"
    check = v["checks"][3]
    assert check["status"] == "pass" and "not one of the retained adducts ([M+H]+)" in check["detail"]
    assert "known positive-mode adduct and the library match confirms it" in check["detail"] and "against C8H10N4O2 [M+Na]+" in check["detail"]
    assert gate(adducts=["[M+Na]+"]).judge(hit())["adduct_source"] == "library"       # caffeine [M+H]+ hit, step 1 kept [M+Na]+ only
    assert gate().judge(hit())["adduct_source"] == "step1"                            # a retained adduct needs no confirmation
    neg = bd.ion_mz("C8H10N4O2", "[M+Cl]-")
    v = gate(precursor_mz=neg, adducts=["[M-H]-"], formulas=[], positive=False).judge(hit(adduct="[M+Cl]-"))
    assert v["passed"] and v["adduct_source"] == "library" and "negative-mode adduct" in v["checks"][3]["detail"]


def test_the_library_match_confirms_an_adduct_only_when_everything_else_holds():
    theo = bd.ion_mz("C8H10N4O2", "[M+Na]+")
    base = {"precursor_mz": theo, "adducts": ["[M+H]+"], "formulas": []}

    def failing(g, **kw):
        v = g.judge(hit(adduct="[M+Na]+", **kw))
        assert not v["passed"] and status(v)["adduct"] == "fail" and v["adduct_source"] == "step1"
        return status(v)

    assert failing(gate(**base), scores={L.PRIMARY_METRIC: 0.75})["similarity"] == "fail"                # a weak match confirms nothing
    assert failing(gate(**base), matches=3)["peaks"] == "fail"
    assert failing(gate(**{**base, "precursor_mz": CAFFEINE_MZ}))["precursor"] == "fail"                # 195.0877 is [M+H]+, not [M+Na]+, of caffeine
    only_library_mz = failing(gate(**base), smiles=None, ppm=0.0)                                       # no structure: no error under its own adduct
    assert only_library_mz["precursor"] == "pass"
    listed = {**base, "formulas": [{"formula": "C7H14O6", "adduct": "[M+Na]+"}], "precursor_mz": theo + 0.01}       # BUDDY disagrees and the mass does not fit
    assert failing(gate(**listed))["formula"] == "fail"
    assert failing(gate(**base, ms1=[(theo, 100.0), (theo + 1.003, 80.0), (theo + 2.006, 60.0)]))["isotope"] == "fail"
    assert not gate(**base, positive=True).judge(hit(adduct="[M-H]-"))["passed"]                         # the other polarity
    assert gate(**base).judge(hit(adduct="[M-H]-"))["checks"][3]["status"] == "fail"
    assert gate(**base).judge(hit(adduct="M-"))["checks"][3]["status"] == "fail"                         # not a known adduct


def test_gate_formula_must_be_one_of_the_formulas_buddy_gave_for_the_hit_adduct():
    g = gate(formulas=[F.FormulaCandidate("C8H10N4O2", "[M+H]+", 1), F.FormulaCandidate("C7H14O6", "[M+H]+", 2),
                       F.FormulaCandidate("C11H12N2", "[M+Na]+", 1)])
    assert status(g.judge(hit()))["formula"] == "pass"
    other = g.judge(hit(smiles="OCC(O)C(O)C(O)C(O)CO"))                                        # C6H14O6 is in no list
    assert status(other)["formula"] == "fail" and "C6H14O6 is not among the BUDDY formulas for [M+H]+ (C8H10N4O2, C7H14O6)" in other["checks"][4]["detail"]
    assert status(g.judge(hit(smiles=None)))["formula"] == "fail"                                # a list exists and the hit has no formula
    assert "unknown" in g.judge(hit(smiles=None))["checks"][4]["detail"]


def test_gate_formula_check_is_not_evaluated_when_buddy_gave_no_formula_for_the_hit_adduct():
    g = gate(adducts=["[M+H]+", "[M+Na]+"], formulas=[{"formula": "C8H10N4O2", "adduct": "[M+H]+"}])
    v = g.judge(hit(adduct="[M+Na]+"))
    assert status(v)["formula"] == "not evaluated" and "BUDDY gave no formula for [M+Na]+" in v["checks"][4]["detail"]
    assert gate(formulas=[]).judge(hit(smiles=None))["checks"][4]["status"] == "not evaluated"      # no list at all: the mass check stands in
    assert status(gate(formulas=["C8H10N4O2"]).judge(hit(adduct="[M+H]+")))["formula"] == "pass"   # a plain string counts for every adduct
    assert status(gate(formulas=["C7H14O6"]).judge(hit()))["formula"] == "not evaluated"             # not listed, but the mass fits (tests/test_library_gate.py)
    assert status(gate(formulas=["C7H14O6"], precursor_mz=CAFFEINE_MZ + 0.01).judge(hit()))["formula"] == "fail"      # not listed and the mass does not fit


def test_gate_formula_check_leaves_what_buddy_cannot_produce_to_the_mass_check():
    diazepam = "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"                                           # C16H13ClN2O: BUDDY searches without Cl by default
    mz = bd.ion_mz("C16H13ClN2O", "[M+H]+")
    v = gate(precursor_mz=mz).judge(hit(smiles=diazepam))
    assert status(v)["formula"] == "not evaluated" and "outside BUDDY's search space" in v["checks"][4]["detail"] and v["passed"]
    v = gate(precursor_mz=mz, halogen=True).judge(hit(smiles=diazepam))                          # with halogens on, BUDDY could have listed it
    assert status(v)["formula"] == "not evaluated" and "not among BUDDY's top 1 formulas" in v["checks"][4]["detail"] and v["passed"]   # listed formulas miss it, its mass fits
    v = gate(precursor_mz=mz + 0.01, halogen=True).judge(hit(smiles=diazepam))
    assert status(v)["formula"] == "fail" and not v["passed"]
    assert F.buddy_covers("C8H10N4O2") and not F.buddy_covers("C16H13ClN2O") and F.buddy_covers("C16H13ClN2O", halogen=True)
    assert not F.buddy_covers("C81H100") and not F.buddy_covers("C4H12Si", halogen=True) and not F.buddy_covers("") and F.buddy_covers("C80H150N20O30P10S15")


def test_gate_formula_comparison_ignores_the_element_order():
    assert L.same_formula("C8H10N4O2", "N4C8H10O2") and not L.same_formula("C8H10N4O2", "C8H10N4O3") and not L.same_formula(None, "C8H10N4O2")
    assert L.same_formula("C6H12O6", "C6H12O6") and not L.same_formula("C6H12O6", "C6H12O")


def test_gate_isotope_check_is_not_evaluated_without_ms1_and_never_blocks():
    v = gate(ms1=None).judge(hit())
    assert v["isotope"] == "not evaluated" and v["isotope_similarity"] is None and v["passed"]
    assert "no usable MS1" in v["checks"][5]["detail"]


def test_gate_isotope_check_compares_the_pattern_of_the_hit_formula_with_ms1():
    good = gate(ms1=CAFFEINE_MS1).judge(hit())
    assert good["isotope"] == "pass" and good["isotope_similarity"] > 0.99 and good["passed"]
    heavy = [(195.0877, 100.0), (196.0910, 80.0), (197.0940, 60.0)]                           # M+1 and M+2 far above anything caffeine has
    bad = gate(ms1=heavy).judge(hit())
    assert bad["isotope"] == "fail" and bad["isotope_similarity"] < L.GATE_ISOTOPE and not bad["passed"]
    chlorine = [(195.0877, 100.0), (196.0910, 10.0), (197.0848, 33.0)]                       # a Cl-like M+2: the cosine does not mind (known limit)
    assert gate(ms1=chlorine).judge(hit())["isotope"] == "pass"
    nothing = gate(ms1=[(250.0, 100.0), (251.0, 10.0)]).judge(hit())                           # no peak where the hit formula puts M
    assert nothing["isotope"] == "fail" and nothing["isotope_similarity"] < L.GATE_ISOTOPE and not nothing["passed"]
    assert status(gate(ms1=CAFFEINE_MS1).judge(hit(smiles=None)))["isotope"] == "not evaluated"


def test_hit_formula_comes_from_the_neutral_largest_fragment():
    assert L.hit_formula(CAFFEINE_SMILES) == "C8H10N4O2"
    assert L.hit_formula("C[NH3+]") == "CH5N" and L.hit_formula("CC(=O)[O-]") == "C2H4O2" and L.hit_formula("[Na+].CC(=O)[O-]") == "C2H4O2"
    assert L.hit_formula("C[N+](C)(C)CCO") == "C5H14NO"                                        # a permanent cation keeps its formula
    assert L.hit_formula("[NH3+]CC([O-])=O") == "C2H5NO2"
    assert L.hit_formula("") is None and L.hit_formula(None) is None and L.hit_formula("not a smiles") is None and L.hit_formula("N/A") is None
    assert "[2H]" in L.hit_formula("[2H]C([2H])([2H])O") and L._ion_mz(L.hit_formula("[2H]C([2H])([2H])O"), "[M+H]+") is None     # labels: no m/z, no BUDDY


# ---------------------------------------------------------------- the search on a tiny index
def record(name, accession, smiles, adduct, precursor, peaks, inchikey="", mode="POSITIVE"):
    lines = [f"Name: {name}", f"DB#: {accession}"] + ([f"InChIKey: {inchikey}"] if inchikey else []) + ([f"SMILES: {smiles}"] if smiles else [])
    lines += ([f"Precursor_type: {adduct}"] if adduct else [])
    lines += ["Spectrum_type: MS2", f"PrecursorMZ: {precursor}", "Instrument_type: LC-ESI-QTOF", "Instrument: fixture", f"Ion_mode: {mode}",
              "Collision_energy: 30", f"Num Peaks: {len(peaks)}"]
    return "\n".join(lines + [f"{mz} {inten}" for mz, inten in peaks]) + "\n\n"


GLUCOSE = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"
GLUCOSE_NA = [(203.0526, 999), (185.0420, 120), (143.0314, 60), (113.0233, 40), (85.0284, 20)]
GLUCOSE_FA = [(225.0612, 80), (179.0561, 999), (89.0244, 300), (59.0139, 100)]
CHOLINE_PEAKS = [(104.1070, 999), (60.0808, 400), (58.0651, 100), (45.0335, 50)]
THEOBROMINE = "Cn1cnc2c1c(=O)[nH]c(=O)n2C"
THEOBROMINE_PEAKS = [(181.0720, 999), (138.0662, 300), (110.0349, 100), (69.0447, 50), (42.0338, 20)]
UNKNOWN_PEAKS = [(300.1, 999), (200.05, 400), (150.02, 300), (100.01, 100)]
UNKNOWN_PEAKS_2 = [(300.1, 900), (200.05, 450), (150.02, 250), (100.01, 120)]            # nearly the same spectrum as UNKNOWN_PEAKS
CAFFEINE_D3_PEAKS = [(198.1065, 100), (141.0851, 999), (113.0901, 40), (86.0, 20)]
CAFFEINE_B = [(69.0447, 20), (83.0604, 30), (110.0713, 400), (138.0663, 999), (195.0878, 300)]      # a second caffeine spectrum
RECORDS = "".join([
    record("Caffeine", "FIX-CAF1", CAFFEINE_SMILES, "[M+H]+", 195.0877, CAFFEINE_MS2, CAFFEINE_KEY),
    record("Caffeine", "FIX-CAF2", CAFFEINE_SMILES, "[M+H]+", 195.0877, CAFFEINE_B, CAFFEINE_KEY),
    record("Caffeine without adduct", "FIX-CAF3", CAFFEINE_SMILES, "", 195.0877, CAFFEINE_B, CAFFEINE_KEY),
    record("Structure-less caffeine", "FIX-NOSMI", "", "[M+H]+", 195.0877, CAFFEINE_B),
    record("Theobromine without adduct", "FIX-THB", THEOBROMINE, "", 181.0720, THEOBROMINE_PEAKS, "YAPQBXQYLJRXSA-UHFFFAOYSA-N"),
    record("Unknown with N/A structure", "FIX-NA1", "N/A", "[M-H]-", 300.1000, UNKNOWN_PEAKS, "N/A", "NEGATIVE"),
    record("Another unknown with N/A structure", "FIX-NA2", "N/A", "[M-H]-", 300.1000, UNKNOWN_PEAKS_2, "N/A", "NEGATIVE"),
    record("Caffeine-d3", "FIX-CAFD3", "[2H]C([2H])([2H])N1C=NC2=C1C(=O)N(C)C(=O)N2C", "[M+H]+", 198.1065, CAFFEINE_D3_PEAKS),
    record("Glucose sodium adduct", "FIX-GLCNA", GLUCOSE, "[M+Na]+", 203.0526, GLUCOSE_NA, "WQZGKKKJIJFFOK-GASJEMHNSA-N"),
    record("Glucose formate adduct", "FIX-GLCFA", GLUCOSE, "[M+FA-H]-", 225.0612, GLUCOSE_FA, "WQZGKKKJIJFFOK-GASJEMHNSA-N", "NEGATIVE"),
    record("Choline", "FIX-CHO", "C[N+](C)(C)CCO", "[M]+", 104.1070, CHOLINE_PEAKS, "OEYIOHPDSNJKLS-UHFFFAOYSA-N"),
])


@pytest.fixture(scope="module")
def tiny(tmp_path_factory):
    """A MassBank index of RECORDS in a temporary library directory (GNPS is not installed there)."""
    root = tmp_path_factory.mktemp("steps_libs")
    mp = pytest.MonkeyPatch()
    mp.setenv("MZDESK_LIBRARY_DIR", str(root))
    lx.clear_cache()
    source = root / "fixture.msp"
    source.write_text(RECORDS)
    setup_mod.build_from("MassBank", source, {"url": "file://fixture", "version": "fixture-1", "date": "2026-01-01"}, progress=lambda *_: None)
    yield root
    mp.undo()
    lx.clear_cache()


@pytest.fixture
def no_libs(monkeypatch, tmp_path):
    monkeypatch.setenv("MZDESK_LIBRARY_DIR", str(tmp_path / "nothing"))
    lx.clear_cache()
    yield
    lx.clear_cache()


CAFFEINE_FORMULAS = [F.FormulaCandidate("C8H10N4O2", "[M+H]+", 1, 0.99, 0.2)]


def run(peaks=CAFFEINE_MS2, mz=CAFFEINE_MZ, positive=True, libraries=("MassBank",), adducts_=("[M+H]+",), formulas=CAFFEINE_FORMULAS, **kw):
    return L.run_library(peaks, mz, positive, list(libraries), adducts=list(adducts_), formulas=formulas, **kw)


def test_a_library_spectrum_of_the_compound_passes_the_gate_and_is_level_2a(tiny):
    layer, summary = run()
    assert layer.layer == "library" and layer.status == "done" and layer.seconds > 0
    assert summary.passed and summary.best_id == "lib1" == summary.top_id and layer.candidates[0].id == "lib1"
    top = layer.candidates[0]
    ev = top.evidence["library"]
    assert (top.name, top.formula, top.adduct, top.smiles) == ("Caffeine", "C8H10N4O2", "[M+H]+", CAFFEINE_SMILES)
    assert (top.inchikey, top.ik14) == (CAFFEINE_KEY, "RYYVLZVUVIJVGH")
    assert ev["gate_passed"] is True and ev["failed_checks"] == [] and ev["metric"] == "SpectralEntropy (weighted)"
    assert ev["similarity"] == pytest.approx(1.0) and ev["matched_peaks"] == 19 and ev["isotope"] == "not evaluated" and ev["isotope_similarity"] is None
    assert top.layer_score == top.layer_score_norm == pytest.approx(ev["similarity"]) and top.hypothesis is False
    assert abs(top.mz_error_ppm) < 1.0 and ev["precursor_error_ppm"] == top.mz_error_ppm and "C8H10N4O2 [M+H]+" in ev["precursor_error_basis"]
    assert (ev["library"], ev["accession"], ev["library_adduct"], ev["n_spectra"]) == ("MassBank", "FIX-CAF1", "[M+H]+", 3)
    assert "Gate passed by Caffeine" in layer.note
    assert [c["status"] for c in summary.checks] == ["pass"] * 5 + ["not evaluated"]
    assert summary.thresholds == {"metric": "SpectralEntropy (weighted)", "metrics": ["SpectralEntropy (weighted)"], "similarity": 0.8, "min_matched_peaks": 4, "precursor_ppm": 10.0,
                                  "isotope": 0.8, "list": 0.7}
    merged = merge_candidates([layer])
    assert merged[0].level == "2a" and merged[0].name == "Caffeine"


def caffeine(layer):
    return next(c for c in layer.candidates if c.name == "Caffeine")


def failed_labels(candidate):
    return [f.split(":")[0] for f in candidate.evidence["library"]["failed_checks"]]


def test_one_candidate_per_compound_passing_ones_first_then_by_similarity_and_a_known_structure(tiny):
    layer, summary = run(peaks=CAFFEINE_B, formulas=[])
    names = [c.name for c in layer.candidates]
    assert sorted(names) == ["Caffeine", "Structure-less caffeine"]                           # three caffeine spectra, one candidate
    assert caffeine(layer).evidence["library"]["n_spectra"] == 3
    assert layer.candidates[0].layer_score == pytest.approx(layer.candidates[1].layer_score) and names[0] == "Caffeine"     # tie: the known structure first
    assert summary.best_id == layer.candidates[0].id
    layer, _ = run(peaks=CAFFEINE_B)                                                          # with a formula list the structure-less hit fails
    flags = [c.evidence["library"]["gate_passed"] for c in layer.candidates]
    assert flags == [True, False] and layer.candidates[1].name == "Structure-less caffeine"
    sims = [c.layer_score_norm for c in run(peaks=[(a, i) for a, i in CAFFEINE_B if a > 100], adducts_=("[M+Na]+",))[0].candidates]
    assert sims == sorted(sims, reverse=True)                                                  # all failing: by similarity


def test_a_perturbed_spectrum_is_listed_with_its_failed_checks(tiny):
    query = [p for p in CAFFEINE_MS2 if abs(p[0] - 138.0663) > 0.01]                           # the base peak removed: similarity about 0.76
    layer, summary = run(peaks=query)
    top = layer.candidates[0]
    ev = top.evidence["library"]
    assert not summary.passed and summary.best_id is None and summary.top_id == top.id
    assert 0.7 <= ev["similarity"] < 0.8 and ev["gate_passed"] is False and len(ev["failed_checks"]) == 1
    assert ev["failed_checks"][0].startswith("Similarity: SpectralEntropy (weighted) 0.") and ev["failed_checks"][0].endswith("< 0.80")
    assert [c["status"] for c in summary.checks][:2] == ["fail", "pass"] and "Gate not passed" in layer.note
    merged = merge_candidates([layer])[0]
    assert merged.level == "3" and "failed: Similarity" in merged.level_reason


def test_a_spectrum_that_matches_nothing_gives_no_candidate(tiny):
    layer, summary = run(peaks=[(a + 0.05, i) for a, i in CAFFEINE_MS2])                        # every peak 0.05 Da off: no match
    assert layer.status == "done" and layer.candidates == [] and not summary.passed and summary.top_id is None and summary.checks == []
    assert "none scored 0.70 or higher" in layer.note
    layer, _ = run(mz=250.0)
    assert layer.candidates == [] and "No library spectrum within 10 ppm of m/z 250.0000" in layer.note


def test_each_requirement_can_fail_on_its_own(tiny):
    layer, summary = run(formulas=[F.FormulaCandidate("C7H14O6", "[M+H]+", 1)])      # not in BUDDY's list, but its mass fits the precursor: not a failure
    assert summary.passed and failed_labels(caffeine(layer)) == [] and summary.checks[4]["status"] == "not evaluated"
    layer, _ = run(formulas=[F.FormulaCandidate("C7H14O6", "[M+H]+", 1)], tol_ppm=0.1)    # ... unless the mass is off too
    assert failed_labels(caffeine(layer)) == ["Precursor m/z", "Formula"] and "C7H14O6" in caffeine(layer).evidence["library"]["failed_checks"][1]
    layer, _ = run(tol_ppm=0.1)                              # the record is inside the 0.1 ppm window, but caffeine [M+H]+ is 0.25 ppm from the query
    assert failed_labels(caffeine(layer)) == ["Precursor m/z"] and caffeine(layer).mz_error_ppm == pytest.approx(0.2459, abs=1e-3)
    layer, _ = run(peaks=[p for p in CAFFEINE_MS2 if abs(p[0] - 138.0663) > 0.01])
    assert failed_labels(caffeine(layer)) == ["Similarity"]
    layer, summary = run(ms1_peaks=[(195.0877, 100.0), (196.0910, 80.0), (197.0940, 60.0)])
    assert failed_labels(caffeine(layer)) == ["Isotope pattern"] and "isotope similarity 0." in caffeine(layer).evidence["library"]["failed_checks"][0]
    assert caffeine(layer).evidence["library"]["isotope"] == "fail" and not summary.passed


def test_a_library_match_confirms_an_adduct_step_1_missed(tiny):
    layer, summary = run(adducts_=("[M+Na]+",))                                             # step 1 kept [M+Na]+; the caffeine spectra say [M+H]+
    top = caffeine(layer)
    ev = top.evidence["library"]
    assert summary.passed and ev["gate_passed"] and ev["adduct_source"] == "library" and top.adduct == "[M+H]+" and ev["failed_checks"] == []
    assert summary.checks[3]["status"] == "pass" and "the library match confirms it" in summary.checks[3]["detail"]
    assert merge_candidates([layer])[0].level == "2a" and merge_candidates([layer])[0].adduct == "[M+H]+"
    assert caffeine(run()[0]).evidence["library"]["adduct_source"] == "step1"
    layer, summary = run(adducts_=("[M+Na]+",), peaks=[p for p in CAFFEINE_MS2 if abs(p[0] - 138.0663) > 0.01])        # ... but not if another check fails
    assert not summary.passed and failed_labels(caffeine(layer)) == ["Similarity", "Adduct"]
    assert caffeine(layer).evidence["library"]["adduct_source"] == "step1"


def test_the_peak_requirement_is_four_or_all_of_a_shorter_query(tiny):
    three = [(138.0663, 999), (110.0713, 101), (195.0878, 847)]                                # all 3 peaks match the library spectrum
    layer, summary = run(peaks=three)
    assert summary.thresholds["min_matched_peaks"] == 3 and summary.passed and caffeine(layer).evidence["library"]["matched_peaks"] == 3
    layer, summary = run(peaks=[(138.0663, 999), (195.0878, 847), (60.0, 1)])                  # 2 of 3 match: not all of them
    assert failed_labels(caffeine(layer)) == ["Matched peaks"] and caffeine(layer).evidence["library"]["matched_peaks"] == 2
    layer, summary = run(peaks=three + [(60.0, 1)])                                           # 3 of 4 match: fewer than 4
    assert summary.thresholds["min_matched_peaks"] == 4 and failed_labels(caffeine(layer)) == ["Matched peaks"]
    assert summary.checks[1]["status"] == "fail" and summary.checks[1]["detail"] == "3 matched, 4 needed"


def test_a_chosen_minimum_of_matched_peaks_replaces_the_four(tiny):
    three_of_four = [(138.0663, 999), (110.0713, 101), (195.0878, 847), (60.0, 1)]
    layer, summary = run(peaks=three_of_four, min_matched_peaks=3)                             # 3 matched peaks are enough now
    assert summary.thresholds["min_matched_peaks"] == 3 and summary.passed and caffeine(layer).evidence["library"]["matched_peaks"] == 3
    layer, summary = run(peaks=three_of_four, min_matched_peaks=6)                             # more than the spectrum has: still a fail, never an error
    assert summary.thresholds["min_matched_peaks"] == 4 and not summary.passed


def test_a_hit_without_an_adduct_fails_the_adduct_check(tiny):
    layer, summary = run(peaks=THEOBROMINE_PEAKS, mz=181.0720, formulas=[])
    top = layer.candidates[0]
    assert top.name == "Theobromine without adduct" and top.adduct is None and top.formula == "C7H8N4O2"
    assert failed_labels(top) == ["Adduct"] and top.evidence["library"]["failed_checks"][0] == "Adduct: the library entry records no adduct"
    assert top.evidence["library"]["library_adduct"] is None and top.mz_error_ppm is None and not summary.passed          # no adduct, no theoretical m/z


def test_a_hit_without_a_structure_is_kept_with_smiles_none(tiny):
    layer, summary = run(peaks=CAFFEINE_B, formulas=[])                                         # no formula list: the library's own precursor m/z stands in
    bare = next(c for c in layer.candidates if c.name == "Structure-less caffeine")
    assert (bare.smiles, bare.inchikey, bare.ik14, bare.formula, bare.mz_error_ppm) == (None, None, None, None, None)
    ev = bare.evidence["library"]
    assert ev["gate_passed"] and ev["precursor_error_basis"] == "the library's recorded precursor m/z" and ev["checks"][4]["status"] == "not evaluated"
    assert merge_candidates([L.LayerResult("library", candidates=[bare])])[0].flags == ["no_structure"]
    layer, _ = run(peaks=CAFFEINE_B)                                                          # with a formula list its formula is unknown: fail
    bare = next(c for c in layer.candidates if c.name == "Structure-less caffeine")
    assert failed_labels(bare) == ["Formula"] and "formula is unknown" in bare.evidence["library"]["failed_checks"][0]


def test_a_structure_the_library_stores_as_n_a_becomes_none_and_does_not_merge_unrelated_hits(tiny):
    layer, summary = run(peaks=UNKNOWN_PEAKS, mz=300.1000, positive=False, adducts_=("[M-H]-",), formulas=[])
    assert sorted(c.name for c in layer.candidates) == ["Another unknown with N/A structure", "Unknown with N/A structure"]      # "N/A" is no shared key
    for c in layer.candidates:
        assert (c.smiles, c.inchikey, c.ik14, c.formula, c.mz_error_ppm) == (None, None, None, None, None)
        assert c.evidence["library"]["gate_passed"] and c.evidence["library"]["precursor_error_basis"] == "the library's recorded precursor m/z"
    assert merge_candidates([layer])[0].flags == ["no_structure"] and len(merge_candidates([layer])) == 2


def test_an_isotope_labelled_hit_is_judged_by_the_library_precursor_not_by_a_formula_list(tiny):
    layer, summary = run(peaks=CAFFEINE_D3_PEAKS, mz=198.1065, formulas=CAFFEINE_FORMULAS)      # the list holds unlabelled caffeine only
    d3 = layer.candidates[0]
    assert d3.name == "Caffeine-d3" and "[2H]" in d3.formula and d3.mz_error_ppm is None and summary.passed
    ev = d3.evidence["library"]
    assert ev["precursor_error_basis"] == "the library's recorded precursor m/z" and [c["status"] for c in ev["checks"]][4] == "not evaluated"
    assert "outside BUDDY's search space" in ev["checks"][4]["detail"]


def test_the_candidate_list_covers_other_adducts_and_charged_structures(tiny):
    na, summary = run(peaks=GLUCOSE_NA, mz=203.0526, adducts_=("[M+Na]+",), formulas=[F.FormulaCandidate("C6H12O6", "[M+Na]+", 1)])
    glucose = na.candidates[0]
    assert summary.passed and glucose.name == "Glucose sodium adduct" and glucose.formula == "C6H12O6" and glucose.adduct == "[M+Na]+"
    fa, summary = run(peaks=GLUCOSE_FA, mz=225.0612, positive=False, adducts_=("[M+HCOO]-",), formulas=[])
    assert summary.passed and fa.candidates[0].adduct == "[M+HCOO]-" and fa.candidates[0].evidence["library"]["library_adduct"] == "[M+HCOO]-"
    choline, summary = run(peaks=CHOLINE_PEAKS, mz=104.1070, adducts_=("[M]+",), formulas=[])
    assert summary.passed and choline.candidates[0].formula == "C5H14NO" and abs(choline.candidates[0].mz_error_ppm) < 5


def test_the_isotope_check_runs_with_a_usable_ms1_pattern(tiny):
    layer, summary = run(ms1_peaks=CAFFEINE_MS1)
    ev = layer.candidates[0].evidence["library"]
    assert ev["isotope"] == "pass" and ev["isotope_similarity"] > 0.99 and summary.checks[5]["status"] == "pass" and summary.passed
    bad, summary = run(ms1_peaks=[(195.0877, 100.0), (196.0910, 80.0), (197.0940, 60.0)])
    ev = bad.candidates[0].evidence["library"]
    assert ev["isotope"] == "fail" and ev["gate_passed"] is False and ev["failed_checks"][0].startswith("Isotope pattern: isotope similarity 0.")
    unusable, summary = run(ms1_peaks=[(250.0, 100.0), (251.0, 10.0)])                         # no peak at the precursor: MS1 is not used
    assert unusable.candidates[0].evidence["library"]["isotope"] == "not evaluated" and summary.passed


def test_without_a_library_the_step_is_skipped_with_a_plain_note(tiny, no_libs):
    layer, summary = run(libraries=())
    assert layer.status == "skipped" and layer.candidates == [] and "No library was selected" in layer.note and not summary.passed
    layer, summary = run(libraries=("MassBank", "GNPS"))
    assert layer.status == "skipped" and "not installed" in layer.note and "library_setup" in layer.note and layer.candidates == []
    assert summary.passed is False and summary.thresholds["similarity"] == 0.8
    assert run(peaks=[])[0].status == "skipped"


def test_a_library_that_is_not_installed_only_adds_a_warning(tiny):
    layer, summary = run(libraries=("MassBank", "GNPS"))
    assert layer.status == "done" and summary.passed
    assert any("GNPS is not installed" in w and "Searching MassBank only" in w for w in summary.warnings)


def test_entropy_that_cannot_be_computed_falls_back_to_cosine(tiny, monkeypatch):
    real = ls.search
    seen = []

    def search(peaks, kind, libraries, metrics, *a, **kw):
        seen.append(metrics[0])
        if metrics[0] == L.PRIMARY_METRIC:
            raise ImportError("No module named 'ms_entropy'")
        return real(peaks, kind, libraries, metrics, *a, **kw)
    monkeypatch.setattr(ls, "search", search)
    layer, summary = run()
    assert seen == [L.PRIMARY_METRIC, L.FALLBACK_METRIC] and summary.passed and summary.thresholds["metric"] == "CosineGreedy"
    ev = layer.candidates[0].evidence["library"]
    assert ev["metric"] == "CosineGreedy" and "CosineGreedy" in ev["checks"][0]["detail"]


def test_the_search_gets_the_chosen_metric_the_precursor_and_the_list_threshold(tiny, monkeypatch):
    real = ls.search
    seen = []

    def search(peaks, kind, libraries, metrics, threshold, **kw):
        seen.append((metrics, threshold, kw["precursor_mz"]))
        return real(peaks, kind, libraries, metrics, threshold, **kw)
    monkeypatch.setattr(ls, "search", search)
    layer, summary = run(metric="ModifiedCosineGreedy", gate_similarity=0.95)
    assert seen == [(["ModifiedCosineGreedy"], 0.85, CAFFEINE_MZ)]
    ev = layer.candidates[0].evidence["library"]
    assert ev["metric"] == "ModifiedCosineGreedy" and summary.passed and "ModifiedCosineGreedy" in ev["checks"][0]["detail"]
    assert summary.thresholds["metric"] == "ModifiedCosineGreedy" and summary.thresholds["similarity"] == 0.95 and summary.thresholds["list"] == 0.85


def test_a_higher_gate_threshold_than_any_hit_scores_keeps_the_hits_but_fails_the_gate(tiny):
    layer, summary = run(metric="CosineGreedy", gate_similarity=1.0)
    assert layer.candidates and summary.thresholds["list"] == 0.9
    assert summary.passed == all(c.evidence["library"]["gate_passed"] for c in layer.candidates[:1])


def test_a_crashing_search_is_a_failed_step_not_an_exception(tiny, monkeypatch):
    def boom(*a, **kw):
        raise RuntimeError("index is corrupt")
    monkeypatch.setattr(ls, "search", boom)
    layer, summary = run()
    assert layer.status == "failed" and "RuntimeError: index is corrupt" in layer.note and layer.candidates == [] and not summary.passed


def test_the_search_is_asked_for_the_primary_metric_the_list_threshold_and_every_adduct(tiny, monkeypatch):
    calls = []
    real = ls.search

    def spy(peaks, kind, libraries, metrics, threshold, **kw):
        calls.append((kind, list(libraries), list(metrics), threshold, kw))
        return real(peaks, kind, libraries, metrics, threshold, **kw)
    monkeypatch.setattr(ls, "search", spy)
    run(tol_ppm=7.0, tol_da=0.01)
    kind, libraries, metrics, threshold, kw = calls[0]
    assert (kind, libraries, metrics, threshold) == ("MS/MS", ["MassBank"], ["SpectralEntropy (weighted)"], 0.7)
    assert kw["adduct"] is None and kw["precursor_tol_ppm"] == 7.0 and kw["tolerance"] == 0.01 and kw["positive"] is True


def test_library_results_are_json_and_feed_the_merge(tiny):
    layer, summary = run()
    json.dumps(layer.to_dict())
    json.dumps(summary.to_dict())
    merged = merge_candidates([layer])
    assert [assign_level(c) for c in merged][0] == "2a" and merged[0].sources[0]["layer"] == "library"


# ================================================================ real engines (opt-in)
needs_real = pytest.mark.skipif(not RUN_REAL, reason="real BUDDY / library runs are opt-in: set MZDESK_RUN_ANNOTATION_STEPS=1 "
                                                     "(the first BUDDY call takes 15-20 s)")


def _real_libraries():
    return [name for name in ("MassBank", "GNPS") if lx.is_installed(name)]


@needs_real
def test_real_caffeine_example_passes_the_gate():
    libs = _real_libraries()
    if not libs:
        pytest.skip("no local library installed (python -m mzdesk.ms2.library_setup)")
    adduct, formula = F.run_adduct_and_formula(CAFFEINE_MZ, True, CAFFEINE_MS2, CAFFEINE_MS1, tol_ppm=10.0)
    assert adduct.retained_adducts[0] == "[M+H]+" and formula.status == "done" and formula.formulas[0].formula == "C8H10N4O2"
    layer, summary = L.run_library(CAFFEINE_MS2, CAFFEINE_MZ, True, libs, adducts=adduct.retained_adducts, formulas=formula.formulas,
                                   ms1_peaks=CAFFEINE_MS1)
    assert summary.passed and layer.candidates[0].name.lower().startswith("caffeine") and layer.candidates[0].evidence["library"]["isotope"] == "pass"


@needs_real
def test_real_perturbed_spectrum_is_listed_but_fails_the_gate():
    libs = _real_libraries()
    if not libs:
        pytest.skip("no local library installed (python -m mzdesk.ms2.library_setup)")
    # the base peak cut to 15 % and two strong foreign peaks added: the real caffeine spectra score about 0.76 (MassBank 2026.03, GNPS 2026-03-29)
    query = [(mz, inten * (0.15 if abs(mz - 138.0663) < 0.01 else 1.0)) for mz, inten in CAFFEINE_MS2] + [(60.0, 300.0), (90.0, 300.0)]
    adduct, formula = F.run_adduct_and_formula(CAFFEINE_MZ, True, query, None, tol_ppm=10.0)
    layer, summary = L.run_library(query, CAFFEINE_MZ, True, libs, adducts=adduct.retained_adducts, formulas=formula.formulas)
    assert layer.candidates and not summary.passed
    assert all(any(f.startswith("Similarity") for f in c.evidence["library"]["failed_checks"]) for c in layer.candidates)


@needs_real
def test_real_library_match_confirms_an_adduct_the_user_did_not_choose():
    libs = _real_libraries()
    if not libs:
        pytest.skip("no local library installed (python -m mzdesk.ms2.library_setup)")
    adduct, formula = F.run_adduct_and_formula(CAFFEINE_MZ, True, CAFFEINE_MS2, None, adduct="[M+Na]+", tol_ppm=10.0)      # wrong on purpose
    layer, summary = L.run_library(CAFFEINE_MS2, CAFFEINE_MZ, True, libs, adducts=adduct.retained_adducts, formulas=formula.formulas)
    assert summary.passed and layer.candidates[0].adduct == "[M+H]+" and layer.candidates[0].evidence["library"]["adduct_source"] == "library"
    assert summary.checks[3]["status"] == "pass" and "confirms it" in summary.checks[3]["detail"]
