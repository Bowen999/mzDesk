"""ESI MS/MS Annotation core, part 2: the re-rank (mzdesk/annotation/rerank.py).

The arithmetic (base score, the three factors, the sort) is tested with the three platform_fit checks replaced by tables, so it needs
no MolGpKa weights. A few tests at the end run the real rules on caffeine, glucose, a phosphatidylcholine and sulfuric acid; they skip
without the MolGpKa weights (python -m mzdesk.platform_fit.pka_setup), the same way the Adducts and Ionization tests do.
"""
import copy
import json
from pathlib import Path

import pytest
from rdkit import Chem

from mzdesk import adducts
from mzdesk.annotation import rerank as rr
from mzdesk.annotation.merge import merge_candidates
from mzdesk.annotation.rerank import RerankConfig, rerank
from mzdesk.annotation.schema import Candidate, LayerResult
from mzdesk.platform_fit import adduct_rules as ar
from mzdesk.platform_fit import pka_engine as pe

needs_pka = pytest.mark.skipif(not Path(pe.weights_path()).is_file(), reason="MolGpKa weights not installed (python -m mzdesk.platform_fit.pka_setup)")

CAFFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
GLUCOSE = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"
POPC = r"CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\CCCCCCCC"
SULFURIC = "OS(=O)(=O)O"                      # no carbon: neither platform
INVALID = "this is not a smiles"
SMILES = [CAFFEINE, GLUCOSE, POPC, "CCO", "CCCO", "CCCCO"]      # distinct valid structures for the table tests


CACHED_CHECKS = (rr._omics_check, rr._ion_check, rr._predicted_ions)          # the real, cached functions (a test may patch the names)


@pytest.fixture(autouse=True)
def fresh_caches():
    """The three checks are cached per structure; a test that stubs the rules must not leave its answers behind."""
    for check in CACHED_CHECKS:
        check.cache_clear()
    yield
    for check in CACHED_CHECKS:
        check.cache_clear()


STUB_SCORE = {"high": 80.0, "medium": 30.0, "low": 10.0, "not expected": 0.0}


@pytest.fixture
def rules(monkeypatch):
    """Replace the platform_fit checks by tables keyed by SMILES (default: metabolome / strong / high). Returns (tables, adduct call log).

    tables["adduct"][smiles] is a likelihood, or a dict {additive: likelihood or (likelihood, rule score)}; an additive missing from the dict is "not expected".
    """
    tables = {"omics": {}, "ion": {}, "adduct": {}}
    log = []
    monkeypatch.setattr(rr, "_omics_check", lambda smiles: (tables["omics"].get(smiles, "metabolome"), "Stub verdict."))
    monkeypatch.setattr(rr, "_ion_check", lambda smiles, sign, ph: (tables["ion"].get(smiles, "strong"), "stub reason"))

    def adduct_check(smiles, adduct, additive, alkali):
        log.append({"smiles": smiles, "adduct": adduct, "additive": additive, "alkali": alkali})
        answer = tables["adduct"].get(smiles, "high")
        if isinstance(answer, dict):
            answer = answer.get(additive, "not expected")
        likelihood, score = answer if isinstance(answer, tuple) else (answer, STUB_SCORE[answer])
        return likelihood, score, "stub reason"
    monkeypatch.setattr(rr, "_adduct_check", adduct_check)
    return tables, log


def cand(cid, smiles, norm, *, layer="library", adduct="[M+H]+", hypothesis=False, sources=None, gate=True):
    c = Candidate(id=cid, layer=layer, smiles=smiles, name=cid, formula="C2H6O", adduct=adduct, layer_score=norm, layer_score_norm=norm,
                  evidence={layer: {"gate_passed": gate} if layer == "library" else {}}, hypothesis=hypothesis)
    c.sources = sources if sources is not None else [{"layer": layer, "layer_score": norm, "layer_score_norm": norm, "rank_in_layer": 1}]
    c.level = {"library": "2a", "msfinder": "3", "frigid": "H"}[layer]
    c.level_reason = "reason"
    return c


def run(cands, **kw):
    kw = {"polarity": "Positive", "adduct": "[M+H]+", "omics": "Either / unknown", **kw}
    return rerank(cands, **kw)


# ---------------------------------------------------------------- the owner's start values
def test_default_constants_are_the_start_values_of_the_spec():
    c = RerankConfig()
    assert (c.omics_match, c.omics_mismatch, c.omics_neither) == (1.0, 0.5, 0.2)
    assert (c.ion_strong, c.ion_moderate, c.ion_weak) == (1.0, 0.8, 0.5)
    assert (c.adduct_high, c.adduct_medium, c.adduct_low, c.adduct_not_expected) == (1.0, 0.75, 0.4, 0.2)
    assert (c.layer_bonus, c.base_cap, c.hypothesis_factor, c.not_evaluated) == (0.05, 1.0, 0.7, 0.8)


# ---------------------------------------------------------------- base score
def test_base_is_the_best_norm_plus_a_bonus_per_extra_layer_capped_at_one(rules):
    def sources(*norms):
        return [{"layer": layer, "layer_score": n, "layer_score_norm": n, "rank_in_layer": 1} for layer, n in zip(("library", "msfinder", "frigid"), norms)]
    one = cand("a", CAFFEINE, 0.70)
    two = cand("b", GLUCOSE, 0.90, sources=sources(0.60, 0.90))
    three = cand("c", POPC, 0.97, sources=sources(0.97, 0.50, 0.40))
    low = cand("d", "CCO", 0.40, sources=sources(0.30, 0.40, 0.20))
    base = {c.id: c.base_score for c in run([one, two, three, low])}
    assert base == {"a": 0.7, "b": 0.95, "c": 1.0, "d": 0.5}


def test_a_pure_frigid_hypothesis_is_multiplied_by_the_hypothesis_factor_and_a_supported_one_is_not(rules):
    only = cand("only", CAFFEINE, 0.80, layer="frigid", hypothesis=True)
    backed = cand("backed", GLUCOSE, 0.80, layer="frigid", hypothesis=False,
                  sources=[{"layer": "frigid", "layer_score_norm": 0.8, "rank_in_layer": 1}, {"layer": "msfinder", "layer_score_norm": 0.5, "rank_in_layer": 1}])
    out = {c.id: c for c in run([only, backed])}
    assert out["only"].base_score == 0.56 and out["backed"].base_score == 0.85


def test_a_candidate_that_was_not_merged_uses_its_own_layer_and_score(rules):
    c = cand("a", CAFFEINE, 0.66, sources=[])
    assert run([c])[0].base_score == 0.66


# ---------------------------------------------------------------- the three factors
@pytest.mark.parametrize("chosen, cls, factor", [
    ("Lipidomics", "lipidome", 1.0), ("Lipidomics", "metabolome", 0.5), ("Lipidomics", "both", 1.0), ("Lipidomics", "neither", 0.2),
    ("Metabolomics", "metabolome", 1.0), ("Metabolomics", "lipidome", 0.5), ("Metabolomics", "both", 1.0), ("Metabolomics", "neither", 0.2),
    ("Either / unknown", "lipidome", 1.0), ("Either / unknown", "metabolome", 1.0), ("Either / unknown", "both", 1.0), ("Either / unknown", "neither", 0.2),
])
def test_omics_factor(rules, chosen, cls, factor):
    tables, _ = rules
    tables["omics"][CAFFEINE] = cls
    out = run([cand("a", CAFFEINE, 0.8)], omics=chosen)[0]
    assert out.detect["omics"] == cls and out.detect["factors"]["omics"] == factor


@pytest.mark.parametrize("level, factor", [("strong", 1.0), ("moderate", 0.8), ("weak", 0.5)])
def test_ionization_factor(rules, level, factor):
    rules[0]["ion"][CAFFEINE] = level
    out = run([cand("a", CAFFEINE, 0.8)])[0]
    assert out.detect["ionizable"] == level and out.detect["factors"]["ionizable"] == factor


@pytest.mark.parametrize("likelihood, factor", [("high", 1.0), ("medium", 0.75), ("low", 0.4), ("not expected", 0.2)])
def test_adduct_factor(rules, likelihood, factor):
    rules[0]["adduct"][CAFFEINE] = likelihood
    out = run([cand("a", CAFFEINE, 0.8)])[0]
    assert out.detect["adduct"] == likelihood and out.detect["factors"]["adduct"] == factor


def test_final_score_is_the_base_times_the_three_factors(rules):
    tables, _ = rules
    tables["omics"][CAFFEINE], tables["ion"][CAFFEINE], tables["adduct"][CAFFEINE] = "lipidome", "weak", "low"
    out = run([cand("a", CAFFEINE, 0.8)], omics="Metabolomics")[0]
    assert out.detect["factor"] == pytest.approx(0.5 * 0.5 * 0.4) and out.base_score == 0.8
    assert out.final_score == pytest.approx(0.8 * 0.1)
    assert len(out.detect["reasons"]) == 3 and out.detect["reasons"][0].startswith("omics:")
    assert out.detect["reasons"][1].startswith("ionization: ESI+ response weak") and out.detect["reasons"][2].startswith("adduct: [M+H]+ low")


def test_the_polarity_picks_the_esi_response(monkeypatch):
    seen = []
    monkeypatch.setattr(rr, "_omics_check", lambda smiles: ("metabolome", "v"))
    monkeypatch.setattr(rr, "_adduct_check", lambda *a: ("high", 80.0, "r"))
    monkeypatch.setattr(rr, "_ion_check", lambda smiles, sign, ph: (seen.append((sign, ph)) or "strong", "r"))
    run([cand("a", CAFFEINE, 0.8)], polarity="negative", adduct="[M-H]-", ph=6.5)
    run([cand("a", CAFFEINE, 0.8)], polarity="+")
    assert seen == [("-", 6.5), ("+", rr.DEFAULT_PH)]


def test_a_config_changes_the_factors(rules):
    rules[0]["omics"][CAFFEINE] = "lipidome"
    out = run([cand("a", CAFFEINE, 0.8)], omics="Metabolomics", config=RerankConfig(omics_mismatch=0.1, base_cap=0.5))[0]
    assert out.detect["factors"]["omics"] == 0.1 and out.base_score == 0.5


# ---------------------------------------------------------------- ranking
def test_a_candidate_impossible_for_the_omics_drops_below_a_plausible_one_in_the_window_but_stays(rules):
    tables, _ = rules
    tables["omics"][CAFFEINE], tables["omics"][POPC] = "neither", "lipidome"
    impossible, plausible = cand("impossible", CAFFEINE, 0.95), cand("plausible", POPC, 0.90, layer="msfinder")
    out = run([impossible, plausible], omics="Lipidomics")
    assert [c.id for c in out] == ["plausible", "impossible"]
    assert out[1].final_score == pytest.approx(0.95 * 0.2) and out[0].final_score == pytest.approx(0.90)
    assert (out[1].layer, out[1].level, out[1].level_reason) == ("library", "2a", "reason")        # nothing is relabelled


def test_reranking_never_deletes_or_relabels_and_does_not_touch_the_input(rules):
    tables, _ = rules
    tables["omics"][GLUCOSE], tables["ion"][POPC], tables["adduct"]["CCO"] = "neither", "weak", "not expected"
    cands = [cand("a", CAFFEINE, 0.9), cand("b", GLUCOSE, 0.8, layer="msfinder"), cand("c", POPC, 0.7, layer="frigid", hypothesis=True),
             cand("d", "CCO", 0.6), cand("e", INVALID, 0.5, layer="msfinder"), cand("f", None, 0.4)]
    before = copy.deepcopy(cands)
    out = run(cands)
    assert cands == before                                                  # the input list and its candidates are untouched
    assert sorted(c.id for c in out) == sorted(c.id for c in cands)
    by_id = {c.id: c for c in cands}
    for c in out:
        old = by_id[c.id]
        assert (c.layer, c.level, c.level_reason, c.sources, c.evidence, c.smiles, c.name, c.hypothesis, c.adduct) == \
               (old.layer, old.level, old.level_reason, old.sources, old.evidence, old.smiles, old.name, old.hypothesis, old.adduct)
    assert [c.rank for c in out] == [1, 2, 3, 4, 5, 6]
    assert [c.base_score for c in out] == sorted((c.base_score for c in out), reverse=True)      # here only "a" is in the window, the rest follows by base score


def test_ties_keep_the_input_order_and_a_higher_base_wins_a_tie_of_the_final_score(rules):
    tables, _ = rules
    first, second = cand("first", CAFFEINE, 0.8), cand("second", GLUCOSE, 0.8)
    assert [c.id for c in run([first, second])] == ["first", "second"]
    assert [c.id for c in run([second, first])] == ["second", "first"]
    tables["adduct"][CAFFEINE] = "medium"
    high_base = cand("high_base", CAFFEINE, 0.8)            # 0.8 x 0.75 = 0.6
    low_base = cand("low_base", POPC, 0.6)                  # 0.6 x 1.0 = 0.6: the same final score, listed first
    out = run([low_base, high_base], omics="Metabolomics")
    assert [c.id for c in out] == ["high_base", "low_base"] and out[0].final_score == out[1].final_score == 0.6


# ---------------------------------------------------------------- the tie-break window
def test_default_window_is_ten_percent():
    assert RerankConfig().tiebreak_window == 0.10


def test_inside_the_window_the_factors_demote_a_candidate(rules):
    tables, _ = rules
    tables["omics"][CAFFEINE] = "neither"                                   # factor 0.2
    top, near, edge = cand("top", CAFFEINE, 0.95), cand("near", GLUCOSE, 0.90), cand("edge", POPC, 0.86)       # 0.855 is the edge of the window
    out = run([top, near, edge])
    assert [c.id for c in out] == ["near", "edge", "top"] and [c.rank for c in out] == [1, 2, 3]
    assert out[2].final_score == pytest.approx(0.19) and out[0].final_score == pytest.approx(0.90)
    assert len(out[2].detect["reasons"]) == 3                               # nobody here is outside the window


def test_a_candidate_far_behind_is_never_promoted_above_the_window(rules):
    tables, _ = rules
    tables["omics"][CAFFEINE], tables["ion"][CAFFEINE], tables["adduct"][CAFFEINE] = "neither", "weak", "not expected"      # 0.2 x 0.5 x 0.2 = 0.02
    bad, far = cand("bad", CAFFEINE, 0.95), cand("far", GLUCOSE, 0.80)       # 0.80 < 0.855
    out = run([bad, far])
    assert [c.id for c in out] == ["bad", "far"]                            # a perfectly detectable 0.80 does not beat a 0.95 that is nearly undetectable
    assert out[0].final_score == pytest.approx(0.019) and out[0].detect["factor"] == pytest.approx(0.02)
    assert out[1].final_score == out[1].base_score == 0.80                  # outside the window: final = base, and it says so
    assert out[1].detect["factors"] == {"omics": 1.0, "ionizable": 1.0, "adduct": 1.0} and len(out[1].detect["reasons"]) == 4
    assert out[1].detect["reasons"][3].startswith("re-rank: base score 0.80") and len(out[0].detect["reasons"]) == 3


def test_outside_the_window_the_order_is_the_base_score_whatever_the_factors_say(rules):
    tables, _ = rules
    tables["omics"][POPC] = "neither"
    out = run([cand("best", CAFFEINE, 0.95), cand("b", GLUCOSE, 0.60), cand("c", POPC, 0.70), cand("d", "CCO", 0.50)])
    assert [c.id for c in out] == ["best", "c", "b", "d"]                   # c is "neither" (x0.2) but 0.70 > 0.60: base order
    assert [c.rank for c in out] == [1, 2, 3, 4]
    assert [c.final_score for c in out[1:]] == [c.base_score for c in out[1:]]
    assert out[1].detect["omics"] == "neither"                              # the details are still computed and shown


def test_window_edge_ties_and_width(rules):
    tables, _ = rules
    tables["omics"][CAFFEINE] = "neither"
    assert [c.id for c in run([cand("edge", CAFFEINE, 0.90), cand("top", GLUCOSE, 1.0)])] == ["top", "edge"]      # 0.90 = (1 - 0.10) x 1.0: inside
    tables["omics"][CAFFEINE] = "metabolome"
    tie = [cand("a", CAFFEINE, 0.8), cand("b", GLUCOSE, 0.8)]
    assert [c.id for c in run(tie)] == ["a", "b"] and [c.id for c in run(tie[::-1])] == ["b", "a"]               # an exact tie keeps the input order
    tables["omics"][CAFFEINE] = "neither"
    pair = [cand("a", CAFFEINE, 0.95), cand("b", GLUCOSE, 0.90)]
    assert [c.id for c in run(pair, config=RerankConfig(tiebreak_window=0.0))] == ["a", "b"]         # window 0: only the best is in it
    assert [c.id for c in run(pair)] == ["b", "a"]


def test_a_single_candidate_is_alone_in_the_window(rules):
    rules[0]["omics"][CAFFEINE] = "neither"
    (only,) = run([cand("only", CAFFEINE, 0.5)])
    assert only.rank == 1 and only.final_score == pytest.approx(0.1) and len(only.detect["reasons"]) == 3


def test_the_window_must_be_a_fraction():
    for bad in (-0.1, 1.0):
        with pytest.raises(ValueError, match="window"):
            run([cand("a", CAFFEINE, 0.5)], config=RerankConfig(tiebreak_window=bad))


def test_no_candidates_give_an_empty_list():
    assert run([]) == []


def test_progress_is_reported_after_each_structure(rules):
    seen = []
    cands = [cand("a", CAFFEINE, 0.9), cand("b", GLUCOSE, 0.8), cand("c", INVALID, 0.7)]
    out = run(cands, progress=lambda done, total: seen.append((done, total)))
    assert seen == [(1, 3), (2, 3), (3, 3)] and len(out) == 3
    seen.clear()
    run([], progress=lambda done, total: seen.append((done, total)))
    assert seen == []
    assert [c.id for c in run(iter(cands))] == ["a", "b", "c"]               # any iterable, not only a list


# ---------------------------------------------------------------- what cannot be scored
def test_a_candidate_without_a_valid_smiles_is_kept_and_counts_the_not_evaluated_factor_for_each_dimension(rules):
    out = run([cand("bad", INVALID, 0.70), cand("none", None, 0.68), cand("blank", "  ", 0.66)])
    assert [c.id for c in out] == ["bad", "none", "blank"]
    for c in out:
        assert (c.detect["omics"], c.detect["ionizable"], c.detect["adduct"]) == ("not evaluated",) * 3
        assert c.detect["factors"] == {"omics": 0.8, "ionizable": 0.8, "adduct": 0.8} and c.detect["factor"] == pytest.approx(0.512)
        assert c.final_score == pytest.approx(c.base_score * 0.512, abs=1e-4)
        assert all("no valid structure" in r for r in c.detect["reasons"])
    assert rules[1] == []                                                   # the adduct rules were never asked


def test_a_rule_that_cannot_run_is_reported_and_counts_the_not_evaluated_factor(monkeypatch):
    def boom(*args, **kwargs):
        raise ValueError("The structure is too large.")
    monkeypatch.setattr(rr.ir, "profile", boom)
    monkeypatch.setattr(rr.ar, "predict", boom)
    out = run([cand("a", CAFFEINE, 0.8)], omics="Metabolomics")[0]
    assert (out.detect["ionizable"], out.detect["adduct"]) == ("not evaluated", "not evaluated")
    assert out.detect["omics"] == "metabolome" and out.detect["factors"] == {"omics": 1.0, "ionizable": 0.8, "adduct": 0.8}
    assert out.final_score == pytest.approx(0.8 * 0.64) and out.base_score == 0.8
    assert "ValueError: The structure is too large." in out.detect["reasons"][1] and "ValueError" in out.detect["reasons"][2]


def test_without_any_adduct_the_adduct_is_not_scored(rules):
    out = run([cand("a", CAFFEINE, 0.8, adduct=None)], adduct=None)[0]
    assert out.detect["adduct"] == "not evaluated" and out.detect["factors"]["adduct"] == 0.8 and rules[1] == []
    assert "no adduct given" in out.detect["reasons"][2]


def test_a_candidates_own_adduct_wins_over_the_run_adduct(rules):
    run([cand("a", CAFFEINE, 0.8, adduct="[M+Na]+"), cand("b", GLUCOSE, 0.7, adduct=None)], adduct="[M+H]+")
    assert [(r["smiles"], r["adduct"]) for r in rules[1]] == [(CAFFEINE, "[M+Na]+"), (GLUCOSE, "[M+H]+")]


def test_adduct_aliases_are_made_canonical_before_the_rules_are_asked(rules):
    own = cand("own", CAFFEINE, 0.8, adduct="[M+FA-H]-")
    assert own.adduct == "[M+HCOO]-"                                       # canonical from the moment the Candidate exists
    run([own, cand("run", GLUCOSE, 0.7, adduct=None)], polarity="Negative", adduct="[M+Hac-H]-")
    assert [r["adduct"] for r in rules[1]] == ["[M+HCOO]-", "[M+CH3COO]-"]
    assert run([cand("auto", CAFFEINE, 0.8, adduct=None)], adduct="Auto")[0].detect["adduct"] == "not evaluated"


# ---------------------------------------------------------------- the mobile-phase context
ALL_ADDITIVES = ["Formic acid", "Ammonium formate", "Ammonium acetate", "Acetic acid", "None"]


@pytest.mark.parametrize("adduct, tried", [
    ("[M+H]+", ALL_ADDITIVES), ("[M-H]-", ALL_ADDITIVES), ("[M+Na]+", ALL_ADDITIVES), ("[M+Cl]-", ALL_ADDITIVES),
    ("[M+NH4]+", ["Ammonium formate", "Ammonium acetate"]),
    ("[M+HCOO]-", ["Formic acid", "Ammonium formate"]),
    ("[M+CH3COO]-", ["Ammonium acetate", "Acetic acid"]),
])
def test_an_unknown_additive_is_searched_over_the_additives_that_can_offer_the_adduct_and_alkali_counts_as_present(rules, adduct, tried):
    rules[0]["adduct"][CAFFEINE] = "low"                                    # never "high", so nothing stops the search
    out = run([cand("a", CAFFEINE, 0.8, adduct=adduct)])[0]
    assert [(r["adduct"], r["additive"], r["alkali"]) for r in rules[1]] == [(adduct, name, True) for name in tried]
    assert out.detect["adduct"] == "low" and out.detect["factors"]["adduct"] == 0.4


def test_an_unknown_additive_stops_at_the_first_high_and_otherwise_keeps_the_best(rules):
    tables, log = rules
    tables["adduct"][CAFFEINE] = {"Formic acid": "low", "Ammonium formate": "medium", "Ammonium acetate": "high", "Acetic acid": "high", "None": "high"}
    out = run([cand("a", CAFFEINE, 0.8)])[0]
    assert [r["additive"] for r in log] == ["Formic acid", "Ammonium formate", "Ammonium acetate"]       # the search stopped at the first high
    assert out.detect["adduct"] == "high" and out.detect["reasons"][2] == "adduct: [M+H]+ high; with Ammonium acetate; rule score 80; stub reason"

    log.clear()
    tables["adduct"][GLUCOSE] = {"Formic acid": ("low", 5.0), "Ammonium formate": ("low", 15.0), "Ammonium acetate": ("low", 12.0), "None": "medium"}
    out = run([cand("b", GLUCOSE, 0.8)])[0]
    assert len(log) == 5 and out.detect["adduct"] == "medium" and "with None; rule score 30" in out.detect["reasons"][2]       # a higher class beats a higher score

    tables["adduct"][POPC] = {"Formic acid": ("low", 5.0), "Ammonium formate": ("low", 15.0), "Ammonium acetate": ("low", 12.0)}
    out = run([cand("c", POPC, 0.8)])[0]
    assert out.detect["adduct"] == "low" and "with Ammonium formate; rule score 15" in out.detect["reasons"][2]               # then the higher rule score

    tables["adduct"]["CCO"] = "not expected"
    out = run([cand("d", "CCO", 0.8)])[0]
    assert out.detect["adduct"] == "not expected" and out.detect["factors"]["adduct"] == 0.2
    assert out.detect["reasons"][2].endswith("(checked with Formic acid, Ammonium formate, Ammonium acetate, Acetic acid, None)")


def test_an_explicit_additive_is_the_only_one_asked_and_alkali_is_passed_through(rules):
    rules[0]["adduct"][CAFFEINE] = {"Formic acid": "low", "Ammonium formate": "high"}
    out = run([cand("a", CAFFEINE, 0.8, adduct="[M+NH4]+")], additive="Formic acid", alkali=False)[0]
    assert rules[1] == [{"smiles": CAFFEINE, "adduct": "[M+NH4]+", "additive": "Formic acid", "alkali": False}]
    assert out.detect["adduct"] == "low" and out.detect["reasons"][2] == "adduct: [M+NH4]+ low; with Formic acid; rule score 10; stub reason"
    rules[1].clear()
    run([cand("a", CAFFEINE, 0.8, adduct="[M+NH4]+")], additive="Unknown")                  # the form's own wording
    run([cand("a", CAFFEINE, 0.8, adduct="[M+NH4]+")], additive=None)
    assert [r["additive"] for r in rules[1]] == ["Ammonium formate", "Ammonium formate"]


def test_the_inputs_are_validated(rules):
    c = [cand("a", CAFFEINE, 0.8)]
    for bad, match in [({"polarity": "Both"}, "Polarity"), ({"omics": "Proteomics"}, "Omics"), ({"ph": 15}, "pH"), ({"additive": "TFA"}, "Additive")]:
        with pytest.raises(ValueError, match=match):
            run(c, **bad)
    assert run(c, omics=None)[0].detect["factors"]["omics"] == 1.0 and run(c, polarity="NEGATIVE", adduct="[M-H]-")


def test_a_structure_is_scored_once_per_additive_even_across_runs(monkeypatch):
    calls = {"profile": 0, "predict": 0}

    class Response:
        level, reason = "strong", "stub"

    class Ion:
        adduct, likelihood, score, reason = "[M+H]+", "high", 100.0, "stub"

    class Prediction:
        all_ions = [Ion()]

    def profile(smiles):
        calls["profile"] += 1
        return object()

    def predict(smiles, polarity, additive, alkali):
        calls["predict"] += 1
        return Prediction()
    monkeypatch.setattr(rr.ir, "profile", profile)
    monkeypatch.setattr(rr.ir, "esi_positive", lambda prof, ph: Response())
    monkeypatch.setattr(rr.ir, "esi_negative", lambda prof, ph: Response())
    monkeypatch.setattr(rr.ar, "predict", predict)
    cands = [cand("a", CAFFEINE, 0.8), cand("b", GLUCOSE, 0.7)]
    run(cands)
    run(cands)
    assert calls == {"profile": 2, "predict": 2}                             # [M+H]+ is high under formic acid: one prediction each, then the cache
    for _ in range(2):                                                       # [M+K]+ is not listed: all five additives are asked, formic acid is cached already
        assert run([cand("c", CAFFEINE, 0.5, adduct="[M+K]+")])[0].detect["adduct"] == "not expected"
    assert calls == {"profile": 2, "predict": 6}


# ---------------------------------------------------------------- the real rules
@needs_pka
def test_real_rules_for_a_metabolite_a_sugar_and_a_lipid_in_positive_mode():
    caf, glc, pc = (cand("caffeine", CAFFEINE, 0.9), cand("glucose", GLUCOSE, 0.9, adduct="[M+Na]+"), cand("pc", POPC, 0.9))
    out = {c.id: c.detect for c in run([caf, glc, pc], omics="Either / unknown")}
    assert out["caffeine"]["omics"] == "metabolome" and out["caffeine"]["adduct"] == "high" and out["caffeine"]["ionizable"] in ("moderate", "strong")
    assert out["glucose"]["omics"] == "metabolome" and out["glucose"]["adduct"] == "high"            # [M+Na]+ of a sugar
    assert (out["pc"]["omics"], out["pc"]["ionizable"], out["pc"]["adduct"]) == ("lipidome", "strong", "high")
    assert out["pc"]["factor"] == pytest.approx(1.0)                                                 # a lipid is detectable by one of the two platforms
    assert "Lipidome" in out["pc"]["reasons"][0] and out["pc"]["reasons"][1].startswith("ionization: ESI+ response strong")


@needs_pka
def test_real_rules_rank_a_lipid_above_a_higher_scored_metabolite_when_the_user_chose_lipidomics():
    out = run([cand("caffeine", CAFFEINE, 0.9), cand("pc", POPC, 0.85)], omics="Lipidomics")
    assert [c.id for c in out] == ["pc", "caffeine"]
    assert out[1].detect["factors"]["omics"] == 0.5
    out = run([cand("caffeine", CAFFEINE, 0.9), cand("pc", POPC, 0.85)], omics="Metabolomics")
    assert [c.id for c in out] == ["caffeine", "pc"]


@needs_pka
def test_real_rules_know_that_a_molecule_without_carbon_is_neither():
    out = run([cand("acid", SULFURIC, 0.9, adduct="[M-H]-")], omics="Metabolomics", polarity="Negative")[0]
    assert out.detect["omics"] == "neither" and out.detect["factors"]["omics"] == 0.2 and "Neither" in out.detect["reasons"][0]


@needs_pka
def test_real_rules_resolve_the_six_casmi_adducts_for_a_sugar_and_a_lipid():
    for smiles in (GLUCOSE, POPC):
        found = {a: run([cand("x", smiles, 0.9, adduct=a)], polarity=adducts.polarity(a))[0].detect["adduct"] for a in adducts.CASMI_ADDUCTS}
        assert set(found) == set(adducts.CASMI_ADDUCTS) and set(found.values()) <= {"high", "medium", "low", "not expected"}
    sugar = {a: run([cand("g", GLUCOSE, 0.9, adduct=a)], polarity=adducts.polarity(a))[0].detect["adduct"] for a in ("[M+H]+", "[M+Na]+", "[M+NH4]+", "[M+HCOO]-")}
    assert sugar == {"[M+H]+": "low", "[M+Na]+": "high", "[M+NH4]+": "high", "[M+HCOO]-": "high"}


@needs_pka
def test_real_rules_an_explicit_additive_is_not_improved_on_but_an_unknown_one_takes_the_best():
    def sugar(adduct, **kw):
        return run([cand("g", GLUCOSE, 0.9, adduct=adduct)], polarity=adducts.polarity(adduct), **kw)[0].detect
    assert sugar("[M+NH4]+")["adduct"] == "high"                                      # unknown: the ammonium additives are searched
    assert sugar("[M+NH4]+", additive="Formic acid")["adduct"] == "not expected"      # formic acid never offers NH4+
    assert sugar("[M+NH4]+", additive="Formic acid")["factors"]["adduct"] == 0.2
    assert sugar("[M+Cl]-")["adduct"] == "medium" and sugar("[M+Cl]-", additive="Formic acid")["adduct"] == "low"     # without an additive anion chloride binds best


@needs_pka
def test_real_rules_never_offer_an_adduct_under_an_additive_the_search_skips():
    """The unknown-additive search skips the additives that cannot supply the adduct's ion partner; this pins that the rules agree."""
    structures = [GLUCOSE, CAFFEINE, POPC, "CCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC", "CCO"]
    for adduct in ("[M+NH4]+", "[M+HCOO]-", "[M+CH3COO]-"):
        skipped = [name for name in ar.ADDITIVES if name not in rr._additives_for(adduct)]
        assert skipped and len(skipped) < len(ar.ADDITIVES)
        for smiles in structures:
            for name in skipped:
                assert adduct not in {i.adduct for i in ar.predict(smiles, "Both", name, True).all_ions}, (adduct, smiles, name)
    assert rr._additives_for("[M+NH4]+") == ["Ammonium formate", "Ammonium acetate"] and rr._additives_for("[M+H]+") == ALL_ADDITIVES


# ---------------------------------------------------------------- the MolGpKa pass is shared
STRUCTURES = {"caffeine": CAFFEINE, "glucose": GLUCOSE, "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "lidocaine": "CCN(CC)CC(=O)Nc1c(C)cccc1C",
              "quercetin": "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12", "metformin": "CN(C)C(=N)N=C(N)N", "carnitine": "C[N+](C)(C)C[C@@H](O)CC(=O)[O-]",
              "POPC": POPC, "sulfate": SULFURIC}


@pytest.mark.parametrize("ph", [2.7, 6.5])
def test_the_esi_verdict_is_the_one_ionization_rules_assess_gives(ph):
    """_ion_check leaves out the EI half of assess(); the ESI levels and reasons must stay exactly what assess() says."""
    for name, smiles in STRUCTURES.items():
        assessed = rr.ir.assess(smiles, ph)
        assert rr._ion_check(smiles, "+", ph) == (assessed.pos.level, assessed.pos.reason), name
        assert rr._ion_check(smiles, "-", ph) == (assessed.neg.level, assessed.neg.reason), name


@pytest.fixture
def counted_passes(monkeypatch):
    """Count the MolGpKa passes (calls of predict_aids: one for the acid sites and one for the base sites per molecule)."""
    pe.predict_sites.cache_clear()
    count = {"n": 0}
    real = pe.predict_aids

    def counting(*args, **kwargs):
        count["n"] += 1
        return real(*args, **kwargs)
    monkeypatch.setattr(pe, "predict_aids", counting)
    yield count
    pe.predict_sites.cache_clear()


@needs_pka
def test_one_molgpka_pass_serves_the_ionization_rules_and_every_additive_of_the_adduct_rules(counted_passes):
    rr._ion_check(GLUCOSE, "+", 2.7)
    first = counted_passes["n"]
    assert first == 2                                                         # acid sites and base sites of one molecule
    for additive in ar.ADDITIVES:
        rr._predicted_ions(GLUCOSE, additive, True)
    rr.ir.assess(GLUCOSE, 6.5)
    rr._ion_check(GLUCOSE, "-", 6.5)
    assert counted_passes["n"] == first


@needs_pka
def test_cached_sites_equal_the_uncached_ones_and_are_copies(counted_passes):
    for smiles in STRUCTURES.values():
        mol = pe.prepare(smiles).mol
        assert pe.predict_sites(mol) == pe.predict_sites.__wrapped__(mol) == pe.predict_sites(mol)
    mol = pe.prepare(CAFFEINE).mol
    sites = pe.predict_sites(mol)
    sites[0].pka, sites[0].kind = 99.0, "none"
    sites.clear()
    again = pe.predict_sites(mol)
    assert again and again == pe.predict_sites.__wrapped__(mol) and all(site.pka != 99.0 for site in again)


@needs_pka
def test_cached_sites_follow_the_atom_numbering_and_the_weights_file(counted_passes, monkeypatch):
    amine_first, amine_last = Chem.MolFromSmiles("NCCO"), Chem.MolFromSmiles("OCCN")           # one molecule, numbered two ways
    base_atom = lambda mol: [s.atom for s in pe.predict_sites(mol) if s.kind == "base"]        # noqa: E731
    assert base_atom(amine_first) == [0] and base_atom(amine_last) == [3] and base_atom(amine_first) == [0]
    monkeypatch.setenv("MZDESK_MOLGPKA_WEIGHTS", "/nonexistent/molgpka_weights.npz")           # without the weights the answer is an error, cached or not
    with pytest.raises(pe.PkaError):
        pe.predict_sites(amine_first)


# ---------------------------------------------------------------- merge and re-rank together
def test_merge_then_rerank_keeps_the_layers_and_levels(rules):
    layers = [LayerResult("library", candidates=[cand("l1", CAFFEINE, 0.93), cand("l2", GLUCOSE, 0.74, gate=False)]),
              LayerResult("msfinder", candidates=[cand("m1", GLUCOSE, 0.8, layer="msfinder")]),
              LayerResult("frigid", candidates=[cand("f1", POPC, 0.9, layer="frigid", hypothesis=True)])]
    for layer in layers:
        for c in layer.candidates:
            c.sources = []                                                 # a layer does not fill these
    merged = merge_candidates(layers)
    out = run(merged, omics="Either / unknown")
    assert {c.id: c.level for c in out} == {"c1": "2a", "c2": "3", "c3": "H"}
    by = {c.id: c for c in out}
    assert by["c2"].base_score == 0.85 and by["c3"].base_score == 0.63 and by["c1"].base_score == 0.93
    assert [c.rank for c in out] == [1, 2, 3]
    rows = json.loads(json.dumps([c.to_dict() for c in out], allow_nan=False))             # the rows of the result are strict JSON
    assert [r["id"] for r in rows] == [c.id for c in out] and rows[0]["detect"]["reasons"][0].startswith("omics:")
