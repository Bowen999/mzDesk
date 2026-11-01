"""ESI MS/MS Annotation with an observed CCS: the form, the CCS check of step 6 and the ranking rule.

  - rerank.apply_ccs is a pure rule (tested on hand-made candidates): measured-inconsistent goes behind all candidates that are not; predicted-inconsistent only
    moves inside the 10% near-tie window; everything else stays; the order is otherwise stable.
  - parse_inputs: the Ion Mobility fields (platform required with a CCS, 1/K0-like values refused).
  - step 6 through pipeline.run_query with the Rig of test_pipeline.py and a fake CCS index (tests/test_ccs.py): measured, predicted, He / FAIMS / TWIMS N2->He,
    index missing, no structure, a failing check.
  - the real CCS index (~/.mzdesk/data/ccs): caffeine's own [M+H]+ value is Consistent, +10 % of it is Inconsistent and moves caffeine behind (skipped when not installed).
  - the same without an observed CCS is tests/test_annotation_ccs_unchanged.py.
"""
import copy
import json
import os

import pytest

from mzdesk import spec
from mzdesk.annotation import ccs_step, pipeline, rerank
from mzdesk.annotation.rerank import apply_ccs
from mzdesk.annotation.schema import Candidate
from mzdesk.api import app as app_module
from mzdesk.api.common import InputError
from mzdesk.methods import compound_index as ci
from mzdesk.orthogonal import ccs
from test_ccs import ASPIRIN, CAFFEINE, IBUPROFEN, fake_ccs, fake_model, needs_real_ccs   # noqa: F401 - fixtures
from test_pipeline import (ESI, ETHER, form, isomer_result, library_candidate, rig,   # noqa: F401 - rig is a fixture
                           structure_candidate)

from fastapi.testclient import TestClient

client = TestClient(app_module.app)
INC, CON, NOREF = "Inconsistent", "Consistent", "No reference"
SINGLE_FIELD = "DTIMS (single field)"
ISOCAFFEINE = "Cn1c(=O)c2ncn(C)c2n(C)c1=O"


# ---------------------------------------------------------------- the rule on hand-made candidates
def cand(cid, base, verdict=None, source=None):
    c = Candidate(id=cid, layer="msfinder", smiles="CCO", name=cid, adduct="[M+H]+", layer_score_norm=base, level="3")
    c.base_score = c.final_score = base
    if verdict:
        c.ccs = {"verdict": verdict, "source": source, "flag": None, "moved_from": None}
    return c


def order(cands):
    return [c.id for c in cands]


def test_the_verdict_words_of_the_rule_are_the_ccs_engines():
    assert (rerank.CCS_INCONSISTENT, rerank.CCS_MEASURED, rerank.CCS_PREDICTED) == (ccs.INCONSISTENT, "measured", "predicted")
    assert ccs_step.NO_STRUCTURE == "No structure" and rerank.FLAG_CCS_MEASURED == "CCS ≠ measured" and rerank.FLAG_CCS_PREDICTED == "CCS ≠ predicted"


def test_measured_inconsistent_goes_behind_every_candidate_that_is_not():
    ranked = [cand("a", 1.0, INC, "measured"), cand("b", 0.95), cand("c", 0.5, CON, "measured"), cand("d", 0.4, NOREF)]
    out = apply_ccs(ranked)
    assert order(out) == ["b", "c", "d", "a"] and [c.rank for c in out] == [1, 2, 3, 4]
    a = out[-1]
    assert a.ccs["flag"] == "CCS ≠ measured" and a.ccs["moved_from"] == 1
    assert all(c.ccs is None or c.ccs["flag"] is None for c in out[:-1])


def test_a_measured_inconsistent_candidate_outside_the_window_is_demoted_too():
    out = apply_ccs([cand("a", 1.0), cand("b", 0.5, INC, "measured"), cand("c", 0.4)])
    assert order(out) == ["a", "c", "b"] and out[2].ccs["moved_from"] == 2


def test_the_measured_inconsistent_ones_keep_their_relative_order_and_so_do_the_rest():
    ranked = [cand("a", 1.0, INC, "measured"), cand("b", 0.99), cand("c", 0.98, INC, "measured"), cand("d", 0.97), cand("e", 0.2), cand("f", 0.1, INC, "measured")]
    assert order(apply_ccs(ranked)) == ["b", "d", "e", "a", "c", "f"]


def test_predicted_inconsistent_only_moves_inside_the_window():
    # window: base >= 0.9 x 1.0 = a, b, c. d and e are outside; e is predicted-inconsistent and stays where it is, without a flag
    ranked = [cand("a", 1.0, INC, "predicted"), cand("b", 0.95, CON, "measured"), cand("c", 0.92, INC, "predicted"), cand("d", 0.5, CON, "predicted"),
              cand("e", 0.45, INC, "predicted")]
    out = apply_ccs(ranked)
    assert order(out) == ["b", "a", "c", "d", "e"]
    by = {c.id: c.ccs for c in out}
    assert (by["a"]["flag"], by["a"]["moved_from"]) == ("CCS ≠ predicted", 1)
    assert (by["c"]["flag"], by["c"]["moved_from"]) == ("CCS ≠ predicted", None)          # the rule applies to it, but it was already behind the others
    assert by["e"]["flag"] is None and by["e"]["moved_from"] is None                       # outside the window: no effect at all


def test_a_flagged_candidate_that_gains_a_place_because_another_left_is_not_moved():
    out = apply_ccs([cand("a", 1.0, INC, "measured"), cand("b", 0.95, INC, "predicted"), cand("c", 0.3)])
    assert order(out) == ["b", "c", "a"] and [c.rank for c in out] == [1, 2, 3]
    assert (out[0].ccs["flag"], out[0].ccs["moved_from"]) == ("CCS ≠ predicted", None) and out[2].ccs["moved_from"] == 1


def test_predicted_inconsistent_stays_ahead_of_everything_outside_the_window():
    out = apply_ccs([cand("a", 1.0, INC, "predicted"), cand("b", 0.95), cand("c", 0.89), cand("d", 0.2)])
    assert order(out) == ["b", "a", "c", "d"]


def test_a_window_of_predicted_inconsistent_candidates_does_not_move():
    ranked = [cand("a", 1.0, INC, "predicted"), cand("b", 0.95, INC, "predicted"), cand("c", 0.3)]
    assert order(apply_ccs(ranked)) == ["a", "b", "c"]


def test_measured_and_predicted_together():
    ranked = [cand("a", 1.0, INC, "predicted"), cand("b", 0.95, INC, "measured"), cand("c", 0.93), cand("d", 0.5)]
    assert order(apply_ccs(ranked)) == ["c", "a", "d", "b"]
    # the only other window member is measured-inconsistent: nothing is ahead of the predicted one, so it stays first
    assert order(apply_ccs([cand("a", 1.0, INC, "predicted"), cand("b", 0.95, INC, "measured")])) == ["a", "b"]


def test_the_window_is_the_one_step_6_used_not_one_recomputed_after_the_demotion():
    # best base 1.0 -> window starts at 0.9. c (0.86) is outside it; had the window been taken from the best remaining candidate (b, 0.9) it would be inside
    out = apply_ccs([cand("a", 1.0, INC, "measured"), cand("b", 0.9), cand("c", 0.86, INC, "predicted")])
    assert order(out) == ["b", "c", "a"] and out[1].ccs["flag"] is None


def test_the_window_is_the_exact_one_even_when_the_rounded_base_score_sits_on_the_line():
    # best 0.999 -> the window starts at 0.8991. b's exact base score is 0.899096: outside, although base_score (rounded to 4 decimals) reads 0.8991.
    # Judged on the rounded value b would be a window member and the predicted-inconsistent a would go behind it, against the rule.
    def real(cid, norm, smiles):
        return Candidate(id=cid, layer="msfinder", smiles=smiles, name=cid, formula="C8H10N4O2", adduct="[M+H]+", layer_score=norm, layer_score_norm=norm, evidence={})
    ranked = rerank.rerank([real("a", 0.999, CAFFEINE), real("b", 0.899096, "CCO")], polarity="Positive", adduct="[M+H]+", omics="Either / unknown")
    assert [(c.id, c.base_score) for c in ranked] == [("a", 0.999), ("b", 0.8991)]
    ranked[0].ccs, ranked[1].ccs = {"verdict": INC, "source": "predicted"}, {"verdict": CON, "source": "predicted"}
    assert order(apply_ccs(ranked)) == ["a", "b"] and apply_ccs(ranked)[0].ccs["flag"] == "CCS ≠ predicted"          # a is alone in the window: flagged, not moved
    assert "_in_tiebreak_window" not in json.dumps(ranked[0].to_dict()) and "in_window" not in ranked[0].to_dict()   # the decision is not part of the result


@pytest.mark.parametrize("verdict, source", [(CON, "measured"), (CON, "predicted"), (NOREF, None), ("No structure", None), ("Not checked", None), (None, None)])
def test_everything_else_is_not_moved(verdict, source):
    ranked = [cand("a", 1.0, verdict, source), cand("b", 0.95, verdict, source), cand("c", 0.3, verdict, source)]
    out = apply_ccs(ranked)
    assert order(out) == ["a", "b", "c"] and all(c.ccs is None or (c.ccs["flag"], c.ccs["moved_from"]) == (None, None) for c in out)


def test_all_measured_inconsistent_keeps_the_order_and_still_flags():
    out = apply_ccs([cand("a", 1.0, INC, "measured"), cand("b", 0.5, INC, "measured")])
    assert order(out) == ["a", "b"] and [c.ccs["flag"] for c in out] == ["CCS ≠ measured"] * 2 and [c.ccs["moved_from"] for c in out] == [None, None]


def test_the_rule_is_pure_and_touches_no_score_or_level():
    ranked = [cand("a", 1.0, INC, "measured"), cand("b", 0.95)]
    before = copy.deepcopy(ranked)
    out = apply_ccs(ranked)
    assert [(c.id, c.rank, c.ccs) for c in ranked] == [(c.id, c.rank, c.ccs) for c in before]            # the input is untouched
    assert [(c.base_score, c.final_score, c.level, c.level_reason) for c in out] == [(c.base_score, c.final_score, c.level, c.level_reason) for c in reversed(before)]
    assert apply_ccs([]) == [] and order(apply_ccs([cand("only", 0.4, INC, "measured")])) == ["only"]


def test_a_wider_window_pulls_more_candidates_into_the_weak_rule():
    ranked = [cand("a", 1.0, INC, "predicted"), cand("b", 0.85), cand("c", 0.2)]
    assert order(apply_ccs(ranked)) == ["a", "b", "c"]                                                    # b is outside the 10 % window
    assert order(apply_ccs(ranked, rerank.RerankConfig(tiebreak_window=0.2))) == ["b", "a", "c"]


# ---------------------------------------------------------------- the form
def ccs_form(observed=140.5, platform=SINGLE_FIELD, **over):
    f = {"ccs_observed": observed, "ccs_platform": platform}
    f.update(over)
    return form(**f)


def test_no_observed_ccs_means_no_ccs_query_and_the_other_fields_are_ignored():
    assert pipeline.parse_inputs(form()).ccs is None
    assert pipeline.parse_inputs(form(ccs_observed="", ccs_platform="TIMS", ccs_gas="He")).ccs is None
    assert pipeline.parse_inputs(form(ccs_observed=None)).ccs is None


def test_the_ion_mobility_fields_are_parsed():
    q = pipeline.parse_inputs(ccs_form("152.3", "TIMS")).ccs
    assert (q.observed, q.platform, q.gas, q.tw_calibration) == (152.3, "TIMS", "N₂", ccs.TW_CALIBRATIONS[0])
    q = pipeline.parse_inputs(ccs_form(140, "TWIMS", ccs_gas="He", ccs_tw_calibration=ccs.TW_CALIBRATIONS[1])).ccs
    assert (q.observed, q.platform, q.gas, q.tw_calibration) == (140.0, "TWIMS", "He", ccs.TW_CALIBRATIONS[1])
    assert q.label == "140 Å² (TWIMS, He)"


@pytest.mark.parametrize("over, message", [
    ({"ccs_platform": ""}, r"Platform is required when an observed CCS is given"),
    ({"ccs_platform": None}, r"Platform is required"),
    ({"ccs_observed": 0.8}, r"too small .* not 1/K0"),                                   # a TIMS 1/K0
    ({"ccs_observed": 19.9}, r"too small"),
    ({"ccs_observed": 1500}, r"too large"),
    ({"ccs_observed": "abc"}, r"Observed CCS .* must be a number"),
    ({"ccs_observed": "152.3, 160.1"}, r"must be a number"),                              # one value only
    ({"ccs_observed": "nan"}, r"finite"),
    ({"ccs_platform": "Magic IMS"}, r"Platform must be one of"),
    ({"ccs_gas": "Ar"}, r"Drift gas must be one of"),
    ({"ccs_platform": "TWIMS", "ccs_tw_calibration": "other"}, r"TWIMS calibration must be one of"),
])
def test_parse_inputs_rejects_a_bad_ion_mobility_section(over, message):
    with pytest.raises(InputError, match=message):
        pipeline.parse_inputs(ccs_form(**over))


def test_the_twims_calibration_only_applies_to_twims():
    q = pipeline.parse_inputs(ccs_form(ccs_tw_calibration=ccs.TW_CALIBRATIONS[1]))
    assert q.ccs.tw_calibration == ccs.TW_CALIBRATIONS[0] and any("only applies to TWIMS" in w for w in q.warnings)
    assert not pipeline.parse_inputs(ccs_form()).warnings


def test_the_api_answers_400_for_a_bad_ion_mobility_section():
    r = client.post("/api/pipelines/esi_ms2_annotation/run", json=ccs_form(ccs_platform=""))
    assert r.status_code == 400 and "Platform is required" in r.json()["detail"]
    r = client.post("/api/pipelines/esi_ms2_annotation/run", json=ccs_form(0.9))
    assert r.status_code == 400 and "1/K0" in r.json()["detail"]


def test_the_spec_has_the_ion_mobility_fields_at_the_end_of_advanced_with_the_ccs_pages_options():
    fields = {f["id"]: f for f in ESI["inputs"]}
    advanced = [f["id"] for f in ESI["inputs"] if f.get("section") == "Advanced"]
    assert advanced == ["additive", "alkali", "ph", "lib_metric", "lib_threshold", "ccs_observed", "ccs_platform", "ccs_gas", "ccs_tw_calibration"]
    sections = list(dict.fromkeys(f["section"] for f in ESI["inputs"]))
    assert sections == ["Spectrum", "Search", "Options", "Advanced"] and "Ion Mobility" not in sections      # Advanced is closed by the page itself
    assert not any(f.get("section_collapsed") for f in ESI["inputs"])
    assert fields["ccs_observed"]["type"] == "number" and fields["ccs_observed"]["value"] is None and not fields["ccs_observed"].get("required")
    assert fields["ccs_platform"]["options"] == ["", *ccs.PLATFORMS] and fields["ccs_platform"]["value"] == "" and fields["ccs_platform"]["placeholder"]
    assert tuple(fields["ccs_gas"]["options"]) == ccs.GASES and fields["ccs_gas"]["value"] == "N₂"
    assert tuple(fields["ccs_tw_calibration"]["options"]) == ccs.TW_CALIBRATIONS
    assert fields["ccs_tw_calibration"]["visible_when"] == {"field": "ccs_platform", "equals": "TWIMS"}
    # one definition for the CCS page and the workflow
    page = {f["id"]: f for f in next(m for m in spec.MODULES if m["id"] == "ccs")["inputs"]}
    assert page["platform"]["options"] == list(ccs.PLATFORMS) == spec.CCS_PLATFORMS and tuple(spec.CCS_GASES) == ccs.GASES
    assert tuple(spec.CCS_TW_CALIBRATIONS) == ccs.TW_CALIBRATIONS
    assert page["gas"]["options"] == spec.CCS_GASES and page["tw_calibration"]["options"] == spec.CCS_TW_CALIBRATIONS
    assert ccs_step.FORM_KEYS == tuple(advanced[-4:])


def test_the_examples_do_not_carry_a_ccs():
    for example in ESI["examples"]:
        assert not set(example["values"]) & set(ccs_step.FORM_KEYS)
        assert pipeline.parse_inputs(example["values"]).ccs is None


# ---------------------------------------------------------------- step 6 with a fake CCS index
def run_step6(inputs):
    """pipeline.run_query on a parsed form -> (result, (status, note, blocks) step 6 ended with)."""
    emitted = {}
    result = pipeline.run_query(pipeline.parse_inputs(inputs), progress=lambda n, status, note, blocks: emitted.update({n: (status, note, blocks)}) if n == 6 and status in ("done", "failed") else None)
    return result, emitted[6]


def follow_the_top(rig):
    """The Rig's isomer check names a fixed id as top-1; a real one is asked about ranked[0]. Make the fake do the same."""
    rig.mp.setattr(pipeline.isomers, "discriminate", lambda pool, *a, **k: isomer_result("single", pool[0]["id"]))


def by_name(result):
    return {c["name"]: c for c in result["candidates"]}


def ranked_names(result):
    return [c["name"] for c in result["candidates"]]


def two_structures(rig, first=(ASPIRIN, "Aspirin", 0.9, "C9H8O4"), second=(CAFFEINE, "Caffeine", 0.85, "C8H10N4O2")):
    rig.msfinder([structure_candidate("msfinder", "m1", first[0], first[2], formula=first[3], name=first[1]),
                  structure_candidate("msfinder", "m2", second[0], second[2], formula=second[3], name=second[1])], good_enough=True)


def test_a_measured_inconsistent_top_goes_behind_a_consistent_one(rig, fake_ccs):
    two_structures(rig)
    plain, _ = run_step6(form())
    assert ranked_names(plain) == ["Aspirin", "Caffeine"]
    result, (status, note, blocks) = run_step6(ccs_form(140.5))
    assert ranked_names(result) == ["Caffeine", "Aspirin"] and [c["rank"] for c in result["candidates"]] == [1, 2]
    caf, asp = by_name(result)["Caffeine"]["ccs"], by_name(result)["Aspirin"]["ccs"]
    assert (caf["verdict"], caf["source"], caf["reference"], caf["tolerance_pct"], caf["ref_label"], caf["dataset"]) == (CON, "measured", 140.5, 1.0, "DTCCS N2", "fakeDT")
    assert caf["delta_pct"] == 0 and caf["flag"] is None and caf["moved_from"] is None and caf["observed"] == 140.5 and caf["adduct"] == "[M+H]+"
    assert (asp["verdict"], asp["source"], asp["reference"], asp["tolerance_pct"], asp["ref_label"], asp["dataset"]) == (INC, "measured", 134.7, 2.0, "TWCCS N2", "fakeTW")
    assert round(asp["delta_pct"], 2) == 4.31 and asp["flag"] == "CCS ≠ measured" and asp["moved_from"] == 1
    # the scores and levels are what they were; only the position changed
    for name in ("Aspirin", "Caffeine"):
        old, new = by_name(plain)[name], by_name(result)[name]
        assert (old["base_score"], old["final_score"], old["level"], old["level_reason"]) == (new["base_score"], new["final_score"], new["level"], new["level_reason"])
    # the isomer check and the answer follow the new order (the Rig's fake isomer check names c1 as its top, so look at what it was asked)
    assert rig.calls["isomers"][-1][0][0][0]["id"] == by_name(result)["Caffeine"]["id"]
    # the step says what happened
    assert status == "done" and "Observed CCS 140.5 Å² (DTIMS (single field), N₂): 1 consistent, 1 inconsistent with a measured value. 1 moved by CCS." in note
    assert any("CCS: Aspirin is +4.3 % from the measured 134.7 Å² (±2.0 %); it was moved behind the other candidates, so Caffeine is now ranked first." == w
               for w in result["warnings"])


def test_the_ccs_check_block_lists_every_candidate_in_the_final_order(rig, fake_ccs):
    two_structures(rig)
    _, (_, _, blocks) = run_step6(ccs_form(140.5))
    assert [b["id"] for b in blocks] == ["ranked_table", "ranked_note", "ccs_table", "ccs_note"]
    table = next(b for b in blocks if b["id"] == "ccs_table")
    assert table["title"] == "CCS check" and table["columns"] == ["Rank", "Candidate", "Adduct", "Verdict", "Δ %", "Tolerance %", "Source"]
    rows = table["data"]["rows"]
    assert [r[0] for r in rows] == [1, 2] and [r[1] for r in rows] == ["Caffeine", "Aspirin"] and [r[2] for r in rows] == ["[M+H]+"] * 2
    assert rows[0][3]["badge"] == CON and rows[0][3]["tone"] == "good" and rows[1][3]["badge"] == INC and rows[1][3]["tone"] == "warn"
    assert "rank 1 to 2" in rows[1][3]["help"]
    assert (rows[0][4], rows[0][5], rows[0][6]) == ("+0.00", "±1.0", "Measured 140.5 Å² · DTCCS N2 · fakeDT")
    assert (rows[1][4], rows[1][5]) == ("+4.31", "±2.0")
    note = next(b for b in blocks if b["id"] == "ccs_note")["data"]["items"]
    assert note == ["Moved by CCS: Aspirin: rank 1 to 2 (+4.3 % from the measured 134.7 Å² (±2.0 %))."]
    assert "Levels do not change" in table["caption"] and "140.5 Å²" in table["caption"]
    ranked = next(b for b in blocks if b["id"] == "ranked_table")["data"]["rows"]
    assert [r[1] for r in ranked] == ["Caffeine", "Aspirin"], "the ranked table shows the final order"


def test_a_predicted_inconsistent_candidate_only_moves_inside_the_window(rig, fake_model):
    # the fake model predicts 150.0 +/- 5.6 % for every compound the index has no value for; 140.5 is 6.3 % below it
    two_structures(rig, first=(IBUPROFEN, "Ibuprofen", 0.9, "C13H18O2"))
    result, (_, note, _) = run_step6(ccs_form(140.5))
    assert ranked_names(result) == ["Caffeine", "Ibuprofen"]
    ibu = by_name(result)["Ibuprofen"]["ccs"]
    assert (ibu["verdict"], ibu["source"], ibu["reference"], ibu["tolerance_pct"], ibu["flag"], ibu["moved_from"]) == (INC, "predicted", 150.0, 5.6, "CCS ≠ predicted", 1)
    assert ibu["ref_label"] == "SVR" and ibu["in_domain"] is not None and ibu["dataset"] is None
    assert "1 inconsistent with a predicted value" in note
    assert any("moved behind the other near-tied candidates" in w for w in result["warnings"])


def test_a_predicted_inconsistent_candidate_outside_the_window_stays(rig, fake_model):
    two_structures(rig, first=(CAFFEINE, "Caffeine", 0.9, "C8H10N4O2"), second=(IBUPROFEN, "Ibuprofen", 0.5, "C13H18O2"))
    result, _ = run_step6(ccs_form(140.5))
    assert ranked_names(result) == ["Caffeine", "Ibuprofen"]
    ibu = by_name(result)["Ibuprofen"]["ccs"]
    assert ibu["verdict"] == INC and ibu["flag"] is None and ibu["moved_from"] is None
    assert not any("CCS:" in w for w in result["warnings"])


def test_a_measured_inconsistent_candidate_is_demoted_even_behind_a_candidate_with_no_reference(rig, fake_ccs):
    # no model installed: the second candidate has no reference at all, and still goes ahead of the one a measured value contradicts
    two_structures(rig, first=(CAFFEINE, "Caffeine", 0.9, "C8H10N4O2"), second=(IBUPROFEN, "Ibuprofen", 0.5, "C13H18O2"))
    result, _ = run_step6(ccs_form(160.0))
    assert ranked_names(result) == ["Ibuprofen", "Caffeine"]
    ibu = by_name(result)["Ibuprofen"]["ccs"]
    assert ibu["verdict"] == NOREF and ibu["source"] is None and ibu["reference"] is None and "No prediction model" in ibu["note"] and ibu["flag"] is None
    assert by_name(result)["Caffeine"]["ccs"]["flag"] == "CCS ≠ measured"


def test_levels_never_change_even_when_a_library_match_is_demoted(rig, fake_ccs):
    rig.library([library_candidate("lib1", CAFFEINE, "Caffeine", 0.93, True), library_candidate("lib2", ETHER, "Dimethyl ether", 0.72, False, formula="C2H6O")], passed=True)
    follow_the_top(rig)
    plain, _ = run_step6(form(libraries=["MassBank"]))
    result, _ = run_step6(ccs_form(170.0, libraries=["MassBank"]))
    assert ranked_names(plain) == ["Caffeine", "Dimethyl ether"] and ranked_names(result) == ["Dimethyl ether", "Caffeine"]
    assert {c["name"]: c["level"] for c in result["candidates"]} == {"Caffeine": "2a", "Dimethyl ether": "3"} == {c["name"]: c["level"] for c in plain["candidates"]}
    assert result["answer"]["level"] == "3"                                           # the answer is the first-ranked candidate, as always
    assert any("re-ranking put a candidate that did not pass the gate first" in w for w in result["warnings"])
    assert any(w.startswith("CCS: Caffeine is ") and "Dimethyl ether is now ranked first" in w for w in result["warnings"])


@pytest.mark.parametrize("platform, gas, tw, phrase", [
    (SINGLE_FIELD, "He", ccs.TW_CALIBRATIONS[0], "measured in He"),
    ("FAIMS / DMS", "N₂", ccs.TW_CALIBRATIONS[0], "FAIMS / DMS gives a compensation voltage"),
    ("TWIMS", "N₂", ccs.TW_CALIBRATIONS[1], "TWCCS N₂→He"),
])
def test_a_ccs_that_cannot_be_compared_is_noted_once_and_changes_nothing(rig, fake_ccs, platform, gas, tw, phrase):
    two_structures(rig)
    plain, _ = run_step6(form())
    result, (status, note, blocks) = run_step6(ccs_form(140.5, platform, ccs_gas=gas, ccs_tw_calibration=tw))
    assert status == "done" and [c["name"] for c in result["candidates"]] == ranked_names(plain)
    assert all("ccs" not in c for c in result["candidates"]) and [b["id"] for b in blocks] == ["ranked_table", "ranked_note"]
    assert note.endswith(f"Observed CCS not used ({ccs.comparability(platform, gas, tw)[1]}).")
    notes = [w for w in result["warnings"] if "observed CCS" in w]
    assert len(notes) == 1 and phrase in notes[0] and notes[0].startswith("The observed CCS (140.5 Å²")
    assert [(c["rank"], c["base_score"], c["final_score"], c["level"]) for c in result["candidates"]] == [(c["rank"], c["base_score"], c["final_score"], c["level"]) for c in plain["candidates"]]


def test_a_missing_index_is_noted_with_the_setup_command_and_changes_nothing(rig, tmp_path, monkeypatch):
    monkeypatch.setenv("MZDESK_CCS_DIR", str(tmp_path))
    two_structures(rig)
    plain, _ = run_step6(form())
    result, (_, note, blocks) = run_step6(ccs_form(140.5))
    assert ranked_names(result) == ranked_names(plain) and all("ccs" not in c for c in result["candidates"]) and len(blocks) == 2
    assert "Observed CCS not used (CCS data not installed)." in note
    [warning] = [w for w in result["warnings"] if "observed CCS" in w]
    assert "python -m mzdesk.orthogonal.ccs_setup" in warning and "not installed" in warning


def test_an_index_in_an_old_format_is_noted_and_changes_nothing(rig, fake_ccs):
    (fake_ccs / ccs.META_FILE).touch()
    ccs._cache.clear()
    import gzip, json
    with gzip.open(fake_ccs / ccs.INDEX_FILE, "wt") as f:
        json.dump({"schema": 0, "columns": [], "rows": [], "sources": {}}, f)
    two_structures(rig)
    result, (_, note, _) = run_step6(ccs_form(140.5))
    assert ranked_names(result) == ["Aspirin", "Caffeine"] and all("ccs" not in c for c in result["candidates"])
    assert "Observed CCS not used (CCS index needs a rebuild)." in note and any("old format" in w and "--force" in w for w in result["warnings"])


def test_no_structure_formula_only_rows_and_unreadable_adducts(rig, fake_ccs):
    # a formula-only result: every row says there is no structure, the step note says the CCS was not used
    result, (_, note, blocks) = run_step6(ccs_form(140.5))
    assert result["answer"]["type"] == "formula_only" and blocks == []
    assert [c["ccs"]["verdict"] for c in result["candidates"]] == ["No structure"] and "The observed CCS was not used: there is no structure to look up." in note
    # check() itself: a candidate without a structure, one with an adduct the engine cannot read
    q = ccs_step.CcsQuery(140.5, SINGLE_FIELD)
    bare = Candidate(id="x1", layer="library", name="1 M+H", adduct="[M+H]+")
    odd = structure_candidate("msfinder", "x2", CAFFEINE, 0.5)
    odd.adduct = "[M+Xx]+"
    noadd = structure_candidate("msfinder", "x3", CAFFEINE, 0.5)
    noadd.adduct = None
    out = ccs_step.check([bare, odd, noadd], q, None)
    assert out.checks["x1"]["verdict"] == "No structure"
    assert out.checks["x2"]["verdict"] == "Not checked" and "Unknown adduct" in out.checks["x2"]["note"]
    assert out.checks["x3"]["verdict"] == "Not checked" and out.checks["x3"]["note"] == "no adduct to look up"
    assert ccs_step.check([noadd], q, "[M+H]+").checks["x3"]["verdict"] == CON            # the run's adduct stands in for a missing one
    assert out.counts == {"No structure": 1, "Not checked": 2}


def test_the_candidates_own_adduct_is_the_one_looked_up(rig, fake_ccs):
    rig.msfinder([structure_candidate("msfinder", "m1", CAFFEINE, 0.9, adduct="[M+Na]+", name="Caffeine")], good_enough=True)
    # the fake index has caffeine [M+Na]+ at 142.9 (DT) only
    result, _ = run_step6(ccs_form(142.9))
    x = result["candidates"][0]["ccs"]
    assert (x["adduct"], x["reference"], x["verdict"]) == ("[M+Na]+", 142.9, CON)


def test_a_failing_ccs_check_leaves_the_ranking_alone(rig, fake_ccs, monkeypatch):
    two_structures(rig)

    def broken(*a, **k):
        raise RuntimeError("index exploded")
    monkeypatch.setattr(ccs_step, "check", broken)
    result, (status, note, blocks) = run_step6(ccs_form(140.5))
    assert status == "done" and ranked_names(result) == ["Aspirin", "Caffeine"] and all("ccs" not in c for c in result["candidates"]) and len(blocks) == 2
    assert "The CCS check failed (RuntimeError: index exploded); the order is the one without CCS." in note
    assert any("CCS check failed" in w for w in result["warnings"])


def test_a_failed_rerank_skips_the_ccs_check(rig, fake_ccs, monkeypatch):
    two_structures(rig)
    monkeypatch.setattr(pipeline, "rerank", lambda *a, **k: (_ for _ in ()).throw(ValueError("rules broke")))
    result, (status, note, _) = run_step6(ccs_form(140.5))
    assert status == "failed" and "The observed CCS was not used." in note and all("ccs" not in c for c in result["candidates"])


def test_the_progress_calls_count_the_ccs_checks(rig, fake_ccs):
    two_structures(rig)
    seen = []
    pipeline.run_query(pipeline.parse_inputs(ccs_form(140.5)), progress=lambda n, status, note, blocks: seen.append(note) if n == 6 else None)
    assert "CCS checked 1/2" in seen and "CCS checked 2/2" in seen


def test_the_summary_and_description_words():
    q = ccs_step.CcsQuery(152.3, "TIMS")
    out = ccs_step.Outcome(checks={}, counts={"Consistent": 2, "Inconsistent measured": 1, "Inconsistent predicted": 3, "No reference": 4, "No structure": 1, "Not checked": 1})
    assert ccs_step.summary(q, out, 2) == ("Observed CCS 152.3 Å² (TIMS, N₂): 2 consistent, 1 inconsistent with a measured value, 3 inconsistent with a predicted value, "
                                           "4 without a reference, 2 not checked. 2 moved by CCS.")
    assert ccs_step.summary(q, ccs_step.Outcome(checks={}, counts={}), 0) == "Observed CCS 152.3 Å² (TIMS, N₂): no candidate checked."


# ---------------------------------------------------------------- through the job API
def test_the_job_api_carries_the_ccs_check(rig, fake_ccs):
    import time
    two_structures(rig)
    job_id = client.post("/api/pipelines/esi_ms2_annotation/run", json=ccs_form(140.5)).json()["job_id"]
    for _ in range(200):
        job = client.get(f"/api/pipelines/jobs/{job_id}").json()
        if job["status"] not in ("queued", "running"):
            break
        time.sleep(0.05)
    assert job["status"] == "done"
    step6 = next(s for s in job["steps"] if s["id"] == 6)
    assert [b["id"] for b in step6["blocks"]] == ["ranked_table", "ranked_note", "ccs_table", "ccs_note"] and "Observed CCS 140.5" in step6["note"]
    assert [c["name"] for c in job["result"]["candidates"]] == ["Caffeine", "Aspirin"] and job["result"]["candidates"][1]["ccs"]["flag"] == "CCS ≠ measured"


# ---------------------------------------------------------------- the real CCS index
def _caffeine_reference():
    """Caffeine's [M+H]+ values in the local index, N2 scale: (the DT value used as the observed CCS, all values)."""
    index = ccs.load()
    rows, _ = index.lookup(ci.normalize(CAFFEINE).key14, "[M+H]+")
    assert rows, "caffeine [M+H]+ is in the CCS index"
    dt = next(r for r in rows if r.ccs_type == "DT")
    return dt.ccs, [r.ccs for r in rows]


def _needs_predictor():
    predictor, why = ccs.load_predictor(None, ccs.load())
    if predictor is None:
        pytest.skip(f"no prediction model: {why}")


def caffeine_and_isocaffeine(rig):
    """The library example's pool in a hand-made form: caffeine (measured values in the index) and isocaffeine (none: predicted) within 10 % of each other, and an outsider."""
    rig.library([library_candidate("lib1", CAFFEINE, "Caffeine", 0.93, True), library_candidate("lib2", ISOCAFFEINE, "Isocaffeine", 0.90, True),
                 library_candidate("lib3", ETHER, "Dimethyl ether", 0.72, False, formula="C2H6O")], passed=True)
    follow_the_top(rig)


@needs_real_ccs
def test_real_index_caffeines_own_value_is_consistent(rig):
    _needs_predictor()
    observed, values = _caffeine_reference()
    caffeine_and_isocaffeine(rig)
    result, (_, note, _) = run_step6(ccs_form(observed, libraries=["MassBank"]))
    print(f"\ncaffeine [M+H]+ in the index: {sorted(values)}; observed {observed}")
    print(" ", note)
    caf = by_name(result)["Caffeine"]["ccs"]
    assert (caf["verdict"], caf["source"], caf["delta_pct"], caf["flag"]) == (CON, "measured", 0.0, None) and caf["n_refs"] == len(values)
    iso = by_name(result)["Isocaffeine"]["ccs"]
    assert iso["source"] == "predicted" and iso["verdict"] in (CON, INC)
    assert ranked_names(result)[:2] == ["Caffeine", "Isocaffeine"] or iso["verdict"] == INC


@needs_real_ccs
def test_real_index_ten_percent_above_caffeines_value_is_inconsistent_and_moves_caffeine_back(rig):
    _needs_predictor()
    observed, values = _caffeine_reference()
    high = round(observed * 1.10, 2)
    caffeine_and_isocaffeine(rig)
    plain, _ = run_step6(form(libraries=["MassBank"]))
    result, (_, note, blocks) = run_step6(ccs_form(high, libraries=["MassBank"]))
    caf, iso = by_name(result)["Caffeine"]["ccs"], by_name(result)["Isocaffeine"]["ccs"]
    print(f"\ncaffeine [M+H]+ measured {sorted(values)}; observed {observed} +10 % = {high}")
    print("  caffeine:", {k: caf[k] for k in ("verdict", "source", "reference", "delta_pct", "tolerance_pct", "ref_label", "dataset", "flag", "moved_from")})
    print("  isocaffeine:", {k: iso[k] for k in ("verdict", "source", "reference", "delta_pct", "tolerance_pct", "ref_label", "flag", "moved_from")})
    print(" ", note)
    print("  order:", ranked_names(plain), "->", ranked_names(result))
    assert ranked_names(plain) == ["Caffeine", "Isocaffeine", "Dimethyl ether"]
    assert (caf["verdict"], caf["source"], caf["flag"], caf["moved_from"]) == (INC, "measured", "CCS ≠ measured", 1) and caf["delta_pct"] > caf["tolerance_pct"]
    assert iso["source"] == "predicted" and iso["verdict"] == INC and iso["flag"] == "CCS ≠ predicted"   # +10 % is far from the predicted ~139 ...
    assert iso["moved_from"] is None                                                                       # ... but it only gained a place: caffeine left, nothing overtook it
    assert ranked_names(result) == ["Isocaffeine", "Dimethyl ether", "Caffeine"]                           # measured demotes globally, the predicted one is first in the window
    assert by_name(result)["Caffeine"]["level"] == "2a" and result["answer"]["top_id"] == by_name(result)["Isocaffeine"]["id"]
    assert next(b for b in blocks if b["id"] == "ccs_table")["data"]["rows"][-1][1] == "Caffeine"


@needs_real_ccs
def test_real_index_he_and_faims_change_nothing(rig):
    observed, _ = _caffeine_reference()
    caffeine_and_isocaffeine(rig)
    plain, _ = run_step6(form(libraries=["MassBank"]))
    for kw in ({"ccs_gas": "He"}, {"ccs_platform": "FAIMS / DMS"}):
        result, _ = run_step6(ccs_form(observed * 1.1, libraries=["MassBank"], **kw))
        assert ranked_names(result) == ranked_names(plain) and all("ccs" not in c for c in result["candidates"])


# ---------------------------------------------------------------- real engines (opt-in, like tests/test_pipeline.py)
@pytest.mark.skipif(os.environ.get("MZDESK_RUN_PIPELINE") != "1" or not ccs.installed(), reason="real runs: set MZDESK_RUN_PIPELINE=1 and install the CCS index")
def test_real_library_example_with_caffeines_own_ccs_and_ten_percent_above():
    example = next(e for e in ESI["examples"] if e["label"].endswith("library match"))["values"]
    observed, _ = _caffeine_reference()
    consistent = pipeline.run({**example, "ccs_observed": observed, "ccs_platform": SINGLE_FIELD})
    shifted = pipeline.run({**example, "ccs_observed": round(observed * 1.1, 2), "ccs_platform": SINGLE_FIELD})
    caf = by_name(consistent)["Caffeine"]["ccs"]
    assert caf["verdict"] == CON and ranked_names(consistent)[0] == "Caffeine"
    caf = by_name(shifted)["Caffeine"]["ccs"]
    assert caf["verdict"] == INC and caf["flag"] == "CCS ≠ measured" and ranked_names(shifted)[-1] == "Caffeine"
