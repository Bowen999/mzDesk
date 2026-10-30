"""ESI MS/MS Annotation workflow: the orchestrator (mzdesk/annotation/pipeline.py) and the job API (/api/pipelines/...).

Unit tests replace the layer functions (BUDDY, library, MS-FINDER, FRIGID, the isomer check) with hand-made answers and cover every branch of the
flow; nothing starts an engine. The re-rank is replaced by a fast stand-in except in one test (the real one takes about a second per structure).
Opt-in real end-to-end runs (BUDDY, local MassBank / GNPS, MS-FINDER, FRIGID; one at a time, minutes): MZDESK_RUN_PIPELINE=1.
"""
import os
import threading
import time

import pytest
from fastapi.testclient import TestClient

from mzdesk.annotation import frigid_step, pipeline
from mzdesk.annotation.adduct_step import AdductHypothesis, AdductResult
from mzdesk.annotation.formula_step import FormulaCandidate, FormulaResult
from mzdesk.annotation.library_step import GateSummary
from mzdesk.annotation.schema import Candidate, LayerResult
from mzdesk.api import app as app_module
from mzdesk.api.common import InputError
from mzdesk.ms2 import library_index as lx
from mzdesk.spec import CAFFEINE_MASSBANK_MS2, CAFFEINE_MS1, CAFFEINE_MSFINDER_MS2, PIPELINES

client = TestClient(app_module.app)
REAL_FRIGID, REAL_RERANK = frigid_step.run_frigid_layer, pipeline.rerank           # the originals, before a Rig replaces them
CAFFEINE, ETHANOL, ETHER, PROPANOL = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "CCO", "COC", "CCCO"       # ethanol and dimethyl ether are isomers


def form(**over):
    """The pipeline form as the page sends it (a positive-mode caffeine spectrum, Auto adduct, no library)."""
    base = {"precursor": {"precursor_mz": "195.0877", "ion_mode": "Positive", "adduct": "Auto"}, "ms2": {"source": "Paste", "value": CAFFEINE_MSFINDER_MS2},
            "ms1": "", "libraries": [], "tol_ppm": 10, "tol_da": 0.02, "omics": "Either / unknown", "halogen": False, "local_isomers": False}
    return {**base, **over}


def form_at(mz, **over):
    return form(precursor={"precursor_mz": str(mz), "ion_mode": "Positive", "adduct": "Auto"}, **over)


def negative(**over):
    return form(precursor={"precursor_mz": "193.0731", "ion_mode": "Negative", "adduct": "Auto"}, **over)


# ---------------------------------------------------------------- hand-made layer answers
def adduct_result(*names):
    hyps = [AdductHypothesis(adduct=n, score=round(1 / (i + 1), 3), rank=i + 1) for i, n in enumerate(names)]
    return AdductResult(mode="auto", hypotheses=hyps, note="adducts scored", seconds=0.01)


def formula_result(*pairs, status="done"):
    formulas = [FormulaCandidate(formula=f, adduct=a, rank=i + 1, prob=0.9, ppm=0.3) for i, (f, a) in enumerate(pairs)]
    return FormulaResult(formulas=formulas, adducts=sorted({a for _, a in pairs}), status=status, note=f"{len(formulas)} formulas", seconds=0.02)


def library_candidate(cid, smiles, name, similarity, passed, formula="C8H10N4O2"):
    ev = {"gate_passed": passed, "failed_checks": [] if passed else ["Similarity: 0.76 is below 0.80"], "similarity": similarity, "matched_peaks": 9,
          "precursor_error_ppm": 0.3}
    return Candidate(id=cid, layer="library", smiles=smiles, name=name, formula=formula, adduct="[M+H]+", layer_score=similarity,
                     layer_score_norm=similarity, evidence={"library": ev})


def structure_candidate(layer, cid, smiles, score_norm, formula="C8H10N4O2", adduct="[M+H]+", name=""):
    return Candidate(id=cid, layer=layer, smiles=smiles, name=name, formula=formula, adduct=adduct, layer_score=score_norm * 10,
                     layer_score_norm=score_norm, evidence={layer: {"score": score_norm}}, hypothesis=layer == "frigid")


def isomer_result(verdict, top_id, *, group=None, challengers=(), notes=("Relative to the 2 candidates considered.",), error=None):
    ids = {"single": [top_id], "indistinguishable": group or [top_id]}.get(verdict, [])
    return {"verdict": verdict, "top_id": top_id, "group": ids, "engine": "MS-FINDER", "log": "engine log", "error": error, "notes": list(notes),
            "isomers": [{"id": i, "name": i, "smiles": ETHER, "evaluated": True, "distinguished": verdict == "single",
                         "distinguishing_peaks": [{"mz": 45.0334, "rel": 30.0, "formula": "C2H5O+", "smiles": "CCO"}]}
                        for i in ("c2",)] if verdict != "insufficient" else [],
            "challengers": [{"id": i, "name": f"isomer {i}", "peaks": [{"mz": 45.0334, "rel": 30.0}]} for i in challengers],
            "candidates": [{"id": top_id, "name": "top", "smiles": ETHANOL, "source": "input", "explained_fraction": 0.8},
                           {"id": "c2", "name": "other", "smiles": ETHER, "source": "input", "explained_fraction": 0.4}] if verdict != "insufficient" else []}


def fast_rerank(candidates, *, progress=None, **kw):
    """Stand-in for rerank.rerank: sorted by the layer score, no detectability rules."""
    out = sorted(candidates, key=lambda c: -c.layer_score_norm)
    for i, c in enumerate(out, start=1):
        c.rank, c.base_score, c.final_score = i, c.layer_score_norm, c.layer_score_norm
        c.detect = {"omics": "both", "ionizable": "strong", "adduct": "high", "factor": 1.0}
        if progress:
            progress(i, len(out))
    return out


class Rig:
    """Replaces the layer functions of pipeline.py and records how each was called."""

    def __init__(self, monkeypatch):
        self.mp, self.calls = monkeypatch, {}
        self.steps12(("[M+H]+",), (("C8H10N4O2", "[M+H]+"),))
        self.library(None)
        self.msfinder([])
        self.frigid([])
        self.isomers(isomer_result("insufficient", "c1", notes=("No isomer to compare with.",)))
        monkeypatch.setattr(pipeline, "rerank", fast_rerank)

    def record(self, name, *args, **kw):
        self.calls.setdefault(name, []).append((args, kw))

    def steps12(self, adducts, formulas, formula_status="done"):
        def fake(*args, **kw):
            self.record("steps12", *args, **kw)
            return adduct_result(*adducts), formula_result(*formulas, status=formula_status)
        self.mp.setattr(pipeline.formula_step, "run_adduct_and_formula", fake)

    def library(self, candidates, passed=False, status="done"):
        """candidates None: the layer must not be called."""
        def fake(*args, **kw):
            self.record("library", *args, **kw)
            assert candidates is not None, "the library step must not run"
            gate = GateSummary(passed=passed, best_id="lib1" if passed else None, checks=[{"label": "Similarity", "status": "pass", "detail": "0.93"}],
                               warnings=["a library warning"])
            return LayerResult("library", status, "library note", 0.1, list(candidates)), gate
        self.mp.setattr(pipeline.library_step, "run_library", fake)

    def msfinder(self, candidates, good_enough=False, error=None):
        def fake(*args, **kw):
            self.record("msfinder", *args, **kw)
            if error:
                raise error
            return LayerResult("msfinder", "done", "msfinder note", 0.1, list(candidates)), {"good_enough": good_enough, "warnings": ["an MS-FINDER warning"]}
        self.mp.setattr(pipeline.msfinder_step, "run_msfinder_layer", fake)

    def frigid(self, candidates, status="done", error=None):
        def fake(*args, **kw):
            self.record("frigid", *args, **kw)
            if error:
                raise error
            return LayerResult("frigid", status, "frigid note", 0.1, list(candidates))
        self.mp.setattr(pipeline.frigid_step, "run_frigid_layer", fake)

    def isomers(self, result, error=None):
        def fake(*args, **kw):
            self.record("isomers", *args, **kw)
            if error:
                raise error
            return result
        self.mp.setattr(pipeline.isomers, "discriminate", fake)


@pytest.fixture
def rig(monkeypatch):
    return Rig(monkeypatch)


# ---------------------------------------------------------------- the branches
def test_gate_passed_skips_steps_4_and_5_but_runs_rerank_and_isomers(rig):
    caffeine = library_candidate("lib1", CAFFEINE, "Caffeine", 0.93, True)
    weak = library_candidate("lib2", ETHER, "Dimethyl ether", 0.72, False, formula="C2H6O")
    rig.library([caffeine, weak], passed=True)
    rig.isomers(isomer_result("single", "c1"))
    result = pipeline.run(form(libraries=["MassBank"]))
    assert [s["status"] for s in result["steps"]] == ["done", "done", "done", "skipped", "skipped", "done", "done", "done"]
    assert "msfinder" not in rig.calls and "frigid" not in rig.calls
    assert [c["name"] for c in result["candidates"]] == ["Caffeine", "Dimethyl ether"]                   # the sub-threshold hit stays in the pool
    assert [c["level"] for c in result["candidates"]] == ["2a", "3"]
    assert result["answer"]["type"] == "single" and result["answer"]["level"] == "2a" and result["answer"]["layer"] == "library"
    assert result["candidates"][0]["isomer"]["role"] == "top" and result["context"]["gate_passed"]
    assert "a library warning" in result["warnings"]
    # the pool for step 7: the top-1 and the same-formula candidates only (dimethyl ether is C2H6O)
    pool = rig.calls["isomers"][0][0][0]
    assert [c["id"] for c in pool] == ["c1"] and rig.calls["isomers"][0][1]["local_isomers"] is False


def test_gate_failed_msfinder_good_enough_skips_frigid(rig):
    rig.library([library_candidate("lib1", CAFFEINE, "Caffeine", 0.75, False)])
    rig.msfinder([structure_candidate("msfinder", "m1", CAFFEINE, 0.8), structure_candidate("msfinder", "m2", ETHER, 0.3, formula="C2H6O")], good_enough=True)
    result = pipeline.run(form(libraries=["MassBank"]))
    assert [s["status"] for s in result["steps"]][2:5] == ["done", "done", "skipped"]
    assert "good enough" in result["steps"][4]["note"] and "frigid" not in rig.calls
    top = result["candidates"][0]
    assert top["name"] == "Caffeine" and top["sources"] and {s["layer"] for s in top["sources"]} == {"library", "msfinder"}      # merged by structure
    assert top["level"] == "3" and result["answer"]["level"] == "3"
    # MS-FINDER got the top adduct, the formula objects and the tolerances
    (args, kw) = rig.calls["msfinder"][0]
    assert args[1] == "[M+H]+" and args[3][0].formula == "C8H10N4O2" and kw["tol_ppm"] == 10 and kw["tol_da"] == 0.02


def test_gate_failed_weak_msfinder_runs_frigid_on_the_top_formula(rig):
    rig.msfinder([structure_candidate("msfinder", "m1", CAFFEINE, 0.3)], good_enough=False)
    rig.frigid([structure_candidate("frigid", "f1", ETHER, 0.2, formula="C8H10N4O2")])
    result = pipeline.run(form())
    assert [s["status"] for s in result["steps"]][2:5] == ["skipped", "done", "done"]              # no library selected
    args, kw = rig.calls["frigid"][0]
    assert args[1] == "[M+H]+" and args[3] == "C8H10N4O2" and kw["top_n"] == 5
    levels = {c["name"] or c["smiles"]: c["level"] for c in result["candidates"]}
    assert levels == {CAFFEINE: "3", ETHER: "H"}
    assert result["candidates"][0]["level"] == "3" and result["answer"]["layer"] == "msfinder" and result["answer"]["level"] == "3"
    assert next(c for c in result["candidates"] if c["level"] == "H")["hypothesis"] is True


def test_negative_mode_skips_frigid_with_its_own_reason(monkeypatch, rig):
    rig.steps12(("[M-H]-",), (("C8H10N4O2", "[M-H]-"),))
    rig.msfinder([structure_candidate("msfinder", "m1", CAFFEINE, 0.3, adduct="[M-H]-")])
    monkeypatch.setattr(pipeline.frigid_step, "run_frigid_layer", REAL_FRIGID)       # the real step: negative ions end before any engine call
    result = pipeline.run(negative())
    frigid = result["steps"][4]
    assert frigid["status"] == "skipped" and "positive ions only" in frigid["note"]
    assert result["answer"]["type"] == "single" and result["context"]["polarity"] == "Negative"


def test_a_layer_that_raises_fails_its_step_and_the_run_goes_on(rig):
    rig.msfinder([], error=RuntimeError("mono crashed"))
    rig.frigid([structure_candidate("frigid", "f1", CAFFEINE, 0.8)])
    result = pipeline.run(form())
    ms = result["steps"][3]
    assert ms["status"] == "failed" and "RuntimeError: mono crashed" in ms["note"]
    assert result["steps"][4]["status"] == "done"                                          # a failed MS-FINDER is not "good enough"
    assert any("Step 4" in w and "mono crashed" in w for w in result["warnings"])
    assert result["answer"]["type"] == "single" and result["answer"]["level"] == "H"


def test_library_that_raises_and_frigid_that_raises(rig):
    rig.library(None)
    rig.mp.setattr(pipeline.library_step, "run_library", lambda *a, **k: (_ for _ in ()).throw(ValueError("index broken")))
    rig.frigid([], error=OSError("no gpu"))
    result = pipeline.run(form(libraries=["GNPS"]))
    assert [s["status"] for s in result["steps"]][2:5] == ["failed", "done", "failed"]
    assert result["answer"]["type"] == "formula_only"                                       # nothing structural, but BUDDY's formula is there


def test_layer_returning_failed_status_is_reported_as_failed(rig):
    rig.library([], status="failed")
    result = pipeline.run(form(libraries=["MassBank"]))
    assert result["steps"][2]["status"] == "failed" and any("Step 3" in w for w in result["warnings"])


def test_no_library_selected_skips_step_3_without_calling_it(rig):
    rig.msfinder([structure_candidate("msfinder", "m1", CAFFEINE, 0.8)], good_enough=True)
    result = pipeline.run(form(libraries=[]))
    assert result["steps"][2]["status"] == "skipped" and "No library was selected" in result["steps"][2]["note"] and "library" not in rig.calls
    assert result["context"]["gate_passed"] is False and result["answer"]["layer"] == "msfinder"


def test_empty_formula_list_goes_to_msfinder_and_frigid_without_one(rig):
    rig.steps12(("[M+H]+",), (), )
    rig.msfinder([], good_enough=False)
    result = pipeline.run(form())
    assert rig.calls["msfinder"][0][0][3] == []                                              # MS-FINDER then searches its own formulas
    assert rig.calls["frigid"][0][0][3] is None                                              # FRIGID skips itself on "no formula"
    assert result["answer"]["type"] == "mass_only" and result["answer"]["level"] == "5"
    assert result["steps"][6]["status"] == "skipped"


def test_no_structure_but_formulas_gives_a_formula_answer(rig):
    rig.steps12(("[M+H]+", "[M+Na]+"), (("C8H10N4O2", "[M+H]+"), ("C7H6N6O2", "[M+H]+")))
    result = pipeline.run(form())
    assert result["answer"]["type"] == "formula_only" and result["answer"]["level"] == "4"
    assert [c["formula"] for c in result["candidates"]] == ["C8H10N4O2", "C7H6N6O2"] and {c["layer"] for c in result["candidates"]} == {"formula_only"}
    assert result["answer"]["ids"] == [c["id"] for c in result["candidates"]]
    assert result["steps"][5]["status"] == "done" and "formulas are listed" in result["steps"][5]["note"] and result["steps"][6]["status"] == "skipped"
    assert result["isomer"] is None and "isomers" not in rig.calls
    assert result["context"]["retained_adducts"] == ["[M+H]+", "[M+Na]+"]


def test_formula_table_rank_matches_result_order_and_preserves_buddy_rank_per_adduct(rig):
    formulas = FormulaResult(formulas=[
        FormulaCandidate("C10H20N2O5", "[M+H]+", 1, prob=0.900004, ppm=0.2),
        FormulaCandidate("C9H16N4O4", "[M+H]+", 3, prob=0.900003, ppm=-0.4),  # rank 2 was outside the tolerance
        FormulaCandidate("C11H16N3O3", "[M+Na]+", 1, prob=0.99, ppm=0.1),
    ], adducts=["[M+H]+", "[M+Na]+"], note="formulas found")
    rig.mp.setattr(pipeline.formula_step, "run_adduct_and_formula", lambda *a, **k: (adduct_result("[M+H]+", "[M+Na]+"), formulas))
    events = []
    result = pipeline.run(form(), progress=lambda *e: events.append(e))
    table = next(b for e in events if e[3] for b in e[3] if b["id"] == "formula_table")
    cols = table["columns"]
    rows = table["data"]["rows"]
    assert [r[cols.index("Rank")] for r in rows] == [c["rank"] for c in result["candidates"]] == [1, 2, 3]
    assert [r[cols.index("Formula")] for r in rows] == [c["formula"] for c in result["candidates"]]
    assert [r[cols.index("BUDDY rank (within adduct)")] for r in rows] == [1, 3, 1]
    assert [r[cols.index("BUDDY probability")] for r in rows] == ["0.900", "0.900", "0.990"]
    assert "adduct" in table["caption"] and "not a joint probability" in table["caption"]
    assert result["answer"]["ids"] == [c["id"] for c in result["candidates"]]
    assert rig.calls["frigid"][0][0][3] == "C10H20N2O5"  # no score-based pooling across adducts


def test_steps_1_and_2_failing_still_gives_a_mass_answer(rig):
    rig.mp.setattr(pipeline.formula_step, "run_adduct_and_formula", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("buddy gone")))
    result = pipeline.run(form())
    assert [s["status"] for s in result["steps"]][:2] == ["failed", "failed"] and "buddy gone" in result["steps"][0]["note"]
    assert result["answer"]["type"] == "mass_only" and result["answer"]["level"] == "5" and result["context"]["adduct"] == "[M+H]+"
    assert rig.calls["msfinder"][0][0][1] == "[M+H]+"                                       # the default adduct of the polarity


def test_user_adduct_is_passed_on(rig):
    pipeline.run(form(precursor={"precursor_mz": "217.0696", "ion_mode": "Positive", "adduct": "[M+Na]+"}))
    assert rig.calls["steps12"][0][1]["adduct"] == "[M+Na]+"


# ---------------------------------------------------------------- the peak caps
def peaks_form(n, **over):
    """n peaks, intensity = their number (so the most intense are the highest m/z), on a 195.0877 precursor."""
    text = "\n".join(f"{100 + i * 0.1:.4f} {i + 1}" for i in range(n))
    return form(ms2={"source": "Paste", "value": text}, **over)


def test_every_engine_step_gets_the_100_most_intense_peaks_and_the_library_all_of_them(rig):
    rig.library([], passed=False)
    rig.msfinder([structure_candidate("msfinder", "m1", CAFFEINE, 0.8)], good_enough=False)
    rig.frigid([structure_candidate("frigid", "f1", ETHANOL, 0.5)])
    result = pipeline.run(peaks_form(150, libraries=["MassBank"]))
    assert pipeline.MAX_ENGINE_PEAKS == 100
    seen = {name: rig.calls[name][0][0][idx] for name, idx in (("steps12", 2), ("msfinder", 2), ("frigid", 2), ("isomers", 1))}
    expected = [(round(100 + i * 0.1, 4), float(i + 1)) for i in range(50, 150)]
    for name, peaks in seen.items():
        assert len(peaks) == 100 and [(round(m, 4), i) for m, i in peaks] == expected, name
    assert len(rig.calls["library"][0][0][0]) == 150                           # the library search is not an engine: it sees every peak
    assert any("150 peaks" in w and "100 most intense" in w for w in result["warnings"])


def test_a_spectrum_within_the_cap_is_passed_on_unchanged_and_without_a_warning(rig):
    result = pipeline.run(peaks_form(100))
    assert len(rig.calls["steps12"][0][0][2]) == 100 and not any("most intense" in w for w in result["warnings"])


def test_a_spectrum_of_more_than_5000_peaks_is_a_plain_400(rig):
    assert pipeline.MAX_INPUT_PEAKS == 5000
    assert len(pipeline.parse_inputs(peaks_form(5000)).ms2) == 5000
    with pytest.raises(InputError, match="5001 peaks; at most 5000"):
        pipeline.parse_inputs(peaks_form(5001))
    r = client.post("/api/pipelines/esi_ms2_annotation/run", json=peaks_form(5001))
    assert r.status_code == 400 and "5001 peaks; at most 5000" in r.json()["detail"] and not app_module._jobs
    assert "steps12" not in rig.calls


# ---------------------------------------------------------------- the isomer check
def two_candidates(rig):
    rig.msfinder([structure_candidate("msfinder", "m1", ETHANOL, 0.8, formula="C2H6O"), structure_candidate("msfinder", "m2", ETHER, 0.5, formula="C2H6O"),
                  structure_candidate("msfinder", "m3", PROPANOL, 0.4, formula="C3H8O")], good_enough=True)
    rig.steps12(("[M+H]+",), (("C2H6O", "[M+H]+"),))


def test_isomer_single_verdict(rig):
    two_candidates(rig)
    rig.isomers(isomer_result("single", "c1"))
    result = pipeline.run(form())
    pool = rig.calls["isomers"][0][0][0]
    assert [c["id"] for c in pool] == ["c1", "c2"] and pool[0]["smiles"] == ETHANOL                  # the same-formula candidates only, ranked first
    assert rig.calls["isomers"][0][1]["engine"] == "msfinder"
    assert result["answer"]["type"] == "single" and result["answer"]["isomer_verdict"] == "single"
    assert [c["isomer"]["role"] for c in result["candidates"]] == ["top", "separated", None]
    assert result["steps"][6]["status"] == "done" and "tell the top-1 apart" in result["steps"][6]["note"]
    assert result["isomer"]["verdict"] == "single" and "log" not in result["isomer"]


def test_isomer_table_shows_each_smiles_with_its_structure_on_hover(rig):
    two_candidates(rig)
    rig.isomers(isomer_result("single", "c1"))
    events = []
    pipeline.run(form(), progress=lambda *e: events.append(e))
    table = next(b for e in events if e[3] for b in e[3] if b["id"] == "isomer_table")
    top, other = table["data"]["rows"]
    assert table["columns"][1] == "SMILES" and top[1]["text"] == ETHANOL and other[1]["text"] == ETHER
    assert "<svg" in other[1]["hover_svg"]                                                    # the 2D structure is drawn on hover
    assert top[5] == "—" and other[5]["text"] == "1"                                          # how many peaks only the top-1 explains
    assert other[5]["help"] == ["m/z 45.0334 · 30.0 % · C2H5O+ · CCO"]                        # m/z, intensity, ion formula, fragment SMILES


def test_isomer_table_survives_a_smiles_that_cannot_be_drawn(rig, monkeypatch):
    two_candidates(rig)
    rig.isomers(isomer_result("single", "c1"))
    monkeypatch.setattr(pipeline, "depiction_svg", lambda *a, **k: (_ for _ in ()).throw(ValueError("drawing broke")))
    events = []
    pipeline.run(form(), progress=lambda *e: events.append(e))
    table = next(b for e in events if e[3] for b in e[3] if b["id"] == "isomer_table")
    assert [row[1] for row in table["data"]["rows"]] == [ETHANOL, ETHER]                      # plain SMILES, no hover structure


def test_isomer_indistinguishable_gives_a_group_answer_and_keeps_the_level(rig):
    two_candidates(rig)
    rig.isomers(isomer_result("indistinguishable", "c1", group=["c1", "c2"]))
    result = pipeline.run(form())
    assert result["answer"]["type"] == "group" and result["answer"]["ids"] == ["c1", "c2"] and "cannot separate" in result["answer"]["statement"]
    assert result["answer"]["level"] == "3"                                                       # the verdict never changes the level
    assert [c["isomer"]["role"] for c in result["candidates"]][:2] == ["top", "not separated"]
    assert "cannot be separated" in result["steps"][6]["note"]


def test_isomer_challengers_are_shown_but_do_not_change_the_verdict(rig):
    two_candidates(rig)
    rig.isomers(isomer_result("single", "c1", challengers=["c2"]))
    events = []
    result = pipeline.run(form(), progress=lambda *e: events.append(e))
    notes = next(b for e in events if e[0] == 7 and e[3] for b in e[3] if b["id"] == "isomer_notes")
    assert notes["data"]["items"][0].startswith("Another isomer explains peaks the top-1 cannot: isomer c2")
    assert result["answer"]["type"] == "single" and result["answer"]["isomer_verdict"] == "single"
    assert "For information only" in result["answer"]["statement"] and "explains observed peaks" in result["answer"]["statement"]
    assert result["candidates"][1]["isomer"]["challenger"] is True and "Another isomer explains peaks the top-1 cannot" in result["steps"][6]["note"]


def test_isomer_insufficient_is_reported_and_the_answer_says_so(rig):
    two_candidates(rig)
    rig.isomers(isomer_result("insufficient", "c1", notes=("MS-FINDER could not run: busy",), error="busy"))
    result = pipeline.run(form())
    assert result["steps"][6]["status"] == "failed" and result["steps"][6]["note"] == "MS-FINDER could not run: busy"
    assert result["answer"]["isomer_verdict"] == "insufficient" and "could not be done" in result["answer"]["statement"]
    rig.isomers(isomer_result("insufficient", "c1", notes=("No isomer to compare with.",)))
    assert pipeline.run(form())["steps"][6]["status"] == "skipped"


def test_isomer_step_that_raises_does_not_stop_the_run(rig):
    two_candidates(rig)
    rig.isomers(None, error=RuntimeError("boom"))
    result = pipeline.run(form())
    assert result["steps"][6]["status"] == "failed" and result["answer"]["type"] == "single" and result["answer"]["isomer_verdict"] == "insufficient"


def test_database_isomer_checkbox_reaches_the_check(rig):
    two_candidates(rig)
    pipeline.run(form(local_isomers=True))
    assert rig.calls["isomers"][0][1]["local_isomers"] is True
    pipeline.run(form())
    assert rig.calls["isomers"][1][1]["local_isomers"] is False


def test_top_candidate_without_a_structure_skips_the_isomer_check(rig):
    rig.library([library_candidate("lib1", None, "Structureless", 0.9, False)])
    result = pipeline.run(form(libraries=["MassBank"]))
    assert result["steps"][6]["status"] == "skipped" and "no structure" in result["steps"][6]["note"] and "isomers" not in rig.calls
    assert result["answer"]["type"] == "single"


# ---------------------------------------------------------------- re-rank, progress, blocks
def test_rerank_failure_keeps_the_candidates(rig):
    rig.msfinder([structure_candidate("msfinder", "m1", CAFFEINE, 0.4), structure_candidate("msfinder", "m2", ETHER, 0.7, formula="C2H6O")], good_enough=True)
    rig.mp.setattr(pipeline, "rerank", lambda *a, **k: (_ for _ in ()).throw(ValueError("rules broke")))
    result = pipeline.run(form())
    assert result["steps"][5]["status"] == "failed" and "rules broke" in result["steps"][5]["note"]
    assert [c["rank"] for c in result["candidates"]] == [1, 2] and result["candidates"][0]["smiles"] == ETHER       # ordered by the layer score
    assert result["answer"]["type"] == "single" and any("Re-rank failed" in w for w in result["warnings"])


def test_real_rerank_runs_on_the_merged_candidates(monkeypatch, rig):
    monkeypatch.setattr(pipeline, "rerank", REAL_RERANK)
    rig.msfinder([structure_candidate("msfinder", "m1", CAFFEINE, 0.7, name="Caffeine")], good_enough=True)
    result = pipeline.run(form(omics="Metabolomics", additive="Formic acid", alkali="Absent", ph=3.0))
    top = result["candidates"][0]
    assert top["rank"] == 1 and top["final_score"] and top["detect"]["omics"] in ("metabolome", "both") and top["detect"]["factors"]


def test_progress_reports_every_step_in_order_with_blocks(rig):
    rig.library([library_candidate("lib1", CAFFEINE, "Caffeine", 0.93, True)], passed=True)
    events = []
    result = pipeline.run(form(libraries=["MassBank"]), progress=lambda *e: events.append(e))
    started = [e[0] for e in events if e[1] == "running" and e[2] is None]
    assert started == [1, 2, 3, 6, 7, 8]                                                       # 4 and 5 are skipped, not started
    ended = {e[0]: e for e in events if e[1] in ("done", "skipped", "failed")}
    assert set(ended) == set(range(1, 9)) and ended[4][1] == "skipped" and ended[4][3] == []
    assert any(e[0] == 6 and e[2] == "scored 1/1" for e in events)
    ids = {b["id"] for e in events if e[3] for b in e[3]}
    assert {"adduct_table", "formula_table", "library_table", "library_gate", "ranked_table", "answer_notes"} <= ids
    table = next(b for e in events if e[3] for b in e[3] if b["id"] == "library_table")
    assert table["type"] == "table" and len(table["columns"]) == len(table["data"]["rows"][0]) and table["data"]["rows"][0][0] == "Caffeine"
    assert [s["id"] for s in result["steps"]] == list(range(1, 9)) and all(s["seconds"] >= 0 for s in result["steps"])
    assert [s["title"] for s in result["steps"]] == [s["title"] for s in PIPELINES[0]["stages"]]


def test_a_failing_progress_listener_does_not_stop_the_run(rig):
    def bad(*e):
        raise RuntimeError("listener down")
    assert pipeline.run(form(), progress=bad)["answer"]["type"] == "formula_only"


def test_result_is_json_safe(rig):
    import json
    rig.msfinder([structure_candidate("msfinder", "m1", CAFFEINE, 0.7)], good_enough=True)
    json.dumps(pipeline.run(form()))


def test_result_carries_what_the_view_button_needs(rig):
    """"View" on the ranked table: the peaks the isomer step considers, and per candidate the adduct MetFrag can use (None = no button)."""
    caffeine = library_candidate("lib1", CAFFEINE, "Caffeine", 0.93, True)
    rig.library([caffeine], passed=True)
    rig.isomers(isomer_result("single", "c1"))
    result = pipeline.run(form(libraries=["MassBank"]))
    query = result["fragment_query"]
    assert query["precursor_mz"] == 195.0877 and query["tol_da"] == 0.02 and query["peaks"]
    assert all(0 < rel <= 100 and mz < 195.0877 for mz, rel in query["peaks"])           # the precursor is not a fragment peak
    assert [c["fragment_adduct"] for c in result["candidates"]] == ["[M+H]+"]


def test_candidates_without_a_structure_or_a_metfrag_adduct_get_no_view_button(rig):
    rig.steps12(("[M+H]+",), (("C8H10N4O2", "[M+H]+"),))
    result = pipeline.run(form())
    assert result["answer"]["type"] == "formula_only"
    assert [c["fragment_adduct"] for c in result["candidates"]] == [None]               # a formula row has no structure to fragment


# ---------------------------------------------------------------- the form
def test_parse_inputs_reads_the_form():
    q = pipeline.parse_inputs(form(ms1=CAFFEINE_MS1, libraries=["MassBank", "MassBank"], omics="Lipidomics", additive="Ammonium formate", alkali="Absent",
                                   ph="3.5", halogen=True, local_isomers=True, tol_ppm="5"))
    assert q.precursor_mz == 195.0877 and q.positive and q.adduct is None and len(q.ms2) == 27 and len(q.ms1) == 3
    assert q.libraries == ["MassBank"] and q.omics == "Lipidomics" and q.additive == "Ammonium formate" and q.alkali is False and q.ph == 3.5
    assert q.halogen and q.local_isomers and q.tol_ppm == 5 and q.tol_da == 0.02
    d = pipeline.parse_inputs({k: v for k, v in form().items() if k in ("precursor", "ms2")})              # the defaults
    assert (d.tol_ppm, d.tol_da, d.omics, d.additive, d.alkali, d.ph, d.halogen, d.local_isomers) == (10.0, 0.02, "Either / unknown", "Unknown", None, 2.7, False, False)


def test_parse_inputs_uses_the_header_precursor():
    ms2 = {"source": "Paste", "value": "BEGIN IONS\nPEPMASS=195.0877\n138.0662 999\n195.0876 58\nEND IONS"}
    assert pipeline.parse_inputs(form(precursor={"precursor_mz": "", "ion_mode": "Positive", "adduct": "Auto"}, ms2=ms2)).precursor_mz == 195.0877


@pytest.mark.parametrize("over, message", [
    ({"precursor": {"precursor_mz": "", "ion_mode": "Positive", "adduct": "Auto"}}, "Precursor m/z is required"),
    ({"precursor": {"precursor_mz": "195.0877", "ion_mode": "Positive", "adduct": "[M-H]-"}}, "not a positive-ion adduct"),
    ({"precursor": {"precursor_mz": "195.0877", "ion_mode": "Positive", "adduct": "nonsense"}}, "Unknown adduct"),
    ({"libraries": ["NIST"]}, "Libraries must be"), ({"libraries": "MassBank"}, "Libraries must be"),
    ({"omics": "Proteomics"}, "Omics must be"), ({"additive": "Water"}, "additive must be"), ({"alkali": "maybe"}, "Sodium / potassium"),
    ({"ph": 15}, "at most 14"), ({"tol_da": 0}, "must be positive"), ({"halogen": "yes"}, "true or false"),
    ({"ms2": {"source": "Paste", "value": "no peaks"}}, "no valid peaks")])
def test_parse_inputs_rejects_bad_forms(over, message):
    with pytest.raises(InputError, match=message):
        pipeline.parse_inputs(form(**over))


def test_spec_form_matches_the_pipeline():
    from mzdesk.platform_fit import adduct_rules
    from mzdesk import spec
    pipe = PIPELINES[0]
    fields = {f["id"]: f for f in pipe["inputs"]}
    assert "structure" not in fields and set(fields) == {"precursor", "ms2", "ms1", "libraries", "tol_ppm", "tol_da", "omics", "halogen", "local_isomers",
                                                         "additive", "alkali", "ph"}
    assert fields["tol_da"]["value"] == 0.02 and fields["omics"]["options"][0] == "Either / unknown" and fields["omics"]["options"] == list(pipeline.OMICS_CHOICES)
    assert fields["additive"]["options"] == ["Unknown", *adduct_rules.ADDITIVES] and spec.PIPELINE_ADDITIVES == list(adduct_rules.ADDITIVES)
    assert fields["halogen"]["label"] == "Allow Cl/Br in formulas" and fields["halogen"]["value"] is False and fields["local_isomers"]["value"] is False
    assert fields["local_isomers"]["label"] == "Also compare with database isomers (slower, stricter)"
    assert {f["section"] for f in pipe["inputs"] if f["id"] in ("additive", "alkali", "ph")} == {"Advanced"}
    assert list(pipeline.STEP_TITLES) == list(range(1, 9)) and [s["id"] for s in pipe["stages"]] == list(range(1, 9))
    assert set(pipeline.ALKALI_CHOICES) == {o.lower() for o in fields["alkali"]["options"]}


def test_both_examples_parse_and_cover_the_library_and_the_no_library_path():
    examples = {e["label"]: e["values"] for e in PIPELINES[0]["examples"]}
    assert len(examples) == 2
    library, plain = (pipeline.parse_inputs(v) for v in examples.values())
    assert library.libraries == ["MassBank", "GNPS"] and len(library.ms1) == 3 and len(library.ms2) == 19
    assert plain.libraries == [] and len(plain.ms2) == 27
    assert set(next(iter(examples.values()))) <= {f["id"] for f in PIPELINES[0]["inputs"]}


# ---------------------------------------------------------------- the job API
@pytest.fixture(autouse=True)
def clean_jobs():
    app_module._jobs.clear()
    yield
    app_module._jobs.clear()


def fake_result():
    return {"answer": {"type": "none"}, "candidates": [], "steps": [{"id": n, "title": t, "status": "done", "seconds": 1.5, "note": "ok"}
                                                                     for n, t in pipeline.STEP_TITLES.items()], "warnings": []}


def wait_done(job_id, timeout=10):
    end = time.time() + timeout
    while time.time() < end:
        body = client.get(f"/api/pipelines/jobs/{job_id}").json()
        if body["status"] not in ("running", "queued"):
            return body
        time.sleep(0.02)
    raise AssertionError("job did not finish")


def test_run_unknown_pipeline_is_404_and_bad_form_is_400():
    assert client.post("/api/pipelines/nope/run", json={}).status_code == 404
    r = client.post("/api/pipelines/esi_ms2_annotation/run", json=form(omics="Proteomics"))
    assert r.status_code == 400 and "Omics" in r.json()["detail"]
    assert client.post("/api/pipelines/esi_ms2_annotation/run", json={}).status_code == 400
    assert client.get("/api/pipelines/jobs/unknown").status_code == 404 and not app_module._jobs


def test_every_pipeline_has_a_runner():
    assert set(app_module._PIPELINE_RUNNERS) == {p["id"] for p in PIPELINES}


def test_job_lifecycle_live_steps_then_result(monkeypatch):
    go, midway = threading.Event(), threading.Event()

    def fake(query, progress):
        progress(1, "running", None, None)
        progress(1, "done", "adducts done", [{"type": "table", "id": "adduct_table", "title": "t", "columns": ["a"], "data": {"rows": [["x"]]}}])
        progress(2, "running", None, None)
        progress(2, "running", "scored 3/9", None)
        midway.set()
        go.wait(5)
        return fake_result()
    monkeypatch.setattr(pipeline, "run_query", fake)
    r = client.post("/api/pipelines/esi_ms2_annotation/run", json=form())
    assert r.status_code == 200 and r.json()["status"] == "running" and r.json()["job_id"]
    job_id = r.json()["job_id"]
    assert midway.wait(5)
    body = client.get(f"/api/pipelines/jobs/{job_id}").json()
    assert body["status"] == "running" and body["result"] is None and body["error"] is None
    steps = {s["id"]: s for s in body["steps"]}
    assert [s["title"] for s in body["steps"]] == list(pipeline.STEP_TITLES.values()) and set(steps[3]) == {"id", "title", "status", "seconds", "note", "blocks"}
    assert steps[1]["status"] == "done" and steps[1]["note"] == "adducts done" and steps[1]["blocks"][0]["id"] == "adduct_table"
    assert steps[2]["status"] == "running" and steps[2]["note"] == "scored 3/9" and steps[3]["status"] == "queued"
    go.set()
    body = wait_done(job_id)
    assert body["status"] == "done" and body["error"] is None and body["result"]["answer"]["type"] == "none"
    assert all(s["status"] == "done" and s["seconds"] == 1.5 for s in body["steps"])
    assert body["steps"][0]["blocks"][0]["id"] == "adduct_table"                                  # blocks survive the end of the run


def test_job_error_is_reported(monkeypatch):
    def boom(query, progress):
        progress(1, "running", None, None)
        raise RuntimeError("engine exploded")
    monkeypatch.setattr(pipeline, "run_query", boom)
    body = wait_done(client.post("/api/pipelines/esi_ms2_annotation/run", json=form()).json()["job_id"])
    assert body["status"] == "error" and "RuntimeError: engine exploded" in body["error"] and body["result"] is None
    assert body["steps"][0]["status"] == "failed" and body["steps"][1]["status"] == "queued"


def test_job_runs_the_real_pipeline_with_mocked_layers(rig):
    rig.msfinder([structure_candidate("msfinder", "m1", CAFFEINE, 0.8)], good_enough=True)
    body = wait_done(client.post("/api/pipelines/esi_ms2_annotation/run", json=form()).json()["job_id"])
    assert body["status"] == "done" and body["result"]["answer"]["layer"] == "msfinder"
    assert [s["status"] for s in body["steps"]] == ["done", "done", "skipped", "done", "skipped", "done", "skipped", "done"]
    assert any(b["id"] == "ranked_table" for s in body["steps"] for b in s["blocks"])


def test_jobs_are_counted_in_the_usage_stats(monkeypatch):
    from mzdesk.api import stats
    seen = []
    monkeypatch.setattr(stats, "record", lambda request, kind, **kw: seen.append((kind, kw["module"], kw["status"])))
    monkeypatch.setattr(pipeline, "run_query", lambda q, progress: fake_result())
    wait_done(client.post("/api/pipelines/esi_ms2_annotation/run", json=form()).json()["job_id"])
    assert seen == [("run", "esi_ms2_annotation", 200)]


def test_at_most_20_jobs_are_kept_and_old_ones_expire(monkeypatch):
    monkeypatch.setattr(pipeline, "run_query", lambda q, progress: fake_result())
    ids = []
    for _ in range(app_module.JOB_LIMIT + 3):
        ids.append(client.post("/api/pipelines/esi_ms2_annotation/run", json=form()).json()["job_id"])
        wait_done(ids[-1])
    assert len(app_module._jobs) == app_module.JOB_LIMIT
    assert client.get(f"/api/pipelines/jobs/{ids[0]}").status_code == 404 and client.get(f"/api/pipelines/jobs/{ids[-1]}").status_code == 200
    app_module._jobs[ids[-1]]["touched"] -= app_module.JOB_TTL_S + 1                               # an hour after the last update
    assert client.get(f"/api/pipelines/jobs/{ids[-1]}").status_code == 404


def test_running_jobs_are_never_dropped_to_make_room(monkeypatch):
    go = threading.Event()
    monkeypatch.setattr(pipeline, "run_query", lambda q, progress: (go.wait(10), fake_result())[1])
    try:
        ids = [client.post("/api/pipelines/esi_ms2_annotation/run", json=form()).json()["job_id"] for _ in range(app_module.JOB_LIMIT)]
        r = client.post("/api/pipelines/esi_ms2_annotation/run", json=form())
        assert r.status_code == 429 and len(app_module._jobs) == app_module.JOB_LIMIT
    finally:
        go.set()
    for job_id in ids:
        wait_done(job_id)


def test_jobs_run_one_at_a_time_and_the_second_waits_in_the_queue(monkeypatch):
    gates, log = {}, []

    def slow(query, progress):
        name = query.precursor_mz                                  # each form below has its own precursor m/z
        gates.setdefault(name, threading.Event())
        log.append(("start", name))
        progress(1, "running", None, None)
        assert gates[name].wait(10)
        log.append(("end", name))
        return fake_result()
    monkeypatch.setattr(pipeline, "run_query", slow)
    first = client.post("/api/pipelines/esi_ms2_annotation/run", json=form_at(195.0877)).json()
    deadline = time.time() + 5
    while ("start", 195.0877) not in log and time.time() < deadline:
        time.sleep(0.01)
    second = client.post("/api/pipelines/esi_ms2_annotation/run", json=form_at(300.0)).json()
    third = client.post("/api/pipelines/esi_ms2_annotation/run", json=form_at(400.0)).json()
    assert first["status"] == "running" and second["status"] == "queued" and third["status"] == "queued"
    time.sleep(0.2)
    assert log == [("start", 195.0877)]                            # the others did not start while the first runs
    one, two = (client.get(f"/api/pipelines/jobs/{j['job_id']}").json() for j in (second, third))
    assert one["status"] == two["status"] == "queued" and one["result"] is None
    assert one["note"] == "waiting for the previous run (position 1 in the queue)" and two["note"].endswith("(position 2 in the queue)")
    assert one["steps"][0]["note"] == one["note"] and all(s["status"] == "queued" for s in one["steps"])
    gates.setdefault(300.0, threading.Event()).set()
    gates.setdefault(400.0, threading.Event()).set()
    gates[195.0877].set()
    for j in (first, second, third):
        assert wait_done(j["job_id"])["status"] == "done"
    assert log == [("start", 195.0877), ("end", 195.0877), ("start", 300.0), ("end", 300.0), ("start", 400.0), ("end", 400.0)]
    assert client.get(f"/api/pipelines/jobs/{second['job_id']}").json()["note"] == ""


def test_a_failed_job_hands_the_turn_to_the_next_one(monkeypatch):
    def run(query, progress):
        if query.precursor_mz == 300.0:
            raise RuntimeError("engine exploded")
        return fake_result()
    monkeypatch.setattr(pipeline, "run_query", run)
    ids = [client.post("/api/pipelines/esi_ms2_annotation/run", json=form_at(mz)).json()["job_id"] for mz in (300.0, 195.0877)]
    assert [wait_done(j)["status"] for j in ids] == ["error", "done"]


def test_job_store_is_safe_under_concurrent_polling(monkeypatch):
    def fake(query, progress):
        for n in range(1, 9):
            progress(n, "running", None, None)
            progress(n, "done", f"step {n}", [{"type": "notes", "id": f"n{n}", "title": "t", "data": {"items": ["x"] * 50}}])
        return fake_result()
    monkeypatch.setattr(pipeline, "run_query", fake)
    errors = []

    def poll(job_id):
        for _ in range(60):
            r = client.get(f"/api/pipelines/jobs/{job_id}")
            if r.status_code != 200:
                errors.append(r.status_code)
    job_id = client.post("/api/pipelines/esi_ms2_annotation/run", json=form()).json()["job_id"]
    threads = [threading.Thread(target=poll, args=(job_id,)) for _ in range(4)]
    [t.start() for t in threads]
    [t.join() for t in threads]
    assert not errors and wait_done(job_id)["status"] == "done"


# ---------------------------------------------------------------- real end-to-end (opt-in: BUDDY, MassBank / GNPS, MS-FINDER, FRIGID)
RUN_REAL = os.environ.get("MZDESK_RUN_PIPELINE") == "1"
needs_real = pytest.mark.skipif(not RUN_REAL, reason="real runs take minutes and up to 5 GB: set MZDESK_RUN_PIPELINE=1")


def _real_libraries():
    libs = [name for name in ("MassBank", "GNPS") if lx.is_installed(name)]
    if not libs:
        pytest.skip("no local library installed (python -m mzdesk.ms2.library_setup)")
    return libs


def _show(title, result):
    print(f"\n=== {title} ===")
    print(result["answer"]["statement"])
    print({k: v for k, v in result["answer"].items() if k != "statement"})
    for s in result["steps"]:
        print(f"  step {s['id']} {s['title']:<22} {s['status']:<8} {s['seconds']:>7.1f} s  {s['note']}")
    for c in result["candidates"][:5]:
        print(f"  #{c['rank']} {c['name'] or c['smiles']} {c['formula']} {c['adduct']} {c['layer_label']} level {c['level']} final {c['final_score']}")
    print("  warnings:", result["warnings"])


@needs_real
def test_real_caffeine_through_the_library_path():
    libs = _real_libraries()
    result = pipeline.run(form(ms2={"source": "Paste", "value": CAFFEINE_MASSBANK_MS2}, ms1=CAFFEINE_MS1, libraries=libs))
    _show("caffeine, library path", result)
    s = {x["id"]: x for x in result["steps"]}
    assert result["context"]["gate_passed"] and s[4]["status"] == s[5]["status"] == "skipped"
    assert result["answer"]["level"] == "2a" and result["answer"]["layer"] == "library" and result["candidates"][0]["name"].lower().startswith("caffeine")
    assert s[6]["status"] == "done" and s[7]["status"] in ("done", "skipped")


@needs_real
def test_real_caffeine_without_a_library_through_msfinder_and_frigid():
    result = pipeline.run(form(ms2={"source": "Paste", "value": CAFFEINE_MSFINDER_MS2}, ms1=CAFFEINE_MS1, libraries=[]))
    _show("caffeine, no library", result)
    s = {x["id"]: x for x in result["steps"]}
    assert s[3]["status"] == "skipped" and s[4]["status"] == "done" and result["candidates"]
    assert s[5]["status"] in ("done", "skipped") and result["answer"]["type"] in ("single", "group") and result["answer"]["level"] in ("3", "H")


@needs_real
def test_real_frigid_runs_when_msfinder_is_not_good_enough(monkeypatch):
    monkeypatch.setattr(pipeline.msfinder_step, "GOOD_ENOUGH", 0.99)                       # caffeine's MS-FINDER score (0.82) now counts as weak
    result = pipeline.run(form(ms2={"source": "Paste", "value": CAFFEINE_MSFINDER_MS2}, ms1=CAFFEINE_MS1, libraries=[]))
    _show("caffeine, no library, MS-FINDER not good enough", result)
    s = {x["id"]: x for x in result["steps"]}
    layers = {src["layer"] for c in result["candidates"] for src in c["sources"]}
    assert s[4]["status"] == "done" and s[5]["status"] == "done" and {"msfinder", "frigid"} <= layers
