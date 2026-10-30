"""Step 4 reruns the formulas a TIME-limited MS-FINDER run left unfinished with TreeDepth=1 (owner's design 2026-10-30); predict_queries is mocked."""
import pytest

from mzdesk.annotation import msfinder_step as step
from mzdesk.ms2 import msfinder_engine as mf
from test_annotation_steps45 import CAFFEINE_KEY, CAFFEINE_PEAKS, answer, struct

TIME = "MS-FINDER did not finish within 90 s and was stopped."
MEMORY = "MS-FINDER used more than 6 GB of memory and was stopped."
F1, F2 = "C8H10N4O2", "C8H18OS2"
UNFINISHED = {"formulas": [], "structures": [], "searches": {}}


@pytest.fixture
def engine(monkeypatch):
    """install(*runs): the n-th predict_queries call does what runs[n] says: a dict {"answers": [...], "stopped": str | None},
    or an MsFinderError to raise. calls collects what the step sent."""
    calls = []

    def install(*runs):
        def fake(queries, **kw):
            calls.append({"queries": queries, **kw})
            run = runs[len(calls) - 1]
            if isinstance(run, Exception):
                raise run
            return {"results": {q["key"]: a for q, a in zip(queries, run["answers"])}, "seconds": 5.0, "peak_rss_mb": 2300.0, "log": "",
                    "stopped": run.get("stopped")}
        monkeypatch.setattr(mf, "predict_queries", fake)
        return calls
    return install


def run(formulas=(F1, F2), **kw):
    return step.run_msfinder_layer(195.0877, "[M+H]+", CAFFEINE_PEAKS, list(formulas), tol_ppm=5.0, tol_da=0.01, **kw)


def test_a_time_stop_reruns_only_the_unfinished_formulas_with_treedepth_1_and_marks_them(engine):
    calls = engine({"answers": [answer(F1, struct("Caffeine", CAFFEINE_KEY, 8.0)), UNFINISHED], "stopped": TIME},
                   {"answers": [answer(F2, struct("Thio", "AAAAAAAAAAAAAA-BBBBBBBBBB-C", 6.0, formula=F2))]})
    layer, summary = run()
    assert len(calls) == 2
    assert "TreeDepth" not in calls[0]["values"] and calls[0]["timeout"] == step.ENGINE_TIMEOUT_S
    assert calls[1]["values"]["TreeDepth"] == 1 and calls[1]["timeout"] == step.RETRY_TIMEOUT_S and calls[1]["values"]["Ms2Tolerance"] == 0.01
    assert [q["key"] for q in calls[1]["queries"]] == ["f1"]                                     # only the unfinished formula
    assert [c.name for c in layer.candidates] == ["Caffeine", "Thio"]                            # merged, best first
    assert [(f["formula"], f["status"], f["tree_depth"], f["note"]) for f in summary["formulas"]] == [
        (F1, "searched", None, None), (F2, "searched", 1, step.RETRY_NOTE)]
    assert [c.evidence["msfinder"].get("tree_depth") for c in layer.candidates] == [None, 1]
    assert summary["retried"] == [F2] and summary["engine_seconds"] == 10.0
    (warning,) = summary["warnings"]
    assert f"Finished at the shipped TreeDepth: {F1}." in warning and f"Finished at TreeDepth=1: {F2}" in warning and "Still not finished" not in warning
    assert f"TreeDepth=1 for {F2} after the TreeDepth=2 search timed out" in layer.note and "stopped early" not in layer.note


def test_a_first_run_that_raises_the_time_stop_is_repeated_for_all_formulas(engine):
    calls = engine(mf.MsFinderError(TIME), {"answers": [answer(F1, struct("Caffeine", CAFFEINE_KEY, 7.5)), answer(F2)]})
    layer, summary = run()
    assert [q["key"] for q in calls[1]["queries"]] == ["f0", "f1"] and calls[1]["values"]["TreeDepth"] == 1
    assert layer.status == "done" and [c.name for c in layer.candidates] == ["Caffeine"]
    assert [(f["status"], f["tree_depth"]) for f in summary["formulas"]] == [("searched", 1), ("no structures", 1)]
    assert summary["retried"] == [F1, F2]


def test_a_rerun_that_is_stopped_too_says_what_is_still_not_finished(engine):
    engine({"answers": [answer(F1, struct("Caffeine", CAFFEINE_KEY, 8.0)), UNFINISHED], "stopped": TIME},
           {"answers": [UNFINISHED], "stopped": "MS-FINDER did not finish within 120 s and was stopped."})
    layer, summary = run()
    assert layer.status == "done" and [c.name for c in layer.candidates] == ["Caffeine"]
    assert [(f["status"], f["note"]) for f in summary["formulas"]][1] == ("not finished", "not finished at TreeDepth=1 either")
    (warning,) = summary["warnings"]
    assert f"Still not finished: {F2} (MS-FINDER did not finish within 120 s and was stopped)" in warning
    assert f"still not finished at TreeDepth=1: {F2}" in layer.note


def test_when_the_rerun_is_stopped_with_nothing_too_the_step_fails_with_both_reasons(engine):
    calls = engine(mf.MsFinderError(TIME), mf.MsFinderError("MS-FINDER did not finish within 120 s and was stopped."))
    layer, _ = run()
    assert len(calls) == 2 and layer.status == "failed" and "within 90 s" in layer.note and "TreeDepth=1 failed too" in layer.note and "within 120 s" in layer.note


@pytest.mark.parametrize("error", [mf.MsFinderError("MS-FINDER is busy with another request. Try again in a few seconds."),
                                   mf.MsFinderError("MS-FINDER is not installed: x missing."), mf.MsFinderError(MEMORY)])
def test_an_error_that_is_not_a_time_stop_is_not_repeated(engine, error):
    calls = engine(error)
    layer, _ = run()
    assert len(calls) == 1 and layer.status == "failed" and layer.note == str(error)


def test_a_memory_stop_with_partial_results_is_not_repeated(engine):
    calls = engine({"answers": [answer(F1, struct("Caffeine", CAFFEINE_KEY, 8.0)), UNFINISHED], "stopped": MEMORY})
    layer, summary = run()
    assert len(calls) == 1 and summary["retried"] is None and "Not finished: C8H18OS2" in summary["warnings"][0] and "stopped early" in layer.note


def test_a_normal_run_is_not_repeated(engine):
    calls = engine({"answers": [answer(F1, struct("Caffeine", CAFFEINE_KEY, 8.0)), answer(F2)]})
    layer, summary = run()
    assert len(calls) == 1 and summary["retried"] is None and summary["warnings"] == [] and all(f["tree_depth"] is None for f in summary["formulas"])


def test_a_time_stop_where_every_formula_wrote_its_result_is_not_repeated(engine):
    calls = engine({"answers": [answer(F1, struct("Caffeine", CAFFEINE_KEY, 8.0)), answer(F2)], "stopped": TIME})
    _, summary = run()
    assert len(calls) == 1 and summary["retried"] is None


def test_without_a_formula_list_a_time_stop_with_nothing_reruns_the_whole_run(engine):
    own = {"formulas": [], "structures": [struct("Caffeine", CAFFEINE_KEY, 8.0)], "searches": {F1: 20.0}}
    calls = engine(mf.MsFinderError(TIME), {"answers": [own]})
    layer, summary = run(())
    assert len(calls) == 2 and calls[1]["values"]["TreeDepth"] == 1 and [q["key"] for q in calls[1]["queries"]] == ["own"]
    assert layer.status == "done" and [c.name for c in layer.candidates] == ["Caffeine"] and summary["retried"] == []
    assert summary["formulas"][0]["tree_depth"] == 1 and summary["formulas"][0]["note"] == step.RETRY_NOTE
    assert "TreeDepth=1 for" in layer.note


def test_without_a_formula_list_a_stop_that_kept_structures_is_not_repeated(engine):
    own = {"formulas": [], "structures": [struct("Caffeine", CAFFEINE_KEY, 8.0)], "searches": {F1: 20.0}}
    calls = engine({"answers": [own], "stopped": TIME})
    run(())
    assert len(calls) == 1


def test_an_empty_search_is_finished_in_a_normal_run_and_unfinished_in_a_stopped_one(engine):
    empty = {"formulas": [], "structures": [], "searches": {F2: 0.0}}            # the engine keeps a 0-byte .sfd in "searches" only when not stopped
    calls = engine({"answers": [answer(F1, struct("Caffeine", CAFFEINE_KEY, 8.0)), {**empty, "formulas": answer(F2)["formulas"]}]})
    _, summary = run()
    assert len(calls) == 1 and [f["status"] for f in summary["formulas"]] == ["searched", "no structures"]
    calls.clear()                                                                # (the fixture's list is shared by both installs)
    engine({"answers": [answer(F1, struct("Caffeine", CAFFEINE_KEY, 8.0)), UNFINISHED], "stopped": TIME}, {"answers": [{**empty, "formulas": answer(F2)["formulas"]}]})
    _, summary = run()                                                           # a search that finished empty is rerun too: cheap (seconds)
    assert len(calls) == 2 and summary["retried"] == [F2] and summary["formulas"][1]["status"] == "no structures"
