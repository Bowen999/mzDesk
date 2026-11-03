"""ESI MS/MS Annotation steps 4-5: the MS-FINDER layer (mzdesk/annotation/msfinder_step.py), the FRIGID layer (frigid_step.py) and the engine
function they use for several MS-FINDER queries in one run (msfinder_engine.predict_queries / parse_sfd_headers).

Unit tests replay recorded engine output (tests/fixtures/annotation/, real runs of 2026-10-28) or hand-made engine answers; nothing starts.
  msfinder/<case>/   raw MS-FINDER result files of a real step run: formula_result.txt, <query key>/<FORMULA>.sfd (the first structures only),
                     input.json (what the step was given). caffeine: MassBank MSBNK-LCSB-LU119905 (the page example, CC BY 4.0), BUDDY's
                     formula list; kaempferol: MassBank MSBNK-LCSB-LU089855 (CC BY, tests/fixtures/isomers); caffeine_own: no formula list.
  frigid/            raw FRIGID driver results (tests/fixtures/frigid/caffeine_quick.result.json is the caffeine one).
Real runs are opt-in, one at a time: MZDESK_RUN_MSFINDER_STEP=1 (about 5-12 s, 2.3 GB) and MZDESK_RUN_FRIGID=1 (20-50 s, 4-5 GB).
"""
import json
import os
from pathlib import Path

import pytest

from mzdesk.annotation import frigid_step, msfinder_step
from mzdesk.annotation.merge import merge_candidates
from mzdesk.annotation.schema import LayerResult
from mzdesk.ms2 import frigid_engine as fe
from mzdesk.ms2 import msfinder_engine as mf
from mzdesk.spec import CAFFEINE_MSFINDER_MS2

FIX = Path(__file__).parent / "fixtures" / "annotation"
CAFFEINE_PEAKS = [tuple(map(float, ln.split())) for ln in CAFFEINE_MSFINDER_MS2.splitlines()]
CAFFEINE_IK14 = "RYYVLZVUVIJVGH"


def _unavailable():
    try:
        mf.check_msfinder()
    except mf.MsFinderError as exc:
        return str(exc)
    return None


NO_MSFINDER = _unavailable()
RUN_MSFINDER = os.environ.get("MZDESK_RUN_MSFINDER_STEP") == "1"
RUN_FRIGID = os.environ.get("MZDESK_RUN_FRIGID") == "1"


# ---------------------------------------------------------------- score map
def test_score_map_is_linear_absolute_and_clipped():
    n = msfinder_step.score_norm
    assert n(5.0) == 0.0 and n(9.0) == 1.0 and n(7.0) == pytest.approx(0.5) and n(7.4) == pytest.approx(msfinder_step.GOOD_ENOUGH)
    assert n(2.0) == 0.0 and n(12.0) == 1.0 and n(None) == 0.0 and n(float("nan")) == 0.0
    scores = [4.0 + 0.1 * k for k in range(70)]
    assert [n(s) for s in scores] == sorted(n(s) for s in scores)                       # monotone
    assert frigid_step.score_norm(0.0) == 0.0 and frigid_step.score_norm(None) == 0.0 and frigid_step.score_norm(1.0) == 1.0
    fs = [k / 50 for k in range(51)]
    assert [frigid_step.score_norm(s) for s in fs] == sorted(frigid_step.score_norm(s) for s in fs)


# ---------------------------------------------------------------- MS-FINDER layer, engine replaced by hand-made answers
CAFFEINE_KEY = "RYYVLZVUVIJVGH-UHFFFAOYSA-N"
CAFFEINE_SMILES = "O=C1C2=C(N=CN2C)N(C(=O)N1C)C"
CAFFEINE_PPM = (195.0877 - 195.08765203118054) / 195.08765203118054 * 1e6          # the page example against C8H10N4O2 [M+H]+
SUBSCORES = {name: 0.5 for _, name in mf.SFD_SCORES}


def struct(name, key, score, formula="C8H10N4O2", smiles="C", dbs=(("HMDB", "HMDB0001847;HMDB0002"), ("PubChem", "2519")), ontology="Xanthines"):
    """One structure as predict_queries returns it."""
    return {"formula": formula, "name": name, "inchikey": key, "smiles": smiles, "score": score, "databases": list(dbs), "ontology": ontology,
            "subscores": dict(SUBSCORES)}


def answer(formula, *structs, accepted=True, seconds=1.5, fscore=4.4):
    """The engine's answer for one query: MS-FINDER's formula list (the formula itself, or nothing), its structures, the search time."""
    return {"formulas": [{"formula": formula, "neutral_mass": 0.0, "mass_error_da": 0.0, "score": fscore}] if accepted else [],
            "structures": list(structs), "searches": {formula: seconds}}


@pytest.fixture
def engine(monkeypatch):
    """Replace mf.predict_queries: install(*answers) answers the queries in order. The returned list collects what the step sent."""
    calls = []

    def install(*answers, stopped=None, error=None, seconds=6.0):
        def fake(queries, **kw):
            calls.append({"queries": queries, **kw})
            if error:
                raise error
            return {"results": {q["key"]: a for q, a in zip(queries, answers)}, "seconds": seconds, "peak_rss_mb": 2300.0, "log": "", "stopped": stopped}
        monkeypatch.setattr(mf, "predict_queries", fake)
        return calls
    return install


def run(formulas=("C8H10N4O2",), *, adduct="[M+H]+", precursor=195.0877, peaks=CAFFEINE_PEAKS, tol_ppm=5.0, tol_da=0.01, **kw):
    return msfinder_step.run_msfinder_layer(precursor, adduct, peaks, list(formulas), tol_ppm=tol_ppm, tol_da=tol_da, **kw)


def test_candidates_fill_the_schema_and_each_formula_is_searched_at_its_own_mz(engine):
    calls = engine(answer("C8H10N4O2", struct("Caffeine", CAFFEINE_KEY, 8.1594, smiles=CAFFEINE_SMILES),
                          struct("1,3,9-Trimethylxanthin (Isocoffein);x", "LPHGQDQBBGAPDZ-UHFFFAOYSA-N", 7.3261, smiles="O=C1C=2N=CN(C2N(C(=O)N1C)C)C")),
                   answer("C8H18OS2", struct("junk", "AAAAAAAAAAAAAA-BBBBBBBBBB-C", 6.0, formula="C8H18OS2")))
    layer, summary = run(("C8H10N4O2", "C8H18OS2"))
    (call,) = calls
    assert [q["key"] for q in call["queries"]] == ["f0", "f1"] and all(q["adduct"] == "[M+H]+" and q["positive"] for q in call["queries"])
    assert [q["precursor_mz"] for q in call["queries"]] == [195.08765, 195.08718]                 # each formula's own theoretical m/z, 5 decimals
    assert call["queries"][0]["ms2"] == sorted(CAFFEINE_PEAKS) and "ms1" not in call["queries"][0]
    assert call["values"] == {"Ms2Tolerance": 0.01, "StructurePredictionTimeOut": 0.5, "Ms1Tolerance": 0.0002, "TryTopNMolecularFormulaSearch": 2}
    assert call["halogen"] is False and call["timeout"] == 90 and call["rss_cap_mb"] == 6144

    assert isinstance(layer, LayerResult) and layer.layer == "msfinder" and layer.status == "done" and layer.seconds > 0
    assert [c.name for c in layer.candidates] == ["Caffeine", "1,3,9-Trimethylxanthin (Isocoffein)", "junk"]          # best first, synonyms cut at ";"
    top = layer.candidates[0]
    assert (top.id, top.layer, top.smiles, top.inchikey, top.ik14, top.formula, top.adduct) == (
        "msfinder_1", "msfinder", CAFFEINE_SMILES, CAFFEINE_KEY, CAFFEINE_IK14, "C8H10N4O2", "[M+H]+")
    assert top.mz_error_ppm == pytest.approx(CAFFEINE_PPM, abs=1e-6) and top.hypothesis is False
    assert top.layer_score == 8.1594 and top.layer_score_norm == pytest.approx((8.1594 - 5.0) / 4.0)
    assert top.evidence == {"msfinder": {"total_score": 8.1594, "subscores": SUBSCORES, "databases": ["HMDB", "PubChem"],
                                         "database_ids": {"HMDB": "HMDB0001847", "PubChem": "2519"}, "ontology": "Xanthines", "formula_score": 4.4,
                                         "rank_in_formula": 1, "formula_structures": 2}}
    assert layer.candidates[2].evidence["msfinder"]["rank_in_formula"] == 1 and layer.candidates[2].formula == "C8H18OS2"
    assert summary["good_enough"] is True and summary["best_norm"] == top.layer_score_norm and summary["best_score"] == 8.1594
    assert summary["margin"] == pytest.approx(8.1594 - 7.3261) and summary["n_candidates"] == 3 and summary["mode"] == "formula_list"
    assert [(f["formula"], f["n_structures"], f["best_score"], f["status"], f["search_seconds"], f["accepted"]) for f in summary["formulas"]] == [
        ("C8H10N4O2", 2, 8.1594, "searched", 1.5, True), ("C8H18OS2", 1, 6.0, "searched", 1.5, True)]
    assert summary["engine_seconds"] == 6.0 and summary["peak_rss_mb"] == 2300.0 and summary["seconds"] == layer.seconds and summary["warnings"] == []
    assert "3 structures from 2 formulas" in layer.note and "best Caffeine (C8H10N4O2)" in layer.note and layer.note.endswith("good enough.")
    json.dumps(summary, allow_nan=False)
    json.dumps(layer.to_dict(), allow_nan=False)


def test_the_layer_merges_like_any_other(engine):
    engine(answer("C8H10N4O2", struct("Caffeine", CAFFEINE_KEY, 8.1594, smiles=CAFFEINE_SMILES)))
    layer, _ = run()
    (cand,) = merge_candidates([layer])
    assert (cand.level, cand.ik14, cand.layer) == ("3", CAFFEINE_IK14, "msfinder") and set(cand.evidence) == {"msfinder"} and cand.hypothesis is False
    assert [(s["layer"], s["rank_in_layer"]) for s in cand.sources] == [("msfinder", 1)]


def test_one_candidate_per_compound_best_first_and_top_n(engine):
    engine(answer("C8H10N4O2", struct("L", "AAAAAAAAAAAAAA-BBBBBBBBBB-C", 6.0), struct("flat", "AAAAAAAAAAAAAA-CCCCCCCCCC-C", 7.0),     # same ik14
                  struct("other", "DDDDDDDDDDDDDD-BBBBBBBBBB-C", 6.5), struct("third", "EEEEEEEEEEEEEE-BBBBBBBBBB-C", 5.5)))
    layer, summary = run(top_n=2)
    assert [(c.id, c.name, c.layer_score) for c in layer.candidates] == [("msfinder_1", "flat", 7.0), ("msfinder_2", "other", 6.5)]
    assert summary["n_candidates"] == 2 and summary["good_enough"] is False and "not good enough (below 0.6)" in layer.note


def test_a_formula_outside_the_precursor_tolerance_is_not_searched(engine):
    calls = engine(answer("C8H10N4O2", struct("Caffeine", CAFFEINE_KEY, 8.0)))
    layer, summary = run(("C8H10N4O2", "C9H12N4O2"))
    assert [q["key"] for q in calls[0]["queries"]] == ["f0"] and layer.status == "done"
    (skipped,) = summary["skipped"]
    assert skipped["formula"] == "C9H12N4O2" and "ppm from the precursor (limit 5 ppm)" in skipped["reason"]
    assert any("C9H12N4O2" in w and "not searched" in w for w in summary["warnings"]) and "1 formula not searched" in layer.note
    calls.clear()
    layer, summary = run(("C9H12N4O2",))
    assert calls == [] and layer.status == "skipped" and "within the precursor tolerance" in layer.note and summary["n_candidates"] == 0


def test_at_most_six_formulas_are_searched_and_the_rest_is_said(engine):
    formulas = [f"C{n}H{n + 4}N2O2" for n in range(6, 14)]                                       # eight formulas, any tolerance
    calls = engine(*[answer(f) for f in formulas[:6]])
    layer, summary = run(formulas, tol_ppm=1e9)
    assert len(calls[0]["queries"]) == 6 and [s["formula"] for s in summary["skipped"]] == formulas[6:]
    assert all("only the first 6 formulas are searched" in s["reason"] for s in summary["skipped"]) and "2 formulas not searched" in layer.note
    (warning,) = summary["warnings"]                                                            # one warning for the whole group, not one per formula
    assert warning == f"Only the first 6 formulas are searched; 2 were not searched: {formulas[6]}, {formulas[7]}."


def test_adducts_go_to_msfinder_in_its_spelling_and_negative_mode_works(engine):
    calls = engine(answer("C6H12O6", struct("glucose", "WQZGKKKJIJFFOK-UHFFFAOYSA-N", 6.0, formula="C6H12O6")))
    layer, _ = run(("C6H12O6",), adduct="[M+HCOO]-", precursor=225.0616)
    (q,) = calls[0]["queries"]
    assert q["adduct"] == "[M+FA-H]-" and q["positive"] is False and q["precursor_mz"] == 225.06159
    assert layer.candidates[0].adduct == "[M+HCOO]-"                                           # canonical on the candidate
    calls.clear()
    layer, _ = run(("C6H12O6",), adduct="[M+FA-H]-", precursor=225.0616)                       # any spelling in
    assert calls[0]["queries"][0]["adduct"] == "[M+FA-H]-" and layer.candidates[0].adduct == "[M+HCOO]-"


def test_formulas_may_be_dicts_or_objects_with_their_own_adduct(engine):
    class Item:
        formula, adduct = "C8H7N3O2", "[M+NH4]+"          # the second adduct of step 1: same ion, another formula

    calls = engine(answer("C8H10N4O2", struct("a", CAFFEINE_KEY, 8.0)), answer("C8H7N3O2", struct("b", "AAAAAAAAAAAAAA-BBBBBBBBBB-C", 7.0, formula="C8H7N3O2")))
    layer, _ = run([{"formula": "C8H10N4O2", "adduct": "[M+H]+", "prob": 0.9}, Item(), "C8H10N4O2", {"nothing": 1}, None, ""])
    qs = calls[0]["queries"]
    assert [(q["key"], q["adduct"]) for q in qs] == [("f0", "[M+H]+"), ("f1", "[M+NH4]+")]               # duplicates and unreadable items are dropped
    assert [(c.formula, c.adduct) for c in layer.candidates] == [("C8H10N4O2", "[M+H]+"), ("C8H7N3O2", "[M+NH4]+")]


def test_no_formula_list_runs_msfinders_own_formula_search_within_tol_ppm(engine):
    ms1 = [(195.0877, 100.0), (196.0911, 9.0)]
    calls = engine({"formulas": [{"formula": "C8H10N4O2", "neutral_mass": 194.08, "mass_error_da": 0.0, "score": 4.4},
                                 {"formula": "C8H14N2O4", "neutral_mass": 202.0954, "mass_error_da": 0.0005, "score": 3.1}],
                    "structures": [struct("Caffeine", CAFFEINE_KEY, 8.1594), struct("wrong formula", "AAAAAAAAAAAAAA-BBBBBBBBBB-C", 6.0, formula="C8H14N2O4")],
                    "searches": {"C8H10N4O2": 2.0, "C8H14N2O4": 0.0}})
    layer, summary = run((), ms1_peaks=ms1)
    (q,) = calls[0]["queries"]
    assert q["key"] == "own" and q["precursor_mz"] == 195.0877 and q["ms1"] == ms1 and calls[0]["values"]["TryTopNMolecularFormulaSearch"] == 5
    assert calls[0]["values"]["Ms1Tolerance"] == 0.001                                       # 5 ppm of 195 is 0.98 mDa: the shipped 1 mDa is the floor
    assert summary["mode"] == "own_formulas"
    assert [c.formula for c in layer.candidates] == ["C8H10N4O2"]                              # C8H14N2O4 is far from the precursor
    assert summary["formulas"][0]["ppm"] == pytest.approx(CAFFEINE_PPM, abs=1e-6) and len(summary["formulas"]) == 1
    calls.clear()
    engine({"formulas": [], "structures": [], "searches": {}})
    layer, _ = run((), tol_ppm=20.0)
    assert calls[0]["values"]["Ms1Tolerance"] == pytest.approx(20e-6 * 195.0877) and layer.status == "done" and "its own formula search" in layer.note


def test_the_isotope_pattern_is_not_sent_when_a_formula_list_exists(engine):
    calls = engine(answer("C8H10N4O2", struct("Caffeine", CAFFEINE_KEY, 8.0)))
    run(ms1_peaks=[(195.0877, 100.0), (196.0911, 9.0)])
    assert all("ms1" not in q for q in calls[0]["queries"])


@pytest.mark.parametrize("res,status,warning", [
    (answer("C8H10N4O2", accepted=False), "not accepted", "element rules do not accept C8H10N4O2"),
    (answer("C8H10N4O2"), "no structures", None),
    (answer("C8H10N4O2", struct("Caffeine", CAFFEINE_KEY, 8.0), seconds=29.0), "time limit", "reached the 30 s limit"),
    (answer("C8H10N4O2", struct("Caffeine", CAFFEINE_KEY, 8.0), seconds=26.9), "searched", None),
])
def test_what_became_of_each_formula_is_reported(engine, res, status, warning):
    engine(res)
    layer, summary = run()
    (f,) = summary["formulas"]
    assert f["status"] == status and layer.status == "done"
    assert (warning is None and not summary["warnings"]) or (warning in " ".join(summary["warnings"]))
    if status != "searched":
        assert layer.note.rstrip(".").endswith(")") or "MS-FINDER" in layer.note
    if status == "time limit":
        assert "time limit reached for C8H10N4O2" in layer.note


def test_no_structure_at_all_is_a_done_layer_with_no_candidates_and_a_plain_note(engine):
    engine(answer("C8H10N4O2"), answer("C8H18OS2"))
    layer, summary = run(("C8H10N4O2", "C8H18OS2"))
    assert layer.status == "done" and layer.candidates == [] and summary["good_enough"] is False and summary["best_norm"] is None
    assert layer.note == "MS-FINDER's local structure databases hold no structure for C8H10N4O2, C8H18OS2."
    engine(answer("C8H10N4O2", accepted=False), answer("C8H18OS2"))
    layer, _ = run(("C8H10N4O2", "C8H18OS2"))
    assert layer.note == "MS-FINDER gave no structure for C8H10N4O2, C8H18OS2: C8H10N4O2 not accepted; C8H18OS2 no structures."


def test_a_run_that_was_stopped_keeps_what_had_finished(engine):
    unfinished = {"formulas": [], "structures": [], "searches": {}}                              # no result file: the search did not end
    engine(answer("C8H10N4O2", struct("Caffeine", CAFFEINE_KEY, 8.0), seconds=3.0), unfinished,
           stopped="MS-FINDER used more than 6 GB of memory and was stopped.")
    layer, summary = run(("C8H10N4O2", "C8H18OS2"))
    assert layer.status == "done" and [c.name for c in layer.candidates] == ["Caffeine"] and "MS-FINDER was stopped early" in layer.note
    assert [(f["formula"], f["status"], f["accepted"]) for f in summary["formulas"]] == [("C8H10N4O2", "searched", None), ("C8H18OS2", "not finished", None)]
    (warning,) = summary["warnings"]                                                            # one warning names both groups
    assert "stopped before it finished" in warning and "Finished: C8H10N4O2." in warning and "Not finished: C8H18OS2 (" in warning


def test_a_stopped_run_tells_a_finished_search_without_structures_from_an_unfinished_one(engine):
    engine(answer("C8H10N4O2"), {"formulas": [], "structures": [], "searches": {}}, stopped="MS-FINDER used more than 6 GB of memory and was stopped.")
    _, summary = run(("C8H10N4O2", "C8H18OS2"))
    assert [f["status"] for f in summary["formulas"]] == ["no structures", "not finished"]
    (warning,) = summary["warnings"]
    assert "Finished: C8H10N4O2." in warning and "Not finished: C8H18OS2" in warning
    engine(answer("C8H10N4O2"), stopped="MS-FINDER used more than 6 GB of memory and was stopped.")
    _, summary = run(("C8H10N4O2",))
    assert "Finished: C8H10N4O2." in summary["warnings"][0] and "Not finished" not in summary["warnings"][0]


@pytest.mark.parametrize("adduct,peaks,why", [("[M+Li-H]-", CAFFEINE_PEAKS, "cannot use the adduct [M+Li-H]-"), ("Auto", CAFFEINE_PEAKS, "cannot use the adduct Auto"),
                                              (None, CAFFEINE_PEAKS, "(none given)"), ("[M+H]+", [], "no MS2 peaks"),
                                              ("[M+H]+", [(0.0, 5.0), (100.0, 0.0), (float("nan"), 1.0)], "no MS2 peaks")])
def test_skipped_when_msfinder_does_not_apply_and_the_engine_is_not_started(engine, adduct, peaks, why):
    calls = engine()
    layer, summary = run(adduct=adduct, peaks=peaks)
    assert layer.status == "skipped" and why in layer.note and layer.candidates == [] and calls == [] and summary["good_enough"] is False


@pytest.mark.parametrize("error,fragment", [(mf.MsFinderError("MS-FINDER is busy with another request (one run at a time).", log="x"), "is busy"),
                                            (mf.MsFinderError("MS-FINDER is not installed: Mono is missing."), "not installed"),
                                            (ValueError("MS-FINDER takes at most 2000 peaks per spectrum"), "ValueError: MS-FINDER takes at most 2000"),
                                            (RuntimeError("boom"), "RuntimeError: boom")])
def test_a_missing_or_failing_engine_gives_a_failed_layer_never_an_exception(engine, error, fragment):
    engine(error=error)
    layer, summary = run()
    assert layer.status == "failed" and fragment in layer.note and layer.candidates == [] and layer.seconds >= 0
    assert summary["good_enough"] is False and summary["best_norm"] is None and summary["n_candidates"] == 0


def test_peaks_are_cleaned_and_the_hundred_most_intense_are_sent(engine):
    calls = engine(answer("C8H10N4O2"))
    peaks = [(50.0 + k, float(k + 1)) for k in range(150)] + [(10.0, 0.0), (float("inf"), 5.0), (-3.0, 4.0)]
    run(peaks=peaks)
    sent = calls[0]["queries"][0]["ms2"]
    assert len(sent) == 100 and sent == sorted(sent) and min(i for _, i in sent) == 51.0 and (10.0, 0.0) not in sent


def test_the_halogen_switch_follows_the_formulas_and_a_formula_msbuddy_cannot_compute_is_skipped(engine):
    calls = engine(answer("C8H9ClN4O2"))
    _, summary = run(("C8H9ClN4O2", "C8H10N4O2Si"), tol_ppm=1e9)               # silicon is not in BUDDY's alphabet: no m/z, no search
    assert calls[0]["halogen"] is True and "silicon" not in calls[0] and len(calls[0]["queries"]) == 1
    assert summary["skipped"] == [{"formula": "C8H10N4O2Si", "reason": "C8H10N4O2Si cannot be searched as [M+H]+"}]
    calls.clear()
    engine(answer("C8H10N4O2"))
    run(("C8H10N4O2",))
    assert calls[0]["halogen"] is False


# ---------------------------------------------------------------- recorded MS-FINDER runs (the real engine code on real result files)
def _case(name):
    return json.loads((FIX / "msfinder" / name / "input.json").read_text())


@pytest.fixture
def replay(monkeypatch, tmp_path):
    """A fake install and a stand-in for the console process that writes the recorded result files of a case into the run's folders.
    replay(case) returns a dict that records the parameter file and the query files of the run."""
    exe, ini = tmp_path / "MsfinderConsoleApp.exe", tmp_path / "x.ini"
    exe.write_text("")
    ini.write_text(mf.INI_PATH.read_text())
    monkeypatch.setattr(mf, "status", lambda: {"mono": "/usr/bin/mono", "dir": tmp_path, "exe": exe, "resources": tmp_path, "ini": ini, "thread_fix": True})
    seen = {}

    def install(case, code=1, error=None):
        def fake_run(cmd, cwd, env, timeout, rss_cap_mb, **kw):
            in_dir, out_dir = Path(cmd[cmd.index("-i") + 1]), Path(cmd[cmd.index("-o") + 1])
            seen.update(ini=Path(cmd[cmd.index("-m") + 1]).read_text(), queries={p.stem: p.read_text() for p in sorted(in_dir.glob("*.mat"))},
                        timeout=timeout, rss_cap_mb=rss_cap_mb)
            src = FIX / "msfinder" / case
            for key in seen["queries"]:
                for sfd in sorted((src / key).glob("*.sfd")) if (src / key).is_dir() else []:
                    (in_dir / key).mkdir(exist_ok=True)
                    (in_dir / key / sfd.name).write_text(sfd.read_text())
            if error is not None:
                raise error
            (out_dir / "Formula result-2099.txt").write_text((src / "formula_result.txt").read_text())
            return code, "log", 2300.0
        monkeypatch.setattr(mf, "run_watched", fake_run)
        return seen
    return install


def _ini_values(text):
    return dict(ln.split("=", 1) for ln in text.splitlines() if "=" in ln and not ln.startswith("#"))


def _query(key, precursor=195.08765, positive=True, adduct="[M+H]+", peaks=CAFFEINE_PEAKS):
    return {"key": key, "precursor_mz": precursor, "positive": positive, "adduct": adduct, "ms2": peaks}


def test_structure_files_give_scores_databases_and_ontology_without_the_fragment_tables():
    rows = mf.parse_sfd_headers((FIX / "msfinder" / "caffeine" / "f0" / "C8H10N4O2.sfd").read_text())
    assert [r["name"] for r in rows] == ["Caffeine", "cyclo(gly-his)", "1,3,9-Trimethylxanthin (Isocoffein)", "UNPD32360"]
    top = rows[0]
    assert (top["inchikey"], top["smiles"], top["ontology"], top["score"]) == (CAFFEINE_KEY, CAFFEINE_SMILES, "Xanthines", 8.1604)
    assert top["subscores"] == {"hr_rules": 0.9385, "bond_cleavage": 0.1184, "mass_accuracy": 0.9934, "fragment_linkage": 1.0, "bond_energy": 0.7227,
                                "database": 0.9167, "substructure": 0.2}
    assert top["databases"][:2] == [("HMDB", "HMDB0001847"), ("KNApSAcK", "C00001492")] and ("PubChem", "2519") in top["databases"]
    assert ("FooDB", "FDB002100;FDB017732") in top["databases"]                                     # an ID list stays one entry
    assert rows[1]["databases"] == [("COCONUT", "CNP0431962")] and rows[1]["ontology"] == "Alpha amino acids and derivatives"
    assert mf.parse_sfd_headers("") == [] and mf.parse_sfd_headers("NAME: a\nTotalScore: x\n")[0]["score"] is None


def test_one_run_takes_every_query_and_returns_each_ones_formulas_structures_and_search_times(replay):
    seen = replay("caffeine")
    res = mf.predict_queries([_query("f0"), _query("f1", 195.08718)], values={"Ms1Tolerance": 0.0002, "TryTopNMolecularFormulaSearch": 2, "Ms2Tolerance": 0.02},
                             timeout=33, rss_cap_mb=5000)
    assert sorted(seen["queries"]) == ["f0", "f1"] and seen["timeout"] == 33 and seen["rss_cap_mb"] == 5000
    assert "NAME: f1" in seen["queries"]["f1"] and "PRECURSORMZ: 195.08718" in seen["queries"]["f1"] and "PRECURSORTYPE: [M+H]+" in seen["queries"]["f0"]
    values = _ini_values(seen["ini"])
    assert (values["Ms1Tolerance"], values["TryTopNMolecularFormulaSearch"], values["Ms2Tolerance"]) == ("0.0002", "2", "0.02")
    assert values["PubChemNeverUse"] == "True" and values["Fcheck"] == "False" and mf.INI_PATH.read_text() == Path(mf.status()["ini"]).read_text()
    f0, f1 = res["results"]["f0"], res["results"]["f1"]
    assert [f["formula"] for f in f0["formulas"]] == ["C8H10N4O2"] and [f["formula"] for f in f1["formulas"]] == ["C8H18OS2"]       # each query's own rows
    assert f0["formulas"][0]["score"] == 4.408 and f1["formulas"][0]["score"] == 3.227
    assert [s["formula"] for s in f0["structures"]] == ["C8H10N4O2"] * 4 and f0["structures"][0]["score"] == 8.1604
    assert [s["score"] for s in f0["structures"]] == sorted((s["score"] for s in f0["structures"]), reverse=True)
    assert set(f0["searches"]) == {"C8H10N4O2"} and set(f1["searches"]) == {"C8H18OS2"}
    assert all(v is None or v >= 0 for v in (*f0["searches"].values(), *f1["searches"].values()))
    assert res["stopped"] is None and res["peak_rss_mb"] == 2300.0 and res["seconds"] >= 0


def test_an_empty_structure_file_is_a_search_without_structures(replay):
    replay("kaempferol")
    res = mf.predict_queries([_query("f1", 285.04115, False, "[M-H]-"), _query("f4", 285.03981, False, "[M-H]-")])
    assert res["results"]["f1"]["structures"] == [] and res["results"]["f4"]["structures"] == []
    assert set(res["results"]["f1"]["searches"]) == {"C8H10N6O4S", "C7H16N2O6P2"} and [f["formula"] for f in res["results"]["f4"]["formulas"]] == ["C7H14N2O8S"]


@pytest.mark.parametrize("queries,match", [([], "No query"), ([_query("a b")], "unique and made of"), ([_query("a"), _query("a")], "unique"),
                                           ([{**_query("a"), "precursor_mz": 0}], "Precursor"), ([{**_query("a"), "adduct": "[M-H]-"}], "positive-mode"),
                                           ([{**_query("a"), "ms2": [(1.0, 1.0)] * 2001}], "at most 2000")])
def test_queries_msfinder_cannot_use_are_refused_before_it_starts(replay, queries, match):
    seen = replay("caffeine")
    with pytest.raises(ValueError, match=match):
        mf.predict_queries(queries)
    assert seen == {}


def test_a_run_that_is_stopped_keeps_the_searches_that_had_finished(replay):
    stop = mf.MsFinderError("MS-FINDER used more than 6 GB of memory and was stopped.", "tail")
    replay("caffeine", error=stop)
    res = mf.predict_queries([_query("f0"), _query("f1", 195.08718)])
    assert res["stopped"] == str(stop) and res["peak_rss_mb"] is None and res["log"] == "tail"
    assert [s["score"] for s in res["results"]["f0"]["structures"]][:1] == [8.1604] and res["results"]["f0"]["formulas"] == []
    replay("kaempferol", error=stop)
    with pytest.raises(mf.MsFinderError, match="was stopped"):                      # nothing had finished: the error is the answer
        mf.predict_queries([_query("f1", 285.04115, False, "[M-H]-")])


@pytest.mark.parametrize("stop,searched", [(True, {"C8H10N4O2"}), (False, {"C8H10N4O2", "C8H18OS2"})])
def test_a_0_byte_structure_file_is_unfinished_in_a_stopped_run_and_empty_in_a_finished_one(replay, monkeypatch, stop, searched):
    replay("caffeine")                                                          # installs the fake install; run_watched is replaced below

    def fake_run(cmd, cwd, env, timeout, rss_cap_mb, **kw):
        in_dir, out_dir = Path(cmd[cmd.index("-i") + 1]), Path(cmd[cmd.index("-o") + 1])
        (in_dir / "f0").mkdir()
        (in_dir / "f0" / "C8H10N4O2.sfd").write_text((FIX / "msfinder" / "caffeine" / "f0" / "C8H10N4O2.sfd").read_text())
        (in_dir / "f0" / "C8H18OS2.sfd").write_text("")                         # created at the start of a search, written at its end
        (out_dir / "Formula result-2099.txt").write_text((FIX / "msfinder" / "caffeine" / "formula_result.txt").read_text())
        if stop:
            raise mf.MsFinderError("MS-FINDER did not finish within 90 s and was stopped.", "tail")
        return 1, "log", 2300.0
    monkeypatch.setattr(mf, "run_watched", fake_run)
    res = mf.predict_queries([_query("f0")])
    assert set(res["results"]["f0"]["searches"]) == searched and bool(res["stopped"]) is stop
    assert res["results"]["f0"]["structures"]


def test_a_run_without_any_result_file_is_an_error(replay, monkeypatch):
    replay("caffeine")
    monkeypatch.setattr(mf, "run_watched", lambda *a, **k: (1, "log", 10.0))
    with pytest.raises(mf.MsFinderError, match="without writing a result file"):
        mf.predict_queries([_query("f0")])


def test_a_second_run_is_turned_away_while_another_heavy_run_holds_the_lock(replay, monkeypatch):
    replay("caffeine")
    monkeypatch.setattr(mf, "BUSY_WAIT_S", 0.05)
    with mf._RUN_LOCK:
        with pytest.raises(mf.MsFinderError, match="busy"):
            mf.predict_queries([_query("f0")])
    assert mf.predict_queries([_query("f0")])["results"]["f0"]["structures"]
    assert mf._RUN_LOCK.acquire(blocking=False)
    mf._RUN_LOCK.release()


# ---------------------------------------------------------------- the step on recorded runs
def _step(case, peaks, **over):
    inp = _case(case)
    args = {"tol_ppm": inp["tol_ppm"], "tol_da": inp["tol_da"], **over}
    return msfinder_step.run_msfinder_layer(inp["precursor_mz"], inp["adduct"], peaks, inp["formulas"], **args)


def test_recorded_caffeine_run_ranks_caffeine_first_for_the_two_formulas_of_the_list(replay):
    seen = replay("caffeine")
    layer, summary = _step("caffeine", CAFFEINE_PEAKS)
    q = seen["queries"]
    assert sorted(q) == ["f0", "f1"] and "PRECURSORMZ: 195.08765" in q["f0"] and "PRECURSORMZ: 195.08718" in q["f1"] and "MSTYPE: MS1" not in q["f0"]
    values = _ini_values(seen["ini"])
    assert values["StructurePredictionTimeOut"] == "0.5" and values["Ms1Tolerance"] == "0.0002" and values["TryTopNMolecularFormulaSearch"] == "2"
    assert layer.status == "done" and layer.note == "6 structures from 2 formulas; best Caffeine (C8H10N4O2), score 8.16, normalised 0.79: good enough."
    assert [(c.name, c.formula, c.layer_score) for c in layer.candidates] == [
        ("Caffeine", "C8H10N4O2", 8.1604), ("cyclo(gly-his)", "C8H10N4O2", 7.4538), ("1,3,9-Trimethylxanthin (Isocoffein)", "C8H10N4O2", 7.3271),
        ("UNPD32360", "C8H10N4O2", 7.2707), (layer.candidates[4].name, "C8H18OS2", 5.0408), (layer.candidates[5].name, "C8H18OS2", 4.8482)]
    top = layer.candidates[0]
    assert (top.ik14, top.smiles, top.adduct, top.hypothesis) == (CAFFEINE_IK14, CAFFEINE_SMILES, "[M+H]+", False)
    assert top.evidence["msfinder"]["subscores"]["mass_accuracy"] == 0.9934 and top.evidence["msfinder"]["ontology"] == "Xanthines"
    assert top.evidence["msfinder"]["database_ids"]["HMDB"] == "HMDB0001847" and top.evidence["msfinder"]["formula_score"] > 4
    assert layer.candidates[4].layer_score_norm == pytest.approx(0.0102, abs=1e-4) and layer.candidates[5].layer_score_norm == 0.0       # below 5.0
    assert summary["good_enough"] is True and summary["best_norm"] == pytest.approx(0.7901) and summary["margin"] == pytest.approx(0.7066)
    assert [(f["formula"], f["accepted"], f["n_structures"], f["status"]) for f in summary["formulas"]] == [("C8H10N4O2", True, 4, "searched"),
                                                                                                           ("C8H18OS2", True, 2, "searched")]
    json.dumps(layer.to_dict(), allow_nan=False)
    json.dumps(summary, allow_nan=False)


def test_recorded_kaempferol_run_in_negative_mode_with_five_formulas(replay):
    peaks = [tuple(p) for p in json.loads((Path(__file__).parent / "fixtures" / "isomers" / "kaempferol_lcsb75.json").read_text())["peaks"]]
    seen = replay("kaempferol")
    layer, summary = _step("kaempferol", peaks)
    assert sorted(seen["queries"]) == ["f0", "f1", "f2", "f3", "f4"]
    assert all("IONMODE: Negative" in t and "PRECURSORTYPE: [M-H]-" in t for t in seen["queries"].values())
    assert len(peaks) == 107 and "MSTYPE: MS2\nNum Peaks: 100\n" in seen["queries"]["f0"]                      # the 100 most intense peaks are sent
    assert [(c.name, c.layer_score, round(c.layer_score_norm, 4)) for c in layer.candidates][:2] == [("Kaempferol", 7.7146, 0.6786), ("Luteolin", 7.3795, 0.5949)]
    assert summary["good_enough"] is True and summary["margin"] == pytest.approx(0.3351) and layer.candidates[0].adduct == "[M-H]-"
    assert [f["status"] for f in summary["formulas"]] == ["searched"] + ["no structures"] * 4 and summary["warnings"] == []
    assert layer.note.startswith("4 structures from 5 formulas; best Kaempferol (C15H10O6), score 7.71, normalised 0.68: good enough")


def test_recorded_run_without_a_formula_list_uses_msfinders_own_formulas(replay):
    seen = replay("caffeine_own")
    ms1 = [(195.0877, 100.0), (196.0911, 9.0)]
    layer, summary = msfinder_step.run_msfinder_layer(195.0877, "[M+H]+", CAFFEINE_PEAKS, [], tol_ppm=5.0, tol_da=0.01, ms1_peaks=ms1)
    assert list(seen["queries"]) == ["own"] and "PRECURSORMZ: 195.08770" in seen["queries"]["own"] and seen["queries"]["own"].count("MSTYPE: MS1") == 1
    values = _ini_values(seen["ini"])
    assert values["TryTopNMolecularFormulaSearch"] == "5" and values["Ms1Tolerance"] == "0.001"
    assert summary["mode"] == "own_formulas" and [f["formula"] for f in summary["formulas"]] == ["C8H10N4O2", "C8H18OS2"]
    assert layer.candidates[0].name == "Caffeine" and layer.candidates[0].layer_score == 8.1594 and summary["good_enough"] is True


# ---------------------------------------------------------------- FRIGID layer, engine replaced by recorded results
CAFFEINE_FRIGID = json.loads((Path(__file__).parent / "fixtures" / "frigid" / "caffeine_quick.result.json").read_text())      # the driver's raw result
LU = json.loads((FIX / "frigid_lu087906.json").read_text())                                                                      # predict()'s result


def caffeine_result(**over):
    """What fe.predict returns for the recorded caffeine driver result (filtered like the engine does)."""
    kept, dropped = fe.filter_candidates(CAFFEINE_FRIGID["candidates"], "[M+H]+")
    res = {"formula": "C8H10N4O2", "adduct": "[M+H]+", "driver_adduct": "[M+H]+", "mode": "quick", "seed": 42, "candidates": kept, "dropped": dropped,
           "n_raw": 7, "n_samples": 16, "n_valid": 16, "n_formula_matched": 9, "n_unique_matched": 4, "n_peaks_assigned": 27, "device": "mps",
           "device_requested": "auto", "device_reason": "auto: mps", "seconds": 21.4, "driver_seconds": 17.3, "peak_footprint_mb": 5047.4,
           "peak_rss_mb": 1586.6, "warnings": [], "length_model": True, "log": ""}
    res.update(over)
    return res


@pytest.fixture
def frigid(monkeypatch):
    """Replace fe.predict: install(result) answers with it; install(error=...) raises. The returned list collects the arguments."""
    calls = []

    def install(result=None, error=None):
        def fake(formula, precursor_mz, adduct, peaks, mode="quick", **kw):
            calls.append({"formula": formula, "precursor_mz": precursor_mz, "adduct": adduct, "peaks": list(peaks), "mode": mode, **kw})
            if error:
                raise error
            return result if result is not None else caffeine_result()
        monkeypatch.setattr(fe, "predict", fake)
        return calls
    return install


def frigid_run(formula="C8H10N4O2", *, adduct="[M+H]+", precursor=195.0877, peaks=CAFFEINE_PEAKS, **kw):
    return frigid_step.run_frigid_layer(precursor, adduct, peaks, formula, **kw)


def test_frigid_structures_are_hypotheses_with_the_schema_fields(frigid):
    calls = frigid()
    layer = frigid_run()
    assert calls == [{"formula": "C8H10N4O2", "precursor_mz": 195.0877, "adduct": "[M+H]+", "peaks": sorted(CAFFEINE_PEAKS), "mode": "quick"}]
    assert isinstance(layer, LayerResult) and layer.layer == "frigid" and layer.status == "done" and layer.seconds >= 0
    assert [c.id for c in layer.candidates] == ["frigid_1", "frigid_2", "frigid_3", "frigid_4"]
    top = layer.candidates[0]
    assert (top.layer, top.smiles, top.ik14, top.inchikey, top.name, top.formula, top.adduct) == (
        "frigid", "Cn1c(=O)c2c(ncn2C)n(C)c1=O", CAFFEINE_IK14, None, "", "C8H10N4O2", "[M+H]+")
    assert top.hypothesis is True and top.layer_score == 1.0 and top.layer_score_norm == 1.0 and top.mz_error_ppm == pytest.approx(CAFFEINE_PPM, abs=1e-3)
    assert top.evidence == {"frigid": {"score": 1.0, "rank": 1, "mode": "quick", "n_samples": 16, "n_formula_matched": 9, "n_unique_matched": 4, "device": "mps"}}
    assert [round(c.layer_score, 3) for c in layer.candidates] == [1.0, 0.515, 0.5, 0.333]
    assert [round(c.layer_score_norm, 3) for c in layer.candidates] == [1.0, 0.03, 0.0, 0.0]                         # linear, 0 at 0.5
    assert layer.note == "4 structures for C8H10N4O2 (Quick mode, 16 samples); best fingerprint similarity 1.00, normalised 1.00."
    json.dumps(layer.to_dict(), allow_nan=False)


def test_frigid_layer_merges_as_a_hypothesis(frigid):
    frigid()
    cands = merge_candidates([frigid_run()])
    assert cands[0].level == "H" and cands[0].hypothesis is True and set(cands[0].evidence) == {"frigid"} and cands[0].ik14 == CAFFEINE_IK14
    assert cands[0].inchikey == "RYYVLZVUVIJVGH-UHFFFAOYSA-N"                                  # merge makes the full key from the SMILES


def test_recorded_run_of_a_second_compound_with_two_structures(frigid):
    frigid(LU["result"])
    peaks = [tuple(p) for p in LU["input"]["peaks"]]
    layer = frigid_run(LU["input"]["formula"], precursor=LU["input"]["precursor_mz"], peaks=peaks)
    assert [(c.ik14, c.layer_score, round(c.layer_score_norm, 3)) for c in layer.candidates] == [("BQGRVFPPZJPWPB", 1.0, 1.0), ("KVBANSJOXDHNIW", 0.842105, 0.684)]
    assert layer.candidates[0].mz_error_ppm == pytest.approx(0.031, abs=0.001) and layer.note.startswith("2 structures for C15H17NO3S (Quick mode, 16 samples)")


def test_top_n_mode_and_the_engine_spelling_of_the_adduct(frigid):
    calls = frigid(caffeine_result(adduct="[M+NH4]+", driver_adduct="[M+H3N+H]+", mode="full"))
    layer = frigid_run(adduct="[M+NH4]+", precursor=212.1142, top_n=2, mode="full")
    assert len(layer.candidates) == 2 and calls[0]["adduct"] == "[M+NH4]+" and calls[0]["mode"] == "full" and layer.candidates[0].adduct == "[M+NH4]+"
    assert layer.note.startswith("2 of 4 structures for C8H10N4O2 (Full mode, 16 samples)")
    assert frigid_run(top_n=0).candidates[0].id == "frigid_1" and len(frigid_run(top_n=0).candidates) == 1                  # at least one
    assert frigid_step.TOP_N == 5 and len(frigid_run().candidates) == 4


def test_frigid_with_no_structure_of_the_formula_is_done_with_no_candidates(frigid):
    frigid(caffeine_result(candidates=[], n_unique_matched=0, n_formula_matched=0))
    layer = frigid_run()
    assert layer.status == "done" and layer.candidates == [] and layer.note == "FRIGID generated no structure with the formula C8H10N4O2 (Quick mode, 16 samples)."


@pytest.mark.parametrize("adduct,formula,peaks,why", [
    ("[M-H]-", "C8H10N4O2", CAFFEINE_PEAKS, "FRIGID supports positive ions only; [M-H]- is a negative-ion adduct."),
    ("[M+FA-H]-", "C8H10N4O2", CAFFEINE_PEAKS, "positive ions only; [M+HCOO]- is a negative-ion adduct"),
    ("[2M+H]+", "C8H10N4O2", CAFFEINE_PEAKS, "does not support the adduct [2M+H]+; it takes [M+H]+, [M+Na]+, [M+K]+, [M+NH4]+"),
    ("Auto", "C8H10N4O2", CAFFEINE_PEAKS, "does not support the adduct Auto"),
    (None, "C8H10N4O2", CAFFEINE_PEAKS, "(none given)"),
    ("[M+H]+", None, CAFFEINE_PEAKS, "no molecular formula to start from"),
    ("[M+H]+", "  ", CAFFEINE_PEAKS, "no molecular formula to start from"),
    ("[M+H]+", "C8H10Zn", CAFFEINE_PEAKS, "FRIGID cannot use this formula: FRIGID cannot use Zn"),
    ("[M+H]+", "H2O", CAFFEINE_PEAKS, "has no carbon"),
    ("[M+H]+", "C8H10N4O2", [], "no MS2 peaks"),
    ("[M+H]+", "C8H10N4O2", [(0.0, 1.0), (100.0, 0.0)], "no MS2 peaks"),
])
def test_skipped_with_a_plain_reason_when_frigid_does_not_apply(frigid, adduct, formula, peaks, why):
    calls = frigid()
    layer = frigid_run(formula, adduct=adduct, peaks=peaks)
    assert layer.status == "skipped" and why.lower() in layer.note.lower() and layer.candidates == [] and calls == []


@pytest.mark.parametrize("error,status,fragment", [
    (fe.FrigidError("FRIGID is busy: another De novo structure run is running.", "log"), "failed", "is busy"),
    (fe.FrigidError("De novo structure runs FRIGID, which is not installed here."), "failed", "not installed"),
    (RuntimeError("boom"), "failed", "RuntimeError: boom"),
    (ValueError("FRIGID rejected the input: bad peaks"), "skipped", "rejected the input"),
])
def test_a_missing_or_failing_frigid_gives_a_failed_layer_never_an_exception(frigid, error, status, fragment):
    frigid(error=error)
    layer = frigid_run()
    assert layer.status == status and fragment in layer.note and layer.candidates == [] and layer.seconds >= 0


def test_frigid_works_without_the_install_when_the_engine_says_so(monkeypatch, tmp_path):
    """The real fe.predict with no FRIGID in place: FrigidError -> a failed layer with the setup command in the note."""
    monkeypatch.setenv("MZDESK_FRIGID_DIR", str(tmp_path / "nothing"))
    layer = frigid_run()
    assert layer.status == "failed" and fe.SETUP_COMMAND in layer.note


# ---------------------------------------------------------------- real runs (opt-in)
@pytest.mark.skipif(not RUN_MSFINDER or NO_MSFINDER is not None,
                    reason=f"a real MS-FINDER step run is opt-in: set MZDESK_RUN_MSFINDER_STEP=1 (MS-FINDER: {NO_MSFINDER or 'ready'})")
def test_real_msfinder_step_finds_caffeine_with_the_formula_list():
    layer, summary = msfinder_step.run_msfinder_layer(195.0877, "[M+H]+", CAFFEINE_PEAKS, ["C8H10N4O2", "C8H18OS2"], tol_ppm=5.0, tol_da=0.01)
    assert layer.status == "done", layer.note
    top = layer.candidates[0]
    assert top.ik14 == CAFFEINE_IK14 and top.layer_score > 7.5 and summary["good_enough"] is True and summary["margin"] > 0.3
    assert all(c.formula in ("C8H10N4O2", "C8H18OS2") for c in layer.candidates) and layer.seconds < 90


@pytest.mark.skipif(not RUN_FRIGID or fe.status()["python"] is None, reason="a real FRIGID step run is opt-in: set MZDESK_RUN_FRIGID=1 and install FRIGID")
def test_real_frigid_step_generates_caffeine():
    layer = frigid_run()
    assert layer.status == "done", layer.note
    assert CAFFEINE_IK14 in [c.ik14 for c in layer.candidates] and all(c.hypothesis and c.formula == "C8H10N4O2" for c in layer.candidates)
