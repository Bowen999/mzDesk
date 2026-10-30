"""Isomer check (mzdesk/exclusivity/isomers.py on MS-FINDER candidate runs) and the Isomers table of Alternative explanations.

Unit tests replay recorded MS-FINDER output or hand-made fragment lists; no MS-FINDER starts. Fixtures in tests/fixtures/isomers/ come from
real runs (2026-10-28): kaempferol / luteolin on MassBank MSBNK-LCSB-LU089855 (CC BY) with only the two candidates, and the page example
(caffeine, MassBank MSBNK-LCSB-LU119905, CC BY 4.0) with MS-FINDER's local databases on.
Real runs are opt-in: MZDESK_RUN_ISOMERS=1 (one MS-FINDER run each, about 5-12 s and 2.3 GB, one at a time). They read their spectra
from the local MassBank / GNPS indexes and skip without them.
"""
import json
import os
import sqlite3
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from mzdesk.api.app import app
from mzdesk.exclusivity import isomers
from mzdesk.ms2 import msfinder_engine as mf
from mzdesk.spec import CAFFEINE, CAFFEINE_MSFINDER_MS2, MODULES

FIX = Path(__file__).parent / "fixtures" / "isomers"
KAEMPFEROL = json.loads((FIX / "kaempferol_lcsb75.json").read_text())
KAEMPFEROL_SFD = (FIX / "kaempferol_luteolin_C15H10O6.sfd").read_text()
CAFFEINE_RUN = json.loads((FIX / "caffeine_local.json").read_text())["explain_candidates"]
CAFFEINE_PEAKS = [tuple(map(float, ln.split())) for ln in CAFFEINE_MSFINDER_MS2.splitlines()]
EXAMPLE = next(m for m in MODULES if m["id"] == "exclusivity")["example"]

LEU, ILE, NLE = "CC(C)C[C@H](N)C(=O)O", "CC[C@H](C)[C@H](N)C(=O)O", "CCCC[C@H](N)C(=O)O"          # C6H13NO2, [M+H]+ 132.1019
LEU_KEY, ILE_KEY, NLE_KEY = "ROHFNLRQFUQHCH-YFKPBYRVSA-N", "AGPKZVBTJJNPAG-WHFBIAKZSA-N", "LRQKBLKVPFOOQJ-YFKPBYRVSA-N"

client = TestClient(app)


def _unavailable():
    try:
        mf.check_msfinder()
    except mf.MsFinderError as exc:
        return str(exc)
    return None


NO_MSFINDER = _unavailable()
RUN_REAL = os.environ.get("MZDESK_RUN_ISOMERS") == "1"


def _structure(inchikey, mzs, score=7.0, name="x"):
    """A hand-made scored structure: it explains exactly the peaks at `mzs`."""
    return {"name": name, "inchikey": inchikey, "smiles": "C", "formula": "C6H13NO2", "score": score,
            "fragments": [{"mz": mz, "intensity": 1.0, "formula": f"F{mz:.0f}", "smiles": f"S{mz:.0f}", "ppm": 0.0} for mz in mzs]}


@pytest.fixture
def engine(monkeypatch):
    """Replace mf.explain_candidates. install(structures) answers with them: a structure without "id" is matched to a candidate by
    InChIKey (a database structure otherwise). The returned list collects the arguments of each call."""
    calls = []

    def install(structures=(), formulas=(), error=None):
        def explain(precursor_mz, positive, adduct, ms2, candidates, **kw):
            calls.append({"precursor_mz": precursor_mz, "positive": positive, "adduct": adduct, "ms2": ms2, "candidates": candidates, **kw})
            if error:
                raise error
            ids = {c["inchikey"]: c["id"] for c in candidates}
            out = [{**s, "id": s["id"] if "id" in s else ids.get(s["inchikey"])} for s in structures]
            return {"adduct": adduct, "structures": out, "formulas": list(formulas), "seconds": 4.2, "peak_rss_mb": 2200.0, "log": ""}
        monkeypatch.setattr(mf, "explain_candidates", explain)
        return calls
    return install


def _kaempferol_structures():
    return [{**s, "formula": "C15H10O6"} for s in mf.parse_sfd(KAEMPFEROL_SFD)]


# ---------------------------------------------------------------- MS-FINDER structure files
def test_structure_file_of_a_real_run_gives_each_candidate_its_explained_peaks():
    rows = mf.parse_sfd(KAEMPFEROL_SFD)
    assert [r["name"] for r in rows] == ["mzdesk_candidate_0", "mzdesk_candidate_1"]
    assert [r["inchikey"][:14] for r in rows] == ["IYRMWMYZSQPJKC", "IQPNAANSBPBGFQ"] and all(r["smiles"] for r in rows)
    assert rows[0]["score"] == pytest.approx(6.6614) and rows[1]["score"] == pytest.approx(6.6048)
    frag = next(f for f in rows[0]["fragments"] if f["mz"] == pytest.approx(195.0451))
    assert frag["formula"] == "C13H7O2" and frag["smiles"] == "OC1=CC=CC(OCC=2C=CC=CC2)=C1"        # C13H12O2 with 5 rearranged H removed
    assert not any(f["mz"] == pytest.approx(195.0451) for f in rows[1]["fragments"])                   # luteolin does not explain it
    assert rows[0]["fragments"][-1]["mz"] == pytest.approx(285.0405)                                    # the precursor row is part of the list


def test_structure_file_formula_follows_the_rearranged_hydrogens():
    head = "NAME: a\nINCHIKEY: K\nSMILES: C\nTotalScore: 1.5\nNum Fragment VS2 (M/Z Intensity Formula RearrangedHydrogen SMILES AdductString): 3\n"
    rows = mf.parse_sfd(head + "100.1\t5\tC7H6O4\t-1\tO\t\n110.2\t6\tCO2\t1\tO=C=O\t\n120.3\t7\tC2H4\t0\tC=C\tNa\nNAME: b\n")
    assert [f["formula"] for f in rows[0]["fragments"]] == ["C7H5O4", "CHO2", "C2H4 (Na)"] and rows[0]["score"] == 1.5
    assert rows[1] == {"name": "b", "inchikey": "", "smiles": "", "score": None, "fragments": []}
    assert mf.parse_sfd("") == []


def test_ini_with_values_sets_any_key_and_keeps_every_other_line():
    text = "#comment\r\nIsUserDefinedDB=False\r\nUserDefinedDbFilePath=\r\nHMDB=True\r\n"
    out = mf.ini_with_values(text, {"IsUserDefinedDB": "True", "UserDefinedDbFilePath": Path("/tmp/db.txt"), "hmdb": False, "TryTopNMolecularFormulaSearch": 20})
    assert out == "#comment\r\nIsUserDefinedDB=True\r\nUserDefinedDbFilePath=/tmp/db.txt\r\nhmdb=False\r\nTryTopNMolecularFormulaSearch=20\n"
    assert mf.ini_with_values(text, {}) == text


@pytest.fixture
def fake_install(monkeypatch, tmp_path):
    """A fake install with the shipped parameter file, and a stand-in for the console run that writes the given result files."""
    exe, ini = tmp_path / "MsfinderConsoleApp.exe", tmp_path / "x.ini"
    exe.write_text("")
    ini.write_text(mf.INI_PATH.read_text())
    monkeypatch.setattr(mf, "status", lambda: {"mono": "/usr/bin/mono", "dir": tmp_path, "exe": exe, "resources": tmp_path, "ini": ini})
    seen = {}

    def install(sfd=KAEMPFEROL_SFD, formula_text="Formula\tFormula score\nC15H10O6\t4.2\n"):
        def fake_run(cmd, cwd, env, timeout, rss_cap_mb, **kw):
            in_dir, out_dir = Path(cmd[cmd.index("-i") + 1]), Path(cmd[cmd.index("-o") + 1])
            seen["ini"] = Path(cmd[cmd.index("-m") + 1]).read_text()
            seen["db_path"] = next(ln.split("=", 1)[1] for ln in seen["ini"].splitlines() if ln.startswith("UserDefinedDbFilePath="))
            seen["db"] = Path(seen["db_path"]).read_text()
            seen["query"] = (in_dir / "query.mat").read_text()
            seen["tmp"] = in_dir.parent
            (out_dir / "Formula result-2099.txt").write_text(formula_text)
            if sfd is not None:
                (in_dir / "query").mkdir()
                (in_dir / "query" / "C15H10O6.sfd").write_text(sfd)
            return 1, "done", 2200.0
        monkeypatch.setattr(mf, "run_watched", fake_run)
        return seen
    return install


CANDS = [{"id": "k", "smiles": "O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12", "inchikey": "IYRMWMYZSQPJKC-UHFFFAOYSA-N", "formula": "C15H10O6", "mass": 286.047738},
         {"id": "l", "smiles": "O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12", "inchikey": "IQPNAANSBPBGFQ-UHFFFAOYSA-N", "formula": "C15H10O6", "mass": 286.047738}]


def _ini_values(text):
    return dict(ln.split("=", 1) for ln in text.splitlines() if "=" in ln and not ln.startswith("#"))


def test_candidate_run_scores_only_the_given_structures_through_a_user_database(fake_install):
    seen = fake_install()
    res = mf.explain_candidates(285.0405, False, "[M-H]-", [(195.0451, 6.7), (285.0405, 100.0)], CANDS, tol_da=0.02)
    values = _ini_values(seen["ini"])
    assert values["IsUserDefinedDB"] == "True" and values["Ms2Tolerance"] == "0.02" and values["TryTopNMolecularFormulaSearch"] == "20"
    assert all(values[k] == "False" for k in mf.STRUCTURE_DB_KEYS) and values["MinesNeverUse"] == "True"
    assert values["Fcheck"] == "False" and values["PubChemNeverUse"] == "True"                  # no halogen in the formulas; nothing online
    db = seen["db"].splitlines()
    assert db[0].split("\t") == list(mf.USER_DB_COLUMNS)
    assert db[1].split("\t") == ["mzdesk_candidate_0", CANDS[0]["inchikey"], "IYRMWMYZSQPJKC", "N/A", "286.047738", "C15H10O6", CANDS[0]["smiles"], "N/A"]
    assert "PRECURSORTYPE: [M-H]-" in seen["query"] and "PRECURSORMZ: 285.04050" in seen["query"]
    assert [(s["id"], s["formula"]) for s in res["structures"]] == [("k", "C15H10O6"), ("l", "C15H10O6")]
    assert res["formulas"][0]["formula"] == "C15H10O6" and res["adduct"] == "[M-H]-"
    assert not seen["tmp"].exists() and mf.INI_PATH.read_text() == Path(mf.status()["ini"]).read_text()     # temp folder gone, shipped file as it was


def test_candidate_run_with_local_databases_keeps_them_on_and_follows_halogens(fake_install):
    seen = fake_install()
    chloro = [{**CANDS[0], "formula": "C15H9ClO6"}]
    mf.explain_candidates(285.0405, False, "[M-H]-", [(195.0451, 6.7)], chloro, local_db=True)
    values = _ini_values(seen["ini"])
    assert values["HMDB"] == "True" and values["MinesNeverUse"] == "False" and values["IsUserDefinedDB"] == "True"
    assert values["ClCheck"] == "True" and values["SiCheck"] == "False" and values["Ms2Tolerance"] == "0.01"


def test_candidate_run_without_structure_file_returns_no_structure(fake_install):
    fake_install(sfd=None)
    assert mf.explain_candidates(285.0405, False, "[M-H]-", [(195.0451, 6.7)], CANDS)["structures"] == []


@pytest.mark.parametrize("args,match", [((285.0405, False, "[M-H]-", [(1.0, 1.0)], []), "No candidate"),
                                        ((0, False, "[M-H]-", [(1.0, 1.0)], CANDS), "Precursor"),
                                        ((285.0405, False, "[M+H]+", [(1.0, 1.0)], CANDS), "negative-mode")])
def test_candidate_run_rejects_input_it_cannot_use(fake_install, args, match):
    fake_install()
    with pytest.raises(ValueError, match=match):
        mf.explain_candidates(*args)


def test_candidate_run_is_turned_away_while_another_heavy_run_holds_the_lock(fake_install, monkeypatch):
    fake_install()
    monkeypatch.setattr(mf, "BUSY_WAIT_S", 0.05)
    with mf._RUN_LOCK:
        with pytest.raises(mf.MsFinderError, match="busy"):
            mf.explain_candidates(285.0405, False, "[M-H]-", [(195.0451, 6.7)], CANDS)
    assert mf.explain_candidates(285.0405, False, "[M-H]-", [(195.0451, 6.7)], CANDS)["structures"]


# ---------------------------------------------------------------- the rule on recorded MS-FINDER output
def test_kaempferol_is_a_single_answer_against_luteolin(engine):
    calls = engine(_kaempferol_structures())
    res = isomers.discriminate(KAEMPFEROL["candidates"], KAEMPFEROL["peaks"], KAEMPFEROL["precursor_mz"], "[M-H]-")
    assert (res["verdict"], res["top_id"], res["group"], res["challengers"]) == ("single", "kaempferol", ["kaempferol"], [])
    (iso,) = res["isomers"]
    assert iso["id"] == "luteolin" and iso["distinguished"] is True and iso["isomer_only_peaks"] == []
    assert iso["distinguishing_peaks"] == [{"mz": 195.0451, "intensity": pytest.approx(0.06707), "rel": 6.71, "formula": "C13H7O2",
                                            "smiles": "OC1=CC=CC(OCC=2C=CC=CC2)=C1"}]
    assert res["peaks_considered"] == 103 and [len(c["explained"]) for c in res["candidates"]] == [40, 39]
    assert res["candidates"][0]["explained_fraction"] == pytest.approx(0.5086, abs=1e-4) and res["candidates"][0]["score"] == pytest.approx(6.6614)
    assert res["notes"][0].startswith("Relative to the 2 candidates considered") and res["engine"] == "MS-FINDER 3.73"
    # what MS-FINDER was given: the theoretical m/z, MS-FINDER's adduct name, peaks >= 1 % (precursor kept), the two prepared structures
    call = calls[0]
    assert call["precursor_mz"] == pytest.approx(285.0405, abs=1e-4) and call["positive"] is False and call["adduct"] == "[M-H]-"
    assert len(KAEMPFEROL["peaks"]) == 107 and len(call["ms2"]) == 104 and all(i >= 0.01 for _, i in call["ms2"])     # 104 of 107 are >= 1 %
    assert call["tol_da"] == 0.01 and call["local_db"] is False
    assert [c["inchikey"][:14] for c in call["candidates"]] == ["IYRMWMYZSQPJKC", "IQPNAANSBPBGFQ"] and call["candidates"][0]["formula"] == "C15H10O6"
    json.dumps(res, allow_nan=False)


def test_the_same_spectrum_with_luteolin_ranked_first_cannot_separate_them(engine):
    engine(_kaempferol_structures())
    cands = KAEMPFEROL["candidates"][::-1]
    res = isomers.discriminate(cands, KAEMPFEROL["peaks"], KAEMPFEROL["precursor_mz"], "[M-H]-")
    assert (res["verdict"], res["group"]) == ("indistinguishable", ["luteolin", "kaempferol"])
    assert [(c["id"], [p["mz"] for p in c["peaks"]]) for c in res["challengers"]] == [("kaempferol", [195.0451])]


def test_page_example_caffeine_against_local_database_isomers(engine):
    calls = engine(CAFFEINE_RUN["structures"], CAFFEINE_RUN["formulas"])
    res = isomers.discriminate([{"id": "compound", "smiles": CAFFEINE, "name": "caffeine"}], CAFFEINE_PEAKS, 195.0877, "[M+H]+", tol_ppm=5,
                               local_isomers=True)
    assert calls[0]["local_db"] is True and res["verdict"] == "indistinguishable"
    assert len(res["isomers"]) == 9 and all(c["source"] == "MS-FINDER database" for c in res["candidates"][1:])
    assert not any(i["id"].startswith("msfinder:RYYVLZVUVIJVGH") for i in res["isomers"])          # caffeine's own database entry is not an isomer
    apart = {i["name"] for i in res["isomers"] if i["distinguished"]}
    assert apart == {"Enprofylline", "UNPD100701", "N-(6-acetamidopyrimidin-4-yl)acetamide"} and len(res["group"]) == 7
    assert res["group"][0] == "compound" and "1,3,9-Trimethylxanthin (Isocoffein)" in {i["name"] for i in res["isomers"] if not i["distinguished"]}
    assert any("local structure databases: the 9 best-scoring" in n for n in res["notes"])


# ---------------------------------------------------------------- the rule on hand-made fragments ([M+H]+ of C6H13NO2 = 132.1019)
PEAKS = [(44.0495, 100.0), (57.0573, 5.0), (69.0699, 6.0), (86.0964, 50.0), (101.0, 3.0), (132.1019, 30.0), (133.1052, 3.0), (150.0, 2.0)]
LEU_ILE = [{"id": "leu", "smiles": LEU, "name": "leucine"}, {"id": "ile", "smiles": ILE, "name": "isoleucine"}]


@pytest.mark.parametrize("rel,verdict", [(5.0, "single"), (4.9, "indistinguishable")])
def test_a_distinguishing_peak_must_be_at_least_5_percent_of_the_base_peak(engine, rel, verdict):
    engine([_structure(LEU_KEY, [44.0495, 57.0573, 86.0964]), _structure(ILE_KEY, [44.0495, 86.0964])])
    peaks = [(mz, rel if mz == 57.0573 else i) for mz, i in PEAKS]
    res = isomers.discriminate(LEU_ILE, peaks, 132.1019, "[M+H]+")
    assert res["verdict"] == verdict and res["group"] == (["leu"] if verdict == "single" else ["leu", "ile"])


def test_precursor_its_isotopes_and_heavier_peaks_never_count(engine):
    engine([_structure(LEU_KEY, [44.0495, 132.1019, 133.1052, 150.0]), _structure(ILE_KEY, [44.0495])])
    res = isomers.discriminate(LEU_ILE, PEAKS, 132.1019, "[M+H]+")
    assert res["peaks_considered"] == 5 and res["verdict"] == "indistinguishable"                   # 44, 57, 69, 86, 101; the rest is out
    assert [p["mz"] for p in res["candidates"][0]["explained"]] == [44.0495]
    assert res["candidates"][0]["explained_fraction"] == round(100 / 164, 4)


def test_challengers_are_information_only(engine):
    engine([_structure(LEU_KEY, [44.0495, 57.0573]), _structure(ILE_KEY, [44.0495, 69.0699])])
    res = isomers.discriminate(LEU_ILE, PEAKS, 132.1019, "[M+H]+")
    assert res["verdict"] == "single" and res["group"] == ["leu"]
    assert res["challengers"] == [{"id": "ile", "name": "isoleucine", "peaks": [{"mz": 69.0699, "intensity": 6.0, "rel": 6.0, "formula": "F69", "smiles": "S69"}]}]


def test_an_isomer_msfinder_did_not_score_is_not_ruled_out(engine):
    engine([_structure(LEU_KEY, [44.0495, 57.0573]), _structure(ILE_KEY, [44.0495])])
    res = isomers.discriminate(LEU_ILE + [{"id": "nle", "smiles": NLE, "name": "norleucine"}], PEAKS, 132.1019, "[M+H]+")
    assert res["verdict"] == "indistinguishable" and res["group"] == ["leu", "nle"]
    nle = res["isomers"][1]
    assert nle["evaluated"] is False and nle["distinguished"] is None and res["candidates"][2]["explained_fraction"] is None
    assert any("did not score norleucine" in n for n in res["notes"])


@pytest.mark.parametrize("structures,formulas,why", [
    ([_structure(ILE_KEY, [44.0495])], [], "C6H13NO2 is not among the 20 formulas"),
    ([_structure(ILE_KEY, [44.0495])], [{"formula": "C6H13NO2"}], "could not fragment"),
    ([_structure(LEU_KEY, [44.0495])], [{"formula": "C6H13NO2"}], "did not score any of the isomers"),
])
def test_missing_scores_give_no_verdict(engine, structures, formulas, why):
    engine(structures, formulas)
    res = isomers.discriminate(LEU_ILE, PEAKS, 132.1019, "[M+H]+")
    assert res["verdict"] == "insufficient" and res["group"] == [] and why in res["notes"][0]


def test_an_engine_failure_gives_no_verdict_and_keeps_the_message(engine):
    engine(error=mf.MsFinderError("MS-FINDER is busy with another request.", log="tail"))
    res = isomers.discriminate(LEU_ILE, PEAKS, 132.1019, "[M+H]+")
    assert res["verdict"] == "insufficient" and res["error"] == "MS-FINDER is busy with another request." and res["log"] == "tail"
    assert res["notes"][0] == "MS-FINDER could not run: MS-FINDER is busy with another request." and res["top_id"] == "leu"


@pytest.mark.parametrize("cands,peaks,precursor,adduct,why", [
    (LEU_ILE, [], 132.1019, "[M+H]+", "No MS2 peaks."),
    (LEU_ILE, [(132.1019, 100.0), (133.1052, 7.0), (44.0, 0.5)], 132.1019, "[M+H]+", "No fragment peak"),
    ([], PEAKS, 132.1019, "[M+H]+", "No candidate structure."),
    (LEU_ILE, PEAKS, 132.1019, "Auto", "cannot use the adduct"),
    (LEU_ILE, PEAKS, 132.1019, None, "cannot use the adduct (none given)"),
    (LEU_ILE, PEAKS, 132.1119, "[M+H]+", "does not match the precursor m/z 132.1119 within 10 ppm"),
    (LEU_ILE[:1], PEAKS, 132.1019, "[M+H]+", "No isomer to compare with"),
    ([{"id": "x", "smiles": "not a smiles"}] + LEU_ILE, PEAKS, 132.1019, "[M+H]+", "The top candidate cannot be used: the SMILES cannot be read."),
    ([{"id": "q", "smiles": "C[N+](C)(C)CCO"}] + LEU_ILE, PEAKS, 132.1019, "[M+H]+", "keeps a net charge"),
])
def test_no_verdict_without_what_the_rule_needs_and_no_engine_run(engine, cands, peaks, precursor, adduct, why):
    calls = engine()
    res = isomers.discriminate(cands, peaks, precursor, adduct)
    assert res["verdict"] == "insufficient" and why in res["notes"][0] and calls == []
    json.dumps(res, allow_nan=False)


def test_only_same_formula_unique_structures_are_compared_up_to_the_cap(engine):
    calls = engine([_structure(LEU_KEY, [44.0495, 57.0573]), _structure(ILE_KEY, [44.0495])])
    cands = LEU_ILE + [{"id": "dup", "smiles": "CC(C)C[C@H](N)C(O)=O", "name": "leucine again"}, {"id": "caf", "smiles": CAFFEINE},
                       {"id": "bad", "smiles": "C1CC", "name": "broken"}, {"id": "nle", "smiles": NLE}]
    res = isomers.discriminate(cands, PEAKS, 132.1019, "[M+H]+", max_candidates=2)
    assert [c["id"] for c in calls[0]["candidates"]] == ["leu", "ile"] and res["verdict"] == "single"
    notes = " ".join(res["notes"])
    assert "leucine again is not compared: the same structure" in notes and "broken is not compared: the SMILES cannot be read" in notes
    assert "1 candidate(s) with another formula than C6H13NO2" in notes and "1 lower-ranked isomer(s) are not compared: at most 2" in notes


def test_the_adduct_is_passed_in_msfinder_spelling_and_formate_works(engine):
    calls = engine([_structure(LEU_KEY, [44.0495, 57.0573]), _structure(ILE_KEY, [44.0495])])
    res = isomers.discriminate(LEU_ILE, PEAKS, None, "[M+HCOO]-")
    assert calls[0]["adduct"] == "[M+FA-H]-" and calls[0]["positive"] is False
    assert calls[0]["precursor_mz"] == pytest.approx(176.0928, abs=1e-4) and res["verdict"] == "single"


def test_only_msfinder_is_an_engine():
    with pytest.raises(ValueError, match="MS-FINDER"):
        isomers.discriminate(LEU_ILE, PEAKS, 132.1019, "[M+H]+", engine="metfrag")


# ---------------------------------------------------------------- the page (engine replaced)
def _run(inputs):
    return client.post("/api/modules/exclusivity/run", json=inputs)


def test_page_is_live_open_and_has_three_tables():
    import re
    mod = next(m for m in MODULES if m["id"] == "exclusivity")
    lock_ts = (Path(__file__).parents[1] / "web" / "src" / "lock.ts").read_text()
    locked_subgroups = re.search(r"LOCKED_SUBGROUPS = new Set\(\[([^\]]*)\]\)", lock_ts).group(1)
    assert mod["live"] and mod["subgroup"] == "tools" and '"tools"' not in locked_subgroups              # open: only DMPK and ICP-MS are locked
    assert [b.get("id") for b in mod["outputs"]] == ["same_formula", "isobaric", "isomers"]


def test_page_example_fills_the_isomers_table_and_says_what_cannot_be_separated(engine):
    calls = engine(CAFFEINE_RUN["structures"], CAFFEINE_RUN["formulas"])
    body = _run(EXAMPLE).json()
    assert body["status"] == "ok" and set(body["outputs"]) == {"same_formula", "isobaric", "isomers"} and calls[0]["local_db"] is True
    assert len(body["outputs"]["same_formula"]["rows"]) == 3 and body["outputs"]["isobaric"]["rows"] == []
    table = body["outputs"]["isomers"]
    rows = table["rows"]
    assert len(rows) == 10 and table["count"] == 9 and rows[0][0] == "caffeine" and rows[0][3] == {"badge": "your compound", "tone": "neutral"}
    assert all(r[1]["text"] and "<svg" in r[1]["hover_svg"] for r in rows)                      # hovering a SMILES draws its 2D structure
    assert [r[3]["badge"] for r in rows[1:]].count("no") == 6 and [r[3]["badge"] for r in rows[1:]].count("yes") == 3
    enprofylline = next(r for r in rows if r[0] == "Enprofylline")
    assert enprofylline[4]["text"] == "83.0603 (31.6 %) · 109.0396 (9.0 %)" and len(enprofylline[4]["help"]) == 2
    assert body["verdict"].startswith("<b>These cannot be separated</b> by this spectrum: <b>caffeine</b> and cyclo(gly-his)")
    assert "and 3 more (of 9 isomers considered)" in body["verdict"] and table["help"][0].startswith("Relative to the 10 candidates")


def test_page_single_answer(engine):
    engine(_kaempferol_structures())
    peaks = "\n".join(f"{mz} {i}" for mz, i in KAEMPFEROL["peaks"])
    body = _run({"structure": {"mode": "SMILES / InChI", "value": KAEMPFEROL["candidates"][0]["smiles"]},
                 "precursor": {"precursor_mz": "285.0405", "ion_mode": "Negative", "adduct": "Auto"},
                 "ms2": {"source": "Paste", "value": peaks}}).json()
    assert body["verdict"] == "<b>C15H10O6</b> is the single answer: the one isomer considered fails to explain at least one of its fragments."


def test_page_without_ms2_gives_no_verdict_and_a_bare_mz_still_runs(engine):
    engine()
    body = _run({"structure": {"mode": "SMILES / InChI", "value": CAFFEINE}, "precursor": {"precursor_mz": "195.0877", "ion_mode": "Positive"}}).json()
    assert body["status"] == "ok" and body["outputs"]["isomers"]["rows"] == []
    assert body["verdict"].endswith("<b>No isomer verdict</b>: No MS2 peaks: the isomer check needs an MS2 spectrum.")
    r = _run({"precursor": {"precursor_mz": "195.0877"}, "ms2": {"source": "Paste", "value": "138.0662 100"}})
    assert r.status_code == 200 and r.json()["outputs"]["same_formula"]["rows"] == []
    r = _run({"ms2": {"source": "Paste", "value": "138.0662 100"}})
    assert r.status_code == 400 and "compound or formula" in r.json()["detail"]


def test_page_engine_failure_is_502_with_the_log(engine):
    engine(error=mf.MsFinderError("MS-FINDER is not installed.", log="no mono"))
    r = _run(EXAMPLE)
    assert r.status_code == 502 and r.json() == {"detail": "MS-FINDER is not installed.", "log": "no mono"}


# ---------------------------------------------------------------- the real thing (opt-in)
real = pytest.mark.skipif(not RUN_REAL or NO_MSFINDER is not None,
                          reason=f"real MS-FINDER isomer runs are opt-in: set MZDESK_RUN_ISOMERS=1 (MS-FINDER: {NO_MSFINDER or 'ready'})")


def _library_spectrum(library, accession):
    """(precursor m/z, peaks, SMILES) of one record of a local library index; skip when the index is not installed."""
    from mzdesk.ms2 import library_index
    if not library_index.is_installed(library):
        pytest.skip(f"{library} index not installed (python -m mzdesk.ms2.library_setup {library})")
    lib = library_index.load(library)
    con = sqlite3.connect(f"{(lib.path / 'meta.sqlite').resolve().as_uri()}?mode=ro", uri=True)
    try:
        row = con.execute("SELECT id, smiles FROM meta WHERE accession = ?", (accession,)).fetchone()
    finally:
        con.close()
    if row is None:
        pytest.skip(f"{accession} is not in the local {library} index")
    mz, inten = lib.peaks(row[0])
    return float(lib.prec_mz[row[0]]), list(zip(mz.tolist(), inten.tolist())), row[1]


@real
def test_real_kaempferol_is_told_apart_from_luteolin_but_not_the_other_way_round():
    precursor, peaks, _ = _library_spectrum("MassBank", "MSBNK-LCSB-LU089855")           # kaempferol [M-H]-, CE 75, CC BY
    res = isomers.discriminate(KAEMPFEROL["candidates"], peaks, precursor, "[M-H]-")
    assert res["verdict"] == "single" and 195.0451 in [round(p["mz"], 4) for p in res["isomers"][0]["distinguishing_peaks"]]
    back = isomers.discriminate(KAEMPFEROL["candidates"][::-1], peaks, precursor, "[M-H]-")
    assert back["verdict"] == "indistinguishable" and [c["id"] for c in back["challengers"]] == ["kaempferol"]


@real
@pytest.mark.parametrize("accession", ["CCMSLIB00012361594", "CCMSLIB00012360432"])               # MSnLib leucine / isoleucine [M+H]+
def test_real_leucine_and_isoleucine_cannot_be_separated(accession):
    precursor, peaks, _ = _library_spectrum("GNPS", accession)
    for cands in (LEU_ILE, LEU_ILE[::-1]):
        res = isomers.discriminate(cands, peaks, precursor, "[M+H]+")
        assert res["verdict"] == "indistinguishable" and set(res["group"]) == {"leu", "ile"}, (accession, cands[0]["id"])


@real
def test_real_page_example():
    body = _run(EXAMPLE).json()
    assert body["status"] == "ok" and body["verdict"] and len(body["outputs"]["isomers"]["rows"]) >= 2
