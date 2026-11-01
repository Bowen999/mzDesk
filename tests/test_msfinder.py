"""Structure prediction（MS-FINDER，Mono）的测试：查询文件写法、结果解析、安装检查、进程看护、runner 外形、安装脚本。

fixtures：tests/fixtures/msfinder/ 是真实运行的结果文件（咖啡因、卡马西平，MassBank MSBNK-LCSB-LU119905 / LU126603，CC BY），路径换成了中性路径。
最后的集成测试（test_real_run_*）真的启动 MS-FINDER（每次约 4-5 s，峰值约 2.3 GB 内存，一次一个）；没有 Mono / MS-FINDER 时它们跳过（python -m mzdesk.ms2.msfinder_setup）。
卤素 / Si 开关的两个集成测试用 MassBank 的真实谱（出处和许可见各常量上方的注释），没有编造峰。
"""
import os
import re
import sys
import threading
import zipfile
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from mzdesk.api.app import app
from mzdesk.ms2 import msfinder_engine as mf
from mzdesk.ms2 import msfinder_setup as su
from mzdesk.spec import CAFFEINE_MS1, CAFFEINE_MSFINDER_MS2, MODULES, MSFINDER_ADDUCTS

FIX = Path(__file__).parent / "fixtures" / "msfinder"
CAFFEINE_FORMULA = (FIX / "caffeine_formula.txt").read_text()
CAFFEINE_STRUCTURE = (FIX / "caffeine_structure.txt").read_text()
CARBAMAZEPINE_FORMULA = (FIX / "carbamazepine_formula.txt").read_text()
CARBAMAZEPINE_STRUCTURE = (FIX / "carbamazepine_structure.txt").read_text()
HEADER_ONLY = CAFFEINE_FORMULA.splitlines()[0] + "\n"          # what MS-FINDER writes when nothing was found
STRUCTURE_HEADER_ONLY = CAFFEINE_STRUCTURE.splitlines()[0] + "\n"

client = TestClient(app)
PEAKS = [(56.0495, 88.0), (110.0712, 590.0), (138.0662, 999.0), (195.0876, 58.0)]
PRECURSOR = {"precursor_mz": "195.0877", "ion_mode": "Positive", "adduct": "Auto"}


def _unavailable():
    try:
        mf.check_msfinder()
    except mf.MsFinderError as exc:
        return str(exc)
    return None


NO_MSFINDER = _unavailable()


def _paste(peaks):
    return {"source": "Paste", "value": peaks}


def _run(inputs):
    return client.post("/api/modules/ms2_structure/run", json=inputs)


def _fake_result(**over):
    f, s = mf.parse_formula_results(CAFFEINE_FORMULA), mf.parse_structure_results(CAFFEINE_STRUCTURE)
    res = {"adduct": "[M+H]+", "precursor_mz": 195.0877, "positive": True, "halogen": False, "silicon": False, "formulas": f["rows"], "structures": s["rows"],
           "ms1_count": 0, "ms2_count": 27, "seconds": 4.0, "peak_rss_mb": 2200.0, "log": ""}
    res.update(over)
    return res


@pytest.fixture
def fake_predict(monkeypatch):
    """Replace the engine call; the returned list collects the arguments of each call."""
    calls = []

    def predict(precursor_mz, positive, adduct, ms2, ms1=None, **kw):
        calls.append({"precursor_mz": precursor_mz, "positive": positive, "adduct": adduct, "ms2": ms2, "ms1": ms1,
                      "halogen": kw.get("halogen"), "silicon": kw.get("silicon")})
        return _fake_result(adduct=mf.resolve_adduct(positive, adduct), positive=positive)

    monkeypatch.setattr(mf, "predict", predict)
    return calls


# ---------------------------------------------------------------- query file
def test_query_writes_precursor_type_ion_mode_and_peaks():
    text = mf.build_query(195.0877, True, "[M+H]+", PEAKS)
    lines = text.splitlines()
    assert "PRECURSORTYPE: [M+H]+" in lines and "IONMODE: Positive" in lines and "PRECURSORMZ: 195.08770" in lines
    assert lines[lines.index("MSTYPE: MS2") + 1] == "Num Peaks: 4"
    assert lines[-4:] == ["56.04950\t88.0881", "110.07120\t590.5906", "138.06620\t1000.0000", "195.08760\t58.0581"]      # base peak = 1000, sorted
    assert "MSTYPE: MS1" not in text


def test_query_ms1_block_comes_first_and_counts_match():
    text = mf.build_query(195.0877, True, "[M+H]+", PEAKS, ms1=[(195.0877, 100), (196.0902, 10.3), (197.0922, 0.9)])
    lines = text.splitlines()
    i1, i2 = lines.index("MSTYPE: MS1"), lines.index("MSTYPE: MS2")
    assert i1 < i2 and lines[i1 + 1] == "Num Peaks: 3" and lines[i2 + 1] == "Num Peaks: 4"
    assert lines[i1 + 2:i2] == ["195.08770\t1000.0000", "196.09020\t103.0000", "197.09220\t9.0000"]
    assert len(lines) == i2 + 2 + 4


def test_query_negative_mode():
    text = mf.build_query(179.0350, False, "[M-H]-", PEAKS)
    assert "IONMODE: Negative" in text and "PRECURSORTYPE: [M-H]-" in text


def test_query_numbers_are_plain_decimals():
    # the reader of the .mat file knows no exponent except "e+"
    text = mf.build_query(195.0877, True, "[M+H]+", [(100.0, 1e-9), (150.0, 3e-12), (195.0877, 1e12)])
    peak_lines = [ln for ln in text.splitlines() if ln[:1].isdigit()]
    assert not any("e" in ln.lower() for ln in peak_lines)


def test_query_name_cannot_inject_header_lines():
    text = mf.build_query(195.0877, True, "[M+H]+", PEAKS, name="x\nPRECURSORTYPE: [M]+.")
    assert text.count("PRECURSORTYPE:") == 1 and text.splitlines()[0].startswith("NAME: x")


def test_query_without_peaks_is_an_error():
    with pytest.raises(ValueError, match="no peaks"):
        mf.build_query(195.0877, True, "[M+H]+", [(100.0, 0.0)])


# ---------------------------------------------------------------- adducts
def test_auto_adduct_follows_the_ion_mode():
    assert mf.resolve_adduct(True, None) == mf.resolve_adduct(True, "Auto") == "[M+H]+"
    assert mf.resolve_adduct(False, None) == mf.resolve_adduct(False, "Auto") == "[M-H]-"
    assert mf.resolve_adduct(True, "[M+Na]+") == "[M+Na]+"


def test_adduct_must_fit_the_mode_and_the_list():
    with pytest.raises(ValueError, match="not a positive-mode adduct"):
        mf.resolve_adduct(True, "[M-H]-")
    with pytest.raises(ValueError, match="not a negative-mode adduct"):
        mf.resolve_adduct(False, "[M+3H]3+")
    with pytest.raises(ValueError, match="does not support"):
        mf.resolve_adduct(True, "[M+3Na]3+")


def test_page_offers_auto_first_then_the_engine_adducts():
    field = next(f for m in MODULES if m["id"] == "ms2_structure" for f in m["inputs"] if f["type"] == "precursor")
    assert field["adducts"] == MSFINDER_ADDUCTS
    for mode in ("Positive", "Negative"):
        assert field["adducts"][mode] == ["Auto"] + mf.ADDUCTS[mode]


@pytest.mark.skipif(NO_MSFINDER is not None, reason="MS-FINDER not installed")
def test_adducts_exist_in_the_msfinder_resource_lists():
    res = mf.install_dir() / "Resources"
    known = {"Positive": (res / "AdductPositives.apf").read_text().split(), "Negative": (res / "AdductNegatives.anf").read_text().split()}
    for mode, names in mf.ADDUCTS.items():
        assert set(names) <= set(known[mode]), mode


# ---------------------------------------------------------------- parameter file: the element switches
SHIPPED_INI = mf.INI_PATH.read_text()
ELEMENT_KEYS = ("Fcheck", "ClCheck", "BrCheck", "Icheck", "SiCheck")


def _ini_values(text):
    return dict(ln.split("=", 1) for ln in text.splitlines() if "=" in ln and not ln.startswith("#"))


def _other_lines(text, keys):
    return [ln for ln in text.splitlines(keepends=True) if ln.split("=", 1)[0] not in keys]


def test_shipped_parameter_file_keeps_halogens_and_si_off():
    values = _ini_values(SHIPPED_INI)
    assert all(values[k] == "False" for k in ELEMENT_KEYS)
    assert (values["CommonRange"], values["ExtendedRange"], values["ExtremeRange"]) == ("TRUE", "FALSE", "FALSE")      # the range is not touched either
    assert values["PubChemNeverUse"] == "True"                                                                           # nothing goes online


def test_no_switch_returns_the_text_unchanged():
    assert mf.ini_with_elements(SHIPPED_INI) == SHIPPED_INI
    assert mf.ini_with_elements(SHIPPED_INI, halogen=False, silicon=False) == SHIPPED_INI


def test_extended_range_switches_only_the_range_keys():
    text = mf.ini_with_elements(SHIPPED_INI, extended=True)
    values = _ini_values(text)
    assert (values["CommonRange"], values["ExtendedRange"], values["ExtremeRange"]) == ("False", "True", "False")
    assert _other_lines(text, ("CommonRange", "ExtendedRange", "ExtremeRange")) == _other_lines(SHIPPED_INI, ("CommonRange", "ExtendedRange", "ExtremeRange"))


@pytest.mark.parametrize("halogen,silicon,on", [(True, False, {"Fcheck", "ClCheck", "BrCheck", "Icheck"}), (False, True, {"SiCheck"}),
                                                (True, True, set(ELEMENT_KEYS))])
def test_ini_rewrite_flips_only_the_chosen_keys(halogen, silicon, on):
    text = mf.ini_with_elements(SHIPPED_INI, halogen=halogen, silicon=silicon)
    values = _ini_values(text)
    assert {k for k in ELEMENT_KEYS if values[k] == "True"} == on                      # halogens and Si are independent of each other
    assert _other_lines(text, ELEMENT_KEYS) == _other_lines(SHIPPED_INI, ELEMENT_KEYS)  # every other line byte for byte
    assert len(text.splitlines()) == len(SHIPPED_INI.splitlines())                     # nothing added or dropped
    assert [ln.split("=")[0] for ln in text.splitlines()] == [ln.split("=")[0] for ln in SHIPPED_INI.splitlines()]      # same order


def test_ini_rewrite_keeps_line_endings_comments_and_adds_a_missing_key():
    crlf = "#Fcheck=False\r\nFcheck=False\r\nOcheck=True\r\nSiCheck=False"
    out = mf.ini_with_elements(crlf, halogen=True)
    assert out.split("\r\n")[:3] == ["#Fcheck=False", "Fcheck=True", "Ocheck=True"] and out.endswith("SiCheck=False\nClCheck=True\nBrCheck=True\nIcheck=True\n")
    assert mf.ini_with_elements("fcheck=False\n", halogen=True).splitlines()[0] == "Fcheck=True"      # the reader ignores case, so do we
    assert mf.ini_with_elements("", silicon=True) == "SiCheck=True\n"


# ---------------------------------------------------------------- result files
def test_formula_results_of_a_real_caffeine_run():
    res = mf.parse_formula_results(CAFFEINE_FORMULA)
    assert [r["formula"] for r in res["rows"]] == ["C8H10N4O2", "C8H18OS2"]
    top = res["rows"][0]
    assert top["score"] == 4.407 and top["neutral_mass"] == pytest.approx(194.08037556) and top["mass_error_da"] == pytest.approx(-4.8e-05)
    assert (res["ms1_count"], res["ms2_count"]) == (0, 27)


def test_structure_results_of_a_real_caffeine_run():
    res = mf.parse_structure_results(CAFFEINE_STRUCTURE)
    rows = res["rows"]
    assert len(rows) == 15 and (res["ms1_count"], res["ms2_count"]) == (0, 27)
    top = rows[0]
    assert (top["name"], top["score"], top["formula"], top["ontology"]) == ("Caffeine", 8.1594, "C8H10N4O2", "Xanthines")
    assert top["inchikey"] == "RYYVLZVUVIJVGH-UHFFFAOYSA-N" and top["smiles"] == "O=C1C2=C(N=CN2C)N(C(=O)N1C)C"
    assert top["databases"][0] == ("HMDB", "HMDB0001847") and ("FooDB", "FDB002100;FDB017732") in top["databases"] and len(top["databases"]) == 20
    scores = [r["score"] for r in rows]
    assert scores == sorted(scores, reverse=True)
    assert any(r["ontology"] == "" for r in rows)                         # a missing ontology stays empty
    assert next(r for r in rows if r["name"] == "cyclo(gly-his)")["databases"] == [("COCONUT", "CNP0431962")]


def test_structure_results_of_a_real_carbamazepine_run():
    rows = mf.parse_structure_results(CARBAMAZEPINE_STRUCTURE)["rows"]
    assert rows[0]["name"] == "Carbamazepine" and rows[0]["inchikey"].startswith("FFGPTBGBLSHEPO")
    assert "spiro[1,3-dihydroindene-2,3'-imidazo[1,2-a]pyridine]-2'-one" in [r["name"] for r in rows]      # the quotes MS-FINDER adds are dropped
    assert mf.parse_formula_results(CARBAMAZEPINE_FORMULA)["rows"][0]["formula"] == "C15H12N2O"


def test_header_only_files_mean_no_candidates():
    assert mf.parse_formula_results(HEADER_ONLY) == {"rows": [], "ms1_count": None, "ms2_count": None}
    assert mf.parse_structure_results(STRUCTURE_HEADER_ONLY) == {"rows": [], "ms1_count": None, "ms2_count": None}
    assert mf.parse_formula_results("")["rows"] == []


def test_results_are_ranked_by_score_whatever_the_file_order():
    head = "\t".join(["File path", "File name", "Title", "MS1 count", "MSMS count", "Precursor mz", "Precursor type", "Structure rank",
                      "Total score", "Databases", "Formula", "Ontology", "InChIKey", "SMILES"])

    def row(name, score, dbs, formula="C1"):
        return "\t".join(["/x", "q", "q", "3", "10", "100", "[M+H]+", name, score, dbs, formula, "", "KEY-" + name, "C"])

    text = "\n".join([head, row("low", "1.5", "A=1"), row("high", "9.25", "A=1,B=x,y,C=3E-2"), row("sci", "1E+1", "")]) + "\n"
    res = mf.parse_structure_results(text)
    assert [r["name"] for r in res["rows"]] == ["sci", "high", "low"]
    assert res["rows"][1]["databases"] == [("A", "1"), ("B", "x,y"), ("C", "3E-2")]            # a comma inside an ID stays in the ID
    assert (res["ms1_count"], res["ms2_count"]) == (3, 10)


# ---------------------------------------------------------------- installation
def test_missing_install_says_what_to_run(monkeypatch, tmp_path):
    monkeypatch.setenv("MZDESK_MSFINDER_DIR", str(tmp_path / "nothing"))
    monkeypatch.setattr(mf, "find_mono", lambda: sys.executable)
    with pytest.raises(mf.MsFinderError) as exc:
        mf.check_msfinder()
    assert mf.SETUP_COMMAND in str(exc.value) and str(tmp_path / "nothing") in str(exc.value) and "MZDESK_MSFINDER_DIR" in str(exc.value)


def test_missing_mono_names_the_local_version(monkeypatch):
    monkeypatch.setattr(mf, "find_mono", lambda: None)
    with pytest.raises(mf.MsFinderError) as exc:
        mf.check_msfinder()
    msg = str(exc.value)
    assert "Mono" in msg and "local version" in msg and mf.SETUP_COMMAND in msg


def test_mono_comes_from_the_environment_first(monkeypatch):
    monkeypatch.setenv("MZDESK_MONO", sys.executable)
    assert mf.find_mono() == sys.executable
    monkeypatch.setenv("MZDESK_MONO", "/nonexistent/mono")
    assert mf.find_mono() != "/nonexistent/mono"


def test_input_errors_come_before_the_install_check(monkeypatch, tmp_path):
    monkeypatch.setenv("MZDESK_MSFINDER_DIR", str(tmp_path / "nothing"))
    with pytest.raises(ValueError, match="not a positive-mode adduct"):
        mf.predict(195.0877, True, "[M-H]-", PEAKS)
    with pytest.raises(ValueError, match="positive"):
        mf.predict(0.0, True, None, PEAKS)
    with pytest.raises(ValueError, match="at most"):
        mf.predict(195.0877, True, None, [(100.0 + i * 0.01, 1.0) for i in range(mf.MAX_PEAKS + 1)])


def test_online_server_without_mono_answers_502_with_a_clear_message(monkeypatch):
    monkeypatch.setattr(mf, "find_mono", lambda: None)
    r = _run({"precursor": PRECURSOR, "ms2": _paste(CAFFEINE_MSFINDER_MS2)})
    assert r.status_code == 502
    detail = r.json()["detail"]
    assert "local version" in detail and mf.SETUP_COMMAND in detail


# ---------------------------------------------------------------- process watching (no MS-FINDER involved)
def _python(code):
    return [sys.executable, "-I", "-c", code]


def test_run_watched_returns_exit_code_and_output():
    rc, log, peak = mf.run_watched(_python("import sys; print('hello'); print('oops', file=sys.stderr); sys.exit(3)"), ".", dict(os.environ), timeout=30)
    assert rc == 3 and "hello" in log and "oops" in log and peak >= 0


def _alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    return True


def test_run_watched_timeout_kills_the_whole_process_group():
    code = "import subprocess, sys, time; p = subprocess.Popen(['sleep', '60']); print(p.pid, flush=True); time.sleep(60)"
    with pytest.raises(mf.MsFinderError, match="did not finish within 1 s") as exc:
        mf.run_watched(_python(code), ".", dict(os.environ), timeout=1, poll_s=0.1)
    child = int(exc.value.log.split()[0])
    for _ in range(50):
        if not _alive(child):
            break
        threading.Event().wait(0.1)
    assert not _alive(child)


def test_run_watched_memory_cap_kills_the_run():
    code = "import time; x = b'x' * (150 << 20); print('allocated', flush=True); time.sleep(60)"
    with pytest.raises(mf.MsFinderError, match="memory") as exc:
        mf.run_watched(_python(code), ".", dict(os.environ), timeout=60, rss_cap_mb=100, poll_s=0.1)
    assert isinstance(exc.value.log, str)


# ---------------------------------------------------------------- predict() plumbing with a stand-in for the process
@pytest.fixture
def fake_install(monkeypatch, tmp_path):
    """A fake install folder and a stand-in for the process that writes the given result files (so no MS-FINDER starts)."""
    exe, ini = tmp_path / "MsfinderConsoleApp.exe", tmp_path / "x.ini"
    exe.write_text("")
    ini.write_text("")
    monkeypatch.setattr(mf, "status", lambda: {"mono": "/usr/bin/mono", "dir": tmp_path, "exe": exe, "resources": tmp_path, "ini": ini})
    seen = {}

    def install(formula=CAFFEINE_FORMULA, structure=CAFFEINE_STRUCTURE, number=77, log="Writing done"):
        def fake_run(cmd, cwd, env, timeout, rss_cap_mb, **kw):
            seen.update(cmd=cmd, cwd=cwd, env=env, timeout=timeout, rss_cap_mb=rss_cap_mb)
            seen["copy_files"] = sorted(q.name for q in Path(cwd).iterdir())                 # what the run's own copy of the install holds
            seen["ini_path"] = Path(cmd[cmd.index("-m") + 1])
            seen["ini_text"] = seen["ini_path"].read_text()
            in_dir, out_dir = Path(cmd[cmd.index("-i") + 1]), Path(cmd[cmd.index("-o") + 1])
            seen["query"] = (in_dir / "query.mat").read_text()
            seen["in_dir"] = in_dir
            if formula is not None:
                (out_dir / f"Formula result-{number}.txt").write_text(formula)
            if structure is not None:
                (out_dir / f"Structure result-{number}.txt").write_text(structure)
            return 1, log, 2300.0                                   # exit code 1 = success for this program
        monkeypatch.setattr(mf, "run_watched", fake_run)
        return seen
    return install


def test_predict_runs_the_console_app_and_parses_both_files(fake_install, tmp_path):
    seen = fake_install()
    res = mf.predict(195.0877, True, None, PEAKS)
    cmd = seen["cmd"]
    assert cmd[:3] == ["/usr/bin/mono", str(Path(seen["cwd"]) / "MsfinderConsoleApp.exe"), "predict"] and cmd[-2:] == ["-m", str(tmp_path / "x.ini")]
    assert Path(seen["cwd"]) != tmp_path and tmp_path not in Path(seen["cwd"]).parents          # the run's own copy of the install, not the shared one
    assert seen["env"]["LC_ALL"] == "en_US.UTF-8" and seen["env"]["LANG"] == "en_US.UTF-8"
    assert "PRECURSORTYPE: [M+H]+" in seen["query"]                                   # Auto -> [M+H]+, written explicitly
    assert res["adduct"] == "[M+H]+" and res["formulas"][0]["formula"] == "C8H10N4O2" and res["structures"][0]["name"] == "Caffeine"
    assert (res["ms1_count"], res["ms2_count"]) == (0, 27) and res["peak_rss_mb"] == 2300.0
    assert not seen["in_dir"].exists() and not Path(seen["cwd"]).exists()             # the temp folder is gone, the copy of the install with it


def test_predict_without_the_switches_uses_the_shipped_parameter_file(fake_install, tmp_path):
    seen = fake_install()
    res = mf.predict(195.0877, True, None, PEAKS)
    assert seen["ini_path"] == tmp_path / "x.ini"                                      # the file itself, no copy
    assert res["halogen"] is False and res["silicon"] is False


@pytest.mark.parametrize("halogen,silicon", [(True, False), (False, True), (True, True)])
def test_predict_with_a_switch_runs_on_a_per_run_copy_and_leaves_the_shipped_file_alone(fake_install, tmp_path, halogen, silicon):
    shipped = tmp_path / "x.ini"
    shipped.write_text(SHIPPED_INI)
    seen = fake_install()
    res = mf.predict(195.0877, True, None, PEAKS, halogen=halogen, silicon=silicon)
    assert seen["ini_path"] != shipped and seen["ini_path"].name == mf.INI_PATH.name           # a copy with the shipped name, in the run's temp folder
    assert tmp_path not in seen["ini_path"].parents                                            # (not in the folder of the shipped file)
    assert not seen["ini_path"].exists()                                                       # gone with the temp folder
    values = _ini_values(seen["ini_text"])
    assert {k for k in ELEMENT_KEYS if values[k] == "True"} == ({"Fcheck", "ClCheck", "BrCheck", "Icheck"} if halogen else set()) | ({"SiCheck"} if silicon else set())
    assert _other_lines(seen["ini_text"], ELEMENT_KEYS) == _other_lines(SHIPPED_INI, ELEMENT_KEYS)
    assert shipped.read_text() == SHIPPED_INI and _ini_values(shipped.read_text())["Fcheck"] == "False"      # the shipped file is untouched
    assert (res["halogen"], res["silicon"]) == (halogen, silicon)


def test_real_shipped_file_is_not_modified_by_a_run_with_both_switches(fake_install, monkeypatch):
    before = mf.INI_PATH.read_bytes()
    seen = fake_install()
    monkeypatch.setattr(mf, "status", lambda: {"mono": "/usr/bin/mono", "dir": mf.INI_PATH.parent, "exe": mf.INI_PATH, "resources": mf.INI_PATH.parent,
                                                "ini": mf.INI_PATH})
    mf.predict(195.0877, True, None, PEAKS, halogen=True, silicon=True)
    assert seen["ini_path"] != mf.INI_PATH and _ini_values(seen["ini_text"])["SiCheck"] == "True"
    assert mf.INI_PATH.read_bytes() == before


def test_predict_finds_result_files_whatever_their_number(fake_install):
    fake_install(number=2123)
    assert mf.predict(195.0877, True, None, PEAKS)["structures"]
    fake_install(number=5)
    assert mf.predict(195.0877, True, None, PEAKS)["structures"]


def test_predict_without_result_files_is_an_engine_error_with_the_log(fake_install):
    fake_install(formula=None, structure=None, log="The program cannot find any query file")
    with pytest.raises(mf.MsFinderError, match="without writing a result file") as exc:
        mf.predict(195.0877, True, None, PEAKS)
    assert "cannot find any query file" in exc.value.log


def test_predict_header_only_results_are_an_empty_answer_not_an_error(fake_install):
    fake_install(formula=HEADER_ONLY, structure=STRUCTURE_HEADER_ONLY)
    res = mf.predict(50.0, True, None, PEAKS, ms1=[(50.0, 1.0)])
    assert res["formulas"] == [] and res["structures"] == [] and (res["ms1_count"], res["ms2_count"]) == (1, 4)


def test_predict_without_structure_file_still_returns_the_formulas(fake_install):
    fake_install(structure=None)
    res = mf.predict(195.0877, True, None, PEAKS)
    assert res["formulas"] and res["structures"] == []


def test_second_request_is_turned_away_while_a_run_is_in_progress(fake_install, monkeypatch):
    fake_install()
    monkeypatch.setattr(mf, "BUSY_WAIT_S", 0.05)
    with mf._RUN_LOCK:
        with pytest.raises(mf.MsFinderError, match="busy"):
            mf.predict(195.0877, True, None, PEAKS)
    assert mf.predict(195.0877, True, None, PEAKS)["structures"]                  # free again afterwards


def test_lock_is_released_after_a_failed_run(fake_install):
    fake_install(formula=None, structure=None)
    with pytest.raises(mf.MsFinderError):
        mf.predict(195.0877, True, None, PEAKS)
    assert mf._RUN_LOCK.acquire(blocking=False)
    mf._RUN_LOCK.release()


# ---------------------------------------------------------------- every run works on its own copy of the install
# MS-FINDER can rewrite Resources/*_bin.* while it runs (concurrent runs once left a 0-byte formula database), so no run uses the shared install.
def _fake_install_tree(root):
    """A small stand-in for the install folder: exe, a DLL and Resources/ with a database file."""
    (root / "Resources").mkdir(parents=True)
    (root / "MsfinderConsoleApp.exe").write_bytes(b"exe")
    (root / "Common.dll").write_bytes(b"dll")
    (root / "Resources" / "MsfinderFormulaDB-VS13_bin.efd").write_bytes(b"database")
    return root / "MsfinderConsoleApp.exe"


@pytest.fixture
def tree_install(monkeypatch, tmp_path):
    """A fake install tree; run_watched records where it ran and, like MS-FINDER's SaveToFile, destroys the database it finds there."""
    live = tmp_path / "live"
    exe = _fake_install_tree(live)
    ini = tmp_path / "x.ini"
    ini.write_text("")
    monkeypatch.setattr(mf, "status", lambda: {"mono": "/usr/bin/mono", "dir": live, "exe": exe, "resources": live / "Resources", "ini": ini})
    seen = {"cwds": []}

    def fake_run(cmd, cwd, env, timeout, rss_cap_mb, **kw):
        seen["cwds"].append(Path(cwd))
        seen["exe"] = Path(cmd[1])
        seen["files"] = sorted(str(q.relative_to(cwd)) for q in Path(cwd).rglob("*") if q.is_file())
        (Path(cwd) / "Resources" / "MsfinderFormulaDB-VS13_bin.efd").write_bytes(b"")           # the damage seen on 2026-10-30
        out_dir = Path(cmd[cmd.index("-o") + 1])
        (out_dir / "Formula result-77.txt").write_text(CAFFEINE_FORMULA)
        (out_dir / "Structure result-77.txt").write_text(CAFFEINE_STRUCTURE)
        return 1, "Writing done", 2300.0
    monkeypatch.setattr(mf, "run_watched", fake_run)
    return live, seen


def test_a_run_works_in_its_own_copy_of_the_install_and_the_shared_install_stays_intact(tree_install):
    live, seen = tree_install
    res = mf.predict(195.0877, True, None, PEAKS)
    cwd = seen["cwds"][0]
    assert cwd != live and live not in cwd.parents and seen["exe"] == cwd / "MsfinderConsoleApp.exe"
    assert seen["files"] == ["Common.dll", "MsfinderConsoleApp.exe", "Resources/MsfinderFormulaDB-VS13_bin.efd"]       # the whole install is there
    assert (live / "Resources" / "MsfinderFormulaDB-VS13_bin.efd").read_bytes() == b"database"                       # the run's damage stayed in its copy
    assert not cwd.exists()                                                                                          # deleted with the temp folder
    assert res["formulas"] and res["install_copy"] in ("clone", "copy")


def test_every_engine_entry_point_uses_a_private_copy_and_removes_it(tree_install, monkeypatch):
    live, seen = tree_install
    monkeypatch.setattr(mf, "parse_sfd", lambda text: [])
    mf.predict(195.0877, True, None, PEAKS)
    mf.explain_candidates(195.0877, True, None, PEAKS, [{"id": "a", "smiles": "C", "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C", "formula": "CH4", "mass": 16.0313}])
    mf.predict_queries([{"key": "q", "precursor_mz": 195.0877, "positive": True, "adduct": "[M+H]+", "ms2": PEAKS}])
    assert len(seen["cwds"]) == 3 and len(set(seen["cwds"])) == 3                                                    # a new copy each run
    assert all(c != live and live not in c.parents and not c.exists() for c in seen["cwds"])
    assert (live / "Resources" / "MsfinderFormulaDB-VS13_bin.efd").read_bytes() == b"database"


def test_the_copy_is_removed_when_the_run_fails(tree_install, monkeypatch):
    live, seen = tree_install
    monkeypatch.setattr(mf, "run_watched", lambda cmd, cwd, *a, **k: (seen["cwds"].append(Path(cwd)), (_ for _ in ()).throw(mf.MsFinderError("stopped")))[1])
    with pytest.raises(mf.MsFinderError, match="stopped"):
        mf.predict(195.0877, True, None, PEAKS)
    assert len(seen["cwds"]) == 1 and not seen["cwds"][0].exists() and (live / "MsfinderConsoleApp.exe").is_file()


def test_failed_clone_falls_back_to_a_plain_copy_and_says_so(tree_install, monkeypatch):
    live, seen = tree_install

    def no_clone(src, dst):
        (Path(dst) / "Resources").mkdir(parents=True)                       # a clone that failed half way
        raise OSError(18, "Invalid cross-device link")
    monkeypatch.setattr(mf, "_clone_tree", no_clone)
    res = mf.predict(195.0877, True, None, PEAKS)
    assert res["install_copy"] == "copy" and "plain copy" in res["log"] and "cross-device" in res["log"] and res["log"].endswith("Writing done")
    assert seen["cwds"][0] != live and not seen["cwds"][0].exists()
    assert (live / "Resources" / "MsfinderFormulaDB-VS13_bin.efd").read_bytes() == b"database"


def test_without_any_copy_the_run_does_not_start_and_the_shared_install_is_not_used(tree_install, monkeypatch):
    live, seen = tree_install
    monkeypatch.setattr(mf, "_clone_tree", lambda src, dst: (_ for _ in ()).throw(OSError("no clones here")))
    monkeypatch.setattr(mf.shutil, "copytree", lambda *a, **k: (_ for _ in ()).throw(OSError(28, "No space left on device")))
    with pytest.raises(mf.MsFinderError, match="private copy"):
        mf.predict(195.0877, True, None, PEAKS)
    assert seen["cwds"] == []                                                # MS-FINDER was never started
    assert mf._RUN_LOCK.acquire(blocking=False)                              # and the lock is free again
    mf._RUN_LOCK.release()


def test_private_install_is_independent_of_the_original(tmp_path):
    exe = _fake_install_tree(tmp_path / "live")
    (tmp_path / "run").mkdir()
    copy_exe, note = mf.private_install(exe, tmp_path / "run" / "install")
    copy = copy_exe.parent
    assert copy_exe == copy / "MsfinderConsoleApp.exe" and copy_exe.read_bytes() == b"exe" and (copy / "Resources" / "MsfinderFormulaDB-VS13_bin.efd").is_file()
    (copy / "Resources" / "MsfinderFormulaDB-VS13_bin.efd").write_bytes(b"")             # write into the copy (as MS-FINDER does) ...
    (copy / "Common.dll").unlink()
    assert (tmp_path / "live" / "Resources" / "MsfinderFormulaDB-VS13_bin.efd").read_bytes() == b"database"   # ... the original is not touched
    assert (tmp_path / "live" / "Common.dll").is_file()
    assert note == "" or "plain copy" in note                                                                # a clone has no note, a plain copy says why


@pytest.mark.skipif(NO_MSFINDER is not None, reason=f"MS-FINDER not available: {NO_MSFINDER}")
def test_real_run_leaves_the_shared_install_untouched_and_uses_a_private_copy(monkeypatch):
    """One real MS-FINDER run (~5 s, ~2.3 GB peak): every file of the install keeps its size and modification time, the copy is gone afterwards."""
    live = mf.status()["exe"].parent

    def snapshot():
        return {str(q.relative_to(live)): (q.stat().st_size, q.stat().st_mtime_ns) for q in live.rglob("*") if q.is_file()}
    before, cwds, real_run = snapshot(), [], mf.run_watched
    monkeypatch.setattr(mf, "run_watched", lambda cmd, cwd, *a, **k: (cwds.append(Path(cwd)), real_run(cmd, cwd, *a, **k))[1])
    ms2 = [tuple(map(float, ln.split())) for ln in CAFFEINE_MSFINDER_MS2.splitlines()]
    res = mf.predict(195.0877, True, None, ms2)
    assert res["formulas"][0]["formula"] == "C8H10N4O2" and res["install_copy"] in ("clone", "copy")
    assert snapshot() == before
    assert len(cwds) == 1 and cwds[0] != live and live not in cwds[0].parents and not cwds[0].exists()      # it ran in a copy, which is gone now


# ---------------------------------------------------------------- runner output (engine replaced)
def test_runner_fills_the_three_blocks(fake_predict):
    r = _run({"precursor": PRECURSOR, "ms2": _paste(CAFFEINE_MSFINDER_MS2)})
    assert r.status_code == 200
    body = r.json()
    out = body["outputs"]
    spec_cols = {b["id"]: b["columns"] for m in MODULES if m["id"] == "ms2_structure" for b in m["outputs"] if b["type"] == "table"}
    assert spec_cols["structures"] == ["Rank", "Structure", "Name", "SMILES", "Formula", "Score", "Ontology", "InChIKey", "Databases"]
    assert spec_cols["formulas"] == ["Rank", "Formula", "Neutral mass", "Mass error (mDa)", "Score"]
    rows = out["structures"]["rows"]
    assert out["structures"]["count"] == 15 and len(rows) == 15 and all(len(r) == 9 for r in rows)
    first = rows[0]
    assert first[0] == "1" and first[1]["img"].startswith("data:image/png;base64,") and first[2:6] == ["Caffeine", "O=C1C2=C(N=CN2C)N(C(=O)N1C)C", "C8H10N4O2", "8.1594"]
    assert first[6] == "Xanthines" and first[7] == "RYYVLZVUVIJVGH-UHFFFAOYSA-N"
    assert first[8]["text"] == "HMDB, KNApSAcK, ChEBI, DrugBank +16" and "HMDB: HMDB0001847" in first[8]["help"]
    assert rows[5][6] == "—"                                                                # empty ontology
    assert out["formulas"]["rows"][:2] == [["1", "C8H10N4O2", "194.0804", "+0.05", "4.407"], ["2", "C8H18OS2", "194.0799", "+0.52", "3.102"]]
    assert out["summary"]["values"]["Top candidate"] == {"value": "Caffeine", "note": "C8H10N4O2 · score 8.16"}
    assert out["summary"]["values"]["Candidates"]["value"] == "15 / 2"
    assert body["verdict"] == "Top candidate: <b>Caffeine</b> (C8H10N4O2), score 8.16 · 15 candidates."
    assert body["warnings"] == []


def test_runner_reports_what_auto_became(fake_predict):
    auto = _run({"precursor": PRECURSOR, "ms2": _paste(CAFFEINE_MSFINDER_MS2)}).json()["outputs"]["summary"]["values"]["Adduct used"]
    assert auto == {"value": "[M+H]+", "note": "Auto: positive mode"}
    neg = _run({"precursor": {**PRECURSOR, "ion_mode": "Negative"}, "ms2": _paste(CAFFEINE_MSFINDER_MS2)}).json()["outputs"]["summary"]["values"]["Adduct used"]
    assert neg == {"value": "[M-H]-", "note": "Auto: negative mode"}
    chosen = _run({"precursor": {**PRECURSOR, "adduct": "[M+Na]+"}, "ms2": _paste(CAFFEINE_MSFINDER_MS2)}).json()["outputs"]["summary"]["values"]["Adduct used"]
    assert chosen == "[M+Na]+"
    assert [c["adduct"] for c in fake_predict] == [None, None, "[M+Na]+"] and [c["positive"] for c in fake_predict] == [True, False, True]


def test_page_has_the_element_boxes_and_extended_range_off_by_default():
    fields = {f["id"]: f for m in MODULES if m["id"] == "ms2_structure" for f in m["inputs"]}
    for fid, label in (("halogen", "Include halogens (F, Cl, Br, I)"), ("silicon", "Include Si"),
                       ("extended", "Extended element range")):
        assert fields[fid]["type"] == "checkbox" and fields[fid]["label"] == label and fields[fid]["value"] is False
    assert not any(k in m.get("example", {}) for m in MODULES if m["id"] == "ms2_structure" for k in ("halogen", "silicon", "extended"))


def test_runner_passes_the_element_boxes_to_the_engine_independently(fake_predict):
    base = {"precursor": PRECURSOR, "ms2": _paste(CAFFEINE_MSFINDER_MS2)}
    _run(base)
    _run({**base, "halogen": True})
    _run({**base, "silicon": True})
    _run({**base, "halogen": True, "silicon": True})
    _run({**base, "halogen": False, "silicon": None})
    assert [(c["halogen"], c["silicon"]) for c in fake_predict] == [(False, False), (True, False), (False, True), (True, True), (False, False)]


@pytest.mark.parametrize("field,value", [("halogen", "yes"), ("halogen", "false"), ("halogen", 1), ("silicon", "true"), ("silicon", [True])])
def test_runner_rejects_a_checkbox_that_is_not_true_or_false(fake_predict, field, value):
    r = _run({"precursor": PRECURSOR, "ms2": _paste(CAFFEINE_MSFINDER_MS2), field: value})
    assert r.status_code == 400 and "must be true or false" in r.json()["detail"]
    assert fake_predict == []                                                          # MS-FINDER is never started for it


def test_runner_states_which_elements_were_searched(fake_predict):
    def card(**boxes):
        return _run({"precursor": PRECURSOR, "ms2": _paste(CAFFEINE_MSFINDER_MS2), **boxes}).json()["outputs"]["summary"]["values"]["Elements searched"]
    assert card() == {"value": "C H N O P S", "note": "halogens off · Si off · common range"}
    assert card(halogen=True) == {"value": "C H N O P S F Cl Br I", "note": "halogens on · Si off · common range"}
    assert card(silicon=True) == {"value": "C H N O P S Si", "note": "halogens off · Si on · common range"}
    assert card(halogen=True, silicon=True) == {"value": "C H N O P S F Cl Br I Si", "note": "halogens on · Si on · common range"}
    assert card(extended=True) == {"value": "C H N O P S", "note": "halogens off · Si off · extended range"}
    labels = [c["label"] for m in MODULES if m["id"] == "ms2_structure" for b in m["outputs"] if b.get("id") == "summary" for c in b["items"]]
    assert "Elements searched" in labels


def test_runner_passes_the_peaks_and_the_ms1_pattern_to_the_engine(fake_predict):
    _run({"precursor": PRECURSOR, "ms2": _paste("100 5\n200 50"), "ms1": CAFFEINE_MS1})
    call = fake_predict[0]
    assert call["precursor_mz"] == 195.0877 and call["ms2"] == [(100.0, 5.0), (200.0, 50.0)]
    assert call["ms1"] == [(195.0877, 100.0), (196.0902, 10.3), (197.0922, 0.9)]


def test_runner_takes_the_precursor_from_the_spectrum_header(fake_predict):
    mgf = "BEGIN IONS\nPEPMASS=195.0877\n138.0662 100\n195.0877 20\nEND IONS"
    body = _run({"precursor": {**PRECURSOR, "precursor_mz": ""}, "ms2": _paste(mgf)}).json()
    assert fake_predict[0]["precursor_mz"] == 195.0877
    assert any("read from the spectrum" in w for w in body["warnings"])


def test_runner_warns_when_ms1_does_not_line_up_with_the_precursor(fake_predict):
    body = _run({"precursor": PRECURSOR, "ms2": _paste(CAFFEINE_MSFINDER_MS2), "ms1": "195.0907 100\n196.0932 10"}).json()
    assert any("0.5 mDa" in w and "isotope" in w for w in body["warnings"])
    ok = _run({"precursor": PRECURSOR, "ms2": _paste(CAFFEINE_MSFINDER_MS2), "ms1": CAFFEINE_MS1}).json()
    assert not any("mDa" in w for w in ok["warnings"])
    junk = _run({"precursor": PRECURSOR, "ms2": _paste(CAFFEINE_MSFINDER_MS2), "ms1": "no peaks here"}).json()
    assert any("No valid MS1 peaks" in w for w in junk["warnings"])


def test_runner_caps_the_tables(monkeypatch, fake_predict):
    res = _fake_result()
    res["structures"] = res["structures"] * 4                           # 60 rows
    res["formulas"] = [dict(res["formulas"][0], formula=f"C{i}") for i in range(1, 61)]
    monkeypatch.setattr(mf, "predict", lambda *a, **k: res)
    body = _run({"precursor": PRECURSOR, "ms2": _paste(CAFFEINE_MSFINDER_MS2)}).json()
    assert len(body["outputs"]["structures"]["rows"]) == 50 and body["outputs"]["structures"]["count"] == 60
    assert len(body["outputs"]["formulas"]["rows"]) == 50
    assert any("best 50" in w for w in body["warnings"])
    assert body["verdict"].endswith("· 60 candidates.")


def test_runner_without_formula_or_structure_says_so(monkeypatch):
    monkeypatch.setattr(mf, "predict", lambda *a, **k: _fake_result(formulas=[], structures=[]))
    body = _run({"precursor": PRECURSOR, "ms2": _paste(CAFFEINE_MSFINDER_MS2)}).json()
    assert body["status"] == "ok" and body["outputs"]["structures"]["rows"] == [] and body["outputs"]["summary"]["values"]["Top candidate"] == "—"
    assert any("no molecular formula" in w for w in body["warnings"])
    assert body["verdict"] == "MS-FINDER found no molecular formula for this precursor."

    only_formulas = _fake_result(structures=[])
    monkeypatch.setattr(mf, "predict", lambda *a, **k: only_formulas)
    body = _run({"precursor": PRECURSOR, "ms2": _paste(CAFFEINE_MSFINDER_MS2)}).json()
    assert any("local structure databases have no candidate" in w for w in body["warnings"])
    assert "No structure in MS-FINDER's local database for <b>C8H10N4O2</b>" in body["verdict"]


def test_runner_hint_names_the_element_box_that_is_off(monkeypatch):
    monkeypatch.setattr(mf, "predict", lambda *a, **k: _fake_result(formulas=[], structures=[]))
    base = {"precursor": PRECURSOR, "ms2": _paste(CAFFEINE_MSFINDER_MS2)}

    def warning(**boxes):
        return next(w for w in _run({**base, **boxes}).json()["warnings"] if "no molecular formula" in w)
    w = warning()
    assert "elements C H N O P S)" in w and "Off: Include halogens (F, Cl, Br, I), Include Si, Extended element range" in w
    w = warning(halogen=True)
    assert "elements C H N O P S F Cl Br I)" in w and "Off: Include Si, Extended element range" in w and "Include halogens" not in w
    w = warning(silicon=True)
    assert "elements C H N O P S Si)" in w and "Off: Include halogens (F, Cl, Br, I), Extended element range" in w
    w = warning(halogen=True, silicon=True)
    assert "elements C H N O P S F Cl Br I Si)" in w and "Off: Extended element range" in w
    w = warning(halogen=True, silicon=True, extended=True)
    assert "tick" not in w                  # every relevant switch is on: nothing left to suggest

    only_formulas = _fake_result(structures=[])
    monkeypatch.setattr(mf, "predict", lambda *a, **k: only_formulas)
    body = _run(base).json()
    assert any("no candidate for C8H10N4O2" in w and "tick a box" in w for w in body["warnings"])
    body = _run({**base, "halogen": True, "silicon": True, "extended": True}).json()
    assert any("no candidate for C8H10N4O2" in w and "tick" not in w for w in body["warnings"])


def test_engine_failures_are_502_with_the_log_and_bad_input_is_400(monkeypatch):
    def broken(*a, **k):
        raise mf.MsFinderError("MS-FINDER did not finish within 300 s and was stopped.", "last lines")
    monkeypatch.setattr(mf, "predict", broken)
    r = _run({"precursor": PRECURSOR, "ms2": _paste(CAFFEINE_MSFINDER_MS2)})
    assert r.status_code == 502 and r.json() == {"detail": "MS-FINDER did not finish within 300 s and was stopped.", "log": "last lines"}

    def bad(*a, **k):
        raise ValueError("Precursor m/z must be positive.")
    monkeypatch.setattr(mf, "predict", bad)
    assert _run({"precursor": PRECURSOR, "ms2": _paste(CAFFEINE_MSFINDER_MS2)}).status_code == 400
    assert _run({"precursor": PRECURSOR, "ms2": _paste("not a spectrum")}).status_code == 400


# ---------------------------------------------------------------- setup command
def test_setup_check_reports_what_is_missing(monkeypatch, tmp_path, capsys):
    monkeypatch.setenv("MZDESK_MSFINDER_DIR", str(tmp_path / "nothing"))
    assert su.check() == 1
    out = capsys.readouterr().out
    assert "NOT FOUND" in out and mf.SETUP_COMMAND in out and str(tmp_path / "nothing") in out


@pytest.mark.skipif(NO_MSFINDER is not None, reason="MS-FINDER not installed")
def test_setup_check_is_ready_on_an_installed_machine(capsys):
    assert su.main(["--check"]) == 0
    assert "Ready            : yes" in capsys.readouterr().out


def test_setup_check_never_starts_the_download_or_msfinder(monkeypatch, tmp_path):
    monkeypatch.setenv("MZDESK_MSFINDER_DIR", str(tmp_path / "nothing"))
    monkeypatch.setattr(su, "_download", lambda *a, **k: pytest.fail("--check must not download"))
    monkeypatch.setattr(mf, "run_watched", lambda *a, **k: pytest.fail("--check must not start MS-FINDER"))
    assert su.main(["--check"]) == 1


def test_setup_pins_nine_sources_with_sha256():
    assert len(su.SOURCES) == 9 and all(re.fullmatch(r"[0-9a-f]{64}", h) for h in su.SOURCES.values())
    assert "Program.cs" in su.SOURCES and "Process/PredictProcess.cs" in su.SOURCES
    assert re.fullmatch(r"[0-9a-f]{64}", su.ZIP_SHA256) and su.ZIP_URL.endswith("/MSFINDER-v3.73/MSFINDER.v3.73.zip")


def test_setup_refuses_a_source_that_does_not_match_its_pin(tmp_path):
    (tmp_path / "Process").mkdir()
    good, other = tmp_path / "a.cs", tmp_path / "Process" / "b.cs"
    good.write_text("class A {}")
    other.write_text("class B {}")
    table = {"a.cs": su.sha256_of(good), "Process/b.cs": su.sha256_of(other)}
    su.verify_sources(tmp_path, table)
    other.write_text("class B { /* patched */ }")
    with pytest.raises(su.SetupError, match="Process/b.cs"):
        su.verify_sources(tmp_path, table)
    other.unlink()
    with pytest.raises(su.SetupError, match="Missing"):
        su.verify_sources(tmp_path, table)


def test_setup_compile_command_matches_the_verified_build(tmp_path):
    src = tmp_path / "src"
    (src / "Process").mkdir(parents=True)
    (src / "Test").mkdir()
    for name in ("Program.cs", "Process/B.cs", "Process/A.cs", "Test/T.cs"):
        (src / name).write_text("")
    cmd = su.compile_command("/m/csc", "/m/Facades", "/bin/dir", src, "/out/x.exe")
    assert cmd[:5] == ["/m/csc", "-nologo", "-target:exe", "-platform:anycpu", "-out:/out/x.exe"]
    refs = [a for a in cmd if a.startswith("-r:")]
    assert "-r:/bin/dir/MsfinderCommon.dll" in refs and "-r:/m/Facades/netstandard.dll" in refs and "-r:System.Numerics.dll" in refs and len(refs) == 11 + 4
    assert [Path(a).name for a in cmd if a.endswith(".cs")] == ["Program.cs", "A.cs", "B.cs", "T.cs"]


def test_setup_finds_mono_facades_from_the_mono_layout(tmp_path):
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "mono").write_text("")
    (tmp_path / "lib" / "mono" / "4.5" / "Facades").mkdir(parents=True)
    assert su.find_facades(tmp_path / "bin" / "mono") == tmp_path.resolve() / "lib" / "mono" / "4.5" / "Facades"


def test_setup_unzip_drops_the_top_folder_and_blocks_path_traversal(tmp_path, monkeypatch):
    z = tmp_path / "r.zip"
    with zipfile.ZipFile(z, "w") as zf:
        zf.writestr("MSFINDER v3.73/", "")
        zf.writestr("MSFINDER v3.73/Common.dll", "dll")
        zf.writestr("MSFINDER v3.73/Resources/AdductPositives.apf", "[M+H]+")
    su.unzip(z, tmp_path / "install")
    assert (tmp_path / "install" / "Common.dll").read_text() == "dll" and (tmp_path / "install" / "Resources" / "AdductPositives.apf").is_file()
    evil = tmp_path / "evil.zip"
    with zipfile.ZipFile(evil, "w") as zf:
        zf.writestr("MSFINDER v3.73/../../escaped.txt", "x")
    with pytest.raises(su.SetupError, match="Unsafe"):
        su.unzip(evil, tmp_path / "install2")
    assert not (tmp_path / "escaped.txt").exists()


# ---------------------------------------------------------------- the real thing
@pytest.mark.skipif(NO_MSFINDER is not None, reason=f"MS-FINDER not available: {NO_MSFINDER}")
def test_real_run_ranks_caffeine_first_with_ms1_isotopes():
    """One real MS-FINDER run: MassBank MSBNK-LCSB-LU119905 (CC BY 4.0) + an MS1 isotope pattern. ~5 s, ~2.3 GB peak."""
    ms2 = [tuple(map(float, ln.split())) for ln in CAFFEINE_MSFINDER_MS2.splitlines()]
    ms1 = [tuple(map(float, ln.split())) for ln in CAFFEINE_MS1.splitlines()]
    res = mf.predict(195.0877, True, None, ms2, ms1)
    assert res["adduct"] == "[M+H]+"
    assert res["formulas"][0]["formula"] == "C8H10N4O2" and abs(res["formulas"][0]["mass_error_da"]) < 0.001
    assert res["structures"][0]["inchikey"].startswith("RYYVLZVUVIJVGH") and res["structures"][0]["name"] == "Caffeine"
    assert (res["ms1_count"], res["ms2_count"]) == (3, 27)                       # the isotope peaks were read
    assert res["peak_rss_mb"] < mf.RSS_CAP_MB and res["seconds"] < mf.TIMEOUT_S


# Halogen / Si spectra for the two real runs below: real MassBank records (peak lists as the MassBank API gives them, "rel" column), nothing invented.
# Diclofenac [M+H]+: MSBNK-CASMI_2016-SM856902, licence CC BY (MassBank), LC-ESI-QFT, MS2, CE 35 NCE, authors Krauss M, Schymanski EL (2016).
DICLOFENAC_PRECURSOR = 296.024
DICLOFENAC_MS2 = [(179.0729, 2), (180.0806, 13), (214.0417, 796), (215.0495, 999), (243.0402, 6), (250.0185, 530), (278.0134, 287), (296.0240, 164)]
# 3,3'-(1,1,3,3-Tetramethyldisiloxane-1,3-diyl)dipropan-1-amine [M+H]+ (two Si, no halogen): MSBNK-EPA-ENTACT_AGILENT000481, licence CC BY (MassBank),
# ESI-QTOF, MS2, CE 40, EPA CCTE and Agilent Technologies (2023).
SILOXANE_PRECURSOR = 249.1813
SILOXANE_MS2 = [(59.0312, 44), (74.0421, 999), (74.0726, 19), (74.0838, 12), (75.0261, 43), (116.0890, 263)]


def _formulas(res):
    return [f["formula"] for f in res["formulas"]]


@pytest.mark.skipif(NO_MSFINDER is not None, reason=f"MS-FINDER not available: {NO_MSFINDER}")
def test_real_run_diclofenac_is_found_only_with_halogens_on():
    """Two real MS-FINDER runs (~5 s, ~2.3 GB peak each, one after the other): the parameter file keeps Cl off, the per-run copy switches it on."""
    off = mf.predict(DICLOFENAC_PRECURSOR, True, None, DICLOFENAC_MS2)
    assert "C14H11Cl2NO2" not in _formulas(off) and off["structures"] == [] and off["halogen"] is False
    on = mf.predict(DICLOFENAC_PRECURSOR, True, None, DICLOFENAC_MS2, halogen=True)
    assert on["halogen"] is True and on["formulas"][0]["formula"] == "C14H11Cl2NO2" and abs(on["formulas"][0]["mass_error_da"]) < 0.001
    assert on["structures"][0]["name"] == "Diclofenac" and on["structures"][0]["inchikey"].startswith("DCOPUUMXTXDBNB")
    assert on["peak_rss_mb"] < mf.RSS_CAP_MB and on["seconds"] < mf.TIMEOUT_S


@pytest.mark.skipif(NO_MSFINDER is not None, reason=f"MS-FINDER not available: {NO_MSFINDER}")
def test_real_run_siloxane_needs_the_si_box_and_not_the_halogen_box():
    """Two real runs: halogens on alone does not allow Si; Si on does (rank 1 formula and structure)."""
    hal = mf.predict(SILOXANE_PRECURSOR, True, None, SILOXANE_MS2, halogen=True)
    assert "C10H28N2OSi2" not in _formulas(hal) and hal["silicon"] is False
    si = mf.predict(SILOXANE_PRECURSOR, True, None, SILOXANE_MS2, silicon=True)
    assert si["silicon"] is True and si["formulas"][0]["formula"] == "C10H28N2OSi2" and abs(si["formulas"][0]["mass_error_da"]) < 0.001
    assert si["structures"][0]["formula"] == "C10H28N2OSi2" and si["structures"][0]["inchikey"].startswith("GPXCORHXFPYJEH")
