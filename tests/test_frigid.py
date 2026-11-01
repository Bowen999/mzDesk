"""De novo structure（FRIGID-base，本机 PyTorch）的测试：加合物映射、任务文件、候选过滤、BUDDY 回退、退出码映射、重任务锁、进程看护、runner 外形、安装脚本。

fixtures：tests/fixtures/frigid/caffeine_quick.result.json 是驱动脚本真实跑出来的结果（咖啡因，MassBank MSBNK-LCSB-LU119905 的谱，CC BY 4.0，Quick 模式，MPS），路径换成了中性路径。
真实运行（test_real_run_*）每次约 20–40 s、峰值约 4–5 GB 内存，所以默认跳过：要设 MZDESK_RUN_FRIGID=1 并且已安装（python -m mzdesk.ms2.frigid_setup）；
内存紧张（空闲 < 20 % 或 swap 余量 < 500 MB）时它们自己也会跳过。其余测试都不启动 FRIGID。
"""
import ast
import hashlib
import importlib.util
import io
import json
import os
import re
import subprocess
import sys
import tarfile
import threading
import time
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from mzdesk.api.app import app
from mzdesk.ms2 import frigid_engine as fe
from mzdesk.ms2 import frigid_setup as su
from mzdesk.ms2 import msfinder_engine as mf
from mzdesk.shared import heavy_lock
from mzdesk.spec import CAFFEINE_MSFINDER_MS2, FRIGID_ADDUCTS, MODULES

ROOT = Path(__file__).resolve().parent.parent
FIX = Path(__file__).parent / "fixtures" / "frigid"
DRIVER_RESULT = json.loads((FIX / "caffeine_quick.result.json").read_text())
client = TestClient(app)
PEAKS = [(56.0495, 88.0), (110.0712, 590.0), (138.0662, 999.0), (195.0876, 58.0)]
PRECURSOR = {"precursor_mz": "195.0877", "ion_mode": "Positive", "adduct": "Auto"}
SMILES_CAFFEINE = "Cn1c(=O)c2c(ncn2C)n(C)c1=O"


def _paste(peaks):
    return {"source": "Paste", "value": peaks}


def _run(inputs):
    return client.post("/api/modules/ms2_denovo/run", json=inputs)


def _inputs(**over):
    return {"precursor": PRECURSOR, "ms2": _paste(CAFFEINE_MSFINDER_MS2), "formula": "C8H10N4O2", **over}


def _unavailable():
    try:
        fe.check_frigid()
    except fe.FrigidError as exc:
        return str(exc)
    return None


NO_FRIGID = _unavailable()


def _module():
    return next(m for m in MODULES if m["id"] == "ms2_denovo")


def _fake_result(**over):
    """What fe.predict returns for the recorded caffeine run (candidates filtered like the engine does)."""
    kept, dropped = fe.filter_candidates(DRIVER_RESULT["candidates"], "[M+H]+")
    res = {"formula": "C8H10N4O2", "adduct": "[M+H]+", "driver_adduct": "[M+H]+", "mode": "quick", "seed": 42, "candidates": kept, "dropped": dropped,
           "n_raw": 7, "n_samples": 16, "n_valid": 16, "n_formula_matched": 9, "n_unique_matched": 4, "n_peaks_assigned": 27, "device": "mps",
           "device_requested": "auto", "device_reason": DRIVER_RESULT["device_decision"]["reason"], "seconds": 21.4, "driver_seconds": 17.3,
           "peak_footprint_mb": 5047.4, "peak_rss_mb": 1586.6, "warnings": [], "length_model": True, "log": ""}
    res.update(over)
    return res


@pytest.fixture
def fake_predict(monkeypatch):
    """Replace the engine call of the runner; returns the dict that records the arguments."""
    seen = {}

    def fake(formula, precursor_mz, adduct, peaks, mode="quick", **kw):
        seen.update(formula=formula, precursor_mz=precursor_mz, adduct=adduct, peaks=list(peaks), mode=mode, kw=kw)
        return seen["result"]
    seen["result"] = _fake_result()
    monkeypatch.setattr(fe, "predict", fake)
    return seen


@pytest.fixture
def fake_install(monkeypatch, tmp_path):
    """A fake install and a stand-in for the driver process that writes a result.json (so no FRIGID starts)."""
    py = tmp_path / "env" / "bin" / "python"
    monkeypatch.setattr(fe, "status", lambda: {"dir": tmp_path, "python": py, "clone": tmp_path / "FRIGID", "mps_patch": True, "safe_patch": True,
                                               "dlm": tmp_path / "d", "mist": tmp_path / "m", "tokenizer": tmp_path / "t", "driver": fe.DRIVER,
                                               "length_model": fe.LENGTH_MODEL})
    seen = {}

    def install(result=None, code=0, log="driver log", memory=None):
        def fake_run(cmd, cwd, env, timeout, **kw):
            seen.update(cmd=cmd, cwd=cwd, env=env, timeout=timeout)
            job_path, out_path = Path(cmd[cmd.index("--input") + 1]), Path(cmd[cmd.index("--output") + 1])
            seen["job"] = json.loads(job_path.read_text())
            seen["tmp"] = job_path.parent
            if result != "none":
                out_path.write_text(json.dumps(DRIVER_RESULT if result is None else result))
            return code, log, memory or {"peak_rss_mb": 1600.0, "peak_footprint_mb": 5100.0}
        monkeypatch.setattr(fe, "run_watched", fake_run)
        return seen
    return install


def _load_driver_module(name):
    spec = importlib.util.spec_from_file_location(f"frigid_driver_{name}", fe.DRIVER.parent / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------- adducts
def test_auto_adduct_is_mh_and_the_driver_names_are_mapped():
    assert fe.resolve_adduct(None) == ("[M+H]+", "[M+H]+") and fe.resolve_adduct("Auto") == ("[M+H]+", "[M+H]+")
    assert fe.resolve_adduct("[M+NH4]+") == ("[M+NH4]+", "[M+H3N+H]+")
    assert fe.resolve_adduct("[M+H-H2O]+") == ("[M+H-H2O]+", "[M-H2O+H]+")
    assert fe.resolve_adduct("[M+H-2H2O]+") == ("[M+H-2H2O]+", "[M-H4O2+H]+")
    assert fe.resolve_adduct("[M]+") == ("[M]+", "[M]+")
    assert fe.resolve_adduct("[M-H2O+H]+")[0] == "[M+H-H2O]+"            # the driver's own spelling is accepted too


def test_seven_positive_adducts_and_the_page_offers_auto_first():
    assert fe.ADDUCTS == ["[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M+H-H2O]+", "[M]+", "[M+H-2H2O]+"]
    field = _module()["inputs"][0]
    assert field["adducts"] == FRIGID_ADDUCTS == {"Positive": ["Auto", *fe.ADDUCTS], "Negative": ["Auto"]}


@pytest.mark.parametrize("adduct", ["[M-H]-", "[M+Cl]-", "[M-H2O-H]-"])
def test_negative_adducts_say_positive_mode_only(adduct):
    with pytest.raises(ValueError, match="positive ion mode only"):
        fe.resolve_adduct(adduct)


def test_unknown_adduct_lists_the_supported_ones():
    with pytest.raises(ValueError, match="does not support the adduct .*\\[M\\+Li\\]\\+.*Supported: \\[M\\+H\\]\\+"):
        fe.resolve_adduct("[M+Li]+")


def test_engine_adduct_table_agrees_with_the_shipped_driver():
    """Every driver name is one the driver accepts, and the added mass is the driver's own ION_TO_MASS."""
    sa = _load_driver_module("subformula_assign")
    for name, (driver_name, delta) in fe.ADDUCT_TABLE.items():
        assert sa.standardize_adduct(driver_name) == driver_name
        assert delta == pytest.approx(sa.ION_TO_MASS[driver_name], abs=1e-9), name
    assert fe.ADDUCT_TABLE["[M+H]+"][1] == pytest.approx(1.007276, abs=1e-6)
    assert fe.ion_mz(194.080376, "[M+H]+") == pytest.approx(195.0877, abs=1e-4)           # caffeine


def test_driver_elements_match_the_engine_element_list():
    sa = _load_driver_module("subformula_assign")
    assert list(fe.ELEMENTS) == sa.VALID_ELEMENTS and len(fe.ELEMENTS) == 18


# ---------------------------------------------------------------- formula
def test_formula_is_checked_and_written_in_hill_order():
    assert fe.check_formula(" H10C8N4O2 ") == "C8H10N4O2"
    assert fe.check_formula("C6H12O6") == "C6H12O6" and fe.check_formula("C2H6OS") == "C2H6OS"
    assert fe.check_formula("C16H15F2N3Si") == "C16H15F2N3Si"
    assert fe.formula_counts("C8H10N4O2") == {"C": 8, "H": 10, "N": 4, "O": 2}


@pytest.mark.parametrize("bad,match", [("", "plain element formula"), ("c8h10", "plain element formula"), ("C8H10N4O2+", "plain element formula"),
                                       ("C8H10Zn", "cannot use Zn"), ("C8H10Xx", "cannot use Xx"), ("H2O", "no carbon"), ("C0H4", "no carbon"),
                                       ("C500H1000N200O200S100", "too large")])
def test_bad_formulas_are_refused_with_a_reason(bad, match):
    with pytest.raises(ValueError, match=match):
        fe.check_formula(bad)


# ---------------------------------------------------------------- mode and device
def test_mode_comes_from_the_page_label_or_the_id():
    assert fe.resolve_mode(None) == "quick" and fe.resolve_mode("") == "quick"
    assert fe.resolve_mode("Quick (16 samples)") == "quick" and fe.resolve_mode("Full (up to 100 samples)") == "full"
    assert fe.resolve_mode("FULL") == "full"
    with pytest.raises(ValueError, match="Mode must be one of"):
        fe.resolve_mode("Slow")


def test_device_defaults_to_auto_and_can_be_overridden_by_the_environment(monkeypatch):
    monkeypatch.delenv("MZDESK_FRIGID_DEVICE", raising=False)
    assert fe.resolve_device() == "auto"
    monkeypatch.setenv("MZDESK_FRIGID_DEVICE", "CPU")
    assert fe.resolve_device() == "cpu" and fe.resolve_device("mps") == "mps"
    monkeypatch.setenv("MZDESK_FRIGID_DEVICE", "cuda")
    with pytest.raises(fe.FrigidError, match="MZDESK_FRIGID_DEVICE"):
        fe.resolve_device()


# ---------------------------------------------------------------- job file and environment
def test_job_json_has_what_the_driver_reads():
    job = fe.build_job("C8H10N4O2", 195.0877, "[M+H3N+H]+", [(138.0662, 999.0), (56.0495, 88.0), (50.0, 0.0), (-1.0, 5.0)], "full",
                       device="cpu", seed=7, name="caffeine; rm -rf /")
    assert job == {"formula": "C8H10N4O2", "precursor_mz": 195.0877, "adduct": "[M+H3N+H]+", "ion_mode": "positive",
                   "peaks": [[56.0495, 88.0], [138.0662, 999.0]], "mode": "full", "device": "cpu", "seed": 7, "name": "caffeine_ rm -rf _",
                   "length_model": str(fe.LENGTH_MODEL)}
    json.dumps(job)


def test_job_passes_the_drivers_own_validation(monkeypatch, tmp_path):
    """frigid_predict.validate() needs only numpy (via subformula_assign); importing the module sets a few defaults, which are undone here."""
    for key in ("HF_HOME", "HF_HUB_OFFLINE", "TOKENIZERS_PARALLELISM", "FRIGID_ROOT"):
        monkeypatch.setenv(key, os.environ.get(key, str(tmp_path)))
    monkeypatch.setattr(sys, "path", list(sys.path))
    monkeypatch.delitem(sys.modules, "subformula_assign", raising=False)
    drv = _load_driver_module("frigid_predict")
    for name, (driver_name, _) in fe.ADDUCT_TABLE.items():
        out = drv.validate(fe.build_job("C8H10N4O2", 195.0877, driver_name, PEAKS, "quick", device="auto"))
        assert out["adduct"] == driver_name and out["mode"] == "quick" and out["device"] == "auto" and out["length_model"] == str(fe.LENGTH_MODEL), name
    with pytest.raises(drv.InputError, match="negative ion mode"):
        drv.validate({**fe.build_job("C8H10N4O2", 195.0877, "[M+H]+", PEAKS), "ion_mode": "negative"})
    with pytest.raises(drv.InputError, match="not one of the positive adducts"):
        drv.validate(fe.build_job("C8H10N4O2", 195.0877, "[M+Li]+", PEAKS))
    monkeypatch.delitem(sys.modules, "subformula_assign", raising=False)


def test_child_environment_is_offline_isolated_and_points_into_the_install(monkeypatch, tmp_path):
    monkeypatch.setenv("PYTHONPATH", "/somewhere")
    monkeypatch.setenv("HF_HUB_CACHE", "/elsewhere")
    monkeypatch.setenv("FRIGID_DEVICE", "cuda")
    monkeypatch.setenv("PATH", "/usr/bin:/bin")
    env = fe.child_env(tmp_path)
    assert "PYTHONPATH" not in env and "HF_HUB_CACHE" not in env and "FRIGID_DEVICE" not in env
    assert env["HF_HOME"] == str(tmp_path / "hf_cache") and env["HF_HUB_OFFLINE"] == "1" and env["PYTHONNOUSERSITE"] == "1"
    assert env["FRIGID_ROOT"] == str(tmp_path / "FRIGID")
    assert env["FRIGID_DLM_CKPT"] == str(tmp_path / "checkpoints" / "DLM.ckpt") and env["FRIGID_MIST_CKPT"].endswith("mist_msg_aug.ckpt")
    assert env["PATH"] == "/usr/bin:/bin"


def test_cpu_threads_follow_the_idle_cores_and_never_drop_below_two():
    assert fe.cpu_threads(load=0.3, ncpu=10) == 10 and fe.cpu_threads(load=4.0, ncpu=10) == 6 and fe.cpu_threads(load=9.6, ncpu=10) == 2
    assert fe.cpu_threads(load=13.7, ncpu=10) == 2 and fe.cpu_threads(load=0.0, ncpu=1) == 2 and fe.cpu_threads(load=30, ncpu=64) == 34


def test_child_environment_limits_the_threads(monkeypatch, tmp_path):
    monkeypatch.delenv("OMP_NUM_THREADS", raising=False)
    monkeypatch.delenv("MZDESK_FRIGID_THREADS", raising=False)
    monkeypatch.setattr(fe, "cpu_threads", lambda: 3)
    assert fe.child_env(tmp_path)["OMP_NUM_THREADS"] == "3"
    monkeypatch.setenv("MZDESK_FRIGID_THREADS", "5")
    assert fe.child_env(tmp_path)["OMP_NUM_THREADS"] == "5"
    monkeypatch.delenv("MZDESK_FRIGID_THREADS")
    monkeypatch.setenv("OMP_NUM_THREADS", "7")                                       # the server's own setting wins over the automatic number
    assert fe.child_env(tmp_path)["OMP_NUM_THREADS"] == "7"
    monkeypatch.delenv("OMP_NUM_THREADS")
    monkeypatch.setattr(fe, "cpu_threads", lambda: None)                              # no load average on this system: PyTorch's default
    assert "OMP_NUM_THREADS" not in fe.child_env(tmp_path)


# ---------------------------------------------------------------- candidate filter
def _cand(smiles, score, **over):
    row = {"rank": 1, "smiles": smiles, "formula": "X", "score": score, "formula_match": True, "inchikey_2d": smiles[:14].ljust(14, "A").upper(),
           "isotope_labeled": False, "radical": False, "charged": False}
    row.update(over)
    return row


def test_filter_keeps_only_formula_matched_and_renumbers():
    kept, dropped = fe.filter_candidates(DRIVER_RESULT["candidates"])
    assert [c["smiles"] for c in kept] == [c["smiles"] for c in DRIVER_RESULT["candidates"][:4]]
    assert [c["rank"] for c in kept] == [1, 2, 3, 4] and dropped["not_formula_match"] == 3 and sum(dropped.values()) == 3
    assert set(kept[0]) == {"rank", "smiles", "formula", "score", "inchikey_2d"}


def test_filter_drops_padding_isotope_labelled_and_radical_rows():
    raw = [_cand("CCO", 0.9), _cand("CC", 0.8, formula_match=False), _cand("[2H]C", 0.7, isotope_labeled=True), _cand("C[O]", 0.6, radical=True),
           _cand("CCN", 0.5)]
    kept, dropped = fe.filter_candidates(raw)
    assert [c["smiles"] for c in kept] == ["CCO", "CCN"]
    assert dropped == {"not_formula_match": 1, "isotope": 1, "radical": 1, "charged": 0, "duplicate": 0}


def test_filter_judges_the_net_charge_not_the_driver_flag():
    """The driver flags any atom with a formal charge: a nitro group or an N-oxide is charged by that flag but a neutral compound."""
    nitro = "O=[N+]([O-])c1ccccc1"
    cation = "C[N+](C)(C)CCO"
    raw = [_cand(nitro, 0.9, charged=True), _cand("CC[NH3+]", 0.8, charged=True), _cand(cation, 0.7, charged=True)]
    kept, dropped = fe.filter_candidates(raw, "[M+H]+")
    assert [c["smiles"] for c in kept] == [nitro] and dropped["charged"] == 2
    kept, dropped = fe.filter_candidates(raw, "[M]+")                              # permanent cation: the formula is the cation's, net charge +1 is right
    assert [c["smiles"] for c in kept] == ["CC[NH3+]", cation] and dropped["charged"] == 1


def test_filter_keeps_one_row_per_2d_structure():
    raw = [_cand("CCO", 0.9, inchikey_2d="AAAAAAAAAAAAAA"), _cand("OCC", 0.8, inchikey_2d="AAAAAAAAAAAAAA"), _cand("CCN", 0.5, inchikey_2d="BBBBBBBBBBBBBB")]
    kept, dropped = fe.filter_candidates(raw)
    assert [c["smiles"] for c in kept] == ["CCO", "CCN"] and dropped["duplicate"] == 1


def test_filter_of_nothing_is_empty():
    assert fe.filter_candidates([]) == ([], {"not_formula_match": 0, "isotope": 0, "radical": 0, "charged": 0, "duplicate": 0})


# ---------------------------------------------------------------- installation
def test_missing_install_says_what_to_run(monkeypatch, tmp_path):
    monkeypatch.setenv("MZDESK_FRIGID_DIR", str(tmp_path / "nothing"))
    with pytest.raises(fe.FrigidError) as exc:
        fe.check_frigid()
    msg = str(exc.value)
    assert fe.SETUP_COMMAND in msg and str(tmp_path / "nothing") in msg and "MZDESK_FRIGID_DIR" in msg and "local version" in msg
    assert "env/" in msg and "checkpoints/DLM.ckpt" in msg


def test_partial_checkpoint_counts_as_missing(monkeypatch, tmp_path):
    (tmp_path / "checkpoints").mkdir()
    (tmp_path / "checkpoints" / "DLM.ckpt").write_bytes(b"x" * 100)
    monkeypatch.setenv("MZDESK_FRIGID_DIR", str(tmp_path))
    assert fe.status()["dlm"] is None


def test_input_errors_come_before_the_install_check(monkeypatch, tmp_path):
    monkeypatch.setenv("MZDESK_FRIGID_DIR", str(tmp_path / "nothing"))
    for args, match in (((("C8H10Zn", 195.0877, None, PEAKS)), "cannot use Zn"), ((("C8H10N4O2", 0, None, PEAKS)), "Precursor"),
                        ((("C8H10N4O2", 195.0877, None, [])), "no peaks"), ((("C8H10N4O2", 195.0877, "[M-H]-", PEAKS)), "positive ion mode only"),
                        ((("C8H10N4O2", 195.0877, None, PEAKS, "slow")), "Mode")):
        with pytest.raises(ValueError, match=match):
            fe.predict(*args)
    with pytest.raises(ValueError, match="at most"):
        fe.predict("C8H10N4O2", 195.0877, None, [(100.0 + i * 0.01, 1.0) for i in range(fe.MAX_PEAKS + 1)])


def test_online_server_without_frigid_answers_502_with_a_clear_message(monkeypatch, tmp_path):
    monkeypatch.setenv("MZDESK_FRIGID_DIR", str(tmp_path / "nothing"))
    r = _run(_inputs())
    assert r.status_code == 502
    detail = r.json()["detail"]
    assert "local version" in detail and fe.SETUP_COMMAND in detail


# ---------------------------------------------------------------- driver result -> engine error
def test_driver_result_ok_passes_through():
    assert fe.parse_result({"ok": True, "candidates": []})["ok"] is True


def test_invalid_input_exit_is_a_value_error():
    with pytest.raises(ValueError, match="rejected the input: adduct 'x' is not one"):
        fe.parse_result({"ok": False, "error_type": "invalid_input", "error": "adduct 'x' is not one of the positive adducts"}, code=2)


def test_memory_guard_exit_is_an_engine_error_with_advice():
    with pytest.raises(fe.FrigidError, match="memory guard.*cpu: process-tree RSS 6200 MB > cap 6144 MB.*Quick mode") as exc:
        fe.parse_result({"ok": False, "error_type": "memory_guard", "error": "cpu: process-tree RSS 6200 MB > cap 6144 MB"}, "tail", 3)
    assert exc.value.log == "tail"


def test_runtime_failure_keeps_the_traceback_as_the_log():
    with pytest.raises(fe.FrigidError, match="FRIGID failed: RuntimeError: boom") as exc:
        fe.parse_result({"ok": False, "error_type": "RuntimeError", "error": "boom", "traceback": "Traceback ...\nRuntimeError: boom"}, "tail", 1)
    assert exc.value.log.startswith("Traceback")


def test_no_result_file_is_an_engine_error_with_the_output(fake_install):
    seen = fake_install(result="none", code=-9, log="Killed: 9")
    with pytest.raises(fe.FrigidError, match="without writing a result \\(exit code -9\\)") as exc:
        fe.predict("C8H10N4O2", 195.0877, None, PEAKS)
    assert exc.value.log == "Killed: 9" and not seen["tmp"].exists()


@pytest.mark.parametrize("result,exc_type,match", [
    ({"ok": False, "error_type": "invalid_input", "error": "no usable peaks"}, ValueError, "no usable peaks"),
    ({"ok": False, "error_type": "memory_guard", "error": "mps: footprint 12400 MB > cap 12288 MB"}, fe.FrigidError, "memory guard"),
    ({"ok": False, "error_type": "FileNotFoundError", "error": "DLM.ckpt", "traceback": "tb"}, fe.FrigidError, "FileNotFoundError: DLM.ckpt")])
def test_predict_maps_driver_failures(fake_install, result, exc_type, match):
    fake_install(result=result, code=2)
    with pytest.raises(exc_type, match=match):
        fe.predict("C8H10N4O2", 195.0877, None, PEAKS)


# ---------------------------------------------------------------- predict() plumbing with a stand-in for the process
def test_predict_runs_the_shipped_driver_with_the_env_python_and_filters(fake_install, tmp_path):
    seen = fake_install()
    res = fe.predict("H10C8N4O2", 195.0877, "[M+NH4]+", PEAKS, "Full (up to 100 samples)", name="caffeine", seed=3, device="cpu")
    cmd = seen["cmd"]
    assert cmd[:4] == [str(tmp_path / "env" / "bin" / "python"), "-I", str(fe.DRIVER), "--input"] and cmd[5] == "--output"
    assert Path(cmd[4]).name == "job.json" and Path(cmd[6]).name == "result.json"
    job = seen["job"]
    assert (job["formula"], job["adduct"], job["mode"], job["device"], job["seed"], job["ion_mode"]) == ("C8H10N4O2", "[M+H3N+H]+", "full", "cpu", 3, "positive")
    assert job["length_model"] == str(fe.LENGTH_MODEL) and Path(job["length_model"]).is_file()
    assert seen["cwd"] == str(seen["tmp"]) and seen["timeout"] == fe.TIMEOUT_S["full"]
    env = seen["env"]
    assert env["HF_HOME"] == str(tmp_path / "hf_cache") and env["HF_HUB_OFFLINE"] == "1" and env["FRIGID_ROOT"] == str(tmp_path / "FRIGID")
    assert res["adduct"] == "[M+NH4]+" and res["driver_adduct"] == "[M+H3N+H]+" and res["formula"] == "C8H10N4O2" and res["mode"] == "full"
    assert [c["rank"] for c in res["candidates"]] == [1, 2, 3, 4] and res["dropped"]["not_formula_match"] == 3 and res["n_raw"] == 7
    assert (res["device"], res["device_requested"]) == ("mps", "cpu")                       # the driver's answer is reported, not the request
    assert res["n_samples"] == 16 and res["n_unique_matched"] == 4 and res["n_peaks_assigned"] == 27 and res["length_model"] is True
    assert res["peak_footprint_mb"] == 5047.4 and res["driver_seconds"] == 17.278 and "auto: mps" in res["device_reason"]
    assert not seen["tmp"].exists()                                                          # the temp folder is gone


def test_predict_quick_has_the_short_timeout_and_auto_device(fake_install, monkeypatch):
    monkeypatch.delenv("MZDESK_FRIGID_DEVICE", raising=False)
    seen = fake_install()
    fe.predict("C8H10N4O2", 195.0877, None, PEAKS)
    assert seen["timeout"] == 300 and seen["job"]["device"] == "auto" and seen["job"]["mode"] == "quick" and seen["job"]["adduct"] == "[M+H]+"
    monkeypatch.setenv("MZDESK_FRIGID_DEVICE", "cpu")
    fe.predict("C8H10N4O2", 195.0877, None, PEAKS)
    assert seen["job"]["device"] == "cpu"


# ---------------------------------------------------------------- the shared heavy-job lock
def test_msfinder_and_frigid_share_one_lock():
    assert mf._RUN_LOCK is heavy_lock.LOCK


def test_frigid_is_turned_away_while_msfinder_runs(fake_install, monkeypatch):
    fake_install()
    monkeypatch.setattr(fe, "BUSY_WAIT_S", 0.05)
    assert heavy_lock.LOCK.acquire(holder="MS-FINDER")
    try:
        with pytest.raises(fe.FrigidError, match="busy: MS-FINDER"):
            fe.predict("C8H10N4O2", 195.0877, None, PEAKS)
    finally:
        heavy_lock.LOCK.release()
    assert fe.predict("C8H10N4O2", 195.0877, None, PEAKS)["candidates"]                     # free again afterwards


def test_frigid_is_turned_away_while_another_frigid_run_is_in_progress(fake_install, monkeypatch):
    fake_install()
    monkeypatch.setattr(fe, "BUSY_WAIT_S", 0.05)
    assert heavy_lock.LOCK.acquire(holder="FRIGID")
    try:
        with pytest.raises(fe.FrigidError, match="busy: another De Novo Structure run"):
            fe.predict("C8H10N4O2", 195.0877, None, PEAKS)
    finally:
        heavy_lock.LOCK.release()


def test_msfinder_waits_for_frigid_and_says_so(monkeypatch):
    monkeypatch.setattr(mf, "check_msfinder", lambda: ("/usr/bin/mono", Path("x.exe"), Path("x.ini")))
    monkeypatch.setattr(mf, "BUSY_WAIT_S", 0.05)
    assert heavy_lock.LOCK.acquire(holder="FRIGID")
    try:
        with pytest.raises(mf.MsFinderError, match="waiting for a De Novo Structure \\(FRIGID\\) run"):
            mf.predict(195.0877, True, None, PEAKS)
    finally:
        heavy_lock.LOCK.release()


def test_msfinder_keeps_its_own_busy_message_when_nobody_is_named(monkeypatch):
    monkeypatch.setattr(mf, "check_msfinder", lambda: ("/usr/bin/mono", Path("x.exe"), Path("x.ini")))
    monkeypatch.setattr(mf, "BUSY_WAIT_S", 0.05)
    with mf._RUN_LOCK:
        with pytest.raises(mf.MsFinderError, match="busy with another request"):
            mf.predict(195.0877, True, None, PEAKS)


def test_lock_is_released_after_a_failed_run(fake_install):
    fake_install(result={"ok": False, "error_type": "RuntimeError", "error": "boom"}, code=1)
    with pytest.raises(fe.FrigidError):
        fe.predict("C8H10N4O2", 195.0877, None, PEAKS)
    assert heavy_lock.LOCK.acquire(blocking=False)
    assert heavy_lock.LOCK.holder is None
    heavy_lock.LOCK.release()


def test_two_requests_never_run_at_the_same_time(fake_install, monkeypatch):
    """Two threads ask at once: the second waits for the first (stand-in process sleeps 0.3 s), so at most one is ever running."""
    fake_install()
    running, peak, guard = [0], [0], threading.Lock()

    def slow_run(cmd, cwd, env, timeout, **kw):
        with guard:
            running[0] += 1
            peak[0] = max(peak[0], running[0])
        time.sleep(0.3)
        Path(cmd[cmd.index("--output") + 1]).write_text(json.dumps(DRIVER_RESULT))
        with guard:
            running[0] -= 1
        return 0, "", {"peak_rss_mb": 1.0, "peak_footprint_mb": 1.0}
    monkeypatch.setattr(fe, "run_watched", slow_run)
    out = []
    threads = [threading.Thread(target=lambda: out.append(fe.predict("C8H10N4O2", 195.0877, None, PEAKS))) for _ in range(2)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=30)
    assert len(out) == 2 and peak[0] == 1


# ---------------------------------------------------------------- process watching (no FRIGID involved)
def _python(code):
    return [sys.executable, "-I", "-c", code]


def _alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    return True


def test_run_watched_returns_exit_code_output_and_peaks():
    rc, log, peaks = fe.run_watched(_python("import sys; print('hello'); print('oops', file=sys.stderr); sys.exit(3)"), ".", dict(os.environ), 30, poll_s=0.1)
    assert rc == 3 and "hello" in log and "oops" in log and peaks["peak_rss_mb"] >= 0 and peaks["peak_footprint_mb"] >= 0


def test_run_watched_timeout_kills_the_whole_process_group():
    code = "import subprocess, sys, time; p = subprocess.Popen(['sleep', '60']); print(p.pid, flush=True); time.sleep(60)"
    with pytest.raises(fe.FrigidError, match="did not finish within 1 s") as exc:
        fe.run_watched(_python(code), ".", dict(os.environ), 1, poll_s=0.1)
    child = int(exc.value.log.split()[0])
    for _ in range(50):
        if not _alive(child):
            break
        time.sleep(0.1)
    assert not _alive(child)


@pytest.mark.skipif(sys.platform != "darwin", reason="the memory guard reads the kernel's per-process numbers (libproc, macOS)")
@pytest.mark.parametrize("caps", [{"rss_cap_mb": 120}, {"footprint_cap_mb": 120}], ids=["rss", "footprint"])
def test_run_watched_memory_caps_kill_the_run(caps):
    code = "import time; x = bytearray(200 << 20); x[::4096] = b'1' * len(x[::4096]); print('allocated', flush=True); time.sleep(60)"
    t0 = time.monotonic()
    with pytest.raises(fe.FrigidError, match="too much memory") as exc:
        fe.run_watched(_python(code), ".", dict(os.environ), 60, poll_s=0.1, **caps)
    assert time.monotonic() - t0 < 30 and isinstance(exc.value.log, str)


# ---------------------------------------------------------------- runner
def test_runner_fills_both_blocks_and_says_where_the_formula_came_from(fake_predict):
    r = _run(_inputs())
    assert r.status_code == 200, r.text
    body = r.json()
    out = body["outputs"]
    values = out["summary"]["values"]
    assert values["Top candidate"]["value"] == SMILES_CAFFEINE and values["Top candidate"]["note"] == "score 1.000 · C8H10N4O2"
    assert values["Formula used"] == {"value": "C8H10N4O2", "note": "entered by you"}
    assert values["Device"]["value"] == "MPS (Apple GPU)" and values["Device"]["note"].startswith("21 s · peak memory 4.9 GB · auto")
    assert any("auto: mps" in line for line in values["Device"]["help"])
    assert values["Samples"]["value"] == "16 / 4" and "generated / shown" in values["Samples"]["note"]
    rows = out["candidates"]["rows"]
    assert out["candidates"]["count"] == 4 and len(rows) == 4
    assert rows[0][0] == "1" and rows[0][1]["img"].startswith("data:image/png;base64,") and rows[0][2:] == [SMILES_CAFFEINE, "C8H10N4O2", "1.000", "RYYVLZVUVIJVGH"]
    assert [r_[4] for r_ in rows] == ["1.000", "0.515", "0.500", "0.333"]
    assert "Hidden" not in " ".join(body["warnings"])
    assert body["verdict"] == f"Top generated candidate: <b>{SMILES_CAFFEINE}</b>, similarity 1.000 · 4 formula-matched candidates (hypotheses)."


def test_runner_passes_formula_adduct_mode_and_peaks_to_the_engine(fake_predict):
    _run(_inputs(precursor={**PRECURSOR, "adduct": "[M+Na]+"}, formula="H10C8N4O2", mode="Full (up to 100 samples)"))
    assert (fake_predict["formula"], fake_predict["adduct"], fake_predict["mode"]) == ("C8H10N4O2", "[M+Na]+", "full")
    assert fake_predict["precursor_mz"] == 195.0877 and len(fake_predict["peaks"]) == 27 and fake_predict["peaks"][0] == (54.0339, 7.0)
    _run(_inputs())
    assert (fake_predict["adduct"], fake_predict["mode"]) == ("[M+H]+", "quick")             # Auto -> [M+H]+, default Quick
    _run({k: v for k, v in _inputs().items() if k != "formula"} | {"formula": "C8H10N4O2", "mode": "Quick (16 samples)"})
    assert fake_predict["mode"] == "quick"


def test_runner_takes_the_precursor_from_the_spectrum_header(fake_predict):
    mgf = "BEGIN IONS\nPEPMASS=195.0877\n56.0495 88\n110.0712 590\n138.0662 999\nEND IONS"
    body = _run(_inputs(precursor={**PRECURSOR, "precursor_mz": ""}, ms2=_paste(mgf))).json()
    assert fake_predict["precursor_mz"] == 195.0877 and any("read from the spectrum" in w for w in body["warnings"])


def test_runner_negative_mode_is_a_400_that_says_positive_only(fake_predict):
    r = _run(_inputs(precursor={**PRECURSOR, "ion_mode": "Negative"}))
    assert r.status_code == 400 and "positive ion mode only" in r.json()["detail"] and "formula" not in fake_predict
    r = _run(_inputs(precursor={**PRECURSOR, "adduct": "[M-H]-"}))
    assert r.status_code == 400 and "positive ion mode only" in r.json()["detail"]


@pytest.mark.parametrize("over,match", [
    ({"formula": "C8H10N4O2Zn"}, "cannot use Zn"), ({"formula": "c8h10n4o2"}, "Invalid formula"), ({"formula": "C8H10N4O2+"}, "Invalid formula"),
    ({"formula": "H2O"}, "no carbon"), ({"formula": "caffeine"}, "Invalid formula"), ({"mode": "Slow"}, "Mode must be one of"),
    ({"precursor": {**PRECURSOR, "adduct": "[M+Li]+"}}, "does not support the adduct"), ({"precursor": {**PRECURSOR, "precursor_mz": ""}}, "Precursor m/z is required"),
    ({"ms2": _paste("not a spectrum")}, "no valid peaks")])
def test_runner_bad_input_is_400(fake_predict, over, match):
    r = _run(_inputs(**over))
    assert r.status_code == 400 and re.search(match, r.json()["detail"]), r.text
    assert "formula" not in fake_predict                                                    # the engine was never reached


def test_runner_warns_when_the_formula_does_not_fit_the_precursor(fake_predict):
    body = _run(_inputs(precursor={**PRECURSOR, "precursor_mz": "215.0000"})).json()
    assert any("C8H10N4O2 as [M+H]+ has m/z 195.0877" in w and "Check the formula and the adduct" in w for w in body["warnings"])
    ok = _run(_inputs(precursor={**PRECURSOR, "adduct": "[M+Na]+", "precursor_mz": "217.0696"})).json()      # caffeine + Na
    assert not any("has m/z" in w for w in ok["warnings"])
    assert not any("has m/z" in w for w in _run(_inputs()).json()["warnings"])


def test_runner_empty_result_is_an_empty_table_and_a_warning(fake_predict):
    dropped = {"not_formula_match": 5, "isotope": 1, "radical": 2, "charged": 0, "duplicate": 0}
    fake_predict["result"] = _fake_result(candidates=[], dropped=dropped, warnings=["no generated molecule matches the requested formula; candidates are padding"])
    body = _run(_inputs()).json()
    out = body["outputs"]
    assert out["candidates"] == {"rows": [], "count": 0} and out["summary"]["values"]["Top candidate"] == "—"
    assert out["summary"]["values"]["Samples"]["value"] == "16 / 0"
    text = " ".join(body["warnings"])
    assert "no usable structure for C8H10N4O2 in 16 samples (1 isotope-labelled, 2 radical hidden)" in text and "Try Full mode" in text
    assert "padding" not in text                                                            # the driver's own wording is replaced
    assert body["verdict"] == "FRIGID generated no formula-matched structure for <b>C8H10N4O2</b> (see the warnings)."


def test_runner_says_what_was_hidden_and_passes_driver_notes_through(fake_predict):
    kept = _fake_result()["candidates"]
    fake_predict["result"] = _fake_result(candidates=kept[:1], dropped={"not_formula_match": 3, "isotope": 0, "radical": 1, "charged": 1, "duplicate": 2},
                                          warnings=["only 2 peaks assigned a sub-formula; fingerprint prediction will be weak"])
    warnings = _run(_inputs()).json()["warnings"]
    assert "Hidden: 1 radical, 1 wrong net charge, 2 duplicate (not valid candidates for this formula and adduct)." in warnings
    assert "Only 2 peaks assigned a sub-formula; fingerprint prediction will be weak" in warnings
    assert any("Only 1 formula-matched structure in 16 samples; Full mode" in w for w in warnings)


def test_runner_warns_when_no_peak_fits_the_formula(fake_predict):
    fake_predict["result"] = _fake_result(n_peaks_assigned=0)
    assert any("No peak could be assigned" in w for w in _run(_inputs()).json()["warnings"])


def test_runner_reports_the_cpu_and_a_long_smiles_is_cut_in_the_card(fake_predict):
    long = "C" * 80
    fake_predict["result"] = _fake_result(device="cpu", candidates=[{"rank": 1, "smiles": long, "formula": "C80H162", "score": 0.5, "inchikey_2d": "AAAAAAAAAAAAAA"}])
    body = _run(_inputs()).json()
    top = body["outputs"]["summary"]["values"]["Top candidate"]
    assert top["value"].endswith("…") and len(top["value"]) == 44 and top["help"] == ["SMILES: " + long]
    assert body["outputs"]["summary"]["values"]["Device"]["value"] == "CPU" and body["outputs"]["candidates"]["rows"][0][2] == long


def test_engine_failures_are_502_with_the_log_and_bad_input_is_400(monkeypatch):
    def broken(*a, **k):
        raise fe.FrigidError("FRIGID did not finish within 300 s and was stopped.", "last lines")
    monkeypatch.setattr(fe, "predict", broken)
    r = _run(_inputs())
    assert r.status_code == 502 and r.json() == {"detail": "FRIGID did not finish within 300 s and was stopped.", "log": "last lines"}

    def bad(*a, **k):
        raise ValueError("FRIGID rejected the input: no usable peaks")
    monkeypatch.setattr(fe, "predict", bad)
    assert _run(_inputs()).status_code == 400


# ---------------------------------------------------------------- BUDDY fallback (BUDDY replaced)
def _buddy(monkeypatch, candidates, error=None):
    from mzdesk.ms1 import msbuddy_engine as bd
    calls = []

    def fake(precursor_mz, positive, adduct, rt, formula_input, ms2_peaks, ms1_peaks, ms_instr="orbitrap", halogen=False):
        calls.append(dict(precursor_mz=precursor_mz, positive=positive, adduct=adduct, rt=rt, formula=formula_input, ms2=ms2_peaks, ms1=ms1_peaks,
                          instr=ms_instr, halogen=halogen))
        if error:
            raise error
        return {"candidates": candidates, "adduct_used": adduct, "from_input": False, "ms1_used": 0}
    monkeypatch.setattr(bd, "annotate_formula", fake)
    return calls


def _no_formula(**over):
    inputs = _inputs(**over)
    inputs.pop("formula")
    return inputs


def test_blank_formula_uses_the_buddy_top_formula_and_says_so(monkeypatch, fake_predict):
    calls = _buddy(monkeypatch, [{"rank": 1, "formula": "C8H10N4O2", "prob": 0.9734}, {"rank": 2, "formula": "C9H14N2O3", "prob": 0.02}])
    body = _run(_no_formula()).json()
    assert fake_predict["formula"] == "C8H10N4O2"
    assert len(calls) == 1 and calls[0]["positive"] is True and calls[0]["adduct"] == "[M+H]+" and calls[0]["formula"] == "" and calls[0]["ms1"] == ()
    assert calls[0]["precursor_mz"] == 195.0877 and len(calls[0]["ms2"]) == 27 and calls[0]["instr"] == "orbitrap" and calls[0]["halogen"] is False
    assert body["outputs"]["summary"]["values"]["Formula used"] == {"value": "C8H10N4O2", "note": "BUDDY top-1 of 2 · probability 0.97"}
    text = " ".join(body["warnings"])
    assert "No formula was entered, so BUDDY's top-1 formula C8H10N4O2 (probability 0.97) was used" in text and "C9H14N2O3" in text
    assert "Enter the formula to override it" in text


def test_a_typed_formula_never_calls_buddy(monkeypatch, fake_predict):
    calls = _buddy(monkeypatch, [], error=AssertionError("BUDDY must not run"))
    assert _run(_inputs()).status_code == 200 and calls == []


def test_buddy_adduct_follows_the_chosen_one(monkeypatch, fake_predict):
    calls = _buddy(monkeypatch, [{"rank": 1, "formula": "C8H10N4O2", "prob": None}])
    body = _run(_no_formula(precursor={**PRECURSOR, "adduct": "[M+Na]+", "precursor_mz": "217.0696"})).json()
    assert calls[0]["adduct"] == "[M+Na]+" and fake_predict["adduct"] == "[M+Na]+"
    assert body["outputs"]["summary"]["values"]["Formula used"]["note"] == "BUDDY top-1 of 1"           # no probability: nothing invented


@pytest.mark.parametrize("adduct", ["[M]+", "[M+H-2H2O]+"])
def test_adducts_buddy_cannot_do_need_a_typed_formula(monkeypatch, fake_predict, adduct):
    calls = _buddy(monkeypatch, [{"rank": 1, "formula": "C8H10N4O2", "prob": 0.9}])
    r = _run(_no_formula(precursor={**PRECURSOR, "adduct": adduct}))
    assert r.status_code == 400 and "BUDDY cannot infer one" in r.json()["detail"] and "Enter the molecular formula" in r.json()["detail"]
    assert calls == [] and "formula" not in fake_predict


def test_buddy_finding_nothing_asks_for_a_formula(monkeypatch, fake_predict):
    _buddy(monkeypatch, [])
    r = _run(_no_formula())
    assert r.status_code == 400 and "BUDDY found none for m/z 195.0877 as [M+H]+" in r.json()["detail"] and "enter the molecular formula" in r.json()["detail"]
    assert "formula" not in fake_predict


def test_buddy_out_of_range_precursor_is_a_400(monkeypatch, fake_predict):
    _buddy(monkeypatch, [], error=ValueError("BUDDY supports precursor m/z between 1 and 1500."))
    r = _run(_no_formula(precursor={**PRECURSOR, "precursor_mz": "2000"}))
    assert r.status_code == 400 and "BUDDY supports precursor m/z between 1 and 1500" in r.json()["detail"]


def test_buddy_formula_that_frigid_cannot_use_is_a_400(monkeypatch, fake_predict):
    _buddy(monkeypatch, [{"rank": 1, "formula": "H2O", "prob": 0.5}])
    r = _run(_no_formula())
    assert r.status_code == 400 and "BUDDY's top formula H2O cannot be used by FRIGID" in r.json()["detail"]


# ---------------------------------------------------------------- page definition
def test_page_definition_is_live_with_the_agreed_inputs_and_blocks():
    m = _module()
    assert m["live"] is True and m["subgroup"] == "ms2" and m["title"] == "De Novo Structure"
    assert [f["id"] for f in m["inputs"]] == ["precursor", "ms2", "ms1", "formula", "mode", "halogen"]
    mode = m["inputs"][4]
    assert mode["options"] == ["Quick (16 samples)", "Full (up to 100 samples)"] and mode["value"] == "Quick (16 samples)"
    assert [b["id"] for b in m["outputs"]] == ["summary", "candidates"]
    assert [i["label"] for i in m["outputs"][0]["items"]] == ["Top candidate", "Formula used", "Device", "Samples"]
    table = m["outputs"][1]
    assert table["columns"] == ["Rank", "Structure", "SMILES", "Formula", "Score", "InChIKey (2D)"] and "Tool" not in table["columns"]
    assert m["example"]["formula"] == "C8H10N4O2" and m["example"]["precursor"]["precursor_mz"] == "195.0877"
    assert m["example"]["ms2"]["value"] == CAFFEINE_MSFINDER_MS2


def test_help_text_is_honest_about_what_the_page_gives():
    text = " ".join(_module()["outputs"][1]["help"])
    for needle in ("hypotheses", "8.21", "9.82", "13.31", "20.15", "likely lower", "Tanimoto", "not a probability", "stereochemistry", "saw it in training",
                   "BUDDY", "frigid_setup"):
        assert needle in text, needle
    assert "Local version only" in _module()["one_liner"] and "hypotheses" in _module()["one_liner"]


# ---------------------------------------------------------------- shipped files
def test_driver_and_length_model_ship_with_the_package():
    drv = fe.DRIVER.parent
    assert {p.name for p in drv.glob("*.py")} == {"frigid_predict.py", "subformula_assign.py", "memguard.py"}
    for p in drv.glob("*.py"):
        ast.parse(p.read_text(encoding="utf-8"))
    assert "samgoldman97/mist" in (drv / "subformula_assign.py").read_text() and "MIT licence" in (drv / "subformula_assign.py").read_text()
    assert fe.LENGTH_MODEL.stat().st_size > 100_000 and json.loads(fe.LENGTH_FEATURES.read_text()) == \
        ["As", "B", "Br", "C", "Cl", "F", "H", "I", "N", "O", "P", "S", "Se", "Si"]
    pyproject = (ROOT / "pyproject.toml").read_text()
    assert '"lib/frigid/*"' in pyproject and '"ms2/frigid_driver/*.py"' in pyproject


def test_memguard_copies_are_identical():
    """mzdesk/shared/memguard.py (used by the engine) and the driver's own copy (runs in the FRIGID env, imports nothing from mzdesk) are one file."""
    a = (ROOT / "mzdesk" / "shared" / "memguard.py").read_bytes()
    assert a == (fe.DRIVER.parent / "memguard.py").read_bytes()


def test_memguard_chooses_a_device_and_reads_this_process():
    from mzdesk.shared import memguard
    assert memguard.choose_device("cpu", mps_available=True)["chosen"] == "cpu"
    assert memguard.choose_device("auto", mps_available=False)["chosen"] == "cpu"
    if sys.platform == "darwin":
        assert memguard.proc_memory(os.getpid())["rss_mb"] > 1


# ---------------------------------------------------------------- setup command
def test_setup_check_reports_what_is_missing(monkeypatch, tmp_path, capsys):
    monkeypatch.setenv("MZDESK_FRIGID_DIR", str(tmp_path / "nothing"))
    assert su.check() == 1
    out = capsys.readouterr().out
    assert "NOT FOUND" in out and "NOT APPLIED" in out and fe.SETUP_COMMAND in out and str(tmp_path / "nothing") in out and "Ready              : no" in out


@pytest.mark.skipif(NO_FRIGID is not None, reason="FRIGID not installed")
def test_setup_check_is_ready_on_an_installed_machine(capsys):
    assert su.main(["--check"]) == 0
    out = capsys.readouterr().out
    assert "Ready              : yes" in out and "Device patch       : applied" in out and "SAFE import patch  : applied" in out


def test_setup_check_never_starts_the_install(monkeypatch, tmp_path):
    monkeypatch.setenv("MZDESK_FRIGID_DIR", str(tmp_path / "nothing"))
    monkeypatch.setattr(su, "setup", lambda *a, **k: pytest.fail("--check must not install"))
    monkeypatch.setattr(su, "download", lambda *a, **k: pytest.fail("--check must not download"))
    assert su.main(["--check"]) == 1 and not (tmp_path / "nothing").exists()


def test_setup_dir_option_sets_the_install_folder(monkeypatch, tmp_path, capsys):
    monkeypatch.setenv("MZDESK_FRIGID_DIR", str(tmp_path / "before"))                 # main() overwrites it; monkeypatch puts the original back
    su.main(["--check", "--dir", str(tmp_path / "other")])
    assert str(tmp_path / "other") in capsys.readouterr().out


def test_pinned_requirements_are_exact_and_leave_out_the_training_and_conda_parts():
    names = {}
    for line in fe.REQUIREMENTS.read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            assert re.fullmatch(r"[A-Za-z0-9_.\-]+==[^\s=]+", line), line
            names[line.split("==")[0].lower()] = line.split("==")[1]
    assert names["torch"] == "2.6.0" and names["transformers"] == "4.56.2" and names["numpy"] == "1.26.4" and "safe-mol" in names and "ngboost" in names
    assert not {"optuna", "alembic", "mako", "sqlalchemy", "colorlog"} & set(names) and not any(n.startswith("antlr") for n in names)


def test_env_commands_name_the_pinned_pieces(tmp_path):
    cmd = su.env_create_command("/opt/conda/bin/conda", tmp_path / "env")
    assert cmd[:5] == ["/opt/conda/bin/conda", "create", "-y", "-p", str(tmp_path / "env")] and "python=3.10" in cmd and "antlr-python-runtime=4.9.3" in cmd
    assert "--override-channels" in cmd and "conda-forge" in cmd
    pip = su.pip_install_command(tmp_path / "env" / "bin" / "python")
    assert pip[1:6] == ["-m", "pip", "install", "--only-binary=:all:", "-r"] and pip[-1] == str(fe.REQUIREMENTS)


def test_setup_without_conda_says_what_to_install(monkeypatch, tmp_path):
    monkeypatch.setattr(su, "find_conda", lambda: None)
    with pytest.raises(su.SetupError, match="conda was not found"):
        su.ensure_env(tmp_path)


def test_pins_match_the_documented_sources():
    assert fe.FRIGID_COMMIT.startswith("9a3e623") and su.DLM_TAR_MD5 == "5a663a3000664ded0cc62fc5a4ca39e0" and su.DLM_MD5 == "717f93448661d9b63c4ceb4f7bcdfab0"
    assert su.MIST_MD5 == "3e92179240e3043733627d720f4b7b86" and fe.DLM_SIZE == 2_676_432_167 and fe.MIST_SIZE == 289_222_437
    assert "records/23543494" in su.DLM_URL and "records/23303289" in su.MIST_URL and su.REPO_URL == "https://github.com/coleygroup/FRIGID"


def _md5(data):
    return hashlib.md5(data).hexdigest()


def test_download_verifies_and_skips_what_is_already_correct(tmp_path):
    src = tmp_path / "src.bin"
    src.write_bytes(b"checkpoint" * 1000)
    dest = tmp_path / "out" / "ckpt.bin"
    assert su.fetch_file(src.as_uri(), dest, 10_000, _md5(src.read_bytes())) is True and dest.read_bytes() == src.read_bytes() and not dest.with_name("ckpt.bin.part").exists()
    assert su.fetch_file("file:///does/not/exist", dest, 10_000, _md5(src.read_bytes())) is False                  # already correct: no download attempted


def test_download_restarts_when_the_server_cannot_resume(tmp_path):
    src = tmp_path / "src.bin"
    src.write_bytes(b"0123456789" * 100)
    dest = tmp_path / "ckpt.bin"
    dest.with_name("ckpt.bin.part").write_bytes(b"junk")                       # an earlier, interrupted attempt; file:// ignores the Range header
    su.download(src.as_uri(), dest)
    assert dest.read_bytes() == src.read_bytes()


def test_download_with_the_wrong_md5_is_refused_and_set_aside(tmp_path):
    src = tmp_path / "src.bin"
    src.write_bytes(b"abc" * 100)
    dest = tmp_path / "ckpt.bin"
    with pytest.raises(su.SetupError, match="md5"):
        su.fetch_file(src.as_uri(), dest, 300, "0" * 32)
    assert not dest.exists() and dest.with_name("ckpt.bin.bad").exists()


def test_a_damaged_existing_file_is_not_silently_replaced(tmp_path):
    dest = tmp_path / "ckpt.bin"
    dest.write_bytes(b"half")
    with pytest.raises(su.SetupError, match="size or md5 is wrong"):
        su.fetch_file("file:///unused", dest, 300, "0" * 32)
    assert dest.read_bytes() == b"half"


def _tar_gz(path, members):
    with tarfile.open(path, "w:gz") as tf:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))


def test_only_the_dlm_checkpoint_is_taken_out_of_the_tar(tmp_path):
    dlm = b"D" * 5000
    tar = tmp_path / "ckpts.tar.gz"
    _tar_gz(tar, {"frigid_pretrained_checkpoints/iceberg_a.ckpt": b"I" * 4000, "frigid_pretrained_checkpoints/DLM.ckpt": dlm, "../evil/DLM.ckpt.bak": b"x"})
    dest = tmp_path / "checkpoints" / "DLM.ckpt"
    su.extract_member(tar, "DLM.ckpt", dest, len(dlm), _md5(dlm))
    assert dest.read_bytes() == dlm and sorted(p.name for p in dest.parent.iterdir()) == ["DLM.ckpt"] and not (tmp_path / "evil").exists()


def test_extraction_refuses_a_checkpoint_that_does_not_match(tmp_path):
    tar = tmp_path / "ckpts.tar.gz"
    _tar_gz(tar, {"x/DLM.ckpt": b"D" * 100})
    with pytest.raises(su.SetupError, match="expected 100 bytes"):
        su.extract_member(tar, "DLM.ckpt", tmp_path / "out" / "DLM.ckpt", 100, "0" * 32)
    assert not (tmp_path / "out" / "DLM.ckpt").exists() and not (tmp_path / "out" / "DLM.ckpt.part").exists()
    with pytest.raises(su.SetupError, match="not in"):
        su.extract_member(tar, "missing.ckpt", tmp_path / "out" / "m.ckpt", 1, "0" * 32)


def test_fetch_checkpoints_extracts_dlm_and_deletes_the_tar(tmp_path, monkeypatch):
    dlm, mist = b"D" * 3000, b"M" * 2000
    monkeypatch.setattr(fe, "DLM_SIZE", len(dlm))
    monkeypatch.setattr(fe, "MIST_SIZE", len(mist))
    monkeypatch.setattr(su, "DLM_MD5", _md5(dlm))
    monkeypatch.setattr(su, "MIST_MD5", _md5(mist))
    (tmp_path / "src").mkdir()
    (tmp_path / "src" / "mist.ckpt").write_bytes(mist)
    _tar_gz(tmp_path / "src" / "pack.tar.gz", {"pack/DLM.ckpt": dlm, "pack/other.ckpt": b"o" * 10})
    tar_bytes = (tmp_path / "src" / "pack.tar.gz").read_bytes()
    monkeypatch.setattr(su, "MIST_URL", (tmp_path / "src" / "mist.ckpt").as_uri())
    monkeypatch.setattr(su, "DLM_URL", (tmp_path / "src" / "pack.tar.gz").as_uri())
    monkeypatch.setattr(su, "DLM_TAR_SIZE", len(tar_bytes))
    monkeypatch.setattr(su, "DLM_TAR_MD5", _md5(tar_bytes))
    install = tmp_path / "install"
    su.fetch_checkpoints(install)
    assert (install / "checkpoints" / "DLM.ckpt").read_bytes() == dlm and (install / "checkpoints" / "mist_msg_aug.ckpt").read_bytes() == mist
    assert not (install / "downloads" / su.DLM_TAR).exists()
    monkeypatch.setattr(su, "DLM_URL", "file:///gone")                               # second run: everything verified, nothing fetched
    monkeypatch.setattr(su, "MIST_URL", "file:///gone")
    su.fetch_checkpoints(install)


def _fake_env_with_safe(install, text="from . import trainer, utils\nfrom .sample import SAFEDesign\n"):
    init = install / "env" / "lib" / "python3.10" / "site-packages" / "safe" / "__init__.py"
    init.parent.mkdir(parents=True)
    init.write_text(text)
    return init


def test_safe_patch_comments_everything_out_once_and_keeps_a_backup(tmp_path):
    original = "from . import trainer, utils\nfrom .sample import SAFEDesign\n"
    init = _fake_env_with_safe(tmp_path, original)
    assert su.apply_safe_patch(tmp_path) == "applied"
    assert init.read_text() == "# from . import trainer, utils\n# from .sample import SAFEDesign\nfrom .converter import SAFEConverter, decode, encode\n"
    assert init.with_name("__init__.py.orig").read_text() == original and fe.safe_patch_applied(init)
    assert su.apply_safe_patch(tmp_path) == "already applied" and init.with_name("__init__.py.orig").read_text() == original


def test_safe_patch_without_the_package_is_a_clear_error(tmp_path):
    with pytest.raises(su.SetupError, match="safe package is not installed"):
        su.apply_safe_patch(tmp_path)


def test_tokenizer_already_cached_is_not_fetched_again(tmp_path, monkeypatch):
    snap = tmp_path / "hf_cache" / "hub" / "models--datamol-io--safe-gpt" / "snapshots" / fe.TOKENIZER_REVISION
    snap.mkdir(parents=True)
    (snap / "tokenizer.json").write_bytes(b"t" * su.TOKENIZER_SIZE)
    monkeypatch.setattr(su, "_run", lambda *a, **k: pytest.fail("must not call the hub"))
    assert su.fetch_tokenizer(tmp_path) is False
    assert (tmp_path / "hf_cache" / "hub" / "models--datamol-io--safe-gpt" / "refs" / "main").read_text() == fe.TOKENIZER_REVISION       # offline lookup of "main"


def test_tokenizer_fetch_runs_in_the_env_python_and_pins_the_revision(tmp_path, monkeypatch):
    calls = []

    def fake_run(cmd, **kw):
        calls.append((cmd, kw))
        snap = tmp_path / "hf_cache" / "hub" / "models--datamol-io--safe-gpt" / "snapshots" / fe.TOKENIZER_REVISION
        snap.mkdir(parents=True)
        (snap / "tokenizer.json").write_bytes(b"t" * su.TOKENIZER_SIZE)
        return subprocess.CompletedProcess(cmd, 0, "", "")
    monkeypatch.setattr(su, "_run", fake_run)
    monkeypatch.setenv("HF_HUB_OFFLINE", "1")
    assert su.fetch_tokenizer(tmp_path) is True
    cmd, kw = calls[0]
    assert cmd[:3] == [str(tmp_path / "env" / "bin" / "python"), "-I", "-c"] and fe.TOKENIZER_REVISION in cmd[3] and "datamol-io/safe-gpt" in cmd[3]
    assert kw["env"]["HF_HOME"] == str(tmp_path / "hf_cache") and "HF_HUB_OFFLINE" not in kw["env"]            # online for this one step


# the next two use the real clone: the first only asks git (nothing is changed); the second works on a throw-away copy of two files
@pytest.mark.skipif(NO_FRIGID is not None, reason="FRIGID not installed")
def test_device_patch_is_recognised_as_applied_on_the_installed_clone():
    assert su.apply_mps_patch(fe.status()["clone"]) == "already applied"
    assert fe.mps_patch_applied(fe.status()["clone"])


@pytest.mark.skipif(NO_FRIGID is not None, reason="FRIGID not installed")
def test_device_patch_applies_to_the_pinned_files_once(tmp_path):
    clone = fe.status()["clone"]
    repo = tmp_path / "repo"
    for rel in ("src/dlm/sampler.py", "src/dlm/utils/spec2mol.py"):
        text = subprocess.run(["git", "-C", str(clone), "show", f"{fe.FRIGID_COMMIT}:{rel}"], capture_output=True, text=True, encoding="utf-8", check=True).stdout
        (repo / rel).parent.mkdir(parents=True, exist_ok=True)
        (repo / rel).write_text(text, encoding="utf-8")
    subprocess.run(["git", "-C", str(repo), "init", "-q"], check=True)
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t", "commit", "-qm", "pinned"], check=True)
    assert not fe.mps_patch_applied(repo)
    assert su.apply_mps_patch(repo) == "applied" and fe.mps_patch_applied(repo)
    assert su.apply_mps_patch(repo) == "already applied"


def test_a_clone_at_another_commit_is_never_moved(tmp_path):
    clone = tmp_path / "FRIGID"
    clone.mkdir()
    subprocess.run(["git", "-C", str(clone), "init", "-q"], check=True)
    subprocess.run(["git", "-C", str(clone), "-c", "user.name=t", "-c", "user.email=t@t", "commit", "-q", "--allow-empty", "-m", "x"], check=True)
    with pytest.raises(su.SetupError, match="mzDesk is verified with 9a3e623"):
        su.clone_frigid(tmp_path)


# ---------------------------------------------------------------- the real thing (opt-in)
def _memory_problem():
    from mzdesk.shared import memguard
    try:
        m = memguard.system_memory()
    except Exception as exc:  # noqa: BLE001
        return f"cannot read memory numbers: {exc}"
    if m["memory_free_pct"] < 20 or m["swap_free_mb"] < 500:
        return f"memory pressure: {m['memory_free_pct']} % free, swap free {m['swap_free_mb']} MB"
    return None


def _real_run_allowed():
    if os.environ.get("MZDESK_RUN_FRIGID") != "1":
        pytest.skip("real FRIGID runs are opt-in: set MZDESK_RUN_FRIGID=1 (about 20-40 s and 4-5 GB each)")
    if NO_FRIGID is not None:
        pytest.skip(f"FRIGID not installed: {NO_FRIGID}")
    problem = _memory_problem()
    if problem:
        pytest.skip(problem)
    assert not heavy_lock.LOCK.locked()


def _report(label, body):
    out = body["outputs"]
    print(f"\n--- {label}")
    print("verdict:", body["verdict"])
    for k, v in out["summary"]["values"].items():
        print(f"  {k}: {v if isinstance(v, str) else (v['value'], v.get('note'))}")
    for row in out["candidates"]["rows"][:5]:
        print("  row:", [c if isinstance(c, str) else "<png>" for c in row])
    print("warnings:", body["warnings"])


def test_real_run_example_through_the_api(monkeypatch):
    """One real run of the page example (caffeine, MassBank MSBNK-LCSB-LU119905, CC BY 4.0), formula given: device auto, Quick, the shipped length model."""
    _real_run_allowed()
    seen = {}
    real_build, real_predict = fe.build_job, fe.predict
    monkeypatch.setattr(fe, "build_job", lambda *a, **k: seen.setdefault("job", real_build(*a, **k)))
    monkeypatch.setattr(fe, "predict", lambda *a, **k: seen.setdefault("res", real_predict(*a, **k)))
    t0 = time.monotonic()
    r = _run(_module()["example"])
    wall = time.monotonic() - t0
    assert r.status_code == 200, r.text
    body = r.json()
    _report(f"example, {wall:.1f} s wall", body)
    res = seen["res"]
    print("engine:", {k: res[k] for k in ("device", "device_requested", "seconds", "driver_seconds", "peak_footprint_mb", "peak_rss_mb", "n_samples", "n_valid",
                                          "n_formula_matched", "n_unique_matched", "n_peaks_assigned", "length_model", "dropped")})
    assert seen["job"]["length_model"] == str(fe.LENGTH_MODEL) and res["length_model"] is True            # the joblib in the repo was loaded by the driver
    rows = body["outputs"]["candidates"]["rows"]
    assert rows and all(row[3] == "C8H10N4O2" for row in rows)
    scores = [float(row[4]) for row in rows]
    assert scores == sorted(scores, reverse=True) and 0 < scores[0] <= 1
    assert res["device"] in ("mps", "cpu") and res["n_samples"] == 16 and res["n_peaks_assigned"] > 0
    assert res["seconds"] < fe.TIMEOUT_S["quick"] and 0 < res["peak_footprint_mb"] < fe.FOOTPRINT_CAP_MB
    assert body["verdict"].startswith("Top generated candidate:") and "(hypotheses)" in body["verdict"]
    assert not heavy_lock.LOCK.locked()


def test_real_run_without_a_formula_takes_buddys_top_formula(monkeypatch):
    """One real run with the formula left blank: BUDDY (real, in this process) gives the formula, FRIGID (real) generates for it."""
    _real_run_allowed()
    from mzdesk.ms1 import msbuddy_engine as bd
    ms2 = [tuple(map(float, ln.split())) for ln in CAFFEINE_MSFINDER_MS2.splitlines()]
    expected = bd.annotate_formula(195.0877, True, "[M+H]+", None, "", tuple(ms2), (), "orbitrap", False)["candidates"]
    assert expected, "BUDDY found no formula for the example spectrum"
    t0 = time.monotonic()
    r = _run(_no_formula())
    wall = time.monotonic() - t0
    assert r.status_code == 200, r.text
    body = r.json()
    _report(f"BUDDY fallback, {wall:.1f} s wall", body)
    used = body["outputs"]["summary"]["values"]["Formula used"]
    assert used["value"] == fe.check_formula(expected[0]["formula"]) and used["note"].startswith("BUDDY top-1")
    assert any("BUDDY's top-1 formula" in w for w in body["warnings"])
    assert all(row[3] == used["value"] for row in body["outputs"]["candidates"]["rows"])


def test_buddy_fallback_gets_ms1_and_the_halogen_switch(monkeypatch, fake_predict):
    calls = _buddy(monkeypatch, [{"rank": 1, "formula": "C8H10N4O2", "prob": 0.9}])
    body = _run({**_no_formula(), "ms1": "195.0877 100\n196.0902 10", "halogen": True}).json()
    assert calls[0]["ms1"] == ((195.0877, 100.0), (196.0902, 10.0)) and calls[0]["halogen"] is True
    assert "no halogens" not in " ".join(body["warnings"])
