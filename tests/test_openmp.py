"""mzdesk.shared.openmp: LightGBM (BUDDY's formula model) on Macs without libomp."""
import pytest

from mzdesk.shared import openmp


@pytest.fixture(autouse=True)
def fresh(monkeypatch, tmp_path):
    monkeypatch.setitem(openmp._state, "ok", False)
    monkeypatch.setattr(openmp.paths, "tools_dir", lambda: tmp_path / "tools")


def failing_then(results):
    calls = []

    def fake():
        calls.append(1)
        outcome = results[len(calls) - 1]
        if outcome:
            raise outcome
    return fake, calls


def test_a_lightgbm_that_loads_is_left_alone(monkeypatch):
    loaded = []
    monkeypatch.setattr(openmp, "_import_lightgbm", lambda: None)
    monkeypatch.setattr(openmp.ctypes, "CDLL", lambda *a, **k: loaded.append(a))
    openmp.ensure_lightgbm()
    assert loaded == [] and openmp._state["ok"]


def test_on_a_mac_without_libomp_mzdesks_copy_is_loaded_first(monkeypatch):
    fake, calls = failing_then([OSError("Library not loaded: @rpath/libomp.dylib"), None])
    loaded = []
    monkeypatch.setattr(openmp, "_import_lightgbm", fake)
    monkeypatch.setattr(openmp.sys, "platform", "darwin")
    monkeypatch.setattr(openmp.ctypes, "CDLL", lambda path, mode=0: loaded.append((path, mode)))
    openmp.dylib().parent.mkdir(parents=True)
    openmp.dylib().write_bytes(b"")
    openmp.ensure_lightgbm()
    assert len(calls) == 2 and loaded == [(str(openmp.dylib()), openmp.ctypes.RTLD_GLOBAL)]


@pytest.mark.parametrize("platform", ["darwin", "linux"])
def test_without_a_copy_the_error_says_how_to_add_it(monkeypatch, platform):
    fake, _ = failing_then([OSError("Library not loaded: @rpath/libomp.dylib")])
    monkeypatch.setattr(openmp, "_import_lightgbm", fake)
    monkeypatch.setattr(openmp.sys, "platform", platform)
    with pytest.raises(RuntimeError, match=r"mzdesk\.shared\.openmp --install"):
        openmp.ensure_lightgbm()
    assert not openmp._state["ok"]


def test_install_is_for_macos_only(monkeypatch):
    monkeypatch.setattr(openmp.sys, "platform", "linux")
    assert openmp.install() is None


def test_buddy_checks_lightgbm_before_it_runs(monkeypatch):
    from mzdesk.ms1 import msbuddy_engine as bd
    def missing():
        raise RuntimeError(openmp.MISSING)
    monkeypatch.setattr(bd.openmp, "ensure_lightgbm", missing)
    bd._annotate_formula.cache_clear()
    with pytest.raises(RuntimeError, match="OpenMP"):
        bd.annotate_formula(195.0877, True, "[M+H]+", None, "", ((138.0662, 100.0),), ())
    bd._annotate_formula.cache_clear()
