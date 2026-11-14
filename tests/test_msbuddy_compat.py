"""msbuddy 0.3 and 1.0 give the same formula counts through msbuddy_engine.read_formula."""
import numpy as np
import pytest

pytest.importorskip("msbuddy")
from mzdesk.ms1 import msbuddy_engine as bd

CAFFEINE = np.array([8, 10, 0, 0, 0, 0, 0, 4, 0, 2, 0, 0], dtype=np.int16)   # C8H10N4O2 in msbuddy's order


def test_the_installed_msbuddy_reads_a_formula_as_the_array():
    assert np.array_equal(bd.read_formula("C8H10N4O2"), CAFFEINE)


@pytest.mark.parametrize("bad", ["C8H10N4O2Xx", "C8(H", "", "**bold**"])
def test_the_installed_msbuddy_gives_none_or_no_counts_for_text_that_is_no_formula(bad):
    assert bd._safe_read_formula(bad) is None


def test_msbuddy_1_dict_becomes_the_array(monkeypatch):
    monkeypatch.setattr(bd, "_buddy_read_formula", lambda s: {"C": 8, "H": 10, "N": 4, "O": 2})
    out = bd.read_formula("C8H10N4O2")
    assert out.dtype == np.int16 and np.array_equal(out, CAFFEINE)


def test_msbuddy_1_error_becomes_none(monkeypatch):
    def raises(s):
        raise ValueError("Unsupported element 'Xx'.")
    monkeypatch.setattr(bd, "_buddy_read_formula", raises)
    assert bd.read_formula("Xx") is None


def test_msbuddy_0_3_array_and_none_pass_through(monkeypatch):
    monkeypatch.setattr(bd, "_buddy_read_formula", lambda s: CAFFEINE.copy())
    assert np.array_equal(bd.read_formula("C8H10N4O2"), CAFFEINE)
    monkeypatch.setattr(bd, "_buddy_read_formula", lambda s: None)
    assert bd.read_formula("Xx") is None


def test_ion_mz_is_the_same_whichever_form_msbuddy_returns(monkeypatch):
    mz_installed = bd.ion_mz_and_formula("C8H10N4O2", "[M+H]+")
    monkeypatch.setattr(bd, "_buddy_read_formula", lambda s: {"C": 8, "H": 10, "N": 4, "O": 2})
    assert bd.ion_mz_and_formula("C8H10N4O2", "[M+H]+") == mz_installed
    assert mz_installed[1] == "C8H11N4O2+" and abs(mz_installed[0] - 195.08765) < 1e-4
