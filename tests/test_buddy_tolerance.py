"""BUDDY must search the precursor window the user asked for (real msbuddy; no network).

msbuddy's presets cover at most +-10 ppm (qtof). With a wider tolerance the ESI MS/MS Annotation used to fall back to qtof and the table said
"within 20 ppm" although BUDDY had only looked within 10 ppm, so a formula between 10 and 20 ppm off never appeared.
"""
import pytest

pytest.importorskip("msbuddy")

from mzdesk import spec                                    # noqa: E402
from mzdesk.annotation import formula_step as F            # noqa: E402
from mzdesk.ms1 import msbuddy_engine as bd                # noqa: E402
from mzdesk.shared.spectrum import parse_spectrum          # noqa: E402

CAFFEINE = "C8H10N4O2"
MS2 = parse_spectrum(spec.CAFFEINE_MASSBANK_MS2)
THEORETICAL = bd.ion_mz_and_formula(CAFFEINE, "[M+H]+")[0]
PRECURSOR_PLUS_14 = THEORETICAL * (1 + 14e-6)               # caffeine measured 14 ppm high


def _formulas(tol_ppm):
    _, result = F.run_adduct_and_formula(PRECURSOR_PLUS_14, True, MS2, adduct="[M+H]+", tol_ppm=tol_ppm)
    return result


def test_a_formula_14_ppm_off_is_found_with_a_20_ppm_tolerance():
    result = _formulas(20.0)
    assert result.status == "done" and result.formulas
    assert CAFFEINE in result.formula_strings and result.formulas[0].formula == CAFFEINE
    assert result.formulas[0].ppm == pytest.approx(14.0, abs=0.1)


def test_a_formula_14_ppm_off_stays_out_of_the_default_10_ppm_table():
    result = _formulas(10.0)
    assert CAFFEINE not in result.formula_strings
    assert all(abs(f.ppm) <= 10.0 for f in result.formulas)                       # the user's window still filters the list


@pytest.mark.parametrize("tol_ppm", [30.0, 50.0])
def test_the_wider_the_tolerance_the_more_formulas_fit(tol_ppm):
    assert CAFFEINE in _formulas(tol_ppm).formula_strings


def test_the_engine_cache_tells_tolerances_apart():
    """Same query, three searches: a cache key without the tolerance would return the first answer for all of them."""
    def top(ms_instr):
        return [c["formula"] for c in bd.annotate_formula(PRECURSOR_PLUS_14, True, "[M+H]+", None, "", tuple(MS2), (), ms_instr, False)["candidates"]]

    assert CAFFEINE not in top("qtof")
    assert CAFFEINE in top((20.0, 20.0))
    assert CAFFEINE not in top((10.0, 20.0))                                      # qtof's own numbers given as a pair: the same search as the preset
    assert CAFFEINE in top((30.0, 20.0))
    assert top("qtof") == top((10.0, 20.0))


def test_a_list_pair_is_accepted_like_a_tuple():
    cands = bd.annotate_formula(PRECURSOR_PLUS_14, True, "[M+H]+", None, "", tuple(MS2), (), [20, 20], False)["candidates"]
    assert CAFFEINE in [c["formula"] for c in cands]
