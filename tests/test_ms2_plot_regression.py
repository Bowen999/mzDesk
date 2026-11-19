"""Regression tests for display only; peak processing and matching stay unchanged."""
import json

import pytest

from mzdesk.api.common import plotly
from mzdesk.api.runners.ms2 import run_clean
from mzdesk.ms2.mirror_plot import build_mirror_figure, match_summary
from mzdesk.shared.ms2_clean import clean_spectrum
from mzdesk.shared.ms2_plot import fmt_mz, spectrum_figure


def _trace(fig, name):
    return next(trace for trace in fig.data if trace.name == name)


def _halves(fig):
    """The two hover traces of a mirror plot (both named "peaks"): query first, then reference."""
    return [trace for trace in fig.data if trace.name == "peaks"]


def _row(card, label):
    return next(value for name, value in card["rows"] if name == label)


def test_mirror_keeps_true_small_stick_heights_and_hover_intensities():
    a = [(100, 100), (120, 0.01), (140, 0.02)]
    b = [(100.01, 100), (120.01, 0.01), (160, 0.02)]
    summary_before = match_summary(a, b, tolerance=0.02)
    fig, cards = build_mirror_figure(a, b, tolerance=0.02)
    assert list(_trace(fig, "Query · matched").y) == pytest.approx([0, 100, None, 0, 0.01, None])
    assert list(_trace(fig, "Reference · matched").y) == pytest.approx([0, -100, None, 0, -0.01, None])
    assert list(_trace(fig, "Query · unmatched").y) == pytest.approx([0, 0.02, None])
    q, r = _halves(fig)
    assert list(q.y) == pytest.approx([100, 0.01, 0.02])
    assert list(r.y) == pytest.approx([-100, -0.01, -0.02])
    assert list(q.customdata) == [0, 1, 2] and list(r.customdata) == [3, 4, 5]
    assert [c["rel"] for c in cards] == pytest.approx([100, 0.01, 0.02, -100, -0.01, -0.02])
    assert _row(cards[1], "Rel. intensity") == "0.01 %"
    assert _row(cards[4], "Rel. intensity") == "0.01 %"
    assert [c["on"] for c in cards] == [True, True, False, True, True, False]
    assert cards[2]["kind"] == "Query · unmatched" and cards[3]["kind"] == "Reference · matched"
    assert summary_before == match_summary(a, b, tolerance=0.02)
    assert summary_before["matched_pairs"] == 2
    json.dumps(plotly(fig, peaks=cards), allow_nan=False)


def test_a_low_mirror_plot_gets_short_y_axis_titles():
    """Each half of the 370 px plot on the workflow page is about 150 px: the long title would run into the other half's."""
    a, b = [(100, 100)], [(100, 90)]
    titles = lambda fig: [ann.text for ann in fig.layout.annotations if "rel." in (ann.text or "")]
    assert titles(build_mirror_figure(a, b, labels=("Experimental", "Library"))[0]) == [
        "Experimental · rel. intensity (%)", "Library · rel. intensity (%)"]
    assert titles(build_mirror_figure(a, b, labels=("Experimental", "Library"), height=370)[0]) == [
        "Experimental · rel. %", "Library · rel. %"]


def test_mirror_hover_and_labels_distinguish_near_masses():
    a, b = [(100.12341, 100)], [(100.12344, 100)]
    fig, cards = build_mirror_figure(a, b, precursor_mz1=200.123456, precursor_mz2=200.1234565)
    q, r = cards
    assert q["mz"] == 100.12341 and r["mz"] == 100.12344
    assert _row(q, "m/z (more precision)") == "100.12341"
    assert _row(r, "m/z (more precision)") == "100.12344"
    assert _row(q, "Reference m/z") == "100.12344" and _row(q, "Δ m/z") == "+3e-05 Da"
    assert _row(r, "Query m/z") == "100.12341" and _row(r, "Δ m/z") == "-3e-05 Da"
    assert "Precursor − peak" in dict(q["rows"])
    assert "Neutral loss" not in dict(q["rows"])
    labels = [ann.text for ann in fig.layout.annotations]
    # distinct precursors that round to the same four decimals keep their full digits
    assert "Query precursor 200.123456" in labels and "Reference precursor 200.1234565" in labels
    assert len(fig.layout.shapes) == 2  # distinct precursor values are not merged


def test_mirror_distinguishes_sub_eight_decimal_peaks_and_nonzero_differences():
    a, b = 100.123456781, 100.123456782
    fig, cards = build_mirror_figure([(a, 100)], [(b, 100)], precursor_mz1=b, precursor_mz2=a)
    q, r = cards
    assert _row(q, "m/z (more precision)") == "100.123456781"
    assert _row(r, "m/z (more precision)") == "100.123456782"
    assert _row(q, "Reference m/z") == "100.123456782"
    assert _row(r, "Query m/z") == "100.123456781"
    assert float(_row(q, "Precursor − peak")) > 0
    assert _row(q, "Δ m/z").startswith("+") and _row(r, "Δ m/z").startswith("-")
    assert len(fig.layout.shapes) == 2


def test_mirror_precursor_label_has_four_decimals():
    fig, _ = build_mirror_figure([(100, 100)], [(100, 100)], precursor_mz1=195.08765203118054, precursor_mz2=195.08765203118054)
    assert "precursor 195.0877" in [ann.text for ann in fig.layout.annotations]


@pytest.mark.parametrize("a,b", [([], []), ([(100, 0)], []), ([], [(100, 10)]), ([(100, 10)], [])])
def test_mirror_empty_zero_and_one_sided_spectra(a, b):
    fig, cards = build_mirror_figure(a, b)
    has_peaks = any(intensity > 0 for _, intensity in a + b)
    assert bool(fig.data) == has_peaks
    assert len(cards) == sum(intensity > 0 for _, intensity in a + b)
    assert match_summary(a, b)["matched_pairs"] == 0
    json.dumps(plotly(fig, peaks=cards), allow_nan=False)


def test_single_precursor_label_and_percent_ticks():
    fig, _ = spectrum_figure([(110.0713, 1), (195.0876, 100)], precursor_mz=195.08765203118054)
    assert [ann.text for ann in fig.layout.annotations] == ["precursor 195.0877"]
    assert not fig.layout.shapes  # a tall precursor peak leaves no stub of dotted line under the top
    assert list(fig.layout.yaxis.tickvals) == [0, 20, 40, 60, 80, 100]


def test_dense_spectra_get_thinner_sticks():
    sparse, _ = spectrum_figure([(100, 100), (150, 50)])
    dense, _ = spectrum_figure([(float(m), 10 + m % 7) for m in range(30, 400)])
    assert _trace(sparse, "kept").line.width == 2
    assert _trace(dense, "kept").line.width == 1


def test_single_zero_and_empty_peaks_keep_cards_and_finite_axis():
    fig, cards = spectrum_figure([])
    assert not fig.data and cards == []
    fig, cards = spectrum_figure([(100, 0), (120, 0)], precursor_mz=120, dim=[0])
    assert fig.layout.yaxis.range == (0, 1)
    assert len(cards) == 2 and not cards[0]["on"]
    assert all(card["rows"][1] == ["Rel. intensity", "—"] for card in cards)
    assert list(_trace(fig, "peaks").customdata) == [0, 1]
    json.dumps(plotly(fig, peaks=cards), allow_nan=False)


@pytest.mark.parametrize("peaks", [[(100, -1)], [(100, float("nan"))], [(float("inf"), 1)], [(0, 1)]])
def test_single_rejects_invalid_peaks_without_changing_card_indices(peaks):
    with pytest.raises(ValueError):
        spectrum_figure(peaks)


def test_single_near_precursor_is_not_asserted_to_be_a_precursor_ion():
    peaks = [(100, 100), (195.04, 50)]
    _, nearby = spectrum_figure(peaks, precursor_mz=195, precursor_tol=0.05)
    _, outside = spectrum_figure(peaks, precursor_mz=195, precursor_tol=0.001)
    assert nearby[1]["kind"] == "Near precursor m/z"
    assert outside[1]["kind"] == "Peak"
    assert nearby[1]["mz"] == outside[1]["mz"] == 195.04


def test_single_exposes_extra_mz_digits_and_small_positive_percentage():
    _, cards = spectrum_figure([(100.12341, 100), (100.12344, 0.01)])
    assert cards[0]["rows"][-1] == ["m/z (more precision)", "100.12341"]
    assert cards[1]["rows"][-1] == ["m/z (more precision)", "100.12344"]
    assert cards[1]["rows"][1] == ["Rel. intensity", "0.01 %"]
    assert fmt_mz(100) == "100.0000" and fmt_mz(100.12345678) == "100.12345678"


def test_single_exposes_very_close_masses_and_small_scientific_notation():
    a, b = 100.123456781, 100.123456782
    _, cards = spectrum_figure([(a, 100), (b, 0.01)])
    assert cards[0]["rows"][-1] == ["m/z (more precision)", "100.123456781"]
    assert cards[1]["rows"][-1] == ["m/z (more precision)", "100.123456782"]
    assert float(fmt_mz(1e-12)) == 1e-12


def test_cleanup_explains_both_bases_and_does_not_change_output():
    inputs = {"spectrum": {"source": "Paste", "value": "100.12341 10\n195 50"},
              "precursor_mz": "195", "remove_precursor": True, "normalize": "Base peak = 100"}
    result = run_clean(inputs)
    out = result["outputs"]
    figure = out["spectrum"]
    assert figure["plotly"]["layout"]["yaxis"]["title"]["text"] == "Intensity (cleanup scale)"
    assert figure["plotly"]["layout"]["yaxis"]["range"] == [0, 575]
    assert figure["peaks"][0]["rows"][1] == ["Rel. int. (input base)", "20.0 %"]
    assert figure["peaks"][1]["rows"][1] == ["Rel. int. (input base)", "100.0 %"]
    assert "input base peak before removal" in figure["help"][1]
    assert out["peaks"]["rows"][0][0] == "100.12341"
    expected = clean_spectrum([(100.12341, 10), (195, 50)], precursor_mz=195,
                              remove_precursor=True, normalize="Base peak = 100")
    assert expected["peaks"] == [(100.12341, 100)]
    assert out["formats"]["formats"]["Peak list"] == "100.12341\t100"
    json.dumps(result, allow_nan=False)


def test_cleanup_uses_selected_precursor_tolerance_and_handles_all_removed():
    common = {"spectrum": {"source": "Paste", "value": "100 100\n195.04 50"},
              "precursor_mz": "195", "precursor_tol": "0.001", "remove_precursor": True}
    kept = run_clean(common)["outputs"]["spectrum"]["peaks"]
    assert kept[1]["kind"] == "Peak" and kept[1]["on"] is True
    removed = run_clean({**common, "spectrum": {"source": "Paste", "value": "195 100\n195.04 50"},
                         "precursor_tol": "0.05", "noise_pct": "100",
                         "normalize": "Base peak = 100"})
    assert "Every peak was removed." in removed["warnings"]
    assert "no kept-peak normalization reference" in removed["outputs"]["spectrum"]["help"][-1]
    assert all(not c["on"] for c in removed["outputs"]["spectrum"]["peaks"])


@pytest.mark.parametrize("normalize", ["None", "Base peak = 100", "Base peak = 999", "Base peak = 1", "Sum = 1"])
def test_cleanup_all_normalizations_keep_axis_and_input_base_explicit(normalize):
    inputs = {"spectrum": {"source": "Paste", "value": "100 10\n150 5\n195 50"},
              "precursor_mz": "195", "remove_precursor": True, "normalize": normalize}
    figure = run_clean(inputs)["outputs"]["spectrum"]
    expected = clean_spectrum([(100, 10), (150, 5), (195, 50)], precursor_mz=195,
                              remove_precursor=True, normalize=normalize)
    tips = next(t for t in figure["plotly"]["data"] if t["name"] == "peaks")
    assert tips["y"] == pytest.approx([r["intensity"] for r in expected["rows"]])
    assert [c["rows"][1] for c in figure["peaks"]] == [
        ["Rel. int. (input base)", "20.0 %"], ["Rel. int. (input base)", "10.0 %"],
        ["Rel. int. (input base)", "100.0 %"],
    ]
    if normalize == "None":
        assert figure["plotly"]["layout"]["yaxis"]["title"]["text"] == "Intensity"
        assert "no normalization was requested" in figure["help"][0]
    else:
        assert figure["plotly"]["layout"]["yaxis"]["title"]["text"] == "Intensity (cleanup scale)"
        assert normalize in figure["help"][0]
        assert "kept peaks" in figure["help"][0]


def test_cleanup_without_normalization_does_not_claim_missing_normalization_reference():
    result = run_clean({"spectrum": {"source": "Paste", "value": "195 50"},
                        "precursor_mz": "195", "remove_precursor": True, "normalize": "None"})
    help_text = " ".join(result["outputs"]["spectrum"]["help"])
    assert "no normalization was requested" in help_text
    assert "normalization reference" not in help_text
