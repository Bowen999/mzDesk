"""MS2 清洗、格式写出、基础可视化的测试。"""
import pytest

from mzdesk.shared.ms2_clean import NORMALIZE, clean_spectrum
from mzdesk.shared.ms2_formats import FORMATS, all_formats, to_format
from mzdesk.shared.ms2_plot import spectrum_figure
from mzdesk.shared.spectrum import first_record, parse_spectrum

PEAKS = [(42.0338, 5.0), (110.0713, 18.5), (138.0662, 100.0), (195.0877, 32.0)]


def _same(a, b, rel=1e-4):
    assert [round(m, 4) for m, _ in a] == [round(m, 4) for m, _ in b]
    assert [i for _, i in a] == pytest.approx([i for _, i in b], rel=rel)


# ---------------------------------------------------------------- writers
@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize("kw", [{}, {"precursor_mz": 195.0877, "name": "caffeine 2:1", "ion_mode": "Negative"}])
def test_every_writer_round_trips(fmt, kw):
    _same(parse_spectrum(to_format(PEAKS, fmt, **kw)), PEAKS)


def test_round_trip_for_raw_and_tiny_intensities():
    for scale in (1e6, 1e-3):
        peaks = [(m, i * scale + 0.123456789 * scale) for m, i in PEAKS]
        for fmt in FORMATS:
            _same(parse_spectrum(to_format(peaks, fmt)), peaks)


def test_all_formats_lists_every_format_in_order():
    out = all_formats(PEAKS, precursor_mz=195.0877)
    assert list(out) == FORMATS and all(out.values())


def test_massbank_relative_intensity_maxes_at_999():
    text = to_format(PEAKS, "MassBank", precursor_mz=195.0877, name="caffeine")
    rel = [int(ln.split()[2]) for ln in text.splitlines() if ln.startswith("  ")]
    assert max(rel) == 999 and min(rel) >= 1
    assert "PK$NUM_PEAK: 4" in text and "PRECURSOR_M/Z 195.0877" in text
    assert first_record(text)["precursor_mz"] == 195.0877


def test_record_headers():
    msp = to_format(PEAKS, "MSP", precursor_mz=195.0877, name="caffeine", ion_mode="Positive")
    assert msp.splitlines()[:4] == ["NAME: caffeine", "PRECURSORMZ: 195.0877", "IONMODE: Positive", "Num Peaks: 4"]
    mgf = to_format(PEAKS, "MGF", precursor_mz=195.0877, name="caffeine", ion_mode="Negative")
    assert "PEPMASS=195.0877" in mgf and "CHARGE=1-" in mgf and mgf.endswith("END IONS")
    assert ">parentmass 195.0877" in to_format(PEAKS, "SIRIUS .ms", precursor_mz=195.0877)


def test_unknown_format_raises():
    with pytest.raises(ValueError):
        to_format(PEAKS, "nope")


# ---------------------------------------------------------------- cleaning
RAW = [(50.0, 0.5), (100.0, 20.0), (150.0, 100.0), (195.0877, 40.0), (195.09, 3.0)]


def test_noise_threshold_is_percent_of_input_base_peak():
    res = clean_spectrum(RAW, noise_pct=1.0)
    assert [r["reason"] for r in res["rows"]] == ["noise", "", "", "", ""]
    assert res["n_in"] == 5 and res["n_out"] == 4 and (50.0, 0.5) not in res["peaks"]


def test_precursor_removal_uses_tolerance_and_needs_a_precursor():
    res = clean_spectrum(RAW, precursor_mz=195.0877, remove_precursor=True, precursor_tol=0.02)
    assert [r["reason"] for r in res["rows"]] == ["", "", "", "precursor", "precursor"]
    assert clean_spectrum(RAW, precursor_mz=195.0877, remove_precursor=True, precursor_tol=0.001)["n_out"] == 4
    assert clean_spectrum(RAW, remove_precursor=True)["n_out"] == 5
    assert clean_spectrum(RAW, precursor_mz=195.0877)["n_out"] == 5


def test_precursor_reason_wins_over_noise_and_kept_flag_matches_reason():
    res = clean_spectrum(RAW, noise_pct=50, precursor_mz=195.0877, remove_precursor=True)
    assert [r["reason"] for r in res["rows"]] == ["noise", "noise", "", "precursor", "precursor"]
    assert all(r["kept"] == (r["reason"] == "") for r in res["rows"])


def test_normalizations_apply_after_removal():
    kept = lambda **kw: [i for _, i in clean_spectrum(RAW, remove_precursor=True, precursor_mz=195.0877, **kw)["peaks"]]
    assert kept() == [0.5, 20.0, 100.0]
    assert kept(normalize="Base peak = 100") == pytest.approx([0.5, 20.0, 100.0])
    assert kept(normalize="Base peak = 999") == pytest.approx([0.5 * 9.99, 199.8, 999.0])
    assert kept(normalize="Base peak = 1") == pytest.approx([0.005, 0.2, 1.0])
    assert sum(kept(normalize="Sum = 1")) == pytest.approx(1.0)
    # base peak removed as precursor: the next peak becomes the base
    res = clean_spectrum([(100.0, 10.0), (195.0877, 50.0)], precursor_mz=195.0877, remove_precursor=True, normalize="Base peak = 100")
    assert res["peaks"] == [(100.0, 100.0)] and res["rows"][1]["rel"] == 100.0 and res["rows"][0]["rel"] == 20.0


def test_sort_orders():
    assert [m for m, _ in clean_spectrum(RAW, sort_by="m/z")["peaks"]] == sorted(m for m, _ in RAW)
    assert [i for _, i in clean_spectrum(RAW, sort_by="Intensity")["peaks"]] == [100.0, 40.0, 20.0, 3.0, 0.5]
    assert [r["mz"] for r in clean_spectrum(RAW, sort_by="Intensity")["rows"]] == sorted(m for m, _ in RAW)


def test_invalid_peaks_dropped_duplicates_collapsed_and_bad_options_rejected():
    res = clean_spectrum([(100.0, 5.0), (100.0, 9.0), (0.0, 5.0), (120.0, -1.0), (130.0, float("nan"))])
    assert res["peaks"] == [(100.0, 9.0)] and res["n_in"] == 1
    assert clean_spectrum([])["peaks"] == []
    with pytest.raises(ValueError):
        clean_spectrum(RAW, normalize="x")
    with pytest.raises(ValueError):
        clean_spectrum(RAW, sort_by="x")
    assert "None" in NORMALIZE


# ---------------------------------------------------------------- figure
def test_spectrum_figure_has_peaks_trace_and_one_card_per_peak():
    fig, cards = spectrum_figure(PEAKS, precursor_mz=195.0877, dim=[0])
    trace = [t for t in fig.data if t.name == "peaks"]
    assert len(trace) == 1 and list(trace[0].customdata) == list(range(len(PEAKS)))
    assert len(cards) == len(PEAKS)
    assert [c["on"] for c in cards] == [False, True, True, True]
    assert cards[0]["kind"] == "Removed" and cards[3]["kind"] == "Precursor ion" and cards[1]["kind"] == "Peak"
    assert cards[2]["rows"] == [["Intensity", "100"], ["Rel. intensity", "100.0 %"]]
    assert fig.layout.yaxis.range[0] == 0 and fig.layout.yaxis.range[1] > 100
    assert any(s.line.dash == "dot" for s in fig.layout.shapes)


def test_spectrum_figure_dim_forms_and_reasons():
    _, by_mask = spectrum_figure(PEAKS, dim=[True, False, False, False])
    _, by_reason = spectrum_figure(PEAKS, dim={3: "precursor"}, precursor_mz=195.0877)
    assert by_mask[0]["on"] is False and by_mask[1]["on"] is True
    assert by_reason[3]["kind"] == "Removed: precursor"


def test_spectrum_figure_raw_intensities():
    raw = [(m, i * 1e4) for m, i in PEAKS]
    fig, cards = spectrum_figure(raw)
    assert fig.layout.yaxis.tickformat == "~s" and len(cards) == 4
    assert cards[2]["rows"][0] == ["Intensity", "1,000,000"]


def test_spectrum_figure_empty():
    fig, cards = spectrum_figure([])
    assert cards == [] and len(fig.data) == 0
