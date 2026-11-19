"""EI marker ions of silyl derivatives (mzdesk/ms2/ei_markers.py): the thresholds, the TMS count bound, the Auto inference and the checks of a chosen derivatization.
All spectra are synthetic; the rules are the ones in the module docstring."""
import pytest

from mzdesk.ms2 import ei_markers as em


def spectrum(**ions):
    """{"mz73": 100, ...} style helper: keyword names like mz73 give the peak list [(73.0, 100.0), ...]."""
    return sorted((float(k[2:]), float(v)) for k, v in ions.items())


def marker(analysis, marker_id):
    return next(m for m in analysis["markers"] if m["id"] == marker_id)


def test_the_constants_are_the_documented_ones():
    assert (em.MIN_73, em.MIN_75, em.MIN_147, em.MIN_147_BLIND, em.MIN_M_ION, em.MIN_M57_TBDMS) == (10.0, 1.0, 2.0, 5.0, 1.0, 20.0)
    assert (em.MIN_M31, em.M31_OF_M15, em.MIN_TBDMS_SI, em.TBDMS_SI_IONS) == (2.0, 0.1, 2.0, (73, 75, 115))
    assert (em.TMS_NOMINAL, em.TBDMS_NOMINAL, em.MEOX_NOMINAL) == (72, 114, 29)       # C3H8Si, C6H14Si, CH3N
    assert [name for name, _ in em.LOSSES] == ["M-15", "M-31", "M-57", "M-90", "M-180"]


def test_relative_intensity_is_a_percentage_of_the_base_peak_within_half_a_dalton():
    peaks = spectrum(mz73=500, mz147=50, mz200=1000)
    assert em.relative_intensity(peaks, 73) == pytest.approx(50.0)
    assert em.relative_intensity(peaks, 73.4) == pytest.approx(50.0) and em.relative_intensity(peaks, 73.6) == 0.0
    assert em.relative_intensity(peaks, 147) == pytest.approx(5.0)
    assert em.relative_intensity([], 73) == 0.0 and em.relative_intensity([(73.0, 0.0)], 73) == 0.0


def test_a_tms_spectrum_without_the_molecular_ion():
    a = em.analyze(spectrum(mz73=100, mz75=20, mz147=30, mz174=75))
    assert a["tms"] == "yes" and a["tms_min"] == 2 and a["derivatization"] == "TMS" and a["meox"] == "undecided" and not a["tbdms"]
    assert marker(a, "73")["present"] is True and marker(a, "75")["present"] is True and marker(a, "147")["present"] is True
    assert all(marker(a, k)["present"] is None and marker(a, k)["rel"] is None for k in ("M", "M-15", "M-31", "M-57", "M-90", "M-180"))
    assert a["m_check"] is None and a["warnings"] == []
    assert [m["id"] for m in a["markers"]] == ["73", "75", "147", "115", "M", "M-15", "M-31", "M-57", "M-90", "M-180"]


def test_one_tms_group_gives_a_lower_bound_of_one():
    a = em.analyze(spectrum(mz73=100, mz75=20, mz117=80, mz147=1))           # 147 under 2 %
    assert a["tms"] == "yes" and a["tms_min"] == 1 and marker(a, "147")["present"] is False


def test_the_thresholds_are_inclusive():
    assert marker(em.analyze(spectrum(mz73=10, mz200=100)), "73")["present"] is True
    assert marker(em.analyze(spectrum(mz73=9.9, mz200=100)), "73")["present"] is False
    assert marker(em.analyze(spectrum(mz73=100, mz147=2)), "147")["present"] is True
    assert marker(em.analyze(spectrum(mz73=100, mz147=1.9)), "147")["present"] is False
    assert marker(em.analyze(spectrum(mz73=100, mz75=1)), "75")["present"] is True and marker(em.analyze(spectrum(mz73=100, mz75=0.9)), "75")["present"] is False
    assert em.analyze(spectrum(mz73=100, mz75=1))["tms"] == "yes" and em.analyze(spectrum(mz73=100, mz75=0.9))["tms"] == "no"
    assert em.analyze(spectrum(mz73=100, mz147=2))["tms"] == "yes" and em.analyze(spectrum(mz73=100, mz147=1.9))["tms"] == "no"
    a = em.analyze(spectrum(mz73=100, mz276=1, mz291=1), 291)
    assert marker(a, "M")["present"] is True and marker(a, "M-15")["present"] is True
    assert marker(em.analyze(spectrum(mz73=100, mz276=0.9), 291), "M-15")["present"] is False


def test_the_molecular_ion_and_its_losses_are_read_at_m_minus_the_loss():
    peaks = spectrum(mz73=100, mz147=20, mz201=4, mz276=8, mz291=2)           # M = 291: M-90 = 201, M-15 = 276
    a = em.analyze(peaks, 291)
    assert marker(a, "M-15")["mz"] == 276 and marker(a, "M-15")["rel"] == pytest.approx(8.0) and marker(a, "M-90")["present"] is True
    assert marker(a, "M-180")["present"] is False and marker(a, "M-31")["present"] is False and marker(a, "M")["rel"] == pytest.approx(2.0)
    assert a["m_check"] == "ok" and a["warnings"] == []


def test_tms_oh_losses_are_shown_but_decide_nothing():
    """M-90 and M-180 do not separate TMS from underivatized spectra in MassBank (module docstring), so they neither raise the count bound nor corroborate."""
    a = em.analyze(spectrum(mz73=100, mz75=20, mz147=1, mz210=3, mz300=3), 390)       # M-90 = 300, M-180 = 210
    assert marker(a, "M-90")["present"] is True and marker(a, "M-180")["present"] is True
    assert a["tms"] == "yes" and a["tms_min"] == 1 and not any("M-90" in r or "M-180" in r for r in a["reasons"])
    alone = em.analyze(spectrum(mz100=100, mz210=3, mz300=3), 390)           # losses of 90 and 180 without m/z 73 or 147: an underivatized spectrum can do that
    assert alone["tms"] == "no" and alone["derivatization"] == "None" and alone["tms_min"] == 0
    assert em.analyze(spectrum(mz73=100, mz210=3, mz300=3), 390)["tms"] == "no"       # even with a strong 73 they are no backing


def test_m_minus_31_makes_it_meox_tms_and_no_m_minus_31_leaves_it_undecided():
    base = {"mz73": 100, "mz147": 25, "mz276": 5}
    assert em.analyze(spectrum(**base), 291)["meox"] == "undecided"
    assert em.analyze(spectrum(**base), 291)["derivatization"] == "TMS"
    with_31 = em.analyze(spectrum(**base, mz260=3), 291)
    assert with_31["meox"] == "yes" and with_31["derivatization"] == "MeOx + TMS" and marker(with_31, "M-31")["present"] is True
    assert em.analyze(spectrum(**base, mz260=3))["meox"] == "undecided"       # without M there is no M-31


def test_m_minus_31_must_be_clear_to_count():
    """MeOx call: M-31 >= 2 % and at least a tenth of M-15 (a lone 1 % M-31 is in many TMS-only spectra, module docstring)."""
    base = {"mz73": 100, "mz147": 25, "mz276": 20}
    assert em.analyze(spectrum(**base, mz260=1.9), 291)["meox"] == "undecided"
    assert em.analyze(spectrum(**base, mz260=2.0), 291)["meox"] == "yes"
    assert em.analyze(spectrum(**{**base, "mz276": 30}, mz260=2.9), 291)["meox"] == "undecided"       # 2.9 < 0.1 x 30
    assert em.analyze(spectrum(**{**base, "mz276": 30}, mz260=3.0), 291)["meox"] == "yes"
    assert marker(em.analyze(spectrum(**base, mz260=1.9), 291), "M-31")["min"] == em.MIN_M31


def test_m_minus_57_makes_it_tbdms_only_when_it_is_the_strong_ion():
    strong = em.analyze(spectrum(mz73=15, mz289=100, mz331=6), 346)           # M-57 = 289 (100 %), M-15 = 331 (6 %), m/z 73 as the silicon ion
    assert marker(strong, "M-57")["mz"] == 289 and marker(strong, "M-15")["mz"] == 331
    assert strong["tbdms"] and strong["derivatization"] == "TBDMS" and strong["meox"] == "n/a" and "tert-butyl" in " ".join(strong["reasons"])
    weak = em.analyze(spectrum(mz73=100, mz75=20, mz289=19, mz331=10), 346)            # M-57 under 20 %
    assert not weak["tbdms"] and weak["derivatization"] == "TMS"
    tms_like = em.analyze(spectrum(mz73=100, mz75=30, mz289=25, mz331=60), 346)        # M-57 over 20 % but M-15 is stronger: still TMS
    assert not tms_like["tbdms"] and tms_like["derivatization"] == "TMS"
    assert not em.analyze(spectrum(mz73=15, mz289=100, mz331=6)).get("tbdms")  # without M, TBDMS cannot be told


def test_m_minus_57_needs_a_silicon_ion_to_make_it_tbdms():
    """A butyl or propionyl loss gives M-57 in 7.5 % of underivatized spectra; with m/z 73, 75 or 115 >= 2 % it is 2.3 % (module docstring)."""
    no_si = em.analyze(spectrum(mz57=30, mz289=100, mz331=6), 346)
    assert not no_si["tbdms"] and no_si["derivatization"] == "None" and marker(no_si, "M-57")["present"] is True
    for ion in ("mz73", "mz75", "mz115"):
        assert em.analyze(spectrum(**{ion: 2}, mz289=100, mz331=6), 346)["tbdms"] is True, ion
        assert em.analyze(spectrum(**{ion: 1.9}, mz289=100, mz331=6), 346)["tbdms"] is False, ion
    assert marker(no_si, "115")["min"] == em.MIN_TBDMS_SI


def test_no_marker_means_underivatized():
    a = em.analyze(spectrum(mz43=100, mz57=80, mz71=40, mz85=20))
    assert a["tms"] == "no" and a["tms_min"] == 0 and a["derivatization"] == "None" and a["meox"] == "n/a" and not a["tbdms"]
    assert "give the molecular ion" in a["reasons"][0] and "0.0 %" in a["reasons"][0]
    assert em.analyze(spectrum(mz43=100, mz57=80), 150)["derivatization"] == "None"


def test_m147_stands_in_for_m73_only_when_m73_cannot_be_seen():
    scan_from_85 = em.analyze(spectrum(mz86=40, mz100=22, mz147=30, mz174=100))
    assert scan_from_85["tms"] == "yes" and scan_from_85["tms_min"] == 2 and "outside the recorded range" in " ".join(scan_from_85["reasons"])
    in_range = em.analyze(spectrum(mz60=40, mz86=40, mz147=30, mz174=100))      # 73 could be seen and is absent
    assert in_range["tms"] == "no" and in_range["derivatization"] == "None"
    with_m15 = em.analyze(spectrum(mz60=40, mz147=30, mz174=100, mz290=2), 305)  # 147 and M-15 together are no stand-in (4.3 points more false positives)
    assert with_m15["tms"] == "no" and with_m15["derivatization"] == "None"


def test_m147_alone_needs_5_percent_where_m73_is_out_of_range():
    assert em.analyze(spectrum(mz86=40, mz147=5, mz174=100))["tms"] == "yes" and em.analyze(spectrum(mz86=40, mz147=4.9, mz174=100))["tms"] == "no"
    weak = em.analyze(spectrum(mz86=40, mz147=3, mz174=100))
    assert "m/z 147 is 3.0 %" in weak["reasons"][0] and "needs 5 %" in weak["reasons"][0]
    assert em.analyze(spectrum(mz86=40, mz174=100))["reasons"][0].endswith("m/z 147 is absent; give the molecular ion to look for TBDMS")


def test_m73_alone_is_not_tms_it_needs_m75_or_m147():
    """Other ions at m/z 73 (C4H9O+, C3H5O2+) are in 7 % of underivatized spectra; with m/z 75 >= 1 % or 147 >= 2 % it is 2.7 % (module docstring)."""
    alone = em.analyze(spectrum(mz45=60, mz73=100, mz87=30))
    assert alone["tms"] == "no" and alone["derivatization"] == "None" and alone["tms_min"] == 0 and marker(alone, "73")["present"] is True
    assert "m/z 73 is 100.0 % of the base peak but m/z 75 (0.0 %) and m/z 147 (0.0 %) do not back it" in alone["reasons"][0]
    backed = em.analyze(spectrum(mz45=60, mz73=100, mz75=1.5, mz87=30))
    assert backed["tms"] == "yes" and backed["tms_min"] == 1 and "backed by m/z 75 (2 %)" in backed["reasons"][0]
    assert em.analyze(spectrum(mz73=100, mz147=2.5))["tms_min"] == 2


def test_a_molecular_ion_that_does_not_fit_the_spectrum_is_reported():
    far = em.analyze(spectrum(mz73=100, mz147=20, mz174=50), 500)
    assert far["m_check"] == "above spectrum" and "above the highest peak" in far["warnings"][0]
    nothing = em.analyze(spectrum(mz73=100, mz147=20, mz174=50, mz300=40), 311)
    assert nothing["m_check"] == "not seen" and "check that value" in nothing["warnings"][0]
    assert em.analyze(spectrum(mz73=100, mz276=5), 291)["m_check"] == "ok"


def test_check_selection_keeps_the_choice_and_explains_a_clash():
    tms = em.analyze(spectrum(mz73=100, mz75=20, mz147=30, mz276=5), 291)
    assert em.check_selection(tms, "TMS") == [] and em.check_selection(tms, "MeOx + TMS") == []
    none_vs_tms = em.check_selection(tms, "None")
    assert len(none_vs_tms) == 1 and "chose None" in none_vs_tms[0] and "Your choice is kept" in none_vs_tms[0]
    none = em.analyze(spectrum(mz43=100, mz57=80))
    assert em.check_selection(none, "None") == []
    weak73 = em.check_selection(none, "TMS")
    assert len(weak73) == 1 and "m/z 73" in weak73[0] and "0.0 %" in weak73[0]
    unbacked = em.check_selection(em.analyze(spectrum(mz45=60, mz73=100, mz87=30)), "TMS")
    assert len(unbacked) == 1 and "m/z 73 (100.0 %) is not backed by m/z 75 or 147" in unbacked[0] and "Your choice is kept" in unbacked[0]
    tbdms = em.analyze(spectrum(mz75=10, mz289=100, mz331=6), 346)
    assert "M-57 points to TBDMS" in em.check_selection(tbdms, "TMS")[0]
    assert "M-57 is only" in em.check_selection(em.analyze(spectrum(mz73=100, mz75=20, mz331=6), 346), "TBDMS")[0]
    assert em.check_selection(tbdms, "TBDMS") == []
    blind = em.analyze(spectrum(mz86=40, mz174=100))                          # 73 not recorded: a TMS choice is not contradicted
    assert em.check_selection(blind, "TMS") == []


def test_state_and_family_names():
    assert em.state_of_spectrum(spectrum(mz73=100, mz147=20)) == "TMS" and em.state_of_spectrum(spectrum(mz43=100)) == "None"
    assert em.state_of_spectrum(spectrum(mz75=10, mz289=100, mz331=5), 346) == "TBDMS"
    assert em.state_of_spectrum(spectrum(mz73=100, mz75=20, mz260=3, mz276=5), 291) == "TMS"        # MeOx + TMS reads as TMS
    assert [em.family(m) for m in ("None", "TMS", "MeOx + TMS", "TBDMS", "MeOx + TBDMS")] == ["None", "TMS", "TMS", "TBDMS", "TBDMS"]
