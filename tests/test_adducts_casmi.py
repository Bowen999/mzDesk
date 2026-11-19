"""CASMI 2022 adducts end to end: notation, polarity, engine names, mass arithmetic and the library adduct filter.

The six adducts ([M+H]+ [M-H]- [M+NH4]+ [M+HCOO]- [M+Na]+ [M+Cl]-) must work in every component the ESI MS/MS workflow uses.
Offline only: no MS-FINDER / MetFrag / FRIGID run, fixture libraries in a temporary MZDESK_LIBRARY_DIR. The CASMI metadata is read
for its adduct column and its precursor arithmetic only; nothing (prior, threshold) is derived from it.
"""
import csv
from collections import Counter
from pathlib import Path

import pytest
from msbuddy.base import Adduct

from mzdesk import adducts, spec
from mzdesk.ms1 import msbuddy_engine as bd
from mzdesk.ms2 import frigid_engine as fe
from mzdesk.ms2 import library_index as lx
from mzdesk.ms2 import library_parse as lp
from mzdesk.ms2 import library_search as ls
from mzdesk.ms2 import metfrag_engine as mf
from mzdesk.ms2 import metric_names as mn
from mzdesk.ms2 import msfinder_engine as mfe
from mzdesk.platform_fit import molecular_properties as mp

CASMI = Path(__file__).resolve().parents[1] / "benchmark" / "casmi2022" / "metadata" / "casmi2022_challenges.csv"
needs_casmi = pytest.mark.skipif(not CASMI.is_file(), reason="CASMI 2022 metadata not in this checkout")
SIX = adducts.CASMI_ADDUCTS
# m/z(ion) - M, electron included (formate = HCOO + e = 44.998203, not 44.997654)
SHIFT = {"[M+H]+": 1.007276452, "[M-H]-": -1.007276452, "[M+NH4]+": 18.033825553, "[M+HCOO]-": 44.998202851,
         "[M+Na]+": 22.989220701, "[M+Cl]-": 34.969401260}
# adduct -> name for (buddy, msfinder, metfrag, frigid); None = the engine does not take it (FRIGID is positive only)
ENGINE_NAMES = {"[M+H]+": ("[M+H]+", "[M+H]+", "[M+H]+", "[M+H]+"),
                "[M-H]-": ("[M-H]-", "[M-H]-", "[M-H]-", None),
                "[M+NH4]+": ("[M+NH4]+", "[M+NH4]+", "[M+NH4]+", "[M+NH4]+"),
                "[M+HCOO]-": ("[M+HCOO]-", "[M+FA-H]-", "[M+HCOO]-", None),
                "[M+Na]+": ("[M+Na]+", "[M+Na]+", "[M+Na]+", "[M+Na]+"),
                "[M+Cl]-": ("[M+Cl]-", "[M+Cl]-", "[M+Cl]-", None)}
FORMATE_SPELLINGS = ("[M+HCOO]-", "[M+FA-H]-", "[M+HCOOH]-", "[M+FA]-")


# ---------------------------------------------------------------- notation
@needs_casmi
def test_casmi_2022_uses_exactly_the_six_adducts():
    counts = Counter(adducts.canonical(r["adduct"]) for r in csv.DictReader(CASMI.open(encoding="utf-8")))
    assert counts == {"[M+H]+": 177, "[M-H]-": 108, "[M+NH4]+": 67, "[M+HCOO]-": 67, "[M+Na]+": 60, "[M+Cl]-": 21}
    assert set(counts) == set(SIX)


@pytest.mark.parametrize("raw,expected", [
    # CASMI metadata (5 rows carry a trailing space)
    ("[M+FA-H]-", "[M+HCOO]-"), ("[M+FA-H]- ", "[M+HCOO]-"),
    # MassBank Precursor_type
    ("[M+HCOOH-H]-", "[M+HCOO]-"), ("[M+HCO2]-", "[M+HCOO]-"), ("[M+Na]+*", "[M+Na]+"),
    # GNPS name suffixes; [M+FA]- (34,611 records) is the formate adduct: precursor = M + 44.998, not M + 46.006
    ("[M+FA]-", "[M+HCOO]-"), ("[M+HCOOH]-", "[M+HCOO]-"), ("M+formate", "[M+HCOO]-"), ("M+FA-", "[M+HCOO]-"),
    ("M+FA-H", "[M+HCOO]-"), ("M+HCOO", "[M+HCOO]-"), ("[M-H+HCOOH]", "[M+HCOO]-"),
    ("M+Cl-", "[M+Cl]-"), ("M+Cl", "[M+Cl]-"), ("[M+Cl]", "[M+Cl]-"), ("M-H-", "[M-H]-"), ("M-H", "[M-H]-"), ("[M-H]", "[M-H]-"),
    ("M+NH4", "[M+NH4]+"), ("[M+NH4]", "[M+NH4]+"), ("M+Na", "[M+Na]+"), ("[M+Na]", "[M+Na]+"), ("M+H", "[M+H]+"), ("[M+H]", "[M+H]+"),
    # not formate adducts / not guessed: a neutral FA cluster on a charged ion, a real loss, sign that contradicts the composition
    ("[M+H+HCOOH]+", "[M+H+HCOOH]+"), ("[M+Cl+HCOOH]-", "[M+Cl+HCOOH]-"), ("[M-HCOO]-", "[M-HCOO]-"),
    ("M-2H-", None), ("M-Ac-H-", None), ("Auto", None), ("", None),
])
def test_library_and_casmi_spellings_fold_to_the_six(raw, expected):
    assert adducts.canonical(raw) == lp.normalize_adduct(raw) == expected


def test_polarity():
    assert [adducts.polarity(a) for a in SIX] == ["Positive", "Negative", "Positive", "Negative", "Positive", "Negative"]
    assert adducts.polarity("[M+FA-H]-") == "Negative" and adducts.polarity("Auto") is None and adducts.polarity("banana") is None


def test_the_six_are_on_the_workflow_and_engine_pages():
    for a in SIX:
        mode = adducts.polarity(a)
        assert a in (spec.ADDUCTS_POS if mode == "Positive" else spec.ADDUCTS_NEG)          # master list of the workflow
        assert a in spec.BUDDY_ADDUCTS[mode] and a in spec.METFRAG_ADDUCTS[mode]
        assert adducts.engine_name(a, "msfinder") in spec.MSFINDER_ADDUCTS[mode]
        assert (adducts.engine_name(a, "frigid") in spec.FRIGID_ADDUCTS[mode]) == (mode == "Positive")


# ---------------------------------------------------------------- engine names
@pytest.mark.parametrize("adduct", SIX)
def test_engine_names_are_accepted_by_each_engine(adduct):
    positive = adducts.polarity(adduct) == "Positive"
    buddy, msfinder, metfrag, frigid = ENGINE_NAMES[adduct]
    assert tuple(adducts.engine_name(adduct, e) for e in adducts.ENGINES) == ENGINE_NAMES[adduct]
    assert adducts.supported_engines(adduct) == tuple(e for e, n in zip(adducts.ENGINES, ENGINE_NAMES[adduct]) if n)
    assert buddy in (bd.ADDUCTS_POSITIVE if positive else bd.ADDUCTS_NEGATIVE)
    ion = bd._parse_adduct(buddy)
    assert ion is not None and (ion.charge > 0) == positive
    assert mfe.resolve_adduct(positive, msfinder) == msfinder
    assert f"PRECURSORTYPE: {msfinder}" in mfe.build_query(300.0, positive, msfinder, [(100.0, 10.0)])
    assert mf.ADDUCTS[metfrag]["positive"] == positive
    if positive:
        assert fe.resolve_adduct(frigid)[0] == frigid
    else:
        assert frigid is None
        with pytest.raises(ValueError):
            fe.resolve_adduct(adduct)


def test_every_formate_spelling_maps_to_the_same_engine_names():
    for engine in adducts.ENGINES:
        assert {adducts.engine_name(s, engine) for s in FORMATE_SPELLINGS} == {ENGINE_NAMES["[M+HCOO]-"][adducts.ENGINES.index(engine)]}
    assert adducts.engine_name("Auto", "buddy") is None and adducts.supported_engines("banana") == ()
    with pytest.raises(ValueError, match="Unknown engine"):
        adducts.engine_name("[M+H]+", "sirius")


# ---------------------------------------------------------------- mass arithmetic
@pytest.mark.parametrize("adduct", SIX)
def test_neutral_mz_neutral_round_trip(adduct):
    shift = SHIFT[adduct]
    ion = Adduct(adduct, adducts.polarity(adduct) == "Positive", report_invalid=True)
    sign = 1 if ion.charge > 0 else -1
    for mass in (60.021129, 194.080376, 342.116212, 1234.56789):
        mz = bd.ion_mz_from_mass(mass, adduct)                    # the repo's calculator (msbuddy parser, float64 masses)
        assert mz - mass == pytest.approx(shift, abs=1e-6)
        assert abs(bd.mass_from_ion_mz(mz, adduct) - mass) < 1e-6
        # msbuddy's own arithmetic (main.py forward, query.py inverse; float32 element masses)
        mz_mb = mass + float(ion.net_formula.mass) - sign * 0.00054858
        assert abs(mz_mb - mz) < 1e-6
        assert abs((mz_mb + sign * 0.00054858 - float(ion.net_formula.mass)) - mass) < 1e-6
    # every spelling of the adduct gives the same m/z once canonicalised
    for spelling in (FORMATE_SPELLINGS if adduct == "[M+HCOO]-" else (adduct,)):
        assert bd.ion_mz_from_mass(0.0, adducts.canonical(spelling)) == pytest.approx(shift, abs=1e-9)
    # the other adduct tables of the repo agree (MetFrag keeps 7-decimal constants: 1.9e-6 Da for formate)
    assert mp.ADDUCTS[adduct][0] == pytest.approx(shift, abs=1e-6)
    assert mf.ADDUCTS[adduct]["shift"] == pytest.approx(shift, abs=5e-6)
    if adduct in fe.ADDUCT_TABLE:
        assert fe.ion_mz(0.0, adduct) == pytest.approx(shift, abs=1e-6)


@needs_casmi
def test_casmi_precursors_match_the_calculator():
    """Arithmetic check only: metadata precursor m/z vs monoisotopic mass + adduct (1 ppm would expose a missing electron)."""
    for r in csv.DictReader(CASMI.open(encoding="utf-8")):
        theo = bd.ion_mz_from_mass(float(r["monoisotopic_mass"]), adducts.canonical(r["adduct"]))
        assert abs(float(r["precursor_mz"]) - theo) / theo * 1e6 < 1.0, r["challenge_id"]


# ---------------------------------------------------------------- library adduct filter
SUCROSE_FA = 387.11441                                  # C12H22O11 as [M+HCOO]-
PEAKS = ["341.1089 100", "179.0561 45", "119.0350 20", "89.0244 30"]


def _massbank(acc, precursor_type):
    rec, why = lp.normalise_massbank({"DB#": acc, "Name": "Sucrose", "Precursor_type": precursor_type, "Spectrum_type": "MS2",
                                      "PrecursorMZ": str(SUCROSE_FA), "Ion_mode": "NEGATIVE", "Instrument_type": "LC-ESI-QFT",
                                      "_peaks": PEAKS})
    assert why is None
    return rec


def _gnps(acc, name):
    rec, why = lp.normalise_gnps({"MSLEVEL": "2", "PEPMASS": str(SUCROSE_FA), "IONMODE": "Negative", "NAME": name,
                                  "SPECTRUMID": acc, "SOURCE_INSTRUMENT": "LC-ESI-qTof", "_peaks": PEAKS})
    assert why is None
    return rec


@pytest.fixture
def formate_libs(monkeypatch, tmp_path):
    """MassBank- and GNPS-style records of one spectrum under every formate spelling, plus [M+Cl]- / [M-H]- decoys at the same m/z.

    GNPS-LEGACY is stored as [M+HCOOH]-, the label the old parser gave GNPS's [M+FA]- (the installed index still has it).
    """
    monkeypatch.setenv("MZDESK_LIBRARY_DIR", str(tmp_path))
    lx.clear_cache()
    legacy = _gnps("GNPS-LEGACY", "Sucrose")
    legacy["adduct"] = "[M+HCOOH]-"
    lx.build("MassBank", [_massbank("MB-FA-H", "[M+FA-H]-"), _massbank("MB-HCOO", "[M+HCOO]-"), _massbank("MB-CL", "[M+Cl]-")],
             {"version": "t1"})
    lx.build("GNPS", [_gnps("GNPS-FA", "Sucrose [M+FA]-"), _gnps("GNPS-FORMATE", "Sucrose M+formate"), _gnps("GNPS-FA-SIGN", "Sucrose M+FA-"),
                      _gnps("GNPS-HCOOH", "Sucrose [M+HCOOH]-"), legacy, _gnps("GNPS-MH", "Sucrose [M-H]-")], {"version": "t1"})
    yield
    lx.clear_cache()


FORMATE_ROWS = {"MB-FA-H", "MB-HCOO", "GNPS-FA", "GNPS-FORMATE", "GNPS-FA-SIGN", "GNPS-HCOOH", "GNPS-LEGACY"}


def _hits(adduct):
    res = ls.search([tuple(map(float, p.split())) for p in PEAKS], "MS/MS", ["MassBank", "GNPS"], [mn.ENTROPY_WEIGHTED], -1,
                    precursor_mz=SUCROSE_FA, positive=False, adduct=adduct)
    return {h["accession"]: h["adduct"] for h in res["hits"]}


def test_formate_filter_finds_every_spelling_including_the_legacy_index_label(formate_libs):
    assert "[M+HCOOH]-" in lx.load("GNPS").meta["adducts"]                  # the legacy vocabulary path is exercised
    for query in FORMATE_SPELLINGS:
        hits = _hits(query)
        assert set(hits) == FORMATE_ROWS, query
        assert set(hits.values()) == {"[M+HCOO]-"}                         # hits are reported in canonical form
    assert set(_hits("[M+Cl]-")) == {"MB-CL"} and set(_hits("[M-H]-")) == {"GNPS-MH"}
    assert set(_hits(None)) == FORMATE_ROWS | {"MB-CL", "GNPS-MH"}
