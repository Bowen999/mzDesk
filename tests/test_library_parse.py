"""Library record parsing: structures that GNPS only gives as a full InChI, and "N/A" placeholders."""
import pytest

from mzdesk.annotation import library_step as step
from mzdesk.ms2 import library_parse as lp

# CCMSLIB00004699083: SMILES=N/A, INCHIAUX=N/A, no INCHIKEY, structure only as INCHI
INCHI = ("InChI=1S/C19H28O10/c20-8-12-14(23)15(24)17(29-18-16(25)13(22)11(21)9-27-18)19(28-12)26-7-6-10-4-2-1-3-5-10/"
         "h1-5,11-25H,6-9H2/t11-,12-,13+,14-,15+,16-,17-,18+,19-/m1/s1")
PEAKS = ["100.0 10", "200.0 100"]


def _gnps(**fields):
    rec = {"MSLEVEL": "2", "PEPMASS": "434.202", "IONMODE": "positive", "NAME": "x [M+NH4]+", "SPECTRUMID": "CCMSLIB00004699083",
           "SMILES": "N/A", "INCHIAUX": "N/A", "_peaks": PEAKS}
    rec.update(fields)
    out, why = lp.normalise_gnps(rec)
    assert why is None
    return out


def test_gnps_inchi_only_record_gets_smiles_and_key():
    r = _gnps(INCHI=INCHI)
    assert r["inchikey"].startswith("CGZNTHDUAFIKOH-") and len(r["inchikey"]) == 27
    assert r["smiles"] and r["smiles"] != "N/A"


@pytest.mark.parametrize("wrapped", [f'"{INCHI}"', f"'{INCHI}'", f"  {INCHI}\t"])
def test_gnps_quoted_inchi(wrapped):
    r = _gnps(INCHI=wrapped)
    assert r["inchikey"].startswith("CGZNTHDUAFIKOH-") and r["smiles"]


def test_gnps_inchi_in_inchiaux():
    assert _gnps(INCHI="N/A", INCHIAUX=INCHI)["inchikey"].startswith("CGZNTHDUAFIKOH-")


@pytest.mark.parametrize("junk", ["InChI=1S/garbage/zzz", "InChI=1S/", "not an inchi", "N/A", ""])
def test_gnps_garbage_inchi_stays_empty(junk):
    r = _gnps(INCHI=junk)
    assert r["smiles"] == "" and r["inchikey"] == ""


def test_gnps_existing_smiles_and_key_untouched():
    r = _gnps(SMILES="CCO", INCHIKEY="LFQSCWFLJHTTHZ-UHFFFAOYSA-N", INCHI=INCHI)
    assert (r["smiles"], r["inchikey"]) == ("CCO", "LFQSCWFLJHTTHZ-UHFFFAOYSA-N")


def test_gnps_smiles_without_key_is_not_keyed_at_parse_time():
    r = _gnps(SMILES="CCO")
    assert r["smiles"] == "CCO" and r["inchikey"] == ""


def test_massbank_na_placeholders_are_empty():
    rec, why = lp.normalise_massbank({"DB#": "MB1", "Name": "x", "Precursor_type": "[M+H]+", "Spectrum_type": "MS2", "PrecursorMZ": "100.1",
                                      "Ion_mode": "POSITIVE", "SMILES": "N/A", "InChIKey": "N/A", "_peaks": PEAKS})
    assert why is None and rec["smiles"] == "" and rec["inchikey"] == ""


def test_hit_without_key_gets_it_from_smiles():
    assert step._inchikey({"inchikey": "", "smiles": "CCO"}) == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    assert step._inchikey({"inchikey": "N/A", "smiles": "CCO"}) == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    assert step._inchikey({"inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C", "smiles": "CCO"}) == "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
    assert step._inchikey({"inchikey": "", "smiles": "N/A"}) is None
    assert step._inchikey({"inchikey": "", "smiles": ""}) is None
