"""共享峰表解析：各种输入格式的识别与解析；MGF / MSP 文件只取第一条谱。"""
import json
from pathlib import Path

from mzdesk.shared.spectrum import (NOTE_DECIMAL_COMMA, NOTE_EXTRA_NUMBERS, NOTE_THOUSANDS, detect_format, fetch_usi, first_record, parse_spectrum,
                                    parse_spectrum_notes)

MGF = """BEGIN IONS
TITLE=caffeine
PEPMASS=195.0877 1000
CHARGE=1+
138.0662 100
110.0713 18
END IONS
BEGIN IONS
TITLE=second
PEPMASS=300
100 1
END IONS"""

MSP = """Name: PG 36:2
PrecursorMZ: 773.5338
Num Peaks: 2
152.99 24
773.53 416

Name: other
PrecursorMZ: 100
Num Peaks: 1
50 1
"""


def test_mgf_first_record():
    rec = first_record(MGF)
    assert rec["n_records"] == 2 and rec["name"] == "caffeine" and rec["precursor_mz"] == 195.0877
    assert parse_spectrum(rec["peaks_text"]) == [(110.0713, 18.0), (138.0662, 100.0)]


def test_msp_first_record():
    rec = first_record(MSP)
    assert rec["n_records"] == 2 and rec["name"] == "PG 36:2" and rec["precursor_mz"] == 773.5338
    assert parse_spectrum(rec["peaks_text"]) == [(152.99, 24.0), (773.53, 416.0)]   # "36:2" in the name is not a peak


def test_plain_peak_list_is_one_record():
    rec = first_record("100 1\n200 2")
    assert rec["n_records"] == 1 and rec["precursor_mz"] is None
    assert parse_spectrum(rec["peaks_text"]) == [(100.0, 1.0), (200.0, 2.0)]


# ---------------------------------------------------------------- formats
PEAKS = [(100.1, 5.0), (200.2, 100.0), (300.3, 40.0)]
INPUTS = {
    "GNPS peaks_json": "[[100.1, 5],\n [200.2, 100],\n [300.3, 40]]",
    "JSON objects": '[{"mz": 100.1, "intensity": 5}, {"m/z": 200.2, "int": 100}, {"mz": 300.3, "abundance": 40}]',
    "Two arrays (matchms)": "mz = [100.1, 200.2, 300.3]\nintensities = [5, 100, 40]",
    "HMDB XML": ("<ms-ms-peaks><ms-ms-peak><id>1</id><mass-charge>100.1</mass-charge><intensity>5.0</intensity></ms-ms-peak>"
                 "<ms-ms-peak><mass-charge>200.2</mass-charge><intensity>100</intensity></ms-ms-peak>"
                 "<ms-ms-peak><mass-charge>300.3</mass-charge><intensity>40</intensity></ms-ms-peak></ms-ms-peaks>"),
    "MassBank record": ("ACCESSION: MSBNK-1\nCH$EXACT_MASS: 194.08\nMS$FOCUSED_ION: PRECURSOR_M/Z 195.0877\nPK$NUM_PEAK: 3\n"
                        "PK$PEAK: m/z int. rel.int.\n  100.1 5 50\n  200.2 100 999\n  300.3 40 400\n//"),
    "MS-DIAL / MoNA (mz:int)": "100.1:5 200.2:100 300.3:40",
    "Peak list (one per line)": "100.1\t5\n200.2\t100\n300.3\t40",
    "MSP": "NAME: PG 36:2\nPRECURSORMZ: 195.0877\nNum Peaks: 3\n100.1\t5\n200.2\t100\n300.3\t40",
    "MGF": "BEGIN IONS\nTITLE=x\nPEPMASS=195.0877 1000\nCHARGE=1+\n100.1 5\n200.2 100\n300.3 40\nEND IONS",
    "CSV": "mz,intensity\n100.1,5\n200.2,100\n300.3,40",
    "SIRIUS .ms": ">compound x\n>parentmass 195.0877\n>ms2peaks\n100.1 5\n200.2 100\n300.3 40",
}


def test_every_input_format_parses_and_is_detected():
    for label, text in INPUTS.items():
        assert parse_spectrum(text) == PEAKS, label
        assert detect_format(text) == label, label


def test_json_variants_without_outer_brackets_or_with_newlines():
    assert parse_spectrum("[100.1, 5], [200.2, 100],\n[300.3, 40]") == PEAKS
    assert parse_spectrum('{"mz": [100.1, 200.2, 300.3], "intensities": [5, 100, 40]}') == PEAKS
    assert parse_spectrum("mz = np.array([100.1, 200.2, 300.3])\nintensity = np.array([5., 100., 40.])") == PEAKS


def test_record_header_numbers_never_become_peaks():
    mb = INPUTS["MassBank record"] + "\nCH$MONOISOTOPIC_MASS: 194.0804"
    assert [mz for mz, _ in parse_spectrum(mb)] == [100.1, 200.2, 300.3]
    assert [mz for mz, _ in parse_spectrum(INPUTS["MSP"])] == [100.1, 200.2, 300.3]      # 36:2 in the name, Num Peaks: 3
    assert [mz for mz, _ in parse_spectrum(INPUTS["MGF"])] == [100.1, 200.2, 300.3]      # PEPMASS 195.0877 1000


def test_detect_unknown_and_empty():
    assert detect_format("") == "Unknown" and detect_format("hello world") == "Unknown"


def test_precursor_from_massbank_and_sirius():
    assert first_record(INPUTS["MassBank record"])["precursor_mz"] == 195.0877
    assert first_record(INPUTS["SIRIUS .ms"])["precursor_mz"] == 195.0877


def test_semicolon_separated_peaks():
    two = [(138.0662, 100.0), (195.0877, 32.0)]
    for text in ("138.0662:100;195.0877:32", "138.0662:100; 195.0877:32;", "138.0662 100;195.0877 32", "138.0662 100; 195.0877 32;\n"):
        assert parse_spectrum(text) == two, text
    assert detect_format("138.0662:100;195.0877:32") == "MS-DIAL / MoNA (mz:int)"
    assert detect_format("138.0662 100;195.0877 32") == "Peak list (one per line)"
    assert parse_spectrum("Name: x\nComment: CE 20; 40 eV\n100 1\n") == [(100.0, 1.0)]


# ---------------------------------------------------------------- commas, ';' columns
CAFFEINE_PEAKS = [(110.0713, 18.0), (138.0662, 100.0), (195.0877, 32.0)]


def test_semicolon_columns_are_mz_and_intensity():
    # Excel's "CSV (semicolon)" export: one peak per line, m/z; intensity. Each number used to become a peak of intensity 1.
    for text in ("195.0877;100\n138.0662;60", "195.0877; 100\r\n138.0662; 60\r\n", "m/z;intensity\n195.0877;100\n138.0662;60\n110.0713;20"):
        assert parse_spectrum_notes(text)[1] == [], text
    assert parse_spectrum("195.0877;100\n138.0662;60") == [(138.0662, 60.0), (195.0877, 100.0)]
    assert parse_spectrum("m/z;intensity\n195.0877;100\n138.0662;60\n110.0713;20") == [(110.0713, 20.0), (138.0662, 60.0), (195.0877, 100.0)]
    assert detect_format("m/z;intensity\n195.0877;100\n138.0662;60") == "CSV"
    assert parse_spectrum("195.0877;100;[M+H]+\n138.0662;60;frag") == [(138.0662, 60.0), (195.0877, 100.0)]      # a text column after the two numbers
    assert parse_spectrum("195.0877;100;") == [(195.0877, 100.0)]                                                 # trailing ';'
    assert parse_spectrum("41;100\n55;85\n69;60") == [(41.0, 100.0), (55.0, 85.0), (69.0, 60.0)]               # an EI table with integer m/z


def test_a_semicolon_list_of_peaks_or_bare_mz_is_still_read_as_before():
    assert parse_spectrum("100.1;200.2;300.3") == [(100.1, 1.0), (200.2, 1.0), (300.3, 1.0)]                      # one line of numbers: m/z only
    assert parse_spectrum("138.0662 100; 195.0877 32;\n") == [(138.0662, 100.0), (195.0877, 32.0)]
    assert parse_spectrum("41,100;55,85") == [(41.0, 100.0), (55.0, 85.0)]                                        # EI "mz,int" pairs separated by ';'
    assert parse_spectrum("195.0877;;100") == [(100.0, 1.0), (195.0877, 1.0)]                                     # an empty column in the middle: not a table row


def test_decimal_commas_next_to_a_tab_space_or_semicolon_are_decimal_points():
    two = [(138.0662, 100.0), (195.0877, 32.0)]
    for text in ("138,0662\t100\n195,0877\t32", "138,0662 100\n195,0877 32", "138,0662;100\n195,0877;32", "m/z;intensity\n138,0662;100\n195,0877;32",
                 "138,0662;100;frag\n195,0877;32;x", "138,0662\t100.0\n195,0877\t32.0"):
        assert parse_spectrum_notes(text) == (two, [NOTE_DECIMAL_COMMA]), text
    # the intensity may use the comma too; the m/z column with a comma decides
    assert parse_spectrum_notes("138,0662;1523,5\n195,0877;32,25") == ([(138.0662, 1523.5), (195.0877, 32.25)], [NOTE_DECIMAL_COMMA])
    assert parse_spectrum_notes("138,066\t100\n195,087\t32") == ([(138.066, 100.0), (195.087, 32.0)], [NOTE_DECIMAL_COMMA])
    assert parse_spectrum("138,0662\t1,200") == [(138.0662, 1.2)]
    # one line of the text shows decimal commas, so a line with nothing next to the comma is one too
    assert parse_spectrum_notes("138,6621\n195,0877") == ([(138.6621, 1.0), (195.0877, 1.0)], [NOTE_DECIMAL_COMMA])


def test_a_comma_alone_still_separates_the_columns():
    # "41,100" = m/z 41, intensity 100 (the usual EI csv); nothing here may be called a decimal comma or flagged
    for text, peaks in (("41,100\n55,85\n69,60", [(41.0, 100.0), (55.0, 85.0), (69.0, 60.0)]), ("195.0877,100", [(195.0877, 100.0)]),
                        ("73,12345\n147,9876", [(73.0, 12345.0), (147.0, 9876.0)]), ("41,100.0\n55,85.3\n69,999.0", [(41.0, 100.0), (55.0, 85.3), (69.0, 999.0)]),
                        ("73,1234.5\n147,999.9", [(73.0, 1234.5), (147.0, 999.9)]), ("mz,intensity\n100.1,5\n200.2,100", [(100.1, 5.0), (200.2, 100.0)]),
                        ("195.0877, 100,\n138.0662 , 60", [(138.0662, 60.0), (195.0877, 100.0)])):
        assert parse_spectrum_notes(text) == (peaks, []), text


def test_thousands_commas_in_the_intensity():
    assert parse_spectrum_notes("100.0 1,200") == ([(100.0, 1200.0)], [NOTE_THOUSANDS])
    assert parse_spectrum_notes("100.0\t1,234.5\n200.1\t12,345.6") == ([(100.0, 1234.5), (200.1, 12345.6)], [NOTE_THOUSANDS])
    assert parse_spectrum_notes("100.0 1,234,567") == ([(100.0, 1234567.0)], [NOTE_THOUSANDS])
    assert parse_spectrum_notes("100.0;1,200\n200.1;500") == ([(100.0, 1200.0), (200.1, 500.0)], [NOTE_THOUSANDS])
    assert parse_spectrum("100 1,200") == [(100.0, 1200.0)]                       # integer m/z: 1,200 is a thousands group (the note says so)


def test_a_row_with_more_than_two_numbers_is_flagged():
    assert parse_spectrum_notes("195.0877 100 5\n138.0662 60 3") == ([(138.0662, 60.0), (195.0877, 100.0)], [NOTE_EXTRA_NUMBERS])
    assert parse_spectrum_notes("100.1,5,50\n200.2,100,100")[1] == [NOTE_EXTRA_NUMBERS]
    assert parse_spectrum_notes("138,0662 100 5") == ([(138.0662, 100.0)], [NOTE_DECIMAL_COMMA, NOTE_EXTRA_NUMBERS])
    assert parse_spectrum_notes("100.1\t5\tC3H5\n200.2\t100\tC4H7")[1] == []         # text after the two numbers is an annotation, not a third number


def test_every_input_format_reads_without_notes():
    for label, text in INPUTS.items():
        assert parse_spectrum_notes(text) == (PEAKS, []), label
    assert parse_spectrum_notes("BEGIN IONS\nPEPMASS=195.0877 1000\n138.0662 100\n110.0713 18\nEND IONS")[1] == []
    assert parse_spectrum_notes("54.03:7 56.04:88 58.1:3")[1] == [] and parse_spectrum_notes("184.07 : 100\n200.1 : 5")[1] == []
    assert parse_spectrum_notes("41\n55\n69") == ([(41.0, 1.0), (55.0, 1.0), (69.0, 1.0)], [])        # a bare integer m/z list (EI) is not suspicious


def test_peak_table_cases():
    """tests/fixtures/peak_table_cases.json: how the pasted tables of different separators / locales are read. The page's peak count
    (web/src/peakCount.ts, tested on the same file by web/tests/peak-count.test.mjs) must agree with this parser."""
    cases = json.loads((Path(__file__).parent / "fixtures" / "peak_table_cases.json").read_text(encoding="utf-8"))
    assert len(cases) >= 30
    for case in cases:
        assert parse_spectrum(first_record(case["text"])["peaks_text"]) == [tuple(p) for p in case["peaks"]], case["name"]


def test_sirius_ms1peaks_are_not_mixed_into_the_ms2_spectrum():
    text = (">compound x\n>parentmass 195.0877\n>ms1peaks\n195.0877 1000\n195.0910 120\n"
            ">ms2peaks\n138.0662 100\n110.0713 20\n")
    assert parse_spectrum(text) == [(110.0713, 20.0), (138.0662, 100.0)]
    assert detect_format(text) == "SIRIUS .ms"


def test_json_record_object_with_peaks_and_precursor():
    text = '{"peaks": [[100, 5], [200, 10]], "precursor_mz": 200.1, "name": "x"}'
    assert parse_spectrum(text) == [(100.0, 5.0), (200.0, 10.0)]
    assert detect_format(text) == "JSON record"
    rec = first_record(text)
    assert rec["precursor_mz"] == 200.1 and rec["name"] == "x"


def test_comment_before_msp_name_is_not_an_empty_first_record():
    text = "# exported by a tool\nName: caffeine\nPrecursorMZ: 195.0877\n100 1\n200 2\n"
    rec = first_record(text)
    assert rec["n_records"] == 1 and rec["name"] == "caffeine" and rec["precursor_mz"] == 195.0877
    assert parse_spectrum(text) == [(100.0, 1.0), (200.0, 2.0)]


def test_usi_resolver_maps_the_gnps_json(monkeypatch):
    from mzdesk.shared import spectrum as sp

    class Response:
        def __enter__(self): return self
        def __exit__(self, *exc): return False
        def read(self, n): return (b'{"n_peaks": 2, "precursor_mz": 195.0877, '
                                   b'"peaks": [[138.0662, 100], [110.0713, 20], [138.0662, 80]]}')

    seen = {}

    def fake_urlopen(request, timeout):
        seen["url"] = request.full_url
        return Response()

    monkeypatch.setattr(sp.urllib.request, "urlopen", fake_urlopen)
    peaks, precursor = fetch_usi("mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00005435822")
    assert peaks == [(110.0713, 20.0), (138.0662, 100.0)]       # duplicate m/z keeps the stronger peak
    assert precursor == 195.0877
    assert "CCMSLIB00005435822" in seen["url"]
