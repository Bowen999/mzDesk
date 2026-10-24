"""共享峰表解析：各种输入格式的识别与解析；MGF / MSP 文件只取第一条谱。"""
from mzdesk.shared.spectrum import detect_format, first_record, parse_spectrum

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
