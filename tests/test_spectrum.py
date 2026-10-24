"""共享峰表解析：MGF / MSP 文件只取第一条谱。"""
from mzdesk.shared.spectrum import first_record, parse_spectrum

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
