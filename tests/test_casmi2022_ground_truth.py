"""CASMI 2022 official key and extracted-spectrum export."""
import csv
import importlib.util
from pathlib import Path

import pytest

from mzdesk import adducts
from mzdesk.annotation.pipeline import parse_inputs


BASE = Path(__file__).resolve().parents[1] / "benchmark" / "casmi2022"
DATA = BASE / "casmi2022_ground_truth.csv"
if not (BASE / "build_ground_truth.py").is_file():
    pytest.skip("the CASMI 2022 benchmark files are not in this checkout", allow_module_level=True)


def load_builder():
    spec = importlib.util.spec_from_file_location("casmi_truth_builder", BASE / "build_ground_truth.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BUILDER = load_builder()


@pytest.fixture(scope="module")
def records():
    """Records loaded from the flat one-row-per-challenge table."""
    return BUILDER.load_records(DATA)


@pytest.fixture(scope="module")
def extracted_records():
    if not (BASE / "sirius_run/ms/challenge_001.ms").exists():
        pytest.skip("blind extracted .ms files are local and not tracked")
    return BUILDER.build_extracted_records()


def test_flat_table_is_one_row_per_challenge_and_round_trips(tmp_path):
    with DATA.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header, rows = next(reader), list(reader)
    assert header == BUILDER.NAMES
    assert len(rows) == 500
    reloaded = tmp_path / "reloaded.csv"
    BUILDER.write_records(records=BUILDER.load_records(DATA), path=reloaded)
    same = reloaded.read_bytes() == DATA.read_bytes()      # a bool: pytest's diff of two 500-row files ran for hours on Windows (CRLF, 2026-11-13)
    assert same, f"{reloaded} differs from {DATA} (line endings? the key must be LF, see .gitattributes)"


def test_all_records_match_curated_key_and_have_expected_fields(records):
    with (BASE / "metadata/casmi2022_challenges.csv").open(newline="", encoding="utf-8") as handle:
        metadata = {int(row["challenge_id"]): row for row in csv.DictReader(handle)}
    assert [r["challenge_id"] for r in records] == list(range(1, 501))
    assert set(metadata) == set(range(1, 501))
    for record in records:
        m = metadata[record["challenge_id"]]
        assert set(record) == set(BUILDER.NAMES)
        for actual, metadata_key in (("formula", "formula"), ("smiles_official", "smiles"),
                                     ("inchikey", "inchikey"), ("adduct_official", "adduct"),
                                     ("priority", "priority"), ("file", "file"),
                                     ("collection", "collection"), ("rack_id", "rack_id"),
                                     ("well_id", "well_id"), ("mix_or_set", "mix_or_set")):
            assert record[actual] == m[metadata_key].strip()
        for actual, metadata_key in (("monoisotopic_mass", "monoisotopic_mass"),
                                     ("rt_min_given", "rt_min"), ("precursor_mz_given", "precursor_mz")):
            assert record[actual] == pytest.approx(float(m[metadata_key]), abs=1e-7)


def test_extracted_records_still_have_blind_valid_inputs(extracted_records):
    for record in extracted_records:
        truth, inputs = record["ground_truth"], record["inputs"]
        assert truth["adduct_mzdesk"] == adducts.canonical(truth["adduct_official"])
        assert inputs["precursor"]["adduct"] == "Auto"
        assert set(inputs) == {"precursor", "ms1", "ms2", "libraries", "tol_ppm", "tol_da", "omics",
                               "halogen", "local_isomers", "additive", "alkali", "ph"}
        query = parse_inputs(inputs)
        assert query.polarity in ("Positive", "Negative")
        assert query.ms1 and query.ms2
        assert record["ms2_by_collision_energy"]


def test_csv_spectra_match_extracted_records(records, extracted_records):
    extracted = {r["challenge_id"]: r for r in extracted_records}
    for record in records:
        inputs = extracted[record["challenge_id"]]["inputs"]
        assert record["ms1_isotope"] == inputs["ms1"]
        assert record["ms2"] == inputs["ms2"]["value"]
        assert record["ms1_isotope"] and record["ms2"]
    assert len(records[250]["ms1_isotope"].splitlines()) == 1


def test_extraction_caveats_are_retained(extracted_records):
    by_id = {r["challenge_id"]: r for r in extracted_records}
    assert {r["challenge_id"] for r in extracted_records if r["provenance"]["extraction_status"] != "OK"} == {18, 393}
    assert "neighbouring scan" in by_id[18]["provenance"]["extraction_notes"]
    assert "widened" in by_id[393]["provenance"]["extraction_notes"]
    assert len(by_id[251]["inputs"]["ms1"].splitlines()) == 1
    assert sum("stereo_missing_in_smiles" in r["ground_truth"]["structure_notes"] for r in extracted_records) == 277


def test_official_final_key_and_corrections(records):
    openpyxl = pytest.importorskip("openpyxl")
    source = BASE / "metadata/source"
    by_id = {r["challenge_id"]: r for r in records}

    def rows(name):
        book = openpyxl.load_workbook(source / name, read_only=True, data_only=True)
        try:
            lines = book.active.iter_rows(values_only=True)
            next(lines)  # title
            headers = next(lines)
            return [dict(zip(headers, row)) for row in lines if row[0] is not None]
        finally:
            book.close()

    key = rows("casmi2022_key_final.xlsx")
    assert len(key) == 500
    for row in key:
        record = by_id[int(row["Compound Number"])]
        for name, official in BUILDER.OFFICIAL_HEADERS.items():
            if BUILDER.KIND[name] in ("int", "float"):
                assert record[name] == pytest.approx(float(row[official]), abs=1e-7)
            else:
                value = row[official]
                if name in ("rack_id", "mix_or_set") and isinstance(value, (int, float)) and value == int(value):
                    value = int(value)
                assert record[name] == str(value).strip()

    corrections = rows("priority_corrections.xlsx") + rows("bonus_corrections.xlsx")
    assert len(corrections) == 4
    for row in corrections:
        record = by_id[int(row["Compound Number"])]
        assert str(row["File"]).strip() == record["file"]
        assert float(row["RT [min]"]) == pytest.approx(record["rt_min_given"], abs=1e-7)
        assert float(row["Precursor m/z (Da)"]) == pytest.approx(record["precursor_mz_given"], abs=1e-7)


def test_rebuild_matches_packaged_data(records, extracted_records):
    pytest.importorskip("openpyxl")
    assert BUILDER.build_records() == records


def test_merge_uses_first_peak_for_window_and_strongest_peak_for_representative():
    blocks = [{"peaks": [[100.0, 2.0], [100.0009, 10.0], [100.0018, 5.0]]}]
    assert BUILDER.merge_peaks(blocks) == [[100.0009, 10.0], [100.0018, 5.0]]
