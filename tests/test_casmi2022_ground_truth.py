"""CASMI 2022 ground-truth package: blind inputs, labels and provenance."""
import csv
import importlib.util
import json
from pathlib import Path

import pytest

from mzdesk import adducts
from mzdesk.annotation.pipeline import parse_inputs


BASE = Path(__file__).resolve().parents[1] / "benchmark" / "casmi2022"
DATA = BASE / "casmi2022_ground_truth.jsonl"


@pytest.fixture(scope="module")
def records():
    with DATA.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle]


def test_all_records_match_curated_key_and_have_blind_valid_inputs(records):
    with (BASE / "metadata/casmi2022_challenges.csv").open(newline="", encoding="utf-8") as handle:
        metadata = {int(row["challenge_id"]): row for row in csv.DictReader(handle)}
    assert [r["challenge_id"] for r in records] == list(range(1, 501))
    assert set(metadata) == set(range(1, 501))
    for record in records:
        m = metadata[record["challenge_id"]]
        truth, inputs = record["ground_truth"], record["inputs"]
        assert truth["formula"] == m["formula"]
        assert truth["smiles_official"] == m["smiles"].strip()
        assert truth["inchikey"] == m["inchikey"]
        assert truth["inchikey14"] == m["inchikey"].split("-")[0]
        assert truth["adduct_official"] == m["adduct"].strip()
        assert truth["adduct_mzdesk"] == adducts.canonical(m["adduct"])
        assert inputs["precursor"]["adduct"] == "Auto"
        assert set(inputs) == {"precursor", "ms1", "ms2", "libraries", "tol_ppm", "tol_da", "omics",
                               "halogen", "local_isomers", "additive", "alkali", "ph"}
        query = parse_inputs(inputs)
        assert query.polarity.lower() == m["polarity"]
        assert query.ms1 and query.ms2
        assert record["ms2_by_collision_energy"]


def test_extraction_caveats_are_retained(records):
    by_id = {r["challenge_id"]: r for r in records}
    assert {r["challenge_id"] for r in records if r["provenance"]["extraction_status"] != "OK"} == {18, 393}
    assert "neighbouring scan" in by_id[18]["provenance"]["extraction_notes"]
    assert "widened" in by_id[393]["provenance"]["extraction_notes"]
    assert len(by_id[251]["inputs"]["ms1"].splitlines()) == 1
    assert sum("stereo_missing_in_smiles" in r["ground_truth"]["structure_notes"] for r in records) == 277


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
        truth, prov = record["ground_truth"], record["provenance"]
        for official, actual in (("Formula", truth["formula"]), ("Adduct", truth["adduct_official"]),
                                 ("InChIKey", truth["inchikey"]), ("SMILES", truth["smiles_official"]),
                                 ("File", prov["file"]), ("Priority/Bonus", prov["priority"])):
            assert str(row[official]).strip() == actual
        for official, actual in (("Monoisotopic Mass", truth["monoisotopic_mass"]),
                                 ("Precursor m/z (Da)", prov["precursor_mz_given"]), ("RT [min]", prov["rt_min_given"])):
            assert float(row[official]) == pytest.approx(actual, abs=1e-7)

    corrections = rows("priority_corrections.xlsx") + rows("bonus_corrections.xlsx")
    assert len(corrections) == 4
    for row in corrections:
        prov = by_id[int(row["Compound Number"])]["provenance"]
        assert str(row["File"]).strip() == prov["file"]
        assert float(row["RT [min]"]) == pytest.approx(prov["rt_min_given"], abs=1e-7)
        assert float(row["Precursor m/z (Da)"]) == pytest.approx(prov["precursor_mz_given"], abs=1e-7)


def test_rebuild_matches_packaged_data_when_extracted_files_are_available(records):
    if not (BASE / "sirius_run/ms/challenge_001.ms").exists():
        pytest.skip("blind extracted .ms files are local and not tracked")
    spec = importlib.util.spec_from_file_location("casmi_truth_builder", BASE / "build_ground_truth.py")
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    assert builder.build_records() == records


def test_merge_uses_first_peak_for_window_and_strongest_peak_for_representative():
    spec = importlib.util.spec_from_file_location("casmi_truth_builder", BASE / "build_ground_truth.py")
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    blocks = [{"peaks": [[100.0, 2.0], [100.0009, 10.0], [100.0018, 5.0]]}]
    assert builder.merge_peaks(blocks) == [[100.0009, 10.0], [100.0018, 5.0]]
