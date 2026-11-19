"""Spectrum in Which Set: the order of the Pan-ReDU enrichment, how the page shows it, the order of the blocks, and 0 versus missing. All offline."""
import csv

from mzdesk import spec
from mzdesk.api.runners import gnps
from mzdesk.gnps import compound_db_check as g


def row(column, value, datasets, files, dsw, **extra):
    return {"column": column, "value": value, "n_hit_datasets": datasets, "n_hit_files": files, "dsw_enrich": dsw, **extra}


def test_enrichment_rows_put_values_from_two_datasets_first_then_single_dataset_values_and_missing_last():
    rows = [row("SampleType", "missing value", 5, 50, 9.0),                       # missing is last in its column however strong
            row("SampleType", "single-strong", 1, 90, 80.0),
            row("SampleType", "multi-weak", 2, 3, 0.5),
            row("SampleType", "multi-none", 3, 9, None),                          # no enrichment: after the ones that have it
            row("SampleType", "multi-strong", 2, 2, 7.0),
            row("SampleType", "single-weak", 1, 1, 1.0),
            row("SampleType", "multi-tie-more-files", 2, 8, 7.0),
            row("SampleType", "single-none", 1, 5, None),
            row("SampleType", "", 2, 4, 3.0),                                     # blank is missing too
            row("SampleTypeSub1", "x", 1, 1, 1.0), row("HealthStatus", "h", 2, 2, 2.0), row("NCBITaxonomy", "9606|Homo sapiens", 2, 2, 2.0)]
    ordered = [(r["column"], r["value"]) for r in sorted(reversed(rows), key=g.redu_row_order)]
    st = [v for c, v in ordered if c == "SampleType"]
    assert st[:7] == ["multi-tie-more-files", "multi-strong", "multi-weak", "multi-none", "single-strong", "single-weak", "single-none"]
    assert sorted(st[7:]) == ["", "missing value"]
    assert [c for c, _ in ordered if c != "SampleType"] == ["SampleTypeSub1", "NCBITaxonomy", "HealthStatus"]          # the fixed column order


def _tsv(path, rows):
    cols = ["search", "hit_set", "column", "value", "n_hits", "n_hit_files", "n_hit_datasets", "bg_files", "bg_datasets", "file_rate", "file_enrich",
            "ds_rate", "ds_enrich", "dsw_hit_share", "dsw_bg_share", "dsw_enrich"]
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({"search": "masst_repository", "hit_set": "all", "n_hits": 3, "bg_files": 10, **r})
    return str(path)


def test_the_page_shows_two_decimals_above_1000_and_dashes_and_hides_constant_columns(tmp_path):
    path = _tsv(tmp_path / "r.tsv", [
        {"column": "SampleType", "value": "plant", "n_hit_files": 4, "n_hit_datasets": 2, "bg_datasets": 0, "file_enrich": "1.23456", "ds_enrich": "1000", "dsw_enrich": "12345.6"},
        {"column": "NCBITaxonomy", "value": "9606|Homo sapiens", "n_hit_files": 0, "n_hit_datasets": 1, "bg_datasets": 7, "file_enrich": "999.999", "ds_enrich": "", "dsw_enrich": "2"},
        {"column": "DOIDCommonName", "value": "missing value", "n_hit_files": 1, "n_hit_datasets": 1, "bg_datasets": 3, "file_enrich": "0.004", "ds_enrich": "3", "dsw_enrich": ""}])
    block = gnps.redu_preview(path)
    assert block["columns"] == ["Property", "Value", "Hit files", "Hit datasets", "Background datasets", "File enrichment", "Dataset enrichment",
                                "Dataset-weighted enrichment"]                                    # Search and Hit set have one value each
    assert block["count"] == 3 and block["rows"] == [
        ["Sample type", "plant", "4", "2", "0", "1.23", ">1000", ">1000"],                         # integers stay integers, 0 stays 0
        ["Organism", "Homo sapiens", "0", "1", "7", ">1000", "—", "2.00"],            # 999.999 rounds to 1000.00: shown as >1000
        ["Disease", "missing value", "1", "1", "3", "0.00", "3.00", "—"]]
    assert tmp_path.joinpath("r.tsv").read_text().count("1.23456") == 1                            # the TSV is untouched


def test_the_page_keeps_search_and_hit_set_when_they_vary_and_names_the_properties_like_where_it_occurs(tmp_path):
    names = {"SampleType": "Sample type", "SampleTypeSub1": "Sample subtype", "NCBITaxonomy": "Organism", "UBERONBodyPartName": "Body part",
             "DOIDCommonName": "Disease", "HealthStatus": "Health status"}
    assert list(names) == g.SUMMARY_COLS and gnps.REDU_PROPERTY == names
    rows = [{"column": c, "value": "v", "n_hit_files": 1, "n_hit_datasets": 1, "bg_datasets": 1, "file_enrich": "2", "ds_enrich": "2", "dsw_enrich": "2"} for c in names]
    rows += [{**rows[0], "hit_set": "diag_pass"}, {**rows[0], "search": "masst_other"}]
    block = gnps.redu_preview(_tsv(tmp_path / "r.tsv", rows))
    assert block["columns"][:4] == ["Search", "Hit set", "Property", "Value"]
    assert [r[2] for r in block["rows"][:6]] == list(names.values()) and [r[0] for r in block["rows"][-2:]] == ["repository", "other"]
    assert [r[1] for r in block["rows"][-2:]] == ["Diagnostic pass", "All"]                        # the words of the Diagnostic verification table
    assert list(gnps.REDU_HIT_SET.values()) == ["All", "Verified", "Diagnostic pass"]
    (v,) = [m for m in spec.MODULES if m["id"] == "gnps_fastmasst"]
    assert "Verified" in str(v) and "Diagnostic pass" in str(v)


def test_the_cli_top_values_are_the_most_frequent_not_the_most_enriched():
    rows = [{"column": "SampleType", "value": "rare-but-enriched", "n_hit_datasets": 2, "n_hit_files": 3, "dsw_enrich": 50.0},
            {"column": "SampleType", "value": "common", "n_hit_datasets": 9, "n_hit_files": 20, "dsw_enrich": 1.2},
            {"column": "SampleType", "value": "common-fewer-files", "n_hit_datasets": 9, "n_hit_files": 12, "dsw_enrich": 1.5},
            {"column": "SampleType", "value": "missing value", "n_hit_datasets": 5, "n_hit_files": 6, "dsw_enrich": None}]
    assert [r["value"] for r in g.redu_most_frequent(rows, 3)] == ["common", "common-fewer-files", "missing value"]   # a missing value ranks by its counts
    tsv = [{k: str(v) for k, v in r.items()} for r in rows]                                        # the batch reads the TSV: strings
    assert [r["value"] for r in g.redu_most_frequent(tsv, 2)] == ["common", "common-fewer-files"]
    assert [r["value"] for r in rows][0] == "rare-but-enriched"                                    # the table's own order is untouched


def test_the_preview_cuts_after_the_order_of_the_tsv_and_an_empty_file_has_no_columns(tmp_path):
    rows = [{"column": "SampleType", "value": f"v{i}", "n_hit_files": 1, "n_hit_datasets": 1, "bg_datasets": 1} for i in range(gnps.PREVIEW_ROWS + 5)]
    block = gnps.redu_preview(_tsv(tmp_path / "r.tsv", rows))
    assert len(block["rows"]) == gnps.PREVIEW_ROWS and block["count"] == gnps.PREVIEW_ROWS + 5 and block["rows"][0][1] == "v0"
    assert gnps.redu_preview(_tsv(tmp_path / "e.tsv", [])) == {"rows": [], "count": 0}
    assert gnps.redu_preview(str(tmp_path / "none.tsv")) == {"rows": [], "count": 0}


def test_the_spec_orders_the_blocks_and_describes_the_enrichment_once():
    (m,) = [m for m in spec.MODULES if m["id"] == "gnps_fastmasst"]
    ids = [b["id"] for b in m["outputs"]]
    assert ids == ["searches", "distribution", "redu", "hits", "reference_hits", "suspect_hits", "verification", "downloads"]
    blocks = {b["id"]: b for b in m["outputs"]}
    redu = blocks["redu"]
    assert redu["columns"][2] == "Property" and redu["columns"][6] == "Background datasets" and redu["columns"][5] == "Hit datasets"
    assert "at least 2 datasets" in redu["help"] and set(redu["column_help"]) <= set(redu["columns"]) and len(redu["column_help"]) == 4
    assert "Coverage is partial" not in str(blocks["distribution"]["help"])
    assert "help" not in blocks["hits"]


def test_zero_is_shown_as_zero_and_only_a_missing_value_is_a_dash():
    assert [gnps._text(v) for v in (0, 0.0, "0", None, "", "x", 12)] == ["0", "0.0", "0", "—", "—", "x", "12"]
    cells = gnps._hit_cells({"USI": "u", "Cosine": "0.9", "Matching Peaks": "0", "CompoundName": "", "Adduct": "[M+H]+"})
    assert cells[1:5] == ["0.900", "0", "—", "[M+H]+"]
