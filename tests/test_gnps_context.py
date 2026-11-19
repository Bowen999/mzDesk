"""Spectrum in Which Set: hits by dataset, "Where it occurs", the bundled domain-MASST table and the GNPS2 Explorer titles. All offline."""
import csv
import gzip
import io
import json
import time

import pytest
from fastapi.testclient import TestClient

from mzdesk import spec
from mzdesk.api.app import app
from mzdesk.gnps import compound_db_check as g, domain_masst_build as build, hit_context as hc, redu_web as web_redu

client = TestClient(app)
REAL_FETCH = hc._fetch_explorer                 # tests/conftest.py replaces hc._fetch_explorer with a stub for every test


def hit(dataset, path, cosine=0.9, search="masst_repository", level="", **extra):
    """One passing hit row as hc.read_hits returns it."""
    row = {c: "" for c in hc.HIT_COLUMNS}
    row.update(search=search, library=g.REPO_LIB, USI=f"mzspec:{dataset}:{path}:scan:1", Dataset=dataset, usi_collection=dataset,
               Cosine=str(cosine), **{"Matching Peaks": "5"}, redu_match=level)
    row.update(extra)
    row["file"] = hc.file_id(row)
    return row


def redu(row, sample=None, taxon=None, part=None, disease=None):
    row.update(redu_match="exact", redu_SampleType=sample or "missing value", redu_NCBITaxonomy=taxon or "missing value",
               redu_UBERONBodyPartName=part or "missing value", redu_DOIDCommonName=disease or "missing value")
    return row


# ---------------------------------------------------------------- files
def test_a_file_is_counted_once_however_its_usi_is_written():
    a = hit("MSV000082718", "ccms_peak/peak/ST000666/EX1/Raw data/POS/S1_Rat.mzXML")
    b = hit("MSV000082718", "ccms_peak/ST000666/EX1/Raw data/POS/S1_Rat.mzXML", cosine=0.8)
    c = hit("MSV000082718", "ccms_peak/peak/ST000666/EX1/Raw data/POS/s2_rat.MZML")
    assert a["file"] == b["file"] == ("MSV000082718", "s1_rat") and c["file"] == ("MSV000082718", "s2_rat")
    (d,) = hc.dataset_rows([a, b, c], joined=False)
    assert (d["hits"], d["files"]) == (3, 2) and [r["Cosine"] for r in d["rows"]] == ["0.9", "0.9", "0.8"]


def test_the_masst_tables_are_joined_by_dataset_and_file_name_without_folder_or_extension(monkeypatch):
    monkeypatch.setattr(hc, "domain_table", lambda: {("MSV000078770", "vern_tomem-1408_ms2-pos_1-3_01_1196"): (("plant", "Vernonia"),)})
    r = hit("MSV000078770", "ccms_peak/VERN_TOMEM-1408_MS2/VERN_TOMEM-1408_MS2-POS_1-3_01_1196.mzXML")
    assert hc.domain_text(hc.domain_of(r)) == "plant: Vernonia"
    assert hc.domain_of(hit("MSV000078770", "ccms_peak/other.mzXML")) == () and hc.domain_of({"USI": "mzspec:GNPS:GNPS-LIBRARY:accession:X"}) == ()


def test_the_bundled_table_matches_its_sidecar_and_the_build_keys(monkeypatch):
    meta = json.loads((hc.DOMAIN_FILE.parent / "domain_masst.json").read_text(encoding="utf-8"))
    assert len(meta["sha"]) == 40 and meta["sha"] == build.SHA and meta["licence"].startswith("MIT")
    with gzip.open(hc.DOMAIN_FILE, "rt", encoding="utf-8", newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    assert rows[0] == build.COLUMNS and len(rows) - 1 == meta["rows"] == sum(v["rows"] for v in meta["domains"].values())
    assert {r[2] for r in rows[1:]} == set(build.TABLES) and all(hc.ACCESSION.match(r[0]) and r[1] == r[1].lower() and r[3] for r in rows[1:])
    dataset, key, domain, category = rows[1]
    assert hc.domain_table()[(dataset, key)][0] in {(r[2], r[3]) for r in rows[1:] if r[:2] == [dataset, key]}
    assert build.file_key("mzspec:MSV000078547:112111_ES129.mzXML") == ("MSV000078547", "112111_es129")
    assert build.file_key("mzspec:MSV1:a/b/C.mzML.gz") == ("MSV1", "c") and build.file_key("not a usi") is None


def test_the_build_drops_blank_and_qc_files_and_keeps_food_names(tmp_path):
    def table(name, header, rows):
        with open(tmp_path / name, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(header)
            w.writerows(rows)
    table("microbe_masst_table.csv", ["Taxaname_file", "Blank", "QC", "file_usi"],
          [["Bacteria", "No", "No", "mzspec:MSV1:a.mzML"], ["Bacteria", "Yes", "No", "mzspec:MSV1:blank.mzML"], ["Fungi", "No", "Yes", "mzspec:MSV1:qc.mzML"],
           ["Bacteria", "No", "No", "mzspec:MSV1:a.mzML"]])
    table("plant_masst_table.csv", ["Taxaname_file", "Blank", "QC", "file_usi"], [["Vernonia", "No", "No", "mzspec:MSV2:P1"]])
    table("tissue_masst_table.csv", ["ID", "file_usi"], [["Mus musculus_feces_healthy_male_unknown", "mzspec:MSV3:T1"], ["blank_qc_blank_qc", "mzspec:MSV3:T2"]])
    table("food_masst_table.csv", ["node_id", "file_usi"], [["coffee", "mzspec:MSV4:F1"], ["", "mzspec:MSV4:F2"]])
    table("personalCareProduct_masst_table.csv", ["Category", "file_usi"], [["Root > Personal Care > Shampoo", "mzspec:MSV5:C1"]])
    meta = build.build(str(tmp_path), str(tmp_path / "out"))
    with gzip.open(tmp_path / "out" / "domain_masst.tsv.gz", "rt", encoding="utf-8", newline="") as f:
        got = [tuple(r) for r in csv.reader(f, delimiter="\t")][1:]
    assert sorted(got) == sorted([("MSV1", "a", "microbe", "Bacteria"), ("MSV2", "p1", "plant", "Vernonia"),
                                  ("MSV3", "t1", "tissue", "Mus musculus_feces_healthy_male_unknown"), ("MSV4", "f1", "food", "coffee"),
                                  ("MSV5", "c1", "personal_care", "Personal Care > Shampoo")])
    assert meta["rows"] == 5 and meta["domains"]["microbe"]["blank_or_qc_dropped"] == 2


# ---------------------------------------------------------------- groups
def test_groups_are_ordered_by_hits_then_best_cosine_and_the_library_comes_last():
    rows = [hit("MSV2", f"a{i}.mzML", cosine=0.8) for i in range(2)] + [hit("MSV1", "b.mzML", cosine=0.99)]
    rows += [hit("MSV3", "c1.mzML", cosine=0.95), hit("MSV3", "c2.mzML", cosine=0.7)]
    rows += [dict(hit("GNPS Library", "x"), search="masst_library", library=g.REF_LIB, USI="mzspec:GNPS:GNPS-LIBRARY:accession:A", Dataset="GNPS Library",
                  CompoundName=f"n{i}", usi_collection="GNPS") for i in range(5)]
    for r in rows:
        r["file"] = hc.file_id(r)
    datasets = hc.dataset_rows(rows, joined=False)
    assert [(d["dataset"], d["hits"]) for d in datasets] == [("MSV3", 2), ("MSV2", 2), ("MSV1", 1), ("", 5)]
    assert [d["best"] for d in datasets[:2]] == [0.95, 0.8]
    assert [r["Cosine"] for r in datasets[0]["rows"]] == ["0.95", "0.7"]
    block = hc.hits_block(datasets, {"MSV3": {"title": "x" * 200}}, False, len(rows), lambda r: [r["USI"]])
    assert [g_["title"] for g_ in block["groups"]] == ["MSV3 · " + "x" * 89 + "…", "MSV2", "MSV1", "GNPS reference library"]
    assert all("url" not in x for x in block["groups"])                          # "MSV3" is not a real accession; the library has no page
    assert block["count"] == len(rows) and "help" not in block
    assert block["groups"][-1]["meta"] == "5 hits · 5 names · best cosine 0.900"


def test_groups_and_rows_are_capped_and_the_block_says_so(monkeypatch):
    monkeypatch.setattr(hc, "MAX_GROUPS", 3)
    monkeypatch.setattr(hc, "MAX_ROWS", 2)
    rows = [hit(f"MSV{i:09d}", f"f{j}.mzML") for i in range(5) for j in range(3 if i == 0 else 1)]
    block = hc.hits_block(hc.dataset_rows(rows, False), {}, False, len(rows), lambda r: [r["USI"]])
    assert len(block["groups"]) == 3 and len(block["groups"][0]["rows"]) == 2 and block["count"] == 7
    assert "first 3 of 5 datasets" in block["note"] and "2 best hits of each dataset" in block["note"] and "help" not in block
    assert block["groups"][0]["url"] == "https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000000000"


def test_group_meta_says_what_pan_redu_knows_and_falls_back_to_explorer_species():
    matched = redu(hit("MSV000000001", "a.mzML"), "plant", "9606|Homo sapiens", "liver")
    matched2 = redu(hit("MSV000000001", "b.mzML", cosine=0.8), "plant", "missing value")
    absent = hit("MSV000000002", "c.mzML", level="no_dataset")
    partial = hit("MSV000000003", "d.mzML", level="dataset_only")
    no_organism = redu(hit("MSV000000004", "e.mzML"), "animal")
    datasets = {d["dataset"]: d for d in hc.dataset_rows([matched, matched2, absent, partial, no_organism], joined=True)}
    info = {a: {"species": "Rattus norvegicus"} for a in datasets}
    meta = {a: hc._meta(d, info, True) for a, d in datasets.items()}
    assert meta["MSV000000001"] == "2 hits · 2 files · best cosine 0.900 · plant · Homo sapiens · liver"
    assert meta["MSV000000002"] == "1 hit · 1 file · best cosine 0.900 · Rattus norvegicus · not in Pan-ReDU"
    assert meta["MSV000000003"] == "1 hit · 1 file · best cosine 0.900 · Rattus norvegicus · no Pan-ReDU file match"
    assert meta["MSV000000004"] == "1 hit · 1 file · best cosine 0.900 · animal · Rattus norvegicus"
    assert hc._meta(datasets["MSV000000002"], {}, False) == "1 hit · 1 file · best cosine 0.900"            # not joined: nothing about Pan-ReDU


def test_dataset_links_are_the_repositories_own_pages():
    assert hc.dataset_url("msv000087728") == "https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000087728"
    assert hc.dataset_url("ST000666").startswith("https://www.metabolomicsworkbench.org/") and hc.dataset_url("MTBLS1").endswith("/metabolights/MTBLS1")
    assert hc.dataset_url("GNPS Library") is None and hc.dataset_url("") is None


def test_tissue_ids_are_shown_as_readable_labels(monkeypatch):
    assert hc.tissue_label("Mus musculus_feces_unknown_unknown_unknown") == "Mus musculus · feces"
    assert hc.tissue_label("Mus musculus_caecum_healthy_male_unknown") == "Mus musculus · caecum · healthy · male"
    assert hc.tissue_label("Homo sapiens_brain_unknown_female_Later Adulthood (>65 yrs)") == "Homo sapiens · brain · female · Later Adulthood (>65 yrs)"
    assert hc.tissue_label("Homo sapiens_unknown_Alzheimer's disease_missing value_") == "Homo sapiens · Alzheimer's disease"
    assert hc.tissue_label("Homo sapiens_unknown_unknown_unknown_unknown") == "Homo sapiens"
    assert hc.tissue_label("not a tissue id") == "not a tissue id" and hc.tissue_label("a_b_c") == "a_b_c"      # not five parts: unchanged
    assert hc.tissue_label("Mus musculus_feces_unknown_unknown_unknown") != hc.tissue_label("Mus musculus_feces_healthy_unknown_unknown")
    ids = ["Mus musculus_feces_unknown_unknown_unknown", "Mus musculus_feces_healthy_male_unknown"]
    monkeypatch.setattr(hc, "domain_table", lambda: {("MSV000000001", "a"): (("tissue", ids[0]), ("plant", "Vernonia")), ("MSV000000001", "b"): (("tissue", ids[1]),),
                                                     ("MSV000000002", "c"): (("tissue", ids[0]),)})
    r = hit("MSV000000001", "a.mzML")
    assert hc.domain_of(r) == (("tissue", "Mus musculus · feces"), ("plant", "Vernonia"))                  # only tissue IDs are reformatted
    assert hc.domain_text(hc.domain_of(r)) == "tissue: Mus musculus · feces; plant: Vernonia"
    rows = [hit("MSV000000001", "a.mzML"), hit("MSV000000001", "b.mzML"), hit("MSV000000002", "c.mzML")]
    groups, table, _ = hc.distribution(rows, False, {})
    (tissue,) = [x for x in groups if x["title"] == "Tissue · tissueMASST"]
    assert tissue["rows"] == [["Mus musculus · feces", "2", "2", "67%"], ["Mus musculus · feces · healthy · male", "1", "1", "33%"]]
    assert tissue["meta"] == "Mus musculus · feces 2 · Mus musculus · feces · healthy · male 1"
    assert sum(int(x[1]) for x in tissue["rows"]) == 3          # the same three files as with the raw IDs


# ---------------------------------------------------------------- where it occurs
def test_distribution_counts_files_and_datasets_and_leaves_out_missing_values(monkeypatch):
    monkeypatch.setattr(hc, "domain_table", lambda: {("MSV000000001", "a"): (("plant", "Vernonia"),), ("MSV000000001", "b"): (("plant", "Vernonia"),),
                                                      ("MSV000000002", "c"): (("food", "coffee"), ("plant", "Coffea"))})
    rows = [redu(hit("MSV000000001", "x/a.mzML"), "plant", "4|Vernonia", "leaf"),
            redu(hit("MSV000000001", "y/a.mzML", cosine=0.8), "plant", "4|Vernonia", "leaf"),          # the same file again: still one file
            redu(hit("MSV000000001", "b.mzML"), "plant", "missing value", "leaf", "healthy"),
            redu(hit("MSV000000002", "c.mzML"), "animal", "9606|Homo sapiens"),
            hit("MSV000000003", "d.mzML", level="no_dataset"),
            dict(hit("GNPS Library", "z"), search="masst_library", library=g.REF_LIB, redu_SampleType="plant")]
    info = {"MSV000000001": {"species": "Viridiplantae"}, "MSV000000002": {"species": "Viridiplantae"}, "MSV000000003": {"species": None}}
    groups, table, note = hc.distribution(rows, True, info)
    by = {x["title"]: x for x in groups}
    assert list(by) == ["Sample type · Pan-ReDU", "Organism · Pan-ReDU", "Body part · Pan-ReDU", "Disease · Pan-ReDU", "Dataset species · GNPS2 Explorer",
                        "Plant · plantMASST", "Food · foodMASST"]                                              # empty groups are skipped
    assert by["Sample type · Pan-ReDU"]["rows"] == [["plant", "2", "1", "50%"], ["animal", "1", "1", "25%"]]          # 4 hit files; the library hit is not a file
    assert by["Organism · Pan-ReDU"]["rows"] == [["Homo sapiens", "1", "1", "25%"], ["Vernonia", "1", "1", "25%"]]    # "4|" is not shown; "missing value" is not counted
    assert by["Body part · Pan-ReDU"]["rows"] == [["leaf", "2", "1", "50%"]] and by["Disease · Pan-ReDU"]["rows"] == [["healthy", "1", "1", "25%"]]
    assert by["Sample type · Pan-ReDU"]["meta"] == "plant 2 · animal 1"                                          # the coverage is in the note, not in every group
    assert by["Dataset species · GNPS2 Explorer"]["rows"] == [["Viridiplantae", "3", "2", "75%"]]
    assert by["Dataset species · GNPS2 Explorer"]["meta"] == "Viridiplantae 2"
    assert by["Plant · plantMASST"]["rows"][0] == ["Vernonia", "2", "1", "50%"] and by["Plant · plantMASST"]["meta"] == "Vernonia 2 · Coffea 1"
    assert note == ("Coverage: Pan-ReDU 3 of 4 hit files · plantMASST 3 of 4 · foodMASST 1 of 4 · GNPS2 Explorer 2 of 3 datasets"
                    " · none in microbeMASST, tissueMASST, personalCareProductMASST")
    assert by["Food · foodMASST"]["rows"] == [["coffee", "1", "1", "25%"]]
    assert ["Sample type · Pan-ReDU", "plant", 2, 1, 0.5] in table and all(len(r) == 5 for r in table)
    assert not any("missing" in row[0] for x in groups for row in x["rows"])


def test_distribution_without_the_join_has_no_pan_redu_groups_and_without_hits_no_groups():
    groups, _, note = hc.distribution([redu(hit("MSV000000001", "a.mzML"), "plant")], False, {"MSV000000001": {"species": "Viridiplantae"}})
    assert [x["title"] for x in groups] == ["Dataset species · GNPS2 Explorer"]
    assert note.startswith("Coverage: GNPS2 Explorer 1 of 1 dataset · none in ") and "Pan-ReDU" not in note
    assert hc.distribution([], True, {}) == ([], [], "")


def test_a_small_share_is_not_rounded_to_zero():
    assert hc._share(1, 1000) == "<1%" and hc._share(1, 3) == "33%" and hc._share(0, 0) == "—"


# ---------------------------------------------------------------- GNPS2 Explorer
def test_explorer_titles_are_cached_for_30_days_and_asked_once(monkeypatch):
    calls = []
    monkeypatch.setattr(hc, "_fetch_explorer", lambda acc: calls.append(acc) or {"title": f"T {acc}", "species": "S", "repository": "MassIVE"})
    warnings = []
    now = 1_000_000.0
    got = hc.explorer_info(["MSV000000001", "MSV000000002", "MSV000000001", "GNPS Library", ""], warnings, now=now)
    assert sorted(calls) == ["MSV000000001", "MSV000000002"] and got["MSV000000002"]["title"] == "T MSV000000002" and not warnings
    assert hc.explorer_info(["MSV000000001"], warnings, now=now + 29 * 86400)["MSV000000001"]["title"] == "T MSV000000001" and len(calls) == 2
    hc.explorer_info(["MSV000000001"], warnings, now=now + 31 * 86400)
    assert len(calls) == 3                                                    # older than 30 days: asked again
    assert (hc.cache_dir() / "explorer_datasets.json").is_file()


def test_explorer_failures_give_no_title_and_one_warning_never_an_error(monkeypatch):
    monkeypatch.setattr(hc, "_fetch_explorer", lambda acc: None if acc.endswith("2") else {"title": "ok", "species": None, "repository": None})
    warnings = []
    got = hc.explorer_info(["MSV000000001", "MSV000000002", "MSV000000003"], warnings)
    assert sorted(got) == ["MSV000000001", "MSV000000003"] and len(warnings) == 1 and "1 of 3" in warnings[0]
    assert "MSV000000002" not in hc.explorer_info(["MSV000000002"], [])                      # a failure is not cached
    assert hc._title("MSV000000002", got) == "MSV000000002" and hc._title("MSV000000001", got) == "MSV000000001 · ok"


def test_explorer_down_stops_waiting_and_a_broken_cache_file_is_ignored(monkeypatch):
    calls = []
    monkeypatch.setattr(hc, "_fetch_explorer", lambda acc: calls.append(acc))
    (hc.cache_dir()).mkdir(parents=True)
    (hc.cache_dir() / "explorer_datasets.json").write_text("{not json", encoding="utf-8")
    warnings = []
    assert hc.explorer_info([f"MSV{i:09d}" for i in range(60)], warnings) == {}
    assert len(calls) < 60 and len(warnings) == 1 and "of 60" in warnings[0]


def test_explorer_is_asked_for_as_many_datasets_as_the_page_shows_and_only_accessions(monkeypatch):
    calls, warnings = [], []
    monkeypatch.setattr(hc, "_fetch_explorer", lambda acc: calls.append(acc) or {"title": acc, "species": None, "repository": None})
    hc.explorer_info([f"MSV{i:09d}" for i in range(250)] + ["not an accession", "GNPS"], warnings)
    assert len(warnings) == 1 and "50 of 250" in warnings[0]                                                    # the 50 over the cap are said, not silent
    assert hc.EXPLORER_MAX == hc.MAX_GROUPS == 200 and len(calls) == 200 and all(c.startswith("MSV") for c in calls)
    assert "MSV000000199" in calls and "MSV000000200" not in calls


def test_explorer_answer_is_read_for_title_species_and_repository(monkeypatch):
    class Reply(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False
    body = {"title": "  A  title\n", "species": "Viridiplantae (NCBITaxon:33090)", "repository": "MassIVE", "pi": "someone", "raw": {}}
    seen = {}

    def urlopen(req, timeout):
        seen.update(url=req.full_url, timeout=timeout)
        return Reply(json.dumps(body).encode())
    monkeypatch.setattr(hc.urllib.request, "urlopen", urlopen)
    assert REAL_FETCH("MSV000087728") == {"title": "A title", "species": "Viridiplantae", "repository": "MassIVE"}
    assert seen == {"url": "https://explorer.gnps2.org/api/datasets/MSV000087728", "timeout": 8}
    body["species"] = "no species"
    assert REAL_FETCH("MSV000087728")["species"] is None
    monkeypatch.setattr(hc.urllib.request, "urlopen", lambda *a, **k: (_ for _ in ()).throw(OSError("offline")))
    assert REAL_FETCH("MSV000087728") is None


def _http_error(code):
    return hc.urllib.error.HTTPError("https://explorer.gnps2.org/x", code, "x", {}, None)


def test_explorer_fetch_tells_a_missing_dataset_from_an_outage(monkeypatch):
    def answer(result):
        def urlopen(req, timeout):
            if isinstance(result, Exception):
                raise result
            return type("R", (io.BytesIO,), {"__enter__": lambda self: self, "__exit__": lambda self, *a: False})(json.dumps(result).encode())
        monkeypatch.setattr(hc.urllib.request, "urlopen", urlopen)
        return REAL_FETCH("MSV000087728")
    assert answer(_http_error(404)) == hc.NOT_FOUND and answer({"pi": "x"}) == hc.NOT_FOUND and answer(["x"]) == hc.NOT_FOUND     # Explorer answered: it does not know
    assert answer({"title": "", "species": "no species"}) == hc.NOT_FOUND
    assert answer(_http_error(500)) is None and answer(_http_error(429)) is None and answer(TimeoutError()) is None and answer(OSError("x")) is None


def test_explorer_reads_at_most_2_mb_of_an_answer(monkeypatch):
    sizes = []

    class Reply(io.BytesIO):
        def read(self, n=-1):
            sizes.append(n)
            return super().read(n)

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False
    monkeypatch.setattr(hc.urllib.request, "urlopen", lambda req, timeout: Reply(b'{"title": "t"}'))
    assert REAL_FETCH("MSV000087728")["title"] == "t" and sizes == [2_000_000]


def test_ordinary_misses_do_not_trip_the_breaker_and_are_cached_for_a_week(monkeypatch):
    calls = []
    monkeypatch.setattr(hc, "_fetch_explorer", lambda acc: calls.append(acc) or (hc.NOT_FOUND if int(acc[3:]) < 20 else
                                                                                 {"title": "T", "species": None, "repository": None}))
    accessions = [f"MSV{i:09d}" for i in range(30)]
    warnings = []
    got = hc.explorer_info(accessions, warnings, now=1_000_000.0)
    assert len(calls) == 30 and sorted(got) == accessions[20:] and warnings == []             # 20 unknown datasets, no outage: all asked, no warning
    hc.explorer_info(accessions, warnings, now=1_000_000.0 + 6 * 86400)
    assert len(calls) == 30                                                                  # misses and answers come from the cache
    hc.explorer_info(accessions, warnings, now=1_000_000.0 + 8 * 86400)
    assert len(calls) == 50                                                                  # a miss is asked again after 7 days (the 10 answers are kept 30)


def test_outage_errors_trip_the_breaker_but_not_after_an_answer(monkeypatch):
    calls = []
    monkeypatch.setattr(hc, "_fetch_explorer", lambda acc: calls.append(acc))                # timeouts
    warnings = []
    assert hc.explorer_info([f"MSV{i:09d}" for i in range(80)], warnings) == {} and len(calls) < 80 and "80 of 80" in warnings[0]
    assert not (hc.cache_dir() / "explorer_datasets.json").exists() or json.loads((hc.cache_dir() / "explorer_datasets.json").read_text()) == {}
    calls.clear()                                                                            # the server answers once, then times out: no breaker
    monkeypatch.setattr(hc, "EXPLORER_WORKERS", 1)
    monkeypatch.setattr(hc, "_fetch_explorer", lambda acc: calls.append(acc) or (hc.NOT_FOUND if not calls[1:] else None))
    hc.explorer_info([f"MSV{i:09d}" for i in range(100, 120)], [])
    assert len(calls) == 20


def test_explorer_has_a_wall_clock_budget(monkeypatch):
    monkeypatch.setattr(hc, "EXPLORER_BUDGET", 0.3)
    monkeypatch.setattr(hc, "_fetch_explorer", lambda acc: time.sleep(0.4) or {"title": "T", "species": None, "repository": None})
    warnings = []
    t = time.monotonic()
    got = hc.explorer_info([f"MSV{i:09d}" for i in range(60)], warnings)
    assert time.monotonic() - t < 1.5 and got == {}                                          # 60 x 0.4 s over 6 workers would be 4 s
    assert len(warnings) == 1 and "60 of 60" in warnings[0]


def test_a_corrupt_cache_entry_is_stale_never_an_error(monkeypatch):
    calls = []
    monkeypatch.setattr(hc, "_fetch_explorer", lambda acc: calls.append(acc) or {"title": "fresh", "species": None, "repository": None})
    now = 5_000_000.0
    good = {"t": now - 86400, "info": {"title": "kept", "species": "S", "repository": None}}
    bad = {"MSV000000001": {"t": "yesterday", "info": {"title": "x"}}, "MSV000000002": {"info": {"title": "x"}}, "MSV000000003": [1],
           "MSV000000004": {"t": now, "info": "text"}, "MSV000000005": {"t": None, "miss": True}, "MSV000000006": {"t": now, "info": {"title": 7}},
           "MSV000000007": {"t": now - 3 * 30 * 86400, "info": {"title": "old"}}, "MSV000000008": good}
    (hc.cache_dir()).mkdir(parents=True)
    (hc.cache_dir() / "explorer_datasets.json").write_text(json.dumps({**bad, "MSV000000010": {"t": now - 61 * 86400, "info": {"title": "very old"}}}), encoding="utf-8")
    got = hc.explorer_info(list(bad), [], now=now)
    assert sorted(calls) == [f"MSV00000000{i}" for i in range(1, 8)] and got["MSV000000008"]["title"] == "kept"
    assert got["MSV000000001"]["title"] == "fresh" and got["MSV000000006"]["title"] == "fresh"
    kept = json.loads((hc.cache_dir() / "explorer_datasets.json").read_text(encoding="utf-8"))
    assert "MSV000000010" not in kept and kept["MSV000000003"]["info"]["title"] == "fresh" and kept["MSV000000008"] == good      # stale entries are replaced or pruned


def test_a_failed_cache_write_leaves_no_temp_file_and_no_error(monkeypatch):
    monkeypatch.setattr(hc.json, "dump", lambda *a, **k: (_ for _ in ()).throw(TypeError("not serialisable")))
    hc._write_cache(hc.cache_dir() / "explorer_datasets.json", {"MSV000000001": {"t": 1.0, "info": {}}}, 2.0)
    assert list(hc.cache_dir().glob("*.tmp")) == [] and not (hc.cache_dir() / "explorer_datasets.json").exists()


def test_species_are_split_cleaned_and_non_taxa_dropped():
    assert hc.species_list("Arabidopsis thaliana (NCBITaxon:3702)::reference compound::Triticum aestivum::Arabidopsis thaliana") == \
        ["Arabidopsis thaliana", "Triticum aestivum"]
    assert hc.species_list("Data from multiple bacteria") == [] and hc.species_list("no species; Not applicable ; blank;;") == []
    assert hc.species_list("Homo sapiens; Mus musculus (NCBITaxon:10090)") == ["Homo sapiens", "Mus musculus"] and hc.species_list(None) == []
    info = {"MSV000000001": {"species": "Arabidopsis thaliana::reference compound::Triticum aestivum"},
            "MSV000000002": {"species": "Triticum aestivum;Zea mays"}, "MSV000000003": {"species": "Data from multiple bacteria"}}
    rows = [hit("MSV000000001", "a.mzML"), hit("MSV000000002", "b.mzML"), hit("MSV000000003", "c.mzML")]
    groups, table, note = hc.distribution(rows, False, info)
    (sp,) = groups
    assert sp["rows"] == [["Triticum aestivum", "2", "2", "67%"], ["Arabidopsis thaliana", "1", "1", "33%"], ["Zea mays", "1", "1", "33%"]]
    assert "reference" not in str(table) and "bacteria" not in str(table) and "GNPS2 Explorer 2 of 3 datasets" in note
    many = {"species": "A a::B b::C c::D d::E e"}
    d = hc.dataset_rows([hit("MSV000000009", "x.mzML")], True)[0]
    assert hc._meta(d, {"MSV000000009": many}, False) == "1 hit · 1 file · best cosine 0.900 · A a, B b, C c +2"       # the fallback of a group Pan-ReDU has no organism for


def test_hits_without_a_dataset_form_a_group_titled_no_dataset():
    r = hit("", "x.mzML")
    r.update(USI="not a usi", Dataset="", usi_collection="")
    r["file"] = hc.file_id(r)
    block = hc.hits_block(hc.dataset_rows([r], False), {}, False, 1, lambda row: [row["USI"]])
    (group,) = block["groups"]
    assert group["title"] == hc.NO_DATASET == "(no dataset)" and "url" not in group and group["meta"] == "1 hit · 0 files · best cosine 0.900"


def test_counts_are_singular_for_one():
    assert hc.count(1, "hit") == "1 hit" and hc.count(0, "file") == "0 files" and hc.count(2, "name") == "2 names"
    lib = dict(hit("GNPS Library", "z"), search="masst_library", library=g.REF_LIB, CompoundName="n", USI="mzspec:GNPS:GNPS-LIBRARY:accession:A")
    lib["file"] = hc.file_id(lib)
    assert hc._meta(hc.dataset_rows([lib], False)[0], {}, False) == "1 hit · 1 name · best cosine 0.900"
    _, _, note = hc.distribution([hit("MSV000000001", "a.mzML")], False, {})
    assert note.startswith("Coverage: GNPS2 Explorer 0 of 1 dataset")


# ---------------------------------------------------------------- the page
RAW = [{"USI": "mzspec:MSV000000001:path/a.mzML:scan:1", "Dataset": "MSV000000001", "Cosine": 0.9, "Matching Peaks": 4, "Delta Mass": 0.0},
       {"USI": "mzspec:MSV000000001:path/b.mzML:scan:2", "Dataset": "MSV000000001", "Cosine": 0.8, "Matching Peaks": 4, "Delta Mass": 0.0},
       {"USI": "mzspec:MSV000000003:x/d.mzML:scan:9", "Dataset": "MSV000000003", "Cosine": 0.75, "Matching Peaks": 3, "Delta Mass": 0.0},
       {"USI": "mzspec:MSV000000009:x/q.mzML:scan:9", "Dataset": "MSV000000009", "Cosine": 0.71, "Matching Peaks": 3, "Delta Mass": 0.0}]
LIBRARY = [{"USI": "mzspec:GNPS:GNPS-LIBRARY:accession:A", "Dataset": "GNPS Library", "Cosine": 0.88, "Matching Peaks": 5, "CompoundName": "Known", "Adduct": "[M+H]+",
            "GNPSLibraryAccession": "A"}]


@pytest.fixture
def fake_search(monkeypatch):
    def search(peaks, precursor, library, cosine, min_peaks, pm_tol, frag_tol, **kw):
        raw = RAW if library == g.REPO_LIB else LIBRARY
        return g.summarise_raw(raw, library, min_peaks, "t"), {"task_id": "t", "raw": raw, "elapsed_s": 0, "n_poll_errors": 0}
    monkeypatch.setattr(g, "fasst_search", search)


def _run(**values):
    r = client.post("/api/modules/gnps_fastmasst/run", json={"spectrum": {"source": "Paste", "value": "100 100\n110 20"}, "precursor_mz": "300", **values})
    assert r.status_code == 200, r.text
    return r.json()


def test_page_lists_hits_by_dataset_with_the_library_apart_and_every_download(fake_search, monkeypatch):
    monkeypatch.setattr(hc, "domain_table", lambda: {("MSV000000001", "a"): (("microbe", "Bacteria"),)})
    body = _run(redu_join=False)
    out = body["outputs"]
    titles = [x["title"] for x in out["hits"]["groups"]]
    assert titles == ["MSV000000001 · Title of MSV000000001", "MSV000000003 · Title of MSV000000003", "MSV000000009 · Title of MSV000000009", "GNPS reference library"]
    first = out["hits"]["groups"][0]
    assert first["meta"] == "2 hits · 2 files · best cosine 0.900 · Testus exemplaris" and first["url"].endswith("accession=MSV000000001")
    assert first["rows"][0][:3] == ["mzspec:MSV000000001:path/a.mzML:scan:1", "0.900", "4"] and first["rows"][0][5:9] == ["—", "—", "—", "microbe: Bacteria"]
    assert out["hits"]["groups"][-1]["rows"][0][3:5] == ["Known", "[M+H]+"]
    assert [x["title"] for x in out["distribution"]["groups"]] == ["Dataset species · GNPS2 Explorer", "Microbe · microbeMASST"]
    assert not any("Pan-ReDU" in w for w in body["warnings"])                                                  # join off by choice: no warning
    files = {f["name"]: f["content"] for f in out["downloads"]["files"]}
    hits = list(csv.DictReader(io.StringIO(files["fastmasst.hits.tsv"]), delimiter="\t"))
    assert len(hits) == 5 and {r["dataset_title"] for r in hits if r["Dataset"] == "MSV000000001"} == {"Title of MSV000000001"}
    assert [r["domain_masst"] for r in hits if r["Dataset"] == "MSV000000001"] == ["microbe: Bacteria", ""]
    datasets = list(csv.DictReader(io.StringIO(files["fastmasst.datasets.tsv"]), delimiter="\t"))
    assert [d["dataset"] for d in datasets] == ["MSV000000001", "MSV000000003", "MSV000000009", "GNPS reference library"]
    assert datasets[0]["title"] == "Title of MSV000000001" and datasets[0]["hits"] == "2" and datasets[0]["best_cosine"] == "0.900" and datasets[0]["redu_status"] == "not joined"
    dist = list(csv.DictReader(io.StringIO(files["fastmasst.distribution.tsv"]), delimiter="\t"))
    assert {r["field"] for r in dist} == {"Dataset species · GNPS2 Explorer", "Microbe · microbeMASST"}


def test_explorer_down_keeps_the_page_with_accessions_and_one_warning(fake_search, monkeypatch):
    monkeypatch.setattr(hc, "_fetch_explorer", lambda acc: None)
    body = _run(redu_join=False)
    assert [x["title"] for x in body["outputs"]["hits"]["groups"]][:3] == ["MSV000000001", "MSV000000003", "MSV000000009"]
    assert sum("GNPS2 Explorer" in w for w in body["warnings"]) == 1


def test_missing_pan_redu_cache_is_a_warning_not_an_error(fake_search, tmp_path, monkeypatch):
    monkeypatch.setattr(g, "REDU_SLIM", str(tmp_path / "nowhere" / "redu_slim.tsv.gz"))
    monkeypatch.setattr(web_redu, "TargetedRedu", lambda *a, **k: (_ for _ in ()).throw(AssertionError("the join must not start")))
    body = _run(redu_join=True)
    warnings = [w for w in body["warnings"] if "Pan-ReDU" in w]
    assert len(warnings) == 1 and "mzdesk setup redu" in warnings[0]
    assert body["outputs"]["redu"]["rows"] == [] and body["outputs"]["hits"]["groups"]
    assert not any(x["title"].endswith("Pan-ReDU") for x in body["outputs"]["distribution"]["groups"])


def test_the_join_with_the_cache_adds_sample_columns_and_the_enrichment_block(fake_search, tmp_path, monkeypatch):
    path = tmp_path / "redu.tsv.gz"
    monkeypatch.setattr(g, "REDU_SLIM", str(path))
    row = {c: "missing value" for c in g.REDU_COLS}
    row.update(ATTRIBUTE_DatasetAccession="MSV000000001", DataSource="MassIVE", USI="mzspec:MSV000000001:other/a.mzML", SampleType="plant", SampleTypeSub1="x",
               NCBITaxonomy="3|Vernonia", UBERONBodyPartName="leaf", IonizationSourceAndPolarity="positive")
    with gzip.open(path, "wt", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(g.REDU_COLS)
        w.writerow([row[c] for c in g.REDU_COLS])
    body = _run(redu_join=True)
    groups = body["outputs"]["hits"]["groups"]
    assert groups[0]["meta"] == "2 hits · 2 files · best cosine 0.900 · plant · Vernonia · leaf"                    # a.mzML matched by file name, b.mzML did not
    assert "not in Pan-ReDU" in groups[1]["meta"]
    dist = {x["title"]: x for x in body["outputs"]["distribution"]["groups"]}
    assert dist["Organism · Pan-ReDU"]["rows"] == [["Vernonia", "1", "1", "25%"]] and dist["Body part · Pan-ReDU"]["meta"] == "leaf 1"
    assert body["outputs"]["redu"]["rows"]
    assert not any("unmatched" in w for w in body["warnings"])                                  # the unmatched fractions stay in fastmasst.json, not on the page
    assert "unmatched_frac_hits" in json.dumps(body["outputs"]["downloads"]["files"][0]["content"])
    note = body["outputs"]["distribution"]["note"]
    assert note.startswith("Coverage: Pan-ReDU 1 of 4 hit files") and "GNPS2 Explorer 3 of 3 datasets" in note
    assert all("covers" not in x["meta"] for x in body["outputs"]["distribution"]["groups"])


def test_the_distribution_note_is_sent_only_when_there_is_something_to_say(fake_search):
    body = _run(redu_join=False, library="library")
    assert body["outputs"]["distribution"] == {"groups": []}


def test_capped_groups_and_rows_are_a_note_under_the_title_not_a_warning(fake_search, monkeypatch):
    monkeypatch.setattr(hc, "MAX_GROUPS", 2)
    monkeypatch.setattr(hc, "MAX_ROWS", 1)
    body = _run(redu_join=False, library="repository")
    hits = body["outputs"]["hits"]
    assert len(hits["groups"]) == 2 and "first 2 of 3 datasets" in hits["note"] and "1 best hits of each dataset" in hits["note"]
    assert "help" not in hits and not [w for w in body["warnings"] if "of 3 datasets" in w]          # the cap is said once, as the note under the title


# ---------------------------------------------------------------- spec
def _module():
    return next(m for m in spec.MODULES if m["id"] == "gnps_fastmasst")


def test_spec_has_the_grouped_blocks_in_order_and_the_join_on_by_default():
    m = _module()
    ids = [b["id"] for b in m["outputs"]]
    assert ids[:3] == ["searches", "distribution", "redu"] and ids[-1] == "downloads" and "datasets" not in ids and ids.count("hits") == 1
    blocks = {b["id"]: b for b in m["outputs"]}
    assert blocks["hits"]["type"] == "grouptable" and blocks["hits"]["title"] == "Passing hits by dataset" and blocks["hits"]["page_size"] == 10
    assert blocks["hits"]["columns"] == ["USI", "Cosine", "Matched peaks", "Name", "Adduct", "Sample type", "Organism", "Body part", "Domain MASST",
                                         "Verification", "Diagnostic pass"]
    assert blocks["distribution"]["type"] == "grouptable" and blocks["distribution"]["title"] == "Where it occurs" and blocks["distribution"]["page_size"] == 10
    assert blocks["distribution"]["columns"] == ["Value", "Hit files", "Hit datasets", "Share of hit files"] and blocks["distribution"]["help"]
    assert next(f for f in m["inputs"] if f["id"] == "redu_join")["value"] is True
    assert "GNPS2 Explorer" in m["network"] and "accessions" in m["network"]


def test_grouptable_helper_carries_page_size_and_compact_only_when_asked():
    plain = spec.grouptable(["A"], "T", id="t")
    assert "page_size" not in plain and "compact" not in plain
    paged = spec.grouptable(["A"], "T", id="t", page_size=10, compact=True)
    assert paged["page_size"] == 10 and paged["compact"] is True and paged["type"] == "grouptable"
