"""GNPS2 migration: offline API, original/streaming parity and real file formats."""
import base64
import copy
import csv
import gzip
import io
import json
import subprocess
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from mzdesk.api.app import app
from mzdesk.api.runners import gnps as web
from mzdesk.gnps import compound_db_check as g, redu_web
from mzdesk.spec import MODULES


client = TestClient(app)
MGF = b"""BEGIN IONS
TITLE=first
SCANS=42
PEPMASS=300
CHARGE=1+
RTINSECONDS=120
IONMODE=Positive
100 100
110 20
END IONS
BEGIN IONS
TITLE=second
SCANS=44
PEPMASS=300
100.3 100
END IONS
"""
QUERY = "QUERY scaninfo(MS2DATA) WHERE MS2PROD=100:TOLERANCEMZ=0.01"


def _file(name, raw):
    return {"name": name, "base64": base64.b64encode(raw).decode()}


def _post(module, **values):
    return client.post(f"/api/modules/{module}/run", json=values)


@pytest.fixture
def small_redu(tmp_path, monkeypatch):
    monkeypatch.setenv("MZDESK_CACHE_DIR", str(tmp_path))
    path = tmp_path / "redu.tsv.gz"
    monkeypatch.setattr(g, "REDU_SLIM", str(path))
    monkeypatch.setattr(g, "BUNDLED_REDU", str(tmp_path / "missing.gz"))
    rows = []

    def add(ds, usi, tax, part, polarity, sample="animal", chromat="LC"):
        row = {c: "missing value" for c in g.REDU_COLS}
        row.update(ATTRIBUTE_DatasetAccession=ds, DataSource="MassIVE", USI=usi, SampleType=sample,
                   SampleTypeSub1="tissue", NCBITaxonomy=tax, UBERONBodyPartName=part,
                   IonizationSourceAndPolarity=polarity, ChromatographyAndPhase=chromat)
        rows.append([row[c] for c in g.REDU_COLS])

    add("MSV1", "mzspec:MSV1:path/a.mzML", "10090|Mus musculus", "brain", "positive")
    add("MSV1", "mzspec:MSV1:path/b.mzML", "10090|Mus musculus", "brain", "negative")
    add("MSV2", "mzspec:MSV2:x/c.mzML", "10090|Mus musculus", "brain", "missing value")
    add("MSV2", "mzspec:MSV2:y/c.mzML", "10090|Mus musculus", "brain", "positive")
    add("MSV3", "mzspec:MSV3:x/d.mzML", "9606|Homo sapiens", "liver", "alternating", "human")
    add("MSV4", "mzspec:MSV4/invalid", "9606|Homo sapiens", "liver", "positive", "human")
    add("MSV5", "mzspec:MSV5:path/E.mzML", "9606|Homo sapiens", "blood", "negative", "human")
    add("MSV5", "mzspec:MSV5:path/e.mzML", "9606|Homo sapiens", "blood", "negative", "human")
    rows.append(rows[0][:])  # Duplicate USI must not contribute twice.
    with gzip.open(path, "wt", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(g.REDU_COLS)
        w.writerows(rows)
    return rows


def test_spec_has_three_live_submodules():
    response = client.get("/api/spec")
    assert response.status_code == 200
    spec = response.json()
    assert [s for s in spec["subgroups"] if s["id"] == "gnps"] == [{"id": "gnps", "title": "Dataset Search"}]
    modules = [m for m in spec["modules"] if m["subgroup"] == "gnps"]
    assert [(m["id"], m["title"]) for m in modules] == [
        ("gnps_fastmasst", "Spectrum in which set"), ("gnps_redu", "Dataset finder"), ("gnps_massql", "Fragment pattern")]
    assert all(m["group"] == "modules" and m["live"] is True and m["example"] for m in modules)
    icons = (Path(__file__).parents[1] / "web" / "src" / "assets" / "icons.svg").read_text()
    assert 'symbol id="gnps"' in icons
    example = next(m for m in spec["modules"] if m["id"] == "gnps_fastmasst")["example"]
    assert example["name"] == "L-leucine" and example["precursor_mz"] == "132.102"
    blocks = {m["id"]: {b["id"]: b for b in m["outputs"] if b.get("id")} for m in spec["modules"] if m["subgroup"] == "gnps"}
    for module, block in (("gnps_fastmasst", "hits"), ("gnps_fastmasst", "redu"), ("gnps_redu", "datasets"), ("gnps_massql", "matches")):
        assert blocks[module][block]["page_size"] == 10 and blocks[module][block]["compact"] is True


def test_redu_streamed_search_matches_original(small_redu):
    filters = [("NCBITaxonomy", ["10090"]), ("UBERONBodyPartName", ["brain"])]
    expected = g.find_datasets(g.REDU_COLS, small_redu, filters, 1, "both")
    result, total, _ = redu_web.find(filters, 1, "both")
    assert result == expected and total == len(small_redu)
    response = _post("gnps_redu", organism="10090", body_part="brain", polarity="both")
    assert response.status_code == 200
    output = response.json()["outputs"]
    assert output["summary"]["values"]["Passed datasets"] == "1"
    assert output["summary"]["values"]["Undetermined"] == "1"
    assert len(output["downloads"]["files"]) == 2


def test_redu_example_runs_offline(small_redu):
    module = next(m for m in MODULES if m["id"] == "gnps_redu")
    response = _post(module["id"], **module["example"])
    assert response.status_code == 200, response.text
    outputs = response.json()["outputs"]
    assert set(outputs) == {block["id"] for block in module["outputs"]}
    assert outputs["summary"]["values"]["Matched files"] == "4"
    assert outputs["summary"]["values"]["Passed datasets"] == "0"


@pytest.mark.parametrize("background", g.REDU_BG_MODES)
def test_redu_targeted_join_and_enrichment_match_original(small_redu, background):
    usis = ["mzspec:MSV1:path/a.mzML:scan:1", "mzspec:MSV1:PATH/B.raw:scan:2",
            "mzspec:MSV2:other/c.mzXML:scan:3", "mzspec:MSV3:x/d.mzML:scan:4",
            "mzspec:MSV5:path/e.mzML:scan:5", "mzspec:MSV9:missing.mzML:scan:1"]
    rows = [{"USI": u, "Dataset": g.parse_usi(u)[0], "file_usi": g.file_usi_of(u),
             "library": g.REPO_LIB, "Cosine": 0.8, "Matching Peaks": 4} for u in usis]
    original = g.Redu(background)
    reduced = redu_web.TargetedRedu(background, {"masst_repository": rows})
    for u in usis:
        oi, ol, on = original.match(u)
        ri, rl, rn = reduced.match(u)
        assert (ol, on) == (rl, rn)
        assert (original.rows[oi] if oi is not None else None) == (reduced.rows[ri] if ri is not None else None)
    left, right = copy.deepcopy(rows), copy.deepcopy(rows)
    g.redu_annotate(original, left)
    g.redu_annotate(reduced, right)
    assert g.redu_summarise(original, left) == g.redu_summarise(reduced, right)
    assert g.redu_by_dataset(original, left) == g.redu_by_dataset(reduced, right)


def test_redu_input_and_lock_errors(small_redu):
    assert _post("gnps_redu").status_code == 400
    assert _post("gnps_redu", organism="mouse", min_files="2.5").status_code == 400
    assert _post("gnps_redu", where="not_a_column=x").status_code == 400
    with web._memory_job():
        response = _post("gnps_redu", organism="10090")
    assert response.status_code == 400
    assert "Another GNPS2" in response.json()["detail"]


def test_massql_mgf_and_download(tmp_path, monkeypatch):
    monkeypatch.setenv("MZDESK_CACHE_DIR", str(tmp_path))
    r = _post("gnps_massql", query=QUERY, input_file=_file("../../input.mgf", MGF))
    assert r.status_code == 200, r.text
    out = r.json()["outputs"]
    assert out["summary"]["values"] == {"Input file": "input.mgf", "MS2 spectra": "2", "Matches": "1"}
    assert out["matches"]["rows"][0][0] == "42"
    downloads = {f["name"]: f["content"] for f in out["downloads"]["files"]}
    assert "first" in downloads["massql-matches.mgf"] and "second" not in downloads["massql-matches.mgf"]
    assert "42\t300" in downloads["massql-matches.tsv"]


def test_massql_example_runs_offline(tmp_path, monkeypatch):
    monkeypatch.setenv("MZDESK_CACHE_DIR", str(tmp_path))
    module = next(m for m in MODULES if m["id"] == "gnps_massql")
    response = _post(module["id"], **module["example"])
    assert response.status_code == 200, response.text
    outputs = response.json()["outputs"]
    assert set(outputs) == {block["id"] for block in module["outputs"]}
    assert outputs["summary"]["values"] == {"Input file": "leucine-example.mgf", "MS2 spectra": "1", "Matches": "1"}
    assert outputs["matches"]["rows"][0][0] == "42"


def test_massql_large_result_count_and_full_download(tmp_path, monkeypatch):
    monkeypatch.setenv("MZDESK_CACHE_DIR", str(tmp_path))
    raw = "".join(f"BEGIN IONS\nTITLE=test {i}\nSCANS={i}\nPEPMASS=300\n100 100\nEND IONS\n" for i in range(205)).encode()
    r = _post("gnps_massql", query=QUERY, input_file=_file("input.mgf", raw))
    assert r.status_code == 200, r.text
    body = r.json()
    matches = body["outputs"]["matches"]
    assert matches["count"] == 205 and len(matches["rows"]) == 200
    assert "200 of 205" in " ".join(body["warnings"])
    files = {f["name"]: f["content"] for f in body["outputs"]["downloads"]["files"]}
    assert len(list(csv.DictReader(io.StringIO(files["massql-matches.tsv"]), delimiter="\t"))) == 205


def test_redu_large_result_count_and_full_download(tmp_path, monkeypatch):
    monkeypatch.setenv("MZDESK_CACHE_DIR", str(tmp_path))
    datasets = [{"dataset": f"MSV{i}", "source": "MassIVE", "n_files": 1,
                 "pos": 1, "neg": 0, "alt": 0, "unk": 0, "other": 0} for i in range(205)]
    result = {"n_files_matched": 205, "passed": datasets, "undetermined": [], "matched_values": {}}
    monkeypatch.setattr(redu_web, "find", lambda *args: (result, 205, "2026-10-12"))
    r = _post("gnps_redu", organism="mouse")
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["outputs"]["datasets"]["count"] == 205 and len(body["outputs"]["datasets"]["rows"]) == 200
    files = {f["name"]: f["content"] for f in body["outputs"]["downloads"]["files"]}
    assert len(list(csv.DictReader(io.StringIO(files["pan-redu-datasets.tsv"]), delimiter="\t"))) == 205


@pytest.mark.parametrize("extension", ["mzML", "mzXML"])
def test_massql_real_mzml_and_mzxml(tmp_path, monkeypatch, extension):
    om = pytest.importorskip("pyopenms")
    np = pytest.importorskip("numpy")
    monkeypatch.setenv("MZDESK_CACHE_DIR", str(tmp_path))
    spec = om.MSSpectrum()
    spec.setMSLevel(2)
    spec.setNativeID("controllerType=0 controllerNumber=1 scan=42")
    spec.setRT(120)
    spec.set_peaks((np.array([100.0, 110.0]), np.array([100.0, 20.0])))
    precursor = om.Precursor()
    precursor.setMZ(300)
    precursor.setCharge(1)
    spec.setPrecursors([precursor])
    experiment = om.MSExperiment()
    experiment.addSpectrum(spec)
    file = tmp_path / f"input.{extension}"
    (om.MzMLFile() if extension == "mzML" else om.MzXMLFile()).store(str(file), experiment)
    r = _post("gnps_massql", query=QUERY, input_file=_file(file.name, file.read_bytes()))
    assert r.status_code == 200, r.text
    output = r.json()["outputs"]
    assert output["summary"]["values"]["MS2 spectra"] == "1"
    assert output["summary"]["values"]["Matches"] == "1"
    if extension == "mzML":
        assert output["matches"]["rows"][0][0] == "42"


def test_massql_rejects_unsupported_and_invalid_input(tmp_path, monkeypatch):
    monkeypatch.setenv("MZDESK_CACHE_DIR", str(tmp_path))
    assert _post("gnps_massql", query="QUERY scanmaxint(MS2DATA) WHERE MS2PROD=100", input_file=_file("x.mgf", MGF)).status_code == 400
    assert _post("gnps_massql", query="QUERY MS2DATA WHERE MS2PREC=300:INTENSITYPERCENT=5", input_file=_file("x.mgf", MGF)).status_code == 400
    assert _post("gnps_massql", query=QUERY, input_file={"name": "x.mgf", "base64": "not base64!"}).status_code == 400
    assert _post("gnps_massql", query=QUERY, input_file=_file("x.xml", MGF)).status_code == 400
    assert _post("gnps_massql", query=QUERY, input_file=_file("x.mzML", b'<!DOCTYPE x [<!ENTITY x "bad">]>')).status_code == 400
    monkeypatch.setattr(web, "MAX_UPLOAD_BYTES", 5)
    assert _post("gnps_massql", query=QUERY, input_file=_file("x.mgf", MGF)).status_code == 400
    assert client.post("/api/modules/gnps_massql/run", content=b"{}", headers={"content-length": "27000001"}).status_code == 413


def test_massql_missing_metadata_is_not_exclusion_evidence():
    q = g.parse_massql("QUERY MS2DATA WHERE MS2PREC=300:EXCLUDED")
    assert not g._massql_eval_query(q, {"precursor_mz": None, "peaks": [(100, 1)]})[0]
    assert g._massql_mgf_spec({"meta": {"CHARGE": "-2-"}, "peaks": [(100, 1)]}, "x", 1)["charge"] == -2


def test_massql_worker_limits_output_file_size(tmp_path):
    src = tmp_path / "input.mgf"
    src.write_text("BEGIN IONS\nSCANS=1\nPEPMASS=300\n" + "".join(f"{100+i/10000:.4f} 100\n" for i in range(1000)) + "END IONS\n")
    query = tmp_path / "query.txt"
    query.write_text(QUERY)
    result = tmp_path / "result.json"
    code = "from mzdesk.gnps import massql_worker as w; import sys; w.MAX_OUTPUT=3000; w.main(*sys.argv[1:])"
    p = subprocess.run([sys.executable, "-c", code, str(src), str(query), str(result)], capture_output=True, text=True, timeout=30)
    assert p.returncode == 0, p.stderr
    out = json.loads(result.read_text())
    assert out["n_matches"] == 1 and "massql-matches.tsv" in out["files"]
    assert "massql-matches.mgf" not in out["files"] and "exceeds" in out["warnings"][0]


def test_request_limit_also_applies_without_content_length(monkeypatch):
    import asyncio
    from mzdesk.api import request_limits

    monkeypatch.setattr(request_limits, "MASSQL_REQUEST_LIMIT", 10)
    messages = iter([{"type": "http.request", "body": b"123456", "more_body": True},
                     {"type": "http.request", "body": b"123456", "more_body": False}])
    sent = []

    async def receive():
        return next(messages)

    async def send(message):
        sent.append(message)

    async def should_not_run(*args):
        raise AssertionError("oversized request reached the app")

    scope = {"type": "http", "method": "POST", "path": "/api/modules/gnps_massql/run", "headers": []}
    asyncio.run(request_limits.MassQLRequestLimit(should_not_run)(scope, receive, send))
    assert sent[0]["status"] == 413


def test_fastmasst_fake_search_diagnostics_and_redu(small_redu, monkeypatch):
    repository = [{"USI": "mzspec:MSV1:path/a.mzML:scan:1", "Dataset": "MSV1", "Cosine": 0.9,
                   "Matching Peaks": 4, "Delta Mass": 0.0},
                  {"USI": "mzspec:MSV1:path/b.mzML:scan:2", "Dataset": "MSV1", "Cosine": 0.8,
                   "Matching Peaks": 2, "Delta Mass": 0.0}]

    def search(peaks, precursor, library, cosine, min_peaks, pm_tol, frag_tol, **kw):
        raw = repository if library == g.REPO_LIB else []
        return g.summarise_raw(raw, library, min_peaks, "task1"), {"task_id": "task1", "raw": raw, "elapsed_s": 1, "n_poll_errors": 0}

    monkeypatch.setattr(g, "fasst_search", search)
    monkeypatch.setattr(g.UsiFetcher, "fetch", lambda self, usi: ({"peaks": [[100.0, 100], [110, 20]], "precursor_mz": 300}, "ok"))
    r = _post("gnps_fastmasst", spectrum={"source": "Paste", "value": "100 100\n110 20"}, precursor_mz="300",
              library="repository", diagnostic="100", redu_join=True)
    assert r.status_code == 200, r.text
    out = r.json()["outputs"]
    assert out["searches"]["rows"][0][1] == "1"
    assert len(out["hits"]["rows"]) == 1
    assert out["verification"]["rows"][0][2:4] == ["1", "1"]
    assert out["redu"]["rows"]
    files = {f["name"]: f["content"] for f in out["downloads"]["files"]}
    assert "fastmasst.hits.tsv" in files and len(list(csv.DictReader(io.StringIO(files["fastmasst.hits.tsv"]), delimiter="\t"))) == 2
    assert "usi_cache" not in json.loads(files["fastmasst.json"])["verification"]


def test_fastmasst_failures_are_not_zero_hits(tmp_path, monkeypatch):
    monkeypatch.setenv("MZDESK_CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(g, "fasst_search", lambda *a, **k: ({"error": "offline"}, {"task_id": None, "raw": None, "elapsed_s": 0, "n_poll_errors": 0}))
    r = _post("gnps_fastmasst", spectrum={"source": "Paste", "value": "100 1"}, precursor_mz="300", library="repository")
    assert r.status_code == 502 and "offline" in r.json()["detail"]


def test_reference_and_suspect_library_hits_are_separate(tmp_path, monkeypatch):
    monkeypatch.setenv("MZDESK_CACHE_DIR", str(tmp_path))
    raw = [{"USI": "mzspec:GNPS:GNPS-LIBRARY:accession:A", "Dataset": "GNPS Library", "Cosine": 0.9,
            "Matching Peaks": 4, "CompoundName": "Known compound", "GNPSLibraryAccession": "A"},
           {"USI": "mzspec:GNPS:GNPS-LIBRARY:accession:B", "Dataset": "GNPS Library", "Cosine": 0.8,
            "Matching Peaks": 4, "CompoundName": "Suspect related to X", "GNPSLibraryAccession": "B"}]
    monkeypatch.setattr(g, "fasst_search", lambda p, m, lib, co, n, pt, ft, **kw:
                        (g.summarise_raw(raw, lib, n, "test"), {"task_id": "test", "raw": raw, "elapsed_s": 0, "n_poll_errors": 0}))
    r = _post("gnps_fastmasst", spectrum={"source": "Paste", "value": "100 1"}, precursor_mz="300", library="library")
    assert r.status_code == 200, r.text
    out = r.json()["outputs"]
    assert out["reference_hits"]["rows"][0][0] == "Known compound"
    assert out["suspect_hits"]["rows"][0][0] == "Suspect related to X"


def test_fastmasst_large_preview_and_missing_download_warning(tmp_path, monkeypatch):
    monkeypatch.setenv("MZDESK_CACHE_DIR", str(tmp_path))
    raw = [{"USI": f"mzspec:MSV{i}:path/a.mzML:scan:1", "Dataset": f"MSV{i}", "Cosine": 0.9, "Matching Peaks": 4} for i in range(205)]
    monkeypatch.setattr(g, "fasst_search", lambda p, m, lib, co, n, pt, ft, **kw:
                        (g.summarise_raw(raw, lib, n, "test"), {"task_id": "test", "raw": raw, "elapsed_s": 0, "n_poll_errors": 0}))
    monkeypatch.setattr(web, "MAX_DOWNLOAD_BYTES", 50)
    r = _post("gnps_fastmasst", spectrum={"source": "Paste", "value": "100 1"}, precursor_mz="300", library="repository")
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["outputs"]["hits"]["count"] == 205 and len(body["outputs"]["hits"]["rows"]) == 200
    assert any("not included in the web download" in w for w in body["warnings"])
    assert not any("use the TSV download" in w for w in body["warnings"])


def test_local_cli_batch_dry_run_and_massql(tmp_path):
    queries = tmp_path / "queries.json"
    queries.write_text(json.dumps([{"name": "example", "precursor": 300, "peaks": [[100, 100], [110, 20]]}]))
    batch = subprocess.run([sys.executable, "-m", "mzdesk.gnps.compound_db_check", "--queries", str(queries),
                            "--out", str(tmp_path / "batch"), "--dry-run"], capture_output=True, text=True, encoding="utf-8", timeout=20)
    assert batch.returncode == 0, batch.stderr
    assert "dry-run" in batch.stdout
    mgf = tmp_path / "input.mgf"
    mgf.write_bytes(MGF)
    output = tmp_path / "result.json"
    result = subprocess.run([sys.executable, "-m", "mzdesk.gnps.compound_db_check", "--massql", QUERY,
                             "--massql-input", str(mgf), "--out", str(output)], capture_output=True, text=True, encoding="utf-8", timeout=20)
    assert result.returncode == 0, result.stderr
    assert json.loads(output.read_text())["n_matches"] == 1


def test_local_cli_replay_does_not_call_fastmasst(tmp_path, monkeypatch):
    raw = [{"USI": "mzspec:MSV1:path/a.mzML:scan:1", "Dataset": "MSV1", "Cosine": 0.9, "Matching Peaks": 4},
           {"USI": "mzspec:MSV1:path/b.mzML:scan:2", "Dataset": "MSV1", "Cosine": 0.8, "Matching Peaks": 2}]
    monkeypatch.setattr(g, "fasst_search", lambda p, m, lib, co, n, pt, ft, **kw:
                        (g.summarise_raw(raw, lib, n, "test"), {"task_id": "test", "raw": raw, "elapsed_s": 0, "n_poll_errors": 0}))
    options = {**g.OPT_DEFAULTS, "name": "test", "precursor": 300, "peaks": [[100, 100]],
               "library": "repository", "keep_all": True}
    first = g.run_query(options, str(tmp_path / "first.json"), {})
    assert first["masst_repository"]["n_matches"] == 1

    def no_network(*args, **kwargs):
        raise AssertionError("replay must not make a fastMASST request")

    monkeypatch.setattr(g, "fasst_search", no_network)
    replay = g.run_query({**options, "from_hits": str(tmp_path / "first.hits.tsv"), "min_peaks": 4},
                         str(tmp_path / "replayed.json"), {})
    assert replay["masst_repository"]["n_matches"] == 1


def test_fastmasst_example_runs_offline_with_fake_service(tmp_path, monkeypatch):
    monkeypatch.setenv("MZDESK_CACHE_DIR", str(tmp_path))

    def search(peaks, precursor, library, cosine, min_peaks, pm_tol, frag_tol, **kw):
        assert precursor == 132.102 and len(peaks) == 31
        return g.summarise_raw([], library, min_peaks, "test"), {"task_id": "test", "raw": [], "elapsed_s": 0, "n_poll_errors": 0}

    monkeypatch.setattr(g, "fasst_search", search)
    mod = next(m for m in MODULES if m["id"] == "gnps_fastmasst")
    r = _post(mod["id"], **mod["example"])
    assert r.status_code == 200, r.text
    assert {b["id"] for b in mod["outputs"]} <= set(r.json()["outputs"])
