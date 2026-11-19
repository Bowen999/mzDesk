"""Web adapters for the migrated GNPS2 tools; local CLI remains available."""
import base64
import binascii
import contextlib
import csv
import io
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

try:
    import fcntl
except ImportError:                  # Windows: msvcrt locks the same file
    fcntl = None
    import msvcrt

from mzdesk.api import common
from mzdesk.api.common import EngineError, InputError, fmt
from mzdesk.gnps import compound_db_check as g, hit_context as hc, redu_web

MAX_UPLOAD_BYTES = 20_000_000
MAX_DOWNLOAD_BYTES = 10_000_000
PREVIEW_ROWS = 200


@contextlib.contextmanager
def _memory_job():
    """One memory-intensive GNPS operation across all API workers."""
    base = hc.cache_dir()
    base.mkdir(parents=True, exist_ok=True)
    fd = os.open(base / "web.lock", os.O_CREAT | os.O_RDWR | getattr(os, "O_NOFOLLOW", 0), 0o600)
    try:
        try:
            if fcntl:
                fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            else:
                msvcrt.locking(fd, msvcrt.LK_NBLCK, 1)
        except (BlockingIOError, PermissionError):      # the second is what msvcrt raises when the lock is taken
            raise InputError("Another GNPS2 data operation is running. Please try again shortly.") from None
        yield
    finally:
        os.close(fd)


def _integer(inputs, key, label, default, minimum, maximum):
    n = common.number(inputs, key, label, default=default, minimum=minimum, maximum=maximum)
    if not float(n).is_integer():
        raise InputError(f"{label} must be a whole number.")
    return int(n)


def _choice(inputs, key, default, choices, label):
    v = inputs.get(key) or default
    if not isinstance(v, str) or v not in choices:
        raise InputError(f"{label} must be one of: {', '.join(choices)}.")
    return v


def _csv_text(columns, rows):
    out = io.StringIO(newline="")
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(columns)
    writer.writerows(rows)
    return out.getvalue()


def _download(name, text):
    return {"name": name, "content": text, "mime": "text/plain;charset=utf-8"}


def _text(value):
    """A cell: a real 0 stays 0, only a missing value (None, empty) is a dash."""
    return "—" if value is None or value == "" else str(value)


REDU_PROPERTY = {"SampleType": "Sample type", "SampleTypeSub1": "Sample subtype", "NCBITaxonomy": "Organism",
                 "UBERONBodyPartName": "Body part", "DOIDCommonName": "Disease", "HealthStatus": "Health status"}
# (TSV column, header, kind); "Search" and "Hit set" are left out of a run in which they have one value
REDU_COLUMNS = [("search", "Search", "search"), ("hit_set", "Hit set", "hitset"), ("column", "Property", "property"), ("value", "Value", "value"),
                ("n_hit_files", "Hit files", "text"), ("n_hit_datasets", "Hit datasets", "text"), ("bg_datasets", "Background datasets", "text"),
                ("file_enrich", "File enrichment", "enrich"), ("ds_enrich", "Dataset enrichment", "enrich"),
                ("dsw_enrich", "Dataset-weighted enrichment", "enrich")]
REDU_CONSTANT = ("search", "hit_set")
# the words of the "Diagnostic verification" table's columns
REDU_HIT_SET = {"all": "All", "verified": "Verified", "diag_pass": "Diagnostic pass"}


def _enrichment(value):
    """2 decimals; 1000 and more is ">1000"; none is a dash."""
    try:
        x = float(value)
    except (TypeError, ValueError):
        return "—"
    if not math.isfinite(x):
        return "—"
    return ">1000" if round(x, 2) >= 1000 else f"{x:.2f}"


def _redu_cell(kind, value):
    if kind == "enrich":
        return _enrichment(value)
    if kind == "search":
        return _text((value or "").replace("masst_", ""))
    if kind == "hitset":
        return _text(REDU_HIT_SET.get(value, value))
    return _text(value)


def redu_preview(path):
    """The "Pan-ReDU enrichment" block from the redu_summary.tsv (already in the page's order): PREVIEW_ROWS rows, the full count, and
    only the columns that vary. The TSV itself keeps every column at full precision."""
    if not os.path.isfile(path):
        return {"rows": [], "count": 0}
    shown, total, seen = [], 0, {k: set() for k in REDU_CONSTANT}
    with open(path, encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            total += 1
            for k in REDU_CONSTANT:
                seen[k].add(row.get(k))
            if len(shown) < PREVIEW_ROWS:
                shown.append(row)
    if not shown:
        return {"rows": [], "count": 0}
    cols = [c for c in REDU_COLUMNS if not (c[0] in REDU_CONSTANT and len(seen[c[0]]) == 1)]
    rows = []
    for row in shown:
        cells = []
        for key, _, kind in cols:
            v = row.get(key)
            if kind == "property":
                v = REDU_PROPERTY.get(v, v)
            elif kind == "value":
                v = (hc.organism(v) or v) if row.get("column") == "NCBITaxonomy" else v
            cells.append(_redu_cell(kind, v))
        rows.append(cells)
    return {"rows": rows, "count": total, "columns": [c[1] for c in cols]}


def _file_downloads(result, warnings):
    files = []
    for name, path in result.get("outputs", {}).items():
        if not os.path.isfile(path):
            continue
        if os.path.getsize(path) > MAX_DOWNLOAD_BYTES:
            warnings.append(f"{name} exceeds 10 MB and is not included in the web download; use the local CLI for complete files.")
            continue
        with open(path, encoding="utf-8") as f:
            files.append(_download(Path(path).name, f.read()))
    return files


def _ions(inputs, key, label):
    try:
        values = g.parse_mz_list(common.text(inputs, key))
    except (TypeError, ValueError):
        raise InputError(f"{label}: use comma-separated positive m/z values.") from None
    if len(values) > 20 or any(not math.isfinite(v) or v <= 0 for v in values):
        raise InputError(f"{label}: use at most 20 positive, finite m/z values.")
    return values


def _cosine(row):
    try:
        return fmt(float(row.get("Cosine")), 3)
    except (TypeError, ValueError):
        return "—"


def _hit_cells(r):
    """One hit of "Passing hits by dataset"; a value the source does not have is a dash."""
    return [_text(r.get("USI")), _cosine(r), _text(r.get("Matching Peaks")), _text(r.get("CompoundName")), _text(r.get("Adduct")),
            hc.redu_value(r, "redu_SampleType") or "—", hc.redu_value(r, "redu_NCBITaxonomy") or "—",
            hc.redu_value(r, "redu_UBERONBodyPartName") or "—", hc.domain_text(hc.domain_of(r)) or "—",
            _text(r.get("verify_status")), _text(r.get("diag_pass"))]


def _add_context_columns(path, info):
    """The hits table of the download gets dataset_title, dataset_species and domain_masst (every row, not only the passing ones)."""
    if not path or not os.path.isfile(path):
        return
    tmp = path + ".tmp"
    with open(path, encoding="utf-8", newline="") as f, open(tmp, "w", encoding="utf-8", newline="") as out:
        reader = csv.DictReader(f, delimiter="\t")
        extra = ["dataset_title", "dataset_species", "domain_masst"]
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow([*reader.fieldnames, *extra])
        for r in reader:
            known = info.get((r.get("Dataset") or "").upper()) or {}
            writer.writerow([*(r.get(c, "") for c in reader.fieldnames), known.get("title") or "", "; ".join(hc.species_list(known.get("species"))),
                             hc.domain_text(hc.domain_of({"USI": r.get("USI")}))])
    os.replace(tmp, path)


def _context_downloads(datasets, info, distribution, joined):
    """Every dataset (not only the 200 on the page) and the full distribution as TSV files."""
    rows = []
    for d in datasets:
        known = info.get(d["dataset"].upper()) or {}
        if d["library"]:
            status = ""
        elif not joined:
            status = "not joined"
        else:
            status = {"matched": f"matched {d['files_matched']} of {d['files']} files", "dataset only": "dataset in Pan-ReDU, no file matched",
                      "not in Pan-ReDU": "not in Pan-ReDU"}.get(d["redu"], "unknown")
        rows.append([d["search"].replace("masst_", ""), "GNPS reference library" if d["library"] else d["dataset"] or hc.NO_DATASET, known.get("title") or "",
                     "; ".join(hc.species_list(known.get("species"))), known.get("repository") or "", d["hits"], d["files"], f"{d['best']:.3f}",
                     d["sample_type"] or "", d["organism"] or "", d["body_part"] or "", status])
    cols = ["search", "dataset", "title", "species", "repository", "hits", "files", "best_cosine", "top_sample_type", "top_organism",
            "top_body_part", "redu_status"]
    return [_download("fastmasst.datasets.tsv", _csv_text(cols, rows)),
            _download("fastmasst.distribution.tsv", _csv_text(["field", "value", "hit_files", "hit_datasets", "share_of_hit_files"], distribution))]


def run_fastmasst(inputs):
    warnings = []
    peaks, header_mz = common.spectrum(inputs, "spectrum", "MS2 spectrum", warnings)
    mz = common.number(inputs, "precursor_mz", "Precursor m/z", positive=True) or header_mz
    if mz is None or not math.isfinite(mz) or mz <= 0:
        raise InputError("Enter a positive precursor m/z or include it in the spectrum.")
    library = _choice(inputs, "library", "both", ("both", "repository", "library"), "Search")
    diag, required = _ions(inputs, "diagnostic", "Diagnostic ions"), _ions(inputs, "diag_required", "Required ions")
    merged = diag + [x for x in required if all(abs(x - d) >= 1e-6 for d in diag)]
    if len(merged) > 20 or len({g.mz_label(x) for x in merged}) != len(merged):
        raise InputError("Use at most 20 distinct diagnostic ions; values must differ at four decimal places.")
    rule = common.text(inputs, "diag_rule") or "all"
    if merged:
        try:
            g.parse_diag_rule(rule, len(merged))
        except ValueError as exc:
            raise InputError(str(exc)) from None
    elif rule != "all":
        raise InputError("A diagnostic rule requires diagnostic ions.")
    join = common.flag(inputs, "redu_join", "Join Pan-ReDU")
    if join and not g.redu_installed():
        # The join is on by default, so a computer that never ran `mzdesk setup redu` must still get its results: one warning, no join.
        warnings.append("The Pan-ReDU metadata is not on this computer, so sample types, organisms and body parts are not shown. "
                        "To add them, run `mzdesk setup redu` (a 1.7 GB download).")
        join = False
    bg = _choice(inputs, "redu_bg", "pos_lc_or_unknown", tuple(g.REDU_BG_MODES), "ReDU background")
    o = {**g.OPT_DEFAULTS, "name": common.text(inputs, "name") or "query", "peaks": peaks, "precursor": mz,
         "library": library, "cosine": common.number(inputs, "cosine", "Cosine threshold", default=0.7, minimum=0, maximum=1),
         "min_peaks": _integer(inputs, "min_peaks", "Minimum matched peaks", 3, 1, 1000),
         "pm_tol": common.number(inputs, "pm_tol", "Precursor tolerance", default=0.05, positive=True, maximum=5),
         "frag_tol": common.number(inputs, "frag_tol", "Fragment tolerance", default=0.05, positive=True, maximum=5),
         "max_wait": 300, "keep_all": True, "diagnostic": diag, "diag_required": required, "diag_rule": rule,
         "diag_tol": common.number(inputs, "diag_tol", "Diagnostic tolerance", default=0.01, positive=True, maximum=5),
         "diag_min_rel": common.number(inputs, "diag_min_rel", "Diagnostic intensity", default=0, minimum=0, maximum=100),
         "verify_max": _integer(inputs, "verify_max", "Maximum spectra to verify per search", 20, 0, 20),
         "query_min_rel": common.number(inputs, "query_min_rel", "Query noise threshold", default=0, minimum=0, maximum=100),
         "redu_join": join, "redu_bg": bg}
    with tempfile.TemporaryDirectory(prefix="mzdesk-gnps-") as directory:
        ctx = {"usi_cache": os.path.join(directory, "usi"), "redu_factory": redu_web.TargetedRedu}
        # Limit web retries and per-request timeout; the CLI retains the original settings.
        ctx["fetcher"] = g.UsiFetcher(ctx["usi_cache"], rate=2.5, retries=1, timeout=12)
        try:
            with _memory_job() if join else contextlib.nullcontext():      # only the Pan-ReDU join is memory-heavy
                result = g.run_query(o, os.path.join(directory, "fastmasst.json"), ctx)
        except (OSError, ValueError, RuntimeError, KeyError) as exc:
            raise EngineError(f"GNPS2 search failed: {exc}") from None
        searches = [(k, result[k]) for k, _ in g.resolve_libraries(library)]
        failures = [(k, r) for k, r in searches if "error" in r]
        if len(failures) == len(searches):
            raise EngineError("fastMASST did not complete: " + "; ".join(str(r["error"]) for _, r in failures))
        warnings.extend(f"{key}: search failed ({r['error']}); results are incomplete." for key, r in failures)
        overview = [[key.replace("masst_", ""), _text(r.get("n_matches")), _text(r.get("n_before_filter")),
                     fmt(r.get("max_cosine"), 3), _text(r.get("task_id"))] for key, r in searches]
        library_result = result.get("masst_library") or {}

        def library_rows(key):
            return [[_text(h.get("name")), fmt(h.get("cosine"), 3), _text(h.get("matching_peaks")),
                     _text(h.get("adduct")), _text(h.get("accession")), str(h.get("n_entries") or 0)]
                    for h in library_result.get(key, [])]

        references, suspects = library_rows("library_hits"), library_rows("suspect_hits")
        hits_path = result.get("outputs", {}).get("hits_tsv", "")
        hits = hc.read_hits(hits_path)
        datasets = hc.dataset_rows(hits, join)
        info = hc.explorer_info([d["dataset"].upper() for d in datasets if not d["library"] and d["dataset"]], warnings)
        hits_block = hc.hits_block(datasets, info, join, len(hits), _hit_cells)
        groups, distribution, coverage = hc.distribution(hits, join, info)
        verification = result.get("verification") or {}
        ver = [[key.replace("masst_", ""), str(v.get("n_candidates", 0)), str(v.get("n_verified_ok", 0)),
                str(v.get("n_diag_pass", 0)), str(v.get("n_not_selected", 0)), str(v.get("n_fetch_failed", 0))]
               for key, v in verification.items() if key.startswith("masst_") and isinstance(v, dict)]
        redu_block = redu_preview(result.get("outputs", {}).get("redu_summary_tsv", ""))
        _add_context_columns(hits_path, info)
        downloads = _file_downloads(result, warnings)
        downloads += _context_downloads(datasets, info, distribution, join)
        public = dict(result)
        public["outputs"] = {k: Path(v).name for k, v in result.get("outputs", {}).items()}
        if public.get("verification"):
            public["verification"] = {k: v for k, v in public["verification"].items() if k != "usi_cache"}
        downloads.insert(0, _download("fastmasst.json", json.dumps(public, ensure_ascii=False, indent=2)))
    where = {"groups": groups, **({"note": coverage} if coverage else {})}
    return {"outputs": {"searches": {"rows": overview}, "distribution": where, "redu": redu_block, "hits": hits_block,
                        "reference_hits": {"rows": references}, "suspect_hits": {"rows": suspects},
                        "verification": {"rows": ver}, "downloads": {"files": downloads}}, "warnings": warnings}


def run_redu(inputs):
    organism, part = common.text(inputs, "organism"), common.text(inputs, "body_part")
    filters = [(k, v.split(",")) for k, v in (("NCBITaxonomy", organism), ("UBERONBodyPartName", part)) if v]
    extra = common.text(inputs, "where")
    for line in extra.splitlines():
        if not line.strip():
            continue
        col, sep, value = line.partition("=")
        col, value = col.strip(), value.strip()
        if not sep or col not in g.REDU_COLS or not value or not any(x.strip() for x in value.split(",")):
            raise InputError("Extra filters must be one column=value per line, using a Pan-ReDU column name.")
        filters.append((col, value.split(",")))
    if not filters or any(not any(x.strip() for x in terms) for _, terms in filters):
        raise InputError("Provide at least one non-empty organism, body part or extra filter.")
    minimum = _integer(inputs, "min_files", "Minimum files", 1, 1, 1_000_000)
    polarity = _choice(inputs, "polarity", "any", ("any", "positive", "negative", "both"), "Polarity")
    try:
        with _memory_job():
            result, total, date = redu_web.find(filters, minimum, polarity)
    except InputError:
        raise
    except (OSError, ValueError, KeyError) as exc:
        raise EngineError(f"Pan-ReDU cache could not be read: {exc}") from None
    cols = ["status", "dataset", "source", "n_files", "pos", "neg", "alt", "unk", "other"]
    datasets = [[status, *(d.get(k, "") for k in cols[1:])] for status, key in (("passed", "passed"), ("undetermined", "undetermined")) for d in result[key]]
    warnings = []
    if len(datasets) > PREVIEW_ROWS:
        warnings.append(f"Only the first {PREVIEW_ROWS} of {len(datasets)} datasets are shown.")
    matched = [[c, v, str(n)] for c, values in result["matched_values"].items() for v, n in values]
    shown = {"rows": datasets[:PREVIEW_ROWS], "count": len(datasets)}
    if result["undetermined"]:
        shown["help"] = "Undetermined: some matched files report no polarity."
    return {"outputs": {"summary": {"values": {"Cache date": date, "Cached rows": str(total), "Matched files": str(result["n_files_matched"]),
                                                   "Passed datasets": str(len(result["passed"])), "Undetermined": str(len(result["undetermined"]))}},
                        "datasets": shown, "values": {"rows": matched, "count": len(matched)},
                        "downloads": {"files": [_download("pan-redu-datasets.tsv", _csv_text(cols, datasets)),
                                               _download("pan-redu.json", json.dumps({"query": {"filters": filters, "min_files": minimum, "polarity": polarity},
                                                                                     "cache_date": date, **result}, ensure_ascii=False, indent=2))]}},
            "warnings": warnings}


def run_massql(inputs):
    query = common.text(inputs, "query")
    try:
        parsed = g.parse_massql(query)
    except (ValueError, TypeError) as exc:
        raise InputError(f"MassQL: {exc}") from None
    if parsed["function"] and parsed["function"].lower() != "scaninfo":
        raise InputError("The web page supports scaninfo(MS2DATA) or MS2DATA only; aggregation functions are not implemented.")
    data = inputs.get("input_file")
    if not isinstance(data, dict) or not isinstance(data.get("base64"), str) or not isinstance(data.get("name"), str):
        raise InputError("Upload one mzML, mzXML or MGF file.")
    name = Path(data["name"]).name
    extension = Path(name).suffix.lower()
    if extension not in (".mzml", ".mzxml", ".mgf"):
        raise InputError("Use an mzML, mzXML or MGF file.")
    payload = data["base64"]
    if len(payload) > 4 * ((MAX_UPLOAD_BYTES + 2) // 3):
        raise InputError("Upload one non-empty file of at most 20 MB.")
    try:
        encoded = base64.b64decode(payload, validate=True)
    except (binascii.Error, ValueError):
        raise InputError("The uploaded file is invalid; choose it again.") from None
    if not encoded or len(encoded) > MAX_UPLOAD_BYTES:
        raise InputError("Upload one non-empty file of at most 20 MB.")
    if extension != ".mgf" and (b"\0" in encoded or b"<!DOCTYPE" in encoded.upper() or b"<!ENTITY" in encoded.upper()):
        raise InputError("XML document type declarations and entities are not accepted.")
    with _memory_job(), tempfile.TemporaryDirectory(prefix="mzdesk-massql-") as directory:
        path = os.path.join(directory, "input" + extension)
        with open(path, "wb") as f:
            f.write(encoded)
        query_path = os.path.join(directory, "query.txt")
        result_path = os.path.join(directory, "result.json")
        with open(query_path, "w", encoding="utf-8") as f:
            f.write(query)
        env = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
        try:
            done = subprocess.run([sys.executable, "-m", "mzdesk.gnps.massql_worker", path, query_path, result_path],
                                  capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=180, env=env)
        except subprocess.TimeoutExpired:
            raise EngineError("MassQL exceeded the 180-second web limit. Use the local CLI for larger jobs.") from None
        if done.returncode or not os.path.isfile(result_path):
            if extension != ".mgf" and "No module named 'pyopenms'" in done.stderr:
                raise EngineError("Reading mzML/mzXML requires pyopenms. Install the updated mzDesk requirements.")
            raise EngineError("MassQL could not read or process the file; check its format and size. "
                              "The web worker also has a memory limit.", (done.stderr or done.stdout)[-2000:])
        with open(result_path, encoding="utf-8") as f:
            result = json.load(f)
        files = []
        for filename in result["files"]:
            with open(os.path.join(directory, filename), encoding="utf-8") as f:
                files.append(_download(filename, f.read()))
    warnings = [*result["warnings"]]
    if result["n_matches"] > PREVIEW_ROWS:
        warnings.append(f"Only the first {PREVIEW_ROWS} of {result['n_matches']} matches are shown.")
    display = [["—" if v is None or v == "" else str(v) for v in row] for row in result["rows"]]
    return {"outputs": {"summary": {"values": {"Input file": name, "MS2 spectra": str(result["n_spectra"]), "Matches": str(result["n_matches"])}},
                        "matches": {"rows": display, "count": result["n_matches"]}, "downloads": {"files": files}}, "warnings": warnings}


RUNNERS = {"gnps_fastmasst": run_fastmasst, "gnps_redu": run_redu, "gnps_massql": run_massql}
