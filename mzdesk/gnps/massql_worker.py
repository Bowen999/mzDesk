"""Isolated, time- and memory-bounded MassQL web execution (not the local CLI)."""
import csv
import json
import os
import resource
import sys
from pathlib import Path

from mzdesk.gnps import compound_db_check as g

MAX_OUTPUT = 10_000_000
PREVIEW_ROWS = 200


def _limits():
    # Native mzML readers may expand compressed peak arrays far beyond the uploaded file size.
    # On Linux cap additional virtual memory after importing the libraries; file size is capped on Unix.
    if sys.platform.startswith("linux"):
        with open("/proc/self/statm", encoding="ascii") as f:
            current = int(f.read().split()[0]) * os.sysconf("SC_PAGE_SIZE")
        cap = current + 768 * 1024 * 1024
        _, hard = resource.getrlimit(resource.RLIMIT_AS)
        if hard != resource.RLIM_INFINITY:
            cap = min(cap, hard)
        resource.setrlimit(resource.RLIMIT_AS, (cap, hard))
    _, hard = resource.getrlimit(resource.RLIMIT_FSIZE)
    cap = min(MAX_OUTPUT, hard) if hard != resource.RLIM_INFINITY else MAX_OUTPUT
    resource.setrlimit(resource.RLIMIT_FSIZE, (cap, hard))


def main(path, query_path, result_path):
    if not path.lower().endswith(".mgf"):
        import pyopenms  # noqa: F401 -- import before measuring the baseline virtual memory
    _limits()
    with open(query_path, encoding="utf-8") as f:
        parsed = g.parse_massql(f.read())
    spectra = g._massql_load_mgf(path) if path.lower().endswith(".mgf") else g._massql_load_mzml(path)
    matches = []
    for spectrum in spectra:
        ok, matched = g._massql_eval_query(parsed, spectrum)
        if ok:
            matches.append({**spectrum, "matched": "; ".join(matched)})
    cols = ["scan", "precursor_mz", "rt_min", "polarity", "charge", "n_peaks", "title", "matched"]

    def values(m):
        return [m.get("scan"), m.get("precursor_mz"), m.get("rt_min"), m.get("polarity"), m.get("charge"),
                len(m["peaks"]), m.get("title"), m.get("matched")]

    # Do not leak the temporary server path when a spectrum has no title.
    for m in matches:
        m["title"] = (m.get("title") or "").replace(path, "input")
    directory = Path(result_path).parent
    files, warnings = [], []
    tsv = directory / "massql-matches.tsv"
    try:
        with tsv.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(cols)
            w.writerows(values(m) for m in matches)
        files.append(tsv.name)
    except OSError:
        tsv.unlink(missing_ok=True)
        warnings.append("The TSV exceeds 10 MB; use the local CLI to export all matches.")
    mgf = directory / "massql-matches.mgf"
    try:
        g._massql_write_mgf(str(mgf), matches)
        files.append(mgf.name)
    except OSError:
        mgf.unlink(missing_ok=True)
        warnings.append("The matched MGF exceeds 10 MB; use the local CLI to export it.")
    result = {"n_spectra": len(spectra), "n_matches": len(matches),
              "rows": [values(m) for m in matches[:PREVIEW_ROWS]], "files": files, "warnings": warnings}
    with open(result_path, "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, allow_nan=False)


if __name__ == "__main__":
    main(*sys.argv[1:])
