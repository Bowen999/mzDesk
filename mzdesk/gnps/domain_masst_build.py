"""python -m mzdesk.gnps.domain_masst_build [--sha SHA] [--source DIR] [--out DIR]

Rebuilds mzdesk/gnps/data/domain_masst.tsv.gz (+ domain_masst.json) from the tables of github.com/robinschmid/microbe_masst (MIT).
The page "Spectrum in Which Set" joins its hit files to this table to say which microbe / plant / tissue / food / personal-care
sample a hit file belongs to. The table ships with mzDesk (no download at run time).

Kept per file: dataset, file_key (lower-case basename without extension: the join key; the MASST tables list a file as
mzspec:<dataset>:<basename>), domain, category. Blank and QC files are dropped (a spectrum in a blank is not an occurrence in a
microbe). The 59 MB microbiome table is not used. The source commit is pinned: a new snapshot is a deliberate change of SHA.

--source DIR  use already downloaded copies of the five tables (same file names) instead of downloading them.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import re
import sys
import tempfile
import time
import urllib.request

REPO = "robinschmid/microbe_masst"
SHA = "5dd2a0f6981697d9a10c26750c45013ca50ab328"        # master on 2026-08-14
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
EXT_RE = re.compile(r"\.(mzml\.gz|mzxml\.gz|mzml|mzxml|mzdata\.xml|mzdata|mzmlb|mz5|raw|cdf|wiff2?|d|mgf)$", re.I)
# domain -> (file in data/, category column, row filter)
TABLES = {
    "microbe": ("microbe_masst_table.csv", "Taxaname_file", lambda r: r.get("Blank") != "Yes" and r.get("QC") != "Yes"),
    "plant": ("plant_masst_table.csv", "Taxaname_file", lambda r: r.get("Blank") != "Yes" and r.get("QC") != "Yes"),
    "tissue": ("tissue_masst_table.csv", "ID", lambda r: not r.get("ID", "").startswith("blank_qc")),
    "food": ("food_masst_table.csv", "node_id", lambda r: True),          # node_id is already the food's name
    "personal_care": ("personalCareProduct_masst_table.csv", "Category", lambda r: True),
}
COLUMNS = ["dataset", "file_key", "domain", "category"]


def file_key(usi):
    """'mzspec:MSV000078547:name.mzXML' -> ('MSV000078547', 'name')"""
    head, _, rest = (usi or "").strip().partition(":")
    if head.lower() != "mzspec":
        return None
    dataset, _, run = rest.partition(":")
    run = EXT_RE.sub("", run.rsplit("/", 1)[-1]).strip().lower()
    return (dataset.strip().upper(), run) if dataset and run else None


def clean_category(domain, value):
    value = (value or "").strip()
    if domain == "personal_care" and value.startswith("Root > "):
        value = value[len("Root > "):]
    return value


def read_table(domain, path):
    _, column, keep = TABLES[domain]
    out, skipped = set(), 0
    with open(path, encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f):
            key = file_key(r.get("file_usi"))
            category = clean_category(domain, r.get(column))
            if not keep(r):
                skipped += 1
                continue
            if key is None or not category:
                continue
            out.add((key[0], key[1], domain, category))
    return out, skipped


def download(directory):
    for name, _, _ in TABLES.values():
        url = f"https://raw.githubusercontent.com/{REPO}/{SHA}/data/{name}"
        print(f"download {url}")
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "mzdesk-build"}), timeout=120) as r, \
                open(os.path.join(directory, name), "wb") as f:
            f.write(r.read())


def build(source, out_dir):
    rows, counts = set(), {}
    for domain, (name, _, _) in TABLES.items():
        got, skipped = read_table(domain, os.path.join(source, name))
        rows |= got
        counts[domain] = {"rows": len(got), "files": len({(d, k) for d, k, _, _ in got}),
                          "categories": len({c for _, _, _, c in got}), "blank_or_qc_dropped": skipped}
        print(f"{domain}: {counts[domain]}")
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "domain_masst.tsv.gz")
    with open(path, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9) as gz, \
            io.TextIOWrapper(gz, encoding="utf-8", newline="") as text:         # mtime 0: the same input gives the same bytes
        w = csv.writer(text, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)
        w.writerows(sorted(rows))
    with open(path, "rb") as f:
        digest = hashlib.sha256(f.read()).hexdigest()
    meta = {"source": f"https://github.com/{REPO}", "sha": SHA, "licence": "MIT (Copyright (c) 2023 Robin Schmid)",
            "built": time.strftime("%Y-%m-%d"), "columns": COLUMNS, "rows": len(rows), "domains": counts,
            "sha256": digest, "file_key": "lower-case basename without extension",
            "note": "Blank and QC files are dropped. The microbiome table is not included."}
    with open(os.path.join(out_dir, "domain_masst.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
        f.write("\n")
    print(f"{path}: {len(rows)} rows, {os.path.getsize(path) / 1e6:.2f} MB")
    return meta


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m mzdesk.gnps.domain_masst_build", description=__doc__.split("\n")[0])
    ap.add_argument("--source", help="folder with the five downloaded tables")
    ap.add_argument("--out", default=OUT_DIR, help="output folder (default: mzdesk/gnps/data)")
    ap.add_argument("--sha", help=f"must equal the pinned SHA ({SHA}); change SHA in this file to move to a new snapshot")
    args = ap.parse_args(argv)
    if args.sha and args.sha != SHA:
        ap.error("the snapshot is pinned in this file; edit SHA there")
    if args.source:
        build(args.source, args.out)
    else:
        with tempfile.TemporaryDirectory(prefix="mzdesk-domain-masst-") as directory:
            download(directory)
            build(directory, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
