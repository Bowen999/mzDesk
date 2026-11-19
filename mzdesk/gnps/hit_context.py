"""Where a fastMASST hit set occurs: hits grouped by dataset, and the sample types / organisms / body parts / diseases / domain-MASST
categories of the hit files. Used by the "Spectrum in Which Set" page (mzdesk/api/runners/gnps.py).

Sources (all joins are by file; a file counted once however many spectra of it hit, and the two USI path forms of one file are one file):
  Pan-ReDU      the redu_* columns that compound_db_check.run_query writes into the hits table (needs the local slim cache)
  domain MASST  mzdesk/gnps/data/domain_masst.tsv.gz (microbeMASST, plantMASST, tissueMASST, foodMASST, personalCareProductMASST file lists;
                github.com/robinschmid/microbe_masst, MIT; built by mzdesk/gnps/domain_masst_build.py). Offline.
  GNPS2 Explorer  title and species of a dataset (https://explorer.gnps2.org/api/datasets/<accession>); only accessions are sent.
Nothing here raises for missing data: a source that is not available leaves its cells and groups out.
"""
import collections
import concurrent.futures
import csv
import functools
import gzip
import json
import os
import re
import sys
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path

from mzdesk.gnps import compound_db_check as g

MISSING = "—"
DOMAIN_FILE = Path(__file__).resolve().parent / "data" / "domain_masst.tsv.gz"
DOMAIN_LABEL = {"microbe": "microbe", "plant": "plant", "tissue": "tissue", "food": "food", "personal_care": "personal care"}
EXPLORER = "https://explorer.gnps2.org/api/datasets/"
EXPLORER_TTL = 30 * 86400          # a title does not change often; 30 days like the lookup caches
EXPLORER_MISS_TTL = 7 * 86400      # a dataset Explorer does not know (HTTP 404 / nothing usable) is asked again after a week
EXPLORER_WORKERS = 6
EXPLORER_TIMEOUT = 8               # seconds per request
EXPLORER_BUDGET = 20               # seconds for all requests of one run; what has not answered by then is listed by accession
EXPLORER_READ = 2_000_000          # bytes of an answer that are read
MAX_GROUPS, MAX_ROWS = 200, 50     # what the page shows; the downloads hold everything
EXPLORER_MAX = MAX_GROUPS          # datasets asked per run: the ones with most hits, the same ones the page shows
NO_DATASET = "(no dataset)"
TITLE_CHARS = 90
ACCESSION = re.compile(r"^(MSV\d{9}|ST\d{6}|MTBLS\d+)$", re.I)
HIT_COLUMNS = ["search", "library", "USI", "Dataset", "Cosine", "Matching Peaks", "CompoundName", "Adduct", "verify_status", "diag_pass",
               "usi_collection", "redu_match", "redu_SampleType", "redu_NCBITaxonomy", "redu_UBERONBodyPartName", "redu_DOIDCommonName"]
# (field, ReDU column, group title); counted by hit files
REDU_FIELDS = [("sample_type", "redu_SampleType", "Sample type · Pan-ReDU"), ("organism", "redu_NCBITaxonomy", "Organism · Pan-ReDU"),
               ("body_part", "redu_UBERONBodyPartName", "Body part · Pan-ReDU"), ("disease", "redu_DOIDCommonName", "Disease · Pan-ReDU")]
DOMAIN_GROUPS = [("microbe", "Microbe · microbeMASST", "microbeMASST"), ("plant", "Plant · plantMASST", "plantMASST"),
                 ("tissue", "Tissue · tissueMASST", "tissueMASST"), ("food", "Food · foodMASST", "foodMASST"),
                 ("personal_care", "Personal care · personalCareProductMASST", "personalCareProductMASST")]


def cache_dir():
    """~/.mzdesk/cache/gnps (MZDESK_CACHE_DIR moves the cache); read at call time."""
    return Path(os.environ.get("MZDESK_CACHE_DIR", Path.home() / ".mzdesk" / "cache")) / "gnps"


# ---------- small helpers ----------
def clean(value):
    """A Pan-ReDU / table value -> text, or None when the table has no value ("missing value", blank)."""
    v = (value or "").strip()
    return None if v.lower() in ("", "missing value") else v


def organism(value):
    """'9606|Homo sapiens' -> 'Homo sapiens'"""
    v = clean(value)
    return re.sub(r"^\d+\|", "", v) if v else None


def file_id(row):
    """(dataset, file) of a hit: the dataset of the USI, and its file name without folder and extension, lower case.
    The MASST tables list a file as mzspec:<dataset>:<basename>, so this is also their join key."""
    coll, run, _, _ = g.parse_usi(row.get("USI"))
    if coll is None:
        return None
    return coll.upper(), g.EXT_RE.sub("", run.rsplit("/", 1)[-1]).strip().lower()


def dataset_url(accession):
    """A page of the dataset's own repository, or None for an accession of another kind."""
    a = (accession or "").strip().upper()
    if re.fullmatch(r"MSV\d{9}", a):
        return f"https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession={a}"
    if re.fullmatch(r"ST\d{6}", a):
        return f"https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID={a}"
    if re.fullmatch(r"MTBLS\d+", a):
        return f"https://www.ebi.ac.uk/metabolights/{a}"
    return None


# ---------- domain MASST ----------
class DomainTable:
    """{(DATASET, file): ((domain, category), ...)} kept small: {dataset: {file: code}}, the code an index into the few distinct
    (domain, category) pairs (an int, or a tuple of ints for a file listed in more than one place). About 4x smaller than a flat dict of tuples."""

    def __init__(self):
        self.pairs, self._index, self._by = [], {}, {}

    def add(self, dataset, key, domain, category):
        pair = (sys.intern(domain), sys.intern(category))
        code = self._index.get(pair)
        if code is None:
            code = self._index[pair] = len(self.pairs)
            self.pairs.append(pair)
        files = self._by.setdefault(sys.intern(dataset), {})
        old = files.get(key)
        files[key] = code if old is None else (old if isinstance(old, tuple) else (old,)) + (code,)

    def get(self, fid, default=()):
        code = self._by.get(fid[0], {}).get(fid[1])
        if code is None:
            return default
        return tuple(self.pairs[c] for c in code) if isinstance(code, tuple) else (self.pairs[code],)

    def __getitem__(self, fid):
        found = self.get(fid, None)
        if found is None:
            raise KeyError(fid)
        return found

    def __len__(self):
        return sum(len(f) for f in self._by.values())


@functools.lru_cache(maxsize=1)
def domain_table():
    """The bundled table as a DomainTable (get / [] by (DATASET, file)); empty when the file is missing."""
    table = DomainTable()
    try:
        with gzip.open(DOMAIN_FILE, "rt", encoding="utf-8", newline="") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader)
            for dataset, key, domain, category in reader:
                table.add(dataset, key, domain, category)
    except (OSError, ValueError):
        return DomainTable()
    return table


TISSUE_FIELDS = 5                   # tissueMASST ID: {species}_{body part}_{health}_{sex}_{life stage}; species names have spaces, never "_"
TISSUE_EMPTY = {"", "unknown", "missing value"}


@functools.lru_cache(maxsize=None)
def tissue_label(category):
    """'Mus musculus_caecum_healthy_male_unknown' -> 'Mus musculus · caecum · healthy · male'; unknown / missing parts are left out.
    An ID that is not five parts is returned as it is. Display only: the table keeps the original ID."""
    parts = category.split("_")
    if len(parts) != TISSUE_FIELDS:
        return category
    kept = [p.strip() for p in parts if p.strip().lower() not in TISSUE_EMPTY]
    return " · ".join(kept) or category


def readable(domain, category):
    return tissue_label(category) if domain == "tissue" else category


def domain_of(row):
    """The (domain, category) pairs of a hit file; tissueMASST IDs are shown as readable labels (see tissue_label)."""
    fid = file_id(row)
    return tuple((d, readable(d, c)) for d, c in domain_table().get(fid, ())) if fid else ()


def domain_text(pairs):
    return "; ".join(f"{DOMAIN_LABEL.get(d, d)}: {c}" for d, c in pairs)


# ---------- GNPS2 Explorer ----------
NOT_FOUND = "not_found"            # _fetch_explorer: Explorer answered and does not know the dataset (a definitive miss, not an outage)
NOT_TAXA = {"reference compound", "reference compounds", "no species", "not applicable", "not available", "n/a", "na", "none", "blank", "blanks",
            "unknown", "missing value", "qc", "quality control", "mixed", "multiple"}


def species_list(value):
    """Explorer's free-text species -> the taxa in it: split on '::' and ';', "(NCBITaxon:...)" removed, non-taxa
    ("reference compound", "no species", "Data from multiple bacteria", ...) dropped, duplicates removed, order kept."""
    out, seen = [], set()
    for part in re.split(r"::|;", str(value or "")):
        v = re.sub(r"\s+", " ", re.sub(r"\(NCBITaxon:[^)]*\)", "", part)).strip(" ,")
        k = v.lower()
        if not v or k in NOT_TAXA or k.startswith("data from") or k in seen:
            continue
        seen.add(k)
        out.append(v)
    return out


def _fetch_explorer(accession):
    """title / species / repository of one dataset; NOT_FOUND when Explorer does not know it (HTTP 404, or an answer without a title or species);
    None when it could not be asked (timeout, connection error, HTTP 5xx or other, unreadable answer)."""
    req = urllib.request.Request(EXPLORER + accession, headers=g.UA)
    try:
        with urllib.request.urlopen(req, timeout=EXPLORER_TIMEOUT) as r:
            d = json.loads(r.read(EXPLORER_READ).decode("utf-8", "replace"))
    except urllib.error.HTTPError as exc:
        return NOT_FOUND if exc.code in (404, 410) else None
    except Exception:
        return None
    if not isinstance(d, dict):
        return NOT_FOUND
    title = re.sub(r"\s+", " ", str(d.get("title") or "")).strip() or None
    species = "; ".join(species_list(d.get("species"))) or None
    if not title and not species:
        return NOT_FOUND
    return {"title": title, "species": species, "repository": str(d.get("repository") or "").strip() or None}


def _text(v):
    return v if isinstance(v, str) and v else None


def _clean_info(info):
    """An answer as kept in the cache -> {"title", "species", "repository"} of str or None; None when it is not that shape."""
    if not isinstance(info, dict):
        return None
    return {k: _text(info.get(k)) for k in ("title", "species", "repository")}


def _entry(entry, now):
    """-> ("info", info) / ("miss", None) for a cache entry that is still fresh, else None. A malformed entry is stale, never an error."""
    if not isinstance(entry, dict):
        return None
    t = entry.get("t")
    if isinstance(t, bool) or not isinstance(t, (int, float)) or t != t:
        return None
    age = now - t
    if entry.get("miss") is True:
        return ("miss", None) if 0 <= age < EXPLORER_MISS_TTL else None
    info = _clean_info(entry.get("info"))
    return ("info", info) if info and (info["title"] or info["species"]) and 0 <= age < EXPLORER_TTL else None


def _read_cache(path):
    try:
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


def _write_cache(path, fresh, now):
    """Merge into the file as it is now (another worker may have written meanwhile), dropping entries older than twice their
    time to live or malformed; a failure to write costs nothing and leaves no temporary file."""
    tmp = None
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        data = _read_cache(path)
        data.update(fresh)
        data = {a: e for a, e in data.items() if isinstance(e, dict) and isinstance(e.get("t"), (int, float)) and not isinstance(e.get("t"), bool)
                and now - e["t"] < 2 * EXPLORER_TTL}
        fd, tmp = tempfile.mkstemp(dir=path.parent, prefix="explorer-", suffix=".tmp")
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False)
        os.replace(tmp, path)
        tmp = None
    except (OSError, TypeError, ValueError):
        pass
    finally:
        if tmp:
            try:
                os.remove(tmp)
            except OSError:
                pass


def explorer_info(accessions, warnings, now=None):
    """{accession: {"title", "species", "repository"}} for the accessions Explorer knows (at most EXPLORER_MAX, in the order given).
    Answers are kept 30 days in <cache>/explorer_datasets.json, a dataset Explorer does not know 7 days. Requests share EXPLORER_BUDGET seconds;
    a dataset left unanswered (outage, too slow, over the cap) is listed by accession and counted in one warning. A dataset Explorer
    does not know is not a failure and adds no warning."""
    now = time.time() if now is None else now
    valid = [a for a in dict.fromkeys(accessions) if ACCESSION.match(a or "")]
    wanted = valid[:EXPLORER_MAX]
    path = cache_dir() / "explorer_datasets.json"
    cache = _read_cache(path)
    out, todo = {}, []
    for a in wanted:
        hit = _entry(cache.get(a), now)
        if hit is None:
            todo.append(a)
        elif hit[0] == "info":
            out[a] = hit[1]
    unanswered = len(valid) - len(wanted)
    if todo:
        deadline = time.monotonic() + EXPLORER_BUDGET
        state = {"ok": 0, "outage": 0}

        def one(acc):
            if time.monotonic() >= deadline or (state["outage"] >= EXPLORER_WORKERS and not state["ok"]):     # over budget, or Explorer is down
                return acc, None
            got = _fetch_explorer(acc)
            if got is None:
                state["outage"] += 1
            else:
                state["ok"] += 1                                  # a 404 is an answer: the server is up
            return acc, got

        fresh = {}
        pool = concurrent.futures.ThreadPoolExecutor(max_workers=EXPLORER_WORKERS)
        try:
            futures = [pool.submit(one, a) for a in todo]
            done, _ = concurrent.futures.wait(futures, timeout=EXPLORER_BUDGET)
        finally:
            pool.shutdown(wait=False, cancel_futures=True)        # requests still running end by their own timeout; their answers are ignored
        for f in done:
            try:
                acc, got = f.result()
            except Exception:
                continue
            if got == NOT_FOUND:
                fresh[acc] = {"t": now, "miss": True}
            elif got is not None and (info := _clean_info(got)) is not None:
                out[acc] = info
                fresh[acc] = {"t": now, "info": info}
            else:
                unanswered += 1
        unanswered += len(futures) - len(done)
        if fresh:
            _write_cache(path, fresh, now)
    if unanswered:
        warnings.append(f"GNPS2 Explorer gave no title for {unanswered} of {len(valid)} datasets (not reachable in time, or over the {EXPLORER_MAX} asked); "
                        "they are listed by accession only.")
    return out


# ---------- reading the hits table ----------
def read_hits(path):
    """Every passing hit of the hits table (the columns the page needs), each with its file and Domain MASST categories."""
    rows = []
    if not path or not os.path.isfile(path):
        return rows
    with open(path, encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("pass_min_peaks") != "1":
                continue
            row = {c: r.get(c) or "" for c in HIT_COLUMNS}
            row["file"] = file_id(row)
            rows.append(row)
    return rows


def _num(row, key):
    try:
        return float(row.get(key))
    except (TypeError, ValueError):
        return 0.0


def redu_value(row, column):
    if row.get("redu_match") not in g.MATCHED_LEVELS:
        return None
    return organism(row.get(column)) if column == "redu_NCBITaxonomy" else clean(row.get(column))


def _top(rows, column):
    """Most frequent value among the files of these rows (a file once), or None."""
    seen, count = set(), collections.Counter()
    for r in rows:
        v = redu_value(r, column)
        if v is None or (r["file"], v) in seen:
            continue
        seen.add((r["file"], v))
        count[v] += 1
    return min(count, key=lambda v: (-count[v], v)) if count else None


def _is_library(row):
    return row.get("library") == g.REF_LIB


def _dataset_of(row):
    return (row.get("Dataset") or row.get("usi_collection") or "").strip()


def dataset_rows(rows, joined):
    """One record per (search, dataset) in the order of the page: most hits first, then best cosine. The reference library is one record."""
    by = collections.defaultdict(list)
    for r in rows:
        by[(r["search"], "" if _is_library(r) else _dataset_of(r))].append(r)
    out = []
    for (search, dataset), rs in by.items():
        lib = _is_library(rs[0])
        files = {r["file"] for r in rs if r["file"]}
        levels = {r.get("redu_match") for r in rs}
        if not joined or lib:
            status = ""
        elif levels & g.MATCHED_LEVELS:
            status = "matched"
        elif "dataset_only" in levels:
            status = "dataset only"
        elif levels == {"no_dataset"}:
            status = "not in Pan-ReDU"
        else:
            status = ""
        out.append({"search": search, "dataset": dataset, "library": lib, "rows": sorted(rs, key=lambda r: -_num(r, "Cosine")),
                    "hits": len(rs), "files": len(files), "best": max(_num(r, "Cosine") for r in rs), "redu": status,
                    "files_matched": len({r["file"] for r in rs if r.get("redu_match") in g.MATCHED_LEVELS and r["file"]}),
                    "sample_type": _top(rs, "redu_SampleType"), "organism": _top(rs, "redu_NCBITaxonomy"),
                    "body_part": _top(rs, "redu_UBERONBodyPartName")})
    out.sort(key=lambda d: (d["library"], -d["hits"], -d["best"], d["dataset"]))   # repository datasets first, the library last
    return out


# ---------- the page's blocks ----------
def count(n, word):
    """'1 hit', '2 hits'"""
    return f"{n} {word}" if n == 1 else f"{n} {word}s"


def _title(accession, info):
    if not accession:
        return NO_DATASET
    t = (info.get(accession.upper()) or {}).get("title")
    if not t:
        return accession
    return f"{accession} · {t if len(t) <= TITLE_CHARS else t[:TITLE_CHARS - 1].rstrip() + '…'}"


def _species_text(value, most=3):
    names = species_list(value)
    return ", ".join(names[:most]) + (f" +{len(names) - most}" if len(names) > most else "")


def _meta(d, info, joined):
    """The line of a group: counts, best cosine, what Pan-ReDU (or Explorer) says about the samples; status words last."""
    species = _species_text((info.get(d["dataset"].upper()) or {}).get("species"))
    if d["library"]:
        names = {r.get("CompoundName") for r in d["rows"] if r.get("CompoundName")}
        return f"{count(d['hits'], 'hit')} · {count(len(names), 'name')} · best cosine {d['best']:.3f}"
    parts = [count(d["hits"], "hit"), count(d["files"], "file"), f"best cosine {d['best']:.3f}"]
    if joined:
        parts += [v for v in (d["sample_type"], d["organism"], d["body_part"]) if v]
    if species and not d["organism"]:
        parts.append(species)
    if joined:
        if d["redu"] == "not in Pan-ReDU":
            parts.append("not in Pan-ReDU")
        elif d["redu"] == "dataset only":
            parts.append("no Pan-ReDU file match")
    return " · ".join(parts)


def hits_block(datasets, info, joined, total_hits, cell):
    """The grouptable data of "Passing hits by dataset". cell(row) -> the list of cells of one hit. When the page shows fewer than
    there are, "note" (the line under the title) says so; it is not repeated as a warning or in the help."""
    groups = []
    for d in datasets[:MAX_GROUPS]:
        group = {"title": "GNPS reference library" if d["library"] else _title(d["dataset"], info),
                 "meta": _meta(d, info, joined), "rows": [cell(r) for r in d["rows"][:MAX_ROWS]]}
        url = None if d["library"] else dataset_url(d["dataset"])
        if url:
            group["url"] = url
        groups.append(group)
    capped = []
    if len(datasets) > MAX_GROUPS:
        capped.append(f"the first {MAX_GROUPS} of {len(datasets)} datasets")
    if any(len(d["rows"]) > MAX_ROWS for d in datasets[:MAX_GROUPS]):
        capped.append(f"the {MAX_ROWS} best hits of each dataset")
    block = {"groups": groups, "count": total_hits}
    if capped:
        block["note"] = "The page shows " + " and ".join(capped) + ". The downloads hold every hit."
    return block


def _share(n, total):
    if not total:
        return MISSING
    pct = 100 * n / total
    return "<1%" if 0 < pct < 1 else f"{pct:.0f}%"


def distribution(rows, joined, info):
    """-> (groups, table, note): groups for the grouptable "Where it occurs", flat rows (field, value, files, datasets, share) for the download,
    and one line saying how much of the hit files / datasets each source covers ("" when there is nothing to show).
    Repository hits only; a file counts once per value and a dataset once per value; missing values are left out."""
    rows = [r for r in rows if not _is_library(r) and r["file"]]
    all_files = {r["file"] for r in rows}
    all_datasets = {r["file"][0] for r in rows}
    fields = []        # (title, {value: files}, {value: datasets}, unit)
    covers = []        # (source, files covered)

    def tally(pairs):
        files, datasets = collections.defaultdict(set), collections.defaultdict(set)
        for fid, value in pairs:
            files[value].add(fid)
            datasets[value].add(fid[0])
        return files, datasets

    if joined:
        covers.append(("Pan-ReDU", len({r["file"] for r in rows if r.get("redu_match") in g.MATCHED_LEVELS})))
        for _, column, title in REDU_FIELDS:
            files, datasets = tally((r["file"], v) for r in rows if (v := redu_value(r, column)) is not None)
            fields.append((title, files, datasets, "file"))
    species = {d: species_list((info.get(d) or {}).get("species")) for d in all_datasets}
    by_species = collections.defaultdict(set)
    for r in rows:
        for sp in species.get(r["file"][0], ()):
            by_species[sp].add(r["file"])
    sp_datasets = {v: {f[0] for f in fs} for v, fs in by_species.items()}
    fields.append(("Dataset species · GNPS2 Explorer", by_species, sp_datasets, "dataset"))
    zero = []
    for domain, title, label in DOMAIN_GROUPS:
        pairs = [(r["file"], c) for r in rows for d, c in domain_of(r) if d == domain]
        files, datasets = tally(pairs)
        fields.append((title, files, datasets, "file"))
        n = len({f for f, _ in pairs})
        (covers.append((label, n)) if n else zero.append(label))
    groups, table = [], []
    for title, files, datasets, unit in fields:
        if not files:
            continue
        key = (lambda v: (-len(datasets[v]), -len(files[v]), v)) if unit == "dataset" else (lambda v: (-len(files[v]), -len(datasets[v]), v))
        order = sorted(files, key=key)
        top = " · ".join(f"{v} {len(datasets[v]) if unit == 'dataset' else len(files[v])}" for v in order[:3])
        groups.append({"title": title, "meta": top,
                       "rows": [[v, str(len(files[v])), str(len(datasets[v])), _share(len(files[v]), len(all_files))] for v in order]})
        table += [[title, v, len(files[v]), len(datasets[v]), round(len(files[v]) / len(all_files), 4)] for v in order]
    note = ""
    if all_files:
        parts = [f"{src} {n} of {len(all_files)}" for src, n in covers]
        if parts:
            parts[0] += " hit file" if len(all_files) == 1 else " hit files"             # the unit once, on the first
        parts.append(f"GNPS2 Explorer {sum(1 for d in all_datasets if species.get(d))} of {count(len(all_datasets), 'dataset')}")
        if zero:
            parts.append("none in " + ", ".join(zero))
        note = "Coverage: " + " · ".join(parts)
    return groups, table, note
