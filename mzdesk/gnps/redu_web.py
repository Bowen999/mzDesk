"""Memory-bounded Pan-ReDU searches for the web app.

Only rows that can match the fastMASST hits are retained. Background counts are
computed in a second streaming pass, using the same definitions as the CLI.
"""
import collections
import csv
import gzip

from mzdesk.gnps import compound_db_check as g


def _reader(path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        yield header
        yield from reader


def find(filters, min_files, polarity):
    path = g.redu_path()
    reader = _reader(path)
    header = next(reader)
    count = 0

    def rows():
        nonlocal count
        for row in reader:
            if len(row) != len(header):
                raise ValueError("The Pan-ReDU cache has a malformed row; refresh the cache.")
            count += 1
            yield row

    result = g.find_datasets(header, rows(), filters, min_files, polarity)
    return result, count, g.redu_cache_date(path)


class TargetedRedu(g.Redu):
    """The original Redu interface, retaining only possible matches to these hits."""

    def __init__(self, bg_mode, all_rows):
        if bg_mode not in g.REDU_BG_MODES:
            raise ValueError("Unknown ReDU background.")
        path = g.redu_path()
        self.cache_date = g.redu_cache_date(path)
        self.bg_mode = bg_mode
        wanted = [r.get("USI") or "" for rows in all_rows.values() for r in rows if r.get("library") != g.REF_LIB]
        exact, noext, basenames = set(), set(), set()
        for usi in wanted:
            coll, run, _, _ = g.parse_usi(usi)
            if coll is None:
                continue
            fu = f"mzspec:{coll}:{run}".lower()
            exact.add(fu)
            noext.add(g.EXT_RE.sub("", fu))
            basenames.add((coll.lower(), g.EXT_RE.sub("", run.lower().rsplit("/", 1)[-1])))

        reader = _reader(path)
        header = next(reader)
        self.col = {c: i for i, c in enumerate(header)}
        c = self.col
        i_usi, i_ds = c["USI"], c["ATTRIBUTE_DatasetAccession"]
        i_pol, i_chr = c["IonizationSourceAndPolarity"], c["ChromatographyAndPhase"]
        self.rows, self.lower, self.noext, self.base, self.datasets = [], {}, {}, {}, set()
        self.case_dup, self.in_bg = {}, []
        self.bg_files = collections.Counter()
        pol_ok = {"pos_lc": ("pos", "alt"), "pos_lc_or_unknown": ("pos", "alt", "unk"),
                  "neg_lc": ("neg", "alt"), "neg_lc_or_unknown": ("neg", "alt", "unk"), "all": None}[bg_mode]

        def in_background(r):
            if pol_ok is None:
                return True
            gc = "gas chromatography" in r[i_chr].lower() or "electron ionization" in r[i_pol].lower()
            return g.polarity_of(r[i_pol]) in pol_ok and not gc

        seen = set()
        for r in reader:
            if len(r) != len(header):
                raise ValueError("The Pan-ReDU cache has a malformed row; refresh the cache.")
            if r[i_usi] in seen:
                continue
            seen.add(r[i_usi])
            lo = r[i_usi].lower()
            coll, run, _, _ = g.parse_usi(lo)
            if coll is not None:
                self.datasets.add(coll)
            self.datasets.add(r[i_ds].lower())
            bg = in_background(r)
            if bg:
                self.bg_files[r[i_ds]] += 1
            base = (coll, g.EXT_RE.sub("", run.rsplit("/", 1)[-1])) if coll is not None else None
            if lo not in exact and g.EXT_RE.sub("", lo) not in noext and base not in basenames:
                continue
            idx = len(self.rows)
            self.rows.append(r)
            self.in_bg.append(bg)
            if lo in self.lower:
                self.case_dup[r[i_usi]] = idx
                self.case_dup.setdefault(self.rows[self.lower[lo]][i_usi], self.lower[lo])
            else:
                self.lower[lo] = idx
            self.noext.setdefault(g.EXT_RE.sub("", lo), []).append(idx)
            if coll is not None:
                self.base.setdefault(base, []).append(idx)
        del seen
        self.n_bg_files = sum(self.bg_files.values())
        self.n_bg_ds = len(self.bg_files)
        selected = {idx for usi in wanted for idx, _, _ in [self.match(usi)] if idx is not None and self.in_bg[idx]}
        values = {k: {self.rows[i][c[k]] for i in selected} for k in g.SUMMARY_COLS}
        self.bg_val_files = {k: collections.Counter() for k in g.SUMMARY_COLS}
        self.bg_val_ds = {k: collections.defaultdict(set) for k in g.SUMMARY_COLS}
        self.bg_ds_val = {k: collections.Counter() for k in g.SUMMARY_COLS}
        reader = _reader(path)
        next(reader)
        seen = set()
        for r in reader:
            if r[i_usi] in seen:
                continue
            seen.add(r[i_usi])
            if not in_background(r):
                continue
            ds = r[i_ds]
            for k in g.SUMMARY_COLS:
                v = r[c[k]]
                if v in values[k]:
                    self.bg_val_files[k][v] += 1
                    self.bg_val_ds[k][v].add(ds)
                    self.bg_ds_val[k][(ds, v)] += 1
        self.dsw_bg = {}
        for k in g.SUMMARY_COLS:
            shares = collections.Counter()
            for (ds, value), n in self.bg_ds_val[k].items():
                shares[value] += n / self.bg_files[ds]
            self.dsw_bg[k] = {v: n / self.n_bg_ds for v, n in shares.items()} if self.n_bg_ds else {}
