#!/usr/bin/env python3
"""Score a SIRIUS / CSI:FingerID / MSNovelist run against CASMI 2022 truth.

Usage:
  python3 -I scripts/04_score.py --summaries results/pilot/sirius_summaries \
      --truth results/truth.tsv --out results/pilot/score [--ids 1,2,3 | --ms-dir ms]

Denominator = all submitted challenges (--ids, else --ms-dir *.ms, else
qc_results_present.tsv in the summaries folder). A challenge with no output
counts as "not found", never dropped.

Ranking rules
  formula   : neutral molecularFormula, SiriusScore. Rows of the same neutral formula
              (different adducts) are collapsed, keeping the best-scored row.
  CSI       : PubChem candidates, ordered by CSI:FingerIDScore (= SIRIUS structurePerIdRank).
              De-duplicated by InChIkey2D, keeping the best-scored row.
  MSNovelist: same, on denovo_structure_identifications_all.tsv (SIRIUS export order,
              i.e. CSI:FingerIDScore, not the MSNovelist ModelScore).
  Ties      : reported rank is the BEST (optimistic) one:
              rank = 1 + number of distinct candidates with a strictly better score.
              *_n_tied = number of other candidates sharing the truth's score.
              The listed (file-order) rank is also written, as *_rank_listed.
"""
import argparse
import csv
import re
import sys
from pathlib import Path

import pandas as pd

TOPS = [1, 5, 10, 20]
ID_RE = re.compile(r"_(challenge_\d+)_challenge")
ID_RE_B = re.compile(rb"_(challenge_\d+)_challenge")


def norm_adduct(a):
    """Normalise adduct strings so SIRIUS and CASMI notation compare equal."""
    a = str(a).replace(" ", "")
    a = a.replace("+H3N+H", "+NH4").replace("+CH2O2-H", "+FA-H")
    return a


def raw_challenge_ids(p):
    """Independent, line-based challenge-id set (no CSV parsing at all)."""
    ids = set()
    with open(p, "rb") as fh:
        for line in fh:
            m = ID_RE_B.search(line)
            if m:
                ids.add(m.group(1).decode())
    return ids


def read(summ, name, usecols, keep_ids=None, chunksize=200_000):
    """Read a SIRIUS summary TSV.

    quoting=QUOTE_NONE: compound names can contain a stray '"', which the default
    CSV quoting would treat as an opening quote and swallow following lines
    (silently dropping whole challenges). Only `usecols` are kept, in chunks, and
    rows are filtered to `keep_ids` so the 1.5 GB structure file stays small.
    Parsed challenge set is checked against a raw line-based scan.
    """
    p = Path(summ) / name
    if not p.exists():
        return pd.DataFrame()
    cols = list(dict.fromkeys(usecols + ["mappingFeatureId"]))
    parts = []
    for ch in pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE,
                          usecols=cols, chunksize=chunksize):
        ch["cid"] = ch["mappingFeatureId"].str.extract(ID_RE)[0]
        if keep_ids is not None:
            ch = ch[ch["cid"].isin(keep_ids)]
        parts.append(ch)
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=cols + ["cid"])
    raw = raw_challenge_ids(p)
    parsed = set(df["cid"].dropna()) if keep_ids is None else set(df["cid"].dropna()) | (raw - keep_ids)
    assert parsed == raw, (f"{name}: parsed challenges != raw line scan "
                           f"(missing {sorted(raw - parsed)[:5]}, extra {sorted(parsed - raw)[:5]})")
    print(f"{name}: {len(raw)} challenges with rows (raw line scan == parsed), "
          f"{len(df)} rows kept", flush=True)
    return df


def rank_by_score(items, truth_key):
    """items: list of (key, score), best first. Returns (optimistic, listed, n_tied)
    for truth_key, or (None, None, 0) if absent. Keys must already be unique."""
    keys = [k for k, _ in items]
    if truth_key not in keys:
        return None, None, 0
    s = dict(items)[truth_key]
    better = sum(1 for _, x in items if x > s)
    tied = sum(1 for k, x in items if x == s and k != truth_key)
    return better + 1, keys.index(truth_key) + 1, tied


def struct_rank(df, truth_key, score_col):
    """df: rows of one challenge. De-dup by InChIkey2D keeping best row."""
    if df.empty:
        return None, None, 0, 0
    df = df.assign(_r=pd.to_numeric(df["structurePerIdRank"]),
                   _s=pd.to_numeric(df[score_col]))
    df = df.sort_values("_r", kind="stable").drop_duplicates("InChIkey2D", keep="first")
    items = list(zip(df["InChIkey2D"], df["_s"]))
    r, rl, t = rank_by_score(items, truth_key)
    return r, rl, t, len(df)


def score_challenge(n, truth, formula, csi, nov):
    cid = f"challenge_{n:03d}"
    t = truth[truth["challenge_id"] == str(n)]
    if t.empty:
        sys.exit(f"challenge {n} not in truth")
    t = t.iloc[0]
    row = dict(challenge_id=n, priority=t["priority"], polarity=t["polarity"],
               true_formula=t["formula"], true_adduct=t["adduct"])
    notes = []

    # formula
    f = formula[formula["cid"] == cid]
    if f.empty:
        notes.append("no formula output")
        row.update(formula_rank="", formula_rank_listed="", formula_n_tied=0,
                   top1_formula="", top1_adduct="", adduct_top1_correct="no")
    else:
        f = f.assign(_r=pd.to_numeric(f["formulaRank"]), _s=pd.to_numeric(f["SiriusScore"]))
        f = f.sort_values("_r", kind="stable")
        top = f.iloc[0]
        f = f.drop_duplicates("molecularFormula", keep="first")  # collapse adducts
        r, rl, nt = rank_by_score(list(zip(f["molecularFormula"], f["_s"])), t["formula"])
        row.update(formula_rank=r or "", formula_rank_listed=rl or "", formula_n_tied=nt,
                   top1_formula=top["molecularFormula"], top1_adduct=top["adduct"],
                   adduct_top1_correct="yes" if norm_adduct(top["adduct"]) == t["adduct"] else "no")
        if r is None:
            notes.append("true formula not in list")
        elif nt:
            notes.append(f"formula tie at rank {r} ({nt} other)")

    # CSI:FingerID (PubChem)
    c = csi[csi["cid"] == cid]
    r, rl, nt, n_c = struct_rank(c, t["inchikey14"], "CSI:FingerIDScore")
    row.update(csi_rank=r or "", csi_rank_listed=rl or "", csi_n_tied=nt)
    if c.empty:
        notes.append("no CSI output")
        row.update(csi_top1_inchikey14="", csi_top1_name="", confidence_exact="", confidence_approx="")
    else:
        top = c.assign(_r=pd.to_numeric(c["structurePerIdRank"])).sort_values("_r").iloc[0]
        row.update(csi_top1_inchikey14=top["InChIkey2D"], csi_top1_name=top["name"],
                   confidence_exact=top["ConfidenceScoreExact"],
                   confidence_approx=top["ConfidenceScoreApproximate"])
        if r is None:
            notes.append("truth not in CSI candidates")
        elif nt:
            notes.append(f"CSI tie at rank {r} ({nt} other)")
    row.update(n_csi_candidates=n_c, truth_in_pubchem_candidates="yes" if r else "no")

    # MSNovelist
    d = nov[nov["cid"] == cid]
    r, rl, nt, n_d = struct_rank(d, t["inchikey14"], "CSI:FingerIDScore")
    row.update(novelist_rank=r or "", novelist_rank_listed=rl or "", novelist_n_tied=nt,
               n_novelist_candidates=n_d)
    if d.empty:
        notes.append("no MSNovelist output")
    elif r is None:
        notes.append("truth not in MSNovelist candidates")
    elif nt:
        notes.append(f"MSNovelist tie at rank {r} ({nt} other)")
    row["notes"] = "; ".join(notes)
    return row


def summarize(pc):
    groups = {"all": pc,
              "Priority": pc[pc["priority"] == "Priority"],
              "Bonus": pc[pc["priority"] == "Bonus"],
              "pos": pc[pc["polarity"] == "positive"],
              "neg": pc[pc["polarity"] == "negative"]}
    methods = {"formula": "formula_rank", "CSI:FingerID": "csi_rank", "MSNovelist": "novelist_rank"}
    tie_cols = {"formula": "formula_n_tied", "CSI:FingerID": "csi_n_tied", "MSNovelist": "novelist_n_tied"}
    rows = []
    for gname, g in groups.items():
        for m, col in methods.items():
            rk = pd.to_numeric(g[col], errors="coerce")  # blank -> NaN = not found
            n = len(g)
            o = dict(method=m, subset=gname, N=n)
            for k in TOPS:
                c = int((rk <= k).sum())
                o[f"top{k}"] = c
                o[f"top{k}_pct"] = round(100 * c / n, 1) if n else ""
            c = int(rk.notna().sum())
            o["any_rank"] = c
            o["any_rank_pct"] = round(100 * c / n, 1) if n else ""
            # found challenges whose truth shares its score with another candidate
            o["n_ties_at_truth"] = int(((g[tie_cols[m]].astype(int) > 0) & rk.notna()).sum())
            rows.append(o)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summaries", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ids", help="comma-separated challenge numbers, e.g. 1,2,3")
    ap.add_argument("--ms-dir", help="folder of challenge_NNN.ms files = submitted set")
    a = ap.parse_args()

    truth = pd.read_csv(a.truth, sep="\t", dtype=str, keep_default_na=False)
    if a.ids:
        ids = [int(x) for x in a.ids.split(",")]
    elif a.ms_dir:
        ids = sorted(int(re.search(r"challenge_(\d+)", p.name).group(1)) for p in Path(a.ms_dir).glob("challenge_*.ms"))
    else:
        qc = Path(a.summaries) / "qc_results_present.tsv"
        if not qc.exists():
            sys.exit("no --ids / --ms-dir and no qc_results_present.tsv")
        ids = sorted(int(re.search(r"challenge_(\d+)", x).group(1))
                     for x in pd.read_csv(qc, sep="\t", dtype=str)["challenge"])
    if not ids:
        sys.exit("empty challenge set")

    keep = {f"challenge_{n:03d}" for n in ids}
    formula = read(a.summaries, "formula_identifications_all.tsv",
                   ["formulaRank", "molecularFormula", "adduct", "SiriusScore"], keep)
    csi = read(a.summaries, "structure_identifications_all.tsv",
               ["structurePerIdRank", "ConfidenceScoreExact", "ConfidenceScoreApproximate",
                "CSI:FingerIDScore", "InChIkey2D", "name"], keep)
    nov = read(a.summaries, "denovo_structure_identifications_all.tsv",
               ["structurePerIdRank", "CSI:FingerIDScore", "InChIkey2D"], keep)
    for name, df in [("formula", formula), ("CSI", csi), ("MSNovelist", nov)]:
        if df.empty:
            df["cid"] = pd.Series(dtype=str)
        extra = set(df["cid"].dropna()) - {f"challenge_{n:03d}" for n in ids}
        if extra:
            print(f"note: {name} output has challenges outside the id set (ignored): {sorted(extra)}")

    cols = ["challenge_id", "priority", "polarity", "true_formula", "true_adduct",
            "formula_rank", "formula_rank_listed", "formula_n_tied", "top1_formula", "top1_adduct",
            "adduct_top1_correct", "csi_rank", "csi_rank_listed", "csi_n_tied",
            "csi_top1_inchikey14", "csi_top1_name", "confidence_exact", "confidence_approx",
            "n_csi_candidates", "truth_in_pubchem_candidates", "novelist_rank", "novelist_rank_listed",
            "novelist_n_tied", "n_novelist_candidates", "notes"]
    pc = pd.DataFrame([score_challenge(n, truth, formula, csi, nov) for n in ids])[cols]
    sm = summarize(pc)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    pc.to_csv(out / "per_challenge.tsv", sep="\t", index=False)
    sm.to_csv(out / "summary.tsv", sep="\t", index=False)

    print(f"Challenges scored: {len(pc)}  ({ids})")
    print("Ranks: best (optimistic) rank under ties; blank = not found.\n")
    print(sm.to_string(index=False))


if __name__ == "__main__":
    main()
