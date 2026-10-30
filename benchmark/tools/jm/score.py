#!/usr/bin/env python3
"""Score JESTR / MIST (and CSI:FingerID on the same candidate files) against CASMI 2022 truth.

Rules (same as casmi2022/sirius_run/scripts/04_score.py):
  * hit = InChIKey first block (2D) of a scored candidate equals the truth
  * rank = 1 + number of candidates with a strictly better score (optimistic under ties);
    rank_pess (truth placed after all tied candidates) is reported too
  * denominator = every challenge of the subset; no output / truth not scored = not found
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

TOPS = [1, 5, 10, 20]
METHODS = {  # name -> (score dir key, column)
    "CSI:FingerID": ("candidates", "csi_score"),
    "JESTR": ("jestr", "score"),
    "MIST-FP": ("mist", "score_fp"),
    "MIST-contrastive": ("mist", "score_contrast"),
}


def rank_of(df: pd.DataFrame, col: str, truth: str):
    s = pd.to_numeric(df[col], errors="coerce")
    ok = s.notna()
    keys, s = df["ikey14"][ok], s[ok]
    hit = keys == truth
    if not hit.any():
        return None, None, 0, int(ok.sum())
    t = s[hit].max()
    better = int((s > t).sum())
    tied = int((s == t).sum()) - 1
    return better + 1, better + 1 + tied, tied, int(ok.sum())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True, help="jm_run dir (inputs/, jestr/, mist/)")
    ap.add_argument("--truth", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    run = Path(a.run_dir)
    chal = pd.read_csv(run / "inputs/challenges.tsv", sep="\t")
    truth = pd.read_csv(a.truth, sep="\t", dtype=str, keep_default_na=False)
    truth["challenge"] = "challenge_" + truth["challenge_id"].astype(int).map("{:03d}".format)
    truth = truth.set_index("challenge")
    dirs = {"candidates": run / "inputs/candidates", "jestr": run / "jestr/scores", "mist": run / "mist/scores"}

    rows = []
    for c in chal.itertuples(index=False):
        t = truth.loc[c.challenge]
        row = dict(challenge=c.challenge, challenge_id=c.challenge_id, priority=t["priority"],
                   polarity=c.polarity, true_adduct=t["adduct"], true_formula=t["formula"],
                   n_candidates=c.n_candidates)
        cand = pd.read_csv(dirs["candidates"] / f"{c.challenge}.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
        row["truth_in_candidates"] = bool((cand["ikey14"] == t["inchikey14"]).any())
        cache = {"candidates": cand}
        for m, (key, col) in METHODS.items():
            if key not in cache:
                p = dirs[key] / f"{c.challenge}.tsv.gz"
                cache[key] = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False) if p.exists() else None
            df = cache[key]
            if df is None or col not in df:
                r = (None, None, 0, 0)
            else:
                r = rank_of(df, col, t["inchikey14"])
            row[f"{m}_rank"], row[f"{m}_rank_pess"], row[f"{m}_n_tied"], row[f"{m}_n_scored"] = r
        rows.append(row)
    pc = pd.DataFrame(rows)

    subsets = {"all": pc, "Priority": pc[pc.priority == "Priority"], "Bonus": pc[pc.priority == "Bonus"],
               "pos": pc[pc.polarity == "positive"], "neg": pc[pc.polarity == "negative"]}
    summ = []
    for sname, g in subsets.items():
        n = len(g)
        for m in METHODS:
            rk = pd.to_numeric(g[f"{m}_rank"], errors="coerce")
            rp = pd.to_numeric(g[f"{m}_rank_pess"], errors="coerce")
            o = dict(method=m, subset=sname, N=n,
                     n_with_output=int((g[f"{m}_n_scored"] > 0).sum()),
                     truth_in_candidates=int(g.truth_in_candidates.sum()))
            for k in TOPS:
                o[f"top{k}"] = int((rk <= k).sum())
                o[f"top{k}_pct"] = round(100 * o[f"top{k}"] / n, 1) if n else np.nan
            o["top1_pess_pct"] = round(100 * (rp <= 1).sum() / n, 1) if n else np.nan
            o["truth_scored"] = int(rk.notna().sum())
            o["median_rank_found"] = float(rk.median()) if rk.notna().any() else np.nan
            summ.append(o)
    sm = pd.DataFrame(summ)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    pc.to_csv(out / "per_challenge.tsv", sep="\t", index=False)
    sm.to_csv(out / "summary.tsv", sep="\t", index=False)

    show = ["method", "subset", "N", "n_with_output", "truth_in_candidates", "truth_scored",
            "top1_pct", "top5_pct", "top10_pct", "top20_pct", "top1_pess_pct"]
    print("Hit = InChIKey 2D block; rank optimistic under ties; no output = not found.")
    print(sm[sm.subset.isin(["all", "pos", "neg"])][show].to_string(index=False))


if __name__ == "__main__":
    main()
