#!/usr/bin/env python3
"""Build the CASMI 2022 summary workbook: truth + top 10 per challenge for every software / situation.

  python3 benchmark/casmi2022/mzdesk_ablation_20261030/write_xlsx.py [--target FILE] [--replace] [--no-check]

Tabs: Summary, Truth, then one wide tab per situation (mzDesk main, mzDesk ablation, single software default settings, comparators).
Reads only the existing outputs (replay/, standalone/, scores.json, jm_run/, sirius_run/) and never changes them; the existing scorers' rules are
reused (IK14 match; ties -> worst rank; challenge without output stays in the denominator). The truth ranks written here are recomputed from the
candidate lists and then compared with scores.json (self-check at the end). Needs pandas, RDKit and openpyxl; does not import mzdesk.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import pickle
import re
import statistics
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import pandas as pd
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import score as S                    # noqa: E402  (only its pure functions are used; importing runs nothing)
import write_report as W             # noqa: E402

ROOT = HERE.parents[2]
CASMI = ROOT / "benchmark/casmi2022"
JM = CASMI / "jm_run"
SIRIUS = CASMI / "sirius_run/results"
CFMID = CASMI / "cfmid_casmi2022_entry"
SUMM = SIRIUS / "sirius_summaries"
TARGET = ROOT / "benchmark/result/CASMI2022_mzdesk_full500.xlsx"
SCRATCH = Path("/private/tmp/claude-501/-Users-bowen-Documents-reasearch-01-Core-Projects-mzDesk/46692548-9c05-45e0-abca-108eb47d0813/scratchpad")
IDS = list(range(1, 501))

# (tab name, group, kind, key). kind: "struct" or "formula". key: replay config / standalone config / comparator / "sirius_formula".
TABS = [
    ("mzDesk full", "mzDesk main", "struct", "full"),
    ("mzDesk leak-free", "mzDesk main", "struct", "leakfree_full"),
    ("Abl all_layers", "mzDesk ablation", "struct", "all_layers"),
    ("Abl library_only", "mzDesk ablation", "struct", "library_only"),
    ("Abl msfinder_only", "mzDesk ablation", "struct", "msfinder_only"),
    ("Abl frigid_only", "mzDesk ablation", "struct", "frigid_only"),
    ("Abl no_library", "mzDesk ablation", "struct", "no_library"),
    ("Abl no_msfinder", "mzDesk ablation", "struct", "no_msfinder"),
    ("Abl no_frigid", "mzDesk ablation", "struct", "no_frigid"),
    ("Abl buddy_only (formula)", "mzDesk ablation", "formula", "buddy_only"),
    ("Default library", "Single software, default settings", "struct", "default_library"),
    ("Default MS-FINDER", "Single software, default settings", "struct", "default_msfinder"),
    ("Default FRIGID", "Single software, default settings", "struct", "default_frigid"),
    ("Default BUDDY (formula)", "Single software, default settings", "formula", "default_buddy"),
    ("CSI FingerID", "Comparators", "struct", "CSI:FingerID"),
    ("MSNovelist", "Comparators", "struct", "MSNovelist"),
    ("JESTR", "Comparators", "struct", "JESTR"),
    ("MIST-contrastive", "Comparators", "struct", "MIST-contrastive"),
    ("CFM-ID (CASMI entry)", "Comparators", "struct", "cfmid_entry"),
    ("SIRIUS formula", "Comparators", "formula", "sirius_formula"),
    ("CFM-ID formula (CASMI entry)", "Comparators", "formula", "cfmid_entry_formula"),
]
# The CFM-ID team's own CASMI 2022 submission: one answer per challenge, so only Top-1 exists; the other cells are "n/a" (not a measured 0).
TOP1_ONLY = {"cfmid_entry", "cfmid_entry_formula"}
NA = "n/a"
POSITIVE_ONLY = set(S.POSITIVE_ONLY) | {"MIST-contrastive"}
NOT_SUPPORTED = "negative: not supported"


# ------------------------------------------------------------------ helpers
def read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE)


def num(x):
    return S.number(x)


def rank_scored(entries, truth):
    """entries: [(key, score, payload)] one per structure. -> (n, worst rank or None, display list best first).
    Same rule as score.scored_rank: rank = (number strictly better) + (number tied with the truth), i.e. the worst rank on ties. The truth is
    placed last among its ties so that its position in the displayed list equals the reported rank."""
    n = len(entries)
    target = next((s for k, s, _ in entries if k == truth), None)
    order = sorted(range(n), key=lambda i: (-entries[i][1], entries[i][0] == truth, i))
    display = [entries[i] for i in order]
    if target is None:
        return n, None, display
    worst = sum(s > target for _, s, _ in entries) + sum(s == target for _, s, _ in entries)
    return n, worst, display


def mk(n, rank, display, status="", cols=None):
    """Result for one challenge: top 10 items as (name, key, score, extra)."""
    return {"n": n, "rank": rank, "status": status, "top": display[:10]}


def counts_key(formula):
    return tuple(sorted(S.counts(formula).items()))


# ------------------------------------------------------------------ truth
def load_truth():
    official = pd.read_csv(CASMI / "casmi2022_ground_truth.csv", dtype=str, keep_default_na=False).set_index("challenge_id")
    old = read_tsv(SIRIUS / "truth.tsv").set_index("challenge_id")
    return official, old


def parentmass(n):
    for line in (CASMI / f"sirius_run/ms/challenge_{n:03d}.ms").read_text().splitlines():
        if line.startswith(">parentmass"):
            return float(line.split()[1])
    raise ValueError(n)


# ------------------------------------------------------------------ mzDesk replay configurations
def mzdesk_result(cfg, n, truth, truth_formula):
    name = f"challenge_{n:03d}"
    run = json.loads((HERE / "replay" / cfg / f"{name}.json").read_text())
    res = run["result"]
    if cfg == "buddy_only":
        return buddy_list_result(res["formulas"], truth_formula, res.get("error"))
    cands = [c for c in res["candidates"] if c["layer"] != "formula_only"]
    try:
        row = S.mzdesk_rank(res, truth)
    except AssertionError:
        row = S.fallback_rank(res, truth)
    display = [(c.get("name") or "", c.get("ik14") or "", c["final_score"], c) for c in cands]
    worst, pos = row["worst"], next((i for i, c in enumerate(cands) if c.get("ik14") == truth), None)
    if worst is not None and pos is not None and worst - 1 > pos:     # tie: put the truth at its worst position
        item = display.pop(pos)
        display.insert(worst - 1, item)
    status = "" if cands else "no output"
    return mk(len(cands), worst, display, status)


def buddy_list_result(formulas, truth_formula, error=None):
    """Distinct formulas in list order (as score.formula_rank); rank = position of the truth formula."""
    seen, distinct = set(), []
    for f in formulas:
        key = counts_key(f["formula"])
        if key not in seen:
            seen.add(key)
            distinct.append(f)
    target = counts_key(truth_formula)
    pos = next((i for i, f in enumerate(distinct, 1) if counts_key(f["formula"]) == target), None)
    display = [(f["formula"], f.get("adduct") or "", f.get("prob"), f) for f in distinct]
    return mk(len(distinct), pos, display, "" if distinct else "no output")


# ------------------------------------------------------------------ single tools, default settings
def standalone_result(cfg, n, truth, truth_formula, pol, truth_adduct):
    name = f"challenge_{n:03d}"
    tool = cfg.replace("default_", "")
    wrong = truth_adduct != S.DEFAULT_ADDUCT[pol]
    if tool == "frigid" and pol != "positive":
        return mk(0, None, [], NOT_SUPPORTED)
    if tool != "library" and wrong:
        return mk(0, None, [], f"not run: true adduct not default (true {truth_adduct}, default {S.DEFAULT_ADDUCT[pol]})")
    path = HERE / "standalone" / tool / f"{name}.json"
    if not path.exists():
        err = HERE / "standalone" / tool / f"{name}.error.json"
        return mk(0, None, [], "no output: run error" if err.exists() else "no output: not run")
    res = json.loads(path.read_text())
    if tool == "library":
        items = res["compounds"]
    elif tool == "buddy":
        if res["status"] != "ok":
            return mk(0, None, [], f"no output (BUDDY status: {res['status']})")
        return buddy_list_result(res["formulas"], truth_formula)
    else:
        if res["status"] != "done":
            return mk(0, None, [], f"no output (status: {res['status']})")
        items = res["candidates"]
    entries, seen = [], set()
    for i, c in enumerate(items):                       # score.candidate_pairs: one entry per IK14 (first), structure-less ones keep their slot
        key = c.get("ik14") or f"?{i}"
        if key not in seen:
            seen.add(key)
            entries.append((key, c["score"], c))
    n_, worst, display = rank_scored(entries, truth)
    shown = [(c.get("name") or "", c.get("ik14") or "", s, c) for _, s, c in display]
    return mk(n_, worst, shown, "" if entries else "no output")


# ------------------------------------------------------------------ comparators
def smiles_map(n):
    p = JM / "inputs/candidates" / f"challenge_{n:03d}.tsv.gz"
    if not p.exists():
        return {}
    df = read_tsv(p)
    out = {}
    for k, s in zip(df["ikey14"], df["smiles"]):
        out.setdefault(k, s)
    return out


def jm_result(cfg, n, truth, pol):
    name = f"challenge_{n:03d}"
    if cfg == "CSI:FingerID":
        path, kc, sc = JM / "inputs/candidates" / f"{name}.tsv.gz", "ikey14", "csi_score"
    elif cfg == "JESTR":
        path, kc, sc = JM / "jestr/scores" / f"{name}.tsv.gz", "ikey14", "score"
    else:
        if pol != "positive":
            return mk(0, None, [], NOT_SUPPORTED)
        path, kc, sc = JM / "mist/scores" / f"{name}.tsv.gz", "ikey14", "score_contrast"
    if not path.exists():
        return mk(0, None, [], "no output")
    df = read_tsv(path)
    entries = [(str(k), num(s), None) for k, s in zip(df[kc], df[sc]) if num(s) is not None]
    assert len({k for k, _, _ in entries}) == len(entries), (cfg, n)
    n_, worst, display = rank_scored(entries, truth)
    ref = S.scored_rank([(k, s) for k, s, _ in entries], truth)
    assert ref["worst"] == worst and ref["n"] == n_, (cfg, n, ref, worst)
    return mk(n_, worst, [("", k, s, None) for k, s, _ in display], "" if entries else "no output")


def novelist_all(truths):
    use = ["mappingFeatureId", "structurePerIdRank", "InChIkey2D", "CSI:FingerIDScore", "name", "smiles"]
    pat = re.compile(r"_(challenge_\d+)_challenge")
    frames = []
    for frame in pd.read_csv(SUMM / "denovo_structure_identifications_all.tsv", sep="\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE,
                             usecols=use, chunksize=100000):
        frame["challenge"] = frame.mappingFeatureId.str.extract(pat, expand=False)
        frames.append(frame)
    rows = pd.concat(frames, ignore_index=True)
    rows["rank"] = pd.to_numeric(rows.structurePerIdRank)
    rows = rows.sort_values("rank", kind="stable").drop_duplicates(["challenge", "InChIkey2D"], keep="first")
    out = {}
    for cid, df in rows.groupby("challenge"):
        out[cid] = [(str(k), num(s), (nm, sm)) for k, s, nm, sm in zip(df.InChIkey2D, df["CSI:FingerIDScore"], df["name"], df["smiles"]) if num(s) is not None]
    return out


def novelist_result(entries, truth):
    if not entries:
        return mk(0, None, [], "no output")
    n_, worst, display = rank_scored(entries, truth)
    ref = S.scored_rank([(k, s) for k, s, _ in entries], truth)
    assert ref["worst"] == worst and ref["n"] == n_
    return mk(n_, worst, [((p[0] or p[1]), k, s, None) for k, s, p in display])


def sirius_formula_all():
    f = read_tsv(SUMM / "formula_identifications_all.tsv")
    f["cid"] = f.mappingFeatureId.str.extract(r"_(challenge_\d+)_challenge", expand=False)
    f["_r"] = pd.to_numeric(f["formulaRank"])
    f["_s"] = pd.to_numeric(f["SiriusScore"])
    f = f.sort_values("_r", kind="stable")
    out = {}
    for cid, df in f.groupby("cid"):
        df = df.drop_duplicates("molecularFormula", keep="first")      # collapse adducts, keep the best row
        out[cid] = [(a, s, adduct) for a, s, adduct in zip(df["molecularFormula"], df["_s"], df["adduct"])]
    return out


def sirius_formula_result(entries, truth_formula):
    if not entries:
        return mk(0, None, [], "no output")
    n_, worst, display = rank_scored(entries, truth_formula)
    status = ""
    target = next((s for f, s, _ in entries if f == truth_formula), None)
    if target is not None:
        tied = sum(s == target for _, s, _ in entries)
        better = sum(s > target for _, s, _ in entries)
        if tied > 1:
            status = f"tie: {tied - 1} other formula(s) share the truth score; best rank {better + 1}"
    return mk(n_, worst, [(f, a, s, None) for f, s, a in display], status)


# ------------------------------------------------------------------ CFM-ID team's own CASMI 2022 submission (cfmid_casmi2022_entry/)
def cfmid_entry_all(official):
    """{challenge id: (submitted formula, submitted adduct, submitted IK14)} from the participant's two submission sheets (Compound Number = challenge id)."""
    frames = [pd.read_excel(CFMID / f, dtype=str, keep_default_na=False) for f in ("Priority_1-250.xlsx", "Bonus_251-500.xlsx")]
    df = pd.concat(frames, ignore_index=True)
    df["id"] = df["Compound Number"].astype(int)
    assert df["id"].tolist() == IDS, "CFM-ID entry: expected one row for each of the 500 challenges, in order"
    assert (df["Priority/Bonus"].values == official.loc[df["id"].astype(str), "priority"].values).all()
    out = {}
    for n, f, ad, ik in zip(df["id"], df["Molecular Formula"], df["Adduct"], df["InChIKey or SMILES"]):
        ik = ik.strip()
        assert re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", ik), (n, ik)       # every answer is an InChIKey (no SMILES)
        out[n] = (f.strip(), ad.strip(), ik[:14])
    return out


def cfmid_entry_result(entry, truth):
    """One answer, n = 1. Truth rank 1 on a hit, blank (truth not in the returned list) otherwise; the score column stays blank (no score submitted)."""
    return mk(1, 1 if entry[2] == truth else None, [("", entry[2], None, None)])


def cfmid_entry_formula_result(entry, truth_formula):
    return mk(1, 1 if counts_key(entry[0]) == counts_key(truth_formula) else None, [(entry[0], entry[1], None, None)])


# ------------------------------------------------------------------ names of PubChem candidates (SIRIUS structure table)
def pubchem_names(needed):
    cache = SCRATCH / "pubchem_names.pkl"
    if cache.exists():
        got = pickle.loads(cache.read_bytes())
        if needed <= set(got) | got.get("__missing__", set()):
            return got
    key_set = set(needed)
    out = {}
    pat = re.compile(r"_(challenge_\d+)_challenge")
    for frame in pd.read_csv(SUMM / "structure_identifications_all.tsv", sep="\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE,
                             usecols=["mappingFeatureId", "InChIkey2D", "name"], chunksize=200000):
        cid = frame.mappingFeatureId.str.extract(pat, expand=False)
        keys = cid + "|" + frame.InChIkey2D
        hit = frame[keys.isin(key_set)]
        for c, k, nm in zip(cid[hit.index], hit.InChIkey2D, hit["name"]):
            out.setdefault(f"{c}|{k}", nm)
    SCRATCH.mkdir(parents=True, exist_ok=True)
    out["__missing__"] = key_set - set(out)
    cache.write_bytes(pickle.dumps(out))
    return out


# ------------------------------------------------------------------ workbook
GREEN = PatternFill("solid", fgColor="C6EFCE")
GREEN_DARK = PatternFill("solid", fgColor="63BE7B")
HEAD = PatternFill("solid", fgColor="D9E1F2")
GREY = PatternFill("solid", fgColor="EDEDED")
BOLD = Font(bold=True)
STATUS_FONT = Font(italic=True, color="7F7F7F")


def write_tab(wb, tab, kind, data, cases):
    ws = wb.create_sheet(tab)
    third = "score"
    if kind == "formula":
        head = ["id", "polarity", "truth rank (worst)", "n candidates", "status"]
        for r in range(1, 11):
            head += [f"#{r} formula", f"#{r} adduct", f"#{r} score"]
    else:
        head = ["id", "polarity", "truth rank (worst)", "n candidates", "status"]
        for r in range(1, 11):
            head += [f"#{r} name", f"#{r} IK14", f"#{r} {third}"]
    ws.append(head)
    for c in ws[1]:
        c.font, c.fill, c.alignment = BOLD, HEAD, Alignment(wrap_text=True, vertical="top")
    for case in cases:
        n = case["id"]
        res = data[n]
        row = [n, case["polarity"], res["rank"], res["n"], res["status"] or ("truth not in list" if res["rank"] is None else "")]
        truth_key = case["truth_formula"] if kind == "formula" else case["truth_ik14"]
        hit_at = None
        for i, (nm, key, sc, _x) in enumerate(res["top"]):
            row += [nm, key, sc if sc is None else float(sc)]
            if (kind == "formula" and counts_key(nm) == counts_key(truth_key)) or (kind != "formula" and key == truth_key):
                hit_at = i
        ws.append(row)
        r = ws.max_row
        if res["status"]:
            ws.cell(r, 5).font = STATUS_FONT
        if hit_at is not None:
            for j in range(3):
                ws.cell(r, 6 + 3 * hit_at + j).fill = GREEN
        if res["rank"] is not None and res["rank"] <= 10:
            ws.cell(r, 3).fill = GREEN_DARK
            ws.cell(r, 3).font = BOLD
        for i in range(len(res["top"])):
            ws.cell(r, 8 + 3 * i).number_format = "0.0000"
    ws.freeze_panes = "B2"
    ws.auto_filter.ref = f"A1:{get_column_letter(len(head))}{ws.max_row}"
    ws.column_dimensions["A"].width = 6
    ws.column_dimensions["B"].width = 10
    ws.column_dimensions["C"].width = 11
    ws.column_dimensions["D"].width = 11
    ws.column_dimensions["E"].width = 34
    for i in range(10):
        a = 6 + 3 * i
        ws.column_dimensions[get_column_letter(a)].width = 14 if kind == "formula" else 30
        ws.column_dimensions[get_column_letter(a + 1)].width = 12 if kind == "formula" else 16
        ws.column_dimensions[get_column_letter(a + 2)].width = 9
    ws.row_dimensions[1].height = 32


def stats_of(rows):
    """rows: [(rank, status)] of the challenges the situation applies to. As write_report.stats."""
    n = len(rows)
    k = lambda t: sum(r is not None and r <= t for r, _ in rows)
    return {"n": n, "wrong": sum(s.startswith("not run: true adduct") for _, s in rows), "no_output": sum(s.startswith("no output") for _, s in rows),
            "present": sum(r is not None for r, _ in rows), "t1": k(1), "t3": k(3), "t10": k(10)}


def tab_rows(data, cases, subset=None, polarity=None):
    out = []
    for case in cases:
        if subset is not None and case["id"] not in subset:
            continue
        if polarity and case["polarity"] != polarity:
            continue
        res = data[case["id"]]
        if res["status"] == NOT_SUPPORTED:
            continue
        out.append((res["rank"], res["status"]))
    return out


def pct(k, n):
    return k / n if n else None


def write_summary(wb, results, cases, unseen, notes_extra):
    ws = wb.active
    ws.title = "Summary"
    ws["A1"] = "CASMI 2022, 500 challenges: mzDesk and comparators"
    ws["A1"].font = Font(bold=True, size=14)
    ws["A2"] = "Top-k = truth IK14 within the first k candidates (ties: worst rank). Challenges without output stay in the denominator."
    head = ["Group", "Situation (tab)", "Level", "n scored", "no output", "not run (true adduct not default)", "truth listed",
            "Top-1", "Top-1 %", "Top-3", "Top-3 %", "Top-10", "Top-10 %",
            "positive n", "positive Top-1", "positive Top-1 %", "negative n", "negative Top-1", "negative Top-1 %",
            "leak-free n", "leak-free Top-1", "leak-free Top-1 %", "leak-free Top-3", "leak-free Top-3 %", "leak-free Top-10", "leak-free Top-10 %"]
    ws.append([])
    ws.append(head)
    hr = ws.max_row
    for c in ws[hr]:
        c.font, c.fill, c.alignment = BOLD, HEAD, Alignment(wrap_text=True, vertical="top")
    for tab, group, kind, key in TABS:
        data = results[tab]
        a, p, ng, lf = (stats_of(tab_rows(data, cases)), stats_of(tab_rows(data, cases, polarity="positive")),
                        stats_of(tab_rows(data, cases, polarity="negative")), stats_of(tab_rows(data, cases, subset=unseen)))
        row = [group, tab, "formula" if kind == "formula" else "structure", a["n"], a["no_output"], a["wrong"], a["present"],
               a["t1"], pct(a["t1"], a["n"]), a["t3"], pct(a["t3"], a["n"]), a["t10"], pct(a["t10"], a["n"]),
               p["n"], p["t1"], pct(p["t1"], p["n"]), ng["n"] or None, ng["t1"] if ng["n"] else None, pct(ng["t1"], ng["n"]),
               lf["n"], lf["t1"], pct(lf["t1"], lf["n"]), lf["t3"], pct(lf["t3"], lf["n"]), lf["t10"], pct(lf["t10"], lf["n"])]
        if key in TOP1_ONLY:                                   # single answer: Top-3 / Top-10 / truth listed are not applicable
            for col in (7, 10, 11, 12, 13, 23, 24, 25, 26):
                row[col - 1] = NA
        ws.append(row)
        r = ws.max_row
        for col in (9, 11, 13, 16, 19, 22, 24, 26):
            ws.cell(r, col).number_format = "0.0%"
        if key in TOP1_ONLY:
            for col in (7, 10, 11, 12, 13, 23, 24, 25, 26):
                ws.cell(r, col).font = STATUS_FONT
                ws.cell(r, col).alignment = Alignment(horizontal="right")
        if kind == "formula":
            for col in range(1, len(head) + 1):
                ws.cell(r, col).fill = GREY
    last = ws.max_row
    ws.freeze_panes = ws.cell(hr + 1, 3)
    ws.auto_filter.ref = f"A{hr}:{get_column_letter(len(head))}{last}"
    widths = [30, 26, 10, 9, 9, 14, 9] + [8, 8] * 3 + [10, 10, 10, 10, 10, 10] + [10, 10, 10, 10, 10, 10, 10]
    for i, w in enumerate(widths, 1):
        ws.column_dimensions[get_column_letter(i)].width = w
    ws.row_dimensions[hr].height = 48
    ws.append([])
    ws.append(["Notes"])
    ws.cell(ws.max_row, 1).font = BOLD
    for line in notes_extra:
        ws.append([line])
    return ws


def write_truth(wb, cases, truth_df, lib, unseen, names, precursor):
    ws = wb.create_sheet("Truth")
    head = ["id", "priority", "polarity", "name", "formula", "adduct", "SMILES", "InChIKey", "IK14", "measured precursor m/z",
            "truth structure has a MassBank/GNPS spectrum (in library)", "in the leak-free subset (202)"]
    ws.append(head)
    for c in ws[1]:
        c.font, c.fill, c.alignment = BOLD, HEAD, Alignment(wrap_text=True, vertical="top")
    for case in cases:
        n = case["id"]
        t = truth_df.loc[str(n)]
        ws.append([n, case["priority"], case["polarity"], names.get(n, ""), case["truth_formula"], case["truth_adduct"], t["smiles_official"], t["inchikey"],
                   case["truth_ik14"], precursor[n], "yes" if n in lib else "no", "yes" if n in unseen else "no"])
        ws.cell(ws.max_row, 10).number_format = "0.00000"
    ws.freeze_panes = "B2"
    ws.auto_filter.ref = f"A1:{get_column_letter(len(head))}{ws.max_row}"
    for i, w in enumerate([6, 9, 10, 36, 14, 11, 50, 29, 16, 14, 18, 14], 1):
        ws.column_dimensions[get_column_letter(i)].width = w
    ws.row_dimensions[1].height = 62


# ------------------------------------------------------------------ main
def build(args):
    scores = json.loads((HERE / "scores.json").read_text())
    cases = scores["cases"]
    assert [c["id"] for c in cases] == IDS
    unseen = set(scores["meta"]["subsets"]["unseen_common"])
    lib = set(json.loads((HERE / "leak_lists/lib_ids.json").read_text())["lib_any"])
    assert len(unseen) == 202
    truth_df, _old = load_truth()
    truth_df = truth_df.assign(inchikey=truth_df.inchikey)

    results = {tab: {} for tab, *_ in TABS}
    novel = novelist_all(IDS)
    sform = sirius_formula_all()
    cfmid = cfmid_entry_all(truth_df)
    for case in cases:
        n, truth, tf, pol, ta = case["id"], case["truth_ik14"], case["truth_formula"], case["polarity"], case["truth_adduct"]
        cid = f"challenge_{n:03d}"
        for tab, _g, kind, key in TABS:
            if key in S.CONFIG_GROUPS:
                r = mzdesk_result(key, n, truth, tf)
                if key in S.POSITIVE_ONLY and pol != "positive":
                    r = mk(0, None, [], NOT_SUPPORTED)
            elif key in S.DEFAULT_CONFIGS:
                r = standalone_result(key, n, truth, tf, pol, ta)
            elif key == "MSNovelist":
                r = novelist_result(novel.get(cid, []), truth)
            elif key == "sirius_formula":
                r = sirius_formula_result(sform.get(cid, []), tf)
            elif key == "cfmid_entry":
                r = cfmid_entry_result(cfmid[n], truth)
            elif key == "cfmid_entry_formula":
                r = cfmid_entry_formula_result(cfmid[n], tf)
            else:
                r = jm_result(key, n, truth, pol)
            results[tab][n] = r
    # names for the PubChem-based comparators (and for the truth)
    needed = set()
    for case in cases:
        n = case["id"]
        needed.add(f"challenge_{n:03d}|{case['truth_ik14']}")
        for tab in ("CSI FingerID", "JESTR", "MIST-contrastive", "CFM-ID (CASMI entry)"):
            for _nm, key, _s, _x in results[tab][n]["top"]:
                needed.add(f"challenge_{n:03d}|{key}")
    pc = pubchem_names(needed)
    for tab in ("CSI FingerID", "JESTR", "MIST-contrastive", "CFM-ID (CASMI entry)"):
        for n in IDS:
            smiles = None
            top = []
            for nm, key, sc, x in results[tab][n]["top"]:
                name = pc.get(f"challenge_{n:03d}|{key}")
                if not name:
                    smiles = smiles if smiles is not None else smiles_map(n)
                    name = smiles.get(key, "")
                top.append((name, key, sc, x))
            results[tab][n]["top"] = top
    # truth names: PubChem name -> a name in any mzDesk / default candidate list -> blank
    names, src = {}, {"pubchem": 0, "mzdesk": 0, "none": 0}
    for case in cases:
        n, truth = case["id"], case["truth_ik14"]
        nm = pc.get(f"challenge_{n:03d}|{truth}")
        if nm:
            names[n] = nm
            src["pubchem"] += 1
            continue
        found = ""
        for tab, _g, kind, key in TABS:
            if kind != "struct" or not (key in S.CONFIG_GROUPS or key in S.DEFAULT_CONFIGS):
                continue
            for nm2, k2, _s, _x in results[tab][n]["top"]:
                if k2 == truth and nm2:
                    found = re.sub(r"\s*\[M[^\]]*\][+-]?$", "", nm2)
                    break
            if found:
                break
        names[n] = found
        src["mzdesk" if found else "none"] += 1

    precursor = {n: parentmass(n) for n in IDS}
    prior = read_tsv(SIRIUS / "score/per_challenge.tsv")
    pr = pd.to_numeric(prior["formula_rank"], errors="coerce")
    sirius_opt = [int((pr <= k).sum()) for k in (1, 3, 10)]
    wb = Workbook()
    n_pos = sum(c["polarity"] == "positive" for c in cases)
    n285 = len(scores["meta"]["subsets"]["default_adduct"])
    notes = [
        "Scoring rule: a hit is an InChIKey first block (IK14) equal to the official key; ties use the worst rank; a challenge without output stays in the denominator. "
        "Formula tabs (BUDDY, SIRIUS formula) rank neutral molecular formulas and are shaded grey in this table.",
        f"Truth rank columns give the real rank even beyond 10 (blank = truth not in the candidate list). The top-10 groups show the highlighted truth candidate (light green); the rank cell is dark green when the truth is within the top 10. "
        "Where candidates tie with the truth, the truth is placed last of its ties, so its position equals its (worst) rank.",
        "mzDesk main: 'full' is the workflow as run (library searched; library entries of the truth structure were not removed). 'leak-free' does no library search "
        "(step 3 empty, gate not passed) and annotates from scratch (MS-FINDER, FRIGID when needed).",
        "mzDesk ablation: every layer listed is merged regardless of gate; step 7 (isomers) is omitted, which does not change ranks. frigid_only is positive ions only. buddy_only ranks the truth formula in the step-2 formula list.",
        f"Single software, default settings: adduct fixed to [M+H]+ (positive) / [M-H]- (negative). Challenges whose true adduct is another one are scored wrong without a run "
        f"('not run: true adduct not default', kept in the denominator; {500 - n285} challenges), so BUDDY / MS-FINDER / FRIGID run on {n285} challenges. "
        "The library search needs no adduct and runs on all 500. FRIGID default uses the BUDDY top formula at [M+H]+. Their ranking is each tool's own score.",
        f"FRIGID (default and frigid_only) and MIST-contrastive support positive ions only: n scored = {n_pos} positive challenges; negative challenges are marked '{NOT_SUPPORTED}' and excluded from the denominator.",
        "BUDDY and SIRIUS formula tabs rank neutral formulas (distinct formulas; BUDDY in its list order, SIRIUS by SiriusScore, adducts collapsed to the best row). "
        f"SIRIUS formula: ties at the truth use the worst rank, as in all other tabs (the status column gives the best rank where tied). The existing SIRIUS table (sirius_run/results/score/per_challenge.tsv) "
        f"uses the best (optimistic) rank: Top-1/3/10 = {sirius_opt[0]}/{sirius_opt[1]}/{sirius_opt[2]} of 500 there, versus the worst-rank counts in this table, because adduct variants of the same fragmentation tree share one SiriusScore.",
        "Comparators: CSI:FingerID ranks the PubChem candidate list by CSI:FingerID score (jm_run/inputs/candidates); MSNovelist is the SIRIUS de novo output ranked by CSI:FingerID score; "
        "JESTR and MIST are scored on the same candidate lists (jm_run/jestr/scores, jm_run/mist/scores). MSNovelist and SIRIUS (CSI, MSNovelist) depend on online services.",
        "CFM-ID (CASMI entry) is not a run of this benchmark: it is the CFM-ID team's own official CASMI 2022 submission (Fei Wang et al., CFM-ID 4; cfmid_casmi2022_entry/), end to end with its own candidate databases "
        "(12 merged compound databases), its own adduct and formula decisions, and one answer per challenge. Only Top-1 exists (Top-3, Top-10 and 'truth listed' are n/a, not 0); the structure tab holds the submitted IK14 "
        "(score blank), the formula tab the submitted neutral formula and adduct. It is not directly comparable to the comparators scored on the same PubChem candidate lists.",
        "Comparator names come from the SIRIUS PubChem table (by IK14); where none exists the name column shows the SMILES. IK14 is the first block of the InChIKey.",
        "Leak-free subset: the 202 challenges whose truth structure is in none of the tools' training / reference sets (SIRIUS fingerprint predictor, FRIGID MSG/4M, MIST, JESTR, local MassBank/GNPS); "
        "MSNovelist and the FRIGID language-model pre-training cannot be excluded. Truth tab 'in library' = the truth IK14 has any MassBank/GNPS spectrum in the local index (leak_lists/lib_ids.json, 'lib_any').",
        f"Truth names: the official key has no compound names; names are taken from the SIRIUS PubChem table ({src['pubchem']}), otherwise from an mzDesk / default candidate list that contains the truth ({src['mzdesk']}); "
        f"{src['none']} are blank.",
    ]
    write_summary(wb, results, cases, unseen, notes)
    write_truth(wb, cases, truth_df, lib, unseen, names, precursor)
    for tab, _g, kind, key in TABS:
        write_tab(wb, tab, kind, results[tab], cases)
    wb.properties.title = "CASMI 2022 mzDesk benchmark"
    return wb, results, scores


# ------------------------------------------------------------------ self-check
def self_check(target, scores):
    wb = load_workbook(target, read_only=False, data_only=True)
    expected = ["Summary", "Truth"] + [t[0] for t in TABS]
    assert wb.sheetnames == expected, wb.sheetnames
    cases = scores["cases"]
    d = W.Data(scores)
    report_path = ROOT / "benchmark/result/CASMI2022_mzdesk_full500.md"
    report = report_path.read_text(encoding="utf-8") if report_path.exists() else None      # the report is checked only when it exists
    official = load_truth()[0]
    cfmid = cfmid_entry_all(official)
    sub_unseen = d.subsets["unseen_common"]
    problems, lines = [], []
    sheet_stats = {}
    for tab, _g, kind, key in TABS:
        ws = wb[tab]
        rows = list(ws.iter_rows(min_row=2, values_only=True))
        assert len(rows) == 500 and [r[0] for r in rows] == IDS, tab
        bypol = {}
        for r in rows:
            bypol[r[0]] = (r[1], r[2], r[3], r[4] or "")
        def sel(subset=None, pol=None):
            return [(v[1], v[3]) for i, v in bypol.items() if v[3] != NOT_SUPPORTED and (subset is None or i in subset) and (pol is None or v[0] == pol)]
        sheet_stats[tab] = {k: stats_of(sel(*a)) for k, a in (("all", ()), ("pos", (None, "positive")), ("neg", (None, "negative")), ("lf", (sub_unseen,)))}
        # --- rank by rank against scores.json
        bad = 0
        for c in cases:
            rank = bypol[c["id"]][1]
            if key in S.CONFIG_GROUPS or key in S.DEFAULT_CONFIGS:
                ref = c["configs"].get(key)
                na = ref is None or (key in POSITIVE_ONLY and c["polarity"] != "positive")
                exp = None if na else ref["worst"]
            elif key in c["comparators"]:
                ref = c["comparators"][key]
                na = ref is None
                exp = None if na else ref["worst"]
            else:
                exp, na = None, False        # SIRIUS formula: compared with per_challenge.tsv below
            if key == "sirius_formula":
                continue
            if key in TOP1_ONLY:         # not in scores.json: recomputed from the submission files (rank 1 on a hit, blank otherwise)
                exp = 1 if ((cfmid[c["id"]][2] == c["truth_ik14"]) if kind == "struct"
                            else (counts_key(cfmid[c["id"]][0]) == counts_key(c["truth_formula"]))) else None
                if rank != exp:
                    bad += 1
                    problems.append(f"{tab} #{c['id']}: workbook rank {rank} vs submission {exp}")
                continue
            if na != (bypol[c["id"]][3] == NOT_SUPPORTED) or (not na and rank != exp):
                bad += 1
                problems.append(f"{tab} #{c['id']}: workbook rank {rank} vs scores.json {exp}")
        # --- aggregate counts against scores.json via the report's own functions
        if key not in ("sirius_formula",) + tuple(TOP1_ONLY):
            for label, sk, subset, pol in (("all", "all", "all", None), ("pos", "pos", "all", "positive"), ("neg", "neg", "all", "negative"), ("lf", "lf", "unseen_common", None)):
                ref = W.stats(d.rows(key, subset, pol))
                mine = sheet_stats[tab][sk]
                for a, b in (("n", "n"), ("no_output", "no_output"), ("present", "present"), ("t1", "t1"), ("t3", "t3"), ("t10", "t10"), ("wrong", "wrong")):
                    if ref[a] != mine[b]:
                        problems.append(f"{tab} [{label}] {a}: workbook {mine[b]} vs scores.json {ref[a]}")
            # --- against the report text (all challenges)
            m = sheet_stats[tab]["all"]
            if key in ("full", "leakfree_full") or key in S.COMPARATORS or key in ("all_layers", "library_only", "msfinder_only", "frigid_only", "no_library", "no_msfinder", "no_frigid", "buddy_only"):
                pattern = f"| {m['n']} | {m['no_output']} | {m['present']} | {W.cell(m['t1'], m['n'])} | {W.cell(m['t3'], m['n'])} | {W.cell(m['t10'], m['n'])} |"
                in_report = report is not None and pattern in report
            else:
                pattern = (f"| {m['n']} | {m['no_output']} | {m['wrong']} | {m['present']} | {W.cell(m['t1'], m['n'])} | {W.cell(m['t3'], m['n'])} | {W.cell(m['t10'], m['n'])} |")
                in_report = report is not None and pattern in report
            if report is not None and not in_report:
                problems.append(f"{tab}: the report has no row '{pattern}'")
            lines.append(f"{tab:26} n={m['n']:3} noout={m['no_output']:3} wrong={m['wrong']:3} listed={m['present']:3} "
                         f"T1={m['t1']:3} T3={m['t3']:3} T10={m['t10']:3}  rank mismatches vs scores.json: {bad}  in report text: {in_report if report is not None else 'report file absent'}")
        elif key in TOP1_ONLY:
            m = sheet_stats[tab]["all"]
            lines.append(f"{tab:26} n={m['n']:3} noout={m['no_output']:3} wrong={m['wrong']:3} listed={m['present']:3} "
                         f"T1={m['t1']:3} (Top-1 only)  rank mismatches vs the submission files: {bad}")
    # SIRIUS formula against the existing SIRIUS score table (best rank + ties = worst rank)
    prior = read_tsv(SIRIUS / "score/per_challenge.tsv").set_index("challenge_id")
    ws = wb["SIRIUS formula"]
    bad = 0
    for r in ws.iter_rows(min_row=2, values_only=True):
        row = prior.loc[str(r[0])]
        best = num(row["formula_rank"])
        exp = None if best is None else int(best) + int(row["formula_n_tied"] or 0)
        if r[2] != exp:
            bad += 1
            problems.append(f"SIRIUS formula #{r[0]}: workbook {r[2]} vs best+ties {exp}")
    pc = pd.to_numeric(prior["formula_rank"], errors="coerce")
    lines.append(f"{'SIRIUS formula':26} rank mismatches vs per_challenge.tsv (formula_rank + n_tied): {bad}; "
                 f"existing optimistic table: Top-1 {int((pc <= 1).sum())}, Top-3 {int((pc <= 3).sum())}, Top-10 {int((pc <= 10).sum())}; "
                 f"this workbook (worst rank on ties): Top-1 {sheet_stats['SIRIUS formula']['all']['t1']}, Top-3 {sheet_stats['SIRIUS formula']['all']['t3']}, Top-10 {sheet_stats['SIRIUS formula']['all']['t10']}")
    # the Summary sheet itself against the same numbers
    sm = wb["Summary"]
    head_row = next(i for i, r in enumerate(sm.iter_rows(values_only=True), 1) if r[0] == "Group")
    for r in sm.iter_rows(min_row=head_row + 1, max_row=head_row + len(TABS), values_only=True):
        s = sheet_stats[r[1]]
        top1_only = next(t[3] for t in TABS if t[0] == r[1]) in TOP1_ONLY
        mine = (s["all"]["n"], s["all"]["no_output"], s["all"]["wrong"], s["all"]["present"], s["all"]["t1"], s["all"]["t3"], s["all"]["t10"],
                s["pos"]["t1"], s["lf"]["n"], s["lf"]["t1"], s["lf"]["t3"], s["lf"]["t10"])
        if top1_only:
            mine = (mine[0], mine[1], mine[2], NA, mine[4], NA, NA, mine[7], mine[8], mine[9], NA, NA)
        got = (r[3], r[4], r[5], r[6], r[7], r[9], r[11], r[14], r[19], r[20], r[22], r[24])
        if mine != got:
            problems.append(f"Summary row {r[1]}: {got} vs tab {mine}")
    return problems, lines


def spot_check(target, scores):
    """3 challenges per group, compared by hand-written code against the source files (not through the functions above)."""
    wb = load_workbook(target, data_only=True)
    out = []

    def row_of(tab, n):
        return list(wb[tab].iter_rows(min_row=n + 1, max_row=n + 1, values_only=True))[0]

    def check(label, ok, detail):
        out.append(f"  {'OK ' if ok else 'BAD'} {label}: {detail}")
    cases = {c["id"]: c for c in scores["cases"]}
    # mzDesk full / leak-free / ablation: first candidate and the truth position straight from replay json
    for tab, cfg, n in (("mzDesk full", "full", 9), ("mzDesk leak-free", "leakfree_full", 116), ("Abl no_frigid", "no_frigid", 227)):
        res = json.loads((HERE / "replay" / cfg / f"challenge_{n:03d}.json").read_text())["result"]
        cs = [c for c in res["candidates"] if c["layer"] != "formula_only"]
        r = row_of(tab, n)
        truth = cases[n]["truth_ik14"]
        pos = next((i + 1 for i, c in enumerate(cs) if c["ik14"] == truth), None)
        ok = r[3] == len(cs) and r[6] == (cs[0]["ik14"] if cs else None) and r[2] == (pos if pos is None or True else None)
        check(f"{tab} #{n}", ok, f"top1 {r[6]} score {r[7]}; n {r[3]} (source {len(cs)}); truth listed pos {pos}, workbook rank {r[2]}")
    # default settings: library and msfinder
    for tab, tool, n in (("Default library", "library", 9), ("Default MS-FINDER", "msfinder", 9), ("Default BUDDY (formula)", "buddy", 9)):
        res = json.loads((HERE / "standalone" / tool / f"challenge_{n:03d}.json").read_text())
        r = row_of(tab, n)
        if tool == "library":
            first = max(res["compounds"], key=lambda c: c["score"])
            check(f"{tab} #{n}", r[6] == first["ik14"] and abs(r[7] - first["score"]) < 1e-9, f"top1 {r[6]} vs {first['ik14']}; n {r[3]} vs {len(res['compounds'])}")
        elif tool == "msfinder":
            first = max(res["candidates"], key=lambda c: c["score"])
            check(f"{tab} #{n}", r[6] == first["ik14"], f"top1 {r[6]} vs {first['ik14']}; rank {r[2]}")
        else:
            check(f"{tab} #{n}", r[5] == res["formulas"][0]["formula"], f"top1 {r[5]} vs {res['formulas'][0]['formula']}; rank {r[2]}")
    # comparators from the raw files
    n = 3
    for tab, path, kc, sc in (("CSI FingerID", JM / "inputs/candidates/challenge_003.tsv.gz", "ikey14", "csi_score"),
                              ("JESTR", JM / "jestr/scores/challenge_003.tsv.gz", "ikey14", "score"),
                              ("MIST-contrastive", JM / "mist/scores/challenge_003.tsv.gz", "ikey14", "score_contrast")):
        rows = [r for r in csv.DictReader(gzip.open(path, "rt"), delimiter="\t")]
        best = max(rows, key=lambda r: float(r[sc]))
        truth = cases[n]["truth_ik14"]
        t = next((float(r[sc]) for r in rows if r[kc] == truth), None)
        exp_rank = None if t is None else sum(float(r[sc]) >= t for r in rows)
        r = row_of(tab, n)
        check(f"{tab} #{n}", r[6] == best[kc] and r[2] == exp_rank, f"top1 {r[6]} vs {best[kc]}; rank {r[2]} vs {exp_rank}; n {r[3]} vs {len(rows)}")
    # SIRIUS formula, MSNovelist, BUDDY ablation
    rows = [r for r in csv.reader(open(SUMM / "formula_identifications_all.tsv"), delimiter="\t", quoting=csv.QUOTE_NONE)]
    h = rows[0]
    fi = {k: h.index(k) for k in ("formulaRank", "molecularFormula", "SiriusScore", "mappingFeatureId")}
    sub = [r for r in rows[1:] if "challenge_007_challenge" in r[fi["mappingFeatureId"]]]
    best = max(sub, key=lambda r: float(r[fi["SiriusScore"]]))
    r = row_of("SIRIUS formula", 7)
    check("SIRIUS formula #7", abs(r[7] - float(best[fi["SiriusScore"]])) < 1e-9 and any(x[fi["molecularFormula"]] == r[5] for x in sub),
          f"top score {r[7]} vs {best[fi['SiriusScore']]} (4-way tie, truth placed last); rank {r[2]}")
    buddy = json.loads((HERE / "replay/buddy_only/challenge_116.json").read_text())["result"]["formulas"]
    r = row_of("Abl buddy_only (formula)", 116)
    check("buddy_only #116", r[5] == buddy[0]["formula"], f"top1 {r[5]} vs {buddy[0]['formula']}; rank {r[2]}")
    r = row_of("MSNovelist", 9)
    check("MSNovelist #9", r[3] == cases[9]["comparators"]["MSNovelist"]["n"], f"n {r[3]} vs scores.json {cases[9]['comparators']['MSNovelist']['n']}")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", default=str(TARGET))
    ap.add_argument("--replace", action="store_true", help="replace an existing workbook")
    ap.add_argument("--no-check", action="store_true")
    args = ap.parse_args()
    target = Path(args.target)
    if target.exists() and not args.replace:
        raise FileExistsError(target)
    wb, results, scores = build(args)
    wb.save(target)
    print("written", target)
    for ws in wb.worksheets:
        print(f"  {ws.title:28} {ws.max_row - 1} data rows" if ws.title != "Summary" else f"  {ws.title:28} {len(TABS)} situations")
    if not args.no_check:
        problems, lines = self_check(target, scores)
        print("\nSelf-check (workbook re-read from disk):")
        print("\n".join(lines))
        print(f"{len(problems)} problems" + "".join("\n  " + p for p in problems[:60]))
        print("\nSpot checks against source files:")
        print("\n".join(spot_check(target, scores)))


if __name__ == "__main__":
    main()
