#!/usr/bin/env python3
"""Score collect.py / replay.py outputs against the official CASMI 2022 key and the comparators -> scores.json.

  python3 benchmark/casmi2022/mzdesk_ablation_20261030/score.py [--dir DIR] [--partial] [--replace]

Judging rules: a hit is an InChIKey first block (IK14) equal to the key's; ties use the worst
rank; a challenge without output stays in the denominator; MIST positive only. Comparator ranks are recomputed from their score files and checked
against the existing score tables. The truth is read here only (never in collect.py or replay.py).
No engine runs. Needs pandas and RDKit; does not import mzdesk.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
from rdkit import Chem

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CASMI = ROOT / "benchmark/casmi2022"
JM = CASMI / "jm_run"
SIRIUS = CASMI / "sirius_run/results"
LEAK = HERE / "leak_lists"
COMPARATORS = ("CSI:FingerID", "MSNovelist", "JESTR", "MIST-FP", "MIST-contrastive")

# configuration -> group in the report; frigid_only is scored on positive challenges only (the engine skips negative adducts)
CONFIG_GROUPS = {
    "full": "workflow", "leakfree_full": "workflow",
    "all_layers": "ablation", "library_only": "ablation", "msfinder_only": "ablation",
    "frigid_only": "ablation", "no_library": "ablation", "no_msfinder": "ablation", "no_frigid": "ablation", "buddy_only": "formula",
}
POSITIVE_ONLY = ("frigid_only", "default_frigid")
# single tools with their default settings (standalone.py): library needs no adduct; BUDDY / MS-FINDER / FRIGID use [M+H]+ / [M-H]- and a
# challenge whose official adduct is another one is scored wrong without a run
DEFAULT_CONFIGS = ("default_library", "default_buddy", "default_msfinder", "default_frigid")
DEFAULT_ADDUCT = {"positive": "[M+H]+", "negative": "[M-H]-"}


def number(value):
    try:
        n = float(value)
        return n if math.isfinite(n) else None
    except (TypeError, ValueError):
        return None


# ------------------------------------------------------------------ comparators (as full500/score.py)
def load_scores(path, key_col, score_col):
    if not path.exists():
        return []
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE)
    return [(str(k), number(s)) for k, s in zip(df[key_col], df[score_col]) if number(s) is not None]


def scored_rank(pairs, truth):
    assert len({k for k, _ in pairs}) == len(pairs), "Duplicate structures in scored candidates"
    if not pairs:
        return {"n": 0, "top": None, "best": None, "worst": None, "ties": 0}
    top = max(pairs, key=lambda p: p[1])[0]
    target = next((s for k, s in pairs if k == truth), None)
    if target is None:
        return {"n": len(pairs), "top": top, "best": None, "worst": None, "ties": 0}
    better = sum(s > target for _, s in pairs)
    tied = sum(s == target for _, s in pairs)
    return {"n": len(pairs), "top": top, "best": better + 1, "worst": better + tied, "ties": tied - 1}


def novelist_pairs(ids):
    p = SIRIUS / "sirius_summaries/denovo_structure_identifications_all.tsv"
    use = ["mappingFeatureId", "structurePerIdRank", "InChIkey2D", "CSI:FingerIDScore"]
    all_rows = []
    pattern = re.compile(r"_(challenge_\d+)_challenge")
    selected = {f"challenge_{n:03d}" for n in ids}
    for frame in pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE, usecols=use, chunksize=100000):
        frame["challenge"] = frame.mappingFeatureId.str.extract(pattern, expand=False)
        all_rows.append(frame[frame.challenge.isin(selected)])
    rows = pd.concat(all_rows, ignore_index=True)
    rows["rank"] = pd.to_numeric(rows.structurePerIdRank)
    rows = rows.sort_values("rank", kind="stable").drop_duplicates(["challenge", "InChIkey2D"], keep="first")
    return {name: [(str(k), number(s)) for k, s in zip(df.InChIkey2D, df["CSI:FingerIDScore"]) if number(s) is not None]
            for name, df in rows.groupby("challenge")}


def legacy_check(rank, row, prefix, *, pess=False):
    expected = number(row.get(prefix + "_rank"))
    assert rank["best"] == expected, (prefix, rank, expected)
    if pess:
        worst = number(row.get(prefix + "_rank_pess"))
    else:
        tied = int(row.get(prefix + "_n_tied") or 0)
        worst = expected + tied if expected is not None else None
    assert rank["worst"] == worst, (prefix, rank, worst)


def _eq(a, b):
    assert a == b, f"{a} != {b}"


# ------------------------------------------------------------------ mzDesk configurations
def mzdesk_rank(result, truth):
    """Rank of the truth in a (compact) run_query result; the ordering keys are rebuilt from the layer scores (ties -> best / worst rank)."""
    candidates = [c for c in result["candidates"] if c["layer"] != "formula_only"]
    assert [c["rank"] for c in candidates] == list(range(1, len(candidates) + 1))
    for c in candidates:
        if c.get("smiles") and c.get("ik14"):
            mol = Chem.MolFromSmiles(c["smiles"])
            assert mol is not None
            if len(Chem.GetMolFrags(mol, asMols=True)) == 1:
                assert Chem.MolToInchiKey(mol)[:14] == c["ik14"], (c["id"], c["ik14"])
    keys = [c.get("ik14") for c in candidates if c.get("ik14")]
    assert all(re.fullmatch(r"[A-Z]{14}", key) for key in keys)
    assert len(keys) == len(set(keys)), "mzDesk candidate structures were not deduplicated"
    unverifiable = [c["rank"] for c in candidates if not c.get("ik14")]
    hit = next((c for c in candidates if c.get("ik14") == truth), None)
    row = {"n": len(candidates), "top": candidates[0].get("ik14") if candidates else None, "best": None, "worst": None, "ties": 0,
           "listed": hit["rank"] if hit else None, "verifiable_n": len(keys), "unverifiable_ranks": unverifiable,
           "hit_layer": hit["layer"] if hit else None, "hit_level": hit["level"] if hit else None}
    if not hit:
        return row
    if any(c.get("detect") is None for c in candidates):
        row["best"] = row["worst"] = hit["rank"]
        row["tie_note"] = "rerank failed; displayed ordering used"
        return row

    def base(c):
        src = c["sources"] or [{"layer": c["layer"], "layer_score_norm": c["layer_score_norm"]}]
        b = min(1.0, max(float(s.get("layer_score_norm") or 0) for s in src) + 0.05 * (len({s["layer"] for s in src}) - 1))
        return max(0.0, b * 0.7 if c["hypothesis"] else b)
    bases = [base(c) for c in candidates]
    best = max(bases)
    ordering = []
    for c, b in zip(candidates, bases):
        assert abs(round(b, 4) - c["base_score"]) < 0.00011
        if b >= 0.9 * best - 1e-9:
            assert abs(round(b * c["detect"]["factor"], 4) - c["final_score"]) < 0.00011
            key = (0, -round(b * c["detect"]["factor"], 6), -round(b, 6))
        else:
            assert abs(round(b, 4) - c["final_score"]) < 0.00011
            key = (1, -round(b, 6))
        ordering.append(key)
    assert ordering == sorted(ordering), "Cannot reproduce mzDesk ordering keys"
    hkey = ordering[hit["rank"] - 1]
    positions = [i for i, key in enumerate(ordering, 1) if key == hkey]
    row.update(best=min(positions), worst=max(positions), ties=len(positions) - 1)
    return row


def fallback_rank(result, truth):
    candidates = [c for c in result["candidates"] if c["layer"] != "formula_only"]
    hit = next((c for c in candidates if c.get("ik14") == truth), None)
    keys = [c.get("ik14") for c in candidates if c.get("ik14")]
    return {"n": len(candidates), "top": candidates[0].get("ik14") if candidates else None, "best": hit["rank"] if hit else None,
            "worst": hit["rank"] if hit else None, "ties": 0, "listed": hit["rank"] if hit else None, "verifiable_n": len(keys),
            "unverifiable_ranks": [c["rank"] for c in candidates if not c.get("ik14")], "tie_note": "ordering not reproduced; displayed rank used"}


def error_rank(message):
    return {"n": 0, "top": None, "best": None, "worst": None, "ties": 0, "listed": None, "verifiable_n": 0, "unverifiable_ranks": [], "error": message}


def counts(formula):
    return Counter({el: int(num or 1) for el, num in re.findall(r"([A-Z][a-z]?)(\d*)", formula)})


def formula_rank(view, truth_formula):
    """buddy_only: rank of the truth's molecular formula in step 2's formula list (distinct formulas, in list order; the list is BUDDY's ranks
    interleaved across adducts). 'best_buddy_rank' = the best BUDDY rank within one adduct that the truth formula has."""
    if view.get("error"):
        return {"n": 0, "top": None, "best": None, "worst": None, "ties": 0, "error": view["error"]}
    target = counts(truth_formula)
    distinct, seen = [], set()
    for f in view["formulas"]:
        key = tuple(sorted(counts(f["formula"]).items()))
        if key not in seen:
            seen.add(key)
            distinct.append(f)
    pos = next((i for i, f in enumerate(distinct, 1) if counts(f["formula"]) == target), None)
    buddy_ranks = [f.get("rank") for f in view["formulas"] if counts(f["formula"]) == target and f.get("rank") is not None]
    return {"n": len(distinct), "top": distinct[0]["formula"] if distinct else None, "best": pos, "worst": pos, "ties": 0,
            "best_buddy_rank": min(buddy_ranks) if buddy_ranks else None, "retained_adducts": view.get("retained_adducts")}


# ------------------------------------------------------------------ single tools with default settings (standalone.py)
def wrong_adduct_row():
    return {"n": 0, "top": None, "best": None, "worst": None, "ties": 0, "wrong_adduct": True}


def candidate_pairs(items):
    """[(key, score)] of a ranked candidate list: one entry per IK14 (the first, i.e. best), candidates without a structure key keep their own slot."""
    pairs, seen = [], set()
    for i, c in enumerate(items):
        key = c.get("ik14") or f"?{i}"
        if key not in seen:
            seen.add(key)
            pairs.append((key, c["score"]))
    return pairs


def standalone_rows(out, n, name, truth, truth_formula, pol, truth_adduct, problems, partial):
    """-> ({config: rank row}, info) for the single tools with default settings of one challenge."""
    base = out / "standalone"
    wrong = truth_adduct != DEFAULT_ADDUCT[pol]
    rows, info = {}, {"default_adduct": DEFAULT_ADDUCT[pol], "wrong_adduct": wrong}

    def load(tool):
        path, err = base / tool / f"{name}.json", base / tool / f"{name}.error.json"
        if path.exists():
            return json.loads(path.read_text())
        if err.exists():
            problems.append(f"{n:03d} default {tool}: error: {json.loads(err.read_text()).get('error')}")
            return {"error": json.loads(err.read_text()).get("error")}
        if not partial:
            raise AssertionError(f"challenge {n}: standalone {tool} has no result")
        return None

    lib = load("library")
    if lib is not None:
        if lib.get("error"):
            rows["default_library"] = error_rank(lib["error"])
        else:
            rows["default_library"] = scored_rank(candidate_pairs(lib["compounds"]), truth)
            rows["default_library"]["n_spectra_kept"] = lib["n_spectra_kept"]
            info["library"] = {"n_spectra_kept": lib["n_spectra_kept"], "n_compounds": len(lib["compounds"]), "seconds": lib["elapsed_seconds"]}
    for tool, key in (("buddy", "default_buddy"), ("msfinder", "default_msfinder"), ("frigid", "default_frigid")):
        if tool == "frigid" and pol != "positive":
            continue
        if wrong:
            rows[key] = wrong_adduct_row()
            continue
        res = load(tool)
        if res is None:
            continue
        if res.get("error"):
            rows[key] = error_rank(res["error"])
        elif tool == "buddy":
            if res["status"] == "ok":
                rows[key] = formula_rank({"formulas": res["formulas"]}, truth_formula)
            else:                                    # e.g. precursor outside BUDDY's m/z range: no output, not an error of the run
                rows[key] = {"n": 0, "top": None, "best": None, "worst": None, "ties": 0, "buddy_status": res["status"], "note": res.get("note")}
            info["buddy_checks"] = res["checks"]
        else:
            if res["status"] not in ("done",):
                rows[key] = {**scored_rank([], truth), "status": res["status"], "note": res.get("note")}
            else:
                rows[key] = scored_rank(candidate_pairs(res["candidates"]), truth)
            info[tool] = {"status": res["status"], "seconds": res.get("seconds"), "elapsed": res["elapsed_seconds"], "source": res.get("source"),
                          "mode": res.get("mode"), "retried": res.get("retried")}
    return rows, info


# ------------------------------------------------------------------ seen sets
def seen_sets(truths_polarity):
    """Per tool: challenges whose truth structure (IK14) the tool is known to have seen. See leak_lists/README.md."""
    over = json.loads((LEAK / "overlap.json").read_text())
    lib = json.loads((LEAK / "lib_ids.json").read_text())["lib_any"]
    pos = {n for n, p in truths_polarity.items() if p == "positive"}
    S = {
        "CSI:FingerID": (set(over["sirius_pos"]) & pos) | (set(over["sirius_neg"]) - pos),   # SIRIUS fingerprint predictor of the matching polarity
        "FRIGID": set(over["msg_train"]) | set(over["msg4m"]),
        "MIST": set(over["mist_real_train_plus_aug"]),
        "JESTR": set(over["jestr_nplib1_train"]),
        "library": set(lib),                                                              # any MassBank / GNPS spectrum of the structure, any polarity
    }
    return {k: sorted(v) for k, v in S.items()}


def subsets(ids, polarity, seen):
    S = set(ids)
    s = {k: set(v) for k, v in seen.items()}
    out = {
        "all": sorted(S), "positive": sorted(n for n in ids if polarity[n] == "positive"), "negative": sorted(n for n in ids if polarity[n] == "negative"),
        "unseen_common": sorted(S - (s["CSI:FingerID"] | s["FRIGID"] | s["MIST"] | s["JESTR"] | s["library"])),
        "unseen_CSI:FingerID": sorted(S - s["CSI:FingerID"]), "unseen_JESTR": sorted(S - s["JESTR"]), "unseen_MIST": sorted(S - s["MIST"]),
        "unseen_FRIGID": sorted(S - s["FRIGID"]), "unseen_library": sorted(S - s["library"]),
        "unseen_mzDesk": sorted(S - (s["library"] | s["FRIGID"])),                       # mzDesk = library + MS-FINDER + FRIGID (+ BUDDY): the two with a known training / reference set
        "unseen_mzDesk_nolib": sorted(S - s["library"]),
    }
    return out


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=str(HERE), help="folder holding manifest.json, store/ and replay/ (default: this folder)")
    ap.add_argument("--partial", action="store_true", help="accept a subset of the 500 challenges / of the configurations (smoke tests)")
    ap.add_argument("--replace", action="store_true", help="replace the previous scores.json")
    ap.add_argument("--target", default=None, help="output file (default: DIR/scores.json)")
    args = ap.parse_args()
    out = Path(args.dir).resolve()
    manifest = json.loads((out / "manifest.json").read_text())
    ids = manifest["ids"]
    if not args.partial:
        assert ids == list(range(1, 501)), "manifest does not cover challenges 1..500"
    problems = []
    official = pd.read_csv(CASMI / "casmi2022_ground_truth.csv", dtype=str, keep_default_na=False).set_index("challenge_id")
    old_truth = pd.read_csv(SIRIUS / "truth.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("challenge_id")
    prior_s = pd.read_csv(SIRIUS / "score/per_challenge.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("challenge_id")
    prior_j = pd.read_csv(JM / "score/per_challenge.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("challenge_id")
    all_ids = list(range(1, 501))
    polarity = {n: old_truth.loc[str(n), "polarity"] for n in all_ids}
    seen = seen_sets(polarity)
    novels = novelist_pairs(ids)

    configs = [c for c in CONFIG_GROUPS if (out / "replay" / c).is_dir()]
    if not args.partial:
        assert len(configs) == len(CONFIG_GROUPS), f"missing replay folders: {set(CONFIG_GROUPS) - set(configs)}"
    std_present = (out / "standalone").is_dir()
    if not args.partial:
        assert std_present, "no standalone/ folder: run standalone.py first"
    cases = []
    for n in ids:
        name = f"challenge_{n:03d}"
        truth = official.loc[str(n), "inchikey"][:14]
        off, old = official.loc[str(n)], old_truth.loc[str(n)]
        assert re.fullmatch(r"[A-Z]{14}", truth) and truth == old["inchikey14"] and off["formula"] == old["formula"]
        pol = old["polarity"]
        meta_path = out / "store" / f"{name}.meta.json"
        assert meta_path.exists(), f"challenge {n} has no collect result"
        meta = json.loads(meta_path.read_text())
        csi = scored_rank(load_scores(JM / "inputs/candidates" / f"{name}.tsv.gz", "ikey14", "csi_score"), truth)
        jestr = scored_rank(load_scores(JM / "jestr/scores" / f"{name}.tsv.gz", "ikey14", "score"), truth)
        nov = scored_rank(novels.get(name, []), truth)
        mist_fp = scored_rank(load_scores(JM / "mist/scores" / f"{name}.tsv.gz", "ikey14", "score_fp"), truth) if pol == "positive" else None
        mist_ct = scored_rank(load_scores(JM / "mist/scores" / f"{name}.tsv.gz", "ikey14", "score_contrast"), truth) if pol == "positive" else None
        for label, fn in (("CSI vs SIRIUS table", lambda: legacy_check(csi, prior_s.loc[str(n)], "csi")),
                          ("CSI vs JM table", lambda: legacy_check(csi, prior_j.loc[str(n)], "CSI:FingerID", pess=True)),
                          ("MSNovelist vs SIRIUS table", lambda: legacy_check(nov, prior_s.loc[str(n)], "novelist")),
                          ("JESTR vs JM table", lambda: legacy_check(jestr, prior_j.loc[str(n)], "JESTR", pess=True)),
                          ("candidate counts", lambda: [_eq(csi["n"], int(prior_s.loc[str(n), "n_csi_candidates"])),
                                                        _eq(nov["n"], int(prior_s.loc[str(n), "n_novelist_candidates"])),
                                                        _eq(csi["n"], int(prior_j.loc[str(n), "CSI:FingerID_n_scored"])),
                                                        _eq(jestr["n"], int(prior_j.loc[str(n), "JESTR_n_scored"]))]),
                          ("MIST vs JM table", lambda: mist_fp is None or [
                              legacy_check(mist_fp, prior_j.loc[str(n)], "MIST-FP", pess=True),
                              legacy_check(mist_ct, prior_j.loc[str(n)], "MIST-contrastive", pess=True),
                              _eq(mist_fp["n"], int(prior_j.loc[str(n), "MIST-FP_n_scored"])),
                              _eq(mist_ct["n"], int(prior_j.loc[str(n), "MIST-contrastive_n_scored"]))])):
            try:
                fn()
            except AssertionError as exc:
                problems.append(f"{n:03d}: comparator check failed ({label}): {exc}")
        cfg_rows, cfg_info = {}, {}
        for c in configs:
            path = out / "replay" / c / f"{name}.json"
            err = out / "replay" / c / f"{name}.error.json"
            if path.exists():
                run = json.loads(path.read_text())
                if run["store_fingerprint"] != meta["code_fingerprint"] or run["code_fingerprint"] != meta["code_fingerprint"]:
                    problems.append(f"{n:03d} {c}: replay code fingerprint differs from the collect run")
                res = run["result"]
                if c == "buddy_only":
                    row = formula_rank(res, off["formula"])
                else:
                    try:
                        row = mzdesk_rank(res, truth)
                    except AssertionError as exc:
                        problems.append(f"{n:03d} {c}: mzDesk ordering check failed ({exc}); displayed rank used")
                        row = fallback_rank(res, truth)
                cfg_rows[c] = row
                info = {"replay_seconds": run["elapsed_seconds"], "step7_cache": run.get("step7_cache"), "copied_from_full": bool(run.get("copied_from_full"))}
                if c in ("full", "leakfree_full"):
                    steps = {s["id"]: s for s in res["steps"]}
                    info.update(answer_type=res["answer"]["type"], level=res["answer"]["level"], gate_passed=res["context"]["gate_passed"],
                                step_status={str(k): v["status"] for k, v in steps.items()}, failed_notes={str(k): v["note"] for k, v in steps.items() if v["status"] == "failed"},
                                top_layer=res["candidates"][0]["layer"] if res["candidates"] else None,
                                step7_engine_seconds=(res.get("isomer") or {}).get("seconds") if steps[7]["status"] != "skipped" else 0.0,
                                step6_seconds=steps[6]["seconds"], n_warnings=len(res["warnings"]))
                cfg_info[c] = info
            elif err.exists():
                cfg_rows[c] = error_rank(json.loads(err.read_text()).get("error"))
                cfg_info[c] = {"error": cfg_rows[c]["error"]}
                problems.append(f"{n:03d} {c}: replay error: {cfg_rows[c]['error']}")
            elif not args.partial:
                raise AssertionError(f"challenge {n}: configuration {c} has no result yet")
        std_info = None
        if std_present:
            std_rows, std_info = standalone_rows(out, n, name, truth, off["formula"], pol, off["adduct_official"], problems, args.partial)
            cfg_rows.update(std_rows)
        # runtime of the workflow as a user would run it: only the steps 'full' did not skip
        seconds = meta["seconds"]
        wf, wf_lf = None, None
        for key, name_ in (("wf", "full"), ("wf_lf", "leakfree_full")):
            if name_ in cfg_info and "step_status" in cfg_info[name_]:
                st = cfg_info[name_]["step_status"]
                w = {"adduct+formula": seconds.get("adduct+formula", 0.0), "library": seconds.get("library", 0.0) if name_ == "full" else 0.0,
                     "msfinder": seconds.get("msfinder", 0.0) if st["4"] != "skipped" else 0.0,
                     "frigid": seconds.get("frigid", 0.0) if st["5"] != "skipped" else 0.0,
                     "rerank": cfg_info[name_]["step6_seconds"], "isomers": cfg_info[name_]["step7_engine_seconds"] or 0.0}
                w["total"] = sum(w.values())
                if key == "wf":
                    wf = w
                else:
                    wf_lf = w
        cases.append({"id": n, "polarity": pol, "priority": off["priority"], "truth_ik14": truth, "truth_formula": off["formula"],
                      "truth_adduct": off["adduct_official"], "comparators": {"CSI:FingerID": csi, "MSNovelist": nov, "JESTR": jestr, "MIST-FP": mist_fp, "MIST-contrastive": mist_ct},
                      "configs": cfg_rows, "config_info": cfg_info, "collect": {"seconds": seconds, "frigid_devices": meta["frigid_devices"],
                      "frigid_runs": meta["frigid_runs"], "layers": meta["layers"], "library_gate_passed": (meta["library_gate"] or {}).get("passed"),
                      "library_n_candidates": (meta["library_gate"] or {}).get("n_candidates"),
                      "msfinder_good_enough": meta["msfinder_good_enough"], "code_fingerprint": meta["code_fingerprint"], "attempts": meta["attempts"],
                      "busy_layers_after_retries": meta.get("busy_layers_after_retries")}, "workflow_seconds": wf, "workflow_seconds_leakfree": wf_lf, "standalone": std_info})
    # leakfree_full == full where the stored library layer has no candidate: checked on every such challenge that was really replayed both ways
    def sig(run):
        r = run["result"]
        return ([(c["id"], c["rank"], c["layer"], c["level"], c["ik14"], c["smiles"], c["base_score"], c["final_score"]) for c in r["candidates"]], r["answer"],
                [(s_["id"], s_["status"]) for s_ in r["steps"] if s_["id"] >= 4],
                (r["isomer"] or {}).get("verdict"), (r["isomer"] or {}).get("top_id"), (r["isomer"] or {}).get("group"))
    cmp_ids, cmp_diff = [], []
    for c in cases:
        info = c["config_info"].get("leakfree_full") or {}
        if info.get("copied_from_full") or "full" not in c["config_info"] or "leakfree_full" not in c["config_info"] or c["collect"]["library_n_candidates"] != 0:
            continue
        nm = f"challenge_{c['id']:03d}.json"
        a, b = json.loads((out / "replay/full" / nm).read_text()), json.loads((out / "replay/leakfree_full" / nm).read_text())
        cmp_ids.append(c["id"])
        if sig(a) != sig(b):
            cmp_diff.append(c["id"])
    fingerprints = sorted({c["collect"]["code_fingerprint"] for c in cases})
    sub = subsets(ids, polarity, seen)
    sub["default_adduct"] = sorted(c["id"] for c in cases if not (c["standalone"] or {}).get("wrong_adduct", True))
    sub["default_adduct_wrong"] = sorted(c["id"] for c in cases if (c["standalone"] or {}).get("wrong_adduct", False))
    verify = [json.loads(p.read_text()) for p in sorted((out / "standalone" / "frigid_verify").glob("challenge_*.json"))] if std_present else []
    standalone_meta = {"configs": list(DEFAULT_CONFIGS) if std_present else [], "default_adduct": DEFAULT_ADDUCT,
                       "wall_seconds": _wall(out / "standalone.log"), "frigid_verify": verify,
                       "buddy_cache_checks": dict(Counter(f"{k}={v}" for c in cases for k, v in ((c["standalone"] or {}).get("buddy_checks") or {}).items())),
                       "frigid_sources": dict(Counter((c["standalone"] or {}).get("frigid", {}).get("source") or "-" for c in cases if (c["standalone"] or {}).get("frigid")))}
    log_wall = lambda name: _wall(out / name)
    meta = {"n": len(cases), "ids": ids, "settings": manifest["settings"], "configs": configs, "config_groups": CONFIG_GROUPS,
            "positive_only_configs": list(POSITIVE_ONLY), "code_fingerprints": fingerprints, "code_root": manifest.get("code_root"),
            "collect_wall_seconds": log_wall("collect.log"), "replay_wall_seconds": log_wall("replay.log"), "problems": problems, "partial": args.partial,
            "leakfree_equivalence": {"compared": len(cmp_ids), "different": cmp_diff},
            "seen": seen, "subsets": sub, "standalone": standalone_meta,
            "seen_definition": {"CSI:FingerID": "overlap.json sirius_pos (positive) / sirius_neg (negative): SIRIUS fingerprint predictor training structures of the matching polarity",
                                "FRIGID": "overlap.json msg_train U msg4m (MIST encoder training + the 4M-structure set of the diffusion model); FRIGID's DLM pre-training corpus is NOT excludable",
                                "MIST": "overlap.json mist_real_train_plus_aug (real training spectra + augmented)", "JESTR": "overlap.json jestr_nplib1_train",
                                "library": "lib_ids.json lib_any: the structure has at least one spectrum (any polarity) in the local MassBank or GNPS index",
                                "MSNovelist": "NOT excludable (training structures not available)"}}
    target = Path(args.target) if args.target else out / "scores.json"
    if target.exists() and not args.replace:
        raise FileExistsError(target)
    target.write_text(json.dumps({"meta": meta, "cases": cases}, indent=1, ensure_ascii=False, allow_nan=False) + "\n")
    print(f"{len(cases)} challenges, {len(configs)} configurations -> {target}")
    for c in configs + (list(DEFAULT_CONFIGS) if std_present else []):
        rows = [x["configs"][c] for x in cases if c in x["configs"] and not (c in POSITIVE_ONLY and x["polarity"] != "positive")]
        top = lambda k: sum(r["worst"] is not None and r["worst"] <= k for r in rows)
        print(f"  {c:18} n={len(rows):3} Top-1 {top(1):3} Top-3 {top(3):3} Top-10 {top(10):3}")
    print(f"{len(problems)} consistency problems" + "".join("\n  " + p for p in problems[:20]))


def _wall(log):
    if not log.exists():
        return None
    return sum(float(line.split()[1]) for line in log.read_text().splitlines() if line.startswith("WALL_SECONDS "))


if __name__ == "__main__":
    main()
