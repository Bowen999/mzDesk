"""Isomer check (workflow step 7) and evidence levels on CASMI 2022, from the saved mzDesk v0.5.3 run.

Reads, and never re-runs anything:
  benchmark/casmi2022/mzdesk_v0.5.3_20261102/scores.json          truth, ranks, answer type and level per challenge
  benchmark/casmi2022/mzdesk_v0.5.3_20261102/replay/full/*.json   workflow output of the default gated workflow

Scoring rules (as in the paper): a structure is correct when the first 14 characters of its InChIKey match the
true structure; tied ranks count at the worst rank (best rank shown for reference); a challenge without output is a
miss; the denominator is always all 500 challenges (or all challenges of the row). 95 % CI: Wilson score interval.

"Answer-level" correct: the true structure is the single answer or one of the structures of the group answer
(set membership, so ties in score do not apply). Formula-only and mass-only answers are never correct.

Usage: python3 validation/isomer_check/analyze_isomer_check.py
"""
import collections
import json
import math
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "benchmark/casmi2022/mzdesk_v0.5.3_20261102"
TYPES = ["single", "group", "formula_only", "mass_only"]
LEVELS = ["2a", "3", "H", "4", "5"]
SPOT = {"Fig. 2c, d": 376, "single answer, separated from every isomer, correct": 86,
        "group answer, truth inside but not at the top": 1, "group answer, truth not inside": 2}


def wilson(k, n, z=1.959964):
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def pct(k, n):
    lo, hi = wilson(k, n)
    return f"{k:4d}/{n:<4d} {100 * k / n:5.1f} % [{100 * lo:4.1f}, {100 * hi:5.1f}]" if n else f"{k:4d}/{n:<4d}   n/a"


def load():
    scores = json.loads((RUN / "scores.json").read_text())
    cases = scores["cases"]
    rows = []
    for c in cases:
        cid = c["id"]
        rep = json.loads((RUN / "replay/full" / f"challenge_{cid:03d}.json").read_text())
        assert rep["challenge_id"] == cid and rep["config"] == "full", cid
        assert rep["options"] == {"layers": ["library", "msfinder", "frigid"], "forced": False, "step7": True,
                                  "kind": "structure"}, (cid, rep["options"])
        res = rep["result"]
        ans, cands, iso = res["answer"], res["candidates"], res.get("isomer")
        by_id = {x["id"]: x for x in cands}
        f, info = c["configs"]["full"], c["config_info"]["full"]
        truth = c["truth_ik14"]
        top = by_id[ans["top_id"]] if ans.get("top_id") in by_id else None
        members = [by_id[i] for i in ans["ids"]] if ans["type"] in ("single", "group") else []
        tcand = [x for x in cands if x.get("ik14") == truth]
        tcand = tcand[0] if tcand else None
        step7 = next(s for s in res["steps"] if s["id"] == 7)
        notes = iso["notes"] if iso else []
        n_considered = next((int(m.group(1)) for n in notes
                             for m in [re.match(r"Relative to the (\d+) candidates considered", n)] if m), None)
        n_not_compared = sum(int(m.group(1)) for n in notes
                             for m in [re.match(r"(\d+) lower-ranked isomer\(s\) are not compared", n)] if m)
        tol = sorted({m.group(1) for n in notes for m in [re.search(r"fragments within ([\d.]+) Da", n)] if m})
        rows.append(dict(
            id=cid, polarity=c["polarity"], truth=truth, truth_formula=c["truth_formula"],
            type=ans["type"], level=ans["level"], verdict=ans.get("isomer_verdict"),
            s_type=info["answer_type"], s_level=info["level"], s_top=f["top"], s_layer=info["top_layer"],
            layer=ans.get("layer"), best=f["best"], worst=f["worst"], ties=f["ties"], n=f["n"],
            top=top, members=members, member_iks=[x["ik14"] for x in members], cands=cands, tcand=tcand,
            step7_status=step7["status"], step7_note=step7["note"], iso=iso, notes=notes,
            n_considered=n_considered, n_not_compared=n_not_compared, tol=tol, statement=ans["statement"],
            formula_ids=[x["formula"] for x in cands if x["layer"] == "formula_only"],
        ))
    return scores, rows


def derive(r):
    r["top1_w"] = r["worst"] == 1
    r["top1_b"] = r["best"] == 1
    r["top10_w"] = r["worst"] is not None and r["worst"] <= 10
    r["top_is_truth"] = r["top"] is not None and r["top"]["ik14"] == r["truth"]
    r["in_answer"] = r["truth"] in r["member_iks"]
    # where the truth sits relative to a structure answer that does not contain it
    if r["type"] in ("single", "group") and not r["in_answer"]:
        t, top = r["tcand"], r["top"]
        if t is None:
            r["miss_where"] = "not among the candidates"
        elif t["formula"] != top["formula"]:
            r["miss_where"] = "candidate with another formula than the top"
        else:
            r["miss_where"] = "isomer of the top, outside the answer"
    else:
        r["miss_where"] = None
    # position of the truth in the step-7 list (top first, then same-formula candidates by rank, at most 10).
    # Inferred from the rank order; the replay keeps the group but not the per-isomer list.
    t = r["tcand"]
    if t is not None and r["top"] is not None and t["formula"] == r["top"]["formula"]:
        r["iso_pos"] = 1 + sum(1 for x in r["cands"] if x["rank"] < t["rank"] and x["formula"] == r["top"]["formula"])
    else:
        r["iso_pos"] = None


def main():
    scores, rows = load()
    for r in rows:
        derive(r)
    N = len(rows)
    out = []
    p = out.append
    meta = scores["meta"]
    p("Isomer check and evidence levels on CASMI 2022 (mzDesk v0.5.3, config 'full' = default gated workflow)")
    p(f"Run: {RUN.relative_to(ROOT)}; code fingerprint {meta['code_fingerprints']}; "
      f"settings.local_isomers = {meta['settings']['local_isomers']}; fragment tol (run) {meta['settings']['tol_da']} Da")
    tols = collections.Counter(t for r in rows for t in r["tol"])
    p(f"Step-7 fragment tolerance recorded in the replay notes: {dict(tols)} (challenges)")
    p("CI: Wilson 95 %. Top-1 at the worst rank of ties unless marked 'best'.")
    p("")

    # ---------------------------------------------------------------- self-checks
    p("== Self-checks ==")
    t1 = sum(r["top1_w"] for r in rows)
    t10 = sum(r["top10_w"] for r in rows)
    p(f"n = {N} (expect 500): {'PASS' if N == 500 else 'FAIL'}")
    p(f"Top-1 (worst) = {t1} (paper 89): {'PASS' if t1 == 89 else 'FAIL'}; "
      f"Top-10 (worst) = {t10} (paper 159): {'PASS' if t10 == 159 else 'FAIL'}; "
      f"Top-1 (best) = {sum(r['top1_b'] for r in rows)}")
    bad = collections.defaultdict(list)
    for r in rows:
        if r["type"] != r["s_type"]:
            bad["answer type"].append(r["id"])
        if r["level"] != r["s_level"]:
            bad["level"].append(r["id"])
        if (r["top"]["ik14"] if r["top"] else None) != r["s_top"]:
            bad["top structure"].append(r["id"])
        if r["type"] in ("single", "group") and r["top"]["id"] not in [m["id"] for m in r["members"]]:
            bad["top in its answer"].append(r["id"])
        if r["type"] == "single" and len(r["members"]) != 1:
            bad["single has 1 member"].append(r["id"])
        if r["type"] == "group" and (len(r["members"]) < 2 or sorted(r["iso"]["group"]) != sorted(m["id"] for m in r["members"])):
            bad["group = step-7 group"].append(r["id"])
        if r["type"] == "group" and len({m["formula"] for m in r["members"]}) != 1:
            bad["group has one formula"].append(r["id"])
        if r["type"] in ("single", "group") and r["level"] != r["top"]["level"]:
            bad["answer level = top level"].append(r["id"])
        rk = r["tcand"]["rank"] if r["tcand"] else None
        if (rk is None) != (r["best"] is None) or (rk is not None and not r["best"] <= rk <= r["worst"]):
            bad["replay rank of truth within [best, worst]"].append(r["id"])
        if r["top1_w"] and not r["in_answer"]:
            bad["Top-1 correct => truth in answer"].append(r["id"])
    for k in ["answer type", "level", "top structure", "top in its answer", "single has 1 member",
              "group = step-7 group", "group has one formula", "answer level = top level",
              "replay rank of truth within [best, worst]", "Top-1 correct => truth in answer"]:
        p(f"replay vs scores.json, {k}: {'PASS' if not bad[k] else 'FAIL ' + str(bad[k][:20])}")
    r376 = next(r for r in rows if r["id"] == 376)
    ok = r376["type"] == "group" and len(r376["members"]) == 4 and r376["in_answer"] and r376["top1_w"]
    p(f"Challenge 376 (Fig. 2d: truth at rank 1 plus 3 isomers not separated): {'PASS' if ok else 'FAIL'}; "
      f"this run: {r376['type']} of {len(r376['members'])} {r376['member_iks']}, truth {r376['truth']} at worst rank "
      f"{r376['worst']}; {r376['notes'][0] if r376['notes'] else ''}")
    fig = ROOT / "paper/figures/figure_2/challenge376_job.json"
    if fig.exists():
        fj = json.loads(fig.read_text())["result"]
        fid = {c["id"]: c for c in fj["candidates"]}
        p(f"  Fig. 2 source run ({fig.relative_to(ROOT)}): group of {len(fj['isomer']['group'])} "
          f"{[fid[i].get('inchikey', '')[:14] for i in fj['isomer']['group']]}, step-7 settings {fj['isomer'].get('settings')}")
    p("")

    # ---------------------------------------------------------------- 1, 2: answer types
    p("== 1-2. Answer types (n = 500) ==")
    p(f"{'type':<13} {'n':>4} {'%':>6}  {'Top-1 (worst rank)':<34} {'Top-1 (best rank)':<34} truth in answer")
    for t in TYPES:
        g = [r for r in rows if r["type"] == t]
        p(f"{t:<13} {len(g):4d} {100 * len(g) / N:5.1f}%  {pct(sum(r['top1_w'] for r in g), len(g)):<34} "
          f"{pct(sum(r['top1_b'] for r in g), len(g)):<34} {pct(sum(r['in_answer'] for r in g), len(g))}")
    p(f"{'all':<13} {N:4d} 100.0%  {pct(t1, N):<34} {pct(sum(r['top1_b'] for r in rows), N):<34} "
      f"{pct(sum(r['in_answer'] for r in rows), N)}")
    p("")
    p("Step-7 verdict behind each answer type:")
    vc = collections.Counter((r["type"], r["verdict"], r["step7_status"]) for r in rows)
    for (t, v, s), k in sorted(vc.items(), key=lambda x: (TYPES.index(x[0][0]), -x[1])):
        g = [r for r in rows if (r["type"], r["verdict"], r["step7_status"]) == (t, v, s)]
        p(f"  {t:<13} verdict={str(v):<18} step 7 {s:<8} n={k:4d}  Top-1 {pct(sum(r['top1_w'] for r in g), k)}")
    why = collections.Counter(re.sub(r"\(.*", "(...)", re.sub(r"formula \S+", "formula X", r["step7_note"])).strip()
                              for r in rows if r["verdict"] == "insufficient")
    p("  why 'insufficient' (no comparison made):")
    for k, v in why.most_common():
        p(f"    {v:3d}  {k}")
    p("")

    # ---------------------------------------------------------------- 3: groups
    p("== 3. Group answers ==")
    G = [r for r in rows if r["type"] == "group"]
    sizes = np.array([len(r["members"]) for r in G])
    q1, med, q3 = np.percentile(sizes, [25, 50, 75])
    p(f"n groups = {len(G)}; size (structures, top included): median {med:g}, IQR {q1:g}-{q3:g}, "
      f"range {sizes.min()}-{sizes.max()}")
    p("  size distribution: " + ", ".join(f"{s}: {k}" for s, k in sorted(collections.Counter(sizes.tolist()).items())))
    p(f"  groups where some same-formula isomers were not compared (10-candidate cap): "
      f"{sum(r['n_not_compared'] > 0 for r in G)}")
    gin = sum(r["in_answer"] for r in G)
    p(f"Groups that contain the true structure: {pct(gin, len(G))}")
    p(f"  of which the true structure is the top (worst rank 1): {sum(r['top1_w'] for r in G)}; "
      f"elsewhere in the group: {gin - sum(r['top1_w'] for r in G)}")
    wrong = [r for r in rows if not r["top1_w"]]
    win = [r for r in wrong if r["in_answer"]]
    p(f"Top-1 wrong (all answer types): {len(wrong)}; true structure inside the reported group: {pct(len(win), len(wrong))}")
    gw = [r for r in G if not r["top1_w"]]
    p(f"  among group answers with Top-1 wrong: {pct(sum(r['in_answer'] for r in gw), len(gw))}")
    p("  of these, true structure tied with the top at rank 1 (best 1, worst > 1): "
      f"{sum(1 for r in win if r['best'] == 1)}")
    ws = np.array([len(r["members"]) for r in win])
    p(f"  size of these {len(win)} groups: median {np.median(ws):g}, range {ws.min()}-{ws.max()}")
    p("Group answers without the true structure, where it is:")
    for k, v in collections.Counter(r["miss_where"] for r in G if not r["in_answer"]).most_common():
        p(f"  {v:4d}  {k}")
    p("")

    # ---------------------------------------------------------------- isomer check, by verdict
    p("== Isomer check: does the verdict match the truth? (structure answers, n = %d) ==" %
      sum(r["type"] in ("single", "group") for r in rows))
    S = [r for r in rows if r["type"] in ("single", "group")]
    p(f"{'verdict':<20} {'n':>4}  {'top = truth (worst rank 1)':<34} truth in answer")
    for v in ["single", "insufficient", "indistinguishable"]:
        g = [r for r in S if r["verdict"] == v]
        p(f"{v:<20} {len(g):4d}  {pct(sum(r['top1_w'] for r in g), len(g)):<34} {pct(sum(r['in_answer'] for r in g), len(g))}")
    iso_truth = [r for r in S if r["tcand"] is not None and r["top"]["ik14"] != r["truth"]
                 and r["tcand"]["formula"] == r["top"]["formula"]]
    p(f"Top is not the truth, truth is a same-formula candidate (an isomer of the top): {len(iso_truth)}")
    p(f"  kept in the answer (check could not separate it from the top): {sum(r['in_answer'] for r in iso_truth)}")
    out_ = [r for r in iso_truth if not r["in_answer"]]
    p(f"  left out of the answer: {len(out_)} -> {[r['id'] for r in out_]}")
    sep = [r for r in out_ if r["n_considered"] and r["iso_pos"] <= r["n_considered"]]
    cap = [r for r in out_ if r["n_considered"] and r["iso_pos"] > r["n_considered"]]
    nochk = [r for r in out_ if not r["n_considered"]]
    p(f"    compared and separated from the top (position in the step-7 list <= candidates considered): "
      f"{len(sep)} {[r['id'] for r in sep]} (in single answers: {[r['id'] for r in sep if r['type'] == 'single']})")
    p(f"    beyond the 10-candidate cap, not compared: {len(cap)} {[r['id'] for r in cap]}")
    p(f"    no comparison made (step 7 skipped): {len(nochk)} {[r['id'] for r in nochk]}")
    for r in out_:
        p(f"    {r['id']}: {r['type']}, verdict {r['verdict']}, truth rank {r['tcand']['rank']}, position in step-7 list "
          f"{r['iso_pos']}, candidates considered {r['n_considered']}, not compared (cap) {r['n_not_compared']}; "
          f"step 7: {r['step7_note']}")
    top_ok = [r for r in S if r["top1_w"]]
    p(f"Top is the truth (worst rank 1, n = {len(top_ok)}), verdict: " + ", ".join(
        f"{k} {v}" for k, v in collections.Counter(r["verdict"] for r in top_ok).most_common()))
    sv = [r for r in S if r["verdict"] == "single"]
    p(f"Verdict 'single' (n = {len(sv)}), where the truth is: correct {sum(r['top1_w'] for r in sv)}; " + ", ".join(
        f"{k} {v}" for k, v in collections.Counter(r["miss_where"] for r in sv if not r["in_answer"]).most_common()))
    p(f"Verdict 'single' with 'another isomer explains peaks the top-1 cannot': "
      f"{sum('Another isomer explains' in r['step7_note'] for r in sv)}; in group answers: "
      f"{sum('Another isomer explains' in r['step7_note'] for r in G)}")
    p("Single answers, where the truth is when it is not the answer:")
    for k, v in collections.Counter(r["miss_where"] for r in rows if r["type"] == "single" and not r["in_answer"]).most_common():
        p(f"  {v:4d}  {k}")
    p("")

    # ---------------------------------------------------------------- 4: answer-level accuracy
    p("== 4. Answer-level accuracy ==")
    a = sum(r["in_answer"] for r in rows)
    p(f"Truth is the single answer or inside the group answer: {pct(a, N)}")
    p(f"Top-1 (worst rank): {pct(t1, N)}; difference {a - t1:+d}")
    p(f"  single answers correct: {sum(r['in_answer'] for r in rows if r['type'] == 'single')}; "
      f"group answers containing the truth: {gin}")
    p("")

    # ---------------------------------------------------------------- 5: levels
    p("== 5. Evidence level of the answer (n = 500) ==")
    p(f"{'level':<6} {'n':>4}  {'Top-1 (worst rank)':<34} {'Top-1 (best rank)':<34} truth in answer")
    for L in LEVELS:
        g = [r for r in rows if r["level"] == L]
        p(f"{L:<6} {len(g):4d}  {pct(sum(r['top1_w'] for r in g), len(g)):<34} "
          f"{pct(sum(r['top1_b'] for r in g), len(g)):<34} {pct(sum(r['in_answer'] for r in g), len(g))}")
    acc = {L: sum(r["top1_w"] for r in rows if r["level"] == L) / max(1, sum(r["level"] == L for r in rows))
           for L in LEVELS}
    order = all(acc[LEVELS[i]] >= acc[LEVELS[i + 1]] for i in range(len(LEVELS) - 1))
    p(f"Top-1 accuracy non-increasing along 2a > 3 > H > 4 > 5: {order} "
      f"({', '.join(f'{L} {100 * acc[L]:.1f} %' for L in LEVELS)})")
    p("Level x answer type (n, Top-1 worst):")
    for L in LEVELS:
        cells = []
        for t in TYPES:
            g = [r for r in rows if r["level"] == L and r["type"] == t]
            if g:
                cells.append(f"{t} {len(g)} ({sum(r['top1_w'] for r in g)})")
        p(f"  {L:<3} " + "; ".join(cells))
    p("Level 2a, top layer: " + ", ".join(
        f"{k} {v}" for k, v in collections.Counter(r["layer"] for r in rows if r["level"] == "2a").items()))
    m2a = [r for r in rows if r["level"] == "2a" and not r["top1_w"]]
    p(f"Level 2a with Top-1 wrong (n = {len(m2a)}): truth in answer {sum(r['in_answer'] for r in m2a)}; "
      f"tied at rank 1 {sum(r['best'] == 1 for r in m2a)}; otherwise " + ", ".join(
          f"{k} {v}" for k, v in collections.Counter(r["miss_where"] for r in m2a if not r["in_answer"]).most_common()))
    p(f"  of the {len(m2a)}, top has the true formula (a library isomer of the truth): "
      f"{sum(r['top']['formula'] == r['truth_formula'] for r in m2a)}")
    p("")

    # ---------------------------------------------------------------- 6: polarity
    p("== 6. By polarity ==")
    for pol in ["positive", "negative"]:
        P = [r for r in rows if r["polarity"] == pol]
        p(f"-- {pol} (n = {len(P)}); Top-1 {pct(sum(r['top1_w'] for r in P), len(P))}; "
          f"answer-level {pct(sum(r['in_answer'] for r in P), len(P))}")
        for t in TYPES:
            g = [r for r in P if r["type"] == t]
            p(f"   type {t:<13} n={len(g):4d}  Top-1 {pct(sum(r['top1_w'] for r in g), len(g)):<34} "
              f"in answer {pct(sum(r['in_answer'] for r in g), len(g))}")
        for L in LEVELS:
            g = [r for r in P if r["level"] == L]
            p(f"   level {L:<12} n={len(g):4d}  Top-1 {pct(sum(r['top1_w'] for r in g), len(g))}")
    p("")

    # ---------------------------------------------------------------- formula-only answers
    F = [r for r in rows if r["type"] == "formula_only"]
    p("== Formula-only answers: rank of the true formula in the listed formulas ==")
    p("  " + ", ".join(f"{r['id']}: {r['formula_ids'].index(r['truth_formula']) + 1 if r['truth_formula'] in r['formula_ids'] else '-'}"
                       f"/{len(r['formula_ids'])}" for r in F))
    p("")

    # ---------------------------------------------------------------- spot checks
    p("== Spot checks (fields as read) ==")
    by = {r["id"]: r for r in rows}
    for label, cid in SPOT.items():
        r = by[cid]
        p(f"-- challenge {cid} ({label})")
        p(f"   scores.json: truth_ik14={r['truth']} polarity={r['polarity']} answer_type={r['s_type']} "
          f"level={r['s_level']} top={r['s_top']} best={r['best']} worst={r['worst']} ties={r['ties']}")
        p(f"   replay: answer.type={r['type']} answer.ids={[m['id'] for m in r['members']]} "
          f"isomer_verdict={r['verdict']} isomer.group={r['iso']['group'] if r['iso'] else None}")
        p(f"   members (id rank ik14 formula layer level): " + "; ".join(
            f"{m['id']} {m['rank']} {m['ik14']}{'*' if m['ik14'] == r['truth'] else ''} {m['formula']} {m['layer']} {m['level']}"
            for m in r["members"]))
        t = r["tcand"]
        p(f"   truth in candidates: {(t['id'], t['rank'], t['formula'], t['layer']) if t else None}; "
          f"truth in answer: {r['in_answer']}")
        p(f"   step 7: {r['step7_status']}: {r['step7_note']}")
        for n in r["notes"]:
            p(f"   note: {n}")
    p("")

    text = "\n".join(out)
    print(text)
    return text


if __name__ == "__main__":
    main()
