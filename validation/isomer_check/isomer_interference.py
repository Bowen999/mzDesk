"""Isomer interference and the isomer check (step 7) on CASMI 2022 (mzDesk v0.5.3, default workflow).

Isomer interference: the final candidate list holds another structure (different first InChIKey block) with the same
molecular formula as the top-1. Correct: first 14 characters of the InChIKey match; tied ranks count at the worst rank.
Challenges without a structure answer (formula only, mass only) are left out of the interference tables and count as
'truth not in answer' in the answer-level table. The answer is the top-1 alone ('single') or the group of structures
the isomer check could not separate ('group').

Usage: python3 validation/isomer_check/isomer_interference.py
"""
import json
import math
import statistics
from collections import Counter
from pathlib import Path

RUN = Path(__file__).resolve().parents[2] / "benchmark" / "casmi2022" / "mzdesk_v0.5.3_20261102"


def wilson(k, n, z=1.959964):
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return 100 * (centre - half), 100 * (centre + half)


def row(label, k, n, width=20):
    lo, hi = wilson(k, n)
    return f"{label:<{width}}{n:>4}{k:>16}   {100 * k / n:.1f} ({lo:.1f}-{hi:.1f})"


def load():
    rows = []
    for case in json.loads((RUN / "scores.json").read_text())["cases"]:
        result = json.loads((RUN / "replay" / "full" / f"challenge_{case['id']:03d}.json").read_text())["result"]
        answer, iso = result["answer"], result["isomer"] or {}
        cands = [c for c in sorted(result["candidates"], key=lambda c: c["rank"]) if c.get("ik14")]
        by_id = {c["id"]: c for c in cands}
        r = dict(id=case["id"], type=answer["type"], verdict=iso.get("verdict"),
                 correct=case["configs"]["full"].get("worst") == 1,
                 in_answer=any(by_id[i]["ik14"] == case["truth_ik14"] for i in answer.get("ids") or [] if i in by_id),
                 size=len(answer.get("ids") or []))
        if answer["type"] in ("single", "group"):
            top = cands[0]
            truth = next((c for c in cands if c["ik14"] == case["truth_ik14"]), None)
            r.update(layer=top["layer"],
                     interference=any(c["formula"] == top["formula"] and c["ik14"] != top["ik14"] for c in cands),
                     truth_where=("top" if truth is top else "absent" if truth is None
                                  else "isomer" if truth["formula"] == top["formula"] else "other formula"))
        rows.append(r)
    return rows


def main():
    rows = load()
    assert len(rows) == 500
    struct = [r for r in rows if r["type"] in ("single", "group")]

    print("1. Top-1 with and without isomer interference (structure answers)")
    print("Isomer interference   n   Top-1 correct   % (95 % CI)")
    for flag, label in ((False, "No"), (True, "Yes")):
        sub = [r for r in struct if r["interference"] == flag]
        print(row(label, sum(r["correct"] for r in sub), len(sub)))

    print("\n2. The same, by the layer of the top-1")
    print("Layer      Interference   n   Top-1 correct   % (95 % CI)")
    for layer in ("library", "msfinder", "frigid"):
        for flag, label in ((False, "No"), (True, "Yes")):
            sub = [r for r in struct if r["layer"] == layer and r["interference"] == flag]
            print(row(f"{layer:<11}{label}", sum(r["correct"] for r in sub), len(sub)))

    print("\n3. Answer level (all 500 challenges)")
    print(f"Top-1 correct: {sum(r['correct'] for r in rows)}; true structure in the answer: {sum(r['in_answer'] for r in rows)}")
    groups = [r["size"] for r in rows if r["type"] == "group"]
    q = statistics.quantiles(groups, n=4)
    print(f"Group answers: {len(groups)}; size median {statistics.median(groups):g} (IQR {q[0]:g}-{q[2]:g}), "
          f"range {min(groups)}-{max(groups)}, at the 10-candidate cap {groups.count(10)}")

    print("\n4. Top-1 wrong, true structure an isomer of the top-1 (same formula, in the list)")
    sub = [r for r in struct if r["truth_where"] == "isomer"]
    kept = [r for r in sub if r["in_answer"]]
    out = [r for r in sub if not r["in_answer"]]
    print(f"n {len(sub)}; kept in the group answer {len(kept)}; left out {len(out)} "
          f"(separated by the check {sum(r['verdict'] in ('single', 'indistinguishable') for r in out)}, "
          f"no comparison made {sum(r['verdict'] not in ('single', 'indistinguishable') for r in out)})")

    print("\n5. Verdict 'single' (the check separated the top-1 from every isomer compared): where the true structure is")
    sub = [r for r in struct if r["verdict"] == "single"]
    print(f"n {len(sub)}: " + ", ".join(f"{k} {v}" for k, v in Counter(r["truth_where"] for r in sub).most_common()))


if __name__ == "__main__":
    main()
