#!/usr/bin/env python3
"""Check that replay == a normal run, for some challenges (the stored outputs come from collect.py).

  MZDESK_CODE_ROOT=... python3 benchmark/casmi2022/mzdesk_ablation_20261030/verify.py --ids 9,16,116,227,370 [--store DIR] [--out DIR] [--no-fresh]

1. replay == normal run on the SAME engine outputs. A normal pipeline.run_query is run on the blind form with the engines running for real,
   but with recording pass-through wrappers (and a step-7 cache); then replay_query is fed what was recorded and must reproduce the whole
   result dictionary exactly (apart from wall-clock seconds and the library step's note, which counts the compounds listed). This is the check of
   the replay mechanism; it cannot be disturbed by engine non-determinism. Also compared, for information only: the collect.py outputs
   against this normal run (MS-FINDER and FRIGID are not bit-reproducible: another FRIGID device, tied MS-FINDER scores). A normal run
   takes up to a few minutes per challenge. --no-fresh skips this part.
2. replay 'full' with step 7 against 'full' with step 7 stubbed out (the ablation setting): the ranked candidate list (ids, order, levels,
   scores) must be identical, i.e. leaving step 7 out of the ablations cannot change a rank.
3. step 7 of 'full' must not have failed (a failed isomer step means MS-FINDER did not run). leakfree_full (library not searched): step 3 skipped, gate not
   passed, no candidate from the library, MS-FINDER step run, step 7 not failed. The truth is never read.
Exit status 1 when a check fails.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import common
from common import HERE

VOLATILE = {"seconds", "log"}


def norm(res, drop_step_notes=()):
    """JSON round trip (tuples -> lists, as the saved files), without wall-clock fields."""
    data = json.loads(json.dumps(res, default=str))

    def strip(o):
        if isinstance(o, dict):
            return {k: strip(v) for k, v in o.items() if k not in VOLATILE}
        if isinstance(o, list):
            return [strip(v) for v in o]
        return o
    data = strip(data)
    for s in data.get("steps", []):
        if s["id"] in drop_step_notes:
            s["note"] = "<ignored>"
    return data


def diff(a, b, path=""):
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                out.append(f"{path}/{k}: only in {'A' if k in a else 'B'}")
            else:
                out += diff(a[k], b[k], f"{path}/{k}")
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append(f"{path}: length {len(a)} vs {len(b)}")
        for i, (x, y) in enumerate(zip(a, b)):
            out += diff(x, y, f"{path}[{i}]")
    elif a != b:
        if isinstance(a, float) and isinstance(b, float) and abs(a - b) < 1e-9:
            return out
        out.append(f"{path}: {str(a)[:90]!r} vs {str(b)[:90]!r}")
    return out



def fresh_recorded(rec, cache):
    """A normal pipeline.run_query (engines run for real, only the ones the workflow decides to run) with recording pass-throughs.
    Returns (result, record) where record has the layout of collect.py's records (engines the workflow did not call are marked 'not called')."""
    import copy
    from mzdesk.annotation import formula_step, frigid_step, library_step, msfinder_step, pipeline
    from mzdesk.exclusivity import isomers
    from replay import patched
    got = {"args": {}}
    real = dict(af=formula_step.run_adduct_and_formula, lib=library_step.run_library, msf=msfinder_step.run_msfinder_layer, fr=frigid_step.run_frigid_layers)

    def keep(key):
        def wrap(*a, **k):
            out = real[key](*a, **k)
            got[key] = ("ok", copy.deepcopy(out))
            if key == "msf":
                got["args"]["msfinder"] = (a[0], a[1], a[2], [f.to_dict() for f in a[3]], k)
            if key == "fr":
                got["args"]["frigid"] = (a[0], a[1], a[2], k)
            return out
        return wrap
    pairs = [(formula_step, "run_adduct_and_formula", keep("af")), (library_step, "run_library", keep("lib")),
             (msfinder_step, "run_msfinder_layer", keep("msf")), (frigid_step, "run_frigid_layers", keep("fr")), (isomers, "discriminate", cache)]
    with patched(*pairs):
        res = pipeline.run_query(pipeline.parse_inputs(rec["form"]))
    for key in ("af", "lib", "msf", "fr"):
        got.setdefault(key, ("exc", "NotCalled", "the normal run did not call this engine"))
    got.update(n=rec["n"], form=rec["form"])
    return res, got


def ranked(res):
    return [(c["id"], c["rank"], c["layer"], c["level"], c["ik14"], c["base_score"], c["final_score"]) for c in res["candidates"]]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ids", required=True)
    ap.add_argument("--out", default=str(HERE), help="folder holding step7_cache/ (default: this folder)")
    ap.add_argument("--store", default=None, help="folder holding store/ (default: --out)")
    ap.add_argument("--no-fresh", action="store_true")
    args = ap.parse_args()
    common.require_code_root()
    os.environ.setdefault("MZDESK_FRIGID_THREADS", "4")
    import replay
    from mzdesk.annotation import pipeline
    out = Path(args.out).resolve()
    store = (Path(args.store).resolve() if args.store else out) / "store"
    cache = replay.IsomerCache(out / "step7_cache")
    failed = 0
    for n in common.parse_ids(args.ids):
        rec = replay.load_store(store, n)
        t = time.perf_counter()
        full = replay.replay_query(rec, iso_cache=cache, display=True)
        t_replay = time.perf_counter() - t
        print(f"== {n}: replay full {t_replay:.1f}s  devices={sorted({r.get('device') or 'none' for r in rec['frigid_runs']})}", flush=True)
        if not args.no_fresh:
            t = time.perf_counter()
            fresh, got = fresh_recorded(rec, cache)
            print(f"   normal run {time.perf_counter() - t:.1f}s", flush=True)
            again = replay.replay_query(got, iso_cache=cache, display=True)
            d = diff(norm(again, drop_step_notes=(3,)), norm(fresh, drop_step_notes=(3,)))
            failed += bool(d)
            print("   replay vs normal run (same engine outputs): " + ("IDENTICAL" if not d else f"{len(d)} DIFFERENCES"))
            for line in d[:25]:
                print("     " + line)
            d2 = diff(norm(full, drop_step_notes=(3,)), norm(fresh, drop_step_notes=(3,)))
            print("   info: collect.py outputs vs this normal run: " + ("identical" if not d2 else f"{len(d2)} differences (engine outputs of two runs differ)"))
        quiet = replay.replay_query(rec, iso_cache=cache, step7=False)
        same = ranked(full) == ranked(quiet)
        failed += not same
        print(f"   ranking with vs without step 7: {'IDENTICAL' if same else 'DIFFERENT'} ({len(full['candidates'])} candidates)")
        st = {x["id"]: x for x in full["steps"]}
        ok7 = st[7]["status"] != "failed"
        failed += not ok7
        print(f"   full: step 7 status {st[7]['status']} (engine {((full.get('isomer') or {}).get('engine'))}, error {((full.get('isomer') or {}).get('error'))})"
              f"{'' if ok7 else '  <-- FAILED'}")
        cfg = replay.CONFIGS["leakfree_full"]
        lf = replay.replay_query(rec, layers=cfg["layers"], forced=cfg["forced"], step7=cfg["step7"], iso_cache=cache, display=True)
        ls = {x["id"]: x for x in lf["steps"]}
        from_library = [c["id"] for c in lf["candidates"] if any(sr["layer"] == "library" for sr in (c["sources"] or [])) or c["layer"] == "library"]
        ok_lf = (ls[3]["status"] == "skipped" and not lf["context"]["gate_passed"] and not from_library and ls[4]["status"] != "skipped"
                 and ls[7]["status"] != "failed")
        failed += not ok_lf
        print(f"   leakfree_full: step3 {ls[3]['status']}, gate_passed {lf['context']['gate_passed']}, library-sourced candidates {len(from_library)}, "
              f"steps 4/5/7 {ls[4]['status']}/{ls[5]['status']}/{ls[7]['status']}, {len(lf['candidates'])} candidates -> {'OK' if ok_lf else 'FAILED'}")
    print("ALL CHECKS PASSED" if not failed else f"{failed} CHECK(S) FAILED")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
