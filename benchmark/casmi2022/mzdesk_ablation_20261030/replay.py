#!/usr/bin/env python3
"""Offline replays of steps 6-8 of the ESI MS/MS Annotation workflow on the outputs stored by collect.py.

  MZDESK_CODE_ROOT=/path/to/the/SAME/frozen/code python3 benchmark/casmi2022/mzdesk_ablation_20261030/replay.py [--workers 3]

A replay calls the workflow's OWN pipeline.run_query with the four engine entry points (adduct+formula, library, MS-FINDER, FRIGID) returning the
stored outputs, so every glue decision (retained adduct, formula list, library gate, unverifiable-hit rule, MS-FINDER "good enough" skip, FRIGID
targets, warnings, merge, re-rank, isomer pool, answer) is the workflow's, not a re-implementation. Step 7 (a second MS-FINDER run on the top-1 and its
same-formula candidates) is real engine work: it is cached by its exact arguments (step7_cache/), shared by all configurations.

Configurations (replay/<name>/challenge_NNN.json):
  full             the workflow as it runs: gate / skip decisions applied by run_query, steps 6-8 incl. step 7
  leakfree_full    as full, but the library is not searched at all: step 3 returns nothing (skipped layer, a gate summary that did not pass), so the
                   workflow goes straight to de novo annotation (steps 4-5 by its own logic: MS-FINDER, FRIGID only when MS-FINDER is not good enough),
                   then merge / re-rank, step 7, answer. No truth is read.
  ablations        ranking only; step 7 is replaced by a stub (verdict "insufficient", step "skipped": the isomer step never changes a rank or a level,
                   see verify.py), every listed layer is used whatever the gate / "good enough" decided (forced), the other layers return nothing:
    all_layers         library + MS-FINDER + FRIGID merged (the reference of the leave-one-out ablations)
    library_only, msfinder_only, frigid_only (positive only: the engine skips negative adducts), no_library, no_msfinder, no_frigid
    buddy_only         formula level: the step-2 formula list, no structure (score.py reports the rank of the truth formula)
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import pickle
import sys
import time
from contextlib import contextmanager
from pathlib import Path

import common
from common import HERE, arm, busy_note, disarm, write_atomic, write_bytes_atomic

ALL = ("library", "msfinder", "frigid")


def _cfg(layers=ALL, forced=False, step7=False, kind="structure"):
    return {"layers": layers, "forced": forced, "step7": step7, "kind": kind}


CONFIGS = {
    "full": _cfg(step7=True),
    "leakfree_full": _cfg(("msfinder", "frigid"), step7=True),      # library not searched: the gate cannot pass, steps 4-5 follow the workflow's own logic
    "all_layers": _cfg(forced=True),
    "library_only": _cfg(("library",), forced=True),
    "msfinder_only": _cfg(("msfinder",), forced=True),
    "frigid_only": _cfg(("frigid",), forced=True),
    "no_library": _cfg(("msfinder", "frigid"), forced=True),
    "no_msfinder": _cfg(("library", "frigid"), forced=True),
    "no_frigid": _cfg(("library", "msfinder"), forced=True),
    "buddy_only": _cfg(layers=(), kind="formula"),
}

# ---------------------------------------------------------------- patching
@contextmanager
def patched(*pairs):
    old = [(obj, name, getattr(obj, name)) for obj, name, _ in pairs]
    try:
        for obj, name, new in pairs:
            setattr(obj, name, new)
        yield
    finally:
        for obj, name, was in old:
            setattr(obj, name, was)


def _no_isomers(*a, **k):
    """Step 7 left out: the engine's own 'cannot do it' answer, so run_query records a skipped step and moves on."""
    return {"verdict": "insufficient", "top_id": None, "group": [], "isomers": [], "challengers": [], "candidates": [], "peaks_considered": 0,
            "engine": "not run", "notes": ["Step 7 was left out of this replay (ranking-only ablation)."], "error": None, "log": "", "seconds": 0.0,
            "settings": {}}


class IsomerCache:
    """Step 7 runs MS-FINDER again for the top-1 and its same-formula candidates, so it cannot be 'offline' for a new top-1 / pool: cache by exact arguments."""

    def __init__(self, directory):
        from mzdesk.exclusivity import isomers
        self.dir = Path(directory)
        self.dir.mkdir(exist_ok=True)
        self.hits = self.misses = 0
        self.miss_seconds = 0.0
        self.real = isomers.discriminate

    def __call__(self, *a, **k):
        key = hashlib.sha256(json.dumps([a, k], sort_keys=True, default=repr).encode()).hexdigest()[:24]
        path = self.dir / f"{key}.pkl"
        if path.exists():
            self.hits += 1
            return copy.deepcopy(pickle.loads(path.read_bytes()))
        t = time.perf_counter()
        res = self.real(*a, **k)
        self.misses += 1
        self.miss_seconds += time.perf_counter() - t
        if not res.get("error"):                       # never cache an engine failure (busy / killed)
            write_bytes_atomic(path, pickle.dumps(res))
        return copy.deepcopy(res)


def _make_exc(entry):
    _, name, msg = entry
    return type(name, (Exception,), {})(msg)


def replay_query(rec, *, layers=ALL, forced=False, step7=True, iso_cache=None, display=False, strict=True, info=None):
    """pipeline.run_query on the stored outputs of `rec` (see the module docstring). Returns the raw result of run_query."""
    from mzdesk.annotation import formula_step, frigid_step, library_step, msfinder_step, pipeline, visualization
    from mzdesk.annotation.schema import LAYER_FRIGID, LAYER_LIBRARY, LAYER_MSFINDER, STATUS_SKIPPED, LayerResult
    from mzdesk.exclusivity import isomers
    from mzdesk.ms2 import fragment_view

    info = info if info is not None else {}
    q = pipeline.parse_inputs(rec["form"])

    def empty(layer):
        return LayerResult(layer, STATUS_SKIPPED, "Layer left out of this replay.", 0.0, [])

    def w_af(*a, **k):
        if rec["af"][0] == "exc":
            raise _make_exc(rec["af"])
        return copy.deepcopy(rec["af"][1])

    def w_lib(*a, **k):
        if "library" not in layers:
            return empty(LAYER_LIBRARY), library_step.GateSummary()
        if rec["lib"][0] == "exc":
            raise _make_exc(rec["lib"])
        layer, summary = copy.deepcopy(rec["lib"][1])
        cands = list(layer.candidates)
        cands = cands[:library_step.MAX_CANDIDATES]
        layer.candidates = cands
        passing = [c for c in cands if c.evidence["library"]["gate_passed"]]
        summary.passed, summary.best_id = bool(passing), (passing[0].id if passing else None)
        summary.top_id = cands[0].id if cands else None
        if forced:
            summary.passed = False
        return layer, summary

    def w_msf(precursor_mz, adduct, ms2, formulas, **k):
        if strict:                                             # steps 4-5 depend on steps 1-2 only: same arguments in every configuration
            assert (precursor_mz, adduct, ms2, [f.to_dict() for f in formulas], k) == rec["args"]["msfinder"], "MS-FINDER arguments drifted"
        if "msfinder" not in layers:
            return empty(LAYER_MSFINDER), {"good_enough": False, "warnings": []}
        if rec["msf"][0] == "exc":
            raise _make_exc(rec["msf"])
        layer, summary = copy.deepcopy(rec["msf"][1])
        if forced:
            summary["good_enough"] = False
        return layer, summary

    def w_fr(precursor_mz, targets, ms2, **k):
        if strict:
            assert (precursor_mz, targets, ms2, k) == rec["args"]["frigid"], "FRIGID arguments drifted"
        if "frigid" not in layers:
            return empty(LAYER_FRIGID), []
        if rec["fr"][0] == "exc":
            raise _make_exc(rec["fr"])
        return copy.deepcopy(rec["fr"][1])

    pairs = [(formula_step, "run_adduct_and_formula", w_af), (library_step, "run_library", w_lib), (msfinder_step, "run_msfinder_layer", w_msf),
             (frigid_step, "run_frigid_layers", w_fr), (isomers, "discriminate", iso_cache if step7 else _no_isomers)]
    if not display:                                            # display-only work (figures, fragment-view query): never affects a rank
        pairs += [(visualization, "spectrum_figures", lambda *a, **k: None), (fragment_view, "query", lambda *a, **k: None)]
    with patched(*pairs):
        return pipeline.run_query(q)


# ---------------------------------------------------------------- compact output
def compact(res):
    """What score.py needs from a run_query result, without the figures, the engine logs and the long evidence blocks."""
    cands = []
    for c in res["candidates"]:
        lib = (c.get("evidence") or {}).get("library") or {}
        cands.append({"id": c["id"], "rank": c["rank"], "layer": c["layer"], "level": c["level"], "smiles": c["smiles"], "ik14": c["ik14"],
                      "name": c["name"], "formula": c["formula"], "adduct": c["adduct"], "mz_error_ppm": c["mz_error_ppm"],
                      "layer_score": c["layer_score"], "layer_score_norm": c["layer_score_norm"], "hypothesis": c["hypothesis"],
                      "sources": c["sources"], "flags": c["flags"], "detect": c["detect"], "base_score": c["base_score"],
                      "final_score": c["final_score"],
                      "library": ({k: lib.get(k) for k in ("similarity", "gate_passed", "library", "accession", "has_structure")} if lib else None)})
    iso = res.get("isomer") or {}
    return {"answer": res["answer"], "context": res["context"], "warnings": res["warnings"], "candidates": cands,
            "steps": [{"id": s["id"], "status": s["status"], "seconds": s["seconds"], "note": s["note"][:300]} for s in res["steps"]],
            "isomer": ({k: iso.get(k) for k in ("verdict", "top_id", "group", "error", "seconds", "notes", "engine")} if iso else None)}


def formula_view(rec):
    """buddy_only: the step-2 formula list as the workflow hands it on (no structure, no ranking beyond BUDDY's)."""
    if rec["af"][0] == "exc":
        return {"error": f"{rec['af'][1]}: {rec['af'][2]}", "formulas": [], "adducts": [], "retained_adducts": []}
    adduct_res, formula_res = rec["af"][1]
    return {"formulas": [f.to_dict() for f in formula_res.formulas], "status": formula_res.status, "note": formula_res.note,
            "adducts": [str(a) for a in formula_res.adducts], "retained_adducts": [str(a) for a in adduct_res.retained_adducts]}


# ---------------------------------------------------------------- one challenge
def load_store(store, n):
    path = Path(store) / f"challenge_{n:03d}.pkl"
    assert path.exists() and (Path(store) / f"challenge_{n:03d}.meta.json").exists(), f"challenge {n}: no complete collect result in {store}"
    return pickle.loads(path.read_bytes())


def run_config(rec, name, cfg, iso_cache, case_timeout, step7_retries=3):
    """-> (payload fields, seconds)."""
    if cfg["kind"] == "formula":
        return {"result": formula_view(rec)}, 0.0
    for attempt in range(1, step7_retries + 2):
        info = {}
        h0, m0, s0 = iso_cache.hits, iso_cache.misses, iso_cache.miss_seconds
        t = time.perf_counter()
        arm(case_timeout)
        try:
            res = replay_query(rec, layers=cfg["layers"], forced=cfg["forced"], step7=cfg["step7"], iso_cache=iso_cache, info=info)
        finally:
            disarm()
        seconds = time.perf_counter() - t
        step7 = next(s for s in res["steps"] if s["id"] == 7)
        if cfg["step7"] and step7["status"] == "failed" and busy_note(step7["note"]) and attempt <= step7_retries:
            time.sleep(30 * attempt)
            continue
        break
    out = {"result": compact(res), "step7_cache": {"hits": iso_cache.hits - h0, "misses": iso_cache.misses - m0,
                                                    "miss_seconds": iso_cache.miss_seconds - s0}, "attempts": attempt}
    return out, seconds


def library_empty(rec):
    """The stored library step ran without error and returned no candidate at all."""
    return rec["lib"][0] == "ok" and not rec["lib"][1][0].candidates


def work(n, out_dir, case_timeout):
    out = Path(out_dir)
    t0 = time.perf_counter()
    selected = os.environ["MZDESK_REPLAY_CONFIGS"].split(",")
    store = os.environ["MZDESK_REPLAY_STORE"]
    allow_mismatch = os.environ.get("MZDESK_REPLAY_ALLOW_CODE_MISMATCH") == "1"
    done, errors, times = [], [], {}
    try:
        rec = load_store(store, n)
        fp = common.STATE.get("fingerprint")
        if rec.get("fingerprint") != fp and not allow_mismatch:
            raise RuntimeError(f"stored outputs come from code {rec.get('fingerprint')}, this replay runs code {fp}; "
                               "use the same MZDESK_CODE_ROOT as collect.py (or --allow-code-mismatch)")
        iso_cache = IsomerCache(out / "step7_cache")
        for name in selected:
            cfg = CONFIGS[name]
            folder = out / "replay" / name
            ok_path, err_path = folder / f"challenge_{n:03d}.json", folder / f"challenge_{n:03d}.error.json"
            if ok_path.exists():
                continue
            c0 = time.perf_counter()
            full_path = out / "replay" / "full" / f"challenge_{n:03d}.json"
            if name == "leakfree_full" and not os.environ.get("MZDESK_REPLAY_NO_COPY") and library_empty(rec) and full_path.exists():
                # No library candidate was stored: leaving the library out changes nothing (the gate fails either way, steps 4-7 see the same inputs),
                # so the result IS the full result. Proven on every challenge of this kind that was replayed both ways (see README).
                payload = json.loads(full_path.read_text())
                payload.update(config=name, options={k: list(v) if isinstance(v, tuple) else v for k, v in cfg.items()}, copied_from_full=True,
                               step7_cache={"hits": 0, "misses": 0, "miss_seconds": 0.0}, elapsed_seconds=0.0)
                write_atomic(ok_path, payload)
                done.append(name)
                times[name] = 0.0
                continue
            try:
                fields, seconds = run_config(rec, name, cfg, iso_cache, case_timeout)
                payload = {"challenge_id": n, "config": name, "options": {k: list(v) if isinstance(v, tuple) else v for k, v in cfg.items()},
                           "elapsed_seconds": time.perf_counter() - c0, "code_fingerprint": fp, "store_fingerprint": rec.get("fingerprint"), **fields}
                write_atomic(ok_path, payload)
                err_path.unlink(missing_ok=True)
                done.append(name)
                times[name] = payload["elapsed_seconds"]
            except BaseException as exc:  # noqa: BLE001 - one configuration must not stop the others
                if isinstance(exc, (KeyboardInterrupt, SystemExit)):
                    raise
                write_atomic(err_path, common.error_payload(n, exc, time.perf_counter() - c0, config=name))
                errors.append(f"{name}:{type(exc).__name__}")
        detail = " ".join(f"{k}={v:.1f}" for k, v in times.items()) + (f" iso(h/m)={iso_cache.hits}/{iso_cache.misses}" if done else "")
        if errors:
            return {"id": n, "status": "error", "seconds": time.perf_counter() - t0, "error": "configurations failed: " + ", ".join(errors)}
        return {"id": n, "status": "ok", "seconds": time.perf_counter() - t0, "detail": detail or "nothing to do"}
    except BaseException as exc:  # noqa: BLE001
        if isinstance(exc, (KeyboardInterrupt, SystemExit)):
            raise
        return {"id": n, "status": "error", "seconds": time.perf_counter() - t0, "error": f"{type(exc).__name__}: {exc}"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=3, help="worker processes (default 3)")
    ap.add_argument("--heavy-slots", type=int, default=2, help="step-7 MS-FINDER jobs allowed at once across all workers (default 2)")
    ap.add_argument("--ids", default="1-500", help='challenges, e.g. "1-500" or "9,116,227" (default all)')
    ap.add_argument("--out", default=str(HERE), help="output directory (default: this folder)")
    ap.add_argument("--store", default=None, help="collect.py output directory holding store/ (default: --out)")
    ap.add_argument("--configs", default=",".join(CONFIGS), help="comma-separated configurations (default all): " + ",".join(CONFIGS))
    ap.add_argument("--case-timeout", type=int, default=1800, help="seconds before one configuration of one challenge is abandoned (default 1800)")
    ap.add_argument("--allow-code-mismatch", action="store_true", help="replay stored outputs made by another code version (not recommended)")
    args = ap.parse_args()
    common.require_code_root()
    out = Path(args.out).resolve()
    store = (Path(args.store).resolve() if args.store else out) / "store"
    configs = args.configs.split(",")
    unknown = [c for c in configs if c not in CONFIGS]
    assert not unknown, f"unknown configurations {unknown}"
    ids = common.parse_ids(args.ids)
    (out / "replay").mkdir(parents=True, exist_ok=True)
    for c in configs:
        (out / "replay" / c).mkdir(exist_ok=True)
    os.environ["MZDESK_REPLAY_CONFIGS"] = ",".join(configs)
    os.environ["MZDESK_REPLAY_STORE"] = str(store)
    os.environ["MZDESK_REPLAY_ALLOW_CODE_MISMATCH"] = "1" if args.allow_code_mismatch else "0"
    missing = [n for n in ids if not (store / f"challenge_{n:03d}.meta.json").exists()]
    assert not missing, f"{len(missing)} challenges have no complete collect result in {store} (first: {missing[:10]}); finish collect.py first."

    def pending(n):
        return [c for c in configs if not (out / "replay" / c / f"challenge_{n:03d}.json").exists()]
    todo = [n for n in ids if pending(n)]
    print(f"{len(ids) - len(todo)} of {len(ids)} challenges complete for {len(configs)} configurations; {len(todo)} to run  (code {common.CODE})", flush=True)

    def dead(n, msg):
        write_atomic(out / "replay" / configs[0] / f"challenge_{n:03d}.error.json", {"challenge_id": n, "error": msg, "elapsed_seconds": None})

    common.run_pool(kind="replay", todo=todo, work=work, out=out, workers=args.workers, heavy_slots=args.heavy_slots,
                    case_timeout=args.case_timeout, on_dead_twice=dead, header=f"configs={','.join(configs)} code={common.CODE}")


if __name__ == "__main__":
    main()
