#!/usr/bin/env python3
"""Heavy pass: steps 1-5 of the ESI MS/MS Annotation workflow on all 500 CASMI 2022 challenges, run ONCE and unconditionally, raw outputs pickled.

  MZDESK_CODE_ROOT=/path/to/frozen/code python3 benchmark/casmi2022/mzdesk_ablation_20261030/collect.py [--workers 3] [--heavy-slots 2]

Per challenge, pipeline.run_query runs with its four engine entry points wrapped (formula_step.run_adduct_and_formula, library_step.run_library,
msfinder_step.run_msfinder_layer, frigid_step.run_frigid_layers). The wrappers record the real outputs and hand run_query a library gate that
"failed" and an MS-FINDER result that was "not good enough", so MS-FINDER (step 4) and FRIGID (step 5) run even when the workflow would skip
them; the FRIGID wrapper aborts the run after recording, so steps 6-8 are not run here. FRIGID is called for every challenge, but the engine
only supports positive adducts: negative ones are recorded as "skipped" by the engine itself. The library layer is untruncated
(library_step.MAX_CANDIDATES raised); replay.py truncates to the real 20 after any filtering. The truth is never read here.

Output (in --out, default this folder):
  store/challenge_NNN.pkl        everything replay needs: adduct+formula result, library layer + GateSummary, MS-FINDER layer + summary,
                                 FRIGID layer + targets, the arguments of the engine calls, the blind form
  store/challenge_NNN.meta.json  readable summary: seconds per step, code fingerprint, FRIGID device(s), layer statuses
  store/challenge_NNN.error.json an error (retried on the next start)
  manifest.json, collect.log (one line per case, 'ERROR' lines, 'ALL DONE n_ok n_err'), worker_logs/
Resumable: a challenge whose pickle and meta exist is skipped.
"""
from __future__ import annotations

import argparse
import copy
import json
import os
import pickle
import sys
import time
from pathlib import Path

import common
from common import HERE, SETTINGS, arm, busy_note, disarm, write_atomic, write_bytes_atomic

BIG = 10 ** 6


class _Abort(BaseException):
    """Raised by the FRIGID recording wrapper: steps 1-5 are done, steps 6-8 are not wanted in the collect pass."""


def _timed(fn, *a, **k):
    t = time.perf_counter()
    try:
        return ("ok", fn(*a, **k)), time.perf_counter() - t
    except Exception as exc:                  # noqa: BLE001 - run_query's own handlers must see the same failure on replay
        return ("exc", type(exc).__name__, str(exc)), time.perf_counter() - t


def collect_one(n):
    """One pass of steps 1-5 -> record dict (not yet pickled)."""
    from mzdesk.annotation import formula_step, frigid_step, library_step, msfinder_step, pipeline
    form, info = common.prepare(n)
    q = pipeline.parse_inputs(form)
    rec = {"n": n, "form": form, "input_info": info, "seconds": {}, "args": {}, "frigid_runs": []}
    real = dict(af=formula_step.run_adduct_and_formula, lib=library_step.run_library, msf=msfinder_step.run_msfinder_layer,
                fr=frigid_step.run_frigid_layers, predict=frigid_step.fe.predict)

    def reraise(rec_entry):
        _, name, msg = rec_entry
        raise type(name, (Exception,), {})(msg)       # run_query reports type(exc).__name__: keep the original class name

    def w_af(*a, **k):
        rec["af"], rec["seconds"]["adduct+formula"] = _timed(real["af"], *a, **k)
        if rec["af"][0] == "exc":
            reraise(rec["af"])
        return copy.deepcopy(rec["af"][1])

    def w_lib(*a, **k):
        old_cap = library_step.MAX_CANDIDATES
        library_step.MAX_CANDIDATES = BIG               # keep every judged compound; replay truncates to the real 20 after filtering
        try:
            rec["lib"], rec["seconds"]["library"] = _timed(real["lib"], *a, **k)
        finally:
            library_step.MAX_CANDIDATES = old_cap
        if rec["lib"][0] == "exc":
            reraise(rec["lib"])
        layer, summary = copy.deepcopy(rec["lib"][1])
        summary.passed = False                          # make run_query run steps 4-5 whatever the gate said
        return layer, summary

    def w_msf(precursor_mz, adduct, ms2, formulas, **k):
        rec["args"]["msfinder"] = (precursor_mz, adduct, ms2, [f.to_dict() for f in formulas], k)
        rec["msf"], rec["seconds"]["msfinder"] = _timed(real["msf"], precursor_mz, adduct, ms2, formulas, **k)
        if rec["msf"][0] == "exc":
            reraise(rec["msf"])
        layer, summary = copy.deepcopy(rec["msf"][1])
        summary["good_enough"] = False                  # make run_query run step 5 whatever MS-FINDER said
        return layer, summary

    def w_predict(*a, **k):                             # the FRIGID device actually used, per engine call
        t = time.perf_counter()
        try:
            res = real["predict"](*a, **k)
        except Exception as exc:                        # noqa: BLE001
            rec["frigid_runs"].append({"error": f"{type(exc).__name__}: {exc}", "seconds": time.perf_counter() - t})
            raise
        keep = ("formula", "adduct", "mode", "seed", "device", "device_requested", "device_reason", "seconds", "driver_seconds",
                "peak_footprint_mb", "peak_rss_mb", "n_samples", "n_valid", "n_formula_matched", "n_unique_matched", "length_model")
        rec["frigid_runs"].append({k_: res.get(k_) for k_ in keep})
        return res

    def w_fr(precursor_mz, targets, ms2, **k):
        rec["args"]["frigid"] = (precursor_mz, targets, ms2, k)
        rec["fr"], rec["seconds"]["frigid"] = _timed(real["fr"], precursor_mz, targets, ms2, **k)
        raise _Abort()

    pairs = ((formula_step, "run_adduct_and_formula", w_af), (library_step, "run_library", w_lib),
             (msfinder_step, "run_msfinder_layer", w_msf), (frigid_step, "run_frigid_layers", w_fr), (frigid_step.fe, "predict", w_predict))
    old = [(o, name, getattr(o, name)) for o, name, _ in pairs]
    t0 = time.perf_counter()
    try:
        for o, name, new in pairs:
            setattr(o, name, new)
        try:
            pipeline.run_query(q)
        except _Abort:
            pass
    finally:
        for o, name, was in old:
            setattr(o, name, was)
    rec["seconds"]["total"] = time.perf_counter() - t0
    missing = [k for k in ("af", "lib", "msf", "fr") if k not in rec]
    assert not missing, f"steps not reached: {missing}"
    return rec


def layer_status(entry):
    if entry[0] == "exc":
        return "failed", f"{entry[1]}: {entry[2]}"
    layer = entry[1][0] if isinstance(entry[1], tuple) else entry[1]
    return layer.status, layer.note


def meta_of(rec, attempts):
    st = {}
    for key, name in (("lib", "library"), ("msf", "msfinder"), ("fr", "frigid")):
        status, note = layer_status(rec[key])
        st[name] = {"status": status, "note": note[:300]}
    af = rec["af"]
    if af[0] == "ok":
        adduct_res, formula_res = af[1]
        st["adduct"] = {"status": adduct_res.status, "retained": [str(a) for a in adduct_res.retained_adducts]}
        st["formula"] = {"status": formula_res.status, "n_formulas": len(formula_res.formulas)}
    else:
        st["adduct"] = st["formula"] = {"status": "failed", "note": f"{af[1]}: {af[2]}"}
    lib = rec["lib"]
    gate = None
    if lib[0] == "ok":
        layer, summary = lib[1]
        gate = {"passed": summary.passed, "n_candidates": len(layer.candidates)}
    msf = rec["msf"]
    return {"challenge_id": rec["n"], "seconds": rec["seconds"], "code_fingerprint": common.STATE.get("fingerprint"),
            "code_root": str(common.CODE), "worker_pid": os.getpid(), "attempts": attempts, "layers": st, "library_gate": gate,
            "msfinder_good_enough": bool(msf[0] == "ok" and msf[1][1].get("good_enough")),
            "frigid_runs": rec["frigid_runs"], "frigid_devices": sorted({r.get("device") or "none" for r in rec["frigid_runs"]}),
            "input": rec["input_info"]}


def busy_layers(rec):
    return [name for key, name in (("msf", "msfinder"), ("fr", "frigid")) if layer_status(rec[key])[0] == "failed" and busy_note(layer_status(rec[key])[1])]


def work(n, out_dir, case_timeout, busy_retries=4):
    store = Path(out_dir) / "store"
    ok_pkl, ok_meta, err = store / f"challenge_{n:03d}.pkl", store / f"challenge_{n:03d}.meta.json", store / f"challenge_{n:03d}.error.json"
    t0 = time.perf_counter()
    try:
        attempts = 0
        while True:
            attempts += 1
            arm(case_timeout)
            try:
                rec = collect_one(n)
            finally:
                disarm()
            busy = busy_layers(rec)
            if not busy or attempts > busy_retries:
                break
            time.sleep(30 * attempts)               # an engine said busy: wait, then run the whole case again
        rec["fingerprint"], rec["code_root"], rec["attempts"] = common.STATE.get("fingerprint"), str(common.CODE), attempts
        rec["busy_after_retries"] = busy
        meta = meta_of(rec, attempts)
        meta["busy_layers_after_retries"] = busy
        write_bytes_atomic(ok_pkl, pickle.dumps(rec, protocol=pickle.HIGHEST_PROTOCOL))
        write_atomic(ok_meta, meta)                   # written last: it marks the case as complete
        err.unlink(missing_ok=True)
        failed = [k for k, v in meta["layers"].items() if v["status"] == "failed"]
        s = rec["seconds"]
        return {"id": n, "status": "ok", "seconds": time.perf_counter() - t0,
                "detail": f"af={s.get('adduct+formula', 0):.0f} lib={s.get('library', 0):.0f} msf={s.get('msfinder', 0):.0f} fr={s.get('frigid', 0):.0f}"
                          f" device={','.join(meta['frigid_devices'])}" + (f" failed={','.join(failed)}" if failed else "")
                          + (f" attempts={attempts}" if attempts > 1 else "")}
    except BaseException as exc:  # noqa: BLE001 - a case must never stop the run
        seconds = time.perf_counter() - t0
        write_atomic(err, common.error_payload(n, exc, seconds))
        if isinstance(exc, (KeyboardInterrupt, SystemExit)):
            raise
        return {"id": n, "status": "error", "seconds": seconds, "error": f"{type(exc).__name__}: {exc}"}


def check_manifest(out, ids, prepared):
    path = out / "manifest.json"
    new = {"seed": None, "sampling": "all 500 challenges", "ids": ids, "settings": SETTINGS,
           "inputs": {str(n): info for n, (_, info) in prepared.items()}, "code_root": str(common.CODE),
           "code_fingerprint": common.code_fingerprint(), "pipeline_sha256": common.sha(common.CODE / "mzdesk/annotation/pipeline.py")}
    if not path.exists():
        path.write_text(json.dumps(new, indent=2, ensure_ascii=False) + "\n")
        return
    old = json.loads(path.read_text())
    assert old["settings"] == SETTINGS, "Existing manifest has other settings; do not mix runs."
    assert old["code_fingerprint"] == new["code_fingerprint"], (
        f"The code root changed since the first start (fingerprint {old['code_fingerprint']} -> {new['code_fingerprint']}); "
        "a store must come from ONE code version. Use a new --out folder.")
    for k, v in new["inputs"].items():
        assert old["inputs"].get(k, v) == v, f"Input of challenge {k} changed since the manifest was written."
    old["inputs"].update(new["inputs"])
    old["ids"] = sorted(set(old["ids"]) | set(ids))
    path.write_text(json.dumps(old, indent=2, ensure_ascii=False) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=3, help="worker processes (default 3)")
    ap.add_argument("--heavy-slots", type=int, default=2, help="MS-FINDER/FRIGID jobs allowed at once across all workers (default 2)")
    ap.add_argument("--frigid-threads", type=int, default=4, help="CPU threads of one FRIGID run (MZDESK_FRIGID_THREADS; default 4)")
    ap.add_argument("--ids", default="1-500", help='challenges, e.g. "1-500" or "9,116,227" (default all)')
    ap.add_argument("--out", default=str(HERE), help="output directory (default: this folder)")
    ap.add_argument("--case-timeout", type=int, default=3600, help="seconds before a case is abandoned as an error (default 3600)")
    args = ap.parse_args()
    assert args.workers >= 1 and args.heavy_slots >= 1
    common.require_code_root()
    out = Path(args.out).resolve()
    (out / "store").mkdir(parents=True, exist_ok=True)
    ids = common.parse_ids(args.ids)
    os.environ["MZDESK_FRIGID_THREADS"] = str(args.frigid_threads)   # OMP_NUM_THREADS=1 would otherwise leak into the FRIGID driver
    todo = [n for n in ids if not ((out / "store" / f"challenge_{n:03d}.pkl").exists() and (out / "store" / f"challenge_{n:03d}.meta.json").exists())]
    print(f"{len(ids) - len(todo)} of {len(ids)} already collected; {len(todo)} to run  (code {common.CODE})", flush=True)
    prepared = {n: common.prepare(n) for n in ids}          # checks every blind input before anything starts
    check_manifest(out, ids, prepared)
    del prepared

    def dead(n, msg):
        write_atomic(out / "store" / f"challenge_{n:03d}.error.json", {"challenge_id": n, "error": msg, "elapsed_seconds": None})

    common.run_pool(kind="collect", todo=todo, work=work, out=out, workers=args.workers, heavy_slots=args.heavy_slots,
                    case_timeout=args.case_timeout, on_dead_twice=dead, header=f"frigid_threads={args.frigid_threads} code={common.CODE}")


if __name__ == "__main__":
    main()
