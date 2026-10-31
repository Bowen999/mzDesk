#!/usr/bin/env python3
"""Record BUDDY outputs (every BUDDY-supported adduct of the polarity) + the baseline run_adduct_and_formula result per CASMI challenge.
Usage: python3 record.py [--workers 4] [--ids 1-500]   (run from anywhere; unmodified current code => baseline)"""
import os
for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[k] = "1"
import argparse, concurrent.futures as cf, multiprocessing as mp, pickle, sys, time
from pathlib import Path
ROOT = Path("/path/to/mzDesk")
HERE = Path(__file__).resolve().parent
CACHE = HERE / "cache"

def load_run():
    sys.path.insert(0, str(HERE))
    import common                    # prepare(): the blind input shared with collect.py
    sys.path.insert(0, str(ROOT))    # current (unmodified) mzdesk code
    return common

def work(n):
    out = CACHE / f"c{n:03d}.pkl"
    if out.exists():
        return n, "cached", 0
    t0 = time.time()
    run = load_run()
    from mzdesk.annotation import formula_step as F, pipeline, msfinder_step
    from mzdesk import adducts
    form, _ = run.prepare(n)
    q = pipeline.parse_inputs(form)
    ms2 = msfinder_step.clean_peaks(q.ms2, pipeline.MAX_ENGINE_PEAKS)
    args = dict(adduct=q.adduct, tol_ppm=q.tol_ppm, tol_da=q.tol_da, halogen=q.halogen)
    seen = {}
    orig = F._Buddy._run
    def rec(self, adduct):
        r = orig(self, adduct); seen[adduct] = r; return r
    F._Buddy._run = rec
    a, f = F.run_adduct_and_formula(q.precursor_mz, q.positive, ms2, q.ms1 or None, **args)
    base = (a.to_dict(), f.to_dict())
    # all supported adducts of the polarity (BUDDY results are lru-cached for the ones already run)
    from mzdesk.annotation import adduct_step as A
    use_ms1, _ = F.ms1_usable(F.clean_peaks(q.ms1), q.precursor_mz, q.tol_ppm)
    b = F._Buddy(q.precursor_mz, q.positive, F.clean_peaks(ms2), F.clean_peaks(q.ms1) if use_ms1 else [], q.tol_ppm, q.halogen)
    for name in A.master_adducts(q.positive):
        b.run(name)
    seen.update(b._runs)
    rec_ = {"id": n, "precursor_mz": q.precursor_mz, "positive": q.positive, "ms2": ms2, "ms1": q.ms1 or None, "args": args,
            "runs": seen, "baseline": base, "seconds": time.time() - t0}
    tmp = out.with_suffix(".tmp")
    tmp.write_bytes(pickle.dumps(rec_)); tmp.rename(out)
    return n, "ok", time.time() - t0

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=4); ap.add_argument("--ids", default="1-500")
    a = ap.parse_args()
    lo, _, hi = a.ids.partition("-"); ids = list(range(int(lo), int(hi or lo) + 1))
    CACHE.mkdir(exist_ok=True)
    ctx = mp.get_context("spawn")
    with cf.ProcessPoolExecutor(a.workers, mp_context=ctx) as ex:
        for n, st, sec in ex.map(work, ids):
            print(n, st, f"{sec:.1f}", flush=True)
