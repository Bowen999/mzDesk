#!/usr/bin/env python3
"""Single tools with their DEFAULT settings on the 500 CASMI 2022 challenges (no mzDesk adduct step, formula step, merge or re-rank).

  MZDESK_CODE_ROOT=/path/to/frozen/code python3 benchmark/casmi2022/mzdesk_ablation_20261030/standalone.py [--workers 3] [--heavy-slots 2]
  ... standalone.py --verify-frigid 16,54,98          # rerun FRIGID for cases whose output would be reused from the collect store and compare

Protocol (agreed with the owner)
  * Default adduct: [M+H]+ (positive) / [M-H]- (negative). A challenge whose TRUE adduct is another one is scored wrong by score.py without any run,
    so BUDDY, MS-FINDER and FRIGID are run only for challenges whose official adduct is the default one (the truth is read here ONLY to select
    those cases, via --all-adducts every tool runs on every challenge). The library search needs no adduct and runs on all 500.
  * library   precursor 10 ppm window, MassBank + GNPS, weighted spectral entropy >= 0.70 (the pool threshold of the workflow), >= min(4, n peaks)
              matched peaks, one hit per IK14 (best similarity), ranked by similarity. No adduct / formula / isotope gate check, no cap of 200 spectra.
              Every compound is stored; score.py also derives the leak-free list (drops IK14 == truth).
  * buddy     no run: the cached BUDDY outputs of buddy_record.py (every adduct of every challenge) -> the default adduct's formulas within 10 ppm
              in BUDDY's own order. Checked against the collect store (same candidates for the adducts both hold).
  * msfinder  msfinder_step.run_msfinder_layer with an EMPTY formula list: MS-FINDER finds its own formulas around the precursor (the step's
              own-formula mode), everything else as in mzDesk's MS-FINDER layer (10 ppm / 0.01 Da, structure databases of the shipped parameter file,
              TreeDepth 2 with the TreeDepth-1 rerun, time limits, top_n 10, MS1 isotope peaks given to MS-FINDER's formula search).
  * frigid    formula = BUDDY's #1 formula of [M+H]+ within 10 ppm (cache), run_frigid_layers with one target, top_n 5, quick mode, seed 42, cpu.
              Positive [M+H]+ only. A run of the collect store for exactly that formula + adduct is reused (identity checked by --verify-frigid).
Output: standalone/<tool>/challenge_NNN.json (+ .error.json), standalone.log, worker_logs/. Resumable: existing results are skipped.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import pickle
import time
from pathlib import Path

import common
from common import CASMI, HERE, arm, busy_note, disarm, write_atomic

CACHE = Path(os.environ.get("MZDESK_BUDDY_CACHE", ""))      # cNNN.pkl made by buddy_record.py
TOOLS = ("library", "buddy", "msfinder", "frigid")
DEFAULT_ADDUCT = {True: "[M+H]+", False: "[M-H]-"}
LIB_THRESHOLD = 0.70
MSF_TOP_N, FRIGID_TOP_N = 10, 5                             # pipeline.MSFINDER_TOP_N / FRIGID_TOP_N (checked against the frozen code in work())
_TRUTH_ADDUCT = {}


def truth_adduct(n):
    """The official adduct of challenge n. Used ONLY to decide which cases to run (a wrong default adduct is scored wrong without a run)."""
    if not _TRUTH_ADDUCT:
        with open(CASMI / "casmi2022_ground_truth.csv", newline="") as f:
            for r in csv.DictReader(f):
                _TRUTH_ADDUCT[int(r["challenge_id"])] = r["adduct_official"]
    return _TRUTH_ADDUCT[n]


def tool_dir(out, tool):
    d = Path(out) / "standalone" / tool
    d.mkdir(parents=True, exist_ok=True)
    return d


def cand_dict(c, ik14):
    return {"smiles": c.smiles, "ik14": ik14, "name": c.name, "formula": c.formula, "adduct": c.adduct, "score": c.layer_score,
            "score_norm": c.layer_score_norm, "mz_error_ppm": c.mz_error_ppm}


def layer_dict(layer, ik14_of):
    return {"status": layer.status, "note": layer.note, "seconds": layer.seconds,
            "candidates": [cand_dict(c, ik14_of(c)) for c in layer.candidates]}


# ---------------------------------------------------------------- tools
def run_library(q, n_hits_info):
    from mzdesk.annotation import library_step as L
    from mzdesk.annotation.adduct_step import clean_peaks
    from mzdesk.annotation.merge import has_structure, structure_keys
    peaks = clean_peaks(q.ms2)
    min_matched = min(L.GATE_MIN_MATCHED, len(peaks))
    old = L.SEARCH_TOP_N
    L.SEARCH_TOP_N = 10 ** 6                       # no cap of the judged spectra: the default search keeps every hit
    try:
        res, metric = L._search(peaks, q.libraries, q.precursor_mz, q.positive, q.tol_ppm, q.tol_da, L.PRIMARY_METRIC, LIB_THRESHOLD)
    finally:
        L.SEARCH_TOP_N = old
    best = {}
    n_low_peaks = 0
    for hit in res["hits"]:
        sim = hit["scores"].get(metric)
        if sim is None or sim < LIB_THRESHOLD:
            continue
        if hit["matches"] < min_matched:
            n_low_peaks += 1
            continue
        ik14 = structure_keys(hit.get("smiles"))[1] or (L._inchikey(hit) or "")[:14] or None
        key = ik14 or (hit.get("name") or "").strip().lower() or (hit["library"], hit["accession"])
        row = {"ik14": ik14, "name": hit.get("name") or "", "score": sim, "matches": hit["matches"], "library": hit["library"],
               "accession": hit["accession"], "has_structure": bool(has_structure(hit.get("smiles")) or ik14), "adduct": hit.get("adduct") or None,
               "ppm": hit["ppm"]}
        if key not in best or (sim, hit["matches"]) > (best[key]["score"], best[key]["matches"]):
            best[key] = row
    compounds = sorted(best.values(), key=lambda r: (-r["score"], -r["matches"]))
    return {"metric": metric, "min_matched_peaks": min_matched, "n_spectra_kept": res["n_kept"], "n_spectra_compared": res["n_candidates"],
            "n_dropped_few_peaks": n_low_peaks, "warnings": list(res["warnings"]), "compounds": compounds}


def buddy_view(n, q, ms2_engine, store_rec=None):
    """The cached BUDDY run of the default adduct, in BUDDY's order within the precursor tolerance; plus the consistency checks."""
    assert CACHE and CACHE.is_dir(), "Set MZDESK_BUDDY_CACHE to the cache folder of buddy_record.py (cNNN.pkl)"
    cache = pickle.loads((CACHE / f"c{n:03d}.pkl").read_bytes())
    adduct = DEFAULT_ADDUCT[q.positive]
    run = cache["runs"].get(adduct)
    checks = {"precursor_equal": cache["precursor_mz"] == q.precursor_mz, "polarity_equal": cache["positive"] == q.positive,
              "ms2_equal": [tuple(p) for p in cache["ms2"]] == [tuple(p) for p in ms2_engine]}
    if run is None:
        return {"adduct": adduct, "status": "missing", "formulas": [], "checks": checks}
    status, payload = run
    out = {"adduct": adduct, "status": status, "formulas": [], "checks": checks}
    if status != "ok":
        out["note"] = str(payload)[:300]
        return out
    out["formulas"] = [{"formula": c["formula"], "rank": c["rank"], "prob": c["prob"], "ppm": c["ppm"]} for c in payload["candidates"]
                       if c["ppm"] is not None and abs(c["ppm"]) <= q.tol_ppm]
    # the collect store holds BUDDY's list of the adducts mzDesk asked for: for the same adduct it must give the same formulas and ranks
    if store_rec is not None and store_rec["af"][0] == "ok":
        stored = {(f.formula, f.rank) for f in store_rec["af"][1][1].formulas if f.adduct == adduct}
        if stored:
            checks["store_formulas_equal"] = stored == {(f["formula"], f["rank"]) for f in out["formulas"]}
            out["checks"] = checks
    return out


def run_msfinder(q, ms2_engine):
    from mzdesk.annotation import msfinder_step
    from mzdesk.annotation.merge import structure_keys
    adduct = DEFAULT_ADDUCT[q.positive]
    attempts = 0
    while True:
        attempts += 1
        layer, summary = msfinder_step.run_msfinder_layer(q.precursor_mz, adduct, ms2_engine, [], tol_ppm=q.tol_ppm, tol_da=q.tol_da,
                                                          top_n=MSF_TOP_N, ms1_peaks=q.ms1 or None)
        if not (layer.status == "failed" and busy_note(layer.note)) or attempts > 4:
            break
        time.sleep(30 * attempts)
    ik14_of = lambda c: structure_keys(c.smiles)[1] or c.ik14       # noqa: E731 - merge's own key
    out = layer_dict(layer, ik14_of)
    out.update(adduct=adduct, attempts=attempts, mode=summary["mode"], retried=summary["retried"], warnings=summary["warnings"],
               formulas=summary["formulas"], engine_seconds=summary["engine_seconds"], peak_rss_mb=summary["peak_rss_mb"])
    return out


def frigid_from_store(store_rec, formula, adduct):
    """The candidates of the collect run for exactly this formula + adduct, or None. The store holds one merged layer; candidates carry formula and adduct.
    Reusable when the run was alone or the other targets have another formula (a structure belongs to one formula, so the merge dropped nothing of this run)."""
    from mzdesk.annotation.merge import structure_keys
    if store_rec["fr"][0] != "ok":
        return None
    targets = store_rec["args"]["frigid"][1]
    mine = [t for t in targets if t == (formula, adduct)]
    others = [t for t in targets if t != (formula, adduct)]
    if len(mine) != 1 or any(t[0] == formula for t in others):
        return None
    layer = store_rec["fr"][1][0]
    runs = store_rec.get("frigid_runs") or []
    if not runs or len(runs) != len(targets) or any("device" not in r or r["device"] != "cpu" for r in runs):
        return None                                            # a run that failed or ran elsewhere than cpu is not reused
    mine_c = [c for c in layer.candidates if c.formula == formula and c.adduct == adduct]
    if layer.status != "done":
        return None
    ik14_of = lambda c: c.ik14 or structure_keys(c.smiles)[1]  # noqa: E731
    return {"status": layer.status, "note": layer.note, "seconds": next((r["seconds"] for r in runs if r.get("formula") == formula), None),
            "candidates": [cand_dict(c, ik14_of(c)) for c in mine_c]}


def run_frigid(q, ms2_engine, formula):
    from mzdesk.annotation import frigid_step
    adduct = DEFAULT_ADDUCT[q.positive]
    runs = []
    real = frigid_step.fe.predict

    def spy(*a, **k):
        res = real(*a, **k)
        runs.append({k_: res.get(k_) for k_ in ("formula", "adduct", "seed", "device", "device_requested", "seconds", "n_samples", "n_valid", "peak_rss_mb")})
        return res
    frigid_step.fe.predict = spy
    attempts = 0
    try:
        while True:
            attempts += 1
            layer, _ = frigid_step.run_frigid_layers(q.precursor_mz, [(formula, adduct)], ms2_engine, top_n=FRIGID_TOP_N)
            if not (layer.status == "failed" and busy_note(layer.note)) or attempts > 4:
                break
            time.sleep(30 * attempts)
    finally:
        frigid_step.fe.predict = real
    out = layer_dict(layer, lambda c: c.ik14)
    out.update(attempts=attempts, runs=runs)
    return out


# ---------------------------------------------------------------- one challenge
def work(n, out_dir, case_timeout):
    out = Path(out_dir)
    t0 = time.perf_counter()
    done, errors, times = [], [], {}
    fp = common.STATE.get("fingerprint")
    verify = os.environ.get("MZDESK_STANDALONE_VERIFY") == "1"
    try:
        from mzdesk.annotation import frigid_step, msfinder_step, pipeline
        assert (pipeline.MSFINDER_TOP_N, pipeline.FRIGID_TOP_N) == (MSF_TOP_N, FRIGID_TOP_N), "top_n of the pipeline changed"
        form, _ = common.prepare(n)
        q = pipeline.parse_inputs(form)
        ms2_engine = msfinder_step.clean_peaks(q.ms2, pipeline.MAX_ENGINE_PEAKS)
        default = DEFAULT_ADDUCT[q.positive]
        run_adduct_tools = os.environ.get("MZDESK_STANDALONE_ALL") == "1" or truth_adduct(n) == default
        store_rec = None
        store_dir = Path(os.environ["MZDESK_STANDALONE_STORE"])
        if (store_dir / f"challenge_{n:03d}.pkl").exists():
            store_rec = pickle.loads((store_dir / f"challenge_{n:03d}.pkl").read_bytes())
        if store_rec is not None and store_rec.get("fingerprint") != fp:
            raise RuntimeError(f"collect store is from code {store_rec.get('fingerprint')}, this run uses {fp}")
        wanted = ["library"] + (["buddy", "msfinder"] + (["frigid"] if q.positive else []) if run_adduct_tools else [])
        wanted = [t for t in wanted if t in os.environ.get("MZDESK_STANDALONE_TOOLS", ",".join(TOOLS)).split(",")]
        buddy = None
        for tool in wanted:
            if verify and tool != "frigid":
                continue
            path = tool_dir(out, tool) / f"challenge_{n:03d}.json"
            if path.exists() and not verify:
                continue
            c0 = time.perf_counter()
            try:
                arm(case_timeout)
                try:
                    if tool == "library":
                        body = run_library(q, None)
                    elif tool == "buddy":
                        body = buddy_view(n, q, ms2_engine, store_rec)
                    elif tool == "msfinder":
                        body = run_msfinder(q, ms2_engine)
                    else:
                        if buddy is None:
                            buddy = buddy_view(n, q, ms2_engine, store_rec)
                        formulas = buddy["formulas"]
                        if not formulas:
                            body = {"status": "skipped", "note": "BUDDY gave no formula within the tolerance for the default adduct.", "candidates": [],
                                    "source": "none"}
                        else:
                            formula = formulas[0]["formula"]
                            reused = frigid_from_store(store_rec, formula, default) if store_rec is not None else None
                            if verify:
                                fresh = run_frigid(q, ms2_engine, formula)
                                same = None if reused is None else (
                                    [(c["ik14"], c["smiles"], c["score"]) for c in reused["candidates"]] ==
                                    [(c["ik14"], c["smiles"], c["score"]) for c in fresh["candidates"]])
                                body = {"formula": formula, "reusable": reused is not None, "identical": same,
                                        "stored": reused and [(c["ik14"], c["score"]) for c in reused["candidates"]],
                                        "fresh": [(c["ik14"], c["score"]) for c in fresh["candidates"]], "fresh_runs": fresh["runs"]}
                                write_atomic(tool_dir(out, "frigid_verify") / f"challenge_{n:03d}.json", {"challenge_id": n, **body})
                                done.append("verify")
                                continue
                            if reused is not None:
                                body = {**reused, "source": "collect store (same formula, adduct, settings)", "formula": formula}
                            else:
                                body = {**run_frigid(q, ms2_engine, formula), "source": "run", "formula": formula}
                finally:
                    disarm()
                if tool == "buddy":
                    buddy = body
                write_atomic(path, {"challenge_id": n, "tool": tool, "default_adduct": default, "code_fingerprint": fp,
                                    "elapsed_seconds": time.perf_counter() - c0, **body})
                path.with_suffix(".error.json").unlink(missing_ok=True)
                done.append(tool)
                times[tool] = time.perf_counter() - c0
            except BaseException as exc:  # noqa: BLE001 - one tool must not stop the others
                if isinstance(exc, (KeyboardInterrupt, SystemExit)):
                    raise
                write_atomic(path.with_suffix(".error.json"), common.error_payload(n, exc, time.perf_counter() - c0, tool=tool))
                errors.append(f"{tool}:{type(exc).__name__}")
        if errors:
            return {"id": n, "status": "error", "seconds": time.perf_counter() - t0, "error": "tools failed: " + ", ".join(errors)}
        detail = " ".join(f"{k}={v:.1f}" for k, v in times.items()) + ("" if run_adduct_tools else " wrong-default-adduct: library only")
        return {"id": n, "status": "ok", "seconds": time.perf_counter() - t0, "detail": detail.strip() or "nothing to do"}
    except BaseException as exc:  # noqa: BLE001
        if isinstance(exc, (KeyboardInterrupt, SystemExit)):
            raise
        return {"id": n, "status": "error", "seconds": time.perf_counter() - t0, "error": f"{type(exc).__name__}: {exc}"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--heavy-slots", type=int, default=2)
    ap.add_argument("--frigid-threads", type=int, default=4)
    ap.add_argument("--ids", default="1-500")
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--store", default=None, help="collect.py output directory holding store/ (default: --out)")
    ap.add_argument("--buddy-cache", default=os.environ.get("MZDESK_BUDDY_CACHE"), help="folder with cNNN.pkl from buddy_record.py (or MZDESK_BUDDY_CACHE)")
    ap.add_argument("--all-adducts", action="store_true", help="run BUDDY / MS-FINDER / FRIGID on every challenge, not only those whose official adduct is the default one")
    ap.add_argument("--verify-frigid", default=None, help='ids: rerun FRIGID where the store output would be reused and compare (writes standalone/frigid_verify/)')
    ap.add_argument("--tools", default=",".join(TOOLS), help="comma-separated subset of: " + ",".join(TOOLS))
    ap.add_argument("--case-timeout", type=int, default=3600)
    args = ap.parse_args()
    common.require_code_root()
    assert args.buddy_cache, "--buddy-cache (or MZDESK_BUDDY_CACHE) is required"
    os.environ["MZDESK_BUDDY_CACHE"] = str(Path(args.buddy_cache).resolve())
    os.environ["MZDESK_FRIGID_DEVICE"] = "cpu"                                 # as the collect run: no GPU / unified-memory swings
    os.environ["MZDESK_FRIGID_THREADS"] = str(args.frigid_threads)
    out = Path(args.out).resolve()
    store = (Path(args.store).resolve() if args.store else out) / "store"
    os.environ["MZDESK_STANDALONE_STORE"] = str(store)
    os.environ["MZDESK_STANDALONE_ALL"] = "1" if args.all_adducts else "0"
    assert set(args.tools.split(",")) <= set(TOOLS), args.tools
    os.environ["MZDESK_STANDALONE_TOOLS"] = args.tools
    os.environ["MZDESK_STANDALONE_VERIFY"] = "1" if args.verify_frigid else "0"
    ids = common.parse_ids(args.verify_frigid or args.ids)
    for t in TOOLS:
        tool_dir(out, t)
    if args.verify_frigid:
        todo, kind = ids, "standalone"
    else:
        def pending(n):
            need = ["library"] + ((["buddy", "msfinder"] + (["frigid"] if truth_adduct(n) in ("[M+H]+",) else [])) if args.all_adducts or truth_adduct(n) in DEFAULT_ADDUCT.values() else [])
            return [t for t in need if t in args.tools.split(",") and not (out / "standalone" / t / f"challenge_{n:03d}.json").exists()]
        todo = [n for n in ids if pending(n)]
    print(f"{len(todo)} of {len(ids)} challenges to run  (code {common.CODE})", flush=True)

    def dead(n, msg):
        write_atomic(out / "standalone" / "library" / f"challenge_{n:03d}.error.json", {"challenge_id": n, "error": msg, "elapsed_seconds": None})

    common.run_pool(kind="standalone", todo=todo, work=work, out=out, workers=args.workers, heavy_slots=args.heavy_slots,
                    case_timeout=args.case_timeout, on_dead_twice=dead,
                    header=f"all_adducts={args.all_adducts} verify={bool(args.verify_frigid)} code={common.CODE}")


if __name__ == "__main__":
    main()
