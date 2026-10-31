#!/usr/bin/env python3
"""Shared pieces of collect.py and replay.py: frozen code root, blind input, file-based heavy lock, parallel resumable pool.

The mzdesk package is imported ONLY from the directory named by the environment variable MZDESK_CODE_ROOT (a `git archive` copy of
one commit), never from the live working tree, which other sessions edit. collect.py and replay.py refuse to start without it.
"""
from __future__ import annotations

import os

for _k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_k] = "1"            # one thread per worker (spawned workers re-import this module, so this runs before numpy there too)

import concurrent.futures as cf
import fcntl
import hashlib
import json
import multiprocessing as mp
import signal
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]                                   # repository root (used for benchmark inputs only, never for mzdesk code)
CASMI = ROOT / "benchmark/casmi2022"
CODE = Path(os.environ["MZDESK_CODE_ROOT"]).resolve() if os.environ.get("MZDESK_CODE_ROOT") else None
if CODE is not None:
    sys.path.insert(0, str(CODE))
sys.path.append(str(ROOT / "benchmark/tools/jm"))       # prep.py / jm_common.py (benchmark helpers; do not import mzdesk)

SETTINGS = {                          # mzDesk settings of every run
    "libraries": ["MassBank", "GNPS"],
    "tol_ppm": 10.0,
    "tol_da": 0.01,
    "omics": "Metabolomics",
    "halogen": True,
    "local_isomers": False,
    "additive": "Formic acid",
    "alkali": "Unknown",
    "ph": 2.7,
}
BUSY_MARKERS = ("is busy", "busy with another", "is waiting for", "Try again in")
STATE = {}


# ------------------------------------------------------------------ code root
def require_code_root():
    """Stop unless MZDESK_CODE_ROOT is set and mzdesk really comes from there."""
    assert CODE is not None, ("Set MZDESK_CODE_ROOT to a frozen copy of the repository (git archive of one commit); "
                              "the live tree is never imported because other sessions edit it.")
    assert (CODE / "mzdesk/annotation/pipeline.py").exists(), f"{CODE} holds no mzdesk/annotation/pipeline.py"
    import mzdesk
    assert Path(mzdesk.__file__).resolve().is_relative_to(CODE), f"mzdesk was imported from {mzdesk.__file__}, not from {CODE}"


def code_fingerprint():
    """Hash of every mzdesk/**/*.py of the code root."""
    h = hashlib.sha256()
    for p in sorted((CODE / "mzdesk").rglob("*.py")):
        h.update(str(p.relative_to(CODE)).encode())
        h.update(p.read_bytes())
    return h.hexdigest()[:16]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ------------------------------------------------------------------ blind input
def prepare(number):
    """Blind form from the SIRIUS-extracted .ms, cross-checked against the JESTR spectrum.

    Accepts challenges with only 1-2 non-empty collision-energy blocks (12 of the 500)."""
    import numpy as np
    from jm_common import report  # noqa: F401 - prep imports it from this directory
    from prep import merge_peaks, parse_ms
    from mzdesk.annotation import pipeline
    name = f"challenge_{number:03d}"
    source = CASMI / "sirius_run/ms" / f"{name}.ms"
    data = parse_ms(source)
    ion = data["meta"]["ionization"]
    assert ion in ("[M+?]+", "[M+?]-"), (name, ion)
    ms1 = [p for label, array in data["blocks"] if label == "ms1" for p in array]
    blocks = [array for label, array in data["blocks"] if label.startswith("collision") and len(array)]
    assert blocks and ms1, name          # at least one non-empty collision-energy block and an MS1 list
    merged = merge_peaks(np.vstack(blocks), 10.0)
    jestr = CASMI / "jm_run/inputs/spectra" / f"{name}.tsv"
    saved = np.loadtxt(jestr, delimiter="\t", skiprows=1, ndmin=2)
    assert merged.shape == saved.shape and np.allclose(merged, saved, rtol=0, atol=0.000051), name
    peaks = lambda rows: "\n".join(f"{m:.6f} {i:.4f}" for m, i in rows)
    form = {
        "precursor": {"precursor_mz": float(data["meta"]["parentmass"]),
                      "ion_mode": "Positive" if ion.endswith("+") else "Negative", "adduct": "Auto"},
        "ms1": peaks(ms1), "ms2": {"source": "Paste", "value": peaks(merged)}, **SETTINGS,
    }
    query = pipeline.parse_inputs(form)
    assert len(query.ms2) == len(saved) and len(query.ms1) == len(ms1)
    return form, {"ms_file": str(source.relative_to(ROOT)), "ms_sha256": sha(source),
                  "jestr_spectrum_sha256": sha(jestr), "ms1_peaks": len(ms1), "ms2_peaks": len(merged),
                  "collision_blocks": len(blocks)}


# ------------------------------------------------------------------ cross-process heavy lock
class SlotLock:
    """Drop-in for heavy_lock.HeavyLock: up to ``slots`` heavy jobs (MS-FINDER / FRIGID) at once across all worker processes.

    acquire() ignores the engine's 20 s timeout and waits (polling) up to ``max_wait`` seconds, so the engines never answer "busy"."""

    def __init__(self, directory, slots, max_wait=6 * 3600):
        self.dir, self.slots, self.max_wait = Path(directory), slots, max_wait
        self.dir.mkdir(parents=True, exist_ok=True)
        self.holder = None
        self._fd = None

    def acquire(self, blocking=True, timeout=-1, holder=None):
        deadline = time.monotonic() + self.max_wait
        while True:
            for i in range(self.slots):
                fd = os.open(self.dir / f"slot{i}.lock", os.O_CREAT | os.O_RDWR, 0o644)
                try:
                    fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except OSError:
                    os.close(fd)
                    continue
                self._fd, self.holder = fd, holder
                return True
            if time.monotonic() > deadline:
                return False
            time.sleep(0.5)

    def release(self):
        self.holder = None
        fd, self._fd = self._fd, None
        if fd is not None:
            fcntl.flock(fd, fcntl.LOCK_UN)
            os.close(fd)

    def locked(self):
        return self._fd is not None

    def __enter__(self):
        self.acquire()
        return self

    def __exit__(self, *exc):
        self.release()


def init_worker(out_dir, heavy_slots, kind):
    """Per-process setup: private log file, file-based heavy lock in place of the in-process one, code fingerprint."""
    out = Path(out_dir)
    logs = out / "worker_logs"
    logs.mkdir(exist_ok=True)
    fd = os.open(logs / f"{kind}_worker_{os.getpid()}.log", os.O_CREAT | os.O_WRONLY | os.O_APPEND, 0o644)
    os.dup2(fd, 1)
    os.dup2(fd, 2)                    # msbuddy / tqdm noise stays out of run.log
    require_code_root()
    from mzdesk.shared import heavy_lock
    from mzdesk.ms2 import msfinder_engine
    lock = SlotLock(out / ".locks", heavy_slots)
    heavy_lock.LOCK = lock
    msfinder_engine._RUN_LOCK = lock  # bound at import time, so patch it too
    if kind in ("collect", "standalone"):
        from mzdesk.ms2 import frigid_engine
        assert frigid_engine.heavy_lock.LOCK is lock
    STATE["fingerprint"] = code_fingerprint()


class CaseTimeout(BaseException):           # BaseException: the workflow's own 'except Exception' must not swallow it
    pass


def _alarm(signum, frame):
    raise CaseTimeout("case exceeded the time limit")


def arm(seconds):
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(int(seconds))


def disarm():
    signal.alarm(0)


def busy_note(note):
    return any(m in str(note) for m in BUSY_MARKERS)


def write_atomic(path, payload):
    path = Path(path)
    tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
    tmp.write_text(json.dumps(payload, indent=1, ensure_ascii=False, allow_nan=False) + "\n")
    tmp.rename(path)


def write_bytes_atomic(path, data):
    path = Path(path)
    tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
    tmp.write_bytes(data)
    tmp.rename(path)


def error_payload(n, exc, seconds, **extra):
    return {"challenge_id": n, "error": f"{type(exc).__name__}: {exc}", "traceback": traceback.format_exc(),
            "elapsed_seconds": seconds, "code_fingerprint": STATE.get("fingerprint"), **extra}


# ------------------------------------------------------------------ parent: parallel, resumable pool
def stamp():
    return datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


def log_line(log, text):
    with open(log, "a") as f:
        f.write(text + "\n")


def parse_ids(text):
    ids = []
    for part in text.split(","):
        a, _, b = part.partition("-")
        ids += range(int(a), int(b or a) + 1)
    assert ids and all(1 <= i <= 500 for i in ids), "ids must lie in 1..500"
    return sorted(set(ids))


def guarded(work, n, out_dir, case_timeout, marker_dir):
    """Runs in the worker: an in-flight marker lets the parent tell which cases a killed worker was holding."""
    marker = Path(marker_dir) / f".inflight_{n:03d}"
    marker.write_text(str(os.getpid()))
    try:
        return work(n, out_dir, case_timeout)
    finally:
        marker.unlink(missing_ok=True)


def run_pool(*, kind, todo, work, out, workers, heavy_slots, case_timeout, on_dead_twice, header=""):
    """Run work(n, out_dir, case_timeout) -> {"id","status": ok|error,"seconds","detail"|"error"} for every n in todo, in worker processes.

    Same behaviour as full500/run.py: spawn context, max_tasks_per_child=None (Python 3.12 ProcessPoolExecutor hangs when it replaces workers),
    a dead pool is restarted and a case that kills a worker twice becomes an error (on_dead_twice(n, message) writes its error file).
    run.log gets one line per case, 'ERROR ...' lines for errors, 'WALL_SECONDS s' and 'ALL DONE n_ok n_err' at the end."""
    out = Path(out)
    log = out / f"{kind}.log"
    marker_dir = out / ".inflight"
    marker_dir.mkdir(exist_ok=True)
    for stale in marker_dir.glob(".inflight_*"):
        stale.unlink()
    log_line(log, f"START {stamp()} kind={kind} workers={workers} heavy_slots={heavy_slots} todo={len(todo)} {header}")
    t_start = time.monotonic()
    n_ok = n_err = breaks = 0
    attempts_in_dead_pool = {}
    pending = list(todo)
    ctx = mp.get_context("spawn")
    while pending:
        pool = cf.ProcessPoolExecutor(max_workers=workers, mp_context=ctx, max_tasks_per_child=None,
                                      initializer=init_worker, initargs=(str(out), heavy_slots, kind))
        futures = {pool.submit(guarded, work, n, str(out), case_timeout, str(marker_dir)): n for n in pending}
        broken = False
        try:
            for fut in cf.as_completed(futures):
                n = futures[fut]
                try:
                    r = fut.result()
                except cf.process.BrokenProcessPool:
                    broken = True
                    break
                except BaseException as exc:  # noqa: BLE001
                    r = {"id": n, "status": "error", "seconds": 0.0, "error": f"{type(exc).__name__}: {exc}"}
                pending.remove(n)
                if r["status"] == "ok":
                    n_ok += 1
                    log_line(log, f"{stamp()} id={n:03d} status=ok seconds={r['seconds']:.1f} {r.get('detail', '')}".rstrip())
                else:
                    n_err += 1
                    log_line(log, f"{stamp()} id={n:03d} status=error seconds={r['seconds']:.1f}")
                    log_line(log, f"ERROR {stamp()} id={n:03d} {r['error']}")
        finally:
            pool.shutdown(wait=not broken, cancel_futures=True)
        if broken:
            markers = list(marker_dir.glob(".inflight_*"))
            running = [int(m.name.split("_")[1]) for m in markers]
            for m in markers:
                m.unlink(missing_ok=True)
            for n in running:
                attempts_in_dead_pool[n] = attempts_in_dead_pool.get(n, 0) + 1
                if attempts_in_dead_pool[n] >= 2 and n in pending:
                    pending.remove(n)
                    n_err += 1
                    msg = "worker process died twice on this case (killed, probably out of memory)"
                    on_dead_twice(n, msg)
                    log_line(log, f"{stamp()} id={n:03d} status=error seconds=0.0")
                    log_line(log, f"ERROR {stamp()} id={n:03d} {msg}")
            breaks += 1
            log_line(log, f"ERROR {stamp()} worker pool broke ({breaks}); restarting with {len(pending)} cases left")
            if breaks >= 20:
                log_line(log, f"ERROR {stamp()} giving up after 20 broken pools; rerun to resume")
                sys.exit(1)
            time.sleep(10)
    wall = time.monotonic() - t_start
    log_line(log, f"WALL_SECONDS {wall:.1f}")
    log_line(log, f"ALL DONE {n_ok} {n_err}")
    print(f"ALL DONE {n_ok} {n_err}  ({wall / 60:.1f} min)", flush=True)
    return n_ok, n_err, wall


def wall_seconds(log):
    """Sum of the WALL_SECONDS lines of a log (one per finished start); None when there is no log."""
    log = Path(log)
    if not log.exists():
        return None
    return sum(float(line.split()[1]) for line in log.read_text().splitlines() if line.startswith("WALL_SECONDS "))
