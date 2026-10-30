"""Fragment view jobs (the "View" pop-up of Alternative explanations and ESI MS/MS Annotation).

POST /api/fragments/run starts one MetFrag run in a background thread and returns its job id; GET /api/fragments/jobs/{job_id} reports
{"status": queued | running | done | error, "percent", "stage", "seconds", "result", "error", "log"} until it is done. The run is made when the
user clicks: nothing is computed ahead of time. The percentage follows the stages of fragment_view.build (MetFrag itself is one call, so the
bar waits at its stage until it returns). The same request twice is answered from a small in-memory cache; jobs are lost on restart.

Request: {"smiles", "adduct" (MetFrag's spelling), "peaks": [[m/z, % of the base peak]], "precursor_mz", "tol_da"?, "tol_ppm"?, "tree_depth"?, "name"?}.
"""
import collections
import hashlib
import json
import logging
import math
import threading
import time
import uuid

from fastapi import APIRouter, HTTPException

from mzdesk.api.common import clean
from mzdesk.ms2 import fragment_view
from mzdesk.ms2 import metfrag_engine as mf

log = logging.getLogger("mzdesk.api")
router = APIRouter()

MAX_RUNNING = 2             # MetFrag runs at once; the next ones wait ("queued")
JOB_LIMIT = 60              # jobs kept in memory
JOB_TTL_S = 900             # a job is dropped 15 minutes after its last update
CACHE_SIZE = 24
MAX_SMILES = 3000
MAX_NAME = 200

_lock = threading.Lock()
_slots = threading.BoundedSemaphore(MAX_RUNNING)
_jobs = {}
_cache = collections.OrderedDict()      # request key -> finished result


def _number(body, key, default, low, high):
    value = body.get(key, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not low <= value <= high:
        raise HTTPException(400, f"{key} must be a number from {low:g} to {high:g}.")
    return float(value)


def parse_request(body):
    """Request JSON -> the arguments of fragment_view.build, or HTTP 400 with a message for the user."""
    smiles = body.get("smiles")
    if not isinstance(smiles, str) or not smiles.strip() or len(smiles) > MAX_SMILES:
        raise HTTPException(400, "smiles must be a structure of at most %d characters." % MAX_SMILES)
    adduct = body.get("adduct")
    if adduct not in mf.ADDUCTS:
        raise HTTPException(400, f"MetFrag cannot use the adduct {adduct!r}. Supported: {', '.join(mf.ADDUCTS)}.")
    raw = body.get("peaks")
    if not isinstance(raw, list) or not raw or len(raw) > fragment_view.MAX_QUERY_PEAKS:
        raise HTTPException(400, f"peaks must be a list of 1 to {fragment_view.MAX_QUERY_PEAKS} [m/z, intensity] pairs.")
    peaks = []
    for p in raw:
        ok = (isinstance(p, (list, tuple)) and len(p) == 2 and all(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) for v in p)
              and p[0] > 0 and p[1] > 0)
        if not ok:
            raise HTTPException(400, "Every peak must be [m/z, intensity] with both numbers above 0.")
        peaks.append((float(p[0]), float(p[1])))
    name = body.get("name") or ""
    if not isinstance(name, str):
        raise HTTPException(400, "name must be text.")
    if body.get("precursor_mz") is None:
        raise HTTPException(400, "precursor_mz is required.")
    return {"smiles": smiles.strip(), "adduct": adduct, "peaks": peaks, "precursor_mz": _number(body, "precursor_mz", None, 1, 20000),
            "tol_da": _number(body, "tol_da", 0.01, 0.0005, 0.5), "tol_ppm": _number(body, "tol_ppm", 0.0, 0, 100),
            "tree_depth": int(_number(body, "tree_depth", 2, 1, 3)), "name": name[:MAX_NAME]}


def _key(req):
    return hashlib.sha1(json.dumps(req, sort_keys=True).encode()).hexdigest()


def _prune(now):
    """Drop expired jobs; call with the lock held."""
    for job_id in [j for j, job in _jobs.items() if now - job["touched"] > JOB_TTL_S]:
        del _jobs[job_id]


def _update(job_id, **fields):
    with _lock:
        job = _jobs.get(job_id)
        if job:
            if "percent" in fields:
                fields["percent"] = max(job["percent"], fields["percent"])      # the bar never goes back
            job.update(fields, touched=time.time())


def _run(job_id, key, req):
    """The job thread: wait for a free slot, run MetFrag, store the result or the error."""
    with _slots:
        _update(job_id, status="running", stage="Starting", percent=1, started=time.time())
        try:
            result = clean(fragment_view.build(
                req["smiles"], req["peaks"], req["adduct"], req["precursor_mz"], tol_da=req["tol_da"], tol_ppm=req["tol_ppm"],
                tree_depth=req["tree_depth"], name=req["name"],
                progress=lambda fraction, text: _update(job_id, percent=round(100 * fraction), stage=text)))
        except ValueError as exc:
            _update(job_id, status="error", error=str(exc))
        except mf.MetFragError as exc:
            _update(job_id, status="error", error=str(exc), log=exc.log)
        except Exception as exc:  # noqa: BLE001 - the job reports it; the thread must not die silently
            log.exception("fragment view failed")
            _update(job_id, status="error", error=f"The fragment view failed: {type(exc).__name__}: {exc}")
        else:
            with _lock:
                _cache[key] = result
                while len(_cache) > CACHE_SIZE:
                    _cache.popitem(last=False)
            _update(job_id, status="done", percent=100, stage="Done", result=result)


@router.post("/api/fragments/run")
def run_fragments(body: dict):
    req = parse_request(body)
    key, now, job_id = _key({**req, "peaks": [list(p) for p in req["peaks"]]}), time.time(), uuid.uuid4().hex[:12]
    with _lock:
        _prune(now)
        if len(_jobs) >= JOB_LIMIT:
            raise HTTPException(503, "Too many fragment views are running or waiting. Try again in a minute.")
        cached = _cache.get(key)
        if cached is not None:
            _cache.move_to_end(key)
        _jobs[job_id] = {"status": "done" if cached is not None else "queued", "percent": 100 if cached is not None else 0,
                         "stage": "Done" if cached is not None else "Waiting for a free slot", "result": cached, "error": None, "log": "",
                         "created": now, "started": now, "touched": now}
    status = "queued" if cached is None else "done"          # read before the thread starts: it may finish first
    if cached is None:
        threading.Thread(target=_run, args=(job_id, key, req), daemon=True, name=f"fragments-{job_id}").start()
    return {"job_id": job_id, "status": status}


@router.get("/api/fragments/jobs/{job_id}")
def get_fragment_job(job_id: str):
    with _lock:
        _prune(time.time())
        job = _jobs.get(job_id)
        if job is None:
            raise HTTPException(404, "unknown or expired job")
        running = job["status"] in ("queued", "running")
        return {"status": job["status"], "percent": job["percent"], "stage": job["stage"], "error": job["error"], "log": job["log"],
                "seconds": round((time.time() if running else job["touched"]) - job["created"], 1), "result": job["result"]}
