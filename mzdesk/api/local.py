"""The local version's API: /api/local/* (first-use downloads) and /api/update (is a newer mzDesk out). Empty on mzdesk.com (MZDESK_HOSTED=1).

  GET  /api/local/status            -> {"local", "platform", "components": [{id, title, size, installed, supported, ask, note}]}
  POST /api/local/needs             {"target": "workflow:<id>" | "module:<id>", "inputs": {form values}} -> {"missing", "blocked", "unsupported"}
  POST /api/local/install           {"ids": [component ids]} -> {"job_id"}     (one job at a time; a second call gets the running job's id)
  GET  /api/local/jobs/{id}         -> {"status", "current", "elapsed_s", "components": [{id, title, status, detail}], "log", "error"}
  POST /api/local/jobs/{id}/cancel  -> {"ok": true}
  GET  /api/update                  -> {"current", "latest", "newer", "released", "notes", "command", "checked_at"}   (mzdesk/update.py)
The POSTs start downloads and programs on this computer, so they refuse a Host that is not localhost / 127.0.0.1 / [::1] (DNS rebinding) and
a cross-site Origin. The request body never reaches a command line: ids are looked up in the registry (mzdesk/local/components.py).
"""
import os
import re
import sys
from urllib.parse import urlsplit

from fastapi import APIRouter, HTTPException, Request

from mzdesk import update
from mzdesk.local import components, jobs, targets
from mzdesk.local.runner import InstallError

router = APIRouter()
LOCAL_HOSTS = ("localhost", "127.0.0.1", "[::1]", "::1")


def hosted():
    """The public website: every route here answers as if it did not exist (read per request, like the other MZDESK_ switches of this package)."""
    return os.environ.get("MZDESK_HOSTED") == "1"


def host_name(host_header):
    """'localhost:8731' -> 'localhost'; '[::1]:8731' -> '[::1]'; '' for nothing usable."""
    host = (host_header or "").strip().lower()
    if host.startswith("["):
        return host.partition("]")[0] + "]"
    return host.rpartition(":")[0] if re.search(r":\d+$", host) else host


def guard(request: Request):
    """403 unless the request comes from a page of this computer's own server: Host is a loopback name, and an Origin header (browsers send one on
    cross-site and on POST requests) names the same host and port. A request without Origin (curl, the terminal) passes the Host check alone."""
    host = request.headers.get("host", "")
    if host_name(host) not in LOCAL_HOSTS:
        raise HTTPException(403, "This address is not allowed.")
    origin = request.headers.get("origin")
    if origin is not None:
        parts = urlsplit(origin)
        if origin == "null" or parts.scheme != request.url.scheme or parts.netloc.lower() != host.strip().lower():
            raise HTTPException(403, "Cross-site requests are not allowed.")


def local_only():
    if hosted():
        raise HTTPException(404, "not found")


@router.get("/api/local/status")
def local_status():
    if hosted():
        return {"local": False}
    return {"local": True, "platform": sys.platform if sys.platform in ("darwin", "win32") else "linux", "components": components.status()}


@router.post("/api/local/needs")
def local_needs(body: dict, request: Request):
    local_only()
    guard(request)
    target = body.get("target")
    if not isinstance(target, str):
        raise HTTPException(400, "target must be a string such as workflow:esi_ms2_annotation or module:ms2_esi.")
    try:
        return targets.needs(target, body.get("inputs"))
    except targets.UnknownTarget as exc:
        raise HTTPException(404, str(exc)) from None


@router.post("/api/local/install")
def local_install(body: dict, request: Request):
    local_only()
    guard(request)
    ids = body.get("ids")
    if not isinstance(ids, list) or not ids or not all(isinstance(i, str) for i in ids):
        raise HTTPException(400, "ids must be a list of component ids.")
    unknown = [i for i in ids if components.get(i) is None]
    if unknown:
        raise HTTPException(404, f"Unknown component: {', '.join(unknown)[:200]}")
    try:
        job = jobs.STORE.start(ids)
    except InstallError as exc:
        raise HTTPException(400, str(exc)) from None
    return {"job_id": job.id}


@router.get("/api/local/jobs/{job_id}")
def local_job(job_id: str):
    local_only()
    job = jobs.STORE.get(job_id)
    if job is None:
        raise HTTPException(404, "unknown or expired job")
    return job.snapshot()


@router.post("/api/local/jobs/{job_id}/cancel")
def local_cancel(job_id: str, request: Request):
    local_only()
    guard(request)
    if not jobs.STORE.cancel(job_id):
        raise HTTPException(404, "unknown or expired job")
    return {"ok": True}


@router.get("/api/update")
def get_update():
    return update.info()          # never raises; on mzdesk.com and with MZDESK_NO_UPDATE_CHECK=1 it answers "no newer version" without fetching anything
