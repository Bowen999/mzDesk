"""FastAPI：/api/spec 给前端页面定义，/api/modules/{id}/run 运行模块，其余路径托管前端静态文件。

已接引擎的模块（spec 里 live=True）由 mzdesk/api/runners 计算；其余模块返回 placeholder。
Pipeline（ESI MS/MS Annotation、EI MS Annotation、Targeted Method）跑得久：POST /api/pipelines/{id}/run 起一个后台线程，GET /api/pipelines/jobs/{job_id} 查每一步的进度和结果。
引擎共用一台机器和一把重任务锁，所以 pipeline 任务一次只跑一个，后到的排队（status "queued"）。
"""
import copy
import functools
import logging
import os
import threading
import time
import uuid
from pathlib import Path

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import FileResponse, JSONResponse, Response
from fastapi.staticfiles import StaticFiles

from mzdesk import __version__
from mzdesk.annotation import pipeline as annotation_pipeline
from mzdesk.api import admin, reports, security
from mzdesk.api import feedback as feedback_api
from mzdesk.api import local as local_api
from mzdesk.api import release as release_api
from mzdesk.api import fragments as fragments_api
from mzdesk.api.request_limits import RequestLimit
from mzdesk.api import stats, verdicts
from mzdesk.api.common import EngineError, InputError, clean
from mzdesk.api.runners import RUNNERS
from mzdesk.shared import versions
from mzdesk.spec import MODULES, PIPELINES, is_deferred, spec
from mzdesk.workflows import ei_ms_annotation, targeted_method

ROOT = Path(__file__).resolve().parents[1]
STATIC = ROOT / "static"
DOCS = ROOT / "docs"                                   # one markdown file per page, written from the long help texts
LIB_FILES = {"3Dmol-min.js": ROOT / "lib" / "3Dmol-min.js"}   # 前端按需加载的第三方脚本
log = logging.getLogger("mzdesk.api")
# the public website (MZDESK_HOSTED=1 in deploy/mzdesk.service): the frontend drops the "Internet search" switch there,
# since what people type reaches this server anyway. The local version leaves it unset. Read when the server starts.
HOSTED = os.environ.get("MZDESK_HOSTED") == "1"
_docs = security.api_docs_enabled()                      # the API browser (/docs, /redoc, /openapi.json) is for local development only
app = FastAPI(title="mzDesk API", version=__version__, docs_url="/docs" if _docs else None, redoc_url="/redoc" if _docs else None,
              openapi_url="/openapi.json" if _docs else None)
# the last one added is the outermost: compression wraps the headers, which wrap the size limit (so a 413 carries them too)
app.add_middleware(RequestLimit)
app.add_middleware(security.SecurityHeaders)
app.add_middleware(GZipMiddleware, minimum_size=1000)   # plotly.js is 4.8 MB raw, 1.5 MB gzipped
app.include_router(fragments_api.router)
app.include_router(local_api.router)        # the local version: first-use downloads, update notice (mzdesk/local/)
app.include_router(release_api.router)      # mzdesk.com only: install scripts, version.json, the release wheel
reports.install_log_handler()                           # log.exception(...) anywhere in mzdesk also files a crash report
_modules = {m["id"]: m for m in MODULES if not is_deferred(m)}     # deferred modules (DMPK, ICP-MS, MRM) answer 404
_pipelines = {p["id"]: p for p in PIPELINES}


@app.get("/api/health")
def health():
    return {"status": "ok"}


@app.get("/api/spec")
def get_spec():
    return {**spec(), "version": __version__, "hosted": HOSTED}


@app.get("/api/versions")
def get_versions():
    """What the results depend on: mzDesk, Python, the engines, models, data and online services (names, versions, dates, counts; no paths).
    Cached for ten minutes; the Session report writes it into report.json (shared/versions.py)."""
    return versions.environment()


# The general Docs pages (Getting Started, Concepts, Reference), served next to the module and workflow notes. Keep equal to GENERAL_PAGES
# in web/src/docs/outline.ts (web/tests/docs-outline.test.mjs checks it).
GENERAL_PAGES = (
    "overview", "quick-start", "local-version",
    "workflows-and-modules", "inputs", "reading-results", "session-report", "privacy",
    "cli", "faq", "changelog", "licenses", "feedback",
)
_DOCS_IMG_CHARS = frozenset("abcdefghijklmnopqrstuvwxyz0123456789-")       # a figure name is [a-z0-9-]+


@app.get("/api/docs/{page_id}")
def get_docs(page_id: str):
    """The text of one Docs page: mzdesk/docs/<id>.md for a general page, a module or a workflow, drawn by the frontend's /docs/<id> page."""
    if page_id not in GENERAL_PAGES and page_id not in _modules and page_id not in _pipelines:
        raise HTTPException(404, "unknown page")
    path = DOCS / f"{page_id}.md"
    if not path.is_file():
        raise HTTPException(404, "no text for this page")
    return {"id": page_id, "markdown": path.read_text(encoding="utf-8")}


@app.get("/api/docs-img/{name}")
def get_docs_img(name: str):
    """A figure of the Docs: mzdesk/docs/img/<name>.svg (the frontend draws it inline). Scripts cannot run if the file is opened on its own."""
    path = DOCS / "img" / f"{name}.svg"
    if not name or not set(name) <= _DOCS_IMG_CHARS or not path.is_file():
        raise HTTPException(404, "no such figure")
    return Response(path.read_bytes(), media_type="image/svg+xml",
                    headers={"Cache-Control": "public, max-age=300", "Content-Security-Policy": "default-src 'none'; style-src 'unsafe-inline'; sandbox"})


@app.post("/api/modules/{module_id}/run")
def run_module(module_id: str, inputs: dict, request: Request):
    """Runs the module and records the run (module, status, time) for the stats page."""
    t0, status = time.perf_counter(), 200
    try:
        out = _run_module(module_id, inputs)
        status = out.status_code if isinstance(out, JSONResponse) else 200
        return out
    except HTTPException as exc:
        status = exc.status_code
        raise
    finally:
        if module_id in _modules:
            stats.record(request, "run", module=module_id, status=status, ms=round(1000 * (time.perf_counter() - t0)))


def _run_module(module_id, inputs):
    if module_id not in _modules:
        raise HTTPException(404, f"unknown module {module_id}")
    runner = RUNNERS.get(module_id)
    if runner is None:
        return {"status": "placeholder", "module": module_id, "outputs": {}, "warnings": []}
    try:
        result = runner(inputs)
    except InputError as exc:
        raise HTTPException(400, str(exc)) from None
    except EngineError as exc:
        return JSONResponse(status_code=502, content={"detail": str(exc), "log": exc.log})
    except Exception as exc:  # noqa: BLE001 — report engine bugs instead of a bare 500
        ref = reports.record_crash(exc, f"module {module_id}", inputs)
        raise HTTPException(500, f"{module_id} failed: {type(exc).__name__}: {security.scrub(exc)} (ref {ref})") from None
    outputs = clean(result.get("outputs", {}))
    return {"status": "ok", "module": module_id, "outputs": outputs, "warnings": result.get("warnings", []),
            "verdict": result.get("verdict") or verdicts.build(module_id, inputs, outputs)}


@app.post("/api/track")
def track(body: dict, request: Request):
    """Page view beacon from the frontend (path inside the app, document.referrer)."""
    path = str(body.get("path") or "")[:200]
    if path and not path.startswith("/stats"):
        stats.record(request, "view", path=path, ref=str(body.get("ref") or ""))
    return {"ok": True}


@app.post("/api/stats")
def get_stats(body: dict, request: Request):
    """Usage numbers for the stats page; needs the site password (checked here, never sent back)."""
    admin.require(body.get("password"), request)
    return stats.summary(body.get("days") or 30)


@app.post("/api/feedback")
async def send_feedback(request: Request):
    """The sidebar's Send feedback form: JSON with a message and up to three images, delivered by email.

    Returns {"ok", "delivered", "note"}: delivered=false means the message was saved as an .eml on the server
    instead (SMTP not configured, or the mail server was unreachable). See mzdesk/api/feedback.py."""
    return await feedback_api.submit(request)


@app.post("/api/crash")
async def send_crash(request: Request):
    """A page that failed to draw in the browser: {"message", "stack", "page"}. Stored as a crash report; the same crash again only counts up."""
    return await feedback_api.submit_crash(request)


@app.post("/api/reports")
def list_reports(body: dict, request: Request):
    """The reports page (admin password): feedback and crashes, newest first, without their screenshots."""
    admin.require(body.get("password"), request)
    return {"reports": reports.listing(), "channels": reports.channels()}


@app.post("/api/reports/{ref}")
def get_report(ref: str, body: dict, request: Request):
    """One report in full, screenshots included as data URLs."""
    admin.require(body.get("password"), request)
    found = reports.detail(ref)
    if found is None:
        raise HTTPException(404, "unknown report")
    return found


@app.post("/api/reports/{ref}/delete")
def delete_report(ref: str, body: dict, request: Request):
    admin.require(body.get("password"), request)
    if not reports.delete(ref):
        raise HTTPException(404, "unknown report")
    return {"ok": True}


@app.exception_handler(Exception)
async def unexpected(request: Request, exc: Exception):
    """Whatever no route caught: stored as a crash report; the caller gets the reference, not the internals."""
    ref = reports.record_crash(exc, f"{request.method} {request.url.path}", log_it=False)      # the server's own log already has the traceback
    return JSONResponse(status_code=500, content={"detail": f"Something went wrong (ref {ref})."})


# ---------------------------------------------------------------- pipeline jobs (in memory, lost on restart)
JOB_LIMIT, JOB_TTL_S = 20, 3600           # at most 20 jobs are kept; a job is dropped an hour after its last update
# pipeline id -> module with parse_inputs() and run_query(); the workflows package modules return the generic result (mzdesk/workflows/common.py)
_PIPELINE_RUNNERS = {m.PIPELINE_ID: m for m in (annotation_pipeline, ei_ms_annotation, targeted_method)}
# The engines (BUDDY, MS-FINDER, FRIGID) share one machine and one heavy-job lock that gives up after 20 s, so two workflow runs at once make
# steps 4-5 fail. Runs therefore go one at a time, first come first served: a job waits with status "queued" until it is first in _order and
# nothing is _active. _turn wakes the waiting threads whenever _active or _order changes.
_jobs_lock = threading.Lock()
_turn = threading.Condition(_jobs_lock)
_jobs, _order, _active = {}, [], None      # job id -> job; ids that have not started, oldest first; id of the run in progress
_FINISHED = ("done", "skipped", "failed")


def _waiting_note(job_id):
    return f"waiting for the previous run (position {_order.index(job_id) + 1} in the queue)"


def _prune_jobs(now, make_room=False):
    """Drop jobs older than JOB_TTL_S; with make_room also the oldest finished ones until a new job fits. Call with the lock held."""
    for job_id in [j for j, job in _jobs.items() if now - job["touched"] > JOB_TTL_S]:
        del _jobs[job_id]
    for job_id in sorted((j for j, job in _jobs.items() if job["status"] in ("done", "error")), key=lambda j: _jobs[j]["touched"]):
        if not make_room or len(_jobs) < JOB_LIMIT:
            break
        del _jobs[job_id]


def _start_job(pipeline_id, runner, query, request):
    """Register the job and start its thread -> (job id, "running" | "queued"): queued when another run is in progress or ahead of it."""
    now = time.time()
    steps = [{"id": s["id"], "title": s["title"], "status": "queued", "seconds": 0.0, "note": "", "blocks": [], "expected_s": s["expected_s"]}
             for s in _pipelines[pipeline_id]["stages"]]
    job_id = uuid.uuid4().hex[:12]
    with _jobs_lock:
        _prune_jobs(now, make_room=True)
        if len(_jobs) >= JOB_LIMIT:
            raise HTTPException(429, f"{JOB_LIMIT} workflow runs are in progress; try again in a minute.")
        queued = _active is not None or bool(_order)
        _order.append(job_id)
        _jobs[job_id] = {"status": "queued" if queued else "running", "steps": steps, "result": None, "error": None, "touched": now, "started": {},
                         "note": _waiting_note(job_id) if queued else ""}
        steps[0]["note"] = _jobs[job_id]["note"]
    threading.Thread(target=_run_job, args=(job_id, pipeline_id, runner, query, request), daemon=True, name=f"job-{job_id}").start()
    return job_id, "queued" if queued else "running"


def _wait_for_turn(job_id):
    """Block until it is this job's turn, then mark it as the active run. False if the job was dropped while it waited."""
    global _active
    with _turn:
        while job_id in _jobs and (_active is not None or _order[0] != job_id):
            job = _jobs[job_id]
            job["note"] = job["steps"][0]["note"] = _waiting_note(job_id)
            _turn.wait()
        _order.remove(job_id)
        if job_id not in _jobs:
            _turn.notify_all()
            return False
        _active = job_id
        job = _jobs[job_id]
        job.update(status="running", note="", touched=time.time())
        job["steps"][0]["note"] = ""
        _turn.notify_all()                              # the others move up one position
        return True


def _release_turn():
    global _active
    with _turn:
        _active = None
        _turn.notify_all()


def _job_progress(job_id):
    """The progress(step_id, status, note, blocks) callback of one job."""
    def progress(step_id, status, note, blocks):
        now = time.time()
        with _jobs_lock:
            job = _jobs.get(job_id)
            step = next((s for s in job["steps"] if s["id"] == step_id), None) if job else None
            if step is None:
                return
            if status == "running":
                job["started"].setdefault(step_id, now)
            step["status"] = status
            if note is not None:
                step["note"] = security.strip_paths(note)
            if blocks is not None:
                step["blocks"] = blocks
            if status in _FINISHED:
                step["seconds"] = round(now - job["started"].get(step_id, now), 3)
            job["touched"] = now
    return progress


def _run_job(job_id, pipeline_id, runner, query, request):
    if not _wait_for_turn(job_id):
        return
    try:
        _run_job_now(job_id, pipeline_id, runner, query, request)
    finally:
        _release_turn()


def _run_job_now(job_id, pipeline_id, runner, query, request):
    t0, error, result = time.perf_counter(), None, None
    try:
        result = clean(runner.run_query(query, _job_progress(job_id)))
    except Exception as exc:  # noqa: BLE001 - the job reports it; the thread must not die silently
        ref = reports.record_crash(exc, f"pipeline {pipeline_id}")
        error = f"{pipeline_id} failed: {type(exc).__name__}: {security.scrub(exc)} (ref {ref})"
    stats.record(request, "run", module=pipeline_id, status=500 if error else 200, ms=round(1000 * (time.perf_counter() - t0)))
    with _jobs_lock:                                      # last: a client that sees "done" knows the run is over
        job = _jobs.get(job_id)
        if job:
            for step in job["steps"]:
                done = next((s for s in result["steps"] if s["id"] == step["id"]), None) if result else None
                if done:                                  # the run's own timing and note are exact
                    step.update(status=done["status"], seconds=done["seconds"], note=security.strip_paths(done["note"]))
                elif step["status"] == "running":
                    step.update(status="failed", note=error or "")
            job.update(status="error" if error else "done", result=result, error=error, touched=time.time())


@app.post("/api/pipelines/{pipeline_id}/run")
def run_pipeline(pipeline_id: str, inputs: dict, request: Request):
    """Check the form, start the run in a background thread and return its job id; progress and result: GET /api/pipelines/jobs/{job_id}."""
    runner = _PIPELINE_RUNNERS.get(pipeline_id)
    if pipeline_id not in _pipelines or runner is None:
        raise HTTPException(404, f"unknown pipeline {pipeline_id}")
    try:
        query = runner.parse_inputs(inputs)
    except InputError as exc:
        raise HTTPException(400, str(exc)) from None
    except Exception as exc:  # noqa: BLE001
        ref = reports.record_crash(exc, f"pipeline {pipeline_id} form", inputs)
        raise HTTPException(500, f"{pipeline_id} failed: {type(exc).__name__}: {security.scrub(exc)} (ref {ref})") from None
    job_id, status = _start_job(pipeline_id, runner, query, request)
    return {"job_id": job_id, "status": status}


@app.get("/api/pipelines/jobs/{job_id}")
def get_pipeline_job(job_id: str):
    """{"status": queued | running | done | error, "note": why it waits (queued) else "",
    "steps": [{id, title, status: queued | running | done | skipped | failed, seconds, note, blocks}],
    "result": the pipeline's Result when done, else null, "error": message or null}.
    Runs go one at a time: a job behind another one is "queued" (note "waiting for the previous run (position n in the queue)"; the same
    text is the note of step 1) until it starts."""
    with _jobs_lock:
        now = time.time()
        _prune_jobs(now)
        job = _jobs.get(job_id)
        if job is None:
            raise HTTPException(404, "unknown or expired job")
        steps = copy.deepcopy(job["steps"])
        for step in steps:                                 # a running step shows its elapsed time
            if step["status"] == "running":
                step["seconds"] = round(now - job["started"].get(step["id"], now), 1)
        return {"status": job["status"], "note": job["note"], "steps": steps, "result": copy.deepcopy(job["result"]), "error": job["error"]}


MAX_DEPICT_SMILES = 1500


@functools.lru_cache(maxsize=1024)
def _depict(smiles: str) -> str | None:
    from mzdesk.api.common import depiction_svg
    try:
        return depiction_svg(smiles, size=(260, 182))
    except Exception:  # noqa: BLE001 - display only; a SMILES typed by a person is not a crash
        log.warning("depiction failed", exc_info=True)
        return None


@app.get("/api/depict")
def depict(smiles: str = ""):
    """2D structure of a SMILES as SVG, for the hover previews of the result pages (404 when there is nothing to draw)."""
    smiles = smiles.strip()
    if not smiles or len(smiles) > MAX_DEPICT_SMILES:
        raise HTTPException(404, "no structure")
    svg = _depict(smiles)
    if not svg:
        raise HTTPException(404, "no structure")
    return Response(svg, media_type="image/svg+xml",
                    headers={"Cache-Control": "public, max-age=3600", "Content-Security-Policy": "default-src 'none'; style-src 'unsafe-inline'; sandbox"})


@app.get("/lib/{name}")
def lib_file(name: str):
    path = LIB_FILES.get(name)
    if path is None or not path.is_file():
        raise HTTPException(404, "not found")
    return FileResponse(path, media_type="text/javascript")


class _ImmutableFiles(StaticFiles):
    """Vite 给 assets 文件名加了内容哈希，浏览器可以一直缓存。"""

    async def get_response(self, path, scope):
        response = await super().get_response(path, scope)
        if response.status_code == 200:
            response.headers["Cache-Control"] = "public, max-age=31536000, immutable"
        return response


if (STATIC / "index.html").exists():
    app.mount("/assets", _ImmutableFiles(directory=STATIC / "assets"), name="assets")

    @app.get("/{path:path}")
    def spa(path: str):
        try:
            candidate = (STATIC / path).resolve()
            is_file = bool(path) and candidate.is_file() and STATIC.resolve() in candidate.parents
        except (ValueError, OSError):                      # a NUL byte or a name too long: it is just another page address
            is_file = False
        if is_file:
            return FileResponse(candidate)
        # the page itself must never be cached: it names the current (hashed) assets
        return FileResponse(STATIC / "index.html", headers={"Cache-Control": "no-cache"})
else:
    @app.get("/")
    def no_frontend():
        return JSONResponse({"error": "frontend not built: run `npm run build` in web/"}, status_code=503)
