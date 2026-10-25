"""FastAPI：/api/spec 给前端页面定义，/api/modules/{id}/run 运行模块，其余路径托管前端静态文件。

已接引擎的模块（spec 里 live=True）由 mzdesk/api/runners 计算；其余模块和 pipeline 返回 placeholder。
"""
import logging
import time
from pathlib import Path

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from mzdesk import __version__
from mzdesk.api import stats, verdicts
from mzdesk.api.common import EngineError, InputError, clean
from mzdesk.api.runners import RUNNERS
from mzdesk.spec import MODULES, PIPELINES, spec

ROOT = Path(__file__).resolve().parents[1]
STATIC = ROOT / "static"
LIB_FILES = {"3Dmol-min.js": ROOT / "lib" / "3Dmol-min.js"}   # 前端按需加载的第三方脚本
log = logging.getLogger("mzdesk.api")
app = FastAPI(title="mzDesk API", version=__version__)
app.add_middleware(GZipMiddleware, minimum_size=1000)   # plotly.js is 4.8 MB raw, 1.5 MB gzipped
_modules = {m["id"]: m for m in MODULES}
_pipelines = {p["id"]: p for p in PIPELINES}


@app.get("/api/health")
def health():
    return {"status": "ok"}


@app.get("/api/spec")
def get_spec():
    return {**spec(), "version": __version__}


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
        log.exception("module %s failed", module_id)
        raise HTTPException(500, f"{module_id} failed: {type(exc).__name__}: {exc}") from None
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
def get_stats(body: dict):
    """Usage numbers for the stats page; needs the site password (checked here, never sent back)."""
    if not stats.check_password(body.get("password")):
        time.sleep(0.6)        # slows guessing; nginx also rate-limits this path
        raise HTTPException(401, "Wrong password.")
    return stats.summary(body.get("days") or 30)


@app.post("/api/pipelines/{pipeline_id}/run")
def run_pipeline(pipeline_id: str, inputs: dict):
    if pipeline_id not in _pipelines:
        raise HTTPException(404, f"unknown pipeline {pipeline_id}")
    stages = [{"title": s["title"], "modules": s["modules"], "status": "queued"} for s in _pipelines[pipeline_id]["stages"]]
    return {"status": "placeholder", "pipeline": pipeline_id, "query": inputs, "stages": stages}


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
        candidate = (STATIC / path).resolve()
        if path and candidate.is_file() and STATIC.resolve() in candidate.parents:
            return FileResponse(candidate)
        # the page itself must never be cached: it names the current (hashed) assets
        return FileResponse(STATIC / "index.html", headers={"Cache-Control": "no-cache"})
else:
    @app.get("/")
    def no_frontend():
        return JSONResponse({"error": "frontend not built: run `npm run build` in web/"}, status_code=503)
