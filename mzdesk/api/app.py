"""FastAPI：/api/spec 给前端页面定义，/api/.../run 是占位接口，其余路径托管前端静态文件。"""
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from mzdesk.spec import MODULES, PIPELINES, spec

STATIC = Path(__file__).resolve().parents[1] / "static"
app = FastAPI(title="mzDesk API", version="0.1.1")
_modules = {m["id"]: m for m in MODULES}
_pipelines = {p["id"]: p for p in PIPELINES}


@app.get("/api/health")
def health():
    return {"status": "ok"}


@app.get("/api/spec")
def get_spec():
    return spec()


@app.post("/api/modules/{module_id}/run")
def run_module(module_id: str, inputs: dict):
    if module_id not in _modules:
        raise HTTPException(404, f"unknown module {module_id}")
    return {"status": "placeholder", "module": module_id, "received": inputs, "outputs": {}}


@app.post("/api/pipelines/{pipeline_id}/run")
def run_pipeline(pipeline_id: str, inputs: dict):
    if pipeline_id not in _pipelines:
        raise HTTPException(404, f"unknown pipeline {pipeline_id}")
    stages = [{"title": s["title"], "module": s["module"], "status": "queued"} for s in _pipelines[pipeline_id]["stages"]]
    return {"status": "placeholder", "pipeline": pipeline_id, "query": inputs, "stages": stages}


if (STATIC / "index.html").exists():
    app.mount("/assets", StaticFiles(directory=STATIC / "assets"), name="assets")

    @app.get("/{path:path}")
    def spa(path: str):
        candidate = STATIC / path
        if path and candidate.is_file():
            return FileResponse(candidate)
        return FileResponse(STATIC / "index.html")
else:
    @app.get("/")
    def no_frontend():
        return JSONResponse({"error": "frontend not built: run `npm run build` in web/"}, status_code=503)
