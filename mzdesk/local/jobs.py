"""Install jobs: one at a time, kept in memory (lost on restart, old ones pruned), with a progress record the dialog polls.

    job = STORE.start(["libraries", "msfinder"])      # a job already running is returned instead of starting a second one
    STORE.get(job.id).snapshot()                      # the shape of GET /api/local/jobs/{id}
    STORE.cancel(job.id)                              # kills the child process tree
"""
import collections
import threading
import time
import uuid

from mzdesk.local import components, progress
from mzdesk.local.runner import Cancelled, Ctx, InstallError

LOG_LINES = 30                  # the tail the dialog shows
KEEP_JOBS, TTL_S = 10, 3600     # at most 10 jobs are kept; a finished job is dropped an hour after it ended


class Job:
    def __init__(self, ids, on_line=None):
        self.id = uuid.uuid4().hex[:12]
        self.status = "running"                       # running | done | error | cancelled
        self.started, self.ended = time.monotonic(), None
        self.rows = [{"id": cid, "title": components.REGISTRY[cid].title, "status": "waiting", "detail": "", "step": "", "progress": None,
                      "got_mb": None, "total_mb": None} for cid in ids]
        self.log = collections.deque(maxlen=LOG_LINES)
        self.error = ""
        self.lock = threading.Lock()
        self.thread = None
        self.on_line = on_line                        # the terminal command prints every line as it comes
        self.ctx = Ctx(self._write)

    def _write(self, line):
        with self.lock:
            self.log.append(line)
            change = progress.parse(line)
            for row in self.rows:
                if row["status"] == "running":
                    row["detail"] = line
                    if change:
                        row.update(change)
        if self.on_line:
            self.on_line(line)

    def _set(self, cid, status, detail=None):
        with self.lock:
            for row in self.rows:
                if row["id"] == cid:
                    row["status"] = status
                    if detail is not None:
                        row["detail"] = detail
                    if status != "running":
                        row.update(step="", progress=None, got_mb=None, total_mb=None)

    def _finish(self, status, error=""):
        with self.lock:
            self.status, self.error, self.ended = status, error, time.monotonic()

    def snapshot(self):
        with self.lock:
            current = next((r["id"] for r in self.rows if r["status"] == "running"), None)
            return {"status": self.status, "current": current, "elapsed_s": round((self.ended or time.monotonic()) - self.started, 1),
                    "components": [dict(r) for r in self.rows], "log": list(self.log), "error": self.error}

    def finished(self):
        return self.status != "running"


def run_job(job, ids):
    """The job's thread: install the components one after the other; the first failure ends the job (the rest stay "waiting")."""
    try:
        for cid in ids:
            comp = components.REGISTRY[cid]
            if comp.installed():
                job._set(cid, "done", "already installed")
                continue
            job._set(cid, "running", "starting")
            try:
                comp.install(job.ctx)
            except Cancelled:
                job._set(cid, "waiting", "cancelled")
                job._finish("cancelled")
                return
            except InstallError as exc:
                job._set(cid, "error", str(exc).splitlines()[0])
                job._finish("error", str(exc))
                return
            except Exception as exc:  # noqa: BLE001 - the thread must report, not die
                job._set(cid, "error", f"{type(exc).__name__}: {exc}")
                job._finish("error", f"{comp.title}: {type(exc).__name__}: {exc}")
                return
            job._set(cid, "done", "installed")
        job._finish("done")
    finally:
        if not job.finished():
            job._finish("error", "The install stopped unexpectedly.")


class JobStore:
    def __init__(self):
        self.lock = threading.Lock()
        self.jobs = {}

    def _prune(self, now):
        for jid in [j for j, job in self.jobs.items() if job.finished() and now - job.ended > TTL_S]:
            del self.jobs[jid]
        while len(self.jobs) > KEEP_JOBS:
            old = min((j for j in self.jobs.values() if j.finished()), key=lambda j: j.ended, default=None)
            if old is None:
                break
            del self.jobs[old.id]

    def running(self):
        return next((j for j in self.jobs.values() if not j.finished()), None)

    def start(self, ids, on_line=None):
        """-> the Job. A job that is still running is returned as it is (one install at a time). InstallError for an unknown or unsupported id.
        on_line(text) gets every log line (the terminal command); job.thread.join() waits for the end."""
        with self.lock:
            self._prune(time.monotonic())
            job = self.running()
            if job is not None:
                return job
            order = components.plan(ids)
            job = Job(order, on_line)
            self.jobs[job.id] = job
            self._prune(time.monotonic())
            job.thread = threading.Thread(target=run_job, args=(job, order), daemon=True, name=f"install-{job.id}")
            job.thread.start()
        return job

    def get(self, job_id):
        with self.lock:
            return self.jobs.get(job_id)

    def cancel(self, job_id):
        job = self.get(job_id)
        if job is None:
            return False
        job.ctx.cancel()
        return True


STORE = JobStore()
