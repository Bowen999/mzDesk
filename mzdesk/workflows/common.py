"""Shared plumbing of the workflows in this package: the step runner and the generic result.

A workflow module (targeted_method.py, ei_ms_annotation.py) has PIPELINE_ID, STEP_TITLES ({step id: title}, from the spec stages), a Query
dataclass, parse_inputs(inputs) -> Query (raises mzdesk.api.common.InputError, HTTP 400) and run_query(query, progress) -> result.
Typical run_query:

    steps = StepRunner(STEP_TITLES, progress)
    got = steps.guarded(1, lambda: (STATUS_DONE, "3 names found", [table_block("names", "Names", ["Name"], rows)], names))
    steps.skip(2, "Not needed for GC-MS input.")
    return make_result(PIPELINE_ID, steps, make_answer("The compound is ...", facts=[["Platform", "LC-MS"]]), headline_blocks, warnings)

Generic result (what the job API returns as job["result"] and PipelinePage.tsx draws for every pipeline except esi_ms2_annotation):

    {"pipeline": PIPELINE_ID,
     "steps": [{"id": int, "title": str, "status": "done" | "skipped" | "failed", "seconds": float, "note": str}, ...],   # all steps, in order
     "answer": {"statement": str,                    # one or two plain sentences, shown large
                "level": str | None,                 # short badge text (a confidence level, a verdict), or None
                "level_label": str | None,           # longer text for the badge tooltip, or None
                "facts": [[label, value], ...]},     # a compact two-column list under the statement; values are strings
     "blocks": [block, ...],                         # the headline blocks, drawn under the answer by OutputList
     "warnings": [str, ...]}

A block is the spec block together with its data, exactly like the blocks a step hands to progress():
    {"type": "table" | "cards" | "chips" | "portal" | "notes" | "figure" | ..., "id": str, "title": str, "columns"?: [...], "data": {...}}
Every block needs a unique "id" (the page keys the data by it; a block without an id is not drawn). The data shapes are the ones in
mzdesk/spec.py and mzdesk/api/blocks.py (table_data, badge, link, row, chips_data, portal_data, ...). A table:

    rows = [["Formic acid", blk.badge("good", "good")], ["Ammonium acetate", "possible"]]
    table_block("additives", "Additives", ["Additive", "Fit"], rows)
    -> {"type": "table", "id": "additives", "title": "Additives", "columns": [...], "data": blk.table_data(rows, count=2)}

The blocks of a step go to progress() through StepRunner.finish / guarded and show inside that step's panel on the page while the run goes
on; result["blocks"] are the ones worth showing first, under the answer (the same block may be in both).
progress(step_id, status, note, blocks): status is running | done | skipped | failed; note / blocks are None when unchanged.
"""
from __future__ import annotations

import logging
import time

from mzdesk.api import blocks as blk

log = logging.getLogger("mzdesk.workflows")

STATUS_DONE, STATUS_SKIPPED, STATUS_FAILED = "done", "skipped", "failed"


class StepRunner:
    """Per-step bookkeeping (status, seconds, note) and the live progress calls; the workflow-package twin of annotation/pipeline.py _Steps.

    records: {step id: {"id", "title", "status", "seconds", "note"}}, every step "queued" until it starts."""

    def __init__(self, step_titles, progress=None):
        self.titles = dict(step_titles)
        self.progress = progress
        self.records = {n: {"id": n, "title": title, "status": "queued", "seconds": 0.0, "note": ""} for n, title in self.titles.items()}
        self.started = {}

    def emit(self, n, status, note=None, blocks=None):
        """One progress(step_id, status, note, blocks) call; a listener that raises never stops the run."""
        if self.progress:
            try:
                self.progress(n, status, note, blocks)
            except Exception:           # noqa: BLE001 - a broken listener must not stop the run
                log.exception("progress callback failed")

    def start(self, n):
        self.started[n] = time.perf_counter()
        self.records[n]["status"] = "running"
        self.emit(n, "running")

    def finish(self, n, status, note, blocks=(), seconds=None):
        if seconds is None:
            seconds = time.perf_counter() - self.started.get(n, time.perf_counter())
        self.records[n].update(status=status, note=note, seconds=round(seconds, 3))
        self.emit(n, status, note, list(blocks))

    def skip(self, n, note):
        self.records[n].update(status=STATUS_SKIPPED, note=note)
        self.emit(n, STATUS_SKIPPED, note, [])

    def guarded(self, n, fn):
        """Run fn() -> (status, note, blocks, value) as step n; an exception makes the step "failed" and the run goes on.
        Returns the value (None when the step failed)."""
        self.start(n)
        try:
            status, note, blocks, value = fn()
        except Exception as exc:        # noqa: BLE001 - a step must never take the workflow down
            log.exception("step %s failed", n)
            status, note, blocks, value = STATUS_FAILED, f"{self.titles[n]} failed: {type(exc).__name__}: {exc}", [], None
        self.finish(n, status, note, blocks)
        return value

    def failed(self):
        """Records of the steps that failed, in order (for the result's warnings)."""
        return [r for r in self.records.values() if r["status"] == STATUS_FAILED]

    def as_list(self):
        """The result's "steps": one record per step, in order."""
        return [dict(r) for r in self.records.values()]


# ---------------------------------------------------------------- blocks and the result
def block(block_type, block_id, title, data, **spec_keys):
    """A spec block with its data: {"type", "id", "title", "data", **spec_keys} (spec_keys: columns, items, style, help, caption, span, ...)."""
    return {"type": block_type, "id": block_id, "title": title, **spec_keys, "data": data}


def table_block(block_id, title, columns, rows, count=None, help=None, **spec_keys):
    """A table block. rows: lists of cells (str | number | blk.badge(...) | blk.link(...) | ...), see mzdesk/api/blocks.py."""
    return block("table", block_id, title, blk.table_data(rows, count=len(rows) if count is None else count, help=help), columns=list(columns), **spec_keys)


def notes_block(block_id, title, items, **spec_keys):
    """A list of short sentences."""
    return block("notes", block_id, title, {"items": [str(i) for i in items]}, **spec_keys)


def make_answer(statement, level=None, level_label=None, facts=()):
    """The result's answer: statement (str), an optional level badge and its longer label, and facts = [[label, value], ...]."""
    return {"statement": str(statement), "level": level, "level_label": level_label, "facts": [[str(k), str(v)] for k, v in facts]}


def make_result(pipeline_id, steps, answer, blocks=(), warnings=()):
    """The generic result (module docstring). steps: a StepRunner or a list of step records; answer: make_answer() or a dict with the same keys.
    A block without type, id or data is dropped with a warning instead of breaking the page."""
    records = steps.as_list() if isinstance(steps, StepRunner) else [dict(r) for r in steps]
    warnings = [str(w) for w in warnings]
    shown, seen = [], set()
    for b in blocks:
        if isinstance(b, dict) and b.get("type") and b.get("id") and isinstance(b.get("data"), dict) and b["id"] not in seen:
            seen.add(b["id"])
            shown.append(b)
        else:
            log.error("workflow %s: dropped a malformed or duplicate block: %r", pipeline_id, (b.get("id"), b.get("type")) if isinstance(b, dict) else b)
            warnings.append("A result block was dropped because it was malformed.")
    a = answer if isinstance(answer, dict) else {}
    return {"pipeline": pipeline_id, "steps": records,
            "answer": {"statement": str(a.get("statement") or ""), "level": a.get("level"), "level_label": a.get("level_label"),
                       "facts": [[str(k), str(v)] for k, v in a.get("facts") or []]},
            "blocks": shown, "warnings": warnings}
