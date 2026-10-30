import assert from "node:assert/strict";
import { mkdtempSync, rmSync } from "node:fs";
import { tmpdir } from "node:os";
import { dirname, join } from "node:path";
import { after, test } from "node:test";
import { fileURLToPath } from "node:url";
import React from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { createServer } from "vite";

const root = dirname(dirname(fileURLToPath(import.meta.url)));
const cache = mkdtempSync(join(tmpdir(), "mzdesk-progress-test-"));
const server = await createServer({ configFile: false, root, cacheDir: cache, server: { middlewareMode: true, hmr: false }, appType: "custom" });
const { WorkflowProgress } = await server.ssrLoadModule("/src/components/WorkflowProgress.tsx");
after(async () => { await server.close(); rmSync(cache, { recursive: true, force: true }); });

const titles = ["Adduct", "Formula", "Library search", "MS-FINDER structures", "FRIGID structures", "Merge and re-rank", "Isomer check", "Answer"];
const steps = statuses => statuses.map((status, i) => ({ id: i + 1, title: titles[i], status, seconds: 0, note: "", blocks: [] }));
const job = (status, statuses, extra = {}) => ({ status, steps: steps(statuses), note: "", error: null, result: null, ...extra });
const view = (job, extra = {}) => renderToStaticMarkup(React.createElement(WorkflowProgress, { job, ...extra }));
const assertPercent = (html, percent) => {
  assert.match(html, new RegExp(`aria-valuenow="${percent}"`));
  assert.match(html, new RegExp(`width:${percent}%`));
  assert.match(html, /not elapsed-time progress/);
  assert.equal((html.match(/role="progressbar"/g) || []).length, 1);
};

test("submission and the first poll have no invented percentage", () => {
  const submitting = view(null);
  assert.match(submitting, /Submitting workflow/);
  assert.match(submitting, /Stage progress unavailable/);
  assert.doesNotMatch(submitting, /aria-valuenow=/);
  assert.match(view(null, { awaitingFirstUpdate: true }), /Waiting for workflow progress/);
});

test("queued state reports the queue position and zero ended stages", () => {
  const html = view(job("queued", Array(8).fill("queued"), { note: "waiting for the previous run (position 2 in the queue)" }));
  assertPercent(html, 0);
  assert.match(html, /0% · 0 \/ 8 stages ended/);
  assert.match(html, /position 2 in the queue/);
  assert.equal((html.match(/position 2 in the queue/g) || []).length, 1);
  assert.match(html, /href="#workflow-step-1"/);
});

test("parallel running stages and their notes appear together", () => {
  const j = job("running", ["running", "running", ...Array(6).fill("queued")]);
  j.steps[1].note = "working";
  const html = view(j);
  assertPercent(html, 0);
  assert.match(html, /Running: Adduct \+ Formula/);
  assert.match(html, /Formula: working/);
  assert.match(html, /Step 2: Formula, running/);
});

test("done, skipped and failed stages count, and current-stage notes remain visible", () => {
  const j = job("running", ["done", "done", "skipped", "failed", "done", "running", "queued", "queued"]);
  j.steps[5].note = "scored 3/9";
  const html = view(j);
  assertPercent(html, 62.5);
  assert.match(html, /5 \/ 8 stages ended/);
  assert.match(html, /Running: Merge and re-rank/);
  assert.match(html, /Merge and re-rank: scored 3\/9/);
  assert.match(html, /Step 4: MS-FINDER structures, failed/);
});

test("all stages can end before finalization; job completion with failures is explicit", () => {
  const statuses = ["done", "done", "skipped", "failed", "done", "done", "done", "done"];
  const finalizing = view(job("running", statuses));
  assertPercent(finalizing, 100);
  assert.match(finalizing, /All stages ended. Finalizing result/);
  assert.doesNotMatch(finalizing, /Workflow completed/);
  const done = view(job("done", statuses));
  assertPercent(done, 100);
  assert.match(done, /Workflow completed with 1 failed stage\./);
  assert.match(view(job("done", Array(8).fill("done"))), /Workflow completed\./);
});

test("job failure and lost polling are not presented as successful runs", () => {
  const failed = view(job("error", ["failed", ...Array(7).fill("queued")]));
  assertPercent(failed, 12.5);
  assert.match(failed, /Workflow failed/);
  const stopped = view(job("running", ["done", "running", ...Array(6).fill("queued")]), { trackingStopped: true });
  assertPercent(stopped, 12.5);
  assert.match(stopped, /Progress tracking stopped/);
  assert.match(stopped, /last update/);
  assert.doesNotMatch(stopped, /Running: Formula/);
  assert.match(view(null, { trackingStopped: true }), /Progress tracking stopped/);
});

test("empty or unknown stage data never produce a percentage", () => {
  const empty = view(job("running", []));
  assert.match(empty, /Stage progress unavailable/);
  assert.match(empty, /Running; stage details unavailable/);
  assert.doesNotMatch(empty, /aria-valuenow=/);
  const unknown = view(job("running", ["done", "unknown"]));
  assert.match(unknown, /Stage progress unavailable/);
  assert.doesNotMatch(unknown, /aria-valuenow=/);
  assertPercent(view(job("running", ["done", "running", "queued"])), 33.3);
});

test("unknown workflow status stops tracking without inventing a running or completed state", () => {
  const html = view(job("unknown", ["running", "queued"]));
  assertPercent(html, 0);
  assert.match(html, /Workflow status unavailable; progress tracking stopped/);
  assert.doesNotMatch(html, /Running: Adduct|Workflow completed|Preparing next stage/);
});

test("missing stage data or status is unavailable, and stopped indeterminate bars are static", () => {
  const missing = view(job("running", [], { steps: null }));
  assert.match(missing, /Stage progress unavailable/);
  assert.doesNotMatch(missing, /aria-valuenow=/);
  const noStatus = view(job("running", ["running", null]));
  assert.match(noStatus, /Step 2: Formula, unknown/);
  assert.doesNotMatch(noStatus, /aria-valuenow=/);
  const stopped = view(null, { trackingStopped: true });
  assert.match(stopped, /is-indeterminate is-idle/);
  assert.doesNotMatch(view(null), /is-idle/);
  assert.match(view(job("error", [])), /is-indeterminate is-idle/);
});
