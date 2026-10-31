import assert from "node:assert/strict";
import { mkdtempSync, readFileSync, rmSync } from "node:fs";
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
const M = await server.ssrLoadModule("/src/components/workflowProgressModel.ts");
after(async () => { await server.close(); rmSync(cache, { recursive: true, force: true }); });

const titles = ["Adduct", "Formula", "Library search", "MS-FINDER structures", "FRIGID structures", "Merge and re-rank", "Isomer check", "Answer"];
const steps = (statuses, seconds = 0) => statuses.map((status, i) => ({ id: i + 1, title: titles[i], status, seconds, note: "", blocks: [] }));
const job = (status, statuses, extra = {}) => ({ status, steps: steps(statuses), note: "", error: null, result: null, ...extra });
const view = (job, extra = {}) => renderToStaticMarkup(React.createElement(WorkflowProgress, { job, ...extra }));
const n8 = (done, running = 0) => [...Array(done).fill("done"), ...Array(running).fill("running"), ...Array(8 - done - running).fill("queued")];
const S = 1000;

/** Drive the pure tracker through [seconds, phase, statuses, seconds-per-step?] samples; returns what the bar shows ("percent"), or the target it heads for.
 *  `expected` (optional) is the stages' expected_s, in order. */
function run(samples, { start = 0, field = "percent", expected = null } = {}) {
  let tracker = M.newTracker(start * S);
  return samples.map(([t, phase, statuses = [], secs = 0]) => {
    const list = steps(statuses, secs);
    if (expected) list.forEach((st, i) => { st.expected_s = expected[i]; });
    const out = M.track(tracker, { phase, steps: list }, t * S);
    tracker = out.tracker;
    return out[field];
  });
}
const targets = (samples, opts = {}) => run(samples, { ...opts, field: "target" });
const nondecreasing = values => values.every((v, i) => i === 0 || v === null || values[i - 1] === null || v >= values[i - 1] - 1e-9);

// ------------------------------------------------------------------ the model

test("the bar is never at 0: submitting, waiting and queued creep from about 2% towards 6%", () => {
  const p = run([[0, "submitting"], [1, "submitting"], [5, "waiting"], [20, "waiting"], [60, "queued", n8(0)], [600, "queued", n8(0)]]);
  assert.ok(p[0] >= 2 && p[0] < 2.5, `starts near 2%, got ${p[0]}`);
  assert.ok(p.every(v => v >= 2 && v < 6), `stays in [2, 6): ${p}`);
  assert.ok(p.every((v, i) => i === 0 || v > p[i - 1]), "creeps while it waits");
  assert.equal(M.creep(0), 2);
  assert.ok(M.creep(1e6) <= 6);
});

test("progress is non-decreasing across any sequence of states and times", () => {
  let seed = 7;
  const rnd = () => (seed = (seed * 1664525 + 1013904223) % 4294967296) / 4294967296;
  const pick = a => a[Math.floor(rnd() * a.length)];
  for (let trial = 0; trial < 40; trial++) {
    let t = 0;
    const samples = [];
    for (let i = 0; i < 80; i++) {
      t += rnd() * 3;
      const phase = i === 0 ? "waiting" : pick(["queued", "running", "running", "running", "stopped", "unknown"]);
      samples.push([t, phase, Array.from({ length: 8 }, () => pick(["queued", "running", "done", "done", "skipped", "failed"])), pick([0, 1, 5, 30])]);
    }
    const p = run(samples);
    assert.ok(nondecreasing(p), `trial ${trial} decreased: ${p.map(v => v?.toFixed(1))}`);
    assert.ok(p.every(v => v >= 2 && v <= 99), "inside [2, 99] until the job is done");
  }
  // a dedicated regression: a poll that reports fewer ended stages than the previous one
  const p = run([[0, "running", n8(4)], [1, "running", n8(1)], [2, "running", n8(0, 1)]]);
  assert.equal(p[1], p[0]);
  assert.equal(p[2], p[0]);
});

test("within a running stage the value rises with time, stays under the next slice floor and never completes the slice", () => {
  const slice = 94 / 8, floor = e => M.PRE_CEIL + slice * e;
  for (const e of [0, 3, 7]) {
    const times = [0, 1, 2, 4, 8, 16, 32, 120, 3600];
    const p = targets(times.map(t => [t, "running", n8(e, 1)]));
    assert.ok(Math.abs(p[0] - floor(e)) < 1e-9, "starts at its own slice floor");
    assert.ok(p.every((v, i) => i === 0 || v > p[i - 1]), `strictly increasing for ${e} ended: ${p}`);
    assert.ok(p.every(v => v < floor(e + 1)), `below the next floor for ${e} ended: ${p}`);
    assert.ok(p.at(-1) > floor(e) + slice * 0.89, "approaches the cap (90% of the slice)");
  }
  assert.equal(M.stageFraction(0, 8), 0);
  assert.ok(M.stageFraction(1e9, 8) <= M.RUN_CAP);
});

test("a stage that ends glides the bar to the next slice floor", () => {
  const p = targets([[0, "running", n8(0, 1)], [30, "running", n8(0, 1)], [31, "running", n8(1, 1)], [31.1, "running", n8(1, 0)], [32, "running", n8(1, 1)]]);
  assert.ok(p[2] > p[1], "next stage started: ahead of the slice floor");
  assert.ok(p[3] >= M.PRE_CEIL + 94 / 8, "at least the new floor");
  assert.ok(nondecreasing(p));
});

test("tau adapts to the finished stages, within sane limits", () => {
  assert.equal(M.stageTau(steps(n8(0, 1))), M.TAU_DEFAULT);
  assert.equal(M.stageTau([{ status: "done", seconds: 4 }, { status: "done", seconds: 6 }, { status: "running", seconds: 99 }]), 5);
  assert.equal(M.stageTau([{ status: "done", seconds: 0 }, { status: "skipped", seconds: 50 }, { status: "failed", seconds: 50 }]), M.TAU_DEFAULT);
  assert.equal(M.stageTau([{ status: "done", seconds: 0.1 }]), M.TAU_MIN);
  assert.equal(M.stageTau([{ status: "done", seconds: 9999 }]), M.TAU_MAX);
  assert.equal(M.stageTau([{ status: "done", seconds: null }, { status: "done", seconds: NaN }]), M.TAU_DEFAULT);
  // slower finished stages -> the running stage fills more slowly
  const fast = targets([[0, "running", n8(2, 1), 2], [10, "running", n8(2, 1), 2]]);
  const slow = targets([[0, "running", n8(2, 1), 40], [10, "running", n8(2, 1), 40]]);
  assert.ok(fast[1] - fast[0] > slow[1] - slow[0]);
});

test("done is 100; every other state stays at 99 or below, even with all stages ended", () => {
  const ended = ["done", "done", "skipped", "failed", "done", "done", "done", "done"];
  const p = targets([[0, "running", n8(0, 1)], [5, "running", ended], [6, "done", ended]]);
  assert.equal(p[1], 99);
  assert.equal(p[2], 100);
  assert.equal(run([[0, "running", n8(0, 1)], [5, "running", ended], [6, "done", ended]])[2], 100, "done is shown at once, not glided to");
  assert.deepEqual(run([[0, "done", []]]), [100]);
});

test("job error and stopped tracking freeze the bar where it is", () => {
  const failedAt4 = ["done", "done", "done", "failed", "queued", "queued", "queued", "queued"];
  const err = run([[0, "running", n8(3, 1)], [12, "running", n8(3, 1)], [13, "error", failedAt4], [90, "error", failedAt4]]);
  assert.equal(err[2], err[1], "an error does not move the bar");
  assert.equal(err[3], err[1], "and neither does time");
  const stopped = run([[0, "running", n8(2, 1)], [9, "running", n8(2, 1)], [10, "stopped", n8(2, 1)], [400, "stopped", n8(2, 1)]]);
  assert.equal(stopped[2], stopped[1]);
  assert.equal(stopped[3], stopped[1]);
  const none = run([[0, "submitting"], [8, "waiting"], [9, "stopped"], [500, "stopped"]]);
  assert.equal(none[2], none[1]);
  assert.equal(none[3], none[1]);
  assert.deepEqual(run([[0, "stopped"]]), [M.PRE_START]);
  // first drawn in a frozen state: the stages that ended
  const cold = (phase, statuses) => run([[0, phase, statuses]])[0];
  assert.equal(cold("unknown", n8(1)), M.PRE_CEIL + 94 / 8);
  assert.equal(cold("error", n8(3, 1)), M.PRE_CEIL + 94 * 3 / 8, "no time-based growth for the running stage");
  assert.equal(cold("error", []), M.PRE_START);
});

test("several running stages add up", () => {
  const t = 8, f = M.stageFraction(t, M.TAU_DEFAULT);
  const one = targets([[0, "running", n8(2, 1)], [t, "running", n8(2, 1)]])[1];
  const two = targets([[0, "running", n8(2, 2)], [t, "running", n8(2, 2)]])[1];
  assert.ok(Math.abs(one - (M.PRE_CEIL + 94 * (2 + f) / 8)) < 1e-9);
  assert.ok(Math.abs(two - (M.PRE_CEIL + 94 * (2 + 2 * f) / 8)) < 1e-9);
  assert.ok(two > one && two < M.PRE_CEIL + 94 * 4 / 8, "under the floor after both end");
});

test("each running stage is timed from when it was first seen", () => {
  const late = targets([[0, "running", n8(0, 1)], [20, "running", ["done", "running", ...Array(6).fill("queued")]], [20, "running", ["done", "running", ...Array(6).fill("queued")]]]);
  const s2 = M.PRE_CEIL + 94 / 8;
  assert.ok(late[1] >= s2 && late[1] < s2 + 0.2, "stage 2 has just been seen: close to its floor, not credited for the 20 s before");
});

test("a new run starts over; a resumed run does not", () => {
  const ended = n8(8);
  const p = run([[0, "running", n8(5, 1)], [3, "done", ended], [10, "submitting"], [11, "waiting"], [12, "queued", n8(0)]]);
  assert.equal(p[1], 100);
  assert.ok(p[2] >= 2 && p[2] < 3, `restarted near 2%: ${p[2]}`);
  assert.ok(p[4] < 6);
  const failedRun = run([[0, "running", n8(1, 1)], [4, "error", [...n8(1), "failed", ...Array(6).fill("queued")]], [20, "running", n8(0, 1)]]);
  assert.ok(failedRun[2] < failedRun[1], "a job that replaces a failed one starts over");
  const stale = run([[0, "running", n8(4, 1)], [5, "stopped", n8(4, 1)], [10, "submitting"]]);
  assert.ok(stale[2] < stale[1], "submitting after a stopped run starts over");
  const resumed = run([[0, "running", n8(4, 1)], [5, "stopped", n8(4, 1)], [10, "running", n8(4, 1)]]);
  assert.ok(resumed[2] >= resumed[1]);
});

test("a running job without usable stage data is indeterminate, not a made-up number", () => {
  assert.deepEqual(run([[0, "running", []]]), [null]);
  assert.deepEqual(run([[0, "running", ["done", "unknown"]]]), [null]);
  const kept = run([[0, "running", n8(2, 1)], [1, "running", []], [2, "running", n8(2, 1)]]);
  assert.equal(kept[1], null);
  assert.ok(kept[2] >= kept[0]);
});

// ------------------------------------------------------------------ expected_s: weighted slices and the shape of a whole run

test("slices are proportional to expected_s; a stage without one gets the mean of the others; no expected_s at all means equal slices", () => {
  assert.deepEqual(M.stageWeights([{ expected_s: 2 }, { expected_s: 6 }, { expected_s: 2 }]), [2, 6, 2]);
  assert.deepEqual(M.stageWeights([{ expected_s: 2 }, {}, { expected_s: 6 }, { expected_s: null }]), [2, 4, 6, 4]);
  assert.deepEqual(M.stageWeights([{}, {}, {}]), [1, 1, 1]);
  for (const bad of [0, -3, NaN, Infinity, "5", null, undefined]) assert.deepEqual(M.stageWeights([{ expected_s: bad }, { expected_s: 3 }]), [3, 3], `${bad} is not an expected_s`);
  const e = [2, 6, 2];
  const p = targets([[0, "running", n8(0, 1).slice(0, 3)], [1, "running", ["done", "running", "queued"]], [30, "running", ["done", "done", "running"]]], { expected: e });
  assert.equal(p[0], M.PRE_CEIL);
  assert.ok(Math.abs(p[1] - (M.PRE_CEIL + 94 * 2 / 10)) < 1e-9, "stage 1 of weight 2 out of 10 is 20% of the band");
  assert.ok(Math.abs(p[2] - (M.PRE_CEIL + 94 * 8 / 10)) < 1e-9, "the first two stages own 80% of the band; stage 3 has only just started");
  // without expected_s the same run has equal thirds
  const eq = targets([[0, "running", ["running", "queued", "queued"]], [1, "running", ["done", "running", "queued"]]]);
  assert.ok(Math.abs(eq[1] - (M.PRE_CEIL + 94 / 3)) < 1e-9);
});

test("a stage that is skipped takes no share of the bar", () => {
  const e = [2, 40, 2];
  const p = targets([[0, "running", ["running", "queued", "queued"]], [1, "running", ["done", "skipped", "running"]]], { expected: e });
  assert.ok(Math.abs(p[1] - (M.PRE_CEIL + 94 * 2 / 4)) < 1e-9, "the 40 s stage is gone from both sides of the fraction");
  assert.equal(targets([[0, "running", ["skipped", "skipped"]]], { expected: [5, 5] })[0], M.ACTIVE_MAX, "nothing left to run");
  assert.ok(nondecreasing(run([[0, "running", ["running", "queued", "queued"]], [1, "running", ["done", "skipped", "running"]], [2, "running", ["done", "skipped", "running"]]], { expected: e })));
});

test("tau of a running stage is its expected_s, scaled by how fast the finished stages ran, within limits", () => {
  const one = { expected_s: 10 };
  assert.equal(M.tauOf(one, [one]), 10);
  assert.equal(M.tauOf({ expected_s: 0.01 }, [{ expected_s: 0.01 }]), M.TAU_EXPECTED_MIN);
  assert.equal(M.tauOf({ expected_s: 5000 }, [{ expected_s: 5000 }]), M.TAU_EXPECTED_MAX);
  assert.equal(M.tauOf({}, [{}]), M.TAU_DEFAULT, "no expected_s: the fallback rule");
  assert.equal(M.tauOf({}, [{ status: "done", seconds: 20 }, {}]), 20);
  // at t = tau a stage has filled RUN_CAP * (1 - 1/e) of its slice
  const e = [10, 10];
  const [a, b] = targets([[0, "running", ["running", "queued"]], [10, "running", ["running", "queued"]]], { expected: e });
  assert.ok(Math.abs((b - a) - 94 / 2 * M.RUN_CAP * (1 - Math.exp(-1))) < 1e-9);
  // finished stages that ran 2x slower than expected slow the next one down; faster ones speed it up
  const slowSteps = [{ status: "done", seconds: 20, expected_s: 10 }, { status: "running", expected_s: 10 }];
  const fastSteps = [{ status: "done", seconds: 5, expected_s: 10 }, { status: "running", expected_s: 10 }];
  assert.ok(M.tauOf(slowSteps[1], slowSteps) > 10 && M.tauOf(fastSteps[1], fastSteps) < 10);
});

test("the speed factor is 1 before anything finished, ignores milliseconds of noise and stays within 0.25-4", () => {
  assert.equal(M.speedFactor([{ status: "running", expected_s: 5 }]), 1);
  assert.equal(M.speedFactor([{ status: "done", seconds: 5, expected_s: 5 }]), 1);
  assert.ok(Math.abs(M.speedFactor([{ status: "done", seconds: 0.01, expected_s: 0.1 }]) - 1) < 0.05, "a 0.1 s stage that took 0.01 s says nothing about the next stage");
  assert.equal(M.speedFactor([{ status: "done", seconds: 9999, expected_s: 1 }]), M.SPEED_MAX);
  assert.equal(M.speedFactor([{ status: "done", seconds: 0, expected_s: 1000 }]), M.SPEED_MIN);
  assert.equal(M.speedFactor([{ status: "skipped", seconds: 0, expected_s: 40 }, { status: "failed", seconds: 99, expected_s: 1 }, { status: "done", seconds: 3 }]), 1, "only finished stages with an expected_s count");
});

test("the bar follows a jump of its target at no more than MAX_RATE points per second", () => {
  const pct = run([[0, "running", n8(0, 1)], [0.5, "running", n8(6, 1)], [1.5, "running", n8(6, 1)], [2, "running", n8(6, 1)], [60, "running", n8(6, 1)]]);
  const tgt = targets([[0, "running", n8(0, 1)], [0.5, "running", n8(6, 1)], [1.5, "running", n8(6, 1)], [2, "running", n8(6, 1)], [60, "running", n8(6, 1)]]);
  assert.equal(pct[0], 6);
  assert.ok(Math.abs(pct[1] - (6 + M.MAX_RATE * 0.5)) < 1e-9);
  assert.ok(Math.abs(pct[2] - (6 + M.MAX_RATE * 1.5)) < 1e-9);
  assert.ok(Math.abs(pct[3] - (6 + M.MAX_RATE * 2)) < 1e-9);
  assert.equal(pct[4], tgt[4], "it gets there");
  assert.ok(tgt[1] > pct[1]);
});

/** A whole run as the page sees it: the server finishes stage i after trueDur[i] seconds (skipped ones take none), the page polls every 1.5 s and
 *  redraws every 0.5 s. Returns [seconds, phase, percent, target] per redraw, from submission until 4 s after completion. */
function simulate(expected, trueDur, skipped = [], pollOffset = 0.3) {
  const n = expected.length, sampleEvery = 0.5, pollEvery = 1.5, jobStart = 0.2;
  const startAt = []; let t = jobStart;
  expected.forEach((_, i) => { startAt[i] = t; if (!skipped.includes(i)) t += trueDur[i]; });
  const end = t;
  const snapshot = at => steps(expected.map((_, i) => skipped.includes(i) ? (at >= startAt[i] ? "skipped" : "queued")
    : at >= startAt[i] + trueDur[i] ? "done" : at >= startAt[i] ? "running" : "queued"), 0).map((st, i) => ({ ...st, expected_s: expected[i], seconds: st.status === "done" ? trueDur[i] : 0 }));
  let tracker = M.newTracker(0);
  const out = [];
  for (let at = 0; at <= end + 4; at += sampleEvery) {
    const k = Math.floor((at - pollOffset) / pollEvery), polled = pollOffset + k * pollEvery;
    const phase = k < 0 || polled < jobStart ? "waiting" : polled >= end ? "done" : "running";
    const r = M.track(tracker, { phase, steps: phase === "waiting" ? [] : snapshot(polled) }, at * S);
    tracker = r.tracker;
    out.push([at, phase, r.percent, r.target]);
  }
  return out;
}
const maxJump = (rows, field = 2) => Math.max(...rows.map((r, i) => i && r[1] !== "done" && r[field] !== null ? r[field] - rows[i - 1][field] : 0));
const GC_EXPECTED = [0.1, 7.0, 0.05, 0.2, 0.2, 0.5, 0.1, 0.05, 0.1];
const ESI_EXPECTED = [6, 0.5, 4, 20, 40, 0.5, 8, 0.1];

test("a GC-MS run dominated by one long stage fills steadily and never jumps before completion", () => {
  for (const offset of [0, 0.3, 0.7, 1.1, 1.4]) {
    for (const lib of [3, 7.5, 15]) {                                   // the library search: faster than, about, and twice as slow as expected
      const rows = simulate(GC_EXPECTED, [0.07, lib, 0, 0.01, 0.02, 0.13, 0.01, 0, 0.05], [7], offset);
      assert.ok(maxJump(rows) <= 15, `offset ${offset}, library ${lib} s: jumped ${maxJump(rows)} points`);
      assert.ok(nondecreasing(rows.map(r => r[2])));
      assert.equal(rows.at(-1)[2], 100);
      const first = rows.find(r => r[1] === "running");
      assert.ok(first[2] >= 2, "never at 0");
      // during the long stage the bar keeps rising: every second of it adds something
      const during = rows.filter(r => r[0] > 1.6 && r[0] < Math.min(lib, 7.5) - 0.5);
      assert.ok(during.every((r, i) => i === 0 || r[2] > during[i - 1][2]), `offset ${offset}, library ${lib} s: stalled while the library search ran`);
    }
  }
  const rows = simulate(GC_EXPECTED, [0.07, 7.5, 0, 0.01, 0.02, 0.13, 0.01, 0, 0.05], [7], 0.3);
  const at = t => rows.find(r => r[0] === t)[2];
  assert.ok(at(2) > 15 && at(4) > at(2) + 10 && at(6) > at(4) + 5, "close to linear while it runs, not stage steps");
  // a warm run (everything done within one poll) can only jump at completion
  const warm = simulate(GC_EXPECTED, [0.01, 0.28, 0, 0.01, 0.02, 0.08, 0.01, 0, 0.02], [7], 0.9);
  assert.ok(maxJump(warm) <= 15);
});

test("an ESI run whose long stages are skipped glides instead of jumping", () => {
  const cases = {
    "library match (MS-FINDER and FRIGID skipped)": [[10.3, 0.1, 7.5, 0, 0, 0.25, 5.07, 0], [3, 4]],
    "no library, MS-FINDER decides (FRIGID skipped)": [[9.6, 0.1, 0, 6.3, 0, 0.5, 4.5, 0], [2, 4]],
    "everything runs": [[12, 1, 3, 40, 70, 1, 10, 0.1], []]};
  for (const [name, [dur, skipped]] of Object.entries(cases)) {
    for (const offset of [0, 0.5, 1.1]) {
      const rows = simulate(ESI_EXPECTED, dur, skipped, offset);
      assert.ok(maxJump(rows) <= 15, `${name}, offset ${offset}: jumped ${maxJump(rows)} points`);
      assert.ok(nondecreasing(rows.map(r => r[2])));
      assert.ok(rows.filter(r => r[1] === "running").every(r => r[2] >= 2 && r[2] <= 99));
      assert.equal(rows.at(-1)[2], 100);
    }
  }
  // without the glide the skipped stages would be a jump of tens of points (this is what the rate limit is for)
  const raw = simulate(ESI_EXPECTED, cases["library match (MS-FINDER and FRIGID skipped)"][0], [3, 4], 0.5);
  assert.ok(maxJump(raw, 3) > 15 && maxJump(raw) <= 15);
});

test("phaseOf maps job state and flags", () => {
  assert.equal(M.phaseOf(null, false, false), "submitting");
  assert.equal(M.phaseOf(null, true, false), "waiting");
  assert.equal(M.phaseOf(null, true, true), "stopped");
  assert.equal(M.phaseOf({ status: "running" }, false, true), "stopped");
  assert.equal(M.phaseOf({ status: "queued" }, false, false), "queued");
  assert.equal(M.phaseOf({ status: "done" }, false, false), "done");
  assert.equal(M.phaseOf({ status: "error" }, false, false), "error");
  assert.equal(M.phaseOf({ status: "weird" }, false, false), "unknown");
  assert.ok(M.isActive("running") && M.isActive("queued") && !M.isActive("done") && !M.isActive("error") && !M.isActive("stopped") && !M.isActive("unknown"));
});

// ------------------------------------------------------------------ the markup

const single = html => {
  assert.equal((html.match(/<section/g) || []).length, 1, "one section");
  assert.equal((html.match(/role="progressbar"/g) || []).length, 1, "one progress bar");
  assert.doesNotMatch(html, /<ol|<li|wf-progress-track|wf-progress-node|Stage-count estimate|not elapsed-time|Workflow progress<\/h2>|<h2/);
};
const valuenow = html => { const m = html.match(/aria-valuenow="(\d+)"/); return m ? Number(m[1]) : null; };

test("submission and the first poll show a small head start, not 0 and not a made-up stage", () => {
  const submitting = view(null);
  single(submitting);
  assert.match(submitting, /Submitting workflow/);
  assert.equal(valuenow(submitting), 2);
  assert.match(submitting, />2%</);
  assert.match(submitting, /aria-valuemin="0"[^>]*aria-valuemax="100"/);
  assert.match(submitting, /aria-valuetext="2%\. Submitting workflow/);
  const waiting = view(null, { awaitingFirstUpdate: true });
  assert.match(waiting, /Waiting for workflow progress/);
  assert.equal(valuenow(waiting), 2);
});

test("queued shows the queue note once and a head start", () => {
  const html = view(job("queued", Array(8).fill("queued"), { note: "waiting for the previous run (position 2 in the queue)" }));
  single(html);
  assert.match(html, /Queued/);
  assert.equal((html.match(/position 2 in the queue/g) || []).length, 1);
  assert.ok(valuenow(html) >= 2 && valuenow(html) < 6);
});

test("a running stage reads 'Step n of N · title' with the percentage on the right and its note below", () => {
  const j = job("running", ["done", "done", "skipped", "failed", "done", "running", "queued", "queued"]);
  j.steps[5].note = "scored 3/9";
  const html = view(j);
  single(html);
  assert.match(html, /Step 6 of 8 · Merge and re-rank/);
  assert.match(html, /aria-valuenow="60"/);                       // 4 of the 7 stages that are not skipped have ended
  assert.match(html, /width:59\.7%/);
  assert.match(html, />60%</);
  assert.match(html, /aria-valuetext="60%\. Step 6 of 8 · Merge and re-rank"/);
  assert.match(html, /<p class="wf-progress-note">scored 3\/9<\/p>/);
  assert.equal((html.match(/wf-progress-note/g) || []).length, 1);
  assert.match(html, /role="status"/);
});

test("parallel stages are listed together with one prefixed note", () => {
  const j = job("running", ["running", "running", ...Array(6).fill("queued")]);
  j.steps[1].note = "working";
  const html = view(j);
  single(html);
  assert.match(html, /Steps 1, 2 of 8 · Adduct \+ Formula/);
  assert.match(html, /Formula: working/);
  assert.equal((html.match(/wf-progress-note/g) || []).length, 1);
  assert.ok(valuenow(html) >= 6);
});

test("expected_s decides how much of the bar a stage owns", () => {
  const j = job("running", ["done", "running"]);
  j.steps[0].expected_s = 1;
  j.steps[1].expected_s = 9;
  const html = view(j);
  single(html);
  assert.equal(valuenow(html), 15, "one tenth of the band after the short stage: 6 + 94 / 10");
  const equal = job("running", ["done", "running"]);
  assert.equal(valuenow(view(equal)), 53, "without expected_s: half of the band");
  assert.equal(valuenow(view(job("running", ["done", "skipped", "running"], { steps: [{ id: 1, title: "A", status: "done", expected_s: 1 }, { id: 2, title: "B", status: "skipped", expected_s: 50 }, { id: 3, title: "C", status: "running", expected_s: 1 }] }))), 53);
});

test("all stages ended but the job still running reads Finalizing at 99%; done reads Completed at 100%", () => {
  const statuses = ["done", "done", "skipped", "failed", "done", "done", "done", "done"];
  const finalizing = view(job("running", statuses));
  single(finalizing);
  assert.match(finalizing, /Finalizing result/);
  assert.equal(valuenow(finalizing), 99);
  const done = view(job("done", statuses));
  single(done);
  assert.match(done, /Completed with 1 failed stage/);
  assert.equal(valuenow(done), 100);
  assert.match(done, /width:100%/);
  assert.match(done, /wf-progress is-done/);
  const clean = view(job("done", Array(8).fill("done")));
  assert.match(clean, />Completed</);
  assert.equal(valuenow(clean), 100);
});

test("job failure and lost polling are stop-styled and say what happened", () => {
  const failed = view(job("error", ["failed", ...Array(7).fill("queued")]));
  single(failed);
  assert.match(failed, /Failed at Adduct/);
  assert.match(failed, /wf-progress is-error/);
  assert.equal(valuenow(failed), 18);
  assert.doesNotMatch(failed, /Completed/);
  assert.equal(valuenow(view(job("error", []))), 2);
  const stopped = view(job("running", ["done", "running", ...Array(6).fill("queued")]), { trackingStopped: true });
  single(stopped);
  assert.match(stopped, /Tracking stopped/);
  assert.match(stopped, /may still be running/);
  assert.match(stopped, /last update/);
  assert.match(stopped, /wf-progress is-stopped/);
  assert.doesNotMatch(stopped, /Step 2 of 8/);
  const early = view(null, { trackingStopped: true });
  assert.match(early, /Tracking stopped/);
  assert.equal(valuenow(early), 2);
});

test("unusable stage data is indeterminate and carries no aria-valuenow", () => {
  for (const j of [job("running", []), job("running", ["done", "unknown"]), job("running", [], { steps: null }), job("running", ["running", null])]) {
    const html = view(j);
    single(html);
    assert.match(html, /Running · stage details unavailable/);
    assert.match(html, /wf-progress-bar is-indeterminate/);
    assert.doesNotMatch(html, /aria-valuenow=/);
    assert.match(html, />—</);
  }
});

test("an unknown job status stops tracking without inventing a running or completed state", () => {
  const html = view(job("unknown", ["running", "queued"]));
  single(html);
  assert.match(html, /Status unavailable/);
  assert.doesNotMatch(html, /Step 1 of 8|Completed|Finalizing/);
  assert.match(html, /wf-progress is-stopped/);
});

test("untitled stages are named", () => {
  const j = job("running", ["running"]);
  j.steps[0].title = "";
  assert.match(view(j), /Step 1 of 1 · Unnamed stage/);
});

// ------------------------------------------------------------------ the stylesheet

test("the stylesheet is flat 2D: a two-tone left-to-right fill, flat track, no shadow, sheen, sweep or vertical gradient; both themes and reduced motion", () => {
  const css = readFileSync(join(root, "src/components/workflow-progress.css"), "utf8");
  assert.match(css, /--wf-from: #5c8bd6; --wf-to: #3262d2/);
  assert.match(css, /:root\[data-theme="dark"\] \.wf-progress \{[^}]*--wf-from: #6282cb; --wf-to: #7b97e0/);
  assert.match(css, /@media \(prefers-color-scheme: dark\) \{\s*:root:not\(\[data-theme="light"\]\) \.wf-progress \{[^}]*--wf-from: #6282cb/);
  assert.match(css, /--wf-track: #e8ecf4/);
  assert.match(css, /linear-gradient\(90deg, var\(--wf-from\), var\(--wf-to\)\)/);
  assert.match(css, /transition: width \.5s/);
  assert.match(css, /@media \(prefers-reduced-motion: reduce\)[\s\S]*transition: none[\s\S]*animation: none/);
  // flat: nothing that reads as 3D
  const bar = css.split("\n").filter(l => /\.wf-progress-bar/.test(l) || /^\s+(display: block; height: 100%|background: linear)/.test(l)).join("\n");
  assert.doesNotMatch(css, /--wf-sheen|--wf-sweep|--wf-track-line|--wf-mid|::before|::after|wf-progress-sweep/, "no sheen, sweep or inner line");
  assert.doesNotMatch(bar, /box-shadow|inset|filter|drop-shadow/, "no shadow on the bar, track or fill");
  assert.doesNotMatch(css, /(?:to bottom|to top|\b(?:0|180|360)deg)\b|repeating-linear-gradient|radial-gradient/, "no vertical, radial or striped gradient");
  assert.deepEqual([...css.matchAll(/linear-gradient\(([^,]+),/g)].map(m => m[1]), ["90deg"], "the only gradient runs left to right");
  assert.match(css, /\.wf-progress\.is-stopped \.wf-progress-bar > span \{ background: var\(--wf-stop\); \}/);
  assert.match(css, /\.wf-progress\.is-error \.wf-progress-bar > span \{ background: var\(--wf-err\); \}/);
  assert.match(css, /\.wf-progress-bar \{[^}]*height: 8px/);
  assert.match(css, /\.wf-progress-bar\.is-indeterminate > span \{[^}]*background: var\(--wf-to\)/, "the sliding segment is one flat colour");
});
