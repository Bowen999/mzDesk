import assert from "node:assert/strict";
import { after, test } from "node:test";
import { startSsr } from "./ssr-helper.mjs";

const ssr = await startSsr("steps");
const { StepList, StepView } = await ssr.load("/src/components/StepList.tsx");
after(() => ssr.close());

const table = { id: "formula", type: "table", title: "Formula detail", columns: ["Formula", "Mass"], data: { rows: [["C8H10N4O2", "194.0804"]] } };
const step = (id, title, status, extra = {}) => ({ id, title, status, seconds: 1.234, note: "", blocks: [], ...extra });
const steps = [
  step(1, "Adduct", "done", { note: "Scored 7 adducts; kept [M+H]+ (0.56).", blocks: [table] }),
  step(2, "Formula", "done", { blocks: [table], seconds: 4.936, note: "BUDDY with the MS1 isotope pattern: 4 formulas for [M+H]+ and 3 formulas for [M+NH4]+ in the first pass of the search." }),
  step(3, "Library search", "running", { seconds: 0 }),
  step(4, "MS-FINDER structures", "skipped", { note: "Library gate passed.", seconds: 0 }),
  step(5, "FRIGID structures", "failed", { note: "FRIGID stopped: the engine returned no candidates for this formula.", blocks: [table] }),
  step(6, "Merge and re-rank", "queued", { seconds: 0 }),
];
const view = list => ssr.render(StepList, { steps: list });
const li = (html, id) => new RegExp(`<li id="workflow-step-${id}".*?</li>`, "s").exec(html)[0];

test("the steps are one timeline: a heading with the totals, an ordered list, a node per step", () => {
  const html = view(steps);
  assert.equal((html.match(/<section class="timeline"/g) ?? []).length, 1);
  assert.match(html, /<h2 class="tl-title">Steps<\/h2><span class="tl-sum">6 steps · 1 skipped · 7\.4 s<\/span>/, "skipped and queued steps add no time");
  assert.equal((html.match(/<li id="workflow-step-/g) ?? []).length, 6);
  assert.equal((html.match(/<span class="tl-node" aria-hidden="true">/g) ?? []).length, 6);
  assert.doesNotMatch(html, /class="badge/, "the status is the node, not a chip");
});

test("no step is open by default, whatever its status or content", () => {
  const html = view(steps);
  assert.equal((html.match(/<details /g) ?? []).length, 3, "the three steps with blocks are details");
  assert.doesNotMatch(html, /<details[^>]*\sopen/);
  assert.doesNotMatch(html, /\sopen=""/);
});

test("each step is one line: title, status word, note, seconds with a bar scaled to the longest step", () => {
  const html = view(steps);
  const first = /<li id="workflow-step-1".*?<\/summary>/s.exec(html)[0];
  assert.match(first, /<span class="step-title">Adduct<\/span><span class="sr-only">done<\/span>/);
  assert.match(first, /<span class="step-note" title="Scored 7 adducts; kept \[M\+H\]\+ \(0\.56\)\.">Scored 7 adducts/);
  assert.match(first, /<span class="step-seconds">1\.2 s<\/span><i class="tl-bar" style="width:25%">/, "1.2 s of the longest 4.9 s");
  assert.match(li(html, 2), /<i class="tl-bar" style="width:100%">/, "the longest step has the full bar");
  assert.match(li(html, 3), /<span class="tl-word">running<\/span>/);
  assert.match(li(html, 4), /<span class="step-title">MS-FINDER structures<\/span><span class="tl-word">skipped<\/span>/);
  assert.doesNotMatch(li(html, 4), /step-seconds|tl-bar/, "no seconds for a skipped step");
  assert.doesNotMatch(li(html, 6), /step-seconds|tl-bar/, "no seconds for a queued step");
});

test("the blocks of a step are inside it and shown only once it is opened", () => {
  const html = view(steps);
  const body = li(html, 1);
  assert.match(body, /<\/summary><div class="step-body">/);
  assert.ok(body.includes("C8H10N4O2"));
  assert.ok(!/<summary.*C8H10N4O2.*<\/summary>/s.test(body), "the table is not in the folded line");
});

test("a long note is cut to one line in the folded step and written out in full inside it", () => {
  const html = view(steps);
  const long = steps[1].note;
  const second = li(html, 2);
  assert.match(second, /<p class="step-full-note">BUDDY with the MS1 isotope pattern/);
  assert.ok(second.includes(`title="${long}"`), "the whole note is the tooltip of the folded line");
  assert.doesNotMatch(li(html, 1), /step-full-note/, "a short note is not repeated");
});

test("a failed step shows its note in the folded line, in the alert colour class", () => {
  const html = view(steps);
  const failed = /<li id="workflow-step-5".*?<\/summary>/s.exec(html)[0];
  assert.match(failed, /data-status="failed"/);
  assert.match(failed, /<span class="tl-word">failed<\/span>/);
  assert.match(failed, /<span class="step-note is-failed"[^>]*>FRIGID stopped: the engine returned no candidates for this formula\.<\/span>/);
  assert.doesNotMatch(failed, /<details[^>]*\sopen/);
});

test("steps with nothing to open are plain lines, and every step keeps its anchor id", () => {
  const html = view(steps);
  for (const s of steps) assert.ok(html.includes(`id="workflow-step-${s.id}"`));
  assert.match(li(html, 3), /<li id="workflow-step-3" class="tl-step" data-status="running"><div class="tl-line is-static">/);
});

test("missing or odd step data does not break the line", () => {
  const odd = { id: 1, title: "", status: null, seconds: null, note: undefined, blocks: undefined };
  const html = ssr.render(StepView, { step: odd, index: 0 });
  assert.match(html, /Unnamed stage/);
  assert.match(html, />unknown</);
  assert.doesNotMatch(html, /step-seconds/);
  assert.equal(view(undefined), "");
  assert.equal(view([]), "", "no steps, no timeline");
});

test("one step reads '1 step', not '1 steps'", () => {
  assert.match(view([step(1, "Only", "done")]), /<span class="tl-sum">1 step · 1\.2 s<\/span>/);
});
