import assert from "node:assert/strict";
import { after, test } from "node:test";
import { startSsr } from "./ssr-helper.mjs";

const ssr = await startSsr("steps");
const { StepList, StepView } = await ssr.load("/src/components/StepList.tsx");
after(() => ssr.close());

const table = { id: "formula", type: "table", title: "Formula detail", columns: ["Formula"], data: { rows: [["C8H10N4O2"]] } };
const step = (id, title, status, extra = {}) => ({ id, title, status, seconds: 1.234, note: "", blocks: [], ...extra });
const steps = [
  step(1, "Adduct", "done", { note: "Scored 7 adducts; kept [M+H]+ (0.56).", blocks: [table] }),
  step(2, "Formula", "done", { blocks: [table], note: "BUDDY with the MS1 isotope pattern: 4 formulas for [M+H]+ and 3 formulas for [M+NH4]+ in the first pass of the search." }),
  step(3, "Library search", "running", { seconds: 0 }),
  step(4, "MS-FINDER structures", "skipped", { note: "Library gate passed.", seconds: 0 }),
  step(5, "FRIGID structures", "failed", { note: "FRIGID stopped: the engine returned no candidates for this formula.", blocks: [table] }),
  step(6, "Merge and re-rank", "queued", { seconds: 0 }),
];
const view = list => ssr.render(StepList, { steps: list });

test("no step is open by default, whatever its status or content", () => {
  const html = view(steps);
  assert.equal((html.match(/<details /g) ?? []).length, 3, "the three steps with blocks are details");
  assert.doesNotMatch(html, /<details[^>]*\sopen/);
  assert.doesNotMatch(html, /\sopen=""/);
});

test("each step is one line: number, title, status, seconds, note", () => {
  const html = view(steps);
  const first = /<details id="workflow-step-1".*?<\/summary>/s.exec(html)[0];
  assert.match(first, /<span class="step-title">1\. Adduct<\/span><span class="badge tone-good">done<\/span><span class="step-seconds">1\.2 s<\/span><span class="step-note" title="Scored 7 adducts; kept \[M\+H\]\+ \(0\.56\)\.">Scored 7 adducts/);
  assert.match(html, /<div id="workflow-step-4" class="stage step-line is-static" data-status="skipped"><span class="step-title">4\. MS-FINDER structures<\/span><span class="badge tone-neutral">skipped<\/span><span class="step-note"[^>]*>Library gate passed\.<\/span>/);
  assert.doesNotMatch(html.slice(html.indexOf('id="workflow-step-4"'), html.indexOf('id="workflow-step-5"')), /step-seconds/, "no seconds for a skipped step");
  assert.doesNotMatch(html.slice(html.indexOf('id="workflow-step-6"')), /step-seconds/, "no seconds for a queued step");
});

test("the blocks of a step are inside it and shown only once it is opened", () => {
  const html = view(steps);
  const body = /<details id="workflow-step-1".*?<\/details>/s.exec(html)[0];
  assert.match(body, /<\/summary><div class="step-body">/);
  assert.ok(body.includes("C8H10N4O2"));
  assert.ok(!/<summary>.*C8H10N4O2.*<\/summary>/s.test(body), "the table is not in the folded line");
});

test("a long note is cut to one line in the folded step and written out in full inside it", () => {
  const html = view(steps);
  const long = steps[1].note;
  const second = /<details id="workflow-step-2".*?<\/details>/s.exec(html)[0];
  assert.match(second, /<p class="step-full-note">BUDDY with the MS1 isotope pattern/);
  assert.ok(second.includes(`title="${long}"`), "the whole note is the tooltip of the folded line");
  assert.doesNotMatch(/<details id="workflow-step-1".*?<\/details>/s.exec(html)[0], /step-full-note/, "a short note is not repeated");
});

test("a failed step shows its note in the folded line", () => {
  const html = view(steps);
  const failed = /<details id="workflow-step-5".*?<\/summary>/s.exec(html)[0];
  assert.match(failed, /<span class="badge tone-bad">failed<\/span>/);
  assert.match(failed, /<span class="step-note is-failed"[^>]*>FRIGID stopped: the engine returned no candidates for this formula\.<\/span>/);
  assert.doesNotMatch(failed, /<details[^>]*\sopen/);
});

test("steps with nothing to open are plain lines, and every step keeps its anchor id", () => {
  const html = view(steps);
  for (const s of steps) assert.ok(html.includes(`id="workflow-step-${s.id}"`));
  assert.match(html, /<div id="workflow-step-3" class="stage step-line is-static" data-status="running">/);
});

test("missing or odd step data does not break the line", () => {
  const odd = { id: 1, title: "", status: null, seconds: null, note: undefined, blocks: undefined };
  const html = ssr.render(StepView, { step: odd, index: 0 });
  assert.match(html, /1\. Unnamed stage/);
  assert.match(html, />unknown</);
  assert.doesNotMatch(html, /step-seconds/);
  assert.equal(view(undefined), "");
});
