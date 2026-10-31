import assert from "node:assert/strict";
import { after, test } from "node:test";
import { startSsr } from "./ssr-helper.mjs";

const ssr = await startSsr("workflow-result");
const { WorkflowResultView } = await ssr.load("/src/pages/PipelinePage.tsx");
after(() => ssr.close());

// The generic result of the other two workflows (EI MS Annotation, Targeted Method): answer, warnings, headline blocks; its steps are folded like ESI's.
const result = {
  pipeline: "targeted_method",
  answer: { statement: "Caffeine fits a reversed-phase LC-MS method.", level: "A", level_label: "A: strong", facts: [["Platform", "RP LC-MS"], ["Ion mode", "ESI+"]] },
  blocks: [{ id: "summary", type: "table", title: "Summary", columns: ["Item"], data: { rows: [["Caffeine"]] } }],
  warnings: ["One input was assumed."],
  steps: [],
};

test("the generic result keeps its answer, level, facts, warnings and blocks", () => {
  const html = ssr.render(WorkflowResultView, { result });
  assert.match(html, /<div class="verdict" aria-live="polite"><div class="verdict-label">Answer<\/div><p class="verdict-text">Caffeine fits a reversed-phase LC-MS method\.<\/p><\/div>/);
  assert.match(html, /<div class="pill-row"><span class="badge" title="A: strong">A<\/span><\/div>/);
  assert.match(html, /<dl class="answer-facts"><div><dt>Platform<\/dt><dd>RP LC-MS<\/dd><\/div><div><dt>Ion mode<\/dt><dd>ESI\+<\/dd><\/div><\/dl>/);
  assert.match(html, /<div class="alert warn">One input was assumed\.<\/div>/);
  assert.match(html, /Summary/);
  assert.ok(html.includes(">Caffeine<") || html.includes("Caffeine</td>"));
  assert.doesNotMatch(html, /answer-headline|answer-details/, "the ESI card is not used here");
});
