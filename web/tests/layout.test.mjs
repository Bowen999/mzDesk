import assert from "node:assert/strict";
import { after, test } from "node:test";
import { React, startSsr } from "./ssr-helper.mjs";

// The layout pieces of the 0.4 redesign: the "have → get" line, the page foot, lists of short values as chips, long tables that fold,
// cards that only repeat the verdict, the Docs markdown and the shared progress curve.
const ssr = await startSsr("layout");
const { Flow, PageFoot, PageHead } = await ssr.load("/src/components/PageHead.tsx");
const { ChipList, isChipList } = await ssr.load("/src/components/Blocks.tsx");
const { OutputBlock, OutputList, blockLayout, coveredCards } = await ssr.load("/src/components/Outputs.tsx");
const { Markdown } = await ssr.load("/src/components/Markdown.tsx");
const { ReportModeContext } = await ssr.load("/src/report/mode.ts");
const curve = await ssr.load("/src/progressCurve.ts");
after(() => ssr.close());

const h = React.createElement;
const count = (html, re) => (html.match(re) ?? []).length;
const report = (component, props) => ssr.render(() => h(ReportModeContext.Provider, { value: { maxRows: 500, track() {} } }, h(component, props)));

// ---------------------------------------------------------------- the description under a title
test("the description is 'what I have → what I get': the arrow splits it, the result side is the strong one", () => {
  const html = ssr.render(Flow, { text: "Compound → names, IDs, classification" });
  assert.match(html, /<p class="one-liner flow"><span class="flow-in">Compound<\/span><span class="flow-to"><span class="flow-arrow" aria-hidden="true">→<\/span>/);
  assert.match(html, /<span class="flow-out">names, IDs, classification<\/span>/);
  assert.match(html, /<span class="sr-only"> gives <\/span>/, "read aloud as 'gives'");
  assert.equal(ssr.render(Flow, { text: "Everyday calculations" }), '<p class="one-liner">Everyday calculations</p>', "a line without an arrow is plain");
});

test("a page head has the icon, the title and the line, and no kicker or subtitle", () => {
  const html = ssr.render(PageHead, { icon: "identity", title: "Identity", flow: "Compound → names" });
  assert.match(html, /<h1>Identity<\/h1>/);
  assert.doesNotMatch(html, /module-kicker|<h2|class="meta"/);
  assert.match(ssr.render(PageHead, { icon: "x", title: "T", flow: "A → B", preview: true }), /<h1>T<span class="badge tone-neutral head-badge">Preview<\/span><\/h1>/);
});

test("the page foot names what the page is built on and links to its Docs; the report has the names only", () => {
  const html = ssr.render(PageFoot, { id: "ms2_fragments", basedOn: ["MetFrag", "mzDesk"] });
  assert.match(html, /<footer class="page-foot"><span>Based on MetFrag, mzDesk<\/span><a href="\/docs\/ms2_fragments">Docs<\/a><\/footer>/);
  const file = ssr.render(PageFoot, { id: "x", basedOn: ["MetFrag"], link: false });
  assert.match(file, /Based on MetFrag/);
  assert.doesNotMatch(file, /<a /);
  assert.equal(ssr.render(PageFoot, { id: "x", link: false }), "", "nothing to say, nothing drawn");
});

// ---------------------------------------------------------------- results that are long or short
const rows = n => Array.from({ length: n }, (_, i) => [`name ${i + 1}`]);

test("a table with one column of short values is a chip list; two columns, long text or images are not", () => {
  assert.ok(isChipList(["Synonym"], rows(3)));
  assert.ok(isChipList(["Synonym"], [[{ text: "aspirin", url: "https://example.org" }]]));
  assert.ok(!isChipList(["Synonym"], []), "an empty table is not a list");
  assert.ok(!isChipList(["A", "B"], [["a", "b"]]));
  assert.ok(!isChipList(["Synonym"], [["x".repeat(65)]]), "a long value stays in a table");
  assert.ok(!isChipList(["Synonym"], [["two\nlines"]]));
  assert.ok(!isChipList(["Structure"], [[{ img: "data:image/png;base64,AA" }]]));
});

test("a short list shows every chip; a long one shows the first rows and '+N more'; two extra chips are not worth a button", () => {
  const short = ssr.render(ChipList, { rows: rows(3) });
  assert.equal(count(short, /class="chip-item"/g), 3);
  assert.doesNotMatch(short, /chip-more/);
  const long = ssr.render(ChipList, { rows: rows(50) });
  assert.equal(count(long, /class="chip-item"/g), 18);
  assert.match(long, /<button type="button" class="chip-more" aria-expanded="false">\+32 more<\/button>/);
  assert.equal(count(ssr.render(ChipList, { rows: rows(20) }), /class="chip-item"/g), 20, "18 + 2 are all shown");
  const all = report(ChipList, { rows: rows(50) });
  assert.equal(count(all, /class="chip-item"/g), 50, "the report shows all of them");
  assert.doesNotMatch(all, /chip-more/);
});

test("the table block draws a one-column table as chips and a long table folded", () => {
  const block = { id: "synonyms", type: "table", title: "Synonyms", columns: ["Synonym"] };
  const html = ssr.render(OutputBlock, { block, data: { rows: rows(4) }, status: "none" });
  assert.match(html, /<div class="chip-list">/);
  assert.doesNotMatch(html, /<table/);

  const wide = { id: "t", type: "table", title: "Hits", columns: ["Name", "Value"] };
  const data = n => ({ rows: Array.from({ length: n }, (_, i) => [`n${i + 1}`, `${i + 1}`]) });
  const folded = ssr.render(OutputBlock, { block: wide, data: data(30), status: "none" });
  assert.equal(count(folded, /<tr><td>/g), 8, "the first eight rows");
  assert.match(folded, /<button type="button" class="fold-more" aria-expanded="false">Show all 30 rows<\/button>/);
  const short = ssr.render(OutputBlock, { block: wide, data: data(14), status: "none" });
  assert.equal(count(short, /<tr><td>/g), 14, "14 rows are not folded");
  assert.doesNotMatch(short, /fold-more/);
  const file = report(OutputBlock, { block: wide, data: data(30), status: "none" });
  assert.equal(count(file, /<tr><td>/g), 30, "the report has every row");
});

test("a short chip list shares a row, a long one takes the whole width", () => {
  const block = { id: "s", type: "table", title: "Synonyms", columns: ["Synonym"] };
  assert.deepEqual(blockLayout(block, { rows: rows(5) }), { span: "half", tier: "quiet" });
  assert.deepEqual(blockLayout(block, { rows: rows(40) }), { span: "full", tier: "quiet" });
  assert.deepEqual(blockLayout({ id: "c", type: "table", columns: ["A", "B", "C", "D"] }, { rows: [["1", "2", "3", "4"]] }), { span: "full", tier: "card" });
});

// ---------------------------------------------------------------- one thing is said once
test("metric cards whose values are all in the verdict are hidden; a card with a new value stays", () => {
  const blocks = [{ id: "resolved", type: "cards", items: [{ label: "Name" }, { label: "Formula" }] }, { id: "summary", type: "cards", items: [{ label: "Formula" }, { label: "Average MW" }] }];
  const outputs = { resolved: { values: { Name: "Aspirin", Formula: "C9H8O4" } }, summary: { values: { Formula: "C9H8O4", "Average MW": "180.1590" } } };
  assert.deepEqual(coveredCards(blocks, outputs, "<b>Aspirin</b> · C9H8O4 · PubChem CID 2244"), ["resolved"]);
  assert.deepEqual(coveredCards(blocks, outputs, "nothing in common"), []);
  assert.deepEqual(coveredCards(blocks, outputs, null), []);
  assert.deepEqual(coveredCards(blocks, { resolved: { values: { Name: "—" } } }, "Aspirin"), [], "blank cards are not 'covered'");
  const html = ssr.render(OutputList, { blocks, outputs, hide: ["resolved"], live: true });
  assert.doesNotMatch(html, /Aspirin/);
  assert.match(html, /180\.1590/);
});

// ---------------------------------------------------------------- Docs markdown
test("the Docs markdown draws headings, bullets with bold and code, and paragraphs", () => {
  const html = ssr.render(Markdown, { text: "# Title\n\nSome <sub>sub</sub> text.\n\n## Inputs\n\n- **Label**: help with `code`\n- second\n\n### Block\n\n- one" });
  assert.match(html, /<h1>Title<\/h1>/);
  assert.match(html, /<p>Some <sub>sub<\/sub> text\.<\/p>/);
  assert.match(html, /<h2>Inputs<\/h2><ul><li><b>Label<\/b>: help with <code>code<\/code><\/li><li>second<\/li><\/ul>/);
  assert.match(html, /<h3>Block<\/h3><ul><li>one<\/li><\/ul>/);
  assert.doesNotMatch(ssr.render(Markdown, { text: "<script>alert(1)</script>" }), /<script/, "other tags stay literal text");
});

// ---------------------------------------------------------------- the curve every bar shares
test("the progress curve is fast at the start, slow at the end, and keeps 0 and 1 fixed", () => {
  assert.equal(curve.ease(0), 0);
  assert.equal(curve.ease(1), 1);
  assert.equal(curve.ease(-3), 0);
  assert.equal(curve.ease(7), 1);
  assert.ok(curve.ease(0.1) > 0.18 && curve.ease(0.5) === 0.75);
  const xs = [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1];
  assert.ok(xs.every((x, i) => i === 0 || curve.ease(x) > curve.ease(xs[i - 1])), "increasing");
});

test("a run of unknown length starts at once, rises with time and never reaches its cap", () => {
  assert.equal(curve.waitingPercent(0), 8);
  const t = [0, 0.5, 1, 2, 4, 8, 16, 60, 3600];
  const p = t.map(s => curve.waitingPercent(s));
  assert.ok(p.every((v, i) => i === 0 || v > p[i - 1]), `increasing: ${p}`);
  assert.ok(p.every(v => v < 92) || Math.abs(p.at(-1) - 92) < 1e-9, "never above the cap");
  assert.ok(p[1] - p[0] > 4, "visible movement in the first half second");
  assert.equal(curve.waitingPercent(-5), 8);
  assert.equal(curve.waitingPercent(1, { start: 20, cap: 50, tau: 1 }).toFixed(2), (20 + 30 * (1 - Math.exp(-1))).toFixed(2));
});

test("a reported stage percent is eased, with a floor so the bar is never empty", () => {
  assert.equal(curve.easedPercent(0), 6);
  assert.equal(curve.easedPercent(100), 100);
  assert.ok(curve.easedPercent(30) > 50, "30% of the stages is more than half of the bar");
  assert.equal(curve.easedPercent(0, 10), 10);
});
