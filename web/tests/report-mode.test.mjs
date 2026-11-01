import assert from "node:assert/strict";
import { after, test } from "node:test";
import { React, startSsr } from "./ssr-helper.mjs";

// Session report (report/mode.ts): the pages are drawn once more, off screen, in report mode, so report.html holds the whole content.
// Outside report mode every component draws what the screen shows.
const ssr = await startSsr("report-mode");
const { MAX_REPORT_ROWS, ReportModeContext } = await ssr.load("/src/report/mode.ts");
const { CandidateTable } = await ssr.load("/src/components/CandidateTable.tsx");
const { FragmentProvider } = await ssr.load("/src/components/FragmentDialog.tsx");
const { GroupTable } = await ssr.load("/src/components/Blocks.tsx");
const { OutputBlock, OutputList } = await ssr.load("/src/components/Outputs.tsx");
const { ReportPage } = await ssr.load("/src/report/ReportPage.tsx");
const { ModulePage } = await ssr.load("/src/pages/ModulePage.tsx");
after(() => ssr.close());

const h = React.createElement;
const mode = (maxRows = MAX_REPORT_ROWS) => ({ maxRows, track() {} });
/** `component` drawn on screen (no context) or in report mode. */
const Wrapped = ({ of, report, ...props }) => (report ? h(ReportModeContext.Provider, { value: report }, h(of, props)) : h(of, props));
const screen = (of, props) => ssr.render(Wrapped, { of, ...props });
const report = (of, props, m = mode()) => ssr.render(Wrapped, { of, report: m, ...props });
const text = html => html.replace(/<[^>]*>/g, " ").replace(/\s+/g, " ");
const count = (html, re) => (html.match(re) ?? []).length;
const heads = html => [...html.matchAll(/<th[ >].*?<\/th>/g)].map(m => m[0].replace(/<[^>]*>/g, ""));

const SMILES = "CN1C=NC2=C1C(=O)N(C)C(=O)N2C";
const candidate = (id, name, extra = {}) => ({
  id, rank: 1, name, smiles: SMILES, formula: "C8H10N4O2", adduct: "[M+H]+", layer: "library", layer_label: "library", level: "2a",
  level_label: "Level 2a: probable structure (library spectrum match)", level_reason: "Library spectrum passed the Level-2 gate.", final_score: 0.8, base_score: 1,
  hypothesis: false, evidence: { library: { gate_passed: true } }, detect: { omics: "metabolome", ionizable: "moderate", adduct: "high" }, fragment_adduct: "[M+H]+",
  isomer: { verdict: "indistinguishable", role: "top", challenger: false, explained_fraction: 1 }, ...extra,
});
const candidates = [
  candidate("c1", "Caffeine"),
  candidate("c2", "Theophylline", { rank: 2, final_score: 0.7, smiles: "Cn1c(=O)c2[nH]cnc2n(C)c1=O", adduct: "[M+Na]+", formula: "C7H8N4O2" }),
];
const query = { peaks: [[69.04, 18]], precursor_mz: 195.0877, tol_da: 0.02 };

// ---------------------------------------------------------------- ESI candidates
test("report mode: every candidate's SMILES, formula, adduct and layer are in the HTML, with no toggle or View column", () => {
  const view = props => h(FragmentProvider, { query }, h(CandidateTable, props));
  const closed = screen(view, { candidates });
  for (const t of [SMILES, "C8H10N4O2", "[M+Na]+", "Detectability"]) assert.ok(!closed.includes(t), `${t} is folded away on screen`);

  const html = report(view, { candidates });
  for (const t of [SMILES, "Cn1c(=O)c2[nH]cnc2n(C)c1=O", "C8H10N4O2", "C7H8N4O2", "[M+H]+", "[M+Na]+", "library", "metabolome / moderate / high", "top"]) {
    assert.ok(html.includes(t), `${t} is in the report`);
  }
  assert.deepEqual(heads(html), ["Rank", "Name", "Level", "Final score"]);
  assert.doesNotMatch(html, /cand-toggle|fv-open|View</, "no toggle or View button");
  assert.equal(count(html, /class="cand-detail"/g), 2, "one open detail line per candidate");
  assert.equal(count(html, /<td colSpan="4">/g), 2, "the detail line spans the 4 columns");
});

test("report mode: the CCS check and the no-candidate row still draw", () => {
  const withCcs = candidate("c3", "Caffeine", { ccs: { verdict: "Consistent", source: "measured", delta_pct: 0.64, tolerance_pct: 1, reference: 139.6, ref_label: "DTCCS N2" } });
  const html = report(CandidateTable, { candidates: [withCcs] });
  assert.match(text(html), /CCS Consistent \+0\.64 % \(±1\.0 %\) · measured 139\.6/);
  assert.match(report(CandidateTable, { candidates: [] }), /<td colSpan="4" class="empty">no candidates<\/td>/);
});

// ---------------------------------------------------------------- folded groups
const groups = [
  { title: "Mass search hits", meta: "2 compounds", rows: [["Caffeine", "CHEBI:27732"], ["Theobromine", "CHEBI:28946"]] },
  { title: "similar", meta: "1 compound", rows: [["Paraxanthine", "CHEBI:4482"]] },
];

test("a group table on screen opens group by group: closed groups have no rows in the HTML", () => {
  const html = screen(GroupTable, { columns: ["Name", "ID"], groups });
  assert.equal(count(html, /class="gt-group"/g), 2);
  assert.equal(count(html, /gt-toggle/g), 2);
  for (const t of ["Caffeine", "CHEBI:27732", "Paraxanthine"]) assert.ok(!html.includes(t), `${t} is not drawn while closed`);
  assert.doesNotMatch(html, /<details/);
});

test("report mode: each group is a native <details> that holds its rows, closed like the screen", () => {
  const html = report(GroupTable, { columns: ["Name", "ID"], groups });
  assert.equal(count(html, /<details class="gt-group">/g), 2);
  assert.equal(count(html, /<summary class="gt-head">/g), 2);
  assert.doesNotMatch(html, /<details[^>]*\sopen/, "closed to start with");
  for (const t of ["Mass search hits", "2 compounds", "Caffeine", "CHEBI:27732", "Theobromine", "CHEBI:28946", "similar", "Paraxanthine", "CHEBI:4482"]) {
    assert.ok(html.includes(t), `${t} is in the report`);
  }
  assert.doesNotMatch(html, /<button/, "no toggle button (a button inside <summary> would not toggle)");
  assert.doesNotMatch(report(GroupTable, { columns: ["Name", "ID"], groups, rowHelp: true }), /help-col/, "no column of ? buttons: the report has no popovers");
  // the "similar" group's rows sit inside its own <details>
  const similar = html.slice(html.lastIndexOf('<details class="gt-group">'));
  assert.ok(similar.includes("Paraxanthine") && !similar.includes("Caffeine"));
});

// ---------------------------------------------------------------- paged tables
const longUsi = i => `mzspec:GNPS:TASK-0123456789abcdef-spectra/specs_ms.mgf:scan:${i}:${"X".repeat(40)}`;
const pagedBlock = { id: "hits", type: "table", title: "GNPS hits", columns: ["Hit", "USI"], page_size: 25, compact: true };
const rowsOf = n => Array.from({ length: n }, (_, i) => [`Hit ${i + 1}`, longUsi(i + 1)]);
const table = (rows, extra = {}) => ({ block: pagedBlock, data: { rows, ...extra }, status: "none" });

test("a paged table on screen shows one page and its controls", () => {
  const html = screen(OutputBlock, table(rowsOf(60)));
  assert.equal(count(html, /<td/g) / 2, 25, "25 rows of the 60");
  assert.match(html, /result-table-pages/);
  assert.match(html, /type="search"/);
  assert.match(html, /Page 1 of 3/);
  assert.ok(!html.includes("Hit 26<") && !html.includes(">Hit 60<"));
});

test("report mode: a paged table lists every row, with no pages and no filter", () => {
  const html = report(OutputBlock, table(rowsOf(60)));
  assert.equal(count(html, /<td/g) / 2, 60, "all 60 rows");
  assert.ok(html.includes(">Hit 1<") || html.includes("Hit 1</"));
  assert.ok(html.includes("Hit 60"));
  assert.doesNotMatch(html, /result-table-pages|type="search"|page-button|Filter preview/);
  assert.doesNotMatch(html, /First \d+ of/, "nothing was cut");
});

test("report mode: a very long table keeps the first rows and says where the rest is", () => {
  const html = report(OutputBlock, table(rowsOf(7)), mode(3));
  assert.equal(count(html, /<td/g) / 2, 3);
  assert.ok(html.includes("Hit 3") && !html.includes("Hit 4"));
  assert.match(html, /<p class="result-table-note">First 3 of 7 rows; all rows are in report\.json\.<\/p>/);
  const big = report(OutputBlock, table(rowsOf(MAX_REPORT_ROWS + 100)));
  assert.equal(count(big, /<td/g) / 2, MAX_REPORT_ROWS);
  assert.match(big, new RegExp(`First ${MAX_REPORT_ROWS} of ${MAX_REPORT_ROWS + 100} rows; all rows are in report\\.json\\.`));
  assert.equal(MAX_REPORT_ROWS, 500);
});

test("report mode: a table without pages and the rows of a group are capped the same way", () => {
  const plain = { id: "t", type: "table", title: "Plain", columns: ["Hit", "USI"] };
  const html = report(OutputBlock, { block: plain, data: { rows: rowsOf(5) }, status: "none" }, mode(2));
  assert.match(html, /First 2 of 5 rows; all rows are in report\.json\./);
  const grouped = report(GroupTable, { columns: ["Name", "ID"], groups: [{ title: "g", rows: Array.from({ length: 4 }, (_, i) => [`n${i}`, `id${i}`]) }] }, mode(2));
  assert.match(grouped, /First 2 of 4 rows/);
});

test("the server's preview limit is still stated in report mode, without the filter sentence", () => {
  const html = report(OutputBlock, table(rowsOf(30), { count: 1200 }));
  assert.match(html, /Preview: first 30 of 1,200 results\./);
  assert.doesNotMatch(html, /Filtering only searches/);
});

test("report mode: a column of View buttons is left out of a table", () => {
  const view = { fragment: { smiles: SMILES, adduct: "[M+H]+" } };
  const block = { id: "f", type: "table", title: "Hits", columns: ["Name", "Fragments", "Score"] };
  const rows = [["Caffeine", view, "0.9"], ["Theobromine", "—", "0.8"]];
  assert.match(screen(OutputBlock, { block, data: { rows }, status: "none" }), /<th>Fragments<\/th>/);
  const html = report(OutputBlock, { block, data: { rows }, status: "none" });
  assert.deepEqual(heads(html), ["Name", "Score"]);
  assert.ok(html.includes("0.9") && html.includes("Theobromine"));
});

// ---------------------------------------------------------------- long text
test("CompactText: cut at 80 characters on screen, whole in report mode", () => {
  const full = "A".repeat(60) + "B".repeat(60) + "END";
  const compact = { id: "c", type: "table", title: "Long", columns: ["Name", "Text"], compact: true };
  const props = { block: compact, data: { rows: [["x", full]] }, status: "none" };
  const on = screen(OutputBlock, props);
  assert.ok(!on.includes("END") && on.includes("…") && on.includes('class="cell-expand"'), "cut, with a More button");
  const html = report(OutputBlock, props);
  assert.ok(html.includes(full), "the whole text");
  assert.doesNotMatch(html, /cell-expand|compact-cell|More/);
});

// ---------------------------------------------------------------- other folded or tabbed content
test("report mode: all formats of every compound, a titled collapsed note, all reaction options, every functional group; no Download buttons", () => {
  const formats = { id: "fm", type: "formats", title: "Formats", label: "Formats" };
  const data = { sets: [{ label: "Caffeine", formats: { SMILES: SMILES, InChI: "InChI=1S/C8H10N4O2" } }, { label: "Aspirin", formats: { SMILES: "CC(=O)Oc1ccccc1C(=O)O", InChI: "InChI=1S/C9H8O4" } }] };
  const one = screen(OutputBlock, { block: formats, data, status: "none" });
  assert.match(one, /role="tablist"/);
  assert.ok(!one.includes("InChI=1S/C9H8O4"), "one tab at a time on screen");
  const html = report(OutputBlock, { block: formats, data, status: "none" });
  for (const t of [SMILES, "InChI=1S/C8H10N4O2", "CC(=O)Oc1ccccc1C(=O)O", "InChI=1S/C9H8O4", "Caffeine", "Aspirin"]) assert.ok(html.includes(t), `${t} is in the report`);
  assert.doesNotMatch(html, /role="tab|class="copy"/);

  const notes = { id: "n", type: "notes", title: "Method notes", collapsed: true };
  assert.doesNotMatch(screen(OutputBlock, { block: notes, data: { items: ["Elutes early."] }, status: "none" }), /Elutes early/, "on screen the items live in the ?");
  const note = report(OutputBlock, { block: notes, data: { items: ["Elutes early."] }, status: "none" });
  assert.match(note, /<details class="out-details"><summary>Method notes<\/summary>.*Elutes early\./);

  const reaction = { id: "rx", type: "reaction", title: "Reactions" };
  const options = [{ label: "Silylation", steps: [{ from: "R-OH", reagent: "BSTFA", to: "R-OTMS" }] }, { label: "Methylation", steps: [{ from: "R-COOH", reagent: "CH2N2", to: "R-COOMe" }] }];
  assert.ok(!screen(OutputBlock, { block: reaction, data: { options }, status: "none" }).includes("R-COOMe"), "one option at a time on screen");
  const rx = report(OutputBlock, { block: reaction, data: { options }, status: "none" });
  for (const t of ["Silylation", "R-OTMS", "Methylation", "R-COOMe", "BSTFA", "CH2N2"]) assert.ok(rx.includes(t), `${t} is in the report`);

  const groupmap = { id: "gm", type: "groupmap", title: "Groups" };
  const found = [
    { family: "Alcohols", name: "Hydroxyl", count: 2, smarts: "[OX2H]", atoms: "1, 2", svg: "<svg/>", all_svg: "<svg/>", sites: [{ atoms: "1", svg: "<svg/>" }, { atoms: "2", svg: "<svg/>" }] },
    { family: "Amines", name: "Primary amine", count: 1, smarts: "[NX3H2]", atoms: "3", svg: "<svg/>", sites: [{ atoms: "3", svg: "<svg/>" }] },
  ];
  const onScreen = screen(OutputBlock, { block: groupmap, data: { svg: "<svg/>", found }, status: "none" });
  assert.ok(onScreen.includes("[OX2H]") && !onScreen.includes("[NX3H2]"), "one group at a time on screen");
  const gm = report(OutputBlock, { block: groupmap, data: { svg: "<svg/>", found }, status: "none" });
  for (const t of ["Hydroxyl", "[OX2H]", "Primary amine", "[NX3H2]"]) assert.ok(gm.includes(t), `${t} is in the report`);
  assert.doesNotMatch(gm, /<button|gm-sites/, "no group or site buttons");

  const blocks = [{ id: "dl", type: "downloads", title: "Files" }, { id: "t", type: "table", title: "T", columns: ["A"] }];
  const outputs = { dl: { files: [{ name: "hits.csv", content: "a,b" }] }, t: { rows: [["1"]] } };
  assert.match(screen(OutputList, { blocks, outputs, live: true }), /Download hits\.csv/);
  const list = report(OutputList, { blocks, outputs, live: true });
  assert.doesNotMatch(list, /Download|hits\.csv/);
  assert.match(list, /<table/);
});

// ---------------------------------------------------------------- whole pages
const block = (id, type, extra = {}) => ({ id, type, title: id, ...extra });
const spec = {
  groups: [], subgroups: [{ id: "ms1", title: "MS1" }], version: "0.0.0",
  modules: [{
    id: "mass_search", title: "Mass search", cls: "x", group: "g", subgroup: "ms1", live: true, one_liner: "Find compounds by mass.", inputs: [], example: null,
    outputs: [block("hits", "grouptable", { columns: ["Name", "ID"] }), block("gnps", "table", { columns: ["Hit", "USI"], page_size: 25, compact: true })],
  }],
  pipelines: [{ id: "esi_ms2_annotation", title: "ESI MS/MS Annotation", one_liner: "Name a <b>compound</b>.", inputs: [], stages: [] }],
};
const stepTable = { id: "formula", type: "table", title: "Formula detail", columns: ["Formula"], data: { rows: [["C8H10N4O2"]] } };
const esi = {
  answer: { type: "single", ids: ["c1"], layer: "library", layer_label: "library", level: "2a", level_label: "Level 2a: probable structure (library spectrum match)",
    statement: "Caffeine is the top-ranked structure (library, Level 2a).", top_id: "c1", isomer_verdict: "single" },
  candidates, warnings: [], steps: [], fragment_query: query, isomer: null, visualization: null,
  context: { precursor_mz: 195.0877, polarity: "Positive", adduct: "[M+H]+", retained_adducts: [], omics: "", libraries: [], gate_passed: true },
};
const entry = (kind, id, result, extra = {}) => ({ key: `${kind}:${id}`, kind, id, title: id, section: "x", inputs: {}, result, seconds: 1.5, finishedAt: "2026-11-01T10:00:00-06:00", ...extra });

test("a module page in the report: title, results, groups and all rows of a long table", () => {
  const result = { status: "ok", warnings: ["One note."], verdict: "Found <b>2</b> hits.", outputs: { hits: { groups }, gnps: { rows: rowsOf(40) } } };
  const html = ssr.render(ReportPage, { spec, entry: entry("module", "mass_search", result), mode: mode() });
  assert.match(html, /<div class="page module-page" data-module="mass_search">/);
  assert.match(html, /<h1>Mass search<\/h1>/);
  assert.match(html, /Found <b>2<\/b> hits\./);
  assert.match(html, /<div class="alert warn">One note\.<\/div>/);
  assert.match(html, /<h2 class="panel-title">Results<\/h2><span class="msg">1\.5 s<\/span>/);
  assert.ok(html.includes("CHEBI:28946") && html.includes("Paraxanthine"), "the folded groups' compounds");
  assert.ok(html.includes("Hit 40"), "the last of 40 rows");
  assert.doesNotMatch(html, /<input|<button|class="panel"|page-button|hp-btn/, "nothing to type into or click");
});

test("the ESI workflow page in the report: answer, open candidate details, steps with their blocks", () => {
  const steps = [{ id: 1, title: "Formula", status: "done", seconds: 1.2, note: "", blocks: [stepTable] }];
  const html = ssr.render(ReportPage, { spec, entry: entry("workflow", "esi_ms2_annotation", esi, { steps, inputs: { precursor: { adduct: "[M+H]+" } } }), mode: mode() });
  assert.match(html, /data-pipeline="esi_ms2_annotation"/);
  assert.match(html, /<p class="one-liner">Name a <b>compound<\/b>\.<\/p>/);
  assert.match(html, /answer-headline/);
  for (const t of [SMILES, "C8H10N4O2", "[M+Na]+", "Formula detail"]) assert.ok(html.includes(t), `${t} is in the page`);
  assert.match(html, /<details id="workflow-step-1" class="stage step-line" data-status="done">/);
  assert.ok(html.indexOf("C8H10N4O2", html.indexOf("workflow-step-1")) > 0, "the step's table is inside its <details>");
  assert.doesNotMatch(html, /cand-toggle|fv-open|<input/);
  assert.match(html, /<span class="msg">2 structure\(s\)<\/span>/);
});

test("a generic workflow page in the report keeps its answer and steps", () => {
  const generic = { pipeline: "targeted_method", answer: { statement: "Fits RP LC-MS.", level: "A", level_label: "A: strong", facts: [["Platform", "RP"]] }, blocks: [], warnings: [], steps: [] };
  const sp = { ...spec, pipelines: [{ ...spec.pipelines[0], id: "targeted_method", title: "Targeted Method" }] };
  const html = ssr.render(ReportPage, { spec: sp, entry: entry("workflow", "targeted_method", generic, { steps: [] }), mode: mode() });
  assert.match(html, /Fits RP LC-MS\./);
  assert.match(html, /<dt>Platform<\/dt><dd>RP<\/dd>/);
  assert.match(html, /<h2 class="panel-title">Steps<\/h2>/);
});

test("a page the spec no longer has draws nothing", () => {
  assert.equal(ssr.render(ReportPage, { spec, entry: entry("module", "gone", { status: "ok", outputs: {}, warnings: [] }), mode: mode() }), "");
});

// ---------------------------------------------------------------- the screen is not changed
test("outside report mode a module page still draws the screen version (inputs, closed groups, paged table)", () => {
  const html = ssr.render(ModulePage, { id: "mass_search", spec, active: true });
  assert.match(html, /<div class="page module-page" data-module="mass_search">/);
  assert.match(html, /<div class="panel">/, "the inputs panel is there");
  assert.match(html, /Results appear here after you run this module\./);
});
