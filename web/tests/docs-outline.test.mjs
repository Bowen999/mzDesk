import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { after, test } from "node:test";
import { startSsr } from "./ssr-helper.mjs";

// The Docs contents (docs/outline.ts) and its tie to the backend list of general pages (GENERAL_PAGES in mzdesk/api/app.py).
const ssr = await startSsr("docs-outline");
const outline = await ssr.load("/src/docs/outline.ts");
after(() => ssr.close());

const spec = {
  groups: [], version: "x",
  subgroups: [{ id: "molecular", title: "Compound" }, { id: "ms2", title: "Spectrum" }, { id: "empty", title: "Empty" }],
  pipelines: [{ id: "esi_ms2_annotation", title: "MS/MS Annotation" }, { id: "ei_ms_annotation", title: "EI MS Annotation" }],
  modules: [
    { id: "identity", title: "Identity", group: "modules", subgroup: "molecular" },
    { id: "ms2_library", title: "Library Search", group: "modules", subgroup: "ms2" },
    { id: "calculator", title: "Calculator", group: "modules", subgroup: null },
  ],
};
const titles = entries => entries.map(e => `${e.id}: ${e.title}`);

test("the contents are Getting Started, Concepts, Workflows, Modules, Reference, with these general pages", () => {
  const o = outline.docsOutline(spec);
  assert.deepEqual(o.map(s => s.title), ["Getting Started", "Concepts", "Workflows", "Modules", "Reference"]);
  assert.deepEqual(titles(o[0].entries), ["overview: Overview", "quick-start: Quick Start", "local-version: Local Version"]);
  assert.deepEqual(titles(o[1].entries), ["workflows-and-modules: Workflows and Modules", "inputs: Inputs and Formats", "reading-results: Reading Results",
    "session-report: Session Report", "privacy: Privacy and Internet Search"]);
  assert.deepEqual(titles(o[4].entries), ["cli: Command Line and Settings", "faq: FAQ", "changelog: Changelog", "licenses: Licenses and Based On", "feedback: Feedback", "contact: Contact Us"]);
});

test("workflows and modules follow the spec: subgroup folders in order, empty ones left out, a module without a subgroup at the end", () => {
  const o = outline.docsOutline(spec);
  assert.deepEqual(titles(o[2].entries), ["esi_ms2_annotation: MS/MS Annotation", "ei_ms_annotation: EI MS Annotation"]);
  assert.ok(o[2].entries.every(e => e.kind === "workflow"));
  const [compound, spectrum, calc] = o[3].entries;
  assert.deepEqual([compound.title, spectrum.title], ["Compound", "Spectrum"]);
  assert.deepEqual(titles(compound.items), ["identity: Identity"]);
  assert.deepEqual(titles(spectrum.items), ["ms2_library: Library Search"]);
  assert.deepEqual([calc.id, calc.kind], ["calculator", "module"]);
  assert.equal(o[3].entries.length, 3);
});

test("every page has a unique id, and Previous / Next walk the pages in reading order", () => {
  const pages = outline.flatPages(outline.docsOutline(spec));
  assert.equal(new Set(pages.map(p => p.id)).size, pages.length);
  assert.equal(pages[0].id, "overview");
  assert.equal(pages.at(-1).id, "contact");
  assert.deepEqual(pages.slice(3, 8).map(p => p.id), ["workflows-and-modules", "inputs", "reading-results", "session-report", "privacy"]);
  assert.equal(pages.findIndex(p => p.id === "esi_ms2_annotation"), 8, "workflows start right after Concepts");
});

test("the general pages of the contents are exactly GENERAL_PAGES of the backend", () => {
  const app = readFileSync(new URL("../../mzdesk/api/app.py", import.meta.url), "utf8");
  const block = /^GENERAL_PAGES = \(([\s\S]*?)\n\)/m.exec(app);
  assert.ok(block, "GENERAL_PAGES is defined in app.py");
  const backend = [...block[1].matchAll(/"([a-z0-9-]+)"/g)].map(m => m[1]);
  const frontend = outline.GENERAL_PAGES.map(p => p.id);
  assert.deepEqual([...backend].sort(), [...frontend].sort());
  const shown = outline.docsOutline(spec).flatMap(s => s.entries).filter(e => e.kind === "page").map(e => e.id);
  assert.deepEqual([...shown].sort(), [...frontend].sort(), "every general page is listed once in the contents");
});
