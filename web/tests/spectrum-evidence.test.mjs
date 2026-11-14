import assert from "node:assert/strict";
import { after, test } from "node:test";
import { startSsr } from "./ssr-helper.mjs";

const ssr = await startSsr("spectrum");
const { SpectrumEvidence } = await ssr.load("/src/components/PipelineVisuals.tsx");
after(() => ssr.close());

const query = { img: "query-spectrum.png" };
const library = { figure: { img: "mirror-spectrum.png" }, library: "MassBank", accession: "MB123", name: "Caffeine" };
const result = (status, note, { libraries = ["MassBank"], visualization = { query, library: null }, steps } = {}) => ({
  context: { precursor_mz: 195.0877, polarity: "Positive", libraries },
  steps: steps ?? [{ id: 1, status: "done" }, { id: 3, status, note }, { id: 4, status: "done" }],
  visualization,
});
const view = data => ssr.render(SpectrumEvidence, { result: data });

test("no library search, no evidence: nothing is drawn when the search was skipped", () => {
  assert.equal(view(result("skipped", "No library was selected, so no library spectrum was searched.", { libraries: [] })), "");
  assert.equal(view(result("skipped", "MassBank is not installed.", { visualization: { query, library } })), "");
  assert.equal(view(result("skipped", "The MS2 spectrum has no peaks, so no library spectrum was searched.")), "");
});

test("missing or unfinished step 3 does not invent evidence from selected libraries", () => {
  assert.equal(view(result("queued", "")), "");
  assert.equal(view(result("running", "")), "");
  assert.equal(view(result("done", "", { steps: [{ id: 1, status: "done" }] })), "");
  assert.equal(view({ ...result("done", ""), steps: undefined }), "");
});

test("a finished search without a library hit to compare draws nothing; the reason stays in the step's note", () => {
  for (const note of ["No library spectrum within 10 ppm of m/z 195.0877.", "8 library spectra compared; none scored 0.70 or higher."]) {
    assert.equal(view(result("done", note)), "");
    assert.equal(view(result("done", note, { visualization: { query: null, library: null } })), "");
    assert.equal(view(result("done", note, { visualization: null })), "");
  }
});

test("a failed search draws nothing", () => {
  assert.equal(view(result("failed", "The library search failed (RuntimeError: index broken).")), "");
  assert.equal(view(result("failed", "failed", { visualization: { query, library } })), "");
});

test("a judged hit displays the comparison, its query tab and a short hit line", () => {
  const html = view(result("done", "1 library spectra scored 0.70 or higher.", { visualization: { query, library } }));
  assert.match(html, /Spectrum evidence/);
  assert.match(html, /Best library hit: <b>Caffeine<\/b> · MassBank · MB123/);
  assert.match(html, /mirror-spectrum\.png/);
  assert.match(html, /Experimental spectrum/);
  assert.match(html, /Library comparison/);
  assert.match(html, /aria-pressed="true"[^>]*>Library comparison/);
});

test("the explanations live behind a help button, not in paragraphs of the section", () => {
  const html = view(result("done", "107 library spectra scored 0.70 or higher.", { visualization: { query, library } }));
  assert.match(html, /<button[^>]*class="hp-btn"[^>]*aria-label="Help: Spectrum evidence"/);
  for (const text of ["Measured m/z", "greedy matching", "may not be the top-ranked", "Library search: done", "Parsed experimental spectrum", "Library comparison figure unavailable"]) {
    assert.ok(!html.includes(text), `${text} is not in the visible markup`);
  }
  assert.doesNotMatch(html, /<p class="workflow-caption">/);
  assert.equal((html.match(/<p /g) ?? []).length, 1, "only the hit line is a paragraph");
});

test("a comparison remains visible when the query figure is unavailable, without tabs", () => {
  const html = view(result("done", "A hit was judged.", { visualization: { query: null, library } }));
  assert.match(html, /mirror-spectrum\.png/);
  assert.match(html, /Best library hit/);
  assert.doesNotMatch(html, /workflow-tabs/);
});
