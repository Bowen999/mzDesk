import assert from "node:assert/strict";
import { mkdtempSync, rmSync } from "node:fs";
import { tmpdir } from "node:os";
import { dirname, join } from "node:path";
import { after, test } from "node:test";
import { fileURLToPath } from "node:url";
import React from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { createServer } from "vite";

// Figure's theme module initializes the document theme on import.
globalThis.document = { documentElement: { removeAttribute() {}, setAttribute() {} } };
const root = dirname(dirname(fileURLToPath(import.meta.url)));
const cache = mkdtempSync(join(tmpdir(), "mzdesk-spectrum-test-"));
const server = await createServer({ configFile: false, root, cacheDir: cache, server: { middlewareMode: true, hmr: false }, appType: "custom" });
const { SpectrumEvidence } = await server.ssrLoadModule("/src/components/PipelineVisuals.tsx");
after(async () => { await server.close(); rmSync(cache, { recursive: true, force: true }); delete globalThis.document; });

const query = { img: "query-spectrum.png" };
const library = { figure: { img: "mirror-spectrum.png" }, library: "MassBank", accession: "MB123", name: "Caffeine" };
const result = (status, note, { libraries = ["MassBank"], visualization = { query, library: null }, steps } = {}) => ({
  context: { precursor_mz: 195.0877, polarity: "Positive", libraries },
  steps: steps ?? [{ id: 1, status: "done" }, { id: 3, status, note }, { id: 4, status: "done" }],
  visualization,
});
const view = data => renderToStaticMarkup(React.createElement(SpectrumEvidence, { result: data }));

test("a generated query spectrum does not show evidence when library search was skipped", () => {
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

test("completed searches with no hits retain the backend's actual reason", () => {
  for (const note of ["No library spectrum within 10 ppm of m/z 195.0877.", "8 library spectra compared; none scored 0.70 or higher."]) {
    const html = view(result("done", note));
    assert.match(html, /Spectrum evidence/);
    assert.match(html, /Library search: done/);
    assert.ok(html.includes(note));
    assert.match(html, /query-spectrum\.png/);
    assert.match(html, /Library comparison figure unavailable/);
  }
});

test("failed search remains distinct from a completed search with no hits", () => {
  const note = "The library search failed (RuntimeError: index broken).";
  const html = view(result("failed", note));
  assert.match(html, /Library search: failed/);
  assert.ok(html.includes(note));
  assert.match(html, /query-spectrum\.png/);
  assert.doesNotMatch(html, /none scored|No library spectrum within/);
});

test("missing comparison figure is not mistaken for no library hits", () => {
  const note = "1 library spectra scored 0.70 or higher (1 compounds listed). Gate passed by Caffeine.";
  const html = view(result("done", note));
  assert.ok(html.includes(note));
  assert.match(html, /Library comparison figure unavailable/);
  assert.doesNotMatch(html, /no hits|none scored|No library spectrum within/i);
});

test("a judged hit displays the comparison and its query tab", () => {
  const html = view(result("done", "1 library spectra scored 0.70 or higher.", { visualization: { query, library } }));
  assert.match(html, /Library search: done/);
  assert.match(html, /Best judged library hit: <b>Caffeine<\/b> · MassBank · MB123/);
  assert.match(html, /mirror-spectrum\.png/);
  assert.match(html, /Experimental spectrum/);
  assert.match(html, /Library comparison/);
  assert.match(html, /aria-pressed="true"[^>]*>Library comparison/);
});

test("failed or completed searches keep their status even if figures failed to render", () => {
  for (const status of ["done", "failed"]) {
    for (const visualization of [null, { query: null, library: null }]) {
      const html = view(result(status, "Backend step detail.", { visualization }));
      assert.match(html, new RegExp(`Library search: ${status}`));
      assert.match(html, /Backend step detail/);
      assert.match(html, /Spectrum figures are unavailable/);
    }
  }
});

test("a comparison remains visible when the query figure is unavailable", () => {
  const html = view(result("done", "A hit was judged.", { visualization: { query: null, library } }));
  assert.match(html, /mirror-spectrum\.png/);
  assert.match(html, /Best judged library hit/);
  assert.doesNotMatch(html, /Spectrum figures are unavailable/);
});
