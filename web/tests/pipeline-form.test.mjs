import assert from "node:assert/strict";
import { after, test } from "node:test";
import { startSsr } from "./ssr-helper.mjs";

const ssr = await startSsr("pipeline-form");
const { PipelinePage } = await ssr.load("/src/pages/PipelinePage.tsx");
const { initialValues, missingRequired, visibleFields } = await ssr.load("/src/components/Fields.tsx");
after(() => ssr.close());

// The form of ESI MS/MS Annotation as the spec sends it (mzdesk/spec.py): a section can start collapsed, and the platform select starts empty.
const PLATFORMS = ["DTIMS (stepped field)", "DTIMS (single field)", "TIMS", "TWIMS", "FAIMS / DMS"];
const inputs = [
  { type: "number", id: "tol_ppm", label: "Precursor tolerance (ppm)", value: 10, section: "Search" },
  { type: "checkbox", id: "halogen", label: "Allow Cl/Br in formulas", value: false, section: "Options" },
  { type: "number", id: "ccs_observed", label: "Observed CCS (Å²)", value: null, placeholder: "e.g. 152.3", section: "Ion Mobility", section_collapsed: true },
  { type: "select", id: "ccs_platform", label: "Platform", options: ["", ...PLATFORMS], value: "", placeholder: "Select the platform", section: "Ion Mobility", section_collapsed: true },
  { type: "select", id: "ccs_gas", label: "Drift gas", options: ["N₂", "He"], value: "N₂", section: "Ion Mobility", section_collapsed: true },
  { type: "select", id: "ccs_tw_calibration", label: "TWIMS calibration", options: ["N₂ reference values", "He reference values (N₂→He)"], value: "N₂ reference values",
    visible_when: { field: "ccs_platform", equals: "TWIMS" }, section: "Ion Mobility", section_collapsed: true },
  { type: "number", id: "ph", label: "Mobile phase pH", value: 2.7, section: "Advanced" },
];
const spec = { groups: [], subgroups: [], modules: [], pipelines: [{ id: "esi_ms2_annotation", title: "ESI MS/MS Annotation", one_liner: "x", inputs, stages: [], examples: [] }] };
const page = () => ssr.render(PipelinePage, { id: "esi_ms2_annotation", spec, active: true });

test("Ion Mobility and Advanced start closed behind their headings; Search and Options are plain groups without a heading", () => {
  const html = page();
  assert.match(html, /<details class="form-section"><summary>Ion Mobility<\/summary>/);
  assert.match(html, /<details class="form-section"><summary>Advanced<\/summary>/);
  assert.equal((html.match(/<div class="form-section"/g) ?? []).length, 2, "Search and Options");
  assert.doesNotMatch(html, /<legend|<fieldset/);
  assert.doesNotMatch(html, />Search<|>Options<|>Inputs</, "no section captions");
  assert.equal((html.match(/<details class="form-section"/g) ?? []).length, 2);
  assert.doesNotMatch(html, /<details class="form-section" open/);
});

test("the platform select starts on an empty choice with a prompt, and the TWIMS calibration is hidden until TWIMS is chosen", () => {
  const html = page();
  const ion = html.slice(html.indexOf("<summary>Ion Mobility"), html.indexOf("<summary>Advanced"));
  assert.match(ion, /<select[^>]*><option value="" selected="">Select the platform<\/option><option value="DTIMS \(stepped field\)">/);
  assert.match(ion, /Drift gas/);
  assert.doesNotMatch(ion, /TWIMS calibration/);
  const values = { ...initialValues(inputs), ccs_platform: "TWIMS" };
  assert.deepEqual(visibleFields(inputs, values).map(f => f.id), ["tol_ppm", "halogen", "ccs_observed", "ccs_platform", "ccs_gas", "ccs_tw_calibration", "ph"]);
  assert.ok(!visibleFields(inputs, initialValues(inputs)).some(f => f.id === "ccs_tw_calibration"));
});

test("nothing in the Ion Mobility section is required on the page (the platform is checked by the server when a CCS is given)", () => {
  assert.deepEqual(missingRequired(inputs, initialValues(inputs)), []);
  assert.doesNotMatch(page(), /aria-required|field-label-text">[^<]*<span class="req"/);
});

test("the observed CCS is an empty number input, a section value that is not sent as a number until typed", () => {
  const init = initialValues(inputs);
  assert.deepEqual([init.ccs_observed, init.ccs_platform, init.ccs_gas], ["", "", "N₂"]);
});

// Lipid Fragment Explanation: the structure field offers a lipid short name first (spec.py LIPID_NAME_MODE), then Draw and SMILES.
const lipidInputs = [
  { type: "structure", id: "lipid", label: "Lipid", modes: ["Lipid name", "Draw", "SMILES / InChI"], network_modes: [], required: true, placeholder: "PC 16:0/18:1", section: "Lipid" },
];
const lipidSpec = { groups: [], subgroups: [], modules: [], pipelines: [{ id: "lipid_fragments", title: "Lipid Fragment Explanation", one_liner: "x", inputs: lipidInputs, stages: [], examples: [] }] };

test("a structure field with modes starts on its first mode and shows only those modes, in order", () => {
  assert.deepEqual(initialValues(lipidInputs).lipid, { mode: "Lipid name", value: "", label: "", buffers: {} });
  const html = ssr.render(PipelinePage, { id: "lipid_fragments", spec: lipidSpec, active: true });
  const labels = [...html.matchAll(/<input type="radio"[^>]*\/>\s*([^<]+)<\/label>/g)].map(m => m[1].trim());
  assert.deepEqual(labels, ["Lipid name", "Draw", "SMILES / InChI"]);
  assert.match(html, /<input type="radio"[^>]*checked=""[^>]*\/>\s*Lipid name/);
  assert.match(html, /placeholder="Lipid short name, e.g. PC 16:0\/18:1"/);
  assert.doesNotMatch(html, /Name \/ ID/);
  assert.deepEqual(missingRequired(lipidInputs, initialValues(lipidInputs)).map(f => f.id), ["lipid"]);
});

test("a structure field without modes keeps Name / ID, Draw and SMILES", () => {
  assert.equal(initialValues([{ type: "structure", id: "s" }]).s.mode, "Name / ID / Formula");
});
