import assert from "node:assert/strict";
import { after, test } from "node:test";
import { React, startSsr } from "./ssr-helper.mjs";

const ssr = await startSsr("candidates");
const { CandidateTable } = await ssr.load("/src/components/CandidateTable.tsx");
const { FragmentProvider } = await ssr.load("/src/components/FragmentDialog.tsx");
const { candidateDetails, ccsNumbers, ccsText } = await ssr.load("/src/components/annotationModel.ts");
after(() => ssr.close());

const SMILES = "CN1C=NC2=C1C(=O)N(C)C(=O)N2C";
const candidate = (id, name, extra = {}) => ({
  id, rank: 1, name, smiles: SMILES, formula: "C8H10N4O2", adduct: "[M+H]+", layer: "library", layer_label: "library", level: "2a",
  level_label: "Level 2a: probable structure (library spectrum match)", level_reason: "Library spectrum passed the Level-2 gate.", final_score: 0.8, base_score: 1,
  hypothesis: false, detect: { omics: "metabolome", ionizable: "moderate", adduct: "high" }, fragment_adduct: "[M+H]+",
  isomer: { verdict: "indistinguishable", role: "top", challenger: false, explained_fraction: 1 }, ...extra,
});
const candidates = [
  candidate("c1", "Caffeine"),
  candidate("c2", "Isocaffeine", { rank: 2, final_score: 0.835, smiles: "Cn1c(=O)c2ncn(C)c2n(C)c1=O", isomer: { verdict: "indistinguishable", role: "not separated", challenger: true, explained_fraction: 1 } }),
  candidate("c3", "1 M+H", { rank: 3, level: "3", level_label: "Level 3: tentative candidate", smiles: null, formula: null, fragment_adduct: null, detect: null, isomer: null, final_score: 0.8132 }),
];
const Table = props => React.createElement(FragmentProvider, { query: { peaks: [[69.04, 18]], precursor_mz: 195.0877, tol_da: 0.02 } }, React.createElement(CandidateTable, props));
const view = (props = {}) => ssr.render(Table, { candidates, ...props });
const heads = html => [...html.matchAll(/<th[^>]*>(.*?)<\/th>/g)].map(m => m[1].replace(/<[^>]*>/g, ""));

test("the visible columns are Rank, Name, Level, Final score and Fragments, with a toggle column in front", () => {
  assert.deepEqual(heads(view()), ["Details", "Rank", "Name", "Level", "Final score", "Fragments"]);
});

test("everything else is folded away: no SMILES, formula, adduct, layer, detectability or isomer role in the closed table", () => {
  const html = view();
  for (const text of ["SMILES", "Formula", "Adduct", "Layer", "Detectability", "Isomer role", "metabolome", SMILES, "C8H10N4O2", "not separated"]) {
    assert.ok(!html.includes(text), `${text} is not in the closed table`);
  }
  assert.doesNotMatch(html, /class="[^"]*\bwide\b/, "no wide, side-scrolling table");
  assert.doesNotMatch(html, /tabindex|scroll sideways|Scroll horizontally/i);
});

test("rows keep the order they came in, even when the final score does not go down the list", () => {
  const html = view();
  const names = [...html.matchAll(/<td class="cand-name">(.*?)<\/td>/g)].map(m => m[1]);
  assert.deepEqual(names, ["Caffeine", "Isocaffeine", "1 M+H"]);
  assert.deepEqual([...html.matchAll(/<td class="num">(.*?)<\/td>/g)].map(m => m[1]), ["0.800", "0.835", "0.813"]);
});

test("the level is a badge with its reason behind a ?, and the View button is where a fragment view exists", () => {
  const html = view();
  assert.match(html, /<span class="badge tone-good">Level 2a<\/span><button[^>]*class="hp-btn inline"/);
  assert.match(html, /<span class="badge tone-jp">Level 3<\/span>/);
  assert.equal((html.match(/class="fv-open"/g) ?? []).length, 2, "the third candidate has no structure: no View button");
  assert.match(html, /<td>—<\/td><\/tr>/);
});

test("each row has a toggle button: aria-expanded, named after its candidate", () => {
  const html = view();
  const buttons = [...html.matchAll(/<button type="button" class="cand-toggle"[^>]*>/g)].map(m => m[0]);
  assert.equal(buttons.length, 3);
  assert.ok(buttons.every(b => b.includes('aria-expanded="false"') && !b.includes("aria-controls")));
  assert.match(buttons[0], /aria-label="Show details of Caffeine"/);
  assert.match(buttons[2], /aria-label="Show details of 1 M\+H"/);
});

test("an open row shows the folded columns in a detail line the toggle points at", () => {
  const html = view({ defaultOpen: ["c2"] });
  const toggles = [...html.matchAll(/<button type="button" class="cand-toggle"[^>]*>/g)].map(m => m[0]);
  assert.ok(toggles[0].includes('aria-expanded="false"') && toggles[1].includes('aria-expanded="true"'));
  const target = /aria-controls="([^"]+)"/.exec(toggles[1])[1];
  assert.match(html, new RegExp(`<tr class="cand-detail" id="${target.replace(/[:]/g, "\\$&")}">`));
  assert.match(toggles[1], /aria-label="Hide details of Isocaffeine"/);
  const detail = html.slice(html.indexOf('<tr class="cand-detail"'));
  const facts = Object.fromEntries([...detail.matchAll(/<dt>(.*?)<\/dt><dd>(.*?)<\/dd>/g)].map(m => [m[1], m[2]]));
  assert.deepEqual(facts, {
    SMILES: "Cn1c(=O)c2ncn(C)c2n(C)c1=O", Formula: "C8H10N4O2", Adduct: "[M+H]+", Layer: "library",
    "Detectability (omics / ionizable / adduct)": "metabolome / moderate / high", "Isomer role": "not separated · challenger",
  });
  assert.equal((html.match(/class="cand-detail"/g) ?? []).length, 1, "only the open row has a detail line");
});

test("missing values read as a dash in the detail line", () => {
  assert.deepEqual(candidateDetails(candidates[2]).map(([, v]) => v), ["—", "—", "[M+H]+", "library", "—", "—"]);
});

test("no candidates: one empty line", () => {
  const html = view({ candidates: [] });
  assert.match(html, /<td colSpan="6" class="empty">no candidates<\/td>|<td colSpan="6" class="empty">no candidates/);
});

// ---------------------------------------------------------------- CCS check (supporting evidence)
const ccs = (over = {}) => ({ verdict: "Consistent", source: "measured", adduct: "[M+H]+", observed: 139.6, reference: 139.6, delta_pct: 0.64, tolerance_pct: 1, ref_label: "DTCCS N2",
  dataset: "Baker2017", n_refs: 6, in_domain: null, note: "", flag: null, moved_from: null, ...over });
const withCcs = [
  candidate("c1", "Isocaffeine", { ccs: ccs({ source: "predicted", reference: 139.3, delta_pct: 0.2, tolerance_pct: 6.5, ref_label: "SVR", dataset: null }) }),
  candidate("c2", "Caffeine", { rank: 2, ccs: ccs({ verdict: "Inconsistent", delta_pct: 9.8, flag: "CCS ≠ measured", moved_from: 1 }) }),
  candidate("c3", "Aspirin", { rank: 3, ccs: ccs({ verdict: "Inconsistent", source: "predicted", reference: 150, delta_pct: -6.33, tolerance_pct: 5.6, ref_label: "SVR", dataset: null, flag: "CCS ≠ predicted" }) }),
  candidate("c4", "Unknown", { rank: 4, smiles: null, ccs: { verdict: "No structure", note: "a CCS can only be looked up for a structure" } }),
];

test("the key columns do not change with a CCS check; a candidate the rule applied to gets a short flag after its name", () => {
  const html = view({ candidates: withCcs });
  assert.deepEqual(heads(html), ["Details", "Rank", "Name", "Level", "Final score", "Fragments"]);
  assert.match(html, /<td class="cand-name">Caffeine<span class="cand-flag" title="Inconsistent · \+9\.80 % \(±1\.0 %\) · measured 139\.6 Å² · DTCCS N2 · Baker2017 · moved from rank 1">CCS ≠ measured<\/span><\/td>/);
  assert.match(html, /<td class="cand-name">Aspirin<span class="cand-flag"[^>]*>CCS ≠ predicted<\/span><\/td>/);
  assert.equal((html.match(/cand-flag/g) ?? []).length, 2, "Consistent and No structure rows carry no flag");
  assert.match(html, /<td class="cand-name">Isocaffeine<\/td>/);
  assert.ok(!html.replace(/ title="[^"]*"/g, "").includes("Baker2017"), "the numbers stay in the detail line (the flag's hover title aside)");
});

test("the detail line shows the CCS: a verdict badge, then delta, tolerance, source and where it moved from", () => {
  const html = view({ candidates: withCcs, defaultOpen: ["c1", "c2", "c3", "c4"] });
  const ccsFacts = [...html.matchAll(/<div><dt>CCS<\/dt><dd>(.*?)<\/dd><\/div>/g)].map(m => m[1]);
  assert.equal(ccsFacts.length, 4);
  assert.equal(ccsFacts[0], '<span class="badge tone-good">Consistent</span><span class="cand-ccs"> +0.20 % (±6.5 %) · predicted 139.3 Å² · SVR</span>');
  assert.equal(ccsFacts[1], '<span class="badge tone-warn">Inconsistent</span><span class="cand-ccs"> +9.80 % (±1.0 %) · measured 139.6 Å² · DTCCS N2 · Baker2017</span><span class="cand-ccs"> · moved from rank 1</span>');
  assert.match(ccsFacts[2], /^<span class="badge tone-warn">Inconsistent<\/span><span class="cand-ccs"> −6\.33 % \(±5\.6 %\) · predicted 150\.0 Å² · SVR<\/span>$/);
  assert.equal(ccsFacts[3], '<span class="badge tone-neutral">No structure</span><span class="cand-ccs"> a CCS can only be looked up for a structure</span>');
});

test("without a CCS check the detail line and the closed table are exactly as before", () => {
  const html = view({ defaultOpen: ["c1"] });
  assert.doesNotMatch(html, /CCS|cand-flag|cand-ccs/);
  assert.equal(candidateDetails(candidates[0]).length, 6);
});

test("ccsNumbers and ccsText: the numbers when there are some, the note when there are not", () => {
  assert.equal(ccsNumbers(ccs({ delta_pct: -3.37, tolerance_pct: 2, source: "measured", reference: 145.4, ref_label: "TWCCS N2", dataset: "moll0218" })), "−3.37 % (±2.0 %) · measured 145.4 Å² · TWCCS N2 · moll0218");
  assert.equal(ccsNumbers({ verdict: "Not checked", note: "Unknown adduct: '[M+Xx]+'." }), "Unknown adduct: '[M+Xx]+'.");
  assert.equal(ccsNumbers({ verdict: "No structure" }), "");
  assert.equal(ccsText({ verdict: "No structure" }), "No structure");
  assert.equal(ccsText(ccs({ moved_from: 3 })), "Consistent · +0.64 % (±1.0 %) · measured 139.6 Å² · DTCCS N2 · Baker2017 · moved from rank 3");
});
