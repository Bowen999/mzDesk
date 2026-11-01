import assert from "node:assert/strict";
import { after, test } from "node:test";
import { startSsr } from "./ssr-helper.mjs";

const ssr = await startSsr("answer");
const { AnswerCard } = await ssr.load("/src/components/AnswerCard.tsx");
const { answerSummary } = await ssr.load("/src/components/annotationModel.ts");
after(() => ssr.close());

const STATEMENT = "The spectrum cannot separate 2 isomers: Caffeine and Isocaffeine. The top-ranked one is Caffeine (library, Level 2a). "
  + "The isomer verdict only holds for the candidates that were considered; an isomer that is not in the candidate list has not been ruled out.";
const candidate = (id, name, extra = {}) => ({
  id, rank: 1, name, smiles: "CN1C=NC2=C1C(=O)N(C)C(=O)N2C", formula: "C8H10N4O2", adduct: "[M+H]+", layer: "library", layer_label: "library",
  level: "2a", level_label: "Level 2a: probable structure (library spectrum match)", level_reason: "", final_score: 0.8, base_score: 1, hypothesis: false,
  evidence: { library: { gate_passed: true } }, detect: null, isomer: null, ...extra,
});
const group = (over = {}) => ({
  answer: { type: "group", ids: ["c1", "c2"], layer: "library", layer_label: "library", level: "2a", level_label: "Level 2a: probable structure (library spectrum match)",
    statement: STATEMENT, top_id: "c1", isomer_verdict: "indistinguishable" },
  candidates: [candidate("c1", "Caffeine"), candidate("c2", "Isocaffeine", { rank: 2, final_score: 0.835 })],
  warnings: [], steps: [], context: { precursor_mz: 195.0877, polarity: "Positive", adduct: "[M+H]+", retained_adducts: [], omics: "", libraries: [], gate_passed: true },
  isomer: { verdict: "indistinguishable", top_id: "c1", group: ["c1", "c2"], challengers: [], notes: [] }, ...over,
});
const view = (result, requestedAdduct = null) => ssr.render(AnswerCard, { result, requestedAdduct });
const text = html => html.replace(/<[^>]*>/g, " ").replace(/\s+/g, " ");

test("the headline is the compound, its level and its score", () => {
  const html = view(group());
  assert.match(html, /<p class="answer-headline"><span class="answer-name"><button type="button" class="cell-depict" aria-expanded="false">Caffeine<\/button><\/span>/);
  assert.match(html, /<span class="badge tone-good" title="Level 2a: probable structure[^"]*">Level 2a<\/span>/);
  assert.match(html, /title="Final score 0.800"><span class="sr-only">Final score <\/span>0\.80<\/span>/);
  assert.equal((html.match(/answer-headline/g) ?? []).length, 1);
});

test("a group gets one short secondary line and no banner repeating the isomer message", () => {
  const html = view(group());
  assert.match(html, /<p class="answer-sub warn">Indistinguishable from Isocaffeine by this spectrum<\/p>/);
  assert.equal((html.match(/answer-sub/g) ?? []).length, 1);
  assert.doesNotMatch(html, /class="alert/);
  assert.doesNotMatch(html, /cannot be separated by the spectrum|The ranking only orders them/);
  assert.doesNotMatch(html, /class="pill-row"|Library gate passed|Isomer check: /, "the chips row is gone");
});

test("the full statement and the other facts stay behind Details", () => {
  const html = view(group());
  assert.match(html, /<details class="answer-details"><summary>Details<\/summary>/);
  assert.ok(html.includes(STATEMENT), "the original statement is kept word for word");
  const details = text(html.slice(html.indexOf("<details")));
  for (const fact of ["Layer library", "Level Level 2a: probable structure (library spectrum match)", "Final score 0.800", "Isomer check indistinguishable", "Library gate passed"]) {
    assert.ok(details.includes(fact), `${fact} is in Details`);
  }
  assert.doesNotMatch(html.slice(0, html.indexOf("<details")), /ruled out|top-ranked one is/, "the long text is not outside Details");
});

test("a single answer says where it comes from, and that the isomers were told apart", () => {
  const result = group({ answer: { ...group().answer, type: "single", ids: ["c1"], isomer_verdict: "single", statement: "Caffeine is the top-ranked structure (library, Level 2a)." },
    isomer: { verdict: "single", group: [], challengers: [], notes: [] } });
  const html = view(result);
  assert.match(html, /<p class="answer-sub">Library match, gate passed · Told apart from the other isomers considered<\/p>/);
  assert.equal(answerSummary(result).warn, false);
});

test("conflicting evidence and challengers are on the secondary line, with their m/z in Details", () => {
  const isomer = { verdict: "single", group: [], notes: [], challengers: [{ id: "c2", name: "Isocaffeine", peaks: [{ mz: 110.0713 }, { mz: 138.0662 }] }] };
  const result = group({ answer: { ...group().answer, type: "single", ids: ["c1"], isomer_verdict: "single", statement: "Caffeine is the top-ranked structure." }, isomer });
  const html = view(result);
  assert.match(html, /<p class="answer-sub warn">Library match, gate passed · Conflicting fragment evidence · Isocaffeine explains peaks it cannot<\/p>/);
  assert.ok(text(html).includes("Isocaffeine explains peaks the top-ranked structure cannot (m/z 110.0713, 138.0662)."));
  assert.doesNotMatch(html, /class="alert/, "no separate banner for the conflict or the challenger");
});

test("an isomer check that was not run or could not be done is named, in a few words", () => {
  const single = (verdict, level = "2a") => group({ answer: { ...group().answer, type: "single", ids: ["c1"], isomer_verdict: verdict, level, level_label: level === "2a" ? group().answer.level_label : "Level 3: tentative candidate", layer: level === "2a" ? "library" : "msfinder" },
    candidates: [candidate("c1", "Caffeine", level === "2a" ? {} : { layer: "msfinder", level: "3", level_label: "Level 3: tentative candidate" })], isomer: null });
  assert.equal(answerSummary(single("not run")).secondary, "Library match, gate passed · Isomer check not run");
  assert.equal(answerSummary(single("insufficient")).secondary, "Library match, gate passed · Isomer check inconclusive");
  assert.equal(answerSummary(single("not run", "3")).secondary, "MS-FINDER candidate · Isomer check not run");
});

test("formula-only, mass-only and empty answers still read as one headline", () => {
  const base = group();
  const formula = { ...base, answer: { type: "formula_only", ids: ["f1"], layer: "formula_only", layer_label: "formula-only", level: "4", level_label: "Level 4: molecular formula only",
    statement: "No structure could be proposed. The data support 1 molecular formula: C8H10N4O2 ([M+H]+). This is formula-level evidence only.", top_id: null, isomer_verdict: null },
    candidates: [candidate("f1", "", { layer: "formula_only", layer_label: "formula-only", level: "4", level_label: "Level 4: molecular formula only", smiles: null, final_score: null })], isomer: null };
  const f = view(formula);
  assert.match(f, /<span class="answer-name">C8H10N4O2 \[M\+H\]\+<\/span>/);
  assert.match(f, /Level 4/);
  assert.match(f, /<p class="answer-sub">Formula only, no structure proposed<\/p>/);
  assert.doesNotMatch(f, /answer-score/);
  const mass = view({ ...base, candidates: [], answer: { type: "mass_only", ids: [], layer: "", layer_label: "", level: "5", level_label: "Level 5: exact mass only", statement: "Only the exact mass is available.", top_id: null, isomer_verdict: null }, isomer: null });
  assert.match(mass, /<span class="answer-name">m\/z 195\.0877<\/span>/);
  assert.match(mass, /Exact mass only/);
  const none = view({ ...base, candidates: [], answer: { type: "none", ids: [], statement: "Workflow result retained", top_id: null, isomer_verdict: "not run" }, isomer: null, context: { precursor_mz: 1, polarity: "Positive", gate_passed: false } });
  assert.match(none, /<p class="verdict-text">Workflow result retained<\/p>/);
  assert.doesNotMatch(none, /answer-headline|<details/);
});

test("adduct notes and the run's warnings are still alerts, and only those", () => {
  const html = view(group({ warnings: ["Peaks above the precursor were ignored."] }), "[M+Na]+");
  assert.equal((html.match(/class="alert warn"/g) ?? []).length, 2);
  assert.match(html, /Input adduct: \[M\+Na\]\+\. A library hit passed with \[M\+H\]\+; review the library checks\./);
  assert.match(html, /<div class="alert warn">Peaks above the precursor were ignored\.<\/div>/);
  // no library hit passed with another adduct, but the top candidate uses one
  const noGate = group({ candidates: [candidate("c1", "Caffeine", { evidence: {} }), candidate("c2", "Isocaffeine", { rank: 2, evidence: {} })] });
  const html2 = view(noGate, "[M+Na]+");
  assert.equal((html2.match(/class="alert warn"/g) ?? []).length, 1);
  assert.match(html2, /Input adduct: \[M\+Na\]\+\. The top candidate uses \[M\+H\]\+; review the adduct evidence\./);
  assert.doesNotMatch(view(noGate, null), /class="alert/);
});

test("a result with the fields missing does not break the card", () => {
  const sparse = { answer: { type: "none", ids: [], statement: "x", top_id: null, isomer_verdict: null }, candidates: [], context: { gate_passed: false } };
  assert.doesNotThrow(() => view(sparse));
  assert.doesNotThrow(() => view({ ...group(), candidates: [], isomer: null }));
});

// ---------------------------------------------------------------- CCS check (supporting evidence)
const ccs = (over = {}) => ({ verdict: "Consistent", source: "measured", adduct: "[M+H]+", observed: 139.6, reference: 139.6, delta_pct: 0, tolerance_pct: 1, ref_label: "DTCCS N2",
  dataset: "Baker2017", n_refs: 6, in_domain: null, note: "", flag: null, moved_from: null, ...over });
const withCcs = x => group({ candidates: [candidate("c1", "Caffeine", { ccs: x }), candidate("c2", "Isocaffeine", { rank: 2, final_score: 0.835 })] });
const detailsOf = html => text(html.slice(html.indexOf("<details")));

test("Details lists the top candidate's CCS line, and the headline does not change", () => {
  const html = view(withCcs(ccs()));
  assert.ok(detailsOf(html).includes("CCS Consistent · +0.00 % (±1.0 %) · measured 139.6 Å² · DTCCS N2 · Baker2017"));
  assert.equal(html.slice(0, html.indexOf("<details")), view(group()).slice(0, view(group()).indexOf("<details")), "the compact card above Details is the one without CCS");
  assert.deepEqual(answerSummary(withCcs(ccs())).facts.at(-1), ["CCS", "Consistent · +0.00 % (±1.0 %) · measured 139.6 Å² · DTCCS N2 · Baker2017"]);
});

test("a demoted top, a predicted value and a candidate without a reference read plainly", () => {
  const moved = ccs({ verdict: "Inconsistent", delta_pct: 9.8, reference: 140.5, flag: "CCS ≠ measured", moved_from: 1 });
  assert.ok(detailsOf(view(withCcs(moved))).includes("CCS Inconsistent · +9.80 % (±1.0 %) · measured 140.5 Å² · DTCCS N2 · Baker2017 · moved from rank 1"));
  const predicted = ccs({ verdict: "Inconsistent", source: "predicted", reference: 150, delta_pct: -6.33, tolerance_pct: 5.6, ref_label: "SVR", dataset: null, in_domain: false });
  assert.ok(detailsOf(view(withCcs(predicted))).includes("CCS Inconsistent · −6.33 % (±5.6 %) · predicted 150.0 Å² · SVR · outside its training domain"));
  const none = ccs({ verdict: "No reference", source: null, reference: null, delta_pct: null, tolerance_pct: null, ref_label: null, dataset: null, note: "No prediction model is installed." });
  assert.ok(detailsOf(view(withCcs(none))).includes("CCS No reference · No prediction model is installed."));
});

test("no CCS check, no CCS fact", () => {
  assert.ok(!answerSummary(group()).facts.some(([k]) => k === "CCS"));
  assert.doesNotMatch(view(group()), /CCS/);
});
