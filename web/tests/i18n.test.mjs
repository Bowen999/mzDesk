import assert from "node:assert/strict";
import { execFileSync } from "node:child_process";
import { readFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { after, afterEach, test } from "node:test";
import { fileURLToPath } from "node:url";
import { React, startSsr } from "./ssr-helper.mjs";

// The language switch: t() with its English fallback, the dictionary store, and localizeSpec (which display text of the page list is
// translated, and what must never be: ids, option lists, values, lookup keys).
const ssr = await startSsr("i18n");
const { LANGS, localizeRunBlocks, localizeSpec, pickLanguage, setDictionary, setLang, specContext, t, useLang, useT } = await ssr.load("/src/i18n/index.ts");
const { OutputBlock } = await ssr.load("/src/components/Outputs.tsx");
const { StepView } = await ssr.load("/src/components/StepList.tsx");
const { FieldInput, summarizeInputs } = await ssr.load("/src/components/Fields.tsx");
const { RichText, emphasize } = await ssr.load("/src/components/Help.tsx");
const { buildHtml } = await ssr.load("/src/report/build.ts");
const { networkUses } = await ssr.load("/src/network.ts");
after(() => ssr.close());
afterEach(() => setLang("en"));

const here = dirname(fileURLToPath(import.meta.url));
const repo = dirname(dirname(here));

test("the six languages, in the order of the menu", () => {
  assert.deepEqual(LANGS.map(l => [l.id, l.short, l.name]), [
    ["en", "EN", "English"], ["zh", "中", "中文"], ["ja", "日", "日本語"], ["de", "DE", "Deutsch"], ["fr", "FR", "Français"], ["es", "ES", "Español"],
  ]);
});

test("the browser's language is used when it is one of the six; zh-* is zh; the first usable one wins", () => {
  assert.equal(pickLanguage(["de-AT", "en"]), "de");
  assert.equal(pickLanguage(["zh-TW"]), "zh");
  assert.equal(pickLanguage(["ZH-hans-CN"]), "zh");
  assert.equal(pickLanguage(["pt-BR", "ja", "fr"]), "ja");
  assert.equal(pickLanguage(["pt-BR", "ko"]), "en");
  assert.equal(pickLanguage([]), "en");
});

test("no translation: the English text comes back, with its {holes} filled", () => {
  assert.equal(t("Run"), "Run");
  assert.equal(t("Download the results of {n} pages", { n: 3 }), "Download the results of 3 pages");
  assert.equal(t("{a} and {a}, not {b}", { a: "x" }), "x and x, not {b}", "a hole without a value stays as it is");
  assert.equal(t("Save {n}"), "Save {n}", "no values given: the text is left alone");
  assert.equal(t("constructor"), "constructor", "inherited object keys are not translations");
});

test("a dictionary translates; a missing or empty entry falls back to English; the choice is shown on <html lang>", () => {
  const lang = [];
  globalThis.document = { documentElement: { setAttribute: (k, v) => lang.push([k, v]), removeAttribute() {} } };
  setDictionary("de", { Run: "Ausführen", "Download the results of {n} pages": "Ergebnisse von {n} Seiten laden", Clear: "" });
  setLang("de");
  assert.equal(t("Run"), "Ausführen");
  assert.equal(t("Download the results of {n} pages", { n: 2 }), "Ergebnisse von 2 Seiten laden");
  assert.equal(t("Clear"), "Clear", "an empty translation counts as missing");
  assert.equal(t("Not in the dictionary"), "Not in the dictionary");
  setLang("zh");
  assert.equal(t("Run"), "Run", "a language whose dictionary has not arrived shows English");
  setLang("en");
  assert.equal(t("Run"), "Run");
  assert.deepEqual(lang.map(l => l[1]), ["de", "zh-CN", "en"], "zh is written zh-CN");
});

test("useT and useLang follow the language", () => {
  setDictionary("fr", { Run: "Lancer" });
  const Probe = () => React.createElement("p", null, `${useLang().lang}:${useT()("Run")}`);
  assert.equal(ssr.render(Probe), "<p>en:Run</p>");
  setLang("fr");
  assert.equal(ssr.render(Probe), "<p>fr:Lancer</p>");
});

test("useT returns a new function when the dictionary changes, so memos that list it re-run", () => {
  const seen = [];
  const Probe = () => { seen.push(useT()); return null; };
  ssr.render(Probe);
  setDictionary("es", { Run: "Ejecutar" });
  setLang("es");
  ssr.render(Probe);
  assert.notEqual(seen[0], seen[1]);
});

test("a context picks the entry 'text|context', then the plain text, then English (without the context)", () => {
  setDictionary("de", { System: "System-Plain", "System|theme": "Design wie System", "Show {n}": "Zeige {n}", "Show {n}|short": "{n} zeigen" });
  setLang("de");
  assert.equal(t("System", { context: "theme" }), "Design wie System");
  assert.equal(t("System", { context: "other" }), "System-Plain", "no entry for that context: the plain translation");
  assert.equal(t("System"), "System-Plain");
  assert.equal(t("Unlisted", { context: "theme" }), "Unlisted", "no entry at all: the English text, not 'Unlisted|theme'");
  assert.equal(t("Show {n}", { context: "short", n: 2 }), "2 zeigen", "{holes} are filled as usual");
  assert.equal(t("Show {n}", { n: 2 }), "Zeige 2");
  setLang("en");
  assert.equal(t("System", { context: "theme" }), "System", "English shows the text alone");
  assert.equal(t("Show {n}", { context: "short", n: 3 }), "Show 3");
});

test("the words the translators flagged have their own entries: System in the theme switch, gives in the one-line flow", () => {
  const dict = lang => JSON.parse(readFileSync(join(here, `../src/i18n/locales/${lang}.json`), "utf8"));
  for (const lang of ["zh", "ja", "de", "fr", "es"]) {
    const words = dict(lang);
    for (const key of ["System|theme", "gives|flow", "Column|metadata", "Species|lipid"]) assert.ok(words[key], `${lang}: ${key}`);
  }
  assert.notEqual(dict("zh")["Column"], dict("zh")["Column|metadata"], "chromatography column / table column");
  assert.notEqual(dict("zh")["Species"], dict("zh")["Species|lipid"], "organism / lipid species");
  assert.notEqual(dict("ja")["gives"], dict("ja")["gives|flow"], "reaction product / the arrow of 'have -> get'");
});

// ---------------------------------------------------------------- localizeSpec
const spec = () => ({
  version: "1.2.3", hosted: true,
  groups: [{ id: "modules", title: "Modules" }],
  subgroups: [{ id: "ms2", title: "Spectrum" }],
  pipelines: [{
    id: "wf", title: "Workflow A", one_liner: "Spectrum → answer", based_on: ["MetFrag", "mzDesk"],
    inputs: [{ type: "select", id: "mode", label: "Mode", options: ["Quick", "Full"], value: "Quick", section: "Advanced", help: ["First tip", "Second tip"] }],
    stages: [{ id: 1, title: "Adduct", modules: [], expected_s: 1 }],
    examples: [{ label: "Caffeine example", values: { mode: "Quick", name: "Caffeine example" } }],
  }],
  modules: [{
    id: "mod", title: "Module A", cls: "A", group: "modules", subgroup: "ms2", live: true, one_liner: "Peaks → clean peaks", network: "Sent to GNPS.",
    based_on: ["RDKit"], example: { omics: "Lipidomics", label: "Compound" },
    inputs: [
      { type: "precursor", id: "precursor", label: "Precursor m/z", adduct: true, auto_label: "Auto (any adduct)", ion_modes: ["Positive"],
        adducts: { Positive: ["Auto", "[M+H]+"], Negative: ["Auto"] }, help: "Required.", section: "Query" },
      { type: "text", id: "t", label: "Compound", placeholder: "e.g. caffeine", value: "Compound", help: "Name or ID" },
      { type: "select", id: "omics", label: "Omics", options: ["Either / unknown", "Lipidomics"], value: "Either / unknown", placeholder: "Pick one" },
      { type: "table", id: "alk", label: "Alkane series", columns: ["Carbon number", "RT"], caption: "One row per line" },
      { type: "heading", label: "Optional part" },
    ],
    outputs: [
      { type: "table", id: "out", title: "Results", columns: ["Name", "Fragments"], caption: "Sorted by score", key: "k", help: ["Row tip"],
        column_help: { Fragments: "View: explained peaks." } },
      { type: "cards", id: "cards", title: "Summary", items: [{ label: "Best hit", unit: "ppm", help: "Item tip" }, { label: "Peaks", unit: "" }] },
      { type: "figure", id: "fig", label: "Spectrum plot", span: "half" },
      { type: "button", label: "Download", key: "b" },
    ],
  }],
});

const mark = s => `«${s}»`;
const everyString = () => {
  const found = new Set();
  const visit = v => { if (typeof v === "string") found.add(v); else if (Array.isArray(v)) v.forEach(visit); else if (v && typeof v === "object") Object.entries(v).forEach(([k, x]) => { found.add(k); visit(x); }); };
  visit(spec());
  return Object.fromEntries([...found].filter(Boolean).map(s => [s, mark(s)]));
};
const deepFreeze = o => { if (o && typeof o === "object") { Object.values(o).forEach(deepFreeze); Object.freeze(o); } return o; };

test("localizeSpec translates the display text, never the ids, values or option lists, and does not touch its input", () => {
  setDictionary("ja", everyString());
  setLang("ja");
  const input = deepFreeze(spec());
  const out = localizeSpec(input);
  assert.deepEqual(input, spec(), "the input is unchanged");
  assert.notEqual(out, input);

  const m = out.modules[0], p = out.pipelines[0];
  assert.deepEqual([out.groups[0].title, out.subgroups[0].title], [mark("Modules"), mark("Spectrum")]);
  assert.deepEqual([m.title, m.one_liner, m.network], [mark("Module A"), mark("Peaks → clean peaks"), mark("Sent to GNPS.")]);
  assert.deepEqual([p.title, p.one_liner, p.examples[0].label], [mark("Workflow A"), mark("Spectrum → answer"), mark("Caffeine example")]);
  assert.deepEqual(p.inputs[0].help, [mark("First tip"), mark("Second tip")], "a help list is translated item by item");
  const [prec, text, omics, table, heading] = m.inputs;
  assert.deepEqual([prec.label, prec.auto_label, prec.help], [mark("Precursor m/z"), mark("Auto (any adduct)"), mark("Required.")]);
  assert.deepEqual([text.label, text.placeholder, text.help], [mark("Compound"), mark("e.g. caffeine"), mark("Name or ID")]);
  assert.equal(omics.placeholder, mark("Pick one"));
  assert.deepEqual([table.label, table.columns, table.caption], [mark("Alkane series"), [mark("Carbon number"), mark("RT")], mark("One row per line")]);
  assert.equal(heading.label, mark("Optional part"));
  const [results, cards, figure, button] = m.outputs;
  assert.deepEqual([results.title, results.columns, results.caption, results.help], [mark("Results"), [mark("Name"), mark("Fragments")], mark("Sorted by score"), [mark("Row tip")]]);
  assert.deepEqual(results.column_help, { Fragments: mark("View: explained peaks."), [mark("Fragments")]: mark("View: explained peaks.") },
    "column help is found under the translated column name and, for a table whose columns come from a run, under the English one");
  assert.equal(cards.title, mark("Summary"));
  assert.deepEqual(cards.items, [{ label: "Best hit", unit: mark("ppm"), help: mark("Item tip") }, { label: "Peaks", unit: "" }],
    "a cards label is the key into the run's values: it stays English; unit and help are translated");
  assert.equal(figure.label, mark("Spectrum plot"));
  assert.equal(button.label, mark("Download"));

  // what must not change: ids, types, flags, numbers, option lists, values, adduct lists, example inputs, section titles, keys, engine names
  assert.deepEqual([m.id, m.cls, m.group, m.subgroup, m.live, m.based_on, m.example], ["mod", "A", "modules", "ms2", true, ["RDKit"], { omics: "Lipidomics", label: "Compound" }]);
  assert.deepEqual([prec.id, prec.type, prec.adduct, prec.ion_modes, prec.adducts, prec.section], ["precursor", "precursor", true, ["Positive"], { Positive: ["Auto", "[M+H]+"], Negative: ["Auto"] }, "Query"]);
  assert.deepEqual([text.value, omics.options, omics.value, omics.id], ["Compound", ["Either / unknown", "Lipidomics"], "Either / unknown", "omics"]);
  assert.deepEqual([results.id, results.key, button.key, figure.span], ["out", "k", "b", "half"]);
  assert.deepEqual(p.inputs[0].options, ["Quick", "Full"]);
  assert.equal(p.inputs[0].section, "Advanced");
  assert.deepEqual(p.examples[0].values, { mode: "Quick", name: "Caffeine example" });
  assert.deepEqual([p.stages, p.based_on], [[{ id: 1, title: "Adduct", modules: [], expected_s: 1 }], ["MetFrag", "mzDesk"]]);
  assert.deepEqual([out.version, out.hosted], ["1.2.3", true]);
});

test("localizeSpec leaves a string alone when the dictionary has no entry for it", () => {
  setDictionary("ja", { "Module A": "モジュールA" });
  setLang("ja");
  const out = localizeSpec(spec());
  assert.equal(out.modules[0].title, "モジュールA");
  assert.equal(out.modules[0].one_liner, "Peaks → clean peaks");
  assert.deepEqual(out.modules[0].outputs[0].columns, ["Name", "Fragments"]);
});

test("localizeSpec gives a string its context in the blocks spec-rules.json names, and a cards label gets it where it is drawn", () => {
  setDictionary("de", { Column: "Säule", "Column|metadata": "Spalte", Value: "Wert", Species: "Spezies", "Species|lipid": "Lipidspezies", Chains: "Ketten", Results: "Ergebnisse" });
  setLang("de");
  const table = (id, columns, extra = {}) => ({ type: "table", id, title: "Results", columns, ...extra });
  const cards = id => ({ type: "cards", id, title: "Results", items: [{ label: "Species", unit: "Chains", help: "Chains" }] });
  const out = localizeSpec({
    groups: [], subgroups: [], pipelines: [],
    modules: [{ id: "m", title: "Results", outputs: [table("redu", ["Column", "Value"], { column_help: { Column: "Chains" } }), table("exact", ["Column"]), cards("lipid_maps"), cards("other")] }],
  });
  const [redu, exact, lipid, other] = out.modules[0].outputs;
  assert.deepEqual(redu.columns, ["Spalte", "Wert"], "a table column in a block with a context");
  assert.deepEqual(redu.column_help, { Column: "Ketten", Spalte: "Ketten" });
  assert.deepEqual(exact.columns, ["Säule"], "the same word in another block: the plain translation");
  assert.equal(lipid.items[0].label, "Species", "a cards label stays English: it is the key into the run's values");
  assert.equal(lipid.items[0].unit, "Ketten");

  assert.deepEqual(specContext("lipid_maps", "Species"), { context: "lipid" });
  assert.equal(specContext("lipid_maps", "Chains"), undefined);
  assert.equal(specContext("other", "Species"), undefined);
  assert.equal(specContext("constructor", "Species"), undefined, "inherited object keys are not blocks");
  assert.equal(specContext("redu", "constructor"), undefined);
  assert.equal(specContext(undefined, "Species"), undefined);

  const label = block => ssr.render(OutputBlock, { block, status: "none" }).match(/class="card-label">([^<]*)</)[1];
  assert.equal(label(cards("lipid_maps")), "Lipidspezies");
  assert.equal(label(cards("other")), "Spezies");
  setLang("en");
  assert.equal(label(cards("lipid_maps")), "Species", "English: the label as written");
});

test("the blocks of a run get their titles translated; notes, cells and everything else stay as the run wrote them", () => {
  setDictionary("de", { "Library hits": "Bibliothekstreffer", Notes: "Hinweise", "Query spectrum": "Abfragespektrum", Name: "Name" });
  const blocks = () => [
    { type: "table", id: "lib", title: "Library hits", columns: ["Name"], data: { rows: [["Library hits"]] } },
    { type: "notes", id: "n", title: "Notes", data: { items: ["Notes"] } },
    { type: "figure", id: "f", title: null, label: "Query spectrum", data: {} },
    { type: "figure", id: "g", title: null, label: "Query vs caffeine (MassBank)", data: {} },
    { type: "table", id: "x", title: "Not in the dictionary", data: {} },
    { type: "button", id: "b", title: null, label: "Library hits", data: {} },
  ];
  const input = deepFreeze(blocks());
  assert.equal(localizeRunBlocks(input), input, "English: the same array");
  setLang("de");
  const out = localizeRunBlocks(input);
  assert.deepEqual(out.map(b => b.title ?? b.label), ["Bibliothekstreffer", "Hinweise", "Abfragespektrum", "Query vs caffeine (MassBank)", "Not in the dictionary", "Library hits"]);
  assert.equal(out[2].title, null);
  assert.equal(out[0].data, input[0].data, "the data is shared, not copied");
  assert.deepEqual(out[0].data.rows, [["Library hits"]], "a cell stays");
  assert.deepEqual(out[1].data.items, ["Notes"], "a note stays");
  assert.equal(out[4], input[4], "an untouched block is the same object");

  const step = { id: 3, title: "Library search", status: "done", seconds: 0.1, note: "Library hits found", blocks: blocks().slice(0, 2) };
  const html = ssr.render(StepView, { step, index: 2 });
  assert.match(html, /<h4 class="out-title">Bibliothekstreffer/);
  assert.match(html, /<h4 class="out-title">Hinweise/);
  assert.match(html, /Library search/, "the step title comes from the run");
  assert.match(html, /Library hits found/, "the step note comes from the run");
  assert.match(html, />Library hits</, "a cell of the table");
});

test("English: the page list comes back as it is", () => {
  const s = spec();
  assert.equal(localizeSpec(s), s);
});

// ---------------------------------------------------------------- the real page list (needs Python; skipped when it cannot run)
function realSpec() {
  const code = "import json; from mzdesk import spec_i18n as s; from mzdesk.spec import spec; sp = spec(); "
    + "print(json.dumps({'spec': sp, 'shown': s.displayed_strings(sp), 'lookup': s.lookup_strings(sp)}))";
  for (const python of [process.env.MZDESK_PYTHON, "/opt/anaconda3/bin/python3", "python3"].filter(Boolean)) {
    try { return JSON.parse(execFileSync(python, ["-c", code], { cwd: repo, encoding: "utf8", maxBuffer: 64e6, stdio: ["ignore", "pipe", "ignore"] })); } catch { /* try the next */ }
  }
  return null;
}

test("localizeSpec and mzdesk/spec_i18n.py agree on the real page list", { timeout: 120_000 }, ctx => {
  const real = realSpec();
  if (!real) { ctx.skip("python with mzdesk is not available"); return; }
  const listed = JSON.parse(readFileSync(join(here, "../src/i18n/source/spec.json"), "utf8"));
  const original = JSON.stringify(real.spec);
  const every = Object.fromEntries(listed.map(s => [s, mark(s)]));
  setDictionary("es", every);
  setLang("es");
  const out = localizeSpec(real.spec);
  assert.equal(JSON.stringify(real.spec), original, "the input is unchanged");

  // a string with a context is looked up as "text|context": the dictionary key is what the mark holds
  const translated = new Set(), bad = [];
  const diff = (a, b, path) => {
    if (typeof a === "string") {
      if (a !== b) { const key = b.slice(1, -1); translated.add(key); if (b !== mark(key) || (key !== a && !key.startsWith(a + "|"))) bad.push(`${path}: ${a} -> ${b}`); if (/\.(options|value|values|example|adducts|ion_modes|based_on|section|id)\b/.test(path)) bad.push(`${path} must not be translated`); }
    } else if (Array.isArray(a)) a.forEach((x, i) => diff(x, b[i], `${path}[${i}]`));
    else if (a && typeof a === "object") for (const k of Object.keys(a)) diff(a[k], b[k], `${path}.${k}`);
  };
  diff(real.spec, out, "spec");
  assert.deepEqual(bad, []);
  assert.deepEqual([...translated].sort(), real.shown, "the strings localizeSpec changes are exactly displayed_strings()");
  const lookup = real.lookup.filter(s => !translated.has(s));
  assert.deepEqual(listed.filter(s => !translated.has(s)), lookup, "the rest of spec.json is the lookup keys the components translate when they draw them");
});

// ---------------------------------------------------------------- option labels, lead-ins, "(optional)", the language of the report

/** The element tree FieldInput returns (its hooks run inside this probe), to find a control and call its handler the way the browser does. */
function treeOf(props) {
  let tree;
  ssr.render(() => { tree = FieldInput(props); return null; });
  return tree;
}
const findAll = (node, test, found = []) => {
  if (Array.isArray(node)) node.forEach(n => findAll(n, test, found));
  else if (node && typeof node === "object" && node.props) { if (test(node)) found.push(node); findAll(node.props.children, test, found); }
  return found;
};
const textOf = node => (Array.isArray(node) ? node.map(textOf).join("") : node && typeof node === "object" ? textOf(node.props?.children) : node == null || typeof node === "boolean" ? "" : String(node));

const german = { Positive: "Positiv", Negative: "Negativ", None: "Keine", "Either / unknown": "Beides / unbekannt", "Formic acid": "Ameisensäure", Lipidomics: "Lipidomik", Metabolomics: "Metabolomik", Pick: "Wählen" };

test("a select shows its options translated and still sends the English value", () => {
  setDictionary("de", german);
  setLang("de");
  const sent = [];
  const field = { type: "select", id: "omics", label: "Omics", options: ["Either / unknown", "Lipidomics", ""], placeholder: "Pick", value: "Either / unknown" };
  const options = findAll(treeOf({ field, value: "Either / unknown", onChange: v => sent.push(v) }), n => n.type === "option");
  assert.deepEqual(options.map(o => [o.props.value, textOf(o)]), [["Either / unknown", "Beides / unbekannt"], ["Lipidomics", "Lipidomik"], ["", "Pick"]],
    "the value is the English option; the empty choice shows the placeholder (the page list translates that one)");
  const html = ssr.render(FieldInput, { field, value: "Lipidomics", onChange() {} });
  assert.match(html, /<option value="Lipidomics" selected="">Lipidomik<\/option>/, "the chosen one is found by its English value");
  assert.match(html, /<option value="Either \/ unknown">Beides \/ unbekannt<\/option>/);
  setLang("en");
  assert.match(ssr.render(FieldInput, { field, value: "Lipidomics", onChange() {} }), /<option value="Either \/ unknown">Either \/ unknown<\/option>/, "English: as it was");
});

test("a radio group shows its options translated, compares and sends the English value", () => {
  setDictionary("de", german);
  setLang("de");
  const sent = [];
  const field = { type: "radio", id: "mode", label: "Ion mode", options: ["Positive", "Negative"], value: "Positive" };
  const tree = treeOf({ field, value: "Positive", onChange: v => sent.push(v) });
  const labels = findAll(tree, n => n.type === "label");
  assert.deepEqual(labels.map(l => textOf(l).trim()), ["Positiv", "Negativ"]);
  findAll(labels[1], n => n.type === "input")[0].props.onChange();
  assert.deepEqual(sent, ["Negative"]);
  const html = ssr.render(FieldInput, { field, value: "Negative", onChange() {} });
  assert.match(html, /<input type="radio"[^>]*\/> Positiv<\/label><label><input type="radio"[^>]* checked=""\/> Negativ<\/label>/, "the checked one is the English value");
});

test("a multiselect (chips and tick boxes) shows its options translated and sends the English values, in the options' order", () => {
  setDictionary("de", german);
  setLang("de");
  for (const checkboxes of [false, true]) {
    const sent = [];
    const field = { type: "multiselect", id: "add", label: "Additive", options: ["Formic acid", "None", "Metabolomics"], value: ["None"], ...(checkboxes ? { checkboxes: true } : {}) };
    const tree = treeOf({ field, value: ["None"], onChange: v => sent.push(v) });
    const controls = checkboxes ? findAll(tree, n => n.type === "label" && n.props.className === "check-line") : findAll(tree, n => n.type === "button" && /chip/.test(n.props.className));
    assert.deepEqual(controls.map(c => textOf(c).trim()), ["Ameisensäure", "Keine", "Metabolomik"]);
    const first = controls[0];
    (checkboxes ? findAll(first, n => n.type === "input")[0].props.onChange : first.props.onClick)();
    assert.deepEqual(sent, [checkboxes ? ["Formic acid", "None"] : ["None", "Formic acid"]], checkboxes ? "tick boxes keep the options' order" : "chips append");
    const html = ssr.render(FieldInput, { field, value: ["None"], onChange() {} });
    assert.match(html, checkboxes ? /<input type="checkbox" checked=""\/> Keine/ : /aria-pressed="true">Keine<\/button>/, "the ticked one is found by its English value");
  }
});

test("the input bar names a chosen option in the page's language and leaves what the person typed alone", () => {
  setDictionary("de", { ...german, "{label}: {value}": "{label}: {value}", "{n} selected": "{n} gewählt" });
  setLang("de");
  const fields = [
    { type: "radio", id: "mode", label: "Ion mode", options: ["Positive", "Negative"], value: "Positive" },
    { type: "select", id: "omics", label: "Omics", options: ["", "Lipidomics"], value: "" },
    { type: "multiselect", id: "add", label: "Additive", options: ["Formic acid", "None", "Lipidomics"], value: [] },
    { type: "text", id: "free", label: "Note" },
  ];
  const values = { mode: "Negative", omics: "Lipidomics", add: ["Formic acid", "None"], free: "None" };
  assert.deepEqual(summarizeInputs(fields, values).items, ["Ion mode: Negativ", "Omics: Lipidomik", "Additive: Ameisensäure, Keine", "Note: None"]);
  assert.deepEqual(values, { mode: "Negative", omics: "Lipidomics", add: ["Formic acid", "None"], free: "None" }, "the values themselves are not changed");
  assert.equal(summarizeInputs(fields, { ...values, add: ["Formic acid", "None", "Lipidomics"] }).items[2], "Additive: 3 gewählt");
  setLang("en");
  assert.deepEqual(summarizeInputs(fields, values).items, ["Ion mode: Negative", "Omics: Lipidomics", "Additive: Formic acid, None", "Note: None"]);
});

test("the confirmation before an online run names the chosen options as the page shows them", () => {
  const field = { type: "multiselect", id: "metrics", label: "Metrics", options: ["CosineGreedy", "MS2DeepScore"], network_options: ["MS2DeepScore"] };
  const values = { metrics: ["CosineGreedy", "MS2DeepScore"] };
  assert.deepEqual(networkUses([field], values), ["MS2DeepScore: model from Zenodo."]);
  setDictionary("fr", { MS2DeepScore: "MS2DeepScore (apprentissage)", "{names}: model from Zenodo.": "{names} : modèle depuis Zenodo." });
  setLang("fr");
  assert.deepEqual(networkUses([field], values), ["MS2DeepScore (apprentissage) : modèle depuis Zenodo."]);
});

test("a lead-in is bold before an ASCII colon and a space, before a full-width colon (with or without a space) and before a colon with a space in front", () => {
  // English, as it was
  assert.equal(emphasize("Yes: expected. Partial: only some"), "<b>Yes</b>: expected. <b>Partial</b>: only some");
  assert.equal(emphasize("Basic site (pKa ≥ 9): protonates"), "<b>Basic site (pKa ≥ 9)</b>: protonates");
  assert.equal(emphasize("Ratio 12:30 stays plain and so does a:b"), "Ratio 12:30 stays plain and so does a:b");
  assert.equal(emphasize("Fires. It is slow, so wait: then read."), "Fires. It is slow, so wait: then read.", "no sentence break inside a lead-in; a later one is at most three words");
  assert.equal(emphasize("Uses [M+H]+: the proton adduct"), "<b>Uses <code>[M+H]+</code></b>: the proton adduct", "ions still become tags");
  // Chinese and Japanese: the full-width colon, no space needed; a one-character lead-in; later sentences after 。and ；
  assert.equal(emphasize("模式：文本"), "<b>模式</b>：文本");
  assert.equal(emphasize("模式： 文本"), "<b>模式</b>： 文本");
  assert.equal(emphasize("是：预期的。部分：只有一部分；否：没有"), "<b>是</b>：预期的。<b>部分</b>：只有一部分；<b>否</b>：没有");
  assert.equal(emphasize("塩基性部位（pKa ≥ 9）：プロトン化されます"), "<b>塩基性部位（pKa ≥ 9）</b>：プロトン化されます");
  assert.equal(emphasize("説明です。次を使います。モード：テキスト"), "説明です。次を使います。<b>モード</b>：テキスト", "a sentence break stops the first lead-in; the later one is found");
  // French: a space before the colon stays outside the bold part
  assert.equal(emphasize("Mode : texte"), "<b>Mode</b> : texte");
  assert.equal(emphasize("Oui : attendu. Partiel : seulement une partie. Écart : aucun"), "<b>Oui</b> : attendu. <b>Partiel</b> : seulement une partie. <b>Écart</b> : aucun");
  assert.equal(emphasize("Oui : attendu"), "<b>Oui</b> : attendu", "a no-break space as well");
  // text with its own bold is left alone
  assert.equal(emphasize("<b>Mode</b>：texte. Other: x"), "<b>Mode</b>：texte. Other: x");
  assert.equal(ssr.render(RichText, { text: emphasize("模式：文本") }), "<b>模式</b>：文本", "and it is drawn as bold");
});

test("'(optional)' is left off a label in the input bar: the English suffix and the one of the page's language, with or without a space", () => {
  const items = (fields, values) => summarizeInputs(fields, values).items;
  const english = [{ type: "text", id: "f", label: "Formula (optional)" }, { type: "number", id: "tol", label: "Mass tolerance (ppm)" }, { type: "checkbox", id: "c", label: "Charge states (optional)" }];
  const values = { f: "C6H6", tol: 5, c: true };
  assert.deepEqual(items(english, values), ["Formula: C6H6", "Mass tolerance: 5", "Charge states"], "English: as it was");

  const pages = {
    zh: ["（可选）", "式", "质量容差"], ja: ["（任意）", "式", "質量許容差"], de: ["(optional)", "Formel", "Massentoleranz"], fr: ["(facultatif)", "Formule", "Tolérance"], es: ["(opcional)", "Fórmula", "Tolerancia"],
  };
  for (const [lang, [suffix, formula, tolerance]] of Object.entries(pages)) {
    setDictionary(lang, { "(optional)": suffix });
    setLang(lang);
    for (const gap of ["", " "]) {
      const fields = [{ type: "text", id: "f", label: `${formula}${gap}${suffix}` }, { type: "number", id: "tol", label: `${tolerance} (ppm)` }, { type: "checkbox", id: "c", label: `${formula}${gap}${suffix}` }];
      assert.deepEqual(items(fields, values), [`${formula}: C6H6`, `${tolerance}: 5`, formula], `${lang}, gap "${gap}"`);
    }
    assert.deepEqual(items(english, values), ["Formula: C6H6", "Mass tolerance: 5", "Charge states"], `${lang}: a label still in English loses its English suffix`);
    setLang("en");
  }
  setDictionary("zh", { "(optional)": "（可选）" });
  setLang("zh");
  assert.deepEqual(items([{ type: "text", id: "f", label: "式（可选）" }], { f: "（可选）" }), ["式: （可选）"], "only the label is cut, never the value");
});

test("the Session report is written in the language of the page: <html lang> follows the UI language", () => {
  const report = () => buildHtml({ version: "1.2.3" }, "run-1", { human: "now" }, [], [], "", "", "light");
  assert.match(report(), /^<!doctype html>\n<html lang="en" data-theme="light">/);
  for (const [lang, tag] of [["zh", "zh-CN"], ["ja", "ja"], ["de", "de"], ["fr", "fr"], ["es", "es"]]) {
    setLang(lang);
    assert.match(report(), new RegExp(`^<!doctype html>\\n<html lang="${tag}" data-theme="light">`), lang);
  }
  setLang("en");
  assert.match(report(), /<html lang="en" /);
});

test("the input bar leaves off the bracket at the end of a label, ASCII or full-width, and keeps one in the middle", () => {
  const label = text => summarizeInputs([{ type: "number", id: "n", label: text }], { n: 5 }).items[0];
  assert.equal(label("Mass tolerance (ppm)"), "Mass tolerance: 5");
  assert.equal(label("質量許容差（ppm）"), "質量許容差: 5");
  assert.equal(label("质量容差 （ppm，仅 MS/MS）"), "质量容差: 5");
  assert.equal(label("Précision (Da) (facultatif)"), "Précision: 5");
  assert.equal(label("Tolérance (ppm) maximale"), "Tolérance (ppm) maximale: 5", "a bracket that is not at the end stays");
  assert.equal(label("質量（ppm）の許容差"), "質量（ppm）の許容差: 5");
});

test("the Session report's own text is in the language of the page; ids, dates, numbers and the product name are not translated", () => {
  const entry = { kind: "module", id: "none", title: "Flow", section: "Workflows", finishedAt: "2026-10-31T15:31:12-06:00", seconds: 2, inputs: { text: "x".repeat(1600) } };
  const report = (pages = [null]) => buildHtml({ version: "1.2.3", modules: [], pipelines: [] }, "run-1", { human: "2026-10-31 15:30:12 MDT" }, [entry], pages, "", "", "light");

  const english = report();
  for (const piece of [
    "<title>mzDesk report run-1</title>", "<h1>mzDesk report</h1>", "<span>Run ID</span><b>run-1</b>", "<span>Generated</span><b>2026-10-31 15:30:12 MDT</b>",
    "<span>Time zone</span><b>America/Edmonton (Edmonton time)</b>", "<span>mzDesk</span><b>V1.2.3</b>", "<span>Pages</span><b>1</b>",
    '<nav class="rp-toc" aria-label="Pages in this report"><b>Pages</b>', "<span>finished 2026-10-31 15:31:12-06:00</span><span>2.0 s</span>",
    "<summary>Inputs as sent (JSON)</summary>", "… (+100 characters, full text in report.json)",
    '<p class="empty-hint">This page could not be drawn for the report; its result is in report.json.</p>',
    '<p class="rp-foot">Made by mzDesk in the browser from the results on screen. report.json in the same folder has the inputs and the raw results.</p>',
  ]) assert.ok(english.includes(piece), piece);

  setDictionary("de", {
    "mzDesk report": "mzDesk-Bericht", "Run ID": "Lauf-ID", Generated: "Erstellt", "Time zone": "Zeitzone", "{zone} (Edmonton time)": "{zone} (Zeit von Edmonton)", Pages: "Seiten",
    "Pages in this report": "Seiten in diesem \"Bericht\"", "finished {when}": "fertig {when}", "Inputs as sent (JSON)": "Eingaben wie gesendet (JSON)",
    "… (+{n} characters, full text in report.json)": "… (+{n} Zeichen, ganzer Text in report.json)",
    "This page could not be drawn for the report; its result is in report.json.": "Diese Seite <konnte> nicht gezeichnet werden.",
    "Made by mzDesk in the browser from the results on screen. report.json in the same folder has the inputs and the raw results.": "Von mzDesk im Browser erstellt.",
  });
  setLang("de");
  const german = report();
  for (const piece of [
    "<html lang=\"de\" ", "<title>mzDesk-Bericht run-1</title>", "<h1>mzDesk-Bericht</h1>", "<span>Lauf-ID</span><b>run-1</b>", "<span>Erstellt</span><b>2026-10-31 15:30:12 MDT</b>",
    "<span>Zeitzone</span><b>America/Edmonton (Zeit von Edmonton)</b>", "<span>mzDesk</span><b>V1.2.3</b>", "<span>Seiten</span><b>1</b>",
    '<nav class="rp-toc" aria-label="Seiten in diesem &quot;Bericht&quot;"><b>Seiten</b>', "<span>fertig 2026-10-31 15:31:12-06:00</span><span>2.0 s</span>",
    "<summary>Eingaben wie gesendet (JSON)</summary>", "… (+100 Zeichen, ganzer Text in report.json)",
    '<p class="empty-hint">Diese Seite &lt;konnte&gt; nicht gezeichnet werden.</p>', '<p class="rp-foot">Von mzDesk im Browser erstellt.</p>',
  ]) assert.ok(german.includes(piece), piece);
  setLang("en");
  assert.equal(report(), english, "back to English: the same file as before");
});
