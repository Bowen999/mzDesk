import assert from "node:assert/strict";
import { after, afterEach, mock, test } from "node:test";
import { startSsr } from "./ssr-helper.mjs";

// Session report, report.json format 2: the versions of what produced the results (GET /api/versions) at the top, and the versions of the
// things each page is based on inside every page entry. report.html is not touched by any of it.
const ssr = await startSsr("report-environment");
const { basedOnEntries, buildHtml, buildJson, loadEnvironment } = await ssr.load("/src/report/build.ts");
const { localizeSpec, setDictionary, setLang } = await ssr.load("/src/i18n/index.ts");
const realFetch = globalThis.fetch;
after(() => ssr.close());
afterEach(() => { globalThis.fetch = realFetch; mock.timers.reset(); setLang("en"); });

const entry = (version, kind, detail = null) => ({ version, kind, detail });
const environment = {
  mzdesk: "0.5.2", python: "3.12.12", checked_at: "2026-11-02T06:35:55+00:00",
  components: {
    mzDesk: entry("0.5.2", "mzdesk"),
    BUDDY: entry("0.3.12", "python-package", "Python package msbuddy"),
    "MS-FINDER": entry("3.73", "engine", "RIKEN MS-FINDER console build; structure DB VS15, formula DB VS13; Mono 6.14.1"),
    "MassBank/GNPS": entry("MassBank 2026.03; GNPS 2026-03-29", "data", "MassBank: 126,507 spectra; GNPS: 955,289 spectra"),
    MassBank: entry("2026.03", "data", "126,507 spectra, built 2026-10-30"),
    GNPS: entry("2026-03-29", "data", "955,289 spectra, built 2026-10-30"),
    "GNPS (online)": entry(null, "online", "online service, not versioned"),
    fastMASST: entry(null, "online", "online service, not versioned"),
    "Pan-ReDU": entry(null, "online", "online service, not versioned"),
    matchms: entry("0.33.1", "python-package"),
    FRIGID: entry(null, "engine", "not installed"),
  },
};

const spec = {
  version: "0.5.2",
  modules: [
    { id: "gnps_fastmasst", title: "fastMASST", based_on: ["fastMASST", "GNPS", "Pan-ReDU"] },
    { id: "ms2_library", title: "Library search", based_on: ["MassBank", "GNPS", "matchms"] },
    { id: "plain", title: "Plain page" },
  ],
  pipelines: [
    { id: "esi_ms2_annotation", title: "MS/MS Annotation", based_on: ["BUDDY", "MS-FINDER", "FRIGID", "MassBank/GNPS"] },
    { id: "gnps_fastmasst", title: "Not the module of that id", based_on: ["PubChem"] },
  ],
};
const run = (kind, id, finishedAt, extra = {}) => ({
  key: `${kind}:${id}`, kind, id, title: id, section: "S", inputs: { a: 1 }, result: { r: id }, seconds: 1.5, finishedAt, ...extra,
});
const entries = [
  run("module", "gnps_fastmasst", "2026-11-02T10:00:00-07:00"),
  run("module", "ms2_library", "2026-11-02T10:05:00-07:00"),
  run("workflow", "esi_ms2_annotation", "2026-11-02T10:10:00-07:00", { steps: [{ n: 1 }] }),
  run("module", "plain", "2026-11-02T10:15:00-07:00"),
];
const when = { iso: "2026-11-02T10:20:00-07:00", abbr: "MST", offset: "-07:00", human: "now" };

// ---------------------------------------------------------------- report.json
test("report.json is format 2: the whole answer of /api/versions at the top, the versions of each page's based_on inside the page", () => {
  const json = buildJson(spec, "run-1", when, entries, { environment });
  assert.equal(json.report_format, 2);
  assert.deepEqual(Object.keys(json), ["report_format", "run_id", "generated_at", "timezone", "timezone_abbreviation", "utc_offset", "app", "environment", "pages"]);
  assert.deepEqual(json.environment, environment);
  assert.deepEqual(json.app, { name: "mzDesk", version: "0.5.2" });
  assert.equal("environment_error" in json, false);

  const [fast, library, flow, plain] = json.pages;
  assert.deepEqual(flow.based_on, [
    { name: "BUDDY", version: "0.3.12", kind: "python-package", detail: "Python package msbuddy" },
    { name: "MS-FINDER", ...environment.components["MS-FINDER"] },
    { name: "FRIGID", version: null, kind: "engine", detail: "not installed" },
    { name: "MassBank/GNPS", ...environment.components["MassBank/GNPS"] },
  ], "a workflow reads the workflow's based_on (not a module of the same id), and a missing version stays null");
  assert.deepEqual(plain.based_on, [], "a page without based_on gets an empty list");
  assert.deepEqual(library.based_on.map(b => b.name), ["MassBank", "GNPS", "matchms"]);
  assert.equal(library.based_on[1].version, "2026-03-29", "GNPS next to MassBank is the local library");
  assert.deepEqual(fast.based_on.map(b => b.name), ["fastMASST", "GNPS", "Pan-ReDU"]);
  assert.equal(fast.based_on[1].kind, "online", "GNPS on a page without MassBank is the online service");
  assert.equal(fast.based_on[1].version, null);
  assert.equal(fast.based_on[1].name, "GNPS", "the name is the one the page lists");
});

test("an online service gets accessed_at = the page's finished_at; nothing else does", () => {
  const json = buildJson(spec, "run-1", when, entries, { environment });
  const [fast, library, flow] = json.pages;
  for (const b of fast.based_on) assert.equal(b.accessed_at, "2026-11-02T10:00:00-07:00", b.name);
  assert.ok(library.based_on.every(b => !("accessed_at" in b)));
  assert.ok(flow.based_on.every(b => !("accessed_at" in b)));
});

test("everything else of a page entry is kept as it was", () => {
  const json = buildJson(spec, "run-1", when, entries, { environment });
  const flow = json.pages[2];
  assert.deepEqual(Object.keys(flow), ["kind", "id", "title", "section", "finished_at", "seconds", "based_on", "inputs", "result", "steps"]);
  assert.deepEqual({ kind: flow.kind, id: flow.id, title: flow.title, section: flow.section, finished_at: flow.finished_at, seconds: flow.seconds, inputs: flow.inputs, result: flow.result, steps: flow.steps },
    { kind: "workflow", id: "esi_ms2_annotation", title: "esi_ms2_annotation", section: "S", finished_at: "2026-11-02T10:10:00-07:00", seconds: 1.5, inputs: { a: 1 }, result: { r: "esi_ms2_annotation" }, steps: [{ n: 1 }] });
  assert.equal("steps" in json.pages[0], false, "steps only where the page has them");
});

test("when the versions could not be had, report.json says so and every page still lists its based_on names", () => {
  const json = buildJson(spec, "run-1", when, entries, { environment: null, error: "502 Bad Gateway" });
  assert.equal(json.environment, null);
  assert.equal(json.environment_error, "502 Bad Gateway");
  assert.equal(json.report_format, 2);
  assert.deepEqual(json.pages[2].based_on, ["BUDDY", "MS-FINDER", "FRIGID", "MassBank/GNPS"].map(name => ({ name, version: null, kind: null, detail: "unknown" })));
  assert.ok(json.pages.flatMap(p => p.based_on).every(b => !("accessed_at" in b)));
  JSON.parse(JSON.stringify(json));
});

test("a name the server does not know is listed as unknown, never dropped", () => {
  const page = { id: "x", based_on: ["Mystery", "GNPS"] };
  const out = basedOnEntries({ modules: [page], pipelines: [] }, run("module", "x", "2026-11-02T10:00:00-07:00"), { ...environment, components: {} });
  assert.deepEqual(out, [{ name: "Mystery", version: null, kind: null, detail: "unknown" }, { name: "GNPS", version: null, kind: null, detail: "unknown" }]);
});

test("the names come from the page list in any language: based_on is not display text", () => {
  setDictionary("de", { "MS/MS Annotation": "MS/MS-Annotation", "Library search": "Bibliotheksuche" });
  setLang("de");
  const german = localizeSpec(spec);
  assert.equal(german.pipelines[0].title, "MS/MS-Annotation", "the title is translated");
  assert.deepEqual(german.pipelines[0].based_on, spec.pipelines[0].based_on);
  const json = buildJson(german, "run-1", when, entries, { environment });
  assert.deepEqual(json.pages[2].based_on.map(b => b.name), ["BUDDY", "MS-FINDER", "FRIGID", "MassBank/GNPS"]);
  assert.equal(json.pages[2].based_on[0].version, "0.3.12");
});

// ---------------------------------------------------------------- report.html
test("report.html does not change: it takes no versions and shows none", () => {
  assert.equal(buildHtml.length, 8, "same arguments as before");
  const html = buildHtml({ version: "1.2.3", modules: [], pipelines: [] }, "run-1", { human: "now" }, [], [], "", "", "light");
  assert.doesNotMatch(html, /environment|based_on|accessed_at|report_format|0\.5\.2/);
});

// ---------------------------------------------------------------- asking the server
const answer = (status, body) => async () => new Response(typeof body === "string" ? body : JSON.stringify(body), { status, statusText: status === 200 ? "OK" : "Bad Gateway" });

test("loadEnvironment returns the answer of /api/versions", async () => {
  let seen;
  globalThis.fetch = async (url, init) => { seen = { url: String(url), signal: init?.signal }; return answer(200, environment)(); };
  assert.deepEqual(await loadEnvironment(), { environment });
  assert.match(seen.url, /\/api\/versions$/);
  assert.ok(seen.signal, "it can be given up on");
});

test("when the request fails the report is still made: environment null and a note", async () => {
  globalThis.fetch = answer(502, { detail: "upstream gone" });
  assert.deepEqual(await loadEnvironment(), { environment: null, error: "upstream gone" });
  globalThis.fetch = async () => { throw new TypeError("Failed to fetch"); };
  assert.deepEqual(await loadEnvironment(), { environment: null, error: "Failed to fetch" });
  globalThis.fetch = answer(200, "<!doctype html><title>mzDesk</title>");       // an older server answers every unknown path with the page itself
  const html = await loadEnvironment();
  assert.equal(html.environment, null);
  assert.ok(html.error);
  globalThis.fetch = answer(200, { components: null });
  assert.deepEqual(await loadEnvironment(), { environment: null, error: "unexpected answer from /api/versions" });
  globalThis.fetch = answer(200, []);
  assert.equal((await loadEnvironment()).environment, null);
  const json = buildJson(spec, "run-1", when, entries, html);
  assert.equal(json.environment, null);
  assert.ok(json.environment_error);
});

test("a server that does not answer is given up on after 20 s", async () => {
  mock.timers.enable({ apis: ["setTimeout"] });
  globalThis.fetch = (url, init) => new Promise((_, reject) => init.signal.addEventListener("abort", () => reject(new DOMException("aborted", "AbortError"))));
  const pending = loadEnvironment();
  mock.timers.tick(19999);
  let settled = false;
  pending.then(() => { settled = true; });
  await new Promise(resolve => setImmediate(resolve));
  assert.equal(settled, false, "still waiting at 19.999 s");
  mock.timers.tick(1);
  assert.deepEqual(await pending, { environment: null, error: "the server did not answer in time" });
});
