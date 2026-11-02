import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { after, afterEach, test } from "node:test";
import { fileURLToPath } from "node:url";
import { React, startSsr } from "./ssr-helper.mjs";

// The first-use dialog (components/SetupDialog.tsx): its four views, the sizes and the bar (setup.ts), and its words in the other languages.
const ssr = await startSsr("setup-dialog");
const { SetupView } = await ssr.load("/src/components/SetupDialog.tsx");
const S = await ssr.load("/src/setup.ts");
const { setDictionary, setLang } = await ssr.load("/src/i18n/index.ts");
after(() => ssr.close());
afterEach(() => setLang("en"));

const here = dirname(fileURLToPath(import.meta.url));
const missing = [
  { id: "libraries", title: "Spectral libraries", size: "2.3 GB", ask: true, why: "MassBank and GNPS spectra" },
  { id: "msfinder", title: "MS-FINDER", size: "0.6 GB", ask: true, why: "Structure candidates from MS/MS" },
];
const state = (phase, extra = {}) => ({
  phase, forRun: true, missing, blocked: [], unsupported: [], job: null, error: "", log: [], cancelling: false, seen: {}, ...extra,
});
const job = (components, extra = {}) => ({ status: "running", current: null, elapsed_s: 0, components, log: [], error: "", ...extra });
const view = (s, props = {}) => ssr.render(SetupView, { state: s, ...props });
const buttons = html => [...html.matchAll(/<button[^>]*>(.*?)<\/button>/g)].map(m => m[1]);

// ---------------------------------------------------------------- confirm
test("confirm: what will be downloaded with its sizes, the total on the button, Not now first", () => {
  const html = view(state("confirm", { unsupported: [{ id: "cfmid", title: "CFM-ID", message: "macOS / Linux only for now." }] }));
  assert.match(html, /role="dialog"/);
  assert.match(html, /aria-modal="true"/);
  assert.match(html, /aria-labelledby="setup-title"/);
  assert.match(html, /<h2 id="setup-title">Download needed<\/h2>/);
  assert.match(html, /mzDesk needs to download:/);
  assert.match(html, /<span class="setup-name">Spectral libraries<\/span><span class="setup-size">2\.3 GB<\/span><span class="setup-why">MassBank and GNPS spectra<\/span>/);
  assert.match(html, /<span class="setup-name">MS-FINDER<\/span><span class="setup-size">0\.6 GB<\/span>/);
  assert.deepEqual(buttons(html), ["Not now", "Download 2.9 GB"]);
  assert.match(html, /class="btn" data-autofocus/, "the primary button takes the focus");
  assert.match(html, /CFM-ID: macOS \/ Linux only for now\./, "an unsupported component is one line of the dialog");
  assert.doesNotMatch(html, /progressbar|setup-rows/);
});

test("confirm: a size that cannot be read leaves the total off the button", () => {
  const html = view(state("confirm", { missing: [{ ...missing[0], size: "see code" }] }));
  assert.deepEqual(buttons(html), ["Not now", "Download"]);
});

// ---------------------------------------------------------------- progress
test("progress: one row per component, the bar, the clock and Cancel; the state is said in words too", () => {
  const s = state("progress", {
    job: job([
      { id: "libraries", title: "Spectral libraries", status: "done", detail: "" },
      { id: "msfinder", title: "MS-FINDER", status: "running", detail: "Unpacking massbank/…" },
      { id: "molgpka", title: "MolGpKa", status: "waiting", detail: "" },
    ], { elapsed_s: 125 }),
    missing: [...missing, { id: "molgpka", title: "MolGpKa", size: "12 MB", ask: false, why: "" }],
  });
  const html = view(s, { percent: 61.4 });
  assert.match(html, /<h2 id="setup-title">First-time setup<\/h2>/);
  assert.match(html, /Keep this page open\. The run starts when this is done\./);
  assert.match(html, /role="progressbar"[^>]*aria-valuemin="0"[^>]*aria-valuemax="100"[^>]*aria-valuenow="61"/);
  assert.match(html, /aria-valuetext="Installing MS-FINDER \(2 of 3\)"/);
  assert.match(html, /<span style="width:61\.4%"><\/span>/);
  assert.match(html, /<span class="setup-clock">2:05<\/span>/);
  assert.match(html, /<li class="setup-row is-done">.*?Spectral libraries<span class="setup-size">2\.3 GB<\/span>.*?<span class="setup-word">Done<\/span><\/li>/);
  assert.match(html, /<li class="setup-row is-running">.*?<span class="setup-word">Installing<\/span><span class="setup-detail" title="Unpacking massbank\/…">Unpacking massbank\/…<\/span><\/li>/);
  assert.match(html, /<li class="setup-row is-waiting">.*?MolGpKa<span class="setup-size">12 MB<\/span>.*?<span class="setup-word">Waiting<\/span><\/li>/);
  assert.match(html, /<p class="sr-only" aria-live="polite">Installing MS-FINDER \(2 of 3\)<\/p>/, "one polite line for screen readers: the component that changed, not the log");
  assert.deepEqual(buttons(html), ["Cancel"]);
  assert.match(html, /tabindex="-1"/, "the card takes the focus, so a stray Enter cannot cancel");
});

test("progress: before the server answers every missing component waits; a tick box has no 'run starts' line; Cancelling… is disabled", () => {
  const html = view(state("progress", { forRun: false, cancelling: true }), { percent: 6 });
  assert.equal((html.match(/class="setup-row is-waiting"/g) ?? []).length, 2);
  assert.match(html, /Keep this page open\.<\/p>/);
  assert.doesNotMatch(html, /run starts/);
  assert.match(html, /<button type="button" class="btn secondary" disabled="">Cancelling…<\/button>/);
});

test("progress: a failed component shows its last line in the alert colour class", () => {
  const html = view(state("progress", { job: job([{ id: "libraries", title: "Spectral libraries", status: "error", detail: "curl: (22) 404" }]) }));
  assert.match(html, /<li class="setup-row is-error">.*?<span class="setup-word">Failed<\/span><span class="setup-detail" title="curl: \(22\) 404">curl: \(22\) 404<\/span>/);
});

// ---------------------------------------------------------------- blocked and error
test("blocked: the message and where to fix it, only Close; a link that is not http(s) is not drawn", () => {
  const b = { id: "cfmid", title: "CFM-ID", message: "Docker is not running. Start Docker Desktop, then run again.", link: "https://docs.docker.com/get-docker/" };
  const html = view(state("blocked", { blocked: [b], missing: [] }));
  assert.match(html, /<h2 id="setup-title">Cannot run yet<\/h2>/);
  assert.match(html, /<span class="setup-name">CFM-ID<\/span>/);
  assert.match(html, /Docker is not running\. Start Docker Desktop, then run again\./);
  assert.match(html, /<a class="setup-link" href="https:\/\/docs\.docker\.com\/get-docker\/" target="_blank" rel="noopener noreferrer">docs\.docker\.com\/get-docker<\/a>/);
  assert.deepEqual(buttons(html), ["Close"]);
  const unsafe = view(state("blocked", { blocked: [{ ...b, link: "javascript:alert(1)" }, { ...b, id: "x", link: null }], missing: [] }));
  assert.doesNotMatch(unsafe, /<a |javascript/);
});

test("error: the message as an alert, the log behind 'Tool log', Close and Retry", () => {
  const html = view(state("error", { error: "Could not download FRIGID.", log: ["Downloading", "curl: (22) 404"] }));
  assert.match(html, /<h2 id="setup-title">Setup failed<\/h2>/);
  assert.match(html, /<div class="alert error" role="alert">Could not download FRIGID\.<\/div>/);
  assert.match(html, /<summary>Tool log<\/summary><pre class="code">Downloading\ncurl: \(22\) 404<\/pre>/);
  assert.deepEqual(buttons(html), ["Close", "Retry"]);
  assert.doesNotMatch(view(state("error", { error: "x" })), /Tool log/, "no log, no disclosure");
});

test("no phase, no dialog", () => {
  assert.equal(view(state(null)), "");
});

// ---------------------------------------------------------------- sizes and the bar
test("sizes: read from the server's text; the total is written as the dialog shows it", () => {
  assert.equal(S.sizeGb("2.3 GB"), 2.3);
  assert.ok(Math.abs(S.sizeGb("600 MB") - 0.6) < 1e-9);
  assert.ok(Math.abs(S.sizeGb("12 MB") - 0.012) < 1e-9);
  assert.equal(S.sizeGb("2,3 GB"), 2.3);
  assert.equal(S.sizeGb("see code"), null);
  assert.equal(S.sizeGb(""), null);
  assert.equal(S.totalSize(missing), "2.9 GB");
  assert.equal(S.totalSize([{ size: "12 MB" }, { size: "8 MB" }]), "20 MB");
  assert.equal(S.totalSize([{ size: "2.3 GB" }, { size: "?" }]), null);
  assert.equal(S.totalSize([]), null);
  assert.equal(S.clock(0), "0:00");
  assert.equal(S.clock(65.9), "1:05");
  assert.equal(S.clock(-3), "0:00");
});

test("the bar: never empty, never 100 before the job ends, fast at first, slower later, a finished component fills its slice", () => {
  const row = (id, size, status) => ({ id, title: id, size, status, detail: "" });
  const t0 = 1_000_000;
  const percent = (rows, seconds) => S.setupPercent(rows, { a: t0, b: t0 }, t0 + seconds * 1000);
  assert.equal(percent([row("a", "2.3 GB", "waiting")], 0), 6, "waiting: the floor");
  assert.equal(percent([], 0), 6);
  const running = [row("a", "2.3 GB", "running")];
  const seconds = [0, 1, 5, 20, 60, 200, 1000, 100000];
  const series = seconds.map(s => percent(running, s));
  assert.ok(series.every((p, i) => p >= 6 && p <= 99 && (i === 0 || p >= series[i - 1])), `bounded and never backwards: ${series}`);
  const rate = (a, b) => (series[b] - series[a]) / (seconds[b] - seconds[a]);
  assert.ok(rate(2, 3) > rate(4, 5) && rate(4, 5) > rate(5, 6), `fast at first, slower later: ${series}`);
  assert.ok(series[7] < 95, "a running component stops short of its end");
  const small = percent([row("a", "12 MB", "running")], 20), big = percent([row("a", "2.3 GB", "running")], 20);
  assert.ok(small > big, "a small download fills faster than a big one");
  const half = percent([row("a", "1 GB", "done"), row("b", "1 GB", "waiting")], 0);
  assert.ok(Math.abs(half - 50) < 1e-9, `the slices are as wide as the downloads: ${half}`);
  assert.equal(percent([row("a", "1 GB", "done"), row("b", "1 GB", "done")], 0), 99, "all done, the job not yet: 99");
  assert.ok(percent([row("a", "1 GB", "done"), row("b", "1 GB", "running")], 30) > half);
  assert.ok(percent([row("a", "9 GB", "running"), row("b", "?", "waiting")], 5) < 100);
});

// ---------------------------------------------------------------- the words
test("the views are written in the page's language: the five dictionaries carry every new string", () => {
  const words = lang => JSON.parse(readFileSync(join(here, `../src/i18n/locales/${lang}.json`), "utf8"));
  const expected = { zh: "需要下载", ja: "ダウンロードが必要です", de: "Download erforderlich", fr: "Téléchargement nécessaire", es: "Descarga necesaria" };
  for (const [lang, title] of Object.entries(expected)) {
    setDictionary(lang, words(lang));
    setLang(lang);
    const html = view(state("confirm"));
    assert.ok(html.includes(`<h2 id="setup-title">${title}</h2>`), `${lang}: ${title}`);
    assert.ok(!html.includes(">Not now<") && !html.includes("Download 2.9 GB"), `${lang}: buttons are translated`);
    assert.ok(html.includes("2.9 GB"), `${lang}: the size is kept`);
    for (const phase of ["progress", "blocked", "error"]) assert.doesNotMatch(view(state(phase, { blocked: [], error: "x" })), />(Cancel|Close|Retry|Setup failed|First-time setup|Cannot run yet)</, `${lang}: ${phase}`);
  }
});
