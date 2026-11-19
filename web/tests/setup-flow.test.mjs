import assert from "node:assert/strict";
import { after, afterEach, beforeEach, test } from "node:test";
import { startSsr } from "./ssr-helper.mjs";

// ensureSetup (setup.ts): what a Run asks the server, when the dialog opens, and what each answer of the person does. The server is a stub of fetch.
const ssr = await startSsr("setup-flow");
const S = await ssr.load("/src/setup.ts");
const realFetch = globalThis.fetch;
after(async () => { globalThis.fetch = realFetch; await ssr.close(); });

let calls, routes;
/** routes: "METHOD /path" -> a body, a function (body) => body | Response-like, or an array answered in order (the last one repeats). */
const answer = (status, body) => ({ ok: status < 400, status, statusText: status < 400 ? "OK" : "Error", json: async () => body });
beforeEach(() => {
  calls = []; routes = {};
  S.setLocalApp(true);
  globalThis.fetch = async (url, init = {}) => {
    const key = `${init.method ?? "GET"} ${String(url).replace(/^https?:\/\/[^/]+/, "")}`;
    const body = init.body ? JSON.parse(init.body) : undefined;
    calls.push({ key, body });
    let r = routes[key];
    if (Array.isArray(r)) r = r.length > 1 ? r.shift() : r[0];
    if (r === undefined) return answer(404, { detail: `no route ${key}` });
    if (typeof r === "function") r = r(body);
    return r && "ok" in r && "json" in r ? r : answer(200, r);
  };
});
afterEach(() => { S.dismissSetup(); });

const until = async (check, what) => {
  for (let i = 0; i < 400; i++) { if (check()) return; await new Promise(r => setTimeout(r, 5)); }
  assert.fail(`timeout: ${what}`);
};
const phase = () => S.getSetupState().phase;
const NEEDS = "POST /api/local/needs", INSTALL = "POST /api/local/install", STATUS = "GET /api/local/status";
const lib = (ask = true) => ({ id: "libraries", title: "Spectral libraries", size: "2.3 GB", ask, why: "MassBank and GNPS spectra" });
const msf = (ask = true) => ({ id: "msfinder", title: "MS-FINDER", size: "0.6 GB", ask, why: "" });
const none = { missing: [], blocked: [], unsupported: [] };
const done = ids => ({ status: "done", current: null, elapsed_s: 3, components: ids.map(id => ({ id, title: id, status: "done", detail: "" })), log: [], error: "" });
const running = ids => ({ status: "running", current: ids[0], elapsed_s: 1, components: ids.map((id, i) => ({ id, title: id, status: i ? "waiting" : "running", detail: "x" })), log: [], error: "" });

test("hosted (or not told yet): true at once, the API is never called", async () => {
  S.setLocalApp(false);
  assert.equal(await S.ensureSetup({ target: "module:ms2_library", inputs: {} }), true);
  assert.equal(await S.ensureSetup({ ids: ["frigid"] }), true);
  assert.deepEqual(calls, []);
  assert.equal(phase(), null);
});

test("nothing missing: true, no dialog; the request carries the target and the form values", async () => {
  routes[NEEDS] = { ...none, unsupported: [{ id: "cfmid", title: "CFM-ID", message: "macOS / Linux only for now." }] };
  assert.equal(await S.ensureSetup({ target: "workflow:esi_ms2_annotation", inputs: { precursor: { precursor_mz: "195.0877" } } }), true);
  assert.deepEqual(calls, [{ key: NEEDS, body: { target: "workflow:esi_ms2_annotation", inputs: { precursor: { precursor_mz: "195.0877" } } } }]);
  assert.equal(phase(), null, "an unsupported component alone never opens the dialog");
});

test("the API fails (HTTP error, no answer, odd body): the run goes on without a dialog", async () => {
  routes[NEEDS] = answer(500, { detail: "boom" });
  assert.equal(await S.ensureSetup({ target: "module:x", inputs: {} }), true);
  globalThis.fetch = async () => { throw new TypeError("Failed to fetch"); };
  assert.equal(await S.ensureSetup({ target: "module:x", inputs: {} }), true);
  assert.equal(phase(), null);
  routes[NEEDS] = answer(404, { local: false });
  globalThis.fetch = async () => answer(404, { detail: "Not Found" });
  assert.equal(await S.ensureSetup({ target: "module:x", inputs: {} }), true);
});

test("something big is missing: the dialog asks; Not now answers false and closes; a double click does not open a second dialog", async () => {
  routes[NEEDS] = { ...none, missing: [lib(), msf(false)] };
  const first = S.ensureSetup({ target: "module:ms2_library", inputs: {} });
  assert.equal(await S.ensureSetup({ target: "module:ms2_library", inputs: {} }), false, "while one flow runs, the next answers false at once");
  await until(() => phase() === "confirm", "confirm view");
  assert.equal(await S.ensureSetup({ target: "module:ms2_library", inputs: {} }), false, "also while the dialog is open");
  assert.equal(S.getSetupState().forRun, true);
  assert.deepEqual(S.getSetupState().missing.map(m => m.id), ["libraries", "msfinder"]);
  S.dismissSetup();
  assert.equal(await first, false);
  assert.equal(phase(), null);
  assert.equal(calls.filter(c => c.key === NEEDS).length, 1, "the double click asked the server once");
  assert.equal(calls.some(c => c.key === INSTALL), false, "nothing was downloaded");
});

test("Download: install the missing ids, show progress, true when the job is done and the server agrees", async () => {
  routes[NEEDS] = [{ ...none, missing: [lib(), msf(false)] }, none];
  routes[INSTALL] = { job_id: "j1" };
  routes["GET /api/local/jobs/j1"] = [done(["libraries", "msfinder"])];
  const p = S.ensureSetup({ target: "workflow:esi_ms2_annotation", inputs: {} });
  await until(() => phase() === "confirm", "confirm");
  S.confirmSetup();
  assert.equal(await p, true);
  assert.deepEqual(calls.map(c => c.key), [NEEDS, INSTALL, "GET /api/local/jobs/j1", NEEDS]);
  assert.deepEqual(calls[1].body, { ids: ["libraries", "msfinder"] });
  assert.equal(phase(), null);
});

test("nothing asks (the MS/MS Annotation case): straight to progress, no confirm view", async () => {
  routes[NEEDS] = [{ ...none, missing: [lib(false), msf(false)] }, none];
  routes[INSTALL] = { job_id: "j2" };
  routes["GET /api/local/jobs/j2"] = [running(["libraries", "msfinder"]), done(["libraries", "msfinder"])];
  const p = S.ensureSetup({ target: "workflow:esi_ms2_annotation", inputs: {} });
  await until(() => phase() === "progress", "progress");
  assert.equal(S.getSetupState().forRun, true);
  await until(() => S.getSetupState().job?.status === "running", "first poll");
  assert.equal(S.getSetupState().seen.libraries > 0, true, "the bar's clock starts when a component is first seen running");
  assert.equal(await p, true);
});

test("Cancel while it downloads: the server is told, false, closed (Escape does the same)", async () => {
  routes[NEEDS] = { ...none, missing: [lib()] };
  routes[INSTALL] = { job_id: "j3" };
  routes["GET /api/local/jobs/j3"] = running(["libraries"]);
  routes["POST /api/local/jobs/j3/cancel"] = { ok: true };
  const p = S.ensureSetup({ target: "module:ms2_library", inputs: {} });
  await until(() => phase() === "confirm", "confirm");
  S.confirmSetup();
  await until(() => S.getSetupState().job?.status === "running", "running");
  S.dismissSetup();
  assert.equal(S.getSetupState().cancelling, true);
  assert.equal(await p, false);
  assert.ok(calls.some(c => c.key === "POST /api/local/jobs/j3/cancel"));
  assert.equal(phase(), null);
});

test("blocked (Docker is not running): the run is cancelled, nothing is downloaded, Close answers false", async () => {
  routes[NEEDS] = { ...none, missing: [lib()], blocked: [{ id: "cfmid", title: "CFM-ID", message: "Docker is not running.", link: "https://docs.docker.com/get-docker/" }] };
  const p = S.ensureSetup({ target: "module:ms2_esi", inputs: {} });
  await until(() => phase() === "blocked", "blocked");
  assert.equal(S.getSetupState().blocked[0].link, "https://docs.docker.com/get-docker/");
  S.dismissSetup();
  assert.equal(await p, false);
  assert.equal(calls.some(c => c.key === INSTALL), false);
});

test("a failed job: the error view with the message and the log tail; Retry installs again; Close answers false", async () => {
  routes[NEEDS] = [{ ...none, missing: [lib(false)] }, none];
  routes[INSTALL] = [{ job_id: "bad" }, { job_id: "good" }];
  routes["GET /api/local/jobs/bad"] = { status: "error", current: null, elapsed_s: 2, components: [{ id: "libraries", title: "Spectral libraries", status: "error", detail: "404" }], log: ["a", "b"], error: "The server answered 404." };
  routes["GET /api/local/jobs/good"] = done(["libraries"]);
  const p = S.ensureSetup({ target: "module:ms2_library", inputs: {} });
  await until(() => phase() === "error", "error");
  assert.equal(S.getSetupState().error, "The server answered 404.");
  assert.deepEqual(S.getSetupState().log, ["a", "b"]);
  S.retrySetup();
  assert.equal(await p, true);
  assert.deepEqual(calls.filter(c => c.key === INSTALL).map(c => c.body), [{ ids: ["libraries"] }, { ids: ["libraries"] }]);

  routes[INSTALL] = { job_id: "bad" };
  routes[NEEDS] = { ...none, missing: [lib(false)] };
  const q = S.ensureSetup({ target: "module:ms2_library", inputs: {} });
  await until(() => phase() === "error", "error again");
  S.dismissSetup();
  assert.equal(await q, false);
});

test("a job that is gone (the server restarted) ends in the error view with the server's words", async () => {
  routes[NEEDS] = { ...none, missing: [lib(false)] };
  routes[INSTALL] = { job_id: "gone" };
  routes["GET /api/local/jobs/gone"] = answer(404, { detail: "Job not found" });
  const p = S.ensureSetup({ target: "module:ms2_library", inputs: {} });
  await until(() => phase() === "error", "error");
  assert.equal(S.getSetupState().error, "Job not found");
  S.dismissSetup();
  assert.equal(await p, false);
});

test("done, but the server still misses the component: an error, never a run on a missing engine", async () => {
  routes[NEEDS] = { ...none, missing: [msf(false)] };
  routes[INSTALL] = { job_id: "j4" };
  routes["GET /api/local/jobs/j4"] = done(["msfinder"]);
  const p = S.ensureSetup({ target: "module:ms2_structure", inputs: {} });
  await until(() => phase() === "error", "error");
  assert.equal(S.getSetupState().error, "MS-FINDER could not be found after the download.");
  S.dismissSetup();
  assert.equal(await p, false);
});

test("a tick box asks for its ids through /api/local/status; installed or unsupported ones need nothing", async () => {
  const comp = (id, extra) => ({ id, title: id.toUpperCase(), size: "3 GB", installed: false, supported: true, ask: true, note: "", ...extra });
  routes[STATUS] = { local: true, platform: "darwin", components: [comp("frigid"), comp("neims", { installed: true }), comp("cfmid", { supported: false, note: "macOS / Linux only for now." })] };
  const p = S.ensureSetup({ ids: ["frigid"] });
  await until(() => phase() === "confirm", "confirm");
  assert.equal(S.getSetupState().forRun, false, "a tick box, not a run");
  assert.deepEqual(S.getSetupState().missing.map(m => m.id), ["frigid"]);
  S.dismissSetup();
  assert.equal(await p, false, "Not now: the box stays unticked");
  assert.equal(await S.ensureSetup({ ids: ["neims"] }), true);
  assert.equal(await S.ensureSetup({ ids: ["cfmid"] }), true, "unsupported never blocks");
  routes[STATUS] = answer(500, {});
  assert.equal(await S.ensureSetup({ ids: ["frigid"] }), true, "no answer: the box is ticked as before");
  assert.ok(calls.every(c => c.key === STATUS), "a tick box never calls /needs");
});
