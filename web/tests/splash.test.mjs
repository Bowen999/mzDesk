import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { test } from "node:test";
import { fileURLToPath } from "node:url";
import vm from "node:vm";

// The opening splash is plain HTML + two inline scripts in index.html (it plays while the app downloads). It plays once per browser session (tab):
// the first load sets a mark in sessionStorage, the later loads of that session drop the splash before the first paint.
const html = readFileSync(join(dirname(dirname(fileURLToPath(import.meta.url))), "index.html"), "utf8");
const [headScript, bodyScript] = [...html.matchAll(/<script>([\s\S]*?)<\/script>/g)].map(m => m[1]);

/** Run both inline scripts against a small fake page. `session` is the sessionStorage content (null = sessionStorage throws, as in some private modes). */
function visit(session) {
  const root = { className: "", setAttribute() {} };
  const splash = { removed: false, classes: [], listeners: [], remove() { this.removed = true; }, classList: { add: c => splash.classes.push(c) }, addEventListener(type) { splash.listeners.push(type); } };
  const timers = [];
  const store = session ?? {};
  const sandbox = {
    document: { documentElement: root, getElementById: id => (id === "splash" ? splash : null), addEventListener() {} },
    localStorage: { getItem: () => null },
    sessionStorage: session === null ? { getItem() { throw new Error("blocked"); }, setItem() { throw new Error("blocked"); } }
      : { getItem: k => store[k] ?? null, setItem: (k, v) => { store[k] = v; } },
    performance: { now: () => 500 },
    setTimeout: (fn, ms) => timers.push({ fn, ms }),
  };
  sandbox.window = sandbox;
  vm.createContext(sandbox);
  vm.runInContext(headScript, sandbox);
  vm.runInContext(bodyScript, sandbox);
  return { root, splash, timers, store, ready: () => sandbox.mzdeskReady() };
}

test("the first load of a session plays the splash and leaves a mark; the app's ready signal ends it after the animation has had its time", () => {
  const v = visit({});
  assert.doesNotMatch(v.root.className, /no-splash/);
  assert.equal(v.store["mzdesk.splash"], "1");
  assert.equal(v.splash.removed, false);
  assert.deepEqual(v.splash.listeners, ["click"]);
  assert.equal(v.timers.length, 1, "only the failsafe is waiting");
  assert.equal(v.timers[0].ms, 15000);
  v.ready();
  assert.equal(v.timers[1].ms, 1500, "2000 ms minimum, 500 ms already gone");
  v.timers[1].fn();
  assert.deepEqual(v.splash.classes, ["sp-out"]);
  v.timers[2].fn();
  assert.equal(v.splash.removed, true);
});

test("a later load of the same session skips the splash: marked in the page before it paints, removed at once, nothing left waiting", () => {
  const v = visit({ "mzdesk.splash": "1" });
  assert.match(v.root.className, /no-splash/);
  assert.equal(v.splash.removed, true);
  assert.equal(v.timers.length, 0);
  v.ready();
  assert.equal(v.timers.length, 0, "the ready signal is a no-op");
  assert.match(html, /html\.no-splash #splash \{ display: none; \}/, "and the stylesheet hides it if the page paints first");
});

test("where sessionStorage is blocked the splash plays on every load rather than never", () => {
  const v = visit(null);
  assert.doesNotMatch(v.root.className, /no-splash/);
  assert.equal(v.splash.removed, false);
});
