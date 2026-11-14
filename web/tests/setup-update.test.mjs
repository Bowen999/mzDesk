import assert from "node:assert/strict";
import { after, beforeEach, test } from "node:test";
import { MemoryRouter } from "react-router-dom";
import { React, startSsr } from "./ssr-helper.mjs";

// The update notice: the sidebar footer item shows only for a newer version the person has not put off (setup.ts, components/UpdateNotice.tsx).
const ssr = await startSsr("setup-update");
const S = await ssr.load("/src/setup.ts");
const { UpdateCard, UpdateNotice } = await ssr.load("/src/components/UpdateNotice.tsx");
const { Sidebar } = await ssr.load("/src/components/Sidebar.tsx");
const { setDictionary, setLang } = await ssr.load("/src/i18n/index.ts");
after(async () => { setLang("en"); await ssr.close(); });
beforeEach(() => { S.resetUpdate(); setLang("en"); });

const info = (extra = {}) => ({ current: "0.5.4", latest: "0.5.5", newer: true, released: "2026-11-05", notes: "https://mzdesk.com/docs/changelog", command: "mzdesk update", checked_at: null, ...extra });
const notice = props => ssr.render(UpdateNotice, props);

test("visibleUpdate: only a newer version that was not put off", () => {
  assert.equal(S.visibleUpdate(null, null), null);
  assert.equal(S.visibleUpdate(info({ newer: false }), null), null);
  assert.equal(S.visibleUpdate(info({ latest: null }), null), null, "offline: no version");
  assert.equal(S.visibleUpdate(info(), "0.5.5"), null, "put off");
  assert.equal(S.visibleUpdate(info({ latest: "0.5.6" }), "0.5.5") !== null, true, "a newer one than the one put off shows again");
  assert.equal(S.visibleUpdate(info(), null).latest, "0.5.5");
});

test("the footer item: nothing until the server says newer; a quiet button with the version; the rail has an icon with a dot", () => {
  assert.equal(notice({}), "");
  S.setUpdateInfo(info({ newer: false, latest: "0.5.4" }));
  assert.equal(notice({}), "", "up to date");
  S.setUpdateInfo(info());
  const html = notice({});
  assert.match(html, /<button[^>]*class="internet-switch update-entry"[^>]*>/);
  assert.match(html, /aria-haspopup="dialog"[^>]*aria-expanded="false"/);
  assert.match(html, /<span class="internet-label">Version 0\.5\.5 is available<\/span>/);
  const rail = notice({ compact: true });
  assert.match(rail, /class="rail-btn update-rail"/);
  assert.match(rail, /aria-label="Version 0\.5\.5 is available"/);
  assert.match(rail, /class="update-dot"/, "the dot");
  assert.doesNotMatch(rail, /internet-label/, "icon only");
});

test("Later: remembered for that version; a newer one shows again", () => {
  S.setUpdateInfo(info());
  assert.notEqual(notice({}), "");
  S.dismissUpdate("0.5.5");
  assert.equal(notice({}), "");
  S.setUpdateInfo(info({ latest: "0.5.6" }));
  assert.match(notice({}), /Version 0\.5\.6 is available/);
});

test("the popover: the version, the changelog link inside the app, the command with Copy, Later", () => {
  const card = ssr.render(UpdateCard, { info: info() });
  assert.match(card, /<p class="update-title">Version 0\.5\.5 is available<\/p>/);
  assert.match(card, /<a class="update-notes" href="\/docs\/changelog">What changed<\/a>/, "no router here: a plain link to the same place");
  assert.match(card, /<code>mzdesk update<\/code><button type="button" class="copy">Copy<\/button>/);
  assert.match(card, /In a terminal:/);
  assert.match(card, /<button type="button" class="btn ghost">Later<\/button>/);
  const routed = ssr.render(() => React.createElement(MemoryRouter, null, React.createElement(UpdateCard, { info: info() })));
  assert.match(routed, /<a class="update-notes" href="\/docs\/changelog">What changed<\/a>/, "a router Link, not mzdesk.com");
  assert.doesNotMatch(routed, /mzdesk\.com/);
});

test("the notice is in the page's language", () => {
  setDictionary("de", { "Version {v} is available": "Version {v} ist verfügbar", Later: "Später", "What changed": "Was sich geändert hat" });
  setLang("de");
  S.setUpdateInfo(info());
  assert.match(notice({}), /Version 0\.5\.5 ist verfügbar/);
  const card = ssr.render(UpdateCard, { info: info() });
  assert.match(card, />Später</);
  assert.match(card, />Was sich geändert hat</);
});

test("the sidebar: the item sits in the footer of the local version, never on the public website", () => {
  S.setUpdateInfo(info());
  const spec = hosted => ({ groups: [], subgroups: [], pipelines: [], modules: [], hosted });
  const side = (hosted, collapsed) => ssr.render(() => React.createElement(MemoryRouter, null, React.createElement(Sidebar, { spec: spec(hosted), collapsed })));
  assert.match(side(false, false), /update-entry/);
  assert.match(side(false, true), /update-rail/);
  assert.doesNotMatch(side(true, false), /update-entry|Version 0\.5\.5/);
  assert.doesNotMatch(side(true, true), /update-rail/);
});
