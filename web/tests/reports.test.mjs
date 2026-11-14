import assert from "node:assert/strict";
import { after, test } from "node:test";
import { startSsr } from "./ssr-helper.mjs";

// The crash screen and the maintainer's reports page.
const ssr = await startSsr("reports");
const { CrashBoundary } = await ssr.load("/src/components/CrashBoundary.tsx");
const { ReportsPage } = await ssr.load("/src/pages/ReportsPage.tsx");
after(() => ssr.close());

test("a page that draws fine is shown as it is", () => {
  const html = ssr.render(CrashBoundary, { label: "Identity", children: "inside" });
  assert.equal(html, "inside");
});

test("a render error turns into the crash screen with a way back", () => {
  assert.equal(CrashBoundary.getDerivedStateFromError(new Error("boom")).error.message, "boom");
  const boundary = new CrashBoundary({ label: "Identity", children: null });
  boundary.state = { error: new Error("boom"), ref: "E-1A2B3C" };
  const html = ssr.render(() => boundary.render());
  assert.match(html, /<h1>Identity stopped<\/h1>/);
  assert.match(html, /sent to the maintainer \(<code>E-1A2B3C<\/code>\)/);
  assert.match(html, />Try again<\/button>/);
  assert.match(html, />Reload<\/button>/);
  boundary.state = { error: new Error("boom"), ref: "" };
  assert.match(ssr.render(() => boundary.render()), /sent to the maintainer\./, "no reference yet, no parentheses");
});

test("the reports page starts with the password form and nothing else", () => {
  const html = ssr.render(ReportsPage, {});
  assert.match(html, /<h1>Reports<\/h1>/);
  assert.match(html, /type="password"/);
  assert.match(html, /Show reports/);
  assert.doesNotMatch(html, /report-list/);
});
