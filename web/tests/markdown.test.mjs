import assert from "node:assert/strict";
import { after, test } from "node:test";
import { MemoryRouter } from "react-router-dom";
import { React, startSsr } from "./ssr-helper.mjs";

// The Docs markdown (components/Markdown.tsx), the figure loader (docs/figure.ts) and the heading scan that "On this page" uses (docs/outline.ts).
const ssr = await startSsr("markdown");
const { Markdown } = await ssr.load("/src/components/Markdown.tsx");
const figures = await ssr.load("/src/docs/figure.ts");
const outline = await ssr.load("/src/docs/outline.ts");
after(() => ssr.close());

const h = React.createElement;
const md = text => ssr.render(() => h(MemoryRouter, null, h(Markdown, { text })));
const count = (html, re) => (html.match(re) ?? []).length;

// ---------------------------------------------------------------- tables
test("a pipe table has a header row, body rows and a frame that scrolls; alignment marks change nothing", () => {
  const html = md("| Tool | Role | Size |\n|:--|:-:|--:|\n| `BUDDY` | **formula** | 1 MB |\n| MS-FINDER | structure | 0.6 GB |");
  assert.match(html, /<div class="doc-table-wrap" tabindex="0"><table><thead><tr><th scope="col">Tool<\/th><th scope="col">Role<\/th><th scope="col">Size<\/th><\/tr><\/thead>/);
  assert.match(html, /<tbody><tr><td><code>BUDDY<\/code><\/td><td><b>formula<\/b><\/td><td>1 MB<\/td><\/tr><tr><td>MS-FINDER<\/td><td>structure<\/td><td>0\.6 GB<\/td><\/tr><\/tbody>/);
  assert.doesNotMatch(html, /align|text-align|style=/);
});

test("table edges: no outer pipes, an escaped pipe, a short row, and a table ends at the blank line", () => {
  const html = md("A | B\n--- | ---\nx \\| y | 1\nonly one |\n\nafter the table");
  assert.match(html, /<th scope="col">A<\/th><th scope="col">B<\/th>/);
  assert.match(html, /<td>x \| y<\/td><td>1<\/td>/);
  assert.match(html, /<td>only one<\/td><td><\/td>/, "a short row is padded");
  assert.match(html, /<\/table><\/div><p>after the table<\/p>/);
});

test("a line with a pipe is not a table without the separator row, and a rule under it stays a rule", () => {
  assert.equal(count(md("a | b\nc | d"), /<table/g), 0);
  const html = md("a | b\n\n---\n\nnext");
  assert.equal(count(html, /<table/g), 0);
  assert.match(html, /<p>a \| b<\/p><hr\/>/);
});

// ---------------------------------------------------------------- code blocks
test("a fenced block keeps its text literally, names its language, and the shell kinds get a Copy button", () => {
  const html = md("```bash\ncurl -fsSL https://mzdesk.com/install.sh | sh\n# not a heading <b>\n```");
  assert.match(html, /<pre tabindex="0"><code class="language-bash">curl -fsSL https:\/\/mzdesk\.com\/install\.sh \| sh\n# not a heading &lt;b&gt;<\/code><\/pre>/);
  assert.match(html, /<button type="button" class="doc-copy" aria-label="Copy code">Copy<\/button>/);
  assert.doesNotMatch(html, /<h1|<b>/, "no heading and no tag came out of the code");
  for (const lang of ["sh", "powershell", "text"]) assert.match(md("```" + lang + "\nx\n```"), /class="doc-copy"/, lang);
  for (const lang of ["json", "python", ""]) assert.doesNotMatch(md("```" + lang + "\nx\n```"), /doc-copy/, lang || "no language");
});

test("a fence ends at its own marker only; ~~~ works too; an unclosed fence runs to the end", () => {
  const html = md("~~~text\n```\ninside\n```\n~~~\n\nafter");
  assert.match(html, /<code class="language-text">```\ninside\n```<\/code>/);
  assert.match(html, /<p>after<\/p>/);
  assert.match(md("```text\nopen\n\nstill code"), /<code class="language-text">open\n\nstill code<\/code>/);
});

test("headings inside a code block get no id and are not in 'On this page'; the ids outside stay in step", () => {
  const text = "# Page\n\n## Install\n\n```bash\n# Install\n## Install\n```\n\n## Install\n\n### Run";
  assert.deepEqual(outline.headingsOf(text), [
    { id: "install", text: "Install", level: 2 }, { id: "install-1", text: "Install", level: 2 }, { id: "run", text: "Run", level: 3 }]);
  const html = md(text);
  assert.deepEqual([...html.matchAll(/<h[123] id="([^"]+)"/g)].map(m => m[1]), ["page", "install", "install-1", "run"]);
  assert.equal(outline.hasTitle(text), true);
  assert.equal(outline.hasTitle("## A\n\n```bash\n# comment\n```"), false, "a # line in code is not the page title");
});

// ---------------------------------------------------------------- lists, callouts, rules, inline
test("an ordered list is an <ol>; numbering that starts elsewhere is kept; a blank line between items does not split it", () => {
  assert.match(md("1. one\n2. two\n3. three"), /^<div class="doc-body"><ol><li>one<\/li><li>two<\/li><li>three<\/li><\/ol><\/div>$/);
  assert.match(md("3. c\n4. d"), /<ol start="3">/);
  const html = md("1. one\n\n2. two\n\n3. three\n\nend");
  assert.equal(count(html, /<ol/g), 1);
  assert.match(html, /<\/ol><p>end<\/p>/);
});

test("lists nest by indent and a wrapped item continues; the old flat bullets are unchanged", () => {
  const html = md("- a\n  - a1\n  - a2\n- b\n  wrapped\n\n1. step\n   - detail\n2. next");
  assert.match(html, /<ul><li>a<ul><li>a1<\/li><li>a2<\/li><\/ul><\/li><li>b wrapped<\/li><\/ul>/);
  assert.match(html, /<ol><li>step<ul><li>detail<\/li><\/ul><\/li><li>next<\/li><\/ol>/);
  assert.match(md("- **Label**: help with `code`\n- second"), /<ul><li><b>Label<\/b>: help with <code>code<\/code><\/li><li>second<\/li><\/ul>/);
});

test("> **Note** and > **Warning** are callouts with the label in front; other quotes are plain", () => {
  const note = md("> **Note** The first run downloads about 3 GB.\n> It takes a few minutes.");
  assert.match(note, /<div class="doc-callout note" role="note"><p><b class="doc-callout-label">Note<\/b> The first run downloads about 3 GB\. It takes a few minutes\.<\/p><\/div>/);
  const warn = md("> **Warning** Do not close the window.");
  assert.match(warn, /<div class="doc-callout warning" role="note"><p><b class="doc-callout-label">Warning<\/b> Do not close the window\.<\/p><\/div>/);
  const rich = md("> **Note** Steps:\n>\n> - one\n> - two\n>\n> Then run `mzdesk`.");
  assert.match(rich, /<p><b class="doc-callout-label">Note<\/b> Steps:<\/p><ul><li>one<\/li><li>two<\/li><\/ul><p>Then run <code>mzdesk<\/code>\.<\/p>/);
  assert.match(md("> **Note**\n> - a"), /<p><b class="doc-callout-label">Note<\/b><\/p><ul><li>a<\/li><\/ul>/);
  assert.match(md("> just a quote"), /^<div class="doc-body"><blockquote><p>just a quote<\/p><\/blockquote><\/div>$/);
  assert.doesNotMatch(md("> # not a page heading"), /<h1 id/, "headings inside a quote get no anchor");
});

test("--- is a horizontal rule", () => {
  assert.match(md("above\n\n---\n\nbelow"), /<p>above<\/p><hr\/><p>below<\/p>/);
  assert.match(md("***"), /<hr\/>/);
});

test("*italic*, **bold**, `code`, and the old <sub> tags work together; code and stars stay literal inside code", () => {
  assert.match(md("An *italic* word, a **bold `code`** one, H<sub>2</sub>O."), /<p>An <i>italic<\/i> word, a <b>bold <code>code<\/code><\/b> one, H<sub>2<\/sub>O\.<\/p>/);
  assert.match(md("Use `a*b*c` and 2 * 3 * 4."), /<code>a\*b\*c<\/code> and 2 \* 3 \* 4\./);
  assert.doesNotMatch(md("<script>alert(1)</script> and <img src=x onerror=1>"), /<script|<img/, "other tags stay literal text");
});

test("links: app pages use the router, https opens a new tab safely, #anchor stays on the page, anything else is not a link", () => {
  const html = md("[Install](/docs/install#windows) · [Paper](https://example.org/p?a=1&b=2) · [Top](#top) · [Mail](mailto:a@b.org) · [Bad](javascript:alert(1)) · [Rel](other.md) · [Script](/install.sh)");
  assert.match(html, /<a href="\/docs\/install#windows">Install<\/a>/);
  assert.match(html, /<a href="https:\/\/example\.org\/p\?a=1&amp;b=2" target="_blank" rel="noopener noreferrer">Paper<\/a>/);
  assert.match(html, /<a href="#top">Top<\/a>/);
  assert.match(html, /<a href="mailto:a@b\.org">Mail<\/a>/);
  assert.doesNotMatch(html, /javascript:/);
  assert.match(html, /<span>Bad<\/span> · <span>Rel<\/span>/);
  assert.match(html, /<a href="\/install\.sh">Script<\/a>/, "a file of the server is a plain link, not a router page");
  assert.match(md("See [the **tour**](/docs/interface-tour)."), /<a href="\/docs\/interface-tour">the <b>tour<\/b><\/a>/);
});

// ---------------------------------------------------------------- figures
const SVG = '<svg viewBox="0 0 720 100" width="720" height="100" xmlns="http://www.w3.org/2000/svg"><title>T</title><desc>D</desc>'
  + '<defs><marker id="a" markerWidth="8" markerHeight="8"><path class="fg-head" d="M0 0L8 4L0 8z"/></marker></defs>'
  + '<rect class="fg-box" x="1" y="1" width="100" height="50"/><path class="fg-arrow" d="M0 0L9 9" marker-end="url(#a)"/><text class="fg-text" x="5" y="9">a &amp; b</text></svg>';

function fetchStub(files) {
  const asked = [];
  globalThis.fetch = async url => {
    asked.push(String(url));
    const body = files[String(url).split("/").pop()];
    return body === undefined ? { ok: false, status: 404, text: async () => "" } : { ok: true, status: 200, text: async () => body };
  };
  return asked;
}

test("a figure shows a skeleton until its SVG is there, then the SVG inline with the alt text as its label, and the caption", async () => {
  figures.clearFigures();
  const text = "![The flow from input to answer](fig:flow-a)\n*One input, one answer.*\n\nAfter.";
  const before = md(text);
  assert.match(before, /<figure class="doc-figure"><div class="doc-fig-skel" aria-hidden="true"><i><\/i><\/div><figcaption>One input, one answer\.<\/figcaption><\/figure><p>After\.<\/p>/);
  const asked = fetchStub({ "flow-a": SVG });
  assert.equal(await figures.loadFigure("flow-a"), figures.cachedFigure("flow-a"));
  await figures.loadFigure("flow-a");
  assert.deepEqual(asked, ["/api/docs-img/flow-a"], "one request, then the module-level cache");
  const html = md(text);
  assert.match(html, /<div class="doc-fig-scroll" role="img" aria-label="The flow from input to answer" tabindex="0"><svg /);
  assert.match(html, /aria-hidden="true" focusable="false"/, "the drawing itself is hidden from screen readers; the wrapper carries the label");
  assert.match(html, /<rect class="fg-box" x="1" y="1" width="100" height="50"\/>/);
  assert.match(html, /id="fg-flow-a-a"/, "ids are prefixed with the figure name");
  assert.match(html, /marker-end="url\(#fg-flow-a-a\)"/);
  assert.match(html, /<text class="fg-text" x="5" y="9">a &amp; b<\/text>/);
  assert.doesNotMatch(html, /width="720"|height="100"|xmlns/, "the page sets the size");
});

test("a figure without a caption line has no <figcaption>; a missing file leaves the description", async () => {
  figures.clearFigures();
  fetchStub({});
  assert.equal(await figures.loadFigure("nothing-here"), null);
  const html = md("![A figure that is gone](fig:nothing-here)\n\nNot a caption line.");
  assert.match(html, /<figure class="doc-figure"><p class="doc-fig-missing">A figure that is gone<\/p><\/figure><p>Not a caption line\.<\/p>/);
  assert.doesNotMatch(html, /figcaption/);
});

test("a caption must be the very next line and a whole *italic* line", () => {
  figures.clearFigures();
  assert.doesNotMatch(md("![A](fig:x)\n\n*Later italic*"), /figcaption/);
  assert.doesNotMatch(md("![A](fig:x)\n**bold line**"), /figcaption/);
  assert.match(md("![A](fig:x)\n*With `code` and **bold**.*"), /<figcaption>With <code>code<\/code> and <b>bold<\/b>\.<\/figcaption>/);
});

// ---------------------------------------------------------------- the sanitizer
test("only text that starts with <svg is accepted", () => {
  for (const bad of ["", "<html><body><svg></svg></body></html>", '{"detail":"no such figure"}', "<!doctype html><svg/>", "plain text", "<svgx></svgx>"])
    assert.equal(figures.sanitizeSvg(bad), null, JSON.stringify(bad));
  assert.match(figures.sanitizeSvg('  \n<svg viewBox="0 0 10 10"><g/></svg>\n'), /^<svg viewBox="0 0 10 10" aria-hidden="true" focusable="false"><g\/><\/svg>$/);
});

test("scripts, foreignObject, event handlers, javascript: links, styles and other elements are removed", () => {
  const evil = '<svg viewBox="0 0 9 9" onload="alert(1)"><script>alert(1)</script><SCRIPT type="x">a<b</SCRIPT>'
    + '<foreignObject><div onclick="x()"><script>alert(2)</script>inside</div></foreignObject>'
    + '<style>*{display:none}</style><image href="http://evil/x.png"/><iframe src="//evil"></iframe><a href="javascript:alert(3)"><text>link</text></a>'
    + '<g onclick="x()" ONMOUSEOVER="y()" style="fill:red" class="fg-box"><use href="javascript:alert(4)"/><use xlink:href=" jAvaScript:alert(5)"/><use href="data:image/svg+xml;base64,AAA"/><use href="#ok"/>'
    + '<rect fill="url(http://evil/x)" stroke="url(#ok)" width="1"/><path d="M0 0" onerror=alert(6) /></g>'
    + '<text>1 &lt; 2</text><!-- <script>x</script> --></svg><script>after()</script>';
  const out = figures.sanitizeSvg(evil, "p-");
  assert.ok(out && out.startsWith("<svg") && out.endsWith("</svg>"));
  assert.doesNotMatch(out, /script|foreignobject|onload|onclick|onmouseover|onerror|javascript|vbscript|data:|<style|style=|<image|<iframe|<a[ >]|evil|inside|alert|after\(|<!--/i);
  assert.match(out, /<g class="fg-box">/);
  assert.match(out, /<use href="#p-ok"\/>/, "an in-page reference stays, with the prefix");
  assert.match(out, /stroke="url\(#p-ok\)"/);
  assert.match(out, /<text>1 &lt; 2<\/text>/);
  assert.equal(count(out, /<rect /g), 1);
});

test("the output is balanced even when the input is not, and nothing after the root is kept", () => {
  assert.equal(figures.sanitizeSvg('<svg><g><text>a</svg>'), '<svg aria-hidden="true" focusable="false"><g><text>a</text></g></svg>');
  assert.equal(figures.sanitizeSvg('<svg><g></g></g></svg><g/>'), '<svg aria-hidden="true" focusable="false"><g></g></svg>');
  assert.equal(figures.sanitizeSvg('<svg><text>a<b</text></svg>'), '<svg aria-hidden="true" focusable="false"><text>a&lt;b</text></svg>');
  const cdata = figures.sanitizeSvg('<svg><text><![CDATA[a <b> & c]]></text></svg>');
  assert.match(cdata, /<text>a &lt;b> &amp; c<\/text>/);
});

test("every figure that is published uses only what the sanitizer keeps", async () => {
  const { readdirSync, readFileSync, existsSync } = await import("node:fs");
  const dir = new URL("../../mzdesk/docs/img/", import.meta.url);
  if (!existsSync(dir)) return;
  for (const f of readdirSync(dir).filter(n => n.endsWith(".svg"))) {
    const src = readFileSync(new URL(f, dir), "utf8");
    const out = figures.sanitizeSvg(src, "t-");
    assert.ok(out, `${f} is an SVG`);
    const tags = s => [...s.matchAll(/<([a-zA-Z][\w:-]*)/g)].map(m => m[1].toLowerCase()).filter(t => t !== "svg");
    assert.deepEqual(tags(out), tags(src.replace(/<!--[\s\S]*?-->/g, "")), `${f}: nothing was dropped`);
    assert.match(src, /<title>[^<]+<\/title>/, `${f} has a title`);
    assert.match(src, /<desc>[^<]+<\/desc>/, `${f} has a description`);
    assert.match(src, /viewBox="0 0 720 \d+"/, `${f} is 720 wide`);
    assert.doesNotMatch(src, /<style|<script|<image|font-family|#[0-9a-fA-F]{3,6}\b|rgb\(/, `${f} uses only the fg-* classes`);
  }
});
