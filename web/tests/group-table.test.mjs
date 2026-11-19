import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { after, test } from "node:test";
import { React, startSsr } from "./ssr-helper.mjs";

// The grouptable block with `page_size` (pages of groups, a filter, "Expand all"), the dataset link on a group, and the report, which shows every group.
// The screen test renders the first state (nothing typed, nothing open); what the filter matches is checked on groupSearchText, the function it filters with.
const ssr = await startSsr("group-table");
const { GroupTable, groupSearchText, groupLinkName, usedColumns, cellKind, numericColumns, RowsTable, trimEmptyColumns } = await ssr.load("/src/components/Blocks.tsx");
const { OutputBlock } = await ssr.load("/src/components/Outputs.tsx");
const { ReportModeContext } = await ssr.load("/src/report/mode.ts");
after(() => ssr.close());

const h = React.createElement;
const count = (html, re) => (html.match(re) ?? []).length;
const screen = (component, props) => ssr.render(component, props);
const report = (component, props) => ssr.render(() => h(ReportModeContext.Provider, { value: { maxRows: 500, track() {} } }, h(component, props)));

const columns = ["USI", "Title"];
const groups = Array.from({ length: 7 }, (_, i) => ({
  title: `MSV00000${i + 1}`, meta: `${i + 1} hits`,
  rows: [[`mzspec:MSV00000${i + 1}:file${i + 1}:scan:1`, i === 3 ? { text: "caffeine standard" } : `sample ${i + 1}`]],
}));

test("a group table without page_size draws as before: every group, no filter, no Expand all, no pages", () => {
  const html = screen(GroupTable, { columns, groups });
  assert.equal(count(html, /class="gt-group"/g), 7);
  assert.doesNotMatch(html, /result-table|<input|Expand all|gt-link|<nav/);
  assert.match(html, /<div class="gt">/);
});

test("with page_size the first page of groups is drawn, with the filter, 'Expand all' and the page controls", () => {
  const html = screen(GroupTable, { columns, groups, pageSize: 3, name: "Passing hits by dataset" });
  assert.equal(count(html, /class="gt-group"/g), 3);
  assert.ok(html.includes("MSV000003") && !html.includes("MSV000004"));
  assert.match(html, /<input type="search"[^>]*aria-label="Filter Passing hits by dataset groups"/);
  assert.match(html, />Expand all<\/button>/);
  assert.match(html, /1–3 of 7 groups/);
  assert.match(html, /<nav class="result-table-pages" aria-label="Passing hits by dataset pages"><span aria-live="polite">Page 1 of 3<\/span>/);
  assert.match(html, /<button type="button" class="page-button" disabled="">Previous<\/button><button type="button" class="page-button">Next<\/button>/);
});

test("page_size at least the number of groups: no filter and no pages, only 'Expand all' (more than one group)", () => {
  const html = screen(GroupTable, { columns, groups, pageSize: 7 });
  assert.equal(count(html, /class="gt-group"/g), 7);
  assert.doesNotMatch(html, /<input|<nav|result-table-count/);
  assert.match(html, />Expand all<\/button>/);
  assert.doesNotMatch(screen(GroupTable, { columns, groups: groups.slice(0, 1), pageSize: 7 }), /Expand all/, "one group has nothing to expand all of");
});

test("the filter matches a group by its title, its meta or the text of its rows (not by hidden markup)", () => {
  const match = (g, q) => groupSearchText(g).includes(q.toLowerCase());
  assert.ok(match(groups[0], "msv000001"), "title");
  assert.ok(match(groups[1], "2 hits"), "meta");
  assert.ok(match(groups[3], "caffeine"), "a cell's text");
  assert.ok(match(groups[3], "file4"), "a USI in a row");
  assert.ok(!match(groups[2], "caffeine"));
  assert.ok(!match({ title: "<b>Set</b> one", rows: [] }, "<b>"), "tags of a title are not searched");
  assert.ok(match({ title: "<b>Set</b> one", rows: [] }, "set one"));
  assert.ok(match({ title: "G", rows: [{ cells: [{ badge: "high", tone: "good" }] }] }, "high"), "badge text of an object row");
  assert.ok(!match({ title: "G", rows: [[{ img: "data:image/png;base64,AA" }]] }, "data:"), "an image has no text");
  assert.ok(match({ title: "G" }, "g"), "a group without rows is fine");
});

test("url: a small external link after the title, opening a new tab; an unsafe url is dropped", () => {
  const html = screen(GroupTable, { columns, groups: [{ ...groups[0], url: "https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000001" }, { ...groups[1], url: "javascript:alert(1)" }, groups[2]] });
  assert.equal(count(html, /class="gt-link"/g), 1);
  assert.match(html, /<\/button><a class="gt-link" href="https:\/\/massive\.ucsd\.edu\/ProteoSAFe\/dataset\.jsp\?accession=MSV000001" target="_blank" rel="noopener noreferrer" title="Open dataset page: MSV000001" aria-label="Open dataset page: MSV000001">/);
  assert.doesNotMatch(html, /javascript:/);
  // the link is outside the toggle button, so the button stays the keyboard target of the group
  assert.doesNotMatch(html, /<button[^>]*gt-toggle[^>]*>(?:(?!<\/button>).)*gt-link/);
});

test("report mode shows every group as <details>, without paging, filter or Expand all, and the url as a link", () => {
  const html = report(GroupTable, { columns, pageSize: 3, name: "Hits", groups: groups.map((g, i) => (i === 0 ? { ...g, url: "https://example.org/d" } : g)) });
  assert.equal(count(html, /<details class="gt-group">/g), 7);
  assert.ok(html.includes("MSV000007") && html.includes("file7"), "the last group and its rows are in the file");
  assert.doesNotMatch(html, /<input|<nav|<button|Expand all|result-table/);
  assert.match(html, /<\/span><a class="gt-link" href="https:\/\/example\.org\/d" target="_blank" rel="noopener noreferrer"[^>]*><svg/);
  assert.equal(count(html, /class="gt-link"/g), 1);
});

test("OutputBlock passes page_size to a grouptable and nothing changes without it", () => {
  const data = { groups, count: 7 };
  const paged = screen(OutputBlock, { block: { id: "g", type: "grouptable", title: "By dataset", columns, page_size: 5 }, data, status: "pending" });
  assert.equal(count(paged, /class="gt-group"/g), 5);
  assert.match(paged, /Filter By dataset groups/);
  const plain = screen(OutputBlock, { block: { id: "g", type: "grouptable", title: "By dataset", columns }, data, status: "pending" });
  assert.equal(count(plain, /class="gt-group"/g), 7);
  assert.doesNotMatch(plain, /<input|Expand all|<nav/);
});

test("a dataset link is named after its group: the accession before ' · ', else the plain title", () => {
  assert.equal(groupLinkName("MSV000087728 · Plant extracts <i>LC-MS</i>"), "MSV000087728");
  assert.equal(groupLinkName("<b>MSV000001</b>"), "MSV000001");
  assert.equal(groupLinkName("GNPS2 Explorer"), "GNPS2 Explorer");
  const html = screen(GroupTable, { columns, groups: [
    { title: "MSV000087728 · Plant extracts", url: "https://example.org/a", rows: [["x", "y"]] },
    { title: "MSV000000002 · Other", url: "https://example.org/b", rows: [["x", "y"]] }] });
  assert.match(html, /title="Open dataset page: MSV000087728" aria-label="Open dataset page: MSV000087728"/);
  assert.match(html, /aria-label="Open dataset page: MSV000000002"/);
  assert.equal(count(html, /aria-label="Open dataset page"/g), 0, "no two links share one label");
});

test("a group's meta is always on a line of its own under the title: one class, no long/short switch, never nowrap", () => {
  const long = "Pan-ReDU 12 of 15 hit files · plantMASST 7 of 15 · GNPS2 Explorer covers 5 of 6 datasets in total";
  const gs = [{ title: "Organism", meta: long, rows: [["a", "b"]] }, { title: "Caffeine", meta: "3 hits", rows: [["a", "b"]] }];
  const html = screen(GroupTable, { columns, groups: gs });
  assert.match(html, new RegExp(`<span class="gt-meta">${long}</span>`));
  assert.match(html, /<span class="gt-meta">3 hits<\/span>/);
  assert.doesNotMatch(html, /gt-meta-long/);
  const inReport = report(GroupTable, { columns, groups: gs });
  assert.doesNotMatch(inReport, /gt-meta-long/);
  assert.match(inReport, /<span class="gt-meta">3 hits<\/span>/);
  const css = readFileSync(new URL("../src/outputs-v2.css", import.meta.url), "utf8");
  const rule = name => css.match(new RegExp(`^\\${name}\\s*\\{([^}]*)\\}`, "m"))?.[1] ?? "";
  assert.doesNotMatch(rule(".gt-meta"), /nowrap/, "a meta is never nowrap (it would be clipped by .gt)");
  assert.match(rule(".gt-meta"), /flex:\s*1 1 100%/, "its own line");
  assert.match(rule(".gt-meta"), /margin-left:\s*0/, "left-aligned, not pushed right");
  assert.match(rule(".gt-meta"), /padding-left:\s*22px/, "past the caret");
  assert.doesNotMatch(css, /gt-meta-long|\.gt-meta[^{]*\{[^}]*white-space:\s*nowrap/);
  assert.doesNotMatch(css.slice(css.indexOf("@media (max-width: 700px)")), /\.gt-meta/, "no width-specific meta layout");
});

test("the meta's colour is --muted (4.5:1 on the page in both themes), not --faint", () => {
  const theme = readFileSync(new URL("../src/theme.css", import.meta.url), "utf8");
  const meta = theme.match(/^\.gt-meta\s*\{([^}]*)\}/m)?.[1] ?? "";
  assert.match(meta, /color:\s*var\(--muted\)/);
  assert.doesNotMatch(meta, /--faint/);
});

test("on a phone the group header toggle, the dataset link and the page buttons are at least 44px; focus stays visible", () => {
  const css = readFileSync(new URL("../src/outputs-v2.css", import.meta.url), "utf8");
  const phone = css.slice(css.indexOf("@media (max-width: 700px)"));
  const body = sel => phone.match(new RegExp(`^\\s*\\${sel}\\s*\\{([^}]*)\\}`, "m"))?.[1] ?? "";
  for (const sel of [".gt-toggle", ".gt-link", ".page-button"]) assert.match(body(sel), /min-height:\s*44px/, `${sel} is 44px tall`);
  for (const sel of [".gt-link", ".page-button"]) assert.match(body(sel), /min-width:\s*44px/, `${sel} is 44px wide`);
  assert.match(css, /\.gt-link:focus-visible/);
  assert.match(css, /\.gt-toggle:focus-visible/);
  assert.match(readFileSync(new URL("../src/theme.css", import.meta.url), "utf8"), /:where\(a, button[^)]*\):focus-visible/, "page-button / expand get the global focus ring");
});

// ---------------------------------------------------------------- numeric columns
test("a numeric cell is a plain number; IDs, adducts, ranges and text are not; blanks are empty", () => {
  for (const n of ["20", "0.800", "-1.5", "+2", "1,234", "12.5%", "40.54", ">1000", "<0.01", "≥2", "≤ 3", "1e-3", "0"]) assert.equal(cellKind(n), "number", n);
  for (const n of [{ text: "0.800" }, { text: ">1000", url: "https://example.org" }, { text: "12.5%", help: "h" }, 5]) assert.equal(cellKind(n), "number", JSON.stringify(n));
  for (const o of ["—", "", "  ", null, undefined, { text: "—" }, { text: "" }]) assert.equal(cellKind(o), "empty", JSON.stringify(o));
  for (const o of ["M+H", "[M+Na]+", "MSV000087728", "CCMSLIB00005435812", "3f2b8c1e-9d4a-4b6e-8f10-2a7c5d9e1b34", "100-200", "100.0–200.0", "caffeine", "1,23", "1.2.3", "12 files",
    "5\n6", "-", ".", "%", { badge: "7", tone: "good" }, { img: "data:image/png;base64,AA" }, { fragment: { smiles: "C", adduct: "M+H" } }, { text: "M+H" }]) assert.equal(cellKind(o), "other", JSON.stringify(o));
});

test("a column is numeric when every non-empty cell is: '—' and blanks do not break it, one text cell does, an all-empty column is not", () => {
  const rows = [
    ["a", "20", "M+H", "3f2b8c1e-9d4a-4b6e-8f10-2a7c5d9e1b34", "1", "—"],
    { cells: ["b", { text: "0.800" }, "M+Na", "3f2b8c1e-9d4a-4b6e-8f10-2a7c5d9e1b35", "—", ""], help: "h" },
    ["c", ">1000", "M+H", "3f2b8c1e-9d4a-4b6e-8f10-2a7c5d9e1b36", "x", null],
    ["d", "—", "M+H", "3f2b8c1e-9d4a-4b6e-8f10-2a7c5d9e1b37", "2", "—"],
  ];
  assert.deepEqual(numericColumns(6, rows), [false, true, false, false, false, false]);
  assert.deepEqual(numericColumns(2, []), [false, false], "no rows, nothing to align");
  assert.deepEqual(numericColumns(2, [["1,234", "12.5%"], ["-1.5", "<0.01"]]), [true, true]);
});

test("a table puts `num` on the th and every td of a numeric column, and only there", () => {
  const html = screen(RowsTable, { columns: ["Compound", "Hits", "Adduct", "Enrichment"], empty: "none", rows: [["Caffeine", "20", "M+H", ">1000"], ["Theobromine", "—", "M+H", "2.50"], ["Xanthine", "3", "M+Na", { text: "0.80" }]] });
  assert.match(html, /<th>Compound<\/th><th class="num">Hits<\/th><th>Adduct<\/th><th class="num">Enrichment<\/th>/);
  assert.equal(count(html, /<td class="num nowrap">/g), 6);
  assert.match(html, /<tr><td class="nowrap">Caffeine<\/td><td class="num nowrap">20<\/td><td class="nowrap">M\+H<\/td><td class="num nowrap">&gt;1000<\/td><\/tr>/);
  assert.match(html, /<td class="nowrap">Xanthine<\/td><td class="num nowrap">3<\/td>/);
  const css = readFileSync(new URL("../src/theme.css", import.meta.url), "utf8");
  const rule = css.match(/^\.out-table th\.num, \.out-table td\.num\s*\{([^}]*)\}/m)?.[1] ?? "";
  assert.match(rule, /text-align:\s*right/);
  assert.match(rule, /tabular-nums/);
});

test("a group table decides the alignment over every row of every group, so the columns line up (screen state and report)", () => {
  const cols = ["Species", "Hits", "Enrichment"];
  const gs = [
    { title: "G1", rows: [["Homo sapiens", "20", "—"], ["Mus musculus", "—", ">1000"]] },       // G1 alone: Hits and Enrichment are numbers
    { title: "G2", rows: [["Danio rerio", "n/a", "2.50"]] },                                    // "n/a" makes Hits a text column for every group
  ];
  assert.deepEqual(numericColumns(3, gs.flatMap(g => g.rows)), [false, false, true]);
  const all = report(GroupTable, { columns: cols, groups: gs });         // the report holds every group's rows
  assert.equal(count(all, /<th class="num">Enrichment<\/th>/g), 2, "both groups' tables right-align Enrichment");
  assert.equal(count(all, /<th>Hits<\/th>/g), 2, "Hits stays left in both, because one group has text in it");
  assert.equal(count(all, /<th class="num">Hits<\/th>/g), 0);
  assert.match(all, /<td class="num nowrap">&gt;1000<\/td>/);
  assert.match(all, /<td class="nowrap">20<\/td>/);
  assert.match(report(GroupTable, { columns: cols, groups: [gs[0]] }), /<th class="num">Hits<\/th>/, "G1 on its own would right-align Hits");
});

test("the report drops a column of 'View' buttons: the numeric flags follow the columns that stay", () => {
  const cols = ["Compound", "Fragments", "Score"];
  const rows = [["A", { fragment: { smiles: "C", adduct: "M+H" } }, "0.9"], ["B", { fragment: { smiles: "CC", adduct: "M+H" } }, "0.8"]];
  const html = report(GroupTable, { columns: cols, groups: [{ title: "G", rows }] });
  assert.doesNotMatch(html, /Fragments/);
  assert.match(html, /<th>Compound<\/th><th class="num">Score<\/th>/);
});

test("columns empty or '—' in every row of every group are hidden, decided over all groups", () => {
  const cols = ["USI", "Name", "Adduct", "Verification", "Notes"];
  const gs = [
    { title: "G1", rows: [["u1", "—", "", { text: "—" }, "—"], { cells: ["u2", "—", null, "", { text: "" }], help: "row help" }] },
    { title: "G2", rows: [["u3", "—", "", "—", "has a note"]] },
  ];
  assert.deepEqual(usedColumns(cols, gs), [0, 4], "Notes is used by one row of the second group");
  const t = trimEmptyColumns(cols, gs);
  assert.deepEqual(t.columns, ["USI", "Notes"]);
  assert.deepEqual(t.groups[0].rows[0], ["u1", "—"]);
  assert.deepEqual(t.groups[0].rows[1], { cells: ["u2", { text: "" }], help: "row help" }, "an object row keeps its other fields");
  assert.deepEqual(t.groups[1].rows[0], ["u3", "has a note"]);
  assert.equal(trimEmptyColumns(["A", "B"], [{ title: "G", rows: [["x", "y"]] }]).groups[0].rows[0].length, 2);
  assert.deepEqual(usedColumns(["A", "B"], [{ title: "G", rows: [["—", ""]] }]), [0, 1], "never every column");
  assert.deepEqual(usedColumns(["A", "B"], [{ title: "G", rows: [] }]), [0, 1], "no rows: nothing to judge");
  // a badge, an image or a link is something even when it has no text
  assert.deepEqual(usedColumns(["A", "B", "C"], [{ title: "G", rows: [["—", { badge: "ok" }, { img: "data:image/png;base64,AA" }]] }]), [1, 2]);
  assert.deepEqual(usedColumns(["A", "B"], [{ title: "G", rows: [["x", { text: "—", url: "https://example.org" }]] }]), [0, 1]);
  // the report draws every group open, so its tables are in the markup (the screen's groups start folded)
  const rendered = report(GroupTable, { columns: cols, groups: gs });
  assert.equal(count(rendered, /<table/g), 2);
  assert.equal(count(rendered, /<th[^>]*>USI/g), 2, "every group has the same two columns");
  assert.equal(count(rendered, /<th[^>]*>Notes/g), 2);
  assert.doesNotMatch(rendered, /<th[^>]*>(Name|Adduct|Verification)/);
  assert.match(rendered, /aria-label="Results table, 2 columns/);
});

test("OutputBlock reaches the column hiding too", () => {
  const cols = ["USI", "Name", "Title"];
  const data = { groups: [{ title: "G", rows: [["u", "—", "t"]] }] };
  const rep = report(OutputBlock, { block: { id: "g", type: "grouptable", title: "T", columns: cols }, data, status: "pending" });
  assert.match(rep, /<th[^>]*>USI/);
  assert.match(rep, /<th[^>]*>Title/);
  assert.doesNotMatch(rep, /<th[^>]*>Name/);
});

test("a block note shows under the title of a table and a grouptable, on screen and in the report; no note, no element", () => {
  const note = "Coverage: Pan-ReDU 8 of 15 hit files · plantMASST 7 of 15 · GNPS2 Explorer 5 of 6 datasets";
  const gblock = { id: "g", type: "grouptable", title: "Where it occurs", columns };
  const tblock = { id: "t", type: "table", title: "Hits", columns };
  const gdata = { groups, note }, tdata = { rows: [["a", "b"]], note };
  for (const [block, data] of [[gblock, gdata], [tblock, tdata]]) {
    for (const render of [screen, report]) {
      const html = render(OutputBlock, { block, data, status: "pending" });
      assert.match(html, new RegExp(`</h4><p class="out-note">${note}</p>`), `${block.type}`);
    }
  }
  const rich = screen(OutputBlock, { block: gblock, data: { groups, note: "<b>8</b> of 15" }, status: "pending" });
  assert.match(rich, /<p class="out-note"><b>8<\/b> of 15<\/p>/);
  for (const [block, data] of [[gblock, { groups }], [tblock, { rows: [["a", "b"]] }]]) {
    assert.doesNotMatch(screen(OutputBlock, { block, data, status: "pending" }), /out-note/);
  }
});

test("a number needs a digit in its mantissa: a bare exponent is text; 1e5, +2, 0.800, >1000, 12.5% are numbers", () => {
  for (const o of ["E1", "e5", "-e5", "+E3", "e", ".e5", "<e2"]) assert.equal(cellKind(o), "other", o);
  for (const n of ["1e5", "1E5", "+2", "0.800", ">1000", "12.5%", ".5", "-.5e-3", "1,234.5"]) assert.equal(cellKind(n), "number", n);
});

test("a plain table decides which columns are numbers over all its rows, not over the page, the fold or the filter", () => {
  const nums = (html, re) => count(html, re);
  // paged: the first page (2 rows) holds only numbers in "Score", a later row holds text -> not numeric on page 1 either
  const mixed = { columns: ["Name", "Score"], rows: [["a", "1"], ["b", "2"], ["c", "3"], ["d", "n/a"], ["e", "5"]] };
  const pagedBlock = { id: "t", type: "table", title: "T", columns: ["Name", "Score"], page_size: 2 };
  const first = screen(OutputBlock, { block: pagedBlock, data: mixed, status: "pending" });
  assert.match(first, /<td class="nowrap">a<\/td><td class="nowrap">1<\/td>/);
  assert.equal(nums(first, /class="num[ "]/g), 0, "a later text row keeps the column left-aligned on page 1");
  // the reverse: the first page shows only blanks in "Score", a later row is a number -> numeric on page 1
  const late = { columns: ["Name", "Score"], rows: [["a", "—"], ["b", ""], ["c", "3"], ["d", "4"]] };
  const lateFirst = screen(OutputBlock, { block: pagedBlock, data: late, status: "pending" });
  assert.match(lateFirst, /<th class="num">Score<\/th>/);
  assert.equal(nums(lateFirst, /<td class="num nowrap">/g), 2, "both visible Score cells of page 1 are right-aligned");
  // no page_size, 20 rows (folded to 8): the same rule
  const many = Array.from({ length: 20 }, (_, i) => ["r" + i, i < 8 ? "—" : String(i)]);
  const folded = screen(OutputBlock, { block: { id: "t", type: "table", title: "T", columns: ["Name", "Score"] }, data: { columns: ["Name", "Score"], rows: many }, status: "pending" });
  assert.match(folded, /Show all 20 rows/);
  assert.match(folded, /<th class="num">Score<\/th>/, "the 12 later numbers make the column numeric from the first fold");
  const text = many.map((r, i) => (i === 15 ? [r[0], "n/a"] : r));
  const foldedText = screen(OutputBlock, { block: { id: "t", type: "table", title: "T", columns: ["Name", "Score"] }, data: { columns: ["Name", "Score"], rows: text }, status: "pending" });
  assert.doesNotMatch(foldedText, /class="num[ "]/, "a text cell hidden in the fold counts");
  // the report shows every row the same way
  assert.doesNotMatch(report(OutputBlock, { block: pagedBlock, data: mixed, status: "pending" }), /class="num[ "]/);
  assert.match(report(OutputBlock, { block: pagedBlock, data: late, status: "pending" }), /<th class="num">Score<\/th>/);
});

test("the empty cell is --muted (readable), and Expand all is a page-button so it gets the phone touch target", () => {
  const theme = readFileSync(new URL("../src/theme.css", import.meta.url), "utf8");
  const empty = theme.match(/^\.out-table td\.empty\s*\{([^}]*)\}/m)?.[1] ?? "";
  assert.match(empty, /color:\s*var\(--muted\)/);
  assert.match(screen(GroupTable, { columns, groups, pageSize: 3 }), /class="page-button gt-expand"/);
});
