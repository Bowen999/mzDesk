import assert from "node:assert/strict";
import { after, test } from "node:test";
import { React, startSsr } from "./ssr-helper.mjs";

// A column whose every cell is a number or a short one-line label never wraps (td.nowrap): the table scrolls sideways anyway, so wrapping "Sample type" only doubles the row.
const ssr = await startSsr("nowrap-columns");
const { RowsTable, nowrapColumns, NOWRAP_MAX } = await ssr.load("/src/components/Blocks.tsx");
const { OutputBlock } = await ssr.load("/src/components/Outputs.tsx");
after(() => ssr.close());

const h = React.createElement;
const screen = (component, props) => ssr.render(component, props);
const long = "x".repeat(NOWRAP_MAX + 1);

test("nowrapColumns: short single-line text and numbers do not wrap; long or multi-line text does", () => {
  const rows = [
    ["Sample type", "blood", "12", long, "a\nb", { text: "short" }, { badge: "hit" }],
    ["Sex", "—", "3.5", "tiny", "c", { text: "x".repeat(NOWRAP_MAX) }, { badge: "y".repeat(NOWRAP_MAX + 1) }],
    { cells: ["Age", null, 7, "", "d", { fragment: { usi: "mzspec:A:b:scan:1" } }, undefined] },
  ];
  assert.deepEqual(nowrapColumns(7, rows), [true, true, true, false, false, true, false]);
});

test("nowrapColumns: the limit is inclusive, numbers always qualify, and an empty column is harmless", () => {
  assert.equal(NOWRAP_MAX, 24);
  assert.deepEqual(nowrapColumns(1, [["x".repeat(24)]]), [true]);
  assert.deepEqual(nowrapColumns(1, [["x".repeat(25)]]), [false]);
  assert.deepEqual(nowrapColumns(2, [["1,234,567.891234567890123456789", "—"]]), [true, true]);       // a number is nowrap through `numeric`
  assert.deepEqual(nowrapColumns(1, []), [true]);
});

const columns = ["Property", "Value", "Hit files", "Note"];
const rows = [
  ["Sample type", "blood", 12, "a long explanation of this property that goes past the limit"],
  ["Disease", "healthy", 3, "short"],
];

test("RowsTable puts nowrap on the td of every short or numeric column, not on a long one, and never on a th", () => {
  const html = screen(RowsTable, { columns, rows, empty: "none" });
  const tds = [...html.matchAll(/<td class="([^"]*)"/g)].map(m => m[1]);
  assert.deepEqual(tds.slice(0, 4), ["nowrap", "nowrap", "num nowrap", "long"]);
  assert.doesNotMatch(html, /<th[^>]*nowrap/);
  assert.match(html, /<th class="num">Hit files/);
});

test("a nowrap array passed in wins over the rows shown, so a page of short rows does not change a column", () => {
  const html = screen(RowsTable, { columns, rows: rows.slice(1), empty: "none", numeric: [false, false, true, false], nowrap: [true, true, true, false] });
  const tds = [...html.matchAll(/<td class="([^"]*)"/g)].map(m => m[1]);
  assert.deepEqual(tds, ["nowrap", "nowrap", "num nowrap"]);
  assert.match(html, /<td>short<\/td>/, "the Note column keeps wrapping: the passed array, not this page, decides");
});

test("a paged table decides over all its rows: page 1 with only short notes keeps the long Note column wrapping", () => {
  const many = [["a", "b", 1, "short"], ["c", "d", 2, "short"], ["e", "f", 3, "short"], ["g", "h", 4, "this note is much longer than the limit allows"]];
  const html = screen(OutputBlock, { block: { type: "table", columns, page_size: 2 }, data: { columns, rows: many }, empty: "none" });
  assert.doesNotMatch(html, /<td class="nowrap">short/, "the Note column is not nowrap, although its first page is all short");
  assert.match(html, /<td class="nowrap">a<\/td>/);
});
