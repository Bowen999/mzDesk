import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { after, test } from "node:test";
import { startSsr } from "./ssr-helper.mjs";

const ssr = await startSsr("peak-count");
const { countPeaks } = await ssr.load("/src/peakCount.ts");
after(() => ssr.close());

// The cases the server's parser is tested on (tests/test_spectrum.py: test_peak_table_cases): the page must count the peaks the run will use.
const cases = JSON.parse(readFileSync(new URL("../../tests/fixtures/peak_table_cases.json", import.meta.url), "utf8"));

test("the peak count of every table in the shared cases is the number of peaks the server reads", () => {
  assert.ok(cases.length >= 30);
  for (const { name, text, peaks } of cases) assert.equal(countPeaks(text), peaks.length, name);
});

test("a semicolon csv is one peak per line, not one per number", () => {
  assert.equal(countPeaks("195.0877;100\n138.0662;60"), 2);
  assert.equal(countPeaks("m/z;intensity\n195.0877;100\n138.0662;60\n110.0713;20"), 3);
  assert.equal(countPeaks("41;100\n55;85\n69;60"), 3);
});

test("decimal commas next to a tab, space or semicolon are decimal points: peaks that share the integer part stay separate", () => {
  assert.equal(countPeaks("138,0662\t100\n138,0877\t50\n138,5\t3"), 3);                  // read as m/z 138 each, this was 1 peak
  assert.equal(countPeaks("138,0662;100\n195,0877;32"), 2);
});

test("a comma alone still separates the columns, and the semicolon forms that were peaks stay peaks", () => {
  assert.equal(countPeaks("41,100\n55,85\n69,60"), 3);
  assert.equal(countPeaks("41,100;55,85;69,60"), 3);
  assert.equal(countPeaks("138.0662 100;195.0877 32"), 2);
  assert.equal(countPeaks("100.1;200.2;300.3"), 3);
  assert.equal(countPeaks("195.0877,100"), 1);
});

test("JSON, XML and array texts are not counted, as before", () => {
  assert.equal(countPeaks("[[100.1, 5], [200.2, 100]]"), 0);
  assert.equal(countPeaks("mz = [100.1, 200.2]\nintensities = [5, 100]"), 0);
  assert.equal(countPeaks(""), 0);
});
