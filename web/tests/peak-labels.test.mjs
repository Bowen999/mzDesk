import assert from "node:assert/strict";
import { after, test } from "node:test";
import { startSsr } from "./ssr-helper.mjs";

const ssr = await startSsr("peak-labels");
const { labelEligible, MIN_LABEL_FRACTION } = await ssr.load("/src/components/peakLabels.ts");
after(() => ssr.close());

const peak = (mz, rel, extra = {}) => ({ mz, rel, on: true, kind: "Peak", rows: [], ...extra });
const all = peaks => peaks.map((_, i) => i);
const eligible = (peaks, precursors = [], from = all(peaks)) => [...labelEligible(peaks, from, precursors)].sort((a, b) => a - b);

test("the threshold is 5 % of the base peak", () => {
  assert.equal(MIN_LABEL_FRACTION, 0.05);
});

test("a single spectrum labels peaks of at least 5 % of the base peak; tiny near-zero peaks stay unlabelled", () => {
  const peaks = [peak(56.05, 2), peak(69.04, 18), peak(97.0, 1), peak(110.07, 101), peak(138.06, 1000), peak(180.06, 49.9), peak(181.07, 50)];
  // base 1000: 5 % = 50; 50 passes (>=), 49.9, 18, 2 and 1 do not
  assert.deepEqual(eligible(peaks), [3, 4, 6]);
});

test("each half of a mirror plot has its own base peak", () => {
  const peaks = [peak(100, 100), peak(110, 4), peak(120, 6),         // upper: 4 is under 5 % of 100, 6 is not
    peak(100, -20), peak(110, -0.9), peak(120, -1.1)];                // lower: base 20, so 1.0 is the line: 0.9 is out, 1.1 is in
  assert.deepEqual(eligible(peaks), [0, 2, 3, 5]);
});

test("the precursor is labelled whatever its height: by its card kind or by its dotted line", () => {
  const peaks = [peak(50, 100), peak(120.5, 0.4, { kind: "Near precursor m/z" }), peak(195.0878, 0.2), peak(130, 0.3)];
  assert.deepEqual(eligible(peaks, [195.0877]), [0, 1, 2]);
  assert.deepEqual(eligible(peaks), [0, 1]);                         // no precursor line: only the card kind marks one
  assert.deepEqual(eligible([peak(50, 100), peak(195.12, 0.2)], [195.0877]), [0], "0.032 Da away is a different peak (window 0.02 Da)");
});

test("coloured fragment peaks are never dropped", () => {
  const peaks = [peak(60, 100), peak(70, 1, { color: "#e4572e" }), peak(80, 1)];
  assert.deepEqual(eligible(peaks), [0, 1]);
});

test("after a zoom the base peak is the tallest peak in view", () => {
  const peaks = [peak(60, 100), peak(150, 4), peak(152, 3), peak(154, 0.1)];
  assert.deepEqual(eligible(peaks), [0]);                            // full view: 4 is under 5 % of 100
  assert.deepEqual(eligible(peaks, [], [1, 2, 3]), [1, 2]);          // zoomed to 145-160: 4 is the base peak there, 3 and 4 both pass, 0.1 does not
});

test("a half without a positive base peak, or an empty view, keeps every peak eligible and does not throw", () => {
  assert.deepEqual(eligible([peak(60, 0), peak(70, 0)]), [0, 1]);
  assert.deepEqual(eligible([]), []);
  assert.deepEqual(eligible([peak(60, 100)], [], []), []);
});
