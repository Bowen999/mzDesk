import assert from "node:assert/strict";
import { mkdtempSync, rmSync } from "node:fs";
import { tmpdir } from "node:os";
import { dirname, join } from "node:path";
import { after, test } from "node:test";
import { fileURLToPath } from "node:url";
import { createServer } from "vite";

const root = dirname(dirname(fileURLToPath(import.meta.url)));
const cache = mkdtempSync(join(tmpdir(), "mzdesk-peak-hover-test-"));
const server = await createServer({ configFile: false, root, cacheDir: cache, server: { middlewareMode: true, hmr: false }, appType: "custom" });
const { hitPeak, hoverCardPosition, stickBaseline } = await server.ssrLoadModule("/src/components/peakHover.ts");
after(async () => { await server.close(); rmSync(cache, { recursive: true, force: true }); });

const peak = (mz, rel) => ({ mz, rel, on: true, kind: "Peak", rows: [] });
const xa = { _offset: 50, _length: 400, l2p: x => x };
const ya = { _offset: 20, _length: 200, l2p: y => 100 - y };

test("hits the whole stick, but not empty space away from it", () => {
  const peaks = [peak(100, 80), peak(140, 10)];
  assert.equal(hitPeak(peaks, xa, ya, 153, 80)?.index, 0);
  assert.equal(hitPeak(peaks, xa, ya, 190, 105)?.index, 1);
  assert.equal(hitPeak(peaks, xa, ya, 150, 25), null); // above the tip
  assert.equal(hitPeak(peaks, xa, ya, 160, 80), null); // beside the stick
  assert.equal(hitPeak(peaks, xa, ya, 45, 100), null); // outside the plot
});

test("mirror halves and closely spaced peaks choose the nearest visible stick", () => {
  const peaks = [peak(100, 80), peak(100, -80), peak(104, 40)];
  assert.equal(hitPeak(peaks, xa, ya, 150, 60)?.index, 0);
  assert.equal(hitPeak(peaks, xa, ya, 150, 180)?.index, 1);
  assert.equal(hitPeak(peaks, xa, ya, 154, 90)?.index, 2);
  assert.equal(hitPeak(peaks, xa, ya, 152, 60)?.index, 0); // third peak does not reach here
});

test("tiny and zero peaks remain selectable at the baseline", () => {
  const peaks = [peak(100, 0.01), peak(140, 0)];
  assert.equal(hitPeak(peaks, xa, ya, 150, 120)?.index, 0);
  assert.deepEqual(hitPeak(peaks, xa, ya, 190, 120), { index: 1, x: 190, y0: 120, y1: 120, tipVisible: true });
});

test("sticks can start below zero, as in MetFrag figures", () => {
  const traces = [{ mode: "markers", y: [50] }, { mode: "lines", hoverinfo: "skip", y: [-1.6, 80, null] }];
  assert.equal(stickBaseline(traces), -1.6);
  assert.equal(stickBaseline([]), 0);
  assert.equal(hitPeak([peak(100, 80)], xa, ya, 150, 121, stickBaseline(traces)).y0, 121.6);
});

test("zoomed axes clip highlights and ignore sticks outside the view", () => {
  const zoomY = { _offset: 20, _length: 40, l2p: y => 80 - y };
  const peaks = [peak(100, 80), peak(140, 10), peak(500, 80)];
  assert.deepEqual(hitPeak(peaks, xa, zoomY, 150, 40), { index: 0, x: 150, y0: 60, y1: 20, tipVisible: true });
  assert.equal(hitPeak(peaks, xa, zoomY, 190, 60), null);
  const clippedTip = hitPeak([peak(100, 100)], xa, zoomY, 150, 40);
  assert.equal(clippedTip.tipVisible, false);
  assert.equal(clippedTip.y1, 20);
});

test("the measured tooltip flips and stays within the viewport", () => {
  assert.deepEqual(hoverCardPosition(100, 100, 250, 160, 800, 600), { left: 116, top: 20 });
  assert.deepEqual(hoverCardPosition(790, 590, 250, 160, 800, 600), { left: 524, top: 434 });
  assert.deepEqual(hoverCardPosition(150, 10, 288, 300, 300, 320), { left: 6, top: 6 });
});
