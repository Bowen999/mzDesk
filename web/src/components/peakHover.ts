import type { PeakCard } from "../types";

type Axis = { _offset: number; _length: number; l2p: (value: number) => number };
export type PeakHit = { index: number; x: number; y0: number; y1: number; tipVisible: boolean };

/** The common stick traces share a baseline; MetFrag starts slightly below zero. */
export function stickBaseline(traces: any[]): number {
  const trace = traces.find(t => t.mode === "lines" && t.hoverinfo === "skip" && Array.isArray(t.y) && Number.isFinite(t.y[0]));
  return trace?.y[0] ?? 0;
}

/** Find the nearest visible stick, including its full length and a small target around its tip. */
export function hitPeak(peaks: PeakCard[], xa: Axis, ya: Axis, x: number, y: number, baseline = 0, radius = 7): PeakHit | null {
  const px = x - xa._offset, py = y - ya._offset;
  if (px < 0 || px > xa._length || py < 0 || py > ya._length) return null;
  const base = ya.l2p(baseline);
  let best: PeakHit | null = null, bestDistance = radius * radius;
  for (let index = 0; index < peaks.length; index++) {
    const peak = peaks[index], sx = xa.l2p(peak.mz), tip = ya.l2p(peak.rel);
    if (!Number.isFinite(sx) || !Number.isFinite(tip) || sx < 0 || sx > xa._length) continue;
    if (Math.max(base, tip) < 0 || Math.min(base, tip) > ya._length) continue;
    const y0 = Math.max(0, Math.min(ya._length, base)), y1 = Math.max(0, Math.min(ya._length, tip));
    const dx = px - sx, dy = py - Math.max(Math.min(y0, y1), Math.min(Math.max(y0, y1), py));
    const distance = dx * dx + dy * dy;
    if (distance <= bestDistance && (!best || distance < bestDistance)) {
      bestDistance = distance;
      best = { index, x: xa._offset + sx, y0: ya._offset + y0, y1: ya._offset + y1, tipVisible: tip >= 0 && tip <= ya._length };
    }
  }
  return best;
}

/** Place the measured card beside the pointer, flipping and clamping it at viewport edges. */
export function hoverCardPosition(x: number, y: number, width: number, height: number, viewportWidth: number, viewportHeight: number) {
  const margin = 6, gap = 16;
  const right = x + gap, left = x - gap - width;
  const preferred = right + width <= viewportWidth - margin ? right : left >= margin ? left : right;
  return {
    left: Math.max(margin, Math.min(preferred, viewportWidth - width - margin)),
    top: Math.max(margin, Math.min(y - height / 2, viewportHeight - height - margin)),
  };
}
