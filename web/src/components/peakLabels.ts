import type { PeakCard } from "../types";

/** A peak gets an m/z label only when it is at least this fraction of the base peak of its half (the upper or the lower half of a mirror plot). */
export const MIN_LABEL_FRACTION = 0.05;
/** A peak this close (Da) to a dotted precursor line is the precursor peak: it is labelled whatever its height. */
export const PRECURSOR_WINDOW = 0.02;

/** The peak at a precursor: its card says so ("Precursor ion", "Near precursor m/z"), or it lies within `window` Da of a dotted precursor line (mirror plots). */
export const isPrecursorPeak = (p: PeakCard, precursorXs: number[], window = PRECURSOR_WINDOW) =>
  /precursor/i.test(p.kind) || precursorXs.some(x => Math.abs(p.mz - x) <= window);

/**
 * Which of `indices` (the peaks in view) may carry an m/z label: the main peaks of each half, and the precursor.
 * This is the single place the label threshold lives (Figure.tsx peakLabels). The hover card of a peak is never affected.
 *  - base peak: the tallest |rel| among `indices` in the same half (rel >= 0 is the upper half, rel < 0 the lower half of a mirror plot);
 *    after a zoom it is the tallest peak in view, so a region of small peaks is still labelled;
 *  - a peak with its own colour (fragment view: every explained fragment is labelled) and a precursor peak (card kind, or a peak within
 *    `precursorWindow` of a precursor line) is always eligible;
 *  - a half without a positive base peak (all zero) has nothing to compare against: its peaks stay eligible.
 */
export function labelEligible(peaks: PeakCard[], indices: number[], precursorXs: number[] = [], minFraction = MIN_LABEL_FRACTION,
  precursorWindow = PRECURSOR_WINDOW): Set<number> {
  const half = (i: number) => (peaks[i].rel < 0 ? "lower" : "upper");
  const base = { upper: 0, lower: 0 };
  for (const i of indices) base[half(i)] = Math.max(base[half(i)], Math.abs(peaks[i].rel));
  const eligible = new Set<number>();
  for (const i of indices) {
    const p = peaks[i], top = base[half(i)];
    if (p.color || isPrecursorPeak(p, precursorXs, precursorWindow) || !(top > 0) || Math.abs(p.rel) >= minFraction * top) eligible.add(i);
  }
  return eligible;
}
