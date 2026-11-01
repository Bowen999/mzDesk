/** Shared curve of every progress bar: fast at the start, slower towards the end. A bar never sits at 0 while work is in flight and never shows
 *  100 before the work ends. Pure functions: the components own the clock. */

export const EASE = 2;                                     // exponent of the ease-out curve: 1 is linear, larger is faster at the start
const clamp01 = (x: number) => Math.min(1, Math.max(0, x));

/** Linear share of the work (0-1) -> shown share (0-1): ease(0.1) = 0.19, ease(0.5) = 0.75. */
export const ease = (share: number) => 1 - Math.pow(1 - clamp01(share), EASE);

export type Waiting = { start?: number; cap?: number; tau?: number };

/** Work of unknown length (one module run): `start` % at once, then towards `cap` % with time constant `tau` seconds, never reaching it. */
export function waitingPercent(seconds: number, { start = 8, cap = 92, tau = 3 }: Waiting = {}): number {
  return start + (cap - start) * (1 - Math.exp(-Math.max(0, seconds) / tau));
}

/** A reported percent (0-100, stage based) -> shown percent: eased, with a floor so the bar is never empty while work runs. */
export const easedPercent = (reported: number, floor = 6) => Math.max(floor, 100 * ease(reported / 100));
