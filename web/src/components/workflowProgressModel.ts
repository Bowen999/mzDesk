/** Pseudo-linear progress for a workflow run. Pure functions only: the component owns the clock, the refs and the render.
 *
 *  The server reports stage statuses, never a percentage, so the bar is an estimate that must never look stuck, never jump and never move backwards:
 *   - the first PRE_CEIL % belong to "submitting / waiting / queued": a slow creep from PRE_START towards PRE_CEIL, so the bar is never at 0;
 *   - the rest is shared by the stages in proportion to their expected duration (`expected_s` of the spec; equal shares when no stage has one).
 *     A stage that is skipped takes no share. A stage that has ended contributes its whole share, a running stage contributes
 *     RUN_CAP * (1 - exp(-t / tau)) of it, where t is the time since the UI first saw it running and tau = TAU_PER_EXPECTED * expected_s * k
 *     (k: how much slower or faster the finished stages ran than expected, clamped; without expected_s the old rule: mean seconds of the finished
 *     stages, TAU_DEFAULT before any). Several running stages add up; a queued stage contributes 0;
 *   - done is 100; every other state is capped at 99 so 100 means "the result is in";
 *   - error / stopped tracking / unknown status freeze the bar where it is;
 *   - the target only ever rises, and the bar follows it at no more than MAX_RATE points per second, so what the estimate cannot know in advance
 *     (a skipped stage that was expected to be long, a stage that ends early) becomes a quick glide instead of a jump. Done is shown at once. */

export type Phase = "submitting" | "waiting" | "queued" | "running" | "done" | "error" | "stopped" | "unknown";
export type StepLike = { id?: number | string | null; title?: string | null; status?: string | null; seconds?: number | null; note?: string | null;
  expected_s?: number | null };
/** since: when this run was first drawn; seen: when each stage was first seen running; max: the highest target so far;
 *  shown / at: the value on screen and when it was set (null before the first sample); last: the previous phase. */
export type Tracker = { since: number; seen: Record<string, number>; max: number; shown: number; at: number | null; last: Phase | null };

export const PRE_START = 2;        // % shown the moment the run is submitted
export const PRE_CEIL = 6;         // % the pre-run creep approaches; stages fill the band between this and 100
export const PRE_TAU = 20;         // seconds: time constant of the pre-run creep
export const RUN_CAP = 0.9;        // share of its slice a running stage can fill before it ends
export const TAU_DEFAULT = 8;      // seconds: stage time constant until a stage has finished
export const TAU_MIN = 2;          // clamp of the fallback tau (no expected_s)
export const TAU_MAX = 60;
export const TAU_PER_EXPECTED = 1;     // tau = expected_s: a stage that takes its expected time has filled 57% of its slice (cap 90%)
export const TAU_EXPECTED_MIN = 0.5;   // clamp of tau when it comes from expected_s
export const TAU_EXPECTED_MAX = 300;
export const SPEED_MIN = 0.25;         // clamp of k, the observed / expected duration ratio of the finished stages
export const SPEED_MAX = 4;
export const SPEED_PRIOR_S = 3;        // seconds added to both sides of that ratio so a few milliseconds of noise do not move it
export const MAX_RATE = 15;            // points per second the bar may advance towards its target
export const ACTIVE_MAX = 99;      // only a finished job shows 100

const ENDED = new Set(["done", "skipped", "failed"]);
const KNOWN = new Set(["queued", "running", ...ENDED]);
const finite = (n: unknown): n is number => typeof n === "number" && Number.isFinite(n);
export const stepKey = (s: StepLike, i: number) => String(s.id ?? i);

export function phaseOf(job: { status?: string | null } | null | undefined, awaitingFirstUpdate: boolean, trackingStopped: boolean): Phase {
  if (trackingStopped) return "stopped";
  if (!job) return awaitingFirstUpdate ? "waiting" : "submitting";
  return job.status === "queued" || job.status === "running" || job.status === "done" || job.status === "error" ? job.status : "unknown";
}

/** A run is moving (and the bar should be redrawn on a timer) while it is submitted, waiting, queued or running. */
export const isActive = (p: Phase) => p === "submitting" || p === "waiting" || p === "queued" || p === "running";

/** Stage data can drive the bar only when it exists and every status is one we know. */
export const measurable = (steps: readonly StepLike[]) => steps.length > 0 && steps.every(s => KNOWN.has(s.status ?? ""));

/** Fallback time constant of a running stage without expected_s: mean seconds of the stages that finished, clamped; TAU_DEFAULT when none has. */
export function stageTau(steps: readonly StepLike[]): number {
  const t = steps.filter(s => s.status === "done" && finite(s.seconds) && s.seconds > 0).map(s => s.seconds as number);
  const mean = t.length ? t.reduce((a, b) => a + b, 0) / t.length : TAU_DEFAULT;
  return Math.min(TAU_MAX, Math.max(TAU_MIN, mean));
}

const expectedOf = (s: StepLike): number | null => (finite(s.expected_s) && s.expected_s > 0 ? s.expected_s : null);

/** Share of each stage: its expected_s; a stage without one gets the mean of the others; equal shares when none has any. */
export function stageWeights(steps: readonly StepLike[]): number[] {
  const known = steps.map(expectedOf);
  const present = known.filter((v): v is number => v !== null);
  if (!present.length) return steps.map(() => 1);
  const mean = present.reduce((a, b) => a + b, 0) / present.length;
  return known.map(v => v ?? mean);
}

/** How much slower (>1) or faster (<1) the finished stages ran than their expected_s, clamped; 1 before any stage has finished. */
export function speedFactor(steps: readonly StepLike[]): number {
  let seen = 0, expected = 0;
  for (const s of steps) {
    const e = expectedOf(s);
    if (s.status === "done" && e !== null && finite(s.seconds) && s.seconds >= 0) { seen += s.seconds; expected += e; }
  }
  return Math.min(SPEED_MAX, Math.max(SPEED_MIN, (seen + SPEED_PRIOR_S) / (expected + SPEED_PRIOR_S)));
}

/** Time constant of one running stage: from its expected_s (scaled by the observed speed), else the fallback. */
export function tauOf(step: StepLike, steps: readonly StepLike[], speed = speedFactor(steps)): number {
  const e = expectedOf(step);
  return e === null ? stageTau(steps) : Math.min(TAU_EXPECTED_MAX, Math.max(TAU_EXPECTED_MIN, TAU_PER_EXPECTED * e * speed));
}

/** Before the first stage runs: PRE_START creeping towards PRE_CEIL, never reaching it. */
export const creep = (elapsedSeconds: number) => PRE_START + (PRE_CEIL - PRE_START) * (1 - Math.exp(-Math.max(0, elapsedSeconds) / PRE_TAU));

/** Share of its slice a stage fills after `elapsedSeconds` of running: 0 at the start, approaching RUN_CAP, never 1. */
export const stageFraction = (elapsedSeconds: number, tau: number) => RUN_CAP * (1 - Math.exp(-Math.max(0, elapsedSeconds) / tau));

export const newTracker = (now: number): Tracker => ({ since: now, seen: {}, max: 0, shown: 0, at: null, last: null });

type Input = { phase: Phase; steps: readonly StepLike[] };

/** Target in percent (0-100), or null when a running job has no usable stage data (indeterminate). Ignores the tracker's max and the slew limit. */
export function estimateProgress(input: Input & { since: number; seen: Record<string, number> }, nowMs: number): number | null {
  const { phase, steps } = input;
  if (phase === "done") return 100;
  if (phase === "submitting" || phase === "waiting" || phase === "queued") return Math.min(PRE_CEIL, creep((nowMs - input.since) / 1000));
  if (!measurable(steps)) return null;
  const weights = stageWeights(steps), speed = speedFactor(steps);
  let total = 0, units = 0;
  steps.forEach((s, i) => {
    if (s.status === "skipped") return;                                   // takes no share of the bar
    total += weights[i];
    if (s.status === "done" || s.status === "failed") units += weights[i];
    else if (phase === "running" && s.status === "running")
      units += weights[i] * stageFraction((nowMs - (input.seen[stepKey(s, i)] ?? nowMs)) / 1000, tauOf(s, steps, speed));
  });
  if (total <= 0) return ACTIVE_MAX;
  return Math.min(ACTIVE_MAX, PRE_CEIL + (100 - PRE_CEIL) * units / total);
}

/** Advance the tracker with the current job state at time `nowMs`. Returns the new tracker, the percent to show (null = indeterminate) and the
 *  target it is heading for. The target never decreases within a run and the percent follows it at most MAX_RATE points per second (done: at once).
 *  A new run (a submission after a finished/stopped one, or a finished job replaced by a queued/running one) starts over. */
export function track(prev: Tracker, input: Input, nowMs: number): { tracker: Tracker; percent: number | null; target: number | null } {
  const { phase, steps } = input;
  const preRun = phase === "submitting" || phase === "waiting";
  const fresh = (preRun && prev.last !== null && prev.last !== "submitting" && prev.last !== "waiting")
    || ((phase === "queued" || phase === "running") && (prev.last === "done" || prev.last === "error"));
  const base: Tracker = fresh ? newTracker(nowMs) : prev;
  const seen = { ...base.seen };
  if (phase === "running") steps.forEach((s, i) => { if (s.status === "running" && seen[stepKey(s, i)] === undefined) seen[stepKey(s, i)] = nowMs; });
  const keep = (shown: number, max: number): Tracker => ({ since: base.since, seen, max, shown, at: nowMs, last: phase });
  if (phase === "done") return { tracker: keep(100, 100), percent: 100, target: 100 };
  const est = estimateProgress({ phase, steps, since: base.since, seen }, nowMs);
  if (!isActive(phase)) {                                                  // error / stopped / unknown: stay where the bar is (first drawn here: the ended stages)
    const shown = Math.min(ACTIVE_MAX, Math.max(base.shown, base.at === null ? est ?? PRE_START : PRE_START));
    return { tracker: keep(shown, Math.max(base.max, shown)), percent: shown, target: shown };
  }
  if (est === null) return { tracker: { ...keep(base.shown, base.max), at: base.at }, percent: null, target: null };
  const target = Math.max(base.max, est);
  const reach = base.at === null ? target : base.shown + MAX_RATE * Math.max(0, nowMs - base.at) / 1000;
  const shown = Math.min(target, Math.max(base.shown, reach));
  return { tracker: keep(shown, target), percent: shown, target };
}
