import { useEffect, useRef, useState } from "react";
import type { Job } from "../types";
import "./workflow-progress.css";
import { isActive, measurable, newTracker, phaseOf, track, type Phase, type StepLike, type Tracker } from "./workflowProgressModel";

const TICK_MS = 200;                 // redraw rate while the run moves; the CSS width transition smooths the steps in between
const ENDED = new Set(["done", "skipped", "failed"]);
const name = (s: StepLike) => s.title || "Unnamed stage";

/** The one-line state ("Step 4 of 9 · Ionization", "Queued", "Failed at …") and the optional short note under the bar. */
function describe(phase: Phase, job: Job | null, steps: StepLike[], usable: boolean): { label: string; note: string } {
  const running = steps.filter(s => s.status === "running");
  const failed = steps.filter(s => s.status === "failed");
  const ended = steps.filter(s => ENDED.has(s.status ?? "")).length;
  switch (phase) {
    case "stopped": return { label: "Tracking stopped", note: "The workflow may still be running; the bar shows the last update." };
    case "submitting": return { label: "Submitting workflow…", note: "" };
    case "waiting": return { label: "Waiting for workflow progress…", note: "" };
    case "queued": return { label: "Queued", note: job?.note ?? "" };
    case "done": return { label: failed.length ? `Completed with ${failed.length} failed ${failed.length === 1 ? "stage" : "stages"}` : "Completed", note: "" };
    case "error": { const at = failed[0] ?? running[0]; return { label: at ? `Failed at ${name(at)}` : "Workflow failed", note: "" }; }
    case "unknown": return { label: "Status unavailable · tracking stopped", note: "" };
    default: {
      if (!usable) return { label: "Running · stage details unavailable", note: "" };
      const noted = running.find(s => s.note);
      const note = noted ? (running.length > 1 ? `${name(noted)}: ${noted.note}` : String(noted.note)) : "";
      if (running.length === 1) return { label: `Step ${steps.indexOf(running[0]) + 1} of ${steps.length} · ${name(running[0])}`, note };
      if (running.length) return { label: `Steps ${running.map(s => steps.indexOf(s) + 1).join(", ")} of ${steps.length} · ${running.map(name).join(" + ")}`, note };
      if (ended === steps.length) return { label: "Finalizing result…", note: "" };
      return { label: ended ? "Preparing next stage…" : "Starting workflow…", note: "" };
    }
  }
}

/** One compact block: current stage and percentage on one line, a pseudo-linear bar below, at most one short note.
 *  The percentage is an estimate (see workflowProgressModel.ts): it creeps while a stage runs, glides to the next slice when one ends,
 *  never goes backwards, and reaches 100 only when the job is done. */
export function WorkflowProgress({ job, awaitingFirstUpdate = false, trackingStopped = false }:
  { job: Job | null; awaitingFirstUpdate?: boolean; trackingStopped?: boolean }) {
  const [, setTick] = useState(0);
  const tracker = useRef<Tracker | null>(null);
  const now = Date.now();
  const steps: StepLike[] = Array.isArray(job?.steps) ? job.steps : [];
  const phase = phaseOf(job, awaitingFirstUpdate, trackingStopped);
  const next = track(tracker.current ?? newTracker(now), { phase, steps }, now);
  tracker.current = next.tracker;                    // refs hold only the monotonic maximum and the first-seen times: re-rendering is idempotent
  const percent = next.percent === null ? null : Math.round(next.percent * 10) / 10;
  const shown = percent === null ? null : Math.round(percent);

  const moving = isActive(phase) && percent !== null;
  useEffect(() => {
    if (!moving) return;
    const id = setInterval(() => setTick(n => n + 1), TICK_MS);
    return () => clearInterval(id);
  }, [moving]);

  const { label, note } = describe(phase, job, steps, measurable(steps));
  const tone = phase === "done" ? " is-done" : phase === "error" ? " is-error" : phase === "stopped" || phase === "unknown" ? " is-stopped" : " is-active";
  return <section className={"wf-progress" + tone} aria-label="Workflow progress">
    <div className="wf-progress-head">
      <span className="wf-progress-stage" role="status" title={label}>{label}</span>
      <span className="wf-progress-percent">{shown === null ? "—" : `${shown}%`}</span>
    </div>
    <div className={"wf-progress-bar" + (percent === null ? " is-indeterminate" : "")} role="progressbar" aria-label="Workflow progress"
      aria-valuemin={0} aria-valuemax={100} aria-valuenow={shown ?? undefined} aria-valuetext={shown === null ? label : `${shown}%. ${label}`}>
      <span style={percent === null ? undefined : { width: `${percent}%` }} />
    </div>
    {note && <p className="wf-progress-note">{note}</p>}
  </section>;
}
