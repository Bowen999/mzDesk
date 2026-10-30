import type { Job } from "../types";
import "./workflow-progress.css";

const SYMBOLS: Record<string, string> = { queued: "–", running: "…", done: "✓", skipped: "–", failed: "!" };
const ENDED = new Set(["done", "skipped", "failed"]);
const KNOWN = new Set(["queued", "running", ...ENDED]);

/** Stage completion is not elapsed-time progress. A failed or skipped stage has ended too. */
export function WorkflowProgress({ job, awaitingFirstUpdate = false, trackingStopped = false }:
  { job: Job | null; awaitingFirstUpdate?: boolean; trackingStopped?: boolean }) {
  const steps = Array.isArray(job?.steps) ? job.steps : [];
  const finished = steps.filter(s => ENDED.has(s.status)).length;
  const measurable = steps.length > 0 && steps.every(s => KNOWN.has(s.status));
  const percent = measurable ? Math.round(finished * 1000 / steps.length) / 10 : null;
  const running = steps.filter(s => s.status === "running");
  const failed = steps.filter(s => s.status === "failed").length;

  let state: string;
  if (trackingStopped) state = "Progress tracking stopped. The workflow may still be running; any stages shown are the last update.";
  else if (!job) state = awaitingFirstUpdate ? "Waiting for workflow progress…" : "Submitting workflow…";
  else if (job.status === "queued") state = "Queued";
  else if (job.status === "error") state = "Workflow failed.";
  else if (job.status === "done") state = failed ? `Workflow completed with ${failed} failed ${failed === 1 ? "stage" : "stages"}.` : "Workflow completed.";
  else if (job.status !== "running") state = "Workflow status unavailable; progress tracking stopped.";
  else if (running.length) state = `Running: ${running.map(s => s.title || "Unnamed stage").join(" + ")}`;
  else if (!measurable) state = "Running; stage details unavailable.";
  else if (finished === steps.length) state = "All stages ended. Finalizing result…";
  else state = "Preparing next stage…";

  const valueText = percent === null ? undefined : `${percent}% of stages ended; ${finished} of ${steps.length}. ${state}`;
  const idle = trackingStopped || (job && !["queued", "running"].includes(job.status));
  return <section className="wf-progress" aria-labelledby="wf-progress-title">
    <div className="wf-progress-heading">
      <h2 id="wf-progress-title">Workflow progress</h2>
      <span>{percent === null ? "Stage progress unavailable" : `${percent}% · ${finished} / ${steps.length} stages ended`}</span>
    </div>
    <div className={"wf-progress-bar" + (percent === null ? " is-indeterminate" : "") + (idle ? " is-idle" : "")} role="progressbar"
      aria-labelledby="wf-progress-title" aria-valuemin={0} aria-valuemax={100} aria-valuenow={percent ?? undefined} aria-valuetext={valueText}>
      <span style={percent === null ? undefined : { width: `${percent}%` }} />
    </div>
    <p className="wf-progress-caption">Stage-count estimate: completed, skipped and failed stages count as ended; this is not elapsed-time progress.</p>
    <p className="wf-progress-state" role="status">{state}</p>
    {!trackingStopped && job?.status === "queued" && job.note && <p className="wf-progress-note">{job.note}</p>}
    {!trackingStopped && job?.status === "running" && running.filter(s => s.note).map(s =>
      <p className="wf-progress-note" key={s.id}>{s.title || "Unnamed stage"}: {s.note}</p>)}
    {steps.length > 0 && <ol className="wf-progress-track">
      {steps.map((s, i) => <li key={s.id} className={`wf-progress-node is-${s.status}`}>
        <a href={`#workflow-step-${s.id}`} aria-label={`Step ${i + 1}: ${s.title || "Unnamed stage"}, ${s.status || "unknown"}`}>
          <span className="wf-progress-symbol" aria-hidden="true">{SYMBOLS[s.status] ?? "–"}</span>
          <span className="wf-progress-node-text"><span>{i + 1}. {s.title || "Unnamed stage"}</span><small>{s.status || "unknown"}</small></span>
        </a>
      </li>)}
    </ol>}
  </section>;
}
