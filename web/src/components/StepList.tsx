import type { JobStep } from "../types";
import { OutputList } from "./Outputs";
import "./pipeline.css";

const STATUS_TONE: Record<string, string> = { done: "good", running: "jp", failed: "bad", skipped: "neutral", queued: "neutral" };
/** A note longer than this is cut to one line in the collapsed step; the whole note is the first line of the opened step. */
const SHORT_NOTE = 70;

export function StatusChip({ status }: { status: string }) {
  return <span className={"badge tone-" + (STATUS_TONE[status] ?? "neutral")}>{status}</span>;
}

/**
 * One step, folded: a single line (number, title, status, seconds, note) that opens on click to the blocks the step has produced (tables with ids,
 * filled as soon as the step ends). A failed step shows its whole note in the folded line. A step with nothing to open (no blocks) is the line alone.
 */
export function StepView({ step, index }: { step: JobStep; index: number }) {
  const blocks = step.blocks ?? [];
  const outputs = Object.fromEntries(blocks.filter(b => b.id).map(b => [b.id as string, b.data]));
  const hasSeconds = typeof step.seconds === "number" && Number.isFinite(step.seconds) && step.seconds >= 0;
  const status = step.status || "unknown";
  const note = step.note || "";
  const line = (
    <>
      <span className="step-title">{index + 1}. {step.title || "Unnamed stage"}</span>
      <StatusChip status={status} />
      {status !== "queued" && status !== "skipped" && hasSeconds && <span className="step-seconds">{step.seconds.toFixed(1)} s</span>}
      {note && <span className={"step-note" + (status === "failed" ? " is-failed" : "")} title={note}>{note}</span>}
    </>);
  const id = `workflow-step-${step.id}`;
  if (!blocks.length) return <div id={id} className="stage step-line is-static" data-status={status}>{line}</div>;
  return (
    <details id={id} className="stage step-line" data-status={status}>
      <summary>{line}</summary>
      <div className="step-body">
        {note.length > SHORT_NOTE && <p className="step-full-note">{note}</p>}
        <OutputList blocks={blocks} outputs={outputs} live />
      </div>
    </details>
  );
}

/** The steps of a job, all folded. */
export function StepList({ steps }: { steps: JobStep[] }) {
  return <>{(Array.isArray(steps) ? steps : []).map((s, i) => <StepView step={s} index={i} key={s.id} />)}</>;
}
