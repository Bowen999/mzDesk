import { localizeRunBlocks, t, useT } from "../i18n";
import type { JobStep } from "../types";
import { OutputList } from "./Outputs";
import "./pipeline.css";

/** A note longer than this is cut to one line in the folded step; the whole note is the first line of the opened step. */
const SHORT_NOTE = 70;
/** Statuses that get a word next to the title; done and queued are told by the node alone. */
const wordOf = (status: string): string | null => (status === "running" ? t("running") : status === "skipped" ? t("skipped") : status === "failed" ? t("failed") : null);
/** The status as text for screen readers (the node tells it by shape only). Other words are shown as received. */
const statusText = (status: string): string => {
  switch (status) {
    case "done": return t("done");
    case "queued": return t("queued");
    case "unknown": return t("unknown");
    default: return status;
  }
};

const secondsOf = (s: JobStep) => (typeof s?.seconds === "number" && Number.isFinite(s.seconds) && s.seconds >= 0 ? s.seconds : null);

/**
 * One step of the timeline: a node (number, check, spinner, dash or "!" by status), the title with its short note under it, and the time with
 * a thin bar scaled to the longest step. A step with blocks opens on click to the tables it produced; a failed step shows its whole note.
 */
export function StepView({ step, index, longest = 0 }: { step: JobStep; index: number; longest?: number }) {
  useT();           // redraw when the language changes
  const blocks = localizeRunBlocks(step.blocks ?? []);          // the titles of the tables and notes; what the run wrote in them stays
  const outputs = Object.fromEntries(blocks.filter(b => b.id).map(b => [b.id as string, b.data]));
  const seconds = secondsOf(step);
  const status = step.status || "unknown";
  const note = step.note || "";
  const showTime = seconds !== null && status !== "queued" && status !== "skipped";
  const line = (
    <>
      <span className="tl-node" aria-hidden="true">{index + 1}</span>
      <span className="tl-main">
        <span className="tl-top">
          <span className="step-title">{step.title || t("Unnamed stage")}</span>
          {wordOf(status) ? <span className="tl-word">{wordOf(status)}</span> : <span className="sr-only">{statusText(status)}</span>}
        </span>
        {note && <span className={"step-note" + (status === "failed" ? " is-failed" : "")} title={note}>{note}</span>}
      </span>
      {showTime && (
        <span className="tl-time">
          <span className="step-seconds">{seconds!.toFixed(1)} s</span>
          {longest > 0 && <i className="tl-bar" style={{ width: `${Math.max(4, Math.min(100, (seconds! / longest) * 100)).toFixed(0)}%` }} />}
        </span>)}
    </>);
  const id = `workflow-step-${step.id}`;
  if (!blocks.length) return <li id={id} className="tl-step" data-status={status}><div className="tl-line is-static">{line}</div></li>;
  return (
    <li id={id} className="tl-step" data-status={status}>
      <details className="tl-details">
        <summary className="tl-line">{line}</summary>
        <div className="step-body">
          {note.length > SHORT_NOTE && <p className="step-full-note">{note}</p>}
          <OutputList blocks={blocks} outputs={outputs} live />
        </div>
      </details>
    </li>
  );
}

/** The steps of a job as one timeline: a rail with a node per step. All folded. */
export function StepList({ steps }: { steps: JobStep[] }) {
  const t = useT();
  const list = Array.isArray(steps) ? steps : [];
  if (!list.length) return null;
  const timed = list.filter(s => s.status !== "queued" && s.status !== "skipped").map(secondsOf).filter((v): v is number => v !== null);
  const total = timed.reduce((a, b) => a + b, 0), longest = Math.max(0, ...timed);
  const skipped = list.filter(s => s.status === "skipped").length;
  return (
    <section className="timeline" aria-label={t("Steps")}>
      <div className="tl-head">
        <h2 className="tl-title">{t("Steps")}</h2>
        <span className="tl-sum">{list.length === 1 ? t("1 step") : t("{n} steps", { n: list.length })}{skipped ? " · " + t("{n} skipped", { n: skipped }) : ""}{total > 0 ? ` · ${total.toFixed(1)} s` : ""}</span>
      </div>
      <ol className="tl-list">{list.map((s, i) => <StepView step={s} index={i} longest={longest} key={s.id} />)}</ol>
    </section>
  );
}
