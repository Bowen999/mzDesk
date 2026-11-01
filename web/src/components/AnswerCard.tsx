import type { AnnotationResult } from "../types";
import { answerSummary } from "./annotationModel";
import { StructText, StructureHover } from "./StructureHover";
import "./pipeline.css";

/**
 * The answer of ESI MS/MS Annotation, compact: one headline ("Caffeine · Level 2a · 0.80"), at most one short line under it, and a "Details"
 * disclosure that keeps the full statement and the facts the headline leaves out. Built from the result's own fields (answerSummary), the
 * statement is the headline only when no structure, formula or mass can be named. The adduct notes and the run's warnings follow as alerts.
 */
export function AnswerCard({ result, requestedAdduct }: { result: AnnotationResult; requestedAdduct: string | null }) {
  const a = result.answer;
  const s = answerSummary(result);
  const candidates = result.candidates ?? [];
  const top = a.top_id ? candidates.find(c => c.id === a.top_id) : undefined;
  const otherAdducts = requestedAdduct ? [...new Set(candidates.filter(c => c.evidence?.library?.gate_passed && c.adduct && c.adduct !== requestedAdduct)
    .map(c => c.adduct))] : [];
  const badge = s.level && <span className={"badge tone-" + s.level.tone} title={s.level.title || undefined}>{s.level.name}</span>;
  return (
    <>
      <div className="verdict answer-card" aria-live="polite">
        <div className="verdict-label">Answer</div>
        {s.headline !== null
          ? <p className="answer-headline">
              <span className="answer-name">{top?.smiles && (a.type === "single" || a.type === "group") ? <StructureHover smiles={top.smiles}>{s.headline}</StructureHover> : s.headline}</span>
              {badge && <><span className="answer-sep" aria-hidden="true">·</span>{badge}</>}
              {s.score !== null && <><span className="answer-sep" aria-hidden="true">·</span>
                <span className="answer-score" title={`Final score ${s.score.toFixed(3)}`}><span className="sr-only">Final score </span>{s.score.toFixed(2)}</span></>}
            </p>
          : <><p className="verdict-text"><StructText text={a.statement} /></p>{badge && <div className="pill-row">{badge}</div>}</>}
        {s.secondary && <p className={"answer-sub" + (s.warn ? " warn" : "")}><StructText text={s.secondary} /></p>}
        {s.headline !== null && (
          <details className="answer-details">
            <summary>Details</summary>
            <p><StructText text={a.statement} /></p>
            {s.challengers.map((t, i) => <p key={i}><StructText text={t} />.</p>)}
            {s.facts.length > 0 && <dl className="answer-facts">{s.facts.map(([k, v]) => <div key={k}><dt>{k}</dt><dd>{v}</dd></div>)}</dl>}
          </details>)}
      </div>
      {otherAdducts.length > 0 && <div className="alert warn" role="status">Input adduct: {requestedAdduct}. A library hit passed with {otherAdducts.join(", ")}; review the library checks.</div>}
      {requestedAdduct && top?.adduct && top.adduct !== requestedAdduct && !otherAdducts.includes(top.adduct) &&
        <div className="alert warn" role="status">Input adduct: {requestedAdduct}. The top candidate uses {top.adduct}; review the adduct evidence.</div>}
      {(result.warnings ?? []).map((w, i) => <div className="alert warn" key={i}><StructText text={w} /></div>)}
    </>
  );
}
