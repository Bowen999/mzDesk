import { Fragment, useId, useState, type MouseEvent } from "react";
import type { AnnotationCandidate, Cell } from "../types";
import { candidateDetails, LEVEL_TONE, levelName } from "./annotationModel";
import { CellView } from "./Blocks";
import "./pipeline.css";

/** "View": the fragment view of this candidate (computed when clicked); a dash when it has no structure or MetFrag cannot use its adduct. */
function fragmentCell(c: AnnotationCandidate): Cell {
  return c.smiles && c.fragment_adduct ? { fragment: { smiles: c.smiles, adduct: c.fragment_adduct, name: c.name || undefined } } : "—";
}

const levelCell = (c: AnnotationCandidate): Cell =>
  c.level ? { badge: levelName(c.level_label), tone: LEVEL_TONE[c.level] ?? "neutral", help: c.level_reason } : "—";

/**
 * The ranked candidates, key columns only: Rank, Name, Level, Final score and the Fragments "View" button. Everything else (SMILES, Formula,
 * Adduct, Layer, Detectability, Isomer role) sits in a detail line that a toggle button (or a click on the row) opens, so the table never
 * scrolls sideways. The order is the one the result came in. `defaultOpen` lists candidate ids whose detail starts open.
 * The "View" buttons need a FragmentProvider above this table.
 */
export function CandidateTable({ candidates, defaultOpen = [] }: { candidates: AnnotationCandidate[]; defaultOpen?: string[] }) {
  const base = useId();
  const [open, setOpen] = useState<Set<string>>(() => new Set(defaultOpen));
  const toggle = (id: string) => setOpen(prev => { const next = new Set(prev); if (!next.delete(id)) next.add(id); return next; });
  // a click on the row opens it, unless it was meant for the button or link inside it (View, "?", the toggle itself)
  const rowClick = (id: string) => (e: MouseEvent) => { if (!(e.target as Element).closest("button, a")) toggle(id); };
  return (
    <div className="table-wrap cand-wrap">
      <table className="out-table cand-table">
        <thead><tr>
          <th className="cand-toggle-col"><span className="sr-only">Details</span></th>
          <th>Rank</th><th>Name</th><th>Level</th><th className="num">Final score</th><th>Fragments</th>
        </tr></thead>
        <tbody>
          {candidates.length ? candidates.map((c, i) => {
            const isOpen = open.has(c.id), detail = `${base}-${i}`, name = c.name || c.smiles || c.id;
            return (
              <Fragment key={c.id ?? i}>
                <tr className={"cand-row" + (isOpen ? " is-open" : "")} onClick={rowClick(c.id)}>
                  <td className="cand-toggle-col">
                    <button type="button" className="cand-toggle" aria-expanded={isOpen} aria-controls={isOpen ? detail : undefined}
                      aria-label={`${isOpen ? "Hide" : "Show"} details of ${name}`} onClick={() => toggle(c.id)}><span className="cand-caret" aria-hidden="true" /></button>
                  </td>
                  <td>{c.rank ?? "—"}</td>
                  <td className="cand-name">{c.name || "—"}</td>
                  <td><CellView cell={levelCell(c)} /></td>
                  <td className="num">{typeof c.final_score === "number" ? c.final_score.toFixed(3) : "—"}</td>
                  <td><CellView cell={fragmentCell(c)} /></td>
                </tr>
                {isOpen && (
                  <tr className="cand-detail" id={detail}>
                    <td colSpan={6}>
                      <dl className="cand-facts">{candidateDetails(c).map(([k, v]) => <div key={k}><dt>{k}</dt><dd>{v}</dd></div>)}</dl>
                    </td>
                  </tr>)}
              </Fragment>);
          }) : <tr><td colSpan={6} className="empty">no candidates</td></tr>}
        </tbody>
      </table>
    </div>
  );
}
