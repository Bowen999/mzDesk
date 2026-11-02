import { Fragment, useId, useState, type MouseEvent } from "react";
import { useT } from "../i18n";
import { useReportMode } from "../report/mode";
import type { AnnotationCandidate, Cell } from "../types";
import { candidateDetails, CCS_TONE, ccsNumbers, ccsText, ccsVerdict, LEVEL_TONE, levelName } from "./annotationModel";
import { CellView } from "./Blocks";
import { StructureHover } from "./StructureHover";
import "./pipeline.css";

/** "View": the fragment view of this candidate (computed when clicked); a dash when it has no structure or MetFrag cannot use its adduct. */
function fragmentCell(c: AnnotationCandidate): Cell {
  return c.smiles && c.fragment_adduct ? { fragment: { smiles: c.smiles, adduct: c.fragment_adduct, name: c.name || undefined } } : "—";
}

const levelCell = (c: AnnotationCandidate): Cell =>
  c.level ? { badge: levelName(c.level_label), tone: LEVEL_TONE[c.level] ?? "neutral", help: c.level_reason } : "—";

/**
 * The ranked candidates, key columns only: Rank, Name, Level, Final score and the Fragments "View" button. Everything else (SMILES, Formula,
 * Adduct, Layer, Detectability, Isomer role, and the CCS check when the run had an observed CCS) sits in a detail line that a toggle button (or a
 * click on the row) opens, so the table never scrolls sideways. A candidate the CCS rule applied to carries a short flag after its name. The order is the
 * one the result came in. `defaultOpen` lists candidate ids whose detail starts open.
 * The "View" buttons need a FragmentProvider above this table.
 * In the Session report (report mode) there is nothing to click: every candidate's detail line is open and the toggle and View columns are left out.
 */
export function CandidateTable({ candidates, defaultOpen = [] }: { candidates: AnnotationCandidate[]; defaultOpen?: string[] }) {
  const t = useT();
  const report = !!useReportMode();
  const base = useId();
  const [open, setOpen] = useState<Set<string>>(() => new Set(defaultOpen));
  const toggle = (id: string) => setOpen(prev => { const next = new Set(prev); if (!next.delete(id)) next.add(id); return next; });
  const cols = report ? 4 : 6;
  // a click on the row opens it, unless it was meant for the button or link inside it (View, "?", the toggle itself)
  const rowClick = (id: string) => (e: MouseEvent) => { if (!(e.target as Element).closest("button, a")) toggle(id); };
  return (
    <div className="table-wrap cand-wrap">
      <table className="out-table cand-table">
        <thead><tr>
          {!report && <th className="cand-toggle-col"><span className="sr-only">{t("Details")}</span></th>}
          <th>{t("Rank")}</th><th>{t("Name")}</th><th>{t("Level")}</th><th className="num">{t("Final score")}</th>{!report && <th>{t("Fragments")}</th>}
        </tr></thead>
        <tbody>
          {candidates.length ? candidates.map((c, i) => {
            const isOpen = report || open.has(c.id), detail = `${base}-${i}`, name = c.name || c.smiles || c.id;
            return (
              <Fragment key={c.id ?? i}>
                <tr className={"cand-row" + (isOpen ? " is-open" : "")} onClick={rowClick(c.id)}>
                  {!report && <td className="cand-toggle-col">
                    <button type="button" className="cand-toggle" aria-expanded={isOpen} aria-controls={isOpen ? detail : undefined}
                      aria-label={isOpen ? t("Hide details of {name}", { name }) : t("Show details of {name}", { name })} onClick={() => toggle(c.id)}><span className="cand-caret" aria-hidden="true" /></button>
                  </td>}
                  <td>{c.rank ?? "—"}</td>
                  <td className="cand-name">{c.name ? (c.smiles ? <StructureHover smiles={c.smiles}>{c.name}</StructureHover> : c.name) : "—"}{c.ccs?.flag && <span className="cand-flag" title={ccsText(c.ccs)}>{c.ccs.flag}</span>}</td>
                  <td><CellView cell={levelCell(c)} /></td>
                  <td className="num">{typeof c.final_score === "number" ? c.final_score.toFixed(3) : "—"}</td>
                  {!report && <td><CellView cell={fragmentCell(c)} /></td>}
                </tr>
                {isOpen && (
                  <tr className="cand-detail" id={detail}>
                    <td colSpan={cols}>
                      <dl className="cand-facts">
                        {candidateDetails(c).map(([k, v]) => <div key={k}><dt>{k}</dt><dd>{k === "SMILES" && c.smiles ? <StructureHover smiles={c.smiles}>{v}</StructureHover> : v}</dd></div>)}
                        {c.ccs && <div key="CCS"><dt>CCS</dt>
                          <dd><span className={"badge tone-" + (CCS_TONE[c.ccs.verdict] ?? "neutral")}>{ccsVerdict(c.ccs.verdict)}</span>
                            {ccsNumbers(c.ccs) && <span className="cand-ccs"> {ccsNumbers(c.ccs)}</span>}
                            {c.ccs.moved_from ? <span className="cand-ccs"> · {t("moved from rank {n}", { n: c.ccs.moved_from })}</span> : null}</dd></div>}
                      </dl>
                    </td>
                  </tr>)}
              </Fragment>);
          }) : <tr><td colSpan={cols} className="empty">{t("no candidates")}</td></tr>}
        </tbody>
      </table>
    </div>
  );
}
