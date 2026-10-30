import { useState } from "react";
import type { AnnotationResult, Job } from "../types";
import { FigureView } from "./Figure";
import "./pipeline.css";

const SYMBOLS: Record<string, string> = { queued: "–", running: "…", done: "✓", skipped: "–", failed: "!" };

/** Counts completed steps, not elapsed time: steps 1 and 2 may run together. */
export function WorkflowOverview({ job }: { job: Job }) {
  const finished = job.steps.filter(s => ["done", "failed", "skipped"].includes(s.status)).length;
  return <section className="workflow-overview" aria-label="Workflow progress">
    <div className="workflow-heading"><h2>Workflow progress</h2><span>{finished} / {job.steps.length} steps finished</span></div>
    {job.status === "queued" && job.note && <p className="workflow-note">{job.note}</p>}
    <ol className="workflow-track">
      {job.steps.map((s, i) => <li key={s.id} className={`workflow-node is-${s.status}`}>
        <a href={`#workflow-step-${s.id}`} aria-label={`Step ${i + 1}: ${s.title}, ${s.status}`}>
          <span className="workflow-symbol" aria-hidden="true">{SYMBOLS[s.status] ?? "–"}</span>
          <span className="workflow-node-text"><span>{i + 1}. {s.title}</span><small>{s.status}</small></span>
        </a>
      </li>)}
    </ol>
  </section>;
}

export function SpectrumEvidence({ result }: { result: AnnotationResult }) {
  const figures = result.visualization;
  const [view, setView] = useState<"query" | "library">(figures?.library ? "library" : "query");
  if (!figures?.query && !figures?.library) return null;
  const library = figures.library;
  const showLibrary = view === "library" && !!library;
  return <section className="workflow-evidence" aria-labelledby="workflow-spectrum-title">
    <div className="workflow-heading"><div><h2 id="workflow-spectrum-title">Spectrum evidence</h2>
      <p>Measured m/z {result.context.precursor_mz.toFixed(4)} · {result.context.polarity}</p></div></div>
    {figures.query && library && <div className="workflow-tabs" role="group" aria-label="Spectrum view">
      <button type="button" className={showLibrary ? "" : "selected"} aria-pressed={!showLibrary} onClick={() => setView("query")}>Experimental spectrum</button>
      <button type="button" className={showLibrary ? "selected" : ""} aria-pressed={showLibrary} onClick={() => setView("library")}>Library comparison</button>
    </div>}
    {showLibrary && library ? <>
      <p className="workflow-caption">Best judged library hit: <b>{library.name || library.accession}</b> · {library.library} · {library.accession}. This may not be the top-ranked final candidate.</p>
      <FigureView data={library.figure} />
      <p className="workflow-caption">Top: experimental. Bottom: library. Dark peaks are paired by greedy matching; the library similarity uses its own metric.</p>
    </> : figures.query && <>
      <FigureView data={figures.query} />
      <p className="workflow-caption">Parsed experimental spectrum. Library search uses all peaks when run; engine steps use at most the 100 most intense peaks.</p>
      {!library && <p className="workflow-caption">No library spectrum is available for comparison.</p>}
    </>}
  </section>;
}

export function IsomerEvidence({ result }: { result: AnnotationResult }) {
  const iso = result.isomer;
  const candidates = iso?.candidates ?? [];
  if (!iso || !candidates.length) return null;
  const group = new Set(iso.group);
  const challengers = new Set(iso.challengers.map(c => c.id));
  const others = new Map((iso.isomers ?? []).map(c => [c.id, c]));
  return <section className="workflow-evidence" aria-labelledby="workflow-isomer-title">
    <div className="workflow-heading"><div><h2 id="workflow-isomer-title">Isomer comparison</h2>
      <p>MS-FINDER fragment explanations for the candidates considered</p></div><span className="badge">{iso.verdict}</span></div>
    <div className="isomer-bars">
      {candidates.map(c => {
        const fraction = c.explained_fraction;
        const other = others.get(c.id);
        const role = c.id === iso.top_id ? "top-ranked" : !c.evaluated ? "not scored" : group.has(c.id) ? "not separated" : other?.distinguished ? "separated" : "not compared";
        return <div className="isomer-bar-row" key={c.id}>
          <div className="isomer-bar-label"><span title={c.smiles}>{c.name || c.smiles || c.id}</span>
            <small>{role}{challengers.has(c.id) ? " · challenger" : ""}{c.source !== "input" ? ` · ${c.source}` : ""}</small></div>
          <div className="isomer-bar-track" aria-hidden="true"><span style={{ width: `${Math.max(0, Math.min(100, (fraction ?? 0) * 100))}%` }} /></div>
          <span className="isomer-bar-value">{fraction === null ? "Not scored" : `${(fraction * 100).toFixed(1)}%`}</span>
          {other?.distinguishing_peaks?.length ? <small className="isomer-peaks">Top-only m/z: {other.distinguishing_peaks.map(p => p.mz.toFixed(4)).join(", ")}</small> : null}
        </div>;
      })}
    </div>
    <p className="workflow-caption">Bars show the share of considered fragment intensity explained by the model, not identification confidence. Unlisted isomers are not ruled out.</p>
  </section>;
}
