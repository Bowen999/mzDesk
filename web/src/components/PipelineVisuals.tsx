import { useState } from "react";
import { useReportMode } from "../report/mode";
import type { AnnotationResult } from "../types";
import { FigureView } from "./Figure";
import { HelpButton } from "./Help";
import { StructText, StructureHover } from "./StructureHover";
import "./pipeline.css";

/**
 * The measured spectrum against the best judged library hit. Only shown when a library search finished and there is a hit to compare:
 * without a search (skipped, not installed, no peaks) or without a hit there is nothing to put next to it, and the reason stays in the step's note.
 * The explanations live behind the "?" next to the title; the plot is lower than on the module pages.
 */
export function SpectrumEvidence({ result }: { result: AnnotationResult }) {
  const report = useReportMode();       // the report has no tabs: the library comparison, then the experimental spectrum
  const figures = result.visualization;
  const [view, setView] = useState<"query" | "library">(figures?.library ? "library" : "query");
  const libraryStep = result.steps?.find(s => s.id === 3);
  const query = figures?.query;
  const library = figures?.library;
  if (!library || (libraryStep && libraryStep.status !== "done")) return null;
  const showLibrary = !query || view === "library";
  const help = [
    `Measured m/z ${result.context.precursor_mz.toFixed(4)} · ${result.context.polarity}.`,
    "Top: experimental. Bottom: library. Dark peaks are paired by greedy matching; the library similarity uses its own metric.",
    "The best judged library hit may not be the top-ranked final candidate.",
    "Only the main peaks (at least 5 % of the base peak) and the precursor are labelled; hover any peak for its details.",
    ...(libraryStep ? [`Library search: ${libraryStep.status}.${libraryStep.note ? " " + libraryStep.note : ""}`] : []),
  ];
  return <section className="workflow-evidence spectrum-evidence" aria-labelledby="workflow-spectrum-title">
    <div className="workflow-heading"><h2 id="workflow-spectrum-title">Spectrum evidence</h2><HelpButton lines={help} label="Help: Spectrum evidence" />
      {query && !report && <div className="workflow-tabs" role="group" aria-label="Spectrum view">
        <button type="button" className={showLibrary ? "" : "selected"} aria-pressed={!showLibrary} onClick={() => setView("query")}>Experimental spectrum</button>
        <button type="button" className={showLibrary ? "selected" : ""} aria-pressed={showLibrary} onClick={() => setView("library")}>Library comparison</button>
      </div>}
    </div>
    {showLibrary
      ? <>
        <p className="workflow-caption spectrum-hit">Best library hit: <b>{library.smiles ? <StructureHover smiles={library.smiles}>{library.name || library.accession}</StructureHover> : <StructText text={library.name || library.accession} />}</b> · {library.library} · {library.accession}</p>
        <FigureView data={library.figure} />
      </>
      : query && <FigureView data={query} />}
    {report && showLibrary && query && <>
      <p className="workflow-caption">Experimental spectrum</p>
      <FigureView data={query} />
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
          <div className="isomer-bar-label"><span>{c.smiles ? <StructureHover smiles={c.smiles}>{c.name || c.smiles || c.id}</StructureHover> : (c.name || c.id)}</span>
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
