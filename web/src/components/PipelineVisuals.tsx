import { useState } from "react";
import { useT } from "../i18n";
import { useReportMode } from "../report/mode";
import type { AnnotationResult } from "../types";
import { FigureView } from "./Figure";
import { isomerWord } from "./annotationModel";
import { HelpButton } from "./Help";
import { StructText, StructureHover } from "./StructureHover";
import "./pipeline.css";

/**
 * The measured spectrum against the best judged library hit. Only shown when a library search finished and there is a hit to compare:
 * without a search (skipped, not installed, no peaks) or without a hit there is nothing to put next to it, and the reason stays in the step's note.
 * The one line of explanation lives behind the "?" next to the title; the plot is lower than on the module pages.
 */
export function SpectrumEvidence({ result }: { result: AnnotationResult }) {
  const t = useT();
  const report = useReportMode();       // the report has no tabs: the library comparison, then the experimental spectrum
  const figures = result.visualization;
  const [view, setView] = useState<"query" | "library">(figures?.library ? "library" : "query");
  const libraryStep = result.steps?.find(s => s.id === 3);
  const query = figures?.query;
  const library = figures?.library;
  if (!library || (libraryStep && libraryStep.status !== "done")) return null;
  const showLibrary = !query || view === "library";
  const help = [t("Top: experimental. Bottom: library. Dark peaks are matched. Precursor m/z {mz} · {polarity}.", { mz: result.context.precursor_mz.toFixed(4), polarity: result.context.polarity })];
  return <section className="workflow-evidence spectrum-evidence" aria-labelledby="workflow-spectrum-title">
    <div className="workflow-heading"><h2 id="workflow-spectrum-title">{t("Spectrum evidence")}</h2><HelpButton lines={help} label={t("Help: Spectrum evidence")} />
      {query && !report && <div className="workflow-tabs" role="group" aria-label={t("Spectrum view")}>
        <button type="button" className={showLibrary ? "" : "selected"} aria-pressed={!showLibrary} onClick={() => setView("query")}>{t("Experimental spectrum")}</button>
        <button type="button" className={showLibrary ? "selected" : ""} aria-pressed={showLibrary} onClick={() => setView("library")}>{t("Library comparison")}</button>
      </div>}
    </div>
    {showLibrary
      ? <>
        <p className="workflow-caption spectrum-hit">{t("Best library hit:")} <b>{library.smiles ? <StructureHover smiles={library.smiles}>{library.name || library.accession}</StructureHover> : <StructText text={library.name || library.accession} />}</b> · {library.library} · {library.accession}</p>
        <FigureView data={library.figure} />
      </>
      : query && <FigureView data={query} />}
    {report && showLibrary && query && <>
      <p className="workflow-caption">{t("Experimental spectrum")}</p>
      <FigureView data={query} />
    </>}
  </section>;
}

export function IsomerEvidence({ result }: { result: AnnotationResult }) {
  const t = useT();
  const iso = result.isomer;
  const candidates = iso?.candidates ?? [];
  if (!iso || !candidates.length) return null;
  const group = new Set(iso.group);
  const challengers = new Set(iso.challengers.map(c => c.id));
  const others = new Map((iso.isomers ?? []).map(c => [c.id, c]));
  return <section className="workflow-evidence" aria-labelledby="workflow-isomer-title">
    <div className="workflow-heading"><h2 id="workflow-isomer-title">{t("Isomer comparison")}<HelpButton lines={[t("Share of the fragment intensity that MS-FINDER explains for each candidate.")]} label={t("Help: Isomer comparison")} /></h2>
      <span className="badge">{isomerWord(iso.verdict)}</span></div>
    <div className="isomer-bars">
      {candidates.map(c => {
        const fraction = c.explained_fraction;
        const other = others.get(c.id);
        const role = c.id === iso.top_id ? t("top-ranked") : !c.evaluated ? t("not scored") : group.has(c.id) ? t("not separated") : other?.distinguished ? t("separated") : t("not compared");
        return <div className="isomer-bar-row" key={c.id}>
          <div className="isomer-bar-label"><span>{c.smiles ? <StructureHover smiles={c.smiles}>{c.name || c.smiles || c.id}</StructureHover> : (c.name || c.id)}</span>
            <small>{role}{challengers.has(c.id) ? " · " + t("challenger") : ""}{c.source !== "input" ? ` · ${c.source}` : ""}</small></div>
          <div className="isomer-bar-track" aria-hidden="true"><span style={{ width: `${Math.max(0, Math.min(100, (fraction ?? 0) * 100))}%` }} /></div>
          <span className="isomer-bar-value">{fraction === null ? t("Not scored") : `${(fraction * 100).toFixed(1)}%`}</span>
          {other?.distinguishing_peaks?.length ? <small className="isomer-peaks">{t("Top-only m/z: {peaks}", { peaks: other.distinguishing_peaks.map(p => p.mz.toFixed(4)).join(", ") })}</small> : null}
        </div>;
      })}
    </div>
  </section>;
}
