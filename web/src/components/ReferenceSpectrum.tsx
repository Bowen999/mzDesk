import { useEffect, useRef, useState, type ReactNode } from "react";
import { getReferenceSpectrum } from "../api";
import { useT } from "../i18n";
import type { ReferenceBlockData, ReferenceRecord, ReferenceView } from "../types";
import { RowsTable } from "./Blocks";
import { FigureView } from "./Figure";
import { FormatsBlock } from "./Formats";
import { helpLines } from "./Help";
import { sameTarget, useReferenceSelection } from "./referenceContext";
import { StructureHover } from "./StructureHover";
import "./fragment-view.css";
import "./reference-view.css";

/**
 * Reference Spectra, the Selected spectrum block under the table (spec.refspectrum). After a run it shows the run's record (the first row, or the
 * one named in Plot record). A click on another row of the table selects it (referenceContext.ts) and this block fetches it
 * (POST /api/reference/spectrum, mzdesk/api/reference.py) with the run's `query`: no new run, the rest of the page stays as it is. While it
 * loads, the previous spectrum stays, dimmed; a click scrolls the block into view when it is out of sight.
 * Shown: the scores against my spectrum (when the page had one), the plot (the reference, or the mirror plot against my spectrum), every field the
 * library index keeps, the peak list and the formats.
 */
const ENTROPY = "SpectralEntropy (weighted)", COSINE = "CosineGreedy";
const score = (v: number | null | undefined) => (typeof v === "number" ? v.toFixed(4) : "—");

function Fact({ label, value, note }: { label: string; value: ReactNode; note?: string }) {
  return <div className="fv-fact"><div className="fv-fact-label">{label}</div><div className="fv-fact-value">{value}</div>{note && <div className="fv-fact-note">{note}</div>}</div>;
}

/** Every field the index keeps for the record, in the order of the table, then the long identifiers. An empty field shows "—". */
function Metadata({ r }: { r: ReferenceRecord }) {
  const t = useT();
  const dash = (v: string | number | null | undefined) => (v === null || v === undefined || v === "" ? "—" : v);
  const items: [string, ReactNode, string?][] = [
    [t("Library"), r.library],
    [t("Accession"), r.url ? <a href={r.url} target="_blank" rel="noopener noreferrer">{r.accession}</a> : dash(r.accession)],
    [t("Name"), dash(r.name)],
    [t("Type"), r.kind],
    [t("Ion mode"), dash(r.ion_mode)],
    [t("Adduct"), dash(r.adduct)],
    [t("Precursor m/z"), typeof r.precursor_mz === "number" ? r.precursor_mz.toFixed(4) : "—"],
    [t("Peaks"), r.n_peaks],
    [t("Instrument"), dash(r.instrument)],
    [t("Collision energy"), dash(r.collision_energy)],
    [t("Stereo"), r.same_stereo ? t("same") : t("2D only")],
    [t("InChIKey"), dash(r.inchikey), "mono"],
    [t("SMILES"), r.smiles ? <StructureHover smiles={r.smiles}>{r.smiles}</StructureHover> : "—", "mono"],
  ];
  if (r.copy_of) items.push([t("Copy of"), r.copy_of, "mono"]);
  return (
    <dl className="rv-meta">
      {items.map(([k, v, cls]) => <div key={k} className={cls}><dt>{k}</dt><dd>{v}</dd></div>)}
    </dl>);
}

function View({ view, scored, tol }: { view: ReferenceView; scored: boolean; tol: number }) {
  const t = useT();
  const r = view.record;
  const caption = helpLines(view.spectrum.help);
  const rows = view.peak_list.map(([mz, rel]) => [mz.toFixed(4), rel.toFixed(1)]);
  const sub = [r.library, r.accession, r.kind, r.adduct, typeof r.precursor_mz === "number" ? r.precursor_mz.toFixed(4) : ""].filter(Boolean).join(" · ");
  return (
    <>
      <div className="rs-head">
        <div className="rs-name" title={r.name || r.accession}>{r.name || r.accession}</div>
        <div className="rs-sub" title={sub}>{sub}</div>
      </div>
      {scored && <div className="fv-facts rv-scores">
        <Fact label={t("Entropy")} value={score(r.scores[ENTROPY])} note={t("weighted spectral entropy")} />
        <Fact label={t("Cosine")} value={score(r.scores[COSINE])} note="CosineGreedy" />
        <Fact label={t("Matched peaks")} value={r.matches ?? "—"} note={t("of {n} reference peaks", { n: r.n_peaks })} />
        <Fact label={t("Fragment tolerance")} value={`${r.kind === "EI" ? 0.5 : tol} Da`} note={r.kind === "EI" ? t("EI always uses 0.5 Da") : undefined} />
      </div>}
      <figure className="fv-spectrum rv-spectrum">
        <FigureView data={view.spectrum} />
        {caption.length > 0 && <figcaption>{caption.join(" ")}</figcaption>}
      </figure>
      <h3 className="rv-heading">{t("Record")}</h3>
      <Metadata r={r} />
      <div className="rv-lower">
        <section className="rv-peaks">
          <h3 className="rv-heading">{t("Peak list")}<span className="count-pill">{rows.length}</span></h3>
          <div className="rv-peaks-box">
            <RowsTable columns={["m/z", "Rel. intensity (%)"]} rows={rows} empty={t("no peaks")} />
          </div>
        </section>
        <section className="rv-formats">
          <FormatsBlock block={{ type: "formats", title: t("Formats"), label: t("Formats") }} data={view.formats} empty={t("empty")} />
        </section>
      </div>
    </>);
}

export function ReferenceSpectrum({ data }: { data: ReferenceBlockData }) {
  const t = useT();
  const selection = useReferenceSelection();
  const box = useRef<HTMLDivElement>(null);
  const target = selection?.current ?? null;
  const own = !target || sameTarget(target, data.record);           // the run's own record: nothing to fetch
  const key = target ? `${target.library}\u0000${target.row}\u0000${target.accession}` : "";
  const [fetched, setFetched] = useState<{ key: string; view: ReferenceView } | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [attempt, setAttempt] = useState(0);
  const query = data.query;
  const scored = !!query?.peaks?.length;

  useEffect(() => {
    setError(null);
    if (own || !target) return;
    const abort = new AbortController();
    getReferenceSpectrum({ ...target, inchikey: query?.inchikey ?? "", compound: query?.compound ?? "", tol_da: query?.tol_da ?? 0.02,
      peaks: scored ? query!.peaks : null, precursor_mz: scored ? query!.precursor_mz : null }, abort.signal)
      .then(view => setFetched({ key, view }))
      .catch(e => { if (!abort.signal.aborted) setError(e?.message || String(e)); });
    return () => abort.abort();
  }, [key, own, data, attempt]);     // eslint-disable-line react-hooks/exhaustive-deps -- `key` stands for the target

  // a click scrolls the block into view when it is not already in sight
  useEffect(() => {
    const el = box.current?.closest<HTMLElement>(".out") ?? box.current;          // the whole block, with its title
    if (!selection?.picked || !el) return;
    const r = el.getBoundingClientRect();
    if (r.top > window.innerHeight * 0.55 || r.bottom < 120) el.scrollIntoView({ behavior: "smooth", block: "start" });
  }, [selection?.picked]);

  const current = own ? data : fetched?.key === key ? fetched.view : null;
  const shown = current ?? fetched?.view ?? data;      // while a row loads, the previous spectrum stays (dimmed)
  const loading = !current && !error;
  return (
    <div ref={box} className={"rs-panel" + (loading ? " is-loading" : "")} aria-busy={loading || undefined}>
      {error && <div className="rs-error">
        <div className="alert error" role="alert">{error}</div>
        <button type="button" className="btn secondary" onClick={() => setAttempt(n => n + 1)}>{t("Try again")}</button>
      </div>}
      {loading && <div className="rs-loading" role="status" aria-live="polite">{t("Loading the spectrum…")}</div>}
      {!error && <View view={shown} scored={scored} tol={query?.tol_da ?? 0.02} />}
    </div>);
}
