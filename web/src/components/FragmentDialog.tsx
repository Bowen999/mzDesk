import { createContext, useContext, useEffect, useMemo, useRef, useState, type ReactNode } from "react";
import { createPortal } from "react-dom";
import { ApiError, getFragmentJob, startFragments } from "../api";
import { useT } from "../i18n";
import { easedPercent, waitingPercent } from "../progressCurve";
import { useReportMode } from "../report/mode";
import type { FragmentItem, FragmentJob, FragmentQuery, FragmentTarget, FragmentView } from "../types";
import { FigureView, ThemedSvg } from "./Figure";
import "./fragment-view.css";

/**
 * Fragment view: the "View" button of a table row and the pop-up it opens (backend: mzdesk/ms2/fragment_view.py, mzdesk/api/fragments.py).
 * The table says which spectrum its rows share (FragmentProvider); a cell says which structure. The pop-up starts a MetFrag job when it opens,
 * polls it, shows its progress, then draws the structure with every explained fragment in its own colour, the spectrum in the same colours
 * and one card per fragment. Nothing is computed before the click.
 * Colour: the whole structure and the unmarked peaks are pure black (pure white in the dark theme); only the fragments carry a colour, and the
 * text around them stays neutral, so nothing competes with the colour blocks.
 */
type Shared = { query: FragmentQuery; note?: string };
const SharedContext = createContext<Shared | null>(null);

/** Gives the "View" cells below it their spectrum. `note` is an extra line shown in the pop-up (what differs from the page's own engine). */
export function FragmentProvider({ query, note, children }: { query?: FragmentQuery | null; note?: string; children: ReactNode }) {
  const value = useMemo(() => (query ? { query, note } : null), [query, note]);
  return <SharedContext.Provider value={value}>{children}</SharedContext.Provider>;
}

const POLL_MS = 350;
const MAX_POLL_FAILURES = 3;
const sleep = (ms: number) => new Promise<void>(resolve => window.setTimeout(resolve, ms));

function Eye() {
  return (
    <svg viewBox="0 0 16 16" width="13" height="13" fill="none" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true" focusable="false">
      <path d="M1.5 8S4 3.5 8 3.5 14.5 8 14.5 8 12 12.5 8 12.5 1.5 8 1.5 8Z" /><circle cx="8" cy="8" r="2" />
    </svg>);
}

/** The "View" cell. Without a spectrum for its table it shows a dash, so a button never opens an empty view. */
export function FragmentButton({ target }: { target: FragmentTarget }) {
  const t = useT();
  const shared = useContext(SharedContext);
  const [open, setOpen] = useState(false);
  if (useReportMode()) return null;       // the Session report is a file: the dialog needs the app
  if (!shared) return <>—</>;
  return (
    <>
      <button type="button" className="fv-open" aria-haspopup="dialog" onClick={() => setOpen(true)}
        title={t("Show which peaks this structure explains, drawn on the structure")}><Eye />{t("View")}</button>
      {open && <FragmentDialog target={target} shared={shared} onClose={() => setOpen(false)} />}
    </>);
}

/** C9H5O3 -> C₉H₅O₃ (subscript digits), with the charge as a superscript sign. */
function Formula({ text, positive }: { text: string; positive?: boolean }) {
  const parts = text.split(/(\d+)/).filter(Boolean);
  return <span className="fv-formula">{parts.map((p, i) => /^\d+$/.test(p) ? <sub key={i}>{p}</sub> : p)}{positive === undefined ? null : <sup>{positive ? "+" : "−"}</sup>}</span>;
}

/** +0.4 / −0.7; a value that rounds to zero is written 0.0, not −0.0. */
const signed = (v: number | null, digits: number) => {
  if (v === null || v === undefined) return "—";
  const text = Math.abs(v).toFixed(digits);
  return Number(text) === 0 ? text : `${v > 0 ? "+" : "−"}${text}`;
};

/** Start the job and poll it until it ends. `attempt` changes to run it again (the Try again button). */
function useFragmentJob(target: FragmentTarget, query: FragmentQuery, attempt: number) {
  const [job, setJob] = useState<FragmentJob | null>(null);
  const [error, setError] = useState<{ message: string; log: string } | null>(null);
  useEffect(() => {
    let alive = true;
    setJob(null);
    setError(null);
    (async () => {
      try {
        const { job_id } = await startFragments({ smiles: target.smiles, adduct: target.adduct, name: target.name ?? "", peaks: query.peaks,
          precursor_mz: query.precursor_mz, tol_da: query.tol_da });
        let failures = 0;
        for (;;) {
          let j: FragmentJob;
          try { j = await getFragmentJob(job_id); failures = 0; }
          catch (e: any) {
            if (!alive) return;
            if (e instanceof ApiError || ++failures >= MAX_POLL_FAILURES) throw e;
            await sleep(POLL_MS);
            continue;
          }
          if (!alive) return;
          setJob(j);
          if (j.status === "done" || j.status === "error") return;
          await sleep(POLL_MS);
          if (!alive) return;
        }
      } catch (e: any) {
        if (alive) setError({ message: e?.message || String(e), log: e instanceof ApiError ? e.log : "" });
      }
    })();
    return () => { alive = false; };
  }, [target.smiles, target.adduct, target.name, query, attempt]);       // values, not identity: the workflow page rebuilds its cells on every render
  return { job, error };
}

/** Stage and elapsed seconds (ticking). The bar is eased: the reported stage percent and time both push it, so it moves at once and keeps moving
 *  while the one long step (MetFrag) runs. */
function Progress({ job, startedAt }: { job: FragmentJob | null; startedAt: number }) {
  const t = useT();
  const [now, setNow] = useState(() => performance.now());
  useEffect(() => { const h = window.setInterval(() => setNow(performance.now()), 200); return () => window.clearInterval(h); }, []);
  const seconds = Math.max(job?.seconds ?? 0, (now - startedAt) / 1000);
  const percent = Math.round(Math.max(easedPercent(job?.percent ?? 0), waitingPercent(seconds, { start: 6, cap: 90, tau: 3 })));
  const stage = !job ? t("Starting") : job.stage;
  return (
    <div className="fv-progress" role="status" aria-live="polite">
      <div className="fv-progress-head"><span>{stage}…</span><span>{seconds.toFixed(1)} s</span></div>
      <div className="fv-bar" role="progressbar" aria-label={t("Fragment view progress")} aria-valuemin={0} aria-valuemax={100} aria-valuenow={percent}
        aria-valuetext={`${stage}`}><span style={{ width: `${percent}%` }} /></div>
    </div>);
}

function FragmentCard({ f, colors }: { f: FragmentItem; colors: string[] }) {
  const t = useT();
  return (
    <article className="fv-frag">
      <header>
        <span className="fv-dot" style={{ background: f.color }} aria-hidden="true" />
        <Formula text={f.ion_formula || "?"} />
        <span className="fv-frag-mz">{f.mz.toFixed(4)}</span>
      </header>
      <div className="fv-frag-sub">{t("Δ {ppm} ppm · {rel} % of base peak", { ppm: signed(f.ppm, 1), rel: f.rel.toFixed(1) })}</div>
      {f.svg ? <ThemedSvg svg={f.svg} keep={colors} className="fv-frag-img" alt={t("Fragment {formula} marked on the structure", { formula: f.ion_formula })} />
        : f.img ? <img className="fv-frag-img raster" src={f.img} alt={t("Fragment {formula}", { formula: f.ion_formula })} /> : null}
      {!f.placed && <p className="fv-frag-note warn">{t("Not placed on the structure: its own structure is drawn.")}</p>}
      {f.placed && f.alternatives > 0 && <p className="fv-frag-note">{f.alternatives === 1 ? t("Also fits 1 other place; the best one is drawn.") : t("Also fits {n} other places; the best one is drawn.", { n: f.alternatives })}</p>}
    </article>);
}

function Fact({ label, value, note }: { label: string; value: ReactNode; note?: string }) {
  return <div className="fv-fact"><div className="fv-fact-label">{label}</div><div className="fv-fact-value">{value}</div>{note && <div className="fv-fact-note">{note}</div>}</div>;
}

function Result({ view, extra }: { view: FragmentView; extra?: string }) {
  const t = useT();
  const { settings: s, coverage } = view;
  const colors = view.fragments.map(f => f.color);
  const lines = [...view.notes];
  if (extra) lines.push(extra);
  lines.push(t("{engine}: in-silico fragmentation, tree depth {depth}, fragment tolerance {tol} Da.", { engine: s.engine, depth: s.tree_depth, tol: s.tol_da }));
  return (
    <>
      <div className="fv-facts">
        <Fact label={t("Explained peaks")} value={`${view.explained} / ${view.peaks}`} note={t("fragment peaks; precursor and isotopes left out")} />
        <Fact label={t("Explained intensity")} value={`${(view.explained_fraction * 100).toFixed(1)} %`} note={t("of those peaks")} />
        <Fact label={t("Atoms covered")} value={`${(coverage.fraction * 100).toFixed(0)} %`}
          note={t("{atoms} of {heavy} non-hydrogen atoms, by the coloured fragments", { atoms: coverage.atoms, heavy: coverage.heavy })} />
        <Fact label={t("Ion formula")} value={<Formula text={view.ion_formula} positive={view.positive} />} note={view.formula} />
      </div>
      {view.fragments.length === 0 && <div className="alert warn">{t("MetFrag explained no peak with this structure.")}</div>}
      <div className="fv-main">
        {view.fragments.length > 0 && <figure className="fv-overview">
          <ThemedSvg svg={view.overview} keep={colors} className="fv-overview-img" alt={t("Structure with every explained fragment in its colour")} />
          <figcaption>{t("One colour per fragment; where fragments overlap, the colours share the atoms.")}</figcaption>
        </figure>}
        <div className="fv-spectrum"><FigureView data={view.figure} /></div>
      </div>
      {view.fragments.length > 0 && <div className="fv-frags">{view.fragments.map(f => <FragmentCard f={f} colors={colors} key={f.mz} />)}</div>}
      <ul className="fv-notes">{lines.map((n, i) => <li key={i}>{n}</li>)}</ul>
    </>);
}

function FragmentDialog({ target, shared, onClose }: { target: FragmentTarget; shared: Shared; onClose: () => void }) {
  const t = useT();
  const [attempt, setAttempt] = useState(0);
  const startedAt = useRef(performance.now());
  const card = useRef<HTMLDivElement>(null);
  const closeRef = useRef(onClose);
  closeRef.current = onClose;
  const { job, error } = useFragmentJob(target, shared.query, attempt);
  const failed = error ?? (job?.status === "error" ? { message: job.error ?? t("The fragment view failed."), log: job.log } : null);
  const view = job?.status === "done" ? job.result : null;

  // Esc closes, Tab stays inside, the page behind does not scroll, focus returns to the button on close
  useEffect(() => {
    const previous = document.activeElement as HTMLElement | null;
    const overflow = document.body.style.overflow;
    document.body.style.overflow = "hidden";
    card.current?.focus();
    const key = (e: KeyboardEvent) => {
      if (e.key === "Escape") { e.preventDefault(); e.stopPropagation(); closeRef.current(); return; }
      if (e.key !== "Tab" || !card.current) return;
      const nodes = Array.from(card.current.querySelectorAll<HTMLElement>("button:not([disabled]), a[href], summary")).filter(n => n.offsetParent !== null);
      if (!nodes.length) { e.preventDefault(); return; }
      const first = nodes[0], last = nodes[nodes.length - 1], active = document.activeElement;
      if (e.shiftKey && (active === first || active === card.current || !card.current.contains(active))) { e.preventDefault(); last.focus(); }
      else if (!e.shiftKey && active === last) { e.preventDefault(); first.focus(); }
    };
    document.addEventListener("keydown", key, true);
    return () => { document.removeEventListener("keydown", key, true); document.body.style.overflow = overflow; previous?.focus?.(); };
  }, []);

  const retry = () => { startedAt.current = performance.now(); setAttempt(n => n + 1); };
  const title = target.name || view?.name || target.smiles;
  return createPortal(
    <div className="dialog-overlay fv-overlay" onMouseDown={e => { if (e.target === e.currentTarget) closeRef.current(); }}>
      <div className="fv-card" role="dialog" aria-modal="true" aria-labelledby="fv-title" tabIndex={-1} ref={card}>
        <header className="fv-head">
          <div className="fv-head-text">
            <div className="fv-kicker">{t("Fragment explanation")}</div>
            <h2 id="fv-title" title={title}>{title}</h2>
            <div className="fv-sub">{view ? <><Formula text={view.formula} /> · {view.adduct} {view.theo_mz.toFixed(4)} · {view.smiles}</> : target.adduct}</div>
          </div>
          <button type="button" className="fv-close" onClick={onClose} aria-label={t("Close")}>×</button>
        </header>
        {failed ? (
          <div className="fv-body">
            <div className="alert error" role="alert">{failed.message}
              {failed.log && <details><summary>{t("Tool log")}</summary><pre className="code">{failed.log}</pre></details>}</div>
            <div className="fv-actions"><button type="button" className="btn secondary" onClick={retry}>{t("Try again")}</button></div>
          </div>
        ) : view ? <div className="fv-body"><Result view={view} extra={shared.note} /></div>
          : <div className="fv-body"><Progress job={job} startedAt={startedAt.current} /></div>}
      </div>
    </div>,
    document.body,
  );
}
