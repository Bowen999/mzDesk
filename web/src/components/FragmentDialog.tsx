import { createContext, useContext, useEffect, useMemo, useRef, useState, type ReactNode } from "react";
import { createPortal } from "react-dom";
import { ApiError, getFragmentJob, startFragments } from "../api";
import { darkColors } from "../darkColors";
import { useTheme } from "../themeMode";
import type { FragmentItem, FragmentJob, FragmentQuery, FragmentTarget, FragmentView } from "../types";
import { FigureView, ThemedSvg } from "./Figure";
import "./fragment-view.css";

/**
 * Fragment view: the "View" button of a table row and the pop-up it opens (backend: mzdesk/ms2/fragment_view.py, mzdesk/api/fragments.py).
 * The table says which spectrum its rows share (FragmentProvider); a cell says which structure. The pop-up starts a MetFrag job when it opens,
 * polls it, shows its progress, then draws the structure with every explained fragment in its own colour, the spectrum in the same colours and
 * one card per fragment. Nothing is computed before the click.
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
  const shared = useContext(SharedContext);
  const [open, setOpen] = useState(false);
  if (!shared) return <>—</>;
  return (
    <>
      <button type="button" className="fv-open" aria-haspopup="dialog" onClick={() => setOpen(true)}
        title="Show which peaks this structure explains, drawn on the structure"><Eye />View</button>
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

/** Start the job and poll it until it ends. `attempt` changes to run it again (the Retry button). */
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

/** Elapsed seconds while `running`, ticking; the job's own figure once it has ended. */
function Progress({ job, startedAt }: { job: FragmentJob | null; startedAt: number }) {
  const [now, setNow] = useState(() => performance.now());
  useEffect(() => { const h = window.setInterval(() => setNow(performance.now()), 200); return () => window.clearInterval(h); }, []);
  const percent = job?.percent ?? 0;
  const seconds = Math.max(job?.seconds ?? 0, (now - startedAt) / 1000);
  const stage = !job ? "Starting" : job.stage;
  return (
    <div className="fv-progress" role="status" aria-live="polite">
      <div className="fv-progress-head"><span>{stage}…</span><span>{percent}% · {seconds.toFixed(1)} s</span></div>
      <div className="fv-bar" role="progressbar" aria-label="Fragment view progress" aria-valuemin={0} aria-valuemax={100} aria-valuenow={percent}
        aria-valuetext={`${percent}%, ${stage}`}><span style={{ width: `${Math.max(percent, 2)}%` }} /></div>
      <p className="fv-progress-note">Stage-based: MetFrag is one step, so the bar waits there until it returns (usually 2–5 s).</p>
    </div>);
}

function FragmentCard({ f }: { f: FragmentItem }) {
  const { dark } = useTheme();
  const tone = dark ? darkColors(f.color) : f.color;
  return (
    <article className="fv-frag" style={{ borderTopColor: tone }}>
      <h4 style={{ color: tone }}><Formula text={f.ion_formula || "?"} /><span className="fv-frag-mz">m/z {f.mz.toFixed(4)}</span></h4>
      <div className="fv-frag-sub">Δ {signed(f.ppm, 1)} ppm · {f.rel.toFixed(1)} %</div>
      {f.svg ? <ThemedSvg svg={f.svg} className="fv-frag-img" alt={`Fragment ${f.ion_formula} marked on the structure`} />
        : f.img ? <img className="fv-frag-img raster" src={f.img} alt={`Fragment ${f.ion_formula}`} /> : null}
      {!f.placed && <p className="fv-frag-note warn">Not placed on the structure: its own structure is drawn.</p>}
      {f.placed && f.alternatives > 0 && <p className="fv-frag-note">Also fits {f.alternatives} other place{f.alternatives === 1 ? "" : "s"}; the best one is drawn.</p>}
    </article>);
}

function Stat({ label, value, note }: { label: string; value: ReactNode; note?: string }) {
  return <div className="fv-stat"><div className="fv-stat-label">{label}</div><div className="fv-stat-value">{value}</div>{note && <div className="fv-stat-note">{note}</div>}</div>;
}

function Result({ view, extra }: { view: FragmentView; extra?: string }) {
  const { settings: s } = view;
  const lines = [...view.notes];
  if (extra) lines.push(extra);
  lines.push(`${s.engine}: in-silico fragmentation, tree depth ${s.tree_depth}, fragment tolerance ${s.tol_da} Da. The fragments are predictions, not proof of identity; an unexplained peak is not evidence against the structure.`);
  return (
    <>
      <div className="fv-stats">
        <Stat label="Explained peaks" value={`${view.explained} / ${view.peaks}`} note="fragment peaks, precursor and isotopes left out" />
        <Stat label="Explained intensity" value={`${(view.explained_fraction * 100).toFixed(1)} %`} note="of those peaks" />
        <Stat label="Precursor m/z" value={view.theo_mz.toFixed(4)} note={view.adduct} />
        <Stat label="Ion formula" value={<Formula text={view.ion_formula} positive={view.positive} />} note={view.formula} />
      </div>
      {view.fragments.length === 0 && <div className="alert warn">MetFrag explained no peak with this structure, so there is nothing to mark. The spectrum is shown below.</div>}
      <div className="fv-main">
        {view.fragments.length > 0 && <figure className="fv-overview">
          <ThemedSvg svg={view.overview} className="fv-overview-img" alt="Structure with every explained fragment in its colour" />
          <figcaption>Each colour is one fragment; where fragments overlap the colours share the atoms.</figcaption>
        </figure>}
        <div className="fv-spectrum"><FigureView data={view.figure} /></div>
      </div>
      {view.fragments.length > 0 && <div className="fv-frags">{view.fragments.map(f => <FragmentCard f={f} key={f.mz} />)}</div>}
      <ul className="fv-notes">{lines.map((n, i) => <li key={i}>{n}</li>)}</ul>
    </>);
}

function FragmentDialog({ target, shared, onClose }: { target: FragmentTarget; shared: Shared; onClose: () => void }) {
  const [attempt, setAttempt] = useState(0);
  const startedAt = useRef(performance.now());
  const card = useRef<HTMLDivElement>(null);
  const closeRef = useRef(onClose);
  closeRef.current = onClose;
  const { job, error } = useFragmentJob(target, shared.query, attempt);
  const failed = error ?? (job?.status === "error" ? { message: job.error ?? "The fragment view failed.", log: job.log } : null);
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
    <div className="dialog-overlay fv-overlay" onMouseDown={e => { if (e.target === e.currentTarget) onClose(); }}>
      <div className="fv-card" role="dialog" aria-modal="true" aria-labelledby="fv-title" tabIndex={-1} ref={card}>
        <header className="fv-head">
          <div className="fv-head-text">
            <div className="fv-kicker">Fragment explanation · {target.adduct}</div>
            <h2 id="fv-title" title={title}>{title}</h2>
            {view && <div className="fv-sub"><Formula text={view.formula} /> · {view.smiles}</div>}
          </div>
          <button type="button" className="fv-close" onClick={onClose} aria-label="Close">×</button>
        </header>
        {failed ? (
          <div className="fv-body">
            <div className="alert error" role="alert">{failed.message}
              {failed.log && <details><summary>Tool log</summary><pre className="code">{failed.log}</pre></details>}</div>
            <div className="fv-actions"><button type="button" className="btn secondary" onClick={retry}>Try again</button></div>
          </div>
        ) : view ? <div className="fv-body"><Result view={view} extra={shared.note} /></div>
          : <div className="fv-body"><Progress job={job} startedAt={startedAt.current} /></div>}
      </div>
    </div>,
    document.body,
  );
}
