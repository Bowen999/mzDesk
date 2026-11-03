import { createPortal } from "react-dom";
import { useEffect, useRef, useState, type ReactNode, type Ref } from "react";
import { useT } from "../i18n";
import {
  cancelSetup, clock, confirmSetup, dismissSetup, isFinishing, retrySetup, setupPercent, setupRows, totalSize, transferText, useSetup,
  type SetupRow, type SetupState,
} from "../setup";
import "./setup.css";

/**
 * The one dialog of first-use downloads (setup.ts): confirm -> progress, or blocked, or error. Mounted once (App.tsx).
 * Escape and the backdrop answer "no" (Not now / Close); while it downloads Escape is Cancel and the backdrop does nothing.
 * Progress: the overall bar (real fraction of the running file when the server reports it, else eased on time), and per row the step it
 * is at, with a thin bar and "412 of 2,300 MB" while a file downloads.
 */
const glyph = (body: ReactNode) => (
  <svg viewBox="0 0 24 24" width="20" height="20" fill="none" stroke="currentColor" strokeWidth="1.7" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true" focusable="false">{body}</svg>
);
const DownloadIcon = () => glyph(<><path d="M12 4v11" /><path d="m7.5 11 4.5 4.5 4.5-4.5" /><path d="M5 19.5h14" /></>);
const AlertIcon = () => glyph(<><circle cx="12" cy="12" r="9" /><path d="M12 7.5v5.5" /><path d="M12 16.4v.1" /></>);
const CheckIcon = () => (
  <svg viewBox="0 0 16 16" width="14" height="14" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true" focusable="false"><path d="m3.5 8.5 3 3 6-6.5" /></svg>
);

/** A link from the server is shown as host + path, and only when it is http(s). */
function linkLabel(link: string | null): { href: string; label: string } | null {
  if (!link || !/^https?:\/\//i.test(link)) return null;
  try { const u = new URL(link); return { href: link, label: (u.host + u.pathname).replace(/\/$/, "") }; } catch { return null; }
}

/** The raw last log line stays as a small mono line, except where a clean `step` already says it and the line is only a path or a link. */
const rawLine = (r: SetupRow) =>
  r.status === "error" ? r.detail : r.status === "running" && r.detail && (!r.step || !/[\\/]|https?:/i.test(r.detail)) && r.detail !== r.step ? r.detail : "";

function Rows({ rows }: { rows: SetupRow[] }) {
  const t = useT();
  const word = { waiting: t("Waiting"), running: t("Installing"), done: t("Done"), error: t("Failed") };
  return (
    <ul className="setup-rows">
      {rows.map(r => {
        const label = r.status === "running" && r.step ? r.step : word[r.status];
        const known = r.status === "running" && typeof r.progress === "number";
        const transfer = known ? transferText(r) : "";
        const raw = rawLine(r);
        return (
          <li key={r.id} className={"setup-row is-" + r.status}>
            <span className="setup-mark" aria-hidden="true">{r.status === "done" ? <CheckIcon /> : r.status === "error" ? "!" : null}</span>
            <span className="setup-name">{r.title}{r.size && <span className="setup-size">{r.size}</span>}</span>
            <span className="setup-word" title={label !== word[r.status] ? label : undefined}>{label}</span>
            {transfer && <span className="setup-transfer">{transfer}</span>}
            {known && <span className="setup-file" aria-hidden="true"><span style={{ width: `${Math.round(Math.min(1, r.progress as number) * 1000) / 10}%` }} /></span>}
            {raw && <span className="setup-detail" title={raw}>{raw}</span>}
          </li>);
      })}
    </ul>
  );
}

export type SetupViewProps = {
  state: SetupState; percent?: number; cardRef?: Ref<HTMLDivElement>;
  onConfirm?: () => void; onDismiss?: () => void; onRetry?: () => void; onCancel?: () => void;
};

/** The card of the dialog for the state it is in (no overlay, no focus handling: SetupDialog adds them). */
export function SetupView({ state, percent = 0, cardRef, onConfirm, onDismiss, onRetry, onCancel }: SetupViewProps) {
  const t = useT();
  const phase = state.phase;
  if (!phase) return null;
  const note = state.unsupported.map(u => (
    <p className="dialog-note setup-unsupported" key={u.id}>{t("{label}: {value}", { label: u.title, value: t(u.message) })}</p>));
  const rows = setupRows(state);
  const running = rows.findIndex(r => r.status === "running");
  const summary = running >= 0 ? t("Installing {name} ({n} of {total})", { name: rows[running].title, n: running + 1, total: rows.length }) : "";
  const total = totalSize(state.missing);
  const finishing = phase === "progress" && isFinishing(rows);
  const shownPercent = finishing ? 100 : percent;
  const tone = phase === "confirm" || phase === "progress" ? " is-info" : "";
  const props = { ref: cardRef, role: "dialog", "aria-modal": true, "aria-labelledby": "setup-title", className: "dialog-card setup-card" + tone, "data-phase": phase } as const;

  if (phase === "confirm") return (
    <div {...props} aria-describedby="setup-desc">
      <div className="dialog-badge"><DownloadIcon /></div>
      <h2 id="setup-title">{t("Download needed")}</h2>
      <p id="setup-desc">{t("mzDesk needs to download:")}</p>
      <ul className="setup-list">
        {state.missing.map(m => (
          <li key={m.id}>
            <span className="setup-name">{m.title}</span><span className="setup-size">{m.size}</span>
            {m.why && <span className="setup-why">{t(m.why)}</span>}
          </li>))}
      </ul>
      <p className="dialog-note">{t("Downloaded once, then kept on this computer.")}</p>
      <p className="dialog-note">{t("You can cancel at any time. Finished downloads are kept.")}</p>
      {note}
      <div className="dialog-actions">
        <button type="button" className="btn secondary" onClick={onDismiss}>{t("Not now")}</button>
        <button type="button" className="btn" data-autofocus onClick={onConfirm}>{total ? t("Download {size}", { size: total }) : t("Download")}</button>
      </div>
    </div>);

  if (phase === "progress") return (
    <div {...props} aria-describedby="setup-desc" tabIndex={-1}>
      <div className="dialog-badge"><DownloadIcon /></div>
      <h2 id="setup-title">{t("Setting up mzDesk")}</h2>
      <p id="setup-desc">{state.forRun ? t("Keep this page open. The run starts when this is done.") : t("Keep this page open.")}</p>
      <div className={"setup-bar-row" + (finishing ? " is-finishing" : "")}>
        <div className="setup-bar" role="progressbar" aria-label={t("Setup progress")} aria-valuemin={0} aria-valuemax={100}
          aria-valuenow={Math.round(shownPercent)} aria-valuetext={finishing ? t("Finishing…") : summary || undefined}>
          <span style={{ width: `${shownPercent}%` }} />
        </div>
        <span className="setup-pct">{finishing ? t("Finishing…") : `${Math.round(shownPercent)}%`}</span>
        <span className="setup-clock">{clock(state.job?.elapsed_s ?? 0)}</span>
      </div>
      <p className="sr-only" aria-live="polite">{finishing ? t("Finishing…") : summary}</p>
      <Rows rows={rows} />
      <p className="dialog-note">{t("You can cancel at any time. Finished downloads are kept.")}</p>
      {note}
      <div className="dialog-actions">
        <button type="button" className="btn secondary" disabled={state.cancelling} onClick={onCancel}>{state.cancelling ? t("Cancelling…") : t("Cancel")}</button>
      </div>
    </div>);

  if (phase === "blocked") return (
    <div {...props}>
      <div className="dialog-badge"><AlertIcon /></div>
      <h2 id="setup-title">{t("Cannot run yet")}</h2>
      <ul className="setup-list">
        {state.blocked.map(b => {
          const link = linkLabel(b.link);
          return (
            <li key={b.id}>
              <span className="setup-name">{b.title}</span>
              <span className="setup-why setup-message">{t(b.message)}</span>
              {link && <a className="setup-link" href={link.href} target="_blank" rel="noopener noreferrer">{link.label}</a>}
            </li>);
        })}
      </ul>
      {note}
      <div className="dialog-actions">
        <button type="button" className="btn" data-autofocus onClick={onDismiss}>{t("Close")}</button>
      </div>
    </div>);

  return (
    <div {...props}>
      <div className="dialog-badge"><AlertIcon /></div>
      <h2 id="setup-title">{t("Setup failed")}</h2>
      <div className="alert error" role="alert">{state.error}</div>
      {state.log.length > 0 && <details className="setup-log"><summary>{t("Tool log")}</summary><pre className="code">{state.log.join("\n")}</pre></details>}
      <div className="dialog-actions">
        <button type="button" className="btn secondary" onClick={onDismiss}>{t("Close")}</button>
        <button type="button" className="btn" data-autofocus onClick={onRetry}>{t("Retry")}</button>
      </div>
    </div>);
}

/** The dialog with its overlay: focus goes to the primary button of each view, Tab stays inside, the page behind does not scroll,
 *  focus returns to where it was on close. */
export function SetupDialog() {
  const state = useSetup();
  const card = useRef<HTMLDivElement>(null);
  const open = state.phase !== null;
  const downloading = state.phase === "progress";
  const [, tick] = useState(0);
  const shown = useRef(0);                         // the bar never goes back

  useEffect(() => {
    if (!downloading) { shown.current = 0; return; }
    const id = window.setInterval(() => tick(n => n + 1), 250);
    return () => window.clearInterval(id);
  }, [downloading]);

  useEffect(() => {
    if (!open) return;
    const previous = document.activeElement as HTMLElement | null;
    const overflow = document.body.style.overflow;
    document.body.style.overflow = "hidden";
    const key = (e: KeyboardEvent) => {
      if (e.key === "Escape") { e.preventDefault(); e.stopPropagation(); dismissSetup(); return; }
      if (e.key !== "Tab" || !card.current) return;
      const nodes = Array.from(card.current.querySelectorAll<HTMLElement>("button:not([disabled]), a[href], summary")).filter(n => n.offsetParent !== null);
      if (!nodes.length) { e.preventDefault(); card.current.focus(); return; }
      const first = nodes[0], last = nodes[nodes.length - 1], active = document.activeElement;
      if (e.shiftKey && (active === first || active === card.current || !card.current.contains(active))) { e.preventDefault(); last.focus(); }
      else if (!e.shiftKey && (active === last || !card.current.contains(active))) { e.preventDefault(); first.focus(); }
    };
    document.addEventListener("keydown", key, true);
    return () => {
      document.removeEventListener("keydown", key, true);
      document.body.style.overflow = overflow;
      previous?.focus?.();
    };
  }, [open]);

  // each view: its primary button (the card itself while it downloads, so a stray Enter cannot cancel)
  useEffect(() => {
    if (!state.phase) return;
    const timer = window.setTimeout(() => (card.current?.querySelector<HTMLElement>("[data-autofocus]") ?? card.current)?.focus(), 30);
    return () => window.clearTimeout(timer);
  }, [state.phase]);

  if (!state.phase) return null;
  if (downloading) shown.current = Math.max(shown.current, setupPercent(setupRows(state), state.seen, Date.now()));
  return createPortal(
    <div className="dialog-overlay setup-overlay" onMouseDown={e => { if (e.target === e.currentTarget && !downloading) dismissSetup(); }}>
      <SetupView state={state} percent={shown.current} cardRef={card} onConfirm={confirmSetup} onDismiss={dismissSetup} onRetry={retrySetup} onCancel={cancelSetup} />
    </div>,
    document.body,
  );
}
