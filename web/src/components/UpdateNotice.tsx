import { useEffect, useLayoutEffect, useRef, useState, type KeyboardEvent } from "react";
import { createPortal } from "react-dom";
import { Link, useInRouterContext } from "react-router-dom";
import { BASE, type UpdateInfo } from "../api";
import { useT } from "../i18n";
import { dismissUpdate, useUpdate } from "../setup";
import "./setup.css";

/**
 * "Version 0.5.5 is available": a quiet item in the sidebar footer (local version only) that opens a small popover with the changelog link,
 * the command `mzdesk update` and [Later]. [Later] is remembered per version (setup.ts). In the collapsed rail: an icon with a dot.
 * The popover is portalled to <body> and placed like the language list (Sidebar.tsx), so the sidebar's overflow never clips it.
 */
const CHANGELOG = "docs/changelog";

function UpdateIcon({ size = 16 }: { size?: number }) {
  return (
    <span className="update-icon">
      <svg viewBox="0 0 24 24" width={size} height={size} fill="none" stroke="currentColor" strokeWidth="1.7" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true" focusable="false">
        <circle cx="12" cy="12" r="9" /><path d="M12 16.5v-8" /><path d="m8.5 11.8 3.5-3.5 3.5 3.5" />
      </svg>
      <span className="update-dot" aria-hidden="true" />
    </span>);
}

function CopyCommand({ text }: { text: string }) {
  const t = useT();
  const [done, setDone] = useState(false);
  const copy = () => { navigator.clipboard?.writeText(text).then(() => { setDone(true); window.setTimeout(() => setDone(false), 1200); }).catch(() => undefined); };
  return <button type="button" className="copy" onClick={copy}>{done ? t("Copied") : t("Copy")}</button>;
}

/** The content of the popover. */
export function UpdateCard({ info, onLater, onNavigate }: { info: UpdateInfo; onLater?: () => void; onNavigate?: () => void }) {
  const t = useT();
  const inRouter = useInRouterContext();
  const notes = t("What changed");
  return (
    <>
      <p className="update-title">{t("Version {v} is available", { v: info.latest ?? "" })}</p>
      {inRouter
        ? <Link className="update-notes" to={`/${CHANGELOG}`} onClick={onNavigate}>{notes}</Link>
        : <a className="update-notes" href={`${BASE}${CHANGELOG}`}>{notes}</a>}
      <p className="update-hint">{t("In a terminal:")}</p>
      <div className="update-cmd"><code>{info.command || "mzdesk update"}</code><CopyCommand text={info.command || "mzdesk update"} /></div>
      <div className="update-actions"><button type="button" className="btn ghost" onClick={onLater}>{t("Later")}</button></div>
    </>);
}

function UpdatePopover({ anchor, side, info, onClose, onLater }:
  { anchor: HTMLElement; side: "up" | "right"; info: UpdateInfo; onClose: () => void; onLater: () => void }) {
  const t = useT();
  const pop = useRef<HTMLDivElement>(null);
  const [pos, setPos] = useState<{ left: number; bottom: number; maxH: number } | null>(null);

  // anchored by its bottom edge, so it grows upwards and never needs its own height to be placed; kept inside the window
  useLayoutEffect(() => {
    const place = () => {
      const a = anchor.getBoundingClientRect(), w = pop.current?.offsetWidth ?? 0, M = 8, GAP = 6;
      const edge = Math.max(a.right, anchor.closest("aside")?.getBoundingClientRect().right ?? 0);
      const left = side === "right" ? edge + GAP : Math.max(M, Math.min(a.left, window.innerWidth - w - M));
      const bottom = side === "right" ? window.innerHeight - a.bottom : window.innerHeight - a.top + GAP;
      setPos({ left, bottom, maxH: Math.max(120, window.innerHeight - bottom - M) });
    };
    place();
    window.addEventListener("resize", place);
    return () => window.removeEventListener("resize", place);
  }, [anchor, side]);

  // focus the first control once it is placed (a hidden element cannot take focus)
  useEffect(() => { if (pos) pop.current?.querySelector<HTMLElement>("a, button")?.focus(); }, [pos === null]);

  useEffect(() => {
    const away = (e: PointerEvent) => { const n = e.target as Node; if (!pop.current?.contains(n) && !anchor.contains(n)) onClose(); };
    // capture phase: Escape closes the popover only, not the navigation drawer behind it
    const key = (e: globalThis.KeyboardEvent) => { if (e.key === "Escape") { e.stopPropagation(); e.preventDefault(); onClose(); anchor.focus(); } };
    const scrolled = (e: Event) => { if (!pop.current?.contains(e.target as Node)) onClose(); };
    document.addEventListener("pointerdown", away);
    document.addEventListener("keydown", key, true);
    document.addEventListener("scroll", scrolled, true);
    return () => { document.removeEventListener("pointerdown", away); document.removeEventListener("keydown", key, true); document.removeEventListener("scroll", scrolled, true); };
  }, [anchor, onClose]);

  // Tab stays inside the popover (it sits at the end of <body>, so the page order would lead out of it)
  const keys = (e: KeyboardEvent<HTMLDivElement>) => {
    if (e.key !== "Tab") return;
    const nodes = [...e.currentTarget.querySelectorAll<HTMLElement>("a, button")];
    const first = nodes[0], last = nodes[nodes.length - 1], active = document.activeElement;
    if (e.shiftKey && active === first) { e.preventDefault(); last.focus(); }
    else if (!e.shiftKey && active === last) { e.preventDefault(); first.focus(); }
  };

  return createPortal(
    <div ref={pop} className="update-pop" role="dialog" aria-label={t("Update available")} onKeyDown={keys}
      style={{ left: pos?.left ?? 0, bottom: pos?.bottom ?? 0, maxHeight: pos?.maxH, visibility: pos ? "visible" : "hidden" }}>
      <UpdateCard info={info} onLater={onLater} onNavigate={onClose} />
    </div>, document.body);
}

/** The sidebar footer item. Nothing when there is no newer version or the person put this one off. */
export function UpdateNotice({ compact = false }: { compact?: boolean }) {
  const t = useT();
  const info = useUpdate();
  const [open, setOpen] = useState(false);
  const button = useRef<HTMLButtonElement>(null);
  if (!info?.latest) return null;
  const label = t("Version {v} is available", { v: info.latest });
  const later = () => {
    // the item goes away: keep the keyboard in the footer
    const foot = button.current?.closest(".sidebar-foot");
    setOpen(false);
    dismissUpdate(info.latest!);
    window.requestAnimationFrame(() => foot?.querySelector<HTMLElement>("button, a")?.focus());
  };
  const props = { ref: button, type: "button" as const, "aria-haspopup": "dialog" as const, "aria-expanded": open, onClick: () => setOpen(o => !o) };
  return (
    <>
      {compact
        ? <button {...props} className="rail-btn update-rail" title={label} aria-label={label}><UpdateIcon size={18} /></button>
        : <button {...props} className="internet-switch update-entry" title={label}><UpdateIcon /><span className="internet-label">{label}</span></button>}
      {open && button.current && (
        <UpdatePopover anchor={button.current} side={compact ? "right" : "up"} info={info} onClose={() => setOpen(false)} onLater={later} />)}
    </>);
}
