import { useEffect, useLayoutEffect, useRef, useState, type KeyboardEvent, type ReactNode } from "react";
import { createPortal } from "react-dom";
import { NavLink, useLocation } from "react-router-dom";
import { LANGS, useLang, useT, type Lang, type TFn } from "../i18n";
import { isLockedModule, lock, useUnlocked } from "../lock";
import { setInternet, useInternet } from "../network";
import { useTheme, type ThemePref } from "../themeMode";
import type { Module, Spec } from "../types";
import { FeedbackButton } from "./Feedback";
import { Icon, moduleIcon, pipelineIcon } from "./Icon";
import { Logo } from "./Logo";
import { LockIcon } from "./LockGate";
import { ReportButton } from "../report/ReportButton";

type Entry = { kind: "link"; m: Module } | { kind: "sub"; id: string; title: string; mods: Module[] };

/** Modules first, then Workflows. Modules: one entry per subgroup in spec.subgroups order, with the modules in spec order;
 *  a module without a subgroup is listed on its own after them. */
function moduleEntries(spec: Spec): Entry[] {
  const mods = spec.modules.filter(m => m.group === "modules");
  const out: Entry[] = spec.subgroups
    .map(s => ({ kind: "sub" as const, id: s.id, title: s.title, mods: mods.filter(m => m.subgroup === s.id) }))
    .filter(e => e.mods.length);
  for (const m of mods) if (!m.subgroup || !spec.subgroups.some(s => s.id === m.subgroup)) out.push({ kind: "link", m });
  return out;
}

function ModuleLink({ m, sub }: { m: Module; sub?: boolean }) {
  const t = useT();
  const unlocked = useUnlocked();
  const showLock = isLockedModule(m) && !unlocked;
  return (
    <NavLink to={`/modules/${m.id}`} title={m.live ? undefined : t("placeholder: format only")}
      className={({ isActive }) => "nav-link" + (sub ? " sub" : "") + (m.live ? "" : " placeholder") + (isActive ? " active" : "")}>
      <Icon name={moduleIcon(m)} className="nav-icon" />
      <span className="nav-link-label">{m.title}</span>
      {showLock && <LockIcon />}
    </NavLink>
  );
}

/** Wordmark: "mzDesk" with "Beta" over the version, stacked to its right. */
export function Brand({ version, compact = false }: { version?: string; compact?: boolean }) {
  const t = useT();
  if (compact) return <div className="brand compact" title={`mzDesk${version ? " V" + version : ""} ${t("Beta")}`}><Logo /></div>;
  return (
    <div className="brand">
      <Logo />
      <span className="brand-name">mzDesk</span>
      <span className="brand-tag"><span>{t("Beta")}</span>{version && <span>V{version}</span>}</span>
    </div>
  );
}

const THEMES: { id: ThemePref; icon: ReactNode }[] = [
  { id: "light", icon: <path d="M12 7a5 5 0 1 0 0 10 5 5 0 0 0 0-10zM12 2v2M12 20v2M4.9 4.9l1.4 1.4M17.7 17.7l1.4 1.4M2 12h2M20 12h2M4.9 19.1l1.4-1.4M17.7 6.3l1.4-1.4" /> },
  { id: "system", icon: <><rect x="3" y="4" width="18" height="12" rx="2" /><path d="M8 20h8M12 16v4" /></> },
  { id: "dark", icon: <path d="M20 14.5A8 8 0 0 1 9.5 4a8 8 0 1 0 10.5 10.5z" /> },
];
const themeName = (t: TFn, id: ThemePref) => (id === "light" ? t("Light") : id === "dark" ? t("Dark") : t("System", { context: "theme" }));      // "System" alone is the classification table of a compound
const glyph = (body: ReactNode) => (
  <svg viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="currentColor" strokeWidth="1.7" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true">{body}</svg>
);

/** Light / System / Dark. In the collapsed rail: one button that cycles. */
function ThemeSwitch({ compact }: { compact: boolean }) {
  const t = useT();
  const { pref, setPref } = useTheme();
  if (compact) {
    const i = THEMES.findIndex(o => o.id === pref), next = THEMES[(i + 1) % THEMES.length];
    const name = themeName(t, THEMES[i].id);
    return <button type="button" className="rail-btn" title={t("Theme: {name} (click for {next})", { name, next: themeName(t, next.id) })} aria-label={t("Theme: {name}", { name })} onClick={() => setPref(next.id)}>{glyph(THEMES[i].icon)}</button>;
  }
  return (
    <div className="theme-switch" role="radiogroup" aria-label={t("Theme")}>
      {THEMES.map(o => {
        const name = themeName(t, o.id);
        return (
          <button key={o.id} type="button" role="radio" aria-checked={pref === o.id} aria-label={name} data-label={name}
            className={pref === o.id ? "on" : ""} onClick={() => setPref(o.id)}>
            {glyph(o.icon)}<span className="theme-label">{name}</span>
          </button>);
      })}
    </div>
  );
}

/** Globe, with a slash when the internet is off. */
function GlobeIcon({ off = false, size = 16 }: { off?: boolean; size?: number }) {
  return (
    <svg viewBox="0 0 24 24" width={size} height={size} fill="none" stroke="currentColor" strokeWidth="1.7"
      strokeLinecap="round" strokeLinejoin="round" aria-hidden="true" focusable="false">
      <circle cx="12" cy="12" r="9" />
      <path d="M3 12h18M12 3a15 15 0 0 1 0 18M12 3a15 15 0 0 0 0 18" />
      {off && <path d="M4 20 20 4" />}
    </svg>
  );
}

/** Internet search on/off. Off asks before a run that needs the internet (see components/NetworkDialog.tsx). Local version only:
 *  the public website (spec.hosted) leaves it out.
 *  In the collapsed rail: a globe button whose slash and dimmed colour say "off". */
function InternetSwitch({ compact = false }: { compact?: boolean }) {
  const t = useT();
  const on = useInternet();
  const toggle = () => setInternet(!on);
  const title = on
    ? t("Internet search is on. Click to turn off.")
    : t("Internet search is off: a run that needs it asks first. Click to turn on.");
  if (compact) return (
    <button type="button" className={"rail-btn internet-rail" + (on ? "" : " off")} onClick={toggle} title={title}
      aria-label={on ? t("Internet search: on") : t("Internet search: off")} aria-pressed={on}>
      <GlobeIcon off={!on} size={18} />
    </button>
  );
  return (
    <button type="button" role="switch" aria-checked={on} aria-label={t("Internet search")} className={"internet-switch" + (on ? " on" : "")}
      onClick={toggle} title={title}>
      <GlobeIcon off={!on} />
      <span className="internet-label">{t("Internet search")}</span>
      <span className="internet-track" aria-hidden="true"><span className="internet-knob" /></span>
      <span className="internet-state" aria-hidden="true">{on ? t("On") : t("Off")}</span>
    </button>
  );
}

/** Open book. */
function BookIcon({ size = 16 }: { size?: number }) {
  return (
    <svg viewBox="0 0 24 24" width={size} height={size} fill="none" stroke="currentColor" strokeWidth="1.7"
      strokeLinecap="round" strokeLinejoin="round" aria-hidden="true" focusable="false">
      <path d="M2.5 5h6a3.5 3.5 0 0 1 3.5 3.5V20a2.5 2.5 0 0 0-2.5-2.5h-7z" />
      <path d="M21.5 5h-6A3.5 3.5 0 0 0 12 8.5V20a2.5 2.5 0 0 1 2.5-2.5h7z" />
    </svg>
  );
}

/** Placeholder: the entry to the general Docs (user guide, methods, how to cite). Not wired yet. */
function DocsButton({ compact = false }: { compact?: boolean }) {
  const t = useT();
  const title = t("Docs: coming soon");
  if (compact) return (
    <button type="button" className="rail-btn soon" aria-disabled="true" title={title} aria-label={title}><BookIcon size={18} /></button>
  );
  return (
    <button type="button" className="internet-switch soon" aria-disabled="true" title={title}>
      <BookIcon />
      <span className="internet-label">{t("Docs")}</span>
      <span className="soon-tag">{t("Soon")}</span>
    </button>
  );
}

/** The language list, opened from the language button. Portalled to <body> with fixed positioning, so the sidebar's overflow never clips it
 *  (the same way the "?" popovers are placed). It opens above the button in the full sidebar and to the right of it in the rail.
 *  Keys: arrows, Home / End move, Enter / Space pick (they click the focused item), Escape and Tab close; a click outside closes too. */
function LanguageMenu({ anchor, side, lang, onPick, onClose }:
  { anchor: HTMLElement; side: "up" | "right"; lang: Lang; onPick: (l: Lang) => void; onClose: () => void }) {
  const t = useT();
  const menu = useRef<HTMLDivElement>(null);
  const [pos, setPos] = useState<{ left: number; bottom: number; maxH: number } | null>(null);

  // anchored by its bottom edge, so the menu grows upwards and never needs its own height to be placed; kept inside the window
  useLayoutEffect(() => {
    const place = () => {
      const a = anchor.getBoundingClientRect(), w = menu.current?.offsetWidth ?? 0, M = 8, GAP = 6;
      const edge = Math.max(a.right, anchor.closest("aside")?.getBoundingClientRect().right ?? 0);       // the rail's own edge, not just the button's
      const left = side === "right" ? edge + GAP : Math.max(M, Math.min(a.left, window.innerWidth - w - M));
      const bottom = side === "right" ? window.innerHeight - a.bottom : window.innerHeight - a.top + GAP;
      setPos({ left, bottom, maxH: Math.max(120, window.innerHeight - bottom - M) });
    };
    place();
    window.addEventListener("resize", place);
    return () => window.removeEventListener("resize", place);
  }, [anchor, side]);

  // focus the current language once it is placed (a hidden element cannot take focus)
  useEffect(() => { if (pos) menu.current?.querySelector<HTMLElement>('[aria-checked="true"]')?.focus(); }, [pos === null]);

  useEffect(() => {
    const away = (e: PointerEvent) => { const n = e.target as Node; if (!menu.current?.contains(n) && !anchor.contains(n)) onClose(); };
    // capture phase: Escape closes the menu only, not the navigation drawer behind it
    const key = (e: globalThis.KeyboardEvent) => { if (e.key === "Escape") { e.stopPropagation(); e.preventDefault(); onClose(); anchor.focus(); } };
    const scrolled = (e: Event) => { if (!menu.current?.contains(e.target as Node)) onClose(); };    // the button moved: the menu would float loose
    document.addEventListener("pointerdown", away);
    document.addEventListener("keydown", key, true);
    document.addEventListener("scroll", scrolled, true);
    return () => { document.removeEventListener("pointerdown", away); document.removeEventListener("keydown", key, true); document.removeEventListener("scroll", scrolled, true); };
  }, [anchor, onClose]);

  const keys = (e: KeyboardEvent<HTMLDivElement>) => {
    const items = [...e.currentTarget.querySelectorAll<HTMLElement>('[role="menuitemradio"]')];
    const at = items.indexOf(document.activeElement as HTMLElement);
    const go = (i: number) => { e.preventDefault(); items[(i + items.length) % items.length]?.focus(); };
    if (e.key === "ArrowDown") go(at + 1);
    else if (e.key === "ArrowUp") go(at <= 0 ? items.length - 1 : at - 1);
    else if (e.key === "Home") go(0);
    else if (e.key === "End") go(items.length - 1);
    else if (e.key === "Tab") { e.preventDefault(); onClose(); anchor.focus(); }
  };

  return createPortal(
    <div ref={menu} className="lang-menu" role="menu" aria-label={t("Language")} onKeyDown={keys}
      style={{ left: pos?.left ?? 0, bottom: pos?.bottom ?? 0, maxHeight: pos?.maxH, visibility: pos ? "visible" : "hidden" }}>
      {LANGS.map(l => (
        <button key={l.id} type="button" role="menuitemradio" aria-checked={l.id === lang} tabIndex={l.id === lang ? 0 : -1}
          lang={l.id === "zh" ? "zh-CN" : l.id} onClick={() => onPick(l.id)}>
          <span>{l.name}</span><span className="lang-check" aria-hidden="true">✓</span>
        </button>))}
    </div>, document.body);
}

/** The interface language: the current one as two letters (or one character); a click lists the six languages. */
function LanguageButton({ compact = false }: { compact?: boolean }) {
  const t = useT();
  const { lang, setLang } = useLang();
  const [open, setOpen] = useState(false);
  const button = useRef<HTMLButtonElement>(null);
  const current = LANGS.find(l => l.id === lang) ?? LANGS[0];
  const label = t("Language: {name}", { name: current.name });
  const close = () => setOpen(false);
  const props = {
    ref: button, type: "button" as const, "aria-haspopup": "menu" as const, "aria-expanded": open, "aria-label": label,
    onClick: () => setOpen(o => !o),
    onKeyDown: (e: KeyboardEvent<HTMLButtonElement>) => { if (e.key === "ArrowDown" || e.key === "ArrowUp") { e.preventDefault(); setOpen(true); } },
  };
  return (
    <>
      {compact
        ? <button {...props} className="rail-btn lang-btn" title={label}>{current.short}</button>
        : <button {...props} className="lang-btn" data-label={label}>{current.short}</button>}
      {open && button.current && (
        <LanguageMenu anchor={button.current} side={compact ? "right" : "up"} lang={lang} onClose={close}
          onPick={l => { setLang(l); close(); button.current?.focus(); }} />)}
    </>
  );
}

/** collapsed = icon rail: one icon per group (a click opens the sidebar on that group) and per pipeline. */
export function Sidebar({ spec, collapsed = false, onCollapse }: { spec: Spec; collapsed?: boolean; onCollapse?: (c: boolean) => void }) {
  const t = useT();
  const { pathname } = useLocation();
  const [open, setOpen] = useState<Record<string, boolean>>({});
  const unlocked = useUnlocked();
  const entries = moduleEntries(spec);
  const collapseLabel = collapsed ? t("Expand sidebar") : t("Collapse sidebar");
  const toggle = onCollapse && (
    <button type="button" className="rail-btn collapse-btn" onClick={() => onCollapse(!collapsed)}
      aria-label={collapseLabel} title={collapseLabel}>
      {glyph(<><rect x="3.5" y="4" width="17" height="16" rx="2.5" /><path d="M9.5 4v16" />{collapsed ? <path d="M13 10l2 2-2 2" /> : <path d="M15.5 10l-2 2 2 2" />}</>)}
    </button>);

  if (collapsed) return (
    <aside className="sidebar rail">
      <div className="rail-top"><Brand version={spec.version} compact />{toggle}</div>
      <nav className="rail-nav" aria-label={t("Modules")}>
        {entries.map(e => {
          if (e.kind === "link") return <NavLink key={e.m.id} to={`/modules/${e.m.id}`} title={e.m.title} className={({ isActive }) => "rail-item" + (isActive ? " active" : "")}><Icon name={moduleIcon(e.m)} className="nav-icon" /></NavLink>;
          const hasActive = e.mods.some(m => pathname === `/modules/${m.id}`);
          return (
            <button key={e.id} type="button" title={e.title} aria-label={e.title} className={"rail-item" + (hasActive ? " active" : "")}
              onClick={() => { setOpen({ ...open, [e.id]: true }); onCollapse?.(false); }}>
              <Icon name={e.id} className="nav-icon" />
            </button>);
        })}
        <span className="rail-sep" aria-hidden="true" />
        {spec.pipelines.map(p => (
          <NavLink key={p.id} to={`/pipelines/${p.id}`} title={p.title} className={({ isActive }) => "rail-item" + (isActive ? " active" : "")}>
            <Icon name={pipelineIcon(p.id)} className="nav-icon" />
          </NavLink>))}
      </nav>
      <div className="sidebar-foot">
        {!spec.hosted && <InternetSwitch compact />}<ReportButton spec={spec} rail /><DocsButton compact />
        <ThemeSwitch compact /><LanguageButton compact /><FeedbackButton rail />
      </div>
    </aside>
  );

  return (
    <aside className="sidebar">
      <div className="lock-brand-row">
        <Brand version={spec.version} />
        <div className="brand-actions">
          {unlocked && <button type="button" className="lock-btn" onClick={lock} title={t("Lock the password-protected pages again")}><LockIcon open />{t("Lock")}</button>}
          {toggle}
        </div>
      </div>
      <div className="nav-group">
        <div className="nav-title">{t("Modules")}</div>
        {entries.map(e => {
          if (e.kind === "link") return <ModuleLink key={e.m.id} m={e.m} />;
          const hasActive = e.mods.some(m => pathname === `/modules/${m.id}`);
          const isOpen = open[e.id] ?? hasActive;
          return (
            <div key={e.id} className="nav-sub">
              <button type="button" className={"nav-sub-title" + (hasActive ? " has-active" : "")} aria-expanded={isOpen}
                onClick={() => setOpen({ ...open, [e.id]: !isOpen })}>
                <Icon name={e.id} className="nav-icon" />
                <span className="nav-sub-label">{e.title}</span>
                {!unlocked && !isOpen && e.mods.every(isLockedModule) && <LockIcon />}
                <span className={"caret" + (isOpen ? " open" : "")} aria-hidden="true">›</span>
              </button>
              {isOpen && <div className="nav-sub-links">{e.mods.map(m => <ModuleLink key={m.id} m={m} sub />)}</div>}
            </div>
          );
        })}
      </div>
      <div className="nav-group">
        <div className="nav-title">{t("Workflows")}</div>
        {spec.pipelines.map(p => (
          <NavLink key={p.id} to={`/pipelines/${p.id}`} className={({ isActive }) => "nav-link" + (isActive ? " active" : "")}>
            <Icon name={pipelineIcon(p.id)} className="nav-icon" />
            <span className="nav-link-label">{p.title}</span>
          </NavLink>
        ))}
      </div>
      <div className="sidebar-foot">
        {!spec.hosted && <InternetSwitch />}
        <ReportButton spec={spec} />
        <DocsButton />
        <div className="sidebar-foot-row"><ThemeSwitch compact={false} /><LanguageButton /><FeedbackButton /></div>
      </div>
    </aside>
  );
}
