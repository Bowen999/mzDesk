import { useState, type ReactNode } from "react";
import { NavLink, useLocation } from "react-router-dom";
import { isLockedModule, lock, useUnlocked } from "../lock";
import { setInternet, useInternet } from "../network";
import { useTheme, type ThemePref } from "../themeMode";
import type { Module, Spec } from "../types";
import { FeedbackButton } from "./Feedback";
import { Icon, moduleIcon } from "./Icon";
import { Logo } from "./Logo";
import { LockIcon } from "./LockGate";

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
  const unlocked = useUnlocked();
  const showLock = isLockedModule(m) && !unlocked;
  return (
    <NavLink to={`/modules/${m.id}`} title={m.live ? undefined : "placeholder: format only"}
      className={({ isActive }) => "nav-link" + (sub ? " sub" : "") + (m.live ? "" : " placeholder") + (isActive ? " active" : "")}>
      <Icon name={moduleIcon(m)} className="nav-icon" />
      <span className="nav-link-label">{m.title}</span>
      {showLock && <LockIcon />}
    </NavLink>
  );
}

/** Wordmark: "mzDesk" with "Beta" over the version, stacked to its right. */
export function Brand({ version, compact = false }: { version?: string; compact?: boolean }) {
  if (compact) return <div className="brand compact" title={`mzDesk${version ? " V" + version : ""} Beta`}><Logo /></div>;
  return (
    <div className="brand">
      <Logo />
      <span className="brand-name">mzDesk</span>
      <span className="brand-tag"><span>Beta</span>{version && <span>V{version}</span>}</span>
    </div>
  );
}

const THEMES: { id: ThemePref; label: string; icon: ReactNode }[] = [
  { id: "light", label: "Light", icon: <path d="M12 7a5 5 0 1 0 0 10 5 5 0 0 0 0-10zM12 2v2M12 20v2M4.9 4.9l1.4 1.4M17.7 17.7l1.4 1.4M2 12h2M20 12h2M4.9 19.1l1.4-1.4M17.7 6.3l1.4-1.4" /> },
  { id: "system", label: "System", icon: <><rect x="3" y="4" width="18" height="12" rx="2" /><path d="M8 20h8M12 16v4" /></> },
  { id: "dark", label: "Dark", icon: <path d="M20 14.5A8 8 0 0 1 9.5 4a8 8 0 1 0 10.5 10.5z" /> },
];
const glyph = (body: ReactNode) => (
  <svg viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="currentColor" strokeWidth="1.7" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true">{body}</svg>
);

/** Light / System / Dark. In the collapsed rail: one button that cycles. */
function ThemeSwitch({ compact }: { compact: boolean }) {
  const { pref, setPref } = useTheme();
  if (compact) {
    const i = THEMES.findIndex(t => t.id === pref), next = THEMES[(i + 1) % THEMES.length];
    return <button type="button" className="rail-btn" title={`Theme: ${THEMES[i].label} (click for ${next.label})`} aria-label={`Theme: ${THEMES[i].label}`} onClick={() => setPref(next.id)}>{glyph(THEMES[i].icon)}</button>;
  }
  return (
    <div className="theme-switch" role="radiogroup" aria-label="Theme">
      {THEMES.map(t => (
        <button key={t.id} type="button" role="radio" aria-checked={pref === t.id} aria-label={t.label} data-label={t.label}
          className={pref === t.id ? "on" : ""} onClick={() => setPref(t.id)}>
          {glyph(t.icon)}<span className="theme-label">{t.label}</span>
        </button>))}
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

/** Internet search on/off. Off asks before a run that needs the internet (see components/NetworkDialog.tsx).
 *  In the collapsed rail: a globe button whose slash and dimmed colour say "off". */
function InternetSwitch({ compact = false }: { compact?: boolean }) {
  const on = useInternet();
  const toggle = () => setInternet(!on);
  const title = on
    ? "Internet search is on: online lookups run without asking. Click to turn off."
    : "Internet search is off: a run that needs the internet asks first. Click to turn on.";
  if (compact) return (
    <button type="button" className={"rail-btn internet-rail" + (on ? "" : " off")} onClick={toggle} title={title}
      aria-label={`Internet search: ${on ? "on" : "off"}`} aria-pressed={on}>
      <GlobeIcon off={!on} size={18} />
    </button>
  );
  return (
    <button type="button" role="switch" aria-checked={on} aria-label="Internet search" className={"internet-switch" + (on ? " on" : "")}
      onClick={toggle} title={title}>
      <GlobeIcon off={!on} />
      <span className="internet-label">Internet search</span>
      <span className="internet-track" aria-hidden="true"><span className="internet-knob" /></span>
      <span className="internet-state" aria-hidden="true">{on ? "On" : "Off"}</span>
    </button>
  );
}

/** collapsed = icon rail: one icon per group (a click opens the sidebar on that group) and per pipeline. */
export function Sidebar({ spec, collapsed = false, onCollapse }: { spec: Spec; collapsed?: boolean; onCollapse?: (c: boolean) => void }) {
  const { pathname } = useLocation();
  const [open, setOpen] = useState<Record<string, boolean>>({});
  const unlocked = useUnlocked();
  const entries = moduleEntries(spec);
  const toggle = onCollapse && (
    <button type="button" className="rail-btn collapse-btn" onClick={() => onCollapse(!collapsed)}
      aria-label={collapsed ? "Expand sidebar" : "Collapse sidebar"} title={collapsed ? "Expand sidebar" : "Collapse sidebar"}>
      {glyph(<><rect x="3.5" y="4" width="17" height="16" rx="2.5" /><path d="M9.5 4v16" />{collapsed ? <path d="M13 10l2 2-2 2" /> : <path d="M15.5 10l-2 2 2 2" />}</>)}
    </button>);

  if (collapsed) return (
    <aside className="sidebar rail">
      <div className="rail-top"><Brand version={spec.version} compact />{toggle}</div>
      <nav className="rail-nav" aria-label="Modules">
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
            <Icon name="pipeline" className="nav-icon" />
          </NavLink>))}
      </nav>
      <div className="sidebar-foot"><InternetSwitch compact /><ThemeSwitch compact /><FeedbackButton rail /></div>
    </aside>
  );

  return (
    <aside className="sidebar">
      <div className="lock-brand-row">
        <Brand version={spec.version} />
        <div className="brand-actions">
          {unlocked && <button type="button" className="lock-btn" onClick={lock} title="Lock the password-protected pages again"><LockIcon open />Lock</button>}
          {toggle}
        </div>
      </div>
      <div className="nav-group">
        <div className="nav-title">Modules</div>
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
        <div className="nav-title">Workflows</div>
        {spec.pipelines.map(p => (
          <NavLink key={p.id} to={`/pipelines/${p.id}`} className={({ isActive }) => "nav-link" + (isActive ? " active" : "")}>
            <Icon name="pipeline" className="nav-icon" />
            <span className="nav-link-label">{p.title}</span>
          </NavLink>
        ))}
      </div>
      <div className="sidebar-foot">
        <InternetSwitch />
        <div className="sidebar-foot-row"><ThemeSwitch compact={false} /><FeedbackButton /></div>
      </div>
    </aside>
  );
}
