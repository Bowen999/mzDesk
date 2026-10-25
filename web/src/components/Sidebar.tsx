import { useState } from "react";
import { NavLink, useLocation } from "react-router-dom";
import { isLockedModule, lock, useUnlocked } from "../lock";
import type { Module, Spec } from "../types";
import { LockIcon } from "./LockGate";

type Entry = { kind: "link"; m: Module } | { kind: "sub"; id: string; title: string; mods: Module[] };

/** 两大类：Pipelines、Modules。Modules 按 spec 顺序列出；有小类的模块收进小类，小类出现在它第一个模块的位置。 */
function moduleEntries(spec: Spec): Entry[] {
  const titles = Object.fromEntries(spec.subgroups.map(s => [s.id, s.title]));
  const mods = spec.modules.filter(m => m.group === "modules");
  const seen = new Set<string>();
  const out: Entry[] = [];
  for (const m of mods) {
    if (!m.subgroup) out.push({ kind: "link", m });
    else if (!seen.has(m.subgroup)) {
      seen.add(m.subgroup);
      out.push({ kind: "sub", id: m.subgroup, title: titles[m.subgroup] ?? m.subgroup, mods: mods.filter(x => x.subgroup === m.subgroup) });
    }
  }
  return out;
}

function ModuleLink({ m, sub }: { m: Module; sub?: boolean }) {
  const unlocked = useUnlocked();
  const showLock = isLockedModule(m) && !unlocked;
  return (
    <NavLink to={`/modules/${m.id}`} title={m.live ? undefined : "placeholder: format only"}
      className={({ isActive }) => "nav-link" + (sub ? " sub" : "") + (m.live ? "" : " placeholder") + (isActive ? " active" : "")}>
      {m.title}{showLock && <LockIcon />}
    </NavLink>
  );
}

export function Sidebar({ spec }: { spec: Spec }) {
  const { pathname } = useLocation();
  const [open, setOpen] = useState<Record<string, boolean>>({});
  const unlocked = useUnlocked();
  return (
    <aside className="sidebar">
      <div className="lock-brand-row">
        <div className="brand">mzDesk</div>
        {unlocked && <button type="button" className="lock-btn" onClick={lock} title="Lock the password-protected pages again"><LockIcon open />Lock</button>}
      </div>
      <div className="nav-group">
        <div className="nav-title">Pipelines</div>
        {spec.pipelines.map(p => (
          <NavLink key={p.id} to={`/pipelines/${p.id}`} className={({ isActive }) => "nav-link placeholder" + (isActive ? " active" : "")}>{p.title}{!unlocked && <LockIcon />}</NavLink>
        ))}
      </div>
      <div className="nav-group">
        <div className="nav-title">Modules</div>
        {moduleEntries(spec).map(e => {
          if (e.kind === "link") return <ModuleLink key={e.m.id} m={e.m} />;
          const hasActive = e.mods.some(m => pathname === `/modules/${m.id}`);
          const isOpen = open[e.id] ?? hasActive;
          return (
            <div key={e.id} className="nav-sub">
              <button type="button" className={"nav-sub-title" + (hasActive ? " has-active" : "")} aria-expanded={isOpen}
                onClick={() => setOpen({ ...open, [e.id]: !isOpen })}>
                <span className={"caret" + (isOpen ? " open" : "")} aria-hidden="true">›</span>{e.title}{!unlocked && !isOpen && e.mods.every(isLockedModule) && <LockIcon />}
              </button>
              {isOpen && e.mods.map(m => <ModuleLink key={m.id} m={m} sub />)}
            </div>
          );
        })}
      </div>
    </aside>
  );
}
