import { NavLink } from "react-router-dom";
import type { Spec } from "../types";

export function Sidebar({ spec }: { spec: Spec }) {
  const groups = [
    { title: "Pipelines", items: spec.pipelines.map(p => ({ to: `/pipelines/${p.id}`, title: p.title })) },
    { title: "Modules", items: spec.modules.filter(m => m.group === "modules").map(m => ({ to: `/modules/${m.id}`, title: m.title })) },
    { title: "DMPK", items: spec.modules.filter(m => m.group === "dmpk").map(m => ({ to: `/modules/${m.id}`, title: m.title })) },
  ];
  return (
    <aside className="sidebar">
      <div className="brand">mzDesk</div>
      <div className="brand-note">Placeholder UI: navigation, inputs and output formats only. Nothing computes.</div>
      {groups.map(g => (
        <div className="nav-group" key={g.title}>
          <div className="nav-title">{g.title}</div>
          {g.items.map(it => <NavLink key={it.to} to={it.to} className={({ isActive }) => "nav-link" + (isActive ? " active" : "")}>{it.title}</NavLink>)}
        </div>
      ))}
    </aside>
  );
}
