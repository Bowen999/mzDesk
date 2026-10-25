import { useEffect, useState, type CSSProperties } from "react";
import { BrowserRouter, Navigate, Route, Routes, useLocation, useParams } from "react-router-dom";
import { getSpec } from "./api";
import { LockGate } from "./components/LockGate";
import { Brand, Sidebar } from "./components/Sidebar";
import { SidebarResizer, useSidebarWidth } from "./components/SidebarResizer";
import { isLockedModule, useUnlocked } from "./lock";
import { DevBlocks } from "./pages/DevBlocks";
import { ModulePage } from "./pages/ModulePage";
import { PipelinePage } from "./pages/PipelinePage";
import type { Spec } from "./types";

/** "/" lands on the first module once unlocked (Modules come first in the sidebar); while locked it lands on the first open
 *  (Detectability) page instead of a blurred wall. */
function Home({ spec }: { spec: Spec }) {
  const unlocked = useUnlocked();
  const first = unlocked ? spec.modules[0] : spec.modules.find(m => !isLockedModule(m)) ?? spec.modules[0];
  return <Navigate to={`/modules/${first.id}`} replace />;
}

function ModuleRoute({ spec }: { spec: Spec }) {
  const { id } = useParams();
  const mod = spec.modules.find(m => m.id === id);
  return <LockGate required={!!mod && isLockedModule(mod)}><ModulePage spec={spec} /></LockGate>;
}

export default function App() {
  const [spec, setSpec] = useState<Spec | null>(null);
  const [err, setErr] = useState("");
  useEffect(() => { getSpec().then(setSpec).catch(e => setErr(String(e))); }, []);
  if (err) return <p className="page">Cannot load /api/spec: {err}</p>;
  if (!spec) return <p className="page">Loading…</p>;
  return <BrowserRouter><Shell spec={spec} /></BrowserRouter>;
}

/** True while the media query matches (re-renders on change). */
function useMedia(query: string) {
  const [on, setOn] = useState(() => typeof window !== "undefined" && window.matchMedia(query).matches);
  useEffect(() => {
    const mq = window.matchMedia(query), f = () => setOn(mq.matches);
    f(); mq.addEventListener("change", f);
    return () => mq.removeEventListener("change", f);
  }, [query]);
  return on;
}

/** Sidebar (resizable on wide screens; a slide-in drawer under 960 px, opened from the top bar) + the page. */
function Shell({ spec }: { spec: Spec }) {
  const { pathname } = useLocation();
  const { width, setWidth, reset } = useSidebarWidth();
  const [drawer, setDrawer] = useState(false);
  const [collapsed, setCollapsed] = useState(() => { try { return localStorage.getItem("mzdesk.sidebarCollapsed") === "1"; } catch { return false; } });
  const narrow = useMedia("(max-width: 959px)");
  const rail = collapsed && !narrow;          // the drawer on small screens is always the full sidebar
  const collapse = (c: boolean) => { setCollapsed(c); try { localStorage.setItem("mzdesk.sidebarCollapsed", c ? "1" : "0"); } catch { /* not remembered */ } };
  useEffect(() => setDrawer(false), [pathname]);
  useEffect(() => {
    if (!drawer) return;
    const key = (e: KeyboardEvent) => { if (e.key === "Escape") setDrawer(false); };
    document.addEventListener("keydown", key);
    return () => document.removeEventListener("keydown", key);
  }, [drawer]);
  return (
    <div className={"layout sidebar-blue" + (drawer ? " drawer-open" : "") + (rail ? " sidebar-collapsed" : "")} style={{ "--sidebar-w": `${rail ? 64 : width}px` } as CSSProperties}>
      <header className="topbar">
        <button type="button" className="menu-btn" aria-label="Open navigation" aria-expanded={drawer} aria-controls="app-nav" onClick={() => setDrawer(true)}>
          <span aria-hidden="true" /><span aria-hidden="true" /><span aria-hidden="true" />
        </button>
        <Brand version={spec.version} />
      </header>
      <div className="sidebar-shell" id="app-nav">
        <Sidebar spec={spec} collapsed={rail} onCollapse={narrow ? undefined : collapse} />
        {!rail && <SidebarResizer width={width} onChange={setWidth} onReset={reset} />}
      </div>
      <div className="drawer-backdrop" aria-hidden="true" onClick={() => setDrawer(false)} />
      <main className="main">
        <Routes>
          <Route path="/" element={<Home spec={spec} />} />
          <Route path="/pipelines/:id" element={<LockGate required><PipelinePage spec={spec} /></LockGate>} />
          <Route path="/modules/:id" element={<ModuleRoute spec={spec} />} />
          {/* hidden dev page: every generic output block with sample data; not in the sidebar */}
          <Route path="/dev/blocks" element={<DevBlocks />} />
          <Route path="*" element={<p className="page">Not found.</p>} />
        </Routes>
      </main>
    </div>
  );
}
