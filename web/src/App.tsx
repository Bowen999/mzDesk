import { useEffect, useState, type CSSProperties } from "react";
import { BrowserRouter, matchPath, Navigate, Route, Routes, useLocation } from "react-router-dom";
import { getSpec, trackView } from "./api";
import { LockGate } from "./components/LockGate";
import { Brand, Sidebar } from "./components/Sidebar";
import { SidebarResizer, useSidebarWidth } from "./components/SidebarResizer";
import { isLockedModule, useUnlocked } from "./lock";
import { DevBlocks } from "./pages/DevBlocks";
import { DocsPage } from "./pages/DocsPage";
import { ModulePage } from "./pages/ModulePage";
import { PipelinePage } from "./pages/PipelinePage";
import { StatsPage } from "./pages/StatsPage";
import type { Spec } from "./types";

/** "/" lands on the first module once unlocked (Modules come first in the sidebar); while locked it lands on the first open
 *  page instead of a blurred wall. */
function Home({ spec }: { spec: Spec }) {
  const unlocked = useUnlocked();
  const first = unlocked ? spec.modules[0] : spec.modules.find(m => !isLockedModule(m)) ?? spec.modules[0];
  return <Navigate to={`/modules/${first.id}`} replace />;
}

export default function App() {
  const [spec, setSpec] = useState<Spec | null>(null);
  const [err, setErr] = useState("");
  useEffect(() => { getSpec().then(setSpec).catch(e => setErr(String(e))); }, []);
  useEffect(() => { if (spec || err) window.mzdeskReady?.(); }, [spec, err]);     // the opening splash (index.html) leaves once the page list is here
  if (err) return <p className="page">Cannot load /api/spec: {err}</p>;
  if (!spec) return <p className="page">Loading…</p>;
  return <BrowserRouter basename={import.meta.env.BASE_URL.replace(/\/$/, "") || "/"}><Shell spec={spec} /></BrowserRouter>;
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
  const moduleId = matchPath("/modules/:id", pathname)?.params.id;
  const pipelineId = matchPath("/pipelines/:id", pathname)?.params.id;
  const activeModule = spec.modules.find(m => m.id === moduleId);
  const activePipeline = spec.pipelines.find(p => p.id === pipelineId);
  const [visitedModules, setVisitedModules] = useState<string[]>([]);
  const [visitedPipelines, setVisitedPipelines] = useState<string[]>([]);
  // Keep visited pages mounted in this tab. Their inputs, results and in-flight runs remain
  // independent; inactive pages are hidden, not reused for another module's route.
  useEffect(() => {
    if (activeModule) setVisitedModules(ids => ids.includes(activeModule.id) ? ids : [...ids, activeModule.id]);
  }, [activeModule]);
  useEffect(() => {
    if (activePipeline) setVisitedPipelines(ids => ids.includes(activePipeline.id) ? ids : [...ids, activePipeline.id]);
  }, [activePipeline]);
  const { width, setWidth, reset } = useSidebarWidth();
  const [drawer, setDrawer] = useState(false);
  const [collapsed, setCollapsed] = useState(() => { try { return localStorage.getItem("mzdesk.sidebarCollapsed") === "1"; } catch { return false; } });
  const narrow = useMedia("(max-width: 959px)");
  const rail = collapsed && !narrow;          // the drawer on small screens is always the full sidebar
  const collapse = (c: boolean) => { setCollapsed(c); try { localStorage.setItem("mzdesk.sidebarCollapsed", c ? "1" : "0"); } catch { /* not remembered */ } };
  useEffect(() => setDrawer(false), [pathname]);
  useEffect(() => { window.scrollTo({ top: 0, left: 0, behavior: "auto" }); }, [pathname]);
  useEffect(() => { if (!pathname.startsWith("/stats")) trackView(pathname); }, [pathname]);   // page views for /stats
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
          <Route path="/pipelines/:id" element={activePipeline ? null : <p className="page">Unknown pipeline.</p>} />
          <Route path="/modules/:id" element={activeModule ? null : <p className="page">Unknown module.</p>} />
          <Route path="/docs/:id" element={<DocsPage spec={spec} />} />
          {/* hidden dev page: every generic output block with sample data; not in the sidebar */}
          <Route path="/dev/blocks" element={<DevBlocks />} />
          {/* usage statistics: not in the sidebar; its own password check runs on the server */}
          <Route path="/stats" element={<StatsPage spec={spec} />} />
          <Route path="*" element={<p className="page">Not found.</p>} />
        </Routes>
        {spec.modules.filter(m => visitedModules.includes(m.id) || m.id === activeModule?.id).map(m => (
          <div key={m.id} hidden={m.id !== activeModule?.id}>
            <LockGate required={isLockedModule(m)}><ModulePage id={m.id} spec={spec} active={m.id === activeModule?.id} /></LockGate>
          </div>))}
        {spec.pipelines.filter(p => visitedPipelines.includes(p.id) || p.id === activePipeline?.id).map(p => (
          <div key={p.id} hidden={p.id !== activePipeline?.id}>
            <PipelinePage id={p.id} spec={spec} active={p.id === activePipeline?.id} />
          </div>))}
      </main>
    </div>
  );
}
