import { useEffect, useState } from "react";
import { BrowserRouter, Navigate, Route, Routes, useParams } from "react-router-dom";
import { getSpec } from "./api";
import { LockGate } from "./components/LockGate";
import { Sidebar } from "./components/Sidebar";
import { isLockedModule, useUnlocked } from "./lock";
import { DevBlocks } from "./pages/DevBlocks";
import { ModulePage } from "./pages/ModulePage";
import { PipelinePage } from "./pages/PipelinePage";
import type { Spec } from "./types";

/** "/" lands on the first pipeline once unlocked; while locked it lands on the first open (Detectability) page instead of a blurred wall. */
function Home({ spec }: { spec: Spec }) {
  const unlocked = useUnlocked();
  const open = spec.modules.find(m => !isLockedModule(m));
  return <Navigate to={!unlocked && open ? `/modules/${open.id}` : `/pipelines/${spec.pipelines[0].id}`} replace />;
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
  return (
    <BrowserRouter>
      <div className="layout">
        <Sidebar spec={spec} />
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
    </BrowserRouter>
  );
}
