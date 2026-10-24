import { useEffect, useState } from "react";
import { BrowserRouter, Navigate, Route, Routes } from "react-router-dom";
import { getSpec } from "./api";
import { Sidebar } from "./components/Sidebar";
import { ModulePage } from "./pages/ModulePage";
import { PipelinePage } from "./pages/PipelinePage";
import type { Spec } from "./types";

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
            <Route path="/" element={<Navigate to={`/pipelines/${spec.pipelines[0].id}`} replace />} />
            <Route path="/pipelines/:id" element={<PipelinePage spec={spec} />} />
            <Route path="/modules/:id" element={<ModulePage spec={spec} />} />
            <Route path="*" element={<p className="page">Not found.</p>} />
          </Routes>
        </main>
      </div>
    </BrowserRouter>
  );
}
