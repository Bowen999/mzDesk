import { useEffect, useState } from "react";
import { useParams } from "react-router-dom";
import { runModule } from "../api";
import { FieldInput, initialValues } from "../components/Fields";
import { OutputList } from "../components/Outputs";
import type { Spec } from "../types";

export function ModulePage({ spec }: { spec: Spec }) {
  const { id } = useParams();
  const mod = spec.modules.find(m => m.id === id);
  const [values, setValues] = useState<Record<string, any>>(() => initialValues(mod ? mod.inputs : []));
  const [msg, setMsg] = useState("");
  useEffect(() => { if (mod) { setValues(initialValues(mod.inputs)); setMsg(""); } }, [id]);
  if (!mod) return <p>Unknown module.</p>;

  const run = async () => {
    setMsg("running…");
    try { const r = await runModule(mod.id, values); setMsg(`${mod.title}: ${r.status}. Nothing computes yet; inputs were received by the API.`); }
    catch (e: any) { setMsg("error: " + e.message); }
  };
  return (
    <div className="page">
      <h1>{mod.title}</h1>
      <div className="meta">Class {mod.cls}</div>
      <p className="one-liner">{mod.one_liner}</p>
      <div className="panel">
        <div className="fields">{mod.inputs.map((f, i) => <FieldInput key={f.id ?? i} field={f} value={f.id ? values[f.id] : undefined} onChange={v => f.id && setValues({ ...values, [f.id]: v })} />)}</div>
        <div className="actions"><button className="btn" onClick={run}>Run {mod.title}</button><span className="msg">{msg}</span></div>
      </div>
      <OutputList blocks={mod.outputs} />
    </div>
  );
}
