import { useEffect, useState } from "react";
import { useParams } from "react-router-dom";
import { runPipeline } from "../api";
import { FieldInput, initialValues } from "../components/Fields";
import { OutputList } from "../components/Outputs";
import type { Module, Spec } from "../types";

export function PipelinePage({ spec }: { spec: Spec }) {
  const { id } = useParams();
  const pipe = spec.pipelines.find(p => p.id === id);
  const [values, setValues] = useState<Record<string, any>>(() => initialValues(pipe ? pipe.inputs : []));
  const [stages, setStages] = useState<{ title: string; status: string }[] | null>(null);
  const [msg, setMsg] = useState("");
  useEffect(() => { if (pipe) { setValues(initialValues(pipe.inputs)); setStages(null); setMsg(""); } }, [id]);
  if (!pipe) return <p>Unknown pipeline.</p>;

  const run = async () => {
    setMsg("running…");
    try { const r = await runPipeline(pipe.id, values); setStages(r.stages); setMsg("placeholder: the pipeline does not run yet; run each module from its own page"); }
    catch (e: any) { setMsg("error: " + e.message); }
  };
  return (
    <div className="page">
      <h1>{pipe.title}</h1>
      <p className="one-liner">{pipe.one_liner}</p>
      <div className="meta">{pipe.stages.map(s => s.title).join("  →  ")}</div>
      <div className="panel">
        <h4 className="panel-title">Input (once)</h4>
        <div className="fields">{pipe.inputs.map((f, i) => <FieldInput key={f.id ?? i} field={f} value={f.id ? values[f.id] : undefined} onChange={v => f.id && setValues({ ...values, [f.id]: v })} />)}</div>
        <div className="actions"><button className="btn" onClick={run}>Run pipeline</button><span className="msg">{msg}</span></div>
        {stages && <ul className="stage-status">{stages.map(s => <li key={s.title}>{s.title}: {s.status}</li>)}</ul>}
        <details><summary>Query object (shared by every stage)</summary><pre className="code">{JSON.stringify(values, null, 2)}</pre></details>
      </div>
      {pipe.stages.map((s, i) => {
        const mods = s.modules.map(id => spec.modules.find(m => m.id === id)).filter(m => m !== undefined);
        const name = (m: Module) => { const g = spec.subgroups.find(x => x.id === m.subgroup)?.title; return g ? `${g} · ${m.title}` : m.title; };
        return (
          <details className="stage" key={s.title} open={i === 0}>
            <summary>{i + 1}. {s.title} <span className="meta">{mods.length > 1 ? "modules" : "module"}: {mods.map(name).join(", ")}</span></summary>
            {mods.map(m => (
              <div key={m.id}>
                {mods.length > 1 && <h4 className="stage-module">{name(m)}</h4>}
                <OutputList blocks={m.outputs} exclude={s.exclude} />
              </div>))}
          </details>
        );
      })}
    </div>
  );
}
