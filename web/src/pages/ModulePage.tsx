import { useEffect, useState } from "react";
import { useParams } from "react-router-dom";
import { ApiError, runModule } from "../api";
import { FieldInput, initialValues } from "../components/Fields";
import { OutputList } from "../components/Outputs";
import type { RunResult, Spec } from "../types";

type Run = { busy: boolean; result: RunResult | null; error: string; log: string; seconds: number | null; signature: string };
const IDLE: Run = { busy: false, result: null, error: "", log: "", seconds: null, signature: "" };

export function ModulePage({ spec }: { spec: Spec }) {
  const { id } = useParams();
  const mod = spec.modules.find(m => m.id === id);
  const [values, setValues] = useState<Record<string, any>>(() => initialValues(mod ? mod.inputs : []));
  const [run, setRun] = useState<Run>(IDLE);
  useEffect(() => { if (mod) { setValues(initialValues(mod.inputs)); setRun(IDLE); } }, [id]);
  if (!mod) return <p className="page">Unknown module.</p>;

  const subgroup = spec.subgroups.find(s => s.id === mod.subgroup)?.title;
  const signature = JSON.stringify(values);
  const go = async () => {
    const t0 = performance.now();
    setRun({ ...IDLE, busy: true });
    try {
      const result = await runModule(mod.id, values);
      setRun({ ...IDLE, result, seconds: (performance.now() - t0) / 1000, signature });
    } catch (e: any) {
      setRun({ ...IDLE, error: e.message ?? String(e), log: e instanceof ApiError ? e.log : "", signature });
    }
  };
  const ok = run.result?.status === "ok";
  return (
    <div className="page">
      <h1>{subgroup ? `${subgroup} · ${mod.title}` : mod.title}</h1>
      <div className="meta">Class {mod.cls}{mod.live ? "" : " · placeholder: shows the page format, nothing computes yet"}</div>
      <p className="one-liner">{mod.one_liner}</p>
      <div className="panel">
        <div className="fields">{mod.inputs.map((f, i) => <FieldInput key={f.id ?? i} field={f} value={f.id ? values[f.id] : undefined} onChange={v => f.id && setValues({ ...values, [f.id]: v })} />)}</div>
        <div className="actions">
          <button className="btn" onClick={go} disabled={run.busy}>{run.busy ? "Running…" : `Run ${mod.title}`}</button>
          {mod.example && <button className="btn secondary" type="button" disabled={run.busy}
            onClick={() => setValues({ ...initialValues(mod.inputs), ...mod.example })}>Load example</button>}
          <span className="msg">
            {ok && run.seconds !== null && `Done in ${run.seconds.toFixed(1)} s.`}
            {run.result?.status === "placeholder" && "Not connected yet: nothing computes on this page."}
            {(ok || run.error) && run.signature !== signature && " Inputs changed: run again to update."}
          </span>
        </div>
        {run.error && (
          <div className="alert error">{run.error}
            {run.log && <details><summary>Tool log</summary><pre className="code">{run.log}</pre></details>}
          </div>)}
        {run.result?.warnings?.map((w, i) => <div className="alert warn" key={i}>{w}</div>)}
      </div>
      <OutputList blocks={mod.outputs} outputs={ok ? run.result!.outputs : null} live={mod.live} />
    </div>
  );
}
