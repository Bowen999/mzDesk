import { useEffect, useRef, useState } from "react";
import { useParams } from "react-router-dom";
import { ApiError, runModule } from "../api";
import { FieldInput, initialValues, summarizeInputs } from "../components/Fields";
import { RichText } from "../components/Help";
import { Icon, moduleIcon } from "../components/Icon";
import { OutputList } from "../components/Outputs";
import type { RunResult, Spec } from "../types";

type Run = { busy: boolean; result: RunResult | null; error: string; log: string; seconds: number | null; signature: string };
const IDLE: Run = { busy: false, result: null, error: "", log: "", seconds: null, signature: "" };

/** Seconds since `busy` turned on, ticking; null when idle. */
function useElapsed(busy: boolean) {
  const [t, setT] = useState<number | null>(null);
  useEffect(() => {
    if (!busy) { setT(null); return; }
    const t0 = performance.now();
    setT(0);
    const h = window.setInterval(() => setT((performance.now() - t0) / 1000), 200);
    return () => window.clearInterval(h);
  }, [busy]);
  return t;
}

export function ModulePage({ spec }: { spec: Spec }) {
  const { id } = useParams();
  const mod = spec.modules.find(m => m.id === id);
  const [values, setValues] = useState<Record<string, any>>(() => initialValues(mod ? mod.inputs : []));
  const [run, setRun] = useState<Run>(IDLE);
  const [editing, setEditing] = useState(true);     // false = inputs folded into a one-line summary (after a successful run)
  const ticket = useRef(0);       // a slow run (ESI MS/MS takes up to minutes) must not land on the page the user has moved to
  const elapsed = useElapsed(run.busy);
  useEffect(() => { ticket.current++; if (mod) { setValues(initialValues(mod.inputs)); setRun(IDLE); setEditing(true); } }, [id]);
  if (!mod) return <p className="page">Unknown module.</p>;

  const subgroup = spec.subgroups.find(s => s.id === mod.subgroup)?.title;
  const detectability = mod.subgroup === "detectability";
  const signature = JSON.stringify(values);
  const go = async () => {
    const t0 = performance.now(), mine = ++ticket.current;
    setRun(r => ({ ...IDLE, busy: true, result: r.result, signature: r.signature }));
    try {
      const result = await runModule(mod.id, values);
      if (ticket.current !== mine) return;
      setRun({ ...IDLE, result, seconds: (performance.now() - t0) / 1000, signature });
      if (result.status === "ok") setEditing(false);
    } catch (e: any) {
      if (ticket.current === mine) { setRun({ ...IDLE, error: e.message ?? String(e), log: e instanceof ApiError ? e.log : "", signature }); setEditing(true); }
    }
  };
  const ok = run.result?.status === "ok" && !run.busy;
  const folded = ok && !editing;
  const summary = summarizeInputs(mod.inputs, values);
  const stale = (ok || run.error) && run.signature !== signature;
  const runLabel = run.busy ? <><span className="spinner" aria-hidden="true" />Running{elapsed !== null && elapsed >= 1 ? ` · ${Math.floor(elapsed)} s` : "…"}</> : `Run ${mod.title}`;

  return (
    <div className={"page module-page" + (detectability ? " detectability-page" : "")} data-module={mod.id}>
      <div className="module-heading">
        {detectability && <span className="module-emblem"><Icon name={moduleIcon(mod)} /></span>}
        <div className="module-heading-copy">
          {detectability && <div className="module-kicker">{subgroup}</div>}
          <h1>{!detectability && subgroup ? `${subgroup} · ${mod.title}` : mod.title}</h1>
          {!mod.live && <div className="meta">Placeholder: shows the page format, nothing computes yet</div>}
        </div>
      </div>

      {folded ? (
        <div className="panel input-summary">
          <div className="input-summary-text">
            {summary.lead && <span className="input-summary-lead">{summary.lead}</span>}
            {summary.items.map((t, i) => <span key={i} className="input-summary-item">{t}</span>)}
          </div>
          <div className="input-summary-actions">
            {run.seconds !== null && <span className="msg">{run.seconds.toFixed(1)} s</span>}
            <button type="button" className="btn secondary" onClick={() => setEditing(true)}>Edit inputs</button>
          </div>
        </div>
      ) : (
        <div className={"panel" + (run.busy ? " is-running" : "")}>
          <div className="fields">{mod.inputs.map((f, i) => <FieldInput key={f.id ?? i} field={f} value={f.id ? values[f.id] : undefined} onChange={v => f.id && setValues({ ...values, [f.id]: v })} />)}</div>
          <div className="actions">
            <button className="btn" onClick={go} disabled={run.busy} aria-busy={run.busy}>{runLabel}</button>
            {mod.example && <button className="btn secondary" type="button" disabled={run.busy}
              onClick={() => setValues({ ...initialValues(mod.inputs), ...mod.example })}>Load example</button>}
            {ok && <button className="btn ghost" type="button" onClick={() => setEditing(false)}>Hide inputs</button>}
            <span className="msg">
              {run.result?.status === "placeholder" && "Not connected yet: nothing computes on this page."}
              {stale && "Inputs changed: run again to update."}
            </span>
          </div>
        </div>
      )}

      {run.busy && <div className="run-progress" role="progressbar" aria-label="Running"><span /></div>}
      {run.error && (
        <div className="alert error">{run.error}
          {run.log && <details><summary>Tool log</summary><pre className="code">{run.log}</pre></details>}
        </div>)}

      {mod.live && !ok && !run.busy && !run.error && (
        <p className="empty-hint">Results appear here. Enter a compound{mod.example ? <>, or <button type="button" className="link-btn" onClick={() => setValues({ ...initialValues(mod.inputs), ...mod.example })}>load the example</button></> : null}, then run.</p>
      )}
      {ok && run.result?.verdict && (
        <div className="verdict" aria-live="polite">
          <div className="verdict-label">Result</div>
          <p className="verdict-text"><RichText text={run.result.verdict} /></p>
        </div>
      )}
      {ok && run.result!.warnings?.map((w, i) => <div className="alert warn" key={i}>{w}</div>)}
      <OutputList blocks={mod.outputs} outputs={ok ? run.result!.outputs : null} live={mod.live} busy={run.busy} />
    </div>
  );
}
