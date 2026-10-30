import { useEffect, useRef, useState } from "react";
import { useParams } from "react-router-dom";
import { ApiError, runModule } from "../api";
import { FieldInput, initialValues, isRequired, missingRequired, summarizeInputs, visibleFields } from "../components/Fields";
import { RichText } from "../components/Help";
import { Icon, moduleIcon } from "../components/Icon";
import { useNetworkGuard } from "../components/NetworkDialog";
import { OutputList } from "../components/Outputs";
import { networkUses } from "../network";
import type { Field, RunResult, Spec } from "../types";

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

type Group = { title: string; fields: Field[] };

/** Group visible fields by `section`, in the order the sections first appear. Fields without a
 *  section (every page except Spectrum) keep the old single flat grid. */
function groupFields(fields: Field[]): Group[] {
  if (!fields.some(f => f.section)) return [{ title: "", fields }];
  const groups: Group[] = [];
  const byTitle = new Map<string, Group>();
  for (const f of fields) {
    const title = f.section ?? "Inputs";
    let g = byTitle.get(title);
    if (!g) { g = { title, fields: [] }; byTitle.set(title, g); groups.push(g); }
    g.fields.push(f);
  }
  return groups;
}

export function ModulePage({ spec }: { spec: Spec }) {
  const { id } = useParams();
  const mod = spec.modules.find(m => m.id === id);
  const [values, setValues] = useState<Record<string, any>>(() => initialValues(mod ? mod.inputs : []));
  const [run, setRun] = useState<Run>(IDLE);
  const [editing, setEditing] = useState(true);     // false = inputs folded into a one-line summary (after a successful run)
  const [errors, setErrors] = useState<Record<string, string>>({});
  const [formError, setFormError] = useState("");
  const [netNote, setNetNote] = useState("");       // "cancelled" note under the Run button
  const guard = useNetworkGuard();                 // asks before an online run while Internet search is off
  const ticket = useRef(0);       // a slow run (ESI MS/MS takes up to minutes) must not land on the page the user has moved to
  const valuesRef = useRef(values);
  valuesRef.current = values;
  const errorRef = useRef<HTMLDivElement | null>(null);
  const elapsed = useElapsed(run.busy);
  useEffect(() => { ticket.current++; if (mod) { setValues(initialValues(mod.inputs)); setRun(IDLE); setEditing(true); setErrors({}); setFormError(""); } }, [id]);
  useEffect(() => { if (run.error) errorRef.current?.focus(); }, [run.error]);
  if (!mod) return <p className="page">Unknown module.</p>;

  const subgroup = spec.subgroups.find(s => s.id === mod.subgroup)?.title;
  const detectability = mod.subgroup === "detectability";
  const gnps = mod.subgroup === "gnps";
  const signature = JSON.stringify(values);
  const ok = run.result?.status === "ok" && !run.busy;
  const stale = (ok || !!run.error) && run.signature !== signature;
  const folded = ok && !editing && !stale;
  const summary = summarizeInputs(mod.inputs, values);
  const shown = visibleFields(mod.inputs, values);
  const groups = groupFields(shown);
  const hasRequired = shown.some(f => isRequired(f, values));
  const runLabel = run.busy ? <><span className="spinner" aria-hidden="true" />Running{elapsed !== null && elapsed >= 1 ? ` · ${Math.floor(elapsed)} s` : "…"}</> : `Run ${mod.title}`;

  const update = (fid: string, v: any) => {
    setValues(prev => ({ ...prev, [fid]: v }));
    if (errors[fid]) setErrors(prev => { const next = { ...prev }; delete next[fid]; return next; });
    if (formError) setFormError("");
    if (netNote) setNetNote("");
  };
  const focusField = (fid: string) => {
    const el = document.querySelector<HTMLElement>(`[data-field="${fid}"]`);
    el?.querySelector<HTMLElement>("input, textarea, select, button")?.focus();
  };
  const loadExample = () => {
    if (!mod.example) return;
    setValues({ ...initialValues(mod.inputs), ...mod.example });
    setRun(IDLE);
    setEditing(true);
    setErrors({});
    setFormError("");
    setNetNote("");
  };
  const go = async () => {
    const missing = missingRequired(mod.inputs, values);
    if (missing.length) {
      const next: Record<string, string> = {};
      for (const f of missing) next[f.id!] = f.type === "structure" ? "Enter or draw a structure." : "Required.";
      setErrors(next);
      setFormError(`Complete ${missing.length === 1 ? "this required input" : "these required inputs"}: ${missing.map(f => f.label).join(", ")}.`);
      focusField(missing[0].id!);
      return;
    }
    setErrors({});
    setFormError("");
    // Internet search off: online runs ask first. Hidden fields are not sent, so they do not trigger the question.
    if (!(await guard.confirm(networkUses(shown, values, mod.network)))) {
      setNetNote("Cancelled: internet search is off.");
      return;
    }
    setNetNote("");
    const t0 = performance.now(), mine = ++ticket.current, started = signature;
    // Hidden controls (for example MS/MS-only fields in EI mode) are neither sent nor validated.
    const payload = Object.fromEntries(shown.filter(f => f.id).map(f => [f.id!, values[f.id!]]));
    setRun(r => ({ ...IDLE, busy: true, result: r.result, signature: r.signature }));
    try {
      const result = await runModule(mod.id, payload);
      if (ticket.current !== mine) return;
      setRun({ ...IDLE, result, seconds: (performance.now() - t0) / 1000, signature: started });
      // Only fold when the inputs did not change while the run was in flight.
      if (result.status === "ok" && JSON.stringify(valuesRef.current) === started) setEditing(false);
    } catch (e: any) {
      if (ticket.current === mine) { setRun({ ...IDLE, error: e.message ?? String(e), log: e instanceof ApiError ? e.log : "", signature: started }); setEditing(true); }
    }
  };

  const field = (f: Field, i: number) => (
    <div className="field-slot" data-field={f.id} key={f.id ?? i}>
      <FieldInput field={f} value={f.id ? values[f.id] : undefined} error={f.id ? errors[f.id] : undefined}
        onChange={v => f.id && update(f.id, v)} />
    </div>
  );

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
      <p className="one-liner"><RichText text={mod.one_liner} /></p>

      {folded ? (
        <div className="panel input-summary">
          <div className="input-summary-text">
            {summary.lead && <span className="input-summary-lead">{summary.lead}</span>}
            {summary.items.map((t, i) => <span key={i} className="input-summary-item">{t}</span>)}
          </div>
          <div className="input-summary-actions">
            {netNote && <span className="msg">{netNote}</span>}
            {run.seconds !== null && <span className="msg">{run.seconds.toFixed(1)} s</span>}
            <button type="button" className="btn secondary" onClick={() => setEditing(true)}>Edit inputs</button>
            <button type="button" className="btn ghost" onClick={go} disabled={run.busy}>Run again</button>
          </div>
        </div>
      ) : (
        <div className={"panel" + (run.busy ? " is-running" : "")}>
          <div className="panel-head">
            <h2 className="panel-title">Inputs</h2>
            <span className="panel-note">{hasRequired ? <><span className="req" aria-hidden="true">*</span> required</> : "all inputs optional"}</span>
          </div>
          {formError && <div className="alert error" role="alert">{formError}</div>}
          {groups.map((g, gi) => {
            // A run of checkboxes at the end of a section reads better as one row; a checkbox
            // in the middle of number/select fields stays next to them.
            const firstCheck = g.fields.findIndex(f => f.type === "checkbox");
            const checkRun = firstCheck >= 0 && g.fields.length - firstCheck >= 2 && g.fields.slice(firstCheck).every(f => f.type === "checkbox");
            const body = checkRun ? <>
              {firstCheck > 0 && <div className="fields">{g.fields.slice(0, firstCheck).map(field)}</div>}
              <div className="check-row">{g.fields.slice(firstCheck).map(field)}</div>
            </> : <div className="fields">{g.fields.map(field)}</div>;
            if (gnps && g.title === "Optional verification") return (
              <details className="gnps-advanced" key={g.title}>
                <summary><span>{g.title}</span><small>Diagnostic ions and verification settings</small></summary>
                <fieldset className="form-section"><legend className="sr-only">{g.title}</legend>{body}</fieldset>
              </details>);
            return g.title
              ? <fieldset className="form-section" key={g.title}><legend>{g.title}</legend>{body}</fieldset>
              : <div className="form-plain" key={gi}>{body}</div>;
          })}
          <div className="actions">
            <button className="btn" onClick={go} disabled={run.busy} aria-busy={run.busy}>{runLabel}</button>
            {mod.example && <button className="btn secondary" type="button" disabled={run.busy} onClick={loadExample}>Load example</button>}
            {ok && !stale && <button className="btn ghost" type="button" onClick={() => setEditing(false)}>Hide inputs</button>}
            <span className="msg" aria-live="polite">
              {netNote || (run.result?.status === "placeholder" && "Not connected yet: nothing computes on this page.")
                || (stale && "Inputs changed since this result. Run again to update it.")}
            </span>
          </div>
        </div>
      )}

      {run.busy && <div className="run-progress" role="progressbar" aria-label="Running"><span /></div>}
      {run.error && (
        <div className="alert error" role="alert" tabIndex={-1} ref={errorRef}>{run.error}
          {run.log && <details><summary>Tool log</summary><pre className="code">{run.log}</pre></details>}
        </div>)}

      {mod.live && !ok && !run.busy && !run.error && (
        <p className="empty-hint">Results appear here after you run this module.{mod.example ? <> Or <button type="button" className="link-btn" onClick={loadExample}>load the example</button>.</> : null}</p>
      )}
      {ok && (
        <div className="results-head">
          <h2 className="panel-title">Results</h2>
          {run.seconds !== null && <span className="msg">{run.seconds.toFixed(1)} s</span>}
          {stale && <span className="msg warn-text">Inputs changed since this run.</span>}
        </div>)}
      {ok && run.result?.verdict && (
        <div className="verdict" aria-live="polite">
          <div className="verdict-label">Result</div>
          <p className="verdict-text"><RichText text={run.result.verdict} /></p>
        </div>
      )}
      {ok && stale && <div className="alert warn">Inputs changed since this run. The results below belong to the previous inputs; run again to update them.</div>}
      {ok && gnps && !!run.result!.warnings?.length && <div className="alert warn gnps-notices" role="note">
        <strong>Notes</strong><ul>{run.result!.warnings.map((w, i) => <li key={i}>{w}</li>)}</ul>
      </div>}
      {ok && !gnps && run.result!.warnings?.map((w, i) => <div className="alert warn" key={i}>{w}</div>)}
      <OutputList blocks={mod.outputs} outputs={ok ? run.result!.outputs : null} live={mod.live} busy={run.busy} />
      {guard.dialog}
    </div>
  );
}
