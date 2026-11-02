import { useEffect, useRef, useState } from "react";
import { ApiError, runModule } from "../api";
import { FieldInput, initialValues, missingRequired, summarizeInputs, visibleFields } from "../components/Fields";
import { RichText } from "../components/Help";
import { moduleIcon } from "../components/Icon";
import { useNetworkGuard } from "../components/NetworkDialog";
import { OutputList, coveredCards } from "../components/Outputs";
import { PageFoot, PageHead } from "../components/PageHead";
import { useT } from "../i18n";
import { RunBar } from "../components/RunBar";
import { networkUses } from "../network";
import { dropRun, recordRun } from "../report/log";
import type { Field, Module, RunResult, Spec } from "../types";

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

export const modulePageClass = (mod: Module) => "page module-page" + (mod.subgroup === "detectability" ? " detectability-page" : "");

/** The title of a module page and its "have → get" line. The Session report draws the same. */
export function ModuleHead({ mod }: { mod: Module; spec?: Spec }) {
  return <PageHead icon={moduleIcon(mod)} title={mod.title} flow={mod.one_liner} preview={!mod.live} />;
}

/** Everything under the inputs once a run has a result (or is running): the verdict, the warnings and the output blocks.
 *  The Session report draws it with ok = true. */
export function ModuleResults({ mod, result, ok, busy, stale }:
  { mod: Module; result: RunResult | null; ok: boolean; busy: boolean; seconds?: number | null; stale: boolean }) {
  const t = useT();
  const gnps = mod.subgroup === "gnps";
  return (
    <>
      {ok && stale && <div className="alert warn">{t("Inputs changed. Run again to update the results.")}</div>}
      {ok && result?.verdict && (
        <div className="verdict" aria-live="polite">
          <p className="verdict-text"><RichText text={result.verdict} /></p>
        </div>
      )}
      {ok && gnps && !!result!.warnings?.length && <div className="alert warn gnps-notices" role="note">
        <ul>{result!.warnings.map((w, i) => <li key={i}>{w}</li>)}</ul>
      </div>}
      {ok && !gnps && result!.warnings?.map((w, i) => <div className="alert warn" key={i}>{w}</div>)}
      <OutputList blocks={mod.outputs} outputs={ok ? result!.outputs : null} hide={ok ? coveredCards(mod.outputs, result!.outputs, result!.verdict) : []} live={mod.live} busy={busy} />
    </>
  );
}

export function ModulePage({ id, spec, active }: { id: string; spec: Spec; active: boolean }) {
  const t = useT();
  const mod = spec.modules.find(m => m.id === id);
  const [values, setValues] = useState<Record<string, any>>(() => initialValues(mod ? mod.inputs : []));
  const [run, setRun] = useState<Run>(IDLE);
  const [editing, setEditing] = useState(true);     // false = inputs folded into a one-line summary (after a successful run)
  const [errors, setErrors] = useState<Record<string, string>>({});
  const [formError, setFormError] = useState("");
  const [netNote, setNetNote] = useState("");       // "cancelled" note under the Run button
  const [formVersion, setFormVersion] = useState(0); // reset uncontrolled file inputs as well
  const guard = useNetworkGuard();                 // asks before an online run while Internet search is off
  const ticket = useRef(0);       // Clear or a newer run invalidates an older response
  const confirming = useRef(false);
  const valuesRef = useRef(values);
  valuesRef.current = values;
  const errorRef = useRef<HTMLDivElement | null>(null);
  const pageRef = useRef<HTMLDivElement | null>(null);
  const elapsed = useElapsed(run.busy);
  useEffect(() => { ticket.current++; dropRun("module", id); if (mod) { setValues(initialValues(mod.inputs)); setRun(IDLE); setEditing(true); setErrors({}); setFormError(""); } }, [id]);
  useEffect(() => { if (!active) guard.cancel(); }, [active, guard.cancel]);
  useEffect(() => { if (run.error) errorRef.current?.focus(); }, [run.error]);
  if (!mod) return <p className="page">{t("Unknown module.")}</p>;

  const subgroup = spec.subgroups.find(s => s.id === mod.subgroup)?.title;
  const gnps = mod.subgroup === "gnps";
  const signature = JSON.stringify(values);
  const ok = run.result?.status === "ok" && !run.busy;
  const stale = (ok || !!run.error) && run.signature !== signature;
  const folded = ok && !editing && !stale;
  const summary = summarizeInputs(mod.inputs, values);
  const shown = visibleFields(mod.inputs, values);
  const groups = groupFields(shown);
  const runLabel = run.busy ? <><span className="spinner" aria-hidden="true" />{elapsed !== null && elapsed >= 1 ? t("Running · {n} s", { n: Math.floor(elapsed) }) : t("Running…")}</> : t("Run");

  const update = (fid: string, v: any) => {
    setValues(prev => ({ ...prev, [fid]: v }));
    if (errors[fid]) setErrors(prev => { const next = { ...prev }; delete next[fid]; return next; });
    if (formError) setFormError("");
    if (netNote) setNetNote("");
  };
  const focusField = (fid: string) => {
    const el = pageRef.current?.querySelector<HTMLElement>(`[data-field="${fid}"]`);
    el?.querySelector<HTMLElement>("input, textarea, select, button")?.focus();
  };
  const loadExample = () => {
    if (!mod.example) return;
    ticket.current++; confirming.current = false; guard.cancel(); dropRun("module", mod.id);
    setValues({ ...initialValues(mod.inputs), ...mod.example });
    setRun(IDLE);
    setEditing(true);
    setErrors({});
    setFormError("");
    setNetNote("");
  };
  const clear = () => {
    ticket.current++;
    confirming.current = false;
    guard.cancel();
    dropRun("module", mod.id);
    setValues(initialValues(mod.inputs));
    setRun(IDLE); setEditing(true); setErrors({}); setFormError("");
    setNetNote(run.busy ? t("Cleared. The previous run may still finish on the server.") : "");
    setFormVersion(n => n + 1);
  };
  const go = async () => {
    if (run.busy || confirming.current) return;
    const missing = missingRequired(mod.inputs, values);
    if (missing.length) {
      const next: Record<string, string> = {};
      for (const f of missing) next[f.id!] = f.type === "structure" ? t("Enter or draw a structure.") : t("Required.");
      setErrors(next);
      const labels = missing.map(f => f.label).join(", ");
      setFormError(missing.length === 1 ? t("Complete this required input: {labels}.", { labels }) : t("Complete these required inputs: {labels}.", { labels }));
      focusField(missing[0].id!);
      return;
    }
    setErrors({});
    setFormError("");
    const mine = ++ticket.current;
    // Internet search off: online runs ask first. Hidden fields are not sent, so they do not trigger the question.
    confirming.current = true;
    const allowed = await guard.confirm(networkUses(shown, values, mod.network));
    if (ticket.current === mine) confirming.current = false;
    if (!allowed) {
      if (ticket.current === mine) setNetNote(t("Cancelled: internet search is off."));
      return;
    }
    if (ticket.current !== mine) return;
    setNetNote("");
    const t0 = performance.now(), started = signature;
    // Hidden controls (for example MS/MS-only fields in EI mode) are neither sent nor validated.
    const payload = Object.fromEntries(shown.filter(f => f.id).map(f => [f.id!, values[f.id!]]));
    setRun(r => ({ ...IDLE, busy: true, result: r.result, signature: r.signature }));
    try {
      const result = await runModule(mod.id, payload);
      if (ticket.current !== mine) return;
      const seconds = (performance.now() - t0) / 1000;
      setRun({ ...IDLE, result, seconds, signature: started });
      // The session report lists the pages that show a result right now.
      if (result.status === "ok") recordRun({ kind: "module", id: mod.id, title: mod.title, section: subgroup ?? "Modules", inputs: payload, result, seconds });
      else dropRun("module", mod.id);
      // Only fold when the inputs did not change while the run was in flight.
      if (result.status === "ok" && JSON.stringify(valuesRef.current) === started) setEditing(false);
    } catch (e: any) {
      if (ticket.current === mine) { dropRun("module", mod.id); setRun({ ...IDLE, error: e.message ?? String(e), log: e instanceof ApiError ? e.log : "", signature: started }); setEditing(true); }
    }
  };

  const field = (f: Field, i: number) => (
    <div className="field-slot" data-field={f.id} key={f.id ?? i}>
      <FieldInput field={f} value={f.id ? values[f.id] : undefined} error={f.id ? errors[f.id] : undefined}
        onChange={v => f.id && update(f.id, v)} />
    </div>
  );

  return (
    <div ref={pageRef} className={modulePageClass(mod)} data-module={mod.id}>
      <ModuleHead mod={mod} />

      {folded ? (
        <div className="panel input-summary">
          <div className="input-summary-text">
            {summary.lead && <span className="input-summary-lead">{summary.lead}</span>}
            {summary.items.map((t, i) => <span key={i} className="input-summary-item">{t}</span>)}
          </div>
          <div className="input-summary-actions">
            {netNote && <span className="msg">{netNote}</span>}
            {run.seconds !== null && <span className="msg">{run.seconds.toFixed(1)} s</span>}
            <button type="button" className="btn secondary" onClick={() => setEditing(true)}>{t("Edit inputs")}</button>
            <button type="button" className="btn ghost" onClick={go} disabled={run.busy}>{t("Run again")}</button>
            <button type="button" className="btn ghost" onClick={clear}>{t("Clear")}</button>
          </div>
        </div>
      ) : (
        <div key={formVersion} className={"panel" + (run.busy ? " is-running" : "")}>
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
            // a section title stays English in the spec (localizeSpec "lookup"): the dictionary translates it where it is drawn, t(g.title)
            if (gnps && g.title === "Optional verification") return (
              <details className="gnps-advanced" key={g.title}>
                <summary>{t(g.title)}</summary>
                <fieldset className="form-section"><legend className="sr-only">{t(g.title)}</legend>{body}</fieldset>
              </details>);
            return <div className="form-section" key={g.title || gi}>{body}</div>;
          })}
          <div className="actions">
            <button className="btn" onClick={go} disabled={run.busy} aria-busy={run.busy}>{runLabel}</button>
            <button className="btn ghost" type="button" onClick={clear}>{t("Clear")}</button>
            {mod.example && <button className="btn secondary" type="button" disabled={run.busy} onClick={loadExample}>{t("Load example")}</button>}
            {ok && !stale && <button className="btn ghost" type="button" onClick={() => setEditing(false)}>{t("Hide inputs")}</button>}
            <span className="msg" aria-live="polite">
              {netNote || (run.result?.status === "placeholder" && t("Not connected yet: nothing computes on this page."))}
            </span>
          </div>
        </div>
      )}

      <RunBar busy={run.busy} />
      {run.error && (
        <div className="alert error" role="alert" tabIndex={-1} ref={errorRef}>{run.error}
          {run.log && <details><summary>{t("Tool log")}</summary><pre className="code">{run.log}</pre></details>}
        </div>)}

      <ModuleResults mod={mod} result={run.result} ok={ok} busy={run.busy} stale={stale} />
      <PageFoot id={mod.id} basedOn={mod.based_on} />
      {guard.dialog}
    </div>
  );
}
