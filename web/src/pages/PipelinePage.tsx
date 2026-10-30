import { useEffect, useRef, useState } from "react";
import { ApiError, getJob, startPipeline } from "../api";
import { RowsTable } from "../components/Blocks";
import { FieldInput, initialValues, missingRequired, visibleFields } from "../components/Fields";
import { Icon, pipelineIcon } from "../components/Icon";
import { FragmentProvider } from "../components/FragmentDialog";
import { useNetworkGuard } from "../components/NetworkDialog";
import { OutputList } from "../components/Outputs";
import { IsomerEvidence, SpectrumEvidence } from "../components/PipelineVisuals";
import { WorkflowProgress } from "../components/WorkflowProgress";
import { networkUses } from "../network";
import type { AnnotationCandidate, AnnotationResult, Cell, Field, Job, JobStep, Row, Spec } from "../types";

const POLL_MS = 1500;
const MAX_POLL_FAILURES = 3;                    // a dropped request is retried; an answer of the server (a job that is gone after a restart) ends the run at once
const sleep = (ms: number) => new Promise<void>(resolve => window.setTimeout(resolve, ms));

const STATUS_TONE: Record<string, string> = { done: "good", running: "jp", failed: "bad", skipped: "neutral", queued: "neutral" };
const LEVEL_TONE: Record<string, string> = { "2a": "good", "3": "jp", H: "warn", "4": "neutral", "5": "neutral" };
const levelName = (label: string) => label.split(":")[0];

/** Fields grouped by `section`, in the order the sections first appear. */
function groupFields(fields: Field[]): { title: string; fields: Field[] }[] {
  const groups: { title: string; fields: Field[] }[] = [];
  for (const f of fields) {
    const title = f.section ?? "";
    let g = groups.find(x => x.title === title);
    if (!g) { g = { title, fields: [] }; groups.push(g); }
    g.fields.push(f);
  }
  return groups;
}

function StatusChip({ status }: { status: string }) {
  return <span className={"badge tone-" + (STATUS_TONE[status] ?? "neutral")}>{status}</span>;
}

/** One step: status chip, seconds, note, and the blocks it has produced so far (tables with ids, filled as soon as the step ends). */
function StepView({ step, index }: { step: JobStep; index: number }) {
  const blocks = step.blocks ?? [];
  const outputs = Object.fromEntries(blocks.filter(b => b.id).map(b => [b.id as string, b.data]));
  const shown = step.status !== "queued";
  const hasSeconds = typeof step.seconds === "number" && Number.isFinite(step.seconds) && step.seconds >= 0;
  return (
    <details id={`workflow-step-${step.id}`} className="stage" open={shown && blocks.length > 0 && step.status !== "skipped"}>
      <summary>
        {index + 1}. {step.title || "Unnamed stage"} <StatusChip status={step.status || "unknown"} />
        {step.status !== "queued" && step.status !== "skipped" && hasSeconds && <span className="meta"> {step.seconds.toFixed(1)} s</span>}
        {step.note && <span className="meta"> {step.note}</span>}
      </summary>
      {blocks.length > 0 && <OutputList blocks={blocks} outputs={outputs} live />}
    </details>
  );
}

const detectText = (c: AnnotationCandidate) => {
  const d = c.detect;
  return d ? `${d.omics ?? "—"} / ${d.ionizable ?? "—"} / ${d.adduct ?? "—"}` : "—";
};

function isomerCell(c: AnnotationCandidate): Cell {
  const iso = c.isomer;
  if (!iso || (!iso.role && !iso.challenger)) return "—";
  const text = iso.role ?? "not compared";
  return iso.challenger ? { badge: `${text} · challenger`, tone: "warn" } : text;
}

/** "View": the fragment view of this candidate (computed when clicked); a dash when it has no structure or MetFrag cannot use its adduct. */
function fragmentCell(c: AnnotationCandidate): Cell {
  return c.smiles && c.fragment_adduct ? { fragment: { smiles: c.smiles, adduct: c.fragment_adduct, name: c.name || undefined } } : "—";
}

const FRAGMENT_NOTE = "The isomer comparison on this page uses MS-FINDER; this view uses MetFrag, so the explained peaks can differ between the two.";

function rankedRows(candidates: AnnotationCandidate[]): Row[] {
  return candidates.map(c => [
    c.rank ?? "—", fragmentCell(c), c.name || "—", c.smiles || "—", c.formula || "—", c.adduct || "—", c.layer_label,
    c.level ? { badge: levelName(c.level_label), tone: LEVEL_TONE[c.level] ?? "neutral", help: c.level_reason } : "—",
    c.final_score === null ? "—" : c.final_score.toFixed(3), detectText(c), isomerCell(c),
  ]);
}

/** The answer: statement, layer, level, the group that cannot be separated, and the isomers that explain peaks the top-1 cannot. */
function AnswerCard({ result, requestedAdduct }: { result: AnnotationResult; requestedAdduct: string | null }) {
  const a = result.answer;
  const byId = new Map(result.candidates.map(c => [c.id, c]));
  const label = (id: string) => byId.get(id)?.name || byId.get(id)?.smiles || id;
  const challengers = result.isomer?.challengers ?? [];
  const conflicting = a.isomer_verdict === "single" && challengers.length > 0;
  const top = a.top_id ? byId.get(a.top_id) : undefined;
  const otherAdducts = requestedAdduct ? [...new Set(result.candidates.filter(c => c.evidence?.library?.gate_passed && c.adduct && c.adduct !== requestedAdduct)
    .map(c => c.adduct))] : [];
  return (
    <>
      <div className="verdict" aria-live="polite">
        <div className="verdict-label">Answer</div>
        <p className="verdict-text">{a.statement}</p>
      </div>
      <div className="pill-row">
        {a.layer_label && <span className="badge">Layer: {a.layer_label}</span>}
        {a.level && <span className={"badge tone-" + (LEVEL_TONE[a.level] ?? "neutral")} title={a.level_label}>{a.level_label ? levelName(a.level_label) : a.level}</span>}
        {a.isomer_verdict && a.isomer_verdict !== "not run" && <span className={"badge" + (conflicting ? " tone-warn" : "")}>Isomer check: {conflicting ? "conflicting evidence" : a.isomer_verdict}</span>}
        {result.context.gate_passed && <span className="badge tone-good">Library gate passed</span>}
      </div>
      {otherAdducts.length > 0 && <div className="alert warn" role="status">Input adduct: {requestedAdduct}. A library hit passed with {otherAdducts.join(", ")}; review the library checks.</div>}
      {requestedAdduct && top?.adduct && top.adduct !== requestedAdduct && !otherAdducts.includes(top.adduct) &&
        <div className="alert warn" role="status">Input adduct: {requestedAdduct}. The top candidate uses {top.adduct}; review the adduct evidence.</div>}
      {conflicting && <div className="alert warn" role="status">Fragment explanations conflict. Treat the single-answer verdict as unresolved.</div>}
      {a.type === "group" && (
        <div className="alert warn" role="status">These cannot be separated by the spectrum: {a.ids.map(label).join(", ")}. The ranking only orders them.</div>
      )}
      {challengers.map(ch => (
        <div className="alert warn" role="status" key={ch.id}>
          Another isomer explains peaks the top-1 cannot: {ch.name || label(ch.id)} (m/z {ch.peaks.map(p => p.mz.toFixed(4)).join(", ")}).
          The verdict above is unchanged; treat the top-1 with care.
        </div>))}
      {result.warnings.map((w, i) => <div className="alert warn" key={i}>{w}</div>)}
    </>
  );
}

export function PipelinePage({ id, spec, active }: { id: string; spec: Spec; active: boolean }) {
  const pipe = spec.pipelines.find(p => p.id === id);
  const [values, setValues] = useState<Record<string, any>>(() => initialValues(pipe ? pipe.inputs : []));
  const [job, setJob] = useState<Job | null>(null);
  const [busy, setBusy] = useState(false);
  const [hasStartedJob, setHasStartedJob] = useState(false);
  const [error, setError] = useState("");
  const [errors, setErrors] = useState<Record<string, string>>({});
  const [requestedAdduct, setRequestedAdduct] = useState<string | null>(null);
  const [netNote, setNetNote] = useState("");     // "cancelled" note under the Run button
  const [formVersion, setFormVersion] = useState(0);
  const guard = useNetworkGuard();               // asks before an online run while Internet search is off
  const ticket = useRef(0);                      // Clear, a newer run or unmount ends polling of an older job
  const confirming = useRef(false);
  useEffect(() => {
    ticket.current++;
    setValues(initialValues(pipe ? pipe.inputs : [])); setJob(null); setBusy(false); setHasStartedJob(false); setError(""); setErrors({}); setRequestedAdduct(null);
    return () => { ticket.current++; };
  }, [id]);
  useEffect(() => { if (!active) guard.cancel(); }, [active, guard.cancel]);
  if (!pipe) return <p className="page">Unknown pipeline.</p>;

  const shown = visibleFields(pipe.inputs, values);
  const update = (fid: string, v: any) => {
    setValues(prev => ({ ...prev, [fid]: v }));
    if (errors[fid]) setErrors(prev => { const next = { ...prev }; delete next[fid]; return next; });
    if (netNote) setNetNote("");
  };
  const loadExample = (example: Record<string, any>) => {
    ticket.current++; confirming.current = false; guard.cancel();
    setValues({ ...initialValues(pipe.inputs), ...example }); setJob(null); setBusy(false); setHasStartedJob(false); setError(""); setErrors({}); setRequestedAdduct(null); setNetNote("");
  };
  const clear = () => {
    ticket.current++;
    confirming.current = false;
    guard.cancel();
    setValues(initialValues(pipe.inputs)); setJob(null); setBusy(false); setHasStartedJob(false); setError("");
    setErrors({}); setRequestedAdduct(null);
    setNetNote(busy ? "Cleared. The previous workflow may still finish on the server." : "");
    setFormVersion(n => n + 1);
  };

  const run = async () => {
    if (busy || confirming.current) return;
    const missing = missingRequired(pipe.inputs, values);
    if (missing.length) {
      setErrors(Object.fromEntries(missing.map(f => [f.id!, "Required."])));
      setError(`Complete ${missing.length === 1 ? "this required input" : "these required inputs"}: ${missing.map(f => f.label).join(", ")}.`);
      return;
    }
    const mine = ++ticket.current;
    // Internet search off: a workflow that would fetch a USI asks first.
    confirming.current = true;
    const allowed = await guard.confirm(networkUses(shown, values));
    if (ticket.current === mine) confirming.current = false;
    if (!allowed) {
      if (ticket.current === mine) setNetNote("Cancelled: internet search is off.");
      return;
    }
    if (ticket.current !== mine) return;
    setNetNote("");
    setBusy(true); setHasStartedJob(false); setError(""); setErrors({}); setJob(null);
    try {
      const payload = Object.fromEntries(shown.filter(f => f.id).map(f => [f.id!, values[f.id!]]));
      const adduct = payload.precursor?.adduct;
      setRequestedAdduct(adduct && adduct !== "Auto" ? adduct : null);
      const { job_id } = await startPipeline(pipe.id, payload);
      if (ticket.current !== mine) return;
      setHasStartedJob(true);
      let failures = 0;
      for (;;) {
        let j: Job;
        try { j = await getJob(job_id); failures = 0; }
        catch (e: any) {
          if (ticket.current !== mine) return;
          if (e instanceof ApiError || ++failures >= MAX_POLL_FAILURES) throw e;
          await sleep(POLL_MS);
          continue;
        }
        if (ticket.current !== mine) return;
        setJob(j);
        if (j.status !== "running" && j.status !== "queued") { if (j.error) setError(j.error); break; }
        await sleep(POLL_MS);
        if (ticket.current !== mine) return;
      }
    } catch (e: any) {
      if (ticket.current === mine) setError(e?.message || String(e) || "An unknown error occurred.");
    } finally {
      if (ticket.current === mine) setBusy(false);
    }
  };

  const field = (f: Field, i: number) => (
    <div className="field-slot" data-field={f.id} key={f.id ?? i}>
      <FieldInput field={f} value={f.id ? values[f.id] : undefined} error={f.id ? errors[f.id] : undefined} onChange={v => f.id && update(f.id, v)} />
    </div>);
  const result = job?.result ?? null;
  return (
    <div className="page module-page" data-pipeline={pipe.id}>
      <div className="module-heading">
        <span className="module-emblem"><Icon name={pipelineIcon(pipe.id)} /></span>
        <div className="module-heading-copy"><div className="module-kicker">Workflows</div><h1>{pipe.title}</h1></div>
      </div>
      <p className="one-liner" dangerouslySetInnerHTML={{ __html: pipe.one_liner }} />
      <div key={formVersion} className={"panel" + (busy ? " is-running" : "")}>
        <div className="panel-head"><h2 className="panel-title">Inputs</h2><span className="panel-note"><span className="req" aria-hidden="true">*</span> required</span></div>
        {groupFields(shown).map(g => {
          const body = g.fields.every(f => f.type === "checkbox")
            ? <div className="check-row">{g.fields.map(field)}</div>
            : <div className="fields">{g.fields.map(field)}</div>;
          return g.title === "Advanced"
            ? <details className="form-section" key={g.title}><summary>Advanced (mobile phase context for the ranking)</summary>{body}</details>
            : g.title ? <fieldset className="form-section" key={g.title}><legend>{g.title}</legend>{body}</fieldset> : <div key="plain">{body}</div>;
        })}
        <div className="actions">
          <button className="btn" onClick={run} disabled={busy} aria-busy={busy}>{busy ? <><span className="spinner" aria-hidden="true" />Running…</> : "Run workflow"}</button>
          <button className="btn ghost" type="button" onClick={clear}>Clear</button>
          {(pipe.examples ?? []).map(ex => (
            <button className="btn secondary" type="button" key={ex.label} disabled={busy} onClick={() => loadExample(ex.values)}>Load example: {ex.label}</button>))}
          <span className="msg" aria-live="polite">{netNote}</span>
        </div>
      </div>
      {error && <div className="alert error" role="alert">{error}</div>}
      {(busy || job || (hasStartedJob && error)) && <WorkflowProgress job={job} awaitingFirstUpdate={hasStartedJob}
        trackingStopped={!!error && !busy && hasStartedJob && (!job || job.status === "running" || job.status === "queued")} />}
      {result && (
        <>
          <AnswerCard result={result} requestedAdduct={requestedAdduct} />
          <SpectrumEvidence result={result} />
          <IsomerEvidence result={result} />
          <div className="results-head"><h2 className="panel-title">Ranked candidates</h2><span className="msg">{result.candidates.length} {result.candidates[0]?.layer === "formula_only" ? "formula(s)" : "structure(s)"}</span></div>
          <FragmentProvider query={result.fragment_query} note={FRAGMENT_NOTE}>
            <RowsTable columns={["Rank", "Fragments", "Name", "SMILES", "Formula", "Adduct", "Layer", "Level", "Final score", "Detectability (omics / ionizable / adduct)", "Isomer role"]}
              rows={rankedRows(result.candidates)} empty="no candidates" />
          </FragmentProvider>
          <div className="meta">Candidates within 10% of the best base score are ordered by base score x detectability; the others follow in base-score order, so the Final score column does not always go down the list.</div>
        </>)}
      {job && (
        <>
          <div className="results-head"><h2 className="panel-title">Steps</h2></div>
          {(Array.isArray(job.steps) ? job.steps : []).map((s, i) => <StepView step={s} index={i} key={s.id} />)}
        </>)}
      {guard.dialog}
    </div>
  );
}
