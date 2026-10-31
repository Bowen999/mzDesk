import { useEffect, useRef, useState } from "react";
import { ApiError, getJob, startPipeline } from "../api";
import { AnswerCard } from "../components/AnswerCard";
import { CandidateTable } from "../components/CandidateTable";
import { FieldInput, initialValues, missingRequired, visibleFields } from "../components/Fields";
import { Icon, pipelineIcon } from "../components/Icon";
import { FragmentProvider } from "../components/FragmentDialog";
import { useNetworkGuard } from "../components/NetworkDialog";
import { HelpButton } from "../components/Help";
import { OutputList } from "../components/Outputs";
import { IsomerEvidence, SpectrumEvidence } from "../components/PipelineVisuals";
import { StepList } from "../components/StepList";
import { WorkflowProgress } from "../components/WorkflowProgress";
import { networkUses } from "../network";
import type { AnnotationResult, Field, Job, Spec, WorkflowResult } from "../types";

const POLL_MS = 1500;
const MAX_POLL_FAILURES = 3;                    // a dropped request is retried; an answer of the server (a job that is gone after a restart) ends the run at once
const sleep = (ms: number) => new Promise<void>(resolve => window.setTimeout(resolve, ms));

const ESI_PIPELINE = "esi_ms2_annotation";       // the one workflow with its own result view; the others return the generic WorkflowResult
const isWorkflowResult = (r: AnnotationResult | WorkflowResult): r is WorkflowResult => "pipeline" in r;

/** Fields grouped by `section`, in the order the sections first appear. A section starts closed (behind its heading) when it is "Advanced" or its fields say section_collapsed. */
function groupFields(fields: Field[]): { title: string; fields: Field[]; collapsed: boolean }[] {
  const groups: { title: string; fields: Field[]; collapsed: boolean }[] = [];
  for (const f of fields) {
    const title = f.section ?? "";
    let g = groups.find(x => x.title === title);
    if (!g) { g = { title, fields: [], collapsed: title === "Advanced" }; groups.push(g); }
    g.fields.push(f);
    if (f.section_collapsed) g.collapsed = true;
  }
  return groups;
}

const FRAGMENT_NOTE = "The isomer comparison on this page uses MS-FINDER; this view uses MetFrag, so the explained peaks can differ between the two.";
const RANK_NOTE = "Candidates within 10% of the best base score are ordered by base score x detectability; the others follow in base-score order, so the Final score column does not always go down the list.";
const CCS_RANK_NOTE = "With an observed CCS, a candidate that disagrees with a measured reference goes behind all the others; one that disagrees with a predicted value only goes behind the other near-ties. Such a candidate carries a CCS flag, and the levels do not change.";

/** The generic result (mzdesk/workflows/common.py): the answer (statement, level badge, facts), the warnings, then the headline blocks. */
export function WorkflowResultView({ result }: { result: WorkflowResult }) {
  const a = result.answer;
  const blocks = result.blocks ?? [];
  const facts = a.facts ?? [];
  const outputs = Object.fromEntries(blocks.filter(b => b.id).map(b => [b.id as string, b.data]));
  return (
    <>
      <div className="verdict" aria-live="polite">
        <div className="verdict-label">Answer</div>
        <p className="verdict-text">{a.statement}</p>
      </div>
      {a.level && <div className="pill-row"><span className="badge" title={a.level_label ?? undefined}>{a.level}</span></div>}
      {facts.length > 0 && <dl className="answer-facts">{facts.map(([label, value], i) => <div key={i}><dt>{label}</dt><dd>{value}</dd></div>)}</dl>}
      {(result.warnings ?? []).map((w, i) => <div className="alert warn" key={i}>{w}</div>)}
      {blocks.length > 0 && <OutputList blocks={blocks} outputs={outputs} live />}
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
      setRequestedAdduct(pipe.id === ESI_PIPELINE && adduct && adduct !== "Auto" ? adduct : null);
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
  const jobResult = job?.result ?? null;
  const result = pipe.id === ESI_PIPELINE && jobResult && !isWorkflowResult(jobResult) ? jobResult : null;          // ESI MS/MS Annotation's own result
  const workflowResult = pipe.id !== ESI_PIPELINE && jobResult && isWorkflowResult(jobResult) ? jobResult : null;   // every other workflow
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
          return g.collapsed
            ? <details className="form-section" key={g.title}><summary>{g.title}</summary>{body}</details>
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
          <div className="results-head"><h2 className="panel-title">Ranked candidates<HelpButton lines={result.candidates.some(c => c.ccs) ? [RANK_NOTE, CCS_RANK_NOTE] : [RANK_NOTE]} label="Help: Ranked candidates" /></h2>
            <span className="msg">{result.candidates.length} {result.candidates[0]?.layer === "formula_only" ? "formula(s)" : "structure(s)"}</span></div>
          <FragmentProvider query={result.fragment_query} note={FRAGMENT_NOTE}>
            <CandidateTable candidates={result.candidates} />
          </FragmentProvider>
        </>)}
      {workflowResult && <WorkflowResultView result={workflowResult} />}
      {job && (
        <>
          <div className="results-head"><h2 className="panel-title">Steps</h2></div>
          <StepList steps={job.steps} />
        </>)}
      {guard.dialog}
    </div>
  );
}
