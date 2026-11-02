import { useEffect, useRef, useState } from "react";
import { ApiError, getJob, startPipeline } from "../api";
import { AnswerCard } from "../components/AnswerCard";
import { CandidateTable } from "../components/CandidateTable";
import { FieldInput, initialValues, missingRequired, visibleFields } from "../components/Fields";
import { pipelineIcon } from "../components/Icon";
import { FragmentProvider } from "../components/FragmentDialog";
import { useNetworkGuard } from "../components/NetworkDialog";
import { HelpButton } from "../components/Help";
import { OutputList } from "../components/Outputs";
import { PageFoot, PageHead } from "../components/PageHead";
import { localizeRunBlocks, useT } from "../i18n";
import { IsomerEvidence, SpectrumEvidence } from "../components/PipelineVisuals";
import { StepList } from "../components/StepList";
import { StructureProvider } from "../components/StructureHover";
import { WorkflowProgress } from "../components/WorkflowProgress";
import { networkUses } from "../network";
import { dropRun, recordRun } from "../report/log";
import { ensureSetup } from "../setup";
import type { AnnotationResult, Field, Job, JobStep, Pipeline, Spec, WorkflowResult } from "../types";

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

/** The generic result (mzdesk/workflows/common.py): the answer (statement, level badge, facts), the warnings, then the headline blocks. */
export function WorkflowResultView({ result }: { result: WorkflowResult }) {
  useT();           // redraw when the language changes
  const a = result.answer;
  const blocks = localizeRunBlocks(result.blocks ?? []);
  const facts = a.facts ?? [];
  const outputs = Object.fromEntries(blocks.filter(b => b.id).map(b => [b.id as string, b.data]));
  return (
    <>
      <div className="verdict" aria-live="polite">
        <p className="verdict-text">{a.statement}</p>
      </div>
      {a.level && <div className="pill-row"><span className="badge" title={a.level_label ?? undefined}>{a.level}</span></div>}
      {facts.length > 0 && <dl className="answer-facts">{facts.map(([label, value], i) => <div key={i}><dt>{label}</dt><dd>{value}</dd></div>)}</dl>}
      {(result.warnings ?? []).map((w, i) => <div className="alert warn" key={i}>{w}</div>)}
      {blocks.length > 0 && <OutputList blocks={blocks} outputs={outputs} live />}
    </>
  );
}

/** The adduct the person asked for (not "Auto"): ESI MS/MS Annotation says when its answer uses another one. */
export const requestedAdductOf = (pipeId: string, payload: Record<string, any>): string | null => {
  const adduct = payload.precursor?.adduct;
  return pipeId === ESI_PIPELINE && adduct && adduct !== "Auto" ? adduct : null;
};

/** The title of a workflow page and its "have → get" line. The Session report draws the same. */
export function PipelineHead({ pipe }: { pipe: Pipeline }) {
  return <PageHead icon={pipelineIcon(pipe.id)} title={pipe.title} flow={pipe.one_liner} />;
}

/** What a workflow page shows under its inputs and progress bar: the result (ESI MS/MS Annotation's own view, or the generic one) and the steps.
 *  The Session report draws the same from the logged result and steps. */
export function PipelineResults({ pipeId, jobResult, job, requestedAdduct }:
  { pipeId: string; jobResult: AnnotationResult | WorkflowResult | null; job: { steps: JobStep[] } | null; requestedAdduct: string | null }) {
  const t = useT();
  const rankNote = t("Candidates within 10% of the best score are ordered by score × detectability.");
  const ccsRankNote = t("A CCS that disagrees with a measured value moves a candidate to the end; with a predicted value, behind its near-ties.");
  const result = pipeId === ESI_PIPELINE && jobResult && !isWorkflowResult(jobResult) ? jobResult : null;          // ESI MS/MS Annotation's own result
  const workflowResult = pipeId !== ESI_PIPELINE && jobResult && isWorkflowResult(jobResult) ? jobResult : null;   // every other workflow
  const content = (
    <>
      {result && (
        <>
          <AnswerCard result={result} requestedAdduct={requestedAdduct} />
          <SpectrumEvidence result={result} />
          <IsomerEvidence result={result} />
          <div className="results-head"><h2 className="panel-title">{t("Ranked candidates")}<HelpButton lines={result.candidates.some(c => c.ccs) ? [rankNote, ccsRankNote] : [rankNote]} label={t("Help: Ranked candidates")} /></h2>
            <span className="msg">{result.candidates[0]?.layer === "formula_only" ? t("{n} formula(s)", { n: result.candidates.length }) : t("{n} structure(s)", { n: result.candidates.length })}</span></div>
          <FragmentProvider query={result.fragment_query} note={t("Isomer comparison uses MS-FINDER, this view uses MetFrag: explained peaks can differ.")}>
            <CandidateTable candidates={result.candidates} />
          </FragmentProvider>
        </>)}
      {workflowResult && <WorkflowResultView result={workflowResult} />}
      {job && <StepList steps={job.steps} />}
    </>
  );
  return pipeId === ESI_PIPELINE ? <StructureProvider result={result}>{content}</StructureProvider> : content;
}

export function PipelinePage({ id, spec, active }: { id: string; spec: Spec; active: boolean }) {
  const t = useT();
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
    ticket.current++; dropRun("workflow", id);
    setValues(initialValues(pipe ? pipe.inputs : [])); setJob(null); setBusy(false); setHasStartedJob(false); setError(""); setErrors({}); setRequestedAdduct(null);
    return () => { ticket.current++; };
  }, [id]);
  useEffect(() => { if (!active) guard.cancel(); }, [active, guard.cancel]);
  if (!pipe) return <p className="page">{t("Unknown pipeline.")}</p>;

  const shown = visibleFields(pipe.inputs, values);
  const update = (fid: string, v: any) => {
    setValues(prev => ({ ...prev, [fid]: v }));
    if (errors[fid]) setErrors(prev => { const next = { ...prev }; delete next[fid]; return next; });
    if (netNote) setNetNote("");
  };
  const loadExample = (example: Record<string, any>) => {
    ticket.current++; confirming.current = false; guard.cancel(); dropRun("workflow", pipe.id);
    setValues({ ...initialValues(pipe.inputs), ...example }); setJob(null); setBusy(false); setHasStartedJob(false); setError(""); setErrors({}); setRequestedAdduct(null); setNetNote("");
  };
  const clear = () => {
    ticket.current++;
    confirming.current = false;
    guard.cancel();
    dropRun("workflow", pipe.id);
    setValues(initialValues(pipe.inputs)); setJob(null); setBusy(false); setHasStartedJob(false); setError("");
    setErrors({}); setRequestedAdduct(null);
    setNetNote(busy ? t("Cleared. The previous workflow may still finish on the server.") : "");
    setFormVersion(n => n + 1);
  };

  const run = async () => {
    if (busy || confirming.current) return;
    const missing = missingRequired(pipe.inputs, values);
    if (missing.length) {
      setErrors(Object.fromEntries(missing.map(f => [f.id!, t("Required.")])));
      const labels = missing.map(f => f.label).join(", ");
      setError(missing.length === 1 ? t("Complete this required input: {labels}.", { labels }) : t("Complete these required inputs: {labels}.", { labels }));
      return;
    }
    const mine = ++ticket.current;
    const payload = Object.fromEntries(shown.filter(f => f.id).map(f => [f.id!, values[f.id!]]));
    // Internet search off: a workflow that would fetch a USI asks first.
    confirming.current = true;
    const allowed = await guard.confirm(networkUses(shown, values, pipe.network));
    if (!allowed) {
      if (ticket.current === mine) { confirming.current = false; setNetNote(t("Cancelled: internet search is off.")); }
      return;
    }
    if (ticket.current !== mine) return;
    // Local version: what the workflow needs and is not installed yet is downloaded first (a dialog asks before a big one). "Not now" stops here, quietly.
    const ready = await ensureSetup({ target: `workflow:${pipe.id}`, inputs: payload });
    if (ticket.current === mine) confirming.current = false;
    if (!ready || ticket.current !== mine) return;
    setNetNote("");
    setBusy(true); setHasStartedJob(false); setError(""); setErrors({}); setJob(null); dropRun("workflow", pipe.id);
    const t0 = performance.now();
    try {
      setRequestedAdduct(requestedAdductOf(pipe.id, payload));
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
        if (j.status !== "running" && j.status !== "queued") {
          if (j.error) setError(j.error);
          else if (j.result) recordRun({ kind: "workflow", id: pipe.id, title: pipe.title, section: "Workflows", inputs: payload, result: j.result, steps: j.steps, seconds: (performance.now() - t0) / 1000 });
          break;
        }
        await sleep(POLL_MS);
        if (ticket.current !== mine) return;
      }
    } catch (e: any) {
      if (ticket.current === mine) setError(e?.message || String(e) || t("An unknown error occurred."));
    } finally {
      if (ticket.current === mine) setBusy(false);
    }
  };

  const field = (f: Field, i: number) => (
    <div className={"field-slot" + (f.new_row ? " new-row" : "")} data-field={f.id} key={f.id ?? i}>
      <FieldInput field={f} value={f.id ? values[f.id] : undefined} error={f.id ? errors[f.id] : undefined} onChange={v => f.id && update(f.id, v)} />
    </div>);
  return (
    <div className="page module-page" data-pipeline={pipe.id}>
      <PipelineHead pipe={pipe} />
      <div key={formVersion} className={"panel" + (busy ? " is-running" : "")}>
        {groupFields(shown).map(g => {
          const body = g.fields.every(f => f.type === "checkbox")
            ? <div className="check-row">{g.fields.map(field)}</div>
            : <div className="fields">{g.fields.map(field)}</div>;
          return g.collapsed
            ? <details className="form-section" key={g.title}><summary>{t(g.title)}</summary>{body}</details>
            : <div className="form-section" key={g.title || "plain"}>{body}</div>;
        })}
        <div className="actions">
          <button className="btn" onClick={run} disabled={busy} aria-busy={busy}>{busy ? <><span className="spinner" aria-hidden="true" />{t("Running…")}</> : t("Run")}</button>
          <button className="btn ghost" type="button" onClick={clear}>{t("Clear")}</button>
          {(pipe.examples ?? []).map(ex => (
            <button className="btn secondary" type="button" key={ex.label} disabled={busy} onClick={() => loadExample(ex.values)}>{ex.label}</button>))}
          <span className="msg" aria-live="polite">{netNote}</span>
        </div>
      </div>
      {error && <div className="alert error" role="alert">{error}</div>}
      {(busy || job || (hasStartedJob && error)) && <WorkflowProgress job={job} awaitingFirstUpdate={hasStartedJob}
        trackingStopped={!!error && !busy && hasStartedJob && (!job || job.status === "running" || job.status === "queued")} />}
      <PipelineResults pipeId={pipe.id} jobResult={job?.result ?? null} job={job} requestedAdduct={requestedAdduct} />
      <PageFoot id={pipe.id} basedOn={pipe.based_on} />
      {guard.dialog}
    </div>
  );
}
