import { ModuleHead, ModuleResults, modulePageClass } from "../pages/ModulePage";
import { PipelineHead, PipelineResults, requestedAdductOf } from "../pages/PipelinePage";
import type { AnnotationResult, JobStep, RunResult, Spec, WorkflowResult } from "../types";
import type { ReportEntry } from "./log";
import { ReportModeContext, type ReportMode } from "./mode";

/**
 * A page as report.html shows it: the same title, result blocks and steps as on screen (the components of ModulePage / PipelinePage), drawn in
 * report mode (report/mode.ts) from the result the page logged. Inputs are not drawn here: build.ts puts them above the page.
 * Returns null for a page the spec no longer has.
 */
export function ReportPage({ spec, entry, mode }: { spec: Spec; entry: ReportEntry; mode: ReportMode }) {
  let page = null;
  if (entry.kind === "module") {
    const mod = spec.modules.find(m => m.id === entry.id);
    if (mod) page = (
      <div className={modulePageClass(mod)} data-module={mod.id}>
        <ModuleHead mod={mod} spec={spec} />
        <ModuleResults mod={mod} result={entry.result as RunResult} ok busy={false} seconds={entry.seconds} stale={false} />
      </div>);
  } else {
    const pipe = spec.pipelines.find(p => p.id === entry.id);
    if (pipe) page = (
      <div className="page module-page" data-pipeline={pipe.id}>
        <PipelineHead pipe={pipe} />
        <PipelineResults pipeId={pipe.id} jobResult={entry.result as AnnotationResult | WorkflowResult} requestedAdduct={requestedAdductOf(pipe.id, entry.inputs)}
          job={Array.isArray(entry.steps) ? { steps: entry.steps as JobStep[] } : null} />
      </div>);
  }
  return <ReportModeContext.Provider value={mode}>{page}</ReportModeContext.Provider>;
}
