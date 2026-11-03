import type { FragmentJob, Job, RunResult, Spec } from "./types";

/** 后端返回的错误：message 给用户看，log 是外部工具（如 MetFrag）的输出末尾。 */
export class ApiError extends Error {
  log: string;
  constructor(message: string, log = "") { super(message); this.log = log; }
}

async function json<T>(url: string, init?: RequestInit): Promise<T> {
  const r = await fetch(url, init);
  if (!r.ok) {
    let detail = `${r.status} ${r.statusText}`, log = "";
    try { const body = await r.json(); detail = typeof body.detail === "string" ? body.detail : JSON.stringify(body.detail); log = body.log ?? ""; }
    catch { /* not JSON */ }
    throw new ApiError(detail, log);
  }
  return r.json();
}
const post = <T,>(url: string, body: unknown) =>
  json<T>(url, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) });

/** Where the app is served: "/" (mzdesk.com and the local version); a sub-path only when built with MZDESK_BASE=/sub/ npm run build. */
export const BASE = import.meta.env.BASE_URL;
export const getSpec = () => json<Spec>(`${BASE}api/spec`);

/** What the results depend on (GET /api/versions): mzDesk, Python, engines, models, data, online services. Names, versions, dates and counts only.
 *  An online service has no version (kind "online"); a version of null with a detail means "not installed here" or "unknown". */
export type VersionEntry = {
  version: string | null;
  kind: "python-package" | "engine" | "model" | "data" | "online" | "mzdesk";
  detail: string | null;
};
export type Versions = { mzdesk: string; python: string; checked_at: string; components: Record<string, VersionEntry> };
export const getVersions = (signal?: AbortSignal) => json<Versions>(`${BASE}api/versions`, { signal });

/** The long notes of a page (Docs page): markdown text. */
export const getDocs = (id: string) => json<{ id: string; markdown: string }>(`${BASE}api/docs/${id}`);
export const runModule = (id: string, inputs: unknown) => post<RunResult>(`${BASE}api/modules/${id}/run`, inputs);
/** A workflow runs as a background job: start it, then poll getJob until its status is not "running". */
export const startPipeline = (id: string, inputs: unknown) => post<{ job_id: string; status: string }>(`${BASE}api/pipelines/${id}/run`, inputs);
export const getJob = (jobId: string) => json<Job>(`${BASE}api/pipelines/jobs/${jobId}`);

/** The fragment view ("View" button): computed when asked. Start the job, then poll it until its status is done or error. */
export const startFragments = (body: unknown) => post<{ job_id: string; status: string }>(`${BASE}api/fragments/run`, body);
export const getFragmentJob = (jobId: string) => json<FragmentJob>(`${BASE}api/fragments/jobs/${jobId}`);

/** The local version only (mzdesk/api/local.py): first-use downloads ("components") and the update check. See setup.ts. */
export type LocalComponent = { id: string; title: string; size: string; installed: boolean; supported: boolean; ask: boolean; note: string };
export type SetupMissing = { id: string; title: string; size: string; ask: boolean; why: string };
export type SetupBlocked = { id: string; title: string; message: string; link: string | null };
export type SetupUnsupported = { id: string; title: string; message: string };
export type SetupNeeds = { missing: SetupMissing[]; blocked: SetupBlocked[]; unsupported: SetupUnsupported[] };
/** One component of an install job. `step`, `progress`, `got_mb` and `total_mb` come from newer servers only (all optional):
 *  step = what the running component does right now ("Unpacking"; "" when unknown); progress = 0..1 of the file downloading right now
 *  (null when unknown, and null again when the next file starts), got_mb / total_mb = the same in megabytes. */
export type SetupJobComponent = {
  id: string; title: string; status: "waiting" | "running" | "done" | "error"; detail: string;
  step?: string; progress?: number | null; got_mb?: number | null; total_mb?: number | null;
};
export type SetupJob = {
  status: "running" | "done" | "error" | "cancelled"; current: string | null; elapsed_s: number;
  components: SetupJobComponent[]; log: string[]; error: string;
};
export type UpdateInfo = { current: string; latest: string | null; newer: boolean; released: string | null; notes: string; command: string; checked_at: string | null };
export const getLocalStatus = () => json<{ local: boolean; platform: string; components: LocalComponent[] }>(`${BASE}api/local/status`);
/** What a Run needs before it can start. `target` is "module:<id>" or "workflow:<id>"; `inputs` the form values the run sends. */
export const getSetupNeeds = (target: string, inputs: unknown) => post<SetupNeeds>(`${BASE}api/local/needs`, { target, inputs });
export const startInstall = (ids: string[]) => post<{ job_id: string }>(`${BASE}api/local/install`, { ids });
export const getInstallJob = (jobId: string) => json<SetupJob>(`${BASE}api/local/jobs/${jobId}`);
export const cancelInstall = (jobId: string) => post<{ ok: boolean }>(`${BASE}api/local/jobs/${jobId}/cancel`, {});
export const getUpdate = () => json<UpdateInfo>(`${BASE}api/update`);

/** Page view beacon for the stats page (fire and forget; never shows an error). */
export function trackView(path: string) {
  try {
    void fetch(`${BASE}api/track`, { method: "POST", headers: { "Content-Type": "application/json" }, keepalive: true,
      body: JSON.stringify({ path, ref: document.referrer }) }).catch(() => undefined);
  } catch { /* ignore */ }
}

export type Stats = {
  days: number; generated: string;
  totals: Record<"today" | "7d" | "range" | "all", { views: number; runs: number; visitors: number }>;
  daily: { day: string; views: number; runs: number; visitors: number }[];
  modules: { module: string; runs: number; errors: number; avg_ms: number }[];
  pages: { path: string; views: number }[];
  referrers: { host: string; views: number }[];
  browsers: { browser: string; visitors: number }[];
  errors: { time: string; module: string; status: number }[];
};
export const getStats = (password: string, days: number) => post<Stats>(`${BASE}api/stats`, { password, days });

/** The sidebar's Send feedback form. Images are data URLs, already scaled down by the browser. */
export type FeedbackResult = { ok: boolean; delivered: boolean; note: string; ref?: string };
export const sendFeedback = (body: {
  message: string; page: string; website: string;
  images: { name: string; data: string }[];
}) => post<FeedbackResult>(`${BASE}api/feedback`, body);

/** A page that failed to draw: the error, its stack and the page path (no form data). Fire and forget; resolves to the report's reference, or "" when it could not be sent. */
export async function reportCrash(error: Error, componentStack = ""): Promise<string> {
  try {
    const r = await post<{ ok: boolean; ref: string }>(`${BASE}api/crash`, {
      message: error.message || String(error), stack: `${error.stack ?? ""}${componentStack ? `\n--- components${componentStack}` : ""}`, page: window.location.pathname,
    });
    return r.ref;
  } catch { return ""; }
}

/** The maintainer's reports page (admin password, checked by the server). */
export type ReportRow = { id: string; time: string; last: string; kind: "feedback" | "crash"; page: string; title: string; version: string; count: number; images: number };
export type ReportFull = Omit<ReportRow, "images"> & { body: string; detail: string; browser: string; images: { name: string; src: string }[] };
export const getReports = (password: string) => post<{ reports: ReportRow[]; channels: string[] }>(`${BASE}api/reports`, { password });
export const getReport = (password: string, id: string) => post<ReportFull>(`${BASE}api/reports/${id}`, { password });
export const deleteReport = (password: string, id: string) => post<{ ok: boolean }>(`${BASE}api/reports/${id}/delete`, { password });
