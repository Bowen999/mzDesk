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

/** Where the app is served: "/" locally, "/mzdesk/" behind bowenyang.org (set at build time: MZDESK_BASE=/mzdesk/ npm run build). */
export const BASE = import.meta.env.BASE_URL;
export const getSpec = () => json<Spec>(`${BASE}api/spec`);
/** The long notes of a page (Docs page): markdown text. */
export const getDocs = (id: string) => json<{ id: string; markdown: string }>(`${BASE}api/docs/${id}`);
export const runModule = (id: string, inputs: unknown) => post<RunResult>(`${BASE}api/modules/${id}/run`, inputs);
/** A workflow runs as a background job: start it, then poll getJob until its status is not "running". */
export const startPipeline = (id: string, inputs: unknown) => post<{ job_id: string; status: string }>(`${BASE}api/pipelines/${id}/run`, inputs);
export const getJob = (jobId: string) => json<Job>(`${BASE}api/pipelines/jobs/${jobId}`);

/** The fragment view ("View" button): computed when asked. Start the job, then poll it until its status is done or error. */
export const startFragments = (body: unknown) => post<{ job_id: string; status: string }>(`${BASE}api/fragments/run`, body);
export const getFragmentJob = (jobId: string) => json<FragmentJob>(`${BASE}api/fragments/jobs/${jobId}`);

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
export type FeedbackResult = { ok: boolean; delivered: boolean; note: string };
export const sendFeedback = (body: {
  message: string; page: string; website: string;
  images: { name: string; data: string }[];
}) => post<FeedbackResult>(`${BASE}api/feedback`, body);
