import type { RunResult, Spec } from "./types";

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

export const getSpec = () => json<Spec>("/api/spec");
export const runModule = (id: string, inputs: unknown) => post<RunResult>(`/api/modules/${id}/run`, inputs);
export const runPipeline = (id: string, inputs: unknown) => post<any>(`/api/pipelines/${id}/run`, inputs);
