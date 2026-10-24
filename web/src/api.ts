import type { Spec } from "./types";

async function json<T>(url: string, init?: RequestInit): Promise<T> {
  const r = await fetch(url, init);
  if (!r.ok) throw new Error(`${r.status} ${r.statusText}`);
  return r.json();
}
const post = (url: string, body: unknown) =>
  json<any>(url, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) });

export const getSpec = () => json<Spec>("/api/spec");
export const runModule = (id: string, inputs: unknown) => post(`/api/modules/${id}/run`, inputs);
export const runPipeline = (id: string, inputs: unknown) => post(`/api/pipelines/${id}/run`, inputs);
