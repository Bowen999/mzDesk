import { useSyncExternalStore } from "react";
import { edmonton } from "./time";

/** One page that has a result in this browser tab. Pages register after a successful run and drop out when they are
 *  cleared, re-run or fail, so the list always matches what the person can see on screen. */
export type ReportEntry = {
  key: string;                       // module:<id> | workflow:<id>
  kind: "module" | "workflow";
  id: string;
  title: string;
  section: string;                   // sidebar group ("MS1", "Workflows", ...)
  inputs: Record<string, unknown>;   // what was sent to the server
  result: unknown;
  steps?: unknown;                   // workflows: the step list shown under the result
  seconds: number | null;
  finishedAt: string;                // ISO with the Edmonton offset
};

const entries = new Map<string, ReportEntry>();
const listeners = new Set<() => void>();
let snapshot: ReportEntry[] = [];

function changed() {
  snapshot = [...entries.values()].sort((a, b) => a.finishedAt.localeCompare(b.finishedAt));
  listeners.forEach(l => l());
}

export const reportKey = (kind: ReportEntry["kind"], id: string) => `${kind}:${id}`;

export function recordRun(e: Omit<ReportEntry, "key" | "finishedAt">) {
  const key = reportKey(e.kind, e.id);
  entries.set(key, { ...e, key, finishedAt: edmonton().iso });
  changed();
}

export function dropRun(kind: ReportEntry["kind"], id: string) {
  if (entries.delete(reportKey(kind, id))) changed();
}

export const listRuns = () => snapshot;

const subscribe = (l: () => void) => { listeners.add(l); return () => { listeners.delete(l); }; };
export const useRuns = () => useSyncExternalStore(subscribe, listRuns, listRuns);
