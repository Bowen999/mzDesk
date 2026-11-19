import { useEffect, useSyncExternalStore } from "react";
import {
  ApiError, cancelInstall, getInstallJob, getLocalStatus, getSetupNeeds, getUpdate, startInstall,
  type SetupBlocked, type SetupJob, type SetupJobComponent, type SetupMissing, type SetupNeeds, type SetupUnsupported, type UpdateInfo,
} from "./api";
import { t } from "./i18n";
import { waitingPercent } from "./progressCurve";

/**
 * First-use downloads and the update notice (local version only; the public website never asks).
 *   ensureSetup({ target, inputs })  before a Run: asks the server what is missing (POST /api/local/needs), shows the one dialog (SetupDialog.tsx)
 *                                    and resolves true when the run may go on, false when the person said "Not now" / Cancel / Close.
 *   ensureSetup({ ids })             a checkbox that needs components (field `needs`): the same dialog for those ids.
 * One flow at a time: while a dialog is open a second call (a double click) answers false at once. A server that cannot be asked (an older
 * build, a dropped request) never blocks a run: ensureSetup then resolves true without a dialog.
 * The dialog reads this store through useSetup(); the buttons call confirmSetup / dismissSetup / retrySetup.
 */
export type SetupRequest = { target: string; inputs: unknown } | { ids: string[] };
export type SetupPhase = "confirm" | "progress" | "blocked" | "error";
export type SetupState = {
  phase: SetupPhase | null;            // null = no dialog
  forRun: boolean;                     // the dialog belongs to a Run (not to a tick box)
  missing: SetupMissing[]; blocked: SetupBlocked[]; unsupported: SetupUnsupported[];
  job: SetupJob | null;
  error: string; log: string[];
  cancelling: boolean;
  seen: Record<string, number>;        // component id -> Date.now() when it was first seen running (the bar's clock)
};
/** A row of the progress view. step / progress / gotMb / totalMb are empty on a server that does not report them (see SetupJobComponent). */
export type SetupRow = {
  id: string; title: string; size: string; status: SetupJobComponent["status"]; detail: string;
  step?: string; progress?: number | null; gotMb?: number | null; totalMb?: number | null;
};

const CLOSED: SetupState = { phase: null, forRun: false, missing: [], blocked: [], unsupported: [], job: null, error: "", log: [], cancelling: false, seen: {} };
const POLL_MS = 1000;
const FLOOR = 6;                                // the bar is never empty while it works
const MAX_POLL_FAILURES = 3;                    // a dropped request is retried; an answer of the server (a job that is gone after a restart) ends it at once
const sleep = (ms: number) => new Promise<void>(resolve => setTimeout(resolve, ms));

type Flow = { request: SetupRequest; resolve: (ok: boolean) => void; jobId: string | null };
let state: SetupState = CLOSED;
let flow: Flow | null = null;
let local = false;
const listeners = new Set<() => void>();
const emit = () => listeners.forEach(l => l());
const subscribe = (l: () => void) => { listeners.add(l); return () => { listeners.delete(l); }; };
const set = (patch: Partial<SetupState>) => { state = { ...state, ...patch }; emit(); };

/** Set once the spec arrives (App.tsx): only the local version (`hosted: false`) downloads components and checks for updates. */
export function setLocalApp(on: boolean) { local = on; }
/** The local version (`hosted: false` in the spec)? False until the spec arrives. */
export const isLocalApp = () => local;
export const getSetupState = () => state;
export const useSetup = () => useSyncExternalStore(subscribe, getSetupState, getSetupState);

function finish(ok: boolean) {
  const f = flow;
  flow = null;
  state = CLOSED;
  emit();
  f?.resolve(ok);
}

// ---------------------------------------------------------------- what is missing
/** The same answer as POST /api/local/needs, for a tick box: the listed ids against GET /api/local/status. */
async function needsOfIds(ids: string[]): Promise<SetupNeeds> {
  const status = await getLocalStatus();
  const wanted = ids.map(id => status.components.find(c => c.id === id)).filter((c): c is NonNullable<typeof c> => !!c);
  return {
    missing: wanted.filter(c => c.supported && !c.installed).map(c => ({ id: c.id, title: c.title, size: c.size, ask: c.ask, why: c.note })),
    blocked: [],
    unsupported: wanted.filter(c => !c.supported).map(c => ({ id: c.id, title: c.title, message: c.note })),
  };
}
const needsOf = (r: SetupRequest) => ("ids" in r ? needsOfIds(r.ids) : getSetupNeeds(r.target, r.inputs));

export async function ensureSetup(request: SetupRequest): Promise<boolean> {
  if (!local) return true;
  if (flow) return false;                                       // one dialog for the whole app: a double click stops quietly
  const mine: Flow = { request, resolve: () => undefined, jobId: null };
  const done = new Promise<boolean>(resolve => { mine.resolve = resolve; });
  flow = mine;
  let needs: SetupNeeds;
  try { needs = await needsOf(request); }
  catch { if (flow === mine) finish(true); return done; }       // cannot ask: the run goes on, as it did before this feature
  if (flow !== mine) return done;
  const missing = needs.missing ?? [], blocked = needs.blocked ?? [], unsupported = needs.unsupported ?? [];
  if (!missing.length && !blocked.length) { finish(true); return done; }
  state = { ...CLOSED, forRun: !("ids" in request), missing, blocked, unsupported };
  if (blocked.length) set({ phase: "blocked" });                // cannot be fixed by downloading: the run is cancelled
  else if (missing.some(m => m.ask)) set({ phase: "confirm" });
  else void install();
  return done;
}

// ---------------------------------------------------------------- download
async function install() {
  const f = flow;
  if (!f) return;
  set({ phase: "progress", error: "", log: [], job: null, cancelling: false, seen: {} });
  try {
    const { job_id } = await startInstall(state.missing.map(m => m.id));
    if (flow !== f) { cancelInstall(job_id).catch(() => undefined); return; }          // cancelled while the request was on its way
    f.jobId = job_id;
    await watch(f);
  } catch (e: any) {
    if (flow === f) fail(e instanceof ApiError ? e.message : String(e?.message ?? e) || t("The download did not finish."), e instanceof ApiError ? e.log.split("\n") : []);
  }
}

function fail(error: string, log: string[]) { set({ phase: "error", error, log: log.filter(Boolean).slice(-30), cancelling: false }); }

async function watch(f: Flow) {
  let failures = 0;
  for (;;) {
    let job: SetupJob;
    try { job = await getInstallJob(f.jobId!); failures = 0; }
    catch (e) {
      if (flow !== f) return;
      if (e instanceof ApiError || ++failures >= MAX_POLL_FAILURES) throw e;
      await sleep(POLL_MS);
      continue;
    }
    if (flow !== f) return;
    const seen = { ...state.seen };
    for (const c of job.components ?? []) if (c.status === "running" && !(c.id in seen)) seen[c.id] = Date.now();
    set({ job, seen });
    if (job.status === "done") return verify(f);
    if (job.status === "error") {
      const bad = (job.components ?? []).find(c => c.status === "error");
      fail(job.error || bad?.detail || t("The download did not finish."), job.log ?? []);
      return;
    }
    if (job.status === "cancelled") { finish(false); return; }   // stopped from somewhere else
    await sleep(POLL_MS);
    if (flow !== f) return;
  }
}

/** The job says done: ask once more, so a run never starts on a component the server does not see (a job another tab started may have installed other ones). */
async function verify(f: Flow) {
  let left: SetupMissing[] = [];
  try { left = (await needsOf(f.request)).missing ?? []; } catch { /* cannot ask: take the job's word */ }
  if (flow !== f) return;
  if (!left.length) { finish(true); return; }
  set({ missing: left });
  fail(t("{names} could not be found after the download.", { names: left.map(m => m.title).join(", ") }), []);
}

// ---------------------------------------------------------------- the dialog's buttons
/** [Download] */
export function confirmSetup() { if (flow && state.phase === "confirm") void install(); }
/** [Retry] */
export function retrySetup() { if (flow && state.phase === "error") void install(); }
/** [Cancel] while it downloads: stops the server's child process, then closes. */
export async function cancelSetup() {
  const f = flow;
  if (!f || state.cancelling) return;
  set({ cancelling: true });
  if (f.jobId) { try { await cancelInstall(f.jobId); } catch { /* it may have ended on its own */ } }
  if (flow === f) finish(false);
}
/** [Not now] / [Close] / Escape: the run is cancelled and its inputs stay. While it downloads this is [Cancel]. */
export function dismissSetup() {
  if (!flow) return;
  if (state.phase === "progress") void cancelSetup(); else finish(false);
}

// ---------------------------------------------------------------- sizes and the bar
const UNIT_GB: Record<string, number> = { KB: 1e-6, MB: 1e-3, GB: 1, TB: 1e3 };
/** "2.3 GB" -> 2.3, "600 MB" -> 0.6; null when the text is not a size. */
export function sizeGb(size: string): number | null {
  const m = /^\s*(\d+(?:[.,]\d+)?)\s*(KB|MB|GB|TB)\b/i.exec(size);
  return m ? parseFloat(m[1].replace(",", ".")) * UNIT_GB[m[2].toUpperCase()] : null;
}

/** The sum of the sizes as the dialog writes it ("2.9 GB", "20 MB"); null when one of them is not a size. */
export function totalSize(items: { size: string }[]): string | null {
  const parts = items.map(i => sizeGb(i.size));
  if (!parts.length || parts.some(p => p === null)) return null;
  const gb = (parts as number[]).reduce((a, b) => a + b, 0);
  return gb >= 0.95 ? `${Math.round(gb * 10) / 10} GB` : `${Math.max(1, Math.round(gb * 1000))} MB`;
}

const num = (v: unknown): number | null => (typeof v === "number" && Number.isFinite(v) && v >= 0 ? v : null);

/** The rows of the progress view: the job's components, with the sizes of the needs; before the first answer of the job, every missing one waits. */
export function setupRows(s: SetupState): SetupRow[] {
  const size = (id: string) => s.missing.find(m => m.id === id)?.size ?? "";
  if (!s.job?.components?.length) return s.missing.map(m => ({ id: m.id, title: m.title, size: m.size, status: "waiting", detail: "", step: "", progress: null, gotMb: null, totalMb: null }));
  return s.job.components.map(c => {
    const gotMb = num(c.got_mb), totalMb = num(c.total_mb);
    let progress = num(c.progress);
    if (progress === null && gotMb !== null && totalMb) progress = gotMb / totalMb;
    return {
      id: c.id, title: c.title, size: size(c.id), status: c.status, detail: c.detail ?? "",
      step: typeof c.step === "string" ? c.step.trim() : "", progress: progress === null ? null : Math.min(1, progress), gotMb, totalMb,
    };
  });
}

/** A size in megabytes without its unit, in the unit it is read in: whole MB with thousands separators, or GB with one decimal from
 *  10 GB and two below (trailing zeros dropped). */
function inUnit(mb: number, unit: "MB" | "GB", gbDecimals: number): string {
  if (unit === "MB") return Math.round(mb).toLocaleString("en-US");
  return String(parseFloat((mb / 1000).toFixed(gbDecimals)));
}

/** Megabytes as the dialog writes them: "412 MB", from 1000 MB on "1.24 GB". */
export function formatMb(mb: number): string {
  const unit = mb >= 1000 ? "GB" : "MB";
  return `${inUnit(mb, unit, mb >= 10000 ? 1 : 2)} ${unit}`;
}

/** "412 of 2,300 MB" for the file that downloads now (both numbers in the unit of the total: "0.41 of 2.3 GB" from 1000 MB on);
 *  "412 MB" when the total is unknown; "18%" when only a fraction is known; "" without a number. */
export function transferText(r: Pick<SetupRow, "gotMb" | "totalMb" | "progress">): string {
  const got = num(r.gotMb), total = num(r.totalMb);
  if (got !== null && total) {
    const unit = total >= 1000 ? "GB" : "MB", decimals = total >= 10000 ? 1 : 2;
    return t("{got} of {total}", { got: inUnit(Math.min(got, total), unit, decimals), total: `${inUnit(total, unit, decimals)} ${unit}` });
  }
  if (got !== null) return formatMb(got);
  const p = num(r.progress);
  return p === null ? "" : `${Math.round(Math.min(p, 1) * 100)}%`;
}

/** Every row is done but the job has not said so yet (it checks the install): the bar is full and says so. */
export const isFinishing = (rows: SetupRow[]) => rows.length > 0 && rows.every(r => r.status === "done");

/** The bar. Each component is a slice as wide as its download (a size we cannot read counts as 0.3 GB). A finished one fills its slice. A
 *  running one fills it by the real fraction of the file that downloads now when the server reports it (progress, capped at 99% of the slice),
 *  otherwise eased on time: fast at first and slower and slower (progressCurve.ts waitingPercent: towards 90%, slower the bigger it is).
 *  Never empty, never 100 before the job ends. (The dialog keeps the highest value it showed: a new file starts at 0 again.) */
export function setupPercent(rows: SetupRow[], seen: Record<string, number>, now: number): number {
  let total = 0, got = 0;
  for (const r of rows) {
    const w = Math.max(sizeGb(r.size) ?? 0.3, 0.02);
    total += w;
    if (r.status === "done") got += w;
    else if (r.status === "running") {
      const real = num(r.progress);
      if (real !== null) got += w * Math.min(real, 0.99);
      else {
        const seconds = Math.max(0, (now - (seen[r.id] ?? now)) / 1000);
        got += w * waitingPercent(seconds, { start: 3, cap: 90, tau: 20 + 60 * w }) / 100;
      }
    }
  }
  return total ? Math.min(99, Math.max(FLOOR, 100 * got / total)) : FLOOR;
}

/** "2:05" */
export function clock(seconds: number): string {
  const s = Math.max(0, Math.floor(seconds));
  return `${Math.floor(s / 60)}:${String(s % 60).padStart(2, "0")}`;
}

// ---------------------------------------------------------------- the update notice
const DISMISSED = "mzdesk.updateLater";           // the version the person put off with [Later]
let update: UpdateInfo | null = null;
let dismissed: string | null | undefined;           // undefined = not read yet
let updateAsked = false;

function dismissedVersion(): string | null {
  if (dismissed === undefined) { try { dismissed = localStorage.getItem(DISMISSED); } catch { dismissed = null; } }
  return dismissed;
}

/** The update to show: a newer version the person has not put off. */
export function visibleUpdate(info: UpdateInfo | null, later: string | null): UpdateInfo | null {
  return info?.newer && info.latest && info.latest !== later ? info : null;
}

/** GET /api/update, once per page load. Any failure means no notice. */
export function checkForUpdate() {
  if (!local || updateAsked) return;
  updateAsked = true;
  getUpdate().then(info => { update = info; emit(); }).catch(() => undefined);
}

/** The update the sidebar footer shows, or null. Asks the server once per page load. */
export function useUpdate(): UpdateInfo | null {
  useEffect(() => { checkForUpdate(); }, []);
  const shown = () => visibleUpdate(update, dismissedVersion());
  return useSyncExternalStore(subscribe, shown, shown);
}

/** [Later]: this version is not shown again (a newer one is). Remembered in localStorage when it works. */
export function dismissUpdate(version: string) {
  dismissed = version;
  try { localStorage.setItem(DISMISSED, version); } catch { /* remembered until the page is closed */ }
  emit();
}

/** The server's answer, for the tests. */
export function setUpdateInfo(info: UpdateInfo | null) { update = info; updateAsked = true; emit(); }
/** Forget the person's choice and the update, for the tests. */
export function resetUpdate() { update = null; dismissed = undefined; updateAsked = false; emit(); }
