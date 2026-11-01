import { useEffect, useState } from "react";
import type { Spec } from "../types";
import { buildReport } from "./build";
import { useRuns } from "./log";
import { download } from "./zip";

function ReportIcon({ size = 16 }: { size?: number }) {
  return (
    <svg viewBox="0 0 24 24" width={size} height={size} fill="none" stroke="currentColor" strokeWidth="1.7"
      strokeLinecap="round" strokeLinejoin="round" aria-hidden="true" focusable="false">
      <path d="M14 3H7a2 2 0 0 0-2 2v14a2 2 0 0 0 2 2h10a2 2 0 0 0 2-2V8z" />
      <path d="M14 3v5h5M12 11v6M9.5 14.5 12 17l2.5-2.5" />
    </svg>
  );
}

/** "Session report" in the sidebar footer: downloads <run id>.zip with report.json + report.html for every page that has
 *  a result in this tab. Disabled until a page has a result. In the collapsed rail: an icon button. */
export function ReportButton({ spec, rail = false }: { spec: Spec; rail?: boolean }) {
  const runs = useRuns();
  const [busy, setBusy] = useState(false);
  const [note, setNote] = useState<{ text: string; bad: boolean } | null>(null);
  useEffect(() => {
    if (!note) return;
    const h = window.setTimeout(() => setNote(null), 12000);
    return () => window.clearTimeout(h);
  }, [note]);
  const n = runs.length;
  const title = n
    ? `Download a report folder (JSON + HTML, zipped) with the ${n} page${n === 1 ? "" : "s"} that have results`
    : "Session report: run a module or workflow first, then download its results here";
  const go = async () => {
    if (busy || !n) return;
    setBusy(true); setNote(null);
    try {
      const r = await buildReport(spec);
      if (!r) return;
      download(r.blob, r.filename);
      setNote({ text: `Saved ${r.filename}`, bad: false });
    } catch (e: any) {
      setNote({ text: `Report failed: ${e?.message ?? e}`, bad: true });
    } finally { setBusy(false); }
  };
  if (rail) return (
    <button type="button" className="rail-btn report-rail" disabled={!n || busy} onClick={go} title={note?.text ?? title} aria-label="Session report">
      <ReportIcon size={18} />
    </button>
  );
  return (
    <>
      <button type="button" className="internet-switch report-btn" disabled={!n || busy} onClick={go} title={title} aria-busy={busy}>
        <ReportIcon />
        <span className="internet-label">{busy ? "Building report…" : "Session report"}</span>
        {n > 0 && <span className="report-count" aria-label={`${n} pages with results`}>{n}</span>}
      </button>
      {note && <div className={"report-note" + (note.bad ? " bad" : "")} role="status">{note.text}</div>}
    </>
  );
}
