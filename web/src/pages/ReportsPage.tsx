import { useRef, useState, type FormEvent } from "react";
import { ApiError, deleteReport, getReport, getReports, type ReportFull, type ReportRow } from "../api";

const CHANNEL_NAMES: Record<string, string> = { resend: "email", smtp: "email (SMTP)", ntfy: "push" };
const when = (iso: string) => new Date(iso).toLocaleString();

/** /reports: feedback and crashes sent from the app. The password is the admin password, checked by the server; it stays in this page's memory. */
export function ReportsPage() {
  const [pw, setPw] = useState("");
  const [rows, setRows] = useState<ReportRow[] | null>(null);
  const [channels, setChannels] = useState<string[]>([]);
  const [kind, setKind] = useState<"all" | "feedback" | "crash">("all");
  const [open, setOpen] = useState<ReportFull | null>(null);
  const [err, setErr] = useState("");
  const [busy, setBusy] = useState(false);
  const request = useRef(0);

  const fail = (e: unknown) => setErr(e instanceof ApiError && e.message ? e.message : String(e));
  const load = async () => {
    const mine = ++request.current;
    setBusy(true); setErr("");
    try {
      const r = await getReports(pw);
      if (mine !== request.current) return;
      setRows(r.reports); setChannels(r.channels);
    } catch (e) { if (mine === request.current) fail(e); }
    finally { if (mine === request.current) setBusy(false); }
  };
  const show = async (id: string) => {
    if (open?.id === id) { setOpen(null); return; }
    setErr("");
    try { setOpen(await getReport(pw, id)); } catch (e) { fail(e); }
  };
  const remove = async (id: string) => {
    try { await deleteReport(pw, id); setOpen(null); await load(); } catch (e) { fail(e); }
  };
  const submit = (e: FormEvent) => { e.preventDefault(); void load(); };

  if (!rows) return (
    <div className="page stats-page">
      <div className="module-heading"><div className="module-heading-copy"><h1>Reports</h1></div></div>
      <form className="panel stats-login" onSubmit={submit}>
        <label className="field">Password<input type="password" autoComplete="current-password" value={pw} onChange={e => setPw(e.target.value)} autoFocus /></label>
        <div className="actions"><button className="btn" disabled={busy || !pw}>{busy ? "Checking…" : "Show reports"}</button>{err && <span className="msg stats-err">{err}</span>}</div>
      </form>
    </div>
  );

  const shown = rows.filter(r => kind === "all" || r.kind === kind);
  const crashes = rows.filter(r => r.kind === "crash").length;
  return (
    <div className="page stats-page reports-page">
      <div className="module-heading"><div className="module-heading-copy"><h1>Reports</h1></div></div>
      <div className="stats-bar">
        <div className="chip-row">
          {([["all", `All ${rows.length}`], ["feedback", `Feedback ${rows.length - crashes}`], ["crash", `Crashes ${crashes}`]] as const).map(([k, label]) =>
            <button key={k} type="button" className={"chip" + (kind === k ? " on" : "")} onClick={() => setKind(k)}>{label}</button>)}
        </div>
        <span className={"msg" + (err ? " stats-err" : "")}>{err || (channels.length ? `Also sent by ${channels.map(c => CHANNEL_NAMES[c] ?? c).join(", ")}` : "Kept here only")}</span>
        <button type="button" className="btn secondary" onClick={() => void load()} disabled={busy}>Refresh</button>
      </div>
      {shown.length === 0 && <p className="stats-note">Nothing yet.</p>}
      <ul className="report-list">
        {shown.map(r => (
          <li key={r.id} className={"report" + (open?.id === r.id ? " open" : "")}>
            <button type="button" className="report-row" aria-expanded={open?.id === r.id} onClick={() => void show(r.id)}>
              <code className="report-ref">{r.id}</code>
              <span className={"badge " + (r.kind === "crash" ? "tone-warn" : "tone-neutral")}>{r.kind === "crash" ? "Crash" : "Feedback"}</span>
              <span className="report-title">{r.title || "(no text)"}</span>
              <span className="report-meta">{r.count > 1 ? `×${r.count} · ` : ""}{r.images ? `${r.images} img · ` : ""}{r.page} · v{r.version} · {when(r.last)}</span>
            </button>
            {open?.id === r.id && (
              <div className="report-detail">
                <pre className="report-body">{open.body}</pre>
                {open.detail && <pre className="report-body report-extra">{open.detail}</pre>}
                {open.images.length > 0 && <div className="report-shots">{open.images.map((im, i) => <a key={i} href={im.src} target="_blank" rel="noreferrer"><img src={im.src} alt={im.name} /></a>)}</div>}
                <p className="report-foot">{open.browser || "unknown browser"} · first {when(open.time)}</p>
                <button type="button" className="btn secondary" onClick={() => void remove(r.id)}>Delete</button>
              </div>)}
          </li>))}
      </ul>
    </div>
  );
}
