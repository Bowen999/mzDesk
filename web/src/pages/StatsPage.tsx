import { useRef, useState, type FormEvent } from "react";
import { ApiError, getStats, type Stats } from "../api";
import type { Spec } from "../types";

/** /stats: usage numbers. The password is the site password, checked by the server; it is kept only in this page's memory. */
export function StatsPage({ spec }: { spec: Spec }) {
  const [pw, setPw] = useState("");
  const [days, setDays] = useState(30);
  const [data, setData] = useState<Stats | null>(null);
  const [err, setErr] = useState("");
  const [busy, setBusy] = useState(false);
  const request = useRef(0);
  const title = (id: string) => spec.modules.find(m => m.id === id)?.title ?? id;

  const load = async (d = days) => {
    const mine = ++request.current;
    setBusy(true); setErr("");
    try {
      const next = await getStats(pw, d);
      if (mine !== request.current) return;
      setData(next);
      setDays(next.days);
    } catch (e) {
      if (mine !== request.current) return;
      setErr(e instanceof ApiError && e.message ? e.message : String(e));
    } finally {
      if (mine === request.current) setBusy(false);
    }
  };
  const submit = (e: FormEvent) => { e.preventDefault(); void load(); };

  if (!data) return (
    <div className="page stats-page">
      <div className="module-heading"><div className="module-heading-copy"><div className="module-kicker">mzDesk</div><h1>Usage statistics</h1></div></div>
      <form className="panel stats-login" onSubmit={submit}>
        <label className="field">Password<input type="password" autoComplete="current-password" value={pw} onChange={e => setPw(e.target.value)} autoFocus /></label>
        <div className="actions"><button className="btn" disabled={busy || !pw}>{busy ? "Checking…" : "Show statistics"}</button>{err && <span className="msg stats-err">{err}</span>}</div>
      </form>
    </div>
  );

  const t = data.totals, max = Math.max(1, ...data.daily.map(d => Math.max(d.views, d.runs)));
  const card = (label: string, v: { views: number; runs: number; visitors: number }) => (
    <div className="card"><div className="card-label">{label}</div>
      <div className="card-value">{v.visitors}</div><div className="card-unit">visitors · {v.views} views · {v.runs} runs</div></div>);
  return (
    <div className="page stats-page">
      <div className="module-heading"><div className="module-heading-copy"><div className="module-kicker">mzDesk</div><h1>Usage statistics</h1></div></div>
      <div className="stats-bar">
        <div className="chip-row">{[7, 30, 90, 365].map(d => <button key={d} type="button" className={"chip" + (d === days ? " on" : "")} onClick={() => void load(d)}>{d} days</button>)}</div>
        <span className={"msg" + (err ? " stats-err" : "")}>{err || (busy ? "Loading…" : `Updated ${new Date(data.generated).toLocaleString()}`)}</span>
        <button type="button" className="btn secondary" onClick={() => void load()} disabled={busy}>Refresh</button>
      </div>
      <div className="outputs has-data">
        <div className="block span-full tier-hero"><div className="out"><div className="cards metric">
          {card("Today", t.today)}{card("Last 7 days", t["7d"])}{card(`Last ${data.days} days`, t.range)}{card("All time", t.all)}
        </div></div></div>
        <div className="block span-full tier-card"><div className="out">
          <h4 className="out-title">Daily activity <span className="stats-legend"><i className="lg-views" />views <i className="lg-runs" />runs</span></h4>
          <div className="stats-chart-scroll">
            <div className="stats-chart-inner" style={{ minWidth: `${data.daily.length * 5}px` }}>
              <div className="stats-chart" role="img" aria-label="Views and runs per day">
                {data.daily.map(d => (
                  <div key={d.day} className="stats-col" title={`${d.day}: ${d.visitors} visitors, ${d.views} views, ${d.runs} runs`}>
                    <span className="bar views" style={{ height: `${(100 * d.views) / max}%` }} />
                    <span className="bar runs" style={{ height: `${(100 * d.runs) / max}%` }} />
                  </div>))}
              </div>
              <div className="stats-axis"><span>{data.daily[0]?.day}</span><span>{data.daily[data.daily.length - 1]?.day}</span></div>
            </div>
          </div>
        </div></div>
        <Table title="Module runs" cols={["Module", "Runs", "Errors", "Avg time"]}
          rows={data.modules.map(m => [title(m.module), m.runs, m.errors, `${(m.avg_ms / 1000).toFixed(1)} s`])} />
        <Table title="Pages" cols={["Page", "Views"]} rows={data.pages.map(p => [p.path, p.views])} half />
        <Table title="Referrers" cols={["Site", "Views"]} rows={data.referrers.map(r => [r.host, r.views])} half />
        <Table title="Browsers" cols={["Browser", "Visitors"]} rows={data.browsers.map(b => [b.browser, b.visitors])} half />
        <Table title="Recent errors" cols={["Time (UTC)", "Module", "Status"]} rows={data.errors.map(e => [e.time.replace("T", " ").slice(0, 16), title(e.module), e.status])} half />
      </div>
      <p className="stats-note">Visitors are counted once per day per browser; no IP addresses are stored. Bots are left out of the totals.</p>
    </div>
  );
}

function Table({ title, cols, rows, half }: { title: string; cols: string[]; rows: (string | number)[][]; half?: boolean }) {
  return (
    <div className={`block ${half ? "span-half tier-quiet" : "span-full tier-card"}`}><div className="out">
      <h4 className="out-title">{title}<span className="count-pill">{rows.length}</span></h4>
      <div className="table-wrap"><table className="out-table">
        <thead><tr>{cols.map(c => <th key={c}>{c}</th>)}</tr></thead>
        <tbody>{rows.length ? rows.map((r, i) => <tr key={i}>{r.map((c, j) => <td key={j}>{c}</td>)}</tr>)
          : <tr><td colSpan={cols.length} className="empty">nothing yet</td></tr>}</tbody>
      </table></div>
    </div></div>);
}
