import { useEffect, useMemo, useRef, useState } from "react";
import { darkColors } from "../darkColors";
import { useTheme } from "../themeMode";
import type { BlockData, PeakCard } from "../types";

const MONO = '"IBM Plex Mono",ui-monospace,"SF Mono",Menlo,Consolas,monospace';
const SANS = '"IBM Plex Sans",-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif';
const INK = "#161A23", MUTED = "#6A7180";     // peak labels; recoloured with the figure in the dark theme
const MAX_LABELS = 40;
const CONFIG = {
  responsive: true, displaylogo: false, displayModeBar: "hover",
  modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"],
  toImageButtonOptions: { format: "svg", filename: "mzdesk_figure" },
};

// plotly.js（约 1.4 MB gzip）和 3Dmol.js 只在第一次画图时加载
let plotlyPromise: Promise<any> | null = null;
const loadPlotly = () => (plotlyPromise ??= import("plotly.js-dist-min").then(m => m.default ?? m));
let mol3dPromise: Promise<any> | null = null;
function load3Dmol(): Promise<any> {
  return (mol3dPromise ??= new Promise((resolve, reject) => {
    const w = window as any;
    if (w.$3Dmol) return resolve(w.$3Dmol);
    const s = document.createElement("script");
    s.src = `${import.meta.env.BASE_URL}lib/3Dmol-min.js`;
    s.onload = () => resolve(w.$3Dmol);
    s.onerror = () => { mol3dPromise = null; reject(new Error("3Dmol.js failed to load")); };
    document.head.appendChild(s);
  }));
}

/** An SVG drawn for a white page (RDKit structures), recoloured in the dark theme. */
export function ThemedSvg({ svg, className = "fig-img svg", alt = "" }: { svg: string; className?: string; alt?: string }) {
  const { dark } = useTheme();
  const src = useMemo(() => `data:image/svg+xml;charset=utf-8,${encodeURIComponent(dark ? darkColors(svg) : svg)}`, [svg, dark]);
  return <img className={className} alt={alt} src={src} />;
}

export function FigureView({ data }: { data: BlockData }) {
  if (data.plotly) return <PlotlyFigure fig={data.plotly} peaks={data.peaks} />;
  if (data.svg) return <ThemedSvg svg={data.svg} />;
  if (data.img) return <img className="fig-img" alt="" src={data.img} />;
  if (data.mol3d) return <Mol3D molblock={data.mol3d} />;
  return null;
}

/** 小数位：全是整数 m/z（EI 单位质量谱）就不写 .0000，否则 4 位。 */
const mzDecimals = (peaks: PeakCard[]) => (peaks.every(p => Number.isInteger(p.mz)) ? 0 : 4);

/** m/z 标签：重要的峰优先（已解释，再按强度）；标签不与别的标签重叠，也不压住别的峰。缩放后重算。 */
function peakLabels(gd: any, peaks: PeakCard[], staticAnn: any[], precursor?: number) {
  const { xaxis: xa, yaxis: ya } = gd._fullLayout;
  const x0 = Math.min(...xa.range), x1 = Math.max(...xa.range), yTop = ya.range[1];
  const vis = peaks.map((_, i) => i).filter(i => peaks[i].mz >= x0 && peaks[i].mz <= x1);
  const sticks = vis.map(i => ({ i, x: xa.l2p(peaks[i].mz), top: ya.l2p(Math.min(peaks[i].rel, yTop)) }));
  if (precursor !== undefined && precursor >= x0 && precursor <= x1) sticks.push({ i: -1, x: xa.l2p(precursor), top: 0 });
  const dec = mzDecimals(peaks);
  const order = [...vis].sort((a, b) => (Number(peaks[b].on) - Number(peaks[a].on)) || (peaks[b].rel - peaks[a].rel));
  const boxes: { l: number; r: number; t: number; b: number }[] = [];
  const anns: any[] = [];
  for (const i of order) {
    if (anns.length >= MAX_LABELS) break;
    const p = peaks[i];
    if (p.rel > yTop) continue;
    const txt = p.mz.toFixed(dec), w = txt.length * 6.1 + 4, px = xa.l2p(p.mz), py = ya.l2p(p.rel);
    const b = { l: px - w / 2, r: px + w / 2, t: py - 16, b: py - 2 };
    if (b.l < 0 || b.r > xa._length || b.t < 0) continue;
    const isPrecursor = /precursor/i.test(p.kind);
    const hit = boxes.some(o => b.l < o.r + 4 && b.r > o.l - 4 && b.t < o.b + 1 && b.b > o.t - 1)
      || sticks.some(s => s.i !== i && !(s.i === -1 && isPrecursor) && s.x > b.l - 1 && s.x < b.r + 1 && s.top < b.b);
    if (hit) continue;
    boxes.push(b);
    anns.push({ x: p.mz, y: p.rel, xref: "x", yref: "y", text: txt, showarrow: false, xanchor: "center", yanchor: "bottom",
                yshift: 2, font: { family: MONO, size: 10, color: p.on ? INK : MUTED } });
  }
  return staticAnn.concat(anns);
}

function PlotlyFigure({ fig: raw, peaks }: { fig: { data: any[]; layout: any }; peaks?: PeakCard[] }) {
  const ref = useRef<HTMLDivElement>(null);
  const { dark } = useTheme();
  // the backend draws for a white page; in the dark theme every colour of the figure (and of the peak labels) is remapped
  const fig = useMemo(() => (dark ? JSON.parse(darkColors(JSON.stringify(raw), true)) : raw), [raw, dark]);
  const [hover, setHover] = useState<{ p: PeakCard; x: number; y: number } | null>(null);
  const [err, setErr] = useState("");
  const [ready, setReady] = useState(false);
  useEffect(() => {
    const el = ref.current as any;
    if (!el) return;
    let Plotly: any = null, alive = true, busy = false, ro: ResizeObserver | null = null;
    loadPlotly().then(P => {
      if (!alive) return;
      Plotly = P;
      const layout = { ...fig.layout, autosize: true, font: { ...(fig.layout.font ?? {}), family: SANS } };
      delete layout.width;
      const staticAnn = [...(layout.annotations ?? [])];
      const precursor = (layout.shapes ?? []).find((s: any) => s.line?.dash === "dot")?.x0;
      return Plotly.newPlot(el, fig.data, layout, CONFIG).then(() => {
        if (alive) setReady(true);
        // the plot follows its container (sidebar drag, drawer, window resize), not only the window
        if (alive && typeof ResizeObserver !== "undefined") {
          let frame = 0;
          ro = new ResizeObserver(() => { cancelAnimationFrame(frame); frame = requestAnimationFrame(() => { if (alive) Plotly.Plots.resize(el); }); });
          ro.observe(el);
        }
        if (!alive || !peaks?.length) return;
        const relabel = () => {
          busy = true;
          const anns = peakLabels(el, peaks, staticAnn, precursor);
          Plotly.relayout(el, { annotations: dark ? JSON.parse(darkColors(JSON.stringify(anns), true)) : anns }).finally(() => { busy = false; });
        };
        relabel();
        el.on("plotly_hover", (ev: any) => {
          const pt = (ev.points ?? []).find((q: any) => q.data.name === "peaks");
          if (pt) setHover({ p: peaks[pt.customdata], x: ev.event?.clientX ?? 0, y: ev.event?.clientY ?? 0 });
        });
        el.on("plotly_unhover", () => setHover(null));
        el.on("plotly_relayout", (e: any) => {
          const keys = Object.keys(e ?? {});
          if (busy || (keys.length && keys.every(k => k.startsWith("annotations")))) return;
          setHover(null);
          relabel();
        });
      });
    }).catch(e => alive && setErr(String(e)));
    return () => { alive = false; ro?.disconnect(); if (Plotly) Plotly.purge(el); };
  }, [fig, peaks]);   // fig changes with the theme too
  if (err) return <div className="out placeholder-box">Figure failed: {err}</div>;
  return (
    <div className="fig" onMouseLeave={() => setHover(null)}>
      {!ready && <div className="fig-loading">Loading figure…</div>}
      <div ref={ref} />
      {hover && <HoverCard {...hover} dec={mzDecimals(peaks ?? [])} />}
    </div>
  );
}

function HoverCard({ p, x, y, dec }: { p: PeakCard; x: number; y: number; dec: number }) {
  const W = 250, H = p.img ? 330 : 170;
  const left = x + 18 + W > window.innerWidth - 6 ? Math.max(6, x - 18 - W) : x + 18;
  const top = Math.max(6, Math.min(y - H / 2, window.innerHeight - H - 6));
  return (
    <div className="hover-card" style={{ left, top, width: W }}>
      <div className={"hc-kind" + (p.on ? "" : " off")}>{p.kind}</div>
      <div className="hc-mz">m/z {p.mz.toFixed(dec)}</div>
      <table><tbody>{p.rows.map(([l, v]) => <tr key={l}><td className="l">{l}</td><td className="v">{v}</td></tr>)}</tbody></table>
      {p.img && <img src={p.img} alt="fragment structure" />}
    </div>
  );
}

/** 3Dmol viewer. `className` lets the compare block make the canvas fill a square box. */
export function Mol3D({ molblock, className = "mol3d" }: { molblock: string; className?: string }) {
  const ref = useRef<HTMLDivElement>(null);
  const [err, setErr] = useState("");
  const { dark } = useTheme();
  useEffect(() => {
    let viewer: any = null, alive = true, ro: ResizeObserver | null = null;
    load3Dmol().then($3Dmol => {
      if (!alive || !ref.current) return;
      ref.current.innerHTML = "";
      viewer = $3Dmol.createViewer(ref.current, { backgroundColor: dark ? "#181B21" : "white" });
      viewer.addModel(molblock, "mol");
      viewer.setStyle({}, { stick: { radius: 0.14 }, sphere: { scale: 0.22 } });
      // a little margin: zoomTo() fits the atom centres, so terminal atoms would touch the box
      const fit = () => { viewer.zoomTo(); viewer.zoom(0.85); viewer.render(); };
      fit();
      // the compare cards follow the sidebar / window: 3Dmol does not resize its canvas by itself
      if (alive && typeof ResizeObserver !== "undefined") {
        ro = new ResizeObserver(() => { if (alive) { viewer.resize(); fit(); } });
        ro.observe(ref.current);
      }
    }).catch(e => alive && setErr(String(e)));
    return () => { alive = false; ro?.disconnect(); viewer?.clear?.(); };
  }, [molblock, dark]);
  if (err) return <div className="out placeholder-box">3D viewer failed: {err}</div>;
  return <div ref={ref} className={className} />;
}
