import { useEffect, useLayoutEffect, useMemo, useRef, useState } from "react";
import { createPortal } from "react-dom";
import { darkColors } from "../darkColors";
import { useReportMode } from "../report/mode";
import { useTheme } from "../themeMode";
import type { BlockData, PeakCard } from "../types";
import { hitPeak, hoverCardPosition, stickBaseline, type PeakHit } from "./peakHover";
import { isPrecursorPeak, labelEligible } from "./peakLabels";

const MONO = '"IBM Plex Mono",ui-monospace,"SF Mono",Menlo,Consolas,monospace';
const SANS = '"IBM Plex Sans",-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif';
const INK = "#161A23", MUTED = "#6A7180";     // peak labels; recoloured with the figure in the dark theme
const LEADER = "#9AA1AE";                      // leader line from a peak to a label set further out: thin and grey, never mistaken for a stick
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
export function ThemedSvg({ svg, className = "fig-img svg", alt = "", keep }: { svg: string; className?: string; alt?: string; keep?: string[] }) {
  const { dark } = useTheme();
  const keepKey = keep?.join(",") ?? "";          // `keep`: colours that stay as drawn in the dark theme (the fragment view's palette)
  const src = useMemo(() => `data:image/svg+xml;charset=utf-8,${encodeURIComponent(dark ? darkColors(svg, false, keep) : svg)}`, [svg, dark, keepKey]);
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

type Box = { l: number; r: number; t: number; b: number };
type Span = { i: number; x: number; lo: number; hi: number };      // a vertical line in px: a peak's stick (i) or a precursor line (-1)
type Precursor = { x: number; y0: number; y1: number };

/** The dotted precursor lines of a figure (ms2_plot.precursor_line), in data units. */
const precursorLines = (layout: any): Precursor[] =>
  (layout.shapes ?? []).filter((s: any) => s.line?.dash === "dot").map((s: any) => ({ x: s.x0, y0: s.y0, y1: s.y1 }));

/** A leader line (from the peak tip to its label) crosses a label box: sampled along the segment. */
const crosses = (x0: number, y0: number, x1: number, y1: number, o: Box) =>
  Array.from({ length: 9 }, (_, k) => k / 8).some(t => {
    const x = x0 + (x1 - x0) * t, y = y0 + (y1 - y0) * t;
    return x > o.l - 1 && x < o.r + 1 && y > o.t - 1 && y < o.b + 1;
  });

/**
 * m/z 标签：只标主峰（每半边 ≥ 5 % 基峰，见 peakLabels.ts）和前体峰；重要的峰优先（彩色碎片，已解释 / 匹配，再按强度）；标签不与别的标签重叠，也不压住别的峰和前体线。缩放后重算。
 * 峰的 rel 可以是负的（镜像图的下半）：标签写在峰尖下方。
 */
function peakLabels(gd: any, peaks: PeakCard[], staticAnn: any[], precursors: Precursor[], dark = false) {
  const { xaxis: xa, yaxis: ya } = gd._fullLayout;
  const x0 = Math.min(...xa.range), x1 = Math.max(...xa.range), yLo = Math.min(...ya.range), yHi = Math.max(...ya.range);
  const clampY = (y: number) => Math.max(yLo, Math.min(yHi, y));
  const vis = peaks.map((_, i) => i).filter(i => peaks[i].mz >= x0 && peaks[i].mz <= x1);
  const span = (i: number, x: number, ya0: number, ya1: number): Span => {
    const p0 = ya.l2p(clampY(ya0)), p1 = ya.l2p(clampY(ya1));
    return { i, x: xa.l2p(x), lo: Math.min(p0, p1), hi: Math.max(p0, p1) };
  };
  const lines = vis.map(i => span(i, peaks[i].mz, 0, peaks[i].rel));
  for (const pr of precursors) if (pr.x >= x0 && pr.x <= x1) lines.push(span(-1, pr.x, pr.y0, pr.y1));
  const dec = mzDecimals(peaks);
  const precursorXs = precursors.map(pr => pr.x);
  const eligible = labelEligible(peaks, vis, precursorXs);      // only the main peaks of each half and the precursor are written; hover cards stay on every peak
  const order = vis.filter(i => eligible.has(i)).sort((a, b) => (Number(!!peaks[b].color) - Number(!!peaks[a].color)) || (Number(peaks[b].on) - Number(peaks[a].on))
    || (Math.abs(peaks[b].rel) - Math.abs(peaks[a].rel)));
  const boxes: Box[] = [];
  const leaders: [number, number, number, number][] = [];
  const anns: any[] = [];
  const coloured = peaks.some(p => p.color);          // fragment view: only the fragments' peaks are labelled, each in its colour
  for (const i of order) {
    if (anns.length >= MAX_LABELS) break;
    const p = peaks[i];
    if (p.rel > yHi || p.rel < yLo || (coloured && !p.color)) continue;
    const down = p.rel < 0;                            // mirror plot, lower half: the label hangs below the tip
    const dir = down ? 1 : -1;                         // px direction away from the baseline
    // a peak with a `label` (fragment view) is written in its own colour on two lines: "m/z<br>Δ ppm"
    const txt = p.label ?? p.mz.toFixed(dec), rows = txt.split("<br>"), w = Math.max(...rows.map(l => l.length)) * 6.1 + 4, h = 12 * rows.length + 2;
    const px = xa.l2p(p.mz), py = ya.l2p(p.rel);
    const isPrecursor = isPrecursorPeak(p, precursorXs);          // its label may sit on the dotted precursor line
    // side: 0 centred over the peak, 1 / -1 beside it (right / left), 2 / 3 one or two label widths further out;
    // lift: px beyond the tip, joined to it by a thin leader line
    // strict 0: the leader passes no stick and no label; 1: it may cross a stick; 2: it may pass a label (a coloured fragment is never left unlabelled)
    const place = (side: number, lift: number, strict = 0): Box | null => {
      const out = (Math.abs(side) - 1) * (w + 6);
      const l = side === 0 ? px - w / 2 : side > 0 ? px + 3 + out : px - 3 - w - out;
      const near = py + dir * (2 + lift), far = near + dir * h;
      const b = { l, r: l + w, t: Math.min(near, far), b: Math.max(near, far) };
      if (b.l < 0 || b.r > xa._length || b.t < 0 || b.b > ya._length) return null;
      const other = (s: Span) => s.i !== i && !(s.i === -1 && isPrecursor);
      const cx = (b.l + b.r) / 2, ly = py + dir * 2;
      const hit = boxes.some(o => b.l < o.r + 4 && b.r > o.l - 4 && b.t < o.b + 1 && b.b > o.t - 1)
        || lines.some(s => other(s) && s.x > b.l - 1 && s.x < b.r + 1 && s.lo < b.b && s.hi > b.t)
        || leaders.some(([ax, ay, bx, by]) => crosses(ax, ay, bx, by, b))
        // the leader itself passes no other label and no other stick
        || (lift > 0 && ((strict < 2 && boxes.some(o => crosses(px, ly, cx, near, o)))
          || (strict < 1 && lines.some(s => {
            if (!other(s) || Math.abs(s.x - px) < 1 || (s.x - px) * (s.x - cx) > 0) return false;
            const y = ly + (near - ly) * (s.x - px) / (cx - px);
            return y > s.lo - 1 && y < s.hi + 1;
          }))));
      return hit ? null : b;
    };
    // a fragment's label that would sit on a taller neighbour goes beside its peak, or further out with a leader:
    // of the free places, the one with the shortest leader
    let lift = 0, b = place(0, 0);
    if (!b && p.color) {
      let best = Infinity;
      for (let strict = 0; strict <= 2 && !b; strict++) {
        for (let up = 0; up <= 168; up += 14) {
          for (const side of up === 0 ? [1, -1] : [0, 1, -1, 2, -2, 3, -3]) {
            const c = place(side, up, strict);
            const cost = c ? up + Math.abs((c.l + c.r) / 2 - px) : Infinity;
            if (cost < best) { best = cost; b = c; lift = up; }
          }
        }
      }
    }
    if (!b) continue;
    boxes.push(b);
    const cx = (b.l + b.r) / 2;
    if (lift) leaders.push([px, py + dir * 2, cx, py + dir * (2 + lift)]);
    // a fragment's label is written in its darker hue on a white page (`ink`), in its own hue on a dark one
    const font = { family: MONO, size: 10, color: p.color ? (dark ? p.color : p.ink ?? p.color) : p.on ? INK : MUTED };
    const base = { x: p.mz, y: p.rel, xref: "x", yref: "y", text: txt, yanchor: down ? "top" : "bottom", font };
    anns.push(lift
      ? { ...base, showarrow: true, arrowhead: 0, arrowwidth: 0.75, arrowcolor: LEADER, standoff: 1, ax: cx - px, ay: dir * (lift + 2), xanchor: "center" }
      : { ...base, showarrow: false, yshift: -dir * 2, xanchor: "center", xshift: cx - px });
  }
  return staticAnn.concat(anns);
}

function PlotlyFigure({ fig: raw, peaks }: { fig: { data: any[]; layout: any }; peaks?: PeakCard[] }) {
  const ref = useRef<HTMLDivElement>(null);
  const { dark } = useTheme();
  // the backend draws for a white page; in the dark theme every colour of the figure (and of the peak labels) is remapped
  const keep = useMemo(() => (peaks ?? []).map(p => p.color).filter((c): c is string => !!c), [peaks]);       // fragment colours are not remapped
  const fig = useMemo(() => (dark ? JSON.parse(darkColors(JSON.stringify(raw), true, keep)) : raw), [raw, dark, keep]);
  const [hover, setHover] = useState<{ p: PeakCard; x: number; y: number; hit: PeakHit } | null>(null);
  const [err, setErr] = useState("");
  const [ready, setReady] = useState(false);
  const report = useReportMode();
  useEffect(() => {
    const el = ref.current as any;
    if (!el) return;
    let Plotly: any = null, alive = true, busy = false, ro: ResizeObserver | null = null;
    let removeHover = () => {};
    let labelled: Promise<unknown> = Promise.resolve();      // the peak labels are a second drawing step
    setHover(null);
    const drawn = loadPlotly().then(P => {
      if (!alive) return;
      Plotly = P;
      const layout = { ...fig.layout, autosize: true, font: { ...(fig.layout.font ?? {}), family: SANS } };
      delete layout.width;
      if (peaks?.length) layout.hovermode = false;  // hit-test the whole stick instead of Plotly's invisible tip marker
      const staticAnn = [...(layout.annotations ?? [])];
      const precursors = precursorLines(layout);
      return Plotly.newPlot(el, fig.data, layout, CONFIG).then(() => {
        if (alive) setReady(true);
        // the plot follows its container (sidebar drag, drawer, window resize), not only the window
        if (alive && typeof ResizeObserver !== "undefined") {
          let frame = 0, width = el.clientWidth, height = el.clientHeight;
          ro = new ResizeObserver(() => {
            // The first observation is not a resize; don't erase a hover just as the plot becomes ready.
            if (width === el.clientWidth && height === el.clientHeight) return;
            width = el.clientWidth; height = el.clientHeight;
            cancelAnimationFrame(frame);
            frame = requestAnimationFrame(() => {
              if (alive) { setHover(null); if (el.clientWidth && el.clientHeight) Plotly.Plots.resize(el); }
            });
          });
          ro.observe(el);
        }
        if (!alive || !peaks?.length) return;
        const relabel = () => {
          busy = true;
          const anns = peakLabels(el, peaks, staticAnn, precursors, dark);
          labelled = Plotly.relayout(el, { annotations: dark ? JSON.parse(darkColors(JSON.stringify(anns), true, keep)) : anns }).finally(() => { busy = false; });
        };
        relabel();
        const clearHover = () => setHover(null);
        const baseline = stickBaseline(fig.data);
        const move = (event: MouseEvent) => {
          if (event.buttons || (event.target as Element).closest?.(".modebar")) { clearHover(); return; }
          const { xaxis, yaxis } = el._fullLayout;
          const rect = el.getBoundingClientRect();
          const hit = hitPeak(peaks, xaxis, yaxis, event.clientX - rect.left, event.clientY - rect.top, baseline);
          setHover(hit ? { p: peaks[hit.index], x: event.clientX, y: event.clientY, hit } : null);
        };
        el.addEventListener("mousemove", move);
        el.addEventListener("mouseleave", clearHover);
        el.addEventListener("mousedown", clearHover);
        window.addEventListener("scroll", clearHover, true);
        removeHover = () => {
          el.removeEventListener("mousemove", move);
          el.removeEventListener("mouseleave", clearHover);
          el.removeEventListener("mousedown", clearHover);
          window.removeEventListener("scroll", clearHover, true);
        };
        el.on("plotly_relayout", (e: any) => {
          const keys = Object.keys(e ?? {});
          if (busy || (keys.length && keys.every(k => k.startsWith("annotations")))) return;
          setHover(null);
          relabel();
        });
      });
    }).catch(e => alive && setErr(String(e)));
    report?.track(drawn.then(() => labelled));                // the Session report waits until the plot and its labels are drawn
    return () => { alive = false; removeHover(); ro?.disconnect(); if (Plotly) Plotly.purge(el); };
  }, [fig, peaks, keep]);   // fig changes with the theme too
  if (err) return <div className="out placeholder-box">Figure failed: {err}</div>;
  return (
    <div className="fig" onMouseLeave={() => setHover(null)}>
      {!ready && <div className="fig-loading">Loading figure…</div>}
      <div className="fig-plot">
        <div ref={ref} />
        {hover && <svg className="fig-peak-highlight" aria-hidden="true">
          <line x1={hover.hit.x} x2={hover.hit.x} y1={hover.hit.y0} y2={hover.hit.y1} className="fig-peak-halo" stroke={hover.p.color ?? "var(--accent)"} />
          <line x1={hover.hit.x} x2={hover.hit.x} y1={hover.hit.y0} y2={hover.hit.y1} className="fig-peak-line" stroke={hover.p.color ?? "var(--accent)"} />
          {hover.hit.tipVisible && <circle cx={hover.hit.x} cy={hover.hit.y1} r="3" fill={hover.p.color ?? "var(--accent)"} />}
        </svg>}
      </div>
      {hover && <HoverCard {...hover} dec={mzDecimals(peaks ?? [])} />}
    </div>
  );
}

function HoverCard({ p, x, y, dec }: { p: PeakCard; x: number; y: number; dec: number }) {
  const ref = useRef<HTMLDivElement>(null);
  useLayoutEffect(() => {
    const card = ref.current;
    if (!card) return;
    const position = () => {
      const { width, height } = card.getBoundingClientRect();
      const { left, top } = hoverCardPosition(x, y, width, height, window.innerWidth, window.innerHeight);
      card.style.left = `${left}px`; card.style.top = `${top}px`;
    };
    position();
    const observer = typeof ResizeObserver === "undefined" ? null : new ResizeObserver(position);
    observer?.observe(card);
    window.addEventListener("resize", position);
    return () => { observer?.disconnect(); window.removeEventListener("resize", position); };
  }, [p, x, y]);
  return createPortal(
    <div ref={ref} className="hover-card" role="tooltip">
      <div className={"hc-kind" + (p.on ? "" : " off")}>{p.kind}</div>
      <div className="hc-mz">m/z {p.mz.toFixed(dec)}</div>
      <table><tbody>{p.rows.map(([l, v]) => <tr key={l}><td className="l">{l}</td><td className="v">{v}</td></tr>)}</tbody></table>
      {p.img && <img src={p.img} alt="fragment structure" />}
    </div>, document.body
  );
}

/**
 * 3Dmol viewer with a reset button (top-right) that goes back to the fitted view: the fitted
 * orientation and scale, because the user can rotate, pan and zoom the model. `className` lets
 * the compare block make the canvas fill a square box.
 */
export function Mol3D({ molblock, className = "mol3d" }: { molblock: string; className?: string }) {
  const ref = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<any>(null);
  const homeRef = useRef<number[] | null>(null);
  const [err, setErr] = useState("");
  const [ready, setReady] = useState(false);
  const { dark } = useTheme();
  const report = useReportMode();
  useEffect(() => {
    let viewer: any = null, alive = true, ro: ResizeObserver | null = null;
    setReady(false);
    viewerRef.current = null;
    homeRef.current = null;
    const drawn = load3Dmol().then($3Dmol => {
      if (!alive || !ref.current) return;
      ref.current.innerHTML = "";
      viewer = $3Dmol.createViewer(ref.current, { backgroundColor: dark ? "#181B21" : "white" });
      viewer.addModel(molblock, "mol");
      viewer.setStyle({}, { stick: { radius: 0.14 }, sphere: { scale: 0.22 } });
      // a little margin: zoomTo() fits the atom centres, so terminal atoms would touch the box
      const fit = () => { viewer.zoomTo(); viewer.zoom(0.85); viewer.render(); homeRef.current = viewer.getView(); };
      fit();
      viewerRef.current = viewer;
      if (alive) setReady(true);
      // the compare cards follow the sidebar / window: 3Dmol does not resize its canvas by itself
      if (alive && typeof ResizeObserver !== "undefined") {
        ro = new ResizeObserver(() => { if (alive) { viewer.resize(); fit(); } });
        ro.observe(ref.current);
      }
    }).catch(e => alive && setErr(String(e)));
    report?.track(drawn);                                     // the Session report waits until the model is drawn
    return () => { alive = false; ro?.disconnect(); viewer?.clear?.(); viewerRef.current = null; homeRef.current = null; };
  }, [molblock, dark]);

  const reset = () => {
    const viewer = viewerRef.current;
    if (!viewer) return;
    if (homeRef.current) { viewer.setView(homeRef.current); viewer.render(); }
    else { viewer.zoomTo(); viewer.zoom(0.85); viewer.render(); }
  };

  if (err) return <div className="out placeholder-box">3D viewer failed: {err}</div>;
  return (
    <div className="mol3d-holder">
      <div ref={ref} className={className} />
      {ready && !report && <button type="button" className="mol3d-reset" onClick={reset} title="Reset the view to the fitted scale" aria-label="Reset the 3D view">
        <svg viewBox="0 0 16 16" width="13" height="13" aria-hidden="true" focusable="false">
          <path d="M13.2 8a5.2 5.2 0 1 1-1.5-3.7" fill="none" stroke="currentColor" strokeWidth="1.7" strokeLinecap="round" />
          <path d="M13.2 2.6v3h-3" fill="none" stroke="currentColor" strokeWidth="1.7" strokeLinecap="round" strokeLinejoin="round" />
        </svg>
      </button>}
    </div>
  );
}
