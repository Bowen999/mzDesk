import { createElement } from "react";
import { flushSync } from "react-dom";
import { createRoot } from "react-dom/client";
import iconsUrl from "../assets/icons.svg?url";
import { initialValues, summarizeInputs } from "../components/Fields";
import { resolvedTheme } from "../themeMode";
import type { Field, Spec } from "../types";
import { listRuns, type ReportEntry } from "./log";
import { MAX_REPORT_ROWS, type ReportMode } from "./mode";
import { ReportPage } from "./ReportPage";
import { edmonton, newRunId, REPORT_TZ } from "./time";
import { makeZip } from "./zip";

/**
 * Session report: one folder (<run id>/) with report.json (inputs + raw results of every page that has a result) and
 * report.html (the same pages, one self-contained file: styles, fonts, icons and images inlined).
 * Each page is drawn once more, off screen, from the result it logged, with the same components in report mode (report/mode.ts): folded
 * blocks are native <details>, tables list their rows instead of paging, long text is whole, and controls that need the app are left out.
 * That second drawing is a separate React root, so the pages on screen keep their open blocks, table pages and filters, and nothing flashes.
 * The folder is delivered as <run id>.zip, because a browser page cannot write a folder to disk.
 */
export type BuiltReport = { runId: string; filename: string; blob: Blob; pages: number };

const esc = (s: unknown) => String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");

function toBase64(buf: ArrayBuffer) {
  const bytes = new Uint8Array(buf);
  let bin = "";
  for (let i = 0; i < bytes.length; i += 0x8000) bin += String.fromCharCode(...bytes.subarray(i, i + 0x8000));
  return btoa(bin);
}

async function toDataUrl(url: string): Promise<string | null> {
  try {
    const r = await fetch(url);
    if (!r.ok) return null;
    const blob = await r.blob();
    const buf = await blob.arrayBuffer();
    return `data:${blob.type || "application/octet-stream"};base64,${toBase64(buf)}`;
  } catch { return null; }
}

/** True when a unicode-range covers "A" (U+0041): the browser normalises "U+0000-00FF" to "U+0-FF", so read the ranges, not the text. */
function coversLatin(range: string) {
  return range.split(",").some(tok => {
    const m = /U\+([0-9A-F?]+)(?:-([0-9A-F]+))?/i.exec(tok.trim());
    if (!m) return false;
    const lo = parseInt(m[1].replace(/\?/g, "0"), 16), hi = m[2] ? parseInt(m[2], 16) : m[1].includes("?") ? parseInt(m[1].replace(/\?/g, "F"), 16) : lo;
    return lo <= 0x41 && 0x41 <= hi;
  });
}

/** A Latin-subset woff2 @font-face becomes a data URL; other subsets are left out (the report is read in English). */
async function fontFace(rule: CSSFontFaceRule, base: string): Promise<string> {
  const st = rule.style, range = st.getPropertyValue("unicode-range");
  if (range && !coversLatin(range)) return "";
  const m = /url\(["']?([^"')]+)["']?\)\s*format\(["']?woff2["']?\)/.exec(st.getPropertyValue("src"));
  if (!m) return "";
  const data = await toDataUrl(new URL(m[1], base).href);
  if (!data) return "";
  const d = (p: string) => st.getPropertyValue(p);
  return `@font-face{font-family:${d("font-family")};font-style:${d("font-style") || "normal"};font-weight:${d("font-weight") || "400"};font-display:swap;src:url(${data}) format("woff2");}`;
}

/** Every stylesheet of the page as text (the app's own, Plotly's), fonts inlined. */
async function collectCss(): Promise<string> {
  const out: string[] = [];
  for (const sheet of Array.from(document.styleSheets)) {
    let rules: CSSRuleList;
    try { rules = sheet.cssRules; } catch { continue; }       // a cross-origin sheet cannot be read
    for (const rule of Array.from(rules)) {
      if (rule instanceof CSSFontFaceRule) out.push(await fontFace(rule, sheet.href ?? document.baseURI));
      else out.push(rule.cssText);
    }
  }
  return out.filter(Boolean).join("\n");
}

/** The icon sprite (icons.svg) is an external <use> reference; the report carries its symbols inline. */
async function spriteSymbols(): Promise<string> {
  try {
    const text = await (await fetch(iconsUrl)).text();
    return (/<svg[^>]*>([\s\S]*)<\/svg>/.exec(text)?.[1] ?? "").trim();
  } catch { return ""; }
}

const STRIP = [":scope > .panel", ":scope > .run-progress", ":scope > .empty-hint", ":scope > .alert.error", ".modebar-container", ".hp-btn", ".copy", "script", "dialog"];

/** Copy of a drawn page, minus the controls that do nothing in a file. */
async function clonePage(src: HTMLElement): Promise<string> {
  const copy = src.cloneNode(true) as HTMLElement;
  for (const sel of STRIP) copy.querySelectorAll(sel).forEach(n => n.remove());
  // the hover previews (names, SMILES) need the app; in the file they are plain text
  copy.querySelectorAll("button.cell-depict").forEach(b => b.replaceWith(...Array.from(b.childNodes)));
  // 3Dmol and other canvases do not survive cloneNode: freeze them as images (an unreadable or empty canvas is dropped)
  const from = src.querySelectorAll("canvas"), to = copy.querySelectorAll("canvas");
  to.forEach((c, i) => {
    let url = "";
    try { url = (from[i] as HTMLCanvasElement).toDataURL("image/png"); } catch { /* tainted */ }
    if (url.length < 2000) { c.remove(); return; }
    const img = document.createElement("img");
    img.src = url; img.className = c.className; img.setAttribute("style", c.getAttribute("style") ?? "");
    img.style.maxWidth = "100%";
    c.replaceWith(img);
  });
  for (const img of Array.from(copy.querySelectorAll("img"))) {
    const s = img.getAttribute("src") ?? "";
    if (s && !s.startsWith("data:")) { const d = await toDataUrl(new URL(s, document.baseURI).href); if (d) img.setAttribute("src", d); }
  }
  copy.querySelectorAll("use").forEach(u => {
    const h = u.getAttribute("href") ?? "";
    const i = h.indexOf("#");
    if (i >= 0) u.setAttribute("href", h.slice(i));
  });
  return copy.outerHTML;
}

const sleep = (ms: number) => new Promise<void>(resolve => setTimeout(resolve, ms));
const FIGURE_WAIT_MS = 20000;        // longest wait for the figures of one page

/** A page is as wide here as the page area on screen (at most the width the report gives it), so a figure is drawn as wide as the person saw it. */
const pageWidth = () => Math.min(1240, Math.max(720, document.querySelector<HTMLElement>("main.main")?.clientWidth || 1100));

/** Draws the page in report mode in a hidden box outside the window, waits for its figures, and returns its HTML. */
async function renderPage(spec: Spec, e: ReportEntry, n: number): Promise<string | null> {
  const host = document.createElement("div");
  host.className = "main";
  host.setAttribute("aria-hidden", "true");
  host.style.cssText = `position:fixed;left:-20000px;top:0;width:${pageWidth()}px;visibility:hidden;pointer-events:none`;
  document.body.appendChild(host);
  const pending: Promise<unknown>[] = [];
  const mode: ReportMode = { maxRows: MAX_REPORT_ROWS, track: p => { pending.push(p); } };
  const root = createRoot(host, { identifierPrefix: `rp${n}-` });      // the ids of this drawing never meet those of the live pages
  try {
    flushSync(() => root.render(createElement(ReportPage, { spec, entry: e, mode })));
    // Plotly and 3Dmol draw after the first render: wait for them, then for what follows (peak labels, the "Loading" note, text fitted to its box)
    const timeout = sleep(FIGURE_WAIT_MS);
    for (let seen = 0; seen < pending.length;) {
      const batch = pending.slice(seen);
      seen = pending.length;
      await Promise.race([Promise.allSettled(batch), timeout]);
    }
    for (let i = 0; i < 20 && host.querySelector(".fig-loading"); i++) await sleep(50);
    await sleep(100);
    const page = host.firstElementChild as HTMLElement | null;
    return page ? await clonePage(page) : null;
  } finally {
    root.unmount();
    host.remove();
  }
}

/** Long strings (spectra, files) are shortened in the HTML view; report.json keeps them whole. */
const shorten = (_k: string, v: unknown) =>
  typeof v === "string" && v.length > 1500 ? `${v.slice(0, 1500)}… (+${v.length - 1500} characters, full text in report.json)` : v;

function fieldsOf(spec: Spec, e: ReportEntry): Field[] {
  return (e.kind === "module" ? spec.modules.find(m => m.id === e.id)?.inputs : spec.pipelines.find(p => p.id === e.id)?.inputs) ?? [];
}

function inputsHtml(spec: Spec, e: ReportEntry) {
  const fields = fieldsOf(spec, e);
  let chips = "";
  try {
    const s = summarizeInputs(fields, { ...initialValues(fields), ...e.inputs });
    chips = [s.lead, ...s.items].filter(Boolean).map(t => `<span class="rp-chip">${esc(t)}</span>`).join("");
  } catch { /* the raw JSON below still shows the inputs */ }
  return `<div class="rp-inputs">${chips}<details><summary>Inputs as sent (JSON)</summary><pre>${esc(JSON.stringify(e.inputs, shorten, 2))}</pre></details></div>`;
}

const REPORT_CSS = `
*,*::before,*::after{animation:none!important;transition:none!important}
.rp-wrap{max-width:1240px;margin:0 auto;padding:clamp(20px,3vw,40px) clamp(16px,3.6vw,52px) 64px}
.rp-head h1{margin:0 0 6px;font-size:26px;letter-spacing:-.02em}
.rp-meta{display:grid;grid-template-columns:max-content 1fr;gap:3px 14px;margin:14px 0 6px;font-size:13px;color:var(--muted)}
.rp-meta b{color:var(--ink);font-weight:500}
.rp-toc{margin:18px 0 8px;padding:12px 16px;background:var(--canvas);border:1px solid var(--line);border-radius:10px}
.rp-toc ol{margin:6px 0 0;padding-left:20px;line-height:1.7}
.rp-toc .rp-when{color:var(--muted);font-size:12px;margin-left:8px}
.rp-entry{margin-top:36px;padding-top:22px;border-top:2px solid var(--line-strong)}
.rp-entry-head{display:flex;flex-wrap:wrap;align-items:baseline;gap:4px 14px;margin-bottom:8px;color:var(--muted);font-size:12.5px}
.rp-entry-head .rp-no{color:var(--accent-ink);font-weight:600}
.rp-chip{display:inline-block;margin:0 6px 6px 0;padding:2px 9px;background:var(--canvas);border:1px solid var(--line);border-radius:99px;font-size:12px}
.rp-inputs{margin:0 0 10px}.rp-inputs pre{max-height:360px;overflow:auto;padding:10px 12px;background:var(--canvas);border:1px solid var(--line);border-radius:8px;font-size:12px}
.rp-entry .page{padding:6px 0 0;max-width:none}
.rp-entry details.gt-group>summary{list-style:none;cursor:pointer}
.rp-entry details.gt-group>summary::-webkit-details-marker{display:none}
.rp-entry details.gt-group[open] .gt-caret{transform:rotate(45deg);margin-bottom:4px}
.rp-entry .cand-row,.rp-entry div.gm-row{cursor:auto}
.rp-entry div.gm-row{box-sizing:border-box}
.rp-entry .cand-table .cand-detail>td{padding-left:12px}
.rp-foot{margin-top:40px;color:var(--faint);font-size:12px}
@media print{.rp-entry{break-before:page}.rp-entry:first-of-type{break-before:auto}}
`;

function buildHtml(spec: Spec, runId: string, when: ReturnType<typeof edmonton>, entries: ReportEntry[], pages: (string | null)[], css: string, sprite: string, theme: string) {
  const toc = entries.map((e, i) => `<li><a href="#rp-${i + 1}">${esc(e.title)}</a><span class="rp-when">${esc(e.section)} · ${esc(e.finishedAt.replace("T", " "))}</span></li>`).join("");
  const body = entries.map((e, i) => `
<section class="rp-entry" id="rp-${i + 1}">
  <div class="rp-entry-head"><span class="rp-no">${i + 1} / ${entries.length}</span><span>${esc(e.section)}</span><span>finished ${esc(e.finishedAt.replace("T", " "))}</span>${e.seconds != null ? `<span>${e.seconds.toFixed(1)} s</span>` : ""}</div>
  ${inputsHtml(spec, e)}
  ${pages[i] ?? `<p class="empty-hint">This page could not be drawn for the report; its result is in report.json.</p>`}
</section>`).join("\n");
  return `<!doctype html>
<html lang="en" data-theme="${theme}">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>mzDesk report ${esc(runId)}</title>
<style>${css}</style>
<style>${REPORT_CSS}</style>
</head>
<body>
<svg width="0" height="0" style="position:absolute" aria-hidden="true">${sprite}</svg>
<div class="rp-wrap">
  <header class="rp-head">
    <h1>mzDesk report</h1>
    <div class="rp-meta">
      <span>Run ID</span><b>${esc(runId)}</b>
      <span>Generated</span><b>${esc(when.human)}</b>
      <span>Time zone</span><b>${REPORT_TZ} (Edmonton time)</b>
      <span>mzDesk</span><b>V${esc(spec.version)}</b>
      <span>Pages</span><b>${entries.length}</b>
    </div>
    <nav class="rp-toc" aria-label="Pages in this report"><b>Pages</b><ol>${toc}</ol></nav>
  </header>
  <main class="main">${body}</main>
  <p class="rp-foot">Made by mzDesk in the browser from the results on screen. report.json in the same folder has the inputs and the raw results.</p>
</div>
</body>
</html>
`;
}

export async function buildReport(spec: Spec): Promise<BuiltReport | null> {
  const entries = listRuns();
  if (!entries.length) return null;
  const when = edmonton(), runId = newRunId(when);
  const theme = resolvedTheme();
  const [css, sprite] = await Promise.all([collectCss(), spriteSymbols()]);
  const pages: (string | null)[] = [];
  for (const [i, e] of entries.entries()) {
    try { pages.push(await renderPage(spec, e, i + 1)); }
    catch (err) { console.error("Session report: page not drawn", e.key, err); pages.push(null); }      // one page that fails does not stop the others
  }
  const json = {
    report_format: 1,
    run_id: runId,
    generated_at: when.iso,
    timezone: REPORT_TZ,
    timezone_abbreviation: when.abbr,
    utc_offset: when.offset,
    app: { name: "mzDesk", version: spec.version },
    pages: entries.map(e => ({
      kind: e.kind, id: e.id, title: e.title, section: e.section, finished_at: e.finishedAt, seconds: e.seconds,
      inputs: e.inputs, result: e.result, ...(e.steps !== undefined ? { steps: e.steps } : {}),
    })),
  };
  const html = buildHtml(spec, runId, when, entries, pages, css, sprite, theme);
  const blob = makeZip([
    { path: `${runId}/report.json`, data: JSON.stringify(json, null, 2) },
    { path: `${runId}/report.html`, data: html },
  ]);
  return { runId, filename: `${runId}.zip`, blob, pages: entries.length };
}
