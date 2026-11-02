import { BASE } from "../api";

/** Figures of the Docs: hand-written SVG files (mzdesk/docs/img/<name>.svg, served by /api/docs-img/<name>) drawn inline, so the theme tokens
 *  (CSS variables) and the manual light / dark switch apply to them. The text is sanitized before it is injected. */

// Elements a figure may use (lower case). Everything else, with its content, is dropped: script, style, foreignObject, image, a, animate, ...
const ALLOWED = new Set(["svg", "g", "defs", "title", "desc", "path", "rect", "circle", "ellipse", "line", "polyline", "polygon", "text", "tspan",
  "marker", "symbol", "use", "clippath", "mask", "lineargradient", "radialgradient", "stop", "pattern"]);
// <script> and <style> hold raw text: skip to their end tag instead of reading the content as markup.
const RAW_TEXT = new Set(["script", "style"]);

const TAG = /<(\/?)([A-Za-z][\w:.-]*)((?:\s+[^\s"'<>\/=]+(?:\s*=\s*(?:"[^"]*"|'[^']*'|[^\s"'=<>`]+))?)*)\s*(\/?)>/y;
const ATTR = /([^\s"'<>\/=]+)(?:\s*=\s*(?:"([^"]*)"|'([^']*)'|([^\s"'=<>`]+)))?/g;

const escAttr = (v: string) => v.replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");

/** One attribute, or null when it is dropped: on* handlers, style, xmlns, a bad name, javascript: / data: values, a link that is not "#id", a url() that is not "#id".
 *  Ids (and the references to them) get `prefix`, so two figures on one page never share an id. */
function cleanAttr(tag: string, name: string, value: string, prefix: string): string | null {
  const key = name.toLowerCase();
  if (!/^[a-z_][\w:.-]*$/.test(key) || key.startsWith("on") || key === "style" || key.startsWith("xmlns")) return null;
  const flat = value.toLowerCase().replace(/[\s\u0000-\u001f]+/g, "");
  if (/(?:javascript|vbscript|data):/.test(flat)) return null;
  if (key === "href" || key === "xlink:href") return value.startsWith("#") && /^#[\w.:-]+$/.test(value) ? `${name}="#${prefix}${value.slice(1)}"` : null;
  if (key === "id") return /^[\w.:-]+$/.test(value) ? `id="${prefix}${value}"` : null;
  if (key === "class") return /^[\w\s-]*$/.test(value) ? `class="${value.trim()}"` : null;
  if (/url\(/i.test(value)) {
    const ref = /^\s*url\(\s*#([\w.:-]+)\s*\)\s*$/.exec(value);
    return ref ? `${name}="url(#${prefix}${ref[1]})"` : null;
  }
  if (tag === "svg" && (key === "width" || key === "height")) return null;      // the page sets the size (CSS), the viewBox sets the shape
  return `${name}="${escAttr(value)}"`;
}

/** Keeps only what a figure needs. Returns the cleaned `<svg>…</svg>` markup, or null when the text is not an SVG (it must start with `<svg`).
 *  Every emitted token is checked on its own (tag allow-list, attribute rules, text with `<` escaped), so the result is safe whatever the input nests. */
export function sanitizeSvg(raw: string, prefix = "fg-"): string | null {
  const src = raw.trim();
  if (!/^<svg[\s/>]/i.test(src)) return null;
  const lower = src.toLowerCase();
  const out: string[] = [];
  const open: string[] = [];          // allowed elements that are open (names as written)
  let skip = 0, i = 0, done = false;  // skip > 0: inside a dropped element
  while (i < src.length && !done) {
    const lt = src.indexOf("<", i);
    if (lt < 0) { if (!skip && open.length) out.push(src.slice(i)); break; }
    if (lt > i && !skip && open.length) out.push(src.slice(i, lt));
    i = lt;
    if (src.startsWith("<!--", i)) { const e = src.indexOf("-->", i + 4); i = e < 0 ? src.length : e + 3; continue; }
    if (src.startsWith("<![CDATA[", i)) {
      const e = src.indexOf("]]>", i + 9), body = src.slice(i + 9, e < 0 ? src.length : e);
      if (!skip && open.length) out.push(body.replace(/&/g, "&amp;").replace(/</g, "&lt;"));
      i = e < 0 ? src.length : e + 3;
      continue;
    }
    if (src.startsWith("<?", i) || src.startsWith("<!", i)) { const e = src.indexOf(">", i); i = e < 0 ? src.length : e + 1; continue; }
    TAG.lastIndex = i;
    const m = TAG.exec(src);
    if (!m) { if (!skip && open.length) out.push("&lt;"); i += 1; continue; }
    i = TAG.lastIndex;
    const [, close, rawName, attrs, selfClose] = m;
    const name = rawName.toLowerCase();
    if (!close && RAW_TEXT.has(name) && !selfClose) {
      const e = lower.indexOf(`</${name}`, i), g = e < 0 ? -1 : src.indexOf(">", e);
      i = g < 0 ? src.length : g + 1;
      continue;
    }
    if (skip) { if (close) skip--; else if (!selfClose) skip++; continue; }
    if (!ALLOWED.has(name)) { if (!close && !selfClose) skip = 1; continue; }
    if (close) {
      const at = open.map(o => o.toLowerCase()).lastIndexOf(name);
      if (at < 0) continue;                                  // a closing tag without its opening tag
      while (open.length > at) out.push(`</${open.pop()}>`);
      if (!open.length) done = true;
      continue;
    }
    if (!open.length && name !== "svg") continue;            // the root is the first element
    const kept: string[] = [];
    for (const a of attrs.matchAll(ATTR)) {
      const clean = cleanAttr(name, a[1], a[2] ?? a[3] ?? a[4] ?? "", prefix);
      if (clean) kept.push(clean);
    }
    if (!open.length) kept.push('aria-hidden="true"', 'focusable="false"');   // the wrapper carries role="img" and the label
    out.push(`<${rawName}${kept.length ? " " + kept.join(" ") : ""}${selfClose ? "/>" : ">"}`);
    if (selfClose) { if (!open.length) done = true; } else open.push(rawName);
  }
  while (open.length) out.push(`</${open.pop()}>`);
  const svg = out.join("");
  return /^<svg[\s/>]/i.test(svg) ? svg : null;
}

// Loaded figures, by name. null = the file is missing or is not an SVG. A network failure is not cached, so a later view tries again.
const loaded = new Map<string, string | null>();
const loading = new Map<string, Promise<string | null>>();

/** The sanitized markup of a figure that was loaded before; undefined when it has not been loaded (yet). */
export const cachedFigure = (name: string) => loaded.get(name);

export function loadFigure(name: string): Promise<string | null> {
  if (loaded.has(name)) return Promise.resolve(loaded.get(name)!);
  let p = loading.get(name);
  if (!p) {
    p = fetch(`${BASE}api/docs-img/${name}`)
      .then(async r => {
        const svg = r.ok ? sanitizeSvg(await r.text(), `fg-${name}-`) : null;
        loaded.set(name, svg);
        return svg;
      })
      .finally(() => loading.delete(name));
    loading.set(name, p);
  }
  return p;
}

/** For tests: forget every loaded figure. */
export const clearFigures = () => { loaded.clear(); loading.clear(); };
