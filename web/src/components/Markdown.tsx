import { useEffect, useRef, useState, type ReactNode } from "react";
import { Link } from "react-router-dom";
import { cachedFigure, loadFigure } from "../docs/figure";
import { fenceClose, fenceOpen, uniqueSlug } from "../docs/outline";
import { RichText } from "./Help";

/** The markdown subset of the Docs files (mzdesk/docs/*.md). Nothing else is interpreted.
 *  Blocks: #, ##, ### headings (#### has no anchor) · paragraphs · "- " / "1. " lists (nested by indent) · ```lang code blocks (a Copy button for
 *  bash, sh, shell, zsh, powershell, pwsh, cmd, console, text) · | pipe | tables · "> **Note** text" and "> **Warning** text" callouts, other quotes ·
 *  --- rules · "![Alt](fig:name)" figures with an optional next line "*Caption*".
 *  Inline: **bold**, *italic*, `code`, [text](/docs/id#anchor) (app link), [text](https://…) (new tab), [text](#anchor), [text](mailto:…),
 *  and the <sub> <sup> <i> <b> tags that RichText draws. */

// ------------------------------------------------------------------ inline
// One pass, leftmost match wins; at the same place code, link, bold, italic in that order.
const INLINE = /`([^`]+)`|\[([^\]]+)\]\(((?:[^()\s]|\([^()\s]*\))+)\)|\*\*(.+?)\*\*|\*([^*\s](?:[^*]*[^*\s])?)\*/g;

function link(label: ReactNode, href: string, key: number): ReactNode {
  if (/^https?:\/\//i.test(href)) return <a key={key} href={href} target="_blank" rel="noopener noreferrer">{label}</a>;
  if (/^\/(docs|modules|pipelines)(?:[/#]|$)/.test(href)) return <Link key={key} to={href}>{label}</Link>;      // pages of the app
  if (href.startsWith("#") || /^mailto:/i.test(href) || (href.startsWith("/") && !href.startsWith("//"))) return <a key={key} href={href}>{label}</a>;
  return <span key={key}>{label}</span>;                            // any other address (javascript:, data:, a relative file) is not a link
}

function inline(text: string): ReactNode {
  const re = new RegExp(INLINE.source, "g");
  const parts: ReactNode[] = [];
  let last = 0, k = 0, m: RegExpExecArray | null;
  const flushText = (to: number) => { if (to > last) parts.push(<RichText key={k++} text={text.slice(last, to)} />); };
  while ((m = re.exec(text))) {
    flushText(m.index);
    if (m[1] !== undefined) parts.push(<code key={k++}>{m[1]}</code>);
    else if (m[2] !== undefined) parts.push(link(inline(m[2]), m[3], k++));
    else if (m[4] !== undefined) parts.push(<b key={k++}>{inline(m[4])}</b>);
    else parts.push(<i key={k++}>{inline(m[5])}</i>);
    last = m.index + m[0].length;
  }
  flushText(text.length);
  return parts.length === 1 ? parts[0] : <>{parts}</>;
}

// ------------------------------------------------------------------ blocks
type Item = { text: string; sub: Block[] };
type Block =
  | { t: "h"; level: 1 | 2 | 3 | 4; text: string; id?: string }
  | { t: "p"; text: string }
  | { t: "list"; ordered: boolean; start: number; items: Item[] }
  | { t: "code"; lang: string; code: string }
  | { t: "table"; head: string[]; rows: string[][] }
  | { t: "quote"; kind: "Note" | "Warning" | null; blocks: Block[] }
  | { t: "hr" }
  | { t: "figure"; alt: string; name: string; caption?: string };

const LIST_ITEM = /^(\s*)([-*+]|\d{1,3}[.)])\s+(.*)$/;
const TABLE_SEP = /^\s*\|?\s*:?-+:?\s*(\|\s*:?-+:?\s*)*\|?\s*$/;
const FIGURE = /^!\[([^\]]*)\]\(fig:([a-z0-9-]+)\)\s*$/;
const CAPTION = /^\*([^*].*[^*]|[^*])\*$/;
const RULE = /^\s*([-*_])(\s*\1){2,}\s*$/;
const COPYABLE = new Set(["bash", "sh", "shell", "zsh", "powershell", "pwsh", "ps1", "cmd", "console", "text"]);

/** The cells of a pipe-table row (`\|` is a literal pipe). */
function cells(line: string): string[] {
  const body = line.trim().replace(/^\|/, "").replace(/(?<!\\)\|$/, "");
  return body.split(/(?<!\\)\|/).map(c => c.replace(/\\\|/g, "|").trim());
}

/** Lines → blocks. `seen` makes the heading anchors (the same ids as headingsOf() in docs/outline.ts); inside a quote it is null: no anchors there. */
function parse(lines: string[], seen: Map<string, number> | null): Block[] {
  const blocks: Block[] = [];
  let para: string[] = [];
  const flush = () => { if (para.length) { blocks.push({ t: "p", text: para.join(" ") }); para = []; } };

  for (let i = 0; i < lines.length; i++) {
    const line = lines[i].trimEnd();

    const fence = fenceOpen(line);
    if (fence) {
      flush();
      const code: string[] = [];
      for (i++; i < lines.length && !fenceClose(lines[i], fence.marker); i++) code.push(lines[i].replace(new RegExp(`^ {0,${fence.indent}}`), ""));
      blocks.push({ t: "code", lang: fence.lang, code: code.join("\n") });
      continue;
    }

    const head = /^(#{1,4})\s+(.*)$/.exec(line);
    if (head) {
      flush();
      const level = head[1].length as 1 | 2 | 3 | 4;
      blocks.push({ t: "h", level, text: head[2], id: seen && level < 4 ? uniqueSlug(head[2], seen) : undefined });
      continue;
    }

    if (!line.trim()) { flush(); continue; }

    const fig = FIGURE.exec(line.trim());
    if (fig) {
      flush();
      const cap = CAPTION.exec((lines[i + 1] ?? "").trim());
      blocks.push({ t: "figure", alt: fig[1], name: fig[2], caption: cap ? cap[0].slice(1, -1) : undefined });
      if (cap) i++;
      continue;
    }

    if (RULE.test(line) && !LIST_ITEM.test(line)) { flush(); blocks.push({ t: "hr" }); continue; }

    if (/^\s*>/.test(line)) {
      flush();
      const quote: string[] = [];
      for (; i < lines.length && /^\s*>/.test(lines[i]); i++) quote.push(lines[i].replace(/^\s*>\s?/, ""));
      i--;
      const kind = /^\*\*(Note|Warning):?\*\*:?\s*/.exec(quote[0]);
      if (kind) quote[0] = quote[0].slice(kind[0].length);
      blocks.push({ t: "quote", kind: kind ? (kind[1] as "Note" | "Warning") : null, blocks: parse(quote, null) });
      continue;
    }

    if (line.includes("|") && i + 1 < lines.length && TABLE_SEP.test(lines[i + 1]) && cells(lines[i + 1]).length === cells(line).length) {
      flush();
      const head = cells(line), rows: string[][] = [];
      for (i += 2; i < lines.length && lines[i].trim() && lines[i].includes("|"); i++) rows.push(cells(lines[i]));
      i--;
      blocks.push({ t: "table", head, rows: rows.map(r => head.map((_, c) => r[c] ?? "")) });
      continue;
    }

    if (LIST_ITEM.test(line)) {
      flush();
      const items: { indent: number; ordered: boolean; num: number; text: string }[] = [];
      for (; i < lines.length; i++) {
        const l = lines[i].trimEnd();
        const m = LIST_ITEM.exec(l);
        if (m && !RULE.test(l)) { items.push({ indent: m[1].length, ordered: /\d/.test(m[2]), num: parseInt(m[2], 10), text: m[3] }); continue; }
        if (!l.trim()) {                                           // a blank line keeps the list going when the next line is another item
          let j = i + 1;
          while (j < lines.length && !lines[j].trim()) j++;
          if (j < lines.length && LIST_ITEM.test(lines[j]) && !RULE.test(lines[j])) { i = j - 1; continue; }
          break;
        }
        if (/^\s{2,}\S/.test(l) && !fenceOpen(l) && items.length) { items[items.length - 1].text += " " + l.trim(); continue; }   // wrapped item text
        break;
      }
      i--;
      let at = 0;
      const build = (): Block => {
        const first = items[at], list: Extract<Block, { t: "list" }> = { t: "list", ordered: first.ordered, start: first.ordered ? first.num : 1, items: [] };
        while (at < items.length && items[at].indent >= first.indent) {
          const it = items[at];
          if (it.indent > first.indent || it.ordered !== first.ordered) break;
          at++;
          const sub: Block[] = [];
          while (at < items.length && items[at].indent > it.indent) sub.push(build());
          list.items.push({ text: it.text, sub });
        }
        return list;
      };
      while (at < items.length) blocks.push(build());
      continue;
    }

    para.push(line.trim());
  }
  flush();
  return blocks;
}

// ------------------------------------------------------------------ drawing
/** Copy to the clipboard: the async API, or a hidden textarea where that is not available (plain http, old browsers). */
async function copyText(text: string): Promise<boolean> {
  try { await navigator.clipboard.writeText(text); return true; } catch { /* fall back */ }
  try {
    const ta = document.createElement("textarea");
    ta.value = text;
    ta.setAttribute("readonly", "");
    ta.style.cssText = "position:fixed;top:0;left:0;opacity:0";
    document.body.appendChild(ta);
    ta.select();
    const ok = document.execCommand("copy");
    document.body.removeChild(ta);
    return ok;
  } catch { return false; }
}

function CodeBlock({ lang, code }: { lang: string; code: string }) {
  const [state, setState] = useState<"idle" | "copied" | "failed">("idle");
  const timer = useRef<number | undefined>(undefined);
  const codeEl = useRef<HTMLElement>(null);
  useEffect(() => () => window.clearTimeout(timer.current), []);
  const copy = async () => {
    const ok = await copyText(code);
    if (!ok && codeEl.current) window.getSelection()?.selectAllChildren(codeEl.current);       // the reader can press Ctrl+C / Cmd+C
    setState(ok ? "copied" : "failed");
    window.clearTimeout(timer.current);
    timer.current = window.setTimeout(() => setState("idle"), 1800);
  };
  return (
    <div className={"doc-code" + (COPYABLE.has(lang) ? " has-copy" : "")}>
      <pre tabIndex={0}><code ref={codeEl} className={lang ? `language-${lang}` : undefined}>{code}</code></pre>
      {COPYABLE.has(lang) && (
        <button type="button" className={"doc-copy" + (state === "copied" ? " done" : "")} onClick={copy} aria-label={state === "copied" ? "Copied" : "Copy code"}>
          {state === "copied" ? "Copied" : state === "failed" ? "Press Ctrl+C" : "Copy"}
        </button>
      )}
      <span className="sr-only" role="status">{state === "copied" ? "Copied to clipboard" : state === "failed" ? "Selected; press Ctrl+C to copy" : ""}</span>
    </div>
  );
}

/** A figure: the SVG file drawn inline (see docs/figure.ts). A skeleton shows until it is there; a missing file leaves the description. */
function DocFigure({ name, alt, caption }: { name: string; alt: string; caption?: string }) {
  const [svg, setSvg] = useState<string | null | undefined>(() => cachedFigure(name));
  useEffect(() => {
    if (svg !== undefined) return;
    let live = true;
    loadFigure(name).then(s => { if (live) setSvg(s); }).catch(() => { if (live) setSvg(null); });
    return () => { live = false; };
  }, [name, svg]);
  return (
    <figure className="doc-figure">
      {svg
        ? <div className="doc-fig-scroll" role="img" aria-label={alt} tabIndex={0} dangerouslySetInnerHTML={{ __html: svg }} />
        : svg === null ? <p className="doc-fig-missing">{alt}</p>
        : <div className="doc-fig-skel" aria-hidden="true"><i /></div>}
      {caption && <figcaption>{inline(caption)}</figcaption>}
    </figure>
  );
}

function draw(b: Block, key: number): ReactNode {
  switch (b.t) {
    case "h": {
      const Tag = (`h${b.level}`) as "h1" | "h2" | "h3" | "h4";
      return <Tag key={key} id={b.id}>{inline(b.text)}</Tag>;
    }
    case "p": return <p key={key}>{inline(b.text)}</p>;
    case "list": {
      const Tag = b.ordered ? "ol" : "ul";
      return (
        <Tag key={key} start={b.ordered && b.start !== 1 ? b.start : undefined}>
          {b.items.map((it, i) => <li key={i}>{inline(it.text)}{it.sub.map(draw)}</li>)}
        </Tag>
      );
    }
    case "code": return <CodeBlock key={key} lang={b.lang} code={b.code} />;
    case "table":
      return (
        <div key={key} className="doc-table-wrap" tabIndex={0}>
          <table>
            <thead><tr>{b.head.map((c, i) => <th key={i} scope="col">{inline(c)}</th>)}</tr></thead>
            <tbody>{b.rows.map((r, i) => <tr key={i}>{r.map((c, j) => <td key={j}>{inline(c)}</td>)}</tr>)}</tbody>
          </table>
        </div>
      );
    case "quote": {
      if (!b.kind) return <blockquote key={key}>{b.blocks.map(draw)}</blockquote>;
      // the label ("Note" / "Warning") leads the first paragraph; a callout that starts with a list or a code block gets it on a line of its own
      const [first, ...more] = b.blocks;
      const label = <b className="doc-callout-label">{b.kind}</b>;
      return (
        <div key={key} className={`doc-callout ${b.kind.toLowerCase()}`} role="note">
          {first?.t === "p" ? <p>{label} {inline(first.text)}</p> : <>{<p>{label}</p>}{first && draw(first, 0)}</>}
          {more.map((x, i) => draw(x, i + 1))}
        </div>
      );
    }
    case "hr": return <hr key={key} />;
    case "figure": return <DocFigure key={key} name={b.name} alt={b.alt} caption={b.caption} />;
  }
}

export function Markdown({ text }: { text: string }) {
  return <div className="doc-body">{parse(text.split("\n"), new Map()).map(draw)}</div>;
}
