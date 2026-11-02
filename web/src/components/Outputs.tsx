import { useEffect, useState } from "react";
import type React from "react";
import { specContext, useT } from "../i18n";
import { useReportMode } from "../report/mode";
import type { Block, BlockData, CardValue, Cell, Help, Row } from "../types";
import "../outputs-v2.css";
import { ChipList, Chips, GroupTable, Portal, Reaction, RowsTable, isChipList } from "./Blocks";
import { CompoundCompare } from "./CompoundCompare";
import { FigureView } from "./Figure";
import { FragmentProvider } from "./FragmentDialog";
import { GroupMap } from "./GroupMap";
import { HelpButton, RichText, emphasize, helpLines } from "./Help";
import { MiddleEllipsis } from "./MiddleEllipsis";

/** 一个输出块。有数据就画结果；没数据时画格式（空表只有表头，卡片是 —，图是虚线框）。
 *  status：pending = 还没运行；none = 运行了但这个块没有结果；off = 这个块还没接引擎。
 *
 *  Generic blocks (data contracts: components/Blocks.tsx header, spec helpers in mzdesk/spec.py):
 *  every block takes `help` (spec) + data.help (runtime) -> a "?" popover after its title; new types: portal, chips, reaction, grouptable;
 *  extended: table (rich cells, row help, count pill, data.columns overriding the spec columns), cards (style "metric", per-item help), notes (collapsed). */
type Status = "pending" | "none" | "off";

/** Block title row: title, optional count pill, "?" (help from spec + runtime + `extra`), optional "not connected yet" tag.
 *  A block with help but no title still gets a bare "?" row. */
function Title({ block, data, status, extra }: { block: Block; data?: BlockData; status?: Status; extra?: string[] }) {
  const t = useT();
  const lines = helpLines(block.help, data?.help, extra);
  if (!block.title && !lines.length) return null;
  const count = data?.count;
  return (
    <h4 className="out-title">
      {block.title}
      {count !== undefined && count !== null && <span className="count-pill">{count}</span>}
      <HelpButton lines={lines} label={block.title ? t("Help: {name}", { name: block.title }) : t("Help")} />
      {status === "off" && <span className="tag">{t("not connected yet")}</span>}
    </h4>);
}

function cardParts(v: CardValue | undefined): [string, string | undefined, string[]] {
  if (v === undefined) return ["—", undefined, []];
  return typeof v === "string" ? [v, undefined, []] : [v.value, v.note, helpLines(v.help)];
}

function CopyButton({ text }: { text: string }) {
  const t = useT();
  const [done, setDone] = useState(false);
  return <button type="button" className="copy" onClick={() => navigator.clipboard?.writeText(text).then(() => { setDone(true); setTimeout(() => setDone(false), 1200); })}>{done ? t("Copied") : t("Copy")}</button>;
}

/** 同一内容的多种格式：一排小标签 + 一个代码框；记住上次选的标签。
 *  Several compounds (Structure & formats) add a compound tab row; the format tab stays selected when switching. */
function FormatsBlock({ block, data, empty }: { block: Block; data?: BlockData; empty: string }) {
  const t = useT();
  const report = useReportMode();
  const [tab, setTab] = useState("");
  const [setIndex, setSetIndex] = useState(0);
  useEffect(() => setSetIndex(0), [data]);
  const sets = data?.sets ?? [];
  const active = sets.length ? (sets[Math.min(setIndex, sets.length - 1)] ?? sets[0]) : null;
  const formats = active?.formats ?? data?.formats ?? {};
  const names = Object.keys(formats);
  const cur = names.includes(tab) ? tab : names[0];
  const text = cur ? formats[cur] : undefined;
  const title = block.title ?? block.label;
  // the report has no tabs: every format of every compound, one after the other
  if (report) return (
    <div className="out">
      <h4 className="out-title">{title}</h4>
      {(sets.length ? sets : [{ label: "", formats }]).map((set, i) => (
        <div key={i}>
          {sets.length > 1 && <div className="fig-label">{set.label}</div>}
          {Object.entries(set.formats ?? {}).map(([name, value]) => <div key={name}><div className="fig-label">{name}</div><pre className="code">{value}</pre></div>)}
        </div>))}
    </div>);
  return (
    <div className="out">
      <h4 className="out-title">{title}<HelpButton lines={helpLines(block.help, data?.help)} label={t("Help: {name}", { name: title ?? "" })} />{text !== undefined && <CopyButton text={text} />}</h4>
      {sets.length > 1 && (
        <div className="tabs compound-tabs" role="tablist" aria-label={t("Compound")}>
          {sets.map((s, i) => <button key={`${s.label}-${i}`} type="button" role="tab" aria-selected={i === Math.min(setIndex, sets.length - 1)}
            className={"tab" + (i === Math.min(setIndex, sets.length - 1) ? " on" : "")} onClick={() => setSetIndex(i)}>
            <MiddleEllipsis text={s.label} maxWidth={216} /></button>)}
        </div>)}
      {names.length > 0 && (
        <div className="tabs" role="tablist">
          {names.map(n => <button key={n} type="button" role="tab" aria-selected={n === cur} className={"tab" + (n === cur ? " on" : "")} onClick={() => setTab(n)}>{n}</button>)}
        </div>)}
      <pre className="code">{text ?? `— ${block.label} (${empty}) —`}</pre>
    </div>);
}

const FOLD_AT = 14;                 // a table longer than this shows its first FOLD_SHOW rows and a button for the rest
const FOLD_SHOW = 8;

/** A long table without its own pages: the first rows, then "Show all N rows". Short tables are never folded. */
function FoldedTable({ rows, ...rest }: { columns: string[]; rows: Row[]; rowHelp?: boolean; columnHelp?: Record<string, Help>; empty: string; compact?: boolean }) {
  const t = useT();
  const [all, setAll] = useState(false);
  useEffect(() => setAll(false), [rows]);
  return (
    <>
      <RowsTable {...rest} rows={all ? rows : rows.slice(0, FOLD_SHOW)} />
      <button type="button" className="fold-more" aria-expanded={all} onClick={() => setAll(v => !v)}>{all ? t("Show fewer rows") : t("Show all {n} rows", { n: rows.length })}</button>
    </>
  );
}

/** Pagination and search operate only on the rows returned for the preview, never on the full download. */
function ResultTable({ block, data, empty }: { block: Block; data?: BlockData; empty: string }) {
  const t = useT();
  const report = useReportMode();
  const rows = data?.rows ?? [];
  const columns = data?.columns ?? block.columns ?? [];
  const size = block.page_size && block.page_size > 0 ? block.page_size : 0;
  const [page, setPage] = useState(0);
  const [filter, setFilter] = useState("");
  useEffect(() => { setPage(0); setFilter(""); }, [data]);
  if (!size) {
    const props = { columns, rows, rowHelp: block.row_help, columnHelp: block.column_help, empty, compact: block.compact };
    return !report && rows.length > FOLD_AT ? <FoldedTable {...props} /> : <RowsTable {...props} />;
  }
  if (report) {         // no pages, no filter: every row (RowsTable caps a very long table and says so)
    const count = typeof data?.count === "number" ? data.count : rows.length;
    return <div className="result-table">
      {count > rows.length && <p className="result-table-note">{t("Preview: first {shown} of {total} results.", { shown: rows.length, total: count.toLocaleString("en-US") })}</p>}
      <RowsTable columns={columns} rows={rows} rowHelp={block.row_help} columnHelp={block.column_help} empty={empty} compact={block.compact} />
    </div>;
  }
  const cellText = (cell: Cell) => {
    if (cell === null || cell === undefined) return "";
    if (typeof cell !== "object") return String(cell);
    return "img" in cell || "fragment" in cell ? "" : String(cell.text ?? cell.badge ?? "");
  };
  const searchable = (row: Row) => (Array.isArray(row) ? row : row.cells).map(cellText).join(" ").toLocaleLowerCase();
  const needle = filter.trim().toLocaleLowerCase();
  const filtered = needle ? rows.filter(row => searchable(row).includes(needle)) : rows;
  const pages = Math.max(1, Math.ceil(filtered.length / size));
  const current = Math.min(page, pages - 1);
  const start = current * size;
  const total = typeof data?.count === "number" ? data.count : rows.length;
  const limited = total > rows.length;
  const controls = rows.length > size;
  const range = { from: start + 1, to: Math.min(start + size, filtered.length), total: filtered.length };
  return <div className="result-table">
    {controls && <div className="result-table-tools">
      <label className="result-table-filter"><span>{t("Filter preview")}</span>
        <input type="search" value={filter} placeholder={t("Search these rows")} aria-label={block.title ? t("Filter {name} preview", { name: block.title }) : t("Filter table preview")}
          onChange={e => { setFilter(e.target.value); setPage(0); }} />
      </label>
      <span className="result-table-count" aria-live="polite">
        {filtered.length ? (needle ? t("{from}–{to} of {total} matching preview rows", range) : t("{from}–{to} of {total} preview rows", range))
          : needle ? t("No matching preview rows") : t("No preview rows")}
      </span>
    </div>}
    {limited && <p className="result-table-note">{t("Preview: first {shown} of {total} results. Filtering only searches this preview.", { shown: rows.length, total: total.toLocaleString() })}</p>}
    <RowsTable key={`${current}:${needle}`} columns={columns} rows={filtered.slice(start, start + size)} rowHelp={block.row_help} columnHelp={block.column_help}
      empty={needle ? t("No preview rows match this filter.") : empty} compact={block.compact} />
    {controls && <nav className="result-table-pages" aria-label={block.title ? t("{name} pages", { name: block.title }) : t("Table pages")}>
      <span aria-live="polite">{t("Page {current} of {pages}", { current: current + 1, pages })}</span>
      <div><button type="button" className="page-button" disabled={current === 0} onClick={() => setPage(current - 1)}>{t("Previous")}</button>
        <button type="button" className="page-button" disabled={current >= pages - 1} onClick={() => setPage(current + 1)}>{t("Next")}</button></div>
    </nav>}
  </div>;
}

export function OutputBlock({ block, data, status }: { block: Block; data?: BlockData; status: Status }) {
  const t = useT();
  const report = useReportMode();
  const empty = status === "none" ? t("no result") : status === "off" ? t("not connected yet") : t("empty");
  const shown = data ? undefined : status;
  switch (block.type) {
    case "heading": return <h4 className="out-title">{block.title}<HelpButton lines={helpLines(block.help, data?.help)} label={t("Help: {name}", { name: block.title ?? "" })} /></h4>;
    case "table": {
      const list = !!data?.rows && !block.page_size && isChipList(data.columns ?? block.columns ?? [], data.rows);
      return (
        <div className="out">
          <Title block={block} data={data} status={shown} />
          {list ? <ChipList rows={data!.rows!} column={(data!.columns ?? block.columns ?? [])[0]} />
            : <FragmentProvider query={data?.fragments}><ResultTable block={block} data={data} empty={empty} /></FragmentProvider>}
          {block.caption && <small>{block.caption}</small>}
        </div>);
    }
    case "grouptable": return (
      <div className="out">
        <Title block={block} data={data} status={shown} />
        {data?.groups?.length ? <GroupTable columns={block.columns!} groups={data.groups} rowHelp={block.row_help} /> : <div className="placeholder-box">{empty}</div>}
        {block.caption && <small>{block.caption}</small>}
      </div>);
    case "cards": return (
      <div className="out">
        <Title block={block} data={data} status={shown} />
        <div className={"cards" + (block.style === "metric" ? " metric" : "")}>{block.items!.map(it => {
          const [value, note, vhelp] = cardParts(data?.values?.[it.label]);
          const lines = helpLines(it.help, vhelp);
          const label = t(it.label, specContext(block.id, it.label));          // a cards label stays English in the spec (it is the key into the run's values): the dictionary translates it here
          return (
            <div className={"card" + (data ? "" : " blank")} key={it.label}>
              <div className="card-label">{label}<HelpButton lines={lines} label={t("Help: {name}", { name: label })} /></div>
              <div className="card-value">{value}</div>
              {(note ?? it.unit) && <div className="card-unit">{note ?? it.unit}</div>}
            </div>);
        })}</div>
      </div>);
    case "portal": return (
      <div className="out">
        <Title block={block} data={data} status={shown} />
        {data?.tiles?.length ? <Portal tiles={data.tiles} /> : <div className="placeholder-box">{empty}</div>}
      </div>);
    case "chips": return (
      <div className="out">
        <Title block={block} data={data} status={shown} />
        {data?.chips?.length ? <Chips chips={data.chips} /> : <div className="placeholder-box">{empty}</div>}
      </div>);
    case "reaction":
      // after a run with no options there is nothing to draw: render nothing, not even the title (before the run the title is the format preview)
      if ((!!data || status === "none") && !data?.options?.length) return null;
      return (
      <div className="out">
        <Title block={block} data={data} status={shown} />
        {data?.options?.length
          // the report has no option buttons: every option, each with its label
          ? report ? data.options.map((o, i) => <div key={i}><div className="fig-label"><RichText text={o.label} /></div><Reaction svg={data.svg} options={[o]} /></div>)
            : <Reaction svg={data.svg} options={data.options} />
          : <div className="placeholder-box">{empty}</div>}
      </div>);
    case "groupmap": return (
      <div className="out">
        <Title block={block} data={data} status={shown} />
        {data?.found?.length && data.svg
          // the report cannot pick a group: one panel per group
          ? report ? data.found.map((g, i) => <GroupMap key={i} svg={data.svg!} found={[g]} />) : <GroupMap svg={data.svg} found={data.found} />
          : <div className="placeholder-box">{empty}</div>}
      </div>);
    case "figure":
      if (data) return (
        <div className="out">
          <div className="fig-label">{block.label}<HelpButton lines={helpLines(block.help, data.help)} label={t("Help: {name}", { name: block.label ?? "" })} /></div>
          <FigureView data={data} />
        </div>);
      return <div className="out placeholder-box">[{empty}] {block.label}</div>;
    case "compare": return (
      <div className="out">
        <Title block={block} data={data} status={shown} />
        {data?.compounds?.length ? <CompoundCompare items={data.compounds} /> : <div className="placeholder-box">{empty}</div>}
      </div>);
    case "code": return (
      <div className="out">
        <h4 className="out-title">{block.title ?? block.label}<HelpButton lines={helpLines(block.help, data?.help)} label={t("Help: {name}", { name: block.title ?? block.label ?? "" })} />{data?.text && <CopyButton text={data.text} />}</h4>
        <pre className="code">{data?.text ?? `— ${block.label} (${empty}) —`}</pre>
      </div>);
    case "formats": return <FormatsBlock block={block} data={data} empty={empty} />;
    case "downloads": return (
      <div className="out">
        <Title block={block} data={data} status={shown} />
        <div className="download-list">{data?.files?.map(file => <button key={file.name} type="button" className="btn secondary" onClick={() => {
          const url = URL.createObjectURL(new Blob([file.content], { type: file.mime ?? "text/plain;charset=utf-8" }));
          const a = document.createElement("a"); a.href = url; a.download = file.name; a.click();
          setTimeout(() => URL.revokeObjectURL(url), 1000);
        }}>{t("Download {name}", { name: file.name })}</button>)}</div>
      </div>);
    case "notes": {
      const items = data?.items ?? [];
      if (block.collapsed) {
        // collapsed: the items live in the "?" next to the title; without a title, a small "Details" disclosure.
        // After a run with no items there is nothing to show: render nothing, not even the title. (Before the run the title is the format preview.)
        const ran = !!data || status === "none";
        if (ran && !items.length) return null;
        if (block.title && !report) return <div className="out"><Title block={block} data={data} status={shown} extra={items} /></div>;
        const lines = helpLines(block.help, data?.help);
        if (!items.length && !lines.length) return null;
        return (
          <div className="out">
            <details className="out-details">
              <summary>{block.title || t("Details")}</summary>
              <ul className="notes">{[...lines, ...items].map((t, i) => <li key={i}><RichText text={emphasize(t)} /></li>)}</ul>
            </details>
          </div>);
      }
      return (
        <div className="out">
          <Title block={block} data={data} status={shown} />
          {items.length ? <ul className="notes">{items.map((t, i) => <li key={i}><RichText text={emphasize(t)} /></li>)}</ul> : <div className="placeholder-box">{empty}</div>}
        </div>);
    }
    case "button": return <button className="btn secondary" disabled>{block.label}</button>;
    default: return <div className="out">unknown block {block.type}</div>;
  }
}

/** Where a block sits in the results grid and how loud it is.
 *  span: "full" = whole width; "half" = shares a row with another half block on wide screens (spec `span` wins).
 *  tier: "hero" = metric cards, no frame (they sit right under the verdict); "card" = framed (figures, wide tables, tiles);
 *        "quiet" = no frame, just a title and a hairline (chips, notes, short tables, lists of short values).
 *  `data` (after a run) lets a list adapt to its length: a short list of chips shares a row, a long one takes the whole width. */
export function blockLayout(b: Block, d?: BlockData): { span: "half" | "full"; tier: "hero" | "card" | "quiet" } {
  const columns = d?.columns ?? b.columns ?? [];
  const list = b.type === "table" && !b.page_size && !!d?.rows && isChipList(columns, d.rows);
  const narrowTable = (b.type === "table") && columns.length <= 3;
  const tier = b.type === "cards" ? "hero" : narrowTable || b.type === "chips" || b.type === "notes" || b.type === "downloads" ? "quiet" : "card";
  const span = b.span ?? (list ? (d!.rows!.length > 24 ? "full" : "half") : narrowTable || b.type === "chips" || b.type === "notes" ? "half" : "full");
  return { span, tier };
}

/** No rows / chips / tiles / items: after a run such a block is hidden instead of showing an empty frame. */
function isEmpty(b: Block, d: BlockData): boolean {
  switch (b.type) {
    case "table": return !(d.rows?.length);
    case "grouptable": return !(d.groups?.length);
    case "chips": return !(d.chips?.length);
    case "portal": return !(d.tiles?.length);
    case "reaction": return !(d.options?.length);
    case "notes": return !(d.items?.length) && !helpLines(d.help).length;
    case "formats": return !Object.keys(d.formats ?? {}).length && !d.sets?.some(s => Object.keys(s.formats ?? {}).length);
    case "compare": return !(d.compounds?.length);
    case "groupmap": return !(d.found?.length);
    case "downloads": return !d.files?.length;
    case "code": return !d.text;
    default: return false;
  }
}

/** Grey stand-in for a block before the run (or while it runs): its title and the rough shape of its content. */
function Skeleton({ block }: { block: Block }) {
  const title = block.title ?? block.label;
  const shape = block.type === "figure" ? <div className="sk sk-figure" />
    : block.type === "cards" ? <div className="sk-cards">{(block.items ?? [1, 2, 3]).slice(0, 5).map((_, i) => <div key={i} className="sk sk-card" />)}</div>
    : block.type === "chips" ? <div className="sk-row">{[56, 72, 64, 88].map((w, i) => <div key={i} className="sk sk-pill" style={{ width: w }} />)}</div>
    : <div className="sk-lines"><div className="sk sk-line" /><div className="sk sk-line" /><div className="sk sk-line short" /></div>;
  return <div className="out skeleton" aria-hidden="true">{title && <div className="out-title">{title}</div>}{shape}</div>;
}

/** 一组输出块。outputs 为 null = 还没运行。live = 模块已接引擎（没 id 的块标成 not connected yet）。
 *  Live module: before the run every block is a skeleton (busy = shimmering); after it, blocks without results are hidden.
 *  Placeholder module (live = false): the blocks show their format, as before. */
export function OutputList({ blocks, exclude = [], hide = [], outputs = null, live = false, busy = false }:
  { blocks: Block[]; exclude?: string[]; hide?: string[]; outputs?: Record<string, BlockData> | null; live?: boolean; busy?: boolean }) {
  const report = useReportMode();
  // the report is a file: "Download" and placeholder buttons would do nothing (the files are in report.json)
  // hide: block ids that say nothing the page does not already say (see coveredCards)
  const shown = blocks.filter(b => (!b.key || !exclude.includes(b.key)) && !(b.id && hide.includes(b.id)) && !(report && (b.type === "downloads" || b.type === "button")));
  const skeleton = live && (busy || !outputs);
  let n = 0;
  return (
    <div className={"outputs" + (skeleton ? " is-skeleton" : "") + (busy ? " is-busy" : "") + (outputs && !busy ? " has-data" : "")}>
      {shown.map((b, i) => {
        const data = b.id && outputs ? outputs[b.id] : undefined;
        const { span, tier } = blockLayout(b, data);
        const cls = `block span-${span} tier-${tier}`;
        if (skeleton) {
          if (b.type === "heading" || b.type === "button" || (b.type === "notes" && b.collapsed)) return null;
          return <div key={b.id ?? i} className={cls}><Skeleton block={b} /></div>;
        }
        const status: Status = live && !b.id ? "off" : outputs ? "none" : "pending";
        if (live && outputs && (status === "none" && !data || data && isEmpty(b, data))) return null;
        const el = <OutputBlock block={b} data={data} status={status} />;
        return <div key={b.id ?? i} className={cls} style={{ "--i": n++ } as React.CSSProperties}>{el}</div>;
      })}
    </div>
  );
}

/** Ids of the cards blocks whose every value is already written in the verdict: the page says it once (the verdict, which is the one that is stressed). */
export function coveredCards(blocks: Block[], outputs: Record<string, BlockData> | null | undefined, verdict: string | null | undefined): string[] {
  if (!verdict || !outputs) return [];
  const plain = verdict.replace(/<[^>]+>/g, "");
  return blocks.filter(b => {
    if (b.type !== "cards" || !b.id) return false;
    const values = Object.values(outputs[b.id]?.values ?? {}).map(v => (typeof v === "string" ? v : v.value)).filter(v => v && v !== "—");
    return values.length > 0 && values.every(v => plain.includes(v));
  }).map(b => b.id as string);
}
