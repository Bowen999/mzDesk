import { useState } from "react";
import type { Block, BlockData, CardValue } from "../types";
import "../outputs-v2.css";
import { Chips, GroupTable, Portal, Reaction, RowsTable } from "./Blocks";
import { FigureView } from "./Figure";
import { HelpButton, helpLines } from "./Help";

/** 一个输出块。有数据就画结果；没数据时画格式（空表只有表头，卡片是 —，图是虚线框）。
 *  status：pending = 还没运行；none = 运行了但这个块没有结果；off = 这个块还没接引擎。
 *
 *  Generic blocks (data contracts: components/Blocks.tsx header, spec helpers in mzdesk/spec.py):
 *  every block takes `help` (spec) + data.help (runtime) -> a "?" popover after its title; new types: portal, chips, reaction, grouptable;
 *  extended: table (rich cells, row help, count pill), cards (style "metric", per-item help), notes (collapsed). */
type Status = "pending" | "none" | "off";

/** Block title row: title, optional count pill, "?" (help from spec + runtime + `extra`), optional "not connected yet" tag.
 *  A block with help but no title still gets a bare "?" row. */
function Title({ block, data, status, extra }: { block: Block; data?: BlockData; status?: Status; extra?: string[] }) {
  const lines = helpLines(block.help, data?.help, extra);
  if (!block.title && !lines.length) return null;
  const count = data?.count;
  return (
    <h4 className="out-title">
      {block.title}
      {count !== undefined && count !== null && <span className="count-pill">{count}</span>}
      <HelpButton lines={lines} label={block.title ? `Help: ${block.title}` : "Help"} />
      {status === "off" && <span className="tag">not connected yet</span>}
    </h4>);
}

function cardParts(v: CardValue | undefined): [string, string | undefined, string[]] {
  if (v === undefined) return ["—", undefined, []];
  return typeof v === "string" ? [v, undefined, []] : [v.value, v.note, helpLines(v.help)];
}

function CopyButton({ text }: { text: string }) {
  const [done, setDone] = useState(false);
  return <button type="button" className="copy" onClick={() => navigator.clipboard?.writeText(text).then(() => { setDone(true); setTimeout(() => setDone(false), 1200); })}>{done ? "Copied" : "Copy"}</button>;
}

/** 同一内容的多种格式：一排小标签 + 一个代码框；记住上次选的标签。 */
function FormatsBlock({ block, data, empty }: { block: Block; data?: BlockData; empty: string }) {
  const [tab, setTab] = useState("");
  const names = Object.keys(data?.formats ?? {});
  const cur = names.includes(tab) ? tab : names[0];
  const text = cur ? data!.formats![cur] : undefined;
  const title = block.title ?? block.label;
  return (
    <div className="out">
      <h4 className="out-title">{title}<HelpButton lines={helpLines(block.help, data?.help)} label={`Help: ${title}`} />{text !== undefined && <CopyButton text={text} />}</h4>
      {names.length > 0 && (
        <div className="tabs" role="tablist">
          {names.map(n => <button key={n} type="button" role="tab" aria-selected={n === cur} className={"tab" + (n === cur ? " on" : "")} onClick={() => setTab(n)}>{n}</button>)}
        </div>)}
      <pre className="code">{text ?? `— ${block.label} (${empty}) —`}</pre>
    </div>);
}

export function OutputBlock({ block, data, status }: { block: Block; data?: BlockData; status: Status }) {
  const empty = status === "none" ? "no result" : status === "off" ? "not connected yet" : "empty";
  const shown = data ? undefined : status;
  switch (block.type) {
    case "heading": return <h4 className="out-title">{block.title}<HelpButton lines={helpLines(block.help, data?.help)} label={`Help: ${block.title}`} /></h4>;
    case "table": return (
      <div className="out">
        <Title block={block} data={data} status={shown} />
        <RowsTable columns={block.columns!} rows={data?.rows ?? []} rowHelp={block.row_help} empty={empty} />
        {block.caption && <small>{block.caption}</small>}
      </div>);
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
          return (
            <div className={"card" + (data ? "" : " blank")} key={it.label}>
              <div className="card-label">{it.label}<HelpButton lines={lines} label={`Help: ${it.label}`} /></div>
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
        {data?.options?.length ? <Reaction svg={data.svg} options={data.options} /> : <div className="placeholder-box">{empty}</div>}
      </div>);
    case "figure":
      if (data) return (
        <div className="out">
          <div className="fig-label">{block.label}<HelpButton lines={helpLines(block.help, data.help)} label={`Help: ${block.label}`} /></div>
          <FigureView data={data} />
        </div>);
      return <div className="out placeholder-box">[{empty}] {block.label}</div>;
    case "code": return (
      <div className="out">
        <h4 className="out-title">{block.title ?? block.label}<HelpButton lines={helpLines(block.help, data?.help)} label={`Help: ${block.title ?? block.label}`} />{data?.text && <CopyButton text={data.text} />}</h4>
        <pre className="code">{data?.text ?? `— ${block.label} (${empty}) —`}</pre>
      </div>);
    case "formats": return <FormatsBlock block={block} data={data} empty={empty} />;
    case "notes": {
      const items = data?.items ?? [];
      if (block.collapsed) {
        // collapsed: the items live in the "?" next to the title; without a title, a small "Details" disclosure.
        // After a run with no items there is nothing to show: render nothing, not even the title. (Before the run the title is the format preview.)
        const ran = !!data || status === "none";
        if (ran && !items.length) return null;
        if (block.title) return <div className="out"><Title block={block} data={data} status={shown} extra={items} /></div>;
        const lines = helpLines(block.help, data?.help);
        if (!items.length && !lines.length) return null;
        return (
          <div className="out">
            <details className="out-details">
              <summary>Details</summary>
              <ul className="notes">{[...lines, ...items].map((t, i) => <li key={i}>{t}</li>)}</ul>
            </details>
          </div>);
      }
      return (
        <div className="out">
          <Title block={block} data={data} status={shown} />
          {items.length ? <ul className="notes">{items.map((t, i) => <li key={i}>{t}</li>)}</ul> : <div className="placeholder-box">{empty}</div>}
        </div>);
    }
    case "button": return <button className="btn secondary" disabled>{block.label}</button>;
    default: return <div className="out">unknown block {block.type}</div>;
  }
}

/** 一组输出块。outputs 为 null = 还没运行。live = 模块已接引擎（没 id 的块标成 not connected yet）。 */
export function OutputList({ blocks, exclude = [], outputs = null, live = false }:
  { blocks: Block[]; exclude?: string[]; outputs?: Record<string, BlockData> | null; live?: boolean }) {
  const shown = blocks.filter(b => !b.key || !exclude.includes(b.key));
  return (
    <div className="outputs">
      {shown.map((b, i) => {
        const data = b.id && outputs ? outputs[b.id] : undefined;
        const status: Status = live && !b.id ? "off" : outputs ? "none" : "pending";
        return <OutputBlock key={b.id ?? i} block={b} data={data} status={status} />;
      })}
    </div>
  );
}
