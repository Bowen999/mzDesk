import { useState } from "react";
import type { Block, BlockData, CardValue, Cell } from "../types";
import { FigureView } from "./Figure";

/** 一个输出块。有数据就画结果；没数据时画格式（空表只有表头，卡片是 —，图是虚线框）。
 *  status：pending = 还没运行；none = 运行了但这个块没有结果；off = 这个块还没接引擎。 */
type Status = "pending" | "none" | "off";

function Title({ block, status }: { block: Block; status?: Status }) {
  if (!block.title) return null;
  return <h4 className="out-title">{block.title}{status === "off" && <span className="tag">not connected yet</span>}</h4>;
}

function CellView({ cell }: { cell: Cell }) {
  if (cell && typeof cell === "object" && "href" in cell) return <a href={cell.href} target="_blank" rel="noopener">{cell.text}</a>;
  if (cell && typeof cell === "object") return <img className="cell-img" src={cell.img} alt="" />;
  return <>{cell}</>;
}

function cardText(v: CardValue | undefined): [string, string | undefined] {
  if (v === undefined) return ["—", undefined];
  return typeof v === "string" ? [v, undefined] : [v.value, v.note];
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
  return (
    <div className="out">
      <h4 className="out-title">{block.title ?? block.label}{text !== undefined && <CopyButton text={text} />}</h4>
      {names.length > 0 && (
        <div className="tabs" role="tablist">
          {names.map(n => <button key={n} type="button" role="tab" aria-selected={n === cur} className={"tab" + (n === cur ? " on" : "")} onClick={() => setTab(n)}>{n}</button>)}
        </div>)}
      <pre className="code">{text ?? `— ${block.label} (${empty}) —`}</pre>
    </div>);
}

export function OutputBlock({ block, data, status }: { block: Block; data?: BlockData; status: Status }) {
  const empty = status === "none" ? "no result" : status === "off" ? "not connected yet" : "empty";
  switch (block.type) {
    case "heading": return <h4 className="out-title">{block.title}</h4>;
    case "table": return (
      <div className="out">
        <Title block={block} status={data ? undefined : status} />
        <div className="table-wrap">
          <table className="out-table">
            <thead><tr>{block.columns!.map(c => <th key={c}>{c}</th>)}</tr></thead>
            <tbody>
              {data?.rows?.length
                ? data.rows.map((r, i) => <tr key={i}>{r.map((c, j) => <td key={j} className={typeof c === "string" && c.length > 40 ? "long" : undefined}><CellView cell={c} /></td>)}</tr>)
                : <tr><td colSpan={block.columns!.length} className="empty">{empty}</td></tr>}
            </tbody>
          </table>
        </div>
        {block.caption && <small>{block.caption}</small>}
      </div>);
    case "cards": return (
      <div className="out">
        <Title block={block} status={data ? undefined : status} />
        <div className="cards">{block.items!.map(it => {
          const [value, note] = cardText(data?.values?.[it.label]);
          return (
            <div className={"card" + (data ? "" : " blank")} key={it.label}>
              <div className="card-label">{it.label}</div>
              <div className="card-value">{value}</div>
              {(note ?? it.unit) && <div className="card-unit">{note ?? it.unit}</div>}
            </div>);
        })}</div>
      </div>);
    case "figure":
      if (data) return <div className="out"><div className="fig-label">{block.label}</div><FigureView data={data} /></div>;
      return <div className="out placeholder-box">[{empty}] {block.label}</div>;
    case "code": return (
      <div className="out">
        <h4 className="out-title">{block.title ?? block.label}{data?.text && <CopyButton text={data.text} />}</h4>
        <pre className="code">{data?.text ?? `— ${block.label} (${empty}) —`}</pre>
      </div>);
    case "formats": return <FormatsBlock block={block} data={data} empty={empty} />;
    case "notes": return (
      <div className="out">
        <Title block={block} status={data ? undefined : status} />
        {data?.items?.length ? <ul className="notes">{data.items.map((t, i) => <li key={i}>{t}</li>)}</ul> : <div className="placeholder-box">{empty}</div>}
      </div>);
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
