import type { Block } from "../types";

/** 输出块只画格式：空表只有表头，卡片是 —，图是虚线框。 */
export function OutputBlock({ block }: { block: Block }) {
  switch (block.type) {
    case "heading": return <h4 className="out-title">{block.title}</h4>;
    case "table": return (
      <div className="out">
        {block.title && <h4 className="out-title">{block.title}</h4>}
        <table className="placeholder-table">
          <thead><tr>{block.columns!.map(c => <th key={c}>{c}</th>)}</tr></thead>
          <tbody><tr><td colSpan={block.columns!.length} className="empty">empty</td></tr></tbody>
        </table>
        {block.caption && <small>{block.caption}</small>}
      </div>);
    case "cards": return (
      <div className="out">
        {block.title && <h4 className="out-title">{block.title}</h4>}
        <div className="cards">{block.items!.map(it => <div className="card" key={it.label}><div className="card-label">{it.label}</div><div className="card-value">—</div>{it.unit && <div className="card-unit">{it.unit}</div>}</div>)}</div>
      </div>);
    case "figure": return <div className="out placeholder-box">[figure placeholder] {block.label}</div>;
    case "code": return (
      <div className="out">
        {block.title && <h4 className="out-title">{block.title}</h4>}
        <pre className="code">— {block.label} —</pre>
      </div>);
    case "button": return <button className="btn secondary" disabled>{block.label}</button>;
    default: return <div className="out">unknown block {block.type}</div>;
  }
}

export function OutputList({ blocks, exclude = [] }: { blocks: Block[]; exclude?: string[] }) {
  const shown = blocks.filter(b => !b.key || !exclude.includes(b.key));
  return <div className="outputs">{shown.map((b, i) => <OutputBlock key={i} block={b} />)}</div>;
}
