import { useId, useState, type ReactNode } from "react";
import type { Cell, Row, RowGroup, RxOption, StatusChip, Tile } from "../types";
import { ThemedSvg } from "./Figure";
import { HelpButton, RichText, helpLines, useHelpPopover } from "./Help";

/**
 * Renderers for the generic output blocks. Data contracts (the spec helpers in mzdesk/spec.py and the builders in
 * mzdesk/api/blocks.py are the reference; help = string | string[] everywhere):
 *   table / grouptable cell : "text" | {"text", "url"} | {"badge", "tone"} | {"text", "help"} | {"img"}   (badge / text may also carry url + help)
 *   table row               : [cells] | {"cells": [...], "help"}   (the "?" column only with spec row_help=True)
 *   table data              : {"rows": [...], "count"?}           count -> pill next to the title
 *   grouptable data         : {"groups": [{"title", "meta"?, "help"?, "rows": [...]}], "count"?}   collapsed by default
 *   portal data             : {"tiles": [{"label", "status": yes|partial|no, "value"?, "chips"?: [str], "help"?}]}
 *   chips data              : {"chips": [{"label", "level": good|possible|poor (strong|moderate|weak, yes|no), "value"?, "help"?}]}
 *   reaction data           : {"svg"?, "options": [{"label", "steps": [{"from", "reagent"?, "to", "count"?}], "formula"?, "mass"?, "heavy"?, "help"?, "svg"?, "product_svg"?}]}
 *   cards (style "metric")  : values {label: "v" | {"value", "note"?, "help"?}}; spec item (label, unit, help?)
 *   notes (collapsed=True)  : {"items": [...]} shown in the "?" next to the title (a "Details" disclosure when the block has no title)
 */

const safeUrl = (u: string) => (/^(https?:|mailto:|\/)/i.test(u) ? u : undefined);
const asText = (v: string | number | undefined) => (typeof v === "number" ? v.toFixed(4) : v);

// ---------------------------------------------------------------- table cells and rows
/** Cell text; "\n" splits it into stacked lines. A line like "A: water + 0.1% formic" gets its label set off;
 *  other lines of a multi-line cell (e.g. "5→95% B in 20 min · 0.4 mL/min") are shown quieter. */
function TextLines({ text }: { text: string }) {
  if (!text.includes("\n")) return <RichText text={text} />;
  return (
    <span className="cell-lines">
      {text.split("\n").map((line, i) => {
        const m = /^([A-Z][\w ]{0,11}): (.*)$/.exec(line);
        return m
          ? <span key={i} className="cell-line"><span className="cell-line-key">{m[1]}</span><span><RichText text={m[2]} /></span></span>
          : <span key={i} className="cell-line rest"><RichText text={line} /></span>;
      })}
    </span>);
}

export function CellView({ cell }: { cell: Cell }) {
  if (cell === null || cell === undefined) return null;
  if (typeof cell !== "object") return typeof cell === "string" ? <TextLines text={cell} /> : <>{String(cell)}</>;
  if ("img" in cell) return <img className="cell-img" src={cell.img} alt="" />;
  const tone = String(cell.tone ?? "neutral").replace(/[^a-z0-9-]/gi, "").toLowerCase() || "neutral";
  const url = safeUrl(cell.url ?? cell.href ?? "");
  let node: ReactNode = cell.badge !== undefined ? <span className={"badge tone-" + tone}>{cell.badge}</span> : <TextLines text={cell.text ?? ""} />;
  if (url) node = <a href={url} target="_blank" rel="noopener noreferrer">{node}</a>;
  const lines = helpLines(cell.help);
  return <>{node}{lines.length > 0 && <HelpButton lines={lines} className="inline" label="Help" />}</>;
}

const cellText = (c: Cell) => (typeof c === "string" ? c : c && typeof c === "object" && "text" in c ? c.text ?? "" : "");
const isLong = (c: Cell) => cellText(c).length > 40;
const isMultiline = (c: Cell) => cellText(c).includes("\n");

/** The shared table body: header + rows (+ trailing "?" column when `rowHelp`). */
export function RowsTable({ columns, rows, rowHelp = false, empty }: { columns: string[]; rows: Row[]; rowHelp?: boolean; empty: string }) {
  return (
    <div className={"table-wrap" + (columns.length > 6 ? " wide" : "")} tabIndex={0}
      role="region" aria-label={`Results table, ${columns.length} columns; scroll horizontally for more`}>
      <table className="out-table">
        <thead><tr>{columns.map(c => <th key={c}>{c}</th>)}{rowHelp && <th className="help-col" aria-label="Help" />}</tr></thead>
        <tbody>
          {rows.length
            ? rows.map((r, i) => {
              const { cells, help } = Array.isArray(r) ? { cells: r, help: undefined } : r;
              const lines = helpLines(help);
              return (
                <tr key={i}>
                  {cells.map((c, j) => <td key={j} className={isMultiline(c) ? "long lines" : isLong(c) ? "long" : undefined}><CellView cell={c} /></td>)}
                  {rowHelp && <td className="help-col">{lines.length > 0 && <HelpButton lines={lines} label="Row help" />}</td>}
                </tr>);
            })
            : <tr><td colSpan={columns.length + (rowHelp ? 1 : 0)} className="empty">{empty}</td></tr>}
        </tbody>
      </table>
    </div>
  );
}

// ---------------------------------------------------------------- grouptable
function GroupItem({ group, columns, rowHelp }: { group: RowGroup; columns: string[]; rowHelp: boolean }) {
  const [open, setOpen] = useState(false);
  const body = useId();
  return (
    <div className={"gt-group" + (open ? " open" : "")}>
      {/* the whole header row toggles; the button inside is the keyboard target (its click bubbles to the row) */}
      <div className="gt-head" onClick={() => setOpen(o => !o)}>
        <button type="button" className="gt-toggle" aria-expanded={open} aria-controls={open ? body : undefined}>
          <span className="gt-caret" aria-hidden="true" />
          <span className="gt-title"><RichText text={group.title} /></span>
        </button>
        <HelpButton lines={helpLines(group.help)} label={`Help: ${group.title}`} />
        {group.meta && <span className="gt-meta">{group.meta}</span>}
      </div>
      {open && <div id={body} className="gt-body"><RowsTable columns={columns} rows={group.rows ?? []} rowHelp={rowHelp} empty="no rows" /></div>}
    </div>
  );
}

export function GroupTable({ columns, groups, rowHelp = false }: { columns: string[]; groups: RowGroup[]; rowHelp?: boolean }) {
  return <div className="gt">{groups.map((g, i) => <GroupItem key={i} group={g} columns={columns} rowHelp={rowHelp} />)}</div>;
}

// ---------------------------------------------------------------- portal
const STATES = { yes: { glyph: "✓", word: "Yes" }, partial: { glyph: "◐", word: "Partial" }, no: { glyph: "–", word: "No" } } as const;

function TileView({ tile }: { tile: Tile }) {
  const st = (tile.status in STATES ? tile.status : "no") as keyof typeof STATES;
  return (
    <div className={"tile st-" + st}>
      <div className="tile-top">
        <span className="tile-label">
          <span className="tile-glyph" aria-hidden="true">{STATES[st].glyph}</span>
          <span className="sr-only">{STATES[st].word}: </span>{tile.label}
        </span>
        <HelpButton lines={helpLines(tile.help)} label={`Help: ${tile.label}`} />
      </div>
      {tile.value && <div className="tile-value"><RichText text={tile.value} /></div>}
      {!!tile.chips?.length && <div className="tile-chips">{tile.chips.map((c, i) => <span className="tile-chip" key={i}>{c}</span>)}</div>}
    </div>
  );
}

export function Portal({ tiles }: { tiles: Tile[] }) {
  return <div className="tile-grid">{tiles.map((t, i) => <TileView key={i} tile={t} />)}</div>;
}

// ---------------------------------------------------------------- chips
const LEVELS: Record<string, "good" | "possible" | "poor"> = {
  good: "good", strong: "good", yes: "good", possible: "possible", moderate: "possible", partial: "possible", poor: "poor", weak: "poor", no: "poor",
};

function Pill({ chip }: { chip: StatusChip }) {
  const lines = helpLines(chip.help);
  const { triggerProps, popover } = useHelpPopover(lines);
  const lvl = LEVELS[String(chip.level).toLowerCase()] ?? "poor";
  const body = (
    <>
      <span className="sr-only">{lvl}: </span>
      <span className="pill-label">{chip.label}</span>
      {chip.value && <span className="pill-value">{chip.value}</span>}
    </>);
  // a chip with help is itself the trigger (click pins its popover)
  return lines.length
    ? <><button type="button" className={"pill has-help lvl-" + lvl} {...triggerProps}>{body}</button>{popover}</>
    : <span className={"pill lvl-" + lvl}>{body}</span>;
}

export function Chips({ chips }: { chips: StatusChip[] }) {
  return <div className="pill-row">{chips.map((c, i) => <Pill key={i} chip={c} />)}</div>;
}

// ---------------------------------------------------------------- reaction
export function Reaction({ svg, options }: { svg?: string; options: RxOption[] }) {
  const [sel, setSel] = useState(0);
  const cur = options[Math.min(sel, options.length - 1)];
  const meta: ReactNode[] = [];
  if (cur.formula) meta.push(<span key="f" className="rx-formula"><RichText text={cur.formula} /></span>);
  if (cur.mass !== undefined && cur.mass !== null) meta.push(<span key="m"><span className="rx-k">mass</span> {asText(cur.mass)}</span>);
  if (cur.heavy !== undefined && cur.heavy !== null) meta.push(<span key="h"><span className="rx-k">heavy mass</span> {asText(cur.heavy)}</span>);
  const help = helpLines(cur.help);
  const currentSvg = cur.svg ?? svg;
  const full = currentSvg && cur.product_svg ? { substrate: currentSvg, product: cur.product_svg } : null;
  const reagent = cur.steps.find((s) => s.reagent)?.reagent;
  return (
    <div className={"rx" + (!full && currentSvg ? " with-parent" : "")}>
      {!full && currentSvg && <div className="rx-parent"><ThemedSvg svg={currentSvg} className="" alt={`Parent structure with reactive groups highlighted for ${cur.label}`} /></div>}
      <div className="rx-main">
        {options.length > 1 && (
          <div className="rx-opts" role="group" aria-label="Reaction options">
            {options.map((o, i) => (
              <button key={i} type="button" className={"rx-opt" + (o === cur ? " on" : "")} aria-pressed={o === cur} onClick={() => setSel(i)}>
                <RichText text={o.label} />
              </button>))}
          </div>)}
        {full ? (
          <div className="rx-full">
            <div className="rx-full-structure">
              <div className="rx-full-label">Substrate</div>
              <ThemedSvg svg={full.substrate} className="" alt={`Substrate with reactive sites highlighted for ${cur.label}`} />
            </div>
            <div className="rx-arrow rx-full-arrow">
              {reagent && <span className="rx-reagent"><span className="sr-only">with </span><RichText text={reagent} /></span>}
              <span className="rx-line" aria-hidden="true" />
            </div>
            <div className="rx-full-structure">
              <div className="rx-full-label">Product</div>
              <ThemedSvg svg={full.product} className="" alt={`Complete product structure for ${cur.label}`} />
            </div>
          </div>
        ) : (
          <div className="rx-steps">
            {(cur.steps ?? []).map((s, i) => (
              <div className="rx-step" key={i}>
                <span className="rx-box"><RichText text={s.from} /></span>
                <span className="rx-arrow">
                  {s.reagent && <span className="rx-reagent"><span className="sr-only">with </span><RichText text={s.reagent} /></span>}
                  <span className="rx-line" aria-hidden="true" />
                </span>
                <span className="sr-only"> gives </span>
                <span className="rx-box to"><RichText text={s.to} /></span>
                {(s.count ?? 1) > 1 && <span className="rx-count" title={`${s.count} groups react`}>×{s.count}</span>}
              </div>))}
          </div>
        )}
        {(meta.length > 0 || help.length > 0) && (
          <div className="rx-meta">
            {meta.map((m, i) => <span key={i} className="rx-item">{m}{i < meta.length - 1 && <span className="rx-dot" aria-hidden="true">·</span>}</span>)}
            <HelpButton lines={help} label={`Help: ${cur.label}`} />
          </div>)}
      </div>
    </div>
  );
}
