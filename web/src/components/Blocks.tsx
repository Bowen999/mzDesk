import { useId, useLayoutEffect, useMemo, useRef, useState, type ReactNode } from "react";
import { useT } from "../i18n";
import { useReportMode } from "../report/mode";
import type { Cell, Help, Row, RowGroup, RxOption, StatusChip, Tile } from "../types";
import { FragmentButton } from "./FragmentDialog";
import { ThemedSvg } from "./Figure";
import { HelpButton, RichText, helpLines, useHelpPopover } from "./Help";
import { StructureHover, structureForCell, useStructures } from "./StructureHover";

/**
 * Renderers for the generic output blocks. Data contracts (the spec helpers in mzdesk/spec.py and the builders in
 * mzdesk/api/blocks.py are the reference; help = string | string[] everywhere):
 *   table / grouptable cell : "text" | {"text", "url"} | {"badge", "tone"} | {"text", "help"} | {"text", "hover_svg"} | {"img"} | {"fragment": {smiles, adduct, name?}}
 *                             (badge / text may also carry url + help; a fragment cell is a "View" button, its spectrum is the table's data.fragments)
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

/** Text (usually a SMILES) whose 2D structure is drawn on hover / keyboard focus, in the same popover as the "?" help; click pins it. */
function Depicted({ svg, children }: { svg: string; children: ReactNode }) {
  const t = useT();
  const { triggerProps, popover } = useHelpPopover([], <ThemedSvg svg={svg} className="pop-struct" alt={t("2D structure")} />);
  return <><button type="button" className="cell-depict" {...triggerProps}>{children}</button>{popover}</>;
}

export function CellView({ cell, column, rowSmiles }: { cell: Cell; column?: string; rowSmiles?: string }) {
  const t = useT();
  const structures = useStructures();      // ESI MS/MS Annotation: a name or SMILES shows its structure on hover (components/StructureHover.tsx)
  if (cell === null || cell === undefined) return null;
  if (typeof cell !== "object") {
    if (typeof cell !== "string") return <>{String(cell)}</>;
    const smiles = structureForCell(structures, cell, column, rowSmiles);
    return smiles ? <StructureHover smiles={smiles}><TextLines text={cell} /></StructureHover> : <TextLines text={cell} />;
  }
  if ("img" in cell) return <img className="cell-img" src={cell.img} alt="" />;
  if ("fragment" in cell) return <><FragmentButton target={cell.fragment} />{helpLines(cell.help).length > 0 && <HelpButton lines={helpLines(cell.help)} className="inline" label={t("Help")} />}</>;
  const tone = String(cell.tone ?? "neutral").replace(/[^a-z0-9-]/gi, "").toLowerCase() || "neutral";
  const url = safeUrl(cell.url ?? cell.href ?? "");
  let node: ReactNode = cell.badge !== undefined ? <span className={"badge tone-" + tone}>{cell.badge}</span> : <TextLines text={cell.text ?? ""} />;
  if (url) node = <a href={url} target="_blank" rel="noopener noreferrer">{node}</a>;
  if (cell.hover_svg) node = <Depicted svg={cell.hover_svg}>{node}</Depicted>;
  else if (!url && cell.badge === undefined && typeof cell.text === "string") {
    const smiles = structureForCell(structures, cell.text, column, rowSmiles);
    if (smiles) node = <StructureHover smiles={smiles}>{node}</StructureHover>;
  }
  const lines = helpLines(cell.help);
  return <>{node}{lines.length > 0 && <HelpButton lines={lines} className="inline" label={t("Help")} />}</>;
}

const cellText = (c: Cell) => (typeof c === "string" ? c : c && typeof c === "object" && "text" in c ? c.text ?? "" : "");
const isLong = (c: Cell) => cellText(c).length > 40;
const isMultiline = (c: Cell) => cellText(c).includes("\n");

/** Keep long identifiers readable without discarding their full value. */
function CompactText({ text, column }: { text: string; column: string }) {
  const t = useT();
  const [open, setOpen] = useState(false);
  if (useReportMode()) return <TextLines text={text} />;       // the report has no button to open it: the whole text
  const short = text.startsWith("mzspec:") ? `${text.slice(0, 46)}…${text.slice(-18)}` : `${text.slice(0, 80)}…`;
  return <span className={"compact-cell" + (open ? " is-open" : "")}>
    <span className="compact-cell-text">{open ? <TextLines text={text} /> : short}</span>
    <button type="button" className="cell-expand" aria-expanded={open} aria-label={open ? t("Collapse {column}", { column }) : t("Show full {column}", { column })}
      onClick={() => setOpen(v => !v)}>{open ? t("Less") : t("More")}</button>
  </span>;
}

/** What a table shows in the Session report: no column of "View" buttons (the dialog needs the app), and at most `max` rows. */
function reportView(columns: string[], rows: Row[], max: number) {
  const cellsOf = (r: Row) => (Array.isArray(r) ? r : r.cells);
  const isView = (c: Cell) => !!c && typeof c === "object" && "fragment" in c;
  const drop = new Set(columns.map((_, j) => j).filter(j => rows.some(r => isView(cellsOf(r)[j]))
    && rows.every(r => { const c = cellsOf(r)[j]; return isView(c) || c === "—" || c === null || c === undefined; })));
  const keep = <T,>(cells: T[]) => (drop.size ? cells.filter((_, j) => !drop.has(j)) : cells);
  return {
    columns: keep(columns), total: rows.length,
    rows: rows.slice(0, max).map(r => (Array.isArray(r) ? keep(r) : { ...r, cells: keep(r.cells) })),
  };
}

/** The shared table body: header + rows (+ trailing "?" column when `rowHelp`; `columnHelp` puts a small "?" after a column's name).
 *  A wide table freezes its first column (theme.css); a "Fragments" column right after it is frozen too (fragment-view.css), so the View button
 *  stays in reach next to the compound while the rest scrolls. */
export function RowsTable({ columns: allColumns, rows: allRows, rowHelp: withHelp = false, empty, compact = false, columnHelp }:
  { columns: string[]; rows: Row[]; rowHelp?: boolean; empty: string; compact?: boolean; columnHelp?: Record<string, Help> }) {
  const t = useT();
  const report = useReportMode();
  const rowHelp = withHelp && !report;       // no "?" column in the report: there are no popovers
  const view = useMemo(() => (report ? reportView(allColumns, allRows, report.maxRows) : null), [report, allColumns, allRows]);
  const columns = view?.columns ?? allColumns, rows = view?.rows ?? allRows;
  const wrap = useRef<HTMLDivElement>(null);
  const wide = columns.length > 6;
  const frozen2 = wide && (columns[1] === "Fragments" || columns[1] === t("Fragments"));      // the spec's own column name, or that name translated
  const smilesCol = columns.findIndex(c => /smiles/i.test(c));      // a name in the same row shows this structure (inside ESI MS/MS Annotation)
  // the second frozen column sticks right after the first, so it needs the first one's width
  useLayoutEffect(() => {
    const el = wrap.current, first = el?.querySelector("th:first-child");
    if (!frozen2 || !el || !first) return;
    const set = () => el.style.setProperty("--col1-w", `${first.getBoundingClientRect().width}px`);
    set();
    const ro = typeof ResizeObserver !== "undefined" ? new ResizeObserver(set) : null;
    ro?.observe(first);
    return () => ro?.disconnect();
  }, [frozen2, columns, rows]);
  const fix = (j: number) => (frozen2 && j === 1 ? " col-frozen" : "");
  // the spec's columns arrive translated; the columns of a run are English, so the header goes through the dictionary: t(c)
  return (
    <>
    {compact && wide && <small className="table-scroll-hint">{t("Scroll sideways to see all columns.")}</small>}
    <div ref={wrap} className={"table-wrap" + (wide ? " wide" : "") + (frozen2 ? " has-frozen2" : "") + (compact ? " compact-table" : "")} tabIndex={0}
      role="region" aria-label={t("Results table, {count} columns; scroll horizontally for more", { count: columns.length })}>
      <table className="out-table">
        <thead><tr>{columns.map((c, j) => <th key={c} className={fix(j).trim() || undefined}>{t(c)}
          {columnHelp?.[c] && <HelpButton lines={helpLines(columnHelp[c])} className="inline" label={t("Help: {name}", { name: t(c) })} />}</th>)}
          {rowHelp && <th className="help-col" aria-label={t("Help")} />}</tr></thead>
        <tbody>
          {rows.length
            ? rows.map((r, i) => {
              const { cells, help } = Array.isArray(r) ? { cells: r, help: undefined } : r;
              const lines = helpLines(help);
              const rowSmiles = smilesCol >= 0 ? cellText(cells[smilesCol]) : "";
              return (
                <tr key={i}>
                  {cells.map((c, j) => <td key={j} className={(isMultiline(c) ? "long lines" : isLong(c) ? "long" : "") + fix(j) || undefined}>
                    {compact && typeof c === "string" && c.length > 80 ? <CompactText text={c} column={columns[j] ?? "text"} /> : <CellView cell={c} column={columns[j]} rowSmiles={rowSmiles} />}
                  </td>)}
                  {rowHelp && <td className="help-col">{lines.length > 0 && <HelpButton lines={lines} label={t("Row help")} />}</td>}
                </tr>);
            })
            : <tr><td colSpan={columns.length + (rowHelp ? 1 : 0)} className="empty">{empty}</td></tr>}
        </tbody>
      </table>
    </div>
    {view && view.total > rows.length && <p className="result-table-note">{t("First {shown} of {total} rows; all rows are in report.json.", { shown: rows.length, total: view.total.toLocaleString("en-US") })}</p>}
    </>
  );
}

// ---------------------------------------------------------------- one short column: a flow of chips
const cellsOfRow = (r: Row) => (Array.isArray(r) ? r : r.cells);
/** A table with one column of short values (Synonyms) reads better as a flow of chips than as a thin column. */
export function isChipList(columns: string[], rows: Row[]): boolean {
  if (columns.length !== 1 || !rows.length) return false;
  return rows.every(r => {
    const c = cellsOfRow(r)[0];
    if (typeof c === "string") return c.length <= 64 && !c.includes("\n");
    return !!c && typeof c === "object" && "text" in c && typeof c.text === "string" && c.text.length <= 64 && !c.hover_svg;
  });
}

/** The values as chips. A long list shows the first rows and "+N more"; the report shows all of them. */
export function ChipList({ rows, limit = 18, column }: { rows: Row[]; limit?: number; column?: string }) {
  const t = useT();
  const report = useReportMode();
  const [all, setAll] = useState(false);
  const folded = !report && !all && rows.length > limit + 2;           // a fold that hides one or two chips is not worth a button
  const shown = folded ? rows.slice(0, limit) : rows;
  return (
    <div className={"chip-list" + (folded ? " is-folded" : "")}>
      {shown.map((r, i) => <span className="chip-item" key={i}><CellView cell={cellsOfRow(r)[0]} column={column} /></span>)}
      {!report && rows.length > limit + 2 && (
        <button type="button" className="chip-more" aria-expanded={all} onClick={() => setAll(v => !v)}>{all ? t("Show fewer") : t("+{n} more", { n: rows.length - limit })}</button>)}
    </div>
  );
}

// ---------------------------------------------------------------- grouptable
function GroupItem({ group, columns, rowHelp }: { group: RowGroup; columns: string[]; rowHelp: boolean }) {
  const t = useT();
  const [open, setOpen] = useState(false);
  const body = useId();
  // the report: a native <details>, closed like the screen's groups, with the rows inside it so they are in the file and open without a script
  if (useReportMode()) return (
    <details className="gt-group">
      <summary className="gt-head">
        <span className="gt-toggle"><span className="gt-caret" aria-hidden="true" /><span className="gt-title"><RichText text={group.title} /></span></span>
        {group.meta && <span className="gt-meta">{group.meta}</span>}
      </summary>
      <div className="gt-body"><RowsTable columns={columns} rows={group.rows ?? []} rowHelp={rowHelp} empty={t("no rows")} /></div>
    </details>);
  return (
    <div className={"gt-group" + (open ? " open" : "")}>
      {/* the whole header row toggles; the button inside is the keyboard target (its click bubbles to the row) */}
      <div className="gt-head" onClick={() => setOpen(o => !o)}>
        <button type="button" className="gt-toggle" aria-expanded={open} aria-controls={open ? body : undefined}>
          <span className="gt-caret" aria-hidden="true" />
          <span className="gt-title"><RichText text={group.title} /></span>
        </button>
        <HelpButton lines={helpLines(group.help)} label={t("Help: {name}", { name: group.title })} />
        {group.meta && <span className="gt-meta">{group.meta}</span>}
      </div>
      {open && <div id={body} className="gt-body"><RowsTable columns={columns} rows={group.rows ?? []} rowHelp={rowHelp} empty={t("no rows")} /></div>}
    </div>
  );
}

export function GroupTable({ columns, groups, rowHelp = false }: { columns: string[]; groups: RowGroup[]; rowHelp?: boolean }) {
  return <div className="gt">{groups.map((g, i) => <GroupItem key={i} group={g} columns={columns} rowHelp={rowHelp} />)}</div>;
}

// ---------------------------------------------------------------- portal
const STATES = { yes: { glyph: "✓" }, partial: { glyph: "◐" }, no: { glyph: "–" } } as const;

function TileView({ tile }: { tile: Tile }) {
  const t = useT();
  const st = (tile.status in STATES ? tile.status : "no") as keyof typeof STATES;
  const word = st === "yes" ? t("Yes") : st === "partial" ? t("Partial") : t("No");
  return (
    <div className={"tile st-" + st}>
      <div className="tile-top">
        <span className="tile-label">
          <span className="tile-glyph" aria-hidden="true">{STATES[st].glyph}</span>
          <span className="sr-only">{word}: </span>{tile.label}
        </span>
        <HelpButton lines={helpLines(tile.help)} label={t("Help: {name}", { name: tile.label })} />
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
  const t = useT();
  const lines = useReportMode() ? [] : helpLines(chip.help);       // the report has no popovers: a plain chip
  const { triggerProps, popover } = useHelpPopover(lines);
  const lvl = LEVELS[String(chip.level).toLowerCase()] ?? "poor";
  const lvlWord = lvl === "good" ? t("good") : lvl === "possible" ? t("possible") : t("poor");
  const body = (
    <>
      <span className="sr-only">{lvlWord}: </span>
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
  const t = useT();
  const [sel, setSel] = useState(0);
  const cur = options[Math.min(sel, options.length - 1)];
  const meta: ReactNode[] = [];
  if (cur.formula) meta.push(<span key="f" className="rx-formula"><RichText text={cur.formula} /></span>);
  if (cur.mass !== undefined && cur.mass !== null) meta.push(<span key="m"><span className="rx-k">{t("mass")}</span> {asText(cur.mass)}</span>);
  if (cur.heavy !== undefined && cur.heavy !== null) meta.push(<span key="h"><span className="rx-k">{t("heavy mass")}</span> {asText(cur.heavy)}</span>);
  const help = helpLines(cur.help);
  const currentSvg = cur.svg ?? svg;
  const full = currentSvg && cur.product_svg ? { substrate: currentSvg, product: cur.product_svg } : null;
  const reagent = cur.steps.find((s) => s.reagent)?.reagent;
  return (
    <div className={"rx" + (!full && currentSvg ? " with-parent" : "")}>
      {!full && currentSvg && <div className="rx-parent"><ThemedSvg svg={currentSvg} className="" alt={t("Parent structure with reactive groups highlighted for {name}", { name: cur.label })} /></div>}
      <div className="rx-main">
        {options.length > 1 && (
          <div className="rx-opts" role="group" aria-label={t("Reaction options")}>
            {options.map((o, i) => (
              <button key={i} type="button" className={"rx-opt" + (o === cur ? " on" : "")} aria-pressed={o === cur} onClick={() => setSel(i)}>
                <RichText text={o.label} />
              </button>))}
          </div>)}
        {full ? (
          <div className="rx-full">
            <div className="rx-full-structure">
              <div className="rx-full-label">{t("Substrate")}</div>
              <ThemedSvg svg={full.substrate} className="" alt={t("Substrate with reactive sites highlighted for {name}", { name: cur.label })} />
            </div>
            <div className="rx-arrow rx-full-arrow">
              {reagent && <span className="rx-reagent"><span className="sr-only">{t("with")} </span><RichText text={reagent} /></span>}
              <span className="rx-line" aria-hidden="true" />
            </div>
            <div className="rx-full-structure">
              <div className="rx-full-label">{t("Product")}</div>
              <ThemedSvg svg={full.product} className="" alt={t("Complete product structure for {name}", { name: cur.label })} />
            </div>
          </div>
        ) : (
          <div className="rx-steps">
            {(cur.steps ?? []).map((s, i) => (
              <div className="rx-step" key={i}>
                <span className="rx-box"><RichText text={s.from} /></span>
                <span className="rx-arrow">
                  {s.reagent && <span className="rx-reagent"><span className="sr-only">{t("with")} </span><RichText text={s.reagent} /></span>}
                  <span className="rx-line" aria-hidden="true" />
                </span>
                <span className="sr-only"> {t("gives")} </span>
                <span className="rx-box to"><RichText text={s.to} /></span>
                {(s.count ?? 1) > 1 && <span className="rx-count" title={t("{count} groups react", { count: s.count ?? 1 })}>×{s.count}</span>}
              </div>))}
          </div>
        )}
        {(meta.length > 0 || help.length > 0) && (
          <div className="rx-meta">
            {meta.map((m, i) => <span key={i} className="rx-item">{m}{i < meta.length - 1 && <span className="rx-dot" aria-hidden="true">·</span>}</span>)}
            <HelpButton lines={help} label={t("Help: {name}", { name: cur.label })} />
          </div>)}
      </div>
    </div>
  );
}
