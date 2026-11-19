import { useEffect, useId, useLayoutEffect, useMemo, useRef, useState, type ReactNode } from "react";
import type React from "react";
import { useT } from "../i18n";
import { useReportMode } from "../report/mode";
import type { Cell, Help, Row, RowGroup, RxOption, StatusChip, Tile } from "../types";
import { FragmentButton } from "./FragmentDialog";
import { ThemedSvg } from "./Figure";
import { HelpButton, RichText, helpLines, useHelpPopover } from "./Help";
import { sameTarget, useReferenceSelection } from "./referenceContext";
import { StructureHover, structureForCell, useStructures } from "./StructureHover";

/**
 * Renderers for the generic output blocks. Data contracts (the spec helpers in mzdesk/spec.py and the builders in
 * mzdesk/api/blocks.py are the reference; help = string | string[] everywhere):
 *   table / grouptable cell : "text" | {"text", "url"} | {"badge", "tone"} | {"text", "help"} | {"text", "hover_svg"} | {"img"} | {"fragment": {smiles, adduct, name?}}
 *                             (badge / text may also carry url + help; a fragment cell is a "View" button, its spectrum is the table's data.fragments)
 *   table row               : [cells] | {"cells": [...], "help", "reference"?}   (the "?" column only with spec row_help=True;
 *                             reference = a library record: a click on the row shows it in the Selected spectrum block, see ReferenceSpectrum.tsx)
 *   table data              : {"rows": [...], "count"?}           count -> pill next to the title
 *   grouptable data         : {"groups": [{"title", "meta"?, "url"?, "help"?, "rows": [...]}], "count"?}   collapsed by default; url = a small link after the title
 *                             (spec page_size = pages of groups + filter + "Expand all"; compact as for a table)
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

// ---------------------------------------------------------------- numeric columns (right-aligned)
/** A plain number as a table shows it: "20", "0.800", "-1.5", "1,234", "12.5%", ">1000", "<0.01", "≥2", "1e-3". Not an ID, an adduct ("M+H"), a range ("100-200") or text. */
const NUMBER_TEXT = /^[<>≤≥~≈±+\-−]?\s?(?:(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?|\.\d+)(?:[eE][+-]?\d+)?\s?%?$/;       // the mantissa has a digit: "E1" is not a number
const isNumberText = (text: string) => /\d/.test(text) && NUMBER_TEXT.test(text.trim());
const isBlankText = (text: string) => text.trim() === "" || text.trim() === "—";

/** What a cell is for the alignment of its column: "empty" (nothing, "", "—"), "number" (a plain number, or an object cell whose text is one) or "other". */
export function cellKind(c: Cell | number | undefined): "empty" | "number" | "other" {
  if (c === null || c === undefined) return "empty";
  if (typeof c === "number") return Number.isFinite(c) ? "number" : "other";
  if (typeof c === "string") return isBlankText(c) ? "empty" : isNumberText(c) ? "number" : "other";
  if (typeof c !== "object") return "other";
  if ("img" in c || "fragment" in c || c.badge !== undefined) return "other";
  const text = c.text ?? "";
  return isBlankText(text) ? "empty" : isNumberText(text) ? "number" : "other";
}

/** Which columns are numeric: every non-empty cell of the column is a plain number, and there is at least one. `rows` may come from several groups, so that
 *  the columns of a group table line up. */
export function numericColumns(columnCount: number, rows: Row[]): boolean[] {
  const kinds = Array.from({ length: columnCount }, () => ({ number: false, other: false }));
  for (const r of rows) {
    const cells = Array.isArray(r) ? r : r.cells;
    for (let j = 0; j < columnCount; j++) {
      const kind = cellKind(cells[j]);
      if (kind !== "empty") kinds[j][kind] = true;
    }
  }
  return kinds.map(k => k.number && !k.other);
}

/** A cell's text up to this long (and on one line) is never wrapped: a table this wide scrolls sideways anyway, so wrapping a short label only doubles the row. */
export const NOWRAP_MAX = 24;
const shownText = (c: Cell | number | undefined) => {
  if (c === null || c === undefined) return "";
  if (typeof c === "number") return String(c);
  if (typeof c === "string") return c;
  if (typeof c !== "object") return String(c);
  if ("img" in c || "fragment" in c) return "";
  return c.badge !== undefined ? String(c.badge) : c.text ?? "";
};

/** Which columns never wrap: numeric ones, and those whose every non-empty cell is short (at most NOWRAP_MAX characters) and on one line. Decided over `rows`, which
 *  should be every row of the table (or of all its groups), so paging, filtering and folding never change a column. */
export function nowrapColumns(columnCount: number, rows: Row[], numeric: boolean[] = numericColumns(columnCount, rows)): boolean[] {
  const short = Array.from({ length: columnCount }, () => true);
  for (const r of rows) {
    const cells = Array.isArray(r) ? r : r.cells;
    for (let j = 0; j < columnCount; j++) {
      if (!short[j]) continue;
      const text = shownText(cells[j]).trim();
      if (text.length > NOWRAP_MAX || text.includes("\n")) short[j] = false;
    }
  }
  return short.map((ok, j) => ok || !!numeric[j]);
}
const classes = (...names: (string | false | undefined)[]) => names.filter(Boolean).join(" ") || undefined;

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
    columns: keep(columns), kept: keep(columns.map((_, j) => j)), total: rows.length,
    rows: rows.slice(0, max).map(r => (Array.isArray(r) ? keep(r) : { ...r, cells: keep(r.cells) })),
  };
}

/** The shared table body: header + rows (+ trailing "?" column when `rowHelp`; `columnHelp` puts a small "?" after a column's name).
 *  A wide table freezes its first column (theme.css); a "Fragments" column right after it is frozen too (fragment-view.css), so the View button
 *  stays in reach next to the compound while the rest scrolls. */
export function RowsTable({ columns: allColumns, rows: allRows, rowHelp: withHelp = false, empty, compact = false, columnHelp, numeric: givenNumeric, nowrap: givenNowrap }:
  { columns: string[]; rows: Row[]; rowHelp?: boolean; empty: string; compact?: boolean; columnHelp?: Record<string, Help>; numeric?: boolean[]; nowrap?: boolean[] }) {
  const t = useT();
  const report = useReportMode();
  const rowHelp = withHelp && !report;       // no "?" column in the report: there are no popovers
  const view = useMemo(() => (report ? reportView(allColumns, allRows, report.maxRows) : null), [report, allColumns, allRows]);
  const columns = view?.columns ?? allColumns, rows = view?.rows ?? allRows;
  // a column whose every filled cell is a number is right-aligned (header too); a group table decides it over all its groups (`numeric`, by the columns as given)
  const numeric = useMemo(
    () => (givenNumeric ? (view ? view.kept.map(j => givenNumeric[j]) : givenNumeric) : numericColumns(columns.length, rows)),
    [givenNumeric, view, columns.length, rows]);
  // a numeric or short-text column keeps each cell on one line (see nowrapColumns); like `numeric`, a table that is paged or grouped passes it, decided over all its rows
  const nowrap = useMemo(
    () => (givenNowrap ? (view ? view.kept.map(j => givenNowrap[j]) : givenNowrap) : nowrapColumns(columns.length, rows, numeric)),
    [givenNowrap, view, columns.length, rows, numeric]);
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
  const fix = (j: number) => (frozen2 && j === 1 ? "col-frozen" : false);
  const selection = useReferenceSelection();
  // a row with a `reference` shows that record in the Selected spectrum block on a click (or Enter), unless the click was meant for a link or
  // button in it or ended a text selection; the record shown is highlighted
  const picker = (r: Row) => {
    const target = !Array.isArray(r) ? r.reference : undefined;
    if (!target || !selection) return {};
    const on = sameTarget(selection.current, target);
    return {
      className: "row-open" + (on ? " is-shown" : ""), tabIndex: 0, "aria-current": on ? ("true" as const) : undefined, title: t("Show this spectrum"),
      onClick: (e: React.MouseEvent) => {
        if ((e.target as Element).closest("a, button, input") || window.getSelection()?.toString()) return;
        selection.select(target);
      },
      onKeyDown: (e: React.KeyboardEvent) => {
        if (e.target === e.currentTarget && (e.key === "Enter" || e.key === " ")) { e.preventDefault(); selection.select(target); }
      },
    };
  };
  // the spec's columns arrive translated; the columns of a run are English, so the header goes through the dictionary: t(c)
  return (
    <>
    {compact && wide && <small className="table-scroll-hint">{t("Scroll sideways to see all columns.")}</small>}
    <div ref={wrap} className={"table-wrap" + (wide ? " wide" : "") + (frozen2 ? " has-frozen2" : "") + (compact ? " compact-table" : "")} tabIndex={0}
      role="region" aria-label={t("Results table, {count} columns; scroll horizontally for more", { count: columns.length })}>
      <table className="out-table">
        <thead><tr>{columns.map((c, j) => <th key={c} className={classes(fix(j), numeric[j] && "num")}>{t(c)}
          {columnHelp?.[c] && <HelpButton lines={helpLines(columnHelp[c])} className="inline" label={t("Help: {name}", { name: t(c) })} />}</th>)}
          {rowHelp && <th className="help-col" aria-label={t("Help")} />}</tr></thead>
        <tbody>
          {rows.length
            ? rows.map((r, i) => {
              const { cells, help } = Array.isArray(r) ? { cells: r, help: undefined } : r;
              const lines = helpLines(help);
              const rowSmiles = smilesCol >= 0 ? cellText(cells[smilesCol]) : "";
              return (
                <tr key={i} {...picker(r)}>
                  {cells.map((c, j) => <td key={j} className={classes(isMultiline(c) ? "long lines" : isLong(c) && "long", fix(j), numeric[j] && "num", nowrap[j] && "nowrap")}>
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

// ---------------------------------------------------------------- table tools (shared by the paged table and the paged group table)
/** The text a filter searches in one cell: what the cell shows (an image or a fragment button has no text). */
export function cellSearchText(cell: Cell): string {
  if (cell === null || cell === undefined) return "";
  if (typeof cell !== "object") return String(cell);
  return "img" in cell || "fragment" in cell ? "" : String(cell.text ?? cell.badge ?? "");
}
/** All the text of one row, lower-cased, for the filter. */
export const rowSearchText = (row: Row) => (Array.isArray(row) ? row : row.cells).map(cellSearchText).join(" ").toLocaleLowerCase();

/** What the filter of a paged group table searches in one group: its title (without markup), its meta line and the text of all its rows; lower-cased. */
export const groupSearchText = (g: RowGroup) => [g.title.replace(/<[^>]*>/g, ""), g.meta ?? "", ...(g.rows ?? []).map(rowSearchText)].join(" ").toLocaleLowerCase();

/** The filter box of a paged result (the label says what is searched). */
export function FilterField({ label, placeholder, ariaLabel, value, onChange }: { label: string; placeholder: string; ariaLabel: string; value: string; onChange: (v: string) => void }) {
  return (
    <label className="result-table-filter"><span>{label}</span>
      <input type="search" value={value} placeholder={placeholder} aria-label={ariaLabel} onChange={e => onChange(e.target.value)} />
    </label>);
}

/** "Page 2 of 5" with Previous / Next. `current` counts from 0. */
export function PageNav({ current, pages, label, onPage }: { current: number; pages: number; label: string; onPage: (page: number) => void }) {
  const t = useT();
  return (
    <nav className="result-table-pages" aria-label={label}>
      <span aria-live="polite">{t("Page {current} of {pages}", { current: current + 1, pages })}</span>
      <div><button type="button" className="page-button" disabled={current === 0} onClick={() => onPage(current - 1)}>{t("Previous")}</button>
        <button type="button" className="page-button" disabled={current >= pages - 1} onClick={() => onPage(current + 1)}>{t("Next")}</button></div>
    </nav>);
}

// ---------------------------------------------------------------- grouptable
const ExternalIcon = () => (
  <svg viewBox="0 0 16 16" width="13" height="13" fill="none" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true" focusable="false">
    <path d="M9 2.5h4.5V7M13.5 2.5 7.5 8.5M6.5 3.5h-3a1 1 0 0 0-1 1v8a1 1 0 0 0 1 1h8a1 1 0 0 0 1-1v-3" />
  </svg>);

/** The text of a group title without markup. */
const plainTitle = (title: string) => title.replace(/<[^>]*>/g, "").trim();

/** What names a group's dataset link: the part of the title before " · " (an accession) or the whole plain title. */
export const groupLinkName = (title: string) => { const text = plainTitle(title); return text.split(" · ")[0].trim() || text; };


/** The small link after a group's title (a dataset's own page). Its click must not toggle the group. */
function GroupLink({ url, title }: { url?: string; title: string }) {
  const t = useT();
  const href = url ? safeUrl(url) : undefined;
  if (!href) return null;
  const label = t("Open dataset page: {name}", { name: groupLinkName(title) });
  return <a className="gt-link" href={href} target="_blank" rel="noopener noreferrer" title={label} aria-label={label} onClick={e => e.stopPropagation()}><ExternalIcon /></a>;
}

/** The meta line of a group: always on a line of its own under the title, left-aligned past the caret, wrapping (outputs-v2.css). */
const GroupMeta = ({ group }: { group: RowGroup }) => group.meta ? <span className="gt-meta">{group.meta}</span> : null;

const isEmptyCell = (c: Cell | undefined) => {
  if (c === null || c === undefined) return true;
  if (typeof c === "string") return c.trim() === "" || c === "—";
  if (typeof c !== "object") return false;
  if ("img" in c || "fragment" in c) return false;
  if (c.badge !== undefined || c.url || c.href || c.hover_svg) return false;
  const text = (c.text ?? "").trim();
  return text === "" || text === "—";
};

/** The columns of a group table that have something in at least one row of ANY group (so the columns line up across groups). A column that is empty or "—" in every
 *  row of every group is left out; if that were every column, all stay. Returns the kept column indexes. */
export function usedColumns(columns: string[], groups: RowGroup[]): number[] {
  const all = columns.map((_, j) => j);
  const rows = groups.flatMap(g => g.rows ?? []);
  if (!rows.length) return all;
  const keep = all.filter(j => rows.some(r => !isEmptyCell((Array.isArray(r) ? r : r.cells)[j])));
  return keep.length ? keep : all;
}

/** The groups and columns without the columns that are empty everywhere (the same object when nothing is dropped). */
export function trimEmptyColumns(columns: string[], groups: RowGroup[]): { columns: string[]; groups: RowGroup[] } {
  const keep = usedColumns(columns, groups);
  if (keep.length === columns.length) return { columns, groups };
  const pick = <T,>(cells: T[]) => keep.map(j => cells[j]);
  return {
    columns: pick(columns),
    groups: groups.map(g => ({ ...g, rows: (g.rows ?? []).map(r => (Array.isArray(r) ? pick(r) : { ...r, cells: pick(r.cells) })) })),
  };
}

function GroupItem({ group, columns, numeric, nowrap, rowHelp, compact, open, onToggle }: { group: RowGroup; columns: string[]; numeric: boolean[]; nowrap: boolean[]; rowHelp: boolean; compact?: boolean; open: boolean; onToggle: () => void }) {
  const t = useT();
  const body = useId();
  // the report: a native <details>, closed like the screen's groups, with the rows inside it so they are in the file and open without a script
  if (useReportMode()) return (
    <details className="gt-group">
      <summary className="gt-head">
        <span className="gt-toggle"><span className="gt-caret" aria-hidden="true" /><span className="gt-title"><RichText text={group.title} /></span></span>
        <GroupLink url={group.url} title={group.title} />
        <GroupMeta group={group} />
      </summary>
      <div className="gt-body"><RowsTable columns={columns} rows={group.rows ?? []} rowHelp={rowHelp} empty={t("no rows")} compact={compact} numeric={numeric} nowrap={nowrap} /></div>
    </details>);
  return (
    <div className={"gt-group" + (open ? " open" : "")}>
      {/* the whole header row toggles; the button inside is the keyboard target (its click bubbles to the row) */}
      <div className="gt-head" onClick={onToggle}>
        <button type="button" className="gt-toggle" aria-expanded={open} aria-controls={open ? body : undefined}>
          <span className="gt-caret" aria-hidden="true" />
          <span className="gt-title"><RichText text={group.title} /></span>
        </button>
        <GroupLink url={group.url} title={group.title} />
        <HelpButton lines={helpLines(group.help)} label={t("Help: {name}", { name: group.title })} />
        <GroupMeta group={group} />
      </div>
      {open && <div id={body} className="gt-body"><RowsTable columns={columns} rows={group.rows ?? []} rowHelp={rowHelp} empty={t("no rows")} compact={compact} numeric={numeric} nowrap={nowrap} /></div>}
    </div>
  );
}

/** Groups, folded by default. With `pageSize` and more groups than that, they come in pages with a text filter (title, meta and the text of the rows);
 *  "Expand all" / "Collapse all" then works on the groups of the current page. Without `pageSize` nothing is added. The report shows every group. */
export function GroupTable({ columns: allColumns, groups: allGroups, rowHelp = false, pageSize = 0, compact, name }: { columns: string[]; groups: RowGroup[]; rowHelp?: boolean; pageSize?: number; compact?: boolean; name?: string | null }) {
  const t = useT();
  // a column empty (or "—") in every row of every group is not shown
  const { columns, groups } = useMemo(() => trimEmptyColumns(allColumns, allGroups), [allColumns, allGroups]);
  const { numeric, nowrap } = useMemo(() => {       // decided over every row of every group, so the columns line up
    const every = groups.flatMap(g => g.rows ?? []);
    const numeric = numericColumns(columns.length, every);
    return { numeric, nowrap: nowrapColumns(columns.length, every, numeric) };
  }, [columns, groups]);
  const report = useReportMode();
  const size = pageSize > 0 ? pageSize : 0;
  const [open, setOpen] = useState<Set<number>>(() => new Set());
  const [page, setPage] = useState(0);
  const [filter, setFilter] = useState("");
  useEffect(() => { setOpen(new Set()); setPage(0); setFilter(""); }, [groups]);
  // lower-cased text of each group, once per data (a run can have 200 groups of many rows)
  const haystacks = useMemo(() => size ? groups.map(groupSearchText) : [], [groups, size]);
  const toggle = (i: number) => setOpen(prev => { const next = new Set(prev); if (!next.delete(i)) next.add(i); return next; });
  if (report || !size) {
    return <div className="gt">{groups.map((g, i) => <GroupItem key={i} group={g} columns={columns} numeric={numeric} nowrap={nowrap} rowHelp={rowHelp} compact={compact} open={!report && open.has(i)} onToggle={() => toggle(i)} />)}</div>;
  }
  const needle = filter.trim().toLocaleLowerCase();
  const matching = groups.map((_, i) => i).filter(i => !needle || haystacks[i].includes(needle));
  const paged = groups.length > size;
  const pages = Math.max(1, Math.ceil(matching.length / size));
  const current = Math.min(page, pages - 1);
  const start = current * size;
  const shown = matching.slice(start, start + size);
  const allOpen = shown.length > 0 && shown.every(i => open.has(i));
  const range = { from: start + 1, to: Math.min(start + size, matching.length), total: matching.length };
  return (
    <div className="result-table">
      {(paged || groups.length > 1) && <div className="result-table-tools">
        {paged && <FilterField label={t("Filter groups")} placeholder={t("Search these groups")} ariaLabel={name ? t("Filter {name} groups", { name }) : t("Filter groups")}
          value={filter} onChange={v => { setFilter(v); setPage(0); }} />}
        {groups.length > 1 && <button type="button" className="page-button gt-expand" disabled={!shown.length}
          onClick={() => setOpen(prev => { const next = new Set(prev); for (const i of shown) allOpen ? next.delete(i) : next.add(i); return next; })}>
          {allOpen ? t("Collapse all") : t("Expand all")}</button>}
        {paged && <span className="result-table-count" aria-live="polite">
          {matching.length ? (needle ? t("{from}–{to} of {total} matching groups", range) : t("{from}–{to} of {total} groups", range)) : t("No matching groups")}
        </span>}
      </div>}
      {shown.length
        ? <div className="gt">{shown.map(i => <GroupItem key={i} group={groups[i]} columns={columns} numeric={numeric} nowrap={nowrap} rowHelp={rowHelp} compact={compact} open={open.has(i)} onToggle={() => toggle(i)} />)}</div>
        : <div className="placeholder-box">{t("No groups match this filter.")}</div>}
      {paged && <PageNav current={current} pages={pages} label={name ? t("{name} pages", { name }) : t("Group pages")} onPage={setPage} />}
    </div>
  );
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
