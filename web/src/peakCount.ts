/** Peaks in a pasted / uploaded spectrum text, counted the way the server reads it (mzdesk/shared/spectrum.py: parse_spectrum on first_record).
 *  Used for the short "MS2 spectrum: 27 peaks" in the collapsed input bar, so it has to agree with what the run will use:
 *   - only the first record of an MGF (BEGIN IONS … END IONS), MSP ("Name:"), SIRIUS (">compound") or MassBank (ends with "//") text;
 *   - a MassBank record counts only the lines under "PK$PEAK:"; in other text, lines that do not start with a number are headers
 *     ("Name: PG 36:2", "PEPMASS=…", "Num Peaks: 13") and never peaks;
 *   - a line that starts with a number is one peak ("mz", "mz int", tab / space / comma separated, any extra columns), or one peak per
 *     "mz:int" token (several per line: "54.03:7 56.04:88"); ";" separates peaks on such a line ("138.0662 100;195.0877 32"), except that
 *     a row of numbers only ("195.0877;100", "195.0877;100;[M+H]+") is one peak with ";" columns (m/z; intensity; …);
 *   - a comma separates the columns only when nothing else does ("195.0877,100", "41,100"); next to a tab, space or ";" it is a decimal comma
 *     ("138,0662<TAB>100") or, after a dot-decimal m/z, a thousands separator ("100.0 1,200"), as in spectrum.py (_row_numbers);
 *   - m/z <= 0 or intensity <= 0 is dropped and a repeated m/z (5 decimals) counts once.
 *  JSON, XML and "mz = [...]" array texts are not counted here (0): the caller then says the spectrum is given, without a number. */

const NUMBER = /^[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$/;
const STARTS_WITH_NUMBER = /^\s*[\d.+-]/;
const COLON_GAP = /([\d.])\s*:\s*(?=[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?(?:[\s,]|$))/g;     // "184.07 : 100" -> "184.07:100"
const MGF_BEGIN = /^\s*BEGIN\s+IONS\s*$/i;
const MGF_END = /^\s*END\s+IONS\s*$/i;
const MSP_NAME = /^\s*(?:NAME\s*:|>compound\b)/i;
const MB_END = /^\s*\/\/\s*$/;
const PK_PEAK = /^\s*PK\$PEAK\s*:/i;
const NOT_LINE_BASED = /^\s*[[{<]|\b(?:m\/?z)\w*["']?\s*[=:]\s*(?:(?:np\.)?array\()?\s*\[/i;

const finite = (tok: string | undefined): number | null => (tok !== undefined && NUMBER.test(tok) && Number.isFinite(Number(tok)) ? Number(tok) : null);

function firstRecord(lines: string[]): string[] {
  if (lines.some(l => MGF_BEGIN.test(l))) {
    let cur: string[] | null = null;
    for (const l of lines) {
      if (MGF_BEGIN.test(l)) cur = [];
      else if (MGF_END.test(l)) { if (cur) return cur; }
      else if (cur) cur.push(l);
    }
    return cur ?? [];                                                  // an unterminated last record still counts
  }
  const cur: string[] = [];
  for (const l of lines) {
    if (MSP_NAME.test(l) && cur.some(x => x.trim())) break;
    cur.push(l);
    if (MB_END.test(l)) break;
  }
  return cur;
}

function massBankPairs(lines: string[]): [number, number][] {
  const out: [number, number][] = [];
  let inside = false;
  for (const l of lines) {
    if (PK_PEAK.test(l)) inside = true;
    else if (inside && STARTS_WITH_NUMBER.test(l)) {
      const toks = l.trim().split(/\s+/);
      const mz = finite(toks[0]), inten = finite(toks[1]);
      if (mz !== null && inten !== null) out.push([mz, inten]);
    } else if (inside) inside = false;
  }
  return out;
}

const PLAIN = NUMBER;
const COMMA_NUM = /^[-+]?\d+,\d+$/;                                  // 138,0662 (decimal comma), 1,200 (thousands) or "mz,int" when nothing else separates the columns
const THOUSANDS = /^[-+]?\d{1,3}(?:,\d{3})+(?:\.\d+)?$/;               // 1,200  1,234,567  1,234.5
const GROUP3 = /^[-+]?\d{1,3},\d{3}(?:\.\d+)?$/;                      // 1,200  1,234.5: a thousands group, or a decimal comma with three decimals, or a csv "mz,int" (41,100.5)
const ZERO_FRAC = /^[-+]?\d+,0\d+$/;                                  // 195,0877: nobody writes an intensity 0877, so this comma is a decimal one

const certainThousands = (c: string) => THOUSANDS.test(c) && c.split(",").length > 2;                   // 1,234,567
const isPlain = (c: string) => PLAIN.test(c) || certainThousands(c);                                     // a number whose commas (if any) cannot be decimal commas
const cellsOf = (text: string) => text.split(/\s+/).filter(c => c.replace(/^,+|,+$/g, ""));            // tab / space separated cells, without "," on its own
const startsNumber = (cells: string[]) => cells.length > 0 && STARTS_WITH_NUMBER.test(cells[0].replace(/^,+/, ""));

function provesDecimalComma(cells: string[]): boolean {
  const commas = cells.filter(c => COMMA_NUM.test(c));
  if (!commas.length) return false;
  if (commas.some(c => ZERO_FRAC.test(c) && !GROUP3.test(c))) return true;
  if (!cells.some(isPlain)) return false;
  return COMMA_NUM.test(cells[0]) || commas.some(c => !GROUP3.test(c));
}

function rowNumbers(cells: string[], decimal: boolean): string[] {
  const plain = cells.some(isPlain), commaMz = COMMA_NUM.test(cells[0]);
  const out: string[] = [];
  cells.forEach((cell, i) => {
    if (!cell.includes(",")) out.push(cell);
    else if (certainThousands(cell) || (i > 0 && plain && !(decimal || commaMz) && GROUP3.test(cell))) out.push(cell.replace(/,/g, ""));
    else if (!COMMA_NUM.test(cell)) out.push(...cell.split(",").filter(Boolean));                       // "195.0877,100": a csv pair
    else if (plain || decimal) out.push(cell.replace(/,/g, "."));
    else out.push(...cell.split(","));                                                                    // "41,100": mz,intensity
  });
  return out;
}

type Row = { kind: "colon"; text: string } | { kind: "cells" | "cols"; cells: string[] };

function linePairs(lines: string[]): [number, number][] {
  const out: [number, number][] = [];
  const nNumeric = lines.filter(l => STARTS_WITH_NUMBER.test(l)).length;
  const rows: Row[] = [];
  for (const raw of lines) {                                           // ";" separates peaks only on lines that start with a number
    const parts = STARTS_WITH_NUMBER.test(raw) ? raw.split(";").map(p => p.trim()) : [raw];
    while (parts.length && !parts[parts.length - 1]) parts.pop();
    // "195.0877;100" and "195.0877;100;[M+H]+" are csv rows with ";" columns. A lone line of numbers only is a ";" separated list of m/z.
    const numeric = parts.map(p => PLAIN.test(p) || COMMA_NUM.test(p) || THOUSANDS.test(p));
    if (parts.length > 1 && numeric[0] && numeric[1] && (!numeric.slice(2).some(Boolean) || (numeric.every(Boolean) && nNumeric > 1))) {
      rows.push({ kind: "cols", cells: parts });
      continue;
    }
    for (const part of parts) {
      const line = part.trim().replace(COLON_GAP, "$1:");
      rows.push(line.includes(":") ? { kind: "colon", text: line } : { kind: "cells", cells: cellsOf(line) });
    }
  }
  const decimal = rows.some(r => r.kind !== "colon" && startsNumber(r.cells) && provesDecimalComma(r.cells));
  for (const row of rows) {
    if (row.kind === "colon") {
      const toks = row.text.split(/[\s,]+/).filter(Boolean);
      if (!toks.length || !/^[\d.+-]/.test(toks[0])) continue;
      const pairs: [number, number][] = [];
      for (const tok of toks) {
        const i = tok.indexOf(":");
        if (i < 0) continue;
        const mz = finite(tok.slice(0, i)), rest = tok.slice(i + 1);
        const inten = rest === "" ? 1 : finite(rest);
        if (mz !== null && inten !== null) pairs.push([mz, inten]);
      }
      if (!pairs.length) {                                             // a plain "mz [intensity] [extra…]" line
        const mz = finite(toks[0]), inten = toks.length > 1 ? finite(toks[1]) : 1;
        if (mz !== null && inten !== null) pairs.push([mz, inten]);
      }
      out.push(...pairs);
      continue;
    }
    // "41,100;55,85" is two "mz,int" peaks, not m/z 41.1 with intensity 55.85
    const split = row.kind === "cols" && !decimal && !row.cells.some(isPlain);
    for (const cells of split ? row.cells.map(cellsOf) : [row.cells]) {
      if (!startsNumber(cells)) continue;
      const nums = rowNumbers(cells, decimal);
      const mz = finite(nums[0]), inten = nums.length > 1 ? finite(nums[1]) : 1;
      if (mz !== null && inten !== null) out.push([mz, inten]);
    }
  }
  return out;
}

export function countPeaks(text: string): number {
  const s = (text ?? "").trim();
  if (!s || NOT_LINE_BASED.test(s)) return 0;
  const lines = firstRecord(s.split(/\r\n|\r|\n/));
  const pairs = lines.some(l => PK_PEAK.test(l)) ? massBankPairs(lines) : linePairs(lines);
  return new Set(pairs.filter(([mz, inten]) => mz > 0 && inten > 0).map(([mz]) => Math.round(mz * 1e5))).size;
}
