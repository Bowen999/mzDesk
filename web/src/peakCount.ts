/** Peaks in a pasted / uploaded spectrum text, counted the way the server reads it (mzdesk/shared/spectrum.py: parse_spectrum on first_record).
 *  Used for the short "MS2 spectrum: 27 peaks" in the collapsed input bar, so it has to agree with what the run will use:
 *   - only the first record of an MGF (BEGIN IONS … END IONS), MSP ("Name:"), SIRIUS (">compound") or MassBank (ends with "//") text;
 *   - a MassBank record counts only the lines under "PK$PEAK:"; in other text, lines that do not start with a number are headers
 *     ("Name: PG 36:2", "PEPMASS=…", "Num Peaks: 13") and never peaks;
 *   - a line that starts with a number is one peak ("mz", "mz int", tab / space / comma separated, any extra columns), or one peak per
 *     "mz:int" token (several per line: "54.03:7 56.04:88"); ";" separates peaks on such a line;
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

function linePairs(lines: string[]): [number, number][] {
  const out: [number, number][] = [];
  const segments: string[] = [];
  for (const raw of lines) segments.push(...(STARTS_WITH_NUMBER.test(raw) ? raw.split(";") : [raw]));
  for (const seg of segments) {
    const toks = seg.trim().replace(COLON_GAP, "$1:").split(/[\s,]+/).filter(Boolean);
    if (!toks.length || !/^[\d.+-]/.test(toks[0])) continue;
    const pairs: [number, number][] = [];
    for (const tok of toks) {
      const i = tok.indexOf(":");
      if (i < 0) continue;
      const mz = finite(tok.slice(0, i)), rest = tok.slice(i + 1);
      const inten = rest === "" ? 1 : finite(rest);
      if (mz !== null && inten !== null) pairs.push([mz, inten]);
    }
    if (!pairs.length) {                                               // a plain "mz [intensity] [extra…]" line
      const mz = finite(toks[0]), inten = toks.length > 1 ? finite(toks[1]) : 1;
      if (mz !== null && inten !== null) pairs.push([mz, inten]);
    }
    out.push(...pairs);
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
