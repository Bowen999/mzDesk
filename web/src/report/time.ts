/** Report timestamps are always Edmonton time (America/Edmonton), whatever the browser's own zone is. */
export const REPORT_TZ = "America/Edmonton";

export type EdmontonTime = {
  stamp: string;      // 20261031-153012
  abbr: string;       // MDT / MST
  offset: string;     // -06:00
  iso: string;        // 2026-10-31T15:30:12-06:00
  human: string;      // 2026-10-31 15:30:12 MDT (UTC-06:00)
};

const pad = (n: number) => String(n).padStart(2, "0");

export function edmonton(d: Date = new Date()): EdmontonTime {
  const parts = new Intl.DateTimeFormat("en-US", {
    timeZone: REPORT_TZ, year: "numeric", month: "2-digit", day: "2-digit",
    hour: "2-digit", minute: "2-digit", second: "2-digit", hourCycle: "h23", timeZoneName: "short",
  }).formatToParts(d);
  const get = (t: string) => parts.find(p => p.type === t)?.value ?? "";
  const [Y, M, D, h, m, s] = ["year", "month", "day", "hour", "minute", "second"].map(t => Number(get(t)));
  // UTC offset = (wall clock read as UTC) - (the real instant), in whole minutes
  const mins = Math.round((Date.UTC(Y, M - 1, D, h, m, s) - Math.floor(d.getTime() / 1000) * 1000) / 60000);
  const sign = mins < 0 ? "-" : "+", a = Math.abs(mins);
  const offset = `${sign}${pad(Math.floor(a / 60))}:${pad(a % 60)}`;
  const abbr = /^[A-Z]{2,5}$/.test(get("timeZoneName")) ? get("timeZoneName") : `UTC${offset}`;
  const date = `${Y}-${pad(M)}-${pad(D)}`, clock = `${pad(h)}:${pad(m)}:${pad(s)}`;
  return {
    stamp: `${Y}${pad(M)}${pad(D)}-${pad(h)}${pad(m)}${pad(s)}`,
    abbr, offset,
    iso: `${date}T${clock}${offset}`,
    human: `${date} ${clock} ${abbr} (UTC${offset})`,
  };
}

const ALPHABET = "abcdefghjkmnpqrstuvwxyz23456789";      // no i, l, o, 0, 1: easy to read out loud

/** Short random suffix (5 characters). */
export function shortRandom(n = 5): string {
  const bytes = new Uint8Array(n);
  crypto.getRandomValues(bytes);
  return Array.from(bytes, b => ALPHABET[b % ALPHABET.length]).join("");
}

/** mzdesk-20261031-153012-MDT-k7q2m: Edmonton timestamp + zone abbreviation + short random string. */
export function newRunId(t: EdmontonTime = edmonton()): string {
  return `mzdesk-${t.stamp}-${t.abbr.replace(/[^A-Za-z0-9+-]/g, "")}-${shortRandom()}`;
}
