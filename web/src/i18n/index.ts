import { useSyncExternalStore } from "react";
import type { Spec } from "../types";
import specRules from "./spec-rules.json";

/**
 * Interface language. English is the source text and needs no dictionary: every visible string in the code is written in English as the
 * literal first argument of t() ("Run", "Download the results of {n} pages") and that text is the key into locales/<lang>.json.
 *   useT()         hook: the t function; the component draws again when the language or its dictionary changes
 *   t(text, vars)  the same lookup outside components (network.ts), with the language at call time
 *   t(text, { context: "theme" })  the same English word with another meaning: the dictionary entry is "text|context" ("System|theme"), the English
 *                  shown is `text`, and a language without that entry uses the plain "text". Keep `context` the first key (extract-i18n.mjs reads it).
 *   localizeSpec   the page list (/api/spec) with its display text translated; ids, option lists, values and flags stay as they are
 * A missing translation shows the English text. Dictionaries are loaded by dynamic import, so the main bundle carries none of them; until the
 * one for the chosen language arrives, English shows.
 * Default language: the browser's (navigator.languages) when it is one of the six, else English. A pick is remembered in localStorage.
 * source/ui.json (scripts/extract-i18n.mjs) and source/spec.json (python -m mzdesk.spec_i18n) list what the dictionaries must cover.
 */
export type Lang = "en" | "zh" | "ja" | "de" | "fr" | "es";
export const LANGS: { id: Lang; short: string; name: string }[] = [
  { id: "en", short: "EN", name: "English" },
  { id: "zh", short: "中", name: "中文" },
  { id: "ja", short: "日", name: "日本語" },
  { id: "de", short: "DE", name: "Deutsch" },
  { id: "fr", short: "FR", name: "Français" },
  { id: "es", short: "ES", name: "Español" },
];

type Dict = Record<string, string>;
type Vars = Record<string, string | number>;
export type TFn = (text: string, vars?: Vars) => string;

const KEY = "mzdesk.lang";
const HTML_LANG: Record<Lang, string> = { en: "en", zh: "zh-CN", ja: "ja", de: "de", fr: "fr", es: "es" };
const LOADERS: Record<Exclude<Lang, "en">, () => Promise<{ default: Dict }>> = {
  zh: () => import("./locales/zh.json"),
  ja: () => import("./locales/ja.json"),
  de: () => import("./locales/de.json"),
  fr: () => import("./locales/fr.json"),
  es: () => import("./locales/es.json"),
};

const isLang = (v: unknown): v is Lang => LANGS.some(l => l.id === v);

/** The first of the browser's languages we have a dictionary for (zh-TW, zh-Hans ... all count as zh), else English. */
export function pickLanguage(tags: readonly string[]): Lang {
  for (const tag of tags) {
    const primary = String(tag).toLowerCase().split("-")[0];
    if (isLang(primary)) return primary;
  }
  return "en";
}

/** The saved pick, else the browser's language. Only a real browser decides: under Node (the tests, server rendering) the answer is English. */
function initialLang(): Lang {
  try { const v = localStorage.getItem(KEY); if (isLang(v)) return v; } catch { /* no storage: fall through */ }
  if (typeof window === "undefined" || typeof navigator === "undefined") return "en";
  return pickLanguage(navigator.languages?.length ? navigator.languages : [navigator.language]);
}

const own = (d: Dict, key: string): string | undefined => (Object.prototype.hasOwnProperty.call(d, key) && d[key] ? d[key] : undefined);
const fill = (text: string, vars?: Vars) =>
  vars ? text.replace(/\{(\w+)\}/g, (hole, name: string) => (Object.prototype.hasOwnProperty.call(vars, name) ? String(vars[name]) : hole)) : text;
const makeT = (d: Dict): TFn => (text, vars) => {
  const context = vars?.context;
  return fill((context !== undefined ? own(d, `${text}|${context}`) : undefined) ?? own(d, text) ?? text, vars);
};

// ---------------------------------------------------------------- the store (same shape as themeMode.ts / network.ts)
const listeners = new Set<() => void>();
const arrived = new Map<Lang, Dict>();          // dictionaries that have loaded (or were set by a test)
const loading = new Set<Lang>();
let lang: Lang = initialLang();
let dict: Dict = {};
let current: TFn = makeT(dict);                 // a new function whenever the dictionary changes, so memos that list it re-run
let version = 0;

function changed() {
  current = makeT(dict);
  version++;
  listeners.forEach(l => l());
}

/** The language as an HTML `lang` value (zh is "zh-CN"): the page's <html lang> and the one of the Session report file. */
export const htmlLang = () => HTML_LANG[lang];

function applyHtmlLang() {
  if (typeof document !== "undefined") document.documentElement?.setAttribute("lang", htmlLang());
}

/** Fetch the dictionary of `l` once; it is used only if `l` is still the chosen language when it arrives. A failed load leaves English (and may be retried). */
function load(l: Lang) {
  if (l === "en" || arrived.has(l) || loading.has(l)) return;
  loading.add(l);
  LOADERS[l]().then(m => { arrived.set(l, m.default ?? {}); if (lang === l) { dict = arrived.get(l)!; changed(); } })
    .catch(() => { /* English stays */ })
    .finally(() => loading.delete(l));
}

export function setLang(l: Lang) {
  lang = l;
  try { localStorage.setItem(KEY, l); } catch { /* not remembered */ }
  dict = arrived.get(l) ?? {};
  applyHtmlLang();
  changed();
  load(l);
}

/** Give a language its dictionary now (tests; the real ones arrive by import). */
export function setDictionary(l: Lang, d: Dict) {
  arrived.set(l, d);
  if (lang === l) { dict = d; changed(); }
}

applyHtmlLang();
load(lang);

const subscribe = (l: () => void) => { listeners.add(l); return () => { listeners.delete(l); }; };
const getLang = () => lang;
const getVersion = () => version;

export function useLang(): { lang: Lang; setLang: (l: Lang) => void } {
  return { lang: useSyncExternalStore(subscribe, getLang, getLang), setLang };
}

export function useT(): TFn {
  useSyncExternalStore(subscribe, getVersion, getVersion);
  return current;
}

export const t: TFn = (text, vars) => current(text, vars);

// ---------------------------------------------------------------- the page list
type Kind = { text?: string[]; keyed?: string[]; children?: Record<string, string>; lookup?: string[] };
type Context = Record<string, Record<string, string>>;      // block id -> { English text -> context }
const RULES = specRules as unknown as { root: Record<string, string>; kinds: Record<string, Kind>; context: Context };

const has = <T>(o: Record<string, T> | undefined, key: string): T | undefined => (o && Object.prototype.hasOwnProperty.call(o, key) ? o[key] : undefined);

/** The `t()` options for a spec string drawn inside block `block`: its context when spec-rules.json gives it one ("Species" in the LIPID MAPS cards). */
export function specContext(block: string | null | undefined, text: string): { context: string } | undefined {
  const context = block == null ? undefined : has(has(RULES.context, block), text);
  return context === undefined ? undefined : { context };
}

/**
 * The page list in the current language. spec-rules.json says which keys are display text (and mzdesk/spec_i18n.py reads the same file); each
 * kind of node is copied with those keys translated, everything else is shared as it is, and the input is never changed.
 * A string with a context in spec-rules.json ("context": the block it is in) is looked up as "text|context", and as plain "text" when that is missing.
 * Section titles, cards item labels and the options of a select / radio / multiselect ("lookup" in the rules) stay English: the pages group and
 * compare by them, a run's values are keyed by them, and an option is the value sent to the server; the component that draws one calls t() on it.
 * English, or a dictionary that has not arrived yet: the same object comes back.
 */
export function localizeSpec(spec: Spec): Spec {
  if (!Object.keys(dict).length) return spec;
  const tr = (s: string, ctx?: Record<string, string>) => {
    const context = has(ctx, s);
    return (context ? own(dict, `${s}|${context}`) : undefined) ?? own(dict, s) ?? s;
  };
  const text = (v: unknown, ctx?: Record<string, string>): unknown => (typeof v === "string" ? tr(v, ctx) : Array.isArray(v) ? v.map(x => (typeof x === "string" ? tr(x, ctx) : x)) : v);
  const walk = (node: any, kind: string, inherited?: Record<string, string>): any => {
    if (!node || typeof node !== "object") return node;
    const rule = RULES.kinds[kind], out: any = { ...node };
    const ctx = (kind === "block" ? has(RULES.context, String(node.id)) : undefined) ?? inherited;      // a block's items share its context
    for (const key of rule.text ?? []) if (key in node) out[key] = text(node[key], ctx);
    // column_help: {column name: help}. The English name stays as a second key, for a table whose columns come from a run.
    for (const key of rule.keyed ?? []) {
      const named = node[key];
      if (!named || typeof named !== "object") continue;
      const copy: Record<string, unknown> = {};
      for (const [name, help] of Object.entries(named)) copy[name] = copy[tr(name, ctx)] = text(help, ctx);
      out[key] = copy;
    }
    for (const [key, child] of Object.entries(rule.children ?? {})) if (Array.isArray(node[key])) out[key] = node[key].map((n: unknown) => walk(n, child, ctx));
    return out;
  };
  const out: any = { ...spec };
  for (const [key, kind] of Object.entries(RULES.root)) if (Array.isArray((spec as any)[key])) out[key] = (spec as any)[key].map((n: unknown) => walk(n, kind));
  return out as Spec;
}

/**
 * The blocks of a run (the tables, notes and figures of a workflow step) come from the backend in English. Their titles are page text with a
 * dictionary entry (source/run.json lists them, written by mzdesk/spec_i18n.py from the Python source), so they are translated here; notes, cells
 * and everything else stay as the run wrote them. A figure's label is its title. English, or a dictionary that has not arrived: the same array.
 */
export function localizeRunBlocks<B extends { type: string; title?: string | null; label?: string | null }>(blocks: B[]): B[] {
  if (!Object.keys(dict).length) return blocks;
  const tr = (s: string) => own(dict, s) ?? s;
  return blocks.map(b => {
    const title = b.title ? tr(b.title) : b.title;
    const label = b.type === "figure" && b.label ? tr(b.label) : b.label;
    return title === b.title && label === b.label ? b : { ...b, title, label };
  });
}
