import type { Spec } from "../types";

/** One page of the Docs. Every page loads its text from GET /api/docs/<id> (mzdesk/docs/<id>.md); `kind` says what it is:
 *  "page" = a general page (Getting Started, Concepts, Reference), "module" / "workflow" = the notes of that app page. */
export type DocItem = { id: string; title: string; kind: "page" | "module" | "workflow" };
export type DocFolder = { id: string; title: string; items: DocItem[] };
export type DocSection = { title: string; entries: (DocItem | DocFolder)[] };

export const isFolder = (e: DocItem | DocFolder): e is DocFolder => "items" in e;

const page = (id: string, title: string): DocItem => ({ id, title, kind: "page" });

/** The general pages in reading order. The backend serves exactly these ids (GENERAL_PAGES in mzdesk/api/app.py; web/tests keeps the two equal). */
export const GENERAL_PAGES: readonly DocItem[] = [
  page("overview", "Overview"),
  page("quick-start", "Quick Start"),
  page("local-version", "Local Version"),
  page("workflows-and-modules", "Workflows and Modules"),
  page("inputs", "Inputs and Formats"),
  page("reading-results", "Reading Results"),
  page("session-report", "Session Report"),
  page("privacy", "Privacy and Internet Search"),
  page("cli", "Command Line and Settings"),
  page("faq", "FAQ"),
  page("changelog", "Changelog"),
  page("licenses", "Licenses and Based On"),
  page("feedback", "Feedback"),
  page("contact", "Contact Us"),
];
const general = (...ids: string[]) => ids.map(id => GENERAL_PAGES.find(p => p.id === id)!);

/** The Docs table of contents. Modules and Workflows follow the app sidebar (spec order); the rest is fixed. */
export function docsOutline(spec: Spec): DocSection[] {
  const mods = spec.modules.filter(m => m.group === "modules");
  const modItem = (m: Spec["modules"][number]): DocItem => ({ id: m.id, title: m.title, kind: "module" });
  const modules: (DocItem | DocFolder)[] = spec.subgroups
    .map(s => ({ id: s.id, title: s.title, items: mods.filter(m => m.subgroup === s.id).map(modItem) }))
    .filter(f => f.items.length);
  for (const m of mods) if (!m.subgroup || !spec.subgroups.some(s => s.id === m.subgroup)) modules.push(modItem(m));

  return [
    { title: "Getting Started", entries: general("overview", "quick-start", "local-version") },
    { title: "Concepts", entries: general("workflows-and-modules", "inputs", "reading-results", "session-report", "privacy") },
    { title: "Workflows", entries: spec.pipelines.map(p => ({ id: p.id, title: p.title, kind: "workflow" as const })) },
    { title: "Modules", entries: modules },
    { title: "Reference", entries: general("cli", "faq", "changelog", "licenses", "feedback", "contact") },
  ];
}

/** Every page in reading order (for Previous / Next). */
export function flatPages(outline: DocSection[]): DocItem[] {
  return outline.flatMap(s => s.entries.flatMap(e => (isFolder(e) ? e.items : [e])));
}

/** Heading text → anchor id ("Results" → "results"). */
export function slug(text: string): string {
  return text.toLowerCase().replace(/<[^>]+>|\*\*|`/g, "").replace(/[^a-z0-9]+/g, "-").replace(/^-|-$/g, "") || "section";
}

/** The opening line of a fenced code block: ``` or ~~~ with an optional language ("```bash"); null for any other line. */
export function fenceOpen(line: string): { marker: string; lang: string; indent: number } | null {
  const m = /^(\s*)(```|~~~)\s*([\w+.-]*)\s*$/.exec(line);
  return m ? { marker: m[2], lang: m[3].toLowerCase(), indent: m[1].length } : null;
}
/** The line that ends a fence opened with `marker`: the bare marker. */
export const fenceClose = (line: string, marker: string) => line.trim() === marker;

type Heading = { id: string; text: string; level: 1 | 2 | 3 };

/** The #, ## and ### headings of a markdown text, ids made unique the same way the Markdown component makes them.
 *  Lines inside fenced code blocks are not headings. */
function scanHeadings(text: string): Heading[] {
  const seen = new Map<string, number>();
  const out: Heading[] = [];
  let fence = "";
  for (const line of text.split("\n")) {
    if (fence) { if (fenceClose(line, fence)) fence = ""; continue; }
    const f = fenceOpen(line);
    if (f) { fence = f.marker; continue; }
    const h = /^(#{1,3})\s+(.*)$/.exec(line.trimEnd());
    if (!h) continue;
    out.push({ id: uniqueSlug(h[2], seen), text: h[2].replace(/\*+|`/g, "").replace(/<[^>]+>/g, ""), level: h[1].length as 1 | 2 | 3 });
  }
  return out;
}

/** The ## and ### headings of a markdown text, for "On this page". */
export function headingsOf(text: string): { id: string; text: string; level: 2 | 3 }[] {
  return scanHeadings(text).filter((h): h is Heading & { level: 2 | 3 } => h.level > 1);
}

/** Does the text start its own page title (a # heading outside code blocks)? */
export const hasTitle = (text: string) => scanHeadings(text).some(h => h.level === 1);

export function uniqueSlug(text: string, seen: Map<string, number>): string {
  const base = slug(text), n = seen.get(base) ?? 0;
  seen.set(base, n + 1);
  return n ? `${base}-${n}` : base;
}
