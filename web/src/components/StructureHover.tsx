import { createContext, useContext, useEffect, useMemo, useState, type ReactNode } from "react";
import { BASE } from "../api";
import type { AnnotationResult } from "../types";
import { ThemedSvg } from "./Figure";
import { useHelpPopover } from "./Help";

/**
 * Hover previews of 2D structures. Anything that shows a SMILES or a compound name can wrap itself in <StructureHover smiles>: hover or
 * keyboard focus draws the structure in the same popover as the "?" help, a click pins it. The drawing is made by the server
 * (GET /api/depict, RDKit) the first time it is needed and kept in this tab, so the result itself stays small.
 *
 * ESI MS/MS Annotation wraps its page in <StructureProvider>: the names and SMILES of the candidates (structureMap) are then recognised
 * in the answer text, the ranked table and the step tables, so a name that appears anywhere on the page shows its structure.
 */
const drawn = new Map<string, Promise<string | null>>();

export function fetchDepiction(smiles: string): Promise<string | null> {
  let p = drawn.get(smiles);
  if (!p) {
    p = fetch(`${BASE}api/depict?smiles=${encodeURIComponent(smiles)}`)
      .then(r => (r.ok ? r.text() : null))
      .catch(() => { drawn.delete(smiles); return null; });      // a dropped request may be tried again
    drawn.set(smiles, p);
  }
  return p;
}

/** The popover content: reserved space first, so the popover is placed once and does not jump when the drawing arrives. */
function StructurePreview({ smiles }: { smiles: string }) {
  const [svg, setSvg] = useState<string | null | undefined>(undefined);
  useEffect(() => { let on = true; setSvg(undefined); fetchDepiction(smiles).then(s => { if (on) setSvg(s); }); return () => { on = false; }; }, [smiles]);
  return (
    <div className="struct-preview">
      {svg ? <ThemedSvg svg={svg} className="pop-struct" alt="2D structure" /> : <span className="struct-msg">{svg === undefined ? "Drawing…" : "No structure to draw"}</span>}
    </div>);
}

/** The text (a name, a SMILES) with its structure on hover. */
export function StructureHover({ smiles, children }: { smiles: string; children: ReactNode }) {
  const { triggerProps, popover } = useHelpPopover([], <StructurePreview smiles={smiles} />);
  return <><button type="button" className="cell-depict" {...triggerProps}>{children}</button>{popover}</>;
}

// ---------------------------------------------------------------- page-wide: names and SMILES of one result
type Structures = { names: Map<string, string> };
const StructureContext = createContext<Structures | null>(null);
export const useStructures = () => useContext(StructureContext);

/** name -> SMILES, and SMILES -> itself, for everything the result names: candidates, isomer-comparison candidates, the library hit. */
export function structureMap(result: AnnotationResult | null): Map<string, string> {
  const m = new Map<string, string>();
  const add = (name?: string | null, smiles?: string | null) => {
    const s = smiles?.trim();
    if (!s) return;
    if (!m.has(s)) m.set(s, s);
    const n = name?.trim();
    if (n && !m.has(n)) m.set(n, s);
  };
  for (const c of result?.candidates ?? []) add(c.name, c.smiles);
  for (const c of result?.isomer?.candidates ?? []) add(c.name, c.smiles);
  const lib = result?.visualization?.library;
  if (lib) add(lib.name, lib.smiles);
  return m;
}

export function StructureProvider({ result, children }: { result: AnnotationResult | null; children: ReactNode }) {
  const value = useMemo(() => ({ names: structureMap(result) }), [result]);
  return <StructureContext.Provider value={value}>{children}</StructureContext.Provider>;
}

const escapeRe = (s: string) => s.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");

/** Plain text in which every known compound name or SMILES (whole words only) shows its structure on hover. Outside a provider: the text as is. */
export function StructText({ text }: { text: string }) {
  const ctx = useStructures();
  const re = useMemo(() => {
    const keys = [...(ctx?.names.keys() ?? [])].filter(k => k.length >= 3).sort((a, b) => b.length - a.length);
    return keys.length ? new RegExp(`(?<![\\p{L}\\p{N}])(${keys.map(escapeRe).join("|")})(?![\\p{L}\\p{N}])`, "gu") : null;
  }, [ctx]);
  if (!ctx || !re || !text) return <>{text}</>;
  const out: ReactNode[] = [];
  let last = 0, m: RegExpExecArray | null;
  re.lastIndex = 0;
  while ((m = re.exec(text))) {
    if (m.index > last) out.push(text.slice(last, m.index));
    out.push(<StructureHover key={m.index} smiles={ctx.names.get(m[1])!}>{m[1]}</StructureHover>);
    last = m.index + m[1].length;
  }
  if (!out.length) return <>{text}</>;
  if (last < text.length) out.push(text.slice(last));
  return <>{out}</>;
}

/** Which structure, if any, a plain table cell stands for: a SMILES column shows its own text; a known name shows its candidate's structure;
 *  a name column falls back to the SMILES of the same row. Only inside a provider. */
export function structureForCell(ctx: Structures | null, text: string, column?: string, rowSmiles?: string): string | null {
  const t = text.trim();
  if (!ctx || !t || t === "—") return null;
  if (column && /smiles/i.test(column)) return t;
  const known = ctx.names.get(t);
  if (known) return known;
  if (column && rowSmiles && /^(name|compound|candidate)/i.test(column)) return rowSmiles;
  return null;
}
