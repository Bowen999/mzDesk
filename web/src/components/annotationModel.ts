import type { AnnotationCandidate, AnnotationResult } from "../types";

/** Pure helpers of the ESI MS/MS Annotation result view (answer card, candidate table): everything is derived from the result's own fields,
 *  nothing is scored or ranked here. */

export const LEVEL_TONE: Record<string, string> = { "2a": "good", "3": "jp", H: "warn", "4": "neutral", "5": "neutral" };
/** "Level 2a: probable structure (library spectrum match)" -> "Level 2a". */
export const levelName = (label: string) => label.split(":")[0];

export type AnswerSummary = {
  headline: string | null;            // the compound (or formula / exact mass) the answer names; null: only the statement can be shown
  level: { name: string; title: string; tone: string } | null;
  score: number | null;               // final score of the answered candidate
  secondary: string;                  // one short line: where it comes from, and what the spectrum cannot separate; "" when there is nothing to add
  warn: boolean;                      // the secondary line carries a caveat (isomers not separated, conflicting evidence, a challenger)
  facts: [string, string][];          // everything the compact card does not say, for the "Details" disclosure
  challengers: string[];              // one sentence per isomer that explains peaks the top-1 cannot (with their m/z)
};

const candidateLabel = (c: AnnotationCandidate) => c.name || c.smiles || c.id;
const scoreOf = (c: AnnotationCandidate) => (typeof c.final_score === "number" && Number.isFinite(c.final_score) ? c.final_score : null);
const listOf = (items: string[], shown = 2) => (items.length <= shown ? items.join(" and ") : `${items.slice(0, shown).join(", ")} and ${items.length - shown} more`);

/** Where the top-ranked structure comes from, in a few words. */
function sourceText(level: string, top: AnnotationCandidate | undefined, layer: string): string {
  if (level === "2a") return "Library match, gate passed";
  if (layer === "library") return top?.evidence?.library?.gate_passed ? "Library hit without a structure" : "Library hit, gate not passed";
  if (layer === "msfinder") return "MS-FINDER candidate";
  if (layer === "frigid") return "FRIGID hypothesis";
  return "";
}

/** The compact answer: "Caffeine · Level 2a · 0.80" and one short line under it. The full statement stays available as the first fact. */
export function answerSummary(result: AnnotationResult): AnswerSummary {
  const a = result.answer;
  const candidates = result.candidates ?? [];
  const byId = new Map(candidates.map(c => [c.id, c]));
  const label = (id: string) => { const c = byId.get(id); return c ? candidateLabel(c) : id; };
  const iso = result.isomer;
  const challengers = iso?.challengers ?? [];
  const verdict = a.isomer_verdict;
  const conflicting = verdict === "single" && challengers.length > 0;
  const top = a.top_id ? byId.get(a.top_id) : undefined;
  const rows = (a.ids ?? []).map(id => byId.get(id)).filter((c): c is AnnotationCandidate => !!c);
  const level = a.level ? { name: a.level_label ? levelName(a.level_label) : a.level, title: a.level_label ?? "", tone: LEVEL_TONE[a.level] ?? "neutral" } : null;

  let headline: string | null = null, score: number | null = null, source = "";
  if ((a.type === "single" || a.type === "group") && top) {
    headline = candidateLabel(top); score = scoreOf(top); source = sourceText(a.level, top, a.layer);
  } else if (a.type === "formula_only" && rows.length) {
    const first = rows[0];
    headline = [first.formula || first.name, first.adduct].filter(Boolean).join(" ") + (rows.length > 1 ? ` +${rows.length - 1} more` : "");
    score = scoreOf(first); source = "Formula only, no structure proposed";
  } else if (a.type === "mass_only" && Number.isFinite(result.context?.precursor_mz)) {
    headline = `m/z ${result.context.precursor_mz.toFixed(4)}`; source = "Exact mass only";
  }

  // the isomer part: what the spectrum can and cannot separate
  const others = (a.ids ?? []).filter(id => id !== a.top_id).map(label);
  const sentence = challengers.map(ch => `${ch.name || label(ch.id)} explains peaks the top-ranked structure cannot (m/z ${(ch.peaks ?? []).map(p => p.mz.toFixed(4)).join(", ")})`);
  const challenger = challengers.length ? `${listOf(challengers.map(ch => ch.name || label(ch.id)))} ${challengers.length === 1 ? "explains" : "explain"} peaks it cannot` : "";
  const parts: string[] = [];
  let warn = false;
  if (a.type === "group" && others.length) { parts.push(`Indistinguishable from ${listOf(others)} by this spectrum`); warn = true; }
  else {
    if (source) parts.push(source);
    if (a.type === "single" || a.type === "group") {
      if (conflicting) { parts.push("Conflicting fragment evidence"); warn = true; }
      else if (verdict === "single") parts.push("Told apart from the other isomers considered");
      else if (verdict === "not run") parts.push("Isomer check not run");
      else if (verdict === "indistinguishable") { parts.push("Isomers not separated"); warn = true; }
      else if (verdict) parts.push("Isomer check inconclusive");
    }
  }
  if (challenger) { parts.push(challenger); warn = true; }

  const facts: [string, string][] = [];
  if (a.layer_label) facts.push(["Layer", a.layer_label]);
  if (a.level_label) facts.push(["Level", a.level_label]);
  if (score !== null) facts.push(["Final score", score.toFixed(3)]);
  if (verdict && verdict !== "not run") facts.push(["Isomer check", conflicting ? "conflicting evidence" : verdict]);
  if (result.context?.gate_passed) facts.push(["Library gate", "passed"]);
  return { headline, level, score, secondary: parts.join(" · "), warn, facts, challengers: sentence };
}

// ---------------------------------------------------------------- candidate table
export const detectText = (c: AnnotationCandidate) => {
  const d = c.detect;
  return d ? `${d.omics ?? "—"} / ${d.ionizable ?? "—"} / ${d.adduct ?? "—"}` : "—";
};

/** The isomer role of a candidate as text ("—" when it was not compared). */
export function isomerText(c: AnnotationCandidate): string {
  const iso = c.isomer;
  if (!iso || (!iso.role && !iso.challenger)) return "—";
  return (iso.role ?? "not compared") + (iso.challenger ? " · challenger" : "");
}

/** The columns of a ranked row that are folded into its detail line: [label, value]. */
export function candidateDetails(c: AnnotationCandidate): [string, string][] {
  return [["SMILES", c.smiles || "—"], ["Formula", c.formula || "—"], ["Adduct", c.adduct || "—"], ["Layer", c.layer_label || "—"],
    ["Detectability (omics / ionizable / adduct)", detectText(c)], ["Isomer role", isomerText(c)]];
}
