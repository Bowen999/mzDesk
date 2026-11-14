import { t } from "../i18n";
import type { AnnotationCandidate, AnnotationResult, CandidateCcs } from "../types";

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
/** "A", "A and B", "A, B and 3 more". */
const listOf = (items: string[]) => (items.length <= 1 ? items.join("") : items.length === 2 ? t("{a} and {b}", { a: items[0], b: items[1] })
  : t("{a}, {b} and {n} more", { a: items[0], b: items[1], n: items.length - 2 }));

/** The words the backend uses for the isomer check (answer.py: verdict single | indistinguishable | insufficient | not run; role top | separated | not separated).
 *  Shown through the dictionary; any other word is shown as received. */
export function isomerWord(word: string): string {
  switch (word) {
    case "single": return t("single");
    case "indistinguishable": return t("indistinguishable");
    case "insufficient": return t("insufficient");
    case "not run": return t("not run");
    case "top": return t("top");
    case "separated": return t("separated");
    case "not separated": return t("not separated");
    default: return word;
  }
}

/** The CCS verdicts of the backend (ccs_step.py); any other word is shown as received. */
export function ccsVerdict(word: string): string {
  switch (word) {
    case "Consistent": return t("Consistent");
    case "Inconsistent": return t("Inconsistent");
    case "No reference": return t("No reference");
    case "No structure": return t("No structure");
    case "Not checked": return t("Not checked");
    default: return word;
  }
}

/** Where the top-ranked structure comes from, in a few words. */
function sourceText(level: string, top: AnnotationCandidate | undefined, layer: string): string {
  if (level === "2a") return t("Library match, gate passed");
  if (layer === "library") return top?.evidence?.library?.gate_passed ? t("Library hit without a structure") : t("Library hit, gate not passed");
  if (layer === "msfinder") return t("MS-FINDER candidate");
  if (layer === "frigid") return t("FRIGID hypothesis");
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
    headline = [first.formula || first.name, first.adduct].filter(Boolean).join(" ") + (rows.length > 1 ? " " + t("+{n} more", { n: rows.length - 1 }) : "");
    score = scoreOf(first); source = t("Formula only, no structure proposed");
  } else if (a.type === "mass_only" && Number.isFinite(result.context?.precursor_mz)) {
    headline = `m/z ${result.context.precursor_mz.toFixed(4)}`; source = t("Exact mass only");
  }

  // the isomer part: what the spectrum can and cannot separate
  const others = (a.ids ?? []).filter(id => id !== a.top_id).map(label);
  const sentence = challengers.map(ch => t("{name} explains peaks the top-ranked structure cannot (m/z {peaks})",
    { name: ch.name || label(ch.id), peaks: (ch.peaks ?? []).map(p => p.mz.toFixed(4)).join(", ") }));
  const names = listOf(challengers.map(ch => ch.name || label(ch.id)));
  const challenger = challengers.length === 1 ? t("{names} explains peaks it cannot", { names }) : challengers.length ? t("{names} explain peaks it cannot", { names }) : "";
  const parts: string[] = [];
  let warn = false;
  if (a.type === "group" && others.length) { parts.push(t("Indistinguishable from {others} by this spectrum", { others: listOf(others) })); warn = true; }
  else {
    if (source) parts.push(source);
    if (a.type === "single" || a.type === "group") {
      if (conflicting) { parts.push(t("Conflicting fragment evidence")); warn = true; }
      else if (verdict === "single") parts.push(t("Told apart from the other isomers considered"));
      else if (verdict === "not run") parts.push(t("Isomer check not run"));
      else if (verdict === "indistinguishable") { parts.push(t("Isomers not separated")); warn = true; }
      else if (verdict) parts.push(t("Isomer check inconclusive"));
    }
  }
  if (challenger) { parts.push(challenger); warn = true; }

  const facts: [string, string][] = [];
  if (a.layer_label) facts.push([t("Layer"), a.layer_label]);
  if (a.level_label) facts.push([t("Level"), a.level_label]);
  if (score !== null) facts.push([t("Final score"), score.toFixed(3)]);
  if (verdict && verdict !== "not run") facts.push([t("Isomer check"), conflicting ? t("conflicting evidence") : isomerWord(verdict)]);
  if (result.context?.gate_passed) facts.push([t("Library gate"), t("passed")]);
  if (top?.ccs) facts.push(["CCS", ccsText(top.ccs)]);
  return { headline, level, score, secondary: parts.join(" · "), warn, facts, challengers: sentence };
}

// ---------------------------------------------------------------- CCS check (supporting evidence; mzdesk/annotation/ccs_step.py)
export const CCS_TONE: Record<string, string> = { Consistent: "good", Inconsistent: "warn" };
const signedPct = (v: number) => `${v < 0 ? "−" : "+"}${Math.abs(v).toFixed(2)} %`;

/** The numbers behind a CCS verdict: "+0.64 % (±1.0 %) · measured 139.6 Å² · DTCCS N2 · MetCCS"; the engine's note when there is no verdict; "" when there is neither. */
export function ccsNumbers(x: CandidateCcs): string {
  if (typeof x.delta_pct === "number" && typeof x.reference === "number" && typeof x.tolerance_pct === "number") {
    return [`${signedPct(x.delta_pct)} (±${x.tolerance_pct.toFixed(1)} %)`, `${x.source ?? t("reference")} ${x.reference.toFixed(1)} Å²`, x.ref_label, x.dataset,
      x.in_domain === false ? t("outside its training domain") : ""].filter(Boolean).join(" · ");
  }
  return x.note ?? "";
}

/** One line for the answer's Details: "Inconsistent · +9.80 % (±1.0 %) · measured 140.5 Å² · DTCCS N2 · tag · moved from rank 1". */
export function ccsText(x: CandidateCcs): string {
  return [ccsVerdict(x.verdict), ccsNumbers(x), x.moved_from ? t("moved from rank {n}", { n: x.moved_from }) : ""].filter(Boolean).join(" · ");
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
  return (iso.role != null ? isomerWord(iso.role) : t("not compared")) + (iso.challenger ? " · " + t("challenger") : "");
}

/** The columns of a ranked row that are folded into its detail line: [label, value]. */
export function candidateDetails(c: AnnotationCandidate): [string, string][] {
  return [["SMILES", c.smiles || "—"], [t("Formula"), c.formula || "—"], [t("Adduct"), c.adduct || "—"], [t("Layer"), c.layer_label || "—"],
    [t("Detectability (omics / ionizable / adduct)"), detectText(c)], [t("Isomer role"), isomerText(c)]];
}
