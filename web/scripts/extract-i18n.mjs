// Lists every interface string: the literal first argument of t("...") / t('...') in web/src/**/*.{ts,tsx}.
// Writes web/src/i18n/source/ui.json (sorted, unique). The translators work from that list and from source/spec.json.
//   node scripts/extract-i18n.mjs          write the file
//   node scripts/extract-i18n.mjs --check  write nothing; exit 1 when the file is out of date
// t("System", { context: "theme" }) is the key "System|theme" (the same English word with another meaning; "context" must be the first key of the object).
// A t() call whose first argument is not one string literal (a variable, a template literal, "a" + b) cannot be found here; the ones
// that are meant (a section title or a cards label from the spec, listed in source/spec.json) are only reported, the others are worth a look.
import { readdirSync, readFileSync, writeFileSync } from "node:fs";
import { dirname, join, relative } from "node:path";
import { fileURLToPath } from "node:url";

const web = dirname(dirname(fileURLToPath(import.meta.url)));
const src = join(web, "src");
const out = join(src, "i18n", "source", "ui.json");

function* files(dir) {
  for (const e of readdirSync(dir, { withFileTypes: true })) {
    const p = join(dir, e.name);
    if (e.isDirectory()) { if (p !== join(src, "i18n")) yield* files(p); }
    else if (/\.(ts|tsx)$/.test(e.name) && !e.name.endsWith(".d.ts")) yield p;
  }
}

/** Comments out of the way, so an example t("...") in a doc comment is not a key. Blanked, not removed, so line numbers stay right.
 *  A "//" after a colon is a URL, not a comment. */
function withoutComments(code) {
  const blank = m => m.replace(/[^\n]/g, " ");
  return code.replace(/\/\*[\s\S]*?\*\//g, blank).replace(/(^|[^:"'\w/])\/\/[^\n]*/g, (m, lead) => lead + blank(m.slice(lead.length)));
}

/** The value of a quoted JavaScript string literal (\" \' \\ \n \t \r \0 \xHH \uHHHH \u{H...}; any other escape stands for the character). */
function unquote(literal) {
  return literal.slice(1, -1).replace(/\\(u\{[0-9a-fA-F]+\}|u[0-9a-fA-F]{4}|x[0-9a-fA-F]{2}|[\s\S])/g, (_, e) => {
    if (e[0] === "u" && e.length > 1) return String.fromCodePoint(parseInt(e.replace(/[u{}]/g, ""), 16));
    if (e[0] === "x" && e.length > 1) return String.fromCharCode(parseInt(e.slice(1), 16));
    return { n: "\n", t: "\t", r: "\r", "0": "\0" }[e] ?? e;
  });
}

const CALL = /(?<![\w$.])t\(\s*/g;
const LITERAL = /^("(?:[^"\\\n]|\\.)*"|'(?:[^'\\\n]|\\.)*')\s*([,)])?/;

const CONTEXT = /^\s*\{\s*context\s*:\s*("(?:[^"\\\n]|\\.)*"|'(?:[^'\\\n]|\\.)*')\s*[,}]/;

const keys = new Set(), odd = [], misplaced = [];
for (const file of files(src)) {
  const code = withoutComments(readFileSync(file, "utf8"));
  const lineOf = at => code.slice(0, at).split("\n").length;
  for (const m of code.matchAll(CALL)) {
    const rest = code.slice(m.index + m[0].length);
    const lit = LITERAL.exec(rest);
    const where = `${relative(web, file)}:${lineOf(m.index)}`;
    if (!lit) { odd.push(`${where}  t(${rest.slice(0, 28).split("\n")[0]}...`); continue; }
    if (!lit[2]) { odd.push(`${where}  t(${lit[1]} + ...) is not one literal`); continue; }
    let key = unquote(lit[1]);
    if (lit[2] === ",") {
      const after = rest.slice(lit[0].length), context = CONTEXT.exec(after);
      if (context) key += "|" + unquote(context[1]);
      else if (/^\s*\{[^}]*\bcontext\s*:/.test(after)) misplaced.push(`${where}  t(${lit[1]}, ...): context must be the first key, with a string literal`);
    }
    keys.add(key);
  }
}
if (misplaced.length) { console.error(misplaced.join("\n")); process.exit(1); }

const list = [...keys].sort();
const text = JSON.stringify(list, null, 2) + "\n";
if (process.argv.includes("--check")) {
  let old = "";
  try { old = readFileSync(out, "utf8"); } catch { /* missing = stale */ }
  if (old !== text) { console.error(`${relative(web, out)} is out of date: run node scripts/extract-i18n.mjs`); process.exit(1); }
  console.log(`${relative(web, out)} is up to date (${list.length} strings)`);
} else {
  writeFileSync(out, text);
  console.log(`${list.length} strings -> ${relative(web, out)}`);
  if (odd.length) console.log(`${odd.length} t() calls without a literal key (fine for spec.json strings):\n  ${odd.join("\n  ")}`);
}
