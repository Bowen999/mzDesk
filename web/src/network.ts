import { useSyncExternalStore } from "react";
import type { Field } from "./types";

/**
 * "Internet search" switch (sidebar footer).
 *   On  = online lookups run without asking.
 *   Off = a run that needs the internet first shows a confirmation dialog; the run continues only if the user agrees.
 * The choice is remembered in localStorage. Both switch shapes (sidebar and rail) and the pages read it through
 * useInternet(), so they always agree. The switch is a preference, not access control: the API stays open.
 */
const KEY = "mzdesk.internetSearch";
let memory: boolean | null = null;
const listeners = new Set<() => void>();

function read(): boolean {
  if (memory === null) {
    try { memory = localStorage.getItem(KEY) !== "0"; } catch { memory = true; }
  }
  return memory;
}

export function isInternetOn(): boolean {
  return read();
}

export function setInternet(on: boolean) {
  memory = on;
  try { localStorage.setItem(KEY, on ? "1" : "0"); } catch { /* memory only */ }
  listeners.forEach(l => l());
}

const subscribe = (listener: () => void) => {
  listeners.add(listener);
  return () => { listeners.delete(listener); };
};

export function useInternet(): boolean {
  return useSyncExternalStore(subscribe, read, read);
}

// The 118 element symbols (mzdesk/api/common.py is_formula() checks the same thing with RDKit).
const ELEMENTS = new Set(("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "
  + "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir "
  + "Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og").split(" "));

/** A formula with only real elements (C8H10N4O2) or an InChI is parsed locally: no prompt for those in Name / ID mode. */
function parsedLocally(text: string): boolean {
  if (text.startsWith("InChI=")) return true;
  if (!/^(?:[A-Z][a-z]?\d*)+$/.test(text)) return false;
  return [...text.matchAll(/[A-Z][a-z]?/g)].every(m => ELEMENTS.has(m[0]));
}

// A candidate list (spec `network_names`) holds one SMILES, InChI or name per line. The backend reads a line as a SMILES when RDKit can and
// otherwise looks the name up online; this is the same question asked without RDKit: only SMILES characters in SMILES order, an atom in
// it, balanced brackets, every ring-closure digit used an even number of times ("F6P", "C1" and "PCB-187" are names, not SMILES), and a
// ring-closure digit when it has aromatic (lower-case) atoms ("Co" is cobalt, not carbon + aromatic oxygen). A line with a space is a name
// (RDKit would read "PC 33:1" as the SMILES "PC" plus a comment; the backend does not).
const SMILES_TOKENS = /^(?:\[[^\]\s]+\]|Br|Cl|[BCNOPSFI]|[bcnops]|%\d\d|\d|[-=#$:/\\().+~*@])+$/;
function wellFormed(line: string): boolean {
  const body = line.replace(/\[[^\]]*\]/g, "[]");                   // digits inside [..] are charges, isotopes, H counts - not ring closures
  let depth = 0;
  for (const ch of body) {
    if (ch === "(") depth++;
    else if (ch === ")" && --depth < 0) return false;
  }
  if (depth !== 0) return false;
  const uses = new Map<string, number>();
  for (const m of body.matchAll(/%\d\d|\d/g)) uses.set(m[0], (uses.get(m[0]) ?? 0) + 1);
  return [...uses.values()].every(n => n % 2 === 0);
}
function looksLikeSmiles(line: string): boolean {
  if (line.startsWith("InChI=")) return true;                       // parsed locally
  if (!SMILES_TOKENS.test(line) || !/[A-Za-z[]/.test(line) || !wellFormed(line)) return false;
  return !/[bcnops]/.test(line.replace(/\[[^\]]*\]/g, "")) || /[\d%]/.test(line);
}

/** The lines of a candidate list that would be looked up by name (blank lines and # comments are ignored, as the backend does). */
export function candidateNames(text: string): string[] {
  return text.split(/\r?\n/).map(l => l.trim()).filter(l => l && !l.startsWith("#") && !looksLikeSmiles(l));
}

/**
 * What this run would send online, as sentences for the confirmation dialog. Empty = the run stays local.
 * The rules come from the spec: module.network (the page is online), and the field flags network_sources /
 * network_modes / network / network_options / network_names that spec.py sets on the inputs that can use the internet.
 */
export function networkUses(fields: Field[], values: Record<string, any>, page?: string): string[] {
  const uses: string[] = page ? [page] : [];
  for (const f of fields) {
    if (!f.id) continue;
    const v = values[f.id];
    if (f.type === "spectrum" && f.network_sources?.length) {
      const source = String(v?.source ?? "");
      if (f.network_sources.includes(source) && String(v?.value ?? "").trim())
        uses.push("The spectrum is fetched from GNPS / MassIVE (USI).");
    } else if (f.type === "structure" && f.network_modes?.length) {
      const text = String(v?.value ?? "").trim();
      // A page-level sentence already covers the page's own lookup; a formula / InChI is parsed locally.
      if (!page && text && f.network_modes.includes(String(v?.mode ?? "")) && !parsedLocally(text))
        uses.push("The compound is looked up in PubChem and UniChem.");
    } else if (f.type === "textarea" && f.network_names) {
      const names = candidateNames(String(v ?? ""));
      if (names.length) uses.push(`${names.length === 1 ? "A candidate name is" : "Candidate names are"} looked up in PubChem and UniChem.`);
    } else if (f.type === "checkbox" && f.network && v === true) {
      uses.push(`“${f.label ?? "This option"}”: models from Zenodo.`);
    } else if (f.type === "multiselect" && f.network_options?.length && Array.isArray(v)) {
      const picked = v.filter((x: string) => f.network_options!.includes(x));
      if (picked.length) uses.push(`${picked.join(" and ")}: ${picked.length === 1 ? "model" : "models"} from Zenodo.`);
    }
  }
  return [...new Set(uses)];
}
