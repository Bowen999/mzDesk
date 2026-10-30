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

/**
 * What this run would send online, as sentences for the confirmation dialog. Empty = the run stays local.
 * The rules come from the spec: module.network (the page is online), and the field flags network_sources /
 * network_modes / network / network_options that spec.py sets on the inputs that can use the internet.
 */
export function networkUses(fields: Field[], values: Record<string, any>, page?: string): string[] {
  const uses: string[] = page ? [page] : [];
  for (const f of fields) {
    if (!f.id) continue;
    const v = values[f.id];
    if (f.type === "spectrum" && f.network_sources?.length) {
      const source = String(v?.source ?? "");
      if (f.network_sources.includes(source) && String(v?.value ?? "").trim())
        uses.push("The spectrum is resolved through the public GNPS / MassIVE USI service.");
    } else if (f.type === "structure" && f.network_modes?.length) {
      const text = String(v?.value ?? "").trim();
      // A page-level sentence already covers the page's own lookup; a formula / InChI is parsed locally.
      if (!page && text && f.network_modes.includes(String(v?.mode ?? "")) && !parsedLocally(text))
        uses.push("The compound is looked up by name or database ID in PubChem and UniChem.");
    } else if (f.type === "checkbox" && f.network && v === true) {
      uses.push(`“${f.label ?? "This option"}” downloads its model from Zenodo on first use.`);
    } else if (f.type === "multiselect" && f.network_options?.length && Array.isArray(v)) {
      const picked = v.filter((x: string) => f.network_options!.includes(x));
      if (picked.length) uses.push(`${picked.join(" and ")} ${picked.length === 1 ? "downloads its model" : "download their models"} from Zenodo on first use.`);
    }
  }
  return [...new Set(uses)];
}
