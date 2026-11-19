import { useSyncExternalStore } from "react";

/** Light / dark theme. The choice ("system" follows the OS) is remembered per browser; the resolved theme sits on <html data-theme>.
 *  CSS reads the attribute (and prefers-color-scheme when the choice is "system"); figures read useTheme().dark to recolour. */
export type ThemePref = "system" | "light" | "dark";
const KEY = "mzdesk.theme";
const listeners = new Set<() => void>();
const mq = typeof window !== "undefined" && window.matchMedia ? window.matchMedia("(prefers-color-scheme: dark)") : null;

let pref: ThemePref = (() => {
  try { const v = localStorage.getItem(KEY); return v === "light" || v === "dark" ? v : "system"; } catch { return "system"; }
})();

export const resolvedTheme = (): "light" | "dark" => (pref === "system" ? (mq?.matches ? "dark" : "light") : pref);

function apply() {
  const root = document.documentElement;
  if (pref === "system") root.removeAttribute("data-theme"); else root.setAttribute("data-theme", pref);
  listeners.forEach(l => l());
}
apply();
mq?.addEventListener?.("change", apply);

export function setThemePref(p: ThemePref) {
  pref = p;
  try { if (p === "system") localStorage.removeItem(KEY); else localStorage.setItem(KEY, p); } catch { /* not remembered */ }
  apply();
}

const subscribe = (l: () => void) => { listeners.add(l); return () => { listeners.delete(l); }; };
const snapshot = () => `${pref}|${resolvedTheme()}`;

export function useTheme() {
  const [p, r] = useSyncExternalStore(subscribe, snapshot, snapshot).split("|") as [ThemePref, "light" | "dark"];
  return { pref: p, dark: r === "dark", setPref: setThemePref };
}
