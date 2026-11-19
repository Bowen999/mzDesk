/** Recolour figures drawn for a white page (RDKit SVG, Plotly JSON) for the dark theme: every "#RRGGBB" is remapped by lightness.
 *  Pure black "#000000" (the fragment view's whole structure and spectrum) becomes pure white; colours listed in `keep` (its fragment colours) stay as drawn.
 *  near-black ink -> light ink; saturated mid colours (atoms, data series) -> lifted; pale tints (highlights) -> deep muted tints;
 *  near-white (grid, hover boxes; Plotly only) -> dark surfaces. Colours given as rgba() are left alone. */
const BG = [24, 27, 33];          // --panel in the dark theme
const INK = "#E3E6EC";

function hsl(r: number, g: number, b: number) {
  const mx = Math.max(r, g, b) / 255, mn = Math.min(r, g, b) / 255, l = (mx + mn) / 2;
  const s = mx === mn ? 0 : l > 0.5 ? (mx - mn) / (2 - mx - mn) : (mx - mn) / (mx + mn);
  return { s, l };
}
const hex = (c: number[]) => "#" + c.map(v => Math.round(Math.max(0, Math.min(255, v))).toString(16).padStart(2, "0")).join("").toUpperCase();
const mix = (a: number[], b: number[], t: number) => a.map((v, i) => v * (1 - t) + b[i] * t);

function remap(h: string, plot: boolean): string {
  if (h.toLowerCase() === "000000") return "#FFFFFF";
  const c = [0, 2, 4].map(i => parseInt(h.slice(i, i + 2), 16));
  const { s, l } = hsl(c[0], c[1], c[2]);
  if (l < 0.28) return INK;
  if (l > 0.9 && s < 0.35) return plot ? (l > 0.985 ? "#1E2128" : "#2A2E38") : "#" + h;
  if (l >= 0.62) return hex(mix(c, BG, 0.5));
  if (s < 0.2) return hex(mix(c, [255, 255, 255], 0.25));          // greys (axis text): a little lighter
  return hex(mix(c, [255, 255, 255], 0.3));
}

export function darkColors(text: string, plot = false, keep: readonly string[] = []): string {
  const kept = new Set(keep.map(c => c.replace("#", "").toLowerCase()));
  return text.replace(/#([0-9a-fA-F]{6})\b/g, (_, h) => (kept.has(h.toLowerCase()) ? "#" + h : remap(h, plot)));
}
