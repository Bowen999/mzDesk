import { OutputList } from "../components/Outputs";
import type { Block, BlockData } from "../types";

/** Hidden dev page (/dev/blocks, not in the sidebar): every generic output block with sample data, then the same blocks before a run.
 *  The blocks and data below have the shapes of the spec helpers in mzdesk/spec.py and the builders in mzdesk/api/blocks.py. */

const PARENT_SVG = `<svg xmlns="http://www.w3.org/2000/svg" width="240" height="170" viewBox="0 0 240 170" font-family="Helvetica, Arial, sans-serif" font-size="15">
<circle cx="190" cy="85" r="19" fill="#fecaca" opacity=".8"/><circle cx="55" cy="34" r="16" fill="#bfdbfe" opacity=".8"/>
<g stroke="#1f2937" stroke-width="1.6" fill="none"><polygon points="140,85 120,119.6 80,119.6 60,85 80,50.4 120,50.4"/>
<line x1="140" y1="85" x2="172" y2="85"/><line x1="80" y1="50.4" x2="64" y2="42"/></g>
<text x="176" y="90" fill="#b91c1c">NH<tspan dy="4" font-size="10">2</tspan></text><text x="40" y="40" fill="#1d4ed8">OH</text>
<text x="120" y="152" text-anchor="middle" fill="#6b7280" font-size="12">4-aminophenol (sample)</text></svg>`;
const CELL_IMG = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(`<svg xmlns="http://www.w3.org/2000/svg" width="60" height="40"><polygon points="30,4 52,16 52,32 30,38 8,32 8,16" fill="none" stroke="#6b7280" stroke-width="2"/></svg>`);

const LONG_HELP = Array.from({ length: 40 }, (_, i) => `Condition line ${i + 1}: column ODS 4.6 mm x 150 mm, mobile phase acetonitrile / phosphate buffer (pH 3.0), flow 1.0 mL/min, detection 273 nm.`);

const BLOCKS: Block[] = [
  { type: "heading", title: "Cards", help: "A heading block can carry help too." },
  { type: "cards", id: "metric", style: "metric", title: "LC-ESI response (metric cards)", help: "Static help from the spec. Compact metric cards: muted label, 20 px value, small note line.",
    items: [{ label: "Polarity", unit: "+ / − / both", help: "Item help from the spec." }, { label: "ESI+ response", unit: "strong / moderate / weak" },
      { label: "ESI− response", unit: "strong / moderate / weak" }, { label: "Best adduct", unit: "" }] },
  { type: "cards", id: "single", style: "metric", title: "Metric cards with a single item (stays one column wide)", items: [{ label: "Polarity", unit: "+ / − / both" }] },
  { type: "cards", id: "plain", title: "Plain cards (unchanged)", items: [{ label: "Name", unit: "" }, { label: "Formula", unit: "" }, { label: "Exact mass", unit: "Da" }] },

  { type: "heading", title: "Portal, chips" },
  { type: "portal", id: "portal", title: "Derivatization routes", help: ["A portal is a grid of tiles.", "yes = teal, partial = amber, no = quiet grey."] },
  { type: "chips", id: "chips", title: "Detectors", help: "Chips: good / possible / poor (aliases strong, moderate, weak, yes, no). Chips with help are clickable." },

  { type: "heading", title: "Reaction" },
  { type: "reaction", id: "rx", title: "Derivatization reaction", help: "Only the functional groups are drawn, not the full product structures." },

  { type: "heading", title: "Tables" },
  { type: "table", id: "tbl", title: "Methods (links, badges, row help)", columns: ["Method", "Source", "Score", "Note"], row_help: true, help: "Rows may carry their own help (trailing ? column).", caption: "Caption under the table." },
  { type: "grouptable", id: "gt", title: "Methods by compound", columns: ["Method", "Source", "Score"], row_help: true, help: "Groups are collapsed by default; click a header or press Enter / Space on it." },
  { type: "table", id: "legacy", title: "Legacy table (plain strings, {text, href}, {img})", columns: ["Name", "Link", "Structure"] },
  { type: "table", id: "depict", title: "Cell with a 2D structure on hover ({text, hover_svg})", columns: ["Name", "SMILES"], help: "Hover or Tab to the SMILES: its 2D structure is drawn in the popover; click pins it." },

  { type: "reaction", id: "rx_empty", title: "Reaction with options: [] (renders nothing after a run)" },
  { type: "reaction", id: "rx_none", title: "Reaction with no data at all (renders nothing after a run)" },
  { type: "heading", title: "Long help (scrolls inside the popover)", help: LONG_HELP },
  { type: "heading", title: "Notes, figure" },
  { type: "notes", id: "n1", title: "Rule hits (collapsed: the items are in the ?)", collapsed: true },
  { type: "notes", id: "n2", title: null, collapsed: true },
  { type: "notes", id: "n3", title: "Rule hits (visible, unchanged)" },
  { type: "notes", id: "n4", title: "Collapsed notes with items: [] (renders nothing after a run)", collapsed: true },
  { type: "notes", id: "n5", title: "Collapsed notes with no data at all (renders nothing after a run)", collapsed: true },
  { type: "figure", id: "fig", label: "Figure with help", help: "Help also works on figures." },
];

const LONG = "Long sample text that wraps inside a table cell so that the layout of the link / help combination can be checked at width.";
const DATA: Record<string, BlockData> = {
  metric: { values: {
    Polarity: { value: "Positive", note: "basic amine, pKa 9.1", help: ["Runtime help for this item.", "A second bullet."] },
    "ESI+ response": "strong", "ESI− response": { value: "weak", note: "no acidic site" }, "Best adduct": "[M+H]+" } },
  single: { values: { Polarity: { value: "Positive", note: "basic amine, pKa 9.1" } } },
  plain: { values: { Name: "Caffeine", Formula: "C8H10N4O2", "Exact mass": { value: "194.0804", note: "Da" } } },
  portal: { tiles: [
    { label: "GC-MS (TMS)", status: "yes", value: "Silylation works", chips: ["–OH", "–NH₂"], help: "Two silylation sites: both groups are derivatized by MSTFA." },
    { label: "CIL LC-MS", status: "partial", value: "Amine only", chips: ["dansyl", "dabsyl"], help: ["Only the amine reacts.", "The phenol needs a second reagent."] },
    { label: "LC-UV / FLD", status: "no", value: "Native absorbance is enough" },
    { label: "Acylation", status: "yes", value: "Acetic anhydride" },
    { label: "Alkylation", status: "no", value: "No acidic carboxyl group", help: "Alkylation targets –COOH, which is absent." } ] },
  chips: { chips: [
    { label: "UV/PDA", level: "good", value: "λmax 273 nm", help: "Strong chromophore (xanthine ring)." },
    { label: "FLD", level: "possible", value: "weak", help: ["Native fluorescence unlikely.", "Derivatization would help."] },
    { label: "CAD", level: "strong" }, { label: "ELSD", level: "moderate" }, { label: "RI", level: "poor", value: "gradient incompatible" },
    { label: "Aliases", level: "weak" }, { label: "yes", level: "yes" }, { label: "no", level: "no", help: "No help text on most chips; this one has some." } ] },
  rx: { svg: PARENT_SVG, options: [
    { label: "CIL · dansyl", steps: [{ from: "R–NH₂", reagent: "Dns-Cl", to: "R–NH–Dns", count: 1 }, { from: "Ar–OH", reagent: "Dns-Cl", to: "Ar–O–Dns" }],
      formula: "C<sub>18</sub>H<sub>18</sub>N<sub>2</sub>O<sub>5</sub>S<sub>2</sub>", mass: 406.0657, heavy: 412.0858, help: ["Dansyl chloride labels amines and phenols.", "Heavy label: <sup>13</sup>C<sub>2</sub>-dansyl."] },
    { label: "GC-MS · TMS (MSTFA)", steps: [{ from: "R–NH₂", reagent: "MSTFA", to: "R–NH–TMS", count: 2 }, { from: "Ar–OH", reagent: "MSTFA", to: "Ar–O–TMS" }],
      formula: "C₁₂H₂₃NOSi₂", mass: "253.1318", help: "Trimethylsilylation of both groups." } ] },
  rx_empty: { options: [] },
  tbl: { count: 4, rows: [
    { cells: [{ text: "JP18 · Caffeine", url: "https://example.org/jp18" }, { badge: "JP", tone: "jp" }, "0.66", { text: "ODS, 273 nm", help: "Column and wavelength as in the monograph." }], help: ["Matches on the structure.", "Score = Tanimoto similarity."] },
    [{ text: "RepoRT 0042", url: "https://example.org/report/42" }, { badge: "RepoRT", tone: "report" }, "0.58", LONG],
    { cells: ["Local note", { badge: "neutral" }, "0.21", { img: CELL_IMG }], help: "Row help without any cell help." },
    ["Plain row", { badge: "ok", tone: "good", help: "A badge with help." }, { badge: "check", tone: "warn" }, { badge: "stop", tone: "bad", url: "https://example.org/stop" }] ] },
  gt: { count: 3, groups: [
    { title: "Caffeine", meta: "0.66 · 3 methods", help: "Best score 0.66 across three methods.", rows: [
      { cells: [{ text: "JP18 · caffeine", url: "https://example.org/a" }, { badge: "JP", tone: "jp" }, "0.66"], help: "Top hit." },
      [{ text: "RepoRT 0007", url: "https://example.org/b" }, { badge: "RepoRT", tone: "report" }, "0.61"], ["RepoRT 0031", { badge: "RepoRT", tone: "report" }, "0.55"]] },
    { title: "Theophylline", meta: "0.52 · 2 methods", rows: [["JP18 · theophylline", { badge: "JP", tone: "jp" }, "0.52"], ["RepoRT 0099", { badge: "RepoRT", tone: "report" }, "0.40"]] },
    { title: "Paracetamol", meta: "0.31 · 1 method", help: ["Weak match.", "Shown for completeness."], rows: [["RepoRT 0123", { badge: "RepoRT", tone: "report" }, "0.31"]] } ] },
  legacy: { rows: [["Caffeine", { text: "PubChem", href: "https://pubchem.ncbi.nlm.nih.gov/compound/2519" }, { img: CELL_IMG }], ["Plain", "no link", "no image"]] },
  depict: { rows: [
    ["Caffeine", { text: "CN1C=NC2=C1C(=O)N(C)C(=O)N2C", hover_svg: PARENT_SVG }],
    ["4-Aminophenol", { text: "NC1=CC=C(O)C=C1", hover_svg: PARENT_SVG, help: "hover_svg combines with help: the structure is drawn above the text." }],
    ["No drawing", "CN1C=NC2=C1C(=O)N(C)C(=O)N2C"] ] },
  n1: { items: ["Basic nitrogen (pKa 9.1) protonates at the chosen pH.", "No acidic site: ESI− weak."] },
  n2: { items: ["Untitled collapsed notes become a small Details disclosure."] },
  n3: { items: ["Visible note one.", "Visible note two."] },
  n4: { items: [] },
  fig: { svg: PARENT_SVG, help: "Runtime help on a figure." },
};

export function DevBlocks() {
  const h2 = { fontSize: 16, margin: "28px 0 12px" } as const;
  return (
    <div className="page">
      <h1>Output blocks (dev)</h1>
      <p className="one-liner">Every generic output block with sample data. Hover or Tab to a ? to preview it, click to pin it, click outside or press Esc to close.</p>
      <h2 style={h2}>With data</h2>
      <OutputList blocks={BLOCKS} outputs={DATA} />
      <h2 style={h2}>Before the first run (no data)</h2>
      <OutputList blocks={BLOCKS} outputs={null} />
    </div>
  );
}
