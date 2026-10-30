export type FieldCondition = { field: string; equals?: any; in?: any[]; empty?: boolean };
export type Field = {
  type: string; id?: string; label?: string; optional?: boolean; adduct?: boolean; mz?: boolean; multi?: boolean;
  value?: any; step?: number | null; placeholder?: string | null; options?: string[]; min?: number; max?: number;
  rows?: number; accept?: string; columns?: string[]; caption?: string | null; adducts?: Record<string, string[]>;
  binary?: boolean; max_bytes?: number;
  // form metadata (Spectrum pages; optional everywhere)
  section?: string; form_order?: number; required?: boolean; required_when?: FieldCondition; visible_when?: FieldCondition;
  header_fallback?: boolean; auto_label?: string; ion_modes?: string[]; allow_formula?: boolean; help?: Help;
  label_field?: boolean;           // structure: an optional label input inside the field names the compound on the result cards
  // internet use (spec.py): the frontend asks before running online while the sidebar switch is off
  network?: boolean;               // checkbox: ticking it uses the internet
  network_sources?: string[];      // spectrum: these input sources are fetched online
  network_modes?: string[];        // structure: these input modes are resolved online
  network_options?: string[];       // multiselect: these options use the internet
};
/** Help text for a "?" popover: one string, or a list shown as bullets. */
export type Help = string | string[];
export type Block = {
  type: string; title?: string | null; columns?: string[]; caption?: string | null;
  items?: { label: string; unit?: string; help?: Help }[]; label?: string; key?: string | null; id?: string | null;
  help?: Help;              // any block: static "?" popover next to the title (runtime data.help is added to it)
  style?: string;           // cards: "metric" = compact metric cards
  row_help?: boolean;       // table / grouptable: trailing "?" column for per-row help
  page_size?: number;       // opt-in pagination and preview-only filtering for large tables
  compact?: boolean;        // long text can be expanded without losing the full value
  collapsed?: boolean;      // notes: items live behind a "?" next to the title (or a "Details" disclosure without a title)
  span?: "half" | "full";   // results grid on wide screens; default chosen by block type (Outputs.tsx blockLayout)
};
export type Module = {
  id: string; title: string; cls: string; group: string; subgroup: string | null; live: boolean;
  one_liner: string; inputs: Field[]; outputs: Block[]; example: Record<string, any> | null;
  network?: string;          // what the page sends online; present only on pages that use the internet
};
export type Stage = { id: number; title: string; modules: string[] };
export type PipelineExample = { label: string; values: Record<string, any> };
export type Pipeline = { id: string; title: string; one_liner: string; inputs: Field[]; stages: Stage[]; examples?: PipelineExample[] };
export type Group = { id: string; title: string };
export type Spec = { groups: Group[]; subgroups: Group[]; pipelines: Pipeline[]; modules: Module[]; version?: string };

/** 后端结果：按输出块 id 填数据（形状见 mzdesk/api/common.py）。 */
/** Table cell: plain text, an image, or an object with text / badge, optional link (url, or href), optional structure
 *  preview on hover (hover_svg: an SVG string), optional help. */
export type Cell = string | number | null
  | { img: string }
  | { fragment: FragmentTarget; help?: Help }
  | { text?: string; badge?: string; tone?: string; url?: string; href?: string; hover_svg?: string; help?: Help };
/** Table row: a list of cells, or {cells, help} when the row has its own "?" (spec row_help=true). */
export type Row = Cell[] | { cells: Cell[]; help?: Help };
export type CardValue = string | { value: string; note?: string; help?: Help };
/** portal tile / chips chip / reaction option: see the helpers in mzdesk/spec.py for the contracts. */
export type Tile = { label: string; status: "yes" | "partial" | "no" | string; value?: string; chips?: string[]; help?: Help };
export type StatusChip = { label: string; level: string; value?: string; help?: Help };
export type RxStep = { from: string; reagent?: string; to: string; count?: number };
export type RxOption = { label: string; steps: RxStep[]; formula?: string; mass?: string | number; heavy?: string | number; help?: Help; svg?: string; product_svg?: string };
export type RowGroup = { title: string; meta?: string; help?: Help; rows: Row[] };
/** color / label: a peak drawn in its fragment's colour with its own two-line label (m/z, Δ ppm) in the fragment view (Figure.tsx). */
export type PeakCard = { mz: number; rel: number; on: boolean; kind: string; rows: [string, string][]; img?: string; color?: string; label?: string };
/** Fragment view (mzdesk/ms2/fragment_view.py): a "View" cell names one structure; the table's `fragments` is the spectrum they share. */
export type FragmentTarget = { smiles: string; adduct: string; name?: string };
export type FragmentQuery = { peaks: [number, number][]; precursor_mz: number; tol_da: number };
export type FragmentItem = {
  color: string; mz: number; rel: number; ion_formula: string; calc_mz: number | null; ppm: number | null; nl: number | null;
  placed: boolean; alternatives: number; svg?: string; img?: string;
};
export type FragmentView = {
  name: string; smiles: string; formula: string; adduct: string; ion_formula: string; precursor_mz: number; theo_mz: number; positive: boolean;
  peaks: number; explained: number; explained_fraction: number; overview: string; fragments: FragmentItem[];
  figure: { plotly: { data: any[]; layout: any }; peaks: PeakCard[] }; notes: string[];
  settings: { tol_da: number; tol_ppm: number; tree_depth: number; engine: string }; seconds: number;
};
export type FragmentJob = { status: "queued" | "running" | "done" | "error"; percent: number; stage: string; seconds: number; error: string | null; log: string;
  result: FragmentView | null };
/** Structure & formats: one card per compound (2D SVG + 3D molblock) and the formats of each compound. */
export type CompareCompound = { label: string; formula?: string; svg?: string; mol3d?: string };
export type FormatSet = { label: string; formats: Record<string, string> };
export type BlockData = {
  rows?: Row[]; values?: Record<string, CardValue>; text?: string; formats?: Record<string, string>; items?: string[];
  files?: { name: string; content: string; mime?: string }[];
  plotly?: { data: any[]; layout: any }; peaks?: PeakCard[]; svg?: string; img?: string; mol3d?: string;
  help?: Help; count?: number | string;                       // any block: runtime "?" text; table / grouptable: count pill next to the title
  fragments?: FragmentQuery;                                   // table: the spectrum behind its "View" (fragment) cells
  columns?: string[];                                          // table: replaces the spec columns for this run (rows then have one cell per data column)
  tiles?: Tile[]; chips?: StatusChip[]; options?: RxOption[]; groups?: RowGroup[];   // portal / chips / reaction (svg above) / grouptable
  compounds?: CompareCompound[];                               // compare: one card per compound
  sets?: FormatSet[];                                          // formats: one tab per compound (single compound = the old shape)
};
/** verdict: the one-line answer shown above the results (may use <b>); none for pages whose output is the answer. */
export type RunResult = { status: "ok" | "placeholder"; outputs: Record<string, BlockData>; warnings: string[]; verdict?: string | null };

/** ESI MS/MS Annotation (mzdesk/annotation/pipeline.py): one step of a running job, with the blocks it has produced (spec block + its data). */
export type JobBlock = Block & { data: BlockData };
export type JobStep = {
  id: number; title: string; status: "queued" | "running" | "done" | "skipped" | "failed"; seconds: number; note: string; blocks: JobBlock[];
};
export type AnnotationAnswer = {
  type: "single" | "group" | "formula_only" | "mass_only" | "none"; ids: string[]; layer: string; layer_label: string; level: string; level_label: string;
  statement: string; top_id: string | null; isomer_verdict: string | null;
};
export type AnnotationCandidate = {
  id: string; rank: number | null; name: string; smiles: string | null; formula: string | null; adduct: string | null; layer: string; layer_label: string;
  level: string | null; level_label: string; level_reason: string; final_score: number | null; base_score: number | null; hypothesis: boolean;
  evidence?: { library?: { gate_passed?: boolean } };
  detect: { omics?: string; ionizable?: string; adduct?: string } | null;
  fragment_adduct?: string | null;     // MetFrag's spelling of the adduct; null / absent = no "View" button
  isomer: { verdict: string; role: string | null; challenger: boolean; explained_fraction: number | null } | null;
};
export type AnnotationResult = {
  answer: AnnotationAnswer; candidates: AnnotationCandidate[]; warnings: string[];
  steps: { id: number; title: string; status: string; seconds: number; note: string }[];
  context: { precursor_mz: number; polarity: string; adduct: string; retained_adducts: string[]; omics: string; libraries: string[]; gate_passed: boolean };
  isomer: { verdict: string; top_id?: string; group: string[]; challengers: { id: string; name: string; peaks: { mz: number; rel?: number }[] }[];
    candidates?: { id: string; name: string; smiles: string; source: string; evaluated: boolean; explained_fraction: number | null }[];
    isomers?: { id: string; distinguished: boolean | null; distinguishing_peaks: { mz: number }[] }[]; notes: string[] } | null;
  visualization?: { query: BlockData | null; library: { figure: BlockData; library: string; accession: string; name: string } | null } | null;
  fragment_query?: FragmentQuery | null;
};
export type Job = { status: "queued" | "running" | "done" | "error"; note: string; steps: JobStep[]; result: AnnotationResult | null; error: string | null };
