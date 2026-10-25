export type Field = {
  type: string; id?: string; label?: string; optional?: boolean; adduct?: boolean; multi?: boolean;
  value?: any; step?: number | null; placeholder?: string | null; options?: string[]; min?: number; max?: number;
  accept?: string; columns?: string[]; caption?: string | null; adducts?: Record<string, string[]>;
};
/** Help text for a "?" popover: one string, or a list shown as bullets. */
export type Help = string | string[];
export type Block = {
  type: string; title?: string | null; columns?: string[]; caption?: string | null;
  items?: { label: string; unit?: string; help?: Help }[]; label?: string; key?: string | null; id?: string | null;
  help?: Help;              // any block: static "?" popover next to the title (runtime data.help is added to it)
  style?: string;           // cards: "metric" = compact metric cards
  row_help?: boolean;       // table / grouptable: trailing "?" column for per-row help
  collapsed?: boolean;      // notes: items live behind a "?" next to the title (or a "Details" disclosure without a title)
  span?: "half" | "full";   // results grid on wide screens; default chosen by block type (Outputs.tsx blockLayout)
};
export type Module = {
  id: string; title: string; cls: string; group: string; subgroup: string | null; live: boolean;
  one_liner: string; inputs: Field[]; outputs: Block[]; example: Record<string, any> | null;
};
export type Stage = { title: string; modules: string[]; exclude: string[] };
export type Pipeline = { id: string; title: string; one_liner: string; inputs: Field[]; stages: Stage[] };
export type Group = { id: string; title: string };
export type Spec = { groups: Group[]; subgroups: Group[]; pipelines: Pipeline[]; modules: Module[]; version?: string };

/** 后端结果：按输出块 id 填数据（形状见 mzdesk/api/common.py）。 */
/** Table cell: plain text, an image, or an object with text / badge, optional link (url, or href), optional help. */
export type Cell = string | number | null
  | { img: string }
  | { text?: string; badge?: string; tone?: string; url?: string; href?: string; help?: Help };
/** Table row: a list of cells, or {cells, help} when the row has its own "?" (spec row_help=true). */
export type Row = Cell[] | { cells: Cell[]; help?: Help };
export type CardValue = string | { value: string; note?: string; help?: Help };
/** portal tile / chips chip / reaction option: see the helpers in mzdesk/spec.py for the contracts. */
export type Tile = { label: string; status: "yes" | "partial" | "no" | string; value?: string; chips?: string[]; help?: Help };
export type StatusChip = { label: string; level: string; value?: string; help?: Help };
export type RxStep = { from: string; reagent?: string; to: string; count?: number };
export type RxOption = { label: string; steps: RxStep[]; formula?: string; mass?: string | number; heavy?: string | number; help?: Help };
export type RowGroup = { title: string; meta?: string; help?: Help; rows: Row[] };
export type PeakCard = { mz: number; rel: number; on: boolean; kind: string; rows: [string, string][]; img?: string };
export type BlockData = {
  rows?: Row[]; values?: Record<string, CardValue>; text?: string; formats?: Record<string, string>; items?: string[];
  plotly?: { data: any[]; layout: any }; peaks?: PeakCard[]; svg?: string; img?: string; mol3d?: string;
  help?: Help; count?: number | string;                       // any block: runtime "?" text; table / grouptable: count pill next to the title
  columns?: string[];                                          // table: replaces the spec columns for this run (rows then have one cell per data column)
  tiles?: Tile[]; chips?: StatusChip[]; options?: RxOption[]; groups?: RowGroup[];   // portal / chips / reaction (svg above) / grouptable
};
/** verdict: the one-line answer shown above the results (may use <b>); none for pages whose output is the answer. */
export type RunResult = { status: "ok" | "placeholder"; outputs: Record<string, BlockData>; warnings: string[]; verdict?: string | null };
