export type Field = {
  type: string; id?: string; label?: string; optional?: boolean; adduct?: boolean; multi?: boolean;
  value?: any; step?: number | null; placeholder?: string | null; options?: string[]; min?: number; max?: number;
  accept?: string; columns?: string[]; caption?: string | null; adducts?: Record<string, string[]>;
};
export type Block = {
  type: string; title?: string | null; columns?: string[]; caption?: string | null;
  items?: { label: string; unit?: string }[]; label?: string; key?: string | null; id?: string | null;
};
export type Module = {
  id: string; title: string; cls: string; group: string; subgroup: string | null; live: boolean;
  one_liner: string; inputs: Field[]; outputs: Block[]; example: Record<string, any> | null;
};
export type Stage = { title: string; modules: string[]; exclude: string[] };
export type Pipeline = { id: string; title: string; one_liner: string; inputs: Field[]; stages: Stage[] };
export type Group = { id: string; title: string };
export type Spec = { groups: Group[]; subgroups: Group[]; pipelines: Pipeline[]; modules: Module[] };

/** 后端结果：按输出块 id 填数据（形状见 mzdesk/api/common.py）。 */
export type Cell = string | { img: string };
export type CardValue = string | { value: string; note?: string };
export type PeakCard = { mz: number; rel: number; on: boolean; kind: string; rows: [string, string][]; img?: string };
export type BlockData = {
  rows?: Cell[][]; values?: Record<string, CardValue>; text?: string; formats?: Record<string, string>; items?: string[];
  plotly?: { data: any[]; layout: any }; peaks?: PeakCard[]; svg?: string; img?: string; mol3d?: string;
};
export type RunResult = { status: "ok" | "placeholder"; outputs: Record<string, BlockData>; warnings: string[] };
