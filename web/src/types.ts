export type Field = {
  type: string; id?: string; label?: string; optional?: boolean; adduct?: boolean; multi?: boolean;
  value?: any; step?: number | null; placeholder?: string | null; options?: string[]; min?: number; max?: number;
  accept?: string; columns?: string[]; caption?: string | null;
};
export type Block = {
  type: string; title?: string | null; columns?: string[]; caption?: string | null;
  items?: { label: string; unit?: string }[]; label?: string; key?: string | null;
};
export type Module = { id: string; title: string; cls: string; group: string; one_liner: string; inputs: Field[]; outputs: Block[] };
export type Stage = { title: string; module: string; exclude: string[] };
export type Pipeline = { id: string; title: string; one_liner: string; inputs: Field[]; stages: Stage[] };
export type Spec = { groups: { id: string; title: string }[]; pipelines: Pipeline[]; modules: Module[] };
