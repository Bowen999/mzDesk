import { lazy, Suspense, useEffect, useId, useRef, useState } from "react";
import type { Field, FieldCondition } from "../types";
import { t, useT } from "../i18n";
import { countPeaks } from "../peakCount";
import { ensureSetup } from "../setup";

// Ketcher 很大，只在选 Draw 时才下载
const StructureEditor = lazy(() => import("./StructureEditor"));

/** The wire value of the first structure mode. The visible label is "Name / ID" on pages that
 *  cannot accept a formula; the value stays the same so the backend contract does not change. */
const NAME_MODE = "Name / ID / Formula";
const LIPID_MODE = "Lipid name";              // Lipid Fragment Explanation: a lipid shorthand, read locally (spec.py LIPID_NAME_MODE)
const POSITIVE = "Positive";
const NEGATIVE = "Negative";

/** The shown name of an ion mode (the value sent stays "Positive" / "Negative"), of a spectrum source and of "Auto" (labels only, never values). */
const ionLabel = (m: string) => (m === POSITIVE ? t("Positive") : m === NEGATIVE ? t("Negative") : m);
const sourceLabel = (s: string) => (s === "Paste" ? t("Paste") : s === "MGF / MSP file" ? t("MGF / MSP file") : s);
/** "label: value" in summaries. */
const pair = (label: string, value: string | number) => t("{label}: {value}", { label, value });

/** 读上传文件的文本内容（多个文件用空行连接）。 */
async function readFiles(files: FileList | null): Promise<{ name: string; text: string }> {
  const list = Array.from(files ?? []);
  const texts = await Promise.all(list.map(f => f.text()));
  return { name: list.map(f => f.name).join(", "), text: texts.join("\n\n") };
}

/** 每个字段的初值。复合字段（structure / precursor / spectrum）是对象。 */
export function initialValues(fields: Field[]): Record<string, any> {
  const v: Record<string, any> = {};
  for (const f of fields) {
    if (!f.id) continue;
    switch (f.type) {
      case "structure": v[f.id] = { mode: f.modes?.[0] ?? NAME_MODE, value: "", label: "", buffers: {} }; break;
      case "precursor": {
        const mode = f.ion_modes?.[0] ?? POSITIVE;
        v[f.id] = { precursor_mz: "", ion_mode: mode, adduct: f.adducts?.[mode]?.[0] ?? "Auto" };
        break;
      }
      case "spectrum": v[f.id] = { source: "Paste", value: "", buffers: {} }; break;
      case "checkbox": v[f.id] = !!f.value; break;
      case "multiselect": v[f.id] = f.value ?? []; break;
      default: v[f.id] = f.value ?? "";
    }
  }
  return v;
}

// ---------------------------------------------------------------- conditions / validation
/** True when a field's `visible_when` / `required_when` condition holds for `values`. */
export function conditionMatches(c: FieldCondition, values: Record<string, any>): boolean {
  const v = values[c.field];
  if (c.empty !== undefined) {
    const empty = v === undefined || v === null || (typeof v === "string" && !v.trim()) || (Array.isArray(v) && !v.length);
    return empty === c.empty;
  }
  if ("equals" in c) return v === c.equals;
  if (c.in) return (c.in as any[]).includes(v);
  return true;
}

/** Fields that the chosen mode actually uses. */
export function visibleFields(fields: Field[], values: Record<string, any>): Field[] {
  return fields.filter(f => !f.visible_when || conditionMatches(f.visible_when, values));
}

/** A preview of the precursor written in an MGF / MSP / MassIVE / SIRIUS header. */
export function headerPrecursor(text: string): number | null {
  if (!text) return null;
  const m = /^\s*(?:(?:PEPMASS|PRECURSORMZ|PRECURSOR_MZ)\s*[=:]|PRECURSOR\s+M\/Z\s*[=:]|MS\$FOCUSED_ION:\s*PRECURSOR_M\/Z|>parentmass)\s*([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)/im.exec(text);
  if (!m) return null;
  const n = Number(m[1]);
  return Number.isFinite(n) && n > 0 ? n : null;
}

export function isRequired(f: Field, values: Record<string, any>): boolean {
  if (f.required) return true;
  return !!f.required_when && conditionMatches(f.required_when, values);
}

function headerSatisfies(field: Field, fields: Field[], values: Record<string, any>): boolean {
  if (!field.header_fallback) return false;
  for (const f of fields) {
    if (f.type !== "spectrum" || !f.id) continue;
    const v = values[f.id] ?? {};
    if (v.source === "USI" && String(v.value ?? "").trim()) return true;
    if (headerPrecursor(String(v.value ?? "")) !== null) return true;
  }
  return false;
}

/** A composite value is missing when the part the backend needs is empty. */
export function isFilled(f: Field, value: any, fields: Field[], values: Record<string, any>): boolean {
  switch (f.type) {
    case "heading": return true;
    case "slider": return true;
    case "structure": return !!String(value?.value ?? "").trim();
    case "precursor":
      if (f.mz === false) return true;                       // adduct-only field
      if (String(value?.precursor_mz ?? "").trim()) return true;
      return headerSatisfies(f, fields, values);
    case "spectrum": return !!String(value?.value ?? "").trim();
    case "file": return !!value?.name && (f.binary ? !!value?.base64 : typeof value?.text === "string");
    case "multiselect": return Array.isArray(value) && value.length > 0;
    case "checkbox": return !!value;
    case "table": return !!String(value ?? "").trim();
    case "ms1": return !!String(value ?? "").trim();
    default: return value !== undefined && value !== null && String(value).trim() !== "";
  }
}

/** Required fields that have no usable value. */
export function missingRequired(fields: Field[], values: Record<string, any>): Field[] {
  return visibleFields(fields, values).filter(f =>
    f.id && f.type !== "heading" && isRequired(f, values) && !isFilled(f, values[f.id] ?? initialValues([f])[f.id], fields, values));
}

// ---------------------------------------------------------------- field chrome
function Req() {
  const t = useT();
  return <span className="req" aria-label={t("required")}>*</span>;
}

function Note({ text, error, id }: { text?: string; error?: string; id?: string }) {
  return (
    <>
      {text && <small className="field-help" id={id}>{text}</small>}
      {error && <span className="field-error" role="alert">{error}</span>}
    </>
  );
}

type P = { field: Field; value: any; onChange: (v: any) => void; error?: string };

export function FieldInput({ field, value, onChange, error }: P) {
  const t = useT();
  const f = field;
  if (f.id && value === undefined) value = initialValues([f])[f.id];   // 首次渲染或切页时 state 还没初始化
  const uid = useId();
  const inputId = `field-${f.id ?? uid}-${uid}`;
  const help = typeof f.help === "string" ? f.help : undefined;
  const required = !!f.required || !!f.required_when;
  const common = {
    id: inputId,
    "aria-required": required || undefined,
    "aria-invalid": error ? true : undefined,
    "aria-describedby": help ? `${inputId}-help` : undefined,
  };

  // Ticking a box whose option needs components (spec `needs`): download what is missing first (setup.ts); "Not now" leaves it unticked.
  const tick = (on: boolean) => {
    if (!on || !f.needs?.length) { onChange(on); return; }
    void ensureSetup({ ids: f.needs }).then(ok => { if (ok) onChange(true); });
  };

  switch (f.type) {
    case "heading": return <h4 className="field-heading">{f.label}</h4>;
    case "structure": return <StructureInput field={f} value={value} onChange={onChange} error={error} />;
    case "precursor": return <PrecursorInput field={f} value={value} onChange={onChange} error={error} />;
    case "spectrum": return <SpectrumInput field={f} value={value} onChange={onChange} error={error} />;
    case "ms1": return (
      <label className="field wide" htmlFor={inputId}><span className="field-label-text">{f.label}{required && <Req />}</span>
        <textarea {...common} rows={4} placeholder={f.placeholder ?? ""} value={value} onChange={e => onChange(e.target.value)} />
        <Note text={help} error={error} id={`${inputId}-help`} />
      </label>);
    case "number": return (
      <label className="field" htmlFor={inputId}><span className="field-label-text">{f.label}{required && <Req />}</span>
        <input {...common} type="number" min={f.min} max={f.max} step={f.step ?? "any"} placeholder={f.placeholder ?? ""}
          value={value ?? ""} onChange={e => onChange(e.target.value)} />
        <Note text={help} error={error} id={`${inputId}-help`} />
      </label>);
    case "text": return (
      <label className="field" htmlFor={inputId}><span className="field-label-text">{f.label}{required && <Req />}</span>
        <input {...common} type="text" placeholder={f.placeholder ?? ""} value={value ?? ""} onChange={e => onChange(e.target.value)} />
        <Note text={help} error={error} id={`${inputId}-help`} />
      </label>);
    case "textarea": return (
      <label className="field wide" htmlFor={inputId}><span className="field-label-text">{f.label}{required && <Req />}</span>
        <textarea {...common} rows={f.rows ?? 5} placeholder={f.placeholder ?? ""} value={value ?? ""} onChange={e => onChange(e.target.value)} />
        <Note text={help} error={error} id={`${inputId}-help`} />
      </label>);
    case "select": return (
      <label className="field" htmlFor={inputId}><span className="field-label-text">{f.label}{required && <Req />}</span>
        <select {...common} value={value} onChange={e => onChange(e.target.value)}>{f.options!.map(o => <option key={o} value={o}>{o === "" ? (f.placeholder ?? "") : t(o)}</option>)}</select>
        <Note text={help} error={error} id={`${inputId}-help`} />
      </label>);
    case "radio": return (
      <div className="field" role="radiogroup" aria-labelledby={`${inputId}-label`}>
        <span id={`${inputId}-label`}><span className="field-label-text">{f.label}{required && <Req />}</span></span>
        <div className="radio-row">
          {f.options!.map(o => <label key={o}><input type="radio" name={inputId} checked={value === o} onChange={() => onChange(o)} /> {t(o)}</label>)}
        </div>
        <Note text={help} error={error} id={`${inputId}-help`} />
      </div>);
    case "multiselect": return (
      <div className="field wide" role="group" aria-labelledby={`${inputId}-label`}>
        <span id={`${inputId}-label`}><span className="field-label-text">{f.label}{required && <Req />}</span></span>
        {f.checkboxes ? (
          <div className="check-grid">
            {f.options!.map(o => {
              const on = (value as string[]).includes(o);
              // keep the options' order whatever order they are ticked in
              return <label className="check-line" key={o}><input type="checkbox" checked={on}
                onChange={() => onChange(f.options!.filter(x => x === o ? !on : (value as string[]).includes(x)))} /> {t(o)}</label>;
            })}
          </div>
        ) : (
        <div className="chip-row">
          {f.options!.map(o => {
            const on = (value as string[]).includes(o);
            return <button type="button" key={o} className={"chip" + (on ? " on" : "")} aria-pressed={on}
              onClick={() => onChange(on ? value.filter((x: string) => x !== o) : [...value, o])}>{t(o)}</button>;
          })}
          {f.options!.length > 8 && (
            <span className="chip-actions">
              <button type="button" className="link-btn" onClick={() => onChange([...f.options!])}>{t("All")}</button>
              <button type="button" className="link-btn" onClick={() => onChange([])}>{t("None")}</button>
            </span>)}
        </div>)}
        <Note text={help} error={error} id={`${inputId}-help`} />
      </div>);
    case "checkbox": return (
      <label className="field check" htmlFor={inputId}>
        <span className="check-line"><input {...common} type="checkbox" checked={!!value} onChange={e => tick(e.target.checked)} /> {f.label}{required && <Req />}</span>
        <Note text={help} error={error} id={`${inputId}-help`} />
      </label>);
    case "slider": return (
      <label className="field" htmlFor={inputId}><span className="field-label-text">{f.label}{required && <Req />}</span>
        <span className="slider-label"><span /> <b>{value}</b></span>
        <input {...common} type="range" min={f.min} max={f.max} step={f.step ?? 1} value={value}
          aria-valuetext={String(value)} onChange={e => onChange(Number(e.target.value))} />
        <Note text={help} error={error} id={`${inputId}-help`} />
      </label>);
    case "file": return <FileInput field={f} value={value} onChange={onChange} error={error} inputId={inputId} required={required} />;
    case "table": return (
      <label className="field wide" htmlFor={inputId}><span className="field-label-text">{f.label}{required && <Req />}</span>
        <div className="grid-head">{f.columns!.join("  |  ")}</div>
        <textarea {...common} rows={5} placeholder={t("One row per line: space, comma or tab between columns")} value={value} onChange={e => onChange(e.target.value)} />
        {f.caption && <small>{f.caption}</small>}
        <Note text={help} error={error} id={`${inputId}-help`} />
      </label>);
    default: return <div className="field">unknown field type {f.type}</div>;
  }
}

function FileInput({ field, value, onChange, error, inputId, required }: P & { inputId: string; required: boolean }) {
  const t = useT();
  const [fileError, setFileError] = useState("");
  const sequence = useRef(0);
  useEffect(() => () => { sequence.current++; }, []); // Ignore file reads that finish after Clear unmounts this field.
  const help = typeof field.help === "string" ? field.help : undefined;
  const select = async (files: FileList | null) => {
    const id = ++sequence.current;
    const file = files?.[0];
    setFileError("");
    onChange("");
    if (!file) return;
    if (field.max_bytes && (file.size === 0 || file.size > field.max_bytes)) {
      setFileError(t("Choose a non-empty file of at most {mb} MB.", { mb: field.max_bytes / 1_000_000 }));
      return;
    }
    try {
      if (!field.binary) {
        const result = await readFiles(files);
        if (id === sequence.current) onChange(result);
        return;
      }
      const data = await new Promise<string>((resolve, reject) => {
        const reader = new FileReader();
        reader.onload = () => resolve(String(reader.result ?? ""));
        reader.onerror = () => reject(reader.error);
        reader.readAsDataURL(file);
      });
      const comma = data.indexOf(",");
      if (comma < 0) throw new Error("Invalid file data");
      if (id === sequence.current) onChange({ name: file.name, base64: data.slice(comma + 1) });
    } catch {
      if (id === sequence.current) setFileError(t("Could not read the file. Please choose it again."));
    }
  };
  return <label className="field" htmlFor={inputId}>
    <span className="field-label-text">{field.label}{required && <Req />}</span>
    <input id={inputId} type="file" accept={field.accept} aria-required={required || undefined}
      aria-invalid={!!(error || fileError) || undefined} onChange={e => void select(e.target.files)} />
    {value?.name && <small>{value.name}</small>}
    <Note text={help} error={fileError || error} />
  </label>;
}

function StructureInput({ field, value, onChange, error }: P) {
  const t = useT();
  const uid = useId();
  const labelId = `structure-${field.id ?? uid}-label`;
  const name = `structure-mode-${field.id ?? "field"}-${uid}`;
  const allowFormula = field.allow_formula === true;
  const nameLabel = allowFormula ? t("Name / ID / Formula") : t("Name / ID");
  const modeLabels: Record<string, string> = { [NAME_MODE]: nameLabel, [LIPID_MODE]: t("Lipid name"), Draw: t("Draw"), "SMILES / InChI": "SMILES / InChI" };
  const modes: [string, string][] = (field.modes ?? [NAME_MODE, "Draw", "SMILES / InChI"]).map(m => [m, modeLabels[m] ?? m]);     // [value sent, label shown]
  const namePlaceholder = value.mode === LIPID_MODE
    ? t("Lipid short name, e.g. {example}", { example: field.placeholder || "PC 16:0/18:1" })
    : allowFormula
      ? t("Name, CAS, formula, InChIKey or database ID, e.g. caffeine")
      : t("Name, CAS, InChIKey or database ID, e.g. caffeine");
  const set = (patch: Partial<typeof value>) => onChange({ ...value, ...patch });
  // Keep one buffer per mode, so Name → Draw → Name does not erase the typed identifier.
  const switchMode = (m: string) => {
    const buffers = { ...(value.buffers ?? {}), [value.mode]: { value: value.value ?? "", molfile: value.molfile ?? "" } };
    const next = buffers[m] ?? { value: "", molfile: "" };
    onChange({ ...value, mode: m, value: next.value ?? "", molfile: next.molfile ?? "", buffers });
  };
  return (
    <div className="field wide composite" role="group" aria-labelledby={labelId}>
      <span id={labelId}>{field.label}{field.required && <Req />}{field.optional ? " " + t("(optional)") : ""}</span>
      <div className="radio-row">{modes.map(([m, label]) => <label key={m}><input type="radio" name={name} checked={value.mode === m} onChange={() => switchMode(m)} /> {label}</label>)}</div>
      {value.mode === "Draw" ? (
        <>
          <Suspense fallback={<div className="placeholder-box editor">{t("Loading structure editor…")}</div>}>
            <StructureEditor smiles={value.value} molfile={value.molfile} onChange={v => set(v)} />
          </Suspense>
          <input type="text" className="mono" readOnly placeholder={t("SMILES from editor")} value={value.value} />
        </>
      ) : value.mode === "SMILES / InChI" ? (
        <input type="text" className="mono" aria-required={field.required || undefined} aria-invalid={error ? true : undefined}
          placeholder="CN1C=NC2=C1C(=O)N(C(=O)N2C)C" value={value.value} onChange={e => set({ value: e.target.value })} />
      ) : (
        <input type="text" aria-required={field.required || undefined} aria-invalid={error ? true : undefined}
          placeholder={namePlaceholder} value={value.value} onChange={e => set({ value: e.target.value })} />
      )}
      {field.label_field && (
        <label className="structure-label">
          <span>{t("Label (optional)")}</span>
          <input type="text" placeholder={t("Name shown on the results")}
            value={value.label ?? ""} onChange={e => set({ label: e.target.value })} />
        </label>
      )}
      <Note text={typeof field.help === "string" ? field.help : undefined} error={error} />
    </div>
  );
}

function PrecursorInput({ field, value, onChange, error }: P) {
  const t = useT();
  const uid = useId();
  const labelId = `precursor-${field.id ?? uid}-label`;
  const name = `ion-mode-${field.id ?? "field"}-${uid}`;
  const set = (patch: Partial<typeof value>) => onChange({ ...value, ...patch });
  const lists = field.adducts ?? { Positive: ["Auto"], Negative: ["Auto"] };
  const modes = field.ion_modes ?? [POSITIVE, NEGATIVE];
  const adducts = lists[value.ion_mode] ?? lists[modes[0]] ?? ["Auto"];
  const label = (a: string) => a === "Auto" ? (field.auto_label || t("Auto")) : a;
  const help = typeof field.help === "string" ? field.help : undefined;
  return (
    <div className="field wide composite row3" role="group" aria-label={field.label ?? t("Precursor and adduct")}>
      {field.mz !== false && (
        <label className="field"><span className="field-label-text">{field.label ?? t("Precursor m/z")}{field.required && <Req />}</span>
          <input type="number" step="0.0001" placeholder="195.0877" value={value.precursor_mz}
            aria-required={field.required || undefined} aria-invalid={error ? true : undefined}
            onChange={e => set({ precursor_mz: e.target.value })} />
        </label>)}
      <div className="field"><span id={labelId}>{t("Ion mode")}</span>
        {modes.length === 1
          ? <span className="field-static">{t("{mode} only", { mode: ionLabel(modes[0]) })}</span>
          : <div className="radio-row">
            {modes.map(m => <label key={m}><input type="radio" name={name} checked={value.ion_mode === m}
              onChange={() => set({ ion_mode: m, adduct: lists[m]?.[0] ?? "Auto" })} /> {ionLabel(m)}</label>)}
          </div>}
      </div>
      {field.adduct !== false && <label className="field">{t("Adduct")}<select value={value.adduct} onChange={e => set({ adduct: e.target.value })}>{adducts.map(a => <option key={a} value={a}>{label(a)}</option>)}</select></label>}
      <div className="field-note-wide"><Note text={help} error={error} /></div>
    </div>
  );
}

function SpectrumInput({ field, value, onChange, error }: P) {
  const t = useT();
  const uid = useId();
  const labelId = `spectrum-${field.id ?? uid}-label`;
  const name = `spectrum-source-${field.id ?? "field"}-${uid}`;
  const fileSequence = useRef(0);
  useEffect(() => () => { fileSequence.current++; }, []);
  const sources = ["Paste", "MGF / MSP file", "USI"];
  const buffers = value.buffers ?? {};
  const switchSource = (s: string) => {
    fileSequence.current++;
    const next = { ...buffers, [value.source]: { value: value.value ?? "", filename: value.filename ?? "" } };
    const b = next[s] ?? { value: "", filename: "" };
    onChange({ ...value, source: s, value: b.value ?? "", filename: b.filename ?? "", buffers: next });
  };
  const set = (patch: Partial<typeof value>) => onChange({ ...value, ...patch });
  const peaks = countPeaks(String(value.value ?? ""));
  const help = typeof field.help === "string" ? field.help : undefined;
  return (
    <div className="field wide composite" role="group" aria-labelledby={labelId}>
      <span id={labelId}>{field.label}{field.required && <Req />}</span>
      <div className="radio-row">{sources.map(s => <label key={s}><input type="radio" name={name} checked={value.source === s} onChange={() => switchSource(s)} /> {sourceLabel(s)}</label>)}</div>
      {value.source === "Paste" && (
        <>
          <textarea rows={5} aria-required={field.required || undefined} aria-invalid={error ? true : undefined}
            placeholder={t("m/z intensity per line, or MGF, MSP, JSON …")} value={value.value}
            onChange={e => set({ value: e.target.value })} />
          {peaks > 0 && <small>{t("{n} peaks read", { n: peaks })}</small>}
        </>
      )}
      {value.source === "MGF / MSP file" && (
        <div className="file-row">
          <input type="file" accept=".mgf,.msp,.txt,.csv,.json,.xml,.ms" multiple={!!field.multi}
            onChange={async e => {
              const mine = ++fileSequence.current, r = await readFiles(e.target.files);
              if (mine === fileSequence.current) set({ value: r.text, filename: r.name });
            }} />
          {value.filename && <small>{value.filename}{peaks > 0 ? " · " + t("{n} peaks", { n: peaks }) : ""}</small>}
        </div>)}
      {value.source === "USI" && (
        <>
          <input type="text" placeholder="mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB…" value={value.value}
            aria-required={field.required || undefined} aria-invalid={error ? true : undefined}
            onChange={e => set({ value: e.target.value })} />
        </>
      )}
      <Note text={help} error={error} />
    </div>
  );
}

/** An option as the input bar shows it: in the language of the page (the value itself, which is what was sent, is not changed). */
const optionText = (f: Field, v: unknown) => (typeof v === "string" && f.options?.includes(v) ? t(v) : String(v));

/** A label as the input bar shows it, without the bracket at its end ("(ppm)", "(optional)"; the full-width "（ppm）" of zh and ja as well). The
 *  page's own way of writing "(optional)" goes too, with or without a space before it: the translated suffix ("（可选）", "（任意）",
 *  "(facultatif)") or the English one. */
function barLabel(label: unknown): string {
  const mine = t("(optional)");
  let s = String(label);
  if (mine !== "(optional)") s = s.replace(new RegExp("\\s*" + mine.replace(/[.*+?^${}()|[\]\\]/g, "\\$&") + "$"), "");
  return s.replace(/\s*[(（].*[)）]$/, "");
}

/** Short text for the collapsed input bar: `lead` = the compound / spectrum / precursor, `items` = the other settings.
 *  Defaults the user did not touch still appear (they shaped the result); empty optional fields are left out. */
export function summarizeInputs(fields: Field[], values: Record<string, any>): { lead: string; items: string[] } {
  const items: string[] = [];
  let lead = "";
  const short = (s: string, n = 48) => (s.length > n ? s.slice(0, n - 1) + "…" : s);
  for (const f of visibleFields(fields, values)) {
    if (!f.id) continue;
    const v = values[f.id] ?? initialValues([f])[f.id];
    let text = "";
    switch (f.type) {
      case "structure":
        text = v?.label ? short(String(v.label)) : v?.value ? (v.mode === "Draw" ? t("drawn structure") : short(String(v.value))) : "";
        if (text && !lead) { lead = text; continue; }
        break;
      case "precursor": {
        const adduct = v?.adduct && v.adduct !== "Auto" ? v.adduct : (f.auto_label ?? "");
        text = [v?.precursor_mz ? `m/z ${v.precursor_mz}` : "", v?.ion_mode ? ionLabel(String(v.ion_mode)) : "", adduct].filter(Boolean).join(" · ");
        break;
      }
      case "spectrum": {
        const n = countPeaks(String(v?.value ?? ""));       // peaks as the server reads them, not text lines (see peakCount.ts)
        text = v?.source === "USI" ? short(String(v.value ?? "")) : v?.filename ? short(v.filename, 32) : v?.value ? pair(f.label ?? "", n ? t("{n} peaks", { n }) : t("given")) : "";
        if (text && !lead) { lead = text; continue; }
        break;
      }
      case "ms1": text = v ? t("MS1 isotopes given") : ""; break;
      case "checkbox": text = v ? short(barLabel(f.label), 40) : ""; break;
      case "multiselect": text = Array.isArray(v) && v.length ? pair(f.label ?? "", v.length <= 2 ? v.map((x: string) => optionText(f, x)).join(", ") : t("{n} selected", { n: v.length })) : ""; break;
      case "file": text = v?.name ? short(v.name, 32) : ""; break;
      case "table": text = v ? pair(f.label ?? "", t("given")) : ""; break;
      case "textarea": text = v ? pair(f.label ?? "", t("{n} lines", { n: String(v).split("\n").filter((l: string) => l.trim()).length })) : ""; break;
      default: text = v !== "" && v !== null && v !== undefined ? pair(barLabel(f.label), short(optionText(f, v), 32)) : "";
    }
    if (text) items.push(text);
  }
  return { lead, items };
}
