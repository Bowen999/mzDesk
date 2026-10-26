import { lazy, Suspense } from "react";
import type { Field } from "../types";

// Ketcher 很大，只在选 Draw 时才下载
const StructureEditor = lazy(() => import("./StructureEditor"));

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
      case "structure": v[f.id] = { mode: "Name / ID / Formula", value: "" }; break;
      case "precursor": v[f.id] = { precursor_mz: "", ion_mode: "Positive", adduct: f.adducts?.Positive?.[0] ?? "Auto" }; break;
      case "spectrum": v[f.id] = { source: "Paste", value: "" }; break;
      case "checkbox": v[f.id] = !!f.value; break;
      case "multiselect": v[f.id] = f.value ?? []; break;
      default: v[f.id] = f.value ?? "";
    }
  }
  return v;
}

type P = { field: Field; value: any; onChange: (v: any) => void };

export function FieldInput({ field, value, onChange }: P) {
  const f = field;
  if (f.id && value === undefined) value = initialValues([f])[f.id];   // 首次渲染或切页时 state 还没初始化
  switch (f.type) {
    case "heading": return <h4 className="field-heading">{f.label}</h4>;
    case "structure": return <StructureInput field={f} value={value} onChange={onChange} />;
    case "precursor": return <PrecursorInput field={f} value={value} onChange={onChange} />;
    case "spectrum": return <SpectrumInput field={f} value={value} onChange={onChange} />;
    case "ms1": return <label className="field wide">{f.label}<textarea rows={4} value={value} onChange={e => onChange(e.target.value)} /></label>;
    case "number": return <label className="field">{f.label}<input type="number" step={f.step ?? "any"} placeholder={f.placeholder ?? ""} value={value ?? ""} onChange={e => onChange(e.target.value)} /></label>;
    case "text": return <label className="field">{f.label}<input type="text" placeholder={f.placeholder ?? ""} value={value ?? ""} onChange={e => onChange(e.target.value)} /></label>;
    case "textarea": return <label className="field wide">{f.label}<textarea rows={f.rows ?? 5} placeholder={f.placeholder ?? ""} value={value ?? ""} onChange={e => onChange(e.target.value)} /></label>;
    case "select": return <label className="field">{f.label}<select value={value} onChange={e => onChange(e.target.value)}>{f.options!.map(o => <option key={o}>{o}</option>)}</select></label>;
    case "radio": return (
      <div className="field"><span>{f.label}</span><div className="radio-row">
        {f.options!.map(o => <label key={o}><input type="radio" checked={value === o} onChange={() => onChange(o)} /> {o}</label>)}
      </div></div>);
    case "multiselect": return (
      <div className="field wide"><span>{f.label}</span><div className="chip-row">
        {f.options!.map(o => {
          const on = (value as string[]).includes(o);
          return <button type="button" key={o} className={"chip" + (on ? " on" : "")}
            onClick={() => onChange(on ? value.filter((x: string) => x !== o) : [...value, o])}>{o}</button>;
        })}
      </div></div>);
    case "checkbox": return <label className="field check"><input type="checkbox" checked={!!value} onChange={e => onChange(e.target.checked)} /> {f.label}</label>;
    case "slider": return <label className="field"><span className="slider-label">{f.label}<b>{value}</b></span><input type="range" min={f.min} max={f.max} step={f.step ?? 1} value={value} onChange={e => onChange(Number(e.target.value))} /></label>;
    case "file": return <label className="field">{f.label}<input type="file" accept={f.accept} onChange={async e => onChange(await readFiles(e.target.files))} /></label>;
    case "table": return (
      <label className="field wide">{f.label}
        <div className="grid-head">{f.columns!.join("  |  ")}</div>
        <textarea rows={5} placeholder="one row per line, columns separated by space, comma or tab" value={value} onChange={e => onChange(e.target.value)} />
        {f.caption && <small>{f.caption}</small>}
      </label>);
    default: return <div className="field">unknown field type {f.type}</div>;
  }
}

function StructureInput({ field, value, onChange }: P) {
  const modes = ["Name / ID / Formula", "Draw", "SMILES / InChI"];
  const set = (patch: Partial<typeof value>) => onChange({ ...value, ...patch });
  // 画布与 SMILES 框共用同一个 SMILES：Draw <-> SMILES 互转时保留，进出 Name 模式时清空
  const switchMode = (m: string) => {
    const keep = (m === "Draw" || m === "SMILES / InChI") && (value.mode === "Draw" || value.mode === "SMILES / InChI");
    onChange({ ...value, mode: m, value: keep ? value.value : "", molfile: "" });
  };
  return (
    <div className="field wide composite">
      <span>{field.label}{field.optional ? " (optional)" : ""}</span>
      <div className="radio-row">{modes.map(m => <label key={m}><input type="radio" checked={value.mode === m} onChange={() => switchMode(m)} /> {m}</label>)}</div>
      {value.mode === "Draw" ? (
        <>
          <Suspense fallback={<div className="placeholder-box editor">Loading structure editor…</div>}>
            <StructureEditor smiles={value.value} molfile={value.molfile} onChange={v => set(v)} />
          </Suspense>
          <input type="text" className="mono" readOnly placeholder="SMILES from editor" value={value.value} />
        </>
      ) : value.mode === "SMILES / InChI" ? (
        <input type="text" className="mono" placeholder="CN1C=NC2=C1C(=O)N(C(=O)N2C)C" value={value.value} onChange={e => set({ value: e.target.value })} />
      ) : (
        <input type="text" placeholder="Name, CAS, formula, InChIKey, PubChem CID or HMDB / KEGG / ChEBI / LIPID MAPS / DrugBank ID — e.g. caffeine" value={value.value} onChange={e => set({ value: e.target.value })} />
      )}
    </div>
  );
}

function PrecursorInput({ field, value, onChange }: P) {
  const set = (patch: Partial<typeof value>) => onChange({ ...value, ...patch });
  const lists = field.adducts ?? { Positive: ["Auto"], Negative: ["Auto"] };
  const adducts = lists[value.ion_mode] ?? lists.Positive;
  return (
    <div className="field wide composite row3">
      {field.mz !== false && <label className="field">{field.label ?? "Precursor m/z"}<input type="number" step="0.0001" placeholder="195.0877" value={value.precursor_mz} onChange={e => set({ precursor_mz: e.target.value })} /></label>}
      <div className="field"><span>Ion mode</span><div className="radio-row">
        {["Positive", "Negative"].map(m => <label key={m}><input type="radio" checked={value.ion_mode === m} onChange={() => set({ ion_mode: m, adduct: lists[m]?.[0] ?? "Auto" })} /> {m}</label>)}
      </div></div>
      {field.adduct !== false && <label className="field">Adduct<select value={value.adduct} onChange={e => set({ adduct: e.target.value })}>{adducts.map(a => <option key={a}>{a}</option>)}</select></label>}
    </div>
  );
}

function SpectrumInput({ field, value, onChange }: P) {
  const set = (patch: Partial<typeof value>) => onChange({ ...value, ...patch });
  const sources = ["Paste", "MGF / MSP file", "USI"];
  return (
    <div className="field wide composite">
      <span>{field.label}</span>
      <div className="radio-row">{sources.map(s => <label key={s}><input type="radio" checked={value.source === s} onChange={() => set({ source: s, value: "", filename: "" })} /> {s}</label>)}</div>
      {value.source === "Paste" && <textarea rows={5} placeholder="m/z intensity per line (or mz:int, MSP, MGF, JSON, MassBank, HMDB XML …)" value={value.value} onChange={e => set({ value: e.target.value })} />}
      {value.source === "MGF / MSP file" && (
        <div className="file-row">
          <input type="file" accept=".mgf,.msp,.txt,.csv,.json,.xml,.ms" multiple={!!field.multi}
            onChange={async e => { const r = await readFiles(e.target.files); set({ value: r.text, filename: r.name }); }} />
          {value.filename && <small>{value.filename}</small>}
        </div>)}
      {value.source === "USI" && <input type="text" placeholder="mzspec:GNPS:..." value={value.value} onChange={e => set({ value: e.target.value })} />}
    </div>
  );
}

/** Short text for the collapsed input bar: `lead` = the compound / spectrum / precursor, `items` = the other settings.
 *  Defaults the user did not touch still appear (they shaped the result); empty optional fields are left out. */
export function summarizeInputs(fields: Field[], values: Record<string, any>): { lead: string; items: string[] } {
  const items: string[] = [];
  let lead = "";
  const short = (s: string, n = 48) => (s.length > n ? s.slice(0, n - 1) + "…" : s);
  for (const f of fields) {
    if (!f.id) continue;
    const v = values[f.id] ?? initialValues([f])[f.id];
    let text = "";
    switch (f.type) {
      case "structure":
        text = v?.value ? (v.mode === "Draw" ? "drawn structure" : short(String(v.value))) : "";
        if (text && !lead) { lead = text; continue; }
        break;
      case "precursor":
        text = [v?.precursor_mz ? `m/z ${v.precursor_mz}` : "", v?.ion_mode, v?.adduct && v.adduct !== "Auto" ? v.adduct : ""].filter(Boolean).join(" · ");
        break;
      case "spectrum": {
        const n = String(v?.value ?? "").split(/\n|(?<=\d)\s+(?=\d+(?:\.\d+)?[:\s])/).filter(l => l.trim()).length;
        text = v?.source === "USI" ? short(String(v.value ?? "")) : v?.filename ? short(v.filename, 32) : v?.value ? `${f.label}: ${n} lines` : "";
        if (text && !lead) { lead = text; continue; }
        break;
      }
      case "ms1": text = v ? "MS1 isotopes given" : ""; break;
      case "checkbox": text = v ? short(String(f.label).replace(/\s*\(.*\)$/, ""), 40) : ""; break;
      case "multiselect": text = Array.isArray(v) && v.length ? `${f.label}: ${v.length <= 2 ? v.join(", ") : v.length + " selected"}` : ""; break;
      case "file": text = v?.name ? short(v.name, 32) : ""; break;
      case "table": text = v ? `${f.label}: given` : ""; break;
      case "textarea": text = v ? `${f.label}: ${String(v).split("\n").filter((l: string) => l.trim()).length} lines` : ""; break;
      default: text = v !== "" && v !== null && v !== undefined ? `${String(f.label).replace(/\s*\(.*\)$/, "")}: ${short(String(v), 32)}` : "";
    }
    if (text) items.push(text);
  }
  return { lead, items };
}
