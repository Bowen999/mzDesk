import type { Field } from "../types";

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
      case "structure": v[f.id] = { mode: "SMILES / InChI", value: "" }; break;
      case "precursor": v[f.id] = { precursor_mz: "", ion_mode: "Positive", adduct: "Auto" }; break;
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
    case "slider": return <label className="field">{f.label}: <b>{value}</b><input type="range" min={f.min} max={f.max} step={f.step ?? 1} value={value} onChange={e => onChange(Number(e.target.value))} /></label>;
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
  const modes = ["SMILES / InChI", "Draw", "Name / ID / Formula"];
  const set = (patch: Partial<typeof value>) => onChange({ ...value, ...patch });
  return (
    <div className="field wide composite">
      <span>{field.label}{field.optional ? " (optional)" : ""}</span>
      <div className="radio-row">{modes.map(m => <label key={m}><input type="radio" checked={value.mode === m} onChange={() => set({ mode: m })} /> {m}</label>)}</div>
      {value.mode === "Draw" ? (
        <>
          <div className="placeholder-box editor">[structure editor placeholder] Ketcher goes here → SMILES appears below</div>
          <input type="text" disabled placeholder="SMILES from editor" value={value.value} />
        </>
      ) : value.mode === "SMILES / InChI" ? (
        <input type="text" placeholder="CN1C=NC2=C1C(=O)N(C(=O)N2C)C" value={value.value} onChange={e => set({ value: e.target.value })} />
      ) : (
        <input type="text" placeholder="caffeine · HMDB0001847 · C8H10N4O2" value={value.value} onChange={e => set({ value: e.target.value })} />
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
      <label className="field">{field.label ?? "Precursor m/z"}<input type="number" step="0.0001" placeholder="195.0877" value={value.precursor_mz} onChange={e => set({ precursor_mz: e.target.value })} /></label>
      <div className="field"><span>Ion mode</span><div className="radio-row">
        {["Positive", "Negative"].map(m => <label key={m}><input type="radio" checked={value.ion_mode === m} onChange={() => set({ ion_mode: m, adduct: "Auto" })} /> {m}</label>)}
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
