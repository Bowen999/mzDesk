import { useEffect, useState } from "react";
import { useT } from "../i18n";
import { useReportMode } from "../report/mode";
import type { Block, BlockData } from "../types";
import { HelpButton, helpLines } from "./Help";
import { MiddleEllipsis } from "./MiddleEllipsis";

/** The formats block (Outputs.tsx) and its Copy button; the Selected spectrum block of Reference Spectra (ReferenceSpectrum.tsx) shows the same block. */
export function CopyButton({ text }: { text: string }) {
  const t = useT();
  const [done, setDone] = useState(false);
  return <button type="button" className="copy" onClick={() => navigator.clipboard?.writeText(text).then(() => { setDone(true); setTimeout(() => setDone(false), 1200); })}>{done ? t("Copied") : t("Copy")}</button>;
}

/** 同一内容的多种格式：一排小标签 + 一个代码框；记住上次选的标签。
 *  Several compounds (Structure & formats) add a compound tab row; the format tab stays selected when switching. */
export function FormatsBlock({ block, data, empty }: { block: Block; data?: BlockData; empty: string }) {
  const t = useT();
  const report = useReportMode();
  const [tab, setTab] = useState("");
  const [setIndex, setSetIndex] = useState(0);
  useEffect(() => setSetIndex(0), [data]);
  const sets = data?.sets ?? [];
  const active = sets.length ? (sets[Math.min(setIndex, sets.length - 1)] ?? sets[0]) : null;
  const formats = active?.formats ?? data?.formats ?? {};
  const names = Object.keys(formats);
  const cur = names.includes(tab) ? tab : names[0];
  const text = cur ? formats[cur] : undefined;
  const title = block.title ?? block.label;
  // the report has no tabs: every format of every compound, one after the other
  if (report) return (
    <div className="out">
      <h4 className="out-title">{title}</h4>
      {(sets.length ? sets : [{ label: "", formats }]).map((set, i) => (
        <div key={i}>
          {sets.length > 1 && <div className="fig-label">{set.label}</div>}
          {Object.entries(set.formats ?? {}).map(([name, value]) => <div key={name}><div className="fig-label">{name}</div><pre className="code">{value}</pre></div>)}
        </div>))}
    </div>);
  return (
    <div className="out">
      <h4 className="out-title">{title}<HelpButton lines={helpLines(block.help, data?.help)} label={t("Help: {name}", { name: title ?? "" })} />{text !== undefined && <CopyButton text={text} />}</h4>
      {sets.length > 1 && (
        <div className="tabs compound-tabs" role="tablist" aria-label={t("Compound")}>
          {sets.map((s, i) => <button key={`${s.label}-${i}`} type="button" role="tab" aria-selected={i === Math.min(setIndex, sets.length - 1)}
            className={"tab" + (i === Math.min(setIndex, sets.length - 1) ? " on" : "")} onClick={() => setSetIndex(i)}>
            <MiddleEllipsis text={s.label} maxWidth={216} /></button>)}
        </div>)}
      {names.length > 0 && (
        <div className="tabs" role="tablist">
          {names.map(n => <button key={n} type="button" role="tab" aria-selected={n === cur} className={"tab" + (n === cur ? " on" : "")} onClick={() => setTab(n)}>{n}</button>)}
        </div>)}
      <pre className="code">{text ?? `— ${block.label} (${empty}) —`}</pre>
    </div>);
}
