import { useMemo, useState } from "react";
import { useT } from "../i18n";
import { useReportMode } from "../report/mode";
import type { FoundGroup } from "../types";
import { ThemedSvg } from "./Figure";
import "./groupmap.css";

/**
 * Functional groups on the input structure (spec block "groupmap", data {"svg", "found"}; contract in mzdesk/spec.py).
 * Left: the groups found, under their family. Right: the input structure with the chosen group tinted, the group
 * alone, its sites (all together, or one at a time with atom numbers) and the SMARTS that was counted.
 */
export function GroupMap({ svg, found }: { svg: string; found: FoundGroup[] }) {
  const t = useT();
  const report = useReportMode();       // the report draws one panel per group, with no site buttons
  const [sel, setSel] = useState(0);        // index into found; -1 = plain input structure
  const [site, setSite] = useState(-1);     // -1 = all sites together
  const families = useMemo(() => {
    const out: { family: string; items: { g: FoundGroup; i: number }[] }[] = [];
    found.forEach((g, i) => {
      let f = out.find(x => x.family === g.family);
      if (!f) out.push(f = { family: g.family, items: [] });
      f.items.push({ g, i });
    });
    return out;
  }, [found]);
  const cur = found[sel] as FoundGroup | undefined;
  const choose = (i: number) => { setSite(-1); setSel(i === sel ? -1 : i); };
  // one site: its picture has the atom numbers; several: the all-sites picture, or the chosen site
  const shown = !cur ? svg : site >= 0 && cur.sites[site] ? cur.sites[site].svg : cur.all_svg ?? cur.sites[0].svg;
  const atoms = !cur ? "" : site >= 0 && cur.sites[site] ? cur.sites[site].atoms : cur.count === 1 ? cur.sites[0].atoms : "";
  return (
    <div className="gm">
      <div className="gm-list" role="group" aria-label={t("Functional groups found")}>
        {families.map(f => (
          <div className="gm-family" key={f.family}>
            <div className="gm-family-name">{f.family}</div>
            {f.items.map(({ g, i }) => {
              const row = <>
                <ThemedSvg svg={g.svg} className="gm-thumb" alt={t("Drawing of {name}", { name: g.name })} />
                <span className="gm-name">{g.name}</span>
                <span className="gm-count" title={g.count === 1 ? t("1 site") : t("{n} sites", { n: g.count })}>{g.count}</span>
              </>;
              return report
                ? <div key={g.name} className={"gm-row" + (i === sel ? " on" : "")}>{row}</div>
                : <button key={g.name} type="button" className={"gm-row" + (i === sel ? " on" : "")} aria-pressed={i === sel} onClick={() => choose(i)}>{row}</button>;
            })}
          </div>))}
      </div>
      <div className="gm-view">
        <div className="gm-head">
          {cur
            ? <><span className="gm-title">{cur.name}</span><span className="gm-sub">{cur.count === 1 ? t("1 site") : t("{n} sites", { n: cur.count })} · {cur.family}</span></>
            : <span className="gm-title">{t("Input structure")}</span>}
        </div>
        {cur && cur.count > 1 && !report && (
          <div className="gm-sites" role="group" aria-label={t("Sites")}>
            <button type="button" className={"gm-site" + (site < 0 ? " on" : "")} aria-pressed={site < 0} onClick={() => setSite(-1)}>{t("All sites")}</button>
            {cur.sites.map((s, i) => (
              <button key={i} type="button" className={"gm-site" + (site === i ? " on" : "")} aria-pressed={site === i} title={t("Atoms {atoms}", { atoms: s.atoms })} onClick={() => setSite(i)}>{t("Site {n}", { n: i + 1 })}</button>))}
            {cur.more ? <span className="gm-more">{t("+{n} more, not drawn", { n: cur.more })}</span> : null}
          </div>)}
        <ThemedSvg svg={shown} className="gm-structure" alt={cur ? t("{name} on the input structure", { name: cur.name }) : t("Input structure")} />
        {cur && (
          <dl className="gm-facts">
            {atoms && <><dt>{t("Atoms")}</dt><dd>{atoms} <span className="gm-note">{t("(numbered on the structure)")}</span></dd></>}
            {!atoms && cur.count > 1 && <><dt>{t("Atoms")}</dt><dd>{cur.atoms}{cur.more ? "" : ""} <span className="gm-note">{t("(one set per site; pick a site to see the numbers)")}</span></dd></>}
            <dt>SMARTS</dt><dd><code>{cur.smarts}</code></dd>
          </dl>)}
        {!cur && <p className="gm-note">{t("Choose a group to see where it sits.")}</p>}
      </div>
    </div>);
}
