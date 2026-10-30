import type { CompareCompound } from "../types";
import { Mol3D, ThemedSvg } from "./Figure";

/**
 * Structure & formats: one card per compound (1-3), for side-by-side comparison.
 * Both boxes are square and take the full card width, so the 2D depiction and the 3D viewer
 * always have the same width and height. One compound: the two boxes sit side by side.
 */
export function CompoundCompare({ items }: { items: CompareCompound[] }) {
  return (
    <div className={"compare-grid" + (items.length === 1 ? " single" : "")}>
      {items.map((c, i) => (
        <div className="compare-col" key={`${c.label}-${i}`}>
          <div className="compare-head">
            <span className="compare-name">{c.label || `Compound ${i + 1}`}</span>
            {c.formula && <span className="compare-formula">{c.formula}</span>}
          </div>
          <div className="compare-views">
            <div className="compare-cell">
              <div className="compare-label">2D depiction</div>
              <div className="compare-box">
                {c.svg
                  ? <ThemedSvg svg={c.svg} className="compare-svg" alt={`${c.label || `Compound ${i + 1}`} 2D depiction`} />
                  : <span className="compare-empty">No 2D depiction</span>}
              </div>
            </div>
            <div className="compare-cell">
              <div className="compare-label">3D viewer</div>
              <div className="compare-box">
                {c.mol3d ? <Mol3D molblock={c.mol3d} /> : <span className="compare-empty">3D viewer unavailable</span>}
              </div>
            </div>
          </div>
        </div>
      ))}
    </div>
  );
}
