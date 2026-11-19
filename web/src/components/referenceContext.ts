import { createContext, useContext } from "react";
import type { ReferenceTarget } from "../types";

/** Reference Spectra: which library record the Selected spectrum block (ReferenceSpectrum.tsx) shows. OutputList holds it; the table's rows with a
 *  `reference` set it on a click (Blocks.tsx) and the one shown is highlighted. A separate file, so that Blocks.tsx and ReferenceSpectrum.tsx (which
 *  draws a table with Blocks.tsx) do not import each other.
 *  current: the record shown (the run's own until a click); picked: counts the clicks, so a click scrolls the block into view even on the same row.
 *  null = nothing to select (no Selected spectrum block on the page, or the Session report, which is a file). */
export type ReferenceSelection = { current: ReferenceTarget; select: (target: ReferenceTarget) => void; picked: number };
export const ReferenceSelectionContext = createContext<ReferenceSelection | null>(null);
export const useReferenceSelection = () => useContext(ReferenceSelectionContext);
export const sameTarget = (a?: ReferenceTarget | null, b?: ReferenceTarget | null) =>
  !!a && !!b && a.library === b.library && a.row === b.row && a.accession === b.accession;
