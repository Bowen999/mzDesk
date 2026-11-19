import { createContext, useContext } from "react";

/**
 * Report mode. The Session report (report/build.ts) draws every page once more, off screen, inside this context; nothing on screen sets it.
 * The same components then show their whole content as a static page: folded blocks are native <details>, a table lists all its rows (up to
 * `maxRows`), long text is shown in full, and controls that need the app (pages, filter, tabs, "?", View, Copy) are left out.
 * Outside the context useReportMode() is null and every component draws exactly what the screen shows.
 */
export type ReportMode = {
  maxRows: number;                           // a table shows at most this many rows, then says the rest is in report.json
  track: (done: Promise<unknown>) => void;   // a figure that draws late (Plotly, 3Dmol) hands its promise in; the report waits for all of them
};

/** Rows of one table in report.html. report.json always has them all. */
export const MAX_REPORT_ROWS = 500;

export const ReportModeContext = createContext<ReportMode | null>(null);
export const useReportMode = () => useContext(ReportModeContext);
