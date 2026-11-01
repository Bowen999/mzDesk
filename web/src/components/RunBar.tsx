import { useEffect, useRef, useState } from "react";
import { waitingPercent } from "../progressCurve";

type State = { pct: number; phase: "idle" | "run" | "done" };

/** The bar under a module's inputs while it runs. The server reports no progress for one module, so the bar is eased on time alone:
 *  it starts at once, slows down, and never reaches 100 by itself. 100 means the result is in; the bar then fades out. */
export function RunBar({ busy }: { busy: boolean }) {
  const [state, setState] = useState<State>({ pct: 0, phase: "idle" });
  const t0 = useRef(0);
  useEffect(() => {
    if (busy) {
      t0.current = performance.now();
      setState({ pct: waitingPercent(0), phase: "run" });
      const id = window.setInterval(
        () => setState(s => (s.phase === "run" ? { pct: waitingPercent((performance.now() - t0.current) / 1000), phase: "run" } : s)), 120);
      return () => window.clearInterval(id);
    }
    setState(s => (s.phase === "run" ? { pct: 100, phase: "done" } : s));
    const id = window.setTimeout(() => setState(s => (s.phase === "done" ? { pct: 0, phase: "idle" } : s)), 700);
    return () => window.clearTimeout(id);
  }, [busy]);
  if (state.phase === "idle") return null;
  return (
    <div className={"run-bar" + (state.phase === "done" ? " is-done" : "")} role="progressbar" aria-label="Running"
      aria-valuemin={0} aria-valuemax={100} aria-valuenow={Math.round(state.pct)}>
      <span style={{ width: `${state.pct}%` }} />
    </div>
  );
}
