import { useEffect, useRef, useState, type KeyboardEvent, type PointerEvent } from "react";
import { useT } from "../i18n";

/** Sidebar width limits (px). The width is remembered per browser; double-click the handle to go back to the default. */
export const SIDEBAR_WIDTH = { min: 200, max: 380, initial: 248 };
const KEY = "mzdesk.sidebarWidth";
const clamp = (w: number) => Math.min(SIDEBAR_WIDTH.max, Math.max(SIDEBAR_WIDTH.min, Math.round(w)));

export function useSidebarWidth() {
  const [width, setWidth] = useState(() => {
    try { const v = Number(localStorage.getItem(KEY)); return v ? clamp(v) : SIDEBAR_WIDTH.initial; } catch { return SIDEBAR_WIDTH.initial; }
  });
  useEffect(() => { try { localStorage.setItem(KEY, String(width)); } catch { /* storage unavailable: width just isn't remembered */ } }, [width]);
  return { width, setWidth: (w: number) => setWidth(clamp(w)), reset: () => setWidth(SIDEBAR_WIDTH.initial) };
}

/** Drag handle on the sidebar's right edge. Keyboard: ←/→ (Shift = bigger steps), Home / End = min / max. */
export function SidebarResizer({ width, onChange, onReset }: { width: number; onChange: (w: number) => void; onReset: () => void }) {
  const t = useT();
  const drag = useRef<{ x: number; w: number } | null>(null);
  const [active, setActive] = useState(false);
  const stop = () => { drag.current = null; setActive(false); document.body.classList.remove("resizing"); };
  useEffect(() => stop, []);

  const down = (e: PointerEvent<HTMLDivElement>) => {
    if (e.button !== 0) return;
    e.preventDefault();
    e.currentTarget.setPointerCapture(e.pointerId);
    drag.current = { x: e.clientX, w: width };
    setActive(true);
    document.body.classList.add("resizing");
  };
  const move = (e: PointerEvent<HTMLDivElement>) => { if (drag.current) onChange(drag.current.w + e.clientX - drag.current.x); };
  const key = (e: KeyboardEvent<HTMLDivElement>) => {
    const step = e.shiftKey ? 40 : 10;
    const next = { ArrowLeft: width - step, ArrowRight: width + step, Home: SIDEBAR_WIDTH.min, End: SIDEBAR_WIDTH.max }[e.key];
    if (next !== undefined) { e.preventDefault(); onChange(next); }
  };

  return (
    <div role="separator" aria-orientation="vertical" aria-label={t("Resize sidebar")} tabIndex={0}
      aria-valuemin={SIDEBAR_WIDTH.min} aria-valuemax={SIDEBAR_WIDTH.max} aria-valuenow={width}
      className={"sidebar-resizer" + (active ? " active" : "")} title={t("Drag to resize · double-click to reset")}
      onPointerDown={down} onPointerMove={move} onPointerUp={stop} onPointerCancel={stop} onDoubleClick={onReset} onKeyDown={key} />
  );
}
