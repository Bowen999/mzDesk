import { createPortal } from "react-dom";
import { useEffect, useId, useLayoutEffect, useRef, useState, type ReactNode, type SyntheticEvent } from "react";
import type { Help } from "../types";

/**
 * Help "?" popovers. Contract (any block, and any item inside cards / chips / portal tiles / table rows):
 *   help = string | string[]    static text from the spec (block.help) and/or runtime text (data.help, item.help); all sources are merged
 *   One line renders as a paragraph, several lines as bullets. Text may use <sub> <sup> <i> <b> (nothing else is interpreted).
 * Behaviour: hover or keyboard focus shows the popover; click pins it until a click outside or Esc; only one popover stays pinned.
 * Accessibility: the trigger is a <button> with aria-expanded and (while open) aria-describedby -> the popover (role="tooltip").
 * The popover is portalled to <body> with fixed positioning, so scrolling table wrappers and collapsed groups never clip it.
 */

/** Merge help sources into a flat list of non-empty strings. */
export function helpLines(...parts: (Help | null | undefined)[]): string[] {
  const out: string[] = [];
  for (const p of parts) for (const s of Array.isArray(p) ? p : [p]) if (typeof s === "string" && s.trim()) out.push(s.trim());
  return out;
}

/** Text with a small safe markup subset (<sub>, <sup>, <i>, <b>); everything else, including any other tag, stays literal text. */
export function RichText({ text }: { text: string }) {
  const re = /<(sub|sup|i|b)>([\s\S]*?)<\/\1>/g;
  const parts: ReactNode[] = [];
  let last = 0, k = 0, m: RegExpExecArray | null;
  while ((m = re.exec(text))) {
    if (m.index > last) parts.push(text.slice(last, m.index));
    const inner = <RichText text={m[2]} />;
    parts.push(m[1] === "sub" ? <sub key={k++}>{inner}</sub> : m[1] === "sup" ? <sup key={k++}>{inner}</sup>
      : m[1] === "i" ? <i key={k++}>{inner}</i> : <b key={k++}>{inner}</b>);
    last = m.index + m[0].length;
  }
  if (last === 0) return <>{text}</>;
  if (last < text.length) parts.push(text.slice(last));
  return <>{parts}</>;
}

let closePinned: (() => void) | null = null;   // the one popover that is currently pinned

/** State + handlers for one popover. Spread `triggerProps` on a <button>; render `popover` anywhere next to it. */
export function useHelpPopover(lines: string[]) {
  const id = useId();
  const trigger = useRef<HTMLButtonElement>(null);
  const pop = useRef<HTMLDivElement>(null);
  const timer = useRef<number | undefined>(undefined);
  const [hover, setHover] = useState(false);
  const [focus, setFocus] = useState(false);
  const [pinned, setPinned] = useState(false);
  const [pos, setPos] = useState<{ left: number; top: number; maxH: number } | null>(null);
  const visible = lines.length > 0 && (hover || focus || pinned);

  useEffect(() => () => window.clearTimeout(timer.current), []);

  // only one pinned popover at a time
  useEffect(() => {
    if (!pinned) return;
    const mine = () => { setPinned(false); setHover(false); };
    if (closePinned && closePinned !== mine) closePinned();
    closePinned = mine;
    return () => { if (closePinned === mine) closePinned = null; };
  }, [pinned]);

  // click outside / Esc close it (also when it was only hovered or focused)
  useEffect(() => {
    if (!visible) return;
    const down = (e: PointerEvent) => {
      const t = e.target as Node;
      if (trigger.current?.contains(t) || pop.current?.contains(t)) return;
      setPinned(false); setHover(false);
    };
    const key = (e: KeyboardEvent) => { if (e.key === "Escape") { setPinned(false); setHover(false); setFocus(false); } };
    document.addEventListener("pointerdown", down);
    document.addEventListener("keydown", key);
    return () => { document.removeEventListener("pointerdown", down); document.removeEventListener("keydown", key); };
  }, [visible]);

  // place it under the trigger (above if there is no room), kept inside the viewport; follow scroll / resize while open
  useLayoutEffect(() => {
    if (!visible) { setPos(null); return; }
    const place = () => {
      const t = trigger.current, p = pop.current;
      if (!t || !p) return;
      const r = t.getBoundingClientRect(), w = p.offsetWidth, M = 8, GAP = 6;
      const cap = window.innerHeight * 0.6;                       // CSS max-height 60vh: long help scrolls inside the popover
      const natural = Math.min(p.scrollHeight + 2, cap);          // height it wants (content + 1px borders), capped
      const below = window.innerHeight - r.bottom - GAP - M, above = r.top - GAP - M;
      const up = natural > below && above > below;                // flip above only when it does not fit below and there is more room above
      const maxH = Math.min(cap, Math.max(up ? above : below, 120));
      const h = Math.min(natural, maxH);
      const left = Math.max(M, Math.min(r.left + r.width / 2 - w / 2, window.innerWidth - w - M));
      const top = up ? r.top - GAP - h : r.bottom + GAP;
      setPos(prev => (prev && prev.left === left && prev.top === top && prev.maxH === maxH ? prev : { left, top, maxH }));
    };
    place();
    window.addEventListener("resize", place);
    window.addEventListener("scroll", place, true);
    return () => { window.removeEventListener("resize", place); window.removeEventListener("scroll", place, true); };
  }, [visible]);

  const enter = () => { window.clearTimeout(timer.current); setHover(true); };
  const leave = () => { window.clearTimeout(timer.current); timer.current = window.setTimeout(() => setHover(false), 120); };

  const triggerProps = {
    ref: trigger,
    "aria-expanded": visible,
    "aria-describedby": visible ? id : undefined,
    onMouseEnter: enter,
    onMouseLeave: leave,
    // mouse clicks also focus a button; only keyboard focus (:focus-visible) should open the popover by itself
    onFocus: (e: SyntheticEvent<HTMLButtonElement>) => { if (e.currentTarget.matches(":focus-visible")) setFocus(true); },
    onBlur: () => setFocus(false),
    onClick: (e: SyntheticEvent) => { e.stopPropagation(); setPinned(p => !p); },
  };

  const popover = visible ? createPortal(
    <div ref={pop} id={id} role="tooltip" className={"hp-pop" + (pinned ? " pinned" : "")}
      style={{ left: pos?.left ?? 0, top: pos?.top ?? 0, maxHeight: pos?.maxH, visibility: pos ? "visible" : "hidden" }}
      onMouseEnter={enter} onMouseLeave={leave}>
      {lines.length === 1
        ? <p><RichText text={lines[0]} /></p>
        : <ul>{lines.map((l, i) => <li key={i}><RichText text={l} /></li>)}</ul>}
    </div>, document.body) : null;

  return { triggerProps, popover };
}

/** The small "?" icon button. Renders nothing when there is no help text. */
export function HelpButton({ lines, label = "More information", className = "" }: { lines: string[]; label?: string; className?: string }) {
  const { triggerProps, popover } = useHelpPopover(lines);
  if (!lines.length) return null;
  return (
    <>
      <button type="button" className={"hp-btn" + (className ? " " + className : "")} aria-label={label} {...triggerProps}>?</button>
      {popover}
    </>
  );
}
