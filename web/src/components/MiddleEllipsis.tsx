import { useEffect, useRef, useState } from "react";

/**
 * Text that keeps its head and tail and drops the middle when its box is too narrow:
 * "Caffeine and theoph…eobromine". The full text stays in the title attribute.
 *
 * `maxWidth` is the character budget for boxes whose own width follows the text (the compound
 * tabs of the formats block); everywhere else the room is the parent's width minus the parent's
 * other children, so the compound header can keep its formula on the same line.
 */
export function MiddleEllipsis({ text, maxWidth, className, title }: {
  text: string; maxWidth?: number; className?: string; title?: string;
}) {
  const ref = useRef<HTMLSpanElement>(null);
  const [shown, setShown] = useState(text);

  useEffect(() => {
    const el = ref.current;
    if (!el) return;
    const ctx = document.createElement("canvas").getContext("2d");
    const room = () => {
      if (maxWidth) return maxWidth;
      const parent = el.parentElement;
      if (!parent) return el.clientWidth;
      const others = Array.from(parent.children).filter(c => c !== el) as HTMLElement[];
      const gap = parseFloat(getComputedStyle(parent).columnGap || "0") || 0;
      return Math.max(0, parent.clientWidth - others.reduce((w, c) => w + c.offsetWidth, 0) - gap * others.length);
    };
    const fit = () => {
      if (!ctx) { setShown(text); return; }
      const s = getComputedStyle(el);
      ctx.font = `${s.fontStyle} ${s.fontWeight} ${s.fontSize} ${s.fontFamily}`;
      const w = room();
      if (!w || ctx.measureText(text).width <= w) { setShown(text); return; }
      // keep 60 % of the kept characters from the head, so both ends stay readable
      const cut = (keep: number) => {
        const head = Math.ceil(keep * 0.6), tail = keep - head;
        return `${text.slice(0, head)}…${tail > 0 ? text.slice(-tail) : ""}`;
      };
      let keep = text.length - 1;
      while (keep > 2 && ctx.measureText(cut(keep)).width > w) keep -= 1;
      setShown(keep > 2 ? cut(keep) : "…");
    };
    fit();
    const ro = typeof ResizeObserver !== "undefined" && !maxWidth ? new ResizeObserver(fit) : null;
    if (ro) ro.observe(el.parentElement ?? el);
    return () => ro?.disconnect();
  }, [text, maxWidth]);

  return <span ref={ref} className={"mid-ellipsis" + (className ? " " + className : "")} title={title ?? text}>{shown}</span>;
}
