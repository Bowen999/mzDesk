import { createPortal } from "react-dom";
import { useCallback, useEffect, useRef, useState } from "react";
import { useT } from "../i18n";
import { isInternetOn } from "../network";

/**
 * Confirmation for a run that needs the internet while "Internet search" is off.
 * Rendered by useNetworkGuard(); Esc, the backdrop and Cancel answer "no", Continue answers "yes".
 * Only the use that triggered it is confirmed: the switch stays off and the next run asks again.
 */

function GlobeSlash() {
  return (
    <svg viewBox="0 0 24 24" width="20" height="20" fill="none" stroke="currentColor" strokeWidth="1.7"
      strokeLinecap="round" strokeLinejoin="round" aria-hidden="true" focusable="false">
      <circle cx="12" cy="12" r="9" />
      <path d="M3 12h18" />
      <path d="M12 3a15 15 0 0 1 0 18M12 3a15 15 0 0 0 0 18" />
      <path d="M4 20 20 4" />
    </svg>
  );
}

function NetworkDialog({ uses, onCancel, onContinue }: { uses: string[]; onCancel: () => void; onContinue: () => void }) {
  const t = useT();
  const card = useRef<HTMLDivElement>(null);
  const continueButton = useRef<HTMLButtonElement>(null);

  // Focus the primary action, keep Tab inside the card, close on Esc, lock the page scroll, restore focus on close.
  useEffect(() => {
    const previous = document.activeElement as HTMLElement | null;
    const overflow = document.body.style.overflow;
    document.body.style.overflow = "hidden";
    const timer = window.setTimeout(() => continueButton.current?.focus(), 30);
    const key = (e: KeyboardEvent) => {
      if (e.key === "Escape") { e.preventDefault(); onCancel(); return; }
      if (e.key !== "Tab" || !card.current) return;
      const nodes = Array.from(card.current.querySelectorAll<HTMLElement>("button:not([disabled])")).filter(n => n.offsetParent !== null);
      if (!nodes.length) return;
      const first = nodes[0], last = nodes[nodes.length - 1], active = document.activeElement;
      if (e.shiftKey && (active === first || !card.current.contains(active))) { e.preventDefault(); last.focus(); }
      else if (!e.shiftKey && active === last) { e.preventDefault(); first.focus(); }
    };
    document.addEventListener("keydown", key, true);
    return () => {
      window.clearTimeout(timer);
      document.removeEventListener("keydown", key, true);
      document.body.style.overflow = overflow;
      previous?.focus?.();
    };
  }, [onCancel]);

  return createPortal(
    <div className="dialog-overlay" onMouseDown={e => { if (e.target === e.currentTarget) onCancel(); }}>
      <div className="dialog-card" role="dialog" aria-modal="true" aria-labelledby="net-dialog-title" aria-describedby="net-dialog-desc" ref={card}>
        <div className="dialog-badge"><GlobeSlash /></div>
        <h2 id="net-dialog-title">{t("Internet search is off")}</h2>
        <p id="net-dialog-desc">{t("This run needs the internet:")}</p>
        <ul className="dialog-list">{uses.map((u, i) => <li key={i}>{u}</li>)}</ul>
        <p className="dialog-note">{t("This run only. The switch stays off.")}</p>
        <div className="dialog-actions">
          <button type="button" className="btn secondary" onClick={onCancel}>{t("Cancel")}</button>
          <button type="button" className="btn" ref={continueButton} onClick={onContinue}>{t("Continue")}</button>
        </div>
      </div>
    </div>,
    document.body,
  );
}

type Pending = { uses: string[]; resolve: (ok: boolean) => void };

/**
 * Ask before an online run while the switch is off. `confirm(uses)` resolves true (run) or false (cancelled);
 * a dialog is already open -> false. Render `dialog` next to the page. Unmounting answers "no".
 */
export function useNetworkGuard() {
  const [pending, setPending] = useState<Pending | null>(null);
  const pendingRef = useRef<Pending | null>(null);
  const answer = useCallback((ok: boolean) => {
    const p = pendingRef.current;
    pendingRef.current = null;
    setPending(null);
    p?.resolve(ok);
  }, []);
  useEffect(() => () => { pendingRef.current?.resolve(false); pendingRef.current = null; }, []);

  const confirm = useCallback((uses: string[]) => {
    if (!uses.length || isInternetOn()) return Promise.resolve(true);
    if (pendingRef.current) return Promise.resolve(false);
    return new Promise<boolean>(resolve => {
      pendingRef.current = { uses, resolve };
      setPending(pendingRef.current);
    });
  }, []);

  const cancel = useCallback(() => answer(false), [answer]);
  const proceed = useCallback(() => answer(true), [answer]);
  const dialog = pending
    ? <NetworkDialog uses={pending.uses} onCancel={cancel} onContinue={proceed} />
    : null;
  return { confirm, cancel, dialog };
}
