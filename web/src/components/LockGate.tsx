import { useEffect, useLayoutEffect, useRef, useState, type FormEvent, type ReactNode } from "react";
import { useLocation } from "react-router-dom";
import { LockUnavailable, unlock, useUnlocked } from "../lock";
import "../lock.css";

/** Small padlock for the sidebar and the unlock card. Colour follows the surrounding text. */
export function LockIcon({ size = 12, open = false, className = "lock-icon" }: { size?: number; open?: boolean; className?: string }) {
  return (
    <svg className={className} width={size} height={size} viewBox="0 0 16 16" fill="none" stroke="currentColor"
      strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true" focusable="false">
      <rect x="3" y="7" width="10" height="7" rx="1.5" />
      <path d={open ? "M5.5 7V4.5a2.5 2.5 0 0 1 4.9-.7" : "M5.5 7V4.5a2.5 2.5 0 0 1 5 0V7"} />
    </svg>
  );
}

function UnlockCard() {
  const [pw, setPw] = useState("");
  const [err, setErr] = useState("");
  const [busy, setBusy] = useState(false);
  const input = useRef<HTMLInputElement>(null);
  useEffect(() => { input.current?.focus(); }, []);

  const submit = async (e: FormEvent) => {
    e.preventDefault();
    if (busy) return;
    setBusy(true); setErr("");
    try {
      if (!(await unlock(pw))) { setErr("Wrong password"); setPw(""); input.current?.focus(); }
    } catch (x) {
      setErr(x instanceof LockUnavailable ? x.message : "Could not check the password in this browser");
    } finally { setBusy(false); }
  };

  return (
    <form className="lock-card" role="dialog" aria-modal="true" aria-labelledby="lock-title" onSubmit={submit}>
      <LockIcon size={28} className="lock-card-icon" />
      <h2 id="lock-title">This page is locked</h2>
      <input ref={input} type="password" autoComplete="off" aria-label="Password" aria-invalid={!!err} placeholder="Password"
        value={pw} onChange={e => { setPw(e.target.value); if (err) setErr(""); }} />
      <div className="lock-error" role="alert">{err}</div>
      <button type="submit" className="btn" disabled={busy}>{busy ? "Checking…" : "Unlock"}</button>
      <p className="lock-note">Password needed only for: DMPK, ICP-MS.</p>
    </form>
  );
}

/**
 * Wraps a route's page. If `required` and the session is still locked, the page is rendered underneath,
 * heavily blurred and inert, with the password card on top; otherwise it renders exactly as before.
 */
export function LockGate({ required, children }: { required: boolean; children: ReactNode }) {
  const unlocked = useUnlocked();
  const locked = required && !unlocked;
  const { pathname } = useLocation();
  const content = useRef<HTMLDivElement>(null);
  // React 18 has no `inert` prop: set it by hand so the blurred page cannot be focused or clicked
  useLayoutEffect(() => {
    const el = content.current as (HTMLDivElement & { inert?: boolean }) | null;
    if (el) el.inert = locked;
  }, [locked]);

  return (
    <div className={"lock-gate" + (locked ? " is-locked" : "")}>
      <div ref={content} className="lock-content" aria-hidden={locked ? true : undefined}>{children}</div>
      {locked && <div className="lock-overlay"><UnlockCard key={pathname} /></div>}
    </div>
  );
}
