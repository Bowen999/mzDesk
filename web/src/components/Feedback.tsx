import { createPortal } from "react-dom";
import { useCallback, useEffect, useRef, useState, type ClipboardEvent, type DragEvent, type FormEvent } from "react";
import { ApiError, sendFeedback } from "../api";
import { t, useT } from "../i18n";

/** Simple beetle: body, head, antennae and six legs. Kept to one stroke width so it reads at 16 px. */
export function BugIcon({ size = 17, className = "" }: { size?: number; className?: string }) {
  return (
    <svg className={className} width={size} height={size} viewBox="0 0 24 24" fill="none" stroke="currentColor"
      strokeWidth="1.7" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true" focusable="false">
      <path d="M9 7.5a3 3 0 0 1 6 0v1" />
      <rect x="7.5" y="8.5" width="9" height="10" rx="4.5" />
      <path d="M12 9v9" />
      <path d="M9 8.5 7 6M15 8.5 17 6" />
      <path d="M7.5 12H4.5M7.5 15H5M9 17.5 7 20M16.5 12h3M16.5 15H19M15 17.5l2 2.5" />
    </svg>
  );
}

const MAX_IMAGES = 3;
const MAX_MESSAGE = 4000;
const MAX_SOURCE_BYTES = 20 * 1024 * 1024;     // before the browser scales the file down
const MAX_IMAGE_BYTES = 2_500_000;            // after scaling, matches the server
const MAX_TOTAL_BYTES = 5_000_000;
const MAX_EDGE = 1600;

type Shot = { id: string; name: string; dataUrl: string; bytes: number };

function dataUrlBytes(dataUrl: string) {
  const b64 = dataUrl.slice(dataUrl.indexOf(",") + 1);
  const pad = b64.endsWith("==") ? 2 : b64.endsWith("=") ? 1 : 0;
  return Math.floor((b64.length * 3) / 4) - pad;
}

/** Scale to at most 1600 px and re-encode (PNG for screenshots, JPEG for photos) so the email stays small. */
async function prepareImage(file: File): Promise<Shot> {
  if (!file.type.startsWith("image/")) throw new Error(t("Only image files can be attached."));
  if (file.size > MAX_SOURCE_BYTES) throw new Error(t("{name} is larger than 20 MB.", { name: file.name }));
  let bitmap: ImageBitmap;
  try {
    bitmap = await createImageBitmap(file);
  } catch {
    throw new Error(t("Could not read {name}. Please use a PNG, JPG, WebP or GIF image.", { name: file.name }));
  }
  try {
    const scale = Math.min(1, MAX_EDGE / Math.max(bitmap.width, bitmap.height));
    const width = Math.max(1, Math.round(bitmap.width * scale));
    const height = Math.max(1, Math.round(bitmap.height * scale));
    const canvas = document.createElement("canvas");
    canvas.width = width;
    canvas.height = height;
    const ctx = canvas.getContext("2d");
    if (!ctx) throw new Error(t("Could not prepare that image."));
    const png = file.type !== "image/jpeg";      // GIF / WebP screenshots keep transparency as PNG
    if (!png) {                                  // JPEG has no alpha: use white, not black
      ctx.fillStyle = "#ffffff";
      ctx.fillRect(0, 0, width, height);
    }
    ctx.drawImage(bitmap, 0, 0, width, height);
    const dataUrl = canvas.toDataURL(png ? "image/png" : "image/jpeg", 0.86);
    const bytes = dataUrlBytes(dataUrl);
    if (bytes > MAX_IMAGE_BYTES) throw new Error(t("{name} is still too large after resizing. Try a smaller screenshot.", { name: file.name }));
    return { id: Math.random().toString(36).slice(2), name: file.name, dataUrl, bytes };
  } finally {
    bitmap.close();
  }
}

/**
 * "Send feedback" in the sidebar footer. Opens a dialog with a message box and up to three screenshots
 * (file picker, drag & drop, or paste). The dialog is portalled to <body> so the scrolling sidebar
 * cannot clip it; Esc, the X and the backdrop close it.
 */
export function FeedbackButton({ rail = false }: { rail?: boolean }) {
  const t = useT();
  const [open, setOpen] = useState(false);
  const trigger = useRef<HTMLButtonElement>(null);
  const close = useCallback(() => setOpen(false), []);
  return (
    <>
      <button ref={trigger} type="button" className="feedback-btn" aria-label={t("Send feedback")}
        data-label={rail ? undefined : t("Send feedback")} title={rail ? t("Send feedback") : undefined}
        aria-haspopup="dialog" aria-expanded={open} onClick={() => setOpen(true)}>
        <BugIcon />
      </button>
      {open && <FeedbackDialog onClose={close} />}
    </>
  );
}

function FeedbackDialog({ onClose }: { onClose: () => void }) {
  const t = useT();
  const [message, setMessage] = useState("");
  const [shots, setShots] = useState<Shot[]>([]);
  const [website, setWebsite] = useState("");          // honeypot; people never see it
  const [error, setError] = useState("");
  const [adding, setAdding] = useState(false);
  const [sending, setSending] = useState(false);
  const [done, setDone] = useState<{ ref: string } | null>(null);
  const [over, setOver] = useState(false);
  const card = useRef<HTMLDivElement>(null);
  const textarea = useRef<HTMLTextAreaElement>(null);
  const fileInput = useRef<HTMLInputElement>(null);
  const closeButton = useRef<HTMLButtonElement>(null);
  const sendingRef = useRef(false);

  useEffect(() => { sendingRef.current = sending; }, [sending]);

  // Focus management, body scroll lock, Esc, Tab trap; restore focus to the trigger when the dialog closes.
  useEffect(() => {
    const previous = document.activeElement as HTMLElement | null;
    const overflow = document.body.style.overflow;
    document.body.style.overflow = "hidden";
    const timer = window.setTimeout(() => textarea.current?.focus(), 40);
    const key = (e: KeyboardEvent) => {
      if (e.key === "Escape" && !sendingRef.current) { e.preventDefault(); onClose(); return; }
      if (e.key !== "Tab" || !card.current) return;
      const nodes = Array.from(card.current.querySelectorAll<HTMLElement>(
        'a[href], button:not([disabled]), textarea:not([disabled]), input:not([disabled]):not([type="hidden"]), [tabindex]:not([tabindex="-1"])',
      )).filter(n => n.tabIndex !== -1 && n.offsetParent !== null);
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
  }, [onClose]);

  useEffect(() => { if (done) closeButton.current?.focus(); }, [done]);

  const addFiles = async (files: File[]) => {
    if (sending || adding) return;
    setError("");
    if (!files.length) return;
    const room = MAX_IMAGES - shots.length;
    if (room <= 0) { setError(t("You can attach up to {max} images.", { max: MAX_IMAGES })); return; }
    if (files.length > room) setError(t("Up to {max} images can be attached; only the first {room} were added.", { max: MAX_IMAGES, room }));
    setAdding(true);
    const added: Shot[] = [];
    let total = shots.reduce((n, s) => n + s.bytes, 0);
    let problem = "";
    try {
      for (const file of files.slice(0, room)) {
        try {
          const shot = await prepareImage(file);
          if (total + shot.bytes > MAX_TOTAL_BYTES) { problem = t("The images are too large together (5 MB total). Remove one and try again."); break; }
          added.push(shot);
          total += shot.bytes;
        } catch (e) {
          problem = e instanceof Error ? e.message : t("Could not add that image.");
        }
      }
    } finally {
      setAdding(false);
    }
    if (added.length) setShots(prev => [...prev, ...added]);
    if (problem) setError(problem);
  };

  const onPaste = (e: ClipboardEvent<HTMLTextAreaElement>) => {
    const files = Array.from(e.clipboardData.items)
      .filter(item => item.kind === "file" && item.type.startsWith("image/"))
      .map(item => item.getAsFile())
      .filter((f): f is File => !!f);
    if (files.length) void addFiles(files);
  };

  const onDrop = (e: DragEvent<HTMLDivElement>) => {
    e.preventDefault();
    setOver(false);
    void addFiles(Array.from(e.dataTransfer.files));
  };

  const submit = async (e: FormEvent) => {
    e.preventDefault();
    if (sending || adding || message.trim().length < 3) return;
    setSending(true);
    setError("");
    try {
      const result = await sendFeedback({
        message: message.trim(),
        page: window.location.pathname,
        website,
        images: shots.map(s => ({ name: s.name, data: s.dataUrl })),
      });
      setDone({ ref: result.ref ?? "" });
    } catch (x) {
      setError(x instanceof ApiError ? x.message : t("Could not reach the server. Please try again."));
    } finally {
      setSending(false);
    }
  };

  const full = shots.length >= MAX_IMAGES;
  const canSend = !sending && !adding && message.trim().length >= 3;

  return createPortal(
    <div className="feedback-overlay" onMouseDown={e => { if (e.target === e.currentTarget && !sendingRef.current) onClose(); }}>
      <div ref={card} className="feedback-card" role="dialog" aria-modal="true" aria-labelledby="feedback-title"
        aria-describedby={done ? undefined : "feedback-desc"}>
        <button type="button" className="feedback-close" aria-label={t("Close feedback")} title={t("Close")} onClick={onClose} disabled={sending}>
          <svg viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="currentColor" strokeWidth="1.8" strokeLinecap="round" aria-hidden="true"><path d="M6 6l12 12M18 6 6 18" /></svg>
        </button>

        {done ? (
          <div className="feedback-done" role="status">
            <span className="feedback-done-icon" aria-hidden="true">
              <svg viewBox="0 0 24 24" width="22" height="22" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round"><path d="m5 12 5 5L20 7" /></svg>
            </span>
            <h2 id="feedback-title">{t("Thanks for the feedback")}</h2>
            <p>{t("Your message reached the mzDesk team.")}{done.ref && <> <code>{done.ref}</code></>}</p>
            <button ref={closeButton} type="button" className="btn" onClick={onClose}>{t("Close")}</button>
          </div>
        ) : (
          <>
            <span className="feedback-badge" aria-hidden="true"><BugIcon size={20} /></span>
            <h2 id="feedback-title">{t("Send feedback")}</h2>
            <p className="feedback-desc" id="feedback-desc">{t("Report a bug, suggest an idea, or share what’s working well.")}</p>

            <form onSubmit={submit}>
              <div className="feedback-field">
                <textarea ref={textarea} value={message} maxLength={MAX_MESSAGE} onChange={e => { setMessage(e.target.value); if (error) setError(""); }}
                  onPaste={onPaste} placeholder={t("What happened? What were you trying to do?")} aria-label={t("Your message")}
                  rows={5} disabled={sending} />
                <span className="feedback-count" aria-hidden="true">{message.length}/{MAX_MESSAGE}</span>
              </div>

              <div className={"feedback-attach" + (over ? " over" : "")}
                onDragOver={e => { e.preventDefault(); if (!sending) setOver(true); }}
                onDragLeave={() => setOver(false)} onDrop={onDrop}>
                <button type="button" className="feedback-add" onClick={() => fileInput.current?.click()}
                  disabled={sending || adding || full}>
                  <svg viewBox="0 0 24 24" width="15" height="15" fill="none" stroke="currentColor" strokeWidth="1.7" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true">
                    <rect x="3" y="3" width="18" height="18" rx="3" /><circle cx="8.5" cy="8.5" r="1.5" /><path d="m21 15-5-5L5 21" /><path d="M16 5h6M19 2v6" />
                  </svg>
                  {adding ? t("Adding…") : full ? t("3 images attached") : shots.length ? t("Add another image") : t("Add images")}
                </button>
                <span className="feedback-attach-hint">{t("or paste a screenshot · PNG, JPG, WebP · up to 3")}</span>
                <input ref={fileInput} type="file" accept="image/png,image/jpeg,image/webp,image/gif" multiple hidden
                  onChange={e => { const picked = Array.from(e.target.files ?? []); e.target.value = ""; void addFiles(picked); }} />
              </div>

              {shots.length > 0 && (
                <ul className="feedback-shots" aria-label={t("Attached images")}>
                  {shots.map(s => (
                    <li key={s.id}>
                      <img src={s.dataUrl} alt={s.name} title={s.name} />
                      <button type="button" className="feedback-shot-x" aria-label={t("Remove {name}", { name: s.name })} title={t("Remove")}
                        onClick={() => setShots(prev => prev.filter(p => p.id !== s.id))} disabled={sending}>×</button>
                    </li>
                  ))}
                </ul>
              )}

              <input className="feedback-hp" type="text" name="website" tabIndex={-1} autoComplete="off"
                aria-hidden="true" value={website} onChange={e => setWebsite(e.target.value)} />

              {error && <p className="feedback-error" role="alert">{error}</p>}

              <div className="feedback-actions">
                <span className="feedback-note">{t("Your message and images go to the mzDesk maintainer.")}</span>
                <button type="button" className="btn ghost" onClick={onClose} disabled={sending}>{t("Cancel")}</button>
                <button type="submit" className="btn" disabled={!canSend} aria-busy={sending}>
                  {sending ? t("Sending…") : t("Send")}
                </button>
              </div>
            </form>
          </>
        )}
      </div>
    </div>,
    document.body,
  );
}
