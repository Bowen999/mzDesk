"""Reports: feedback from the sidebar and crashes, kept on the server and forwarded to the maintainer.

Every report is written to a small SQLite file first (MZDESK_FEEDBACK_DIR/reports.sqlite, default ~/.mzdesk/feedback), so
nothing is lost when a mail service is down. The maintainer reads them at /reports (admin password, api/admin.py). Then a
background thread forwards the report through whatever is configured, each channel independent of the others:

    email via Resend's HTTPS API   MZDESK_RESEND_API_KEY  (+ MZDESK_MAIL_FROM, MZDESK_FEEDBACK_TO)
    email via SMTP                 MZDESK_SMTP_HOST / _PORT / _USER / _PASSWORD / _FROM / _SECURITY
    push via ntfy                  MZDESK_NTFY_URL        (e.g. https://ntfy.sh/<a-long-random-topic>; title only, no message text)

A crash is recorded by record_crash(): the same crash again (same place, same error, same line) adds to a counter of the
first report instead of making a new one, and it is forwarded at most once an hour. deploy/README.md lists the setup.
"""
from __future__ import annotations

import base64
import hashlib
import json
import logging
import os
import secrets
import smtplib
import sqlite3
import ssl
import threading
import time
import traceback
import urllib.request
from datetime import datetime, timezone
from email.message import EmailMessage
from email.utils import formatdate, make_msgid, parseaddr
from pathlib import Path

from mzdesk import __version__

log = logging.getLogger("mzdesk.reports")

DEFAULT_TO = "feedback@example.org"
RESEND_URL = "https://api.resend.com/emails"
DEFAULT_FROM = "mzDesk <onboarding@resend.dev>"        # Resend's shared sender: it only delivers to the address the account was opened with
CRASH_WINDOW_S = 24 * 3600.0                           # the same crash within a day counts up the first report
CRASH_NOTIFY_EVERY_S = 3600.0
KEEP = 500                                             # reports kept; the oldest go first
MAX_IMAGE_BYTES = 300_000_000                          # screenshots kept in all; past that the oldest reports lose theirs (the text stays)
CRASH_FORWARDS_PER_HOUR = 20                           # crash emails / pushes in all, however many different crashes there are
_lock = threading.Lock()
_con: sqlite3.Connection | None = None
_con_path: Path | None = None


# ---------------------------------------------------------------- store
def spool_dir() -> Path:
    return Path(os.environ.get("MZDESK_FEEDBACK_DIR") or Path.home() / ".mzdesk" / "feedback")


def _db() -> sqlite3.Connection:
    """The shared connection; reopened when MZDESK_FEEDBACK_DIR changes (tests). Call with the lock held."""
    global _con, _con_path
    path = spool_dir() / "reports.sqlite"
    if _con is None or _con_path != path:
        path.parent.mkdir(parents=True, exist_ok=True)
        con = sqlite3.connect(path, check_same_thread=False, timeout=5)
        con.execute("PRAGMA journal_mode=WAL")
        con.execute("""CREATE TABLE IF NOT EXISTS reports (
                       id TEXT PRIMARY KEY, t REAL, last REAL, kind TEXT, page TEXT, title TEXT, body TEXT, detail TEXT,
                       version TEXT, agent TEXT, sig TEXT, count INTEGER DEFAULT 1, notified REAL DEFAULT 0)""")
        con.execute("CREATE INDEX IF NOT EXISTS reports_sig ON reports(sig)")
        con.execute("CREATE TABLE IF NOT EXISTS images (report TEXT, n INTEGER, name TEXT, mime TEXT, data BLOB)")
        con.execute("CREATE INDEX IF NOT EXISTS images_report ON images(report)")
        try:
            path.chmod(0o600)
        except OSError:
            pass
        _con, _con_path = con, path
    return _con


def _new_id(kind: str) -> str:
    return ("F-" if kind == "feedback" else "E-") + secrets.token_hex(3).upper()


def _trim(con: sqlite3.Connection) -> None:
    con.execute("DELETE FROM reports WHERE id NOT IN (SELECT id FROM reports ORDER BY last DESC LIMIT ?)", (KEEP,))
    con.execute("DELETE FROM images WHERE report NOT IN (SELECT id FROM reports)")
    total = con.execute("SELECT COALESCE(SUM(LENGTH(data)), 0) FROM images").fetchone()[0]
    for report, size in con.execute("SELECT i.report, SUM(LENGTH(i.data)) FROM images i JOIN reports r ON r.id = i.report GROUP BY i.report ORDER BY r.last").fetchall():
        if total <= MAX_IMAGE_BYTES:
            break
        con.execute("DELETE FROM images WHERE report = ?", (report,))
        total -= size


def add_feedback(message: str, page: str, agent: str, images: list[tuple[str, str, bytes]]) -> str:
    """Store a feedback message with its screenshots [(file name, mime, bytes)] -> its reference."""
    now = time.time()
    first = " ".join(message.strip().splitlines()[:1])[:100]
    with _lock:
        con = _db()
        ref = _new_id("feedback")
        con.execute("INSERT INTO reports (id, t, last, kind, page, title, body, detail, version, agent, sig) VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                    (ref, now, now, "feedback", page[:200], first, message.strip(), "", __version__, agent[:300], ""))
        for n, (name, mime, data) in enumerate(images, 1):
            con.execute("INSERT INTO images VALUES (?,?,?,?,?)", (ref, n, name, mime, data))
        _trim(con)
        con.commit()
    return ref


def _signature(where: str, kind: str, frame: str) -> str:
    return hashlib.sha1(f"{where}|{kind}|{frame}".encode()).hexdigest()[:16]


def add_crash(where: str, error: str, trace: str, frame: str, page: str = "", agent: str = "", context: str = "") -> tuple[str, bool]:
    """Store a crash -> (reference, forward it). A repeat of a crash from the last day only counts up the first report."""
    now = time.time()
    sig = _signature(where, error.split(":", 1)[0], frame)
    with _lock:
        con = _db()
        row = con.execute("SELECT id, notified FROM reports WHERE sig=? AND last > ? ORDER BY last DESC LIMIT 1", (sig, now - CRASH_WINDOW_S)).fetchone()
        if row:
            ref, notified = row
            con.execute("UPDATE reports SET count = count + 1, last = ? WHERE id = ?", (now, ref))
            forward = now - (notified or 0) > CRASH_NOTIFY_EVERY_S
        else:
            ref = _new_id("crash")
            forward = con.execute("SELECT COUNT(*) FROM reports WHERE kind = 'crash' AND notified > ?", (now - 3600,)).fetchone()[0] < CRASH_FORWARDS_PER_HOUR
            con.execute("INSERT INTO reports (id, t, last, kind, page, title, body, detail, version, agent, sig) VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                        (ref, now, now, "crash", page[:200] or where[:200], f"{where}: {error}"[:140], trace[-6000:], context[:2000], __version__, agent[:300], sig))
            _trim(con)
        if forward:
            con.execute("UPDATE reports SET notified = ? WHERE id = ?", (now, ref))
        con.commit()
    return ref, forward


def record_crash(exc: BaseException, where: str, inputs: dict | None = None, log_it: bool = True) -> str:
    """Log and store an exception from the server's own code -> its reference ('E-1A2B3C') for the message the user sees.
    `inputs` (the form of the run) is described by field names and sizes only; the values never leave the run."""
    trace = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
    tb = traceback.extract_tb(exc.__traceback__)
    frame = f"{Path(tb[-1].filename).name}:{tb[-1].name}:{tb[-1].lineno}" if tb else ""
    error = f"{type(exc).__name__}: {exc}"
    shape = ""
    if isinstance(inputs, dict):
        shape = "inputs: " + ", ".join(f"{k} ({len(json.dumps(v, default=str))} chars)" for k, v in list(inputs.items())[:40])
    try:
        ref, forward = add_crash(where, error, trace, frame, context=shape)
    except Exception:  # noqa: BLE001 - reporting must never turn one failure into two
        log.exception("could not store a crash report")
        return "E-000000"
    if log_it:
        log.error("crash %s in %s: %s\n%s", ref, where, error, trace)
    if forward:
        forward_async(ref)
    return ref


class _CrashHandler(logging.Handler):
    """Every log.exception(...) of the mzdesk code (a step that failed, a figure that could not be drawn) becomes a crash report too."""

    def __init__(self):
        super().__init__(level=logging.ERROR)
        self._busy = threading.local()

    def emit(self, record):
        if not record.exc_info or record.exc_info[1] is None or record.name.startswith(("mzdesk.reports", "mzdesk.feedback")) or getattr(self._busy, "on", False):
            return
        self._busy.on = True                              # reporting a crash must not report its own failure
        try:
            record_crash(record.exc_info[1], f"{record.name}: {record.msg}"[:200], log_it=False)    # the template, not the text: ids in it would make every crash look new
        except Exception:  # noqa: BLE001
            pass
        finally:
            self._busy.on = False


def install_log_handler() -> None:
    logger = logging.getLogger("mzdesk")
    if not any(isinstance(h, _CrashHandler) for h in logger.handlers):
        logger.addHandler(_CrashHandler())


def add_browser_crash(message: str, stack: str, page: str, agent: str) -> tuple[str, bool]:
    """A crash seen by the browser (a page that failed to draw). The frame is the first stack line that is not the message."""
    lines = [ln.strip() for ln in stack.splitlines() if ln.strip()]
    frame = lines[1] if len(lines) > 1 else (lines[0] if lines else "")
    frame = frame.split("(", 1)[-1].rsplit("/", 1)[-1][:120]
    return add_crash("browser", message[:300], stack, frame, page=page, agent=agent)


# ---------------------------------------------------------------- reading (admin)
def listing(limit: int = 200) -> list[dict]:
    with _lock:
        con = _db()
        rows = con.execute("""SELECT r.id, r.t, r.last, r.kind, r.page, r.title, r.version, r.count,
                              (SELECT COUNT(*) FROM images i WHERE i.report = r.id) FROM reports r ORDER BY r.last DESC LIMIT ?""", (limit,)).fetchall()
    return [{"id": i, "time": _iso(t), "last": _iso(last), "kind": k, "page": p, "title": title, "version": v, "count": c, "images": n}
            for i, t, last, k, p, title, v, c, n in rows]


def detail(ref: str) -> dict | None:
    with _lock:
        con = _db()
        row = con.execute("SELECT id, t, last, kind, page, title, body, detail, version, agent, count FROM reports WHERE id=?", (ref,)).fetchone()
        if row is None:
            return None
        images = con.execute("SELECT name, mime, data FROM images WHERE report=? ORDER BY n", (ref,)).fetchall()
    i, t, last, kind, page, title, body, extra, version, agent, count = row
    return {"id": i, "time": _iso(t), "last": _iso(last), "kind": kind, "page": page, "title": title, "body": body, "detail": extra,
            "version": version, "browser": agent, "count": count,
            "images": [{"name": n, "src": f"data:{m};base64,{base64.b64encode(d).decode()}"} for n, m, d in images]}


def delete(ref: str) -> bool:
    with _lock:
        con = _db()
        gone = con.execute("DELETE FROM reports WHERE id=?", (ref,)).rowcount
        con.execute("DELETE FROM images WHERE report=?", (ref,))
        con.commit()
    return bool(gone)


def _iso(t: float) -> str:
    return datetime.fromtimestamp(t, timezone.utc).isoformat(timespec="seconds")


# ---------------------------------------------------------------- forwarding
def _to_address() -> str:
    return (os.environ.get("MZDESK_FEEDBACK_TO") or DEFAULT_TO).strip() or DEFAULT_TO


def _smtp_config():
    """None when SMTP is not configured; a half-finished configuration counts as missing."""
    host = (os.environ.get("MZDESK_SMTP_HOST") or "").strip()
    if not host:
        return None
    try:
        port = int((os.environ.get("MZDESK_SMTP_PORT") or "587").strip())
    except ValueError:
        port = 587
    user = (os.environ.get("MZDESK_SMTP_USER") or "").strip()
    sender = (os.environ.get("MZDESK_SMTP_FROM") or user).strip()
    security = (os.environ.get("MZDESK_SMTP_SECURITY") or "starttls").strip().lower()
    if not sender or security not in ("starttls", "ssl", "none"):
        log.error("reports: incomplete SMTP configuration (need MZDESK_SMTP_FROM or MZDESK_SMTP_USER)")
        return None
    return {"host": host, "port": port, "user": user, "password": os.environ.get("MZDESK_SMTP_PASSWORD", ""), "sender": sender, "security": security}


def compose(report: dict, sender: str) -> EmailMessage:
    """The report as an email: text body, screenshots attached."""
    crash = report["kind"] == "crash"
    msg = EmailMessage()
    msg["Subject"] = f"[mzDesk] {'Crash' if crash else 'Feedback'} {report['id']}: {report['title']}"[:200]
    msg["From"] = sender
    msg["To"] = _to_address()
    msg["Date"] = formatdate(localtime=True)
    address = parseaddr(sender)[1]                       # the sender's own domain in the Message-ID, not a fake TLD
    msg["Message-ID"] = make_msgid(domain=address.rsplit("@", 1)[-1] if "@" in address else "mzdesk.local")
    lines = [f"{'Crash' if crash else 'Feedback'} {report['id']} from mzDesk.", "", report["body"].strip()]
    if report.get("detail"):
        lines += ["", report["detail"]]
    lines += ["", "---", f"Page: {report['page'] or '(unknown)'}", f"Version: {report['version']}", f"Submitted: {report['time']}",
              f"Browser: {report.get('browser') or '(unknown)'}", f"Seen: {report['count']} time(s)" if crash else f"Images: {len(report['images'])}"]
    msg.set_content("\n".join(lines) + "\n")
    for image in report["images"]:
        header, payload = image["src"].split(",", 1)
        maintype, subtype = header[5:].split(";")[0].split("/", 1)
        msg.add_attachment(base64.b64decode(payload), maintype=maintype, subtype=subtype, filename=image["name"])
    return msg


def _send_resend(msg: EmailMessage, key: str) -> None:
    payload = {"from": msg["From"], "to": [msg["To"]], "subject": msg["Subject"],
               "text": msg.get_body(preferencelist=("plain",)).get_content(),
               "attachments": [{"filename": p.get_filename(), "content": base64.b64encode(p.get_content()).decode()} for p in msg.iter_attachments()]}
    request = urllib.request.Request(os.environ.get("MZDESK_RESEND_URL") or RESEND_URL, data=json.dumps(payload).encode(), method="POST",
                                     headers={"Authorization": f"Bearer {key}", "Content-Type": "application/json", "User-Agent": f"mzDesk/{__version__}"})
    with urllib.request.urlopen(request, timeout=20) as response:
        response.read()


def _send_smtp(cfg: dict, msg: EmailMessage) -> None:
    context = ssl.create_default_context()
    server = smtplib.SMTP_SSL(cfg["host"], cfg["port"], timeout=20, context=context) if cfg["security"] == "ssl" else smtplib.SMTP(cfg["host"], cfg["port"], timeout=20)
    try:
        if cfg["security"] == "starttls":
            server.starttls(context=context)
        if cfg["user"]:
            server.login(cfg["user"], cfg["password"])
        server.send_message(msg)
    finally:
        try:
            server.quit()
        except Exception:  # noqa: BLE001 - the message is sent, or the connection is gone
            pass


def _send_ntfy(url: str, report: dict) -> None:
    """A push with the reference, the kind and the page. The message text stays on the server."""
    crash = report["kind"] == "crash"
    text = f"{report['id']} · {report['page'] or 'mzDesk'}" + (f" · {report['title']}" if crash else "")
    request = urllib.request.Request(url, data=text.encode("utf-8"), method="POST",
                                     headers={"Title": f"mzDesk {'crash' if crash else 'feedback'}", "Tags": "rotating_light" if crash else "speech_balloon",
                                              "User-Agent": f"mzDesk/{__version__}"})
    with urllib.request.urlopen(request, timeout=15) as response:
        response.read()


def channels() -> list[str]:
    """The forwarding channels that are configured (names only)."""
    found = []
    if (os.environ.get("MZDESK_RESEND_API_KEY") or "").strip():
        found.append("resend")
    if _smtp_config():
        found.append("smtp")
    if (os.environ.get("MZDESK_NTFY_URL") or "").strip():
        found.append("ntfy")
    return found


def forward(ref: str) -> list[str]:
    """Send one stored report through every configured channel -> the channels that worked. Never raises."""
    report = detail(ref)
    if report is None:
        return []
    key = (os.environ.get("MZDESK_RESEND_API_KEY") or "").strip()
    smtp = _smtp_config()
    ntfy = (os.environ.get("MZDESK_NTFY_URL") or "").strip()
    attempts = []
    if key:
        attempts.append(("resend", lambda: _send_resend(compose(report, (os.environ.get("MZDESK_MAIL_FROM") or DEFAULT_FROM).strip()), key)))
    if smtp:
        attempts.append(("smtp", lambda: _send_smtp(smtp, compose(report, smtp["sender"]))))
    if ntfy:
        attempts.append(("ntfy", lambda: _send_ntfy(ntfy, report)))
    done = []
    for name, send in attempts:
        try:
            send()
            done.append(name)
        except Exception:  # noqa: BLE001 - one channel failing must not stop the next; the report is stored anyway
            log.exception("reports: forwarding %s by %s failed", ref, name)
    return done


def forward_async(ref: str) -> None:
    if channels():
        threading.Thread(target=forward, args=(ref,), daemon=True, name=f"report-{ref}").start()
