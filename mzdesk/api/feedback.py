"""Feedback from the app to the maintainer's inbox.

The browser posts JSON (so the API needs no multipart dependency):
    {"message": "...", "page": "/modules/ms1_search", "website": "",
     "images": [{"name": "screenshot.png", "data": "data:image/png;base64,..."}]}
Images are already scaled down in the browser; the server checks the real file type and the size, so only images are
ever attached. `website` is a honeypot: people never see it, bots fill it in.

Delivery is SMTP. Configure it with environment variables (deploy/mail.env.example lists them all). When
MZDESK_SMTP_HOST is missing, or delivery fails, the message is written as an .eml file under
MZDESK_FEEDBACK_DIR (default ~/.mzdesk/feedback/) instead of being lost; the API then says delivered=false.
"""
from __future__ import annotations

import base64
import binascii
import json
import logging
import os
import re
import smtplib
import ssl
import tempfile
import threading
import time
import uuid
from datetime import datetime, timezone
from email.message import EmailMessage
from email.utils import formatdate, make_msgid, parseaddr
from pathlib import Path

from fastapi import HTTPException, Request
from starlette.concurrency import run_in_threadpool

from mzdesk import __version__

log = logging.getLogger("mzdesk.feedback")

DEFAULT_TO = "feedback@example.org"
MAX_BODY_BYTES = 8 * 1024 * 1024          # nginx also caps the request at 8 MB
MAX_MESSAGE = 4000
MIN_MESSAGE = 3
MAX_IMAGES = 3
MAX_IMAGE_BYTES = 2_500_000               # per image, after the browser compressed it
MAX_TOTAL_IMAGE_BYTES = 5_000_000

# One sender per IP: a short pause between messages and a small hourly budget. The app runs two workers, so this is
# per worker; nginx adds a stricter shared limit on the same path.
_RATE_LOCK = threading.Lock()
_HITS: dict[str, list[float]] = {}
COOLDOWN_S, HOURLY_LIMIT, HOUR_S = 30.0, 8, 3600.0


def _to_address() -> str:
    return (os.environ.get("MZDESK_FEEDBACK_TO") or DEFAULT_TO).strip() or DEFAULT_TO


def spool_dir() -> Path:
    return Path(os.environ.get("MZDESK_FEEDBACK_DIR") or Path.home() / ".mzdesk" / "feedback")


def _rate_limited(ip: str) -> str:
    """Empty when the request may go through, else the message to show."""
    now = time.time()
    with _RATE_LOCK:
        hits = [t for t in _HITS.get(ip, []) if now - t < HOUR_S]
        if hits and now - hits[-1] < COOLDOWN_S:
            _HITS[ip] = hits
            return f"Please wait {max(1, int(COOLDOWN_S - (now - hits[-1])))} s before sending another message."
        if len(hits) >= HOURLY_LIMIT:
            _HITS[ip] = hits
            return "Too many messages from this address. Please try again later."
        hits.append(now)
        _HITS[ip] = hits
        if len(_HITS) > 4096:                              # many distinct addresses: drop the inactive ones
            for key in [k for k, v in _HITS.items() if not v or v[-1] < now - HOUR_S]:
                del _HITS[key]
    return ""


def _sniff_image(data: bytes):
    """(mime, extension) of a real image, or None. Trusts the bytes, never the browser's content type."""
    if data.startswith(b"\x89PNG\r\n\x1a\n"):
        return "image/png", "png"
    if data.startswith(b"\xff\xd8\xff"):
        return "image/jpeg", "jpg"
    if data[:6] in (b"GIF87a", b"GIF89a"):
        return "image/gif", "gif"
    if data[:4] == b"RIFF" and data[8:12] == b"WEBP":
        return "image/webp", "webp"
    return None


_DATA_URL = re.compile(r"^data:image/[a-z0-9.+-]+;base64,(.+)$", re.I | re.S)


def _decode_images(raw) -> list[tuple[str, str, bytes]]:
    """Validate the images list -> [(filename, mime, bytes)]. Only real PNG / JPEG / WebP / GIF pass."""
    if raw is None:
        return []
    if not isinstance(raw, list):
        raise HTTPException(400, "images must be a list.")
    if len(raw) > MAX_IMAGES:
        raise HTTPException(400, f"Attach at most {MAX_IMAGES} images.")
    out, total = [], 0
    for i, item in enumerate(raw, 1):
        data_url = item.get("data") if isinstance(item, dict) else None
        if not isinstance(data_url, str):
            raise HTTPException(400, "Every image needs a data URL.")
        m = _DATA_URL.match(data_url.strip())
        if not m:
            raise HTTPException(400, "Images must be PNG, JPEG, WebP or GIF files.")
        payload = m.group(1)
        if len(payload) > (MAX_IMAGE_BYTES * 4) // 3 + 8:       # do not decode something that cannot fit
            raise HTTPException(400, "Each image must be under 2.5 MB.")
        try:
            data = base64.b64decode(payload, validate=True)
        except (binascii.Error, ValueError):
            raise HTTPException(400, "That image could not be read.") from None
        if len(data) > MAX_IMAGE_BYTES:
            raise HTTPException(400, "Each image must be under 2.5 MB.")
        kind = _sniff_image(data)
        if kind is None:
            raise HTTPException(400, "Images must be PNG, JPEG, WebP or GIF files.")
        total += len(data)
        if total > MAX_TOTAL_IMAGE_BYTES:
            raise HTTPException(400, "Images must be under 5 MB together.")
        mime, ext = kind
        wanted = str(item.get("name") or "")[:80]
        stem = re.sub(r"[^A-Za-z0-9._-]+", "-", Path(wanted).stem).strip(".-") or f"screenshot-{i}"
        out.append((f"{stem[:40]}.{ext}", mime, data))
    return out


def _smtp_config():
    """None when mail is not configured; otherwise a dict. A half-finished configuration is treated as missing."""
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
        log.error("feedback: incomplete SMTP configuration (need MZDESK_SMTP_FROM or MZDESK_SMTP_USER)")
        return None
    return {"host": host, "port": port, "user": user, "password": os.environ.get("MZDESK_SMTP_PASSWORD", ""),
            "sender": sender, "security": security}


def _compose(message: str, page: str, images, request: Request, sender: str) -> EmailMessage:
    first = " ".join(message.strip().splitlines()[:1])[:72]
    msg = EmailMessage()
    msg["Subject"] = f"[mzDesk] Feedback: {first}" if first else "[mzDesk] Feedback"
    msg["From"] = sender
    msg["To"] = _to_address()
    msg["Date"] = formatdate(localtime=True)
    address = parseaddr(sender)[1]                       # use the sender's domain, not a fake TLD, in the Message-ID
    msg["Message-ID"] = make_msgid(domain=address.rsplit("@", 1)[-1] if "@" in address else "mzdesk.local")
    agent = (request.headers.get("user-agent") or "")[:300]
    submitted = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    msg.set_content(
        "New feedback from mzDesk.\n\n"
        f"{message.strip()}\n\n"
        "---\n"
        f"Page: {page or '(unknown)'}\n"
        f"Version: {__version__}\n"
        f"Submitted: {submitted}\n"
        f"Browser: {agent or '(unknown)'}\n"
        f"Images: {len(images)}\n"
    )
    for filename, mime, data in images:
        maintype, subtype = mime.split("/", 1)
        msg.add_attachment(data, maintype=maintype, subtype=subtype, filename=filename)
    return msg


def _spool(msg: EmailMessage) -> Path:
    folder = spool_dir()
    try:
        folder.mkdir(parents=True, exist_ok=True)
    except OSError:
        # A read-only home (or a sandbox) must not lose the message: fall back to the system temp folder.
        folder = Path(tempfile.gettempdir()) / "mzdesk-feedback"
        folder.mkdir(parents=True, exist_ok=True)
        log.warning("feedback: cannot write %s; using %s", spool_dir(), folder)
    path = folder / f"{datetime.now(timezone.utc):%Y%m%dT%H%M%SZ}-{uuid.uuid4().hex[:8]}.eml"
    path.write_bytes(bytes(msg))
    try:
        path.chmod(0o600)
    except OSError:
        pass
    return path


def _send(cfg, msg: EmailMessage):
    context = ssl.create_default_context()
    if cfg["security"] == "ssl":
        server = smtplib.SMTP_SSL(cfg["host"], cfg["port"], timeout=20, context=context)
    else:
        server = smtplib.SMTP(cfg["host"], cfg["port"], timeout=20)
    try:
        if cfg["security"] == "starttls":
            server.starttls(context=context)
        if cfg["user"]:
            server.login(cfg["user"], cfg["password"])
        server.send_message(msg)
    finally:
        try:
            server.quit()
        except Exception:  # noqa: BLE001 - the message is already sent (or the connection is gone)
            pass


def _deliver(msg: EmailMessage, cfg) -> tuple[bool, str]:
    """-> (delivered, note for the user). Never raises for a delivery problem: the .eml copy is the fallback."""
    if cfg is None:
        try:
            path = _spool(msg)
        except Exception:  # noqa: BLE001
            log.exception("feedback: could not save the message")
            raise HTTPException(503, "Could not save your message right now. Please try again later.") from None
        log.warning("feedback: MZDESK_SMTP_HOST is not set; saved the message to %s", path)
        return False, "Your message was saved on the server; email delivery is not configured yet."
    try:
        _send(cfg, msg)
        return True, ""
    except Exception:  # noqa: BLE001 - SMTP errors must not lose the message
        log.exception("feedback: SMTP delivery failed")
    try:
        path = _spool(msg)
    except Exception:  # noqa: BLE001
        log.exception("feedback: could not save the message after SMTP failed")
        raise HTTPException(502, "The message could not be sent. Please try again later.") from None
    log.warning("feedback: saved the message to %s after SMTP failed", path)
    return False, "Your message was saved on the server, but the email could not be sent right now."


async def submit(request: Request) -> dict:
    """POST /api/feedback: validate, rate-limit, then hand the SMTP call to a thread."""
    length = request.headers.get("content-length")
    if length:
        try:
            if int(length) > MAX_BODY_BYTES:
                raise HTTPException(413, "The message is too large. Attach fewer or smaller images.")
        except ValueError:
            pass
    raw = await request.body()
    if len(raw) > MAX_BODY_BYTES:
        raise HTTPException(413, "The message is too large. Attach fewer or smaller images.")
    try:
        body = json.loads(raw or b"{}")
    except (ValueError, UnicodeDecodeError):
        raise HTTPException(400, "Expected a JSON body.") from None
    if not isinstance(body, dict):
        raise HTTPException(400, "Expected a JSON object.")

    if str(body.get("website") or "").strip():          # honeypot: pretend success, send nothing
        return {"ok": True, "delivered": True, "note": ""}

    message = str(body.get("message") or "").strip()
    if len(message) < MIN_MESSAGE:
        raise HTTPException(400, "Please write a short message.")
    if len(message) > MAX_MESSAGE:
        raise HTTPException(400, f"Please keep the message under {MAX_MESSAGE} characters.")
    page = str(body.get("page") or "")[:200]
    if page and not page.startswith("/"):
        page = "/" + page
    images = _decode_images(body.get("images"))

    ip = request.client.host if request.client else "unknown"
    wait = _rate_limited(ip)                              # only valid messages use up the sender's budget
    if wait:
        raise HTTPException(429, wait)

    cfg = _smtp_config()
    sender = cfg["sender"] if cfg else "mzDesk feedback <mzdesk@localhost>"
    msg = _compose(message, page, images, request, sender)
    delivered, note = await run_in_threadpool(_deliver, msg, cfg)
    return {"ok": True, "delivered": delivered, "note": note}
