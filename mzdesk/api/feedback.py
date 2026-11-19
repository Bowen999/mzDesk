"""Feedback and browser crash reports from the app to the maintainer.

The browser posts JSON (so the API needs no multipart dependency):
    {"message": "...", "page": "/modules/ms1_search", "website": "",
     "images": [{"name": "screenshot.png", "data": "data:image/png;base64,..."}]}
Images are already scaled down in the browser; the server checks the real file type and the size, so only images are
ever stored. `website` is a honeypot: people never see it, bots fill it in.

A report is stored first (mzdesk/api/reports.py) and then forwarded by email / push when that is configured; the form
reports success once it is stored. A page that crashed in the browser posts {"message", "stack", "page"} to /api/crash.
"""
from __future__ import annotations

import base64
import binascii
import json
import logging
import re
import threading
import time
from pathlib import Path

from fastapi import HTTPException, Request

from mzdesk.api import reports

log = logging.getLogger("mzdesk.feedback")

MAX_BODY_BYTES = 8 * 1024 * 1024          # nginx also caps the request at 8 MB
MAX_MESSAGE = 4000
MIN_MESSAGE = 3
MAX_IMAGES = 3
MAX_IMAGE_BYTES = 2_500_000               # per image, after the browser compressed it
MAX_TOTAL_IMAGE_BYTES = 5_000_000

# One sender per IP: a short pause between messages and a small hourly budget. nginx adds a stricter shared limit on the same path.
_RATE_LOCK = threading.Lock()
_HITS: dict[str, list[float]] = {}
COOLDOWN_S, HOURLY_LIMIT, HOUR_S = 30.0, 8, 3600.0


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


async def _json_body(request: Request, too_big: str) -> dict:
    length = request.headers.get("content-length")
    if length:
        try:
            if int(length) > MAX_BODY_BYTES:
                raise HTTPException(413, too_big)
        except ValueError:
            pass
    raw = await request.body()
    if len(raw) > MAX_BODY_BYTES:
        raise HTTPException(413, too_big)
    try:
        body = json.loads(raw or b"{}")
    except (ValueError, UnicodeDecodeError):
        raise HTTPException(400, "Expected a JSON body.") from None
    if not isinstance(body, dict):
        raise HTTPException(400, "Expected a JSON object.")
    return body


def _page(value) -> str:
    page = str(value or "")[:200]
    return "/" + page if page and not page.startswith("/") else page


async def submit(request: Request) -> dict:
    """POST /api/feedback: validate, rate-limit, store, forward in the background."""
    body = await _json_body(request, "The message is too large. Attach fewer or smaller images.")
    if str(body.get("website") or "").strip():          # honeypot: pretend success, store nothing
        return {"ok": True, "delivered": True, "note": ""}

    message = str(body.get("message") or "").strip()
    if len(message) < MIN_MESSAGE:
        raise HTTPException(400, "Please write a short message.")
    if len(message) > MAX_MESSAGE:
        raise HTTPException(400, f"Please keep the message under {MAX_MESSAGE} characters.")
    page = _page(body.get("page"))
    images = _decode_images(body.get("images"))

    wait = _rate_limited(request.client.host if request.client else "unknown")     # only valid messages use up the sender's budget
    if wait:
        raise HTTPException(429, wait)

    try:
        ref = reports.add_feedback(message, page, request.headers.get("user-agent") or "", images)
    except Exception:  # noqa: BLE001
        log.exception("feedback: could not store the message")
        raise HTTPException(503, "Could not save your message right now. Please try again later.") from None
    reports.forward_async(ref)
    return {"ok": True, "delivered": True, "note": "", "ref": ref}


_CRASH_HITS: dict[str, list[float]] = {}
CRASH_HOURLY_LIMIT = 30


def _crash_budget(ip: str) -> bool:
    now = time.time()
    with _RATE_LOCK:
        hits = [t for t in _CRASH_HITS.get(ip, []) if now - t < HOUR_S]
        if len(hits) >= CRASH_HOURLY_LIMIT:
            _CRASH_HITS[ip] = hits
            return False
        hits.append(now)
        _CRASH_HITS[ip] = hits
        if len(_CRASH_HITS) > 4096:
            for key in [k for k, v in _CRASH_HITS.items() if not v or v[-1] < now - HOUR_S]:
                del _CRASH_HITS[key]
    return True


async def submit_crash(request: Request) -> dict:
    """POST /api/crash: a page that failed to draw in the browser. No form data is sent: the error, its stack and the page path."""
    body = await _json_body(request, "Too large.")
    message = str(body.get("message") or "").strip()[:300]
    if not message:
        raise HTTPException(400, "Missing error message.")
    if not _crash_budget(request.client.host if request.client else "unknown"):
        raise HTTPException(429, "Too many reports from this address.")
    stack = str(body.get("stack") or "")[:6000]
    try:
        ref, forward = reports.add_browser_crash(message, stack, _page(body.get("page")), request.headers.get("user-agent") or "")
    except Exception:  # noqa: BLE001
        log.exception("crash: could not store the report")
        raise HTTPException(503, "Could not save the report.") from None
    if forward:
        reports.forward_async(ref)
    return {"ok": True, "ref": ref}
