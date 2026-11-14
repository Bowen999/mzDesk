"""Update notice of the local version: is a newer mzDesk published?

Fetches MZDESK_UPDATE_URL (default $MZDESK_BASE_URL/version.json, which is https://mzdesk.com/version.json) in a background thread, with a
3 s timeout, at most once per 24 h, and keeps the answer in ~/.mzdesk/cache/update.json. The request is a plain GET with
`User-Agent: mzdesk/<version>` and nothing else: no user data, no ids. A failed check (offline) is tried again after an hour; it never raises
and never delays the app. Off on mzdesk.com (MZDESK_HOSTED=1) and with MZDESK_NO_UPDATE_CHECK=1.
`mzdesk update` (mzdesk/cli.py) installs the newer version by running the installer again.
"""
import json
import os
import re
import threading
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from mzdesk import __version__

DEFAULT_BASE_URL = "https://mzdesk.com"
CHANGELOG = "https://mzdesk.com/docs/changelog"
COMMAND = "mzdesk update"
TIMEOUT_S = 3
CHECK_EVERY_S = 24 * 3600
RETRY_AFTER_FAILURE_S = 3600
MAX_BYTES = 65536
FIRST_ANSWER_WAIT_S = 1.5          # GET /api/update waits this long for a check that is already running and has nothing cached yet

_lock = threading.Lock()
_thread = None


def base_url():
    return (os.environ.get("MZDESK_BASE_URL") or DEFAULT_BASE_URL).rstrip("/")


def update_url():
    return os.environ.get("MZDESK_UPDATE_URL") or base_url() + "/version.json"


def enabled():
    return os.environ.get("MZDESK_HOSTED") != "1" and os.environ.get("MZDESK_NO_UPDATE_CHECK") != "1"


def cache_file():
    return Path.home() / ".mzdesk" / "cache" / "update.json"


def parse_version(text):
    """'0.5.10' -> (0, 5, 10); None for anything that is not dotted numbers. Compared as tuples, so 0.5.10 > 0.5.9."""
    if not isinstance(text, str) or not re.fullmatch(r"\d{1,6}(\.\d{1,6}){0,3}", text.strip()):
        return None
    return tuple(int(part) for part in text.strip().split("."))


def is_newer(latest, current=None):
    a, b = parse_version(latest), parse_version(current or __version__)
    return a is not None and b is not None and a > b


def _now():
    return time.time()


def _iso(stamp):
    return datetime.fromtimestamp(stamp, timezone.utc).isoformat(timespec="seconds")


def read_cache():
    """The cached answer {"latest", "released", "notes", "checked_at" (epoch), "tried_at" (epoch)} or {} (missing, damaged, wrong shape)."""
    try:
        data = json.loads(cache_file().read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def _write_cache(data):
    try:
        path = cache_file()
        path.parent.mkdir(parents=True, exist_ok=True)
        part = path.with_name(path.name + ".part")
        part.write_text(json.dumps(data), encoding="utf-8")
        part.replace(path)
    except OSError:
        pass                                   # the notice is a convenience; a read-only home must not matter


def fetch(timeout=TIMEOUT_S):
    """GET the version file -> {"latest", "released", "notes"}; None when it cannot be read (offline, bad answer, not a newer-style file)."""
    request = urllib.request.Request(update_url(), headers={"User-Agent": f"mzdesk/{__version__}", "Accept": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as resp:
            body = json.loads(resp.read(MAX_BYTES).decode("utf-8"))
    except (urllib.error.URLError, OSError, ValueError, TimeoutError):
        return None
    latest = body.get("version") if isinstance(body, dict) else None
    if parse_version(latest) is None:
        return None
    notes = body.get("notes")
    released = body.get("released")
    return {"latest": latest.strip(), "released": released if isinstance(released, str) and re.fullmatch(r"\d{4}-\d{2}-\d{2}", released) else None,
            "notes": notes if isinstance(notes, str) and notes.startswith(("https://", "http://")) else CHANGELOG}


def refresh(timeout=TIMEOUT_S):
    """Fetch now and write the cache (a failure is remembered as tried_at, so it is not retried before an hour has passed). Returns the cache."""
    data = read_cache()
    got = fetch(timeout)
    if got:
        data = {**got, "checked_at": _now(), "tried_at": _now()}
    else:
        data = {**data, "tried_at": _now()}
    _write_cache(data)
    return data


def _num(value):
    return value if isinstance(value, (int, float)) and not isinstance(value, bool) else 0


def _due(data):
    """A check is due when the last good answer is older than 24 h and the last attempt (good or not) older than an hour."""
    now = _now()
    return now - _num(data.get("checked_at")) >= CHECK_EVERY_S and now - _num(data.get("tried_at")) >= RETRY_AFTER_FAILURE_S


def check_in_background():
    """Start the check in a daemon thread when one is due and none is running. Returns the thread or None."""
    global _thread
    if not enabled():
        return None
    with _lock:
        if (_thread is not None and _thread.is_alive()) or not _due(read_cache()):
            return None
        _thread = threading.Thread(target=refresh, daemon=True, name="update-check")
        _thread.start()
        return _thread


def info():
    """The body of GET /api/update. Never raises. Offline, disabled or nothing cached yet: latest null, newer false."""
    off = {"current": __version__, "latest": None, "newer": False, "released": None, "notes": CHANGELOG, "command": COMMAND, "checked_at": None}
    if not enabled():
        return off
    try:
        thread = check_in_background() or _thread
        data = read_cache()
        if not data.get("latest") and thread is not None and thread.is_alive():
            thread.join(FIRST_ANSWER_WAIT_S)             # the first page load of a fresh install gets the answer too, if it is quick
            data = read_cache()
        latest = data.get("latest") if parse_version(data.get("latest")) else None
        checked = data.get("checked_at")
        return {**off, "latest": latest, "newer": bool(latest and is_newer(latest)), "released": data.get("released") if latest else None,
                "notes": (data.get("notes") or CHANGELOG) if latest else CHANGELOG, "checked_at": _iso(checked) if _num(checked) else None}
    except Exception:  # noqa: BLE001 - the promise is: never raises
        return off


def terminal_line():
    """The line `mzdesk` prints at start when the cache holds a newer version, else None. Reads the cache only; no network."""
    if not enabled():
        return None
    latest = read_cache().get("latest")
    if parse_version(latest) and is_newer(latest):
        return f"mzDesk {latest} is available (you have {__version__}). Update: {COMMAND}"
    return None
