"""Usage statistics: page views and module runs, kept in a small SQLite file, readable only with the site password.

What is stored per event: time, kind ("view" | "run"), page path, module id, HTTP status, run time (ms), browser family,
referrer host, and a visitor id. No IP address is stored: the visitor id is a hash of (secret, day, IP, user agent), so the
same person counts once per day and cannot be followed from one day to the next.

Where: MZDESK_STATS_DB, default ~/.mzdesk/stats.sqlite. MZDESK_STATS=0 turns recording off.
Password: checked here on the server so the numbers never reach a browser without it (api/admin.py). A server sets its own
(MZDESK_ADMIN_SALT / MZDESK_ADMIN_HASH); without them it is the page lock's (web/src/lock.ts holds the same salt / hash; keep both in sync).
"""
import hashlib
import hmac
import os
import re
import secrets
import sqlite3
import threading
import time
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from urllib.parse import urlparse

# same values as web/src/lock.ts
SALT_HEX = "7f4282fbfe80e0a447687b7fd2817b09"
HASH_HEX = "9bcbf15d44d26180b6d179732218c949367dae9a2e5be513211ce742c6d343b0"
ITERATIONS = 200_000

_lock = threading.Lock()
_con = None
_secret = None


def enabled():
    return os.environ.get("MZDESK_STATS", "1") != "0"


def db_path():
    return Path(os.environ.get("MZDESK_STATS_DB") or Path.home() / ".mzdesk" / "stats.sqlite")


def _db():
    """Shared connection (opened once; WAL so the page can read while runs write)."""
    global _con, _secret
    if _con is None:
        p = db_path()
        p.parent.mkdir(parents=True, exist_ok=True)
        con = sqlite3.connect(p, check_same_thread=False, timeout=5)
        con.execute("PRAGMA journal_mode=WAL")
        con.execute("""CREATE TABLE IF NOT EXISTS events (t REAL, day TEXT, kind TEXT, path TEXT, module TEXT, status INTEGER,
                       ms INTEGER, visitor TEXT, browser TEXT, ref TEXT)""")
        con.execute("CREATE INDEX IF NOT EXISTS events_day ON events(day)")
        con.execute("CREATE TABLE IF NOT EXISTS meta (k TEXT PRIMARY KEY, v TEXT)")
        row = con.execute("SELECT v FROM meta WHERE k='secret'").fetchone()
        if row is None:
            con.execute("INSERT INTO meta VALUES ('secret', ?)", (secrets.token_hex(16),))
            con.commit()
            row = con.execute("SELECT v FROM meta WHERE k='secret'").fetchone()
        _secret = row[0]
        _con = con
    return _con


def _browser(ua):
    ua = ua or ""
    if re.search(r"bot|crawl|spider|curl|wget|python|httpx|scan", ua, re.I):
        return "bot"
    for name, pat in (("Edge", r"Edg/"), ("Chrome", r"Chrome/"), ("Firefox", r"Firefox/"), ("Safari", r"Safari/")):
        if re.search(pat, ua):
            return name
    return "other"


def _ref_host(ref):
    try:
        host = urlparse(ref or "").hostname or ""
    except ValueError:
        return ""
    return "" if host in ("", "mzdesk.com", "www.mzdesk.com", "bowenyang.org", "www.bowenyang.org", "localhost", "127.0.0.1") else host[:80]


def record(request, kind, *, path="", module="", status=0, ms=0, ref=""):
    """Add one event. Never raises: statistics must not break a page or a run."""
    if not enabled():
        return
    try:
        ip = request.client.host if request.client else ""
        ua = request.headers.get("user-agent", "")
        day = date.today().isoformat()
        with _lock:
            con = _db()
            visitor = hashlib.sha256(f"{_secret}|{day}|{ip}|{ua}".encode()).hexdigest()[:16]
            con.execute("INSERT INTO events VALUES (?,?,?,?,?,?,?,?,?,?)",
                        (time.time(), day, kind, path[:200], module[:60], int(status), int(ms), visitor, _browser(ua), _ref_host(ref)))
            con.commit()
    except Exception:  # noqa: BLE001
        pass


def _credentials():
    """(salt, hash) in hex: the server's own (MZDESK_ADMIN_SALT / MZDESK_ADMIN_HASH, see api/admin.py) or the page-lock pair."""
    salt, digest = os.environ.get("MZDESK_ADMIN_SALT", "").strip(), os.environ.get("MZDESK_ADMIN_HASH", "").strip()
    return (salt, digest) if salt and digest else (SALT_HEX, HASH_HEX)


def check_password(password):
    salt, digest = _credentials()
    got = hashlib.pbkdf2_hmac("sha256", str(password or "").encode(), bytes.fromhex(salt), ITERATIONS)
    return hmac.compare_digest(got, bytes.fromhex(digest))


def summary(days=30):
    """Everything the stats page shows, for the last `days` days (today included)."""
    days = max(1, min(int(days), 365))
    first = (date.today() - timedelta(days=days - 1)).isoformat()
    with _lock:
        con = _db()
        q = lambda sql, *a: con.execute(sql, a).fetchall()   # noqa: E731
        human = "browser != 'bot'"

        def totals(since):
            v, r, u = q(f"""SELECT SUM(kind='view'), SUM(kind='run'), COUNT(DISTINCT visitor) FROM events WHERE day >= ? AND {human}""", since)[0]
            return {"views": v or 0, "runs": r or 0, "visitors": u or 0}

        today = date.today()
        out = {
            "days": days,
            "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "totals": {"today": totals(today.isoformat()), "7d": totals((today - timedelta(days=6)).isoformat()), "range": totals(first),
                       "all": totals("0000")},
            "daily": [],
            "modules": [{"module": m, "runs": n, "errors": e or 0, "avg_ms": round(a or 0)} for m, n, e, a in
                        q(f"""SELECT module, COUNT(*), SUM(status >= 400), AVG(ms) FROM events WHERE kind='run' AND day >= ? AND {human}
                              GROUP BY module ORDER BY COUNT(*) DESC""", first)],
            "pages": [{"path": p, "views": n} for p, n in
                      q(f"SELECT path, COUNT(*) FROM events WHERE kind='view' AND day >= ? AND {human} GROUP BY path ORDER BY COUNT(*) DESC LIMIT 15", first)],
            "referrers": [{"host": h, "views": n} for h, n in
                          q(f"SELECT ref, COUNT(*) FROM events WHERE kind='view' AND ref != '' AND day >= ? AND {human} GROUP BY ref ORDER BY COUNT(*) DESC LIMIT 10", first)],
            "browsers": [{"browser": b, "visitors": n} for b, n in
                         q("SELECT browser, COUNT(DISTINCT visitor) FROM events WHERE day >= ? GROUP BY browser ORDER BY 2 DESC", first)],
            "errors": [{"time": datetime.fromtimestamp(t, timezone.utc).isoformat(timespec="seconds"), "module": m, "status": s} for t, m, s in
                       q("SELECT t, module, status FROM events WHERE kind='run' AND status >= 400 AND day >= ? ORDER BY t DESC LIMIT 20", first)],
        }
        by_day = {d: (v or 0, r or 0, u or 0) for d, v, r, u in
                  q(f"SELECT day, SUM(kind='view'), SUM(kind='run'), COUNT(DISTINCT visitor) FROM events WHERE day >= ? AND {human} GROUP BY day", first)}
        for i in range(days):
            d = (today - timedelta(days=days - 1 - i)).isoformat()
            v, r, u = by_day.get(d, (0, 0, 0))
            out["daily"].append({"day": d, "views": v, "runs": r, "visitors": u})
    return out
