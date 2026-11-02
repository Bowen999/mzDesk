"""The maintainer's password for the stats and reports pages, checked on the server with a limit on wrong guesses.

Password hash: PBKDF2-HMAC-SHA256, 200 000 rounds. A server sets its own in MZDESK_ADMIN_SALT / MZDESK_ADMIN_HASH (hex), so the
hash that ships to browsers in web/src/lock.ts (the page lock) is not the one that guards the data. Without them the page-lock
hash is used (local version).
    python -m mzdesk.api.admin            # asks for a password, prints the two lines for /etc/mzdesk/mzdesk.env
"""
import hashlib
import secrets
import threading
import time

from fastapi import HTTPException, Request

ITERATIONS = 200_000
WINDOW_S, MAX_WRONG = 900.0, 8          # 8 wrong passwords per address in 15 minutes, then a pause until the oldest has aged out
_lock = threading.Lock()
_wrong: dict[str, list[float]] = {}


def hash_password(password: str, salt_hex: str | None = None) -> tuple[str, str]:
    """-> (salt_hex, hash_hex)."""
    salt = bytes.fromhex(salt_hex) if salt_hex else secrets.token_bytes(16)
    return salt.hex(), hashlib.pbkdf2_hmac("sha256", password.encode(), salt, ITERATIONS).hex()


def client_ip(request: Request) -> str:
    return request.client.host if request.client else "unknown"


def require(password, request: Request) -> None:
    """Return when the password is right; 401 when it is wrong; 429 after too many wrong ones from the same address."""
    from mzdesk.api import stats          # late: tests replace stats.check_password
    ip, now = client_ip(request), time.time()
    with _lock:
        recent = [t for t in _wrong.get(ip, []) if now - t < WINDOW_S]
        _wrong[ip] = recent
        if len(recent) >= MAX_WRONG:
            raise HTTPException(429, "Too many wrong passwords. Try again in a few minutes.")
    if stats.check_password(password):
        return
    with _lock:
        _wrong.setdefault(ip, []).append(time.time())
        if len(_wrong) > 2048:
            for key in [k for k, v in _wrong.items() if not v or v[-1] < now - WINDOW_S]:
                del _wrong[key]
    time.sleep(0.6)                        # slows guessing; nginx also rate-limits these paths
    raise HTTPException(401, "Wrong password.")


def reset_throttle() -> None:
    with _lock:
        _wrong.clear()


if __name__ == "__main__":
    import getpass
    first = getpass.getpass("New admin password: ")
    if len(first) < 12 or first != getpass.getpass("Again: "):
        raise SystemExit("Use at least 12 characters, typed the same twice.")
    salt, digest = hash_password(first)
    print(f"MZDESK_ADMIN_SALT={salt}\nMZDESK_ADMIN_HASH={digest}")
