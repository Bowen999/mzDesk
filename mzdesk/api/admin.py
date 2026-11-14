"""The maintainer's password for the stats and reports pages, checked on the server with a limit on wrong guesses.

Password hash: PBKDF2-HMAC-SHA256, 200 000 rounds. A server sets its own in MZDESK_ADMIN_SALT / MZDESK_ADMIN_HASH (hex), so the
hash that ships to browsers in web/src/lock.ts (the page lock) is not the one that guards the data. Without them the page-lock
hash is used (local version).
    python -m mzdesk.api.admin            # asks for a password, prints the two lines for /etc/mzdesk/mzdesk.env
    deploy/set-admin-password.sh          # the same, and installs them on the server
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


MIN_LENGTH = 1           # the maintainer's choice; online guesses are limited (8 per address per 15 minutes) and the hash stays on the server


def main(argv=None) -> None:
    """Print the two MZDESK_ADMIN_* lines for a password typed twice (hidden), or read from stdin with --stdin (deploy/set-admin-password.sh)."""
    import getpass
    import sys
    argv = sys.argv[1:] if argv is None else argv
    if "--stdin" in argv:
        password = sys.stdin.read().rstrip("\r\n")
    else:
        password = getpass.getpass("New admin password: ")
        if password != getpass.getpass("Again: "):
            raise SystemExit("The two passwords differ.")
    if len(password) < MIN_LENGTH:
        raise SystemExit(f"Use at least {MIN_LENGTH} characters.")
    if len(password) < 8:
        print("note: under 8 characters can be guessed by many machines working together; 12 or more is safe.", file=sys.stderr)
    elif len(password) < 12:
        print("note: 12 characters or more is safer.", file=sys.stderr)
    salt, digest = hash_password(password)
    print(f"MZDESK_ADMIN_SALT={salt}\nMZDESK_ADMIN_HASH={digest}")


if __name__ == "__main__":
    main()
