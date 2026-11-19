"""Security headers for every response, and a scrub for error text that is about to leave the server.

The headers live in the app (not only in nginx) so they hold behind any proxy and in the local version too. HSTS is the
exception: it is a transport promise, so nginx sets it (deploy/nginx-mzdesk.conf).

Content-Security-Policy: scripts only from this origin plus the inline scripts of index.html (allowed by hash, computed
from the built file, so a rebuild needs no edit here); nothing is allowed to embed the app; no plugins, no <base> tricks;
requests only to this origin. Styles may be inline (React style attributes, the opening splash). 'unsafe-eval' stays on because
the structure editor (Ketcher, PaperScript) builds functions with new Function; the app itself has no eval and renders no
untrusted HTML.
"""
import base64
import hashlib
import os
import re
from pathlib import Path

from starlette.datastructures import MutableHeaders

STATIC_INDEX = Path(__file__).resolve().parents[1] / "static" / "index.html"

HEADERS = {
    "X-Content-Type-Options": "nosniff",
    "X-Frame-Options": "SAMEORIGIN",
    "Referrer-Policy": "strict-origin-when-cross-origin",
    "Permissions-Policy": "camera=(), microphone=(), geolocation=(), payment=(), usb=(), interest-cohort=()",
    "Cross-Origin-Opener-Policy": "same-origin",
}

_INLINE_SCRIPT = re.compile(r"<script(?![^>]*\bsrc=)[^>]*>(.*?)</script>", re.S | re.I)


def inline_script_hashes(html: str) -> list[str]:
    """CSP source tokens ('sha256-…') for the inline <script> blocks of a page."""
    return [f"'sha256-{base64.b64encode(hashlib.sha256(body.encode()).digest()).decode()}'" for body in _INLINE_SCRIPT.findall(html) if body.strip()]


def content_security_policy(index_html: str = "") -> str:
    scripts = " ".join(["'self'", *inline_script_hashes(index_html)])
    return "; ".join([
        "default-src 'self'",
        f"script-src {scripts} 'unsafe-eval' 'wasm-unsafe-eval'",
        "style-src 'self' 'unsafe-inline'",
        "img-src 'self' data: blob:",
        "font-src 'self' data:",
        "connect-src 'self'",
        "worker-src 'self' blob:",
        "object-src 'none'",
        "base-uri 'self'",
        "form-action 'self'",
        "frame-ancestors 'self'",
    ])


def _read_index() -> str:
    try:
        return STATIC_INDEX.read_text(encoding="utf-8")
    except OSError:
        return ""


class SecurityHeaders:
    """ASGI middleware: adds the headers above to every HTTP response. A header the route set itself is kept; API answers are not cached."""

    def __init__(self, app):
        self.app = app
        self.headers = {**HEADERS, "Content-Security-Policy": content_security_policy(_read_index())}

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http":
            return await self.app(scope, receive, send)
        api = scope.get("path", "").startswith("/api/")

        async def add_headers(message):
            if message["type"] == "http.response.start":
                headers = MutableHeaders(scope=message)
                for name, value in self.headers.items():
                    headers.setdefault(name, value)
                if api:
                    headers.setdefault("Cache-Control", "no-store")
            await send(message)

        return await self.app(scope, receive, add_headers)


# ---------------------------------------------------------------- error text
_ROOTS = "opt|var|usr|home|root|Users|tmp|private|etc|srv|mnt|Library|Applications|proc|sys|nix|snap"
_PATH = re.compile(rf"(?<![\w.])/(?:{_ROOTS})/(?:[\w.@+~=,-]+/)*([\w.@+~=,-]+)")        # /opt/mzdesk/x/y.py -> y.py
_WINDOWS_PATH = re.compile(r"\b[A-Za-z]:\\(?:[^\\\s'\"]+\\)+([^\\\s'\"]+)")


def strip_paths(text) -> str:
    """The text with the server's folders removed: only the file names stay."""
    return _WINDOWS_PATH.sub(r"\1", _PATH.sub(r"\1", str(text)))


def scrub(text, limit: int = 240) -> str:
    """Error text for a user: no folders, one line, at most `limit` characters."""
    out = " ".join(strip_paths(text).split())
    return out if len(out) <= limit else out[: limit - 1] + "…"


def api_docs_enabled() -> bool:
    """FastAPI's /docs, /redoc and /openapi.json stay off unless MZDESK_API_DOCS=1 (local development)."""
    return os.environ.get("MZDESK_API_DOCS") == "1"
