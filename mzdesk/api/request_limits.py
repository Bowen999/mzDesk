"""Bound the size of every request body, also when the app runs without nginx."""
from fastapi.responses import JSONResponse

DEFAULT_REQUEST_LIMIT = 8_000_000     # pasted spectra, structures and up to three feedback screenshots; nginx caps the same
MASSQL_REQUEST_LIMIT = 27_000_000     # 20 MB raw file, base64-encoded, plus form fields
MASSQL_PATH = "/api/modules/gnps_massql/run"


def limit_for(path: str) -> int:
    return MASSQL_REQUEST_LIMIT if path == MASSQL_PATH else DEFAULT_REQUEST_LIMIT


def too_large(path: str) -> JSONResponse:
    if path == MASSQL_PATH:
        return JSONResponse({"detail": "The MassQL request is too large; upload one file of at most 20 MB."}, status_code=413)
    return JSONResponse({"detail": "The request is too large."}, status_code=413)


class RequestLimit:
    """ASGI middleware: a POST / PUT / PATCH body over the limit gets 413, by Content-Length or, without it, by the bytes received."""

    def __init__(self, app):
        self.app = app

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http" or scope.get("method") not in ("POST", "PUT", "PATCH"):
            return await self.app(scope, receive, send)
        path = scope.get("path", "")
        limit = limit_for(path)
        headers = dict(scope.get("headers", []))
        try:
            length = int(headers.get(b"content-length", b"0"))
        except ValueError:
            length = limit + 1
        if length > limit:
            return await too_large(path)(scope, receive, send)
        messages, size = [], 0
        while True:
            message = await receive()
            if message["type"] == "http.disconnect":
                return
            size += len(message.get("body", b""))
            if size > limit:
                return await too_large(path)(scope, receive, send)
            messages.append(message)
            if not message.get("more_body", False):
                break
        iterator = iter(messages)

        async def replay():
            try:
                return next(iterator)
            except StopIteration:
                return await receive()

        return await self.app(scope, replay, send)


MassQLRequestLimit = RequestLimit      # the name this module had when it only guarded MassQL
