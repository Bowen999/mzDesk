"""Bound the MassQL JSON request even when the app runs without nginx."""
from fastapi.responses import JSONResponse

MASSQL_REQUEST_LIMIT = 27_000_000  # 20 MB raw file, base64-encoded, plus form fields


class MassQLRequestLimit:
    def __init__(self, app):
        self.app = app

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http" or scope.get("method") != "POST" or scope.get("path") != "/api/modules/gnps_massql/run":
            return await self.app(scope, receive, send)
        headers = dict(scope.get("headers", []))
        try:
            length = int(headers.get(b"content-length", b"0"))
        except ValueError:
            length = MASSQL_REQUEST_LIMIT + 1
        if length > MASSQL_REQUEST_LIMIT:
            return await JSONResponse({"detail": "The MassQL request is too large; upload one file of at most 20 MB."}, status_code=413)(scope, receive, send)
        messages, size = [], 0
        while True:
            message = await receive()
            if message["type"] == "http.disconnect":
                return
            size += len(message.get("body", b""))
            if size > MASSQL_REQUEST_LIMIT:
                return await JSONResponse({"detail": "The MassQL request is too large; upload one file of at most 20 MB."}, status_code=413)(scope, receive, send)
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
