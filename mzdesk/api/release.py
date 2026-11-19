"""mzdesk.com only (MZDESK_HOSTED=1): /install.sh, /install.ps1, /version.json and /download/<wheel>, so the local version installs with one command.

The scripts are the files in mzdesk/install/. version.json, the wheels and their .sha256 files are read from MZDESK_RELEASES_DIR
(default /var/lib/mzdesk/releases) at each request: a deploy only has to put new files there. On the local version
(MZDESK_HOSTED unset) every route answers 404.
"""
import os
import re
from pathlib import Path

from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse

HOSTED = os.environ.get("MZDESK_HOSTED") == "1"          # same convention as app.py: read when the server starts
INSTALL_DIR = Path(__file__).resolve().parents[1] / "install"
DEFAULT_RELEASES_DIR = "/var/lib/mzdesk/releases"
# The only names /download/ serves: the wheel of a release and its checksum file (written by deploy/build-release.sh).
# [0-9], not \d: \d also matches digits of other scripts. fullmatch, not $: $ accepts a trailing newline.
DOWNLOAD_NAME = re.compile(r"mzdesk-[0-9]+\.[0-9]+\.[0-9]+-py3-none-any\.whl(\.sha256)?")
NO_CACHE = {"Cache-Control": "no-cache"}

router = APIRouter()


def releases_dir() -> Path:
    return Path(os.environ.get("MZDESK_RELEASES_DIR") or DEFAULT_RELEASES_DIR)


def _file(path: Path, media_type: str, headers=None, filename=None) -> FileResponse:
    if not HOSTED or not path.is_file():
        raise HTTPException(404, "not found")
    return FileResponse(path, media_type=media_type, headers={**NO_CACHE, **(headers or {})}, filename=filename)


@router.api_route("/install.sh", methods=["GET", "HEAD"], include_in_schema=False)
def install_sh():
    return _file(INSTALL_DIR / "install.sh", "text/plain")


@router.api_route("/install.ps1", methods=["GET", "HEAD"], include_in_schema=False)
def install_ps1():
    return _file(INSTALL_DIR / "install.ps1", "text/plain")


@router.api_route("/version.json", methods=["GET", "HEAD"], include_in_schema=False)
def version_json():
    return _file(releases_dir() / "version.json", "application/json")


@router.api_route("/download/{name}", methods=["GET", "HEAD"], include_in_schema=False)
def download(name: str):
    if not DOWNLOAD_NAME.fullmatch(name):
        raise HTTPException(404, "not found")
    if name.endswith(".sha256"):
        return _file(releases_dir() / name, "text/plain")
    # A wheel is a zip file, already compressed: gzip would only cost CPU and drop Content-Length. GZipMiddleware (app.py) leaves
    # a response alone when it already names a Content-Encoding, so say "identity" (no encoding) and the bytes pass unchanged.
    return _file(releases_dir() / name, "application/octet-stream", headers={"Content-Encoding": "identity"}, filename=name)
