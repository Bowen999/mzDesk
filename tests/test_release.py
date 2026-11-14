"""mzdesk.com's release routes (mzdesk/api/release.py): the two install scripts, version.json and the wheel download.
They answer only when MZDESK_HOSTED=1 was set at start; the tests switch that on in the module and run the real app, with its middleware."""
import hashlib
import json

import pytest
from fastapi.testclient import TestClient

from mzdesk.api import release
from mzdesk.api.app import app

WHEEL = "mzdesk-9.8.7-py3-none-any.whl"
client = TestClient(app)


@pytest.fixture
def releases(tmp_path, monkeypatch):
    """A releases folder with one fake wheel (compressible, over GZipMiddleware's size limit), its checksum and version.json."""
    data = b"mzdesk wheel stand-in " * 250_000                     # 5.5 MB; gzip would shrink this to a few KB
    (tmp_path / WHEEL).write_bytes(data)
    digest = hashlib.sha256(data).hexdigest()
    (tmp_path / f"{WHEEL}.sha256").write_text(f"{digest}  {WHEEL}\n")
    info = {"version": "9.8.7", "released": "2026-11-02", "wheel": WHEEL, "sha256": digest, "size": len(data), "python": ">=3.10",
            "notes": "https://mzdesk.com/docs/changelog"}
    (tmp_path / "version.json").write_text(json.dumps(info, indent=2) + "\n")
    (tmp_path / "secret.txt").write_text("not a release file")
    monkeypatch.setenv("MZDESK_RELEASES_DIR", str(tmp_path))
    monkeypatch.setattr(release, "HOSTED", True)
    return {"dir": tmp_path, "data": data, "digest": digest, "info": info}


def test_scripts_are_the_files_of_the_package(releases):
    for name in ("install.sh", "install.ps1"):
        r = client.get(f"/{name}")
        assert r.status_code == 200
        assert r.headers["content-type"].startswith("text/plain")
        assert r.content == (release.INSTALL_DIR / name).read_bytes()
        assert r.headers["x-content-type-options"] == "nosniff"          # the app's middleware


def test_version_json_is_read_at_each_request(releases):
    r = client.get("/version.json")
    assert r.status_code == 200 and r.headers["content-type"].startswith("application/json")
    assert r.json() == releases["info"]
    (releases["dir"] / "version.json").write_text(json.dumps({**releases["info"], "version": "9.8.8"}))
    assert client.get("/version.json").json()["version"] == "9.8.8"


def test_wheel_arrives_unchanged_even_when_the_client_accepts_gzip(releases):
    r = client.get(f"/download/{WHEEL}", headers={"Accept-Encoding": "gzip"})
    assert r.status_code == 200
    assert r.headers["content-type"] == "application/octet-stream"
    assert r.headers.get("content-encoding") != "gzip"                  # GZipMiddleware must not touch it
    assert r.headers["content-length"] == str(len(releases["data"]))    # and the size is known up front (a progress bar)
    assert hashlib.sha256(r.content).hexdigest() == releases["digest"]
    assert f'filename="{WHEEL}"' in r.headers["content-disposition"]
    plain = client.get(f"/download/{WHEEL}")
    assert plain.content == r.content


def test_wheel_checksum_file(releases):
    r = client.get(f"/download/{WHEEL}.sha256")
    assert r.status_code == 200 and r.headers["content-type"].startswith("text/plain")
    assert r.text.split()[0] == releases["digest"]


@pytest.mark.parametrize("path", ["/install.sh", "/install.ps1", "/version.json", f"/download/{WHEEL}", f"/download/{WHEEL}.sha256"])
def test_head_answers_with_headers_and_no_body(releases, path):
    r = client.head(path)
    assert r.status_code == 200 and r.content == b""
    assert int(r.headers["content-length"]) == len(client.get(path).content)


@pytest.mark.parametrize("name", [
    "version.json",                                       # inside the folder, but not a release file
    "secret.txt",
    "mzdesk-9.8.7.tar.gz",
    "mzdesk-9.8.7-py3-none-any.whl.exe",
    "mzdesk-9.8.7-py3-none-any.whl.sha256.bak",
    "mzdesk-9.8-py3-none-any.whl",
    "mzdesk-9.8.7-py3-none-any.whl%0A",                   # a trailing newline
    "mzdesk-%D9%A9.8.7-py3-none-any.whl",                 # an Arabic-Indic digit
    "other-9.8.7-py3-none-any.whl",
    "..%2Fversion.json",
    "%2e%2e%2fversion.json",
    "..%5Cversion.json",
    "mzdesk-9.8.7-py3-none-any.whl/../secret.txt",
    ".mzdesk-9.8.7-py3-none-any.whl",
])
def test_download_refuses_every_other_name(releases, name):
    r = client.get(f"/download/{name}")
    # a name with a slash in it is not a /download/ route at all and falls through to the app's own page (text/html) or a 404
    assert r.status_code == 404 or r.headers["content-type"].startswith("text/html")
    assert b"not a release file" not in r.content and b'"sha256"' not in r.content and b"stand-in" not in r.content


def test_missing_files_are_404(tmp_path, monkeypatch):
    monkeypatch.setenv("MZDESK_RELEASES_DIR", str(tmp_path))                  # empty: nothing released yet
    monkeypatch.setattr(release, "HOSTED", True)
    assert client.get("/version.json").status_code == 404
    assert client.get(f"/download/{WHEEL}").status_code == 404
    assert client.get("/install.sh").status_code == 200                       # the scripts do not depend on a release


@pytest.mark.parametrize("path", ["/install.sh", "/install.ps1", "/version.json", f"/download/{WHEEL}", f"/download/{WHEEL}.sha256"])
def test_local_version_has_none_of_it(releases, monkeypatch, path):
    monkeypatch.setattr(release, "HOSTED", False)
    assert client.get(path).status_code == 404 and client.head(path).status_code == 404


def test_releases_folder_comes_from_the_environment_at_request_time(tmp_path, monkeypatch):
    monkeypatch.delenv("MZDESK_RELEASES_DIR", raising=False)
    assert release.releases_dir().as_posix() == "/var/lib/mzdesk/releases"        # the server folder (\\var\\lib\\... when run on Windows)
    monkeypatch.setenv("MZDESK_RELEASES_DIR", str(tmp_path))
    assert release.releases_dir() == tmp_path
