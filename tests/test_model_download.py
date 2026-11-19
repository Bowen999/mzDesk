"""ms2_similarity._download_file: a stuck server must not block forever, two downloads must not share a file, a short file is never kept.

A throw-away HTTP server on 127.0.0.1 plays the model host (nothing leaves the machine).
"""
import os
import socket
import threading
import time

import pytest

from mzdesk.ms2 import ms2_similarity as sim

BODY = bytes(range(256)) * 8192                                    # 2 MiB


class Host:
    """One-connection-per-request HTTP/1.0-style server. mode: ok | slow | stall | no_headers | truncated."""

    def __init__(self, mode="ok", body=BODY):
        self.mode, self.body = mode, body
        self.release = threading.Event()                           # lets stalled handlers go at teardown
        self.requests = 0
        self.sock = socket.socket()
        self.sock.bind(("127.0.0.1", 0))
        self.sock.listen(8)
        self.url = f"http://127.0.0.1:{self.sock.getsockname()[1]}/model.bin"
        threading.Thread(target=self._accept, daemon=True).start()

    def _accept(self):
        while True:
            try:
                conn, _ = self.sock.accept()
            except OSError:
                return
            self.requests += 1
            threading.Thread(target=self._handle, args=(conn,), daemon=True).start()

    def _handle(self, conn):
        try:
            conn.recv(65536)
            if self.mode == "no_headers":
                self.release.wait(60)
                return
            n = len(self.body)
            conn.sendall(f"HTTP/1.1 200 OK\r\nContent-Length: {n}\r\nConnection: close\r\n\r\n".encode())
            if self.mode == "stall":                               # headers and a part of the data, then silence
                conn.sendall(self.body[:1000])
                self.release.wait(60)
            elif self.mode == "truncated":                         # promises n bytes, sends 1000, closes
                conn.sendall(self.body[:1000])
            elif self.mode == "slow":                              # the whole file, in pieces, so two downloads overlap
                for i in range(0, n, 1 << 16):
                    conn.sendall(self.body[i:i + (1 << 16)])
                    time.sleep(0.02)
            else:
                conn.sendall(self.body)
        except OSError:
            pass
        finally:
            conn.close()

    def close(self):
        self.release.set()
        self.sock.close()


@pytest.fixture
def host(monkeypatch):
    for key in ("http_proxy", "https_proxy", "HTTP_PROXY", "HTTPS_PROXY", "all_proxy", "ALL_PROXY"):
        monkeypatch.delenv(key, raising=False)
    monkeypatch.setenv("no_proxy", "127.0.0.1")
    made = []

    def make(mode="ok", body=BODY):
        made.append(Host(mode, body))
        return made[-1]
    yield make
    for h in made:
        h.close()


def test_a_download_gives_the_whole_file_and_leaves_nothing_else(host, tmp_path):
    dest = tmp_path / "model.bin"
    assert sim._download_file(host("ok").url, str(dest)) is True
    assert dest.read_bytes() == BODY and os.listdir(tmp_path) == ["model.bin"]


@pytest.mark.skipif(os.name == "nt", reason="POSIX permission bits")
@pytest.mark.parametrize("umask, mode", [(0o022, 0o644), (0o077, 0o600), (0o002, 0o664)])
def test_the_downloaded_model_gets_the_permissions_of_any_new_file(host, tmp_path, umask, mode):
    """0o666 minus the umask (usually 0644), not always the 0600 of a mkstemp file."""
    url, dest = host("ok").url, tmp_path / "model.bin"
    old = os.umask(umask)
    try:
        assert sim._download_file(url, str(dest)) is True
    finally:
        os.umask(old)
    assert dest.stat().st_mode & 0o777 == mode


@pytest.mark.parametrize("mode", ["stall", "no_headers"])
def test_a_stuck_server_ends_in_failure_after_the_timeout_not_never(host, tmp_path, monkeypatch, mode):
    """Headers and 1000 bytes (or nothing at all) and then silence: the old urlretrieve had no timeout and the thread waited for ever."""
    monkeypatch.setattr(sim, "DOWNLOAD_TIMEOUT_S", 1)
    dest = tmp_path / "model.bin"
    started = time.monotonic()
    done = []
    worker = threading.Thread(target=lambda: done.append(sim._download_file(host(mode).url, str(dest))), daemon=True)
    worker.start()
    worker.join(10)
    assert not worker.is_alive(), "the download is still blocked after 10 s with a 1 s timeout"
    assert done == [False] and time.monotonic() - started < 6
    assert os.listdir(tmp_path) == []                              # no model, no temporary file


def test_a_file_shorter_than_its_content_length_is_not_kept(host, tmp_path):
    dest = tmp_path / "model.bin"
    assert sim._download_file(host("truncated").url, str(dest)) is False
    assert os.listdir(tmp_path) == []


def test_an_empty_file_is_not_kept(host, tmp_path):
    dest = tmp_path / "model.bin"
    assert sim._download_file(host("ok", body=b"").url, str(dest)) is False
    assert os.listdir(tmp_path) == []


def test_a_failed_download_leaves_the_cached_model_alone(host, tmp_path):
    dest = tmp_path / "model.bin"
    dest.write_bytes(b"the good model")
    assert sim._download_file(host("truncated").url, str(dest)) is False
    assert dest.read_bytes() == b"the good model" and os.listdir(tmp_path) == ["model.bin"]


def test_a_refused_connection_is_a_plain_failure(tmp_path):
    with socket.socket() as s:                                     # a port nobody listens on
        s.bind(("127.0.0.1", 0))
        port = s.getsockname()[1]
    assert sim._download_file(f"http://127.0.0.1:{port}/model.bin", str(tmp_path / "model.bin")) is False
    assert os.listdir(tmp_path) == []


def test_two_first_time_downloads_at_once_do_not_overwrite_each_other(host, tmp_path):
    """Both used to write dest + '.part': one os.replace then found its file gone (or the file held both writers' bytes)."""
    url, dest = host("slow").url, str(tmp_path / "model.bin")
    results = []
    threads = [threading.Thread(target=lambda: results.append(sim._download_file(url, dest))) for _ in range(2)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(30)
    assert results == [True, True]
    assert open(dest, "rb").read() == BODY and os.listdir(tmp_path) == ["model.bin"]
