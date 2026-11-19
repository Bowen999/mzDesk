"""Update notice: version compare, the fetch (one plain GET with a User-Agent), the 24 h cache, offline, the switches. A local HTTP server
stands in for mzdesk.com; HOME is a temporary folder. No network."""
import http.server
import json
import threading
import time

import pytest

from mzdesk import __version__, update


@pytest.fixture(autouse=True)
def home(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("USERPROFILE", str(tmp_path))
    for key in ("MZDESK_HOSTED", "MZDESK_NO_UPDATE_CHECK", "MZDESK_UPDATE_URL", "MZDESK_BASE_URL"):
        monkeypatch.delenv(key, raising=False)
    monkeypatch.setattr(update, "_thread", None)
    return tmp_path


@pytest.fixture
def server(monkeypatch):
    """A version.json server: .body (dict, bytes or None = 404), .requests (headers of each GET)."""
    class State:
        body = {"version": "99.0.0", "released": "2099-01-02", "notes": "https://mzdesk.com/docs/changelog"}
        requests = []
        delay = 0.0

    class Handler(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            State.requests.append({"path": self.path, "headers": dict(self.headers)})
            time.sleep(State.delay)
            body = State.body
            if body is None:
                self.send_response(404)
                self.end_headers()
                return
            data = body if isinstance(body, bytes) else json.dumps(body).encode()
            self.send_response(200)
            self.end_headers()
            self.wfile.write(data)

        def log_message(self, *a):
            pass

    State.requests = []
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    monkeypatch.setenv("MZDESK_UPDATE_URL", f"http://127.0.0.1:{srv.server_port}/version.json")
    yield State
    srv.shutdown()


# ---------------------------------------------------------------- compare
def test_versions_compare_as_numbers_not_text():
    assert update.parse_version("0.5.10") > update.parse_version("0.5.9")
    assert update.is_newer("0.5.10", "0.5.9") and update.is_newer("1.0.0", "0.99.99")
    assert not update.is_newer("0.5.3", "0.5.3") and not update.is_newer("0.5.2", "0.5.3") and not update.is_newer("0.5", "0.5.0")
    for bad in (None, "", "abc", "0.5.x", "1..2", "v1.2", 5, "0.5.3-beta", "1.2.3.4.5"):
        assert update.parse_version(bad) is None
        assert not update.is_newer(bad, "0.5.3")


def test_urls_follow_the_environment(monkeypatch):
    assert update.update_url() == "https://mzdesk.com/version.json" and update.base_url() == "https://mzdesk.com"
    monkeypatch.setenv("MZDESK_BASE_URL", "http://localhost:9000/")
    assert update.base_url() == "http://localhost:9000" and update.update_url() == "http://localhost:9000/version.json"
    monkeypatch.setenv("MZDESK_UPDATE_URL", "http://example.test/v.json")
    assert update.update_url() == "http://example.test/v.json"


# ---------------------------------------------------------------- fetch and cache
def test_fetch_sends_one_plain_get_with_the_user_agent(server):
    got = update.fetch()
    assert got == {"latest": "99.0.0", "released": "2099-01-02", "notes": "https://mzdesk.com/docs/changelog"}
    assert len(server.requests) == 1 and server.requests[0]["path"] == "/version.json"
    headers = {k.lower(): v for k, v in server.requests[0]["headers"].items()}
    assert headers["user-agent"] == f"mzdesk/{__version__}"
    assert not {"cookie", "authorization", "referer"} & set(headers) and set(headers) <= {"host", "user-agent", "accept", "accept-encoding", "connection"}


def test_fetch_survives_bad_answers(server):
    for body in (None, b"not json", b"[]", {"version": "abc"}, {"nothing": 1}, b""):
        server.body = body
        assert update.fetch() is None
    server.body = {"version": "1.2.3", "released": "yesterday", "notes": "javascript:alert(1)"}
    assert update.fetch() == {"latest": "1.2.3", "released": None, "notes": update.CHANGELOG}                  # odd extras are dropped


def test_fetch_offline_returns_none(monkeypatch):
    monkeypatch.setenv("MZDESK_UPDATE_URL", "http://127.0.0.1:9/version.json")                                  # nothing listens there
    assert update.fetch(timeout=1) is None


def test_info_checks_once_and_then_uses_the_cache(server, home):
    first = update.info()
    assert first["latest"] == "99.0.0" and first["newer"] is True and first["current"] == __version__ and first["released"] == "2099-01-02"
    assert first["command"] == "mzdesk update" and first["notes"] == "https://mzdesk.com/docs/changelog" and first["checked_at"].endswith("+00:00")
    assert update.cache_file() == home / ".mzdesk/cache/update.json" and update.cache_file().is_file()
    for _ in range(3):
        again = update.info()
        assert again["latest"] == "99.0.0"
    update._thread.join(5) if update._thread else None
    assert len(server.requests) == 1


def test_a_cache_older_than_24_hours_is_checked_again(server):
    update.refresh()
    data = update.read_cache()
    data["checked_at"] = time.time() - 25 * 3600
    data["tried_at"] = time.time() - 25 * 3600
    update._write_cache(data)
    server.body = {"version": "99.1.0"}
    update._thread = None
    body = update.info()
    update._thread.join(5)
    assert update.info()["latest"] == "99.1.0" and body["latest"] in ("99.0.0", "99.1.0")
    assert len(server.requests) == 2


def test_info_waits_briefly_for_the_first_answer(server):
    server.delay = 0.3
    assert update.info()["latest"] == "99.0.0"


def test_info_is_not_slowed_by_a_server_that_hangs(server, monkeypatch):
    server.delay = 5
    monkeypatch.setattr(update, "TIMEOUT_S", 1)
    start = time.monotonic()
    body = update.info()
    assert time.monotonic() - start < update.FIRST_ANSWER_WAIT_S + 1 and body["latest"] is None and body["newer"] is False


def test_offline_means_no_newer_version_and_a_retry_after_an_hour(monkeypatch):
    monkeypatch.setenv("MZDESK_UPDATE_URL", "http://127.0.0.1:9/version.json")
    body = update.info()
    update._thread.join(10)
    assert body["newer"] is False and body["latest"] is None and body["checked_at"] is None
    cache = update.read_cache()
    assert "latest" not in cache and cache["tried_at"] > 0
    assert not update._due(cache)                                                      # tried a moment ago: no new attempt yet
    assert update._due({"tried_at": time.time() - 3700})
    assert not update._due({"checked_at": time.time() - 3600, "tried_at": time.time() - 7200})        # a good answer less than 24 h old


def test_a_failed_check_keeps_the_last_good_answer(server):
    update.refresh()
    server.body = None
    update.refresh()
    cache = update.read_cache()
    assert cache["latest"] == "99.0.0" and cache["checked_at"] > 0
    assert update.info()["latest"] == "99.0.0"


def test_a_damaged_cache_is_ignored(home):
    path = update.cache_file()
    path.parent.mkdir(parents=True)
    for text in ("{", "[]", '"x"', json.dumps({"latest": {"x": 1}, "checked_at": "soon"}), json.dumps({"latest": "99.0.0", "checked_at": True})):
        path.write_text(text)
        body = update.info()
        assert set(body) == {"current", "latest", "newer", "released", "notes", "command", "checked_at"}


def test_the_cache_is_not_written_when_home_is_read_only(server, home, monkeypatch):
    monkeypatch.setattr(update.Path, "mkdir", lambda *a, **k: (_ for _ in ()).throw(OSError("read-only")))
    update.refresh()                                                                   # must not raise
    assert update.read_cache() == {}


# ---------------------------------------------------------------- switches
def test_switches_turn_the_check_off(server, monkeypatch):
    monkeypatch.setenv("MZDESK_NO_UPDATE_CHECK", "1")
    body = update.info()
    assert body["latest"] is None and body["newer"] is False and not server.requests and update.terminal_line() is None
    monkeypatch.delenv("MZDESK_NO_UPDATE_CHECK")
    monkeypatch.setenv("MZDESK_HOSTED", "1")
    assert update.info()["latest"] is None and update.check_in_background() is None and not server.requests


# ---------------------------------------------------------------- terminal line
def test_terminal_line_comes_from_the_cache_only(server):
    assert update.terminal_line() is None and not server.requests                       # nothing cached, no network
    update.refresh()
    assert update.terminal_line() == f"mzDesk 99.0.0 is available (you have {__version__}). Update: mzdesk update"
    update._write_cache({"latest": __version__, "checked_at": time.time()})
    assert update.terminal_line() is None
    update._write_cache({"latest": "0.0.1", "checked_at": time.time()})
    assert update.terminal_line() is None
