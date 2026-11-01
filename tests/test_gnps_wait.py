"""fastMASST waiting (compound_db_check.fasst_search): max_wait_s is a deadline on the clock, not a count of sleeps.

The old loop added only the sleeps to its `waited` counter, so a status request that took 30 s never counted: max_wait_s=300 with 30 s requests
waited about 3000 s. Time here is virtual (a fake clock the sleeps and the slow requests advance), so the tests take no real time.
"""
from types import SimpleNamespace

import pytest

from mzdesk.gnps import compound_db_check as g

PEAKS = [(100.0, 100.0), (110.0, 20.0)]
REQUEST_S = 30.0


class Clock:
    def __init__(self):
        self.now = 1000.0

    def sleep(self, seconds):
        self.now += seconds

    def monotonic(self):
        return self.now

    time = monotonic


@pytest.fixture
def clock(monkeypatch):
    c = Clock()
    monkeypatch.setattr(g, "time", SimpleNamespace(sleep=c.sleep, monotonic=c.monotonic, time=c.time))
    monkeypatch.setattr(g, "post_json", lambda url, payload, timeout=60: {"id": "task-1"})
    return c


def _status_costs(clock, monkeypatch, reply, seconds=REQUEST_S):
    """Every get_json call takes `seconds` and returns reply(call number)."""
    calls = []

    def get_json(url, timeout=40):
        calls.append(url)
        clock.now += seconds
        return reply(len(calls))
    monkeypatch.setattr(g, "get_json", get_json)
    return calls


def test_slow_status_requests_count_toward_max_wait(clock, monkeypatch):
    calls = _status_costs(clock, monkeypatch, lambda n: {"status": "RUNNING"})
    started = clock.now
    out = g.fasst_search(PEAKS, 300.0, g.REF_LIB, max_wait_s=300)
    waited = clock.now - started
    assert "timeout after" in out["error"] and out["task_id"] == "task-1"
    assert waited <= 300 + REQUEST_S, f"waited {waited:.0f} s for max_wait_s=300 (one request is {REQUEST_S:.0f} s)"
    assert waited >= 300                                         # and it did use the whole allowance
    assert len(calls) <= 10                                      # the old loop made 38 requests (about 1100 s)


def test_the_timeout_message_gives_the_time_really_waited(clock, monkeypatch):
    _status_costs(clock, monkeypatch, lambda n: {"status": "RUNNING"})
    started = clock.now
    out = g.fasst_search(PEAKS, 300.0, g.REF_LIB, max_wait_s=300)
    said = int(out["error"].split("timeout after ")[1].split("s")[0])
    assert said == pytest.approx(clock.now - started, abs=1) and said >= 300


def test_elapsed_in_the_run_info_is_the_whole_run(clock, monkeypatch):
    _status_costs(clock, monkeypatch, lambda n: {"status": "RUNNING"})
    started = clock.now
    _, info = g.fasst_search(PEAKS, 300.0, g.REF_LIB, max_wait_s=60, return_raw=True)
    assert info["elapsed_s"] == pytest.approx(clock.now - started, abs=0.1) and info["task_id"] == "task-1"


def test_fast_requests_still_wait_the_whole_allowance(clock, monkeypatch):
    calls = _status_costs(clock, monkeypatch, lambda n: {"status": "RUNNING"}, seconds=0.2)
    started = clock.now
    out = g.fasst_search(PEAKS, 300.0, g.REF_LIB, max_wait_s=80)
    assert "timeout after" in out["error"] and 80 <= clock.now - started <= 80.5
    assert len(calls) >= 9                                       # polled about every 8 s


def test_a_search_that_finishes_after_slow_requests_is_returned(clock, monkeypatch):
    raw = [{"USI": "mzspec:X:1", "Dataset": "MSV1", "Cosine": 0.9, "Matching Peaks": 5}]

    def reply(n):                                                # RUNNING three times, then COMPLETED, then the result
        return {"status": "RUNNING"} if n < 4 else {"status": "COMPLETED"} if n == 4 else {"results": raw}
    _status_costs(clock, monkeypatch, reply)
    out = g.fasst_search(PEAKS, 300.0, g.REF_LIB, max_wait_s=300)
    assert out["n_matches"] == 1 and "error" not in out


def test_a_failed_status_still_ends_the_wait_at_once(clock, monkeypatch):
    calls = _status_costs(clock, monkeypatch, lambda n: {"status": "FAILED"}, seconds=1.0)
    out = g.fasst_search(PEAKS, 300.0, g.REF_LIB, max_wait_s=300)
    assert out["error"] == {"status": "FAILED"} and len(calls) == 1
