"""The terminal view of `mzdesk setup` (mzdesk/local/terminal.py)."""
import io
import re

from mzdesk.local import terminal as tm

ANSI = re.compile(r"\x1b\[[0-9;]*m")


class Utf8(io.StringIO):
    encoding = "utf-8"


def snap(rows, status="running", elapsed=75, **extra):
    return {"status": status, "elapsed_s": elapsed, "components": rows, "log": extra.get("log", []), "error": extra.get("error", "")}


def row(cid, status, **kw):
    return {"id": cid, "title": kw.pop("title", cid.upper()), "status": status, "step": "", "progress": None, "got_mb": None, "total_mb": None, **kw}


def test_sizes_and_bar():
    assert tm.fmt_mb(412) == "412 MB" and tm.fmt_mb(2300) == "2.2 GB"
    assert tm.bar(0.5, 10) == tm.FILL * 5 + tm.EMPTY * 5 and tm.bar(2, 4, fancy=False) == "####" and tm.bar(-1, 4) == tm.EMPTY * 4
    assert tm.step_text(row("a", "running", step="Downloading x.zip", progress=0.4, got_mb=136.0, total_mb=325.0)) == "Downloading x.zip · 136 MB of 325 MB"
    assert tm.step_text(row("a", "running")) == "Working"


def test_a_terminal_gets_colour_a_live_line_and_the_component_lines():
    out = Utf8()
    ui = tm.Terminal(out, color=True, width=100)
    assert ui.live and ui.fancy
    ui.title([row("a", "waiting"), row("b", "waiting")])
    ui.render(snap([row("a", "running", title="MS-FINDER", step="Downloading z.zip", progress=0.25, got_mb=81.0, total_mb=325.0), row("b", "waiting")]))
    ui.render(snap([row("a", "running", title="MS-FINDER", step="Downloading z.zip", progress=0.5, got_mb=162.0, total_mb=325.0), row("b", "waiting")]))
    ui.render(snap([row("a", "done", title="MS-FINDER"), row("b", "running", title="NEIMS")]))
    ui.finish(snap([row("a", "done"), row("b", "done")], status="done"))
    raw = out.getvalue()
    text = ANSI.sub("", raw)
    assert "\x1b[38;2;123;156;240m" in raw                                   # the accent of mzDesk
    assert "[1/2] MS-FINDER" in text and "[1/2] ✓ MS-FINDER" in text and "[2/2] NEIMS" in text
    assert " 50% Downloading z.zip · 162 MB of 325 MB" in text and "All set." in text
    assert text.count("[1/2] MS-FINDER") == 1                                  # a component line is printed once, however often it is rendered
    assert "\r\x1b[2K" in raw                                                  # the live line is redrawn in place


def test_no_colour_when_asked_for_and_never_on_a_file(monkeypatch):
    out = Utf8()
    monkeypatch.setenv("NO_COLOR", "1")
    assert tm.Terminal(out).live is False                                      # StringIO is not a TTY anyway
    ui = tm.Terminal(out, color=False)
    ui.render(snap([row("a", "running", title="MS-FINDER", step="Unpacking")]))
    assert "\x1b" not in out.getvalue() and "[1/1] MS-FINDER" in out.getvalue()


def test_a_failure_shows_the_last_log_lines_and_the_first_line_of_the_error():
    out = Utf8()
    ui = tm.Terminal(out, color=True, width=100)
    ui.finish(snap([row("a", "error")], status="error", log=[f"line {i}" for i in range(20)], error="Docker is not running.\nmore"))
    text = ANSI.sub("", out.getvalue())
    assert "Setup failed." in text and "line 19" in text and "line 11" not in text and "Docker is not running." in text and "more" not in text


def test_an_unknown_amount_shows_a_walking_bar_not_a_percentage():
    out = Utf8()
    ui = tm.Terminal(out, color=True, width=100)
    ui.render(snap([row("a", "running", step="Compiling")], elapsed=3))
    assert "Compiling" in ANSI.sub("", out.getvalue()) and "%" not in ANSI.sub("", out.getvalue())
