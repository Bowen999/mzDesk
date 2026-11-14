"""How `mzdesk setup` looks in a terminal: one line per component, a live progress bar in mzDesk's blue, plain English. Standard library only.

    ui = Terminal(sys.stdout)
    ui.render(job.snapshot())          # call it a few times a second while the job runs; it prints what changed
    ui.finish(job.snapshot())          # the closing lines (done / failed)

On a terminal (a TTY, colours allowed) the current step is one line that is redrawn in place. When the output goes to a file or a pipe nothing is redrawn:
the log lines are printed as they come by the caller (cli.run_setup) and this class only prints the component lines. NO_COLOR turns the colours off.
"""
import os
import shutil
import sys

ACCENT = (123, 156, 240)          # --accent of mzDesk's dark theme: readable on a dark and on a light terminal
ALERT = (229, 72, 77)
CHECK, CROSS, DOT = "✓", "✗", "·"
FILL, EMPTY = "━", "─"


def enable_ansi():
    """Windows 10+: turn on escape sequences in the console of this process (a no-op elsewhere). True when they work."""
    if sys.platform != "win32":
        return True
    try:
        import ctypes
        k32 = ctypes.windll.kernel32
        handle = k32.GetStdHandle(-11)
        mode = ctypes.c_ulong()
        if not k32.GetConsoleMode(handle, ctypes.byref(mode)):
            return False
        return bool(k32.SetConsoleMode(handle, mode.value | 0x0004))          # ENABLE_VIRTUAL_TERMINAL_PROCESSING
    except Exception:  # noqa: BLE001 - an old console: plain text
        return False


def fmt_mb(mb):
    """412.0 -> '412 MB', 2300 -> '2.2 GB'."""
    return f"{mb / 1024:.1f} GB" if mb >= 1000 else f"{mb:.0f} MB"


def bar(fraction, width, fancy=True):
    fraction = max(0.0, min(1.0, fraction))
    full = round(width * fraction)
    return (FILL if fancy else "#") * full + (EMPTY if fancy else "-") * (width - full)


def step_text(row, limit=None):
    """What the running component is doing, in one short phrase: 'Downloading MSFINDER.v3.73.zip · 136 MB of 325 MB'.
    With ``limit`` the step is shortened (with an ellipsis) so that the whole phrase fits; the numbers are never cut."""
    text = row.get("step") or "Working"
    tail = f" {DOT} {fmt_mb(row['got_mb'])} of {fmt_mb(row['total_mb'])}" if row.get("progress") is not None and row.get("total_mb") else ""
    if limit is not None and len(text) + len(tail) > limit:
        text = text[: max(6, limit - len(tail) - 1)].rstrip() + "…"
    return text + tail


class Terminal:
    def __init__(self, stream=None, color=None, width=None):
        self.stream = stream or sys.stdout
        tty = bool(getattr(self.stream, "isatty", lambda: False)())
        allowed = tty and not os.environ.get("NO_COLOR") and os.environ.get("TERM") != "dumb"
        self.live = bool(color if color is not None else allowed) and (enable_ansi() if color is None else True)
        enc = (getattr(self.stream, "encoding", None) or "ascii").lower()
        try:
            (FILL + CHECK + CROSS + DOT).encode(enc)
            self.fancy = True
        except (UnicodeEncodeError, LookupError):
            self.fancy = False
        self.width = width or shutil.get_terminal_size((80, 20)).columns
        self.seen = {}                    # component id -> last status printed ("running" | "done" | "error")
        self.drawn = False                # a live line is on screen

    # ---- pieces
    def paint(self, text, rgb=None, bold=False, dim=False):
        if not self.live:
            return text
        codes = ([1] if bold else []) + ([2] if dim else []) + ([38, 2, *rgb] if rgb else [])
        return f"\x1b[{';'.join(map(str, codes))}m{text}\x1b[0m" if codes else text

    def _write(self, text):
        self.stream.write(text)
        self.stream.flush()

    def _clear(self):
        if self.live and self.drawn:
            self._write("\r\x1b[2K")
            self.drawn = False

    def say(self, text=""):
        self._clear()
        self._write(text + "\n")

    # ---- the job
    def title(self, rows):
        self.say(self.paint("mzDesk setup", ACCENT, bold=True) + self.paint(f"  {len(rows)} component{'s' if len(rows) != 1 else ''}", dim=True))

    def render(self, snap):
        """Print the component lines that changed and redraw the live line of the running one."""
        rows = snap["components"]
        total = len(rows)
        for i, row in enumerate(rows, 1):
            status = row["status"]
            if status == "waiting" or self.seen.get(row["id"]) == status:
                continue
            self.seen[row["id"]] = status
            head = self.paint(f"[{i}/{total}]", dim=True)
            if status == "running":
                self.say(f"{head} {self.paint(row['title'], bold=True)}")
            elif status == "done":
                self.say(f"{head} {self.paint(CHECK if self.fancy else 'ok', ACCENT, bold=True)} {row['title']}")
            elif status == "error":
                self.say(f"{head} {self.paint(CROSS if self.fancy else 'x', ALERT, bold=True)} {row['title']}")
        running = next((r for r in rows if r["status"] == "running"), None)
        if self.live and running:
            self._live_line(running, snap.get("elapsed_s", 0))

    def _live_line(self, row, elapsed):
        frac = row.get("progress")
        clock = f"{int(elapsed) // 60}:{int(elapsed) % 60:02d}"
        room = max(10, min(32, self.width // 3))
        label = step_text(row, self.width - room - len(clock) - 16)              # the whole line fits the window: it is redrawn in place
        if frac is None:                  # no number to show: a short bar that walks, so it is clear that something is happening
            span = 6
            pos = int(elapsed * 4) % (room + span)
            cells = "".join((FILL if self.fancy else "#") if pos - span <= k < pos else (EMPTY if self.fancy else "-") for k in range(room))
            gauge, pct = self.paint(cells, ACCENT), "    "
        else:
            filled = bar(frac, room, self.fancy)
            n = filled.count(FILL if self.fancy else "#")
            gauge = self.paint(filled[:n], ACCENT) + self.paint(filled[n:], dim=True)
            pct = f"{round(frac * 100):>3}%"
        line = f"      {gauge} {pct} {label}  {self.paint(clock, dim=True)}"
        self._clear()
        self._write(line)
        self.drawn = True

    def finish(self, snap):
        self._clear()
        if snap["status"] == "done":
            self.say(self.paint(f"{CHECK if self.fancy else 'ok'} All set.", ACCENT, bold=True) + self.paint(f"  {int(snap['elapsed_s'])} s", dim=True))
        elif snap["status"] == "cancelled":
            self.say(self.paint("Cancelled. Files that were already downloaded are kept.", dim=True))
        else:
            self.say(self.paint(f"{CROSS if self.fancy else 'x'} Setup failed.", ALERT, bold=True))
            for line in (snap.get("log") or [])[-8:]:
                self.say("  " + self.paint(line, dim=True))
            if snap.get("error"):
                self.say(self.paint(snap["error"].splitlines()[0], ALERT))
