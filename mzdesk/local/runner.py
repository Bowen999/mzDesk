"""Run the child processes of an install and stream their output; kill the whole process tree on cancel. Standard library only."""
import collections
import os
import re
import signal
import subprocess
import sys
import threading
import time

ANSI = re.compile(r"\x1b\[[0-9;?]*[A-Za-z]")
TAIL_LINES = 40                    # kept for the error message
KILL_AFTER_S = 5                   # SIGTERM first, SIGKILL when the tree is still there after this long


class InstallError(RuntimeError):
    """An install step failed; the message is shown to the user as it is."""


class Cancelled(Exception):
    """The user cancelled the install."""


class Ctx:
    """What an install step may do: write a log line, run a program (its output goes to the log), notice a cancel.

    ``log`` gets every output line of every child process (one call per line, progress lines that end in a carriage return included).
    One Ctx belongs to one install job; ``cancel()`` may be called from another thread."""

    def __init__(self, log=print):
        self._log = log
        self.tail = collections.deque(maxlen=TAIL_LINES)
        self.cancelled = threading.Event()
        self._proc = None
        self._lock = threading.Lock()

    def log(self, text):
        for line in str(text).splitlines():
            line = ANSI.sub("", line).rstrip()
            if line.strip():
                self.tail.append(line)
                self._log(line)

    def check(self):
        if self.cancelled.is_set():
            raise Cancelled()

    def cancel(self):
        self.cancelled.set()
        with self._lock:
            proc = self._proc
        if proc is not None:
            kill_tree(proc)

    def run(self, cmd, env=None, cwd=None, what=None, check=True):
        """Run ``cmd`` (a list; nothing goes through a shell), log its stdout and stderr line by line, return the exit code.

        check=True: a non-zero exit raises InstallError("<what> failed (exit code N)" + the last lines). Cancel raises Cancelled."""
        self.check()
        self.tail.clear()
        kwargs = {"creationflags": subprocess.CREATE_NEW_PROCESS_GROUP} if sys.platform == "win32" else {"start_new_session": True}
        try:
            proc = subprocess.Popen([str(c) for c in cmd], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                    env=env, cwd=cwd, **kwargs)
        except OSError as exc:
            raise InstallError(f"Could not start {cmd[0]}: {exc}") from None
        with self._lock:
            self._proc = proc
        if self.cancelled.is_set():              # cancel() ran between check() and the assignment
            kill_tree(proc)
        try:
            self._pump(proc)
            code = proc.wait()
        finally:
            with self._lock:
                self._proc = None
            if proc.poll() is None:
                kill_tree(proc)
        self.check()
        if check and code != 0:
            last = "\n".join(list(self.tail)[-12:])
            raise InstallError(f"{what or os.path.basename(str(cmd[0]))} failed (exit code {code})." + (f"\n{last}" if last else ""))
        return code

    def _pump(self, proc):
        """Read the child's output; both \\n and \\r end a line (pip and uv redraw one line with \\r)."""
        buf = b""
        for chunk in iter(lambda: proc.stdout.read1(4096), b""):
            buf += chunk
            parts = re.split(rb"[\r\n]", buf)
            buf = parts.pop()
            for part in parts:
                self.log(part.decode("utf-8", "replace"))
        if buf:
            self.log(buf.decode("utf-8", "replace"))


def kill_tree(proc):
    """Terminate ``proc`` and everything it started (the child runs in its own process group / session)."""
    if proc.poll() is not None:
        return
    if sys.platform == "win32":
        subprocess.run(["taskkill", "/PID", str(proc.pid), "/T", "/F"], capture_output=True)
        return
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except (ProcessLookupError, PermissionError):
        return
    deadline = time.monotonic() + KILL_AFTER_S
    while time.monotonic() < deadline:
        if proc.poll() is not None:
            break
        time.sleep(0.05)
    try:
        os.killpg(proc.pid, signal.SIGKILL)        # also takes grandchildren that ignored SIGTERM
    except (ProcessLookupError, PermissionError):
        pass
