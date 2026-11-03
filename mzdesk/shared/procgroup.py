"""Start a child in its own process group, measure its memory, kill the whole group: one place for macOS, Linux and Windows. Standard library only.

    proc = subprocess.Popen(cmd, **procgroup.popen_kwargs())
    procgroup.group_rss_mb(proc.pid)     # resident memory of the child and everything it started, in MB (0.0 when it cannot be read)
    procgroup.kill_group(proc)           # SIGKILL to the group (POSIX) / taskkill /T /F (Windows), then wait up to 10 s for the child
"""
import os
import signal
import subprocess
import sys

WINDOWS = sys.platform == "win32"


def popen_kwargs():
    """Popen arguments that put the child in a process group (POSIX: its own session; Windows: new group, and no console window)."""
    if WINDOWS:
        return {"creationflags": subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.CREATE_NO_WINDOW}
    return {"start_new_session": True}


def group_rss_mb(pid):
    """Resident memory (MB) of ``pid`` and its descendants. POSIX: `ps` for the process group, else the kernel's number for the leader."""
    if WINDOWS:
        try:
            from mzdesk.shared.memguard import tree_memory
            return float(tree_memory(int(pid))["rss_mb"])
        except Exception:  # noqa: BLE001 - a missing probe must not kill the run
            return 0.0
    out = ""
    try:
        out = subprocess.run(["ps", "-axo", "pgid=,rss="], capture_output=True, text=True, timeout=5).stdout
    except (OSError, subprocess.TimeoutExpired):
        out = ""
    total = 0
    for line in out.splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0] == str(pid) and parts[1].isdigit():
            total += int(parts[1])
    if total:
        return total / 1024.0
    try:
        from mzdesk.shared.memguard import proc_memory
        return float((proc_memory(int(pid)) or {}).get("rss_mb", 0.0))
    except Exception:  # noqa: BLE001
        return 0.0


def kill_group(proc, wait_s=10):
    """Kill the child's whole process group and wait for the child."""
    if WINDOWS:
        subprocess.run(["taskkill", "/PID", str(proc.pid), "/T", "/F"], capture_output=True)
    else:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            pass
    try:
        proc.wait(timeout=wait_s)
    except subprocess.TimeoutExpired:
        pass
