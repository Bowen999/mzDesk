"""memguard.py - small, stdlib-only memory helpers for running heavy ML jobs on macOS (Apple silicon) and on Linux (CPU).

Linux: the numbers come from /proc (MemAvailable, SwapFree, VmRSS / VmHWM per process); there is no MPS, so "phys_footprint" is the RSS.

Copy this single file into another project (mzDesk etc.); it has no third-party imports.

Why it exists: on Apple silicon the MPS (GPU) backend allocates *unified memory* that is NOT shown in a process's RSS
(FRIGID on MPS: RSS ~1.5 GB but physical footprint 8-10 GB). So the cap must be on the physical footprint for MPS and
may be on RSS for CPU. Both are read from the kernel with proc_pid_rusage() (ctypes/libproc, ~microseconds, same
number the `footprint` CLI prints), which also exposes the kernel-tracked *lifetime peak* footprint, so peaks need no sampling.

Typical use
-----------
    import memguard

    # 1) decide the device at job start ("auto": MPS only if enough memory AND swap headroom, else CPU)
    d = memguard.choose_device("auto", mps_available=torch.backends.mps.is_available())
    device = d["chosen"]            # "mps" | "cpu";  d also holds the numbers + a human-readable reason for logs/JSON

    # 2) guard the job (in the job process itself, or from a parent supervising a child pid)
    guard = memguard.MemoryGuard(device, pid=None,              # None = this process (+ its children)
                                 on_exceed=my_handler)          # my_handler(reason: str) -> None; default kills the pid tree
    guard.start()
    ... heavy work ...
    summary = guard.stop()          # {"metric": "rss"|"phys_footprint", "cap_mb":..., "peak_mb":..., "exceeded": None|str, ...}

    # supervising a child from a parent process:
    p = subprocess.Popen([...]);  g = memguard.MemoryGuard("mps", pid=p.pid); g.start(); p.wait(); g.stop()

Caps (defaults): CPU -> process-tree RSS 6144 MB; MPS -> process-tree physical footprint 12288 MB.
Self-test:  python memguard.py --selftest
"""
from __future__ import annotations

import ctypes
import ctypes.util
import os
import re
import signal
import subprocess
import sys
import threading
import time

MB = 1024 * 1024
DEFAULT_CAPS_MB = {"cpu": 6144, "mps": 12288}
AUTO_MIN_AVAILABLE_MB = 12 * 1024  # MPS only if available memory >= 12 GB ...
AUTO_MIN_SWAP_FREE_MB = 1024       # ... AND swap free > 1 GB


# ----------------------------------------------------------------------------------------------------------------------
# system-level numbers
# ----------------------------------------------------------------------------------------------------------------------
def _sysctl(name: str) -> str:
    return subprocess.run(["sysctl", "-n", name], capture_output=True, text=True, timeout=5).stdout.strip()


def _to_mb(num: str, unit: str) -> float:
    return float(num) * {"K": 1 / 1024, "M": 1, "G": 1024, "T": 1024 * 1024}[unit.upper()]


IS_LINUX = sys.platform.startswith("linux")


def parse_meminfo(text: str) -> dict:
    """Linux: the numbers choose_device() uses, from the text of /proc/meminfo (kB values)."""
    kb = {key: int(rest.split()[0]) for key, rest in (line.split(":", 1) for line in text.splitlines() if ":" in line) if rest.split()}
    total, avail = kb["MemTotal"] / 1024, kb["MemAvailable"] / 1024
    strict = (kb.get("MemFree", 0) + kb.get("Inactive(file)", 0)) / 1024
    swap_total, swap_free = kb.get("SwapTotal", 0) / 1024, kb.get("SwapFree", 0) / 1024
    return dict(total_mb=round(total), memory_free_pct=round(100 * avail / total), available_mb_kernel=round(avail), available_mb_strict=round(strict),
                swap_total_mb=round(swap_total), swap_used_mb=round(swap_total - swap_free), swap_free_mb=round(swap_free))


def parse_proc_status(text: str) -> dict | None:
    """Linux: {rss_mb, phys_footprint_mb, lifetime_max_footprint_mb} from the text of /proc/<pid>/status (VmRSS, VmHWM = the peak RSS)."""
    kb = {key: int(rest.split()[0]) for key, rest in (line.split(":", 1) for line in text.splitlines() if ":" in line) if rest.split() and rest.split()[0].isdigit()}
    if "VmRSS" not in kb:
        return None
    rss = kb["VmRSS"] / 1024
    return dict(rss_mb=rss, phys_footprint_mb=rss, lifetime_max_footprint_mb=kb.get("VmHWM", kb["VmRSS"]) / 1024)


def system_memory() -> dict:
    """Numbers used by choose_device(). All sizes in MB.

    available_mb_kernel : total * kern.memorystatus_level/100 - the same "System-wide memory free percentage" that
                          `memory_pressure` prints (counts reclaimable + compressible memory; best predictor of whether a
                          big allocation will be granted without trouble)
    available_mb_strict : free + inactive + speculative + purgeable pages from vm_stat (no compression credit; much lower
                          on a busy Mac)
    """
    if IS_LINUX:
        with open("/proc/meminfo") as f:
            return parse_meminfo(f.read())
    total_mb = int(_sysctl("hw.memsize")) / MB
    level = int(_sysctl("kern.memorystatus_level"))
    out = subprocess.run(["vm_stat"], capture_output=True, text=True, timeout=5).stdout
    page = int(re.search(r"page size of (\d+) bytes", out).group(1))
    pages = {k: int(v) for k, v in re.findall(r"Pages ([a-z ]+):\s+(\d+)\.", out)}
    strict = sum(pages.get(k, 0) for k in ("free", "inactive", "speculative", "purgeable")) * page / MB
    sw = re.search(r"total = ([\d.]+)([KMGT]).*?used = ([\d.]+)([KMGT]).*?free = ([\d.]+)([KMGT])", _sysctl("vm.swapusage"))
    return dict(total_mb=round(total_mb), memory_free_pct=level, available_mb_kernel=round(total_mb * level / 100),
                available_mb_strict=round(strict),
                swap_total_mb=round(_to_mb(sw.group(1), sw.group(2))), swap_used_mb=round(_to_mb(sw.group(3), sw.group(4))),
                swap_free_mb=round(_to_mb(sw.group(5), sw.group(6))))


def choose_device(requested: str = "auto", mps_available: bool = True, min_available_mb: int = AUTO_MIN_AVAILABLE_MB,
                  min_swap_free_mb: int = AUTO_MIN_SWAP_FREE_MB, avail_metric: str = "kernel") -> dict:
    """requested: "auto" | "mps" | "cpu". "auto" -> MPS only if MPS exists AND available >= min_available_mb AND
    swap_free > min_swap_free_mb, else CPU. avail_metric: "kernel" (default, see system_memory) or "strict".
    Explicit "mps"/"cpu" are honoured (mps falls back to cpu only if MPS does not exist). Never raises on probe failure:
    if the numbers cannot be read, "auto" picks CPU."""
    res = dict(requested=requested, mps_available=bool(mps_available), avail_metric=avail_metric,
               thresholds=dict(min_available_mb=min_available_mb, min_swap_free_mb=min_swap_free_mb))
    try:
        res.update(system_memory())
    except Exception as e:  # noqa: BLE001
        res["probe_error"] = f"{type(e).__name__}: {e}"
    avail = res.get(f"available_mb_{avail_metric}")
    if requested == "cpu":
        res.update(chosen="cpu", reason="cpu requested")
    elif requested == "mps":
        ok = bool(mps_available)
        res.update(chosen="mps" if ok else "cpu", reason="mps requested" if ok else "mps requested but unavailable -> cpu")
    elif requested == "auto":
        if not mps_available:
            res.update(chosen="cpu", reason="auto: MPS not available")
        elif avail is None or res.get("swap_free_mb") is None:
            res.update(chosen="cpu", reason="auto: memory probe failed -> cpu")
        else:
            why = []
            if avail < min_available_mb:
                why.append(f"available {avail} MB < {min_available_mb} MB")
            if res["swap_free_mb"] <= min_swap_free_mb:
                why.append(f"swap free {res['swap_free_mb']} MB <= {min_swap_free_mb} MB")
            res.update(chosen="cpu" if why else "mps",
                       reason=("auto: cpu because " + " and ".join(why)) if why else
                       f"auto: mps (available {avail} MB >= {min_available_mb}, swap free {res['swap_free_mb']} MB > {min_swap_free_mb})")
    else:
        raise ValueError("requested must be 'auto', 'mps' or 'cpu'")
    return res


# ----------------------------------------------------------------------------------------------------------------------
# per-process numbers (kernel proc_pid_rusage, flavor 4)
# ----------------------------------------------------------------------------------------------------------------------
class _RUsageV4(ctypes.Structure):  # struct rusage_info_v4 from <sys/resource.h>
    _fields_ = [("ri_uuid", ctypes.c_uint8 * 16)] + [(n, ctypes.c_uint64) for n in (
        "ri_user_time", "ri_system_time", "ri_pkg_idle_wkups", "ri_interrupt_wkups", "ri_pageins", "ri_wired_size",
        "ri_resident_size", "ri_phys_footprint", "ri_proc_start_abstime", "ri_proc_exit_abstime", "ri_child_user_time",
        "ri_child_system_time", "ri_child_pkg_idle_wkups", "ri_child_interrupt_wkups", "ri_child_pageins",
        "ri_child_elapsed_abstime", "ri_diskio_bytesread", "ri_diskio_byteswritten", "ri_cpu_time_qos_default",
        "ri_cpu_time_qos_maintenance", "ri_cpu_time_qos_background", "ri_cpu_time_qos_utility", "ri_cpu_time_qos_legacy",
        "ri_cpu_time_qos_user_initiated", "ri_cpu_time_qos_user_interactive", "ri_billed_system_time",
        "ri_serviced_system_time", "ri_logical_writes", "ri_lifetime_max_phys_footprint", "ri_instructions", "ri_cycles",
        "ri_billed_energy", "ri_serviced_energy", "ri_interval_max_phys_footprint", "ri_runnable_time")]


_libproc = None


def proc_memory(pid: int) -> dict | None:
    """{rss_mb, phys_footprint_mb, lifetime_max_footprint_mb} for one pid, or None if it is gone / not readable."""
    global _libproc
    if IS_LINUX:
        try:
            with open(f"/proc/{int(pid)}/status") as f:
                return parse_proc_status(f.read())
        except OSError:
            return None
    if _libproc is None:
        _libproc = ctypes.CDLL(ctypes.util.find_library("proc"), use_errno=True)
    r = _RUsageV4()
    if _libproc.proc_pid_rusage(int(pid), 4, ctypes.byref(r)) != 0:
        return None
    return dict(rss_mb=r.ri_resident_size / MB, phys_footprint_mb=r.ri_phys_footprint / MB,
                lifetime_max_footprint_mb=r.ri_lifetime_max_phys_footprint / MB)


def tree_pids(root: int) -> list[int]:
    """root + all descendants (via ps). Falls back to [root] when ps is unavailable or blocked."""
    root = int(root)
    try:
        out = subprocess.run(["ps", "-axo", "pid=,ppid="], capture_output=True, text=True, timeout=5).stdout
    except (OSError, subprocess.TimeoutExpired):
        return [root]
    kids: dict[int, list[int]] = {}
    for line in out.splitlines():
        p = line.split()
        if len(p) == 2:
            try:
                kids.setdefault(int(p[1]), []).append(int(p[0]))
            except ValueError:
                continue
    seen, stack = [], [root]
    while stack:
        x = stack.pop()
        if x not in seen:
            seen.append(x)
            stack.extend(kids.get(x, []))
    return seen or [root]


def tree_memory(root: int, include_children: bool = True) -> dict:
    """Sum of rss / phys_footprint over the process tree. lifetime_max is the ROOT process's kernel-tracked peak."""
    pids = tree_pids(root) if include_children else [int(root)]
    rss = fp = 0.0
    root_mem = None
    for p in pids:
        m = proc_memory(p)
        if m:
            rss += m["rss_mb"]
            fp += m["phys_footprint_mb"]
            if p == root:
                root_mem = m
    return dict(rss_mb=rss, phys_footprint_mb=fp, n_procs=len(pids),
                root_lifetime_max_footprint_mb=(root_mem or {}).get("lifetime_max_footprint_mb"))


# ----------------------------------------------------------------------------------------------------------------------
# the guard
# ----------------------------------------------------------------------------------------------------------------------
def _kill_tree(pid: int) -> None:
    for p in sorted(tree_pids(pid), reverse=True):  # children first
        try:
            os.kill(p, signal.SIGKILL)
        except ProcessLookupError:
            pass


class MemoryGuard:
    """Background thread that polls memory of `pid` (default: this process) + children and fires `on_exceed(reason)` once
    when the cap for the device is crossed. cpu -> RSS cap; mps -> physical-footprint cap (RSS misses unified memory).
    Default on_exceed: SIGKILL the whole process tree (if pid is this process that includes ourselves)."""

    def __init__(self, device: str, pid: int | None = None, cap_mb: int | None = None, poll_s: float = 0.5,
                 include_children: bool = True, on_exceed=None):
        if device not in DEFAULT_CAPS_MB:
            raise ValueError("device must be 'cpu' or 'mps'")
        self.device, self.pid = device, int(pid or os.getpid())
        self.metric = "phys_footprint_mb" if device == "mps" else "rss_mb"
        self.cap_mb = int(cap_mb or DEFAULT_CAPS_MB[device])
        self.poll_s, self.include_children = poll_s, include_children
        self.on_exceed = on_exceed or (lambda reason: _kill_tree(self.pid))
        self.peak_mb = 0.0
        self.peak_rss_mb = 0.0
        self.peak_footprint_mb = 0.0
        self.exceeded: str | None = None
        self.n_polls = 0
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None

    def _poll(self) -> None:
        m = tree_memory(self.pid, self.include_children)
        self.n_polls += 1
        self.peak_rss_mb = max(self.peak_rss_mb, m["rss_mb"])
        self.peak_footprint_mb = max(self.peak_footprint_mb, m["phys_footprint_mb"])
        v = m[self.metric]
        self.peak_mb = max(self.peak_mb, v)
        if v > self.cap_mb and self.exceeded is None:
            self.exceeded = (f"{self.device}: process-tree {'physical footprint' if self.device == 'mps' else 'RSS'} "
                             f"{v:.0f} MB > cap {self.cap_mb} MB")
            try:
                self.on_exceed(self.exceeded)
            except Exception:  # noqa: BLE001  a broken handler must not leave the job running past its cap
                _kill_tree(self.pid)

    def _run(self) -> None:
        while not self._stop.is_set():
            try:
                self._poll()
            except Exception:  # noqa: BLE001  (a failed probe must never kill the job)
                pass
            self._stop.wait(self.poll_s)

    def start(self) -> "MemoryGuard":
        self._thread = threading.Thread(target=self._run, name="memguard", daemon=True)
        self._thread.start()
        return self

    def stop(self) -> dict:
        self._stop.set()
        if self._thread and self._thread is not threading.current_thread():
            self._thread.join(timeout=2)
        try:
            self._poll()
        except Exception:  # noqa: BLE001
            pass
        return self.snapshot()

    def snapshot(self) -> dict:
        """Current summary; never joins the polling thread, so it is safe to call from inside `on_exceed`."""
        life = proc_memory(self.pid)
        return dict(device=self.device, metric="phys_footprint" if self.device == "mps" else "rss", cap_mb=self.cap_mb,
                    peak_mb=round(self.peak_mb, 1), peak_rss_mb=round(self.peak_rss_mb, 1),
                    peak_footprint_mb=round(self.peak_footprint_mb, 1),
                    kernel_lifetime_peak_footprint_mb=round(life["lifetime_max_footprint_mb"], 1) if life else None,
                    exceeded=self.exceeded, n_polls=self.n_polls, poll_s=self.poll_s, include_children=self.include_children)

    def __enter__(self):
        return self.start()

    def __exit__(self, *a):
        self.stop()


# ----------------------------------------------------------------------------------------------------------------------
def _selftest() -> int:
    ok = True
    print("system_memory:", system_memory())
    for req in ("auto", "cpu", "mps"):
        d = choose_device(req, mps_available=True)
        print(f"choose_device({req!r}) -> {d['chosen']}: {d['reason']}")
    # synthetic decisions with thresholds forced both ways
    assert choose_device("auto", True, min_available_mb=1, min_swap_free_mb=-1)["chosen"] == "mps"
    assert choose_device("auto", True, min_available_mb=10**9)["chosen"] == "cpu"
    assert choose_device("auto", False, min_available_mb=1, min_swap_free_mb=-1)["chosen"] == "cpu"
    # footprint reading vs `footprint` CLI, and the guard firing on a tiny cap
    blob = bytearray(300 * MB)
    for i in range(0, len(blob), 4096):
        blob[i] = 1
    mine = proc_memory(os.getpid())
    cli = subprocess.run(["footprint", str(os.getpid())], capture_output=True, text=True).stdout
    m = re.search(r"Footprint: ([\d.]+) ([KMG]B)", cli)
    cli_mb = _to_mb(m.group(1), m.group(2)[0]) if m else float("nan")
    print(f"phys_footprint via proc_pid_rusage {mine['phys_footprint_mb']:.0f} MB vs footprint CLI {cli_mb:.0f} MB")
    ok &= abs(mine["phys_footprint_mb"] - cli_mb) < 0.05 * cli_mb + 10
    fired = []
    with MemoryGuard("cpu", cap_mb=100, poll_s=0.1, on_exceed=fired.append) as g:
        time.sleep(0.5)
    ok &= bool(fired) and g.exceeded is not None
    with MemoryGuard("mps", cap_mb=10**6, poll_s=0.1, on_exceed=fired.append) as g2:
        time.sleep(0.3)
    ok &= g2.exceeded is None and g2.peak_mb >= 300
    print("guard fired on tiny cap:", g.exceeded, "| quiet under big cap, peak", round(g2.peak_mb), "MB")
    del blob
    # parent supervising a child pid, with a deliberately broken handler: the guard must still kill the child
    child = subprocess.Popen([sys.executable, "-I", "-c",
                              "import time\nb=bytearray(400*1024*1024)\nfor i in range(0,len(b),4096): b[i]=1\ntime.sleep(30)"])
    def broken(reason):
        raise RuntimeError("handler bug")
    g3 = MemoryGuard("cpu", pid=child.pid, cap_mb=150, poll_s=0.2, on_exceed=broken).start()
    try:
        rc = child.wait(timeout=15)
    except subprocess.TimeoutExpired:
        child.kill()
        rc = None
    g3.stop()
    print("child supervised by guard with broken handler -> exit code", rc, "(negative = killed by signal)")
    ok &= rc is not None and rc < 0 and g3.exceeded is not None
    print("SELFTEST", "OK" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_selftest() if "--selftest" in sys.argv else print(__doc__) or 0)
