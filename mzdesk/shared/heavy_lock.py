"""One heavy local job at a time: MS-FINDER (Mono, about 2.3 GB) and FRIGID (PyTorch, about 4-5 GB) share this lock.

Both engines run an external program that needs gigabytes of memory, so the server never runs two of them together, whichever
page started them. A second request waits (MS-FINDER: 20 s, FRIGID: 20 s) and is then turned away as busy. The lock is per server
process; the local version runs one uvicorn worker (`mzdesk` command). Standard library only.

    if not LOCK.acquire(timeout=20, holder="FRIGID"):
        ... busy (LOCK.holder says who is running)
    try: ...
    finally: LOCK.release()

`with LOCK:` still works (blocking, no holder name).
"""
import threading


class HeavyLock:
    """A threading.Lock that remembers the name of the job holding it, so a turned-away request can say what it is waiting for."""

    def __init__(self):
        self._lock = threading.Lock()
        self.holder = None

    def acquire(self, blocking=True, timeout=-1, holder=None):
        got = self._lock.acquire(blocking, timeout)
        if got:
            self.holder = holder
        return got

    def release(self):
        self.holder = None
        self._lock.release()

    def locked(self):
        return self._lock.locked()

    def __enter__(self):
        self.acquire()
        return self

    def __exit__(self, *exc):
        self.release()


LOCK = HeavyLock()
