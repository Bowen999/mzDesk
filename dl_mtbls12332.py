#!/usr/bin/env python3
"""Download MTBLS12332 positive-mode raw files (.d.zip) with a live progress bar."""
import os, re, sys, time, subprocess, threading, urllib.request
from concurrent.futures import ThreadPoolExecutor

BASE = "https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS12332/FILES/RAW_FILES"
DEST = "/Volumes/G-TECH/MetabolomicsWorkbench/MTBLS12332"
DIRS = ["RP18_pos"]
PREFIX = "B9_Comp"   # only this condition
WORKERS = 24

def listing(d):
    h = urllib.request.urlopen(f"{BASE}/{d}/", timeout=60).read().decode("utf-8", "ignore")
    out = []
    for m in re.finditer(r'<a href="([^"]+\.zip)">.*?</td><td[^>]*>([^<]*)</td><td[^>]*>([^<]*)</td>', h):
        name, _date, size = m.groups()
        mm = re.match(r"([\d.]+)([KMGT])", size.strip().upper())
        b = int(float(mm.group(1)) * {"K": 1 << 10, "M": 1 << 20, "G": 1 << 30, "T": 1 << 40}[mm.group(2)]) if mm else 0
        out.append((name, b))
    return out

def human(n):
    for u in ["B", "KB", "MB", "GB", "TB"]:
        if n < 1024 or u == "TB":
            return f"{n:.1f} {u}" if u != "B" else f"{int(n)} B"
        n /= 1024

def fetch(url, path, expected):
    part = path + ".part"
    for attempt in range(3):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if have == expected:
            os.replace(part, path)
            return True
        r = subprocess.run(["curl", "-sS", "-L", "--fail", "--retry", "5", "--retry-delay", "5",
                            "-C", "-", "-o", part, url])
        if r.returncode == 0 or (os.path.exists(part) and os.path.getsize(part) == expected):
            os.replace(part, path)
            return True
    return False

def main():
    os.makedirs(DEST, exist_ok=True)
    jobs = []
    for d in DIRS:
        os.makedirs(os.path.join(DEST, d), exist_ok=True)
        for name, size in listing(d):
            if not name.startswith(PREFIX):
                continue
            jobs.append((f"{BASE}/{d}/{name}", os.path.join(DEST, d, name), size))
    total = sum(s for _, _, s in jobs)
    print(f"{len(jobs)} files, {human(total)} -> {DEST}", flush=True)

    done_bytes = [0] * len(jobs)
    done_flags = [False] * len(jobs)
    lock = threading.Lock()
    t0 = time.time()
    hist = []

    # sizes already present / partial on disk
    for i, (_, p, s) in enumerate(jobs):
        for cand in (p, p + ".part"):
            if os.path.exists(cand):
                done_bytes[i] = min(os.path.getsize(cand), s or os.path.getsize(cand))
                break
        if os.path.exists(p) and os.path.getsize(p) == s:
            done_flags[i] = True

    def worker(i):
        url, path, size = jobs[i]
        ok = fetch(url, path, size)
        with lock:
            done_flags[i] = ok
            if ok:
                done_bytes[i] = size or os.path.getsize(path)

    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = [ex.submit(worker, i) for i in range(len(jobs)) if not done_flags[i]]

        while any(not f.done() for f in futs):
            cur = 0
            for i, (_, p, s) in enumerate(jobs):
                if done_flags[i]:
                    cur += done_bytes[i]
                else:
                    part = p + ".part"
                    if os.path.exists(part):
                        done_bytes[i] = min(os.path.getsize(part), s or os.path.getsize(part))
                    cur += done_bytes[i]
            now = time.time()
            hist.append((now, cur))
            hist = [(t, b) for t, b in hist if now - t <= 20]
            speed = (cur - hist[0][1]) / max(now - hist[0][0], 1e-6) if len(hist) > 1 else 0
            pct = 100 * cur / total if total else 0
            filled = int(pct / 100 * 30)
            bar = "#" * filled + "-" * (30 - filled)
            eta = (total - cur) / speed if speed > 0 else 0
            eta_s = f"{int(eta // 3600)}h{int(eta % 3600 // 60):02d}m" if eta > 3600 else f"{int(eta // 60)}m{int(eta % 60):02d}s"
            n_done = sum(done_flags)
            sys.stdout.write(f"\r[{bar}] {pct:5.1f}%  {human(cur)}/{human(total)}  "
                             f"{human(speed)}/s  ETA {eta_s}  files {n_done}/{len(jobs)}")
            sys.stdout.flush()
            time.sleep(5)

    ok = sum(done_flags)
    dt = time.time() - t0
    sys.stdout.write("\n")
    print(f"done: {ok}/{len(jobs)} files, {human(sum(done_bytes))}, {dt/60:.1f} min", flush=True)
    bad = [os.path.basename(jobs[i][1]) for i in range(len(jobs)) if not done_flags[i]]
    if bad:
        print("FAILED:", *bad[:20], sep="\n  ", flush=True)
        return 1
    return 0

if __name__ == "__main__":
    sys.exit(main())
