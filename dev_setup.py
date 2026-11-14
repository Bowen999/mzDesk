"""Get a git checkout of mzDesk ready for development and the whole test suite, in one command:

    python dev_setup.py                      everything
    python dev_setup.py --skip frigid,libraries   leave out large components (the tests that need them then fail or skip)
    python dev_setup.py --no-pip --no-web    only the MetFrag jar and the components
    python dev_setup.py --no-pip --no-web --skip all   only the MetFrag jar

Steps (each one is skipped when its result is already there, so running it again is cheap):
  1. Python packages: requirements.txt plus pytest and httpx.
  2. MetFrag jar: not in git (68 MB); downloaded from the MetFrag release on GitHub and checked by SHA-256.
  3. Web interface: `npm ci` + `npm run build` in web/ -> mzdesk/static/ (needs Node.js; skipped with a note when there is none).
  4. macOS: the OpenMP runtime LightGBM needs (python -m mzdesk.shared.openmp --install), then the components: the same installs as
     `mzdesk setup --yes <id>` (the ones the local version downloads on first use). CFM-ID is skipped when Docker is not running, and a
     component this computer cannot run (FRIGID on an Intel Mac) is left out with its reason.
  5. Summary of what is in place.
The progress shown is a rough estimate (weighted by download size), not a measurement. Exit code 1 when a step failed; later steps still run.
"""
import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
JAR = ROOT / "mzdesk" / "lib" / "MetFragCommandLine-2.5.2.jar"
JAR_URL = "https://github.com/ipb-halle/MetFragRelaunched/releases/download/v2.5.2/MetFragCommandLine-2.5.2.jar"
JAR_SHA256 = "aefe802cd964bc88142636d67306f53ef3700038dfdea4efca2b5a45c8bbc796"       # the v2.5.2 release file, unmodified
STATIC_INDEX = ROOT / "mzdesk" / "static" / "index.html"

# rough weight of each piece in the progress line (about its download size in 100 MB; small steps get 1)
WEIGHTS = {"pip": 5, "jar": 1, "web": 3, "libraries": 23, "msfinder": 5, "frigid": 30, "cfmid": 2, "neims": 10, "molgpka": 2,
           "ccs": 1, "repo": 2, "redu": 17}

T0 = time.monotonic()
done_weight = 0
total_weight = 1


def say(step, text):
    pct = min(99, int(100 * done_weight / total_weight))
    mins, secs = divmod(int(time.monotonic() - T0), 60)
    print(f"\n[{step}] about {pct}% done, {mins}:{secs:02d} elapsed: {text}", flush=True)


def finish(key):
    global done_weight
    done_weight += WEIGHTS.get(key, 1)


def run(cmd, cwd=None):
    print("  $ " + " ".join(str(c) for c in cmd), flush=True)
    return subprocess.run([str(c) for c in cmd], cwd=cwd).returncode == 0


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def step_pip():
    return run([sys.executable, "-m", "pip", "install", "-r", ROOT / "requirements.txt", "pytest", "httpx"])


def step_jar():
    if JAR.is_file() and sha256(JAR) == JAR_SHA256:
        print("  already there")
        return True
    part = JAR.with_suffix(".jar.part")
    print(f"  downloading {JAR_URL}")
    with urllib.request.urlopen(JAR_URL, timeout=60) as resp, open(part, "wb") as out:
        size, got, shown = int(resp.headers.get("Content-Length") or 0), 0, 0
        while block := resp.read(1 << 20):
            out.write(block)
            got += len(block)
            if size and got * 4 // size > shown:              # a line at about 25 %, 50 %, 75 %, 100 %
                shown = got * 4 // size
                print(f"  {got >> 20} / {size >> 20} MB", flush=True)
    digest = sha256(part)
    if digest != JAR_SHA256:
        part.unlink()
        print(f"  the download has SHA-256 {digest}, expected {JAR_SHA256}; not used")
        return False
    os.replace(part, JAR)
    return True


def step_web():
    npm = shutil.which("npm")
    if not npm:
        print("  Node.js (npm) is not installed: the web interface is not built. Install Node.js 20 or newer and run this again.")
        return False
    web = ROOT / "web"
    if not (web / "node_modules").is_dir() and not run([npm, "ci", "--no-audit", "--no-fund"], cwd=web):
        return False
    return run([npm, "run", "build"], cwd=web) and STATIC_INDEX.is_file()


def docker_running():
    docker = shutil.which("docker")
    return bool(docker) and subprocess.run([docker, "info"], capture_output=True).returncode == 0


def main(argv=None):
    global total_weight, done_weight
    ap = argparse.ArgumentParser(description="Get a git checkout of mzDesk ready for development and the tests.")
    ap.add_argument("--skip", default="", help="comma-separated component ids to leave out, e.g. frigid,libraries; `all` leaves out every component")
    ap.add_argument("--no-pip", action="store_true", help="do not install the Python packages")
    ap.add_argument("--no-web", action="store_true", help="do not build the web interface")
    args = ap.parse_args(argv)

    sys.path.insert(0, str(ROOT))
    skip = {s.strip() for s in args.skip.split(",") if s.strip()}
    plan = [] if args.no_pip else ["pip"]
    plan += ["jar"] + ([] if args.no_web else ["web"])
    skip_all = "all" in skip
    total_weight = sum(WEIGHTS.get(k, 1) for k in plan) + sum(v for k, v in WEIGHTS.items()
                                                               if k not in ("pip", "jar", "web") and k not in skip and not skip_all)
    results = {}

    n_steps = 5 - args.no_pip - args.no_web
    i = 0
    for key, title, fn in (("pip", "Python packages", step_pip), ("jar", "MetFrag jar", step_jar), ("web", "Web interface", step_web)):
        if key not in plan:
            continue
        i += 1
        say(f"step {i}/{n_steps}", title)
        try:
            results[title] = fn()
        except Exception as exc:                        # a network error in one step does not stop the others
            print(f"  failed: {exc}")
            results[title] = False
        finish(key)

    if sys.platform == "darwin":                        # BUDDY's LightGBM model needs libomp, which most Macs lack (mzdesk/shared/openmp.py)
        say(f"step {i}/{n_steps}", "OpenMP runtime for LightGBM (macOS)")
        results["OpenMP runtime"] = run([sys.executable, "-m", "mzdesk.shared.openmp", "--install"], cwd=ROOT)

    i += 1
    from mzdesk.local import components
    ids = [] if skip_all else [cid for cid in components.REGISTRY if cid not in skip]
    if "cfmid" in ids and not docker_running():
        ids.remove("cfmid")
        print("\nCFM-ID is left out: Docker is not running (start Docker Desktop and run this again to add it).")
        results["cfmid (no Docker)"] = None
    elif "cfmid" in ids:
        from mzdesk.local import prereq
        if prereq.docker_os_type() == "windows":
            ids.remove("cfmid")
            print("\nCFM-ID is left out: " + prereq.WINDOWS_CONTAINERS)
            results["cfmid (Windows containers)"] = None
    for k, cid in enumerate(ids, 1):
        comp = components.REGISTRY[cid]
        say(f"step {i}/{n_steps}, component {k}/{len(ids)}", f"{comp.title} ({comp.size} download)")
        if not comp.supported():
            print(f"  left out: {comp.unsupported_message()}")
            results[f"{cid} (not on this computer)"] = None
        elif comp.installed():
            print("  already there")
            results[cid] = True
        else:
            results[cid] = run([sys.executable, "-m", "mzdesk", "setup", "--yes", cid], cwd=ROOT)
        finish(cid)

    i += 1
    done_weight = total_weight
    say(f"step {i}/{n_steps}", "Summary")
    for name, ok in results.items():
        print(f"  {'ok     ' if ok else ('skipped' if ok is None else 'FAILED ')}  {name}")
    for cid in sorted(skip):
        print(f"  skipped  {cid} (--skip)")
    mins, secs = divmod(int(time.monotonic() - T0), 60)
    failed = [n for n, ok in results.items() if ok is False]
    print(f"\n{'Done' if not failed else 'Done with failures: ' + ', '.join(failed)} in {mins}:{secs:02d}. Run the tests with: python -m pytest tests")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
