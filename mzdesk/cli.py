"""The mzdesk command.

Usage:
  mzdesk                      start mzDesk at http://localhost:8731 (if that port is busy, the next free one up to 8750; the address is printed)
  mzdesk --open               also open the browser (the launcher the installer writes passes this)
  mzdesk --port 8599          use this port only (no search; for development: python -m mzdesk --port 8599)
  mzdesk --host 0.0.0.0       listen on all addresses (Docker)
  mzdesk --reload             development mode: restart when the backend code changes
  mzdesk setup [--list] [--check] [--yes] [id ... | all]
                              install, in the terminal, the components that are otherwise downloaded on first use
                              (libraries msfinder frigid cfmid neims molgpka ccs repo); asks y/N for each one, --yes does not ask,
                              --list / --check only show the status
  mzdesk update [--check]     --check only shows the newest version; otherwise downloads $MZDESK_BASE_URL/install.sh (Windows: install.ps1) and runs it
  mzdesk version

Environment: MZDESK_BASE_URL (default https://mzdesk.com: the server of the installer scripts and the version file), MZDESK_UPDATE_URL (the version file,
default <BASE>/version.json), MZDESK_NO_UPDATE_CHECK=1 (do not look for a new version). The frontend is in web/; `npm run build` writes it to
mzdesk/static/, which the backend serves.
"""
import argparse
import itertools
import os
import socket
import subprocess
import sys
import tempfile
import threading
import time
import urllib.error
import urllib.request
import webbrowser
from pathlib import Path

from mzdesk import __version__

DEFAULT_PORT, LAST_PORT = 8731, 8750
SUBCOMMANDS = ("setup", "update", "version")


# ---------------------------------------------------------------- run
def port_is_free(host, port):
    """True when a server could listen on (host, port) now."""
    try:
        family = socket.AF_INET6 if ":" in host else socket.AF_INET
        with socket.socket(family, socket.SOCK_STREAM) as sock:
            if sys.platform != "win32":                # like uvicorn: a port in TIME_WAIT is free; on Windows SO_REUSEADDR would share a live port
                sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            sock.bind((host, port))
        return True
    except OSError:
        return False


def find_port(host, first=DEFAULT_PORT, last=LAST_PORT):
    """First free port from ``first`` to ``last`` (inclusive); None when all are taken."""
    return next((p for p in range(first, last + 1) if port_is_free(host, p)), None)


def browser_url(host, port):
    return f"http://{'localhost' if host in ('0.0.0.0', '::', '') else host}:{port}"


def _open_when_up(host, port, wait_s=60.0):
    """Open the browser as soon as the server accepts connections (called in a thread; gives up quietly after ``wait_s``)."""
    target = "127.0.0.1" if host in ("0.0.0.0", "") else ("::1" if host == "::" else host)
    deadline = time.monotonic() + wait_s
    while time.monotonic() < deadline:
        try:
            with socket.create_connection((target, port), timeout=0.5):
                break
        except OSError:
            time.sleep(0.2)
    else:
        return
    try:
        webbrowser.open(browser_url(host, port))
    except Exception:  # noqa: BLE001 - no browser on this machine: the URL was printed
        pass


def run_server(args):
    from mzdesk import update

    port = args.port
    if port is None:
        port = find_port(args.host)
        if port is None:
            print(f"Ports {DEFAULT_PORT}-{LAST_PORT} are all in use. Close another mzdesk, or choose a port: mzdesk --port 9731", file=sys.stderr)
            return 1
    line = update.terminal_line()
    if line:
        print(line, flush=True)
    update.check_in_background()
    url = browser_url(args.host, port)
    print(f"mzDesk {__version__} is running at {url}  (Ctrl+C to stop)", flush=True)
    if args.open:
        threading.Thread(target=_open_when_up, args=(args.host, port), daemon=True).start()
    import uvicorn
    uvicorn.run("mzdesk.api.app:app", host=args.host, port=port, reload=args.reload)
    return 0


# ---------------------------------------------------------------- setup
def _status_line(comp):
    if not comp.supported():
        return f"not available: {comp.unsupported_message()}"
    return "installed" if comp.installed() else "not installed"


def run_setup(args, ask=input, out=print):
    """`mzdesk setup`: 0 when everything asked for is installed, 1 otherwise. ``ask`` and ``out`` are the terminal (replaced in the tests)."""
    from mzdesk.local import components, jobs, terminal
    from mzdesk.local.runner import InstallError

    registry = components.REGISTRY
    names = list(args.ids)
    unknown = [n for n in names if n != "all" and n not in registry]
    if unknown:
        out(f"Unknown component: {', '.join(unknown)}. Components: {', '.join(registry)}.")
        return 2
    chosen = list(registry) if (not names or "all" in names) else list(dict.fromkeys(names))
    if args.list or args.check:
        out("Components (size is the download):")
        for cid in chosen:
            comp = registry[cid]
            out(f"  {cid:<10} {comp.size:>7}  {_status_line(comp):<30} {comp.title}")
        missing = [c for c in chosen if registry[c].supported() and not registry[c].installed()]
        return 1 if (args.check and missing) else 0
    if args.yes and not names:
        out("Name the components to install, or `all`: mzdesk setup --yes all")
        return 2

    todo = []
    for cid in chosen:
        comp = registry[cid]
        if not comp.supported():
            out(f"{cid}: not available on this computer ({comp.unsupported_message()}).")
        elif comp.installed():
            out(f"{cid}: already installed.")
        elif args.yes or ask(f"Install {comp.title} ({comp.size} download)? [y/N] ").strip().lower() in ("y", "yes"):
            todo.append(cid)
        else:
            out(f"{cid}: skipped.")
    if not todo:
        return 0
    ui = terminal.Terminal(sys.stdout) if out is print else None
    live = bool(ui and ui.live)                 # a terminal: the themed live view; a file or a pipe: the plain log lines, as before
    try:
        job = jobs.JobStore().start(todo, on_line=None if live else out)
    except InstallError as exc:
        out(str(exc))
        return 1
    try:
        if live:
            ui.title(job.snapshot()["components"])
            while job.thread.is_alive():
                ui.render(job.snapshot())
                job.thread.join(0.25)
            ui.render(job.snapshot())
        else:
            job.thread.join()
    except KeyboardInterrupt:
        if live:
            ui.say()
        out("Cancelling ...")
        job.ctx.cancel()
        job.thread.join()
        return 130
    snap = job.snapshot()
    if live:
        ui.finish(snap)
        return 0 if snap["status"] == "done" else 1
    for row in snap["components"]:
        out(f"{row['id']}: {row['status']}")
    if snap["status"] != "done":
        out(snap["error"] or f"The install ended with status {snap['status']}.")
        return 1
    return 0


# ---------------------------------------------------------------- update
def run_update(args, out=print):
    from mzdesk import update

    got = update.refresh(timeout=8)
    latest = got.get("latest") if update.parse_version(got.get("latest")) else None
    if latest is None:
        out(f"mzDesk {__version__}. Could not read the newest version from {update.update_url()}.")
        return 1
    if args.check:
        out(f"Current {__version__}, latest {latest}" + (f". Update: {update.COMMAND}" if update.is_newer(latest) else " (up to date)."))
        return 0
    if not update.is_newer(latest):
        out(f"mzDesk {__version__} is up to date.")
        return 0
    return run_installer(update.base_url(), out)


def run_installer(base, out=print):
    """Download the installer script of the server and run it with MZDESK_UPDATING=1 (it installs the new version over this one).

    POSIX: run it and wait. Windows: the running mzdesk holds .pyd files open, so the script starts in a new console window, detached, and this
    process ends right away (not verified on Windows)."""
    windows = sys.platform == "win32"
    name = "install.ps1" if windows else "install.sh"
    url = f"{base}/{name}"
    out(f"Updating mzDesk with {url}")
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": f"mzdesk/{__version__}"}), timeout=30) as resp:
            script = resp.read(2_000_000)
    except (urllib.error.URLError, OSError) as exc:
        out(f"Could not download {url} ({exc}). Check the network connection and try again.")
        return 1
    env = {**os.environ, "MZDESK_UPDATING": "1"}
    if windows:
        path = Path(tempfile.mkdtemp(prefix="mzdesk-update-")) / name          # stays: the detached script runs after this process has ended
        path.write_bytes(script)
        flags = getattr(subprocess, "CREATE_NEW_CONSOLE", 0x10) | getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0x200)
        try:
            subprocess.Popen(["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File", str(path)], env=env, creationflags=flags, close_fds=True)
        except OSError as exc:
            out(f"Could not run the installer ({exc}).")
            return 1
        out("Updating in a new window; mzDesk will be closed.")
        return 0
    with tempfile.TemporaryDirectory(prefix="mzdesk-update-") as tmp:
        path = Path(tmp) / name
        path.write_bytes(script)
        try:
            return subprocess.run(["sh", str(path)], env=env).returncode
        except OSError as exc:
            out(f"Could not run the installer ({exc}).")
            return 1


# ---------------------------------------------------------------- entry point
def build_parsers():
    run = argparse.ArgumentParser(prog="mzdesk", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    run.add_argument("--port", type=int, default=None, help=f"port (default {DEFAULT_PORT}, or the next free one up to {LAST_PORT})")
    run.add_argument("--host", default="127.0.0.1")
    run.add_argument("--open", action="store_true", help="open the browser when the server is up")
    run.add_argument("--reload", action="store_true", help="development: restart on code changes")
    run.add_argument("-V", "--version", action="store_true", help="print the version and exit")
    setup = argparse.ArgumentParser(prog="mzdesk setup", description="Install the components that pages download on first use.")
    setup.add_argument("ids", nargs="*", metavar="id", help="component ids, or `all` (default: all, asking for each)")
    setup.add_argument("--list", action="store_true", help="show the components, their sizes and whether they are installed")
    setup.add_argument("--check", action="store_true", help="like --list; exit code 1 when one is missing")
    setup.add_argument("--yes", "-y", action="store_true", help="do not ask")
    update = argparse.ArgumentParser(prog="mzdesk update", description="Install the newest mzDesk (runs the installer again).")
    update.add_argument("--check", action="store_true", help="only print the newest version")
    return run, setup, update


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    run, setup, update = build_parsers()
    lead = len(argv) - len(list(itertools.dropwhile(lambda a: a == "--open", argv)))
    if lead and argv[lead:lead + 1] and argv[lead] in SUBCOMMANDS:        # the launcher always adds --open: `mzdesk version` arrives as `--open version`
        argv = argv[lead:]
    if argv and argv[0] in SUBCOMMANDS:
        command, rest = argv[0], argv[1:]
        if command == "version":
            print(f"mzDesk {__version__}")
            return 0
        if command == "setup":
            return run_setup(setup.parse_args(rest))
        return run_update(update.parse_args(rest))
    args = run.parse_args(argv)
    if args.version:
        print(f"mzDesk {__version__}")
        return 0
    return run_server(args)


if __name__ == "__main__":
    sys.exit(main())
