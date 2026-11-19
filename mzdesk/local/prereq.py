"""What an install provisions before it runs a setup command: micromamba, Mono (with its C# compiler), git, a throw-away Python environment.

Every function takes the install's Ctx (mzdesk.local.runner): its log lines reach the dialog, a cancel stops the child process.
Nothing here edits a shell profile, installs system-wide or asks for administrator rights; everything lands in ~/.mzdesk.
Windows: MS-FINDER needs no Mono there (the .NET Framework runs it); micromamba, git from conda-forge and the throw-away environment are written
for it and are exercised by the Windows CI job (.github/workflows/windows.yml).
"""
import io
import os
import platform
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path

from mzdesk.local import paths
from mzdesk.local.runner import Ctx, InstallError  # noqa: F401  (InstallError is re-exported for the callers)

MICROMAMBA_URL = "https://micro.mamba.pm/api/micromamba/{platform}/latest"       # redirects to conda-forge's micromamba package (tar.bz2, about 7 MB)
MICROMAMBA_PLATFORMS = {("darwin", "arm64"): "osx-arm64", ("darwin", "x86_64"): "osx-64", ("linux", "x86_64"): "linux-64",
                        ("linux", "aarch64"): "linux-aarch64", ("win32", "amd64"): "win-64"}
MONO_SPEC = "mono=6.12.0.199"            # conda-forge; the version MS-FINDER's console front end was compiled and tested with
PYTHON_FOR_UV = "3.12"                   # TensorFlow and PyTorch both have wheels for it
DOCKER_LINK = "https://docs.docker.com/get-docker/"
ROSETTA_COMMAND = "softwareupdate --install-rosetta --agree-to-license"
INSTALL_COMMAND = "curl -fsSL https://mzdesk.com/install.sh | sh"
INSTALL_COMMAND_WINDOWS = "irm https://mzdesk.com/install.ps1 | iex"


# ---------------------------------------------------------------- environment of the child processes
def clean_env(extra=None):
    """The environment for a child process: ours, unbuffered output, UTF-8; no activated Python environment leaks in."""
    env = {k: v for k, v in os.environ.items() if k not in ("VIRTUAL_ENV", "PYTHONHOME", "CONDA_PREFIX", "CONDA_DEFAULT_ENV")}
    env.update(PYTHONUNBUFFERED="1", PYTHONIOENCODING="utf-8")
    env.update({k: str(v) for k, v in (extra or {}).items() if v is not None})
    return env


def prepend_path(env, *folders):
    """Put ``folders`` at the front of PATH in ``env`` (a dict, changed in place)."""
    env["PATH"] = os.pathsep.join([str(f) for f in folders if f] + [env.get("PATH", "")])
    return env


def cache_dir():
    d = paths.home_dir() / "cache"
    d.mkdir(parents=True, exist_ok=True)
    return d


# ---------------------------------------------------------------- micromamba
def micromamba_platform():
    found = MICROMAMBA_PLATFORMS.get((sys.platform, platform.machine().lower()))
    if not found:
        raise InstallError(f"There is no micromamba build for {sys.platform} / {platform.machine()}.")
    return found


def find_micromamba():
    return paths.find_executable("micromamba")


def ensure_micromamba(ctx):
    """Path of micromamba: the one in ~/.mzdesk/bin or on PATH, else the official build downloaded into ~/.mzdesk/bin (about 7 MB)."""
    found = find_micromamba()
    if found:
        return found
    plat = micromamba_platform()
    url = MICROMAMBA_URL.format(platform=plat)
    ctx.log(f"Downloading micromamba ({plat}) from {url}")
    ctx.check()
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "mzdesk-setup"}), timeout=60) as resp:
            data = resp.read()
    except (urllib.error.URLError, OSError) as exc:
        raise InstallError(f"Could not download micromamba ({exc}). Check the network connection and try again.") from None
    member = "Library/bin/micromamba.exe" if plat == "win-64" else "bin/micromamba"
    target = paths.bin_dir() / paths.exe("micromamba")
    try:
        with tarfile.open(fileobj=io.BytesIO(data), mode="r:bz2") as tf:
            src = tf.extractfile(member)
            if src is None:
                raise InstallError(f"The micromamba package from {url} has no {member}.")
            paths.bin_dir().mkdir(parents=True, exist_ok=True)
            part = target.with_name(target.name + ".part")
            part.write_bytes(src.read())
    except (tarfile.TarError, KeyError, OSError) as exc:
        raise InstallError(f"The micromamba download from {url} is damaged ({exc}).") from None
    part.chmod(part.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    part.replace(target)
    ctx.log(f"micromamba installed: {target}")
    return str(target)


def mamba_env(extra=None):
    """Environment for micromamba: its package cache goes to ~/.mzdesk/cache/mamba (cleaned after each environment is made)."""
    return clean_env({"MAMBA_ROOT_PREFIX": cache_dir() / "mamba", "MAMBA_NO_BANNER": "1", **(extra or {})})


def conda_create(ctx, prefix, specs, subdir=None, what="conda environment"):
    """micromamba create -p ``prefix`` from conda-forge only; ``subdir`` ('osx-64') makes it for another platform. Frees the package cache afterwards."""
    mamba = ensure_micromamba(ctx)
    env = mamba_env({"CONDA_SUBDIR": subdir} if subdir else None)
    Path(prefix).parent.mkdir(parents=True, exist_ok=True)
    ctx.log(f"Creating {what} ({', '.join(specs)}) with micromamba")
    try:
        ctx.run([mamba, "create", "-y", "-p", prefix, "--override-channels", "-c", "conda-forge", *specs], env=env, what=f"Creating the {what}")
    except InstallError:
        shutil.rmtree(prefix, ignore_errors=True)               # a half-made environment would look like an installed one
        raise
    ctx.run([mamba, "clean", "-a", "-y"], env=env, check=False)
    return mamba


# ---------------------------------------------------------------- Mono (MS-FINDER)
def rosetta_ready():
    """True when x86_64 programs run on this Apple Silicon Mac (Rosetta 2 is installed)."""
    try:
        return subprocess.run(["arch", "-x86_64", "/usr/bin/true"], capture_output=True, timeout=30).returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False


def mono_works(mono):
    """(True, "") when this mono comes with a C# compiler that starts and the facade assemblies MS-FINDER's front end is built against."""
    from mzdesk.ms2 import msfinder_setup as setup

    csc, facades = setup.find_csc(mono), setup.find_facades(mono)
    if not csc or not facades:
        return False, "no C# compiler (csc) or no facade assemblies next to it"
    try:
        done = subprocess.run([csc, "-version"], capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return False, f"csc does not start ({exc})"
    return (True, "") if done.returncode == 0 else (False, "csc does not start")


def ensure_mono(ctx, system=True):
    """-> (mono path, folder to put first on PATH or None). A Mono that works is used: the one this installer made earlier, then (``system``) a Mono on
    this computer whose csc starts. Else Mono 6.12 is installed from conda-forge into ~/.mzdesk/tools/mono with micromamba.

    conda-forge has no osx-arm64 build: on Apple Silicon the environment is made for osx-64 and Mono runs under Rosetta 2."""
    from mzdesk.ms2 import msfinder_engine as mf

    if os.name == "nt":                         # MS-FINDER is a Windows program: the .NET Framework that comes with Windows runs it, no Mono
        if mf.find_mono():
            ctx.log("Runtime: .NET Framework (built into Windows)")
            return mf.NATIVE, None
        raise InstallError("MS-FINDER needs the .NET Framework 4.7.2 or newer. Install it from https://dotnet.microsoft.com/download/dotnet-framework, "
                           "then try again.")
    mine = paths.mono_in_prefix()
    candidates = [str(mine)] if mine.is_file() else []
    if system:
        found = mf.find_mono()
        if found and found not in candidates:
            candidates.append(found)
    for mono in candidates:
        ok, why = mono_works(mono)
        if ok:
            ctx.log(f"Mono: {mono}")
            return mono, _own_bin(mono)
        ctx.log(f"The Mono at {mono} cannot be used ({why}).")
    subdir = None
    if sys.platform == "darwin" and platform.machine().lower() == "arm64":
        if not rosetta_ready():
            raise InstallError("Mono for Apple Silicon runs under Rosetta 2, which is not installed on this Mac. Install it once in Terminal with "
                               f"`{ROSETTA_COMMAND}`, then try again.")
        subdir = "osx-64"
    prefix = paths.mono_prefix()
    if prefix.exists():
        shutil.rmtree(prefix, ignore_errors=True)
    conda_create(ctx, prefix, [MONO_SPEC], subdir=subdir, what="Mono environment")
    mono = str(paths.mono_in_prefix())
    ok, why = mono_works(mono) if Path(mono).is_file() else (False, "mono is missing in the new environment")
    if not ok:
        shutil.rmtree(prefix, ignore_errors=True)
        raise InstallError(f"Mono was installed but does not work here ({why}).")
    ctx.log(f"Mono installed: {mono}")
    return mono, _own_bin(mono)


def _own_bin(mono):
    """The bin folder of the Mono this installer made (its csc must come first on PATH: msfinder_setup looks there), None for any other Mono."""
    return Path(mono).parent if paths.mono_prefix() in Path(mono).parents else None


# ---------------------------------------------------------------- git (FRIGID)
def ensure_git(ctx):
    """Path of a git that works: the system one, else git from conda-forge in ~/.mzdesk/tools/git."""
    git = paths.working_git()
    if git:
        return git
    prefix = paths.tools_dir() / "git"
    if prefix.exists():
        shutil.rmtree(prefix, ignore_errors=True)
    conda_create(ctx, prefix, ["git"], what="git environment")
    git = paths.working_git()
    if not git:
        raise InstallError("git was installed from conda-forge but does not start.")
    return git


# ---------------------------------------------------------------- uv: the throw-away environment
def find_uv():
    return paths.find_executable("uv", extra=(Path.home() / ".local" / "bin", Path.home() / ".cargo" / "bin"))


def throwaway_python(ctx, requirements, script, what, env=None, args=()):
    """Run the standard-library-only setup script ``script`` with `uv run --with <requirements>` in a throw-away environment.

    The packages (TensorFlow, PyTorch) are only needed once, to convert the downloaded weights, so they never touch mzDesk's own environment: uv
    builds an ephemeral one with its own cache folder, which is deleted afterwards. `python -I` keeps the script's folder and every PYTHON*
    variable out of sys.path. The script must import nothing from mzdesk."""
    uv = find_uv()
    if not uv:
        raise InstallError("uv was not found. mzDesk's installer puts it in ~/.mzdesk/bin; run the installer again: "
                           + (INSTALL_COMMAND_WINDOWS if os.name == "nt" else INSTALL_COMMAND))
    work = Path(tempfile.mkdtemp(prefix="uv-", dir=cache_dir()))
    cmd = [uv, "run", "--no-project", "--python", PYTHON_FOR_UV]
    if sys.platform.startswith("linux"):
        # the default PyTorch wheel on Linux is the CUDA build (2 GB+); the CPU build is all the conversion needs
        cmd += ["--extra-index-url", "https://download.pytorch.org/whl/cpu", "--index-strategy", "unsafe-best-match"]
    for req in requirements:
        cmd += ["--with", req]
    cmd += ["python", "-I", str(script), *args]
    child = clean_env({"UV_CACHE_DIR": work, "UV_NO_CONFIG": "1", **(env or {})})
    child.pop("PYTHONPATH", None)
    try:
        ctx.log(f"Starting a throw-away Python environment with {', '.join(requirements)} (deleted afterwards)")
        ctx.run(cmd, env=child, what=what)
    finally:
        shutil.rmtree(work, ignore_errors=True)


# ---------------------------------------------------------------- Docker (CFM-ID)
_docker_cache = (0.0, None)
DOCKER_TTL_S = 10


def docker_state(image, use_cache=True):
    """'no-docker' (not installed) | 'down' (installed, the daemon does not answer) | 'no-image' (running, the image is not pulled) | 'ready'.

    One `docker image inspect` call (a few ms when Docker runs, a timeout of 8 s when it hangs); the answer is kept for 10 s."""
    global _docker_cache
    stamp, state = _docker_cache
    if use_cache and state and time.monotonic() - stamp < DOCKER_TTL_S:
        return state
    state = _docker_probe(image)
    _docker_cache = (time.monotonic(), state)
    return state


def _docker_probe(image):
    if shutil.which("docker") is None:
        return "no-docker"
    try:
        done = subprocess.run(["docker", "image", "inspect", image], capture_output=True, text=True, timeout=8)
    except subprocess.TimeoutExpired:
        return "down"
    except OSError:
        return "no-docker"
    if done.returncode == 0:
        return "ready"
    return "no-image" if "No such image" in (done.stderr or "") else "down"


def docker_os_type():
    """'linux' or 'windows': the kind of container the running Docker makes ('' when it does not answer). Docker Desktop on Windows can be
    switched to Windows containers (and GitHub's Windows machines only have those); CFM-ID's image is a Linux one."""
    try:
        done = subprocess.run(["docker", "info", "--format", "{{.OSType}}"], capture_output=True, text=True, timeout=15)
    except (OSError, subprocess.TimeoutExpired):
        return ""
    return done.stdout.strip().lower() if done.returncode == 0 else ""


WINDOWS_CONTAINERS = ("Docker is set to Windows containers; CFM-ID needs Linux containers. Right-click the Docker icon in the taskbar, "
                      "choose \"Switch to Linux containers...\", and try again.")


def forget_docker():
    global _docker_cache
    _docker_cache = (0.0, None)
