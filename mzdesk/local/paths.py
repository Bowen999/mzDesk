"""Folders of the local version and the discovery of Java. Standard library only (the engines import this at the top).

Layout on the user's machine:
  ~/.mzdesk/          data/ models/ cache/ (existing), venv/ (the Python environment), bin/ (launcher, uv, micromamba), tools/ (engines)
  ~/.mzdesk/tools/    msfinder/ frigid/ mono/ git/ (installed on first use)
Everything is looked up when it is needed (never at import), so a test or a second user session with another HOME sees its own folders and an engine
that was installed a moment ago is found without a restart.
"""
import functools
import os
import shutil
import subprocess
import sys
from pathlib import Path

WINDOWS = sys.platform == "win32"


def home_dir():
    return Path.home() / ".mzdesk"


def bin_dir():
    return home_dir() / "bin"


def tools_dir():
    return home_dir() / "tools"


def exe(name):
    """Program file name on this platform (micromamba -> micromamba.exe on Windows)."""
    return name + ".exe" if WINDOWS else name


def tool_folder(env_var, legacy, name):
    """Folder of an engine: the environment variable if set, else the old default ``legacy`` (a path under ~/tools) when it already exists
    (the machine of a developer who installed the engine before), else ``~/.mzdesk/tools/<name>``."""
    value = os.environ.get(env_var)
    if value:
        return Path(value).expanduser()
    legacy = Path(legacy)
    return legacy if legacy.is_dir() else tools_dir() / name


def mono_prefix():
    """Conda environment with Mono, made by the installer when the computer has no usable Mono."""
    return tools_dir() / "mono"


def mono_in_prefix():
    return mono_prefix() / "bin" / "mono"


def find_executable(name, extra=()):
    """First existing, executable file for ``name``: ~/.mzdesk/bin, PATH, then the folders in ``extra``."""
    for cand in (bin_dir() / exe(name), shutil.which(name), *(Path(d) / exe(name) for d in extra)):
        if cand and Path(cand).is_file() and os.access(cand, os.X_OK):
            return str(cand)
    return None


# ---------------------------------------------------------------- Java
def _runs(java):
    try:
        done = subprocess.run([str(java), "-version"], capture_output=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired):
        return False
    return done.returncode == 0


@functools.lru_cache(maxsize=1)
def find_java():
    """Path of a Java that starts: ``java`` on PATH, else the runtime of the jdk4py package (part of the normal install), else None.

    macOS has a /usr/bin/java stub that only prints "Unable to locate a Java Runtime" when no JDK is installed, so a java on PATH is
    started once (`java -version`) before it is trusted. The answer is cached for the life of the process."""
    found = shutil.which("java")
    if found and _runs(found):
        return found
    try:
        from jdk4py import JAVA
    except ImportError:
        return None
    return str(JAVA) if Path(JAVA).is_file() else None


# ---------------------------------------------------------------- git
def working_git():
    """Path of a git that can be used: the one the installer made (~/.mzdesk/tools/git), else the one on PATH, else None.

    On macOS /usr/bin/git is a stub that opens the "install the command line developer tools" dialog when they are missing, so it only counts
    when `xcode-select -p` finds them."""
    own = tools_dir() / "git"
    for cand in (own / "bin" / exe("git"), own / "Library" / "cmd" / "git.exe", own / "Library" / "bin" / "git.exe", shutil.which("git")):
        if not cand or not Path(cand).is_file():
            continue
        if sys.platform == "darwin" and str(cand) == "/usr/bin/git":
            try:
                if subprocess.run(["xcode-select", "-p"], capture_output=True, timeout=30).returncode != 0:
                    continue
            except (OSError, subprocess.TimeoutExpired):
                continue
        return str(cand)
    return None
