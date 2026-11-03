"""The first-use downloads ("components") of the local version: id, size, installed(), install().

installed() looks at files only (the engines' own status functions); the one exception is cfmid, which asks Docker (one `docker image inspect`,
kept for 10 s). install(ctx) provisions what the component needs first (mzdesk.local.prereq), then runs the existing setup command as a child
process (`python -m mzdesk.x_setup`, or a throw-away `uv run --with ...` environment for the TensorFlow / PyTorch conversions). Folders come from
the engines (env var, else an existing ~/tools/..., else ~/.mzdesk/tools/...) and are read at call time, so a component installed a moment
ago works without a restart. Commands come from this file only; nothing a request says reaches a command line.
"""
import os
import shlex
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from mzdesk.local import paths, prereq
from mzdesk.local.runner import InstallError

ALL_PLATFORMS = ("darwin", "linux", "win32")
POSIX_ONLY = ("darwin", "linux")
POSIX_ONLY_MESSAGE = "macOS / Linux only for now"
LIBRARY_NAMES = ("MassBank", "GNPS")


def current_platform():
    return sys.platform


@dataclass(frozen=True)
class Component:
    id: str
    title: str
    size: str                       # download, as shown to the user
    ask: bool                       # the user confirms before it is downloaded (the MS/MS Annotation workflow overrides this, targets.py)
    platforms: tuple
    check: Callable[[], bool]       # installed? files only
    run: Callable                   # run(ctx): provision the prerequisites, run the setup command
    note: str = ""
    refresh: Callable = None        # in-process clean-up after a successful install
    first: tuple = field(default=())  # components that are installed before this one when they can be (ordering, see plan())

    def installed(self):
        try:
            return bool(self.check())
        except Exception:           # noqa: BLE001 - a broken status check must not break the status endpoint
            return False

    def supported(self, platform=None):
        return (platform or current_platform()) in self.platforms

    def unsupported_message(self):
        return POSIX_ONLY_MESSAGE if self.platforms == POSIX_ONLY else "not available on this platform"

    def install(self, ctx):
        """Run the install; raises InstallError / Cancelled. The component must be installed afterwards."""
        if not self.supported():
            raise InstallError(f"{self.title}: {self.unsupported_message()}.")
        self.run(ctx)
        ctx.check()
        if not self.installed():
            raise InstallError(f"{self.title}: the setup finished but the files are not complete. Run it again; the log above says what was missing.")
        if self.refresh:
            self.refresh()

    def info(self):
        ok = self.supported()
        return {"id": self.id, "title": self.title, "size": self.size, "installed": self.installed() if ok else False, "supported": ok,
                "ask": self.ask, "note": self.note if ok else self.unsupported_message()}


# ---------------------------------------------------------------- installed()
def _libraries_installed():
    from mzdesk.ms2 import library_index as lx
    return all((lx.index_dir(name) / "meta.json").is_file() for name in LIBRARY_NAMES)


def library_installed(name):
    """One library (MassBank | GNPS): is its index there?"""
    from mzdesk.ms2 import library_index as lx
    return (lx.index_dir(name) / "meta.json").is_file()


def _msfinder_installed():
    from mzdesk.ms2 import msfinder_engine as mf
    s = mf.status()
    return all(s[k] for k in ("mono", "exe", "resources", "ini", "thread_fix"))     # no thread fix: the setup runs again and applies it


def _frigid_installed():
    from mzdesk.ms2 import frigid_engine as fe
    return fe.installed()


def _cfmid_installed():
    from mzdesk.ms2 import cfmid_engine as cf
    return prereq.docker_state(cf.IMAGE) == "ready"


def _neims_installed():
    from mzdesk.ms2 import neims_setup
    return Path(os.environ.get("MZDESK_NEIMS_WEIGHTS") or neims_setup.default_weights()).is_file()


def _molgpka_installed():
    from mzdesk.platform_fit import pka_setup
    return Path(os.environ.get("MZDESK_MOLGPKA_WEIGHTS") or pka_setup.default_weights()).is_file()


def _ccs_installed():
    from mzdesk.orthogonal import ccs
    return ccs.installed()


def _repo_installed():
    from mzdesk.methods import report
    return report.installed()


# ---------------------------------------------------------------- install
def child_env(extra=None):
    """Environment of a setup command. mzDesk must be importable from any folder (a development checkout is not always installed)."""
    import mzdesk
    root = str(Path(mzdesk.__file__).resolve().parent.parent)
    env = prereq.clean_env(extra)
    env["PYTHONPATH"] = os.pathsep.join([root] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else []))
    return env


def run_setup(ctx, module, args=(), env=None, what=None):
    """`python -u -m <module> args` in a child process; the working folder is ~/.mzdesk so that no folder named mzdesk can shadow the package."""
    home = paths.home_dir()
    home.mkdir(parents=True, exist_ok=True)
    ctx.run([sys.executable, "-u", "-m", module, *args], env=env or child_env(), cwd=home, what=what or module)


def _install_libraries(ctx):
    run_setup(ctx, "mzdesk.ms2.library_setup", what="The library download")


def _install_msfinder(ctx):
    from mzdesk.ms2 import msfinder_engine as mf
    mono, own_bin = prereq.ensure_mono(ctx)
    env = child_env({"MZDESK_MONO": mono, "MZDESK_MSFINDER_DIR": mf.install_dir()})
    if own_bin:
        prereq.prepend_path(env, own_bin)
    run_setup(ctx, "mzdesk.ms2.msfinder_setup", env=env, what="The MS-FINDER setup")


def _install_frigid(ctx):
    from mzdesk.ms2 import frigid_engine as fe
    from mzdesk.ms2 import frigid_setup
    git = prereq.ensure_git(ctx)
    env = child_env({"MZDESK_FRIGID_DIR": fe.install_dir()})
    prereq.prepend_path(env, Path(git).parent)
    if not frigid_setup.find_conda():          # no conda on this computer: micromamba makes the Python 3.10 environment
        env.update(CONDA_EXE=prereq.ensure_micromamba(ctx), MAMBA_ROOT_PREFIX=str(prereq.cache_dir() / "mamba"), MAMBA_NO_BANNER="1")
    run_setup(ctx, "mzdesk.ms2.frigid_setup", env=env, what="The FRIGID setup")


def _install_cfmid(ctx):
    from mzdesk.ms2 import cfmid_engine as cf
    state = prereq.docker_state(cf.IMAGE, use_cache=False)
    if state == "no-docker":
        raise InstallError(f"Docker is not installed. CFM-ID runs in Docker; install Docker Desktop ({prereq.DOCKER_LINK}), start it and try again.")
    if state == "down":
        raise InstallError(cf.DOCKER_DOWN)
    if state == "no-image":
        ctx.run(shlex.split(cf.PULL_CMD), what="docker pull")
    prereq.forget_docker()


def _install_neims(ctx):
    from mzdesk.ms2 import neims_setup
    prereq.throwaway_python(ctx, ["tensorflow", "numpy"], neims_setup.__file__, "The NEIMS weight conversion", env={k: os.environ[k] for k in ("MZDESK_NEIMS_DIR",) if k in os.environ})


def _install_molgpka(ctx):
    from mzdesk.platform_fit import pka_setup
    prereq.throwaway_python(ctx, ["torch", "numpy"], pka_setup.__file__, "The MolGpKa weight conversion", env={k: os.environ[k] for k in ("MZDESK_MOLGPKA_DIR",) if k in os.environ})


def _install_ccs(ctx):
    run_setup(ctx, "mzdesk.orthogonal.ccs_setup", what="The CCS data setup")


def _install_repo(ctx):
    run_setup(ctx, "mzdesk.methods.report_setup", what="The RepoRT setup")


def _refresh_molgpka():
    """Pages that ran before the weights existed remember "no MolGpKa"; forget it (only modules that are loaded are touched)."""
    ci = sys.modules.get("mzdesk.methods.compound_index")
    if ci is not None:
        ci.reset_molgpka_check()
    rr = sys.modules.get("mzdesk.annotation.rerank")
    if rr is not None:
        rr._ion_check.cache_clear()
        rr._predicted_ions.cache_clear()


def _refresh_cfmid():
    prereq.forget_docker()


# ---------------------------------------------------------------- the registry (order = the order installs run in)
REGISTRY = {c.id: c for c in (
    Component("libraries", "Spectral libraries (MassBank, GNPS)", "2.3 GB", True, ALL_PLATFORMS, _libraries_installed, _install_libraries,
              note="MassBank 137 MB and GNPS 2.2 GB, then a local index of about 1.3 GB on disk."),
    Component("msfinder", "MS-FINDER structure prediction", "0.5 GB", True, POSIX_ONLY, _msfinder_installed, _install_msfinder,
              note="MS-FINDER 325 MB; adds Mono (130 MB, from conda-forge) when this computer has no usable Mono."),
    Component("frigid", "FRIGID de novo structures", "3 GB", True, POSIX_ONLY, _frigid_installed, _install_frigid,
              note="A PyTorch environment, the FRIGID code and two model files; 6 GB on disk. Adds conda and git when missing."),
    Component("cfmid", "CFM-ID 4 and MSRB (Docker image)", "0.2 GB", True, POSIX_ONLY, _cfmid_installed, _install_cfmid, refresh=_refresh_cfmid,
              note="A Docker image, about 1 GB on disk. Needs Docker Desktop running."),
    Component("neims", "NEIMS EI spectrum prediction", "1 GB", True, POSIX_ONLY, _neims_installed, _install_neims,
              note="The official weights are 466 MB; converting them needs TensorFlow (0.2-0.6 GB) once, in a throw-away environment."),
    Component("molgpka", "MolGpKa pKa model", "0.2 GB", False, ALL_PLATFORMS, _molgpka_installed, _install_molgpka, refresh=_refresh_molgpka,
              note="The weights are 12 MB; converting them needs PyTorch once, in a throw-away environment that is deleted afterwards."),
    Component("ccs", "CCS reference data and model", "8 MB", False, ALL_PLATFORMS, _ccs_installed, _install_ccs,
              note="Downloads the data files of one commit (4 MB) and trains the prediction model, a few minutes."),
    Component("repo", "RepoRT retention data", "0.2 GB", False, POSIX_ONLY, _repo_installed, _install_repo, first=("molgpka",),
              note="100 MB with git, 190 MB as a GitHub archive without it; builds a 40 MB index."),
)}


def get(component_id):
    return REGISTRY.get(component_id)


def status():
    return [c.info() for c in REGISTRY.values()]


def plan(ids):
    """Component ids to install, in install order. ``first`` components that are missing and supported join the list ahead of their dependant
    (RepoRT is classified with MolGpKa when the weights exist). Unknown id or a platform the component cannot have -> InstallError."""
    wanted = list(dict.fromkeys(ids))
    for cid in wanted:
        comp = REGISTRY.get(cid)
        if comp is None:
            raise InstallError(f"Unknown component {cid!r}.")
        if not comp.supported():
            raise InstallError(f"{comp.title}: {comp.unsupported_message()}.")
    for cid in list(wanted):
        for before in REGISTRY[cid].first:
            if before not in wanted and REGISTRY[before].supported() and not REGISTRY[before].installed():
                wanted.append(before)
    order = list(REGISTRY)
    return sorted(wanted, key=order.index)
