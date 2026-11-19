"""macOS: the OpenMP runtime (libomp) that LightGBM needs, so that BUDDY's formula model loads on every Mac.

msbuddy's formula model is a LightGBM model. PyPI's lightgbm wheels for macOS link ``@rpath/libomp.dylib`` and look for it only in the
Homebrew and MacPorts folders, which most Macs do not have: without it every BUDDY run fails ("Library not loaded: @rpath/libomp.dylib").
mzDesk installs conda-forge's llvm-openmp into ~/.mzdesk/tools/libomp (install.sh, ``python -m mzdesk.shared.openmp --install``,
dev_setup.py). ensure_lightgbm() loads it before LightGBM only when LightGBM cannot load by itself: LightGBM then finds it, already loaded,
by its install name. A Mac that has libomp already (Homebrew, or PyTorch loaded first) is left alone, so two OpenMP runtimes never meet.
Other systems: nothing to do (the Linux and Windows wheels bring their OpenMP runtime).

    python -m mzdesk.shared.openmp [--check | --install]
"""
import ctypes
import sys

from mzdesk.local import paths

SPEC = "llvm-openmp"                     # conda-forge; libomp.dylib with the install name @rpath/libomp.dylib
MISSING = ("BUDDY needs the OpenMP runtime (libomp) on macOS, which this Mac does not have. Run `python -m mzdesk.shared.openmp --install` "
           "(about 1 MB from conda-forge; `mzdesk update` does it too), or `brew install libomp`, then restart mzDesk.")
_state = {"ok": False}


def prefix():
    return paths.tools_dir() / "libomp"


def dylib():
    return prefix() / "lib" / "libomp.dylib"


def _import_lightgbm():
    import lightgbm  # noqa: F401  (loads lib_lightgbm.dylib, which needs libomp)


def ensure_lightgbm():
    """Make sure LightGBM loads; raises RuntimeError(MISSING) on a Mac without libomp. Cheap after the first call."""
    if _state["ok"]:
        return
    try:
        _import_lightgbm()
    except OSError:
        if sys.platform != "darwin" or not dylib().is_file():
            raise RuntimeError(MISSING) from None
        ctypes.CDLL(str(dylib()), mode=ctypes.RTLD_GLOBAL)
        try:
            _import_lightgbm()
        except OSError as exc:
            raise RuntimeError(f"{MISSING} (mzDesk's {dylib()} did not help: {exc})") from None
    _state["ok"] = True


def install(ctx=None):
    """Install llvm-openmp from conda-forge into prefix() with micromamba (macOS only; elsewhere nothing to do). Returns the dylib path or None."""
    if sys.platform != "darwin":
        return None
    if dylib().is_file():
        return dylib()
    from mzdesk.local import prereq
    from mzdesk.local.runner import Ctx

    prereq.conda_create(ctx or Ctx(), prefix(), [SPEC], what="OpenMP runtime for LightGBM")
    if not dylib().is_file():
        raise prereq.InstallError(f"{SPEC} was installed but {dylib()} is not there.")
    return dylib()


def main(argv=None):
    args = sys.argv[1:] if argv is None else argv
    if sys.platform != "darwin":
        print("Nothing to do: only macOS needs a separate OpenMP runtime for LightGBM.")
        return 0
    if "--install" in args:
        print(f"OpenMP runtime: {install()}")
    try:
        ensure_lightgbm()
    except RuntimeError as exc:
        print(exc)
        return 1
    print("LightGBM loads." + (f" mzDesk's OpenMP runtime: {dylib()}" if dylib().is_file() else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
