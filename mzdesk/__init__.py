"""mzDesk：小分子质谱一站式工作台。版本号只在这里定义（pyproject.toml 和 API 都从这里读）。"""
import os
import tempfile

# Matchms -> pynndescent -> numba compiles with cache=True. When numba cannot write a cache
# next to the installed package (read-only site-packages, a package-manager prefix, or a
# sandbox), importing matchms fails with "no locator available" instead of a useful error.
# Point numba at a writable, disposable cache before anything imports matchms. setdefault is
# deliberate: an explicit NUMBA_CACHE_DIR in the environment wins.
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(tempfile.gettempdir(), "mzdesk-numba"))
# Matplotlib builds a font cache on first import; the default ~/.matplotlib is not writable in
# some sandboxes and then every process rebuilds it in a fresh temp folder. Keep one cache.
os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "mzdesk-matplotlib"))

__version__ = "0.5.2"
