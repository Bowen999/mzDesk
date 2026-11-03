"""Versions of everything a result depends on: mzDesk, Python, the Python packages, the engines, the models, the data and the online services.

GET /api/versions (api/app.py) returns environment(); the browser puts it into report.json of the Session report (web/src/report/build.ts).
The keys of "components" are the short names of the pages' `based_on` (spec.py), plus a few more key packages. Each value is
    {"version": "..." | None, "kind": "python-package" | "engine" | "model" | "data" | "online" | "mzdesk", "detail": "..." | None}
A version of None means "not installed here" or "could not be found" (detail says which). An online service has no version (kind "online"):
the report records when it was asked instead.

Names, versions, dates and counts only: never a folder, a user name, a host name or an environment variable (every text goes through strip_paths).
Nothing here raises: a resolver that fails gives {"version": None, "detail": "unknown"}. The answer is computed lazily, under a lock, and kept for
TTL_S seconds, so a request does not wait for the probes. The only slow things are three FIXED commands (docker image inspect, mono --version,
java -version), each with a short time limit; no user input ever reaches a command.
"""
import base64
import copy
import gzip
import hashlib
import importlib.metadata as metadata
import importlib.util
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import threading
import time
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path

from mzdesk import __version__
from mzdesk.api.security import strip_paths

TTL_S = 600                 # a model downloaded on first use, or Docker started later, shows up within ten minutes
DOCKER_TIMEOUT_S = 4
RUNTIME_TIMEOUT_S = 3       # mono --version, java -version
KINDS = ("python-package", "engine", "model", "data", "online", "mzdesk")
ONLINE = "online service, not versioned"


# ---------------------------------------------------------------- building blocks
def _entry(version, kind, detail=None):
    return {"version": version, "kind": kind, "detail": detail}


def _absent(kind, why="not installed"):
    return _entry(None, kind, why)


def _day(path):
    """Date a file was last written, as YYYY-MM-DD (the server's clock, like the dates in the library indexes)."""
    return time.strftime("%Y-%m-%d", time.localtime(Path(path).stat().st_mtime))


def _release(text):
    """'2024.09.6' -> (2024, 9, 6): the numbers of a version, to compare two spellings of it."""
    found = re.match(r"\d+(?:\.\d+)*", str(text))
    return tuple(int(n) for n in found.group(0).split(".")) if found else ()


def _module_file(module):
    """The file an import of `module` runs (its __init__.py), found without importing it; None when there is none."""
    try:
        loaded = sys.modules.get(module)
        origin = getattr(loaded, "__file__", None) or getattr(importlib.util.find_spec(module), "origin", None)
    except Exception:                   # not importable, or a broken package
        return None
    return Path(origin) if origin and origin not in ("built-in", "frozen") else None


def _record_matches(dist, origin):
    """True when the RECORD of this metadata folder lists `origin` with the hash the file has now, i.e. this is the release that wrote it."""
    try:
        for entry in dist.files or []:
            if entry.hash and entry.hash.mode == "sha256" and Path(dist.locate_file(entry)).resolve() == origin.resolve():
                digest = base64.urlsafe_b64encode(hashlib.sha256(origin.read_bytes()).digest()).rstrip(b"=").decode()
                return digest == entry.hash.value
    except Exception:
        pass
    return False


def _dist_version(*names, module=None):
    """Version of a Python package, or None when it is not installed. Read from the installed metadata (importlib.metadata); nothing is imported.

    The package's own __version__ is not used: some packages never update it (msbuddy 0.3.12 says "0.3.9"). An upgrade can leave the old metadata
    folder beside the new one (this happens in conda environments), so with several versions listed the one whose RECORD matches the file `module`
    imports wins, and without a match the highest version."""
    found = []
    for name in names:
        try:
            found += [d for d in metadata.distributions(name=name) if d.version]
        except Exception:               # a broken metadata folder
            continue
    if len({str(d.version) for d in found}) > 1 and module:
        origin = _module_file(module)
        written = [d for d in found if origin and _record_matches(d, origin)]
        found = written or found
    return str(max(found, key=lambda d: _release(d.version)).version) if found else None


def python_package(*dists, module=None, detail=None):
    """Entry of a Python package; `module` is its import name (see _dist_version)."""
    version = _dist_version(*dists, module=module)
    return _entry(version, "python-package", detail) if version else _absent("python-package")


def _online(detail=ONLINE):
    return _entry(None, "online", detail)


def _run(cmd, timeout):
    """(exit code, output) of a FIXED command, or None when it cannot start or does not finish in time. cmd never holds user input."""
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, stdin=subprocess.DEVNULL)
    except (OSError, subprocess.SubprocessError, ValueError):
        return None
    return proc.returncode, (proc.stdout or "") + (proc.stderr or "")


@lru_cache(maxsize=1)
def _docker_image():
    """('ok', short image id) when the CFM-ID image is local, else ('missing', why). One `docker image inspect`, never a pull."""
    from mzdesk.ms2 import cfmid_engine
    if shutil.which("docker") is None:
        return "missing", "Docker not found"
    out = _run(["docker", "image", "inspect", "--format", "{{.Id}}", cfmid_engine.IMAGE], DOCKER_TIMEOUT_S)
    if out is None:
        return "missing", "Docker did not answer"
    code, text = out
    if code == 0:
        found = re.search(r"sha256:([0-9a-f]{12})", text)
        return "ok", found.group(1) if found else ""
    return "missing", "Docker image not pulled" if "No such image" in text else "Docker is not running"


@lru_cache(maxsize=1)
def _mono_version():
    """Version of the Mono that runs MS-FINDER (MZDESK_MONO, then PATH), or None."""
    from mzdesk.ms2 import msfinder_engine
    mono = msfinder_engine.find_mono()
    out = _run([mono, "--version"], RUNTIME_TIMEOUT_S) if mono and mono != msfinder_engine.NATIVE else None
    found = re.search(r"compiler version ([\d.]+)", out[1]) if out and out[0] == 0 else None
    return found.group(1) if found else None


@lru_cache(maxsize=1)
def _java_version():
    """Version of the Java that runs MetFrag (java on PATH, else the one inside the jdk4py package), or None."""
    out = _run(["java", "-version"], RUNTIME_TIMEOUT_S) if shutil.which("java") else None
    if not (out and out[0] == 0):                        # no usable java on PATH (macOS has a stub that only complains): the one in jdk4py
        from mzdesk.local.paths import find_java
        bundled = find_java()
        out = _run([bundled, "-version"], RUNTIME_TIMEOUT_S) if bundled else None
    found = re.search(r'version "([^"]+)"', out[1]) if out and out[0] == 0 else None
    return found.group(1) if found else None


def _git_head(clone):
    """Commit a git clone is on, read from .git/HEAD (a detached hash or a branch name), or None."""
    try:
        head = (Path(clone) / ".git" / "HEAD").read_text(encoding="utf-8").strip()
        if head.startswith("ref:"):
            head = (Path(clone) / ".git" / head[4:].strip()).read_text(encoding="utf-8").strip()
    except OSError:
        return None
    return head if re.fullmatch(r"[0-9a-f]{40}", head) else None


# ---------------------------------------------------------------- mzDesk and the Python packages
def _mzdesk():
    return _entry(__version__, "mzdesk")


def _massql():
    # the page runs a MassQL subset that is part of mzDesk (gnps/compound_db_check.py); the massql package is not used
    from mzdesk.gnps import compound_db_check
    return _entry(compound_db_check.TOOL_VERSION, "engine", "MassQL subset built into mzDesk; the massql package is not used")


def _zenodo_record(url):
    found = re.search(r"/records/(\d+)/", str(url))
    return found.group(1) if found else "?"


def _ms2deepscore():
    from mzdesk.ms2 import ms2_similarity as sim
    version = _dist_version("ms2deepscore", module="ms2deepscore")
    if version is None:
        return _absent("python-package")
    model = Path(sim.CACHE_DIR) / "ms2deepscore_model.pt"
    record = _zenodo_record(sim.MS2DEEPSCORE_MODEL_URLS[0])
    if model.is_file() and model.stat().st_size > 0:
        return _entry(version, "python-package", f"model {model.name} (Zenodo record {record}), cached {_day(model)}")
    return _entry(version, "python-package", f"model not downloaded yet (Zenodo record {record}); it is fetched on first use")


def _spec2vec():
    # mzDesk scores Spec2Vec straight from the pretrained gensim Word2Vec model; the spec2vec package is not used
    from mzdesk.ms2 import ms2_similarity as sim
    gensim = _dist_version("gensim", module="gensim")
    if gensim is None:
        return _absent("model", "not installed (needs gensim)")
    files = [Path(sim.CACHE_DIR) / name for name in sim.SPEC2VEC_MODEL_FILES]
    main = next(f for f in files if f.name.endswith(".model"))
    if not all(f.is_file() and f.stat().st_size > 0 for f in files):
        return _absent("model", "not installed (the model is fetched on first use)")
    record = _zenodo_record(next(iter(sim.SPEC2VEC_MODEL_FILES.values())))
    name = main.name.removeprefix("spec2vec_").removesuffix(".model")
    return _entry(name, "model", f"pretrained Word2Vec model (Zenodo record {record}), scored with gensim {gensim}, cached {_day(main)}")


def _lipid_plus():
    from mzdesk import lipid_name_map
    return _entry(lipid_name_map.__version__, "engine", "lipid_name_map of the LIPID+ project, bundled with mzDesk")


# ---------------------------------------------------------------- engines
def _msfinder():
    from mzdesk.ms2 import msfinder_engine as me
    s = me.status()
    if not s["exe"] or not s["resources"]:
        return _absent("engine")
    if not s["mono"]:
        return _absent("engine", "not installed (Mono not found)")
    found = re.search(r"(\d+(?:\.\d+)+)$", s["dir"].name)                       # the release folder is called MSFINDER_v3.73
    dbs = []
    for label, pattern in (("structure DB", r"MsfinderStructureDB-(VS\d+)\.esd"), ("formula DB", r"MsfinderFormulaDB-(VS\d+)\.efd")):
        for f in sorted(s["resources"].iterdir()):
            hit = re.fullmatch(pattern, f.name)
            if hit:
                dbs.append(f"{label} {hit.group(1)}")
                break
    mono = _mono_version()
    detail = "RIKEN MS-FINDER console build" + (f"; {', '.join(dbs)}" if dbs else "") + (f"; Mono {mono}" if mono else "")
    return _entry(found.group(1) if found else me.VERSION, "engine", detail)


def _frigid():
    from mzdesk.ms2 import frigid_engine as fe
    try:
        s = fe.check_frigid()
    except fe.FrigidError:
        return _absent("engine")
    commit = _git_head(s["clone"]) or fe.FRIGID_COMMIT
    dist_info = sorted(Path(s["dir"], "env").glob("lib/python3*/site-packages/torch-*.dist-info"))      # its own PyTorch, in its own environment
    torch = re.fullmatch(r"torch-(.+)\.dist-info", dist_info[0].name) if dist_info else None
    detail = "FRIGID-base (DLM and MIST checkpoints, SAFE length model trained for mzDesk)" + (f"; PyTorch {torch.group(1)}" if torch else "")
    return _entry(commit[:7], "engine", detail)


def _cfmid_part(group):
    """CFM-ID (group 1) or MSRB (group 2): both versions are in the tag of the Docker image that mzDesk runs."""
    from mzdesk.ms2 import cfmid_engine as ce
    state, info = _docker_image()
    if state != "ok":
        return _absent("engine", f"not installed ({info})")
    found = re.search(r"msml([\d.]+)_msrb([\d.]+)$", ce.IMAGE)
    if not found:
        return _entry(None, "engine", "unknown")
    return _entry(found.group(group), "engine", f"Docker image {ce.IMAGE}" + (f" (image id {info})" if info else ""))


def _cfmid():
    return _cfmid_part(1)


def _msrb():
    return _cfmid_part(2)


def _metfrag():
    from mzdesk.ms2 import metfrag_engine
    jar = Path(metfrag_engine.JAR)
    if not jar.is_file():
        return _absent("engine")
    java = _java_version()
    if not java:
        return _absent("engine", "not installed (Java not found)")
    found = re.fullmatch(r"MetFragCommandLine-([\d.]+)\.jar", jar.name)
    return _entry(found.group(1), "engine", f"bundled {jar.name}; runs on Java {java}") if found else _entry(None, "engine", "unknown")


def _mono():
    version = _mono_version()
    return _entry(version, "engine", "runtime for MS-FINDER") if version else _absent("engine")


def _java():
    version = _java_version()
    return _entry(version, "engine", "runtime for MetFrag") if version else _absent("engine")


# ---------------------------------------------------------------- models
def _neims():
    from mzdesk.ms2 import neims_setup
    weights = Path(os.environ.get("MZDESK_NEIMS_WEIGHTS") or neims_setup.default_weights())
    if not weights.is_file():
        return _absent("model")
    return _entry(neims_setup.CHECKPOINT, "model", f"weights of Wei et al. 2019 converted to numpy ({weights.name}), file dated {_day(weights)}")


def _molgpka():
    from mzdesk.platform_fit import pka_setup
    weights = Path(os.environ.get("MZDESK_MOLGPKA_WEIGHTS") or pka_setup.default_weights())
    if not weights.is_file():
        return _absent("model")
    commit = pka_setup.COMMIT[:7]
    return _entry(commit, "model", f"weights of Xundrug/MolGpKa at that commit, converted to numpy ({weights.name}), file dated {_day(weights)}")


# ---------------------------------------------------------------- data
def _library(name):
    """Entry of one local spectral library, from the meta.json its index carries (the index itself is not opened)."""
    from mzdesk.ms2 import library_index as lx
    meta_file = lx.index_dir(name) / "meta.json"
    if not meta_file.is_file():
        return _absent("data")
    meta = json.loads(meta_file.read_text(encoding="utf-8"))
    if meta.get("format") != lx.FORMAT:
        return _absent("data", "not installed (the index has an old format)")
    version, built, released = str(meta.get("version") or ""), str(meta.get("built") or "")[:10], str(meta.get("date") or "")
    if not version:
        return _entry(None, "data", "unknown")
    detail = f"{int(meta['n_spectra']):,} spectra" + (f", built {built}" if built else "") + (f", source dated {released}" if released and released != version else "")
    return _entry(version, "data", detail)


def _massbank():
    return _library("MassBank")


def _gnps():
    return _library("GNPS")


def _libraries():
    """MassBank and GNPS together (the MS/MS Annotation workflow searches both through one index)."""
    found = {name: _library(name) for name in ("MassBank", "GNPS")}
    parts = [f"{name} {e['version']}" for name, e in found.items() if e["version"]]
    details = [f"{name}: {e['detail']}" for name, e in found.items()]
    return _entry("; ".join(parts) or None, "data", "; ".join(details) if parts else "not installed")


def _ccsbase():
    from mzdesk.orthogonal import ccs
    folder = ccs.ccs_dir()
    if not ccs.installed(folder):
        return _absent("data")
    meta = json.loads((folder / ccs.META_FILE).read_text(encoding="utf-8"))
    commit = str(meta.get("commit") or "")[:7]
    if not commit:
        return _entry(None, "data", "unknown")
    sklearn = (meta.get("model") or {}).get("sklearn")
    if (folder / ccs.MODEL_FILE).is_file():
        model_text = "; prediction model trained here" + (f" with scikit-learn {sklearn}" if sklearn else "")
    else:
        model_text = "; no prediction model"
    detail = (f"dylanhross/c3sdb at that commit: {int(meta['rows_kept']):,} CCS values, {int(meta['compounds']):,} compounds, "
              f"{int(meta['datasets'])} datasets, built {meta.get('built', '?')}{model_text}")
    return _entry(commit, "data", detail)


def _report():
    from mzdesk.methods import report
    folder = report.report_dir()
    if not report.installed(folder):
        return _absent("data")
    payload = json.loads((folder / report.DATASETS_FILE).read_text(encoding="utf-8"))
    commit = str(payload.get("commit") or "")[:7]
    if not commit:
        return _entry(None, "data", "unknown")
    return _entry(commit, "data", f"michaelwitting/RepoRT at that commit, {len(payload.get('datasets') or []):,} datasets, indexed {_day(folder / report.DATASETS_FILE)}")


def _jp18():
    from mzdesk.methods import jp18
    if not jp18.DATA_FILE.is_file():
        return _absent("data")
    with gzip.open(jp18.DATA_FILE, "rt", encoding="utf-8") as f:
        n = sum(1 for line in f if line.strip())
    return _entry("JP18", "data", f"Japanese Pharmacopoeia 18th edition with Supplements I and II, {n:,} LC/GC methods, bundled with mzDesk")


def _nist_asd():
    return _entry(None, "data", "unknown (first ionization energies copied into mzDesk; the ASD release is not recorded)")


def _pan_redu():
    from mzdesk.gnps import compound_db_check as g
    # redu_path() would download the whole table when nothing is cached: look at the two files only
    path = g.REDU_SLIM if os.path.isfile(g.REDU_SLIM) else g.BUNDLED_REDU if os.path.isfile(g.BUNDLED_REDU) else None
    return _online(f"{ONLINE}; the Dataset Finder page reads a local metadata snapshot dated {g.redu_cache_date(path)}" if path else ONLINE)


# ---------------------------------------------------------------- the list
def _package(*dists, module=None, detail=None):
    return lambda: python_package(*dists, module=module, detail=detail)


# name -> (kind, resolver). The first block is every name that a page's `based_on` uses (tests/test_versions.py checks it against spec.py);
# "GNPS (online)" is the online GNPS2 services: build.ts uses it for "GNPS" on a page that does not list the MassBank library beside it.
RESOLVERS = {
    "mzDesk": ("mzdesk", _mzdesk),
    "RDKit": ("python-package", _package("rdkit", "rdkit-pypi", module="rdkit")),
    "IsoSpecPy": ("python-package", _package("IsoSpecPy", module="IsoSpecPy")),
    "matchms": ("python-package", _package("matchms", module="matchms")),
    "ms_entropy": ("python-package", _package("ms_entropy", module="ms_entropy")),
    "pyOpenMS": ("python-package", _package("pyopenms", module="pyopenms")),
    "BUDDY": ("python-package", _package("msbuddy", module="msbuddy", detail="Python package msbuddy")),
    "MS2DeepScore": ("python-package", _ms2deepscore),
    "Spec2Vec": ("model", _spec2vec),
    "MassQL": ("engine", _massql),
    "MS-FINDER": ("engine", _msfinder),
    "FRIGID": ("engine", _frigid),
    "CFM-ID": ("engine", _cfmid),
    "MSRB": ("engine", _msrb),
    "MetFrag": ("engine", _metfrag),
    "LIPID+": ("engine", _lipid_plus),
    "NEIMS": ("model", _neims),
    "MolGpKa": ("model", _molgpka),
    "MassBank": ("data", _massbank),
    "GNPS": ("data", _gnps),
    "MassBank/GNPS": ("data", _libraries),
    "CCSbase": ("data", _ccsbase),
    "RepoRT": ("data", _report),
    "Japanese Pharmacopoeia": ("data", _jp18),
    "NIST ASD": ("data", _nist_asd),
    "fastMASST": ("online", _online),
    "Pan-ReDU": ("online", _pan_redu),
    "PubChem": ("online", _online),
    "UniChem": ("online", _online),
    "LOTUS": ("online", _online),
    "ClassyFire": ("online", _online),
    # not named by a page, but the results depend on them
    "GNPS (online)": ("online", _online),
    "Mono": ("engine", _mono),
    "Java": ("engine", _java),
    "numpy": ("python-package", _package("numpy", module="numpy")),
    "scipy": ("python-package", _package("scipy", module="scipy")),
    "scikit-learn": ("python-package", _package("scikit-learn", module="sklearn")),
    "gensim": ("python-package", _package("gensim", module="gensim")),
    "torch": ("python-package", _package("torch", module="torch")),
    "plotly": ("python-package", _package("plotly", module="plotly")),
    "matplotlib": ("python-package", _package("matplotlib", module="matplotlib")),
}


def _text(value, limit):
    return " ".join(strip_paths(value).split())[:limit]


def _guard(kind, resolve):
    """The resolver's entry in the fixed shape, with no path in it; {"version": None, "detail": "unknown"} when it fails for any reason."""
    try:
        entry = resolve()
        version = entry.get("version")
        version = _text(version, 120) if version not in (None, "") else None
        detail = _text(entry["detail"], 300) if entry.get("detail") else None
        kind = entry.get("kind") if entry.get("kind") in KINDS else kind
        return _entry(version, kind, detail or ("unknown" if version is None else None))
    except Exception:                   # a probe must never break the report
        return _entry(None, kind, "unknown")


def _build():
    for probe in (_docker_image, _mono_version, _java_version):
        probe.cache_clear()
    return {"mzdesk": __version__, "python": platform.python_version(),
            "checked_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "components": {name: _guard(kind, resolve) for name, (kind, resolve) in RESOLVERS.items()}}


_LOCK = threading.Lock()
_CACHE = {"at": 0.0, "data": None}


def environment():
    """The body of GET /api/versions (see the top of this file). Cached for TTL_S seconds; the caller gets its own copy."""
    with _LOCK:
        now = time.monotonic()
        if _CACHE["data"] is None or now - _CACHE["at"] > TTL_S:
            _CACHE["data"] = _build()
            _CACHE["at"] = now
        return copy.deepcopy(_CACHE["data"])


def reset():
    """Forget the cached answer (the next call looks again). For tests."""
    with _LOCK:
        _CACHE["data"], _CACHE["at"] = None, 0.0
