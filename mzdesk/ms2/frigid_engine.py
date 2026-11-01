"""FRIGID-base engine for the De novo structure page: generate candidate structures for one MS/MS spectrum with a known molecular formula.

FRIGID (Coley group, github.com/coleygroup/FRIGID, Apache-2.0) predicts a molecular fingerprint from the spectrum (a MIST encoder, MIT) and lets a
masked diffusion language model write SMILES for that formula and fingerprint. Positive ion mode only. The result is a set of hypotheses, not an
identification: the score is the Tanimoto similarity of a generated structure to the *predicted* fingerprint.

How it runs: on this machine only. `python -m mzdesk.ms2.frigid_setup` builds a Python 3.10 + PyTorch environment, a pinned clone of FRIGID and the two
checkpoints in ~/tools/frigid (MZDESK_FRIGID_DIR). A run starts the driver that ships with mzDesk (mzdesk/ms2/frigid_driver/frigid_predict.py, FRIGID's
own code plus a small sub-formula assignment) with that environment's Python, in its own process group:
    <dir>/env/bin/python -I frigid_predict.py --input job.json --output result.json
job.json / result.json sit in a fresh temp folder. The driver files are byte-for-byte the verified ones, so they find the install through the
environment variables that child_env() sets (FRIGID_ROOT, FRIGID_MIST_CKPT, FRIGID_DLM_CKPT, HF_HOME) and the length model through the job's "length_model". The driver chooses the device itself (device "auto": Apple GPU (MPS) only when the Mac has
>= 12 GB available and swap to spare, else CPU); MZDESK_FRIGID_DEVICE = auto | mps | cpu overrides it. The SAFE token-length model that the
driver uses (mzdesk/lib/frigid/best_ngboost_MSG.joblib) was trained by us on MIT-licensed MassSpecGym labels (train split) and ships with mzDesk.

Resource rules (this runs gigabytes of PyTorch): one heavy job at a time across the whole app (mzdesk.shared.heavy_lock, shared with MS-FINDER);
a run is stopped after 300 s (Quick) / 1200 s (Full); an outer memory guard (mzdesk.shared.memguard, watching the child's process tree) kills the
run above 6.5 GB resident memory or 13 GB physical footprint (the driver's own guard sits at 6 / 12 GB and normally fires first). The driver gets only
as many CPU threads as there are idle cores (cpu_threads), so a busy machine slows a run down instead of stalling it.

Driver exit codes: 0 ok, 2 invalid input, 1 runtime failure, 3 its memory guard. Top-level imports are the standard library only (spec.py reads
ADDUCTS and MODE_LABELS); RDKit is imported inside the functions that need it.
"""
import json
import os
import re
import signal
import subprocess
import sys
import tempfile
import time
from math import prod
from pathlib import Path

from mzdesk.shared import heavy_lock, memguard

FRIGID_COMMIT = "9a3e623007dc5771b0e4488cd8c91b7e49490a59"
DEFAULT_DIR = Path.home() / "tools" / "frigid"
SETUP_COMMAND = "python -m mzdesk.ms2.frigid_setup"
HERE = Path(__file__).resolve().parent
DRIVER = HERE / "frigid_driver" / "frigid_predict.py"
LIB_DIR = HERE.parent / "lib" / "frigid"
LENGTH_MODEL = LIB_DIR / "best_ngboost_MSG.joblib"
LENGTH_FEATURES = LIB_DIR / "best_ngboost_MSG_features.json"
MPS_PATCH = LIB_DIR / "mps_patch.diff"
REQUIREMENTS = LIB_DIR / "requirements.txt"

# Checkpoint sizes (bytes): a cheap guard against a half-downloaded file at run time; the md5 sums are checked by the setup command.
DLM_SIZE = 2_676_432_167
MIST_SIZE = 289_222_437
TOKENIZER_REPO, TOKENIZER_REVISION = "datamol-io/safe-gpt", "3d5fa0988383e898d5ac5db7cd52bf715bc37061"

MODES = {"quick": "Quick (16 samples)", "full": "Full (up to 100 samples)"}
MODE_LABELS = list(MODES.values())
DEFAULT_MODE = "quick"
TIMEOUT_S = {"quick": 300, "full": 1200}
RSS_CAP_MB = 6656               # 6.5 GB resident memory of the process tree (driver's own CPU cap: 6144 MB)
FOOTPRINT_CAP_MB = 13312        # 13 GB physical footprint (driver's own MPS cap: 12288 MB); the Apple GPU's memory is not in RSS
POLL_S = 0.5
BUSY_WAIT_S = 20                # how long a second request waits for the run in progress before it is turned away
MAX_PEAKS = 2000
DEVICES = ("auto", "mps", "cpu")

# Elements the MIST/FRIGID featuriser knows (subformula_assign.VALID_ELEMENTS) and the cap on its sub-formula lattice.
ELEMENTS = ("C", "H", "As", "B", "Br", "Cl", "Co", "F", "Fe", "I", "K", "N", "Na", "O", "P", "S", "Se", "Si")
MAX_SUBSETS = 5_000_000

# Same masses the driver uses (RDKit's most common isotopes) and the electron mass it subtracts.
_H, _N, _O, _NA, _K, _E = 1.007825032, 14.003074, 15.99491462, 22.98976928, 38.96370668, 0.00054858
# mzDesk adduct name -> (name the driver wants, mass added to the neutral formula's monoisotopic mass to get the ion's m/z). Positive mode only.
# [M]+ is the permanent cation: the formula is the cation's own.
ADDUCT_TABLE = {
    "[M+H]+": ("[M+H]+", _H - _E),
    "[M+Na]+": ("[M+Na]+", _NA - _E),
    "[M+K]+": ("[M+K]+", _K - _E),
    "[M+NH4]+": ("[M+H3N+H]+", _N + 4 * _H - _E),
    "[M+H-H2O]+": ("[M-H2O+H]+", -_O - _H - _E),
    "[M]+": ("[M]+", -_E),
    "[M+H-2H2O]+": ("[M-H4O2+H]+", -2 * _O - 3 * _H - _E),
}
ADDUCTS = list(ADDUCT_TABLE)                    # "Auto" on the page = the first entry
_ALIASES = {"[M-H2O+H]+": "[M+H-H2O]+", "[M-2H2O+H]+": "[M+H-2H2O]+", "[M+H3N+H]+": "[M+NH4]+", "[M-H4O2+H]+": "[M+H-2H2O]+"}
NEGATIVE_MESSAGE = ("De Novo Structure works in positive ion mode only: FRIGID and its MIST encoder were trained on positive-mode adducts. "
                    "Choose Positive, or use Structure Prediction (MS-FINDER) for negative-mode spectra.")


class FrigidError(RuntimeError):
    """FRIGID is missing, busy, or failed; log is the tail of its output."""

    def __init__(self, message, log=""):
        super().__init__(message)
        self.log = log


# ---------------------------------------------------------------- installation
def install_dir():
    return Path(os.environ.get("MZDESK_FRIGID_DIR") or DEFAULT_DIR).expanduser()


def env_python(d=None):
    return (d or install_dir()) / "env" / "bin" / "python"


def tokenizer_file(d=None):
    """The cached SAFE tokenizer (hf_cache/hub/models--datamol-io--safe-gpt/snapshots/<revision>/tokenizer.json) or None."""
    d = d or install_dir()
    snaps = d / "hf_cache" / "hub" / "models--datamol-io--safe-gpt" / "snapshots"
    hits = sorted(snaps.glob("*/tokenizer.json")) if snaps.is_dir() else []
    return hits[0] if hits else None


def safe_init(d=None):
    """<env>/lib/python3.x/site-packages/safe/__init__.py or None."""
    hits = sorted(((d or install_dir()) / "env").glob("lib/python3*/site-packages/safe/__init__.py"))
    return hits[0] if hits else None


def mps_patch_applied(clone):
    """True when FRIGID's sampler has the local device patch (cuda -> mps -> cpu, FRIGID_DEVICE override)."""
    try:
        return "FRIGID_DEVICE" in (Path(clone) / "src" / "dlm" / "sampler.py").read_text(encoding="utf-8")
    except OSError:
        return False


def safe_patch_applied(init):
    """True when the SAFE package's __init__ is the trimmed one (FRIGID's env/fix_safe_imports.sh edit)."""
    try:
        lines = Path(init).read_text(encoding="utf-8").splitlines()
    except OSError:
        return False
    return bool(lines) and lines[-1] == "from .converter import SAFEConverter, decode, encode" and all(ln.startswith("# ") for ln in lines[:-1])


def _sized(path, size):
    return path if path.is_file() and path.stat().st_size == size else None


def status():
    """What is in place: path or None for each part (never starts FRIGID). Keys: dir, python, clone, mps_patch, safe_patch, dlm, mist, tokenizer,
    driver, length_model."""
    d = install_dir()
    clone = d / "FRIGID"
    has_clone = (clone / "src" / "dlm" / "sampler.py").is_file() and (clone / "src" / "mist").is_dir()
    py = env_python(d)
    init = safe_init(d)
    return {"dir": d, "python": py if py.is_file() and os.access(py, os.X_OK) else None,
            "clone": clone if has_clone else None, "mps_patch": has_clone and mps_patch_applied(clone),
            "safe_patch": bool(init) and safe_patch_applied(init),
            "dlm": _sized(d / "checkpoints" / "DLM.ckpt", DLM_SIZE), "mist": _sized(d / "checkpoints" / "mist_msg_aug.ckpt", MIST_SIZE),
            "tokenizer": tokenizer_file(d),
            "driver": DRIVER if DRIVER.is_file() else None, "length_model": LENGTH_MODEL if LENGTH_MODEL.is_file() and LENGTH_FEATURES.is_file() else None}


def check_frigid():
    """Return the status dict when everything is in place, else raise FrigidError with a message that says what to do."""
    s = status()
    if not s["driver"] or not s["length_model"]:
        raise FrigidError(f"The mzDesk files for De Novo Structure are missing ({DRIVER.parent} or {LIB_DIR}); reinstall mzDesk.")
    missing = [label for key, label in (("python", "the Python environment (env/)"), ("clone", "the FRIGID code (FRIGID/)"),
                                        ("mps_patch", "FRIGID's device patch"), ("safe_patch", "the SAFE import patch"),
                                        ("dlm", "the DLM checkpoint (checkpoints/DLM.ckpt)"), ("mist", "the MIST checkpoint (checkpoints/mist_msg_aug.ckpt)"),
                                        ("tokenizer", "the SAFE tokenizer (hf_cache/)")) if not s[key]]
    if missing:
        raise FrigidError("De Novo Structure runs FRIGID, which is part of the local version of mzDesk only and "
                          f"is not installed here: missing {', '.join(missing)} in {s['dir']}. "
                          f"Run {SETUP_COMMAND} (about 3 GB; MZDESK_FRIGID_DIR changes the folder).")
    return s


# ---------------------------------------------------------------- input
def resolve_adduct(adduct):
    """-> (mzDesk adduct name, name for the driver); None / "Auto" = [M+H]+. ValueError for an adduct FRIGID was not trained on."""
    if not adduct or adduct == "Auto":
        adduct = ADDUCTS[0]
    adduct = _ALIASES.get(str(adduct).replace(" ", ""), str(adduct).replace(" ", ""))
    if adduct not in ADDUCT_TABLE:
        if adduct.endswith("-") or adduct.endswith("-]"):
            raise ValueError(NEGATIVE_MESSAGE)
        raise ValueError(f"FRIGID does not support the adduct {adduct}. Supported: {', '.join(ADDUCTS)}.")
    return adduct, ADDUCT_TABLE[adduct][0]


def resolve_mode(mode):
    """'quick' | 'full' from an id or a page label ('Quick (16 samples)'); None = the default. ValueError otherwise."""
    if mode is None or mode == "":
        return DEFAULT_MODE
    for key, label in MODES.items():
        if str(mode).strip().lower() in (key, label.lower()):
            return key
    raise ValueError(f"Mode must be one of: {', '.join(MODE_LABELS)}.")


def formula_counts(formula):
    """{'C': 8, 'H': 10, ...} for a plain neutral formula; ValueError when it is not one FRIGID can take."""
    f = (formula or "").strip()
    if not f or not re.fullmatch(r"(?:[A-Z][a-z]?\d*)+", f):
        raise ValueError(f"{f or 'The formula'} is not a plain element formula such as C8H10N4O2.")
    counts = {}
    for sym, num in re.findall(r"([A-Z][a-z]?)(\d*)", f):
        counts[sym] = counts.get(sym, 0) + (int(num) if num else 1)
    counts = {el: n for el, n in counts.items() if n > 0}
    bad = [el for el in counts if el not in ELEMENTS]
    if bad:
        raise ValueError(f"FRIGID cannot use {', '.join(bad)} in {f}. Elements it knows: {' '.join(ELEMENTS)}.")
    if "C" not in counts:
        raise ValueError(f"{f} has no carbon; FRIGID generates organic molecules.")
    if prod(n + 1 for n in counts.values()) > MAX_SUBSETS:
        raise ValueError(f"{f} is too large for FRIGID's sub-formula assignment (more than {MAX_SUBSETS:,} sub-formulae).")
    return counts


def hill(counts):
    """Hill order: C, H, then the rest alphabetically."""
    order = [el for el in ("C", "H") if el in counts] + sorted(el for el in counts if el not in ("C", "H"))
    return "".join(el + (str(counts[el]) if counts[el] > 1 else "") for el in order)


def check_formula(formula):
    """Formula in Hill order, or ValueError."""
    return hill(formula_counts(formula))


def ion_mz(neutral_mass, adduct):
    """m/z of the ion of the given monoisotopic mass (a formula's) as the given mzDesk adduct (singly charged)."""
    return neutral_mass + ADDUCT_TABLE[resolve_adduct(adduct)[0]][1]


def resolve_device(device=None):
    """Device to ask the driver for: the argument, else MZDESK_FRIGID_DEVICE, else auto."""
    d = str(device or os.environ.get("MZDESK_FRIGID_DEVICE") or "auto").strip().lower()
    if d not in DEVICES:
        raise FrigidError(f"MZDESK_FRIGID_DEVICE must be one of {', '.join(DEVICES)}, not {d!r}.")
    return d


def build_job(formula, precursor_mz, adduct, peaks, mode="quick", *, device="auto", seed=42, name="query", length_model=None):
    """The job.json content. adduct: the driver's name (see resolve_adduct); peaks: [(m/z, intensity)] with positive values."""
    clean = sorted((float(mz), float(i)) for mz, i in peaks if i > 0 and mz > 0)
    return {"formula": formula, "precursor_mz": float(precursor_mz), "adduct": adduct, "ion_mode": "positive",
            "peaks": [[mz, i] for mz, i in clean], "mode": mode, "device": device, "seed": int(seed),
            "name": re.sub(r"[^A-Za-z0-9 _.-]", "_", str(name)) or "query", "length_model": str(length_model or LENGTH_MODEL)}


def cpu_threads(load=None, ncpu=None):
    """Threads PyTorch may use on the CPU: the cores that are free right now (cores minus the 1-minute load average), at least 2.

    PyTorch's OpenMP threads spin-wait for each other, so on a machine whose cores are busy with something else a run with one thread per core
    does not just slow down, it stalls (measured: a Quick run that takes 20 s on an idle Mac took more than 300 s next to a SIRIUS job, and 50 s with 2 threads).
    MZDESK_FRIGID_THREADS sets the number, an OMP_NUM_THREADS in the server's own environment is respected. None = leave PyTorch's default.
    """
    try:
        load = os.getloadavg()[0] if load is None else load
    except (OSError, AttributeError):
        return None
    ncpu = ncpu or os.cpu_count() or 1
    return max(2, min(ncpu, int(ncpu - load + 0.5)))


def child_env(d):
    """Environment of the driver process: offline Hugging Face cache inside the install folder, no user site, no inherited Python variables."""
    drop = ("HF_HUB_CACHE", "HUGGINGFACE_HUB_CACHE", "TRANSFORMERS_CACHE", "FRIGID_DEVICE", "FRIGID_TOKEN_MODEL")
    env = {k: v for k, v in os.environ.items() if not k.startswith("PYTHON") and k not in drop}
    env.update(HF_HOME=str(d / "hf_cache"), HF_HUB_OFFLINE="1", PYTHONNOUSERSITE="1", TOKENIZERS_PARALLELISM="false",
               FRIGID_ROOT=str(d / "FRIGID"), FRIGID_MIST_CKPT=str(d / "checkpoints" / "mist_msg_aug.ckpt"),
               FRIGID_DLM_CKPT=str(d / "checkpoints" / "DLM.ckpt"))
    forced = os.environ.get("MZDESK_FRIGID_THREADS", "").strip()
    threads = int(forced) if forced.isdigit() and int(forced) > 0 else (None if "OMP_NUM_THREADS" in os.environ else cpu_threads())
    if threads:
        env["OMP_NUM_THREADS"] = str(threads)
    return env


# ---------------------------------------------------------------- output
def filter_candidates(raw, adduct="[M+H]+"):
    """Keep only what the page shows: formula-matched, not isotope-labelled, not radical, right charge, one row per 2D structure.

    raw: the driver's candidates[] (it appends non-matching "padding" after the formula-matched ones). The driver's `charged` flag means
    "some atom carries a formal charge", which is also true for a nitro group, an N-oxide or a zwitterion; so the NET charge is recomputed
    from the SMILES: it must be 0, or +1 for the permanent cation [M]+ (a charged-looking formula match is expected there).
    Returns (kept rows with fresh ranks, {"not_formula_match", "isotope", "radical", "charged", "duplicate"} counts).
    """
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    want_charge = 1 if resolve_adduct(adduct)[0] == "[M]+" else 0
    dropped = {"not_formula_match": 0, "isotope": 0, "radical": 0, "charged": 0, "duplicate": 0}
    kept, seen = [], set()
    for c in sorted(raw, key=lambda c: -(c.get("score") or 0.0)):          # stable: the driver's order is already best first
        if not c.get("formula_match"):
            dropped["not_formula_match"] += 1
        elif c.get("isotope_labeled"):
            dropped["isotope"] += 1
        elif c.get("radical"):
            dropped["radical"] += 1
        elif _net_charge(Chem, c.get("smiles"), c.get("charged")) != want_charge:
            dropped["charged"] += 1
        elif c.get("inchikey_2d") and c["inchikey_2d"] in seen:
            dropped["duplicate"] += 1
        else:
            seen.add(c.get("inchikey_2d"))
            kept.append({"rank": len(kept) + 1, "smiles": c["smiles"], "formula": c.get("formula"), "score": c.get("score"),
                         "inchikey_2d": c.get("inchikey_2d")})
    return kept, dropped


def _net_charge(Chem, smiles, flagged):
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return None if flagged else 0           # unreadable: only keep it when the driver saw no charge
    return Chem.GetFormalCharge(mol)


# ---------------------------------------------------------------- process
def _kill_group(proc):
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        pass
    try:
        proc.wait(timeout=10)
    except subprocess.TimeoutExpired:
        pass


def _group_rss_mb(pgid):
    """Resident memory (MB) of the process group from `ps`; 0 when ps is not available. Stand-in for the guard on systems without libproc."""
    try:
        out = subprocess.run(["ps", "-axo", "pgid=,rss="], capture_output=True, text=True, timeout=5).stdout
    except (OSError, subprocess.TimeoutExpired):
        return 0.0
    return sum(int(p[1]) for p in (ln.split() for ln in out.splitlines()) if len(p) == 2 and p[0] == str(pgid) and p[1].isdigit()) / 1024.0


def _stop_guard(guard):
    """guard.stop() summary; zeros when the platform has no libproc (the snapshot reads the kernel's per-process numbers)."""
    try:
        return guard.stop()
    except Exception:  # noqa: BLE001
        return {"peak_rss_mb": 0.0, "peak_footprint_mb": 0.0}


def run_watched(cmd, cwd, env, timeout, rss_cap_mb=RSS_CAP_MB, footprint_cap_mb=FOOTPRINT_CAP_MB, poll_s=POLL_S, label="FRIGID"):
    """Run ``cmd`` in its own process group (stdout + stderr to a temp file) -> (exit code, output tail, {"peak_rss_mb", "peak_footprint_mb"}).

    Two memguard.MemoryGuard threads watch the child's process tree: resident memory (the CPU device) and physical footprint (the Apple GPU's
    unified memory is not in RSS). Past ``timeout`` seconds or either cap the whole group is killed and FrigidError raised (with the output so far).
    """
    with tempfile.TemporaryFile(mode="w+b") as log_file:
        proc = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=log_file, stderr=subprocess.STDOUT, start_new_session=True)
        reasons = []

        def on_exceed(reason):
            reasons.append(reason)
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                pass

        def tail():
            log_file.seek(0)
            return log_file.read().decode("utf-8", "replace").strip()[-4000:]

        guards = [memguard.MemoryGuard("cpu", pid=proc.pid, cap_mb=rss_cap_mb, poll_s=poll_s, on_exceed=on_exceed),
                  memguard.MemoryGuard("mps", pid=proc.pid, cap_mb=footprint_cap_mb, poll_s=poll_s, on_exceed=on_exceed)]
        for g in guards:
            g.start()
        peak_rss = 0.0
        deadline = time.monotonic() + timeout
        try:
            while True:
                try:
                    proc.wait(timeout=poll_s)
                    break
                except subprocess.TimeoutExpired:
                    pass
                if sys.platform != "darwin":                                       # no libproc: the guards read nothing, check RSS here
                    rss = _group_rss_mb(proc.pid)
                    peak_rss = max(peak_rss, rss)
                    if rss > rss_cap_mb:
                        on_exceed(f"process-group RSS {rss:.0f} MB > cap {rss_cap_mb} MB")
                if reasons:
                    _kill_group(proc)
                    break
                if time.monotonic() > deadline:
                    _kill_group(proc)
                    raise FrigidError(f"{label} did not finish within {timeout:g} s and was stopped.", tail())
        finally:
            if proc.poll() is None:                     # an exception other than the ones above (for example KeyboardInterrupt)
                _kill_group(proc)
            snaps = [_stop_guard(g) for g in guards]
        peak_rss = max([peak_rss] + [s["peak_rss_mb"] for s in snaps])
        peak_fp = max([0.0] + [s["peak_footprint_mb"] for s in snaps])
        if reasons:
            raise FrigidError(f"{label} used too much memory and was stopped ({reasons[0]}). Try Quick mode, a smaller spectrum, or close other applications.", tail())
        return proc.returncode, tail(), {"peak_rss_mb": peak_rss, "peak_footprint_mb": peak_fp}


def _read_result(path, log, code):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        raise FrigidError(f"FRIGID stopped without writing a result (exit code {code}).", log) from None


def parse_result(result, log="", code=0):
    """Check the driver's result JSON. Returns it when ok; raises ValueError (the driver rejected the input) or FrigidError."""
    if result.get("ok"):
        return result
    kind, msg = result.get("error_type"), str(result.get("error") or "")
    if kind == "invalid_input":
        raise ValueError(f"FRIGID rejected the input: {msg}")
    if kind == "memory_guard":
        raise FrigidError(f"FRIGID was stopped by its memory guard ({msg}). Try Quick mode, a smaller spectrum, or close other applications.", log)
    raise FrigidError(f"FRIGID failed: {kind or 'error'}: {msg}".strip(), result.get("traceback") or log)


def predict(formula, precursor_mz, adduct, peaks, mode="quick", *, name="query", seed=42, device=None, timeout=None):
    """Candidate structures for one positive-mode MS/MS spectrum with a known neutral molecular formula.

    adduct: mzDesk name (see ADDUCTS) or None / "Auto" = [M+H]+. peaks: [(m/z, intensity)]. mode: "quick" | "full" (or the page label).
    device: auto | mps | cpu (None = MZDESK_FRIGID_DEVICE, else auto). Raises ValueError for input FRIGID cannot use and FrigidError when it is
    missing, busy or fails. Returns
    {"formula", "adduct", "driver_adduct", "mode", "seed", "candidates": [{"rank", "smiles", "formula", "score", "inchikey_2d"}] (filtered),
     "dropped": {...}, "n_raw", "n_samples", "n_valid", "n_formula_matched", "n_unique_matched", "n_peaks_assigned", "device", "device_requested",
     "device_reason", "seconds", "driver_seconds", "peak_footprint_mb", "peak_rss_mb", "warnings", "length_model", "log"}.
    """
    formula = check_formula(formula)
    adduct, driver_adduct = resolve_adduct(adduct)
    mode = resolve_mode(mode)
    if not precursor_mz or precursor_mz <= 0:
        raise ValueError("Precursor m/z must be positive.")
    peaks = [(float(mz), float(i)) for mz, i in peaks if i > 0 and mz > 0]
    if not peaks:
        raise ValueError("The MS2 spectrum has no peaks.")
    if len(peaks) > MAX_PEAKS:
        raise ValueError(f"FRIGID takes at most {MAX_PEAKS} peaks per spectrum; clean the spectrum first.")
    requested = resolve_device(device)
    inst = check_frigid()
    job = build_job(formula, precursor_mz, driver_adduct, peaks, mode, device=requested, seed=seed, name=name)

    if not heavy_lock.LOCK.acquire(timeout=BUSY_WAIT_S, holder="FRIGID"):
        who = "MS-FINDER (Structure Prediction)" if heavy_lock.LOCK.holder == "MS-FINDER" else "another De Novo Structure run"
        raise FrigidError(f"FRIGID is busy: {who} is running (one heavy run at a time; a FRIGID run needs about 4-5 GB of memory). Try again in a minute.")
    try:
        t0 = time.monotonic()
        with tempfile.TemporaryDirectory(prefix="mzdesk-frigid-") as tmp:
            job_path, out_path = Path(tmp) / "job.json", Path(tmp) / "result.json"
            job_path.write_text(json.dumps(job), encoding="utf-8")
            cmd = [str(inst["python"]), "-I", str(DRIVER), "--input", str(job_path), "--output", str(out_path)]
            code, log, peaks_mb = run_watched(cmd, tmp, child_env(inst["dir"]), timeout or TIMEOUT_S[mode])
            result = parse_result(_read_result(out_path, log, code), log, code)
        seconds = time.monotonic() - t0
    finally:
        heavy_lock.LOCK.release()

    kept, dropped = filter_candidates(result.get("candidates") or [], adduct)
    mem, lm = result.get("memory") or {}, result.get("length_model")
    decision = result.get("device_decision") or {}
    return {"formula": formula, "adduct": adduct, "driver_adduct": driver_adduct, "mode": mode, "seed": job["seed"],
            "candidates": kept, "dropped": dropped, "n_raw": len(result.get("candidates") or []),
            "n_samples": result.get("n_samples_generated"), "n_valid": result.get("n_valid"),
            "n_formula_matched": result.get("n_formula_matched_with_duplicates"), "n_unique_matched": result.get("n_candidates_formula_match"),
            "n_peaks_assigned": (result.get("subformula") or {}).get("n_peaks_assigned"),
            "device": result.get("device_used"), "device_requested": requested, "device_reason": decision.get("reason"),
            "seconds": seconds, "driver_seconds": (result.get("timings_s") or {}).get("process_total_incl_imports"),
            "peak_footprint_mb": mem.get("peak_phys_footprint_mb") or peaks_mb["peak_footprint_mb"], "peak_rss_mb": peaks_mb["peak_rss_mb"],
            "warnings": list(result.get("warnings") or []), "length_model": isinstance(lm, dict) and bool(lm.get("used")), "log": log}
