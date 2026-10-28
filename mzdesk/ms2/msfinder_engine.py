"""MS-FINDER engine for the Structure prediction page: formula and structure candidates from one MS/MS spectrum.

Runs the official MS-FINDER v3.73 release (RIKEN, Tsugawa lab) through Mono, on this machine only. The release is the Windows
GUI build; `python -m mzdesk.ms2.msfinder_setup` adds a console front end compiled from the project's own sources (nothing patched).
One run = one spectrum: write a .mat query file into a fresh temp folder, call
    mono MsfinderConsoleApp.exe predict -i <in> -o <out> -m mzdesk_msfinder.ini
and read `Formula result-N.txt` / `Structure result-N.txt` (N is a sum of date parts, so the files are found by glob).
The shipped parameter file (mzdesk/lib/msfinder/mzdesk_msfinder.ini) keeps every online source off: PubChem is never queried, the
structure candidates come from MS-FINDER's local databases only. It also keeps F, Cl, Br, I and Si off; "Include halogens" and
"Include Si" on the page switch them on for ONE run through a copy of the file in that run's temp folder (the shipped file is never changed).

Things worth knowing (all checked against the sources of tag MSFINDER-v3.73 and real runs):
  * Exit code is 1 on success and 255 on failure, so success is judged by the result files.
  * PRECURSORTYPE must be written; without it the reader silently uses "[M]+." and the formulas are wrong.
  * The number reader of the query file knows no exponent except "e+": peaks are written as plain decimals.
  * MS-FINDER writes .fgt / .sfd side files into the input folder, and resolves Resources/ from the exe folder.
  * A run loads about 2.3 GB of databases and takes 4-5 s, so runs are serialised and watched (time and memory).
Top-level imports are the standard library only (spec.py reads ADDUCTS).
"""
import os
import re
import shutil
import signal
import subprocess
import tempfile
import threading
import time
from pathlib import Path

VERSION = "3.73"
EXE_NAME = "MsfinderConsoleApp.exe"
DEFAULT_DIR = Path.home() / "tools" / "msfinder" / f"MSFINDER_v{VERSION}"
INI_PATH = Path(__file__).resolve().parent.parent / "lib" / "msfinder" / "mzdesk_msfinder.ini"
SETUP_COMMAND = "python -m mzdesk.ms2.msfinder_setup"

TIMEOUT_S = 300             # one run is 4-5 s; a very large spectrum can take longer
RSS_CAP_MB = 4096           # a run peaks at about 2.3 GB; anything near double is a runaway
POLL_S = 0.5
BUSY_WAIT_S = 20            # how long a second request waits for the run in progress before it is turned away
MAX_PEAKS = 2000

# MS-FINDER's own adduct names (Resources/AdductPositives.apf, AdductNegatives.anf): the common singly charged ones, each checked with a
# real run. FA = formate, Hac = acetate. [M]+ is the permanent cation (the formula is the cation's). "Auto" on the page = the first entry.
ADDUCTS = {"Positive": ["[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M+Li]+", "[M+H-H2O]+", "[M+ACN+H]+", "[M]+"],
           "Negative": ["[M-H]-", "[M+Cl]-", "[M+FA-H]-", "[M+Hac-H]-", "[M+Br]-", "[M-H2O-H]-", "[M+Na-2H]-"]}
# The isotope peaks are matched to the precursor within half of Ms1Tolerance of the parameter file (0.001 Da): 0.5 mDa.
MS1_MATCH_DA = 0.0005
# Element switches of the parameter file (MS-FINDER's INI reader takes a value that contains "T" as on, anything else as off).
# CommonRange (the shipped setting) allows, per carbon atom: F 1.5, Cl 1, Br 1, I 0.5, Si 0.5 (and H 4, N 2, O 2.5, P 0.5, S 1; O/P >= 2 when P is present).
HALOGEN_KEYS = ("Fcheck", "ClCheck", "BrCheck", "Icheck")
SILICON_KEYS = ("SiCheck",)

_RUN_LOCK = threading.Lock()


class MsFinderError(RuntimeError):
    """MS-FINDER is missing, busy, or failed; log is the tail of its output."""

    def __init__(self, message, log=""):
        super().__init__(message)
        self.log = log


# ---------------------------------------------------------------- installation
def install_dir():
    return Path(os.environ.get("MZDESK_MSFINDER_DIR") or DEFAULT_DIR).expanduser()


def find_mono():
    """Path of the mono executable (MZDESK_MONO, then PATH, then Homebrew) or None."""
    for cand in (os.environ.get("MZDESK_MONO"), shutil.which("mono"), "/opt/homebrew/bin/mono"):
        if cand and Path(cand).is_file() and os.access(cand, os.X_OK):
            return str(cand)
    return None


def status():
    """What is in place: {"mono", "dir", "exe", "resources", "ini"} with paths or None (never starts MS-FINDER)."""
    d = install_dir()
    exe = d / EXE_NAME
    res = d / "Resources"
    has_db = res.is_dir() and any(res.glob("MsfinderStructureDB-*")) and any(res.glob("MsfinderFormulaDB-*"))
    return {"mono": find_mono(), "dir": d, "exe": exe if exe.is_file() else None,
            "resources": res if has_db else None, "ini": INI_PATH if INI_PATH.is_file() else None}


def check_msfinder():
    """Return (mono, exe, ini) or raise MsFinderError with a message that says what to do."""
    s = status()
    if not s["mono"]:
        raise MsFinderError("Structure prediction runs MS-FINDER through Mono, which is not installed here. The module is part of the local "
                            "version of mzDesk only (the online version has no Mono). On your own machine: install Mono (for example "
                            f"`brew install mono`), then run `{SETUP_COMMAND}`.")
    if not s["exe"] or not s["resources"]:
        missing = f"{EXE_NAME} is" if not s["exe"] else "the MS-FINDER databases (Resources/) are"
        raise MsFinderError(f"MS-FINDER is not installed: {missing} missing in {s['dir']}. "
                            f"Run once in a terminal: {SETUP_COMMAND} (MZDESK_MSFINDER_DIR changes the folder).")
    if not s["ini"]:
        raise MsFinderError(f"The mzDesk parameter file {INI_PATH} is missing; reinstall mzDesk.")
    return s["mono"], s["exe"], s["ini"]


# ---------------------------------------------------------------- parameter file
def ini_with_elements(text, halogen=False, silicon=False):
    """Text of the parameter file with the chosen element switches set to True; every other line stays byte for byte.

    halogen -> Fcheck, ClCheck, BrCheck, Icheck; silicon -> SiCheck (independent of each other). A switch that is not in the text is
    appended, so an older copy of the file still works. Neither flag: the text comes back unchanged.
    """
    wanted = {k.lower(): k for flag, keys in ((halogen, HALOGEN_KEYS), (silicon, SILICON_KEYS)) if flag for k in keys}
    done, out = set(), []
    for line in text.splitlines(keepends=True):
        key, eq, _ = line.partition("=")
        name = wanted.get(key.strip().lower()) if eq else None         # "#Fcheck=..." is a comment: its key starts with "#"
        if name:
            eol = line[len(line.rstrip("\r\n")):]                      # keep the line ending as it was
            out.append(f"{name}=True{eol}")
            done.add(name)
        else:
            out.append(line)
    extra = [f"{name}=True\n" for name in wanted.values() if name not in done]
    if extra and out and not out[-1].endswith("\n"):
        out[-1] += "\n"
    return "".join(out + extra)


# ---------------------------------------------------------------- query file
def resolve_adduct(positive, adduct):
    """Adduct name for MS-FINDER; None / "Auto" = [M+H]+ (positive) or [M-H]- (negative). ValueError when it does not fit the mode."""
    pool = ADDUCTS["Positive" if positive else "Negative"]
    if not adduct or adduct == "Auto":
        return pool[0]
    if adduct not in pool:
        other = ADDUCTS["Negative" if positive else "Positive"]
        raise ValueError(f"{adduct} is not a {'positive' if positive else 'negative'}-mode adduct MS-FINDER can use here."
                         if adduct in other else f"MS-FINDER does not support {adduct} here. Supported: {', '.join(pool)}.")
    return adduct


def _peak_lines(peaks):
    """Plain-decimal 'm/z<TAB>intensity' lines, base peak = 1000 (MS-FINDER normalises anyway). Peaks without intensity are dropped."""
    peaks = [(float(mz), float(i)) for mz, i in peaks if i > 0 and mz > 0]
    top = max((i for _, i in peaks), default=0.0)
    return [f"{mz:.5f}\t{i / top * 1000:.4f}" for mz, i in sorted(peaks)]


def build_query(precursor_mz, positive, adduct, ms2, ms1=None, name="query"):
    """Text of the .mat query file. ms2 / ms1: lists of (m/z, intensity); ms1 = isotope peaks around the precursor (optional)."""
    ms2_lines = _peak_lines(ms2)
    if not ms2_lines:
        raise ValueError("The MS2 spectrum has no peaks.")
    head = [f"NAME: {re.sub(r'[^A-Za-z0-9 _.-]', '_', name) or 'query'}", f"PRECURSORMZ: {precursor_mz:.5f}", f"PRECURSORTYPE: {adduct}",
            f"IONMODE: {'Positive' if positive else 'Negative'}", "SPECTRUMTYPE: Centroid"]
    ms1_lines = _peak_lines(ms1 or [])
    if ms1_lines:
        head += ["MSTYPE: MS1", f"Num Peaks: {len(ms1_lines)}"] + ms1_lines
    return "\n".join(head + ["MSTYPE: MS2", f"Num Peaks: {len(ms2_lines)}"] + ms2_lines) + "\n"


# ---------------------------------------------------------------- result files
def _num(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def _table(text):
    """Rows of a tab-separated result file as dicts keyed by the header line."""
    lines = [ln.rstrip("\r\n") for ln in text.splitlines() if ln.strip()]
    if not lines:
        return []
    header = lines[0].split("\t")
    rows = []
    for ln in lines[1:]:
        cells = ln.split("\t")
        cells += [""] * (len(header) - len(cells))
        rows.append(dict(zip(header, cells)))
    return rows


def _databases(text):
    """'HMDB=HMDB0001847,FooDB=FDB1;FDB2,PubChem=2519' -> [('HMDB', 'HMDB0001847'), ('FooDB', 'FDB1;FDB2'), ('PubChem', '2519')]."""
    out = []
    for part in (text or "").split(","):
        name, eq, ident = part.partition("=")
        if eq:
            out.append((name.strip(), ident.strip()))
        elif part.strip() and out:                      # an ID that itself contains a comma
            out[-1] = (out[-1][0], out[-1][1] + "," + part.strip())
    return out


def _by_score(row):
    return float("inf") if row["score"] is None else -row["score"]          # best first, unscored last (the sort is stable)


def _counts(rows):
    first = rows[0] if rows else {}
    ms1, ms2 = _num(first.get("MS1 count")), _num(first.get("MSMS count"))
    return (None if ms1 is None else int(ms1)), (None if ms2 is None else int(ms2))


def parse_formula_results(text):
    """-> {"rows": [{"formula", "neutral_mass", "mass_error_da", "score"}], "ms1_count", "ms2_count"}; best score first."""
    rows = _table(text)
    out = [{"formula": r.get("Formula", ""), "neutral_mass": _num(r.get("Theoretical mass")), "mass_error_da": _num(r.get("Mass error")),
            "score": _num(r.get("Formula score"))} for r in rows if r.get("Formula")]
    out.sort(key=_by_score)
    ms1, ms2 = _counts(rows)
    return {"rows": out, "ms1_count": ms1, "ms2_count": ms2}


def parse_structure_results(text):
    """-> {"rows": [{"name", "score", "databases": [(db, id)], "formula", "ontology", "inchikey", "smiles"}], ...}; best score first.

    The column called "Structure rank" holds the candidate NAME (names with a comma may come wrapped in quotes).
    """
    rows = _table(text)
    out = []
    for r in rows:
        name = r.get("Structure rank", "").strip()
        if len(name) > 1 and name[0] == name[-1] == '"':
            name = name[1:-1]
        if not (r.get("SMILES") or r.get("InChIKey")):
            continue
        out.append({"name": name, "score": _num(r.get("Total score")), "databases": _databases(r.get("Databases")),
                    "formula": r.get("Formula", ""), "ontology": r.get("Ontology", "").strip(),
                    "inchikey": r.get("InChIKey", "").strip(), "smiles": r.get("SMILES", "").strip()})
    out.sort(key=_by_score)
    ms1, ms2 = _counts(rows)
    return {"rows": out, "ms1_count": ms1, "ms2_count": ms2}


# ---------------------------------------------------------------- process
def _group_rss_mb(pgid):
    """Resident memory (MB) of every process in the process group, from `ps`; 0 when ps is not available."""
    try:
        out = subprocess.run(["ps", "-axo", "pgid=,rss="], capture_output=True, text=True, timeout=5).stdout
    except (OSError, subprocess.TimeoutExpired):
        return 0.0
    total = 0
    for line in out.splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0] == str(pgid) and parts[1].isdigit():
            total += int(parts[1])
    return total / 1024.0


def _kill_group(proc):
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        pass
    try:
        proc.wait(timeout=10)
    except subprocess.TimeoutExpired:
        pass


def run_watched(cmd, cwd, env, timeout=TIMEOUT_S, rss_cap_mb=RSS_CAP_MB, poll_s=POLL_S, label="MS-FINDER"):
    """Run ``cmd`` in its own process group, stdout+stderr to a temp file -> (exit code, output tail, peak RSS in MB).

    Past ``timeout`` seconds or ``rss_cap_mb`` of resident memory the whole group is killed and MsFinderError raised (with the output so far).
    """
    with tempfile.TemporaryFile(mode="w+b") as log_file:
        proc = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=log_file, stderr=subprocess.STDOUT, start_new_session=True)

        def tail():
            log_file.seek(0)
            return log_file.read().decode("utf-8", "replace").strip()[-4000:]

        peak, deadline = 0.0, time.monotonic() + timeout
        try:
            while True:
                try:
                    proc.wait(timeout=poll_s)
                    break
                except subprocess.TimeoutExpired:
                    pass
                rss = _group_rss_mb(proc.pid)
                peak = max(peak, rss)
                if rss > rss_cap_mb:
                    _kill_group(proc)
                    raise MsFinderError(f"{label} used more than {rss_cap_mb / 1024:g} GB of memory and was stopped. Try a smaller spectrum.", tail())
                if time.monotonic() > deadline:
                    _kill_group(proc)
                    raise MsFinderError(f"{label} did not finish within {timeout:g} s and was stopped.", tail())
        finally:
            if proc.poll() is None:                     # an exception other than the two above (for example KeyboardInterrupt)
                _kill_group(proc)
        return proc.returncode, tail(), peak


def _read(path):
    return path.read_text(encoding="utf-8", errors="replace")       # the files are ASCII: MS-FINDER writes '?' for any other character


def predict(precursor_mz, positive, adduct, ms2, ms1=None, *, halogen=False, silicon=False, name="query", timeout=TIMEOUT_S, rss_cap_mb=RSS_CAP_MB):
    """Formula and structure candidates for one MS/MS spectrum.

    adduct: MS-FINDER name or None / "Auto". ms2 / ms1: [(m/z, intensity)]. halogen / silicon: also search formulas with F, Cl, Br, I /
    with Si (a per-run copy of the parameter file; both off = the shipped file is used as it is). Raises ValueError for input MS-FINDER
    cannot use and MsFinderError when it is missing, busy or fails. Returns
    {"adduct", "precursor_mz", "positive", "halogen", "silicon", "formulas": [...], "structures": [...], "ms1_count", "ms2_count",
     "seconds", "peak_rss_mb", "log"}.
    """
    adduct = resolve_adduct(positive, adduct)
    if not precursor_mz or precursor_mz <= 0:
        raise ValueError("Precursor m/z must be positive.")
    if len(ms2) > MAX_PEAKS or len(ms1 or []) > MAX_PEAKS:
        raise ValueError(f"MS-FINDER takes at most {MAX_PEAKS} peaks per spectrum; clean the spectrum first.")
    query = build_query(precursor_mz, positive, adduct, ms2, ms1, name)
    mono, exe, ini = check_msfinder()

    if not _RUN_LOCK.acquire(timeout=BUSY_WAIT_S):
        raise MsFinderError("MS-FINDER is busy with another request (one run at a time, each needs about 2.3 GB of memory). Try again in a few seconds.")
    try:
        t0 = time.monotonic()
        with tempfile.TemporaryDirectory(prefix="mzdesk-msfinder-") as tmp:
            in_dir, out_dir = Path(tmp) / "in", Path(tmp) / "out"
            in_dir.mkdir()
            out_dir.mkdir()
            (in_dir / "query.mat").write_text(query, encoding="ascii")
            if halogen or silicon:                                  # this run only: a copy in the run's temp folder, the shipped file stays as it is
                run_ini = Path(tmp) / INI_PATH.name
                run_ini.write_bytes(ini_with_elements(Path(ini).read_bytes().decode("utf-8"), halogen, silicon).encode("utf-8"))
                ini = run_ini
            env = {**os.environ, "LC_ALL": "en_US.UTF-8", "LANG": "en_US.UTF-8"}
            cmd = [mono, str(exe), "predict", "-i", str(in_dir), "-o", str(out_dir), "-m", str(ini)]
            _, log, peak = run_watched(cmd, str(exe.parent), env, timeout, rss_cap_mb)         # exit code 1 = success: judge by the files
            formula_files = sorted(out_dir.glob("Formula result-*.txt"))
            structure_files = sorted(out_dir.glob("Structure result-*.txt"))
            if not formula_files:
                raise MsFinderError("MS-FINDER finished without writing a result file.", log)
            formulas = parse_formula_results(_read(formula_files[0]))
            structures = parse_structure_results(_read(structure_files[0])) if structure_files else {"rows": [], "ms1_count": None, "ms2_count": None}
        seconds = time.monotonic() - t0
    finally:
        _RUN_LOCK.release()

    ms1_count = formulas["ms1_count"] if formulas["ms1_count"] is not None else structures["ms1_count"]
    ms2_count = formulas["ms2_count"] if formulas["ms2_count"] is not None else structures["ms2_count"]
    return {"adduct": adduct, "precursor_mz": precursor_mz, "positive": positive, "halogen": bool(halogen), "silicon": bool(silicon),
            "formulas": formulas["rows"], "structures": structures["rows"],
            "ms1_count": ms1_count if ms1_count is not None else len(ms1 or []), "ms2_count": ms2_count if ms2_count is not None else len(ms2),
            "seconds": seconds, "peak_rss_mb": peak, "log": log}
