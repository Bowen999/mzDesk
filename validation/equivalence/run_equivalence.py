#!/usr/bin/env python3
"""Equivalence experiment (Supplementary S5): do BUDDY, the library search, MS-FINDER and FRIGID give the same answer inside mzDesk as when
each is run directly, on the same inputs and settings, without mzDesk?

One command does everything and can be left alone (a progress bar shows where it is; the rates are printed and written at the end):

  python3 validation/equivalence/run_equivalence.py run [--out DIR] [--commit REV] [--workers 2] [--ids 1 2 51] [--limit N]
                                                         [--library-cache DIR] [--status FILE]
  python3 validation/equivalence/run_equivalence.py summary --out DIR        # rates again from the finished cases (the library
                                                                             # comparison is redone from the stored outputs)

`run` is resumable: started again on the same --out it skips the cases that have a result.json. Ctrl-C stops the running cases.
--status FILE: one line in FILE is kept current (cases done / total) for a progress board.

Steps of `run`
  1. code     the mzDesk code under test is frozen in OUT/code (git archive of mzdesk/ at --commit, default HEAD; MZDESK_CODE_ROOT overrides
              it); OUT/code/FROZEN.json names the commit, and the sha256 of this script and of msfinder_equivalence.py (whose sample draw and
              MS-FINDER comparison are reused). FRIGID's driver files of that commit are copied to OUT/direct/frigid_driver for the direct route.
  2. library  the raw MassBank 2026.03 (.msp) and GNPS 2026-03-29 (.mgf) files are read once by this script's own parser (no mzDesk code):
              every MS/MS record with a precursor m/z -> OUT/direct_library/<library>/ (precursor, ion mode, raw peaks, accession, name,
              InChIKey, InChI, SMILES). The direct library search reads only this.
  3. sample   100 MassBank 2026.03 spectra (seed 20261102, 50 positive + 50 negative compounds; msfinder_equivalence.cmd_sample, the same
              rules and draw as the v0.5.3 run of 2 Nov 2026 and the MS-FINDER run of 3 Nov 2026) -> OUT/inputs/sample.json.
  4. cases    one process per spectrum (--workers at once), MZDESK_FRIGID_DEVICE=cpu:
              mzDesk route  one ESI MS/MS Annotation run (precursor m/z, ion mode, MS2 peaks; adduct Auto; every other field at its default)
                            with Force merge on (force_merge=True: the library gate and MS-FINDER's "good enough" do not stop steps 4-5)
                            and Run FRIGID on (FRIGID: positive mode only). The run is ended when step 6 starts (steps 6-8 call no engine
                            that is tested here; step 7 would start another MS-FINDER). Every call is recorded: BUDDY (msbuddy_engine.
                            annotate_formula: inputs and output), the library search (library_search.search: inputs and hits), every
                            MsfinderConsoleApp run (query files, parameter file, result and structure files) and every FRIGID driver run
                            (job.json, the FRIGID part of its environment, result.json); and each layer as the workflow built it.
              direct route  the recorded inputs and settings, without any mzDesk code (each tool in its own process):
                            BUDDY    msbuddy (Msbuddy / MsbuddyConfig / MetaFeature) on the recorded peaks, adduct and tolerances;
                            library  ms_entropy.calculate_entropy_similarity of the recorded query against every raw MassBank and GNPS
                                     MS/MS spectrum of the same ion mode within the recorded precursor tolerance (|query - library| /
                                     library), at the recorded fragment tolerance; spectra at or above the recorded threshold are kept;
                            MS-FINDER  `mono MsfinderConsoleApp.exe predict` on the same query and parameter files, same time limit, in a
                                     private copy of the same install (msfinder_equivalence.run_direct; the install carries mzDesk's thread fix);
                            FRIGID   FRIGID's driver (the files of that commit) with FRIGID's own Python on the same job.json, the same
                                     FRIGID environment variables (thread count included) and the same time limit.
  5. summary  OUT/summary.csv (one row per case), OUT/summary.json and OUT/report.md, and the rates printed.

What is compared (the criteria of S5, Table S4)
  BUDDY      per call: the same formulas in the same order (BUDDY's top 5, what mzDesk keeps), probabilities within 1e-6.
             layer: every formula of mzDesk's formula list (formula + adduct + probability) is in the direct output of its adduct.
  library    per call: the same compounds (first 14 InChIKey characters: the record's InChIKey, else the key of its InChI, else of its
             SMILES; the name when there is no structure), best score per compound within 0.1 (mzDesk hides the GNPS copies of MassBank
             spectra, and the best-scoring copy of a compound can be a hidden one). When mzDesk kept more spectra than it returns (top 200),
             only direct compounds scoring at least mzDesk's last returned spectrum are asked for.
             layer: a candidate of mzDesk's library layer counts as found only if the direct output has the compound with the same best
             score (within 1e-6); per spectrum: every spectrum the mzDesk search returned is in the direct output with the same score.
  MS-FINDER  per run: formula table, structure table and structure files (InChIKey, SMILES, total score) the same, scores within 1e-6,
             rows sorted by score (msfinder_equivalence.compare_run). Runs that reached a time limit (run stopped, or a formula search
             >= 27 s of the 30 s limit) are counted apart: their output depends on machine speed.
             layer: every candidate of the MS-FINDER layer (InChIKey block, total score) is in a direct structure file with the same score.
  FRIGID     per run: the driver's candidate list the same (SMILES in the same order, scores within 1e-6). A run stopped by the time limit
             (either route) is counted apart. layer: every candidate of the FRIGID layer is in the direct output with the same score.

Operational burden (per spectrum, derived from the recorded calls; the counting rules of the manuscript's 2 Nov run are lost, these are
this script's): tools = tools with at least one call; tool runs = calls; files written = MS-FINDER query files + its parameter file per run,
plus FRIGID's job file per run; parameters set = BUDDY 6-8 per call (precursor m/z, charge, adduct, MS2 peaks, [MS1 peaks], instrument
preset or ppm + MS1 + MS2 tolerance, halogens), library 7 per call (libraries, metric, precursor m/z, ion mode, precursor and fragment
tolerance, threshold), MS-FINDER per run the parameter-file lines that differ from the shipped MSFINDER.INI plus the header fields of each
query file, FRIGID the keys of each job file; hand-offs = values carried from one tool's output into another's input (each MS-FINDER query
file carries a BUDDY formula and adduct, each FRIGID run a formula and adduct).
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import hashlib
import json
import math
import os
import re
import shutil
import signal
import statistics
import subprocess
import sys
import tempfile
import threading
import time
from datetime import datetime
from functools import lru_cache
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
SCRIPT = Path(__file__).resolve()
if os.environ.get("MZDESK_CODE_ROOT"):
    sys.path.insert(0, str(Path(os.environ["MZDESK_CODE_ROOT"]).resolve()))

RAW = Path.home() / ".mzdesk/data/libraries/raw"
MASSBANK_MSP = RAW / "MassBank_NISTformat.msp"
GNPS_MGF = RAW / "ALL_GNPS_NO_PROPOGATED.mgf"
REFERENCE_SAMPLE = HERE / "msfinder_20261103/inputs/sample.json"     # the draw of the MS-FINDER run (and of the v0.5.3 run)
SCORE_TOL = 1e-6
LIBRARY_TOL = 0.1                  # best score per compound; 0.05 in the v1.0.0 run (S5), widened after case 066 of the 1.2.0 run (0.05003)
BUDDY_TOP = 5                       # mzDesk keeps BUDDY's top 5 (msbuddy_engine.MAX_CANDIDATES)
CASE_TIMEOUT = 7200
FRIGID_ENV_KEYS = ("HF_HOME", "HF_HUB_OFFLINE", "PYTHONNOUSERSITE", "TOKENIZERS_PARALLELISM", "FRIGID_ROOT", "FRIGID_MIST_CKPT",
                   "FRIGID_DLM_CKPT", "OMP_NUM_THREADS")
TOOLS = ("buddy", "library", "msfinder", "frigid")
TOOL_LABELS = {"buddy": "BUDDY", "library": "Library search", "msfinder": "MS-FINDER", "frigid": "FRIGID"}
_IK = re.compile(r"[A-Z]{14}-[A-Z]{10}-[A-Z]")
_NA = ("", "N/A", "NA", "NONE", "NULL")


def stamp():
    return datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


def write_json(path, payload):
    path = Path(path)
    tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
    tmp.write_text(json.dumps(payload, indent=1, ensure_ascii=False, default=str) + "\n")
    tmp.rename(path)


def read_json(path):
    return json.loads(Path(path).read_text())


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def hms(seconds):
    seconds = int(max(0, seconds))
    return f"{seconds // 3600}h{seconds % 3600 // 60:02d}m" if seconds >= 3600 else f"{seconds // 60}m{seconds % 60:02d}s"


def _mfe():
    """msfinder_equivalence (same folder): the sample draw, the MS-FINDER direct run and comparison, the progress bar."""
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))
    import msfinder_equivalence
    return msfinder_equivalence


def direct_env():
    """Environment of a direct-route process: no mzDesk or Python variables."""
    return {k: v for k, v in os.environ.items() if not k.startswith(("MZDESK_", "PYTHON"))}


def no_mzdesk():
    assert not any(m == "mzdesk" or m.startswith("mzdesk.") for m in sys.modules), "the direct route imported mzDesk code"


# ================================================================== 1. frozen code
def freeze_code(out, commit):
    """OUT/code holds the mzdesk package of `commit` (or MZDESK_CODE_ROOT is used as it is). Returns the code root."""
    if os.environ.get("MZDESK_CODE_ROOT"):
        return Path(os.environ["MZDESK_CODE_ROOT"]).resolve()
    git = lambda *a: subprocess.run(["git", "-C", str(REPO), *a], capture_output=True, text=True, check=True).stdout.strip()   # noqa: E731
    full = git("rev-parse", f"{commit}^{{commit}}")
    code = out / "code"
    if (code / "FROZEN.json").exists():
        frozen = read_json(code / "FROZEN.json")
        if frozen["commit"] != full:
            sys.exit(f"{code} holds commit {frozen['commit'][:7]}, not {full[:7]}; use another --out or the same --commit.")
        return code
    dirty = git("status", "--porcelain", "--", "mzdesk") if full == git("rev-parse", "HEAD") else ""
    if dirty:
        print(f"WARNING: mzdesk/ has uncommitted changes; the run tests commit {full[:7]} without them:\n{dirty}", flush=True)
    tmp = out / "code.partial"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    archive = subprocess.run(["git", "-C", str(REPO), "archive", "--format=tar", full, "mzdesk", "pyproject.toml"], capture_output=True, check=True)
    subprocess.run(["tar", "-x", "-C", str(tmp)], input=archive.stdout, check=True)
    (tmp / "scripts").mkdir()
    for p in (SCRIPT, HERE / "msfinder_equivalence.py"):
        shutil.copy2(p, tmp / "scripts" / p.name)
    write_json(tmp / "FROZEN.json", {"commit": full, "subject": git("log", "-1", "--format=%s", full), "frozen_at": stamp(),
                                     "uncommitted_changes_left_out": dirty or None,
                                     "scripts": {p.name: sha256(p) for p in (SCRIPT, HERE / "msfinder_equivalence.py")}})
    tmp.rename(code)
    return code


def child_env(code):
    env = {**os.environ, "MZDESK_CODE_ROOT": str(code), "PYTHONUNBUFFERED": "1", "MZDESK_FRIGID_DEVICE": "cpu"}
    env.pop("PYTHONPATH", None)
    return env


def copy_frigid_driver(out, code):
    """OUT/direct/frigid_driver: FRIGID's driver files of the commit under test, run by the direct route outside the mzdesk package."""
    dst = out / "direct" / "frigid_driver"
    src = code / "mzdesk/ms2/frigid_driver"
    if not dst.exists():
        dst.mkdir(parents=True)
        for p in sorted(src.glob("*.py")):
            shutil.copy2(p, dst / p.name)
    return dst


# ================================================================== 2. raw libraries, read by this script (direct route)
def _na(value):
    v = str(value or "").strip()
    return "" if v.upper() in _NA else v


def _first_number(text):
    try:
        v = float(str(text).split()[0])
    except (TypeError, ValueError, IndexError):
        return None
    return v if math.isfinite(v) and v > 0 else None


def _peak_array(lines):
    import numpy as np
    rows = []
    for ln in lines:
        tok = ln.replace(",", " ").split()
        if len(tok) >= 2:
            try:
                rows.append((float(tok[0]), float(tok[1])))
            except ValueError:
                continue
    return np.array(rows, dtype=np.float64).reshape(-1, 2)


def iter_msp(path):
    """MassBank NIST-format .msp -> (fields (first value of each key), peak lines)."""
    rec, peaks, in_peaks = {}, [], False
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if not line.strip():
                if rec:
                    yield rec, peaks
                rec, peaks, in_peaks = {}, [], False
                continue
            if in_peaks:
                peaks.append(line)
                continue
            key, sep, value = line.partition(":")
            if not sep:
                continue
            key = key.strip()
            if key == "Num Peaks":
                in_peaks = True
            elif key not in rec:
                rec[key] = value.strip()
    if rec:
        yield rec, peaks


def iter_mgf(path):
    """GNPS .mgf -> (fields with upper-case keys (first value), peak lines)."""
    rec, peaks, inside = None, [], False
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not inside:
                if line == "BEGIN IONS":
                    rec, peaks, inside = {}, [], True
                continue
            if line == "END IONS":
                yield rec, peaks
                inside = False
                continue
            if not line:
                continue
            if line[0].isdigit() or line[0] in ".-":
                peaks.append(line)
            else:
                key, sep, value = line.partition("=")
                key = key.strip().upper()
                if sep and key not in rec:
                    rec[key] = value.strip()


def _mode(text, adduct=""):
    t = str(text or "").strip().lower()
    if t.startswith("pos"):
        return 1
    if t.startswith("neg"):
        return -1
    a = _na(adduct).rstrip("*")
    return 1 if a.endswith("+") else -1 if a.endswith("-") else 0


def select_massbank(rec):
    """MS/MS record with a precursor m/z -> (precursor, ion mode, meta) or None."""
    if (rec.get("Spectrum_type") or "").strip().upper() != "MS2":
        return None
    prec = _first_number(rec.get("PrecursorMZ"))
    if prec is None:
        return None
    inchi = _na(rec.get("InChI"))
    return prec, _mode(rec.get("Ion_mode"), rec.get("Precursor_type")), {
        "accession": rec.get("DB#", ""), "name": rec.get("Name", ""), "inchikey": _na(rec.get("InChIKey")), "smiles": _na(rec.get("SMILES")),
        "inchi": inchi if inchi.startswith("InChI=") else "", "adduct": _na(rec.get("Precursor_type"))}


def select_gnps(rec):
    if str(rec.get("MSLEVEL", "")).strip() != "2":
        return None
    prec = _first_number(rec.get("PEPMASS"))
    if prec is None:
        return None
    adduct = _na(rec.get("ADDUCT")) or _na(rec.get("PRECURSOR_TYPE"))
    fields = [_na(rec.get(k)).strip("\"'` ") for k in ("INCHIKEY", "INCHI", "INCHIAUX")]
    inchikey = next((v for v in fields if _IK.fullmatch(v)), "")
    inchi = next((v for v in fields[1:] if v.startswith("InChI=")), "")
    return prec, _mode(rec.get("IONMODE"), adduct), {
        "accession": rec.get("SPECTRUMID", ""), "name": _na(rec.get("NAME")), "inchikey": inchikey, "smiles": _na(rec.get("SMILES")),
        "inchi": inchi, "adduct": adduct}


RAW_LIBRARIES = {"MassBank": (MASSBANK_MSP, iter_msp, select_massbank), "GNPS": (GNPS_MGF, iter_mgf, select_gnps)}


def cmd_prepare_library(cache):
    """Read the raw files once into cache/<library>/: prec.npy, mode.npy, offsets.npy, mz.f64, intensity.f64 (raw values, peak lines as
    written), order.npy (rows by precursor), meta.jsonl + meta_offsets.npy, info.json. No mzDesk code."""
    import numpy as np
    for name, (path, reader, select) in RAW_LIBRARIES.items():
        dst = cache / name
        if (dst / "info.json").exists():
            continue
        tmp = cache / f"{name}.partial"
        shutil.rmtree(tmp, ignore_errors=True)
        tmp.mkdir(parents=True)
        t0 = time.monotonic()
        prec, mode, offsets, meta_off = [], [], [0], []
        n_read = n_no_peaks = 0
        with open(tmp / "mz.f64", "wb") as fm, open(tmp / "intensity.f64", "wb") as fi, open(tmp / "meta.jsonl", "wb") as fmeta:
            for rec, lines in reader(path):
                n_read += 1
                got = select(rec)
                if got is None:
                    continue
                peaks = _peak_array(lines)
                if not len(peaks):
                    n_no_peaks += 1
                    continue
                peaks[:, 0].tofile(fm)
                peaks[:, 1].tofile(fi)
                offsets.append(offsets[-1] + len(peaks))
                prec.append(got[0])
                mode.append(got[1])
                meta_off.append(fmeta.tell())
                fmeta.write((json.dumps(got[2], ensure_ascii=False) + "\n").encode("utf-8"))
                if n_read % 100000 == 0:
                    print(f"  {name}: {n_read:,} records read, {len(prec):,} MS/MS kept ({hms(time.monotonic() - t0)})", flush=True)
        prec = np.asarray(prec, dtype=np.float64)
        np.save(tmp / "prec.npy", prec)
        np.save(tmp / "mode.npy", np.asarray(mode, dtype=np.int8))
        np.save(tmp / "offsets.npy", np.asarray(offsets, dtype=np.int64))
        np.save(tmp / "meta_offsets.npy", np.asarray(meta_off, dtype=np.int64))
        np.save(tmp / "order.npy", np.argsort(prec, kind="stable"))
        st = path.stat()
        write_json(tmp / "info.json", {"library": name, "source": str(path), "source_size": st.st_size, "source_mtime": st.st_mtime,
                                       "records_read": n_read, "msms_kept": len(prec), "msms_without_peaks": n_no_peaks,
                                       "peaks": int(offsets[-1]), "built": stamp(), "seconds": round(time.monotonic() - t0, 1),
                                       "rule": "MS2 record (MassBank Spectrum_type MS2, GNPS MSLEVEL 2) with a precursor m/z; ion mode from the "
                                               "ion-mode field, else the sign of the adduct"})
        shutil.rmtree(dst, ignore_errors=True)
        tmp.rename(dst)
        print(f"  {name}: {len(prec):,} MS/MS spectra of {n_read:,} records ({hms(time.monotonic() - t0)})", flush=True)
    no_mzdesk()


class RawLibrary:
    def __init__(self, folder):
        import numpy as np
        self.folder = Path(folder)
        self.prec = np.load(self.folder / "prec.npy", mmap_mode="r")
        self.mode = np.load(self.folder / "mode.npy", mmap_mode="r")
        self.offsets = np.load(self.folder / "offsets.npy", mmap_mode="r")
        self.order = np.load(self.folder / "order.npy", mmap_mode="r")
        self.meta_off = np.load(self.folder / "meta_offsets.npy", mmap_mode="r")
        self.mz = np.memmap(self.folder / "mz.f64", dtype=np.float64, mode="r")
        self.intensity = np.memmap(self.folder / "intensity.f64", dtype=np.float64, mode="r")
        self.prec_sorted = np.asarray(self.prec)[np.asarray(self.order)]
        self._meta = open(self.folder / "meta.jsonl", "rb")

    def window(self, precursor_mz, tol_ppm, sign):
        """Rows of the ion mode `sign` with |query - library| / library <= tol_ppm (the window mzDesk's index uses)."""
        import numpy as np
        t = tol_ppm * 1e-6
        lo = int(np.searchsorted(self.prec_sorted, precursor_mz / (1 + t), side="left"))
        hi = int(np.searchsorted(self.prec_sorted, precursor_mz / (1 - t), side="right"))
        rows = np.asarray(self.order[lo:hi])
        return np.sort(rows[np.asarray(self.mode)[rows] == sign])

    def peaks(self, row):
        import numpy as np
        a, b = int(self.offsets[row]), int(self.offsets[row + 1])
        return np.column_stack((self.mz[a:b], self.intensity[a:b]))

    def meta(self, row):
        self._meta.seek(int(self.meta_off[row]))
        return json.loads(self._meta.readline())


def cmd_direct_library(folder, cache):
    """Direct library search of every recorded call in folder/calls.json -> folder/direct.json. No mzDesk code."""
    import ms_entropy
    import numpy as np
    calls = read_json(folder / "calls.json")
    libs = {}
    out = []
    for c in calls:
        t0 = time.perf_counter()
        q = np.asarray(c["peaks"], dtype=np.float64).reshape(-1, 2)
        sign = 1 if c["positive"] else -1
        kept, n_scored = [], {}
        for name in c["libraries"]:
            if name not in libs:
                libs[name] = RawLibrary(cache / name)
            lib = libs[name]
            rows = lib.window(c["precursor_mz"], c["precursor_tol_ppm"], sign)
            n_scored[name] = int(rows.size)
            for r in rows:
                s = float(ms_entropy.calculate_entropy_similarity(q, lib.peaks(r), ms2_tolerance_in_da=c["tolerance"]))
                if s >= c["threshold"]:
                    kept.append({"library": name, "row": int(r), "score": s, "precursor_mz": float(lib.prec[r]), **lib.meta(r)})
        kept.sort(key=lambda h: -h["score"])
        out.append({"n_scored": n_scored, "kept": kept, "seconds": time.perf_counter() - t0,
                    "function": "ms_entropy.calculate_entropy_similarity", "ms_entropy": getattr(ms_entropy, "__version__", None)})
    write_json(folder / "direct.json", out)
    no_mzdesk()


# ================================================================== BUDDY, direct route
def cmd_direct_buddy(folder):
    """msbuddy on every recorded call in folder/calls.json -> folder/direct.json. No mzDesk code."""
    import numpy as np
    import msbuddy
    from msbuddy import Msbuddy, MsbuddyConfig
    from msbuddy.base import MetaFeature, Spectrum

    def spectrum(peaks):
        if not peaks:
            return None
        return Spectrum(mz_array=np.array([p[0] for p in peaks], dtype=float), int_array=np.array([p[1] for p in peaks], dtype=float))

    out = []
    for c in read_json(folder / "calls.json"):
        t0 = time.perf_counter()
        try:
            if not 1.0 < float(c["precursor_mz"]) < 1500.0:
                raise ValueError("msbuddy annotates 1 < m/z < 1500 only")
            if isinstance(c["ms_instr"], str):
                config = MsbuddyConfig(ms_instr=c["ms_instr"], halogen=c["halogen"])
            else:
                config = MsbuddyConfig(ms_instr=None, ppm=True, ms1_tol=float(c["ms_instr"][0]), ms2_tol=float(c["ms_instr"][1]), halogen=c["halogen"])
            feature = MetaFeature(identifier=0, mz=float(c["precursor_mz"]), charge=1 if c["positive"] else -1, rt=c["rt"], adduct=c["adduct"],
                                  ms1=spectrum(c["ms1"]), ms2=spectrum(c["ms2"]))
            engine = Msbuddy(config)
            engine.add_data([feature])
            engine.annotate_formula()
            meta = engine.data[0]
            cands = [{"rank": i + 1, "formula": str(x.formula), "prob": None if x.estimated_prob is None else float(x.estimated_prob),
                      "normed_prob": None if x.normed_estimated_prob is None else float(x.normed_estimated_prob),
                      "fdr": None if x.estimated_fdr is None else float(x.estimated_fdr),
                      "isotope_sim": None if x.ms1_isotope_similarity is None else float(x.ms1_isotope_similarity)}
                     for i, x in enumerate(meta.candidate_formula_list or [])]
            engine.clear_data()
            out.append({"status": "ok", "adduct_used": str(meta.adduct.string), "candidates": cands, "seconds": time.perf_counter() - t0})
        except Exception as exc:        # noqa: BLE001 - recorded and compared
            out.append({"status": f"error: {type(exc).__name__}: {exc}", "candidates": [], "seconds": time.perf_counter() - t0})
    write_json(folder / "direct.json", {"msbuddy": getattr(msbuddy, "__version__", None), "calls": out})
    no_mzdesk()


# ================================================================== 4a. one case, mzDesk route
class _Abort(BaseException):
    """Raised when step 6 starts: steps 6-8 call no engine tested here."""


def _peaks(peaks):
    return [[float(m), float(i)] for m, i in (peaks or [])]


def run_mzdesk(rec, case_dir):
    """Steps 1-5 of the workflow, Force merge and Run FRIGID on, every engine call and every layer recorded."""
    mfe = _mfe()
    from mzdesk import adducts
    from mzdesk.annotation import formula_step, frigid_step, library_step, msfinder_step, pipeline
    from mzdesk.ms1 import msbuddy_engine as bd
    from mzdesk.ms2 import frigid_engine as fe
    from mzdesk.ms2 import library_search as ls
    from mzdesk.ms2 import msfinder_engine as mf

    calls = {t: [] for t in TOOLS}
    layers, steps, patches = {}, {}, []

    def patch(obj, name, fn):
        patches.append((obj, name, getattr(obj, name)))
        setattr(obj, name, fn)

    real = {"bd": bd.annotate_formula, "ls": ls.search, "mf": mf.run_watched, "fe": fe.run_watched, "formula": formula_step.run_adduct_and_formula,
            "library": library_step.run_library, "msfinder": msfinder_step.run_msfinder_layer, "frigid": frigid_step.run_frigid_layers}

    def w_buddy(precursor_mz, positive, adduct, rt, formula_input, ms2_peaks, ms1_peaks, ms_instr="orbitrap", halogen=False):
        entry = {"precursor_mz": float(precursor_mz), "positive": bool(positive), "adduct": adduct, "rt": rt, "formula_input": formula_input,
                 "ms2": _peaks(ms2_peaks), "ms1": _peaks(ms1_peaks), "ms_instr": ms_instr if isinstance(ms_instr, str) else [float(x) for x in ms_instr],
                 "halogen": bool(halogen), "status": "ok"}
        t = time.perf_counter()
        try:
            res = real["bd"](precursor_mz, positive, adduct, rt, formula_input, ms2_peaks, ms1_peaks, ms_instr, halogen)
            entry["output"] = res
            return res
        except Exception as exc:                                    # noqa: BLE001 - recorded, then handed back to mzDesk unchanged
            entry["status"] = f"error: {type(exc).__name__}: {exc}"
            raise
        finally:
            entry["seconds"] = time.perf_counter() - t
            calls["buddy"].append(entry)

    def w_formula(*a, **k):
        ar, fr = real["formula"](*a, **k)
        layers["formula"] = {"mode": ar.mode, "retained_adducts": list(ar.retained_adducts), "result": fr.to_dict(),
                             "buddy_names": {f.adduct: adducts.engine_name(f.adduct, "buddy") for f in fr.formulas}}
        return ar, fr

    def w_search(peaks, spectrum_type, libraries, metrics, threshold, **k):
        entry = {"peaks": _peaks(peaks), "spectrum_type": spectrum_type, "libraries": list(libraries), "metrics": list(metrics),
                 "threshold": float(threshold), **{key: v for key, v in k.items()}, "status": "ok"}
        t = time.perf_counter()
        try:
            res = real["ls"](peaks, spectrum_type, libraries, metrics, threshold, **k)
            keys = ("rank", "library", "row", "accession", "name", "adduct", "precursor_mz", "ppm", "matches", "scores", "inchikey", "smiles")
            entry["output"] = {"hits": [{key: h.get(key) for key in keys} for h in res["hits"]],
                               **{key: res.get(key) for key in ("n_kept", "n_candidates", "hidden_copies", "rank_metric", "n_metrics_scored",
                                                                "warnings", "notes", "libraries")}}
            return res
        except Exception as exc:                                    # noqa: BLE001
            entry["status"] = f"error: {type(exc).__name__}: {exc}"
            raise
        finally:
            entry["seconds"] = time.perf_counter() - t
            calls["library"].append(entry)

    def w_library(*a, **k):
        layer, summary = real["library"](*a, **k)
        layers["library"] = {"layer": layer.to_dict(), "gate_passed": summary.passed, "best_id": summary.best_id, "thresholds": summary.thresholds}
        return layer, summary

    def w_mf_watch(cmd, cwd, env, *a, **k):
        d = case_dir / "msfinder" / f"{len(calls['msfinder']):02d}"
        in_dir, out_dir, ini = Path(cmd[cmd.index("-i") + 1]), Path(cmd[cmd.index("-o") + 1]), Path(cmd[cmd.index("-m") + 1])
        shutil.copytree(in_dir, d / "input")                        # the query files as MS-FINDER gets them
        shutil.copy2(ini, d / "param.ini")
        timeout = a[0] if a else k.get("timeout", mf.TIMEOUT_S)
        rss_cap = a[1] if len(a) > 1 else k.get("rss_cap_mb")
        entry = {"dir": str(d.relative_to(case_dir)), "timeout": timeout, "rss_cap_mb": rss_cap, "status": "ok"}
        t = time.perf_counter()
        try:
            return real["mf"](cmd, cwd, env, *a, **k)
        except Exception as exc:                                    # noqa: BLE001
            entry["status"] = f"stopped: {type(exc).__name__}: {exc}"
            raise
        finally:
            entry["seconds"] = time.perf_counter() - t
            entry["search_seconds"] = mfe.search_seconds(in_dir)
            side = d / "mzdesk"
            if out_dir.exists():
                shutil.copytree(out_dir, side / "out")
            shutil.copytree(in_dir, side / "in")                    # .fgt and <query>/<FORMULA>.sfd written by MS-FINDER
            calls["msfinder"].append(entry)

    def w_msfinder(*a, **k):
        res, summary = real["msfinder"](*a, **k)
        layers["msfinder"] = {"layer": mfe._jsonable(res), "summary": mfe._jsonable(summary)}
        return res, summary

    def w_fe_watch(cmd, cwd, env, timeout, *a, **k):
        d = case_dir / "frigid" / f"{len(calls['frigid']):02d}"
        d.mkdir(parents=True)
        job, result = Path(cmd[cmd.index("--input") + 1]), Path(cmd[cmd.index("--output") + 1])
        shutil.copy2(job, d / "job.json")
        entry = {"dir": str(d.relative_to(case_dir)), "timeout": timeout, "python": cmd[0], "driver": cmd[2], "driver_sha256": sha256(cmd[2]),
                 "env": {key: env[key] for key in FRIGID_ENV_KEYS if key in env}, "status": "ok"}
        t = time.perf_counter()
        try:
            code, log, mem = real["fe"](cmd, cwd, env, timeout, *a, **k)
            entry["exit_code"], entry["memory"] = code, mem
            return code, log, mem
        except Exception as exc:                                    # noqa: BLE001
            entry["status"] = f"stopped: {type(exc).__name__}: {exc}"
            raise
        finally:
            entry["seconds"] = time.perf_counter() - t
            if result.exists():
                shutil.copy2(result, d / "mzdesk_result.json")
            calls["frigid"].append(entry)

    def w_frigid(*a, **k):
        layer, warns = real["frigid"](*a, **k)
        layers["frigid"] = {"layer": layer.to_dict(), "warnings": warns}
        return layer, warns

    def stop_at_merge(*a, **k):
        raise _Abort()

    def progress(n, status, note, blocks):
        if status != "running":
            steps[n] = {"status": status, "note": note}

    form = mfe.form_of(rec)
    form.update(force_merge=True, run_frigid=True)
    t0 = time.perf_counter()
    for obj, name, fn in ((bd, "annotate_formula", w_buddy), (ls, "search", w_search), (mf, "run_watched", w_mf_watch), (fe, "run_watched", w_fe_watch),
                          (formula_step, "run_adduct_and_formula", w_formula), (library_step, "run_library", w_library),
                          (msfinder_step, "run_msfinder_layer", w_msfinder), (frigid_step, "run_frigid_layers", w_frigid),
                          (pipeline, "merge_candidates", stop_at_merge)):
        patch(obj, name, fn)
    try:
        query = pipeline.parse_inputs(form)
        assert query.force_merge and query.run_frigid and query.libraries, (query.force_merge, query.run_frigid, query.libraries)
        try:
            pipeline.run_query(query, progress)
        except _Abort:
            pass
    finally:
        for obj, name, fn in reversed(patches):
            setattr(obj, name, fn)
    return {"calls": calls, "layers": layers, "steps": steps, "seconds": time.perf_counter() - t0, "form": form,
            "query": {"tol_ppm": query.tol_ppm, "tol_da": query.tol_da, "halogen": query.halogen, "libraries": query.libraries,
                      "lib_metrics": list(query.lib_metrics), "force_merge": query.force_merge, "run_frigid": query.run_frigid},
            "frigid_device_env": os.environ.get("MZDESK_FRIGID_DEVICE")}


# ================================================================== 4b. one case, direct route
def _run_child(args, log, timeout=3600):
    with open(log, "w") as f:
        p = subprocess.run([sys.executable, "-I", str(SCRIPT), *args], stdout=f, stderr=subprocess.STDOUT, env=direct_env(), timeout=timeout)
    if p.returncode != 0:
        raise RuntimeError(f"{' '.join(args[:1])} failed (exit {p.returncode}); see {log}")


def run_frigid_direct(case_dir, call, driver_dir):
    """FRIGID's driver with FRIGID's Python on the recorded job.json, same FRIGID environment and time limit, in a private temp folder."""
    d = case_dir / call["dir"]
    with tempfile.TemporaryDirectory(prefix="eq-frigid-") as tmp:
        tmp = Path(tmp)
        shutil.copy2(d / "job.json", tmp / "job.json")
        driver = driver_dir / "frigid_predict.py"
        cmd = [call["python"], "-I", str(driver), "--input", str(tmp / "job.json"), "--output", str(tmp / "result.json")]
        env = {**direct_env(), **call["env"]}
        t = time.perf_counter()
        with open(tmp / "log.txt", "w") as log:
            p = subprocess.Popen(cmd, cwd=str(tmp), env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                p.wait(timeout=call["timeout"])
                status = "ok"
            except subprocess.TimeoutExpired:
                os.killpg(p.pid, signal.SIGKILL)
                p.wait()
                status = f"stopped: did not finish within {call['timeout']} s"
        seconds = time.perf_counter() - t
        if (tmp / "result.json").exists():
            shutil.copy2(tmp / "result.json", d / "direct_result.json")
        shutil.copy2(tmp / "log.txt", d / "direct_log.txt")
    return {"status": status, "exit_code": p.returncode, "seconds": seconds, "driver_sha256": sha256(driver)}


# ================================================================== 4c. comparison
@lru_cache(maxsize=None)
def compound_key(inchikey, inchi, smiles, name):
    """First 14 InChIKey characters: the record's InChIKey, else the key of its InChI, else of its SMILES (largest fragment); "name:<name>"
    when there is no structure; None when there is nothing."""
    from rdkit import Chem
    from rdkit.Chem import inchi as rd_inchi
    from rdkit.rdBase import BlockLogs
    ik = _na(inchikey)
    if _IK.fullmatch(ik):
        return ik[:14]
    with BlockLogs():
        if _na(inchi).startswith("InChI="):
            try:
                mol = Chem.MolFromInchi(inchi)
                key = rd_inchi.InchiToInchiKey(inchi) if mol is not None and mol.GetNumAtoms() else ""
                if key and _IK.fullmatch(key):
                    return key[:14]
            except Exception:       # noqa: BLE001 - an odd InChI falls through to the SMILES
                pass
        if _na(smiles):
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None and mol.GetNumHeavyAtoms():
                    mol = max(Chem.GetMolFrags(mol, asMols=True), key=lambda m: m.GetNumHeavyAtoms())
                    key = Chem.MolToInchiKey(mol)
                    if key and _IK.fullmatch(key):
                        return key[:14]
            except Exception:       # noqa: BLE001
                pass
    n = _na(name).strip().lower()
    return f"name:{n}" if n else None


def _same(a, b, tol=SCORE_TOL):
    if a is None or b is None:
        return a is None and b is None
    return abs(a - b) <= tol


def compare_buddy(mz, direct):
    calls, out = mz["calls"]["buddy"], direct["calls"]
    runs = []
    for i, (c, d) in enumerate(zip(calls, out)):
        a = [(x["formula"], x["prob"]) for x in (c.get("output") or {}).get("candidates", [])]
        b = [(x["formula"], x["prob"]) for x in d["candidates"]][:BUDDY_TOP]
        ok_a, ok_b = c["status"] == "ok", d["status"] == "ok"
        same = (len(a) == len(b) and all(x[0] == y[0] and _same(x[1], y[1]) for x, y in zip(a, b))) if ok_a and ok_b else (not ok_a and not ok_b)
        diff = max((abs(x[1] - y[1]) for x, y in zip(a, b) if x[1] is not None and y[1] is not None), default=0.0)
        first = None if same else f"{c['adduct']}: mzDesk {a[:3]} vs direct {b[:3]} ({c['status']} / {d['status']})"
        runs.append({"n": i, "adduct": c["adduct"], "identical": same, "formulas_mzdesk": len(a), "formulas_direct": len(d["candidates"]),
                     "max_prob_diff": diff, "status_mzdesk": c["status"], "status_direct": d["status"], "first_difference": first})
    layer = mz["layers"].get("formula") or {"result": {"formulas": []}, "buddy_names": {}}
    missing = []
    formulas = layer["result"]["formulas"]
    for f in formulas:
        name = layer["buddy_names"].get(f["adduct"])
        idx = next((i for i, c in enumerate(calls) if c["adduct"] == name), None)
        found = idx is not None and any(x["formula"] == f["formula"] and _same(x["prob"], f["prob"]) for x in out[idx]["candidates"])
        if not found:
            missing.append({"formula": f["formula"], "adduct": f["adduct"], "prob": f["prob"]})
    return {"runs": runs, "layer": {"candidates": len(formulas), "found": len(formulas) - len(missing), "missing": missing}}


def compare_library(mz, direct):
    runs, layer_missing, n_layer = [], [], 0
    spectra_total = spectra_found = 0
    for i, (c, d) in enumerate(zip(mz["calls"]["library"], direct)):
        metric = c["metrics"][0] if c["metrics"] else None
        hits = (c.get("output") or {}).get("hits", [])
        n_kept = (c.get("output") or {}).get("n_kept", 0)
        a = {}
        for h in hits:
            k = compound_key(h.get("inchikey"), None, h.get("smiles"), h.get("name")) or f"spectrum:{h['library']}:{h['accession']}"
            s = (h.get("scores") or {}).get(metric)
            if s is not None and (k not in a or s > a[k]):
                a[k] = s
        b = {}
        for h in d["kept"]:
            k = compound_key(h.get("inchikey"), h.get("inchi"), h.get("smiles"), h.get("name")) or f"spectrum:{h['library']}:{h['accession']}"
            if k not in b or h["score"] > b[k]:
                b[k] = h["score"]
        truncated = n_kept > len(hits)
        floor = min((h["scores"][metric] for h in hits if (h.get("scores") or {}).get(metric) is not None), default=0.0) if truncated else None
        b_asked = {k: v for k, v in b.items() if floor is None or v >= floor}
        only_a, only_b = sorted(set(a) - set(b)), sorted(set(b_asked) - set(a))
        diffs = {k: abs(a[k] - b[k]) for k in set(a) & set(b)}
        comparable = metric == "SpectralEntropy (weighted)" and c["status"] == "ok"
        same = comparable and not only_a and not only_b and all(v <= LIBRARY_TOL for v in diffs.values())
        exact = same and all(v <= SCORE_TOL for v in diffs.values())
        # per spectrum: every spectrum the mzDesk search returned, in the direct output with the same score
        by_acc = {(h["library"], h["accession"]): h["score"] for h in d["kept"]}
        sp_missing = [{"library": h["library"], "accession": h["accession"], "score": h["scores"].get(metric),
                       "direct": by_acc.get((h["library"], h["accession"]))}
                      for h in hits if not _same(by_acc.get((h["library"], h["accession"])), h["scores"].get(metric))]
        spectra_total += len(hits)
        spectra_found += len(hits) - len(sp_missing)
        first = None
        if not same:
            first = (f"not comparable: metric {metric}, status {c['status']}" if not comparable else
                     f"only mzDesk {only_a[:5]}, only direct {only_b[:5]}, largest best-score difference {max(diffs.values(), default=0):.4f}")
        mz_acc = {(h["library"], h["accession"]) for h in hits}
        runs.append({"n": i, "identical": same, "identical_scores_1e-6": exact, "comparable": comparable, "truncated": truncated,
                     "compounds_mzdesk": len(a), "compounds_direct": len(b), "compounds_direct_asked": len(b_asked),
                     "only_mzdesk": only_a, "only_direct": only_b, "max_best_score_diff": max(diffs.values(), default=0.0),
                     "spectra_mzdesk": len(hits), "spectra_kept_mzdesk": n_kept, "spectra_direct": len(d["kept"]),
                     "spectra_direct_not_returned_by_mzdesk": sum((h["library"], h["accession"]) not in mz_acc for h in d["kept"]),
                     "hidden_copies_mzdesk": (c.get("output") or {}).get("hidden_copies"),
                     "spectra_not_found": sp_missing, "status_mzdesk": c["status"], "seconds_mzdesk": c["seconds"], "seconds_direct": d["seconds"],
                     "first_difference": first})
        # layer: each library candidate's compound has the same best score in the direct output
        hit_of = {(h["library"], h["accession"]): h for h in hits}
        cands = ((mz["layers"].get("library") or {}).get("layer") or {}).get("candidates") or []
        for cand in cands:
            ev = (cand.get("evidence") or {}).get("library") or {}
            h = hit_of.get((ev.get("library"), ev.get("accession")))
            k = (compound_key(h.get("inchikey"), None, h.get("smiles"), h.get("name")) if h else None) or cand.get("ik14")
            n_layer += 1
            if not (k in b and _same(b[k], cand.get("layer_score"))):
                layer_missing.append({"key": k, "name": cand.get("name"), "accession": ev.get("accession"), "library": ev.get("library"),
                                      "score": cand.get("layer_score"), "direct_best": b.get(k)})
    return {"runs": runs, "layer": {"candidates": n_layer, "found": n_layer - len(layer_missing), "missing": layer_missing,
                                    "spectra": spectra_total, "spectra_found": spectra_found}}


def compare_frigid(case_dir, mz, direct):
    runs, raw_direct = [], []
    for i, (c, d) in enumerate(zip(mz["calls"]["frigid"], direct)):
        p = case_dir / c["dir"]
        ra = read_json(p / "mzdesk_result.json") if (p / "mzdesk_result.json").exists() else {}
        rb = read_json(p / "direct_result.json") if (p / "direct_result.json").exists() else {}
        a, b = ra.get("candidates") or [], rb.get("candidates") or []
        raw_direct += b
        stopped = c["status"] != "ok" or d["status"] != "ok"
        same = (bool(ra.get("ok")) == bool(rb.get("ok")) and len(a) == len(b)
                and all(x.get("smiles") == y.get("smiles") and _same(x.get("score"), y.get("score")) for x, y in zip(a, b)))
        first = None
        if not same:
            k = next((j for j, (x, y) in enumerate(zip(a, b)) if x.get("smiles") != y.get("smiles") or not _same(x.get("score"), y.get("score"))), None)
            first = (f"{len(a)} vs {len(b)} candidates" + (f"; first difference at {k + 1}: {a[k].get('smiles')} {a[k].get('score')} vs "
                                                         f"{b[k].get('smiles')} {b[k].get('score')}" if k is not None else "")
                     + f" (ok {ra.get('ok')} / {rb.get('ok')})")
        job = read_json(p / "job.json")
        runs.append({"dir": c["dir"], "identical": same, "time_limited": stopped, "formula": job.get("formula"), "adduct": job.get("adduct"),
                     "candidates_mzdesk": len(a), "candidates_direct": len(b), "device_mzdesk": ra.get("device_used"), "device_direct": rb.get("device_used"),
                     "threads": c["env"].get("OMP_NUM_THREADS"), "status_mzdesk": c["status"], "status_direct": d["status"],
                     "same_driver": c["driver_sha256"] == d["driver_sha256"], "seconds_mzdesk": c["seconds"], "seconds_direct": d["seconds"],
                     "first_difference": first})
    cands = ((mz["layers"].get("frigid") or {}).get("layer") or {}).get("candidates") or []
    missing = []
    for cand in cands:
        found = any(((cand.get("ik14") and (x.get("inchikey_2d") or "")[:14] == cand["ik14"]) or x.get("smiles") == cand.get("smiles"))
                    and _same(x.get("score"), cand.get("layer_score")) for x in raw_direct)
        if not found:
            missing.append({"ik14": cand.get("ik14"), "smiles": cand.get("smiles"), "score": cand.get("layer_score")})
    return {"runs": runs, "layer": {"candidates": len(cands), "found": len(cands) - len(missing), "missing": missing,
                                    "status": ((mz["layers"].get("frigid") or {}).get("layer") or {}).get("status"),
                                    "note": ((mz["layers"].get("frigid") or {}).get("layer") or {}).get("note")}}


def _ini_pairs(path):
    pairs = {}
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if line and not line.startswith("#") and "=" in line:
            k, _, v = line.partition("=")
            pairs[k.strip()] = v.strip()
    return pairs


def burden(case_dir, mz):
    """Operational burden of the direct route for this spectrum (counting rules in the module docstring)."""
    mfe = _mfe()
    c = mz["calls"]
    default_ini = _ini_pairs(mfe.MSFINDER_DIR / "MSFINDER.INI")
    files = params = handoffs = 0
    for call in c["buddy"]:
        params += 4 + bool(call["ms1"]) + (1 if isinstance(call["ms_instr"], str) else 3) + 1
    params += 7 * len(c["library"])
    for call in c["msfinder"]:
        d = case_dir / call["dir"]
        queries = sorted(p for p in (d / "input").iterdir() if p.is_file())
        files += len(queries) + 1
        handoffs += len(queries)
        mine = _ini_pairs(d / "param.ini")
        params += sum(1 for k, v in mine.items() if default_ini.get(k) != v)
        for q in queries:
            head = q.read_text(encoding="utf-8", errors="replace").split("Num Peaks")[0]
            params += sum(1 for ln in head.splitlines() if ":" in ln)
    for call in c["frigid"]:
        files += 1
        handoffs += 1
        params += len(read_json(case_dir / call["dir"] / "job.json"))
    return {"tools": sum(1 for t in TOOLS if c[t]), "tool_runs": sum(len(c[t]) for t in TOOLS), "files_written": files,
            "parameters_set": params, "handoffs": handoffs, "per_tool_runs": {t: len(c[t]) for t in TOOLS}}


def cmd_case(out, n, cache):
    mfe = _mfe()
    version = mfe.require_code_root()
    assert os.environ.get("MZDESK_FRIGID_DEVICE") == "cpu", "MZDESK_FRIGID_DEVICE must be cpu (the run command sets it)"
    from mzdesk.ms2 import msfinder_engine as mf
    assert mf.install_dir().resolve() == mfe.MSFINDER_DIR.resolve(), f"mzDesk uses {mf.install_dir()}, the direct route {mfe.MSFINDER_DIR}"
    assert mfe.thread_fix_ok(mfe.MSFINDER_DIR), f"{mfe.MSFINDER_DIR} does not carry mzDesk's thread fix; run mzDesk's MS-FINDER setup first"
    rec = {r["id"]: r for r in mfe.load_sample(out)}[n]
    case_dir = out / "cases" / f"{n:03d}"
    shutil.rmtree(case_dir, ignore_errors=True)
    case_dir.mkdir(parents=True)
    mz = run_mzdesk(rec, case_dir)
    write_json(case_dir / "mzdesk.json", mz)

    t = time.perf_counter()
    secs = {}
    (case_dir / "buddy").mkdir()
    write_json(case_dir / "buddy/calls.json", [{k: v for k, v in c.items() if k not in ("output", "status", "seconds")} for c in mz["calls"]["buddy"]])
    t1 = time.perf_counter()
    _run_child(["direct-buddy", "--dir", str(case_dir / "buddy")], case_dir / "buddy/direct.log")
    secs["buddy"] = time.perf_counter() - t1
    (case_dir / "library").mkdir()
    write_json(case_dir / "library/calls.json", [{k: v for k, v in c.items() if k not in ("output", "status", "seconds")} for c in mz["calls"]["library"]])
    t1 = time.perf_counter()
    _run_child(["direct-library", "--dir", str(case_dir / "library"), "--cache", str(cache)], case_dir / "library/direct.log")
    secs["library"] = time.perf_counter() - t1
    t1 = time.perf_counter()
    msf_direct = [mfe.run_direct(case_dir, call) for call in mz["calls"]["msfinder"]]
    secs["msfinder"] = time.perf_counter() - t1
    t1 = time.perf_counter()
    driver_dir = out / "direct" / "frigid_driver"
    frigid_direct = [run_frigid_direct(case_dir, call, driver_dir) for call in mz["calls"]["frigid"]]
    secs["frigid"] = time.perf_counter() - t1
    direct_seconds = time.perf_counter() - t

    msf_layer = {"calls": mz["calls"]["msfinder"], "layer": (mz["layers"].get("msfinder") or {}).get("layer")}
    result = {"id": n, "accession": rec["accession"], "name": rec["name"], "inchikey": rec["inchikey"], "adduct": rec["adduct"],
              "positive": rec["positive"], "mzdesk_version": version, "code_root": os.environ["MZDESK_CODE_ROOT"], "script_sha256": sha256(SCRIPT),
              "seconds": {"mzdesk": mz["seconds"], "direct": direct_seconds, "direct_by_tool": secs},
              "steps": mz["steps"], "retained_adducts": (mz["layers"].get("formula") or {}).get("retained_adducts"),
              "buddy": compare_buddy(mz, read_json(case_dir / "buddy/direct.json")),
              "library": compare_library(mz, read_json(case_dir / "library/direct.json")),
              "msfinder": {"runs": [mfe.compare_run(case_dir, call, d) for call, d in zip(mz["calls"]["msfinder"], msf_direct)],
                           "layer": mfe.compare_layer(case_dir, msf_layer)},
              "frigid": compare_frigid(case_dir, mz, frigid_direct),
              "burden": burden(case_dir, mz)}
    write_json(case_dir / "result.json", result)


# ================================================================== 4d. all cases with a progress bar
def case_verdict(r):
    """identical (every call of every tool identical and every layer candidate found) / differs."""
    ok = all(x["identical"] for t in TOOLS for x in r[t]["runs"]) and all(r[t]["layer"]["found"] == r[t]["layer"]["candidates"] for t in TOOLS)
    return "identical" if ok else "differs"


def case_line(r):
    parts = []
    for t in TOOLS:
        runs = r[t]["runs"]
        lim = sum(x.get("time_limited", False) for x in runs)
        parts.append(f"{TOOL_LABELS[t].split()[0]} {sum(x['identical'] for x in runs)}/{len(runs)}{f'({lim} tl)' if lim else ''}"
                     f" L{r[t]['layer']['found']}/{r[t]['layer']['candidates']}")
    return "  ".join(parts)


def write_status(path, text):
    if path:
        try:
            Path(path).write_text(text + "\n")
        except OSError:
            pass


def run_cases(out, code, cache, ids, workers, status_file):
    mfe = _mfe()
    log = out / "run.log"
    (out / "case_logs").mkdir(exist_ok=True)
    todo = [n for n in ids if not (out / "cases" / f"{n:03d}" / "result.json").exists()]
    with open(log, "a") as f:
        f.write(f"START {stamp()} todo={len(todo)} of {len(ids)} workers={workers} code={code} cache={cache}\n")
    if not todo:
        print(f"All {len(ids)} cases already have results.")
        return
    prog = mfe.Progress(len(ids), workers)
    for n in ids:                               # cases finished by an earlier start count towards the bar
        if n not in todo:
            prog.counts[case_verdict(read_json(out / "cases" / f"{n:03d}" / "result.json"))] += 1
    if len(todo) < len(ids):
        print(f"{len(ids) - len(todo)} of {len(ids)} cases already have results; running the other {len(todo)}.", flush=True)

    def board():
        c = prog.counts
        write_status(status_file, f"equivalence {out.name}: {prog.done()}/{len(ids)} cases done (identical {c['identical']}, differs {c['differs']}, "
                                  f"error {c['error']}), running {len(prog.running)}")
    board()
    procs, stop = {}, threading.Event()

    def one(n):
        if stop.is_set():
            return n, "error", "not started", 0.0
        prog.start(n)
        board()
        t = time.monotonic()
        with open(out / "case_logs" / f"{n:03d}.log", "w") as f:
            p = subprocess.Popen([sys.executable, str(SCRIPT), "case", "--out", str(out), "--id", str(n), "--library-cache", str(cache)],
                                 stdout=f, stderr=subprocess.STDOUT, env=child_env(code), start_new_session=True)
            procs[n] = p
            try:
                p.wait(timeout=CASE_TIMEOUT)
            except subprocess.TimeoutExpired:
                os.killpg(p.pid, signal.SIGKILL)
                p.wait()
            procs.pop(n, None)
        res = out / "cases" / f"{n:03d}" / "result.json"
        sec = time.monotonic() - t
        if p.returncode != 0 or not res.exists():
            return n, "error", f"{n:03d}  ERROR (exit {p.returncode}; see case_logs/{n:03d}.log)", sec
        r = read_json(res)
        v = case_verdict(r)
        return n, v, f"{n:03d}  {v:<9} {case_line(r)}  [{r['accession']} {r['adduct']}, {hms(sec)}]", sec

    ticker_stop = threading.Event()

    def ticker():
        while not ticker_stop.wait(1.0):
            prog.draw()

    th = threading.Thread(target=ticker, daemon=True)
    th.start()
    try:
        with cf.ThreadPoolExecutor(max_workers=workers) as pool:
            futs = [pool.submit(one, n) for n in todo]
            for fut in cf.as_completed(futs):
                n, v, text, sec = fut.result()
                prog.finish(n, v, text, sec)
                board()
                with open(log, "a") as f:
                    f.write(f"{stamp()} id={n:03d} status={v} seconds={sec:.1f}\n")
    except KeyboardInterrupt:
        stop.set()
        for p in list(procs.values()):
            try:
                os.killpg(p.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        ticker_stop.set()
        print("\nStopped. Run the same command again to continue with the cases that have no result yet.")
        with open(log, "a") as f:
            f.write(f"INTERRUPTED {stamp()}\n")
        write_status(status_file, f"equivalence {out.name}: interrupted at {prog.done()}/{len(ids)} cases")
        sys.exit(130)
    ticker_stop.set()
    th.join()
    if prog.tty:
        sys.stdout.write("\n")
    with open(log, "a") as f:
        f.write(f"DONE {stamp()} {json.dumps(prog.counts)} wall={time.monotonic() - prog.t0:.1f}\n")


# ================================================================== 5. summary
def _spread(values):
    if not values:
        return {"median": None, "min": None, "max": None, "text": "n/a"}
    med = statistics.median(values)
    return {"median": med, "min": min(values), "max": max(values), "text": f"{med:g} ({min(values)}-{max(values)})"}


def cmd_summary(out, ids=None):
    mfe = _mfe()
    rate = mfe.rate
    results = []
    for path in sorted((out / "cases").glob("*/result.json")):
        r = read_json(path)
        if ids is None or r["id"] in ids:
            mz_path, direct_path = path.parent / "mzdesk.json", path.parent / "library/direct.json"
            if mz_path.exists() and direct_path.exists():                 # redo the library comparison at the current LIBRARY_TOL
                lib = compare_library(read_json(mz_path), read_json(direct_path))
                if lib != r["library"]:
                    r["library"] = lib
                    write_json(path, r)
            results.append(r)
    expected = len(ids) if ids is not None else len(mfe.load_sample(out))
    frozen = read_json(out / "code/FROZEN.json") if (out / "code/FROZEN.json").exists() else {}
    runs = {t: [x for r in results for x in r[t]["runs"]] for t in TOOLS}
    layer = {t: (sum(r[t]["layer"]["found"] for r in results), sum(r[t]["layer"]["candidates"] for r in results)) for t in TOOLS}
    tool_level = {}
    for t in TOOLS:
        xs = runs[t]
        clean = [x for x in xs if not x.get("time_limited")]
        tool_level[t] = {"all": rate(sum(x["identical"] for x in xs), len(xs)),
                         "without_time_limited": rate(sum(x["identical"] for x in clean), len(clean)),
                         "time_limited": rate(sum(x["identical"] for x in xs if x.get("time_limited")), sum(bool(x.get("time_limited")) for x in xs))}
    lib_sp = (sum(r["library"]["layer"]["spectra_found"] for r in results), sum(r["library"]["layer"]["spectra"] for r in results))
    bur = [r["burden"] for r in results]
    tools_count = {}
    for b in bur:
        tools_count[b["tools"]] = tools_count.get(b["tools"], 0) + 1
    s = {
        "question": "Do BUDDY, the library search, MS-FINDER and FRIGID give the same result inside mzDesk (Force merge, every step run) as when run "
                    "directly on the recorded inputs and settings?",
        "mzdesk_version": results[0]["mzdesk_version"] if results else None, "commit": frozen.get("commit"), "scripts": frozen.get("scripts"),
        "cases_expected": expected, "cases_finished": len(results), "cases_missing_or_failed": expected - len(results),
        "cases_identical": rate(sum(case_verdict(r) == "identical" for r in results), len(results)),
        "tool_level": tool_level,
        "layer_level": {t: rate(*layer[t]) for t in TOOLS},
        "library_spectra": rate(*lib_sp),
        "library_identical_scores_1e-6": rate(sum(x["identical_scores_1e-6"] for x in runs["library"]), len(runs["library"])),
        "library_max_best_score_diff": max((x["max_best_score_diff"] for x in runs["library"]), default=None),
        "library_direct_spectra_not_returned_by_mzdesk": sum(x["spectra_direct_not_returned_by_mzdesk"] for x in runs["library"]),
        "library_calls_truncated": sum(x["truncated"] for x in runs["library"]),
        "frigid_devices": sorted({str(x.get("device_mzdesk")) for x in runs["frigid"]} | {str(x.get("device_direct")) for x in runs["frigid"]}),
        "burden": {"tools": tools_count, "tool_runs": _spread([b["tool_runs"] for b in bur]), "files_written": _spread([b["files_written"] for b in bur]),
                   "parameters_set": _spread([b["parameters_set"] for b in bur]), "handoffs": _spread([b["handoffs"] for b in bur])},
        "differences": [{"id": r["id"], "tool": t, "where": x.get("dir", x.get("n")), "time_limited": x.get("time_limited", False),
                         "first_difference": x.get("first_difference")} for r in results for t in TOOLS for x in r[t]["runs"] if not x["identical"]],
        "layer_missing": [{"id": r["id"], "tool": t, **m} for r in results for t in TOOLS for m in r[t]["layer"]["missing"]],
        "seconds_per_case": _spread([round(r["seconds"]["mzdesk"] + r["seconds"]["direct"]) for r in results]),
    }
    write_json(out / "summary.json", s)
    rows = []
    for r in results:
        row = {"id": r["id"], "accession": r["accession"], "adduct": r["adduct"], "mode": "positive" if r["positive"] else "negative",
               "verdict": case_verdict(r)}
        for t in TOOLS:
            row[f"{t}_calls"] = len(r[t]["runs"])
            row[f"{t}_identical"] = sum(x["identical"] for x in r[t]["runs"])
            row[f"{t}_time_limited"] = sum(bool(x.get("time_limited")) for x in r[t]["runs"])
            row[f"{t}_layer_found"] = r[t]["layer"]["found"]
            row[f"{t}_layer_candidates"] = r[t]["layer"]["candidates"]
        row.update({"library_spectra_found": r["library"]["layer"]["spectra_found"], "library_spectra": r["library"]["layer"]["spectra"],
                    **{f"burden_{k}": v for k, v in r["burden"].items() if k != "per_tool_runs"},
                    "seconds_mzdesk": round(r["seconds"]["mzdesk"], 1), "seconds_direct": round(r["seconds"]["direct"], 1),
                    "first_difference": next((f"{t}: {x['first_difference']}" for t in TOOLS for x in r[t]["runs"] if not x["identical"]), "")})
        rows.append(row)
    if rows:
        import csv
        with open(out / "summary.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

    def tool_cell(t):
        tl = tool_level[t]
        if t in ("msfinder", "frigid") and tl["time_limited"]["n"]:
            return (f"{tl['without_time_limited']['k']} of {tl['without_time_limited']['n']} runs without a time limit "
                    f"({tl['all']['k']} of {tl['all']['n']} in all)")
        unit = "runs" if t in ("msfinder", "frigid") else "calls"
        return f"{tl['all']['k']} of {tl['all']['n']} {unit}"

    b = s["burden"]
    lines = [
        "# Equivalence experiment: mzDesk route vs direct route",
        "",
        f"mzDesk {s['mzdesk_version']} (commit {(s['commit'] or '?')[:7]}), Force merge on, Run FRIGID on, FRIGID on the CPU. "
        f"Cases: {len(results)} of {expected} finished; {s['cases_missing_or_failed']} missing or failed; "
        f"{s['cases_identical']['k']} with every call and every layer candidate identical.",
        "",
        "| | Tool level, identical | Layer level, candidates found |",
        "|---|---|---|",
        f"| BUDDY | {tool_cell('buddy')} | {layer['buddy'][0]} of {layer['buddy'][1]} |",
        f"| Library search | {tool_cell('library')} | {layer['library'][0]} of {layer['library'][1]} (spectra: {lib_sp[0]} of {lib_sp[1]}) |",
        f"| MS-FINDER | {tool_cell('msfinder')} | {layer['msfinder'][0]} of {layer['msfinder'][1]} |",
        f"| FRIGID | {tool_cell('frigid')} | {layer['frigid'][0]} of {layer['frigid'][1]} |",
        "",
        f"Library search: compounds by the first 14 InChIKey characters, best score per compound within {LIBRARY_TOL}; "
        f"{s['library_identical_scores_1e-6']['k']} of {s['library_identical_scores_1e-6']['n']} calls also within 1e-6 "
        f"(largest best-score difference {s['library_max_best_score_diff'] if s['library_max_best_score_diff'] is not None else 'n/a'}); "
        f"the direct search kept {s['library_direct_spectra_not_returned_by_mzdesk']} spectra that mzDesk did not return; "
        f"{s['library_calls_truncated']} calls had more than mzDesk's 200 returned spectra.",
        f"FRIGID devices: {', '.join(s['frigid_devices']) or 'n/a'}.",
        "",
        f"Operational burden per spectrum (direct route; this script's counting rules): tools {', '.join(f'{k}: {v} spectra' for k, v in sorted(b['tools'].items()))}; "
        f"tool runs {b['tool_runs']['text']}; files written {b['files_written']['text']}; parameters set {b['parameters_set']['text']}; "
        f"hand-offs {b['handoffs']['text']} (median, range).",
        f"Wall time per case (both routes): {s['seconds_per_case']['text']} s.",
        "",
    ]
    if s["differences"]:
        lines += ["## Calls that differ", ""] + [
            f"- case {d['id']:03d} {TOOL_LABELS[d['tool']]} {d['where']}{' (time-limited)' if d['time_limited'] else ''}: {d['first_difference']}"
            for d in s["differences"]] + [""]
    if s["layer_missing"]:
        lines += ["## Layer candidates not found", ""] + [
            f"- case {m['id']:03d} {TOOL_LABELS[m['tool']]}: " + ", ".join(f"{k} {v}" for k, v in m.items() if k not in ("id", "tool"))
            for m in s["layer_missing"]] + [""]
    (out / "report.md").write_text("\n".join(lines) + "\n")
    print("\n" + "\n".join(lines))
    print(f"Files: {out / 'report.md'}, {out / 'summary.csv'}, {out / 'summary.json'}")
    return s


# ================================================================== main
def cmd_run(out, commit, workers, ids, limit, cache, status_file):
    mfe = _mfe()
    out.mkdir(parents=True, exist_ok=True)
    print(f"Equivalence experiment (four tools) -> {out}", flush=True)
    write_status(status_file, f"equivalence {out.name}: freezing code")
    code = freeze_code(out, commit)
    frozen = read_json(code / "FROZEN.json") if (code / "FROZEN.json").exists() else {}
    print(f"code: {code}  {('commit ' + frozen['commit'][:7] + ' ' + frozen.get('subject', '')) if frozen else ''}", flush=True)
    if not mfe.thread_fix_ok(mfe.MSFINDER_DIR):
        sys.exit(f"{mfe.MSFINDER_DIR} does not carry mzDesk's thread fix; run mzDesk's MS-FINDER setup first.")
    copy_frigid_driver(out, code)
    if not all((cache / name / "info.json").exists() for name in RAW_LIBRARIES):
        write_status(status_file, f"equivalence {out.name}: reading the raw MassBank and GNPS files")
        print(f"library: reading the raw MassBank and GNPS files into {cache} (about 10 minutes, once) ...", flush=True)
        subprocess.run([sys.executable, "-I", str(SCRIPT), "prepare-library", "--cache", str(cache)], env=direct_env(), check=True)
    if not (out / "inputs/sample.json").exists():
        write_status(status_file, f"equivalence {out.name}: drawing the sample")
        print("sample: reading MassBank (a few minutes) ...", flush=True)
        subprocess.run([sys.executable, str(SCRIPT), "sample", "--out", str(out)], env=child_env(code), check=True)
    sample = mfe.load_sample(out)
    if REFERENCE_SAMPLE.exists():
        same = [r["accession"] for r in sample] == [r["accession"] for r in read_json(REFERENCE_SAMPLE)["records"]]
        print(f"sample: {'the same' if same else 'NOT the same'} {len(sample)} spectra as {REFERENCE_SAMPLE.relative_to(REPO)}", flush=True)
    all_ids = [r["id"] for r in sample]
    chosen = [n for n in all_ids if n in set(ids)] if ids else all_ids
    if limit:
        chosen = chosen[:limit]
    run_cases(out, code, cache, chosen, workers, status_file)
    s = cmd_summary(out, set(chosen) if (ids or limit) else None)
    write_status(status_file, f"equivalence {out.name}: finished {s['cases_finished']}/{s['cases_expected']} cases, "
                              f"{s['cases_identical']['k']} fully identical; report in {out / 'report.md'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["run", "summary", "sample", "case", "prepare-library", "direct-buddy", "direct-library"])
    ap.add_argument("--out", default=None, help="output folder (default: validation/equivalence/run_<YYYYMMDD>)")
    ap.add_argument("--commit", default="HEAD", help="mzDesk commit under test (frozen into OUT/code)")
    ap.add_argument("--workers", type=int, default=2, help="cases at once (each runs one engine at a time; MS-FINDER ~2.3-4 GB, FRIGID ~5 GB)")
    ap.add_argument("--ids", type=int, nargs="*", help="only these case ids (1-50 positive, 51-100 negative)")
    ap.add_argument("--limit", type=int, help="only the first N cases")
    ap.add_argument("--library-cache", default=None, help="raw-library cache of the direct route (default: OUT/direct_library)")
    ap.add_argument("--status", default=None, help="file whose one line is kept current (progress board)")
    ap.add_argument("--id", type=int, help=argparse.SUPPRESS)
    ap.add_argument("--dir", default=None, help=argparse.SUPPRESS)
    ap.add_argument("--cache", default=None, help=argparse.SUPPRESS)
    args = ap.parse_args()
    if args.command == "prepare-library":
        return cmd_prepare_library(Path(args.cache))
    if args.command == "direct-buddy":
        return cmd_direct_buddy(Path(args.dir))
    if args.command == "direct-library":
        return cmd_direct_library(Path(args.dir), Path(args.cache))
    out = Path(args.out).resolve() if args.out else HERE / f"run_{datetime.now():%Y%m%d}"
    cache = Path(args.library_cache).resolve() if args.library_cache else out / "direct_library"
    if args.command == "run":
        cmd_run(out, args.commit, args.workers, args.ids, args.limit, cache, args.status)
    elif args.command == "summary":
        cmd_summary(out)
    elif args.command == "sample":
        _mfe().cmd_sample(out)
    else:
        cmd_case(out, args.id, cache)


if __name__ == "__main__":
    main()
