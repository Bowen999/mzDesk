#!/usr/bin/env python3
"""MS-FINDER result consistency: does MS-FINDER give the same answer inside mzDesk as when it is run on its own?

One command does everything and can be left alone (a progress bar shows where it is; the rates are printed and written at the end):

  python3 validation/equivalence/msfinder_equivalence.py run [--out DIR] [--workers 2] [--ids 1 2 51] [--limit N]
  python3 validation/equivalence/msfinder_equivalence.py summary --out DIR        # rates again from the finished cases

`run` is resumable: started again on the same --out it skips the cases that have a result.json. Ctrl-C stops the running cases.

Steps of `run`
  1. code     the mzDesk code under test is frozen in OUT/code (git archive of HEAD's mzdesk/; MZDESK_CODE_ROOT overrides it);
              OUT/code/FROZEN.json names the commit and says whether mzdesk/ had uncommitted changes.
  2. sample   100 MassBank 2026.03 spectra (seed 20261102, 50 positive + 50 negative compounds; the same rules and draw as the
              v0.5.3 four-layer run of 2 Nov 2026) -> OUT/inputs/sample.json.
  3. cases    one process per spectrum (--workers at once):
              mzDesk side  the ESI MS/MS Annotation workflow with the default form (precursor m/z, ion mode, MS2 peaks; adduct Auto) up to
                           step 4 (MS-FINDER). The library search is switched off (it does not feed MS-FINDER) and the run stops after
                           step 4. Every MsfinderConsoleApp run mzDesk starts is recorded: its query files, parameter file, result files
                           and the .fgt / .sfd files MS-FINDER writes next to the queries (what mzDesk reads its structures from).
              direct side  `mono MsfinderConsoleApp.exe predict -i in -o out -m param.ini` on the same query and parameter files, in a
                           private copy of the same install (~/tools/msfinder/MSFINDER_v3.73, which carries mzDesk's thread fix; the
                           script refuses another state), with the same time limit, without any mzDesk code.
  4. summary  OUT/summary.csv (one row per case), OUT/summary.json and OUT/report.md, and the rates printed.

What is compared (per MsfinderConsoleApp run, then per case)
  run identical   the formula table (Formula result-*.txt), the structure table (Structure result-*.txt) and every structure file
                  (<query>/<FORMULA>.sfd: candidate InChIKey, SMILES and Total score) are the same on both sides; scores within 1e-6;
                  the "File path" column (a temp folder) is left out; rows are sorted by query, score, then their keys, because MS-FINDER
                  writes equal scores in thread-finish order. "byte identical" additionally needs every .sfd structure block (subscores and
                  fragment tables included) to be the same text.
  layer           every candidate mzDesk shows in its MS-FINDER layer (InChIKey block + Total score) is found in the direct output.
  time limited    a run that mzDesk or the direct run stopped (whole-run time or memory limit), or in which a formula's structure search
                  reached MS-FINDER's StructurePredictionTimeOut (search wall time >= 27 s of 30 s): its candidate list depends on the
                  machine's speed at that moment, so the rates are given with and without these runs.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import copy
import csv
import dataclasses
import json
import math
import os
import random
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
if os.environ.get("MZDESK_CODE_ROOT"):
    sys.path.insert(0, str(Path(os.environ["MZDESK_CODE_ROOT"]).resolve()))

MASSBANK_MSP = Path.home() / ".mzdesk/data/libraries/raw/MassBank_NISTformat.msp"
MSFINDER_DIR = Path(os.environ.get("MZDESK_MSFINDER_DIR") or Path.home() / "tools/msfinder/MSFINDER_v3.73")
SEED = 20261102
N_PER_MODE = 50
SCORE_TOL = 1e-6
CASE_TIMEOUT = 3600
SEARCH_LIMIT_S = 27.0               # 0.9 x StructurePredictionTimeOut (0.5 min), as msfinder_step marks "time limit"
HIGHRES = re.compile(r"TOF|FT|Orbitrap", re.I)
ALLOWED_ELEMENTS = {"C", "H", "N", "O", "P", "S", "F", "Cl", "Br", "I"}


def stamp():
    return datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


def write_json(path, payload):
    path = Path(path)
    tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
    tmp.write_text(json.dumps(payload, indent=1, ensure_ascii=False, default=str) + "\n")
    tmp.rename(path)


def hms(seconds):
    seconds = int(max(0, seconds))
    return f"{seconds // 3600}h{seconds % 3600 // 60:02d}m" if seconds >= 3600 else f"{seconds // 60}m{seconds % 60:02d}s"


# ================================================================== 1. frozen code
def freeze_code(out):
    """OUT/code holds the mzdesk package of HEAD (or MZDESK_CODE_ROOT is used as it is). Returns the code root."""
    if os.environ.get("MZDESK_CODE_ROOT"):
        return Path(os.environ["MZDESK_CODE_ROOT"]).resolve()
    code = out / "code"
    if (code / "FROZEN.json").exists():
        return code
    git = lambda *a: subprocess.run(["git", "-C", str(REPO), *a], capture_output=True, text=True, check=True).stdout.strip()   # noqa: E731
    commit = git("rev-parse", "HEAD")
    dirty = git("status", "--porcelain", "--", "mzdesk")
    if dirty:
        print(f"WARNING: mzdesk/ has uncommitted changes; the run tests commit {commit[:7]} without them:\n{dirty}", flush=True)
    tmp = out / "code.partial"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    archive = subprocess.run(["git", "-C", str(REPO), "archive", "--format=tar", commit, "mzdesk", "pyproject.toml"], capture_output=True, check=True)
    subprocess.run(["tar", "-x", "-C", str(tmp)], input=archive.stdout, check=True)
    write_json(tmp / "FROZEN.json", {"commit": commit, "subject": git("log", "-1", "--format=%s", commit), "frozen_at": stamp(),
                                     "uncommitted_changes_left_out": dirty or None})
    tmp.rename(code)
    return code


def child_env(code):
    env = {**os.environ, "MZDESK_CODE_ROOT": str(code), "PYTHONUNBUFFERED": "1"}
    env.pop("PYTHONPATH", None)
    return env


def require_code_root():
    code = os.environ.get("MZDESK_CODE_ROOT")
    assert code, "MZDESK_CODE_ROOT is not set (the `run` command sets it for its children)."
    import mzdesk
    assert Path(mzdesk.__file__).resolve().is_relative_to(Path(code).resolve()), f"mzdesk imported from {mzdesk.__file__}, not {code}"
    return mzdesk.__version__


# ================================================================== 2. sample (same rules and draw as the v0.5.3 four-layer run)
def parse_msp(path):
    """MassBank NIST-format .msp -> dicts (first value of each key) with "peaks": [(mz, intensity)]."""
    rec, peaks, in_peaks = {}, [], False
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if not line.strip():
                if rec:
                    rec["peaks"] = peaks
                    yield rec
                rec, peaks, in_peaks = {}, [], False
                continue
            if in_peaks:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        peaks.append((float(parts[0]), float(parts[1])))
                    except ValueError:
                        pass
                continue
            key, _, value = line.partition(":")
            key = key.strip()
            if key.lower() == "num peaks":
                in_peaks = True
                continue
            rec.setdefault(key, value.strip())
    if rec:
        rec["peaks"] = peaks
        yield rec


def eligible(r, adducts_by_mode):
    """A MassBank record usable as a query, or None."""
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    from rdkit.rdBase import BlockLogs
    from mzdesk.adducts import canonical
    from mzdesk.ms1 import msbuddy_engine as bd
    if r.get("Spectrum_type") != "MS2":
        return None
    mode = r.get("Ion_mode", "").upper()
    if mode not in ("POSITIVE", "NEGATIVE"):
        return None
    positive = mode == "POSITIVE"
    itype = r.get("Instrument_type", "")
    if "ESI" not in itype.upper() or not HIGHRES.search(itype):
        return None
    adduct = canonical(r.get("Precursor_type", "")) or ""
    if adduct not in adducts_by_mode[positive]:
        return None
    ik, smi = r.get("InChIKey", ""), r.get("SMILES", "")
    if not re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", ik) or not smi:
        return None
    try:
        pmz = float(r.get("PrecursorMZ", ""))
    except ValueError:
        return None
    peaks = [(m, i) for m, i in r["peaks"] if m > 0 and i > 0]
    if not 6 <= len(peaks) <= 5000:
        return None
    _ = BlockLogs()
    mol = Chem.MolFromSmiles(smi)
    if mol is None or len(Chem.GetMolFrags(mol)) != 1 or Chem.GetFormalCharge(mol) != 0:
        return None
    if any(a.GetSymbol() not in ALLOWED_ELEMENTS or a.GetNumRadicalElectrons() or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    if not 100 <= Descriptors.ExactMolWt(mol) <= 800:
        return None
    if Chem.MolToInchiKey(mol)[:14] != ik[:14]:
        return None
    formula = rdMolDescriptors.CalcMolFormula(mol)
    theo = bd.ion_mz(formula, adduct)
    if theo is None or abs((pmz - theo) / theo * 1e6) > 10:
        return None
    return {"accession": r.get("DB#", ""), "name": r.get("Name", ""), "inchikey": ik, "smiles": smi, "formula": formula,
            "adduct": adduct, "precursor_mz": pmz, "positive": positive, "instrument_type": itype, "instrument": r.get("Instrument", ""),
            "collision_energy": r.get("Collision_energy", ""), "precursor_ppm": (pmz - theo) / theo * 1e6, "peaks": peaks}


def cmd_sample(out):
    require_code_root()
    from mzdesk import spec
    from mzdesk.adducts import canonical
    adducts_by_mode = {True: {canonical(a) for a in spec.ADDUCTS_POS if canonical(a)},
                       False: {canonical(a) for a in spec.ADDUCTS_NEG if canonical(a)}}
    by_mode = {True: {}, False: {}}
    n_read = 0
    for r in parse_msp(MASSBANK_MSP):
        n_read += 1
        e = eligible(r, adducts_by_mode)
        if e:
            by_mode[e["positive"]].setdefault(e["inchikey"][:14], []).append(e)
    rng = random.Random(SEED)
    chosen = []
    for positive in (True, False):
        pool = sorted(k for k in by_mode[positive] if k not in {c["inchikey"][:14] for c in chosen})
        for k in rng.sample(pool, N_PER_MODE):
            recs = sorted(by_mode[positive][k], key=lambda e: e["accession"])
            chosen.append(rng.choice(recs))
    for i, c in enumerate(chosen, 1):
        c["id"] = i
    meta = {"source": str(MASSBANK_MSP), "massbank_version": "2026.03", "seed": SEED, "records_read": n_read,
            "eligible_compounds": {"positive": len(by_mode[True]), "negative": len(by_mode[False])},
            "rules": "MS2, ESI, TOF/FT/Orbitrap instrument, adduct in mzDesk's Auto list of the polarity, InChIKey + SMILES that agree, "
                     "one neutral component of C H N O P S F Cl Br I, 100-800 Da, 6-5000 peaks, precursor within 10 ppm of formula + adduct; "
                     "one random spectrum per compound; 50 positive compounds, then 50 negative compounds not already chosen"}
    (out / "inputs").mkdir(parents=True, exist_ok=True)
    write_json(out / "inputs/sample.json", {"meta": meta, "records": chosen})
    print(f"sample: {len(chosen)} spectra ({meta['eligible_compounds']} eligible compounds)", flush=True)


def load_sample(out):
    return json.loads((out / "inputs/sample.json").read_text())["records"]


# ================================================================== 3a. one case, mzDesk side
class _Abort(BaseException):
    """Raised after step 4 to end the workflow run (steps 5-8 are not part of this test)."""


def _jsonable(x):
    if dataclasses.is_dataclass(x) and not isinstance(x, type):
        return x.to_dict() if hasattr(x, "to_dict") else {k: _jsonable(v) for k, v in dataclasses.asdict(x).items()}
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, set)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (str, int, float, bool)) or x is None:
        return x
    return repr(x)


def search_seconds(in_dir):
    """{"<query>/<FORMULA>": wall seconds of that structure search} from the .sfd files' birth and modification times (macOS); a 0-byte .sfd
    (a search that was still running when the run was stopped) is None. Must be read before the files are copied (a copy has a new birth time)."""
    out = {}
    for p in sorted(Path(in_dir).glob("*/*.sfd")):
        st = p.stat()
        born = getattr(st, "st_birthtime", None)
        out[f"{p.parent.name}/{p.stem}"] = None if born is None or st.st_size == 0 else max(0.0, st.st_mtime - born)
    return out


def form_of(rec):
    """The workflow form as the page sends it: every field at its default (spec); the user fills precursor m/z, ion mode and MS2 peaks."""
    from mzdesk import spec
    pipe = next(p for p in spec.PIPELINES if p["id"] == "esi_ms2_annotation")
    form = {f["id"]: copy.deepcopy(f["value"]) for f in pipe["inputs"] if f.get("value") is not None}
    peaks = "\n".join(f"{m:.6f} {i:.6g}" for m, i in rec["peaks"])
    form.update({"precursor": {"precursor_mz": rec["precursor_mz"], "ion_mode": "Positive" if rec["positive"] else "Negative", "adduct": "Auto"},
                 "ms2": {"source": "Paste", "value": peaks}, "ms1": ""})
    return form


def run_mzdesk(rec, case_dir):
    """Steps 1-4 of the workflow with every MsfinderConsoleApp run recorded under case_dir/runs/NN/mzdesk."""
    from mzdesk.annotation import msfinder_step, pipeline
    from mzdesk.ms2 import msfinder_engine as mf
    calls, layer = [], {}
    real_watch, real_layer = mf.run_watched, msfinder_step.run_msfinder_layer

    def w_watch(cmd, cwd, env, *a, **k):
        d = case_dir / "runs" / f"{len(calls):02d}"
        in_dir, out_dir, ini = Path(cmd[cmd.index("-i") + 1]), Path(cmd[cmd.index("-o") + 1]), Path(cmd[cmd.index("-m") + 1])
        shutil.copytree(in_dir, d / "input")                        # the query files as MS-FINDER gets them
        shutil.copy2(ini, d / "param.ini")
        timeout = a[0] if a else k.get("timeout", mf.TIMEOUT_S)
        rss_cap = a[1] if len(a) > 1 else k.get("rss_cap_mb")
        entry = {"dir": str(d.relative_to(case_dir)), "timeout": timeout, "rss_cap_mb": rss_cap, "status": "ok"}
        t = time.perf_counter()
        try:
            return real_watch(cmd, cwd, env, *a, **k)
        except Exception as exc:                                    # noqa: BLE001 - recorded, then handed back to mzDesk unchanged
            entry["status"] = f"stopped: {type(exc).__name__}: {exc}"
            raise
        finally:
            entry["seconds"] = time.perf_counter() - t
            entry["search_seconds"] = search_seconds(in_dir)
            side = d / "mzdesk"
            if out_dir.exists():
                shutil.copytree(out_dir, side / "out")
            shutil.copytree(in_dir, side / "in")                    # .fgt and <query>/<FORMULA>.sfd written by MS-FINDER
            calls.append(entry)

    def w_layer(*a, **k):
        res, summary = real_layer(*a, **k)
        layer["layer"], layer["summary"] = _jsonable(res), _jsonable(summary)
        raise _Abort()

    form = form_of(rec)
    t0 = time.perf_counter()
    mf.run_watched, msfinder_step.run_msfinder_layer = w_watch, w_layer
    try:
        query = pipeline.parse_inputs(form)
        query.libraries = []                    # step 3 skipped: the library search does not feed MS-FINDER, and the gate can then not skip step 4
        try:
            pipeline.run_query(query)
        except _Abort:
            pass
    finally:
        mf.run_watched, msfinder_step.run_msfinder_layer = real_watch, real_layer
    return {"calls": calls, "layer": layer.get("layer"), "summary": layer.get("summary"), "seconds": time.perf_counter() - t0,
            "query": {"tol_ppm": query.tol_ppm, "tol_da": query.tol_da, "halogen": query.halogen}}


# ================================================================== 3b. one case, direct side (no mzDesk code)
def thread_fix_ok(install):
    """The install's StructureFinder.dll is the one named by the marker mzDesk's setup wrote (the fixed build)."""
    import hashlib
    dll, marker = install / "StructureFinder.dll", install / "StructureFinder.dll.mzdesk-threadfix.json"
    if not dll.is_file() or not marker.is_file():
        return False
    digest = hashlib.sha256(dll.read_bytes()).hexdigest()
    return json.loads(marker.read_text()).get("fixed_sha256") == digest


def clone_install(dst):
    """Copy-on-write clone of the install (APFS: instant), else a plain copy."""
    if subprocess.run(["cp", "-cR", str(MSFINDER_DIR), str(dst)], capture_output=True).returncode != 0:
        shutil.rmtree(dst, ignore_errors=True)
        shutil.copytree(MSFINDER_DIR, dst, symlinks=True)
    return dst / "MsfinderConsoleApp.exe"


def run_direct(case_dir, call):
    """MsfinderConsoleApp.exe on the recorded query and parameter files, same time limit, in a private copy of the install."""
    d = case_dir / call["dir"]
    mono = os.environ.get("MZDESK_MONO") or shutil.which("mono") or "/opt/homebrew/bin/mono"
    with tempfile.TemporaryDirectory(prefix="eq-msfinder-") as tmp:
        tmp = Path(tmp)
        exe = clone_install(tmp / "install")
        shutil.copytree(d / "input", tmp / "in")
        (tmp / "out").mkdir()
        shutil.copy2(d / "param.ini", tmp / "param.ini")
        cmd = [mono, str(exe), "predict", "-i", str(tmp / "in"), "-o", str(tmp / "out"), "-m", str(tmp / "param.ini")]
        env = {**os.environ, "LC_ALL": "en_US.UTF-8", "LANG": "en_US.UTF-8"}
        t = time.perf_counter()
        with open(tmp / "log.txt", "w") as log:
            p = subprocess.Popen(cmd, cwd=str(exe.parent), env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                p.wait(timeout=call["timeout"])
                status = "ok"                                       # exit code 1 is MS-FINDER's success code: judged by the files
            except subprocess.TimeoutExpired:
                os.killpg(p.pid, signal.SIGKILL)
                p.wait()
                status = f"stopped: did not finish within {call['timeout']} s"
        seconds = time.perf_counter() - t
        searches = search_seconds(tmp / "in")
        side = d / "direct"
        shutil.copytree(tmp / "out", side / "out")
        shutil.copytree(tmp / "in", side / "in")
        shutil.copy2(tmp / "log.txt", side / "log.txt")
    return {"status": status, "exit_code": p.returncode, "seconds": seconds, "search_seconds": searches}


# ================================================================== 3c. comparison
def _rows(path):
    lines = [ln.rstrip("\r\n") for ln in Path(path).read_text(encoding="utf-8", errors="replace").splitlines() if ln.strip()]
    if not lines:
        return []
    head = lines[0].split("\t")
    return [dict(zip(head, ln.split("\t") + [""] * (len(head) - len(ln.split("\t"))))) for ln in lines[1:]]


def _score(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def _sorted(rows, score_key):
    """Query, best score first, then the rest of the row (ties are written in thread-finish order)."""
    def key(r):
        s = _score(r.get(score_key))
        return (r.get("File name", ""), float("inf") if s is None else -s, json.dumps(r, sort_keys=True))
    return sorted(rows, key=key)


def table(side_dir, pattern, score_key):
    """Rows of the result table(s), "File path" left out (a temp folder), sorted."""
    rows = [{k: v for k, v in r.items() if k != "File path"} for f in sorted((side_dir / "out").glob(pattern)) for r in _rows(f)]
    return _sorted(rows, score_key)


def sfd_blocks(path):
    """Structure blocks of one .sfd file: [{"inchikey", "smiles", "score", "text"}] (text = the whole block, fragments included)."""
    blocks, cur = [], None
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("NAME: "):
            cur = {"inchikey": "", "smiles": "", "score": None, "lines": []}
            blocks.append(cur)
        if cur is None:
            continue
        cur["lines"].append(line)
        key, _, value = line.partition(": ")
        if key == "INCHIKEY":
            cur["inchikey"] = value.strip()
        elif key == "SMILES":
            cur["smiles"] = value.strip()
        elif key == "TotalScore":
            cur["score"] = _score(value)
    out = [{"inchikey": b["inchikey"], "smiles": b["smiles"], "score": b["score"], "text": "\n".join(b["lines"]).rstrip()} for b in blocks]
    return sorted(out, key=lambda b: (float("inf") if b["score"] is None else -b["score"], b["inchikey"], b["smiles"], b["text"]))


def sfd_files(side_dir):
    """{"<query>/<FORMULA>": blocks} of the finished searches (a 0-byte file is a search stopped before it ended)."""
    return {f"{p.parent.name}/{p.stem}": sfd_blocks(p) for p in sorted((side_dir / "in").glob("*/*.sfd")) if p.stat().st_size > 0}


def _same_score(a, b):
    if a is None or b is None:
        return a is None and b is None
    return abs(a - b) <= SCORE_TOL


def same_rows(a, b, score_key):
    if len(a) != len(b):
        return False
    for x, y in zip(a, b):
        if {k: v for k, v in x.items() if k != score_key} != {k: v for k, v in y.items() if k != score_key}:
            return False
        sx, sy = _score(x.get(score_key)), _score(y.get(score_key))
        if (sx is None or sy is None) and x.get(score_key) != y.get(score_key):
            return False
        if not _same_score(sx, sy):
            return False
    return True


def compare_run(case_dir, call, direct):
    """One MsfinderConsoleApp run, mzDesk's vs the direct one -> dict of verdicts and the first difference found."""
    d = case_dir / call["dir"]
    mz, dr = d / "mzdesk", d / "direct"
    fa, fb = table(mz, "Formula result-*.txt", "Formula score"), table(dr, "Formula result-*.txt", "Formula score")
    ta, tb = table(mz, "Structure result-*.txt", "Total score"), table(dr, "Structure result-*.txt", "Total score")
    sa, sb = sfd_files(mz), sfd_files(dr)
    formulas_same, table_same = same_rows(fa, fb, "Formula score"), same_rows(ta, tb, "Total score")
    sfd_same, sfd_bytes, first = set(sa) == set(sb), set(sa) == set(sb), None
    if not sfd_same:
        first = f"structure files differ: only mzDesk {sorted(set(sa) - set(sb))}, only direct {sorted(set(sb) - set(sa))}"
    for key in sorted(set(sa) & set(sb)):
        a, b = sa[key], sb[key]
        same = len(a) == len(b) and all(x["inchikey"] == y["inchikey"] and x["smiles"] == y["smiles"] and _same_score(x["score"], y["score"])
                                        for x, y in zip(a, b))
        if not same:
            sfd_same = False
            if first is None:
                first = (f"{key}: {len(a)} vs {len(b)} structures; best mzDesk "
                         f"{[(x['inchikey'][:14], x['score']) for x in a[:3]]} vs direct {[(y['inchikey'][:14], y['score']) for y in b[:3]]}")
        if [x["text"] for x in a] != [y["text"] for y in b]:
            sfd_bytes = False
    if first is None and not formulas_same:
        first = f"formula table: {len(fa)} vs {len(fb)} rows"
    if first is None and not table_same:
        first = f"structure table: {len(ta)} vs {len(tb)} rows"
    secs = [v for v in list(call["search_seconds"].values()) + list(direct["search_seconds"].values()) if v is not None]
    stopped = call["status"] != "ok" or direct["status"] != "ok"
    at_limit = sorted({k for side in (call["search_seconds"], direct["search_seconds"]) for k, v in side.items() if v is not None and v >= SEARCH_LIMIT_S})
    identical = formulas_same and table_same and sfd_same
    return {"dir": call["dir"], "identical": identical, "byte_identical": identical and sfd_bytes, "formula_table_identical": formulas_same,
            "structure_table_identical": table_same, "structure_files_identical": sfd_same, "structure_files": len(sa),
            "structures_mzdesk": sum(map(len, sa.values())), "structures_direct": sum(map(len, sb.values())),
            "time_limited": stopped or bool(at_limit), "stopped": stopped, "searches_at_limit": at_limit,
            "longest_search_s": max(secs) if secs else None, "status_mzdesk": call["status"], "status_direct": direct["status"],
            "seconds_mzdesk": call["seconds"], "seconds_direct": direct["seconds"], "first_difference": first}


def compare_layer(case_dir, mz):
    """Every candidate of mzDesk's MS-FINDER layer (InChIKey block or SMILES, Total score) must be in a direct .sfd with the same score."""
    direct = [(b["inchikey"][:14] if b["inchikey"] else None, b["smiles"], b["score"])
              for call in mz["calls"] for blocks in sfd_files(case_dir / call["dir"] / "direct").values() for b in blocks]
    cands = (mz["layer"] or {}).get("candidates") or []
    missing = []
    for c in cands:
        found = any((c.get("ik14") and k == c["ik14"] or not c.get("ik14") and s == c.get("smiles")) and _same_score(sc, c.get("layer_score"))
                    for k, s, sc in direct)
        if not found:
            missing.append({"ik14": c.get("ik14"), "formula": c.get("formula"), "score": c.get("layer_score")})
    return {"candidates": len(cands), "found": len(cands) - len(missing), "missing": missing,
            "status": (mz["layer"] or {}).get("status"), "note": (mz["layer"] or {}).get("note")}


def cmd_case(out, n):
    version = require_code_root()
    from mzdesk.ms2 import msfinder_engine as mf
    assert mf.install_dir().resolve() == MSFINDER_DIR.resolve(), f"mzDesk uses {mf.install_dir()}, the direct side {MSFINDER_DIR}"
    assert thread_fix_ok(MSFINDER_DIR), f"{MSFINDER_DIR} does not carry mzDesk's thread fix; run mzDesk's MS-FINDER setup first"
    rec = {r["id"]: r for r in load_sample(out)}[n]
    case_dir = out / "cases" / f"{n:03d}"
    shutil.rmtree(case_dir, ignore_errors=True)
    case_dir.mkdir(parents=True)
    mz = run_mzdesk(rec, case_dir)
    write_json(case_dir / "mzdesk.json", mz)
    t = time.perf_counter()
    direct = [run_direct(case_dir, call) for call in mz["calls"]]
    runs = [compare_run(case_dir, call, d) for call, d in zip(mz["calls"], direct)]
    result = {"id": n, "accession": rec["accession"], "name": rec["name"], "inchikey": rec["inchikey"], "adduct": rec["adduct"],
              "positive": rec["positive"], "mzdesk_version": version, "code_root": os.environ["MZDESK_CODE_ROOT"],
              "seconds": {"mzdesk": mz["seconds"], "direct": time.perf_counter() - t}, "runs": runs, "layer": compare_layer(case_dir, mz)}
    write_json(case_dir / "result.json", result)


# ================================================================== 4. all cases with a progress bar
def case_verdict(r):
    """identical / differs / no run (MS-FINDER was not started for this spectrum)."""
    if not r["runs"]:
        return "no run"
    ok = all(x["identical"] for x in r["runs"]) and r["layer"]["found"] == r["layer"]["candidates"]
    return "identical" if ok else "differs"


class Progress:
    """One live line at the bottom of the terminal; finished cases are printed above it. Without a terminal: one line per case."""

    def __init__(self, total, workers):
        self.total, self.workers, self.tty = total, workers, sys.stdout.isatty()
        self.counts = {"identical": 0, "differs": 0, "no run": 0, "error": 0}
        self.running, self.durations, self.t0, self.lock = {}, [], time.monotonic(), threading.Lock()

    def done(self):
        return sum(self.counts.values())

    def line(self):
        n, width = self.done(), max(10, min(40, shutil.get_terminal_size((100, 20)).columns - 90))
        fill = int(width * n / self.total) if self.total else width
        el = time.monotonic() - self.t0
        eta = ""
        if self.durations and n < self.total:
            left = (self.total - n) * (sum(self.durations) / len(self.durations)) / self.workers
            eta = f" ~{hms(left)} left"
        now = time.monotonic()
        run = " ".join(f"{k:03d}({hms(now - v)})" for k, v in sorted(self.running.items()))
        c = self.counts
        return (f"[{'#' * fill}{'.' * (width - fill)}] {n}/{self.total}  identical {c['identical']}  differs {c['differs']}  "
                f"no-run {c['no run']}  error {c['error']} | {hms(el)}{eta} | now: {run or '-'}")

    def draw(self):
        if self.tty:
            with self.lock:
                sys.stdout.write("\r\x1b[2K" + self.line()[:shutil.get_terminal_size((100, 20)).columns - 1])
                sys.stdout.flush()

    def start(self, n):
        with self.lock:
            self.running[n] = time.monotonic()

    def finish(self, n, verdict, text, seconds):
        with self.lock:
            self.running.pop(n, None)
            self.counts[verdict] += 1
            self.durations.append(seconds)
            prefix = "\r\x1b[2K" if self.tty else ""
            sys.stdout.write(f"{prefix}  {text}\n")
            if not self.tty:
                sys.stdout.write(f"  {self.line()}\n")
            sys.stdout.flush()
        self.draw()


def run_cases(out, code, ids, workers):
    log = out / "run.log"
    (out / "case_logs").mkdir(exist_ok=True)
    todo = [n for n in ids if not (out / "cases" / f"{n:03d}" / "result.json").exists()]
    with open(log, "a") as f:
        f.write(f"START {stamp()} todo={len(todo)} of {len(ids)} workers={workers} code={code}\n")
    if not todo:
        print(f"All {len(ids)} cases already have results.")
        return
    prog = Progress(len(ids), workers)
    for n in ids:                               # cases finished by an earlier start count towards the bar
        if n not in todo:
            prog.counts[case_verdict(json.loads((out / "cases" / f"{n:03d}" / "result.json").read_text()))] += 1
    if len(todo) < len(ids):
        print(f"{len(ids) - len(todo)} of {len(ids)} cases already have results; running the other {len(todo)}.", flush=True)
    procs, stop = {}, threading.Event()

    def one(n):
        if stop.is_set():
            return n, "error", "not started", 0.0
        prog.start(n)
        t = time.monotonic()
        with open(out / "case_logs" / f"{n:03d}.log", "w") as f:
            p = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), "case", "--out", str(out), "--id", str(n)],
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
        r = json.loads(res.read_text())
        v = case_verdict(r)
        ok_runs = sum(x["identical"] for x in r["runs"])
        lim = sum(x["time_limited"] for x in r["runs"])
        text = (f"{n:03d}  {v:<9} runs {ok_runs}/{len(r['runs'])} identical"
                f"{f' ({lim} time-limited)' if lim else ''}, layer {r['layer']['found']}/{r['layer']['candidates']}  "
                f"[{r['accession']} {r['adduct']}, {hms(sec)}]")
        return n, v, text, sec

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
        sys.exit(130)
    ticker_stop.set()
    th.join()
    if prog.tty:
        sys.stdout.write("\n")
    with open(log, "a") as f:
        f.write(f"DONE {stamp()} {json.dumps(prog.counts)} wall={time.monotonic() - prog.t0:.1f}\n")


# ================================================================== 5. summary
def wilson(k, n, z=1.96):
    if n == 0:
        return None
    p = k / n
    den = 1 + z * z / n
    mid = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, mid - half), min(1.0, mid + half)


def rate(k, n):
    if n == 0:
        return {"k": 0, "n": 0, "rate": None, "ci95": None, "text": "n/a (0)"}
    lo, hi = wilson(k, n)
    return {"k": k, "n": n, "rate": k / n, "ci95": [lo, hi], "text": f"{k}/{n} = {100 * k / n:.1f}%  (95% CI {100 * lo:.1f}-{100 * hi:.1f}%)"}


def cmd_summary(out, ids=None):
    results = []
    for path in sorted((out / "cases").glob("*/result.json")):
        r = json.loads(path.read_text())
        if ids is None or r["id"] in ids:
            results.append(r)
    expected = len(ids) if ids is not None else len(load_sample(out))
    runs = [x for r in results for x in r["runs"]]
    clean = [x for x in runs if not x["time_limited"]]
    cases_run = [r for r in results if r["runs"]]
    cases_clean = [r for r in cases_run if not any(x["time_limited"] for x in r["runs"])]
    cands = sum(r["layer"]["candidates"] for r in results)
    found = sum(r["layer"]["found"] for r in results)
    frozen = {}
    for p in (out / "code/FROZEN.json",):
        if p.exists():
            frozen = json.loads(p.read_text())
    s = {
        "question": "Does MS-FINDER give the same result inside mzDesk as when run on its own (same install, same query and parameter files)?",
        "mzdesk_version": results[0]["mzdesk_version"] if results else None, "commit": frozen.get("commit"),
        "cases_expected": expected, "cases_finished": len(results), "cases_without_msfinder_run": len(results) - len(cases_run),
        "cases_missing_or_failed": expected - len(results),
        "runs": {"all": rate(sum(x["identical"] for x in runs), len(runs)),
                 "without_time_limited": rate(sum(x["identical"] for x in clean), len(clean)),
                 "time_limited": rate(sum(x["identical"] for x in runs if x["time_limited"]), sum(x["time_limited"] for x in runs)),
                 "byte_identical_all": rate(sum(x["byte_identical"] for x in runs), len(runs)),
                 "formula_table_identical": rate(sum(x["formula_table_identical"] for x in runs), len(runs)),
                 "structure_files_identical": rate(sum(x["structure_files_identical"] for x in runs), len(runs))},
        "cases": {"all_runs_identical": rate(sum(all(x["identical"] for x in r["runs"]) for r in cases_run), len(cases_run)),
                  "all_runs_identical_without_time_limited": rate(sum(all(x["identical"] for x in r["runs"]) for r in cases_clean), len(cases_clean)),
                  "layer_fully_reproduced": rate(sum(r["layer"]["found"] == r["layer"]["candidates"] for r in cases_run), len(cases_run))},
        "layer_candidates": rate(found, cands),
        "differences": [{"id": r["id"], "accession": r["accession"], "dir": x["dir"], "time_limited": x["time_limited"],
                         "searches_at_limit": x["searches_at_limit"], "status_mzdesk": x["status_mzdesk"], "status_direct": x["status_direct"],
                         "first_difference": x["first_difference"]} for r in results for x in r["runs"] if not x["identical"]],
        "layer_missing": [{"id": r["id"], **m} for r in results for m in r["layer"]["missing"]],
    }
    write_json(out / "summary.json", s)
    rows = [{"id": r["id"], "accession": r["accession"], "adduct": r["adduct"], "mode": "positive" if r["positive"] else "negative",
             "verdict": case_verdict(r), "msfinder_runs": len(r["runs"]), "runs_identical": sum(x["identical"] for x in r["runs"]),
             "runs_byte_identical": sum(x["byte_identical"] for x in r["runs"]), "runs_time_limited": sum(x["time_limited"] for x in r["runs"]),
             "structures_mzdesk": sum(x["structures_mzdesk"] for x in r["runs"]), "structures_direct": sum(x["structures_direct"] for x in r["runs"]),
             "layer_candidates": r["layer"]["candidates"], "layer_found": r["layer"]["found"], "layer_status": r["layer"]["status"],
             "seconds_mzdesk": round(r["seconds"]["mzdesk"], 1), "seconds_direct": round(r["seconds"]["direct"], 1),
             "first_difference": next((x["first_difference"] for x in r["runs"] if not x["identical"]), "")} for r in results]
    if rows:
        with open(out / "summary.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    lines = [
        f"# MS-FINDER consistency: mzDesk vs standalone",
        "",
        f"mzDesk {s['mzdesk_version']} (commit {(s['commit'] or '?')[:7]}), MS-FINDER v3.73 with mzDesk's thread fix ({MSFINDER_DIR}).",
        f"Cases: {len(results)} of {expected} finished; {s['cases_without_msfinder_run']} without an MS-FINDER run; "
        f"{s['cases_missing_or_failed']} missing or failed.",
        "",
        "| Measure | Consistency |",
        "|---|---|",
        f"| MS-FINDER runs, identical (formulas, structures, scores) | {s['runs']['all']['text']} |",
        f"| ... without time-limited runs | {s['runs']['without_time_limited']['text']} |",
        f"| ... time-limited runs only | {s['runs']['time_limited']['text']} |",
        f"| MS-FINDER runs, byte-identical structure files | {s['runs']['byte_identical_all']['text']} |",
        f"| Cases, every run identical | {s['cases']['all_runs_identical']['text']} |",
        f"| ... without time-limited runs | {s['cases']['all_runs_identical_without_time_limited']['text']} |",
        f"| Cases, mzDesk's MS-FINDER layer fully found in the standalone output | {s['cases']['layer_fully_reproduced']['text']} |",
        f"| Layer candidates found with the same score | {s['layer_candidates']['text']} |",
        "",
    ]
    if s["differences"]:
        lines += ["## Runs that differ", ""] + [
            f"- case {d['id']:03d} {d['dir']}{' (time-limited)' if d['time_limited'] else ''}: {d['first_difference']}" for d in s["differences"]] + [""]
    (out / "report.md").write_text("\n".join(lines) + "\n")
    print("\n" + "\n".join(lines))
    print(f"Files: {out / 'report.md'}, {out / 'summary.csv'}, {out / 'summary.json'}")


# ================================================================== main
def cmd_run(out, workers, ids, limit):
    out.mkdir(parents=True, exist_ok=True)
    print(f"MS-FINDER consistency test -> {out}", flush=True)
    code = freeze_code(out)
    frozen = json.loads((code / "FROZEN.json").read_text()) if (code / "FROZEN.json").exists() else {}
    print(f"code: {code}  {('commit ' + frozen['commit'][:7] + ' ' + frozen.get('subject', '')) if frozen else ''}", flush=True)
    if not thread_fix_ok(MSFINDER_DIR):
        sys.exit(f"{MSFINDER_DIR} does not carry mzDesk's thread fix; run mzDesk's MS-FINDER setup first.")
    if not (out / "inputs/sample.json").exists():
        print("sample: reading MassBank (a few minutes) ...", flush=True)
        subprocess.run([sys.executable, str(Path(__file__).resolve()), "sample", "--out", str(out)], env=child_env(code), check=True)
    all_ids = [r["id"] for r in load_sample(out)]
    chosen = [n for n in all_ids if n in set(ids)] if ids else all_ids
    if limit:
        chosen = chosen[:limit]
    run_cases(out, code, chosen, workers)
    cmd_summary(out, set(chosen) if (ids or limit) else None)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["run", "summary", "sample", "case"])
    ap.add_argument("--out", default=None, help="output folder (default: validation/equivalence/msfinder_<YYYYMMDD>)")
    ap.add_argument("--workers", type=int, default=2, help="cases at once (each runs one MS-FINDER at a time, ~2.3-4 GB)")
    ap.add_argument("--ids", type=int, nargs="*", help="only these case ids (1-50 positive, 51-100 negative)")
    ap.add_argument("--limit", type=int, help="only the first N cases")
    ap.add_argument("--id", type=int, help=argparse.SUPPRESS)
    args = ap.parse_args()
    out = Path(args.out).resolve() if args.out else HERE / f"msfinder_{datetime.now():%Y%m%d}"
    if args.command == "run":
        cmd_run(out, args.workers, args.ids, args.limit)
    elif args.command == "summary":
        cmd_summary(out)
    elif args.command == "sample":
        cmd_sample(out)
    else:
        cmd_case(out, args.id)


if __name__ == "__main__":
    main()
