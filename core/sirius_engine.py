"""SIRIUS 6 molecular formula identification engine (optional second engine).

Wraps the external `sirius` CLI (Bright Giant) for the Structure Annotation tab.
SIRIUS is OPTIONAL: the tab degrades gracefully when it is not installed.

Input contract (per https://v6.docs.sirius-ms.io):
  - builds a single-compound `.ms` peak-list file (>compound/>parentmass/>ion/
    optional >formula/>ms1/>ms2 blocks) from the tab's inputs;
  - runs `sirius --input x.ms --project <tmp> formulas [-c N]
    [--formulas F] write-summaries --output <dir>`.

Output: ranked formula candidates parsed from the console listing when
available, otherwise the top hit from `formula_identifications.tsv`.
"""
import csv
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import streamlit as st

# sirius instrument profiles: qtof (default, 10 ppm) / orbitrap (5 ppm);
# fticr data should use the orbitrap profile with a lower ppm-max (docs).
SIRIUS_PROFILE = {"qtof": "qtof", "orbitrap": "orbitrap", "fticr": "orbitrap"}

MAX_CANDIDATES = 5

# a Hill-formula token used when parsing the console candidate listing
_FORMULA_TOKEN = re.compile(r"^((?:[A-Z][a-z]?\d*)+)$")
_RANK_LINE = re.compile(r"^\s*(\d{1,2})[.)\]]?\s+(.+)$")


def sirius_available():
    """Path to the `sirius` executable, or None if not found.

    Lookup order: $SIRIUS_EXE -> PATH -> common macOS app-bundle locations.
    """
    import os
    exe = os.environ.get("SIRIUS_EXE")
    if exe and Path(exe).is_file():
        return exe
    on_path = shutil.which("sirius")
    if on_path:
        return on_path
    for pattern in (
            str(Path.home() / "tools" / "sirius*.app" / "Contents" / "MacOS" / "sirius"),
            "/Applications/sirius*.app/Contents/MacOS/sirius",
            "/Applications/SIRIUS*.app/Contents/MacOS/sirius",
    ):
        import glob as _glob
        hits = sorted(_glob.glob(pattern))
        if hits:
            return hits[0]
    return None


def auto_adduct(positive):
    """SIRIUS wildcard adduct when the user leaves Adduct on Auto."""
    return "[M+?]+" if positive else "[M+?]-"


def _build_ms_file(path, precursor_mz, positive, adduct, formula_input,
                   ms2_peaks, ms1_peaks):
    lines = [
        ">compound query",
        ">feature_id query",
        f">parentmass {precursor_mz:.6f}",
        f">ion {adduct or auto_adduct(positive)}",
    ]
    if formula_input:
        lines.append(f">formula {formula_input}")
    if ms1_peaks:
        lines.append(">ms1")
        lines += [f"{mz:.6f} {inten:.6f}" for mz, inten in ms1_peaks]
    if ms2_peaks:
        lines.append(">ms2")
        lines += [f"{mz:.6f} {inten:.6f}" for mz, inten in ms2_peaks]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _looks_like_formula(tok):
    tok = tok.strip("[],")
    return bool(_FORMULA_TOKEN.match(tok)) and any(c.isdigit() for c in tok)


def _parse_stdout_candidates(stdout_text):
    """Parse the ranked candidate list SIRIUS prints to the console.

    Returns [(rank, formula, raw_score_str), ...] or [].
    """
    out = []
    for line in stdout_text.splitlines():
        m = _RANK_LINE.match(line)
        if not m:
            continue
        rest = m.group(2).split()
        if not rest:
            continue
        # first token of a rank line is expected to be a formula like C8H10N4O2
        cand = rest[0].strip("[],")
        if _looks_like_formula(cand):
            score = rest[1] if len(rest) > 1 else ""
            out.append((int(m.group(1)), cand, score))
    return sorted(out)[: MAX_CANDIDATES]


def _parse_summary_tsv(summary_path):
    """Parse formula_identifications.tsv (the authoritative result source).

    SIRIUS 6 exports only the TOP-ranked formula here; columns are rich:
    molecularFormula, adduct, SiriusScore, SiriusScoreNormalized, TreeScore,
    IsotopeScore, numExplainedPeaks, explainedIntensity, massErrorPrecursor...
    Returns a de-duplicated list of candidate dicts.
    """
    rows = []
    with open(summary_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            f = (row.get("molecularFormula") or "").strip()
            if not f:
                continue
            cand = {
                "rank": int((row.get("formulaRank") or "1").strip() or 1),
                "formula": f,
                "adduct": (row.get("adduct") or "").strip(),
                "score": (row.get("SiriusScore") or "").strip(),
                "score_norm": (row.get("SiriusScoreNormalized") or "").strip(),
                "tree_score": (row.get("TreeScore") or "").strip(),
                "iso_score": (row.get("IsotopeScore") or "").strip(),
                "explained_peaks": (row.get("numExplainedPeaks") or "").strip(),
                "explained_int": (row.get("explainedIntensity") or "").strip(),
            }
            if not any(c["rank"] == cand["rank"] and c["formula"] == cand["formula"]
                       and c["score"] == cand["score"] for c in rows):
                rows.append(cand)
    return rows


@st.cache_data(show_spinner=False)
def run_sirius(precursor_mz, positive, adduct, formula_input,
               ms2_peaks, ms1_peaks,
               ms_instr="orbitrap", timeout_secs=600, n_candidates=5):
    """Run SIRIUS `formulas` for one query. Returns dict:

      candidates : [{rank, formula, score, adduct}]
      source     : "console" | "summary"
      log_tail   : last lines of console output (for error reporting)
    Raises RuntimeError on failure so the caller can degrade gracefully.
    """
    exe = sirius_available()
    if exe is None:
        raise RuntimeError("SIRIUS executable not found on PATH.")

    with tempfile.TemporaryDirectory(prefix="sirius_tab_") as tmp:
        tmp = Path(tmp)
        ms_file = tmp / "query.ms"
        project = tmp / "project"
        outdir = tmp / "summary"
        _build_ms_file(ms_file, float(precursor_mz), bool(positive),
                       adduct, formula_input, ms2_peaks, ms1_peaks)
        cmd = [
            exe, "--input", str(ms_file), "--project", str(project),
            "formulas",
            "-p", SIRIUS_PROFILE.get(ms_instr, "orbitrap"),
            "-c", str(max(int(n_candidates), 1)),
        ]
        if formula_input:
            cmd += ["--formulas", formula_input]
        cmd += ["write-summaries", "--output", str(outdir)]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True,
                                  timeout=int(timeout_secs))
        except subprocess.TimeoutExpired:
            raise RuntimeError(f"SIRIUS timed out after {timeout_secs}s.")
        log_tail = "\n".join((proc.stdout or "").splitlines()[-12:])
        if proc.returncode != 0:
            raise RuntimeError(
                f"SIRIUS exited with code {proc.returncode}.\n{log_tail}")

        candidates = []
        summary_file = outdir / "formula_identifications.tsv"
        if summary_file.exists():
            candidates = _parse_summary_tsv(summary_file)
        if not candidates:
            # last resort: console listing (older SIRIUS versions)
            for rank, formula, score in _parse_stdout_candidates(proc.stdout or ""):
                candidates.append({"rank": rank, "formula": formula,
                                   "score": score, "adduct": ""})

        if not candidates:
            raise RuntimeError(
                "SIRIUS produced no parsable formula result.\n" + log_tail)
        return {"candidates": candidates[:MAX_CANDIDATES], "source": "summary",
                "log_tail": log_tail}
