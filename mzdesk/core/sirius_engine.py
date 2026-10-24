"""SIRIUS 6 molecular formula identification engine (optional second engine).

Wraps the external `sirius` CLI (Bright Giant) for the Structure Annotation tab.
SIRIUS is OPTIONAL: the tab degrades gracefully when it is missing or not logged in.

Input contract (SIRIUS 6.3, https://v6.docs.sirius-ms.io):
  - builds a single-compound `.ms` peak-list file (>compound/>feature_id/>parentmass/
    >ion, optional >formula/>ms1/>ms2 blocks) from the tab's inputs;
  - runs `sirius --input x.ms --project <tmp> formulas -p <profile> -c N
    [--formulas F] write-summaries --output <dir>`.
  - SIRIUS 6 refuses to compute without an account login (`sirius login`): it then exits
    with code 0, prints the login help on stdout and "Login ERROR" on stderr. That case
    is raised as SiriusLoginError with instructions.

Output: the top-ranked formula from `formula_identifications.tsv` (SIRIUS only writes
the top hit there by default); a ranked listing on stdout (older SIRIUS versions) is
used as a last resort.
"""
import csv
import glob
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import streamlit as st

# sirius instrument profiles: qtof (default, 10 ppm) / orbitrap (5 ppm). SIRIUS 6.3.5
# lists 'fticr' in `formulas --help` but ships no fticr profile ("Could not find profile
# JSON"), so FT-ICR data use the orbitrap profile.
SIRIUS_PROFILE = {"qtof": "qtof", "orbitrap": "orbitrap", "fticr": "orbitrap"}

MAX_CANDIDATES = 5

LOGIN_HINT = (
    "SIRIUS 6 only computes for a logged-in SIRIUS account (free for academic use). "
    "Run `sirius login -u <email> -p` once in a terminal, then click Annotate formula again."
)

# a Hill-formula token used when parsing the console candidate listing
_FORMULA_TOKEN = re.compile(r"^((?:[A-Z][a-z]?\d*)+)$")
_RANK_LINE = re.compile(r"^\s*(\d{1,2})[.)\]]?\s+(.+)$")
_LOGIN_ERROR = re.compile(r"Login ERROR|Please Login", re.IGNORECASE)
_LOG_LEVEL = re.compile(r"^\s*(SEVERE|ERROR|WARNING)\b")


class SiriusLoginError(RuntimeError):
    """SIRIUS 6 refused to run because no valid account login is stored."""


def sirius_available():
    """Path to the `sirius` executable, or None if not found.

    Lookup order: $SIRIUS_EXE -> PATH -> common macOS app-bundle locations.
    """
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
        hits = sorted(glob.glob(pattern))
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
    """Parse a ranked candidate list printed to the console (older SIRIUS versions).

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


def _num(value):
    """float(value) for numeric TSV cells, None for empty / N/A / NaN."""
    try:
        v = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return v if v == v else None


def _parse_summary_tsv(summary_path):
    """Parse formula_identifications.tsv (the authoritative result source).

    SIRIUS 6 exports only the TOP-ranked formula here. Columns used (checked against
    the SIRIUS 6.3.5 summary writer): formulaRank, molecularFormula, adduct,
    SiriusScore, SiriusScoreNormalized, TreeScore, IsotopeScore, numExplainedPeaks,
    explainedIntensity (a 0-1 ratio), massErrorPrecursor(ppm).
    Returns a de-duplicated list of candidate dicts with numeric fields as floats.
    """
    rows = []
    with open(summary_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            f = (row.get("molecularFormula") or "").strip()
            if not f:
                continue
            cand = {
                "rank": int(_num(row.get("formulaRank")) or 1),
                "formula": f,
                "adduct": (row.get("adduct") or "").replace(" ", ""),
                "score": _num(row.get("SiriusScore")),
                "score_norm": _num(row.get("SiriusScoreNormalized")),
                "tree_score": _num(row.get("TreeScore")),
                "iso_score": _num(row.get("IsotopeScore")),
                "explained_peaks": _num(row.get("numExplainedPeaks")),
                "explained_int": _num(row.get("explainedIntensity")),
                "mass_error_ppm": _num(row.get("massErrorPrecursor(ppm)")),
            }
            if not any(c["rank"] == cand["rank"] and c["formula"] == cand["formula"]
                       and c["adduct"] == cand["adduct"] for c in rows):
                rows.append(cand)
    return sorted(rows, key=lambda c: c["rank"])


def _diagnostics(stdout_text, stderr_text, n_lines=8):
    """Last SEVERE/ERROR/WARNING log lines (stderr), else the tail of stdout."""
    lines = [ln.strip() for ln in (stderr_text or "").splitlines() if _LOG_LEVEL.match(ln)]
    if not lines:
        lines = [ln.rstrip() for ln in (stdout_text or "").splitlines() if ln.strip()]
    return "\n".join(lines[-n_lines:])


@st.cache_data(show_spinner=False)
def run_sirius(precursor_mz, positive, adduct, formula_input,
               ms2_peaks, ms1_peaks,
               ms_instr="orbitrap", timeout_secs=600, n_candidates=MAX_CANDIDATES):
    """Run SIRIUS `formulas` for one query (``timeout_secs`` = wall-clock limit).

    Returns dict:
      candidates : [{rank, formula, adduct, score, score_norm, tree_score, iso_score,
                     explained_peaks, explained_int, mass_error_ppm}] (console fallback:
                     rank/formula/score/adduct only)
      source     : "summary" | "console"
      log_tail   : last warning/error lines (for troubleshooting)
    Raises SiriusLoginError when SIRIUS needs a login, RuntimeError on other failures,
    so the caller can degrade gracefully.
    """
    exe = sirius_available()
    if exe is None:
        raise RuntimeError("SIRIUS executable not found (set $SIRIUS_EXE or put `sirius` on PATH).")

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
            raise RuntimeError(f"SIRIUS timed out after {int(timeout_secs)} s "
                               "(Advanced → SIRIUS timeout).")
        stdout, stderr = proc.stdout or "", proc.stderr or ""
        log_tail = _diagnostics(stdout, stderr)
        login_failed = bool(_LOGIN_ERROR.search(stderr) or _LOGIN_ERROR.search(stdout))
        if proc.returncode != 0:
            if login_failed:
                raise SiriusLoginError(LOGIN_HINT)
            raise RuntimeError(f"SIRIUS exited with code {proc.returncode}.\n{log_tail}")

        candidates, source = [], "summary"
        summary_file = outdir / "formula_identifications.tsv"
        if summary_file.exists():
            candidates = _parse_summary_tsv(summary_file)
        if not candidates and login_failed:
            raise SiriusLoginError(LOGIN_HINT)
        if not candidates:
            # last resort: console listing (older SIRIUS versions)
            source = "console"
            for rank, formula, score in _parse_stdout_candidates(stdout):
                candidates.append({"rank": rank, "formula": formula,
                                   "score": _num(score), "adduct": ""})

        if not candidates:
            raise RuntimeError("SIRIUS produced no formula result (no candidate passed its "
                               "filters, or the output could not be parsed).\n" + log_tail)
        return {"candidates": candidates[:MAX_CANDIDATES], "source": source,
                "log_tail": log_tail}
