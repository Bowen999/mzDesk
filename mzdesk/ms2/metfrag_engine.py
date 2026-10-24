"""MetFrag engine for the MS² Interpretation tab.

Runs the bundled MetFrag CLI (mzdesk/lib/MetFragCommandLine-2.5.2.jar) on ONE user-supplied
candidate structure and maps MetFrag's explained peaks back onto the user's peak list
(ion formula, theoretical m/z, mass error, neutral loss, fragment structure).

The shared peak-list parser ``parse_spectrum`` lives in mzdesk/shared/spectrum.py and is
re-exported here. Also hosts the interactive spectrum chart (Plotly inside a components.html iframe with a custom
hover card that shows the fragment structure).
"""
import base64
import json
import math
import re
import subprocess
import tempfile
from pathlib import Path

import plotly
import plotly.graph_objects as go
import streamlit as st
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.MolStandardize import rdMolStandardize

from mzdesk.shared.spectrum import parse_spectrum  # noqa: F401  re-exported: ui_fragments uses engine.parse_spectrum

RDLogger.DisableLog('rdApp.*')

JAR = Path(__file__).parent.parent / "lib" / "MetFragCommandLine-2.5.2.jar"
METFRAG_TIMEOUT_S = 1200


@st.cache_resource
def get_plotlyjs():
    return plotly.offline.get_plotlyjs()


# NOTE: the ~4.9 MB of Plotly below is re-sent with the iframe on every rerun. Serving it
# from ./static (server.enableStaticServing) does NOT work: Streamlit returns app static
# files as text/plain, so the browser refuses to execute the script ("Plotly is not
# defined"). A CDN would work but would break offline use; switching to st.plotly_chart
# would drop the hover card with the fragment structure image.


# --------------------------------------------------------------------------
# Masses, formulas, adducts
# --------------------------------------------------------------------------
# Monoisotopic (most abundant isotope) masses.
ATOMIC = {'C': 12.0, 'H': 1.00782503, 'O': 15.99491462, 'P': 30.97376199, 'N': 14.0030740,
          'S': 31.9720707, 'F': 18.99840322, 'Cl': 34.96885268, 'Br': 78.9183376,
          'I': 126.904473, 'Na': 22.98976928, 'K': 38.96370649, 'Si': 27.9769265,
          'B': 11.0093054, 'Se': 79.9165218, 'D': 2.01410178}
ELECTRON = 0.000548579909
DARK = "#1B3A8F"
LIGHT = "#A9BEEA"

# PrecursorIonMode codes are MetFrag's adduct nominal masses (Constants.ADDUCT_NOMINAL_MASSES);
# "shift" = m/z(ion) - M (electron mass included); "delta" = ion composition minus M.
ADDUCTS = {
    "[M-H]-":      {"mode": -1, "positive": False, "shift": -1.0072764, "delta": {"H": -1}},
    "[M+Cl]-":     {"mode": 35, "positive": False, "shift": 34.9694011, "delta": {"Cl": 1}},
    "[M+HCOO]-":   {"mode": 45, "positive": False, "shift": 44.9982010, "delta": {"C": 1, "H": 1, "O": 2}},
    "[M+CH3COO]-": {"mode": 59, "positive": False, "shift": 59.0138510, "delta": {"C": 2, "H": 3, "O": 2}},
    "[M+H]+":      {"mode":  1, "positive": True,  "shift":  1.0072764, "delta": {"H": 1}},
    "[M+Na]+":     {"mode": 23, "positive": True,  "shift": 22.9892214, "delta": {"Na": 1}},
    "[M+NH4]+":    {"mode": 18, "positive": True,  "shift": 18.0338254, "delta": {"N": 1, "H": 4}},
    "[M+K]+":      {"mode": 39, "positive": True,  "shift": 38.9631578, "delta": {"K": 1}},
    # [M]+ : the drawn structure IS the observed ion (permanently charged cation, or a
    # radical cation from a neutral molecule). run_metfrag drops the electron shift when the
    # input already carries the charge, because RDKit's ExactMolWt accounts for it.
    "[M]+":        {"mode":  0, "positive": True,  "shift": -ELECTRON,  "delta": {}},
}

# Adduct groups MetFrag may append to a fragment ion (Constants.ADDUCT_NAMES) use plain
# formulas except acetonitrile.
_GROUP_ALIASES = {"ACN": {"C": 2, "H": 3, "N": 1}}


def _element_mass(el):
    mass = ATOMIC.get(el)
    if mass is None:
        try:
            mass = Chem.GetPeriodicTable().GetMostCommonIsotopeMass(el)
        except Exception:
            mass = 0.0
    return mass


def _formula_counts(formula):
    """'C18H33O2' -> {'C': 18, 'H': 33, 'O': 2} (no brackets / charges)."""
    counts = {}
    for el, n in re.findall(r'([A-Z][a-z]?)(\d*)', formula or ""):
        counts[el] = counts.get(el, 0) + (int(n) if n else 1)
    return counts


def _add_counts(counts, delta, factor=1):
    out = dict(counts)
    for el, n in delta.items():
        out[el] = out.get(el, 0) + factor * n
    return out


def decode_metfrag_ion(s):
    """Decode a MetFrag fragment-ion string -> (ion element counts, H shift inside the brackets).

    MetFrag writes '[' + fragment formula + H shift + ']' + adduct + charge sign, e.g.
    '[C3H7O4P+H]-', '[C18H33O2-H]-H-', '[C2H4O4P+H]+Na+', '[C39H73O8P]+Na+', '[C10H20O2]+2ACN+H+'.
    The H shift ('+H', '+2H', '-H', '-2H') and every adduct group (-H, +H, +NH4, +Na, +K, +Cl,
    +HCOO, +CH3COO, -2H, -2H+Na, -2H+K, +CH3OH+H, +ACN+H, +ACN+Na, +2ACN+H) are folded into
    the returned composition, which is therefore the composition of the detected ion.
    """
    m = re.match(r'\s*\[([^\]]+)\](.*)$', s or "")
    if not m:
        return {}, 0
    inner, tail = m.group(1).strip(), m.group(2).strip()
    im = re.match(r'^(.*?)([+-])(\d*)H$', inner)
    hshift = 0
    if im and im.group(1):
        inner = im.group(1)
        hshift = (1 if im.group(2) == '+' else -1) * (int(im.group(3)) if im.group(3) else 1)
    counts = _add_counts(_formula_counts(inner), {'H': hshift})
    if tail and tail[-1] in '+-':
        tail = tail[:-1]                                   # trailing charge sign
    for sign, mult, name in re.findall(r'([+-])(\d*)([A-Za-z][A-Za-z0-9]*)', tail):
        group = _GROUP_ALIASES.get(name) or _formula_counts(name)
        counts = _add_counts(counts, group, (1 if sign == '+' else -1) * (int(mult) if mult else 1))
    return counts, hshift


def ion_mz(counts, positive):
    """m/z of a singly-charged ion. Cation: -1 electron; anion: +1 electron."""
    mass = sum(_element_mass(el) * n for el, n in counts.items())
    return mass + (-ELECTRON if positive else ELECTRON)


def formula_str(counts):
    """Hill-order formula string (C, H, then alphabetical) — same convention as RDKit's CalcMolFormula."""
    els = [el for el, n in counts.items() if n]
    if 'C' in els:
        order = ['C'] + (['H'] if 'H' in els else []) + sorted(e for e in els if e not in ('C', 'H'))
    else:
        order = sorted(els)
    return ''.join(el + (str(counts[el]) if counts[el] != 1 else '') for el in order)


# --------------------------------------------------------------------------
# Structure rendering
# --------------------------------------------------------------------------
def render_frag_png(smi, size=(240, 170)):
    """Best-effort 2D depiction of a (possibly artifact-laden) MetFrag fragment SMILES -> base64 PNG str."""
    if not smi:
        return ""
    m = Chem.MolFromSmiles(smi)
    if m is None:
        m = Chem.MolFromSmiles(smi, sanitize=False)
        if m is None:
            return ""
        try:
            m.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(m, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^
                             Chem.SanitizeFlags.SANITIZE_PROPERTIES ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
        except Exception:
            pass
    try:
        d = rdMolDraw2D.MolDraw2DCairo(*size)
        rdMolDraw2D.PrepareAndDrawMolecule(d, m)
        d.FinishDrawing()
        return base64.b64encode(d.GetDrawingText()).decode()
    except Exception:
        return ""


@st.cache_data(show_spinner=False, max_entries=32)
def render_input_png(smiles, size=(560, 360)):
    """Render the INPUT molecule to PNG bytes for display (None if the SMILES is invalid)."""
    if not smiles or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    d = rdMolDraw2D.MolDraw2DCairo(*size)
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    return d.GetDrawingText()


# --------------------------------------------------------------------------
# Input structure / precursor
# --------------------------------------------------------------------------
def _prepare_molecule(smiles, allow_charge=False):
    """Parse the candidate structure -> (mol, dropped_smiles, None) or (None, dropped, error).

    A salt / mixture keeps its largest fragment (heavy atoms, then exact mass) and the dropped
    parts are returned so the caller can say what was ignored. Unless ``allow_charge`` (the
    [M]+ adduct, where the drawn structure is itself the ion) the structure is neutralised and
    a remaining net charge is an error — MetFrag adds the adduct itself.
    """
    smiles = (smiles or "").strip()
    if not smiles:
        return None, [], "Enter the SMILES of the candidate structure."
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, [], "Invalid SMILES string. Please check the structure."
    dropped = []
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        ranked = sorted(frags, key=lambda f: (f.GetNumHeavyAtoms(), Descriptors.ExactMolWt(f)),
                        reverse=True)
        mol = ranked[0]
        dropped = [Chem.MolToSmiles(f) for f in ranked[1:]]
    if not allow_charge:
        try:
            mol = rdMolStandardize.Uncharger().uncharge(mol)
        except Exception:
            pass
        charge = Chem.GetFormalCharge(mol)
        if charge != 0:
            return None, dropped, (
                f"The input structure carries a net charge that could not be neutralized "
                f"(net {charge:+d}). Provide the NEUTRAL molecule's SMILES — MetFrag adds the "
                "adduct itself — or pick the [M]+ adduct if the molecule is permanently charged.")
    return mol, dropped, None


def _adduct_shift(add, mol):
    """m/z(ion) - ExactMolWt(structure). Zero for [M]+ of an already-charged structure:
    RDKit's exact mass has the missing electron accounted for already."""
    if add["mode"] == 0 and Chem.GetFormalCharge(mol) > 0:
        return 0.0
    return add["shift"]


@st.cache_data(show_spinner=False, max_entries=32)
def precursor_info(smiles, adduct):
    """Neutral formula / theoretical precursor m/z of a structure for an adduct (no MetFrag run).

    Returns {"error": str|None, "formula", "mono_mass", "theo_mz", "ion_formula"}.
    """
    mol, dropped, err = _prepare_molecule(smiles, allow_charge=(adduct == "[M]+"))
    if mol is None:
        return {"error": err, "dropped": dropped}
    formula = rdMolDescriptors.CalcMolFormula(mol)
    mono = Descriptors.ExactMolWt(mol)
    add = ADDUCTS[adduct]
    return {"error": None, "formula": formula, "mono_mass": mono, "dropped": dropped,
            "theo_mz": mono + _adduct_shift(add, mol),
            "ion_formula": formula_str(_add_counts(_formula_counts(formula), add["delta"]))}


# --------------------------------------------------------------------------
# MetFrag run
# --------------------------------------------------------------------------
class MetFragError(RuntimeError):
    """MetFrag could not be run or failed; ``log`` holds the tail of its console output."""

    def __init__(self, message, log=""):
        super().__init__(message)
        self.log = log


def _psv_mapping(field):
    out = {}
    for entry in (field or "").split(";"):
        k, sep, v = entry.partition(":")
        if not sep:
            continue
        try:
            out[float(k)] = v.strip()
        except ValueError:
            continue
    return out


# 'mz:value' entry, and a complete MetFrag ion string such as '[C18H33O2-H]-H-' / '[C2H4O4P+H]+Na+'.
_PSV_ENTRY = re.compile(r'^\d+(?:\.\d+)?(?:[eE][-+]?\d+)?:\S')
_PSV_ION = re.compile(r'^\[[^\[\]]+\](?:[+-]\d*[A-Za-z][A-Za-z0-9]*)*[+-]$')


def parse_metfrag_psv(psv_path):
    """Parse MetFrag's FragmentSmilesPSV output (first candidate row).

    Returns (formulas, smiles): dicts keyed by explained peak m/z (float) -> MetFrag ion-formula
    string (FormulasOfExplPeaks) / fragment SMILES (SmilesOfExplPeaks). Both are empty when
    MetFrag explained nothing or returned no candidate.

    Columns are recognised by content, not by header name: MetFrag 2.5.2 writes one more data
    field than there are header names (everything after 'Identifier' is shifted by one).
    The formula column is the field whose every entry is 'mz:<ion string>'; the SMILES column
    is the first other field made of 'mz:value' entries.
    """
    lines = [ln for ln in Path(psv_path).read_text(encoding='utf-8').splitlines() if ln.strip()]
    if len(lines) < 2:
        return {}, {}
    formulas, smiles = {}, {}
    for field in lines[1].split('|'):
        entries = [e.strip() for e in field.split(';') if e.strip()]
        if not entries or not all(_PSV_ENTRY.match(e) for e in entries):
            continue
        mapping = _psv_mapping(field)
        if all(_PSV_ION.match(v) for v in mapping.values()):
            formulas = formulas or mapping
        else:
            smiles = smiles or mapping
    return formulas, smiles


def _nearest_key(keys, mz, tol=1e-4):
    """Explained-peak key (MetFrag echoes the m/z we wrote, 5 dp) closest to mz within tol."""
    best, best_d = None, tol
    for k in keys:
        d = abs(k - mz)
        if d <= best_d:
            best, best_d = k, d
    return best


DEFAULT_SPECTRUM = """139.00977 3
152.99622 24
168.04965 3
279.22757 48
280.2258 27
281.24707 1101
282.24652 24
283.26199 14
491.25638 10
509.2948 17
599.81622 3
713.46429 10
773.53589 416"""


def run_metfrag(smiles, peaks, adduct, precursor_mz, tree_depth, tol_abs, tol_ppm):
    """Run MetFrag on one candidate and return per-peak annotations.

    ``precursor_mz`` is the measured precursor (None = use the theoretical m/z of the
    structure + adduct); it is used for neutral losses and a mismatch warning. MetFrag
    itself only receives the neutral formula and the ion mode. A peak is matched with
    MetFrag's own rule |Δm| <= tol_abs + tol_ppm·m/z·1e-6.
    """
    if adduct not in ADDUCTS:
        raise ValueError(f"Unsupported adduct {adduct!r}.")
    peaks = [(float(mz), float(i)) for mz, i in peaks
             if math.isfinite(mz) and math.isfinite(i) and mz > 0 and i > 0]
    if not peaks:
        raise ValueError("No valid peaks found in the spectrum. Check the peak list format.")
    mol, dropped, err = _prepare_molecule(smiles, allow_charge=(adduct == "[M]+"))
    if mol is None:
        raise ValueError(err)
    neutral_smiles = Chem.MolToSmiles(mol)
    inchi = Chem.MolToInchi(mol)
    inchikey = Chem.InchiToInchiKey(inchi) if inchi else ""
    mol_formula = rdMolDescriptors.CalcMolFormula(mol)
    mono_mass = Descriptors.ExactMolWt(mol)
    add = ADDUCTS[adduct]
    positive = add["positive"]
    theo_mz = mono_mass + _adduct_shift(add, mol)
    ion_formula = formula_str(_add_counts(_formula_counts(mol_formula), add["delta"]))

    precursor_given = precursor_mz is not None and precursor_mz != ""
    precursor_mz = float(precursor_mz) if precursor_given else theo_mz
    warnings = []
    if dropped:
        warnings.append(
            f"The SMILES contained {len(dropped) + 1} disconnected parts; only the largest was "
            f"annotated. Ignored: {', '.join(dropped)}.")
    if precursor_given and abs(precursor_mz - theo_mz) > 0.5:
        warnings.append(
            f"Precursor m/z {precursor_mz:.4f} does not match the structure's {adduct} m/z "
            f"{theo_mz:.4f} (Δ {precursor_mz - theo_mz:+.3f} Da). "
            "Check the SMILES, adduct, or precursor value — annotation may be poor.")

    if not JAR.exists():
        raise MetFragError(f"MetFrag JAR not found at {JAR}.")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        (tmp / "results").mkdir(exist_ok=True)
        candidate_path = tmp / "candidate.csv"
        with open(candidate_path, "w") as f:
            f.write('"Identifier","InChI","MonoisotopicMass","MolecularFormula","SMILES","InChIKey","CompoundName"\n')
            f.write(f'"candidate","{inchi}","{mono_mass:.6f}","{mol_formula}","{neutral_smiles}","{inchikey}","candidate"\n')
        peaks_path = tmp / "peaks.txt"
        with open(peaks_path, "w") as f:
            for mz, inten in peaks:
                f.write(f"{mz:.5f} {inten:.10g}\n")
        param_path = tmp / "param.txt"
        params = [
            f"PeakListPath = {peaks_path}",
            "MetFragDatabaseType = LocalCSV",
            f"LocalDatabasePath = {candidate_path}",
            f"NeutralPrecursorMolecularFormula = {mol_formula}",
            f"PrecursorIonMode = {add['mode']}",
            f"IsPositiveIonMode = {positive}",
            "MetFragCandidateWriter = FragmentSmilesPSV",
            "SampleName = sample",
            f"ResultsPath = {tmp / 'results'}",
            f"MaximumTreeDepth = {tree_depth}",
            f"FragmentPeakMatchAbsoluteMassDeviation = {tol_abs}",
            f"FragmentPeakMatchRelativeMassDeviation = {tol_ppm}",
            "DatabaseSearchRelativeMassDeviation = 10",
            "MetFragPreProcessingCandidateFilter = UnconnectedCompoundFilter",
        ]
        with open(param_path, "w") as f:
            f.write("\n".join(params))
        try:
            proc = subprocess.run(
                ["java", "-jar", str(JAR), str(param_path)],
                capture_output=True, text=True, timeout=METFRAG_TIMEOUT_S)
        except FileNotFoundError:
            raise MetFragError("Java was not found on PATH — MetFrag needs a Java 8+ runtime.") from None
        except subprocess.TimeoutExpired:
            raise MetFragError(
                f"MetFrag timed out after {METFRAG_TIMEOUT_S} s — try a lower tree depth.") from None
        log = ((proc.stdout or "") + (proc.stderr or "")).strip()
        if proc.returncode != 0:
            raise MetFragError(f"MetFrag failed (exit code {proc.returncode}).", log[-4000:])
        psv_path = tmp / "results" / "sample.psv"
        if not psv_path.exists():
            raise MetFragError("MetFrag finished without writing a result file.", log[-4000:])
        n_rows = sum(1 for ln in psv_path.read_text(encoding="utf-8").splitlines() if ln.strip())
        formulas, smiles_dict = parse_metfrag_psv(psv_path)
    if n_rows < 2:
        warnings.append("MetFrag returned no result for this structure, so no fragment could be "
                        "annotated.")

    max_inten = max(i for _, i in peaks)
    tol_prec = tol_abs + theo_mz * tol_ppm / 1e6        # MetFrag's peak window (absolute + relative)
    keys = sorted(formulas)
    precursor_png = None
    annotation = []
    for mz, inten in peaks:
        row = {"mz": mz, "intensity": inten, "rel": inten / max_inten * 100.0,
               "annotated": False, "kind": "", "ion_formula": "", "metfrag_ion": "",
               "calc_mz": None, "ppm": None, "nl": None, "smiles": "", "b64": ""}
        if abs(mz - theo_mz) <= tol_prec:
            if precursor_png is None:
                precursor_png = render_frag_png(neutral_smiles)
            row.update(annotated=True, kind="precursor", ion_formula=ion_formula,
                       calc_mz=theo_mz, ppm=(mz - theo_mz) / theo_mz * 1e6, nl=0.0,
                       smiles=neutral_smiles, b64=precursor_png)
        else:
            key = _nearest_key(keys, mz)
            if key is not None:
                ion = formulas[key]
                counts, _ = decode_metfrag_ion(ion)
                calc = ion_mz(counts, positive) if counts else None
                frag_smiles = smiles_dict.get(key, "")
                row.update(annotated=True, kind="fragment", metfrag_ion=ion,
                           ion_formula=formula_str(counts) if counts else ion,
                           calc_mz=calc, ppm=((mz - calc) / calc * 1e6) if calc else None,
                           nl=precursor_mz - mz, smiles=frag_smiles,
                           b64=render_frag_png(frag_smiles))
        annotation.append(row)
    return {
        "mol_formula": mol_formula, "mono_mass": mono_mass, "ion_formula": ion_formula,
        "adduct": adduct, "theo_precursor_mz": theo_mz,
        "precursor_mz": precursor_mz, "precursor_given": precursor_given,
        "positive": positive, "tree_depth": tree_depth, "tol_abs": tol_abs, "tol_ppm": tol_ppm,
        "annotation": annotation, "warnings": warnings,
        "warning": warnings[0] if warnings else None}


# --------------------------------------------------------------------------
# Interactive spectrum (Plotly in a components.html iframe)
# --------------------------------------------------------------------------
FONT = '"Helvetica Neue",Helvetica,Arial,sans-serif'   # same stack as theme.PLOT_FONT
MONO = 'ui-monospace,"SF Mono",Menlo,Consolas,monospace'
CHART_HEIGHT = 460
Y_FLOOR = 1.6      # sticks start slightly below 0 % so even tiny peaks stay visible above the axis line
Y_TOP = 116.0      # headroom for the m/z labels above the base peak


def build_figure(result, dark=DARK, light=LIGHT, ink="#0A1633", muted="#5C6B8A",
                 grid="#E9EDF7", font=FONT, mono=MONO, height=CHART_HEIGHT):
    """Stick spectrum: dark = annotated, light = unannotated, dotted line = theoretical precursor.

    Sticks are drawn as lines (constant pixel width), so neighbouring peaks never merge when
    zooming and every peak stays visible when zoomed out. An invisible marker per peak (trace
    'peaks', customdata = peak index) drives the hover card in hovermode 'x'.
    """
    ann = result["annotation"]
    fig = go.Figure()
    if not ann:
        return fig
    theo = result.get("theo_precursor_mz") or result.get("precursor_mz")
    xs = [a["mz"] for a in ann]
    lo, hi = min(xs), max(xs)
    if theo:
        lo, hi = min(lo, theo), max(hi, theo)
    pad = max((hi - lo) * 0.05, 5.0)
    x0, x1 = max(0.0, lo - pad), hi + pad

    for name, color, rows in (("unannotated", light, [a for a in ann if not a["annotated"]]),
                              ("annotated", dark, [a for a in ann if a["annotated"]])):
        if not rows:
            continue
        sx, sy = [], []
        for a in rows:
            sx += [a["mz"], a["mz"], None]
            sy += [-Y_FLOOR, a["rel"], None]
        fig.add_trace(go.Scatter(x=sx, y=sy, mode="lines", name=name, hoverinfo="skip",
                                 line=dict(color=color, width=2)))
    fig.add_trace(go.Scatter(x=xs, y=[a["rel"] for a in ann], mode="markers", name="peaks",
                             marker=dict(size=10, color="rgba(0,0,0,0)"),
                             customdata=list(range(len(ann))), hoverinfo="none"))

    shapes, annotations = [], []
    if theo and x0 <= theo <= x1:
        # dotted marker above the precursor peak (and its m/z label) — or full height if absent
        prec_rel = max((a["rel"] for a in ann if a["kind"] == "precursor"), default=None)
        y_start = -Y_FLOOR if prec_rel is None else prec_rel + 9.0
        if y_start < Y_TOP:
            shapes.append(dict(type="line", xref="x", yref="y", x0=theo, x1=theo, y0=y_start,
                               y1=Y_TOP, layer="below", line=dict(color=muted, width=1, dash="dot")))
        frac = (theo - x0) / (x1 - x0)
        annotations.append(dict(
            x=theo, y=1, xref="x", yref="paper", yanchor="bottom", yshift=4, showarrow=False,
            xanchor="right" if frac > 0.8 else ("left" if frac < 0.2 else "center"),
            text=f"precursor {result.get('adduct', '')} {theo:.4f}".replace("  ", " "),
            font=dict(family=mono, size=10, color=muted)))
    axis = dict(showline=True, linecolor=ink, linewidth=1.5, ticks="outside", tickcolor=ink,
                ticklen=5, tickfont=dict(color=muted, size=11), zeroline=False)
    fig.update_layout(
        height=height, margin=dict(l=62, r=22, t=30, b=50),
        paper_bgcolor="#FFFFFF", plot_bgcolor="#FFFFFF",
        font=dict(family=font, color=ink, size=12),
        showlegend=False, hovermode="x", hoverdistance=16, dragmode="zoom",
        shapes=shapes, annotations=annotations,
        xaxis=dict(axis, title=dict(text="m/z", font=dict(color=muted, size=12)),
                   range=[x0, x1], showgrid=False),
        yaxis=dict(axis, title=dict(text="Relative intensity (%)", font=dict(color=muted, size=12)),
                   range=[-Y_FLOOR, Y_TOP], tickvals=[0, 20, 40, 60, 80, 100],
                   showgrid=True, gridcolor=grid))
    return fig


CHART_TEMPLATE = r'''<!doctype html><html><head><meta charset="utf-8">
<style>
 html,body{margin:0;padding:0;background:#FFFFFF;overflow:hidden;font-family:__FONT__;color:__INK__}
 #plot{width:100%;height:__HEIGHT__px}
 #hl{position:fixed;display:none;width:2px;margin-left:-1px;background:__HL__;pointer-events:none;z-index:5}
 #card{position:fixed;display:none;z-index:10;box-sizing:border-box;width:250px;background:#FFFFFF;
   border:1px solid __LINE__;border-radius:6px;box-shadow:none;padding:9px 11px 10px;
   pointer-events:none}
 #card .k{font-family:__MONO__;font-size:10px;font-weight:600;color:__BLUE__}
 #card .k.off{color:__MUTED__}
 #card .mz{font-family:__MONO__;font-size:17px;font-weight:600;color:__INK__;margin:3px 0 6px}
 #card table{border-collapse:collapse;width:100%}
 #card td{padding:2px 0;vertical-align:baseline}
 #card td.l{font-size:10px;font-weight:500;color:__MUTED__;white-space:nowrap;padding-right:10px}
 #card td.v{font-family:__MONO__;font-size:12px;color:__INK__;text-align:right;word-break:break-all}
 #card img{display:block;width:224px;height:159px;margin:7px auto 0;border-top:1px solid __LINE__;
   padding-top:6px}
</style>
<script>__PLOTLYJS__</script></head>
<body><div id="plot"></div><div id="hl"></div><div id="card"></div>
<script>
(function(){
 var fig = __FIGJSON__, PEAKS = __PEAKSJSON__, CFG = __CFGJSON__;
 var gd = document.getElementById('plot'), card = document.getElementById('card'),
     hl = document.getElementById('hl');
 var STATIC_ANN = (fig.layout.annotations || []).slice(), busy = false;
 var ESC = {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'};
 function esc(s){ return String(s).replace(/[&<>"']/g, function(c){ return ESC[c]; }); }
 function num(v){ return v !== null && v !== undefined && isFinite(v); }
 function fmt(v, d){ return num(v) ? Number(v).toFixed(d) : '—'; }
 function signed(v, d){
   if (!num(v)) return '—';
   var r = Number(Number(v).toFixed(d));          // no "−0.0" for tiny negative errors
   return (r >= 0 ? '+' : '−') + Math.abs(r).toFixed(d);
 }

 function cardHtml(p){
   var kind = p.kind === 'precursor' ? 'Precursor ion · ' + CFG.adduct
            : p.kind === 'fragment' ? 'MetFrag fragment' : 'Not annotated';
   var rows = [];
   if (p.ion) rows.push(['Ion formula', esc(p.ion)]);
   if (num(p.calc)) rows.push(['Theor. m/z', fmt(p.calc, 4)]);
   if (num(p.ppm)) rows.push(['Mass error', signed(p.ppm, 1) + ' ppm']);
   if (p.kind === 'fragment' && num(p.nl)) rows.push(['Neutral loss', fmt(p.nl, 4)]);
   rows.push(['Rel. intensity', fmt(p.rel, 1) + ' %']);
   return '<div class="k' + (p.kind ? '' : ' off') + '">' + esc(kind) + '</div>'
     + '<div class="mz">m/z ' + fmt(p.mz, 4) + '</div><table>'
     + rows.map(function(r){ return '<tr><td class="l">' + r[0] + '</td><td class="v">' + r[1] + '</td></tr>'; }).join('')
     + '</table>' + (p.img ? '<img alt="fragment structure" src="data:image/png;base64,' + p.img + '">' : '');
 }
 function geom(){
   var fl = gd._fullLayout;
   return {r: gd.getBoundingClientRect(), xa: fl.xaxis, ya: fl.yaxis};
 }
 function hide(){ card.style.display = 'none'; hl.style.display = 'none'; }
 function show(p, ev){
   var g = geom(), px = g.r.left + g.xa._offset + g.xa.l2p(p.mz);
   var yTop = g.r.top + g.ya._offset + g.ya.l2p(Math.min(p.rel, g.ya.range[1]));
   var yBot = g.r.top + g.ya._offset + g.ya._length;
   if (px >= g.r.left + g.xa._offset && px <= g.r.left + g.xa._offset + g.xa._length) {
     hl.style.left = px + 'px'; hl.style.top = yTop + 'px';
     hl.style.height = Math.max(3, yBot - yTop) + 'px'; hl.style.display = 'block';
   } else { hl.style.display = 'none'; }
   card.innerHTML = cardHtml(p);
   card.style.display = 'block';
   var e = ev || {}, cx = num(e.clientX) ? e.clientX : px,
       cy = num(e.clientY) ? e.clientY : (yTop + yBot) / 2;
   var bw = card.offsetWidth, bh = card.offsetHeight, vw = window.innerWidth,
       vh = window.innerHeight, m = 6, gap = 18;
   var left = cx + gap;
   if (left + bw > vw - m) left = cx - gap - bw;
   if (left < m) left = Math.max(m, vw - bw - m);
   var top = Math.max(m, Math.min(cy - bh / 2, vh - bh - m));
   card.style.left = left + 'px'; card.style.top = top + 'px';
 }

 // m/z labels: most important peaks first (annotated, then by intensity); a label is only
 // drawn where it neither overlaps another label nor crosses another peak's stick.
 function placeLabels(){
   var g = geom(), xa = g.xa, ya = g.ya;
   var x0 = Math.min(xa.range[0], xa.range[1]), x1 = Math.max(xa.range[0], xa.range[1]);
   var yTop = ya.range[1], vis = [];
   PEAKS.forEach(function(p, i){ if (p.mz >= x0 && p.mz <= x1) vis.push(i); });
   var sticks = vis.map(function(i){
     return {i: i, x: xa.l2p(PEAKS[i].mz), top: ya.l2p(Math.min(PEAKS[i].rel, yTop))};
   });
   if (num(CFG.precursor) && CFG.precursor >= x0 && CFG.precursor <= x1)
     sticks.push({i: -1, x: xa.l2p(CFG.precursor), top: 0});
   var order = vis.slice().sort(function(a, b){
     var pa = PEAKS[a], pb = PEAKS[b];
     return ((pb.kind ? 1 : 0) - (pa.kind ? 1 : 0)) || (pb.rel - pa.rel);
   });
   var boxes = [], anns = [];
   for (var n = 0; n < order.length && anns.length < CFG.maxLabels; n++){
     var i = order[n], p = PEAKS[i];
     if (p.rel > yTop) continue;
     var txt = p.mz.toFixed(4), w = txt.length * 6.1 + 4, px = xa.l2p(p.mz), py = ya.l2p(p.rel);
     var b = {l: px - w / 2, r: px + w / 2, t: py - 16, b: py - 2};
     if (b.l < 0 || b.r > xa._length || b.t < 0) continue;
     var hit = boxes.some(function(o){ return b.l < o.r + 4 && b.r > o.l - 4 && b.t < o.b + 1 && b.b > o.t - 1; })
       || sticks.some(function(s){
            if (s.i === i || (s.i === -1 && p.kind === 'precursor')) return false;
            return s.x > b.l - 1 && s.x < b.r + 1 && s.top < b.b; });
     if (hit) continue;
     boxes.push(b);
     anns.push({x: p.mz, y: p.rel, xref: 'x', yref: 'y', text: txt, showarrow: false,
                xanchor: 'center', yanchor: 'bottom', yshift: 2,
                font: {family: CFG.mono, size: 10, color: p.kind ? CFG.ink : CFG.muted}});
   }
   busy = true;
   Plotly.relayout(gd, {annotations: STATIC_ANN.concat(anns)})
     .then(function(){ busy = false; }, function(){ busy = false; });
 }

 Plotly.newPlot(gd, fig.data, fig.layout, {responsive: true, displayModeBar: false})
   .then(placeLabels);
 gd.on('plotly_hover', function(ev){
   var pt = (ev.points || []).filter(function(q){ return q.data.name === 'peaks'; })[0];
   if (!pt) return;
   show(PEAKS[pt.customdata], ev.event);
 });
 gd.on('plotly_unhover', hide);
 gd.addEventListener('mouseleave', hide);
 gd.on('plotly_relayout', function(e){
   if (busy) return;
   var keys = Object.keys(e || {});
   if (keys.length && keys.every(function(k){ return k.indexOf('annotations') === 0; })) return;
   hide(); placeLabels();
 });
})();
</script></body></html>'''


def _js_json(obj):
    """JSON for inline <script> embedding (no '</script>' breakout)."""
    return json.dumps(obj, allow_nan=False).replace("</", "<\\/")


def _clean(v, nd):
    return None if v is None or not math.isfinite(v) else round(float(v), nd)


@st.cache_data(show_spinner=False, max_entries=3)
def build_chart_html(result, tokens=None, height=CHART_HEIGHT):
    """Self-contained HTML (Plotly + hover card) for components.html.

    ``tokens`` are the design-system colours (keys: ink, muted, line, blue, blue_bright,
    dark, light, grid); defaults match theme.py.
    """
    t = {"ink": "#0A1633", "muted": "#5C6B8A", "line": "#C7D0E4", "blue": "#1B3A8F",
         "blue_bright": "#3B5BDB", "dark": DARK, "light": LIGHT, "grid": "#E9EDF7"}
    t.update(tokens or {})
    fig = build_figure(result, dark=t["dark"], light=t["light"], ink=t["ink"], muted=t["muted"],
                       grid=t["grid"], height=height)
    peaks = [{"mz": a["mz"], "rel": _clean(a["rel"], 3), "kind": a["kind"],
              "ion": a["ion_formula"], "calc": _clean(a["calc_mz"], 6), "ppm": _clean(a["ppm"], 3),
              "nl": _clean(a["nl"], 6), "img": a["b64"]} for a in result["annotation"]]
    cfg = {"adduct": result.get("adduct", ""), "precursor": _clean(result.get("theo_precursor_mz"), 6),
           "ink": t["ink"], "muted": t["muted"], "mono": MONO, "maxLabels": 40}
    html = CHART_TEMPLATE
    for key, val in (("__FONT__", FONT), ("__MONO__", MONO), ("__INK__", t["ink"]),
                     ("__MUTED__", t["muted"]), ("__LINE__", t["line"]), ("__BLUE__", t["blue"]),
                     ("__HL__", t["blue_bright"]), ("__HEIGHT__", str(int(height)))):
        html = html.replace(key, val)
    return (html.replace("__FIGJSON__", fig.to_json().replace("</", "<\\/"))
                .replace("__PEAKSJSON__", _js_json(peaks))
                .replace("__CFGJSON__", _js_json(cfg))
                .replace("__PLOTLYJS__", get_plotlyjs()))
