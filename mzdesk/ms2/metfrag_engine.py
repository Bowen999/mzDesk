"""MetFrag engine for the MS2 fragment-explanation page.

Runs the bundled MetFrag CLI (mzdesk/lib/MetFragCommandLine-2.5.2.jar) on ONE user-supplied
candidate structure and maps MetFrag's explained peaks back onto the user's peak list
(ion formula, theoretical m/z, mass error, neutral loss, fragment structure).

The shared peak-list parser ``parse_spectrum`` lives in mzdesk/shared/spectrum.py and is
re-exported here. Also builds the stick-spectrum Plotly figure (on mzdesk/shared/ms2_plot.py); the web
page adds the hover card with the fragment structure (web/src/components/Figure.tsx).
"""
import base64
import math
import re
import subprocess
import tempfile
from functools import lru_cache
from pathlib import Path

import plotly.graph_objects as go
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.MolStandardize import rdMolStandardize

from mzdesk.shared import ms2_plot
from mzdesk.shared.spectrum import parse_spectrum  # noqa: F401  re-exported for callers of the engine

RDLogger.DisableLog('rdApp.*')

JAR = Path(__file__).parent.parent / "lib" / "MetFragCommandLine-2.5.2.jar"
METFRAG_TIMEOUT_S = 1200


# --------------------------------------------------------------------------
# Masses, formulas, adducts
# --------------------------------------------------------------------------
# Monoisotopic (most abundant isotope) masses.
ATOMIC = {'C': 12.0, 'H': 1.00782503, 'O': 15.99491462, 'P': 30.97376199, 'N': 14.0030740,
          'S': 31.9720707, 'F': 18.99840322, 'Cl': 34.96885268, 'Br': 78.9183376,
          'I': 126.904473, 'Na': 22.98976928, 'K': 38.96370649, 'Si': 27.9769265,
          'B': 11.0093054, 'Se': 79.9165218, 'D': 2.01410178}
ELECTRON = 0.000548579909

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


def precursor_info(smiles, adduct):
    """Neutral formula / theoretical precursor m/z of a structure for an adduct (no MetFrag run).

    Returns {"error": str|None, "formula", "mono_mass", "theo_mz", "ion_formula"}.
    """
    return dict(_precursor_info(smiles, adduct))


@lru_cache(maxsize=32)
def _precursor_info(smiles, adduct):
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
# Stick spectrum (Plotly)
# --------------------------------------------------------------------------
Y_FLOOR = 1.6      # sticks start slightly below 0 % so even tiny peaks stay visible above the axis line
Y_TOP = 116.0      # headroom for the m/z labels above the base peak


def build_figure(result, height=ms2_plot.HEIGHT):
    """Stick spectrum: dark = annotated, light = unannotated, dotted line = theoretical precursor.

    Built from the shared MS2 plot helpers (mzdesk/shared/ms2_plot.py); an invisible marker per peak
    (trace 'peaks', customdata = peak index) drives the hover card in hovermode 'x'.
    """
    ann = result["annotation"]
    fig = go.Figure()
    if not ann:
        return fig
    theo = result.get("theo_precursor_mz") or result.get("precursor_mz")
    xs = [a["mz"] for a in ann]
    x_range = ms2_plot.mz_range(xs + ([theo] if theo else []))

    for name, color, rows in (("unannotated", ms2_plot.LIGHT, [a for a in ann if not a["annotated"]]),
                              ("annotated", ms2_plot.DARK, [a for a in ann if a["annotated"]])):
        if rows:
            fig.add_trace(ms2_plot.stick_trace([a["mz"] for a in rows], [a["rel"] for a in rows], color, name, y0=-Y_FLOOR))
    fig.add_trace(ms2_plot.hover_markers(xs, [a["rel"] for a in ann]))

    shapes, annotations = [], []
    if theo and x_range[0] <= theo <= x_range[1]:
        # dotted marker above the precursor peak (and its m/z label) — or full height if absent
        prec_rel = max((a["rel"] for a in ann if a["kind"] == "precursor"), default=None)
        y_start = -Y_FLOOR if prec_rel is None else prec_rel + 9.0
        if y_start < Y_TOP:
            shapes.append(ms2_plot.precursor_line(theo, y_start, Y_TOP))
        annotations.append(ms2_plot.precursor_label(
            theo, f"precursor {result.get('adduct', '')} {theo:.4f}".replace("  ", " "), x_range))
    fig.update_layout(**ms2_plot.base_layout(
        height, x_range=x_range, y_range=[-Y_FLOOR, Y_TOP], y_title="Relative intensity (%)",
        shapes=shapes, annotations=annotations, yaxis=dict(tickvals=[0, 20, 40, 60, 80, 100])))
    return fig


def _clean(v, nd):
    return None if v is None or not math.isfinite(v) else round(float(v), nd)
