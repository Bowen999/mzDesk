"""Runner 共用的输入解析与输出格式。

输入：前端按 spec 字段发来的 JSON（复合字段 structure / precursor / spectrum 是对象，数字可能是字符串）。
输出：{块 id: 数据}。数据形状：
  table   {"rows": [[cell, ...], ...]}；cell 是字符串，或 {"img": data-URL}，或 {"text": str, "href": URL}（链接）
  cards   {"values": {label: 字符串}}
  figure  {"plotly": {data, layout}, "peaks": [...]?} | {"svg": "..."} | {"mol3d": molblock}
  code    {"text": "..."}
  formats {"formats": {format label: text}}
  notes   {"items": ["...", ...]}
"""
import json
import math
import re

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

from mzdesk.shared.spectrum import first_record, parse_spectrum

RDLogger.DisableLog("rdApp.*")

MISSING = "—"
_FORMULA = re.compile(r"^(?:[A-Z][a-z]?\d*)+$")


class InputError(ValueError):
    """用户输入不对 → HTTP 400，消息直接给用户看。"""


class EngineError(RuntimeError):
    """外部工具运行失败 → HTTP 502；log 是工具输出的末尾。"""

    def __init__(self, message, log=""):
        super().__init__(message)
        self.log = log


# ---------------------------------------------------------------- inputs
def text(inputs, key):
    value = inputs.get(key)
    return "" if value is None else str(value).strip()


def number(inputs, key, label, *, default=None, required=False, positive=False, minimum=None, maximum=None):
    """Optional float field; blank -> default. Raises InputError with the field label."""
    raw = inputs.get(key)
    if raw is None or (isinstance(raw, str) and not raw.strip()):
        if required:
            raise InputError(f"{label} is required.")
        return default
    try:
        value = float(raw)
    except (TypeError, ValueError):
        raise InputError(f"{label} must be a number.") from None
    if not math.isfinite(value):
        raise InputError(f"{label} must be a finite number.")
    if positive and value <= 0:
        raise InputError(f"{label} must be positive.")
    if minimum is not None and value < minimum:
        raise InputError(f"{label} must be at least {minimum:g}.")
    if maximum is not None and value > maximum:
        raise InputError(f"{label} must be at most {maximum:g}.")
    return value


def is_formula(s):
    """True for a plain element formula such as C8H10N4O2 (every symbol a real element)."""
    if not s or not _FORMULA.match(s):
        return False
    table = Chem.GetPeriodicTable()
    for el in re.findall(r"[A-Z][a-z]?", s):
        try:
            if table.GetAtomicNumber(el) <= 0:
                return False
        except Exception:
            return False
    return True


def structure(inputs, key="structure", label="Structure", *, required=True, allow_formula=False):
    """Structure field -> {"kind": "smiles", "smiles": ...} | {"kind": "formula", "formula": ...} | None.

    SMILES / InChI mode accepts SMILES or an InChI string. Name / ID / Formula mode accepts a formula
    where the page allows one; anything else is resolved to a SMILES by mzdesk.identity.resolve.
    """
    value = inputs.get(key) or {}
    mode = value.get("mode") or "SMILES / InChI"
    s = str(value.get("value") or "").strip()
    if mode == "Draw" and not s:
        raise InputError("The structure editor is not connected yet. Use SMILES / InChI.")
    if not s:
        if required:
            raise InputError(f"{label} is required.")
        return None
    if mode == "Name / ID / Formula":
        from mzdesk.identity import resolve
        # C00001 之类的 KEGG ID 长得像分子式：只有识别为 Name 才当分子式
        if allow_formula and is_formula(s) and resolve.detect(s) == "Name":
            return {"kind": "formula", "formula": s}
        try:
            return {"kind": "smiles", "smiles": resolve.to_smiles(s)}
        except resolve.ResolveError as exc:
            raise InputError(f"{label}: {exc}") from None
    if s.startswith("InChI="):
        mol = Chem.MolFromInchi(s)
        if mol is None:
            raise InputError(f"{label}: the InChI could not be parsed.")
        return {"kind": "smiles", "smiles": Chem.MolToSmiles(mol)}
    from mzdesk.platform_fit.molecular_properties import parse_error
    err = parse_error(s)
    if err:
        hint = " For a formula, choose Name / ID / Formula." if allow_formula and is_formula(s) else ""
        raise InputError(f"{label}: invalid SMILES. {err}{hint}")
    return {"kind": "smiles", "smiles": s}


def record_name(inputs, smiles, key="structure"):
    """Compound name for exported spectrum records: what the user typed in Name / ID mode, else the molecular formula of ``smiles``."""
    value = inputs.get(key) or {}
    typed = str(value.get("value") or "").strip()
    if value.get("mode") == "Name / ID / Formula" and typed:
        return typed
    return rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smiles))


def precursor(inputs, key="precursor", label="Precursor m/z", *, required=False):
    """Precursor field -> (m/z or None, positive: bool, adduct or None for Auto)."""
    value = inputs.get(key) or {}
    mz = number(value, "precursor_mz", label, required=required, positive=True)
    positive = value.get("ion_mode", "Positive") != "Negative"
    adduct = value.get("adduct") or "Auto"
    return mz, positive, (None if adduct == "Auto" else adduct)


def spectrum(inputs, key, label, warnings, *, required=True):
    """Spectrum field -> (peaks [(mz, intensity)], precursor m/z read from the text or None).

    The text may be any format parse_spectrum reads; MGF / MSP / MassBank / SIRIUS headers give the
    precursor (and only the first spectrum of an MGF / MSP file is used).
    """
    value = inputs.get(key) or {}
    source = value.get("source") or "Paste"
    raw = str(value.get("value") or "")
    if source == "USI":
        if raw.strip():
            raise InputError(f"{label}: USI lookup is not connected yet. Paste the peak list or upload an MGF / MSP file.")
        raw = ""
    file_precursor = None
    if raw.strip():
        rec = first_record(raw)
        if rec["n_records"] > 1:
            warnings.append(f"{label}: the input holds {rec['n_records']} spectra; only the first"
                            + (f" ({rec['name']})" if rec["name"] else "") + " is used.")
        raw, file_precursor = rec["peaks_text"], rec["precursor_mz"]
    peaks = parse_spectrum(raw)
    if not peaks and (required or raw.strip()):
        raise InputError(f"{label}: no valid peaks. Use one 'm/z intensity' pair per line.")
    return peaks, file_precursor


# ---------------------------------------------------------------- outputs
def fmt(value, nd=4):
    if value is None or not isinstance(value, (int, float)) or not math.isfinite(value):
        return MISSING
    return f"{value:.{nd}f}"


def signed(value, nd=1):
    if value is None or not math.isfinite(value):
        return MISSING
    r = round(float(value), nd)
    return f"{'+' if r >= 0 else '−'}{abs(r):.{nd}f}"


def plotly(fig, **extra):
    return {"plotly": json.loads(fig.to_json()), **extra}


def link(text, href):
    """Link cell for a table; MISSING when there is no text."""
    return {"text": str(text), "href": href} if text else MISSING


def png(b64):
    return {"img": f"data:image/png;base64,{b64}"} if b64 else ""


def depiction_svg(smiles, size=(420, 300)):
    """Plain 2D depiction as SVG (None if RDKit cannot parse it)."""
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return None
    drawer = rdMolDraw2D.MolDraw2DSVG(*size)
    drawer.drawOptions().clearBackground = False
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def clean(obj):
    """JSON-safe copy: NaN / inf -> None, numpy scalars -> Python, tuples -> lists."""
    if isinstance(obj, dict):
        return {str(k): clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean(v) for v in obj]
    if hasattr(obj, "item") and not isinstance(obj, (str, bytes)):
        obj = obj.item()
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    return obj
