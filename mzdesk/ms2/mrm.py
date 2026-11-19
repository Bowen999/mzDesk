"""MRM transition selection for a targeted method: ranking, roles, collision energy, CSV.

Pure functions (no I/O, no external tools). The runner supplies the candidate product ions --
merged CFM-ID energies, a MassBank spectrum or a pasted peak list -- and this module shapes them
into the transition table: distinct product m/z, quantifier / qualifier roles, relative
intensities, the eV -> NCE conversion and the method-table CSV.

One transition = one precursor m/z -> product m/z pair.
"""
import csv
import io
import math
import re

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import inchi as rd_inchi

RDLogger.DisableLog("rdApp.*")

# Product ions closer than this (Da) are the same ion; also the distance a product must keep
# from the precursor m/z (an ion at the precursor m/z is the precursor, not a fragment).
TOLERANCE_DA = 0.01
# Generic NCE conversion used by the Calculator page: NCE = CE * 500 / (precursor m/z * z).
NCE_SCALE = 500.0
_MORE = re.compile(r"\s*\(\+\d+ more\)\s*$")


# ---------------------------------------------------------------- fragment annotations
def fragment_smiles(annotation):
    """Fragment SMILES from a CFM-ID peak annotation, or None.

    CFM-ID annotations look like ``C6H8N3O+ · C=[N+](C=O)C1=CN(C)C=N1`` and end in
    ``(+N more)`` when the peak has several explanations (only the first is used).
    An annotation that is already a bare SMILES is returned as it is.
    """
    text = _MORE.sub("", str(annotation or "").strip())
    if not text:
        return None
    if " · " in text:
        text = text.split(" · ", 1)[1].strip()
    return text or None


def fragment_formula(annotation):
    """Molecular formula (RDKit ``CalcMolFormula``, charge included) of the annotated fragment.

    None when the annotation is empty or its SMILES cannot be parsed, so the caller can show a
    dash instead of a wrong formula.
    """
    smiles = fragment_smiles(annotation)
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        return rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        return None


def structure_key(smiles):
    """First 14 characters of the InChIKey (the 2D structure block) of ``smiles``.

    A salt / mixture keeps its largest part, the same rule the engines use. Raises ValueError
    when the SMILES cannot be parsed or no InChIKey can be computed.
    """
    mol = Chem.MolFromSmiles(str(smiles or "").strip())
    if mol is None:
        raise ValueError("Invalid SMILES string. Please check the structure.")
    parts = Chem.GetMolFrags(mol, asMols=True)
    if len(parts) > 1:
        mol = max(parts, key=lambda part: (part.GetNumHeavyAtoms(), Descriptors.ExactMolWt(part)))
    key = rd_inchi.MolToInchiKey(mol)
    if not key:
        raise ValueError("The InChIKey of the structure could not be computed.")
    return key[:14]


# ---------------------------------------------------------------- predicted spectra
def merge_energies(spectra, tolerance=TOLERANCE_DA):
    """Merge per-energy CFM-ID peak lists into one product-ion list.

    ``spectra``: [{"energy", "ev", "peaks": [(m/z, relative intensity)], "fragments": [annotation]}].
    Peaks whose m/z differ by at most ``tolerance`` are one ion. A merged ion keeps the highest
    intensity over the energies, and the m/z, annotation and energy (``ev``) of that occurrence.
    Returns [{"mz", "intensity", "ev", "annotation"}], sorted by m/z.
    """
    spectra = list(spectra or [])
    if not spectra:
        raise ValueError("There is no predicted spectrum to build transitions from.")
    records = []
    for spectrum in spectra:
        peaks = list(spectrum.get("peaks") or [])
        annotations = list(spectrum.get("fragments") or [])
        for index, peak in enumerate(peaks):
            records.append({"mz": float(peak[0]), "intensity": float(peak[1]), "ev": spectrum.get("ev"),
                            "annotation": annotations[index] if index < len(annotations) else None})
    if not records:
        raise ValueError("The predicted spectra have no peaks to build transitions from.")
    records.sort(key=lambda record: (record["mz"], -record["intensity"]))   # stable: equal peaks keep the lower energy
    merged = []
    for record in records:
        if merged and record["mz"] - merged[-1]["anchor"] <= tolerance:
            current = merged[-1]
            if record["intensity"] > current["intensity"]:
                current.update(record)                  # the most intense occurrence defines m/z, annotation and energy
        else:
            merged.append({**record, "anchor": record["mz"]})
    return [{key: ion[key] for key in ("mz", "intensity", "ev", "annotation")} for ion in merged]


# ---------------------------------------------------------------- ranking and roles
def select_transitions(products, precursor_mz, limit, tolerance=TOLERANCE_DA):
    """Rank candidate products and assign the quantifier / qualifier roles.

    ``products``: [{"mz", "intensity", ...}]; extra keys are kept in the result. Peaks at or above
    the precursor m/z are dropped (the precursor ion is not a product), as are duplicate m/z
    values. Products are ranked by descending intensity (ties: lower m/z first) and the most
    intense is the quantifier, the rest are qualifiers. ``relative`` normalises to the base peak
    of ``products`` -- the spectrum the page chose, precursor included -- so the quantifier is
    100 only when it is also that base peak.

    Returns {"transitions": [...], "n_precursor_skipped": int, "n_requested": int}.
    Raises ValueError for a missing / non-positive precursor m/z, a limit below 1, or when no
    product ion is left.
    """
    precursor_mz = _positive(precursor_mz, "precursor m/z")
    try:
        if int(limit) != limit or limit < 1:
            raise ValueError
        limit = int(limit)
    except (TypeError, ValueError):
        raise ValueError("At least one transition is required.") from None
    candidates, intensities, n_precursor = [], [], 0
    for product in products or []:
        try:
            mz, intensity = float(product["mz"]), float(product["intensity"])
        except (KeyError, TypeError, ValueError):
            raise ValueError("Every product ion needs a numeric m/z and intensity.") from None
        if not (math.isfinite(mz) and math.isfinite(intensity)) or mz <= 0 or intensity <= 0:
            continue
        intensities.append(intensity)
        if mz >= precursor_mz - tolerance:
            n_precursor += 1
            continue
        candidates.append({**product, "mz": mz, "intensity": intensity})
    candidates.sort(key=lambda product: (-product["intensity"], product["mz"]))
    selected = []
    for product in candidates:
        if all(abs(product["mz"] - chosen["mz"]) > tolerance for chosen in selected):
            selected.append(product)
            if len(selected) == limit:
                break
    if not selected:
        raise ValueError("No product ion below the precursor m/z. Check the precursor m/z and the spectrum.")
    base = max(intensities)
    transitions = [{**product, "relative": product["intensity"] / base * 100.0,
                    "role": "quantifier" if index == 0 else "qualifier"} for index, product in enumerate(selected)]
    return {"transitions": transitions, "n_precursor_skipped": n_precursor, "n_requested": limit}


# ---------------------------------------------------------------- collision energy
def nce_from_ev(ev, precursor_mz, charge=1):
    """Normalized collision energy from CE in eV: NCE = CE * 500 / (precursor m/z * z), rounded.

    The generic conversion also used by the Calculator page; instruments define NCE differently.
    Returns None for an unknown CE. Raises ValueError for a non-positive precursor m/z or charge.
    """
    if ev is None:
        return None
    try:
        ev, precursor_mz, charge = float(ev), float(precursor_mz), float(charge)
        if not all(math.isfinite(value) for value in (ev, precursor_mz, charge)):
            raise TypeError
    except (TypeError, ValueError):
        raise ValueError("Collision energy, precursor m/z and charge must be numbers.") from None
    if ev <= 0:
        raise ValueError("Collision energy must be positive.")
    if precursor_mz <= 0 or charge <= 0:
        raise ValueError("A positive precursor m/z and charge are needed for the NCE conversion.")
    return int(round(ev * NCE_SCALE / (precursor_mz * charge)))


def ce_label(ev, precursor_mz, unit="eV"):
    """CE cell of a predicted fragment: ``20 eV`` or ``51 NCE``; None when the energy is unknown."""
    if ev is None:
        return None
    if unit == "eV":
        try:
            return f"{float(ev):g} eV"
        except (TypeError, ValueError):
            raise ValueError("Collision energy must be a number.") from None
    if unit == "NCE":
        return f"{nce_from_ev(ev, precursor_mz)} NCE"
    raise ValueError(f"Unknown CE unit {unit!r}; choose eV or NCE.")


# ---------------------------------------------------------------- method table
def method_table_csv(header, rows):
    """CSV text of the transition table (header + rows, quoted when needed, no trailing newline)."""
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(list(header))
    writer.writerows([list(row) for row in rows])
    return buffer.getvalue().rstrip("\n")


def _positive(value, label):
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"A positive {label} is required.") from None
    if not math.isfinite(number) or number <= 0:
        raise ValueError(f"A positive {label} is required.")
    return number
