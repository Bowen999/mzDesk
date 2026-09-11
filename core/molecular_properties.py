"""molecular_properties.py

Basic chemical-identity descriptors and adduct exact-mass helpers for the
Lipidome / Metabolome tab.

Kept independent of the lipid/metabolite heuristic so the interface can show the
same identifiers (InChI, InChI Key, formula) and monoisotopic descriptors next to
the classifier verdict, without coupling the two concerns.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Callable, Dict, List, Optional, Tuple, TypeVar

from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, rdMolDescriptors


T = TypeVar("T")


# (Δ monoisotopic mass, charge) for commonly observed ion forms.
#
# Deltas follow the widely used metabolomics adduct tables (H+, Na+, K+, NH4+,
# water loss, chloride, formate, acetate and the common doubly charged forms).
# Exact m/z = (neutral exact mass + Δ) / |charge|.
ADDUCTS: Dict[str, Tuple[float, int]] = {
    "[M+H]+": (1.007276, 1),
    "[M+Na]+": (22.989218, 1),
    "[M+K]+": (38.963158, 1),
    "[M+NH4]+": (18.033823, 1),
    "[M+H-H2O]+": (-17.003289, 1),
    "[M+2H]2+": (2.014552, 2),
    "[M-H]-": (-1.007276, -1),
    "[M+Cl]-": (34.968852, -1),
    "[M+HCOO]-": (44.997655, -1),
    "[M+CH3COO]-": (59.013305, -1),
    "[M-2H]2-": (-2.014552, -2),
}


@dataclass
class MolecularProperties:
    """Neutral-molecule identity and descriptor set shown in the UI."""
    smiles: str
    valid: bool
    canonical_smiles: Optional[str] = None
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    formula: Optional[str] = None
    exact_mass: Optional[float] = None
    mw: Optional[float] = None
    logp: Optional[float] = None
    tpsa: Optional[float] = None
    hbd: Optional[int] = None
    hba: Optional[int] = None
    rotatable_bonds: Optional[int] = None
    rings: Optional[int] = None
    aromatic_rings: Optional[int] = None
    heavy_atoms: Optional[int] = None


def _try(call: Callable[[], T]) -> Optional[T]:
    """Run a descriptor call and return None if RDKit raises for exotic input."""
    try:
        return call()
    except Exception:
        return None


def from_smiles(smiles: str) -> MolecularProperties:
    """Parse SMILES and compute the basic-property block.

    Returns a MolecularProperties with valid=False when the string is empty or
    cannot be parsed, so callers can branch before rendering values.
    """
    s = (smiles or "").strip()
    if not s:
        return MolecularProperties(smiles=s, valid=False)

    mol = Chem.MolFromSmiles(s)
    if mol is None:
        return MolecularProperties(smiles=s, valid=False)

    canonical = _try(lambda: Chem.MolToSmiles(mol)) or s
    inchi = _try(lambda: Chem.MolToInchi(mol))
    inchikey = _try(lambda: Chem.InchiToInchiKey(inchi)) if inchi else None
    formula = _try(lambda: rdMolDescriptors.CalcMolFormula(mol))

    return MolecularProperties(
        smiles=s,
        valid=True,
        canonical_smiles=canonical,
        inchi=inchi,
        inchikey=inchikey,
        formula=formula,
        exact_mass=_try(lambda: float(Descriptors.ExactMolWt(mol))),
        mw=_try(lambda: float(Descriptors.MolWt(mol))),
        logp=_try(lambda: float(Crippen.MolLogP(mol))),
        tpsa=_try(lambda: float(rdMolDescriptors.CalcTPSA(mol))),
        hbd=_try(lambda: int(rdMolDescriptors.CalcNumHBD(mol))),
        hba=_try(lambda: int(rdMolDescriptors.CalcNumHBA(mol))),
        rotatable_bonds=_try(lambda: int(rdMolDescriptors.CalcNumRotatableBonds(mol))),
        rings=_try(lambda: int(rdMolDescriptors.CalcNumRings(mol))),
        aromatic_rings=_try(lambda: int(rdMolDescriptors.CalcNumAromaticRings(mol))),
        heavy_atoms=_try(lambda: int(mol.GetNumHeavyAtoms())),
    )


def adduct_spec(label: str) -> Tuple[float, int]:
    """Return (Δ mass, charge) for a known adduct label."""
    if label not in ADDUCTS:
        raise KeyError(f"Unknown adduct: {label}")
    return ADDUCTS[label]


def adduct_exact_mz(exact_mass: float, label: str) -> float:
    """Exact m/z of the selected adduct of a neutral exact mass."""
    delta, charge = adduct_spec(label)
    denom = abs(charge) if charge != 0 else 1
    return (exact_mass + delta) / denom


def property_rows(props: MolecularProperties) -> List[Tuple[str, str]]:
    """Format the descriptor grid for `theme.stat_grid`, tolerating None."""

    def _fmt(value: Optional[float], fmt: str) -> str:
        return "N/A" if value is None else fmt.format(value)

    def _int(value: Optional[int], fmt: str = "{:d}") -> str:
        return "N/A" if value is None else fmt.format(value)

    return [
        ("Exact mass (Da)", _fmt(props.exact_mass, "{:.4f}")),
        ("Average MW (Da)", _fmt(props.mw, "{:.4f}")),
        ("logP", _fmt(props.logp, "{:.2f}")),
        ("TPSA (Å²)", _fmt(props.tpsa, "{:.1f}")),
        ("H-bond donors", _int(props.hbd, "{:d}")),
        ("H-bond acceptors", _int(props.hba, "{:d}")),
        ("Rotatable bonds", _int(props.rotatable_bonds, "{:d}")),
        ("Rings", _int(props.rings, "{:d}")),
        ("Aromatic rings", _int(props.aromatic_rings, "{:d}")),
        ("Heavy atoms", _int(props.heavy_atoms, "{:d}")),
    ]


@lru_cache(maxsize=128)
def structure_3d(smiles: str) -> Optional[str]:
    """Generate a 3D conformer and return it as an MDL molblock (with Hs).

    Mirrors MolProbe's ``backend/chem.structure_3d``: AddHs -> ETKDGv3 embed ->
    MMFF optimization (UFF fallback). Returns None when 3D embedding fails so
    the UI can fall back to the 2D depiction.
    """
    mol = Chem.MolFromSmiles(smiles.strip()) if smiles else None
    if mol is None:
        return None

    try:
        m = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xC0FFEE
        if AllChem.EmbedMolecule(m, params) != 0:
            params.useRandomCoords = True
            if AllChem.EmbedMolecule(m, params) != 0:
                return None
        try:
            AllChem.MMFFOptimizeMolecule(m, maxIters=500)
        except Exception:
            try:
                AllChem.UFFOptimizeMolecule(m, maxIters=500)
            except Exception:
                pass
        return Chem.MolToMolBlock(m)
    except Exception:
        return None
