"""molecular_properties.py

Basic chemical-identity descriptors and adduct exact-mass helpers for the
Basic tab.

Kept independent of the lipid/metabolite heuristic so the interface can show the
same identifiers (InChI, InChI Key, formula) and monoisotopic descriptors next to
the classifier verdict, without coupling the two concerns.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Callable, Dict, Optional, Tuple, TypeVar

from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, rdMolDescriptors
from rdkit.Geometry import Point3D


T = TypeVar("T")


# Monoisotopic masses (u): electron and proton from CODATA 2018, nuclides from
# AME2020. The electron mass (0.55 mDa) is visible at the 4th decimal of m/z.
ELECTRON = 0.000548579909
PROTON = 1.007276466621
_H, _C, _N, _O = 1.00782503223, 12.0, 14.00307400443, 15.99491461957
_NA, _K, _CL = 22.9897692820, 38.9637064864, 34.968852682

# (Δ mass, charge) for commonly observed ion forms of a neutral molecule M:
# exact m/z = (M + Δ) / |charge|. Δ is the mass of the ion minus M, so it
# includes the electron: a cation has lost one (Na+ = Na - e), an anion has
# gained one (Cl- = Cl + e), e.g. [M+Cl]- = M + 34.969401, not M + 34.968853.
ADDUCTS: Dict[str, Tuple[float, int]] = {
    "[M+H]+": (PROTON, 1),
    "[M+Na]+": (_NA - ELECTRON, 1),
    "[M+K]+": (_K - ELECTRON, 1),
    "[M+NH4]+": (_N + 4 * _H - ELECTRON, 1),
    "[M+H-H2O]+": (PROTON - (2 * _H + _O), 1),
    "[M+2H]2+": (2 * PROTON, 2),
    "[M-H]-": (-PROTON, -1),
    "[M+Cl]-": (_CL + ELECTRON, -1),
    "[M+HCOO]-": (_H + _C + 2 * _O + ELECTRON, -1),
    "[M+CH3COO]-": (2 * _C + 3 * _H + 2 * _O + ELECTRON, -1),
    "[M-2H]2-": (-2 * PROTON, -2),
}

# Per-attempt time limit (s) for 3D embedding, so a pathological cage such as
# C60 falls back to the 2D depiction instead of blocking the page for ~15 s.
# No ETKDG wall-clock cap: measured, params.timeout does not bound one attempt, it keeps the
# embedder retrying until the cap expires (a C100 alkane went from 7 s to 25 s, PEG50 to 46 s).
# Plain ETKDG already gives up quickly on big flexible molecules (~1-6 s) and the
# random-coordinate retry is the attempt that succeeds for them (C100 4.4 s, PEG50 6.6 s);
# the worst case is a molecule that fails both (C60 fullerene, ~19 s -> 2D depiction).


@dataclass
class MolecularProperties:
    """Identity and descriptor set shown in the UI (for the input as written)."""
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
    formal_charge: int = 0   # net formal charge; != 0 means the input is already an ion
    n_fragments: int = 0     # disconnected components (> 1 for salts / mixtures)


def _try(call: Callable[[], T]) -> Optional[T]:
    """Run a descriptor call and return None if RDKit raises for exotic input."""
    try:
        return call()
    except Exception:
        return None


def parse_error(smiles: str) -> Optional[str]:
    """Explain why RDKit rejects a SMILES, or return None if it parses.

    Separates a syntax error from a chemistry (sanitization) error such as a
    hypervalent atom or an aromatic ring that cannot be kekulized.
    """
    s = (smiles or "").strip()
    if not s:
        return "The SMILES is empty."
    mol = _try(lambda: Chem.MolFromSmiles(s, sanitize=False))
    if mol is None:
        return "Syntax error: check brackets, ring-closure digits and atom symbols."
    if mol.GetNumAtoms() == 0:
        return "The SMILES contains no atoms."
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:  # AtomValenceException, KekulizeException, ...
        detail = " ".join(str(exc).split())
        return f"Chemistry error: {detail}." if detail else "Chemistry error."
    return None


def from_smiles(smiles: str) -> MolecularProperties:
    """Parse SMILES and compute the basic-property block.

    Returns a MolecularProperties with valid=False when the string is empty or
    cannot be parsed, so callers can branch before rendering values.
    """
    s = (smiles or "").strip()
    if not s:
        return MolecularProperties(smiles=s, valid=False)

    mol = _try(lambda: Chem.MolFromSmiles(s))
    if mol is None or mol.GetNumAtoms() == 0:
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
        formal_charge=_try(lambda: int(Chem.GetFormalCharge(mol))) or 0,
        n_fragments=_try(lambda: len(Chem.GetMolFrags(mol))) or 1,
    )


def adduct_spec(label: str) -> Tuple[float, int]:
    """Return (Δ mass, charge) for a known adduct label."""
    if label not in ADDUCTS:
        raise KeyError(f"Unknown adduct: {label}")
    return ADDUCTS[label]


def adduct_exact_mz(exact_mass: Optional[float], label: str) -> Optional[float]:
    """Exact m/z of an adduct of a neutral exact mass (None if the mass is unknown)."""
    delta, charge = adduct_spec(label)
    if exact_mass is None:
        return None
    return (exact_mass + delta) / (abs(charge) or 1)


def _spread_fragments(mol: Chem.Mol, gap: float = 3.0) -> None:
    """Line disconnected fragments up along x, `gap` Å apart (in place).

    ETKDG embeds every fragment of a salt or mixture around the origin, so
    without this the ions / molecules would be drawn on top of each other.
    """
    frags = Chem.GetMolFrags(mol)
    if len(frags) < 2:
        return
    conf = mol.GetConformer()
    edge = None
    for frag in frags:
        pts = [conf.GetAtomPosition(i) for i in frag]
        cx = sum(p.x for p in pts) / len(pts)
        cy = sum(p.y for p in pts) / len(pts)
        cz = sum(p.z for p in pts) / len(pts)
        shift = (0.0 if edge is None else edge + gap) - (min(p.x for p in pts) - cx)
        for i, p in zip(frag, pts):
            conf.SetAtomPosition(i, Point3D(p.x - cx + shift, p.y - cy, p.z - cz))
        edge = max(p.x for p in pts) - cx + shift


@lru_cache(maxsize=128)
def structure_3d(smiles: str) -> Optional[str]:
    """Generate a 3D conformer and return it as an MDL molblock (with Hs).

    Mirrors MolProbe's ``backend/chem.structure_3d``: AddHs -> ETKDGv3 embed ->
    MMFF optimization (UFF fallback). Fragments of salts / mixtures are laid
    out side by side. Returns None when 3D embedding fails so the UI can fall
    back to the 2D depiction.
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
        _spread_fragments(m)
        # fewer iterations for big molecules: 500 iters cost ~8 s on a 350-atom PEG chain
        iters = 500 if m.GetNumAtoms() <= 200 else 200
        try:
            AllChem.MMFFOptimizeMolecule(m, maxIters=iters)
        except Exception:
            try:
                AllChem.UFFOptimizeMolecule(m, maxIters=iters)
            except Exception:
                pass
        return Chem.MolToMolBlock(m)
    except Exception:
        return None
