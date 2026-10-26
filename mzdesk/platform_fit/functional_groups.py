"""Transparent SMARTS-based functional-group search for the Properties page.

This is a curated library of common small-molecule groups, not a claim that
chemistry has one universally accepted, exhaustive functional-group ontology.
Specific and generic matches may overlap (for example, lactone and ester).
Atom indices are zero-based here; the API formats them as one-based numbers.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List

from rdkit import Chem


@dataclass(frozen=True)
class GroupMatch:
    name: str
    smarts: str
    matches: tuple

    @property
    def count(self) -> int:
        return len(self.matches)


# Ordered roughly by chemical family. Keep names stable: they are displayed in
# the result table and used by tests. Each pattern is intentionally visible in
# the UI so users can audit exactly what was counted.
RULES = (
    ("Carboxylic acid", "[CX3](=[OX1])[OX2H1,OX1-]"),
    ("Ester", "[CX3;!$(C-[NX3]);!$(C-[OX2]-[CX3]=[OX1])](=[OX1])[OX2H0][#6]"),
    ("Lactone", "[CX3;R](=[OX1])[OX2H0;R][#6;R]"),
    ("Amide", "[NX3][CX3](=[OX1])[#6]"),
    ("Imide", "[CX3](=[OX1])[NX3][CX3]=[OX1]"),
    ("Carbamate", "[NX3][CX3](=[OX1])[OX2]"),
    ("Urea", "[NX3][CX3](=[OX1])[NX3]"),
    ("Carbonate", "[OX2][CX3](=[OX1])[OX2]"),
    ("Acyl halide", "[CX3](=[OX1])[F,Cl,Br,I]"),
    ("Carboxylic anhydride", "[CX3](=[OX1])[OX2][CX3]=[OX1]"),
    ("Aldehyde", "[CX3H1](=[OX1])[#6]"),
    ("Ketone", "[#6][CX3](=[OX1])[#6]"),
    ("Thioester", "[CX3](=[OX1])[SX2H0][#6]"),
    ("Thioamide", "[NX3][CX3](=[SX1])[#6]"),
    ("Alcohol", "[OX2H1][CX4;!$(C=O)]"),
    ("Phenol", "[OX2H1][c]"),
    ("Enol", "[OX2H1][CX3]=[CX3]"),
    ("Ether", "[OD2H0;!$(O-[CX3]=[OX1])]([#6])[#6]"),
    ("Epoxide", "[O;r3]1[C;r3][C;r3]1"),
    ("Acetal / ketal", "[OX2H0][CX4]([OX2H0])"),
    ("Hemiacetal / hemiketal", "[OX2H1][CX4][OX2H0]"),
    ("Peroxide", "[OX2H0][OX2H0]"),
    ("Hydroperoxide", "[OX2H1][OX2]"),
    ("Thiol", "[SX2H1]"),
    ("Thioether", "[SX2H0]([#6])[#6]"),
    ("Disulfide", "[SX2H0][SX2H0]"),
    ("Sulfoxide", "[#6][SX3](=[OX1])[#6]"),
    ("Sulfone", "[#6][SX4](=[OX1])(=[OX1])[#6]"),
    ("Sulfonic acid", "[SX4](=[OX1])(=[OX1])[OX2H1,OX1-]"),
    ("Sulfonamide", "[NX3][SX4](=[OX1])=[OX1]"),
    ("Primary amine", "[NX3;H2;D1;!$(N-[CX3]=[OX1]);!$(N-[SX4]=[OX1])][#6]"),
    ("Secondary amine", "[NX3;H1;D2;!$(N-[CX3]=[OX1]);!$(N-[SX4]=[OX1])]([#6])[#6]"),
    ("Tertiary amine", "[NX3;H0;D3;!$(N-[CX3]=[OX1]);!$(N-[SX4]=[OX1])]([#6])([#6])[#6]"),
    ("Quaternary ammonium", "[NX4+;D4]([#6])([#6])([#6])[#6]"),
    ("Aromatic NH", "[nH]"),
    ("Pyridine-type aromatic N", "[nX2;H0]"),
    ("Amine N-oxide", "[NX4+;D4]([OX1-])([#6])([#6])[#6]"),
    ("Imine", "[NX2]=[CX3]"),
    ("Oxime", "[OX2H1][NX2]=[CX3]"),
    ("Hydrazone", "[NX3][NX2]=[CX3]"),
    ("Hydrazine", "[NX3][NX3]"),
    ("Azo", "[#6][NX2]=[NX2][#6]"),
    ("Azide", "[#6][NX2]=[NX2+]=[NX1-]"),
    ("Enamine", "[NX3;!$(N-[CX3]=[OX1])][CX3]=[CX3]"),
    ("Nitrile", "[CX2]#[NX1]"),
    ("Isonitrile", "[#6][NX1+]#[CX1-]"),
    ("Isocyanate", "[#6][NX2]=[CX2]=[OX1]"),
    ("Isothiocyanate", "[#6][NX2]=[CX2]=[SX1]"),
    ("Cyanate", "[OX2][CX2]#[NX1]"),
    ("Thiocyanate", "[SX2][CX2]#[NX1]"),
    ("Nitro", "[N;H0;$(N-[#6]);D3](=[O;D1])~[O;D1]"),
    ("Nitroso", "[#6][NX2]=[OX1]"),
    ("Amidine", "[NX3][CX3](=[NX2])[#6]"),
    ("Guanidine", "[NX3][CX3](=[NX2,NX2+])([NX3])"),
    ("Alkene", "[CX3]=[CX3]"),
    ("Alkyne", "[CX2]#[CX2]"),
    ("Allene", "[CX3]=[CX2]=[CX3]"),
    ("Alkyl halide", "[F,Cl,Br,I][CX4]"),
    ("Aryl halide", "[F,Cl,Br,I][c]"),
    ("Phosphate / phosphonate", "[PX4](=[OX1])([OX2,OX1-])([OX2,OX1-])[#6,OX2,OX1-]"),
    ("Boronic acid", "[#6][BX3]([OX2H1])[OX2H1]"),
    ("Boronate ester", "[#6][BX3]([OX2H0])[OX2H0]"),
    ("Silanol", "[Si][OX2H1]"),
    ("Silyl ether", "[Si][OX2H0][#6]"),
)

_PATTERNS = tuple((name, smarts, Chem.MolFromSmarts(smarts)) for name, smarts in RULES)
assert all(pattern is not None for _, _, pattern in _PATTERNS)
assert len({name for name, _ in RULES}) == len(RULES)


def find_functional_groups(smiles: str) -> List[GroupMatch]:
    """Return every matching rule, in RULES order, with unique atom tuples.

    Raises ValueError for an empty or invalid SMILES. Different rules may cover
    the same atoms; within one rule RDKit's uniquify removes symmetry duplicates.
    """
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError("RDKit could not parse the SMILES.")
    found = []
    for name, smarts, pattern in _PATTERNS:
        matches = tuple(mol.GetSubstructMatches(pattern, uniquify=True))
        if matches:
            found.append(GroupMatch(name=name, smarts=smarts, matches=matches))
    return found


def atom_text(matches: tuple, *, shown: int = 8) -> str:
    """One-based atom indices; semicolons separate individual group matches."""
    parts = [", ".join(str(index + 1) for index in match) for match in matches[:shown]]
    if len(matches) > shown:
        parts.append("…")
    return "; ".join(parts)
