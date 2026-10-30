"""Transparent SMARTS-based functional-group search for the Properties page.

This is a curated library of common small-molecule groups, ring systems and biomolecule motifs, not a claim that
chemistry has one universally accepted, exhaustive functional-group ontology.
Specific and generic matches may overlap (for example, lactone and ester, or primary alcohol and alcohol).
Atom indices are zero-based here; the API formats them as one-based numbers.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import List

from rdkit import Chem

from mzdesk.shared import colors as theme, depict


@dataclass(frozen=True)
class GroupMatch:
    family: str
    name: str
    smarts: str
    matches: tuple

    @property
    def count(self) -> int:
        return len(self.matches)


# Each rule: (family, name, SMARTS, drawing SMILES). Families are listed in display order. Keep names stable: they are
# shown in the result and used by tests. The SMARTS is shown in the UI so users can audit exactly what was counted.
# The drawing SMILES is the group as a chemist would sketch it ("*" is drawn as R); a test checks that every SMARTS
# matches its own drawing (with "*" read as a methyl group), so the picture and the rule cannot drift apart.
# (family, name, SMARTS, drawing SMILES; "*" is drawn as R)
RULES = (
    # ---- Carbonyl
    ("Carbonyl", "Carboxylic acid", "[CX3](=[OX1])[OX2H1,OX1-]", "*C(=O)O"),
    ("Carbonyl", "Carboxylate anion", "[CX3](=[OX1])[OX1-]", "*C(=O)[O-]"),
    ("Carbonyl", "Ester", "[CX3;!$(C-[NX3]);!$(C-[OX2]-[CX3]=[OX1])](=[OX1])[OX2H0][#6]", "*C(=O)O*"),
    ("Carbonyl", "Lactone", "[CX3;R](=[OX1])[OX2H0;R][#6;R]", "O=C1CCCO1"),
    ("Carbonyl", "Amide", "[NX3][CX3](=[OX1])[#6]", "*C(=O)N"),
    ("Carbonyl", "Lactam", "[NX3]@[CX3;R](=[OX1])[#6]", "O=C1CCCN1"),
    ("Carbonyl", "β-Lactam", "[CX3;r4](=[OX1])[NX3;r4]", "O=C1CCN1"),
    ("Carbonyl", "Imide", "[CX3](=[OX1])[NX3][CX3]=[OX1]", "*C(=O)NC(=O)*"),
    ("Carbonyl", "Carbamate", "[NX3][CX3](=[OX1])[OX2]", "*OC(=O)N"),
    ("Carbonyl", "Urea", "[NX3][CX3](=[OX1])[NX3]", "NC(=O)N"),
    ("Carbonyl", "Carbonate", "[OX2][CX3](=[OX1])[OX2]", "*OC(=O)O*"),
    ("Carbonyl", "Acyl halide", "[CX3](=[OX1])[F,Cl,Br,I]", "*C(=O)Cl"),
    ("Carbonyl", "Carboxylic anhydride", "[CX3](=[OX1])[OX2][CX3]=[OX1]", "*C(=O)OC(=O)*"),
    ("Carbonyl", "Aldehyde", "[CX3H1](=[OX1])[#6]", "*C=O"),
    ("Carbonyl", "Ketone", "[#6][CX3](=[OX1])[#6]", "*C(=O)*"),
    ("Carbonyl", "1,3-Dicarbonyl", "[CX3](=[OX1])[CX4,CX3][CX3](=[OX1])", "*C(=O)CC(=O)*"),
    ("Carbonyl", "α,β-Unsaturated carbonyl", "[CX3]=[CX3][CX3]=[OX1]", "*C=CC(=O)*"),
    ("Carbonyl", "Heteroaromatic C=O (pyridone-type)", "[c]=[OX1]", "O=c1cccc[nH]1"),
    ("Carbonyl", "Quinone", "[#6X3]1(=[OX1])[#6X3]=,:[#6X3][#6X3](=[OX1])[#6X3]=,:[#6X3]1", "O=C1C=CC(=O)C=C1"),
    ("Carbonyl", "Thioester", "[CX3](=[OX1])[SX2H0][#6]", "*C(=O)S*"),
    ("Carbonyl", "Thioamide", "[NX3][CX3](=[SX1])[#6]", "*C(=S)N"),
    ("Carbonyl", "Thiourea", "[NX3][CX3](=[SX1])[NX3]", "NC(=S)N"),
    ("Carbonyl", "Dithiocarbamate", "[NX3][CX3](=[SX1])[SX2]", "NC(=S)S*"),
    ("Carbonyl", "Thioketone", "[#6][CX3](=[SX1])[#6]", "*C(=S)*"),
    ("Carbonyl", "Hydroxamic acid", "[CX3](=[OX1])[NX3][OX2H1]", "*C(=O)NO"),
    ("Carbonyl", "Hydrazide", "[CX3](=[OX1])[NX3][NX3]", "*C(=O)NN"),
    ("Carbonyl", "Sulfonylurea", "[NX3][CX3](=[OX1])[NX3][SX4](=[OX1])=[OX1]", "NC(=O)NS(=O)(=O)*"),
    # ---- Oxygen
    ("Oxygen", "Alcohol", "[OX2H1][CX4;!$(C=O)]", "*CO"),
    ("Oxygen", "Primary alcohol", "[OX2H1][CH2X4]", "*CO"),
    ("Oxygen", "Secondary alcohol", "[OX2H1][CH1X4]", "*C(*)O"),
    ("Oxygen", "Tertiary alcohol", "[OX2H1][CX4H0]", "*C(*)(*)O"),
    ("Oxygen", "Vicinal diol", "[OX2H1][CX4][CX4][OX2H1]", "*C(O)C(O)*"),
    ("Oxygen", "Phenol", "[OX2H1][c]", "Oc1ccccc1"),
    ("Oxygen", "Catechol", "[OX2H1]c:c[OX2H1]", "Oc1ccccc1O"),
    ("Oxygen", "Hydroquinone", "[OX2H1]c1ccc([OX2H1])cc1", "Oc1ccc(O)cc1"),
    ("Oxygen", "Enol", "[OX2H1][CX3]=[CX3]", "*C(O)=C*"),
    ("Oxygen", "Ether", "[OD2H0;!$(O-[CX3]=[OX1])]([#6])[#6]", "*O*"),
    ("Oxygen", "Aryl ether", "[OD2H0;!$(O-[CX3]=[OX1])](c)[#6]", "*Oc1ccccc1"),
    ("Oxygen", "Methoxy", "[CH3][OD2H0][#6]", "*OC"),
    ("Oxygen", "Methylenedioxy", "c[OX2][CH2][OX2]c", "c1ccc2OCOc2c1"),
    ("Oxygen", "Epoxide", "[O;r3]1[C;r3][C;r3]1", "C1CO1"),
    ("Oxygen", "Acetal / ketal", "[OX2H0][CX4]([OX2H0])", "*OC(*)O*"),
    ("Oxygen", "Hemiacetal / hemiketal", "[OX2H1][CX4][OX2H0]", "*OC(*)O"),
    ("Oxygen", "Peroxide", "[OX2H0][OX2H0]", "*OO*"),
    ("Oxygen", "Hydroperoxide", "[OX2H1][OX2]", "*OO"),
    # ---- Nitrogen
    ("Nitrogen", "Primary amine", "[NX3;H2;D1;!$(N-[CX3]=[OX1]);!$(N-[SX4]=[OX1])][#6]", "*N"),
    ("Nitrogen", "Secondary amine", "[NX3;H1;D2;!$(N-[CX3]=[OX1]);!$(N-[SX4]=[OX1])]([#6])[#6]", "*N*"),
    ("Nitrogen", "Tertiary amine", "[NX3;H0;D3;!$(N-[CX3]=[OX1]);!$(N-[SX4]=[OX1])]([#6])([#6])[#6]", "*N(*)*"),
    ("Nitrogen", "Aryl amine (aniline-type)", "[NX3;!$(N-[CX3]=[OX1,SX1]);!$(N-[SX4]=[OX1]);!$(N-[N,O]);!$(N=*)]-c", "Nc1ccccc1"),
    ("Nitrogen", "Protonated amine", "[NX4+;H1,H2,H3]", "*[NH3+]"),
    ("Nitrogen", "Quaternary ammonium", "[NX4+;D4]([#6])([#6])([#6])[#6]", "*[N+](*)(*)*"),
    ("Nitrogen", "Aromatic NH", "[nH]", "c1cc[nH]c1"),
    ("Nitrogen", "Pyridine-type aromatic N", "[nX2;H0]", "c1ccncc1"),
    ("Nitrogen", "Aromatic N⁺ (pyridinium-type)", "[n+;!$([n+][O-])]", "C[n+]1ccccc1"),
    ("Nitrogen", "Aromatic N-oxide", "[n+][O-]", "[O-][n+]1ccccc1"),
    ("Nitrogen", "Amine N-oxide", "[NX4+;D4]([OX1-])([#6])([#6])[#6]", "*[N+]([O-])(*)*"),
    ("Nitrogen", "Aziridine", "[NX3;r3]", "C1CN1"),
    ("Nitrogen", "Imine", "[NX2]=[CX3]", "*C=N"),
    ("Nitrogen", "Oxime", "[OX2H1][NX2]=[CX3]", "*C=NO"),
    ("Nitrogen", "Oxime ether", "[CX3]=[NX2][OX2H0]", "*C=NO*"),
    ("Nitrogen", "Hydroxylamine", "[NX3;!$(N-[CX3]=[OX1])][OX2H1]", "*NO"),
    ("Nitrogen", "Hydrazone", "[NX3][NX2]=[CX3]", "*NN=C*"),
    ("Nitrogen", "Hydrazine", "[NX3][NX3]", "*NN*"),
    ("Nitrogen", "Azo", "[#6][NX2]=[NX2][#6]", "*N=N*"),
    ("Nitrogen", "Azide", "[#6][NX2]=[NX2+]=[NX1-]", "*N=[N+]=[N-]"),
    ("Nitrogen", "Diazo", "[#6]=[N+]=[N-]", "*C=[N+]=[N-]"),
    ("Nitrogen", "Diazonium", "[#6][N+]#[N]", "*[N+]#N"),
    ("Nitrogen", "Enamine", "[NX3;!$(N-[CX3]=[OX1])][CX3]=[CX3]", "*NC=C*"),
    ("Nitrogen", "Nitrile", "[CX2]#[NX1]", "*C#N"),
    ("Nitrogen", "Isonitrile", "[#6][NX2+]#[CX1-]", "*[N+]#[C-]"),
    ("Nitrogen", "Isocyanate", "[#6][NX2]=[CX2]=[OX1]", "*N=C=O"),
    ("Nitrogen", "Isothiocyanate", "[#6][NX2]=[CX2]=[SX1]", "*N=C=S"),
    ("Nitrogen", "Cyanate", "[OX2][CX2]#[NX1]", "*OC#N"),
    ("Nitrogen", "Thiocyanate", "[SX2][CX2]#[NX1]", "*SC#N"),
    ("Nitrogen", "Cyanamide", "[NX3][CX2]#[NX1]", "*NC#N"),
    ("Nitrogen", "Carbodiimide", "[NX2]=[CX2]=[NX2]", "*N=C=N*"),
    ("Nitrogen", "Nitro", "[N;H0;$(N-[#6]);D3](=[O;D1])~[O;D1]", "*[N+](=O)[O-]"),
    ("Nitrogen", "Nitroso", "[#6][NX2]=[OX1]", "*N=O"),
    ("Nitrogen", "Nitrate ester", "[OX2][N+](=[OX1])[OX1-]", "*O[N+](=O)[O-]"),
    ("Nitrogen", "Nitrite ester", "[OX2][NX2]=[OX1]", "*ON=O"),
    ("Nitrogen", "Amidine", "[NX3][CX3](=[NX2])[#6]", "*C(=N)N"),
    ("Nitrogen", "Guanidine", "[NX3][CX3](=[NX2,NX2+])([NX3])", "NC(=N)N"),
    # ---- Sulfur
    ("Sulfur", "Thiol", "[SX2H1]", "*S"),
    ("Sulfur", "Thioether", "[SX2H0]([#6])[#6]", "*S*"),
    ("Sulfur", "Disulfide", "[SX2H0][SX2H0]", "*SS*"),
    ("Sulfur", "Sulfoxide", "[#6][SX3](=[OX1])[#6]", "*S(=O)*"),
    ("Sulfur", "Sulfone", "[#6][SX4](=[OX1])(=[OX1])[#6]", "*S(=O)(=O)*"),
    ("Sulfur", "Sulfinic acid", "[#6][SX3](=[OX1])[OX2H1,OX1-]", "*S(=O)O"),
    ("Sulfur", "Sulfonic acid", "[SX4](=[OX1])(=[OX1])[OX2H1,OX1-]", "*S(=O)(=O)O"),
    ("Sulfur", "Sulfonate ester", "[#6][SX4](=[OX1])(=[OX1])[OX2H0][#6]", "*S(=O)(=O)O*"),
    ("Sulfur", "Sulfate ester", "[OX2]([#6])[SX4](=[OX1])(=[OX1])[OX2H1,OX1-]", "*OS(=O)(=O)O"),
    ("Sulfur", "Sulfonamide", "[NX3][SX4](=[OX1])=[OX1]", "*S(=O)(=O)N"),
    ("Sulfur", "Primary sulfonamide", "[NX3H2][SX4](=[OX1])(=[OX1])[#6]", "*S(=O)(=O)N"),
    ("Sulfur", "Sulfonyl halide", "[SX4](=[OX1])(=[OX1])[F,Cl,Br,I]", "*S(=O)(=O)Cl"),
    ("Sulfur", "Sulfonium", "[SX3+]([#6])([#6])[#6]", "*[S+](*)*"),
    # ---- P / B / Si / Se / As
    ("P / B / Si / Se / As", "Phosphate / phosphonate", "[PX4](=[OX1])([OX2,OX1-])([OX2,OX1-])[#6,OX2,OX1-]", "*OP(=O)(O)O"),
    ("P / B / Si / Se / As", "Phosphate monoester", "[OX2]([#6])[PX4](=[OX1])([OX2H1,OX1-])[OX2H1,OX1-]", "*OP(=O)(O)O"),
    ("P / B / Si / Se / As", "Phosphate diester", "[OX2]([#6])[PX4](=[OX1])([OX2H1,OX1-])[OX2][#6]", "*OP(=O)(O)O*"),
    ("P / B / Si / Se / As", "Phosphate triester", "[OX2]([#6])[PX4](=[OX1])([OX2][#6])[OX2][#6]", "*OP(=O)(O*)O*"),
    ("P / B / Si / Se / As", "Phosphonic acid", "[#6][PX4](=[OX1])([OX2H1,OX1-])[OX2H1,OX1-]", "*P(=O)(O)O"),
    ("P / B / Si / Se / As", "Diphosphate (pyrophosphate)", "[PX4](=[OX1])[OX2][PX4]=[OX1]", "*OP(=O)(O)OP(=O)(O)O"),
    ("P / B / Si / Se / As", "Phosphoramidate", "[PX4](=[OX1])[NX3]", "*NP(=O)(O)O"),
    ("P / B / Si / Se / As", "Thiophosphate (P=S)", "[PX4]=[SX1]", "*OP(=S)(O)O"),
    ("P / B / Si / Se / As", "Phosphine", "[PX3]([#6])([#6])[#6]", "*P(*)*"),
    ("P / B / Si / Se / As", "Phosphine oxide", "[#6][PX4](=[OX1])([#6])[#6]", "*P(=O)(*)*"),
    ("P / B / Si / Se / As", "Phosphonium", "[PX4+]", "*[P+](*)(*)*"),
    ("P / B / Si / Se / As", "Boronic acid", "[#6][BX3]([OX2H1])[OX2H1]", "*B(O)O"),
    ("P / B / Si / Se / As", "Boronate ester", "[#6][BX3]([OX2H0])[OX2H0]", "*B(O*)O*"),
    ("P / B / Si / Se / As", "Silanol", "[Si][OX2H1]", "*[Si](*)(*)O"),
    ("P / B / Si / Se / As", "Silyl ether", "[Si][OX2H0][#6]", "*[Si](*)(*)O*"),
    ("P / B / Si / Se / As", "Trimethylsilyl", "[Si]([CH3])([CH3])[CH3]", "C[Si](C)(C)*"),
    ("P / B / Si / Se / As", "Organoselenium", "[Se]", "*[Se]*"),
    ("P / B / Si / Se / As", "Arsenic", "[As]", "*[As](*)*"),
    # ---- Halogen
    ("Halogen", "Alkyl halide", "[F,Cl,Br,I][CX4]", "*CCl"),
    ("Halogen", "Aryl halide", "[F,Cl,Br,I][c]", "Clc1ccccc1"),
    ("Halogen", "Vinyl halide", "[F,Cl,Br,I][CX3]=[CX3]", "*C(Cl)=C*"),
    ("Halogen", "Organofluorine (C–F)", "[F][#6]", "*F"),
    ("Halogen", "Organochlorine (C–Cl)", "[Cl][#6]", "*Cl"),
    ("Halogen", "Organobromine (C–Br)", "[Br][#6]", "*Br"),
    ("Halogen", "Organoiodine (C–I)", "[I][#6]", "*I"),
    ("Halogen", "Trifluoromethyl", "[CX4](F)(F)F", "*C(F)(F)F"),
    ("Halogen", "Trichloromethyl", "[CX4](Cl)(Cl)Cl", "*C(Cl)(Cl)Cl"),
    ("Halogen", "Halide ion", "[F-,Cl-,Br-,I-]", "[Cl-]"),
    # ---- Hydrocarbon
    ("Hydrocarbon", "Methyl", "[CH3X4]", "*C"),
    ("Hydrocarbon", "Isopropyl", "[CH3][CH1X4]([CH3])[#6]", "*C(C)C"),
    ("Hydrocarbon", "tert-Butyl", "[CX4]([CH3])([CH3])([CH3])[#6]", "*C(C)(C)C"),
    ("Hydrocarbon", "Alkene", "[CX3]=[CX3]", "*C=C*"),
    ("Hydrocarbon", "Terminal alkene (vinyl)", "[CH2]=[CX3]", "C=C*"),
    ("Hydrocarbon", "Conjugated diene", "[CX3]=[CX3][CX3]=[CX3]", "C=CC=C"),
    ("Hydrocarbon", "Alkyne", "[CX2]#[CX2]", "*C#C*"),
    ("Hydrocarbon", "Terminal alkyne", "[CH1]#[CX2]", "C#C*"),
    ("Hydrocarbon", "Allene", "[CX3]=[CX2]=[CX3]", "C=C=C"),
    ("Hydrocarbon", "Cyclopropane", "[CX4;r3]1[CX4;r3][CX4;r3]1", "C1CC1"),
    ("Hydrocarbon", "Cyclobutane", "[CX4;r4]1[CX4;r4][CX4;r4][CX4;r4]1", "C1CCC1"),
    # ---- Rings
    ("Rings", "Benzene ring", "c1ccccc1", "c1ccccc1"),
    ("Rings", "Biphenyl", "c1ccccc1-c1ccccc1", "c1ccc(cc1)-c1ccccc1"),
    ("Rings", "Naphthalene", "c1ccc2ccccc2c1", "c1ccc2ccccc2c1"),
    ("Rings", "Pyridine", "n1ccccc1", "c1ccncc1"),
    ("Rings", "Pyrimidine", "c1cncnc1", "c1cncnc1"),
    ("Rings", "Pyrazine", "c1cnccn1", "c1cnccn1"),
    ("Rings", "Pyridazine", "c1ccnnc1", "c1ccnnc1"),
    ("Rings", "Triazine", "c1ncncn1", "c1ncncn1"),
    ("Rings", "Pyrrole", "n1cccc1", "c1cc[nH]c1"),
    ("Rings", "Furan", "o1cccc1", "c1ccoc1"),
    ("Rings", "Thiophene", "s1cccc1", "c1ccsc1"),
    ("Rings", "Imidazole", "c1cncn1", "c1c[nH]cn1"),
    ("Rings", "Pyrazole", "c1ccnn1", "c1cn[nH]c1"),
    ("Rings", "Oxazole", "c1cocn1", "c1cocn1"),
    ("Rings", "Isoxazole", "c1ccon1", "c1cnoc1"),
    ("Rings", "Thiazole", "c1cscn1", "c1cscn1"),
    ("Rings", "Triazole", "[#7]1:[#7]:[#6]:[#7,#6]:[#6,#7]:1", "c1nc[nH]n1"),
    ("Rings", "Tetrazole", "[#6]1:[#7]:[#7]:[#7]:[#7]:1", "c1nnn[nH]1"),
    ("Rings", "Indole", "n1ccc2ccccc12", "c1ccc2[nH]ccc2c1"),
    ("Rings", "Benzimidazole", "c1nc2ccccc2n1", "c1nc2ccccc2[nH]1"),
    ("Rings", "Benzofuran", "o1ccc2ccccc12", "c1ccc2occc2c1"),
    ("Rings", "Benzothiophene", "s1ccc2ccccc12", "c1ccc2sccc2c1"),
    ("Rings", "Quinoline", "n1cccc2ccccc12", "c1ccc2ncccc2c1"),
    ("Rings", "Isoquinoline", "c1nccc2ccccc12", "c1ccc2cnccc2c1"),
    ("Rings", "Purine", "c1ncc2ncnc2n1", "c1ncc2[nH]cnc2n1"),
    ("Rings", "Coumarin", "O=c1ccc2ccccc2o1", "O=c1ccc2ccccc2o1"),
    ("Rings", "Chromone", "O=c1ccoc2ccccc12", "O=c1ccoc2ccccc12"),
    ("Rings", "Flavone core", "O=c1cc(-c2ccccc2)oc2ccccc12", "O=c1cc(-c2ccccc2)oc2ccccc12"),
    ("Rings", "Tetrahydrofuran-type ring", "[CX4;r5]1[CX4;r5][CX4;r5][OX2;r5][CX4;r5]1", "C1CCOC1"),
    ("Rings", "Tetrahydropyran-type ring", "[CX4;r6]1[CX4;r6][CX4;r6][CX4;r6][OX2;r6][CX4;r6]1", "C1CCOCC1"),
    ("Rings", "Piperidine-type ring", "[CX4;r6]1[CX4;r6][CX4;r6][NX3;r6][CX4;r6][CX4;r6]1", "C1CCNCC1"),
    ("Rings", "Pyrrolidine-type ring", "[CX4;r5]1[CX4;r5][CX4;r5][NX3;r5][CX4;r5]1", "C1CCNC1"),
    ("Rings", "Piperazine-type ring", "[CX4;r6]1[CX4;r6][NX3;r6][CX4;r6][CX4;r6][NX3;r6]1", "C1CNCCN1"),
    ("Rings", "Morpholine-type ring", "[CX4;r6]1[CX4;r6][NX3;r6][CX4;r6][CX4;r6][OX2;r6]1", "C1COCCN1"),
    # ---- Biomolecule
    ("Biomolecule", "α-Amino acid", "[NX3;!$(N-[CX3]=[OX1])][CX4][CX3](=[OX1])[OX2H1,OX1-]", "NC(*)C(=O)O"),
    ("Biomolecule", "O-Glycoside (cyclic acetal)", "[CX4;R]([OX2;R])[OX2;!R][#6]", "C1CCC(O1)O*"),
    ("Biomolecule", "N-Glycoside (nucleoside-type)", "[CX4;R]([OX2;R])[nX3]", "C1CCC(O1)n1ccnc1"),
    ("Biomolecule", "Xanthine core", "O=c1nc(=O)c2ncnc2n1", "O=c1[nH]c(=O)c2[nH]cnc2[nH]1"),
    ("Biomolecule", "Uracil-type ring", "O=c1nc(=O)ccn1", "O=c1cc[nH]c(=O)[nH]1"),
    ("Biomolecule", "Steroid core", "[#6]1~[#6]~[#6]~[#6]2~[#6](~[#6]~1)~[#6]~[#6]~[#6]1~[#6]~2~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]~[#6]~1~2", "C1CCC2C(C1)CCC1C2CCC2CCCC12"),
    ("Biomolecule", "Phosphocholine head group", "[CH3][N+]([CH3])([CH3])[CH2][CH2][OX2][PX4](=[OX1])([OX1-,OX2H1])[OX2]", "C[N+](C)(C)CCOP(=O)([O-])O*"),
    ("Biomolecule", "Phosphoethanolamine head group", "[NX3H2,NX4H3+][CH2][CH2][OX2][PX4](=[OX1])([OX1-,OX2H1])[OX2]", "NCCOP(=O)(O)O*"),
)

FAMILIES = tuple(dict.fromkeys(family for family, *_ in RULES))
_PATTERNS = tuple((family, name, smarts, Chem.MolFromSmarts(smarts)) for family, name, smarts, _ in RULES)
_DRAWING = {name: drawing for _, name, _, drawing in RULES}
assert all(pattern is not None for *_, pattern in _PATTERNS)
assert len({name for _, name, _, _ in RULES}) == len(RULES)


def find_functional_groups(smiles: str) -> List[GroupMatch]:
    """Return every matching rule, in RULES order, with unique atom tuples.

    Raises ValueError for an empty or invalid SMILES. Different rules may cover
    the same atoms; within one rule RDKit's uniquify removes symmetry duplicates.
    """
    return _find(_parse(smiles))


def _parse(smiles: str):
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError("RDKit could not parse the SMILES.")
    return mol


def _find(mol) -> List[GroupMatch]:
    found = []
    for family, name, smarts, pattern in _PATTERNS:
        matches = tuple(mol.GetSubstructMatches(pattern, uniquify=True))
        if matches:
            found.append(GroupMatch(family=family, name=name, smarts=smarts, matches=matches))
    return found


def atom_text(matches: tuple, *, shown: int = 8) -> str:
    """One-based atom indices; semicolons separate individual group matches."""
    parts = [", ".join(str(index + 1) for index in match) for match in matches[:shown]]
    if len(matches) > shown:
        parts.append("…")
    return "; ".join(parts)


# ---- pictures ----
# soft tints of the data colours: one per site when all sites of a group are drawn together
_SITE_COLORS = tuple(theme._rgb(c, .5) for c in (theme.DATA_TEAL, theme.DATA_CORAL, theme.DATA_VIOLET, theme.DATA_AMBER, theme.DATA_BLUE, theme.DATA_GREY))
MAX_SITE_PICTURES = 6          # a group with more sites gets the all-sites picture and these first sites only


@lru_cache(maxsize=None)
def pattern_svg(name: str) -> str:
    """The group alone, drawn from its drawing SMILES with R for the rest of the molecule."""
    mol = Chem.MolFromSmiles(_DRAWING[name])
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetProp("atomLabel", "R")
    drawer = depict.svg_drawer(150, 110)
    drawer.drawOptions().padding = 0.14
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def _sites_svg(mol, sites, *, numbered: bool, size=(420, 320)) -> str:
    """The input structure with each site tinted (sites[i] gets colour i); numbered adds the one-based atom number on every tinted atom."""
    work = Chem.Mol(mol)
    atom_colors, bond_colors = {}, {}
    for i, site in enumerate(sites):
        color = _SITE_COLORS[i % len(_SITE_COLORS)]
        members = set(site)
        for index in site:
            atom_colors[index] = color
            if numbered:
                work.GetAtomWithIdx(index).SetProp("atomNote", str(index + 1))
        for bond in work.GetBonds():
            if bond.GetBeginAtomIdx() in members and bond.GetEndAtomIdx() in members:
                bond_colors[bond.GetIdx()] = color
    drawer = depict.svg_drawer(*size)
    drawer.drawOptions().annotationFontScale = 0.85
    drawer.DrawMolecule(work, highlightAtoms=list(atom_colors), highlightAtomColors=atom_colors,
                        highlightBonds=list(bond_colors), highlightBondColors=bond_colors)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def group_views(smiles: str) -> dict:
    """Everything the Properties page draws: the input structure plus, per matching group, the group picture, its sites and where they sit.

    Returns {"svg": input structure, "groups": [{"family", "name", "count", "smarts", "atoms", "svg": group picture,
    "all_svg"?: every site on the input (only with several sites), "sites": [{"atoms": "3, 2", "svg": input with that site tinted and numbered}], "more"?: sites not drawn}]}.
    """
    mol = _parse(smiles)
    base = depict.svg_drawer(420, 320)
    base.DrawMolecule(mol)
    base.FinishDrawing()
    groups = []
    for g in _find(mol):
        shown = g.matches[:MAX_SITE_PICTURES]
        item = {"family": g.family, "name": g.name, "count": g.count, "smarts": g.smarts, "atoms": atom_text(g.matches),
                "svg": pattern_svg(g.name),
                "sites": [{"atoms": ", ".join(str(i + 1) for i in site), "svg": _sites_svg(mol, [site], numbered=True)} for site in shown]}
        if g.count > 1:
            item["all_svg"] = _sites_svg(mol, g.matches, numbered=False)
        if g.count > len(shown):
            item["more"] = g.count - len(shown)
        groups.append(item)
    return {"svg": base.GetDrawingText(), "groups": groups}
