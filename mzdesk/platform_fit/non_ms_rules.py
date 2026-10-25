"""Non-MS detection from the structure: UV / PDA, native fluorescence, universal detectors (CAD / ELSD / RI). Rules only, no model.

How it works
  1. pi units: aromatic ring systems (with their exocyclic C=O / C=S) and single multiple bonds (C=C, C=O, C=N, N=N, N=O, C=S, C#C, C#N).
  2. Units joined by a single bond form one conjugated system (a chromophore). Open-chain sugars are not read as C=O (cyclic hemiacetal).
     strong = aromatic or >= 2 conjugated units;  weak = one isolated unit (C=O, C=C, amide, acid ...);  none = no unit.
  3. lambda_max per system: Woodward-Fieser (dienes, enones; Fieser-Kuhn from 5 C=C), Scott (aromatic carbonyls), else the class table below.
     The report value is the largest lambda of the strongest tier. It is approximate (+-10-20 nm); solvent and pH are not modelled.
  4. Fluorescence per aromatic system: class table plus quenchers (nitro, azo, Br / I, conjugated ketone).
  5. Detectors: UV from 2-3, FLD from 4, CAD / ELSD from the GC verdict of derivatization.gc_assess and MW, RI always "possible".
All thresholds and tables are in this file.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import rdMolDraw2D

from mzdesk.platform_fit import derivatization as deriv


# ---------------------------------------------------------------- thresholds
@dataclass(frozen=True)
class Thresholds:
    mw_cad_min: float = 150.0          # CAD / ELSD: good above this MW (Da)
    end_absorption_nm: float = 210.0   # no chromophore: absorbs below this
    usable_uv_min_nm: float = 200.0    # a weak chromophore below this is not usable (treated as none)
    low_uv_margin_nm: float = 10.0     # strong chromophore is "good" when lambda_max >= mobile phase cutoff + this
    low_uv_max_cutoff_nm: float = 205.0  # end-absorption detection needs a mobile phase cutoff at or below this
    polyene_min: int = 5               # Fieser-Kuhn from this many conjugated C=C (Woodward-Fieser below)
    max_wf_ring: int = 7               # ring terms (homoannular, exocyclic, endocyclic) only for rings up to this size, not macrocycles
    max_wf_units: int = 12             # skip the path search for bigger conjugated systems
    generic_conj_nm: float = 210.0     # unclassified conjugated system: this + step per unit beyond 2
    generic_conj_step_nm: float = 30.0
    ring_system_nm: Tuple[float, ...] = (255.0, 275.0, 330.0, 360.0)   # unclassified aromatic system by fused ring count (1, 2, 3, 4+)
    extension_nm: float = 25.0         # per conjugated C=C / C#C / C=N that no class-table hit covers
    extension_kinds: Tuple[str, ...] = ("C=C", "C#C", "C=N")
    max_label_atoms: int = 40          # atom numbers are drawn in the figure up to this many highlighted atoms
    max_fluor_ring: int = 6            # only 5- and 6-membered sp2 rings count as rigid planar fused rings


TH = Thresholds()
LAMBDA_ERROR = "±10–20 nm"
# mobile phase hint -> UV cutoff (nm) of the solvent / additive
MOBILE_PHASES = {"Water/acetonitrile + phosphate": 200.0, "Water/methanol + formic acid": 215.0}
DEFAULT_MOBILE_PHASE = "Water/acetonitrile + phosphate"
STRONG, WEAK, NONE = "strong", "weak", "none"
RANK = {NONE: 0, WEAK: 1, STRONG: 2}


# ---------------------------------------------------------------- chromophore class table
# (name, SMARTS or tuple of SMARTS, lambda_max nm of the main band, strength, text for the table)
# Patterns of substituted benzenes contain the whole ring, so they win over the plain benzene ring (a pattern whose atoms are inside another is dropped).
# Weak "azole" rings only count when they are the whole aromatic system (whole=True); otherwise the system is strong.
_NITRO = "$([NX3+](=O)[O-]),$([NX3](=O)=O)"
CLASSES = [
    # fused aromatics
    ("Pyrene", "c1cc2ccc3cccc4ccc(c1)c2c34", 335, STRONG, "~335 nm"),
    ("Anthracene", "c1ccc2cc3ccccc3cc2c1", 356, STRONG, "~356 nm"),
    ("Phenanthrene", "c1ccc2c(c1)ccc1ccccc12", 293, STRONG, "~293 nm (also ~250)"),
    ("Naphthalene", "c1ccc2ccccc2c1", 275, STRONG, "~275 nm (strong band ~220)"),
    ("Acridine / phenazine", ("c1ccc2nc3ccccc3cc2c1", "c1ccc2nc3ccccc3nc2c1"), 360, STRONG, "~360 nm"),
    ("Carbazole", "c1ccc2c(c1)[nX3]c1ccccc12", 292, STRONG, "~292 nm (also ~330)"),
    ("Fluorene / dibenzofuran type", "c1ccc2c(c1)[#6,#8,#16]c1ccccc12", 265, STRONG, "~265 nm"),
    ("Xanthene dye (fluorescein type)", "c1ccc2c(c1)[#6]c1ccccc1O2", 490, STRONG, "~490 nm (anion; ~440 neutral)"),
    ("Phenothiazine", "c1ccc2c(c1)Sc1ccccc1N2", 255, STRONG, "~255 nm (also ~305)"),
    ("Coumarin (7-oxygenated)", "O=c1ccc2ccc([OX2])cc2o1", 325, STRONG, "~325 nm"),
    ("Coumarin", "O=c1ccc2ccccc2o1", 310, STRONG, "~275 and ~310 nm"),
    ("Flavonol (3-hydroxyflavone)", "O=c1c([OX2])c(-c2ccccc2)oc2ccccc12", 365, STRONG, "~250 and ~360-370 nm"),
    ("Flavone", "O=c1cc(-c2ccccc2)oc2ccccc12", 330, STRONG, "~250-270 and ~300-350 nm"),
    ("Isoflavone", "O=c1c(-c2ccccc2)coc2ccccc12", 260, STRONG, "~260 nm"),
    ("Chromone", "O=c1ccoc2ccccc12", 300, STRONG, "~300 nm"),
    ("Anthraquinone", "O=C1c2ccccc2C(=O)c2ccccc12", 255, STRONG, "~255 nm (also ~325)"),
    ("Hydroxy- / amino-anthraquinone (emodin type)", "O=C1c2ccccc2C(=O)c2c([OX2H1,NX3])cccc12", 430, STRONG, "~430 nm (also ~250-290; coloured)"),
    ("Aminoquinone (mitomycin type)", "[#6;R]1(=O)[#6]([NX3])=,:[#6][#6;R](=O)[#6]=,:[#6]1", 360, STRONG, "~360 nm (also ~250-310)"),
    ("Quinone (benzo- / naphtho- / anthraquinone)", "[#6;R]1(=O)[#6]=,:[#6][#6;R](=O)[#6]=,:[#6]1", 250, STRONG, "~250 nm (also ~330)"),
    ("Anthrone / dianthrone (aloin, sennoside type)", "O=C1c2ccccc2[CX4]c2ccccc12", 355, STRONG, "~355 nm (also ~255, ~290)"),
    ("1,4-Dihydropyridine (nifedipine type)", "[#6]1=[#6](C(=O)[#8,#7])[#6;X4][#6](C(=O)[#8,#7])=[#6][#7]1", 238, STRONG, "~238 nm (also ~350)"),
    ("Cephem (cephalosporin)", "[#16,#8]1[#6][#6]=[#6](C(=O)[#8,#7])[#7]2[#6]1[#6][#6]2=O", 260, STRONG, "~260 nm"),
    ("Carbapenem / penem", "O=C1[#6][#6]2[#7]1C(C(=O)[#8,#7])=C[#6,#16]2", 300, STRONG, "~300 nm"),
    ("Nitroenamine (ranitidine type)", f"[NX3]-[CX3]=[CX3]-[{_NITRO}]", 320, STRONG, "~320 nm"),
    ("Quinolone (4-quinolone)", "[#8]=c1cc[nX3]c2ccccc12", 285, STRONG, "~275-295 nm (also ~320-330)"),
    ("Isoalloxazine (flavin)", "[#6]1~[#6]~[#6]2~[#7]~[#6]3~[#6](=[#8])~[#7]~[#6](=[#8])~[#7]~[#6]~3~[#7]~[#6]~2~[#6]~[#6]~1", 445, STRONG, "~445 and ~370 nm (also ~267; coloured)"),
    ("Quinoline / isoquinoline / diazanaphthalene",
     ("c1ccc2c(c1)naaa2", "c1ccc2c(c1)anaa2", "c1ccc2c(c1)aana2", "c1ccc2c(c1)aaan2"), 315, STRONG, "~315 nm (also ~225, ~275)"),
    ("Indole", "c1ccc2[nX3]ccc2c1", 280, STRONG, "~280 nm (also ~220)"),
    ("Benzofuran", "c1ccc2occc2c1", 280, STRONG, "~245 and ~280 nm"),
    ("Benzothiophene", "c1ccc2sccc2c1", 297, STRONG, "~260 and ~297 nm"),
    ("Benzazole (benzimidazole / benzoxazole / benzothiazole)", "c1ccc2[#7,#8,#16]c[#7,#8,#16]c2c1", 280, STRONG, "~275-285 nm"),
    ("Indazole / benzotriazole", ("c1ccc2[nX3]ncc2c1", "c1ccc2[nX3]nnc2c1"), 285, STRONG, "~280-295 nm"),
    ("Xanthine (caffeine, theophylline)", "[#8]=c1[nX3]c(=[#8])c2[n]c[n]c2[nX3]1", 273, STRONG, "~273 nm"),
    ("Purine (adenine, guanine)", "[#6]1:[#7]:[#6]:[#6]2:[#7]:[#6]:[#7]:[#6]:2:[#7]:1", 263, STRONG, "~260-275 nm"),
    ("Pyrimidine (uracil, cytosine, thymine)", "c1cncnc1", 260, STRONG, "~245-270 nm"),
    ("Pyrazine / pyridazine", ("c1cnccn1", "c1ccnnc1"), 255, STRONG, "~245-260 nm"),
    ("2-Pyridone", "[#8]=c1cccc[nX3]1", 295, STRONG, "~295 nm"),
    ("2-Pyranone (bufadienolide)", "[#8]=c1ccc[c,n]o1", 300, STRONG, "~295-300 nm"),
    ("Pyridine", "c1ccncc1", 257, STRONG, "~257 nm"),
    ("Thiophene", "c1ccsc1", 235, STRONG, "~235 nm"),
    ("Thiazole", "c1cscn1", 240, STRONG, "~240 nm"),
    # weak azole rings (low-UV only)
    ("Imidazole", "c1c[nX3]c[nX2]1", 210, WEAK, "~210 nm (low UV)"),
    ("Pyrazole", "c1cc[nX2][nX3]1", 210, WEAK, "~210 nm (low UV)"),
    ("Triazole / tetrazole", ("c1nc[nX3]n1", "c1c[nX3]nn1", "c1nnn[nX3]1", "c1nnc[nX3]1"), 205, WEAK, "<210 nm (low UV)"),
    ("Pyrrole", "c1cc[nX3]c1", 210, WEAK, "~210 nm (low UV)"),
    ("Furan", "c1ccoc1", 205, WEAK, "~205 nm (low UV)"),
    ("Oxazole / isoxazole / oxadiazole", ("c1cocn1", "c1ccno1", "c1nnco1", "c1ncon1"), 208, WEAK, "~205-210 nm (low UV)"),
    # benzene and its substituted forms (whole ring in every pattern)
    ("Curcuminoid (bis-cinnamoyl β-diketone)", "c-[CX3]=[CX3]-[CX3](=O)-[CX4]-[CX3](=O)-[CX3]=[CX3]-c", 425, STRONG, "~420-430 nm (enol form; coloured)"),
    ("Hydroxy- / alkoxy-cinnamoyl", "[OX2]c1ccc(cc1)-[CX3]=[CX3]-[CX3]=O", 315, STRONG, "~310-325 nm (p-coumaric, ferulic)"),
    ("Chalcone", "c1ccccc1-[CX3]=[CX3]-[CX3](=O)-c", 310, STRONG, "~310 nm"),
    ("Cinnamoyl (aryl-C=C-C=O)", "c1ccccc1-[CX3]=[CX3]-[CX3]=O", 280, STRONG, "~275-285 nm"),
    ("Stilbene", "c1ccccc1-[CX3]=[CX3]-c", 295, STRONG, "~295 nm"),
    ("Styrene (aryl alkene)", "c1ccccc1-[CX3]=[CX3]", 250, STRONG, "~250 nm"),
    ("Azobenzene (aryl azo)", "c-N=N-c", 330, STRONG, "~320-350 nm (n→π* ~440, visible)"),
    ("Acylaminophenol / alkoxyacetanilide (paracetamol type)", "[OX2]-c1ccc(cc1)-[NX3;$(N-[CX3]=O)]", 245, STRONG, "~245 nm"),
    ("Push-pull nitroarene (nitroaniline)", f"[NX3;!$(N=*);!$(N-[CX3]=O)]-c1ccc(cc1)-[{_NITRO}]", 380, STRONG, "~380 nm (4-nitroaniline)"),
    ("Push-pull nitroarene (nitrophenol / nitroanisole)", f"[OX2]-c1ccc(cc1)-[{_NITRO}]", 320, STRONG, "~310-320 nm (4-nitrophenol; ~400 as phenolate)"),
    ("Nitro-heteroarene (5-ring)", f"[c;r5]-[{_NITRO}]", 320, STRONG, "~320 nm (nitroimidazole), ~365 (nitrofuran)"),
    ("Nitroarene", f"c1ccccc1-[{_NITRO}]", 265, STRONG, "~260-270 nm"),
    ("Diaryl ketone", "c1ccccc1-[CX3](=O)-c", 255, STRONG, "~255 nm (benzophenone)"),
    ("Aryl ketone (benzoyl)", "c1ccccc1-[CX3](=O)[#6]", 245, STRONG, "~245 nm (also ~280)"),
    ("Aryl aldehyde", "c1ccccc1-[CX3H1]=O", 250, STRONG, "~250 nm"),
    ("Benzoic acid / ester", "c1ccccc1-[CX3](=O)[OX2]", 230, STRONG, "~230 nm (also ~273)"),
    ("Benzamide", "c1ccccc1-[CX3](=O)[NX3]", 230, STRONG, "~225-230 nm (also ~270)"),
    ("Benzonitrile", "c1ccccc1-C#N", 225, STRONG, "~224 nm (also ~271)"),
    ("Anilide", "[NX3;$(N-[CX3]=O)]-c1ccccc1", 245, STRONG, "~240-250 nm"),
    ("Aniline", "[NX3;!$(N=*);!$(N-[CX3]=[O,S,N]);!$(N-S(=O)=O);!$(N-[#7,#8])]-c1ccccc1", 280, STRONG, "~280 nm (also ~230; ~255 when protonated)"),
    ("Phenol / aryl ether", "[OX2]-c1ccccc1", 270, STRONG, "~270 nm (also ~220; ~290 as phenolate)"),
    ("Biaryl", "a-!@a", 255, STRONG, "~250-260 nm"),
    ("Benzene ring", "c1ccccc1", 255, STRONG, "~255 nm (weak band; strong band ~205)"),
    # isolated groups (weak)
    ("Biguanide (metformin type)", "[CX3](=[NX2])([NX3])[NX2]=[CX3]([NX3])[NX3]", 232, STRONG, "~232 nm"),
    ("Thiocarbonyl / thioamide (C=S, thiourea, thiobarbiturate)", ("[#6]=[SX1]", "[#6](=[NX2])[SX2H1,SX1-]"), 265, STRONG, "~235-290 nm (also weak n→π* ~330+)"),
    ("Aliphatic azo", "[CX4]N=N[CX4]", 350, WEAK, "~350 nm (very weak)"),
    ("Aliphatic nitro", "[CX4][NX3+](=O)[O-]", 270, WEAK, "~270 nm (very weak)"),
    ("Nitroso", "[#6][NX2]=O", 300, WEAK, "~300 nm (weak)"),
    ("Aldehyde (isolated C=O)", "[CX3H1](=O)[#6]", 290, WEAK, "~290 nm (n→π*, ε < 30); end absorption < 210"),
    ("Ketone (isolated C=O)", "[#6][CX3](=O)[#6]", 280, WEAK, "~280 nm (n→π*, ε < 30); end absorption < 210"),
    ("Acid / ester / lactone (isolated C=O)", "[CX3](=O)[OX2]", 205, WEAK, "~205 nm (end absorption)"),
    ("Amide / urea / carbamate (isolated C=O)", "[CX3](=O)[NX3]", 210, WEAK, "~205-215 nm (end absorption)"),
    ("Imine / oxime / amidine (isolated C=N)", "[CX3]=[NX2,NX3+]", 210, WEAK, "~210 nm (end absorption)"),
    ("Isolated C=C", "[CX3]=[CX3]", 200, WEAK, "~190-205 nm (end absorption)"),
    ("Isolated C≡C / C≡N", ("[CX2]#[CX2]", "[CX2]#[NX1]"), 190, WEAK, "< 200 nm (not usable)"),
]
CHROMOPHORES = [(n, tuple(s) if isinstance(s, (tuple, list)) else (s,), lam, st, text) for n, s, lam, st, text in CLASSES]
_PATTERNS = [[Chem.MolFromSmarts(s) for s in smarts] for _, smarts, *_ in CHROMOPHORES]
assert all(p is not None for ps in _PATTERNS for p in ps)
WHOLE_RING = {n for n, _, _, st, _ in CHROMOPHORES if st == WEAK and n in ("Imidazole", "Pyrazole", "Triazole / tetrazole", "Pyrrole", "Furan", "Oxazole / isoxazole / oxadiazole")}

# open-chain sugars (aldose / ketose): the C=O is really a cyclic hemiacetal, no chromophore
_OPEN_SUGAR = [Chem.MolFromSmarts(s) for s in ("[CX3;!R](=O)[CX4;!R]([OX2H1])[CX4;!R]([OX2H1])[CX4;!R][OX2H1]",
                                               "O=[CX3;!R]([CX4;!R][OX2H1])[CX4;!R]([OX2H1])[CX4;!R][OX2H1]")]

# ---------------------------------------------------------------- Woodward-Fieser / Nielsen / Fieser-Kuhn
WF_ENONE_BASE = {"ketone": 215.0, "ketone 5-ring": 202.0, "aldehyde": 207.0, "acid / ester": 195.0}
WF_EXTRA_CC = 30.0          # each further conjugated C=C
WF_EXOCYCLIC = 5.0          # per ring the C=C is exocyclic to
# substituent increments by position: alpha, beta, gamma, delta, epsilon and beyond
WF_ENONE_SUB = {"alkyl": (10, 12, 18, 18, 18), "OH": (35, 30, 50, 50, 50), "OR": (35, 30, 17, 31, 31), "OAc": (6, 6, 6, 6, 6),
                "Cl": (15, 12, 0, 0, 0), "Br": (25, 30, 0, 0, 0), "SR": (0, 85, 0, 0, 0), "NR2": (0, 95, 0, 0, 0)}
WF_DIENE_BASE = {"heteroannular / acyclic": 214.0, "homoannular": 253.0}
WF_DIENE_SUB = {"alkyl": 5, "OR": 6, "OH": 6, "OAc": 0, "SR": 30, "NR2": 60, "Cl": 5, "Br": 5}
FK_BASE, FK_PER_ALKYL, FK_N, FK_N2, FK_ENDO, FK_EXO = 114.0, 5.0, 48.0, 1.7, 16.5, 10.0

# Scott rules for aromatic carbonyls: base by carbonyl type, then increments by position (ortho, meta, para)
SCOTT_BASE = {"ketone": 246.0, "ketone 5-ring": 244.0, "ketone 7-ring": 238.0, "aldehyde": 250.0, "acid / ester": 230.0}
SCOTT_SUB = {"alkyl": (3, 3, 10), "OR": (7, 7, 25), "Cl": (0, 0, 10), "Br": (2, 2, 15),
             "NH2": (13, 13, 58), "NHAc": (20, 20, 45), "NHR": (20, 20, 73), "NR2": (20, 20, 85)}
SCOTT_NAME = {"ketone": "Aryl ketone (benzoyl)", "ketone 5-ring": "Aryl ketone (benzoyl)", "ketone 7-ring": "Aryl ketone (benzoyl)",
              "aldehyde": "Aryl aldehyde", "acid / ester": "Benzoic acid / ester"}

# ---------------------------------------------------------------- fluorescence
LIKELY, POSSIBLE, UNLIKELY = "likely", "possible", "unlikely"
FRANK = {UNLIKELY: 0, POSSIBLE: 1, LIKELY: 2}
# (name, SMARTS, level, level with a donor on the ring system (or None), typical ex / em)
FLUOR = [
    ("coumarin", "O=c1ccc2ccccc2o1", LIKELY, None, "ex ~320 / em ~400-450 nm"),
    ("indole (tryptophan-like)", "c1ccc2[nX3]ccc2c1", LIKELY, None, "ex ~280 / em ~350 nm"),
    ("4-quinolone (fluoroquinolone-like)", "[#8]=c1cc[nX3]c2ccccc12", LIKELY, None, "ex ~280 / em ~450 nm"),
    ("xanthene dye", "c1ccc2c(c1)[#6]c1ccccc1O2", LIKELY, None, "ex ~490 / em ~520 nm"),
    ("bridged biaryl (fluorene, carbazole, dibenzofuran)", "c1ccc2c(c1)[#6,#7,#8,#16]c1ccccc12", LIKELY, None, "ex ~300 / em ~320-360 nm"),
    ("quinoline-type (benzo-fused pyridine / diazine)", ("c1ccc2c(c1)naaa2", "c1ccc2c(c1)anaa2", "c1ccc2c(c1)aana2", "c1ccc2c(c1)aaan2"), POSSIBLE, LIKELY,
     "ex ~330-350 / em ~400-450 nm with a donor (quinine)"),
    ("naphthalene", "c1ccc2ccccc2c1", POSSIBLE, LIKELY, "ex ~280-330 / em ~330-350 nm"),
    ("benzazole", "c1ccc2[#7,#8,#16]c[#7,#8,#16]c2c1", POSSIBLE, None, "ex ~290 / em ~350 nm"),
    ("benzofuran / benzothiophene", ("c1ccc2occc2c1", "c1ccc2sccc2c1"), POSSIBLE, None, "ex ~290 / em ~350 nm"),
    ("flavonoid (3-OH / 5-OH patterns vary)", ("O=c1cc(-c2ccccc2)oc2ccccc12", "O=c1c(-c2ccccc2)coc2ccccc12"), POSSIBLE, None, "weak; depends on OH pattern and metal chelation"),
    ("stilbene / biaryl", ("c1ccccc1-[CX3]=[CX3]-c", "c-!@c"), POSSIBLE, None, "ex ~300 / em ~350 nm"),
    ("benzene with a donor (OH, OR, N)", "[OX2,$([NX3;!$(N-[#7,#8]);!$(N-S(=O)=O);!$([N+])])]-c1ccccc1", POSSIBLE, None, "ex ~270-290 / em ~300-350 nm (tyrosine-like)"),
]
_FLUOR_PATTERNS = [[Chem.MolFromSmarts(s) for s in ((smarts,) if isinstance(smarts, str) else smarts)] for _, smarts, *_ in FLUOR]
assert all(p is not None for ps in _FLUOR_PATTERNS for p in ps)
_DONOR = Chem.MolFromSmarts("[OX2,$([NX3;!$(N-[#7,#8]);!$([N+])])]-[a]")
_KETONE_CONJ = Chem.MolFromSmarts("[#6][CX3;!a](=O)[#6,#1]")
_BENZENE = Chem.MolFromSmarts("c1ccccc1")

# ---------------------------------------------------------------- result types
@dataclass
class Chromophore:
    name: str
    strength: str
    lam: float
    method: str                 # "Woodward–Fieser", "Fieser–Kuhn", "Scott", "class table", "generic"
    text: str                   # lambda text for the table
    atoms: List[int]            # 0-based
    count: int = 1
    hit: frozenset = frozenset()   # atoms of the class-table hit (one hit spanning several systems counts once)


@dataclass
class Analysis:
    smiles: str
    mw: float
    strength: str = NONE
    lam: Optional[float] = None
    lam_method: str = ""
    chromophores: List[Chromophore] = field(default_factory=list)
    fluor: str = UNLIKELY
    fluor_reason: str = ""
    fluor_note: str = ""
    gc: str = ""
    gc_volatility: str = ""
    reagents: Dict[str, int] = field(default_factory=dict)      # OPA / FMOC / DNPH -> number of reactive groups
    detectors: List[list] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    svg: str = ""
    mobile_phase: str = DEFAULT_MOBILE_PHASE
    cutoff: float = 200.0
    mol: Optional[Chem.Mol] = None


class NonMsError(ValueError):
    """The structure cannot be analysed."""


# ---------------------------------------------------------------- pi units and conjugated systems
_ELEMENT_ORDER = {"C": 0, "N": 1, "O": 2, "S": 3}


def _aromatic_units(mol, active):
    ri = mol.GetRingInfo()
    rings = [set(r) for r in ri.AtomRings() if set(r) <= active and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    groups: List[dict] = []
    for r in rings:
        hit = [g for g in groups if len(g["ring_atoms"] & r) >= 2]
        for g in hit:
            groups.remove(g)
        merged = {"ring_atoms": set(r), "n_rings": 1}
        for g in hit:
            merged["ring_atoms"] |= g["ring_atoms"]
            merged["n_rings"] += g["n_rings"]
        groups.append(merged)
    units = []
    for g in groups:
        atoms = set(g["ring_atoms"])
        for i in g["ring_atoms"]:                    # exocyclic double bonds (c=O of pyridone, caffeine ...) belong to the ring system
            for b in mol.GetAtomWithIdx(i).GetBonds():
                j = b.GetOtherAtomIdx(i)
                if b.GetBondType() == Chem.BondType.DOUBLE and j not in g["ring_atoms"] and j in active:
                    atoms.add(j)
        units.append({"kind": "aromatic", "atoms": atoms, "disp": set(atoms), "ring_atoms": g["ring_atoms"], "n_rings": g["n_rings"]})
    return units


def _multiple_bond_units(mol, active, aromatic_atoms):
    units = []
    for b in mol.GetBonds():
        if b.GetBondType() not in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE):
            continue
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i not in active or j not in active or i in aromatic_atoms or j in aromatic_atoms:
            continue
        a, c = mol.GetAtomWithIdx(i), mol.GetAtomWithIdx(j)
        syms = sorted((a.GetSymbol(), c.GetSymbol()), key=lambda s: _ELEMENT_ORDER.get(s, 9))
        if any(s not in _ELEMENT_ORDER for s in syms) or ("S" in syms and syms != ["C", "S"]):
            continue                                 # S=O, S=N, P=O ... are not pi chromophores (C=S is)
        kind = syms[0] + ("=" if b.GetBondType() == Chem.BondType.DOUBLE else "#") + syms[1]
        atoms, disp = {i, j}, {i, j}
        for k in (i, j):                             # for display: the hetero atoms around a C=O, the other O of a nitro group
            at = mol.GetAtomWithIdx(k)
            if kind == "C=O" and at.GetSymbol() == "C":
                disp |= {n.GetIdx() for n in at.GetNeighbors() if n.GetSymbol() in ("N", "O") and n.GetIdx() in active}
            if kind == "N=O" and at.GetSymbol() == "N":
                disp |= {n.GetIdx() for n in at.GetNeighbors() if n.GetSymbol() == "O" and n.GetIdx() in active}
        units.append({"kind": kind, "atoms": atoms, "disp": disp, "ring_atoms": set(), "n_rings": 0})
    return units


def _merge_touching(units):
    """Units that share an atom (cumulated bonds, allenes) are one unit."""
    out: List[dict] = []
    for u in units:
        hit = [v for v in out if v["atoms"] & u["atoms"]]
        for v in hit:
            out.remove(v)
            u = {"kind": "aromatic" if "aromatic" in (u["kind"], v["kind"]) else "cumulated", "atoms": u["atoms"] | v["atoms"], "disp": u["disp"] | v["disp"],
                 "ring_atoms": u["ring_atoms"] | v["ring_atoms"], "n_rings": u["n_rings"] + v["n_rings"]}
        out.append(u)
    return out


@dataclass
class Comp:
    """A conjugated system: units joined by single bonds."""
    units: List[dict]
    atoms: set = field(default_factory=set)
    disp: set = field(default_factory=set)

    @property
    def aromatic(self):
        return any(u["kind"] == "aromatic" for u in self.units)

    @property
    def kinds(self):
        return Counter(u["kind"] for u in self.units)

    @property
    def n_rings(self):
        return max([u["n_rings"] for u in self.units if u["kind"] == "aromatic"] or [0])


def conjugated_systems(mol, active) -> List[Comp]:
    arom = _aromatic_units(mol, active)
    arom_atoms = set().union(*[u["atoms"] for u in arom]) if arom else set()
    units = _merge_touching(arom + _multiple_bond_units(mol, active, arom_atoms))
    owner = {a: k for k, u in enumerate(units) for a in u["atoms"]}
    parent = list(range(len(units)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in owner and j in owner and owner[i] != owner[j] and b.GetBondType() == Chem.BondType.SINGLE:
            if units[owner[i]]["kind"] == "C=O" and units[owner[j]]["kind"] == "C=O":
                continue                             # 1,2-dicarbonyl: no conjugation to speak of
            parent[find(owner[i])] = find(owner[j])
    groups: Dict[int, List[dict]] = {}
    for k, u in enumerate(units):
        groups.setdefault(find(k), []).append(u)
    return [Comp(us, set().union(*[u["atoms"] for u in us]), set().union(*[u["disp"] for u in us])) for us in groups.values()]


# ---------------------------------------------------------------- substituent helpers
def _substituent_class(mol, atom, origin):
    """Class of the group on a chain / ring atom: alkyl, OH, OR, OAc, SR, NR2, Cl, Br, or None (ignored)."""
    sym = atom.GetSymbol()
    if sym == "C":
        return "alkyl"
    if sym == "O":
        if atom.GetTotalNumHs():
            return "OH"
        if any(n.GetIdx() != origin and _is_carbonyl_c(mol, n) for n in atom.GetNeighbors()):
            return "OAc"
        return "OR"
    return {"S": "SR", "N": "NR2", "Cl": "Cl", "Br": "Br"}.get(sym)


def _is_carbonyl_c(mol, atom):
    return atom.GetSymbol() == "C" and any(b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(atom).GetSymbol() == "O" for b in atom.GetBonds())


def _exocyclic_rings(mol, pair):
    """Rings the C=C (atom pair) is exocyclic to: exactly one of its atoms is in the ring."""
    return {tuple(sorted(r)) for r in mol.GetRingInfo().AtomRings() if len(r) <= TH.max_wf_ring and (pair[0] in r) != (pair[1] in r)}


def _endocyclic_rings(mol, pair):
    return {tuple(sorted(r)) for r in mol.GetRingInfo().AtomRings() if len(r) <= TH.max_wf_ring and pair[0] in r and pair[1] in r}


def _heavy_neighbors(atom, exclude):
    return [n for n in atom.GetNeighbors() if n.GetIdx() not in exclude]


# ---------------------------------------------------------------- lambda_max rules
def woodward_fieser(mol, comp: Comp):
    """Dienes / polyenes and enones made of C=C (and one C=O) only -> (name, lambda, method) or None."""
    kinds = comp.kinds
    n_cc, n_co = kinds.get("C=C", 0), kinds.get("C=O", 0)
    if set(kinds) - {"C=C", "C=O"} or n_cc < 1 or n_co > 1 or len(comp.units) > TH.max_wf_units or comp.aromatic:
        return None
    if n_co == 0:
        return _diene(mol, comp) if n_cc >= 2 else None
    if n_cc >= TH.polyene_min:
        return _diene(mol, Comp([u for u in comp.units if u["kind"] == "C=C"]), end_group=True)
    return _enone(mol, comp)


def _cc_pairs(comp):
    return [tuple(sorted(u["atoms"])) for u in comp.units if u["kind"] == "C=C"]


def _diene(mol, comp, end_group=False):
    pairs = _cc_pairs(comp)
    atoms = {a for p in pairs for a in p}
    n = len(pairs)
    # substituents on the sp2 carbons: alkyl / ring residues, OR, SR, NR2, halogen
    subs = Counter()
    for a in atoms:
        for nb in _heavy_neighbors(mol.GetAtomWithIdx(a), atoms):
            cls = _substituent_class(mol, nb, a)
            if cls:
                subs[cls] += 1
    exo = set().union(*[_exocyclic_rings(mol, p) for p in pairs])
    endo = set().union(*[_endocyclic_rings(mol, p) for p in pairs])
    if n >= TH.polyene_min:
        m_alkyl = subs["alkyl"]
        lam = FK_BASE + FK_PER_ALKYL * m_alkyl + n * (FK_N - FK_N2 * n) - FK_ENDO * len(endo) - FK_EXO * len(exo)
        return ("Conjugated polyene" + (" (end group ignored)" if end_group else ""), lam, "Fieser–Kuhn")
    homo = any(len(_endocyclic_rings(mol, p) & _endocyclic_rings(mol, q)) for p in pairs for q in pairs if p < q and _linked(mol, p, q))
    lam = WF_DIENE_BASE["homoannular" if homo else "heteroannular / acyclic"] + WF_EXTRA_CC * (n - 2)
    lam += sum(WF_DIENE_SUB[k] * v for k, v in subs.items()) + WF_EXOCYCLIC * len(exo)
    return ("Conjugated diene" if n == 2 else "Conjugated polyene", lam, "Woodward–Fieser")


def _linked(mol, p, q):
    return any(mol.GetBondBetweenAtoms(a, b) is not None for a in p for b in q)


def _enone(mol, comp):
    co = next(u for u in comp.units if u["kind"] == "C=O")
    c1 = next(mol.GetAtomWithIdx(i) for i in co["atoms"] if mol.GetAtomWithIdx(i).GetSymbol() == "C")
    o_idx = next(i for i in co["atoms"] if i != c1.GetIdx())
    cc = {a: u for u in comp.units if u["kind"] == "C=C" for a in u["atoms"]}
    other = [n for n in c1.GetNeighbors() if n.GetIdx() != o_idx]
    # carbonyl type from the substituent that is not the chain
    chain_starts = [n for n in other if n.GetIdx() in cc and mol.GetBondBetweenAtoms(c1.GetIdx(), n.GetIdx()).GetBondType() == Chem.BondType.SINGLE]
    if not chain_starts:
        return None
    best = None
    for start in chain_starts:
        rest = [n for n in other if n.GetIdx() != start.GetIdx()]
        if len(rest) > 1:
            continue
        r = rest[0] if rest else None
        if r is None:
            ctype = "aldehyde"
        elif r.GetSymbol() == "C" and not r.GetIsAromatic():
            ctype = "ketone 5-ring" if c1.IsInRingSize(5) else "ketone"
        elif r.GetSymbol() in ("O", "N"):
            ctype = "acid / ester"                   # acid, ester, lactone; amides follow the ester rule
        else:
            continue
        for path in _enone_paths(mol, cc, c1.GetIdx(), start.GetIdx()):
            lam = _enone_lambda(mol, c1.GetIdx(), path, WF_ENONE_BASE[ctype], ctype == "acid / ester")
            if lam is not None and (best is None or lam > best[1]):
                best = (ctype, lam, len(path) // 2)
    if best is None:
        return None
    ctype, lam, n = best
    name = "α,β-Unsaturated acid / ester / amide" if ctype == "acid / ester" else ("α,β-Unsaturated aldehyde" if ctype == "aldehyde" else "α,β-Unsaturated ketone (enone)")
    if n > 1:
        name = name.replace("α,β-Unsaturated", "Conjugated dienyl" if n == 2 else "Conjugated polyenyl")
    return (name, lam, "Woodward–Fieser")


def _enone_paths(mol, cc, c1, start):
    """Chains [alpha, beta, gamma, delta ...] of conjugated C=C bonds leading away from the carbonyl carbon."""
    paths = []

    def partner(a):
        return next(i for i in cc[a]["atoms"] if i != a)

    def walk(chain, used):
        paths.append(list(chain))
        for nb in mol.GetAtomWithIdx(chain[-1]).GetNeighbors():
            k = nb.GetIdx()
            if k in cc and id(cc[k]) not in used and mol.GetBondBetweenAtoms(chain[-1], k).GetBondType() == Chem.BondType.SINGLE:
                walk(chain + [k, partner(k)], used | {id(cc[k])})
    walk([start, partner(start)], {id(cc[start])})
    return paths


def _enone_lambda(mol, c1, chain, base, acid_like=False):
    lam = base + WF_EXTRA_CC * (len(chain) // 2 - 1)
    for pos, idx in enumerate(chain):
        prev = c1 if pos == 0 else chain[pos - 1]
        nxt = chain[pos + 1] if pos + 1 < len(chain) else None
        for nb in _heavy_neighbors(mol.GetAtomWithIdx(idx), {prev, nxt}):
            if pos == 0 and nb.GetSymbol() == "N":
                return None                          # alpha amino / enamide acids: not covered by the rules
            cls = _substituent_class(mol, nb, idx)
            if cls in ("NR2", "SR") and acid_like:
                return None                          # beta amino / thio acrylates: the increments are for ketones only
            if cls:
                lam += WF_ENONE_SUB[cls][min(pos, 4)]
    for k in range(0, len(chain), 2):
        lam += WF_EXOCYCLIC * len(_exocyclic_rings(mol, (chain[k], chain[k + 1])))
    return lam


def scott(mol, comp: Comp):
    """One benzene ring carrying aromatic carbonyl(s) and nothing else conjugated -> (name, lambda, method) or None."""
    kinds = comp.kinds
    arom = [u for u in comp.units if u["kind"] == "aromatic"]
    if len(arom) != 1 or set(kinds) - {"aromatic", "C=O"} or not kinds.get("C=O"):
        return None
    ring_atoms = arom[0]["ring_atoms"]
    if arom[0]["n_rings"] != 1 or len(ring_atoms) != 6 or any(mol.GetAtomWithIdx(i).GetSymbol() != "C" for i in ring_atoms):
        return None
    ring = next(list(r) for r in mol.GetRingInfo().AtomRings() if set(r) == ring_atoms)
    best = None
    for u in comp.units:
        if u["kind"] != "C=O":
            continue
        c = next(mol.GetAtomWithIdx(i) for i in u["atoms"] if mol.GetAtomWithIdx(i).GetSymbol() == "C")
        o_idx = next(i for i in u["atoms"] if i != c.GetIdx())
        anchor = [n.GetIdx() for n in c.GetNeighbors() if n.GetIdx() in ring_atoms]
        if len(anchor) != 1:
            continue
        rest = [n for n in c.GetNeighbors() if n.GetIdx() not in (o_idx, anchor[0])]
        if len(rest) > 1:
            continue
        r = rest[0] if rest else None
        if r is None:
            ctype = "aldehyde"
        elif r.GetSymbol() == "C" and not r.GetIsAromatic() and not any(b.GetBondType() != Chem.BondType.SINGLE for b in r.GetBonds()):
            ctype = "ketone 5-ring" if c.IsInRingSize(5) else ("ketone 7-ring" if c.IsInRingSize(7) else "ketone")
        elif r.GetSymbol() == "O":
            ctype = "acid / ester"
        else:
            continue                                 # amides, vinyl / aryl ketones: class table
        lam = SCOTT_BASE[ctype]
        pos0 = ring.index(anchor[0])
        for i in ring:
            if i == anchor[0]:
                continue
            d = abs(ring.index(i) - pos0)
            d = min(d, 6 - d)                        # 1 ortho, 2 meta, 3 para
            for nb in mol.GetAtomWithIdx(i).GetNeighbors():
                if nb.GetIdx() in ring_atoms or nb.GetIdx() == c.GetIdx() or nb.GetSymbol() == "H":
                    continue
                cls = _scott_class(mol, nb, i)
                if cls:
                    lam += SCOTT_SUB[cls][d - 1]
        if best is None or lam > best[1]:
            best = (SCOTT_NAME[ctype], lam, "Scott")
    return best


def _scott_class(mol, atom, origin):
    sym = atom.GetSymbol()
    if sym == "C":
        return "alkyl" if not atom.GetIsAromatic() and all(b.GetBondType() == Chem.BondType.SINGLE for b in atom.GetBonds()) else None
    if sym == "O":
        return "OR"
    if sym == "Cl":
        return "Cl"
    if sym == "Br":
        return "Br"
    if sym == "N":
        if any(_is_carbonyl_c(mol, n) for n in atom.GetNeighbors() if n.GetIdx() != origin):
            return "NHAc"
        if atom.GetFormalCharge() or any(b.GetBondType() != Chem.BondType.SINGLE for b in atom.GetBonds()) or any(n.GetSymbol() == "O" for n in atom.GetNeighbors()):
            return None                              # nitro and similar
        hs = atom.GetTotalNumHs()
        return {2: "NH2", 1: "NHR"}.get(hs, "NR2")
    return None


# ---------------------------------------------------------------- classify the systems
def _candidates(comp, matches):
    """Class-table hits inside one system, after dropping hits whose atoms are inside another hit -> [(entry, atoms)]."""
    hits = []
    for atoms, k in matches:
        if not (atoms & comp.atoms):
            continue
        if CHROMOPHORES[k][0] in WHOLE_RING and not (len(comp.units) == 1 and comp.aromatic and comp.units[0]["ring_atoms"] <= atoms):
            continue
        hits.append((atoms, k))
    hits.sort(key=lambda h: (-len(h[0]), h[1]))
    kept = []
    for atoms, k in hits:
        if not any(atoms <= other for other, _ in kept):
            kept.append((atoms, k))
    return [(CHROMOPHORES[k], atoms) for atoms, k in kept]


def _extension(comp, covered):
    """Extra nm for conjugated C=C / C#C / C=N units that no class-table hit covers (they lengthen the chromophore)."""
    n = sum(1 for u in comp.units if u["kind"] in TH.extension_kinds and not (u["atoms"] <= covered))
    return n, n * TH.extension_nm


def classify(mol, comp, matches) -> Optional[Chromophore]:
    atoms = sorted(comp.disp)
    strong_graph = comp.aromatic or len(comp.units) >= 2
    calc = woodward_fieser(mol, comp) or scott(mol, comp)
    if calc:
        name, lam, method = calc
        return Chromophore(name, STRONG, lam, method, f"~{lam:.0f} nm ({method})", atoms)
    cands = _candidates(comp, matches)
    if strong_graph:
        pool = [c for c in cands if c[0][3] == STRONG] or ([c for c in cands if c[0][0] in WHOLE_RING] if len(comp.units) == 1 else [])
    else:
        pool = cands
    covered = set().union(*[at for e, at in cands if e[3] == STRONG] or [set()])
    if pool:
        (name, _, lam, strength, text), hit_atoms = max(pool, key=lambda c: (RANK[c[0][3]], c[0][2]))
        n_ext, add = _extension(comp, covered) if strength == STRONG else (0, 0.0)
        if n_ext:
            lam, text = lam + add, f"~{lam + add:.0f} nm (class value {lam:.0f} + {add:.0f} for {n_ext} more conjugated bond{'s' if n_ext > 1 else ''})"
        return Chromophore(name, strength, lam, "class table", text, atoms, hit=frozenset(hit_atoms))
    if comp.aromatic:
        n = min(comp.n_rings, len(TH.ring_system_nm)) - 1
        base = TH.ring_system_nm[max(n, 0)]
        n_ext, add = _extension(comp, set())
        lam = base + add
        label = "heteroaromatic" if any(mol.GetAtomWithIdx(i).GetSymbol() != "C" for u in comp.units for i in u["ring_atoms"]) else "aromatic"
        name = f"Fused {label} system ({comp.n_rings} rings)" if comp.n_rings > 1 else f"{label.capitalize()} ring"
        if n_ext:
            name += f" + {n_ext} conjugated bond{'s' if n_ext > 1 else ''}"
        return Chromophore(name, STRONG, lam, "generic", f"~{lam:.0f} nm (generic: by ring count" + (f", +{add:.0f} for extra conjugation)" if n_ext else ")"), atoms)
    if len(comp.units) >= 2:
        lam = TH.generic_conj_nm + TH.generic_conj_step_nm * (len(comp.units) - 2)
        return Chromophore(f"Conjugated π system ({len(comp.units)} multiple bonds)", STRONG, lam, "generic", f"~{lam:.0f} nm (generic: by length)", atoms)
    return Chromophore("Isolated π bond (" + comp.units[0]["kind"] + ")", WEAK, 210.0, "generic", "~210 nm (end absorption)", atoms)


def _is_open_sugar_carbonyl(mol):
    return {m[0] if mol.GetAtomWithIdx(m[0]).GetSymbol() == "C" else m[1] for p in _OPEN_SUGAR for m in mol.GetSubstructMatches(p)}


def find_chromophores(mol, active) -> List[Chromophore]:
    comps = conjugated_systems(mol, active)
    sugar_c = _is_open_sugar_carbonyl(mol)
    if sugar_c:
        # an open-chain sugar C=O is a hemiacetal in solution: drop its unit (alone) so nothing is reported for it
        comps = [c for c in comps if not (len(c.units) == 1 and c.units[0]["kind"] == "C=O" and c.units[0]["atoms"] & sugar_c)]
    matches = []
    for k, patterns in enumerate(_PATTERNS):
        seen = set()
        for p in patterns:
            for m in mol.GetSubstructMatches(p, uniquify=True, maxMatches=2000):
                s = frozenset(m)
                if s <= active and s not in seen:
                    seen.add(s)
                    matches.append((set(s), k))
    found = [classify(mol, c, matches) for c in comps]
    return _group([f for f in found if f])


def _group(found: List[Chromophore]) -> List[Chromophore]:
    """One row per class: same name, strength and lambda are merged (atoms joined, count)."""
    rows: Dict[tuple, Chromophore] = {}
    seen_hits: Dict[tuple, set] = {}
    for f in found:
        key = (f.name, f.strength, round(f.lam))
        if key in rows:
            rows[key].atoms = sorted(set(rows[key].atoms) | set(f.atoms))
            if not f.hit or f.hit not in seen_hits.setdefault(key, set()):
                rows[key].count += 1
        else:
            rows[key] = Chromophore(f.name, f.strength, f.lam, f.method, f.text, list(f.atoms))
        seen_hits.setdefault(key, set()).add(f.hit)
    return sorted(rows.values(), key=lambda r: (-RANK[r.strength], -r.lam))


def headline(chromophores: List[Chromophore]):
    """-> (strength, lambda or None, method, name): the largest lambda of the strongest tier; a weak lambda below the usable window counts as none."""
    for tier in (STRONG, WEAK):
        rows = [c for c in chromophores if c.strength == tier]
        if rows:
            best = max(rows, key=lambda c: c.lam)
            if tier == WEAK and best.lam < TH.usable_uv_min_nm:
                return NONE, None, "", best.name
            return tier, best.lam, best.method, best.name
    return NONE, None, "", ""


# ---------------------------------------------------------------- fluorescence
def _fused_sp2_rings(mol, comp, active):
    """Largest group of fused rings that are all sp2 (aromatic, or carbonyl / alkene ring atoms), rings of 5 or 6 atoms only."""
    ri = mol.GetRingInfo()

    def sp2(i):
        a = mol.GetAtomWithIdx(i)
        return a.GetIsAromatic() or a.GetHybridization() in (Chem.HybridizationType.SP2, Chem.HybridizationType.SP)
    rings = [set(r) for r in ri.AtomRings() if len(r) <= TH.max_fluor_ring and set(r) <= active and all(sp2(i) for i in r) and set(r) & comp.atoms]
    best, seen = 0, set()
    for start in range(len(rings)):
        if start in seen:
            continue
        group, stack = {start}, [start]
        while stack:
            k = stack.pop()
            for m in range(len(rings)):
                if m not in group and len(rings[k] & rings[m]) >= 2:
                    group.add(m)
                    stack.append(m)
        seen |= group
        best = max(best, len(group))
    return best


def fluorescence(mol, comps, active):
    """-> (level, reason, ex/em note) over all aromatic systems."""
    best = (UNLIKELY, "", "")
    reasons = []
    for comp in comps:
        if not comp.aromatic:
            continue
        level, why, note = UNLIKELY, "", ""
        n_fused = _fused_sp2_rings(mol, comp, active)
        if n_fused >= 3:
            level, why, note = LIKELY, f"{n_fused} fused sp2 rings: rigid and planar", "ex / em depend on the ring size"
        for (name, _, lvl, lvl_donor, ex_em), patterns in zip(FLUOR, _FLUOR_PATTERNS):
            hit = [set(m) for p in patterns for m in mol.GetSubstructMatches(p) if set(m) & comp.atoms and set(m) <= active]
            if not hit:
                continue
            cand = lvl
            if lvl_donor and any(_has_donor(mol, h) for h in hit):
                cand = lvl_donor
            if FRANK[cand] > FRANK[level]:
                level, why, note = cand, name + (" with a donor" if cand == lvl_donor and lvl_donor else ""), ex_em
        quench = _quenchers(mol, comp, active)
        if level != UNLIKELY and quench:
            hard = [q for q in quench if q[1]]
            if hard:
                reasons.append(f"{why} quenched by {hard[0][0]}")
                level, why = UNLIKELY, ""
            else:
                level = POSSIBLE if level == LIKELY else UNLIKELY
                why += f"; reduced by {quench[0][0]}"
        if FRANK[level] > FRANK[best[0]] or (level != UNLIKELY and level == best[0] and not best[1]):
            best = (level, why, note)
    if best[0] == UNLIKELY and not best[1]:
        if reasons:
            best = (UNLIKELY, reasons[0], "")
        elif any(c.aromatic for c in comps):
            benzene = mol.HasSubstructMatch(_BENZENE)
            best = (UNLIKELY, "benzene ring without a donor or a fused fluorophore" if benzene else "heteroaromatic ring without a known fluorophore pattern", "")
        else:
            best = (UNLIKELY, "no aromatic system", "")
    return best


def _has_donor(mol, atoms):
    """Donor (OH, OR, amine) bonded to an atom of the matched ring system."""
    return any(m[0] not in atoms and m[1] in atoms for m in mol.GetSubstructMatches(_DONOR))


def _quenchers(mol, comp, active):
    """[(name, hard)]: hard quenchers (nitro, azo, Br, I) switch fluorescence off; a conjugated ketone only lowers it."""
    out = []
    kinds = comp.kinds
    if kinds.get("N=O"):
        out.append(("nitro group", True))
    if kinds.get("N=N"):
        out.append(("azo group", True))
    heavy = {n.GetSymbol() for u in comp.units if u["kind"] == "aromatic" for i in u["ring_atoms"]
             for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetSymbol() in ("Br", "I")}
    if heavy:
        out.append(("heavy atom (" + " / ".join(sorted(heavy)) + ")", True))
    if any(set(m) & comp.atoms for m in mol.GetSubstructMatches(_KETONE_CONJ)) and any(u["kind"] == "C=O" for u in comp.units):
        out.append(("conjugated ketone (n→π* state)", False))
    return out


# ---------------------------------------------------------------- detectors and notes
_LABEL_REAGENTS = {"OPA": "UV/FLD · OPA + 2-mercaptoethanol (primary amine)", "FMOC": "UV/FLD · FMOC-Cl (amine)", "DNPH": "UV/FLD · DNPH (carbonyl)"}


def _reagents(smiles):
    res = deriv.derive(smiles, list(_LABEL_REAGENTS.values()))
    by_opt = {r["reagent"]: r for r in res["results"]}
    out = {}
    for short, opt in _LABEL_REAGENTS.items():
        if by_opt[opt]["n"] > 0:
            out[short] = by_opt[opt]["n"]
    return out


def deriv_hint(reagents):
    """Sentence naming the OPA / FMOC / DNPH routes that add a chromophore or fluorophore ("" when no group reacts)."""
    return "Derivatization adds one: " + " / ".join(f"{k} ({v}×)" for k, v in reagents.items()) + " (Derivatization page)." if reagents else ""


def _uv_row(strength, lam, name, cutoff):
    if strength == STRONG:
        if lam >= cutoff + TH.low_uv_margin_nm:
            return ["good", f"Strong chromophore ({name}), λmax ~{lam:.0f} nm: detect there; PDA gives the spectrum for peak purity."]
        return ["possible", f"Strong chromophore but absorbing near the mobile phase cutoff (~{cutoff:.0f} nm): λmax ~{lam:.0f} nm; "
                f"use a transparent mobile phase (acetonitrile, phosphate) and detect at 200–{lam + 10:.0f} nm."]
    if strength == WEAK:
        if lam >= cutoff + TH.low_uv_margin_nm:
            return ["possible", f"Weak chromophore ({name}); band at ~{lam:.0f} nm has a low absorptivity (ε < ~100): poor sensitivity."]
        if cutoff <= TH.low_uv_max_cutoff_nm:
            return ["possible", f"End absorption only ({name}): detect at 200–215 nm; needs a transparent mobile phase (cutoff ~{cutoff:.0f} nm), high purity solvents, isocratic or a flat baseline."]
        return ["poor", f"End absorption only ({name}) at 200–215 nm, but this mobile phase cuts off at ~{cutoff:.0f} nm: switch to acetonitrile / phosphate, or use CAD / ELSD."]
    return ["poor", "No chromophore: no absorption above ~210 nm."]


def analyze(smiles: str, mobile_phase: Optional[str] = None) -> Analysis:
    mobile_phase = mobile_phase or DEFAULT_MOBILE_PHASE
    if mobile_phase not in MOBILE_PHASES:
        raise NonMsError(f"Unknown mobile phase: {mobile_phase}.")
    mol0 = Chem.MolFromSmiles(smiles)
    if mol0 is None:
        raise NonMsError("RDKit could not parse the structure.")
    mol, neutralized = deriv.neutralize(mol0)
    frags = Chem.GetMolFrags(mol)
    main = max(frags, key=lambda f: (sum(mol.GetAtomWithIdx(i).GetAtomicNum() > 1 for i in f), sum(mol.GetAtomWithIdx(i).GetSymbol() == "C" for i in f)))
    active = set(main)
    frag = Chem.GetMolFrags(mol, asMols=True)[frags.index(main)]
    a = Analysis(smiles=smiles, mw=Descriptors.MolWt(frag), mobile_phase=mobile_phase, cutoff=MOBILE_PHASES[mobile_phase], mol=mol)
    if neutralized:
        a.warnings.append("Input neutralized (charges removed).")
    if len(frags) > 1:
        a.warnings.append(f"{len(frags)} disconnected parts (salt or mixture): only the largest part ({len(main)} atoms) is analysed; atom numbers are those of the input.")

    a.chromophores = find_chromophores(mol, active)
    a.strength, a.lam, a.lam_method, top = headline(a.chromophores)
    comps = conjugated_systems(mol, active)
    a.fluor, a.fluor_reason, a.fluor_note = fluorescence(mol, comps, active)

    groups = deriv.find_groups(frag)
    gc = deriv.gc_assess(frag, groups)
    a.gc, a.gc_volatility = gc["gc"], gc["volatility"]
    frag_smiles = Chem.MolToSmiles(frag)
    reagents = a.reagents = _reagents(frag_smiles)
    a.detectors = _detector_rows(a, top, reagents)
    a.notes = _notes(a, top, reagents, bool(_is_open_sugar_carbonyl(mol)))
    a.svg = highlight_svg(mol, a.chromophores) or ""
    return a


def _detector_rows(a: Analysis, top: str, reagents: dict) -> List[list]:
    uv = _uv_row(a.strength, a.lam or 0.0, top, a.cutoff)
    if a.strength == NONE and top:
        uv[1] = f"Only {top} (absorbs below {TH.usable_uv_min_nm:.0f} nm): not usable for UV."
    deriv_text = " " + deriv_hint(reagents) if reagents else ""
    fld_suit = {LIKELY: "good", POSSIBLE: "possible", UNLIKELY: "poor"}[a.fluor]
    fld = [fld_suit, {LIKELY: f"Native fluorescence likely: {a.fluor_reason} ({a.fluor_note}).",
                      POSSIBLE: f"Native fluorescence possible: {a.fluor_reason} ({a.fluor_note}); verify with a standard.",
                      UNLIKELY: f"Native fluorescence unlikely: {a.fluor_reason}."}[a.fluor] + (deriv_text if a.fluor != LIKELY else "")]
    volatile = a.gc == "yes"
    if volatile:
        cad = ["poor", f"GC verdict yes ({a.gc_volatility}): treated as volatile, so it largely evaporates in the CAD / ELSD drying step and the response is lost or very low."]
    elif a.mw >= TH.mw_cad_min:
        cad = ["good", f"Non-volatile (GC verdict {a.gc}), MW {a.mw:.0f} ≥ {TH.mw_cad_min:.0f}: response does not depend on a chromophore. Mobile phase must be volatile (formic / acetic acid, ammonium formate / acetate; no phosphate)."]
    else:
        cad = ["possible", f"Non-volatile but small (MW {a.mw:.0f} < {TH.mw_cad_min:.0f}): response is low and the baseline noise from the mobile phase matters. Mobile phase must be volatile."]
    elsd = [cad[0], cad[1].replace("response does not depend on a chromophore", "response does not depend on a chromophore (less sensitive than CAD, non-linear)")
            if cad[0] == "good" else cad[1]]
    ri_extra = " UV or FLD is far more sensitive here." if a.strength == STRONG else ""
    ri = ["possible", "Universal but insensitive (µg level): needs an isocratic method (no gradient), a thermostatted detector and a high concentration." + ri_extra]
    return [["UV/PDA"] + uv, ["FLD"] + fld, ["CAD"] + cad, ["ELSD"] + elsd, ["RI"] + ri]


def _notes(a: Analysis, top: str, reagents: dict, sugar: bool) -> List[str]:
    notes = []
    if a.chromophores:
        parts = []
        for tier in (STRONG, WEAK):
            rows = [c for c in a.chromophores if c.strength == tier]
            if rows:
                text = "; ".join(f"{c.name}{' ×' + str(c.count) if c.count > 1 else ''}: {c.text}" for c in rows[:6])
                parts.append(text if tier == a.strength else f"Also {tier}: {text}")
        notes.append(f"UV rule: {a.strength} — " + ". ".join(parts) + ".")
    elif sugar:
        notes.append("UV rule: only an open-chain sugar C=O, which is a cyclic hemiacetal in solution: no chromophore.")
    else:
        notes.append("UV rule: no aromatic ring or multiple bond: no chromophore (end absorption below ~210 nm).")
    if a.lam is not None:
        notes.append(f"λmax rule: {a.lam_method} gives ~{a.lam:.0f} nm for {top} (approximate, {LAMBDA_ERROR}).")
    elif a.chromophores:
        notes.append(f"λmax rule: no usable band above {TH.usable_uv_min_nm:.0f} nm ({top} absorbs below).")
    else:
        notes.append(f"λmax rule: no band above {TH.end_absorption_nm:.0f} nm expected.")
    near_cutoff = (a.lam or 999) < a.cutoff + TH.low_uv_margin_nm
    if a.strength == WEAK and near_cutoff:
        if a.cutoff <= TH.low_uv_max_cutoff_nm:
            notes.append(f"Low-UV rule: with {a.mobile_phase} (cutoff ~{a.cutoff:.0f} nm) detect the end absorption at 200–215 nm; use a phosphate buffer, not formate or acetate.")
        else:
            notes.append(f"Low-UV rule: {a.mobile_phase} cuts off at ~{a.cutoff:.0f} nm, above the end absorption (200–215 nm): use acetonitrile with phosphate for UV.")
    elif a.strength == STRONG and near_cutoff:
        notes.append(f"Low-UV rule: λmax ~{a.lam:.0f} nm is near the cutoff of {a.mobile_phase} (~{a.cutoff:.0f} nm).")
    notes.append(f"Fluorescence rule: {a.fluor}" + (f" — {a.fluor_reason}" if a.fluor_reason else "") + (f" ({a.fluor_note})." if a.fluor_note else "."))
    verdict = "volatile as is, so poor" if a.gc == "yes" else ("non-volatile and MW ≥ " + format(TH.mw_cad_min, ".0f") + ", so good" if a.mw >= TH.mw_cad_min
                                                             else "non-volatile but MW < " + format(TH.mw_cad_min, ".0f") + ", so possible")
    notes.append(f"CAD / ELSD rule: GC verdict “{a.gc}” ({a.gc_volatility}), MW {a.mw:.0f}: {verdict}. The mobile phase must be volatile.")
    notes.append("RI rule: universal, but detection limits are in the µg range (100–1000× worse than UV for a strong chromophore); isocratic only, no gradient.")
    if a.strength in (NONE, WEAK):
        if reagents:
            names = {"OPA": "OPA + 2-mercaptoethanol (primary amine, fluorescence)", "FMOC": "FMOC-Cl (amine / phenol, UV and fluorescence)", "DNPH": "DNPH (carbonyl, UV ~360 nm)"}
            notes.append("Derivatization option: " + "; ".join(names[k] for k in reagents) + ". See the Derivatization page.")
        else:
            notes.append("Derivatization option: no group for OPA / FMOC / DNPH (no amine or carbonyl): use CAD / ELSD, RI or MS instead.")
    notes.append(f"Caveat: λmax values are class rules, not predictions ({LAMBDA_ERROR}); solvent and pH shift them (phenols, anilines, acids and pyridines move 20–50 nm with ionization). "
                 "Concentration, path length and slit width decide the real sensitivity.")
    notes.append("Caveat: fluorescence is a structural alert only; quantum yield depends on solvent, pH and quenchers.")
    return notes


# ---------------------------------------------------------------- figure
STRONG_COLOR = (0.30, 0.43, 0.88)    # indigo-blue
WEAK_COLOR = (0.96, 0.66, 0.15)      # amber


def highlight_svg(mol, chromophores: List[Chromophore], size=(560, 420)) -> Optional[str]:
    """The molecule with the atoms of each chromophore highlighted (strong: blue, weak: amber; atom numbers as in the table)."""
    colors: Dict[int, tuple] = {}
    for c in sorted(chromophores, key=lambda c: RANK[c.strength]):      # strong drawn last
        for i in c.atoms:
            colors[i] = STRONG_COLOR if c.strength == STRONG else WEAK_COLOR
    drawer = rdMolDraw2D.MolDraw2DSVG(*size)
    opts = drawer.drawOptions()
    opts.clearBackground = False
    opts.addStereoAnnotation = False
    opts.annotationFontScale = 0.7
    m = Chem.Mol(mol)
    if colors and len(colors) <= TH.max_label_atoms:
        for i in colors:
            m.GetAtomWithIdx(i).SetProp("atomNote", str(i + 1))
    bonds = [b.GetIdx() for b in m.GetBonds() if b.GetBeginAtomIdx() in colors and b.GetEndAtomIdx() in colors]
    bond_colors = {b: colors[m.GetBondWithIdx(b).GetBeginAtomIdx()] for b in bonds}
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, m, highlightAtoms=list(colors), highlightAtomColors=colors, highlightBonds=bonds, highlightBondColors=bond_colors)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()
