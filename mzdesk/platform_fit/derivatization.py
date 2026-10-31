"""Derivatization：找反应基团，按试剂算衍生物的分子式和单同位素质量，并给 Platform 页判断 CIL / GC-MS 适用性。

母体取画出来的中性分子。质量全部用精确同位素质量（不用平均质量）。
每个试剂还带一份官能团反应式（SCHEMES）；四个 CIL Channel 均生成完整产物结构，Platform 和 Derivatization 两页的 reaction 块共用。
"""
from __future__ import annotations

import re
from collections import Counter
from typing import Dict, List, Optional

from rdkit import Chem
from rdkit.Chem import Descriptors, rdCoordGen
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.rdBase import BlockLogs

from mzdesk.shared import colors as _theme, depict as _depict

C13_SHIFT = 1.0033548          # 13C - 12C
PT = Chem.GetPeriodicTable()

# 基团：key -> (名称, SMARTS, 代表原子在 SMARTS 里的序号)
GROUPS = {
    "acid": ("carboxylic acid", "[CX3](=O)[OX2H1]", 2),
    "phenol": ("phenol OH", "[OX2H1][c]", 0),
    "alcohol": ("alcohol OH", "[OX2H1][CX4]", 0),
    "amine1": ("primary amine", "[NX3;H2;!$(N[#6]=[O,S,N]);!$(NS(=O)=O)]", 0),
    "amine2": ("secondary amine", "[NX3;H1;!$(N[#6]=[O,S,N]);!$(NS(=O)=O);!$(N=*)]", 0),
    "amide": ("amide / imide NH", "[NX3;H1,H2;$(N[#6]=[O,S,N])]", 0),
    "arom_nh": ("aromatic NH", "[nH]", 0),                                    # 只用于 TMS 的活泼氢
    "arom_nh_cil": ("aromatic NH (not next to C=O)", "[nX3;H1,H2;!$(nc=O)]", 0),   # CIL 通道（Zhao 2019）
    "thiol": ("thiol", "[SX2H1]", 0),
    "enol": ("enol OH", "[CX3,NX2]=[CX3][OX2H1]", 2),
    "aldehyde": ("aldehyde", "[CX3H1](=O)[#6]", 0),
    "ketone": ("ketone", "[#6][CX3](=O)[#6]", 1),
    "ketone_ar": ("ketone (aromatic C=O)", "[#6][cX3](=O)[#6]", 1),
}
_PATTERNS = {k: Chem.MolFromSmarts(v[1]) for k, v in GROUPS.items()}
_HEAD = ("acid", "phenol", "alcohol", "enol", "amine1", "amine2", "amide", "arom_nh", "arom_nh_cil", "thiol",
         "aldehyde", "ketone", "ketone_ar")
# TMS 取代的 H；每个原子只属于其中一类，arom_nh_cil 是 arom_nh 的子集，所以不放进来
ACTIVE_H = ("acid", "phenol", "alcohol", "enol", "amine1", "amine2", "amide", "arom_nh", "thiol")
CARBONYL = ("aldehyde", "ketone", "ketone_ar")
_TERTIARY_AMINE = Chem.MolFromSmarts("[NX3;H0;!$(N[#6]=[O,S,N]);!$(NS(=O)=O);!$(N-[a])]")
_HEMIACETAL = Chem.MolFromSmarts("[OX2H1][CX4;R][OX2;R]")

TMS, MEOX, TBDMS = {"C": 3, "H": 8, "Si": 1}, {"C": 1, "H": 3, "N": 1}, {"C": 6, "H": 14, "Si": 1}
TFAA = {"C": 2, "H": -1, "F": 3, "O": 1}                    # CF3CO 取代 H
DANSYL = {"C": 12, "H": 11, "N": 1, "O": 2, "S": 1}         # Dns(C12H12NO2S) 取代 H
DMPA = {"C": 10, "H": 11, "N": 1, "O": 1}
DNS_HYDRAZINE = {"C": 12, "H": 13, "N": 3, "O": 1, "S": 1}  # C12H15N3O2S - H2O
# OPA(C8H6O2) + 2-巯基乙醇(C2H6OS) + RNH2 - 2 H2O：净增 C10H12O3S - H4O2 = C10H8OS（异吲哚）
OPA = {"C": 10, "H": 8, "O": 1, "S": 1}
FMOC = {"C": 15, "H": 10, "O": 2}                           # C15H11ClO2 - HCl
DNPH = {"C": 6, "H": 4, "N": 4, "O": 3}                     # C6H6N4O4 - H2O

# 试剂：(名称, 反应基团 {key: "H" 按氢数 / "atom" 按原子数 / int 固定每原子标记数}, delta, 是否 CIL 同位素标记)
_ALL_H = {k: "H" for k in ACTIVE_H}
REAGENTS = {
    "GC · TMS (MSTFA / BSTFA)": ("TMS", _ALL_H, TMS, False),
    "GC · MeOx + TMS": ("MeOx + TMS", _ALL_H, TMS, False),
    "GC · TBDMS (MTBSTFA)": ("TBDMS", _ALL_H, TBDMS, False),
    "GC · Acylation (TFAA)": ("TFAA", {"alcohol": "atom", "phenol": "atom", "amine1": 1, "amine2": "atom", "thiol": "atom"}, TFAA, False),
    "CIL · Dansyl chloride (amine / phenol)": ("Dansyl chloride", {"amine1": 1, "amine2": "atom", "phenol": "atom", "arom_nh_cil": "atom"}, DANSYL, True),
    "CIL · Dansyl chloride, base-activated (hydroxyl)": ("Dansyl chloride (base)", {"alcohol": "atom", "enol": "atom", "thiol": "atom", "phenol": "atom",
                                                                              "amine1": 1, "amine2": "atom", "arom_nh_cil": "atom"}, DANSYL, True),
    "CIL · DmPA bromide (carboxyl)": ("DmPA bromide", {"acid": "atom"}, DMPA, True),
    "CIL · Dansylhydrazine (carbonyl)": ("Dansylhydrazine", {k: "atom" for k in CARBONYL}, DNS_HYDRAZINE, True),
    "UV/FLD · OPA + 2-mercaptoethanol (primary amine)": ("OPA + 2-ME", {"amine1": 1}, OPA, False),
    "UV/FLD · FMOC-Cl (amine)": ("FMOC-Cl", {"amine1": 1, "amine2": "atom", "phenol": "atom"}, FMOC, False),
    "UV/FLD · DNPH (carbonyl)": ("DNPH", {k: "atom" for k in CARBONYL}, DNPH, False),
}
DEFAULT_REAGENTS = [k for k in REAGENTS if k.startswith("CIL")] + list(REAGENTS)[:2]


def platform_of(option: str) -> str:
    """"GC" / "CIL" / "UV/FLD": the part of the reagent key before the first " · "."""
    return option.split(" · ", 1)[0]


def reagent_label(option: str) -> str:
    """Short label of a reagent, shared by the Derivatization table and the reaction options: "GC · MeOx + TMS", "CIL · Dansyl chloride (base)"."""
    return f"{platform_of(option)} · {REAGENTS[option][0]}"


# ---------------------------------------------------------------- functional-group reaction schemes
# SCHEMES[reagent name] = {"abbr": what the abbreviations on the right stand for, "steps": [(group key, reagent label, from, to), ...]} in display order.
# Functional groups only, never full product structures. Group keys are the GROUPS keys plus "carbonyl" (aldehyde + ketone + aromatic ketone
# together) and "hemiacetal" (cyclic hemiacetal of a reducing sugar, opened by MeOx). reaction_steps() keeps the steps whose group the reagent
# labels (the REAGENTS rule) and attaches the group counts; one step = one group type. Counts are numbers of groups, not of labels
# (a primary amine counts once but takes two TMS: "R–NH₂ → R–N(TMS)₂").
def _silyl(tag):
    return [("acid", "R–COOH", f"R–COO–{tag}"), ("phenol", "Ar–OH", f"Ar–O–{tag}"), ("alcohol", "R–OH", f"R–O–{tag}"),
            ("enol", "C=C–OH", f"C=C–O–{tag}"), ("amine1", "R–NH₂", f"R–N({tag})₂"), ("amine2", "R₂NH", f"R₂N–{tag}"),
            ("amide", "R–CO–NH–", f"R–CO–N({tag})–"), ("arom_nh", "Ar–NH", f"Ar–N–{tag}"), ("thiol", "R–SH", f"R–S–{tag}")]


def _with(reagent, rows):
    return [(key, reagent, frm, to) for key, frm, to in rows]


_DNS = [("amine1", "R–NH₂", "R–NH–Dns"), ("amine2", "R₂NH", "R₂N–Dns"), ("phenol", "Ar–OH", "Ar–O–Dns"), ("arom_nh_cil", "Ar–NH", "Ar–N–Dns")]
_DNS_BASE = [("alcohol", "R–OH", "R–O–Dns"), ("enol", "C=C–OH", "C=C–O–Dns"), ("thiol", "R–SH", "R–S–Dns")] + _DNS
_DNS_DEF = "Dns = dansyl, 5-(dimethylamino)naphthalene-1-sulfonyl"
_MEOX = [("carbonyl", "MeOx", "R₂C=O", "R₂C=N–OCH₃"), ("hemiacetal", "MeOx", "ring hemiacetal", "R–CH=N–OCH₃")]

SCHEMES = {
    "TMS": {"abbr": "TMS = –Si(CH₃)₃", "steps": _with("MSTFA / BSTFA", _silyl("TMS"))},
    "MeOx + TMS": {"abbr": "MeOx = methoxyamine (CH₃O–NH₂), TMS = –Si(CH₃)₃; a reducing sugar's ring opens",
                   "steps": _MEOX + _with("MSTFA / BSTFA", _silyl("TMS"))},
    "TBDMS": {"abbr": "TBDMS = –Si(CH₃)₂C(CH₃)₃", "steps": _with("MTBSTFA", _silyl("TBDMS"))},
    "TFAA": {"abbr": "COCF₃ = trifluoroacetyl",
             "steps": _with("TFAA", [("phenol", "Ar–OH", "Ar–O–COCF₃"), ("alcohol", "R–OH", "R–O–COCF₃"), ("amine1", "R–NH₂", "R–NH–COCF₃"),
                                     ("amine2", "R₂NH", "R₂N–COCF₃"), ("thiol", "R–SH", "R–S–COCF₃")])},
    "Dansyl chloride": {"abbr": _DNS_DEF, "steps": _with("Dns-Cl", _DNS)},
    "Dansyl chloride (base)": {"abbr": _DNS_DEF + "; the base also lets it label amines and phenols", "steps": _with("Dns-Cl + base", _DNS_BASE)},
    "DmPA bromide": {"abbr": "DmPA = –CH₂–CO–C₆H₄–N(CH₃)₂ (4-dimethylaminophenacyl ester)", "steps": _with("DmPA-Br", [("acid", "R–COOH", "R–COO–DmPA")])},
    "Dansylhydrazine": {"abbr": _DNS_DEF, "steps": _with("Dns-NHNH₂", [("carbonyl", "R₂C=O", "R₂C=N–NH–Dns")])},
    "OPA + 2-ME": {"abbr": "OPA = o-phthalaldehyde, 2-ME = 2-mercaptoethanol; the product is a fluorescent isoindole",
                   "steps": _with("OPA + 2-ME", [("amine1", "R–NH₂", "R–N–isoindole")])},
    "FMOC-Cl": {"abbr": "Fmoc = 9-fluorenylmethoxycarbonyl",
                "steps": _with("FMOC-Cl", [("amine1", "R–NH₂", "R–NH–Fmoc"), ("amine2", "R₂NH", "R₂N–Fmoc"), ("phenol", "Ar–OH", "Ar–O–Fmoc")])},
    "DNPH": {"abbr": "DNP = 2,4-dinitrophenyl", "steps": _with("DNPH", [("carbonyl", "R₂C=O", "R₂C=N–NH–DNP")])},
}
assert set(SCHEMES) == {v[0] for v in REAGENTS.values()}

HIGHLIGHT = _theme.HL_REACTIVE


# ---------------------------------------------------------------- 分子式 / 质量
def parent_formula(mol) -> Counter:
    mh = Chem.AddHs(mol)
    return Counter(a.GetSymbol() for a in mh.GetAtoms())


def formula_string(f) -> str:
    f = {k: v for k, v in f.items() if v}
    order = ([e for e in ("C", "H") if e in f] + sorted(e for e in f if e not in ("C", "H"))) if "C" in f else sorted(f)
    return "".join(e + (str(f[e]) if f[e] != 1 else "") for e in order)


def formula_mass(f) -> float:
    return sum(PT.GetMostCommonIsotopeMass(e) * n for e, n in f.items())


def parse_formula(s: str) -> Counter:
    return Counter({e: int(n or 1) for e, n in re.findall(r"([A-Z][a-z]?)(\d*)", s)})


def _add(f, delta, n):
    out = Counter(f)
    for e, v in delta.items():
        out[e] += v * n
    return out


# ---------------------------------------------------------------- 基团
def neutralize(mol):
    """去掉羧酸盐、两性离子等可中和的电荷（季铵等不可中和的保留）-> (中性分子, 是否改动)。"""
    neutral = rdMolStandardize.Uncharger().uncharge(mol)
    return neutral, Chem.MolToSmiles(neutral) != Chem.MolToSmiles(mol)


def find_groups(mol) -> Dict[str, dict]:
    """key -> {"atoms": [0 基序号], "h": 活泼氢总数}（只含 count > 0 的基团）。"""
    found = {}
    taken = set()          # 烯醇 OH 不和酸 / 酚 / 醇重复计
    for key, patt in _PATTERNS.items():
        idx = GROUPS[key][2]
        atoms = sorted({m[idx] for m in mol.GetSubstructMatches(patt)} - (taken if key == "enol" else set()))
        if key in ("acid", "phenol", "alcohol"):
            taken.update(atoms)
        if atoms:
            h = sum(mol.GetAtomWithIdx(i).GetTotalNumHs() for i in atoms)
            found[key] = {"atoms": atoms, "h": h}
    return found


def group_text(groups) -> List[str]:
    """["alcohol OH ×5", "primary amine ×1", ...] in display order: every reactive group found."""
    return [f"{GROUPS[key][0]} ×{len(groups[key]['atoms'])}" for key in _HEAD if key in groups]


def reacting_groups(results, groups) -> Dict[str, dict]:
    """The part of ``groups`` that at least one of the derivatization ``results`` labels: what a structure picture should highlight."""
    used = set()
    for r in results:
        if r["n"]:
            used.update(r["labels"])
            if r["name"] == "MeOx + TMS":
                used.update(CARBONYL)
    return {k: g for k, g in groups.items() if k in used}


def reaction_steps(name, rule, groups, n_hemi=0) -> List[tuple]:
    """[(from, reagent, to, count)]: the functional-group scheme of reagent ``name`` for the groups it labels (``rule`` = its REAGENTS rule)."""
    counts = {}
    for key, g in groups.items():
        if key in rule:
            sk = "carbonyl" if key in CARBONYL else key
            counts[sk] = counts.get(sk, 0) + len(g["atoms"])
    if name == "MeOx + TMS":                      # MeOx takes the carbonyls (not part of the TMS rule) and opens ring hemiacetals
        counts["carbonyl"] = sum(len(groups[k]["atoms"]) for k in CARBONYL if k in groups)
        counts["hemiacetal"] = n_hemi
    return [(frm, reagent, to, counts[key]) for key, reagent, frm, to in SCHEMES[name]["steps"] if counts.get(key)]


def _labels(rule, groups):
    """-> {group key: 标记数}"""
    out = {}
    for key, how in rule.items():
        g = groups.get(key)
        if not g:
            continue
        out[key] = g["h"] if how == "H" else len(g["atoms"]) * (1 if how == "atom" else how)
    return out


def highlight_svg(mol, groups, size=(560, 420)) -> Optional[str]:
    atoms = sorted({i for g in groups.values() for i in g["atoms"]})
    drawer = _depict.svg_drawer(*size)
    drawer.DrawMolecule(mol, highlightAtoms=atoms, highlightAtomColors={i: HIGHLIGHT for i in atoms})
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def _product_sites(result, groups, reagent):
    targets = [i for key in result["labels"] for i in groups[key]["atoms"]]
    if not targets or len(targets) != result["n"]:
        raise ValueError(f"{reagent} product sites do not agree with the label count.")
    return targets


def _without_site_hydrogens(mol, targets, reagent):
    """Explicit Hs make the leaving H unambiguous for N / O / S and aromatic [nH]; heavy-atom indices stay unchanged."""
    work = Chem.RWMol(Chem.AddHs(mol))
    leaving = []
    for i in targets:
        hs = [n.GetIdx() for n in work.GetAtomWithIdx(i).GetNeighbors() if n.GetAtomicNum() == 1]
        if not hs:
            raise ValueError(f"{reagent} reaction site {i} has no removable hydrogen.")
        leaving.append(hs[0])
    for i in sorted(leaving, reverse=True):
        work.RemoveAtom(i)
    return work.GetMol()


def _attach_at_sites(product, targets, tag, anchor, bond_type):
    for target in targets:
        offset = product.GetNumAtoms()
        product = Chem.CombineMols(product, tag)
        edit = Chem.RWMol(product)
        edit.AddBond(target, offset + anchor, bond_type)
        product = edit.GetMol()
    return product


def _clean_product(product):
    Chem.SanitizeMol(product)
    product = Chem.RemoveHs(product)
    Chem.SanitizeMol(product)
    return product


def _dansyl_chloride_product(mol, result, groups):
    targets = _product_sites(result, groups, "Dansyl chloride")
    product = _without_site_hydrogens(mol, targets, "Dansyl chloride")
    # PubChem CID 11801. Remove Cl, then join every amine / phenol / hydroxyl / thiol / aromatic NH site to sulfur.
    tag = Chem.RWMol(Chem.MolFromSmiles("CN(C)C1=CC=CC2=C1C=CC=C2S(=O)(=O)Cl"))
    cl = next(a.GetIdx() for a in tag.GetAtoms() if a.GetAtomicNum() == 17)
    sulfur = next(a.GetIdx() for a in tag.GetAtoms() if a.GetAtomicNum() == 16)
    tag.RemoveAtom(cl)
    if cl < sulfur:
        sulfur -= 1
    product = _attach_at_sites(product, targets, tag.GetMol(), sulfur, Chem.BondType.SINGLE)
    return _clean_product(product)


def _dmpa_product(mol, result, groups):
    targets = _product_sites(result, groups, "DmPA bromide")
    product = _without_site_hydrogens(mol, targets, "DmPA bromide")
    # PubChem CID 142225. Remove Br, then join every carboxyl O to the phenacyl CH2 carbon, forming the complete ester.
    tag = Chem.RWMol(Chem.MolFromSmiles("CN(C)C1=CC=C(C=C1)C(=O)CBr"))
    bromine = next(a.GetIdx() for a in tag.GetAtoms() if a.GetAtomicNum() == 35)
    anchor = next(n.GetIdx() for n in tag.GetAtomWithIdx(bromine).GetNeighbors() if n.GetAtomicNum() == 6)
    tag.RemoveAtom(bromine)
    if bromine < anchor:
        anchor -= 1
    product = _attach_at_sites(product, targets, tag.GetMol(), anchor, Chem.BondType.SINGLE)
    return _clean_product(product)


def _dansylhydrazine_product(mol, result, groups):
    targets = _product_sites(result, groups, "Dansylhydrazine")
    product = Chem.AddHs(mol)
    carbonyl_o = []
    for i in targets:
        oxygen = [n.GetIdx() for n in product.GetAtomWithIdx(i).GetNeighbors()
                  if n.GetAtomicNum() == 8 and product.GetBondBetweenAtoms(i, n.GetIdx()).GetBondType() == Chem.BondType.DOUBLE]
        if len(oxygen) != 1:
            raise ValueError(f"Dansylhydrazine reaction site {i} is not a carbonyl carbon.")
        carbonyl_o.append(oxygen[0])

    # PubChem CID 94442: Dns-SO2-NH-NH2. Terminal NH2 loses two H; carbonyl O leaves as water after C=N formation.
    tag = Chem.RWMol(Chem.AddHs(Chem.MolFromSmiles("CN(C)C1=CC=CC2=C1C=CC=C2S(=O)(=O)NN")))
    terminal = next(a for a in tag.GetAtoms() if a.GetAtomicNum() == 7
                    and sum(n.GetAtomicNum() == 1 for n in a.GetNeighbors()) == 2
                    and any(n.GetAtomicNum() == 7 for n in a.GetNeighbors()))
    anchor = terminal.GetIdx()
    leaving = [n.GetIdx() for n in terminal.GetNeighbors() if n.GetAtomicNum() == 1]
    for i in sorted(leaving, reverse=True):
        tag.RemoveAtom(i)
    product = _attach_at_sites(product, targets, tag.GetMol(), anchor, Chem.BondType.DOUBLE)
    edit = Chem.RWMol(product)
    for i in sorted(carbonyl_o, reverse=True):
        edit.RemoveAtom(i)
    return _clean_product(edit.GetMol())


def reaction_product(mol, result, groups):
    """Complete product for the four CIL channels (A / H / C / K); None for reagents without a full-structure renderer."""
    if not result["n"]:
        return None
    if result["name"] in ("Dansyl chloride", "Dansyl chloride (base)"):
        return _dansyl_chloride_product(mol, result, groups)
    if result["name"] == "DmPA bromide":
        return _dmpa_product(mol, result, groups)
    if result["name"] == "Dansylhydrazine":
        return _dansylhydrazine_product(mol, result, groups)
    return None


def reaction_product_svg(mol, result, groups, size=(640, 480)) -> Optional[str]:
    """Full-product depiction for the four CIL channels; other reagents keep their functional-group schemes."""
    product = reaction_product(mol, result, groups)
    if product is None:
        return None
    rdCoordGen.AddCoords(product)       # keeps multiple large labels separated
    drawer = _depict.svg_drawer(*size)
    drawer.DrawMolecule(product)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


# ---------------------------------------------------------------- 衍生化
def count_hemiacetal(mol) -> int:
    """环状半缩醛碳（还原糖的异头碳）个数；甲氧胺化会把环打开成 C=N-OCH3。"""
    return len({m[1] for m in mol.GetSubstructMatches(_HEMIACETAL)})


def _notes(name, labels, groups, mol, n_carbonyl, n_h_labels):
    """-> (key EI ions, caveats), each a list of short strings."""
    n, notes = [], []
    if name in ("TMS", "MeOx + TMS", "TBDMS"):
        if name == "TBDMS":
            n.append("M−57 (base); M−15")
        else:
            n.append("m/z 73; M−15" + ("; 147 (≥2 TMS)" if n_h_labels >= 2 else ""))
            if name == "MeOx + TMS" and n_carbonyl:
                n.append("M−31")
        if "amine1" in labels:
            notes.append("primary amine may take only 1 silyl group")
        if any(k in labels for k in ("amide", "arom_nh")):
            notes.append("amide / imide / aromatic NH often derivatizes incompletely")
        if name == "MeOx + TMS" and n_carbonyl:
            notes.append("MeOx gives syn / anti double peaks per carbonyl")
        if n_carbonyl and name != "MeOx + TMS":
            notes.append("carbonyl not protected (use MeOx first)")
    elif name == "TFAA":
        n.append("m/z 69 (CF3+)")
        if "acid" in groups:
            notes.append("COOH not acylated")
    else:
        if mol.HasSubstructMatch(_TERTIARY_AMINE):
            notes.append("tertiary amine not labelled")
        if "amide" in groups:
            notes.append("amide not labelled")
    n_oh = len(groups.get("alcohol", {"atoms": []})["atoms"])
    sugar = (n_carbonyl and n_oh >= 3) or mol.HasSubstructMatch(_HEMIACETAL)
    if sugar and name in ("TMS", "MeOx + TMS", "TBDMS"):
        notes.append("reducing sugar: ring / open-chain forms give several peaks")
    return n, notes


def derive(smiles, reagents=None) -> dict:
    """-> {"groups", "results", "n", "neutralized", "mol"}. ``results`` has one dict per known reagent, in order (``n`` is the list of their label counts):
    reagent (key), name, platform, n (labels), labels, formula / mass / heavy (None when n = 0), text (groups labelled), ions (EI key ions),
    warns (caveats), steps (functional-group scheme, reaction_steps)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("RDKit could not parse the structure.")
    mol, neutralized = neutralize(mol)
    groups = find_groups(mol)
    base = parent_formula(mol)
    n_carbonyl = sum(len(groups[k]["atoms"]) for k in CARBONYL if k in groups)
    n_hemi = count_hemiacetal(mol)
    ns, results = [], []
    for option in (reagents if reagents is not None else DEFAULT_REAGENTS):
        if option not in REAGENTS:
            continue
        name, rule, delta, cil = REAGENTS[option]
        labels = _labels(rule, groups)
        n = sum(labels.values())
        res = {"reagent": option, "name": name, "platform": platform_of(option), "n": n, "labels": labels, "formula": None, "mass": None, "heavy": None,
               "text": "", "ions": [], "warns": [], "steps": []}
        text = "; ".join(f"{GROUPS[k][0]} ×{v}" for k, v in labels.items())
        if name == "MeOx + TMS":
            n_meox = n_carbonyl + n_hemi        # 开环后的活泼氢数不变：异头 OH 变 C=O，环 O 变 OH
            res["n_meox"], res["n_tms"] = n_meox, n
            n_total = n + n_meox
            meox_text = ([f"carbonyl ×{n_carbonyl} (MeOx)"] if n_carbonyl else []) + \
                        ([f"hemiacetal ×{n_hemi} (ring opened, MeOx)"] if n_hemi else [])
            text = "; ".join(meox_text + ([text] if text else []))
            f = _add(_add(base, MEOX, n_meox), delta, n)
            res["n"] = n = n_total
        else:
            f = _add(base, delta, n)
        if n:
            m = formula_mass(f)
            heavy = m + 2 * C13_SHIFT * n if cil else None
            res.update(formula=formula_string(f), mass=m, heavy=heavy, text=text,
                       steps=reaction_steps(name, rule, groups, n_hemi if name == "MeOx + TMS" else 0))
            res["ions"], res["warns"] = _notes(name, labels, groups, mol, n_carbonyl + (n_hemi if name == "MeOx + TMS" else 0), res.get("n_tms", n))
        ns.append(res["n"])
        results.append(res)
    return {"groups": groups, "results": results, "n": ns, "neutralized": neutralized, "mol": mol}


# ---------------------------------------------------------------- Platform 页：CIL / GC-MS
# CIL channels (Zhao 2019): key, short name, long name (reagent in brackets), reagent key, the group keys that open the channel
CHANNELS = [
    ("amine_phenol", "amine / phenol", "amine / phenol (dansyl)", "CIL · Dansyl chloride (amine / phenol)", ("amine1", "amine2", "phenol", "arom_nh_cil")),
    ("hydroxyl", "hydroxyl", "hydroxyl (base-activated dansyl)", "CIL · Dansyl chloride, base-activated (hydroxyl)", ("alcohol", "enol", "thiol")),   # 论文把硫醇归入羟基
    ("carboxyl", "carboxyl", "carboxyl (DmPA)", "CIL · DmPA bromide (carboxyl)", ("acid",)),
    ("carbonyl", "carbonyl", "carbonyl (dansylhydrazine)", "CIL · Dansylhydrazine (carbonyl)", CARBONYL),
]


def cil_assess(groups) -> dict:
    """-> {"cil": bool, "channels": [long names], "short": [short names], "options": [reagent keys]}, the three lists aligned (one entry per channel that applies)."""
    open_ = [c for c in CHANNELS if any(k in groups for k in c[4])]
    return {"cil": bool(open_), "channels": [c[2] for c in open_], "short": [c[1] for c in open_], "options": [c[3] for c in open_]}


def gc_option(mol, groups) -> str:
    """The GC reagent whose product the Platform page shows: MeOx + TMS when there is a carbonyl or a reducing-sugar ring, else TMS."""
    n_meox = sum(len(groups[k]["atoms"]) for k in CARBONYL if k in groups) + count_hemiacetal(mol)
    return "GC · MeOx + TMS" if n_meox else "GC · TMS (MSTFA / BSTFA)"


_CHARGED_OK = Chem.MolFromSmarts("[$([N+](=O)[O-]),$([O-][N+]=O)]")      # 硝基不算离子
_LABILE = [("quaternary ammonium", "[NX4+]"), ("phosphate / phosphonate", "P(=O)(O)O"),
           ("sulfate / sulfonate", "S(=O)(=O)[OX2H1,O-]")]


def gc_assess(mol, groups, n_glycosyl=0) -> dict:
    """mol 应是中和后的分子。-> {"gc": "yes" / "after derivatization" / "no", "volatility", "stability", "reasons"}"""
    reasons, bad = [], []
    for name, sma in _LABILE:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(sma)):
            bad.append(name)
    nitro = {i for m in mol.GetSubstructMatches(_CHARGED_OK) for i in m}
    if any(a.GetFormalCharge() and a.GetIdx() not in nitro for a in mol.GetAtoms()):
        bad.append("formal charge (ionic)")
    if n_glycosyl >= 2:
        bad.append(f"{n_glycosyl} glycosidic sugar rings")
    mw = Descriptors.MolWt(mol)
    if bad:
        why = ", ".join(dict.fromkeys(bad))
        reasons.append(f"Rule: {why} -> not volatile, thermally labile")
        return {"gc": "no", "volatility": "not volatile", "stability": f"labile / ionic: {why}", "reasons": reasons}
    n_h = sum(groups[k]["h"] for k in ACTIVE_H if k in groups)
    n_carbonyl = sum(len(groups[k]["atoms"]) for k in CARBONYL if k in groups)
    ok = "OK"
    if n_h == 0 and not n_carbonyl:
        reasons.append(f"Rule: no active H, no carbonyl, MW {mw:.0f}")
        if mw <= 500:
            return {"gc": "yes", "volatility": "volatile as is", "stability": ok, "reasons": reasons}
        return {"gc": "no", "volatility": f"too heavy for GC (MW {mw:.0f})", "stability": ok, "reasons": reasons}
    if mw <= 250 and n_h <= 1:
        reasons.append(f"Rule: small molecule (MW {mw:.0f}), {n_h} active H -> often run directly; derivatization also works")
        return {"gc": "yes", "volatility": "volatile as is (derivatization optional)", "stability": ok, "reasons": reasons}
    base = parent_formula(mol)
    n_meox = n_carbonyl + count_hemiacetal(mol)         # 还原糖环会被 MeOx 打开
    f = _add(_add(base, MEOX, n_meox), TMS, n_h)
    mass = formula_mass(f)
    label = "MeOx + TMS" if n_meox else "TMS"
    reasons.append(f"Rule: {n_h} active H" + (f", {n_meox} carbonyl / opened hemiacetal" if n_meox else "") + f" -> {label} product {mass:.2f} Da")
    if mass <= 600:
        return {"gc": "after derivatization", "stability": ok, "reasons": reasons,
                "volatility": f"needs derivatization ({label}, {n_h} TMS → {mass:.2f} Da)"}
    return {"gc": "no", "volatility": f"derivative too heavy ({mass:.0f} Da)", "stability": ok, "reasons": reasons}


# ---------------------------------------------------------------- GC derivative structures (EI MS Annotation workflow)
# derive() above only counts labels and gives a formula; the EI workflow needs the derivative STRUCTURE (a SMILES for the NEIMS prediction) and the
# way back (a library entry stored as its derivative -> the original compound). Same sites, same counts as derive(): the formula of the structure
# built here equals derive()'s formula for the same reagent (asserted in tests/test_derivatization.py).
EI_MODES = {"TMS": "GC · TMS (MSTFA / BSTFA)", "MeOx + TMS": "GC · MeOx + TMS", "TBDMS": "GC · TBDMS (MTBSTFA)"}
_SILYL_TAG = {"TMS": "C[Si](C)C", "MeOx + TMS": "C[Si](C)C", "TBDMS": "C[Si](C)C(C)(C)C"}      # open valence on Si: the bond to the O / N / S site
_OXIME_ETHER = Chem.MolFromSmarts("[#6X3]=[NX2]-[OX2]-[CH3]")                                    # C=N-OCH3, what MeOx makes of C=O
_HETERO = (7, 8, 16)


def _bond_neighbors(atom, atomic_num, bond_type=None):
    """Indices of the neighbours of ``atom`` with this atomic number (and, if given, joined by this bond type)."""
    mol = atom.GetOwningMol()
    return [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == atomic_num
            and (bond_type is None or mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondType() == bond_type)]


def _meox_edit(work, groups, mol, sites):
    """MeOx on the editable molecule ``work`` (heavy-atom indices as in ``mol``): every C=O becomes C=N-OCH3 and every cyclic hemiacetal (the
    anomeric carbon of a reducing sugar) is opened to C=N-OCH3, its ring oxygen becoming an OH. ``sites`` ({atom: silyl count}) and the
    hydrogens to delete are updated: the anomeric OH is gone, the ring oxygen is a new silylation site. -> (list of H atoms to delete, new sites)"""
    delete_h, new_sites = [], {}
    carbonyl_c = [i for key in CARBONYL if key in groups for i in groups[key]["atoms"]]
    kekulize = any(work.GetAtomWithIdx(i).GetIsAromatic() for i in carbonyl_c)
    if kekulize:
        Chem.Kekulize(work, clearAromaticFlags=True)
    targets = []                                                                    # (nitrogen atom) that gets the O-CH3
    for c in carbonyl_c:
        oxygen = _bond_neighbors(work.GetAtomWithIdx(c), 8, Chem.BondType.DOUBLE)
        if len(oxygen) != 1:
            raise ValueError(f"MeOx reaction site {c} is not a carbonyl carbon.")
        work.GetAtomWithIdx(oxygen[0]).SetAtomicNum(7)
        targets.append(oxygen[0])
    for match in mol.GetSubstructMatches(_HEMIACETAL):
        oh, c, ring_o = match[0], match[1], match[2]
        if oh not in sites:
            continue
        delete_h += _bond_neighbors(work.GetAtomWithIdx(oh), 1)
        work.GetAtomWithIdx(oh).SetAtomicNum(7)
        work.GetBondBetweenAtoms(c, oh).SetBondType(Chem.BondType.DOUBLE)
        work.RemoveBond(c, ring_o)
        targets.append(oh)
        new_sites[ring_o] = 1
        sites.pop(oh)
    for n in targets:
        o_new, c_new = work.AddAtom(Chem.Atom(8)), work.AddAtom(Chem.Atom(6))
        work.AddBond(n, o_new, Chem.BondType.SINGLE)
        work.AddBond(o_new, c_new, Chem.BondType.SINGLE)
    return delete_h, new_sites


def derivatize(smiles, reagent) -> dict:
    """The derivative of ``smiles`` with a GC reagent (a key of EI_MODES.values(), e.g. "GC · TMS (MSTFA / BSTFA)").

    -> {"smiles", "formula", "mass" (monoisotopic), "n_silyl", "n_meox", "n_labels" (what derive() reports as ``n``: silyl groups, plus MeOx for
    MeOx + TMS), "neutralized"}. Sites and counts are derive()'s: every active-H atom of ACTIVE_H takes as many silyl groups as it has hydrogens
    (a primary amine N(TMS)2, a primary amide too); MeOx first turns every C=O into C=N-OCH3 and opens the ring of a reducing sugar (the open-chain
    oxime is the product, as derive() assumes); the silyl groups go on afterwards. No site: the (neutralised) parent comes back unchanged.
    ValueError: unparsable SMILES or an unknown reagent."""
    mode = next((m for m, option in EI_MODES.items() if option == reagent), None)
    if mode is None:
        raise ValueError(f"{reagent!r} is not a GC silylation reagent; use one of {', '.join(EI_MODES.values())}.")
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        raise ValueError("RDKit could not parse the structure.")
    with BlockLogs():                                                               # the Uncharger logs "Running Uncharger" at INFO
        mol, neutralized = neutralize(mol)
    groups = find_groups(mol)
    rule = REAGENTS[reagent][1]
    labels = _labels(rule, groups)
    sites = {}                                                                      # atom -> silyl groups it takes (= its hydrogens)
    for key in labels:
        for i in groups[key]["atoms"]:
            sites[i] = mol.GetAtomWithIdx(i).GetTotalNumHs()
    if sum(sites.values()) != sum(labels.values()):
        raise ValueError(f"{reagent}: the silylation sites do not agree with the label count.")
    work = Chem.RWMol(Chem.AddHs(mol))
    delete_h, n_meox = [], 0
    attach = dict(sites)
    if mode == "MeOx + TMS":
        n_meox = sum(len(groups[k]["atoms"]) for k in CARBONYL if k in groups) + count_hemiacetal(mol)
        extra_h, new_sites = _meox_edit(work, groups, mol, attach)
        delete_h += extra_h
        attach.update(new_sites)
        sites = {i: n for i, n in sites.items() if i in attach}                     # the anomeric OH is a nitrogen now: nothing to delete there
    for i, n in sites.items():
        hydrogens = _bond_neighbors(work.GetAtomWithIdx(i), 1)
        if len(hydrogens) < n:
            raise ValueError(f"{reagent} reaction site {i} has no removable hydrogen.")
        delete_h += hydrogens[:n]
    tag = Chem.MolFromSmiles(_SILYL_TAG[mode])
    anchor = next(a.GetIdx() for a in tag.GetAtoms() if a.GetAtomicNum() == 14)
    targets = [i for i, n in attach.items() for _ in range(n)]
    product = work.GetMol()
    if targets:
        for i in sorted(set(delete_h), reverse=True):                               # hydrogens come last in the atom list: heavy-atom indices stay valid
            edit = Chem.RWMol(product)
            edit.RemoveAtom(i)
            product = edit.GetMol()
        product = _attach_at_sites(product, targets, tag, anchor, Chem.BondType.SINGLE)
    elif delete_h:
        raise ValueError("Hydrogens were removed without a silyl group to replace them.")
    product = _clean_product(product)
    formula = parent_formula(product)
    return {"smiles": Chem.MolToSmiles(product), "formula": formula_string(formula), "mass": formula_mass(formula), "n_silyl": len(targets),
            "n_meox": n_meox, "n_labels": len(targets) + n_meox, "neutralized": neutralized}


def derivative_smiles(smiles, reagent):
    """-> (derivative SMILES, number of labels): see derivatize()."""
    d = derivatize(smiles, reagent)
    return d["smiles"], d["n_labels"]


def _silyl_groups(mol):
    """[(Si atom, 'TMS' | 'TBDMS', heteroatom)] for every Si(CH3)3 / Si(CH3)2C(CH3)3 that sits on an O, N or S."""
    found = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 14 or atom.GetDegree() != 4:
            continue
        hetero = [n for n in atom.GetNeighbors() if n.GetAtomicNum() in _HETERO]
        carbons = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(hetero) != 1 or len(carbons) != 3:
            continue
        if not any(n.GetAtomicNum() not in (1, 14) for n in hetero[0].GetNeighbors()):
            continue                                                                # Si-O-Si, Si-N-Si, Si-OH: a siloxane, silazane or silanol (reagent, column bleed), not a derivative of water or ammonia
        methyl = [c for c in carbons if c.GetDegree() == 1]
        tert = [c for c in carbons if c.GetDegree() == 4 and sum(1 for n in c.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetDegree() == 1) == 3]
        if len(methyl) == 3:
            found.append((atom, "TMS", hetero[0]))
        elif len(methyl) == 2 and len(tert) == 1:
            found.append((atom, "TBDMS", hetero[0]))
    return found


def underivatize_info(smiles) -> dict:
    """Undo the GC derivatization of a SMILES -> {"parent" (SMILES), "n_tms", "n_tbdms", "n_meox", "mode" (None | "TMS" | "TBDMS" | "MeOx" |
    "MeOx + TMS" | "MeOx + TBDMS"), "changed"}. Every O/N/S-Si(CH3)3 and O/N/S-Si(CH3)2C(CH3)3 becomes O/N/S-H, every C=N-OCH3 becomes C=O. A cyclic sugar
    that MeOx opened comes back as its open chain (a different InChIKey from the ring form). A plain organosilane (Si-C bonds only) is not a
    derivative and stays as it is. ValueError: unparsable SMILES."""
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        raise ValueError("RDKit could not parse the structure.")
    silyl = _silyl_groups(mol)
    oximes = mol.GetSubstructMatches(_OXIME_ETHER)
    work = Chem.RWMol(mol)
    drop = set()
    for si, kind, hetero in silyl:
        drop.add(si.GetIdx())
        for c in si.GetNeighbors():
            if c.GetAtomicNum() == 6:
                drop.add(c.GetIdx())
                if c.GetDegree() == 4:                                              # the tert-butyl carbon takes its methyls with it
                    drop.update(n.GetIdx() for n in c.GetNeighbors() if n.GetIdx() != si.GetIdx())
        h = work.GetAtomWithIdx(hetero.GetIdx())
        h.SetNumExplicitHs(h.GetNumExplicitHs() + 1)
    for c, n, o, me in oximes:
        work.GetAtomWithIdx(n).SetAtomicNum(8)
        drop.update((o, me))
    for i in sorted(drop, reverse=True):
        work.RemoveAtom(i)
    parent = work.GetMol()
    Chem.SanitizeMol(parent)
    n_tms, n_tbdms, n_meox = sum(1 for s in silyl if s[1] == "TMS"), sum(1 for s in silyl if s[1] == "TBDMS"), len(oximes)
    kind = "TBDMS" if n_tbdms else "TMS" if n_tms else None
    mode = " + ".join(x for x in ("MeOx" if n_meox else None, kind) if x) or None
    return {"parent": Chem.MolToSmiles(parent), "n_tms": n_tms, "n_tbdms": n_tbdms, "n_meox": n_meox, "mode": mode,
            "changed": bool(silyl or oximes)}


def underivatize(smiles) -> str:
    """The original compound behind a GC derivative SMILES (see underivatize_info())."""
    return underivatize_info(smiles)["parent"]


def nominal_mass(formula) -> int:
    """Nominal mass of a formula ("C11H29NO2Si3" or a Counter): the sum of the integer masses of the most abundant isotopes (what an EI molecular ion is quoted as)."""
    counts = parse_formula(formula) if isinstance(formula, str) else formula
    return sum(PT.GetMostCommonIsotope(element) * n for element, n in counts.items())


def carbonyl_count(smiles) -> int:
    """C=O groups MeOx would turn into C=N-OCH3 (aldehydes, ketones, and the cyclic hemiacetal carbon of a reducing sugar); 0 for an unparsable SMILES."""
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return 0
    with BlockLogs():
        mol, _ = neutralize(mol)
    groups = find_groups(mol)
    return sum(len(groups[k]["atoms"]) for k in CARBONYL if k in groups) + count_hemiacetal(mol)
