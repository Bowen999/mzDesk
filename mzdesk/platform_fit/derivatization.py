"""Derivatization：找反应基团，按试剂算衍生物的分子式和单同位素质量，并给 Platform 页判断 CIL / GC-MS 适用性。

母体取画出来的中性分子。质量全部用精确同位素质量（不用平均质量）。
每个试剂还带一份官能团反应式（SCHEMES）；Dansyl chloride 另生成完整产物结构，Platform 和 Derivatization 两页的 reaction 块共用。
"""
from __future__ import annotations

import re
from collections import Counter
from typing import Dict, List, Optional

from rdkit import Chem
from rdkit.Chem import Descriptors, rdCoordGen
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Draw import rdMolDraw2D

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


def dansyl_product(mol, result, groups):
    """Complete Dansyl chloride product, with one tag attached at every labelled atom; None for other reagents."""
    if result["name"] not in ("Dansyl chloride", "Dansyl chloride (base)"):
        return None
    targets = [i for key in result["labels"] for i in groups[key]["atoms"]]
    if not targets or len(targets) != result["n"]:
        raise ValueError("Dansyl product sites do not agree with the label count.")

    # Explicit Hs make the leaving H unambiguous for N / O / S and aromatic [nH]. Heavy-atom indices stay unchanged.
    work = Chem.RWMol(Chem.AddHs(mol))
    leaving = []
    for i in targets:
        hs = [n.GetIdx() for n in work.GetAtomWithIdx(i).GetNeighbors() if n.GetAtomicNum() == 1]
        if not hs:
            raise ValueError(f"Dansyl reaction site {i} has no removable hydrogen.")
        leaving.append(hs[0])
    for i in sorted(leaving, reverse=True):
        work.RemoveAtom(i)
    product = work.GetMol()

    # PubChem CID 11801: 5-(dimethylamino)naphthalene-1-sulfonyl chloride. Remove Cl, then join each site to sulfur.
    tag = Chem.RWMol(Chem.MolFromSmiles("CN(C)C1=CC=CC2=C1C=CC=C2S(=O)(=O)Cl"))
    cl = next(a.GetIdx() for a in tag.GetAtoms() if a.GetAtomicNum() == 17)
    sulfur = next(a.GetIdx() for a in tag.GetAtoms() if a.GetAtomicNum() == 16)
    tag.RemoveAtom(cl)
    if cl < sulfur:
        sulfur -= 1
    tag_mol = tag.GetMol()
    for target in targets:
        offset = product.GetNumAtoms()
        product = Chem.CombineMols(product, tag_mol)
        edit = Chem.RWMol(product)
        edit.AddBond(target, offset + sulfur, Chem.BondType.SINGLE)
        product = edit.GetMol()
    Chem.SanitizeMol(product)
    product = Chem.RemoveHs(product)
    Chem.SanitizeMol(product)
    return product


def dansyl_product_svg(mol, result, groups, size=(640, 480)) -> Optional[str]:
    """Full-product depiction for Dansyl chloride options; other reagents keep their functional-group schemes."""
    product = dansyl_product(mol, result, groups)
    if product is None:
        return None
    rdCoordGen.AddCoords(product)       # keeps multiple large Dansyl tags separated
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
