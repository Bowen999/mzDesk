"""lipid_metabolite_classifier.py

Classify a molecule (given as SMILES) as detected by lipidomics, metabolomics,
both, or neither, from RDKit-computed structural features and a tunable
heuristic. Scientific rationale documented in each detector.
"""
from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen


# --------------------------------------------------------------------------
# Tunable thresholds. Every heuristic cut-off lives here so users can retune.
# --------------------------------------------------------------------------
@dataclass
class Thresholds:
    mw_min: float = 50.0            # below -> outside typical LC-MS window
    mw_max: float = 2000.0          # above -> outside typical LC-MS window
    acyl_min_carbons: int = 12      # "long" fatty/acyl chain length (>=C12)
    logp_polar_max: float = 3.0     # logP below this (no lipid scaffold) -> polar/metabolite
    medium_chain_max: int = 14      # C12-C14 counts as "medium/boundary" chain
    boundary_logp_low: float = 1.0  # boundary FFA logP lower bound (both possible)
    boundary_logp_high: float = 3.0 # boundary FFA logP upper bound
    boundary_mw_max: float = 400.0  # boundary FFA size ceiling
    isoprene_min_units: int = 2     # >=2 isoprene units -> isoprenoid chain


# SMARTS for polar head groups, ordered specific -> generic. Matches are claimed
# greedily so a phosphate inside phosphocholine is not double-counted.
POLAR_HEAD_SMARTS: List[Tuple[str, str]] = [
    ("phosphocholine",    "OP(=O)([O-,OX2H,OX1])OCC[N+](C)(C)C"),
    ("phosphoethanolamine","OP(=O)([O-,OX2H,OX1])OCC[NX3;!$(NC=O)]"),
    ("phosphoserine",     "OP(=O)([O-,OX2H,OX1])OCC([NX3])C(=O)[O-,OX2H,OX1]"),
    ("phosphoinositol",   "OP(=O)([O-,OX2H,OX1])OC1C(O)C(O)C(O)C(O)C1O"),
    ("phosphoglycerol",   "OP(=O)([O-,OX2H,OX1])OCC(O)CO"),
    ("phosphate",         "P(=O)([O-,OX2H,OX1])([O-,OX2H,OX1])[O-,OX2H,OX1]"),
    ("sulfate/sulfo",     "S(=O)(=O)[O-,OX2H,OX1]"),
    ("glycerol",          "[OX2H][CH2;!R][CH;!R]([OX2H])[CH2;!R][OX2H]"),
    ("quaternary_ammonium","[NX4+;!$([NX4+]~[!#6])]"),  # choline/carnitine style
]

# Isoprenoid detection (see detect_isoprenoid). A chain is treated as isoprenoid
# if it carries >= ISOPRENE_ALKENE_METHYL_MIN "methyl-on-alkene" motifs
# (unsaturated terpenes/prenols/retinoids/carotenoids) OR >= 1 saturated 2-unit
# isoprene motif (phytanyl / dolichol / tocopherol tail). Empirically validated
# against a labeled positive/negative panel.
ISOPRENE_ALKENE_METHYL_SMARTS = "[CH3][CX3]=[CX3]"   # a methyl on a C=C carbon
ISOPRENE_SATURATED_2UNIT_SMARTS = "[CH3][CHX4]([CX4])[CH2][CH2][CH2][CHX4][CH3]"
ISOPRENE_ALKENE_METHYL_MIN = 2                       # >= this -> isoprenoid


@dataclass
class Features:
    """All structural features extracted from one molecule."""
    smiles: str
    valid: bool
    mw: float = 0.0
    logp: float = 0.0
    has_carbon: bool = False
    has_long_acyl: bool = False          # has a >=C(acyl_min_carbons) aliphatic chain
    longest_chain_carbons: int = 0       # longest contiguous non-ring aliphatic C chain
    is_acyl_linked: bool = False         # that long chain sits on an ester/amide/acid carbonyl
    has_isoprenoid: bool = False         # sterol nucleus OR >= isoprene_min_units isoprene units
    isoprenoid_type: Optional[str] = None  # 'sterol' | 'terpenoid/prenol' | None
    n_polar_head: int = 0
    polar_head_types: List[str] = field(default_factory=list)
    n_glycosyl: int = 0
    # atom-index lists for highlighting in visualization
    highlight: Dict[str, List[int]] = field(default_factory=dict)  # keys: acyl,isoprenoid,head,glycosyl


@dataclass
class Result:
    category: str            # '脂质组' | '代谢组' | '都有可能' | '都不可能' | '无效SMILES'
    reasons: List[str]       # human-readable rule hits (Chinese)
    features: Features


# --------------------------------------------------------------------------
# Feature detectors
# --------------------------------------------------------------------------
def parse_smiles(smiles: str) -> Optional[Chem.Mol]:
    """Return an RDKit Mol or None if the SMILES is invalid/empty."""
    if smiles is None:
        return None
    s = smiles.strip()
    if not s:
        return None
    try:
        mol = Chem.MolFromSmiles(s)
    except Exception:
        return None
    return mol


def longest_aliphatic_carbon_chain(mol: Chem.Mol) -> Tuple[int, List[int]]:
    """Longest contiguous chain of non-aromatic, non-ring carbon atoms.

    Returns (num_carbons_on_longest_path, atom_indices_on_that_path).
    ALGORITHM:
      - eligible atoms: GetAtomicNum()==6 and not GetIsAromatic() and not IsInRing().
      - build adjacency among eligible atoms (bond between two eligible atoms).
      - eligible subgraph is a forest (ring atoms excluded); for each connected
        component find the longest simple path (tree diameter) via double-BFS,
        returning the actual atom indices on that path.
      - return the global maximum over components (count = number of atoms).
      - empty -> (0, []).
    """
    if mol is None:
        return (0, [])

    eligible = []
    eligible_set = set()
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic() and not atom.IsInRing():
            eligible.append(idx)
            eligible_set.add(idx)

    if not eligible:
        return (0, [])

    adj = {i: [] for i in eligible}
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if i in eligible_set and j in eligible_set:
            adj[i].append(j)
            adj[j].append(i)

    def bfs(start: int) -> Tuple[List[int], set]:
        parent = {start: None}
        dist = {start: 0}
        q = deque([start])
        seen = {start}
        farthest = start

        while q:
            cur = q.popleft()
            for nb in adj[cur]:
                if nb not in seen:
                    seen.add(nb)
                    parent[nb] = cur
                    dist[nb] = dist[cur] + 1
                    if dist[nb] > dist[farthest]:
                        farthest = nb
                    q.append(nb)

        parent2 = {farthest: None}
        dist2 = {farthest: 0}
        q2 = deque([farthest])
        seen2 = {farthest}
        farthest2 = farthest

        while q2:
            cur = q2.popleft()
            for nb in adj[cur]:
                if nb not in seen2:
                    seen2.add(nb)
                    parent2[nb] = cur
                    dist2[nb] = dist2[cur] + 1
                    if dist2[nb] > dist2[farthest2]:
                        farthest2 = nb
                    q2.append(nb)

        path = []
        node = farthest2
        while node is not None:
            path.append(node)
            node = parent2[node]
        return path, seen

    best_path: List[int] = []
    visited = set()
    for start in eligible:
        if start in visited:
            continue
        path, component = bfs(start)
        visited.update(component)
        if len(path) > len(best_path):
            best_path = path

    return len(best_path), best_path


def detect_long_acyl(mol: Chem.Mol, min_carbons: int = 12) -> Tuple[bool, int, bool, List[int]]:
    """Detect a long (>=min_carbons) fatty/acyl chain.

    Returns (has_long, longest_carbons, is_acyl_linked, atom_indices).
    ALGORITHM:
      - (longest_carbons, path) = longest_aliphatic_carbon_chain(mol).
      - has_long = longest_carbons >= min_carbons.
      - is_acyl_linked: any atom on `path` is a carbonyl carbon
        C(=O) that is part of ester/amide/thioester/carboxylic acid, i.e. matches
        one of: "[CX3](=O)[OX2H0][#6]", "[CX3](=O)[OX2H]", "[CX3](=O)[NX3]",
        "[CX3](=O)[SX2]". (Substructure match atom sets intersect path.)
      - atom_indices = path if has_long else [].
    """
    if mol is None:
        return (False, 0, False, [])

    longest_carbons, path = longest_aliphatic_carbon_chain(mol)
    has_long = longest_carbons >= min_carbons
    is_acyl = False

    if has_long:
        path_set = set(path)
        acyl_smarts = [
            "[CX3](=O)[OX2H0][#6]",
            "[CX3](=O)[OX2H]",
            "[CX3](=O)[NX3]",
            "[CX3](=O)[SX2]",
        ]
        for smarts in acyl_smarts:
            patt = Chem.MolFromSmarts(smarts)
            if patt is None:
                continue
            for match in mol.GetSubstructMatches(patt, uniquify=True):
                if any(idx in path_set for idx in match):
                    is_acyl = True
                    break
            if is_acyl:
                break

    return has_long, longest_carbons, is_acyl, (path if has_long else [])


def detect_isoprenoid(mol: Chem.Mol, min_units: int = 2) -> Tuple[bool, Optional[str], List[int]]:
    """Detect a sterol nucleus or an isoprenoid (terpenoid/prenol) chain.

    Returns (has_isoprenoid, type_str, atom_indices).
    ALGORITHM:
      1) STEROL NUCLEUS via ring fusion analysis:
         - ri = mol.GetRingInfo(); rings = ri.AtomRings().
         - build fused-ring graph: two rings are fused if they share >=2 atoms.
         - for each connected component of fused rings, take union of atoms and
           the multiset of ring sizes.
         - a steroid nucleus = component with >=4 rings whose size-multiset
           contains three 6-membered rings and one 5-membered ring, whose atom
           union has 16..18 atoms, and where >=15 of those atoms are carbon.
         - if found -> return (True, 'sterol', sorted(atom_union)).
      2) ISOPRENOID CHAIN (empirically validated against a labeled panel):
         - count matches of ISOPRENE_ALKENE_METHYL_SMARTS (a methyl on a C=C carbon)
           -- these tally isoprene units in unsaturated terpenes/prenols/retinoids/
           carotenoids.
         - count matches of ISOPRENE_SATURATED_2UNIT_SMARTS (a saturated 2-unit
           isoprene motif) -- for phytanyl/dolichol/tocopherol-type saturated tails.
         - isoprenoid if (alkene-methyls >= max(min_units, 1)) OR (saturated motif >= 1).
         - return (True, 'terpenoid/prenol', sorted(union of matched atoms)).
      3) else (False, None, []).
    """
    if mol is None:
        return (False, None, [])

    ri = mol.GetRingInfo()
    atom_rings = [set(r) for r in ri.AtomRings()]

    if atom_rings:
        n = len(atom_rings)
        ring_adj = [[] for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                if len(atom_rings[i] & atom_rings[j]) >= 2:
                    ring_adj[i].append(j)
                    ring_adj[j].append(i)

        visited_rings = [False] * n
        for i in range(n):
            if visited_rings[i]:
                continue
            component = []
            stack = [i]
            visited_rings[i] = True
            while stack:
                cur = stack.pop()
                component.append(cur)
                for nb in ring_adj[cur]:
                    if not visited_rings[nb]:
                        visited_rings[nb] = True
                        stack.append(nb)

            atom_union = set()
            ring_sizes = Counter()
            for ring_idx in component:
                atom_union |= atom_rings[ring_idx]
                ring_sizes[len(atom_rings[ring_idx])] += 1

            if len(component) >= 4 and ring_sizes[6] >= 3 and ring_sizes[5] >= 1:
                if 16 <= len(atom_union) <= 18:
                    carbon_count = sum(
                        1 for idx in atom_union if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6
                    )
                    if carbon_count >= 15:
                        return True, 'sterol', sorted(atom_union)

    # 2) Isoprenoid chain: count methyl-on-alkene motifs (unsaturated terpenes/
    #    prenols/retinoids/carotenoids) and saturated 2-unit isoprene motifs
    #    (phytanyl/dolichol/tocopherol). Either signal (enough alkene methyls, or
    #    one saturated 2-unit motif) marks the molecule as isoprenoid.
    matched_atoms = set()

    alkene_methyl = Chem.MolFromSmarts(ISOPRENE_ALKENE_METHYL_SMARTS)
    n_alkene_methyl = 0
    if alkene_methyl is not None:
        matches = mol.GetSubstructMatches(alkene_methyl, uniquify=True)
        n_alkene_methyl = len(matches)
        for match in matches:
            matched_atoms.update(match)

    sat_2unit = Chem.MolFromSmarts(ISOPRENE_SATURATED_2UNIT_SMARTS)
    n_sat = 0
    if sat_2unit is not None:
        matches = mol.GetSubstructMatches(sat_2unit, uniquify=True)
        n_sat = len(matches)
        for match in matches:
            matched_atoms.update(match)

    if n_alkene_methyl >= max(min_units, 1) or n_sat >= 1:
        return True, 'terpenoid/prenol', sorted(matched_atoms)

    return False, None, []


def detect_polar_heads(mol: Chem.Mol) -> Tuple[int, List[str], List[int]]:
    """Count polar head groups using POLAR_HEAD_SMARTS with greedy claiming.

    Returns (count, type_names, atom_indices).
    ALGORITHM:
      - claimed = set(); count=0; names=[]; atoms=set().
      - iterate POLAR_HEAD_SMARTS in order; for each, GetSubstructMatches(uniquify=True);
        for each match: if match atoms are disjoint from `claimed`, claim them,
        count += 1, append name, atoms |= match.
      - return count, names, sorted(atoms).
    """
    if mol is None:
        return (0, [], [])

    claimed = set()
    count = 0
    names = []
    atoms = set()

    for name, smarts in POLAR_HEAD_SMARTS:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue
        for match in mol.GetSubstructMatches(patt, uniquify=True):
            match_set = set(match)
            if match_set.isdisjoint(claimed):
                claimed.update(match_set)
                count += 1
                names.append(name)
                atoms.update(match_set)

    return count, names, sorted(atoms)


def detect_glycosyl(mol: Chem.Mol) -> Tuple[int, List[int]]:
    """Count sugar (glycosyl) rings and return their atoms.

    Returns (count, atom_indices).
    ALGORITHM:
      - for each ring in GetRingInfo().AtomRings():
          size in (5,6); exactly ONE ring atom is oxygen; all other ring atoms
          are carbon; and >=2 ring carbons carry an exocyclic oxygen neighbor
          (an O outside the ring, via single bond) -> count it as a glycosyl ring.
      - return (count, sorted(union of qualifying ring atoms)).
    """
    if mol is None:
        return (0, [])

    ri = mol.GetRingInfo()
    count = 0
    atoms = set()

    for ring in ri.AtomRings():
        if len(ring) not in (5, 6):
            continue
        ring_set = set(ring)

        oxygen_indices = [
            idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8
        ]
        if len(oxygen_indices) != 1:
            continue

        carbon_indices = [
            idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6
        ]
        if len(carbon_indices) != len(ring) - 1:
            continue

        exo_oxygen_carbons = 0
        for cidx in carbon_indices:
            atom = mol.GetAtomWithIdx(cidx)
            has_exo_o = False
            for nbr in atom.GetNeighbors():
                nidx = nbr.GetIdx()
                if nidx in ring_set:
                    continue
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(cidx, nidx)
                    if bond is not None and bond.GetBondType() == Chem.BondType.SINGLE:
                        has_exo_o = True
                        break
            if has_exo_o:
                exo_oxygen_carbons += 1

        if exo_oxygen_carbons >= 2:
            count += 1
            atoms.update(ring_set)

    return count, sorted(atoms)


def extract_features(smiles: str, thr: Optional[Thresholds] = None) -> Features:
    """Run all detectors and assemble a Features object.

    ALGORITHM:
      - thr = thr or Thresholds()
      - mol = parse_smiles(smiles); if None -> Features(smiles, valid=False).
      - mw = Descriptors.MolWt(mol); logp = Crippen.MolLogP(mol).
      - has_carbon = any atom is carbon.
      - (has_long, longest, acyl_linked, acyl_atoms) = detect_long_acyl(mol, thr.acyl_min_carbons)
      - (has_iso, iso_type, iso_atoms) = detect_isoprenoid(mol, thr.isoprene_min_units)
      - (n_head, head_types, head_atoms) = detect_polar_heads(mol)
      - (n_gly, gly_atoms) = detect_glycosyl(mol)
      - fill Features incl. highlight={'acyl':acyl_atoms,'isoprenoid':iso_atoms,
        'head':head_atoms,'glycosyl':gly_atoms}.
    """
    thr = thr or Thresholds()
    mol = parse_smiles(smiles)

    if mol is None:
        return Features(smiles=smiles, valid=False)

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    has_carbon = any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())

    has_long, longest, acyl_linked, acyl_atoms = detect_long_acyl(
        mol, thr.acyl_min_carbons
    )
    has_iso, iso_type, iso_atoms = detect_isoprenoid(mol, thr.isoprene_min_units)
    n_head, head_types, head_atoms = detect_polar_heads(mol)
    n_gly, gly_atoms = detect_glycosyl(mol)

    return Features(
        smiles=smiles,
        valid=True,
        mw=mw,
        logp=logp,
        has_carbon=has_carbon,
        has_long_acyl=has_long,
        longest_chain_carbons=longest,
        is_acyl_linked=acyl_linked,
        has_isoprenoid=has_iso,
        isoprenoid_type=iso_type,
        n_polar_head=n_head,
        polar_head_types=head_types,
        n_glycosyl=n_gly,
        highlight={
            "acyl": acyl_atoms,
            "isoprenoid": iso_atoms,
            "head": head_atoms,
            "glycosyl": gly_atoms,
        },
    )


def classify(features: Features, thr: Optional[Thresholds] = None) -> Result:
    """Apply the tunable decision tree; return a Result with Chinese reasons.

    DECISION TREE (thr defaults in Thresholds):
      - if not features.valid: category='无效SMILES', reason 'SMILES无法解析'.
      - if not features.has_carbon OR mw<mw_min OR mw>mw_max:
            category='都不可能'  (outside typical LC-MS detection window / 非有机小分子)
      - lipid_scaffold = has_long_acyl or has_isoprenoid
      - amphipathic = n_polar_head + n_glycosyl
      - if lipid_scaffold:
            boundary = (has_long_acyl and not has_isoprenoid
                        and longest_chain_carbons <= medium_chain_max
                        and boundary_logp_low <= logp <= boundary_logp_high
                        and mw <= boundary_mw_max
                        and amphipathic == 0)
            category = '都有可能' if boundary else '脂质组'
      - else:  # no lipid scaffold
            category = '代谢组' if logp < logp_polar_max else '都有可能'
      - Build `reasons` listing: logP, MW & window status, long-chain result
        (with carbon count), isoprenoid result (+type), polar head count(+types),
        glycosyl count, and a final one-sentence rationale for the category.
    """
    thr = thr or Thresholds()

    if not features.valid:
        return Result("无效SMILES", ["SMILES无法解析"], features)

    reasons = [f"logP={features.logp:.2f}，分子量={features.mw:.1f} Da"]
    details = []

    if features.has_long_acyl:
        if features.is_acyl_linked:
            linkage_desc = "与酯/酰胺/羧酸等羰基相连"
        else:
            linkage_desc = "未连接酯/酰胺/羧酸羰基"
        details.append(
            f"含≥{thr.acyl_min_carbons}碳长脂肪链"
            f"（最长C{features.longest_chain_carbons}），{linkage_desc}"
        )
    else:
        details.append(
            f"未检测到≥{thr.acyl_min_carbons}碳长脂肪链"
            f"（最长C{features.longest_chain_carbons}）"
        )

    if features.has_isoprenoid:
        details.append(f"检测到类异戊二烯骨架（{features.isoprenoid_type}）")
    else:
        details.append("未检测到类异戊二烯骨架")

    head_str = "、".join(features.polar_head_types) if features.polar_head_types else "无"
    details.append(f"极性头基数量={features.n_polar_head}（{head_str}）")

    if features.n_glycosyl:
        details.append(f"糖基环数量={features.n_glycosyl}（有）")
    else:
        details.append("糖基环数量=0（无）")

    if features.mw < thr.mw_min or features.mw > thr.mw_max:
        reasons.append(
            f"分子量{features.mw:.1f} Da不在典型LC-MS窗口"
            f"（{thr.mw_min}-{thr.mw_max} Da）内"
        )
    else:
        reasons.append(
            f"分子量{features.mw:.1f} Da在典型LC-MS窗口"
            f"（{thr.mw_min}-{thr.mw_max} Da）内"
        )

    if not features.has_carbon:
        reasons.extend(details)
        reasons.append("不含碳原子，非有机小分子；判为都不可能")
        return Result("都不可能", reasons, features)

    if features.mw < thr.mw_min or features.mw > thr.mw_max:
        reasons.extend(details)
        reasons.append("分子量超出典型LC-MS检测窗口；判为都不可能")
        return Result("都不可能", reasons, features)

    reasons.extend(details)

    lipid_scaffold = features.has_long_acyl or features.has_isoprenoid
    amphipathic = features.n_polar_head + features.n_glycosyl

    if lipid_scaffold:
        boundary = (
            features.has_long_acyl
            and not features.has_isoprenoid
            and features.longest_chain_carbons <= thr.medium_chain_max
            and thr.boundary_logp_low <= features.logp <= thr.boundary_logp_high
            and features.mw <= thr.boundary_mw_max
            and amphipathic == 0
        )
        if boundary:
            reasons.append("中等长度脂链、低极性且无极性头基，脂质组/代谢组均可能检出")
            return Result("都有可能", reasons, features)
        reasons.append("存在明确脂质骨架，判为脂质组")
        return Result("脂质组", reasons, features)

    if features.logp < thr.logp_polar_max:
        reasons.append(
            f"无脂质骨架且logP={features.logp:.2f}<{thr.logp_polar_max}，"
            "倾向于极性代谢物，归为代谢组"
        )
        return Result("代谢组", reasons, features)

    reasons.append(
        f"无脂质骨架但logP={features.logp:.2f}≥{thr.logp_polar_max}，"
        "两类检测均可能"
    )
    return Result("都有可能", reasons, features)


def analyze(smiles: str, thr: Optional[Thresholds] = None) -> Result:
    """Top-level convenience: extract_features then classify."""
    return classify(extract_features(smiles, thr), thr)
