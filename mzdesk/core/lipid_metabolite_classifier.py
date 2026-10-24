"""lipid_metabolite_classifier.py

Classify a molecule (given as SMILES) as detected by lipidomics, metabolomics,
both, or neither, from RDKit-computed structural features and a tunable
heuristic. Scientific rationale documented in each detector.

WHAT THE HEURISTIC MODELS
-------------------------
Which platform reports a compound is decided by three things, and the decision
tree below follows them in that order:

  1. LIPID-CLASS MEMBERSHIP (structural, decisive). A lipidomics workflow targets
     the defined LIPID MAPS classes with class-specific precursor/neutral-loss
     scans and class-based annotation libraries, so anything carrying a canonical
     lipid backbone is reported *regardless* of how borderline its
     physicochemistry is (PC 10:0/10:0 and LPC 8:0 are sold as lipid standards).
     `assign_lipid_class` therefore looks for the backbone, not for size.
  2. EXTRACTION PARTITIONING (physicochemical). Lipidomics keeps the organic
     phase of a biphasic extraction, metabolomics the aqueous/methanolic phase,
     so logP separates the two for everything without a class backbone.
  3. THE INSTRUMENT WINDOW. Outside ~50-2000 Da neither platform reports it.

Every cut-off lives in `Thresholds` so the whole tree can be retuned.
"""
from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen


# --------------------------------------------------------------------------
# Tunable thresholds. Every heuristic cut-off lives here so users can retune.
# --------------------------------------------------------------------------
@dataclass
class Thresholds:
    mw_min: float = 50.0            # below -> outside typical LC-MS window
    mw_max: float = 2000.0          # above -> outside typical LC-MS window
    acyl_min_carbons: int = 12      # "long" apolar fatty/acyl tail (>=C12)
    logp_polar_max: float = 3.0     # logP below this (no lipid class) -> polar/metabolite
    medium_chain_max: int = 14      # C12-C14 counts as "medium/boundary" chain
    boundary_logp_low: float = 1.0  # boundary FFA logP lower bound (both possible)
    boundary_logp_high: float = 5.0 # boundary FFA logP upper bound (covers C12-C14 FFAs)
    boundary_mw_max: float = 400.0  # boundary FFA size ceiling
    isoprene_min_units: int = 2     # >=2 isoprene units -> isoprenoid chain
    # A lipid tail is a *mostly* hydrocarbon run. Hydroxy/oxo/epoxy fatty acids
    # (HETEs, leukotrienes, resolvins) carry a few oxygens and are still lipids,
    # so the test is the DENSITY of substituted carbons, not their absence -- at
    # 0.3 a C20 oxylipin with three oxygenated carbons passes while a fully
    # hydroxylated C16 polyol (logP ~ -10) does not.
    chain_max_polar_fraction: float = 0.3
    fatty_acid_max_carbons: int = 40  # above this a carboxylic acid is not a fatty acid
    # Terpenes, prenols, retinoids, quinones, carotenoids, waxes and fatty
    # alcohols are all hydrocarbon-rich, and those two classes are assigned on a
    # greasy motif alone, so both have to stay heteroatom-poor. Capping
    # heteroatoms per carbon stops the methyl-on-alkene motif from claiming
    # tetrapyrroles (bilirubin) and keeps short prenyl pyrophosphates (FPP,
    # GGPP) -- which are mevalonate-pathway intermediates, not membrane lipids --
    # out of the lipid classes, while long polyprenyl phosphates still qualify.
    apolar_max_hetero_ratio: float = 0.2
    # Size floor for evidence that is only physicochemical (a free fatty acyl or
    # an isoprenoid chain). No lipidomics workflow reports a ~100 Da methylated
    # alkene as a lipid; a recognised backbone (GP/GL/SP/ST) bypasses this.
    lipid_min_mw: float = 200.0
    # A steroid nucleus with no lipophilic tail and a logP at or below this is a
    # hormone, bile acid or oxysterol: water-soluble enough for the aqueous phase,
    # and reported by steroid-metabolomics panels as well as by lipidomics. Above
    # it sit the free sterols and steryl esters, which partition into the organic
    # phase only.
    polar_lipid_logp_max: float = 5.0
    # Shortest apolar tail that makes a 2-amino-1,3-diol a sphingoid base rather
    # than a small polar amino alcohol (serinol and friends).
    sphingoid_min_tail: int = 8
    # Shortest acyl chain on a glycerol backbone that makes it a glycerolipid
    # rather than a small glycerol ester.
    glycerolipid_min_tail: int = 8


# SMARTS for polar head groups, ordered specific -> generic. Matches are claimed
# greedily so a phosphate inside phosphocholine is not double-counted.
POLAR_HEAD_SMARTS: List[Tuple[str, str]] = [
    ("phosphocholine",    "OP(=O)([O-,OX2H,OX1])OCC[N+](C)(C)C"),
    ("phosphoethanolamine","OP(=O)([O-,OX2H,OX1])OCC[NX3;!$(NC=O)]"),
    ("phosphoserine",     "OP(=O)([O-,OX2H,OX1])OCC([NX3])C(=O)[O-,OX2H,OX1]"),
    ("phosphoinositol",   "OP(=O)([O-,OX2H,OX1])OC1C(O)C(O)C(O)C(O)C1O"),
    # acyclic carbons only: a sugar phosphate (glucose-6-phosphate, AMP, ATP) is not a
    # phosphoglycerol head group, it is a plain phosphate ester.
    ("phosphoglycerol",   "OP(=O)([O-,OX2H,OX1])O[C;!R][C;!R](O)[C;!R]O"),
    # two free oxygens + any third: matches free phosphate AND phosphate mono-esters
    # (sugar phosphates, nucleotides, phosphatidic acid), which are polar head groups too.
    ("phosphate",         "P(=O)([O-,OX2H,OX1])([O-,OX2H,OX1])O"),
    ("sulfate/sulfo",     "S(=O)(=O)[O-,OX2H,OX1]"),
    ("glycerol",          "[OX2H][CH2;!R][CH;!R]([OX2H])[CH2;!R][OX2H]"),
    ("quaternary_ammonium","[NX4+;!$([NX4+]~[!#6])]"),  # choline/carnitine style
]

# Head groups that, on a glycerol-3-phospho backbone, define a glycerophospholipid.
PHOSPHO_HEAD_NAMES = {
    "phosphocholine", "phosphoethanolamine", "phosphoserine",
    "phosphoinositol", "phosphoglycerol", "phosphate",
}

# Carbonyls that terminate a fatty-acyl tail. The carboxylate form is listed
# explicitly: HMDB/LIPID MAPS routinely store fatty acids deprotonated, and
# without it the same molecule would lose its acyl linkage depending only on
# whether the depositor drew "C(=O)O" or "C(=O)[O-]".
ACYL_CARBONYL_SMARTS = [
    "[CX3](=O)[OX2H0][#6]",   # ester
    "[CX3](=O)[OX2H1]",       # carboxylic acid
    "[CX3](=O)[OX1H0-]",      # carboxylate
    "[CX3](=O)[NX3]",         # amide
    "[CX3](=O)[SX2]",         # thioester
]

# Canonical lipid backbones ------------------------------------------------
# sn-glycerol-3-phosphate: three acyclic carbons each bearing an O, the third
# O going to phosphorus. Matches PC/PE/PS/PI/PG/PA, their lyso forms, the ether
# and plasmalogen forms, and cardiolipin.
GLYCEROPHOSPHO_SMARTS = "[CX4;!R]([OX2])[CX4;!R]([OX2])[CX4;!R][OX2][PX4](=O)"
# Plain glycerol backbone (no phosphorus needed): MG/DG/TG and ether glycerolipids.
GLYCEROL_BACKBONE_SMARTS = "[CX4;!R]([OX2])[CX4;!R]([OX2])[CX4;!R][OX2]"
# O-acyl ester onto an acyclic carbon: the "acyl" half of a glycerolipid.
O_ACYL_ESTER_SMARTS = "[CX3](=O)[OX2H0][CX4;!R]"
# Sphingoid base: 2-amino-1,3-diol. C1-O covers the free base, the
# phosphate (S1P/SM) and the glycosyl (HexCer) forms; the N covers both the
# free amine and the ceramide amide.
SPHINGOID_SMARTS = "[OX2,OX1][CX4H2;!R][CX4H;!R]([NX3])[CX4H;!R][OX2,OX1]"

# Isoprenoid detection (see detect_isoprenoid).
# An isoprene UNIT is counted as one non-aromatic C=C bond carrying at least one
# methyl. Counting methyl *groups* (the previous behaviour) double-counts a
# single isopropylidene terminus, which made every gem-dimethyl alkene -- tiglic
# acid, mesityl oxide, 2-methyl-3-buten-2-ol -- look like a two-unit isoprenoid.
ISOPRENE_ALKENE_METHYL_SMARTS = "[CH3][CX3]=[CX3]"   # a methyl on a C=C carbon
ISOPRENE_SATURATED_2UNIT_SMARTS = "[CH3][CHX4]([CX4])[CH2][CH2][CH2][CHX4][CH3]"

# Fused-ring nucleus geometry (see detect_isoprenoid).
STEROL_NUCLEUS_ATOMS = (16, 18)       # cyclopentanoperhydrophenanthrene = 17
TRITERPENE_NUCLEUS_ATOMS = (16, 23)   # hopane/oleanane/lupane skeletons
TRITERPENE_MIN_RING_METHYLS = 3


@dataclass
class Features:
    """All structural features extracted from one molecule."""
    smiles: str
    valid: bool
    mw: float = 0.0
    logp: float = 0.0
    has_carbon: bool = False
    n_carbons: int = 0                   # carbon atoms in the main organic fragment
    n_heteroatoms: int = 0               # non-carbon heavy atoms in that fragment
    has_long_acyl: bool = False          # lipid tail >= thr.acyl_min_carbons
    longest_chain_carbons: int = 0       # longest contiguous non-ring aliphatic C chain
    lipid_tail_carbons: int = 0          # longest mostly-hydrocarbon run within that chain
    is_acyl_linked: bool = False         # that apolar tail ends in an ester/amide/acid carbonyl
    has_isoprenoid: bool = False         # sterol/triterpene nucleus OR >= min_units isoprene units
    isoprenoid_type: Optional[str] = None  # 'sterol' | 'triterpenoid' | 'terpenoid/prenol' | None
    n_isoprene_units: int = 0            # methyl-bearing C=C bonds (see detect_isoprenoid)
    n_polar_head: int = 0
    polar_head_types: List[str] = field(default_factory=list)
    n_glycosyl: int = 0
    lipid_class: Optional[str] = None    # canonical LIPID MAPS-style backbone, if any
    n_fragments: int = 1                 # salts/mixtures: mw & logp use the largest organic part
    # atom-index lists for highlighting in visualization
    highlight: Dict[str, List[int]] = field(default_factory=dict)  # keys: acyl,isoprenoid,head,glycosyl


@dataclass
class Result:
    category: str            # '脂质组' | '代谢组' | '都有可能' | '都不可能' | '无效SMILES'
    reasons: List[str]       # human-readable rule hits (Chinese)
    features: Features
    rule: str = ""           # stable key of the branch that decided `category`


# Chinese names of the lipid classes, for the `reasons` text.
LIPID_CLASS_CN = {
    "glycerophospholipid": "甘油磷脂",
    "glycerolipid": "甘油酯",
    "sphingolipid": "鞘脂",
    "sterol lipid": "甾醇脂",
    "fatty acyl": "脂肪酰",
    "prenol lipid": "异戊烯醇脂",
    "hydrocarbon/fatty alcohol": "烃类/脂肪醇",
}


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


def _longest_paths(adj: Dict[int, List[int]]) -> List[List[int]]:
    """Longest simple path of every connected component of a forest `adj`.

    `adj` must be acyclic (a forest); each component's longest path is its
    diameter, found by double BFS. Returns one atom-index path per component.
    """
    def bfs(start: int) -> Tuple[List[int], Set[int]]:
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

    paths: List[List[int]] = []
    visited: Set[int] = set()
    for start in adj:
        if start in visited:
            continue
        path, component = bfs(start)
        visited.update(component)
        paths.append(path)
    return paths


def _forest_adjacency(mol: Chem.Mol, eligible: Set[int]) -> Dict[int, List[int]]:
    """Adjacency dict restricted to `eligible` atoms."""
    adj: Dict[int, List[int]] = {i: [] for i in eligible}
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if i in eligible and j in eligible:
            adj[i].append(j)
            adj[j].append(i)
    return adj


def aliphatic_chain_paths(mol: Chem.Mol) -> List[List[int]]:
    """Longest path of every connected run of non-aromatic, non-ring carbons.

    Returns one atom-index path per connected component (in atom order).
    ALGORITHM:
      - eligible atoms: GetAtomicNum()==6 and not GetIsAromatic() and not IsInRing().
      - build adjacency among eligible atoms (bond between two eligible atoms).
      - eligible subgraph is a forest (ring atoms excluded); for each connected
        component find the longest simple path (tree diameter) via double-BFS,
        returning the actual atom indices on that path.
      - empty -> [].
    """
    if mol is None:
        return []

    eligible = {
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic() and not atom.IsInRing()
    }
    if not eligible:
        return []
    return _longest_paths(_forest_adjacency(mol, eligible))


def longest_aliphatic_carbon_chain(mol: Chem.Mol) -> Tuple[int, List[int]]:
    """Longest contiguous chain of non-aromatic, non-ring carbon atoms.

    Returns (num_carbons_on_longest_path, atom_indices_on_that_path): the
    longest of aliphatic_chain_paths() (first one on ties); (0, []) if none.
    """
    best_path: List[int] = max(aliphatic_chain_paths(mol), key=len, default=[])
    return len(best_path), best_path


def _acyl_carbonyl_indices(mol: Chem.Mol) -> Set[int]:
    """Indices of every carbonyl carbon that anchors an acyl chain.

    Ester, carboxylic acid, carboxylate, amide and thioester. The carboxylate
    form matters: HMDB and LIPID MAPS routinely store fatty acids deprotonated,
    and without it the same molecule gained or lost its acyl linkage depending
    only on whether the depositor drew "C(=O)O" or "C(=O)[O-]".
    """
    out: Set[int] = set()
    if mol is None:
        return out
    for smarts in ACYL_CARBONYL_SMARTS:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue
        for match in mol.GetSubstructMatches(patt, uniquify=True):
            out.add(match[0])
    return out


def _polar_flags(mol: Chem.Mol, path: List[int]) -> List[bool]:
    """Per chain carbon: does it carry a heavy neighbour that is not carbon?"""
    return [
        any(nb.GetAtomicNum() not in (1, 6)
            for nb in mol.GetAtomWithIdx(idx).GetNeighbors())
        for idx in path
    ]


def _longest_lipid_window(flags: List[bool], max_fraction: float) -> Tuple[int, int]:
    """Longest slice of `flags` whose True-fraction is <= max_fraction.

    Returns (start, stop) as a half-open range; (0, 0) when nothing qualifies.
    """
    n = len(flags)
    prefix = [0] * (n + 1)
    for i, flag in enumerate(flags):
        prefix[i + 1] = prefix[i] + (1 if flag else 0)

    best = (0, 0)
    for i in range(n):
        for j in range(n, i, -1):
            if j - i <= best[1] - best[0]:
                break  # cannot beat the incumbent from here on
            if prefix[j] - prefix[i] <= max_fraction * (j - i):
                best = (i, j)
                break
    return best


def lipid_chain_paths(mol: Chem.Mol, max_polar_fraction: float = 0.3) -> List[List[int]]:
    """Per acyclic-carbon component, the longest run that is mostly hydrocarbon.

    A lipid tail is a greasy run, so counting bare carbons is not enough: on
    carbon count alone a fully hydroxylated C16 polyol, or any sugar alcohol,
    reads as a C16 "fatty chain". Trimming each aliphatic path down to its
    longest window with at most `max_polar_fraction` substituted carbons keeps
    oxygenated fatty acids (which are lipids) and drops polyols (which are not).
    """
    if mol is None:
        return []
    out: List[List[int]] = []
    for path in aliphatic_chain_paths(mol):
        i, j = _longest_lipid_window(_polar_flags(mol, path), max_polar_fraction)
        if j > i:
            out.append(path[i:j])
    return out


def _chain_is_acyl_linked(mol: Chem.Mol, path: List[int], carbonyls: Set[int]) -> bool:
    """Does this tail terminate in (or hang directly off) an acyl carbonyl?"""
    if not path or not carbonyls:
        return False
    for end in (path[0], path[-1]):
        if end in carbonyls:
            return True
        if any(nb.GetIdx() in carbonyls for nb in mol.GetAtomWithIdx(end).GetNeighbors()):
            return True
    return False


def detect_long_acyl(mol: Chem.Mol, min_carbons: int = 12,
                     max_polar_fraction: float = 0.3) -> Tuple[bool, int, bool, List[int]]:
    """Detect a long (>=min_carbons) fatty/acyl tail.

    Returns (has_long, tail_carbons, is_acyl_linked, atom_indices).
    ALGORITHM:
      - paths = lipid_chain_paths(mol, max_polar_fraction): each acyclic carbon
        run trimmed to its longest mostly-hydrocarbon window.
      - tail_carbons = length of the longest such run; has_long = tail_carbons >= min_carbons.
      - is_acyl_linked: ANY run that clears min_carbons ends on an ester /
        carboxylic acid / carboxylate / amide / thioester carbonyl. Checking
        every qualifying run (not only the longest) is what keeps a ceramide's
        N-acyl chain visible next to its longer sphingoid base.
      - atom_indices = every qualifying run, so a phospholipid highlights both
        of its acyl chains.
    """
    if mol is None:
        return (False, 0, False, [])

    paths = lipid_chain_paths(mol, max_polar_fraction)
    tail_carbons = max((len(p) for p in paths), default=0)
    has_long = tail_carbons >= min_carbons

    qualifying = [p for p in paths if len(p) >= min_carbons]
    carbonyls = _acyl_carbonyl_indices(mol)
    is_acyl = any(_chain_is_acyl_linked(mol, p, carbonyls) for p in qualifying)
    atoms = sorted({idx for p in qualifying for idx in p})
    return has_long, tail_carbons, is_acyl, atoms


def _fused_ring_components(mol: Chem.Mol) -> List[Tuple[List[int], Set[int], Counter]]:
    """Group SSSR rings into fused systems.

    Two rings are fused when they share >=2 atoms. Returns, per system,
    (ring indices, union of atom indices, Counter of ring sizes).
    """
    atom_rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    n = len(atom_rings)
    if not n:
        return []

    ring_adj: List[List[int]] = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if len(atom_rings[i] & atom_rings[j]) >= 2:
                ring_adj[i].append(j)
                ring_adj[j].append(i)

    out = []
    seen = [False] * n
    for i in range(n):
        if seen[i]:
            continue
        component = []
        stack = [i]
        seen[i] = True
        while stack:
            cur = stack.pop()
            component.append(cur)
            for nb in ring_adj[cur]:
                if not seen[nb]:
                    seen[nb] = True
                    stack.append(nb)

        atom_union: Set[int] = set()
        sizes: Counter = Counter()
        for ring_idx in component:
            atom_union |= atom_rings[ring_idx]
            sizes[len(atom_rings[ring_idx])] += 1
        out.append((component, atom_union, sizes))
    return out


def detect_isoprenoid(mol: Chem.Mol, min_units: int = 2) -> Tuple[bool, Optional[str], List[int]]:
    """Detect a sterol/triterpene nucleus or an isoprenoid (terpenoid/prenol) chain.

    Returns (has_isoprenoid, type_str, atom_indices).
    ALGORITHM:
      1) FUSED ALL-CARBON NUCLEUS. Group SSSR rings into fused systems (>=2
         shared atoms). A system qualifies only when every atom in it is carbon
         and at least three of its rings are non-aromatic -- the steroid nucleus
         is an all-carbon cyclopentanoperhydrophenanthrene, so requiring that
         rejects fused alkaloids (morphine's system carries the furan O and the
         piperidine N) while still accepting the aromatic A-ring of estrogens.
           - 'sterol': >=4 rings, >=3 six-membered and >=1 five-membered, 16-18 atoms.
           - 'triterpenoid': >=4 rings, 16-23 atoms and >=3 methyls on the
             nucleus -- hopane/lupane/oleanane skeletons are 21-atom systems that
             the sterol atom-count test alone throws away.
      2) ISOPRENOID CHAIN:
           - units = number of non-aromatic C=C bonds carrying >=1 methyl. One
             bond per unit, so an isopropylidene terminus counts once rather
             than twice.
           - a saturated 2-unit isoprene motif (phytanyl/dolichol/tocopherol
             tail) counts on its own.
           - isoprenoid if units >= min_units, or one saturated motif is present.
      3) else (False, None, []).
    """
    if mol is None:
        return (False, None, [])

    for component, atom_union, ring_sizes in _fused_ring_components(mol):
        if len(component) < 4:
            continue
        if not all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in atom_union):
            continue
        n_aliphatic_rings = sum(
            1 for r in mol.GetRingInfo().AtomRings()
            if set(r) <= atom_union and not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)
        )
        if n_aliphatic_rings < 3:
            continue

        size = len(atom_union)
        if (ring_sizes[6] >= 3 and ring_sizes[5] >= 1
                and STEROL_NUCLEUS_ATOMS[0] <= size <= STEROL_NUCLEUS_ATOMS[1]):
            return True, 'sterol', sorted(atom_union)

        if TRITERPENE_NUCLEUS_ATOMS[0] <= size <= TRITERPENE_NUCLEUS_ATOMS[1]:
            n_ring_methyl = sum(
                1 for atom in mol.GetAtoms()
                if atom.GetAtomicNum() == 6 and atom.GetTotalNumHs() == 3
                and any(nb.GetIdx() in atom_union for nb in atom.GetNeighbors())
            )
            if n_ring_methyl >= TRITERPENE_MIN_RING_METHYLS:
                return True, 'triterpenoid', sorted(atom_union)

    units, matched_atoms = count_isoprene_units(mol)

    sat_2unit = Chem.MolFromSmarts(ISOPRENE_SATURATED_2UNIT_SMARTS)
    n_sat = 0
    if sat_2unit is not None:
        matches = mol.GetSubstructMatches(sat_2unit, uniquify=True)
        n_sat = len(matches)
        for match in matches:
            matched_atoms.update(match)

    if units >= max(min_units, 1) or n_sat >= 1:
        return True, 'terpenoid/prenol', sorted(matched_atoms)

    return False, None, []


def count_isoprene_units(mol: Chem.Mol) -> Tuple[int, Set[int]]:
    """Count isoprene units as methyl-bearing non-aromatic C=C bonds.

    Returns (n_units, atoms involved). One unit per double bond: the previous
    implementation counted methyl groups, so a single isopropylidene terminus
    scored 2 and every gem-dimethyl alkene (tiglic acid, mesityl oxide) cleared
    a "2 isoprene units" threshold on its own.
    """
    atoms: Set[int] = set()
    if mol is None:
        return 0, atoms

    patt = Chem.MolFromSmarts(ISOPRENE_ALKENE_METHYL_SMARTS)
    if patt is None:
        return 0, atoms

    bonds: Set[int] = set()
    for methyl, c1, c2 in mol.GetSubstructMatches(patt, uniquify=True):
        bond = mol.GetBondBetweenAtoms(c1, c2)
        if bond is None:
            continue
        bonds.add(bond.GetIdx())
        atoms.update((methyl, c1, c2))
    return len(bonds), atoms


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


def _has_match(mol: Chem.Mol, smarts: str) -> bool:
    patt = Chem.MolFromSmarts(smarts)
    return bool(patt is not None and mol.GetSubstructMatch(patt))


def is_fatty_acid_skeleton(mol: Chem.Mol, thr: Thresholds) -> bool:
    """A carboxylic acid / ester / amide on a long, mostly-hydrocarbon skeleton.

    The chain-length test alone measures *acyclic* carbon runs, so it misses the
    prostanoids: PGE2 is a C20 fatty acid whose chain is cut in half by a
    cyclopentane ring, leaving two runs of eight. Measuring the skeleton instead
    of the longest run recovers that whole family (prostaglandins, thromboxanes,
    isoprostanes) without letting drug-like acids in, because an aromatic ring or
    a second ring disqualifies the molecule and the substituted-carbon density
    still has to stay below thr.chain_max_polar_fraction.
    """
    if mol is None or not _acyl_carbonyl_indices(mol):
        return False

    carbons = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 6]
    if not (thr.acyl_min_carbons <= len(carbons) <= thr.fatty_acid_max_carbons):
        return False
    if any(a.GetIsAromatic() for a in mol.GetAtoms()):
        return False

    rings = mol.GetRingInfo().AtomRings()
    if len(rings) > 1:
        return False
    if rings and (len(rings[0]) > 7
                  or any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in rings[0])):
        return False

    n_substituted = sum(
        1 for a in carbons
        if any(nb.GetAtomicNum() not in (1, 6) for nb in a.GetNeighbors())
    )
    return n_substituted <= thr.chain_max_polar_fraction * len(carbons)


def assign_lipid_class(mol: Chem.Mol, f: Features, thr: Thresholds) -> Optional[str]:
    """Name the canonical lipid backbone, or None.

    This is the decisive test in `classify`: a lipidomics run annotates by class,
    so a recognised backbone means "lipidome" no matter how small or polar the
    molecule is. That is exactly why the head-group features had to start
    carrying weight -- PC 10:0/10:0 and LPC 8:0 are sold as lipid standards, yet a
    chain-length-only rule called them "both" and "metabolome".

    Backbone classes (glycerophospholipid / glycerolipid / sphingolipid / sterol
    lipid) carry no size floor. The three that rest on a chain or a methylation
    count alone -- fatty acyl, prenol, hydrocarbon/fatty alcohol -- must clear
    thr.lipid_min_mw, which is what stops a ~100 Da methylated alkene from being
    called a lipid, and the latter two must also be heteroatom-poor.
    """
    if mol is None:
        return None

    if _has_match(mol, GLYCEROPHOSPHO_SMARTS) and any(
        h in PHOSPHO_HEAD_NAMES for h in f.polar_head_types
    ):
        return "glycerophospholipid"

    if (_has_match(mol, GLYCEROL_BACKBONE_SMARTS)
            and _has_match(mol, O_ACYL_ESTER_SMARTS)
            and f.lipid_tail_carbons >= thr.glycerolipid_min_tail):
        return "glycerolipid"

    if _has_match(mol, SPHINGOID_SMARTS) and f.lipid_tail_carbons >= thr.sphingoid_min_tail:
        return "sphingolipid"

    if f.isoprenoid_type in ("sterol", "triterpenoid"):
        return "sterol lipid"

    if f.mw >= thr.lipid_min_mw:
        if (f.has_long_acyl and f.is_acyl_linked) or is_fatty_acid_skeleton(mol, thr):
            return "fatty acyl"
        hydrocarbon_rich = (
            f.n_heteroatoms <= thr.apolar_max_hetero_ratio * max(f.n_carbons, 1)
        )
        if f.has_isoprenoid and hydrocarbon_rich:
            return "prenol lipid"
        # A bare greasy tail with no carbonyl is still a lipid: alkanes, waxes and
        # fatty alcohols (LIPID MAPS FA05/FA11) have no acyl linkage to find.
        if f.has_long_acyl and hydrocarbon_rich and f.logp >= thr.logp_polar_max:
            return "hydrocarbon/fatty alcohol"

    return None


def _main_organic_fragment(mol: Chem.Mol) -> Tuple[Chem.Mol, int]:
    """(largest carbon-containing fragment, number of fragments).

    Counterions and co-crystallised acids must not move the analyte's MW or
    logP: sodium palmitate and palmitic acid are the same lipid, and a TFA salt
    of glutamate is still glutamate.
    """
    try:
        frags = Chem.GetMolFrags(mol, asMols=True)
    except Exception:
        return mol, 1
    if len(frags) <= 1:
        return mol, max(1, len(frags))
    organic = [m for m in frags if any(a.GetAtomicNum() == 6 for a in m.GetAtoms())]
    if not organic:
        return mol, len(frags)
    return max(organic, key=lambda m: m.GetNumHeavyAtoms()), len(frags)


def extract_features(smiles: str, thr: Optional[Thresholds] = None) -> Features:
    """Run all detectors and assemble a Features object.

    ALGORITHM:
      - thr = thr or Thresholds()
      - mol = parse_smiles(smiles); if None -> Features(smiles, valid=False).
      - mw/logp are measured on the largest carbon-containing fragment so a
        counterion cannot shift them; every structural detector runs on the whole
        input, keeping the atom indices valid for highlighting.
      - (has_long, tail, acyl_linked, acyl_atoms) = detect_long_acyl(mol, thr.acyl_min_carbons)
      - (has_iso, iso_type, iso_atoms) = detect_isoprenoid(mol, thr.isoprene_min_units)
      - (n_head, head_types, head_atoms) = detect_polar_heads(mol)
      - (n_gly, gly_atoms) = detect_glycosyl(mol)
      - lipid_class = assign_lipid_class(mol, features, thr)
      - fill Features incl. highlight={'acyl':acyl_atoms,'isoprenoid':iso_atoms,
        'head':head_atoms,'glycosyl':gly_atoms}.
    """
    thr = thr or Thresholds()
    mol = parse_smiles(smiles)

    if mol is None:
        return Features(smiles=smiles, valid=False)

    main, n_fragments = _main_organic_fragment(mol)
    mw = Descriptors.MolWt(main)
    logp = Crippen.MolLogP(main)
    n_carbons = sum(1 for atom in main.GetAtoms() if atom.GetAtomicNum() == 6)
    n_heteroatoms = sum(
        1 for atom in main.GetAtoms() if atom.GetAtomicNum() not in (1, 6)
    )
    has_carbon = n_carbons > 0

    longest, _ = longest_aliphatic_carbon_chain(mol)
    has_long, tail, acyl_linked, acyl_atoms = detect_long_acyl(
        mol, thr.acyl_min_carbons, thr.chain_max_polar_fraction
    )
    has_iso, iso_type, iso_atoms = detect_isoprenoid(mol, thr.isoprene_min_units)
    n_units, _ = count_isoprene_units(mol)
    n_head, head_types, head_atoms = detect_polar_heads(mol)
    n_gly, gly_atoms = detect_glycosyl(mol)

    features = Features(
        smiles=smiles,
        valid=True,
        mw=mw,
        logp=logp,
        has_carbon=has_carbon,
        n_carbons=n_carbons,
        n_heteroatoms=n_heteroatoms,
        has_long_acyl=has_long,
        longest_chain_carbons=longest,
        lipid_tail_carbons=tail,
        is_acyl_linked=acyl_linked,
        has_isoprenoid=has_iso,
        isoprenoid_type=iso_type,
        n_isoprene_units=n_units,
        n_polar_head=n_head,
        polar_head_types=head_types,
        n_glycosyl=n_gly,
        n_fragments=n_fragments,
        highlight={
            "acyl": acyl_atoms,
            "isoprenoid": iso_atoms,
            "head": head_atoms,
            "glycosyl": gly_atoms,
        },
    )
    features.lipid_class = assign_lipid_class(mol, features, thr)
    return features


def classify(features: Features, thr: Optional[Thresholds] = None) -> Result:
    """Apply the tunable decision tree; return a Result with Chinese reasons.

    DECISION TREE (thr defaults in Thresholds). `Result.rule` names the branch
    that fired, so a display layer can render the verdict without re-deriving
    (and drifting from) the logic here.

      - not features.valid                      -> '无效SMILES'   rule='invalid'
      - not has_carbon                          -> '都不可能'     rule='no_carbon'
      - mw outside [mw_min, mw_max]             -> '都不可能'     rule='outside_window'
      - lipid_class is not None:
          * 'sterol lipid', no long tail, logp <= polar_lipid_logp_max
                                                -> '都有可能'     rule='steroid_both'
            (hormones, bile acids, oxysterols: steroid panels and lipidomics
             both report them; free sterols and steryl esters sit above the
             cut-off and stay lipidome-only)
          * boundary medium-chain free acyl     -> '都有可能'     rule='boundary_both'
            ('fatty acyl', tail and total carbons <= medium_chain_max, no polar
             head or sugar, logp in [boundary_logp_low, boundary_logp_high],
             mw <= boundary_mw_max)
          * otherwise                           -> '脂质组'       rule='lipid_class'
      - no lipid class:
          * logp < logp_polar_max               -> '代谢组'       rule='polar_metabolite'
          * otherwise                           -> '都有可能'     rule='nonpolar_both'
    """
    thr = thr or Thresholds()

    if not features.valid:
        return Result("无效SMILES", ["SMILES无法解析"], features, rule="invalid")

    reasons = [f"logP={features.logp:.2f}，分子量={features.mw:.1f} Da"]
    if features.n_fragments > 1:
        reasons.append(
            f"输入含{features.n_fragments}个片段（盐/混合物），"
            "logP与分子量按最大的含碳片段计算"
        )
    details = []

    if features.has_long_acyl:
        linkage_desc = ("与酯/酰胺/羧酸等羰基相连" if features.is_acyl_linked
                        else "未连接酯/酰胺/羧酸羰基")
        details.append(
            f"含≥{thr.acyl_min_carbons}碳疏水脂链"
            f"（最长疏水链C{features.lipid_tail_carbons}，最长碳链C{features.longest_chain_carbons}），"
            f"{linkage_desc}"
        )
    else:
        details.append(
            f"未检测到≥{thr.acyl_min_carbons}碳疏水脂链"
            f"（最长疏水链C{features.lipid_tail_carbons}，最长碳链C{features.longest_chain_carbons}）"
        )

    if features.has_isoprenoid:
        details.append(
            f"检测到类异戊二烯骨架（{features.isoprenoid_type}，"
            f"异戊二烯单元数={features.n_isoprene_units}）"
        )
    else:
        details.append(f"未检测到类异戊二烯骨架（异戊二烯单元数={features.n_isoprene_units}）")

    head_str = "、".join(features.polar_head_types) if features.polar_head_types else "无"
    details.append(f"极性头基数量={features.n_polar_head}（{head_str}）")

    if features.n_glycosyl:
        details.append(f"糖基环数量={features.n_glycosyl}（有）")
    else:
        details.append("糖基环数量=0（无）")

    if features.lipid_class:
        cn = LIPID_CLASS_CN.get(features.lipid_class, features.lipid_class)
        details.append(f"匹配脂质大类：{cn}（{features.lipid_class}）")
    else:
        details.append("未匹配到任何脂质大类骨架")

    in_window = thr.mw_min <= features.mw <= thr.mw_max
    reasons.append(
        f"分子量{features.mw:.1f} Da{'在' if in_window else '不在'}典型LC-MS窗口"
        f"（{thr.mw_min}-{thr.mw_max} Da）内"
    )

    if not features.has_carbon:
        reasons.extend(details)
        reasons.append("不含碳原子，非有机小分子；判为都不可能")
        return Result("都不可能", reasons, features, rule="no_carbon")

    if not in_window:
        reasons.extend(details)
        reasons.append("分子量超出典型LC-MS检测窗口；判为都不可能")
        return Result("都不可能", reasons, features, rule="outside_window")

    reasons.extend(details)

    amphipathic = features.n_polar_head + features.n_glycosyl

    if features.lipid_class:
        if (features.lipid_class == "sterol lipid"
                and not features.has_long_acyl
                and features.logp <= thr.polar_lipid_logp_max):
            reasons.append(
                f"甾核但无长疏水尾且logP={features.logp:.2f}≤{thr.polar_lipid_logp_max}，"
                "属激素/胆汁酸/氧化甾醇一类，脂质组与代谢组均常规检出"
            )
            return Result("都有可能", reasons, features, rule="steroid_both")

        # "Medium chain" has to mean the whole molecule is a C12-C14 fatty acid.
        # Testing only the longest acyclic run let PGE2 in: it is a C20 fatty
        # acid whose chain is split into two runs of eight by its cyclopentane.
        boundary = (
            features.lipid_class == "fatty acyl"
            and not features.has_isoprenoid
            and features.lipid_tail_carbons <= thr.medium_chain_max
            and features.n_carbons <= thr.medium_chain_max
            and thr.boundary_logp_low <= features.logp <= thr.boundary_logp_high
            and features.mw <= thr.boundary_mw_max
            and amphipathic == 0
        )
        if boundary:
            reasons.append("中等长度脂链、低极性且无极性头基，脂质组/代谢组均可能检出")
            return Result("都有可能", reasons, features, rule="boundary_both")

        cn = LIPID_CLASS_CN.get(features.lipid_class, features.lipid_class)
        reasons.append(f"存在明确的{cn}骨架，判为脂质组")
        return Result("脂质组", reasons, features, rule="lipid_class")

    if features.logp < thr.logp_polar_max:
        reasons.append(
            f"无脂质大类骨架且logP={features.logp:.2f}<{thr.logp_polar_max}，"
            "倾向于极性代谢物，归为代谢组"
        )
        return Result("代谢组", reasons, features, rule="polar_metabolite")

    reasons.append(
        f"无脂质大类骨架但logP={features.logp:.2f}≥{thr.logp_polar_max}，"
        "两类检测均可能"
    )
    return Result("都有可能", reasons, features, rule="nonpolar_both")


def analyze(smiles: str, thr: Optional[Thresholds] = None) -> Result:
    """Top-level convenience: extract_features then classify."""
    return classify(extract_features(smiles, thr), thr)
