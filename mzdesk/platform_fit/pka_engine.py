"""pKa 预测与 pH 行为引擎（MolGpKa），只依赖 numpy 和 RDKit。

模型：Pan, Wang, Yang, Zhang, "MolGpKa: A Web Server for Small Molecule pKa Prediction Using a Graph-Convolutional
Neural Network", J. Chem. Inf. Model. 2021, 61, 3159-3165（https://doi.org/10.1021/acs.jcim.1c00075）。
原代码 github.com/Xundrug/MolGpKa（MIT，PyTorch + torch_geometric），权重是作者放在仓库 models/ 里的 weight_acid.pth /
weight_base.pth，由 pka_setup.py 一次性转成 numpy。这里按原代码的网络用 numpy 重写前向计算（不需要 torch、torch_geometric），
特征和位点匹配与原代码一致，并和原代码逐位点对照验证过（validation/molgpka/）。

做法：
1. 位点：143 条 SMARTS（molgpka_smarts.py）在加了显式 H 的分子上匹配出酸性位点（被失去的 H）和碱性位点（接收质子的重原子）。
2. 每个位点一次前向：29 维原子特征 → 5 层 GCN（对称归一化邻接，ReLU，BatchNorm）→ 注意力池化 → 3 层全连接 → pKa。
   酸、碱各一套权重。特征里有一维是各原子到该位点的距离，另一维标出该位点，所以每个位点要单独算一遍。
3. 预测之后去掉化学上不呈碱性的 N 的碱位点（酰胺 / 磺酰胺 / 吡咯型芳香 N，规则见 NONBASIC_N；酸位点保留）。
4. 各位点当作互相独立，用 Henderson–Hasselbalch 得到每个 pH 下的净电荷分布、完全不带电的分子比例，再得到 logD(pH)。

和原代码的差别（有意的）：原代码先把分子转成 SMILES 再读回来（原子编号会变），这里保留输入的原子顺序，位点编号才能对上输入的结构；
原代码的氢键供体 / 受体两个特征因为 `int in list-of-tuples` 永远为 False，训练时也是如此，所以这里这两列恒为 0。
"""
import os
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

from mzdesk.platform_fit import pka_setup
from mzdesk.platform_fit.molgpka_smarts import PATTERNS

RDLogger.DisableLog("rdApp.*")

N_FEATURES = 29
BN_EPS = 1e-5                                    # torch BatchNorm1d 默认
SYMBOLS = ("C", "H", "O", "N", "S", "Cl", "F", "Br", "P", "I")      # 原代码的 one-hot；其余元素落到最后一项（I）
HYBRIDIZATION = (Chem.HybridizationType.SP, Chem.HybridizationType.SP2, Chem.HybridizationType.SP3,
                 Chem.HybridizationType.SP3D, Chem.HybridizationType.SP3D2)
RING_SIZES = (3, 4, 5, 6, 7, 8)
WARN_HEAVY_ATOMS = 100                           # 超过这个大小提醒（MolGpKa 在类药小分子上训练）
MAX_HEAVY_ATOMS = 300                            # 超过就拒绝（稠密邻接矩阵的开销）
PH_GRID = np.round(np.arange(0, 141) * 0.1, 1)   # pH 0-14，步长 0.1


class PkaError(RuntimeError):
    """权重缺失或损坏（不是用户输入的问题）。"""


# ---------------------------------------------------------------- weights
def weights_path():
    return os.environ.get("MZDESK_MOLGPKA_WEIGHTS") or str(pka_setup.default_weights())


@lru_cache(maxsize=2)
def _load(path):
    """npz -> {"acid": model, "base": model}；BatchNorm（推理时）折成 x * scale + shift。"""
    if not os.path.isfile(path):
        raise PkaError(f"MolGpKa weights not found at {path}. Run once: {pka_setup.SETUP_COMMAND}")
    try:
        with np.load(path) as z:
            w = {k: z[k] for k in z.files}
        models = {}
        for kind in pka_setup.KINDS:
            m = {"conv": [], "bn": []}
            for i in range(1, 6):
                m["conv"].append((w[f"{kind}.conv{i}.weight"], w[f"{kind}.conv{i}.bias"]))
                scale = w[f"{kind}.bn{i}.weight"] / np.sqrt(w[f"{kind}.bn{i}.running_var"] + np.float32(BN_EPS))
                m["bn"].append((scale.astype(np.float32), (w[f"{kind}.bn{i}.bias"] - w[f"{kind}.bn{i}.running_mean"] * scale).astype(np.float32)))
            m["att"] = (w[f"{kind}.att.gate_nn.weight"][0], w[f"{kind}.att.gate_nn.bias"][0])
            m["fc"] = [(w[f"{kind}.fc{i}.weight"], w[f"{kind}.fc{i}.bias"]) for i in (2, 3, 4)]
            if m["conv"][0][0].shape != (N_FEATURES, 1024) or m["fc"][2][0].shape != (1, 16):
                raise ValueError("unexpected tensor shapes")
            models[kind] = m
    except (KeyError, ValueError, OSError) as exc:
        raise PkaError(f"MolGpKa weights at {path} are unreadable ({exc!r}). Rebuild with: {pka_setup.SETUP_COMMAND} --force") from None
    return models


# ---------------------------------------------------------------- sites
@lru_cache(maxsize=1)
def _patterns():
    out = []
    for _id, smarts, index, kind in PATTERNS:
        patt = Chem.MolFromSmarts(smarts)
        out.append((patt, [int(i) for i in index.split(",")], kind))
    return out


def find_site_atoms(mol_h):
    """Acid sites (index of the H atom that is lost) and base sites (index of the heavy atom), both sorted, on a molecule with explicit H."""
    acid, base = set(), set()
    for patt, indices, kind in _patterns():
        for match in mol_h.GetSubstructMatches(patt):
            (acid if kind == "A" else base).update(match[i] for i in indices)
    return sorted(acid), sorted(base)


# ---------------------------------------------------------------- network
def static_features(mol_h):
    """Per-atom features that do not depend on the site (27 columns): N x 27 float32.

    Same order as the original descriptor.py: element one-hot (10), degree, hybridization one-hot (5), implicit valence,
    aromatic, in a ring of size 3-8 (6), H-bond donor, H-bond acceptor, formal charge.  The donor and acceptor columns are
    always 0 (see the module docstring), exactly as in the original code and its training data.
    """
    ring = mol_h.GetRingInfo()
    rows = []
    for atom in mol_h.GetAtoms():
        idx = atom.GetIdx()
        symbol = atom.GetSymbol() if atom.GetSymbol() in SYMBOLS else SYMBOLS[-1]
        hyb = atom.GetHybridization() if atom.GetHybridization() in HYBRIDIZATION else HYBRIDIZATION[-1]
        row = [s == symbol for s in SYMBOLS]
        row.append(atom.GetDegree())
        row += [hyb == h for h in HYBRIDIZATION]
        row.append(atom.GetImplicitValence())
        row.append(atom.GetIsAromatic())
        row += [ring.IsAtomInRingOfSize(idx, n) for n in RING_SIZES]
        row += [False, False]
        row.append(atom.GetFormalCharge())
        rows.append(row)
    return np.asarray(rows, dtype=np.float32)


def normalized_adjacency(mol_h):
    """D^-1/2 (A + I) D^-1/2 as a dense float32 matrix (the GCNConv propagation of the original code)."""
    n = mol_h.GetNumAtoms()
    a = Chem.GetAdjacencyMatrix(mol_h).astype(np.float32) + np.eye(n, dtype=np.float32)
    d = a.sum(axis=1) ** -0.5
    return (a * d[:, None] * d[None, :]).astype(np.float32)


def _forward(x, a_hat, m):
    h = x
    for (w, b), (scale, shift) in zip(m["conv"], m["bn"]):
        h = np.maximum(a_hat @ (h @ w) + b, 0, dtype=np.float32)        # GCNConv -> ReLU -> BatchNorm (the original order)
        h = h * scale + shift
    gate = h @ m["att"][0] + m["att"][1]                                # global attention pooling: softmax over the atoms
    alpha = np.exp(gate - gate.max())
    h = (alpha / alpha.sum()) @ h
    for (w, b), last in zip(m["fc"], (False, False, True)):
        h = w @ h + b
        if not last:
            h = np.maximum(h, 0)
    return float(h[0])


def predict_aids(mol_h, kind, aids, distances=None, static=None, a_hat=None):
    """{aid: pKa} for the given site atoms of ``mol_h`` (explicit H) with the "acid" or "base" weights."""
    model = _load(weights_path())[kind]
    static = static_features(mol_h) if static is None else static
    a_hat = normalized_adjacency(mol_h) if a_hat is None else a_hat
    dist = Chem.GetDistanceMatrix(mol_h) if distances is None else distances
    out = {}
    for aid in aids:
        x = np.zeros((static.shape[0], N_FEATURES), dtype=np.float32)
        x[:, :27] = static
        x[:, 27] = np.where(dist[aid] < 1e7, dist[aid] + 1, 0)           # original: len(shortest path) = bonds + 1; 0 for the site itself
        x[aid, 27] = 0                                                   # and for atoms in another fragment (no path); never happens here
        x[aid, 28] = 1
        out[aid] = _forward(x, a_hat, model)
    return out


# ---------------------------------------------------------------- molecule preparation
@dataclass
class Prepared:
    mol: Chem.Mol                  # neutral form of the largest fragment, heavy atoms only (atom order of the input when it has one fragment)
    notes: list                    # what was changed, for the page warnings
    fixed_charge: int              # net formal charge that cannot be neutralized (quaternary N+ ...); nitro / N-oxide pairs cancel
    permanent: bool                # has a charged atom with no opposite charge next to it (quaternary N+, N-alkylpyridinium ...)


def _neutralize_anions(mol):
    """Protonate O- / S- / N- that are not next to a cation: carboxylate, phosphate, phenolate, sulfonate ...

    RDKit's Uncharger leaves such an anion alone when a quaternary N+ is in the molecule (betaine, phosphatidylcholine), which
    would hide an ionizable site. Charge-separated groups (nitro, N-oxide, azide: anion next to a cation) are left as they are.
    """
    rw = Chem.RWMol(mol)
    changed = False
    for atom in rw.GetAtoms():
        if atom.GetFormalCharge() == -1 and atom.GetSymbol() in ("O", "S", "N") and not any(n.GetFormalCharge() > 0 for n in atom.GetNeighbors()):
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
            changed = True
    if not changed:
        return mol, False
    try:
        Chem.SanitizeMol(rw)
    except Exception:  # noqa: BLE001 — keep the molecule as it was
        return mol, False
    return rw.GetMol(), True


def _largest_fragment(mol):
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) == 1:
        return mol, 1
    def key(f):
        return (any(a.GetSymbol() == "C" for a in f.GetAtoms()), f.GetNumHeavyAtoms(), Descriptors.MolWt(f))
    return max(frags, key=key), len(frags)


def prepare(smiles):
    """SMILES -> Prepared. Raises ValueError for a structure that cannot be used."""
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        raise ValueError("RDKit could not parse the structure.")
    notes = []
    mol, n_frag = _largest_fragment(mol)
    if n_frag > 1:
        notes.append(f"{n_frag} disconnected parts (salt or mixture): only the largest is analysed.")
    before = Chem.MolToSmiles(mol)
    mol, _ = _neutralize_anions(mol)
    mol = rdMolStandardize.Uncharger().uncharge(mol)
    if Chem.MolToSmiles(mol) != before:
        notes.append("Input neutralized (charges that can be removed by adding or losing H); pKa is predicted for the neutral form.")
    if mol.GetNumHeavyAtoms() > MAX_HEAVY_ATOMS:
        raise ValueError(f"Structure has {mol.GetNumHeavyAtoms()} heavy atoms; MolGpKa is limited to {MAX_HEAVY_ATOMS}.")
    if mol.GetNumHeavyAtoms() > WARN_HEAVY_ATOMS:
        notes.append(f"{mol.GetNumHeavyAtoms()} heavy atoms: far larger than the drug-like molecules MolGpKa was trained on; treat pKa values with caution.")
    odd = sorted({a.GetSymbol() for a in mol.GetAtoms()} - set(SYMBOLS))
    if odd:
        notes.append(f"Element(s) {', '.join(odd)} not in the MolGpKa training set (C H N O S P F Cl Br I): pKa values are unreliable.")
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    fixed = sum(a.GetFormalCharge() for a in charged)
    permanent = any(not any(n.GetFormalCharge() * a.GetFormalCharge() < 0 for n in a.GetNeighbors()) for a in charged)
    if permanent:
        notes.append(f"Permanent charge (net {fixed:+d}, e.g. quaternary N+): counted as fixed; there is no uncharged form, so logD is not defined.")
    return Prepared(mol=mol, notes=notes, fixed_charge=fixed, permanent=permanent)


# ---------------------------------------------------------------- sites and pKa
@dataclass
class Site:
    atom: int          # 0-based index of the heavy atom (the one that holds / takes the proton)
    element: str
    kind: str          # "acid" | "base"
    pka: float

    @property
    def label(self):
        return f"{self.element}{self.atom + 1}"          # 1-based, like the atom numbers of the Derivatization page


def predict_sites(mol):
    """Ionizable sites of the neutral molecule ``mol`` (no explicit H), sorted by pKa."""
    mol_h = Chem.AddHs(mol)                                # heavy atoms keep their indices; the H atoms come last
    acid_h, base_atoms = find_site_atoms(mol_h)
    static, a_hat, dist = static_features(mol_h), normalized_adjacency(mol_h), Chem.GetDistanceMatrix(mol_h)
    heavy_of = {h: mol_h.GetAtomWithIdx(h).GetNeighbors()[0].GetIdx() for h in acid_h}
    first_h = {}                                           # NH2 and the like: the H atoms are equivalent, predict one
    for h in acid_h:
        first_h.setdefault(heavy_of[h], h)
    acid = predict_aids(mol_h, "acid", sorted(first_h.values()), dist, static, a_hat)
    base = predict_aids(mol_h, "base", base_atoms, dist, static, a_hat)
    sites = [Site(heavy_of[h], mol_h.GetAtomWithIdx(heavy_of[h]).GetSymbol(), "acid", p) for h, p in acid.items()]
    sites += [Site(a, mol_h.GetAtomWithIdx(a).GetSymbol(), "base", p) for a, p in base.items()]
    return sorted(sites, key=lambda s: (s.pka, s.atom))


# ---------------------------------------------------------------- chemically non-basic N
# Textbook rules for N whose lone pair is not available for protonation. The MolGpKa patterns still match such N as "base"
# sites (e.g. the anilide N of lidocaine gets pKa 5.2), so these base sites are dropped after prediction. Acid sites on the
# same atoms stay (sulfonamide NH, imide NH ...). First matching rule gives the reason.
NONBASIC_N = (
    ("amide N", "[NX3][CX3]=[O,S]"),                 # amide / imide / carbamate / urea / thioamide: N on a C=O or C=S carbon
    ("sulfonamide N", "[NX3]S(=O)=O"),               # sulfonamide / sulfonyl N
    ("pyrrole-type aromatic N", "[nX3]"),            # aromatic N with three connections: lone pair is part of the aromatic sextet
    ("pyrrole-type aromatic N", "[nH]"),             # (indole NH, N-alkyl imidazole / purine ring N); pyridine-type [nX2] stays a base
)


@lru_cache(maxsize=1)
def _nonbasic_patterns():
    return [(reason, Chem.MolFromSmarts(smarts)) for reason, smarts in NONBASIC_N]


def drop_nonbasic_sites(mol, sites):
    """Remove the base sites on chemically non-basic N (``NONBASIC_N``) -> (kept sites, [(dropped site, reason), ...]).

    ``mol`` is the neutral molecule the sites were predicted on (same atom indices). Acid sites and base sites on other atoms are kept.
    """
    reason_of = {}
    for reason, patt in _nonbasic_patterns():
        for match in mol.GetSubstructMatches(patt):
            reason_of.setdefault(match[0], reason)
    kept, dropped = [], []
    for s in sites:
        if s.kind == "base" and s.element == "N" and s.atom in reason_of:
            dropped.append((s, reason_of[s.atom]))
        else:
            kept.append(s)
    return kept, dropped


# ---------------------------------------------------------------- Henderson-Hasselbalch, independent sites
def charged_fraction(site, ph):
    """Fraction of molecules with this site ionized: acid -> A- (charge -1), base -> BH+ (charge +1)."""
    ph = np.asarray(ph, dtype=float)
    return 1.0 / (1.0 + 10.0 ** ((site.pka - ph) if site.kind == "acid" else (ph - site.pka)))


def charge_states(sites, ph, fixed_charge=0):
    """Net-charge distribution at each pH: (charges, probabilities [n_pH x n_charges], fraction fully uncharged [n_pH]).

    Sites are independent. The net charge is the fixed charge plus -1 per ionized acid and +1 per ionized base, so a
    zwitterion (COO- and NH3+) counts as net 0; the fully uncharged fraction is the product of the neutral fractions.
    """
    ph = np.atleast_1d(np.asarray(ph, dtype=float))
    n_acid = sum(s.kind == "acid" for s in sites)
    n_base = len(sites) - n_acid
    qmin = fixed_charge - n_acid
    dist = np.zeros((len(ph), n_acid + n_base + 1))
    dist[:, fixed_charge - qmin] = 1.0
    uncharged = np.ones(len(ph))
    for s in sites:
        p = charged_fraction(s, ph)[:, None]
        uncharged = uncharged * (1.0 - p[:, 0])
        stay, moved = dist * (1.0 - p), dist * p
        dist = stay
        if s.kind == "acid":
            dist[:, :-1] += moved[:, 1:]
        else:
            dist[:, 1:] += moved[:, :-1]
    return np.arange(qmin, qmin + dist.shape[1]), dist, uncharged


def crippen_logp(mol):
    return float(Crippen.MolLogP(mol))


@dataclass
class Analysis:
    mol: Chem.Mol
    notes: list
    sites: list                    # ionizable sites after the non-basic-N filter, sorted by pKa
    ignored: list                  # [(Site, reason)]: base sites on non-basic N that were predicted but dropped
    fixed_charge: int
    permanent: bool
    logp: float
    tpsa: float
    mw: float
    n_stereo: int                  # stereocentres, counted even when the input does not specify them
    n_unassigned: int              # how many of them the input leaves unspecified
    meso: bool                     # all stereocentres are specified and the molecule equals its mirror image (achiral)
    charges: np.ndarray            # possible net charges (columns of ``dist``)
    ph: np.ndarray                 # PH_GRID
    dist: np.ndarray               # P(net charge) at each grid pH
    mean_charge: np.ndarray        # expected net charge on the grid
    f_neutral: np.ndarray          # fraction fully uncharged on the grid (0 when permanent)
    logd: np.ndarray               # logP + log10(f_neutral) on the grid; NaN when permanent or f_neutral underflows

    def at(self, ph):
        """State at one pH: {"charge", "p_charge", "f_neutral", "logd", "mean_charge", "dist"} (logd None when not defined)."""
        charges, dist, uncharged = charge_states(self.sites, [ph], self.fixed_charge)
        f = 0.0 if self.permanent else float(uncharged[0])
        k = int(np.argmax(dist[0]))
        return {"charge": int(charges[k]), "p_charge": float(dist[0, k]), "f_neutral": f,
                "logd": (self.logp + float(np.log10(f))) if f > 0 else None,
                "mean_charge": float(dist[0] @ charges), "charges": charges, "dist": dist[0]}


def _stereo(mol):
    """(number of stereocentres, how many unassigned, meso flag): meso = every centre assigned and the mirror image equals the molecule."""
    centres = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    if not centres:
        return 0, 0, False
    unassigned = sum(c == "?" for _, c in centres)
    if unassigned:
        return len(centres), unassigned, False
    mirror = Chem.Mol(mol)
    for a in mirror.GetAtoms():
        if a.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        elif a.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    return len(centres), 0, Chem.MolToSmiles(mirror) == Chem.MolToSmiles(mol)


def analyze(smiles):
    """SMILES -> Analysis (pKa sites, curves over pH 0-14, logP / TPSA / MW, stereocentres). Raises ValueError / PkaError."""
    prep = prepare(smiles)
    sites, ignored = drop_nonbasic_sites(prep.mol, predict_sites(prep.mol))
    charges, dist, uncharged = charge_states(sites, PH_GRID, prep.fixed_charge)
    logp = crippen_logp(prep.mol)
    f = np.zeros_like(uncharged) if prep.permanent else uncharged
    with np.errstate(divide="ignore"):
        logd = np.where(f > 0, logp + np.log10(np.where(f > 0, f, 1.0)), np.nan)
    n_stereo, n_unassigned, meso = _stereo(prep.mol)
    return Analysis(mol=prep.mol, notes=prep.notes, sites=sites, ignored=ignored, fixed_charge=prep.fixed_charge, permanent=prep.permanent,
                    logp=logp, tpsa=float(rdMolDescriptors.CalcTPSA(prep.mol)), mw=float(Descriptors.MolWt(prep.mol)),
                    n_stereo=n_stereo, n_unassigned=n_unassigned, meso=meso, charges=charges, ph=PH_GRID, dist=dist,
                    mean_charge=dist @ charges, f_neutral=f, logd=logd)
