"""Compound index behind the Known methods page: normalisation, fingerprints, descriptors, charge classes, vectorised search.

One ``CompoundIndex`` holds the structures of one source (JP18 or RepoRT) as numpy arrays, so a query is a few array operations
over ~100k compounds (Tanimoto on packed Morgan bits by popcount, property scores by broadcasting), not a Python loop.

Normalisation (query and library alike): largest fragment, neutralise (RDKit MolStandardize), InChIKey. The "connectivity key" is the
first 14 characters of that InChIKey (formula + skeleton + H layers; no stereo, no isotopes).

Every weight, scale and cut-off of the similarity score is in ``Config``.
"""
import json
import math
import os
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

SCHEMA = 2                                   # bump when the arrays or their meaning change (forces a rebuild of the JP cache)
METHOD_MOLGPKA = "MolGpKa"
METHOD_FG = "functional groups"
CHARGE_NAMES = {-1: "negative", 0: "neutral", 1: "positive"}


# ---------------------------------------------------------------- the one place for weights, scales and cut-offs
@dataclass(frozen=True)
class Config:
    # similar-compound score = w_structure * Tanimoto + w_properties * mean(property terms)
    w_structure: float = 0.6
    w_properties: float = 0.4
    # property terms, in this order: exp(-|dlogP| / scale_logp), exp(-|dTPSA| / scale_tpsa), exp(-|dMW| / scale_mw),
    # same charge class at ph_acidic, same charge class at ph_neutral (1 or 0); prop_weights weigh them in the mean
    scale_logp: float = 1.0
    scale_tpsa: float = 30.0
    scale_mw: float = 100.0
    ph_acidic: float = 2.7                   # typical acidic LC-MS mobile phase
    ph_neutral: float = 7.4
    prop_weights: tuple = (1.0, 1.0, 1.0, 1.0, 1.0)
    # Morgan fingerprint (fixed when an index is built; an index made with other values is refused)
    fp_radius: int = 2
    fp_bits: int = 2048
    # what the page shows
    min_score: float = 0.5                   # similar compounds below this score are not listed
    min_structure: float = 0.25              # floor on the Tanimoto term alone: weak structural similarity cannot be made up by matching properties (0 = none)
    max_similar: int = 10                    # similar compounds listed
    max_methods_per_similar: int = 3         # methods listed per similar compound
    max_similar_per_dataset: int = 2         # similar hits taken from one RepoRT dataset (SMRT alone holds 80k compounds)
    max_exact_per_source: int = 25           # methods listed per source for the exact compound


CONFIG = Config()


# ---------------------------------------------------------------- directories
def methods_dir():
    """Local data folder of the Known methods page (RepoRT index, JP18 index cache). MZDESK_METHODS_DIR overrides it."""
    return Path(os.environ.get("MZDESK_METHODS_DIR") or Path.home() / ".mzdesk" / "data" / "methods")


# ---------------------------------------------------------------- normalisation
_tools = {}


def _chooser():
    if "largest" not in _tools:
        _tools["largest"] = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
        _tools["uncharger"] = rdMolStandardize.Uncharger()
    return _tools["largest"], _tools["uncharger"]


@dataclass
class Normalized:
    mol: Chem.Mol
    smiles: str
    inchikey: str
    notes: list

    @property
    def key14(self):
        return self.inchikey[:14]

    @property
    def stereo_block(self):
        return self.inchikey[15:25]


def normalize(smiles):
    """SMILES -> Normalized (largest fragment, neutralised, InChIKey). Raises ValueError for a structure that cannot be used."""
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None or mol.GetNumHeavyAtoms() == 0:
        raise ValueError("RDKit could not parse the structure." if mol is None else "The structure has no heavy atoms.")
    largest, uncharger = _chooser()
    notes = []
    if len(Chem.GetMolFrags(mol)) > 1:
        mol = largest.choose(mol)
        notes.append("disconnected parts (salt or mixture): only the largest fragment is used")
    if any(a.GetFormalCharge() for a in mol.GetAtoms()):
        neutral = uncharger.uncharge(mol)
        if Chem.MolToSmiles(neutral) != Chem.MolToSmiles(mol):
            mol = neutral
            notes.append("charges removed (neutralised)")
    ik = Chem.MolToInchiKey(mol)
    if not ik:
        raise ValueError("RDKit could not make an InChIKey for this structure.")
    return Normalized(mol=mol, smiles=Chem.MolToSmiles(mol), inchikey=ik, notes=notes)


# ---------------------------------------------------------------- fingerprints
_fpgens = {}


def fingerprint(mol, cfg=CONFIG):
    """Morgan bit fingerprint packed to bytes (fp_bits / 8)."""
    key = (cfg.fp_radius, cfg.fp_bits)
    if key not in _fpgens:
        _fpgens[key] = rdFingerprintGenerator.GetMorganGenerator(radius=cfg.fp_radius, fpSize=cfg.fp_bits)
    return np.packbits(_fpgens[key].GetFingerprintAsNumPy(mol)).tobytes()


_M1, _M2, _M4, _H01 = (np.uint64(v) for v in (0x5555555555555555, 0x3333333333333333, 0x0F0F0F0F0F0F0F0F, 0x0101010101010101))


def popcount(x):
    """Set bits per uint64 element (np.bitwise_count on numpy >= 2, bit tricks on older numpy)."""
    if hasattr(np, "bitwise_count"):
        return np.bitwise_count(x).astype(np.uint8)
    x = x - ((x >> np.uint64(1)) & _M1)
    x = (x & _M2) + ((x >> np.uint64(2)) & _M2)
    x = (x + (x >> np.uint64(4))) & _M4
    return ((x * _H01) >> np.uint64(56)).astype(np.uint8)


def tanimoto(fps, pops, q_fp, q_pop):
    """Tanimoto of every row of ``fps`` (N x words, uint64) with the query ``q_fp`` (words,)."""
    inter = popcount(fps & q_fp).sum(axis=1, dtype=np.int32)
    union = pops.astype(np.int32) + int(q_pop) - inter
    out = np.zeros(len(fps), dtype=np.float32)
    np.divide(inter, union, out=out, where=union > 0)
    return out


# ---------------------------------------------------------------- charge classes
# Class = sign of the dominant net charge at a pH: -1 / 0 / +1.
_FG = {name: Chem.MolFromSmarts(s) for name, s in {
    "carboxylic": "[CX3](=O)[OX2H1,OX1-]",
    "sulfonic": "[SX4](=O)(=O)[OX2H1,OX1-]",
    "phosphoric": "[PX4](=O)[OX2H1,OX1-]",
    "amine": "[NX3;!$(N-[!#6;!#1]);!$(N-[#6]=[#6,#7,#8,#16]);!$(N-[#6]#*);!$(N-a)]",
    "amidine": "[NX3;!$(N-C=O)]-[CX3;!$(C-[OX2,SX2])]=[NX2;!$(N-[!#6;!#1]);!$(N-C=O)]",
    "quaternary": "[NX4+]",
    "pyridinium": "[n+;!$([n+][O-])]",
}.items()}


def functional_group_classes(mol):
    """Fallback charge classes from functional groups (no pKa prediction): (class at pH 2.7, class at pH 7.4).

    Sulfonic and phosphoric acids are ionised at both pH, carboxylic acids only at the neutral pH (pKa about 4); a phosphoric acid with two
    acidic OH loses both protons at pH 7.4. Aliphatic amines, amidines / guanidines and quaternary N are positive at both.
    Net charge = positives - negatives, so an amino acid is neutral at pH 7.4.
    """
    def count(name):
        return len(mol.GetSubstructMatches(_FG[name]))
    amidines = len({m[1] for m in mol.GetSubstructMatches(_FG["amidine"])})
    pos = count("amine") + amidines + count("quaternary") + count("pyridinium")
    sulfonic, carboxylic = count("sulfonic"), count("carboxylic")
    per_p = Counter(m[0] for m in mol.GetSubstructMatches(_FG["phosphoric"]))
    low = pos - (sulfonic + len(per_p))
    neutral = pos - (sulfonic + carboxylic + sum(min(n, 2) for n in per_p.values()))
    return int(np.sign(low)), int(np.sign(neutral))


_molgpka = {}


def molgpka_available():
    """(True, "") when the MolGpKa weights can be used, else (False, reason)."""
    if "state" not in _molgpka:
        from mzdesk.platform_fit import pka_engine as pe
        try:
            pe._load(pe.weights_path())
            _molgpka["state"] = (True, "")
        except pe.PkaError as exc:
            _molgpka["state"] = (False, str(exc))
    return _molgpka["state"]


def reset_molgpka_check():
    _molgpka.clear()


def molgpka_classes(smiles, phs):
    """MolGpKa charge classes of the (already normalised) SMILES at the pH values ``phs``. Raises ValueError for an unusable structure."""
    from mzdesk.platform_fit import pka_engine as pe
    prep = pe.prepare(smiles)
    sites, _dropped = pe.drop_nonbasic_sites(prep.mol, pe.predict_sites(prep.mol))
    charges, dist, _ = pe.charge_states(sites, list(phs), prep.fixed_charge)
    return tuple(int(np.sign(charges[k])) for k in dist.argmax(axis=1))


def charge_classes(norm, method, cfg=CONFIG):
    """-> (class at ph_acidic, class at ph_neutral, method actually used). ``method`` is METHOD_MOLGPKA or METHOD_FG; a structure
    MolGpKa cannot handle (too large ...) falls back to the functional groups. The fallback is defined for pH 2.7 and 7.4."""
    if method == METHOD_MOLGPKA:
        try:
            q = molgpka_classes(norm.smiles, (cfg.ph_acidic, cfg.ph_neutral))
            return q[0], q[1], METHOD_MOLGPKA
        except Exception:  # noqa: BLE001 - too large, unparsable for the engine ...: use the functional groups for this structure
            pass
    q = functional_group_classes(norm.mol)
    return q[0], q[1], METHOD_FG


def best_charge_method():
    return METHOD_MOLGPKA if molgpka_available()[0] else METHOD_FG


# ---------------------------------------------------------------- one structure -> everything the index stores
def describe(smiles, method, cfg=CONFIG):
    """SMILES -> (inchikey, normalised smiles, fingerprint bytes, logP, TPSA, MW, class low pH, class neutral pH, method used) or None."""
    try:
        norm = normalize(smiles)
        q = charge_classes(norm, method, cfg)
        mol = norm.mol
        return (norm.inchikey, norm.smiles, fingerprint(mol, cfg), float(Crippen.MolLogP(mol)), float(rdMolDescriptors.CalcTPSA(mol)),
                float(Descriptors.MolWt(mol)), q[0], q[1], q[2])
    except Exception:  # noqa: BLE001 - one odd structure must not stop an index build
        return None


def _worker_init():
    from rdkit import RDLogger as _R
    _R.DisableLog("rdApp.*")


def _describe_chunk(args):
    smiles_list, method = args
    return [describe(s, method) for s in smiles_list]


def describe_many(smiles_list, method, workers=1, progress=None):
    """describe() for each SMILES, in order; ``workers`` > 1 uses a process pool (MolGpKa takes ~8 ms per structure)."""
    smiles_list = list(smiles_list)
    if workers <= 1 or len(smiles_list) < 2000:
        out = []
        for i, s in enumerate(smiles_list):
            out.append(describe(s, method))
            if progress and (i + 1) % 5000 == 0:
                progress(i + 1, len(smiles_list))
        return out
    from concurrent.futures import ProcessPoolExecutor
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(var, "1")           # one BLAS thread per worker; inherited by the spawned processes
    size = 500
    chunks = [(smiles_list[i:i + size], method) for i in range(0, len(smiles_list), size)]
    out = []
    with ProcessPoolExecutor(max_workers=workers, initializer=_worker_init) as pool:
        for part in pool.map(_describe_chunk, chunks):
            out.extend(part)
            if progress:
                progress(len(out), len(smiles_list))
    return out


# ---------------------------------------------------------------- display names
_ID_LIKE = re.compile(r"^(HMDB\d+|CHEBI:\d+|C\d{5}|LM[A-Z]{2}\d+|[A-Z]{2,}\d{4,}|\d{2,7}-\d\d-\d|\d+)$", re.I)


def is_id_like(name):
    """True for a name that is just a database ID (HMDB0000123, CHEBI:1234, C00031, LMFA..., a CAS number, digits)."""
    return bool(_ID_LIKE.match(name))


def best_name(counter):
    """Display name of a compound from the names its entries carry: the most common one that is not just a database ID; ALL-CAPS words are lower-cased."""
    names = [n for n, _ in counter.most_common()]
    if not names:
        return ""
    name = next((n for n in names if not _ID_LIKE.match(n)), names[0])
    return name.capitalize() if name.isupper() and len(name) > 5 and not _ID_LIKE.match(name) else name


# ---------------------------------------------------------------- packed strings
class Packed:
    """List of UTF-8 strings as one uint8 blob + offsets (no pickling in the npz)."""

    def __init__(self, blob, off):
        self.blob, self.off = blob, off

    @classmethod
    def from_list(cls, items):
        data = [s.encode("utf-8") for s in items]
        off = np.zeros(len(data) + 1, dtype=np.int64)
        np.cumsum([len(b) for b in data], out=off[1:])
        return cls(np.frombuffer(b"".join(data), dtype=np.uint8).copy(), off)

    def __getitem__(self, i):
        return bytes(self.blob[self.off[i]:self.off[i + 1]]).decode("utf-8")

    def __len__(self):
        return len(self.off) - 1


# ---------------------------------------------------------------- the index
@dataclass
class Query:
    """The query structure with everything the score needs (same descriptors as the index rows)."""
    norm: Normalized
    fp: np.ndarray
    pop: int
    logp: float
    tpsa: float
    mw: float
    q_low: int
    q_neutral: int
    method: str


def make_query(smiles, cfg=CONFIG):
    """SMILES -> Query. The charge class is MolGpKa whenever its weights are installed (also when an index used the functional-group
    fallback), else the functional-group fallback. Raises ValueError for an unusable structure."""
    norm = normalize(smiles)
    q_low, q_neutral, used = charge_classes(norm, best_charge_method(), cfg)
    mol = norm.mol
    fp = np.frombuffer(fingerprint(mol, cfg), dtype=np.uint8)
    return Query(norm=norm, fp=fp.view(np.uint64), pop=int(popcount(fp.view(np.uint64)).sum()), logp=float(Crippen.MolLogP(mol)),
                 tpsa=float(rdMolDescriptors.CalcTPSA(mol)), mw=float(Descriptors.MolWt(mol)), q_low=q_low, q_neutral=q_neutral, method=used)


@dataclass
class Scores:
    structure: np.ndarray          # Tanimoto
    properties: np.ndarray         # mean of the property terms
    total: np.ndarray              # w_structure * structure + w_properties * properties


class CompoundIndex:
    """Structures of one source. Row i of every per-compound array is compound i; the entries (method references) are sorted by compound."""

    ARRAYS = ("ik14", "fp", "pop", "logp", "tpsa", "mw", "q_low", "q_neutral", "ent_comp", "ent_ref", "ent_rt", "ent_ik2", "ent_start")
    OPTIONAL = {"ent_role": np.int8}     # per entry: 0 = the method measures the compound, 1 = only the compound's monograph (JP18); absent = all 0

    def __init__(self, meta, arrays, smiles, names):
        self.meta = meta
        self.__dict__.update(arrays)
        self.smiles, self.names = smiles, names
        self.fp = self.fp.view(np.uint64)                 # (N, fp_bits / 64)

    def __len__(self):
        return len(self.ik14)

    @property
    def charge_method(self):
        return self.meta.get("charge_method")

    def check(self, cfg=CONFIG):
        want = {"fp": [cfg.fp_radius, cfg.fp_bits], "ph": [cfg.ph_acidic, cfg.ph_neutral]}
        got = {k: self.meta.get(k) for k in want}
        if got != want:
            raise ValueError(f"index built with {got}, the page now uses {want}: rebuild the index")

    # -- search
    def exact(self, key14):
        return np.flatnonzero(self.ik14 == key14.encode("ascii"))

    def score(self, q, cfg=CONFIG):
        tani = tanimoto(self.fp, self.pop, q.fp, q.pop)
        w = cfg.prop_weights
        terms = (w[0] * np.exp(-np.abs(self.logp - q.logp) / cfg.scale_logp),
                 w[1] * np.exp(-np.abs(self.tpsa - q.tpsa) / cfg.scale_tpsa),
                 w[2] * np.exp(-np.abs(self.mw - q.mw) / cfg.scale_mw),
                 w[3] * (self.q_low == q.q_low),
                 w[4] * (self.q_neutral == q.q_neutral))
        prop = sum(terms) / float(sum(w))
        return Scores(structure=tani, properties=prop.astype(np.float32),
                      total=(cfg.w_structure * tani + cfg.w_properties * prop).astype(np.float32))

    def candidates(self, scores, exclude_key14, cfg=CONFIG):
        """Compound rows with total >= min_score (and structure >= min_structure), best first, without the compound ``exclude_key14``."""
        ok = (scores.total >= cfg.min_score) & (scores.structure >= cfg.min_structure) & (self.ik14 != exclude_key14.encode("ascii"))
        rows = np.flatnonzero(ok)
        order = np.lexsort((self.ik14[rows], -scores.structure[rows], -scores.total[rows]))
        return rows[order]

    def entries(self, i):
        """(refs, rts, stereo blocks) of the method references of compound i."""
        a, b = self.ent_start[i], self.ent_start[i + 1]
        return self.ent_ref[a:b], self.ent_rt[a:b], self.ent_ik2[a:b]

    def roles(self, i):
        """Role of each method reference of compound i (see OPTIONAL)."""
        a, b = self.ent_start[i], self.ent_start[i + 1]
        return self.ent_role[a:b]

    # -- save / load
    def save(self, path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(path.name + ".part")
        arrays = {k: getattr(self, k) for k in self.ARRAYS + tuple(self.OPTIONAL)}
        arrays["fp"] = self.fp.view(np.uint8)
        for name, packed in (("smiles", self.smiles), ("names", self.names)):
            arrays[name + "_blob"], arrays[name + "_off"] = packed.blob, packed.off
        arrays["meta_json"] = np.frombuffer(json.dumps(self.meta).encode("utf-8"), dtype=np.uint8)
        with open(tmp, "wb") as f:
            np.savez(f, **arrays)
        os.replace(tmp, path)

    @classmethod
    def load(cls, path):
        with np.load(path, allow_pickle=False) as z:
            data = {k: z[k] for k in z.files}
        meta = json.loads(bytes(data.pop("meta_json")).decode("utf-8"))
        smiles = Packed(data.pop("smiles_blob"), data.pop("smiles_off"))
        names = Packed(data.pop("names_blob"), data.pop("names_off"))
        arrays = {k: data[k] for k in cls.ARRAYS}
        for k, dtype in cls.OPTIONAL.items():
            arrays[k] = data[k] if k in data else np.zeros(len(arrays["ent_comp"]), dtype=dtype)
        return cls(meta, arrays, smiles, names)

    # -- build
    @classmethod
    def build(cls, rows, meta, method, cfg=CONFIG, workers=1, progress=None, dedupe_refs=False):
        """Index from ``rows`` = [(smiles, ref, rt, name[, role])] (one per method reference; rt may be NaN, name may be "", role 0 or 1).

        Compounds are merged by connectivity key (the first row of a key supplies the structure). With ``dedupe_refs`` a compound has one
        entry per reference (the lowest role wins). Returns (index, number of rows that could not be processed, number of structures that
        fell back from MolGpKa to the functional groups).
        """
        uniq = list(dict.fromkeys(r[0] for r in rows if r[0]))
        described = dict(zip(uniq, describe_many(uniq, method, workers, progress)))
        comp_of, comps, skipped = {}, [], 0
        names_of = []
        refs = []
        for smiles, ref, rt, name, *rest in rows:
            d = described.get(smiles)
            if d is None:
                skipped += 1
                continue
            key = d[0][:14]
            if key not in comp_of:
                comp_of[key] = len(comps)
                comps.append(d)
                names_of.append(Counter())
            c = comp_of[key]
            if name:
                names_of[c][name] += 1
            refs.append((c, ref, rt, d[0][15:25], rest[0] if rest else 0))
        if dedupe_refs:
            best_of = {}
            for r in refs:
                if (r[0], r[1]) not in best_of or r[4] < best_of[(r[0], r[1])][4]:
                    best_of[(r[0], r[1])] = r
            refs = list(best_of.values())
        fallbacks = sum(1 for d in comps if method == METHOD_MOLGPKA and d[8] != METHOD_MOLGPKA)
        order = sorted(range(len(refs)), key=lambda i: (refs[i][0], refs[i][1]))      # entries sorted by compound, then reference
        refs = [refs[i] for i in order]
        n = len(comps)
        ent_start = np.zeros(n + 1, dtype=np.int64)
        np.cumsum(np.bincount([r[0] for r in refs], minlength=n), out=ent_start[1:])
        arrays = {
            "ik14": np.array([d[0][:14] for d in comps], dtype="S14"),
            "fp": np.frombuffer(b"".join(d[2] for d in comps), dtype=np.uint8).reshape(n, cfg.fp_bits // 8).view(np.uint64).copy(),
            "pop": np.array([int(popcount(np.frombuffer(d[2], dtype=np.uint64)).sum()) for d in comps], dtype=np.uint16),
            "logp": np.array([d[3] for d in comps], dtype=np.float32),
            "tpsa": np.array([d[4] for d in comps], dtype=np.float32),
            "mw": np.array([d[5] for d in comps], dtype=np.float32),
            "q_low": np.array([d[6] for d in comps], dtype=np.int8),
            "q_neutral": np.array([d[7] for d in comps], dtype=np.int8),
            "ent_comp": np.array([r[0] for r in refs], dtype=np.int32),
            "ent_ref": np.array([r[1] for r in refs], dtype=np.int32),
            "ent_rt": np.array([r[2] for r in refs], dtype=np.float32),
            "ent_ik2": np.array([r[3] for r in refs], dtype="S10"),
            "ent_role": np.array([r[4] for r in refs], dtype=np.int8),
            "ent_start": ent_start,
        }
        meta = dict(meta, schema=SCHEMA, charge_method=method, fp=[cfg.fp_radius, cfg.fp_bits], ph=[cfg.ph_acidic, cfg.ph_neutral],
                    n_compounds=n, n_entries=len(refs), charge_fallbacks=fallbacks)
        best = [best_name(c) for c in names_of]
        index = cls(meta, arrays, Packed.from_list([d[1] for d in comps]), Packed.from_list(best))
        return index, skipped, fallbacks


def rt_text(rts, limit=3):
    """Retention times (min) of one compound in one dataset: one value, or the distinct values (several entries)."""
    vals = sorted({round(float(r), 2) for r in rts if r is not None and math.isfinite(float(r))})
    if not vals:
        return None
    text = " / ".join(f"{v:g}" for v in vals[:limit])
    return text + (" / …" if len(vals) > limit else "")
