"""CCS（碰撞截面）页的引擎：本地索引、同化合物同加合物的比较、没有参考值时的模型预测。

数据：CCSbase 的源数据（dylanhross/c3sdb，固定提交号，MIT 许可，26 个数据集，全部在 N2 中测得，DT / TW / TIMS 三种平台）。
`python -m mzdesk.orthogonal.ccs_setup` 一次性下载、建索引、训练预测模型（ccs_setup.py）；这里读索引，回答页面。
引用：Ross, Cho, Xu. Anal. Chem. 2020, 92, 4548-4557. https://doi.org/10.1021/acs.analchem.9b05772

比较规则（用户定的，不要改）：
1. 只比较同一化合物（规范化后 InChIKey 前 14 位）的同一加合物。别的加合物只列出来作参考。
2. 只比较 N2 标度的 CCS。He 中测的值、用 He 参考值校准的 TWIMS（TWCCS N2→He）：不可比（参考数据只有 N2）。
3. FAIMS / DMS 给的是补偿电压，不是 CCS：不比较。
4. TIMS 填供应商软件报告的 CCS（Å²），不收 1/K0。
5. 没有任何参考值时，才用模型预测；预测的容差是测试集误差的 95 分位。

所有数字（容差、领域阈值、训练划分）都在 Config 里，不要散落在别处。
"""
import gzip
import json
import math
import os
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from mzdesk import adducts as adduct_names
from mzdesk.methods import compound_index as ci

SETUP_COMMAND = "python -m mzdesk.orthogonal.ccs_setup"
REPO_URL = "https://github.com/dylanhross/c3sdb"
LICENCE = "MIT"
LICENCE_URL = "https://github.com/dylanhross/c3sdb/blob/master/LICENSE.txt"
CITATION = ("D. H. Ross, J. J. Cho, L. Xu. Breaking Down Structural Diversity for Comprehensive Prediction of Ion-Neutral Collision Cross Sections. "
            "Anal. Chem. 92, 4548-4557 (2020). https://doi.org/10.1021/acs.analchem.9b05772")
INDEX_FILE = "ccs_index.json.gz"
META_FILE = "ccs_meta.json"
MODEL_FILE = "ccs_model.joblib"
DOMAIN_FILE = "ccs_domain.npz"
SCHEMA = 1                                 # 索引或特征的含义变了就加 1（旧索引会被拒绝，提示重建）
COLUMNS = ("key", "name", "smiles", "adduct", "mz", "z", "ccs", "ccs_type", "ccs_method", "src_tag")


# ---------------------------------------------------------------- 全部数字都在这里
@dataclass(frozen=True)
class Config:
    # ---- 容差（±%）：实测值 vs 参考值，N2 标度；按"观测平台 / 参考行的 ccs_type"这一对选
    tol_same: float = 1.0          # 同一类型。DTIMS 实验室间 RSD 0.29 %（Stow et al., Anal. Chem. 2017, 89, 9048），单场法偏差 0.54 %
    tol_dt_tims: float = 1.0       # DT <-> TIMS。Feuerstein et al., J. Am. Soc. Mass Spectrom. 2022（类固醇）：95 % 的离子在 ±1 % 内
    tol_dt_tw: float = 2.0         # DT <-> TW。同一研究：在 ±2 % 内
    tol_tims_tw: float = 3.0       # TIMS <-> TW。推导值：两个平台各自相对 DT 的偏差相加，没有直接研究
    # 没有参考值时的预测：容差 = 测试集绝对相对误差的 95 分位（建索引时算好，向上取整到 tol_round）
    tol_round: float = 0.1
    # ---- 预测的适用范围
    # 查询与训练集的最大 Tanimoto（Morgan r=2, 2048 位）低于它：标"训练范围之外"。只提示，不改结论。
    # 0.5 的依据（commit 5f1393e 的数据，测试集按化合物留出）：最近邻 Tanimoto < 0.5 的测试行（n=288）误差中位 2.3 %、95 分位 9.7 %；
    # ≥ 0.5（n=1779）中位 1.1 %、95 分位 5.5 %，在整体容差（6.5 %）以内。0.3 只能分出 75 行（95 分位 11 %）。重建索引后各区间的数字在 ccs_meta.json 的 model.domain.test_bins。
    domain_min_tanimoto: float = 0.5
    min_adduct_rows: int = 100           # 训练集里某加合物至少这么多行才预测；否则 "no prediction for this adduct"
    # ---- 训练划分（按连接键分组，同一化合物不会同时在训练集和测试集里）
    seed: int = 20241203
    test_fraction: float = 0.15
    val_fraction: float = 0.15
    # ---- 观测值输入
    observed_min: float = 20.0           # 小于它多半是 1/K0（Vs/cm²）或单位错了：拒绝
    observed_max: float = 1000.0
    max_observed: int = 10
    # ---- 行内提示：索引里的 m/z 与按这个结构和加合物算出的 m/z 差超过它（Da），说明来源的名字 / 结构 / 加合物可能对不上
    mz_check_da: float = 0.05


CONFIG = Config()

class QueryError(ValueError):
    """查询的输入不能用（结构、加合物）。消息直接给用户看。"""


class IndexProblem(RuntimeError):
    """索引过旧或读不了（不是用户的输入）。"""


# ---------------------------------------------------------------- 平台、气体、可比性
PLATFORMS = {"DTIMS (stepped field)": "DT", "DTIMS (single field)": "DT", "TIMS": "TIMS", "TWIMS": "TW", "FAIMS / DMS": None}
GASES = ("N₂", "He")
TW_CALIBRATIONS = ("N₂ reference values", "He reference values (N₂→He)")
TYPE_LABEL = {"DT": "DTCCS N2", "TW": "TWCCS N2", "TIMS": "TIMSCCS N2"}
TYPE_ORDER = {"DT": 0, "TIMS": 1, "TW": 2}
CONSISTENT, INCONSISTENT, NOT_COMPARABLE, NO_OBSERVED, NO_REFERENCE = "Consistent", "Inconsistent", "Not comparable", "No observed value", "No reference"


def tolerance(obs_type, ref_type, cfg=CONFIG):
    """观测平台类型 vs 参考行类型 -> 容差（±%）。"""
    if obs_type == ref_type:
        return cfg.tol_same
    pair = {obs_type, ref_type}
    if pair == {"DT", "TIMS"}:
        return cfg.tol_dt_tims
    if pair == {"DT", "TW"}:
        return cfg.tol_dt_tw
    if pair == {"TIMS", "TW"}:
        return cfg.tol_tims_tw
    raise ValueError(f"Unknown platform type pair: {obs_type!r} / {ref_type!r}")


def comparability(platform, gas, tw_calibration):
    """(能比较?, 简短原因, 完整原因)。FAIMS 最优先；然后是 He；然后是用 He 参考值校准的 TWIMS。"""
    if platform not in PLATFORMS:
        raise ValueError(f"Unknown platform: {platform!r}")
    if PLATFORMS[platform] is None:
        return False, "FAIMS / DMS gives no CCS", "FAIMS / DMS gives a compensation voltage, not a CCS, so nothing is compared."
    if gas == "He":
        return False, "measured in He", "A CCS measured in He is not comparable with N₂ values (the reference data are N₂ only)."
    if platform == "TWIMS" and tw_calibration == TW_CALIBRATIONS[1]:
        return False, "TWCCS N₂→He", "A TWIMS CCS calibrated with He reference values (TWCCS N₂→He) is not comparable with N₂ values."
    return True, "", ""


def parse_observed(raw, cfg=CONFIG):
    """"152.3, 160.1" -> [152.3, 160.1]；空 -> []。逗号、分号或空白分隔。错误抛 ValueError（消息直接给用户看）。"""
    if raw is None:
        return []
    if isinstance(raw, (int, float)):
        parts = [str(raw)]
    else:
        parts = [p for p in re.split(r"[,;\s]+", str(raw).strip()) if p]
    values = []
    for p in parts:
        try:
            v = float(p)
        except ValueError:
            raise ValueError(f"Observed CCS: '{p}' is not a number.") from None
        if not math.isfinite(v):
            raise ValueError(f"Observed CCS: '{p}' is not a finite number.")
        if v < cfg.observed_min:
            raise ValueError(f"Observed CCS: {p} is too small for a CCS in Å². For TIMS, enter the CCS your vendor software reports, not 1/K0.")
        if v > cfg.observed_max:
            raise ValueError(f"Observed CCS: {p} Å² is too large for a small molecule (limit {cfg.observed_max:g} Å²).")
        values.append(v)
    if len(values) > cfg.max_observed:
        raise ValueError(f"Observed CCS: at most {cfg.max_observed} values.")
    return values


# ---------------------------------------------------------------- 加合物
_C3SDB_FIX = {"[M+]+": "[M]+", "M+NH4]+": "[M+NH4]+", "[M+H]+*": "[M+H]+", "[M+Na]+*": "[M+Na]+", "[M+H20-H]-": "[M+H2O-H]-"}   # c3sdb 加载数据时同样修正
_SYNONYM = {"[M+H3C2O2]-": "[M+CH3COO]-"}           # H3C2O2 就是乙酸根
ELECTRON = 0.000548579909
_ADDUCT = re.compile(r"^\[(\d*)M((?:[+-]\d*[A-Z][A-Za-z0-9]*)*)\](\d*)([+-])$")
_TERM = re.compile(r"([+-])(\d*)([A-Z][A-Za-z0-9]*)")
_ELEMENT = re.compile(r"([A-Z][a-z]?)(\d*)")


def map_adduct(raw):
    """c3sdb 或用户的加合物写法 -> mzDesk 的规范写法（[M+HCOO]-、[M-H-H2O]- ……）；认不出的保持原样（去掉空白）。"""
    s = _C3SDB_FIX.get(str(raw).strip(), str(raw).strip())
    s = _SYNONYM.get(s, s)
    return adduct_names.canonical(s) or s


def adduct_charge(adduct):
    """[M+H]+ -> 1，[M-2H]2- -> -2；认不出 -> 0。"""
    m = re.search(r"\](\d*)([+-])$", adduct or "")
    return 0 if not m else (int(m.group(1) or 1)) * (1 if m.group(2) == "+" else -1)


_ADDUCT_ELEMENTS = {"H", "C", "N", "O", "F", "Li", "Na", "P", "S", "Cl", "K", "Br", "Rb", "I", "Cs"}      # 加合物里会出现的元素（"OAcO" 里的 Ac 不是锕）


def _formula_mass(formula):
    table = Chem.GetPeriodicTable()
    pos, total = 0, 0.0
    for m in _ELEMENT.finditer(formula):
        if m.start() != pos or m.group(1) not in _ADDUCT_ELEMENTS:
            return None
        pos = m.end()
        total += table.GetMostCommonIsotopeMass(m.group(1)) * int(m.group(2) or 1)
    return total if pos == len(formula) and total else None


def adduct_mz(neutral_mass, adduct):
    """中性分子精确质量 -> 加合物离子的 m/z；认不出的加合物返回 None。m/z = (n·M + Σ加减的基团 − z·电子质量) / |z|。"""
    m = _ADDUCT.match(adduct or "")
    if not m:
        return None
    n_m, terms = int(m.group(1) or 1), m.group(2)
    z = adduct_charge(adduct)
    if not z:
        return None
    total = n_m * neutral_mass
    for sign, count, formula in _TERM.findall(terms):
        mass = _formula_mass(formula)
        if mass is None:
            return None
        total += (1 if sign == "+" else -1) * int(count or 1) * mass
    return (total - z * ELECTRON) / abs(z)


# ---------------------------------------------------------------- 本地文件
def ccs_dir():
    """CCS 索引所在目录。MZDESK_CCS_DIR 可改（默认 ~/.mzdesk/data/ccs）。"""
    return Path(os.environ.get("MZDESK_CCS_DIR") or Path.home() / ".mzdesk" / "data" / "ccs")


def installed(path=None):
    d = Path(path) if path else ccs_dir()
    return (d / INDEX_FILE).is_file() and (d / META_FILE).is_file()


@dataclass(frozen=True)
class Row:
    """索引的一行：一个化合物、一个加合物、一个实测 CCS。"""
    key: str            # 连接键：规范化后 InChIKey 前 14 位
    name: str
    smiles: str         # 规范化后的 SMILES（最大片段、中性化）
    adduct: str
    mz: float
    z: int
    ccs: float
    ccs_type: str       # DT | TW | TIMS
    ccs_method: str
    src_tag: str


class CcsIndex:
    """所有行 + 数据集信息（来源、引用）+ meta。按连接键查。"""

    def __init__(self, rows, sources, meta):
        self.rows, self.sources, self.meta = rows, sources, meta
        self.by_key = defaultdict(list)
        for r in rows:
            self.by_key[r.key].append(r)

    def __len__(self):
        return len(self.rows)

    def lookup(self, key, adduct):
        """(同一加合物的行, 其他加合物的行)；行按 (类型, 数据集, CCS) 排序。"""
        rows = sorted(self.by_key.get(key, ()), key=lambda r: (TYPE_ORDER[r.ccs_type], r.src_tag, r.ccs))
        return [r for r in rows if r.adduct == adduct], [r for r in rows if r.adduct != adduct]

    def source(self, row):
        """数据集信息 {reference, url, ccs_type, ccs_method, n}。"""
        return self.sources.get(row.src_tag) or {}

    def compound_name(self, key):
        return ci.best_name(Counter(r.name for r in self.by_key.get(key, ())))


def save_index(path, rows, sources):
    """索引写成 gzip 的 JSON（列表式，体积小）；先写 .part 再改名。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".part")
    payload = {"schema": SCHEMA, "columns": list(COLUMNS), "rows": [[getattr(r, c) for c in COLUMNS] for r in rows], "sources": sources}
    with gzip.open(tmp, "wt", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, separators=(",", ":"))
    os.replace(tmp, path)


_cache = {}


def load(path=None):
    """已安装的索引（带缓存，文件变了就重读）；没装返回 None；格式过旧抛 IndexProblem。"""
    d = Path(path) if path else ccs_dir()
    if not installed(d):
        return None
    stamp = ((d / INDEX_FILE).stat().st_mtime_ns, (d / META_FILE).stat().st_mtime_ns)
    hit = _cache.get(str(d))
    if hit and hit[0] == stamp:
        return hit[1]
    with gzip.open(d / INDEX_FILE, "rt", encoding="utf-8") as f:
        payload = json.load(f)
    if payload.get("schema") != SCHEMA or payload.get("columns") != list(COLUMNS):
        raise IndexProblem(f"The CCS index in {d} has an old format. Run {SETUP_COMMAND} --force to rebuild it.")
    with open(d / META_FILE, encoding="utf-8") as f:
        meta = json.load(f)
    index = CcsIndex([Row(*r) for r in payload["rows"]], payload["sources"], meta)
    _cache[str(d)] = (stamp, index)
    return index


# ---------------------------------------------------------------- 特征（训练和预测共用）
def descriptor_vector(mol):
    """分子描述符：42 个 MQN（rdMolDescriptors.MQNs_，和 CCSbase 一样）。mol 是规范化后的中性最大片段。
    试过再加 TPSA / logP / Fsp3：验证集中位误差 1.20 % 对 1.16 %，没有好处，所以不加。"""
    return np.array(list(rdMolDescriptors.MQNs_(mol)), dtype=np.float64)


def feature_vector(desc, mz, adduct, adduct_list):
    """[m/z, 加合物 one-hot…, 描述符…]。"""
    onehot = np.array([1.0 if adduct == a else 0.0 for a in adduct_list])
    return np.concatenate([[mz], onehot, desc])


# ---------------------------------------------------------------- 预测
@dataclass
class Prediction:
    ccs: float
    tolerance_pct: float        # 测试集 95 分位，向上取整到 0.1 %
    median_pct: float           # 测试集绝对相对误差的中位数 / 95 分位 / 行数
    q95_pct: float
    n_test: int
    model: str
    max_tanimoto: float
    in_domain: bool


class Predictor:
    """训练好的模型 + 训练集指纹（判断是否在训练范围内）。estimator 只需要 predict(X)。"""

    def __init__(self, estimator, info, train_fps):
        self.estimator, self.info, self.fps = estimator, info, train_fps
        self.pops = ci.popcount(train_fps).sum(axis=1).astype(np.uint16) if len(train_fps) else np.zeros(0, dtype=np.uint16)

    @property
    def adducts(self):
        return self.info.get("adducts") or {}

    def predict(self, norm, adduct, mz, cfg=CONFIG):
        """Prediction，或 None（这个加合物训练行不够）。"""
        if adduct not in self.adducts:
            return None
        x = feature_vector(descriptor_vector(norm.mol), mz, adduct, list(self.adducts)).reshape(1, -1)
        value = float(self.estimator.predict(x)[0])
        fp = np.frombuffer(ci.fingerprint(norm.mol), dtype=np.uint8).view(np.uint64)
        sim = float(ci.tanimoto(self.fps, self.pops, fp, int(ci.popcount(fp).sum())).max()) if len(self.fps) else 0.0
        test = self.info["test"]
        return Prediction(ccs=value, tolerance_pct=self.info["tolerance_pct"], median_pct=test["median_pct"], q95_pct=test["q95_pct"], n_test=test["n"],
                          model=self.info["name"], max_tanimoto=sim, in_domain=sim >= cfg.domain_min_tanimoto)


def load_predictor(path=None, index=None):
    """已安装的预测模型，或 None（没训练 / scikit-learn 版本不一致 / 文件坏了）。第二项是不可用的原因。"""
    d = Path(path) if path else ccs_dir()
    meta = index.meta if index is not None else None
    if meta is None:
        try:
            meta = json.loads((d / META_FILE).read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return None, "The CCS index has no readable meta file."
    info = meta.get("model")
    if not info or not (d / MODEL_FILE).is_file() or not (d / DOMAIN_FILE).is_file():
        return None, f"No prediction model is installed. Run {SETUP_COMMAND} --force."
    key = (str(d), (d / MODEL_FILE).stat().st_mtime_ns)
    if key in _models:
        return _models[key], ""
    try:
        import joblib
        import sklearn
        if info.get("sklearn") != sklearn.__version__:
            return None, (f"The prediction model was trained with scikit-learn {info.get('sklearn')}, this is {sklearn.__version__}. "
                          f"Run {SETUP_COMMAND} --force to retrain it.")
        estimator = joblib.load(d / MODEL_FILE)           # 本机训练生成的文件，不是下载来的
        with np.load(d / DOMAIN_FILE, allow_pickle=False) as z:
            fps = z["fps"].view(np.uint64)
    except Exception as exc:  # noqa: BLE001 - 缺 scikit-learn、文件损坏：页面照常工作，只是没有预测
        return None, f"The prediction model could not be loaded ({type(exc).__name__}). Run {SETUP_COMMAND} --force."
    _models[key] = Predictor(estimator, info, fps)
    return _models[key], ""


_models = {}


# ---------------------------------------------------------------- 一次查询
@dataclass
class Comparison:
    """一个参考值（实测行或预测）和观测值的比较。"""
    ccs: float
    ref_type: str | None          # DT / TW / TIMS；预测为 None
    tolerance: float | None
    observed: float | None        # 最接近参考值的那个观测值
    delta_pct: float | None       # (观测 − 参考) / 参考 × 100
    verdict: str
    row: Row | None = None


@dataclass
class Result:
    installed: bool
    key: str = ""
    name: str = ""
    smiles: str = ""
    adduct: str = ""
    mz: float | None = None
    observed: list = field(default_factory=list)
    platform: str = ""
    gas: str = "N₂"
    tw_calibration: str = ""
    obs_type: str | None = None
    comparable: bool = True
    reason: str = ""                                      # 不可比的完整原因 / 简短原因
    reason_short: str = ""
    rows: list = field(default_factory=list)              # 同化合物同加合物的参考行
    others: list = field(default_factory=list)            # 同化合物的其他加合物
    comparisons: list = field(default_factory=list)       # 与 rows 一一对应
    prediction: Prediction | None = None
    prediction_note: str = ""                             # 为什么没有预测 / 为什么没试
    prediction_comparison: Comparison | None = None
    used: Comparison | None = None                        # 汇总卡片用的那一个
    ref_ccs: float | None = None
    source: str = ""                                      # "measured" | "predicted" | ""
    verdict: str = NO_REFERENCE
    meta: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)             # 规范化时的提示（盐、中性化）
    index: CcsIndex | None = None
    cfg: Config = CONFIG


def mz_mismatch(row, cfg=CONFIG):
    """(索引里的 m/z, 按结构和加合物算出的 m/z)，两者差超过 mz_check_da 时返回；一致或算不出返回 None。"""
    mol = Chem.MolFromSmiles(row.smiles)
    calc = adduct_mz(Descriptors.ExactMolWt(mol), row.adduct) if mol is not None else None
    return (row.mz, calc) if calc is not None and abs(row.mz - calc) > cfg.mz_check_da else None


def compare(ccs, ref_type, obs_type, observed, tol, row=None):
    """一个参考值 vs 一组观测值：取最接近的观测值；|Δ| ≤ 容差 才算一致。tol 为 None 表示这一对不可比。"""
    if not observed:
        return Comparison(ccs, ref_type, tol, None, None, NO_OBSERVED, row)
    best = min(observed, key=lambda v: abs(v - ccs))
    delta = (best - ccs) / ccs * 100.0
    if tol is None:
        return Comparison(ccs, ref_type, None, best, None, NOT_COMPARABLE, row)
    return Comparison(ccs, ref_type, tol, best, delta, CONSISTENT if abs(delta) <= tol + 1e-9 else INCONSISTENT, row)


def median(values):
    s = sorted(values)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2


def analyse(smiles, adduct, observed, platform, gas="N₂", tw_calibration=TW_CALIBRATIONS[0], *, cfg=CONFIG, path=None):
    """查询：结构 + 加合物 + 观测 CCS（列表，可空）+ 平台 / 气体 -> Result。

    Raises QueryError：结构不能用、加合物认不出；IndexProblem：索引过旧。索引没装：Result.installed = False。
    """
    adduct = map_adduct(adduct)
    if not adduct_mz(100.0, adduct):
        raise QueryError(f"Unknown adduct: {adduct!r}.")
    try:
        norm = ci.normalize(smiles)
    except ValueError as exc:
        raise QueryError(f"Structure: {exc}") from None
    comparable, short, reason = comparability(platform, gas, tw_calibration)
    res = Result(installed=False, key=norm.key14, smiles=norm.smiles, adduct=adduct, observed=list(observed), platform=platform, gas=gas,
                 tw_calibration=tw_calibration, obs_type=PLATFORMS[platform], comparable=comparable, reason=reason, reason_short=short,
                 notes=list(norm.notes), cfg=cfg)
    res.mz = adduct_mz(Descriptors.ExactMolWt(norm.mol), adduct)
    index = load(path)
    if index is None:
        return res
    res.installed, res.index, res.meta = True, index, index.meta
    res.rows, res.others = index.lookup(norm.key14, adduct)
    res.name = index.compound_name(norm.key14)

    def tol_for(ref_type):
        return tolerance(res.obs_type, ref_type, cfg) if comparable else None

    res.comparisons = [compare(r.ccs, r.ccs_type, res.obs_type, res.observed, tol_for(r.ccs_type), r) for r in res.rows]
    if res.rows:
        res.source = "measured"
        pool = [c for c in res.comparisons if c.verdict in (CONSISTENT, INCONSISTENT)]
        if pool:                                           # 与观测值比较时，用相对容差最接近的参考值（任一参考值一致就算一致）
            res.used = min(pool, key=lambda c: (abs(c.delta_pct) / c.tolerance, abs(c.delta_pct)))
            res.ref_ccs = res.used.ccs
        else:
            res.ref_ccs = median([r.ccs for r in res.rows])
    else:
        predictor, why = load_predictor(path, index)
        if predictor is None:
            res.prediction_note = why
        else:
            res.prediction = predictor.predict(norm, adduct, res.mz, cfg)
            if res.prediction is None:
                res.prediction_note = f"No prediction for this adduct (fewer than {cfg.min_adduct_rows} training rows)."
        if res.prediction is not None:
            p = res.prediction
            res.source, res.ref_ccs = "predicted", p.ccs
            tol = p.tolerance_pct if comparable else None
            res.prediction_comparison = compare(p.ccs, None, res.obs_type, res.observed, tol)
            res.used = res.prediction_comparison if res.prediction_comparison.verdict in (CONSISTENT, INCONSISTENT) else None
    if res.ref_ccs is None:
        res.verdict = NO_REFERENCE
    elif not res.observed:
        res.verdict = NO_OBSERVED
    elif not comparable:
        res.verdict = NOT_COMPARABLE
    else:
        res.verdict = res.used.verdict
    return res
