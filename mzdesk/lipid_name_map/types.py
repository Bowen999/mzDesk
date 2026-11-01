"""共享数据结构（SYNTHESIS §2.1 + §3.1，W0 冻结；改动须架构师批准）。

所有模块只通过本文件的数据结构互相通信，互不 import 对方内部函数。全部是 ``frozen`` dataclass
（``slots=True``）；``dataclasses.replace`` 可用。运行时只依赖标准库。

数据流（谁填、谁读）
--------------------
M2 ``Segment`` -> M5a ``split_head`` -> M4a ``HeadInfo(rows: AliasRow)`` -> M5b ``Lexed(chains: RawChain,
feats: Features)`` -> M4b ``ClassRes`` -> M6 ``Interp(alts: Struct(chains: Chain))`` -> M7 ``Labels``。
``ClassInfo`` / ``AliasRow`` 是类表 CSV 的行（A5 的 loader 产出）；``ClassRes`` 是选定后的类。

常量（本文件定义，各模块读取，不得各自再抄一份）
-------------------------------------------------
``TOKEN_ORDER`` ``TOKEN_VALENCE``      L3 链 token 的输出顺序与氧价（COOH/OOH/OO=2，OH/oxo/Ep/O=1，Me/cPr/P=0）
``ATOM_KINDS`` ``CONJ_NAMES``          ``Atom.kind`` 取值
``OXYLIPIN_ABBREVS``                   7 个介质缩写（词法识别用；展开由语义层读 ``data/oxylipin.csv``）
``FLAGS`` ``FLAG_GROUPS``              flags 词表（SYNTHESIS §1.5）
``ALIAS_COLUMNS`` ``CLASS_INFO_COLUMNS``  两种 CSV 的列序

关键语义约定（写进契约，避免各模块各自理解）
--------------------------------------------
1. **DB 口径**：``RawChain.DB`` 是名字里**字面写的**双键数；``Chain.DB`` 是**标签口径**
   （烯基醚 ``P-`` 已 +1；DB = C=C + 环 + C=O 的并集口径，见 label_rules §5.1）。字面 -> 标签口径的换算只在
   ``semantics/conventions.py`` 一处做。
2. **未分型氧**累加进 ``Struct.tail_o``，不进 ``Chain.toks``；``Chain.toks`` 只含已分型 token
   （``TOKEN_ORDER`` 去掉 ``P``；``P`` 由 ``Chain.ether == "P"`` 表达）。
3. **介质缩写链**（``RawChain.oxylipin != ""``，如 ``12-HETE``）：词法只做识别，**不展开**——
   ``C == DB == 0``、``atoms == ()``、``placeholder == False``，位置前缀只留在 ``raw`` 里；
   ``db_ge_c`` 检查对其跳过。展开（读 ``data/oxylipin.csv``、检查 ``Config.allow_oxylipin_abbrev``、
   打 ``oxylipin_expanded`` flag）由语义层（B0/B4）做。
4. **括号 O-酰基前缀**（``AHexCer (O-24:6)18:1…`` 的 ``(O-24:6)``）放 ``Lexed.pre_acyl``，其 ``pre == ""``
   （它不是醚；不进 ``Features.ether_marks``）。
5. **is_sum 是纯句法**：恰好写了一个链组（没有任何分隔符）。它是不是"和式"取决于类的链数，由语义层决定
   （单链类 -> 链，多链类 -> 和式，PNNL 例外 D4）。
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Final

from .reasons import REASONS, describe, module_of, status_of

__all__ = [
    "Reject", "Segment", "Atom", "RawChain", "Features", "Lexed", "AliasRow", "HeadInfo", "ClassInfo",
    "ClassRes", "Chain", "Struct", "Interp", "Labels",
    "TOKEN_ORDER", "TOKEN_VALENCE", "ATOM_KINDS", "CONJ_NAMES", "AA_CONJ_NAMES", "DOUBLE_CONJ", "TAIL_CONJ_VALUES", "conj_key",
    "OXYLIPIN_ABBREVS", "LIPIDMATCH_OX_MODIFIERS", "SEGMENT_KINDS",
    "PRE_VALUES", "BASE_VALUES", "SUF_VALUES", "ETHER_MARKS", "CLASS_SUFFIXES", "FORMS", "FAMILIES",
    "IMPLIED_MODES", "SLASH_ROLES", "OXY_POLICIES", "ALIAS_STATUSES", "ETHER_KINDS", "ETHER_LINKAGES",
    "FLAGS", "FLAG_GROUPS", "ALIAS_COLUMNS", "CLASS_INFO_COLUMNS",
]

# ==========================================================================
# 常量
# ==========================================================================
#: L3 链上 token 的输出顺序（SYNTHESIS §4.1-7）
TOKEN_ORDER: Final = ("OH", "oxo", "Ep", "COOH", "OOH", "OO", "O", "Me", "cPr", "P")
#: 每个 token 的氧价（X2：COOH/OOH/OO=2，OH/oxo/Ep/O=1，Me/cPr/P=0）
TOKEN_VALENCE: Final = MappingProxyType(
    {"OH": 1, "oxo": 1, "Ep": 1, "COOH": 2, "OOH": 2, "OO": 2, "O": 1, "Me": 0, "cPr": 0, "P": 0})
_CHAIN_TOKS: Final = frozenset(TOKEN_ORDER) - {"P"}          # Chain.toks 允许的键

#: 氨基酸共轭名（0.4.0，词表 v0.2.11 的 ``ST-<AA>``：LIPID MAPS 写法 ``ST 24:1;O4;Phe``，数据集谱库名同形）。``Ile`` / ``Leu`` 是旧拼写，
#: 与 ``Xle`` 同归 ``ST-Xle``（数据集 v0.2.8 把两者并成 IUPAC 的 "Leu 或 Ile" 符号）
AA_CONJ_NAMES: Final = ("Ala", "Arg", "Asn", "Asp", "Cit", "Cys", "Gln", "Glu", "His", "Ile", "Leu", "Lys", "Met", "Orn",
                        "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "Xle")
#: 共轭名（``;Gly ;Tau ;S ;Hex ;G ;T`` 与 ``AA_CONJ_NAMES``）；``Atom.kind == "conj_" + name``
CONJ_NAMES: Final = ("Gly", "Tau", "S", "Hex", "G", "T") + AA_CONJ_NAMES
#: 双共轭（D17，如 ``ST 24:1;O3;G;S`` = ST-Gly-sulf、``ST 24:1;O3;T;S`` = ST-Tau-sulf）：``Features.tail_conj`` 的合法取值
DOUBLE_CONJ: Final = ("G+S", "T+S")
#: ``Features.tail_conj`` 的全部合法非空取值：单共轭 ``CONJ_NAMES`` 或双共轭 ``DOUBLE_CONJ``
TAIL_CONJ_VALUES: Final = CONJ_NAMES + DOUBLE_CONJ
#: ``Atom.kind`` 的全部取值，语义见 ``Atom``
ATOM_KINDS: Final = frozenset(
    {"O_count", "OH_pos", "OH", "oxo", "Ep", "COOH", "OOH", "OO", "O_other", "Me", "cPr", "P"}
    | {"conj_" + c for c in CONJ_NAMES})

#: LipidMatch 氧化磷脂链后括号修饰 -> 对应的 ``Atom.kind``（``Atom.via`` 的合法取值）。``Ke`` 酮、``CHO`` 醛都收成 ``oxo``（醛基也是 C=O，+O −2H），
#: ``COOH`` 羧基（+2O −2H）。库 m/z 核对见 docs/INPUT_FIXES.md
LIPIDMATCH_OX_MODIFIERS: Final = MappingProxyType({"Ke": "oxo", "CHO": "oxo", "COOH": "COOH"})

#: Riken 介质缩写（封闭白名单，C18）
OXYLIPIN_ABBREVS: Final = ("HETE", "HEPE", "HDoHE", "HODE", "EET", "EpETE", "EpDPE")

SEGMENT_KINDS: Final = ("plain", "lipidin_chain", "alias")
PRE_VALUES: Final = ("", "O-", "P-", "dO-")
BASE_VALUES: Final = ("", "d", "t", "m")
SUF_VALUES: Final = ("", "e", "p")
ETHER_MARKS: Final = ("O-", "P-", "e", "p", "dO-")
CLASS_SUFFIXES: Final = ("", "-O", "-P", "-d", "-t", "-m")
FORMS: Final = ("space", "paren")
FAMILIES: Final = ("gpl", "gl", "sl", "fa")
IMPLIED_MODES: Final = ("none", "add", "max")
SLASH_ROLES: Final = ("none", "base_first", "sterol_first", "hfa_second")
OXY_POLICIES: Final = ("typed_only", "sphingo_OH", "class_OH")
ALIAS_STATUSES: Final = ("verified", "inferred", "unverified", "reject")
ETHER_KINDS: Final = ("", "O", "P")                         # Chain.ether
ETHER_LINKAGES: Final = ("P", "O", "unspecified", "")       # meta.ether_linkage


def conj_key(names: Iterable[str]) -> str | None:
    """一个名字里出现的共轭名集合 -> ``Features.tail_conj`` 的取值（D17）。

    * 空 -> ``""``；一种 -> 该名原样（``Gly`` 与 ``G`` 是两种写法，单共轭时各自保留）
    * 两种且恰为 {甘氨酸, 硫酸} / {牛磺酸, 硫酸} -> ``"G+S"`` / ``"T+S"``（与书写顺序无关；``Gly``=``G``、``Tau``=``T``）
    * 其余组合（三种、重复写法 ``Gly;G``、``Hex;S`` …）-> ``None``（词法抛 ``malformed``，``Lexed.problems()`` 报告）

    >>> conj_key([]), conj_key(["Gly"]), conj_key(["S", "G"]), conj_key(["Tau", "S"]), conj_key(["Hex", "S"])
    ('', 'Gly', 'G+S', 'T+S', None)
    """
    got = set(names)
    if len(got) <= 1:
        return next(iter(got), "")
    canon = {"Gly": "G", "Tau": "T"}
    c = {canon.get(n, n) for n in got}
    if len(got) == 2 and len(c) == 2 and "S" in c:
        other = next(iter(c - {"S"}))
        if other in ("G", "T"):
            return other + "+S"
    return None

#: SYNTHESIS §1.5 flags 词表（分组）；家族可另行 ``contract.register_family_flags`` 注册 ``FLAGS_<fam>``
FLAG_GROUPS: Final = MappingProxyType({
    "decoration": ("adduct_stripped", "alias_split", "lipidin_chain_used", "lipidin_sum_O_mismatch",
                   "sn_suffix_stripped", "ox_prefix_ignored", "branch_dropped", "placeholder_dropped"),
    "convention": ("P_db_plus1", "implied_O", "implied_sterol_chain", "subclass_implied_O", "pnnl_suffix_O",
                   "lipidblast_eos_rule", "oxylipin_expanded", "source_hint_lyso"),
    "reach_limit": ("sum_only", "ether_chain_unknown", "p_landing_ambiguous", "base_chain_unknown",
                    "hfa_chain_unknown", "sterol_chain_unknown", "oxygen_untyped"),
    "hypothesis": ("oxygen_assumed_OH", "class_implied_OH", "alias_inferred",
                   "harmonized_name_oxygen_placement", "convention_unverified", "base_oxygen_assumed"),
    "downgrade": ("selfcheck_downgraded_L3", "selfcheck_downgraded_L2"),
})
FLAGS: Final = frozenset(f for grp in FLAG_GROUPS.values() for f in grp)

#: alias CSV 的列序（SYNTHESIS §3.1）
ALIAS_COLUMNS: Final = ("head", "ether", "zero", "fa", "when", "class_full", "status", "reject_reason",
                        "inject", "dialects", "evidence", "note")
#: class_info CSV 的列序（SYNTHESIS §3.1）
CLASS_INFO_COLUMNS: Final = ("class_full", "nucleus", "mod", "family", "nchain", "ether", "lyso", "implied_O",
                             "implied_mode", "slash_role", "oxy_policy", "pnnl_O", "notes")


# ==========================================================================
# 异常
# ==========================================================================
class Reject(Exception):
    """解析中途的拒绝/不支持信号。只传 ``code`` 与 ``detail``；``status`` 由 ``reasons`` 查表。

    >>> r = Reject("chain_count", "TG needs 3 chains, got 2")
    >>> (r.code, r.status, r.module)
    ('chain_count', 'rejected', 'M6')

    未知 ``code`` 在构造时抛 ``ValueError``（即拼写错误当场暴露；``normalize`` 会把它当 ``internal_error``）。
    可 pickle（``args == (code, detail)``）。
    """

    def __init__(self, code: str, detail: str = "") -> None:
        if not isinstance(code, str) or code not in REASONS:
            raise ValueError(f"unknown reason code {code!r}")
        super().__init__(code, detail)
        self.code: str = code
        self.detail: str = detail

    @property
    def status(self) -> str:
        """``"rejected"`` | ``"unsupported"``（reasons 查表）"""
        return status_of(self.code)

    @property
    def module(self) -> str:
        return module_of(self.code)

    @property
    def description(self) -> str:
        return describe(self.code)

    def __str__(self) -> str:
        return f"{self.code}: {self.detail}" if self.detail else self.code


# ==========================================================================
# M2 / M5：词法层
# ==========================================================================
@dataclass(frozen=True, slots=True)
class Segment:
    """M2 输出：一个别名段（已剥加合物/装饰）。

    text     剥离后的文本，交给 M3/M5a
    adduct   被剥掉的加合物原文（无则 ``""``）
    kind     ``"plain"`` | ``"lipidin_chain"``（取自 ``X | Chain: Y`` 的 Y 段）| ``"alias"``（``A | B`` 的一段）
    notes    被剥掉的其他装饰（诊断用；每项 ``"<类别>:<原文>"``，如 ``("sn_suffix:SN1",)``、``("ox_prefix:Ox",)``；
             flag 由 ``decorations.segment_flags(seg)`` 取）；进 ``meta.stripped = list(seg.notes)``
    """
    text: str
    adduct: str = ""
    kind: str = "plain"
    notes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.kind not in SEGMENT_KINDS:
            raise ValueError(f"Segment.kind must be one of {SEGMENT_KINDS}, got {self.kind!r}")


@dataclass(frozen=True, slots=True)
class Atom:
    """链上一个氧/修饰记号（词法产出，语义层消费）。``n`` 的含义依 ``kind``：

    =============  ==============================================================================
    kind           含义（``n`` = 个数）
    =============  ==============================================================================
    ``O_count``    **未分型**氧个数：``;O``(1) ``;On``/``;nO``/``;O2``/``;2O``(n)；LipidMatch ``(OO)(OOH)(OOO)``
                   等括号氧（D10，只计数不分型）；``(nO)``
    ``OH_pos``     位置式羟基 ``;(2OH)``：每个括号项 **1 个** OH（``n`` = 项数，``;(9OH,13OH)`` -> 2）
    ``OH`` ``oxo`` ``Ep`` ``COOH`` ``OOH`` ``OO``
                   已分型 token 的个数：``;OH`` -> 1，``;OH2``/``;(OH)2`` -> 2，LipidMatch ``(OH)`` -> 1
    ``O_other``    名字明确写成数据集通用 ``O`` token 的氧（v1 词法没有已知来源，保留给语义层映射到 ``toks["O"]``）
    ``Me`` ``cPr`` 甲基/环丙基支链（价 0，默认不输出）
    ``P``          ``;P`` 缩醛标记（数据集原生 L3 写法；其 DB 已是标签口径）
    ``conj_<X>``   共轭：``conj_Gly conj_Tau conj_S conj_Hex conj_G conj_T`` 与 ``conj_<氨基酸>``（``AA_CONJ_NAMES``；``n`` 恒 1）
    =============  ==============================================================================

    注意 ``;(2OH)``(位置，1 个) ≠ ``;(O2)``(个数 2) ≠ ``;(OH)2``(个数 2)（S13）。

    ``via``：LipidMatch 氧化磷脂的链后括号修饰 ``(Ke)`` ``(CHO)`` ``(COOH)``（2026-10-31）写成什么。非空时 ``kind`` 是 ``oxo``（Ke、CHO）或 ``COOH``，
    且**这条链的字面 DB 不含这个 C=O**（``18:1(Ke)`` 是 1 个 C=C 加 1 个酮基）；标签口径的 DB 含 C=O（naming_spec_v2「双键数」），
    所以语义层要给每个这样的 C=O 加 1 个 DB（``18:1(Ke)`` -> ``18:2;oxo``）。原生写法 ``18:2;oxo`` 的 ``via`` 恒为 ``""``，DB 已是标签口径。
    """
    kind: str
    n: int = 1
    via: str = ""

    def __post_init__(self) -> None:
        if self.kind not in ATOM_KINDS:
            raise ValueError(f"unknown Atom.kind {self.kind!r}")
        if self.n < 1:
            raise ValueError(f"Atom.n must be >= 1, got {self.n}")
        if self.via and (self.via not in LIPIDMATCH_OX_MODIFIERS or LIPIDMATCH_OX_MODIFIERS[self.via] != self.kind):
            raise ValueError(f"Atom.via {self.via!r} does not match kind {self.kind!r} (allowed: {dict(LIPIDMATCH_OX_MODIFIERS)})")


@dataclass(frozen=True, slots=True)
class RawChain:
    """M5 输出：一条链的纯词法结果，**无类语义**（不知道哪条是醚、哪条是碱基）。

    raw          该链在输入中的原文片段（含前后缀与 ``;`` 记号；诊断用）
    pre          链前缀 ``""`` | ``"O-"`` | ``"P-"`` | ``"dO-"``
    base         鞘氨醇碱基前缀 ``""`` | ``"d"`` | ``"t"`` | ``"m"``
    C, DB        **字面**碳数、双键数（见模块文档"DB 口径"）；``0:0`` 占位为 ``0, 0``
    suf          链后缀 ``""`` | ``"e"``（D5：等价 ``O-``）| ``"p"``（零实例 -> 语义层 ``unsupported_syntax``）
    atoms        按书写顺序的氧/修饰记号（``Atom``）
    n_acyl       链带 ``N-`` 前缀（LNAPE ``N-9:0``）
    branch       支链记号原文（``"iso"`` ``"anteiso"`` ``"15me"`` ``"methyl"`` …）；非空 ⇒ C、DB 不变、不输出 ``;Me``（S17）
    oxylipin     介质缩写（``OXYLIPIN_ABBREVS`` 之一）或 ``""``；非空时 ``C==DB==0`` 且 ``atoms==()``（模块文档第 3 条）
    placeholder  是 ``0:0`` 占位链
    """
    raw: str
    pre: str
    base: str
    C: int
    DB: int
    suf: str
    atoms: tuple[Atom, ...] = ()
    n_acyl: bool = False
    branch: str = ""
    oxylipin: str = ""
    placeholder: bool = False

    def __post_init__(self) -> None:
        if self.pre not in PRE_VALUES:
            raise ValueError(f"RawChain.pre must be one of {PRE_VALUES}, got {self.pre!r}")
        if self.base not in BASE_VALUES:
            raise ValueError(f"RawChain.base must be one of {BASE_VALUES}, got {self.base!r}")
        if self.suf not in SUF_VALUES:
            raise ValueError(f"RawChain.suf must be one of {SUF_VALUES}, got {self.suf!r}")
        if self.C < 0 or self.DB < 0:
            raise ValueError(f"RawChain C/DB must be >= 0, got {self.C}:{self.DB}")
        if self.oxylipin and self.oxylipin not in OXYLIPIN_ABBREVS:
            raise ValueError(f"RawChain.oxylipin must be one of {OXYLIPIN_ABBREVS}, got {self.oxylipin!r}")


@dataclass(frozen=True, slots=True)
class Features:
    """词法层对整个名字做的特征摘要。**A3 词法负责填，A5 的命名谓词 / `select_class` 负责读。**

    has_ether_prefix  任一主链的 ``pre`` ∈ {``O-``, ``P-``, ``dO-``}（和式 ``O-34:1`` 也算）。
                      仅指**链前缀**；类头里的 ``Plasmanyl-``/``Ether`` 属于别名表 head，不在这里
    ether_marks       每个带醚标记的主链各贡献其标记，按书写顺序；取值 ``'O-' 'P-' 'e' 'p' 'dO-'``
                      （一条链同时有前缀和后缀则两个都记）。``has_ether_prefix`` 必须与之一致。
                      不含 ``pre_acyl``（括号 O-酰基不是醚）
    n_placeholder     主链中 ``0:0`` 占位链的个数
    has_fa_paren      有嵌套 ``(FA c:d)``（即 ``Lexed.nested_fa`` 非空）
    has_pre_acyl      有括号 O-酰基前缀（``Lexed.pre_acyl`` 非 None）
    class_suffix      头部末尾的类后缀，纯词法：``head`` 以 ``-O -P -d -t -m`` 结尾（且去掉后非空）则取之，否则 ``""``。
                      是不是真后缀由 ``select_class`` 结合别名表决定（``SQDG-O`` 这类整体就是类名的情况由
                      ``resolve_head`` 先按整个 head 查表）
    sphingoid_prefix  第一条带 ``base`` 的主链的 ``base``（``'d' 't' 'm'``），无则 ``""``
    eos2O             LipidBlast EOS 判据：第 2 条主链（占位链不计）原文里有**数字在前**的 ``;2O`` 记号
                      （``HexCer 14:1;2O/58:1;2O``）。LM 写法 ``;O2`` 不算（那是字面 2 个 OH，S6）。
                      **纯词法判据，对任何类头都可能为真**（``PC 18:0_18:0;2O`` 也是）——只有与 ``Cer`` / ``HexCer`` 头联用
                      才有 EOS 含义（D18）：别名表里 ``HexCer`` + ``eos2O`` -> ``Cer-Hex-oacyl``，``Cer`` + ``eos2O`` -> 拒绝
                      ``unsupported_subclass``；其他类头不读它
    tail_conj         共轭名：``'' 'Gly' 'Tau' 'S' 'Hex' 'G' 'T'``、氨基酸三字母 ``AA_CONJ_NAMES``（``ST 24:1;O4;Phe``）或双共轭 ``'G+S' 'T+S'``（D17；``ST 24:1;O3;G;S`` =
                      ST-Gly-sulf、``ST 24:1;O3;T;S`` = ST-Tau-sulf）。取自主链 ``atoms`` 与 ``tail_atoms`` 中出现的 ``conj_*``，
                      与书写顺序无关；``Gly``/``G``、``Tau``/``T`` 在双共轭里归一为 ``G``/``T``。``conj_key`` 是唯一的换算处；
                      其余组合（三种、``Hex;S``、``Gly;G``）词法抛 ``malformed``
    n_chains          主链数，**不含**占位链、``pre_acyl``、``nested_fa``（含 ``n_acyl`` 链）
    is_sum            纯句法：恰好写了一个链组且无分隔符（等价 ``len(Lexed.chains) == 1``）
    seps              主链之间出现过的分隔符集合，⊆ ``{'/', '_'}``
    """
    has_ether_prefix: bool = False
    ether_marks: tuple[str, ...] = ()
    n_placeholder: int = 0
    has_fa_paren: bool = False
    has_pre_acyl: bool = False
    class_suffix: str = ""
    sphingoid_prefix: str = ""
    eos2O: bool = False
    tail_conj: str = ""
    n_chains: int = 0
    is_sum: bool = False
    seps: frozenset[str] = frozenset()

    def __post_init__(self) -> None:
        bad = [m for m in self.ether_marks if m not in ETHER_MARKS]
        if bad:
            raise ValueError(f"Features.ether_marks has unknown marks {bad}")
        has_prefix = any(m in ("O-", "P-", "dO-") for m in self.ether_marks)
        if has_prefix != self.has_ether_prefix:
            raise ValueError("Features.has_ether_prefix inconsistent with ether_marks")
        if self.class_suffix not in CLASS_SUFFIXES:
            raise ValueError(f"Features.class_suffix must be one of {CLASS_SUFFIXES}, got {self.class_suffix!r}")
        if self.sphingoid_prefix not in BASE_VALUES:
            raise ValueError(f"Features.sphingoid_prefix must be one of {BASE_VALUES}")
        if self.tail_conj and self.tail_conj not in TAIL_CONJ_VALUES:
            raise ValueError(f"Features.tail_conj must be '' or one of {TAIL_CONJ_VALUES}")
        if self.n_chains < 0 or self.n_placeholder < 0:
            raise ValueError("Features counts must be >= 0")
        if not self.seps <= {"/", "_"}:
            raise ValueError(f"Features.seps must be a subset of {{'/', '_'}}, got {set(self.seps)}")

    @property
    def has_ether(self) -> bool:
        """任意醚标记（前缀或 ``e``/``p`` 后缀）。别名表 ``ether`` 列与谓词 ``has_ether`` 用这个。"""
        return bool(self.ether_marks)


@dataclass(frozen=True, slots=True)
class Lexed:
    """M5b 输出：整个名字的词法结果。

    head       头部原文（``split_head`` 的第一返回值，**含**类后缀，如 ``"PE-O"``、``"Hex(4)-HexNAc(2)-Cer"``）
    form       ``"space"``（``PC 16:0/18:1``）| ``"paren"``（``PC(16:0/18:1)``）
    chains     主链，按书写顺序，**含** ``0:0`` 占位链与 ``n_acyl`` 链；不含 ``pre_acyl`` / ``nested_fa``
    seps       相邻主链之间的分隔符 ``'/'`` 或 ``'_'``；``len(seps) == max(len(chains) - 1, 0)``
    pre_acyl   括号 O-酰基前缀（``AHexCer (O-24:6)18:1…``），其 ``pre == ""``；无则 None
    nested_fa  嵌套 ``(FA c:d)`` 链
    tail_atoms 位于最后一条链之后的记号（LipidMatch 氧括号、共轭 ``;G`` 等）
    is_sum     同 ``feats.is_sum``（便捷副本，必须一致）
    feats      ``Features``
    """
    head: str
    form: str
    chains: tuple[RawChain, ...]
    seps: tuple[str, ...] = ()
    pre_acyl: RawChain | None = None
    nested_fa: tuple[RawChain, ...] = ()
    tail_atoms: tuple[Atom, ...] = ()
    is_sum: bool = False
    feats: Features = field(default_factory=Features)

    def __post_init__(self) -> None:
        if self.form not in FORMS:
            raise ValueError(f"Lexed.form must be one of {FORMS}, got {self.form!r}")
        if len(self.seps) != max(len(self.chains) - 1, 0):
            raise ValueError(f"Lexed.seps length {len(self.seps)} != chains-1 ({len(self.chains)} chains)")
        if self.is_sum != self.feats.is_sum:
            raise ValueError("Lexed.is_sum must equal Lexed.feats.is_sum")

    def problems(self) -> list[str]:
        """核对 ``feats`` 与 ``chains/seps/...`` 是否一致，返回问题清单（空 = 一致）。A3 的测试对每个输出调用它。"""
        out: list[str] = []
        f = self.feats
        n_ph = sum(1 for c in self.chains if c.placeholder)
        if f.n_placeholder != n_ph:
            out.append(f"n_placeholder {f.n_placeholder} != {n_ph}")
        if f.n_chains != len(self.chains) - n_ph:
            out.append(f"n_chains {f.n_chains} != {len(self.chains) - n_ph}")
        if f.has_pre_acyl != (self.pre_acyl is not None):
            out.append("has_pre_acyl inconsistent with pre_acyl")
        if f.has_fa_paren != bool(self.nested_fa):
            out.append("has_fa_paren inconsistent with nested_fa")
        marks: list[str] = []
        for c in self.chains:
            if c.pre:
                marks.append(c.pre)
            if c.suf:
                marks.append(c.suf)
        if tuple(marks) != f.ether_marks:
            out.append(f"ether_marks {f.ether_marks} != derived {tuple(marks)}")
        first_base = next((c.base for c in self.chains if c.base), "")
        if f.sphingoid_prefix != first_base:
            out.append(f"sphingoid_prefix {f.sphingoid_prefix!r} != {first_base!r}")
        conj = {a.kind[5:] for c in self.chains for a in c.atoms if a.kind.startswith("conj_")}
        conj |= {a.kind[5:] for a in self.tail_atoms if a.kind.startswith("conj_")}
        want_conj = conj_key(conj)
        if want_conj is None:
            out.append(f"unsupported conjugate combination {sorted(conj)}")
        elif f.tail_conj != want_conj:
            out.append(f"tail_conj {f.tail_conj!r} != {want_conj!r} (conjugates {sorted(conj)})")
        if f.seps != frozenset(self.seps):
            out.append(f"feats.seps {set(f.seps)} != set(seps) {set(self.seps)}")
        if f.is_sum != (len(self.chains) == 1):
            out.append(f"is_sum {f.is_sum} but {len(self.chains)} chain groups")
        for c in self.chains:
            if c.placeholder and (c.C, c.DB) != (0, 0):
                out.append(f"placeholder chain {c.raw!r} is not 0:0")
        return out


# ==========================================================================
# M4：类表
# ==========================================================================
def _parse_set(cell: str) -> frozenset[int]:
    cell = cell.strip()
    return frozenset(int(x) for x in cell.split("|") if x.strip()) if cell else frozenset()


def _flag01(cell: str, what: str) -> bool:
    cell = cell.strip()
    if cell not in ("0", "1"):
        raise ValueError(f"{what} must be '0' or '1', got {cell!r}")
    return cell == "1"


@dataclass(frozen=True, slots=True)
class AliasRow:
    """别名表一行（``data/alias/*.csv``，列序 ``ALIAS_COLUMNS``，SYNTHESIS §3.1）。

    head           类头文本（空白规范化，大小写敏感）
    ether          ``"0"|"1"|"*"``  是否带醚标记（``Features.has_ether``）；``*`` 不限
    zero           ``"0"|"1"|"2"|"3"|"*"``  ``0:0`` 占位个数（``Features.n_placeholder``）；``*`` 不限
    fa             ``"0"|"1"|"*"``  是否带 ``(FA c:d)``（``Features.has_fa_paren``）
    when           命名谓词表达式（``classes.PREDICATES``，如 ``eos2O``、``suffix=='-d'``）；空 = 无
    class_full     目标类（``status=="reject"`` 行可空）
    status         ``verified | inferred | unverified | reject``
    reject_reason  ``status=="reject"`` 时的原因码（``reasons.REASONS`` 的键）
    inject         类带入的隐含链，``|`` 分隔（如 ``chain:27:1``）
    dialects       适用方言，``|`` 分隔
    evidence       证据（参考表行数 / 来源 / 示例）
    note           备注

    选择规则（A5 实现）：按 head 取候选 -> ``matches_shape`` 过滤 -> 谓词 -> ``specificity`` 最大者胜 -> 仍多个即表 bug。
    """
    head: str
    ether: str = "*"
    zero: str = "*"
    fa: str = "*"
    when: str = ""
    class_full: str = ""
    status: str = "verified"
    reject_reason: str = ""
    inject: str = ""
    dialects: str = ""
    evidence: str = ""
    note: str = ""

    def __post_init__(self) -> None:
        if not self.head:
            raise ValueError("AliasRow.head must be non-empty")
        if self.ether not in ("0", "1", "*"):
            raise ValueError(f"AliasRow.ether must be 0|1|*, got {self.ether!r} (head={self.head!r})")
        if self.zero not in ("0", "1", "2", "3", "*"):
            raise ValueError(f"AliasRow.zero must be 0|1|2|3|*, got {self.zero!r} (head={self.head!r})")
        if self.fa not in ("0", "1", "*"):
            raise ValueError(f"AliasRow.fa must be 0|1|*, got {self.fa!r} (head={self.head!r})")
        if self.status not in ALIAS_STATUSES:
            raise ValueError(f"AliasRow.status must be one of {ALIAS_STATUSES}, got {self.status!r}")
        if self.status == "reject":
            if self.reject_reason not in REASONS:
                raise ValueError(f"reject row {self.head!r} needs a valid reject_reason, got {self.reject_reason!r}")
        else:
            if not self.class_full:
                raise ValueError(f"AliasRow {self.head!r}: class_full is required unless status == 'reject'")
            if self.reject_reason:
                raise ValueError(f"AliasRow {self.head!r}: reject_reason only allowed with status == 'reject'")

    @property
    def specificity(self) -> int:
        """非通配条件的个数（``ether/zero/fa`` 各 1，``when`` 非空 1）；最具体者胜。"""
        return (self.ether != "*") + (self.zero != "*") + (self.fa != "*") + bool(self.when)

    def matches_shape(self, has_ether: bool, n_zero: int, has_fa: bool) -> bool:
        """``ether/zero/fa`` 三列是否匹配（不含 ``when``）。"""
        if self.ether != "*" and (self.ether == "1") != has_ether:
            return False
        if self.zero != "*" and int(self.zero) != n_zero:
            return False
        if self.fa != "*" and (self.fa == "1") != has_fa:
            return False
        return True

    @classmethod
    def from_row(cls, row: Mapping[str, str]) -> AliasRow:
        """CSV 行（``csv.DictReader`` 的 dict）-> ``AliasRow``；空单元格取默认值；未知列抛 ``ValueError``。"""
        extra = set(row) - set(ALIAS_COLUMNS)
        if extra:
            raise ValueError(f"unknown alias columns {sorted(extra)}")
        g = lambda k, d="": (row.get(k) or d).strip()  # noqa: E731
        return cls(head=g("head"), ether=g("ether", "*"), zero=g("zero", "*"), fa=g("fa", "*"), when=g("when"),
                   class_full=g("class_full"), status=g("status", "verified"), reject_reason=g("reject_reason"),
                   inject=g("inject"), dialects=g("dialects"), evidence=g("evidence"), note=g("note"))

    def to_row(self) -> dict[str, str]:
        return {k: getattr(self, k) for k in ALIAS_COLUMNS}


@dataclass(frozen=True, slots=True)
class HeadInfo:
    """M4a 输出：头部命中的别名行集合（尚未定类）。``rows`` 非空（空即 ``resolve_head`` 抛 ``class_unknown``）。

    ``resolve_head`` 约定：先按整个 head 查表；查不到且 head 带 ``Features.class_suffix`` 同形的后缀时，
    再按去后缀的 head 查（``rows[i].head`` 告诉 ``select_class`` 命中的是哪一种键）。
    """
    head: str
    rows: tuple[AliasRow, ...] = ()


@dataclass(frozen=True, slots=True)
class ClassInfo:
    """类信息表一行（``data/class_info/<fam>.csv``，每个 class_full 一行，共 305；列序 ``CLASS_INFO_COLUMNS``）。

    class_full    ``nucleus`` 或 ``nucleus-mod``（nucleus/mod 显式给出，不靠字符串切分；``As-PL``/``CDP-DAG`` 含连字符）
    nucleus, mod  ``mod`` 可为 ``""``
    family        ``gpl | gl | sl | fa``（X12：类名全局一张表，家族只由这一列决定）
    nchain        名义链数集合（CSV ``"2"`` 或 ``"2|3"``）；**空集 = 豁免**（``ST``、``WD``，不检查链数）
    ether         该类是醚类（``class_full`` 含 ``ether`` 段）
    lyso          该类是溶血类
    implied_O     类隐含的氧个数（FAHFA/NAGly/SE/CE/ST-Gly…）；``implied_mode == "none"`` ⇔ 0
    implied_mode  ``none | add | max``：``add`` = 名字氧 + implied_O；``max`` = max(名字氧, implied_O)
    slash_role    ``none | base_first | sterol_first | hfa_second``（S14，``/`` 的角色）
    oxy_policy    ``typed_only | sphingo_OH | class_OH``（§4.3）
    pnnl_O        PNNL 类后缀 ``-d`` 对该类贡献的氧数（S4；``-t`` 的氧数归家族钩子）；无则 None
    notes         备注
    """
    class_full: str
    nucleus: str
    mod: str
    family: str
    nchain: frozenset[int]
    ether: bool = False
    lyso: bool = False
    implied_O: int = 0
    implied_mode: str = "none"
    slash_role: str = "none"
    oxy_policy: str = "typed_only"
    pnnl_O: int | None = None
    notes: str = ""

    def __post_init__(self) -> None:
        expect = f"{self.nucleus}-{self.mod}" if self.mod else self.nucleus
        if not self.nucleus or self.class_full != expect:
            raise ValueError(f"ClassInfo.class_full {self.class_full!r} != nucleus[-mod] ({expect!r})")
        if self.family not in FAMILIES:
            raise ValueError(f"ClassInfo.family must be one of {FAMILIES}, got {self.family!r} ({self.class_full})")
        if any(n < 1 for n in self.nchain):
            raise ValueError(f"ClassInfo.nchain must hold ints >= 1 ({self.class_full})")
        if self.implied_mode not in IMPLIED_MODES:
            raise ValueError(f"ClassInfo.implied_mode must be one of {IMPLIED_MODES} ({self.class_full})")
        if (self.implied_mode == "none") != (self.implied_O == 0) or self.implied_O < 0:
            raise ValueError(f"ClassInfo: implied_mode={self.implied_mode!r} inconsistent with implied_O={self.implied_O} "
                             f"({self.class_full})")
        if self.slash_role not in SLASH_ROLES:
            raise ValueError(f"ClassInfo.slash_role must be one of {SLASH_ROLES} ({self.class_full})")
        if self.oxy_policy not in OXY_POLICIES:
            raise ValueError(f"ClassInfo.oxy_policy must be one of {OXY_POLICIES} ({self.class_full})")

    @classmethod
    def from_row(cls, row: Mapping[str, str]) -> ClassInfo:
        """CSV 行 -> ``ClassInfo``。``ether/lyso`` 必须是 ``0/1``；空 ``nchain`` = 豁免；空 ``pnnl_O`` = None。"""
        extra = set(row) - set(CLASS_INFO_COLUMNS)
        if extra:
            raise ValueError(f"unknown class_info columns {sorted(extra)}")
        g = lambda k, d="": (row.get(k) or d).strip()  # noqa: E731
        cf = g("class_full")
        pnnl = g("pnnl_O")
        return cls(class_full=cf, nucleus=g("nucleus"), mod=g("mod"), family=g("family"),
                   nchain=_parse_set(g("nchain")), ether=_flag01(g("ether"), f"{cf}.ether"),
                   lyso=_flag01(g("lyso"), f"{cf}.lyso"), implied_O=int(g("implied_O", "0")),
                   implied_mode=g("implied_mode", "none"), slash_role=g("slash_role", "none"),
                   oxy_policy=g("oxy_policy", "typed_only"), pnnl_O=int(pnnl) if pnnl else None,
                   notes=g("notes"))

    def to_row(self) -> dict[str, str]:
        return {
            "class_full": self.class_full, "nucleus": self.nucleus, "mod": self.mod, "family": self.family,
            "nchain": "|".join(str(n) for n in sorted(self.nchain)),
            "ether": "1" if self.ether else "0", "lyso": "1" if self.lyso else "0",
            "implied_O": str(self.implied_O), "implied_mode": self.implied_mode, "slash_role": self.slash_role,
            "oxy_policy": self.oxy_policy, "pnnl_O": "" if self.pnnl_O is None else str(self.pnnl_O),
            "notes": self.notes,
        }

    def to_class_res(self, alias_status: str = "verified", inject: Iterable[str] = ()) -> ClassRes:
        """用本行 + 命中的别名行状态/inject 生成 ``ClassRes``（字段一一对应，避免各处手抄）。"""
        return ClassRes(
            class_full=self.class_full, nucleus=self.nucleus, mod=self.mod, family=self.family,
            nchain=self.nchain, implied_O=self.implied_O, implied_mode=self.implied_mode,
            oxy_policy=self.oxy_policy, slash_role=self.slash_role, alias_status=alias_status,
            inject=tuple(inject), ether=self.ether, lyso=self.lyso, pnnl_O=self.pnnl_O)


@dataclass(frozen=True, slots=True)
class ClassRes:
    """M4b 输出：选定的类（``ClassInfo`` + 命中别名行的 ``status`` / ``inject``）。

    前 11 个字段与 SYNTHESIS §2.1 一致；``ether / lyso / pnnl_O`` 是 W0 补的可选字段（取自 ``ClassInfo``），
    供语义层判断 ``ether_not_allowed``、``source_hint_lyso``、PNNL 后缀氧。

    nchain       名义链数集合；**空集 = 不检查**（``ST``、``WD``）
    alias_status ``verified | inferred | unverified``（reject 行不会走到这里，它在 M4 抛 ``Reject``）
    inject       类带入的隐含链（如 ``("chain:27:1",)``）
    """
    class_full: str
    nucleus: str
    mod: str
    family: str
    nchain: frozenset[int]
    implied_O: int
    implied_mode: str
    oxy_policy: str
    slash_role: str
    alias_status: str
    inject: tuple[str, ...] = ()
    ether: bool = False
    lyso: bool = False
    pnnl_O: int | None = None


# ==========================================================================
# M6 / M7：标签层
# ==========================================================================
@dataclass(frozen=True, slots=True)
class Chain:
    """标签级的一条链（label_rules §5.1）。**DB 是标签口径**（``ether == "P"`` 时已含 +1）。

    toks   已分型 token ``((name, n), ...)``，name ∈ ``TOKEN_ORDER`` 去掉 ``P``，``n >= 1``；顺序不限（M7 负责排序）
    ether  ``""`` | ``"O"``（输出 ``O-x:y``，无 ``;P``）| ``"P"``（输出 ``O-x:(y+1);P``，DB 已 +1）
    role   链角色：``""`` 或家族自定义（已知：``base`` 鞘氨醇碱基、``hfa`` 羟基脂肪酸链、``sterol`` 甾醇链、``acyl`` N-酰基链）；
           仅供诊断/元数据，M7 排序**不用**角色
    oh_read  这条链的 ``OH`` token 里有几个是**输入里显式写出的羟基**（``;OH`` ``;OH2`` ``;(2OH)`` 以及 ``d/t/m`` 前缀换算），
           其余（``;O`` / ``;On`` 按鞘脂政策补成 ``OH``、PNNL 后缀 / 缺氧数补的碱基氧、类隐含氧、家族钩子加的 ``OH``）都是包假设出来的。
           ``None`` = 不知道来源（直接构造的 ``Chain``）。**只有 ``lm.py`` 读它**（鞘脂碱基按羟基数推断，见 ``lm._base_by_hydroxyl``）；
           不参与相等 / 哈希 / 假设合并，M7 与 L1-L3 不看它
    """
    C: int
    DB: int
    toks: tuple[tuple[str, int], ...] = ()
    ether: str = ""
    role: str = ""
    oh_read: int | None = field(default=None, compare=False)

    def __post_init__(self) -> None:
        if self.ether not in ETHER_KINDS:
            raise ValueError(f"Chain.ether must be one of {ETHER_KINDS}, got {self.ether!r}")
        for name, n in self.toks:
            if name not in _CHAIN_TOKS or n < 1:
                raise ValueError(f"bad Chain token {(name, n)!r}; allowed names {sorted(_CHAIN_TOKS)}")


@dataclass(frozen=True, slots=True)
class Struct:
    """一个完整假设（结构）。

    class_full  该假设下的类（通常等于 ``ClassRes.class_full``）
    chains      标签级链列表；**None = 只有总量**（level 1，不拆链）
    tail_o      未分型氧个数（不在任何链上）
    total       ``(C, DB, O)`` 总量；``chains is None`` 时必填，否则由链推出可省略（None）
    """
    class_full: str
    chains: tuple[Chain, ...] | None
    tail_o: int = 0
    total: tuple[int, int, int] | None = None

    def __post_init__(self) -> None:
        if self.chains is None and self.total is None:
            raise ValueError("Struct with chains=None needs total=(C, DB, O)")
        if self.total is not None and len(self.total) != 3:
            raise ValueError("Struct.total must be (C, DB, O)")
        if self.tail_o < 0:
            raise ValueError("Struct.tail_o must be >= 0")


@dataclass(frozen=True, slots=True)
class Interp:
    """M6 输出。``alts`` 是 ≥1、≤6 个完整假设；某一层只有所有假设给出同一字符串才输出（D1，§4.1-8）。

    flags          语义层打的 flags（词表 ``FLAGS`` 或已注册的家族 flags）
    ether_linkage  ``"P" | "O" | "unspecified" | ""``（label_rules Q4）
    """
    cls: ClassRes
    alts: tuple[Struct, ...]
    flags: tuple[str, ...] = ()
    ether_linkage: str = ""

    def __post_init__(self) -> None:
        if not 1 <= len(self.alts) <= 6:
            raise ValueError(f"Interp.alts must hold 1..6 hypotheses, got {len(self.alts)}")
        if self.ether_linkage not in ETHER_LINKAGES:
            raise ValueError(f"Interp.ether_linkage must be one of {ETHER_LINKAGES}")


@dataclass(frozen=True, slots=True)
class Labels:
    """M7 输出（benchmark 去 sn 写法）。层级不足为 ``""``。

    level        ``3 if L3 else 2 if L2 else 1``（构造时核对；``dataclasses.replace`` 降级时必须同时改 ``level``）
    chains_meta  与 L2/L3 同序的 ``[{"C","DB","toks":{"OH":2},"ether":""|"O"|"P"}]``（进 ``meta.chains``）
    sum_meta     ``{"C","DB","O"}``（进 ``meta.sum``）
    """
    L1: str
    L2: str
    L3: str
    level: int
    chains_meta: list = field(default_factory=list)
    sum_meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.L1:
            raise ValueError("Labels.L1 must be non-empty")
        if self.L3 and not self.L2:
            raise ValueError("Labels: L3 without L2")
        want = 3 if self.L3 else 2 if self.L2 else 1
        if self.level != want:
            raise ValueError(f"Labels.level {self.level} != {want}")
