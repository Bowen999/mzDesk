"""M4 类别名表：加载、校验、头部解析（M4a）与类选择（M4b）。SYNTHESIS §3，A5。

对外接口
--------
``resolve_head(head: str) -> HeadInfo``                 M4a。按头部文本查别名表，返回命中的候选行（尚未定类）
``select_class(hi: HeadInfo, lx: Lexed, cfg: Config) -> ClassRes``   M4b。用 ``Lexed.feats`` 在候选行里选唯一行
``class_info(class_full: str) -> ClassInfo``            一个 class_full 的类信息（nucleus/mod/家族/链数/隐含氧…）
``all_class_info() -> Mapping[str, ClassInfo]``         305 个类信息（只读）
``PREDICATES``                                          命名谓词注册表（13 个，见下）；``eval_when`` / ``parse_when`` 执行 ``when`` 表达式
``class_sibling(class_full, *, ether=None, lyso=None) -> str | None``   同核、仅醚/溶血不同的兄弟类（如 PC -> PC-lyso、PE-nacyl -> PE-ether-nacyl）
``head_ether_marker(head) -> str``                      类头部自带的醚种类：``"O"``（Plasmanyl-/-O）、``"P"``（Plasmenyl-/-P/plasmalogen）、``""``
``alias_rows() -> tuple[AliasRow, ...]``                加载后的全部别名行（诊断/测试用）
``parse_glycan_head(head) -> str | None``               LIPID MAPS 糖链组成式 ``Hex(n)-HexNAc(m)-…-Cer`` -> 规范组成 key
``check_tables() -> list[str]``                         重跑表自检（唯一性穷举）；加载时已跑过，冲突即加载错误
``load_class_info / load_alias_rows / load_glycan / build_tables``   纯加载函数（可指向任意目录，测试用）
``ClassTableError``                                     表数据有误（加载错误）

数据文件（``lipid_name_map/data/``，按家族分文件，加载时拼接）
--------------------------------------------------------------
* ``class_info/<fam>.csv``  每个 class_full 一行（列序 ``types.CLASS_INFO_COLUMNS``）；文件名 = 家族；305 类须全覆盖、不得重复
* ``alias/*.csv``           别名行（列序 ``types.ALIAS_COLUMNS``）；重复键 ``(head, ether, zero, fa, when)`` = 加载错误；
                            每行 ``evidence`` 必须非空；``inject`` 只能是 ``chain:C:DB`` 列表
* ``glycan.csv``            糖链组成 key -> class_full（``kind,key,class_full,status,evidence,note``）；class_full 必须在词表内

选择规则（SYNTHESIS §3.1）
--------------------------
1. ``resolve_head``：整个 head（空白折叠、大小写敏感）先查；查不到且 head 以 ``-O -P -d -t -m`` 结尾，再查去后缀的 head
   （此时只保留显式接受该后缀的行，即 ``when`` 里的 ``suffix`` 项对该后缀成立，例如 ``suffix=='-O'``；没有这样的行 ->
   ``class_unknown``，所以 ``Cer-O``、``PC-d`` 不会被悄悄当成 ``Cer`` / ``PC``）；再不行，试 LM 糖链组成式（``glycan.csv``）；
   都没有 -> ``class_unknown``。
   整个 head 的行全是同一原因的 reject 行 -> 直接抛那个原因码（``class_not_in_vocab`` / ``polymer`` / …）。
2. ``select_class``：候选行按 ``ether/zero/fa`` 三列（对 ``Features.has_ether / n_placeholder / has_fa_paren``）与 ``when``
   过滤，**最具体者胜**（``AliasRow.specificity``：非 ``*`` 的列数 + ``when`` 非空）。最具体的行若指向不同结果（类/状态/原因/inject）
   -> ``class_chain_conflict``（表有 bug；加载期的穷举检查保证不会发生）；没有行匹配 -> ``class_chain_conflict``。
   选中 reject 行 -> 抛其 ``reject_reason``；选中 unverified 行且 ``cfg.allow_unverified_alias`` 为假 -> ``class_alias_unverified``。
3. 加载期穷举：对每个 head，把它的行引用到的特征轴（醚、0:0 个数、(FA)、谓词读的字段）取遍值域，逐个状态选行，
   最具体者指向不同结果即加载错误（``Features`` 的约束不在此处检查，只看轴取值）。分两种命中方式各穷举一遍：
   直接命中（``class_suffix`` 由键本身的结尾决定）与去后缀命中（``class_suffix`` 取遍 ``-O -P -d -t -m``，只有接受该后缀的行参与）。
4. 家族 agent 往 ``alias/<fam>.csv`` 追加行即可（文件名随意，加载器 glob 整个目录；键不得与已有行重复）。``native.csv`` /
   ``tools.csv`` / ``ref_main.csv`` 自 A8 起手工维护（D21）：属于本家族 ``class_full`` 的行允许小范围 Edit 修改，不要整文件脚本重写；
   ``tools/build_alias_table.py`` 已冻结（``--init-class-info`` 在 ``class_info`` 已存在时拒绝执行，``--write`` 需显式强制）。

``when`` 表达式（``&`` 连接的项；项 = ``name`` | ``!name`` | ``name==值`` | ``name!=值``；字符串值用单引号）
-----------------------------------------------------------------------------------------------------------
=================  =====  ==================================================================
谓词               类型   读取的 ``Features`` 字段
=================  =====  ==================================================================
``has_ether``      bool   ``has_ether``（任意醚标记）
``zero``           int    ``n_placeholder``（``0:0`` 个数）
``nested_fa``      bool   ``has_fa_paren``
``pre_acyl``       bool   ``has_pre_acyl``
``eos2O``          bool   ``eos2O``（LipidBlast ``;2O/;2O``；纯词法，只有 ``Cer`` / ``HexCer`` 头的行才读它，D18）
``tail_Gly``       bool   ``tail_conj in {Gly, G, G+S}``（名字带甘氨酸共轭；双共轭 ``G+S`` 也算）
``tail_Tau``       bool   ``tail_conj in {Tau, T, T+S}``（名字带牛磺酸共轭）
``tail_S``         bool   ``tail_conj in {S, G+S, T+S}``（名字带硫酸酯）
``tail_Hex``       bool   ``tail_conj == Hex``
``tail_conj``      str    ``tail_conj`` 原文（``tail_conj=='Phe'``，氨基酸共轭 ``ST …;Phe``，0.4.0；取值见 ``types.TAIL_CONJ_VALUES``）
``chains``         int    ``n_chains``（不含占位链）
``suffix``         str    ``class_suffix``（``'' -O -P -d -t -m``）
``sphingoid_prefix``  str ``sphingoid_prefix``（``'' d t m``）
=================  =====  ==================================================================

双共轭（D17）
-------------
``ST …;G;S`` / ``ST …;T;S``（ST-Gly-sulf / ST-Tau-sulf，参考表各 3 行）：``tail_conj`` 取 ``'G+S'`` / ``'T+S'``，三个共轭谓词
按"名字里有没有这种共轭"取值，所以 ``ST`` 的各共轭行要写成互斥（``tail_Gly & !tail_S`` 与 ``tail_Gly & tail_S`` …），否则加载期穷举报歧义。

已知限制
--------
* 类头部里的醚种类标记（``Plasmanyl-`` ``Plasmenyl-`` ``plasmalogen`` 后缀 ``-O/-P``）不进 ``ClassRes``；语义层从 ``Lexed.head`` 读。
"""
from __future__ import annotations

import csv
import re
import threading
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, field
from functools import lru_cache
from itertools import product
from pathlib import Path
from types import MappingProxyType, SimpleNamespace
from typing import Any, Final

from .config import Config
from .contract import load_vocab
from .types import (ALIAS_COLUMNS, BASE_VALUES, CLASS_INFO_COLUMNS, CLASS_SUFFIXES, FAMILIES, TAIL_CONJ_VALUES,
                    AliasRow, ClassInfo, ClassRes, Features, HeadInfo, Lexed, Reject)

__all__ = [
    "resolve_head", "select_class", "class_info", "all_class_info", "alias_rows", "PREDICATES", "Predicate", "Term",
    "parse_when", "eval_when", "parse_glycan_head", "normalize_head", "check_tables", "find_dead_rows", "class_sibling",
    "head_ether_marker",
    "load_class_info", "load_alias_rows", "load_glycan", "build_tables", "tables", "Tables", "GlycanRow",
    "ClassTableError", "DATA_DIR", "GLYCAN_ORDER",
]

DATA_DIR: Final = Path(__file__).resolve().parent / "data"
_SUFFIXES: Final = tuple(s for s in CLASS_SUFFIXES if s)
GLYCAN_ORDER: Final = ("Hex", "HexNAc", "Fuc", "NeuAc", "NeuGc", "KDN")
GLYCAN_COLUMNS: Final = ("kind", "key", "class_full", "status", "evidence", "note")
GLYCAN_KINDS: Final = ("comp",)


class ClassTableError(ValueError):
    """类表数据有误（格式、重复键、未知类、无证据、冲突…）。加载期抛出，等同于 bug，不是用户输入问题。"""


# ==========================================================================
# 命名谓词
# ==========================================================================
@dataclass(frozen=True, slots=True)
class Predicate:
    """一个命名谓词：``pred(feats)`` 返回值（bool/int/str）；``axis`` = 它读取的 ``Features`` 字段（穷举检查按轴取值）。"""
    name: str
    kind: str                      # 'bool' | 'int' | 'str'
    axis: str
    fn: Callable[[Any], Any]
    values: tuple[str, ...] = ()   # str 型的合法字面值
    doc: str = ""

    def __call__(self, feats: Any) -> Any:
        return self.fn(feats)


def _mk(name: str, kind: str, axis: str, fn: Callable[[Any], Any], values: tuple[str, ...] = (), doc: str = "") -> Predicate:
    return Predicate(name, kind, axis, fn, values, doc)


PREDICATES: Final[Mapping[str, Predicate]] = MappingProxyType({p.name: p for p in (
    _mk("has_ether", "bool", "has_ether", lambda f: bool(f.has_ether), doc="任意醚标记（前缀或 e/p 后缀）"),
    _mk("zero", "int", "n_placeholder", lambda f: f.n_placeholder, doc="0:0 占位链个数"),
    _mk("nested_fa", "bool", "has_fa_paren", lambda f: bool(f.has_fa_paren), doc="有嵌套 (FA c:d)"),
    _mk("pre_acyl", "bool", "has_pre_acyl", lambda f: bool(f.has_pre_acyl), doc="有括号 O-酰基前缀 (O-c:d)"),
    _mk("eos2O", "bool", "eos2O", lambda f: bool(f.eos2O), doc="LipidBlast EOS 判据：第 2 条主链有数字在前的 ;2O"),
    _mk("tail_Gly", "bool", "tail_conj", lambda f: f.tail_conj in ("Gly", "G", "G+S"), doc="甘氨酸共轭 ;Gly / ;G（含 ;G;S）"),
    _mk("tail_Tau", "bool", "tail_conj", lambda f: f.tail_conj in ("Tau", "T", "T+S"), doc="牛磺酸共轭 ;Tau / ;T（含 ;T;S）"),
    _mk("tail_S", "bool", "tail_conj", lambda f: f.tail_conj in ("S", "G+S", "T+S"), doc="硫酸酯 ;S（含 ;G;S / ;T;S）"),
    _mk("tail_Hex", "bool", "tail_conj", lambda f: f.tail_conj == "Hex", doc="己糖苷 ;Hex"),
    _mk("tail_conj", "str", "tail_conj", lambda f: f.tail_conj, values=("",) + TAIL_CONJ_VALUES,
        doc="共轭名原文（``tail_conj=='Phe'``：氨基酸共轭 ``ST …;Phe``；单共轭原样，``Gly``/``G`` 是两种写法）"),
    _mk("chains", "int", "n_chains", lambda f: f.n_chains, doc="主链数（不含占位链）"),
    _mk("suffix", "str", "class_suffix", lambda f: f.class_suffix, values=CLASS_SUFFIXES, doc="类后缀 -O -P -d -t -m"),
    _mk("sphingoid_prefix", "str", "sphingoid_prefix", lambda f: f.sphingoid_prefix, values=BASE_VALUES,
        doc="第一条带碱基前缀链的 d/t/m"),
)})

#: 穷举检查里每个特征轴的取值域（int 轴会再并上 ``when`` 里引用到的值）
_AXIS_DOMAIN: Final = MappingProxyType({
    "has_ether": (False, True),
    "n_placeholder": (0, 1, 2, 3),
    "has_fa_paren": (False, True),
    "has_pre_acyl": (False, True),
    "eos2O": (False, True),
    "tail_conj": ("",) + TAIL_CONJ_VALUES,
    "n_chains": (0, 1, 2, 3, 4, 5),
    "class_suffix": CLASS_SUFFIXES,
    "sphingoid_prefix": BASE_VALUES,
})


@dataclass(frozen=True, slots=True)
class Term:
    """``when`` 表达式里的一项。"""
    name: str
    neg: bool = False
    op: str = ""               # '' | '==' | '!='
    value: Any = None

    def holds(self, feats: Any) -> bool:
        v = PREDICATES[self.name](feats)
        r = (v == self.value) if self.op == "==" else (v != self.value) if self.op == "!=" else bool(v)
        return (not r) if self.neg else r


_TERM_RE: Final = re.compile(r"^(!?)\s*([A-Za-z_][A-Za-z0-9_]*)\s*(?:(==|!=)\s*(.+?))?\s*$")


@lru_cache(maxsize=4096)
def parse_when(expr: str) -> tuple[Term, ...]:
    """解析 ``when``；空串 -> ``()``（恒真）。非法表达式抛 ``ClassTableError``。"""
    expr = expr.strip()
    if not expr:
        return ()
    terms: list[Term] = []
    for part in expr.split("&"):
        m = _TERM_RE.match(part.strip())
        if not m:
            raise ClassTableError(f"bad when term {part.strip()!r} in {expr!r}")
        bang, name, op, raw = m.groups()
        pred = PREDICATES.get(name)
        if pred is None:
            raise ClassTableError(f"unknown predicate {name!r} in {expr!r}; known: {sorted(PREDICATES)}")
        if pred.kind == "bool":
            if op:
                raise ClassTableError(f"boolean predicate {name!r} takes no comparison in {expr!r}")
            terms.append(Term(name, bool(bang)))
            continue
        if not op:
            raise ClassTableError(f"predicate {name!r} needs a comparison (== or !=) in {expr!r}")
        if bang:
            raise ClassTableError(f"use != instead of ! with {name!r} in {expr!r}")
        if pred.kind == "int":
            if not re.fullmatch(r"\d+", raw):
                raise ClassTableError(f"predicate {name!r} compares with an integer, got {raw!r} in {expr!r}")
            value: Any = int(raw)
        else:
            if len(raw) >= 2 and raw[0] == raw[-1] and raw[0] in "'\"":
                raw = raw[1:-1]
            if raw not in pred.values:
                raise ClassTableError(f"predicate {name!r}: value {raw!r} not in {pred.values} ({expr!r})")
            value = raw
        terms.append(Term(name, False, op, value))
    return tuple(terms)


def eval_when(expr: str, feats: Any) -> bool:
    """``when`` 表达式在 ``feats``（``Features`` 或同名属性对象）上是否成立。"""
    return all(t.holds(feats) for t in parse_when(expr))


def _mentions_suffix(row: AliasRow) -> bool:
    return any(t.name == "suffix" for t in parse_when(row.when))


def _accepts_suffix(row: AliasRow, suffix: str) -> bool:
    """靠去后缀命中时，这一行是否显式接受该后缀（``when`` 里有 ``suffix`` 项，且所有 ``suffix`` 项对该后缀成立）。"""
    ns = SimpleNamespace(class_suffix=suffix)
    terms = [t for t in parse_when(row.when) if t.name == "suffix"]
    return bool(terms) and all(t.holds(ns) for t in terms)


# ==========================================================================
# 加载
# ==========================================================================
@dataclass(frozen=True, slots=True)
class GlycanRow:
    kind: str
    key: str
    class_full: str
    status: str
    evidence: str
    note: str = ""


@dataclass(frozen=True, slots=True)
class Tables:
    info: Mapping[str, ClassInfo]
    rows: tuple[AliasRow, ...]
    index: Mapping[str, tuple[AliasRow, ...]]
    glycan: Mapping[str, GlycanRow]
    sources: Mapping[tuple[str, str, str, str, str], str] = field(default_factory=dict)


def _read_csv(path: Path, columns: tuple[str, ...]) -> list[tuple[int, dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as fh:      # utf-8-sig：容忍编辑器加的 BOM
        rdr = csv.DictReader(fh)
        if tuple(rdr.fieldnames or ()) != tuple(columns):
            raise ClassTableError(f"{path.name}: header {rdr.fieldnames} != expected {list(columns)}")
        out = []
        for r in rdr:
            line = rdr.line_num
            if None in r or any(v is None for v in r.values()):
                raise ClassTableError(f"{path.name}:{line}: wrong number of cells")
            out.append((line, r))
    return out


def normalize_head(head: str) -> str:
    """空白折叠 + strip（大小写敏感，其余原样）。"""
    return " ".join(str(head).split())


def _vocab_sets(vocab: Mapping[str, Any]) -> tuple[frozenset[str], frozenset[str], frozenset[str]]:
    return frozenset(vocab["class_full"]), frozenset(vocab["class_nucleus"]), frozenset(vocab["class_mod"])


def load_class_info(directory: Path | None = None, *, vocab: Mapping[str, Any] | None = None,
                    require_complete: bool = True) -> dict[str, ClassInfo]:
    """读 ``class_info/*.csv``。文件名须是家族名（行的 ``family`` 与之一致）；重复、未知类、nucleus/mod 不在词表 -> 错误。"""
    directory = Path(directory) if directory is not None else DATA_DIR / "class_info"
    vocab = vocab if vocab is not None else load_vocab()
    classes, nuclei, mods = _vocab_sets(vocab)
    info: dict[str, ClassInfo] = {}
    where: dict[str, str] = {}
    files = sorted(directory.glob("*.csv"))
    if not files:
        raise ClassTableError(f"no class_info CSV in {directory}")
    for path in files:
        for line, r in _read_csv(path, CLASS_INFO_COLUMNS):
            loc = f"{path.name}:{line}"
            try:
                ci = ClassInfo.from_row(r)
            except (ValueError, TypeError) as e:
                raise ClassTableError(f"{loc}: {e}") from None
            if path.stem in FAMILIES and ci.family != path.stem:
                raise ClassTableError(f"{loc}: {ci.class_full} has family {ci.family!r} but lives in {path.name}")
            if ci.class_full not in classes:
                raise ClassTableError(f"{loc}: {ci.class_full!r} is not in the v0.2.11 vocabulary")
            if ci.nucleus not in nuclei or (ci.mod and ci.mod not in mods):
                raise ClassTableError(f"{loc}: nucleus/mod of {ci.class_full} not in the vocabulary lists")
            if ci.class_full in info:
                raise ClassTableError(f"{loc}: duplicate class {ci.class_full!r} (first at {where[ci.class_full]})")
            info[ci.class_full] = ci
            where[ci.class_full] = loc
    if require_complete and set(info) != set(classes):
        miss, extra = sorted(classes - set(info)), sorted(set(info) - classes)
        raise ClassTableError(f"class_info must cover the vocabulary exactly: missing {miss[:8]}, extra {extra[:8]}")
    return info


_INJECT_RE: Final = re.compile(r"^chain:\d+:\d+(\|chain:\d+:\d+)*$")


def load_alias_rows(directory: Path | None = None, *, classes: Iterable[str] | None = None
                    ) -> tuple[list[AliasRow], dict[tuple[str, str, str, str, str], str]]:
    """读 ``alias/*.csv``（按文件名排序）。返回 (行, 键->位置)。重复键、非法 when、无证据、class 不在词表 -> 错误。"""
    directory = Path(directory) if directory is not None else DATA_DIR / "alias"
    classes = frozenset(classes) if classes is not None else _vocab_sets(load_vocab())[0]
    rows: list[AliasRow] = []
    where: dict[tuple[str, str, str, str, str], str] = {}
    files = sorted(directory.glob("*.csv"))
    if not files:
        raise ClassTableError(f"no alias CSV in {directory}")
    for path in files:
        for line, r in _read_csv(path, ALIAS_COLUMNS):
            loc = f"{path.name}:{line}"
            try:
                row = AliasRow.from_row(r)
            except (ValueError, TypeError) as e:
                raise ClassTableError(f"{loc}: {e}") from None
            if row.head != normalize_head(row.head):
                raise ClassTableError(f"{loc}: head {row.head!r} is not whitespace-normalised")
            parse_when(row.when)
            if row.status != "reject" and row.class_full not in classes:
                raise ClassTableError(f"{loc}: class_full {row.class_full!r} is not in the v0.2.11 vocabulary")
            if row.status == "reject" and row.class_full:
                raise ClassTableError(f"{loc}: reject row {row.head!r} must leave class_full empty")
            if not row.evidence:
                raise ClassTableError(f"{loc}: evidence is required ({row.head!r})")
            if row.inject and not _INJECT_RE.match(row.inject):
                raise ClassTableError(f"{loc}: inject {row.inject!r} must look like chain:C:DB[|chain:C:DB]")
            key = (row.head, row.ether, row.zero, row.fa, row.when)
            if key in where:
                raise ClassTableError(f"{loc}: duplicate alias key {key} (first at {where[key]})")
            where[key] = loc
            rows.append(row)
    return rows, where


_GLY_KEY_RE: Final = re.compile(r"^Hex(?:\d+)?(?:-(?:HexNAc|Fuc|NeuAc|NeuGc|KDN)(?:\d+)?)*$")
_LM_HEAD_RE: Final = re.compile(r"^Hex(?:\(\d+\))?(?:-(?:HexNAc|Fuc|NeuAc|NeuGc|KDN)(?:\(\d+\))?)*-Cer$")
_GLY_TOKEN_RE: Final = re.compile(r"(HexNAc|Hex|Fuc|NeuAc|NeuGc|KDN)(?:\((\d+)\))?")


def parse_glycan_head(head: str) -> str | None:
    """LM 糖链组成式 ``Hex(2)-HexNAc-Fuc-NeuAc(2)-Cer`` -> 规范 key ``Hex2-HexNAc-Fuc-NeuAc2``（顺序 Hex,HexNAc,Fuc,NeuAc,NeuGc,KDN，
    个数 1 省略）。不是组成式、有重复糖或个数为 0 返回 ``None``。"""
    if not _LM_HEAD_RE.match(head):
        return None
    counts: dict[str, int] = {}
    for tok in head[:-4].split("-"):
        m = _GLY_TOKEN_RE.fullmatch(tok)
        if m is None:
            return None
        name, n = m.group(1), int(m.group(2) or 1)
        if n < 1 or name in counts:
            return None
        counts[name] = n
    return "-".join(k + (str(counts[k]) if counts[k] > 1 else "") for k in GLYCAN_ORDER if k in counts)


def _canonical_glycan_key(key: str) -> bool:
    """规范组成 key：糖按 GLYCAN_ORDER 严格递增、个数 1 省略、个数 >= 2 不带前导 0，且以 Hex 开头。"""
    if not _GLY_KEY_RE.match(key):
        return False
    last = -1
    for tok in key.split("-"):
        m = re.fullmatch(r"(HexNAc|Hex|Fuc|NeuAc|NeuGc|KDN)(\d*)", tok)
        if m is None:
            return False
        idx = GLYCAN_ORDER.index(m.group(1))
        n = m.group(2)
        if idx <= last or (n and (n.startswith("0") or n == "1")):
            return False
        last = idx
    return True


def load_glycan(path: Path | None = None, *, classes: Iterable[str] | None = None) -> dict[str, GlycanRow]:
    """读 ``glycan.csv``：组成 key -> 行。class_full 必须在词表内且等于 ``Cer-<key>``；key 须是规范写法且唯一。"""
    path = Path(path) if path is not None else DATA_DIR / "glycan.csv"
    classes = frozenset(classes) if classes is not None else _vocab_sets(load_vocab())[0]
    out: dict[str, GlycanRow] = {}
    for line, r in _read_csv(path, GLYCAN_COLUMNS):
        loc = f"{path.name}:{line}"
        g = GlycanRow(r["kind"].strip(), r["key"].strip(), r["class_full"].strip(), r["status"].strip(),
                      r["evidence"].strip(), r["note"].strip())
        if g.kind not in GLYCAN_KINDS:
            raise ClassTableError(f"{loc}: unknown glycan kind {g.kind!r}; known {GLYCAN_KINDS}")
        if g.status not in ("verified", "inferred", "unverified"):
            raise ClassTableError(f"{loc}: glycan status must be verified|inferred|unverified, got {g.status!r}")
        if not g.evidence:
            raise ClassTableError(f"{loc}: evidence is required ({g.key!r})")
        if g.class_full not in classes:
            raise ClassTableError(f"{loc}: class_full {g.class_full!r} is not in the v0.2.11 vocabulary")
        if not _canonical_glycan_key(g.key):
            raise ClassTableError(f"{loc}: composition key {g.key!r} is not canonical "
                                  f"(order {GLYCAN_ORDER}, count 1 omitted)")
        if g.class_full != "Cer-" + g.key:
            raise ClassTableError(f"{loc}: class_full {g.class_full!r} != 'Cer-' + key {g.key!r}")
        if g.key in out:
            raise ClassTableError(f"{loc}: duplicate composition key {g.key!r}")
        out[g.key] = g
    return out


# ==========================================================================
# 穷举唯一性检查
# ==========================================================================
def _outcome(r: AliasRow) -> tuple[str, str, str, str]:
    return (r.class_full, r.status, r.reject_reason, r.inject)


def _axes_of(rows: list[AliasRow]) -> dict[str, tuple]:
    axes: dict[str, set] = {}
    for r in rows:
        if r.ether != "*":
            axes.setdefault("has_ether", set())
        if r.zero != "*":
            axes.setdefault("n_placeholder", set())
        if r.fa != "*":
            axes.setdefault("has_fa_paren", set())
        for t in parse_when(r.when):
            pred = PREDICATES[t.name]
            axes.setdefault(pred.axis, set())
            if pred.kind == "int" and t.op:
                axes[pred.axis].add(t.value)
    out = {}
    for ax, extra in axes.items():
        dom = list(_AXIS_DOMAIN[ax])
        out[ax] = tuple(dom + sorted(v for v in extra if v not in dom))
    return out


def _head_suffix(head: str) -> str:
    """词法层对整个 head 取的类后缀（与 ``Features.class_suffix`` 同口径）。"""
    for suf in _SUFFIXES:
        if head.endswith(suf) and len(head) > len(suf):
            return suf
    return ""


def _states(rows: list[AliasRow], head: str):
    """穷举一个 head 的行所引用的特征轴；产出 (state_dict, via_suffix, 最具体的匹配行)。

    直接命中（via_suffix=False）时 ``class_suffix`` 由键本身决定；靠去后缀命中（via_suffix=True）时只有写了
    ``suffix`` 的行参与，且 ``class_suffix`` 取遍非空值。
    """
    axes = _axes_of(rows)
    axes.pop("class_suffix", None)
    names = sorted(axes)
    modes = [(False, _head_suffix(head))]
    if any(_mentions_suffix(r) for r in rows):
        modes += [(True, s) for s in _SUFFIXES]
    for combo in product(*(axes[n] for n in names)):
        base = dict(zip(names, combo))
        he, nz, fa = base.get("has_ether", False), base.get("n_placeholder", 0), base.get("has_fa_paren", False)
        for via, suf in modes:
            state = dict(base, class_suffix=suf)
            ns = SimpleNamespace(**state)
            cands = [r for r in rows
                     if (not via or _accepts_suffix(r, suf)) and r.matches_shape(he, nz, fa)
                     and all(t.holds(ns) for t in parse_when(r.when))]
            if not cands:
                yield state, via, []
                continue
            top = max(r.specificity for r in cands)
            yield state, via, [r for r in cands if r.specificity == top]


def check_alias_uniqueness(rows: Iterable[AliasRow], *, max_per_head: int = 3) -> list[str]:
    """穷举每个 head 的特征状态，最具体的行指向不同结果即报告。返回问题清单（空 = 唯一）。"""
    by_head: dict[str, list[AliasRow]] = {}
    for r in rows:
        by_head.setdefault(r.head, []).append(r)
    problems: list[str] = []
    for head, rs in by_head.items():
        n = 0
        for state, via, win in _states(rs, head):
            if len({_outcome(r) for r in win}) > 1:
                n += 1
                if n <= max_per_head:
                    desc = "; ".join(f"{r.class_full or 'REJECT:' + r.reject_reason}[when={r.when!r}]" for r in win)
                    problems.append(f"head {head!r}: ambiguous at {state} via_suffix={via}: {desc}")
    return problems


def find_dead_rows(rows: Iterable[AliasRow]) -> list[AliasRow]:
    """在任何特征状态下都不会胜出的行（被更具体的行完全遮蔽，或条件不可满足）。"""
    by_head: dict[str, list[AliasRow]] = {}
    for r in rows:
        by_head.setdefault(r.head, []).append(r)
    dead: list[AliasRow] = []
    for head, rs in by_head.items():
        alive: set[int] = set()
        for _state, _via, win in _states(rs, head):
            alive.update(id(r) for r in win)
        dead.extend(r for r in rs if id(r) not in alive)
    return dead


def build_tables(data_dir: Path | None = None, *, vocab: Mapping[str, Any] | None = None,
                 require_complete: bool = True) -> Tables:
    """读全部类表并做自检；任何问题抛 ``ClassTableError``。"""
    data_dir = Path(data_dir) if data_dir is not None else DATA_DIR
    vocab = vocab if vocab is not None else load_vocab()
    classes = _vocab_sets(vocab)[0]
    info = load_class_info(data_dir / "class_info", vocab=vocab, require_complete=require_complete)
    rows, where = load_alias_rows(data_dir / "alias", classes=classes)
    glycan = load_glycan(data_dir / "glycan.csv", classes=classes) if (data_dir / "glycan.csv").exists() else {}
    problems = check_alias_uniqueness(rows)
    if problems:
        raise ClassTableError("alias table is ambiguous:\n  " + "\n  ".join(problems[:12]))
    index: dict[str, list[AliasRow]] = {}
    for r in rows:
        index.setdefault(r.head, []).append(r)
    return Tables(info=MappingProxyType(info), rows=tuple(rows),
                  index=MappingProxyType({h: tuple(v) for h, v in index.items()}),
                  glycan=MappingProxyType(glycan), sources=MappingProxyType(where))


_LOCK: Final = threading.Lock()
_TABLES: Tables | None = None


def tables(*, reload: bool = False) -> Tables:
    """进程内缓存的默认类表（首次调用时加载并自检）。"""
    global _TABLES
    if _TABLES is None or reload:
        with _LOCK:
            if _TABLES is None or reload:
                _TABLES = build_tables()
    return _TABLES


def check_tables() -> list[str]:
    """重跑默认类表的唯一性穷举；返回问题清单（空 = 通过）。"""
    return check_alias_uniqueness(tables().rows)


# ==========================================================================
# 公共接口
# ==========================================================================
def class_info(class_full: str) -> ClassInfo:
    """一个 class_full 的类信息；不在词表 -> ``KeyError``。"""
    try:
        return tables().info[class_full]
    except KeyError:
        raise KeyError(f"unknown class_full {class_full!r}") from None


def all_class_info() -> Mapping[str, ClassInfo]:
    return tables().info


def alias_rows() -> tuple[AliasRow, ...]:
    return tables().rows


def class_sibling(class_full: str, *, ether: bool | None = None, lyso: bool | None = None) -> str | None:
    """同核、``mod`` 里仅"醚 / 溶血"不同的兄弟类；没有则 ``None``。``ether`` / ``lyso`` 为 ``None`` 表示保持不变。

    按 ``mod`` 的部件集合比较（``ether``、``lyso``；``lyso2`` 是另一个层级，不当作 ``lyso``）：
    ``PC -> PC-lyso``、``PC-ether -> PC-ether-lyso``、``PE-nacyl -> PE-ether-nacyl``、``Cer-Hex -> None``。
    """
    info = tables().info
    ci = info[class_full]
    parts = set(ci.mod.split("-")) if ci.mod else set()
    for flag, name in ((ether, "ether"), (lyso, "lyso")):
        if flag is True:
            parts.add(name)
        elif flag is False:
            parts.discard(name)
    hits = [cf for cf, o in info.items()
            if o.nucleus == ci.nucleus and (set(o.mod.split("-")) if o.mod else set()) == parts]
    return hits[0] if len(hits) == 1 else None


#: 后缀 -O 指 O-酰基链而不是醚的头部（MS-DIAL 规范化名，S16）
_O_ACYL_SUFFIX_HEADS: Final = frozenset({"AHexCer-O"})


#: 位置号 + 单字母取代基（``-1-P`` 1-磷酸、``-3-O``）：这种 ``-P`` / ``-O`` 不是醚类后缀
_LOCANT_SUFFIX: Final = re.compile(r"-\d+-[OP]$")


def head_ether_marker(head: str) -> str:
    """类头部文本自带的醚种类标记：``"O"``（Plasmanyl- / plasmanyl- / 后缀 -O）、``"P"``（Plasmenyl- / plasmenyl- / 后缀 -P /
    含 plasmalogen 的旧名）、``""``（无标记；``Ether*`` 只说明是醚、种类未定，也返回 ``""``；``AHexCer-O`` 的 -O 是 O-酰基，
    也返回 ``""``）。纯字符串规则，不查表；是不是真的醚类由 ``ClassRes.ether`` 与 ``select_class`` 决定。
    位置号之后的 ``-P`` / ``-O``（``Sph-1-P`` 的 1-磷酸）不是醚后缀，返回 ``""``。"""
    h = normalize_head(head)
    low = h.lower()
    if h in _O_ACYL_SUFFIX_HEADS:
        return ""
    locant = _LOCANT_SUFFIX.search(h) is not None                  # ``Sph-1-P`` 的 -P 是 1-磷酸（位置号之后），不是缩醛磷脂后缀
    if low.startswith("plasmenyl-") or "plasmalogen" in low or (h.endswith("-P") and len(h) > 2 and not locant):
        return "P"
    if low.startswith("plasmanyl-") or (h.endswith("-O") and len(h) > 2 and not locant):
        return "O"
    return ""


def resolve_head(head: str) -> HeadInfo:
    """M4a：头部文本 -> 命中的别名行（见模块文档"选择规则"1）。

    抛出：``class_unknown``（表里没有）；头部的行全是同一原因的 reject 行时抛该原因（``class_not_in_vocab`` / ``polymer`` /
    ``common_name`` / ``unsupported_subclass``）；LM 糖链组成式不在 glycan.csv -> ``glycan_unknown``。
    ``HeadInfo.head`` 是规范化后的整个头部；若靠去后缀命中，``rows[i].head`` 是去后缀的键。
    """
    t = tables()
    h = normalize_head(head)
    if not h:
        raise Reject("class_unknown", "empty head")
    rows = t.index.get(h)
    if rows is None:
        for suf in _SUFFIXES:
            if h.endswith(suf) and len(h) > len(suf):
                base = t.index.get(h[:-len(suf)])
                if base:
                    rows = tuple(r for r in base if _accepts_suffix(r, suf))
                    if not rows:
                        raise Reject("class_unknown", f"{h}: class suffix {suf} is not accepted for {h[:-len(suf)]}")
                    break
    if rows is None and _LM_HEAD_RE.match(h):
        key = parse_glycan_head(h)
        g = t.glycan.get(key) if key else None
        if g is None:
            raise Reject("glycan_unknown", f"{h}: composition {key!r} is not in the vocabulary")
        rows = (AliasRow(head=h, class_full=g.class_full, status=g.status, dialects="lipidmaps", evidence=g.evidence,
                         note="glycan composition rule"),)
    if rows is None:
        raise Reject("class_unknown", h)
    if all(r.status == "reject" for r in rows):
        reasons = {r.reject_reason for r in rows}
        if len(reasons) == 1:
            raise Reject(next(iter(reasons)), h)
    return HeadInfo(head=h, rows=rows)


def select_class(hi: HeadInfo, lx: Lexed, cfg: Config) -> ClassRes:
    """M4b：在 ``hi.rows`` 里按 ``lx.feats`` 选唯一行并给出 ``ClassRes``（见模块文档"选择规则"2）。"""
    f: Features = lx.feats
    he, nz, fa = f.has_ether, f.n_placeholder, f.has_fa_paren
    cands = [r for r in hi.rows if r.matches_shape(he, nz, fa) and all(t.holds(f) for t in parse_when(r.when))]
    if not cands:
        raise Reject("class_chain_conflict",
                     f"{hi.head}: no class row fits ether={he}, placeholders={nz}, (FA)={fa}, suffix={f.class_suffix!r}")
    top = max(r.specificity for r in cands)
    win = [r for r in cands if r.specificity == top]
    if len({_outcome(r) for r in win}) > 1:
        raise Reject("class_chain_conflict", f"{hi.head}: ambiguous alias rows {[r.class_full or r.reject_reason for r in win]}")
    r = win[0]
    if r.status == "reject":
        raise Reject(r.reject_reason, hi.head)
    if r.status == "unverified" and not cfg.allow_unverified_alias:
        raise Reject("class_alias_unverified", f"{hi.head} -> {r.class_full}")
    return class_info(r.class_full).to_class_res(r.status, r.inject.split("|") if r.inject else ())
