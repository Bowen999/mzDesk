"""M7a：去 sn 写法的链排序、链渲染、``to_nosn``（SYNTHESIS §2.2 M7 / §4.1，依据 ``label_rules`` §2、§4）。

对外接口
--------
``sort_key(chain, mecpr=False) -> tuple``
    一条标签级链的排序键 ``(C, DB, token种类数, 首token首字母, 含前缀的链全串)``；**全降序**使用。
    ``chain.ether == "P"`` 时 ``P`` 记作一个 token（算入种类数、可作首字母），``O-`` 前缀只出现在最后一项
    （同 (C, DB)、同 token 时 ``O-18:1`` > ``18:1``，故醚链排在前；醚链不置首）。
    ``mecpr=False``（默认）= benchmark 口径：``Me``/``cPr`` 不参与键，也不出现在链串里。
``render_chain(chain, mecpr=False) -> str``
    L3 链串：``[O-]C:DB[;tok…][;P]``。token 顺序 ``OH oxo Ep COOH OOH OO O [Me cPr] P``，计数 n>=2 紧跟数字
    （``OH2``，永不 ``OH1``）；``ether`` 非空 -> ``O-`` 前缀，``"P"`` -> 末尾 ``;P``。``Chain.DB`` 已是标签口径，
    本函数**不再** +1。
``order_chains(chains, mecpr=False) -> list[Chain]``
    按 ``sort_key`` 全降序的稳定排序（L2、L3 共用同一次排序）。
``to_nosn(label) -> str``
    sn 版（参考表写法）L2/L3 字符串 -> 去 sn 写法：链按同一排序键全降序、分隔符统一 ``_``、括号占位（``;(O2)``）
    不参与排序并按原顺序挂回末尾。空串/None 原样返回；单链恒等。无法解析 -> ``ValueError``。
    与 ``lipid_plus_v2/.../spectra_pipeline/nosn_labels.nosn_label`` 在参考表全部 793,906 个唯一 L2/L3 组合上逐字符一致。
``parse_label(label) -> (class_full, [ParsedChain], placeholders)``  ``to_nosn`` 的解析层
``is_ether_class(class_full) -> bool``  ``class_full`` 含 ``ether`` 段（``(?:^|-)ether(?:-|$)``）
``chain_valence(chain) -> int``         一条链上已分型 token 的氧价之和（COOH/OOH/OO=2，OH/oxo/Ep/O=1，Me/cPr=0）
``tok_fmt(name, n) / tok_split(text)``   token 的文本形式 <-> (名, 个数)

只依赖 ``types``（常量 ``TOKEN_ORDER`` ``TOKEN_VALENCE`` 与 ``Chain``）；不 import 任何其他包内模块。
"""
from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass

from .types import TOKEN_ORDER, TOKEN_VALENCE, Chain

__all__ = [
    "sort_key", "render_chain", "order_chains", "chain_tokens", "chain_valence", "to_nosn", "parse_label",
    "ParsedChain", "is_ether_class", "tok_fmt", "tok_split",
]

#: 链上可出现的 token（``P`` 由 ``Chain.ether == "P"`` 表达，不在 ``Chain.toks`` 里）
_L3_TOKENS: tuple[str, ...] = tuple(t for t in TOKEN_ORDER if t != "P")
#: benchmark 口径默认丢弃的支链/环丙烷 token
_MECPR: frozenset[str] = frozenset({"Me", "cPr"})

_ETHER_CLASS_RE = re.compile(r"(?:^|-)ether(?:-|$)")
_TOK_RE = re.compile(r"^([A-Za-z]+?)(\d*)$")
# 一条链：可选前缀 + C:DB + 若干 ``;token`` 或 ``;(括号占位)``。token 须以字母开头。
_CHAIN_RE = re.compile(r"^(O-|P-)?(\d+):(\d+)((?:;(?:[A-Za-z][A-Za-z0-9]*|\([^()]*\)))*)$")
_MOD_RE = re.compile(r";(\([^()]*\)|[A-Za-z][A-Za-z0-9]*)")
_SEP_RE = re.compile(r"[/_]")


# ==========================================================================
# token 小工具
# ==========================================================================
def tok_fmt(name: str, n: int) -> str:
    """``("OH", 1) -> "OH"``，``("OH", 2) -> "OH2"``（n=1 不写数字，永不 ``OH1``）。"""
    if n < 1:
        raise ValueError(f"token count must be >= 1, got {n} for {name!r}")
    return name if n == 1 else f"{name}{n}"


def tok_split(text: str) -> tuple[str, int]:
    """``"OH2" -> ("OH", 2)``，``"oxo" -> ("oxo", 1)``。无法解析抛 ``ValueError``。"""
    m = _TOK_RE.match(text)
    if not m:
        raise ValueError(f"bad token {text!r}")
    return m.group(1), int(m.group(2) or 1)


def is_ether_class(class_full: str) -> bool:
    """``class_full`` 含 ``ether`` 段（``PE-ether``、``DG-Hex-ether``、``PE-ether-nacyl-lyso``…）。"""
    return bool(_ETHER_CLASS_RE.search(class_full))


# ==========================================================================
# 标签级链：渲染与排序键
# ==========================================================================
def _merged_counts(chain: Chain) -> dict[str, int]:
    counts: dict[str, int] = {}
    for name, n in chain.toks:
        counts[name] = counts.get(name, 0) + n
    return counts


def chain_tokens(chain: Chain, mecpr: bool = False) -> list[str]:
    """链上 token 的文本列表（已按 ``TOKEN_ORDER`` 排好、同族合并、``ether == "P"`` 时末尾追加 ``"P"``）。

    ``mecpr=False`` 时 ``Me``/``cPr`` 被丢弃。
    """
    counts = _merged_counts(chain)
    out = [tok_fmt(name, counts[name]) for name in _L3_TOKENS
           if counts.get(name, 0) and (mecpr or name not in _MECPR)]
    if chain.ether == "P":
        out.append("P")
    return out


def _render(chain: Chain, toks: list[str]) -> str:
    return ("O-" if chain.ether else "") + f"{chain.C}:{chain.DB}" + "".join(";" + t for t in toks)


def render_chain(chain: Chain, mecpr: bool = False) -> str:
    """L3 链串，如 ``O-16:0``、``O-18:1;P``、``18:1;OH2``、``20:4;OH;oxo2;COOH``。"""
    return _render(chain, chain_tokens(chain, mecpr))


def sort_key(chain: Chain, mecpr: bool = False) -> tuple[int, int, int, str, str]:
    """排序键 ``(C, DB, token种类数, 首token首字母, 含前缀的链全串)``，**全降序**使用（label_rules §2.1）。

    * ``len`` 是 token **种类**数（``OH2`` 算 1 种；``;P`` 算 1 种）；
    * 首字母区分大小写，降序下 ``o``(oxo) > ``c``(cPr) > ``P`` > ``O`` > ``M``(Me) > ``E``(Ep) > ``C``(COOH)；
    * 链全串含 ``O-`` 前缀：同 (C, DB)、同 token 时醚链在前；前三项按数值比（``9:0`` 在 ``18:0`` 之后）。
    """
    toks = chain_tokens(chain, mecpr)
    return (chain.C, chain.DB, len(toks), toks[0][0] if toks else "", _render(chain, toks))


def order_chains(chains: Iterable[Chain], mecpr: bool = False) -> list[Chain]:
    """按 ``sort_key`` 全降序排列（稳定；键完全相同的链渲染相同，次序无所谓）。"""
    keyed = [(sort_key(c, mecpr), c) for c in chains]
    keyed.sort(key=lambda kc: kc[0], reverse=True)
    return [c for _, c in keyed]


def chain_valence(chain: Chain) -> int:
    """链上已分型 token 的氧价之和（不含 ``tail_o``）。"""
    return sum(TOKEN_VALENCE[name] * n for name, n in chain.toks)


# ==========================================================================
# to_nosn：sn 版 L2/L3 字符串 -> 去 sn 字符串
# ==========================================================================
@dataclass(frozen=True, slots=True)
class ParsedChain:
    """``parse_label`` 的一条链：``prefix`` ``""|"O-"|"P-"``；``mods`` 是不含括号占位的 ``;`` 记号（原样文本）。"""
    prefix: str
    C: int
    DB: int
    mods: tuple[str, ...] = ()

    def render(self) -> str:
        return f"{self.prefix}{self.C}:{self.DB}" + "".join(";" + m for m in self.mods)

    def key(self) -> tuple[int, int, int, str, str]:
        """与 ``sort_key`` 同构的键（官方 ``nosn_label`` 的键）。"""
        return (self.C, self.DB, len(self.mods), self.mods[0][0] if self.mods else "", self.render())


def parse_label(label: str) -> tuple[str, list[ParsedChain], tuple[str, ...]]:
    """``"<cf> 链<sep>链…[;(O..)]"`` -> ``(class_full, [ParsedChain], 括号占位元组)``。

    ``sep`` 为 ``/`` 或 ``_``；括号占位（``(O2)``）可出现在任一链后，按出现顺序收集（不含前导 ``;``）。
    无法解析抛 ``ValueError``。
    """
    if not isinstance(label, str):
        raise ValueError(f"label must be str, got {type(label).__name__}")
    cf, sp, body = label.partition(" ")
    if not sp or not cf or not body:
        raise ValueError(f"label has no '<class> <chains>' shape: {label!r}")
    chains: list[ParsedChain] = []
    placeholders: list[str] = []
    for part in _SEP_RE.split(body):
        m = _CHAIN_RE.match(part)
        if not m:
            raise ValueError(f"chain token does not match: {part!r} in {label!r}")
        mods: list[str] = []
        for mod in _MOD_RE.findall(m.group(4)):
            (placeholders if mod.startswith("(") else mods).append(mod)
        chains.append(ParsedChain(m.group(1) or "", int(m.group(2)), int(m.group(3)), tuple(mods)))
    return cf, chains, tuple(placeholders)


def to_nosn(label: str) -> str:
    """sn 版 L2/L3 -> 去 sn 版（见模块文档）。空串/None 原样返回；无法解析抛 ``ValueError``。"""
    if label is None or label == "":
        return label
    cf, chains, placeholders = parse_label(label)
    chains.sort(key=ParsedChain.key, reverse=True)
    return (f"{cf} " + "_".join(c.render() for c in chains)
            + "".join(f";{ph}" for ph in placeholders))
