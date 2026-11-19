"""M6c 约定层：醚种类 / P-DB+1 / 鞘氨醇碱基氧 / 原子 -> token / 价数 / 介质缩写字典 / 类带入链（纯函数）。

**唯一实现点**：``P-`` 的 DB+1（``P_DB_BONUS``、``chain_mark``、``ether_kind``）只在本文件。core 与家族钩子都不得自己 +1。

对外接口
--------
``ether_kind(lx, cr) -> (kind, flags)``
    把链标记（``O-`` ``P-`` ``e``、原生 ``O-x:y;P``）与类头标记（``Plasmanyl-`` ``Plasmenyl-`` ``plasmalogen`` ``-O`` ``-P``）
    合并成醚种类 ``kind ∈ {"", "O", "P"}``：

    * 链标记与类标记都有且不同 -> ``class_chain_conflict``；``dO-`` 前缀 / ``p`` 后缀 -> ``unsupported_syntax``（D18）
    * 有醚标记但 ``cr.ether`` 为假 -> ``ether_not_allowed``
    * ``cr.ether`` 为真而没有任何标记（``Ether*``、``PC-ether 16:0_18:1``）-> ``"O"``（链落点未知，linkage 未表态）
    * ``flags`` 含 ``P_db_plus1`` 当且仅当 ``kind == "P"`` 且 +1 确实要加（``P-`` 前缀或类标记；原生 ``;P`` 标记的 DB 已是标签口径，不加）
    * ``flags`` 含 ``ether_prefix_inferred_from_P``：类是醚类、某条链写了 ``;P`` 却没有 ``O-`` / ``e``（LIPID+ 模型解码输出 ``PE-ether 22:6_20:1;P``）
      -> 当作 ``O-20:1;P``（见 ``chain_mark`` 的 ``ether_class``）
``chain_mark(raw, *, ether_class=False) -> (kind, db_extra)``      一条链自带的醚标记；``db_extra`` 是待加的 DB（``P-`` 为 1，其余 0）
``label_db(db, db_extra) -> int``            字面 DB -> 标签口径（**唯一**的 P 加一）
``class_marker(head) -> str``                类头自带醚种类 ``"O" | "P" | ""``（``classes.head_ether_marker`` 的薄封装）
``base_oxygen(base) -> int``                 鞘氨醇碱基前缀 ``d/t/m`` -> 2/3/1（参数可为 ``RawChain`` 或字符串）
``valence(toks) -> int``                     token 氧价之和（``toks`` 是 ``{name: n}`` 或 ``((name, n), ...)``）
``atoms_to_tokens(atoms, policy, *, assume_oh=True, base="") -> (toks: dict, tail_o: int)``
    一条链上的原子 -> 已分型 token 与"未分型氧"：

    * 已分型原子（``OH oxo Ep COOH OOH OO``、``OH_pos`` 每项 1 个 OH、``O_other`` -> ``O``、``Me`` ``cPr``）-> token
    * 碱基前缀 ``base`` ∈ d/t/m -> ``OH`` += 2/3/1（结构事实，不是假设）；同一条链上再带未分型氧 -> ``class_chain_conflict``
    * ``O_count``（``;O`` ``;2O`` ``;O2``）：policy 为 ``sphingo_OH`` 且 ``assume_oh`` -> ``OH``（调用方打 ``oxygen_assumed_OH``）；
      否则若本链已有携带氧价的 token -> ``O`` token（"其他型 O"）；否则 -> 返回的 ``tail_o``（未分型，L3 空）
    * ``P`` / ``conj_*`` 原子不在这里处理
``demote_overlimit(toks, DB) -> bool``       token 是否违反标签合法性（DB < oxo+Ep+COOH，或某 token 超过词表上限）
``expand_oxylipin(abbr) -> (C, DB, {token: n})``   读 ``data/oxylipin.csv``（7 个介质缩写）；未知缩写 ``KeyError``
``parse_inject(inject) -> [(C, DB), ...]``    ``ClassRes.inject`` 的 ``chain:C:DB`` 项
``adjacent_groups(lx) -> bool``               名字是相邻括号组写法 ``CL(28:0)(28:0)``（词法把各组当 ``_`` 连接的链，``RawChain.raw`` 保留括号）
``is_harmonized_name(lx, works, separators) -> bool``   MS-DIAL 规范化名（S16）：``AHexCer-O`` 头，或 ``_`` 分隔、无碱基前缀 / 类后缀、氧
                                              （``;O`` / ``;On``）贴在链上。core 对**所有** ``sphingo_OH`` 策略的类用它（含家族是 ``fa`` 的 ``SL``）
``L3_TOKEN_MAX``                             L3 词表各 token 的个数上限（fam_fa 4.3）
``FLAGS_input``                              本模块打的 flags（经 ``register_family_flags("input", …)`` 登记，没有改 ``types.FLAGS``）：
                                              ``ether_prefix_inferred_from_P``、``aldehyde_as_oxo``
``lipidmatch_ox_db_extra(atoms) -> int``             LipidMatch 括号形式 ``(Ke) (CHO) (COOH)`` 的 C=O 个数 = 字面 DB 之外要加的 DB（标签口径的 DB 含 C=O）
"""
from __future__ import annotations

import csv
from collections.abc import Iterable, Mapping
from functools import lru_cache
from pathlib import Path
from typing import Final

from ..classes import head_ether_marker
from ..contract import register_family_flags
from ..types import TOKEN_VALENCE, Atom, ClassRes, Lexed, RawChain, Reject

__all__ = [
    "P_DB_BONUS", "BASE_OXYGEN", "L3_TOKEN_MAX", "ether_kind", "chain_mark", "label_db", "class_marker", "base_oxygen",
    "valence", "atoms_to_tokens", "demote_overlimit", "expand_oxylipin", "oxylipin_table", "parse_inject",
    "adjacent_groups", "is_harmonized_name", "FLAGS_input", "bare_p_chain", "lipidmatch_ox_db_extra",
]

#: 输入端（2026-10-31 LM_INPUT_AUDIT 修复）新增的 flags；同 ``lm.FLAGS_lm``，经 ``register_family_flags`` 登记而不改 ``types.FLAGS``（那个词表有测试钉住）
FLAGS_input: frozenset[str] = frozenset({"ether_prefix_inferred_from_P", "aldehyde_as_oxo"})
register_family_flags("input", FLAGS_input)

#: 烯基醚（``P-``）的 DB 加一量；整个包只有这一个常数
P_DB_BONUS: Final = 1
#: 鞘氨醇碱基前缀 -> 碱基上的羟基个数（C1-O 含在内；S1）
BASE_OXYGEN: Final = {"d": 2, "t": 3, "m": 1}
#: L3 词表各 token 的个数上限（超限的链不出 L3）
L3_TOKEN_MAX: Final = {"OH": 9, "oxo": 6, "Ep": 2, "COOH": 4, "OOH": 2, "OO": 2, "O": 5, "Me": 9, "cPr": 3}

_DATA: Final = Path(__file__).resolve().parents[1] / "data"
_TYPED: Final = frozenset({"OH", "oxo", "Ep", "COOH", "OOH", "OO"})
_ZERO_VAL: Final = frozenset({"Me", "cPr"})


# ==========================================================================
# 价数 / 碱基
# ==========================================================================
def valence(toks: Mapping[str, int] | Iterable[tuple[str, int]]) -> int:
    """token 氧价之和（``Me``/``cPr`` 为 0）。"""
    items = toks.items() if isinstance(toks, Mapping) else toks
    return sum(TOKEN_VALENCE[k] * n for k, n in items)


def base_oxygen(base: str | RawChain) -> int:
    """``d``->2 ``t``->3 ``m``->1，其余 0。"""
    b = base.base if isinstance(base, RawChain) else base
    return BASE_OXYGEN.get(b, 0)


# ==========================================================================
# 醚种类
# ==========================================================================
def class_marker(head: str) -> str:
    """类头自带的醚种类标记 ``"O"`` / ``"P"`` / ``""``（A5 的 ``head_ether_marker``）。"""
    return head_ether_marker(head)


def label_db(db: int, db_extra: int) -> int:
    """字面 DB -> 标签口径。``db_extra`` 来自 ``chain_mark``（P- 前缀）或假设枚举里的 P 落点（``P_DB_BONUS``）。"""
    return db + db_extra


def bare_p_chain(rc: RawChain) -> bool:
    """链上有 ``;P`` 原子，却没有 ``O-`` / ``P-`` / ``dO-`` 前缀，也没有 ``e`` / ``p`` 后缀（LIPID+ 解码输出 ``20:1;P``）。"""
    return not rc.pre and not rc.suf and any(a.kind == "P" for a in rc.atoms)


def chain_mark(rc: RawChain, *, ether_class: bool = False) -> tuple[str, int]:
    """一条链自带的醚标记 -> ``(kind, db_extra)``。

    ``O-`` / ``e`` -> ``("O", 0)``；``P-`` -> ``("P", 1)``；原生缩醛 ``O-x:y;P`` -> ``("P", 0)``（DB 已是标签口径）；无标记 -> ``("", 0)``。
    ``dO-`` 前缀与 ``p`` 后缀 -> ``unsupported_syntax``（D18）；``;P`` 与 ``P-`` 并存 -> ``class_chain_conflict``；
    ``;P`` 落在 DB 为 0 的链上（``O-16:0;P``，结构上不存在）-> ``malformed``（对所有家族通用，原来只有 GL 钩子查）。

    ``;P`` 贴在没有 ``O-`` 的链上（LIPID+ 模型的解码输出 ``PE-ether 22:6_20:1;P``）：``;P`` 只能出现在烯基醚链上，语义上等于 ``O-20:1;P``。
    调用方说明这是醚类（``ether_class``）时就当 ``O-x:y;P`` 收（``("P", 0)``，flag ``ether_prefix_inferred_from_P`` 由 ``ether_kind`` 打）；
    不是醚类、或链上还有 ``N-`` / ``d t m`` 前缀（不可能是醚链）-> 照旧 ``class_chain_conflict``，不放宽别的冲突。
    """
    if rc.pre == "dO-":
        raise Reject("unsupported_syntax", f"dialkyl-ether prefix dO- in {rc.raw!r}")
    if rc.suf == "p":
        raise Reject("unsupported_syntax", f"p (plasmalogen) chain suffix in {rc.raw!r}")
    has_p_atom = False
    for a in rc.atoms:
        if a.kind == "P":
            has_p_atom = True
            break
    if has_p_atom:
        if rc.pre != "O-" and rc.suf != "e":
            if not (ether_class and bare_p_chain(rc) and not rc.n_acyl and not rc.base):
                raise Reject("class_chain_conflict", f";P acetal marker needs an O- ether chain: {rc.raw!r}")
        if rc.DB < 1:                                      # 原生 ;P 的 DB 是标签口径，已含乙烯醚那个双键（至少 1）：O-16:0;P 不存在
            raise Reject("malformed", f";P acetal marker on a chain without a double bond (its label DB includes the vinyl-ether "
                                      f"one, so DB >= 1): {rc.raw!r}")
        return "P", 0
    if rc.pre == "P-":
        return "P", P_DB_BONUS
    if rc.pre == "O-" or rc.suf == "e":
        return "O", 0
    return "", 0


def ether_kind(lx: Lexed, cr: ClassRes) -> tuple[str, tuple[str, ...]]:
    """链标记 + 类头标记 -> ``(kind, flags)``，见模块文档。"""
    mark = ""
    extra = 0
    inferred = False
    for rc in lx.chains:
        if rc.placeholder:
            continue
        k, e = chain_mark(rc, ether_class=cr.ether)
        if k:
            if mark:
                raise Reject("multi_ether", "more than one chain carries an ether mark")
            mark, extra = k, e
            inferred = bare_p_chain(rc)
    ck = class_marker(lx.head)
    if mark and ck and mark != ck:
        raise Reject("class_chain_conflict", f"class head {lx.head!r} says {ck!r}-ether but the chain says {mark!r}")
    kind = mark or ck
    if kind and not cr.ether:
        raise Reject("ether_not_allowed", f"{cr.class_full} does not allow an ether chain (mark {kind!r})")
    if not kind and cr.ether:
        kind = "O"
    flags: tuple[str, ...] = ()
    if kind == "P" and (not mark or extra):
        flags = ("P_db_plus1",)
    if inferred:
        flags += ("ether_prefix_inferred_from_P",)
    return kind, flags


# ==========================================================================
# 原子 -> token
# ==========================================================================
def atoms_to_tokens(atoms: Iterable[Atom], policy: str, *, assume_oh: bool = True,
                    base: str = "") -> tuple[dict[str, int], int]:
    """一条链上的原子 -> ``(toks, tail_o)``，规则见模块文档。"""
    toks: dict[str, int] = {}
    untyped = 0
    for a in atoms:
        k = a.kind
        if k == "O_count":
            untyped += a.n
        elif k in _TYPED or k in _ZERO_VAL:
            toks[k] = toks.get(k, 0) + a.n
        elif k == "OH_pos":
            toks["OH"] = toks.get("OH", 0) + a.n
        elif k == "O_other":
            toks["O"] = toks.get("O", 0) + a.n
        # "P" 与 "conj_*" 由调用方处理
    bo = BASE_OXYGEN.get(base, 0)
    if bo:
        if untyped:
            raise Reject("class_chain_conflict", f"sphingoid prefix {base!r} together with an oxygen count on the same chain")
        toks["OH"] = toks.get("OH", 0) + bo
    tail = 0
    if untyped:
        if policy == "sphingo_OH" and assume_oh:
            toks["OH"] = toks.get("OH", 0) + untyped
        elif any(TOKEN_VALENCE[k] > 0 for k in toks):
            toks["O"] = toks.get("O", 0) + untyped
        else:
            tail = untyped
    return toks, tail


def lipidmatch_ox_db_extra(atoms: Iterable[Atom]) -> int:
    """LipidMatch 括号修饰 ``(Ke) (CHO) (COOH)`` 的 C=O 个数（``Atom.via`` 非空的原子）：字面 DB 不含它，标签口径的 DB 含（每个 +1）。"""
    return sum(a.n for a in atoms if a.via)


def demote_overlimit(toks: Mapping[str, int], db: int) -> bool:
    """True = 这条链的 token 不合法（``DB < oxo+Ep+COOH`` 或某 token 超词表上限），调用方应把它的氧退回未分型。"""
    need = toks.get("oxo", 0) + toks.get("Ep", 0) + toks.get("COOH", 0)
    if need and db < need:
        return True
    for k, n in toks.items():
        lim = L3_TOKEN_MAX.get(k)
        if lim is not None and n > lim:
            return True
    return False


# ==========================================================================
# 介质缩写 / 类带入链
# ==========================================================================
@lru_cache(maxsize=1)
def oxylipin_table() -> Mapping[str, tuple[int, int, str, int]]:
    """``data/oxylipin.csv`` -> ``{abbr: (C, DB, token, n)}``（7 行）。"""
    out: dict[str, tuple[int, int, str, int]] = {}
    with (_DATA / "oxylipin.csv").open(newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            out[r["abbr"].strip()] = (int(r["C"]), int(r["DB"]), r["token"].strip(), int(r["n"]))
    return out


def expand_oxylipin(abbr: str) -> tuple[int, int, dict[str, int]]:
    """介质缩写 -> ``(C, DB, {token: n})``（DB 已是标签口径，环氧环已含在内）。"""
    c, db, tok, n = oxylipin_table()[abbr]
    return c, db, {tok: n}


def parse_inject(inject: Iterable[str]) -> list[tuple[int, int]]:
    """``("chain:27:1",)`` -> ``[(27, 1)]``；非 ``chain:`` 项忽略。"""
    out: list[tuple[int, int]] = []
    for item in inject:
        parts = item.split(":")
        if len(parts) == 3 and parts[0] == "chain":
            out.append((int(parts[1]), int(parts[2])))
    return out


def adjacent_groups(lx: Lexed) -> bool:
    """名字是相邻括号组写法 ``(a)(b)``（LipidMatch ``CL(28:0)(28:0)``）？词法 ``lexer._adjacent_groups`` 把各组当作 ``_`` 连接的链，
    并让每条链的 ``RawChain.raw`` 保留括号 -> 全部链 ``raw`` 以 ``(`` 开头、``)`` 结尾（普通 ``PC(16:0_18:1)`` 的链 ``raw`` 不带括号）。"""
    return len(lx.chains) >= 2 and all(c.raw.startswith("(") and c.raw.endswith(")") for c in lx.chains)


#: MS-DIAL 规范化名里 ``-O`` 后缀代表 O-酰基的头（SYNTHESIS S16）
_HARMONIZED_HEADS: Final = frozenset({"AHexCer-O"})


def is_harmonized_name(lx: Lexed, works: Iterable, separators: Iterable[str]) -> bool:
    """S16：MS-DIAL 规范化名 = ``AHexCer-O`` 头，或 ``_`` 分隔、无碱基前缀、无类后缀、氧（``;O`` / ``;On``）贴在链上的鞘脂名。
    基准库规范化名有 13% 的 L3 氧换位，名字本身看不出，所以 core 对这类名字打 ``harmonized_name_oxygen_placement``。
    ``works`` 是 ``Ctx.works``（用 ``.base``），``separators`` 是 ``Ctx.separators``；判据只看词法事实与氧策略，与家族无关。"""
    if lx.head in _HARMONIZED_HEADS:
        return True
    if "_" not in separators or any(w.base for w in works) or lx.feats.class_suffix:
        return False
    return any(not rc.placeholder and any(a.kind == "O_count" for a in rc.atoms) for rc in lx.chains)
