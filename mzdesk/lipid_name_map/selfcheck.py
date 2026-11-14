"""M8：包内轻量自检与降级（SYNTHESIS §2.2 M8，依据 ``label_rules`` §5.4）。

对外接口
--------
``check_labels(lb: Labels, it: Interp) -> list[str]``
    返回违规码列表（空 = 通过；去重，顺序固定）。只检查**非空层**。码名前缀表示层：``L1_* L2_* L3_*`` 单层格式/
    与 ``it`` 的一致性，``X12_* X13_* X23_*`` 跨层，``meta_*`` 元数据，``level``。
    格式：L1 ``^<cf> C:DB(;O(n>=2)?)?$``（无 ``;O1``/``;O0``）；L2 ``C:DB(_C:DB)*`` + 至多一组尾括号 ``;(O)``/``;(On)``，
    无 token、无 ``O-``；L3 每链 ``[O-]C:DB(;tok)*``，token 顺序 ``OH oxo Ep COOH OOH OO O P``、同族不重复、n>=2 带数字
    （无 ``OH1``）、无 Me/cPr、``;P`` 至多一个且在 ``O-`` 链上、无尾括号、分隔符只有 ``_``。
    跨层：L1 的 (C,DB) = ΣL2 = ΣL3；L1 的 O = L2 括号氧 = ΣL3 氧价；L2/L3 链数与 (C,DB) 序列逐位相同；
    ``to_nosn(L2/L3)`` 幂等且 (C,DB) 非增；醚类 L3 恰一条 ``O-``、非醚类零条。
    与 ``it`` 一致：L1 的类与 (C,DB,O) = 每个假设；L2/L3 的链多重集 = 每个假设的链多重集；链数在名义链数内；
    L3 非空时每个假设 ``tail_o == 0``。
    ``chains_meta/sum_meta`` 须与已输出层一致。
``check_and_downgrade(lb: Labels, it: Interp) -> Labels``
    有违规 -> 先丢 L3，再（仍违规）丢 L2；L1 本身违规（``L1_*`` 码）或丢到只剩 L1 仍违规 -> ``Reject("internal_error")``。
    **返回值是 ``Labels``，没有 flags 槽**；要拿 ``selfcheck_downgraded_L3/L2`` 用 ``check_and_downgrade_flags``。
``check_and_downgrade_flags(lb, it) -> (Labels, tuple[str, ...])``
    同上，另返回被打的 flags（``("selfcheck_downgraded_L3",)``、``(…_L3, …_L2)`` 或 ``()``）。
``downgrade_flags(before: Labels, after: Labels) -> tuple[str, ...]``
    由降级前后的两个 ``Labels`` 推出应打的 flags（调用方只保留了两个对象时用）。
``drop_idle_flags(flags, lb: Labels) -> tuple[str, ...]``
    最终标签里没有 L3 时，去掉那些只对 L3 字符串作断言的 flag（``L3_ONLY_FLAGS``：``oxygen_assumed_OH`` ``class_implied_OH``
    ``harmonized_name_oxygen_placement``）。这三个 flag 说的是"L3 里的这个 OH / 这个氧位置是政策假设出来的"；结果只到 L1 / L2
    （和式、``_`` 分隔的 FAHFA、碱基氧放不下的 Cer 亚类……）时没有任何标签带着这个假设，flag 就是空转的（flag ⇔ 假设真的写进了 L3，
    检查器才有东西可核对）。api 在 ``check_and_downgrade_flags`` 之后调用，所以自检降级丢掉 L3 时同样生效。

``it`` 里没有 ``alts`` 能对上的情形（语义层 bug）也体现为违规码，不抛异常；本模块只在 L1 无法挽救时抛 ``internal_error``。
"""
from __future__ import annotations

import re
from collections import Counter
from collections.abc import Iterable
from dataclasses import replace
from typing import Final

from .errors import Reject
from .labels import chains_meta_l2, nominal_chain_counts
from .nosn import is_ether_class, to_nosn
from .types import TOKEN_ORDER, TOKEN_VALENCE, Interp, Labels, Struct

_L1Parsed = tuple[str, int, int, int]                              # (cf, C, DB, O)
_L2Parsed = tuple[str, list[tuple[int, int]], int]                 # (cf, [(C, DB)], 尾括号 O)
_L3Chain = tuple[bool, int, int, dict[str, int], bool]             # (O- 前缀, C, DB, {token: n}, 有 ;P)
_L3Parsed = tuple[str, list[_L3Chain]]

__all__ = ["check_labels", "check_and_downgrade", "check_and_downgrade_flags", "downgrade_flags", "drop_idle_flags", "L3_ONLY_FLAGS"]

#: 只对 L3 字符串作断言的 flag（见模块文档 ``drop_idle_flags``）
L3_ONLY_FLAGS: Final = frozenset({"oxygen_assumed_OH", "class_implied_OH", "harmonized_name_oxygen_placement"})

_L1_RE = re.compile(r"^(?P<cf>\S+) (?P<C>\d+):(?P<DB>\d+)(?:;O(?P<n>\d*))?$")
_L2_RE = re.compile(r"^(?P<cf>\S+) (?P<chains>\d+:\d+(?:_\d+:\d+)*)(?:;\(O(?P<n>\d*)\))?$")
_L3_HEAD_RE = re.compile(r"^(?P<cf>\S+) (?P<body>\S+)$")
_L3_CHAIN_RE = re.compile(r"^(?P<eth>O-)?(?P<C>\d+):(?P<DB>\d+)(?P<mods>(?:;[A-Za-z]+\d*)*)$")
_TOK_RE = re.compile(r"^([A-Za-z]+?)(\d*)$")
_TOK_INDEX = {t: i for i, t in enumerate(TOKEN_ORDER)}
_MECPR = frozenset({"Me", "cPr"})


def _count_ok(digits: str) -> bool:
    """尾缀里的氧数字：无数字 = 1；有数字则 >= 2 且无前导零。"""
    return digits == "" or (digits == str(int(digits)) and int(digits) >= 2)


def _o_from_digits(digits: str) -> int:
    return 1 if digits == "" else int(digits)


# ==========================================================================
# 各层解析（严格文法；失败返回 None 并由调用方记 ``*_format``）
# ==========================================================================
def _parse_l1(s: str) -> _L1Parsed | None:
    m = _L1_RE.match(s)
    if not m:
        return None
    has_o = ";O" in s
    n = m.group("n") or ""
    if has_o and not _count_ok(n):
        return None
    return m.group("cf"), int(m.group("C")), int(m.group("DB")), (_o_from_digits(n) if has_o else 0)


def _parse_l2(s: str) -> _L2Parsed | None:
    m = _L2_RE.match(s)
    if not m:
        return None
    n = m.group("n") or ""
    has_tail = s.endswith(")")
    if has_tail and not _count_ok(n):
        return None
    chains = [(int(c.split(":")[0]), int(c.split(":")[1])) for c in m.group("chains").split("_")]
    return m.group("cf"), chains, (_o_from_digits(n) if has_tail else 0)


def _parse_l3(s: str, errs: list[str]) -> _L3Parsed | None:
    """-> ``(cf, [(has_O_prefix, C, DB, toks{name:n}, has_P)], ...)``；文法问题记入 ``errs`` 并返回 None。"""
    m = _L3_HEAD_RE.match(s)
    if not m:
        errs.append("L3_format")
        return None
    chains: list[_L3Chain] = []
    for part in m.group("body").split("_"):
        cm = _L3_CHAIN_RE.match(part)
        if not cm:
            errs.append("L3_format")
            return None
        toks: dict[str, int] = {}
        has_p = False
        last_idx = -1
        names = [t for t in cm.group("mods").split(";") if t]
        for t in names:
            tm = _TOK_RE.match(t)
            if tm is None:                                # 不可达（_L3_CHAIN_RE 已保证），防御
                errs.append("L3_format")
                return None
            name, digits = tm.group(1), tm.group(2)
            if name not in _TOK_INDEX:
                errs.append("L3_token_unknown")
                return None
            if name in _MECPR:
                errs.append("L3_mecpr")
            if digits and not _count_ok(digits):
                errs.append("L3_token_count")
            if name == "P" and digits:
                errs.append("L3_token_count")
            idx = _TOK_INDEX[name]
            if idx <= last_idx:
                errs.append("L3_token_order")             # 乱序或同族重复
            last_idx = max(last_idx, idx)
            if name == "P":
                has_p = True
                if not cm.group("eth"):
                    errs.append("L3_P_misplaced")
            else:
                toks[name] = toks.get(name, 0) + (int(digits) if digits else 1)
        chains.append((bool(cm.group("eth")), int(cm.group("C")), int(cm.group("DB")), toks, has_p))
    if sum(1 for c in chains if c[4]) > 1:
        errs.append("L3_P_misplaced")
    return m.group("cf"), chains


def _alt_totals(st: Struct) -> tuple[int, int, int]:
    if st.chains is None:
        if st.total is None:                               # Struct.__post_init__ 已保证；防御：永不等于 L1 -> L1_vs_interp
            return (-1, -1, -1)
        c, db, o = st.total
        return c, db, o
    return (sum(c.C for c in st.chains), sum(c.DB for c in st.chains),
            sum(TOKEN_VALENCE[k] * n for c in st.chains for k, n in c.toks) + st.tail_o)


def _non_increasing(pairs: list[tuple[int, int]]) -> bool:
    return all(a >= b for a, b in zip(pairs, pairs[1:]))


def _to_nosn_safe(label: str) -> str:
    """``to_nosn`` 的容错版：无法解析 -> 返回 None（与任何合法串不等，由调用方记 order 违规）。"""
    try:
        return to_nosn(label)
    except ValueError:
        return None  # type: ignore[return-value]


# ==========================================================================
# check_labels
# ==========================================================================
def check_labels(lb: Labels, it: Interp) -> list[str]:  # noqa: C901 - 一条条列出的检查，拆开反而难读
    errs: list[str] = []
    want_level = 3 if lb.L3 else 2 if lb.L2 else 1
    if lb.level != want_level:
        errs.append("level")

    # ---------------- L1 ----------------
    p1 = _parse_l1(lb.L1)
    if p1 is None:
        errs.append("L1_format")
        cf, C1, DB1, O1 = None, None, None, None
    else:
        cf, C1, DB1, O1 = p1
        if any(st.class_full != cf for st in it.alts):
            errs.append("L1_class")
        if any(_alt_totals(st) != (C1, DB1, O1) for st in it.alts):
            errs.append("L1_vs_interp")

    # ---------------- L2 ----------------
    p2: _L2Parsed | None = None
    if lb.L2:
        p2 = _parse_l2(lb.L2)
        if p2 is None:
            errs.append("L2_format")
        else:
            cf2, ch2, O2 = p2
            seq2 = list(ch2)
            if cf is not None and cf2 != cf:
                errs.append("L2_class")
            if not _non_increasing(seq2) or _to_nosn_safe(lb.L2) != lb.L2:
                errs.append("L2_order")
            if cf is not None:
                if (sum(c for c, _ in seq2), sum(d for _, d in seq2)) != (C1, DB1):
                    errs.append("X12_CDB")
                if O2 != O1:
                    errs.append("X12_O")
            for st in it.alts:
                if st.chains is None:
                    errs.append("L2_vs_interp")
                    break
                if sorted(seq2) != sorted((c.C, c.DB) for c in st.chains):
                    errs.append("L2_vs_interp")
                    break
            for st in it.alts:
                nom = nominal_chain_counts(it, st.class_full)
                if nom and len(seq2) not in nom:
                    errs.append("L2_chain_count")
                    break

    # ---------------- L3 ----------------
    if lb.L3:
        l3errs: list[str] = []
        p3 = _parse_l3(lb.L3, l3errs)
        errs.extend(l3errs)
        if p3 is not None:
            cf3, ch3 = p3
            seq3 = [(c, d) for _, c, d, _, _ in ch3]
            if cf is not None and cf3 != cf:
                errs.append("L3_class")
            if not _non_increasing(seq3) or _to_nosn_safe(lb.L3) != lb.L3:
                errs.append("L3_order")
            if cf is not None:
                if (sum(c for c, _ in seq3), sum(d for _, d in seq3)) != (C1, DB1):
                    errs.append("X13_CDB")
                o3 = sum(TOKEN_VALENCE[k] * n for _, _, _, toks, _ in ch3 for k, n in toks.items())
                if o3 != O1:
                    errs.append("X13_O")
            if p2 is not None and [tuple(c) for c in p2[1]] != seq3:
                errs.append("X23_seq")
            n_eth = sum(1 for e, *_ in ch3 if e)
            if n_eth != (1 if is_ether_class(cf3) else 0):
                errs.append("L3_ether_count")
            for st in it.alts:
                if st.chains is None or st.tail_o != 0:
                    errs.append("L3_vs_interp")
                    break
                alt = Counter(_chain_sig(c.C, c.DB, dict(_merge(c.toks)), c.ether) for c in st.chains)
                got = Counter(_chain_sig(c, d, toks, "P" if p else ("O" if e else ""))
                              for e, c, d, toks, p in ch3)
                if alt != got:
                    errs.append("L3_vs_interp")
                    break

    # ---------------- meta ----------------
    if p1 is not None and lb.sum_meta != {"C": C1, "DB": DB1, "O": O1}:
        errs.append("meta_sum")
    if not _meta_ok(lb, p2):
        errs.append("meta_chains")

    return list(dict.fromkeys(errs))


def _merge(toks: Iterable[tuple[str, int]]) -> list[tuple[str, int]]:
    d: dict[str, int] = {}
    for k, n in toks:
        if k not in _MECPR:
            d[k] = d.get(k, 0) + n
    return sorted(d.items())


def _chain_sig(C: int, DB: int, toks: dict[str, int], ether: str) -> tuple[int, int, tuple[tuple[str, int], ...], str]:
    return (C, DB, tuple(sorted((k, n) for k, n in toks.items() if k not in _MECPR)), ether)


def _meta_ok(lb: Labels, p2: _L2Parsed | None) -> bool:
    cm = lb.chains_meta
    if not lb.L2:
        return cm == []
    if not isinstance(cm, list) or not all(isinstance(m, dict) for m in cm):
        return False
    seq = [(m.get("C"), m.get("DB")) for m in cm]
    if p2 is None:
        return True                                       # L2 本身已记 L2_format
    if seq != [tuple(c) for c in p2[1]]:
        return False
    if not lb.L3:
        return all(m.get("toks") == {} and m.get("ether") == "" for m in cm)
    parsed: list[str] = []
    p3 = _parse_l3(lb.L3, parsed)
    if p3 is None or len(p3[1]) != len(cm):
        return True                                       # L3 文法问题已另记
    for (e, c, d, toks, p), m in zip(p3[1], cm):
        if m.get("toks") != toks or m.get("ether") != ("P" if p else ("O" if e else "")):
            return False
    return True


# ==========================================================================
# 降级
# ==========================================================================
def check_and_downgrade_flags(lb: Labels, it: Interp) -> tuple[Labels, tuple[str, ...]]:
    """见模块文档。"""
    viol = check_labels(lb, it)
    if not viol:
        return lb, ()
    if any(v.startswith("L1_") for v in viol):
        raise Reject("internal_error", f"L1 self-check failed {viol}: {lb.L1!r}")
    cur, flags = lb, []
    if cur.L3:
        cur = replace(cur, L3="", level=2, chains_meta=chains_meta_l2(cur.chains_meta))
        flags.append("selfcheck_downgraded_L3")
        viol = check_labels(cur, it)
    if viol and cur.L2:
        cur = replace(cur, L2="", L3="", level=1, chains_meta=[])
        flags.append("selfcheck_downgraded_L2")
        viol = check_labels(cur, it)
    if viol:
        raise Reject("internal_error", f"self-check failed after downgrade {viol}: {cur.L1!r}")
    return cur, tuple(flags)


def check_and_downgrade(lb: Labels, it: Interp) -> Labels:
    """有违规 -> 先丢 L3 再丢 L2；``L1`` 违规或救不回 -> ``Reject("internal_error")``。要 flags 用 ``check_and_downgrade_flags``。"""
    return check_and_downgrade_flags(lb, it)[0]


def downgrade_flags(before: Labels, after: Labels) -> tuple[str, ...]:
    """降级前后的两个 ``Labels`` -> ``selfcheck_downgraded_*`` flags。"""
    out: list[str] = []
    if before.L3 and not after.L3:
        out.append("selfcheck_downgraded_L3")
    if before.L2 and not after.L2:
        out.append("selfcheck_downgraded_L2")
    return tuple(out)


def drop_idle_flags(flags: Iterable[str], lb: Labels) -> tuple[str, ...]:
    """去掉空转的 L3 断言 flag：最终标签没有 L3 时丢 ``L3_ONLY_FLAGS``，其余原样（排序后返回）。"""
    out = sorted(set(flags))
    if lb.L3:
        return tuple(out)
    return tuple(f for f in out if f not in L3_ONLY_FLAGS)
