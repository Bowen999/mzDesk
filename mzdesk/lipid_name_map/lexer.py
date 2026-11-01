"""M5a / M5b：切头体与整名词法（纯句法，不做类语义：类选择、醚种类合并、氧分配归 M4/M6）。

对外接口
--------
``split_head(text: str) -> tuple[str, str]``
    找第一个 ``C:DB``（连同紧贴它的链前缀 ``N- O- P- dO- d t m iso anteiso ai methyl-``），其前为 ``head``
    （去掉末尾空白与形式 B 的 ``(``；保留 ``Hex(4)-HexNAc(2)-Cer``、``(Fuc)iGb3Cer`` 这类类括号），其余为 ``body``。
    形式 B（括号）的 ``body`` 以 ``(`` 开头：``"PC(16:0/18:1)" -> ("PC", "(16:0/18:1)")``；
    形式 A 的 ``body`` 即链文本：``"PC 16:0/18:1" -> ("PC", "16:0/18:1")``；
    ``AHexCer (O-24:6)18:1;2O/26:0;O -> ("AHexCer", "(O-24:6)18:1;2O/26:0;O")``（括号 O-酰基前缀由 ``lex_body`` 识别）。
    前缀只在其前一字符为 行首/空格/``(``/``/``/``_`` 时才算链前缀（``Cer-d 24:1`` 的 ``d`` 属于类名）。
    抛 ``malformed``（无 ``C:DB``、head 为空、第一条链落在嵌套括号里而整串括号其实配对）、``dangling_sep``（head 以分隔符结尾）、
    ``unbalanced_parens``（整串括号真的不配对）。
``lex_body(head: str, body: str, cfg: Config | None = None) -> Lexed``
    括号/空格两种形式；分隔符 ``/ _``（两侧允许空白）；``0:0`` 占位链；括号 O-酰基前缀 -> ``Lexed.pre_acyl``；
    相邻单链括号组 ``CL(28:0)(28:0)``（LipidMatch 的 CL 对记法，D24）= 以 ``_`` 连接的链列表，每条链的 ``RawChain.raw`` 保留括号
    （语义层靠它认出这种写法；只有 CL 对记法有定义，其余类 core 拒绝 ``malformed``）；组内含分隔符 / ``(FA c:d)`` / 氧括号的不算；
    嵌套 ``(FA c:d)`` -> ``Lexed.nested_fa``；末尾的 ``;(On)``（数据集 L2 写法）与外层括号之后的 ``;atom`` ->
    ``Lexed.tail_atoms``；``Lexed.feats`` 的每个字段（``Lexed.problems()`` 恒为空）。
``lex_chain(tok: str) -> RawChain``  见 ``chains.py``（此处 re-export）。
``lex(text: str, cfg: Config | None = None) -> Lexed``  便利函数 = ``split_head`` + ``lex_body``。

所有失败抛 ``Reject``（``malformed glued_tokens dangling_sep unbalanced_parens db_ge_c dup_atom multi_ether
all_zero_chains class_in_chain unsupported_syntax unsupported_chain_modifier``），不返回半成品。``cfg`` 目前不被词法使用
（介质缩写开关 ``allow_oxylipin_abbrev`` 由语义层检查）。

约定（与 ``types.py`` 一致）：``is_sum`` 纯句法 = 恰好一个链组；``n_chains`` 不含占位/``pre_acyl``/``nested_fa``；
``ether_marks`` 按书写顺序逐链收 ``pre`` 再收 ``suf``，总数 >= 2 即 ``multi_ether``；``tail_conj`` = 共轭名（D17）：一种原样，
``;G;S`` / ``;T;S``（任意顺序）-> ``'G+S'`` / ``'T+S'``，其他多共轭组合抛 ``malformed``（``types.conj_key``）；
``eos2O`` = 第 2 条非占位链原文含数字在前的 ``;2O``，**且**第 1、2 条非占位链之间的分隔符全是 ``/``（P1：``Cer 24:1_18:1;2O``
不触发；纯词法，只有 Cer/HexCer 头联用才有 EOS 含义，D18）。
"""
from __future__ import annotations

import dataclasses
import re
from typing import Final

from .chains import lex_chain, scan_chain, scan_tail
from .config import Config
from .types import Atom, Features, Lexed, RawChain, Reject, conj_key

__all__ = ["split_head", "lex_body", "lex_chain", "lex"]

_CORE_SEARCH: Final = re.compile(r"(?<![\d:.])(\d{1,3}):(\d{1,3})", re.ASCII)
_BOUNDARY: Final = " (/_\\"
_PRE_END: Final = re.compile(
    r"(?<![^ (/_\\])(?:N-|dO-|tO-|eO-|O-|P-|anteiso-?|iso-?|ai-?|methyl-|[dtm])+$", re.ASCII)
_PRE_ACYL: Final = re.compile(
    r"\(O-(\d{1,3}):(\d{1,3})\)(?=\d|[dtm]\d|dO-|[OPN]-|iso|anteiso|ai\d|methyl-)", re.ASCII)
_SPECIAL: Final = re.compile(r"[()/_]", re.ASCII)
#: 相邻括号组 ``(a)(b)`` 里的一个组：单条链（``[前缀]C:DB…``），不含 ``FA c:d`` 嵌套、``(OH)`` 氧括号
_CHAIN_GROUP: Final = re.compile(r"(?:dO-|tO-|eO-|[OPN]-|[dtm])?\d{1,3}:\d{1,3}", re.ASCII)
_PAREN: Final = re.compile(r"[()]", re.ASCII)
_EOS2O: Final = re.compile(r";\d{1,3}O(?![A-Za-z])", re.ASCII)
_CLASS_SUFFIXES: Final = frozenset({"-O", "-P", "-d", "-t", "-m"})
_PREFIX_MARKS: Final = frozenset({"O-", "P-", "dO-"})


def split_head(text: str) -> tuple[str, str]:
    """切头体，见模块文档。"""
    s = text.strip()
    m = _CORE_SEARCH.search(s)
    if m is None:
        raise Reject("malformed", "no C:DB chain found")
    p = m.start()
    if p == 0 or s[p - 1] in _BOUNDARY:
        q = p
    else:
        pm = _PRE_END.search(s, 0, p)
        q = pm.start() if pm else p
    hr = s[:q].rstrip()
    if hr.endswith("("):
        head = hr[:-1].rstrip()
        body = s[len(hr) - 1:]
    else:
        head = hr
        body = s[q:]
    if not head:
        raise Reject("malformed", "empty class head")
    if head[-1] in "/_\\":
        raise Reject("dangling_sep", f"class head {head!r} ends with a separator")
    if head.count("(") != head.count(")"):
        if _parens_balanced(s):
            # 整串括号是配对的，只是第一条链落在嵌套括号里（LipidMatch ``OxCL((34:2(OH))_(44:9))+NH4``、``CL(1'-[18:1/18:2],3'-[…])``）：
            # 括号没有不配对，是头与链体切不开（``malformed``：整体语法不成形，无法切分头部与链体）
            raise Reject("malformed", f"the first chain sits inside a nested bracket, the class head cannot be cut: {head!r}")
        raise Reject("unbalanced_parens", f"unbalanced parentheses in class head {head!r}")
    return head, body.strip()


def _parens_balanced(s: str) -> bool:
    """``(`` / ``)`` 整串配对（深度从不为负、最后为 0）。"""
    depth = 0
    for m in _PAREN.finditer(s):
        depth += 1 if m.group() == "(" else -1
        if depth < 0:
            return False
    return depth == 0


def _match_close(s: str, i: int) -> int:
    """``s[i] == '('`` 的配对 ``)``；无则抛 ``unbalanced_parens``。"""
    depth = 0
    for m in _PAREN.finditer(s, i):
        if m.group() == "(":
            depth += 1
        else:
            depth -= 1
            if depth == 0:
                return m.start()
    raise Reject("unbalanced_parens", f"unclosed '(' in {s!r}")


def _adjacent_groups(s: str) -> tuple[list[str], str]:
    """``s`` 以 ``(`` 开头：取开头连续的"单链括号组" ``(a)(b)…``，返回 ``(各组内文, 其后剩余文本)``。

    只有 **>= 2 个** 紧挨着的单链组才算（CL 对记法 ``CL(28:0)(28:0)``，D24）；一个组 / 组内含 ``/ _`` / ``(FA c:d)`` 嵌套 / ``(OH)``
    氧括号都不算，返回 ``([], s)``，由调用方走老路径。括号不配对抛 ``unbalanced_parens``。
    """
    if s[_match_close(s, 0) + 1:][:1] != "(":                # 快速路径：第一个组之后不是紧挨着的 "("（绝大多数括号写法）
        return [], s
    groups: list[str] = []
    i = 0
    while i < len(s) and s[i] == "(":
        j = _match_close(s, i)
        g = s[i + 1:j]
        if _CHAIN_GROUP.match(g) is None or len(_split_top(g)[0]) != 1:
            break
        groups.append(g)
        i = j + 1
    if len(groups) < 2:
        return [], s
    return groups, s[i:].strip()


def _split_top(inner: str) -> tuple[list[str], list[str]]:
    """在括号深度 0 处按 ``/``、``_`` 切开；括号不配对抛 ``unbalanced_parens``。"""
    pieces: list[str] = []
    seps: list[str] = []
    depth = 0
    start = 0
    for m in _SPECIAL.finditer(inner):
        c = m.group()
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
            if depth < 0:
                raise Reject("unbalanced_parens", f"stray ')' in {inner!r}")
        elif depth == 0:
            pieces.append(inner[start:m.start()])
            seps.append(c)
            start = m.end()
    if depth:
        raise Reject("unbalanced_parens", f"unclosed '(' in {inner!r}")
    pieces.append(inner[start:])
    return pieces, seps


def lex_body(head: str, body: str, cfg: Config | None = None) -> Lexed:
    """整名词法，见模块文档。"""
    del cfg                                                    # 词法层不读配置
    s = body.strip()
    pre_acyl: RawChain | None = None
    rest = ""
    adjacent = False
    m = _PRE_ACYL.match(s)
    if m is not None:
        c, db = int(m.group(1)), int(m.group(2))
        if c < 1:
            raise Reject("malformed", "O-acyl prefix without carbons")
        if db > 0 and db >= c:
            raise Reject("db_ge_c", f"{m.group()}: DB {db} >= C {c}")
        pre_acyl = RawChain(m.group(), "", "", c, db, "")
        form, inner = "space", s[m.end():]
    elif s.startswith("("):
        groups, after = _adjacent_groups(s)
        if groups:                                             # D24：``(a)(b)`` 相邻组 = 以 ``_`` 连接的链列表（raw 保留括号作标记）
            form, inner, rest = "paren", "_".join(groups), after
            adjacent = True
        else:
            j = _match_close(s, 0)
            form, inner, rest = "paren", s[1:j], s[j + 1:].strip()
    else:
        form, inner = "space", s
    if not inner.strip():
        raise Reject("malformed", "empty chain list")

    pieces, seps = _split_top(inner)
    chains: list[RawChain] = []
    nested: list[RawChain] = []
    tail: list[Atom] = []
    last = len(pieces) - 1
    for k, piece in enumerate(pieces):
        if not piece.strip():
            raise Reject("dangling_sep", f"empty chain slot in {inner!r}")
        rc, nf, n_tail = scan_chain(piece)
        if adjacent:
            rc = dataclasses.replace(rc, raw=f"({rc.raw})")
        if nf:
            nested.extend(nf)
        if k == last and n_tail:
            tail.extend(rc.atoms[len(rc.atoms) - n_tail:])
            rc = dataclasses.replace(rc, atoms=rc.atoms[:len(rc.atoms) - n_tail])
        chains.append(rc)
    if rest:
        t_atoms, t_nested = scan_tail(rest)
        tail.extend(t_atoms)
        nested.extend(t_nested)

    return _assemble(head, form, chains, seps, pre_acyl, nested, tail)


def _assemble(head: str, form: str, chains: list[RawChain], seps: list[str], pre_acyl: RawChain | None,
              nested: list[RawChain], tail: list[Atom]) -> Lexed:
    n_ph = 0
    marks: list[str] = []
    first_base = ""
    conj: set[str] = set()
    for c in chains:
        if c.placeholder:
            n_ph += 1
        if c.pre:
            marks.append(c.pre)
        if c.suf:
            marks.append(c.suf)
        if c.base and not first_base:
            first_base = c.base
        for a in c.atoms:
            if a.kind[:5] == "conj_":
                conj.add(a.kind[5:])
    for a in tail:
        if a.kind[:5] == "conj_":
            conj.add(a.kind[5:])
    if n_ph == len(chains):
        raise Reject("all_zero_chains", "every chain is a 0:0 placeholder")
    if len(marks) > 1:
        raise Reject("multi_ether", f"ether markers {marks}")
    tail_conj = conj_key(conj)
    if tail_conj is None:
        raise Reject("malformed", f"unsupported conjugate combination {sorted(conj)}")
    if tail:
        seen: set[str] = set()
        for a in tail:
            if a.kind in seen:
                raise Reject("dup_atom", f"repeated {a.kind} in trailing atoms")
            seen.add(a.kind)

    eos2O = False
    real = [i for i, c in enumerate(chains) if not c.placeholder]
    if len(real) >= 2:
        # P1：LipidBlast EOS 只写斜线 ``a;2O/b;2O``；第 1、2 条非占位主链之间的分隔符只要有一个 ``_``，
        # ``;2O`` 就只是 ``;O2`` 的另一种写法（MS-DIAL 规范化名 ``Cer 24:1_18:1;2O``），不触发 EOS。
        i1, i2 = real[0], real[1]
        eos2O = all(sp == "/" for sp in seps[i1:i2]) and _EOS2O.search(chains[i2].raw) is not None

    suffix = head[-2:] if len(head) > 2 and head[-2:] in _CLASS_SUFFIXES else ""
    is_sum = len(chains) == 1
    feats = Features(
        has_ether_prefix=any(mk in _PREFIX_MARKS for mk in marks),
        ether_marks=tuple(marks),
        n_placeholder=n_ph,
        has_fa_paren=bool(nested),
        has_pre_acyl=pre_acyl is not None,
        class_suffix=suffix,
        sphingoid_prefix=first_base,
        eos2O=eos2O,
        tail_conj=tail_conj,
        n_chains=len(chains) - n_ph,
        is_sum=is_sum,
        seps=frozenset(seps),
    )
    return Lexed(head=head, form=form, chains=tuple(chains), seps=tuple(seps), pre_acyl=pre_acyl,
                 nested_fa=tuple(nested), tail_atoms=tuple(tail), is_sum=is_sum, feats=feats)


def lex(text: str, cfg: Config | None = None) -> Lexed:
    """``split_head`` + ``lex_body``。"""
    head, body = split_head(text)
    return lex_body(head, body, cfg)
