"""M5b（链词法）：把单条链的文本词法成 ``RawChain``（纯句法，无类语义）。

对外接口
--------
``lex_chain(tok: str) -> RawChain``
    一条链（不含分隔符）。覆盖：链前缀 ``O- P- dO- N- d t m iso anteiso ai methyl-``；``C:DB``；双键位置
    ``(9Z,12Z)``（丢弃）；醚后缀 ``e``/``p``；支链 ``-15me`` ``(15-me)`` ``(methyl)``（只记 ``branch``，C/DB 不变）；
    氧/修饰记号 ``;atom`` 与 LipidMatch/SwissLipids 的括号形式；``0:0`` 占位；Riken 介质链 ``12-HETE``。
    链上若带嵌套酰基 ``(FA c:d)``/``;FA c:d``，``lex_chain`` 抛 ``malformed``（嵌套链要由 ``lex_chain_ex`` 取回）。
``lex_chain_ex(tok: str) -> tuple[RawChain, tuple[RawChain, ...]]``
    同上，另返回该链带的嵌套酰基链（``Lexed.nested_fa`` 的来源）。
``scan_chain(tok) -> tuple[RawChain, tuple[RawChain, ...], int]``
    ``lexer.lex_body`` 用：第三项 = 链末尾以 ``;(On)`` 书写的未分型氧记号个数（数据集 L2 写法，不挂链，
    由 ``lex_body`` 对**最后一条链**挪进 ``Lexed.tail_atoms``）。
``scan_tail(rest: str) -> tuple[tuple[Atom, ...], tuple[RawChain, ...]]``
    外层括号之后的尾部记号（``;O2;G`` 之类）。
``MAX_DIGITS``  ``C``、``DB``、各记号计数最多 3 位数字（防 ``int()`` 长数字、同时让 ``16:017:0`` 判成粘连）。

记号表（``;`` 后；n 缺省 1）——``;(2OH)``（位置，1 个） ≠ ``;(O2)``/``;O2``/``;2O``（个数） ≠ ``;(OH)2``（个数）
----------------------------------------------------------------------------------------------
``;O`` ``;On`` ``;nO`` ``;(On)``            -> ``O_count`` n（未分型氧个数）
``;OH`` ``;OHn`` ``;(OH)`` ``;(OH)n``         -> ``OH`` n（已分型羟基个数）
``;(kOH)`` ``;kOH`` ``(kOH)`` ``-kOH`` ``;(aOH,bOH)`` -> ``OH_pos``（每项 1 个，n = 项数）
``;oxo[n]`` ``;Ep[n]`` ``;COOH[n]`` ``;OOH[n]`` ``;OO[n]``   -> 同名 kind，n 个
``;Me[n]`` ``;cPr[n]``                        -> ``Me`` / ``cPr``
``;P``                                        -> ``P``（缩醛）
``;Gly ;Tau ;S ;Hex ;G ;T`` 与氨基酸 ``;Phe`` …（``types.AA_CONJ_NAMES``）-> ``conj_*``
``;FA c:d``                                   -> 嵌套酰基（不进 atoms）
LipidMatch 链后括号（无 ``;``）：``(OH)``->OH；``(O)`` ``(OO)`` ``(OOH)`` ``(OOO)`` ``(nO)``->``O_count``（只计氧数）；
``(Ke)`` -> ``oxo``、``(CHO)`` -> ``oxo``、``(COOH)`` -> ``COOH``，``Atom.via`` 记原写法（2026-10-31：这个 C=O 不在字面 DB 里，语义层加 1 个 DB；
库 m/z 核对见 docs/INPUT_FIXES.md；只认链后括号，``;(Ke)`` 仍不支持）；
``(KeKe)`` ``(Ke,OH)`` ``(KeO)`` ``(1O(1Cyc))`` 等 -> ``unsupported_chain_modifier``（D10）；
炔 ``(12Y)`` -> ``unsupported_syntax``。

抛出（全部是 ``Reject``，不返回半成品）：``glued_tokens dangling_sep unbalanced_parens db_ge_c dup_atom multi_ether
class_in_chain malformed unsupported_syntax unsupported_chain_modifier``。``db_ge_c``：非占位链 ``DB>0 且 DB>=C``。
"""
from __future__ import annotations

import re
from typing import Final, NamedTuple

from .types import CONJ_NAMES, LIPIDMATCH_OX_MODIFIERS, OXYLIPIN_ABBREVS, Atom, RawChain, Reject

__all__ = ["lex_chain", "lex_chain_ex", "scan_chain", "scan_tail", "MAX_DIGITS"]

MAX_DIGITS: Final = 3

_CORE: Final = re.compile(r"(\d{1,3}):(\d{1,3})", re.ASCII)
_OXY: Final = re.compile(r"\d{1,3}(?:,\d{1,3})*-(" + "|".join(OXYLIPIN_ABBREVS) + r")", re.ASCII)
_PRE_TOK: Final = re.compile(r"N-|dO-|tO-|eO-|O-|P-|anteiso-?|iso-?|ai-?|methyl-|[dtm]", re.ASCII)
_PAREN: Final = re.compile(r"[()]", re.ASCII)
_BRANCH_NAME: Final = {"iso": "iso", "iso-": "iso", "anteiso": "anteiso", "anteiso-": "anteiso",
                       "ai": "anteiso", "ai-": "anteiso", "methyl-": "methyl"}

# ---- 括号组内容 ------------------------------------------------------------------------------------------
_DBPOS: Final = re.compile(r"\d{1,3}[EZ]?(?:,\d{1,3}[EZ]?)*", re.ASCII)
_ALKYNE: Final = re.compile(r"\d{1,3}[EZY]?(?:,\d{1,3}[EZY]?)*", re.ASCII)
_FA: Final = re.compile(r"FA (\d{1,3}):(\d{1,3})(?:\(\d{1,3}[EZ]?(?:,\d{1,3}[EZ]?)*\))?", re.ASCII)   # 0.4.0：LIPID MAPS 的 (FA 18:1(9Z)) 带双键位置（丢弃）
_OHPOS: Final = re.compile(r"\d{1,3}OH(?:,\d{1,3}OH)*", re.ASCII)
_MEPOS: Final = re.compile(r"(\d{1,3})-?[Mm]e", re.ASCII)
_O_RUN: Final = re.compile(r"(O+)H?", re.ASCII)                      # (O) (OO) (OOH) (OOO) (OOOO)
_O_N_FIRST: Final = re.compile(r"(\d{1,3})O", re.ASCII)              # (3O)
_O_N_LAST: Final = re.compile(r"O(\d{1,3})", re.ASCII)               # (O2)
_WORDY: Final = re.compile(r"[A-Za-z0-9,\[\]\-]+", re.ASCII)
_CYC: Final = re.compile(r"[A-Za-z0-9,]*O\(\d{0,3}Cyc\)", re.ASCII)  # (1O(1Cyc)) (NO(NCyc))

# ---- ';' 记号与后缀 ----------------------------------------------------------------------------------------
_ATOM_ALPHA: Final = re.compile(r"([A-Za-z]+)(\d{0,3})", re.ASCII)
_ATOM_NUM: Final = re.compile(r"(\d{1,3})([A-Za-z]+)", re.ASCII)
_SUF_DASH: Final = re.compile(r"-(\d{1,3})(?:(me)(?![A-Za-z])|OH(?![A-Za-z]))", re.ASCII)
_COUNT: Final = re.compile(r"\d{1,3}", re.ASCII)
#: 括号氧组（``(OH)`` ``(OOH)`` ``(OOOO)``）之后紧跟的个数后缀 ``(OOH)2``：数字后面不再跟 ``:`` / 数字（跟 ``:`` 才是粘连的下一条链 ``(OH)16:0``）
_COUNT_SUFFIX: Final = re.compile(r"\d{1,3}(?![\d:])", re.ASCII)
_COUNTABLE: Final = frozenset({"OH", "O_count"})
#: 链后面解析不了的大写缩写修饰词（LipidMatch AcCar 的 ``14:1-OH`` 的 ``-OH``、``14:1DC`` 的 ``DC``）：未知修饰记号，口径同 ``;XYZ``；
#: 小写 / 混杂的残片（``18:1abc``）不算修饰记号，仍是 ``malformed``
_MODIFIER_WORD: Final = re.compile(r"-?[A-Z]{1,4}", re.ASCII)
_TYPED: Final = frozenset({"OH", "oxo", "Ep", "COOH", "OOH", "OO", "Me", "cPr"})
_CONJ: Final = frozenset(CONJ_NAMES)

# ---- 失败分类 ----------------------------------------------------------------------------------------------
_CLASS_IN_CHAIN: Final = re.compile(
    r"[A-Za-z][A-Za-z0-9+\-]*(?:\s+[A-Za-z][A-Za-z0-9+\-]*)*\s+(?:[dtm]|[OPN]-)?\d{1,3}:\d{1,3}", re.ASCII)
_OACYL_TAIL: Final = re.compile(r"-O-\d{1,3}:\d{1,3}", re.ASCII)


class _Item(NamedTuple):
    """一个括号组的分类结果。``tag`` ∈ ``dbpos`` ``fa``（n=C, m=DB）``branch``（text）``atom``（kind, n；``via`` = LipidMatch 氧化修饰的原写法）。"""
    tag: str
    kind: str = ""
    n: int = 0
    m: int = 0
    text: str = ""
    via: str = ""


class _Acc:
    """一条链扫描过程中的累加器（只在本模块内用）。"""
    __slots__ = ("atoms", "tail", "nested", "suf", "branch")

    def __init__(self, branch: str = "") -> None:
        self.atoms: list[Atom] = []
        self.tail: list[bool] = []          # 与 atoms 同长：该记号是否以 ';(On)' 书写（整分子氧，不挂链）
        self.nested: list[RawChain] = []
        self.suf = ""
        self.branch = branch

    def add(self, kind: str, n: int = 1, *, paren_o: bool = False, via: str = "") -> None:
        self.atoms.append(Atom(kind, n, via))
        self.tail.append(paren_o)

    def set_branch(self, new: str) -> None:
        if self.branch:
            raise Reject("malformed", f"two branch markers ({self.branch!r}, {new!r})")
        self.branch = new


def _is_digit(c: str) -> bool:
    return "0" <= c <= "9"


def _paren_atom(content: str) -> tuple[str, int] | None:
    """括号内容 -> ``(kind, n)``（不含 dbpos / FA / 支链）；识别不了返回 None。"""
    if content == "OH":
        return "OH", 1
    if _OHPOS.fullmatch(content):
        return "OH_pos", content.count("OH")
    m = _O_RUN.fullmatch(content)
    if m:
        return "O_count", len(m.group(1))
    m = _O_N_FIRST.fullmatch(content) or _O_N_LAST.fullmatch(content)
    if m:
        n = int(m.group(1))
        if n < 1:
            raise Reject("malformed", f"zero oxygen count in ({content})")
        return "O_count", n
    return None


def _paren_item(content: str) -> _Item:
    """分类一个括号组；不支持/畸形直接抛。"""
    if _DBPOS.fullmatch(content):
        return _Item("dbpos")
    m = _FA.fullmatch(content)
    if m:
        return _Item("fa", n=int(m.group(1)), m=int(m.group(2)))
    if content == "methyl":
        return _Item("branch", text="methyl")
    m = _MEPOS.fullmatch(content)
    if m:
        return _Item("branch", text=f"{int(m.group(1))}me")
    at = _paren_atom(content)
    if at is not None:
        return _Item("atom", kind=at[0], n=at[1])
    if content in LIPIDMATCH_OX_MODIFIERS:                               # (Ke) (CHO) (COOH)：LipidMatch 氧化磷脂的 C=O 修饰
        return _Item("atom", kind=LIPIDMATCH_OX_MODIFIERS[content], n=1, via=content)
    if "Y" in content and _ALKYNE.fullmatch(content):
        raise Reject("unsupported_syntax", f"alkyne position ({content})")
    if _WORDY.fullmatch(content) or _CYC.fullmatch(content):
        raise Reject("unsupported_chain_modifier", f"({content})")
    raise Reject("malformed", f"cannot parse parenthesised group ({content})")


def _close_paren(s: str, i: int) -> int:
    """``s[i] == '('`` 的配对 ``)`` 下标；不配对抛 ``unbalanced_parens``。"""
    depth = 0
    for m in _PAREN.finditer(s, i):
        if m.group() == "(":
            depth += 1
        else:
            depth -= 1
            if depth == 0:
                return m.start()
    raise Reject("unbalanced_parens", f"unclosed '(' in {s!r}")


def _nested(c: int, db: int, raw: str) -> RawChain:
    if c < 1:
        raise Reject("malformed", f"nested acyl chain {raw!r} has no carbons")
    if db > 0 and db >= c:
        raise Reject("db_ge_c", f"{raw}: DB {db} >= C {c}")
    return RawChain(raw, "", "", c, db, "")


def _fail_rest(rest: str, whole: str) -> None:
    """链上解析不了的尾巴 -> 对应的拒绝（永远抛出）。"""
    c0 = rest[0]
    if _is_digit(c0) or c0 == ":":
        raise Reject("glued_tokens", f"{rest!r} glued to a chain in {whole!r}")
    if c0 == ")":
        raise Reject("unbalanced_parens", f"stray ')' in {whole!r}")
    if _OACYL_TAIL.match(rest):
        raise Reject("unsupported_syntax", f"omega-O-acyl tail {rest!r}")
    if _CLASS_IN_CHAIN.match(rest):
        raise Reject("class_in_chain", f"class name inside chain: {rest!r}")
    if _MODIFIER_WORD.fullmatch(rest):
        raise Reject("unsupported_chain_modifier", f"unknown modifier {rest!r} after the chain in {whole!r}")
    raise Reject("malformed", f"cannot parse {rest!r} in chain {whole!r}")


def _scan_atom(s: str, i: int, acc: _Acc) -> int:
    """解析 ``s[i] == ';'`` 开头的一个记号，追加到 ``acc``；返回新下标。"""
    j = i + 1
    if j >= len(s):
        raise Reject("malformed", f"dangling ';' in {s!r}")
    c = s[j]
    if c == "(":
        k = _close_paren(s, j)
        content = s[j + 1:k]
        item = _paren_item(content)
        k += 1
        if item.tag == "atom" and item.via:                               # ``;(Ke)``：只认链后括号 ``18:1(Ke)``，不放宽
            raise Reject("unsupported_chain_modifier", f";({content})")
        if item.tag == "atom":
            n = item.n
            if item.kind == "OH" and n == 1:                              # ;(OH)2 -> 个数
                m = _COUNT.match(s, k)
                if m:
                    n = int(m.group())
                    k = m.end()
                    if n < 1:
                        raise Reject("malformed", "zero count after (OH)")
            acc.add(item.kind, n, paren_o=item.kind == "O_count")        # ;(On) = 整分子氧
        elif item.tag == "fa":
            acc.nested.append(_nested(item.n, item.m, f"({content})"))
        else:
            raise Reject("malformed", f"';({content})' is not an atom")
        return k
    if _is_digit(c):
        m = _ATOM_NUM.match(s, j)
        if not m:
            raise Reject("malformed", f"bad atom after ';' in {s!r}")
        n, name = int(m.group(1)), m.group(2)
        if n < 1:
            raise Reject("malformed", f"zero count in ;{m.group()}")
        if name == "O":
            acc.add("O_count", n)
        elif name == "OH":
            acc.add("OH_pos")
        else:
            raise Reject("unsupported_syntax", f"position-prefixed modifier ;{m.group()}")
        return m.end()
    m = _ATOM_ALPHA.match(s, j)
    if not m:
        raise Reject("malformed", f"bad atom after ';' in {s!r}")
    name, digits = m.group(1), m.group(2)
    end = m.end()
    if end < len(s) and s[end] == ":":
        raise Reject("glued_tokens", f"{s[j:]!r} looks like a glued chain")
    if name == "FA":
        mf = _FA.match(s, j)
        if not mf:
            raise Reject("malformed", f"bad ';FA' in {s!r}")
        acc.nested.append(_nested(int(mf.group(1)), int(mf.group(2)), mf.group()))
        return mf.end()
    if name == "O" or name in _TYPED:
        cnt = int(digits) if digits else 1
        if cnt < 1:
            raise Reject("malformed", f"zero count in ;{name}{digits}")
        acc.add("O_count" if name == "O" else name, cnt)
    elif name == "P" or name in _CONJ:
        if digits:
            raise Reject("malformed", f";{name} takes no count")
        acc.add(name if name == "P" else "conj_" + name)
    else:
        raise Reject("unsupported_chain_modifier", f"unknown modifier ;{name}{digits}")
    return end


def _scan_items(s: str, i: int, acc: _Acc, allow_suffix: bool) -> int:
    """从 ``s[i]`` 起吃掉 ``( ; - e p`` 开头的修饰项，直到吃不动；返回停下的下标。

    ``allow_suffix=False``（外层括号之后的尾部）时只认括号组与 ``;`` 记号。
    """
    n = len(s)
    while i < n:
        c = s[i]
        if c == "(":
            j = _close_paren(s, i)
            content = s[i + 1:j]
            item = _paren_item(content)
            i = j + 1
            if item.tag == "atom":
                if item.via and _COUNT_SUFFIX.match(s, i):
                    raise Reject("unsupported_chain_modifier", f"count suffix after ({content}): {s[j:]!r}")   # (Ke)2 / (COOH)2：没核实 -> 不支持
                if item.kind in _COUNTABLE and _COUNT_SUFFIX.match(s, i):
                    # LipidMatch ``18:2(OOH)2`` / ``(OH)2``：氧括号后面的个数后缀（不是粘连的链记号）。语义没核实 -> 不支持，不猜
                    raise Reject("unsupported_chain_modifier", f"count suffix after the bracketed oxygen group ({content}): {s[j:]!r}")
                acc.add(item.kind, item.n, via=item.via)
            elif item.tag == "fa":
                acc.nested.append(_nested(item.n, item.m, f"({content})"))
            elif item.tag == "branch":
                acc.set_branch(item.text)
        elif c == ";":
            i = _scan_atom(s, i, acc)
        elif not allow_suffix:
            break
        elif c == "-":
            m = _SUF_DASH.match(s, i)
            if not m:
                break
            if m.group(2):                                               # -15me
                acc.set_branch(f"{int(m.group(1))}me")
            else:                                                        # -2OH
                acc.add("OH_pos")
            i = m.end()
        elif c == "e" or c == "p":
            if i + 1 < n and s[i + 1] not in ";()-, ":
                break
            if acc.suf:
                raise Reject("malformed", f"two ether suffixes in {s!r}")
            acc.suf = c
            i += 1
        else:
            break
    return i


def _check_dups(atoms: list[Atom], where: str) -> None:
    if len(atoms) > 1:
        seen: set[str] = set()
        for a in atoms:
            if a.kind in seen:
                raise Reject("dup_atom", f"repeated {a.kind} {where}")
            seen.add(a.kind)


def scan_chain(tok: str) -> tuple[RawChain, tuple[RawChain, ...], int]:
    """词法一条链；返回 ``(RawChain, 嵌套酰基链, 末尾 ;(On) 记号个数)``。见模块文档。"""
    s = tok.strip()
    if not s:
        raise Reject("malformed", "empty chain")
    # ---- 前缀 ----
    pre = base = branch = ""
    n_acyl = False
    i = 0
    if not _is_digit(s[0]):
        while True:
            m = _PRE_TOK.match(s, i)
            if not m:
                break
            t = m.group()
            if t == "N-":
                if n_acyl:
                    raise Reject("malformed", f"two N- prefixes in {s!r}")
                n_acyl = True
            elif t in ("O-", "P-", "dO-"):
                if pre:
                    raise Reject("multi_ether", f"two ether prefixes in {s!r}")
                pre = t
            elif t in ("tO-", "eO-"):
                raise Reject("unsupported_syntax", f"ether prefix {t!r}")
            elif t in ("d", "t", "m"):
                if base:
                    raise Reject("malformed", f"two sphingoid prefixes in {s!r}")
                base = t
            else:
                if branch:
                    raise Reject("malformed", f"two branch markers in {s!r}")
                branch = _BRANCH_NAME[t]
            i = m.end()
    # ---- C:DB ----
    m = _CORE.match(s, i)
    if m is None:
        mo = _OXY.fullmatch(s) if i == 0 else None
        if mo is not None:
            return RawChain(s, "", "", 0, 0, "", oxylipin=mo.group(1)), (), 0
        if _CLASS_IN_CHAIN.match(s):
            raise Reject("class_in_chain", f"class name inside chain: {s!r}")
        raise Reject("malformed", f"bad chain {s!r}")
    C, DB = int(m.group(1)), int(m.group(2))
    end_core = m.end()
    acc = _Acc(branch)
    i = _scan_items(s, end_core, acc, True)
    if i < len(s):
        _fail_rest(s[i:], s)
    # ---- 校验 ----
    if C == 0 and DB == 0:
        if pre or base or n_acyl or acc.branch or acc.suf or acc.atoms or acc.nested or len(s) != end_core:
            raise Reject("malformed", f"placeholder 0:0 cannot carry modifiers: {s!r}")
        return RawChain(s, "", "", 0, 0, "", placeholder=True), (), 0
    if DB > 0 and DB >= C:
        raise Reject("db_ge_c", f"{s}: DB {DB} >= C {C}")
    n_co = sum(a.n for a in acc.atoms if a.via)                        # (Ke) (CHO) (COOH) 的 C=O 也占不饱和度：字面 DB + C=O 个数 < C
    if n_co and len(acc.atoms) > 1:                                    # 一条链上 (Ke) 之外还有别的氧记号（(Ke)(OH)、(CHO);O）：库 m/z 只核对了单个修饰 -> 不支持
        raise Reject("unsupported_chain_modifier", f"bracketed Ke/CHO/COOH combined with another oxygen mark on one chain: {s!r}")
    if n_co and DB + n_co >= C:
        raise Reject("db_ge_c", f"{s}: DB {DB} + {n_co} C=O from the bracketed modifier >= C {C}")
    _check_dups(acc.atoms, f"on one chain: {s!r}")
    n_tail = 0
    while n_tail < len(acc.tail) and acc.tail[len(acc.tail) - 1 - n_tail]:
        n_tail += 1
    rc = RawChain(s, pre, base, C, DB, acc.suf, tuple(acc.atoms), n_acyl, acc.branch)
    return rc, tuple(acc.nested), n_tail


def lex_chain_ex(tok: str) -> tuple[RawChain, tuple[RawChain, ...]]:
    """见模块文档：``(RawChain, 嵌套酰基链)``。"""
    rc, nested, _ = scan_chain(tok)
    return rc, nested


def lex_chain(tok: str) -> RawChain:
    """词法一条链（见模块文档）。带嵌套 ``(FA c:d)`` 的链用 ``lex_chain_ex``。"""
    rc, nested, _ = scan_chain(tok)
    if nested:
        raise Reject("malformed", f"{tok!r} carries a nested (FA) chain; use lex_chain_ex")
    return rc


def scan_tail(rest: str) -> tuple[tuple[Atom, ...], tuple[RawChain, ...]]:
    """外层括号之后的尾部（以 ``;`` 或括号氧开头）-> ``(atoms, nested)``；解析不了抛 ``malformed``。"""
    s = rest.strip()
    acc = _Acc()
    i = _scan_items(s, 0, acc, False)
    if i < len(s) or acc.branch:
        raise Reject("malformed", f"cannot parse tail {s[i:]!r}")
    _check_dups(acc.atoms, f"in tail: {s!r}")
    return tuple(acc.atoms), tuple(acc.nested)
