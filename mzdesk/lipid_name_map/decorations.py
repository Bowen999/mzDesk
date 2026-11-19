"""M2 装饰 / 加合物 / 别名拆分（SYNTHESIS §2.2 M2、§3.4 S12/S16/S17、D6/D7）。

对外接口
--------
``split_decorations(s: str, cfg: Config) -> list[Segment]``
    输入是 ``preclean`` 之后的字符串；输出 >=1 个 ``Segment``（``types.Segment``），每段已剥掉装饰，交给 M3/M5a。
    不改大小写，不碰链内部写法（``;2O`` ``(FA 16:0)`` ``O-`` ``d18:1`` ``-O`` ``/N-`` ``(methyl)`` 都原样留给词法层）。

``segment_flags(seg: Segment) -> tuple[str, ...]``
    由 ``Segment`` 推出应打的装饰类 flags（``types.FLAGS`` 词表内）：
    ``kind=="alias"`` -> ``alias_split``；``kind=="lipidin_chain"`` -> ``lipidin_chain_used``；``adduct`` 非空 ->
    ``adduct_stripped``；``notes`` 里按 ``NOTE_FLAGS`` 映射（``sn_suffix`` -> ``sn_suffix_stripped``，
    ``ox_prefix`` -> ``ox_prefix_ignored``，``lipidin_sum_O_mismatch`` -> 同名）。排序去重。
    （``api._merge_segments`` 可直接用它；``meta.adduct = seg.adduct``，``meta.stripped = list(seg.notes)``。）

``adduct_info(raw: str) -> AdductInfo | None``
    把一个加合物写法（``[M+H]+`` / ``+NH4-CO`` / ``M-H``）解成结构；认不出返回 ``None``。
    ``AdductInfo.delta_mass``：``m/z = (multimer*M + delta_mass) / |charge|``（单电荷时即 m/z 增量）。

``NOTE_FLAGS``  ``Segment.notes`` 键 -> flag 的映射；``ADDUCT_SPECIES``  加合物物种表（``data/adducts.csv``）。

处理顺序（每个别名段内）
------------------------
1. ``" | "`` 拆别名。含 ``Chain:`` 的段 -> LipidIN ``X | Chain: Y``：取 Y（``kind="lipidin_chain"``）；X 与 Y 的 C:DB 总和
   （不计 ``(FA c:d)`` 嵌套酰基，LipidIN 标题不含它）不同 -> ``Reject("sum_detail_conflict")``（D6）；仅氧数不同 ->
   notes 里加 ``lipidin_sum_O_mismatch``，信 Y。其余多段 -> ``kind="alias"`` 各自独立处理（D7 的合并在 M9）。
2. **先判同位素**：``-d7``、``(d9)``、``[D9]``、``d7-`` 前缀、``[2H7]`` ``[13C6]`` ... -> ``Reject("isotope_label")``。
   这一步在剥加合物之前，所以 ``PC(16:0/18:1)-d7`` 不会被当成 ``)-d7`` 加合物。
   不会误判的写法：PNNL 类后缀 ``Cer-d 24:1_18:1``（``-d`` 后是空格）、鞘氨醇 ``d18:1``、``(d18:1/24:0)``、
   ``iso-d17:1``、``GM3(d26:0)``、``dO-``；维生素 D 名 ``Vitamin D3`` / ``25-hydroxyvitamin D3`` / 裸 ``D3``（编号不是氘标记，
   后续由 M3 拒绝为俗名；``Vitamin D3-d3`` 的 ``-d3`` 仍判同位素）。
3. 自段尾向前反复剥（顺序无关，直到不再变化）：
   * 加合物（至多一个；再出现 -> ``unrecognized_tail``）：
     - 括号式 ``[M+H]+ [M-2H]2- [M+H-H2O]+ [M+2NH4]2+ [M+OAc]- [M-CH3]- [M]+``（LipidDetective / MS-DIAL / SwissLipids）；
     - PNNL mgf ``; [M-H]- C42H82N1O3 M-H``（括号 + 分子式 + 裸码；裸码必须与括号一致）；
     - LipidMatch 后缀 ``)+H )-H )+NH4 )+Na )-2H )+HCO2 )+C2H3O2 )+NH4-CO )-2H+3Na``（紧跟 ``)``）；
     - 裸后缀 ``PC 16:0/18:1+H`` ``TG 16:0_18:1_18:2+NH4`` ``PS 18:0_20:4-H``（2026-10-31：紧跟链末的数字、没有括号；同样每项都要在物种表里）；
     - 裸码 `` M+H``（空格 + M + 加减项）。
     加合物写法由物种表按"项"解析：每项 ``[+-] 个数 物种``，物种须在 ``data/adducts.csv``；认不出的物种 -> ``unrecognized_tail``。
     裸后缀**只认物种表**：``+O`` ``+2O``（LipidIN 的氧化写法）、``-2OH``、``-15me``、``-SN1`` 都不在表里，原样留给后面（不剥、不当加合物）。
     裸后缀只认紧跟数字的（``)`` 的见上一条）；链以字母结尾（``;O+H``、``16:0e+H``）仍不认。（原边界 R9 / B45：``PC 34:1+H`` 曾拒绝，现在剥。）
   * ``-SN1/-SN2/-SN3``、``[sn1]``（notes ``sn_suffix:SN1``）。
   * LIPID MAPS 立体/异构体括号 ``[iso2] [iso3] [iso6] [rac] [R] [S] [U]``（notes ``stereo:iso3``）、
     SwissLipids 变体 ``(a)``..``(g)``（notes ``variant:a``）。仅当前面紧跟 ``)`` 或数字时才剥（``Cer[NS]`` 不动）。
   * 其余形如 ``[...]`` 的尾部括号：若前面已出现 ``C:DB`` -> ``unrecognized_tail``；否则原样保留（如 ``LipidA-PP [10_10_(40)]``）。
4. 段首 ``Ox`` / ``O2x`` 前缀（LipidIN/LipidMatch 的氧化变体，``OxPC`` ``O2xPE`` ``OxFA``）：**仅当**其后有氧标记
   （``;O`` ``;2O`` ``;(2OH)`` 或链后括号 ``(1O) (OH) (Ke) (CHO) (COOH)``）时才剥掉前缀并记 notes ``ox_prefix:Ox``，
   氧数以后缀为准；没有氧标记的 ``OxPC 36:3`` 原样放行，由类表拒绝（D20：``alias/tools.csv`` 里 ``Ox*``/``O2x*`` 头全是
   ``reject unsupported_syntax``，前缀是唯一的氧化信息，丢掉就是静默给错标签）。

不做的事（交给别人）
    S16（MS-DIAL 规范化名的 ``harmonized_name_oxygen_placement``）需要类/链信息，M2 无法判定；见交付说明。

任何剥离都不会把段剥成空串：剥完为空则该步不动，原文交给 M3 拒绝。
所有正则只用 ``[0-9]``，不依赖 Unicode 数字。纯函数、线程安全。
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Final

from .config import Config
from .types import Reject, Segment

__all__ = ["split_decorations", "segment_flags", "adduct_info", "AdductInfo", "Species",
           "NOTE_FLAGS", "ADDUCT_SPECIES"]

# ==========================================================================
# 加合物物种表
# ==========================================================================
ADDUCT_CSV: Final = Path(__file__).resolve().parent / "data" / "adducts.csv"
_ELECTRON: Final = 0.000548579909


@dataclass(frozen=True, slots=True)
class Species:
    """加合物里的一个"项"：``mass`` 是该物种**作为离子/中性分子加入时**的质量（带电物种已含电子修正）。"""
    name: str
    formula: str
    charge: int
    mass: float


def _load_species() -> dict[str, Species]:
    table: dict[str, Species] = {}
    with ADDUCT_CSV.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            sp = Species(row["species"], row["formula"], int(row["charge"]), float(row["mass"]))
            for key in (sp.name, *(a for a in row["aliases"].split("|") if a)):
                if key in table:
                    raise ValueError(f"adducts.csv: duplicate species key {key!r}")
                table[key] = sp
    return table


#: 写法（含别名）-> Species
ADDUCT_SPECIES: Final = MappingProxyType(_load_species())


@dataclass(frozen=True, slots=True)
class AdductInfo:
    """解析后的加合物。``terms`` = ``((符号, 个数, 物种规范名), ...)``。"""
    raw: str
    style: str                       # bracket（[M+H]+）| suffix（+H，LipidMatch）| code（M+H）
    multimer: int
    terms: tuple[tuple[str, int, str], ...]
    charge: int                      # 括号式写明电荷时取之，否则取物种电荷之和
    delta_mass: float                # m/z = (multimer*M + delta_mass) / |charge|


_TERM_RE: Final = re.compile(r"([+-])([0-9]*)([A-Za-z][A-Za-z0-9]*)")


def _parse_terms(body: str) -> list[tuple[str, int, Species]] | None:
    """``+H-H2O`` -> [('+',1,H),('-',1,H2O)]；有任何一项认不出就返回 None。"""
    body = body.replace(" ", "")
    pos, out = 0, []
    while pos < len(body):
        m = _TERM_RE.match(body, pos)
        if m is None:
            return None
        sp = ADDUCT_SPECIES.get(m.group(3))
        if sp is None:
            return None
        out.append((m.group(1), int(m.group(2) or 1), sp))
        pos = m.end()
    return out


def _make_info(raw: str, style: str, multimer: int, terms: list[tuple[str, int, Species]],
               declared: int | None) -> AdductInfo:
    delta = 0.0
    z_species = 0
    for sign, n, sp in terms:
        k = n if sign == "+" else -n
        delta += k * sp.mass
        z_species += k * sp.charge
    z = z_species if declared is None else declared
    delta += (z_species - z) * _ELECTRON          # 补/减电子（如 [M]+ 少一个电子）
    return AdductInfo(raw, style, multimer, tuple((s, n, sp.name) for s, n, sp in terms), z, round(delta, 6))


_BRACKET_BODY: Final = r"\[\s*([0-9]*)M(?=[\]+\-\s])([^\[\]]*)\]\s*([0-9]*)([+-]?)"
_BRACKET_FULL: Final = re.compile(_BRACKET_BODY)


def _info_bracket(token: str) -> AdductInfo | None:
    m = _BRACKET_FULL.fullmatch(token.strip())
    if m is None:
        return None
    digits, sign = m.group(3), m.group(4)
    if digits and not sign:
        return None
    terms = _parse_terms(m.group(2))
    if terms is None:
        return None
    declared = None
    if sign:
        declared = (int(digits) if digits else 1) * (1 if sign == "+" else -1)
    return _make_info(token.strip(), "bracket", int(m.group(1) or 1), terms, declared)


def _info_terms(token: str, style: str, multimer: int = 1, raw: str | None = None) -> AdductInfo | None:
    terms = _parse_terms(token)
    if not terms:
        return None
    return _make_info(raw if raw is not None else token, style, multimer, terms, None)


_CODE_FULL: Final = re.compile(r"([0-9]*)M((?:[+-][0-9]*[A-Za-z][A-Za-z0-9]*)+)")


def adduct_info(raw: str) -> AdductInfo | None:
    """解析单个加合物写法；认不出返回 ``None``。

    >>> adduct_info("[M+H]+").delta_mass
    1.007276
    >>> adduct_info("+NH4-CO").style
    'suffix'
    >>> adduct_info("-d7") is None
    True
    """
    raw = raw.strip()
    if raw.startswith("["):
        return _info_bracket(raw)
    if raw[:1] in ("+", "-"):
        return _info_terms(raw, "suffix")
    m = _CODE_FULL.fullmatch(raw)
    if m is not None:
        return _info_terms(m.group(2), "code", int(m.group(1) or 1), raw)
    return None


# ==========================================================================
# 正则
# ==========================================================================
_CDB: Final = re.compile(r"(?<![0-9:])([0-9]+):([0-9]+)(?![0-9:])")
_NESTED_FA: Final = re.compile(r"\(\s*FA\s+[^)]*\)")

# --- 同位素 ---
_ISOTOPE_RES: Final = (
    # (d7) [D9] [13C6] (15N) ... —— 括号里只有同位素标记
    re.compile(r"[(\[]\s*(?:[dD][0-9]+|13C[0-9]*|15N[0-9]*|18O[0-9]*)\s*[)\]]"),
    # 氘：方括号 [2H7] [2H]（同位素修饰记法，个数可省）；圆括号 (2H7) 必须带个数——不带个数的 ``(2H)`` ``(12H)`` 是 IUPAC 的"指示氢"
    # （``benzofuran-3(2H)-one``、``6(12H)-one``），不是氘标记（植物次生代谢物系统名，原来被误报 isotope_label）
    re.compile(r"\[\s*[0-9]*2H[0-9]*\s*\]"),
    re.compile(r"\(\s*[0-9]*2H[0-9]+\s*\)"),
    # 后缀 -d7 / -D9（后面不能是数字或冒号：排除 iso-d17:1、Cer-d18:1；PNNL 的 "-d " 后是空格，也不匹配）
    re.compile(r"-[dD][0-9]+(?![0-9:])"),
    # 前缀 d7-PC ...
    re.compile(r"(?:^|[\s(\[/_])[dD][0-9]{1,2}-(?=[A-Za-z0-9])"),
    # 前缀标记 13C16-PC / 15N-PC / 2H9-PC：``2H`` 必须带个数——``-2H-pyran`` ``2H-chromen`` 的 ``2H`` 是指示氢
    re.compile(r"(?:^|[\s(\[/_-])(?:13C[0-9]*|15N[0-9]*|2H[0-9]+)-(?=[A-Za-z])"),
    # 独立的 d7 记号（前后是空白/串端）："PC 16:0 d7"、"d7 PC 16:0"
    re.compile(r"(?:^|\s)[dD][0-9]{1,2}(?=\s|$)"),
)

#: 维生素 D 名（``Vitamin D3``、``25-hydroxyvitamin D3``、``vitamin-D2``）：``D3`` 是维生素编号，不是氘标记。同位素检测前先遮掉
#: （``Vitamin D3-d3`` 遮掉后剩 ``-d3``，仍判同位素）。``[D9]`` ``(d9)`` ``-d7`` 这些带括号 / 连字符的真氘标记不受影响
_VITAMIN_D: Final = re.compile(r"(?i:vitamin)[\s-]*[dD][0-9]{1,2}(?![0-9:])")
#: 整个名字只有一个 ``D3``：没有任何脂质内容，是维生素缩写而不是氘标记（交给 M3 拒绝）
_BARE_D: Final = re.compile(r"[dD][0-9]{1,2}")

# --- 加合物 ---
_PNNL_TAIL: Final = re.compile(
    r"\s*;\s*(\[[^\[\]]*\]\s*[0-9]*[+-]?)(?:\s+((?:[A-Z][a-z]?[0-9]*)+)(?:\s+(M(?:[+-][A-Za-z0-9]+)*))?)?\s*$")
_BRACKET_TAIL: Final = re.compile(r"(?:\s*[;,])?\s*(" + _BRACKET_BODY + r")\s*$")
_SUFFIX_TAIL: Final = re.compile(r"(?<=[)\]])((?:[+-][0-9]*[A-Za-z][A-Za-z0-9]*)+)$")
_CODE_TAIL: Final = re.compile(r"\s(M(?:[+-][0-9]*[A-Za-z][A-Za-z0-9]*)+)$")
_NOPAREN_TAIL: Final = re.compile(r"(?<=[0-9])((?:[+-][0-9]*[A-Za-z][A-Za-z0-9]*)+)$")

# --- 其它尾部装饰 ---
_SN_TAIL: Final = re.compile(r"(?:-sn([123])|\[sn([123])\])$", re.IGNORECASE)
_STEREO_TAIL: Final = re.compile(r"\s*\[(iso[0-9]*|rac|R|S|U)\]$")
_VARIANT_TAIL: Final = re.compile(r"\s*\(([a-g])\)$")
_ANY_BRACKET_TAIL: Final = re.compile(r"\s*\[[^\[\]]*\]\s*[0-9]*[+-]?$")

# --- Ox 前缀 ---
_OX_PREFIX: Final = re.compile(r"(O2?x)(?=[A-Z])")
_OXY_MARK: Final = re.compile(r";\s*\(?[0-9]*(?:[Oo]|Ep|COOH)|[0-9]\([0-9]*(?:O|Ke|CHO|COOH)")

# --- LipidIN ---
_CHAIN_LABEL: Final = re.compile(r"chain\s*:\s*", re.IGNORECASE)
_O_ATOM: Final = re.compile(
    r";\s*(?:\(\s*([0-9]+)\s*OH\s*\)|([0-9]+)O(?![A-Za-z])|O([0-9]+)(?![A-Za-z0-9])|(OH|oxo|Ep|O)(?![A-Za-z0-9])"
    r"|(COOH|OOH|OO)(?![A-Za-z0-9]))")

#: ``Segment.notes`` 的键（冒号前部分）-> flag
NOTE_FLAGS: Final = MappingProxyType({
    "sn_suffix": "sn_suffix_stripped",
    "ox_prefix": "ox_prefix_ignored",
    "lipidin_sum_O_mismatch": "lipidin_sum_O_mismatch",
})


# ==========================================================================
# 同位素
# ==========================================================================
def _check_isotope(text: str) -> None:
    probe = text.strip()
    if _BARE_D.fullmatch(probe):
        return
    probe = _VITAMIN_D.sub(lambda m: "V" * len(m.group()), probe)
    for rx in _ISOTOPE_RES:
        m = rx.search(probe)
        if m is not None:
            raise Reject("isotope_label", f"isotope label {m.group(0).strip()!r} in {text!r}")


# ==========================================================================
# 加合物
# ==========================================================================
def _has_chain(text: str) -> bool:
    return _CDB.search(text) is not None


def _pop_adduct(text: str) -> tuple[str, str, tuple[str, ...]] | None:
    """认出并剥掉段尾**已识别**的加合物。返回 ``(剩余文本, 加合物原文, notes)``；没有则 ``None``。
    剥完为空时不剥。形似加合物但认不出的写法留给 ``_reject_adduct_like``。"""
    # 1) PNNL mgf：; [M-H]- C42H82N1O3 M-H
    m = _PNNL_TAIL.search(text)
    if m is not None:
        rest = text[:m.start()].rstrip()
        if rest and _info_bracket(m.group(1)) is not None:
            code = m.group(3)
            if code is None or code == re.sub(r"\s+", "", m.group(1))[1:].split("]")[0]:
                notes = (f"pnnl_formula:{m.group(2)}",) if m.group(2) else ()
                return rest, m.group(1).strip(), notes
    # 2) 括号式 [M+H]+
    m = _BRACKET_TAIL.search(text)
    if m is not None:
        rest = text[:m.start()].rstrip()
        if rest and _info_bracket(m.group(1)) is not None:
            return rest, m.group(1).strip(), ()
    # 3) LipidMatch 后缀 )+H
    m = _SUFFIX_TAIL.search(text)
    if m is not None and text[:m.start()] and _info_terms(m.group(1), "suffix") is not None:
        return text[:m.start()], m.group(1), ()
    # 3b) 裸后缀 "PC 16:0/18:1+H"：紧跟链末数字、没有括号。必须有链、剩下的非空、每一项都在物种表里（+O / +2O 这类不在表里 -> 不剥）
    m = _NOPAREN_TAIL.search(text)
    if (m is not None and text[:m.start()].strip() and _has_chain(text[:m.start()]) and _SN_TAIL.search(text[:m.start()]) is None
            and _info_terms(m.group(1), "suffix") is not None):                  # ``LPC(18:1)-SN1+H``：sn 后缀夹在链与加合物之间，不猜（同 ``)-SN1+H``）
        return text[:m.start()], m.group(1), ()
    # 4) 裸码 " M+H"
    m = _CODE_TAIL.search(text)
    if m is not None:
        rest = text[:m.start()].rstrip()
        if rest and _info_terms(m.group(1)[1:], "code") is not None:
            return rest, m.group(1), ()
    return None


def _reject_adduct_like(text: str) -> None:
    """段尾剥完后若还剩"形似加合物/尾部括号"的东西 -> ``Reject("unrecognized_tail")``（boundary R9）。"""
    m = _PNNL_TAIL.search(text)
    if m is not None and text[:m.start()].strip():
        raise Reject("unrecognized_tail", f"unrecognized PNNL adduct tail {m.group(0).strip()!r} in {text!r}")
    m = _BRACKET_TAIL.search(text)
    if m is not None and text[:m.start()].strip():
        raise Reject("unrecognized_tail", f"unrecognized adduct {m.group(1).strip()!r} in {text!r}")
    m = _SUFFIX_TAIL.search(text)
    if m is not None and _has_chain(text[:m.start()]):
        raise Reject("unrecognized_tail", f"unrecognized adduct suffix {m.group(1)!r} in {text!r}")
    m = _CODE_TAIL.search(text)
    if m is not None and text[:m.start()].strip():
        raise Reject("unrecognized_tail", f"unrecognized adduct code {m.group(1)!r} in {text!r}")
    m = _NOPAREN_TAIL.search(text)
    if m is not None and _has_chain(text[:m.start()]) and _info_terms(m.group(1), "suffix") is not None:
        raise Reject("unrecognized_tail", f"adduct {m.group(1)!r} not directly after ')' in {text!r}")
    m = _ANY_BRACKET_TAIL.search(text)
    if m is not None and _has_chain(text[:m.start()]):
        raise Reject("unrecognized_tail", f"unrecognized bracket tail {m.group(0).strip()!r} in {text!r}")


def _strip_tails(text: str) -> tuple[str, str, list[str]]:
    """反复剥段尾装饰，返回 ``(文本, 加合物, notes)``。"""
    adduct = ""
    notes: list[str] = []
    for _ in range(8):                       # 每轮至少缩短一次；上限只是保险
        r = _pop_adduct(text)
        if r is not None:
            if adduct:
                raise Reject("unrecognized_tail", f"more than one adduct in {text!r}")
            text, adduct, extra = r
            notes.extend(extra)
            continue
        m = _SN_TAIL.search(text)
        if m is not None and text[:m.start()].strip():
            notes.append(f"sn_suffix:SN{m.group(1) or m.group(2)}")
            text = text[:m.start()].rstrip()
            continue
        m = _STEREO_TAIL.search(text)
        if m is not None:
            rest = text[:m.start()]
            if rest[-1:] in ")0123456789" and rest.strip():
                notes.append(f"stereo:{m.group(1)}")
                text = rest.rstrip()
                continue
        m = _VARIANT_TAIL.search(text)
        if m is not None:
            rest = text[:m.start()]
            if rest[-1:] in ")0123456789" and rest.strip():
                notes.append(f"variant:{m.group(1)}")
                text = rest.rstrip()
                continue
        break
    _reject_adduct_like(text)
    return text, adduct, notes


# ==========================================================================
# Ox 前缀
# ==========================================================================
def _strip_ox_prefix(text: str) -> tuple[str, list[str]]:
    m = _OX_PREFIX.match(text)
    if m is None:
        return text, []
    rest = text[m.end():]
    if _OXY_MARK.search(rest) is None:        # 没有氧标记：前缀是唯一的氧化信息，不能丢
        return text, []
    return rest, [f"ox_prefix:{m.group(1)}"]


# ==========================================================================
# 段构造
# ==========================================================================
def _segment(text: str, kind: str) -> Segment:
    _check_isotope(text)
    text, adduct, notes = _strip_tails(text)
    text, ox = _strip_ox_prefix(text)
    return Segment(text=text, adduct=adduct, kind=kind, notes=tuple(notes + ox))


def _sum_cdb(text: str) -> tuple[int, int] | None:
    """文本里所有 ``C:DB`` 的总和（不计 ``(FA c:d)``）；一个也没有返回 ``None``。"""
    pairs = _CDB.findall(_NESTED_FA.sub("", text))
    if not pairs:
        return None
    return sum(int(c) for c, _ in pairs), sum(int(d) for _, d in pairs)


def _o_count(text: str) -> int:
    """名字里写出的氧个数的粗估（仅用于 ``lipidin_sum_O_mismatch`` 这个提示 flag）。"""
    total = 0
    for m in _O_ATOM.finditer(_NESTED_FA.sub("", text)):
        if m.group(1):
            total += 1
        elif m.group(2):
            total += int(m.group(2))
        elif m.group(3):
            total += int(m.group(3))
        elif m.group(4):
            total += 1
        else:
            total += 2
    return total


def _lipidin(title: str, detail: str) -> Segment:
    """LipidIN ``标题 | Chain: 详情``：``detail`` 已去掉 ``Chain:`` 标签。"""
    if not detail:
        raise Reject("unrecognized_tail", f"empty 'Chain:' segment after {title!r}")
    _check_isotope(title)
    seg = _segment(detail, "lipidin_chain")
    sx, sy = _sum_cdb(title), _sum_cdb(seg.text)
    if sx is not None and sy is not None and sx != sy:
        raise Reject("sum_detail_conflict",
                     f"title C:DB {sx[0]}:{sx[1]} != chain C:DB {sy[0]}:{sy[1]} ({title!r} | {detail!r})")
    notes = seg.notes
    if _o_count(title) != _o_count(seg.text):
        notes = notes + ("lipidin_sum_O_mismatch",)
    return Segment(text=seg.text, adduct=seg.adduct, kind="lipidin_chain", notes=notes)


# ==========================================================================
# 入口
# ==========================================================================
def split_decorations(s: str, cfg: Config) -> list[Segment]:
    """见模块文档。``cfg`` 目前不影响 M2（保留以便与 M1/M3 同签名）。

    >>> from lipid_name_map.config import DEFAULT
    >>> [(g.kind, g.text, g.adduct) for g in split_decorations("PC(16:0_18:1)+H", DEFAULT)]
    [('plain', 'PC(16:0_18:1)', '+H')]
    >>> [g.text for g in split_decorations("Cer 26:7 | Chain: Cer 8:5;2O/18:2", DEFAULT)]
    ['Cer 8:5;2O/18:2']
    """
    if not isinstance(s, str):
        raise Reject("not_string", f"expected str, got {type(s).__name__}")
    text = s.strip()
    if not text:
        raise Reject("empty", "input is empty or whitespace only")
    if "|" not in text:
        return [_segment(text, "plain")]
    parts = [p.strip() for p in text.split("|")]
    if not all(parts):
        raise Reject("unrecognized_tail", f"empty alias segment around '|' in {text!r}")
    labels = [_CHAIN_LABEL.match(p) for p in parts]
    if any(labels):
        label = labels[1] if len(parts) == 2 and labels[0] is None else None
        if label is None:
            raise Reject("unrecognized_tail", f"unexpected 'Chain:' layout in {text!r}")
        return [_lipidin(parts[0], parts[1][label.end():].strip())]
    return [_segment(p, "alias") for p in parts]


def segment_flags(seg: Segment) -> tuple[str, ...]:
    """由 ``Segment`` 推出装饰类 flags（排序去重）；见模块文档。"""
    out: set[str] = set()
    if seg.kind == "alias":
        out.add("alias_split")
    elif seg.kind == "lipidin_chain":
        out.add("lipidin_chain_used")
    if seg.adduct:
        out.add("adduct_stripped")
    for note in seg.notes:
        flag = NOTE_FLAGS.get(note.split(":", 1)[0])
        if flag is not None:
            out.add(flag)
    return tuple(sorted(out))
