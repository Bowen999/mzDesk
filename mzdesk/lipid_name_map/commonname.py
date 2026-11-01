"""M3 common-name / non-shorthand rejection（SYNTHESIS §2.2 M3；依据 ref_census §5 的 ``is_shorthand_like`` 原型）。

这是"拒绝俗名"的第一道门：**纯语法层**，不认识类表（类名认不认识是 M4a 的事），宁可放过、不可误拒——
放过的东西下游还会被类表/词法/语义拒掉，而误拒会直接丢掉合法名字。

对外接口
--------
``check_shorthand(seg, cfg=DEFAULT) -> None``
    ``seg`` 是 M2 产出的 ``Segment``（也接受裸 ``str``，方便工具/测试）。合格返回 ``None``；
    不合格抛 ``Reject(code, detail)``，``code`` 只会是 ``REASON_CODES`` 之一::

        id_or_number  数据库编号/纯数字（LMGP01010001、SLM:…、CHEBI:…、InChIKey、CAS、760.585）
        class_only    只有类名、没有 C:DB（PC、CoQ9、GM1、Cer-NS）
        no_chain      没有 C:DB 且不是类名（Stearamide、Cholesterol、palmitic acid、系统名、PEG n=…、``PC aa C36:2``）
        common_name   有链，但类位置/链位置是俗名：长式类名（Phosphatidylcholine 16:0/18:1）、文字链、系统名片段、
                      SPLASH_* 内标名；单独出现的介质缩写（12-HETE）；``allow_oxylipin_abbrev=False`` 时链位里的介质缩写
        bare_chain    有链但无类（16:0_16:0、18:1;O2/22:0;O、(16:0/18:1)、;O2/18:0）
        bad_chars     含 shorthand 里不可能出现的字符（非 ASCII、``{ } " # % & ? !`` …，残留的 ``|``）
``diagnose(text, cfg=DEFAULT) -> tuple[str, str] | None``
    不抛异常版：返回 ``(code, detail)`` 或 ``None``（供统计/测试用）。
``is_shorthand_like(text, cfg=DEFAULT) -> tuple[bool, str]``
    与普查原型 ``(bool, reason)`` 同形：``(True, "ok")`` 或 ``(False, code)``。
``REASON_CODES``
    本模块可能抛出的 6 个原因码（均在 ``reasons.REASONS`` 注册，status=rejected）。

判定顺序（每个 `` | `` 分段各自判定，任一段不过即整体不过；段内先剥 ``Chain:`` 前缀与尾部加合物）
-------------------------------------------------------------------------------------------
1. ID/数字                                         → id_or_number
2. 字符集（``A-Za-z0-9 ()[],.'+-/;:_=*`` 与空格之外；``< >`` 只在 IUPAC 糖链 ``(2->3)`` 里出现）→ bad_chars
3. ``reject_prefixes.csv`` 的前缀（SPLASH）           → 表中给定的码
4. 找第一个 ``C:DB``：其前不得紧贴字母/数字（``C18:1``、``PC34:1``、``PC aa C36:2`` 不是链）
   - 没有链（含 ``a:b`` 紧跟在 ``数字,`` 后的 IUPAC 位置号列表）：单独的介质缩写 → common_name；整串是"类名形"（≤3 个词、无长式词、括号配对、不是 IUPAC 系统名）→ class_only；否则 no_chain
   - 有链但类位置为空/无字母/不以字母数字或括号起头（``;O2/18:0`` 这类残片）  → bare_chain
   - 类位置含 ≥6 个连续小写字母（Phosphatidylcholine）或 ≥9 个连续大写字母、或多词（Sterol ester）→ common_name
   - 链位置含 ≥6 字母的文字词（palmitoyl）或系统名片段（glycero/phospho/sn-）→ common_name
   - 链位置含介质缩写且 ``cfg.allow_oxylipin_abbrev`` 为假                        → common_name
白名单（``data/commonname/allow_words.csv``，整词、区分大小写）：Forssman Branched Plasmenyl Plasmanyl plasmenyl plasmanyl anteiso methyl，
再加旧 v26 类名里的 plasmalogen（``PC-plasmalogen``，§3.2 第 4 步的别名行要靠它才到得了 M4a）。拒绝前缀表 ``reject_prefixes.csv``：SPLASH。
注意：糖链 IUPAC 记法（``Galbeta-Cer(d18:1/16:0)``、``NeuAcalpha2-3Galbeta1-4Glcbeta-Cer(…)``）是 LIPID MAPS 的系统名写法，
含 ≥6 连续小写字母，按验收（SYSTEMATIC_NAME 0% 通过）在这里以 common_name 拒掉，到不了 M4a 的 ``glycan_unknown``。

不在本模块范围（放过，交下游）：未知类名（LysoPC、TAG、PtdCho、pc）→ M4a；链数/醚/氧 → M5/M6；同位素/加合物 → M2。
"""
from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Final

from .config import DEFAULT, Config
from .types import OXYLIPIN_ABBREVS, Reject, Segment

__all__ = ["REASON_CODES", "check_shorthand", "diagnose", "is_shorthand_like"]

#: 本模块可能抛出的原因码（SYNTHESIS §1.4 的 M3 行）
REASON_CODES: Final = ("id_or_number", "class_only", "no_chain", "common_name", "bare_chain", "bad_chars")

_DATA_DIR: Final = Path(__file__).resolve().parent / "data" / "commonname"
_DETAIL_MAX: Final = 200       # Reject.detail 的最大长度（输入最长 1000 字符，别把整串回显进诊断）


# --------------------------------------------------------------------------
# 数据文件（白名单词、拒绝前缀）
# --------------------------------------------------------------------------
def _read_csv(name: str) -> list[dict[str, str]]:
    with open(_DATA_DIR / name, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(ln for ln in fh if ln.strip() and not ln.lstrip().startswith("#")))


def _word_mask(words: list[str]) -> re.Pattern[str] | None:
    """整词（不被字母夹住；其后不得再接小写字母）匹配这些白名单词。"""
    if not words:
        return None
    alt = "|".join(re.escape(w) for w in sorted(set(words), key=len, reverse=True))
    return re.compile(rf"(?<![A-Za-z])(?:{alt})(?![a-z])")


def _load_allow() -> tuple[re.Pattern[str] | None, re.Pattern[str] | None]:
    head: list[str] = []
    tail: list[str] = []
    for r in _read_csv("allow_words.csv"):
        scope = r["scope"].strip()
        if scope not in ("head", "tail", "both"):
            raise ValueError(f"allow_words.csv: bad scope {scope!r} for {r['word']!r}")
        if scope in ("head", "both"):
            head.append(r["word"].strip())
        if scope in ("tail", "both"):
            tail.append(r["word"].strip())
    return _word_mask(head), _word_mask(tail)


def _load_reject_prefixes() -> tuple[tuple[str, str], ...]:
    rows = []
    for r in _read_csv("reject_prefixes.csv"):
        code = r["code"].strip()
        if code not in REASON_CODES:
            raise ValueError(f"reject_prefixes.csv: code {code!r} is not an M3 reason code")
        rows.append((r["prefix"].strip().casefold(), code))
    return tuple(rows)


_ALLOW_HEAD, _ALLOW_TAIL = _load_allow()
_REJECT_PREFIXES: Final = _load_reject_prefixes()


# --------------------------------------------------------------------------
# 正则
# --------------------------------------------------------------------------
# 1. ID / 纯数字
_ID_RE: Final = re.compile(
    r"""^(?:
        LM[A-Z]{2}\d{6,}[A-Za-z0-9]*                 # LIPID MAPS
      | SLM:\d+                                      # SwissLipids
      | (?i:chebi):\d+ | HMDB\d+ | [CDGKR]\d{5}      # ChEBI / HMDB / KEGG
      | [A-Z]{14}-[A-Z]{10}-[A-Z]                    # InChIKey
      | InChI(?:Key)?=\S+                            # InChI / InChIKey=…
      | \d{2,7}-\d{2}-\d                             # CAS
      | (?:CID|PubChem|KEGG|LIPIDBANK|CAS)[:_\s-]*\d{3,}
      | [+-]?\d+(?:[.,]\d+)*                         # 纯数字 / m/z
    )$""",
    re.X,
)

# 2. 字符集（shorthand 里可能出现的全部字符；``|`` 不在内——分段后残留的 ``|`` 算 bad_chars）
_BAD_CHAR_RE: Final = re.compile(r"[^A-Za-z0-9()\[\],.'+\-/ ;:_=*]")

# 剥离：``Chain:`` 前缀、尾部加合物（M2 通常已剥；这里只为独立使用时判得准）。加合物只认元素式（H Na NH4 HCO2 C2H3O2…），
# 不能放宽成"任意大写开头词"：否则 ``5(S)-HETE`` 的 ``-HETE`` 会被当加合物剥掉，绕过介质缩写开关。
_CHAIN_LABEL_RE: Final = re.compile(r"^Chain:\s*")
_ADDUCT_ELEM: Final = r"(?:Na|Li|Cl|Br|Mg|Ca|Fe|Cu|Zn|Ag|H|C|N|O|K|F|I|S|P|D)\d*"
_ADDUCT_TAIL_RE: Final = re.compile(rf"\s*(?:\[M[^\]]*\]\d*[+-]*|(?<=\))(?:[+-]\d*(?:{_ADDUCT_ELEM})+)+)\s*$")
_PIPE_SPLIT_RE: Final = re.compile(r"\s+\|\s+")

# 3. 链
_CORE_RE: Final = re.compile(r"\d+:\d+")
_CHAIN_PREFIX: Final = r"(?:[dte]?[OP]-|N-|(?:anteiso|iso|ai)-?|methyl-?|[dtem](?=\d))"
# 紧跟在 ``数字,`` 后面的 ``a:b`` 是 IUPAC 位置号列表（环氧桥 ``5,6:5',6'-diepoxy``、``16,24:20,24-diepoxy``），不是链记号
_CHAIN_START_RE: Final = re.compile(rf"(?<![A-Za-z0-9])(?<![0-9],){_CHAIN_PREFIX}?\d+:\d+")

_HAS_LETTER_RE: Final = re.compile(r"[A-Za-z]")
_HEAD_START_RE: Final = re.compile(r"[A-Za-z0-9(\[]")                      # 类位置的合法首字符（参考表 357 个类头都满足）
_LOWER_RUN_RE: Final = re.compile(r"[a-z]{6,}")                         # ≥6 个连续小写字母 = 长式词
_UPPER_RUN_RE: Final = re.compile(r"[A-Z]{9,}")                         # ≥9 个连续大写字母 = 全大写长式名
_HEAD_EXTRA_WORD_RE: Final = re.compile(r"(?:alpha|beta|gamma)?(?:\([^)]*\))?")    # 类名后允许的附加词
_CLASS_TOKEN_CHARS_RE: Final = re.compile(r"[A-Za-z0-9()\[\]+,.'\-_ ]+")

# IUPAC 系统名（``Octan-2-one`` ``(2E)-undec-2-en-1-ol`` ``Hexane-2,3-dione`` ``pregn-4-ene-3,20-dione`` ``Decane``）不是类缩写：
#   ① 位置号 + 连字符 + 带 ≥2 个小写字母的词（``2-one`` ``22(29)-ene`` ``2E-Hexene``；位置号不在括号里起头，``Hex(4)-HexNAc-Cer`` 的 ``(4)-`` 不算，
#      ``AcCar(14-OH)`` 的 ``-OH`` 没有小写字母也不算）  ② 位置号列表 ``3beta,6alpha``  ③ 烷 / 烯 / 炔 / 醇 / 酮 / 醛 词尾（``Decane`` ``Nonene``）。
# 类缩写（``PC`` ``CoQ9`` ``GM1`` ``Cer-NS`` ``AcCar(4-OH)`` ``Acetyl-CoA`` ``Hex(4)-HexNAc(3)-Fuc-NeuAc-Cer``）一个都不命中；
# 命中的名字没有链且不是类名 -> no_chain（reasons.py：``no_chain`` 的例子就是"系统名"），不是 class_only。
_SYSTEMATIC_RE: Final = re.compile(
    r"(?<![(A-Za-z0-9])\d+(?:\(\d+\))?[A-Za-z]{0,5}\)?-[A-Za-z]*[a-z]{2}"
    r"|\d[A-Za-z]{0,5},\d"
    r"|(?<![A-Za-z])[A-Za-z]*[a-z]{2}(?:an|en|yn)(?:e|ol|one|al|diol|triol|dione|olide)(?![a-z])"
)

_LETTER_RUN_RE: Final = re.compile(r"[A-Za-z]+")
_FRAGMENT_RE: Final = re.compile(r"oyl|yl\)|glycero|phospho|sn-")        # 系统名片段（enoyl 已含在 oyl 内）
# LipidMatch 的氧括号：``(O)`` ``(OOH)`` ``(OOOOOO)`` ``(OOOOOOOOH)`` 以及复合修饰 ``(KeOKeO)`` ``(OOOKeO)``——括号里只由 ``O`` / ``H`` / ``Ke``
# 拼成，是氧修饰记号（词法 ``chains._O_RUN`` 接受任意长度的 ``O+H?``，``Ke`` 之类走 ``unsupported_chain_modifier``），不是文字词；
# 长度不能成为拒绝理由，所以先把这种括号整体屏蔽再做"≥6 字母文字词"检测。
_OXY_BRACKET_RE: Final = re.compile(r"\(H*(?:O|Ke)(?:O|H|Ke)*\)")

# 4. 介质缩写（单独出现 / 链位）
_OXY_ALT: Final = "|".join(sorted((re.escape(a) for a in OXYLIPIN_ABBREVS), key=len, reverse=True))
_OXY_WORD_RE: Final = re.compile(rf"(?<![A-Za-z])(?:{_OXY_ALT})(?![A-Za-z])")
_OXY_POS: Final = r"\d+[RSEZ]?(?:\([RSEZ]\))?"
_OXY_ALONE_RE: Final = re.compile(rf"(?:{_OXY_POS}(?:\s*,\s*{_OXY_POS})*\s*-\s*)?(?:{_OXY_ALT})")


# --------------------------------------------------------------------------
# 判定
# --------------------------------------------------------------------------
def _mask(text: str, pat: re.Pattern[str] | None) -> str:
    return pat.sub("", text) if pat is not None else text


def _head_problem(head: str) -> str | None:
    """类位置（链之前的部分）是否像长式/多词名；合格返回 None。"""
    h = _mask(head, _ALLOW_HEAD)
    for w in _LETTER_RUN_RE.findall(h):
        if len(w) >= 6 and _LOWER_RUN_RE.search(w):
            return f"long-form word {w!r} in class position"
        if len(w) >= 9 and _UPPER_RUN_RE.search(w):
            return f"long upper-case word {w!r} in class position"
    words = head.split()
    if len(words) == 2 and words[1].startswith(";") and words[1].endswith(("/", "_")):
        # ``SM ;O2/22:0``：LipidDetective 丢了第一条链，类名后面只剩 ``;O2`` 残片和分隔符。残片是链残留，不是第二个类词；
        # 这是悬空分隔符，交词法报 ``dangling_sep``（长式词的判定在上面，照旧；残片多于一个词的仍按多词类位置拒绝）
        return None
    if len(words) > 3 or any(not _HEAD_EXTRA_WORD_RE.fullmatch(w) for w in words[1:]):
        return f"multi-word class position {head!r}"
    return None


def _tail_problem(tail: str, cfg: Config) -> str | None:
    """链位置（第一个 C:DB 起的部分）是否含文字链/系统名片段/被关掉的介质缩写；合格返回 None。"""
    m = _OXY_WORD_RE.search(tail)
    if m:
        if not cfg.allow_oxylipin_abbrev:
            return f"oxylipin abbreviation {m.group()!r} in chain position (Config.allow_oxylipin_abbrev=False)"
        tail = _OXY_WORD_RE.sub("", tail)
    t = _OXY_BRACKET_RE.sub("()", _mask(tail, _ALLOW_TAIL))
    for w in _LETTER_RUN_RE.findall(t):
        if len(w) >= 6:
            return f"text word {w!r} in chain part"
    m = _FRAGMENT_RE.search(t)
    if m:
        return f"systematic-name fragment {m.group()!r} in chain part"
    return None


def _is_class_token(s: str) -> bool:
    """整串是否"类名形"（无链时用来区分 class_only 与 no_chain）：≤3 词、字符集受限、括号配对、无长式词、不是系统名。"""
    return bool(
        _HAS_LETTER_RE.search(s)
        and s.count("(") == s.count(")")
        and s.count("[") == s.count("]")
        and _CLASS_TOKEN_CHARS_RE.fullmatch(s)
        and _head_problem(s) is None
        and _SYSTEMATIC_RE.search(s) is None
    )


def _diagnose_part(part: str, cfg: Config) -> tuple[str, str] | None:
    s = _ADDUCT_TAIL_RE.sub("", _CHAIN_LABEL_RE.sub("", part.strip())).strip()
    if not s:
        return "no_chain", "empty segment"
    if _ID_RE.match(s):
        return "id_or_number", f"{s!r} looks like a database id or a number"
    bad = _BAD_CHAR_RE.search(s)
    if bad:
        return "bad_chars", f"character {bad.group()!r} cannot occur in a shorthand"
    low = s.casefold()
    for prefix, code in _REJECT_PREFIXES:
        if low.startswith(prefix):
            return code, f"{prefix!r}-prefixed product/mixture name"

    core = _CORE_RE.search(s)
    chain = _CHAIN_START_RE.search(s)
    if core is None or chain is None or core.start() < chain.start():
        # 没有（合法起头的）链
        if core is not None:
            return "no_chain", f"{core.group()!r} is glued to a preceding letter/digit, so it is not a chain token"
        if _OXY_ALONE_RE.fullmatch(s):
            return "common_name", f"{s!r} is a free oxylipin mediator name, not a lipid shorthand"
        if _is_class_token(s):
            return "class_only", f"{s!r} has a class-like token but no C:DB chain"
        return "no_chain", f"no C:DB chain and not a class abbreviation: {s!r}"

    head = s[: chain.start()].rstrip()
    if head.endswith("("):
        head = head[:-1].rstrip()
    if not _HAS_LETTER_RE.search(head) or not _HEAD_START_RE.match(head):
        return "bare_chain", f"chain without a class: {s!r}"      # 含 ``;O2/18:0`` 这类只剩氧记号的残片
    problem = _head_problem(head)
    if problem:
        return "common_name", problem
    problem = _tail_problem(s[chain.start():], cfg)
    if problem:
        return "common_name", problem
    return None


def diagnose(text: str, cfg: Config = DEFAULT) -> tuple[str, str] | None:
    """不抛异常版：``(code, detail)`` 或 ``None``（合格）。多段（`` | ``）时第一个不合格段决定结果。"""
    if not isinstance(text, str):
        return "bad_chars", f"non-string input of type {type(text).__name__}"
    parts = [p for p in _PIPE_SPLIT_RE.split(text.strip()) if p.strip()]
    if not parts:
        return "no_chain", "empty segment"
    for i, p in enumerate(parts):
        res = _diagnose_part(p, cfg)
        if res is not None:
            code, detail = res
            if len(parts) > 1:
                detail = f"alias part {i + 1}: {detail}"
            return code, (detail if len(detail) <= _DETAIL_MAX else detail[: _DETAIL_MAX - 1] + "…")
    return None


def check_shorthand(seg: Segment | str, cfg: Config = DEFAULT) -> None:
    """M3 入口：合格返回 ``None``，否则抛 ``Reject``（码见模块文档）。"""
    res = diagnose(seg.text if isinstance(seg, Segment) else seg, cfg)
    if res is not None:
        raise Reject(res[0], res[1])


def is_shorthand_like(text: str, cfg: Config = DEFAULT) -> tuple[bool, str]:
    """普查原型同形接口：``(True, "ok")`` 或 ``(False, 原因码)``。"""
    res = diagnose(text, cfg)
    return (True, "ok") if res is None else (False, res[0])
