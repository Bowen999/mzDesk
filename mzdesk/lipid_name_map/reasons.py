"""原因码注册表（SYNTHESIS §1.4，W0 冻结；新增原因码须架构师批准）。

对外接口
--------
``REASONS``            ``dict[str, Reason]``；``Reason`` 是 ``(status, module, text)`` 三元 NamedTuple
                       - status: ``"rejected"`` | ``"unsupported"``
                       - module: 抛出该码的模块（``"M3"``；多个抛出方用 ``/`` 连接，如 ``"M6/M4b"``）
                       - text:   白话描述（给人读，也进文档）
``status_of(code)``    原因码 -> status（未知码抛 ``ValueError``）
``module_of(code)`` ``describe(code)``
``ALL_CODES`` ``REJECTED_CODES`` ``UNSUPPORTED_CODES``    frozenset
``STATUS_OK`` ``STATUS_REJECTED`` ``STATUS_UNSUPPORTED``  状态字符串常量
``check_integrity()``  注册表自检，返回问题清单（空 = 通过）

规则
----
* ok 行的 ``reason`` 恒为 ``""``，不在本表里。
* 本模块不 import 包内任何其他模块（被 types / contract 依赖，避免循环）。
* 调用方只给 code：``raise Reject("chain_count", detail)``；status 由本表查，不得在调用处自拟。
"""
from __future__ import annotations

import re
from typing import Final, NamedTuple

STATUS_OK: Final = "ok"
STATUS_REJECTED: Final = "rejected"
STATUS_UNSUPPORTED: Final = "unsupported"
STATUSES: Final = (STATUS_OK, STATUS_REJECTED, STATUS_UNSUPPORTED)

#: 允许出现在 ``Reason.module`` 里的模块记号（SYNTHESIS §2.2；粗粒度 M1..M9 与细粒度子模块都允许）
MODULES: Final = frozenset({
    "M1", "M2", "M3", "M4", "M4a", "M4b", "M5", "M5a", "M5b",
    "M6", "M6c", "M6f", "M7", "M8", "M9",
})


class Reason(NamedTuple):
    status: str
    module: str
    text: str


def _block(status: str, rows: list[tuple[str, str, str]]) -> dict[str, Reason]:
    return {code: Reason(status, module, text) for code, module, text in rows}


# --------------------------------------------------------------------------
# rejected：不是可接受的 shorthand（俗名、乱码、链数不符…）
# --------------------------------------------------------------------------
_REJECTED = _block(STATUS_REJECTED, [
    # M1 预清洗
    ("not_string", "M1", "输入不是字符串（None、数字、bytes、列表等）"),
    ("empty", "M1", "输入为空串或只有空白"),
    ("too_long", "M1", "输入长度超过 Config.max_len"),
    # M3 common-name 拒绝
    ("id_or_number", "M3", "像数据库编号或纯数字（如 LMGP01010001、CHEBI:12345），不是 shorthand"),
    ("class_only", "M3", "只有类名、没有链信息（如 PC、GM1、SPB）"),
    ("no_chain", "M3", "没有 C:DB 链信息，也不是类名（如 Stearamide、系统名）"),
    ("common_name", "M3", "俗名/长式类名/文字链/系统名片段（如 Phosphatidylcholine 16:0/18:1）"),
    ("bare_chain", "M3", "只有链、没有类名（如 16:0_16:0）"),
    ("bad_chars", "M3", "含有 shorthand 里不可能出现的字符"),
    # M2 装饰/加合物/别名拆分
    ("unrecognized_tail", "M2", "名字末尾有无法识别的装饰（既不是加合物，也不是已知后缀）"),
    ("sum_detail_conflict", "M2", "LipidIN 'X | Chain: Y' 中和式与链段的 C:DB 互相矛盾"),
    # M9 段合并
    ("alias_conflict", "M9", "'A | B' 别名各段的类或 L1 互相矛盾，无法选出一致的段"),
    # M5 词法/语法
    ("glued_tokens", "M5", "链记号粘连（如 16:017:0）"),
    ("dangling_sep", "M5", "分隔符悬空（开头、结尾或连续分隔符）"),
    ("unbalanced_parens", "M5", "括号不配对"),
    ("db_ge_c", "M5", "双键数不小于碳数（C>=2 时 DB>=C）"),
    ("dup_atom", "M5", "同一条链上同种氧/修饰记号重复出现"),
    ("multi_ether", "M5", "出现两条以上醚链（或一条链上多个醚标记）"),
    ("all_zero_chains", "M5", "所有链都是 0:0 占位，没有实际链"),
    ("class_in_chain", "M5", "链位置里混入了类名"),
    ("malformed", "M5", "整体语法不成形，无法切分头部与链体"),
    # M6 语义
    ("chain_count", "M6", "链数与该类要求的链数不符（且不能解释为占位/和式）"),
    ("missing_oxygen_count", "M6", "鞘脂名字缺少氧数（如 SM 30:1、SPB 18:1）；默认按 d 型补 ;O2，只在 Config.assume_base_oxygen=False 时拒绝"),
    ("single_token_base", "M6", "只有一个 d/t/m 记号且 C<26，分不清碱基链还是和式（如 Cer d18:1）"),
    ("class_chain_conflict", "M6/M4b", "类名携带的信息与链上的标记互相矛盾（如醚种类不一致、无唯一类行）"),
    ("ether_not_allowed", "M6", "该类不允许醚链，却出现了醚标记"),
    ("o_not_expressible", "M6", "氧总数在 L1 词表里无法表达（O12 或 >O13）"),
    # M9 兜底
    ("internal_error", "M9", "内部异常（永远是 bug，检查器按失败处理）"),
])

# --------------------------------------------------------------------------
# unsupported：语法像 shorthand，但本版本有意不处理
# --------------------------------------------------------------------------
_UNSUPPORTED = _block(STATUS_UNSUPPORTED, [
    # M4 类解析
    ("class_unknown", "M4a", "类名不在别名表里"),
    ("class_not_in_vocab", "M4a/M4b", "类名已知，但对应的类不在 v0.2.11 词表（如 BA、TG-OAHFA、SQDG-O、总量写法的 TG-EST；CoQ 只收 CoQn 与原生 CoQ C:DB;O4）"),
    ("class_alias_unverified", "M4a", "别名未经验证，默认不收（Config.allow_unverified_alias=False）"),
    ("polymer", "M4a", "聚合物/非脂质（PEG、LipidA-PP 等）"),
    ("glycan_unknown", "M4a", "糖链组成不在词表里（组成式生成后校验失败）"),
    # M2
    ("isotope_label", "M2", "同位素标记（-d7、[D9] 等），本版本不处理"),
    # M5/M6
    ("unsupported_subclass", "M6", "不支持的亚类（Cer 的 EOS/EODS/EBDS/HS/HDS 等）"),
    ("unsupported_syntax", "M5", "语法像 shorthand 但本版本不处理（炔位置、dO-、p 后缀等）"),
    ("unsupported_chain_modifier", "M5", "链修饰记号本版本不处理（如 (Ke)(CHO)(COOH)(KeKe)）"),
    ("oxylipin_disabled", "M6", "介质缩写（HETE、EET 等）被 Config.allow_oxylipin_abbrev 关闭"),
])

REASONS: Final[dict[str, Reason]] = {**_REJECTED, **_UNSUPPORTED}

ALL_CODES: Final = frozenset(REASONS)
REJECTED_CODES: Final = frozenset(_REJECTED)
UNSUPPORTED_CODES: Final = frozenset(_UNSUPPORTED)


def status_of(code: str) -> str:
    """原因码 -> ``"rejected"`` / ``"unsupported"``；未知码抛 ``ValueError``。"""
    try:
        return REASONS[code].status
    except (KeyError, TypeError):
        raise ValueError(f"unknown reason code {code!r}") from None


def module_of(code: str) -> str:
    try:
        return REASONS[code].module
    except (KeyError, TypeError):
        raise ValueError(f"unknown reason code {code!r}") from None


def describe(code: str) -> str:
    """白话描述。"""
    try:
        return REASONS[code].text
    except (KeyError, TypeError):
        raise ValueError(f"unknown reason code {code!r}") from None


#: SYNTHESIS §1.4 的规模（rejected 28 + unsupported 10）。改动须同步架构师批准。
EXPECTED_COUNTS: Final = {STATUS_REJECTED: 28, STATUS_UNSUPPORTED: 10}

_CODE_RE = re.compile(r"[a-z][a-z0-9_]*")


def check_integrity() -> list[str]:
    """注册表自检；返回问题描述列表，空列表 = 通过。

    检查：码名合法且唯一、status 合法、module 记号合法、描述非空、
    ok 不在表里、两个 status 的集合不相交且并集 == 全集、各 status 数量等于冻结规模。
    """
    problems: list[str] = []
    for code, r in REASONS.items():
        if not _CODE_RE.fullmatch(code):
            problems.append(f"{code!r}: illegal code name")
        if r.status not in (STATUS_REJECTED, STATUS_UNSUPPORTED):
            problems.append(f"{code}: bad status {r.status!r}")
        mods = r.module.split("/")
        if not mods or any(m not in MODULES for m in mods):
            problems.append(f"{code}: bad module {r.module!r}")
        if not r.text.strip():
            problems.append(f"{code}: empty description")
    if "ok" in REASONS or "" in REASONS:
        problems.append("'ok' / '' must not be reason codes")
    if REJECTED_CODES & UNSUPPORTED_CODES:
        problems.append(f"codes in both statuses: {sorted(REJECTED_CODES & UNSUPPORTED_CODES)}")
    if REJECTED_CODES | UNSUPPORTED_CODES != ALL_CODES:
        problems.append("status sets do not cover REASONS")
    if len(REASONS) != len(_REJECTED) + len(_UNSUPPORTED):
        problems.append("duplicate code across status blocks")
    for status, n in EXPECTED_COUNTS.items():
        got = sum(1 for r in REASONS.values() if r.status == status)
        if got != n:
            problems.append(f"{status}: expected {n} codes, found {got}")
    if REASONS.get("internal_error", Reason("", "", "")).status != STATUS_REJECTED:
        problems.append("internal_error must be 'rejected'")
    return problems
