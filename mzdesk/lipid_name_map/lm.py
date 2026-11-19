"""M10：LIPID MAPS 名 ``LM_species`` / ``LM_molecular``（0.3.0，只加字段，不碰 L1/L2/L3）。

两个字段
--------
``LM_species``     LIPID MAPS 数据库 ABBREVIATION 风格的**总量级**名：``<lm_head> [O-]<C>:<DB>[;On]``，如 ``PC 34:1``、``LPC 18:1``、
                   ``PE O-34:2``（缩醛 ``P-`` 在总量级也写 ``O-``，DB 取 L1 口径，已含 +1）、``Cer 42:1;O2``、``ST 24:1;O4;Gly``。
                   数值取 ``Labels.sum_meta``（= L1），按 ``lm_heads.csv`` 的 ``@adj`` / ``@suffix`` 调整。
``LM_molecular``   数据库类名 + LIPID MAPS 2020 链写法的**链级**名：``PC 16:0_18:1``、``PE P-18:0_16:1``、``Cer 18:1;O2/24:0``、
                   ``TG 16:0_18:1_18:2``、``SE 27:1/18:1``。只有下面全部成立才出，否则 ``""`` 并在 ``meta.lm.molecular`` 记原因码。

两个字段都不从 L2/L3 字符串反解析：输入是内部 ``Interp``（``alts[*].chains``）+ ``Labels.sum_meta`` + 输入的类头文本 + 分隔符是否有序（``api`` 传入）。
固定结构类 ``CoQ`` 没有 ``Interp``，按 ``lm_heads.csv`` 的 ``@render`` 写（两个字段相同）。

数据：``data/lm_heads.csv``（``class_full -> lm_head``，另有 ``source`` ``fmt`` 与 ``note`` 里的机器记号；读法见 ``docs/LM_HEADS.md``）。
本模块认的记号：``@adj=C-27,DB-1,O-1``（LM 值 = L1 值 + 偏移）、``@suffix=;Gly``、``@render=CoQ{n}``、``@blank=<原因码>``（0.3.0 新增：
该类两个字段都留空，原因码进 ``meta.lm``——用于默认配置下读不回的类）。``@alt_head`` 一律不用（只有胆固醇酯 CE/SE 见下）。

``source``：``lmsd``（ABBREVIATION 里有）/ ``inferred``（数据库没有，用推断的类名；出名时打 flag ``lm_head_inferred``）/ ``none``（两个字段都留空）。

LM_molecular 的条件（逐条对每个假设 ``alts[i]`` 检查，**所有假设渲染出同一个字符串**才出；否则 ``hypotheses_differ``）
----------------------------------------------------------------------------------------------------------------
a. ``chains`` 非 None（否则 ``sum_only``）且链数 ∈ 类的名义链数（``chain_count``）；
b. 所有氧都在链上（``Struct.tail_o == 0``，否则 ``oxygen_off_chain``）；唯一例外：单链类，全部氧算在那条链上；
c. 醚类恰好一条醚链（类型 O 或 P 已确定），否则 ``ether_chain_unknown``；非醚类没有醚链；
d. ``fmt=sphingo``：恰好 2 链，其中能确定哪条是碱基，否则 ``sphingo_base_unknown``。按顺序认三种依据：
     ① 链的 ``role == "base"``（输入有 ``d/t/m`` 前缀）；
     ② 没有任何 ``role`` 但输入是有序写法——分隔符全是 ``/``——且类的 ``slash_role == "base_first"``，此时第 1 条是碱基；
     ③ （2026-10-31，用户已批准）没有任何 ``role`` 又是 ``_`` 写法：**输入里显式写出的羟基严格更多的那条是碱基**（``_base_by_hydroxyl``；LIPID+ 原生
       ``Cer 24:0_18:1;OH2``、``SM 18:1;OH2_16:0``）。只数 ``Chain.oh_read``（``;OH`` ``;OH2`` ``;(2OH)``）；任何一条链上有包补出来的羟基
       （MS-DIAL 的 ``;O2`` 贴链写法按鞘脂政策补成 ``OH``、PNNL 后缀 / 缺氧数补的碱基氧、亚类字母的类隐含氧）、有泛化的 ``O`` token、或来源不明
       （``oh_read is None``）都不推断；两链显式 OH 数相同（含都为 0）也不推断。推断成功时 ``LmResult.base_inferred`` 为真，api 加 flag ``sphingo_base_inferred``。
   其余情形**不猜**；
e. ``fmt ∈ {species_only, n_acyl, o_acyl}`` 不出链级（``fmt_species_only`` / ``fmt_n_acyl`` / ``fmt_o_acyl``）：N-酰基 / O-酰基的 2020 写法
   本包读不回，LIPID+ 的 L2/L3 也不标哪条链是 N-/O-酰基；
f. 各链 C / DB / 氧之和 = L1（``sum_mismatch``，防御）。

渲染：每条链 ``C:DB`` + 氧尾缀（链上氧价 1 -> ``;O``，n -> ``;On``；不写 OH/oxo 等官能团，那是更高一级的信息）；
``ether=="O"`` -> ``O-C:DB``；``ether=="P"`` -> ``P-C:(DB-1)``（``Chain.DB`` 已含缩醛 +1，写回字面）。
``standard``：醚链在前，其余按 (C, DB, 氧) 升序，``_`` 连接；溶血类单链、不写 ``0:0``。``sphingo``：碱基 ``/`` N-酰基。

胆固醇酯 ``fmt=sterol_ester``（用户规则 4，**按输入写**）
    输入类头是 ``CE`` / ``ChE`` -> 两级都写 ``CE <酰基链>``（``CE 18:1``）；总量级按 ``@adj``（C-27, DB-1, O-1）从 L1 换算；
        链级要求甾醇链是 27:1 且恰好 1 个氧（胆固醇 3-OH），酰基链 = 另一条。
    其他输入（``SE``、原生 ``ST-oacyl``…）-> 两级都写 ``SE <甾醇链>/<酰基链>``（``SE 27:1/18:1``），**只有总量时两个字段都留空**
        （``sterol_ester_sum_only``）。甾醇链的 ``;O`` 记 氧价-1（3-O 已含在 ``SE`` 里）。

公共接口
--------
``lm_names(cf, lb, it, head="", ordered=True, methyl_branch=False) -> LmResult``
                                                                  永不抛异常（内部错误 -> 两个字段为空、原因码 ``lm_internal_error``）；
                                                                  ``methyl_branch``：有链带位置甲基支链 ``(15Me)`` -> 两个字段留空（``branch_methyl_carbon``）
``heads() -> dict[str, LmHead]``                                  读 ``lm_heads.csv``（进程内缓存）
``LmResult``（``species molecular head_source species_reason molecular_reason base_inferred``；``.inferred``、``.to_meta()``）
``LM_REASONS``                                                    原因码全集；``FLAGS_lm``  本模块的 flag（``lm_head_inferred`` ``sphingo_base_inferred``）
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Final

from .contract import register_family_flags
from .labels import fmt_l1_o
from .nosn import chain_valence, is_ether_class
from .semantics.families.fa import CHOLESTEROL, STEROL_CHAINS
from .types import Chain, Interp, Labels, Struct

__all__ = ["lm_names", "heads", "LmHead", "LmResult", "LM_REASONS", "FLAGS_lm", "SOURCES", "FORMATS", "CE_HEADS"]

#: ``lm_head_inferred``：类名是推断的（数据库没有）；``sphingo_base_inferred``：鞘脂链级名的碱基是按"显式羟基更多的那条"推断的（条件 d③）
FLAGS_lm: frozenset[str] = frozenset({"lm_head_inferred", "sphingo_base_inferred"})
register_family_flags("lm", FLAGS_lm)

CSV_PATH: Final = Path(__file__).resolve().parent / "data" / "lm_heads.csv"
SOURCES: Final = ("lmsd", "inferred", "none")
FORMATS: Final = ("standard", "sphingo", "species_only", "n_acyl", "o_acyl", "sterol_ester", "fixed")
#: 输入类头是这些时，胆固醇酯写 ``CE <酰基链>``（用户规则 4：输入写 CE / ChE -> CE）；其他输入写 ``SE <甾醇链>/<酰基链>``。
#: ``ChE``（LipidSearch / LipidView 的胆固醇酯写法）2026-10-31 随别名行 ``alias/fa.csv`` 加入，别名行与 ``CE`` 同（inferred、带入胆固醇骨架）
CE_HEADS: Final = frozenset({"CE", "ChE"})

#: 原因码（``meta.lm.species`` / ``meta.lm.molecular``）。``@blank=`` 记号的取值也必须在这里
LM_REASONS: Final = frozenset({
    # 两个字段都空
    "no_lm_head",                     # source=none：数据库和推断都给不出类名
    "readback_unverified_alias",      # 该类的 LM 名在默认配置下被别名表拒绝（class_alias_unverified），读不回
    "readback_class_not_in_vocab",    # 该类的 LM 名读成词表外的类（NAGly 等单链 N-酰基氨基酸）
    "adj_out_of_range",               # @adj 换算出负数
    "sterol_ester_sum_only",          # SE / 原生 ST-oacyl 只有总量：LM 名要写甾醇链，写不出
    "branch_methyl_carbon",           # 位置甲基支链 (15Me)：L1 不计甲基碳，LM 要计，对不上
    "lm_internal_error",              # 本模块的 bug（测试里恒为 0）
    # 只有 LM_molecular 空
    "fmt_species_only", "fmt_n_acyl", "fmt_o_acyl",
    "sum_only", "chain_count", "oxygen_off_chain", "ether_chain_unknown", "ether_chain_in_non_ether_class",
    "p_chain_without_db", "sphingo_base_unknown", "sphingo_chain_count", "hypotheses_differ", "sum_mismatch",
    "sterol_unknown", "sterol_not_cholesterol", "sterol_oxygen_missing", "no_structure",
})

_TOKEN: Final = re.compile(r"@(\w+)=(\S+)")
_ADJ: Final = re.compile(r"(C|DB|O)([+-]\d+)")


# ==========================================================================
# 类名表
# ==========================================================================
@dataclass(frozen=True, slots=True)
class LmHead:
    """``lm_heads.csv`` 的一行（machine 记号已解析）。``head`` 为空 <=> ``source == "none"``。"""
    class_full: str
    head: str
    source: str
    fmt: str
    adj: tuple[int, int, int] = (0, 0, 0)       # (C, DB, O) 偏移：LM 值 = L1 值 + adj
    suffix: str = ""
    render: str = ""
    blank: str = ""


def _parse_adj(s: str) -> tuple[int, int, int]:
    d = {"C": 0, "DB": 0, "O": 0}
    for part in s.split(",") if s else ():
        m = _ADJ.fullmatch(part)
        if m is None:
            raise ValueError(f"lm_heads.csv: bad @adj part {part!r}")
        d[m.group(1)] += int(m.group(2))
    return d["C"], d["DB"], d["O"]


@lru_cache(maxsize=1)
def heads() -> dict[str, LmHead]:
    """读 ``data/lm_heads.csv``。格式有误（未知 source / fmt / 原因码、重复类）抛 ``ValueError``。用 ``csv`` 模块（``FA-amide`` 的头就是字面 ``NA``）。"""
    out: dict[str, LmHead] = {}
    with CSV_PATH.open(encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh):
            cf = r["class_full"]
            if cf in out:
                raise ValueError(f"lm_heads.csv: duplicate class {cf!r}")
            if r["source"] not in SOURCES or r["fmt"] not in FORMATS:
                raise ValueError(f"lm_heads.csv: {cf}: bad source/fmt {r['source']!r}/{r['fmt']!r}")
            if (r["source"] == "none") != (r["lm_head"] == ""):
                raise ValueError(f"lm_heads.csv: {cf}: source=none <=> empty lm_head violated")
            tok = dict(_TOKEN.findall(r["note"].split(" | ", 1)[0]))
            blank = tok.get("blank", "")
            if blank and blank not in LM_REASONS:
                raise ValueError(f"lm_heads.csv: {cf}: unknown @blank reason {blank!r}")
            out[cf] = LmHead(cf, r["lm_head"], r["source"], r["fmt"], _parse_adj(tok.get("adj", "")),
                             tok.get("suffix", ""), tok.get("render", ""), blank)
    return out


# ==========================================================================
# 结果
# ==========================================================================
@dataclass(frozen=True, slots=True)
class LmResult:
    """``lm_names`` 的返回。``species_reason`` / ``molecular_reason`` 只在对应字段为空时有意义；``base_inferred`` 只在 ``molecular`` 非空时有意义。"""
    species: str = ""
    molecular: str = ""
    head_source: str = "none"
    species_reason: str = ""
    molecular_reason: str = ""
    base_inferred: bool = False           # 链级名的鞘脂碱基是按显式羟基数推断的（条件 d③）；api 据此加 flag ``sphingo_base_inferred``

    @property
    def inferred(self) -> bool:
        """类名是推断的（数据库没有）。"""
        return self.head_source == "inferred"

    def to_meta(self) -> dict[str, str]:
        """``meta.lm``：``head_source`` 恒有；字段为空时附原因码。"""
        m = {"head_source": self.head_source}
        if not self.species and self.species_reason:
            m["species"] = self.species_reason
        if not self.molecular and self.molecular_reason:
            m["molecular"] = self.molecular_reason
        return m


_EMPTY: Final = LmResult()


def _blank(src: str, why: str) -> LmResult:
    return LmResult("", "", src, why, why)


# ==========================================================================
# 总量级
# ==========================================================================
def _species(h: LmHead, cf: str, lb: Labels) -> tuple[str, str]:
    """``(名字, 原因码)``：``<head> [O-]C:DB[;On]<suffix>``，数值 = L1 + adj。"""
    sm = lb.sum_meta
    c, db, o = sm["C"] + h.adj[0], sm["DB"] + h.adj[1], sm["O"] + h.adj[2]
    if min(c, db, o) < 0:
        return "", "adj_out_of_range"
    eth = "O-" if is_ether_class(cf) else ""
    return f"{h.head} {eth}{c}:{db}{fmt_l1_o(o)}{h.suffix}", ""


# ==========================================================================
# 链级
# ==========================================================================
def _chain_text(c: Chain, n_o: int) -> str | None:
    """一条链的 2020 写法 ``[O-|P-]C:DB[;On]``；``P`` 写回字面 DB（``Chain.DB`` 已 +1）。``P`` 链 DB<1 返回 None。"""
    if c.ether == "P":
        if c.DB < 1:
            return None
        core = f"P-{c.C}:{c.DB - 1}"
    elif c.ether == "O":
        core = f"O-{c.C}:{c.DB}"
    else:
        core = f"{c.C}:{c.DB}"
    return core + fmt_l1_o(n_o)


def _base_index(st: Struct, cslash: str, ordered: bool) -> int | None:
    """鞘脂碱基链的下标：恰一条 ``role == "base"``；或没有任何 role 时，有序写法 + ``slash_role == base_first`` -> 第 1 条。其余 None。"""
    chains = st.chains or ()
    idx = [i for i, c in enumerate(chains) if c.role == "base"]
    if len(idx) == 1:
        return idx[0]
    if idx:
        return None
    if ordered and cslash == "base_first" and len(chains) == 2 and not any(c.role for c in chains):
        return 0
    return None


def _base_by_hydroxyl(st: Struct) -> int | None:
    """鞘脂碱基的下标，依据"输入里**显式写出**的羟基严格更多的那条是碱基"（条件 d③）；不能可靠地比较 -> None。

    要求（任何一条不满足都不推断，宁可留空）：
    * 恰好两条链、都没有 ``role``（有 ``d/t/m`` 前缀 / ``N-`` 前缀的名字不走这里）、没有链外氧；
    * 每条链的 ``OH`` token 都是输入显式写出的（``Chain.oh_read == OH 个数``；``oh_read is None``＝来源不明，也算不满足）：
      包按鞘脂政策把 ``;O`` / ``;On`` 补成 ``OH``（MS-DIAL ``Cer 18:1;O2_24:0``，基准库里约 13% 氧位置不可靠）、PNNL 后缀 / 缺氧数补的碱基氧、
      亚类字母的类隐含氧、家族钩子加的 ``OH`` 都会让 ``oh_read`` 小于 ``OH`` 个数；
    * 没有泛化的 ``O`` token（可能是羟基也可能不是）；
    * 两链显式 ``OH`` 数不同（相同含都为 0 -> 留空）。
    其他氧（``oxo`` ``Ep`` …）不参与比较。证据：参考表两链鞘脂 157,928 个唯一 L3 里碱基的 OH 数严格大于酰基链，0 反例、0 平局（``docs/CONTROVERSIES.md`` LM-2）。"""
    chains = st.chains or ()
    if len(chains) != 2 or st.tail_o or any(c.role for c in chains):
        return None
    n_oh: list[int] = []
    for c in chains:
        t = dict(c.toks)
        oh = t.get("OH", 0)
        if t.get("O", 0) or c.oh_read != oh:
            return None
        n_oh.append(oh)
    if n_oh[0] == n_oh[1]:
        return None
    return 0 if n_oh[0] > n_oh[1] else 1


def _sterol_index(st: Struct) -> int | None:
    """甾醇链的下标：恰一条 ``role == "sterol"``；没有 role 时（原生写法）恰一条 (C, DB) 属于常见甾醇集合（同 ``fa.role_candidates`` 的口径）。"""
    chains = st.chains or ()
    idx = [i for i, c in enumerate(chains) if c.role == "sterol"]
    if len(idx) == 1:
        return idx[0]
    if idx or any(c.role for c in chains):
        return None
    cand = [i for i, c in enumerate(chains) if (c.C, c.DB) in STEROL_CHAINS]
    return cand[0] if len(cand) == 1 else None


def _render(h: LmHead, cf: str, st: Struct, nominal: frozenset[int], lb: Labels, cslash: str, ordered: bool,
            ce_form: bool) -> tuple[str, str, bool]:
    """一个假设 -> ``(链级名, 原因码, 碱基是推断的)``（见模块文档的条件 a-f；第三项只有鞘脂按条件 d③ 推断碱基时为真）。"""
    chains = st.chains
    if chains is None:
        return "", "sum_only", False
    if nominal and len(chains) not in nominal:
        return "", "chain_count", False
    single = nominal == frozenset({1})
    if st.tail_o and not single:
        return "", "oxygen_off_chain", False
    eth = [c for c in chains if c.ether]
    if is_ether_class(cf):
        if len(eth) != 1:
            return "", "ether_chain_unknown", False
    elif eth:
        return "", "ether_chain_in_non_ether_class", False
    ox = [chain_valence(c) + (st.tail_o if single else 0) for c in chains]
    sm = lb.sum_meta
    if (sum(c.C for c in chains), sum(c.DB for c in chains), sum(ox)) != (sm["C"], sm["DB"], sm["O"]):
        return "", "sum_mismatch", False

    fmt = h.fmt
    if fmt == "sphingo":
        if len(chains) != 2:
            return "", "sphingo_chain_count", False
        b = _base_index(st, cslash, ordered)
        by_oh = False
        if b is None:
            b = _base_by_hydroxyl(st)
            by_oh = b is not None
        if b is None:
            return "", "sphingo_base_unknown", False
        tb, ta = _chain_text(chains[b], ox[b]), _chain_text(chains[1 - b], ox[1 - b])
        if tb is None or ta is None:
            return "", "p_chain_without_db", False
        return f"{h.head} {tb}/{ta}", "", by_oh
    if fmt == "sterol_ester":
        if len(chains) != 2:
            return "", "chain_count", False
        s = _sterol_index(st)
        if s is None:
            return "", "sterol_unknown", False
        if ox[s] < 1:
            return "", "sterol_oxygen_missing", False          # 3-O 必须在甾醇链上（SE 里已含）
        a = 1 - s
        ta = _chain_text(chains[a], ox[a])
        if ta is None:
            return "", "p_chain_without_db", False
        if ce_form:
            if (chains[s].C, chains[s].DB) != CHOLESTEROL or ox[s] != 1:
                return "", "sterol_not_cholesterol", False
            return f"{h.head} {ta}", "", False
        ts = _chain_text(chains[s], ox[s] - 1)
        return f"SE {ts}/{ta}", "", False
    # standard
    order = sorted(range(len(chains)), key=lambda i: (0 if chains[i].ether else 1, chains[i].C, chains[i].DB, ox[i]))
    parts = [_chain_text(chains[i], ox[i]) for i in order]
    if None in parts:
        return "", "p_chain_without_db", False
    return f"{h.head} " + "_".join(parts), "", False


def _molecular(h: LmHead, cf: str, it: Interp | None, lb: Labels, ordered: bool, ce_form: bool = False) -> tuple[str, str, bool]:
    """``(链级名, 原因码, 碱基是推断的)``；所有假设都推断才算推断（实际上同一个输入的各假设对碱基的判断一致）。"""
    if h.fmt in ("species_only", "n_acyl", "o_acyl"):
        return "", f"fmt_{h.fmt}", False
    if it is None:
        return "", "no_structure", False
    nominal, cslash = it.cls.nchain, it.cls.slash_role
    out: set[str] = set()
    inferred: list[bool] = []
    for st in it.alts:
        s, why, inf = _render(h, cf, st, nominal, lb, cslash, ordered, ce_form)
        if why:
            return "", why, False
        out.add(s)
        inferred.append(inf)
    if len(out) != 1:
        return "", "hypotheses_differ", False
    return next(iter(out)), "", any(inferred)


# ==========================================================================
# 入口
# ==========================================================================
def _coq(h: LmHead, lb: Labels) -> LmResult:
    sm = lb.sum_meta
    c, db, o = sm["C"], sm["DB"], sm["O"]
    n = db - 2
    if h.render and o == 4 and 1 <= n <= 13 and c == 9 + 5 * n:
        s = h.render.replace("{n}", str(n))
    else:
        s = f"{h.head} {c}:{db}{fmt_l1_o(o)}"
    return LmResult(s, s, h.source)


def _names(h: LmHead, cf: str, lb: Labels, it: Interp | None, head: str, ordered: bool) -> LmResult:
    src = h.source
    if h.fmt == "fixed":
        return _coq(h, lb)
    if h.fmt == "sterol_ester":
        if head in CE_HEADS:
            sp, why = _species(h, cf, lb)
            if why:
                return _blank(src, why)
            mol, mwhy, _ = _molecular(h, cf, it, lb, ordered, ce_form=True)
            return LmResult(sp, mol, src, "", mwhy)
        mol, why, _ = _molecular(h, cf, it, lb, ordered, ce_form=False)
        if why == "sum_only":
            why = "sterol_ester_sum_only"
        return LmResult(mol, mol, src, why, why) if mol else _blank(src, why)
    sp, why = _species(h, cf, lb)
    if why:
        return _blank(src, why)
    mol, mwhy, by_oh = _molecular(h, cf, it, lb, ordered)
    return LmResult(sp, mol, src, "", mwhy, by_oh and bool(mol))


def lm_names(cf: str, lb: Labels, it: Interp | None, head: str = "", ordered: bool = True,
             methyl_branch: bool = False) -> LmResult:
    """``cf`` = 输出的 ``class_full``；``lb`` = 该段的 ``Labels``；``it`` = 该段的 ``Interp``（固定结构类为 None）；
    ``head`` = 输入的类头（``HeadInfo.head``，规范化后、大小写敏感；固定结构类传 ``""``）；
    ``ordered`` = 主链之间没有 ``_`` 分隔符（有序写法，``/`` 才能指认鞘脂碱基）；
    ``methyl_branch`` = 有链带位置甲基支链 ``(15Me)``：L1 不计这个甲基的 C（``branch_dropped``，LIPID+ 口径），LIPID MAPS 的总量名要计，
    两边对不上 -> 两个字段都留空（``branch_methyl_carbon``）。``iso`` / ``anteiso`` 的 C 已含支链碳，不受影响。永不抛异常。"""
    src = "none"
    try:
        h = heads().get(cf)
        if h is None or h.source == "none":
            return _blank("none", "no_lm_head")
        src = h.source
        if h.blank:
            return _blank(src, h.blank)
        if methyl_branch:
            return _blank(src, "branch_methyl_carbon")
        return _names(h, cf, lb, it, head, ordered)
    except Exception:                                                    # noqa: BLE001 - 永不抛；测试里要求这个码出现 0 次
        return _blank(src, "lm_internal_error")
