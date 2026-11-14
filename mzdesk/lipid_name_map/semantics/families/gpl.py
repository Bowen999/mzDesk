"""甘油磷脂（GPL）家族钩子（B1）：``Handler`` 只处理本家族的例外，通用规则在 ``semantics/core.py`` 与类表。

家族范围
--------
``class_info/gpl.csv`` 的 70 个 ``class_full``（词表 v0.2.11 的 PC PE PG PI PS PA CL BMP CDP-DAG As-PL PT PnC PnE CPA 核；
``lyso / ether / nacyl / oacyl / phos / Hex / me / OMe / OEt`` 变体）。``AMP`` 与 ``WD`` 归 ``fa`` 家族。

分工（谁负责什么；本家族的规则几乎都是**表驱动**的，钩子只补两处）
-------------------------------------------------------------------
类表（``data/class_info/gpl.csv``、``data/alias/*.csv``）
    * 变体选择：头部 × 醚标记（``ether``）× ``0:0`` 个数（``zero``）-> ``class_full``：
      ``PC(0:0/x)`` ≡ ``LPC x`` -> ``PC-lyso``；NAPE 四变体（醚 × ``0:0``）-> ``PE-nacyl`` / ``PE-ether-nacyl`` / ``PE-nacyl-lyso``
      / ``PE-ether-nacyl-lyso``；``CL`` + 1/2 个 ``0:0`` -> ``CL-lyso`` / ``CL-lyso2``（4 槽）；``BMP`` 4 槽（2 个 ``0:0``）；
      ``PIP/PIP2/PIP3`` -> ``PI-phos*``；``PIMn / Ac2PIMn / LPIMn``；``PEtOH PMeOH MMPE DMPE``；v26 旧名（``lyso1`` ``plasmalogen``
      ``PA-me`` ``BMP-nacyl`` …）。
    * 链数 ``nchain``（丢掉 ``0:0`` 后的槽数）：CL 4、CL-lyso 3、CL-lyso2 2、NAPE 3、BMP-oacyl 3、PI-Hex2-oacyl 3 / oacyl2 4、其余 2、溶血 1。
``core``（``semantics/core.py`` / ``conventions.py``，B0）
    * 醚种类合并（链 ``O- P- e`` / 类头 ``Plasmanyl- Plasmenyl- plasmalogen -O -P``），``P`` 的 DB+1 只加一次；
    * 醚链位置未知（``-O`` / ``-P`` / ``Plasmenyl-X`` 无链标记 / ``PE-ether`` 无标记）-> 假设枚举（``ether_chain_unknown`` /
      ``p_landing_ambiguous``）：L1 恒可出，L2 只有所有落点同串才出，L3 同理；
    * 未分型氧 ``;O`` ``;2O`` -> ``tail_o``（L2 ``;(On)``，L3 空）；已分型（``;OH`` ``(OH)`` Riken 介质链）才进 L3；
    * 介质缩写（``12-HETE`` 等 7 个）展开；``0:0`` 占位丢弃；总量名（一个 C:DB 而类需要更多链）-> ``chains=None``（level 1）。

本模块的钩子
------------
``Handler.select_hint(hi, lx, cr, cfg, source) -> ClassRes | None``      （D13）
    ``Config.single_token_is_lyso(source)``（默认仅 ``source="pnnl"``，D4）为真，且名字只有一个 C:DB 链组、类是两链类
    （``nchain == {2}``）而"仅溶血不同"的兄弟类是单链类（``classes.class_sibling(cf, lyso=True)``：``PC -> PC-lyso``、
    ``PC-ether -> PC-ether-lyso``）-> 换成溶血类（core 打 ``source_hint_lyso``）。
    PNNL 里 ``PC 16:0`` 是 LPC（1,509 谱 = 3.3%），LIPID MAPS 里是总量 ``PC 8:0/8:0``，所以默认仍当总量。
    ``single_token_as_lyso="auto"``（默认）时只对 PNNL 实测过的类生效：核 PC PE PG PI PS PA 的普通类与醚类
    （PNNL 数据里没有 PIM / PE-Hex / PA-Hex 等，它们的单 C:DB 仍当总量）；显式 ``True`` 则对所有有单链兄弟的类生效。
    没有溶血兄弟（``CL`` ``BMP`` ``PE-nacyl`` …）或带 ``0:0`` / ``(FA)`` / O-酰基前缀的名字不换类。
``Handler.build(lx, cr, cfg, ctx) -> None``
    N-酰基标记一致性（``N-`` 前缀只在 ``*-nacyl*`` 类里才有意义；core 对 ``N-`` 只用来把该链排除出醚链候选，不会拒绝）：
    * 非 nacyl 类里出现 ``N-`` 链（``PC 18:1_N-16:0``）-> ``class_chain_conflict``，否则 ``N-`` 会被静默丢掉；
    * nacyl 类里 ``N-`` 链多于一条，或 ``N-`` 链同时带醚标记（``N-O-18:1``，N-酰基不可能是醚链）-> ``class_chain_conflict``。
    其余返回 ``None``（core 默认装配）。

``Handler.allows_adjacent_groups(cr)``
    ``CL`` -> True：相邻括号组 ``CL(28:0)(28:0)`` 只有 CL 对记法有定义（其余类 core 拒绝 ``malformed``）。
``Handler.build`` 的 CL 对记法（D24）
    ``CL 36:4_36:2``（MS-DIAL MSP）、``CL(28:0)(28:0)``（LipidMatch，``+NH4`` 等加合物由 M2 先剥掉）：``CL`` 头 + 恰好 2 个链组
    （``_`` 分隔或相邻括号组；斜线 ``CL 36:4/36:2`` 同收——两段都 >= C28 时它只可能是两个磷脂酰基，检查器的 B4 元变换 ``_`` <-> ``/`` 也要求
    一致）、没有 ``0:0`` / 氧 / 醚标记 / 嵌套链，每个组的 C >= ``PAIR_SEGMENT_MIN_C``：**每组是一个磷脂酰基
    （二酰，两条酰基链）的总组成**，按和式收，只到 level 1：``CL <ΣC>:<ΣDB>``，L2/L3 留空，flag ``pair_sum`` + ``sum_only``。
    不满足这些条件（组太短 -> 更像两条酰基链的写错，如 ``CL 18:2_18:2``；带氧）-> 不当对记法，走默认路径
    拒绝 ``chain_count``（相邻括号组则拒绝 ``malformed``）。``MLCL`` / ``DLCL`` 头（``CL-lyso`` / ``CL-lyso2``）与 4 链 ``CL`` 不受影响。

``FLAGS_gpl``   本家族自带 flag：``pair_sum``（CL 对记法，D24）。GPL 需要的其余 flag 全在基础词表里：``P_db_plus1 p_landing_ambiguous
ether_chain_unknown placeholder_dropped source_hint_lyso oxylipin_expanded oxygen_untyped alias_inferred`` …

不处理（有意的，见 ``docs/CONTROVERSIES.md`` B1 小节）
------------------------------------------------------
* PNNL ``PE-P a_b`` 的 +1 落点（``Config.pnnl_ether_heuristic`` 才套 benchmark 的"放最饱和链"）；
* ``PS-NAc`` / ``Glc-GP 18:0`` 等参考表标签来自结构的行（名字不可知）。
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Final

from ...classes import class_info, class_sibling
from ...config import Config
from ...contract import register_family_flags
from ...types import ClassRes, HeadInfo, Lexed, Reject, Struct
from .. import conventions as cv
from . import BaseHandler

if TYPE_CHECKING:  # pragma: no cover
    from ..model import Ctx

__all__ = ["Handler", "FLAGS_gpl", "is_nacyl_class", "PAIR_SEGMENT_MIN_C"]

#: ``pair_sum``：CL 对记法（D24）——名字里 CL 的每一段都是一个磷脂酰基的总组成
FLAGS_gpl: frozenset[str] = frozenset({"pair_sum"})
register_family_flags("gpl", FLAGS_gpl)

#: CL 对记法每一段（一个二酰磷脂酰基）的 C 下限：两条酰基链、每条至少 C14（参考表 CL / MLCL / DLCL 里最短的酰基链是 14；
#: LipidMatch / MSP 库里最短的段是 28:0）。更短的段（``CL 18:2_18:2``）更可能是两条酰基链写成了 ``_``，不当对记法收
PAIR_SEGMENT_MIN_C: Final = 28

#: PNNL 的"单 C:DB = 溶血"约定实测过的核与 mod（fam_gpl §5.1、pnnl_name_to_benchmark_gt.csv）；其余类 auto 模式下不套用
_PNNL_LYSO_NUCLEI = frozenset({"PC", "PE", "PG", "PI", "PS", "PA"})
_PNNL_LYSO_MODS = frozenset({"", "ether"})


def is_nacyl_class(cr: ClassRes) -> bool:
    """类名的 ``mod`` 里有 ``nacyl`` 段（``PE-nacyl``、``PE-ether-nacyl-lyso``、``PE-nacyl-lyso2``、``PS-nacyl``…）。

    ``oacyl``（``BMP-oacyl``、``PI-Hex2-oacyl``）是 O-酰基，不是 N-酰基。
    """
    return "nacyl" in cr.mod.split("-")


class Handler(BaseHandler):
    family = "gpl"

    # ------------------------------------------------------------------ D13：PNNL 单链 -> 溶血
    def select_hint(self, hi: HeadInfo, lx: Lexed, cr: ClassRes, cfg: Config, source: str | None) -> ClassRes | None:
        if not lx.is_sum or cr.lyso or cr.nchain != frozenset({2}):
            return None
        if lx.feats.n_placeholder or lx.nested_fa or lx.pre_acyl is not None:
            return None
        if not cfg.single_token_is_lyso(source):
            return None
        if cfg.single_token_as_lyso == "auto" and not (cr.nucleus in _PNNL_LYSO_NUCLEI and cr.mod in _PNNL_LYSO_MODS):
            return None
        sib = class_sibling(cr.class_full, lyso=True)
        if sib is None:
            return None
        ci = class_info(sib)
        if ci.nchain != frozenset({1}):
            return None
        return ci.to_class_res(cr.alias_status, cr.inject)

    def allows_adjacent_groups(self, cr: ClassRes) -> bool:
        return cr.class_full == "CL"

    # ------------------------------------------------------------------ D24：CL 对记法 / N-酰基标记一致性
    def build(self, lx: Lexed, cr: ClassRes, cfg: Config, ctx: "Ctx") -> Struct | list[Struct] | None:
        if cr.class_full == "CL":
            pair = self._cl_pair_sum(lx, cr, ctx)
            if pair is not None:
                return pair
            if cv.adjacent_groups(lx):
                raise Reject("malformed", "adjacent parenthesised chain groups (a)(b) are only defined for the CL pair "
                                          f"notation (two plain diacyl totals, each C >= {PAIR_SEGMENT_MIN_C})")
        marked = [w for w in ctx.works if w.n_acyl]
        if not marked:
            return None
        if not is_nacyl_class(cr):
            raise Reject("class_chain_conflict",
                         f"N- (N-acyl) chain marker {marked[0].raw.raw if marked[0].raw else ''!r} in {cr.class_full}, "
                         f"which has no N-acyl chain")
        if len(marked) > 1:
            raise Reject("class_chain_conflict", f"{len(marked)} N-acyl chains in {cr.class_full} (it has one)")
        if marked[0].mark:
            raise Reject("class_chain_conflict", "an N-acyl chain cannot be the ether chain")
        return None

    # ------------------------------------------------------------------
    @staticmethod
    def _cl_pair_sum(lx: Lexed, cr: ClassRes, ctx: "Ctx") -> Struct | None:
        """D24：CL 对记法 -> level 1 和式（见模块文档）；不是对记法返回 None。"""
        ws = ctx.works
        if (lx.feats.n_placeholder or lx.nested_fa or lx.pre_acyl is not None or len(lx.chains) != 2 or len(ws) != 2
                or lx.seps not in (("_",), ("/",)) or ctx.n_groups != 2 or ctx.n_inject or ctx.ether_kind):
            return None
        if ctx.tail_o or ctx.pending or ctx.name_o:
            return None
        for w, rc in zip(ws, lx.chains):
            if (w.toks or w.mark or w.n_acyl or w.base or w.db_extra or rc.atoms or rc.pre or rc.suf or rc.branch
                    or rc.oxylipin or rc.base):
                return None
            if w.C < PAIR_SEGMENT_MIN_C:
                return None
        ctx.flags.update(("pair_sum", "sum_only"))
        return Struct(cr.class_full, None, total=(ws[0].C + ws[1].C, ws[0].DB + ws[1].DB, 0))

