"""甘油脂（GL）家族钩子（B2）：MG / DG / TG、DGCC / DGTS / DGTA（及溶血）、糖基 / 硫酸 / 酰基甘油脂（共 23 个 class_full）。

设计结论（``docs/census/fam_gl.md`` 的实测：纯规则原型在 GT 299 名与参考表 305,828 行上逐字复现）：本家族的绝大多数规则
**是表驱动的**，不需要代码例外——

* 类 token -> class_full：``data/alias/{ref_main,tools,native,gl}.csv``（MGDG = DG-Hex、DGDG = DG-Hex2、SQDG = DG-sq、
  DGGA / GlcADG = DG-HexA、ADGGA = DG-HexA-oacyl、LDGCC = DGCC-lyso、LDGTS / MGTS = DGTS-lyso、DGMG = MG-Hex2 ……；
  糖 / 硫酸 / 酰基 / 溶血全由**类 token** 触发，不由链语法触发）。词表外一律在 M4 拒绝：``TG-EST`` / ``TG-OAHFA``、``SQDG-O``、
  ``ADGGA-O``、``SMGDG``（无醚）、``DGTSA``（歧义）、带 ``O-`` 的 ``SQDG`` / ``ADGGA``。
* 类的名义链数、醚 / 溶血标志：``data/class_info/gl.csv``。类自带的氧（TG 6、DG 5、MGDG 10 ……）从不写进 L1/L2/L3，所以
  ``implied_O`` 全为 0、氧策略全为 ``typed_only``。
* 醚种类（``O-`` / ``P-`` / ``e`` / 类头 ``plasmanyl- plasmenyl- -O -P``）、``P`` 的 DB+1、``0:0`` 占位丢弃（1,2- / 1,3- / 2,3-sn
  去 sn 后同一标签）、链数核对、``-O`` 后缀的"每条链各当一次醚链"假设枚举（L2 可达、L3 仅当各假设一致，如链全同）、
  未分型 ``;O`` -> 只到 L2、sum 级 -> 只到 L1：全部在 ``semantics/core.py`` / ``conventions.py``。

本钩子只处理四类 GL 例外（``BaseHandler`` 接口见 ``families/__init__.py``）：

``select_hint``（D13 / D4）
    0. ``TG-oacyl`` 带 2 个嵌套酰基 -> 换成 ``TG-oacyl2``（别名表按 ``(FA c:d)`` 有无选类，数不了个数）。
    ``Config.single_token_is_lyso(source)``（默认仅 ``source="pnnl"``）为真，名字只有一个 C:DB 链组，类是两链类而同核溶血
    兄弟类是单链类（``classes.class_sibling(cf, lyso=True)``；GL 里只有 ``DGCC -> DGCC-lyso``、``DGTS -> DGTS-lyso``）
    -> 换成溶血类，core 打 ``source_hint_lyso``。与 GPL 的 PNNL 溶血提示同一口径；PNNL 名字表里没有这两类的实例（只有
    ``DGTSA``，被拒绝），所以这条是 D4 的外推，未被 PNNL 数据验证。
``build``（在默认装配前做 GL 专属的拒绝 / 改写）
    1. 嵌套 ``(FA c:d)`` 或括号 ``(O-c:d)`` 酰基前缀 -> ``unsupported_syntax``，**除了** ``TG-oacyl`` / ``TG-oacyl2``（见 4）。core 的通用闸只拒绝 *没有* oacyl 的类，
       但 GL 里的 oacyl 类 ``DG-HexA-oacyl``（``ADGGA``）的第三条链在参考表里是普通链（``ADGGA 16:0_16:0_18:1``），嵌套 / 括号前缀写法
       没有来源，仍拒绝（``TG 16:0_18:1(FA 16:0)`` 没有宿主氧，本来就被拒绝）。
       不用 ``class_not_in_vocab``：检查器认为头部是词表类（``TG``）的名字不可能是词表外的类。
    2. ``N-`` 酰胺链前缀 -> ``class_chain_conflict``（GL 没有 N-酰基链）。
    4. ``TG-oacyl`` / ``TG-oacyl2``（0.4.0，词表 v0.2.7 新增：TG 的某条链上的羟基被另一个脂肪酸酯化，即 estolide / O-酰基 TG）的嵌套写法
       ``TG 16:0;5O(FA 16:0)/16:0/16:0``（LIPID MAPS，数据集参考表 17 个名字）与 ``TG-EST 8:4_8:0_10:0;O(FA 16:0)``（LipidIN）：宿主链上的 ``;nO`` 在这里
       是**酯化位置**（LIPID MAPS 的 ``5O(FA 16:0)`` = 5 位氧上挂 FA 16:0），不是氧个数；每个嵌套酰基恰有一个宿主氧，所以宿主链记一个 ``OH``
       （数据集口径：被酯化的羟基氧留在所在链上，L3 记 ``OH``；flag ``class_implied_OH``），O-酰基链是另一条无氧链。要求 3 条主链、
       宿主链数 = 嵌套酰基数、宿主链上除 ``;nO`` 外没有别的氧记号，否则 ``unsupported_syntax``。两个嵌套酰基 = ``TG-oacyl2``（``select_hint`` 换类）。
    其余一律 ``return None``，交给 core 的默认装配。缩醛 ``;P`` 落在 DB 为 0 的链上 -> ``malformed`` 原来也在这里，现在是 core
    （``conventions.chain_mark``）对所有家族通用的规则：``;P`` 链的 DB 已含乙烯醚那个双键，``O-16:0;P`` 结构上不存在。

新增 flag ``ester_position_read``（0.4.0）：TG-oacyl 宿主链上的 ``;nO``（n != 1）按酯化位置读、不是氧个数。醚链未知 / P 落点未知 / 氧未分型等状态用已有全局 flag
（``ether_chain_unknown`` ``p_landing_ambiguous`` ``oxygen_untyped`` ``sum_only``）。
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ...classes import class_info, class_sibling
from ...config import Config
from ...contract import register_family_flags
from ...types import ClassRes, HeadInfo, Lexed, Reject, Struct
from .. import conventions as cv
from . import BaseHandler

if TYPE_CHECKING:  # pragma: no cover
    from ..model import Ctx

__all__ = ["FLAGS_gl", "Handler"]

#: ``ester_position_read``：``TG-oacyl`` 宿主链上的 ``;nO``（n != 1）按酯化位置读、不是氧个数（L1 的氧只有 1 个）
FLAGS_gl: frozenset[str] = frozenset({"ester_position_read"})
register_family_flags("gl", FLAGS_gl)

#: 嵌套 ``(FA c:d)`` 酰基的个数 -> TG 的 O-酰基类（``TG-oacyl`` 一个、``TG-oacyl2`` 两个）；主链数恒为 3
TG_OACYL: dict[str, int] = {"TG-oacyl": 1, "TG-oacyl2": 2}
TG_OACYL_MAIN_CHAINS = 3


class Handler(BaseHandler):
    """GL 家族钩子，见模块文档。"""

    family = "gl"

    # ---- D13：类提示 ----------------------------------------------------------------------
    def allows_nested_acyl(self, cr: ClassRes) -> bool:
        return cr.class_full in TG_OACYL or super().allows_nested_acyl(cr)

    def select_hint(self, hi: HeadInfo, lx: Lexed, cr: ClassRes, cfg: Config, source: str | None) -> ClassRes | None:
        if cr.class_full == "TG-oacyl" and len(lx.nested_fa) == 2:                 # 两个嵌套酰基：TG-oacyl2
            return class_info("TG-oacyl2").to_class_res(cr.alias_status, cr.inject)
        if not lx.is_sum or cr.lyso or cr.nchain != frozenset({2}):
            return None
        if lx.feats.n_placeholder or lx.nested_fa or lx.pre_acyl is not None:
            return None
        if not cfg.single_token_is_lyso(source):
            return None
        sib = class_sibling(cr.class_full, lyso=True)
        if sib is None:
            return None
        ci = class_info(sib)
        if ci.family != self.family or ci.nchain != frozenset({1}):
            return None
        return ci.to_class_res(cr.alias_status, cr.inject)

    # ---- 装配前的 GL 专属拒绝 -------------------------------------------------------------
    def build(self, lx: Lexed, cr: ClassRes, cfg: Config, ctx: "Ctx") -> Struct | list[Struct] | None:
        if cr.class_full in TG_OACYL and lx.nested_fa and lx.pre_acyl is None:
            self._tg_oacyl_hosts(lx, cr, cfg, ctx)
        elif lx.nested_fa or lx.pre_acyl is not None:
            raise Reject("unsupported_syntax",
                         f"{cr.class_full}: an (FA ..) / (O-..) acyl group has no source in the glycerolipid names "
                         "(only the TG-oacyl classes take one, on a chain that carries the esterified oxygen; "
                         "the O-acyl prefix is the AHexCer notation)")
        for w in ctx.works:
            if w.n_acyl:
                raise Reject("class_chain_conflict", f"{cr.class_full} has no N-acyl (amide) chain: {w.raw.raw if w.raw else ''!r}")
        return None

    # ---- TG-oacyl / TG-oacyl2 的宿主氧 ----------------------------------------------------
    @staticmethod
    def _tg_oacyl_hosts(lx: Lexed, cr: ClassRes, cfg: Config, ctx: "Ctx") -> None:
        """嵌套 ``(FA c:d)`` 的宿主链：``;nO`` 是酯化位置，记一个 ``OH``（模块文档第 4 条）。不满足形状 -> ``unsupported_syntax``。"""
        want = TG_OACYL[cr.class_full]
        if len(lx.nested_fa) != want or len(ctx.works) != TG_OACYL_MAIN_CHAINS or ctx.sum_mode:
            raise Reject("unsupported_syntax", f"{cr.class_full}: needs {TG_OACYL_MAIN_CHAINS} main chains and {want} nested (FA ..) group(s), "
                                               f"got {len(ctx.works)} and {len(lx.nested_fa)}")
        hosts = [w for w in ctx.works if w.raw is not None and any(a.kind == "O_count" for a in w.raw.atoms)]
        if len(hosts) != want:
            raise Reject("unsupported_syntax", f"{cr.class_full}: {want} nested (FA ..) group(s) need {want} chain(s) carrying the esterified "
                                               f"oxygen (';O(FA ..)'), got {len(hosts)}")
        for w in hosts:
            assert w.raw is not None
            if any(a.kind != "O_count" for a in w.raw.atoms) or w.mark or w.base or w.n_acyl:
                raise Reject("unsupported_syntax", f"{cr.class_full}: the chain carrying the esterified oxygen has other modifiers: {w.raw.raw!r}")
            _toks, written = cv.atoms_to_tokens(w.raw.atoms, cr.oxy_policy)            # typed_only：整条链的氧都在 tail_o 里
            ctx.tail_o -= written
            if written != 1:
                ctx.flags.add("ester_position_read")                                    # ``;5O(FA ..)``：5 是位置，L1 只有 1 个氧
            w.toks["OH"] = w.toks.get("OH", 0) + 1
            w.role = "oacyl_host"
        ctx.flags.add("class_implied_OH")
