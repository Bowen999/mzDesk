"""脂肪酰 / 固醇 / 其他（FA）家族钩子（B4）。依据 ``docs/census/fam_fa.md`` 与 ``data/class_info/fa.csv``。

家族范围：``FA FAL FOH CAR CoA VAE AMP WE WD``、N-酰基胺/氨基酸（``FA-amide FA-etn FA-Tau FA-<AA>``，0.4.0 起还有 ``FA-Xle``（原 Ile + Leu）与半胱氨酰白三烯 ``FA-GSH FA-SCys FA-SCysGly``）、两链 oacyl 类
（``FAHFA FA-Gly/GlySer/Orn-oacyl``）、固醇族（``ST`` 及其共轭 ``ST-Gly ST-Tau ST-sulf ST-Hex ST-Gly-sulf ST-Tau-sulf`` 与 0.4.0 起的氨基酸结合物 ``ST-<AA>``（20 个）、
``ST-oacyl`` = CE/SE、``ST-oacyl2``）、磺脂 ``SL``。

core 已经表驱动完成的规则（本文件**不**重复）
-----------------------------------------------
* 氧类型策略 ``typed_only``：只有类型化词（``OH oxo Ep COOH OOH OO``、``;(2OH)``、Riken 展开）才进 L3，裸 ``;On`` 只到 L2（C34）。
* ``class_OH``（FAHFA / ``*-oacyl`` / ``ST-oacyl``）：类隐含的那个氧落成角色链上的 ``OH``，``/`` 才能定位（FAHFA、N-酰基 oacyl 的第 2 链
  是羟基链；``ST-oacyl`` 的第 1 条是甾醇链），``_`` 枚举（L3 空）。``implied_mode``：FAHFA 与 N-酰基 oacyl 是 ``max``（名字氧已够则不加），
  ``ST-oacyl`` / ``ST-Gly`` / ``ST-Tau`` 是 ``add``。
* ``CE`` 的类带入链（甾醇 27:1，``inject``）、介质缩写（``data/oxylipin.csv``）展开、``0:0`` 占位、链数核对。

本文件的钩子（只处理 FA 家族的例外）
------------------------------------
``allows_conjugate(cr)``     只有 26 个 ``ST`` 共轭类（6 个原有 + 20 个 ``ST-<AA>``）允许链上的共轭记号；``SE 24:1;O4;T``、``CE 18:1;Gly`` 一类（共轭挂在酯类上）
                             -> core 抛 ``unsupported_chain_modifier``（以前会被悄悄丢掉共轭）。
``role_candidates(ctx,role)`` ``ST-oacyl`` 的 ``SE`` 两链写法：按甾醇集合 ``STEROL_CHAINS`` 判甾醇链。``/`` 保持"甾醇居首"（LipidBlast、
                             SwissLipids 200/200）；``_`` 时恰有一个 token 在集合里 -> 它是甾醇（MS-DIAL 规范化写法 FA 在前）；
                             两个都在 / 都不在 -> 两种都枚举（L3 空，L2 不受影响）。
``build(lx,cr,cfg,ctx)``     按序做（任一条拒绝即抛 ``Reject``）：
    1. 共轭一致性：链上共轭记号必须与类相符（``ST-Tau`` + ``;Gly`` -> ``class_chain_conflict``）；显式共轭 + 原生类头 = LM 写法，补上
       ``implied_O``（native 的幂等豁免不适用）。
    2. ``WD`` 单个链组（占位链不计）= 总量（参考表 ``WD`` 的 L2/L3 为空，unresolved）-> ``sum_mode``。
    3. ``CE`` + 单条 FA 链 C>=40 且带氧（``CE 50:3;O2`` 一类，参考表 4 名标 ``ST-oacyl2`` 却违反 CE 规则）-> ``unsupported_subclass``。
    4. 多链 ``ST``：只有"氧已完整表达"（有类型化词，或氧在尾括号 ``;(On)``，或没有氧）才收（数据集原生 L2/L3 写法与 GT 的两链 ST）；
       链上带未分型 ``;On``（LipidBlast / LipidDetective 的 ``ST 24:1;O4/13:1`` 胆汁酸酰基写法）-> ``unsupported_subclass``；
       超过 2 个链组 -> ``chain_count``。
    5. ``SE`` 总量+FA 写法 ``SE 45:2_18:1``（benchmark ``normalize_sterol_ester`` 读成甾醇 27:1 + FA 18:1；字面读法 45:2 不是甾醇）：
       读成 27:1 要凭空造出名字里没写的链（BRIEF 判据 3、检查器 P1/B4/B7 都不允许），字面读法的 L1 又可能是错的 -> 拒绝
       ``unsupported_syntax``（判据：较大的 token 不在甾醇集合，且两 token 之差恰为 27:1；``/`` 与 ``_`` 同样处理，B4 元变换一致）。
    6. 类隐含羟基被名字写成裸 ``;O``：``max`` 型 oacyl 类（FAHFA / ``FA-*-oacyl``）``/`` 形式、名字的氧恰好只有第 2 链（羟基链）上的
       一个裸 ``;O``（LipidIN 链级 ``FAHFA 14:0/18:1;O``）-> 它就是类隐含的酯化 OH，记为该链 ``OH``（``class_implied_OH``）。
       ``_`` 看不出羟基链、``;O`` 落在第 1 链或氧不止这一个 -> 不动。
    7. 名字级类型化方言：名字里**任一**链有类型化氧词时，其余链上的裸 ``;On`` 记为 ``O`` 型（GT ``ST 30:5;OH2;oxo2;COOH2/2:0;O2``）；
       只有裸 ``;On`` 的名字不动（C34）。

对外接口
--------
``Handler``（``family = "fa"``）、``FLAGS_fa``（空）、常量 ``STEROL_CHAINS CHOLESTEROL CE_MAX_C_WITH_O CONJUGATE_CLASSES MAX_ST_CHAINS NAME_LEVEL_OTHER_OXYGEN``。
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Final

from ...config import Config
from ...contract import register_family_flags
from ...types import AA_CONJ_NAMES, TOKEN_VALENCE, ClassRes, Lexed, RawChain, Reject, Struct
from .. import conventions as cv
from . import BaseHandler

if TYPE_CHECKING:  # pragma: no cover
    from ..model import Ctx

__all__ = ["Handler", "FLAGS_fa", "STEROL_CHAINS", "CHOLESTEROL", "CE_MAX_C_WITH_O", "CONJUGATE_CLASSES", "MAX_ST_CHAINS",
           "NAME_LEVEL_OTHER_OXYGEN"]

#: 本家族没有自带 flag（全局词表够用）；保留登记以便将来扩展
FLAGS_fa: frozenset[str] = frozenset()
register_family_flags("fa", FLAGS_fa)

#: 常见植物/动物甾醇的 (C, DB) 集合（fam_fa §2：胆固醇 27:1、菜油甾醇 28:1、菜籽甾醇 28:2、谷甾醇 29:1、豆甾醇 29:2）
STEROL_CHAINS: Final = frozenset({(27, 1), (28, 1), (28, 2), (29, 1), (29, 2)})
#: 胆固醇骨架 (C, DB)；``CE`` 的类带入链与总量+FA 写法都用它
CHOLESTEROL: Final = (27, 1)
#: ``CE`` 的单条 FA 链带氧时 C 的上限（含）：>= 此值 -> 拒绝（fam_fa §7：LM 的 ``CE n:m;O2`` n=46..52 共 4 名，参考表标 ST-oacyl2）
CE_MAX_C_WITH_O: Final = 40
#: ``ST`` 多链（类型化方言）时允许的最多链组数（甾体链 + 一个小酰基；GT 只有 2 个）
MAX_ST_CHAINS: Final = 2
#: 名字级类型化方言开关（模块文档第 7 条）。检查器 P3 目前不认 ``O`` 型氧（同链的也不认），关掉它只会让 GT 的
#: ``ST 30:5;OH2;oxo2;COOH2/2:0;O2`` 少一层 L3（core 同链规则不受影响）
NAME_LEVEL_OTHER_OXYGEN: Final = True
#: 允许链上共轭记号的类 -> 与之相符的 ``Features.tail_conj`` 取值（D17：双共轭 ``G+S`` / ``T+S``）
CONJUGATE_CLASSES: Final = {
    "ST-Gly": frozenset({"Gly", "G"}),
    "ST-Tau": frozenset({"Tau", "T"}),
    "ST-sulf": frozenset({"S"}),
    "ST-Hex": frozenset({"Hex"}),
    "ST-Gly-sulf": frozenset({"G+S"}),
    "ST-Tau-sulf": frozenset({"T+S"}),
    # 0.4.0：词表 v0.2.11 的胆汁酸-氨基酸结合物 ST-<AA>（LIPID MAPS ``ST 24:1;O4;Phe``）；Leu / Ile 是旧拼写，同归 ST-Xle
    **{f"ST-{aa}": frozenset({aa}) for aa in AA_CONJ_NAMES if aa not in ("Ile", "Leu", "Xle")},
    "ST-Xle": frozenset({"Xle", "Ile", "Leu"}),
}

_TYPED_KINDS: Final = frozenset({"OH", "oxo", "Ep", "COOH", "OOH", "OO", "OH_pos", "O_other"})
_STEROL_ESTER: Final = "ST-oacyl"


def _chain_atoms(lx: Lexed):
    for rc in lx.chains:
        if not rc.placeholder:
            yield from rc.atoms


def _name_has_typed_oxygen(lx: Lexed) -> bool:
    """名字里任一主链写了类型化氧词（``OH oxo Ep COOH OOH OO``、``;(2OH)``、``O``）。Riken 展开的链不算（它不是名字里的词）。"""
    return any(a.kind in _TYPED_KINDS for a in _chain_atoms(lx))


def _chain_untyped_oxygen(rc: RawChain) -> int:
    return sum(a.n for a in rc.atoms if a.kind == "O_count")


def _name_oxygen(lx: Lexed) -> int:
    """名字字面写出的氧总数（链上 + 尾部记号；不含类隐含氧、不含碱基/类后缀氧）。"""
    n = 0
    for a in (*_chain_atoms(lx), *lx.tail_atoms):
        if a.kind in ("O_count", "OH_pos", "O_other"):
            n += a.n
        elif a.kind in TOKEN_VALENCE:
            n += TOKEN_VALENCE[a.kind] * a.n
    return n


class Handler(BaseHandler):
    family = "fa"

    # ------------------------------------------------------------------ 共轭
    def allows_conjugate(self, cr: ClassRes) -> bool:
        return cr.class_full in CONJUGATE_CLASSES

    # ------------------------------------------------------------------ 角色链
    def role_candidates(self, ctx: "Ctx", role: str) -> list[int] | None:
        """``SE`` 的甾醇链：``/`` 甾醇居首；``_`` 按甾醇集合判（见模块文档）。其余情形走默认。"""
        if role != "sterol" or ctx.cr.class_full != _STEROL_ESTER or ctx.injected or len(ctx.works) != 2:
            return None
        is_sterol = [(w.C, w.DB) in STEROL_CHAINS for w in ctx.works]
        if "_" not in ctx.separators:                                   # '/'：甾醇居首（实测 200/200）
            if not is_sterol[0] and is_sterol[1]:                       # 与甾醇集合矛盾 -> 不下判断
                return [0, 1]
            return [0]
        if is_sterol[0] != is_sterol[1]:                                # '_'：恰有一个像甾醇
            return [is_sterol.index(True)]
        return None                                                     # 都像 / 都不像 -> 默认（两种都枚举）

    # ------------------------------------------------------------------ 装配前的例外
    def build(self, lx: Lexed, cr: ClassRes, cfg: Config, ctx: "Ctx") -> Struct | list[Struct] | None:
        cf = cr.class_full
        self._check_conjugate(lx, cr, ctx)
        if cf == "WD" and ctx.n_groups == 1:
            ctx.sum_mode = True                                         # WD：参考表 L2/L3 为空（unresolved），单 token 只是总量
        if cf == _STEROL_ESTER:
            self._check_cholesteryl_ester(lx, cr, ctx)
            self._reject_sterol_total_notation(lx, ctx)
        elif cf == "ST":
            self._check_multichain_sterol(lx, ctx)
        self._implied_oxygen_written_as_bare_o(lx, cr, cfg, ctx)
        self._name_level_other_oxygen(cr, ctx)
        return None

    # ---- 1 共轭
    @staticmethod
    def _check_conjugate(lx: Lexed, cr: ClassRes, ctx: "Ctx") -> None:
        tc = lx.feats.tail_conj
        if not tc:
            return
        want = CONJUGATE_CLASSES.get(cr.class_full)
        if want is None or tc not in want:
            raise Reject("class_chain_conflict", f"conjugate ;{tc.replace('+', ';')} does not match class {cr.class_full}")
        if ctx.native and cr.implied_mode == "add" and cr.implied_O:    # 显式共轭 = LM 写法：隐含氧要加（native 的幂等豁免不适用）
            ctx.tail_o += cr.implied_O
            ctx.flags.add("implied_O")

    # ---- 3 CE
    @staticmethod
    def _check_cholesteryl_ester(lx: Lexed, cr: ClassRes, ctx: "Ctx") -> None:
        if not ctx.injected or len(ctx.works) != 1:
            return
        w = ctx.works[0]
        if w.C >= CE_MAX_C_WITH_O and _name_oxygen(lx) > 0:
            raise Reject("unsupported_subclass",
                         f"{cr.class_full} fatty-acyl chain with C={w.C} >= {CE_MAX_C_WITH_O} and oxygen "
                         "(not one acyl chain; the reference files these as the three-chain sterol ester class)")

    # ---- 4 多链 ST
    @staticmethod
    def _check_multichain_sterol(lx: Lexed, ctx: "Ctx") -> None:
        if ctx.n_groups < 2:
            return
        if ctx.n_groups > MAX_ST_CHAINS:
            raise Reject("chain_count", f"ST takes one steroid chain (plus at most one small acyl): {ctx.n_groups} given")
        if _name_has_typed_oxygen(lx):
            return
        if any(_chain_untyped_oxygen(rc) for rc in lx.chains if not rc.placeholder) or (
                lx.nested_fa or lx.pre_acyl is not None):
            raise Reject("unsupported_subclass",
                         "ST with an acyl chain and chain-attached untyped oxygen: bile-acid / sterol acyl conjugate "
                         "(the dropped acyl classes of the reference)")

    # ---- 5 SE 总量+FA
    @staticmethod
    def _reject_sterol_total_notation(lx: Lexed, ctx: "Ctx") -> None:
        if lx.head != "SE" or ctx.injected or len(ctx.works) != 2:
            return
        a, b = ctx.works
        big, small = (a, b) if (a.C, a.DB) >= (b.C, b.DB) else (b, a)
        if (big.C, big.DB) in STEROL_CHAINS:
            return
        if (big.C - small.C, big.DB - small.DB) == CHOLESTEROL:
            raise Reject("unsupported_syntax",
                         f"'SE {big.C}:{big.DB} + {small.C}:{small.DB}' reads as the total (cholesterol skeleton plus one fatty "
                         "acid); the total notation is not supported and the literal reading has no sterol chain")

    # ---- 6 类隐含羟基被写成裸 ;O
    @staticmethod
    def _implied_oxygen_written_as_bare_o(lx: Lexed, cr: ClassRes, cfg: Config, ctx: "Ctx") -> None:
        if cr.implied_mode != "max" or cr.oxy_policy != "class_OH" or cr.slash_role != "hfa_second":
            return
        if not cfg.class_implied_OH_in_L3 or len(ctx.works) != 2 or ctx.separators != ("/",):
            return
        w = ctx.works[1]
        if w.raw is None or any(a.kind in _TYPED_KINDS for a in w.raw.atoms):
            return
        n = _chain_untyped_oxygen(w.raw)
        if n != cr.implied_O or _name_oxygen(lx) != n or ctx.tail_o < n or w.toks:
            return
        w.toks["OH"] = n
        ctx.tail_o -= n
        ctx.flags.add("class_implied_OH")

    # ---- 7 名字级方言
    @staticmethod
    def _name_level_other_oxygen(cr: ClassRes, ctx: "Ctx") -> None:
        if not NAME_LEVEL_OTHER_OXYGEN or cr.oxy_policy == "sphingo_OH" or not _name_has_typed_oxygen(ctx.lx):
            return
        for w in ctx.works:
            rc = w.raw
            if rc is None or rc.oxylipin or any(a.kind in _TYPED_KINDS for a in rc.atoms):
                continue
            n = _chain_untyped_oxygen(rc)
            if n and not w.toks.get("O") and n <= cv.L3_TOKEN_MAX["O"] and ctx.tail_o >= n and not any(
                    TOKEN_VALENCE[k] > 0 for k in w.toks):
                w.toks["O"] = n
                ctx.tail_o -= n
