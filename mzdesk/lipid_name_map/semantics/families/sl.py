"""鞘脂与糖鞘脂（SL）家族钩子（B3，``docs/census/fam_sl.md``）。

core（``semantics/core.py``）已经做了的鞘脂通用规则：``d/t/m`` 碱基氧（S1）、``sphingo_OH`` 策略（``;nO`` -> ``OH``，打
``oxygen_assumed_OH``）、PNNL 类后缀 ``-d/-t/-m``（碱基链位置未知的假设枚举，``ClassRes.pnnl_O`` 超出 2 的部分进 ``tail_o``，S4）、
缺氧数按 d 型补碱基氧（C26 改判，``base_oxygen_assumed``；``Config.assume_base_oxygen=False`` 时拒绝）/ 单记号碱基拒绝（C33）、
链数核对（D16）。本模块只处理**本家族的例外**。

对外接口
--------
``Handler``（``family == "sl"``，注册表按 ``ClassRes.family`` 分派）::

    select_hint(hi, lx, cr, cfg, source) -> ClassRes | None      D13：换类 / 拒绝
    build(lx, cr, cfg, ctx)              -> Struct | None         在 core 默认装配前改 ctx，或直接给总量结构
``FLAGS_sl``    本家族自带 flags（现为空：用到的 flag 都在全局 ``types.FLAGS`` 里）

``select_hint``（先于装配，可换类或拒绝）
    1. **LM 糖链组成式必须服从组成规则**（S9）：``Hex(n)-HexNAc(m)-…-Cer`` 的头，若别名表的显式行指向的类不是
       ``Cer-<组成 key>``（参考表按结构重标的 ``Hex(4)-HexNAc(3)-Fuc-NeuAc-Cer`` 8 行，规则给出的类不在词表），
       -> ``glycan_unknown``（census §9-④ / ``fam_sl_examples`` 的 REJECT 行）。
    2. **游离碱基**：``Cer`` + 单链 + 类后缀 ``-d/-t``（PNNL ``Cer-d 18:1``，与 ``source`` 无关：后缀本身就是 PNNL 写法）-> ``SPB``；
       ``Cer(d18:0/0:0)``（碱基居首，第 2 链是 ``0:0``）-> ``SPB``（S4、census §3）。后缀写法要求 C < ``BASE_ONLY_MAX_C``（26，C33：
       ``Cer-d 42:2`` 是和式，不是 42 碳的碱基），占位写法不限。
    3. **溶血鞘磷脂**：``SM`` + 单链（C < 26）+ 类后缀 ``-d/-t``，或 ``SM(d18:1/0:0)``，或 ``Config.single_token_is_lyso(source)``（D4，
       ``source="pnnl"``）为真的单个 C:DB -> ``SM-lyso``（core 再打 ``source_hint_lyso``）。其它两链类（``GalCer-d 18:1``）的单链后缀写法是
       溶血形式（psychosine 型），词表没有对应类 -> ``class_not_in_vocab``（不当 18 碳的总量）。
    4. **EOS 规则只认 ``a;2O/b;2O`` 的斜线写法**：别名表按 ``eos2O`` 谓词把 ``HexCer`` 归到 ``Cer-Hex-oacyl``。P1 之后词法 ``eos2O`` 本身
       就只在第 1、2 条非占位链之间的分隔符全是 ``/`` 时为真（``HexCer 18:0;O2_16:0;2O`` 的 ``;2O`` 只是氧个数的另一种写法，
       ``eos2O`` 为假，别名表直接给 ``Cer-Hex``；``Cer`` 头的 ``Cer 24:1_18:1;2O`` 同理不再被 ``unsupported_subclass`` 拒绝），
       所以这里不再需要把 ``Cer-Hex-oacyl`` 换回 ``Cer-Hex`` 的补丁（原来的第 4 条钩子已删，``has_slash_pair`` 保留作 ``build`` 的防御性核对）。
``build``（先改 ``ctx`` 再交 core 默认装配；或直接返回总量结构）
    1. **LipidBlast ``HexCer a;2O/b;2O``**（S6、D9）：``Cer-Hex-oacyl``，总量 ``C=ΣC, DB=ΣDB-1, O=3``，level 1，flag ``lipidblast_eos_rule``
       （参考表 3,759/3,759）。
    2. **LipidMatch 亚类字母**（S11，连字符写法 ``Cer-AS`` ``Cer-BS`` ``Cer-ADS`` ``Cer-BDS`` ``Cer-AP`` ``HexCer-AP``）：首字母 ``A/B`` 的酰基
       链比名字多 1 个 OH（打 ``subclass_implied_O``；能唯一定位酰基链时记成链上的 ``OH`` 并打 ``class_implied_OH``，否则进 ``tail_o``，
       L3 留空）；字母 ``P``（植物鞘氨醇）要求碱基是 ``t``，``d/m`` -> ``class_chain_conflict``；名字没给碱基氧、core 按 d 型补了 2 个
       （``base_oxygen_assumed``）时，字母 ``P`` 改补 3 个（t 型）。MS-DIAL 下划线写法 ``Cer_AS`` 的字母
       **仅元数据**（A8：名字自带氧标记，不补氧）。
    （MS-DIAL 规范化名 S16 的 ``harmonized_name_oxygen_placement`` 已移到 core，按氧策略 ``sphingo_OH`` 判——``SL`` 的家族是 ``fa``，
    钩子管不到它，原来这里按家族打会漏。）

不在这里做的（数据或别的模块已处理）
    * ``Cer`` + ``eos2O`` -> ``unsupported_subclass``（D18，``alias/tools.csv`` 数据行）；``Cer-EOS/EODS/EBDS/HS/HDS`` 同（别名表 reject 行）
    * ``dO-`` / ``p`` 链标记、``Ox*`` 前缀（B0 / 别名表）
    * 糖链类名生成：LM 组成式 -> ``data/glycan.csv``（``classes.resolve_head``）；Svennerholm / SwissLipids 俗名 -> ``alias/*.csv``
    * ``SL``（磺基脂）类的家族是 ``fa``（``class_info/fa.csv``，B4），这里的钩子管不到它
"""
from __future__ import annotations

import re
from typing import Final

from ...classes import class_info, class_sibling, parse_glycan_head
from ...config import Config
from ...contract import register_family_flags
from ...types import ClassRes, HeadInfo, Lexed, Reject, Struct
from .. import conventions as cv
from ..model import Ctx
from . import BaseHandler

__all__ = ["Handler", "FLAGS_sl", "SUBCLASS_RE", "BASE_ONLY_MAX_C", "has_slash_pair"]

#: 单个 C:DB 的 C 低于它是"一条游离碱基 / 溶血链"，不低于它是和式（与 core 的 ``SINGLE_TOKEN_BASE_MIN_C`` 同值，C33）
BASE_ONLY_MAX_C: Final = 26

FLAGS_sl: frozenset[str] = frozenset()
register_family_flags("sl", FLAGS_sl)

#: LipidMatch 连字符亚类写法：首字母 N（普通）/A（α-羟基）/B（β-羟基），次字母 S / DS（鞘氨醇 / 二氢鞘氨醇）/ P（植物鞘氨醇）。
#: 带 ``-`` 才是 LipidMatch（A8：``Cer_AS`` 下划线的 MS-DIAL ontology 名只是元数据）。H / E 开头的由别名表 reject（unsupported_subclass）。
SUBCLASS_RE: Final = re.compile(r"^(?P<core>Cer|HexCer)-(?P<acyl>[NAB])(?P<base>DS|S|P)$")
_LYSO_SUFFIXES: Final = ("-d", "-t")


def has_slash_pair(lx: Lexed) -> bool:
    """第 1、2 条非占位主链之间的分隔符全是 ``/``（LipidBlast EOS 的 ``a;2O/b;2O`` 写法）。少于两条链 -> False。"""
    real = [i for i, c in enumerate(lx.chains) if not c.placeholder]
    if len(real) < 2:
        return False
    return all(s == "/" for s in lx.seps[real[0]:real[1]])


def _no_extra_chains(lx: Lexed) -> bool:
    return not lx.nested_fa and lx.pre_acyl is None


def _eos_head(lx: Lexed, cr: ClassRes) -> bool:
    """别名表按 ``eos2O`` 把 ``HexCer`` 归到 ``Cer-Hex-oacyl``（ref_main.csv）：头是 ``HexCer``、无 ``(FA)``/``(O-)`` 嵌套链。
    ``AHexCer`` / ``AHexCer-O`` 本来就是 ``Cer-Hex-oacyl``（3 条链），与 EOS 规则无关。"""
    return (lx.head == "HexCer" and cr.class_full == "Cer-Hex-oacyl" and lx.feats.eos2O and _no_extra_chains(lx))


def _swap(cr: ClassRes, cf: str) -> ClassRes:
    """用 ``class_info(cf)`` 造替换类，沿用别名状态与带入链。"""
    return class_info(cf).to_class_res(cr.alias_status, cr.inject)


class Handler(BaseHandler):
    family = "sl"

    # ------------------------------------------------------------------------------------------
    # select_hint
    # ------------------------------------------------------------------------------------------
    def select_hint(self, hi: HeadInfo, lx: Lexed, cr: ClassRes, cfg: Config, source: str | None) -> ClassRes | None:
        f = lx.feats
        cf = cr.class_full
        head = hi.head if hi is not None else lx.head

        # 1. LM 糖链组成式：显式行与组成规则不一致 -> 拒绝（S9）
        key = parse_glycan_head(head)
        if key is not None and cf != "Cer-" + key:
            raise Reject("glycan_unknown",
                         f"{head}: the composition rule gives Cer-{key}, which is not the class {cf} the table row says")

        single = f.n_chains == 1 and _no_extra_chains(lx)
        # 单个 C:DB 且 C < 26 才是"一条游离链"（C33：C >= 26 是和式 Cer-d 42:2，不是 42 碳的碱基）
        one_chain = single and f.n_placeholder == 0 and lx.chains[0].C < BASE_ONLY_MAX_C
        free_base = single and f.n_placeholder == 1 and not lx.chains[0].placeholder       # d18:0/0:0 ：游离碱基
        suffix_single = one_chain and f.class_suffix in _LYSO_SUFFIXES                     # PNNL Cer-d 18:1

        # 2. 游离鞘氨醇碱基
        if cf == "Cer" and (free_base or suffix_single):
            return _swap(cr, "SPB")

        # 3. 溶血鞘磷脂
        if cf == "SM" and (free_base or suffix_single
                           or (one_chain and not f.class_suffix and cfg.single_token_is_lyso(source))):
            return _swap(cr, "SM-lyso")
        if suffix_single and cr.nchain == frozenset({2}) and class_sibling(cf, lyso=True) is None:
            # PNNL 的 GalCer-d 18:1 之类：单链 = 溶血形式（psychosine 型），词表里没有这些溶血类；不能当 18 碳的 HexCer 总量
            raise Reject("class_not_in_vocab", f"{head}: a single chain with a -d/-t suffix is the lyso form, which has no class in the vocabulary")
        return None

    # ------------------------------------------------------------------------------------------
    # build
    # ------------------------------------------------------------------------------------------
    def build(self, lx: Lexed, cr: ClassRes, cfg: Config, ctx: Ctx) -> Struct | list[Struct] | None:
        cf = cr.class_full

        # 1. LipidBlast HexCer a;2O/b;2O -> Cer-Hex-oacyl 总量（DB-1，O=3，只到 L1）
        if _eos_head(lx, cr) and ctx.n_groups == 2 and has_slash_pair(lx):
            c = sum(w.C for w in ctx.works)
            db = sum(w.DB for w in ctx.works) - 1
            if db < 0:
                raise Reject("malformed", "the LipidBlast EOS rule (DB-1) needs at least one double bond")
            ctx.flags.add("lipidblast_eos_rule")
            ctx.flags.discard("oxygen_assumed_OH")                  # 只有 L1，没有链级 OH，这个假设标志没有意义
            return Struct(cf, None, total=(c, db, 3))

        # 2. LipidMatch 亚类字母
        m = SUBCLASS_RE.match(lx.head)
        if m is not None:
            self._subclass(m, lx, cfg, ctx)

        return None

    # ------------------------------------------------------------------------------------------
    def _subclass(self, m: re.Match[str], lx: Lexed, cfg: Config, ctx: Ctx) -> None:
        """S11：``P`` 要求 ``t``；首字母 ``A/B`` -> 酰基链 +1 OH。"""
        bases = {w.base for w in ctx.works if w.base}
        if m.group("base") == "P" and bases - {"t"}:
            raise Reject("class_chain_conflict",
                         f"subclass {lx.head} (phytosphingosine) needs a t-prefixed base, the name has {sorted(bases)}")
        if m.group("base") == "P" and "base_oxygen_assumed" in ctx.flags:   # core 按 d 型补的碱基氧：字母 P（植物鞘氨醇）是 t 型，补 3 个
            ctx.pending = [(r, cv.base_oxygen("t") if r == "base" else n) for r, n in ctx.pending]
        if m.group("acyl") not in ("A", "B"):
            return
        if _name_oxygen(ctx) == 0:                                   # 字母只补酰基的 1 个 OH，碱基的氧名字里必须有（core 默认已按 d 型补，
                                                                     # 只有 Config.assume_base_oxygen=False 才走到这里）
            raise Reject("missing_oxygen_count", f"{lx.head}: the name gives no oxygen for the sphingoid base (d/t/m prefix or ;On)")
        ctx.flags.add("subclass_implied_O")
        if ctx.sum_mode:                                             # 只有总量：氧加在总量上
            ctx.tail_o += 1
            return
        acyl = _acyl_candidates(ctx)
        if len(acyl) != 1 or not cfg.class_implied_OH_in_L3:         # 酰基链定位不了 / 配置关掉类隐含 OH 出 L3：未分型氧（L3 留空）
            ctx.tail_o += 1
            return
        w = ctx.works[acyl[0]]
        if cv.valence(w.toks) >= 1:                                  # 名字已经给了这条链的氧：字母不再补
            ctx.flags.discard("subclass_implied_O")
            return
        w.toks["OH"] = w.toks.get("OH", 0) + 1
        ctx.flags.add("class_implied_OH")


def _name_oxygen(ctx: Ctx) -> int:
    """名字自己给出的氧个数（链上 token + 尾部未分型氧 + 待放置羟基）。"""
    return sum(cv.valence(w.toks) for w in ctx.all_chains()) + ctx.tail_o + sum(n for _, n in ctx.pending)


def _acyl_candidates(ctx: Ctx) -> list[int]:
    """酰基链在 ``ctx.works`` 里的下标：有碱基前缀的链之外的链；没有任何前缀时按 ``/`` 的角色（碱基居首，第 2 条是酰基）。"""
    works = ctx.works
    if any(w.base for w in works):
        return [i for i, w in enumerate(works) if not w.base and not w.n_acyl]
    if len(works) == 2 and "_" not in ctx.separators:
        return [1]
    return []
