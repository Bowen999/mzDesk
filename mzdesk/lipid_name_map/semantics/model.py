"""M6 的工作对象：``WChain``（工作链）与 ``Ctx``（一次 ``interpret`` 的上下文）。

``core.interpret`` 先把 ``Lexed`` 加工成 ``Ctx``（丢 ``0:0``、展开介质缩写、算醚种类、按氧策略把原子变成 token、
落类隐含氧与 PNNL 后缀氧……），再调家族钩子 ``Handler.build(lx, cr, cfg, ctx)``：

* 钩子可以**改 ctx**（``works[i].toks``、``tail_o``、``pending``、``flags``、``sum_mode``……）后返回 ``None``，由 core 用默认装配；
* 也可以自己造 ``Struct`` / ``list[Struct]`` 返回（core 只做一致性检查：同类、≤6 个假设、O 总数可表达）。
* 钩子要默认装配结果再改，可以调 ``ctx.default_structs()``。

``Ctx`` 的字段（可变；除 ``lx cr cfg source`` 外都可由钩子修改）
------------------------------------------------------------
lx, cr, cfg, source   输入（只读约定）
works         主链 ``WChain``（书写顺序，``0:0`` 占位已丢，介质缩写已展开）；``kind == "main"``
extras        O-酰基链：``Lexed.pre_acyl``、``Lexed.nested_fa``（``kind == "oacyl"``，无 token，不能当醚链/碱基/角色链）
injected      类带入的链（``ClassRes.inject``，如 CE 的甾醇 ``27:1``；``kind == "inject"``）
n_groups      ``len(works) + len(extras)``（D16 的链数，**不含** injected）
n_inject      ``len(injected)``
sum_mode      True = 名字只有总量（单个链组而类需要更多链）-> ``chains=None``，不拆链
ether_kind    ``"" | "O" | "P"``（``conventions.ether_kind``；``"P"`` 的 DB+1 由 ``label_db`` 在装配时加，且只加一次）
tail_o        未分型氧个数（不在任何链上；L3 空）
pending       待放置的"角色羟基" ``[(role, n)]``：``base`` = PNNL 类后缀 ``-d/-t/-m`` 的碱基羟基（碱基链未知），
              ``sterol`` / ``hfa`` = 类隐含羟基（``class_OH`` 策略）；装配时按 ``role_candidates`` 枚举落点
flags         ``set[str]``，最终进 ``Interp.flags``
ether_linkage ``"P" | "O" | "unspecified" | ""``
handler       本类的家族 Handler
name_o        名字**自己写出**的氧个数（链上 token + 尾部记号，不含类隐含氧 / PNNL 后缀氧）；用于判"显式氧数"：D27 的单 token 多链鞘脂和式
              （``Cer 18:1;O2``）打 ``convention_unverified``，以及 gpl 的 CL 对记守卫
native        True = 输入是数据集原生写法（头部就是 ``class_full``，或别名行的 ``dialects`` 含 ``native``）：标签空间里已含类隐含氧，
              ``implied_mode == "add"`` 的类不再加一次（幂等）
separators    主链之间的分隔符（占位链已去掉）``"/"`` 或 ``"_"``
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

from ..config import Config
from ..types import ClassRes, Lexed, RawChain, Struct

if TYPE_CHECKING:  # pragma: no cover
    from .families import BaseHandler

__all__ = ["WChain", "Ctx"]


@dataclass(slots=True)
class WChain:
    """一条工作链（可变）。

    C, DB      碳数 / 双键数。DB 是**标签口径，但不含待加的 P 加一**（见 ``db_extra``）；介质缩写展开后的 DB 已是标签口径
    toks       已分型 token ``{name: n}``（``TOKEN_ORDER`` 去掉 ``P``，含 ``Me`` ``cPr``）
    untyped    暂存的未分型氧（``Ctx`` 构造完成后恒为 0：已按策略转成 token 或进了 ``Ctx.tail_o``）
    mark       这条链自带的醚种类 ``"" | "O" | "P"``
    db_extra   待加 DB（字面 ``P-`` 前缀为 1；原生 ``;P`` 为 0）；装配时 ``Chain.DB = label_db(DB, db_extra)``
    role       ``"" | "base" | "hfa" | "sterol" | "acyl" | "oacyl" | "inject"``（诊断用；M7 排序不看角色）
    kind       ``"main" | "oacyl" | "inject"``
    n_acyl     词法 ``N-`` 前缀（不能当醚链）
    base       词法碱基前缀 ``"" | "d" | "t" | "m"``
    raw        来源 ``RawChain``（介质展开/类带入时为 None）
    oh_read    ``toks["OH"]`` 里输入显式写出的羟基个数（``;OH`` ``;OH2`` ``;(2OH)`` 与 ``d/t/m`` 换算）；``;O`` / ``;On`` 按鞘脂政策补成的 ``OH``、
               钩子（亚类字母等）后加的 ``OH``、待放置羟基落点都不算。**默认 0 = 全部视为假设**（失败保守）。装配时写进 ``Chain.oh_read``，只有
               ``lm.py`` 的鞘脂碱基推断读它
    """
    C: int
    DB: int
    toks: dict[str, int] = field(default_factory=dict)
    untyped: int = 0
    mark: str = ""
    db_extra: int = 0
    role: str = ""
    kind: str = "main"
    n_acyl: bool = False
    base: str = ""
    raw: RawChain | None = None
    oh_read: int = 0


@dataclass(slots=True)
class Ctx:
    """一次 ``interpret`` 的上下文，见模块文档。"""
    lx: Lexed
    cr: ClassRes
    cfg: Config
    source: str | None
    handler: "BaseHandler"
    works: list[WChain] = field(default_factory=list)
    extras: list[WChain] = field(default_factory=list)
    injected: list[WChain] = field(default_factory=list)
    n_groups: int = 0
    n_inject: int = 0
    sum_mode: bool = False
    ether_kind: str = ""
    tail_o: int = 0
    pending: list[tuple[str, int]] = field(default_factory=list)
    flags: set[str] = field(default_factory=set)
    ether_linkage: str = ""
    native: bool = False
    name_o: int = 0
    separators: tuple[str, ...] = ()
    scratch: dict[str, Any] = field(default_factory=dict)      # 钩子之间/钩子内部传递的草稿区，core 不读

    # ---- 便捷视图 ------------------------------------------------------------------------
    def all_chains(self) -> list[WChain]:
        """``works + extras + injected``（装配顺序；候选下标都指这个列表）。"""
        return [*self.works, *self.extras, *self.injected]

    def role_candidates(self, role: str) -> list[int]:
        """角色 ``role``（``base`` / ``hfa`` / ``sterol``）的候选链下标（指 ``all_chains()``）。

        先问家族钩子 ``Handler.role_candidates``（返回 None = 用默认）；默认：``sterol`` 有带入的甾醇链则只有它；
        否则主链只有一条 -> 它；分隔符全是 ``/`` -> 按 ``slash_role`` 定位（``base_first`` / ``sterol_first`` 首链，
        ``hfa_second`` 第二链）；有 ``_`` -> 所有主链（``n_acyl`` 链除外）。
        """
        got = self.handler.role_candidates(self, role)
        if got is not None:
            return list(got)
        return default_role_candidates(self, role)

    def default_structs(self) -> list[Struct]:
        """core 的默认装配（链数核对 -> 总量/拆链 -> 假设枚举 -> ``Struct``）；会抛 ``chain_count`` 等拒绝。"""
        from . import core
        return core.default_structs(self)


def default_role_candidates(ctx: Ctx, role: str) -> list[int]:
    """见 ``Ctx.role_candidates``。"""
    n = len(ctx.works)
    if role == "sterol" and ctx.injected:
        for j, w in enumerate(ctx.injected):
            if w.role == "sterol":
                return [n + len(ctx.extras) + j]
    if n <= 1:
        return [0] if n else []
    ordered = "_" not in ctx.separators
    if ordered:
        slash = ctx.cr.slash_role
        if slash in ("base_first", "sterol_first"):
            return [0]
        if slash == "hfa_second":
            return [1]
    cands = [i for i, w in enumerate(ctx.works) if not w.n_acyl]
    return cands or list(range(n))
