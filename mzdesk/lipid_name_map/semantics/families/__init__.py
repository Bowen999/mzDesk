"""家族钩子基类 ``BaseHandler``（M6f）。四个家族文件 ``gpl gl sl fa`` 各自 ``class Handler(BaseHandler)``。

core（``semantics/core.py``）把绝大多数规则做成**表驱动的通用装配**（``ClassRes`` 的 ``nchain implied_O implied_mode slash_role
oxy_policy pnnl_O inject`` 字段），家族钩子只处理例外。钩子的调用点与权限：

``select_hint(hi, lx, cr, cfg, source) -> ClassRes | None``
    D13：core 在装配前先调。返回非 None 就用它**替换** ``cr``（例：``source="pnnl"`` 的多链类单 C:DB -> 溶血类；
    ``Cer-d`` 单链 -> ``SPB``；``Cer(d18:0/0:0)`` -> ``SPB``）。新类用 ``classes.class_info(cf).to_class_res(cr.alias_status, cr.inject)``
    造（可 import ``classes`` 的公开 API：``class_info class_sibling head_ether_marker``，D15）。
    返回的家族可以与原来不同（core 会改用新类家族的 Handler 装配）。默认 None。
``allows_conjugate(cr) -> bool``
    链上出现共轭记号（``;Gly ;Tau ;S ;Hex ;G ;T``）时，core 问这个。默认 False -> ``unsupported_chain_modifier``。
    ``fa`` 家族对 ``ST`` 核返回 True。
``allows_nested_acyl(cr) -> bool``
    名字带嵌套 ``(FA c:d)`` / 括号 ``(O-c:d)`` 酰基前缀时，core 问这个。默认：类的 ``mod`` 里有 ``oacyl`` 段（``Cer-oacyl`` ``SM-oacyl`` ``TG-oacyl`` …）才允许，
    否则 ``unsupported_syntax``；``gl`` 对 ``TG-oacyl2`` 也返回 True（两个嵌套酰基，0.4.0）。
``allows_adjacent_groups(cr) -> bool``
    名字写成相邻括号组 ``CL(28:0)(28:0)``（词法当 ``_`` 连接的链列表，``conventions.adjacent_groups``）时，core 问这个。默认 False ->
    ``malformed``（这种写法只有 CL 对记法有定义，D24）；``gpl`` 对 ``CL`` 返回 True。
``role_candidates(ctx, role) -> list[int] | None``
    角色链（``"base"`` ``"hfa"`` ``"sterol"``）的候选下标（指 ``ctx.all_chains()``）。None = 用默认
    （``model.default_role_candidates``）。例：SE 的 ``_`` 写法按甾醇集合 ``{27:1, 28:1, 28:2, 29:1, 29:2}`` 判甾醇。
``build(lx, cr, cfg, ctx) -> Struct | list[Struct] | None``
    core 构造完 ``Ctx`` 之后、做任何"判决"（链数、缺氧数、O 不可表达）之前调用。
    * 返回 ``None``：core 用默认装配（钩子可先改 ``ctx``：``ctx.works[i].toks``、``ctx.tail_o``、``ctx.pending``、``ctx.flags``、
      ``ctx.sum_mode``、``ctx.ether_linkage`` …）。
    * 返回 ``Struct`` / ``list[Struct]``：直接作为假设（1..6 个，``class_full`` 必须一致）；``ctx.flags`` 与 ``ctx.ether_linkage``
      仍会进 ``Interp``。钩子自己负责链数、氧数的合法性（抛 ``Reject``）。
    * 想在默认结果上改：``structs = ctx.default_structs()``，再 ``dataclasses.replace``。

家族可以声明 ``FLAGS_<fam> = frozenset({...})`` 并 ``contract.register_family_flags("<fam>", FLAGS_<fam>)``（模块导入时）。
钩子**不得** import 其它家族模块或 core 内部函数。

配方（``ctx`` 的字段见 ``semantics/model.py``；``conventions`` 是公开的纯函数集合，可以 import）
------------------------------------------------------------------------------------------
* **换类**（PNNL 溶血、``Cer-d`` 单链 -> ``SPB``）：``select_hint`` 里返回 ``class_info(cf).to_class_res(cr.alias_status, cr.inject)``；
  之后 core 按新类装配（``PNNL`` 类后缀 ``-d`` 的碱基羟基在候选唯一时直接落到那条链上）。
* **强制总量并改总数**（LipidBlast ``HexCer a;2O/b;2O``：DB-1、O=3、level 1）：``build`` 里
  ``ctx.flags.add("lipidblast_eos_rule")`` 后返回 ``Struct(cr.class_full, None, total=(C, DB - 1, 3))``。
* **改某条链的氧**（Cer 子类字母 A/B 的酰基 +1 OH）：``build`` 里
  ``w = ctx.works[<酰基链下标>]; w.toks["OH"] = w.toks.get("OH", 0) + 1``，加 flag 后 ``return None``（默认装配会重算总氧）。
* **角色链判别**（SE 的 ``_`` 写法按甾醇集合）：覆盖 ``role_candidates(ctx, "sterol")``，返回 ``ctx.all_chains()`` 的下标列表。
* **在默认结果上改**：``structs = ctx.default_structs()``（会抛 ``chain_count`` 等），再 ``dataclasses.replace``。
* 想拒绝：``raise Reject("<code>", detail)``（码必须在 ``reasons.REASONS``）。
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from ...config import Config
from ...types import ClassRes, HeadInfo, Lexed, Struct

if TYPE_CHECKING:  # pragma: no cover
    from ..model import Ctx

__all__ = ["BaseHandler"]


class BaseHandler:
    """所有钩子的默认实现 = "没有例外"。"""

    family: str = ""

    def select_hint(self, hi: HeadInfo, lx: Lexed, cr: ClassRes, cfg: Config, source: str | None) -> ClassRes | None:
        return None

    def allows_conjugate(self, cr: ClassRes) -> bool:
        return False

    def allows_nested_acyl(self, cr: ClassRes) -> bool:
        return "oacyl" in cr.mod.split("-")

    def allows_adjacent_groups(self, cr: ClassRes) -> bool:
        return False

    def role_candidates(self, ctx: "Ctx", role: str) -> list[int] | None:
        return None

    def build(self, lx: Lexed, cr: ClassRes, cfg: Config, ctx: "Ctx") -> Struct | list[Struct] | None:
        return None
