"""M6 语义层：把词法结果 + 选定的类变成"结构假设"（``Interp``），再交给 M7 造标签。

对外接口
--------
``interpret(lx, cr, cfg, source=None, hi=None) -> Interp``     见 ``core``
``conventions``     纯函数：醚种类 / P 的 DB+1（唯一实现点）/ 碱基氧 / 原子->token / 价数 / 介质缩写字典
``registry``        ``handler(family)`` 按 ``ClassRes.family`` 分派家族钩子
``families``        ``BaseHandler`` 与四个家族模块 ``gpl gl sl fa``（B1-B4 各自负责）
``model``           ``Ctx`` / ``WChain``（家族钩子的工作对象）
"""
from __future__ import annotations

from .core import interpret

__all__ = ["interpret"]
