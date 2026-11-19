"""家族 Handler 注册表（M6f 分派）。

``handler(family) -> BaseHandler``     按 ``ClassRes.family``（``gpl gl sl fa``）取（首次使用时才 import 家族模块，线程安全）
``register(family, handler)``          测试/扩展用：覆盖某家族的 Handler（返回旧的，可还原）
``unregister(family)``                 去掉覆盖，下次重新按模块加载
``FAMILIES``                           ``types.FAMILIES``

家族模块的约定：``lipid_name_map/semantics/families/<fam>.py`` 里有 ``class Handler(BaseHandler)``；模块导入时用
``contract.register_family_flags("<fam>", FLAGS_<fam>)`` 登记自带 flags。未知家族 -> ``KeyError``。
"""
from __future__ import annotations

import importlib
import threading

from ..types import FAMILIES
from .families import BaseHandler

__all__ = ["handler", "register", "unregister", "FAMILIES"]

_LOCK = threading.Lock()
_HANDLERS: dict[str, BaseHandler] = {}


def handler(family: str) -> BaseHandler:
    """``family`` 的 Handler 单例。"""
    h = _HANDLERS.get(family)
    if h is not None:
        return h
    if family not in FAMILIES:
        raise KeyError(f"unknown family {family!r}")
    with _LOCK:
        h = _HANDLERS.get(family)
        if h is None:
            mod = importlib.import_module(f"{__package__}.families.{family}")
            h = mod.Handler()
            _HANDLERS[family] = h
    return h


def register(family: str, h: BaseHandler) -> BaseHandler | None:
    """覆盖 ``family`` 的 Handler；返回被覆盖的（没有则 None）。"""
    if family not in FAMILIES:
        raise KeyError(f"unknown family {family!r}")
    with _LOCK:
        old = _HANDLERS.get(family)
        _HANDLERS[family] = h
    return old


def unregister(family: str) -> None:
    with _LOCK:
        _HANDLERS.pop(family, None)
