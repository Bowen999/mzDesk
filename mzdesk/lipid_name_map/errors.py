"""错误类型的统一出口（只做 re-export）。

``Reject``      解析过程中的拒绝/不支持信号：``raise Reject("chain_count", "TG needs 3 chains")``；
                status 由 ``reasons`` 查表（调用方只给 code 与 detail）。
``status_of``   原因码 -> status。
``REASONS``     原因码注册表。
"""
from __future__ import annotations

from .reasons import REASONS, status_of
from .types import Reject

__all__ = ["Reject", "REASONS", "status_of"]
