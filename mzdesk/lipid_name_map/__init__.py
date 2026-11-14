"""lipid_name_map：脂质 shorthand 名 -> 标准化 L1 / L2 / L3（benchmark 去 sn 写法）。

对外接口
--------
``normalize(name, *, source=None, config=None) -> dict``     单个名字；永不抛异常；返回字段见 ``contract``（0.3.0 起含 LIPID MAPS 名 ``LM_species`` / ``LM_molecular``，见 ``lm``）
``normalize_many(names, *, source=None, config=None) -> list[dict]``
``Config`` ``DEFAULT``                                        运行配置（``config.py``）
``__version__``

``normalize`` / ``normalize_many`` 在首次调用时才 import ``.api``（W3 实现）。``.api`` 尚不存在时，调用抛
``NotImplementedError``，但 ``import lipid_name_map`` 本身总是成功。

契约与数据结构：``contract``（返回字典构造/检查）、``types``（共享 dataclass）、``reasons``（原因码）、
``errors``（``Reject``）、``config``。详见 ``docs/CONTRACT.md``。
"""
from __future__ import annotations

import importlib
from collections.abc import Iterable
from typing import Any

from .config import DEFAULT, Config

__version__ = "0.4.0"

__all__ = ["normalize", "normalize_many", "Config", "DEFAULT", "__version__"]

_api: Any = None


def _load_api() -> Any:
    """惰性 import ``.api``；它不存在时抛 ``NotImplementedError``（api 内部缺依赖仍按原样抛 ImportError）。"""
    global _api
    if _api is None:
        try:
            _api = importlib.import_module(".api", __name__)
        except ModuleNotFoundError as e:
            if e.name == f"{__name__}.api":
                raise NotImplementedError("lipid_name_map.api is not implemented yet (owned by W3)") from None
            raise
    return _api


def normalize(name: object, *, source: str | None = None, config: Config | None = None) -> dict:
    """把一个脂质 shorthand 名标准化；见 ``lipid_name_map.api.normalize``（W3）。"""
    return _load_api().normalize(name, source=source, config=config)


def normalize_many(names: Iterable[object], *, source: str | None = None, config: Config | None = None) -> list[dict]:
    """批量版 ``normalize``；见 ``lipid_name_map.api.normalize_many``（W3）。"""
    return _load_api().normalize_many(names, source=source, config=config)
