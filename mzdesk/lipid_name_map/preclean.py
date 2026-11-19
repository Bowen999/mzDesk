"""M1 预清洗（SYNTHESIS §2.2 M1）：把任意输入变成"可交给 M2 的单行 ASCII 化文本"。

对外接口
--------
``preclean(s: object, cfg: Config) -> str``

做什么（按顺序）
    1. ``s`` 不是 ``str`` -> ``Reject("not_string")``（``bytes``、``None``、数字、列表……一律不收）。
    2. 原始长度 > ``cfg.max_len`` -> ``Reject("too_long")``。**先于**任何 Unicode 处理，避免对超长输入白做 NFKC。
    3. （仅非 ASCII 串）NFKC；破折号/减号族 ``‐ ‑ ‒ – — ― − ﹘ ﹣ －`` -> ``-``；除号/分数线 ``∕ ⁄ ／`` -> ``/``；
       比号 ``∶`` -> ``:``；Unicode 格式字符（类别 Cf：零宽空格/连接符、BOM、软连字符、双向控制符……）-> 空格；
       NBSP 等各种空白由 NFKC 或 ``str.split`` 归为普通空格。上述步骤迭代到不再变化（NFKC 本身可能产生新的破折号，
       例如 ``﹘`` -> ``—``），因此 ``preclean`` 是幂等的。
    4. 折叠所有空白串为单个空格并 strip。
    5. 结果为空串 -> ``Reject("empty")``；归一化后长度仍 > ``max_len``（NFKC 可展开，如 ``ﷺ``）-> ``Reject("too_long")``。

不做什么
    * **不改大小写**（``pc`` 与 ``PC`` 对类名查表是不同的输入，大小写敏感由 M4 负责）。
    * 不删/不移动任何可见的 ASCII 字符（``( ) [ ] / _ ; : | + -`` 原样保留，剥装饰是 M2 的事）。
    * 零宽字符变成**空格**而不是被删除：夹在记号中间的零宽字符因此会让后续词法拒绝该名字，而不是悄悄拼出另一个名字。
    * 其余控制字符（``\\x00`` 等）原样保留，由 M3 的 ``bad_chars`` 拒绝。

纯函数、线程安全、确定性。
"""
from __future__ import annotations

import unicodedata

from .config import Config
from .types import Reject

__all__ = ["preclean"]

#: NFKC 之后仍需手工归一的字符（NFKC 不处理，或把它们换成另一个破折号/斜线）
_DASHES = "‐‑‒–—―−⁃﹘﹣－⁻₋"
_SLASHES = "⁄∕⧸／"
_TRANSLATE = {
    **{ord(c): "-" for c in _DASHES},
    **{ord(c): "/" for c in _SLASHES},
    ord("∶"): ":",
}

_MAX_ROUNDS = 4  # NFKC -> translate 的不动点迭代上限（实际 1-2 轮即稳定）


def _normalize_unicode(t: str) -> str:
    """非 ASCII 串的 Unicode 归一（不含空白折叠）。迭代到不动点。"""
    for _ in range(_MAX_ROUNDS):
        u = unicodedata.normalize("NFKC", t).translate(_TRANSLATE)
        if not u.isascii() and any(unicodedata.category(c) == "Cf" for c in u):
            u = "".join(" " if unicodedata.category(c) == "Cf" else c for c in u)
        if u == t:
            return u
        t = u
    return t


def preclean(s: object, cfg: Config) -> str:
    """见模块文档。返回折叠空白后的 ``str``（从不返回 ``str`` 子类实例）。

    >>> from lipid_name_map.config import DEFAULT
    >>> preclean("ＰＣ　１６：０／１８：１", DEFAULT)
    'PC 16:0/18:1'
    >>> preclean("PC\\u00a016:0\\u2212", DEFAULT)
    'PC 16:0-'
    """
    if not isinstance(s, str):
        raise Reject("not_string", f"expected str, got {type(s).__name__}")
    max_len = cfg.max_len
    if len(s) > max_len:
        raise Reject("too_long", f"length {len(s)} > max_len {max_len}")
    t = s if s.isascii() else _normalize_unicode(s)
    t = " ".join(t.split())
    if not t:
        raise Reject("empty", "input is empty or whitespace only")
    if len(t) > max_len:
        raise Reject("too_long", f"normalized length {len(t)} > max_len {max_len}")
    return t
