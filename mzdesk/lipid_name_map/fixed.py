"""固定结构类：名字本身就给出整个分子（目前只有辅酶 Q，ubiquinone）。在每段流水线的最前面处理（``api._run_segment``，M3 之前）。

CoQ 不在 v0.2.11 词表里（该词表只收脂质；``alias/tools.csv`` 的 ``CoQ`` / ``CoQ1``–``CoQ13`` 行是 ``class_not_in_vocab``）。数据集 v0.2.12
的 PR_PK 词表有类 ``CoQ``，三层标签相同：``CoQ <9+5n>:<n+2>;O4``——CoQn = C(9+5n)H(10+8n)O4，DB = 侧链 n 个 C=C + 醌环 2 个 C=C
（数据集 CoQ6 / CoQ8 / CoQ9 / CoQ10 = ``39:8`` / ``49:10`` / ``54:11`` / ``59:12``，逐个核对过）。2026-10-20 起按这个写法收，出到 L3
（名字给出了整个结构，与数据集一致）；类 ``CoQ`` 记在 ``contract.EXTRA_CLASSES``。

对外接口
--------
``fixed_labels(text) -> tuple[str, Labels] | None``    一段文本（已预清洗、装饰与加合物已剥）-> ``(class_full, 标签)``；不归本模块 -> None
``COQ_N``                                            认的异戊二烯单元数 n（1..13，与别名表原来列出的 ``CoQ1``–``CoQ13`` 一致）

认的写法
    ``CoQn``           MS-DIAL / LipidBlast / LipidIN（``CoQ6 | Chain: CoQ6``）/ LipiDetective（``CoQ9 [M+H]+``）
    ``CoQ C:DB;O4``    数据集原生写法（回灌幂等；MS-DIAL 式 ``;4O`` 同）；C、DB 必须对得上某个 n
其余（``CoQ`` 单独、``CoQ 10:0``、``CoQ 39:8``、``CoQ14``、LipidMatch 式 ``Co(Q10)``）返回 None，交回流水线，照旧被拒
（``class_only`` / ``class_not_in_vocab``）。
"""
from __future__ import annotations

import re
from typing import Final

from .types import Labels

__all__ = ["fixed_labels", "COQ_N"]

COQ_N: Final = range(1, 14)
_COQ_NAME: Final = re.compile(r"CoQ(?P<n>[1-9][0-9]?)")
_COQ_LABEL: Final = re.compile(r"CoQ (?P<C>[1-9][0-9]*):(?P<DB>[1-9][0-9]*);(?:O4|4O)")


def _coq(n: int) -> tuple[str, Labels]:
    c, db = 9 + 5 * n, n + 2
    lab = f"CoQ {c}:{db};O4"
    return "CoQ", Labels(lab, lab, lab, 3, sum_meta={"C": c, "DB": db, "O": 4})


def fixed_labels(text: str) -> tuple[str, Labels] | None:
    """见模块文档。"""
    m = _COQ_NAME.fullmatch(text)
    if m:
        n = int(m["n"])
        return _coq(n) if n in COQ_N else None
    m = _COQ_LABEL.fullmatch(text)
    if m:
        n = int(m["DB"]) - 2
        if n in COQ_N and int(m["C"]) == 9 + 5 * n:
            return _coq(n)
    return None
