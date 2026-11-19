"""M7：从解析结果 ``Interp`` 构造去 sn 的 L1 / L2 / L3（SYNTHESIS §4.1–4.2，依据 ``label_rules`` §4、§5）。

对外接口
--------
``build_labels(it: Interp, cfg: Config) -> Labels``
    对 ``it.alts`` 的每个假设（``Struct``，≤6 个）各算一遍 L1/L2/L3，再合并（见下），``level = 3 if L3 else 2 if L2 else 1``。
    ``cfg`` 目前不影响标签字符串（``Config.emit_mecpr_field`` 只影响 ``l3_with_mecpr``）。
``l3_with_mecpr(it: Interp, cfg: Config) -> str``
    ``meta.L3_with_MeCpr`` 的值：仅当 ``cfg.emit_mecpr_field`` 且默认 L3 非空且某链带 ``Me``/``cPr`` 时非空，
    否则 ``""``。``Labels`` 没有放它的字段，由 ``api`` 调用本函数后写进 ``meta``。
``struct_totals(st: Struct) -> (C, DB, O)``   一个假设的总量（``chains is None`` 取 ``st.total``）
``nominal_chain_counts(it: Interp, class_full: str) -> frozenset[int]``  名义链数集合；空集 = 不检查（ST、WD）
``fmt_l1_o(n) / fmt_l2_tail(n)``              L1 的 ``;O``/``;On`` 与 L2 的 ``;(O)``/``;(On)`` 尾缀（n=0 为 ``""``，永不 ``;O1``）

构造算法（每个假设）
--------------------
1. 总量 ``C = ΣC, DB = ΣDB, O = Σ氧价(toks) + tail_o``（``chains is None`` 取 ``total``）。
2. ``L1 = "<cf> C:DB" + ;O/;On``。
3. ``chains is None``，或链数不在名义链数集合内（非空时） -> ``L2 = L3 = ""``。
4. 排序：L2、L3 共用一次 ``nosn.order_chains``（键 ``(C, DB, token种类数, 首token首字母, 链全串)`` 全降序，``O-`` 不置首）。
5. ``L2 = "<cf> C:DB_C:DB…" + ;(O)/;(On)``（至多一组尾括号，链上无 token、无 ``O-``）。
6. ``L3`` 门槛：``tail_o == 0``；醚类（``class_full`` 含 ``ether`` 段）恰一条醚链、非醚类零条；否则 ``""``。
   通过则 ``L3 = "<cf> " + "_".join(render_chain)``。``Me``/``cPr`` 默认丢弃。

假设合并（D1）
--------------
* L1 必须在所有假设上相同，否则 ``Reject("internal_error")``（semantics 给出了总量冲突的假设 = bug）。
* L2 只有所有假设的 L2 串相同（且非空）才输出；L3 同理（每个假设都过门槛且串相同）。
* ``Labels.sum_meta = {"C","DB","O"}``（所有假设一致）。``Labels.chains_meta`` 只反映**已输出层**所断言的内容：
  level 3 -> 每链 ``{"C","DB","toks":{...},"ether":""|"O"|"P"}``（不含 Me/cPr，与 L3 串一致）；
  level 2 -> 每链 ``{"C","DB","toks":{},"ether":""}``（L2 不断言 token 与醚）；level 1 -> ``[]``。
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .config import Config
from .contract import load_vocab
from .errors import Reject
from .nosn import chain_valence, is_ether_class, order_chains, render_chain
from .types import Chain, Interp, Labels, Struct

__all__ = [
    "build_labels", "l3_with_mecpr", "struct_totals", "nominal_chain_counts", "fmt_l1_o", "fmt_l2_tail",
    "chains_meta_l2", "chains_meta_l3",
]

_MECPR = frozenset({"Me", "cPr"})


# ==========================================================================
# 小工具
# ==========================================================================
def fmt_l1_o(n: int) -> str:
    """L1 的氧尾缀：0 -> ``""``，1 -> ``";O"``，n>=2 -> ``";On"``（永不 ``;O1``）。"""
    if n < 0:
        raise ValueError(f"oxygen count must be >= 0, got {n}")
    return "" if n == 0 else ";O" if n == 1 else f";O{n}"


def fmt_l2_tail(n: int) -> str:
    """L2 的尾括号氧：0 -> ``""``，1 -> ``";(O)"``，n>=2 -> ``";(On)"``。"""
    if n < 0:
        raise ValueError(f"oxygen count must be >= 0, got {n}")
    return "" if n == 0 else ";(O)" if n == 1 else f";(O{n})"


def struct_totals(st: Struct) -> tuple[int, int, int]:
    """一个假设的总量 ``(C, DB, O)``；``chains`` 为 None 时取 ``st.total``（``Struct`` 构造时保证非空）。"""
    if st.chains is None:
        if st.total is None:                                   # Struct.__post_init__ 已保证；防御
            raise Reject("internal_error", "Struct with chains=None has no total")
        c, db, o = st.total
        return int(c), int(db), int(o)
    if not st.chains:
        raise Reject("internal_error", "Struct.chains is an empty tuple")
    return (sum(c.C for c in st.chains), sum(c.DB for c in st.chains),
            sum(chain_valence(c) for c in st.chains) + st.tail_o)


def nominal_chain_counts(it: Interp, class_full: str) -> frozenset[int]:
    """类的名义链数集合。假设的类就是选定的类时取 ``it.cls.nchain``；否则查词表 ``class_to_numchain``。
    空集 = 不检查链数（``ST``、``WD`` 等豁免类，或词表里没有）。"""
    if class_full == it.cls.class_full:
        return frozenset(it.cls.nchain)
    n = load_vocab()["class_to_numchain"].get(class_full)
    return frozenset({int(n)}) if n is not None else frozenset()


def chains_meta_l3(chains: tuple[Chain, ...] | list[Chain]) -> list[dict[str, Any]]:
    """level 3 的 ``meta.chains``：与 L3 串一致（不含 Me/cPr）。"""
    out: list[dict[str, Any]] = []
    for c in chains:
        toks: dict[str, int] = {}
        for name, n in c.toks:
            if name not in _MECPR:
                toks[name] = toks.get(name, 0) + n
        out.append({"C": c.C, "DB": c.DB, "toks": toks, "ether": c.ether})
    return out


def chains_meta_l2(chains_meta: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """把某层 ``chains_meta`` 压成 level 2 的形态（只保留 C、DB；L2 不断言 token 与醚）。"""
    return [{"C": m["C"], "DB": m["DB"], "toks": {}, "ether": ""} for m in chains_meta]


# ==========================================================================
# 单个假设
# ==========================================================================
@dataclass(frozen=True, slots=True)
class _Hyp:
    L1: str
    L2: str
    L3: str
    totals: tuple[int, int, int]
    ordered: tuple[Chain, ...]          # 排好序的链（chains is None 时为空）


def _build_one(st: Struct, nominal: frozenset[int], mecpr: bool) -> _Hyp:
    cf = st.class_full
    C, DB, n_o = struct_totals(st)
    L1 = f"{cf} {C}:{DB}{fmt_l1_o(n_o)}"
    if st.chains is None:
        return _Hyp(L1, "", "", (C, DB, n_o), ())
    if nominal and len(st.chains) not in nominal:           # 链清单不全/过多 -> 只到 L1
        return _Hyp(L1, "", "", (C, DB, n_o), ())
    ordered = tuple(order_chains(st.chains, mecpr))
    L2 = f"{cf} " + "_".join(f"{c.C}:{c.DB}" for c in ordered) + fmt_l2_tail(n_o)
    n_eth = sum(1 for c in ordered if c.ether)
    gate = st.tail_o == 0 and (n_eth == 1 if is_ether_class(cf) else n_eth == 0)
    L3 = f"{cf} " + "_".join(render_chain(c, mecpr) for c in ordered) if gate else ""
    return _Hyp(L1, L2, L3, (C, DB, n_o), ordered)


# ==========================================================================
# 合并
# ==========================================================================
def _assemble(it: Interp, mecpr: bool) -> Labels:
    hyps = [_build_one(st, nominal_chain_counts(it, st.class_full), mecpr) for st in it.alts]
    h0 = hyps[0]
    l1s = {h.L1 for h in hyps}
    if len(l1s) > 1:
        raise Reject("internal_error", f"hypotheses disagree on L1: {sorted(l1s)}")
    L2 = h0.L2 if h0.L2 and all(h.L2 == h0.L2 for h in hyps) else ""
    L3 = h0.L3 if L2 and h0.L3 and all(h.L3 == h0.L3 for h in hyps) else ""
    C, DB, n_o = h0.totals
    if L3:
        meta = chains_meta_l3(h0.ordered)
    elif L2:
        meta = chains_meta_l2(chains_meta_l3(h0.ordered))
    else:
        meta = []
    return Labels(L1=h0.L1, L2=L2, L3=L3, level=3 if L3 else 2 if L2 else 1,
                  chains_meta=meta, sum_meta={"C": C, "DB": DB, "O": n_o})


def build_labels(it: Interp, cfg: Config) -> Labels:
    """``Interp`` -> ``Labels``（见模块文档）。假设的 L1 不一致 -> ``Reject("internal_error")``。"""
    return _assemble(it, mecpr=False)


def l3_with_mecpr(it: Interp, cfg: Config) -> str:
    """``meta.L3_with_MeCpr``：``cfg.emit_mecpr_field`` 开着、默认 L3 非空、且至少一条链带 Me/cPr 时，
    返回保留 Me/cPr 的 L3（排序键也含 Me/cPr，与参考表 ``L3_full_nosn_MeCpr`` 写法一致）；否则 ``""``。"""
    if not cfg.emit_mecpr_field:
        return ""
    plain = _assemble(it, mecpr=False).L3
    if not plain:
        return ""
    full = _assemble(it, mecpr=True).L3
    return full if full != plain else ""
