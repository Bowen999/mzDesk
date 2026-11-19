"""M9 对外入口：``normalize`` / ``normalize_many``（SYNTHESIS §2.2 伪代码 + AUDIT_A 的 D13 / D14 / D7）。

对外接口
--------
``normalize(name, *, source=None, config=None) -> dict``
    返回字典的契约见 ``lipid_name_map.contract``（key 恒全、JSON 可序列化、``check_result_contract`` 恒为空）。**永不抛异常**：
    ``Reject`` -> ``rejected`` / ``unsupported`` 字典；其它异常 -> ``rejected/internal_error``（永远是 bug，检查器按失败处理）。
    ``source`` 只有 ``"pnnl"`` 影响解析（D4：多链类单 C:DB 当溶血）；``SOURCES`` 里的其它值只写进 ``meta.source_used``；
    不认识的值当 None。``config`` 不是 ``Config`` 时当默认配置。
``normalize_many(names, *, source=None, config=None) -> list[dict]``    逐个 ``normalize``（单个 ``str`` 当一个名字）
``_merge_segments(outs) -> _Out``      D7 段合并（见下）

每段的流水线（``_run_segment``）::

    check_shorthand(M3) -> split_head(M5a) -> resolve_head(M4a) -> lex_body(M5b) -> select_class(M4b)
      -> interpret(M6, 传入 hi，D13) -> build_labels(M7) -> check_and_downgrade_flags(M8) -> drop_idle_flags(M8) -> lm_names(M10, 0.3.0)

``lm_names``（``lm.py``）只读 ``Interp`` / ``Labels``，给出 LIPID MAPS 风格的 ``LM_species`` / ``LM_molecular``；它不影响 L1-L3 / level，
flags 只加两个：写出了名字且类名是推断的 -> ``lm_head_inferred``；链级名的鞘脂碱基是按"显式羟基更多的那条"推断的 -> ``sphingo_base_inferred``。
自身出错也只会让两个字段留空（``meta.lm`` 里原因码 ``lm_internal_error``）。

固定结构类（``fixed.fixed_labels``：``CoQn`` 与原生 ``CoQ C:DB;O4``）在 M3 之前直接给出三层标签，不走上面的流水线；
flags 只有装饰类，``class_nucleus / class_mod`` 取 ``contract.EXTRA_CLASSES``。

D14 清单：降级 flags 取 ``selfcheck.check_and_downgrade_flags``；装饰类 flags 取 ``decorations.segment_flags(seg)``；
``meta.adduct = seg.adduct``，``meta.stripped = list(seg.notes)``；``Config.emit_mecpr_field`` 时 ``meta.L3_with_MeCpr =
labels.l3_with_mecpr(it, cfg)``；输出的 ``class_nucleus / class_mod / class_full`` 取 ``it.alts[0].class_full`` 再查 ``classes.class_info``
（家族钩子可能改了类，不能直接用 ``select_class`` 的结果）。

D7 段合并（``A | B`` 别名）：每段独立走完整条流水线；各段 ``L1``（含类）必须相同、且任意两段在**双方都非空**的层（L2、L3）上必须
相同，否则 ``alias_conflict``；通过后取层级最深的一段（并列取靠前者）。任何一段被拒绝 -> 整体被拒绝（取第一个拒绝的原因）。
非 ok 结果的 ``meta.class_hint`` 放已识别的类头文本（诊断用）。
"""
from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any

from .classes import class_info, resolve_head, select_class
from .commonname import check_shorthand
from .config import DEFAULT, SOURCES, Config
from .contract import EXTRA_CLASSES, internal_error_dict, ok_dict, rejected_dict
from .decorations import segment_flags, split_decorations
from .fixed import fixed_labels
from .labels import build_labels, l3_with_mecpr
from .lexer import lex_body, split_head
from .lm import LmResult, lm_names
from .preclean import preclean
from .selfcheck import check_and_downgrade_flags, drop_idle_flags
from .semantics import interpret
from .types import Interp, Labels, Reject, Segment

__all__ = ["normalize", "normalize_many"]

_POS_METHYL = re.compile(r"\d+me")          # 位置甲基支链 ``(15Me)``（RawChain.branch）；iso / anteiso 的 C 已含支链碳


@dataclass(frozen=True, slots=True)
class _Out:
    """一段的完整结果（固定结构类没有 ``Interp``：``it`` 为 None）。"""
    seg: Segment
    cf: str
    it: Interp | None
    lb: Labels
    flags: tuple[str, ...]
    mecpr: str = ""
    lm: LmResult = LmResult()


def _run_segment(seg: Segment, cfg: Config, source: str | None, trace: list[str]) -> _Out:
    fx = fixed_labels(seg.text)
    if fx is not None:                                                 # 固定结构类（CoQn）：名字就是整个分子
        cf, lb = fx
        trace.append(cf)
        return _Out(seg, cf, None, lb, tuple(sorted(segment_flags(seg))), lm=lm_names(cf, lb, None))
    check_shorthand(seg, cfg)
    head, body = split_head(seg.text)
    trace.append(head)
    hi = resolve_head(head)
    lx = lex_body(head, body, cfg)
    cr = select_class(hi, lx, cfg)
    it = interpret(lx, cr, cfg, source, hi)
    lb, dflags = check_and_downgrade_flags(build_labels(it, cfg), it)
    flags = set(it.flags)
    flags.update(dflags)
    flags.update(segment_flags(seg))
    flags = set(drop_idle_flags(flags, lb))                            # 没有 L3 就没有"L3 里的 OH / 氧位置是假设"可言
    mecpr = l3_with_mecpr(it, cfg) if (cfg.emit_mecpr_field and lb.L3) else ""
    cf = it.alts[0].class_full
    me = any(_POS_METHYL.fullmatch(rc.branch) for rc in lx.chains)
    lm = lm_names(cf, lb, it, hi.head, "_" not in lx.seps, me)           # LIPID MAPS 名：只读 Interp / Labels，不影响 L1-L3
    return _Out(seg, cf, it, lb, tuple(sorted(flags)), mecpr, lm)


def _merge_segments(outs: list[_Out]) -> _Out:
    """D7 段合并，见模块文档。"""
    if len(outs) == 1:
        return outs[0]
    first = outs[0]
    for o in outs[1:]:
        if o.lb.L1 != first.lb.L1:
            raise Reject("alias_conflict", f"alias segments disagree on L1: {first.lb.L1!r} vs {o.lb.L1!r}")
    for i, a in enumerate(outs):
        for b in outs[i + 1:]:
            if a.lb.L2 and b.lb.L2 and a.lb.L2 != b.lb.L2:
                raise Reject("alias_conflict", f"alias segments disagree on L2: {a.lb.L2!r} vs {b.lb.L2!r}")
            if a.lb.L3 and b.lb.L3 and a.lb.L3 != b.lb.L3:
                raise Reject("alias_conflict", f"alias segments disagree on L3: {a.lb.L3!r} vs {b.lb.L3!r}")
    best = first
    for o in outs[1:]:
        if o.lb.level > best.lb.level:
            best = o
    return best


def _ok_result(name: object, best: _Out, outs: list[_Out], cfg: Config, source: str | None) -> dict:
    if best.cf in EXTRA_CLASSES:
        nucleus, mod = EXTRA_CLASSES[best.cf]
    else:
        ci = class_info(best.cf)
        nucleus, mod = ci.nucleus, ci.mod
    lb = best.lb
    flags = set(best.flags)
    meta: dict[str, Any] = {"sum": dict(lb.sum_meta)}
    if lb.chains_meta:
        meta["chains"] = [dict(m) for m in lb.chains_meta]
    if best.it is not None and best.it.ether_linkage:
        meta["ether_linkage"] = best.it.ether_linkage
    adduct = best.seg.adduct or next((o.seg.adduct for o in outs if o.seg.adduct), "")
    if adduct:
        meta["adduct"] = adduct
    if best.seg.notes:
        meta["stripped"] = list(best.seg.notes)
    if best.mecpr:
        meta["L3_with_MeCpr"] = best.mecpr
    if source is not None:
        meta["source_used"] = source
    lm = best.lm
    meta["lm"] = lm.to_meta()
    if lm.inferred and (lm.species or lm.molecular):
        flags.add("lm_head_inferred")
    if lm.base_inferred and lm.molecular:
        flags.add("sphingo_base_inferred")                               # 鞘脂碱基按"显式羟基更多的那条"推断（lm.py 条件 d③）；只在链级名写出时打
    if len(outs) > 1:
        meta["segments"] = [{"text": o.seg.text, "kind": o.seg.kind, "class_full": o.cf, "L1": o.lb.L1,
                             "level": o.lb.level} for o in outs]
        for o in outs:
            flags.update(segment_flags(o.seg))
    return ok_dict(name, class_nucleus=nucleus, class_mod=mod, class_full=best.cf, L1=lb.L1, L2=lb.L2, L3=lb.L3,
                   LM_species=lm.species, LM_molecular=lm.molecular, flags=flags, meta=meta)


def normalize(name: object, *, source: str | None = None, config: Config | None = None) -> dict:
    """标准化一个脂质 shorthand 名，见模块文档。"""
    cfg = config if isinstance(config, Config) else DEFAULT
    src = source if isinstance(source, str) and source in SOURCES else None
    trace: list[str] = []
    try:
        segs = split_decorations(preclean(name, cfg), cfg)
        outs = [_run_segment(sg, cfg, src, trace) for sg in segs]
        best = _merge_segments(outs)
        return _ok_result(name, best, outs, cfg, src)
    except Reject as r:
        meta: dict[str, Any] = {}
        if trace and trace[-1]:
            meta["class_hint"] = trace[-1][:80]
        if src is not None:
            meta["source_used"] = src
        try:
            return rejected_dict(name, r, meta=meta)
        except Exception as e:                                           # noqa: BLE001 - 永不抛
            return internal_error_dict(name, e)
    except Exception as e:                                               # noqa: BLE001 - 永不抛
        return internal_error_dict(name, e)


def normalize_many(names: Iterable[object], *, source: str | None = None, config: Config | None = None) -> list[dict]:
    """逐个 ``normalize``。"""
    if isinstance(names, (str, bytes)):
        names = [names]
    return [normalize(n, source=source, config=config) for n in names]
