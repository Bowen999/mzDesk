"""输出契约（SYNTHESIS §1，W0 冻结）：返回字典的构造器、契约检查器、flags 词表注册。

对外接口
--------
``ok_dict(name, /, *, class_nucleus, class_mod, class_full, L1, L2="", L3="", LM_species="", LM_molecular="", flags=(), meta=None) -> dict``
``rejected_dict(name, reason, detail="", /, *, flags=(), meta=None) -> dict``
    ``reason`` 是 ``Reject`` 实例（或其 ``code`` 字符串）；status 由 ``reasons`` 查表
``internal_error_dict(name, exc=None) -> dict``      ``rejected/internal_error``（永远是 bug）
    三个构造器返回的字典 **key 恒全**（``RESULT_KEYS``）、**JSON 可序列化**、``flags`` 排序去重、
    ``level`` 按 ``3 if L3 else 2 if L2 else 1``（非 ok 为 0）。``LM_species`` / ``LM_molecular`` 是 LIPID MAPS 名（0.3.0，``lm.py``）：
    非 ok 恒为 ``""``；``LM_molecular`` 非空 => ``LM_species`` 非空。违反契约的参数（空 L1、L3 无 L2、LM 链级名无总量级名、未知 flag、
    未知 meta 键、非 JSON 值…）抛 ``ValueError``/``TypeError``——``api.normalize`` 的兜底会把它变成
    ``internal_error``，所以这些是"当场暴露 bug"，不是用户可见行为。
``check_result_contract(d, original=<未给>) -> list[str]``
    检查 §1.3 的第 1、2、5 条与全部字段类型；返回违规描述（空 = 合规）。给了 ``original`` 才查第 5 条
    与 ``input`` 回显。（第 3、4 条——确定性/线程安全/幂等——需要调用 ``normalize``，见 ``tests/test_contract.py``。）
``level_of(L1, L2, L3) -> int``
``FLAGS``  ``all_flags()``  ``register_family_flags(fam, flags)``  ``registered_family_flags()``
``RESULT_KEYS`` ``META_KEYS`` ``STATUSES``
``load_vocab()`` ``class_full_set()`` ``VOCAB_PATH``     ``data/class_vocab_v0211.json`` 的只读访问
``EXTRA_CLASSES`` ``output_class_set()``                v0.2.11 之外允许输出的类（固定结构类 ``CoQ``，``fixed.py``）；ok 结果的类 ∈ 两者之并

flags 词表：基础词表 ``FLAGS``（SYNTHESIS §1.5，定义在 ``types.FLAGS``，此处同名再导出）；家族在自己模块里声明
``FLAGS_<fam> = frozenset({...})`` 并在导入时调用 ``register_family_flags("<fam>", FLAGS_<fam>)``，之后
``check_result_contract`` 才认这些 flag。名字冲突（与基础词表或别的家族）抛 ``ValueError``。
"""
from __future__ import annotations

import json
import math
import re
import threading
from collections.abc import Iterable, Mapping
from functools import lru_cache
from pathlib import Path
from types import MappingProxyType
from typing import Any, Final

from .reasons import REASONS, STATUS_OK, STATUS_REJECTED, STATUSES, status_of
from .types import ETHER_KINDS, ETHER_LINKAGES, FLAG_GROUPS, FLAGS

__all__ = [
    "RESULT_KEYS", "META_KEYS", "META_LM_KEYS", "LM_HEAD_SOURCES", "STATUSES", "FLAGS", "FLAG_GROUPS", "VOCAB_PATH",
    "ok_dict", "rejected_dict", "internal_error_dict", "check_result_contract", "level_of",
    "all_flags", "register_family_flags", "registered_family_flags",
    "load_vocab", "class_full_set", "EXTRA_CLASSES", "output_class_set",
]

#: §1.2 返回字段，固定顺序，key 恒全
RESULT_KEYS: Final = ("input", "status", "reason", "level", "class_nucleus", "class_mod", "class_full",
                      "L1", "L2", "L3", "LM_species", "LM_molecular", "flags", "meta")
#: §1.2 ``meta`` 允许的顶层键（全部可选）
META_KEYS: Final = frozenset({
    "dialect", "source_used", "segments", "sum", "chains", "ether_linkage", "adduct", "stripped",
    "class_hint", "detail", "L3_with_MeCpr", "lm",
})
#: ``meta.lm`` 允许的键与 ``head_source`` 取值（``lm.py``）
META_LM_KEYS: Final = frozenset({"head_source", "species", "molecular"})
LM_HEAD_SOURCES: Final = ("lmsd", "inferred", "none")

VOCAB_PATH: Final = Path(__file__).resolve().parent / "data" / "class_vocab_v0211.json"
#: v0.2.11 之外允许输出的类 ``{class_full: (nucleus, mod)}``。词表文件是数据集脂质词表的逐字拷贝（sha 锁定），不往里加；
#: 这里只放固定结构类（``fixed.py``）：``CoQ`` 是数据集 v0.2.12 PR_PK 词表的类（2026-10-20 起支持）
EXTRA_CLASSES: Final = MappingProxyType({"CoQ": ("CoQ", "")})


# ==========================================================================
# 词表
# ==========================================================================
@lru_cache(maxsize=1)
def load_vocab() -> dict[str, Any]:
    """读 ``data/class_vocab_v0211.json``（进程内缓存；调用方不得修改返回值）。"""
    with VOCAB_PATH.open(encoding="utf-8") as fh:
        return json.load(fh)


@lru_cache(maxsize=1)
def class_full_set() -> frozenset[str]:
    """v0.2.11 词表的 305 个 ``class_full``。"""
    return frozenset(load_vocab()["class_full"])


@lru_cache(maxsize=1)
def output_class_set() -> frozenset[str]:
    """ok 结果允许的 ``class_full``：v0.2.11 词表的 305 个 + ``EXTRA_CLASSES``。"""
    return class_full_set() | frozenset(EXTRA_CLASSES)


# ==========================================================================
# flags 注册
# ==========================================================================
_FAMILY_FLAGS: dict[str, frozenset[str]] = {}
_FLAGS_LOCK = threading.Lock()
_NAME_RE = re.compile(r"[A-Za-z][A-Za-z0-9_]*")
_FAM_RE = re.compile(r"[a-z][a-z0-9_]*")


def register_family_flags(fam: str, flags: Iterable[str]) -> frozenset[str]:
    """登记模块自带的 flags（``FLAGS_<fam>``）。幂等（同一 ``fam`` 重复登记同一集合无副作用）。

    校验：``fam`` 小写标识符；每个 flag 是标识符；不得与基础词表 ``FLAGS`` 重名；不得与别的 ``fam`` 已登记的重名；
    同一 ``fam`` 再次登记必须是同一集合。违反抛 ``ValueError``。返回登记的 frozenset。
    """
    if not isinstance(fam, str) or not _FAM_RE.fullmatch(fam):
        raise ValueError(f"family name must match {_FAM_RE.pattern}, got {fam!r}")
    fs = frozenset(flags)
    bad = sorted(f for f in fs if not isinstance(f, str) or not _NAME_RE.fullmatch(f))
    if bad:
        raise ValueError(f"illegal flag names {bad}")
    clash = sorted(fs & FLAGS)
    if clash:
        raise ValueError(f"family {fam!r}: flags already in the global vocabulary: {clash}")
    with _FLAGS_LOCK:
        for other, ofs in _FAMILY_FLAGS.items():
            if other != fam and (ofs & fs):
                raise ValueError(f"family {fam!r}: flags {sorted(ofs & fs)} already registered by {other!r}")
        prev = _FAMILY_FLAGS.get(fam)
        if prev is not None and prev != fs:
            raise ValueError(f"family {fam!r} already registered with a different flag set")
        _FAMILY_FLAGS[fam] = fs
    return fs


def registered_family_flags() -> dict[str, frozenset[str]]:
    """已登记的家族 flags 快照 ``{fam: flags}``。"""
    with _FLAGS_LOCK:
        return dict(_FAMILY_FLAGS)


def all_flags() -> frozenset[str]:
    """基础词表 ∪ 已登记的家族 flags。"""
    with _FLAGS_LOCK:
        extra = frozenset().union(*_FAMILY_FLAGS.values()) if _FAMILY_FLAGS else frozenset()
    return FLAGS | extra


def _clear_family_flags_for_tests() -> None:
    with _FLAGS_LOCK:
        _FAMILY_FLAGS.clear()


# ==========================================================================
# 构造器
# ==========================================================================
def level_of(L1: str, L2: str = "", L3: str = "") -> int:
    """``3 if L3 else 2 if L2 else 1 if L1 else 0``（只看哪些层非空；层间包含关系由调用方保证）。"""
    return 3 if L3 else 2 if L2 else 1 if L1 else 0


def _input_str(name: object) -> str:
    if isinstance(name, str):
        return name
    try:
        return repr(name)
    except Exception:  # 自定义 __repr__ 抛错
        return f"<unrepresentable {type(name).__name__}>"


def _jsonable(x: Any, path: str = "meta") -> Any:
    """转成 JSON 原生类型：tuple/list -> list，set/frozenset -> 排序 list；其余非 JSON 值抛 ``TypeError``。"""
    if x is None or isinstance(x, (str, bool, int)):
        return x
    if isinstance(x, float):
        if not math.isfinite(x):
            raise TypeError(f"{path}: non-finite float")
        return x
    if isinstance(x, (list, tuple)):
        return [_jsonable(v, f"{path}[{i}]") for i, v in enumerate(x)]
    if isinstance(x, (set, frozenset)):
        return [_jsonable(v, f"{path}{{}}") for v in sorted(x, key=repr)]
    if isinstance(x, Mapping):
        out = {}
        for k, v in x.items():
            if not isinstance(k, str):
                raise TypeError(f"{path}: non-str key {k!r}")
            out[k] = _jsonable(v, f"{path}.{k}")
        return out
    raise TypeError(f"{path}: {type(x).__name__} is not JSON-serializable")


def _clean_flags(flags: Iterable[str]) -> list[str]:
    if isinstance(flags, str):
        raise TypeError("flags must be an iterable of str, not a bare str")
    fs = set(flags)
    unknown = sorted(f for f in fs if f not in all_flags())
    if unknown:
        raise ValueError(f"unknown flags {unknown} (add to the vocabulary or register_family_flags)")
    return sorted(fs)


def _clean_meta(meta: Mapping[str, Any] | None) -> dict[str, Any]:
    if not meta:
        return {}
    m = _jsonable(dict(meta))
    unknown = sorted(set(m) - META_KEYS)
    if unknown:
        raise ValueError(f"unknown meta keys {unknown}; allowed {sorted(META_KEYS)}")
    return m


def ok_dict(name: object, /, *, class_nucleus: str, class_mod: str, class_full: str, L1: str,
            L2: str = "", L3: str = "", LM_species: str = "", LM_molecular: str = "", flags: Iterable[str] = (),
            meta: Mapping[str, Any] | None = None) -> dict:
    """``status="ok"`` 的返回字典。``level`` 由 L1/L2/L3 推出（LM 名不影响 ``level``）。"""
    if not L1:
        raise ValueError("ok result needs a non-empty L1")
    if L3 and not L2:
        raise ValueError("L3 without L2 is not allowed (no L2 stand-in for L3)")
    if LM_molecular and not LM_species:
        raise ValueError("LM_molecular without LM_species is not allowed")
    if not class_nucleus or not class_full:
        raise ValueError("ok result needs class_nucleus and class_full")
    if class_full != (f"{class_nucleus}-{class_mod}" if class_mod else class_nucleus):
        raise ValueError(f"class_full {class_full!r} != nucleus[-mod] ({class_nucleus!r}, {class_mod!r})")
    if class_full not in output_class_set():
        raise ValueError(f"class_full {class_full!r} is not in the v0.2.11 vocabulary (nor in EXTRA_CLASSES)")
    return {
        "input": _input_str(name), "status": STATUS_OK, "reason": "", "level": level_of(L1, L2, L3),
        "class_nucleus": class_nucleus, "class_mod": class_mod, "class_full": class_full,
        "L1": L1, "L2": L2, "L3": L3, "LM_species": LM_species, "LM_molecular": LM_molecular,
        "flags": _clean_flags(flags), "meta": _clean_meta(meta),
    }


def rejected_dict(name: object, reason: object, detail: str = "", /, *, flags: Iterable[str] = (),
                  meta: Mapping[str, Any] | None = None) -> dict:
    """``rejected`` / ``unsupported`` 的返回字典。

    ``reason`` 是 ``Reject`` 实例（取其 ``code``/``detail``）或原因码字符串；``status`` 由 ``reasons`` 查表。
    非空 ``detail`` 写进 ``meta["detail"]``（调用方已在 ``meta`` 里给了 ``detail`` 则保留）。
    非 ok 结果的 ``class_*``、``L1..L3``、``LM_*`` 恒为空，``level`` 恒为 0（D8）；诊断信息放 ``meta.class_hint``。
    """
    code = getattr(reason, "code", reason)
    if not isinstance(code, str) or code not in REASONS:
        raise ValueError(f"unknown reason code {code!r}")
    if not detail:
        detail = getattr(reason, "detail", "") or ""
    m = _clean_meta(meta)
    if detail:
        m.setdefault("detail", str(detail))
    return {
        "input": _input_str(name), "status": status_of(code), "reason": code, "level": 0,
        "class_nucleus": "", "class_mod": "", "class_full": "", "L1": "", "L2": "", "L3": "",
        "LM_species": "", "LM_molecular": "", "flags": _clean_flags(flags), "meta": m,
    }


def internal_error_dict(name: object, exc: BaseException | None = None) -> dict:
    """``rejected/internal_error``。只记异常类型与截断后的消息，不带 traceback。"""
    detail = ""
    if exc is not None:
        detail = f"{type(exc).__name__}: {exc}"[:200]
    return rejected_dict(name, "internal_error", detail)


# ==========================================================================
# 检查器
# ==========================================================================
_UNSET: Final = object()


def _is_int(x: object) -> bool:
    return isinstance(x, int) and not isinstance(x, bool)


def _check_meta(meta: dict, errs: list[str]) -> None:
    unknown = sorted(set(meta) - META_KEYS)
    if unknown:
        errs.append(f"meta has unknown keys {unknown}")
    for k in ("dialect", "source_used", "adduct", "class_hint", "detail", "L3_with_MeCpr"):
        if k in meta and not isinstance(meta[k], str):
            errs.append(f"meta.{k} must be str, got {type(meta[k]).__name__}")
    if "stripped" in meta and not (isinstance(meta["stripped"], list) and all(isinstance(s, str) for s in meta["stripped"])):
        errs.append("meta.stripped must be list[str]")
    if "segments" in meta and not (isinstance(meta["segments"], list) and all(isinstance(s, dict) for s in meta["segments"])):
        errs.append("meta.segments must be list[dict]")
    if "ether_linkage" in meta and meta["ether_linkage"] not in ETHER_LINKAGES:
        errs.append(f"meta.ether_linkage must be one of {ETHER_LINKAGES}, got {meta['ether_linkage']!r}")
    if "sum" in meta:
        s = meta["sum"]
        if not (isinstance(s, dict) and set(s) == {"C", "DB", "O"} and all(_is_int(v) for v in s.values())):
            errs.append("meta.sum must be {'C': int, 'DB': int, 'O': int}")
    if "lm" in meta:
        lm = meta["lm"]
        if not (isinstance(lm, dict) and set(lm) <= META_LM_KEYS and all(isinstance(v, str) for v in lm.values())
                and lm.get("head_source") in LM_HEAD_SOURCES):
            errs.append(f"meta.lm must be {{'head_source': one of {LM_HEAD_SOURCES}, 'species'?: str, 'molecular'?: str}}, got {lm!r}")
    if "chains" in meta:
        ch = meta["chains"]
        if not isinstance(ch, list):
            errs.append("meta.chains must be a list")
        else:
            for i, c in enumerate(ch):
                if not (isinstance(c, dict) and _is_int(c.get("C")) and _is_int(c.get("DB"))
                        and isinstance(c.get("toks", {}), dict) and c.get("ether", "") in ETHER_KINDS):
                    errs.append(f"meta.chains[{i}] must be {{'C': int, 'DB': int, 'toks': dict, 'ether': ''|'O'|'P'}}")
                    break


def check_result_contract(d: object, original: object = _UNSET) -> list[str]:
    """检查一个 ``normalize`` 结果字典是否满足 §1.3 的第 1、2、5 条与字段类型；返回违规描述（空 = 合规）。

    * 第 1 条：``status=="ok"`` ⇔ ``level>=1`` ⇔ ``L1!=""``；非 ok ⇒ ``level==0``、``L1==L2==L3==""``、
      ``class_*==""``、``reason!=""``（且是已注册的码、status 与码一致）
    * 第 2 条：``L3!="" ⇒ L2!="" ⇒ L1!=""``；``level == 3 if L3 else 2 if L2 else 1``
    * LM 名（0.3.0）：非 ok ⇒ ``LM_species == LM_molecular == ""``；``LM_molecular != "" ⇒ LM_species != "" 且 L1 != ""``
    * 第 5 条（给了 ``original`` 才查）：非 str -> ``rejected/not_string``；``""`` -> ``rejected/empty``；
      全空白串 -> 非 ok（rejected）；``input`` 回显 = 原串（非 str 则 ``repr``）
    * 字段类型：key 恒全且无多余、flags 为排序去重且在词表内、meta 顶层键在 ``META_KEYS`` 内、整个字典可 JSON 序列化
    * ok 时：``class_full == nucleus[-mod]`` 且 ∈ v0.2.11 词表（305）或 ``EXTRA_CLASSES``；非空的各层都以 ``class_full + " "`` 开头
    """
    if not isinstance(d, dict):
        return [f"result is {type(d).__name__}, not dict"]
    errs: list[str] = []
    missing = [k for k in RESULT_KEYS if k not in d]
    extra = [k for k in d if k not in RESULT_KEYS]
    if missing:
        errs.append(f"missing keys {missing}")
    if extra:
        errs.append(f"unexpected keys {extra}")
    if missing:
        return errs

    # ---- 类型 ----
    for k in ("input", "reason", "class_nucleus", "class_mod", "class_full", "L1", "L2", "L3", "LM_species", "LM_molecular"):
        if not isinstance(d[k], str):
            errs.append(f"{k} must be str, got {type(d[k]).__name__}")
    if d["status"] not in STATUSES:
        errs.append(f"status must be one of {STATUSES}, got {d['status']!r}")
    if not _is_int(d["level"]) or not 0 <= d["level"] <= 3:
        errs.append(f"level must be int in 0..3, got {d['level']!r}")
    flags = d["flags"]
    if not (isinstance(flags, list) and all(isinstance(f, str) for f in flags)):
        errs.append("flags must be list[str]")
    if not isinstance(d["meta"], dict):
        errs.append("meta must be dict")
    if errs:
        return errs

    if flags != sorted(set(flags)):
        errs.append("flags must be sorted and de-duplicated")
    known = all_flags()
    unknown = sorted(set(flags) - known)
    if unknown:
        errs.append(f"unknown flags {unknown}")
    _check_meta(d["meta"], errs)
    try:
        json.dumps(d, allow_nan=False)
    except (TypeError, ValueError) as e:
        errs.append(f"not JSON-serializable: {e}")

    status, reason, level = d["status"], d["reason"], d["level"]
    L1, L2, L3 = d["L1"], d["L2"], d["L3"]

    # ---- §1.3-2 ----
    if L3 and not L2:
        errs.append("L3 non-empty but L2 empty")
    if L2 and not L1:
        errs.append("L2 non-empty but L1 empty")

    # ---- LM 名（0.3.0）：非 ok => 两个字段空；LM_molecular 非空 => LM_species 非空且 L1 非空 ----
    LMs, LMm = d["LM_species"], d["LM_molecular"]
    if status != STATUS_OK and (LMs or LMm):
        errs.append(f"status {status} but LM_species/LM_molecular not both empty")
    if LMm and not LMs:
        errs.append("LM_molecular non-empty but LM_species empty")
    if (LMs or LMm) and not L1:
        errs.append("LM_species/LM_molecular non-empty but L1 empty")

    # ---- §1.3-1 ----
    if status == STATUS_OK:
        if level < 1:
            errs.append("status ok but level < 1")
        if not L1:
            errs.append("status ok but L1 empty")
        if level != (3 if L3 else 2 if L2 else 1):
            errs.append(f"level {level} inconsistent with filled layers (L2={bool(L2)}, L3={bool(L3)})")
        if reason != "":
            errs.append(f"status ok but reason {reason!r}")
        nuc, mod, cf = d["class_nucleus"], d["class_mod"], d["class_full"]
        if not nuc or not cf:
            errs.append("status ok but class_nucleus/class_full empty")
        else:
            if cf != (f"{nuc}-{mod}" if mod else nuc):
                errs.append(f"class_full {cf!r} != nucleus[-mod] ({nuc!r}, {mod!r})")
            if cf not in output_class_set():
                errs.append(f"class_full {cf!r} not in the v0.2.11 vocabulary (nor in EXTRA_CLASSES)")
            for lab, s in (("L1", L1), ("L2", L2), ("L3", L3)):
                if s and not s.startswith(cf + " "):
                    errs.append(f"{lab} {s!r} does not start with class_full {cf!r}")
    else:
        if level != 0:
            errs.append(f"status {status} but level {level} != 0")
        if L1 or L2 or L3:
            errs.append(f"status {status} but L1/L2/L3 not all empty")
        if d["class_nucleus"] or d["class_mod"] or d["class_full"]:
            errs.append(f"status {status} but class_* not empty")
        if not reason:
            errs.append(f"status {status} but reason empty")
        elif reason not in REASONS:
            errs.append(f"unknown reason code {reason!r}")
        elif REASONS[reason].status != status:
            errs.append(f"reason {reason!r} belongs to status {REASONS[reason].status!r}, not {status!r}")

    # ---- §1.3-5 + input 回显 ----
    if original is not _UNSET:
        if not isinstance(original, str):
            if not (status == STATUS_REJECTED and reason == "not_string"):
                errs.append(f"non-str input must give rejected/not_string, got {status}/{reason}")
            if d["input"] != _input_str(original):
                errs.append("input must echo repr() of a non-str argument")
        else:
            if d["input"] != original:
                errs.append("input must echo the original string")
            if original == "" and not (status == STATUS_REJECTED and reason == "empty"):
                errs.append(f"empty string must give rejected/empty, got {status}/{reason}")
            elif original.strip() == "" and status != STATUS_REJECTED:
                errs.append(f"whitespace-only string must be rejected, got {status}/{reason}")
    return errs
