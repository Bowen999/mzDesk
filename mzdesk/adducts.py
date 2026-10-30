"""加合物写法的统一入口：规范写法、离子极性、各引擎的写法（ESI MS/MS 注释工作流和各模块共用）。

规范写法 = library_parse.normalize_adduct 的输出（"[M+H]+"、"[M+HCOO]-" ……），谱库索引、加合物下拉框、工作流都按它比较。
甲酸根写法最多：[M+HCOO]- / [M+FA-H]- / [M+HCOOH-H]- / [M+HCO2]- / [M+FA]- / [M+HCOOH]- / M+formate，规范写法都是 [M+HCOO]-，
离子 m/z = M + 44.998203（HCOO 加一个电子）。
引擎写法取自引擎模块自己的加合物表（不另抄一份）：引擎页面能用的，这里就认。引擎模块第一次用到时才导入
（msbuddy_engine 约 1 s，metfrag_engine 会拉进 RDKit）。
"""
from functools import lru_cache

from mzdesk.ms2.library_parse import adduct_charge_sign, normalize_adduct

# CASMI 2022 的六种加合物（benchmark/casmi2022/metadata/casmi2022_challenges.csv，500 条：177 / 108 / 67 / 67 / 60 / 21）
CASMI_ADDUCTS = ("[M+H]+", "[M-H]-", "[M+NH4]+", "[M+HCOO]-", "[M+Na]+", "[M+Cl]-")
ENGINES = ("buddy", "msfinder", "metfrag", "frigid")


def canonical(name):
    """任意写法 -> 规范写法；"Auto"、空、认不出 -> None。"""
    return normalize_adduct(name)


def polarity(name):
    """"Positive" / "Negative"；认不出 -> None。"""
    return {1: "Positive", -1: "Negative"}.get(adduct_charge_sign(canonical(name)))


@lru_cache(maxsize=None)
def _engine_names(engine):
    """规范写法 -> 该引擎接受的写法。"""
    if engine == "buddy":
        from mzdesk.ms1 import msbuddy_engine as m
        names = m.ADDUCTS_POSITIVE + m.ADDUCTS_NEGATIVE
    elif engine == "msfinder":
        from mzdesk.ms2 import msfinder_engine as m
        names = m.ADDUCTS["Positive"] + m.ADDUCTS["Negative"]
    elif engine == "metfrag":
        from mzdesk.ms2 import metfrag_engine as m
        names = list(m.ADDUCTS)
    elif engine == "frigid":
        from mzdesk.ms2 import frigid_engine as m
        names = list(m.ADDUCTS)                     # 只有正离子
    else:
        raise ValueError(f"Unknown engine {engine!r}; choose one of {', '.join(ENGINES)}.")
    return {c: n for c, n in ((canonical(n), n) for n in names) if c}


def engine_name(name, engine):
    """name 在该引擎里的写法，即 msbuddy_engine.annotate_formula / msfinder_engine.predict / metfrag_engine.run_metfrag /
    frigid_engine.predict 的 adduct 参数；引擎不支持、name 认不出或是 "Auto" -> None。engine：buddy | msfinder | metfrag | frigid。

    例：engine_name("[M+HCOO]-", "msfinder") == "[M+FA-H]-"；engine_name("[M+FA-H]-", "buddy") == "[M+HCOO]-"；
    engine_name("[M-H]-", "frigid") is None。
    """
    return _engine_names(engine).get(canonical(name))


def supported_engines(name):
    """支持 name 的引擎，按 ENGINES 的顺序（tuple，可能为空）。"""
    return tuple(e for e in ENGINES if engine_name(name, e))
