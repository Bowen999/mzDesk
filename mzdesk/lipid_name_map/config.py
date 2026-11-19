"""运行配置（SYNTHESIS §1.6，W0 冻结；字段名不可改）。

对外接口
--------
``Config``      frozen dataclass，字段见下；``dataclasses.replace(cfg, x=...)`` 派生新配置
``DEFAULT``     ``Config()``
``SOURCES``     ``normalize(source=...)`` 允许的来源提示集合（目前只有 ``"pnnl"`` 影响解析，D4）
``Config.single_token_is_lyso(source)``  D4：多链类只写一个 C:DB 时是否当溶血

字段
----
allow_oxylipin_abbrev      True   展开 7 个介质缩写 HETE/HEPE/HDoHE/HODE/EET/EpETE/EpDPE（C18）
allow_unverified_alias     False  D2：unverified 别名默认 ``unsupported(class_alias_unverified)``
assume_sphingolipid_OH     True   §4.3 sphingo_OH 策略：鞘脂未分型氧按羟基出 L3
assume_base_oxygen         True   C26（2026-10-20 改判）：鞘脂名字一个氧都没给（``SM 30:1``、``SPB 18:1``）时按 d 型碱基补 2 个羟基
                                  （flag ``base_oxygen_assumed``）；False 退回原设计，拒绝 ``missing_oxygen_count``
class_implied_OH_in_L3     True   §4.3 class_OH 策略：类隐含氧作为角色链 ``;OH`` 出 L3
single_token_as_lyso       "auto" D4：``"auto"``=仅 source=="pnnl" 时真；也可显式 True/False
emit_mecpr_field           False  名字显式写了 Me/cPr 支链时，另给 ``meta.L3_with_MeCpr``
pnnl_ether_heuristic       False  开则退回 benchmark 的"放饱和链"启发式（不推荐）
max_len                    1000   超长输入 -> ``rejected/too_long``
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Final

#: normalize(source=...) 的合法提示值（None 也合法）；目前只有 "pnnl" 影响解析（D4），其余只写进 meta.source_used
SOURCES: Final = frozenset({
    "pnnl", "lipidin", "lipidmatch", "lipiddetective", "elgordo",
    "msdial", "lipidblast", "lipidmaps", "swisslipids", "riken",
})


@dataclass(frozen=True, slots=True)
class Config:
    allow_oxylipin_abbrev: bool = True
    allow_unverified_alias: bool = False
    assume_sphingolipid_OH: bool = True
    assume_base_oxygen: bool = True
    class_implied_OH_in_L3: bool = True
    single_token_as_lyso: str | bool = "auto"
    emit_mecpr_field: bool = False
    pnnl_ether_heuristic: bool = False
    max_len: int = 1000

    def __post_init__(self) -> None:
        stl = self.single_token_as_lyso
        if not (isinstance(stl, bool) or stl == "auto"):
            raise ValueError(f"single_token_as_lyso must be 'auto', True or False, got {self.single_token_as_lyso!r}")
        if isinstance(self.max_len, bool) or not isinstance(self.max_len, int) or self.max_len < 1:
            raise ValueError(f"max_len must be a positive int, got {self.max_len!r}")

    def single_token_is_lyso(self, source: str | None) -> bool:
        """D4：多链类只写一个 C:DB 时，是否当溶血（True）还是总量（False）。"""
        if self.single_token_as_lyso == "auto":
            return source == "pnnl"
        return bool(self.single_token_as_lyso)


DEFAULT: Final = Config()
