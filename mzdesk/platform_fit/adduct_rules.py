"""Adducts 页的规则引擎：结构 + ESI 极性 + 流动相添加剂 → 会出现哪些离子（加合物、常见源内丢失、二聚体），按预期强度排序，每个离子给出理由。

没有可靠的通用预测工具，所以这里是透明的打分规则，不是模型。每个候选离子得一个 0-100 的分数（大致是"相对一个理想离子化分子的预期响应"），
分数 → likelihood：high >= 50，medium >= 20，其余 low。所有分数和阈值都在 Rules 里（一处）。

1. 结构：pka_engine.prepare（最大片段、能中和的电荷先中和）。M 是中性分子；有永久正电荷（季铵等）且带酸性基团时 M 取内盐（两性离子，
   [M+H]+ 就是那个阳离子），没有酸性基团时 M 就是阳离子本身（离子 [M]+）。m/z 和离子式用 msbuddy_engine.ion_mz_and_formula，
   与 Theoretical m/z 页同一套代码，数值一致。
2. 脂质：lipid_metabolite_classifier 给出大类，这里按头基 / 链数细分到 PC、PE、PS、PG、PI、PIP、PA、CL、SM、Cer、HexCer、TG、DG、MG、
   糖基甘油二酯、CE（甾醇酯）、游离甾醇、FA，按各类常用的加合物约定打分（LIPID MAPS、MS-DIAL、LipidSearch 的做法）。分不出亚类的脂质走小分子规则。
3. 小分子：MolGpKa 的位点 pKa（权重缺失或失败时用官能团的典型 pKa，Rule hits 里说明）+ 几组 SMARTS：
   - 碱性位点 pKa >= 3 → [M+H]+ high；没有这样的碱但有能接受质子的基团（酰胺 / 内酰胺 C=O、芳香 N、共轭羰基、P=O、S=O、N-氧化物、
     亚胺、间位二氧代苯环）→ [M+H]+ 仍然 high（液滴 / 气相质子化）；只有 O 的中性分子 → [M+Na]+ / [M+NH4]+ 高于 [M+H]+；
     没有 N / O / P / S → 都是 low（ESI 响应差）。
   - 永久正电荷 → [M]+，没有 [M+H]+；负模式下不出峰。
   - 强酸（pKa < 5：羧酸、磺酸 / 硫酸酯、磷酸、四唑、酰基磺酰胺）→ [M-H]- high；弱酸（酚、酰亚胺、磺酰胺 NH，pKa < 10.5）→ [M-H]- high；
     多元醇 OH（pKa < 13）→ [M-H]- medium；中性极性分子在负模式下 → [M+HCOO]- / [M+CH3COO]-（看添加剂），[M+Cl]- low。
   - 一致性：只有强酸、没有碱性 / 质子受体的分子，正模式所有离子乘 0.3（有质子受体乘 0.6），并且封顶在 medium 以下。
   - 源内丢失：叔 / 苄位 / 烯丙位醇、β-羟基羰基 → [M+H-H2O]+；糖 / 多元醇 → [M+H-H2O]+ low；芳环或羰基旁的伯胺 → [M+H-NH3]+ low；
     芳基乙酸、水杨酸型、β-酮酸 → [M-H-CO2]- low。
   - 二聚体 [2M+H]+ / [2M+Na]+ / [2M+NH4]+ / [2M-H]- low（与浓度有关），每个极性最多一个，只给 MW <= 500 的小分子。
   - 多电荷：>= 2 个 pKa >= 6 的碱且 MW >= 300 → [M+2H]2+；>= 2 个强酸位点且 MW >= 300 → [M-2H]2-。
4. 流动相 pH（水相，添加剂之外什么都不加）：甲酸 0.1 % ≈ 2.7，乙酸 0.1 % ≈ 3.2，甲酸铵 5-10 mM（不调 pH）≈ 6.5，乙酸铵 5-10 mM ≈ 6.8，
   无添加剂 ≈ 7（不缓冲，实际 pH 不确定）。溶液里的电离比例只把分数在 0.8-1 倍之间调整：ESI 里酸在酸性流动相中照样出 [M-H]-，
   碱在中性流动相中照样出 [M+H]+（"wrong-way-round" 电离），所以起决定作用的是有没有可电离的位点。

依据（简要）：Cech & Enke, Mass Spectrom. Rev. 2001（ESI 响应与碱性 / 表面活性）；Kebarle & Verkerk, Mass Spectrom. Rev. 2009（ESI 机理）；
Kruve & Kaupmees, JASMS 2017（流动相添加剂与 Na+ / NH4+ 加合物）；Liigand et al., JASMS 2014（流动相对电离效率的影响）；
Xu et al., Anal. Chem. 2015（源内碎裂造成的误注释）；Hsu & Turk 的磷脂 ESI 系列工作、Murphy《Tandem Mass Spectrometry of Lipids》2015、
MS-DIAL / LipidSearch 的脂质类别加合物表（各类脂质的约定离子）。
"""
import re
from dataclasses import dataclass, field
from functools import lru_cache

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from mzdesk.platform_fit import lipid_metabolite_classifier as clf
from mzdesk.platform_fit import pka_engine as pe

ADDITIVES = ("Formic acid", "Ammonium formate", "Ammonium acetate", "Acetic acid", "None")
POLARITIES = ("Positive", "Negative", "Both")


def _lipid_ions():
    """Class conventions: code -> [(polarity, ion, score with an ammonium additive, score without, reason)].

    "[M+A]-" is the additive anion adduct ([M+HCOO]- with formic acid / ammonium formate, [M+CH3COO]- with acetic acid / ammonium acetate);
    it is dropped with no additive. [M+NH4]+ needs an ammonium additive (score without = 0).
    """
    choline = "choline headgroup carries a fixed + charge"
    no_h = "zwitterion with no free acidic H: seen as the anion adduct of the additive (PC / SM convention), not [M-H]-"
    return {
        "PC": [("+", "[M+H]+", 100, 100, f"{choline}; the phosphate takes H+ (standard PC / SM positive ion)"),
               ("+", "[M+Na]+", 25, 40, "Na+ replaces the phosphate H (common second PC / SM ion)"),
               ("+", "[M+K]+", 8, 12, "K+ adduct, minor"),
               ("-", "[M+A]-", 100, 100, no_h),
               ("-", "[M-CH3]-", 10, 10, "demethylated choline, formed in the source from the anion adduct"),
               ("-", "[M+Cl]-", 10, 15, "chloride traces")],
        "PE": [("+", "[M+H]+", 100, 100, "ethanolamine amine takes H+ (standard PE positive ion)"),
               ("+", "[M+Na]+", 20, 35, "Na+ adduct"),
               ("-", "[M-H]-", 100, 100, "phosphodiester loses H+ (standard PE negative ion)")],
        "PS": [("-", "[M-H]-", 100, 100, "net negative headgroup (phosphate + carboxylate); [M-H]- is the standard PS ion"),
               ("+", "[M+H]+", 45, 45, "serine amine takes H+, weaker than negative mode"),
               ("+", "[M+Na]+", 15, 25, "Na+ adduct")],
        "PG": [("-", "[M-H]-", 100, 100, "acidic phosphodiester; [M-H]- is the standard PG ion"),
               ("+", "[M+NH4]+", 45, 0, "NH4+ adduct (used for PG in positive mode with ammonium additives)"),
               ("+", "[M+Na]+", 15, 30, "Na+ adduct"),
               ("+", "[M+H]+", 15, 20, "weak: acidic lipid without a basic group")],
        "PI": [("-", "[M-H]-", 100, 100, "acidic phosphodiester; [M-H]- is the standard PI ion"),
               ("+", "[M+NH4]+", 45, 0, "NH4+ adduct (used for PI in positive mode with ammonium additives)"),
               ("+", "[M+Na]+", 15, 25, "Na+ adduct"),
               ("+", "[M+H]+", 10, 12, "weak: acidic lipid without a basic group")],
        "PIP": [("-", "[M-H]-", 100, 100, "polyphosphoinositide: several acidic phosphates"),
                ("-", "[M-2H]2-", 60, 60, "two acidic phosphates: doubly deprotonated ion is common"),
                ("+", "[M+NH4]+", 30, 0, "NH4+ adduct, weak"),
                ("+", "[M+H]+", 5, 5, "very weak: strongly acidic lipid")],
        "PA": [("-", "[M-H]-", 100, 100, "phosphomonoester loses H+; [M-H]- is the standard PA ion"),
               ("+", "[M+NH4]+", 40, 0, "NH4+ adduct (PA in positive mode with ammonium additives)"),
               ("+", "[M+Na]+", 15, 25, "Na+ adduct"),
               ("+", "[M+H]+", 10, 15, "weak: acidic lipid without a basic group")],
        "CL": [("-", "[M-H]-", 100, 100, "two acidic phosphodiesters; [M-H]- is the usual CL ion"),
               ("-", "[M-2H]2-", 70, 70, "two phosphates: [M-2H]2- is a characteristic CL ion"),
               ("+", "[M+NH4]+", 30, 0, "NH4+ adduct, weak"),
               ("+", "[M+Na]+", 10, 15, "Na+ adduct, weak")],
        "TG": [("+", "[M+NH4]+", 100, 0, "neutral lipid with six ester O: NH4+ adduct is the standard TG ion with ammonium additives"),
               ("+", "[M+Na]+", 40, 100, "Na+ adduct (the main TG ion without ammonium)"),
               ("+", "[M+K]+", 12, 30, "K+ adduct"),
               ("+", "[M+H]+", 8, 15, "little [M+H]+: no proton-accepting group")],
        "DG": [("+", "[M+NH4]+", 100, 0, "neutral lipid: NH4+ adduct is the standard DG ion with ammonium additives"),
               ("+", "[M+Na]+", 40, 90, "Na+ adduct"),
               ("+", "[M+H-H2O]+", 25, 30, "free OH: DG loses water in the source"),
               ("+", "[M+K]+", 10, 25, "K+ adduct"),
               ("+", "[M+H]+", 8, 15, "little [M+H]+"),
               ("-", "[M+A]-", 15, 15, "weak anion adduct; DG is rarely measured in negative mode")],
        "MG": [("+", "[M+NH4]+", 90, 0, "NH4+ adduct (MG with ammonium additives)"),
               ("+", "[M+Na]+", 50, 90, "Na+ adduct"),
               ("+", "[M+H]+", 30, 40, "[M+H]+ seen for MG"),
               ("+", "[M+H-H2O]+", 25, 30, "two free OH: loses water in the source"),
               ("-", "[M+A]-", 25, 25, "anion adduct, weak")],
        "GlyDG": [("+", "[M+NH4]+", 100, 0, "glycosyl diacylglycerol (MGDG / DGDG): NH4+ adduct with ammonium additives"),
                  ("+", "[M+Na]+", 50, 100, "Na+ adduct (sugar and ester O)"),
                  ("-", "[M+A]-", 80, 80, "sugar OH groups bind the additive anion (MGDG / DGDG negative ion)"),
                  ("-", "[M-H]-", 15, 15, "weak: no acidic group")],
        "Cer": [("+", "[M+H]+", 100, 100, "amide and OH groups: [M+H]+ is the standard ceramide positive ion"),
                ("+", "[M+H-H2O]+", 55, 55, "sphingoid 1,3-diol loses water readily in the source"),
                ("+", "[M+Na]+", 25, 40, "Na+ adduct"),
                ("-", "[M+A]-", 100, 100, "OH / NH groups bind the additive anion (standard ceramide negative ion)"),
                ("-", "[M-H]-", 45, 45, "amide / OH deprotonation, weaker than the anion adduct"),
                ("-", "[M+Cl]-", 10, 15, "chloride traces")],
        "CE": [("+", "[M+NH4]+", 100, 0, "steryl ester: NH4+ adduct is the standard CE ion with ammonium additives"),
               ("+", "[M+Na]+", 40, 100, "Na+ adduct (the main CE ion without ammonium)"),
               ("+", "[M+K]+", 10, 25, "K+ adduct")],
        "ST": [("+", "[M+H-H2O]+", 45, 45, "3β-OH sterol: protonation is followed by water loss; weak in ESI (APCI is the usual choice)"),
               ("+", "[M+NH4]+", 40, 0, "NH4+ adduct (used for sterols with ammonium additives); weak in ESI"),
               ("+", "[M+Na]+", 12, 20, "Na+ adduct, weak"),
               ("+", "[M+H]+", 5, 5, "intact [M+H]+ is barely seen"),
               ("-", "[M+A]-", 8, 8, "very weak: free sterols ionize poorly in negative ESI")],
        "FA": [("-", "[M-H]-", 100, 100, "free carboxylic acid: [M-H]- is the standard fatty acid ion"),
               ("-", "[2M-H]-", 10, 10, "acid dimer at high concentration"),
               ("+", "[M+NH4]+", 10, 0, "weak: fatty acids ionize poorly in positive mode"),
               ("+", "[M+Na]+", 8, 12, "weak Na+ adduct"),
               ("+", "[M+H]+", 3, 3, "barely protonated")],
    }


@dataclass(frozen=True)
class Rules:
    # mobile-phase pH implied by the additive (aqueous phase; 0.1 % acid or 5-10 mM salt, nothing else added)
    ph_formic: float = 2.7
    ph_acetic: float = 3.2
    ph_ammonium_formate: float = 6.5
    ph_ammonium_acetate: float = 6.8
    ph_none: float = 7.0
    # site classes from pKa
    base_pka: float = 3.0             # base at or above this protonates readily -> [M+H]+ high
    multi_base_pka: float = 6.0       # bases counted for [M+2H]2+
    strong_acid_pka: float = 5.0      # carboxylic, sulfonic, phosphate, tetrazole, acylsulfonamide
    weak_acid_pka: float = 10.5       # phenol, imide, sulfonamide NH
    very_weak_acid_pka: float = 13.0  # sugar / polyol OH (MolGpKa gives about 12)
    inner_salt_pka: float = 7.0       # permanent cation + an acid below this -> M is the inner salt
    # likelihood from the score
    high: float = 50.0
    medium: float = 20.0
    min_shown: float = 2.0
    ph_floor: float = 0.8             # solution ionization moves a score only between 0.8x and 1x
    # positive mode
    h_base: float = 100.0
    h_weak: float = 65.0              # proton-accepting group, no base with pKa >= base_pka
    h_polar: float = 12.0             # O-only neutral
    h_none: float = 3.0               # no N, O, P, S
    radical: float = 5.0              # [M]+. of aromatic hydrocarbons
    na_base: float = 15.0
    na_weak: float = 30.0
    na_polar: float = 40.0
    na_polyoxy: float = 75.0
    polyoxy_min_o: int = 3            # this many O atoms -> "polyoxygenated" (sugars, glycosides, esters)
    k_factor: float = 0.3             # [M+K]+ relative to [M+Na]+
    nh4_weak: float = 25.0
    nh4_polar: float = 40.0
    nh4_polyoxy: float = 85.0
    na_with_ammonium: float = 0.6     # NH4+ competes with Na+
    na_with_acid: float = 0.85        # H+ competes a little
    na_with_none: float = 1.3         # no other cation source: Na+ / K+ rise
    na_absent: float = 0.1            # "Na / K present" unticked: trace Na only, no K
    weak_with_ammonium: float = 0.85  # fewer free protons than with an acid
    weak_with_none: float = 0.75
    acid_pos: float = 0.3             # strong acid, no proton acceptor: every positive ion x this
    acid_pos_weak: float = 0.6        # strong acid + proton acceptor
    # negative mode
    mh_strong: float = 100.0
    mh_weak: float = 75.0
    mh_very_weak: float = 30.0
    mh_none: float = 4.0
    anion_donors2: float = 70.0       # >= 2 OH / NH donors bind HCOO- / CH3COO-
    anion_donors1: float = 40.0
    anion_donors0: float = 15.0
    anion_with_acid: float = 5.0      # acids lose H+ instead
    anion_with_weak_acid: float = 0.3
    anion_with_base: float = 0.3
    cl_factor: float = 0.25           # [M+Cl]- relative to the anion adduct score
    cl_factor_none: float = 0.35
    lipid_cl_none: float = 1.5        # lipid [M+Cl]- with no additive anion
    # inner salts (betaine, carnitine)
    quat_na: float = 20.0
    quat_anion: float = 20.0
    dication_mh: float = 30.0         # [M-H]+ of a permanent dication
    # in-source losses, dimers, multiple charges
    water_frac: float = 0.35          # [M+H-H2O]+ = this x [M+H]+ for loss-prone OH
    water_min: float = 10.0
    water_sugar: float = 12.0
    nh3_loss: float = 12.0
    co2_loss: float = 10.0
    dimer: float = 8.0
    dimer_max_mw: float = 500.0
    z2_min_mw: float = 300.0
    z2_mid_mw: float = 700.0
    z2_high_mw: float = 1200.0
    z2_low: float = 10.0
    z2_mid: float = 35.0
    z2_high: float = 60.0
    z2_acid_mid_mw: float = 600.0
    z2_acid_mid: float = 40.0
    max_sites_listed: int = 8
    lipid_ions: dict = field(default_factory=_lipid_ions)


RULES = Rules()


@dataclass(frozen=True)
class Additive:
    name: str
    ph: float
    ammonium: bool
    anion: str            # "HCOO" | "CH3COO" | ""
    acid: bool
    assumption: str


def additive(name, r=RULES):
    table = {
        "Formic acid": (r.ph_formic, False, "HCOO", True, "0.1 % formic acid in water"),
        "Acetic acid": (r.ph_acetic, False, "CH3COO", True, "0.1 % acetic acid in water"),
        "Ammonium formate": (r.ph_ammonium_formate, True, "HCOO", False, "5-10 mM ammonium formate, pH not adjusted"),
        "Ammonium acetate": (r.ph_ammonium_acetate, True, "CH3COO", False, "5-10 mM ammonium acetate, pH not adjusted"),
        "None": (r.ph_none, False, "", False, "water / organic without additive; unbuffered, so the real pH is uncertain"),
    }
    if name not in table:
        raise ValueError(f"Unknown additive {name!r}; choose one of {', '.join(ADDITIVES)}.")
    return Additive(name, *table[name])


# ---------------------------------------------------------------- functional groups (site atom first in every SMARTS)
ACID_GROUPS = (
    ("carboxylic acid", "[OX2H1][CX3]=O"),
    ("phosphate / phosphonate", "[OX2H1][PX4]=O"),
    ("sulfonic acid / sulfate", "[OX2H1][SX4](=O)=O"),
    ("tetrazole NH", "[#7H1]1:[#7]:[#7]:[#7]:[#6]:1"),
    ("tetrazole NH", "[#7H1]1:[#7]:[#7]:[#6]:[#7]:1"),
    ("acylsulfonamide NH", "[NX3H1]([CX3]=O)[SX4](=O)=O"),
    ("sulfonamide NH", "[NX3;H1,H2][SX4](=O)=O"),
    ("imide NH", "[#7H1]([#6]=O)[#6]=O"),
    ("phenol", "[OX2H1]c"),
    ("amide / lactam NH", "[#7;H1,H2][#6]=O"),
    ("aromatic NH", "[nH]"),
    ("thiol", "[SX2H1]"),
    ("alcohol", "[OX2H1][CX4]"),
    ("enol / other OH", "[OX2H1]"),
)
BASE_GROUPS = (
    ("guanidine", "[#7]~[CX3](~[#7])~[#7]"),
    ("amidine", "[#7;!$(N~[#6]=O)]~[CX3;!$(C=O);!$(C[OX2,SX2])]~[#7]"),
    ("aliphatic amine", "[NX3;!$(N-a);!$(N-[#6]=[#6,#7,#8,#16]);!$(N-[#6]#*);!$(N-[!#6;!#1])]"),
    ("aniline-type N", "[NX3]a"),
    ("imidazole-type N", "[nX2;r5]:[#6]:[nX3]"),
    ("pyridine-type N", "[nX2]"),
    ("imine N", "[NX2]=[#6]"),
)
# fallback when MolGpKa cannot run: typical pKa of the group (textbook values)
FALLBACK_SITES = (
    ("sulfonic acid / sulfate", "[OX2H1][SX4](=O)=O", "acid", -1.0),
    ("phosphate / phosphonate", "[OX2H1][PX4]=O", "acid", 1.5),
    ("carboxylic acid", "[OX2H1][CX3]=O", "acid", 4.5),
    ("tetrazole NH", "[#7H1]1:[#7]:[#7]:[#7]:[#6]:1", "acid", 4.9),
    ("tetrazole NH", "[#7H1]1:[#7]:[#7]:[#6]:[#7]:1", "acid", 4.9),
    ("acylsulfonamide NH", "[NX3H1]([CX3]=O)[SX4](=O)=O", "acid", 4.5),
    ("imide NH", "[#7H1]([#6]=O)[#6]=O", "acid", 9.5),
    ("sulfonamide NH", "[NX3;H1,H2][SX4](=O)=O", "acid", 10.0),
    ("phenol", "[OX2H1]c", "acid", 10.0),
    ("guanidine", "[NX2]=[CX3]([NX3])[NX3]", "base", 13.0),
    ("amidine", "[NX2;!$(N~[#6]=O)]=[CX3;!$(C[OX2,SX2]);!$(C([NX3])[NX3])][NX3;!$(N[#6]=O)]", "base", 12.0),
    ("aliphatic amine", "[NX3;!$(N-a);!$(N-[#6]=[#6,#7,#8,#16]);!$(N-[#6]#*);!$(N-[!#6;!#1])]", "base", 9.5),
    ("imidazole-type N", "[nX2;r5]:[#6]:[nX3]", "base", 7.0),
    ("pyridine-type N", "[nX2;r6]", "base", 5.0),
    ("aniline-type N", "[NX3;!$(N[#6]=[O,S,N]);!$(N[SX4])]c", "base", 4.0),
)
# groups that accept a proton in the ESI droplet / gas phase without being bases in water
ACCEPTORS = (
    ("amide / lactam / urea C=O", "[#6X3](=O)[#7]"),
    ("aromatic N", "[n]"),
    ("conjugated carbonyl", "[#6X3;!$([#6][OX2H1]);!$([#6][OX1-])](=O)[$([#6X3]=[#6X3]),$(c),$([#6X3]=O)]"),
    ("aromatic carbonyl (pyrone / coumarin / flavone)", "c=O"),
    ("phosphoryl P=O", "[PX4;!$(P[OX2H1]);!$(P[OX1-])]=O"),
    ("sulfoxide", "[SX3]=O"),
    ("N-oxide", "[#7+][OX1-]"),
    ("imine / oxime / hydrazone N", "[NX2;!a]=[#6,#7]"),
    ("electron-rich (meta-dioxygenated) aromatic ring", "[OX2]c:c:c[OX2]"),
)
WATER_LOSS = (
    ("tertiary alcohol", "[OX2H1][CX4;H0]([#6])([#6])[#6]"),
    ("benzylic alcohol", "[OX2H1][CX4]c"),
    ("allylic alcohol", "[OX2H1][CX4][CX3]=[CX3]"),
    ("β-hydroxy ketone / aldehyde", "[OX2H1][CX4][CX4][CX3;!$(C[OX2,OX1-,#7])]=O"),
)
NH3_LOSS = ("[NX3H2][CX4]c", "[NX3H2][CX4][CX4]c", "[NX3H2][CX4][CX3]=O")
CO2_LOSS = ("[CX3](=O)([OX2H1])[CX4]c", "[CX3](=O)([OX2H1])c:c[OX2H1]", "[CX3](=O)([OX2H1])[CX4][CX3]=O")
STRONG_ACID_SMARTS = ("[OX2H1][SX4](=O)=O", "[OX2H1][PX4]=O", "[OX2H1][CX3]=O")   # inner-salt anion, most acidic first

# lipid headgroups on the phosphate (the classifier's own head names mix PS / PE and PA / PG, so they are re-detected here)
PL_HEADS = (
    ("PC", "[PX4](=O)(~[OX1,OX2])O[CH2][CH2][NX4+]([CH3])([CH3])[CH3]"),
    ("PS", "[PX4](=O)(~[OX1,OX2])O[CH2][CX4H1]([NX3])[CX3](=O)[OX2H1,OX1-]"),
    ("PE", "[PX4](=O)(~[OX1,OX2])O[CH2][CH2][NX3;!$(N[#6]=O)]"),
    ("PI", "[PX4](=O)(~[OX1,OX2])O[CX4;R]1[CX4;R]([OX2])[CX4;R]([OX2])[CX4;R]([OX2])[CX4;R]([OX2])[CX4;R]1[OX2]"),
    ("PA", "[PX4](=O)([OX2H1,OX1-])([OX2H1,OX1-])O[CX4]"),
    ("PG", "[PX4](=O)(~[OX1,OX2])(O[CH2][CX4H1]([OX2])[CH2][OX2])O[CH2][CX4H1]([OX2])[CH2][OX2]"),
)
_GLYCEROL = "[CX4;!R]([OX2])[CX4;!R]([OX2])[CX4;!R][OX2]"
_AMIDE_NH = "[CX3](=O)[NX3H1][CX4]"
_STERYL_ESTER = "[CX3](=O)[OX2][CX4;R]"


@lru_cache(maxsize=None)
def _patt(smarts):
    p = Chem.MolFromSmarts(smarts)
    if p is None:
        raise ValueError(f"bad SMARTS {smarts}")
    return p


def _has(mol, smarts):
    return mol.HasSubstructMatch(_patt(smarts))


def _first_atoms(mol, smarts):
    return {m[0] for m in mol.GetSubstructMatches(_patt(smarts))}


# ---------------------------------------------------------------- sites
@dataclass
class Site:
    atom: int
    label: str        # element + 1-based atom number, as on the Separation page
    kind: str         # "acid" | "base"
    pka: float
    group: str


def _group(mol, atom, kind):
    for name, smarts in (ACID_GROUPS if kind == "acid" else BASE_GROUPS):
        if atom in _first_atoms(mol, smarts):
            return name
    return "acidic site" if kind == "acid" else "basic site"


def _fallback_sites(mol):
    seen, out = set(), []
    for name, smarts, kind, pka in FALLBACK_SITES:
        for a in sorted(_first_atoms(mol, smarts)):
            if (a, kind) not in seen:
                seen.add((a, kind))
                out.append(Site(a, f"{mol.GetAtomWithIdx(a).GetSymbol()}{a + 1}", kind, pka, name))
    return out


def find_sites(mol):
    """(sites, method, notes): MolGpKa sites of the neutral molecule, or the functional-group fallback when MolGpKa cannot run."""
    try:
        sites, dropped = pe.drop_nonbasic_sites(mol, pe.predict_sites(mol))
        out = [Site(s.atom, s.label, s.kind, s.pka, _group(mol, s.atom, s.kind)) for s in sites]
        notes = [f"{s.label} ({why}): MolGpKa base site dropped as not basic" for s, why in dropped]
        return out, "MolGpKa", notes
    except pe.PkaError:
        why = f"MolGpKa weights are not installed ({pe.pka_setup.SETUP_COMMAND})"
    except Exception as exc:  # noqa: BLE001 - odd structure for the network: use the functional groups
        why = f"MolGpKa failed on this structure ({exc})"
    return _fallback_sites(mol), "functional groups", [
        f"pKa fallback: {why}; typical pKa values of the functional groups are used instead (polyol OH is not counted as acidic)."]


def protonated(pka, ph):
    return 1.0 / (1.0 + 10.0 ** (ph - pka))


def deprotonated(pka, ph):
    return 1.0 / (1.0 + 10.0 ** (pka - ph))


def _pct(f):
    return "< 1 %" if f < 0.01 else ("> 99 %" if f > 0.99 else f"{100 * f:.0f} %")


# ---------------------------------------------------------------- structure profile
@dataclass
class Profile:
    sites: list
    method: str
    bases: list            # pKa >= base_pka, strongest first
    weak_bases: list       # base sites below base_pka
    strong_acids: list     # most acidic first
    weak_acids: list
    very_weak_acids: list
    acceptors: list        # names of proton-accepting groups
    permanent: bool
    fixed_charge: int
    inner_salt: bool
    n_o: int
    n_n: int
    polar: bool            # has N, O, P or S
    aromatic: bool
    donors: int            # alcohol OH + amide NH
    mw: float
    water: list            # in-source water-loss groups
    sugar: bool
    nh3: bool
    co2: bool


def profile(mol, sites, method, prep, r=RULES):
    bases = sorted([s for s in sites if s.kind == "base" and s.pka >= r.base_pka], key=lambda s: -s.pka)
    weak_bases = sorted([s for s in sites if s.kind == "base" and s.pka < r.base_pka], key=lambda s: -s.pka)
    acids = sorted([s for s in sites if s.kind == "acid"], key=lambda s: s.pka)
    strong = [s for s in acids if s.pka < r.strong_acid_pka]
    weak = [s for s in acids if r.strong_acid_pka <= s.pka < r.weak_acid_pka]
    very_weak = [s for s in acids if r.weak_acid_pka <= s.pka < r.very_weak_acid_pka]
    acceptors = [name for name, smarts in ACCEPTORS if _has(mol, smarts)]
    acceptors += [f"weakly basic {s.group} {s.label} (pKa {s.pka:.1f})" for s in weak_bases]
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    alcohols = len(_first_atoms(mol, "[OX2H1][CX4]"))
    donors = alcohols + len(_first_atoms(mol, "[#7;H1,H2][#6]=O"))
    sugar = clf.detect_glycosyl(mol)[0] > 0 or alcohols >= 3
    inner = prep.permanent and prep.fixed_charge > 0 and any(s.pka < r.inner_salt_pka for s in acids)
    return Profile(
        sites=sites, method=method, bases=bases, weak_bases=weak_bases, strong_acids=strong, weak_acids=weak, very_weak_acids=very_weak,
        acceptors=acceptors, permanent=prep.permanent, fixed_charge=prep.fixed_charge, inner_salt=inner,
        n_o=symbols.count("O"), n_n=symbols.count("N"), polar=any(s in ("N", "O", "P", "S") for s in symbols),
        aromatic=any(a.GetIsAromatic() for a in mol.GetAtoms()), donors=donors, mw=Descriptors.MolWt(mol),
        water=[name for name, smarts in WATER_LOSS if _has(mol, smarts)], sugar=sugar,
        nh3=any(_has(mol, s) for s in NH3_LOSS), co2=any(_has(mol, s) for s in CO2_LOSS))


# ---------------------------------------------------------------- candidates
class Candidates:
    """{(polarity, ion): [score, reason]}; put keeps the higher score."""

    def __init__(self):
        self.d = {}

    def put(self, pol, ion, score, reason):
        if score <= 0:
            return
        old = self.d.get((pol, ion))
        if old is None or score > old[0]:
            self.d[(pol, ion)] = [float(score), reason]

    def score(self, pol, ion):
        return self.d.get((pol, ion), [0.0])[0]

    def scale(self, pol, factor, cap, note):
        for (p, _), v in self.d.items():
            if p == pol:
                v[0] = min(v[0] * factor, cap)
                v[1] += f"; {note}"


def _blend(f, r):
    return r.ph_floor + (1 - r.ph_floor) * f


def _anion_score(p, r):
    return r.anion_donors2 if p.donors >= 2 else (r.anion_donors1 if p.donors == 1 else r.anion_donors0)


def small_molecule(p, add, r=RULES):
    """Candidate ions from ionizable sites and functional groups -> (Candidates, rule hits)."""
    c, hits = Candidates(), []
    ph = add.ph
    na_factor = r.na_with_ammonium if add.ammonium else (r.na_with_acid if add.acid else r.na_with_none)
    weak_factor = r.weak_with_ammonium if add.ammonium else (r.weak_with_none if add.name == "None" else 1.0)
    cation = p.permanent and p.fixed_charge > 0 and not p.inner_salt
    anion_only = p.permanent and p.fixed_charge < 0

    # ---- positive
    cls = "base" if p.bases else "weak" if p.acceptors else "polar" if p.polar else "none"
    if cation:
        if p.fixed_charge == 1:
            c.put("+", "[M]+", r.h_base, "permanent cation: the molecule is already an ion, so no [M+H]+")
        else:
            c.put("+", "[M]2+", r.h_base, f"permanent {p.fixed_charge}+ cation: seen as the multiply charged ion itself")
            c.put("+", "[M-H]+", r.dication_mh, "dication loses H+ to a singly charged ion")
        hits.append("Permanent cation (quaternary N+, sulfonium, phosphonium ...) without an acidic group: [M]+ only; nothing in negative mode.")
    elif not anion_only:
        if p.inner_salt:
            h, why = r.h_base, "permanent cation with an acidic group (inner salt): the acid takes H+, so it is seen as [M+H]+ (the cation itself)"
            hits.append("Inner salt (betaine / carnitine type): M is the zwitterion; [M+H]+ in positive, only the anion adduct in negative.")
        elif cls == "base":
            b = p.bases[0]
            f = protonated(b.pka, ph)
            h = r.h_base * _blend(f, r)
            why = f"basic {b.group} {b.label} (pKa {b.pka:.1f}; {_pct(f)} protonated at pH {ph:g})"
            hits.append(f"Basic site (pKa ≥ {r.base_pka:g}): [M+H]+ dominates; Na+ adducts stay minor.")
        elif cls == "weak":
            h = r.h_weak * weak_factor
            why = f"no basic site with pKa ≥ {r.base_pka:g}, but a proton-accepting {p.acceptors[0]} (protonated in the droplet / gas phase)"
            hits.append(f"No basic site, proton-accepting group(s): {', '.join(p.acceptors[:4])} → [M+H]+ still the main positive ion.")
        elif cls == "polar":
            h, why = r.h_polar, "no basic or proton-accepting group: protonation is weak"
            hits.append("Neutral, O-only polar molecule: cation adducts ([M+Na]+, [M+NH4]+ with ammonium) rank above [M+H]+.")
        else:
            h, why = r.h_none, "no N, O, P or S: ESI barely ionizes it"
            hits.append("No N / O / P / S: poor ESI response in both polarities; use APCI or APPI (or GC-MS).")
        c.put("+", "[M+H]+", h, why)
        if cls == "none" and p.aromatic:
            c.put("+", "[M]+", r.radical, "radical cation [M]+• by charge transfer / oxidation in the source (APPI / APCI do this much better)")
        # alkali adducts: O / N lone pairs chelate Na+; more O, more Na+ (Kruve & Kaupmees 2017)
        if p.n_o + p.n_n > 0:
            if p.inner_salt:
                na, why = r.quat_na, "Na+ replaces the acidic H of the inner salt"
            elif cls == "base":
                na, why = r.na_base, "Na+ adduct, minor next to [M+H]+ of a base"
            elif cls == "weak":
                na, why = r.na_weak, "Na+ binds the carbonyl / heteroatom lone pairs"
            elif p.n_o >= r.polyoxy_min_o:
                na, why = r.na_polyoxy, f"{p.n_o} O atoms chelate Na+ (polyoxygenated neutral: sugars, glycosides, esters)"
            else:
                na, why = r.na_polar, f"neutral with {p.n_o} O atom(s): Na+ adduct"
            if add.ammonium:
                why += "; NH4+ from the additive competes"
            elif add.name == "None":
                why += "; no additive cation, so Na+ / K+ rise"
            c.put("+", "[M+Na]+", na * na_factor, why)
            c.put("+", "[M+K]+", na * na_factor * r.k_factor, "K+ adduct, usually a fraction of [M+Na]+")
        # ammonium adducts of neutrals (proton transfer wins for bases)
        if add.ammonium and p.n_o > 0 and cls != "base" and not p.inner_salt:
            if cls == "weak":
                nh4, why = r.nh4_weak, "NH4+ from the additive binds the carbonyl O"
            elif p.n_o >= r.polyoxy_min_o:
                nh4, why = r.nh4_polyoxy, f"NH4+ from the additive binds the {p.n_o} O atoms of a neutral molecule (usually its main positive ion)"
            else:
                nh4, why = r.nh4_polar, f"NH4+ from the additive binds the {p.n_o} O atom(s)"
            c.put("+", "[M+NH4]+", nh4, why)
        # in-source losses
        h = c.score("+", "[M+H]+")
        if p.water:
            c.put("+", "[M+H-H2O]+", max(r.water_min, r.water_frac * h), f"{p.water[0]}: loses water in the source")
            hits.append(f"In-source water loss ({', '.join(p.water)}) → [M+H-H2O]+.")
        elif p.sugar:
            c.put("+", "[M+H-H2O]+", r.water_sugar, "sugar / polyol: [M+H]+ loses water in the source")
            hits.append("Sugar / polyol → [M+H-H2O]+ (low).")
        if p.nh3 and (cls == "base" or p.inner_salt):
            c.put("+", "[M+H-NH3]+", r.nh3_loss, "primary amine next to an aromatic ring / carbonyl: loses NH3 in the source")
            hits.append("Primary amine next to an aromatic ring or carbonyl → [M+H-NH3]+ (low).")
        n2 = len({b.atom for b in p.bases if b.pka >= r.multi_base_pka})
        if n2 >= 2 and p.mw >= r.z2_min_mw:
            z2 = r.z2_high if (n2 >= 3 and p.mw >= r.z2_high_mw) else r.z2_mid if p.mw >= r.z2_mid_mw else r.z2_low
            c.put("+", "[M+2H]2+", z2, f"{n2} basic sites with pKa ≥ {r.multi_base_pka:g}, MW {p.mw:.0f}")
            hits.append("Several strong bases and enough size → [M+2H]2+.")
        # an acid in positive mode never gets a high-likelihood ion
        if p.strong_acids and cls != "base" and not p.inner_salt:
            a = p.strong_acids[0]
            factor = r.acid_pos_weak if cls == "weak" else r.acid_pos
            c.scale("+", factor, r.high - 1, f"scaled down: acidic {a.group} (negative mode is the better choice)")
            hits.append(f"Consistency: acidic {a.group} {a.label} (pKa {a.pka:.1f}) and no basic site → positive ions × {factor:g}, never high.")

    # ---- negative
    if cation:
        pass
    elif anion_only:
        c.put("-", "[M]-", r.mh_strong, "permanent anion: the molecule is already an ion")
    elif not p.polar:
        pass                             # hydrocarbons: no deprotonation, nothing for an anion to bind
    else:
        if not p.inner_salt:
            if p.strong_acids:
                a = p.strong_acids[0]
                f = deprotonated(a.pka, ph)
                s, why = r.mh_strong * _blend(f, r), f"acidic {a.group} {a.label} (pKa {a.pka:.1f}; {_pct(f)} ionized at pH {ph:g})"
                hits.append(f"Acidic site (pKa < {r.strong_acid_pka:g}): [M-H]- dominates; anion adducts stay minor.")
            elif p.weak_acids:
                a = p.weak_acids[0]
                f = deprotonated(a.pka, ph)
                s, why = r.mh_weak * _blend(f, r), f"weakly acidic {a.group} {a.label} (pKa {a.pka:.1f}; {_pct(f)} ionized at pH {ph:g}; deprotonates in ESI−)"
                hits.append(f"Weakly acidic site (pKa {r.strong_acid_pka:g}–{r.weak_acid_pka:g}, e.g. phenol): [M-H]- in negative mode.")
            elif p.very_weak_acids:
                a = p.very_weak_acids[0]
                s, why = r.mh_very_weak, f"only very weakly acidic {a.group} {a.label} (pKa {a.pka:.1f}): [M-H]- forms, but weakly"
            else:
                s, why = r.mh_none, "no acidic site: deprotonation is unlikely"
                if p.bases:
                    hits.append("Base without an acidic site: negative mode is a poor choice.")
            c.put("-", "[M-H]-", s, why)
        anion_name = {"HCOO": "formate", "CH3COO": "acetate"}.get(add.anion, "")
        if add.anion:
            ion = f"[M+{add.anion}]-"
            if p.inner_salt:
                s, why = r.quat_anion, f"inner salt has no free acidic H: seen as the {anion_name} adduct"
            elif p.strong_acids:
                s, why = r.anion_with_acid, f"{anion_name} adduct, minor: acids lose H+ instead"
            else:
                s = _anion_score(p, r)
                why = f"{p.donors} OH / NH group(s) hydrogen-bond the {anion_name} of the additive" if p.donors else f"{anion_name} adduct of a neutral without OH / NH donors, weak"
                if p.weak_acids:
                    s *= r.anion_with_weak_acid
                    why += f"; the {p.weak_acids[0].group} deprotonates instead"
                elif p.bases:
                    s *= r.anion_with_base
                    why += "; bases form anion adducts poorly"
                if not p.weak_acids and s >= r.medium:
                    hits.append(f"Neutral with OH / NH donors in negative mode → [M+{add.anion}]- ({anion_name} from the additive).")
            c.put("-", ion, s, why)
        if not p.strong_acids and (p.donors or p.inner_salt):
            s = r.quat_anion if p.inner_salt else _anion_score(p, r)
            if p.weak_acids:
                s *= r.anion_with_weak_acid
            elif p.bases:
                s *= r.anion_with_base
            s *= r.cl_factor_none if add.name == "None" else r.cl_factor
            c.put("-", "[M+Cl]-", s, "chloride traces (solvents, salts, matrix) give [M+Cl]- with neutral molecules")
        if p.co2 and p.strong_acids and not p.inner_salt:
            c.put("-", "[M-H-CO2]-", r.co2_loss, "aryl-acetic / salicylic / β-keto acid: decarboxylates in the source")
            hits.append("Easily decarboxylated acid → [M-H-CO2]- (low).")
        if p.strong_acids and not p.inner_salt and "β-hydroxy ketone / aldehyde" in p.water:
            c.put("-", "[M-H-H2O]-", r.water_min, "β-hydroxy ketone: loses water in the source")
        n2 = len({a.atom for a in p.strong_acids})
        if n2 >= 2 and p.mw >= r.z2_min_mw:
            s = r.z2_acid_mid if p.mw >= r.z2_acid_mid_mw else r.z2_low
            c.put("-", "[M-2H]2-", s, f"{n2} acidic sites with pKa < {r.strong_acid_pka:g}, MW {p.mw:.0f}")
            hits.append("Several strong acid sites and enough size → [M-2H]2-.")
    return c, hits


# ---------------------------------------------------------------- lipids
@dataclass
class Lipid:
    code: str          # rule key in Rules.lipid_ions
    name: str          # display name (LPC, PE, HexCer ...)
    category: str      # the classifier's lipid class
    detail: str


def _chain_on_o(mol, o, backbone_c):
    """True if glycerol O ``o`` carries an acyl (ester) or a long alkyl / alkenyl (ether) chain."""
    for nb in mol.GetAtomWithIdx(o).GetNeighbors():
        if nb.GetIdx() == backbone_c:
            continue
        if nb.GetSymbol() != "C":
            return False
        if any(b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(nb).GetSymbol() == "O" for b in nb.GetBonds()):
            return True
        seen, todo = {nb.GetIdx()}, [nb]            # ether: count the acyclic carbons reachable without passing the O
        while todo:
            a = todo.pop()
            for x in a.GetNeighbors():
                if x.GetSymbol() == "C" and not x.IsInRing() and x.GetIdx() not in seen:
                    seen.add(x.GetIdx())
                    todo.append(x)
        return len(seen) >= 6 and not nb.IsInRing()
    return False


def glycerol_chains(mol):
    """Most acyl / alkyl chains on one glycerol unit (TG 3, DG / PC 2, MG / LPC 1)."""
    best = 0
    for m in mol.GetSubstructMatches(_patt(_GLYCEROL)):
        best = max(best, sum(_chain_on_o(mol, o, c) for c, o in ((m[0], m[1]), (m[2], m[3]), (m[4], m[5]))))
    return best


def lipid_subclass(mol):
    """(Lipid with class-specific ion rules or None, the classifier's lipid class or None)."""
    f = clf.extract_features(Chem.MolToSmiles(mol))
    cat = f.lipid_class if f.valid else None
    if not cat:
        return None, None
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    n_p = symbols.count("P")
    strong_acid = any(_has(mol, s) for s in STRONG_ACID_SMARTS)
    chains = glycerol_chains(mol)
    lyso = "L" if chains == 1 else ""
    if cat == "glycerophospholipid":
        head = next((code for code, smarts in PL_HEADS if _has(mol, smarts)), None)
        if head == "PI" and n_p > 1:
            return Lipid("PIP", "PIP", cat, "phosphoinositol head with extra phosphate(s)"), cat
        if head == "PG" and n_p >= 2:
            return Lipid("CL", "CL", cat, "two phosphates bridged by glycerol"), cat
        if head:
            return Lipid(head, lyso + head, cat, f"{head} headgroup, {chains} chain(s)"), cat
        return None, cat
    if cat == "sphingolipid":
        if _has(mol, PL_HEADS[0][1]):
            return Lipid("PC", "SM", cat, "phosphocholine on a sphingoid base"), cat
        if not _has(mol, _AMIDE_NH) or n_p or strong_acid:
            return None, cat             # sphingoid base, S1P, CerP, sulfatide, ganglioside: small-molecule rules
        if f.n_glycosyl:
            return Lipid("Cer", "HexCer", cat, f"N-acyl sphingoid base with {f.n_glycosyl} sugar(s), no acidic group"), cat
        return Lipid("Cer", "Cer", cat, "N-acyl sphingoid base"), cat
    if cat == "glycerolipid":
        if n_p or strong_acid:
            return None, cat             # SQDG and the like: small-molecule acid rules
        if f.n_glycosyl:
            return Lipid("GlyDG", "MGDG / DGDG", cat, f"glycosyl glycerolipid ({f.n_glycosyl} sugar(s))"), cat
        code = {3: "TG", 2: "DG", 1: "MG"}.get(chains)
        return (Lipid(code, code, cat, f"{chains} chain(s) on glycerol"), cat) if code else (None, cat)
    if cat == "sterol lipid":
        if "N" in symbols or strong_acid or f.n_glycosyl:
            return None, cat
        if _has(mol, _STERYL_ESTER) and f.has_long_acyl:
            return Lipid("CE", "CE / steryl ester", cat, "sterol esterified with a fatty acid"), cat
        carbonyl = _has(mol, "[CX3]=O")
        aromatic = any(a.GetIsAromatic() for a in mol.GetAtoms())
        if f.isoprenoid_type == "sterol" and not carbonyl and not aromatic and symbols.count("O") <= 2:
            return Lipid("ST", "free sterol", cat, "sterol with OH only (no carbonyl, acid or phenol)"), cat
        return None, cat                 # hormones, bile acids, oxysterols with C=O: small-molecule rules
    if cat == "fatty acyl":
        if _has(mol, "[OX2H1][CX3]=O") and not any(s in symbols for s in ("N", "P", "S")):
            return Lipid("FA", "FA", cat, "free fatty acid"), cat
        return None, cat
    return None, cat


def lipid_candidates(lipid, add, r=RULES):
    c = Candidates()
    for pol, ion, s_amm, s_no, why in r.lipid_ions[lipid.code]:
        s = s_amm if add.ammonium else s_no
        if ion in ("[M+A]-", "[M-CH3]-") and not add.anion:
            continue
        if ion == "[M+A]-":
            ion = f"[M+{add.anion}]-"
        if ion == "[M+Cl]-" and add.name == "None":
            s *= r.lipid_cl_none
        c.put(pol, ion, s, f"{lipid.name}: {why}")
    return c


# ---------------------------------------------------------------- M, m/z
def neutral_m(prep, p_sites, r=RULES):
    """(formula of M without a charge sign, mol of M, note): the inner salt for a permanent cation with an acid, else the prepared molecule."""
    mol = prep.mol
    if not (prep.permanent and prep.fixed_charge > 0):
        return _formula(mol), mol, ""
    acid_atoms = [s.atom for s in sorted(p_sites, key=lambda s: s.pka) if s.kind == "acid" and s.pka < r.inner_salt_pka] if p_sites is not None else []
    if p_sites is None:
        for smarts in STRONG_ACID_SMARTS:
            acid_atoms += sorted(_first_atoms(mol, smarts))
    if not acid_atoms:
        return _formula(mol), mol, f"M is the cation itself (permanent charge {prep.fixed_charge:+d}, no acidic group)."
    rw = Chem.RWMol(mol)
    atom = rw.GetAtomWithIdx(acid_atoms[0])
    atom.SetFormalCharge(atom.GetFormalCharge() - 1)
    atom.SetNumExplicitHs(0)
    atom.SetNoImplicit(False)
    inner = rw.GetMol()
    Chem.SanitizeMol(inner)
    return _formula(inner), inner, "M is the inner salt (zwitterion): the permanent cation balanced by its deprotonated acid, as in PC / SM / carnitine."


def _formula(mol):
    return re.sub(r"[+-]\d*$", "", rdMolDescriptors.CalcMolFormula(mol))


def ion_mz(formula, ion):
    """(m/z, ion formula with charge) from the Theoretical m/z page's code; (None, None) when it cannot be computed."""
    from mzdesk.ms1 import msbuddy_engine as bd
    if not bd.formula_is_valid(formula):
        return None, None
    if ion == "[M]2+":
        mz1, f1 = bd.ion_mz_and_formula(formula, "[M]+")
        return ((mz1 - ELECTRON) / 2, f1[:-1] + "2+") if mz1 else (None, None)
    mz, f = bd.ion_mz_and_formula(formula, ion)
    if mz is not None and re.search(r"\]2[+-]$", ion):
        f = f[:-1] + "2" + f[-1]          # msbuddy writes only the sign
    return mz, f


ELECTRON = 0.00054857990946


# ---------------------------------------------------------------- top level
@dataclass
class Ion:
    polarity: str      # "+" | "-"
    adduct: str
    score: float
    likelihood: str
    reason: str
    mz: float = None
    formula: str = None


@dataclass
class Prediction:
    ions: list                 # chosen polarities, each ranked by score
    best: dict                 # "+" / "-" -> top Ion (or None), for both polarities
    formula: str               # M
    lipid: object              # Lipid or None
    notes: list                # rule hits
    warnings: list
    all_ions: list = field(default_factory=list)   # both polarities


def likelihood(score, r=RULES):
    return "high" if score >= r.high else "medium" if score >= r.medium else "low"


_DIMER = {"[M+H]+": "[2M+H]+", "[M+Na]+": "[2M+Na]+", "[M+NH4]+": "[2M+NH4]+", "[M-H]-": "[2M-H]-"}


def predict(smiles, polarity="Positive", additive_name="Formic acid", alkali=True, r=RULES):
    """Ranked ions for a structure, ESI polarity ("Positive" / "Negative" / "Both") and mobile-phase additive. Raises ValueError."""
    if polarity not in POLARITIES:
        raise ValueError(f"Unknown ion mode {polarity!r}; choose one of {', '.join(POLARITIES)}.")
    add = additive(additive_name, r)
    prep = pe.prepare(smiles)                      # largest fragment, neutralized; raises ValueError
    mol = prep.mol
    # keep the fragment / size / element notes; the neutralization and logD notes are about the Separation page
    warnings = [n for n in prep.notes if not n.startswith(("Input neutralized", "Permanent charge"))]
    notes = [f"Mobile phase: {add.name} → pH ≈ {add.ph:g} assumed ({add.assumption}). Solution ionization only fine-tunes the scores "
             "(× 0.8–1): in ESI, acids still give [M-H]- from acidic eluents and bases [M+H]+ from neutral ones."]

    lipid, category = lipid_subclass(mol)
    sites = None
    if lipid:
        c = lipid_candidates(lipid, add, r)
        notes.append(f"Lipid class: {category} → {lipid.name} ({lipid.detail}). Class conventions (LIPID MAPS / MS-DIAL / LipidSearch practice) "
                     "are used instead of the pKa rules.")
        if lipid.code == "ST":
            notes.append("Free sterols ionize poorly by ESI; APCI (or derivatization) is the usual choice.")
        if lipid.code in ("TG", "CE"):
            notes.append(f"{lipid.name} is not detected in negative ESI.")
        mw = Descriptors.MolWt(mol)
        permanent = prep.permanent
    else:
        if category:
            notes.append(f"Lipid class: {category}, no class-specific rules for this structure; small-molecule rules used.")
        sites, method, site_notes = find_sites(mol)
        p = profile(mol, sites, method, prep, r)
        c, hits = small_molecule(p, add, r)
        if method == "MolGpKa":
            notes.append(f"pKa: MolGpKa, {len(sites)} ionizable site(s) after dropping non-basic N.")
        notes += site_notes
        for s in sorted(sites, key=lambda s: s.pka)[: r.max_sites_listed]:
            f = protonated(s.pka, add.ph) if s.kind == "base" else deprotonated(s.pka, add.ph)
            state = "protonated" if s.kind == "base" else "ionized"
            notes.append(f"Site {s.label}: {s.group}, {s.kind}, pKa {s.pka:.1f} ({method}); {_pct(f)} {state} at pH {add.ph:g}.")
        if len(sites) > r.max_sites_listed:
            notes.append(f"... and {len(sites) - r.max_sites_listed} more site(s).")
        notes += hits
        mw, permanent = p.mw, p.permanent

    formula, m_mol, m_note = neutral_m(prep, sites if sites is not None else None, r)
    if m_note:
        notes.append(m_note)
    notes.append(f"M = {formula} (monoisotopic {rdMolDescriptors.CalcExactMolWt(m_mol):.4f}); m/z and ion formulas use the Theoretical m/z page's code.")
    drawn, _ = pe._largest_fragment(Chem.MolFromSmiles(smiles))
    if _formula(drawn) != formula:
        warnings.append(f"The input is drawn as an ion ({_formula(drawn)}, charge {Chem.GetFormalCharge(drawn):+d}); M is taken as {formula}. "
                        "The Theoretical m/z page uses the formula as drawn: enter the neutral structure there to get the same m/z.")

    # alkali switch
    if not alkali:
        for key in list(c.d):
            if "Na" in key[1]:
                c.d[key][0] *= r.na_absent
                c.d[key][1] += "; Na / K marked absent: trace Na only"
            elif "+K]" in key[1]:
                del c.d[key]
        notes.append("Sodium / potassium marked absent: Na+ adducts cut to traces, K+ adducts removed.")

    # one dimer per polarity, for strongly ionized small molecules
    if mw <= r.dimer_max_mw and not permanent and (lipid is None or lipid.code == "FA"):
        for pol in ("+", "-"):
            mono = sorted(((v[0], ion) for (pp, ion), v in c.d.items() if pp == pol and not ion.startswith("[2M")), reverse=True)
            have = any(pp == pol and ion.startswith("[2M") for pp, ion in c.d)
            if mono and mono[0][0] >= r.high and mono[0][1] in _DIMER and not have:
                c.put(pol, _DIMER[mono[0][1]], r.dimer, f"dimer of {mono[0][1]}; grows with concentration")
        if any(ion.startswith("[2M") for _, ion in c.d):
            notes.append(f"Dimers ([2M+X]) are concentration dependent; only one per polarity is listed, for molecules up to {r.dimer_max_mw:g} Da.")

    all_ions = []
    mz_missing = False
    for (pol, ion), (score, reason) in c.d.items():
        if score < r.min_shown:
            continue
        mz, f = ion_mz(formula, ion)
        mz_missing |= mz is None
        all_ions.append(Ion(pol, ion, round(score, 1), likelihood(score, r), reason, mz, f))
    all_ions.sort(key=lambda i: (i.polarity == "-", -i.score, i.adduct))
    if mz_missing:
        warnings.append(f"m/z not computed for some ions of {formula} (supported elements: C H N O P S F Cl Br I Na K).")
    best = {pol: next((i for i in all_ions if i.polarity == pol), None) for pol in ("+", "-")}
    notes.append(_polarity_note(best))
    notes.append(f"Scores are rule-based ranks (0–100), not measured intensities: high ≥ {r.high:g}, medium ≥ {r.medium:g}. "
                 "Matrix suppression, source settings and concentration are not modelled.")
    chosen = {"Positive": ("+",), "Negative": ("-",), "Both": ("+", "-")}[polarity]
    return Prediction(ions=[i for i in all_ions if i.polarity in chosen], best=best, formula=formula, lipid=lipid,
                      notes=notes, warnings=warnings, all_ions=all_ions)


def _polarity_note(best):
    def text(i):
        return "no ion expected" if i is None else f"{i.adduct} ({i.likelihood})"
    pos, neg = best["+"], best["-"]
    sp, sn = (pos.score if pos else 0), (neg.score if neg else 0)
    verdict = "positive" if sp > sn * 1.2 else "negative" if sn > sp * 1.2 else "either polarity"
    return f"Polarity check: positive best {text(pos)}, negative best {text(neg)} → {verdict}."
