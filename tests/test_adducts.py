"""Adducts 页的测试：规则引擎（mzdesk/platform_fit/adduct_rules.py）、页面输出、与 Theoretical m/z 页的 m/z 一致性，
以及一组公认主离子的参考化合物（每个一行来源 / 理由）。

参考集的判据：排第一的离子 = 公认的主离子，而且它的 likelihood 至少是给定等级。规则没有为单个化合物调过；
原则上解释不了的真实失败用 xfail(strict=True) 记下（理由写在 MISSES 里），规则改进后会变成 XPASS 提醒更新。
小分子需要 MolGpKa 权重（python -m mzdesk.platform_fit.pka_setup），没有时跳过；脂质走类别规则，不需要权重。
"""
from functools import lru_cache
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from mzdesk.api.app import app
from mzdesk.platform_fit import adduct_rules as ar
from mzdesk.platform_fit import pka_engine as pe
from mzdesk.spec import ADDUCTS_NEG, ADDUCTS_POS, MODULES

client = TestClient(app)
NO_MOLGPKA = not Path(pe.weights_path()).is_file()
needs_weights = pytest.mark.skipif(NO_MOLGPKA, reason="MolGpKa weights not installed (python -m mzdesk.platform_fit.pka_setup)")
FA, AF, AA, AC = "Formic acid", "Ammonium formate", "Ammonium acetate", "Acetic acid"
RANK = {"low": 0, "medium": 1, "high": 2}

S = dict(
    caffeine="CN1C=NC2=C1C(=O)N(C(=O)N2C)C", metformin="CN(C)C(=N)N=C(N)N", propranolol="CC(C)NCC(O)COc1cccc2ccccc12",
    amitriptyline="CN(C)CCC=C1c2ccccc2CCc2ccccc12", nicotine="CN1CCC[C@H]1c1cccnc1", nicotinamide="NC(=O)c1cccnc1",
    lidocaine="CCN(CC)CC(=O)Nc1c(C)cccc1C", carbamazepine="NC(=O)N1c2ccccc2C=Cc2ccccc21", acetaminophen="CC(=O)Nc1ccc(O)cc1",
    diazepam="CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21", testosterone="C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@@]34C)[C@@H]1CC[C@@H]2O",
    progesterone="CC(=O)[C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@H]3CC[C@]12C", spermine="NCCCNCCCCNCCCN",
    atrazine="CCNc1nc(Cl)nc(NC(C)C)n1", tphp="O=P(Oc1ccccc1)(Oc1ccccc1)Oc1ccccc1", dbp="CCCCOC(=O)c1ccccc1C(=O)OCCCC",
    sulfamethoxazole="Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1", losartan="CCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1",
    adenosine="Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O", creatinine="CN1CC(=O)N=C1N", serotonin="NCCc1c[nH]c2ccc(O)cc12",
    catechin="O[C@@H]1Cc2c(O)cc(O)cc2O[C@@H]1c1ccc(O)c(O)c1",
    ibuprofen="CC(C)Cc1ccc(cc1)C(C)C(=O)O", salicylic="OC(=O)c1ccccc1O", diclofenac="OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl",
    palmitic="CCCCCCCCCCCCCCCC(=O)O", arachidonic=r"CCCCC/C=C\C/C=C\C/C=C\C/C=C\CCCC(=O)O",
    taurocholic="C[C@H](CCC(=O)NCCS(=O)(=O)O)[C@H]1CC[C@H]2[C@@H]3[C@H](O)C[C@@H]4C[C@H](O)CC[C@]4(C)[C@H]3C[C@H](O)[C@]12C",
    cholic="C[C@H](CCC(=O)O)[C@H]1CC[C@H]2[C@@H]3[C@H](O)C[C@@H]4C[C@H](O)CC[C@]4(C)[C@H]3C[C@H](O)[C@]12C",
    pfoa="OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    pfos="OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    citric="OC(=O)CC(O)(CC(=O)O)C(=O)O", bpa="CC(C)(c1ccc(O)cc1)c1ccc(O)cc1",
    estradiol="C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc34)[C@@H]1CC[C@@H]2O", quercetin="O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
    uric="O=c1[nH]c(=O)c2[nH]c(=O)[nH]c2[nH]1", amp="Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]1O",
    phe="N[C@@H](Cc1ccccc1)C(=O)O", glu="N[C@@H](CCC(=O)O)C(=O)O", lys="NCCCC[C@H](N)C(=O)O", glycine="NCC(=O)O", taurine="NCCS(=O)(=O)O",
    carnitine="C[N+](C)(C)C[C@@H](O)CC(=O)[O-]", betaine="C[N+](C)(C)CC(=O)[O-]", choline="C[N+](C)(C)CCO",
    acetylcholine="CC(=O)OCC[N+](C)(C)C", paraquat="C[n+]1ccc(-c2cc[n+](C)cc2)cc1",
    glucose="OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
    sucrose="OC[C@H]1O[C@@](CO)(O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@@H]1O",
    cortisol="C[C@]12CCC(=O)C=C1CC[C@@H]1[C@@H]2[C@@H](O)C[C@@]2(C)[C@H]1CC[C@]2(O)C(=O)CO",
    salicin="OC[C@H]1O[C@@H](Oc2ccccc2CO)[C@H](O)[C@@H](O)[C@@H]1O", pyrene="c1cc2ccc3cccc4ccc(c1)c2c34",
    dhea="C[C@]12CC[C@H]3[C@@H](CC=C4C[C@@H](O)CC[C@]34C)[C@@H]1CCC2=O", retinol="CC1=C(C(C)(C)CCC1)/C=C/C(C)=C/C=C/C(C)=C/CO",
    chloramphenicol="OC[C@@H](NC(=O)C(Cl)Cl)[C@H](O)c1ccc(cc1)[N+](=O)[O-]",
    # lipids (sn-1 16:0 / sn-2 18:1 where it applies)
    PC=r"CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\CCCCCCCC",
    LPC="CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(=O)([O-])OCC[N+](C)(C)C",
    PE=r"CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OCCN)OC(=O)CCCCCCC/C=C\CCCCCCCC",
    PS=r"CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OC[C@H](N)C(=O)O)OC(=O)CCCCCCC/C=C\CCCCCCCC",
    PG=r"CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OCC(O)CO)OC(=O)CCCCCCC/C=C\CCCCCCCC",
    PI=r"CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O)OC(=O)CCCCCCC/C=C\CCCCCCCC",
    PA=r"CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)O)OC(=O)CCCCCCC/C=C\CCCCCCCC",
    CL="CCCCCCCCCCCCCCCC(=O)OCC(COP(=O)(O)OCC(O)COP(=O)(O)OCC(COC(=O)CCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC",
    SM="CCCCCCCCCCCCC/C=C/[C@@H](O)[C@H](COP(=O)([O-])OCC[N+](C)(C)C)NC(=O)CCCCCCCCCCCCCCC",
    Cer="CCCCCCCCCCCCC/C=C/[C@@H](O)[C@H](CO)NC(=O)CCCCCCCCCCCCCCC",
    HexCer="CCCCCCCCCCCCC/C=C/[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)NC(=O)CCCCCCCCCCCCCCC",
    TG=r"CCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCC/C=C\CCCCCCCC)OC(=O)CCCCCCC/C=C\CCCCCCCC",
    DG=r"CCCCCCCCCCCCCCCC(=O)OCC(CO)OC(=O)CCCCCCC/C=C\CCCCCCCC", MG="CCCCCCCCCCCCCCCC(=O)OCC(O)CO",
    MGDG=r"CCCCCCCCCCCCCCCC(=O)OC[C@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)OC(=O)CCCCCCC/C=C\CCCCCCCC",
    CE=r"CCCCCCCC/C=C\CCCCCCCC(=O)O[C@H]1CC[C@@]2(C)C(=CC[C@H]3[C@@H]4CC[C@H]([C@H](C)CCCC(C)C)[C@@]4(C)CC[C@@H]32)C1",
    cholesterol="CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
)
LIPIDS = {"PC", "LPC", "PE", "PS", "PG", "PI", "PA", "CL", "SM", "Cer", "HexCer", "TG", "DG", "MG", "MGDG", "CE", "cholesterol"}

# (compound, polarity, additive, expected top ion, its minimum likelihood, source / rationale)
REFERENCE = [
    # ---- bases and other N-containing drugs: [M+H]+
    ("caffeine", "Positive", FA, "[M+H]+", "high", "standard ESI+ test compound; m/z 195.088 is its reference ion (MassBank / mzCloud)"),
    ("caffeine", "Positive", AA, "[M+H]+", "high", "protonated with ammonium acetate too; NH4+ adduct only minor"),
    ("metformin", "Positive", FA, "[M+H]+", "high", "biguanide (pKa ~12); bioanalysis uses m/z 130 [M+H]+"),
    ("propranolol", "Positive", FA, "[M+H]+", "high", "secondary amine beta-blocker; MRM 260 [M+H]+"),
    ("amitriptyline", "Positive", FA, "[M+H]+", "high", "tertiary amine TCA; MRM 278 [M+H]+"),
    ("nicotine", "Positive", AA, "[M+H]+", "high", "pyrrolidine N (pKa ~8); MRM 163 [M+H]+ with ammonium buffers"),
    ("nicotinamide", "Positive", FA, "[M+H]+", "high", "pyridine N + amide; m/z 123 [M+H]+ in metabolomics"),
    ("lidocaine", "Positive", FA, "[M+H]+", "high", "tertiary amine; MRM 235 [M+H]+"),
    ("carbamazepine", "Positive", FA, "[M+H]+", "high", "no basic N (urea), still m/z 237 [M+H]+ in every TDM method"),
    ("acetaminophen", "Positive", FA, "[M+H]+", "high", "anilide; m/z 152 [M+H]+ in ESI+"),
    ("diazepam", "Positive", FA, "[M+H]+", "high", "benzodiazepine imine N; MRM 285 [M+H]+"),
    ("testosterone", "Positive", FA, "[M+H]+", "high", "3-keto-4-ene steroid; MRM 289 [M+H]+"),
    ("progesterone", "Positive", AA, "[M+H]+", "high", "enone steroid; m/z 315 [M+H]+ even with ammonium acetate"),
    ("spermine", "Positive", FA, "[M+H]+", "high", "polyamine, MW 202: singly protonated ion dominates"),
    ("atrazine", "Positive", FA, "[M+H]+", "high", "triazine herbicide; MRM 216 [M+H]+"),
    ("tphp", "Positive", AA, "[M+H]+", "high", "organophosphate ester (P=O); m/z 327 [M+H]+ in flame-retardant methods"),
    ("dbp", "Positive", FA, "[M+H]+", "medium", "phthalate ester; m/z 279 [M+H]+ (149 fragment)"),
    ("sulfamethoxazole", "Positive", FA, "[M+H]+", "high", "aniline N; MRM 254 [M+H]+ in antibiotic residue methods"),
    ("losartan", "Positive", FA, "[M+H]+", "high", "imidazole base (tetrazole acid also present); MRM 423 [M+H]+"),
    ("losartan", "Negative", AA, "[M-H]-", "high", "tetrazole acid; m/z 421 [M-H]- in negative methods"),
    ("adenosine", "Positive", FA, "[M+H]+", "high", "purine nucleoside; m/z 268 [M+H]+"),
    ("creatinine", "Positive", FA, "[M+H]+", "high", "guanidine-type base; m/z 114 [M+H]+ in clinical methods"),
    ("serotonin", "Positive", FA, "[M+H]+", "high", "primary amine; m/z 177 [M+H]+ (160 = in-source NH3 loss)"),
    ("catechin", "Positive", FA, "[M+H]+", "high", "flavan-3-ol; m/z 291 [M+H]+ widely reported (ring protonation of the phloroglucinol A-ring)"),
    # ---- acids: [M-H]-
    ("ibuprofen", "Negative", FA, "[M-H]-", "high", "carboxylic acid; m/z 205 [M-H]- (161 = in-source CO2 loss)"),
    ("ibuprofen", "Negative", AA, "[M-H]-", "high", "same with ammonium acetate (higher response)"),
    ("salicylic", "Negative", FA, "[M-H]-", "high", "m/z 137 [M-H]-"),
    ("diclofenac", "Negative", AA, "[M-H]-", "high", "m/z 294 [M-H]- (250 = CO2 loss)"),
    ("palmitic", "Negative", AA, "[M-H]-", "high", "free fatty acid: [M-H]- (LIPID MAPS / MS-DIAL FA)"),
    ("arachidonic", "Negative", FA, "[M-H]-", "high", "m/z 303 [M-H]- in eicosanoid methods"),
    ("taurocholic", "Negative", FA, "[M-H]-", "high", "taurine-conjugated bile acid (sulfonic acid); m/z 514 [M-H]-"),
    ("cholic", "Negative", AA, "[M-H]-", "high", "unconjugated bile acid; m/z 407 [M-H]-"),
    ("pfoa", "Negative", AA, "[M-H]-", "high", "EPA 537 / 1633 PFAS methods: m/z 413 [M-H]-"),
    ("pfos", "Negative", AA, "[M-H]-", "high", "EPA PFAS methods: m/z 499 [M-H]-"),
    ("citric", "Negative", FA, "[M-H]-", "high", "TCA acid; m/z 191 [M-H]-"),
    ("bpa", "Negative", AA, "[M-H]-", "high", "phenol; m/z 227 [M-H]- is the standard BPA ion"),
    ("estradiol", "Negative", AA, "[M-H]-", "medium", "phenolic estrogen; underivatized methods use m/z 271 [M-H]-"),
    ("quercetin", "Negative", FA, "[M-H]-", "high", "flavonol; m/z 301 [M-H]-"),
    ("uric", "Negative", FA, "[M-H]-", "high", "acidic N-H (pKa 5.4); m/z 167 [M-H]-"),
    ("amp", "Negative", AA, "[M-H]-", "high", "nucleotide phosphate; m/z 346 [M-H]-"),
    ("amp", "Positive", FA, "[M+H]+", "high", "adenine base: m/z 348 [M+H]+ in positive mode"),
    # ---- zwitterions and permanent cations
    ("phe", "Positive", FA, "[M+H]+", "high", "amino acid; m/z 166 [M+H]+"),
    ("phe", "Negative", AA, "[M-H]-", "high", "amino acid; m/z 164 [M-H]-"),
    ("glu", "Negative", FA, "[M-H]-", "high", "acidic amino acid; m/z 146 [M-H]-"),
    ("glu", "Positive", FA, "[M+H]+", "high", "m/z 148 [M+H]+"),
    ("lys", "Positive", FA, "[M+H]+", "high", "basic amino acid, MW 146: m/z 147 [M+H]+ (no 2+)"),
    ("glycine", "Negative", AA, "[M-H]-", "high", "m/z 74 [M-H]-"),
    ("taurine", "Negative", FA, "[M-H]-", "high", "sulfonic acid zwitterion; m/z 124 [M-H]-"),
    ("carnitine", "Positive", FA, "[M+H]+", "high", "inner salt: m/z 162 = [M+H]+ of the zwitterion (= the cation)"),
    ("betaine", "Positive", FA, "[M+H]+", "high", "inner salt: m/z 118 [M+H]+"),
    ("choline", "Positive", FA, "[M]+", "high", "quaternary ammonium: m/z 104 [M]+"),
    ("acetylcholine", "Positive", AA, "[M]+", "high", "quaternary ammonium: m/z 146 [M]+"),
    ("paraquat", "Positive", FA, "[M]2+", "high", "bipyridinium dication: m/z 93 [M]2+ (185 [M-H]+ second)"),
    # ---- neutrals, sugars, glycosides
    ("glucose", "Positive", FA, "[M+Na]+", "high", "sugars: Na+ adduct (m/z 203) dominates without ammonium"),
    ("glucose", "Positive", AA, "[M+NH4]+", "high", "with ammonium: m/z 198 [M+NH4]+"),
    ("glucose", "Negative", FA, "[M+HCOO]-", "high", "neutral sugar with formate: m/z 225 [M+HCOO]-"),
    ("sucrose", "Positive", FA, "[M+Na]+", "high", "m/z 365 [M+Na]+"),
    ("sucrose", "Negative", AA, "[M+CH3COO]-", "high", "with acetate: m/z 401 [M+CH3COO]-"),
    ("salicin", "Positive", AA, "[M+NH4]+", "high", "phenolic glycoside, no basic group: NH4+ adduct with ammonium"),
    ("cortisol", "Positive", FA, "[M+H]+", "high", "enone corticosteroid; MRM 363 [M+H]+"),
    ("cortisol", "Negative", FA, "[M+HCOO]-", "high", "corticosteroids give formate / acetate adducts in ESI-"),
    # ---- lipids (class conventions: LIPID MAPS, MS-DIAL, LipidSearch, Hsu & Turk)
    ("PC", "Positive", FA, "[M+H]+", "high", "PC 34:1 m/z 760.585 [M+H]+"),
    ("PC", "Negative", FA, "[M+HCOO]-", "high", "PC in ESI- with formate: [M+HCOO]-, not [M-H]-"),
    ("PC", "Negative", AA, "[M+CH3COO]-", "high", "PC in ESI- with acetate: [M+CH3COO]-"),
    ("LPC", "Positive", AF, "[M+H]+", "high", "LPC 16:0 m/z 496.340 [M+H]+"),
    ("SM", "Positive", FA, "[M+H]+", "high", "SM d18:1/16:0 m/z 703.575 [M+H]+"),
    ("SM", "Negative", AA, "[M+CH3COO]-", "high", "SM in ESI- with acetate: [M+CH3COO]-"),
    ("PE", "Negative", FA, "[M-H]-", "high", "PE: [M-H]- in negative mode"),
    ("PE", "Positive", FA, "[M+H]+", "high", "PE: [M+H]+ in positive mode"),
    ("PS", "Negative", AA, "[M-H]-", "high", "PS: [M-H]-"),
    ("PG", "Negative", AA, "[M-H]-", "high", "PG: [M-H]-"),
    ("PG", "Positive", AA, "[M+NH4]+", "medium", "PG positive with ammonium: [M+NH4]+ (MS-DIAL)"),
    ("PI", "Negative", FA, "[M-H]-", "high", "PI: [M-H]-"),
    ("PA", "Negative", AA, "[M-H]-", "high", "PA: [M-H]-"),
    ("CL", "Negative", AA, "[M-H]-", "high", "CL: [M-H]- (with [M-2H]2-)"),
    ("TG", "Positive", AA, "[M+NH4]+", "high", "TG with ammonium: [M+NH4]+"),
    ("TG", "Positive", FA, "[M+Na]+", "high", "TG without ammonium: [M+Na]+"),
    ("DG", "Positive", AF, "[M+NH4]+", "high", "DG with ammonium: [M+NH4]+"),
    ("MG", "Positive", AA, "[M+NH4]+", "high", "MG with ammonium: [M+NH4]+"),
    ("MGDG", "Positive", AA, "[M+NH4]+", "high", "galactolipid with ammonium: [M+NH4]+"),
    ("MGDG", "Negative", AF, "[M+HCOO]-", "high", "galactolipid in ESI-: [M+HCOO]-"),
    ("Cer", "Positive", FA, "[M+H]+", "high", "Cer d18:1/16:0: [M+H]+ (with [M+H-H2O]+)"),
    ("Cer", "Negative", AF, "[M+HCOO]-", "high", "Cer in ESI- with formate: [M+HCOO]-"),
    ("HexCer", "Positive", FA, "[M+H]+", "high", "HexCer: [M+H]+"),
    ("CE", "Positive", AA, "[M+NH4]+", "high", "CE with ammonium: [M+NH4]+"),
    ("cholesterol", "Positive", AA, "[M+H-H2O]+", "medium", "free cholesterol in ESI / APCI: m/z 369.35 [M+H-H2O]+"),
]
# genuine misses: the rule principle does not explain the known ion; recorded, not tuned away
MISSES = {
    ("dhea", "Positive", FA): ("[M+H-H2O]+", "DHEA shows m/z 271 [M+H-H2O]+ (homoallylic 3β-OH next to Δ5); the rules only flag allylic / benzylic / "
                                             "tertiary OH and rank [M+Na]+ first for an O-only neutral"),
    ("retinol", "Positive", FA): ("[M+H-H2O]+", "retinol is seen as m/z 269 [M+H-H2O]+ (very stable polyene carbocation); the water-loss rule "
                                                "scales with the weak [M+H]+ of an O-only neutral, so [M+Na]+ ranks first"),
    ("chloramphenicol", "Negative", FA): ("[M-H]-", "residue methods use m/z 321 [M-H]-: the dichloroacetamide NH is far more acidic in the gas phase "
                                                    "than its solution pKa suggests; the rules rank the formate adduct (3 OH / NH donors) first"),
}


@lru_cache(maxsize=None)
def _predict(name, polarity, additive, alkali=True):
    return ar.predict(S[name], polarity, additive, alkali)


def _check(name, polarity, additive, expected, min_lk):
    res = _predict(name, polarity, additive)
    assert res.ions, f"{name}: no ion predicted"
    got = res.ions[0]
    assert got.adduct == expected, f"{name} {polarity} {additive}: top {got.adduct}, expected {expected}; " + str([(i.adduct, i.score) for i in res.ions[:4]])
    assert RANK[got.likelihood] >= RANK[min_lk], f"{name}: {expected} is only {got.likelihood}"


def _params():
    out = []
    for name, pol, add, exp, lk, why in REFERENCE:
        marks = [] if name in LIPIDS else [needs_weights]
        out.append(pytest.param(name, pol, add, exp, lk, id=f"{name}-{pol[:3]}-{add.replace(' ', '_')}", marks=marks))
    return out


@pytest.mark.parametrize("name,polarity,additive,expected,min_lk", _params())
def test_reference_set(name, polarity, additive, expected, min_lk):
    _check(name, polarity, additive, expected, min_lk)


@needs_weights
@pytest.mark.parametrize("key", list(MISSES), ids=[k[0] for k in MISSES])
def test_recorded_misses(key):
    name, pol, add = key
    expected, why = MISSES[key]
    try:
        _check(name, pol, add, expected, "low")
    except AssertionError:
        pytest.xfail(why)
    pytest.fail(f"{name} now passes: move it from MISSES to REFERENCE")


def test_reference_set_size_and_coverage():
    names = {r[0] for r in REFERENCE}
    assert len(names) >= 40
    assert {r[1] for r in REFERENCE} == {"Positive", "Negative"} and {FA, AA} <= {r[2] for r in REFERENCE}
    assert set(S) >= names | {k[0] for k in MISSES}


# ---------------------------------------------------------------- rule properties
ACIDS = ["ibuprofen", "salicylic", "palmitic", "cholic", "pfoa", "citric", "taurocholic"]


@needs_weights
@pytest.mark.parametrize("name", ACIDS)
@pytest.mark.parametrize("additive", [FA, AA])
def test_acid_in_positive_mode_is_never_high(name, additive):
    res = _predict(name, "Positive", additive)
    assert all(i.likelihood != "high" for i in res.ions)
    assert _predict(name, "Negative", additive).ions[0].adduct == "[M-H]-"


@needs_weights
def test_base_in_negative_mode_is_only_low():
    res = _predict("amitriptyline", "Negative", FA)
    assert all(i.likelihood == "low" for i in res.ions)


@needs_weights
def test_hydrocarbon_is_all_low_with_a_poor_response_note():
    pos, neg = _predict("pyrene", "Positive", FA), _predict("pyrene", "Negative", FA)
    assert pos.ions and all(i.likelihood == "low" for i in pos.ions)
    assert not neg.ions
    assert any("poor ESI response" in n for n in pos.notes)


@pytest.mark.parametrize("name", ["glucose", "PC", "TG", "Cer"])
def test_no_additive_means_no_formate_acetate_or_ammonium_adducts(name):
    if name not in LIPIDS and NO_MOLGPKA:
        pytest.skip("MolGpKa weights not installed")
    res = _predict(name, "Both", "None")
    assert not any(t in i.adduct for i in res.ions for t in ("HCOO", "CH3COO", "NH4"))


@needs_weights
def test_no_additive_raises_sodium_relative_to_formic_acid():
    assert _predict("glucose", "Positive", "None").ions[0].score > _predict("glucose", "Positive", FA).ions[0].score
    na = {a: next(i.score for i in _predict("caffeine", "Positive", a).ions if i.adduct == "[M+Na]+") for a in ("None", FA, AA)}
    assert na["None"] > na[FA] > na[AA]


def test_additive_picks_the_anion_adduct():
    assert _predict("PC", "Negative", AF).ions[0].adduct == "[M+HCOO]-"
    assert _predict("PC", "Negative", AC).ions[0].adduct == "[M+CH3COO]-"
    assert not any(i.adduct == "[M-H]-" for i in _predict("PC", "Negative", FA).ions)


def test_ammonium_only_with_ammonium_additives():
    assert _predict("TG", "Positive", AF).ions[0].adduct == "[M+NH4]+"
    assert not any(i.adduct == "[M+NH4]+" for i in _predict("TG", "Positive", AC).ions)


def test_alkali_absent_removes_potassium_and_cuts_sodium():
    on, off = _predict("TG", "Positive", FA, True), _predict("TG", "Positive", FA, False)
    assert any(i.adduct == "[M+K]+" for i in on.ions) and not any(i.adduct == "[M+K]+" for i in off.ions)
    na_on = next(i.score for i in on.ions if i.adduct == "[M+Na]+")
    na_off = next((i.score for i in off.ions if i.adduct == "[M+Na]+"), 0)
    assert na_off < na_on * 0.2 and off.ions[0].likelihood == "low"


def test_permanent_cation_has_no_negative_ions_and_no_mh():
    res = _predict("choline", "Both", FA)
    assert [i.adduct for i in res.ions] == ["[M]+"]
    assert res.best["-"] is None


def test_inner_salt_m_is_the_zwitterion():
    assert _predict("PC", "Positive", FA).formula == "C42H82NO8P"
    assert _predict("carnitine", "Positive", FA).formula == "C7H15NO3"
    assert _predict("choline", "Positive", FA).formula == "C5H14NO"


@needs_weights
def test_at_most_one_dimer_per_polarity_and_none_for_lipids():
    for name in ("caffeine", "ibuprofen", "glucose", "phe"):
        res = _predict(name, "Both", FA)
        for pol in "+-":
            assert sum(i.adduct.startswith("[2M") for i in res.ions if i.polarity == pol) <= 1
    assert not any(i.adduct.startswith("[2M") for i in _predict("PC", "Both", FA).ions)


@needs_weights
def test_multiply_charged_ions():
    cl = _predict("CL", "Negative", AA)
    assert next(i for i in cl.ions if i.adduct == "[M-2H]2-").likelihood in ("high", "medium")
    assert not any(i.adduct == "[M+2H]2+" for i in _predict("lys", "Positive", FA).ions)
    gentamicin = "CN[C@@H]1[C@@H](O)[C@@H](O[C@H]2[C@H](N)C[C@H](N)[C@@H](O[C@H]3O[C@H](C(C)NC)CC[C@H]3N)[C@@H]2O)OC[C@]1(C)O"
    res = ar.predict(gentamicin, "Positive", FA)
    z2 = next(i for i in res.ions if i.adduct == "[M+2H]2+")
    assert z2.formula.endswith("2+") and abs(z2.mz - (next(i.mz for i in res.ions if i.adduct == "[M+H]+") + 1.007276) / 2) < 1e-4


@needs_weights
def test_in_source_losses():
    assert any(i.adduct == "[M-H-CO2]-" for i in _predict("ibuprofen", "Negative", FA).ions)
    assert any(i.adduct == "[M+H-NH3]+" for i in _predict("serotonin", "Positive", FA).ions)
    assert any(i.adduct == "[M+H-H2O]+" for i in _predict("cortisol", "Positive", FA).ions)
    assert next(i for i in _predict("Cer", "Positive", FA).ions if i.adduct == "[M+H-H2O]+").likelihood in ("high", "medium")


def test_lipid_subclasses():
    from rdkit import Chem
    names = {k: ar.lipid_subclass(ar.pe.prepare(S[k]).mol)[0] for k in LIPIDS}
    got = {k: (v.name if v else None) for k, v in names.items()}
    assert got == {"PC": "PC", "LPC": "LPC", "PE": "PE", "PS": "PS", "PG": "PG", "PI": "PI", "PA": "PA", "CL": "CL", "SM": "SM",
                   "Cer": "Cer", "HexCer": "HexCer", "TG": "TG", "DG": "DG", "MG": "MG", "MGDG": "MGDG / DGDG",
                   "CE": "CE / steryl ester", "cholesterol": "free sterol"}
    # lipids without class rules fall back to the small-molecule rules
    for smi in (S["estradiol"], S["cholic"], "CCCCCCCCCCCCC/C=C/[C@@H](O)[C@@H](N)CO"):
        assert ar.lipid_subclass(Chem.MolFromSmiles(smi))[0] is None


def test_every_ion_has_a_valid_notation_and_mz():
    used = {ion for rows in ar.RULES.lipid_ions.values() for _, ion, *_ in rows}
    used |= {"[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M]+", "[M+H-H2O]+", "[M+H-NH3]+", "[M+2H]2+", "[2M+H]+", "[2M+Na]+", "[2M+NH4]+",
             "[M-H]-", "[M+HCOO]-", "[M+CH3COO]-", "[M+Cl]-", "[M-H-CO2]-", "[M-H-H2O]-", "[M-2H]2-", "[2M-H]-", "[M-CH3]-", "[M-H]+", "[M]2+"}
    used = {u.replace("[M+A]-", "[M+HCOO]-") for u in used}
    for ion in used:
        mz, formula = ar.ion_mz("C20H40NO8P", ion)
        assert mz and formula, ion
    # notation shared with the rest of the app where it exists
    assert {"[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M]+", "[M+H-H2O]+", "[2M+H]+"} <= set(ADDUCTS_POS)
    assert {"[M-H]-", "[M+Cl]-", "[M+HCOO]-", "[M+CH3COO]-", "[M-H-H2O]-", "[2M-H]-"} <= set(ADDUCTS_NEG)


def test_rules_live_in_one_dataclass():
    r = ar.Rules()
    assert r.high > r.medium > r.min_shown and r.base_pka == 3.0
    assert ar.additive("Formic acid", r).ph == 2.7 and ar.additive("Ammonium acetate", r).ph == 6.8
    with pytest.raises(ValueError):
        ar.additive("TFA", r)


def test_functional_group_fallback_without_molgpka(monkeypatch):
    def broken(mol):
        raise pe.PkaError("weights missing")
    monkeypatch.setattr(pe, "predict_sites", broken)
    for name, pol, add, top in (("caffeine", "Positive", FA, "[M+H]+"), ("ibuprofen", "Negative", FA, "[M-H]-"),
                                ("phe", "Positive", FA, "[M+H]+"), ("glucose", "Negative", FA, "[M+HCOO]-")):
        res = ar.predict(S[name], pol, add)
        assert res.ions[0].adduct == top, name
        assert any("pKa fallback" in n for n in res.notes)


# ---------------------------------------------------------------- page
def _run(inputs):
    return client.post("/api/modules/adducts/run", json=inputs)


def _smiles(s):
    return {"mode": "SMILES / InChI", "value": s}


def _rows(out):
    """[(cells, help)] of the ions table; every row carries its reason as help."""
    return [(r["cells"], r["help"]) for r in out["ions"]["rows"]]


def test_page_spec_blocks():
    mod = next(m for m in MODULES if m["id"] == "adducts")
    blocks = {b["id"]: b for b in mod["outputs"]}
    assert [b["id"] for b in mod["outputs"]] == ["base", "ions", "rule_hits"]
    assert blocks["base"]["style"] == "metric" and [i["label"] for i in blocks["base"]["items"]] == ["Dominant ion", "m/z"]
    assert blocks["ions"]["columns"] == ["Ion", "Formula", "m/z", "Likelihood"] and blocks["ions"]["row_help"] is True and blocks["ions"]["help"]
    assert blocks["rule_hits"]["type"] == "notes" and blocks["rule_hits"]["collapsed"] is True


def test_page_example_and_blocks():
    mod = next(m for m in MODULES if m["id"] == "adducts")
    body = _run(mod["example"]).json()
    assert body["status"] == "ok"
    out = body["outputs"]
    dominant = out["base"]["values"]["Dominant ion"]
    assert dominant["value"] == "[M+HCOO]-" and dominant["note"] == "high likelihood" and "PC" in dominant["help"]       # the reason is the help
    assert out["base"]["values"]["m/z"]["value"] == "804.5760"
    cells, reason = _rows(out)[0]
    assert cells == ["[M+HCOO]-", "C43H83NO10P-", "804.5760", {"badge": "high", "tone": "good"}] and reason == dominant["help"]
    assert out["ions"]["count"] == len(out["ions"]["rows"]) and all(len(c) == 4 for c, _ in _rows(out))
    assert any("PC" in t for t in out["rule_hits"]["items"]) and any("pH ≈ 6.5 assumed" in t for t in out["rule_hits"]["items"])


def test_likelihood_badge_tones():
    out = _run({"structure": _smiles(S["caffeine"]), "ion_mode": "Positive", "additive": FA}).json()["outputs"]
    tones = {c[3]["badge"]: c[3]["tone"] for c, _ in _rows(out)}
    assert tones == {"high": "good", "medium": "warn", "low": "neutral"}


def test_page_accepts_both_and_reports_both_polarities():
    body = _run({"structure": _smiles(S["PE"]), "ion_mode": "Both", "additive": FA}).json()
    assert {c[0][-1] for c, _ in _rows(body["outputs"])} == {"+", "-"}
    dominant = body["outputs"]["base"]["values"]["Dominant ion"]
    assert dominant["value"] == "[M+H]+ / [M-H]-" and [h.split(":")[0] for h in dominant["help"]] == ["positive", "negative"]


def test_page_permanent_cation_negative_says_none():
    body = _run({"structure": _smiles(S["choline"]), "ion_mode": "Negative", "additive": AA}).json()
    dominant = body["outputs"]["base"]["values"]["Dominant ion"]
    assert dominant["value"] == "none expected" and "nothing forms in negative ESI" in dominant["help"]
    assert body["outputs"]["ions"] == {"rows": [], "count": 0}


def test_page_bad_inputs():
    assert _run({"structure": _smiles("CCO"), "ion_mode": "Sideways", "additive": FA}).status_code == 400
    assert _run({"structure": _smiles("CCO"), "ion_mode": "Positive", "additive": "TFA"}).status_code == 400
    assert _run({"structure": _smiles("C1CC"), "ion_mode": "Positive", "additive": FA}).status_code == 400


def test_page_warns_when_input_is_drawn_as_an_ion():
    body = _run({"structure": _smiles("CC(C)Cc1ccc(cc1)C(C)C(=O)[O-].[Na+]"), "ion_mode": "Negative", "additive": AA}).json()
    assert any("drawn as an ion" in w for w in body["warnings"]) and any("disconnected" in w for w in body["warnings"])


# m/z must equal the Theoretical m/z page for the same neutral structure and adduct
@pytest.mark.parametrize("name,polarity,additive,adduct", [
    ("caffeine", "Positive", FA, "[M+H]+"), ("glucose", "Positive", FA, "[M+Na]+"), ("PC", "Negative", FA, "[M+HCOO]-")])
def test_mz_matches_the_theoretical_mz_page(name, polarity, additive, adduct):
    if name not in LIPIDS and NO_MOLGPKA:
        pytest.skip("MolGpKa weights not installed")
    ours = {c[0]: c for c, _ in _rows(_run({"structure": _smiles(S[name]), "ion_mode": polarity, "additive": additive}).json()["outputs"])}
    theirs = {}
    for a in ours:
        r = client.post("/api/modules/ms1_theory/run", json={"structure": _smiles(S[name]),
                                                             "precursor": {"precursor_mz": "", "ion_mode": polarity, "adduct": a}})
        if r.status_code == 200:
            v = r.json()["outputs"]["theory"]["values"]
            theirs[a] = (v["Ion formula"], v["Precursor m/z"]["value"])
    assert ours[adduct][2] == theirs[adduct][1] and ours[adduct][1] == theirs[adduct][0]
    shared = set(ours) & set(theirs)
    assert len(shared) >= 2 and all(ours[a][2] == theirs[a][1] for a in shared)
