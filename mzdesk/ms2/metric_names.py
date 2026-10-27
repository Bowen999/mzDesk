"""相似度指标的名字与类别（只用标准库）。

Similarity 页的得分表、Library search 页的指标多选、ms2_similarity.py 的计算都从这里取名字，
spec.py 也导入它（所以这里不能引入 numpy / matchms），名字因此不会漂移。顺序 = 得分表的顺序。
"""
CAT_MATCHMS = "Classic peak matching (matchms)"
CAT_ENTROPY = "Spectral entropy (ms_entropy)"
CAT_CUSTOM = "Custom (numpy/scipy)"
CAT_DEEP = "Deep learning"

COSINE_GREEDY = "CosineGreedy"
COSINE_HUNGARIAN = "CosineHungarian"
COSINE_LINEAR = "CosineLinear"
MOD_COSINE_GREEDY = "ModifiedCosineGreedy"
MOD_COSINE_HUNGARIAN = "ModifiedCosineHungarian"
NEUTRAL_LOSSES = "NeutralLossesCosine"
BLINK = "BlinkCosine"
ENTROPY_WEIGHTED = "SpectralEntropy (weighted)"
ENTROPY_UNWEIGHTED = "SpectralEntropy (unweighted)"
NIST_DOT = "NIST weighted dot product (mz^3*I^0.6)"
REVERSE_DOT = "Reverse dot product"
BINNED_COSINE = "Binned cosine"
BINNED_PEARSON = "Binned Pearson correlation"
BINNED_SPEARMAN = "Binned Spearman correlation"
MS2DEEPSCORE = "MS2DeepScore"
SPEC2VEC = "Spec2Vec"

MATCHMS = (COSINE_GREEDY, COSINE_HUNGARIAN, COSINE_LINEAR, MOD_COSINE_GREEDY, MOD_COSINE_HUNGARIAN, NEUTRAL_LOSSES, BLINK)
ENTROPY = (ENTROPY_WEIGHTED, ENTROPY_UNWEIGHTED)
CUSTOM = (NIST_DOT, REVERSE_DOT, BINNED_COSINE, BINNED_PEARSON, BINNED_SPEARMAN)
DEEP = (MS2DEEPSCORE, SPEC2VEC)

CLASSIC = list(MATCHMS + ENTROPY + CUSTOM)          # 14
ALL = CLASSIC + list(DEEP)                          # 16，Similarity 页得分表的顺序
# Library search 页指标多选的选项顺序：谱熵排第一（Similarity 页的主分数），所以默认按谱熵排序；其余项顺序不变。排序用"所选指标里最靠前的一项"
LIBRARY_ORDER = list(ENTROPY + MATCHMS + CUSTOM + DEEP)
assert sorted(LIBRARY_ORDER) == sorted(ALL)
CATEGORY = {**{n: CAT_MATCHMS for n in MATCHMS}, **{n: CAT_ENTROPY for n in ENTROPY},
            **{n: CAT_CUSTOM for n in CUSTOM}, **{n: CAT_DEEP for n in DEEP}}

# 需要两张谱都有 precursor m/z 的指标（缺了就是 N/A）
NEEDS_PRECURSOR = (MOD_COSINE_GREEDY, MOD_COSINE_HUNGARIAN, NEUTRAL_LOSSES, MS2DEEPSCORE)
