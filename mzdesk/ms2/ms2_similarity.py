#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ms2_similarity.py — 计算两张 MS/MS 谱图之间的多种相似度得分。

输入：两张 MS2 谱（JSON 峰列表，precursor m/z 可选）。
输出：终端打印按类别分组的所有相似度得分表。

用法（CLI）：
    python ms2_similarity.py \
        --spec1 '{"peaks": [[100.05, 1000], [200.10, 500]], "precursor_mz": 500.20}' \
        --spec2 spec2.json

    --spec1/--spec2   JSON 字符串、JSON 文件路径、或 @文件路径
    --tolerance       碎片匹配容差（Da），默认 0.02
    --no-deep         跳过深度学习方法（MS2DeepScore / Spec2Vec）
    --ms2deepscore-model PATH   本地 MS2DeepScore 模型路径（默认自动下载并缓存）
    --spec2vec-model PATH       本地 Spec2Vec 模型路径（默认自动下载并缓存）

作为模块导入：
    from ms2_similarity import compute_all_similarities
    results = compute_all_similarities(peaks1, peaks2, precursor_mz1=500.2,
                                       precursor_mz2=516.2, tolerance=0.02)

依赖：
    必需: matchms, ms_entropy, numpy, scipy
    可选: ms2deepscore (+torch), gensim
    注: Spec2Vec 打分直接基于 gensim Word2Vec 模型实现（与官方预训练模型兼容），
        无需安装 spec2vec 包（其与新版 matchms 不兼容）。

覆盖指标（共 14 种）：
    经典峰匹配 (matchms): CosineGreedy, CosineHungarian, CosineLinear,
        ModifiedCosineGreedy, ModifiedCosineHungarian, NeutralLossesCosine, BlinkCosine
    谱熵 (ms_entropy): 加权 / 未加权谱熵相似度
    自定义 (numpy/scipy): NIST 加权点积, Reverse dot product,
        分箱余弦, 分箱 Pearson 相关, 分箱 Spearman 相关
    深度学习: MS2DeepScore, Spec2Vec
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import urllib.request
import uuid
from functools import lru_cache

import numpy as np

from mzdesk.ms2 import metric_names as mn   # 指标名字与类别（与 spec.py 共用）

# --------------------------------------------------------------------------- #
# 模型下载地址（首次使用时自动下载并缓存到 ~/.cache/ms2_similarity/）
# --------------------------------------------------------------------------- #
CACHE_DIR = os.path.join(os.path.expanduser("~"), ".cache", "ms2_similarity")
DOWNLOAD_TIMEOUT_S = 30        # 下载时连接、每次读取的等待上限（秒）；不是总时长上限
DOWNLOAD_CHUNK = 1 << 20

MS2DEEPSCORE_MODEL_URLS = [
    "https://zenodo.org/records/17826815/files/ms2deepscore_model.pt?download=1",
]
# Spec2Vec 预训练模型（GNPS AllPositive, iter 15; Zenodo record 3978070）
# gensim 加载需要主 .model 文件及两个分离保存的 .npy 数组
SPEC2VEC_MODEL_FILES = {
    "spec2vec_AllPositive_ratio05_filtered_iter_15.model":
        "https://zenodo.org/records/3978070/files/spec2vec_AllPositive_ratio05_filtered_iter_15.model?download=1",
    "spec2vec_AllPositive_ratio05_filtered_iter_15.model.wv.vectors.npy":
        "https://zenodo.org/records/3978070/files/spec2vec_AllPositive_ratio05_filtered_iter_15.model.wv.vectors.npy?download=1",
    "spec2vec_AllPositive_ratio05_filtered_iter_15.model.trainables.syn1neg.npy":
        "https://zenodo.org/records/3978070/files/spec2vec_AllPositive_ratio05_filtered_iter_15.model.trainables.syn1neg.npy?download=1",
}


# --------------------------------------------------------------------------- #
# 输入解析与预处理
# --------------------------------------------------------------------------- #
def load_spectrum_input(arg: str) -> dict:
    """接受 JSON 字符串、JSON 文件路径或 @文件路径，返回 dict。"""
    if arg.startswith("@"):
        with open(arg[1:], "r", encoding="utf-8") as fh:
            return json.load(fh)
    if os.path.isfile(arg):
        with open(arg, "r", encoding="utf-8") as fh:
            return json.load(fh)
    return json.loads(arg)


def normalize_peaks_input(peaks) -> tuple[np.ndarray, np.ndarray]:
    """接受 [[mz, intensity], ...] 或 (mz_array, intensity_array)，返回两个数组。"""
    if isinstance(peaks, tuple) and len(peaks) == 2:
        mz = np.asarray(peaks[0], dtype=float)
        intensities = np.asarray(peaks[1], dtype=float)
    else:
        arr = np.asarray(peaks, dtype=float)
        if arr.ndim != 2 or arr.shape[1] != 2:
            raise ValueError("peaks must be in [[mz, intensity], ...] format")
        mz, intensities = arr[:, 0], arr[:, 1]
    return mz, intensities


def clean_peaks(mz: np.ndarray, intensities: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """去除非法/零强度峰，按 m/z 排序，强度归一化到基峰（=100）。"""
    mask = (intensities > 0) & np.isfinite(mz) & np.isfinite(intensities) & (mz > 0)
    mz, intensities = mz[mask], intensities[mask]
    if mz.size == 0:
        raise ValueError("spectrum has no valid peaks")
    order = np.argsort(mz)
    mz, intensities = mz[order], intensities[order]
    intensities = intensities / intensities.max() * 100.0
    return mz, intensities


def to_matchms_spectrum(mz, intensities, precursor_mz=None, ionmode=None):
    """Build a matchms Spectrum (only intensity normalization; no filtering chain needs a precursor)."""
    try:
        from matchms import Spectrum
        from matchms.filtering import normalize_intensities
    except Exception as exc:  # noqa: BLE001 - a broken optional dependency must not break the page
        raise MatchmsUnavailable(f"{type(exc).__name__}: {exc}") from exc

    metadata = {}
    if precursor_mz is not None:
        metadata["precursor_mz"] = float(precursor_mz)
    if ionmode is not None:
        metadata["ionmode"] = ionmode
    spec = Spectrum(mz=np.asarray(mz, dtype=float),
                    intensities=np.asarray(intensities, dtype=float),
                    metadata=metadata)
    return normalize_intensities(spec)


class MatchmsUnavailable(RuntimeError):
    """matchms cannot be imported (the rest of the page can still run)."""


# --------------------------------------------------------------------------- #
# 自定义相似度实现（numpy / scipy）
# --------------------------------------------------------------------------- #
def _greedy_match(mz1, w1, mz2, w2, tolerance):
    """贪心峰匹配：候选峰对按权重乘积降序，每个峰最多匹配一次。"""
    candidates = []
    for i, m in enumerate(mz1):
        lo = np.searchsorted(mz2, m - tolerance)
        hi = np.searchsorted(mz2, m + tolerance, side="right")
        for j in range(lo, hi):
            candidates.append((w1[i] * w2[j], i, j))
    candidates.sort(key=lambda x: -x[0])
    used1, used2, matches = set(), set(), []
    for _, i, j in candidates:
        if i not in used1 and j not in used2:
            used1.add(i)
            used2.add(j)
            matches.append((i, j))
    return matches


def greedy_cosine(mz1, i1, mz2, i2, tolerance):
    """CosineGreedy-compatible score and matched-peak count without matchms.

    Used when matchms cannot be imported (for example a numba/pynndescent cache problem):
    the rest of the Similarity page and the library search still work, and the 7 matchms
    metrics are N/A. ``i1`` / ``i2`` are already normalized to the base peak by clean_peaks().
    """
    w1 = np.asarray(i1, dtype=float) / 100.0
    w2 = np.asarray(i2, dtype=float) / 100.0
    matches = _greedy_match(mz1, w1, mz2, w2, tolerance)
    num = sum(w1[i] * w2[j] for i, j in matches)
    den = float(np.linalg.norm(w1) * np.linalg.norm(w2))
    return (float(num / den) if den > 0 else 0.0), len(matches)


def _nist_weights(mz, intensities, mz_power=3.0, int_power=0.6):
    """NIST 加权：w = mz^mz_power * I^int_power（Stein & Scott 1994）。"""
    return (mz ** mz_power) * (intensities ** int_power)


def nist_dot_product(mz1, i1, mz2, i2, tolerance):
    """Stein & Scott (1994) 复合点积：(Σ W₁W₂)² / (ΣW₁²·ΣW₂²)，W = mz³·I⁰·⁶。

    分子只累加匹配上的峰对，分母用两张谱的全部峰，未匹配的峰因此会拉低得分。
    """
    w1 = _nist_weights(mz1, i1)
    w2 = _nist_weights(mz2, i2)
    matches = _greedy_match(mz1, w1, mz2, w2, tolerance)
    num = sum(w1[i] * w2[j] for i, j in matches)
    den = float(np.dot(w1, w1) * np.dot(w2, w2))
    return (float(num * num / den) if den > 0 else 0.0), len(matches)


def reverse_dot_product(mz1, i1, mz2, i2, tolerance):
    """Reverse dot product：同样的 Stein & Scott 复合点积，但查询谱只用匹配上的峰归一化
    （spec1 = 查询谱，spec2 = 库谱），所以查询谱里多出来的峰不扣分。
    """
    w1 = _nist_weights(mz1, i1)
    w2 = _nist_weights(mz2, i2)
    matches = _greedy_match(mz1, w1, mz2, w2, tolerance)
    num = sum(w1[i] * w2[j] for i, j in matches)
    den = float(sum(w1[i] * w1[i] for i, _ in matches) * np.dot(w2, w2))
    return (float(num * num / den) if den > 0 else 0.0), len(matches)


# 分箱向量长度上限（约 40 MB/向量），防止极小 bin 宽度导致内存爆炸
MAX_BINS = 5_000_000


def _binned_vectors(mz1, i1, mz2, i2, bin_width):
    if not (np.isfinite(bin_width) and bin_width > 0):
        raise ValueError("bin width must be a positive number")
    lo = np.floor(min(mz1.min(), mz2.min()) / bin_width) * bin_width
    hi = max(mz1.max(), mz2.max()) + 2 * bin_width
    n_bins = int(np.ceil((hi - lo) / bin_width))
    if n_bins > MAX_BINS:
        raise ValueError(f"{n_bins:,} bins would be needed at {bin_width:g} Da; "
                         "use a larger tolerance")
    edges = np.arange(lo, hi, bin_width)
    v1, _ = np.histogram(mz1, bins=edges, weights=i1)
    v2, _ = np.histogram(mz2, bins=edges, weights=i2)
    return v1.astype(float), v2.astype(float)


def binned_cosine(mz1, i1, mz2, i2, bin_width):
    v1, v2 = _binned_vectors(mz1, i1, mz2, i2, bin_width)
    den = np.linalg.norm(v1) * np.linalg.norm(v2)
    return float(np.dot(v1, v2) / den) if den > 0 else 0.0


def binned_correlation(mz1, i1, mz2, i2, bin_width, method="pearson"):
    """两张谱分箱后的相关系数，只统计至少一方有信号的 bin。

    整段 m/z 范围里绝大多数 bin 是空的（默认例子有几千个），把它们算进去会让
    相关系数趋近 1（Spearman 因大量并列的 0 甚至正好是 1.0000）。
    """
    from scipy import stats
    v1, v2 = _binned_vectors(mz1, i1, mz2, i2, bin_width)
    used = (v1 > 0) | (v2 > 0)
    v1, v2 = v1[used], v2[used]
    if v1.size < 3:
        return None
    if method == "pearson":
        r = stats.pearsonr(v1, v2).statistic
    else:
        r = stats.spearmanr(v1, v2).statistic
    return float(r) if np.isfinite(r) else None


# --------------------------------------------------------------------------- #
# 各指标计算（统一封装，失败时写入 note，不中断整体运行）
# --------------------------------------------------------------------------- #
def _add_result(results, category, method, score=None, matches=None, note=""):
    results.append({"category": category, "method": method,
                    "score": score, "matches": matches, "note": note})


def _unpack_matchms_score(out):
    """matchms .pair() 返回带 ('score','matches') 字段的 0 维结构化数组。"""
    names = getattr(getattr(out, "dtype", None), "names", None)
    if names:
        score = float(out["score"])
        matches = int(out["matches"]) if "matches" in names else None
        return score, matches
    return float(out), None


def _matchms_scores(spec1, spec2, tolerance, results):
    from matchms import similarity

    cat = mn.CAT_MATCHMS
    has_precursor = (spec1.get("precursor_mz") is not None
                     and spec2.get("precursor_mz") is not None)

    entries = [(name, name, name in mn.NEEDS_PRECURSOR) for name in mn.MATCHMS]   # 名字 = matchms 类名
    for name, cls_name, needs_precursor in entries:
        cls = getattr(similarity, cls_name, None)
        if cls is None:
            _add_result(results, cat, name,
                        note="skipped: not provided by this matchms version")
            continue
        if needs_precursor and not has_precursor:
            _add_result(results, cat, name, note="N/A: missing precursor m/z")
            continue
        try:
            score, matches = _unpack_matchms_score(cls(tolerance=tolerance).pair(spec1, spec2))
            note = ""
            if name == mn.BLINK:
                # BLINK 默认预过滤：去掉 <1% 基峰及高于 precursor 的峰；matches 为窗口内近似计数
                note = ("BLINK approximation: drops peaks <1% of base peak and above the "
                        "precursor; match count is approximate")
            _add_result(results, cat, name, score=score, matches=matches, note=note)
        except Exception as exc:  # noqa: BLE001
            _add_result(results, cat, name, note=f"Failed: {exc}")


def _entropy_scores(mz1, i1, mz2, i2, tolerance, results):
    cat = mn.CAT_ENTROPY
    try:
        import ms_entropy as me
    except ImportError:
        _add_result(results, cat, "SpectralEntropy", note="skipped: ms_entropy is not installed")
        return

    peaks1 = np.column_stack([mz1, i1])
    peaks2 = np.column_stack([mz2, i2])
    for name, fn_name in [(mn.ENTROPY_WEIGHTED, "calculate_entropy_similarity"),
                          (mn.ENTROPY_UNWEIGHTED, "calculate_unweighted_entropy_similarity")]:
        fn = getattr(me, fn_name, None)
        if fn is None:
            _add_result(results, cat, name, note="skipped: ms_entropy does not provide this function")
            continue
        try:
            score = fn(peaks1, peaks2, ms2_tolerance_in_da=tolerance)
            _add_result(results, cat, name, score=float(score))
        except Exception as exc:  # noqa: BLE001
            _add_result(results, cat, name, note=f"Failed: {exc}")


def _custom_scores(mz1, i1, mz2, i2, tolerance, results):
    cat = mn.CAT_CUSTOM
    try:
        score, n = nist_dot_product(mz1, i1, mz2, i2, tolerance)
        _add_result(results, cat, mn.NIST_DOT, score=score, matches=n,
                    note="Stein & Scott (1994) composite dot product")
    except Exception as exc:  # noqa: BLE001
        _add_result(results, cat, mn.NIST_DOT, note=f"Failed: {exc}")
    try:
        score, n = reverse_dot_product(mz1, i1, mz2, i2, tolerance)
        _add_result(results, cat, mn.REVERSE_DOT, score=score, matches=n,
                    note="Directional: query spectrum vs. reference library")
    except Exception as exc:  # noqa: BLE001
        _add_result(results, cat, mn.REVERSE_DOT, note=f"Failed: {exc}")
    for name, fn in [(mn.BINNED_COSINE, lambda: binned_cosine(mz1, i1, mz2, i2, tolerance)),
                     (mn.BINNED_PEARSON, lambda: binned_correlation(mz1, i1, mz2, i2, tolerance, "pearson")),
                     (mn.BINNED_SPEARMAN, lambda: binned_correlation(mz1, i1, mz2, i2, tolerance, "spearman"))]:
        try:
            score = fn()
            note = f"bin width = {tolerance:g} Da"
            if score is None:
                note += "; correlation undefined for these binned vectors"
            _add_result(results, cat, name, score=score, note=note)
        except Exception as exc:  # noqa: BLE001
            _add_result(results, cat, name, note=f"Failed: {exc}")


def _download_file(url, dest):
    """下载单个文件，成功返回 True。

    连接和每次读取最多等 DOWNLOAD_TIMEOUT_S 秒（服务器停住时放弃，返回 False）；先写到同目录下本次独有的临时文件，
    大小对得上（有 Content-Length 时）才 os.replace 成 dest，所以两个同时开始的下载不会互相覆盖，失败也不留文件。
    """
    tmp = None
    try:
        print(f"  Downloading: {url}", file=sys.stderr)
        with urllib.request.urlopen(url, timeout=DOWNLOAD_TIMEOUT_S) as resp:
            expected = resp.headers.get("Content-Length")
            tmp = f"{dest}.{uuid.uuid4().hex[:12]}.part"
            # 0o666 & ~umask (like any file the process creates): mkstemp would make the cached model 0600
            fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_BINARY", 0), 0o666)
            size = 0
            with os.fdopen(fd, "wb") as out:
                while True:
                    chunk = resp.read(DOWNLOAD_CHUNK)
                    if not chunk:
                        break
                    out.write(chunk)
                    size += len(chunk)
        if size == 0:
            raise OSError("the server sent an empty file")
        if expected is not None and expected.strip().isdigit() and size != int(expected):
            raise OSError(f"incomplete download: got {size} of {int(expected)} bytes")
        os.replace(tmp, dest)   # a partial download must never become the cached model
        return True
    except Exception as exc:  # noqa: BLE001
        print(f"  Download failed: {exc}", file=sys.stderr)
        return False
    finally:
        if tmp and os.path.exists(tmp):
            os.remove(tmp)


def _get_ms2deepscore_model():
    os.makedirs(CACHE_DIR, exist_ok=True)
    dest = os.path.join(CACHE_DIR, "ms2deepscore_model.pt")
    if os.path.isfile(dest) and os.path.getsize(dest) > 0:
        return dest
    for url in MS2DEEPSCORE_MODEL_URLS:
        if _download_file(url, dest):
            return dest
    return None


def _get_spec2vec_model():
    os.makedirs(CACHE_DIR, exist_ok=True)
    main_dest = os.path.join(CACHE_DIR, "spec2vec_AllPositive_ratio05_filtered_iter_15.model")
    for fname, url in SPEC2VEC_MODEL_FILES.items():
        dest = os.path.join(CACHE_DIR, fname)
        if os.path.isfile(dest) and os.path.getsize(dest) > 0:
            continue
        if not _download_file(url, dest):
            return None
    return main_dest


@lru_cache(maxsize=1)
def _load_ms2deepscore_model(model_path: str):
    """Load and cache the MS2DeepScore model in memory for the process lifetime."""
    from ms2deepscore.models import load_model
    return load_model(model_path)


def _ms2deepscore_score(spec1, spec2, model_path, results):
    cat, name = mn.CAT_DEEP, mn.MS2DEEPSCORE
    if spec1.get("precursor_mz") is None or spec2.get("precursor_mz") is None:
        _add_result(results, cat, name, note="N/A: missing precursor m/z")
        return
    try:
        from ms2deepscore import MS2DeepScore
    except ImportError:
        _add_result(results, cat, name, note="skipped: ms2deepscore is not installed")
        return
    if model_path is None:
        model_path = _get_ms2deepscore_model()
    if model_path is None or not os.path.isfile(model_path):
        _add_result(results, cat, name,
                    note="skipped: model unavailable; provide a local path via --ms2deepscore-model")
        return
    try:
        model = _load_ms2deepscore_model(model_path)
        sim = MS2DeepScore(model, progress_bar=False)
        score, _ = _unpack_matchms_score(sim.pair(spec1, spec2))
        _add_result(results, cat, name, score=score, note="Predicted structural Tanimoto similarity")
    except Exception as exc:  # noqa: BLE001
        _add_result(results, cat, name, note=f"Failed: {exc}")


class _DummyBase:
    """占位基类：实例可调用、可迭代、属性可任意访问。"""

    def __call__(self, *args, **kwargs):
        return None

    def __iter__(self):
        return iter(())

    def __getattr__(self, attr):
        return _DummyBase()


def _install_custom_functions_shim():
    """官方 Spec2Vec 预训练模型（gensim 3.x 保存）的 pickle 引用了训练仓库中的
    custom_functions 模块（及其子模块）。注册一个 meta path finder，为
    custom_functions[.*] 提供占位模块，使模型可在新版 gensim 中加载。"""
    import importlib.abc
    import importlib.machinery
    import types

    def _make_dummy_module(name):
        mod = types.ModuleType(name)
        mod.__path__ = []

        def _getattr(attr):
            if attr.startswith("__") and attr.endswith("__"):
                raise AttributeError(attr)
            # 返回动态类而非实例：同时兼容函数引用与 pickle NEWOBJ 类实例化
            return type(attr, (_DummyBase,), {})

        mod.__getattr__ = _getattr
        return mod

    class _DummyLoader(importlib.abc.Loader):
        def create_module(self, spec):
            return _make_dummy_module(spec.name)

        def exec_module(self, module):
            pass

    class _CustomFunctionsFinder(importlib.abc.MetaPathFinder):
        def find_spec(self, fullname, path=None, target=None):
            if fullname == "custom_functions" or fullname.startswith("custom_functions."):
                return importlib.machinery.ModuleSpec(fullname, _DummyLoader(), is_package=True)
            return None

    if not any(type(f).__name__ == "_CustomFunctionsFinder" for f in sys.meta_path):
        sys.meta_path.insert(0, _CustomFunctionsFinder())


def _spectrum_to_words(mz, intensities, precursor_mz, n_decimals):
    """把谱图转成 spec2vec 词序列（peak@xxx / loss@xxx）及强度权重。"""
    words = [f"peak@{m:.{n_decimals}f}" for m in mz]
    weights = list(intensities ** 0.5)
    if precursor_mz is not None:
        for m, inten in zip(mz, intensities):
            loss = float(precursor_mz) - m
            if loss > 0:
                words.append(f"loss@{loss:.{n_decimals}f}")
                weights.append(inten ** 0.5)
    return words, weights


@lru_cache(maxsize=1)
def _load_spec2vec_model(model_path: str):
    """Load and cache the Spec2Vec gensim model in memory for the process lifetime."""
    from gensim.models import Word2Vec
    _install_custom_functions_shim()
    return Word2Vec.load(model_path)


def spec2vec_decimals(wv, mz, intensities):
    """根据（查询谱前 20 个峰的）词表覆盖率自动选择 n_decimals（1 或 2）。"""
    best_nd, best_cov = 2, -1.0
    for nd in (2, 1):
        words, _ = _spectrum_to_words(mz[:20], intensities[:20], None, nd)
        cov = float(np.mean([w in wv for w in words])) if words else 0.0
        if cov > best_cov:
            best_nd, best_cov = nd, cov
    return best_nd


def spec2vec_embedding(wv, mz, intensities, precursor_mz, n_decimals):
    """返回 (嵌入向量, 峰词权重的词表覆盖率)；词表外的词被忽略，一个词都没命中返回 (None, 0.0)。"""
    words, weights = _spectrum_to_words(mz, intensities, precursor_mz, n_decimals)
    vecs, ws = [], []
    peak_total = peak_found = 0.0
    for w, wt in zip(words, weights):
        is_peak = w.startswith("peak@")
        if is_peak:
            peak_total += wt
        if w in wv:
            vecs.append(wv[w] * wt)
            ws.append(wt)
            if is_peak:
                peak_found += wt
    if not vecs:
        return None, 0.0
    coverage = peak_found / peak_total if peak_total > 0 else 0.0
    return np.sum(vecs, axis=0) / sum(ws), coverage


def _spec2vec_score(mz1, i1, pmz1, mz2, i2, pmz2, model_path, results, ionmode=None):
    """基于 gensim Word2Vec 模型直接计算 Spec2Vec 嵌入余弦（无需 spec2vec 包）。"""
    cat, name = mn.CAT_DEEP, mn.SPEC2VEC
    if ionmode == "negative":
        _add_result(results, cat, name,
                    note="N/A: the pretrained model covers positive-mode spectra only")
        return
    try:
        import gensim  # noqa: F401
    except ImportError:
        _add_result(results, cat, name, note="skipped: gensim is not installed")
        return
    if model_path is None:
        model_path = _get_spec2vec_model()
    if model_path is None or not os.path.isfile(model_path):
        _add_result(results, cat, name,
                    note="skipped: model unavailable; provide a local path via --spec2vec-model")
        return
    try:
        model = _load_spec2vec_model(model_path)
        wv = model.wv

        best_nd = spec2vec_decimals(wv, mz1, i1)

        def embedding(mz, ii, pmz):
            return spec2vec_embedding(wv, mz, ii, pmz, best_nd)

        # neutral-loss words only when BOTH spectra have a precursor m/z: adding them on one
        # side only compares different word sets and depresses the score
        use_losses = pmz1 is not None and pmz2 is not None
        v1, cov1 = embedding(mz1, i1, pmz1 if use_losses else None)
        v2, cov2 = embedding(mz2, i2, pmz2 if use_losses else None)
        if v1 is None or v2 is None:
            _add_result(results, cat, name, note="Failed: spectrum words are not in the model vocabulary")
            return
        den = np.linalg.norm(v1) * np.linalg.norm(v2)
        score = float(np.dot(v1, v2) / den) if den > 0 else 0.0
        note = (f"GNPS AllPositive pretrained model, n_decimals={best_nd}; vocabulary covers "
                f"{cov1:.0%} / {cov2:.0%} of query / ref peak weight")
        if (pmz1 is None) != (pmz2 is None):
            note += ("; only one spectrum has a precursor m/z, so neutral-loss words are "
                     "left out on both sides")
        _add_result(results, cat, name, score=score, note=note)
    except Exception as exc:  # noqa: BLE001
        _add_result(results, cat, name, note=f"Failed: {exc}")


# --------------------------------------------------------------------------- #
# 主入口
# --------------------------------------------------------------------------- #
def compute_all_similarities(peaks1, peaks2, precursor_mz1=None, precursor_mz2=None,
                             tolerance=0.02, use_deep=True, ionmode="positive",
                             ms2deepscore_model_path=None, spec2vec_model_path=None):
    """
    计算两张 MS2 谱之间的所有相似度得分。

    参数
    ----
    peaks1, peaks2 : [[mz, intensity], ...] 或 (mz_array, intensity_array)
    precursor_mz1, precursor_mz2 : float 或 None（可选）
    tolerance : 碎片匹配容差（Da）
    use_deep : 是否计算深度学习方法（MS2DeepScore / Spec2Vec）
    ionmode : "positive" 或 "negative"（MS2DeepScore 需要），默认 "positive"

    返回
    ----
    list[dict]，字段：category, method, score, matches, note
    """
    tolerance = float(tolerance)
    if not (np.isfinite(tolerance) and tolerance > 0):
        raise ValueError("tolerance must be a positive number (Da)")
    precursor_mz1 = _check_precursor(precursor_mz1, "precursor m/z of spectrum 1")
    precursor_mz2 = _check_precursor(precursor_mz2, "precursor m/z of spectrum 2")

    mz1, i1 = clean_peaks(*normalize_peaks_input(peaks1))
    mz2, i2 = clean_peaks(*normalize_peaks_input(peaks2))

    results: list[dict] = []

    matchms_error = None
    try:
        spec1 = to_matchms_spectrum(mz1, i1, precursor_mz1, ionmode)
        spec2 = to_matchms_spectrum(mz2, i2, precursor_mz2, ionmode)
    except MatchmsUnavailable as exc:
        matchms_error = f"matchms is unavailable ({exc})"
        spec1 = spec2 = None
        for name in mn.MATCHMS:
            _add_result(results, mn.CAT_MATCHMS, name, note=f"N/A: {matchms_error}")

    if matchms_error is None:
        _matchms_scores(spec1, spec2, tolerance, results)
    _entropy_scores(mz1, i1, mz2, i2, tolerance, results)
    _custom_scores(mz1, i1, mz2, i2, tolerance, results)

    if use_deep:
        if matchms_error is None:
            _ms2deepscore_score(spec1, spec2, ms2deepscore_model_path, results)
        else:
            _add_result(results, mn.CAT_DEEP, mn.MS2DEEPSCORE, note=f"N/A: {matchms_error}")
        _spec2vec_score(mz1, i1, precursor_mz1, mz2, i2, precursor_mz2,
                        spec2vec_model_path, results, ionmode=ionmode)

    return results


def _check_precursor(value, label):
    """None 保持 None；否则必须是有限正数（NaN/负数会让 matchms 静默给出 0 分）。"""
    if value is None:
        return None
    value = float(value)
    if not (np.isfinite(value) and value > 0):
        raise ValueError(f"{label} must be a positive number")
    return value


def print_results(results, info_lines=None):
    """终端打印按类别分组的得分表。"""
    if info_lines:
        for line in info_lines:
            print(line)
        print()

    headers = ["Method", "Score", "Matched peaks", "Note"]
    rows = []
    for r in results:
        rows.append([
            r["method"],
            f"{r['score']:.4f}" if r["score"] is not None else "N/A",
            str(r["matches"]) if r["matches"] is not None else "-",
            r["note"],
        ])
    widths = [len(h) for h in headers]
    for row in rows:
        for k, cell in enumerate(row):
            widths[k] = max(widths[k], len(str(cell)))

    def fmt(row):
        return "  ".join(str(c).ljust(widths[k]) for k, c in enumerate(row))

    current_cat = None
    for r, row in zip(results, rows):
        if r["category"] != current_cat:
            current_cat = r["category"]
            print(f"\n-- {current_cat} " + "-" * max(2, 60 - len(current_cat)))
            print(fmt(headers))
            print("  ".join("-" * w for w in widths))
        print(fmt(row))
    print()


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Compute multiple similarity scores between two MS2 spectra",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example:\n"
               '  python ms2_similarity.py --spec1 \'{"peaks": [[100.05, 1000]], "precursor_mz": 500.2}\''
               ' --spec2 spec2.json\n')
    parser.add_argument("--spec1", required=True, help="Spectrum 1: JSON string / JSON file path / @file path")
    parser.add_argument("--spec2", required=True, help="Spectrum 2: same as above")
    parser.add_argument("--tolerance", type=float, default=0.02, help="Fragment matching tolerance (Da), default 0.02")
    parser.add_argument("--no-deep", action="store_true", help="Skip deep-learning methods")
    parser.add_argument("--ionmode", default="positive", choices=["positive", "negative"],
                        help="Ion mode (required by MS2DeepScore), default positive")
    parser.add_argument("--ms2deepscore-model", default=None, help="Local MS2DeepScore model path")
    parser.add_argument("--spec2vec-model", default=None, help="Local Spec2Vec model path")
    args = parser.parse_args(argv)

    s1 = load_spectrum_input(args.spec1)
    s2 = load_spectrum_input(args.spec2)

    results = compute_all_similarities(
        s1["peaks"], s2["peaks"],
        precursor_mz1=s1.get("precursor_mz"),
        precursor_mz2=s2.get("precursor_mz"),
        tolerance=args.tolerance,
        use_deep=not args.no_deep,
        ionmode=args.ionmode,
        ms2deepscore_model_path=args.ms2deepscore_model,
        spec2vec_model_path=args.spec2vec_model,
    )

    mz1, i1 = clean_peaks(*normalize_peaks_input(s1["peaks"]))
    mz2, i2 = clean_peaks(*normalize_peaks_input(s2["peaks"]))
    info = [
        f"Spectrum 1: {len(mz1)} peaks, precursor_mz={s1.get('precursor_mz', 'N/A')}",
        f"Spectrum 2: {len(mz2)} peaks, precursor_mz={s2.get('precursor_mz', 'N/A')}",
        f"Fragment matching tolerance: {args.tolerance} Da",
    ]
    print_results(results, info)
    return 0


if __name__ == "__main__":
    sys.exit(main())
