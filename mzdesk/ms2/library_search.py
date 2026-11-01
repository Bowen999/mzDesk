"""谱库检索引擎：一张查询谱 vs 本地 MassBank / GNPS 索引，只算用户选的相似度指标。

两种谱：
  MS/MS  前体 m/z 窗口（ppm）+ 同离子模式 +（可选）同加合物 -> 候选；每个候选算全部所选指标（精确，不做预筛）。
         加合物按规范写法比（mzdesk/adducts.py）：[M+FA-H]- / [M+HCOO]- / GNPS 的 [M+FA]- 都是 [M+HCOO]-。
         同时选了 MassBank 和 GNPS 时，GNPS 里转存的 MassBank 记录（副本）只要原件也在候选里就不算（只留 MassBank 原件）。
         候选太多时按前体偏差保留最近的 MAX_MSMS_CANDIDATES 个并警告。
  EI     前体不用，碎片容差固定 0.5 Da。第一步：对库里全部 EI 谱算加权谱熵相似度（ms_entropy 原函数，逐对，约 10 µs/对），
         留前 EI_PREFILTER_N 个；第二步：只对这些算所选指标。库里 EI 谱不超过 EI_PREFILTER_N 条时第二步覆盖全部（精确）。
保留规则：任一所选指标 >= 阈值（N/A 忽略）。排序：按所选指标里排在最前（顺序 = metric_names.LIBRARY_ORDER）且对命中有分数的那个降序
（某个指标对所有命中都是 N/A 就跳过它），并列比匹配峰数（贪心 1:1，同 CosineGreedy），再比前体偏差。

打分复用 ms2_similarity.py 的实现：matchms 各类、ms_entropy、NIST / Reverse 点积、Spec2Vec 嵌入；
省掉的是逐对开销（不建 matchms Spectrum 对象；三个分箱指标共用一次分箱，Pearson / Spearman 用 rankdata + 相关系数）。
深度指标整批算：模型只载入一次，嵌入按批算。服务器没装 torch / gensim 时这两项返回 N/A 加说明，不会报错。
"""
from __future__ import annotations

import math
import os

import numpy as np

from mzdesk import adducts
from mzdesk.ms2 import library_index as lx
from mzdesk.ms2 import metric_names as mn
from mzdesk.ms2 import ms2_similarity as sim
from mzdesk.ms2.library_parse import adduct_charge_sign

SPECTRUM_TYPES = ("MS/MS", "EI")
EI_TOLERANCE = 0.5                 # EI 按单位质量匹配
EI_PREFILTER_N = 5000              # EI 第一步（谱熵）后留给第二步的候选数
MAX_MSMS_CANDIDATES = 2000         # MS/MS 前体窗口内最多打分的候选数
EI_NA = {mn.MOD_COSINE_GREEDY: "needs a precursor m/z, which EI spectra do not use",
         mn.MOD_COSINE_HUNGARIAN: "needs a precursor m/z, which EI spectra do not use",
         mn.NEUTRAL_LOSSES: "needs a precursor m/z, which EI spectra do not use",
         mn.MS2DEEPSCORE: "trained on ESI MS/MS spectra, not EI",
         mn.SPEC2VEC: "trained on ESI MS/MS spectra, not EI"}


def order_metrics(selected):
    """所选指标按选项列表（metric_names.LIBRARY_ORDER）的顺序排好（不是点选顺序）；未知名字抛 ValueError。"""
    names = list(dict.fromkeys(selected))
    unknown = [n for n in names if n not in mn.LIBRARY_ORDER]
    if unknown:
        raise ValueError(f"Unknown metric: {', '.join(unknown)}.")
    if not names:
        raise ValueError("Select at least one metric.")
    return [m for m in mn.LIBRARY_ORDER if m in names]


# ---------------------------------------------------------------- 逐对打分
class _Fragments:
    """matchms Fragments 的最小替身：只给相似度类用到的 mz / intensities / to_numpy。"""
    __slots__ = ("mz", "intensities")

    def __init__(self, mz, intensities):
        self.mz, self.intensities = mz, intensities

    @property
    def to_numpy(self):
        return np.vstack((self.mz, self.intensities)).T          # 与 matchms 同样的形状 / 内存布局

    def __len__(self):
        return self.mz.size


class _Spectrum:
    """matchms Spectrum 的最小替身（建真 Spectrum 约 0.2-0.4 ms/条，比大部分指标本身还慢）。"""
    __slots__ = ("peaks", "_precursor")

    def __init__(self, mz, intensities, precursor_mz):
        self.peaks = _Fragments(mz, intensities / intensities.max())      # 与 to_matchms_spectrum 一样归一到基峰 = 1
        self._precursor = precursor_mz

    def get(self, key, default=None):
        if key == "precursor_mz" and self._precursor is not None:
            return float(self._precursor)
        return default


def _binned_scores(mz1, i1, mz2, i2, bin_width, want):
    """Binned cosine / Pearson / Spearman，共用一次分箱（ms2_similarity._binned_vectors）；数值与 Similarity 页的单指标函数一致。

    单指标函数每个都重新分箱一遍，Pearson / Spearman 还走 scipy.stats（约 1-2.5 ms/对）；这里只分一次，相关系数用 rankdata + numpy。
    """
    v1, v2 = sim._binned_vectors(mz1, i1, mz2, i2, bin_width)
    out = {}
    if mn.BINNED_COSINE in want:
        den = np.linalg.norm(v1) * np.linalg.norm(v2)
        out[mn.BINNED_COSINE] = float(np.dot(v1, v2) / den) if den > 0 else 0.0
    if mn.BINNED_PEARSON in want or mn.BINNED_SPEARMAN in want:
        used = (v1 > 0) | (v2 > 0)                      # 只统计至少一方有信号的箱（见 ms2_similarity.binned_correlation）
        x, y = v1[used], v2[used]
        if x.size < 3:
            out.update({n: None for n in (mn.BINNED_PEARSON, mn.BINNED_SPEARMAN) if n in want})
        else:
            if mn.BINNED_PEARSON in want:
                out[mn.BINNED_PEARSON] = _pearson(x, y)
            if mn.BINNED_SPEARMAN in want:
                from scipy.stats import rankdata
                out[mn.BINNED_SPEARMAN] = _pearson(rankdata(x), rankdata(y))
    return out


def _pearson(x, y):
    xm, ym = x - x.mean(), y - y.mean()
    den = np.sqrt(np.dot(xm, xm) * np.dot(ym, ym))
    if not den > 0:
        return None
    r = float(np.dot(xm, ym) / den)
    return max(-1.0, min(1.0, r))


class _Scorer:
    """一张查询谱对多个候选的逐对打分器（候选的峰已是 clean_peaks 同样的形式：m/z 升序、基峰 = 100）。"""

    def __init__(self, q_mz, q_int, q_precursor, tol, metrics, ei):
        self.q_mz, self.q_int, self.q_prec, self.tol, self.ei = q_mz, q_int, q_precursor, tol, ei
        self.metrics = metrics
        self.pair_metrics = [m for m in metrics if m not in mn.DEEP]
        self.matchms = {}
        self.matchms_error = None
        try:
            from matchms import similarity          # 贪心匹配峰数总要算（平局排序、"Matched peaks" 列），所以 CosineGreedy 总建
        except Exception as exc:  # noqa: BLE001 - fall back to the local greedy cosine
            self.matchms_error = f"{type(exc).__name__}: {exc}"
        else:
            for name in dict.fromkeys([mn.COSINE_GREEDY] + [m for m in self.pair_metrics if m in mn.MATCHMS]):
                cls = getattr(similarity, name, None)
                if cls is not None:
                    self.matchms[name] = cls(tolerance=tol)
        self.q_spec = _Spectrum(q_mz, q_int, q_precursor)
        self.q_peaks = np.column_stack([q_mz, q_int])
        self.binned = {m for m in self.pair_metrics if m in (mn.BINNED_COSINE, mn.BINNED_PEARSON, mn.BINNED_SPEARMAN)}
        self.me = None
        if any(m in mn.ENTROPY for m in self.pair_metrics):
            import ms_entropy
            self.me = ms_entropy

    def na_reason(self, metric):
        if self.matchms_error and metric in mn.MATCHMS and metric != mn.COSINE_GREEDY:
            return f"N/A: matchms is unavailable ({self.matchms_error})"
        if self.ei and metric in EI_NA:
            return EI_NA[metric]
        return None

    def score(self, mz2, i2, p2):
        """-> ({指标: 分数或 None}, 贪心匹配峰数)。"""
        out = {}
        cand = _Spectrum(mz2, i2, p2)
        has_prec = self.q_prec is not None and p2 is not None
        matches = 0
        g = self.matchms.get(mn.COSINE_GREEDY)
        if g is not None:
            try:
                s, matches = sim._unpack_matchms_score(g.pair(self.q_spec, cand))
            except Exception:  # noqa: BLE001
                s, matches = None, 0
            if mn.COSINE_GREEDY in self.pair_metrics:
                out[mn.COSINE_GREEDY] = s
        else:
            # matchms cannot be imported: use the same greedy 1:1 pairing and cosine formula.
            s, matches = sim.greedy_cosine(self.q_mz, self.q_int, mz2, i2, self.tol)
            if mn.COSINE_GREEDY in self.pair_metrics:
                out[mn.COSINE_GREEDY] = s
        for name in self.pair_metrics:
            if name in out:
                continue
            if name in mn.MATCHMS:
                inst = self.matchms.get(name)
                if inst is None or (name in mn.NEEDS_PRECURSOR and not has_prec):
                    out[name] = None
                    continue
                try:
                    out[name] = sim._unpack_matchms_score(inst.pair(self.q_spec, cand))[0]
                except Exception:  # noqa: BLE001
                    out[name] = None
            elif name in mn.ENTROPY:
                try:
                    fn = self.me.calculate_entropy_similarity if name == mn.ENTROPY_WEIGHTED else self.me.calculate_unweighted_entropy_similarity
                    out[name] = float(fn(self.q_peaks, np.column_stack([mz2, i2]), ms2_tolerance_in_da=self.tol))
                except Exception:  # noqa: BLE001
                    out[name] = None
            elif name == mn.NIST_DOT:
                out[name] = self._safe(lambda: sim.nist_dot_product(self.q_mz, self.q_int, mz2, i2, self.tol)[0])
            elif name == mn.REVERSE_DOT:
                out[name] = self._safe(lambda: sim.reverse_dot_product(self.q_mz, self.q_int, mz2, i2, self.tol)[0])
        if self.binned:
            try:
                out.update(_binned_scores(self.q_mz, self.q_int, mz2, i2, self.tol, self.binned))
            except Exception:  # noqa: BLE001
                out.update({n: None for n in self.binned})
        return out, matches

    @staticmethod
    def _safe(fn):
        try:
            return float(fn())
        except Exception:  # noqa: BLE001
            return None


# ---------------------------------------------------------------- 深度指标（整批）
def _deep_scores(metrics, q_mz, q_int, q_prec, positive, cands, tol, notes):
    """{指标: 分数数组(NaN = N/A)}。cands = [(mz, 强度(基峰 100), 前体 m/z)]。失败写 notes，不抛异常。"""
    out = {}
    n = len(cands)
    for name in metrics:
        out[name] = np.full(n, np.nan)
    if not n:
        return out
    ionmode = "positive" if positive else "negative"
    if mn.MS2DEEPSCORE in metrics:
        _ms2deepscore_batch(out[mn.MS2DEEPSCORE], q_mz, q_int, q_prec, ionmode, cands, notes)
    if mn.SPEC2VEC in metrics:
        _spec2vec_batch(out[mn.SPEC2VEC], q_mz, q_int, q_prec, ionmode, cands, notes)
    return out


def _ms2deepscore_batch(dest, q_mz, q_int, q_prec, ionmode, cands, notes):
    name = mn.MS2DEEPSCORE
    if q_prec is None:
        notes[name] = "N/A: the query has no precursor m/z"
        return
    try:
        from ms2deepscore import MS2DeepScore
        from ms2deepscore.vector_operations import cosine_similarity_matrix
    except ImportError:
        notes[name] = "N/A: ms2deepscore (PyTorch) is not installed on this server"
        return
    try:
        path = sim._get_ms2deepscore_model()
        if path is None:
            notes[name] = "N/A: the MS2DeepScore model is unavailable (download failed)"
            return
        scorer = MS2DeepScore(sim._load_ms2deepscore_model(path), progress_bar=False)
        specs = [sim.to_matchms_spectrum(q_mz, q_int, q_prec, ionmode)]
        specs += [sim.to_matchms_spectrum(mz, ii, pm, ionmode) for mz, ii, pm in cands]
        emb = scorer.get_embedding_array(specs)
        dest[:] = cosine_similarity_matrix(emb[:1], emb[1:])[0]
    except Exception as exc:  # noqa: BLE001
        notes[name] = f"N/A: {type(exc).__name__}: {exc}"


def _spec2vec_batch(dest, q_mz, q_int, q_prec, ionmode, cands, notes):
    name = mn.SPEC2VEC
    if ionmode == "negative":
        notes[name] = "N/A: the pretrained model covers positive-mode spectra only"
        return
    try:
        import gensim  # noqa: F401
    except ImportError:
        notes[name] = "N/A: gensim is not installed on this server"
        return
    try:
        path = sim._get_spec2vec_model()
        if path is None or not os.path.isfile(path):
            notes[name] = "N/A: the Spec2Vec model is unavailable (download failed)"
            return
        wv = sim._load_spec2vec_model(path).wv
        nd = sim.spec2vec_decimals(wv, q_mz, q_int)
        v1_loss, _ = sim.spec2vec_embedding(wv, q_mz, q_int, q_prec, nd)       # 候选都有前体时带 neutral loss 词
        v1_plain, _ = sim.spec2vec_embedding(wv, q_mz, q_int, None, nd)
        for k, (mz, ii, pm) in enumerate(cands):
            use_losses = q_prec is not None and pm is not None
            v1 = v1_loss if use_losses else v1_plain
            v2, _ = sim.spec2vec_embedding(wv, mz, ii, pm if use_losses else None, nd)
            if v1 is None or v2 is None:
                continue
            den = np.linalg.norm(v1) * np.linalg.norm(v2)
            dest[k] = float(np.dot(v1, v2) / den) if den > 0 else 0.0
    except Exception as exc:  # noqa: BLE001
        notes[name] = f"N/A: {type(exc).__name__}: {exc}"


# ---------------------------------------------------------------- 检索
def _cand_arrays(lib, ids):
    """行号 -> [(mz, 强度(基峰 100), 前体 m/z 或 None)]，顺序同 ids。"""
    mz_flat, rel_flat, lens = lib.gather(ids)
    off = np.concatenate([[0], np.cumsum(lens)])
    inten = rel_flat.astype(np.float64) * 100.0
    prec = np.asarray(lib.prec_mz[np.asarray(ids, dtype=np.int64)])
    return [(mz_flat[off[k]:off[k + 1]], inten[off[k]:off[k + 1]], float(prec[k]) if prec[k] > 0 else None) for k in range(len(ids))]


def _entropy_step(q_peaks, lib, tol):
    """EI 第一步：库里全部 EI 谱的加权谱熵相似度（ms_entropy 原函数，与 SpectralEntropy (weighted) 指标一致）。"""
    import ms_entropy
    ids = lib.ei_ids()
    if not ids.size:
        return ids, np.zeros(0)
    mz_flat, rel_flat, lens = lib.gather(ids)
    off = np.concatenate([[0], np.cumsum(lens)])
    inten = rel_flat.astype(np.float64) * 100.0
    fn = ms_entropy.calculate_entropy_similarity
    scores = np.empty(ids.size)
    for k in range(ids.size):
        a, b = off[k], off[k + 1]
        scores[k] = fn(q_peaks, np.column_stack((mz_flat[a:b], inten[a:b])), ms2_tolerance_in_da=tol)
    return ids, scores


def _msms_candidates(lib, precursor_mz, tol_ppm, positive, canonical):
    """lib.msms_candidates，但加合物按现在的规范写法比：索引词表里规范化后等于 canonical 的写法都算。

    旧解析器把 GNPS 的 [M+FA]-（甲酸根，2026-03 版 34,611 条）存成了 [M+HCOOH]-；索引重建之前靠这里找回。
    """
    ids, ppm = lib.msms_candidates(precursor_mz, tol_ppm, positive)
    if canonical is None:
        return ids, ppm
    wanted = [lib.adduct_id_of(a) for a in lib.meta.get("adducts", ()) if adducts.canonical(a) == canonical]
    keep = np.isin(np.asarray(lib.adduct_id[ids]), wanted)
    return ids[keep], ppm[keep]


def _hide_massbank_copies(found):
    """[(库, 行号, Δppm)] -> (去掉 GNPS 中 MassBank 副本后的同样结构, 隐藏条数)。

    只有 MassBank 也在候选里才处理。副本要在候选里找得到原件（同一个短 ID）才隐藏：原件被前体窗口 / 离子模式 / 加合物筛掉时，
    副本留着，不然这条命中就丢了。
    """
    mb = next((f for f in found if f[0].name == lx.COPY_SOURCE), None)
    if mb is None or not len(mb[1]):
        return found, 0
    mb_rows = set(mb[1].tolist())
    out, hidden = [], 0
    for lib, ids, ppm in found:
        if lib is mb[0] or not ids.size or not lib.meta.get("counts", {}).get("massbank_copies"):
            out.append((lib, ids, ppm))
            continue
        flagged = np.flatnonzero(np.asarray(lib.copy_flag[ids]) == 1)
        drop = np.zeros(ids.size, dtype=bool)
        if flagged.size:
            infos = lib.info(ids[flagged])
            rows_of = mb[0].short_id_rows({i["copy_of"] for i in infos if i["copy_of"]})
            for pos, info in zip(flagged, infos):
                if info["copy_of"] and mb_rows.intersection(rows_of.get(info["copy_of"], ())):
                    drop[pos] = True
        hidden += int(drop.sum())
        out.append((lib, ids[~drop], ppm[~drop]))
    return out, hidden


def search(peaks, spectrum_type, libraries, metrics, threshold, *, precursor_mz=None, positive=True, adduct=None,
           precursor_tol_ppm=10.0, tolerance=0.02, top_n=100):
    """查询谱 -> 结果 dict。

    peaks [(mz, intensity)]；adduct 是用户写法或 None（= Auto，不按加合物筛）。
    返回 {"hits": 前 top_n 个保留的命中（已排序）, "n_kept", "n_candidates"（第二步打分的数量）, "hidden_copies"（隐藏的 GNPS 副本数）,
          "rank_metric"（排序用的指标）, "libraries": [...], "warnings", "notes", "metrics",
          "kind", "tolerance", "prefilter": EI 第一步信息或 None, "top_peaks": 第一名的库谱峰 (mz, 相对强度) 或 None}。
    库没装 -> lx.LibraryMissing（所选库全没装）；ValueError = 输入不对。
    """
    if spectrum_type not in SPECTRUM_TYPES:
        raise ValueError(f"Spectrum type must be one of {', '.join(SPECTRUM_TYPES)}.")
    ei = spectrum_type == "EI"
    metrics = order_metrics(metrics)
    tol = EI_TOLERANCE if ei else float(tolerance)
    if not (math.isfinite(tol) and tol > 0):
        raise ValueError("Fragment tolerance must be a positive number.")
    q_mz, q_int = sim.clean_peaks(*sim.normalize_peaks_input(peaks))
    q_prec = None
    canonical = None
    if not ei:
        if precursor_mz is None or not (math.isfinite(precursor_mz) and precursor_mz > 0):
            raise ValueError("Precursor m/z is required for an MS/MS library search.")
        q_prec = float(precursor_mz)
        if not (math.isfinite(precursor_tol_ppm) and 0 < precursor_tol_ppm < 1e5):
            raise ValueError("Precursor tolerance must be a positive number of ppm.")
        if adduct is not None:
            canonical = adducts.canonical(adduct)
            if canonical is None:
                raise ValueError(f"Unknown adduct notation: {adduct}.")
            if adduct_charge_sign(canonical) != (1 if positive else -1):
                raise ValueError(f"{adduct} does not match the {'positive' if positive else 'negative'} ion mode.")

    warnings, loaded, missing = [], [], []
    for name in dict.fromkeys(libraries):
        try:
            loaded.append(lx.load(name))
        except lx.LibraryMissing as exc:
            missing.append(str(exc))
    if not loaded:
        raise lx.LibraryMissing(" ".join(missing) if missing else "Select at least one library.")
    warnings += [f"{msg} Searching {', '.join(l.name for l in loaded)} only." for msg in missing]

    notes = {}
    kind = "EI" if ei else "MS/MS"
    prefilter = None
    hidden = 0
    # ---- 候选
    cand = []                 # (library, row, ppm, mz, 强度, 前体)
    if ei:
        q_peaks = np.column_stack([q_mz, q_int])
        parts = []
        for lib in loaded:
            ids, ent = _entropy_step(q_peaks, lib, tol)
            if not ids.size:
                warnings.append(f"{lib.name} has no EI spectra.")
            parts.append((lib, ids, ent))
        total = sum(len(p[1]) for p in parts)
        if total == 0:
            return _empty(metrics, kind, tol, loaded, warnings, notes, None)
        lib_idx = np.concatenate([np.full(len(p[1]), i) for i, p in enumerate(parts)])
        rows = np.concatenate([p[1] for p in parts])
        ent = np.concatenate([p[2] for p in parts])
        order = np.lexsort((rows, lib_idx, -ent))
        keep = order[:EI_PREFILTER_N]
        prefilter = {"total": int(total), "kept": int(len(keep)), "dropped": total > len(keep)}
        if prefilter["dropped"]:
            warnings.append(f"EI: all {total:,} library EI spectra were ranked by weighted entropy similarity; "
                            f"the top {EI_PREFILTER_N:,} were scored with the selected metrics.")
        for li in range(len(parts)):
            sel = np.sort(keep[lib_idx[keep] == li])
            if sel.size:
                lib = parts[li][0]
                for k, c in zip(sel, _cand_arrays(lib, rows[sel])):
                    cand.append((lib, int(rows[k]), None, *c))
    else:
        found = []
        for lib in loaded:
            ids, ppm = _msms_candidates(lib, q_prec, precursor_tol_ppm, positive, canonical)
            found.append((lib, ids, ppm))
        found, hidden = _hide_massbank_copies(found)
        total = sum(len(f[1]) for f in found)
        if total > MAX_MSMS_CANDIDATES:
            allppm = np.concatenate([np.abs(f[2]) for f in found])
            thr = np.partition(allppm, MAX_MSMS_CANDIDATES - 1)[MAX_MSMS_CANDIDATES - 1]
            warnings.append(f"{total:,} library spectra are within {precursor_tol_ppm:g} ppm of m/z {q_prec:.4f}; "
                            f"only the {MAX_MSMS_CANDIDATES:,} closest by precursor Δ were scored. Lower the precursor tolerance to see all of them.")
            kept = []
            for lib, ids, ppm in found:
                m = np.abs(ppm) <= thr
                kept.append((lib, ids[m], ppm[m]))
            found = kept
        for lib, ids, ppm in found:
            for i, p, c in zip(ids, ppm, _cand_arrays(lib, ids)):
                cand.append((lib, int(i), float(p), *c))
        if len(cand) > MAX_MSMS_CANDIDATES:        # 恰在阈值上的并列：截到上限
            cand.sort(key=lambda c: abs(c[2]))
            cand = cand[:MAX_MSMS_CANDIDATES]
        if not cand:
            mode = "positive" if positive else "negative"
            warnings.append(f"No library spectrum within {precursor_tol_ppm:g} ppm of m/z {q_prec:.4f}"
                            + (f" with adduct {canonical}" if canonical else "") + f" in {mode} mode.")

    n = len(cand)
    result = _empty(metrics, kind, tol, loaded, warnings, notes, prefilter)
    result["n_candidates"] = n
    result["hidden_copies"] = hidden
    counts = {}
    for c in cand:
        counts[c[0].name] = counts.get(c[0].name, 0) + 1
    for entry in result["libraries"]:
        entry["candidates"] = counts.get(entry["name"], 0)
    if not n:
        return result

    # ---- 打分
    scorer = _Scorer(q_mz, q_int, q_prec, tol, metrics, ei)
    S = np.full((n, len(metrics)), np.nan)
    M = np.zeros(n, dtype=np.int64)
    col = {m: j for j, m in enumerate(metrics)}
    for k, (lib, row, ppm, mz2, i2, p2) in enumerate(cand):
        got, M[k] = scorer.score(mz2, i2, p2)
        for name, v in got.items():
            if name in col and v is not None and math.isfinite(v):
                S[k, col[name]] = v
    deep = [m for m in metrics if m in mn.DEEP]
    for name in deep:
        reason = scorer.na_reason(name)
        if reason:
            notes[name] = f"N/A: {reason}"
    run_deep = [m for m in deep if m not in notes]
    if run_deep:
        for name, arr in _deep_scores(run_deep, q_mz, q_int, q_prec, positive, [c[3:] for c in cand], tol, notes).items():
            S[:, col[name]] = arr
    for name in metrics:
        if name in mn.MATCHMS or name in mn.ENTROPY or name in mn.CUSTOM:
            reason = scorer.na_reason(name)
            if reason:
                notes[name] = f"N/A: {reason}"
        if name not in notes and np.isnan(S[:, col[name]]).all():
            notes[name] = "N/A for every candidate"

    # How many selected metrics produced at least one number; 0 means nothing could be scored at
    # all (for example only MS2DeepScore when torch/model is missing), which is not a "no hit".
    result["n_metrics_scored"] = int(np.isfinite(S).any(axis=0).sum())

    # ---- 保留、排序
    kept = np.flatnonzero((S >= threshold).any(axis=1))          # NaN 比较为 False：N/A 不参与
    # 排序用的指标：所选指标（选项顺序）里第一个对命中有分数的；第一个若对所有命中都是 N/A（如 EI 的 ModifiedCosine）就跳过
    rank_col = next((j for j in range(len(metrics)) if np.isfinite(S[kept, j]).any()), 0)
    result["rank_metric"] = metrics[rank_col] if kept.size else None
    primary = np.where(np.isnan(S[:, rank_col]), -np.inf, S[:, rank_col])
    ppm_abs = np.array([abs(c[2]) if c[2] is not None else 0.0 for c in cand])
    lib_order = {l.name: i for i, l in enumerate(loaded)}
    lib_no = np.array([lib_order[c[0].name] for c in cand])
    rows_a = np.array([c[1] for c in cand])
    order = kept[np.lexsort((rows_a[kept], lib_no[kept], ppm_abs[kept], -M[kept], -primary[kept]))]
    result["n_kept"] = int(order.size)
    shown = order[:top_n]
    meta = {}
    for lib in loaded:
        mine = [int(k) for k in shown if cand[k][0] is lib]
        for k, info in zip(mine, lib.info([cand[k][1] for k in mine])):
            meta[k] = info
    hits = []
    for rank, k in enumerate(shown, 1):
        lib, row, ppm, *_ = cand[k]
        info = meta[int(k)]
        hits.append({"rank": rank, "library": lib.name, "row": row, "accession": info["accession"], "name": info["name"],
                     "adduct": adducts.canonical(info["adduct"]) or info["adduct"],      # 老索引的 [M+HCOOH]- 报成 [M+HCOO]-
                     "precursor_mz": float(lib.prec_mz[row]) if not ei else None, "ppm": ppm,
                     "matches": int(M[k]), "scores": {m: (None if np.isnan(S[k, col[m]]) else float(S[k, col[m]])) for m in metrics},
                     "instrument": info["instrument"], "collision_energy": info["collision_energy"],
                     "inchikey": info["inchikey"], "smiles": info["smiles"]})
    result["hits"] = hits
    if hits:
        lib0 = next(l for l in loaded if l.name == hits[0]["library"])
        result["top_peaks"] = lib0.peaks(hits[0]["row"])
        result["top_precursor"] = hits[0]["precursor_mz"]
    return result


def _empty(metrics, kind, tol, loaded, warnings, notes, prefilter):
    return {"hits": [], "n_kept": 0, "n_candidates": 0, "hidden_copies": 0, "rank_metric": None, "metrics": metrics, "kind": kind, "tolerance": tol,
            "n_metrics_scored": 0,
            "libraries": [{"name": l.name, "version": l.version, "label": l.label(), "candidates": 0} for l in loaded],
            "warnings": warnings, "notes": notes, "prefilter": prefilter, "top_peaks": None, "top_precursor": None}
