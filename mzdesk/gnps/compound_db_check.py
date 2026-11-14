#!/usr/bin/env python3
"""
compound_db_check.py — 用 fastMASST (GNPS2) 检索一张 MS/MS 谱图在公共数据中的分布

检索内容:
  1. fastMASST        (谱图全库检索 library=metabolomicspanrepo_index_latest:
                       GNPS/MassIVE + MetaboLights + Metabolomics Workbench + NORMAN)
  2. GNPS 参考谱库    (fastMASST library=gnpslibrary, MSI Level-2 身份确认;
                       命中按 cosine 排序, Suspect 条目单列)

两项检索都只保留匹配峰数 >= --min-peaks (默认 3) 的命中

用法示例:
  # 直接给 MS/MS 谱图
  python compound_db_check.py --name "MG 16:0" \
      --precursor 331.284 \
      --peaks "[[60.081,98618],[71.0856,16836],[85.02837,860576]]"
  (注: 此示例的 60.081 / 85.028 是酰基肉碱型碎片, 并非 MG 诊断碎片; 仅作命令格式示例)

  # 从 in_house.parquet 按行号自动提取 (需要 duckdb)
  python compound_db_check.py --from-parquet /path/to/in_house.parquet --row 87

  # 不检索谱图, 只按 Pan-ReDU 样本元数据找公开数据集
  # (例: 小鼠大脑, 符合条件的文件数 > 12, 同时有正负离子模式)
  python compound_db_check.py --find-datasets --organism "Mus musculus" \
      --body-part brain --min-files 13 --polarity both

输出: 终端汇总表 + JSON 结果文件 (--out)

---------------- 2026-10-15 新增功能 (旧命令的默认值与 JSON 字段不变) ----------------
旧备份: compound_db_check.py.orig_20261015

  # 只查全库, 自定义容差, 保存全部命中 (<out>.hits.tsv)
  python compound_db_check.py --name "MG 18:1 [M+H]+" --precursor 357.2999 \
      --peaks "[[247.242,100],[265.252,75],[339.289,60],[121.101,72]]" \
      --library repository --pm-tol 0.02 --frag-tol 0.02 --min-peaks 3 \
      --keep-all --out mg181/mg181.json

  # 诊断离子核对 (取回命中谱图, 判断诊断离子是否存在) + Pan-ReDU 元数据关联与富集
  python compound_db_check.py ... --diagnostic "265.2526,247.2420,339.2894" \
      --diag-required 265.2526 --diag-rule 2-of-3 --diag-tol 0.01 \
      --verify-max 300 --redu-join

  # 用 USI 作查询谱图 (仓库谱图或 GNPS 参考谱库条目, 经 metabolomics-usi.gnps2.org 取峰表)
  python compound_db_check.py --name "MG 20:4 lib" \
      --usi mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00000000000 --library repository

  # 重放: 不再调用 fastMASST, 用上次保存的全部命中换一套诊断离子规则 / --min-peaks / ReDU 背景重新分析
  python compound_db_check.py --from-hits mg181/mg181.hits.tsv --min-peaks 4 \
      --diagnostic "265.2526,247.2420,339.2894" --diag-rule all --diag-min-rel 5 --redu-join \
      --out mg181/mg181_strict.json

  # 批量模式: 依次检索多个查询 (不并发), 每个查询一个子目录, 汇总表 <out>/batch_summary.tsv
  python compound_db_check.py --queries queries.json --out mg_batch --keep-all --redu-join
  python compound_db_check.py --queries queries.json --out mg_batch --dry-run   # 只检查查询文件

  # MassQL 本地检索: 在 mzML/mzXML/MGF 中按碎片模式找谱 (不需要 fastMASST)
  python compound_db_check.py \
      --massql "QUERY scaninfo(MS2DATA) WHERE MS2PROD=265.2526:TOLERANCEMZ=0.01 AND MS2PROD=247.2420:TOLERANCEMZ=0.01" \
      --massql-input raw_data/ --out mg_massql.json --massql-mgf mg_massql.mgf

新增选项 (语义):
  --library       both (默认, 与旧版相同: 全库 + 参考谱库) | repository | library
                  | 其他 fastMASST 索引名 (逗号分隔可多个)
  --pm-tol / --frag-tol   前体 / 碎片 m/z 容差 (Da, 默认 0.05 / 0.05, 与旧版相同)
  --max-wait      单次检索最长等待秒数 (默认 600); 轮询间隔固定 8 s
  --usi           查询谱图来自 USI (未给 --precursor 时用 USI 的 precursor_mz)
  --from-hits     重放上次 --keep-all 的 <out>.hits.tsv (需同名 JSON): 不调用 fastMASST; 可改 --min-peaks、
                  调高 --cosine、诊断离子与 ReDU 参数; 容差 (--pm-tol / --frag-tol) 沿用原检索, 无法事后修改
  --query-min-rel 查询谱图预处理: 去掉低于基峰 X% 的峰 (默认不过滤)
  --keep-all      保存 API 返回的每一条命中 (全部字段) 到 <out>.hits.tsv (含未达 --min-peaks 的,
                  pass_min_peaks 列标记); --diagnostic / --redu-join 自动打开
  --diagnostic    诊断离子 m/z 列表; 对命中 (仅 pass_min_peaks) 经 USI 取回谱图, 记录每个诊断离子
                  的相对强度 (% 基峰) 与实测 m/z, 并按 --diag-rule 给出 diag_pass
  --diag-required 必需离子 (自动并入 --diagnostic), 缺任何一个即 fail
  --diag-tol      诊断离子 m/z 容差 (Da, 默认 0.01)
  --diag-min-rel  诊断离子最低相对强度 (% 基峰, 默认 0 = 有峰即算)
  --diag-rule     all (默认) | any | k (至少 k 个) | k-of-n (n 必须等于诊断离子数)
  --verify-max    每次检索最多取回的谱图数 (默认 200); 超出时按数据集分层轮转抽取可复现子集
                  (--verify-seed, 默认 0), 未核对的命中写入 <out>.not_verified.tsv (原因: not_selected /
                  fetch_failed / invalid_scan). 取谱失败 (5xx/网络) 退避重试 3, 8, 20 s, 再整体暂停 20 s
                  重试一轮; USI 扫描号为 4294967295 (索引无效值) 的命中无法取谱, 不参与抽样 (invalid_scan)
  --usi-cache     谱图缓存目录 (默认: 输出 JSON 同目录下 usi_cache/; 批量模式共用 <out>/usi_cache)
  --usi-rate      USI 取谱限速 (次/秒, 默认 2.5, 上限 3)
  --redu-join     命中文件关联本地 Pan-ReDU 缓存, 写逐命中元数据 (hits.tsv 的 redu_* 列)、
                  <out>.redu_summary.tsv (按 SampleType / SampleTypeSub1 / NCBITaxonomy /
                  UBERONBodyPartName / DOIDCommonName / HealthStatus 汇总, 文件级与数据集级富集)、
                  <out>.redu_by_dataset.tsv
  --redu-bg       富集背景 (命中文件也只统计落在背景内的): pos_lc_or_unknown (默认: 正离子 / 正负交替 /
                  极性未标注, 且非 GC/EI) | pos_lc (只要已标注正离子或交替, 非 GC/EI) | neg_lc_or_unknown
                  | neg_lc | all. 默认不用严格 pos_lc 的原因: Pan-ReDU 中 46% 文件极性为 missing value,
                  GNPS/MassIVE 数据集尤甚 (2026-10-15 测试: MG 16:0 示例 312 个已匹配命中文件只有 24 个在
                  pos_lc 内; MG 18:0;O 测试 14352 个中 6744 个). 已标注负离子的命中文件数记在 JSON
                  (matched_hit_files_polarity)
  --queries       批量模式 JSON; --out 为输出目录; --pause 两次检索之间的间隔秒数 (默认 5);
                  --resume 跳过已成功完成的查询; --dry-run 只解析不检索

fastMASST 命中字段 (2026-10-15 实测, hits.tsv 原样保留): USI, Dataset, Cosine, Matching Peaks, Delta Mass
  (= 查询前体 - 命中前体, 2 位小数), Unit Delta Mass, Charge, Status (HasID/NoID), Query Scan, Index UnitPM,
  Index IdxInUnitPM; 参考谱库另有 GNPSLibraryAccession, CompoundName, Adduct. API 不返回命中前体 m/z 与
  扫描号: hits.tsv 的 scan 取自 USI, est_hit_precursor_mz = 查询前体 - Delta Mass (约 ±0.005),
  核对后 hit_precursor_mz 为 resolver 返回的真实值.

USI 归一化 (命中 -> Pan-ReDU 文件): 去掉 USI 末尾的 :scan:N / :index:N / :nativeId:... 得到文件级 USI
  mzspec:<数据集>:<文件路径>, 然后依次尝试: 1 完全相同 (exact) 2 忽略大小写 (case_insensitive)
  3 去扩展名后相同 (no_extension, 如 .mzML/.mzXML/.raw/.CDF/.d) 4 同数据集内同文件名 (忽略目录;
  唯一 basename_unique, 多个取路径后缀最长相同者 basename_ambiguous). 仍失败: 数据集在 ReDU 中
  (dataset_only) / 数据集不在 ReDU (no_dataset); GNPS 参考谱库命中记为 library (不适用).
  数据集前缀 MSV / MTBLS / ST / NORMAN-<uuid> 都取自 USI 第二段, 比较时不区分大小写.

ReDU 覆盖率是主要限制: 很多数据集 (尤其 2025-2026 新上传的 MSV) 不在 Pan-ReDU 中, 或只登记了部分文件.
  2026-10-15 测试 (MG 18:0;O [M+H-H2O]+): 61% 命中 / 60% 文件 / 63% 数据集未匹配, 其中 no_dataset 占大头.
  未匹配的命中不进入元数据汇总, 比例写在 JSON (unmatched_frac_*) 与终端输出.

富集指标 (redu_summary.tsv, hit_set = all 全部命中 / verified 已核对 / diag_pass 诊断通过;
  只用落在背景内的命中文件):
  file_rate   = 含该值且有命中的文件数 / 背景中含该值的文件数
  file_enrich = (该值命中文件 / 全部命中文件) / (该值背景文件 / 全部背景文件)
  ds_rate     = 该值有命中的数据集数 / 背景中有该值的数据集数
  ds_enrich   = ds_rate / (有命中的数据集数 / 背景数据集数)
  dsw_*       = 每个数据集权重为 1 的组成比例 (命中 vs 背景), dsw_enrich = 两者之比
"""
import argparse, ast, collections, csv, gzip, hashlib, io, json, os, re, time, sys, urllib.error, urllib.parse, urllib.request

UA = {"User-Agent": "Mozilla/5.0 (compound-db-check; research)"}
FASST = "https://api.fasst.gnps2.org"
USI_JSON = "https://metabolomics-usi.gnps2.org/json/?usi1="
REDU_URL = "https://redu.gnps2.org/dump"  # Pan-ReDU 全表 (逐文件元数据, 约 1.7 GB, 每日更新)
# Pan-ReDU metadata is not shipped with mzDesk (no redistribution licence): it is downloaded on first use
# (local version: the "redu" component, `mzdesk setup redu`; server: deploy/install-engines.sh) into the writable cache.
REDU_SLIM = os.environ.get("MZDESK_GNPS_REDU_CACHE", os.path.join(
    os.environ.get("MZDESK_CACHE_DIR", os.path.expanduser("~/.mzdesk/cache")), "gnps", "redu_slim.tsv.gz"))
REDU_MISSING = ("The Pan-ReDU metadata is not downloaded yet. Local version: run the page again and choose Download, or run "
                "`mzdesk setup redu` (about 1.7 GB is streamed from GNPS2; 8 MB is kept).")
REDU_COLS = ["ATTRIBUTE_DatasetAccession", "DataSource", "USI", "SampleType", "SampleTypeSub1",
             "NCBITaxonomy", "UBERONBodyPartName", "DOIDCommonName", "HealthStatus", "LifeStage",
             "BiologicalSex", "Country", "MassSpectrometer", "IonizationSourceAndPolarity",
             "ChromatographyAndPhase", "SampleExtractionMethod"]
TOOL_VERSION = "2026-10-15"


def get(url, timeout=40):
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read().decode("utf-8", "replace")
    except Exception as e:
        return f"__ERR__ {e}"


def get_json(url, timeout=40):
    t = get(url, timeout)
    if t.startswith("__ERR__"):
        return {"error": t}
    try:
        return json.loads(t)
    except Exception:
        return {"error": "not-json", "raw": t[:300]}


def post_json(url, payload, timeout=60):
    req = urllib.request.Request(url, data=json.dumps(payload).encode(),
                                 headers={**UA, "Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return json.loads(r.read().decode())
    except Exception as e:
        return {"error": str(e)}


# ---------- fastMASST ----------
SUSPECT = "Suspect related to"  # GNPS suspect 谱库条目 (传播注释, 非参考谱)
REPO_LIB, REF_LIB = "metabolomicspanrepo_index_latest", "gnpslibrary"
# fastMASST 命中字段 (2026-10-15 实测; 参考谱库命中另有 GNPSLibraryAccession / CompoundName / Adduct)
# Delta Mass = 查询前体 m/z - 命中前体 m/z (保留 2 位小数); Unit Delta Mass = 其整数部分
API_FIELDS = ["USI", "Dataset", "Cosine", "Matching Peaks", "Delta Mass", "Unit Delta Mass", "Charge",
              "Status", "Query Scan", "Index UnitPM", "Index IdxInUnitPM",
              "GNPSLibraryAccession", "CompoundName", "Adduct"]


def library_hits(results, top=5):
    best = {}  # 同名条目合并, 保留 cosine 最高的一条
    for x in results:
        n = x.get("CompoundName")
        if not n:
            continue
        b = best.setdefault(n, {"name": n, "cosine": -1, "n_entries": 0})
        b["n_entries"] += 1
        if x.get("Cosine", 0) > b["cosine"]:
            b.update(cosine=x.get("Cosine", 0), matching_peaks=x.get("Matching Peaks"),
                     adduct=x.get("Adduct"), accession=x.get("GNPSLibraryAccession"))
    hits = sorted(best.values(), key=lambda b: -b["cosine"])
    ref = [b for b in hits if not b["name"].startswith(SUSPECT)]
    sus = [b for b in hits if b["name"].startswith(SUSPECT)]
    return {"n_reference": sum(b["n_entries"] for b in ref), "n_suspect": sum(b["n_entries"] for b in sus),
            "library_hits": ref[:top], "suspect_hits": sus[:top]}


def summarise_raw(raw, library, min_peaks, tid):
    """API 原始命中 -> 汇总 dict (旧版 fasst_search 的输出字段)"""
    results = [x for x in raw if (x.get("Matching Peaks") or 0) >= min_peaks]
    out = {"task_id": tid, "n_matches": len(results), "n_before_filter": len(raw), "min_peaks": min_peaks}
    if results:
        ds = collections.Counter(x.get("Dataset") for x in results)
        out["top_datasets"] = ds.most_common(10)
        out["max_cosine"] = max(x.get("Cosine", 0) for x in results)
        if library == "gnpslibrary":
            out.update(library_hits(results))
    return out


def fasst_search(peaks, precursor, library, cosine=0.7, min_peaks=3, pm_tol=0.05, frag_tol=0.05,
                 max_wait_s=600, poll_s=8, return_raw=False, max_poll_errors=3):
    """返回汇总 dict (字段与旧版相同). return_raw=True 时返回 (汇总 dict, info),
    info = {"task_id", "raw" (API 原始命中列表, 失败时 None), "elapsed_s", "n_poll_errors"}.
    轮询不快于 8 s; 状态/结果请求网络瞬断时最多连续容忍 max_poll_errors 次.
    max_wait_s 是提交后等结果的截止时间 (time.monotonic, 请求本身的耗时也算在内); 最后一次请求可能越过截止时间, 所以总等待至多是
    max_wait_s 加一次请求的时间."""
    poll_s = max(poll_s, 8)
    t0 = time.monotonic()
    info = {"task_id": None, "raw": None, "elapsed_s": None, "n_poll_errors": 0}

    def done(out):
        info["elapsed_s"] = round(time.monotonic() - t0, 1)
        return (out, info) if return_raw else out

    qs = json.dumps({"peaks": [[float(m), float(i)] for m, i in peaks],
                     "precursor_mz": float(precursor)})
    r = post_json(f"{FASST}/search", {"library": library, "query_spectrum": qs, "analog": "No",
                                      "pm_tolerance": pm_tol, "fragment_tolerance": frag_tol,
                                      "cosine_threshold": cosine})
    if "id" not in r:
        return done({"error": r})
    tid = r["id"]
    info["task_id"] = tid
    submitted = time.monotonic()
    deadline = submitted + max_wait_s
    n_err = 0
    while time.monotonic() < deadline:
        time.sleep(min(poll_s, max(deadline - time.monotonic(), 0)))   # the last look is at the deadline
        st = get_json(f"{FASST}/search/status/{tid}")
        s = st.get("status")
        if s == "COMPLETED":
            res = get_json(f"{FASST}/search/result/{tid}", timeout=120)
            if not isinstance(res, dict) or not isinstance(res.get("results"), list):
                if n_err < max_poll_errors:  # 取结果瞬断: 稍后重试
                    n_err += 1
                    info["n_poll_errors"] += 1
                    continue
                return done({"error": res, "task_id": tid})  # 取结果失败, 不能当作 0 匹配
            raw = res["results"]
            info["raw"] = raw
            return done(summarise_raw(raw, library, min_peaks, tid))
        if s in ("PENDING", "RUNNING"):
            n_err = 0
            continue
        if "error" in st and s is None and n_err < max_poll_errors:  # 网络瞬断, 继续等
            n_err += 1
            info["n_poll_errors"] += 1
            continue
        return done({"error": st})
    return done({"error": f"timeout after {time.monotonic() - submitted:.0f}s waiting for the search (limit {max_wait_s}s)", "task_id": tid})


# ---------- Pan-ReDU 数据集检索 ----------
REDU_TRIES = 6        # a dropped connection is resumed from where it stopped (HTTP Range), up to this many times


class _ResumableBody(io.RawIOBase):
    """The Pan-ReDU response as one byte stream: when GNPS2 drops the connection early (it does, after a few minutes), it is opened again
    with Range from the byte where it stopped, as long as the file is the same (If-Range with its ETag). pos/total give the progress."""

    def __init__(self, url, timeout=120, tries=REDU_TRIES, wait=5):
        self.url, self.timeout, self.tries, self.wait = url, timeout, tries, wait
        self.r = urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout)
        self.total = int(self.r.headers.get("Content-Length") or 0) or None
        self.etag = self.r.headers.get("ETag")
        self.pos = 0
        self.resumed = 0

    def readable(self):
        return True

    def readinto(self, b):
        while True:
            try:
                n = self.r.readinto(b)
            except OSError:                                # a reset or a timeout in the middle: resumed like an early end
                n = 0
            if n:
                self.pos += n
                return n
            if self.total is None or self.pos >= self.total:
                return 0
            self._resume()

    def _resume(self):
        if self.resumed >= self.tries:
            raise IOError(f"下载不完整, 还差 {self.total - self.pos} 字节 (已续传 {self.resumed} 次), 请重试")
        self.resumed += 1
        self.r.close()
        print(f"[ReDU]   连接在 {self.pos >> 20} MB 处断开, 从这里接着下 (第 {self.resumed} 次) ...", flush=True)
        time.sleep(self.wait)
        headers = dict(UA, Range=f"bytes={self.pos}-", **({"If-Range": self.etag} if self.etag else {}))
        try:
            self.r = urllib.request.urlopen(urllib.request.Request(self.url, headers=headers), timeout=self.timeout)
        except OSError:
            return                                         # the next read finds nothing and tries again (counted)
        start = (self.r.headers.get("Content-Range") or "").removeprefix("bytes ").split("-")[0]
        if self.r.status != 206 or start != str(self.pos):  # GNPS2 has a newer table now (it is updated daily): start again
            self.r.close()
            raise IOError("Pan-ReDU 表在下载途中更新了, 无法续传, 请重试")

    def close(self):
        self.r.close()
        super().close()


def redu_download():
    """流式下载 Pan-ReDU 全表, 只保留 REDU_COLS 写成精简缓存 (原表不落盘); 断线时从断点续传"""
    os.makedirs(os.path.dirname(REDU_SLIM), exist_ok=True)
    tmp = REDU_SLIM + ".part"
    with _ResumableBody(REDU_URL) as body, gzip.open(tmp, "wt", newline="", encoding="utf-8") as f:
        size = f"{body.total / 1e9:.1f} GB" if body.total else "about 1.7 GB"
        print(f"[ReDU] 下载 {REDU_URL} ({size}, 只保留 {len(REDU_COLS)} 列) ...", flush=True)
        rd = csv.reader(io.TextIOWrapper(io.BufferedReader(body, 1 << 20), encoding="utf-8", errors="replace", newline=""), delimiter="\t")
        header = next(rd)
        idx = [header.index(c) for c in REDU_COLS]
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(REDU_COLS)
        n = bad = 0
        for row in rd:
            if len(row) != len(header):
                bad += 1
                continue
            w.writerow([row[i] for i in idx])
            n += 1
            if n % 200000 == 0:
                pct = f", 约 {100 * body.pos // body.total}%" if body.total else ""
                print(f"[ReDU]   已处理 {n} 行 ({body.pos >> 20} MB{pct}) ...", flush=True)
    os.replace(tmp, REDU_SLIM)
    print(f"[ReDU] 完成: {n} 行" + (f", 跳过 {bad} 行格式异常" if bad else ""))


def redu_installed():
    return os.path.isfile(REDU_SLIM)


def redu_path(refresh=False, download=True):
    """The downloaded Pan-ReDU cache. Missing: the command line downloads it (download=True); the web app never starts a 1.7 GB download
    inside a request (download=False) and raises FileNotFoundError(REDU_MISSING) instead."""
    if refresh:
        redu_download()
    if os.path.isfile(REDU_SLIM):
        return REDU_SLIM
    if not download:
        raise FileNotFoundError(REDU_MISSING)
    redu_download()
    return REDU_SLIM


def redu_cache_date(path=None):
    path = path or redu_path()
    return time.strftime("%Y-%m-%d", time.localtime(os.path.getmtime(path)))


def redu_load(refresh=False):
    with gzip.open(redu_path(refresh), "rt", newline="", encoding="utf-8") as f:
        rd = csv.reader(f, delimiter="\t")
        header = next(rd)
        return header, list(rd)


def polarity_of(v):
    v = v.lower()
    for key, p in [("alternating", "alt"), ("positive", "pos"), ("negative", "neg")]:
        if key in v:
            return p
    return "unk" if v in ("", "missing value") else "other"  # other: 无正负之分的离子源, 如 GC-EI


def find_datasets(header, rows, filters, min_files=1, polarity="any"):
    """filters: [(列名, [关键词, ...])], 列之间 AND, 同列关键词之间 OR, 不区分大小写的子串匹配
    (NCBITaxonomy 的纯数字关键词按 taxid 精确匹配). 文件数与极性都只统计符合 filters 的文件."""
    col = {c: i for i, c in enumerate(header)}
    flt = [(col[c], c == "NCBITaxonomy", [t.strip().lower() for t in terms if t.strip()]) for c, terms in filters]
    i_ds, i_src, i_usi, i_pol = (col[c] for c in ["ATTRIBUTE_DatasetAccession", "DataSource", "USI",
                                                  "IonizationSourceAndPolarity"])
    seen, ds = set(), {}
    matched = {c: collections.Counter() for c, _ in filters}
    for r in rows:
        if r[i_usi] in seen:
            continue
        ok = True
        for i, is_tax, terms in flt:
            v = r[i].lower()
            if not any(v.startswith(t + "|") if is_tax and t.isdigit() else t in v for t in terms):
                ok = False
                break
        if not ok:
            continue
        seen.add(r[i_usi])
        for c in matched:
            matched[c][r[col[c]]] += 1
        d = ds.setdefault(r[i_ds], {"dataset": r[i_ds], "source": r[i_src], "n_files": 0,
                                    "pos": 0, "neg": 0, "alt": 0, "unk": 0, "other": 0})
        d["n_files"] += 1
        d[polarity_of(r[i_pol])] += 1

    out = {"n_files_matched": len(seen), "n_datasets_matched": len(ds),
           "matched_values": {c: m.most_common(15) for c, m in matched.items()},
           "passed": [], "undetermined": [], "n_failed_min_files": 0, "n_failed_polarity": 0}
    for d in sorted(ds.values(), key=lambda d: -d["n_files"]):
        if d["n_files"] < min_files:
            out["n_failed_min_files"] += 1
            continue
        has = {"positive": d["pos"] or d["alt"], "negative": d["neg"] or d["alt"]}
        need = {"any": [], "both": ["positive", "negative"]}.get(polarity, [polarity])
        if all(has[p] for p in need):
            out["passed"].append(d)
        elif d["unk"]:  # 有文件未标注极性, 无法判断
            out["undetermined"].append(d)
        else:
            out["n_failed_polarity"] += 1
    return out


def run_find_datasets(args, ap):
    filters = [(c, v.split(",")) for c, v in [("NCBITaxonomy", args.organism),
                                              ("UBERONBodyPartName", args.body_part)] if v]
    for w in args.where:
        c, _, v = w.partition("=")
        if c not in REDU_COLS or not v:
            ap.error(f"--where 格式为 列=值, 可用列: {', '.join(REDU_COLS)}")
        filters.append((c, v.split(",")))
    if not filters:
        ap.error("--find-datasets 至少需要一个条件: --organism / --body-part / --where")

    header, rows = redu_load(args.redu_refresh)
    cache_date = redu_cache_date()
    res = find_datasets(header, rows, filters, args.min_files, args.polarity)
    res = {"query": {"filters": filters, "min_files": args.min_files, "polarity": args.polarity},
           "redu_cache_date": cache_date, "redu_n_files": len(rows), **res}

    def show(ds, top):
        for d in ds[:top]:
            print(f"  {d['dataset']:<13}| {d['source']:<12} | {d['n_files']:>5} 文件 | "
                  f"正 {d['pos']} 负 {d['neg']} 交替 {d['alt']} 未标注 {d['unk']} 其他 {d['other']}")
        if len(ds) > top:
            print(f"  ... 其余 {len(ds) - top} 个见 JSON")

    print("\n" + "=" * 60)
    print(f"  Pan-ReDU 数据集检索 (缓存 {cache_date}, 共 {len(rows)} 个文件)")
    print("  条件: " + "; ".join(f"{c} 含 {'/'.join(t)}" for c, t in filters)
          + f"; 文件数 >= {args.min_files}; 极性 {args.polarity}")
    print("=" * 60)
    print(f"符合条件的文件: {res['n_files_matched']} 个, 分布在 {res['n_datasets_matched']} 个数据集")
    for c, m in res["matched_values"].items():
        print(f"  {c} 实际匹配到: " + ", ".join(f"{v} ({n})" for v, n in m[:6]))
    if not res["n_files_matched"]:
        print("  无匹配: 物种请用拉丁名或 taxid (Mus musculus / 10090), 组织请用 UBERON 英文名 (brain)")
    print(f"通过: {len(res['passed'])} 个数据集")
    show(res["passed"], 30)
    print(f"待定 (文件数达标, 但有文件未标注极性, 无法判断): {len(res['undetermined'])} 个数据集")
    show(res["undetermined"], 10)
    print(f"未通过: 文件数不足 {res['n_failed_min_files']} 个, 极性不符 {res['n_failed_polarity']} 个")
    print("=" * 60)

    out = args.out or "redu_datasets.json"
    with open(out, "w") as f:
        json.dump(res, f, indent=2, ensure_ascii=False)
    print(f"\n结果已保存: {out}")


# ---------- 汇总打印 ----------
SEARCH_LABEL = {"masst_repository": "fastMASST 全库", "masst_library": "GNPS 参考谱库"}


def print_summary(name, results):
    print("\n" + "=" * 60)
    print(f"  化合物: {name}")
    print("=" * 60)
    keys = [k for k in ("masst_repository", "masst_library") if k in results]
    keys += [k for k in results if k.startswith("masst_") and k not in keys]  # --library 自定义索引
    for key in keys:
        label = SEARCH_LABEL.get(key, f"fastMASST {key[6:]}")
        m = results[key]
        if "error" in m:
            print(f"{label}: 失败 {m['error']}")
            continue
        head = f"{label}: {m['n_matches']} 匹配 (匹配峰 >={m['min_peaks']}, 过滤前 {m['n_before_filter']})"
        if key != "masst_library":
            print(f"{head}, max cosine {m.get('max_cosine', '-')}")
            continue
        print(head)
        for hk, nk, title in [("library_hits", "n_reference", "参考条目"),
                              ("suspect_hits", "n_suspect", "Suspect 条目 (传播注释, 非参考谱)")]:
            print(f"  {title}: {m.get(nk, 0)} 条")
            for h in m.get(hk, []):
                print(f"    cos {h['cosine']:.2f} | {h['matching_peaks']} 峰 | {h['name'][:70]} (x{h['n_entries']})")
    print("=" * 60)


# ======================================================================================
# 2026-10-15 新增: 全部命中 / 诊断离子核对 / Pan-ReDU 关联 / 批量模式
# ======================================================================================
# 新选项的内置默认值 (命令行未给时用这里的值; 旧选项 --cosine / --min-peaks 的默认值不变)
OPT_DEFAULTS = {"library": "both", "cosine": 0.7, "min_peaks": 3, "pm_tol": 0.05, "frag_tol": 0.05,
                "max_wait": 600, "query_min_rel": None, "keep_all": False, "diagnostic": None,
                "diag_required": None, "diag_tol": 0.01, "diag_min_rel": 0.0, "diag_rule": "all",
                "verify_max": 200, "verify_seed": 0, "redu_join": False, "redu_bg": "pos_lc_or_unknown",
                "usi_rate": 2.5}
NEW_OPTS = ["library", "pm_tol", "frag_tol", "max_wait", "query_min_rel", "keep_all", "diagnostic",
            "diag_required", "diag_tol", "diag_min_rel", "diag_rule", "verify_max", "verify_seed",
            "redu_join", "redu_bg", "usi_rate", "usi_cache", "usi", "from_hits"]
QUERY_KEYS = {"name", "inchikey", "precursor", "peaks", "mgf", "usi", "from_hits", "notes"} | set(OPT_DEFAULTS)
SUMMARY_COLS = ["SampleType", "SampleTypeSub1", "NCBITaxonomy", "UBERONBodyPartName", "DOIDCommonName",
                "HealthStatus"]
REDU_BG_MODES = ["pos_lc_or_unknown", "pos_lc", "neg_lc_or_unknown", "neg_lc", "all"]
EXT_RE = re.compile(r"\.(mzml\.gz|mzxml\.gz|mzml|mzxml|mzdata\.xml|mzdata|mzmlb|mz5|raw|cdf|wiff2?|d|mgf)$",
                    re.I)
INDEX_RE = re.compile(r":(scan|index|nativeId|accession):")
INVALID_SCAN = "4294967295"  # fastMASST 索引中的无效扫描号 (2^32-1), resolver 无法取谱 (HTTP 404/500)


def resolve_libraries(spec):
    """--library -> [(结果键, fastMASST 索引名)]"""
    spec = (spec or "both").strip()
    if spec == "both":
        return [("masst_repository", REPO_LIB), ("masst_library", REF_LIB)]
    out = []
    for s in [x.strip() for x in spec.split(",") if x.strip()]:
        if s in ("repository", "repo", REPO_LIB):
            out.append(("masst_repository", REPO_LIB))
        elif s in ("library", "lib", REF_LIB):
            out.append(("masst_library", REF_LIB))
        else:
            out.append(("masst_" + re.sub(r"[^A-Za-z0-9]+", "_", s), s))
    if len({k for k, _ in out}) != len(out):
        raise ValueError(f"--library 有重复: {spec}")
    return out


def parse_usi(usi):
    """'mzspec:<集合>:<文件>[:scan|index|nativeId|accession:<值>][:<注释>]' -> (集合, 文件, 索引类型, 索引)"""
    s = (usi or "").strip()
    if not s.lower().startswith("mzspec:"):
        return None, None, None, None
    coll, _, rest = s[7:].partition(":")
    last = None
    for last in INDEX_RE.finditer(rest):
        pass
    if last is None:
        return coll, rest, None, None
    return coll, rest[:last.start()], last.group(1), rest[last.end():].split(":")[0]


def file_usi_of(usi):
    coll, run, _, _ = parse_usi(usi)
    return f"mzspec:{coll}:{run}" if coll is not None else None


def write_tsv(path, rows, cols):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            out = []
            for c in cols:
                v = r.get(c)
                out.append("" if v is None else (int(v) if isinstance(v, bool) else v))
            w.writerow(out)


def mz_label(x):
    return ("%.4f" % x).rstrip("0").rstrip(".")


def parse_mz_list(v):
    if v is None or v == "":
        return []
    if isinstance(v, (int, float)):
        return [float(v)]
    if isinstance(v, str):
        return [float(t) for t in v.replace(";", ",").split(",") if t.strip()]
    return [float(t) for t in v]


def parse_diag_rule(rule, n):
    """'all' | 'any' | 'k' | 'k-of-n' -> 至少需要的诊断离子数 k"""
    r = str(rule).strip().lower()
    if r == "all":
        k = n
    elif r == "any":
        k = 1
    elif re.fullmatch(r"\d+", r):
        k = int(r)
    elif re.fullmatch(r"\d+-of-\d+", r):
        k, nn = (int(t) for t in r.split("-of-"))
        if nn != n:
            raise ValueError(f"--diag-rule {rule}: n={nn} 与诊断离子数 {n} 不一致")
    else:
        raise ValueError(f"--diag-rule 无法解析: {rule} (可用 all / any / k / k-of-n)")
    if not 1 <= k <= n:
        raise ValueError(f"--diag-rule {rule}: k={k} 超出 1..{n}")
    return k


class UsiFetcher:
    """经 metabolomics-usi.gnps2.org/json 取谱图: 限速 (<= 3 次/秒), 重试, 本地缓存 (只缓存成功结果)"""

    BACKOFF = (3, 8, 20)  # 秒; 2026-10-15 实测 resolver 偶发 HTTP 500, 几秒后即可恢复

    def __init__(self, cache_dir, rate=2.5, retries=4, timeout=60):
        self.cache_dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)
        self.min_interval = 1.0 / min(float(rate), 3.0)
        self.retries, self.timeout = retries, timeout
        self.last = 0.0
        self.stats = collections.Counter()

    def _path(self, usi):
        return os.path.join(self.cache_dir, hashlib.sha1(usi.encode()).hexdigest() + ".json")

    def fetch(self, usi):
        """-> (谱图 dict 或 None, 状态字符串)"""
        p = self._path(usi)
        if os.path.exists(p):
            try:
                with open(p) as f:
                    d = json.load(f)
                if d.get("usi") == usi and isinstance(d.get("data", {}).get("peaks"), list):
                    self.stats["cache"] += 1
                    return d["data"], "ok"
            except Exception:
                pass
        url = USI_JSON + urllib.parse.quote(usi, safe="")
        err = None
        for attempt in range(self.retries):
            wait = self.min_interval - (time.time() - self.last)
            if wait > 0:
                time.sleep(wait)
            self.last = time.time()
            self.stats["requests"] += 1
            try:
                req = urllib.request.Request(url, headers=UA)
                with urllib.request.urlopen(req, timeout=self.timeout) as r:
                    data = json.loads(r.read().decode("utf-8", "replace"))
                if not isinstance(data.get("peaks"), list):
                    raise ValueError("no peaks in response")
                tmp = p + ".part"
                with open(tmp, "w") as f:
                    json.dump({"usi": usi, "fetched": time.strftime("%Y-%m-%d %H:%M:%S"), "data": data}, f)
                os.replace(tmp, p)
                self.stats["network_ok"] += 1
                return data, "ok"
            except urllib.error.HTTPError as e:
                err = f"HTTP {e.code}"
                self.stats["request_error " + err] += 1
                if 400 <= e.code < 500 and e.code != 429:  # 请求本身有误, 重试无用
                    break
            except Exception as e:
                err = f"{type(e).__name__}: {str(e)[:120]}"
                self.stats["request_error " + type(e).__name__] += 1
            if attempt < self.retries - 1:
                time.sleep(self.BACKOFF[min(attempt, len(self.BACKOFF) - 1)])  # 退避
        self.stats["calls_failed_after_retries"] += 1  # 可能在第 2 轮重试中恢复; 最终失败数见 n_fetch_failed
        return None, f"fetch_failed: {err}"


def select_for_verification(cands, verify_max, seed):
    """超过 verify_max 时按数据集分层: 数据集与数据集内命中都按 sha1(seed|键) 排序 (可复现的伪随机),
    再在数据集之间轮转抽取, 保证先覆盖尽可能多的数据集"""
    if verify_max is None or verify_max < 0 or len(cands) <= verify_max:
        return list(cands)

    def h(s):
        return hashlib.sha1(f"{seed}|{s}".encode()).hexdigest()

    by_ds = collections.defaultdict(list)
    for r in cands:
        by_ds[r.get("Dataset") or ""].append(r)
    queues = [sorted(v, key=lambda r: h(r["USI"])) for _, v in sorted(by_ds.items(), key=lambda kv: h(kv[0]))]
    sel, i = [], 0
    while len(sel) < verify_max:
        progressed = False
        for q in queues:
            if i < len(q):
                sel.append(q[i])
                progressed = True
                if len(sel) >= verify_max:
                    break
        if not progressed:
            break
        i += 1
    return sel


def diag_eval(peaks, diag, tol, min_rel):
    """每个诊断离子: 容差内最强峰的 (实测 m/z, 相对强度 % 基峰); 相对强度 >= min_rel 才算存在"""
    base = max((float(i) for _, i in peaks), default=0.0)
    res = []
    for d in diag:
        best = None
        for m, i in peaks:
            if abs(float(m) - d) <= tol and (best is None or float(i) > best[1]):
                best = (float(m), float(i))
        if best is None or base <= 0:
            res.append((None, None, False))
        else:
            rel = 100.0 * best[1] / base
            res.append((best[0], rel, rel >= min_rel))
    return res, base


class Redu:
    """本地 Pan-ReDU 缓存: 命中文件匹配 + 背景统计 (一次加载, 批量模式共用)"""

    def __init__(self, bg_mode="pos_lc_or_unknown"):
        if bg_mode not in REDU_BG_MODES:
            raise ValueError(f"--redu-bg 只能是 {REDU_BG_MODES}")
        t0 = time.time()
        header, rows = redu_load()
        self.cache_date = redu_cache_date()
        self.col = {c: i for i, c in enumerate(header)}
        c = self.col
        i_usi, i_ds, i_pol, i_chr = c["USI"], c["ATTRIBUTE_DatasetAccession"], \
            c["IonizationSourceAndPolarity"], c["ChromatographyAndPhase"]
        self.rows, self.lower, self.noext, self.base, self.datasets = [], {}, {}, {}, set()
        self.case_dup = {}  # 只差大小写的不同文件 (2026-10-12 缓存中 3 对, MetaboLights): 原样 USI -> 行号
        seen = set()
        for r in rows:
            if r[i_usi] in seen:  # 完全相同的重复 USI 只保留第一行
                continue
            seen.add(r[i_usi])
            lo = r[i_usi].lower()
            idx = len(self.rows)
            self.rows.append(r)
            if lo in self.lower:
                self.case_dup[r[i_usi]] = idx
                self.case_dup.setdefault(self.rows[self.lower[lo]][i_usi], self.lower[lo])
            else:
                self.lower[lo] = idx
            self.noext.setdefault(EXT_RE.sub("", lo), []).append(idx)
            coll, run, _, _ = parse_usi(lo)
            if coll is not None:
                self.datasets.add(coll)
                self.base.setdefault((coll, EXT_RE.sub("", run.rsplit("/", 1)[-1])), []).append(idx)
            self.datasets.add(r[i_ds].lower())
        del rows, seen
        self.bg_mode = bg_mode
        pol_ok = {"pos_lc": ("pos", "alt"), "pos_lc_or_unknown": ("pos", "alt", "unk"),
                  "neg_lc": ("neg", "alt"), "neg_lc_or_unknown": ("neg", "alt", "unk"),
                  "all": None}[bg_mode]
        self.in_bg = []
        for r in self.rows:
            if pol_ok is None:
                self.in_bg.append(True)
                continue
            gc = "gas chromatography" in r[i_chr].lower() or "electron ionization" in r[i_pol].lower()
            self.in_bg.append(polarity_of(r[i_pol]) in pol_ok and not gc)
        # 背景统计
        self.bg_files = collections.Counter()  # 数据集 -> 背景文件数
        self.bg_val_files = {k: collections.Counter() for k in SUMMARY_COLS}
        self.bg_val_ds = {k: collections.defaultdict(set) for k in SUMMARY_COLS}
        self.bg_ds_val = {k: collections.Counter() for k in SUMMARY_COLS}  # (数据集, 值) -> 文件数
        for idx, r in enumerate(self.rows):
            if not self.in_bg[idx]:
                continue
            ds = r[i_ds]
            self.bg_files[ds] += 1
            for k in SUMMARY_COLS:
                v = r[c[k]]
                self.bg_val_files[k][v] += 1
                self.bg_val_ds[k][v].add(ds)
                self.bg_ds_val[k][(ds, v)] += 1
        self.n_bg_files = sum(self.bg_files.values())
        self.n_bg_ds = len(self.bg_files)
        self.dsw_bg = {}
        for k in SUMMARY_COLS:
            share = collections.Counter()
            for (ds, v), n in self.bg_ds_val[k].items():
                share[v] += n / self.bg_files[ds]
            self.dsw_bg[k] = {v: s / self.n_bg_ds for v, s in share.items()} if self.n_bg_ds else {}
        print(f"[ReDU] 载入 {len(self.rows)} 个文件 (缓存 {self.cache_date}); 背景 {bg_mode}: "
              f"{self.n_bg_files} 个文件 / {self.n_bg_ds} 个数据集 ({time.time() - t0:.0f} s)")

    def match(self, usi):
        """命中 USI -> (ReDU 行号或 None, 匹配级别, 候选数)"""
        coll, run, _, _ = parse_usi(usi)
        if coll is None:
            return None, "unparsed", 0
        if coll.upper() == "GNPS" and run.upper() == "GNPS-LIBRARY":
            return None, "library", 0
        fu = f"mzspec:{coll}:{run}"
        if fu in self.case_dup:
            return self.case_dup[fu], "exact", 1
        lo = fu.lower()
        i = self.lower.get(lo)
        if i is not None:
            return i, ("exact" if self.rows[i][self.col["USI"]] == fu else "case_insensitive"), 1
        cands = self.noext.get(EXT_RE.sub("", lo))
        level = "no_extension"
        if not cands:
            cl, rl = coll.lower(), run.lower()
            cands = self.base.get((cl, EXT_RE.sub("", rl.rsplit("/", 1)[-1])))
            level = "basename_unique"
        if cands:
            if len(cands) == 1:
                return cands[0], level, 1

            def suffix_len(i):
                a = EXT_RE.sub("", parse_usi(self.rows[i][self.col["USI"]].lower())[1] or "")[::-1]
                b = EXT_RE.sub("", run.lower())[::-1]
                n = 0
                while n < min(len(a), len(b)) and a[n] == b[n]:
                    n += 1
                return n

            best = max(cands, key=lambda i: (suffix_len(i), -i))
            return best, ("no_extension_ambiguous" if level == "no_extension" else "basename_ambiguous"), len(cands)
        if coll.lower() in self.datasets:
            return None, "dataset_only", 0
        return None, "no_dataset", 0

    def value(self, idx, colname):
        return self.rows[idx][self.col[colname]]


def redu_annotate(R, rows):
    """给每条命中加 redu_* 列"""
    for r in rows:
        idx, lvl, nc = R.match(r["USI"])
        r["_redu_idx"] = idx
        r["redu_match"], r["redu_n_candidates"] = lvl, nc
        if idx is not None:
            for c in REDU_COLS:
                r["redu_" + c] = R.value(idx, c)
            r["redu_in_background"] = R.in_bg[idx]


def redu_summarise(R, hits):
    """hits: 一组命中 (已 redu_annotate). 返回 (汇总 dict, 汇总表行列表)"""
    lv = collections.Counter(r["redu_match"] for r in hits)
    files = {}
    for r in hits:
        files.setdefault(r["file_usi"], r)
    f_lv = collections.Counter(r["redu_match"] for r in files.values())
    ds_all = {r.get("Dataset") for r in hits}
    ds_matched = {r.get("Dataset") for r in hits if r["_redu_idx"] is not None}
    inb = [r for r in hits if r["_redu_idx"] is not None and R.in_bg[r["_redu_idx"]]]
    hit_files = {}  # ReDU 行号 -> 命中谱图数
    for r in inb:
        hit_files[r["_redu_idx"]] = hit_files.get(r["_redu_idx"], 0) + 1
    i_ds = R.col["ATTRIBUTE_DatasetAccession"]
    hit_ds_files = collections.Counter(R.rows[i][i_ds] for i in hit_files)
    n_hf, n_hd = len(hit_files), len(hit_ds_files)
    unmatched = lambda c: sum(n for k, n in c.items() if k not in MATCHED_LEVELS)  # noqa: E731
    summ = {"n_hits": len(hits), "match_levels_hits": dict(lv), "match_levels_files": dict(f_lv),
            "n_files": len(files), "n_datasets": len(ds_all), "n_datasets_with_matched_file": len(ds_matched),
            "unmatched_frac_hits": round(unmatched(lv) / len(hits), 4) if hits else None,
            "unmatched_frac_files": round(unmatched(f_lv) / len(files), 4) if files else None,
            "unmatched_frac_datasets": round(1 - len(ds_matched) / len(ds_all), 4) if ds_all else None,
            "n_hits_in_background": len(inb), "n_hit_files_in_background": n_hf,
            "n_hit_datasets_in_background": n_hd,
            "n_matched_hit_files_outside_background": len({r["_redu_idx"] for r in hits if r["_redu_idx"] is not None
                                                           and not R.in_bg[r["_redu_idx"]]}),
            "matched_hit_files_polarity": dict(collections.Counter(
                polarity_of(R.value(i, "IonizationSourceAndPolarity"))
                for i in {r["_redu_idx"] for r in hits if r["_redu_idx"] is not None})),
            "matched_hit_files_gc": sum(1 for i in {r["_redu_idx"] for r in hits if r["_redu_idx"] is not None}
                                        if "gas chromatography" in R.value(i, "ChromatographyAndPhase").lower()),
            "background": R.bg_mode, "n_bg_files": R.n_bg_files, "n_bg_datasets": R.n_bg_ds}
    table = []
    for k in SUMMARY_COLS:
        ci = R.col[k]
        v_hits, v_files, v_ds = collections.Counter(), collections.Counter(), collections.defaultdict(set)
        ds_v = collections.Counter()
        for i, n in hit_files.items():
            v, ds = R.rows[i][ci], R.rows[i][i_ds]
            v_hits[v] += n
            v_files[v] += 1
            v_ds[v].add(ds)
            ds_v[(ds, v)] += 1
        dsw = collections.Counter()
        for (ds, v), n in ds_v.items():
            dsw[v] += n / hit_ds_files[ds]
        for v in v_files:
            bgf, bgd = R.bg_val_files[k].get(v, 0), len(R.bg_val_ds[k].get(v, ()))
            hs = dsw[v] / n_hd if n_hd else None
            bs = R.dsw_bg[k].get(v)
            row = {"column": k, "value": v, "n_hits": v_hits[v], "n_hit_files": v_files[v],
                   "n_hit_datasets": len(v_ds[v]), "bg_files": bgf, "bg_datasets": bgd,
                   "file_rate": v_files[v] / bgf if bgf else None,
                   "file_enrich": (v_files[v] / n_hf) / (bgf / R.n_bg_files) if bgf and n_hf else None,
                   "ds_rate": len(v_ds[v]) / bgd if bgd else None,
                   "ds_enrich": (len(v_ds[v]) / bgd) / (n_hd / R.n_bg_ds) if bgd and n_hd else None,
                   "dsw_hit_share": hs, "dsw_bg_share": bs,
                   "dsw_enrich": hs / bs if hs is not None and bs else None}
            for kk in ("file_rate", "file_enrich", "ds_rate", "ds_enrich", "dsw_hit_share", "dsw_bg_share",
                       "dsw_enrich"):
                if row[kk] is not None:
                    row[kk] = round(row[kk], 6)
            table.append(row)
    table.sort(key=lambda r: (SUMMARY_COLS.index(r["column"]), -r["n_hit_datasets"], -r["n_hit_files"]))
    return summ, table


MATCHED_LEVELS = {"exact", "case_insensitive", "no_extension", "no_extension_ambiguous", "basename_unique",
                  "basename_ambiguous"}


def redu_by_dataset(R, hits):
    """每个数据集一行: 命中数 / 命中文件数 / 背景文件数 / 主要元数据值"""
    by = collections.defaultdict(list)
    for r in hits:
        by[r.get("Dataset") or ""].append(r)
    out = []
    for ds, rs in by.items():
        files = {r["file_usi"] for r in rs}
        mf = {r["_redu_idx"] for r in rs if r["_redu_idx"] is not None}
        bf = {i for i in mf if R.in_bg[i]}
        acc = R.value(next(iter(mf)), "ATTRIBUTE_DatasetAccession") if mf else ds
        row = {"dataset": ds, "redu_source": R.value(next(iter(mf)), "DataSource") if mf else "",
               "n_hits": len(rs), "n_hit_files": len(files), "n_hit_files_matched": len(mf),
               "n_hit_files_in_background": len(bf), "bg_files_in_dataset": R.bg_files.get(acc, 0),
               "frac_bg_files_hit": round(len(bf) / R.bg_files[acc], 4) if R.bg_files.get(acc) else None,
               "max_cosine": max(r.get("Cosine") or 0 for r in rs),
               "max_matching_peaks": max(r.get("Matching Peaks") or 0 for r in rs),
               "n_verified": sum(1 for r in rs if r.get("verify_status") == "ok"),
               "n_diag_pass": sum(1 for r in rs if r.get("diag_pass") is True),
               "match_levels": ";".join(f"{k}:{n}" for k, n in collections.Counter(
                   r["redu_match"] for r in rs).most_common())}
        for k in SUMMARY_COLS:
            c = collections.Counter(R.value(i, k) for i in mf)
            row["top_" + k] = "; ".join(f"{v} ({n})" for v, n in c.most_common(2))
        out.append(row)
    out.sort(key=lambda r: (-r["n_hit_files"], r["dataset"]))
    return out


HIT_TYPES = {"Cosine": float, "Matching Peaks": int, "Delta Mass": float, "Unit Delta Mass": int, "Charge": int,
             "Query Scan": int, "Index UnitPM": int, "Index IdxInUnitPM": int}


def load_hits_tsv(path):
    """读取 --keep-all 写出的 <out>.hits.tsv 与同名 JSON -> {"prev": JSON, "raw": {结果键: [API 命中 dict]}}"""
    pj = re.sub(r"\.hits\.tsv$", ".json", path)
    if pj == path or not os.path.exists(pj):
        raise ValueError(f"--from-hits: 找不到同名 JSON {pj}")
    with open(pj) as f:
        prev = json.load(f)
    if "query_params" not in prev:
        raise ValueError(f"--from-hits: {pj} 没有 query_params (不是 --keep-all 的输出)")
    raw = collections.defaultdict(list)
    with open(path, newline="", encoding="utf-8") as f:
        rd = csv.DictReader(f, delimiter="\t")
        api_cols = [c for c in rd.fieldnames if c in API_FIELDS or c[:1].isupper()]
        for r in rd:
            x = {}
            for c in api_cols:
                v = r[c]
                if v == "":
                    continue
                x[c] = HIT_TYPES[c](v) if c in HIT_TYPES else v
            raw[r["search"]].append(x)
    return {"prev": prev, "raw": raw}


# ---------- 单个查询 ----------
def load_query_spectrum(spec, fetcher_factory):
    """spec 中 peaks / mgf / usi 三选一 -> (peaks, precursor, 来源说明)"""
    precursor = spec.get("precursor")
    src = None
    if spec.get("mgf"):
        pm, pk = None, []
        for line in open(spec["mgf"]):
            line = line.strip()
            if line.startswith("PEPMASS"):
                pm = float(line.split("=")[1].split()[0])
            elif line and line[0].isdigit():
                a, b = line.split()[:2]
                pk.append([float(a), float(b)])
        precursor = precursor or pm
        peaks, src = pk, {"type": "mgf", "path": spec["mgf"]}
    elif spec.get("usi"):
        data, st = fetcher_factory().fetch(spec["usi"])
        if data is None:
            raise RuntimeError(f"查询 USI 取谱失败 ({st}): {spec['usi']}")
        peaks = [[float(m), float(i)] for m, i in data["peaks"]]
        src = {"type": "usi", "usi": spec["usi"], "usi_precursor_mz": data.get("precursor_mz"),
               "usi_charge": data.get("precursor_charge")}
        precursor = precursor or data.get("precursor_mz")
    else:
        pk = spec.get("peaks")
        peaks = json.loads(pk) if isinstance(pk, str) else pk
        src = {"type": "peaks"}
    if not peaks or not precursor:
        raise ValueError("需要提供谱图: peaks 与 precursor, 或 mgf, 或 usi")
    peaks = [[float(m), float(i)] for m, i in peaks]
    n0 = len(peaks)
    qmr = spec.get("query_min_rel")
    if qmr:
        base = max(i for _, i in peaks)
        peaks = [p for p in peaks if p[1] >= base * float(qmr) / 100.0]
    src["n_peaks_input"], src["n_peaks_used"] = n0, len(peaks)
    return peaks, float(precursor), src


def run_query(spec, out_json, ctx, new_mode=True):
    """执行一个查询. spec: 合并后的参数 dict. 返回结果 dict (同时写 JSON 与表格)"""
    o = {k: spec.get(k, OPT_DEFAULTS.get(k)) for k in OPT_DEFAULTS}
    name = spec.get("name") or spec.get("inchikey") or "unknown"
    out_dir = os.path.dirname(os.path.abspath(out_json))
    os.makedirs(out_dir, exist_ok=True)
    stem = out_json[:-5] if out_json.endswith(".json") else out_json
    cache_dir = ctx.get("usi_cache") or os.path.join(out_dir, "usi_cache")

    def fetcher():
        if ctx.get("fetcher") is None or ctx["fetcher"].cache_dir != cache_dir:
            ctx["fetcher"] = UsiFetcher(cache_dir, o["usi_rate"])
        return ctx["fetcher"]

    replay = None
    if spec.get("from_hits"):  # 重放: 不调用 fastMASST, 用上次 --keep-all 保存的全部命中重新分析
        replay = load_hits_tsv(spec["from_hits"])
        qp = replay["prev"]["query_params"]
        peaks, precursor = qp["peaks"], qp["precursor"]
        src = {"type": "from_hits", "path": spec["from_hits"], "orig_spectrum_source": qp.get("spectrum_source"),
               "n_peaks_used": len(peaks)}
        if o["cosine"] < qp["cosine"] - 1e-9:
            raise ValueError(f"--from-hits: cosine {o['cosine']} 低于原检索的 {qp['cosine']}, 无法重放")
        for k in ("pm_tol", "frag_tol"):
            o[k] = qp[k]  # 容差无法事后改变, 沿用原检索
        name = spec.get("name") if spec.get("name") not in (None, "unknown") else replay["prev"].get("name")
    else:
        peaks, precursor, src = load_query_spectrum(spec, fetcher)
    diag = parse_mz_list(o["diagnostic"])
    req = parse_mz_list(o["diag_required"])
    for m in req:
        if not any(abs(m - d) < 1e-6 for d in diag):
            diag.append(m)
    k_need = parse_diag_rule(o["diag_rule"], len(diag)) if diag else None
    keep_all = bool(o["keep_all"] or diag or o["redu_join"])
    searches = resolve_libraries(o["library"]) if replay is None else \
        [(x["key"], x["library"]) for x in replay["prev"]["query_params"]["searches"]]

    results = {"name": name, "inchikey": spec.get("inchikey") if replay is None else replay["prev"].get("inchikey")}
    run_info = {}
    raw_by_key = {}
    both = [k for k, _ in searches] == ["masst_repository", "masst_library"]
    for n, (key, lib) in enumerate(searches, 1):
        if replay is not None:
            prev_m = replay["prev"][key]
            prev_ri = replay["prev"].get("run_info", {}).get(key, {})
            if "error" in prev_m:
                results[key], raw = prev_m, None
            else:
                raw = [x for x in replay["raw"].get(key, []) if (x.get("Cosine") or 0) >= o["cosine"] - 1e-9]
                results[key] = summarise_raw(raw, lib, o["min_peaks"], prev_m.get("task_id"))
            raw_by_key[key] = (lib, raw)
            run_info[key] = {"library": lib, "task_id": prev_m.get("task_id") or prev_ri.get("task_id"),
                             "replayed_from": spec["from_hits"], "elapsed_s": 0, "n_poll_errors": 0,
                             "n_datasets": len({x.get("Dataset") for x in raw or []
                                                if (x.get("Matching Peaks") or 0) >= o["min_peaks"]})}
            print(f"[重放] {key}: {len(raw or [])} 条命中 (来自 {spec['from_hits']}, 未调用 fastMASST)")
            continue
        if both and key == "masst_repository":
            print(f"[1/2] fastMASST 全库检索 ({len(peaks)} peaks, 可能需要 1-2 分钟) ...")
        elif both:
            print("[2/2] GNPS 参考谱库检索 ...")
        else:
            print(f"[{n}/{len(searches)}] fastMASST 检索 library={lib} ({len(peaks)} peaks) ...")
        if n > 1 and ctx.get("pause"):
            time.sleep(ctx["pause"])
        r, info = fasst_search(peaks, precursor, lib, o["cosine"], o["min_peaks"], o["pm_tol"], o["frag_tol"],
                               max_wait_s=o["max_wait"], return_raw=True)
        results[key] = r
        raw_by_key[key] = (lib, info["raw"])
        run_info[key] = {"library": lib, "task_id": info["task_id"], "elapsed_s": info["elapsed_s"],
                         "n_poll_errors": info["n_poll_errors"],
                         "n_datasets": len({x.get("Dataset") for x in info["raw"] or []
                                            if (x.get("Matching Peaks") or 0) >= o["min_peaks"]})}
    print_summary(name, results)

    if new_mode:
        results["query_params"] = {
            "tool_version": TOOL_VERSION, "run_time": time.strftime("%Y-%m-%d %H:%M:%S"),
            "precursor": precursor, "peaks": peaks, "spectrum_source": src, "notes": spec.get("notes"),
            "searches": [{"key": k, "library": l} for k, l in searches],
            **{k: o[k] for k in ("cosine", "min_peaks", "pm_tol", "frag_tol", "max_wait", "query_min_rel")},
            "analog": "No"}
        results["run_info"] = run_info

    if keep_all:
        cols_api = list(API_FIELDS)
        all_rows = {}
        for key, (lib, raw) in raw_by_key.items():
            rows = []
            for x in raw or []:
                for f in x:
                    if f not in cols_api:
                        cols_api.append(f)  # API 将来新增的字段也保留
                coll, run, it, iv = parse_usi(x.get("USI"))
                dm = x.get("Delta Mass")
                r = {"search": key, "library": lib, "task_id": run_info[key]["task_id"], **x,
                     "file_usi": file_usi_of(x.get("USI")), "usi_collection": coll, "index_type": it,
                     "scan": iv,
                     "est_hit_precursor_mz": round(precursor - dm, 4) if isinstance(dm, (int, float)) else None,
                     "pass_min_peaks": (x.get("Matching Peaks") or 0) >= o["min_peaks"]}
                rows.append(r)
            all_rows[key] = rows
        outputs = {"hits_tsv": stem + ".hits.tsv"}
        extra_cols = []

        # ---- 诊断离子核对 ----
        if diag:
            labels = [mz_label(d) for d in diag]
            f = fetcher()
            ver = {}
            not_ver = []
            for key, rows in all_rows.items():
                cands = [r for r in rows if r["pass_min_peaks"]]
                resolvable = [r for r in cands if r.get("scan") != INVALID_SCAN]
                sel = select_for_verification(resolvable, o["verify_max"], o["verify_seed"])
                sel_ids = {id(r) for r in sel}
                for r in rows:
                    if not r["pass_min_peaks"]:
                        r["verify_status"] = "below_min_peaks"
                    elif r.get("scan") == INVALID_SCAN:
                        r["verify_status"] = "invalid_scan"
                    elif id(r) not in sel_ids:
                        r["verify_status"] = "not_selected"
                if sel:
                    print(f"[核对] {key}: 取回 {len(sel)}/{len(cands)} 张命中谱图 (缓存 {cache_dir}) ...")
                ion_counts = collections.Counter()
                todo = list(sel)
                for npass in (1, 2):  # 第 2 轮: 暂停 20 s 后重试第 1 轮失败的 (非 4xx)
                    if npass == 2:
                        todo = [r for r in sel if r["verify_status"].startswith("fetch_failed")
                                and not re.search(r"HTTP 4(?!29)", r["verify_status"])]
                        if not todo:
                            break
                        print(f"[核对]   第 2 轮重试 {len(todo)} 个失败的谱图 (20 s 后) ...")
                        time.sleep(20)
                    for j, r in enumerate(todo, 1):
                        data, st = f.fetch(r["USI"])
                        r["verify_status"] = st
                        if data is not None:
                            r["_spec"] = data
                        if j % 25 == 0 or j == len(todo):
                            print(f"[核对]   {j}/{len(todo)}")
                for r in sel:
                    data = r.pop("_spec", None)
                    if data is not None:
                        pk = data.get("peaks") or []
                        ev, base = diag_eval(pk, diag, o["diag_tol"], o["diag_min_rel"])
                        r["hit_precursor_mz"] = data.get("precursor_mz")
                        r["hit_n_peaks"] = len(pk)
                        r["hit_base_peak_mz"] = round(max(pk, key=lambda p: p[1])[0], 4) if pk else None
                        npres = 0
                        req_ok = True
                        for lab, d, (mz, rel, pres) in zip(labels, diag, ev):
                            r[f"diag_{lab}_rel"] = round(rel, 2) if rel is not None else None
                            r[f"diag_{lab}_mz"] = round(mz, 4) if mz is not None else None
                            npres += pres
                            ion_counts[lab] += pres
                            if not pres and any(abs(d - q) < 1e-6 for q in req):
                                req_ok = False
                        r["diag_n_present"] = npres
                        r["diag_pass"] = bool(npres >= k_need and req_ok)
                for r in rows:
                    if r["pass_min_peaks"] and r["verify_status"] != "ok":
                        not_ver.append({"search": key, "USI": r["USI"], "Dataset": r.get("Dataset"),
                                        "Cosine": r.get("Cosine"), "Matching Peaks": r.get("Matching Peaks"),
                                        "reason": r["verify_status"]})
                okr = [r for r in cands if r.get("verify_status") == "ok"]
                pas = [r for r in okr if r.get("diag_pass")]
                ver[key] = {"n_candidates": len(cands), "n_candidate_datasets": len({r.get("Dataset") for r in cands}),
                            "n_invalid_scan": len(cands) - len(resolvable),
                            "n_selected": len(sel), "n_not_selected": len(resolvable) - len(sel),
                            "n_verified_ok": len(okr), "n_fetch_failed": len(sel) - len(okr),
                            "n_verified_datasets": len({r.get("Dataset") for r in okr}),
                            "n_diag_pass": len(pas), "n_diag_pass_datasets": len({r.get("Dataset") for r in pas}),
                            "diag_pass_rate": round(len(pas) / len(okr), 4) if okr else None,
                            "ion_present_counts": {lab: ion_counts[lab] for lab in labels}}
                print(f"[核对] {key}: 成功 {len(okr)}, 失败 {len(sel) - len(okr)}, 未抽中 {len(resolvable) - len(sel)}, "
                      f"无效扫描号 {len(cands) - len(resolvable)}; "
                      f"诊断通过 {len(pas)} ({ver[key]['n_diag_pass_datasets']} 个数据集)")
            extra_cols += ["verify_status", "hit_precursor_mz", "hit_n_peaks", "hit_base_peak_mz"]
            for lab in labels:
                extra_cols += [f"diag_{lab}_rel", f"diag_{lab}_mz"]
            extra_cols += ["diag_n_present", "diag_pass"]
            outputs["not_verified_tsv"] = stem + ".not_verified.tsv"
            write_tsv(outputs["not_verified_tsv"], not_ver,
                      ["search", "USI", "Dataset", "Cosine", "Matching Peaks", "reason"])
            results["verification"] = {"diagnostic": labels, "required": [mz_label(m) for m in req],
                                       "diag_tol": o["diag_tol"], "diag_min_rel": o["diag_min_rel"],
                                       "diag_rule": str(o["diag_rule"]), "k_required": k_need,
                                       "verify_max": o["verify_max"], "verify_seed": o["verify_seed"],
                                       "usi_cache": cache_dir, "fetch_stats": dict(f.stats), **ver}

        # ---- Pan-ReDU 关联 ----
        if o["redu_join"]:
            if ctx.get("redu") is None or ctx["redu"].bg_mode != o["redu_bg"]:
                factory = ctx.get("redu_factory")
                ctx["redu"] = factory(o["redu_bg"], all_rows) if factory else Redu(o["redu_bg"])
            R = ctx["redu"]
            rsum, rtab, rds = {}, [], []
            for key, rows in all_rows.items():
                if raw_by_key[key][0] == REF_LIB:
                    continue
                redu_annotate(R, rows)
                hits = [r for r in rows if r["pass_min_peaks"]]
                sets = [("all", hits)]
                if diag:
                    sets += [("verified", [r for r in hits if r.get("verify_status") == "ok"]),
                             ("diag_pass", [r for r in hits if r.get("diag_pass") is True])]
                rsum[key] = {"redu_cache_date": R.cache_date}
                for sname, hs in sets:
                    s, t = redu_summarise(R, hs)
                    rsum[key][sname] = s
                    rtab += [{"search": key, "hit_set": sname, **row} for row in t]
                rds += [{"search": key, **row} for row in redu_by_dataset(R, hits)]
                a = rsum[key]["all"]
                print(f"[ReDU] {key}: {a['n_hits']} 命中 / {a['n_files']} 文件 / {a['n_datasets']} 数据集; "
                      f"未匹配比例 命中 {a['unmatched_frac_hits']} 文件 {a['unmatched_frac_files']} "
                      f"数据集 {a['unmatched_frac_datasets']}; 背景内 {a['n_hit_files_in_background']} 文件 / "
                      f"{a['n_hit_datasets_in_background']} 数据集")
                for col in ("SampleType", "NCBITaxonomy", "UBERONBodyPartName"):
                    top = [r for r in rtab if r["search"] == key and r["hit_set"] == "all" and r["column"] == col][:5]
                    print(f"  {col}: " + ", ".join(f"{r['value']} ({r['n_hit_datasets']} 数据集/{r['n_hit_files']} 文件)"
                                                   for r in top))
            extra_cols += ["redu_match", "redu_n_candidates"] + ["redu_" + c for c in REDU_COLS] + \
                          ["redu_in_background"]
            outputs["redu_summary_tsv"] = stem + ".redu_summary.tsv"
            outputs["redu_by_dataset_tsv"] = stem + ".redu_by_dataset.tsv"
            write_tsv(outputs["redu_summary_tsv"], rtab,
                      ["search", "hit_set", "column", "value", "n_hits", "n_hit_files", "n_hit_datasets", "bg_files",
                       "bg_datasets", "file_rate", "file_enrich", "ds_rate", "ds_enrich", "dsw_hit_share",
                       "dsw_bg_share", "dsw_enrich"])
            ds_cols = ["search", "dataset", "redu_source", "n_hits", "n_hit_files", "n_hit_files_matched",
                       "n_hit_files_in_background", "bg_files_in_dataset", "frac_bg_files_hit", "max_cosine",
                       "max_matching_peaks", "n_verified", "n_diag_pass", "match_levels"] + \
                      ["top_" + k for k in SUMMARY_COLS]
            write_tsv(outputs["redu_by_dataset_tsv"], rds, ds_cols)
            results["redu"] = rsum

        cols = ["search", "library", "task_id"] + cols_api + \
               ["file_usi", "usi_collection", "index_type", "scan", "est_hit_precursor_mz", "pass_min_peaks"] + \
               extra_cols
        write_tsv(outputs["hits_tsv"], [r for rows in all_rows.values() for r in rows], cols)
        results["outputs"] = outputs
        print("[输出] " + ", ".join(outputs.values()))

    with open(out_json, "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\n结果已保存: {out_json}")
    return results


# ---------- MassQL 本地检索 (MS2 常用子集) ----------
# 完整 MassQL 由官方 massql 包/GNPS2 MassQL 应用执行; 这里实现常用 MS2 子集, 不引入第三方依赖,
# 便于直接筛本地 mzML/mzXML/MGF。遇到不支持的语法会明确报错, 不会静默给出错误结果。
ELEMENT_MASS = {
    "C": 12.0, "H": 1.00782503223, "N": 14.00307400443, "O": 15.99491461957,
    "P": 30.97376199842, "S": 31.9720711744, "Na": 22.9897692820, "K": 38.9637064864,
    "Cl": 34.968852682, "Br": 78.9183376, "F": 18.99840316273, "Si": 27.97692653465,
    "B": 11.00930536, "Se": 79.9165218, "I": 126.9044719,
}
MASSQL_MASS_KEYS = {"MS2PROD", "MS2NL", "MS2PREC"}
MASSQL_OTHER_KEYS = {"RTMIN", "RTMAX", "SCANMIN", "SCANMAX", "POLARITY", "CHARGE"}
MASSQL_QUALS = {"TOLERANCEMZ", "TOLERANCEPPM", "INTENSITYPERCENT", "INTENSITYVALUE",
                "EXCLUDED", "CARDINALITY"}


def _massql_formula_mass(formula):
    """'H2O' / 'C6H12O6' -> 单同位素质量; 仅支持简单分子式"""
    s = (formula or "").strip()
    if not s:
        raise ValueError("formula() 为空")
    total, pos = 0.0, 0
    for m in re.finditer(r"([A-Z][a-z]?)(\d*)", s):
        if m.start() != pos:
            raise ValueError(f"无法解析分子式: {formula}")
        el, n = m.group(1), m.group(2)
        if el not in ELEMENT_MASS:
            raise ValueError(f"分子式中有未知元素: {el}")
        total += ELEMENT_MASS[el] * (int(n) if n else 1)
        pos = m.end()
    if pos != len(s):
        raise ValueError(f"无法解析分子式: {formula}")
    return total


def _massql_eval_expr(expr):
    """质量表达式: 数字、+ - * /、括号、formula(...), 如 100+14 / formula(H2O)"""
    expr = (expr or "").strip()
    if not expr:
        raise ValueError("空的质量表达式")
    if re.search(r"\bX\b", expr):
        raise ValueError("暂不支持 MassQL 变量 X")
    try:
        tree = ast.parse(expr, mode="eval")
    except SyntaxError as e:
        raise ValueError(f"无法解析质量表达式: {expr}") from e

    def ev(node):
        if isinstance(node, ast.Expression):
            return ev(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return float(node.value)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            v = ev(node.operand)
            return v if isinstance(node.op, ast.UAdd) else -v
        if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div)):
            a, b = ev(node.left), ev(node.right)
            if isinstance(node.op, ast.Add):
                return a + b
            if isinstance(node.op, ast.Sub):
                return a - b
            if isinstance(node.op, ast.Mult):
                return a * b
            if b == 0:
                raise ValueError(f"质量表达式除以 0: {expr}")
            return a / b
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "formula":
            if len(node.args) != 1:
                raise ValueError("formula() 需要一个参数, 如 formula(H2O)")
            arg = node.args[0]
            if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                formula = arg.value
            elif isinstance(arg, ast.Name):
                formula = arg.id
            else:
                raise ValueError("formula() 需要一个分子式, 如 formula(H2O)")
            return _massql_formula_mass(formula)
        raise ValueError(f"不支持的质量表达式: {expr}")

    return ev(tree)


def _massql_split_top(text, sep):
    """按顶层分隔符切分, 忽略括号内的 AND/OR/冒号"""
    out, buf, depth, i = [], [], 0, 0
    n = len(text)
    while i < n:
        c = text[i]
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        if depth == 0 and text[i:i + len(sep)].upper() == sep.upper():
            if sep.isalpha():
                before = text[i - 1] if i else " "
                after = text[i + len(sep)] if i + len(sep) < n else " "
                if before.isalnum() or after.isalnum() or before == "_" or after == "_":
                    buf.append(c)
                    i += 1
                    continue
            out.append("".join(buf).strip())
            buf = []
            i += len(sep)
            continue
        buf.append(c)
        i += 1
    out.append("".join(buf).strip())
    return [x for x in out if x]


def _massql_cardinality(value, n):
    v = "" if value is True else str(value or "").strip().lower()
    if not v:
        return (1, 1) if n == 1 else (1, n)
    m = re.fullmatch(r"range\s*\(\s*min\s*=\s*(\d+)(?:\s*,\s*max\s*=\s*(\d+))?\s*\)", v)
    if m:
        lo, hi = int(m.group(1)), int(m.group(2)) if m.group(2) else n
    elif re.fullmatch(r"\d+", v):
        lo = hi = int(v)
    else:
        raise ValueError(f"CARDINALITY 无法解析: {value}")
    if not 0 <= lo <= hi <= n:
        raise ValueError(f"CARDINALITY 超出范围 (候选数 {n}): {value}")
    return lo, hi


def _massql_parse_condition(text):
    if "=" not in text:
        raise ValueError(f"MassQL 条件缺少 '=': {text}")
    key, rest = text.split("=", 1)
    key = key.strip().upper()
    if key == "MS2MZ":
        key = "MS2PROD"
    if key in ("MS1MZ", "MOBILITY"):
        raise ValueError(f"暂不支持条件 {key}")
    if key not in MASSQL_MASS_KEYS | MASSQL_OTHER_KEYS:
        raise ValueError(f"暂不支持条件 {key}")

    pieces = _massql_split_top(rest, ":")
    raw_value = pieces[0].strip()
    quals = {}
    for p in pieces[1:]:
        p = p.strip()
        if not p:
            continue
        if "=" in p:
            qk, qv = p.split("=", 1)
            qk, qv = qk.strip().upper(), qv.strip()
        else:
            qk, qv = p.upper(), True
        if qk not in MASSQL_QUALS:
            raise ValueError(f"暂不支持限定符 {qk}")
        quals[qk] = qv

    if key in MASSQL_MASS_KEYS:
        if raw_value.upper() == "ANY":
            raise ValueError("暂不支持 ANY")
        if raw_value.startswith("("):
            if not raw_value.endswith(")"):
                raise ValueError(f"括号不匹配: {text}")
            vals = [_massql_eval_expr(x) for x in _massql_split_top(raw_value[1:-1], "OR")]
        else:
            vals = [_massql_eval_expr(raw_value)]
    elif key == "POLARITY":
        p = raw_value.lower()
        if "pos" in p:
            vals = ["Positive"]
        elif "neg" in p:
            vals = ["Negative"]
        else:
            raise ValueError(f"POLARITY 只能是 Positive/Negative: {raw_value}")
    elif key == "CHARGE":
        vals = [int(float(raw_value))]
    else:  # RTMIN/RTMAX/SCANMIN/SCANMAX
        vals = [float(raw_value)]

    tol_da = float(quals.get("TOLERANCEMZ", 0.1))
    tol_ppm = float(quals["TOLERANCEPPM"]) if "TOLERANCEPPM" in quals else None
    int_pct = float(quals["INTENSITYPERCENT"]) if "INTENSITYPERCENT" in quals else None
    int_val = float(quals["INTENSITYVALUE"]) if "INTENSITYVALUE" in quals else None
    if tol_da < 0 or (tol_ppm is not None and tol_ppm < 0):
        raise ValueError(f"质量容差不能为负: {text}")
    if key not in MASSQL_MASS_KEYS and any(q in quals for q in
                                           ("TOLERANCEMZ", "TOLERANCEPPM", "INTENSITYPERCENT",
                                            "INTENSITYVALUE", "CARDINALITY")):
        raise ValueError(f"{key} 不支持质量/强度限定符: {text}")
    if key == "MS2PREC" and any(q in quals for q in ("INTENSITYPERCENT", "INTENSITYVALUE")):
        raise ValueError("MS2PREC 不支持碎片强度限定符")

    lo, hi = _massql_cardinality(quals.get("CARDINALITY"), len(vals)) \
        if key in MASSQL_MASS_KEYS else (1, 1)
    return {"key": key, "values": vals, "tol_da": tol_da, "tol_ppm": tol_ppm,
            "int_pct": int_pct, "int_val": int_val, "min": lo, "max": hi,
            "excluded": "EXCLUDED" in quals, "raw": text.strip()}


def parse_massql(query):
    """解析常用 MS2 MassQL 子集 -> {query, conditions, ...}"""
    q = " ".join((query or "").split()).strip().rstrip(";")
    if not q:
        raise ValueError("--massql 为空")
    if not re.match(r"^QUERY\b", q, re.I):
        raise ValueError("MassQL 查询必须以 QUERY 开头")
    if re.search(r"\bFILTER\b", q, re.I):
        raise ValueError("暂不支持 FILTER 子句")
    m = re.match(r"^QUERY\s+(?:([A-Za-z_]\w*)\s*\(\s*(MS[12]DATA)\s*\)|(MS[12]DATA))"
                 r"\s*(?:WHERE\s+(.+))?$", q, re.I)
    if not m:
        raise ValueError("无法解析 MassQL 查询; 例: QUERY scaninfo(MS2DATA) WHERE MS2PROD=100:TOLERANCEMZ=0.1")
    func, datatype, where = m.group(1), (m.group(2) or m.group(3)).upper(), m.group(4)
    if datatype != "MS2DATA":
        raise ValueError("当前只支持 MS2DATA")
    if func and func.lower() != "scaninfo":
        raise ValueError(f"暂不支持查询函数 {func}()")
    if not where:
        raise ValueError("MassQL 查询需要 WHERE 条件")
    conditions = [_massql_parse_condition(p) for p in _massql_split_top(where, "AND")]
    if not conditions:
        raise ValueError("MassQL 查询没有有效条件")
    return {"query": q, "datatype": datatype, "function": func, "conditions": conditions}


def _massql_tol(cond, target):
    if cond["tol_ppm"] is not None:
        return abs(target) * cond["tol_ppm"] / 1e6
    return cond["tol_da"]


def _massql_peak_ok(cond, intensity, base):
    if cond["int_pct"] is not None:
        if base <= 0 or intensity / base * 100.0 < cond["int_pct"]:
            return False
    if cond["int_val"] is not None and intensity < cond["int_val"]:
        return False
    return True


def _massql_match_condition(cond, spec):
    key = cond["key"]
    # Missing metadata is unknown, not evidence for an excluded condition.
    required = {"RTMIN": "rt_min", "RTMAX": "rt_min", "SCANMIN": "scan", "SCANMAX": "scan",
                "POLARITY": "polarity", "CHARGE": "charge", "MS2PREC": "precursor_mz", "MS2NL": "precursor_mz"}
    if key in required and spec.get(required[key]) is None:
        return False
    if key == "RTMIN":
        ok = spec.get("rt_min") is not None and spec["rt_min"] >= cond["values"][0]
    elif key == "RTMAX":
        ok = spec.get("rt_min") is not None and spec["rt_min"] <= cond["values"][0]
    elif key in ("SCANMIN", "SCANMAX"):
        m = re.search(r"\d+", str(spec.get("scan") or ""))
        if not m:
            ok = False
        else:
            scan = int(m.group())
            ok = scan >= cond["values"][0] if key == "SCANMIN" else scan <= cond["values"][0]
    elif key == "POLARITY":
        ok = (spec.get("polarity") or "").lower() == cond["values"][0].lower()
    elif key == "CHARGE":
        try:
            ok = int(spec.get("charge")) == int(cond["values"][0])
        except (TypeError, ValueError):
            ok = False
    else:
        peaks = spec.get("peaks") or []
        base = max((i for _, i in peaks), default=0)
        prec = spec.get("precursor_mz")
        n = 0
        for target in cond["values"]:
            tol = _massql_tol(cond, target)
            hit = False
            if key == "MS2PREC":
                hit = prec is not None and abs(prec - target) <= tol
            elif key == "MS2NL":
                if prec is not None:
                    for mz, inten in peaks:
                        if abs((prec - mz) - target) <= tol and _massql_peak_ok(cond, inten, base):
                            hit = True
                            break
            else:  # MS2PROD
                for mz, inten in peaks:
                    if abs(mz - target) <= tol and _massql_peak_ok(cond, inten, base):
                        hit = True
                        break
            if hit:
                n += 1
        ok = cond["min"] <= n <= cond["max"]
    return (not ok) if cond["excluded"] else ok


def _massql_eval_query(q, spec):
    matched = [c["raw"] for c in q["conditions"] if _massql_match_condition(c, spec)]
    return len(matched) == len(q["conditions"]), matched


def _massql_mgf_spec(cur, path, idx):
    peaks = cur["peaks"]
    if not peaks:
        return None
    meta = cur["meta"]
    pm = None
    if "PEPMASS" in meta:
        try:
            pm = float(meta["PEPMASS"].split()[0])
        except (ValueError, IndexError):
            pass
    rt = None
    for k in ("RTINSECONDS", "RTINMINUTES", "RT"):
        if k in meta:
            try:
                v = float(meta[k].split()[0])
                rt = v / 60.0 if k == "RTINSECONDS" else v
            except (ValueError, IndexError):
                pass
            break
    charge = None
    if "CHARGE" in meta:
        m = re.search(r"-?\d+", meta["CHARGE"])
        if m:
            charge = int(m.group())
            if meta["CHARGE"].strip().endswith("-"):
                charge = -abs(charge)
    scan = None
    for k in ("SCANS", "SCAN", "SCAN_NUMBER"):
        if k in meta:
            m = re.search(r"\d+", meta[k])
            if m:
                scan = m.group()
            break
    if scan is None:
        scan = str(idx)
    pol = None
    for k in ("IONMODE", "POLARITY", "IONIZATION"):
        if k in meta:
            v = meta[k].lower()
            if "neg" in v:
                pol = "Negative"
            elif "pos" in v:
                pol = "Positive"
            break
    title = meta.get("TITLE") or meta.get("NAME") or f"{os.path.basename(path)} scan {scan}"
    return {"file": path, "scan": scan, "precursor_mz": pm, "rt_min": rt, "polarity": pol,
            "charge": charge, "peaks": peaks, "title": title}


def _massql_load_mgf(path):
    out, cur = [], None
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            u = s.upper()
            if u == "BEGIN IONS":
                cur = {"meta": {}, "peaks": []}
                continue
            if u == "END IONS":
                if cur is not None:
                    spec = _massql_mgf_spec(cur, path, len(out) + 1)
                    if spec:
                        out.append(spec)
                cur = None
                continue
            if cur is None:
                continue
            if "=" in s:
                k, v = s.split("=", 1)
                cur["meta"][k.strip().upper()] = v.strip()
            else:
                a = s.split()
                if len(a) >= 2:
                    try:
                        cur["peaks"].append((float(a[0]), float(a[1])))
                    except ValueError:
                        pass
    return out


def _massql_load_mzml(path):
    try:
        from pyopenms import MSExperiment, MzMLFile, MzXMLFile
    except ImportError as e:
        raise RuntimeError("读取 mzML/mzXML 需要 pyopenms: pip install pyopenms") from e
    exp = MSExperiment()
    loader = MzXMLFile() if path.lower().endswith(".mzxml") else MzMLFile()
    loader.load(path, exp)
    out = []
    for idx, spec in enumerate(exp):
        if spec.getMSLevel() != 2:
            continue
        mz, inten = spec.get_peaks()
        if len(mz) == 0:
            continue
        pm, charge, pol = None, None, None
        precs = spec.getPrecursors()
        if precs:
            pm = precs[0].getMZ()
            charge = precs[0].getCharge()
        try:
            p = spec.getInstrumentSettings().getPolarity()
            if p in (1, 2):
                pol = "Positive" if p == 1 else "Negative"
        except Exception:
            pass
        native = spec.getNativeID() or ""
        scan_match = re.search(r"(?:^|\s)scan=(\d+)(?:\s|$)", native, re.I)
        out.append({"file": path, "scan": scan_match.group(1) if scan_match else None, "precursor_mz": pm,
                     "rt_min": spec.getRT() / 60.0, "polarity": pol, "charge": charge,
                     "peaks": [(float(a), float(b)) for a, b in zip(mz, inten)],
                     "title": native or f"{os.path.basename(path)} spectrum {idx + 1}"})
    return out


def _massql_collect_files(paths):
    exts = (".mzml", ".mzxml", ".mgf")
    files, missing = [], []
    for p in paths or []:
        if os.path.isdir(p):
            for root, dirs, names in os.walk(p):
                dirs[:] = sorted(d for d in dirs if not d.startswith("."))
                for n in sorted(names):
                    if n.lower().endswith(exts):
                        files.append(os.path.join(root, n))
        elif os.path.isfile(p):
            files.append(p)
        else:
            missing.append(p)
    if missing:
        raise ValueError("找不到 --massql-input: " + ", ".join(missing))
    seen, out = set(), []
    for f in files:
        a = os.path.abspath(f)
        if a not in seen:
            seen.add(a)
            out.append(f)
    return out


def _massql_write_mgf(path, matches):
    with open(path, "w", encoding="utf-8") as f:
        for m in matches:
            f.write("BEGIN IONS\n")
            f.write(f"TITLE={m.get('title') or 'massql'}\n")
            if m.get("precursor_mz") is not None:
                f.write(f"PEPMASS={m['precursor_mz']}\n")
            if m.get("charge") is not None:
                c = int(m["charge"])
                f.write(f"CHARGE={abs(c)}{'-' if c < 0 else '+'}\n")
            if m.get("scan") is not None:
                f.write(f"SCANS={m['scan']}\n")
            if m.get("rt_min") is not None:
                f.write(f"RTINSECONDS={m['rt_min'] * 60}\n")
            if m.get("polarity"):
                f.write(f"IONMODE={m['polarity']}\n")
            for mz, inten in m.get("peaks") or []:
                f.write(f"{mz} {inten}\n")
            f.write("END IONS\n\n")


def run_massql(args, ap):
    try:
        q = parse_massql(args.massql)
        files = _massql_collect_files(args.massql_input)
    except ValueError as e:
        ap.error(f"--massql: {e}")
    if not files:
        ap.error("--massql-input 没有找到 mzML/mzXML/MGF 文件")

    print("\n" + "=" * 60)
    print("  MassQL 本地检索")
    print("=" * 60)
    print(f"查询: {q['query']}")
    print(f"文件: {len(files)} 个")
    matches, errors = [], []
    n_spectra = 0
    for n, path in enumerate(files, 1):
        try:
            if path.lower().endswith((".mzml", ".mzxml")):
                specs = _massql_load_mzml(path)
            elif path.lower().endswith(".mgf"):
                specs = _massql_load_mgf(path)
            else:
                raise ValueError(f"不支持的文件格式: {path}")
        except Exception as e:
            msg = f"{type(e).__name__}: {e}"
            errors.append({"file": path, "error": msg})
            print(f"[MassQL] {n}/{len(files)} 跳过 {os.path.basename(path)}: {msg}")
            continue
        n_hit = 0
        for spec in specs:
            n_spectra += 1
            ok, matched = _massql_eval_query(q, spec)
            if ok:
                n_hit += 1
                matches.append({**spec, "n_peaks": len(spec["peaks"]), "matched": "; ".join(matched)})
        print(f"[MassQL] {n}/{len(files)} {os.path.basename(path)}: {len(specs)} 张 MS2, 命中 {n_hit}")

    if errors and n_spectra == 0:
        ap.error("所有输入文件都读取失败: " + errors[0]["error"])

    out_json = args.out or "massql_results.json"
    os.makedirs(os.path.dirname(os.path.abspath(out_json)), exist_ok=True)
    stem = out_json[:-5] if out_json.lower().endswith(".json") else out_json
    rows = []
    for m in matches:
        rows.append({"file": m["file"], "scan": m["scan"],
                     "precursor_mz": round(m["precursor_mz"], 5) if m.get("precursor_mz") is not None else None,
                     "rt_min": round(m["rt_min"], 3) if m.get("rt_min") is not None else None,
                     "polarity": m.get("polarity"), "charge": m.get("charge"), "n_peaks": m["n_peaks"],
                     "title": m.get("title"), "matched": m.get("matched")})
    tsv_path = stem + ".massql.tsv"
    write_tsv(tsv_path, rows, ["file", "scan", "precursor_mz", "rt_min", "polarity", "charge",
                               "n_peaks", "title", "matched"])
    outputs = {"massql_tsv": tsv_path}
    if args.massql_mgf:
        os.makedirs(os.path.dirname(os.path.abspath(args.massql_mgf)), exist_ok=True)
        _massql_write_mgf(args.massql_mgf, matches)
        outputs["massql_mgf"] = args.massql_mgf

    res = {"tool": "compound_db_check.py", "tool_version": TOOL_VERSION, "mode": "massql",
           "query": q["query"], "conditions": [c["raw"] for c in q["conditions"]],
           "inputs": files, "n_files": len(files), "n_spectra": n_spectra,
           "n_matches": len(matches), "errors": errors, "outputs": outputs, "matches": rows}
    with open(out_json, "w") as f:
        json.dump(res, f, indent=2, ensure_ascii=False)

    print("=" * 60)
    print(f"MS2 谱图: {n_spectra}; 命中: {len(matches)}; 文件错误: {len(errors)}")
    print(f"TSV: {tsv_path}")
    if args.massql_mgf:
        print(f"MGF: {args.massql_mgf}")
    print(f"JSON: {out_json}")
    return res


# ---------- 批量模式 ----------
BATCH_COLS = ["name", "search", "library", "status", "error", "task_id", "elapsed_s", "precursor", "n_query_peaks",
              "cosine", "min_peaks", "pm_tol", "frag_tol", "n_before_filter", "n_matches", "n_datasets",
              "max_cosine", "top_library_hit", "top_library_cosine", "n_verify_candidates", "n_verified_ok",
              "n_fetch_failed", "n_not_selected", "n_diag_pass", "n_diag_pass_datasets", "diag_pass_rate",
              "redu_unmatched_frac_hits", "redu_unmatched_frac_files", "redu_hit_files_in_bg",
              "redu_hit_datasets_in_bg", "top_SampleType", "top_NCBITaxonomy", "top_UBERONBodyPartName",
              "diag_pass_top_SampleType", "out_json"]


def batch_rows(res, out_json):
    rows = []
    qp = res.get("query_params", {})
    for key in [k for k in res if k.startswith("masst_")]:
        m, ri = res[key], res.get("run_info", {}).get(key, {})
        v = res.get("verification", {}).get(key, {})
        rd = res.get("redu", {}).get(key, {})
        r = {"name": res.get("name"), "search": key, "library": ri.get("library"),
             "status": "error" if "error" in m else "ok", "error": json.dumps(m["error"])[:300] if "error" in m else None,
             "task_id": m.get("task_id") or ri.get("task_id"), "elapsed_s": ri.get("elapsed_s"),
             "precursor": qp.get("precursor"), "n_query_peaks": len(qp.get("peaks") or []),
             **{k: qp.get(k) for k in ("cosine", "min_peaks", "pm_tol", "frag_tol")},
             "n_before_filter": m.get("n_before_filter"), "n_matches": m.get("n_matches"),
             "n_datasets": ri.get("n_datasets"), "max_cosine": m.get("max_cosine"),
             "n_verify_candidates": v.get("n_candidates"), "n_verified_ok": v.get("n_verified_ok"),
             "n_fetch_failed": v.get("n_fetch_failed"), "n_not_selected": v.get("n_not_selected"),
             "n_diag_pass": v.get("n_diag_pass"), "n_diag_pass_datasets": v.get("n_diag_pass_datasets"),
             "diag_pass_rate": v.get("diag_pass_rate"), "out_json": out_json}
        if m.get("library_hits"):
            r["top_library_hit"], r["top_library_cosine"] = m["library_hits"][0]["name"], m["library_hits"][0]["cosine"]
        a = rd.get("all")
        if a:
            r.update(redu_unmatched_frac_hits=a["unmatched_frac_hits"], redu_unmatched_frac_files=a["unmatched_frac_files"],
                     redu_hit_files_in_bg=a["n_hit_files_in_background"],
                     redu_hit_datasets_in_bg=a["n_hit_datasets_in_background"])
            stem = out_json[:-5]
            try:
                with open(stem + ".redu_summary.tsv") as f:
                    t = [x for x in csv.DictReader(f, delimiter="\t") if x["search"] == key]
                for col in ("SampleType", "NCBITaxonomy", "UBERONBodyPartName"):
                    top = [x for x in t if x["hit_set"] == "all" and x["column"] == col][:3]
                    r["top_" + col] = "; ".join(f"{x['value']} ({x['n_hit_datasets']}ds)" for x in top)
                top = [x for x in t if x["hit_set"] == "diag_pass" and x["column"] == "SampleType"][:3]
                r["diag_pass_top_SampleType"] = "; ".join(f"{x['value']} ({x['n_hit_datasets']}ds)" for x in top)
            except OSError:
                pass
        rows.append(r)
    return rows


def safe_name(s):
    return re.sub(r"[^A-Za-z0-9.+-]+", "_", s).strip("_") or "query"


def run_batch(args, ap, cli_opts):
    with open(args.queries) as f:
        qj = json.load(f)
    defaults, queries = ({}, qj) if isinstance(qj, list) else (qj.get("defaults", {}), qj.get("queries"))
    if not isinstance(queries, list) or not queries:
        ap.error("--queries: 需要 JSON 列表, 或 {\"defaults\": {...}, \"queries\": [...]}")
    bad = [k for k in defaults if k not in QUERY_KEYS]
    specs, names = [], set()
    for i, q in enumerate(queries):
        bad += [f"queries[{i}].{k}" for k in q if k not in QUERY_KEYS and not k.startswith("_")]
        if not q.get("name"):
            ap.error(f"--queries: 第 {i} 个查询缺少 name")
        sn = safe_name(q["name"])
        if sn in names:
            ap.error(f"--queries: 查询名重复 (按文件名规范化后): {q['name']}")
        names.add(sn)
        s = {**OPT_DEFAULTS, **cli_opts, **defaults, **{k: v for k, v in q.items() if not k.startswith("_")}}
        if sum(bool(s.get(k)) for k in ("peaks", "mgf", "usi", "from_hits")) != 1:
            ap.error(f"--queries: 查询 {q['name']} 需要 peaks / mgf / usi / from_hits 四者之一")
        diag = parse_mz_list(s.get("diagnostic")) + [m for m in parse_mz_list(s.get("diag_required"))
                                                     if not any(abs(m - d) < 1e-6 for d in parse_mz_list(s.get("diagnostic")))]
        try:
            if diag:
                parse_diag_rule(s.get("diag_rule", "all"), len(diag))
            resolve_libraries(s.get("library"))
            if s.get("redu_bg") not in REDU_BG_MODES:
                raise ValueError(f"redu_bg 只能是 {REDU_BG_MODES}")
        except ValueError as e:
            ap.error(f"--queries: 查询 {q['name']}: {e}")
        specs.append((sn, s))
    if bad:
        ap.error(f"--queries: 未知键 {bad}; 可用键: {sorted(QUERY_KEYS)}")
    out_dir = args.out or "masst_batch"
    print(f"[批量] {len(specs)} 个查询 -> {out_dir}/ (依次执行, 检索间隔 {args.pause} s)")
    for sn, s in specs:
        if s.get("usi") or s.get("mgf") or s.get("from_hits"):
            src = " ".join(f"{k} {s[k]}" for k in ("usi", "mgf", "from_hits") if s.get(k))
        else:
            src = f"{len(json.loads(s['peaks']) if isinstance(s['peaks'], str) else s['peaks'])} peaks"
        print(f"  {sn}: precursor {s.get('precursor')}, {src}, library {s['library']}, cos {s['cosine']}, "
              f"min_peaks {s['min_peaks']}, pm {s['pm_tol']}, frag {s['frag_tol']}, diag {s.get('diagnostic')}, "
              f"rule {s.get('diag_rule')}, verify_max {s.get('verify_max')}, redu {s.get('redu_join')}")
    if args.dry_run:
        print("[批量] --dry-run: 只检查, 未检索")
        return
    os.makedirs(out_dir, exist_ok=True)
    ctx = {"usi_cache": args.usi_cache or os.path.join(out_dir, "usi_cache"), "pause": args.pause}
    allrows = []
    summary_path = os.path.join(out_dir, "batch_summary.tsv")
    for n, (sn, s) in enumerate(specs, 1):
        oj = os.path.join(out_dir, sn, sn + ".json")
        print(f"\n######## [{n}/{len(specs)}] {s['name']} ########")
        if args.resume and os.path.exists(oj):
            try:
                with open(oj) as f:
                    prev = json.load(f)
                if all("error" not in prev[k] for k in prev if k.startswith("masst_")):
                    print(f"[批量] --resume: 已完成, 跳过 ({oj})")
                    allrows += batch_rows(prev, oj)
                    write_tsv(summary_path, allrows, BATCH_COLS)
                    continue
            except Exception:
                pass
        if n > 1 and args.pause:
            time.sleep(args.pause)
        try:
            res = run_query(s, oj, ctx, new_mode=True)
            allrows += batch_rows(res, oj)
        except Exception as e:
            print(f"[批量] 查询失败: {type(e).__name__}: {e}")
            allrows.append({"name": s["name"], "status": "error", "error": f"{type(e).__name__}: {e}"[:300],
                            "out_json": oj})
        write_tsv(summary_path, allrows, BATCH_COLS)  # 每个查询后都更新, 中断也不丢
    print(f"\n[批量] 汇总表: {summary_path}")


def main():
    ap = argparse.ArgumentParser(description="用 fastMASST 检索谱图在 GNPS 公共数据与参考谱库中的分布")
    ap.add_argument("--name", default=None, help="化合物名 (用于展示)")
    ap.add_argument("--inchikey", default=None)
    ap.add_argument("--precursor", type=float, default=None, help="前体 m/z")
    ap.add_argument("--peaks", default=None, help='MS2 峰列表 JSON, 如 "[[60.08,98618],[85.03,860576]]"')
    ap.add_argument("--mgf", default=None, help="或提供 MGF 文件路径 (单谱图)")
    ap.add_argument("--cosine", type=float, default=0.7)
    ap.add_argument("--min-peaks", type=int, default=3, help="最小匹配峰数 (低于此值的命中丢弃)")
    ap.add_argument("--from-parquet", default=None, help="in_house.parquet 路径 (自动提取某行)")
    ap.add_argument("--row", type=int, default=None, help="parquet 行号 (0 起)")
    ap.add_argument("--out", default=None, help="JSON 结果输出路径 (批量模式为输出目录)")
    g = ap.add_argument_group("Pan-ReDU 数据集检索 (不需要谱图)")
    g.add_argument("--find-datasets", action="store_true", help="按样本元数据找公开数据集")
    g.add_argument("--organism", default=None, help='物种名或 NCBI taxid, 如 "Mus musculus" 或 10090')
    g.add_argument("--body-part", default=None, help='组织 (UBERON 名), 逗号分隔表示任一, 如 "brain,cerebellum"')
    g.add_argument("--where", action="append", default=[], metavar="列=值",
                   help="其他元数据条件, 可重复, 如 SampleType=animal")
    g.add_argument("--min-files", type=int, default=1, help="数据集内符合条件的最少文件数 (大于 12 即填 13)")
    g.add_argument("--polarity", choices=["any", "positive", "negative", "both"], default="any",
                   help="符合条件的文件须包含的离子模式 (both = 正负都有; alternating 算两者都有)")
    g.add_argument("--redu-refresh", action="store_true", help="重新下载 Pan-ReDU 全表")
    n = ap.add_argument_group("检索参数与新功能 (2026-10-15 新增; 不给则与旧版行为相同)")
    n.add_argument("--library", default=None,
                   help="both (默认: 全库 + 参考谱库) | repository | library | 其他 fastMASST 索引名 (逗号分隔)")
    n.add_argument("--pm-tol", type=float, default=None, help="前体 m/z 容差 Da (默认 0.05)")
    n.add_argument("--frag-tol", type=float, default=None, help="碎片 m/z 容差 Da (默认 0.05)")
    n.add_argument("--max-wait", type=int, default=None, help="单次检索最长等待秒数 (默认 600)")
    n.add_argument("--usi", default=None, help="查询谱图来自 USI (经 metabolomics-usi.gnps2.org 取峰表)")
    n.add_argument("--from-hits", default=None,
                   help="重放: 读上次 --keep-all 的 <out>.hits.tsv (与同名 JSON), 不调用 fastMASST, "
                        "按新的 --min-peaks / 更高的 --cosine / 诊断离子 / ReDU 参数重新分析")
    n.add_argument("--query-min-rel", type=float, default=None, help="去掉查询谱图中低于基峰 X%% 的峰")
    n.add_argument("--keep-all", action="store_true", default=None,
                   help="保存 API 返回的全部命中 (全部字段) 到 <out>.hits.tsv")
    n.add_argument("--diagnostic", default=None, help='诊断离子 m/z, 逗号分隔, 如 "265.2526,247.2420"')
    n.add_argument("--diag-required", default=None, help="必需诊断离子 m/z (逗号分隔, 自动并入 --diagnostic)")
    n.add_argument("--diag-tol", type=float, default=None, help="诊断离子容差 Da (默认 0.01)")
    n.add_argument("--diag-min-rel", type=float, default=None, help="诊断离子最低相对强度 %% 基峰 (默认 0)")
    n.add_argument("--diag-rule", default=None, help="all (默认) | any | k | k-of-n")
    n.add_argument("--verify-max", type=int, default=None, help="每次检索最多核对的命中谱图数 (默认 200)")
    n.add_argument("--verify-seed", type=int, default=None, help="分层抽样的种子 (默认 0)")
    n.add_argument("--usi-cache", default=None, help="谱图缓存目录 (默认: 输出 JSON 同目录 usi_cache/)")
    n.add_argument("--usi-rate", type=float, default=None, help="USI 取谱限速 次/秒 (默认 2.5, 上限 3)")
    n.add_argument("--redu-join", action="store_true", default=None, help="命中关联 Pan-ReDU 元数据并做富集")
    n.add_argument("--redu-bg", choices=REDU_BG_MODES, default=None,
                   help="富集背景 (默认 pos_lc_or_unknown: 排除已标注负离子与 GC/EI; pos_lc: 只要已标注正离子)")
    mq = ap.add_argument_group("MassQL 本地检索 (常用 MS2 子集; 不需要 fastMASST)")
    mq.add_argument("--massql", default=None,
                    help='MassQL 查询, 如 "QUERY scaninfo(MS2DATA) WHERE MS2PROD=265.2526:TOLERANCEMZ=0.01 '
                         'AND MS2PROD=247.2420:TOLERANCEMZ=0.01" (支持 MS2PROD/MS2NL/MS2PREC/RT/POLARITY 等)')
    mq.add_argument("--massql-input", action="append", default=[], metavar="PATH",
                    help="mzML/mzXML/MGF 文件或目录, 可重复")
    mq.add_argument("--massql-mgf", default=None, help="把命中的 MS2 谱写成 MGF (可选)")
    b = ap.add_argument_group("批量模式")
    b.add_argument("--queries", default=None, help="查询 JSON 文件 (列表, 或 {defaults, queries})")
    b.add_argument("--pause", type=float, default=5.0, help="两次 fastMASST 检索之间的间隔秒数 (批量模式, 默认 5)")
    b.add_argument("--resume", action="store_true", help="跳过已成功完成的查询")
    b.add_argument("--dry-run", action="store_true", help="只解析查询文件并打印计划, 不检索")
    args = ap.parse_args()

    if args.massql:
        return run_massql(args, ap)
    if args.massql_input or args.massql_mgf:
        ap.error("--massql-input / --massql-mgf 需要配合 --massql")

    if args.find_datasets:
        return run_find_datasets(args, ap)
    if args.organism or args.body_part or args.where or args.redu_refresh:
        ap.error("--organism / --body-part / --where / --redu-refresh 需要配合 --find-datasets")

    # 命令行显式给出的选项 (None = 未给)
    cli_opts = {k: getattr(args, k) for k in NEW_OPTS
                if k not in ("usi_cache", "usi", "from_hits") and getattr(args, k) is not None}
    if args.cosine != 0.7:
        cli_opts["cosine"] = args.cosine
    if args.min_peaks != 3:
        cli_opts["min_peaks"] = args.min_peaks
    if cli_opts.get("usi_rate", 2.5) > 3:
        ap.error("--usi-rate 不能超过 3 次/秒")
    if args.queries:
        return run_batch(args, ap, cli_opts)
    if args.resume or args.dry_run:
        ap.error("--resume / --dry-run 需要配合 --queries")
    new_mode = any(getattr(args, k) is not None for k in NEW_OPTS)

    # 从 parquet 提取
    if args.from_parquet is not None and args.row is not None:
        import duckdb
        df = duckdb.sql(f"SELECT * FROM '{args.from_parquet}'").df()
        r = df.iloc[args.row]
        args.name = args.name or r.get("name_hp")
        args.inchikey = args.inchikey or r.get("inchikey")
        args.precursor = args.precursor or float(r["precursor_mz"])
        if args.peaks is None:
            args.peaks = r["ms2"] if isinstance(r["ms2"], str) else json.dumps(r["ms2"])
        print(f"[parquet] 行 {args.row}: {args.name}, precursor={args.precursor}")

    if not ((args.peaks or args.mgf) and (args.precursor or args.mgf) or args.usi or args.from_hits):
        ap.error("需要提供谱图: --peaks 与 --precursor, 或 --mgf, 或 --from-parquet 与 --row, 或 --usi, "
                 "或 --from-hits")

    name = args.name or args.inchikey or "unknown"
    spec = {**OPT_DEFAULTS, **cli_opts, "cosine": args.cosine, "min_peaks": args.min_peaks,
            "name": name, "inchikey": args.inchikey, "precursor": args.precursor,
            "peaks": args.peaks if not args.mgf else None, "mgf": args.mgf, "usi": args.usi,
            "from_hits": args.from_hits}
    out = args.out or f"dbcheck_{name.replace(' ', '_').replace(':', '')}.json"
    try:
        run_query(spec, out, {"usi_cache": args.usi_cache}, new_mode=new_mode)
    except (ValueError, RuntimeError) as e:
        ap.error(str(e))


if __name__ == "__main__":
    main()
