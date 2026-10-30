"""MS1 精确质量检索（页面 ms1_search）：观测 m/z + 加合物 -> 中性质量窗口 -> 在 PubChem 查单同位素质量落在窗口内的化合物。

全部在线、不下载。所有库都走 PubChem：
  PubChem                                 不加来源过滤；
  HMDB / DrugBank / ChEBI / KEGG / LIPID MAPS / NORMAN
                                          同一个 Entrez 检索加 [SourceName] 过滤，即 PubChem 里这些库的副本
                                          （HMDB、DrugBank 拒绝脚本下载，这是唯一可脚本化的路）。
流程：每个库一次 esearch（各加合物窗口 OR 在一个 term 里）-> PUG REST 批量取属性 -> 用与 Theoretical m/z 页相同的加合物定义
(msbuddy_engine.ion_mz_from_mass) 重算 m/z 和 Δ ppm，窗口外的丢掉 -> PUG REST xrefs 取各库自己的 ID。

实测过的 PubChem 行为：
  - MonoisotopicMass 对永久带电的离子（胆碱阳离子 C5H14NO+ 等）是原子质量之和，不扣电子；所以 [M]+ 的 m/z = 质量 - 电子，
    反推质量 = m/z + 电子，与 Theoretical m/z 页的 [M]+ 一致；带电化合物只匹配 [M]+ / [M]-。
  - esearch 默认新 CID 在前，没有可用的 sort；DCNT（DepositorCount）可做范围过滤，所以超过上限时按 DCNT 分档取，
    每档内取最小的 CID（老 CID 多半是知名化合物）。
  - xrefs RegistryID 不带来源（各家编号混在一起），而且冷缓存时很慢（20 个 CID 可达 5 秒）：只给确认属于 HMDB / DrugBank / ChEBI /
    KEGG / LIPID MAPS 的化合物取，分批并发，超过 ID_BUDGET 秒就不等（没取到的显示库名，请求继续跑进缓存，再查一次就有）。
    KEGG 的 C 号（化合物）/ D 号（药物）会和别家编号撞，所以必须先确认属于 KEGG；PubChem 里 NORMAN 没有自己的编号，只记成员关系。
  - NCBI 限速：不带 API key 最多 3 次/秒（带 tool=mzdesk，不发邮箱）；PUG REST 最多 5 次/秒。
"""
import concurrent.futures as cf
import json
import re
import socket
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from functools import lru_cache

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid"
FASTFORMULA = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula"
TIMEOUT = 20.0
MAX_PER_DB = 100        # 每个库最多列出多少个化合物
FORMULA_CAP = 100       # 同分子式表最多考察多少个 PubChem 结构（再多 xrefs 太慢）
ISOBARIC_CAP = 50       # 同量异位表最多列多少个分子式
SOURCE_BUDGET = 8.0     # 同分子式表最多等多少秒拿库归属
BATCH = 200             # PUG REST 属性每次 POST 的 CID 数
XREF_BATCH = 15         # xrefs 每次 POST 的 CID 数（小批并发，冷缓存时每个 CID 都要服务器现算）
XREF_WORKERS = 4
ID_BUDGET = 10.0        # 最多等多少秒拿原生 ID
MIN_GAP = {"eutils": 0.35, "pug": 0.22}      # 两次请求的最小间隔（秒）
MIN_NEUTRAL_MASS = 1.0  # 小于它的窗口没有意义（也会出现负数质量）

# 页面上的库名 -> PubChem 的来源名（已逐个核对：每个都能查到咖啡因 / PC 34:1 等已知化合物）；PubChem 本身不加过滤
SOURCES = {         # 顺序与 spec.DATABASES 一致
    "HMDB": "Human Metabolome Database (HMDB)",
    "PubChem": None,
    "LIPID MAPS": "LIPID MAPS",
    "KEGG": "KEGG",
    "ChEBI": "ChEBI",
    "DrugBank": "DrugBank",
    "NORMAN SusDat": "NORMAN Suspect List Exchange",   # PubChem 里的 NORMAN 来源是整个 SLE，SusDat 是其中的核心清单
}
PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/compound/{}"
DB_URL = {      # 各库自己的页面；NORMAN 没有，用 PubChem 页面
    "HMDB": "https://hmdb.ca/metabolites/{}",
    "DrugBank": "https://go.drugbank.com/drugs/{}",
    "ChEBI": "https://www.ebi.ac.uk/chebi/searchId.do?chebiId={}",
    "KEGG": "https://www.kegg.jp/entry/{}",
    "LIPID MAPS": "https://www.lipidmaps.org/databases/lmsd/{}",
}
# xrefs RegistryID 里认得出的原生 ID（形状已在真实数据上核对）
_ID_PATTERNS = {
    "HMDB": re.compile(r"HMDB\d{7}", re.I),
    "DrugBank": re.compile(r"DB\d{5}"),
    "ChEBI": re.compile(r"CHEBI:\d+", re.I),
    "LIPID MAPS": re.compile(r"LM[A-Z]{2}\d{8}"),
}
_KEGG_PATTERNS = (re.compile(r"C\d{5}"), re.compile(r"D\d{5}"))     # 化合物优先，其次药物；只在确认属于 KEGG 时用
_ID_DBS = set(_ID_PATTERNS) | {"KEGG"}      # 能取到原生 ID 的库（NORMAN 没有）
# 超过上限时的 DCNT 分档（含两端）：先取多家机构收录的，档内取最小的 CID
POPULARITY_TIERS = ((100, 1_000_000), (30, 99), (10, 29), (3, 9), (1, 2))


def _bd():
    """加合物定义所在的引擎。延后导入：msbuddy 要几秒才载入，而联网部分（和测试里替换它）用不到。"""
    from mzdesk.ms1 import msbuddy_engine
    return msbuddy_engine


class ExactMassError(RuntimeError):
    """PubChem 联网失败或回复异常。消息直接给用户看。"""


@dataclass
class Hit:
    cid: int
    name: str
    formula: str
    mass: float            # PubChem 单同位素质量
    adduct: str
    mz: float              # 该化合物作为这个加合物的理论 m/z
    ppm: float             # (观测 - 理论) / 理论 * 1e6
    dbs: list = field(default_factory=list)       # 确认收录它的库（不含 PubChem），按 SOURCES 的顺序
    ids: dict = field(default_factory=dict)       # 库名 -> 原生 ID
    tier: int = 0                                 # 超过上限时它所在的 DCNT 档（0 最好）
    refs: int = 0                                 # PubChem 里这个化合物的注册号个数（知名度的粗略指标；只对取过 ID 的化合物有）

    @property
    def order(self):
        """同一组里的排序：先专库里有的，再 DCNT 档靠前的，再被引用多的，最后 CID 小的"""
        return (0 if self.dbs else 1, self.tier, -self.refs, self.cid)


@dataclass
class Result:
    hits: list             # 去重（按 CID）后的命中，已按 |Δ ppm| 排序
    totals: dict           # 库名 -> PubChem 报告的窗口内化合物总数（未封顶）
    shown: dict            # 库名 -> 实际取回的数量
    windows: list          # [(加合物, 质量下限, 质量上限)]
    warnings: list = field(default_factory=list)


# ---------------------------------------------------------------- network
_slot_lock = threading.Lock()
_next_slot = {"eutils": 0.0, "pug": 0.0}


def _wait_turn(kind):
    """占下一个请求时间槽（距上一个至少 MIN_GAP[kind]）并等到那一刻。多线程共用，所以不会超速。"""
    with _slot_lock:
        now = time.monotonic()
        start = max(now, _next_slot[kind])
        _next_slot[kind] = start + MIN_GAP[kind]
    if start > now:
        time.sleep(start - now)


def _http(url, data=None):
    """唯一的联网入口（测试时 monkeypatch 它）：GET，或带 data（表单编码的 bytes）时 POST。PUG 的 404（没有这些 CID）按空结果处理。"""
    kind = "eutils" if url.startswith(ESEARCH) else "pug"
    headers = {"User-Agent": "mzDesk", "Accept": "application/json"}
    if data:
        headers["Content-Type"] = "application/x-www-form-urlencoded"
    for attempt in (1, 2):
        _wait_turn(kind)
        try:
            with urllib.request.urlopen(urllib.request.Request(url, data=data, headers=headers), timeout=TIMEOUT) as r:
                return r.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return "{}"
            if exc.code in (429, 503):
                if attempt == 1:
                    time.sleep(1.5)
                    continue
                raise ExactMassError(f"PubChem is busy (HTTP {exc.code}). Try again in a minute.") from None
            raise ExactMassError(f"PubChem answered with HTTP {exc.code}. Try again later.") from None
        except (urllib.error.URLError, socket.timeout, TimeoutError, OSError):
            raise ExactMassError("PubChem could not be reached (timeout or no network). Try again later.") from None


@lru_cache(maxsize=256)
def _get_json(url, data=None):
    """解析后的 JSON，按请求缓存（失败不缓存）。返回值只读。"""
    try:
        return json.loads(_http(url, data))
    except ValueError:
        raise ExactMassError("PubChem sent an unexpected reply. Try again later.") from None


def _esearch(term, retmax, retstart=0):
    """(命中总数, 本页 CID 列表)。默认顺序是新 CID 在前。"""
    query = urllib.parse.urlencode({"db": "pccompound", "term": term, "retmode": "json", "retmax": retmax,
                                    "retstart": retstart, "tool": "mzdesk"})
    res = _get_json(f"{ESEARCH}?{query}").get("esearchresult") or {}
    if res.get("ERROR") or res.get("errorlist"):
        raise ExactMassError("PubChem rejected the search. Try again later.")
    for phrase in (res.get("warninglist") or {}).get("quotedphrasesnotfound") or []:
        if "[SourceName]" in phrase:      # 来源名写错会静默返回 0 条，所以这里报错
            raise ExactMassError(f"PubChem does not know the source {phrase.split('[')[0]}.")
    try:
        return int(res["count"]), [int(c) for c in res.get("idlist") or []]
    except (KeyError, TypeError, ValueError):
        raise ExactMassError("PubChem sent an unexpected reply. Try again later.") from None


def _pug_post(path, cids):
    """PUG REST 的批量 POST（每批 BATCH 个 CID），逐批给出解析后的回复。"""
    for i in range(0, len(cids), BATCH):
        body = urllib.parse.urlencode({"cid": ",".join(str(c) for c in cids[i:i + BATCH])}).encode()
        yield _get_json(f"{PUG}/{path}", body)


def _properties(cids):
    """{cid: {"formula", "mass", "name", "charge", "inchikey"}}"""
    out = {}
    for reply in _pug_post("property/MolecularFormula,MonoisotopicMass,Title,Charge,InChIKey/JSON", cids):
        for p in (reply.get("PropertyTable") or {}).get("Properties") or []:
            try:
                out[int(p["CID"])] = {"formula": p.get("MolecularFormula") or "", "mass": float(p["MonoisotopicMass"]),
                                      "name": p.get("Title") or "", "charge": int(p.get("Charge") or 0),
                                      "inchikey": p.get("InChIKey") or ""}
            except (KeyError, TypeError, ValueError):
                continue
    return out


def _registry_batch(cids):
    out = {}
    body = urllib.parse.urlencode({"cid": ",".join(str(c) for c in cids)}).encode()
    for info in (_get_json(f"{PUG}/xrefs/RegistryID/JSON", body).get("InformationList") or {}).get("Information") or []:
        try:
            out[int(info["CID"])] = list(info.get("RegistryID") or [])
        except (KeyError, TypeError, ValueError):
            continue
    return out


def _registry_ids(cids, budget=None):
    """({cid: [注册号, ...]}, 没取到的 CID 数)。分批并发（请求间隔由 _wait_turn 保证），cids 靠前的先取；
    超过 budget 秒就不等了：没跑完的请求继续跑、结果进缓存，所以再查一次通常就齐了。"""
    chunks = [cids[i:i + XREF_BATCH] for i in range(0, len(cids), XREF_BATCH)]
    pool = cf.ThreadPoolExecutor(max_workers=XREF_WORKERS)
    futures = [pool.submit(_registry_batch, c) for c in chunks]
    done, _ = cf.wait(futures, timeout=ID_BUDGET if budget is None else budget)
    pool.shutdown(wait=False, cancel_futures=True)
    out, missing = {}, 0
    for future, chunk in zip(futures, chunks):
        if future in done and future.exception() is None:
            out.update(future.result())
        else:
            missing += len(chunk)
    return out, missing


def native_ids(registry_ids, member_dbs):
    """{库名: 原生 ID}：只看化合物确认属于的库。注册号里各家编号混在一起，所以按 ID 的形状挑。"""
    ids = {}
    for db in member_dbs:
        pats = _KEGG_PATTERNS if db == "KEGG" else (_ID_PATTERNS[db],) if db in _ID_PATTERNS else ()
        for pat in pats:       # KEGG：化合物 (C) 优先，其次药物 (D)
            found = next((r for r in registry_ids if pat.fullmatch(r)), None)
            if found:
                ids[db] = found.upper() if db in ("HMDB", "ChEBI") else found
                break
    return ids


def _source_batch(cids):
    out = {}
    body = urllib.parse.urlencode({"cid": ",".join(str(c) for c in cids)}).encode()
    reply = _get_json(f"{PUG}/xrefs/SourceName/JSON", body)
    for info in (reply.get("InformationList") or {}).get("Information") or []:
        try:
            out[int(info["CID"])] = set(info.get("SourceName") or [])
        except (KeyError, TypeError, ValueError):
            continue
    return out


def _source_names(cids, budget=SOURCE_BUDGET):
    """({cid: {PubChem 来源名, ...}}, 没取到的 CID 数)；与 _registry_ids 相同的分批并发与时间预算。"""
    chunks = [cids[i:i + XREF_BATCH] for i in range(0, len(cids), XREF_BATCH)]
    pool = cf.ThreadPoolExecutor(max_workers=XREF_WORKERS)
    futures = [pool.submit(_source_batch, c) for c in chunks]
    done, _ = cf.wait(futures, timeout=budget)
    pool.shutdown(wait=False, cancel_futures=True)
    out, missing = {}, 0
    for future, chunk in zip(futures, chunks):
        if future in done and future.exception() is None:
            out.update(future.result())
        else:
            missing += len(chunk)
    return out, missing


# ---------------------------------------------------------------- masses
def mass_windows(mz, adducts, tol_ppm):
    """[(加合物, 质量下限, 质量上限)]：理论 m/z 离观测 m/z 不超过 tol_ppm 时，中性质量的取值范围。

    |观测 - 理论| / 理论 <= tol  <=>  理论 m/z 在 [obs / (1 + tol), obs / (1 - tol)]；再按加合物的定义换成质量
    ([2M+H]+ 等 m = 2 的加合物，窗口自然只有一半宽）。质量小于 MIN_NEUTRAL_MASS 的窗口丢掉。
    """
    bd, tol, out = _bd(), tol_ppm * 1e-6, []
    for adduct in adducts:
        lo = bd.mass_from_ion_mz(mz / (1 + tol), adduct)
        hi = bd.mass_from_ion_mz(mz / (1 - tol), adduct)
        if lo is None:
            raise ValueError(f"Unknown adduct {adduct}.")
        if hi >= MIN_NEUTRAL_MASS:
            out.append((adduct, max(lo, 0.0), hi))
    return out


def _mass_term(windows):
    """Entrez 检索式：所有窗口 OR 在一起。窗口外扩 5e-7 Da（6 位小数的取整余量），精确的筛选在取回属性后做。"""
    return "(" + " OR ".join(f"{lo - 5e-7:.6f}:{hi + 5e-7:.6f}[MonoisotopicMass]" for _, lo, hi in windows) + ")"


def _best_adduct(prop, mz, adducts):
    """(加合物, 理论 m/z, Δ ppm)：Δ 最小的那个。永久带电的化合物只配 [M]+ / [M]-。"""
    bd, best = _bd(), None
    for adduct in adducts:
        if prop["charge"] and adduct != ("[M]+" if prop["charge"] > 0 else "[M]-"):
            continue
        theo = bd.ion_mz_from_mass(prop["mass"], adduct)
        ppm = bd.ppm_error(mz, theo)
        if best is None or abs(ppm) < abs(best[2]):
            best = (adduct, theo, ppm)
    return best


# ---------------------------------------------------------------- search
def _pick(term, cap):
    """(PubChem 报告的总数, [(cid, 档)] 按显示顺序)。

    不超过 cap：全部，按 CID 升序（档 0）。超过：按 DCNT 分档（POPULARITY_TIERS，多家机构收录的先取），档内取最小的 CID；
    DCNT 取不到东西时退回全体里最小的 CID。
    """
    total, cids = _esearch(term, cap)
    if total <= cap:
        return total, [(c, 0) for c in sorted(cids)]
    picked = []
    for tier, (lo, hi) in enumerate(POPULARITY_TIERS):
        need = cap - len(picked)
        sub = f"{term} AND {lo}:{hi}[DCNT]"
        n, ids = _esearch(sub, need)
        if n > need:       # 默认新 CID 在前，要最老的就从末尾取
            _, ids = _esearch(sub, need, retstart=n - need)
        picked += [(c, tier) for c in sorted(ids)]
        if len(picked) >= cap:
            break
    if not picked:
        _, ids = _esearch(term, cap, retstart=total - cap)
        picked = [(c, 0) for c in sorted(ids)]
    return total, picked


def search(mz, adducts, tol_ppm, databases, cap=MAX_PER_DB, with_ids=True):
    """观测 m/z 在给定加合物、容差下的数据库命中。

    adducts: 要试的加合物（页面上的 Auto 由调用方展开成该离子模式的全部加合物）。databases: SOURCES 里的库名。
    每个化合物只出现一次（按 CID），归到 |Δ ppm| 最小的那个加合物；每个库最多取 cap 个。
    with_ids=False 时不取各库的原生 ID（只需要库名时省掉最慢的一步）。
    """
    unknown = [d for d in databases if d not in SOURCES]
    if unknown:
        raise ValueError(f"Unknown database: {', '.join(unknown)}.")
    windows = mass_windows(mz, adducts, tol_ppm)
    if not windows:
        raise ValueError("The m/z is too small for the chosen adducts.")
    term = _mass_term(windows)
    chosen = [d for d in SOURCES if d in databases]      # 固定顺序，请求次序也就固定

    totals, picks = {}, {}
    for db in chosen:
        src = SOURCES[db]
        totals[db], picks[db] = _pick(term if src is None else f'{term} AND "{src}"[SourceName]', cap)
    shown = {db: len(p) for db, p in picks.items()}
    tier_of = {}
    for pairs in picks.values():
        for cid, tier in pairs:
            tier_of[cid] = min(tier, tier_of.get(cid, tier))
    members = {db: {cid for cid, _ in pairs} for db, pairs in picks.items()}

    props = _properties(sorted(tier_of)) if tier_of else {}
    hits = []
    for cid, prop in props.items():
        best = _best_adduct(prop, mz, adducts)
        if best is None or abs(best[2]) > tol_ppm + 1e-9:       # 窗口边缘（取整余量、电荷不符）的丢掉
            continue
        dbs = [db for db in chosen if db != "PubChem" and cid in members[db]]
        hits.append(Hit(cid, prop["name"] or prop["formula"], prop["formula"], prop["mass"], best[0], best[1], best[2], dbs,
                        tier=tier_of[cid]))

    rank = {a: i for i, a in enumerate(adducts)}

    def arrange():
        # 同一个离子式（caffeine [M+H]+ 与 C8H12N4O3 [M+H-H2O]+）的 |Δ| 只差 PubChem 质量取整的 1e-5 ppm：按 0.001 ppm 取整后并列，先排靠前的加合物
        hits.sort(key=lambda h: (round(abs(h.ppm), 3), rank[h.adduct], h.formula) + h.order)

    arrange()
    warnings = []
    want = [h for h in hits if any(d in _ID_DBS for d in h.dbs)] if with_ids else []
    if want:
        regs, missing = _registry_ids([h.cid for h in want])
        for h in want:
            h.ids = native_ids(regs.get(h.cid, []), h.dbs)
            h.refs = len(regs.get(h.cid, []))
        arrange()
        if missing:
            warnings.append(f"PubChem was slow: the database IDs of {missing} compounds are not shown yet (the rows give the database name). "
                            "Run the search again in a minute.")
    for db in chosen:
        if totals[db] > cap:
            warnings.append(f"{db}: {totals[db]:,} compounds in the mass window; the {cap} best known are listed. "
                            "Narrow the tolerance or pick one adduct to see fewer.")
    return Result(hits, totals, shown, windows, warnings)


# ---------------------------------------------------------------- alternative explanations
_FORMULA_RE = re.compile(r"^(?:[A-Z][a-z]?\d*)+$")
_SOURCE_DB = {name: db for db, name in SOURCES.items() if name}


def _neutral_formula(formula):
    """去掉末尾的电荷（C5H14NO+ -> C5H14NO），只用于按分子式比较。"""
    return re.sub(r"[+-]\d*$", "", (formula or "").strip())


def _fastformula(formula, cap):
    """(CID 列表, 是否超过 cap)。fastformula 是按分子式精确检索，不像 esearch 的文本检索那样把式量当字符串。"""
    if not _FORMULA_RE.fullmatch(formula or ""):
        raise ValueError(f"Invalid molecular formula {formula!r}.")
    url = f"{FASTFORMULA}/{urllib.parse.quote(formula, safe='')}/cids/JSON?MaxRecords={int(cap) + 1}"
    reply = _get_json(url)
    try:
        cids = [int(c) for c in (reply.get("IdentifierList") or {}).get("CID") or []]
    except (TypeError, ValueError):
        raise ExactMassError("PubChem sent an unexpected reply. Try again later.") from None
    return cids[:cap], len(cids) > cap


def formula_groups(formula, databases, cap=FORMULA_CAP):
    """同一分子式的候选结构，按 InChIKey 前 14 位（连接方式）合并。

    返回 {"formula", "rows": [{"key", "name", "names", "databases", "cids", "stereoisomers"}], "examined",
          "truncated", "warnings"}。databases 是所选库里收录它的库；只选了 PubChem 时写 PubChem。
    """
    unknown = [d for d in databases if d not in SOURCES]
    if unknown:
        raise ValueError(f"Unknown database: {', '.join(unknown)}.")
    formula = _neutral_formula(formula)
    if not formula:
        raise ValueError("Enter a molecular formula.")
    chosen = [d for d in SOURCES if d in databases]
    cids, truncated = _fastformula(formula, cap)
    props = {cid: p for cid, p in _properties(cids).items() if _neutral_formula(p["formula"]) == formula}
    sources, missing = _source_names(sorted(props)) if props else ({}, 0)

    groups = {}
    for cid in sorted(props):
        p = props[cid]
        key = (p["inchikey"] or f"cid:{cid}")[:14]
        g = groups.setdefault(key, {"key": key, "cids": [], "names": [], "dbs": set(), "stereo": set()})
        g["cids"].append(cid)
        if p["name"] and p["name"] not in g["names"]:
            g["names"].append(p["name"])
        if p["inchikey"]:
            g["stereo"].add(p["inchikey"])
        g["dbs"].update(_SOURCE_DB[s] for s in sources.get(cid, ()) if s in _SOURCE_DB)

    rows = []
    for g in groups.values():
        member = [db for db in chosen if db != "PubChem" and db in g["dbs"]]
        if "PubChem" not in chosen and not member:
            continue
        rows.append({"key": g["key"], "name": g["names"][0] if g["names"] else "", "names": g["names"][1:],
                     "databases": member or (["PubChem"] if "PubChem" in chosen else []),
                     "cids": g["cids"], "stereoisomers": max(0, len(g["stereo"]) - 1)})
    rows.sort(key=lambda g: (-len([d for d in g["databases"] if d != "PubChem"]), min(g["cids"])))

    warnings = []
    if truncated:
        warnings.append(f"PubChem lists more than {cap} structures with {formula}; only the first {cap} it returned were examined.")
    if missing:
        warnings.append(f"PubChem was slow: the database membership of {missing} structures is incomplete. Run again in a minute.")
    return {"formula": formula, "rows": rows, "examined": len(props), "truncated": truncated, "warnings": warnings}


def isobaric(mz, adducts, tol_ppm, databases, exclude_formula=None, cap=ISOBARIC_CAP):
    """同一 m/z 窗口里分子式不同的化合物，按分子式合并（每式一行，取 |Δ ppm| 最小的化合物）。

    exclude_formula: 查询物自己的（中性）分子式；同式的命中不列在这里。with_ids 永远是 False：这里只要库名。
    """
    res = search(mz, adducts, tol_ppm, databases, with_ids=False)
    exclude = _neutral_formula(exclude_formula) if exclude_formula else None
    groups, same = {}, 0
    for h in res.hits:
        f = _neutral_formula(h.formula)
        if exclude and f == exclude:
            same += 1
            continue
        g = groups.setdefault(f, {"formula": h.formula, "hits": [], "names": []})
        g["hits"].append(h)
        if h.name and h.name not in g["names"]:
            g["names"].append(h.name)

    rows = []
    for g in groups.values():
        best = min(g["hits"], key=lambda h: (abs(h.ppm), h.order))
        dbs = [db for db in SOURCES if db in {d for h in g["hits"] for d in h.dbs}]
        rows.append({"formula": g["formula"], "name": best.name, "names": [n for n in g["names"] if n != best.name],
                     "databases": dbs or (["PubChem"] if "PubChem" in databases else []),
                     "mz": best.mz, "ppm": best.ppm, "adduct": best.adduct, "compounds": len(g["hits"])})
    rows.sort(key=lambda r: (round(abs(r["ppm"]), 3), r["formula"]))
    warnings = list(res.warnings)
    if len(rows) > cap:
        warnings.append(f"{len(rows)} isobaric formulas are within {tol_ppm:g} ppm; the {cap} closest are listed.")
    return {"rows": rows[:cap], "warnings": warnings, "same_formula": same, "total": len(rows)}
