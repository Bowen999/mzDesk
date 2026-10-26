"""resolve.py

名称 / ID / 结构 的互转与跨库链接。输入一个字符串（名称、CAS、SMILES、InChI、
InChIKey、PubChem CID、HMDB、KEGG、ChEBI、LIPID MAPS、DrugBank），返回结构与其他标识符。

结构互转走 RDKit（本地）；名称与跨库 ID 走 PubChem PUG REST 和 UniChem。
所有网络访问经过唯一的 `_http`（带 sqlite 磁盘缓存），测试时只需 monkeypatch 它。
"""
from __future__ import annotations

import json
import os
import re
import socket
import sqlite3
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import inchi as rd_inchi

KINDS = ["Auto", "Name", "SMILES", "InChI", "InChIKey", "CAS", "PubChem CID",
         "HMDB", "KEGG", "ChEBI", "LIPID MAPS", "DrugBank"]

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
UNICHEM = "https://www.ebi.ac.uk/unichem/api/v1"
CACHE_TTL = 30 * 86400
MAX_CANDIDATES = 10
MAX_SYNONYMS = 50
_MIN_INTERVAL = 0.2  # PubChem: 5 req/s

_RE_KEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_RE_CAS = re.compile(r"^\d{2,7}-\d{2}-\d$")
_RE_HMDB = re.compile(r"^HMDB\d{5,7}$")
_RE_KEGG = re.compile(r"^C\d{5}$")
_RE_CHEBI = re.compile(r"^CHEBI:\d+$", re.I)
_RE_LM = re.compile(r"^LM[A-Z]{2}\d{8,12}$")
_RE_DB = re.compile(r"^DB\d{5}$")

# UniChem shortName -> 本模块的 ID 标签
_UNICHEM_MAP = {"hmdb": "HMDB", "chebi": "ChEBI", "drugbank": "DrugBank", "comptox": "CompTox",
                "chembl": "ChEMBL", "lipidmaps": "LIPID MAPS", "kegg_ligand": "KEGG"}
# 反向查询用的 UniChem 数据源名
_KIND_SOURCE = {"HMDB": "hmdb", "KEGG": "kegg_ligand", "ChEBI": "chebi",
                "LIPID MAPS": "lipidmaps", "DrugBank": "drugbank"}

_LINKS = {
    "PubChem": "https://pubchem.ncbi.nlm.nih.gov/compound/{}",
    "HMDB": "https://hmdb.ca/metabolites/{}",
    "KEGG": "https://www.kegg.jp/entry/{}",
    "ChEBI": "https://www.ebi.ac.uk/chebi/searchId.do?chebiId={}",
    "LIPID MAPS": "https://www.lipidmaps.org/databases/lmsd/{}",
    "DrugBank": "https://go.drugbank.com/drugs/{}",
    "CompTox": "https://comptox.epa.gov/dashboard/chemical/details/{}",
    "ChEMBL": "https://www.ebi.ac.uk/chembl/compound_report_card/{}/",
    "CAS": "https://commonchemistry.cas.org/detail?cas_rn={}",
}


class ResolveError(ValueError):
    """User-facing: not found / network unavailable / cannot parse."""


class _NotFound(Exception):
    """HTTP 404 from the remote service."""


@dataclass
class Resolved:
    query: str
    kind: str
    smiles: Optional[str] = None
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    cid: Optional[int] = None
    name: Optional[str] = None
    iupac: Optional[str] = None
    ids: dict = field(default_factory=dict)
    synonyms: list = field(default_factory=list)
    candidates: list = field(default_factory=list)
    stereo_fallback: bool = False
    links: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)


# ---------------------------------------------------------------- detect

def _cas_ok(s: str) -> bool:
    if not _RE_CAS.match(s):
        return False
    digits = s.replace("-", "")
    body, check = digits[:-1], int(digits[-1])
    total = sum((i + 1) * int(d) for i, d in enumerate(reversed(body)))
    return total % 10 == check


def detect(query: str) -> str:
    """Guess the input kind (never returns 'Auto')."""
    q = (query or "").strip()
    if q.startswith("InChI="):
        return "InChI"
    if _RE_KEY.match(q):
        return "InChIKey"
    if _RE_HMDB.match(q):
        return "HMDB"
    if _RE_KEGG.match(q):
        return "KEGG"
    if _RE_CHEBI.match(q):
        return "ChEBI"
    if _RE_LM.match(q):
        return "LIPID MAPS"
    if _RE_DB.match(q):
        return "DrugBank"
    if _cas_ok(q):
        return "CAS"
    if q.isdigit():
        return "PubChem CID"
    if q and " " not in q:
        RDLogger.DisableLog("rdApp.*")  # 名称会被当 SMILES 试解析，别刷屏
        try:
            mol = Chem.MolFromSmiles(q)
        finally:
            RDLogger.EnableLog("rdApp.*")
        if mol is not None and (len(q) <= 3 or re.search(r"[()=#@\[\]\d]", q)):
            return "SMILES"
    return "Name"


# ---------------------------------------------------------------- network + cache

_last_call = 0.0
_sources_cache: Optional[dict] = None


def _cache_path() -> Path:
    base = os.environ.get("MZDESK_CACHE_DIR") or os.path.join(os.path.expanduser("~"), ".cache", "mzdesk")
    return Path(base) / "identity.sqlite"


def _cache_open():
    """Open the cache db; None if unwritable (cache is optional)."""
    try:
        p = _cache_path()
        p.parent.mkdir(parents=True, exist_ok=True)
        con = sqlite3.connect(str(p), timeout=5)
        con.execute("CREATE TABLE IF NOT EXISTS cache (k TEXT PRIMARY KEY, v TEXT, t REAL)")
        return con
    except Exception:
        return None


def _cache_get(key: str) -> Optional[str]:
    con = _cache_open()
    if con is None:
        return None
    try:
        row = con.execute("SELECT v, t FROM cache WHERE k=?", (key,)).fetchone()
        if row and time.time() - row[1] < CACHE_TTL:
            return row[0]
    except Exception:
        pass
    finally:
        con.close()
    return None


def _cache_put(key: str, value: str) -> None:
    con = _cache_open()
    if con is None:
        return
    try:
        with con:
            con.execute("INSERT OR REPLACE INTO cache VALUES (?,?,?)", (key, value, time.time()))
    except Exception:
        pass
    finally:
        con.close()


def _http(url: str, data: Optional[bytes] = None, timeout: float = 15.0) -> str:
    """The only network entry point. GET, or POST when data is given (JSON if it starts with '{')."""
    global _last_call
    key = url + "\n" + (data.decode("utf-8", "replace") if data else "")
    hit = _cache_get(key)
    if hit is not None:
        return hit
    wait = _MIN_INTERVAL - (time.time() - _last_call)
    if wait > 0:
        time.sleep(wait)
    headers = {"User-Agent": "mzDesk", "Accept": "application/json, text/plain"}
    if data:
        if "reactome.org/AnalysisService/" in url:
            headers["Content-Type"] = "text/plain"
        else:
            headers["Content-Type"] = ("application/json" if data[:1] == b"{"
                                       else "application/x-www-form-urlencoded")
    req = urllib.request.Request(url, data=data, headers=headers)
    try:
        _last_call = time.time()
        with urllib.request.urlopen(req, timeout=timeout) as r:
            text = r.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as e:
        if e.code == 404:
            raise _NotFound(url) from e
        raise ResolveError(f"Network lookup failed (PubChem/UniChem unreachable): HTTP {e.code}") from e
    except (urllib.error.URLError, socket.timeout, TimeoutError, OSError) as e:
        raise ResolveError(f"Network lookup failed (PubChem/UniChem unreachable): {e}") from e
    _cache_put(key, text)
    return text


def _get_text(url: str, data: Optional[bytes] = None, timeout: float = 15.0) -> str:
    return _http(url, data, timeout)


def _get_json(url: str, data: Optional[bytes] = None, timeout: float = 15.0):
    text = _http(url, data, timeout)
    try:
        return json.loads(text)
    except ValueError as e:
        raise ResolveError(f"Unexpected response from {url.split('/')[2]}.") from e


# ---------------------------------------------------------------- structure helpers

def _mol_from(query: str, kind: str):
    if kind == "InChI":
        mol = Chem.MolFromInchi(query)
        if mol is None:
            raise ResolveError("The InChI could not be parsed.")
        return mol
    mol = Chem.MolFromSmiles(query)
    if mol is None:
        raise ResolveError("The SMILES could not be parsed.")
    return mol


def _fill_structure(res: Resolved, mol) -> None:
    res.smiles = Chem.MolToSmiles(mol)
    try:
        res.inchi = rd_inchi.MolToInchi(mol)
        res.inchikey = rd_inchi.InchiToInchiKey(res.inchi) if res.inchi else None
    except Exception:
        pass


def _fill_links(res: Resolved) -> None:
    if res.cid:
        res.links["PubChem"] = _LINKS["PubChem"].format(res.cid)
    for label, vals in res.ids.items():
        if vals and label in _LINKS:
            res.links[label] = _LINKS[label].format(vals[0])


# ---------------------------------------------------------------- PubChem

def _pc_cids(path: str, data: Optional[bytes], query: str, timeout: float) -> list:
    try:
        text = _get_text(f"{PUG}/{path}", data, timeout)
    except _NotFound:
        raise ResolveError(f"'{query}' not found in PubChem.")
    cids = [int(x) for x in text.split() if x.isdigit() and int(x) > 0]
    if not cids:
        raise ResolveError(f"'{query}' not found in PubChem.")
    return cids


def _pc_props(cid: int, timeout: float) -> dict:
    fields = "Title,IUPACName,SMILES,ConnectivitySMILES,IsomericSMILES,CanonicalSMILES,InChI,InChIKey"
    try:
        js = _get_json(f"{PUG}/compound/cid/{cid}/property/{fields}/JSON", None, timeout)
    except _NotFound:
        raise ResolveError(f"PubChem CID {cid} not found.")
    try:
        return js["PropertyTable"]["Properties"][0]
    except (KeyError, IndexError, TypeError):
        raise ResolveError(f"PubChem CID {cid} not found.")


def _pc_title(cid: int, timeout: float) -> str:
    try:
        js = _get_json(f"{PUG}/compound/cid/{cid}/property/Title/JSON", None, timeout)
        return js["PropertyTable"]["Properties"][0].get("Title") or ""
    except (_NotFound, ResolveError, KeyError, IndexError, TypeError):
        return ""


def _pc_synonyms(cid: int, timeout: float) -> list:
    try:
        text = _get_text(f"{PUG}/compound/cid/{cid}/synonyms/TXT", None, timeout)
    except _NotFound:
        return []
    return [s.strip() for s in text.splitlines() if s.strip()]


def _is_id(s: str) -> bool:
    return bool(_cas_ok(s) or _RE_KEGG.match(s) or _RE_CHEBI.match(s) or _RE_HMDB.match(s)
                or _RE_DB.match(s) or re.match(r"^LM[A-Z]{2}\d{8,}$", s))


def _ids_from_synonyms(syns: list) -> dict:
    ids: dict = {}

    def add(label, v):
        lst = ids.setdefault(label, [])
        if v not in lst:
            lst.append(v)

    for s in syns:
        if _cas_ok(s):
            add("CAS", s)
        elif _RE_KEGG.match(s):
            add("KEGG", s)
        elif _RE_CHEBI.match(s):
            add("ChEBI", "CHEBI:" + s.split(":")[1])
        elif _RE_HMDB.match(s):
            add("HMDB", s)
        elif _RE_DB.match(s):
            add("DrugBank", s)
        elif re.match(r"^LM[A-Z]{2}\d{8,}$", s):
            add("LIPID MAPS", s)
    return ids


# ---------------------------------------------------------------- UniChem

def _unichem_sources(timeout: float) -> dict:
    """shortName(name) -> sourceID, fetched once."""
    global _sources_cache
    if _sources_cache is None:
        js = _get_json(f"{UNICHEM}/sources", None, timeout)
        out = {}
        for s in js.get("sources", []):
            nm = s.get("name") or s.get("shortName")
            if nm and s.get("sourceID") is not None:
                out[str(nm).lower()] = s["sourceID"]
        _sources_cache = out
    return _sources_cache


def _unichem_post(body: dict, timeout: float):
    try:
        return _get_json(f"{UNICHEM}/compounds", json.dumps(body).encode(), timeout)
    except _NotFound:
        return {}


def _unichem_by_key(key: str, timeout: float) -> list:
    js = _unichem_post({"type": "inchikey", "compound": key}, timeout)
    comps = js.get("compounds") or []
    return comps[0].get("sources", []) if comps else []


def _merge_unichem(sources: list, ids: dict) -> Optional[int]:
    """Merge UniChem sources into ids; return a PubChem CID if present."""
    cid = None
    for s in sources:
        short = str(s.get("shortName", "")).lower()
        val = str(s.get("compoundId", "")).strip()
        if not val:
            continue
        if short == "pubchem":
            if val.isdigit() and cid is None:
                cid = int(val)
            continue
        label = _UNICHEM_MAP.get(short)
        if not label:
            continue
        if label == "ChEBI" and not val.upper().startswith("CHEBI:"):
            val = "CHEBI:" + val
        lst = ids.setdefault(label, [])
        if val not in lst:
            lst.append(val)
    return cid


def _key_from_source_id(kind: str, query: str, timeout: float) -> Optional[str]:
    src = _unichem_sources(timeout).get(_KIND_SOURCE[kind])
    if src is None:
        return None
    js = _unichem_post({"type": "sourceID", "compound": query, "sourceID": src}, timeout)
    comps = js.get("compounds") or []
    return comps[0].get("standardInchiKey") if comps else None


# ---------------------------------------------------------------- resolve

def _norm_query(q: str, kind: str) -> str:
    if kind == "ChEBI" and q.lower().startswith("chebi:"):
        return "CHEBI:" + q.split(":", 1)[1]
    if kind in ("HMDB", "DrugBank", "LIPID MAPS", "KEGG"):
        return q.upper()
    return q


def _find_cids(q: str, kind: str, timeout: float) -> list:
    """Query -> PubChem CIDs for non-structure kinds (first one is used)."""
    if kind == "PubChem CID":
        return [int(q)]
    if kind == "InChIKey":
        return _pc_cids(f"compound/inchikey/{urllib.parse.quote(q)}/cids/TXT", None, q, timeout)
    try:
        return _pc_cids(f"compound/name/{urllib.parse.quote(q, safe='')}/cids/TXT", None, q, timeout)
    except ResolveError as first:
        if kind in _KIND_SOURCE:  # 数据库 ID：PubChem 没认出来就走 UniChem
            key = _key_from_source_id(kind, q, timeout)
            if key:
                return _pc_cids(f"compound/inchikey/{key}/cids/TXT", None, q, timeout)
        raise first


def resolve(query: str, kind: str = "Auto", *, want_ids: bool = True, timeout: float = 15.0) -> Resolved:
    q = (query or "").strip()
    if not q:
        raise ResolveError("Enter a name, ID, SMILES, InChI or InChIKey.")
    if kind not in KINDS:
        raise ResolveError(f"Unknown input type '{kind}'.")
    if kind == "Auto":
        kind = detect(q)
    q = _norm_query(q, kind)
    res = Resolved(query=query.strip(), kind=kind)

    cid: Optional[int] = None
    if kind in ("SMILES", "InChI"):
        mol = _mol_from(q, kind)
        _fill_structure(res, mol)
        if not want_ids:
            return res
        try:
            if kind == "SMILES":
                cids = _pc_cids("compound/smiles/cids/TXT",
                                urllib.parse.urlencode({"smiles": q}).encode(), q, timeout)
            else:
                cids = _pc_cids("compound/inchi/cids/TXT",
                                urllib.parse.urlencode({"inchi": q}).encode(), q, timeout)
            cid = cids[0]
        except ResolveError as e:
            res.warnings.append(f"{e} Showing the structure only.")
            return res
    else:
        cids = _find_cids(q, kind, timeout)
        cid = cids[0]
        if len(cids) > 1 and kind not in ("InChIKey", "PubChem CID"):
            res.candidates = [(c, _pc_title(c, timeout)) for c in cids[:MAX_CANDIDATES]]
            others = ", ".join(f"{c} ({t})" for c, t in res.candidates[1:])
            res.warnings.append(f"{len(cids)} PubChem compounds match '{q}'; showing CID {cid} "
                                f"({res.candidates[0][1]}). Others: {others}")

    props = _pc_props(cid, timeout)
    res.cid = cid
    res.name = props.get("Title") or None
    res.iupac = props.get("IUPACName") or None
    if res.smiles is None:
        smi = next((props.get(k) for k in ("SMILES", "IsomericSMILES", "ConnectivitySMILES", "CanonicalSMILES")
                    if props.get(k)), None)
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is not None:
            _fill_structure(res, mol)
        else:
            res.smiles = smi
            res.inchi = props.get("InChI")
            res.inchikey = props.get("InChIKey")
            if smi:
                res.warnings.append("RDKit could not parse the PubChem structure; using PubChem's strings.")
    if not res.inchikey:
        res.inchikey = props.get("InChIKey")

    if want_ids:
        _fill_ids(res, timeout)
    _fill_links(res)
    return res


def _fill_ids(res: Resolved, timeout: float) -> None:
    syns = _pc_synonyms(res.cid, timeout)
    res.ids = _ids_from_synonyms(syns)
    res.synonyms = [s for s in syns if not _is_id(s)][:MAX_SYNONYMS]
    if res.inchikey:
        try:
            sources = _unichem_by_key(res.inchikey, timeout)
            if not sources and res.smiles:
                m = Chem.MolFromSmiles(res.smiles)
                if m is not None:
                    Chem.RemoveStereochemistry(m)
                    flat = rd_inchi.MolToInchiKey(m)
                    if flat and flat != res.inchikey:
                        sources = _unichem_by_key(flat, timeout)
                        if sources:
                            res.stereo_fallback = True
                            res.warnings.append("Cross-database IDs matched by structure without stereochemistry "
                                                "(InChIKey first block); they may refer to a stereoisomer.")
            cid = _merge_unichem(sources, res.ids)
            if res.cid is None and cid:
                res.cid = cid
        except ResolveError as e:
            res.warnings.append(f"UniChem lookup failed: {e}")


def to_smiles(query: str, kind: str = "Auto") -> str:
    """Structure only: SMILES/InChI are parsed offline, everything else needs PubChem."""
    q = (query or "").strip()
    if kind == "Auto":
        kind = detect(q)
    r = resolve(q, kind, want_ids=False)
    if not r.smiles:
        raise ResolveError(f"No structure found for '{query}'.")
    return r.smiles
