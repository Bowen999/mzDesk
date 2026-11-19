"""Structure mapping: analyte -> PubChem CID -> SMILES / InChIKey / parent (salt-stripped, neutralised) InChIKey.

Order per analyte: CAS RN (monograph header) -> English name -> name with salt / hydrate words stripped -> none.
All PubChem responses are cached in a JSON file so that reruns are offline. Politeness: <= 4 requests / second.
Biologics, polymers and mixtures are not looked up (map_method 'none' with a reason).
"""
from __future__ import annotations

import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
MIN_INTERVAL = 0.25
PROPS = "Title,SMILES,ConnectivitySMILES,InChIKey,MolecularFormula,MolecularWeight"

SALT_WORDS = ["hydrochloride", "hydrobromide", "hydroiodide", "sulfate", "sodium", "potassium", "calcium", "magnesium", "hydrate",
              "anhydrous", "maleate", "mesilate", "mesylate", "besilate", "tartrate", "citrate", "phosphate", "acetate", "fumarate",
              "succinate", "nitrate", "oxalate", "lactate", "gluconate", "dried", "monohydrate", "dihydrate", "trihydrate",
              "hemihydrate", "tosilate", "napadisilate", "bitartrate", "bromide", "chloride", "iodide", "hydrogen", "disodium",
              "dipotassium", "ammonium", "benzoate", "palmitate", "stearate", "propionate", "valerate", "pivoxil", "axetil",
              "cilexetil", "medoxomil", "etabonate", "fine granules", "n-hydrate"]
GENERIC_ANALYTES = {"aconitum monoester alkaloids standard ts": "reagent solution (mixture), not one compound","related substances": "class of impurities, not one compound", "residual solvents": "mixture of solvents",
                    "bhc isomers": "mixture of isomers", "ddt isomers": "mixture of isomers", "total alkaloids": "class of compounds",
                    "total saponin": "class of compounds"}
NO_STRUCTURE = [
    (r"genetical recombination|\(genetical", "biologic (recombinant protein / peptide)"),
    (r"\b(insulin|glucagon|calcitonin|heparin|globulin|interferon|epoetin|filgrastim|lenograstim|celmoleukin|teceleukin|nartograstim|"
     r"vaccine|antitoxin|toxin|toxoid|hyaluronate|gonadotrophin|gonadotropin|somatropin|pancreatin|elcatonin|oxytocin|vasopressin|"
     r"albumin|aprotinin|bacitracin|polymyxin|nystatin|colistin|actinomycin|peptide)\b|"
     r"galactosidase|amylase|lipase|protease|diastase|pepsin|trypsin|lysozyme|kallidinogenase|urokinase|ulinastatin|thrombin|"
     r"hyaluronidase|asparaginase|streptokinase|serratiopeptidase|lactase|invertase|cellulase", "biologic / large biomolecule"),
    (r"dried yeast|creosote|parnaparin|enoxaparin|opium alkaloids|royal jelly|bezoar|bear bile|"
     r"\btotal (?:alkaloid|capsaicin|saponin|flavonoid|curcuminoid)s?\b", "natural product / mixture / class of compounds"),
    (r"macrogol|polyethylene glycol|polysorbate|povidone|dextran|hypromellose|cellulose|carmellose|starch|gelatin|polyvinyl|carboxyvinyl|"
     r"sodium alginate|pectin|chondroitin|hyaluronic|lanolin|beeswax|paraffin|petrolatum|colestimide|colestilan|polystyrene sulfonate|"
     r"\boil\b|essential oil", "polymer / mixture / natural product"),
]


class Cache:
    def __init__(self, path: Path | None, offline: bool = False):
        self.path, self.offline = path, offline
        self.data: dict = {}
        self.dirty = 0
        self.last = 0.0
        self.requests = 0
        self.errors = 0
        if path and Path(path).exists():
            try:
                self.data = json.loads(Path(path).read_text(encoding="utf-8"))
            except ValueError:
                self.data = {}

    def save(self) -> None:
        if self.path and self.dirty:
            tmp = Path(str(self.path) + ".tmp")
            tmp.write_text(json.dumps(self.data, ensure_ascii=False), encoding="utf-8")
            tmp.replace(self.path)
            self.dirty = 0

    def get(self, key: str, fetch):
        """cached value, or fetch() -> value | None (None = not found). Network failures are not cached."""
        if key in self.data:
            return self.data[key]
        if self.offline:
            return "__offline__"
        try:
            val = fetch()
        except _NetError:
            self.errors += 1
            return "__error__"
        self.data[key] = val
        self.dirty += 1
        if self.dirty >= 25:
            self.save()
        return val


class _NetError(Exception):
    pass


def _http_json(cache: Cache, url: str):
    """GET -> parsed JSON; None on 404; raises _NetError on persistent failure. Rate limited."""
    for attempt in range(5):
        wait = MIN_INTERVAL - (time.time() - cache.last)
        if wait > 0:
            time.sleep(wait)
        cache.last = time.time()
        cache.requests += 1
        req = urllib.request.Request(url, headers={"User-Agent": "mzDesk-known-methods/1.0", "Accept": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.loads(r.read().decode("utf-8", "replace"))
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if e.code in (429, 500, 502, 503, 504):
                time.sleep(2 ** attempt)
                continue
            if e.code == 400:
                return None
            raise _NetError(str(e))
        except (urllib.error.URLError, TimeoutError, OSError, ValueError):
            time.sleep(2 ** attempt)
    raise _NetError(url)


def cids_by_name(cache: Cache, name: str):
    key = "name:" + name.lower()
    url = f"{PUG}/compound/name/{urllib.parse.quote(name, safe='')}/cids/JSON"

    def fetch():
        js = _http_json(cache, url)
        return (js or {}).get("IdentifierList", {}).get("CID", []) if js else []
    return cache.get(key, fetch)


def props(cache: Cache, cid: int):
    key = f"cid:{cid}"
    url = f"{PUG}/compound/cid/{cid}/property/{PROPS}/JSON"

    def fetch():
        js = _http_json(cache, url)
        try:
            return js["PropertyTable"]["Properties"][0] if js else None
        except (KeyError, IndexError, TypeError):
            return None
    return cache.get(key, fetch)


def cas_ok(s: str) -> bool:
    m = re.fullmatch(r"(\d{2,7})-(\d{2})-(\d)", s or "")
    if not m:
        return False
    digits = (m.group(1) + m.group(2))[::-1]
    return sum((i + 1) * int(d) for i, d in enumerate(digits)) % 10 == int(m.group(3))


def standardize(smiles: str) -> dict:
    """Largest organic fragment, neutralised -> parent SMILES / InChIKey (None when RDKit cannot handle the molecule)."""
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem.MolStandardize import rdMolStandardize
        RDLogger.DisableLog("rdApp.*")
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"parent_smiles": None, "parent_inchikey": None, "inchikey_calc": None}
        full_key = Chem.MolToInchiKey(mol)
        m = rdMolStandardize.LargestFragmentChooser(preferOrganic=True).choose(mol)
        m = rdMolStandardize.Uncharger().uncharge(m)
        return {"parent_smiles": Chem.MolToSmiles(m), "parent_inchikey": Chem.MolToInchiKey(m), "inchikey_calc": full_key}
    except Exception:
        return {"parent_smiles": None, "parent_inchikey": None, "inchikey_calc": None}


GREEK_SPELL = {"α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta", "ε": "epsilon", "ω": "omega", "μ": "mu", "′": "'"}


def _spell(n: str) -> str:
    for g, w in GREEK_SPELL.items():
        n = n.replace(g, w)
    return n


def name_variants(name: str) -> list[tuple[str, str]]:
    """(query, method) in the order they are tried after the CAS attempt."""
    base = re.sub(r"\s+", " ", name).strip()
    base = re.sub(r"\[\s*", "[", re.sub(r"\s*\]", "]", base)).replace(", ", ",") if "[" in base else base
    out = [(base, "name")]
    cands: list[str] = []
    paren = re.match(r"^(.+?)\s+\((.+)\)$", base)                        # "alpha-bhc (alpha-hexachlorocyclohexane)": both names
    if paren and re.search(r"[a-z]{5,}", paren.group(2)):
        cands += [_spell(paren.group(2)), paren.group(2), _spell(paren.group(1))]
    if re.search(r"[a-z] (?=sulfon|carbox|phosph)", base):                  # "o-toluene sulfonamide"
        cands.append(re.sub(r"(?<=[a-z]) (?=sulfon|carbox|phosph)", "", base))
    no_paren = re.sub(r"\s*\([^)]*\)\s*$", "", base).strip()
    if no_paren and no_paren != base and len(no_paren) > 3:
        cands.append(no_paren)
    spelled = _spell(base)
    if spelled != base:
        cands.append(spelled)
    nostereo = re.sub(r"\((?:e|z|r|s|\+|-|rs|dl)\)-?", "", base, flags=re.I).strip(" -")
    if nostereo != base and len(nostereo) > 3:
        cands.append(nostereo)
    if re.search(r"\[\d+\]-", base):                                      # [6]-shogaol -> 6-shogaol
        cands.append(re.sub(r"\[(\d+)\]-", r"\1-", base))
    if "n-caprylate" in base.lower():
        cands.append(re.sub(r"n-caprylate", "octanoate", base, flags=re.I))
    for c in cands:
        out.append((c, "name-stripped"))
    words = no_paren.split()
    stripped = list(words)
    changed = True
    while changed and len(stripped) > 1:
        changed = False
        if stripped[-1].lower() in SALT_WORDS:
            stripped.pop()
            changed = True
        elif stripped[0].lower() in SALT_WORDS:
            stripped.pop(0)
            changed = True
    s = " ".join(stripped)
    if s and s.lower() != no_paren.lower() and len(s) > 3:
        out.append((s, "name-stripped"))
    return out


def no_structure_reason(name: str, rec: dict) -> str | None:
    n = name.lower()
    if n in GENERIC_ANALYTES:
        return GENERIC_ANALYTES[n]
    for pat, why in NO_STRUCTURE:
        if re.search(pat, n):
            return why
    return None


def map_one(cache: Cache, name: str, cas_list: list[str], rec: dict | None = None) -> dict:
    out = {"name": name, "cid": None, "smiles": None, "inchikey": None, "parent_inchikey": None, "parent_smiles": None,
           "formula": None, "map_method": "none", "pubchem_title": None, "reason": None}
    why = no_structure_reason(name, rec or {})
    if why:
        out["reason"] = why
        return out
    tries: list[tuple[str, str]] = [(c, "cas") for c in cas_list if cas_ok(c)] + name_variants(name)
    seen = set()
    err = False
    for q, method in tries:
        if q.lower() in seen:
            continue
        seen.add(q.lower())
        cids = cids_by_name(cache, q)
        if cids in ("__offline__", "__error__"):
            err = True
            continue
        if not cids:
            continue
        p = props(cache, cids[0])
        if p in ("__offline__", "__error__"):
            err = True
            continue
        if not p:
            continue
        if method == "name-stripped":                                   # derived query: the hit must look related (guards "d-(+)-..." -> daidzein)
            qt = {t for t in re.findall(r"[a-z0-9]{3,}", q.lower())}
            tt = {t for t in re.findall(r"[a-z0-9]{3,}", (p.get("Title") or "").lower())}
            if qt and tt and not (qt & tt) and not any(x in (p.get("Title") or "").lower() for x in qt):
                continue
        smi = p.get("SMILES") or p.get("ConnectivitySMILES")
        std = standardize(smi) if smi else {}
        out.update(cid=int(cids[0]), smiles=smi, inchikey=p.get("InChIKey"), formula=p.get("MolecularFormula"),
                   pubchem_title=p.get("Title"), map_method=method, parent_inchikey=std.get("parent_inchikey"),
                   parent_smiles=std.get("parent_smiles"))
        if len(cids) > 1:
            out["reason"] = f"{len(cids)} PubChem CIDs matched; first used"
        return out
    out["reason"] = "network error / not cached" if err else "not found in PubChem"
    return out


def map_records(records: list[dict], cache: Cache, log=print) -> dict:
    """Add `analyte_structures` (+ primary-analyte convenience fields) to every record. Returns stats."""
    memo: dict[tuple, dict] = {}
    todo = []
    for r in records:
        for a in list(r["analytes"]) + list(r.get("substance_names") or []):
            cas = tuple((r.get("analyte_cas") or {}).get(a, ()))
            todo.append((a, cas))
    uniq = list(dict.fromkeys(todo))
    log(f"mapping {len(uniq)} distinct analyte queries ({len(todo)} analyte slots) ...")
    t0 = time.time()
    for n, (a, cas) in enumerate(uniq, 1):
        memo[(a, cas)] = map_one(cache, a, list(cas))
        if n % 100 == 0:
            log(f"  {n}/{len(uniq)}  ({cache.requests} requests, {time.time() - t0:.0f}s)")
            cache.save()
    cache.save()
    stats = {"unique": len(uniq), "requests": cache.requests, "errors": cache.errors}
    for r in records:
        sts = []
        for i, a in enumerate(r["analytes"]):
            cas = tuple((r.get("analyte_cas") or {}).get(a, ()))
            sts.append(dict(memo[(a, cas)]))
        r["analyte_structures"] = sts
        r["substances"] = [dict(memo[(a, tuple((r.get("analyte_cas") or {}).get(a, ())))]) for a in r.get("substance_names") or []]
        prim = next((s for s in sts if s["cid"]), sts[0] if sts else None)
        # primary = first analyte that is a real compound
        first = sts[0] if sts else None
        prim = first if first and first["cid"] else next((s for s in sts if s["cid"]), first)
        for k in ("cid", "smiles", "inchikey", "parent_inchikey", "map_method"):
            r[k] = prim[k] if prim else (None if k != "map_method" else "none")
        r["map_reason"] = prim["reason"] if prim else "no analyte"
    return stats
