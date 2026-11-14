"""A small stand-in for the PubChem calls of the exact-mass and alternative-explanation searches.

FakePubChem(compounds) replaces mzdesk.ms1.exact_mass._http: it answers from an in-memory list, in the reply shapes PubChem uses
(checked against the real service), and records every request in .calls. Esearch terms are read with the same few patterns
the engine writes: mass windows, one [SourceName] and one [DCNT] range; the default order is newest CID first, like NCBI.
fastformula answers by exact formula; xrefs answer RegistryID or SourceName.
"""
import json
import re
import time
import urllib.parse

HMDB, DRUGBANK, CHEBI, KEGG, LIPID, NORMAN = ("Human Metabolome Database (HMDB)", "DrugBank", "ChEBI", "KEGG", "LIPID MAPS",
                                              "NORMAN Suspect List Exchange")      # PubChem source names
KNOWN_SOURCES = {HMDB, DRUGBANK, CHEBI, KEGG, LIPID, NORMAN}


def compound(cid, title, formula, mass, sources=(), dcnt=1, charge=0, regids=(), inchikey=""):
    return {"cid": cid, "title": title, "formula": formula, "mass": mass, "sources": set(sources), "dcnt": dcnt,
            "charge": charge, "regids": list(regids), "inchikey": inchikey}


CAFFEINE = compound(2519, "Caffeine", "C8H10N4O2", 194.08037557, [HMDB, DRUGBANK, CHEBI, KEGG, NORMAN], 300,
                    regids=["58-08-2", "HMDB0001847", "CHEBI:27732", "DB00201", "C07481", "D00528", "D96575", "G12345"] + [f"X{i}" for i in range(40)],
                    inchikey="RYYVLZVUVIJVGH-UHFFFAOYSA-N")
ENPROFYLLINE = compound(1676, "Enprofylline", "C8H10N4O2", 194.08037557, [HMDB, DRUGBANK], 120,
                        regids=["HMDB0014962", "DB00824", "G56441", "D04006"], inchikey="SIQPXVQCUCHWDI-UHFFFAOYSA-N")
ISOCAFFEINE = compound(1326, "Isocaffeine", "C8H10N4O2", 194.08037557, [], 40, inchikey="LPHGQDQBBGAPDZ-UHFFFAOYSA-N")
METFORMIN = compound(4091, "Metformin", "C4H11N5", 129.10144537, [HMDB, DRUGBANK, CHEBI, KEGG], 200,
                     regids=["HMDB0001921", "CHEBI:6801", "DB00331", "C07151", "657-24-9"], inchikey="XZWYZXLIPXDOLR-UHFFFAOYSA-N")
CHOLINE = compound(305, "Choline", "C5H14NO+", 104.107539070, [HMDB, CHEBI], 200, charge=1,
                   regids=["HMDB0000097", "CHEBI:15354"], inchikey="OEYIOHPDSNJKLS-UHFFFAOYSA-N")
PC_34_1 = compound(5497103, "1-Palmitoyl-2-oleoyl-sn-glycero-3-phosphocholine", "C42H82NO8P", 759.57780, [LIPID], 80,
                   regids=["LMGP01010005"], inchikey="WTJKGGKOPKCXLL-UHFFFAOYSA-N")


def default_compounds():
    """What every test sees unless it builds its own fake (conftest): enough for the page example."""
    return [CAFFEINE, ENPROFYLLINE, ISOCAFFEINE, CHOLINE, PC_34_1, METFORMIN]


class FakePubChem:
    def __init__(self, compounds=None, slack=0.0, xref_delay=0.0):
        self.compounds = {c["cid"]: c for c in (default_compounds() if compounds is None else compounds)}
        self.slack = slack              # widen every mass window by this many Da (like PubChem's rounding), to test the exact filter
        self.xref_delay = xref_delay
        self.calls = []                 # (kind, esearch term or PUG url, POST body)
        self.urls = []                  # every full URL requested

    # the engine's single network entry point
    def __call__(self, url, data=None):
        body = urllib.parse.unquote_plus(data.decode()) if data else ""
        self.urls.append(url + ("?" + data.decode() if data else ""))
        if "esearch.fcgi" in url:
            q = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)
            self.calls.append(("esearch", q["term"][0], body))
            return self._esearch(q)
        if "/fastformula/" in url:
            formula = urllib.parse.unquote(re.search(r"/fastformula/([^/]+)/", url).group(1))
            maxrec = int((urllib.parse.parse_qs(urllib.parse.urlparse(url).query).get("MaxRecords") or ["100000"])[0])
            self.calls.append(("fastformula", formula, body))
            ids = sorted(c["cid"] for c in self.compounds.values() if c["formula"] == formula)[:maxrec]
            return json.dumps({"IdentifierList": {"CID": ids}})
        kind = "xrefs" if "/xrefs/" in url else "property"
        self.calls.append((kind, url, body))
        cids = [int(c) for c in re.search(r"cid=([\d,]+)", body).group(1).split(",")]
        known = [self.compounds[c] for c in cids if c in self.compounds]
        if not known:
            return "{}"                 # PUG answers 404, which _http turns into an empty result
        if kind == "xrefs":
            time.sleep(self.xref_delay)
            if "/xrefs/SourceName/" in url:
                return json.dumps({"InformationList": {"Information": [{"CID": c["cid"], "SourceName": sorted(c["sources"])} for c in known]}})
            return json.dumps({"InformationList": {"Information": [{"CID": c["cid"], "RegistryID": c["regids"]} for c in known]}})
        return json.dumps({"PropertyTable": {"Properties": [
            {"CID": c["cid"], "MolecularFormula": c["formula"], "MonoisotopicMass": f"{c['mass']:.8f}", "Title": c["title"],
             **({"Charge": c["charge"]} if c["charge"] else {}), **({"InChIKey": c["inchikey"]} if c["inchikey"] else {})} for c in known]}})

    def _esearch(self, q):
        term = q["term"][0]
        retmax, retstart = int(q["retmax"][0]), int(q.get("retstart", ["0"])[0])
        windows = [(float(a), float(b)) for a, b in re.findall(r"(-?[\d.]+):(-?[\d.]+)\[MonoisotopicMass\]", term)]
        hits = [c for c in self.compounds.values() if any(a - self.slack <= c["mass"] <= b + self.slack for a, b in windows)]
        result = {}
        if (src := re.search(r'"([^"]+)"\[SourceName\]', term)):
            if src.group(1) not in KNOWN_SOURCES:
                result["warninglist"] = {"quotedphrasesnotfound": [f'"{src.group(1)}"[SourceName]'], "outputmessages": ["No items found."]}
            hits = [c for c in hits if src.group(1) in c["sources"]]
        if (dc := re.search(r"(\d+):(\d+)\[DCNT\]", term)):
            hits = [c for c in hits if int(dc.group(1)) <= c["dcnt"] <= int(dc.group(2))]
        ids = sorted((c["cid"] for c in hits), reverse=True)
        page = ids[retstart:retstart + retmax]
        result.update({"count": str(len(ids)), "retmax": str(len(page)), "retstart": str(retstart), "idlist": [str(i) for i in page]})
        return json.dumps({"esearchresult": result})

    @property
    def esearch_terms(self):
        return [c[1] for c in self.calls if c[0] == "esearch"]
