"""Identity annotations from independent public sources and an optional local HMDB XML.

Every remote source is optional: a failure adds a warning and leaves the other
annotations intact. Rows describe source records, not a unified truth claim.
All network access goes through identity.resolve's cached HTTP entry point.
"""
from __future__ import annotations

import os
import re
import time
import urllib.parse
import xml.etree.ElementTree as ET

from mzdesk.identity import resolve as resolver

OUTPUT_KEYS = ("classification", "origin", "pathways", "biofluids", "species", "standards", "ghs")
MAX_SPECIES = 50
MAX_IDS = 3
PUBCHEM_VIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{}/JSON"
LOTUS_SEARCH = "https://lotus.naturalproducts.net/api/search/simple?query={}"
NPCLASSIFIER = "https://npclassifier.gnps2.org/classify?smiles={}"
CLASSYFIRE = "https://api.naturalproducts.net/latest/chem/classyfire/classify?smiles={}"
CLASSYFIRE_RESULT = "https://api.naturalproducts.net/latest/chem/classyfire/{}/result"
PUBCHEM_CATEGORIES = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/categories/compound/{}/JSON"
PUBCHEM_SAME_CONNECTIVITY = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastidentity/cid/{}/cids/JSON?identity_type=same_connectivity"
MAX_ISOMERS = 5
# Vendors that sell catalog products (standards houses included); case-insensitive substring of the PubChem source name.
CATALOG_VENDORS = (
    "sigma-aldrich", "merck", "supelco", "tci", "thermo fisher", "alfa aesar", "acros", "cayman", "avanti", "larodan",
    "lgc standards", "accustandard", "toronto research", "santa cruz", "matreya", "medchemexpress", "selleck", "biosynth",
    "chem-impex", "spectrum chemical", "fujifilm", "wako", "combi-blocks", "ambeed", "bld pharm", "apollo scientific",
    "oakwood", "fluorochem", "ak scientific", "abcam", "tocris", "glentham", "carbosynth", "hpc standards", "chiron",
    "dr. ehrenstorfer", "cambridge isotope", "cdn isotopes", "clearsynth", "aablocks", "angene", "j&k", "aladdin",
    "macklin", "nacalai", "kanto")
# Aggregators, screening libraries and make-on-demand sources: a listing is not stock.
SCREENING_VENDORS = (
    "molport", "zinc", "mcule", "chemspace", "enamine", "aurora fine", "vitas-m", "chemdiv", "otava", "life chemicals",
    "akos", "ambinter", "princeton biomolecular", "specs", "chembridge", "asinex", "interbioscreen", "key organics",
    "uorsy", "chemfaces", "aronis", "zerenex", "emolecules", "abamachem", "chemieliva")
REACTOME = ("https://reactome.org/AnalysisService/identifiers/projection?"
            "pageSize=100&page=1&interactors=false&includeDisease=true")


def empty_outputs():
    """Fresh empty blocks, including sources that were not queried."""
    return {key: {"rows": [], "count": 0} for key in OUTPUT_KEYS}


def _link(text, url=None, help=None):
    cell = {"text": str(text)}
    if url:
        cell["url"] = url
    if help:
        cell["help"] = help
    return cell


def _clean(value):
    return re.sub(r"\s+", " ", str(value or "")).strip()


def _term(value):
    """Text from the string / named-object forms used by the source APIs."""
    if isinstance(value, dict):
        for key in ("name", "label", "value", "text", "term"):
            if value.get(key) not in (None, ""):
                return _term(value[key])
        return ""
    if isinstance(value, (list, tuple)):
        return "; ".join(dict.fromkeys(filter(None, (_term(v) for v in value))))
    text = _clean(value)
    if "|" in text:                              # LOTUS may repeat a term once per taxonomy provider
        return "; ".join(dict.fromkeys(part for part in (_clean(v) for v in text.split("|")) if part))
    return text


def _items(value):
    return list(value) if isinstance(value, (list, tuple)) else ([] if value in (None, "") else [value])


def _add(rows, row):
    if row not in rows:
        rows.append(row)


def _attempt(label, function, warnings, default=None):
    try:
        return function()
    except Exception as exc:                       # one optional source must not fail the whole Identity run
        message = _clean(exc) or exc.__class__.__name__
        warnings.append(f"{label} lookup failed: {message}")
        return default


def _local_classification(smiles):
    from mzdesk.platform_fit import display_en
    from mzdesk.platform_fit import lipid_metabolite_classifier as classifier

    result = classifier.analyze(smiles)
    category = display_en.category_en(result.category)
    rows = [["mzDesk rules", "Workflow", category]]
    if result.features.lipid_class:
        rows.append(["mzDesk rules", "Lipid class", display_en.lipid_class_en(result.features)])
    return rows, category in ("Lipidome", "Both")


# ---------------------------------------------------------------- PubChem PUG View

def _pc_references(record):
    return {str(ref.get("ReferenceNumber")): ref for ref in record.get("Reference", [])
            if ref.get("ReferenceNumber") is not None}


def _pc_info_refs(info, references):
    numbers = _items(info.get("ReferenceNumber"))
    refs = [references.get(str(number)) for number in numbers]
    refs = [ref for ref in refs if ref]
    return refs or [{}]


def _pc_texts(info):
    value = info.get("Value") or {}
    out = []
    for item in value.get("StringWithMarkup", []):
        text = _clean(item.get("String"))
        if text:
            urls = [m.get("URL") for m in item.get("Markup", []) if m.get("URL")]
            out.append((text, urls))
    for text in _items(value.get("String")):
        if _clean(text):
            out.append((_clean(text), []))
    return out


def _pc_entries(payload):
    record = payload.get("Record") or {}
    references = _pc_references(record)

    def walk(section, path=()):
        heading = _clean(section.get("TOCHeading"))
        current = path + ((heading,) if heading else ())
        for info in section.get("Information", []):
            yield {"path": current, "name": _clean(info.get("Name")),
                   "texts": _pc_texts(info), "refs": _pc_info_refs(info, references)}
        for child in section.get("Section", []):
            yield from walk(child, current)

    for section in record.get("Section", []):
        yield from walk(section)


def _pc_source(ref):
    return _clean(ref.get("SourceName") or ref.get("Name")) or "PubChem", ref.get("URL")


def _pubchem_annotations(res, timeout):
    cid = res.cid
    payload = resolver._get_json(PUBCHEM_VIEW.format(cid), None, timeout)
    entries = list(_pc_entries(payload))
    compound_url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
    origin, pathways, standards, ghs = [], [], [], []
    sections = {heading for entry in entries for heading in entry["path"]}

    origin_specs = [
        ("Drug and Medication Information", "Drug", "PubChem Drug and Medication Information record"),
        ("Food Additives and Ingredients", "Food", "PubChem Food Additives and Ingredients record"),
    ]
    for heading, label, evidence in origin_specs:
        if heading in sections:
            _add(origin, [label, evidence, _link("PubChem", compound_url)])
    environmental = next((h for h in sections if "environment" in h.lower()
                          and ("contamin" in h.lower() or "pollut" in h.lower())), None)
    if environmental:
        _add(origin, ["Environmental pollutant", f"PubChem {environmental} record", _link("PubChem", compound_url)])

    reports = {}
    for entry in entries:
        path = entry["path"]
        if "Human Metabolite Information" in path and "Metabolite Pathways" in path:
            for text, urls in entry["texts"]:
                url = next((u for u in urls if "smpdb" in u.lower() or "hmdb" in u.lower()), None)
                url = url or (urls[0] if urls else None)
                _add(pathways, [text, "HMDB / SMPDB", _link("Open", url) if url else "—"])
        if "GHS Classification" not in path:
            continue
        kind = entry["name"].lower()
        is_signal = kind == "signal"
        is_hazard = "hazard statement" in kind
        if not (is_signal or is_hazard):
            continue
        for ref in entry["refs"]:
            source, url = _pc_source(ref)
            key = str(ref.get("ReferenceNumber") or source)
            report = reports.setdefault(key, {"source": _link(source, url) if url else source,
                                              "signal": [], "hazards": []})
            field = "signal" if is_signal else "hazards"
            for text, _ in entry["texts"]:
                if text not in report[field]:
                    report[field].append(text)
    for report in reports.values():
        signal = "; ".join(report["signal"]) or "—"
        hazards = "; ".join(report["hazards"])
        if signal != "—" or hazards:
            ghs.append([report["source"], signal, hazards or "—"])
    return {"origin": origin, "pathways": pathways, "standards": standards, "ghs": ghs}


# ---------------------------------------------------------------- commercial standards

def _vendor_names(cid, timeout):
    try:
        payload = resolver._get_json(PUBCHEM_CATEGORIES.format(cid), None, timeout)
    except resolver._NotFound:
        return []
    for category in (payload or {}).get("SourceCategories", {}).get("Categories", []):
        if category.get("Category") == "Chemical Vendors":
            return sorted({_clean(s.get("SourceName")) for s in category.get("Sources", []) if s.get("SourceName")})
    return []


def _split_vendors(names):
    catalog, screening, other = [], [], []
    for name in names:
        low = name.lower()
        if any(v in low for v in CATALOG_VENDORS):
            catalog.append(name)
        elif any(v in low for v in SCREENING_VENDORS):
            screening.append(name)
        else:
            other.append(name)
    return catalog, screening, other


def _standards_row(res, timeout):
    """One row: PubChem vendor listings, a grade and the vendor counts.
    A: a catalog vendor lists this structure. B: only screening / other vendors. B-iso: none for this structure, but a catalog
    vendor lists a stereoisomer / isotopologue. C: no vendor listing. A listing is not a certified standard."""
    catalog, screening, other = _split_vendors(_vendor_names(res.cid, timeout))
    isomers = []
    if not catalog:
        try:
            payload = resolver._get_json(PUBCHEM_SAME_CONNECTIVITY.format(res.cid), None, timeout)
            cids = [c for c in (payload or {}).get("IdentifierList", {}).get("CID", []) if c and c != res.cid]
        except resolver._NotFound:
            cids = []
        for cid in cids[:MAX_ISOMERS]:
            hit, _, _ = _split_vendors(_vendor_names(cid, timeout))
            if hit:
                isomers.append(f"CID {cid}")
    if catalog:
        grade = "A · catalog vendor lists this structure"
    elif screening or other:
        grade = "B · screening or other vendors only"
    elif isomers:
        grade = "B-iso · catalog vendor lists an isomer"
    else:
        grade = "C · no vendor listing"
    counts = f"Catalog {len(catalog)} · Screening {len(screening)} · Other {len(other)}"
    if isomers:
        counts += " · Isomer catalog hits: " + ", ".join(isomers)
    url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{res.cid}#section=Chemical-Vendors"
    return [_link("PubChem vendor listings", url), grade, counts]


# ---------------------------------------------------------------- LOTUS and structure classifiers

def _lotus_record(inchikey, timeout):
    url = LOTUS_SEARCH.format(urllib.parse.quote(inchikey, safe=""))
    payload = resolver._get_json(url, None, timeout)
    products = payload.get("naturalProducts", []) if isinstance(payload, dict) else payload
    return {"record": products[0] if products else {}, "url": url}


def _taxonomy_dicts(value):
    if isinstance(value, dict):
        if value.get("species"):
            yield value
        for child in value.values():
            if isinstance(child, (dict, list, tuple)):
                yield from _taxonomy_dicts(child)
    elif isinstance(value, (list, tuple)):
        for child in value:
            yield from _taxonomy_dicts(child)


def _lotus_species_origin(lotus):
    record, url = lotus["record"], lotus["url"]
    taxa = list(_taxonomy_dicts(record.get("taxonomyReferenceObjects", [])))
    species = {_clean(_term(t.get("species"))) for t in taxa}
    if not species:
        species = {_clean(_term(t)) for t in _items(record.get("allTaxa"))}
    species.discard("")
    species_rows = [[name, _link("LOTUS", url)] for name in sorted(species, key=str.casefold)[:MAX_SPECIES]]

    origin = []
    kingdoms = {_term(t.get("kingdom")) for t in taxa if _term(t.get("kingdom"))}
    for kingdom in sorted(kingdoms):
        low = kingdom.casefold()
        if low in ("plant", "plantae", "plants"):
            label = "Plant"
        elif low in ("bacteria", "archaea", "fungi", "microbial", "microorganism"):
            label = "Microbial"
        else:
            continue
        _add(origin, [label, f"LOTUS occurrence record; kingdom: {kingdom}", _link("LOTUS", url)])
    return species_rows, origin


def _lotus_classification(record):
    np_levels = (("Pathway", "chemicalTaxonomyNPclassifierPathway"),
                 ("Superclass", "chemicalTaxonomyNPclassifierSuperclass"),
                 ("Class", "chemicalTaxonomyNPclassifierClass"))
    cf_levels = (("Kingdom", "chemicalTaxonomyClassyfireKingdom"),
                 ("Superclass", "chemicalTaxonomyClassyfireSuperclass"),
                 ("Class", "chemicalTaxonomyClassyfireClass"),
                 ("Direct parent", "chemicalTaxonomyClassyfireDirectParent"))
    np_rows = [["LOTUS / NPClassifier", level, _term(record.get(key))] for level, key in np_levels if _term(record.get(key))]
    cf_rows = [["LOTUS / ClassyFire", level, _term(record.get(key))] for level, key in cf_levels if _term(record.get(key))]
    return np_rows, cf_rows


def _npclassification(smiles, timeout, present):
    url = NPCLASSIFIER.format(urllib.parse.quote(smiles, safe=""))
    payload = resolver._get_json(url, None, timeout)
    levels = (("Pathway", "pathway_results"), ("Superclass", "superclass_results"), ("Class", "class_results"))
    return [["NPClassifier", level, _term(term)] for level, key in levels if level not in present
            for term in _items(payload.get(key)) if _term(term)]


def _normal_key(key):
    return re.sub(r"[^a-z]", "", str(key).lower())


def _classy_dicts(value):
    if isinstance(value, dict):
        keys = {_normal_key(key) for key in value}
        if keys & {"kingdom", "superclass", "class", "subclass", "directparent"}:
            yield value
        for child in value.values():
            if isinstance(child, (dict, list, tuple)):
                yield from _classy_dicts(child)
    elif isinstance(value, (list, tuple)):
        for child in value:
            yield from _classy_dicts(child)


def _classy_rows(payload, present=()):
    levels = (("Kingdom", "kingdom"), ("Superclass", "superclass"), ("Class", "class"),
              ("Subclass", "subclass"), ("Direct parent", "directparent"))
    rows = []
    for candidate in _classy_dicts(payload):
        values = {_normal_key(key): value for key, value in candidate.items()}
        for level, key in levels:
            term = _term(values.get(key))
            if term and level not in present:
                _add(rows, ["ClassyFire", level, term])
        if rows:
            break
    return rows


def _classyfire(smiles, timeout, present):
    payload = resolver._get_json(CLASSYFIRE.format(urllib.parse.quote(smiles, safe="")), None, timeout)
    rows = _classy_rows(payload, present)
    if rows:
        return rows
    job = payload.get("id") if isinstance(payload, dict) else None
    if job is None:
        return []
    for attempt in range(3):
        time.sleep(0.5)
        payload = resolver._get_json(CLASSYFIRE_RESULT.format(urllib.parse.quote(str(job))), None, timeout)
        rows = _classy_rows(payload, present)
        if rows:
            return rows
    raise LookupError("classification job was not ready after three polls")


def _lipidmaps(res, timeout):
    lm_ids = res.ids.get("LIPID MAPS", [])
    if lm_ids:
        url = "https://www.lipidmaps.org/rest/compound/lm_id/{}/classification".format(
            urllib.parse.quote(str(lm_ids[0]), safe=""))
    else:
        url = f"https://www.lipidmaps.org/rest/compound/pubchem_cid/{res.cid}/classification"
    payload = resolver._get_json(url, None, timeout)
    candidates = [payload]
    if isinstance(payload, dict):
        candidates.extend(payload.values())
    for candidate in candidates:
        if not isinstance(candidate, dict):
            continue
        values = {_normal_key(key): value for key, value in candidate.items()}
        rows = []
        for level, key in (("Category", "core"), ("Main class", "mainclass"), ("Subclass", "subclass")):
            term = _term(values.get(key))
            if term:
                rows.append(["LIPID MAPS", level, term])
        if rows:
            return rows
    return []


# ---------------------------------------------------------------- Pathways

def _kegg_pathways(kegg_ids, timeout, warnings):
    map_ids = []
    for kegg_id in kegg_ids[:MAX_IDS]:
        text = _attempt("KEGG", lambda kid=kegg_id: resolver._get_text(
            f"https://rest.kegg.jp/link/pathway/{urllib.parse.quote(kid)}", None, timeout), warnings, "") or ""
        for line in text.splitlines():
            fields = line.strip().split("\t")
            if len(fields) == 2 and re.fullmatch(r"path:map\d{5}", fields[1]):
                map_id = fields[1].split(":", 1)[1]
                if map_id not in map_ids:
                    map_ids.append(map_id)
    if not map_ids:
        return []
    listing = _attempt("KEGG pathway names", lambda: resolver._get_text(
        "https://rest.kegg.jp/list/pathway", None, timeout), warnings, "") or ""
    names = {}
    for line in listing.splitlines():
        fields = line.strip().split("\t", 1)
        if len(fields) == 2:
            map_id = fields[0].split(":", 1)[-1]       # KEGG currently returns "map..."; older responses may use "path:map..."
            if re.fullmatch(r"map\d{5}", map_id):
                names[map_id] = _clean(fields[1])
    return [[names.get(map_id) or map_id, "KEGG", _link(map_id, f"https://www.kegg.jp/pathway/{map_id}")]
            for map_id in map_ids]


def _reactome_pathways(chebi_ids, timeout, warnings):
    rows = []
    for chebi in chebi_ids[:MAX_IDS]:
        numeric = str(chebi).split(":", 1)[-1]
        if not numeric.isdigit():
            continue
        payload = _attempt("Reactome", lambda body=numeric: resolver._get_json(
            REACTOME, body.encode(), timeout), warnings, {}) or {}
        for pathway in payload.get("pathways", []):
            stable_id = _clean(pathway.get("stId"))
            name = _clean(pathway.get("name"))
            if stable_id and name:
                _add(rows, [name, "Reactome", _link(stable_id, f"https://reactome.org/content/detail/{stable_id}")])
    return rows


# ---------------------------------------------------------------- Optional local HMDB XML

def _tag(element):
    return element.tag.rsplit("}", 1)[-1]


def _children(element, name):
    return [child for child in element if _tag(child) == name]


def _child(element, name):
    return next((child for child in element if _tag(child) == name), None)


def _child_text(element, name):
    child = _child(element, name)
    return _clean(child.text) if child is not None and child.text else ""


def _hmdb_concentrations(path, accessions, links):
    wanted = set(accessions)
    rows, found = [], set()
    root = None
    for event, element in ET.iterparse(path, events=("start", "end")):
        if root is None and event == "start":
            root = element
        if event != "end" or _tag(element) != "metabolite":
            continue
        accession = _child_text(element, "accession")
        if accession in wanted:
            found.add(accession)
            source_url = links.get("HMDB")
            container = _child(element, "normal_concentrations")
            for concentration in _children(container, "normal_concentration") if container is not None else []:
                value = _child_text(concentration, "concentration_value")
                if not value:
                    continue
                references = []
                refs = _child(concentration, "references")
                for reference in _children(refs, "reference") if refs is not None else []:
                    text = _child_text(reference, "reference_text")
                    pmid = _child_text(reference, "pubmed_id")
                    if text:
                        references.append(text)
                    elif pmid:
                        references.append(f"PMID: {pmid}")
                source = _link("HMDB XML", source_url, references or None) if source_url else "HMDB XML"
                rows.append([_child_text(concentration, "biospecimen") or "—", value,
                             _child_text(concentration, "concentration_units") or "—", source])
        element.clear()
        if root is not None:
            root.clear()
        if found == wanted:
            break
    return rows


# ---------------------------------------------------------------- Top level

def enrich(res, *, timeout=12.0):
    """Return all seven annotation blocks plus non-fatal source warnings."""
    outputs = empty_outputs()
    warnings = []

    local_rows, is_lipid = _attempt("mzDesk classification", lambda: _local_classification(res.smiles), warnings, ([], False))
    outputs["classification"]["rows"].extend(local_rows)

    lotus = None
    if res.inchikey:
        lotus = _attempt("LOTUS", lambda: _lotus_record(res.inchikey, timeout), warnings)
    np_rows, cf_rows = [], []
    if lotus and lotus["record"]:
        np_rows, cf_rows = _lotus_classification(lotus["record"])
        species, origin = _lotus_species_origin(lotus)
        outputs["species"]["rows"].extend(species)
        outputs["origin"]["rows"].extend(origin)
    if res.smiles:
        np_present = {row[1] for row in np_rows}
        if len(np_present) < 3:
            extra = _attempt("NPClassifier", lambda: _npclassification(res.smiles, timeout, np_present), warnings, []) or []
            np_rows.extend(extra)
        cf_present = {row[1] for row in cf_rows}
        if len(cf_present) < 4:
            extra = _attempt("ClassyFire", lambda: _classyfire(res.smiles, timeout, cf_present), warnings, []) or []
            cf_rows.extend(extra)
    outputs["classification"]["rows"].extend(np_rows + cf_rows)

    lm_ids = res.ids.get("LIPID MAPS", [])
    if (lm_ids or is_lipid) and (lm_ids or res.cid):
        rows = _attempt("LIPID MAPS", lambda: _lipidmaps(res, timeout), warnings, []) or []
        outputs["classification"]["rows"].extend(rows)

    if res.cid:
        annotations = _attempt("PubChem annotations", lambda: _pubchem_annotations(res, timeout), warnings, {}) or {}
        for key in ("origin", "pathways", "standards", "ghs"):
            for row in annotations.get(key, []):
                _add(outputs[key]["rows"], row)

    if res.cid:
        row = _attempt("PubChem vendors", lambda: _standards_row(res, timeout), warnings)
        if row:
            _add(outputs["standards"]["rows"], row)

    if res.ids.get("HMDB"):
        hmdb_id = res.ids["HMDB"][0]
        url = res.links.get("HMDB")
        _add(outputs["origin"]["rows"], [
            "Human-associated", "HMDB metabolite record",
            _link(hmdb_id, url) if url else hmdb_id])

    if res.ids.get("KEGG"):
        for row in _kegg_pathways(res.ids["KEGG"], timeout, warnings):
            _add(outputs["pathways"]["rows"], row)
    if res.ids.get("ChEBI"):
        for row in _reactome_pathways(res.ids["ChEBI"], timeout, warnings):
            _add(outputs["pathways"]["rows"], row)

    hmdb_xml = os.environ.get("MZDESK_HMDB_XML")
    if hmdb_xml and res.ids.get("HMDB"):
        rows = _attempt("HMDB XML", lambda: _hmdb_concentrations(hmdb_xml, res.ids["HMDB"], res.links), warnings, []) or []
        outputs["biofluids"]["rows"].extend(rows)

    for key in OUTPUT_KEYS:
        rows = outputs[key]["rows"]
        outputs[key]["count"] = len(rows)
    return {"outputs": outputs, "warnings": warnings}
