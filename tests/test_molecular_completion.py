"""Offline coverage for the completed Molecular Identity and Properties blocks."""
import json

import pytest

from mzdesk.identity import enrich as E
from mzdesk.identity import resolve as R
from mzdesk.platform_fit import functional_groups as FG


@pytest.fixture(autouse=True)
def isolated_identity(monkeypatch, tmp_path):
    monkeypatch.setenv("MZDESK_CACHE_DIR", str(tmp_path))
    monkeypatch.delenv("MZDESK_HMDB_XML", raising=False)


def resolved(**changes):
    values = dict(
        query="caffeine", kind="Name", smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        inchikey="RYYVLZVUVIJVGH-UHFFFAOYSA-N", cid=2519, name="Caffeine",
        ids={"CAS": ["58-08-2"], "HMDB": ["HMDB0001847"], "KEGG": ["C07481"], "ChEBI": ["CHEBI:27732"]},
        links={"HMDB": "https://hmdb.ca/metabolites/HMDB0001847"}, warnings=[],
    )
    values.update(changes)
    return R.Resolved(**values)


def as_json(payload):
    return json.dumps(payload)


def pc_payload():
    def info(name, text, ref=1, url=None):
        markup = [{"URL": url}] if url else []
        return {"ReferenceNumber": ref, "Name": name,
                "Value": {"StringWithMarkup": [{"String": text, "Markup": markup}]}}

    return {"Record": {
        "Reference": [
            {"ReferenceNumber": 1, "SourceName": "Reporter A", "URL": "https://example.test/a"},
            {"ReferenceNumber": 2, "SourceName": "Reporter B", "URL": "https://example.test/b"},
        ],
        "Section": [
            {"TOCHeading": "Drug and Medication Information", "Information": [info("Drug record", "Caffeine")]},
            {"TOCHeading": "Food Additives and Ingredients", "Information": [info("Food record", "Caffeine")]},
            {"TOCHeading": "Environmental Information", "Section": [
                {"TOCHeading": "Environmental Contaminants", "Information": [info("Occurrence", "Reported")]}]},
            {"TOCHeading": "Human Metabolite Information", "Section": [
                {"TOCHeading": "Metabolite Pathways", "Information": [
                    info("Pathway", "Purine metabolism", url="https://smpdb.ca/view/SMP00001")]}]},
            {"TOCHeading": "Chemical Vendors", "Information": [info("Availability", "Available")]},
            {"TOCHeading": "Safety and Hazards", "Section": [
                {"TOCHeading": "GHS Classification", "Information": [
                    info("Signal", "Danger", 1), info("GHS Hazard Statements", "H301", 1),
                    info("GHS Hazard Statements", "H315", 1), info("Signal", "Warning", 2),
                    info("GHS Hazard Statements", "H302", 2)]}]},
        ],
    }}


def lotus_payload():
    return {"naturalProducts": [{
        "lotus_id": "LTS0001",
        "chemicalTaxonomyNPclassifierPathway": "Alkaloids",
        "chemicalTaxonomyNPclassifierSuperclass": "Pseudoalkaloids",
        "chemicalTaxonomyNPclassifierClass": "Purine alkaloids|Purine alkaloids",
        "chemicalTaxonomyClassyfireKingdom": "Organic compounds",
        "chemicalTaxonomyClassyfireSuperclass": "Organoheterocyclic compounds",
        "chemicalTaxonomyClassyfireClass": "Imidazopyrimidines",
        "chemicalTaxonomyClassyfireDirectParent": "Purines and purine derivatives",
        "taxonomyReferenceObjects": {
            "NCBI": [{"species": "Coffea arabica", "kingdom": "Plantae"}],
            "GBIF": [{"species": "Coffea arabica", "kingdom": "Plantae"},
                     {"species": "Aspergillus niger", "kingdom": "Fungi"}],
        },
    }]}


def test_identity_annotations_from_all_sources(monkeypatch):
    calls = []

    def fake_http(url, data=None, timeout=15.0):
        calls.append((url, data))
        if "pug_view/categories" in url:
            return as_json({"SourceCategories": {"Categories": [{"Category": "Chemical Vendors", "Sources": [
                {"SourceName": "Sigma-Aldrich"}, {"SourceName": "MolPort"}, {"SourceName": "Some Shop"}]}]}})
        if "pug_view" in url:
            return as_json(pc_payload())
        if "lotus.naturalproducts.net" in url:
            return as_json(lotus_payload())
        if "rest.kegg.jp/link/pathway/C07481" in url:
            return "cpd:C07481\tpath:map00230\ncpd:C07481\tpath:map01100\n"
        if "rest.kegg.jp/list/pathway" in url:
            return "map00230\tPurine metabolism\nmap01100\tMetabolic pathways\n"
        if "reactome.org/AnalysisService" in url:
            assert data == b"27732"
            return as_json({"pathways": [{"stId": "R-HSA-8936459", "name": "Caffeine metabolism"}]})
        if "npclassifier" in url or "classyfire" in url:
            raise AssertionError("complete LOTUS classifications should be reused")
        raise R._NotFound(url)

    monkeypatch.setattr(R, "_http", fake_http)
    result = E.enrich(resolved())
    out = result["outputs"]

    assert set(out) == set(E.OUTPUT_KEYS)
    assert all(block["count"] == len(block["rows"]) for block in out.values())
    assert ["LOTUS / NPClassifier", "Class", "Purine alkaloids"] in out["classification"]["rows"]
    assert ["LOTUS / ClassyFire", "Kingdom", "Organic compounds"] in out["classification"]["rows"]
    assert [row[0] for row in out["species"]["rows"]] == ["Aspergillus niger", "Coffea arabica"]
    origins = {row[0] for row in out["origin"]["rows"]}
    assert origins == {"Drug", "Food", "Environmental pollutant", "Plant", "Microbial", "Human-associated"}
    pathways = {(row[0], row[1]) for row in out["pathways"]["rows"]}
    assert ("Purine metabolism", "HMDB / SMPDB") in pathways
    assert ("Metabolic pathways", "KEGG") in pathways
    assert ("Caffeine metabolism", "Reactome") in pathways
    assert out["standards"]["rows"][0][1:] == ["A · catalog vendor lists this structure", "Catalog 1 · Screening 1 · Other 1"]
    assert len(out["ghs"]["rows"]) == 2
    reporter_a = next(row for row in out["ghs"]["rows"] if row[0]["text"] == "Reporter A")
    assert reporter_a[1:] == ["Danger", "H301; H315"]
    assert out["biofluids"]["rows"] == []
    assert result["warnings"] == []


def test_source_failures_degrade_and_classyfire_job_is_polled(monkeypatch):
    monkeypatch.setattr(E.time, "sleep", lambda _: None)

    def fake_http(url, data=None, timeout=15.0):
        if "npclassifier" in url:
            return as_json({"class_results": ["Purine alkaloids"], "superclass_results": [], "pathway_results": ["Alkaloids"]})
        if "classyfire/classify" in url:
            return as_json({"id": "job-1"})
        if "classyfire/job-1/result" in url:
            return as_json({"entities": [{"classification": {
                "kingdom": {"name": "Organic compounds"}, "superclass": {"name": "Alkaloids"},
                "class": {"name": "Purine alkaloids"}, "direct_parent": {"name": "Purines"},
            }}]})
        raise R.ResolveError("offline")

    monkeypatch.setattr(R, "_http", fake_http)
    result = E.enrich(resolved())
    out = result["outputs"]
    assert set(out) == set(E.OUTPUT_KEYS)
    assert ["NPClassifier", "Class", "Purine alkaloids"] in out["classification"]["rows"]
    assert ["ClassyFire", "Direct parent", "Purines"] in out["classification"]["rows"]
    assert ["Human-associated", out["origin"]["rows"][0][1], out["origin"]["rows"][0][2]] == out["origin"]["rows"][0]
    assert out["species"]["rows"] == [] and out["ghs"]["rows"] == []
    assert all(block["count"] == len(block["rows"]) for block in out.values())
    assert any(w.startswith("LOTUS lookup failed:") for w in result["warnings"])
    assert any(w.startswith("PubChem annotations lookup failed:") for w in result["warnings"])
    assert any(w.startswith("KEGG lookup failed:") for w in result["warnings"])
    assert any(w.startswith("Reactome lookup failed:") for w in result["warnings"])


def test_lipidmaps_classification_is_independent(monkeypatch):
    def fake_http(url, data=None, timeout=15.0):
        assert "lipidmaps.org/rest/compound/lm_id/LMFA01010001/classification" in url
        return as_json({"core": "Fatty Acyls", "main_class": "Fatty acids and conjugates", "sub_class": "Straight chain fatty acids"})

    monkeypatch.setattr(R, "_http", fake_http)
    monkeypatch.setattr(E, "_local_classification", lambda _: ([], False))
    result = E.enrich(resolved(smiles=None, inchikey=None, cid=None,
                               ids={"LIPID MAPS": ["LMFA01010001"]}, links={}))
    assert result["outputs"]["classification"]["rows"] == [
        ["LIPID MAPS", "Category", "Fatty Acyls"],
        ["LIPID MAPS", "Main class", "Fatty acids and conjugates"],
        ["LIPID MAPS", "Subclass", "Straight chain fatty acids"],
    ]


def test_local_hmdb_xml_concentrations(monkeypatch, tmp_path):
    xml = """<?xml version="1.0"?>
    <hmdb xmlns="http://www.hmdb.ca">
      <metabolite><accession>HMDB9999999</accession></metabolite>
      <metabolite>
        <accession>HMDB0001847</accession>
        <normal_concentrations>
          <normal_concentration>
            <biospecimen>Saliva</biospecimen>
            <concentration_value>0.20-1.10</concentration_value>
            <concentration_units>uM</concentration_units>
            <references><reference><pubmed_id>12345</pubmed_id></reference></references>
          </normal_concentration>
          <normal_concentration>
            <biospecimen>Urine</biospecimen>
            <concentration_value>1.0</concentration_value>
            <concentration_units>umol/mmol creatinine</concentration_units>
            <references><reference><reference_text>Example cohort</reference_text></reference></references>
          </normal_concentration>
        </normal_concentrations>
      </metabolite>
    </hmdb>"""
    path = tmp_path / "hmdb.xml"
    path.write_text(xml)
    monkeypatch.setenv("MZDESK_HMDB_XML", str(path))
    monkeypatch.setattr(E, "_local_classification", lambda _: ([], False))
    result = E.enrich(resolved(smiles=None, inchikey=None, cid=None,
                               ids={"HMDB": ["HMDB0001847"]}))
    rows = result["outputs"]["biofluids"]["rows"]
    assert [row[:3] for row in rows] == [["Saliva", "0.20-1.10", "uM"],
                                         ["Urine", "1.0", "umol/mmol creatinine"]]
    assert rows[0][3]["help"] == ["PMID: 12345"]
    assert rows[1][3]["url"] == "https://hmdb.ca/metabolites/HMDB0001847"
    assert result["outputs"]["biofluids"]["count"] == 2


@pytest.mark.parametrize("smiles,expected,absent", [
    ("CC(=O)O", {"Carboxylic acid": 1}, {"Ester", "Ketone"}),
    ("CCOC(C)=O", {"Ester": 1}, {"Carboxylic acid", "Ketone", "Ether"}),
    ("CCO", {"Alcohol": 1}, {"Phenol"}),
    ("Oc1ccccc1", {"Phenol": 1}, {"Alcohol"}),
    ("CC(C)=O", {"Ketone": 1}, {"Aldehyde", "Carboxylic acid"}),
    ("CC(N)=O", {"Amide": 1}, {"Primary amine", "Ketone"}),
    ("NC(N)=O", {"Urea": 1}, {"Primary amine"}),
    ("O=[N+]([O-])c1ccccc1", {"Nitro": 1}, set()),
    ("O=C1OCCCC1", {"Lactone": 1, "Ester": 1}, {"Carboxylic acid"}),
    ("CC(=O)OC(C)=O", {"Carboxylic anhydride": 1}, {"Ester"}),
    ("CC(=O)NC(C)=O", {"Imide": 1}, set()),
    ("CP(=O)(O)O", {"Phosphate / phosphonate": 1}, set()),
])
def test_functional_groups_known_molecules(smiles, expected, absent):
    found = {group.name: group.count for group in FG.find_functional_groups(smiles)}
    for name, count in expected.items():
        assert found.get(name) == count
    for name in absent:
        assert name not in found


def test_functional_group_atoms_are_one_based_and_capped():
    groups = {group.name: group for group in FG.find_functional_groups("CCO")}
    assert groups["Alcohol"].matches == ((2, 1),)
    assert FG.atom_text(groups["Alcohol"].matches) == "3, 2"
    assert FG.atom_text(tuple((i,) for i in range(10)), shown=8).endswith("; …")


def test_all_functional_group_smarts_compile_and_invalid_input_fails():
    from rdkit import Chem

    assert len(FG.RULES) >= 170
    assert all(Chem.MolFromSmarts(smarts) is not None for _, _, smarts, _ in FG.RULES)
    with pytest.raises(ValueError):
        FG.find_functional_groups("not a SMILES")


def test_every_rule_matches_its_own_drawing():
    """The picture shown for a group and the SMARTS that is counted must agree ("*" is the rest of the molecule, read as methyl here)."""
    from rdkit import Chem

    for family, name, smarts, drawing in FG.RULES:
        assert family in FG.FAMILIES
        assert Chem.MolFromSmiles(drawing) is not None, name
        stand_in = Chem.MolFromSmiles(drawing.replace("*", "C"))
        assert stand_in.HasSubstructMatch(Chem.MolFromSmarts(smarts)), f"{name}: {smarts} does not match {drawing}"
        assert "<svg" in FG.pattern_svg(name)


@pytest.mark.parametrize("smiles,expected", [
    ("Cn1cnc2c1c(=O)n(C)c(=O)n2C", {"Xanthine core": 1, "Purine": 1, "Methyl": 3, "Heteroaromatic C=O (pyridone-type)": 2}),
    ("O=c1c(O)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12", {"Phenol": 4, "Flavone core": 1, "Chromone": 1, "Benzene ring": 2}),
    ("CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O", {"Steroid core": 1, "Secondary alcohol": 1, "α,β-Unsaturated carbonyl": 1}),
    ("OCC1OC(Oc2ccccc2)C(O)C(O)C1O", {"O-Glycoside (cyclic acetal)": 1, "Primary alcohol": 1, "Tetrahydropyran-type ring": 1}),
    ("NC(Cc1c[nH]c2ccccc12)C(=O)O", {"α-Amino acid": 1, "Indole": 1, "Carboxylic acid": 1, "Primary amine": 1}),
    ("C[N+]#[C-]", {"Isonitrile": 1}),
])
def test_functional_groups_rings_and_motifs(smiles, expected):
    found = {group.name: group.count for group in FG.find_functional_groups(smiles)}
    for name, count in expected.items():
        assert found.get(name) == count, (name, found)


def test_group_views_pictures_and_site_cap():
    views = FG.group_views("CCCCCCCCCCO")                      # 1 primary alcohol, 1 methyl, ...
    by_name = {g["name"]: g for g in views["groups"]}
    assert by_name["Primary alcohol"]["sites"][0]["atoms"] == "11, 10"
    assert "all_svg" not in by_name["Primary alcohol"]
    many = {g["name"]: g for g in FG.group_views("CC(C)(C)CC(C)(C)CC(C)(C)C")["groups"]}["Methyl"]
    assert many["count"] > FG.MAX_SITE_PICTURES and len(many["sites"]) == FG.MAX_SITE_PICTURES
    assert many["more"] == many["count"] - FG.MAX_SITE_PICTURES and "all_svg" in many
    assert FG.group_views("CCO")["svg"].startswith("<?xml")
