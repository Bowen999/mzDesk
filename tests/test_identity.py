import json
import os
from pathlib import Path

import pytest

from mzdesk.identity import resolve as R

FIX = Path(__file__).parent / "fixtures" / "identity_caffeine.json"


@pytest.fixture(autouse=True)
def tmp_cache(monkeypatch, tmp_path):
    monkeypatch.setenv("MZDESK_CACHE_DIR", str(tmp_path))


@pytest.fixture
def canned(monkeypatch):
    data = json.loads(FIX.read_text())
    calls = []

    def fake(url, data_=None, timeout=15.0):
        key = url + "\n" + (data_.decode() if data_ else "")
        calls.append(key)
        if key not in data:
            raise R._NotFound(url)
        return data[key]

    monkeypatch.setattr(R, "_http", fake)
    return calls


@pytest.mark.parametrize("q,kind", [
    ("caffeine", "Name"),
    ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "SMILES"),
    ("CO", "SMILES"),
    ("water", "Name"),
    ("58-08-2", "CAS"),
    ("58-08-3", "Name"),
    ("2519", "PubChem CID"),
    ("HMDB0001847", "HMDB"),
    ("C07481", "KEGG"),
    ("CHEBI:27732", "ChEBI"),
    ("DB00201", "DrugBank"),
    ("DTXSID7020182", "CompTox"),
    ("dtxsid7020182", "CompTox"),
    ("CHEMBL25", "ChEMBL"),
    ("chembl25", "ChEMBL"),
    ("LMFA01010001", "LIPID MAPS"),
    ("RYYVLZVUVIJVGH-UHFFFAOYSA-N", "InChIKey"),
    ("InChI=1S/CH4O/c1-2/h2H,1H3", "InChI"),
    ("C8H10N4O2", "Name"),
])
def test_detect(q, kind):
    assert R.detect(q) == kind


def test_comptox_and_chembl_are_known_id_kinds():
    assert R.KINDS[-2:] == ["CompTox", "ChEMBL"]
    assert R._norm_query("dtxsid7020182", "CompTox") == "DTXSID7020182"
    assert R._norm_query("chembl25", "ChEMBL") == "CHEMBL25"
    assert R._is_id("DTXSID7020182") and R._is_id("CHEMBL25")
    ids = R._ids_from_synonyms(["DTXSID7020182", "chembl25", "58-08-2"])
    assert ids == {"CompTox": ["DTXSID7020182"], "ChEMBL": ["CHEMBL25"], "CAS": ["58-08-2"]}
    assert R._KIND_SOURCE["CompTox"] == "comptox" and R._KIND_SOURCE["ChEMBL"] == "chembl"


def test_resolve_comptox_reads_the_id_kind(monkeypatch):
    seen = {}

    def find(q, kind, timeout):
        seen.update(q=q, kind=kind)
        return [6623]

    monkeypatch.setattr(R, "_find_cids", find)
    monkeypatch.setattr(R, "_pc_props", lambda cid, timeout: {
        "Title": "Bisphenol A", "SMILES": "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1",
        "InChIKey": "IISBACLAFKSPIT-UHFFFAOYSA-N"})
    r = R.resolve("dtxsid7020182", "CompTox", want_ids=False)
    assert seen == {"q": "DTXSID7020182", "kind": "CompTox"}
    assert r.kind == "CompTox" and r.cid == 6623
    assert r.links["PubChem"].endswith("/6623")


def test_resolve_caffeine_offline(canned):
    r = R.resolve("caffeine")
    assert r.kind == "Name" and r.cid == 2519 and r.name == "Caffeine"
    assert r.inchikey == "RYYVLZVUVIJVGH-UHFFFAOYSA-N"
    assert r.smiles == "Cn1c(=O)c2c(ncn2C)n(C)c1=O"
    assert "58-08-2" in r.ids["CAS"]
    assert r.ids["HMDB"] == ["HMDB0001847"] and r.ids["KEGG"] == ["C07481"]
    assert "CHEBI:27732" in r.ids["ChEBI"] and "DB00201" in r.ids["DrugBank"]
    for label in ("PubChem", "CAS", "HMDB", "KEGG", "ChEBI", "DrugBank", "CompTox", "ChEMBL"):
        assert label in r.links
    assert r.synonyms and len(r.synonyms) <= 50
    assert not any(R._is_id(s) for s in r.synonyms)


def test_structure_inputs_no_network(monkeypatch):
    def boom(*a, **k):
        raise AssertionError("network used")
    monkeypatch.setattr(R, "_http", boom)
    assert R.to_smiles("OCC") == "CCO"
    assert R.to_smiles("InChI=1S/CH4O/c1-2/h2H,1H3") == "CO"
    r = R.resolve("CCO", want_ids=False)
    assert r.inchikey == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"


def test_bad_structure_raises():
    with pytest.raises(R.ResolveError):
        R.resolve("C(C", "SMILES", want_ids=False)


def test_not_found(canned):
    with pytest.raises(R.ResolveError, match="not found"):
        R.resolve("zzzz-not-a-compound")


def test_cache_roundtrip(tmp_path):
    R._cache_put("k", "v")
    assert R._cache_get("k") == "v"


def test_cache_unwritable(monkeypatch):
    monkeypatch.setenv("MZDESK_CACHE_DIR", "/proc/nonexistent/x")
    R._cache_put("k", "v")
    assert R._cache_get("k") is None


@pytest.mark.network
@pytest.mark.skipif(os.environ.get("MZDESK_NETWORK_TESTS") != "1", reason="set MZDESK_NETWORK_TESTS=1")
def test_real_network_caffeine():
    r = R.resolve("58-08-2")
    assert r.cid == 2519 and "HMDB0001847" in r.ids["HMDB"]
