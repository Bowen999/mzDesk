"""Alternative explanations: bounded PubChem candidates and local MetFrag (offline)."""
import json
import shutil
from types import SimpleNamespace

import fake_pubchem as fp
import pytest
from fastapi.testclient import TestClient
from rdkit import Chem

from mzdesk.api.app import app
from mzdesk.exclusivity import pubchem_isomers as check
from mzdesk.ms1 import exact_mass as em
from mzdesk.ms2 import metfrag_engine as mf
from mzdesk.spec import EXAMPLE_STRUCTURES, MODULES

FORMULA = "C3H8O"
TOP = "CC(O)C"
OTHER = "CCCO"
ETHER = "CCOC"
PEAKS = [(30.0, 4.0), (40.0, 100.0), (50.0, 20.0), (61.0653, 10.0)]
TOP_NAME = Chem.MolToSmiles(Chem.MolFromSmiles(TOP))
OTHER_NAME = Chem.MolToSmiles(Chem.MolFromSmiles(OTHER))
ETHER_NAME = Chem.MolToSmiles(Chem.MolFromSmiles(ETHER))
client = TestClient(app)


def pubchem(monkeypatch, rows):
    compounds = []
    structures = {}
    for cid, name, smiles, sources in rows:
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        compounds.append(fp.compound(cid, name, FORMULA, 60.0575149, sources,
                                     inchikey=Chem.MolToInchiKey(mol) if mol else ""))
        structures[cid] = smiles
    fake = fp.FakePubChem(compounds)

    def answer(url, data=None):
        reply = json.loads(fake(url, data))
        if "CanonicalSMILES,IsomericSMILES" in url:
            for row in reply.get("PropertyTable", {}).get("Properties", []):
                if structures[row["CID"]]:
                    # Current PubChem property names, not just the legacy aliases.
                    row["SMILES"] = structures[row["CID"]]
        return json.dumps(reply)

    monkeypatch.setattr(em, "_http", answer)
    em._get_json.cache_clear()
    return fake


def fake_metfrag(monkeypatch, fail=(), no_result=()):
    calls = []

    def run(smiles, peaks, adduct, precursor, depth, tol_da, tol_ppm, **kw):
        calls.append((smiles, list(peaks), adduct, depth, tol_da, tol_ppm, kw))
        if smiles in fail:
            raise mf.MetFragError("engine unavailable", "log")
        explained = {40.0, 50.0} if smiles == TOP_NAME else {40.0}
        return {"scored": smiles not in no_result,
                "annotation": [{"mz": mz, "annotated": mz in explained, "kind": "fragment" if mz in explained else "",
                                "ion_formula": "C2H5", "smiles": "CC"} for mz, _ in peaks]}

    monkeypatch.setattr(mf, "run_metfrag", run)
    return calls


def inputs(smiles=TOP):
    return {"structure": {"mode": "SMILES / InChI", "value": smiles},
            "precursor": {"ion_mode": "Positive", "adduct": "[M+H]+"},
            "databases": ["HMDB"], "ms2": {"source": "Paste", "value": "\n".join(f"{x} {y}" for x, y in PEAKS)}}


def test_candidate_properties_are_bounded_and_missing_is_reported(monkeypatch):
    pubchem(monkeypatch, [(1, "top", TOP, []), (2, "other", OTHER, []), (3, "missing", "", [])])
    res = em.formula_structures(FORMULA, cap=2)
    assert res["retrieved"] == 2 and res["truncated"] and [r["cid"] for r in res["candidates"]] == [1, 2]
    assert res["candidates"][1]["smiles"] == OTHER and res["candidates"][1]["source"] == "PubChem CID 2"
    res = em.formula_structures(FORMULA, cap=3)
    assert not res["truncated"] and res["missing"] == 1


def test_page_uses_pubchem_even_when_only_hmdb_selected_and_reports_scope(monkeypatch):
    pubchem(monkeypatch, [(1, "top", TOP, []), (2, "other", OTHER, []), (3, "ether", ETHER, [])])
    calls = fake_metfrag(monkeypatch)
    response = client.post("/api/modules/exclusivity/run", json=inputs())
    assert response.status_code == 200, response.text
    body = response.json()
    rows = body["outputs"]["isomers"]["rows"]
    assert len(rows) == 3 and rows[1][2] == "PubChem CID 2" and rows[2][2] == "PubChem CID 3"
    assert [r[5]["badge"] for r in rows] == ["your compound", "yes", "yes"]
    assert "can be distinguished from 2 tested isomers." in body["verdict"] and "Not exhaustive" not in body["verdict"]
    assert "3 CIDs retrieved (limit 100; not truncated" in body["verdict"]
    assert "2 other connectivities selected (limit 9, not truncated)" in body["verdict"]
    assert len(calls) == 3 and all(call[2:6] == ("[M+H]+", 2, 0.01, 0.0) for call in calls)
    assert body["outputs"]["same_formula"]["rows"] == []  # neither PubChem isomer is in selected HMDB


def test_candidate_and_comparison_truncation_are_reported(monkeypatch):
    pubchem(monkeypatch, [(i, f"record {i}", OTHER if i % 2 else ETHER, []) for i in range(1, 102)])
    monkeypatch.setattr(check, "MAX_CANDIDATES", 2)
    fake_metfrag(monkeypatch)
    body = client.post("/api/modules/exclusivity/run", json=inputs()).json()
    assert "100 CIDs retrieved (limit 100; truncated)" in body["verdict"]
    assert "1 other connectivities selected (limit 1, truncated)" in body["verdict"]
    assert len(body["outputs"]["isomers"]["rows"]) == 2
    assert "Not exhaustive" not in body["verdict"]


@pytest.mark.parametrize("rows", [[], [(1, "same compound", TOP, [])], [(1, "missing", "", [])]])
def test_no_candidate_never_produces_a_unique_verdict(monkeypatch, rows):
    pubchem(monkeypatch, rows)
    calls = fake_metfrag(monkeypatch)
    body = client.post("/api/modules/exclusivity/run", json=inputs()).json()
    assert "No isomer verdict" in body["verdict"] and "No other usable same-formula connectivity" in body["verdict"]
    assert body["outputs"]["isomers"]["rows"] == [] and calls == []
    if rows and not rows[0][2]:
        assert any("missing properties" in w for w in body["warnings"])


def test_network_failure_is_not_a_unique_answer(monkeypatch):
    def offline(url, data=None):
        raise em.ExactMassError("offline")
    monkeypatch.setattr(em, "_http", offline)
    em._get_json.cache_clear()
    calls = fake_metfrag(monkeypatch)
    body = client.post("/api/modules/exclusivity/run", json=inputs()).json()
    assert "PubChem isomer candidate search failed" in body["verdict"]
    assert "offline" in body["verdict"] and "not ruled out" not in body["verdict"] and calls == []
    assert not body["outputs"]["isomers"]["rows"]


@pytest.mark.parametrize("failed", [OTHER_NAME, ETHER_NAME])
def test_partial_metfrag_failure_is_not_a_unique_answer(monkeypatch, failed):
    pubchem(monkeypatch, [(1, "other", OTHER, []), (2, "ether", ETHER, [])])
    fake_metfrag(monkeypatch, fail=[failed])
    body = client.post("/api/modules/exclusivity/run", json=inputs()).json()
    assert "Not distinguished or not scored" in body["verdict"]
    assert any(w.startswith("Not scored") and "not ruled out" not in w for w in body["warnings"])
    assert "not scored" in [r[5]["badge"] for r in body["outputs"]["isomers"]["rows"]]


def test_total_comparison_timeout_leaves_the_unscored_candidate_in_the_group(monkeypatch):
    pubchem(monkeypatch, [(1, "other", OTHER, []), (2, "ether", ETHER, [])])
    calls = fake_metfrag(monkeypatch)
    ticks = iter([0, 1, 2, 121])
    monkeypatch.setattr(check, "time", SimpleNamespace(monotonic=lambda: next(ticks)))
    body = client.post("/api/modules/exclusivity/run", json=inputs()).json()
    assert len(calls) == 2 and calls[0][-1]["timeout_s"] == 119 and calls[1][-1]["timeout_s"] == 118
    assert "Not distinguished or not scored" in body["verdict"]
    assert any("time limit reached" in w and w.startswith("Not scored") for w in body["warnings"])


@pytest.mark.parametrize("adduct,mz,why", [
    ("[M+H]+", "61.1", "does not match the precursor"),
    ("[2M+H]+", "", "MetFrag cannot use the adduct"),
])
def test_mismatch_or_unsupported_adduct_cannot_give_a_verdict(monkeypatch, adduct, mz, why):
    pubchem(monkeypatch, [(1, "other", OTHER, [])])
    calls = fake_metfrag(monkeypatch)
    form = inputs()
    form["precursor"].update(adduct=adduct, precursor_mz=mz)
    body = client.post("/api/modules/exclusivity/run", json=form).json()
    assert "No isomer verdict" in body["verdict"] and why in body["verdict"]
    assert "PubChem same-formula search" in body["verdict"] and "Not exhaustive" not in body["verdict"] and calls == []


@pytest.mark.parametrize("fail,no_result", [([TOP_NAME], []), ([], [TOP_NAME]), ([], [OTHER_NAME])])
def test_no_engine_result_never_produces_a_unique_verdict(monkeypatch, fail, no_result):
    pubchem(monkeypatch, [(1, "other", OTHER, [])])
    fake_metfrag(monkeypatch, fail=fail, no_result=no_result)
    body = client.post("/api/modules/exclusivity/run", json=inputs()).json()
    assert "No isomer verdict" in body["verdict"] and "no isomer verdict is possible" in body["verdict"]


def test_lipid_skips_isomers_even_with_spectrum(monkeypatch):
    fake = pubchem(monkeypatch, [])
    calls = fake_metfrag(monkeypatch)
    body = client.post("/api/modules/exclusivity/run", json=inputs(EXAMPLE_STRUCTURES["DOPG"])).json()
    assert body["outputs"]["isomers"]["rows"] == [] and calls == []
    assert "lipid" in body["verdict"] and "C=C" in body["verdict"] and "sn-positions" in body["verdict"]
    assert "specific methods or additional evidence" not in body["verdict"]                      # the hedge is in the Docs page, not on the result
    assert not any("CanonicalSMILES,IsomericSMILES" in url for url in fake.urls)


def test_network_confirmation_mentions_pubchem_isomer_search():
    mod = next(m for m in MODULES if m["id"] == "exclusivity")
    assert "isomer candidate search query PubChem" in mod["network"]
    assert "MetFrag runs locally" in mod["network"]
    table = next(b for b in mod["outputs"] if b.get("id") == "isomers")
    assert "Distinguishable by MS2" in table["column_help"] and "Atoms covered" in table["column_help"]


@pytest.mark.skipif(not shutil.which("java") or not mf.JAR.exists(), reason="MetFrag Java runtime or JAR missing")
def test_real_metfrag_scoring_contract():
    result = mf.run_metfrag(TOP, [(43.0542, 100.0), (45.0335, 30.0)], "[M+H]+", None, 2, 0.01, 0.0, timeout_s=30)
    assert result["scored"] is True and len(result["annotation"]) == 2


@pytest.mark.skipif(not shutil.which("java") or not mf.JAR.exists(), reason="MetFrag Java runtime or JAR missing")
def test_real_metfrag_compares_two_pubchem_candidates_without_network():
    search = {"formula": FORMULA, "retrieved": 1, "cap": 100, "truncated": False, "missing": 0, "wrong_formula": 0,
              "candidates": [{"id": "pubchem:1", "name": "propanol", "smiles": OTHER, "source": "PubChem CID 1"}]}
    result = check.discriminate({"id": "compound", "name": "isopropanol", "smiles": TOP}, search,
                                [(43.0542, 100.0), (45.0335, 30.0)], None, "[M+H]+")
    assert result["verdict"] in ("single", "indistinguishable")
    assert len(result["candidates"]) == 2 and all(c["evaluated"] for c in result["candidates"])
    assert result["candidates"][1]["source"] == "PubChem CID 1" and "PubChem same-formula search" in result["source_summary"] and "Not exhaustive" not in result["source_summary"]
