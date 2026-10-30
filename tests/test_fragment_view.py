"""Fragment view (the "View" pop-up of Alternative explanations and ESI MS/MS Annotation): placing fragments on the structure, the peaks it
considers, the job API and the two pages' wiring. MetFrag is replaced by a stub except in the two tests that need Java."""
import re
import shutil
import time

import fake_pubchem as fp
import pytest
from fastapi.testclient import TestClient
from rdkit import Chem

from mzdesk.api import fragments as api
from mzdesk.api.app import app
from mzdesk.exclusivity import pubchem_isomers as check
from mzdesk.ms1 import exact_mass as em
from mzdesk.ms2 import fragment_view as fv
from mzdesk.ms2 import metfrag_engine as mf

client = TestClient(app)
HAS_JAVA = shutil.which("java") is not None
# a B-type flavan dimer (the kind of structure the feature is made for) and the fragments MetFrag explains for it in negative mode
DIMER = "Oc1ccc(cc1)C1Oc2cc(O)cc(O)c2C(C1O)c1c(O)cc(O)c2CC(O)C(Oc12)c1ccc(O)c(O)c1"
DIMER_PEAKS = [(109.0295, 29), (125.0244, 65), (161.0243, 100), (271.0613, 98), (289.0718, 47), (409.0928, 70), (177.0, 22), (245.0, 55)]


# ---------------------------------------------------------------- placing a fragment on the structure
def test_a_fragment_is_placed_on_the_atoms_with_its_elements():
    parent = Chem.MolFromSmiles(DIMER)
    place = fv.fragment_atoms(parent, "C1=CC(=C(C=C1)O)O", "C6H5O2")               # the catechol ring of the upper unit
    assert place and place["alternatives"] == 0 and len(place["atoms"]) == 8
    symbols = sorted(parent.GetAtomWithIdx(i).GetSymbol() for i in place["atoms"])
    assert symbols == ["C"] * 6 + ["O"] * 2


def test_a_fragment_that_fits_several_places_reports_the_alternatives():
    parent = Chem.MolFromSmiles("OCCOCCO")           # HO-CH2CH2-O-CH2CH2-OH: C2H5O2 fits at either end
    place = fv.fragment_atoms(parent, "OCCO", "C2H5O2")
    assert place and place["alternatives"] >= 1


def test_the_hydrogen_counts_decide_between_places():
    parent = Chem.MolFromSmiles("CC(=O)OC")          # methyl acetate: C-O bonded to C: "CO" fits the methoxy and ... the acetyl C–O too
    place = fv.fragment_atoms(parent, "[CH3]O", "CH3O")
    assert place is not None
    assert parent.GetAtomWithIdx(place["atoms"][0]).GetTotalNumHs() == 3 or parent.GetAtomWithIdx(place["atoms"][-1]).GetTotalNumHs() == 3


@pytest.mark.parametrize("smiles, ion", [
    ("", "C2H5O"),                       # MetFrag gave no structure
    ("CCCCCC", "C6H13"),                 # no such chain in the structure
    ("OCCO", "C3H5O2"),                  # atoms do not add up to the ion formula
])
def test_a_fragment_that_cannot_be_placed_is_none(smiles, ion):
    assert fv.fragment_atoms(Chem.MolFromSmiles("OCCOCCO"), smiles, ion) is None


def test_a_metfrag_smiles_rdkit_cannot_sanitize_still_places():
    parent = Chem.MolFromSmiles(DIMER)
    place = fv.fragment_atoms(parent, "C=CC[OH]C1=CC(=CC(=C1)O)O", "C9H5O3")        # as MetFrag wrote it for m/z 161.0243
    assert place and len(place["atoms"]) == 12


@pytest.mark.parametrize("ion, formula", [
    ("[C5H10O4]+HCOO-", "C5H10O4"),            # the adduct the fragment carries is not part of the fragment
    ("[C3H7O4P+H]-", "C3H7O4P"),               # nor is a hydrogen it moved
    ("[C18H33O2-H]-H-", "C18H33O2"),
    ("[CH2]+HCOO-", "CH2"), ("", ""), ("not an ion", ""),
])
def test_the_fragments_own_formula_is_read_from_metfrags_ion_string(ion, formula):
    assert fv.fragment_formula(ion) == formula


def test_a_fragment_carrying_an_adduct_is_checked_against_its_own_formula():
    parent = Chem.MolFromSmiles("OCC(O)C(O)C(O)C(O)CO")
    place = fv.fragment_atoms(parent, "C(CC(C(CO)O)O)O", fv.fragment_formula("[C5H10O4]+HCOO-"))
    assert place and len(place["atoms"]) == 9


# ---------------------------------------------------------------- the peaks
def test_the_peaks_are_the_ones_the_isomer_check_considers():
    peaks = [(30.0, 4.0), (40.0, 100.0), (50.0, 20.0), (60.0, 0.5), (61.0653, 10.0), (62.07, 2.0), (80.0, 5.0)]    # precursor 61.0653 [M+H]+ of C3H8O
    got = fv.considered_peaks(peaks, 61.0653, "[M+H]+", 0.01, 5)
    assert [mz for mz, _ in got] == [30.0, 40.0, 50.0]               # 0.5 % is too weak; the precursor, its M+1 (62.07) and what lies above it are out
    assert got[1][1] == pytest.approx(100.0) and got[0][1] == pytest.approx(4.0)


def test_the_view_peaks_match_the_isomer_table(monkeypatch):
    peaks = [(30.0, 4.0), (40.0, 100.0), (50.0, 20.0), (60.0, 0.5), (61.0648, 10.0), (62.07, 2.0)]       # C3H8O [M+H]+ = 61.0648

    def run(smiles, pk, *a, **kw):
        return {"scored": True, "annotation": [{"mz": mz, "annotated": False, "kind": ""} for mz, _ in pk]}
    monkeypatch.setattr(mf, "run_metfrag", run)
    search = {"formula": "C3H8O", "retrieved": 2, "cap": 100, "truncated": False, "missing": 0, "wrong_formula": 0,
              "candidates": [{"id": 2, "smiles": "CCCO", "name": "other", "source": "PubChem CID 2"}]}
    res = check.discriminate({"id": "compound", "smiles": "CC(O)C", "name": "top"}, search, peaks, 61.0648, "[M+H]+", tol_da=0.01, tol_ppm=5)
    assert res["peaks_considered"] == 3 == len(fv.considered_peaks(peaks, 61.0648, "[M+H]+", 0.01, 5))


def test_query_is_json_able_and_none_without_fragment_peaks():
    q = fv.query([(40.0, 100.0), (61.0653, 10.0)], 61.0653, "[M+H]+", tol_da=0.01, tol_ppm=5)
    assert q == {"peaks": [[40.0, 100.0]], "precursor_mz": 61.0653, "tol_da": 0.01}
    assert fv.query([(61.0653, 10.0)], 61.0653, "[M+H]+", tol_da=0.01, tol_ppm=5) is None
    assert fv.metfrag_adduct("[M+H]+") == "[M+H]+" and fv.metfrag_adduct("[M-H2O+H]+") is None and fv.metfrag_adduct(None) is None


# ---------------------------------------------------------------- the view itself (Java)
@pytest.mark.skipif(not HAS_JAVA, reason="MetFrag needs Java")
def test_real_view_of_the_flavan_dimer():
    q = fv.query(DIMER_PEAKS + [(561.14, 5)], 561.14, "[M-H]-", tol_da=0.01, tol_ppm=5)
    ticks = []
    view = fv.build(DIMER, q["peaks"], "[M-H]-", q["precursor_mz"], tol_da=0.01, name="dimer", progress=lambda f, t: ticks.append((f, t)))
    assert view["explained"] == 6 and view["peaks"] == 8 and view["formula"] == "C30H26O11" and view["ion_formula"] == "C30H25O11"
    assert [f["ion_formula"] for f in view["fragments"]] == ["C6H5O2", "C6H5O3", "C9H5O3", "C15H11O5", "C15H13O6", "C22H17O8"]
    assert all(f["placed"] and f["svg"].startswith("<?xml") or f["svg"].startswith("<svg") for f in view["fragments"])
    assert all(abs(f["ppm"]) < 1 for f in view["fragments"])                        # the mass errors in the attached figure are all below 1 ppm
    assert len({f["color"] for f in view["fragments"]}) == 6 and view["overview"].lstrip().startswith(("<?xml", "<svg"))
    assert [f for f, _ in ticks] == sorted(f for f, _ in ticks) and ticks[0][0] < 0.1 and ticks[-1] == (1.0, "Done")
    labelled = [c for c in view["figure"]["peaks"] if c.get("label")]
    assert len(labelled) == 6 and all("Δ" in c["label"] for c in labelled)
    assert len(view["figure"]["peaks"]) == 8                                       # one card per peak: hover works on every peak


@pytest.mark.skipif(not HAS_JAVA, reason="MetFrag needs Java")
def test_a_structure_with_nothing_explained_still_gives_a_view():
    view = fv.build("CCO", [(30.0, 100.0), (99.0, 50.0)], "[M+H]+", 47.0491, tol_da=0.01)
    assert view["explained"] == 0 and view["fragments"] == [] and view["figure"]["peaks"]


@pytest.mark.parametrize("value, text", [(-0.027, "0.0"), (0.04, "0.0"), (0.4, "+0.4"), (-0.73, "-0.7"), (12.34, "+12.3")])
def test_the_mass_error_is_written_without_a_negative_zero(value, text):
    assert fv._ppm(value) == text


@pytest.mark.skipif(not HAS_JAVA, reason="MetFrag needs Java")
@pytest.mark.parametrize("smiles, adduct, precursor, peaks", [
    ("C[N+](C)(C)CCO", "[M]+", 104.1070, [(60.0808, 60.0), (58.0651, 40.0), (45.0335, 20.0)]),                  # a permanent cation: MetFrag rejects "C5H14NO+"
    ("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "[M+HCOO]-", 225.0616, [(179.0561, 80.0), (89.0244, 40.0)]),       # fragments that carry the formate
    ("Cn1c(=O)c2c(ncn2C)n(C)c1=O", "[M+Na]+", 217.0696, [(189.0747, 30.0), (161.0798, 25.0)]),
])
def test_charged_and_adduct_structures_are_explained_and_placed(smiles, adduct, precursor, peaks):
    view = fv.build(smiles, peaks, adduct, precursor, tol_da=0.01)
    assert view["explained"] == len(peaks) and view["fragments"] and all(f["placed"] for f in view["fragments"])


def test_a_bad_structure_is_a_value_error():
    with pytest.raises(ValueError):
        fv.build("not a smiles", [(30.0, 100.0)], "[M+H]+", 47.0491)


# ---------------------------------------------------------------- the job API
def request(**over):
    return {"smiles": "CCO", "adduct": "[M+H]+", "peaks": [[30.0, 100.0], [31.0, 20.0]], "precursor_mz": 47.0491, **over}


def poll(job_id, timeout=20):
    deadline = time.time() + timeout
    while time.time() < deadline:
        job = client.get(f"/api/fragments/jobs/{job_id}").json()
        if job["status"] in ("done", "error"):
            return job
        time.sleep(0.05)
    raise AssertionError("the job did not finish")


@pytest.fixture(autouse=True)
def fresh_jobs():
    api._jobs.clear()
    api._cache.clear()
    yield
    api._jobs.clear()
    api._cache.clear()


@pytest.mark.parametrize("over, text", [
    ({"smiles": ""}, "smiles"), ({"smiles": "C" * 4000}, "smiles"), ({"adduct": "[M+X]+"}, "MetFrag cannot use"), ({"peaks": []}, "peaks"),
    ({"peaks": [[30, -1]]}, "above 0"), ({"peaks": [["a", 1]]}, "above 0"), ({"precursor_mz": None}, "precursor_mz is required"),
    ({"tol_da": 0}, "tol_da"), ({"tree_depth": 9}, "tree_depth"), ({"name": 5}, "name"),
])
def test_a_bad_request_is_a_400_with_the_reason(over, text):
    r = client.post("/api/fragments/run", json=request(**over))
    assert r.status_code == 400 and text in r.json()["detail"]


def test_unknown_job_is_404():
    assert client.get("/api/fragments/jobs/nope").status_code == 404


def test_job_reports_progress_and_the_result(monkeypatch):
    def build(smiles, peaks, adduct, precursor_mz, **kw):
        kw["progress"](0.06, "MetFrag is fragmenting the structure")
        kw["progress"](0.6, "Matching peaks to fragments")
        return {"name": kw["name"], "fragments": []}
    monkeypatch.setattr(api.fragment_view, "build", build)
    r = client.post("/api/fragments/run", json=request(name="ethanol"))
    assert r.status_code == 200 and re.fullmatch(r"[0-9a-f]{12}", r.json()["job_id"])
    job = poll(r.json()["job_id"])
    assert job["status"] == "done" and job["percent"] == 100 and job["result"] == {"name": "ethanol", "fragments": []} and job["error"] is None


def test_job_percent_never_goes_back(monkeypatch):
    seen = []

    def build(smiles, peaks, adduct, precursor_mz, **kw):
        for f in (0.5, 0.2, 0.7):
            kw["progress"](f, "x")
            seen.append(client.get(f"/api/fragments/jobs/{next(iter(api._jobs))}").json()["percent"])
        return {}
    monkeypatch.setattr(api.fragment_view, "build", build)
    job_id = client.post("/api/fragments/run", json=request()).json()["job_id"]
    poll(job_id)
    assert seen == sorted(seen) and seen[-1] == 70


@pytest.mark.parametrize("exc, message, log", [
    (ValueError("The structure is not valid."), "The structure is not valid.", ""),
    (mf.MetFragError("MetFrag failed (exit code 1).", "tail of the log"), "MetFrag failed (exit code 1).", "tail of the log"),
    (RuntimeError("boom"), "RuntimeError: boom", ""),
])
def test_errors_come_back_with_their_message(monkeypatch, exc, message, log):
    def build(*a, **kw):
        raise exc
    monkeypatch.setattr(api.fragment_view, "build", build)
    job = poll(client.post("/api/fragments/run", json=request()).json()["job_id"])
    assert job["status"] == "error" and message in job["error"] and job["log"] == log and job["result"] is None


def test_the_same_request_is_answered_from_the_cache(monkeypatch):
    calls = []

    def build(*a, **kw):
        calls.append(1)
        return {"n": len(calls)}
    monkeypatch.setattr(api.fragment_view, "build", build)
    first = poll(client.post("/api/fragments/run", json=request()).json()["job_id"])
    second = client.post("/api/fragments/run", json=request())
    assert second.json()["status"] == "done" and len(calls) == 1
    assert client.get(f"/api/fragments/jobs/{second.json()['job_id']}").json()["result"] == first["result"]
    other = client.post("/api/fragments/run", json=request(tol_da=0.02)).json()      # another tolerance is another request
    assert other["status"] != "done" and poll(other["job_id"])["status"] == "done" and len(calls) == 2


def test_too_many_jobs_are_refused(monkeypatch):
    monkeypatch.setattr(api, "JOB_LIMIT", 1)
    monkeypatch.setattr(api.fragment_view, "build", lambda *a, **kw: {})
    assert client.post("/api/fragments/run", json=request()).status_code == 200
    assert client.post("/api/fragments/run", json=request(tol_da=0.02)).status_code == 503


@pytest.mark.skipif(not HAS_JAVA, reason="MetFrag needs Java")
def test_real_job_end_to_end():
    r = client.post("/api/fragments/run", json={"smiles": DIMER, "adduct": "[M-H]-", "precursor_mz": 561.14, "name": "dimer",
                                                "peaks": [[mz, rel] for mz, rel in DIMER_PEAKS], "tol_da": 0.01})
    job = poll(r.json()["job_id"], timeout=60)
    assert job["status"] == "done", job
    assert job["result"]["explained"] == 6 and len(job["result"]["fragments"]) == 6 and job["seconds"] > 0


# ---------------------------------------------------------------- the two pages
FORMULA = "C3H8O"


def _pubchem(monkeypatch):
    compounds = [fp.compound(1, "top", FORMULA, 60.0575149, [], inchikey=Chem.MolToInchiKey(Chem.MolFromSmiles("CC(O)C"))),
                 fp.compound(2, "other", FORMULA, 60.0575149, [], inchikey=Chem.MolToInchiKey(Chem.MolFromSmiles("CCCO")))]
    fake = fp.FakePubChem(compounds)
    smiles = {1: "CC(O)C", 2: "CCCO"}

    def answer(url, data=None):
        import json
        reply = json.loads(fake(url, data))
        if "CanonicalSMILES,IsomericSMILES" in url:
            for row in reply.get("PropertyTable", {}).get("Properties", []):
                row["SMILES"] = smiles[row["CID"]]
        return json.dumps(reply)
    monkeypatch.setattr(em, "_http", answer)
    em._get_json.cache_clear()


def test_isomers_table_has_a_view_column_and_the_spectrum_for_it(monkeypatch):
    from mzdesk.spec import MODULES
    _pubchem(monkeypatch)
    monkeypatch.setattr(mf, "run_metfrag", lambda smiles, peaks, *a, **kw: {
        "scored": True, "annotation": [{"mz": mz, "annotated": mz == 40.0, "kind": "fragment" if mz == 40.0 else "", "ion_formula": "C2H5", "smiles": "CC"}
                                       for mz, _ in peaks]})
    peaks = [(30.0, 4.0), (40.0, 100.0), (50.0, 20.0), (61.0653, 10.0)]
    body = client.post("/api/modules/exclusivity/run", json={
        "structure": {"mode": "SMILES / InChI", "value": "CC(O)C"}, "precursor": {"ion_mode": "Positive", "adduct": "[M+H]+"}, "databases": ["HMDB"],
        "ms2": {"source": "Paste", "value": "\n".join(f"{x} {y}" for x, y in peaks)}}).json()
    block = next(b for b in next(m for m in MODULES if m["id"] == "exclusivity")["outputs"] if b.get("id") == "isomers")
    assert block["columns"][-1] == "Fragments"
    data = body["outputs"]["isomers"]
    assert all(len(r) == len(block["columns"]) for r in data["rows"])
    cells = [r[-1]["fragment"] for r in data["rows"]]
    assert [c["adduct"] for c in cells] == ["[M+H]+", "[M+H]+"] and {c["smiles"] for c in cells} == {"CC(C)O", "CCCO"}
    query = data["fragments"]                                           # no precursor given: the theoretical m/z of C3H8O [M+H]+ is used
    assert query["peaks"] == [[30.0, 4.0], [40.0, 100.0], [50.0, 20.0]] and query["precursor_mz"] == pytest.approx(61.0648, abs=1e-4) and query["tol_da"] == 0.01


def test_no_view_buttons_without_an_isomer_check(monkeypatch):
    body = client.post("/api/modules/exclusivity/run", json={"precursor": {"ion_mode": "Positive", "adduct": "[M+H]+", "precursor_mz": "100.1"},
                                                              "databases": ["HMDB"]}).json()
    assert "fragments" not in body["outputs"]["isomers"] and body["outputs"]["isomers"]["rows"] == []
