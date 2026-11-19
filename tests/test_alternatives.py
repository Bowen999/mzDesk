"""Alternative explanations: the same-formula and isobaric tables (PubChem answered by tests/fake_pubchem.py) and the page wiring.

The isomer table itself is covered by tests/test_isomers.py; here only its place on the page (three output blocks) and the new
online tables are tested. No real PubChem, no MS-FINDER.
"""
import fake_pubchem as fp
import pytest
from fastapi.testclient import TestClient

from mzdesk.api.app import app
from mzdesk.ms1 import exact_mass as em

client = TestClient(app)
FORMULA = "C8H10N4O2"
CAFFEINE = em._neutral_formula(fp.CAFFEINE["formula"])


@pytest.fixture(autouse=True)
def clear_cache():
    em._get_json.cache_clear()
    yield
    em._get_json.cache_clear()


def install(monkeypatch, compounds):
    fake = fp.FakePubChem(compounds)
    monkeypatch.setattr(em, "_http", fake)
    em._get_json.cache_clear()
    return fake


def run(inputs):
    return client.post("/api/modules/exclusivity/run", json=inputs)


# ---------------------------------------------------------------- same formula
def test_same_formula_collapses_stereoisomers_and_keeps_database_membership(monkeypatch):
    compounds = [
        fp.compound(10, "First", FORMULA, 194.08037557, [fp.HMDB], inchikey="RYYVLZVUVIJVGH-UHFFFAOYSA-N"),
        fp.compound(11, "First stereo", FORMULA, 194.08037557, [fp.HMDB], inchikey="RYYVLZVUVIJVGH-QHFFFAOYSA-N"),
        fp.compound(12, "Other", FORMULA, 194.08037557, [], inchikey="SIQPXVQCUCHWDI-UHFFFAOYSA-N"),
    ]
    install(monkeypatch, compounds)
    res = em.formula_groups(FORMULA, ["HMDB", "PubChem"])
    assert res["examined"] == 3 and not res["truncated"]
    assert [(r["key"], r["name"], r["databases"], r["stereoisomers"], len(r["cids"])) for r in res["rows"]] == [
        ("RYYVLZVUVIJVGH", "First", ["HMDB"], 1, 2),
        ("SIQPXVQCUCHWDI", "Other", ["PubChem"], 0, 1),
    ]


def test_same_formula_follows_the_database_filter(monkeypatch):
    install(monkeypatch, [
        fp.compound(10, "In HMDB", FORMULA, 194.08037557, [fp.HMDB], inchikey="RYYVLZVUVIJVGH-UHFFFAOYSA-N"),
        fp.compound(12, "PubChem only", FORMULA, 194.08037557, [], inchikey="SIQPXVQCUCHWDI-UHFFFAOYSA-N"),
    ])
    res = em.formula_groups(FORMULA, ["HMDB"])
    assert [(r["name"], r["databases"]) for r in res["rows"]] == [("In HMDB", ["HMDB"])]


def test_same_formula_cap_and_invalid_formula(monkeypatch):
    install(monkeypatch, [fp.compound(10 + i, f"C{i}", FORMULA, 194.08037557,
                                      inchikey=f"RYYVLZVUVIJV{i:02d}-UHFFFAOYSA-N") for i in range(3)])
    res = em.formula_groups(FORMULA, ["PubChem"], cap=2)
    assert res["examined"] == 2 and res["truncated"] and "only the first 2" in res["warnings"][0]
    with pytest.raises(ValueError, match="Invalid molecular formula"):
        em.formula_groups("not a formula", ["PubChem"])


# ---------------------------------------------------------------- isobaric
def test_isobaric_groups_by_formula_and_hides_the_query_formula(monkeypatch):
    install(monkeypatch, [
        fp.compound(2519, "Caffeine", FORMULA, 194.08037557, [fp.HMDB], inchikey="RYYVLZVUVIJVGH-UHFFFAOYSA-N"),
        fp.compound(9001, "Isobar", "C7H14N4O3", 194.08037557, [fp.HMDB], inchikey="AAAAAAAAAAAAAA-UHFFFAOYSA-N"),
    ])
    mz = 195.0877
    res = em.isobaric(mz, ["[M+H]+"], 5, ["HMDB", "PubChem"], exclude_formula=FORMULA)
    assert res["same_formula"] == 1 and res["total"] == 1
    (row,) = res["rows"]
    assert row["formula"] == "C7H14N4O3" and row["name"] == "Isobar" and row["databases"] == ["HMDB"]
    assert row["mz"] == pytest.approx(195.0877, abs=1e-3) and abs(row["ppm"]) < 5
    assert len(em.isobaric(mz, ["[M+H]+"], 5, ["HMDB", "PubChem"])["rows"]) == 2


# ---------------------------------------------------------------- page
def test_page_accepts_a_formula_without_a_structure(monkeypatch):
    install(monkeypatch, [fp.CAFFEINE, fp.ENPROFYLLINE, fp.ISOCAFFEINE])
    body = run({"structure": {"mode": "Name / ID / Formula", "value": FORMULA},
                "precursor": {"precursor_mz": "195.0877", "ion_mode": "Positive"},
                "databases": ["HMDB", "PubChem"]}).json()
    assert body["status"] == "ok"
    assert len(body["outputs"]["same_formula"]["rows"]) == 3 and body["outputs"]["isomers"]["rows"] == []
    assert "Formula input: the isomer check needs a structure." in body["verdict"]


def test_page_with_only_a_precursor_mz_runs_the_isobaric_table(monkeypatch):
    install(monkeypatch, [fp.compound(9001, "Isobar", "C7H14N4O3", 194.08037557, [fp.HMDB],
                                      inchikey="AAAAAAAAAAAAAA-UHFFFAOYSA-N")])
    body = run({"precursor": {"precursor_mz": "195.0877", "ion_mode": "Positive", "adduct": "[M+H]+"},
                "databases": ["HMDB"]}).json()
    assert body["status"] == "ok" and body["outputs"]["same_formula"]["rows"] == []
    assert len(body["outputs"]["isobaric"]["rows"]) == 1 and body["outputs"]["isomers"]["rows"] == []
    assert "No isomer verdict" in body["verdict"]


def test_page_needs_a_compound_or_an_mz():
    r = run({"databases": ["HMDB"]})
    assert r.status_code == 400 and "compound or formula" in r.json()["detail"]


def test_page_survives_pubchem_failure(monkeypatch):
    def down(url, data=None):
        raise em.ExactMassError("PubChem could not be reached (timeout or no network). Try again later.")
    monkeypatch.setattr(em, "_http", down)
    em._get_json.cache_clear()
    body = run({"structure": {"mode": "Name / ID / Formula", "value": "C9H10N2O3"},
                "precursor": {"precursor_mz": "195.0877", "ion_mode": "Positive", "adduct": "[M+H]+"},
                "databases": ["HMDB"]}).json()
    assert body["status"] == "ok" and body["outputs"]["same_formula"]["rows"] == []
    assert any("PubChem could not be reached" in w for w in body["warnings"])
    assert any("PubChem could not be reached" in h for h in body["outputs"]["isobaric"]["help"])
