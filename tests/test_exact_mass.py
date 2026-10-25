"""Exact-mass search (MS1 · ms1_search): mass windows, adduct round trips, [M]+ electron, grouping, caps, errors.

Offline: every PubChem call is answered by tests/fake_pubchem.py (conftest patches exact_mass._http for all tests).
Real PubChem: the tests marked `network` (MZDESK_NETWORK_TESTS=1) check the source names, the ID patterns and the example page.
"""
import os
import re
import socket
import time
import types
import urllib.error

import pytest
from fastapi.testclient import TestClient

import fake_pubchem as fp
from mzdesk.api.app import app
from mzdesk.api.common import EngineError, InputError
from mzdesk.api.runners.ms1 import run_search, run_theory
from mzdesk.ms1 import exact_mass as em
from mzdesk.ms1 import msbuddy_engine as bd
from mzdesk.spec import ADDUCTS_NEG, ADDUCTS_POS, CAFFEINE, MODULES

REAL_HTTP = em._http        # taken at import, before conftest swaps in the fake
ELECTRON = 0.00054857990946
ALL_ADDUCTS = ADDUCTS_POS[1:] + ADDUCTS_NEG[1:]
CAFFEINE_MZ = bd.ion_mz_and_formula("C8H10N4O2", "[M+H]+")[0]
client = TestClient(app)
needs_network = pytest.mark.skipif(os.environ.get("MZDESK_NETWORK_TESTS") != "1", reason="set MZDESK_NETWORK_TESTS=1")


def install(monkeypatch, compounds=None, **kw):
    fake = fp.FakePubChem(compounds, **kw)
    monkeypatch.setattr(em, "_http", fake)
    em._get_json.cache_clear()
    return fake


def inputs(mz, ion_mode="Positive", adduct="Auto", tol="5", dbs=("HMDB", "PubChem")):
    return {"precursor": {"precursor_mz": str(mz), "ion_mode": ion_mode, "adduct": adduct}, "tol_ppm": tol, "databases": list(dbs)}


def mass_at(mz, ppm, adduct="[M+H]+"):
    """Neutral mass whose `adduct` ion sits `ppm` away from the observed `mz` (observed = theoretical * (1 + ppm))."""
    return bd.mass_from_ion_mz(mz / (1 + ppm * 1e-6), adduct)


def cell_text(cell):
    return cell["text"] if isinstance(cell, dict) else cell


def rows_of(out):
    return [(g["title"], g["meta"], r) for g in out["hits"]["groups"] for r in g["rows"]]


# ---------------------------------------------------------------- mass windows and adducts
def test_window_edges_are_plus_and_minus_tol_in_mz():
    for adduct in ALL_ADDUCTS:
        (_, lo, hi), = em.mass_windows(400.0, [adduct], 5)
        assert bd.ppm_error(400.0, bd.ion_mz_from_mass(lo, adduct)) == pytest.approx(5.0, abs=1e-6), adduct
        assert bd.ppm_error(400.0, bd.ion_mz_from_mass(hi, adduct)) == pytest.approx(-5.0, abs=1e-6), adduct


def test_dimer_window_is_half_as_wide():
    (_, lo1, hi1), = em.mass_windows(389.168, ["[M+H]+"], 5)
    (_, lo2, hi2), = em.mass_windows(389.168, ["[2M+H]+"], 5)
    (_, lo3, hi3), = em.mass_windows(389.168, ["[2M-H]-"], 5)
    assert hi1 - lo1 == pytest.approx(2 * 5e-6 * 389.168, rel=1e-3)
    assert hi2 - lo2 == pytest.approx((hi1 - lo1) / 2, rel=1e-9)
    assert hi3 - lo3 == pytest.approx((hi1 - lo1) / 2, rel=1e-9)


def test_adduct_arithmetic_is_the_theoretical_mz_page():
    """ion_mz_from_mass / mass_from_ion_mz are the same arithmetic as ion_mz_and_formula, for every adduct of both modes."""
    mass = bd.ion_mz_and_formula("C8H10N4O2", "[M]+")[0] + ELECTRON            # sum of the atom masses
    for adduct in ALL_ADDUCTS:
        mz = bd.ion_mz_and_formula("C8H10N4O2", adduct)[0]
        assert bd.ion_mz_from_mass(mass, adduct) == pytest.approx(mz, abs=1e-9), adduct
        assert bd.mass_from_ion_mz(mz, adduct) == pytest.approx(mass, abs=1e-9), adduct
    assert bd.ion_mz_from_mass(mass, "[M+Zz]+") is None


@pytest.mark.parametrize("adduct", ALL_ADDUCTS)
def test_round_trip_every_adduct(adduct):
    """The theoretical m/z of caffeine for any adduct finds caffeine again: with only that adduct and in Auto mode."""
    mz = bd.ion_mz_and_formula("C8H10N4O2", adduct)[0]
    pool = (ADDUCTS_POS if adduct in ADDUCTS_POS else ADDUCTS_NEG)[1:]
    for adducts in ([adduct], pool):
        res = em.search(mz, adducts, 5, ["PubChem"])
        caffeine = [h for h in res.hits if h.cid == 2519]
        assert len(caffeine) == 1, adducts
        assert caffeine[0].adduct == adduct and abs(caffeine[0].ppm) < 0.1


def test_round_trip_from_the_theoretical_mz_page():
    out = run_theory({"structure": {"mode": "SMILES / InChI", "value": CAFFEINE},
                      "precursor": {"precursor_mz": "", "ion_mode": "Positive", "adduct": "[M+H]+"}})["outputs"]
    shown = out["theory"]["values"]["Precursor m/z"]["value"]
    assert float(shown) == round(CAFFEINE_MZ, 4)
    exact = run_search(inputs(CAFFEINE_MZ, adduct="[M+H]+"))["outputs"]
    assert rows_of(exact)[0][2][0] == "Caffeine" and re.search(r"Δ [+−]0\.0 ppm", rows_of(exact)[0][1])
    rounded = run_search(inputs(shown, adduct="[M+H]+"))["outputs"]       # the 4-decimal text a user would paste back
    assert rows_of(rounded)[0][2][0] == "Caffeine" and "Δ +0.2 ppm" in rows_of(rounded)[0][1]


# ---------------------------------------------------------------- [M]+ and charged compounds
def test_permanent_cation_mass_has_no_electron_correction():
    """PubChem's MonoisotopicMass of choline (C5H14NO+) is the atom sum; the ion m/z is one electron (5.3 ppm here) lighter."""
    mz = bd.ion_mz_and_formula("C5H14NO", "[M]+")[0]
    assert mz == pytest.approx(104.107539070 - ELECTRON, abs=2e-8)
    (_, lo, hi), = em.mass_windows(mz, ["[M]+"], 5)
    assert lo < 104.107539070 < hi              # PubChem's number is inside the window
    assert mz < lo                              # taking the m/z as the mass would miss it
    res = em.search(mz, ["[M]+"], 5, ["PubChem", "HMDB"])
    choline, = res.hits
    assert (choline.cid, choline.adduct, choline.dbs) == (305, "[M]+", ["HMDB"]) and abs(choline.ppm) < 0.001
    assert choline.ids == {"HMDB": "HMDB0000097"}


def test_charged_compounds_match_only_m_plus():
    as_h = bd.ion_mz_from_mass(104.107539070, "[M+H]+")      # a cation cannot take another proton in these lists
    assert [h.cid for h in em.search(as_h, ["[M+H]+"], 5, ["PubChem"]).hits] == []
    assert [h.cid for h in em.search(as_h, ["[M+H]+", "[M]+"], 5, ["PubChem"]).hits] == []


def test_neutral_compound_can_be_the_radical_cation():
    mz = bd.ion_mz_and_formula("C8H10N4O2", "[M]+")[0]
    res = em.search(mz, ["[M+H]+", "[M]+"], 5, ["PubChem"])
    assert {(h.cid, h.adduct) for h in res.hits} == {(2519, "[M]+"), (1676, "[M]+"), (1326, "[M]+")}


def test_anions_are_left_out(monkeypatch):
    anion = fp.compound(107689, "an anion of the same mass", "C3H5O3-", 90.031694, [fp.HMDB], 50, charge=-1)     # inside the [M-H]- window
    lactic_acid = fp.compound(612, "Lactic acid", "C3H6O3", 90.031694, [fp.HMDB], 200)
    install(monkeypatch, [anion, lactic_acid])
    mz = bd.ion_mz_and_formula("C3H6O3", "[M-H]-")[0]
    assert [h.cid for h in em.search(mz, ["[M-H]-"], 5, ["PubChem"]).hits] == [612]


# ---------------------------------------------------------------- hits: dedupe, filter, IDs
def test_one_row_per_cid_across_databases():
    res = em.search(CAFFEINE_MZ, ["[M+H]+"], 5, ["HMDB", "PubChem", "KEGG", "NORMAN SusDat"])
    cids = [h.cid for h in res.hits]
    assert sorted(cids) == [1326, 1676, 2519]
    caffeine = next(h for h in res.hits if h.cid == 2519)
    assert caffeine.dbs == ["HMDB", "KEGG", "NORMAN SusDat"]
    assert caffeine.ids == {"HMDB": "HMDB0001847", "KEGG": "C07481"}      # the C number, not D00528; PubChem has no NORMAN ID
    assert next(h for h in res.hits if h.cid == 1676).ids == {"HMDB": "HMDB0014962"}     # its D / G numbers are not KEGG here
    assert res.totals == {"HMDB": 2, "PubChem": 3, "KEGG": 1, "NORMAN SusDat": 1}


def test_native_id_patterns():
    regs = ["58-08-2", "hmdb0001847", "HMDB01847", "CHEBI:27732", "chebi:99", "DB00201", "DB-123456", "LMGP01010005", "C07481", "D00528", "G12345"]
    assert em.native_ids(regs, ["HMDB", "ChEBI", "DrugBank", "LIPID MAPS", "KEGG"]) == {
        "HMDB": "HMDB0001847", "ChEBI": "CHEBI:27732", "DrugBank": "DB00201", "LIPID MAPS": "LMGP01010005", "KEGG": "C07481"}
    assert em.native_ids(regs, []) == {}                                    # a database counts only when the compound is known to be in it
    assert em.native_ids(["G56441", "D04006"], ["KEGG"]) == {"KEGG": "D04006"}      # a drug without a compound entry; G is not KEGG here
    assert em.native_ids(regs, ["NORMAN SusDat"]) == {}


def test_hits_outside_the_ppm_window_are_dropped(monkeypatch):
    near, far = (fp.compound(900001 + i, f"c{i}", "C9H9", mass_at(CAFFEINE_MZ, ppm)) for i, ppm in enumerate((4.9, -6.0)))
    install(monkeypatch, [near, far], slack=0.002)            # the fake returns both: PubChem's own rounding can do the same
    res = em.search(CAFFEINE_MZ, ["[M+H]+"], 5, ["PubChem"])
    assert [h.cid for h in res.hits] == [900001] and res.hits[0].ppm == pytest.approx(4.9, abs=1e-3)
    assert res.totals["PubChem"] == 2


def test_best_adduct_wins_when_windows_overlap():
    res = em.search(CAFFEINE_MZ, ["[M]+", "[M+H]+"], 6000, ["PubChem"])      # 6000 ppm (~1.2 Da): the [M]+ and [M+H]+ windows overlap
    assert [(h.cid, h.adduct) for h in res.hits if h.cid == 2519] == [(2519, "[M+H]+")]       # one row, under the closer adduct


def test_request_hygiene(monkeypatch):
    fake = install(monkeypatch)
    em.search(CAFFEINE_MZ, ADDUCTS_POS[1:], 5, ["HMDB", "PubChem"])
    esearch = [u for u in fake.urls if "esearch.fcgi" in u]
    assert esearch and all("tool=mzdesk" in u and "db=pccompound" in u for u in esearch)
    assert not any("email" in u.lower() for u in fake.urls)
    assert all(re.search(r"%5BMonoisotopicMass%5D", u) for u in esearch)
    assert "OR" in fake.esearch_terms[0]                      # all adduct windows in one term


def test_unknown_input_is_refused():
    with pytest.raises(ValueError, match="Unknown database"):
        em.search(CAFFEINE_MZ, ["[M+H]+"], 5, ["Nope"])
    with pytest.raises(ValueError, match="Unknown adduct"):
        em.search(CAFFEINE_MZ, ["[M+Zz]+"], 5, ["PubChem"])
    with pytest.raises(ValueError, match="too small"):
        em.search(0.5, ["[M+H]+"], 5, ["PubChem"])


# ---------------------------------------------------------------- the page: groups, order, cells
def test_groups_sorted_by_abs_ppm_with_ties_by_adduct(monkeypatch):
    mz = CAFFEINE_MZ
    install(monkeypatch, [fp.CAFFEINE, fp.ENPROFYLLINE, fp.ISOCAFFEINE,
                          fp.compound(901, "Caffeine monohydrate", "C8H12N4O3", bd.mass_from_ion_mz(mz, "[M+H-H2O]+"), [], 20),   # same ion as [M+H]+
                          fp.compound(902, "minus one", "C9H9", mass_at(mz, -1.0)),
                          fp.compound(903, "plus three", "C10H10", mass_at(mz, 3.0))])
    out = run_search(inputs(mz))["outputs"]
    titles = [g["title"] for g in out["hits"]["groups"]]
    assert titles == ["C8H10N4O2 · [M+H]+", "C8H12N4O3 · [M+H-H2O]+", "C9H9 · [M+H]+", "C10H10 · [M+H]+"]
    metas = [g["meta"] for g in out["hits"]["groups"]]
    assert re.fullmatch(r"m/z 195\.0877 · Δ [+−]0\.0 ppm · 3 compounds", metas[0])
    assert metas[2].endswith("Δ −1.0 ppm · 1 compound") and metas[3].endswith("Δ +3.0 ppm · 1 compound")
    assert out["hits"]["count"] == 6
    names = [r[0] for r in out["hits"]["groups"][0]["rows"]]
    assert names == ["Caffeine", "Enprofylline", "Isocaffeine"]            # database members first, the better-referenced one before the other


def test_cells_ids_and_links():
    out = run_search(inputs(CAFFEINE_MZ))["outputs"]
    rows = out["hits"]["groups"][0]["rows"]
    assert rows[0] == ["Caffeine", {"text": "HMDB0001847", "url": "https://hmdb.ca/metabolites/HMDB0001847"},
                       {"text": "2519", "url": "https://pubchem.ncbi.nlm.nih.gov/compound/2519"}]
    assert rows[2][1] == "—"                                                # Isocaffeine: in no chosen database
    every = run_search(inputs(CAFFEINE_MZ, dbs=["NORMAN SusDat", "DrugBank", "ChEBI", "KEGG", "HMDB", "LIPID MAPS"]))["outputs"]["hits"]["groups"][0]["rows"]
    assert every[0][1] == {"text": "HMDB0001847 · KEGG C07481 · CHEBI:27732 · DB00201 · NORMAN SusDat", "url": "https://hmdb.ca/metabolites/HMDB0001847"}
    norman = run_search(inputs(CAFFEINE_MZ, dbs=["NORMAN SusDat"]))["outputs"]["hits"]["groups"][0]["rows"]
    assert norman == [["Caffeine", {"text": "NORMAN SusDat", "url": "https://pubchem.ncbi.nlm.nih.gov/compound/2519"},
                       {"text": "2519", "url": "https://pubchem.ncbi.nlm.nih.gov/compound/2519"}]]


def test_native_urls():
    cases = {"HMDB": ("HMDB0001847", "https://hmdb.ca/metabolites/HMDB0001847"),
             "DrugBank": ("DB00201", "https://go.drugbank.com/drugs/DB00201"),
             "ChEBI": ("CHEBI:27732", "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:27732"),
             "KEGG": ("C07481", "https://www.kegg.jp/entry/C07481"),
             "LIPID MAPS": ("LMGP01010005", "https://www.lipidmaps.org/databases/lmsd/LMGP01010005")}
    assert {db: em.DB_URL[db].format(i) for db, (i, _) in cases.items()} == {db: url for db, (_, url) in cases.items()}
    assert "NORMAN SusDat" not in em.DB_URL and set(em.SOURCES) == {"PubChem", "NORMAN SusDat"} | set(cases)


def test_sources_cover_the_page_databases():
    from mzdesk.spec import DATABASES
    assert list(em.SOURCES) == DATABASES


def test_summary_cards_and_verdict():
    r = client.post("/api/modules/ms1_search/run", json=next(m for m in MODULES if m["id"] == "ms1_search")["example"])
    assert r.status_code == 200, r.text
    body = r.json()
    values = body["outputs"]["summary"]["values"]
    assert values["Compounds"] == "3" and values["Formulas"] == "1"
    assert values["Closest hit"]["value"] == "Caffeine" and values["Closest hit"]["note"].startswith("C8H10N4O2 [M+H]+ · ")
    assert values["Databases searched"] == {"value": "HMDB, PubChem", "note": "HMDB 2 · PubChem 3"}
    assert body["verdict"].startswith("<b>3 compounds</b> (1 formula) within 5 ppm; closest: <b>Caffeine</b> C8H10N4O2 [M+H]+ (")
    assert body["warnings"] == []


def test_nothing_matched():
    r = client.post("/api/modules/ms1_search/run", json=inputs(50.0))
    body = r.json()
    assert body["outputs"]["hits"] == {"groups": [], "count": 0}
    assert body["outputs"]["summary"]["values"]["Closest hit"] == "—"
    assert body["verdict"] == "Nothing matched within 5 ppm in HMDB and PubChem."


# ---------------------------------------------------------------- the 100 limit
def _crowd(dcnt_by_range):
    out = []
    for (lo, hi), dcnt in dcnt_by_range.items():
        out += [fp.compound(cid, f"c{cid}", "C8H10N4O2", 194.08037557, [], dcnt) for cid in range(lo, hi)]
    return out


def test_cap_takes_the_best_known_first(monkeypatch):
    fake = install(monkeypatch, _crowd({(1000, 1020): 150, (1020, 1080): 50, (1080, 1180): 15, (1180, 1250): 1}))      # 250 compounds
    res = em.search(CAFFEINE_MZ, ["[M+H]+"], 5, ["PubChem"])
    assert res.totals["PubChem"] == 250 and res.shown["PubChem"] == 100 and len(res.hits) == 100
    assert sorted(h.cid for h in res.hits) == list(range(1000, 1100))        # tiers 150, 50, then the 20 lowest CIDs of the DCNT 10-29 tier
    assert [h.tier for h in res.hits][:1] == [0] and [h.cid for h in res.hits][:3] == [1000, 1001, 1002]
    assert len(fake.esearch_terms) <= 6
    assert "250" in res.warnings[0] and "100 best known" in res.warnings[0]


def test_cap_falls_back_to_lowest_cids_without_dcnt(monkeypatch):
    install(monkeypatch, _crowd({(2000, 2300): 0}))
    res = em.search(CAFFEINE_MZ, ["[M+H]+"], 5, ["PubChem"], cap=100)
    assert sorted(h.cid for h in res.hits) == list(range(2000, 2100))


def test_cap_is_per_database_and_small_lists_are_complete(monkeypatch):
    crowd = _crowd({(3000, 3150): 150})
    for c in crowd[:5]:
        c["sources"].add(fp.HMDB)
    install(monkeypatch, crowd)
    res = em.search(CAFFEINE_MZ, ["[M+H]+"], 5, ["HMDB", "PubChem"])
    assert res.totals == {"HMDB": 5, "PubChem": 150} and res.shown == {"HMDB": 5, "PubChem": 100}
    assert len(res.hits) == 100                                            # the five HMDB compounds are among the first 100 PubChem ones


def test_capped_page_says_so(monkeypatch):
    install(monkeypatch, _crowd({(1000, 1300): 150}))
    out = run_search(inputs(CAFFEINE_MZ, adduct="[M+H]+", dbs=["PubChem"]))
    values = out["outputs"]["summary"]["values"]
    assert values["Compounds"] == {"value": "100", "note": "capped at 100 per database"}
    assert values["Databases searched"] == {"value": "PubChem", "note": "PubChem 300 (100 listed)"}
    assert "300" in out["warnings"][0]
    from mzdesk.api import verdicts
    assert "list capped" in verdicts.build("ms1_search", {"tol_ppm": 5}, out["outputs"])


# ---------------------------------------------------------------- slow IDs
def test_slow_ids_do_not_hold_the_page(monkeypatch):
    install(monkeypatch, xref_delay=1.0)
    monkeypatch.setattr(em, "ID_BUDGET", 0.05)
    t0 = time.monotonic()
    out = run_search(inputs(CAFFEINE_MZ, dbs=["HMDB"]))
    assert time.monotonic() - t0 < 0.9                      # did not wait for the 1 s answer
    rows = {r[0]: r for r in out["outputs"]["hits"]["groups"][0]["rows"]}
    assert rows["Caffeine"][1] == {"text": "HMDB", "url": "https://pubchem.ncbi.nlm.nih.gov/compound/2519"}       # no ID yet: the database name
    assert "database IDs of 2 compounds" in out["warnings"][0]
    time.sleep(1.3)             # let the background request finish before the cache is cleared


# ---------------------------------------------------------------- errors
def test_network_failure_is_an_engine_error(monkeypatch):
    def down(url, data=None):
        raise em.ExactMassError("PubChem could not be reached (timeout or no network). Try again later.")
    monkeypatch.setattr(em, "_http", down)
    with pytest.raises(EngineError, match="could not be reached"):
        run_search(inputs(CAFFEINE_MZ))
    r = client.post("/api/modules/ms1_search/run", json=inputs(CAFFEINE_MZ))
    assert r.status_code == 502 and "could not be reached" in r.json()["detail"]


def test_failure_after_the_search_step(monkeypatch):
    fake = fp.FakePubChem()
    def flaky(url, data=None):
        if "/property/" in url:
            raise em.ExactMassError("PubChem is busy (HTTP 429). Try again in a minute.")
        return fake(url, data)
    monkeypatch.setattr(em, "_http", flaky)
    with pytest.raises(EngineError, match="busy"):
        run_search(inputs(CAFFEINE_MZ))


@pytest.fixture
def raw_http(monkeypatch):
    """The real _http with urlopen replaced and no waiting."""
    monkeypatch.setattr(em, "_http", REAL_HTTP)
    monkeypatch.setattr(em, "_wait_turn", lambda kind: None)
    monkeypatch.setattr(em, "time", types.SimpleNamespace(sleep=lambda s: None, monotonic=time.monotonic))
    answers = []

    class Reply:
        def __init__(self, text): self.text = text
        def read(self): return self.text.encode()
        def __enter__(self): return self
        def __exit__(self, *a): return False

    def urlopen(req, timeout):
        answers.append((req.full_url, req.get_method(), timeout))
        a = state.pop(0)
        if isinstance(a, BaseException):
            raise a
        return Reply(a)

    state = []
    monkeypatch.setattr(em.urllib.request, "urlopen", urlopen)
    return state, answers


def _http_error(code):
    return urllib.error.HTTPError("https://x", code, "msg", {}, None)


def test_http_failures_become_short_messages(raw_http):
    state, calls = raw_http
    for failure, message in [(urllib.error.URLError("dns"), "could not be reached"), (socket.timeout(), "could not be reached"),
                             (_http_error(500), "HTTP 500"), (_http_error(403), "HTTP 403")]:
        state[:] = [failure]
        with pytest.raises(em.ExactMassError, match=message):
            em._http("https://pubchem.ncbi.nlm.nih.gov/rest/pug/x")
    assert all(timeout == 20.0 for _, _, timeout in calls)


def test_busy_server_is_retried_once(raw_http):
    state, calls = raw_http
    state[:] = [_http_error(429), '{"ok": 1}']
    assert em._http("https://pubchem.ncbi.nlm.nih.gov/rest/pug/x") == '{"ok": 1}' and len(calls) == 2
    state[:] = [_http_error(503), _http_error(503)]
    with pytest.raises(em.ExactMassError, match="busy"):
        em._http("https://pubchem.ncbi.nlm.nih.gov/rest/pug/x")


def test_pug_404_is_an_empty_result_and_posts_are_posts(raw_http):
    state, calls = raw_http
    state[:] = [_http_error(404)]
    assert em._http(em.PUG + "/property/Title/JSON", b"cid=1") == "{}"
    assert calls[0][1] == "POST"


def test_bad_replies_are_not_cached(raw_http):
    state, calls = raw_http
    state[:] = ["<html>maintenance</html>", '{"a": 1}']
    with pytest.raises(em.ExactMassError, match="unexpected reply"):
        em._get_json("https://pubchem.ncbi.nlm.nih.gov/rest/pug/y")
    assert em._get_json("https://pubchem.ncbi.nlm.nih.gov/rest/pug/y") == {"a": 1}
    assert em._get_json("https://pubchem.ncbi.nlm.nih.gov/rest/pug/y") == {"a": 1} and len(calls) == 2      # the good one is cached


def test_unknown_source_name_is_an_error():
    with pytest.raises(em.ExactMassError, match="does not know the source"):
        em._esearch('(194.0:194.1[MonoisotopicMass]) AND "NORMAN SusDat"[SourceName]', 5)


def test_pubchem_error_reply(monkeypatch):
    monkeypatch.setattr(em, "_http", lambda url, data=None: '{"esearchresult": {"ERROR": "Invalid query"}}')
    with pytest.raises(em.ExactMassError, match="rejected"):
        em._esearch("x", 5)
    monkeypatch.setattr(em, "_http", lambda url, data=None: '{"header": {}}')
    with pytest.raises(em.ExactMassError, match="unexpected reply"):
        em._esearch("y", 5)


def test_requests_are_spaced(monkeypatch):
    clock, slept = [100.0], []
    def sleep(s):
        slept.append(round(s, 3)); clock[0] += s
    monkeypatch.setattr(em, "time", types.SimpleNamespace(monotonic=lambda: clock[0], sleep=sleep))
    monkeypatch.setattr(em, "_next_slot", {"eutils": 0.0, "pug": 0.0})
    for _ in range(4):
        em._wait_turn("eutils")
    em._wait_turn("pug")
    assert slept == [0.35, 0.35, 0.35]                   # <= 3 requests per second for NCBI; the PUG clock is separate
    assert 1 / em.MIN_GAP["eutils"] <= 3 and 1 / em.MIN_GAP["pug"] <= 5


# ---------------------------------------------------------------- input checks
@pytest.mark.parametrize("bad, message", [
    ({"precursor": {"precursor_mz": "", "ion_mode": "Positive", "adduct": "Auto"}}, "required"),
    ({"databases": []}, "at least one database"),
    ({"databases": ["HMDB", "Nope"]}, "Unknown database: Nope"),
    ({"tol_ppm": "0"}, "positive"),
    ({"tol_ppm": "500"}, "at most 100"),
    ({"tol_ppm": "x"}, "number"),
    ({"precursor": {"precursor_mz": "195.0877", "ion_mode": "Positive", "adduct": "[M-H]-"}}, "not a positive-mode adduct"),
    ({"precursor": {"precursor_mz": "20", "ion_mode": "Negative", "adduct": "[M+Cl]-"}}, "too small"),
])
def test_bad_inputs(bad, message):
    with pytest.raises(InputError, match=message):
        run_search({**inputs(195.0877), **bad})


def test_api_returns_400_for_bad_input():
    assert client.post("/api/modules/ms1_search/run", json={**inputs(195.0877), "databases": []}).status_code == 400


def test_negative_mode_runs():
    mz = bd.ion_mz_and_formula("C8H10N4O2", "[M+Cl]-")[0]
    out = run_search(inputs(mz, ion_mode="Negative"))["outputs"]
    assert out["hits"]["groups"][0]["title"] == "C8H10N4O2 · [M+Cl]-"


# ---------------------------------------------------------------- real PubChem
@pytest.mark.parametrize("db, key, native", [("HMDB", "HMDB", "HMDB0001847"), ("DrugBank", "DrugBank", "DB00201"), ("ChEBI", "ChEBI", "CHEBI:27732"),
                                             ("KEGG", "KEGG", "C07481"), ("NORMAN SusDat", None, None)])
@pytest.mark.network
@needs_network
def test_real_caffeine_in_each_database(db, key, native):
    res = em.search(CAFFEINE_MZ, ["[M+H]+"], 1, [db])
    caffeine = next((h for h in res.hits if h.cid == 2519), None)
    assert caffeine is not None, f"caffeine not found in {db}: {res.totals}"
    assert db in caffeine.dbs and abs(caffeine.ppm) < 0.1
    if key:
        assert caffeine.ids[key] == native


@pytest.mark.network
@needs_network
def test_real_lipid_maps_pc_34_1():
    res = em.search(760.5851, ["[M+H]+"], 5, ["LIPID MAPS"])
    popc = next((h for h in res.hits if h.cid == 5497103), None)
    assert popc is not None and popc.formula == "C42H82NO8P" and abs(popc.ppm) < 0.5
    assert popc.ids["LIPID MAPS"] == "LMGP01010005" and "LIPID MAPS" in popc.dbs


@pytest.mark.network
@needs_network
def test_real_choline_as_m_plus_has_no_electron_error():
    mz = bd.ion_mz_and_formula("C5H14NO", "[M]+")[0]
    choline = next(h for h in em.search(mz, ["[M]+"], 5, ["HMDB"]).hits if h.cid == 305)
    assert abs(choline.ppm) < 0.001 and choline.ids["HMDB"].startswith("HMDB")


@pytest.mark.network
@needs_network
def test_real_example_page_and_the_pubchem_cap():
    example = next(m for m in MODULES if m["id"] == "ms1_search")["example"]
    result = run_search(example)
    values = result["outputs"]["summary"]["values"]
    assert values["Closest hit"]["value"] == "Caffeine" and "capped" in values["Compounds"]["note"]
    assert any("PubChem" in w and "100" in w for w in result["warnings"])
    only_pubchem = run_search({**example, "databases": ["PubChem"], "precursor": {**example["precursor"], "adduct": "[M+H]+"}})
    assert 2519 in {int(cell_text(r[2])) for _, _, r in rows_of(only_pubchem["outputs"])}      # caffeine survives the 100 limit
