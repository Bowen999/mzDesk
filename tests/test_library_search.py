"""Library search (Spectrum · Library search) tests.

Everything runs offline on two tiny fixture libraries (tests/fixtures/library_massbank.msp, library_gnps.mgf) that are indexed into a
temporary MZDESK_LIBRARY_DIR; the real MassBank / GNPS indexes are never touched. The real-index example runs in tests/test_api.py (skipped
when the index is not installed).
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from fastapi.testclient import TestClient

from mzdesk import spec
from mzdesk.api import verdicts
from mzdesk.api.app import app
from mzdesk.ms2 import library_index as lx
from mzdesk.ms2 import library_parse as lp
from mzdesk.ms2 import library_search as ls
from mzdesk.ms2 import library_setup as setup_mod
from mzdesk.ms2 import metric_names as mn
from mzdesk.ms2 import ms2_similarity as sim
from mzdesk.shared.spectrum import parse_spectrum

FIX = Path(__file__).parent / "fixtures"
MASSBANK, GNPS = FIX / "library_massbank.msp", FIX / "library_gnps.mgf"
GNPS_COPIES = FIX / "library_gnps_copies.mgf"
CAFFEINE = spec.CAFFEINE_MASSBANK_MS2                  # = MSBNK-FIX-0001
EI_QUERY = "55 25\n67 30\n82 150\n109 200\n137 100\n165 120\n194 999"        # = MSBNK-FIX-E001
ENT, COS = mn.ENTROPY_WEIGHTED, mn.COSINE_GREEDY
client = TestClient(app)


def _paste(text):
    return {"source": "Paste", "value": text}


def _prec(mz="195.0877", mode="Positive", adduct="Auto"):
    return {"precursor_mz": mz, "ion_mode": mode, "adduct": adduct}


def _run(**inputs):
    return client.post("/api/modules/ms2_library/run", json=inputs)


def _ms(**kw):
    body = {"spectrum_type": "MS/MS", "ms2": _paste(CAFFEINE), "precursor": _prec()}
    body.update(kw)
    return body


def _build(names=("MassBank", "GNPS"), gnps=GNPS):
    for name in names:
        setup_mod.build_from(name, MASSBANK if name == "MassBank" else gnps,
                             {"url": "file://fixture", "version": "fixture-1", "date": "2026-01-01"}, progress=lambda *_: None)


@pytest.fixture(scope="module")
def libs(tmp_path_factory):
    """Both fixture libraries indexed into a temporary library directory (module-wide)."""
    mp = pytest.MonkeyPatch()
    mp.setenv("MZDESK_LIBRARY_DIR", str(tmp_path_factory.mktemp("libs")))
    lx.clear_cache()
    _build()
    yield
    mp.undo()
    lx.clear_cache()


@pytest.fixture
def copies_libs(monkeypatch, tmp_path):
    """MassBank fixture + a GNPS fixture that holds copies of MassBank records."""
    monkeypatch.setenv("MZDESK_LIBRARY_DIR", str(tmp_path / "with_copies"))
    lx.clear_cache()
    _build(gnps=GNPS_COPIES)
    yield
    lx.clear_cache()


@pytest.fixture
def massbank_only(monkeypatch, tmp_path):
    monkeypatch.setenv("MZDESK_LIBRARY_DIR", str(tmp_path / "only_massbank"))
    lx.clear_cache()
    _build(("MassBank",))
    yield
    lx.clear_cache()


@pytest.fixture
def no_libs(monkeypatch, tmp_path):
    monkeypatch.setenv("MZDESK_LIBRARY_DIR", str(tmp_path / "nothing"))
    lx.clear_cache()
    yield
    lx.clear_cache()


def _accessions(lib, ids):
    return {r["accession"] for r in lib.info(ids)}


def _search(peaks_text, kind="MS/MS", libraries=("MassBank",), metrics=(ENT, COS), threshold=0.7, **kw):
    kw.setdefault("precursor_mz", 195.0877 if kind == "MS/MS" else None)
    return ls.search(parse_spectrum(peaks_text), kind, list(libraries), list(metrics), threshold, **kw)


# ---------------------------------------------------------------- names / spec
def test_metric_names_are_the_similarity_page_names():
    res = sim.compute_all_similarities([(100, 1), (150, 3), (200, 5)], [(100, 2), (151, 3), (200, 4)], 250, 251, use_deep=False)
    assert [r["method"] for r in res] == mn.CLASSIC and len(mn.CLASSIC) == 14
    assert [r["category"] for r in res] == [mn.CATEGORY[n] for n in mn.CLASSIC]
    assert mn.ALL == mn.CLASSIC + [mn.MS2DEEPSCORE, mn.SPEC2VEC] and len(mn.ALL) == 16
    assert sorted(spec.LIBRARY_METRICS) == sorted(mn.ALL)


def test_deep_names_come_from_the_same_constants(monkeypatch):
    monkeypatch.setitem(sys.modules, "ms2deepscore", None)
    monkeypatch.setitem(sys.modules, "gensim", None)
    res = sim.compute_all_similarities([(100, 1), (150, 3), (200, 5)], [(100, 2), (151, 3), (200, 4)], 250, 251, use_deep=True)
    assert [r["method"] for r in res] == mn.ALL


def test_library_page_definition():
    mod = next(m for m in spec.MODULES if m["id"] == "ms2_library")
    assert mod["live"] and mod["example"]
    assert [f["id"] for f in mod["inputs"] if f.get("id")] == ["spectrum_type", "ms2", "precursor", "tol_ppm", "tol_da", "libraries", "metrics", "threshold"]
    f = {x["id"]: x for x in mod["inputs"] if x.get("id")}
    assert f["spectrum_type"]["options"] == ["MS/MS", "EI"] and f["spectrum_type"]["value"] == "MS/MS"
    assert f["libraries"]["options"] == ["MassBank", "GNPS"] and f["libraries"]["value"] == ["MassBank", "GNPS"]
    assert f["metrics"]["options"] == spec.LIBRARY_METRICS and f["metrics"]["value"] == [ENT, COS]
    assert f["tol_ppm"]["value"] == 10.0 and f["threshold"]["value"] == 0.7
    assert "MS/MS only" in f["precursor"]["label"] and "MS/MS only" in f["tol_da"]["label"]
    assert f["precursor"]["required_when"] == {"field": "spectrum_type", "equals": "MS/MS"}
    assert f["precursor"]["visible_when"] == {"field": "spectrum_type", "equals": "MS/MS"}
    assert f["tol_da"]["visible_when"] == {"field": "spectrum_type", "equals": "MS/MS"}
    assert f["precursor"]["auto_label"] == "Auto (any adduct)"
    assert "MoNA" not in spec.LIBRARIES and "MoNA" not in mod["one_liner"]
    assert set(setup_mod.NAMES) == set(spec.LIBRARIES)
    hits = next(b for b in mod["outputs"] if b.get("id") == "hits")
    assert hits["columns"] == ["Rank", "Library", "Accession", "Name", "Adduct", "Precursor m/z", "Δ ppm", "Matched peaks", ENT, COS, "Instrument / CE"]
    text = " ".join(hits["help"])
    assert "first selected metric" in text and "threshold" in text
    assert "first selected metric" in f["metrics"]["help"] and "header" in f["precursor"]["help"]       # the ranking rule and the header fallback stay on the form
    assert [b["id"] for b in mod["outputs"]] == ["summary", "hits", "mirror"]
    assert [lbl for lbl, *_ in [(i["label"],) for i in mod["outputs"][0]["items"]]] == [
        "Best hit", "Rank score", "Hits ≥ threshold", "Candidates scored", "Libraries"]


def test_spec_stays_light():
    code = ("import sys, mzdesk.spec\n"
            "bad = {'matchms', 'ms_entropy', 'scipy', 'gensim', 'torch', 'ms2deepscore'} & set(sys.modules)\n"
            "assert not bad, bad")
    done = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert done.returncode == 0, done.stderr


# ---------------------------------------------------------------- adducts
@pytest.mark.parametrize("raw,expected", [
    ("[M+H]+", "[M+H]+"), ("M+H", "[M+H]+"), ("[M+H]1+", "[M+H]+"), ("[M-H]-", "[M-H]-"), ("M-H", "[M-H]-"), ("M+Na", "[M+Na]+"),
    ("[M+2H]++", "[M+2H]2+"), ("[M+2H]2+", "[M+2H]2+"), ("[M-2H]--", "[M-2H]2-"),
    ("[M-H2O+H]+", "[M+H-H2O]+"), ("[M+H-H2O]+", "[M+H-H2O]+"), ("[M-H2O-H]-", "[M-H-H2O]-"),
    ("[M+FA-H]-", "[M+HCOO]-"), ("[M+HCOO]-", "[M+HCOO]-"), ("[M+CH3COOH-H]-", "[M+CH3COO]-"), ("[M+Hac-H]-", "[M+CH3COO]-"),
    ("[M]+*", "[M]+"), ("[M]+", "[M]+"), ("[2M+H]+", "[2M+H]+"), ("2M+H", "[2M+H]+"), ("[M+NH4]+", "[M+NH4]+"),
    ("[M+CH3COO]-/[M-CH3]-", None), ("N/A", None), ("", None), (None, None), ("something", None),
])
def test_adduct_notations_are_normalised(raw, expected):
    assert lp.normalize_adduct(raw) == expected


def test_every_ui_adduct_has_a_canonical_form():
    for name in spec.ADDUCTS_POS[1:] + spec.ADDUCTS_NEG[1:]:
        canon = lp.normalize_adduct(name)
        assert canon is not None and lp.adduct_charge_sign(canon) == (1 if name.endswith("+") else -1), name


# ---------------------------------------------------------------- parsers
def test_massbank_parser_keeps_ms2_and_ei_and_says_why_it_skips():
    kept, skipped = {}, {}
    for rec, why in lp.iter_records("MassBank", MASSBANK):
        if rec is None:
            skipped[why] = skipped.get(why, 0) + 1
        else:
            kept[rec["accession"]] = rec
    assert len(kept) == 15 and sum(skipped.values()) == 4
    assert skipped["MS/MS without a precursor m/z"] == 1 and sum(v for k, v in skipped.items() if k.startswith("not MS/MS")) == 3
    r = kept["MSBNK-FIX-0001"]
    assert (r["kind"], r["ion_mode"], r["adduct"], r["precursor_mz"]) == (lp.MSMS, 1, "[M+H]+", 195.0877)
    assert r["inchikey"] == "RYYVLZVUVIJVGH-UHFFFAOYSA-N" and r["collision_energy"] == "45" and "Q Exactive" in r["instrument"] and r["instrument"].startswith("LC-ESI-QFT")
    assert len(r["mz"]) == 19 and r["intensity"].dtype == np.float32 and r["intensity"].max() == 1.0 and np.all(np.diff(r["mz"]) >= 0)
    assert kept["MSBNK-FIX-0006"]["adduct"] == "[M+H]+"                 # "M+H" written without brackets
    assert kept["MSBNK-FIX-0008"]["adduct"] is None
    assert kept["MSBNK-FIX-0007"]["ion_mode"] == -1
    ei = {a: r for a, r in kept.items() if r["kind"] == lp.EI}
    assert sorted(ei) == [f"MSBNK-FIX-E00{i}" for i in range(1, 6)]
    assert all(r["precursor_mz"] is None for r in ei.values())              # E002 has a PrecursorMZ line; EI ignores it
    assert "MSBNK-FIX-E006" not in kept                                      # CI-B is not EI


def test_gnps_parser_handles_ion_mode_adduct_in_name_and_inchikey_field():
    kept, skipped = {}, {}
    for rec, why in lp.iter_records("GNPS", GNPS):
        if rec is None:
            skipped[why] = skipped.get(why, 0) + 1
        else:
            kept[rec["accession"]] = rec
    assert len(kept) == 8 and skipped == {"MS/MS without a precursor m/z": 1, "no valid peaks": 1, "not MS/MS or EI": 1}
    r = kept["CCMSLIB00000FIX001"]
    assert (r["adduct"], r["ion_mode"], r["precursor_mz"], r["kind"]) == ("[M+H]+", 1, 195.0877, lp.MSMS)
    assert r["inchikey"] == "RYYVLZVUVIJVGH-UHFFFAOYSA-N"                    # GNPS puts the key in the INCHI field
    assert kept["CCMSLIB00000FIX002"]["adduct"] == "[M+H]+" and kept["CCMSLIB00000FIX002"]["ion_mode"] == 1      # "[M+H]+" at the end of the name, "positive"
    assert kept["CCMSLIB00000FIX003"]["adduct"] == "[M+Na]+" and kept["CCMSLIB00000FIX003"]["ion_mode"] == 1     # " Positive" with a space
    assert kept["CCMSLIB00000FIX004"]["ion_mode"] == 0 and kept["CCMSLIB00000FIX004"]["adduct"] is None
    assert kept["CCMSLIB00000FIX006"]["ion_mode"] == -1
    assert kept["CCMSLIB00000FIX011"]["kind"] == lp.EI and kept["CCMSLIB00000FIX011"]["precursor_mz"] is None


# ---------------------------------------------------------------- index
def test_index_layout_counts_and_metadata(libs):
    lib = lx.load("MassBank")
    assert lx.load("MassBank") is lib                                         # cached in the process
    assert (lib.n, lib.n_msms) == (15, 10) and lib.version == "fixture-1" and lib.label() == "MassBank fixture-1"
    assert lib.meta["counts"] == {"msms_positive": 8, "msms_negative": 2, "msms_unknown_mode": 0, "ei": 5, "msms_without_adduct": 1, "massbank_copies": 0}
    assert int(np.asarray(lib.copy_flag).sum()) == 0
    assert lib.meta["source_url"] == "file://fixture" and lib.meta["date"] == "2026-01-01" and lib.meta["format"] == lx.FORMAT
    assert sum(lib.meta["skipped"].values()) == 4
    for f in lx.ARRAYS + ("meta",):
        assert (lx.index_dir("MassBank") / (f + (".sqlite" if f == "meta" else ".npy"))).is_file()
    assert (lx.index_dir("MassBank") / "meta.json").is_file()
    assert isinstance(lib.peak_mz, np.memmap) and isinstance(lib.peak_int, np.memmap) and isinstance(lib.prec_mz, np.memmap)
    assert np.all(np.diff(lib.prec_mz[:lib.n_msms]) >= 0) and np.all(np.asarray(lib.kind[:lib.n_msms]) == lp.MSMS) and np.all(np.asarray(lib.kind[lib.n_msms:]) == lp.EI)
    row = int(np.flatnonzero([r["accession"] == "MSBNK-FIX-0001" for r in lib.info(range(lib.n))])[0])
    mz, rel = lib.peaks(row)
    assert len(mz) == 19 and rel.max() == pytest.approx(1.0) and mz[0] == pytest.approx(56.0496)
    info = lib.info([row])[0]
    assert info["name"] == "Caffeine" and info["adduct"] == "[M+H]+" and info["collision_energy"] == "45" and info["smiles"].startswith("CN1")
    assert lx.load("GNPS").meta["counts"]["msms_unknown_mode"] == 1 and lx.load("GNPS").meta["counts"]["ei"] == 1


def test_gather_gives_the_same_peaks_as_single_reads(libs):
    lib = lx.load("MassBank")
    ids = np.array([3, 1, 7])
    mz, rel, lens = lib.gather(ids)
    off = np.concatenate([[0], np.cumsum(lens)])
    for k, i in enumerate(ids):
        a, b = lib.peaks(i)
        assert np.array_equal(mz[off[k]:off[k + 1]], a) and np.array_equal(rel[off[k]:off[k + 1]], b)
    mz2, _, lens2 = lib.gather(lib.ei_ids())                                  # contiguous rows take the slice path
    assert lens2.sum() == len(mz2) == int(lib.offsets[lib.n] - lib.offsets[lib.n_msms])


def test_index_in_a_folder_with_spaces_and_a_hash(tmp_path, monkeypatch):
    monkeypatch.setenv("MZDESK_LIBRARY_DIR", str(tmp_path / "my libs #1"))
    lx.clear_cache()
    _build(("MassBank",))
    assert lx.load("MassBank").info([0])[0]["accession"].startswith("MSBNK-FIX")
    lx.clear_cache()


# ---------------------------------------------------------------- MS/MS candidates
def test_precursor_window_ion_mode_and_adduct_filters(libs):
    lib = lx.load("MassBank")
    ids, ppm = lib.msms_candidates(195.0877, 10, True)
    assert _accessions(lib, ids) == {f"MSBNK-FIX-000{i}" for i in (1, 2, 3, 5, 6, 8)}      # not 0004 (12 ppm), 0007 (negative), 0009 / 0010 (other masses)
    assert np.all(np.abs(ppm) <= 10)
    info = {r["accession"]: p for r, p in zip(lib.info(ids), ppm)}
    assert info["MSBNK-FIX-0003"] == pytest.approx((195.0877 - 195.0880) / 195.0880 * 1e6)       # query minus library
    assert "MSBNK-FIX-0004" in _accessions(lib, lib.msms_candidates(195.0877, 100, True)[0])
    assert _accessions(lib, lib.msms_candidates(195.0877, 10, False)[0]) == {"MSBNK-FIX-0007"}   # negative mode
    # a specific adduct keeps only that precursor type (and the "M+H" notation counts as [M+H]+); unknown adduct = none
    assert _accessions(lib, lib.msms_candidates(195.0877, 10, True, "[M+H]+")[0]) == {f"MSBNK-FIX-000{i}" for i in (1, 2, 3, 6)}
    assert _accessions(lib, lib.msms_candidates(195.0877, 10, True, "[M+NH4]+")[0]) == {"MSBNK-FIX-0005"}
    assert lib.msms_candidates(195.0877, 10, True, "[M+Na]+")[0].size == 0
    assert lx.load("GNPS").msms_candidates(195.0877, 10, True)[0].size == 3                   # FIX001-003; FIX004 has no ion mode, FIX007 is negative
    assert lx.load("GNPS").msms_candidates(195.0877, 10, False)[0].size == 1


def test_search_filters_adduct_by_normalised_notation(libs):
    res = _search(CAFFEINE, metrics=[ENT], threshold=-1, adduct="M+H")
    assert {h["accession"] for h in res["hits"]} == {f"MSBNK-FIX-000{i}" for i in (1, 2, 3, 6)} and res["n_candidates"] == 4
    res = _search(CAFFEINE, metrics=[ENT], threshold=-1)
    assert res["n_candidates"] == 6
    with pytest.raises(ValueError, match="does not match"):
        _search(CAFFEINE, metrics=[ENT], adduct="[M-H]-")
    with pytest.raises(ValueError, match="Unknown adduct"):
        _search(CAFFEINE, metrics=[ENT], adduct="banana")


def test_too_many_candidates_keeps_the_closest_and_warns(libs, monkeypatch):
    monkeypatch.setattr(ls, "MAX_MSMS_CANDIDATES", 3)
    res = _search(CAFFEINE, metrics=[ENT], threshold=-1)
    assert res["n_candidates"] == 3
    assert {h["accession"] for h in res["hits"]} == {"MSBNK-FIX-0001", "MSBNK-FIX-0002", "MSBNK-FIX-0008"}      # exactly 0.0 ppm away
    assert any("6 library spectra" in w and "3 closest" in w for w in res["warnings"])


def test_no_candidate_is_a_warning_not_an_error(libs):
    res = _search(CAFFEINE, metrics=[ENT], precursor_mz=500.0)
    assert res["n_candidates"] == 0 and res["hits"] == [] and res["n_kept"] == 0
    assert any("No library spectrum within 10 ppm of m/z 500.0000" in w for w in res["warnings"])


# ---------------------------------------------------------------- scoring
def test_precursor_is_required_for_ms_ms(libs):
    with pytest.raises(ValueError, match="Precursor m/z is required"):
        _search(CAFFEINE, precursor_mz=None)
    r = _run(**_ms(precursor=_prec("")))
    assert r.status_code == 400 and "Precursor m/z is required" in r.json()["detail"]
    # a precursor in an MGF / MSP header counts (like the other pages)
    header = f"BEGIN IONS\nPEPMASS=195.0877\n{CAFFEINE}\nEND IONS"
    r = _run(**_ms(ms2=_paste(header), precursor=_prec("")))
    assert r.status_code == 200 and any("195.0877 is read from the spectrum" in w for w in r.json()["warnings"])


def test_scores_equal_the_similarity_page_for_every_candidate(libs):
    lib = lx.load("MassBank")
    q = parse_spectrum(CAFFEINE)
    res = _search(CAFFEINE, metrics=mn.CLASSIC, threshold=-1)
    assert res["metrics"] == [m for m in mn.LIBRARY_ORDER if m in mn.CLASSIC] and len(res["hits"]) == 6
    for h in res["hits"]:
        mz, rel = lib.peaks(h["row"])
        ref = {r["method"]: r["score"] for r in sim.compute_all_similarities(q, list(zip(mz, rel.astype(float) * 100)), 195.0877, h["precursor_mz"], use_deep=False)}
        for m in mn.CLASSIC:
            got, want = h["scores"][m], ref[m]
            assert (got is None) == (want is None), (h["accession"], m)
            if got is not None:
                assert got == pytest.approx(want, abs=1e-9), (h["accession"], m)
    best = res["hits"][0]
    assert best["accession"] == "MSBNK-FIX-0001" and best["scores"][ENT] == pytest.approx(1.0)


def test_matched_peaks_are_the_greedy_pairs(libs):
    from mzdesk.ms2 import mirror_plot
    res = _search(CAFFEINE, metrics=[ENT], threshold=-1)
    lib = lx.load("MassBank")
    q = parse_spectrum(CAFFEINE)
    for h in res["hits"]:
        mz, rel = lib.peaks(h["row"])
        assert h["matches"] == mirror_plot.match_summary(q, list(zip(mz, rel)), 0.02)["matched_pairs"], h["accession"]


def test_binned_fast_path_equals_the_page_functions():
    rng = np.random.default_rng(7)
    for tol in (0.02, 0.5):
        for _ in range(25):
            mz1 = np.sort(rng.uniform(50, 400, rng.integers(3, 40)).round(3))
            mz2 = np.sort(np.concatenate([mz1[: rng.integers(1, 3)], rng.uniform(50, 400, rng.integers(2, 40)).round(3)]))
            i1, i2 = rng.uniform(1, 100, mz1.size), rng.uniform(1, 100, mz2.size)
            got = ls._binned_scores(mz1, i1, mz2, i2, tol, {mn.BINNED_COSINE, mn.BINNED_PEARSON, mn.BINNED_SPEARMAN})
            assert got[mn.BINNED_COSINE] == pytest.approx(sim.binned_cosine(mz1, i1, mz2, i2, tol), abs=1e-12)
            for name, method in ((mn.BINNED_PEARSON, "pearson"), (mn.BINNED_SPEARMAN, "spearman")):
                want = sim.binned_correlation(mz1, i1, mz2, i2, tol, method)
                assert (got[name] is None) == (want is None)
                if want is not None:
                    assert got[name] == pytest.approx(want, abs=1e-10)


def _stub_scores(monkeypatch, table, matches):
    """Replace the per-pair scorer: scores / matches come from tables keyed by the library accession."""
    lib = lx.load("MassBank")
    by_peaks = {}
    for row in range(lib.n_msms):
        mz, _ = lib.peaks(row)
        by_peaks[tuple(np.round(mz, 4))] = lib.info([row])[0]["accession"]

    def fake(self, mz2, i2, p2):
        acc = by_peaks[tuple(np.round(mz2, 4))]
        return {m: table[acc].get(m) for m in self.pair_metrics}, matches.get(acc, 0)

    monkeypatch.setattr(ls._Scorer, "score", fake)


def test_any_metric_at_or_above_the_threshold_keeps_a_hit_and_na_is_ignored(libs, monkeypatch):
    table = {"MSBNK-FIX-0001": {ENT: 0.60, COS: 0.95}, "MSBNK-FIX-0002": {ENT: 0.90, COS: 0.50}, "MSBNK-FIX-0003": {ENT: 0.90, COS: 0.40},
             "MSBNK-FIX-0006": {ENT: None, COS: 0.99}, "MSBNK-FIX-0005": {ENT: 0.65, COS: None}, "MSBNK-FIX-0008": {ENT: 0.20, COS: 0.20}}
    matches = {"MSBNK-FIX-0001": 5, "MSBNK-FIX-0002": 3, "MSBNK-FIX-0003": 7, "MSBNK-FIX-0006": 9, "MSBNK-FIX-0005": 9, "MSBNK-FIX-0008": 1}
    _stub_scores(monkeypatch, table, matches)
    res = _search(CAFFEINE, metrics=[ENT, COS], threshold=0.7)
    # kept: 0001 (cosine), 0002 / 0003 (entropy), 0006 (cosine; entropy N/A). Not kept: 0005 (0.65 and N/A), 0008
    assert res["n_kept"] == 4 and res["n_candidates"] == 6
    # ranked by the first metric (entropy); 0002 and 0003 tie on entropy -> more matched peaks first (0003 has 7); the N/A entropy goes last
    assert [h["accession"] for h in res["hits"]] == ["MSBNK-FIX-0003", "MSBNK-FIX-0002", "MSBNK-FIX-0001", "MSBNK-FIX-0006"]
    assert res["hits"][-1]["scores"][ENT] is None and [h["rank"] for h in res["hits"]] == [1, 2, 3, 4]
    # a single metric: only that metric decides
    res = _search(CAFFEINE, metrics=[COS], threshold=0.7)
    assert [h["accession"] for h in res["hits"]] == ["MSBNK-FIX-0006", "MSBNK-FIX-0001"]
    # the threshold is inclusive
    assert _search(CAFFEINE, metrics=[COS], threshold=0.95)["n_kept"] == 2
    assert _search(CAFFEINE, metrics=[COS], threshold=0.96)["n_kept"] == 1


def test_ranking_follows_the_option_list_order_not_the_click_order(libs, monkeypatch):
    assert mn.LIBRARY_ORDER.index(ENT) < mn.LIBRARY_ORDER.index(COS)
    assert ls.order_metrics([COS, ENT]) == [ENT, COS] == ls.order_metrics([ENT, COS])
    assert ls.order_metrics([mn.SPEC2VEC, mn.BINNED_COSINE, COS]) == [COS, mn.BINNED_COSINE, mn.SPEC2VEC]
    table = {"MSBNK-FIX-0001": {ENT: 0.80, COS: 0.71}, "MSBNK-FIX-0002": {ENT: 0.72, COS: 0.99}, "MSBNK-FIX-0003": {ENT: 0.75, COS: 0.80}}
    _stub_scores(monkeypatch, {**table, **{a: {ENT: 0.0, COS: 0.0} for a in ("MSBNK-FIX-0005", "MSBNK-FIX-0006", "MSBNK-FIX-0008")}}, {})
    clicked_cos_first = _search(CAFFEINE, metrics=[COS, ENT], threshold=0.7)
    clicked_ent_first = _search(CAFFEINE, metrics=[ENT, COS], threshold=0.7)
    expect = ["MSBNK-FIX-0001", "MSBNK-FIX-0003", "MSBNK-FIX-0002"]                # by entropy 0.80 > 0.75 > 0.72
    assert [h["accession"] for h in clicked_cos_first["hits"]] == expect == [h["accession"] for h in clicked_ent_first["hits"]]
    assert clicked_cos_first["metrics"][0] == ENT
    # the real scorer also returns entropy-sorted hits for [COS, ENT]
    monkeypatch.undo()
    real = _search(CAFFEINE, metrics=[COS, ENT], threshold=-1)
    ent = [h["scores"][ENT] for h in real["hits"]]
    assert ent == sorted(ent, reverse=True)


def test_unknown_or_empty_metric_selection_is_an_error(libs):
    with pytest.raises(ValueError, match="Unknown metric"):
        _search(CAFFEINE, metrics=["Nope"])
    with pytest.raises(ValueError, match="at least one metric"):
        _search(CAFFEINE, metrics=[])
    assert _run(**_ms(metrics=["Nope"])).status_code == 400
    assert _run(**_ms(metrics=[])).status_code == 400
    assert _run(**_ms(libraries=[])).status_code == 400
    assert _run(**_ms(libraries=["MoNA"])).status_code == 400


# ---------------------------------------------------------------- deep metrics
def test_deep_metrics_are_na_with_a_note_when_the_packages_are_missing(libs, monkeypatch):
    monkeypatch.setitem(sys.modules, "ms2deepscore", None)         # what the cloud server looks like: no torch / gensim
    monkeypatch.setitem(sys.modules, "gensim", None)
    res = _search(CAFFEINE, metrics=[ENT, mn.MS2DEEPSCORE, mn.SPEC2VEC], threshold=0.7)
    assert res["n_kept"] > 0 and all(h["scores"][mn.MS2DEEPSCORE] is None and h["scores"][mn.SPEC2VEC] is None for h in res["hits"])
    assert "not installed" in res["notes"][mn.MS2DEEPSCORE] and "not installed" in res["notes"][mn.SPEC2VEC]
    r = _run(**_ms(metrics=[ENT, mn.MS2DEEPSCORE, mn.SPEC2VEC]))
    body = r.json()
    assert r.status_code == 200
    assert any(w.startswith(f"{mn.MS2DEEPSCORE}: N/A") for w in body["warnings"]) and any(w.startswith(f"{mn.SPEC2VEC}: N/A") for w in body["warnings"])
    row = body["outputs"]["hits"]["rows"][0]
    assert row[-3:-1] == ["—", "—"]                                 # the two deep columns


def test_deep_metrics_are_na_for_ei_and_spec2vec_for_negative_mode(libs):
    res = _search(EI_QUERY, kind="EI", metrics=[ENT, mn.MS2DEEPSCORE, mn.SPEC2VEC], threshold=0.5)
    assert res["hits"] and all(h["scores"][mn.MS2DEEPSCORE] is None for h in res["hits"])
    assert "EI" in res["notes"][mn.MS2DEEPSCORE] or "not EI" in res["notes"][mn.MS2DEEPSCORE]
    neg = _search("110.0713 300\n135.0451 999\n179.0350 500", metrics=[ENT, mn.SPEC2VEC], precursor_mz=179.0350, positive=False, threshold=-1)
    assert "positive-mode" in neg["notes"][mn.SPEC2VEC] and all(h["scores"][mn.SPEC2VEC] is None for h in neg["hits"])


# ---------------------------------------------------------------- EI
def test_ei_ignores_the_precursor_uses_half_a_dalton_and_never_mixes_with_ms_ms(libs):
    assert ls.EI_TOLERANCE == 0.5
    res = _search(EI_QUERY, kind="EI", libraries=("MassBank", "GNPS"), metrics=[ENT, COS], threshold=0.5, tolerance=0.02, precursor_mz=999.0)
    assert res["tolerance"] == 0.5 and res["kind"] == "EI"
    accs = [h["accession"] for h in res["hits"]]
    assert all("-E00" in a or a == "CCMSLIB00000FIX011" for a in accs) and "MSBNK-FIX-0001" not in accs     # EI only
    assert {accs[0], accs[1]} == {"MSBNK-FIX-E001", "CCMSLIB00000FIX011"}                                    # the exact spectrum from both libraries
    assert all(h["precursor_mz"] is None and h["ppm"] is None for h in res["hits"])
    assert res["n_candidates"] == 6 and res["prefilter"] == {"total": 6, "kept": 6, "dropped": False}       # 5 MassBank + 1 GNPS EI spectra
    # no MS/MS spectrum shows up in an MS/MS search either
    ms = _search(CAFFEINE, libraries=("MassBank", "GNPS"), threshold=-1)
    assert all("-E00" not in h["accession"] and h["accession"] != "CCMSLIB00000FIX011" for h in ms["hits"])


def test_ei_fragments_match_within_half_a_dalton(libs):
    near = "\n".join(f"{float(m) + 0.3} {i}" for m, i in (ln.split() for ln in EI_QUERY.splitlines()))
    far = "\n".join(f"{float(m) + 0.6} {i}" for m, i in (ln.split() for ln in EI_QUERY.splitlines()))
    top = lambda text: next(h for h in _search(text, kind="EI", metrics=[ENT, mn.BINNED_COSINE], threshold=-1)["hits"] if h["accession"] == "MSBNK-FIX-E001")
    assert top(near)["scores"][ENT] == pytest.approx(1.0) and top(near)["matches"] == 7
    assert top(far)["scores"][ENT] < 0.01 and top(far)["matches"] == 0


def test_ei_na_metrics_need_a_precursor(libs):
    res = _search(EI_QUERY, kind="EI", metrics=[ENT, mn.MOD_COSINE_GREEDY, mn.NEUTRAL_LOSSES, COS], threshold=0.9)
    assert res["hits"][0]["scores"][mn.MOD_COSINE_GREEDY] is None and res["hits"][0]["scores"][mn.NEUTRAL_LOSSES] is None
    assert "precursor" in res["notes"][mn.MOD_COSINE_GREEDY] and mn.NEUTRAL_LOSSES in res["notes"]


def test_ei_step_one_entropy_is_ms_entropys_own_score(libs):
    import ms_entropy
    lib = lx.load("MassBank")
    q = parse_spectrum(EI_QUERY)
    q_mz, q_int = sim.clean_peaks(*sim.normalize_peaks_input(q))
    ids, ent = ls._entropy_step(np.column_stack([q_mz, q_int]), lib, 0.5)
    assert len(ids) == 5 and ids.tolist() == lib.ei_ids().tolist()
    for i, score in zip(ids, ent):
        mz, rel = lib.peaks(i)
        direct = ms_entropy.calculate_entropy_similarity(np.column_stack([q_mz, q_int]), np.column_stack([mz, rel.astype(np.float64) * 100.0]), ms2_tolerance_in_da=0.5)
        assert score == direct
    # ... and the SpectralEntropy (weighted) metric scored in step two is the same number
    res = _search(EI_QUERY, kind="EI", metrics=[ENT], threshold=-1)
    by_row = {h["row"]: h["scores"][ENT] for h in res["hits"]}
    assert all(by_row[int(i)] == pytest.approx(s, abs=1e-12) for i, s in zip(ids, ent))


def test_ei_prefilter_keeps_the_top_entropy_spectra_and_says_so(libs, monkeypatch):
    default = ls.EI_PREFILTER_N
    assert default == 5000
    monkeypatch.setattr(ls, "EI_PREFILTER_N", 2)
    res = _search(EI_QUERY, kind="EI", metrics=[COS], threshold=-1)
    assert res["n_candidates"] == 2 and res["prefilter"] == {"total": 5, "kept": 2, "dropped": True}
    assert {h["accession"] for h in res["hits"]} == {"MSBNK-FIX-E001", "MSBNK-FIX-E002"}       # the two caffeine EI spectra have the highest entropy similarity
    assert any("top 2" in w and "weighted entropy similarity" in w and "low entropy" not in w for w in res["warnings"])
    body = _run(spectrum_type="EI", ms2=_paste(EI_QUERY), libraries=["MassBank"], metrics=[COS], threshold="-1").json()
    assert "top 2 of 5" in body["outputs"]["summary"]["values"]["Candidates scored"]["note"]
    monkeypatch.setattr(ls, "EI_PREFILTER_N", default)                                        # everything fits: exact, and no warning
    res = _search(EI_QUERY, kind="EI", metrics=[COS], threshold=-1)
    assert res["n_candidates"] == 5 and not res["prefilter"]["dropped"] and not any("entropy" in w for w in res["warnings"])


def test_ei_runner_warns_about_a_given_precursor_and_drops_the_ms_ms_columns(libs):
    r = _run(spectrum_type="EI", ms2=_paste(EI_QUERY), precursor=_prec("194.0"), libraries=["MassBank"], tol_da="0.02", threshold="0.5")
    body = r.json()
    assert r.status_code == 200 and any("precursor m/z is not used for EI" in w for w in body["warnings"])
    hits = body["outputs"]["hits"]
    assert hits["columns"] == ["Rank", "Library", "Accession", "Name", "Matched peaks", ENT, COS, "Instrument / CE"]
    assert all(len(row) == len(hits["columns"]) for row in hits["rows"])
    r = _run(spectrum_type="EI", ms2=_paste(EI_QUERY), libraries=["MassBank"])             # no precursor needed, none warned about
    assert r.status_code == 200 and not any("not used" in w for w in r.json()["warnings"])


def test_gnps_has_no_ei_spectra_message(libs, monkeypatch):
    # make the GNPS fixture index look like the real one (no EI)
    lib = lx.load("GNPS")
    monkeypatch.setattr(lib, "n_msms", lib.n)
    res = _search(EI_QUERY, kind="EI", libraries=("MassBank", "GNPS"), metrics=[ENT], threshold=0.5)
    assert "GNPS has no EI spectra." in res["warnings"] and res["n_candidates"] == 5


# ---------------------------------------------------------------- runner
def test_runner_output_shapes_for_ms_ms(libs):
    r = _run(**_ms(libraries=["MassBank", "GNPS"], metrics=[COS, mn.NIST_DOT, ENT], threshold="0.9"))
    assert r.status_code == 200
    body = r.json()
    out = body["outputs"]
    assert set(out) == {"summary", "hits", "mirror"} and body["status"] == "ok"
    s = out["summary"]["values"]
    assert list(s) == ["Best hit", "Rank score", "Hits ≥ threshold", "Candidates scored", "Libraries"]
    assert s["Best hit"]["value"] == "Caffeine" and s["Best hit"]["note"].startswith("MassBank MSBNK-FIX-0001") or s["Best hit"]["note"].startswith("GNPS")
    assert s["Rank score"] == {"value": "1.0000", "note": ENT}                       # first metric of the list order, whatever the click order
    assert s["Libraries"] == {"value": "MassBank · GNPS", "note": "fixture-1 · fixture-1"}
    assert s["Candidates scored"]["note"] == "within 10 ppm" and int(s["Candidates scored"]["value"]) == 9
    hits = out["hits"]
    assert hits["columns"] == ["Rank", "Library", "Accession", "Name", "Adduct", "Precursor m/z", "Δ ppm", "Matched peaks", ENT, COS, mn.NIST_DOT, "Instrument / CE"]
    assert hits["count"] == int(s["Hits ≥ threshold"]["value"]) and len(hits["rows"]) == hits["count"] <= 100
    for row in hits["rows"]:
        assert len(row) == len(hits["columns"])
    first = hits["rows"][0]
    assert first[0] == "1" and first[1] in ("MassBank", "GNPS") and first[3] == "Caffeine" and first[4] == "[M+H]+" and first[5] == "195.0877"
    assert first[8] == {"badge": "1.0000", "tone": "good"}                           # >= threshold -> good badge
    mb = next(r for r in hits["rows"] if r[1] == "MassBank")
    assert mb[2] == {"text": "MSBNK-FIX-0001", "url": "https://massbank.eu/MassBank/RecordDisplay?id=MSBNK-FIX-0001"}
    gn = next(r for r in hits["rows"] if r[1] == "GNPS")
    assert gn[2]["url"] == f"https://gnps.ucsd.edu/ProteoSAFe/gnpslibraryspectrum.jsp?SpectrumID={gn[2]['text']}"
    assert mb[-1].startswith("LC-ESI-QFT · Q Exactive") and mb[-1].endswith("· 45")
    assert "plotly" in out["mirror"] and out["mirror"]["plotly"]["data"]
    assert body["verdict"].startswith("Best hit <b>Caffeine</b> (") and "entropy 1.0000" in body["verdict"] and "hits ≥ 0.9" in body["verdict"]


def test_below_threshold_scores_are_plain_text_and_the_table_is_capped_at_100(libs):
    body = _run(**_ms(metrics=[ENT, COS], threshold="-1")).json()
    cells = [c for row in body["outputs"]["hits"]["rows"] for c in row[8:10]]
    assert all(isinstance(c, dict) and c["tone"] == "good" for c in cells)           # threshold -1: everything is a hit
    body = _run(**_ms(metrics=[mn.BINNED_PEARSON], threshold="0.999")).json()
    rows = body["outputs"]["hits"]["rows"]
    assert rows and all(isinstance(row[8], dict) for row in rows)
    # count pill is the total kept, rows are at most 100
    many = ls.search(parse_spectrum(CAFFEINE), "MS/MS", ["MassBank"], [ENT], -1.0, precursor_mz=195.0877, top_n=2)
    assert many["n_kept"] == 6 and len(many["hits"]) == 2


def test_no_hits_gives_an_empty_table_a_sentence_and_no_mirror(libs):
    body = _run(**_ms(ms2=_paste("60.0 999\n61.0 500"), threshold="0.9", libraries=["MassBank"])).json()
    out = body["outputs"]
    assert out["hits"]["rows"] == [] and out["hits"]["count"] == 0 and "mirror" not in out
    assert out["summary"]["values"]["Best hit"] == "—" and out["summary"]["values"]["Hits ≥ threshold"]["value"] == "0"
    assert body["verdict"].startswith("No hit: none of the 6 candidates reaches 0.9")
    body = _run(**_ms(precursor=_prec("500.0"))).json()
    assert body["verdict"] == "No library spectrum has a matching precursor m/z, so nothing was scored." and body["outputs"]["summary"]["values"]["Candidates scored"]["value"] == "0"
    assert any("No library spectrum within 10 ppm" in w for w in body["warnings"])


def test_runner_adduct_and_ion_mode(libs):
    body = _run(**_ms(precursor=_prec("195.0877", "Positive", "[M+NH4]+"), libraries=["MassBank"], threshold="-1")).json()
    assert [row[2]["text"] for row in body["outputs"]["hits"]["rows"]] == ["MSBNK-FIX-0005"]
    body = _run(**_ms(precursor=_prec("195.0877", "Negative", "Auto"), libraries=["MassBank"], threshold="-1")).json()
    assert [row[2]["text"] for row in body["outputs"]["hits"]["rows"]] == ["MSBNK-FIX-0007"]
    r = _run(**_ms(precursor=_prec("195.0877", "Positive", "[M-H]-")))
    assert r.status_code == 400 and "does not match" in r.json()["detail"]
    assert _run(**_ms(tol_ppm="0")).status_code == 400
    assert _run(**_ms(threshold="2")).status_code == 400
    assert _run(spectrum_type="IR", ms2=_paste(CAFFEINE)).status_code == 400


# ---------------------------------------------------------------- missing libraries
def test_nothing_installed_names_the_setup_command(no_libs):
    r = _run(**_ms())
    assert r.status_code == 502
    detail = r.json()["detail"]
    assert "MassBank is not installed" in detail and "GNPS is not installed" in detail and "python -m mzdesk.ms2.library_setup MassBank" in detail
    r = _run(**_ms(libraries=["GNPS"]))
    assert "python -m mzdesk.ms2.library_setup GNPS" in r.json()["detail"]


def test_one_missing_library_is_searched_around_with_a_warning(massbank_only):
    body = _run(**_ms(libraries=["MassBank", "GNPS"], threshold="0.9")).json()
    assert body["status"] == "ok" and body["outputs"]["hits"]["rows"]
    assert any("GNPS is not installed" in w and "library_setup GNPS" in w and "Searching MassBank only" in w for w in body["warnings"])
    assert body["outputs"]["summary"]["values"]["Libraries"]["value"] == "MassBank"
    assert {row[1] for row in body["outputs"]["hits"]["rows"]} == {"MassBank"}


def test_an_old_or_broken_index_asks_for_a_rebuild(massbank_only):
    meta_file = lx.index_dir("MassBank") / "meta.json"
    import json
    meta = json.loads(meta_file.read_text())
    meta_file.write_text(json.dumps({**meta, "format": 0}))
    lx.clear_cache()
    with pytest.raises(lx.LibraryMissing, match="--rebuild"):
        lx.load("MassBank")
    meta_file.write_text(json.dumps(meta))
    (lx.index_dir("MassBank") / "peak_mz.npy").unlink()
    lx.clear_cache()
    with pytest.raises(lx.LibraryMissing, match="--force"):
        lx.load("MassBank")


# ---------------------------------------------------------------- ranking metric
def test_ranking_skips_a_first_metric_that_is_na_for_every_hit(libs):
    # EI has no precursor: ModifiedCosineGreedy is N/A for every spectrum, so the NIST dot product ranks (and is the Rank score)
    res = _search(EI_QUERY, kind="EI", libraries=("MassBank", "GNPS"), metrics=[mn.NIST_DOT, mn.MOD_COSINE_GREEDY], threshold=0.5)
    assert res["metrics"] == [mn.MOD_COSINE_GREEDY, mn.NIST_DOT] and res["rank_metric"] == mn.NIST_DOT
    scores = [h["scores"][mn.NIST_DOT] for h in res["hits"]]
    assert len(scores) >= 3 and scores == sorted(scores, reverse=True) and scores[0] == pytest.approx(1.0)
    # a different first metric that does apply is the ranking metric
    assert _search(EI_QUERY, kind="EI", metrics=[mn.MOD_COSINE_GREEDY, ENT], threshold=0.5)["rank_metric"] == ENT
    assert _search(EI_QUERY, kind="EI", metrics=[ENT, mn.MOD_COSINE_GREEDY], threshold=0.5)["rank_metric"] == ENT
    # MS/MS: both apply, the first of the list order ranks
    assert _search(CAFFEINE, metrics=[mn.NIST_DOT, COS], threshold=0.5)["rank_metric"] == COS               # CosineGreedy comes before the NIST dot product in the list
    assert _search(CAFFEINE, metrics=[ENT, COS], threshold=0.99999999)["rank_metric"] == ENT


def test_runner_rank_score_card_and_verdict_use_the_ranking_metric(libs):
    body = _run(spectrum_type="EI", ms2=_paste(EI_QUERY), libraries=["MassBank", "GNPS"], metrics=[mn.NIST_DOT, mn.MOD_COSINE_GREEDY], threshold="0.5").json()
    s = body["outputs"]["summary"]["values"]
    assert s["Rank score"] == {"value": "1.0000", "note": mn.NIST_DOT}                      # not "—"
    assert s["Best hit"] != "—" and mn.NIST_DOT in body["verdict"] and "1.0000" in body["verdict"]
    assert any(w.startswith(f"Ranked by {mn.NIST_DOT}, the first selected metric that applies ({mn.MOD_COSINE_GREEDY} is N/A for every hit)") for w in body["warnings"])
    rows = body["outputs"]["hits"]["rows"]
    cols = body["outputs"]["hits"]["columns"]
    j = cols.index(mn.NIST_DOT)
    val = lambda c: float(c["badge"] if isinstance(c, dict) else c)
    nist = [val(r[j]) for r in rows]
    assert nist == sorted(nist, reverse=True)
    # nothing special when the first metric applies
    body = _run(**_ms(metrics=[ENT, COS])).json()
    assert not any(w.startswith("Ranked by") for w in body["warnings"]) and body["outputs"]["summary"]["values"]["Rank score"]["note"] == ENT


# ---------------------------------------------------------------- GNPS copies of MassBank records
def test_gnps_parser_flags_copies_of_massbank_records():
    got = {r["accession"]: r["copy_of"] for r, _ in lp.iter_records("GNPS", GNPS_COPIES) if r}
    assert got == {"CCMSLIB00000COPY01": "0001", "CCMSLIB00000COPY02": "0002", "CCMSLIB00000COPY03": "ZZ9999", "CCMSLIB00000COPY04": None,
                   "CCMSLIB00000COPY05": "0009", "CCMSLIB00000COPY06": ""}             # MassbankEU: prefix too; ORGANISM only -> flagged without an ID
    assert all(r["copy_of"] is None for r, _ in lp.iter_records("GNPS", GNPS) if r)
    assert all(r["copy_of"] is None for r, _ in lp.iter_records("MassBank", MASSBANK) if r)


def test_copy_flags_short_ids_and_the_massbank_lookup(copies_libs):
    assert lx.short_id("MSBNK-LCSB-LU119903") == "LU119903" and lx.short_id("MSBNK-FIX-0001") == "0001" and lx.short_id("odd") == "odd"
    gn, mb = lx.load("GNPS"), lx.load("MassBank")
    assert gn.meta["counts"]["massbank_copies"] == 5 and int(np.asarray(gn.copy_flag).sum()) == 5
    flagged = np.flatnonzero(np.asarray(gn.copy_flag) == 1)
    assert sorted(i["copy_of"] for i in gn.info(flagged)) == ["", "0001", "0002", "0009", "ZZ9999"]
    rows = mb.short_id_rows(["0001", "ZZ9999"])
    assert [r["accession"] for r in mb.info(rows["0001"])] == ["MSBNK-FIX-0001"] and rows["ZZ9999"] == []


def test_copies_are_hidden_only_when_the_massbank_original_is_a_candidate(copies_libs):
    both = _search(CAFFEINE, libraries=("MassBank", "GNPS"), metrics=[ENT], threshold=-1)
    accs = {h["accession"] for h in both["hits"]}
    assert both["hidden_copies"] == 2 and both["n_candidates"] == 10 and len(both["hits"]) == 10
    assert {"CCMSLIB00000COPY01", "CCMSLIB00000COPY02"}.isdisjoint(accs)                     # copies of 0001 / 0002, whose originals are in the results
    assert {"MSBNK-FIX-0001", "MSBNK-FIX-0002"} <= accs                                      # the MassBank originals stay
    assert {"CCMSLIB00000COPY03", "CCMSLIB00000COPY04", "CCMSLIB00000COPY05", "CCMSLIB00000COPY06"} <= accs   # original not in MassBank / not a candidate / no ID / not a copy
    assert [l["candidates"] for l in both["libraries"]] == [6, 4]
    # GNPS alone, or MassBank missing: nothing is hidden
    gnps_only = _search(CAFFEINE, libraries=("GNPS",), metrics=[ENT], threshold=-1)
    assert gnps_only["hidden_copies"] == 0 and gnps_only["n_candidates"] == 6
    # the original is filtered out by the adduct -> its copy stays
    nh4 = _search(CAFFEINE, libraries=("MassBank", "GNPS"), metrics=[ENT], threshold=-1, adduct="[M+NH4]+")
    assert nh4["hidden_copies"] == 0 and {h["accession"] for h in nh4["hits"]} == {"MSBNK-FIX-0005"}
    # EI has no copies
    assert _search(EI_QUERY, kind="EI", libraries=("MassBank", "GNPS"), metrics=[ENT], threshold=-1)["hidden_copies"] == 0


def test_runner_says_how_many_copies_were_hidden(copies_libs):
    body = _run(**_ms(libraries=["MassBank", "GNPS"], metrics=[ENT], threshold="-1")).json()
    note = body["outputs"]["summary"]["values"]["Candidates scored"]
    assert note["value"] == "10" and note["note"] == "within 10 ppm · 2 GNPS copies of MassBank spectra hidden"
    body = _run(**_ms(libraries=["GNPS"], metrics=[ENT], threshold="-1")).json()
    assert body["outputs"]["summary"]["values"]["Candidates scored"]["note"] == "within 10 ppm"


# ---------------------------------------------------------------- verdict
def test_verdict_sentences():
    out = {"summary": {"values": {"Best hit": {"value": "Caffeine", "note": "MassBank MSBNK-X"}, "Rank score": {"value": "0.9312", "note": ENT},
                                  "Hits ≥ threshold": {"value": "4", "note": "any of 2 metrics ≥ 0.7"}, "Candidates scored": {"value": "85"}}}}
    assert verdicts.build("ms2_library", {}, out) == "Best hit <b>Caffeine</b> (MassBank MSBNK-X) · entropy 0.9312 · 4 hits ≥ 0.7."
    out["summary"]["values"]["Hits ≥ threshold"] = {"value": "1", "note": "any of 2 metrics ≥ 0.5"}
    out["summary"]["values"]["Rank score"] = {"value": "0.8", "note": COS}
    assert verdicts.build("ms2_library", {}, out) == f"Best hit <b>Caffeine</b> (MassBank MSBNK-X) · {COS} 0.8 · 1 hit ≥ 0.5."
    out["summary"]["values"]["Hits ≥ threshold"] = {"value": "0", "note": "any of 2 metrics ≥ 0.7"}
    assert verdicts.build("ms2_library", {}, out) == "No hit: none of the 85 candidates reaches 0.7 on any selected metric."
    out["summary"]["values"]["Candidates scored"] = {"value": "0"}
    assert "nothing was scored" in verdicts.build("ms2_library", {"spectrum_type": "MS/MS"}, out)
    assert verdicts.build("ms2_library", {"spectrum_type": "EI"}, out) == "The selected libraries have no EI spectra."
    assert verdicts.build("ms2_library", {}, {}) is None


def test_verdict_builder_sits_right_after_ms2_esi():
    names = list(verdicts.BUILDERS)
    assert names[names.index("ms2_esi") + 1] == "ms2_library"


# ---------------------------------------------------------------- setup (offline)
def test_setup_downloads_builds_and_records_the_source(tmp_path, monkeypatch):
    import shutil
    monkeypatch.setenv("MZDESK_LIBRARY_DIR", str(tmp_path))
    lx.clear_cache()
    info = {"url": "https://example.test/MassBank_NISTformat.msp", "file": setup_mod.MASSBANK_ASSET, "version": "2099.01", "date": "2099-01-02",
            "size": MASSBANK.stat().st_size}
    calls = []
    monkeypatch.setattr(setup_mod, "source_info", lambda name: info)
    monkeypatch.setattr(setup_mod, "_download", lambda url, dest, size=None: calls.append(url) or shutil.copy(MASSBANK, dest))
    meta = setup_mod.setup("MassBank")
    assert calls == [info["url"]] and meta["version"] == "2099.01" and meta["date"] == "2099-01-02" and meta["source_url"] == info["url"]
    assert meta["n_spectra"] == 15 and (tmp_path / "MassBank" / "meta.json").is_file()
    assert (tmp_path / "raw" / setup_mod.MASSBANK_ASSET).is_file() and (tmp_path / "raw" / (setup_mod.MASSBANK_ASSET + ".json")).is_file()
    assert lx.load("MassBank").label() == "MassBank 2099.01"
    # already installed: nothing is downloaded or rebuilt
    assert setup_mod.setup("MassBank")["version"] == "2099.01" and len(calls) == 1
    # --rebuild reuses the downloaded file without touching the network
    monkeypatch.setattr(setup_mod, "source_info", lambda name: pytest.fail("no network on --rebuild"))
    assert setup_mod.setup("MassBank", rebuild=True)["n_spectra"] == 15 and len(calls) == 1
    setup_mod.setup("MassBank", rebuild=True, delete_raw=True)
    assert not (tmp_path / "raw" / setup_mod.MASSBANK_ASSET).exists()
    lx.clear_cache()


def test_setup_rejects_unknown_libraries():
    with pytest.raises(SystemExit, match="Unknown library"):
        setup_mod.setup("MoNA")


def test_a_failed_build_leaves_the_old_index_alone(tmp_path, monkeypatch):
    monkeypatch.setenv("MZDESK_LIBRARY_DIR", str(tmp_path))
    lx.clear_cache()
    _build(("MassBank",))
    good = lx.load("MassBank").meta["built"]
    empty = tmp_path / "empty.msp"
    empty.write_text("")
    with pytest.raises(ValueError, match="no usable spectra"):
        setup_mod.build_from("MassBank", empty, {"version": "x"}, progress=lambda *_: None)
    assert lx.load("MassBank").meta["built"] == good and not (tmp_path / ".MassBank.building").exists()
    lx.clear_cache()


def test_no_applicable_metric_is_not_a_no_hit(libs, monkeypatch):
    monkeypatch.setitem(sys.modules, "ms2deepscore", None)
    monkeypatch.setitem(sys.modules, "gensim", None)
    res = _search(CAFFEINE, metrics=[mn.MS2DEEPSCORE], threshold=0.7)
    assert res["n_candidates"] > 0 and res["n_kept"] == 0 and res["n_metrics_scored"] == 0
    assert "N/A" in res["notes"][mn.MS2DEEPSCORE]
    r = _run(**_ms(metrics=[mn.MS2DEEPSCORE]))
    body = r.json()
    assert body["outputs"]["summary"]["values"]["Hits ≥ threshold"]["note"].endswith("no applicable metric")
    assert "No hit could be scored" in body["verdict"]


def test_ei_ignores_hidden_msms_only_values(libs):
    r = _run(spectrum_type="EI", ms2=_paste(EI_QUERY), precursor=_prec("0"), tol_ppm="0", tol_da="10", threshold="0.5")
    assert r.status_code == 200, r.text
    body = r.json()
    assert any("precursor m/z is not used for EI" in w for w in body["warnings"])
    assert "Adduct" not in body["outputs"]["hits"]["columns"]


def test_on_an_intel_mac_ms2deepscore_says_why_it_is_left_out(monkeypatch):
    """Intel Macs stop at ms2deepscore 2.5.5, which cannot read the model ('model_params', CI 2026-11-13): say so instead of failing.
    sys.platform itself is left alone: torch loads its Windows DLLs only when it reads "win32" (a first import here failed on Windows)."""
    import types
    fake = types.ModuleType("ms2deepscore")
    fake.MS2DeepScore = object
    fake.vector_operations = types.ModuleType("ms2deepscore.vector_operations")
    fake.vector_operations.cosine_similarity_matrix = None
    monkeypatch.setitem(sys.modules, "ms2deepscore", fake)
    monkeypatch.setitem(sys.modules, "ms2deepscore.vector_operations", fake.vector_operations)
    monkeypatch.setattr(sim, "_this_machine", lambda: ("darwin", "x86_64"))
    monkeypatch.setattr(sim, "_get_ms2deepscore_model", lambda: pytest.fail("no model download on an Intel Mac"))
    res = sim.compute_all_similarities([(100, 1), (150, 3), (200, 5)], [(100, 2), (151, 3), (200, 4)], 250, 251, use_deep=True)
    row = next(r for r in res if r["method"] == mn.MS2DEEPSCORE)
    assert row["score"] is None and row["note"] == f"skipped: {sim.MS2DS_INTEL_MAC}"
    notes = {}
    ls._ms2deepscore_batch(np.zeros(1), [100.0], [1.0], 250.0, "positive", [([100.0], [1.0], 250.0)], notes)
    assert notes[mn.MS2DEEPSCORE] == f"N/A: {sim.MS2DS_INTEL_MAC}"
    monkeypatch.setattr(sim, "_this_machine", lambda: ("darwin", "arm64"))
    assert sim.ms2deepscore_unsupported() is None
