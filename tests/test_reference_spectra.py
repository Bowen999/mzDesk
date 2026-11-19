"""Reference Spectra (Spectrum · Reference Spectra): a structure -> the recorded MassBank / GNPS spectra of that compound.

Offline, on the same two fixture libraries as tests/test_library_search.py, indexed into a temporary MZDESK_LIBRARY_DIR.
Caffeine in the fixtures: MassBank MSBNK-FIX-0001 / 0002 (MS/MS, [M+H]+) and E001 / E002 (EI); GNPS CCMSLIB00000FIX001.
The other caffeine-named decoys carry no InChIKey, so a structure lookup must not find them.
"""
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from mzdesk import spec
from mzdesk.api import verdicts
from mzdesk.api.app import app
from mzdesk.ms2 import library_index as lx
from mzdesk.ms2 import library_setup as setup_mod
from mzdesk.ms2 import metric_names as mn
from mzdesk.ms2 import reference_lookup as rl

FIX = Path(__file__).parent / "fixtures"
MASSBANK, GNPS = FIX / "library_massbank.msp", FIX / "library_gnps.mgf"
CAFFEINE = {"mode": "SMILES / InChI", "value": spec.CAFFEINE}
CAFFEINE_MSMS = {"MSBNK-FIX-0001", "MSBNK-FIX-0002", "CCMSLIB00000FIX001"}
CAFFEINE_EI = {"MSBNK-FIX-E001", "MSBNK-FIX-E002"}
client = TestClient(app)


def _run(**inputs):
    return client.post("/api/modules/ms2_reference/run", json=inputs)


def _accessions(res):
    return {r["accession"] for r in res["records"]}


@pytest.fixture(scope="module")
def libs(tmp_path_factory):
    mp = pytest.MonkeyPatch()
    mp.setenv("MZDESK_LIBRARY_DIR", str(tmp_path_factory.mktemp("libs")))
    lx.clear_cache()
    for name in ("MassBank", "GNPS"):
        setup_mod.build_from(name, MASSBANK if name == "MassBank" else GNPS,
                             {"url": "file://fixture", "version": "fixture-1", "date": "2026-01-01"}, progress=lambda *_: None)
    yield
    mp.undo()
    lx.clear_cache()


@pytest.fixture
def no_libs(monkeypatch, tmp_path):
    monkeypatch.setenv("MZDESK_LIBRARY_DIR", str(tmp_path / "nothing"))
    lx.clear_cache()
    yield
    lx.clear_cache()


def test_page_is_in_the_spectrum_group_after_library_search():
    ids = [m["id"] for m in spec.MODULES if m["subgroup"] == "ms2"]
    assert ids.index("ms2_reference") == ids.index("ms2_library") + 1
    mod = next(m for m in spec.MODULES if m["id"] == "ms2_reference")
    assert mod["live"] and mod["example"] == {"structure": {"mode": "Name / ID / Formula", "value": "caffeine"}}
    assert {f["id"] for f in mod["inputs"]} == {"structure", "kinds", "ion_mode", "libraries", "ms2", "tol_da", "plot"}


def test_structure_rows_match_the_first_inchikey_block_only(libs):
    lib = lx.load("MassBank")
    rows = lib.structure_rows("RYYVLZVUVIJVGH")
    assert {r["accession"] for r in lib.info(rows)} == {"MSBNK-FIX-0001", "MSBNK-FIX-0002"} | CAFFEINE_EI
    assert lib.structure_rows("rYYVLZVUVIJVGH").size == rows.size          # case does not matter
    assert lib.structure_rows("").size == 0 and lib.structure_rows("RYYVLZVUVIJVGH-UHFFFAOYSA-N").size == 0


def test_caffeine_finds_every_record_with_its_inchikey_and_no_decoy(libs):
    res = rl.lookup(spec.CAFFEINE, ["MassBank", "GNPS"])
    assert _accessions(res) == CAFFEINE_MSMS | CAFFEINE_EI
    assert res["key"] == "RYYVLZVUVIJVGH-UHFFFAOYSA-N" and res["key14"] == "RYYVLZVUVIJVGH"
    assert [(l["name"], l["found"]) for l in res["libraries"]] == [("MassBank", 4), ("GNPS", 1)]
    assert all(r["same_stereo"] and r["n_peaks"] > 0 and not r["scores"] for r in res["records"])
    # sorted: library, MS/MS before EI
    assert [r["library"] for r in res["records"]] == ["MassBank"] * 4 + ["GNPS"]
    assert [r["kind"] for r in res["records"][:4]] == ["MS/MS", "MS/MS", "EI", "EI"]
    ei = next(r for r in res["records"] if r["kind"] == "EI")
    assert ei["precursor_mz"] is None and ei["adduct"] == ""


def test_filters_spectrum_type_ion_mode_and_library(libs):
    assert _accessions(rl.lookup(spec.CAFFEINE, ["MassBank", "GNPS"], kinds=["EI"])) == CAFFEINE_EI
    assert _accessions(rl.lookup(spec.CAFFEINE, ["MassBank", "GNPS"], kinds=["MS/MS"])) == CAFFEINE_MSMS
    # the ion-mode filter applies to MS/MS only: EI spectra stay
    assert _accessions(rl.lookup(spec.CAFFEINE, ["MassBank", "GNPS"], ion_mode="Negative")) == CAFFEINE_EI
    assert _accessions(rl.lookup(spec.CAFFEINE, ["GNPS"])) == {"CCMSLIB00000FIX001"}
    with pytest.raises(ValueError, match="spectrum type"):
        rl.lookup(spec.CAFFEINE, ["MassBank"], kinds=[])
    with pytest.raises(ValueError, match="Ion mode"):
        rl.lookup(spec.CAFFEINE, ["MassBank"], ion_mode="Both")


def test_a_salt_is_looked_up_by_its_largest_part(libs):
    assert _accessions(rl.lookup(spec.CAFFEINE + ".Cl", ["MassBank", "GNPS"])) == CAFFEINE_MSMS | CAFFEINE_EI


def test_my_spectrum_scores_every_reference_and_ranks_by_entropy(libs):
    query = [tuple(map(float, line.split())) for line in spec.CAFFEINE_MASSBANK_MS2.splitlines()]     # = MSBNK-FIX-0001
    res = rl.lookup(spec.CAFFEINE, ["MassBank", "GNPS"], query=query)
    best = res["records"][0]
    assert best["accession"] == "MSBNK-FIX-0001" and best["scores"][mn.ENTROPY_WEIGHTED] == pytest.approx(1.0)
    entropies = [r["scores"][mn.ENTROPY_WEIGHTED] for r in res["records"] if r["scores"].get(mn.ENTROPY_WEIGHTED) is not None]
    assert entropies == sorted(entropies, reverse=True) and len(entropies) == len(res["records"])
    assert all(set(r["scores"]) == set(rl.METRICS) and r["matches"] is not None for r in res["records"])


def test_page_lists_the_spectra_plots_the_first_and_gives_the_formats(libs):
    r = _run(structure=CAFFEINE)
    assert r.status_code == 200, r.text
    body = r.json()
    out = body["outputs"]
    cards = out["summary"]["values"]
    assert cards["Compound"] == {"value": "Caffeine", "note": "RYYVLZVUVIJVGH-UHFFFAOYSA-N"}
    assert cards["Reference spectra"] == {"value": "5", "note": "3 MS/MS · 2 EI"}
    assert cards["Libraries"]["value"] == "MassBank 4 · GNPS 1" and cards["Best match"] == "—"
    rows = out["spectra"]["rows"]
    assert out["spectra"]["count"] == 5 and len(rows) == 5 and "Entropy" not in out["spectra"]["columns"]
    assert rows[0]["cells"][2]["url"] == "https://massbank.eu/MassBank/RecordDisplay?id=MSBNK-FIX-0001"
    assert rows[0]["reference"]["library"] == "MassBank" and rows[0]["reference"]["accession"] == "MSBNK-FIX-0001"
    shown = out["shown"]
    assert shown["query"] == {"inchikey": "RYYVLZVUVIJVGH-UHFFFAOYSA-N", "compound": "Caffeine", "tol_da": 0.02, "peaks": None, "precursor_mz": None}
    assert shown["record"]["accession"] == "MSBNK-FIX-0001" and shown["record"]["url"] == rows[0]["cells"][2]["url"]
    assert len(shown["peak_list"]) == shown["record"]["n_peaks"] and "spectrum" not in out and "formats" not in out
    assert "plotly" in shown["spectrum"] and "Peak list" in shown["formats"]["formats"]
    assert "MSBNK-FIX-0001" in str(shown["formats"])
    assert body["verdict"].startswith("<b>5 recorded spectra</b>: 3 MS/MS · 2 EI (MassBank 4 · GNPS 1)")


def test_page_with_my_spectrum_ranks_scores_and_draws_a_mirror(libs):
    r = _run(structure=CAFFEINE, ms2={"source": "Paste", "value": spec.CAFFEINE_MASSBANK_MS2}, kinds=["MS/MS"])
    assert r.status_code == 200, r.text
    out = r.json()["outputs"]
    assert out["spectra"]["columns"][-3:] == ["Entropy", "Cosine", "Matched peaks"]
    first = out["spectra"]["rows"][0]["cells"]
    assert first[2]["text"] == "MSBNK-FIX-0001" and first[-3]["tone"] == "good"
    assert out["shown"]["query"]["peaks"] and out["shown"]["query"]["tol_da"] == 0.02
    assert out["shown"]["record"]["scores"][mn.ENTROPY_WEIGHTED] == pytest.approx(1.0) and "Up: my spectrum" in out["shown"]["spectrum"]["help"][0]
    assert out["summary"]["values"]["Best match"]["value"] == "1.0000"
    assert "Best match with your spectrum" in r.json()["verdict"]


def test_plot_record_picks_a_row_and_an_unknown_one_says_so(libs):
    r = _run(structure=CAFFEINE, plot="msbnk-fix-e002")
    assert r.status_code == 200 and r.json()["outputs"]["shown"]["record"]["accession"] == "MSBNK-FIX-E002"
    r = _run(structure=CAFFEINE, plot="NOPE-1")
    assert any("NOPE-1 is not among the listed spectra" in w for w in r.json()["warnings"])


def test_a_compound_without_spectra_says_so(libs):
    r = _run(structure={"mode": "SMILES / InChI", "value": "CCCCCCCCCCCCCCCCO"})
    assert r.status_code == 200
    body = r.json()
    assert body["outputs"]["summary"]["values"]["Reference spectra"]["value"] == "0"
    assert "shown" not in body["outputs"]
    assert any("No MS/MS or EI spectrum with InChIKey block" in w and "Library Search finds them" in w for w in body["warnings"])
    assert body["verdict"] == "No recorded spectrum of this compound in MassBank or GNPS."


def test_bad_inputs_are_400(libs):
    assert _run(structure=CAFFEINE, kinds=[]).status_code == 400
    assert _run(structure=CAFFEINE, kinds=["MS3"]).status_code == 400
    assert _run(structure=CAFFEINE, libraries=["MoNA"]).status_code == 400
    assert _run(structure=CAFFEINE, ion_mode="Both").status_code == 400
    assert _run(structure={"mode": "SMILES / InChI", "value": ""}).status_code == 400


def test_nothing_installed_names_the_setup_command(no_libs):
    r = _run(structure=CAFFEINE)
    assert r.status_code == 502 and "python -m mzdesk.ms2.library_setup MassBank" in r.json()["detail"]


# ---------------------------------------------------------------- a clicked row: the Selected spectrum block fetches it (POST /api/reference/spectrum)
def _open(**body):
    return client.post("/api/reference/spectrum", json=body)


def _rows(**inputs):
    return {r["reference"]["accession"]: r["reference"] for r in _run(structure=CAFFEINE, **inputs).json()["outputs"]["spectra"]["rows"]}


def test_a_clicked_row_gives_its_record_plot_peaks_and_formats(libs):
    ref = _rows()["MSBNK-FIX-0001"]
    r = _open(**ref, inchikey="RYYVLZVUVIJVGH-UHFFFAOYSA-N", compound="Caffeine")
    assert r.status_code == 200, r.text
    body = r.json()
    rec = body["record"]
    assert rec["accession"] == "MSBNK-FIX-0001" and rec["kind"] == "MS/MS" and rec["adduct"] == "[M+H]+" and rec["same_stereo"]
    assert rec["url"] == "https://massbank.eu/MassBank/RecordDisplay?id=MSBNK-FIX-0001" and rec["smiles"]
    assert rec["scores"] == {} and rec["matches"] is None
    assert len(body["peak_list"]) == rec["n_peaks"] > 0 and max(p[1] for p in body["peak_list"]) == 100
    assert "plotly" in body["spectrum"] and "Peak list" in body["formats"]["formats"]
    assert _open(**ref).json()["record"]["same_stereo"] is False          # no query InChIKey: nothing to compare with


def test_a_clicked_row_with_my_spectrum_is_scored_and_mirrored(libs):
    query = [[float(v) for v in line.split()] for line in spec.CAFFEINE_MASSBANK_MS2.splitlines()]
    body = _open(**_rows()["MSBNK-FIX-0001"], peaks=query, precursor_mz=195.0877, tol_da=0.02).json()
    assert body["record"]["scores"][mn.ENTROPY_WEIGHTED] == pytest.approx(1.0) and body["record"]["matches"] > 0
    assert "Up: my spectrum" in body["spectrum"]["help"][0]
    ei = _open(**_rows()["MSBNK-FIX-E001"], peaks=query).json()["record"]           # EI: scored with the EI scorer, no precursor
    assert ei["kind"] == "EI" and ei["precursor_mz"] is None and set(ei["scores"]) == set(rl.METRICS)


def test_a_clicked_row_checks_its_request(libs):
    ref = _rows()["MSBNK-FIX-0001"]
    assert _open(**{**ref, "accession": "MSBNK-OTHER"}).status_code == 409          # the library changed between the run and the click
    assert _open(**{**ref, "row": 10 ** 9}).status_code == 409
    assert _open(**{**ref, "library": "MoNA"}).status_code == 400
    assert _open(**{**ref, "row": -1}).status_code == 400 and _open(**{**ref, "row": "3"}).status_code == 400
    assert _open(**ref, peaks=[[0, 10]]).status_code == 400 and _open(**ref, peaks="138 100").status_code == 400
    assert _open(**ref, tol_da=50).status_code == 400


def test_a_clicked_row_without_the_library_is_503(no_libs):
    r = _open(library="MassBank", row=0, accession="MSBNK-FIX-0001")
    assert r.status_code == 503 and "MassBank" in r.json()["detail"]


def test_verdict_builder_is_registered():
    assert "ms2_reference" in verdicts.BUILDERS
