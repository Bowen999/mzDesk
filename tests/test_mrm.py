"""MRM transitions page: ranking / roles, CE conversion, the method-table CSV and the three paths.

Offline: the engine unit tests and the fake CFM-ID result. Real runs are marked: the predicted
path needs Docker with the CFM-ID image (`needs_docker` pattern of tests/test_api.py); the library
path needs the MassBank index (`python -m mzdesk.ms2.library_setup MassBank`).
"""
import pytest
from fastapi.testclient import TestClient

from mzdesk import spec
from mzdesk.api import verdicts
from mzdesk.api.app import app
from mzdesk.api.common import InputError
from mzdesk.api.runners import mrm as runner
from mzdesk.ms2 import cfmid_engine, library_index
from mzdesk.ms2 import mrm as engine

client = TestClient(app)


def _module():
    return next(module for module in spec.MODULES if module["id"] == "mrm")


def _run(inputs):
    return client.post("/api/modules/mrm/run", json=inputs)


def _paste(peaks):
    return {"source": "Paste", "value": peaks}


def _my_spectrum(peaks, precursor=195.0877):
    return _paste(f"BEGIN IONS\nPEPMASS={precursor}\n{peaks}\nEND IONS")


def _docker_unavailable():
    """Why the CFM-ID image cannot run here (same check as tests/test_api.py), or None."""
    try:
        cfmid_engine.check_docker()
    except cfmid_engine.CFMIDError as exc:
        return str(exc)
    return None


NO_DOCKER = _docker_unavailable()
NO_MASSBANK = not library_index.is_installed("MassBank")


# ---------------------------------------------------------------- fake engine result
def _fake_cfmid(smiles, adduct):
    """Two products plus the precursor per energy: 138 is strongest at 20 eV, 42 at 40 eV."""
    precursor = (195.0877, "C8H11N4O2+ · Cn1cnc2c1c(=O)n(C)c(=O)n2C")
    product_138 = (138.0662, "C6H8N3O+ · C=[N+](C=O)C1=CN(C)C=N1")
    product_110 = (110.0713, "C5H8N3+ · C=[NH+]C1=CN(C)C=N1")
    product_42 = (42.0338, "C2H4N+ · C=[N+]=C")
    spectra = [
        {"energy": "low", "ev": 10, "peaks": [(product_138[0], 30.0), (product_110[0], 5.0), (precursor[0], 100.0)],
         "fragments": [product_138[1], product_110[1], precursor[1]]},
        {"energy": "medium", "ev": 20, "peaks": [(product_138[0], 100.0), (product_110[0], 20.0), (precursor[0], 80.0)],
         "fragments": [product_138[1], product_110[1], precursor[1]]},
        {"energy": "high", "ev": 40, "peaks": [(product_42[0], 90.0), (product_110[0], 50.0), (precursor[0], 70.0)],
         "fragments": [product_42[1], product_110[1], precursor[1]]},
    ]
    return {"engine": "CFM-ID 4.4.10", "smiles": smiles, "adduct": adduct, "precursor_mz": 195.08765,
            "ion_formula": "C8H11N4O2+", "spectra": spectra, "warnings": [], "log": ""}


STRUCTURE_INPUTS = {"start_from": "Structure (predicted fragments)",
                    "structure": {"mode": "SMILES / InChI", "value": spec.CAFFEINE},
                    "ion_mode": "Positive", "n_transitions": 3, "ce_unit": "eV"}


# ---------------------------------------------------------------- engine: annotations
def test_fragment_annotation_parsing():
    assert engine.fragment_smiles("C2H4N+ · C=[N+]=C") == "C=[N+]=C"
    assert engine.fragment_smiles("C2H4N+ · C=[N+]=C (+1 more)") == "C=[N+]=C"
    assert engine.fragment_smiles("C=[N+]=C") == "C=[N+]=C"
    assert engine.fragment_smiles("") is None
    assert engine.fragment_formula("C6H8N3O+ · C=[N+](C=O)C1=CN(C)C=N1") == "C6H8N3O+"
    assert engine.fragment_formula("C2H4N+ · C=[N+]=C") == "C2H4N+"
    assert engine.fragment_formula("not a smiles") is None


def test_structure_key_uses_the_largest_part():
    assert engine.structure_key(spec.CAFFEINE) == "RYYVLZVUVIJVGH"
    assert engine.structure_key("CC(=O)[O-].[Na+]") == "QTBSBXVTEAMEQO"
    with pytest.raises(ValueError, match="Invalid SMILES"):
        engine.structure_key("not a smiles")


# ---------------------------------------------------------------- engine: merging and ranking
def test_merge_keeps_the_most_intense_energy():
    spectra = [
        {"ev": 10, "peaks": [(100.0, 10.0), (195.0, 100.0)], "fragments": ["a", "p"]},
        {"ev": 20, "peaks": [(100.0, 50.0), (120.0, 20.0)], "fragments": ["b", "c"]},
        {"ev": 40, "peaks": [(100.004, 30.0), (100.02, 99.0)], "fragments": ["d", "e"]},
    ]
    merged = engine.merge_energies(spectra)
    assert [(round(ion["mz"], 3), ion["intensity"], ion["ev"], ion["annotation"]) for ion in merged] == [
        (100.0, 50.0, 20, "b"), (100.02, 99.0, 40, "e"), (120.0, 20.0, 20, "c"), (195.0, 100.0, 10, "p")]


def test_merge_needs_peaks():
    with pytest.raises(ValueError, match="no predicted spectrum"):
        engine.merge_energies([])
    with pytest.raises(ValueError, match="no peaks"):
        engine.merge_energies([{"ev": 10, "peaks": [], "fragments": []}])


def test_ranking_assigns_quantifier_and_qualifiers():
    products = [{"mz": 110.0, "intensity": 5.0}, {"mz": 138.0, "intensity": 100.0}, {"mz": 83.0, "intensity": 20.0}]
    result = engine.select_transitions(products, 195.0, 3)
    assert [round(t["mz"], 3) for t in result["transitions"]] == [138.0, 83.0, 110.0]
    assert [t["role"] for t in result["transitions"]] == ["quantifier", "qualifier", "qualifier"]
    assert [round(t["relative"], 1) for t in result["transitions"]] == [100.0, 20.0, 5.0]
    assert result["n_precursor_skipped"] == 0 and result["n_requested"] == 3


def test_ranking_is_deterministic_on_ties_and_limits():
    products = [{"mz": 138.0, "intensity": 50.0}, {"mz": 110.0, "intensity": 50.0}, {"mz": 80.0, "intensity": 50.0}]
    assert [t["mz"] for t in engine.select_transitions(products, 195.0, 2)["transitions"]] == [80.0, 110.0]
    with pytest.raises(ValueError, match="At least one transition"):
        engine.select_transitions(products, 195.0, 0)


def test_ranking_drops_the_precursor_and_duplicate_mz():
    products = [{"mz": 195.0877, "intensity": 100.0}, {"mz": 195.0878, "intensity": 99.0},
                {"mz": 138.0, "intensity": 20.0}, {"mz": 138.005, "intensity": 15.0}, {"mz": 138.02, "intensity": 10.0}]
    result = engine.select_transitions(products, 195.0877, 3)
    assert [round(t["mz"], 3) for t in result["transitions"]] == [138.0, 138.02]
    assert result["n_precursor_skipped"] == 2
    with pytest.raises(ValueError, match="No product ion below the precursor"):
        engine.select_transitions([{"mz": 195.0877, "intensity": 10.0}], 195.0877, 1)


# ---------------------------------------------------------------- engine: CE and CSV
@pytest.mark.parametrize("ev, expected", [(10, 26), (20, 51), (40, 103)])
def test_nce_conversion(ev, expected):
    assert engine.nce_from_ev(ev, 195.0877) == expected
    assert engine.nce_from_ev(ev, 193.0731) == int(round(ev * 500 / 193.0731))
    assert engine.nce_from_ev(None, 195.0877) is None


def test_nce_conversion_rejects_bad_values():
    with pytest.raises(ValueError, match="positive precursor"):
        engine.nce_from_ev(20, 0)
    with pytest.raises(ValueError, match="positive"):
        engine.nce_from_ev(0, 195.0877)


def test_ce_labels():
    assert engine.ce_label(20, 195.0877, "eV") == "20 eV"
    assert engine.ce_label(20, 195.0877, "NCE") == "51 NCE"
    assert engine.ce_label(None, 195.0877, "eV") is None
    with pytest.raises(ValueError, match="Unknown CE unit"):
        engine.ce_label(20, 195.0877, "V")


def test_method_table_csv_quotes_when_needed():
    assert engine.method_table_csv(["A", "B"], [["x,y", 'q"z']]) == 'A,B\n"x,y","q""z"'
    assert engine.method_table_csv(["A"], []) == "A"


def test_columns_match_the_spec_table():
    table = next(block for block in _module()["outputs"] if block["type"] == "table")
    assert table["columns"] == runner.COLUMNS


# ---------------------------------------------------------------- the example (My spectrum)
def test_example_is_the_my_spectrum_path():
    response = _run(_module()["example"])
    assert response.status_code == 200, response.text
    body = response.json()
    assert body["status"] == "ok"
    out = body["outputs"]
    assert out["summary"]["values"] == {"Polarity": "Positive", "Transitions": "3"}
    assert out["transitions"]["rows"] == [
        ["195.0877", "138.0663", "—", "Quantifier", "—", "My spectrum", "100.0 %"],
        ["195.0877", "110.0713", "—", "Qualifier", "—", "My spectrum", "10.1 %"],
        ["195.0877", "69.0447", "—", "Qualifier", "—", "My spectrum", "1.8 %"]]
    assert out["method_table"]["formats"]["Method table"].splitlines()[0] == ",".join(runner.COLUMNS)
    assert "pasted MS2 spectrum" in " ".join(out["notes"]["items"])
    assert body["verdict"] == "3 transitions for [M+H]+ <b>195.0877</b> → <b>138.0663</b> (quantifier)."


def test_my_spectrum_reads_the_header_and_ignores_a_structure():
    body = _run({"start_from": "My spectrum", "ion_mode": "Positive", "n_transitions": 2,
                 "structure": {"mode": "SMILES / InChI", "value": spec.CAFFEINE},
                 "ms2": _my_spectrum("195.0878 80\n138.0663 100\n110.0713 40\n69.0447 5")}).json()
    rows = body["outputs"]["transitions"]["rows"]
    assert [row[1] for row in rows] == ["138.0663", "110.0713"]
    assert rows[0][0] == "195.0877" and rows[0][6] == "100.0 %" and rows[1][6] == "40.0 %"
    assert any("structure is ignored" in warning for warning in body["warnings"])
    assert any("read from the spectrum header" in warning for warning in body["warnings"])


def test_relative_intensity_uses_the_spectrum_base_peak():
    # the precursor peak is the base peak of the pasted spectrum, so the quantifier is 50 %
    body = _run({"start_from": "My spectrum", "n_transitions": 2,
                 "ms2": _my_spectrum("195.0877 1000\n138.0663 500\n110.0713 100")}).json()
    rows = body["outputs"]["transitions"]["rows"]
    assert [row[1] for row in rows] == ["138.0663", "110.0713"]
    assert [row[6] for row in rows] == ["50.0 %", "10.0 %"]
    assert rows[0][3] == "Quantifier"


def test_fewer_products_than_requested_warns():
    body = _run({"start_from": "My spectrum", "ms2": _my_spectrum("100 100\n150 50", precursor=200),
                 "n_transitions": 3}).json()
    assert len(body["outputs"]["transitions"]["rows"]) == 2
    assert any("Fewer product ions" in warning for warning in body["warnings"])


def test_my_spectrum_without_a_header_is_400():
    response = _run({"start_from": "My spectrum", "ms2": _paste("138.0663 100\n110.0713 40")})
    assert response.status_code == 400
    assert "PEPMASS" in response.json()["detail"]


# ---------------------------------------------------------------- the structure path (fake CFM-ID)
def test_structure_path_merges_energies(monkeypatch):
    monkeypatch.setattr(cfmid_engine, "run_cfmid", _fake_cfmid)
    result = runner.run_mrm(STRUCTURE_INPUTS)
    rows = result["outputs"]["transitions"]["rows"]
    assert rows == [
        ["195.0877", "138.0662", "20 eV", "Quantifier", "C6H8N3O+", "CFM-ID", "100.0 %"],
        ["195.0877", "42.0338", "40 eV", "Qualifier", "C2H4N+", "CFM-ID", "90.0 %"],
        ["195.0877", "110.0713", "40 eV", "Qualifier", "C5H8N3+", "CFM-ID", "50.0 %"]]
    assert any("precursor ion is not a product" in warning for warning in result["warnings"])
    csv_text = result["outputs"]["method_table"]["formats"]["Method table"]
    assert csv_text.splitlines()[0] == ",".join(runner.COLUMNS)
    assert csv_text.splitlines()[1] == "195.0877,138.0662,20 eV,Quantifier,C6H8N3O+,CFM-ID,100.0 %"
    assert result["outputs"]["notes"]["items"][0].startswith("Rows come from CFM-ID")
    assert any("model output" in item for item in result["outputs"]["notes"]["items"])
    assert verdicts.build("mrm", STRUCTURE_INPUTS, result["outputs"]) == \
        "3 transitions for [M+H]+ <b>195.0877</b> → <b>138.0662</b> (quantifier), CE <b>20 eV</b>."


def test_structure_path_nce_units(monkeypatch):
    monkeypatch.setattr(cfmid_engine, "run_cfmid", _fake_cfmid)
    rows = runner.run_mrm({**STRUCTURE_INPUTS, "ce_unit": "NCE"})["outputs"]["transitions"]["rows"]
    assert rows[0][2] == "51 NCE" and rows[1][2] == "103 NCE"


def test_negative_mode_uses_m_h_minus(monkeypatch):
    def fake(smiles, adduct):
        return {**_fake_cfmid(smiles, adduct), "adduct": adduct,
                "spectra": [{**spectrum, "fragments": []} for spectrum in _fake_cfmid(smiles, adduct)["spectra"]]}
    monkeypatch.setattr(cfmid_engine, "run_cfmid", fake)
    result = runner.run_mrm({**STRUCTURE_INPUTS, "ion_mode": "Negative"})
    assert result["outputs"]["transitions"]["rows"][0][0] == "193.0731"
    assert "for [M-H]-" in verdicts.build("mrm", {}, result["outputs"])
    assert "[M-H]-" in " ".join(result["outputs"]["notes"]["items"])


@pytest.mark.docker
@pytest.mark.skipif(NO_DOCKER is not None, reason=f"Docker unavailable: {NO_DOCKER}")
def test_real_cfmid_structure_path():
    result = runner.run_mrm(STRUCTURE_INPUTS)
    rows = result["outputs"]["transitions"]["rows"]
    assert len(rows) == 3
    assert all(row[0] == "195.0877" and row[5] == "CFM-ID" for row in rows)
    assert all(row[2].endswith(" eV") for row in rows)
    assert all(float(row[1]) < 195.0877 for row in rows)
    assert rows[0][3] == "Quantifier" and rows[0][4] != "—" and rows[0][6] == "100.0 %"


# ---------------------------------------------------------------- the library path
@pytest.mark.skipif(NO_MASSBANK, reason="MassBank index not installed (python -m mzdesk.ms2.library_setup MassBank)")
def test_library_path_picks_the_same_compound_and_adduct():
    result = runner.run_mrm({"start_from": "Library spectrum",
                             "structure": {"mode": "SMILES / InChI", "value": spec.CAFFEINE},
                             "ion_mode": "Positive", "n_transitions": 3, "ce_unit": "eV"})
    rows = result["outputs"]["transitions"]["rows"]
    assert len(rows) == 3
    assert all(row[5].startswith("MassBank MSBNK-") for row in rows)
    assert rows[0][1] == "138.0662" and rows[0][3] == "Quantifier" and rows[0][4] == "—"
    assert rows[0][2] == "60"                      # the [M+H]+ record; the MALDI record has no CE
    assert rows[0][6] == "100.0 %"
    assert any("different InChIKey" in warning for warning in result["warnings"])
    assert "MassBank MSBNK-" in result["outputs"]["notes"]["items"][0]


def test_missing_massbank_index_says_how_to_install_it(monkeypatch):
    monkeypatch.setattr(library_index, "is_installed", lambda name: False)
    with pytest.raises(InputError) as error:
        runner.run_mrm({"start_from": "Library spectrum",
                        "structure": {"mode": "SMILES / InChI", "value": spec.CAFFEINE}, "ion_mode": "Positive"})
    assert str(error.value) == "MassBank index is not installed. Run: python -m mzdesk.ms2.library_setup MassBank"


# ---------------------------------------------------------------- bad input
@pytest.mark.parametrize("inputs", [
    {"start_from": "Structure (predicted fragments)", "ion_mode": "Positive"},
    {"start_from": "Library spectrum", "ion_mode": "Positive"},
    {"start_from": "My spectrum", "ion_mode": "Positive"},
    {"start_from": "Nope", "ion_mode": "Positive"},
    {"start_from": "My spectrum", "ion_mode": "Positive", "ms2": _my_spectrum("138.0663 100"), "n_transitions": 0},
    {"start_from": "My spectrum", "ion_mode": "Positive", "ms2": _my_spectrum("138.0663 100"), "n_transitions": 1.5},
    {"start_from": "My spectrum", "ion_mode": "Positive", "ms2": _my_spectrum("138.0663 100"), "ce_unit": "V"},
    {"start_from": "My spectrum", "ion_mode": "Both", "ms2": _my_spectrum("138.0663 100")},
])
def test_bad_input_is_400(inputs):
    response = _run(inputs)
    assert response.status_code == 400
    assert response.json()["detail"]


def test_missing_structure_does_not_touch_docker(monkeypatch):
    def boom(*args, **kwargs):
        raise AssertionError("CFM-ID must not run without a structure")
    monkeypatch.setattr(cfmid_engine, "run_cfmid", boom)
    with pytest.raises(InputError, match="needs a structure"):
        runner.run_mrm({"start_from": "Structure (predicted fragments)"})
