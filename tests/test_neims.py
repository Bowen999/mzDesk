"""NEIMS（EI 质谱预测）的测试：指纹和谱对照原代码的参考值，runner 的输出，错误输入。

参考值 tests/fixtures/neims_reference.json 来自原代码 + 原权重在 Docker 里的运行（validation/neims/）。
权重没安装（python -m mzdesk.ms2.neims_setup）时，需要权重的测试跳过。
"""
import json
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from mzdesk.api.common import EngineError, InputError
from mzdesk.api.runners import ei_ms
from mzdesk.ms2 import neims_engine as ne
from mzdesk.ms2 import neims_setup
from mzdesk.shared.spectrum import parse_spectrum

REF = json.loads((Path(__file__).parent / "fixtures" / "neims_reference.json").read_text())["molecules"]
CAFFEINE = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
needs_weights = pytest.mark.skipif(not Path(ne.weights_path()).is_file(),
                                   reason="NEIMS weights not installed (python -m mzdesk.ms2.neims_setup)")


def _smiles_input(smiles):
    return {"structure": {"mode": "SMILES / InChI", "value": smiles}}


def _dense(sparse):
    a = np.zeros(ne.FP_LENGTH, dtype=np.float32)
    a[sparse["idx"]] = sparse["cnt"]
    return a


# ---------------------------------------------------------------- features and spectrum vs the original
@pytest.mark.parametrize("ref", REF, ids=[r["name"] for r in REF])
def test_features_match_original(ref):
    mol = Chem.MolFromSmiles(ref["smiles"])
    assert np.array_equal(ne.fingerprint(mol), _dense(ref["fp"]))
    assert ne.model_mass(mol) == ref["model_mass"]


@needs_weights
@pytest.mark.parametrize("ref", REF, ids=[r["name"] for r in REF])
def test_spectrum_matches_original(ref):
    expected = np.zeros(ne.N_BINS)
    expected[ref["mz"]] = ref["value"]
    got = ne.spectrum_array(Chem.MolFromSmiles(ref["smiles"]))
    assert got.max() == pytest.approx(999.0, rel=1e-5)
    np.testing.assert_allclose(got, expected, atol=0.01)         # 999 刻度；float32 数值噪声约 0.002


# ---------------------------------------------------------------- engine
def test_nominal_mass_versus_model_mass():
    caffeine = Chem.MolFromSmiles(CAFFEINE)
    assert ne.nominal_mass(caffeine) == ne.model_mass(caffeine) == 194
    assert ne.nominal_mass(Chem.MolFromSmiles("Clc1ccccc1")) == 112         # 35Cl，不是 114
    stearin = Chem.MolFromSmiles("CCCCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCCCC")
    assert (ne.nominal_mass(stearin), ne.model_mass(stearin)) == (890, 891)  # 精确质量 890.83 四舍五入比标称质量多 1


@needs_weights
def test_predict_caffeine():
    res = ne.predict(CAFFEINE)
    mz = [m for m, _ in res["peaks"]]
    rel = dict(res["peaks"])
    assert res["formula"] == "C8H10N4O2" and res["molecular_ion_mz"] == res["model_mass"] == 194 and res["warnings"] == []
    assert mz == sorted(mz) and max(rel.values()) == pytest.approx(100.0, rel=1e-5)
    assert min(rel.values()) >= ne.MIN_REL and rel[194] == pytest.approx(100.0, rel=1e-5)
    assert all(rel.get(m, 0) > 20 for m in (109, 82, 67, 55))     # 咖啡因的主要碎片
    assert max(mz) <= 194 + ne.ABOVE_MASS


@needs_weights
def test_no_peak_above_molecular_ion_window():
    res = ne.predict("CCO")
    assert res["molecular_ion_mz"] == 46 and max(m for m, _ in res["peaks"]) <= 46 + ne.ABOVE_MASS


@needs_weights
def test_largest_fragment_is_used_for_salts():
    salt = ne.predict(CAFFEINE + ".Cl")
    assert any("several components" in w for w in salt["warnings"])
    assert salt["peaks"] == ne.predict(CAFFEINE)["peaks"] and salt["formula"] == "C8H10N4O2"


@needs_weights
def test_warnings_for_charge_rare_elements_and_size():
    assert any("net charge" in w for w in ne.predict("C[N+](C)(C)C")["warnings"])
    assert any("Ga" in w for w in ne.predict("C[Ga](C)C")["warnings"])
    assert not ne.predict("[O-][N+](=O)c1ccccc1")["warnings"]       # 硝基是电荷平衡的，不算带电
    assert any("heavy atoms" in w for w in ne._prepare("C" * 101)[1])


def test_charged_input_is_neutralised():
    mol, warnings = ne._prepare("[O-]C(=O)c1ccccc1.[Na+]")                  # sodium benzoate -> benzoic acid
    assert Chem.MolToSmiles(mol) == "O=C(O)c1ccccc1" and Chem.GetFormalCharge(mol) == 0
    assert len(warnings) == 2 and "several components" in warnings[0]
    assert "(-1) was neutralized" in warnings[1] and "O=C(O)c1ccccc1" in warnings[1]
    mol, warnings = ne._prepare("CC[NH3+].[Cl-]")                           # the cation is neutralised to the amine
    assert Chem.MolToSmiles(mol) == "CCN" and any("(+1) was neutralized" in w for w in warnings)
    assert not any("neutralized" in w for w in ne._prepare("[O-][N+](=O)c1ccccc1")[1])        # the nitro group is not a charge
    mol, warnings = ne._prepare("C[N+](C)(C)C")                             # a quaternary ammonium stays charged: only the net charge warning
    assert Chem.GetFormalCharge(mol) == 1 and len(warnings) == 1 and "net charge of +1" in warnings[0]


@needs_weights
def test_charged_input_predicts_the_neutral_molecule():
    res = ne.predict("[O-]C(=O)c1ccccc1.[Na+]")
    acid = ne.predict("OC(=O)c1ccccc1")
    assert res["smiles"] == acid["smiles"] and res["formula"] == "C7H6O2" and res["molecular_ion_mz"] == 122
    assert res["peaks"] == acid["peaks"] and acid["warnings"] == []
    assert not any("net charge" in w for w in res["warnings"])


@needs_weights
@pytest.mark.parametrize("smiles", ["", "   ", "C1CC", "not a molecule"])
def test_invalid_smiles_is_a_value_error(smiles):
    with pytest.raises(ValueError, match="Invalid"):
        ne.predict(smiles)


@needs_weights
def test_too_heavy_is_a_value_error():
    with pytest.raises(ValueError, match="above the model's maximum"):
        ne.predict("C" * 80)                                        # C80H162，M = 1122


def test_missing_weights_error_names_the_setup_command(tmp_path, monkeypatch):
    monkeypatch.setenv("MZDESK_NEIMS_WEIGHTS", str(tmp_path / "none.npz"))
    with pytest.raises(ne.NEIMSError, match="neims_setup"):
        ne.predict(CAFFEINE)


def test_incomplete_weights_file_is_reported(tmp_path, monkeypatch):
    bad = tmp_path / "bad.npz"
    np.savez(bad, **{"dense.kernel": np.zeros((2, 2), dtype=np.float32)})
    monkeypatch.setenv("MZDESK_NEIMS_WEIGHTS", str(bad))
    with pytest.raises(ne.NEIMSError, match="unreadable"):
        ne.predict(CAFFEINE)


def test_setup_paths_and_checksum(tmp_path, monkeypatch):
    monkeypatch.setenv("MZDESK_NEIMS_DIR", str(tmp_path))
    assert neims_setup.default_weights() == tmp_path / "neims_weights.npz"
    f = tmp_path / "x"
    f.write_bytes(b"abc")
    assert neims_setup.sha256_of(f) == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    assert len(neims_setup.ZIP_SHA256) == 64


# ---------------------------------------------------------------- runner
@needs_weights
def test_runner_fills_every_block():
    out = ei_ms.run_ei(_smiles_input(CAFFEINE))
    blocks = out["outputs"]
    assert set(blocks) == {"summary", "spectrum", "peaks", "formats"} and out["warnings"] == []
    json.dumps(out)
    assert blocks["summary"]["values"] == {
        "Ion": {"value": "[M]+•", "note": "70 eV EI"}, "Formula": "C8H10N4O2",
        "Molecular ion m/z": {"value": "194", "note": "nominal"}, "Base peak m/z": "194"}
    rows = blocks["peaks"]["rows"]
    mz = [int(r[0]) for r in rows]
    assert mz == sorted(mz) and len(set(mz)) == len(mz) and ["194", "100.0"] in rows
    cards = blocks["spectrum"]["peaks"]
    assert len(cards) == len(rows) and blocks["spectrum"]["plotly"]["data"]
    assert [c["kind"] for c in cards].count("Molecular ion [M]+•") == 1
    assert next(c for c in cards if c["mz"] == 194)["kind"] == "Molecular ion [M]+•"
    fmts = blocks["formats"]["formats"]
    assert "MSP" in fmts and "IONMODE: Positive" in fmts["MSP"]
    assert "NAME: C8H10N4O2" in fmts["MSP"]                       # SMILES input -> the molecular formula names the record
    assert len(parse_spectrum(fmts["Peak list"])) == len(rows)


@needs_weights
def test_runner_uses_the_typed_name_in_formats(monkeypatch):
    monkeypatch.setattr("mzdesk.identity.resolve.to_smiles", lambda s: CAFFEINE)
    out = ei_ms.run_ei({"structure": {"mode": "Name / ID / Formula", "value": "caffeine"}})
    assert "NAME: caffeine" in out["outputs"]["formats"]["formats"]["MSP"]


@needs_weights
def test_runner_passes_engine_warnings_through():
    out = ei_ms.run_ei(_smiles_input(CAFFEINE + ".Cl"))
    assert len(out["warnings"]) == 1 and "largest" in out["warnings"][0]


@needs_weights
@pytest.mark.parametrize("smiles", ["", "C1CC", "C" * 80])
def test_runner_input_errors(smiles):
    with pytest.raises(InputError):
        ei_ms.run_ei(_smiles_input(smiles))


def test_runner_reports_missing_weights_as_engine_error(tmp_path, monkeypatch):
    monkeypatch.setenv("MZDESK_NEIMS_WEIGHTS", str(tmp_path / "none.npz"))
    with pytest.raises(EngineError, match="neims_setup"):
        ei_ms.run_ei(_smiles_input(CAFFEINE))
