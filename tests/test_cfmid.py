"""ESI MS/MS（CFM-ID 4 + MSRB）测试：解析器、加合物、runner 外形离线测；真实 Docker 测试在没有 Docker / 镜像时跳过。

fixtures：tests/fixtures/cfmid_*.log 是镜像里真实跑出来的输出（MSRB：PC 16:0/18:1 [M+H]+、心磷脂 [M-2H]2-；cfm-predict：咖啡因 [M+H]+）。
"""
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

from mzdesk.api.common import EngineError, InputError
from mzdesk.api.runners import esi_msms
from mzdesk.ms2 import cfmid_engine as cf
from mzdesk.shared.spectrum import first_record, parse_spectrum

FIX = Path(__file__).parent / "fixtures"
MSRB_PC = (FIX / "cfmid_msrb_pc.log").read_text()
MSRB_CL = (FIX / "cfmid_msrb_cl.log").read_text()
CFM_CAFFEINE = (FIX / "cfmid_cfm_caffeine.log").read_text()

CAFFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
PC = "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\\CCCCCCCC"      # PC 16:0/18:1
CARDIOLIPIN = next(ln for ln in MSRB_CL.splitlines() if ln.startswith("#SMILES=")).split("=", 1)[1]


@pytest.fixture(autouse=True)
def fresh_cache():
    cf._predict.cache_clear()


def _docker_unavailable():
    try:
        cf.check_docker()
    except cf.CFMIDError as exc:
        return str(exc)
    return None


NO_DOCKER = _docker_unavailable()


def real_docker(fn):
    """Marks a test as a real Docker run and skips it when Docker or the image is missing."""
    skip = pytest.mark.skipif(NO_DOCKER is not None, reason=f"Docker unavailable: {NO_DOCKER}")
    return pytest.mark.docker(skip(fn))


def _has_peak(res, mz, energy=0, tol=0.005):
    return any(abs(m - mz) < tol for m, _ in res["spectra"][energy]["peaks"])


# ---------------------------------------------------------------- adducts and structure
def test_adduct_list():
    assert cf.ADDUCTS[:2] == ["[M+H]+", "[M-H]-"]
    assert len(set(cf.ADDUCTS)) == len(cf.ADDUCTS) == 7
    assert set(cf.ML_ADDUCTS) <= set(cf.ADDUCTS) and set(cf._MSRB_NAME) <= set(cf.ADDUCTS)


@pytest.mark.parametrize("adduct, mz, ion", [
    ("[M+H]+", 195.0877, "C8H11N4O2+"), ("[M-H]-", 193.0731, "C8H9N4O2-"), ("[M+Na]+", 217.0696, "C8H10N4NaO2+"),
    ("[M+NH4]+", 212.1142, "C8H14N5O2+"), ("[M+Li]+", 201.0958, "C8H10LiN4O2+"), ("[M+Cl]-", 229.0498, "C8H10ClN4O2-"),
    ("[M-2H]2-", 96.0329, "[C8H8N4O2]2-")])
def test_precursor_and_ion_formula(adduct, mz, ion):
    from rdkit import Chem
    got_mz, got_ion = cf.ion_info(Chem.MolFromSmiles(CAFFEINE), adduct)
    assert got_mz == pytest.approx(mz, abs=1e-4) and got_ion == ion


def test_ion_formula_rejects_unsupported_elements():
    from rdkit import Chem
    with pytest.raises(ValueError, match="cannot form"):
        cf.ion_info(Chem.MolFromSmiles("C[Si](C)C"), "[M+H]+")


def test_prepare_keeps_largest_part_and_neutralises():
    smiles, warnings = cf.prepare("CC(=O)[O-].[Na+]")
    assert smiles == "CC(=O)O" and len(warnings) == 2 and "Ignored: [Na+]" in warnings[0]
    assert cf.prepare(CAFFEINE)[1] == []
    zwitterion = cf.prepare(PC)                          # a net-neutral zwitterion stays as it is
    assert "[N+](C)(C)C" in zwitterion[0] and zwitterion[1] == []


@pytest.mark.parametrize("smiles, message", [("", "Enter"), ("not a smiles", "Invalid SMILES"), ("C[N+](C)(C)C", "net charge")])
def test_prepare_rejects(smiles, message):
    with pytest.raises(ValueError, match=message):
        cf.prepare(smiles)


def test_unknown_adduct_is_a_user_error():
    with pytest.raises(ValueError, match="Unsupported adduct"):
        cf.run_cfmid(CAFFEINE, "[M+K]+")


# ---------------------------------------------------------------- parsers
def test_parse_msrb_log():
    p = cf.parse_log(MSRB_PC)
    assert p["adduct"] == "[M+H]+" and p["by"] == "MSRB 1.1.15"
    assert [len(e) for e in p["energies"]] == [1, 2, 12]
    assert all(max(y for _, y, _ in e) == 100.0 for e in p["energies"])
    mz, rel, frags = p["energies"][1][0]
    assert (mz, rel) == (184.07332, 66.7) and "[N+](C)(C)C" in frags[0]
    assert [e[0] for e in p["energies"][2]] == sorted(e[0] for e in p["energies"][2])


def test_parse_msrb_log_with_several_fragments_per_peak():
    p = cf.parse_log(MSRB_CL)
    assert p["adduct"] == "[M-2H](2H)-"                  # MSRB's own spelling of [M-2H]2-
    peak = next(x for x in p["energies"][1] if abs(x[0] - 602.34835) < 1e-5)
    assert len(peak[2]) == 2
    spectra = cf._spectra(p)
    i = [round(m, 5) for m, _ in spectra[1]["peaks"]].index(602.34835)
    assert spectra[1]["fragments"][i].endswith("(+1 more)")


def test_parse_cfm_log():
    p = cf.parse_log(CFM_CAFFEINE)
    assert p["adduct"] == "[M+H]+" and p["by"] == "CFM-ID 4.4.10"
    assert [len(e) for e in p["energies"]] == [2, 5, 10]
    low = dict((round(mz, 5), (y, f)) for mz, y, f in p["energies"][0])
    assert low[195.08765][0] == 100.0 and low[110.07127][0] == pytest.approx(1.0)
    assert len(low[195.08765][1]) == 9 and low[195.08765][1][0].startswith("Cn1c(=[OH+])")


def test_spectra_carry_energy_labels_and_annotation():
    spectra = cf._spectra(cf.parse_log(CFM_CAFFEINE))
    assert [(s["energy"], s["ev"]) for s in spectra] == [("low", 10), ("medium", 20), ("high", 40)]
    assert all(len(s["peaks"]) == len(s["fragments"]) for s in spectra)
    assert spectra[2]["fragments"][0] == "C2H4N+ · C=[N+]=C"          # m/z 42.0338


def test_fewer_than_three_energies_is_an_error():
    with pytest.raises(cf.CFMIDError, match="2 spectra"):
        cf._spectra(cf.parse_log("energy0\n100.0 100.0 0 (1)\nenergy1\n100.0 100.0 0 (1)\n\n0 100.0 C\n"))


def test_fragment_text():
    assert cf.fragment_text("C=[N+]=C") == "C2H4N+ · C=[N+]=C"
    assert cf.fragment_text("not smiles") == "not smiles"


# ---------------------------------------------------------------- docker plumbing (no Docker needed)
def test_docker_missing_or_down(monkeypatch):
    monkeypatch.setattr(cf.shutil, "which", lambda _: None)
    with pytest.raises(cf.CFMIDError) as err:
        cf.check_docker()
    assert str(err.value) == "Docker is not running. Start Docker Desktop and try again."

    monkeypatch.setattr(cf.shutil, "which", lambda _: "/usr/bin/docker")
    down = subprocess.CompletedProcess([], 1, "", "Cannot connect to the Docker daemon. Is the docker daemon running?")
    monkeypatch.setattr(cf.subprocess, "run", lambda *a, **k: down)
    with pytest.raises(cf.CFMIDError, match="Docker is not running"):
        cf.check_docker()


def test_missing_image_gives_the_pull_command(monkeypatch):
    monkeypatch.setattr(cf.shutil, "which", lambda _: "/usr/bin/docker")
    gone = subprocess.CompletedProcess([], 1, "", f"Error response from daemon: No such image: {cf.IMAGE}")
    monkeypatch.setattr(cf.subprocess, "run", lambda *a, **k: gone)
    with pytest.raises(cf.CFMIDError) as err:
        cf.check_docker()
    assert "docker pull --platform linux/amd64 wishartlab/cfmid:msml4.4.10_msrb1.1.15" in str(err.value)


def test_timeout_removes_the_container(monkeypatch, tmp_path):
    calls = []

    def fake_run(cmd, **kw):
        calls.append(cmd)
        if cmd[1] == "run":
            raise subprocess.TimeoutExpired(cmd, kw["timeout"])
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(cf.subprocess, "run", fake_run)
    with pytest.raises(cf.CFMIDError, match="did not finish"):
        cf._docker_run(["true"], tmp_path)
    run, rm = calls
    assert "--pull" in run and "never" in run and run[run.index("--name") + 1] == rm[-1] and rm[:3] == ["docker", "rm", "-f"]


class FakeDocker:
    """Stands in for _docker_run: writes the given tool outputs into the work dir and records the commands."""

    def __init__(self, monkeypatch, msrb=None, cfm=None):
        self.calls = []
        monkeypatch.setattr(cf, "check_docker", lambda: None)
        monkeypatch.setattr(cf, "_docker_run", self)
        self.msrb, self.cfm = msrb, cfm

    def __call__(self, args, workdir):
        self.calls.append(args)
        if args[0] == "java" and self.msrb:
            (workdir / "msrb_x").write_text(self.msrb)
        if args[0] == "cfm-predict" and self.cfm:
            (workdir / "cfm.log").write_text(self.cfm)
        return "tool output", 0


def test_msrb_is_tried_first_and_never_returns_another_adduct(monkeypatch):
    fake = FakeDocker(monkeypatch, msrb=MSRB_CL, cfm=CFM_CAFFEINE)       # MSRB answers with [M-2H](2H)-, not what was asked
    res = cf.run_cfmid(CAFFEINE, "[M+H]+")
    assert fake.calls[0][0] == "java" and "-na" in fake.calls[0] and fake.calls[0][fake.calls[0].index("-a") + 1] == "[M+H]+"
    assert fake.calls[1][0] == "cfm-predict" and fake.calls[1][3].endswith("[M+H]+/param_output.log")
    assert res["engine"] == "CFM-ID 4.4.10 (ML)" and res["adduct"] == "[M+H]+"


def test_msrb_result_wins_and_dianion_uses_msrb_spelling(monkeypatch):
    fake = FakeDocker(monkeypatch, msrb=MSRB_CL)
    res = cf.run_cfmid(CARDIOLIPIN, "[M-2H]2-")
    assert [c[0] for c in fake.calls] == ["java"] and fake.calls[0][fake.calls[0].index("-a") + 1] == "[M-2H](2H)-"
    assert res["engine"] == "MSRB 1.1.15" and _has_peak(res, 719.44754, tol=1e-4)
    assert res["precursor_mz"] == pytest.approx(719.4475, abs=1e-3) and res["ion_formula"].endswith("]2-")


def test_msrb_only_adduct_not_covered_is_a_user_error(monkeypatch):
    fake = FakeDocker(monkeypatch)                                       # MSRB writes nothing
    with pytest.raises(ValueError, match=r"MSRB rules do not cover this structure for \[M\+Na\]\+; .* only \[M\+H\]\+ and \[M-H\]-"):
        cf.run_cfmid(CAFFEINE, "[M+Na]+")
    assert [c[0] for c in fake.calls] == ["java"]                        # no ML fallback for this adduct


def test_tool_without_output_is_an_engine_error(monkeypatch):
    FakeDocker(monkeypatch)
    with pytest.raises(cf.CFMIDError, match="no spectrum"):
        cf.run_cfmid(CAFFEINE, "[M-H]-")


def test_results_are_cached_and_copies(monkeypatch):
    fake = FakeDocker(monkeypatch, cfm=CFM_CAFFEINE)
    first = cf.run_cfmid(CAFFEINE, "[M+H]+")
    first["spectra"].clear()
    second = cf.run_cfmid("CN1C=NC2=C1C(=O)N(C)C(=O)N2C", "[M+H]+")      # same molecule, other SMILES
    assert len(second["spectra"]) == 3 and len(fake.calls) == 2


# ---------------------------------------------------------------- runner
def _canned(adduct="[M+H]+"):
    spectra = cf._spectra(cf.parse_log(CFM_CAFFEINE))
    return {"engine": "CFM-ID 4.4.10 (ML)", "smiles": CAFFEINE, "adduct": adduct, "precursor_mz": 195.087652,
            "ion_formula": "C8H11N4O2+", "spectra": spectra, "warnings": ["a warning"], "log": ""}


def _run(monkeypatch, structure=None, adduct="[M+H]+"):
    monkeypatch.setattr(cf, "run_cfmid", lambda smiles, adduct: _canned(adduct))
    return esi_msms.run_esi({"structure": structure or {"mode": "SMILES / InChI", "value": CAFFEINE}, "adduct": adduct})


def test_runner_fills_every_block(monkeypatch):
    out = _run(monkeypatch)
    assert out["warnings"] == ["a warning"]
    o = out["outputs"]
    assert list(o) == ["summary", "spectrum_low", "spectrum_medium", "spectrum_high", "peaks", "formats"]
    assert o["summary"]["values"] == {
        "Engine": "CFM-ID 4.4.10 (ML)", "Precursor m/z": {"value": "195.0877", "note": "[M+H]+"}, "Ion formula": "C8H11N4O2+",
        "Peaks": {"value": "2 · 5 · 10", "note": "10 eV · 20 eV · 40 eV"}}
    for energy, n in (("low", 2), ("medium", 5), ("high", 10)):
        block = o[f"spectrum_{energy}"]
        assert block["plotly"]["data"] and len(block["peaks"]) == n
        assert block["peaks"][0]["rows"][-1][0] == "Fragment"
    rows = o["peaks"]["rows"]
    assert len(rows) == 17 and rows[0][:3] == ["low (10 eV)", "110.0713", "1.0"] and rows[-1][0] == "high (40 eV)"
    assert all(len(r) == 4 and r[3] for r in rows)


def test_runner_formats_hold_three_readable_records(monkeypatch):
    out = _run(monkeypatch, adduct="[M-H]-")
    formats = out["outputs"]["formats"]["formats"]
    assert list(formats) == ["MSP", "MGF", "MassBank", "SIRIUS .ms"]
    for name, text in formats.items():
        rec = first_record(text)
        first = parse_spectrum(rec["peaks_text"])
        assert rec["n_records"] == 3, name                                   # one record per energy
        assert rec["name"] == "C8H10N4O2 [M-H]- | low (10 eV)", name        # SMILES input -> molecular formula as the compound
        assert rec["precursor_mz"] == pytest.approx(195.08765), name
        assert [m for m, _ in first] == pytest.approx([110.07127, 195.08765]), name   # the first record only
    assert "IONMODE: Negative" in formats["MSP"] and "CHARGE=1-" in formats["MGF"]


def test_runner_names_records_after_a_typed_name(monkeypatch):
    monkeypatch.setattr(esi_msms.common, "structure", lambda inputs, label="", **kw: {"kind": "smiles", "smiles": CAFFEINE})
    out = _run(monkeypatch, structure={"mode": "Name / ID / Formula", "value": "caffeine"})
    assert first_record(out["outputs"]["formats"]["formats"]["MSP"])["name"] == "caffeine [M+H]+ | low (10 eV)"


def test_runner_error_mapping(monkeypatch):
    def broken(smiles, adduct):
        raise cf.CFMIDError("Docker is not running. Start Docker Desktop and try again.", "daemon log")

    monkeypatch.setattr(cf, "run_cfmid", broken)
    with pytest.raises(EngineError) as err:
        esi_msms.run_esi({"structure": {"mode": "SMILES / InChI", "value": CAFFEINE}, "adduct": "[M+H]+"})
    assert err.value.log == "daemon log"
    for inputs in ({"structure": {"mode": "SMILES / InChI", "value": "not a smiles"}, "adduct": "[M+H]+"},
                   {"structure": {"mode": "SMILES / InChI", "value": ""}, "adduct": "[M+H]+"},
                   {"structure": {"mode": "SMILES / InChI", "value": CAFFEINE}, "adduct": "[M+K]+"}):
        monkeypatch.undo()
        with pytest.raises(InputError):
            esi_msms.run_esi(inputs)


def test_runner_registers_the_page_id():
    assert list(esi_msms.RUNNERS) == ["ms2_esi"]


# ---------------------------------------------------------------- the first records of the multi-record exports
def test_massbank_and_sirius_records_are_split():
    one = "\n".join(["CH$NAME: a", "MS$FOCUSED_ION: PRECURSOR_M/Z 100.5", "PK$PEAK: m/z int. rel.int.", "  50.1 10 999", "//"])
    two = one.replace("a", "b").replace("50.1", "70.2").replace("100.5", "200.5")
    rec = first_record(one + "\n\n" + two)
    assert rec["n_records"] == 2 and rec["name"] == "a" and rec["precursor_mz"] == 100.5
    assert parse_spectrum(rec["peaks_text"]) == [(50.1, 10.0)]
    sirius = ">compound a\n>parentmass 100.5\n>ms2peaks\n50.1 10\n\n>compound b\n>parentmass 200.5\n>ms2peaks\n70.2 20"
    rec = first_record(sirius)
    assert rec["n_records"] == 2 and rec["precursor_mz"] == 100.5 and parse_spectrum(rec["peaks_text"]) == [(50.1, 10.0)]


# ---------------------------------------------------------------- real Docker (slow; skipped without Docker + image)
def _timed(label, fn, *args):
    t = time.perf_counter()
    res = fn(*args)
    print(f"\n[docker] {label}: {time.perf_counter() - t:.1f} s")
    return res


@real_docker
def test_real_caffeine_protonated_falls_back_to_ml():
    res = _timed("caffeine [M+H]+ (MSRB no, ML)", cf.run_cfmid, CAFFEINE, "[M+H]+")
    assert res["engine"] == "CFM-ID 4.4.10 (ML)" and res["adduct"] == "[M+H]+"
    assert [s["energy"] for s in res["spectra"]] == ["low", "medium", "high"] and all(s["peaks"] for s in res["spectra"])
    assert _has_peak(res, res["precursor_mz"]) and res["precursor_mz"] == pytest.approx(195.0877, abs=1e-4)
    assert any("C5H8N3+" in f for f in res["spectra"][0]["fragments"])


@real_docker
@pytest.mark.parametrize("adduct, mz", [("[M+H]+", 760.5851), ("[M+Na]+", 782.5670)])
def test_real_pc_goes_through_msrb(adduct, mz):
    res = _timed(f"PC 16:0/18:1 {adduct} (MSRB)", cf.run_cfmid, PC, adduct)
    assert res["engine"] == "MSRB 1.1.15" and res["adduct"] == adduct
    assert res["precursor_mz"] == pytest.approx(mz, abs=1e-4) and _has_peak(res, res["precursor_mz"])
    assert len(res["spectra"][2]["peaks"]) > len(res["spectra"][0]["peaks"])         # more fragments at 40 eV


@real_docker
def test_real_cardiolipin_dianion_matches_msrb_precursor():
    res = _timed("cardiolipin [M-2H]2- (MSRB)", cf.run_cfmid, CARDIOLIPIN, "[M-2H]2-")
    assert res["engine"] == "MSRB 1.1.15" and _has_peak(res, res["precursor_mz"])     # our m/z = MSRB's own 719.4475


@real_docker
def test_real_caffeine_sodium_is_not_covered():
    with pytest.raises(ValueError, match=r"do not cover this structure for \[M\+Na\]\+"):
        _timed("caffeine [M+Na]+ (error)", cf.run_cfmid, CAFFEINE, "[M+Na]+")
    with pytest.raises(InputError, match="do not cover"):
        esi_msms.run_esi({"structure": {"mode": "SMILES / InChI", "value": CAFFEINE}, "adduct": "[M+Na]+"})


@real_docker
def test_real_runner_end_to_end():
    out = _timed("runner: PC [M+H]+", esi_msms.run_esi, {"structure": {"mode": "SMILES / InChI", "value": PC}, "adduct": "[M+H]+"})
    assert out["outputs"]["summary"]["values"]["Engine"] == "MSRB 1.1.15"
    rec = first_record(out["outputs"]["formats"]["formats"]["MassBank"])
    assert rec["n_records"] == 3 and rec["precursor_mz"] == pytest.approx(760.58508)


def test_runner_writes_the_charge_of_a_multiply_charged_adduct(monkeypatch):
    out = _run(monkeypatch, adduct="[M-2H]2-")
    assert "CHARGE=2-" in out["outputs"]["formats"]["formats"]["MGF"]


def test_the_container_gets_no_network_no_capabilities_and_runs_as_the_caller(monkeypatch):
    monkeypatch.delenv("MZDESK_CFMID_MEMORY", raising=False)
    flags = cf.container_flags()
    assert flags[:8] == ["--network", "none", "--cap-drop", "ALL", "--security-opt", "no-new-privileges", "--pids-limit", "256"]
    if sys.platform.startswith("linux"):
        assert flags[-2:] == ["--user", f"{os.getuid()}:{os.getgid()}"]
    monkeypatch.setenv("MZDESK_CFMID_MEMORY", "1g")
    assert ["--memory", "1g"] == cf.container_flags()[8:10]
    seen = []
    monkeypatch.setattr(cf.subprocess, "run", lambda cmd, **k: seen.append(cmd) or subprocess.CompletedProcess(cmd, 0, "ok", ""))
    cf._docker_run(["true"], Path("/work"))
    assert seen[0][:6] == ["docker", "run", "--rm", "--name", seen[0][4], "--platform"] and "--network" in seen[0] and seen[0][-2:] == [cf.IMAGE, "true"]


def test_the_working_folder_can_be_set_for_a_service_with_a_private_tmp(monkeypatch, tmp_path):
    monkeypatch.delenv("MZDESK_CFMID_WORKDIR", raising=False)
    assert cf._workdir() is None
    target = tmp_path / "shared" / "cfmid"
    monkeypatch.setenv("MZDESK_CFMID_WORKDIR", str(target))
    assert cf._workdir() == str(target) and target.is_dir()
