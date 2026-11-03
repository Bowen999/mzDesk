"""ESI MS/MS 预测引擎：在 Docker 里跑 wishartlab/cfmid 镜像（CFM-ID 4.4.10 + MSRB 1.1.15）。

先试 MSRB（脂质规则碎裂，约 20 类脂质、7 种加合物）；MSRB 不覆盖该结构的这个加合物时，回退 cfm-predict
（CFM-ID 4 的 ML 模型，只有 [M+H]+ / [M-H]-）。永远不返回别的加合物。三个碰撞能量 = low / medium / high = 10 / 20 / 40 eV。
本模块顶层只用标准库（spec.py 要读 ADDUCTS）；rdkit 等在函数里才导入。
"""
import copy
import os
import re
import shutil
import subprocess
import sys
import tempfile
import uuid
from functools import lru_cache
from pathlib import Path

IMAGE = "wishartlab/cfmid:msml4.4.10_msrb1.1.15"
PLATFORM = "linux/amd64"
PULL_CMD = f"docker pull --platform {PLATFORM} {IMAGE}"
TIMEOUT_S = 900                     # per docker call; emulated cfm-predict on a big lipid takes minutes
MSRB_JAR = "/cfmid/msrb-fragmenter.jar"
MODEL_DIR = "/trained_models_cfmid4.0"

# 页面给用户选的加合物（[M+H]+ / [M-H]- 在前）：MSRB 的 7 种；ML 模型只覆盖前两种
ADDUCTS = ["[M+H]+", "[M-H]-", "[M+Na]+", "[M+NH4]+", "[M+Li]+", "[M+Cl]-", "[M-2H]2-"]
ML_ADDUCTS = ("[M+H]+", "[M-H]-")
_MSRB_NAME = {"[M-2H]2-": "[M-2H](2H)-"}            # MSRB 的写法
# energy0 / 1 / 2：MSRB jar 的用法说明写明 10 / 20 / 40 eV，CFM-ID 4 论文同
ENERGIES = (("low", 10), ("medium", 20), ("high", 40))

DOCKER_DOWN = "Docker is not running. Start Docker Desktop and try again."


class CFMIDError(RuntimeError):
    """Docker / CFM-ID / MSRB 运行失败；log 是工具输出的末尾。"""

    def __init__(self, message, log=""):
        super().__init__(message)
        self.log = log


# ---------------------------------------------------------------- docker
def check_docker():
    """Raise CFMIDError unless Docker answers and the image is already local (never pulls)."""
    if shutil.which("docker") is None:
        raise CFMIDError(DOCKER_DOWN)
    try:
        proc = subprocess.run(["docker", "image", "inspect", IMAGE], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired):
        raise CFMIDError(DOCKER_DOWN) from None
    if proc.returncode != 0:
        err = (proc.stderr or "").strip()
        if "No such image" in err:
            raise CFMIDError(f"The CFM-ID image is not installed. Run: {PULL_CMD}", err)
        raise CFMIDError(DOCKER_DOWN, err)


def _workdir():
    """Where the run's files go: MZDESK_CFMID_WORKDIR when set, else the system temp. The folder must be one the Docker daemon sees too: a service
    with a private /tmp (the server's systemd unit) has to set it."""
    folder = os.environ.get("MZDESK_CFMID_WORKDIR")
    if folder:
        Path(folder).mkdir(parents=True, exist_ok=True)
    return folder or None


def container_flags():
    """What the container may do: no network, no capabilities, no privilege gain, a cap on processes (and on memory, MZDESK_CFMID_MEMORY, e.g. 1g).
    On Linux it also runs as the caller's user, so the files it writes into the mounted folder can be read and removed by the caller."""
    flags = ["--network", "none", "--cap-drop", "ALL", "--security-opt", "no-new-privileges", "--pids-limit", "256"]
    memory = (os.environ.get("MZDESK_CFMID_MEMORY") or "").strip()
    if memory:
        flags += ["--memory", memory]
    if sys.platform.startswith("linux"):
        flags += ["--user", f"{os.getuid()}:{os.getgid()}"]
    return flags


def _docker_run(args, workdir):
    """Run ``args`` in the image with ``workdir`` mounted at /data (the cwd) -> (stdout + stderr, exit code)."""
    name = f"mzdesk-cfmid-{uuid.uuid4().hex[:10]}"
    cmd = ["docker", "run", "--rm", "--name", name, "--platform", PLATFORM, "--pull", "never", *container_flags(),
           "--mount", f"type=bind,source={workdir},target=/data", "-w", "/data", IMAGE, *args]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT_S)
    except subprocess.TimeoutExpired:
        subprocess.run(["docker", "rm", "-f", name], capture_output=True)       # killing the client leaves the container running
        raise CFMIDError(f"CFM-ID did not finish within {TIMEOUT_S // 60} minutes. Try a smaller molecule.") from None
    except OSError:
        raise CFMIDError(DOCKER_DOWN) from None
    return ((proc.stdout or "") + (proc.stderr or "")).strip(), proc.returncode


# ---------------------------------------------------------------- structure and ion
def prepare(smiles):
    """SMILES -> (canonical neutral SMILES, warnings). ValueError for input the engine cannot use.

    A salt / mixture keeps its largest part; a charged part is neutralised; a charge that stays
    (quaternary ammonium cation, ...) is an error, because the adduct is added by the engine.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors
    from rdkit.Chem.MolStandardize import rdMolStandardize

    RDLogger.DisableLog("rdApp.*")
    smiles = (smiles or "").strip()
    if not smiles:
        raise ValueError("Enter the SMILES of the structure.")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string. Please check the structure.")
    warnings = []
    parts = Chem.GetMolFrags(mol, asMols=True)
    if len(parts) > 1:
        parts = sorted(parts, key=lambda f: (f.GetNumHeavyAtoms(), Descriptors.ExactMolWt(f)), reverse=True)
        mol = parts[0]
        warnings.append(f"The SMILES has {len(parts)} disconnected parts; only the largest is predicted. "
                        f"Ignored: {', '.join(Chem.MolToSmiles(p) for p in parts[1:])}.")
    charge = Chem.GetFormalCharge(mol)
    mol = rdMolStandardize.Uncharger().uncharge(mol)
    if Chem.GetFormalCharge(mol) != 0:
        raise ValueError(f"The structure carries a net charge ({Chem.GetFormalCharge(mol):+d}). "
                         "Give the neutral molecule; the adduct is added for you.")
    if charge:
        warnings.append(f"The input charge ({charge:+d}) was neutralized; the neutral molecule is predicted.")
    return Chem.MolToSmiles(mol), warnings


def ion_info(mol, adduct):
    """(precursor m/z, ion formula with charge) of the neutral ``mol`` as ``adduct``; ValueError if it cannot form."""
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    from mzdesk.ms1 import msbuddy_engine as bd          # the MS1 theory page's calculation

    formula = rdMolDescriptors.CalcMolFormula(mol)
    if adduct == "[M+Li]+":                              # msbuddy has no Li: add Li+ to the [M]+ ion
        from mzdesk.ms2 import metfrag_engine as mf
        mz, ion = bd.ion_mz_and_formula(formula, "[M]+")
        if mz is not None:
            counts = mf._formula_counts(ion[:-1])
            counts["Li"] = 1
            mz += Chem.GetPeriodicTable().GetMostCommonIsotopeMass("Li")
            ion = mf.formula_str(counts) + "+"
    else:
        mz, ion = bd.ion_mz_and_formula(formula, adduct)
        charge = adduct[adduct.rindex("]") + 1:]          # the ion formula carries the sign only; keep a "2-"
        if ion and len(charge) > 1:
            ion = f"[{ion[:-1]}]{charge}"
    if mz is None:
        raise ValueError(f"{formula} cannot form {adduct}: the formula has elements other than "
                         "C H N O P S F Cl Br I Na K, or the adduct removes atoms it does not have.")
    return mz, ion


# ---------------------------------------------------------------- output parsing
_HEADER = re.compile(r"^#In-silico ESI-MS/MS (\S+)\s+Spectra")
_BY = re.compile(r"^#PREDICTED BY (.+?)\s*$")
_ENERGY = re.compile(r"^energy(\d+)\s*$")
_PEAK = re.compile(r"^(\d+\.\d+)\s+(\d+(?:\.\d+)?)\s+(\d[\d ]*?)\s*\(")     # mz, intensity, fragment ids, (scores)
_FRAGMENT = re.compile(r"^(\d+)\s+\d+\.\d+\s+(\S+)\s*$")                      # id, mass, SMILES


def parse_log(text):
    """Parse an MSRB or cfm-predict log -> {"adduct", "by", "energies": [[(mz, rel, [fragment SMILES])]]}.

    Both tools write the same layout: ``#`` header, ``energyN`` blocks of ``mz intensity ids (scores)``,
    then a fragment table ``id mass SMILES``. Intensities are scaled so each energy's base peak is 100.
    """
    adduct = by = None
    energies, fragments = [], {}
    for line in (text or "").splitlines():
        m = _HEADER.match(line)
        if m:
            adduct = m.group(1)
        m = _BY.match(line)
        if m:
            by = m.group(1)
        if _ENERGY.match(line):
            energies.append([])
            continue
        m = _PEAK.match(line)
        if m and energies:
            energies[-1].append((float(m.group(1)), float(m.group(2)), [int(i) for i in m.group(3).split()]))
            continue
        m = _FRAGMENT.match(line)
        if m:
            fragments[int(m.group(1))] = m.group(2)
    out = []
    for peaks in energies:
        top = max((p[1] for p in peaks), default=0.0) or 1.0
        out.append([(mz, y / top * 100.0, [fragments[i] for i in ids if i in fragments]) for mz, y, ids in sorted(peaks)])
    return {"adduct": adduct, "by": by, "energies": out}


@lru_cache(maxsize=2048)
def fragment_text(smiles):
    """'C5H9N2+ · SMILES' for a fragment (formula from RDKit; just the SMILES when RDKit cannot read it)."""
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    mol = Chem.MolFromSmiles(smiles)
    return f"{rdMolDescriptors.CalcMolFormula(mol)} · {smiles}" if mol is not None else smiles


def _annotation(fragments):
    if not fragments:
        return ""
    more = f" (+{len(fragments) - 1} more)" if len(fragments) > 1 else ""
    return fragment_text(fragments[0]) + more


def _spectra(parsed):
    if len(parsed["energies"]) != len(ENERGIES):
        raise CFMIDError(f"The tool returned {len(parsed['energies'])} spectra instead of {len(ENERGIES)}.")
    return [{"energy": label, "ev": ev, "peaks": [(mz, y) for mz, y, _ in peaks],
             "fragments": [_annotation(f) for _, _, f in peaks]}
            for (label, ev), peaks in zip(ENERGIES, parsed["energies"])]


# ---------------------------------------------------------------- the two tools
def _msrb(tmp, smiles, adduct):
    """(MSRB log text, tool output); the text is None when MSRB does not cover ``adduct`` (-na: never another adduct)."""
    name = _MSRB_NAME.get(adduct, adduct)
    out, _ = _docker_run(["java", "-jar", MSRB_JAR, "-a", name, "-na", "-ismi", smiles, "-o", "msrb"], tmp)
    for path in tmp.glob("msrb_*"):
        text = path.read_text(errors="replace")
        if parse_log(text)["adduct"] == name:
            return text, out
    return None, out


def _ml(tmp, smiles, adduct):
    model = f"{MODEL_DIR}/{adduct}"
    out, _ = _docker_run(["cfm-predict", smiles, "0.001", f"{model}/param_output.log", f"{model}/param_config.txt",
                          "1", "/data/cfm.log", "1"], tmp)
    path = tmp / "cfm.log"
    if not path.exists() or not path.read_text().strip():
        raise CFMIDError("CFM-ID returned no spectrum for this structure.", out[-4000:])
    return path.read_text(errors="replace"), out


@lru_cache(maxsize=32)
def _predict(smiles, adduct):
    from rdkit import Chem

    precursor_mz, ion_formula = ion_info(Chem.MolFromSmiles(smiles), adduct)      # user errors before Docker starts
    check_docker()
    with tempfile.TemporaryDirectory(prefix="mzdesk-cfmid-", dir=_workdir()) as d:
        tmp = Path(d).resolve()
        text, log = _msrb(tmp, smiles, adduct)
        engine = "MSRB 1.1.15"                                                   # overwritten by the log's own header
        if text is None:
            if adduct not in ML_ADDUCTS:
                raise ValueError(f"MSRB rules do not cover this structure for {adduct}; "
                                 "the CFM-ID ML model supports only [M+H]+ and [M-H]-.")
            text, ml_log = _ml(tmp, smiles, adduct)
            log, engine = f"{log}\n{ml_log}", "CFM-ID 4.4.10 (ML)"
    parsed = parse_log(text)
    if parsed["by"]:
        engine = parsed["by"] + (" (ML)" if parsed["by"].startswith("CFM-ID") else "")
    return {"engine": engine, "smiles": smiles, "adduct": adduct, "precursor_mz": precursor_mz,
            "ion_formula": ion_formula, "spectra": _spectra(parsed), "log": log[-4000:]}


def run_cfmid(smiles, adduct):
    """Predict the ESI-MS/MS spectrum of ``smiles`` as ``adduct`` at low / medium / high energy.

    Returns {"engine", "smiles" (neutral, canonical), "adduct", "precursor_mz", "ion_formula",
    "spectra": [{"energy", "ev", "peaks": [(mz, rel 0-100)], "fragments": [str per peak]}],
    "warnings", "log"}. ValueError = input problem (invalid SMILES, charged, adduct not covered);
    CFMIDError = Docker or the tool failed. Results are cached per (canonical SMILES, adduct).
    """
    if adduct not in ADDUCTS:
        raise ValueError(f"Unsupported adduct {adduct!r}. Choose one of {', '.join(ADDUCTS)}.")
    canonical, warnings = prepare(smiles)
    result = copy.deepcopy(_predict(canonical, adduct))
    result["warnings"] = warnings
    return result
