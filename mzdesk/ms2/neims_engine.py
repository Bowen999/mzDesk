"""EI 质谱预测引擎（NEIMS），只依赖 numpy 和 RDKit。

模型：Wei, Belanger, Adams, Sculley, "Rapid Prediction of Electron-Ionization Mass Spectrometry Using
Neural Networks", ACS Cent. Sci. 2019, 5, 700-708（https://doi.org/10.1021/acscentsci.9b00085）。
原代码 github.com/brain-research/deep-molecular-massspec（Apache-2.0，已归档，TensorFlow 1.13）。
权重是原作者的官方 zip（从 Wayback Machine 存档下载），由 neims_setup.py 一次性转成 numpy；
这里按原代码的 MLP 结构用 numpy 重写前向计算，特征和后处理与原代码一致，并和原代码逐分子对照验证过（validation/neims/）。

输入一个中性分子的结构，输出 70 eV 的单位质量 EI 谱（m/z 0-999）。
网络：Morgan 计数指纹（半径 2，4096 位）→ 全连接 2000 → 7 个残差块（BN-ReLU-1000-BN-ReLU-2000）→ BN-ReLU，
再分出三个线性头：正向谱、以分子离子质量为锚点的反向谱、门控（sigmoid）；正向谱只保留 m/z <= M+5，
反向谱 rev[j] = bwd[M-j+5]，两者按门控加权后取 ReLU，最大峰归一化。
"""
import os
from functools import lru_cache

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

from mzdesk.ms2 import neims_setup

RDLogger.DisableLog("rdApp.*")

FP_RADIUS, FP_LENGTH = 2, 4096
N_BINS = 1000                   # 输出 m/z 0-999
MAX_MZ = N_BINS - 1
ABOVE_MASS = 5                  # 原代码 max_prediction_above_molecule_mass
N_BLOCKS = 7
BN_EPS = 1e-3
MIN_REL = 0.1                   # 显示阈值（相对强度 %）：原代码不设阈值，只是把强度按 999 取整写出；
                                # 0.1 % 约等于取整后不为 0 的最小值（1/999）
MAX_HEAVY_ATOMS = 100           # 原代码训练集的原子数上限（训练集中最大 91）
# 训练集（260k 个分子）里至少出现 100 次的元素（H 只作显式氢 / 氘出现）；其余元素训练样本极少或没有
COVERED_ELEMENTS = {"C", "H", "N", "O", "F", "S", "Cl", "Si", "Br", "P", "I", "B", "Se", "Co", "Fe", "Sn", "Ni", "As", "Ge", "Cr"}


class NEIMSError(RuntimeError):
    """权重缺失或损坏（不是用户输入的问题）。"""


# ---------------------------------------------------------------- weights
def weights_path():
    return os.environ.get("MZDESK_NEIMS_WEIGHTS") or str(neims_setup.default_weights())


@lru_cache(maxsize=2)
def _load(path):
    if not os.path.isfile(path):
        raise NEIMSError(f"NEIMS weights not found at {path}. Run: {neims_setup.SETUP_COMMAND}")
    try:
        with np.load(path) as z:
            w = {k: z[k] for k in z.files}
        layers = [f"dense{'' if i == 0 else f'_{i}'}" for i in range(18)]
        bns = [f"batch_normalization{'' if i == 0 else f'_{i}'}" for i in range(2 * N_BLOCKS + 1)]
        model = {"dense": [(w[f"{n}.kernel"], w[f"{n}.bias"]) for n in layers], "bn": []}
        for n in bns:                        # 推理时 BN 等价于 x * scale + shift
            scale = w[f"{n}.gamma"] / np.sqrt(w[f"{n}.moving_variance"] + np.float32(BN_EPS))
            model["bn"].append((scale.astype(np.float32), (w[f"{n}.beta"] - w[f"{n}.moving_mean"] * scale).astype(np.float32)))
        if model["dense"][0][0].shape != (FP_LENGTH, 2000) or model["dense"][17][0].shape != (2000, N_BINS):
            raise ValueError("unexpected tensor shapes")
    except (KeyError, ValueError, OSError) as exc:
        raise NEIMSError(f"NEIMS weights at {path} are unreadable ({exc!r}). Rebuild with: {neims_setup.SETUP_COMMAND} --force") from None
    return model


def _relu(x):
    return np.maximum(x, 0, out=x)


def _bn(x, layer):
    scale, shift = layer
    return x * scale + shift


def forward(fp, mass, model):
    """Raw network output (1000 values, before the max normalisation) for one fingerprint and nominal mass."""
    dense, bn = model["dense"], model["bn"]
    h = _relu(fp @ dense[0][0] + dense[0][1])
    for k in range(N_BLOCKS):
        a = _relu(_bn(h, bn[2 * k]))
        a = a @ dense[1 + 2 * k][0] + dense[1 + 2 * k][1]
        a = _relu(_bn(a, bn[2 * k + 1]))
        h = h + a @ dense[2 + 2 * k][0] + dense[2 + 2 * k][1]
    h = _relu(_bn(h, bn[2 * N_BLOCKS]))
    head = [h @ dense[i][0] + dense[i][1] for i in (15, 16, 17)]
    forward_pred, backward_pred, gate = head[0], head[1], 1.0 / (1.0 + np.exp(-head[2]))
    j = np.arange(N_BINS)
    forward_pred = np.where(j > mass + ABOVE_MASS, 0.0, forward_pred)
    src = mass - j + ABOVE_MASS
    reverse_pred = np.where((src >= 0) & (src < N_BINS), backward_pred[np.clip(src, 0, N_BINS - 1)], 0.0)
    return _relu(gate * forward_pred + (1.0 - gate) * reverse_pred).astype(np.float32)


# ---------------------------------------------------------------- features
_FP = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_LENGTH)


def fingerprint(mol):
    """Hashed Morgan count fingerprint (radius 2, 4096), float32 - the model input."""
    return _FP.GetCountFingerprintAsNumPy(mol).astype(np.float32)


def model_mass(mol):
    """The mass the network is anchored on: exact monoisotopic mass rounded as the original code does it
    (float32, round half to even). Differs from the nominal mass by 1 for some large H-rich or polybrominated molecules."""
    return int(np.round(np.float32(rdMolDescriptors.CalcExactMolWt(mol))))


def nominal_mass(mol):
    """Nominal mass of M+•: sum of the integer masses of the most abundant isotopes (implicit H included)."""
    table = Chem.GetPeriodicTable()
    return sum(a.GetIsotope() or table.GetMostCommonIsotope(a.GetAtomicNum()) for a in Chem.AddHs(mol).GetAtoms())


def spectrum_array(mol, model=None):
    """1000-bin prediction scaled so the largest peak is 999 - the array the original code writes to its SDF."""
    pred = forward(fingerprint(mol), model_mass(mol), model or _load(weights_path()))
    top = pred.max()
    return pred / top * np.float32(999.0) if top > 0 else pred


# ---------------------------------------------------------------- input checks
def _prepare(smiles):
    """SMILES -> (mol, warnings); ValueError for input that cannot be predicted."""
    smiles = (smiles or "").strip()
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError("Invalid or empty SMILES.")
    warnings = []
    if len(Chem.GetMolFrags(mol)) > 1:
        mol = rdMolStandardize.LargestFragmentChooser(preferOrganic=True).choose(mol)
        warnings.append(f"The structure has several components; the largest ({Chem.MolToSmiles(mol)}) is used, "
                        "as EI measures the neutral molecule.")
    charge = Chem.GetFormalCharge(mol)
    if charge:                                  # a salt's ion (benzoate, ammonium ...) -> the neutral acid / base EI would measure
        mol = rdMolStandardize.Uncharger().uncharge(mol)
        left = Chem.GetFormalCharge(mol)
        if left != charge:
            warnings.append(f"The input charge ({charge:+d}) was neutralized; the neutral molecule ({Chem.MolToSmiles(mol)}) is predicted.")
        if left:                                # e.g. a quaternary ammonium cannot be neutralized
            warnings.append(f"The structure has a net charge of {left:+d}. NEIMS predicts EI of neutral molecules, "
                            "so this prediction is unreliable.")
    if mol.GetNumAtoms() > MAX_HEAVY_ATOMS:
        warnings.append(f"{mol.GetNumAtoms()} heavy atoms is beyond the model's training range (up to {MAX_HEAVY_ATOMS}).")
    rare = sorted({a.GetSymbol() for a in mol.GetAtoms()} - COVERED_ELEMENTS)
    if rare:
        warnings.append(f"{', '.join(rare)} {'is' if len(rare) == 1 else 'are'} rare or absent in the NEIMS training data; "
                        "the prediction is unreliable.")
    return mol, warnings


# ---------------------------------------------------------------- predict
def predict(smiles):
    """EI spectrum of one structure.

    Returns {"smiles", "formula", "molecular_ion_mz" (nominal mass of M+•), "model_mass" (what the network is
    anchored on, see model_mass()), "peaks": [(mz, rel 0-100)] sorted by m/z (peaks below MIN_REL % are dropped),
    "warnings": [str]}. ValueError = unusable input; NEIMSError = no weights.
    """
    model = _load(weights_path())
    mol, warnings = _prepare(smiles)
    mass = nominal_mass(mol)
    if mass > MAX_MZ:
        raise ValueError(f"The molecular ion m/z {mass} is above the model's maximum (m/z {MAX_MZ}).")
    spec = spectrum_array(mol, model)
    if not spec.max() > 0:
        warnings.append("The model predicts no peaks for this structure.")
    rel = spec / np.float32(999.0) * 100.0
    peaks = [(int(j), float(rel[j])) for j in np.flatnonzero(rel >= MIN_REL)]
    return {"smiles": Chem.MolToSmiles(mol), "formula": rdMolDescriptors.CalcMolFormula(mol),
            "molecular_ion_mz": mass, "model_mass": model_mass(mol), "peaks": peaks, "warnings": warnings}
