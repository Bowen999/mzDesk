"""在 Docker 里用原代码（TF 1.13 + RDKit 2019.03）和原权重预测，结果写成 JSON。

Python 3.6。由 run_original.sh 调用，不要在宿主机直接运行。
输入：molecules.csv（name,smiles,note）和 extra_test_set_smiles.txt（每行一个 SMILES）。
输出：每个分子的指纹（稀疏）、原代码用的分子量、原代码输出的 1000 维谱（最大峰 = 999）、原代码写进 SDF 的峰文本。
另外对前 N 个训练集分子只算指纹，用来比较不同 RDKit 版本。
"""
import csv
import json
import sys

sys.path.insert(0, "/orig")      # 原仓库（只读挂载）

import numpy as np                # noqa: E402
import rdkit                      # noqa: E402
import tensorflow as tf           # noqa: E402
from rdkit import Chem            # noqa: E402

import feature_utils              # noqa: E402
import spectra_predictor          # noqa: E402


def sparse(fp):
    idx = np.nonzero(fp)[0]
    return {"idx": idx.tolist(), "cnt": [int(fp[i]) for i in idx]}


def main(out_path, weights_dir, n_fp_only):
    mols, rows = [], []
    with open("/work/molecules.csv") as f:
        for r in csv.DictReader(f):
            rows.append((r["name"], r["smiles"]))
    with open("/work/extra_test_set_smiles.txt") as f:
        rows += [("test_set_%03d" % i, s.strip()) for i, s in enumerate(f) if s.strip()]
    mols = [Chem.MolFromSmiles(s) for _, s in rows]
    assert all(m is not None for m in mols)

    predictor = spectra_predictor.NeimsSpectraPredictor(model_checkpoint_dir=weights_dir)
    fps, weights = predictor.get_inputs_for_model_from_mol_list(mols)
    spectra = predictor.make_spectra_prediction(fps, weights)
    # 原代码写进 SDF 的文本
    updated = spectra_predictor.update_mols_with_spectra(mols, spectra)
    out = []
    for (name, smi), fp, w, spec, mol in zip(rows, fps, weights, spectra, updated):
        out.append({"name": name, "smiles": smi, "fp": sparse(fp), "mol_weight": float(w),
                    "spectrum": spec.astype(float).tolist(),
                    "sdf_text": mol.GetProp(spectra_predictor.PREDICTED_SPECTRA_PROP_NAME)})

    # 只比指纹：前 n_fp_only 个训练集分子
    fp_only = []
    key = predictor._fingerprint_key
    with open("/work/fp_smiles.txt") as f:
        for line in f:
            s = line.strip()
            m = Chem.MolFromSmiles(s) if s else None
            if m is not None:
                fp_only.append({"smiles": s, "fp": sparse(feature_utils.make_circular_fingerprint(m, key))})
            if len(fp_only) >= n_fp_only:
                break

    versions = {"python": sys.version.split()[0], "rdkit": rdkit.__version__, "tensorflow": tf.__version__,
                "numpy": np.__version__}
    with open(out_path, "w") as f:
        json.dump({"versions": versions, "molecules": out, "fp_only": fp_only}, f)
    print("wrote", out_path, versions, len(out), "molecules,", len(fp_only), "fingerprint-only")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
