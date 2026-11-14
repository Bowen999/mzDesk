"""把原代码（Docker，original_output.json）和 mzdesk 的 numpy 引擎在同一批分子上逐项对比。

运行：python validation/neims/compare.py          （先运行 run_original.sh，并已执行 python -m mzdesk.ms2.neims_setup）
输出：comparison.csv（每个分子一行）、summary.json（汇总），并打印汇总；同时从原代码的输出里取 5 个分子写成
tests/fixtures/neims_reference.json（单元测试用的参考值）。
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger, __version__ as rdkit_version

from mzdesk.ms2 import neims_engine as E

RDLogger.DisableLog("rdApp.*")
HERE = Path(__file__).parent
FIXTURE = HERE.parent.parent / "tests" / "fixtures" / "neims_reference.json"
FIXTURE_NAMES = ["caffeine", "pentachlorobenzene", "malathion", "testosterone", "tristearin"]


def dense(sp):
    a = np.zeros(E.FP_LENGTH, dtype=np.float32)
    a[sp["idx"]] = sp["cnt"]
    return a


def sdf_text(arr):
    """原代码 feature_utils.convert_spectrum_array_to_string 的写法：非零 bin，强度按 %d 截断。"""
    return "\n".join("%d %d" % (p, arr[p]) for p in np.nonzero(arr)[0].tolist())


def display(arr):
    """显示用峰表：相对强度 >= MIN_REL %，保留 1 位小数（与页面一致）。"""
    rel = arr / np.float32(999.0) * 100.0
    return [(int(j), round(float(rel[j]), 1)) for j in np.flatnonzero(rel >= E.MIN_REL)]


def main():
    data = json.loads((HERE / "original_output.json").read_text())
    model = E._load(E.weights_path())
    rows = []
    for rec in data["molecules"]:
        mol = Chem.MolFromSmiles(rec["smiles"])
        fp_o, w_o = dense(rec["fp"]), rec["mol_weight"]
        fp_e = E.fingerprint(mol)
        mass_o = int(np.round(np.float32(w_o)))
        mass_e = E.model_mass(mol)
        orig = np.asarray(rec["spectrum"], dtype=np.float32)
        eng = E.spectrum_array(mol, model)
        # 把原代码的指纹和质量喂给 numpy 前向：只剩数值差异（排除 RDKit 版本的影响）
        eng_same_input = E.forward(fp_o, mass_o, model)
        top = eng_same_input.max()
        eng_same_input = eng_same_input / top * np.float32(999.0) if top > 0 else eng_same_input
        text_o, text_e = rec["sdf_text"], sdf_text(eng)
        lo, le = text_o.split("\n") if text_o else [], text_e.split("\n") if text_e else []
        po, pe = display(orig), display(eng)
        do, de = dict(po), dict(pe)
        common = set(do) & set(de)
        rows.append({
            "name": rec["name"], "smiles": rec["smiles"], "mass": mass_o,
            "fp_equal": bool(np.array_equal(fp_o, fp_e)), "mass_equal": mass_o == mass_e,
            "max_abs_diff": float(np.abs(orig - eng).max()),                 # 999 刻度
            "max_abs_diff_same_input": float(np.abs(orig - eng_same_input).max()),
            "sdf_text_identical": text_o == text_e,
            "n_lines_orig": len(lo), "n_lines_engine": len(le),
            "n_int_intensity_differ": sum(a != b for a, b in zip(lo, le)) + abs(len(lo) - len(le)),
            "display_peaks_orig": len(po), "display_peaks_identical": po == pe,
            "display_peaks_only_in_one": len(set(do) ^ set(de)),
            "display_max_rel_diff": max((abs(do[k] - de[k]) for k in common), default=0.0),
        })
    # 指纹：不同 RDKit 版本
    fp_only = data["fp_only"]
    diff_fp = []
    for rec in fp_only:
        mol = Chem.MolFromSmiles(rec["smiles"])
        if mol is None:
            diff_fp.append((rec["smiles"], "host RDKit cannot parse"))
        elif not np.array_equal(dense(rec["fp"]), E.fingerprint(mol)):
            diff_fp.append((rec["smiles"], "different counts"))

    named = [r for r in rows if not r["name"].startswith("test_set_")]
    extra = [r for r in rows if r["name"].startswith("test_set_")]

    def stats(group):
        same = [r for r in group if r["fp_equal"]]         # 指纹一致的分子：只剩浮点数值差异
        return {
            "molecules": len(group),
            "fingerprints_equal": len(same),
            "model_mass_equal": sum(r["mass_equal"] for r in group),
            "max_abs_diff_999_scale_fp_equal": max(r["max_abs_diff"] for r in same),
            "max_abs_diff_999_scale_all": max(r["max_abs_diff"] for r in group),
            "max_abs_diff_999_scale_with_original_fp": max(r["max_abs_diff_same_input"] for r in group),
            "sdf_peak_text_identical": sum(r["sdf_text_identical"] for r in group),
            "display_peak_lists_identical": sum(r["display_peaks_identical"] for r in group),
            "display_peak_lists_identical_fp_equal": sum(r["display_peaks_identical"] for r in same),
            "display_peaks_only_in_one_list_fp_equal": sum(r["display_peaks_only_in_one"] for r in same),
            "display_max_rel_diff_percent_points_fp_equal": round(max(r["display_max_rel_diff"] for r in same), 3),
        }

    summary = {"original": data["versions"], "host": {"rdkit": rdkit_version, "numpy": np.__version__},
               "named_molecules": stats(named), "extra_test_set_molecules": stats(extra),
               "fingerprint_only": {"molecules": len(fp_only), "different": len(diff_fp), "smiles": [s for s, _ in diff_fp]},
               "min_rel_percent": E.MIN_REL}
    with open(HERE / "comparison.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    (HERE / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    for r in rows:
        if not (r["fp_equal"] and r["sdf_text_identical"] and r["display_peaks_identical"]):
            print("DIFF", r["name"], r["smiles"][:50], {k: r[k] for k in ("fp_equal", "max_abs_diff", "n_int_intensity_differ", "display_peaks_only_in_one", "display_max_rel_diff")})
    write_fixture(data)


def write_fixture(data):
    """Reference values from the ORIGINAL run (fingerprint, nominal mass, 999-scale spectrum) for a few molecules."""
    by_name = {r["name"]: r for r in data["molecules"]}
    out = {"source": "original deep-molecular-massspec f82884a + official weights, run in Docker " + json.dumps(data["versions"]),
           "scale": "spectrum values are the original's output with the largest peak = 999", "molecules": []}
    for name in FIXTURE_NAMES:
        r = by_name[name]
        spec = np.asarray(r["spectrum"])
        keep = [j for j in np.flatnonzero(spec) if round(spec[j], 3) > 0]
        out["molecules"].append({"name": name, "smiles": r["smiles"], "model_mass": int(np.round(np.float32(r["mol_weight"]))),
                                 "fp": r["fp"], "mz": [int(j) for j in keep], "value": [round(float(spec[j]), 3) for j in keep]})
    FIXTURE.write_text(json.dumps(out, separators=(",", ":")) + "\n")
    print("wrote", FIXTURE, FIXTURE.stat().st_size, "bytes")


if __name__ == "__main__":
    sys.exit(main())
