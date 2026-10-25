"""Compare the original MolGpKa (original_output.json, from run_original.sh) with the numpy engine in mzdesk on the same molecules.

Run:  python validation/molgpka/compare.py          (after run_original.sh and python -m mzdesk.platform_fit.pka_setup)
Out:  comparison.csv (one row per molecule), summary.json; prints the summary. Also writes tests/fixtures/molgpka_reference.json
      (a few molecules with the original code's per-site pKa; used by tests/test_separation.py).

Three levels, from strict to practical:
  all_atoms   original network vs numpy network with EVERY atom of the molecule as the site (acid and base weights):
              exercises the feature code (distance, site flag, ring sizes, hybridization ...) on all kinds of atoms
  same_graph  site matching (143 SMARTS) and pKa on AddHs(MolFromSmiles(smiles)), input atom order, per atom index
  pipeline    the whole mzdesk path (largest fragment, neutralization, sites, one prediction per heavy atom) vs the original
              predict() (Uncharger, SMILES round trip); sites are matched by canonical atom rank. Only compared when both
              ended up with the same neutral molecule; the rest are the intended differences (salts, betaines ...), listed.
"""
import json
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger

from mzdesk.platform_fit import pka_engine as E

RDLogger.DisableLog("rdApp.*")
HERE = Path(__file__).parent
FIXTURE = HERE.parent.parent / "tests" / "fixtures" / "molgpka_reference.json"
FIXTURE_NAMES = ["ibuprofen", "lidocaine", "glycine", "morphine", "metformin", "phosphoric acid dimethyl ester"]


def mine_pipeline(smiles):
    """-> (canonical SMILES of the molecule used, {(kind, rank): [pKa, ...]}) in the same keys as run_original.py."""
    prep = E.prepare(smiles)
    sites = E.predict_sites(prep.mol)
    mol_h = Chem.AddHs(prep.mol)
    ranks = list(Chem.CanonicalRankAtoms(mol_h, breakTies=False))
    out = {}
    for s in sites:
        out.setdefault((s.kind, ranks[s.atom]), []).append(s.pka)
    return Chem.MolToSmiles(prep.mol), out


def main():
    data = json.loads((HERE / "original_output.json").read_text())
    rows, fixture = [], {}
    worst = {"all_atoms": 0.0, "same_graph": 0.0, "pipeline": 0.0}
    n_all = n_all_atoms = n_sites_same = n_site_sets_equal = 0
    pipe_compared, pipe_equal_sites, pipe_differ = 0, 0, []
    for rec in data["molecules"]:
        smiles = rec["smiles"]
        mol_h = Chem.AddHs(Chem.MolFromSmiles(smiles))
        row = {"name": rec["name"], "smiles": smiles}
        # ---- all atoms as the site
        if "all_atoms" in rec:
            d = 0.0
            for kind in E.pka_setup.KINDS:
                aids = sorted(int(a) for a in rec["all_atoms"][kind])
                mine = E.predict_aids(mol_h, kind, aids)
                d = max([d] + [abs(mine[a] - rec["all_atoms"][kind][str(a)]) for a in aids])
                n_all_atoms += len(aids)
            row["all_atoms_max_abs_diff"] = f"{d:.2e}"
            worst["all_atoms"] = max(worst["all_atoms"], d)
        # ---- site matching + pKa on the same graph
        acid_h, base_atoms = E.find_site_atoms(mol_h)
        sets_equal = (sorted(int(a) for a in rec["same_graph"]["acid"]) == acid_h and
                      sorted(int(a) for a in rec["same_graph"]["base"]) == base_atoms)
        d = 0.0
        for kind, aids in (("acid", acid_h), ("base", base_atoms)):
            mine = E.predict_aids(mol_h, kind, aids)
            for a in aids:
                if str(a) in rec["same_graph"][kind]:
                    d = max(d, abs(mine[a] - rec["same_graph"][kind][str(a)]))
                    n_sites_same += 1
        n_all += 1
        n_site_sets_equal += sets_equal
        worst["same_graph"] = max(worst["same_graph"], d)
        row.update(site_sets_equal=sets_equal, same_graph_max_abs_diff=f"{d:.2e}", n_sites_original=len(rec["same_graph"]["acid"]) + len(rec["same_graph"]["base"]))
        # ---- whole pipeline
        if "pipeline" in rec:
            smi_m, mine = mine_pipeline(smiles)
            if smi_m == rec["pipeline_smiles"]:
                orig = {}
                for s in rec["pipeline"]:
                    orig.setdefault((s["kind"], s["rank"]), []).append(s["pka"])
                same_keys = set(orig) == set(mine)
                d = 0.0
                for k in set(orig) & set(mine):
                    d = max(d, max(abs(mine[k][0] - v) for v in orig[k]))
                pipe_compared += 1
                pipe_equal_sites += same_keys
                worst["pipeline"] = max(worst["pipeline"], d)
                row.update(pipeline="compared", pipeline_sites_equal=same_keys, pipeline_max_abs_diff=f"{d:.2e}")
            else:
                pipe_differ.append(f"{rec['name']}: original {rec['pipeline_smiles']} / mzdesk {smi_m}")
                row.update(pipeline="different molecule (intended)")
        rows.append(row)
        if rec["name"] in FIXTURE_NAMES:
            fixture[rec["name"]] = {"smiles": smiles, "same_graph": rec["same_graph"]}
    cols = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("name", "smiles"), k))
    import csv
    with open(HERE / "comparison.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    summary = {
        "versions_original": data["versions"], "molecules": n_all,
        "all_atoms": {"atoms_compared": n_all_atoms, "max_abs_diff_pka": worst["all_atoms"]},
        "same_graph": {"site_sets_equal": f"{n_site_sets_equal} / {n_all}", "sites_compared": n_sites_same, "max_abs_diff_pka": worst["same_graph"]},
        "pipeline": {"compared": pipe_compared, "site_sets_equal": pipe_equal_sites, "max_abs_diff_pka": worst["pipeline"],
                     "different_molecule": len(pipe_differ), "different_molecule_list": pipe_differ},
    }
    (HERE / "summary.json").write_text(json.dumps(summary, indent=1))
    FIXTURE.parent.mkdir(exist_ok=True)
    FIXTURE.write_text(json.dumps({"versions": data["versions"], "molecules": fixture}, indent=1))
    print(json.dumps({k: v for k, v in summary.items() if k != "versions_original"}, indent=1))


if __name__ == "__main__":
    sys.exit(main())
