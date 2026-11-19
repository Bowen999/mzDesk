"""Run the ORIGINAL MolGpKa code (Xundrug/MolGpKa, torch + torch_geometric) on molecules.csv and extra_test_set_smiles.txt.

Run it with a Python that has torch, torch_geometric, pandas and RDKit (run_original.sh builds a throw-away venv for that):

    python run_original.py <MolGpKa repo> <output.json>

The original code is used unmodified: its src/ directory is put on sys.path and its own predict_pka / descriptor / net /
ionization_group modules do the work. Only two things are adapted around it, neither changes any computation:
  * torch_scatter is not installable everywhere; the original gcn_conv.py only uses scatter_add, so a 3-line stand-in
    (index_add_) is registered as the module torch_scatter before the import.
  * ionization_group.py finds smarts_pattern.tsv relative to the current directory; its module variable is pointed at the repo file.
Also, load_model() is cached so that the 300+ molecules do not reload the weights each time (same weights, same eval mode).

For every molecule it records
  pipeline   what the original predict() returns for the SMILES (uncharge, SMILES round trip, AddHs, all sites), with every
             site identified by canonical atom ranks (the round trip renumbers the atoms); pipeline_smiles = the molecule it used
  same_graph the original network and site matching on AddHs(MolFromSmiles(smiles)) with the input atom order, keyed by atom index
  all_atoms  (first 40 molecules of molecules.csv, up to 90 atoms) the original network with EVERY atom as the site, both weight sets
"""
import csv
import json
import os
import sys
import time
import types
from functools import lru_cache
from pathlib import Path

import torch


def _install_torch_scatter_standin():
    mod = types.ModuleType("torch_scatter")

    def scatter_add(src, index, dim=0, dim_size=None):
        size = list(src.shape)
        size[dim] = int(index.max()) + 1 if dim_size is None else dim_size
        return torch.zeros(size, dtype=src.dtype, device=src.device).index_add_(dim, index, src)

    mod.scatter_add = scatter_add
    sys.modules["torch_scatter"] = mod


def main(repo, out_path):
    here = Path(__file__).parent
    repo = Path(repo).resolve()
    _install_torch_scatter_standin()
    sys.path.insert(0, str(repo / "src"))
    import warnings
    warnings.filterwarnings("ignore")
    from rdkit import Chem, RDLogger, __version__ as rdkit_version
    RDLogger.DisableLog("rdApp.*")
    import torch_geometric
    import predict_pka as orig                                  # the original module
    from utils import ionization_group as ig
    from utils.descriptor import mol2vec

    ig.smarts_file = str(repo / "src" / "utils" / "smarts_pattern.tsv")
    orig.load_model = lru_cache(maxsize=None)(orig.load_model)
    models = {"acid": orig.load_model(str(repo / "models" / "weight_acid.pth")),
              "base": orig.load_model(str(repo / "models" / "weight_base.pth"))}
    # predict_acid / predict_base call load_model with their own path string, so they load (and then cache) identical copies

    with open(here / "molecules.csv", newline="") as f:
        named = [(r["name"], r["smiles"]) for r in csv.DictReader(f)]
    extra = [("extra_%03d" % i, s) for i, s in enumerate((here / "extra_test_set_smiles.txt").read_text().split())]

    def site_records(mol_h, base_dict, acid_dict):
        ranks = list(Chem.CanonicalRankAtoms(mol_h, breakTies=False))
        recs = []
        for aid, pka in acid_dict.items():
            heavy = mol_h.GetAtomWithIdx(aid).GetNeighbors()[0].GetIdx()
            recs.append({"kind": "acid", "rank": ranks[heavy], "element": mol_h.GetAtomWithIdx(heavy).GetSymbol(), "pka": float(pka)})
        for aid, pka in base_dict.items():
            recs.append({"kind": "base", "rank": ranks[aid], "element": mol_h.GetAtomWithIdx(aid).GetSymbol(), "pka": float(pka)})
        return recs

    out = {"versions": {"torch": torch.__version__, "torch_geometric": torch_geometric.__version__, "rdkit": rdkit_version,
                        "python": sys.version.split()[0], "repo_commit": os.environ.get("MOLGPKA_COMMIT", "")},
           "molecules": []}
    t0 = time.time()
    for n, (name, smiles) in enumerate(named + extra):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        rec = {"name": name, "smiles": smiles}
        try:
            base_dict, acid_dict, mol_h = orig.predict_for_protonate(mol)
            rec["pipeline"] = site_records(mol_h, base_dict, acid_dict)
            rec["pipeline_smiles"] = Chem.MolToSmiles(Chem.RemoveHs(mol_h))
        except Exception as exc:  # noqa: BLE001 — record, the comparison reports it
            rec["pipeline_error"] = repr(exc)
        mol_h0 = Chem.AddHs(mol)
        same = {}
        for kind in ("acid", "base"):
            aids = ig.get_ionization_aid(mol_h0, acid_or_base=kind)
            same[kind] = {}
            for aid in sorted(set(aids)):
                with torch.no_grad():
                    same[kind][aid] = float(models[kind](mol2vec(mol_h0, aid))[0][0])
        rec["same_graph"] = same
        if n < 40 and mol_h0.GetNumAtoms() <= 90:
            rec["all_atoms"] = {}
            for kind in ("acid", "base"):
                rec["all_atoms"][kind] = {}
                for aid in range(mol_h0.GetNumAtoms()):
                    with torch.no_grad():
                        rec["all_atoms"][kind][aid] = float(models[kind](mol2vec(mol_h0, aid))[0][0])
        out["molecules"].append(rec)
    out["seconds"] = round(time.time() - t0, 1)
    Path(out_path).write_text(json.dumps(out))
    print(f"{len(out['molecules'])} molecules in {out['seconds']} s -> {out_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
