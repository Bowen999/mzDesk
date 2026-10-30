#!/usr/bin/env python3
"""JESTR: score every SIRIUS PubChem candidate of every CASMI 2022 challenge.

Same model path as JESTR/smoke_test/run_jestr_smoke.py (NPLIB1 weights, 10-element
atom featurizer = 74-dim, intensity scaled to 999, m/z <= 1000, 1 Da bins,
log10over3), cosine(mol embedding, spectrum embedding). The spectrum embedding
is computed once per challenge; candidates are featurised in blocks to bound memory.

Candidates JESTR cannot featurise (disconnected SMILES, single atom, elements
outside H C O N P S Cl F Br I) get no score and therefore no rank.

Output per challenge (atomic write; existing files are skipped = resumable):
  <out>/scores/challenge_NNN.tsv.gz   ikey14, formula, adduct, score
  <out>/meta/challenge_NNN.tsv        n_candidates, n_scored, seconds
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("DGLBACKEND", "pytorch")
import torch.utils.data.datapipes.utils.common as _td_common  # noqa: E402

if not hasattr(_td_common, "DILL_AVAILABLE"):
    _td_common.DILL_AVAILABLE = False  # torchdata 0.7.1 vs torch>=2.3

import dgl  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402
import yaml  # noqa: E402
from rdkit import Chem, RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from jm_common import assign_shards, atomic_write_tsv, report  # noqa: E402

ELEMENTS_74 = ["H", "C", "O", "N", "P", "S", "Cl", "F", "Br", "I"]


def norm_mzi(mz, inten):
    """JESTR training-time normalisation: intensity to 999, keep m/z <= 1000."""
    mz = np.asarray(mz, dtype=float).reshape(-1, 1)
    inten = np.asarray(inten, dtype=float)
    inten = (inten / inten.max() * 999).reshape(-1, 1)
    mzi = np.hstack((inten, mz))
    return mzi[mzi[:, 1] <= 1000]


class Jestr:
    def __init__(self, repo: Path, device=torch.device("cpu")):
        sys.path.insert(0, str(repo))
        from dataset import mol_to_graph
        from train_contr import train_contr
        from utils import DatasetBuilder, collate_spectra_data

        with open(repo / "params.yaml") as fh:
            params = yaml.safe_load(fh)
        params["load_dicts"] = False
        # Released NPLIB1 weights are 74-dim: 'full' atom featurizer with 10 elements
        # (params.yaml lists 14, which gives 78 dims and fails to load).
        params["element_list"] = ELEMENTS_74
        params["pretrained_mol_enc_model"] = str(repo / "data/NPLIB1/pretrained_mol_enc_model_1707829192911_best.pt")
        params["pretrained_spec_enc_model"] = str(repo / "data/NPLIB1/pretrained_spec_enc_model_1707829192911_best.pt")
        self.params, self.device = params, device
        self.mol_to_graph = mol_to_graph

        db = DatasetBuilder(params["exp"])
        for attr in ("data_dict", "mol_dict", "pair_dict", "split_dict", "fp_dict",
                     "in_to_id_dict", "in_to_id_dict_wneg"):
            setattr(db, attr, [])
        db.data_dir = str(repo / "data/NPLIB1/")
        self.graphs = {}
        db.molgraph_dict = self.graphs
        # train_contr infers the node feature size from one graph in the dict
        assert mol_to_graph(Chem.MolFromSmiles("CCO"), "_probe", self.graphs, params, device)
        self.probe = self.graphs["_probe"]
        with open(os.devnull, "w") as quiet:
            self.mol_enc, self.spec_enc, _ = train_contr(db, self.graphs, params, quiet, device, db.data_dir, True)
        self.mol_enc.eval()
        self.spec_enc.eval()
        self.collate = collate_spectra_data(self.graphs, params)
        self.graphs.clear()

    @torch.no_grad()
    def spec_embedding(self, mzi: np.ndarray) -> torch.Tensor:
        """Spectrum branch only; the probe graph just satisfies the collate function."""
        self.graphs["_probe"] = self.probe
        _, mz_b, int_b, pad, _, _, lengths, _ = self.collate([("_probe", mzi, [0.0] * 4096, 1)])
        self.graphs.clear()
        return self.spec_enc(mz_b, int_b, pad, lengths)

    @torch.no_grad()
    def score(self, spec_e: torch.Tensor, smiles: list[str], batch_size: int, block: int) -> np.ndarray:
        scores = np.full(len(smiles), np.nan, dtype=np.float32)
        for b0 in range(0, len(smiles), block):
            self.graphs.clear()
            keys, idx = [], []
            for i in range(b0, min(b0 + block, len(smiles))):
                mol = Chem.MolFromSmiles(smiles[i])
                if mol is None:
                    continue
                key = str(i)
                try:
                    ok = self.mol_to_graph(mol, key, self.graphs, self.params, self.device)
                except Exception:
                    ok = False
                if ok:
                    keys.append(key)
                    idx.append(i)
            for s in range(0, len(keys), batch_size):
                g = dgl.batch([self.graphs[k] for k in keys[s:s + batch_size]])
                mol_e = self.mol_enc(g, g.ndata["h"])
                cos = torch.nn.functional.cosine_similarity(mol_e, spec_e.expand_as(mol_e), dim=-1)
                scores[idx[s:s + batch_size]] = cos.numpy()
        self.graphs.clear()
        return scores


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--repo", required=True, help="JESTR repo dir")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--batch-size", type=int, default=256)
    ap.add_argument("--block", type=int, default=4096, help="candidates featurised at once")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--progress-file", default=None)
    args = ap.parse_args()

    torch.set_num_threads(args.threads)
    inp, out = Path(args.inputs), Path(args.out)
    (out / "scores").mkdir(parents=True, exist_ok=True)
    (out / "meta").mkdir(parents=True, exist_ok=True)

    chal = pd.read_csv(inp / "challenges.tsv", sep="\t")
    chal = chal[(chal.n_candidates > 0) & (chal.n_peaks_merged > 0)]
    mine = assign_shards(chal, args.nshards)[args.shard]
    units_total = int(mine.n_candidates.sum())
    units_done, items_done = 0, 0

    todo = []
    for r in mine.itertuples(index=False):
        if (out / "scores" / f"{r.challenge}.tsv.gz").exists():
            units_done += r.n_candidates
            items_done += 1
        else:
            todo.append(r)
    report(args.progress_file, units_done, units_total, items_done, len(mine))
    if not todo:
        return 0

    model = Jestr(Path(args.repo))
    for r in todo:
        t0 = time.time()
        spec = pd.read_csv(inp / "spectra" / f"{r.challenge}.tsv", sep="\t")
        mzi = norm_mzi(spec["mz"].values, spec["intensity"].values)
        cand = pd.read_csv(inp / "candidates" / f"{r.challenge}.tsv.gz", sep="\t", dtype=str,
                           keep_default_na=False)
        if len(mzi):
            spec_e = model.spec_embedding(mzi)
            cand["score"] = model.score(spec_e, cand["smiles"].tolist(), args.batch_size, args.block)
        else:
            cand["score"] = np.nan
        scored = cand.dropna(subset=["score"])
        atomic_write_tsv(scored[["ikey14", "formula", "adduct", "score"]], out / "scores" / f"{r.challenge}.tsv.gz")
        atomic_write_tsv(pd.DataFrame([dict(challenge=r.challenge, n_candidates=len(cand), n_scored=len(scored),
                                            n_peaks_used=len(mzi), seconds=round(time.time() - t0, 2))]),
                         out / "meta" / f"{r.challenge}.tsv")
        units_done += r.n_candidates
        items_done += 1
        report(args.progress_file, units_done, units_total, items_done, len(mine))
        print(f"{r.challenge}: {len(scored)}/{len(cand)} scored, {time.time() - t0:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
