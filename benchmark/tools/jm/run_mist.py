#!/usr/bin/env python3
"""MIST: score SIRIUS PubChem candidates of the positive-mode CASMI 2022 challenges.

For every (challenge, formula, adduct) hypothesis that SIRIUS proposed (mist_hypotheses.tsv)
the official MIST CLIs are run, in chunks of challenges:
  1. subformulae/assign_subformulae.py   (20 ppm, max 50 formulae = MIST training settings)
  2. pred_fp.py                          (mist_fp_canopus_pretrain.ckpt  -> 4096-bit Morgan r2 probs)
  3. embed_contrast.py                   (mist_contrastive_canopus_pretrain.ckpt -> 256-d embedding)
Each candidate is then scored with the hypothesis matching its own formula + adduct
(the way CSI:FingerID scores candidates per formula):
  score_fp       = cosine(predicted fingerprint, candidate Morgan4096 r2)      (MIST retrieval_fp.py, cosine)
  score_contrast = cosine(spectrum embedding, mol_aux_encoder(Morgan4096 r2))  (MIST retrieval_contrast.py)

Output per challenge (atomic; existing files are skipped = resumable):
  <out>/scores/challenge_NNN.tsv.gz   ikey14, formula, adduct, score_fp, score_contrast
  <out>/meta/challenge_NNN.tsv        n_candidates, n_scored, n_hyp, n_hyp_ok, seconds
"""
from __future__ import annotations

import argparse
import os
import pickle
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from jm_common import assign_shards, atomic_write_tsv, report  # noqa: E402

MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=4096)


def morgan4096(smiles: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Morgan r=2, 4096 bits (= MIST 'morgan4096'). Returns (fps float32, ok mask)."""
    fps = np.zeros((len(smiles), 4096), dtype=np.float32)
    ok = np.zeros(len(smiles), dtype=bool)
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        DataStructs.ConvertToNumpyArray(MORGAN.GetFingerprint(m), fps[i])
        ok[i] = True
    return fps, ok


def cosine_rows(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """cosine(a[i], b) for every row of a; b is one vector."""
    den = np.linalg.norm(a, axis=1) * np.linalg.norm(b)
    with np.errstate(invalid="ignore", divide="ignore"):
        return (a @ b) / den


def run(cmd: list[str], log, env) -> None:
    log.write("$ " + " ".join(cmd) + "\n")
    log.flush()
    subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env, check=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--repo", required=True, help="MIST repo dir")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--chunk", type=int, default=10, help="challenges per MIST CLI call")
    ap.add_argument("--mass-diff-thresh", type=float, default=20)
    ap.add_argument("--max-formulae", type=int, default=50)
    ap.add_argument("--keep-work", action="store_true")
    ap.add_argument("--progress-file", default=None)
    args = ap.parse_args()

    inp, out, repo = Path(args.inputs).resolve(), Path(args.out).resolve(), Path(args.repo).resolve()
    for sub in ("scores", "meta", "logs"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    py = str(repo / ".venv/bin/python")
    fp_ckpt = repo / "pretrained_models/mist_fp_canopus_pretrain.ckpt"
    con_ckpt = repo / "pretrained_models/mist_contrastive_canopus_pretrain.ckpt"
    env = dict(os.environ, PYTHONPATH=str(repo / "src"), OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

    # contrastive mol encoder = one linear layer 4096 -> 256 on the Morgan fingerprint
    sd = torch.load(con_ckpt, map_location="cpu")["state_dict"]
    w_mol = sd["mol_aux_encoder.weight"].numpy().astype(np.float32)
    b_mol = sd["mol_aux_encoder.bias"].numpy().astype(np.float32)

    chal = pd.read_csv(inp / "challenges.tsv", sep="\t")
    chal = chal[(chal.polarity == "positive") & (chal.n_hypotheses > 0)]
    hyps = pd.read_csv(inp / "mist_hypotheses.tsv", sep="\t", keep_default_na=False)
    hyps = hyps[hyps.supported.astype(str) == "True"]
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

    log = open(out / "logs" / f"mist_cli_shard{args.shard}.log", "a")
    for c0 in range(0, len(todo), args.chunk):
        chunk = todo[c0:c0 + args.chunk]
        t0 = time.time()
        names = [r.challenge for r in chunk]
        work = out / "work" / f"shard{args.shard}_{names[0]}"
        shutil.rmtree(work, ignore_errors=True)
        (work / "subforms").mkdir(parents=True)
        h = hyps[hyps.challenge.isin(names)]
        prec = dict(zip(chal.challenge, chal.precursor_mz))
        labels = pd.DataFrame({"spec": h.spec, "formula": h.formula, "ionization": h.ionization,
                               "dataset": "casmi2022", "compound": h.spec,
                               "parentmass": h.challenge.map(prec),  # precursor m/z (MIST reads it from the .ms)
                               "instrument": "Orbitrap (LCMS)"})
        labels.to_csv(work / "labels.tsv", sep="\t", index=False)
        tag = work.name
        run([py, str(repo / "src/mist/subformulae/assign_subformulae.py"),
             "--spec-files", str(inp / "mist_ms"), "--labels-file", str(work / "labels.tsv"),
             "--mass-diff-thresh", str(args.mass_diff_thresh), "--max-formulae", str(args.max_formulae),
             "--output-dir", str(work / "subforms"), "--num-workers", "1"], log, env)
        run([py, str(repo / "src/mist/pred_fp.py"), "--num-workers", "0",
             "--labels-file", str(work / "labels.tsv"), "--subform-folder", str(work / "subforms"),
             "--spec-folder", str(inp / "mist_ms"), "--dataset-name", tag,
             "--model-ckpt", str(fp_ckpt), "--save-dir", str(work)], log, env)
        run([py, str(repo / "src/mist/embed_contrast.py"), "--num-workers", "0",
             "--labels-file", str(work / "labels.tsv"), "--subform-folder", str(work / "subforms"),
             "--spec-folder", str(inp / "mist_ms"), "--dataset-name", tag,
             "--model-ckpt", str(con_ckpt), "--save-dir", str(work), "--out-name", "contrast.p"], log, env)
        fp = pickle.load(open(work / f"fp_preds_{tag}.p", "rb"))
        ce = pickle.load(open(work / "contrast.p", "rb"))
        fp_by = dict(zip(map(str, fp["names"]), np.asarray(fp["preds"], dtype=np.float32)))
        ce_by = dict(zip(map(str, ce["names"]), np.asarray(ce["embeds"], dtype=np.float32)))

        for r in chunk:
            t1 = time.time()
            cand = pd.read_csv(inp / "candidates" / f"{r.challenge}.tsv.gz", sep="\t", dtype=str,
                               keep_default_na=False)
            cand["score_fp"], cand["score_contrast"] = np.nan, np.nan
            hc = h[h.challenge == r.challenge]
            n_ok = 0
            for hr in hc.itertuples(index=False):
                if hr.spec not in fp_by or hr.spec not in ce_by:
                    continue
                n_ok += 1
                sel = np.flatnonzero(((cand.formula == hr.formula) & (cand.adduct == hr.sirius_adduct)).values)
                for s0 in range(0, len(sel), 8192):
                    rows = sel[s0:s0 + 8192]
                    fps, ok = morgan4096(cand.smiles.values[rows].tolist())
                    rows, fps = rows[ok], fps[ok]
                    if not len(rows):
                        continue
                    cand.loc[cand.index[rows], "score_fp"] = cosine_rows(fps, fp_by[hr.spec])
                    mol_e = fps @ w_mol.T + b_mol
                    cand.loc[cand.index[rows], "score_contrast"] = cosine_rows(mol_e, ce_by[hr.spec])
            scored = cand.dropna(subset=["score_fp", "score_contrast"], how="all")
            atomic_write_tsv(scored[["ikey14", "formula", "adduct", "score_fp", "score_contrast"]],
                             out / "scores" / f"{r.challenge}.tsv.gz")
            atomic_write_tsv(pd.DataFrame([dict(challenge=r.challenge, n_candidates=len(cand),
                                                n_scored=len(scored), n_hyp=len(hc), n_hyp_ok=n_ok,
                                                seconds=round(time.time() - t1, 2))]),
                             out / "meta" / f"{r.challenge}.tsv")
            units_done += r.n_candidates
            items_done += 1
            report(args.progress_file, units_done, units_total, items_done, len(mine))
        if not args.keep_work:
            shutil.rmtree(work, ignore_errors=True)
        print(f"chunk {names[0]}..{names[-1]} ({len(chunk)} challenges, {len(h)} hypotheses): "
              f"{time.time() - t0:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
