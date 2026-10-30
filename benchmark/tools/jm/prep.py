#!/usr/bin/env python3
"""Build JESTR / MIST inputs for the full CASMI 2022 run (blind).

Reads only model inputs, never the answer key:
  * spectra   : casmi2022/sirius_run/ms/challenge_NNN.ms (blind extraction,
                measured precursor m/z, 3 collision energies + MS1 block)
  * candidates: SIRIUS PubChem candidates (structure_identifications_all.tsv),
                i.e. the same list CSI:FingerID ranked. Formula and adduct of
                each candidate come from SIRIUS, not from CASMI.

Outputs (--out):
  challenges.tsv              one row per challenge
  spectra/challenge_NNN.tsv   merged MS2 (3 CEs, 10 ppm clusters, max intensity) -> JESTR
  candidates/challenge_NNN.tsv.gz
                              ikey14, smiles, formula, adduct, csi_score, csi_rank
                              (dedup by InChIKey first block, best SIRIUS rank kept,
                              same rule as sirius_run/scripts/04_score.py)
  mist_ms/<hyp>.ms            one SIRIUS-style .ms per (challenge, formula, adduct)
                              hypothesis, MS2 blocks only (positive mode only)
  mist_hypotheses.tsv         hyp name -> challenge, formula, MIST ionization
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from jm_common import report  # noqa: E402

ID_RE = r"_(challenge_\d+)_challenge"
USECOLS = ["mappingFeatureId", "structurePerIdRank", "CSI:FingerIDScore",
           "molecularFormula", "adduct", "InChIkey2D", "smiles"]
MIST_IONS = {"[M+H]+", "[M+Na]+", "[M+K]+", "[M-H2O+H]+", "[M+H3N+H]+", "[M]+", "[M-H4O2+H]+"}
MERGE_PPM = 10.0  # same MS2 tolerance as the SIRIUS run


def parse_ms(path: Path) -> dict:
    """SIRIUS .ms -> header dict + list of (block name, peaks array)."""
    meta, blocks, cur_name, cur = {}, [], None, []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            key, _, val = line[1:].partition(" ")
            if key in ("ms1", "ms2", "collision"):
                if cur_name is not None:
                    blocks.append((cur_name, np.array(cur, dtype=float).reshape(-1, 2)))
                cur_name, cur = line[1:], []
            else:
                meta[key] = val.strip()
        else:
            mz, it = line.split()[:2]
            cur.append((float(mz), float(it)))
    if cur_name is not None:
        blocks.append((cur_name, np.array(cur, dtype=float).reshape(-1, 2)))
    return {"meta": meta, "blocks": blocks}


def merge_peaks(peaks: np.ndarray, ppm: float = MERGE_PPM) -> np.ndarray:
    """Cluster peaks within `ppm`; keep the m/z of the strongest peak and its intensity."""
    peaks = peaks[peaks[:, 1] > 0]
    if len(peaks) == 0:
        return peaks
    peaks = peaks[np.argsort(peaks[:, 0])]
    out, start = [], 0
    for i in range(1, len(peaks) + 1):
        if i == len(peaks) or peaks[i, 0] - peaks[start, 0] > peaks[start, 0] * ppm * 1e-6:
            grp = peaks[start:i]
            out.append(grp[np.argmax(grp[:, 1])])
            start = i
    return np.vstack(out)


def mist_ion(adduct: str) -> str:
    return str(adduct).replace(" ", "")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ms-dir", required=True)
    ap.add_argument("--summaries", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ids", default="", help="comma-separated challenge numbers (default: all .ms)")
    ap.add_argument("--progress-file", default=None)
    args = ap.parse_args()

    out = Path(args.out)
    for sub in ("spectra", "candidates", "mist_ms"):
        (out / sub).mkdir(parents=True, exist_ok=True)

    ms_files = sorted(Path(args.ms_dir).glob("challenge_*.ms"))
    if args.ids:
        keep = {int(x) for x in args.ids.split(",") if x.strip()}
        ms_files = [p for p in ms_files if int(p.stem.split("_")[1]) in keep]
    names = [p.stem for p in ms_files]
    name_set = set(names)

    # ---- spectra
    chal = {}
    for p in ms_files:
        ms = parse_ms(p)
        ion = ms["meta"].get("ionization", "")
        pol = "negative" if ion.endswith("-") else "positive"
        ms2 = [(n, a) for n, a in ms["blocks"] if not n.startswith("ms1") and len(a)]
        merged = merge_peaks(np.vstack([a for _, a in ms2])) if ms2 else np.zeros((0, 2))
        np.savetxt(out / "spectra" / f"{p.stem}.tsv", merged, fmt="%.6f\t%.4f",
                   header="mz\tintensity", comments="")
        chal[p.stem] = dict(challenge=p.stem, challenge_id=int(p.stem.split("_")[1]), polarity=pol,
                            precursor_mz=float(ms["meta"]["parentmass"]), n_ce=len(ms2),
                            n_peaks_merged=len(merged), ms=ms, ms2=ms2)

    # ---- candidates (1.4 GB file, chunked; QUOTE_NONE because names contain stray quotes)
    src = Path(args.summaries) / "structure_identifications_all.tsv"
    total_bytes = src.stat().st_size
    parts = []
    with open(src, "rb") as fh:
        reader = pd.read_csv(fh, sep="\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE,
                             usecols=USECOLS, chunksize=200_000)
        for ch in reader:
            ch["challenge"] = ch["mappingFeatureId"].str.extract(ID_RE)[0]
            parts.append(ch[ch["challenge"].isin(name_set)].drop(columns="mappingFeatureId"))
            report(args.progress_file, 600 * fh.tell() / total_bytes, 1000, 0, len(names))  # read = 60%
    cand = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=USECOLS + ["challenge"])
    cand["csi_rank"] = pd.to_numeric(cand["structurePerIdRank"])
    cand = (cand.sort_values(["challenge", "csi_rank"], kind="stable")
                .drop_duplicates(["challenge", "InChIkey2D"], keep="first")
                .rename(columns={"InChIkey2D": "ikey14", "molecularFormula": "formula",
                                 "CSI:FingerIDScore": "csi_score"}))
    cols = ["ikey14", "smiles", "formula", "adduct", "csi_score", "csi_rank"]
    groups = dict(tuple(cand.groupby("challenge", sort=False)))

    # ---- write per challenge
    hyp_rows = []
    for k, name in enumerate(names, 1):
        c = chal[name]
        g = groups.get(name)
        g = g[cols] if g is not None else pd.DataFrame(columns=cols)
        g.to_csv(out / "candidates" / f"{name}.tsv.gz", sep="\t", index=False)
        c["n_candidates"] = len(g)
        c["n_hypotheses"] = 0
        if c["polarity"] == "positive" and len(g) and c["ms2"]:
            hyps = g.groupby(["formula", "adduct"], sort=False).size().reset_index(name="n")
            for j, h in enumerate(hyps.itertuples(index=False)):
                ion = mist_ion(h.adduct)
                spec = f"{name}__h{j:02d}"
                supported = ion in MIST_IONS
                hyp_rows.append(dict(spec=spec, challenge=name, formula=h.formula, ionization=ion,
                                     sirius_adduct=h.adduct, n_candidates=int(h.n), supported=supported))
                if not supported:
                    continue
                lines = [f">compound {spec}", f">formula {h.formula}",
                         f">parentmass {c['precursor_mz']:.5f}", f">ionization {ion}", ""]
                for bname, arr in c["ms2"]:
                    lines.append(f">{bname}")
                    lines += [f"{m:.6f} {i:.2f}" for m, i in arr]
                    lines.append("")
                (out / "mist_ms" / f"{spec}.ms").write_text("\n".join(lines))
                c["n_hypotheses"] += 1
        report(args.progress_file, 600 + 400 * k / len(names), 1000, k, len(names))  # write = 40%

    tab = pd.DataFrame([{k: v for k, v in c.items() if k not in ("ms", "ms2")} for c in chal.values()])
    tab.sort_values("challenge_id").to_csv(out / "challenges.tsv", sep="\t", index=False)
    hyp_cols = ["spec", "challenge", "formula", "ionization", "sirius_adduct", "n_candidates", "supported"]
    pd.DataFrame(hyp_rows, columns=hyp_cols).to_csv(out / "mist_hypotheses.tsv", sep="\t", index=False)

    pos = tab[tab.polarity == "positive"]
    print(f"challenges: {len(tab)} (pos {len(pos)}, neg {len(tab) - len(pos)}); "
          f"candidates: {int(tab.n_candidates.sum())} (0 for {int((tab.n_candidates == 0).sum())} challenges); "
          f"MIST hypotheses: {int(tab.n_hypotheses.sum())}")


if __name__ == "__main__":
    main()
