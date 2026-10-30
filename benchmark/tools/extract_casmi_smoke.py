#!/usr/bin/env python3
"""Extract 10 CASMI2022 challenges into an MGF + labels table for smoke tests.

Selection: 5 x [M+H]+ (positive) and 5 x [M-H]- (negative), one challenge per
source file, monoisotopic mass < 500, distinct molecular formulae.

For each challenge the mzML run is scanned once for MS2 scans; the scan whose
(precursor m/z, retention time) is closest to the challenge annotation wins.

Outputs (in --out-dir):
  smoke.mgf        MGF with 10 MS2 spectra (FEATURE_ID = challenge id)
  labels.tsv       MIST-style label table
  selected.csv     provenance table (challenge id, scan, rt, precursor, ...)
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd
from pyteomics import mzml, mgf


def pick_challenges(meta: pd.DataFrame, n_pos: int = 5, n_neg: int = 5, min_mass: float = 150.0) -> pd.DataFrame:
    """Pick small, singly-charged, formula-distinct challenges, one per file."""
    meta = meta.dropna(subset=["smiles", "formula", "precursor_mz", "rt_min"])
    meta = meta[(meta["monoisotopic_mass"] >= min_mass) & (meta["monoisotopic_mass"] < 500)]

    picked = []
    for polarity, adduct, n in (
        ("positive", "[M+H]+", n_pos),
        ("negative", "[M-H]-", n_neg),
    ):
        pool = meta[(meta.polarity == polarity) & (meta.adduct.str.strip() == adduct)]
        pool = pool.drop_duplicates(subset=["file"]).drop_duplicates(subset=["formula"])
        pool = pool.sort_values(["monoisotopic_mass", "challenge_id"])
        picked.append(pool.head(n))
    return pd.concat(picked).reset_index(drop=True)


def scan_mzml(path: Path, target_mz: float, target_rt: float):
    """Return (scan_dict, mz_array, intensity_array) closest to the target."""
    best = None
    with mzml.MzML(str(path), use_index=False) as fh:
        for spec in fh:
            if spec.get("ms level") != 2:
                continue
            try:
                prec = spec["precursorList"]["precursor"][0]["selectedIonList"][
                    "selectedIon"
                ][0]["selected ion m/z"]
                rt = float(spec["scanList"]["scan"][0]["scan start time"])
            except (KeyError, IndexError, TypeError):
                continue
            dmz = abs(float(prec) - target_mz)
            if dmz > 0.02:  # 20 mDa window keeps the search honest
                continue
            score = dmz / 0.02 + abs(rt - target_rt) / 0.5
            if best is None or score < best[0]:
                best = (score, spec, prec, rt)
    if best is None:
        return None
    _, spec, prec, rt = best
    mz = spec["m/z array"]
    intens = spec["intensity array"]
    keep = intens > 0
    return (
        {"scan": spec.get("index"), "rt": rt, "precursor": float(prec)},
        mz[keep],
        intens[keep],
    )


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--casmi", required=True, help="casmi2022 directory")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--n-pos", type=int, default=5)
    ap.add_argument("--n-neg", type=int, default=5)
    ap.add_argument("--min-mass", type=float, default=150.0)
    args = ap.parse_args()

    casmi = Path(args.casmi)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(casmi / "metadata/casmi2022_challenges.csv")
    picked = pick_challenges(meta, args.n_pos, args.n_neg, args.min_mass)
    print(f"selected {len(picked)} challenges")

    selected_rows, mgf_specs, label_rows = [], [], []
    for _, row in picked.iterrows():
        mzml_path = casmi / row.mzml_path
        hit = scan_mzml(mzml_path, float(row.precursor_mz), float(row.rt_min))
        if hit is None:
            print(f"  !! challenge {row.challenge_id}: no matching MS2 scan")
            continue
        scan, mz, intens = hit
        feat = f"casmi{int(row.challenge_id):03d}"
        params = {
            "FEATURE_ID": feat,
            "TITLE": feat,
            "PEPMASS": float(row.precursor_mz),
            "CHARGE": "1+",
            "IONMODE": row.polarity,
            "RTINSECONDS": round(scan["rt"] * 60, 2),
            "EQID": row.file,
            "COMPOUND": str(row.formula),
        }
        mgf_specs.append({"params": params, "m/z array": mz, "intensity array": intens})
        label_rows.append(
            {
                "spec": feat,
                "formula": row.formula,
                "ionization": row.adduct,
                "dataset": "casmi2022",
                "compound": feat,
                "parentmass": round(float(row.monoisotopic_mass), 6),
                "instrument": "unknown",
            }
        )
        selected_rows.append(
            {
                "challenge_id": int(row.challenge_id),
                "feature_id": feat,
                "file": row.file,
                "polarity": row.polarity,
                "adduct": row.adduct,
                "formula": row.formula,
                "neutral_mass": float(row.monoisotopic_mass),
                "precursor_mz": float(row.precursor_mz),
                "rt_min": float(row.rt_min),
                "matched_scan": scan["scan"],
                "matched_rt_min": round(scan["rt"], 3),
                "matched_precursor_mz": round(scan["precursor"], 5),
                "n_peaks": int(len(mz)),
                "smiles": row.smiles,
                "inchikey": row.inchikey,
                "mzml_path": str(mzml_path),
            }
        )
        print(
            f"  {feat}: {row.formula:12s} {row.adduct:8s} scan={scan['scan']} "
            f"rt={scan['rt']:.3f}/{row.rt_min} prec={scan['precursor']:.4f}"
            f"/{row.precursor_mz} peaks={len(mz)}"
        )

    with open(out / "smoke.mgf", "w") as fh:
        mgf.write(list(mgf_specs), fh)
    pd.DataFrame(label_rows).to_csv(out / "labels.tsv", sep="\t", index=False)
    pd.DataFrame(selected_rows).to_csv(out / "selected.csv", index=False)
    print(f"wrote {out/'smoke.mgf'} ({len(mgf_specs)} spectra), {out/'labels.tsv'}")


if __name__ == "__main__":
    main()
