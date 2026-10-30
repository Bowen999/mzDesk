#!/usr/bin/env python3
"""Convert the smoke-test MGF into SIRIUS-style .ms files (one per spectrum).

MIST's prediction CLI loads spectra from a folder of ``<spec>.ms`` files where
``<spec>`` matches the ``spec`` column of the labels table.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from pyteomics import mgf


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mgf", required=True)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    labels = pd.read_csv(args.labels, sep="\t").astype(str)
    meta = {r.spec: r for r in labels.itertuples(index=False)}
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    n = 0
    for spec in mgf.read(args.mgf):
        params = spec["params"]
        name = params["feature_id"]
        row = meta[name]
        mz = spec["m/z array"]
        inten = spec["intensity array"]
        order = inten.argsort()[::-1]
        mz, inten = mz[order], inten[order]

        lines = [
            f">compound {name}",
            f">parentmass {row.parentmass}",
            f">ionization {row.ionization}",
            f">formula {row.formula}",
            f">instrument {row.instrument}",
            f">feature_id {name}",
            "",
            f">spectrum {name}",
        ]
        lines += [f"{m:.6f} {i:.4f}" for m, i in zip(mz, inten)]
        (out / f"{name}.ms").write_text("\n".join(lines) + "\n")
        n += 1
    print(f"wrote {n} .ms files to {out}")


if __name__ == "__main__":
    main()
