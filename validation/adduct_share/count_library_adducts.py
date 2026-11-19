"""Share of MS/MS spectra in MassBank and GNPS whose adduct is not [M+H]+ or [M-H]-.

Reads mzDesk's local library indexes (~/.mzdesk/data/libraries/<name>/, or MZDESK_LIBRARY_DIR).
Counts MS/MS records only (kind == 0) and leaves out GNPS copies of MassBank records (copy_flag == 1).
Adduct names are the canonical names of the index vocabulary; -1 means the record gives no usable adduct.

Two shares are printed:
  among spectra with a known adduct;
  a lower bound over all spectra, counting every spectrum with an unknown adduct as [M+H]+ or [M-H]-.

Usage: python3 validation/adduct_share/count_library_adducts.py
"""
import collections
import json
import os
from pathlib import Path

import numpy as np

PROTONATED = ("[M+H]+", "[M-H]-")


def library_dir():
    return Path(os.environ.get("MZDESK_LIBRARY_DIR") or Path.home() / ".mzdesk" / "data" / "libraries")


def count(lib):
    d = library_dir() / lib
    meta = json.loads((d / "meta.json").read_text())
    vocab = meta["adducts"]
    kind = np.load(d / "kind.npy", mmap_mode="r")
    adduct = np.load(d / "adduct_id.npy", mmap_mode="r")
    copy = np.load(d / "copy_flag.npy", mmap_mode="r")
    keep = (kind == 0) & (copy == 0)
    ids = np.asarray(adduct[keep])
    counts = collections.Counter(vocab[i] for i in ids[ids >= 0].tolist())
    return meta.get("version"), int(keep.sum()), int((ids < 0).sum()), counts


def main():
    total = collections.Counter()
    n_all = n_unknown = 0
    for lib in ("MassBank", "GNPS"):
        version, n, unknown, counts = count(lib)
        total.update(counts)
        n_all += n
        n_unknown += unknown
        print(f"{lib} {version}: {n} MS/MS spectra, {unknown} with unknown adduct")
    known = sum(total.values())
    protonated = sum(total[a] for a in PROTONATED)
    print(f"All: {n_all} MS/MS spectra, {known} with a known adduct")
    print(f"Not [M+H]+ or [M-H]-, among known: {100 * (known - protonated) / known:.1f} %")
    print(f"Not [M+H]+ or [M-H]-, lower bound over all: {100 * (known - protonated) / n_all:.1f} %")
    print("Most common:", ", ".join(f"{a} {c}" for a, c in total.most_common(10)))


if __name__ == "__main__":
    main()
