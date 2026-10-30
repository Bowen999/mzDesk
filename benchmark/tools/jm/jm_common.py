"""Helpers shared by the JESTR / MIST runners (no third-party imports beyond pandas)."""
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd


def report(path: str | None, done: float, total: float, items_done: int, items_total: int) -> None:
    """Write '<units done> <units total> <items done> <items total>' for j_m_run.sh's progress bar."""
    if not path:
        return
    tmp = f"{path}.tmp"
    with open(tmp, "w") as fh:
        fh.write(f"{done:.0f} {total:.0f} {items_done} {items_total}\n")
    os.replace(tmp, path)


def assign_shards(chal: pd.DataFrame, nshards: int, weight: str = "n_candidates") -> list[pd.DataFrame]:
    """Deterministic greedy balancing by candidate count (largest first)."""
    order = chal.sort_values([weight, "challenge_id"], ascending=[False, True])
    loads, members = [0] * nshards, [[] for _ in range(nshards)]
    for idx, w in zip(order.index, order[weight]):
        k = loads.index(min(loads))
        loads[k] += int(w) + 50  # +50: fixed per-challenge cost
        members[k].append(idx)
    return [chal.loc[sorted(m, key=lambda i: chal.at[i, "challenge_id"])] for m in members]


def atomic_write_tsv(df: pd.DataFrame, path: Path) -> None:
    path = Path(path)
    suffix = ".tsv.gz" if path.name.endswith(".gz") else ".tsv"
    tmp = path.with_name(path.name + ".part" + suffix)
    df.to_csv(tmp, sep="\t", index=False, compression="gzip" if suffix.endswith("gz") else None)
    os.replace(tmp, path)
