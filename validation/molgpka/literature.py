"""Predicted pKa (numpy MolGpKa engine) against literature pKa for ten common drugs / amino acids.

Run:  python validation/molgpka/literature.py            Out: literature.csv; prints the table and the summary.

The literature values are those given in the task brief (approximate, macroscopic, mostly 25 C in water); they were not
re-checked against the primary sources. Each literature value is assigned to a site by a SMARTS for the functional group
(site = atom 0 of the match), so the comparison is per group, not "closest value". Nothing here was used to tune anything:
the engine has no parameter that depends on these molecules.

Caveat: MolGpKa was trained on ChEMBL pKa data, which very likely contains several of these drugs, so agreement here is not
a held-out accuracy estimate.
"""
import csv
import sys
from pathlib import Path

from rdkit import Chem

from mzdesk.platform_fit import pka_engine as E

HERE = Path(__file__).parent

# (compound, SMILES, group SMARTS, kind, literature low, literature high, note)
CASES = [
    ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "[OX2H1][CX3]=O", "acid", 4.4, 4.9, "carboxylic acid"),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O", "[OX2H1][CX3]=O", "acid", 3.5, 3.5, "carboxylic acid"),
    ("lidocaine", "CCN(CC)CC(=O)Nc1c(C)cccc1C", "[NX3;H0;!$(N-C=O);!$(N-c)]", "base", 7.9, 7.9, "tertiary amine"),
    ("morphine (amine)", "CN1CC[C@]23[C@@H]4[C@H]1Cc5ccc(O)c(O[C@H]2[C@@H](O)C=C4)c35", "[NX3;H0;!$(N-C=O);!$(N-c)]", "base", 8.2, 8.2, "tertiary amine"),
    ("morphine (phenol)", "CN1CC[C@]23[C@@H]4[C@H]1Cc5ccc(O)c(O[C@H]2[C@@H](O)C=C4)c35", "[OX2H1]c", "acid", 9.9, 9.9, "phenol"),
    ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "[nX2;r5]", "base", 0.6, 0.6, "imidazole N (protonation site)"),
    ("glycine (acid)", "NCC(=O)O", "[OX2H1][CX3]=O", "acid", 2.3, 2.3, "carboxylic acid"),
    ("glycine (amine)", "NCC(=O)O", "[NX3;H2][CX4]", "base", 9.6, 9.6, "primary amine"),
    ("sulfamethoxazole (aniline)", "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1", "[NX3;H2]c", "base", 1.8, 1.8, "aniline N"),
    ("sulfamethoxazole (sulfonamide)", "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1", "[NX3;H1](S(=O)=O)", "acid", 5.7, 5.7, "sulfonamide NH"),
    ("metformin", "CN(C)C(=N)NC(N)=N", "[NX3,NX2;$(N~C(~N)~N)]", "base", 12.4, 12.4, "biguanide, strongest base site"),
    ("salicylic acid", "OC(=O)c1ccccc1O", "[OX2H1][CX3]=O", "acid", 3.0, 3.0, "carboxylic acid"),
    ("procainamide", "CCN(CC)CCNC(=O)c1ccc(N)cc1", "[NX3;H0;!$(N-C=O);!$(N-c)]", "base", 9.2, 9.2, "tertiary amine"),
]


def predict(smiles, smarts, kind):
    """pKa of the site of ``kind`` on the atoms matched by ``smarts`` (the highest one when several atoms match, e.g. biguanide)."""
    prep = E.prepare(smiles)
    atoms = {m[0] for m in prep.mol.GetSubstructMatches(Chem.MolFromSmarts(smarts))}
    sites, _ = E.drop_nonbasic_sites(prep.mol, E.predict_sites(prep.mol))          # the sites the page uses
    hits = [s for s in sites if s.atom in atoms and s.kind == kind]
    if not hits:
        return None
    return max(hits, key=lambda s: s.pka)


def error(pred, lo, hi):
    """Signed error against the literature value; for a range, 0 inside it and the distance to the nearest end outside."""
    if lo <= pred <= hi:
        return 0.0
    return pred - lo if pred < lo else pred - hi


def main():
    rows = []
    for name, smiles, smarts, kind, lo, hi, note in CASES:
        site = predict(smiles, smarts, kind)
        pred = None if site is None else site.pka
        lit = f"{lo:g}" if lo == hi else f"{lo:g}–{hi:g}"
        rows.append({"compound": name, "site": "" if site is None else site.label, "type": kind, "group": note,
                     "predicted": "" if pred is None else f"{pred:.2f}", "literature": lit,
                     "error": "" if pred is None else f"{error(pred, lo, hi):+.2f}"})
    with open(HERE / "literature.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"{'compound':32s} {'site':5s} {'type':5s} {'predicted':>9s} {'literature':>10s} {'error':>6s}")
    for r in rows:
        print(f"{r['compound']:32s} {r['site']:5s} {r['type']:5s} {r['predicted']:>9s} {r['literature']:>10s} {r['error']:>6s}")
    errs = [abs(float(r["error"])) for r in rows if r["error"]]
    print(f"\n{len(errs)} sites: MAE {sum(errs) / len(errs):.2f}, median {sorted(errs)[len(errs) // 2]:.2f}, max {max(errs):.2f}; "
          f"within 1.0: {sum(e <= 1.0 for e in errs)} / {len(errs)}")


if __name__ == "__main__":
    sys.exit(main())
