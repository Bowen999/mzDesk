# MolGpKa (pKa prediction) validation: numpy rewrite vs the original code

`mzdesk/platform_fit/pka_engine.py` re-implements the forward pass of MolGpKa in numpy (no PyTorch, no torch_geometric at run time).
This folder compares it with the **original code and weights** on the same molecules, and with literature pKa values.

- Original code: <https://github.com/Xundrug/MolGpKa>, commit `4dc8352` (2024-01-11), **MIT licence** (the LICENSE.md copyright line reads "2017-Present OpenNMT", copied from a template; the grant is MIT)
- Weights: `models/weight_acid.pth`, `models/weight_base.pth` in that repository (5.97 MB each; sha256 pinned in `pka_setup.py`)
- Paper: Pan, Wang, Yang, Zhang, J. Chem. Inf. Model. 2021, 61, 3159-3165

## Reproduce

```bash
python -m mzdesk.platform_fit.pka_setup       # download the two .pth, verify sha256, convert to ~/.mzdesk/models/molgpka/molgpka_weights.npz (needs torch once)
PY=/opt/anaconda3/bin/python3 validation/molgpka/run_original.sh     # original code in a throw-away venv (torch_geometric 2.4.0 only there) -> original_output.json
python validation/molgpka/compare.py          # -> comparison.csv, summary.json, tests/fixtures/molgpka_reference.json
python validation/molgpka/literature.py       # -> literature.csv
<python> -I validation/molgpka/rdkit_versions.py       # RDKit-version check, run with other Pythons / conda envs
```

`run_original.sh` clones the repository, builds `venv --system-site-packages` (borrows torch, pandas and RDKit of the Python you give in `PY`; only `torch_geometric==2.4.0` is installed, inside the venv) and runs the **unmodified** original modules (`predict_pka`, `utils.descriptor`, `utils.net`, `utils.ionization_group`). Two adaptations around them, neither changes a computation: a 3-line stand-in for `torch_scatter.scatter_add` (the original `gcn_conv.py` only uses that function), and the path of `smarts_pattern.tsv` (the original finds it relative to the current directory).

## What the original does (the numpy code follows it)

| Step | Original |
|---|---|
| Sites | 143 SMARTS (`smarts_pattern.tsv`, 54 acid / 89 base) on the molecule with explicit H. Acid site = the H atom that is lost; base site = the heavy atom that takes the proton |
| Features | 29 per atom: element one-hot (C H O N S Cl F Br P I; others count as I), degree, hybridization one-hot (5), implicit valence, aromatic, ring size 3-8 (6), H-bond donor, H-bond acceptor, formal charge, topological distance to the site (path length in atoms; 0 for the site), site flag |
| Network | 5 x [GCNConv (D^-1/2 (A+I) D^-1/2 X W + b) -> ReLU -> BatchNorm], global attention pooling, Linear 1024-128-16-1 with ReLU. One forward pass per site, separate weights for acid and base |
| `predict()` | Uncharger, SMILES round trip (renumbers the atoms), AddHs, all sites |

Two things in the original that are easy to miss, and which the rewrite reproduces on purpose:

- The donor / acceptor columns are **always 0**. `descriptor.py` builds lists of 1-tuples from `GetSubstructMatches` and then tests `atom_idx in <that list>`, which is never true. The same code made the training data, so the model learned with those two columns at 0.
- Atoms in another fragment get distance 0 (no path), not "far". Only matters for salts, which the engine reduces to the largest fragment anyway.

## Deliberate differences of the engine (not bugs)

- Atom order of the input is kept (the original's SMILES round trip renumbers atoms), so site labels such as `N12` point to the atom numbers of the user's SMILES.
- The largest fragment is analysed (salts, solvates); anions not next to a cation (carboxylate, phosphate, phenolate ...) are protonated even when a quaternary N+ is present (RDKit's Uncharger would keep them charged, hiding the site in betaines and phosphatidylcholines); then the Uncharger removes the other removable charges. Quaternary N+, N-alkylpyridinium etc. are fixed charges; nitro / N-oxide pairs are not.
- One prediction per heavy atom for acids (NH2 has two equivalent H atoms; the original lists both with the same value).
- After prediction, base sites on amide / sulfonamide / pyrrole-type aromatic N are dropped (see "What the page does" below). `predict_sites` still returns every pattern hit; `analyze` applies the filter.

## Results (349 molecules: the 49 in `molecules.csv` + 300 random from the original training file, seed 20241024)

Original run: torch 2.6.0, torch_geometric 2.4.0, RDKit 2024.09.6, Python 3.12.12.

| Level | Result |
|---|---|
| Network, every atom as the site, acid and base weights (1910 atom predictions, first 40 molecules) | max abs diff **2.1e-05** pKa |
| Site matching on `AddHs(MolFromSmiles(smiles))`, input atom order (349 molecules) | site sets identical **349 / 349** |
| pKa of those sites (1331 predictions) | max abs diff **1.9e-05** pKa |
| Whole engine path vs original `predict()` (Uncharger, round trip; sites matched by canonical atom rank), molecules for which both ended on the same neutral molecule | 332 molecules, site sets identical **332 / 332**, max abs diff **1.9e-05** |
| Not comparable (the engine takes the largest fragment / protonates anions next to a quaternary N+) | 17 molecules: 16 salts or charged inputs, carnitine, phosphatidylcholine; listed in `summary.json` |

The 300 extra molecules come from the original training file, so they test equivalence, not accuracy.

### RDKit version

The model was trained years before the RDKit used here. The engine's predictions were recomputed under other RDKit versions and compared with the original code's values (`rdkit_versions.py`):

| RDKit | molecules (unparsable) | site sets equal | max abs diff pKa |
|---|---|---|---|
| 2021.09.2 | 345 (4) | 345 / 345 | 9.5e-06 |
| 2023.09.6 | 349 (0) | 349 / 349 | 1.9e-05 |
| 2025.09.4 | 349 (0) | 349 / 349 | 9.5e-06 |
| 2026.03.2 | 349 (0) | 349 / 349 | 1.9e-05 |
| 2024.09.6 (host) | 349 (0) | 349 / 349 | 1.9e-05 |

Not tested: RDKit before 2021 (the original environment used 2019-2020).

### Speed

Import 0.3 s; weights load 6 ms; 4 ms for ibuprofen or lidocaine, 15 ms for a phosphatidylcholine (50 heavy atoms), 64 ms for a 7-residue peptide (64 heavy atoms, 16 sites), 1.6 s for a 30-residue peptide (281 heavy atoms, 76 sites). Peak memory about 110 MB.

## Literature pKa (`literature.py`, `literature.csv`)

Literature values are those given in the task brief (approximate, macroscopic; not re-checked against primary sources). Each is assigned to a site by the functional group (SMARTS), not by "closest value". Nothing in the engine depends on these molecules; no threshold was tuned to them. The non-basic N filter was added afterwards for chemical reasons (it removes sites that are not among these groups), so the table is the same before and after it.

| Compound | Site | Type | Predicted | Literature | Error |
|---|---|---|---|---|---|
| ibuprofen | O15 (COOH) | acid | 4.40 | 4.4-4.9 | 0.00 (inside the range) |
| aspirin | O13 (COOH) | acid | 3.39 | 3.5 | -0.11 |
| lidocaine | N3 (tertiary amine) | base | 7.60 | 7.9 | -0.30 |
| morphine | N2 (tertiary amine) | base | 8.03 | 8.2 | -0.17 |
| morphine | O13 (phenol) | acid | 9.43 | 9.9 | -0.47 |
| caffeine | N4 (imidazole N) | base | 2.06 | ~0.6 | +1.46 |
| glycine | O5 (COOH) | acid | 2.09 | 2.3 | -0.21 |
| glycine | N1 (NH2) | base | 9.64 | 9.6 | +0.04 |
| sulfamethoxazole | N13 (aniline N) | base | 4.13 | ~1.8 | +2.33 |
| sulfamethoxazole | N5 (sulfonamide NH) | acid | 7.02 | 5.7 | +1.32 |
| metformin | N5 (biguanide, highest base site) | base | 9.57 | ~12.4 | -2.83 |
| salicylic acid | O1 (COOH) | acid | 2.97 | 3.0 | -0.03 |
| procainamide | N3 (tertiary amine) | base | 9.09 | 9.2 | -0.11 |

13 groups: mean absolute error 0.72, median 0.21, max 2.83; 9 of 13 within 1.0. All the common acids and amines are within 0.5; the misses are weakly basic heteroaromatic / aniline N (caffeine, sulfamethoxazole) and the strongly basic biguanide (metformin; the model's 12.42 is for its N-H *acid* site, a different quantity). MolGpKa was trained on ChEMBL pKa data that probably includes several of these drugs, so this is not a held-out accuracy estimate.

## What the page does with the pKa values, and where it misleads

Per-site Henderson-Hasselbalch with independent sites; logD = Crippen logP + log10(fraction fully uncharged); ions are assumed not to partition.

**Non-basic N filter.** The 143 patterns also match N atoms that are not basic, and MolGpKa gives them a base pKa (lidocaine's anilide N: 5.2; caffeine's three N-methyl / ring N: 2.0-2.4; sulfamethoxazole's sulfonamide N: 4.3). After prediction, and before the charge states, `drop_nonbasic_sites` (`NONBASIC_N` in `pka_engine.py`) removes **base** sites on N that are

- amide / imide / carbamate / urea / thioamide N: `[NX3][CX3]=[O,S]`
- sulfonamide N: `[NX3]S(=O)=O`
- pyrrole-type aromatic N (lone pair in the aromatic sextet): `[nX3]`, `[nH]`; pyridine-type `[nX2]` stays

Acid sites on the same atoms stay (sulfonamide NH, imide NH, anilide NH). The dropped sites are listed in the page's "Rule hits". The table here is not affected (`literature.py` applies the filter; the groups compared are all kept), and the equivalence checks above compare the unfiltered MolGpKa output.

| Molecule | Base sites before | After the filter |
|---|---|---|
| lidocaine | N3 7.6, N9 (anilide) 5.2 | N3 7.6 (net +1, logD -2.3 at pH 2.7; before: +2, -4.8) |
| caffeine | N4 2.06 and three more at 1.99-2.38 | N4 2.06 (the pyridine-type N, atom index 3) |
| sulfamethoxazole | N13 (aniline) 4.13, N5 (sulfonamide) 4.30 | N13 4.13; N5 stays as an acid, 7.02 |

What still misleads:

- **Weakly basic sites the rules do not cover** stay (aniline N of sulfamethoxazole: predicted 4.13, literature 1.8). Other pattern hits are unchanged too: glucose's five OH groups get acid pKa 11.8-12.2 (harmless at LC pH, but they narrow the "robust pH" ranges above pH 11).
- **Polyprotic groups** (phosphate pKa1 / pKa2) are two independent sites with similar pKa; the real second pKa is 4-5 units higher.
- **Zwitterions** count as charged, so logD of an amino acid is very negative (glycine at pH 6: about -8.5; the measured octanol/water logP is about -3). The page says so in "Rule hits" whenever the dominant species is a zwitterion.
- **Permanent charges** (quaternary N+) leave no uncharged form: logD is not defined, and the page falls back to documented special cases in `separation_rules.py`.

## Files

| File | Content |
|---|---|
| `molecules.csv`, `extra_test_set_smiles.txt` | input molecules |
| `run_original.sh`, `run_original.py` | run the original code (throw-away venv) |
| `compare.py` | comparison; writes `comparison.csv`, `summary.json` and `tests/fixtures/molgpka_reference.json` (6 molecules with the original per-site pKa, used by `tests/test_separation.py`) |
| `literature.py`, `literature.csv` | literature check |
| `rdkit_versions.py` | RDKit-version check |
| `original_output.json` | output of the original code (0.25 MB, git-ignored like the NEIMS one) |
