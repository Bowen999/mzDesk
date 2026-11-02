# MS/MS Prediction

Predict the ESI MS/MS spectrum of a structure at 10 / 20 / 40 eV (CFM-ID 4; MSRB rules for lipids). Needs Docker and the wishartlab/cfmid image; the container has no network and runs with no extra privileges.

## Inputs

- **Compound**: A salt keeps its largest component; a charged structure is neutralised before the adduct is added.
- **Adduct**: MSRB covers seven lipid adducts; the CFM-ID ML fallback only supports [M+H]+ and [M-H]-. Other adducts return a clear error instead of a different ion.

## Results

### Peaks

- Fragment: formula · SMILES from the engine
