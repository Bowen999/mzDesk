# MS/MS Prediction

Give a structure and an adduct. Get the predicted ESI MS/MS spectrum at 10 / 20 / 40 eV (CFM-ID 4; MSRB rules for lipids). Needs Docker.

## Inputs

- **Compound**: A salt keeps its largest component; a charged structure is neutralised before the adduct is added.
- **Adduct**: MSRB covers seven lipid adducts; the CFM-ID ML fallback only supports [M+H]+ and [M-H]-. Other adducts return a clear error instead of a different ion.

## Results

### Peaks

- Fragment: formula · SMILES from the engine

### Run notes

- Needs Docker. mzdesk.com has it ready. In the local version, install and start Docker Desktop; mzDesk pulls the CFM-ID image (about 1 GB) the first time.
