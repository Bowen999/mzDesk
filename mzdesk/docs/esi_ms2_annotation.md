# MS/MS Annotation

I have an ESI MS/MS spectrum and a precursor m/z. Find the adduct and formula, search the libraries and, if no library spectrum matches, propose structures; then let detectability order the near-ties (candidates within 10% of the best score) and check which isomers the spectrum can tell apart. If you also have an observed CCS, each candidate is checked against its reference or predicted CCS.

## Inputs

- **Mobile phase additive**: Unknown = the best additive counts for each candidate's adduct.
- **Sodium / potassium present**: Glassware, salts, biological matrix. Unknown counts as present.
- **Mobile phase pH**: For the ionization check of each candidate (default: 0.1 % formic acid).
- **Library similarity algorithms**: Tick one or more. Each library hit is scored with every ticked algorithm and passes the Level-2 gate when at least one reaches the similarity threshold; the first ticked algorithm orders the hits. One that cannot be computed is left out (CosineGreedy replaces all of them if none can be).
- **Similarity threshold**: A library hit needs this similarity, on at least one ticked algorithm, to pass the Level-2 gate. Hits join the candidate pool from this value minus 0.10.
- **Observed CCS (Å²)**: Optional. One value measured in N₂. Each candidate is compared with its CCS in the local reference data (same compound and adduct), or with a predicted value if there is none. A candidate that disagrees with a measured value goes behind the others; a predicted value only reorders near-ties. For TIMS enter the CCS your vendor software reports, not 1/K0.
- **Platform**: Required when a CCS is given: the allowed difference depends on the platform.

## Results

### Run notes

- The isomer verdict only holds for the candidates that were considered; an isomer that is not in the candidate list has not been ruled out. When the isomer check was not run or could not be done, other isomers are not ruled out either.
- A formula answer (no structure proposed) is formula-level evidence only.
- FRIGID candidates are hypotheses, not identifications.
- Run FRIGID (Advanced): on by default in the local version, off on mzdesk.com, where a run takes about two minutes and 3 GB of the server's 4 GB. Tick it to run FRIGID when the library and MS-FINDER do not settle the structure.
- In the Formula table, Rank is the list position, not a joint probability: the list is interleaved by BUDDY rank across the adducts (best adduct first within each rank). Filtered formulas can leave gaps in BUDDY rank.
- Candidates within 10 % of the best base score are ordered by base score × detectability; the others follow in base-score order, so the Final score column does not always go down the list.
- CCS check: a candidate that disagrees with a measured reference goes behind all the others; one that disagrees with a predicted value only goes behind the other near-ties (within 10 % of the best base score). Levels do not change.
