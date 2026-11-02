# MS/MS Annotation

Give an ESI MS/MS spectrum and its precursor m/z. Get the adduct, the formula and ranked structure candidates, each with a confidence level. An observed CCS is optional and checks every candidate.

## Inputs

- **Mobile phase additive**: Unknown = the best additive counts for each candidate's adduct.
- **Sodium / potassium present**: Glassware, salts, biological matrix. Unknown counts as present.
- **Mobile phase pH**: For the ionization check of each candidate (default: 0.1 % formic acid).
- **Library similarity algorithms**: Tick one or more. Each library hit is scored with every ticked algorithm and passes the Level-2 gate when at least one reaches the similarity threshold; the first ticked algorithm orders the hits. One that cannot be computed is left out (CosineGreedy replaces all of them if none can be).
- **Similarity threshold**: A library hit needs this similarity, on at least one ticked algorithm, to pass the Level-2 gate. Hits join the candidate pool from this value minus 0.10.
- **Observed CCS (Å²)**: Optional. One value measured in N₂. Each candidate is compared with its CCS in the local reference data (same compound and adduct), or with a predicted value if there is none. A candidate that disagrees with a measured value goes behind the others; a predicted value only reorders near-ties. For TIMS enter the CCS your vendor software reports, not 1/K0.
- **Platform**: Required when a CCS is given: the allowed difference depends on the platform.
- **Minimum matched peaks**: A library hit needs at least this many matched peaks to pass the Level-2 gate (default 4; all of them if the spectrum has fewer).
- **Allow Cl/Br in formulas**: Widens the formula search to chlorine and bromine.
- **Compare with database isomers**: The isomer check also compares the top structure with the best-scoring same-formula structures of MS-FINDER's local databases.

## Results

### Run notes

- The isomer verdict covers only the candidates that were compared. When the check did not run, no isomer is ruled out.
- A formula answer (no structure proposed) is formula-level evidence only.
- FRIGID candidates are hypotheses, not identifications.
- Run FRIGID (Advanced): runs FRIGID when the library and MS-FINDER do not settle the structure. In the local version it is ticked once FRIGID is installed; ticking it earlier offers the download. On mzdesk.com it starts unticked, because a run takes about two minutes.
- In the Formula table, Rank is the list position, not a joint probability: the list is interleaved by BUDDY rank across the adducts (best adduct first within each rank). Filtered formulas can leave gaps in BUDDY rank.
- Candidates within 10 % of the best base score are ordered by base score × detectability; the others follow in base-score order, so the Final score column does not always go down the list.
- CCS check: a candidate that disagrees with a measured reference goes behind all the others; one that disagrees with a predicted value only goes behind the other near-ties (within 10 % of the best base score). Levels do not change.
