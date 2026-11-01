# Mass Search

Which database compounds match the precursor m/z, optionally filtered by the isotope pattern.

## Results

### Summary

- **Compounds**: One row per PubChem compound (CID). Stereoisomers and salts have their own CIDs, so they count separately.
- **Formulas**: Different molecular formulas among the hits.
- **Closest hit**: The compound with the smallest |Δ ppm|; among equals, one listed in a chosen database comes first.
- **Databases searched**: Compounds PubChem finds in the mass window, per database (before the 100-per-database limit).

### Exact-mass hits

- Every database is searched through PubChem, so hits come from PubChem's copy of each database (HMDB, DrugBank, ChEBI, KEGG, LIPID MAPS; NORMAN is its Suspect List Exchange, which includes SusDat).
- One group per formula and adduct, closest first. Δ ppm = (observed m/z − theoretical m/z) / theoretical m/z × 10<sup>6</sup>; the theoretical m/z comes from the compound's monoisotopic mass and the adduct, as on the Theoretical page.
- Auto tries every adduct of the ion mode. A compound is listed once, under its best-fitting adduct. Permanently charged compounds (e.g. choline) only match [M]+.
- At most 100 compounds are listed per database. PubChem alone can have thousands: compounds deposited by more PubChem sources come first (depositor count), then lower CIDs (older, usually better known).
- IDs are the databases' own IDs, taken from PubChem cross-references; PubChem has no ID for NORMAN.
- With an isotope MS, each group's theoretical pattern (IsoSpecPy, merged to the chosen resolution) is compared with it: every theoretical peak takes the most intense observed peak within the mass tolerance (Unit mass: 0.5 Da), each observed peak once; the score is the cosine between theoretical abundances and matched intensities, and a theoretical peak with no match counts as 0.
- Observed peaks that match no theoretical peak are ignored (they may be other ions). Groups below the threshold are hidden; groups whose formula cannot be scored (e.g. isotope labels) are kept and marked N/A.
