# Alternative Explanations

Give a compound or formula and a precursor m/z. Get same-formula compounds, isobaric interferences and isomers. The isomer check also needs an MS2 spectrum.

## Inputs

- **Compound or formula**: Draw or paste a structure, or choose Name / ID / Formula and type a formula (e.g. C8H10N4O2). With only a precursor m/z, only isobaric interferences are searched.
- **Precursor m/z**: Blank = the theoretical m/z of the structure or formula (used for the isobaric search and the isomer check).
- **Databases**: Used by the same-formula and isobaric tables (PubChem). The isomer check searches PubChem independently of this selection.
- **MS2 spectrum (optional, for isomer check)**: Optional. Without it the isomer check is skipped; the two database tables still run.
- **Fragment tolerance (Da, isomer check)**: Absolute fragment tolerance for MetFrag's in-silico fragments (isomer table only; relative tolerance 0 ppm).

## Results

### Same-formula candidates

- Structures PubChem lists with exactly this formula, in its copies of the selected databases (at most the first 100 are examined). One row per connectivity (InChIKey first 14 characters), so stereoisomers and isotopologues fold into one row.
- Databases: which of the selected databases list the structure. PubChem is the source of every row; it is shown when no other selected database has the structure. A run takes a few seconds.

### Isobaric interferences

- Compounds with a different formula whose m/z falls within the tolerance at the chosen adduct, in PubChem's copies of the selected databases; one row per formula, closest first, at most 50.
- Exact m/z and Δ ppm use the compound's closest adduct. Same-formula hits are in the table above, not here.
- At most 100 compounds per database are examined, so a crowded m/z window is not exhaustive.

### Isomers

- PubChem same-formula candidates, independent of the Databases field: at most its first 100 CIDs are retrieved, and at most 9 other usable connectivities are compared with your compound. The run reports truncation; the set is not exhaustive. Structures with the same connectivity are collapsed, so stereoisomers are not assessed.
- Distinguishable: MetFrag predicts an observed peak for your compound but not the isomer (local in-silico fragmentation, tree depth 2). This is model-dependent and does not establish a unique identity; not-scored candidates are not ruled out. Needs Java (the local install includes one). The comparison is limited to 120 seconds.
- Peaks: at least 1 % of the base peak, without the precursor and its isotope peaks; a distinguishing peak must be at least 5 %. Explained peaks: how many of those peaks the structure explains · their share of the intensity.
- For structures recognised as lipids this check is skipped: conventional MS/MS usually cannot reliably determine C=C or sn-positions; specific methods or additional evidence may distinguish them.
- Fragments: View opens that structure with the peaks MetFrag explains for it, each fragment in its own colour on the structure and in the spectrum. It is computed when you click (a few seconds, with a progress bar), with the same peaks and tolerance as this table. Where a fragment fits several places on the structure, one is drawn and the card says how many others fit equally well.
- **Explained peaks**: Peaks MetFrag explains for this structure, and their share of the intensity.
- **Distinguishable by MS2**: yes: a peak of at least 5 % is explained for your compound but not for this isomer. A model result, not proof.
- **Distinguishing fragments**: Those peaks, each with the fragment that explains it.
- **Atoms covered**: Share of this structure's non-hydrogen atoms inside the coloured fragments (the 8 strongest explained peaks): how much of it the MS2 accounts for.
- **Fragments**: View: the explained peaks drawn on the structure. Computed when you click.

### Run notes

- PubChem structures that are unavailable or unusable, and candidates MetFrag could not score, are not ruled out. The same holds when the PubChem search fails.
