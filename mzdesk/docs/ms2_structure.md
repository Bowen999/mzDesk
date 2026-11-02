# Structure Prediction

Ranked candidate structures from MS-FINDER (offline, local structure DB). Needs Mono (python -m mzdesk.ms2.msfinder_setup).

## Inputs

- **Precursor m/z**: MS-FINDER searches the formula from this m/z and adduct. Auto = [M+H]+ in positive mode, [M-H]- in negative mode.
- **MS1 isotope pattern (m/z intensity per line, optional)**: Optional isotope pattern; the first peak within 0.5 mDa of the precursor is used.
- **Include halogens (F, Cl, Br, I)**: F, Cl, Br and I are off by default because they widen the formula search.
- **Include Si**: Si is off by default; tick it for silicon-containing compounds.
- **Extended element range (PFAS / unusual compositions)**: Switches MS-FINDER from CommonRange to ExtendedRange, which allows more F per carbon (PFAS) and other unusual compositions. Slower and more false formulas.

## Results

### Structure candidates

- Candidates come from MS-FINDER's local structure databases (HMDB, ChEBI, PubChem-derived, COCONUT, UNPD and others). Nothing is sent online, and a compound that is not in them cannot be returned.
- Score: MS-FINDER's total score from in-silico fragmentation (hydrogen-rearrangement rules, tree depth 2) against the MS2 peaks; higher is better. Candidates of different formulas are listed together; the Formula column shows which one each belongs to.
- Structures are searched for the 5 best formulas only. The table shows the best 50. Databases: where the structure is listed; the ? shows the IDs.

### Formula candidates

- Formula score: mass accuracy, the isotope pattern (when MS1 is given) and how well the MS2 fragments and neutral losses fit. Only formulas with a positive score are listed (best 50).
- Search limits of the shipped parameter file: precursor tolerance 1 mDa (isotope peaks must sit within 0.5 mDa of the precursor m/z), fragment tolerance 10 mDa, elements C H N O P S. F, Cl, Br, I and Si are searched only when you tick Include halogens / Include Si (each off by default; the Elements searched card shows what was used).
- Mass error (mDa) is observed − theoretical (the same sign as Δ ppm elsewhere); positive means the measured m/z is heavier. Atoms per carbon atom allowed by MS-FINDER's common range: H 4, N 2, O 2.5, S 1, P 0.5, and with the boxes ticked F 1.5, Cl 1, Br 1, I 0.5, Si 0.5. A compound with more fluorine than 1.5 per carbon (for example a perfluorinated acid) is outside it.
- [M]+ is a permanent cation: the formula is the cation's own.

### Run notes

- When the local databases hold no candidate for the top formulas, the compound may simply not be in them: PubChem is not used.
- The extended element range is for PFAS and unusual compositions.
