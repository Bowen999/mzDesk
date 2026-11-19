# Targeted Method

Give one compound. Get a method sheet: can it be measured, platform, column, pH, polarity, main ion, derivatization and published methods.

## Inputs

- **Compound**: Name, database ID, SMILES or InChI.
- **pH**: Mobile-phase pH. Blank = 2.7 (0.1 % formic acid).
- **Mobile-phase additive**: Formic acid, ammonium formate, ammonium acetate, acetic acid or none.
- **Instrument**: Any = the best platform overall. LC-MS or GC-MS = judge only that platform.

## Results

### Run notes

- The method sheet is rule-based and built from the nine steps.
- A structure without carbon is outside the rules, which are written for organic compounds: the sheet is not reliable for it. Inorganic ions are usually measured by ion chromatography or ICP-MS.
