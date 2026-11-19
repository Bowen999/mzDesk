# Non-MS Detection

Give a compound. Get how well UV, fluorescence and universal detectors (CAD / ELSD / RI) suit it. Rule-based.

## Results

### Detectors

- yes = good, partial = possible, no = poor
- UV/PDA value: chromophore strength and estimated λmax. strong: aromatic or two or more conjugated π units; weak: one isolated C=O, C=C or amide; none: no π system
- FLD value: native fluorescence (likely / possible / unlikely)
- Rules and thresholds: mzdesk/platform_fit/non_ms_rules.py

### Chromophores

- blue: strong, amber: weak
- numbers = atoms (as in the SMILES, from 1)

### Run notes

- Fluorescence is a structural alert only: the quantum yield depends on solvent, pH and quenchers.
- λmax values are class rules, not predictions (±10–20 nm). Solvent and pH shift them: phenols, anilines, acids and pyridines move 20–50 nm with ionization. Concentration, path length and slit width decide the real sensitivity.
- The RI detection limit is in the µg range, 100–1000× worse than UV for a strong chromophore.
- CAD and ELSD need a volatile mobile phase: formic or acetic acid, ammonium formate or acetate, no phosphate.
- A possible native fluorescence should be checked with a standard.
- The ELSD response is non-linear and less sensitive than CAD.
