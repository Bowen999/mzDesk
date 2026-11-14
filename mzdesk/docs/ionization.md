# Ionization

Give a compound and a mobile-phase pH. Get the better ESI polarity for LC-MS and whether EI works for GC-MS.

## Results

### Ionization

- Polarity: the better ESI mode (+, −, or both); poor = weak in both
- Simple rules at the chosen pH: pKa sites (MolGpKa; textbook group values without its weights), O/N counts, GC amenability and the NEIMS molecular ion
- Thresholds: mzdesk/platform_fit/ionization_rules.py

### Ionizable sites

- Coral: acid site. Teal: base site. Violet: both.
- Label: atom and pKa.
- acids above pKa 10 are left out

### Run notes

- The pKa of a site carries an error of about one unit.
- ESI+ and ESI- rest on the solution pKa. Gas-phase proton affinity and droplet chemistry can move a compound up or down, mostly for weak bases (pyridine-type N, carbonyls) and for acids at low pH.
- Response strength is qualitative (strong / moderate / weak, relative to other small molecules). It says nothing about absolute sensitivity, matrix effects or ion source settings.
- The M+• intensity is a neural-network prediction (NEIMS). The intensity of a single peak is far less certain than the peak pattern.
- For a molecule that needs derivatization, NEIMS predicts the underivatized molecule: the M+• and spectrum are not those of the TMS or MeOx-TMS derivative that is measured (masses: Derivatization page).
- A permanent positive charge makes ESI+ strong at any pH; an acid in the same molecule is its counterion.
- Molecules much larger than the drug-like training set of MolGpKa get less reliable pKa values.
- Without MolGpKa weights the pKa values are textbook values of the functional groups, not a prediction; the rules still run.
- The MolGpKa patterns do not cover N-oxides, so their pKa is a textbook value (not predicted).
- Where NEIMS weights are missing or the structure is outside its range, the M+• line is a textbook rule of thumb.
