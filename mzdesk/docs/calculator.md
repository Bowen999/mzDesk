# Calculator

Give masses, concentrations or retention times. Get the ppm error, labelled masses, a concentration conversion, a dilution, a collision-energy conversion and a Kovats RI. Each block has its own inputs.

## Results

### Run notes

- Isotope labelling assumes complete incorporation of the label. The table shows the natural isotope distribution of the remaining atoms, normalised to the most abundant labelled isotopologue.
- Crosstalk is the unlabelled analyte's natural M+n abundance relative to its M peak, where n is the label's nominal mass shift.
- Mass ↔ molar concentration conversion needs a molar mass; 1 µg/mL = 1 mg/L. Dilution uses the converted stock concentration and C1V1 = C2V2.
- NCE ↔ eV uses CE = NCE × precursor m/z × charge / 500. Instrument vendors define NCE differently, so this is a generic conversion.
- Kovats RI is linearly interpolated between adjacent n-alkanes. Values outside the measured range are extrapolated and flagged.
