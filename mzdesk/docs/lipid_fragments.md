# Lipid Fragment Explanation

Give a lipid, as a short name or a structure, and its MS2 spectrum. Get the LIPID+ and LIPID MAPS names, an illustrative structure and the peaks that structure explains (MetFrag).

## Inputs

- **Lipid**: a short name such as `PC 16:0/18:1`, a SMILES or a drawing. Only a short name gives the names.
- **Precursor m/z (optional)**: Blank = theoretical m/z. **Auto** picks the adduct whose m/z fits the precursor within 10 ppm or 0.01 Da; with no precursor, or none that fits, it uses [M+H]+ or [M-H]-.
- **MS2 spectrum**: the measured spectrum. A precursor in the file header is used when the field is blank.
- **Tree depth (bond cleavages)**: default 3. Higher is slower and can over-explain.

## Results

### LIPID+ and LIPID MAPS

- The same names as [Lipid Short Name](/docs/lipid_short_name). A SMILES or a drawing gives — here.

### Illustrative structure

- One structure that fits the name. The shaded parts are assumed: the C=C positions, group positions or chain split the name does not give.
- A name at the sum level (PC 34:1) is split into chains by the default rule. The spectrum does not choose the split.
- A class without a template (free sterols, sulfatides and others) has no structure, so the fragments are not explained. Enter the structure as SMILES or draw it.

### Fragment explanation

- The same results as [Fragment Explanation](/docs/ms2_fragments), on the structure above.

### Run notes

- The explanation is only as good as the structure. A fragment that depends on an assumed C=C position or chain split is not evidence for that position or split.
