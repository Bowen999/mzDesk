# Reference Spectra

Enter a compound. Get its recorded spectra from MassBank and GNPS. Add your own spectrum to score them.

## Inputs

- **Compound**: a name, ID, SMILES, InChI or drawing. For a salt, the largest part is used.
- **Spectrum types**: MS/MS, EI or both.
- **Ion mode (MS/MS)**: filters MS/MS only. EI spectra are always positive.
- **Libraries**: In the local version, a library that is not installed is downloaded the first time (2.3 GB; it asks first).
- **My spectrum (optional)**: every reference is scored against it.
- **Fragment tolerance (Da, MS/MS)**: used for scoring MS/MS. EI always uses 0.5 Da.
- **Plot record (accession, optional)**: the record shown under the table after the run. Empty = the first row.

## Results

### Summary

- **Ion types**: the 3 most common adducts.
- **Best match**: only with your spectrum.

### Reference spectra

- Same compound = same first 14 characters of the InChIKey. Stereoisomers are listed too; **Stereo** says if the full InChIKey also matches.
- Without your spectrum: sorted by library, type, ion mode and adduct, then most peaks. With it: highest entropy first.
- Adduct "—": the record does not give one.
- Entropy and Cosine are scored as in Library Search. Scores ≥ 0.7 are marked.
- At most 1000 rows.
- Click a row to show that spectrum below. The row shown is highlighted.

### Selected spectrum

- After a run: the first row, or the Plot record. A click on another row shows that one here, without a new run. The page scrolls to it when it is out of sight.
- **Scores**: only with your spectrum. Entropy, Cosine, matched peaks and the fragment tolerance used.
- **Plot**: the spectrum. With your spectrum: the mirror plot. Up: yours. Down: the reference. Dark: matched peaks.
- **Record**: the fields mzDesk keeps for each record: library, accession (links to the record on MassBank or GNPS), name, type, ion mode, adduct, precursor m/z, peaks, instrument, collision energy, stereo, InChIKey and SMILES. **Copy of**: a GNPS record that copies a MassBank one.
- **Peak list**: m/z and intensity of every peak (base peak = 100).
- **Formats**: the spectrum as text, to paste into another page or program.

### Run notes

- Records without an InChIKey cannot be found by structure (about 1 in 5 in GNPS). Library Search finds them by spectrum.
