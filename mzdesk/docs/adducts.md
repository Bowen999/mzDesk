# Adducts

Which ions appear (adducts, in-source losses, dimers) for the chosen polarity and mobile phase additive, and which is strongest. Rule-based: pKa sites (MolGpKa) and lipid class.

## Results

### Ions

- Ranked by expected intensity: a rule score (high / medium / low), not a measured intensity.
- m/z computed as on the Theoretical m/z page.

### Run notes

- Solution ionization only fine-tunes the scores (× 0.8–1): in ESI, acids still give [M-H]- from acidic eluents and bases give [M+H]+ from neutral ones.
- Lipid classes use class conventions (LIPID MAPS / MS-DIAL / LipidSearch practice) instead of the pKa rules.
- m/z and ion formulas use the same code as the Theoretical m/z page. That page uses the formula as drawn, so for a structure drawn as an ion enter the neutral structure there to get the same m/z.
- Scores are rule-based ranks (0–100), not measured intensities. Matrix suppression, source settings and concentration are not modelled.
- Dimers ([2M+X]) are concentration dependent.
- Without an additive the mobile phase is unbuffered, so the real pH is uncertain.
