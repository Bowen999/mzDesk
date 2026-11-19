# Formula Annotation

Give a precursor m/z. Add an MS1 isotope pattern and an MS2 spectrum to sharpen the ranking. Get ranked molecular formulas (BUDDY).

## Inputs

- **Precursor m/z**: Required unless the MS2 spectrum header carries it. BUDDY supports m/z 1–1500.
- **MS1 isotope pattern (m/z intensity per line, optional)**: Optional isotope pattern around the precursor (m/z intensity per line). A pattern that misses the precursor is ignored with a warning.
- **Formula (optional: skip ranking, explain the MS2 with it)**: Optional: if given, ranking is skipped and the MS2 is explained with this formula. It must be a neutral formula (no charge); clear it to use MS1 and rank candidates.
- **Include halogens (F, Cl, Br, I)**: Only used when Formula is blank; it widens BUDDY's element search.

## Results

### MS2 fragment formulas · top candidate

- Hover a peak for its formula, Δ ppm and neutral loss.
