# Theoretical

Give a compound or formula and an adduct. Get the precursor m/z and the isotope pattern.

## Results

### Theory

- **Precursor m/z**: m/z of the ion made of the most abundant isotope of each element (M), as the Mass Search page uses it.

### Isotope pattern

- Computed with IsoSpecPy from the ion formula (isotope fine structure), then merged to the chosen resolution: peaks closer than the FWHM (m/z ÷ resolution) become one peak (intensity-weighted m/z, summed abundance); Unit mass merges by nominal mass.
- M is the ion with the most abundant isotope of every element (<sup>11</sup>B, <sup>79</sup>Br, <sup>35</sup>Cl …). Lighter peaks are M−1 (e.g. <sup>10</sup>B).
- Peaks that stay apart at the chosen resolution carry (a), (b) after the label. Peaks below 0.1 % of the largest are dropped.

### Theoretical isotope pattern

- Hover a peak for its m/z and relative abundance.
