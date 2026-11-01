# Fragment Explanation

Which peaks the candidate structure explains, and with which fragment (MetFrag).

## Inputs

- **Precursor m/z (optional; blank = theoretical)**: Blank = theoretical m/z from the structure and adduct. A measured value is used for neutral losses and a mismatch warning.
- **MS2 spectrum**: The measured MS2 spectrum; MetFrag annotates its peaks.
- **Tree depth (bond cleavages)**: Bond cleavages considered by MetFrag; higher is slower and can over-explain.

## Results

### Observed spectrum

- Hover a peak: fragment formula, Δ ppm, neutral loss and fragment structure.

### Run notes

- Where the same fragment fits several places (symmetric or repeated units), the drawing marks the one whose hydrogen counts agree best; the card says how many other places fit equally well.
- A fragment is not placed on the structure when it has no fit or its atoms do not match the ion formula; the fragment's own structure is drawn instead.
