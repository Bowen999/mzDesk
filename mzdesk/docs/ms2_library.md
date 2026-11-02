# Library Search

Give an MS/MS or EI spectrum. Get ranked matches from MassBank and GNPS, scored with the similarity metrics you choose.

## Inputs

- **Precursor m/z (MS/MS only)**: Required for MS/MS unless the spectrum header carries it. Auto = any adduct; a chosen adduct must also match the library record.
- **Fragment tolerance (Da, MS/MS only)**: Peak-matching tolerance for MS/MS. EI always uses unit mass (0.5 Da) and hides this field.
- **Libraries**: In the local version, a library that is not installed is downloaded the first time (2.3 GB; it asks first).
- **Metrics**: Hits are ranked by the first selected metric in this list order, not by the order you clicked.

## Results

### Summary

- **Best hit**: The first row of the table: the hit with the highest score on the ranking metric.

### Library hits

- A spectrum is a hit when ANY selected metric is at or above the threshold (N/A scores are ignored). The table shows the best 100.
- Hits are ranked by the first selected metric that applies, in the order of the Metrics list (not the order you clicked them); a metric that is N/A for every hit is skipped. Ties: more matched peaks first.
- Matched peaks: peak pairs matched one-to-one within the fragment tolerance (greedy, as in CosineGreedy). Δ ppm: query minus library precursor.
- MS/MS: the precursor m/z is required (a header in an MGF / MSP spectrum counts). Candidates are the library MS/MS spectra within the precursor tolerance, in the same ion mode (and with the chosen adduct, unless Auto); every candidate is scored with every selected metric. With both libraries selected, GNPS copies of MassBank spectra are hidden when the MassBank original is also a candidate.
- EI: the precursor is not used. Fragments are matched at unit mass (0.5 Da). Candidates = top 5000 by entropy similarity; a spectrum with low entropy but high other scores can be missed.
- Deep-learning scores need torch / gensim; they are N/A where those are missing.

### Mirror: query vs top hit

- Up: query. Down: the top-ranked hit.
- Dark: matched peaks.

### Run notes

- MS2DeepScore and Spec2Vec are N/A for a negative-mode query: the pretrained Spec2Vec model is GNPS AllPositive, i.e. positive-mode spectra only, so a negative-mode score would be meaningless.
