# Similarity

Give two MS2 spectra. Get 14 classic similarity scores, optional deep-learning scores and a mirror plot.

## Inputs

- **Query spectrum**: Paste a peak list or upload an MGF / MSP file. A header precursor m/z is used when that field is blank.
- **Ion mode**: Used only by MS2DeepScore / Spec2Vec; the classic scores do not need it.

## Results

### Mirror plot

- Up: query. Down: reference.
- Dark: matched peaks.

### Run notes

- The Spec2Vec score is N/A in negative mode: the pretrained model is GNPS AllPositive, i.e. positive-mode spectra only, so a negative-mode score would be meaningless.
