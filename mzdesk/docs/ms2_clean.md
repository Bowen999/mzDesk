# Spectrum Cleanup

Give an MS2 spectrum in any common format. Get it cleaned, plotted and written in every common format (MS-DIAL, GNPS, MassBank, HMDB, MSP, MGF …).

## Inputs

- **MS2 spectrum (any common format)**: Paste any common format, upload a file, or resolve a USI. The precursor m/z and name in an MGF / MSP header are used when the fields are blank.
- **Precursor m/z (optional)**: Blank = read from the MGF / MSP header when present.

## Results

### Cleaned spectrum

- Grey: removed peaks.
- Hover a peak for details.

### Peaks

- Rel. intensity: % of the input base peak. Reason: noise, precursor

### Run notes

- The axis and Intensity use the cleanup output scale. When possible, normalization is based on kept peaks; removed peaks use the same scale and can exceed the largest kept peak.
