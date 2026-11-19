# CCS

Give a compound, an adduct and an observed CCS. Get the reference CCS of the same compound and adduct, or a model prediction if there is none, and whether your value agrees.

## Results

### Summary

- Observed CCS against the N₂ reference CCS of the same compound and adduct.
- **Reference CCS**: The reference value used: the one closest to your observed value, or the median of the measured values when there is nothing to compare. If the index has no value for this compound and adduct, the model prediction.
- **Source**: Measured: a published N₂ value from the CCSbase source data (CCS type and dataset in the note). Predicted: a model, weaker than a measurement.
- **Observed**: Your value(s). With several, the one closest to each reference is used.
- **Δ**: (observed − reference) / reference × 100, with the tolerance for this platform pair.
- **Verdict**: Consistent: within the tolerance. Inconsistent: outside it; this never excludes the compound. Not comparable: He, TWIMS calibrated with He, or FAIMS / DMS. No observed value, or no reference value for this compound and adduct.

### Reference values

- Every measured N₂ value in the index for this compound and adduct. Δ is taken against the observed value closest to that reference.
- The tolerance depends on the platform of your measurement and the CCS type of the row. The ? on a row gives the publication and the name in the source.

### Other ions of this compound

- The same compound as other adducts. Information only: a CCS depends on the adduct, so these are not compared.

### Prediction

- Only when the index has no value for this compound and adduct. A model prediction, weaker than a measured value; the ? gives the model and its test error.

### Assumptions

- Comparison rules, tolerances with their sources, the prediction model, and the data source with its licence and citation.

### Run notes

- Comparison rule: the same compound (first 14 characters of the InChIKey after taking the largest fragment and neutralising; stereo is ignored) and the same adduct. Other adducts of the compound are listed for information only.
- When several measured values exist, the summary uses the one closest to your value (relative to its tolerance): one agreeing measurement is enough. Every value is listed in the table. With several observed values (for example two conformers) the one closest to each reference is used.
- N₂ scale only: all reference values were measured in N₂ (DTCCS N2; TWCCS N2, calibrated with N₂ values; TIMSCCS N2). A value measured in He, or a TWIMS value calibrated with He reference values (TWCCS N₂→He), is reported as not comparable.
- FAIMS / DMS gives a compensation voltage, not a CCS, so nothing is compared. For TIMS, enter the CCS your vendor software reports (Å²), not 1/K0.
- Tolerance basis: DTIMS interlaboratory RSD 0.29 % and single-field bias 0.54 % (Stow et al., Anal. Chem. 2017); DT vs TIMS within ±1 % and DT vs TW within ±2 % for 95 % of ions (Feuerstein et al., J. Am. Soc. Mass Spectrom. 2022, steroids). TIMS vs TW is derived (the two biases against DT added up), not measured.
- A value outside the tolerance is flagged inconsistent and is never used to exclude a compound: in the same study fewer than 1.5 % of ions differed by up to 7 % between platforms (conformers, protomers).
- Prediction: only when the index has no value for the compound and adduct, and only for adducts with at least 100 training rows. The model uses m/z, adduct and 42 MQN descriptors and was chosen against HistGradientBoosting on a validation set. The tolerance is the 95th percentile of the test error, rounded up to 0.1 %. A prediction is weaker than a measured value.
- A compound whose nearest training compound has a Tanimoto (Morgan radius 2, 2048 bits) below 0.5 is flagged as outside the training domain. In the test set such compounds had a larger 95th percentile error than the rest; treat their predicted CCS with extra care.
- A row marked 'Check this row' lists an m/z that differs from the m/z calculated for the structure and adduct. The source may list average masses, or its name, structure or adduct may not match.
- Data: CCSbase source data (dylanhross/c3sdb, MIT licence), values from 26 published datasets (DT, TW, TIMS), all measured in N₂. Structures come from the source or were assigned to the compound names by the c3sdb authors (lipid structures can be generic); values without a structure were skipped.
- Citation: D. H. Ross, J. J. Cho, L. Xu. Breaking Down Structural Diversity for Comprehensive Prediction of Ion-Neutral Collision Cross Sections. Anal. Chem. 92, 4548-4557 (2020). https://doi.org/10.1021/acs.analchem.9b05772
