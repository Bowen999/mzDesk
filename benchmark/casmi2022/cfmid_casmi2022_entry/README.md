# CFM-ID team's official CASMI 2022 submission

The participant submission of the CFM-ID team (Fei Wang et al.) to CASMI 2022, kept here as an external comparator. **Nothing in this folder was run by us; CFM-ID was not run.**

## Files

| File | Content |
|---|---|
| `Priority_1-250.xlsx`, `Bonus_251-500.xlsx` | The submission as downloaded: one row per challenge (500 in total), columns Compound Number (= challenge id), File, RT [min], Precursor m/z (Da), Priority/Bonus, Molecular Formula, Adduct, InChIKey or SMILES, Compound Class. Every answer is an InChIKey. |
| `Abstract_Fei.txt` | The participant's method abstract. |

## Source

- CASMI 2022 results page: https://fiehnlab.ucdavis.edu/casmi/casmi-2022-results (Wayback Machine snapshot 20260729183139), which links the public Google Drive folder "8_Participant Abstracts & Submissions".
- Fei Wang subfolder: https://drive.google.com/drive/folders/1cO4nzUP3MuLW_hldvS9NKwxv7aJqsXm8
- Downloaded 2026-10-31.

## What it is (and is not)

- The participant's own end-to-end entry: its own candidate databases, its own adduct and formula decisions. Not a run on the candidate lists used for the other comparators (PubChem via SIRIUS), so it is **not directly comparable** to CSI:FingerID / JESTR / MIST.
- CASMI 2022 accepted one answer per challenge, so only **Top-1** exists. Top-3, Top-10 and "truth listed" are not applicable (shown as n/a in the workbook, not 0).
- Version: CFM-ID 4 (minor version not stated).
- Method (paraphrased from the abstract): [M+H]+ and [M-H]- spectra were predicted with CFM-ID 4 for the compounds of 12 merged databases (ChEBI, DrugBank, MassBankEU, HMDB, DSSTox, ECMDB, FooDB, KEGG, LMSD, NP-MRD, STOFFIDENT, YMDB). CFM-ID predicts QToF spectra only, so the predictions at 10, 20 and 40 eV were stacked into one spectrum and the stacked CASMI spectra were searched against them (peaks below 1 % intensity removed, 0.5 Da window around each peak). Candidates were selected by precursor mass (10 ppm) and ranked by Dice coefficient plus dot product. When nothing matched under [M+H]+ / [M-H]-, other adducts were tried: candidates from the mass computed for the adduct of choice, then matched against the predicted [M+H]+ / [M-H]- spectra.

## Official CASMI scores

Structure: Priority 23/249, Bonus 10/247 (the official evaluation excluded challenges 81, 282, 432 and 476).

Our workbook keeps all 500 challenges in the denominator, like every other sheet. The submission has no hit among those four challenges, so the counts are unchanged: 23/250 Priority, 10/250 Bonus, 33/500 overall (positive 26/304, negative 7/196).

## Use in the workbook

`mzdesk_ablation_20261030/write_xlsx.py` reads the two xlsx files here and writes the tabs `CFM-ID (CASMI entry)` (structure; hit = submitted InChIKey first block equals the official one) and `CFM-ID formula (CASMI entry)` (submitted neutral molecular formula equals the official one).

Challenges 81 and 432: the `File` column of the submission differs from the official file name of the challenge (A_M24_posPFP_02 vs _01; A_M16_posPFP_01 vs _02) and in 81 also the precursor m/z; the rows are matched by Compound Number as in every other sheet.
