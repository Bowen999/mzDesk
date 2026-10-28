# SIRIUS settings for the blind CASMI 2022 run

SIRIUS 6.3.5 (lib 5.8.0), CSI:FingerID lib 3.0.14, web service via academic login.
Script: `scripts/02_run_sirius.sh [ids...]` (all settings are in the script, in the `CFG` and `CHAIN` arrays).
Same settings for every compound. No adduct, formula or structure is given per compound.

## Input (`01_extract_ms.py`)

| Item | Value | Why |
|---|---|---|
| `>parentmass` | measured MS1 monoisotopic m/z (`--parentmass ms1`, default) | CASMI-given m/z is rounded and up to 8 ppm off; the measured value is data, not an answer. Over the 500 files the measured value is a median 1.3 ppm (max 8.1) from the given one, mean -1.5 ppm. Given value stays in `ms/extraction_report.tsv` (`precursor_given`; `parentmass_used`, `parentmass_source` added). `--parentmass given` restores the old behaviour. |
| `>ionization [M+?]+/-` | polarity only | adduct is left to SIRIUS |
| MS2 | 3 `>collision 35/45/65` blocks, MS1 isotope block | checked with `mgf-export --write-ms1`: one compound, 3 MS2 + MS1 parsed, PEPMASS = measured m/z |

## SIRIUS config (`config ...`)

| Setting | Value | Why |
|---|---|---|
| `AlgorithmProfile` | `orbitrap` | Thermo Q Exactive HF |
| `MS1MassDeviation.allowedMassDeviation` | 10 ppm | parentmass is measured, but up to ~8 ppm off the theoretical mass at the tail; 10 ppm is the SIRIUS default and safe |
| `MS2MassDeviation.allowedMassDeviation` | 10 ppm | MS2 at 15k resolution, fragments are low intensity |
| `AdductSettings.enforced` | pos `[M+H]+,[M+Na]+,[M+NH4]+,[M+K]+,[M+H-H2O]+`; neg `[M-H]-,[M+Cl]-,[M+FA-H]-,[M+HAc-H]-,[M-H2O-H]-` (both lists in one string) | `enforced` = always considered, no MS1 evidence needed. SIRIUS only uses entries of the polarity in the `.ms` file. Verified on the pilot: formate and NH4 adducts appear in the ranked list. SIRIUS spells them `[M + CH2O2 - H]-` (formate), `[M + C2H4O2 - H]-` (acetate), `[M + H3N + H]+` (NH4+), `[M - H2O + H]+` |
| `AdductSettings.fallback` | same list | used if no adduct is detected; kept identical so nothing changes |
| `AdductSettings.detectable` | default | not needed, all wanted adducts are enforced |
| `FormulaSettings.enforced` | `C,H,N,O,P,F,I` | default is CHNOP only. F and I have no isotope signature so they cannot be auto-detected; without them every fluorinated or iodinated compound would be unsolvable (the first pilot compound has an F). This equals SIRIUS `--elements-extended-organic`. Cost: more candidate formulas, slightly slower |
| `FormulaSettings.detectable` | `S,Br,Cl,B,Se,Fe,Zn,Mg` (default) | added only when the MS1 isotope pattern indicates them |
| `FormulaSettings.fallback` | `S` (default) | used when there is no usable isotope pattern |
| `StructureSearchDB` and `structures -d` | `PUBCHEM` | user decision. (SIRIUS 6 default is BIO.) Both the config key and the subtool flag are set |
| Formula candidates (`formulas -c`, `NumberOfCandidates`) | 10 | default. `NumberOfCandidatesPerIonization=1` (default) keeps at least one formula per adduct |
| Formula threshold before fingerprint/structure search | default (soft threshold on) | default. Only the formula candidates passing the threshold are sent to CSI:FingerID. Not changed, so some low-ranked formulas have no structures |
| `NumberOfStructureCandidates` | 10000 (default) | cap per compound; exported in full |
| MSNovelist (`denovo-structures -c`, `NumberOfMsNovelistCandidates`) | 128 | default |
| Heuristic tree computation | default (>300 m/z heuristic ranking, >650 m/z heuristic only) | default; large compounds use heuristic trees, which is much faster |
| `Timeout.secondsPerInstance` | 300 s | default. A compound that times out is retried in the completion pass (and every later re-run) |
| ZODIAC | not run | single compounds, no network; adds nothing for the formula/structure ranking |
| CANOPUS (`classes`) | runs but its output is not used | the CLI chain `formulas fingerprints classes structures` only accepts `structures` after `classes`; cost is small |
| Spectral library search | runs by default (`SpectralSearchDB=ALL`), found no matches in the pilot | not turned off; SIRIUS libraries are not the CASMI reference |
| `--cores` | 8 | machine has 10 logical CPUs (5 physical) |

## Tool chain

`sirius --cores 8 -i <ms dir> -o project.sirius config <CFG> formulas -c 10 fingerprints classes structures -d PUBCHEM denovo-structures -c 128`

Then `summaries --full-summary --top-hit-summary --top-hit-adduct-summary --format=TSV --output <dir>`.

## Resume / re-run

* The `.sirius` project is kept. SIRIUS does not dedupe imports: importing the same `.ms` into an existing project creates a duplicate compound. The script therefore keeps a ledger `imported_ids.txt` next to the project and imports each ID once (in chunks of `CHUNK`, default 25).
* After the chunks, a completion pass `-i project -o project` runs the same chain; SIRIUS skips results that already exist (verified: 8 finished compounds, 0.3 s of compute). It retries compounds that timed out or failed on the web service.
* Do not change the settings on an existing project (SIRIUS recomputes and drops downstream results). To start over, delete the project dir.
* The export is rewritten each time (written to `*.tmp`, then swapped in).

## Exported files (`results/.../sirius_summaries/`)

Challenge ID is in `mappingFeatureId`, e.g. `903346292385728429_challenge_003_challenge_003` (regex `_(challenge_\d+)_challenge`).

| File | Content |
|---|---|
| `formula_identifications_all.tsv` | ranked formula + adduct (`formulaRank, molecularFormula, adduct, precursorFormula, SiriusScore, SiriusScoreNormalized, TreeScore, IsotopeScore, numExplainedPeaks, ...`). Several adducts of the same ion can tie with identical scores |
| `formula_identifications.tsv`, `..._adducts.tsv` | rank 1 per compound; best formula per adduct |
| `structure_identifications_all.tsv` | all PubChem candidates over all formulas: `structurePerIdRank, formulaRank, ConfidenceScoreExact, ConfidenceScoreApproximate, CSI:FingerIDScore, molecularFormula, adduct, InChIkey2D, InChI, smiles, name, pubchemids, ...` (confidence is per compound, repeated on each row; `InChIkey2D` is the 14-character first block only, full InChI is given) |
| `structure_identifications.tsv` | rank 1 per compound |
| `denovo_structure_identifications_all.tsv` | all MSNovelist candidates: `structurePerIdRank, CSI:FingerIDScore, ModelScore, molecularFormula, adduct, InChIkey2D, InChI, smiles, ...` (no confidence score for MSNovelist) |
| `qc_results_present.tsv` | per requested ID: has formula / PubChem structure / MSNovelist result (written by the script) |

## Pilot findings (8 compounds: 1 2 3 4 18 253 300 450)

* Compounds without structure/MSNovelist output exist by design (SIRIUS drops formulas whose fragmentation tree has < 2 fragments). In the pilot, challenge_450 (m/z 1335) has no MSNovelist candidates; in 5 of 24 other test compounds MSNovelist was empty and 3 had no PubChem structures.
* Web-service errors: none that stopped a run. A harmless `InterruptedIOException` is logged at shutdown. One compound hit the 300 s timeout in a 24-compound batch and was completed by the completion pass.
