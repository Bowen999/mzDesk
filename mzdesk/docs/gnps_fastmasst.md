# Spectrum in Which Set

Give an MS2 spectrum. Get the public datasets and GNPS reference-library entries that match it (GNPS2 fastMASST), and where the spectrum occurs: which sample types, organisms, body parts, diseases and domain-MASST categories its hit files come from.

## Inputs

- **Join Pan-ReDU metadata**: On by default. Without the Pan-ReDU metadata on your computer, the page runs without the join and says so in a warning. The local version adds the metadata with `mzdesk setup redu`.
- **Background**: Choose the background for your ion mode. The default includes positive, alternating and unknown polarity but excludes negative and GC/EI.

## Results

### Where it occurs

One group per property. Open a group to see every value with its **Hit files**, **Hit datasets** and **Share of hit files**.

- Pan-ReDU: sample type, organism, body part, disease.
- GNPS2 Explorer: the species the dataset reports, counted by datasets. A dataset that lists several species counts under each; entries such as "reference compound" or "Data from multiple bacteria" are not species and are left out.
- Domain MASST: microbe, plant, tissue, food and personal care product, from microbeMASST, plantMASST, tissueMASST, foodMASST and personalCareProductMASST.

A file counts once, however many of its spectra hit. The group line gives the three largest values. One line under the title says how many hit files (Explorer: datasets) each source covers, and which domain-MASST tables have none of them. Blank and QC files are not counted. Values the source marks as missing are left out.

### Pan-ReDU enrichment

Which Pan-ReDU values (sample type, sample subtype, organism, body part, disease, health status) are over-represented among your hit files compared with the background you chose. Within each property, values found in at least 2 datasets come first, ranked by dataset-weighted enrichment; single-dataset values and missing values come last.

- **Hit files**, **Hit datasets**: how many hit files and hit datasets have the value. **Background datasets**: how many datasets of the background have it.
- **File enrichment**: the share of hit files with the value divided by the share of background files with it.
- **Dataset enrichment**: the share of the background datasets with the value that are hit datasets, divided by the same share over all background datasets.
- **Dataset-weighted enrichment**: like file enrichment, but each dataset weighs 1, however many files it has, so a large dataset does not dominate.

The page shows two decimals, and ">1000" for 1000 or more; the TSV in the downloads keeps the full precision and every column. The page leaves out **Search** and **Hit set** when they have one value; the TSV always has them. A search with a diagnostic rule adds the hit sets "Verified" and "Diagnostic pass" (the TSV writes them as `verified` and `diag_pass`; the whole run is `all`, shown as "All"). In the command line summary and the batch columns `top_*`, the top values are the most frequent ones (hit datasets, then hit files), not the most enriched.

### Passing hits by dataset

One group per dataset, folded. The group line gives the hit count, the number of files, the best cosine, the main sample type, organism and body part, and, last, whether Pan-ReDU has the dataset or its files. The link opens the dataset page of its repository. The GNPS reference library has its own group. The page shows up to 200 datasets and 50 hits each; when it cuts, a note under the title says so, and the downloads hold all of them. Hits without a dataset are grouped as "(no dataset)".

### Top suspect names

- Suspect entries are propagated annotations, not reference standards.

### Run notes

- The query spectrum is sent to GNPS2. An MS/MS match alone does not establish compound identity.
- Pan-ReDU enrichment describes the selected background and the matched metadata. It is not a statistical significance test.

### Sources and limits

- **Pan-ReDU** (GNPS2): per-file sample metadata, read from a copy on your computer.
- **Domain MASST tables**: the file lists of [microbe_masst](https://github.com/robinschmid/microbe_masst) (MIT), shipped with mzDesk as a snapshot of 14 August 2026 (commit `5dd2a0f`). They are matched by dataset and file name, so mzDesk works offline for this part. They are not updated between releases.
- **GNPS2 Explorer**: the dataset title and species. mzDesk asks for the datasets shown on the page (at most 200 per run), only by accession, within about 20 seconds in all, and keeps the answers for 30 days (7 days for a dataset Explorer does not know). A dataset without an answer is listed by accession only; a warning says how many.

Coverage is partial. Pan-ReDU and the domain MASST tables do not cover every public dataset, and a file that is in none of them adds nothing to the counts. Read the numbers as "in what is annotated", not as a survey of all data. The coverage line says how many hit files each source covers.

See [Licenses](/docs/licenses) for the licenses and citations.
