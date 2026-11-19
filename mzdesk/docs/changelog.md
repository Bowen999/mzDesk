# Changelog

This page lists what changed in each version of mzDesk, newest first.

## v1.2.0 (2026-11-17)

- [ESI MS/MS Annotation](/docs/esi_ms2_annotation): new **Force merge** under Advanced. MS-FINDER and FRIGID also run after a library hit, and all candidates are ranked together. Ticked in the local version, unticked on mzdesk.com, because a run takes longer.
- One score instead of Base and Final: **Score** is the evidence score, lowered by up to 10 % when the structure is unlikely to be detected. Candidates are ordered by it, so rank 1 has the highest score (an observed CCS can still move a candidate). The evidence score is in the detail line.
- On mzdesk.com and in the local version (`mzdesk update`).

## v1.1.5 (2026-11-16)

- [ESI MS/MS Annotation](/docs/esi_ms2_annotation): the ranked-candidates table shows the final score only for candidates within 10 % of the best base score, where detectability set the order. The others show "—" and their base score is in the detail line. The ranking is unchanged.
- On mzdesk.com now. The local version gets it in a later update.

## v1.1.4 (2026-11-14)

- [Spectrum in Which Set](/docs/gnps_fastmasst): the blocks now run Searches, Where it occurs, Pan-ReDU enrichment, Passing hits by dataset, names, Diagnostic verification, and the downloads last.
- **Pan-ReDU enrichment**: within each property (sample type, sample subtype, organism, body part, disease, health status), values found in at least 2 datasets come first, ranked by dataset-weighted enrichment. Single-dataset and missing values come last. Values show 2 decimals; the TSV keeps full precision. A new Background datasets column and column help were added, and columns that are the same in every row are hidden.
- **Passing hits by dataset**: a cap on the hits shown is now a note under the title, not a warning. Each group line reads counts, best cosine, sample type, organism, body part, then status. tissueMASST categories read like "Mus musculus · caecum · healthy · male".
- All tables: numbers are right-aligned, short labels stay on one line, group lines sit under the title, empty cells and group lines are easier to read, and rows are easier to tap on a phone.
- On mzdesk.com now. The local version gets it in a later update.

## v1.1.3 (2026-11-14)

- [Spectrum in Which Set](/docs/gnps_fastmasst): the passing hits now come in one folded group per dataset, with the dataset title, a link to its page, the hit count, the number of files, the best cosine and the main sample type, organism and body part. The flat Top datasets table is gone.
- New block **Where it occurs**: how many hit files and datasets fall under each sample type, organism, body part and disease (Pan-ReDU), each species (GNPS2 Explorer), and each domain-MASST category (microbe, plant, tissue, food, personal care product). Two new downloads hold the dataset list and these counts.
- **Join Pan-ReDU metadata** is now on by default. Without the Pan-ReDU metadata on your computer, the page warns and still runs.
- On mzdesk.com now. The local version gets it in a later update.

## v1.1.2 (2026-11-14)

- [Reference Spectra](/docs/ms2_reference): click a row of the table to show that spectrum below it, without a new run. You see the record, the plot (the mirror plot and the scores when you gave your spectrum), the peak list and the spectrum as text.
- On mzdesk.com now. The local version gets it in a later update.

## v1.1.1 (2026-11-14)

- Fixed: a hard scroll past the top or bottom of a page no longer drags the sidebar along and shows a blank strip. The sidebar and the Docs contents list now scroll on their own.
- On mzdesk.com now. The local version gets it in a later update.

## v1.1.0 (2026-11-14)

- New page: [Reference Spectra](/docs/ms2_reference), under Spectrum. Enter a compound and get its recorded MassBank and GNPS spectra. Add your own spectrum to score them and see a mirror plot.
- On mzdesk.com now. The local version gets it in a later update.

## v1.0.0 (2026-11-14)

- First public release. mzDesk brings the everyday tools of small-molecule mass spectrometry into one place: compound identity, MS1 formula and database search, MS2 annotation, spectral library search and similarity, structure elucidation, separation, ionization and adduct prediction, GNPS searches, and more. See [Workflows and Modules](/docs/workflows-and-modules).
- Use it on mzdesk.com, or install the local version on macOS (Apple silicon and Intel), Linux or Windows. See [Local Version](/docs/local-version).
