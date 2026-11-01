# EI MS Annotation

Enter one EI spectrum (usually from GC-MS). Get candidate metabolites, named as the original underivatized compounds, with evidence.

## Inputs

- **Derivatization**: Auto reads the marker ions: m/z 73 backed by 75 or 147 for TMS (m/z 147 alone when the scan starts above m/z 73); M-31 (MeOx) and M-57 (TBDMS) need the molecular ion.
- **Molecular ion m/z**: Optional. The M+• of the derivative when derivatized.
- **Candidate structures**: Optional.
- **Sample RT**: Same unit as the alkane RTs. Give it together with the alkane series.

## Results

### Run notes

- The Kovats RI is reported only: the libraries have no retention index, so it does not change the ranking.
- A compound with an experimental library spectrum ranks above one with only a NEIMS-predicted spectrum; inside a group the score decides.
- EI library search: all library EI spectra are ranked by weighted entropy similarity and the top 5,000 are scored with the selected metrics. A spectrum with low entropy but high other scores can be missed.
