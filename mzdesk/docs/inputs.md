# Inputs and Formats

This page lists what the input fields accept.

## Compounds

A compound field has three modes:

| Mode | You enter |
|---|---|
| Name / ID | A name, CAS number, InChIKey, PubChem CID, or an ID from HMDB, KEGG, ChEBI, LIPID MAPS, DrugBank, CompTox or ChEMBL. |
| Draw | A structure drawn in the editor. |
| SMILES / InChI | A SMILES, or an InChI that starts with `InChI=`. |

- A name or ID is looked up in PubChem and UniChem, so it needs the internet. SMILES, InChI and formulas are read on your computer.
- Two pages also accept a formula such as `C8H10N4O2` in the first mode: Theoretical and Alternative Explanations. On the pages that need a structure, a formula is refused.
- A structure with isotope labels keeps them. A page whose values ignore the labels says so.

## Spectra

A spectrum field has three sources:

- **Paste**: paste text. The box shows how many peaks were read.
- **MGF / MSP file**: choose a file (`.mgf`, `.msp`, `.txt`, `.csv`, `.json`, `.xml`, `.ms`). Only the first spectrum is used, and the page warns if there are more.
- **USI**: a Universal Spectrum Identifier such as `mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB…`. The spectrum is fetched from GNPS / MassIVE, so it needs the internet.

mzDesk detects the format by itself. These are read:

| Format | Looks like |
|---|---|
| Peak list | One `m/z intensity` pair per line, separated by space, tab, comma or `;`. A lone m/z gets intensity 1. |
| MS-DIAL / MoNA | `m/z:intensity` tokens, several per line. |
| MGF | `BEGIN IONS` to `END IONS`; the precursor comes from `PEPMASS`. |
| MSP | `Name:` and `Num Peaks:` blocks. |
| MassBank record | The lines under `PK$PEAK:`. |
| SIRIUS .ms | The `>ms2peaks` block; the precursor comes from `>parentmass`. |
| HMDB XML | `<ms-ms-peak>` elements. |
| GNPS peaks_json | `[[m/z, intensity], ...]`. |
| JSON | A list of `{mz, intensity}` objects, or a record with `peaks`. |
| Two arrays | `mz = [...]` and `intensities = [...]`. |
| CSV | A header such as `mz,intensity`. |

- Lines that do not start with a number are headers and are skipped.
- Intensities can be on any scale.
- Peaks with m/z or intensity of 0 or less are dropped. Peaks with the same m/z (to 5 decimals) are merged, keeping the stronger.
- Decimal commas (`138,0662`) and thousands separators (`1,200`) are read, and the page tells you.
- The workflows accept at most 5000 peaks per spectrum.
- [Spectrum Cleanup](/docs/ms2_clean) writes a spectrum in 11 formats.

## Precursor, ion mode and adduct

- **Precursor m/z**: a number. Leave it blank to use the m/z in the spectrum header (MGF, MSP, MassBank, SIRIUS or USI).
- **Ion mode**: Positive or Negative.
- **Adduct**: **Auto** or a choice such as `[M+H]+`, `[M+Na]+`, `[M-H]-`, `[M+HCOO]-`. The list depends on the engine. What Auto does is written next to it: MS/MS Annotation ranks the adducts, Library Search accepts any adduct, Structure Prediction uses `[M+H]+` or `[M-H]-`.

## MS1 isotope pattern

Optional on some pages. One `m/z intensity` pair per line, around the precursor.

## Tables and files

- **Alkane series** (EI MS Annotation, Calculator): one row per line, the carbon number and the RT, for example `C8 8.00`. At least two rows.
- **Fragment Pattern**: one `.mzml`, `.mzxml` or `.mgf` file of at most 20 MB, and a MassQL query.
- **Candidate structures** (EI MS Annotation): one SMILES or name per line. Names are looked up online.

## Tolerances

Each tolerance field names its unit, ppm or Da. A request to the server is limited to 8 MB.
