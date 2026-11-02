# De Novo Structure

Give a precursor m/z and an MS2 spectrum. Get generated candidate structures from FRIGID, for compounds that are in no database. They are hypotheses.

## Inputs

- **Precursor m/z (positive mode only)**: FRIGID only generates positive-mode structures; the negative option is hidden.
- **MS1 isotope pattern (m/z intensity per line, optional)**: Optional isotope pattern. It is used only when Formula is blank, to improve BUDDY's top-1 formula.
- **Formula (optional; blank = BUDDY's top-1 formula)**: Blank = BUDDY's top-1 formula (not available for [M]+ or [M+H-2H2O]+; enter a formula there). FRIGID only knows C H N O P S F Cl Br I B Si Se As Fe Co Na K.
- **Mode**: Quick = 16 samples (a first look). Full = up to 100 samples and takes longer.
- **Include halogens (BUDDY formula fallback)**: Used only when Formula is blank: widen BUDDY's formula search to F, Cl, Br and I.

## Results

### Generated candidates

- These are hypotheses, not identifications. FRIGID does not search a database: it writes structures for the formula you give and the molecular fingerprint it predicts from the MS/MS, so a candidate can be a molecule that exists nowhere. Check any candidate against a reference standard, a library spectrum (Library Search) or a predicted spectrum (MS/MS Prediction).
- How often the right structure is among them (the authors' benchmark of FRIGID-base, revised 2026-10-21, with the true formula given): MassSpecGym top-1 8.21 %, top-10 9.82 %; NPLIB1 top-1 13.31 %, top-10 20.15 %. Quick mode draws 16 samples, so it is a first look; Full mode draws until 10 different formula-matched structures or 100 samples.
- Score: Tanimoto similarity (Morgan fingerprint, radius 2, 4096 bits) between the candidate and the fingerprint FRIGID predicted from the spectrum. It is a similarity, not a probability, and scores of different spectra are not comparable.
- Only candidates with the right formula are shown; isotope-labelled, radical, wrongly charged and duplicate structures are hidden. SMILES carry no stereochemistry (FRIGID does not predict it), so the InChIKey is the 2D block only.
- A familiar compound can come back because the model saw it in training, not because it read this spectrum.
- Formula: if you leave it blank the page takes the top-1 formula of BUDDY (Orbitrap settings; halogens are off unless you tick Include halogens) and says so under Formula used; a wrong formula gives wrong structures. FRIGID knows C H N O P S F Cl Br I B Si Se As Fe Co Na K, and only positive-mode adducts.
- In the local version, FRIGID is downloaded the first time (about 3 GB; it always asks first). On an Apple-silicon Mac a run takes roughly 20-40 s and 4-5 GB of memory. On mzdesk.com a Quick run takes about two minutes. One heavy job (this page, Structure Prediction, or the MS-FINDER and FRIGID steps of MS/MS Annotation) runs at a time.

### Run notes

- If the formula is not entered and the BUDDY formula is wrong, every structure is wrong too. Enter the formula to override it.
- When no usable structure is generated, the compound may be outside what FRIGID can generate.
- Device auto picks the Apple GPU (MPS) only when the Mac has at least 12 GB of memory available and swap to spare; otherwise the CPU. MZDESK_FRIGID_DEVICE = auto | mps | cpu overrides it.
- The generated structures are hypotheses, not identifications.
