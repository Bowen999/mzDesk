# Workflows and Modules

A module does one job on one page. A workflow runs several steps in order and ends with one answer.

## Workflows

Each workflow shows its steps on the page. A step can be skipped, and a step that fails does not stop the run.

### MS/MS Annotation

![The steps of MS/MS Annotation](fig:workflow-esi)

- **You give**: an ESI MS/MS spectrum and its precursor m/z. Optional: an MS1 isotope pattern, an observed CCS.
- **You get**: the adduct, the formula and ranked structure candidates, each with a confidence level.

The library search comes first. MS-FINDER runs only if no library hit passes the Level-2 gate. FRIGID runs only if MS-FINDER is not good enough and **Run FRIGID** is ticked. [Notes](/docs/esi_ms2_annotation)

### EI MS Annotation

![The steps of EI MS Annotation](fig:workflow-ei)

- **You give**: an EI spectrum, usually from GC-MS. Optional: the derivatization, the molecular ion, candidate structures, an alkane series with the sample RT.
- **You get**: candidate metabolites, named as the original compounds without the derivative, with a level and a Kovats RI.

[Notes](/docs/ei_ms_annotation)

### Targeted Method

![The steps of Targeted Method](fig:workflow-targeted)

- **You give**: one compound. Optional: pH, mobile-phase additive, instrument (Any, LC-MS, GC-MS).
- **You get**: a method sheet with platform, column, mobile phase and pH, polarity, main ion, theoretical m/z, derivatization, non-MS alternative and published methods.

[Notes](/docs/targeted_method)

## Modules

### Compound

- [Identity](/docs/identity): names, IDs, classification.
- [Structure & Formats](/docs/structure_formats): up to 3 structures in 2D and 3D, as SMILES, InChI, InChIKey, MOL.
- [Properties & Groups](/docs/properties): formula, masses, descriptors, functional groups.
- [Lipid Short Name](/docs/lipid_short_name): lipid shorthand to LIPID+ and LIPID MAPS names.

### Detectability

- [Platform](/docs/platform_fit): lipidomics, metabolomics, CIL LC-MS or GC-MS.
- [Separation](/docs/separation): charge, logD and separation modes at a pH.
- [Ionization](/docs/ionization): ESI polarity and EI usability.
- [Adducts](/docs/adducts): expected ions, strongest first.
- [Non-MS Detection](/docs/non_ms): UV, fluorescence, CAD / ELSD / RI.
- [Derivatization](/docs/derivatization): reactive groups, derivative formula and mass.
- [Known Methods](/docs/known_methods): published methods, or those of similar compounds.

### MS1

- [Mass Search](/docs/ms1_search): precursor m/z to database compounds.
- [Theoretical](/docs/ms1_theory): precursor m/z and isotope pattern.

### Spectrum

- [Spectrum Cleanup](/docs/ms2_clean): clean a spectrum, write it in 11 formats.
- [Similarity](/docs/ms2_similarity): scores and a mirror plot for two spectra.
- [Library Search](/docs/ms2_library): MS/MS or EI spectrum against MassBank and GNPS.
- [Formula Annotation](/docs/ms1_formula): ranked molecular formulas.
- [Structure Prediction](/docs/ms2_structure): ranked structures from MS-FINDER.
- [De Novo Structure](/docs/ms2_denovo): generated structures from FRIGID.
- [EI-MS Prediction](/docs/ms2_ei): predicted 70 eV spectrum.
- [MS/MS Prediction](/docs/ms2_esi): predicted spectra at 10, 20 and 40 eV.
- [Fragment Explanation](/docs/ms2_fragments): which peaks a structure explains.
- [Alternative Explanations](/docs/exclusivity): same-formula compounds, isobars, isomers.

### Orthogonal

- [CCS](/docs/ccs): reference or predicted CCS and a match check.

### Dataset Search

- [Spectrum in Which Set](/docs/gnps_fastmasst): public datasets and library hits for a spectrum.
- [Dataset Finder](/docs/gnps_redu): public datasets by sample metadata.
- [Fragment Pattern](/docs/gnps_massql): MS2 spectra in a file that match a fragment query.

### Calculator

- [Calculator](/docs/calculator): ppm error, labelled masses, concentration, dilution, collision energy, Kovats RI.
