# Changelog

This page lists what changed in each version of mzDesk, newest first.

## v0.5 (2026-11-01)

- v0.5.8: The Platform tile **CIL LC-MS** is now called **CIL Metabolomics**. Every derivatization reaction now shows the complete product structure next to its functional-group scheme: the GC reagents (TMS, MeOx + TMS, TBDMS, TFAA) and the LC-UV / FLD reagents (OPA, FMOC-Cl, DNPH), as the four CIL channels already did.
- v0.5.7: New workflow **Lipid Fragment Explanation**. Give a lipid short name (or a SMILES or a drawing) and its MS2 spectrum; get the LIPID+ and LIPID MAPS names, an illustrative structure with its assumed parts shaded, and the peaks that structure explains (MetFrag, tree depth 3 by default). **Auto** picks the adduct that fits the precursor m/z.
- v0.5.6: MS-FINDER, FRIGID, CFM-ID, NEIMS and RepoRT now install and run on Windows (MS-FINDER uses the .NET Framework that comes with Windows; CFM-ID needs Docker Desktop). First-use downloads show a real progress bar with the file and step, in the app and in the terminal. mzDesk now opens at `http://localhost:8731` (it was 8000, which other programs often use).
- v0.5.5: MS-FINDER gives the same structure scores every time it runs on the same spectrum. Before, MS-FINDER's parallel scoring mixed values between candidates, so scores changed slightly from run to run. Run `mzdesk setup` once if Structure Prediction asks for it; nothing is downloaded.
- mzdesk.com is live, with every engine installed: spectral libraries, CCS, MS-FINDER, FRIGID, CFM-ID and the deep-learning scores.
- The bug button sends feedback to the maintainer. A page that fails to draw, or a run that fails, shows a reference code such as `E-1A2B3C`.
- MS/MS Annotation: **Run FRIGID** under **Advanced**. It is on in the local version and off on mzdesk.com.
- The interface comes in six languages: English, 中文, 日本語, Deutsch, Français and Español. The Docs stay in English.
- The Session report lists the versions of the tools, models and data behind each result.
- mzdesk.com has no Internet search switch: its online lookups always run.

## v0.4 (2026-11-01)

- New look: one blue theme, one warning colour, no subtitles.
- Every page opens with one line, what you give and what you get, and ends with **Based on** and a link to its Docs.
- Workflows show a **Steps** timeline with the status, note and time of each step.
- Smoother progress bars, and an opening splash that plays once per browser session.
- Tidier results: one-column tables become chips and long tables fold. The layout works on phones and tablets.
- MS/MS Annotation: **Minimum matched peaks** under **Advanced**.

## v0.3 (2026-10-30)

- New workflows: **EI MS Annotation** and **Targeted Method**.
- MS/MS Annotation shows a compact answer, shows a structure when you hover a name, and lets you tick several of ten library similarity algorithms.
- MS/MS Annotation can check each candidate against an observed CCS.
- New module: **Lipid Short Name**.
- **Session report**: one download with the results of every page on screen.
- DMPK, ICP-MS and MRM Transitions are no longer offered. Alternative Explanations moved into Spectrum, and Calculator stands alone.

## v0.2 (2026-10-24)

- New modules: Identity, Structure & Formats, Properties & Groups, Spectrum Cleanup, Mass Search, Theoretical, Library Search (MassBank and GNPS), Structure Prediction (MS-FINDER), De Novo Structure (FRIGID), EI-MS Prediction, MS/MS Prediction, Fragment Explanation and Alternative Explanations.
- Detectability pages: Platform, Separation, Ionization, Adducts, Non-MS Detection, Derivatization and Known Methods.
- New pages: CCS, Calculator and the GNPS2 Dataset Search pages.
- One workflow, MS/MS Annotation, replaced the four earlier workflows.
- Any page can take a structure as a name or an ID. The Internet search switch asks before an online run.
- A **View** button shows the fragments that explain each peak. Properties & Groups draws 176 functional groups on the structure.

## v0.1 (2026-10-24)

- First release: 25 module pages, 6 of them computing.
- Computing: theoretical m/z, BUDDY formulas, spectrum similarity with a mirror plot, MetFrag fragment explanation, 2D and 3D structures with descriptors, and lipidomics or metabolomics platform fit.
- Reads MGF and MSP files. **Load example** fills a page.
