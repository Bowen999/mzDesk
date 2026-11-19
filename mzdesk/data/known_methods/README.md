# Known methods: JP18 chromatography conditions

`jp18_methods.jsonl.gz` has one JSON object per "Operating conditions" block of the Japanese Pharmacopoeia, 18th edition (English
version) that describes a liquid-chromatography (LC) or gas-chromatography (GC) method. It is the dataset behind the
"Known methods" page (compound -> published methods). Rebuild: `python -m mzdesk.methods.jp_build <folder with the PDFs>`.

* Edition: JP18 (Ministerial Notification No. 220, 2021), English version, with Supplement I (Notification No. 355, 2022) and Supplement II
  (Notification No. 238, 2024) merged into a "current" set. Errata (正誤表) published later by MHLW are NOT applied.
* Build date: 2026-10-24.  Records: 1924 (1780 from the main volume, 78 from
  Supplement I, 66 from Supplement II); LC 1796, GC 128.
  NMR (qNMR) blocks and non-chromatographic blocks (viscometer, capillary electrophoresis, thermogravimetry, polarimeter) are dropped.
* Structure mapping: 1815 records have a PubChem structure for their primary analyte (see `map_method`).

## Sources (official MHLW PDFs)

| edition | content | URL |
|---|---|---|
| JP18 | General Notices + General Tests | https://www.mhlw.go.jp/content/11120000/000945683.pdf |
| JP18 | Official Monographs A-L | https://www.mhlw.go.jp/content/11120000/000904440.pdf |
| JP18 | Official Monographs M-Z | https://www.mhlw.go.jp/content/11120000/000904449.pdf |
| JP18 | Crude Drugs and Related Drugs | https://www.mhlw.go.jp/content/11120000/000904450.pdf |
| JP18 | General Information | https://www.mhlw.go.jp/content/11120000/000904454.pdf |
| JP18 Supplement I | Supplement I | https://www.mhlw.go.jp/content/11120000/001109611.pdf |
| JP18 Supplement II | Supplement II | https://www.mhlw.go.jp/content/11120000/001354067.pdf |

Index page: https://www.mhlw.go.jp/stf/seisakunitsuite/bunya/0000066530.html

## Licence

The Japanese Pharmacopoeia is a Ministerial Notification (厚生労働省告示). Under Article 13 of the Japanese Copyright Act, laws, notifications
and the like issued by national or local authorities are not protected by copyright. This dataset therefore stores the extracted
parameters plus the text of each condition block (`condition_text`) and page references (`file`, `page`, `page_printed`); the PDFs are not redistributed.
The English version is a convenience translation: where it differs from the Japanese original, the Japanese is authoritative. Check the
printed page before relying on a value.

## Records

One object per block. Cross-referenced blocks ("Proceed as directed in the operating conditions in the Assay") are resolved: the missing fields are
copied from the referenced block (`resolved`, `resolved_from`, `inherited_fields`); the original reference stays in `xref_target` / `xref_fields`.
Supplement II overrides Supplement I overrides the main volume (matched by monograph, section and item); superseded blocks are not in this file
(`revised_in`, `replaces` mark the winning records). `null` = the block does not state the value.

Fields (in order):

- `block_id`
- `source`
- `file`
- `page`
- `page_printed`
- `monograph`
- `entry`
- `section`
- `item`
- `item_no`
- `test_label`
- `kind`
- `jp_name`
- `formula`
- `mw`
- `cas`
- `chemical_name`
- `analytes`
- `analyte_source`
- `technique`
- `block_type`
- `xref_target`
- `xref_fields`
- `xref_time_span`
- `xref_except`
- `resolved`
- `resolved_from`
- `inherited_fields`
- `resolve_note`
- `detector_type`
- `wavelength_nm`
- `column_count`
- `column_material`
- `column_id_mm`
- `column_length_mm`
- `packing_class`
- `particle_or_film_um`
- `column_temp_c`
- `injection_temp_c`
- `detector_temp_c`
- `elution_mode`
- `n_mobile_phases`
- `mobile_phase_ph`
- `organic_solvents`
- `gradient_table`
- `carrier_gas`
- `split_ratio`
- `flow_rate_value`
- `flow_rate_unit`
- `flow_target_compound`
- `flow_target_rt_min`
- `run_time_min`
- `run_time_rt_multiple`
- `cid`
- `smiles`
- `inchikey`
- `parent_inchikey`
- `map_method`
- `map_reason`
- `analyte_structures`
- `substances`
- `revised_in`
- `replaces`
- `cond_label`
- `notes`
- `qa`
- `condition_text`

Notes on selected fields

* `file`, `page`, `page_printed`: PDF file name, page index in that PDF (1-based) and the page number printed in the running header (differs by a constant per file).
* `kind`: monograph / crude drug / extract / general test / reagent / general information.
* `item`, `item_no`, `test_label`: e.g. `(3) Related substances`, 3, `Purity (3) Related substances`. Sub-procedure headings are kept when they name the block:
  `(3) Related substances (iii)`, `(ii) Determination of molecular mass`; long titles that wrap over two printed lines are joined.
* `cas`: CAS RN(s) in brackets in the monograph header; preparations (tablets ...) inherit the CAS of the bulk-drug monograph.
* `analytes`: lower-case names of what the test MEASURES. Assay / dissolution / content uniformity: the monograph substance (bulk monograph = its title, a
  preparation = the bulk monograph it is made from, a combination preparation = the active(s) the test text talks about). A Purity item whose title names
  compounds ("Methanol", "Ethylene glycol, diethylene glycol and related substances", "1-Vinyl-2-pyrrolidone") measures those compounds and not the substance
  (the substance is kept only when the test compares the impurity peak with the substance's own peak). "Residual solvents" take the solvents named in the
  quantity formulae. `related substances` is added for related-substance tests (monographs and reagent entries). Crude drugs / extracts: marker compounds from
  the item title, else the quantity formula or the peak statement. Reagent entries: the reagent compound. Unnamed impurity classes ("Diastereomer",
  "Multimers") keep the substance.
* `detector_type`: UV, PDA, visible, fluorescence, RI, conductivity, MS, FID, TCD, ECD, other. `packing_class`: C18, C8, C4, phenyl, CN, NH2, silica, diol, PFP, anion-exchange, cation-exchange, SEC, chiral, polymer-RP, GC-PEG, GC-6%cyanopropylphenyl, GC-5%phenyl, GC-dimethylpolysiloxane, GC-other, other. `silica` is only the plain silica gel for liquid chromatography; any bonded or modified gel (carbamoyl, trimethylsilanized,
  sulfonamide / hexadecylsilanized, ...) is `other`, albumin / glycoprotein / cellulose / cyclodextrin phases are `chiral`, "sulfo" alone is not a cation exchanger.
* `column_temp_c`: number, `"program"` (temperature programme) or `"room"`. `elution_mode`: isocratic / gradient (a gradient table, or several mobile phases switched).
* `gradient_table`: rows `[t_start, t_end, A_start, A_end, B_start, B_end]` in min / vol%; further mobile phases append `[C_start, C_end]` ...
* `flow_rate_value` + `flow_rate_unit` (mL/min or cm/s) when a flow is stated; `flow_target_compound` + `flow_target_rt_min` when the text says "adjust so that the retention time of X is about N minutes" (midpoint of a range).
* `xref_time_span`: the cross-reference also covers the time span of measurement. `xref_except`: text after "except".
* `cid`, `smiles` (isomeric), `inchikey`, `parent_inchikey` (largest fragment, neutralised; use it for matching), `map_method` (`cas`, `name`, `name-stripped`, `none`) describe the primary analyte; `analyte_structures` has one entry per analyte. `substances` (drug monographs only) lists the structure of the monograph substance the
  sample consists of, which differs from `analytes` for impurity tests (the Captopril Purity (4) test measures captopril disulfide; `substances` still names captopril). Biologics, polymers and mixtures have no structure (`map_reason`).
* `qa`: automatic QA flags (missing fields, out-of-range values, unresolved references ...); see QA_REPORT.md in the build work folder.

## Known limitations

* Only blocks introduced by an "Operating conditions" heading are covered. A few monographs print the conditions under another heading (the Assay of
  Propiverine Hydrochloride has "System suitability" there), so references to them stay unresolved (`resolved: false`, listed in the QA report).
* The mobile-phase composition, buffer concentrations, injection volume and system-suitability limits are not parsed; read `condition_text`.
* Alternatives take the first value ("0.32 mm (or 0.53 mm)"); ranges of retention time or flow take the midpoint; for a guard / pre-column plus a
  separation column the analytical column is reported (`notes` says so). Stepwise switching of several mobile phases counts as `gradient`.
* Cross-references that name a sub-procedure only by its position (e.g. "prescription 1)", "Procedure (ii)") are resolved when the position matches a
  block of the same test, otherwise left unresolved rather than guessed.
* Text comes from the PDF text layer; glyph fixes (percent, degree, prime, Greek letters, fi / fl ligatures) are in `extract.py`. Printed page numbers are
  those of the running header.
