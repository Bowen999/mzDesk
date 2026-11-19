# Identity

Give a compound. Get its names, database IDs, synonyms, classification and background.

## Results

### Classification

- ClassyFire and LOTUS records are source annotations; NPClassifier is a natural-product-oriented structure classifier.
- LIPID MAPS rows appear for LMSD records. mzDesk rows come from the local lipid / metabolite rules. Different systems need not agree.

### Origin

- Evidence of association or occurrence, not proof of exclusive origin. Human-associated records do not by themselves prove endogenous biosynthesis.

### Pathways

- KEGG, Reactome and HMDB / SMPDB annotations. Reactome uses the numeric ChEBI identifier and human pathway projection.

### Biofluid concentrations

- Read from a local HMDB XML only; values are never inferred. Set MZDESK_HMDB_XML to the unzipped HMDB metabolites XML path.
- Normal ranges depend on cohort, matrix, method and units. Check the linked HMDB record before reuse.

### Species

- Reported species-level LOTUS occurrences, deduplicated across taxonomic providers. Absence of a row is not evidence of absence.

### Standards & availability

Is there a commercial standard for this compound? The grade comes from the vendors PubChem lists for the structure:

| Grade | Meaning |
|---|---|
| A | A catalog vendor sells this exact structure. Best chance of a real standard. |
| B | Only screening-library or other vendors list it. Often make-on-demand, not in stock. |
| B-iso | No vendor for this structure, but a catalog vendor sells an isomer (stereoisomer or isotopologue; up to 5 checked). Check it is the one you need. |
| C | No vendor listing. |

- The counts show catalog, screening and other vendors. **PubChem vendor listings** opens the full list.
- A listing is a link-out, not a certified-standard claim. Check grade, purity and stock with the vendor.

### GHS safety

- Reporter-specific GHS classifications from PubChem. Reports can differ by jurisdiction, impurities and additives.

### Run notes

- Human-associated means an HMDB metabolite record exists.
- When no exact stereoisomer matches, UniChem IDs are matched on the first InChIKey block (connectivity only) and may belong to a stereoisomer.
