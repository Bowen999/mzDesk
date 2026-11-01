# Identity

Who is this compound: names, database IDs, synonyms, classification and background.

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

- PubChem vendor availability is a link-out, not a certified-standard claim. Product grade and purity must be checked with each vendor.

### GHS safety

- Reporter-specific GHS classifications from PubChem. Reports can differ by jurisdiction, impurities and additives; this is not a safety data sheet.

### Run notes

- Human-associated means an HMDB metabolite record exists. This alone does not prove endogenous biosynthesis.
- Vendor purity is not standardised: verify it with each vendor.
- When no exact stereoisomer matches, UniChem IDs are matched on the first InChIKey block (connectivity only) and may belong to a stereoisomer.
