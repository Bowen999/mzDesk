"""ESI MS/MS Annotation workflow: the steps, the data contract they share and the orchestrator that runs them.

schema         Candidate and LayerResult (what every layer returns and what the result rows look like)
levels         Schymanski confidence levels and their assignment
adduct_step    step 1: adduct hypotheses (prior x MS2 signals x BUDDY formula probability x isotope fit)
formula_step   steps 1+2: runs BUDDY per retained adduct, feeds step 1 and gives the shared formula list
library_step   step 3: library search and the Level-2 gate
msfinder_step  step 4: MS-FINDER structures for the shared formulas (Level 3) and the "good enough" test
frigid_step    step 5: FRIGID structure hypotheses for the top formula (only when steps 3-4 gave nothing good enough)
merge          step 6a: dedupe the layers' candidates by InChIKey first block
rerank         step 6b: the layers' base score decides the order; detectability factors (omics, ionization, adduct) only break near-ties
answer         step 8: the answer block of the result, and the isomer annotation of the rows (step 7)
pipeline       the orchestrator: parses the form, runs steps 1-8 in order (step 7 is mzdesk.exclusivity.isomers) and builds the Result JSON
"""
