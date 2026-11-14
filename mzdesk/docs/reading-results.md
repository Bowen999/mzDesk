# Reading Results

This page explains the answer card, the confidence levels, the scores and the tables you see after a run.

## The answer card

MS/MS Annotation starts with one line, for example **Caffeine · Level 2a · 1.00**.

- **Name**: the top-ranked structure. Without a structure, the formula or the exact m/z.
- **Badge**: the confidence level.
- **Number**: the final score, from 0 to 1. Higher is better.
- **Line below**: where the answer comes from, and what the spectrum cannot separate.
- **Details**: the full statement, and the layer, level, score, isomer check, library gate and CCS check.

EI MS Annotation, Targeted Method and Lipid Fragment Explanation start with a statement, a badge when there is one, and a short list of facts.

## Confidence levels

![The confidence levels](fig:evidence-levels)

The levels follow Schymanski et al., Environ. Sci. Technol. 2014. Level 1 needs a reference standard, so mzDesk never reports it.

| Level | Meaning | Comes from |
|---|---|---|
| Level 2a | Probable structure | A library spectrum passed the Level-2 gate, and the entry has a structure. |
| Level 3 | Tentative candidate | An MS-FINDER candidate, or a library hit that failed the gate or has no structure. |
| Hypothesis | Formula-level evidence only | A FRIGID structure. |
| Level 4 | Molecular formula only | No structure, only formulas. |
| Level 5 | Exact mass only | Only the m/z. |

- A structure found by several tools gets the best level among them.
- The isomer check and the re-ranking never change a level.
- EI MS Annotation uses Level 2a, Level 3 and **none** (no identification).

### The Level-2 gate

A library hit reaches Level 2a only if it passes these checks. Open **Steps**, then **Library search**, to see each result.

1. **Similarity** reaches the threshold (default 0.80) on at least one ticked algorithm.
2. **Matched peaks**: at least 4, or all peaks if the spectrum has fewer.
3. **Precursor** is within the ppm tolerance of the hit's formula and adduct.
4. **Adduct** is one of the adducts found in step 1.
5. **Formula** is in the formula list of step 2.
6. **Isotope pattern** matches at 0.80 or more.

A check that cannot be evaluated does not block the hit. For example, the isotope check needs an MS1 pattern. Change the threshold and the minimum peaks under **Advanced**.

## Hypothesis

A FRIGID structure is a hypothesis. FRIGID writes structures for a formula and a fingerprint predicted from the spectrum. It does not search a database, so a candidate can be a molecule that exists nowhere. Check it against a reference standard, a library spectrum (Library Search) or a predicted spectrum (MS/MS Prediction).

## Ranked candidates

| Column | Meaning |
|---|---|
| Rank | The order in the list. |
| Name | Hover to see the structure. |
| Level | The confidence level. Hover for the reason. |
| Final score | 0 to 1. |
| Fragments | **View** opens the [fragment view](#fragment-view). |

Click a row for its SMILES, formula, adduct, layer (library, MS-FINDER, FRIGID or formula), detectability and isomer role.

The order is not by final score alone. Candidates within 10 % of the best base score are ordered by base score × detectability. The others follow by base score. So the final score does not always fall down the list.

## Isomer check

The check compares the top structure with its isomers using the fragments of your spectrum.

| Word | Meaning |
|---|---|
| single | The fragments tell the top structure apart from every isomer compared. |
| indistinguishable | The spectrum cannot separate it from the listed isomers. |
| insufficient | The check could not be done. |
| not run | The check was skipped. |
| challenger | An isomer that explains peaks the top structure cannot. The card calls this conflicting evidence. |

The result holds only for the isomers that were compared.

## Fragment view

Click **View** in the Fragments column of the ranked candidates, or of the Isomers table in [Alternative Explanations](/docs/exclusivity). A window opens, computed when you click:

- the structure, with each explained fragment in its own colour;
- the spectrum in the same colours, with m/z and mass error on each peak;
- one card per fragment.

The fragment view uses MetFrag. The isomer comparison uses MS-FINDER, so the explained peaks can differ.

## Module results

- **Cards**: single values. A `?` next to a label or title gives a short explanation.
- **Tables**: long tables fold or page. Hover over a name or SMILES to see the structure.
- **Chips and tiles**: good, possible or poor (yes, partial or no). Hover for the reason.
- **Figures**: hover a peak for its m/z and, where known, its fragment.
- **Warnings**: a box above the results.
