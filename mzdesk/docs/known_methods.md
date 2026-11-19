# Known Methods

Give a compound. Get its published methods (Japanese Pharmacopoeia, RepoRT), or, if there are none, the methods of the most similar compounds.

## Results

### Summary

- **Match**: Exact: the compound itself has a published method. Similar: only similar compounds do. None: neither.
- **Methods found**: Exact matches plus the methods listed under the similar compounds.
- **Sources**: Which sources were searched; licence and citation.

### Exact matches

- The compound itself: same connectivity key (first 14 characters of the InChIKey after largest fragment and neutralisation); stereo is noted in the row's ?.
- The ? on a row gives match details, the reference and the printed JP conditions or the RepoRT gradient. RepoRT is LC-MS: RT is the compound's retention time in that dataset (min). JP18 has no RT.

### Similar compounds

- Click a compound to see its methods. RT is that compound's retention time.

### Run notes

- JP18 source: Japanese Pharmacopoeia, 18th edition (English version) with Supplements I and II, Ministry of Health, Labour and Welfare (Ministerial Notification; not protected by copyright under Article 13 of the Japanese Copyright Act).
- JP18 parameters were extracted from the PDF text layer of the 'Operating conditions' blocks, with the printed page. The English version is a convenience translation: check the printed page before relying on a value. Errata published after the supplements are not applied. Page links open the MHLW PDFs.
- RepoRT citation: F. Kretschmer, E.-M. Harrieder, M. A. Hoffmann, S. Böcker, M. Witting. RepoRT: a comprehensive repository for small molecule retention times. Nat. Methods 21, 153-155 (2024). https://doi.org/10.1038/s41592-023-02143-z
- RepoRT data are CC BY-SA 4.0 (https://github.com/michaelwitting/RepoRT). Each dataset keeps its own origin (the link of its row). Column, eluent and gradient texts are summaries of the dataset metadata. RepoRT is LC-MS; retention times are in minutes. Reuse of RepoRT-derived values requires attribution and the same licence.
- Similarity score = 0.6 × Tanimoto (Morgan radius 2, 2048 bits) + 0.4 × the mean of five property terms: exp(−|ΔlogP| / 1), exp(−|ΔTPSA| / 30), exp(−|ΔMW| / 100), and 1 / 0 for the same / a different charge class (sign of the dominant net charge) at pH 2.7 and at pH 7.4. Δ = library compound − query.
- A similar compound is listed when its score is ≥ 0.5 and its Tanimoto is ≥ 0.25, whatever the property score. At most 10 compounds, 3 methods per compound and 2 per RepoRT dataset are shown.
- The query's logP is the Crippen value. Molecules are reduced to the largest fragment and neutralised first; charge classes come from MolGpKa.
- Without MolGpKa weights, the charge class comes from functional groups: carboxylic / sulfonic / phosphoric acid is negative; aliphatic amine / amidine / guanidine / quaternary N is positive.
