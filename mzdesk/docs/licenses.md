# Licenses and Based On

This page gives the license of mzDesk and of the tools, models and data that its pages rely on.

## mzDesk

mzDesk is free software under the MIT License. Copyright (c) 2026 mzDesk contributors.

- You may use, copy, change and share it, also commercially.
- Keep the copyright line and the license text in every copy.
- The full text is the `LICENSE` file of the [source code](https://github.com/Bowen999/mzDesk).

## Based on

The foot of every page reads **Based on**, then the names of what the page relies on, most important first, at most four. **mzDesk** in that list means the page also uses mzDesk's own rules.

- The [Session report](/docs/session-report) records the version of each name.
- The tables below give the license of each name.
- When you reuse data from a source, for example a library spectrum, the license of that source applies.

## Tools and models

| Name | Used by | License | Link |
|---|---|---|---|
| BUDDY (msbuddy) | Formula Annotation, De Novo Structure, MS/MS Annotation | Apache-2.0 | [GitHub](https://github.com/Philipbear/msbuddy) |
| MS-FINDER | Structure Prediction, MS/MS Annotation | Program CC BY 4.0, source code LGPL-3.0 | [GitHub](https://github.com/systemsomicslab/MsdialWorkbench) |
| FRIGID | De Novo Structure, MS/MS Annotation | Code and DLM checkpoint Apache-2.0; MIST checkpoint MIT | [GitHub](https://github.com/coleygroup/FRIGID) |
| MetFrag | Fragment Explanation, Alternative Explanations, Lipid Fragment Explanation, the fragment view | | |
| CFM-ID and MSRB | MS/MS Prediction | | |
| NEIMS | EI-MS Prediction, Ionization, EI MS Annotation | Original code Apache-2.0 | [GitHub](https://github.com/brain-research/deep-molecular-massspec) |
| MolGpKa | Separation, Ionization, Adducts, Targeted Method | MIT | [GitHub](https://github.com/Xundrug/MolGpKa) |
| matchms | Similarity, Library Search | Apache-2.0 | [GitHub](https://github.com/matchms/matchms) |
| ms_entropy | Similarity, Library Search | Apache-2.0 | [GitHub](https://github.com/YuanyueLi/MSEntropy) |
| MS2DeepScore | Similarity | Apache-2.0 | [GitHub](https://github.com/iomega/ms2deepscore) |
| Spec2Vec | Similarity | | [GitHub](https://github.com/iomega/spec2vec) |
| RDKit | Structure & Formats, Properties & Groups, Lipid Short Name, Platform, Separation, Non-MS Detection, Derivatization, Known Methods, EI MS Annotation | BSD-3-Clause | [GitHub](https://github.com/rdkit/rdkit) |
| IsoSpecPy | Mass Search, Theoretical, Calculator | BSD-2-Clause | [GitHub](https://github.com/MatteoLacki/IsoSpec) |
| pyOpenMS | Fragment Pattern (reads mzML and mzXML) | BSD-3-Clause | [Project page](https://pyopenms.readthedocs.io) |
| MassQL | Fragment Pattern | mzDesk runs its own code for the MS2 part of the query language. | |
| LIPID+ | Lipid Short Name, Lipid Fragment Explanation | Ships inside mzDesk | |

- MetFrag ships inside mzDesk.
- In the local version, MS-FINDER, FRIGID, NEIMS and MolGpKa are downloaded from their official sources the first time a page needs them. CFM-ID and MSRB come as a Docker image from their authors. See [First-use downloads](/docs/local-version#first-use-downloads).

## Data and online services

| Name | Used by | License or terms | Link |
|---|---|---|---|
| MassBank | Library Search, MS/MS Annotation, EI MS Annotation | Set per record. The records used in the examples are CC BY 4.0 or CC BY-SA. | [MassBank](https://massbank.eu) |
| GNPS spectral libraries | Library Search, MS/MS Annotation | | [GNPS](https://gnps.ucsd.edu) |
| fastMASST, Pan-ReDU (GNPS2) | Spectrum in Which Set, Dataset Finder | | [Pan-ReDU](https://redu.gnps2.org) |
| PubChem | Identity, Mass Search, Alternative Explanations, Targeted Method, name and ID lookup | | [PubChem](https://pubchem.ncbi.nlm.nih.gov) |
| UniChem | Identity, name and ID lookup | | [UniChem](https://www.ebi.ac.uk/unichem/) |
| ClassyFire | Identity | | |
| LOTUS | Identity | | [LOTUS](https://lotus.naturalproducts.net) |
| CCSbase source data (c3sdb) | CCS; MS/MS Annotation with an observed CCS | MIT | [GitHub](https://github.com/dylanhross/c3sdb) |
| RepoRT | Known Methods, Targeted Method | CC BY-SA 4.0 | [GitHub](https://github.com/michaelwitting/RepoRT) |
| Japanese Pharmacopoeia, 18th edition | Known Methods | Not protected by copyright (Japanese Copyright Act, Article 13). | [MHLW](https://www.mhlw.go.jp) |

## In the interface

| Name | For | License |
|---|---|---|
| Ketcher | Drawing structures | Apache-2.0 |
| Plotly.js | Spectra and figures | MIT |
| 3Dmol.js | 3D structures | |
| React | The interface | MIT |

## Other packages

mzDesk also uses FastAPI (MIT), Uvicorn (BSD-3-Clause), NumPy and SciPy (BSD), scikit-learn (BSD-3-Clause, the CCS model), gensim (LGPL-2.1, for Spec2Vec), Plotly for Python (MIT) and Matplotlib (Matplotlib license). The local version gets a Java runtime from the jdk4py package, used by MetFrag.
