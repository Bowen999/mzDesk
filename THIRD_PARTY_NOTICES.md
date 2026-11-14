# Third-party notices

mzDesk is MIT-licensed (see `LICENSE`). This file lists the third-party software, code and data that are **distributed with mzDesk**:
in this repository, in the Python package (wheel) and in the built web interface. The full licence texts are in `third_party_licenses/`.

Software and data that mzDesk downloads on first use (MS-FINDER, FRIGID and its checkpoints, NEIMS, MolGpKa, the MassBank and GNPS
spectral libraries, CCSbase, RepoRT, Pan-ReDU, the CFM-ID Docker image) are not distributed with mzDesk; they come from their authors'
official sources under their own licences. Python packages are installed by pip under their own licences. A full inventory with
sources is kept in `paper/supplementary/licenses/` (checked 2026-11-12).

## Java library: MetFrag

`mzdesk/lib/MetFragCommandLine-2.5.2.jar`, unmodified.

- **MetFrag CommandLine 2.5.2**: MetFrag authors, Leibniz Institute of Plant Biochemistry (IPB Halle).
  LGPL-2.1-or-later (`third_party_licenses/LGPL-2.1-or-later.txt`). Source: https://github.com/ipb-halle/MetFragRelaunched;
  licence statement: https://ipb-halle.github.io/MetFrag/. You may replace this JAR with your own build of MetFrag.
- Libraries contained in the JAR (as published by the MetFrag authors):

| Library | Version | Licence | Licence text |
|---|---|---|---|
| Chemistry Development Kit (CDK) | 2.9 | LGPL-2.1 | `LGPL-2.1-or-later.txt` |
| JNA-InChI | 1.2 | LGPL-2.1 | `LGPL-2.1-or-later.txt` |
| vecmath (javax.vecmath) | 1.5.2 | GPL-2.0 with Classpath Exception | `GPL-2.0-only.txt`, `Classpath-exception-2.0.txt` |
| javax.annotation-api | 3 | CDDL-1.1 or GPL-2.0 with Classpath Exception | `CDDL-1.1.txt` |
| jakarta.json (Eclipse JSON-P) | 2.0.1 | EPL-2.0 | `EPL-2.0.txt` |
| FreeHEP vector graphics (CERN, SLAC) | 2.4 | Apache-2.0 or LGPL-2.0 | `Apache-2.0.txt` |
| Beam (EMBL-EBI) | 1.3.5 | BSD-style | `BSD-3-Clause.txt` |
| jsonld-java | 0.13.4 | BSD-3-Clause | `BSD-3-Clause.txt` |
| protobuf-java (Google Inc.) | 3.21.9 | BSD-3-Clause | `BSD-3-Clause.txt` |
| SLF4J; signatures | 1.7.36; 1.1 | MIT | `MIT.txt` |
| Apache Jena, Commons, HttpComponents, Log4j; Guava; Gson; Jackson; Joda-Time; json-simple; titanium-json-ld; cmlxom; euclid; hierarchical-clustering | various | Apache-2.0 | `Apache-2.0.txt` |

## JavaScript library: 3Dmol.js

`mzdesk/lib/3Dmol-min.js`, the unmodified `build/3Dmol-min.js` of **3Dmol.js 2.4.2** (David Koes and contributors; Copyright (c) 2014,
University of Pittsburgh and contributors). BSD-3-Clause; it includes code from GLmol (MIT or LGPL-3.0), Three.js (MIT) and jQuery (MIT).
The full notice is `mzdesk/lib/3Dmol-min.js.LICENSE.txt`. Source: https://github.com/3dmol/3Dmol.js.

## MS-FINDER files for Windows

`mzdesk/lib/msfinder/win/`. MS-FINDER v3.73, RIKEN (Tsugawa lab), https://github.com/systemsomicslab/MsdialWorkbench; details and
SHA-256 in `mzdesk/lib/msfinder/win/NOTICE.md`.

- `MsfinderConsoleApp.exe`: compiled by mzDesk from the 9 unmodified source files of `src/MSFINDER/MsfinderConsoleAppCore` at tag
  `MSFINDER-v3.73`. LGPL-3.0 (`LGPL-3.0-only.txt`, which supplements `GPL-3.0-only.txt`). The source is at the tag above.
- `StructureFinder.fixed.dll`: a **modified** copy of the `StructureFinder.dll` of the official MS-FINDER v3.73 release. The only change
  marks the field `FragmentPeakMatcher.maxBondDissociationEnergy` as `[ThreadStatic]` (`mzdesk/lib/msfinder/ThreadStaticPatch.cs`).
  The MS-FINDER assemblies are CC BY 4.0 (`CC-BY-4.0.txt`).
- Citation: Tsugawa H. et al. Anal. Chem. 2016, 88, 7946-7958.

## FRIGID patch

`mzdesk/lib/frigid/mps_patch.diff` changes `src/dlm/sampler.py` of FRIGID (Coley group, https://github.com/coleygroup/FRIGID), which is
Apache-2.0 (`Apache-2.0.txt`). The patch says what it changes; mzDesk applies it to the FRIGID code that the user downloads.
`mzdesk/lib/frigid/best_ngboost_MSG.joblib` was trained by the mzDesk authors on MassSpecGym labels (MIT) and is part of mzDesk (MIT).

## LIPID MAPS

Some tables of `mzdesk/lipid_name_map/data/` (`alias/lmsd.csv`, `lm_heads.csv`) contain information derived from the LIPID MAPS®
Structure Database (LMSD), https://www.lipidmaps.org, licensed under CC BY 4.0 (`CC-BY-4.0.txt`).

## Japanese Pharmacopoeia

`mzdesk/data/known_methods/` holds chromatography conditions from the Japanese Pharmacopoeia, 18th edition (Ministry of Health, Labour and
Welfare, Japan), which is not protected by copyright (Japanese Copyright Act, Article 13). Sources: `mzdesk/data/known_methods/README.md`.

## Web interface

The built interface (`mzdesk/static/`) contains, among others:

| Component | Version | Copyright | Licence |
|---|---|---|---|
| React, React DOM | 18.3.1 | Facebook, Inc. and its affiliates | MIT |
| React Router | 6.30.6 | React Training LLC; Remix Software Inc. | MIT |
| Plotly.js | 3.3.1 | Plotly Technologies Inc. | MIT |
| Ketcher | 3.18.0 | EPAM Systems, Inc. | Apache-2.0 |
| Indigo (indigo-ketcher, WebAssembly) | 1.46.0 | EPAM Systems, Inc. | Apache-2.0 |
| IBM Plex Sans, IBM Plex Mono (fonts) | 5.3.0 | IBM Corp. | OFL-1.1 (`OFL-1.1.txt`) |

All 330 packages of the web interface (`web/package-lock.json`, non-development) are under MIT, ISC, Apache-2.0, BSD-2-Clause,
BSD-3-Clause, BlueOak-1.0.0, OFL-1.1, 0BSD or Unlicense; the list is `paper/supplementary/licenses/dependencies_npm.csv`.
