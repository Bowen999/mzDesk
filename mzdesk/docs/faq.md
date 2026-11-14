# FAQ

This page answers the questions people ask most, with a link to the page that has the details.

## General

### Do I need an account? Does it cost anything?

No account and no login. mzDesk is free, and its code is open source under the MIT license. See [Licenses and Based On](/docs/licenses).

### Which languages does mzDesk speak?

The interface comes in English, 中文, 日本語, Deutsch, Français and Español. Change it with the language button in the sidebar footer. The Docs are in English only.

## Inputs

### Can I upload an mzML or other raw data file?

Not for annotation. Export the spectrum or a peak list from your instrument or processing software, then paste it or load an MGF / MSP file. The one exception is [Fragment Pattern](/docs/gnps_massql), which takes one `.mzml`, `.mzxml` or `.mgf` file of at most 20 MB.

### Which spectrum formats can I paste?

Peak lists, MGF, MSP, MassBank records, SIRIUS `.ms`, JSON, CSV and more. mzDesk detects the format by itself. See [Inputs and Formats](/docs/inputs#spectra).

### Can I run many spectra at once?

No. Each run takes one spectrum or one compound. A file with several spectra uses only the first, and the page warns you.

### Does it work in negative ion mode?

Yes, apart from FRIGID: De Novo Structure and the FRIGID step of MS/MS Annotation take positive ions only. In Library Search, the MS2DeepScore and Spec2Vec scores are N/A for a negative-mode spectrum.

## Results

### Why is the answer never Level 1?

Level 1 needs a reference standard, so mzDesk never reports it. The best it reports is Level 2a. See [confidence levels](/docs/reading-results#confidence-levels).

### What is the difference between Level 2a and Level 3?

Level 2a means a library spectrum matched and passed every check of [the Level-2 gate](/docs/reading-results#the-level-2-gate). Level 3 is a tentative candidate: an MS-FINDER structure, or a library hit that failed the gate or has no structure.

### Why were the MS-FINDER and FRIGID steps skipped?

MS-FINDER runs only when no library hit passes the Level-2 gate. FRIGID runs only when MS-FINDER is not good enough, **Run FRIGID** is ticked and the ion is positive. The note of each step in **Steps** says why it was skipped.

### Why is Run FRIGID off on mzdesk.com?

A FRIGID run on the server takes about two minutes and 3 GB of its 4 GB of memory. You can tick it under **Advanced**. In the local version it is on once FRIGID is installed. See [Local Version](/docs/local-version#first-use-downloads).

### Why are MassBank and GNPS the only libraries? Where is NIST?

The spectral libraries are MassBank and GNPS. mzDesk includes no NIST or other commercial library.

### My result differs from another program, or a page shows an error. What should I do?

Click the bug button in the sidebar footer and describe what you saw. A screenshot helps. See [Feedback](/docs/feedback). The **Based on** line at the foot of each page names the tools behind the result.

### How do I save my results?

Click **Session report** in the sidebar. It downloads a `.zip` with `report.html` to read and `report.json` with the inputs, raw results and tool versions. See [Session Report](/docs/session-report).

## Privacy and the local version

### Does mzdesk.com keep my data?

Your inputs are sent to the server, computed in memory and not saved. Only compound lookups are cached, for 30 days, then deleted. No IP address is written to disk. If your data must stay on your computer, use the local version. See [Privacy and Internet Search](/docs/privacy).

### Does the local version work without the internet?

Mostly. SMILES, InChI, formulas and spectra are read on your computer. Names, IDs and USIs are looked up online, and engines download the first time a page needs them. With the Internet search switch off, mzDesk asks before a run that needs the internet.

### Which systems does the local version run on?

macOS and Linux on x86_64 or arm64, and 64-bit Windows. Every component has a Windows version, but the engines (MS-FINDER, FRIGID, NEIMS, CFM-ID, RepoRT) are new there and less tested than on macOS and Linux. See [Local Version](/docs/local-version#install).
