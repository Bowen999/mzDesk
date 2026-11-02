# Overview

mzDesk helps chemists and biologists identify and check small molecules from mass-spectrometry data, without writing code.

![What mzDesk does](fig:overview-scene)
*Give a spectrum, a compound or an m/z value; mzDesk runs a workflow or a module and returns the answer, on mzdesk.com or on your own computer.*

## What you can do

- Annotate an MS/MS or EI spectrum: adduct, formula, structure candidates, confidence level.
- Search a spectrum against MassBank and GNPS, compare two spectra, clean a spectrum and convert its format.
- Look up a compound: names, IDs, properties, functional groups.
- Check whether a compound can be measured, and how: platform, column, pH, polarity, ions.
- Predict EI and MS/MS spectra, and explain a spectrum with a structure.

## Workflows and modules

A **module** does one job. A **workflow** chains several modules to answer one question.

![From input to answer](fig:overview-flow)
*You give an input. A module or a workflow gives the answer.*

mzDesk has 3 workflows and 28 modules: six groups and the Calculator. Every page says in one line what it needs and what it returns.

![The modules by group](fig:module-map)

| Group | For |
|---|---|
| Compound | names, IDs, structures, properties, lipid names |
| Detectability | platform, separation, ionization, adducts, other detectors, derivatization, published methods |
| MS1 | match a precursor m/z to compounds; theoretical m/z and isotopes |
| Spectrum | clean, compare, search and annotate spectra; predict spectra |
| Orthogonal | CCS check |
| Dataset Search | find public datasets and spectra (GNPS2) |
| Calculator | ppm error, dilution, collision energy, Kovats RI |

Details: [Workflows and Modules](/docs/workflows-and-modules).

## Two ways to use it

- **mzdesk.com**: open it in a browser. Nothing to install.
- **Local version**: install it once and run it on your own computer. See [Local Version](/docs/local-version).

## Start here

1. [Quick Start](/docs/quick-start): your first result, and what is on the screen.
2. [Reading Results](/docs/reading-results): what the answer means.
