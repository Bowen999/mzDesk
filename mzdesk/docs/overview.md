# Overview

mzDesk helps chemists and biologists identify and check small molecules from mass-spectrometry data, without writing code.

![What mzDesk does](fig:overview-scene)
*Give a spectrum, a compound or an m/z value; mzDesk runs a workflow or a module and returns the answer, on mzdesk.com or on your own computer.*

## Workflows and modules

A **module** does one job. A **workflow** chains several modules to answer one question.

![From input to answer](fig:overview-flow)
*You give an input. A module or a workflow gives the answer.*

mzDesk has 28 modules (six groups and the Calculator) and 3 workflows. Every page says in one line what it needs and what it returns.

![28 modules and 3 workflows](fig:catalog)
*Each module group has its own colour. The workflows are solid blue.*

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
