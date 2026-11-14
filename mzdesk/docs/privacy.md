# Privacy and Internet Search

This page lists what leaves your computer, where it goes and what mzdesk.com keeps.

## mzdesk.com

- Your inputs go to the mzDesk server over HTTPS.
- The server computes in memory. A workflow result stays in memory for up to one hour, so the page can collect it. Temporary files are deleted when the run ends.
- The [online lookups](#what-goes-online) run from the server. Your browser talks only to mzdesk.com.
- There is no account, no login and no cookie.

The server keeps:

| What | Content |
|---|---|
| Usage counts | Time, page, module, status, run time, browser family, the referring site and a visitor ID. No IP address and no inputs. |
| Lookup cache | The answers of compound lookups, filed under the request, which holds the name, ID, InChIKey or SMILES looked up. An answer is reused for 30 days, then deleted. |
| Crash reports | The error and where it happened. Input fields are listed by name and size, not by value. |
| Feedback | What you write and the screenshots you attach. See [Feedback](/docs/feedback). |

No IP address is written to disk: the web server keeps no access log. The visitor ID is a hash of a secret key, the day, the IP address and the browser. A visitor counts once a day and cannot be followed from one day to the next.

## Local version

- mzDesk runs on your computer, and your inputs stay there.
- Three things leave it: the [online lookups](#what-goes-online), the [first-use downloads](#first-use-downloads) and the [version check](#version-check).
- Usage counts and crash reports are saved in `~/.mzdesk` and never sent.
- The bug button opens the feedback form of mzdesk.com. What you write there goes to the maintainer. See [Feedback](/docs/feedback).
- Lookup answers are cached in `~/.cache/mzdesk`. An answer is reused for 30 days, then deleted.

## Internet search switch

The switch is in the sidebar footer of the local version.

- **On**: online lookups run without asking.
- **Off**: a run that needs the internet first opens a dialog that lists what it would send. **Continue** runs it this once, and the switch stays off. **Cancel** stops the run and keeps your inputs.
- The choice is remembered in this browser.

mzdesk.com has no switch. Your inputs reach the server anyway, so its lookups always run.

## What goes online

| Input or page | What is sent | Sent to |
|---|---|---|
| A compound in **Name / ID** mode, on any page | The name or ID | PubChem, UniChem |
| A spectrum from a **USI** | The USI | GNPS / MassIVE USI service |
| Candidate names in **EI MS Annotation** | Each line that is not a SMILES | PubChem, UniChem |
| [Identity](/docs/identity) | The compound's InChIKey, SMILES and IDs | PubChem, UniChem, LOTUS, NPClassifier, ClassyFire, LIPID MAPS, KEGG, Reactome |
| [Targeted Method](/docs/targeted_method) | The compound's InChIKey, or its SMILES | PubChem, UniChem |
| [Mass Search](/docs/ms1_search) | The mass window | PubChem (NCBI E-utilities and PUG REST) |
| [Alternative Explanations](/docs/exclusivity) | The formula and the mass | PubChem |
| [Spectrum in Which Set](/docs/gnps_fastmasst) | The spectrum and the precursor m/z | GNPS2 fastMASST; hits are checked through the GNPS USI service |
| **Deep-learning scores** in Similarity and Library Search | Nothing; the models are downloaded once | Zenodo |

- SMILES, InChI and formulas are read on your computer.
- Targeted Method falls back to values computed on your computer when the lookup fails.
- [Dataset Finder](/docs/gnps_redu) reads a copy of the Pan-ReDU metadata that mzDesk downloads from GNPS2 once, on first use; your search terms are not sent. [Fragment Pattern](/docs/gnps_massql) reads your file on mzDesk and does not send it to GNPS2.
- Every other page computes on mzDesk alone.

## First-use downloads

The local version downloads engines, libraries and models the first time a page needs them. See [First-use downloads](/docs/local-version#first-use-downloads). The files come from the sites of their projects, such as GitHub, Zenodo, GNPS and MassBank. These requests fetch files and carry none of your data.

## Version check

The local version asks mzdesk.com for the newest version number once a day, by reading `https://mzdesk.com/version.json`. The request carries no inputs and no ID. mzdesk.com sees your IP address, as with any web request.

- When a newer version exists, the sidebar shows **Version X is available**. **Later** hides it until the next version.
- To turn the check off, set `MZDESK_NO_UPDATE_CHECK=1`. See [Command Line and Settings](/docs/cli#settings).
