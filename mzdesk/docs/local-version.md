# Local Version

The local version runs mzDesk on your own computer and opens it in your browser. This page says how it differs from mzdesk.com, and how to install, update and remove it.

## Local or mzdesk.com?

![Where each version computes](fig:local-vs-web)

The pages, the results and the Docs are the same. These things differ:

| | mzdesk.com | Local version |
|---|---|---|
| Install | None. Open it in a browser. | One command, about 1.6 GB. |
| Where it computes | On the mzDesk server. | On your computer. |
| Your inputs | Sent to the server over HTTPS. | Stay on your computer. Only online lookups leave it. |
| Internet search switch | None. Online lookups always run. | In the sidebar. Off asks before a run that needs the internet. |
| Engines and libraries | All installed. | Downloaded on first use. Large ones ask first. |
| Run FRIGID (MS/MS Annotation, Advanced) | Off. | On once FRIGID is installed. |
| Speed | Shared with other users. A De Novo Structure run takes about two minutes. | Your hardware. On an Apple-silicon Mac, about 20 to 40 s. |
| CFM-ID (MS/MS Prediction) | Ready. | You install Docker Desktop. |
| Updates | Automatic. | `mzdesk update`. |

Use the local version when your data must stay on your computer. For everything that leaves it, see [Privacy and Internet Search](/docs/privacy).

## Install

Run one command in a terminal.

macOS and Linux:

```bash
curl -fsSL https://mzdesk.com/install.sh | sh
```

Windows, in PowerShell:

```powershell
irm https://mzdesk.com/install.ps1 | iex
```

- It needs no `sudo` and no administrator rights, an internet connection and about 1.6 GB of disk space.
- It ends with `Run: mzdesk`. If a step fails, it names the step and gives one sentence of reason.

![What the installer does](fig:install-flow)

Supported systems:

- macOS and Linux, on x86_64 or arm64. On Apple Silicon you need macOS 14 or newer. Linux needs glibc (for example Ubuntu or Debian), not musl (Alpine); on arm64, glibc 2.34 or newer.
- Windows, 64-bit Intel or AMD, with PowerShell 5.1 or 7. Windows has the core, the libraries, CCS and MolGpKa only. A workflow skips an engine it does not have and says so in the step note.

> **Note** The Windows installer is not verified on Windows yet.

### What gets installed

Everything goes into one folder, `~/.mzdesk` (on Windows, `.mzdesk` in your user folder).

| Path | Content |
|---|---|
| `venv/` | A private Python 3.12 and the mzDesk packages. Your system Python is not used. |
| `bin/` | The `mzdesk` launcher. Also `uv`, and `micromamba` when an engine needs it. |
| `tools/` | Engines installed on first use: `msfinder/`, `frigid/`, `mono/`. |
| `data/`, `models/`, `cache/` | Libraries, model weights and caches. |

The installer puts `~/.mzdesk/bin` on your PATH: a link in `~/.local/bin` if that is on it, otherwise one `export PATH` line in your shell profile. On Windows it edits your user PATH. A Java runtime comes with the install; mzDesk uses the `java` on your PATH if there is one.

## First run

```bash
mzdesk
```

mzDesk starts and opens `http://localhost:8000` in your browser. If port 8000 is busy, it takes the next free port up to 8020 and prints the address. Press `Ctrl+C` in the terminal to stop it. If `mzdesk` is not found, open a new terminal window.

## First-use downloads

The first time you run a page that needs a missing component, mzDesk downloads it into `~/.mzdesk`. A dialog lists what it needs, with sizes. **Download** starts it. **Not now** cancels the run and keeps your inputs. Small components download without asking and show a progress bar.

| Component | Used by | Download | Asks first | Systems |
|---|---|---|---|---|
| `libraries`: MassBank and GNPS spectra | Library Search, MS/MS Annotation, EI MS Annotation | 2.3 GB | Yes, except inside MS/MS Annotation | macOS, Linux, Windows |
| `msfinder`: MS-FINDER with Mono | Structure Prediction, MS/MS Annotation | 0.5 GB | Yes, except inside MS/MS Annotation | macOS, Linux |
| `frigid`: FRIGID model | De Novo Structure; MS/MS Annotation with **Run FRIGID** | 3 GB | Always | macOS, Linux |
| `cfmid`: CFM-ID 4 and MSRB, a Docker image | MS/MS Prediction | 0.2 GB (about 1 GB on disk) | Yes | macOS, Linux |
| `neims`: NEIMS weights | EI-MS Prediction, EI MS Annotation | about 1 GB | Yes | macOS, Linux |
| `molgpka`: MolGpKa weights | Separation, Ionization, Adducts, Known Methods, Targeted Method, MS/MS Annotation | 0.2 GB | No | macOS, Linux, Windows |
| `ccs`: CCSbase data and model | CCS; MS/MS Annotation when you give an observed CCS | 8 MB | No | macOS, Linux, Windows |
| `repo`: RepoRT retention data | Known Methods, Targeted Method | 0.2 GB | No | macOS, Linux |

- `neims` and `molgpka` convert their weights once with TensorFlow or PyTorch in a temporary environment that is deleted afterwards. That is most of their download size.
- `cfmid` needs Docker Desktop. mzDesk cannot install it. If Docker is missing or stopped, the page says so and links to [docs.docker.com/get-docker](https://docs.docker.com/get-docker/).
- The MS2DeepScore (416 MB) and Spec2Vec (about 280 MB) models download themselves the first time you tick **Deep-learning scores**.
- To download ahead of time, use [`mzdesk setup`](/docs/cli#download-components).

## Update

```bash
mzdesk update
```

This runs the installer again. When a newer version exists, the sidebar shows **Version X is available**. See [Version check](/docs/privacy#version-check).

## Uninstall

Delete two folders: `~/.mzdesk` holds the program, the downloaded components and your local data; `~/.cache/mzdesk` holds the cached lookup answers.

```bash
rm -rf ~/.mzdesk ~/.cache/mzdesk
```

On Windows:

```powershell
Remove-Item -Recurse -Force $HOME\.mzdesk, $HOME\.cache\mzdesk
```

The PATH entry stays. Remove the `~/.local/bin/mzdesk` link or the `export PATH` line in your shell profile if you want it gone. On Windows, remove `.mzdesk\bin` from your user PATH.
