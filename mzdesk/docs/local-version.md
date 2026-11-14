# Local Version

The local version runs mzDesk on your own computer and opens it in your browser. This page says how to install, update and remove it.

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

mzDesk starts and opens `http://localhost:8731` in your browser. If port 8731 is busy, it takes the next free port up to 8750 and prints the address. Press `Ctrl+C` in the terminal to stop it. If `mzdesk` is not found, open a new terminal window.

## First-use downloads

The first time you run a page that needs a missing component, mzDesk downloads it into `~/.mzdesk`. A dialog lists what it needs, with sizes. **Download** starts it. **Not now** cancels the run and keeps your inputs. Small components download without asking and show a progress bar.

| Component | Used by | Download | Asks first | Systems |
|---|---|---|---|---|
| `libraries`: MassBank and GNPS spectra | Library Search, MS/MS Annotation, EI MS Annotation | 2.3 GB | Yes, except inside MS/MS Annotation | macOS, Linux, Windows |
| `msfinder`: MS-FINDER (with Mono on macOS and Linux) | Structure Prediction, MS/MS Annotation | 0.5 GB | Yes, except inside MS/MS Annotation | macOS, Linux, Windows |
| `frigid`: FRIGID model | De Novo Structure; MS/MS Annotation with **Run FRIGID** | 3 GB | Always | macOS, Linux, Windows |
| `cfmid`: CFM-ID 4 and MSRB, a Docker image | MS/MS Prediction | 0.2 GB (about 1 GB on disk) | Yes | macOS, Linux, Windows |
| `neims`: NEIMS weights | EI-MS Prediction, EI MS Annotation | about 1 GB | Yes | macOS, Linux, Windows |
| `molgpka`: MolGpKa weights | Separation, Ionization, Adducts, Known Methods, Targeted Method, MS/MS Annotation | 0.2 GB | No | macOS, Linux, Windows |
| `ccs`: CCSbase data and model | CCS; MS/MS Annotation when you give an observed CCS | 8 MB | No | macOS, Linux, Windows |
| `repo`: RepoRT retention data | Known Methods, Targeted Method | 0.2 GB | No | macOS, Linux, Windows |
| `redu`: Pan-ReDU sample metadata from GNPS2 | Dataset Finder; Spectrum in Which Set with **Join Pan-ReDU** | 1.7 GB (8 MB kept on disk) | Yes | macOS, Linux, Windows |

- On Windows, `msfinder` needs no Mono: it runs MS-FINDER with the .NET Framework 4.7.2 or newer that comes with Windows 10 (version 1803 and later) and Windows 11. `cfmid` needs Docker Desktop in Linux-container mode. `frigid` builds a 6 GB Python environment; if its setup stops on a very long file path, turn on long paths in Windows or install mzDesk closer to the drive root.
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
