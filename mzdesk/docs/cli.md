# Command Line and Settings

This page lists the `mzdesk` commands of the local version and the environment variables that change how it runs.

## Start

```bash
mzdesk
```

mzDesk starts and opens `http://localhost:8000` in your browser. If port 8000 is busy, it takes the next free port up to 8020 and prints the address. Press `Ctrl+C` in the terminal to stop it.

| Option | Effect |
|---|---|
| `--port 8599` | Use this port. |
| `--host 0.0.0.0` | Listen on every network address, so other computers can open mzDesk. The default, `127.0.0.1`, is this computer only. |
| `--reload` | Restart when the code changes. For development. |

## Download components

```bash
mzdesk setup <component> ...
```

This downloads components ahead of time, for example before you work offline. The ids are `libraries`, `msfinder`, `frigid`, `cfmid`, `neims`, `molgpka`, `ccs` and `repo`. For sizes and the pages that use them, see [First-use downloads](/docs/local-version#first-use-downloads).

```bash
mzdesk setup libraries msfinder
```

## Update

```bash
mzdesk update
```

This runs the installer again and installs the newest version. `mzdesk update --check` only prints the newest version.

To install a given version, run the installer with `MZDESK_VERSION`:

```bash
curl -fsSL https://mzdesk.com/install.sh | MZDESK_VERSION=0.5.4 sh
```

On Windows, in PowerShell:

```powershell
$env:MZDESK_VERSION = "0.5.4"; irm https://mzdesk.com/install.ps1 | iex
```

## Settings

A setting is an environment variable. Set it when you start mzDesk:

```bash
MZDESK_FRIGID_DEVICE=cpu mzdesk
```

On Windows, in PowerShell:

```powershell
$env:MZDESK_FRIGID_DEVICE = "cpu"; mzdesk
```

To keep a setting, add a line such as `export MZDESK_FRIGID_DEVICE=cpu` to your shell profile.

### Folders and files

| Variable | What | Default |
|---|---|---|
| `MZDESK_LIBRARY_DIR` | MassBank and GNPS spectra | `~/.mzdesk/data/libraries` |
| `MZDESK_CCS_DIR` | CCSbase data and model | `~/.mzdesk/data/ccs` |
| `MZDESK_METHODS_DIR` | Known Methods data | `~/.mzdesk/data/methods` |
| `MZDESK_REPORT_DIR` | RepoRT retention data | `report/` in the Known Methods folder |
| `MZDESK_NEIMS_DIR` | NEIMS weights | `~/.mzdesk/models/neims` |
| `MZDESK_NEIMS_WEIGHTS` | The NEIMS weights file | `neims_weights.npz` in the NEIMS folder |
| `MZDESK_MOLGPKA_DIR` | MolGpKa weights | `~/.mzdesk/models/molgpka` |
| `MZDESK_MOLGPKA_WEIGHTS` | The MolGpKa weights file | `molgpka_weights.npz` in the MolGpKa folder |
| `MZDESK_MSFINDER_DIR` | MS-FINDER | `~/.mzdesk/tools/msfinder/MSFINDER_v3.73` |
| `MZDESK_MONO` | The `mono` program that runs MS-FINDER | `~/.mzdesk/tools/mono`, then the one on your PATH |
| `MZDESK_FRIGID_DIR` | FRIGID | `~/.mzdesk/tools/frigid` |
| `MZDESK_CACHE_DIR` | Lookup answers and GNPS working files | `~/.cache/mzdesk` for lookups, `~/.mzdesk/cache` for GNPS |
| `MZDESK_GNPS_REDU_CACHE` | The Pan-ReDU metadata snapshot of Dataset Finder | `gnps/redu_slim.tsv.gz` in `~/.mzdesk/cache` if it exists, else the snapshot that comes with mzDesk |
| `MZDESK_HMDB_XML` | An HMDB metabolites XML file you downloaded. Identity then adds biofluid concentrations. | None |
| `MZDESK_STATS_DB` | The usage counts | `~/.mzdesk/stats.sqlite` |
| `MZDESK_FEEDBACK_DIR` | Crash reports | `~/.mzdesk/feedback` |

### Engines

| Variable | Effect | Default |
|---|---|---|
| `MZDESK_FRIGID_DEVICE` | Where FRIGID computes: `auto`, `mps` (the Apple GPU) or `cpu`. `auto` takes the Apple GPU when enough memory is free. | `auto` |
| `MZDESK_FRIGID_THREADS` | CPU threads for FRIGID. | The cores that are free, at least 2 |
| `MZDESK_ANNOTATION_FRIGID` | `off` leaves **Run FRIGID** in MS/MS Annotation unticked, even when FRIGID is installed. | Ticked once FRIGID is installed |
| `MZDESK_CFMID_MEMORY` | Memory limit of the CFM-ID container, for example `4g`. | No limit |
| `MZDESK_CFMID_WORKDIR` | Folder for the files of a CFM-ID run. Docker must be able to see it. | The system temporary folder |
| `MZDESK_USI_ENDPOINT` | The service that turns a USI into a spectrum. | The GNPS USI service |

### Privacy

| Variable | Effect |
|---|---|
| `MZDESK_STATS=0` | Do not record usage counts. |
| `MZDESK_NO_UPDATE_CHECK=1` | Do not check for a new version. See [Version check](/docs/privacy#version-check). |

### Installer

These are read by the install command, not by `mzdesk`.

| Variable | Effect | Default |
|---|---|---|
| `MZDESK_VERSION` | Install this version instead of the newest. | The newest |
| `MZDESK_BASE_URL` | Where the installer reads `version.json` and the package. | `https://mzdesk.com` |
