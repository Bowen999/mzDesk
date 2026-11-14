# Files for the Windows install of MS-FINDER

MS-FINDER v3.73 (RIKEN, Tsugawa lab) is released for Windows, so on Windows mzDesk runs it directly (no Mono). Three small files are added to the
official release zip. Both were made on macOS / Linux by `python -m mzdesk.ms2.msfinder_setup`, the same way as on those systems, and are checked
by SHA-256 before they are used (`mzdesk/ms2/msfinder_setup.py`).

| File | What it is | SHA-256 |
|---|---|---|
| `MsfinderConsoleApp.exe` | Console front end compiled from the 9 unmodified source files of `src/MSFINDER/MsfinderConsoleAppCore` at tag `MSFINDER-v3.73` (LGPL-3.0) | `15d38b986b4139c8e7398bbf25f3afbec847525c5af71b7c4c0620ed6714f931` |
| `StructureFinder.fixed.dll` | The release `StructureFinder.dll` with ONE change: the field `FragmentPeakMatcher.maxBondDissociationEnergy` is marked `[ThreadStatic]` (`mzdesk/lib/msfinder/ThreadStaticPatch.cs`), so that identical input gives identical structure scores on every run | `2bd9f959afb9457f40ad6911b67fa29ac58d9aee07745381fa43ece4712d9820` |
| `MsfinderConsoleApp.exe.config` | Binding redirects written by mzDesk (no MS-FINDER code): the .NET Framework otherwise refuses the newer `System.Threading.Tasks.Extensions` and `System.Buffers` builds that the release ships. Copied next to the exe in each run's private copy | not pinned (plain text) |

The assemblies of MS-FINDER are CC BY 4.0, the source LGPL-3.0, from github.com/systemsomicslab/MsdialWorkbench. `StructureFinder.fixed.dll` is a
modified copy of the release file (the change is described above). Citation: Tsugawa H. et al. Anal. Chem. 2016, 88, 7946-7958 (MS-FINDER).
