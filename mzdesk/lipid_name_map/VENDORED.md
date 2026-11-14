# lipid_name_map (vendored)

- What: the lipid name normalization package of the LIPID+ project, version 0.4.0 (`__version__`), copied in on 2026-10-31.
- Changes: none. The folder is the source package as released, with its `data/` files; only `__pycache__` / `.pyc` files are left out.
  Submodules use relative imports and find their data with `Path(__file__).parent / "data"`, so the package works under `mzdesk.lipid_name_map`.
- Python: standard library only. The source targets Python >= 3.12; the code uses no syntax newer than 3.10, and it was tested on 3.10 and 3.12 with the same results.
- Use: `from mzdesk.lipid_name_map import normalize`. The Lipid Short Name page (`mzdesk/api/runners/lipid.py`) and the illustrative structures (`mzdesk/lipid/`) use it.
- Packaging: the package data in `pyproject.toml` includes `lipid_name_map/data/**/*`.
- Upgrade: copy the new source folder over this one, update the version and date above, and run `tests/test_lipid_short_name.py`
  (its version check and golden comparison flag changes; new reason codes or flags need texts in `mzdesk/lipid/texts.py`).
