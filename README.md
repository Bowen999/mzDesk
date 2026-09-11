# Metabolomics Workbench

An integrated Streamlit app for small-molecule annotation, unified under a single
**neo-brutalist × minimalist, cool-white + indigo** design system.

## Tabs

| Tab | Status | What it does |
|-----|--------|--------------|
| **Basic** | ✅ live | Two output panels: **basic properties** (interactive 3D structure, InChI Key, InChI, formula, molecular descriptors, exact m/z for common adducts via a dropdown, copy buttons for long identifiers) and the **lipidome/metabolome** classification (logP, acyl chains, isoprenoid/sterol scaffold, polar head groups, glycosyl rings) with substructure highlights on the 2D depiction. |
| **MS² Interpretation** | ✅ live | Annotate an MS/MS spectrum with **MetFrag** in-silico fragments. Interactive spectrum — hover an annotated peak for m/z, neutral loss, ion formula and the fragment structure. |
| **MS2 Similarity** | ✅ live | Compare two MS/MS spectra across classic peak-matching, spectral-entropy, custom and (optionally) deep-learning metrics, with an interactive **Mirror MS² plot** that highlights matched peaks within tolerance. |
| **Structure Annotation** | ✅ live | Dual-engine molecular formula annotation: **BUDDY** (always available) + optional **SIRIUS 6** CLI if installed. Engines are picked automatically from the inputs — SIRIUS only runs when an MS¹ isotope pattern and/or MS² spectrum is provided. BUDDY gives Top-5 formulas with ML probability + estimated FDR; per-peak MS² fragment formulas on an interactive spectrum.
| *(was)* **BUDDY** molecular formula annotation (Xing et al., Nat. Methods 2023). Required: precursor m/z + ion mode; optional: MS², MS¹ isotope pattern, retention time, adduct, or a known formula (skips ranking). Returns Top-5 candidate formulas (ML probability + estimated FDR) with one-click switching, and explains every MS² peak by its fragment formula on an interactive spectrum. |
| **Substructure Search** | 🚧 under construction | Placeholder. |

## Run

```bash
cd "~/Documents/实用小工具/metabolite"
python3 -m streamlit run app.py
```

Then open the URL Streamlit prints (default http://localhost:8501).

## Requirements

```bash
pip install -r requirements.txt
```

Plus a **Java 8+** runtime on `PATH` — the MS² tab shells out to the bundled
`lib/MetFragCommandLine-2.5.2.jar`.

## Layout

```
metabolite/
├── app.py                 # entry: page config, global CSS, hero, st.tabs wiring
├── theme.py               # unified design system (CSS tokens + render helpers)
├── tabs/
│   ├── tab_classifier.py  # Basic tab: properties + classifier UI
│   ├── tab_ms2.py         # MS² MetFrag annotation UI (redesigned)
│   ├── tab_ms2_similarity.py  # MS2 similarity computation + results UI
│   ├── _ms2_mirror.py     # interactive Mirror MS² plot helper
│   └── tab_structure.py   # Structure Annotation (BUDDY + SIRIUS)
│   └── tab_placeholders.py# Substructure Search (WIP)
├── core/
│   ├── metfrag_engine.py            # MetFrag subprocess + spectrum/annotation logic
│   ├── ms2_similarity.py            # 14 MS² similarity metrics (classic/entropy/custom/deep)
│   ├── molecular_properties.py      # InChI / formula / descriptors + adduct exact m/z
│   ├── lipid_metabolite_classifier.py
│   ├── lipid_metabolite_viz.py      # RDKit substructure-highlight rendering
│   └── display_en.py                # English display layer for the classifier
├── lib/
│   ├── 3Dmol-min.js                    # bundled 3Dmol.js for the structure viewer
│   └── MetFragCommandLine-2.5.2.jar
└── .streamlit/config.toml           # light indigo base theme
```

`app.py` puts the project root and `core/` on `sys.path`, so tab modules use flat
imports (`import theme`, `import metfrag_engine as engine`, `from lipid_metabolite_classifier import ...`).

## Design system (`theme.py`)

- Colour tokens: `INK`, `BG`, `SURFACE`, `MUTED`, `BLUE`, `BLUE_SOFT`, … and
  `CATEGORY_ACCENT` / `PLOT_DARK` / `PLOT_LIGHT` / `FRAG_HL`.
- Render helpers: `inject_css`, `hero`, `tab_header`, `section`, `metric_cards`,
  `stat_grid`, `verdict`, `chips`, `under_construction`, `divider`, `note`.

To extend a tab, keep all styling in `theme.py` — no per-tab CSS.

## Provenance

Integrated from two existing tools:
- `lipid_or_metabolite/` → classifier core + visualization.
- `formula_assign/metfrag_ms2_app.py` → MetFrag engine (refactored) + redesigned UI.

Tab/engine modules were drafted by a `deepseek-v4-flash` executor against written specs,
then reviewed, corrected and integrated. The `theme.py` design system and `app.py` shell
were authored directly.
