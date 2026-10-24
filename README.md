# mzDesk

An integrated Streamlit app for small-molecule annotation, unified under a single
**neo-brutalist × minimalist, cool-white + indigo** design system.

## Tabs

| Tab | Status | What it does |
|-----|--------|--------------|
| **Basic** | ✅ live | One framed query panel — examples picker + SMILES beside the live 2D/3D structure, with the classifier thresholds underneath — then the results as full-width bands: the **lipidome/metabolome** classification (logP, acyl chains, isoprenoid/sterol scaffold, polar head groups, glycosyl rings, MW detection window) with substructure highlights on the 2D depiction, the rule-by-rule "why", and the **basic properties** (InChI Key, InChI, formula, molecular descriptors, ESI+ / ESI− adduct m/z tables, copy buttons for identifiers and m/z values). |
| **MS² Interpretation** | ✅ live | Annotate an MS/MS spectrum with **MetFrag** in-silico fragments. Interactive stick spectrum with m/z labels and the theoretical precursor marked — hover any peak for its ion formula, mass error, neutral loss and fragment structure; sortable annotation table. Adducts: [M+H]+, [M+Na]+, [M+K]+, [M+NH4]+, [M]+ and [M−H]−, [M+Cl]−, [M+HCOO]−, [M+CH3COO]−; a salt / mixture SMILES is annotated by its largest fragment. |
| **MS² Similarity** | ✅ live | Compare two MS/MS spectra across classic peak-matching, spectral-entropy, custom and (optionally) deep-learning metrics. The headline reports the weighted spectral-entropy similarity plus how many metrics pass the score threshold; the interactive **mirror MS² plot** marks peaks matched within tolerance (same greedy matching as matchms CosineGreedy) and labels the strongest m/z. |
| **Structure Annotation** | ✅ live | Molecular formula annotation with **BUDDY** (msbuddy; Xing et al., Nat. Methods 2023), plus the optional **SIRIUS 6** CLI when installed and signed in (`sirius login`). Required: precursor m/z + ion mode; optional: MS² spectrum, MS¹ isotope pattern (ignored when it does not contain the precursor peak), adduct, or a known formula (skips ranking; checked against the precursor m/z). Engines are picked automatically — SIRIUS only runs when an MS¹ isotope pattern and/or MS² spectrum is provided. BUDDY returns up to 5 candidate formulas (ML probability + estimated FDR) with one-click switching and explains every MS² peak by its fragment formula on an interactive spectrum. A caffeine [M+H]+ demo is pre-filled. |
| **Substructure Search** | 🚧 under construction | Placeholder. |

## Run

```bash
python3 -m streamlit run mzdesk/app.py   # from the repository root
```

Then open the URL Streamlit prints (default http://localhost:8501).

## Requirements

```bash
pip install -r requirements.txt
```

Plus:

- **Java 8+** on `PATH` — the MS² Interpretation tab shells out to the bundled
  `mzdesk/lib/MetFragCommandLine-2.5.2.jar`.
- *Optional:* the **SIRIUS 6** CLI for the Structure Annotation tab. It is looked up via
  `$SIRIUS_EXE`, then `PATH`, then the usual macOS app-bundle locations. SIRIUS 6 needs a
  signed-in SIRIUS account (`sirius login`); if the CLI cannot run, the tab shows a warning
  and keeps the BUDDY result.
- *Optional:* deep-learning scores in the MS² Similarity tab need `ms2deepscore` (+ torch)
  and `gensim`; the pretrained MS2DeepScore / Spec2Vec models are downloaded and cached on
  first use.

## Layout

```
mzDesk/
├── mzdesk/                           # application source
│   ├── app.py                        # entry: page config, global CSS, hero, st.tabs wiring
│   ├── theme.py                      # unified design system (CSS tokens + render helpers)
│   ├── tabs/
│   │   ├── tab_classifier.py         # Basic tab: properties + classifier UI
│   │   ├── tab_ms2.py                # MS² Interpretation: MetFrag annotation UI
│   │   ├── tab_ms2_similarity.py     # MS² Similarity: computation + results UI
│   │   ├── _ms2_mirror.py            # interactive mirror MS² plot helper
│   │   ├── tab_structure.py          # Structure Annotation: BUDDY + SIRIUS UI
│   │   └── tab_placeholders.py       # Substructure Search (placeholder)
│   ├── core/
│   │   ├── metfrag_engine.py         # MetFrag subprocess + spectrum/annotation logic
│   │   ├── ms2_similarity.py         # MS² similarity metrics (classic / entropy / custom / deep)
│   │   ├── msbuddy_engine.py         # BUDDY formula annotation + MS² fragment explanation
│   │   ├── sirius_engine.py          # optional SIRIUS 6 CLI wrapper
│   │   ├── molecular_properties.py   # InChI / formula / descriptors + adduct exact m/z
│   │   ├── lipid_metabolite_classifier.py
│   │   ├── lipid_metabolite_viz.py   # RDKit substructure-highlight rendering
│   │   └── display_en.py             # English display layer for the classifier
│   └── lib/
│       ├── 3Dmol-min.js              # bundled 3Dmol.js for the structure viewer
│       └── MetFragCommandLine-2.5.2.jar
├── tests/                            # unit tests (all modules) + equivalence tests (evidence modules)
├── examples/                         # bundled example spectra / queries used by the tutorial
├── validation/                       # benchmarks for the paper: equivalence, integration gain, case studies
├── paper/
│   ├── manuscript/
│   ├── figures/
│   └── literature/
├── admin/                            # internal docs (Chinese)
│   ├── README.md                     # project plan: goals, six-question protocol, scope, target journals
│   ├── version_change.html           # version history page
│   └── notes/
├── requirements.txt
└── .streamlit/config.toml            # Streamlit base theme matching mzdesk/theme.py
```

`mzdesk/app.py` puts `mzdesk/` and `mzdesk/core/` on `sys.path`, so tab modules use flat
imports (`import theme`, `import metfrag_engine as engine`, `from lipid_metabolite_classifier import ...`).

## Design system (`mzdesk/theme.py`)

Thick navy-ink borders, hard offset shadows, near-square corners, uppercase mono
labels, a restrained indigo palette. Two things sit on top of that look and are
deliberately not decoration:

**Two heading levels.** `section()` is a page section (indigo square marker + label);
`subsection()` is a block inside one (no marker, hairline rule, lighter). They are not
interchangeable — using one level for both is what made a card title read with exactly
the same weight as a page section.

**Motion.** Switching tabs, and a result arriving after a run, are the two moments a
reader can lose their place, so both are animated (`nb-rise` / `nb-fade`, ~0.3 s). A
warning banner also gets one pass of the bright accent through its shadow so a stale
result is noticed. Nothing else moves; the button press stays instant on purpose.
Everything collapses under `prefers-reduced-motion: reduce`.

`st.tabs` switches panels on the client, so a revealed panel is not rebuilt and a plain
CSS animation would never replay — `MOTION_JS` (injected by `inject_css` through
`st.html`, since `st.markdown` strips `<script>`) re-adds the `.nb-enter` class on each
switch to retrigger it.

- Colour tokens: `INK`, `BG`, `SURFACE`, `SURFACE_2`, `MUTED`, `LINE`, `BLUE`, `BLUE_BRIGHT`,
  `BLUE_SOFT`, `BLUE_GHOST`, plus `CATEGORY_ACCENT` / `PLOT_DARK` / `PLOT_LIGHT` / `FRAG_HL`.
- Render helpers: `inject_css`, `hero`, `tab_header`, `section`, `subsection`, `panel`,
  `metric_cards`, `stat_grid`, `verdict`, `chips`, `chip_line`, `reasons`, `notes`,
  `under_construction`, `divider`, `note`, and the tab-specific `copy_fields` /
  `molecule_3d` (Basic) and `similarity_score_cards` (MS² Similarity). The module
  docstring lists the full API.
- **Filling a full-width band.** A section whose body is one narrow column of text leaves
  more than half the canvas empty, so `reasons()` (bulleted conclusions, last entry = the
  verdict) and `notes()` (muted paragraphs) flow their items through CSS multi-column and
  let the browser balance the columns. Both collapse to one column under 820 px.
- **Stat grids.** `.nb-statgrid` sets `grid-auto-rows:1fr`, so every row takes the height of
  the tallest one — otherwise a card whose label wraps ("Isoprenoid / sterol") leaves its
  neighbours short and the grid reads as misaligned. `stat_grid(..., wide=True, max_cols=5)`
  gives ten cards a 5 + 5 layout instead of a ragged 4 + 4 + 2; `min_px` lowers the
  auto-fit floor when the container is narrow.
- Bordered input panels: `with theme.panel("setup"): ...` — a keyed
  `st.container(border=True)`; Streamlit 1.5x gives bordered containers no CSS hook of
  their own, so the key class is what `GLOBAL_CSS` styles. Every tab frames its query in
  one, with the run button on the last row beside the advanced options.
- `.streamlit/config.toml` holds the matching Streamlit base theme (colours, 2px
  `baseRadius`, indigo `st.info`); `GLOBAL_CSS` adds the ink borders and hard shadows.

**Peak-list heights.** `st.text_area(height=H)` gives the textarea `H - 30` px (Streamlit's
label strip); with the 20 px line box `GLOBAL_CSS` sets, a whole number of rows is
`H = rows * 20 + 39`. Use that formula for any new peak-list input.

The `+39` (not `+48`) matters: the textarea's 9 px top and bottom padding are *inside* the
box, and Chrome lets a line's ascent render into the bottom padding. With `+48` the box
shows `rows + 0.45` lines, so the next row is sliced horizontally across the middle — a
peak list with 13 rows in an 8-row box visibly clipped its 9th line.

**Adduct rows are a grid, not flex.** An m/z is one unbreakable token, so in a narrow card
a flex row lets it overflow on top of the Δ note instead of shrinking; grid tracks cannot
collide. Below 250px of card width the Δ note is dropped rather than crushed.

To extend a tab, keep all styling in `theme.py` — no per-tab CSS.

## Provenance

Integrated from two existing tools:
- `lipid_or_metabolite/` → classifier core + visualization.
- `formula_assign/metfrag_ms2_app.py` → MetFrag engine (refactored) + redesigned UI.

Tab/engine modules were drafted by a `deepseek-v4-flash` executor against written specs,
then reviewed, corrected and integrated. The `theme.py` design system and `app.py` shell
were authored directly.
