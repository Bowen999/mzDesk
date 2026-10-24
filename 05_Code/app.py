"""Metabolomics Workbench — integrated Streamlit app.

Five tabs under one neo-brutalist / minimalist blue design system:
  1. Basic                 — SMILES classifier + basic properties
                             (core/lipid_metabolite_*, core/molecular_properties)
  2. MS² Interpretation    — MetFrag fragment annotation (core/metfrag_engine)
  3. MS² Similarity        — multi-metric MS² comparison + interactive mirror plot
                             (core/ms2_similarity)
  4. Structure Annotation  — molecular formula annotation with BUDDY (core/msbuddy_engine),
                             plus the SIRIUS CLI when installed (core/sirius_engine)
  5. Substructure Search   — under construction

This shell owns the single set_page_config, the global CSS, the app hero, and the
tab wiring. Each tab module exposes a render function and never calls set_page_config.
"""
import sys
from pathlib import Path

import streamlit as st

# Make the app root and core/ importable with flat names, so the copied core
# modules (which import each other flatly) and the tab modules resolve cleanly.
ROOT = Path(__file__).parent.resolve()
CORE = ROOT / "core"
for _p in (str(ROOT), str(CORE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import theme  # noqa: E402
from tabs.tab_classifier import render as render_classifier  # noqa: E402
from tabs.tab_ms2 import render as render_ms2  # noqa: E402
from tabs.tab_ms2_similarity import render as render_ms2_similarity  # noqa: E402
from tabs.tab_structure import render as render_structure  # noqa: E402
from tabs.tab_placeholders import render_substructure  # noqa: E402


def main():
    st.set_page_config(
        page_title="Metabolomics Workbench",
        page_icon="🔷",
        layout="wide",
    )
    theme.inject_css()

    theme.hero(
        'Metabolomics <span class="accent">Workbench</span>',
        subtitle="An integrated bench for small-molecule annotation — classify lipids vs "
                 "metabolites, interpret MS² spectra with MetFrag, compare MS² spectra, and "
                 "annotate molecular formulas with BUDDY (plus SIRIUS when installed).",
        kicker="lipidomics · metabolomics · ms² · structure",
    )

    tab_clf, tab_ms2, tab_sim, tab_struct, tab_sub = st.tabs([
        "Basic",
        "MS² Interpretation",
        "MS² Similarity",
        "Structure Annotation",
        "Substructure Search",
    ])

    with tab_clf:
        render_classifier()
    with tab_ms2:
        render_ms2()
    with tab_sim:
        render_ms2_similarity()
    with tab_struct:
        render_structure()
    with tab_sub:
        render_substructure()


# `streamlit run 05_Code/app.py` executes the script with __name__ == "__main__".
if __name__ == "__main__":
    main()
