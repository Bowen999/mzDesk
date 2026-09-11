import streamlit as st

import theme
import molecular_properties
from lipid_metabolite_classifier import analyze, Thresholds
from lipid_metabolite_viz import draw_structure_highlight, HIGHLIGHT_COLORS
from display_en import category_en, category_subtitle_en, feature_rows_en, english_reasons


EXAMPLES = {
    "— pick an example —": "",
    "Cholesterol (lipid)": "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
    "Palmitic acid (lipid)": "CCCCCCCCCCCCCCCC(=O)O",
    "POPC phospholipid (lipid)": "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
    "Retinol / Vitamin A (lipid)": "OC/C=C(C)/C=C/C=C(C)/C=C/C1=C(C)CCCC1(C)C",
    "Glucose (metabolite)": "OCC1OC(O)C(O)C(O)C1O",
    "Caffeine (metabolite)": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Citric acid (metabolite)": "OC(=O)CC(O)(C(=O)O)CC(=O)O",
    "Ibuprofen (both)": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "Anthracene (both)": "c1ccc2cc3ccccc3cc2c1",
    "Water (neither)": "O",
}

LEGEND_LABELS = {
    "acyl": ">=C12 acyl chain",
    "isoprenoid": "Isoprenoid / sterol",
    "head": "Polar head group",
    "glycosyl": "Glycosyl ring",
}

LEGEND_KEYS = ["acyl", "isoprenoid", "head", "glycosyl"]


def _rgb(t):
    r, g, b = t
    return f"rgb({int(r*255)},{int(g*255)},{int(b*255)})"


def _on_example_change():
    sel = st.session_state.get("clf_example")
    if sel and sel != "— pick an example —":
        st.session_state["clf_smiles"] = EXAMPLES[sel]


def render():
    theme.tab_header(
        "Basic",
        subtitle="Chemical identity, molecular descriptors, adduct m/z and a "
                 "lipidome/metabolome classification from a SMILES string.",
        kicker="SMILES classifier",
    )

    if "clf_smiles" not in st.session_state:
        st.session_state["clf_smiles"] = EXAMPLES["Cholesterol (lipid)"]

    c1, c2, c3 = st.columns([1.1, 2, 0.7])
    c1.selectbox("Examples", list(EXAMPLES.keys()), key="clf_example", on_change=_on_example_change)
    c2.text_input("SMILES", key="clf_smiles")
    c3.markdown("<div style='height:1.75rem'></div>", unsafe_allow_html=True)
    c3.button("Analyze", key="clf_run")

    smiles = st.session_state["clf_smiles"].strip()
    if not smiles:
        st.info("Enter a SMILES string or choose an example.")
        return

    with st.expander("Thresholds (heuristic cut-offs)"):
        t1, t2, t3 = st.columns(3)
        with t1:
            mw_min = st.number_input("MW min", min_value=0.0, max_value=2000.0, value=50.0, step=10.0, key="clf_mw_min")
            acyl_min_carbons = st.number_input("Acyl min carbons", min_value=4, max_value=30, value=12, step=1, key="clf_acyl_min_carbons")
            medium_chain_max = st.number_input("Medium chain max", min_value=4, max_value=30, value=14, step=1, key="clf_medium_chain_max")
        with t2:
            mw_max = st.number_input("MW max", min_value=100.0, max_value=10000.0, value=2000.0, step=100.0, key="clf_mw_max")
            logp_polar_max = st.number_input("LogP polar max", min_value=-10.0, max_value=20.0, value=3.0, step=0.1, key="clf_logp_polar_max")
            boundary_logp_low = st.number_input("Boundary logP low", min_value=-10.0, max_value=20.0, value=1.0, step=0.1, key="clf_boundary_logp_low")
        with t3:
            boundary_logp_high = st.number_input("Boundary logP high", min_value=-10.0, max_value=20.0, value=3.0, step=0.1, key="clf_boundary_logp_high")
            boundary_mw_max = st.number_input("Boundary MW max", min_value=50.0, max_value=5000.0, value=400.0, step=10.0, key="clf_boundary_mw_max")
            isoprene_min_units = st.number_input("Isoprene min units", min_value=1, max_value=20, value=2, step=1, key="clf_isoprene_min_units")

        thr = Thresholds(
            mw_min=mw_min,
            mw_max=mw_max,
            acyl_min_carbons=int(acyl_min_carbons),
            logp_polar_max=logp_polar_max,
            medium_chain_max=int(medium_chain_max),
            boundary_logp_low=boundary_logp_low,
            boundary_logp_high=boundary_logp_high,
            boundary_mw_max=boundary_mw_max,
            isoprene_min_units=int(isoprene_min_units),
        )

    result = analyze(smiles, thr)
    props = molecular_properties.from_smiles(smiles)
    en = category_en(result.category)
    accent = theme.CATEGORY_ACCENT.get(en, theme.BLUE)

    if not props.valid or not result.features.valid:
        theme.section("Basic Properties")
        st.error("The SMILES could not be parsed by RDKit. Please check the string.")
        theme.divider()
        theme.section("Lipidome / Metabolome")
        theme.verdict(en, category_subtitle_en(en), accent)
        return

    # --- Section 1: Basic properties ------------------------------------
    theme.section("Basic Properties")

    viewer_col, info_col = st.columns([1.15, 1])
    with viewer_col:
        st.markdown(theme.section_label_html("3D Structure"), unsafe_allow_html=True)
        molblock = molecular_properties.structure_3d(smiles)
        if molblock:
            theme.molecule_3d(molblock, height=430)
        else:
            fallback_img = draw_structure_highlight(result, size=(600, 430))
            st.image(fallback_img, width="stretch")
            theme.note("3D embedding unavailable — showing the 2D depiction.")

    with info_col:
        theme.copy_fields(
            [
                ("InChI Key", props.inchikey or "—", True),
                ("Formula", props.formula or "—", False),
                ("Canonical SMILES", props.canonical_smiles or "—", True),
                ("InChI", props.inchi or "—", True),
            ],
            heading="Molecular Identity",
        )
        st.markdown(theme.section_label_html("Common Adducts"), unsafe_allow_html=True)
        adduct_label = st.selectbox(
            "Adduct", list(molecular_properties.ADDUCTS), key="clf_adduct"
        )
        adduct_delta, _ = molecular_properties.adduct_spec(adduct_label)
        adduct_mz = molecular_properties.adduct_exact_mz(props.exact_mass, adduct_label)
        theme.metric_cards([
            (f"{adduct_label} exact m/z", f"{adduct_mz:.4f}"),
            ("Adduct Δ (Da)", f"{adduct_delta:+.4f}"),
        ])
        theme.note("Exact m/z = (neutral monoisotopic mass + adduct Δ), adjusted for the adduct's charge state.")

    st.markdown(theme.section_label_html("Molecular Properties"), unsafe_allow_html=True)
    theme.stat_grid(molecular_properties.property_rows(props), wide=True)

    theme.divider()

    # --- Section 2: Lipidome / Metabolome classification -----------------
    theme.section("Lipidome / Metabolome")
    theme.verdict(en, category_subtitle_en(en), accent)

    left, right = st.columns([1.05, 1])

    with left:
        st.markdown(theme.section_label_html("Structure"), unsafe_allow_html=True)
        img = draw_structure_highlight(result, size=(640, 520))
        st.image(img, width="stretch")
        theme.chips([(LEGEND_LABELS[k], _rgb(HIGHLIGHT_COLORS[k])) for k in LEGEND_KEYS])

    with right:
        theme.section("Features")
        theme.stat_grid(feature_rows_en(result.features))

        theme.section("Detection window (MW)")

        def _window_bar(mw, lo, hi):
            denom = hi - lo
            pct = 0.0 if denom <= 0 else max(0.0, min(100.0, (mw - lo) / denom * 100.0))
            inside = lo <= mw <= hi
            mark = theme.BLUE if inside else "#C0392B"
            html = (
                f'<div style="margin-top:4px;">'
                f'<div style="position:relative;height:22px;border:3px solid {theme.INK};'
                f'background:{theme.SURFACE_2};border-radius:2px;">'
                f'<div style="position:absolute;top:-6px;left:calc({pct:.2f}% - 7px);width:13px;'
                f'height:32px;border:3px solid {theme.INK};background:{mark};border-radius:2px;"></div></div>'
                f'<div style="display:flex;justify-content:space-between;font-size:.7rem;'
                f'font-weight:700;margin-top:5px;color:{theme.MUTED};">'
                f'<span>{lo:.0f}</span><span style="color:{theme.INK};font-weight:800;">MW {mw:.1f}</span>'
                f'<span>{hi:.0f}</span></div></div>'
            )
            st.markdown(html, unsafe_allow_html=True)

        _window_bar(result.features.mw, float(thr.mw_min), float(thr.mw_max))

        theme.section("Why")
        for r in english_reasons(result, thr):
            st.markdown(f"- {r}")

    theme.divider()
    theme.note("Heuristic tool — lipid vs metabolite has no hard boundary. Tune thresholds above.")
