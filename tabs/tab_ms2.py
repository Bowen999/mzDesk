import streamlit as st
import pandas as pd
import streamlit.components.v1 as components
import theme
import metfrag_engine as engine


DEFAULT_SMILES = r"CCCCCCCC\C=C/CCCCCCCC(=O)OCC(COP(O)(=O)OCC(O)CO)OC(=O)CCCCCCC\C=C/CCCCCCCC"


def render():
    theme.tab_header(
        "MS² Interpretation",
        subtitle=(
            "Annotate an MS/MS spectrum with MetFrag in-silico fragment structures. "
            "Hover an annotated peak to see its ion formula, neutral loss and fragment."
        ),
        kicker="MetFrag · fragment annotation",
    )

    left, right = st.columns([1.1, 1])

    with left:
        smiles = st.text_input(
            "SMILES",
            value=DEFAULT_SMILES,
            key="ms2_smiles",
        )

        spectrum_text = st.text_area(
            "Spectrum (peak list)",
            value=engine.DEFAULT_SPECTRUM,
            height=180,
            help="`mz intensity` per line, OR space-separated `mz:intensity` pairs. Intensity optional.",
            key="ms2_spectrum",
        )

        c1, c2 = st.columns(2)
        adduct = c1.selectbox(
            "Adduct",
            list(engine.ADDUCTS.keys()),
            index=0,
            key="ms2_adduct",
        )
        precursor_input = c2.text_input(
            "Precursor m/z",
            value="",
            help="Blank = auto from neutral mass + adduct",
            key="ms2_precursor",
        )

        tree_depth = st.slider(
            "Tree depth",
            1,
            4,
            value=2,
            key="ms2_depth",
        )

        with st.expander("Advanced tolerances"):
            tol_abs = st.number_input(
                "Absolute mass tolerance (Da)",
                value=0.005,
                min_value=0.0,
                step=0.001,
                format="%.3f",
                key="ms2_tolabs",
            )
            tol_ppm = st.number_input(
                "Relative tolerance (ppm)",
                value=15.0,
                min_value=0.0,
                step=1.0,
                key="ms2_tolppm",
            )

        run = st.button("Run MetFrag", type="primary", key="ms2_run")

    with right:
        theme.section("Input structure")
        if smiles.strip():
            png = engine.render_input_png(smiles)
            if png:
                st.image(png, width="stretch")
            else:
                st.warning("Invalid SMILES — cannot draw structure.")
        else:
            theme.note("Enter a SMILES to preview the structure.")

    if run:
        try:
            peaks = engine.parse_spectrum(spectrum_text)
            if not peaks:
                st.error("No valid peaks found in the spectrum. Check the peak list format.")
                return

            precursor_val = None
            if precursor_input.strip():
                try:
                    precursor_val = float(precursor_input)
                except ValueError:
                    st.error("Precursor m/z must be a number or left blank.")
                    return

            with st.spinner("Running MetFrag…"):
                st.session_state["ms2_result"] = engine.run_metfrag(
                    smiles,
                    peaks,
                    adduct,
                    precursor_val,
                    tree_depth,
                    tol_abs,
                    tol_ppm,
                )
        except Exception as e:
            st.error(str(e))
            return

    if "ms2_result" in st.session_state:
        result = st.session_state["ms2_result"]

        theme.divider()

        if result.get("warning"):
            st.warning(result["warning"])

        annotation = result["annotation"]
        n_peaks = len(annotation)
        n_ann = sum(1 for a in annotation if a["annotated"])

        theme.section("Summary")
        theme.metric_cards([
            ("Peaks", str(n_peaks)),
            ("Annotated", str(n_ann)),
            ("Precursor m/z", f"{result['precursor_mz']:.4f}"),
            ("Formula", result["mol_formula"]),
        ])

        theme.chips([
            ("Annotated", theme.PLOT_DARK),
            ("Unannotated", theme.PLOT_LIGHT),
        ])
        theme.note(
            "Hover an annotated (dark) peak for m/z · neutral loss · ion formula · fragment structure."
        )

        theme.section("Spectrum")

        fig = engine.build_figure(
            result,
            dark=theme.PLOT_DARK,
            light=theme.PLOT_LIGHT,
        )
        html = (
            engine.HOVER_TEMPLATE
            .replace("__PLOTLYJS__", engine.get_plotlyjs())
            .replace("__FIGJSON__", fig.to_json())
        )
        components.html(html, height=500, scrolling=False)

        with st.expander("Annotation table"):
            rows = []
            for a in annotation:
                rows.append({
                    "m/z": a["mz"],
                    "Intensity": a["intensity"],
                    "Rel. %": round(a["rel"], 1),
                    "Annotated": "✓" if a["annotated"] else "✗",
                    "Ion formula": a["ion_formula"],
                    "NL": f"{a['nl']:.4f}" if a["nl"] is not None else "",
                    "ppm": f"{a['ppm']:.1f}" if a["ppm"] is not None else "",
                })
            st.dataframe(pd.DataFrame(rows), width="stretch")
    else:
        st.info("Set your SMILES and spectrum above, then click Run MetFrag.")
