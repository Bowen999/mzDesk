import math

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
import theme
import metfrag_engine as engine


DEFAULT_SMILES = r"CCCCCCCC\C=C/CCCCCCCC(=O)OCC(COP(O)(=O)OCC(O)CO)OC(=O)CCCCCCC\C=C/CCCCCCCC"

_KIND_LABEL = {"precursor": "Precursor ion", "fragment": "MetFrag fragment"}


def _chart_tokens():
    """Design-system colours handed to the chart iframe (it cannot see the global CSS)."""
    return {"ink": theme.INK, "muted": theme.MUTED, "line": theme.LINE, "blue": theme.BLUE,
            "blue_bright": theme.BLUE_BRIGHT, "dark": theme.PLOT_DARK,
            "light": theme.PLOT_LIGHT, "grid": theme.SURFACE_2}


def render():
    theme.tab_header(
        "MS² Interpretation",
        subtitle=(
            "Annotate an MS² spectrum with MetFrag in-silico fragment structures. "
            "Hover a peak to see its ion formula, mass error, neutral loss and fragment."
        ),
        kicker="MetFrag · fragment annotation",
    )

    # The setup lives in one bordered panel like the other tabs, so the four tabs
    # read as the same kind of page: a framed query on top, results below.
    with theme.panel("ms2-query"):
        theme.section("Query")

        col_in, col_fig = st.columns([1.35, 1], gap="large")

        with col_in:
            smiles = st.text_input(
                "SMILES",
                value=DEFAULT_SMILES,
                help="Neutral candidate structure — MetFrag adds the adduct itself.",
                key="ms2_smiles",
            )

            spectrum_text = st.text_area(
                "MS² peak list",
                value=engine.DEFAULT_SPECTRUM,
                height=219,  # 9 rows: 9*20 + 39 — see "Peak-list heights" in the README
                help="`mz intensity` per line (space, tab or comma separated), OR space-separated "
                     "`mz:intensity` pairs. Intensity optional; malformed lines and peaks with "
                     "intensity ≤ 0 are skipped.",
                key="ms2_spectrum",
            )

        # The preview sits before the adduct picker, so read the current adduct out of
        # the session to report a structure/adduct mismatch next to the picture.
        adduct_now = st.session_state.get("ms2_adduct") or next(iter(engine.ADDUCTS))
        info = engine.precursor_info(smiles, adduct_now) if smiles.strip() else None

        with col_fig:
            theme.subsection("Input structure")
            if smiles.strip():
                png = engine.render_input_png(smiles)
                if png:
                    st.image(png, width="stretch")
                    if info and info["error"]:
                        theme.note(info["error"])
                else:
                    st.warning("Invalid SMILES — cannot draw structure.")
            else:
                theme.note("Enter a SMILES to preview the structure.")

        c_adduct, c_prec, c_depth = st.columns([1, 1, 1.25], vertical_alignment="bottom")
        adduct = c_adduct.selectbox(
            "Adduct",
            list(engine.ADDUCTS.keys()),
            index=0,
            key="ms2_adduct",
        )
        precursor_input = c_prec.text_input(
            "Precursor m/z",
            value="",
            placeholder="auto",
            help="Measured precursor m/z (optional) — used for neutral losses and a mismatch "
                 "check. Blank = theoretical m/z of the structure + adduct.",
            key="ms2_precursor",
        )
        tree_depth = c_depth.slider(
            "Tree depth",
            1,
            4,
            value=2,
            help="Maximum number of consecutive bond cleavages MetFrag explores "
                 "(MetFrag default 2). Higher finds more, less specific fragments.",
            key="ms2_depth",
        )

        c_adv, c_run = st.columns([2.4, 1], vertical_alignment="bottom")
        with c_adv:
            with st.expander("Advanced tolerances"):
                t_abs, t_ppm = st.columns(2)
                tol_abs = t_abs.number_input(
                    "Absolute mass tolerance (Da)",
                    value=0.002,
                    min_value=0.0,
                    step=0.001,
                    format="%.4f",
                    help="MetFrag matches a fragment to a peak when |Δm| ≤ absolute + relative "
                         "tolerance (the two are added).",
                    key="ms2_tolabs",
                )
                tol_ppm = t_ppm.number_input(
                    "Relative tolerance (ppm)",
                    value=10.0,
                    min_value=0.0,
                    step=1.0,
                    help="Added to the absolute tolerance: window = abs + ppm × m/z × 10⁻⁶.",
                    key="ms2_tolppm",
                )
        run = c_run.button("Run MetFrag", type="primary", key="ms2_run", width="stretch")

        if info and not info["error"]:
            theme.note(f"Theoretical {adduct} m/z {info['theo_mz']:.4f} · ion {info['ion_formula']}")

    peaks = engine.parse_spectrum(spectrum_text)
    signature = (smiles.strip(), tuple(peaks), adduct, precursor_input.strip(),
                 int(tree_depth), float(tol_abs), float(tol_ppm))

    if run:
        # A new run replaces the previous result; on failure only the error is shown.
        st.session_state.pop("ms2_result", None)
        st.session_state.pop("ms2_result_sig", None)
        _run(smiles, peaks, adduct, precursor_input, tree_depth, tol_abs, tol_ppm, signature)

    result = st.session_state.get("ms2_result")
    if result is None:
        if not run:
            st.info("Set your SMILES and MS² peak list above, then click Run MetFrag.")
        return
    _render_result(result, stale=st.session_state.get("ms2_result_sig") != signature)


def _run(smiles, peaks, adduct, precursor_input, tree_depth, tol_abs, tol_ppm, signature):
    if not peaks:
        st.error("No valid peaks found in the MS² peak list. Check the format "
                 "(`mz intensity` per line, or `mz:intensity` pairs).")
        return

    precursor_val = None
    if precursor_input.strip():
        try:
            precursor_val = float(precursor_input)
        except ValueError:
            precursor_val = float("nan")
        if not math.isfinite(precursor_val) or precursor_val <= 0:
            st.error("Precursor m/z must be a positive number or left blank.")
            return

    try:
        with st.spinner("Running MetFrag…"):
            result = engine.run_metfrag(
                smiles,
                peaks,
                adduct,
                precursor_val,
                tree_depth,
                tol_abs,
                tol_ppm,
            )
    except engine.MetFragError as e:
        st.error(str(e))
        if e.log:
            with st.expander("MetFrag log (last lines)"):
                st.code(e.log, language=None)
        return
    except ValueError as e:
        st.error(str(e))
        return
    except Exception as e:  # unexpected failure: show it instead of a traceback
        st.error(f"Unexpected error while running MetFrag: {e}")
        return
    st.session_state["ms2_result"] = result
    st.session_state["ms2_result_sig"] = signature


def _render_result(result, stale):
    theme.divider()

    if stale:
        st.info("Inputs changed since this result was computed — click **Run MetFrag** to update it.")
    for msg in result.get("warnings", []):
        st.warning(msg)

    annotation = result["annotation"]
    n_ann = sum(1 for a in annotation if a["annotated"])
    theo = result["theo_precursor_mz"]

    theme.section("Summary")
    theme.metric_cards([
        ("Annotated peaks", f"{n_ann} / {len(annotation)}"),
        ("Precursor (entered)" if result["precursor_given"] else "Precursor (theor.)",
         f"{result['precursor_mz']:.4f}"),
        ("Neutral formula", result["mol_formula"]),
    ])
    parts = [f"{result['adduct']} ion {result['ion_formula']}, theoretical m/z {theo:.4f}"]
    if result["precursor_given"]:
        parts.append(f"entered precursor Δ {(result['precursor_mz'] - theo) / theo * 1e6:+.1f} ppm")
    parts.append(f"tree depth {result['tree_depth']}")
    parts.append(f"match tolerance {result['tol_abs']:g} Da + {result['tol_ppm']:g} ppm")
    theme.note(" · ".join(parts))

    theme.section("Spectrum")
    theme.chips([
        ("Annotated", theme.PLOT_DARK),
        ("Unannotated", theme.PLOT_LIGHT),
    ])
    theme.note(
        "Hover any peak for its ion formula, mass error, neutral loss and MetFrag fragment. "
        "Drag to zoom · double-click to reset · dotted line = theoretical precursor m/z."
    )
    components.html(
        engine.build_chart_html(result, _chart_tokens()),
        height=engine.CHART_HEIGHT + 8,
        scrolling=False,
    )

    with st.expander("Annotation table"):
        rows = []
        for a in annotation:
            is_frag = a["kind"] == "fragment"
            rows.append({
                "m/z": a["mz"],
                "Rel. int. (%)": a["rel"],
                "Intensity": a["intensity"],
                "Annotation": _KIND_LABEL.get(a["kind"], "—"),
                "Ion formula": a["ion_formula"],
                "Theor. m/z": a["calc_mz"],
                "Error (ppm)": a["ppm"],
                "Neutral loss": a["nl"] if is_frag else None,
                "Fragment SMILES": a["smiles"] if is_frag else "",
            })
        st.dataframe(
            pd.DataFrame(rows),
            width="stretch",
            hide_index=True,
            column_config={
                "m/z": st.column_config.NumberColumn(format="%.4f"),
                "Rel. int. (%)": st.column_config.NumberColumn(format="%.1f"),
                "Theor. m/z": st.column_config.NumberColumn(format="%.4f"),
                "Error (ppm)": st.column_config.NumberColumn(format="%.1f"),
                "Neutral loss": st.column_config.NumberColumn(format="%.4f"),
            },
        )
