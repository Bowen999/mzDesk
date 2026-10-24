import math

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from mzdesk.shared import theme
from mzdesk.ms2 import metfrag_engine as engine

DEFAULT_SMILES = r"CCCCCCCC\C=C/CCCCCCCC(=O)OCC(COP(O)(=O)OCC(O)CO)OC(=O)CCCCCCC\C=C/CCCCCCCC"
_KIND_LABEL = {"precursor": "Precursor ion", "fragment": "MetFrag fragment"}


def _chart_tokens():
    return {
        "ink": theme.INK,
        "muted": theme.MUTED,
        "line": theme.LINE,
        "blue": theme.BLUE,
        "blue_bright": theme.BLUE_BRIGHT,
        "dark": theme.PLOT_DARK,
        "light": theme.PLOT_LIGHT,
        "grid": theme.SURFACE_2,
    }


def render():
    with theme.panel("ms2-query"):
        col_in, col_fig = st.columns([1.35, 1], gap="large")
        with col_in:
            smiles = st.text_input(
                "SMILES",
                value=DEFAULT_SMILES,
                help="Neutral structure; MetFrag adds the adduct.",
                key="ms2_smiles",
            )
            spectrum_text = st.text_area(
                "MS² peak list",
                value=engine.DEFAULT_SPECTRUM,
                height=219,
                help="Lines: `m/z intensity`; or pairs: `m/z:intensity`.",
                key="ms2_spectrum",
            )

        adduct_now = st.session_state.get("ms2_adduct") or next(iter(engine.ADDUCTS))
        info = engine.precursor_info(smiles, adduct_now) if smiles.strip() else None

        with col_fig:
            theme.subsection("Structure")
            if smiles.strip():
                png = engine.render_input_png(smiles)
                if png:
                    st.image(png, width="stretch")
                    if info and info["error"]:
                        theme.note(info["error"])
                else:
                    st.warning("Invalid SMILES.")

        c_adduct, c_prec, c_depth = st.columns([1, 1, 1.25], vertical_alignment="bottom")
        adduct = c_adduct.selectbox("Adduct", list(engine.ADDUCTS), index=0, key="ms2_adduct")
        precursor_input = c_prec.text_input(
            "Precursor m/z",
            value="",
            placeholder="Auto",
            help="Optional; blank uses theoretical m/z.",
            key="ms2_precursor",
        )
        tree_depth = c_depth.slider(
            "Tree depth", 1, 4, value=2,
            help="Maximum consecutive bond cleavages.",
            key="ms2_depth",
        )

        c_adv, c_run = st.columns([2.4, 1], vertical_alignment="bottom")
        with c_adv:
            with st.expander("Tolerances"):
                t_abs, t_ppm = st.columns(2)
                tol_abs = t_abs.number_input(
                    "Absolute (Da)", value=0.002, min_value=0.0, step=0.001, format="%.4f",
                    help="Match window = absolute + relative tolerance.",
                    key="ms2_tolabs",
                )
                tol_ppm = t_ppm.number_input(
                    "Relative (ppm)", value=10.0, min_value=0.0, step=1.0,
                    help="Relative window = ppm × m/z × 10⁻⁶.",
                    key="ms2_tolppm",
                )
        run = c_run.button("Run MetFrag", type="primary", key="ms2_run", width="stretch")

        if info and not info["error"]:
            theme.note(f"Theoretical {adduct}: m/z {info['theo_mz']:.4f} · {info['ion_formula']}")

    peaks = engine.parse_spectrum(spectrum_text)
    signature = (smiles.strip(), tuple(peaks), adduct, precursor_input.strip(),
                 int(tree_depth), float(tol_abs), float(tol_ppm))

    if run:
        st.session_state.pop("ms2_result", None)
        st.session_state.pop("ms2_result_sig", None)
        _run(smiles, peaks, adduct, precursor_input, tree_depth, tol_abs, tol_ppm, signature)

    result = st.session_state.get("ms2_result")
    if result is None:
        if not run:
            st.info("Click Run MetFrag.")
        return
    _render_result(result, stale=st.session_state.get("ms2_result_sig") != signature)


def _run(smiles, peaks, adduct, precursor_input, tree_depth, tol_abs, tol_ppm, signature):
    if not peaks:
        st.error("No valid MS² peaks. Use `m/z intensity` lines or `m/z:intensity` pairs.")
        return

    precursor_val = None
    if precursor_input.strip():
        try:
            precursor_val = float(precursor_input)
        except ValueError:
            precursor_val = float("nan")
        if not math.isfinite(precursor_val) or precursor_val <= 0:
            st.error("Precursor m/z must be positive or blank.")
            return

    try:
        with st.spinner("Running MetFrag…"):
            result = engine.run_metfrag(
                smiles, peaks, adduct, precursor_val, tree_depth, tol_abs, tol_ppm,
            )
    except engine.MetFragError as exc:
        st.error(str(exc))
        if exc.log:
            with st.expander("MetFrag log"):
                st.code(exc.log, language=None)
        return
    except ValueError as exc:
        st.error(str(exc))
        return
    except Exception as exc:
        st.error(f"MetFrag failed: {exc}")
        return
    st.session_state["ms2_result"] = result
    st.session_state["ms2_result_sig"] = signature


def _render_result(result, stale):
    theme.divider()
    if stale:
        st.info("Inputs changed. Run MetFrag to update.")
    for message in result.get("warnings", []):
        st.warning(message)

    annotation = result["annotation"]
    n_annotated = sum(1 for item in annotation if item["annotated"])
    theoretical = result["theo_precursor_mz"]

    theme.section("Summary")
    theme.metric_cards([
        ("Annotated peaks", f"{n_annotated} / {len(annotation)}"),
        ("Precursor (entered)" if result["precursor_given"] else "Precursor (theoretical)",
         f"{result['precursor_mz']:.4f}"),
        ("Neutral formula", result["mol_formula"]),
    ])
    parts = [f"{result['adduct']} · {result['ion_formula']} · theoretical m/z {theoretical:.4f}"]
    if result["precursor_given"]:
        parts.append(f"precursor Δ {(result['precursor_mz'] - theoretical) / theoretical * 1e6:+.1f} ppm")
    parts.append(f"depth {result['tree_depth']}")
    parts.append(f"tolerance {result['tol_abs']:g} Da + {result['tol_ppm']:g} ppm")
    theme.note(" · ".join(parts))

    theme.section("Spectrum")
    theme.chips([("Annotated", theme.PLOT_DARK), ("Unannotated", theme.PLOT_LIGHT)])
    theme.note("Hover for details; drag to zoom. Dotted line: theoretical precursor.")
    components.html(engine.build_chart_html(result, _chart_tokens()),
                    height=engine.CHART_HEIGHT + 8, scrolling=False)

    with st.expander("Annotations"):
        rows = []
        for item in annotation:
            is_fragment = item["kind"] == "fragment"
            rows.append({
                "m/z": item["mz"],
                "Rel. int. (%)": item["rel"],
                "Intensity": item["intensity"],
                "Annotation": _KIND_LABEL.get(item["kind"], "—"),
                "Ion formula": item["ion_formula"],
                "Theor. m/z": item["calc_mz"],
                "Error (ppm)": item["ppm"],
                "Neutral loss": item["nl"] if is_fragment else None,
                "Fragment SMILES": item["smiles"] if is_fragment else "",
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
