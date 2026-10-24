import math

import pandas as pd
import streamlit as st
from mzdesk.shared import theme
from mzdesk.ms1 import msbuddy_engine as engine

DEMO_PRECURSOR = "195.0877"
DEFAULT_MS2 = """42.0338 5
56.0495 4
69.0447 6
83.0604 9
110.0713 18
138.0662 100
195.0877 32"""
DEFAULT_ISO = """195.0877 100
196.0902 10.3
197.0922 0.9"""

ION_MODES = ["Positive", "Negative"]
AUTO = "Auto"
ELEMENTS_NOTE = "Supported: C H N O P S F Cl Br I Na K"


def buddy_status(precursor_mz=None):
    lo, hi = engine.PRECURSOR_MZ_RANGE
    if precursor_mz is not None and not lo < precursor_mz < hi:
        return False, f"BUDDY supports precursor m/z {lo:.0f}–{hi:.0f} only."
    return True, ""


def _keep_choice(key, shadow_key, options):
    value = st.session_state.get(key)
    if value is None:
        previous = st.session_state.get(shadow_key)
        st.session_state[key] = previous if previous in options else options[0]
    else:
        st.session_state[shadow_key] = value


def _parse_precursor(text):
    text = (text or "").strip()
    if not text:
        return None, "Precursor m/z is required."
    try:
        value = float(text)
    except ValueError:
        return None, "Precursor m/z must be a number."
    if not math.isfinite(value) or value <= 0:
        return None, "Precursor m/z must be positive."
    return value, None


def _fmt(value, digits=4):
    return "—" if value is None else f"{value:.{digits}f}"


def _formula_mismatch(formula, precursor_mz, positive, adduct_choice, tol_ppm):
    adduct = engine.default_adduct(positive) if adduct_choice == AUTO else adduct_choice
    theoretical = engine.ion_mz(formula, adduct)
    if theoretical is None:
        return None
    error = engine.ppm_error(precursor_mz, theoretical)
    if abs(error) <= tol_ppm:
        return None
    others = engine.ADDUCTS_POSITIVE if positive else engine.ADDUCTS_NEGATIVE
    matches = []
    for other in others:
        other_mz = engine.ion_mz(formula, other) if other != adduct else None
        if other_mz is not None and abs(engine.ppm_error(precursor_mz, other_mz)) <= tol_ppm:
            matches.append(other)
    message = (
        f"{formula} as {adduct}: theoretical m/z {theoretical:.4f}; measured {precursor_mz:.4f} "
        f"(Δ {precursor_mz - theoretical:+.4f} Da, {error:+,.1f} ppm; tolerance {tol_ppm:g} ppm)."
    )
    if matches:
        return message + f" Use {matches[0]}."
    return message + " Check formula, adduct, and precursor."


def _input_checks(precursor_text, precursor_mz, precursor_error, positive, adduct_choice,
                  formula, ms2_text, ms2_peaks, iso_text, ms1_peaks, ms1_tol, halogen):
    checks = []
    if precursor_text.strip() and precursor_error:
        checks.append(("warning", precursor_error))
    if formula and not engine.formula_is_valid(formula):
        checks.append(("warning", f"Invalid formula `{formula}`. {ELEMENTS_NOTE}."))
    elif formula and precursor_mz is not None:
        mismatch = _formula_mismatch(formula, precursor_mz, positive, adduct_choice, ms1_tol)
        if mismatch:
            checks.append(("warning", mismatch))
    if ms2_text.strip() and not ms2_peaks:
        checks.append(("warning", "No valid MS² peaks; spectrum ignored."))
    if iso_text.strip() and not ms1_peaks:
        checks.append(("warning", "No valid MS¹ peaks; isotope pattern ignored."))
    if ms1_peaks and precursor_mz is not None:
        m0_mz, m0_int = min(ms1_peaks, key=lambda peak: abs(peak[0] - precursor_mz))
        error = engine.ppm_error(m0_mz, precursor_mz)
        if abs(error) > ms1_tol:
            checks.append((
                "warning",
                f"MS¹ has no peak within {ms1_tol:g} ppm of the precursor "
                f"(closest {m0_mz:.4f}, {error:+,.1f} ppm). BUDDY will ignore it.",
            ))
        elif not halogen and not formula:
            m2 = [intensity for mz, intensity in ms1_peaks if 1.98 <= mz - m0_mz <= 2.02]
            if m2 and max(m2) >= 0.2 * m0_int:
                checks.append((
                    "info",
                    f"M+2 intensity is {max(m2) / m0_int:.0%}; enable Include halogens for Cl/Br.",
                ))
    return checks


def _signature(precursor_text, positive, adduct_choice, formula, ms2_peaks, ms1_peaks,
               ms_instr, halogen):
    return (precursor_text.strip(), bool(positive), adduct_choice, formula,
            tuple(ms2_peaks), tuple(ms1_peaks), ms_instr, bool(halogen))


def render():
    with theme.panel("struct-query"):
        r1 = st.columns([1, 1.25, 1, 1.3])
        precursor_text = r1[0].text_input(
            "Precursor m/z *",
            value=DEMO_PRECURSOR,
            placeholder="e.g. 195.0877",
            help="Required; singly charged; BUDDY below 1500.",
            key="struct_precursor",
        )
        mode = r1[1].segmented_control(
            "Ion mode *", ION_MODES, default=ION_MODES[0], key="struct_mode",
            on_change=_keep_choice, args=("struct_mode", "_struct_mode_last", ION_MODES),
        )
        positive = mode != "Negative"
        adduct_choice = r1[2].selectbox(
            "Adduct",
            [AUTO] + (engine.ADDUCTS_POSITIVE if positive else engine.ADDUCTS_NEGATIVE),
            index=0,
            help="Auto uses [M+H]+ or [M-H]-.",
            key="struct_adduct",
        )
        formula_input = r1[3].text_input(
            "Formula (optional)",
            value="",
            placeholder="e.g. C8H10N4O2",
            help="Known neutral formula; skips ranking. " + ELEMENTS_NOTE + ".",
            key="struct_formula",
        )

        r2 = st.columns(2)
        ms2_text = r2[0].text_area(
            "MS² spectrum (optional)",
            value=DEFAULT_MS2,
            height=179,
            help="Lines: `m/z intensity`; or pairs: `m/z:intensity`.",
            key="struct_ms2",
        )
        iso_text = r2[1].text_area(
            "MS¹ isotope pattern (optional)",
            value=DEFAULT_ISO,
            height=179,
            help="M, M+1, M+2 …; first peak is the precursor; same format as MS².",
            key="struct_iso",
        )

        with st.expander("Advanced"):
            a1, a2 = st.columns(2)
            ms_instr = a1.selectbox(
                "Instrument", list(engine.INSTRUMENT_TOLS), index=0,
                help="MS¹/MS² ppm: orbitrap 5/10; qtof 10/20; fticr 2/5.",
                key="struct_instr",
            )
            halogen = a2.checkbox(
                "Include halogens", value=False,
                help="Adds F, Cl, Br, and I to BUDDY.",
                key="struct_halogen",
            )

        precursor_mz, precursor_error = _parse_precursor(precursor_text)
        formula_clean = formula_input.strip()
        ms2_peaks = tuple(engine.parse_peak_list(ms2_text))
        ms1_peaks = tuple(engine.parse_peak_list(iso_text))
        ms1_tol, ms2_tol = engine.INSTRUMENT_TOLS.get(ms_instr, (5, 10))
        for level, message in _input_checks(
                precursor_text, precursor_mz, precursor_error, positive, adduct_choice,
                formula_clean, ms2_text, ms2_peaks, iso_text, ms1_peaks, ms1_tol, halogen):
            (st.warning if level == "warning" else st.info)(message)

        buddy_ok, buddy_reason = buddy_status(precursor_mz)
        run = st.button("Annotate", type="primary", key="struct_run")
        if not buddy_ok:
            theme.note(buddy_reason)

    signature = _signature(precursor_text, positive, adduct_choice, formula_clean,
                           ms2_peaks, ms1_peaks, ms_instr, halogen)
    if run:
        _handle_run(precursor_mz, precursor_error, positive, adduct_choice, formula_clean,
                    ms2_peaks, ms1_peaks, ms_instr, ms2_tol, halogen, signature)

    result = st.session_state.get("struct_result")
    if result is not None and ("signature" not in result or "ran" in result):
        st.session_state.pop("struct_result", None)
        result = None
    if result is None:
        st.info("Caffeine demo loaded. Click Annotate.")
    else:
        _render_results(result, signature)


def _handle_run(precursor_mz, precursor_error, positive, adduct_choice, formula,
                ms2_peaks, ms1_peaks, ms_instr, ms2_tol, halogen, signature):
    if precursor_error:
        st.error(precursor_error)
        return
    if formula and not engine.formula_is_valid(formula):
        st.error(f"Invalid formula `{formula}`. {ELEMENTS_NOTE}.")
        return
    buddy_ok, buddy_reason = buddy_status(precursor_mz)
    if not buddy_ok:
        st.error(buddy_reason)
        return

    adduct = None if adduct_choice == AUTO else adduct_choice
    combined = {
        "precursor_mz": precursor_mz,
        "ms2_peaks": ms2_peaks,
        "ms2_tol_ppm": ms2_tol,
        "signature": signature,
        "buddy": None,
        "buddy_error": None,
    }

    try:
        with st.spinner("Running BUDDY…"):
            combined["buddy"] = engine.annotate_formula(
                precursor_mz=precursor_mz,
                positive=positive,
                adduct=adduct,
                rt=None,
                formula_input=formula,
                ms2_peaks=ms2_peaks,
                ms1_peaks=ms1_peaks,
                ms_instr=ms_instr,
                halogen=halogen,
            )
    except Exception as exc:
        combined["buddy_error"] = f"BUDDY failed: {exc}"

    st.session_state["struct_result"] = combined
    st.session_state.pop("struct_selection", None)
    st.session_state.pop("_struct_selection_last", None)


def _warn(message):
    st.warning(message.replace("\n", "  \n"))


def _render_results(result, signature):
    theme.divider()
    if result["signature"] != signature:
        st.warning("Inputs changed. Annotate again to update.")

    chosen = None
    theme.section("BUDDY")
    buddy = result.get("buddy")
    if result.get("buddy_error"):
        _warn(result["buddy_error"])
    elif not buddy["candidates"]:
        st.warning("No BUDDY candidates. Check precursor/adduct or Advanced settings.")
    else:
        chosen = _render_buddy(result, buddy)

    peaks = result["ms2_peaks"]
    if chosen and peaks:
        formula, adduct, source = chosen
        _render_fragments(peaks, formula, adduct, source,
                          result["precursor_mz"], result["ms2_tol_ppm"])
    elif chosen:
        theme.note("No MS² spectrum; formula annotation only.")


def _render_buddy(result, buddy):
    candidates = buddy["candidates"]
    adduct = buddy["adduct_used"]
    if buddy["from_input"]:
        labels = [f"Input · {candidates[0]['formula']}"]
    else:
        labels = [f"#{candidate['rank']} {candidate['formula']}" for candidate in candidates]

    index = 0
    if len(candidates) > 1:
        selection = st.segmented_control(
            "Candidate", labels, default=labels[0], key="struct_selection",
            on_change=_keep_choice,
            args=("struct_selection", "_struct_selection_last", labels),
        )
        index = labels.index(selection) if selection in labels else 0
    chosen = candidates[index]

    fdr = chosen["fdr"]
    if buddy["from_input"]:
        subtitle = "Input formula; ranking skipped"
    else:
        subtitle = (f"Rank #{chosen['rank']} of {len(candidates)}"
                    + (f" · estimated FDR {fdr:.1%}" if fdr is not None else ""))
    theme.formula_verdict(chosen["formula"], subtitle, theme.BLUE)

    stats = [
        ("Theoretical m/z", _fmt(chosen["mz_theo"])),
        ("Mass error", "—" if chosen["ppm"] is None else f"{chosen['ppm']:+.2f} ppm"),
        ("Adduct", adduct),
    ]
    if not buddy["from_input"]:
        stats += [
            ("ML probability", _fmt(chosen["prob"], 3)),
            ("Normalised prob.", _fmt(chosen["normed_prob"], 3)),
            ("Est. FDR", "—" if fdr is None else f"{fdr:.1%}"),
        ]
    theme.stat_grid(stats, wide=True, max_cols=3)
    if buddy.get("ms1_used") == 0 and result["signature"][5]:
        theme.note("MS¹ ignored: no peak matched the precursor within tolerance.")

    if len(candidates) > 1:
        with st.expander(f"Candidates (top {len(candidates)})"):
            st.dataframe(
                pd.DataFrame([{
                    "Rank": candidate["rank"],
                    "Formula": candidate["formula"],
                    "Theor. m/z": candidate["mz_theo"],
                    "Δ ppm": candidate["ppm"],
                    "ML probability": candidate["prob"],
                    "Norm. prob.": candidate["normed_prob"],
                    "Est. FDR": candidate["fdr"],
                    "Iso. similarity": candidate["isotope_sim"],
                } for candidate in candidates]),
                width="stretch",
                hide_index=True,
                column_config={
                    "Theor. m/z": st.column_config.NumberColumn(format="%.4f"),
                    "Δ ppm": st.column_config.NumberColumn(format="%+.2f"),
                    "ML probability": st.column_config.NumberColumn(format="%.4f"),
                    "Norm. prob.": st.column_config.NumberColumn(format="%.4f"),
                    "Est. FDR": st.column_config.NumberColumn(format="%.4f"),
                    "Iso. similarity": st.column_config.NumberColumn(format="%.4f"),
                },
            )

    source = "user formula" if buddy["from_input"] else f"BUDDY rank #{chosen['rank']}"
    return chosen["formula"], adduct, source


def _render_fragments(peaks, formula, adduct, source, precursor_mz, ms2_tol_ppm):
    theme.section("Fragment formulas")
    try:
        explained = engine.explain_fragments(
            tuple(peak[0] for peak in peaks), formula, adduct, ms2_tol_ppm)
    except Exception as exc:
        st.warning(f"Could not assign fragment formulas: {exc}")
        return
    rows = engine.build_annotation_rows(peaks, explained, precursor_mz)

    n_annotated = sum(1 for row in rows if row["annotated"])
    total = sum(row["intensity"] for row in rows)
    explained_intensity = sum(row["intensity"] for row in rows if row["annotated"])
    theme.metric_cards([
        ("Peaks", str(len(rows))),
        ("Explained", f"{n_annotated}/{len(rows)}"),
        ("Explained intensity", f"{explained_intensity / total:.0%}" if total else "—"),
        ("Adduct", adduct),
        ("MS² tolerance", f"{ms2_tol_ppm:g} ppm"),
    ])
    theme.chips([("Explained", theme.PLOT_DARK), ("Unexplained", theme.PLOT_LIGHT)])
    theme.note(f"{formula} {adduct} · {source}. Hover for details; drag to zoom.")

    figure = engine.build_figure(
        rows, precursor_mz, adduct=adduct,
        dark=theme.PLOT_DARK, light=theme.PLOT_LIGHT, ink=theme.INK,
        muted=theme.MUTED, grid=theme.SURFACE_2, highlight=theme.BLUE_BRIGHT,
    )
    st.plotly_chart(figure, width="stretch", theme=None, config={"displayModeBar": False})

    with st.expander("Fragment table"):
        st.dataframe(
            pd.DataFrame([{
                "m/z": row["mz"],
                "Rel. int. (%)": row["rel"],
                "Intensity": row["intensity"],
                "Explained": "✓" if row["annotated"] else "✗",
                "Ion formula": row["formula"] + (" (precursor)" if row["precursor"] else ""),
                "Δ ppm": row["ppm"],
                "Neutral loss": row["nl"],
                "Other fits": ", ".join(row["alternatives"]),
            } for row in rows]),
            width="stretch",
            hide_index=True,
            column_config={
                "m/z": st.column_config.NumberColumn(format="%.4f"),
                "Rel. int. (%)": st.column_config.NumberColumn(format="%.1f"),
                "Δ ppm": st.column_config.NumberColumn(format="%+.1f"),
                "Neutral loss": st.column_config.NumberColumn(format="%.4f"),
            },
        )
