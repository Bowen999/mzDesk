import math

import pandas as pd
import streamlit as st
import theme
import msbuddy_engine as engine
import sirius_engine


# Demo feature, runnable out of the box: caffeine (C8H10N4O2) as [M+H]+, a well-known
# standard with a known answer. MS² = the classic caffeine CID fragment ions plus the
# precursor; MS¹ = the theoretical isotope envelope of C8H11N4O2+.
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
ELEMENTS_NOTE = "supported elements: C H N O P S F Cl Br I Na K"


def decide_engines(has_ms1, has_ms2, precursor_mz=None):
    """Which annotation engines can run given the current inputs.

    BUDDY: needs only precursor m/z (within msbuddy's range) + ion mode.
    SIRIUS: additionally requires at least one spectrum (MS¹ isotopes or MS²)
    and must be installed.
    """
    decision = {
        "buddy": {"ok": True, "reason": "runs on m/z + ion mode"},
        "sirius": {"ok": False, "reason": ""},
    }
    lo, hi = engine.PRECURSOR_MZ_RANGE
    if precursor_mz is not None and not lo < precursor_mz < hi:
        decision["buddy"] = {"ok": False,
                             "reason": f"supports precursor m/z {lo:.0f}–{hi:.0f} only"}
    if not (has_ms1 or has_ms2):
        decision["sirius"]["reason"] = "needs an MS¹ isotope pattern and/or MS² spectrum"
    elif sirius_engine.sirius_available() is None:
        decision["sirius"]["reason"] = "sirius executable not found"
    else:
        decision["sirius"] = {"ok": True, "reason": "spectra available"}
    return decision


def _keep_choice(key, shadow_key, options):
    """on_change for single-select segmented controls: clicking the active option
    re-selects it instead of clearing the control (Streamlit allows deselection)."""
    value = st.session_state.get(key)
    if value is None:
        previous = st.session_state.get(shadow_key)
        st.session_state[key] = previous if previous in options else options[0]
    else:
        st.session_state[shadow_key] = value


def _parse_precursor(text):
    """(precursor m/z or None, error message or None)."""
    text = (text or "").strip()
    if not text:
        return None, "Precursor m/z is required."
    try:
        value = float(text)
    except ValueError:
        return None, "Precursor m/z must be a number."
    if not math.isfinite(value) or value <= 0:
        return None, "Precursor m/z must be a positive number."
    return value, None


def _fmt(value, digits=4):
    return "—" if value is None else f"{value:.{digits}f}"


def _formula_mismatch(formula, precursor_mz, positive, adduct_choice, tol_ppm):
    """Warning text when the formula + adduct does not reproduce the precursor m/z."""
    adduct = engine.default_adduct(positive) if adduct_choice == AUTO else adduct_choice
    theo = engine.ion_mz(formula, adduct)
    if theo is None:
        return None
    err = engine.ppm_error(precursor_mz, theo)
    if abs(err) <= tol_ppm:
        return None
    others = engine.ADDUCTS_POSITIVE if positive else engine.ADDUCTS_NEGATIVE
    matches = []
    for other in others:
        other_mz = engine.ion_mz(formula, other) if other != adduct else None
        if other_mz is not None and abs(engine.ppm_error(precursor_mz, other_mz)) <= tol_ppm:
            matches.append(other)
    msg = (f"Formula `{formula}` as {adduct} has m/z {theo:.4f}, but the precursor is "
           f"{precursor_mz:.4f} (Δ {precursor_mz - theo:+.4f} Da / {err:+,.1f} ppm; "
           f"MS¹ tolerance {tol_ppm} ppm).")
    if matches:
        return msg + f" It matches {matches[0]} — select that adduct."
    return msg + (" Check the formula, adduct and precursor — otherwise the MS² peaks are "
                  "explained with a formula that does not fit the precursor.")


def _input_checks(precursor_text, precursor_mz, precursor_error, positive, adduct_choice,
                  formula, ms2_text, ms2_peaks, iso_text, ms1_peaks, ms1_tol, halogen):
    """Live validation of the query inputs -> [(level, message)]."""
    checks = []
    if precursor_text.strip() and precursor_error:
        checks.append(("warning", precursor_error))
    if formula and not engine.formula_is_valid(formula):
        checks.append(("warning", f"`{formula}` is not a valid formula ({ELEMENTS_NOTE})."))
    elif formula and precursor_mz is not None:
        mismatch = _formula_mismatch(formula, precursor_mz, positive, adduct_choice, ms1_tol)
        if mismatch:
            checks.append(("warning", mismatch))
    if ms2_text.strip() and not ms2_peaks:
        checks.append(("warning", "No valid peaks could be read from the MS² spectrum, so it "
                                  "will be ignored — use one `m/z intensity` pair per line."))
    if iso_text.strip() and not ms1_peaks:
        checks.append(("warning", "No valid peaks could be read from the MS¹ isotope pattern, "
                                  "so it will be ignored."))
    if ms1_peaks and precursor_mz is not None:
        m0_mz, m0_int = min(ms1_peaks, key=lambda p: abs(p[0] - precursor_mz))
        err = engine.ppm_error(m0_mz, precursor_mz)
        if abs(err) > ms1_tol:
            checks.append(("warning",
                           f"The MS¹ isotope pattern has no peak within {ms1_tol} ppm of the "
                           f"precursor m/z (closest: {m0_mz:.4f}, {err:+,.1f} ppm). BUDDY needs "
                           "the monoisotopic peak there to use the pattern."))
        elif not halogen and not formula:  # the hint only matters for BUDDY ranking
            m2 = [inten for mz, inten in ms1_peaks if 1.98 <= mz - m0_mz <= 2.02]
            if m2 and max(m2) >= 0.2 * m0_int:
                checks.append(("info",
                               f"The M+2 isotope peak is {max(m2) / m0_int:.0%} of the "
                               "monoisotopic peak — typical of Cl or Br. Enable **Include "
                               "halogens** in Advanced so BUDDY considers them."))
    return checks


def _signature(precursor_text, positive, adduct_choice, formula, ms2_peaks, ms1_peaks,
               ms_instr, halogen):
    """Everything that changes the annotation result (RT and timeout do not)."""
    return (precursor_text.strip(), bool(positive), adduct_choice, formula,
            tuple(ms2_peaks), tuple(ms1_peaks), ms_instr, bool(halogen))


def render():
    theme.tab_header(
        "Structure Annotation",
        subtitle=(
            "Annotate the molecular formula of an unknown feature with BUDDY and SIRIUS — "
            "engines are picked automatically from your inputs."
        ),
        kicker="msbuddy · sirius · molecular formula",
    )

    with theme.panel("struct-query"):
        theme.section("Query feature")
        r1 = st.columns([1, 1.25, 1, 1.3])
        precursor_text = r1[0].text_input(
            "Precursor m/z *",
            value=DEMO_PRECURSOR,
            placeholder="e.g. 195.0877",
            help="Required. Singly charged precursor ion m/z (BUDDY: below 1500).",
            key="struct_precursor",
        )
        mode = r1[1].segmented_control(
            "Ion mode *",
            ION_MODES,
            default=ION_MODES[0],
            key="struct_mode",
            on_change=_keep_choice,
            args=("struct_mode", "_struct_mode_last", ION_MODES),
        )
        positive = (mode != "Negative")
        adduct_choice = r1[2].selectbox(
            "Adduct",
            [AUTO] + (engine.ADDUCTS_POSITIVE if positive else engine.ADDUCTS_NEGATIVE),
            index=0,
            help="Auto: BUDDY assumes [M+H]+ / [M-H]-, SIRIUS detects the adduct itself "
                 "([M+?] wildcard).",
            key="struct_adduct",
        )
        formula_input = r1[3].text_input(
            "Formula (optional)",
            value="",
            placeholder="e.g. C8H10N4O2",
            help="Known neutral formula (" + ELEMENTS_NOTE + "). Skips BUDDY ranking: the "
                 "MS² peaks are explained with this formula and SIRIUS is restricted to it.",
            key="struct_formula",
        )

        r2 = st.columns(2)
        ms2_text = r2[0].text_area(
            "MS² spectrum (optional)",
            value=DEFAULT_MS2,
            height=179,  # 7 rows: 7*20 + 39 — see "Peak-list heights" in the README
            help="`m/z intensity` per line, OR space-separated `m/z:intensity` pairs "
                 "(e.g. `138.0662:100 110.0713:18`). Used by BUDDY's bottom-up search, "
                 "SIRIUS fragmentation trees and the fragment-formula plot.",
            key="struct_ms2",
        )
        iso_text = r2[1].text_area(
            "MS¹ isotope pattern (optional)",
            value=DEFAULT_ISO,
            height=179,  # 7 rows: 7*20 + 39 — see "Peak-list heights" in the README
            help="Isotope peaks of the precursor (M, M+1, M+2 …) for isotope-similarity "
                 "scoring; the first peak should be the precursor m/z. Same formats as MS². "
                 "Either MS¹ or MS² unlocks SIRIUS.",
            key="struct_iso",
        )

        with st.container():
            with st.expander("Advanced"):
                a1, a2, a3 = st.columns(3)
                ms_instr = a1.selectbox(
                    "Instrument", list(engine.INSTRUMENT_TOLS.keys()), index=0,
                    help="Mass tolerances (MS¹ / MS²): orbitrap 5 / 10 ppm · qtof 10 / 20 ppm "
                         "· fticr 2 / 5 ppm. SIRIUS has no FT-ICR profile and uses its "
                         "orbitrap profile for fticr.",
                    key="struct_instr",
                )
                halogen = a2.checkbox(
                    "Include halogens (F Cl Br I)", value=False,
                    help="Adds F, Cl, Br and I to BUDDY's formula space. SIRIUS detects Cl/Br "
                         "from the MS¹ isotope pattern on its own.",
                    key="struct_halogen")
                timeout_secs = a3.number_input(
                    "SIRIUS timeout (s)", min_value=60, value=600, step=60,
                    help="Wall-clock limit for the SIRIUS run. BUDDY single queries take "
                         "seconds and have no timeout (msbuddy only enforces it in batch mode).",
                    key="struct_timeout")

        # ---- derived inputs + live validation ----
        precursor_mz, precursor_error = _parse_precursor(precursor_text)
        formula_clean = formula_input.strip()
        ms2_peaks = tuple(engine.parse_peak_list(ms2_text))
        ms1_peaks = tuple(engine.parse_peak_list(iso_text))
        ms1_tol, ms2_tol = engine.INSTRUMENT_TOLS.get(ms_instr, (5, 10))
        for level, msg in _input_checks(precursor_text, precursor_mz, precursor_error,
                                        positive, adduct_choice, formula_clean,
                                        ms2_text, ms2_peaks, iso_text, ms1_peaks,
                                        ms1_tol, halogen):
            (st.warning if level == "warning" else st.info)(msg)

        # ---- live engine preview + run ----
        # No rule above this row: the subsection label drew a hairline, and together
        # with the divider that was ~80 px of air between the expander and two chips.
        decision = decide_engines(bool(ms1_peaks), bool(ms2_peaks), precursor_mz)
        c_eng, c_run = st.columns([2.4, 1], vertical_alignment="center")
        with c_eng:
            theme.chip_line("Engines to run", [
                (("✓ " if d["ok"] else "✗ ") + name,
                 theme.BLUE if d["ok"] else theme.CATEGORY_ACCENT["Invalid"])
                for name, d in (("BUDDY", decision["buddy"]), ("SIRIUS", decision["sirius"]))
            ])
        run = c_run.button("Annotate formula", type="primary", key="struct_run",
                           width="stretch")
        skipped = [k.upper() + ": " + v["reason"] for k, v in decision.items() if not v["ok"]]
        if skipped:
            theme.note("Skipped — " + "; ".join(skipped) + ".")

    signature = _signature(precursor_text, positive, adduct_choice, formula_clean,
                           ms2_peaks, ms1_peaks, ms_instr, halogen)
    if run:
        _handle_run(precursor_mz, precursor_error, positive, adduct_choice, formula_clean,
                    ms2_peaks, ms1_peaks, decision, ms_instr, ms2_tol, halogen,
                    timeout_secs, signature)

    result = st.session_state.get("struct_result")
    if result is not None and "signature" not in result:
        result = None  # stored by an older version of this tab
    if result is None:
        st.info("A demo feature (caffeine, [M+H]+) is pre-filled — click **Annotate formula**, "
                "or enter your own precursor m/z and spectra above.")
    else:
        _render_results(result, signature)

    theme.divider()
    theme.section("How it works")
    theme.notes([
        "BUDDY enumerates candidate formulas for the precursor m/z, shrinks the space "
        "by explaining MS² peaks bottom-up, ranks candidates with machine learning and "
        "reports an estimated FDR.",
        "SIRIUS 6 scores the MS¹ isotope pattern and MS² fragmentation trees; it only runs "
        "when an MS¹ isotope pattern and/or MS² spectrum is provided, sirius is installed "
        "and you are logged in (`sirius login`).",
        "Provide a Formula to skip ranking — it is fed to both engines directly. "
        "Otherwise the BUDDY Top-1 is shown and can be switched with one click.",
        "MS² fragment-ion formulas are subformulas of the chosen formula + adduct within the "
        "instrument's MS² tolerance (msbuddy); when several fit, the smallest mass error is "
        "shown and the others are listed. Mass errors are (observed − theoretical) / "
        "theoretical.",
    ])
    theme.note("Xing et al. BUDDY. Nat Methods 20, 881-890 (2023) · "
               "Dührkop et al. SIRIUS 4. Nat Methods 16, 299-302 (2019).")


def _handle_run(precursor_mz, precursor_error, positive, adduct_choice, formula,
                ms2_peaks, ms1_peaks, decision, ms_instr, ms2_tol, halogen,
                timeout_secs, signature):
    if precursor_error:
        st.error(precursor_error)
        return
    if formula and not engine.formula_is_valid(formula):
        st.error(f"`{formula}` is not a valid formula ({ELEMENTS_NOTE}).")
        return
    if not any(d["ok"] for d in decision.values()):
        reasons = "; ".join(f"{k.upper()} {v['reason']}" for k, v in decision.items())
        st.error(f"No engine can run with the current inputs ({reasons}).")
        return

    adduct = None if adduct_choice == AUTO else adduct_choice
    combined = {
        "precursor_mz": precursor_mz,
        "positive": positive,
        "ms2_peaks": ms2_peaks,
        "ms2_tol_ppm": ms2_tol,
        "adduct_explain": adduct or engine.default_adduct(positive),
        "signature": signature,
        "ran": {"buddy": decision["buddy"]["ok"], "sirius": decision["sirius"]["ok"]},
        "buddy": None, "buddy_error": None,
        "sirius": None, "sirius_error": None,
    }

    if decision["buddy"]["ok"]:
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
        except Exception as e:
            combined["buddy_error"] = f"BUDDY failed: {e}"

    if decision["sirius"]["ok"]:
        try:
            with st.spinner("Running SIRIUS… (this can take a few minutes)"):
                combined["sirius"] = sirius_engine.run_sirius(
                    precursor_mz=precursor_mz,
                    positive=positive,
                    adduct=(adduct or sirius_engine.auto_adduct(positive)),
                    formula_input=formula,
                    ms2_peaks=ms2_peaks,
                    ms1_peaks=ms1_peaks,
                    ms_instr=ms_instr,
                    timeout_secs=int(timeout_secs),
                )
        except sirius_engine.SiriusLoginError as e:
            combined["sirius_error"] = str(e)
        except Exception as e:
            combined["sirius_error"] = f"SIRIUS failed: {e}"

    st.session_state["struct_result"] = combined
    st.session_state.pop("struct_selection", None)
    st.session_state.pop("_struct_selection_last", None)


def _warn(message):
    st.warning(message.replace("\n", "  \n"))  # keep multi-line engine logs readable


def _render_results(result, signature):
    theme.divider()
    if result["signature"] != signature:
        st.warning("The inputs above changed after these results were computed — click "
                   "**Annotate formula** to update them.")

    chosen = None  # (formula, adduct, source label) used for the fragment plot
    if result["ran"]["buddy"]:
        theme.section("BUDDY")
        buddy = result.get("buddy")
        if result.get("buddy_error"):
            _warn(result["buddy_error"])
        elif not buddy["candidates"]:
            st.warning("BUDDY returned no candidate formulas — check the precursor m/z and "
                       "adduct, or try Include halogens / another instrument in Advanced.")
        else:
            chosen = _render_buddy(result, buddy)

    if result["ran"]["sirius"]:
        theme.section("SIRIUS")
        if result.get("sirius_error"):
            _warn(result["sirius_error"])
        else:
            sirius_choice = _render_sirius(result, result["sirius"])
            if chosen is None:
                chosen = sirius_choice  # fall back to the SIRIUS top hit for the fragments

    peaks = result["ms2_peaks"]
    if chosen and peaks:
        formula, adduct, source = chosen
        _render_fragments(peaks, formula, adduct, source,
                          result["precursor_mz"], result["ms2_tol_ppm"])
    elif chosen:
        theme.note("No MS² spectrum was provided — showing the formula annotation only.")


def _render_buddy(result, buddy):
    """Render the BUDDY block. Returns (formula, adduct, source) for the fragment plot."""
    candidates = buddy["candidates"]
    adduct = buddy["adduct_used"]
    if buddy["from_input"]:
        labels = [f"Input · {candidates[0]['formula']}"]
    else:
        labels = [f"#{c['rank']} {c['formula']}" for c in candidates]

    idx = 0
    if len(candidates) > 1:
        selection = st.segmented_control(
            "Switch candidate", labels, default=labels[0], key="struct_selection",
            on_change=_keep_choice, args=("struct_selection", "_struct_selection_last", labels))
        idx = labels.index(selection) if selection in labels else 0
    chosen = candidates[idx]

    fdr = chosen["fdr"]
    if buddy["from_input"]:
        sub = "user-provided formula · BUDDY ranking skipped"
    else:
        sub = (f"rank #{chosen['rank']} of {len(candidates)}"
               + (f" · estimated FDR {fdr:.1%}" if fdr is not None else ""))
    theme.formula_verdict(chosen["formula"], sub, theme.BLUE)

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
    if not buddy["from_input"]:
        theme.note("ML probability is BUDDY's raw score for this formula; the normalised "
                   "probability spreads it over the candidate list and is what the estimated "
                   "FDR is derived from.")
    if buddy.get("ms1_used") == 0 and result["signature"][5]:
        theme.note("The MS¹ isotope pattern was not used: none of its peaks matched the "
                   "precursor m/z within the MS¹ tolerance.")

    if len(candidates) > 1:
        with st.expander(f"All candidates (top {len(candidates)})"):
            st.dataframe(
                pd.DataFrame([{
                    "Rank": c["rank"],
                    "Formula": c["formula"],
                    "Theor. m/z": c["mz_theo"],
                    "Δ ppm": c["ppm"],
                    "ML probability": c["prob"],
                    "Norm. prob.": c["normed_prob"],
                    "Est. FDR": c["fdr"],
                    "Iso. similarity": c["isotope_sim"],
                } for c in candidates]),
                width="stretch", hide_index=True,
                column_config={
                    "Theor. m/z": st.column_config.NumberColumn(format="%.4f"),
                    "Δ ppm": st.column_config.NumberColumn(format="%+.2f"),
                    "ML probability": st.column_config.NumberColumn(format="%.4f"),
                    "Norm. prob.": st.column_config.NumberColumn(format="%.4f"),
                    "Est. FDR": st.column_config.NumberColumn(format="%.4f"),
                    "Iso. similarity": st.column_config.NumberColumn(format="%.4f"),
                })

    source = "user formula" if buddy["from_input"] else f"BUDDY rank #{chosen['rank']}"
    return chosen["formula"], adduct, source


def _render_sirius(result, siri):
    """Render the SIRIUS block. Returns (formula, adduct, source) of the top hit."""
    cands = siri["candidates"]
    top = cands[0]
    sub = "SIRIUS top hit"
    if top.get("score") is not None:
        sub += f" · SIRIUS score {top['score']:.2f}"
    if top.get("score_norm") is not None:
        sub += f" · normalized {top['score_norm']:.3f}"
    theme.formula_verdict(top["formula"], sub, theme.CATEGORY_ACCENT["Metabolome"])

    explained_int = top.get("explained_int")
    mass_error = top.get("mass_error_ppm")
    peaks_expl = top.get("explained_peaks")
    theme.stat_grid([
        ("Adduct", top.get("adduct") or "—"),
        ("Tree score", _fmt(top.get("tree_score"), 2)),
        ("Isotope score", _fmt(top.get("iso_score"), 2)),
        ("Explained peaks", "—" if peaks_expl is None else f"{peaks_expl:.0f}"),
        ("Explained intensity", "—" if explained_int is None else f"{explained_int:.0%}"),
        ("Precursor error", "—" if mass_error is None else f"{mass_error:+.2f} ppm"),
    ], wide=True, max_cols=3)
    if len(cands) > 1:
        with st.expander(f"SIRIUS candidates (top {len(cands)})"):
            st.dataframe(pd.DataFrame([{
                "Rank": c.get("rank"),
                "Formula": c["formula"],
                "Adduct": c.get("adduct") or "—",
                "SIRIUS score": c.get("score"),
            } for c in cands]), width="stretch", hide_index=True)

    adduct = top.get("adduct") or ""
    if engine.ion_mz(top["formula"], adduct) is None:
        adduct = result["adduct_explain"]
    return top["formula"], adduct, "SIRIUS top hit"


def _render_fragments(peaks, formula, adduct, source, precursor_mz, ms2_tol_ppm):
    theme.section("MS² fragment formulas")
    try:
        explained = engine.explain_fragments(
            tuple(p[0] for p in peaks), formula, adduct, ms2_tol_ppm)
    except Exception as e:
        st.warning(f"Fragment formulas could not be assigned: {e}")
        return
    rows = engine.build_annotation_rows(peaks, explained, precursor_mz)

    n_ann = sum(1 for r in rows if r["annotated"])
    total = sum(r["intensity"] for r in rows)
    explained_int = sum(r["intensity"] for r in rows if r["annotated"])
    theme.metric_cards([
        ("Peaks", str(len(rows))),
        ("Explained", f"{n_ann}/{len(rows)}"),
        ("Explained intensity", f"{explained_int / total:.0%}" if total else "—"),
        ("Adduct", adduct),
        ("MS² tolerance", f"{ms2_tol_ppm:g} ppm"),
    ])
    theme.chips([
        ("Explained", theme.PLOT_DARK),
        ("Unexplained", theme.PLOT_LIGHT),
    ])
    theme.note(
        f"Fragment-ion formulas of {formula} {adduct} ({source}). Labels give the m/z of the "
        "strongest peaks and, when explained, the ion formula — hover near any peak for its "
        "formula, neutral loss and mass error; drag to zoom, double-click to reset."
    )

    fig = engine.build_figure(rows, precursor_mz, adduct=adduct,
                              dark=theme.PLOT_DARK, light=theme.PLOT_LIGHT, ink=theme.INK,
                              muted=theme.MUTED, grid=theme.SURFACE_2,
                              highlight=theme.BLUE_BRIGHT)
    st.plotly_chart(fig, width="stretch", theme=None, config={"displayModeBar": False})

    with st.expander("Fragment formula table"):
        st.dataframe(
            pd.DataFrame([{
                "m/z": r["mz"],
                "Rel. int. (%)": r["rel"],
                "Intensity": r["intensity"],
                "Explained": "✓" if r["annotated"] else "✗",
                "Ion formula": r["formula"] + (" (precursor)" if r["precursor"] else ""),
                "Δ ppm": r["ppm"],
                "Neutral loss": r["nl"],
                "Other fits": ", ".join(r["alternatives"]),
            } for r in rows]),
            width="stretch", hide_index=True,
            column_config={
                "m/z": st.column_config.NumberColumn(format="%.4f"),
                "Rel. int. (%)": st.column_config.NumberColumn(format="%.1f"),
                "Δ ppm": st.column_config.NumberColumn(format="%+.1f"),
                "Neutral loss": st.column_config.NumberColumn(format="%.4f"),
            })
