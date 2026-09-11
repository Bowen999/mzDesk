import pandas as pd
import streamlit as st
import theme
import msbuddy_engine as engine
import sirius_engine


DEFAULT_MS2 = """107.04965 30
123.04422 12
149.02332 55
166.04988 18
177.05410 8
221.06680 40
249.08720 6"""

DEFAULT_ISO = """291.0954 100
292.0987 11
293.1010 1.2"""


def decide_engines(has_ms1, has_ms2):
    """Which annotation engines can run given the current inputs.

    BUDDY: needs only precursor m/z + ion mode (both mandatory fields).
    SIRIUS: additionally requires at least one spectrum (MS1 isotopes or MS2)
    and must be installed.
    """
    decision = {
        "buddy": {"ok": True, "reason": "runs on m/z + ion mode"},
        "sirius": {"ok": False, "reason": ""},
    }
    if not (has_ms1 or has_ms2):
        decision["sirius"]["reason"] = "needs MS1 isotopes and/or MS2 spectrum"
    elif sirius_engine.sirius_available() is None:
        decision["sirius"]["reason"] = "sirius executable not on PATH"
    else:
        decision["sirius"] = {"ok": True, "reason": "spectra available"}
    return decision


def render():
    theme.tab_header(
        "Structure Annotation",
        subtitle=(
            "Annotate the molecular formula of an unknown feature with BUDDY "
            "and SIRIUS - engines are picked automatically from your inputs."
        ),
        kicker="msbuddy · sirius · molecular formula",
    )

    # ---- inputs: full-width compact grid ----
    r1 = st.columns([1.2, 1, 1, 1.4])
    precursor_text = r1[0].text_input(
        "Precursor m/z *",
        value="",
        placeholder="e.g. 291.0954",
        help="Required. Singly charged precursor ion m/z.",
        key="struct_precursor",
    )
    mode = r1[1].segmented_control(
        "Ion mode *",
        ["Positive", "Negative"],
        default="Positive",
        key="struct_mode",
    )
    positive = (mode != "Negative")
    adduct_choice = r1[2].selectbox(
        "Adduct",
        ["Auto"] + (engine.ADDUCTS_POSITIVE if positive else engine.ADDUCTS_NEGATIVE),
        index=0,
        help="Auto: BUDDY uses [M+H]+ / [M-H]-, SIRIUS uses its [M+?]+ wildcard.",
        key="struct_adduct",
    )
    formula_input = r1[3].text_input(
        "Formula (optional)",
        value="",
        placeholder="e.g. C15H16O5 - skips ranking, explains MS2 directly",
        key="struct_formula",
    )

    r2 = st.columns(2)
    ms2_text = r2[0].text_area(
        "MS2 spectrum (optional)",
        value=DEFAULT_MS2,
        height=150,
        help="`mz intensity` per line, OR space-separated `mz:intensity` pairs "
             "(e.g. `141.0148:92 327.23431:5 ...`).",
        key="struct_ms2",
    )
    iso_text = r2[1].text_area(
        "MS1 isotope pattern (optional)",
        value=DEFAULT_ISO,
        height=150,
        help="Isotope peaks of the precursor; used for isotope-similarity scoring. "
             "Same formats as MS2. Either MS1 or MS2 unlocks SIRIUS.",
        key="struct_iso",
    )

    r3 = st.columns([1, 3])
    rt = r3[0].number_input(
        "Retention time (s)",
        min_value=0.0,
        value=0.0,
        step=1.0,
        help="Optional; 0 = not provided.",
        key="struct_rt",
    )
    with r3[1]:
        with st.expander("Advanced"):
            a1, a2 = st.columns(2)
            ms_instr = a1.selectbox(
                "Instrument", list(engine.INSTRUMENT_TOLS.keys()), index=0,
                help="orbitrap 5/10 ppm - qtof 10/20 ppm - fticr 2/5 ppm "
                     "(SIRIUS fticr uses its orbitrap profile).",
                key="struct_instr",
            )
            halogen = a2.checkbox(
                "Include halogens (F Cl Br I)", value=False, key="struct_halogen")
            timeout_secs = a1.number_input(
                "Timeout / query (s)", min_value=30, value=300, step=30, key="struct_timeout")

    # ---- live engine preview based on current inputs ----
    has_ms2 = bool(engine.parse_peak_list(ms2_text))
    has_ms1 = bool(engine.parse_peak_list(iso_text))
    decision = decide_engines(has_ms1, has_ms2)

    theme.section("Engines to run")
    theme.chips([
        ("✓ BUDDY" if decision["buddy"]["ok"] else "✗ BUDDY",
         theme.BLUE if decision["buddy"]["ok"] else "#C0392B"),
        ("✓ SIRIUS" if decision["sirius"]["ok"] else "✗ SIRIUS",
         theme.BLUE if decision["sirius"]["ok"] else "#C0392B"),
    ])
    skipped = [k.upper() + ": " + v["reason"] for k, v in decision.items() if not v["ok"]]
    if skipped:
        theme.note("Skipped - " + "; ".join(skipped) + ".")

    run = st.button("Annotate formula", type="primary", key="struct_run")

    if run:
        _handle_run(precursor_text, positive, adduct_choice, formula_input,
                    ms2_text, iso_text, rt, decision, ms_instr, halogen,
                    timeout_secs)

    result = st.session_state.get("struct_result")

    # ---- engine messages ----
    if result:
        for msg in result["messages"]:
            st.warning(msg)

    if result is None:
        st.info("Enter a precursor m/z above, then click Annotate formula.")
    else:
        buddy = result.get("buddy")
        siri = result.get("sirius")
        chosen = None

        if buddy:
            chosen = _render_buddy(result, buddy)
        if siri:
            _render_sirius(siri)
            if chosen is None:
                # fall back to SIRIUS top hit for fragment plotting
                chosen = (siri["candidates"][0]["formula"],
                          result["adduct_explain"])

        # ---- shared MS2 fragment explanation for whichever formula is active ----
        peaks = result["ms2_peaks"]
        if peaks and chosen:
            active_formula, active_adduct = chosen
            _render_fragments(peaks, active_formula, active_adduct,
                              result["precursor_mz"], result["ms2_tol_ppm"])
        elif not peaks and result is not None:
            theme.note("No MS2 spectrum provided - showing the formula annotation only.")

    theme.divider()
    theme.section("How it works")
    theme.note(
        "BUDDY enumerates candidate formulas for the precursor m/z, shrinks the space "
        "by explaining MS2 peaks bottom-up, ranks candidates with machine learning and "
        "reports an estimated FDR."
    )
    theme.note(
        "SIRIUS scores isotope pattern and fragmentation trees; it only runs when an "
        "MS1 isotope pattern and/or MS2 spectrum is provided (and sirius is installed)."
    )
    theme.note(
        "Provide a Formula to skip ranking - it is fed to both engines directly. "
        "Otherwise Top-1 is shown and can be switched with one click."
    )
    theme.note("Xing et al. BUDDY. Nat Methods 20, 881-890 (2023) · "
               "Dührkop et al. SIRIUS 4. Nat Methods 16, 299-302 (2019).")


def _handle_run(precursor_text, positive, adduct_choice, formula_input,
                ms2_text, iso_text, rt, decision, ms_instr, halogen, timeout_secs):
    error = None
    precursor_mz = None
    try:
        if precursor_text.strip():
            precursor_mz = float(precursor_text)
    except ValueError:
        error = "Precursor m/z must be a number."
    if precursor_mz is None:
        error = error or "Precursor m/z is required."
    if precursor_mz is not None and precursor_mz <= 0:
        error = "Precursor m/z must be positive."

    formula_clean = formula_input.strip()
    if formula_clean and not engine.formula_is_valid(formula_clean):
        error = (f"'{formula_clean}' is not a valid formula "
                 "(supported elements: C H N O P S F Cl Br I Na K).")

    adduct = None if adduct_choice == "Auto" else adduct_choice

    if error:
        st.error(error)
        return
    if not any(d["ok"] for d in decision.values()):
        st.error("No engine can run with the current inputs.")
        return

    ms2_peaks = tuple(map(tuple, engine.parse_peak_list(ms2_text)))
    ms1_peaks = tuple(map(tuple, engine.parse_peak_list(iso_text)))
    combined = {
        "precursor_mz": precursor_mz,
        "positive": positive,
        "ms2_peaks": ms2_peaks,
        "ms2_tol_ppm": engine.INSTRUMENT_TOLS.get(ms_instr, (5, 10))[1],
        "adduct_explain": adduct or ("[M+H]+" if positive else "[M-H]-"),
        "buddy": None,
        "sirius": None,
        "messages": [],
    }

    if decision["buddy"]["ok"]:
        try:
            with st.spinner("Running BUDDY…"):
                buddy = engine.annotate_formula(
                    precursor_mz=precursor_mz,
                    positive=positive,
                    adduct=adduct,
                    rt=(rt if rt > 0 else None),
                    formula_input=formula_clean,
                    ms2_peaks=ms2_peaks,
                    ms1_peaks=ms1_peaks,
                    ms_instr=ms_instr,
                    halogen=halogen,
                    timeout_secs=int(timeout_secs),
                )
            buddy["adduct_explain"] = (
                buddy["adduct_used"] or ("[M+H]+" if positive else "[M-H]-"))
            combined["buddy"] = buddy
            if not buddy["candidates"]:
                combined["messages"].append(
                    "BUDDY returned no candidate formulas - try enabling halogens "
                    "or widening tolerances in Advanced.")
        except Exception as e:
            combined["messages"].append(f"BUDDY failed: {e}")

    if decision["sirius"]["ok"]:
        try:
            with st.spinner("Running SIRIUS…"):
                combined["sirius"] = sirius_engine.run_sirius(
                    precursor_mz=precursor_mz,
                    positive=positive,
                    adduct=(adduct or sirius_engine.auto_adduct(positive)),
                    formula_input=formula_clean,
                    ms2_peaks=ms2_peaks,
                    ms1_peaks=ms1_peaks,
                    ms_instr=ms_instr,
                    timeout_secs=int(timeout_secs) * 2,
                )
        except Exception as e:
            combined["messages"].append(f"SIRIUS failed: {e}")

    if not combined["buddy"] and not combined["sirius"]:
        combined["messages"].append("Both engines failed - no results available.")
    st.session_state["struct_result"] = combined
    st.session_state.pop("struct_selection", None)


def _render_buddy(result, buddy):
    """Render the BUDDY block. Returns (chosen_formula, adduct) for fragment plotting."""
    candidates = buddy["candidates"]
    if not candidates:
        return None

    sel_key = "struct_selection"
    theme.section("BUDDY")
    if buddy["from_input"]:
        labels = [f"Input · {candidates[0]['formula']}"]
    else:
        labels = [f"#{c['rank']} {c['formula']}" for c in candidates]
    prev = st.session_state.get(sel_key)
    default_idx = labels.index(prev) if prev in labels else 0

    if len(candidates) > 1:
        selection = st.segmented_control(
            "Switch candidate", labels, default=labels[default_idx], key=sel_key)
        idx = labels.index(selection) if selection in labels else 0
    else:
        idx = 0

    chosen = candidates[idx]
    fdr = chosen["fdr"]
    sub = "user-provided formula" if buddy["from_input"] else (
        f"rank #{chosen['rank']} of {len(candidates)}"
        + (f" · estimated FDR {fdr:.1%}" if fdr is not None else ""))
    theme.verdict(chosen["formula"], sub, theme.BLUE)

    stats = [("Precursor m/z", f"{result['precursor_mz']:.4f}"),
             ("Ion mode", "Positive" if result["positive"] else "Negative")]
    if not buddy["from_input"]:
        stats += [
            ("ML probability", "—" if chosen["prob"] is None else f"{chosen['prob']:.3f}"),
            ("Est. FDR", "—" if chosen["fdr"] is None else f"{chosen['fdr']:.1%}"),
        ]
    theme.stat_grid(stats, wide=True)

    if len(candidates) > 1:
        def _fmt(v, digits=4):
            return "—" if v is None else f"{v:.{digits}f}"
        with st.expander("All candidates (Top 5)"):
            st.dataframe(pd.DataFrame([{
                "Rank": c["rank"],
                "Formula": c["formula"],
                "Probability": _fmt(c["prob"]),
                "Norm. prob": _fmt(c["normed_prob"]),
                "Est. FDR": _fmt(c["fdr"]),
                "Iso similarity": _fmt(c["isotope_sim"]),
            } for c in candidates]), width="stretch", hide_index=True)

    return chosen["formula"], buddy.get("adduct_explain") or result["adduct_explain"]


def _render_sirius(siri):
    theme.section("SIRIUS")
    cands = siri["candidates"]
    top = cands[0]
    sub = f"top hit · SIRIUS score {top.get('score') or '—'}" + (
        f" · normalized {top['score_norm']}" if top.get("score_norm") else "")
    theme.verdict(top["formula"], sub, "#0E7C9B")

    with st.expander("SIRIUS details"):
        st.dataframe(pd.DataFrame([{
            "Formula": c["formula"],
            "Adduct": c.get("adduct") or "—",
            "SIRIUS score": c.get("score") or "—",
            "Norm.": c.get("score_norm") or "—",
            "Tree score": c.get("tree_score") or "—",
            "Iso score": c.get("iso_score") or "—",
            "Explained peaks": c.get("explained_peaks") or "—",
        } for c in cands]), width="stretch", hide_index=True)


def _render_fragments(peaks, formula, adduct, precursor_mz, ms2_tol_ppm):
    theme.section("MS² fragment formulas")
    theme.chips([
        ("Explained", theme.PLOT_DARK),
        ("Unexplained", theme.PLOT_LIGHT),
    ])
    theme.note("Hover an explained (dark) peak for its fragment formula, neutral loss and mass error.")

    explained = engine.explain_fragments(
        tuple(p[0] for p in peaks), formula, adduct, ms2_tol_ppm)
    rows = engine.build_annotation_rows(peaks, explained, precursor_mz)

    fig = engine.build_figure(rows, dark=theme.PLOT_DARK, light=theme.PLOT_LIGHT)
    st.plotly_chart(fig, width="stretch", config={"displayModeBar": False})

    n_ann = sum(1 for r in rows if r["annotated"])
    theme.metric_cards([
        ("Peaks", str(len(rows))),
        ("Explained", str(n_ann)),
        ("Adduct", adduct),
        ("Tolerance", f"{ms2_tol_ppm:.0f} ppm"),
    ])

    with st.expander("Fragment formula table"):
        st.dataframe(pd.DataFrame([{
            "m/z": r["mz"],
            "Intensity": r["intensity"],
            "Rel. %": round(r["rel"], 1),
            "Explained": "✓" if r["annotated"] else "✗",
            "Formula": r["formula"],
            "Neutral loss": f"{r['nl']:.4f}" if r["nl"] is not None else "",
            "ppm": f"{r['ppm']:.1f}" if r["ppm"] is not None else "",
        } for r in rows]), width="stretch", hide_index=True)
