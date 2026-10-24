from __future__ import annotations

import math

import streamlit as st

from mzdesk.shared import theme
from mzdesk.ms2 import ms2_similarity as engine
from mzdesk.shared.spectrum import parse_spectrum

try:
    from mzdesk.ms2.mirror_plot import match_summary, render_mirror
except Exception:
    match_summary = None
    render_mirror = None

DEFAULT_SPEC_A = "197.04643:10 255.22775:92 256.23184:31 277.60822:10 283.24814:20"
DEFAULT_SPEC_B = "197.04643:8 255.22775:100 256.23184:28 277.60822:7 283.24814:22"
PRIMARY_METHOD = "SpectralEntropy (weighted)"
FALLBACK_METHOD = "CosineGreedy"


def _to_float(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def parse_peaks(text):
    return parse_spectrum(text)


def _parse_optional_float(value, label):
    value = (value or "").strip()
    if not value:
        return None
    number = _to_float(value)
    if number is None or not math.isfinite(number) or number <= 0:
        raise ValueError(f"{label} must be positive or blank.")
    return number


def _peak_list_note(text):
    peaks = parse_peaks(text)
    if peaks:
        theme.note(f"{len(peaks)} peaks · m/z {peaks[0][0]:.4f}–{peaks[-1][0]:.4f}")
    else:
        theme.note("No valid peaks. Use `m/z intensity` lines or `m/z:intensity` pairs.")


def _is_finite_score(result):
    score = result.get("score")
    return isinstance(score, (int, float)) and math.isfinite(score)


def _fmt_threshold(value):
    return f"{value:.2f}" if round(value, 2) == value else f"{value:g}"


def _ordered_categories(results):
    ordered = []
    for result in results:
        category = result.get("category", "Other")
        if category not in ordered:
            ordered.append(category)
    return ordered


def _primary_result(results):
    by_name = {result.get("method"): result for result in results}
    primary = by_name.get(PRIMARY_METHOD)
    if primary is not None and _is_finite_score(primary):
        return primary, None
    reason = (primary.get("note") if primary else "") or "unavailable"
    fallback = by_name.get(FALLBACK_METHOD)
    if fallback is not None and _is_finite_score(fallback):
        return fallback, reason
    return None, reason


def _render_verdict(results, threshold):
    finite = [result for result in results if _is_finite_score(result)]
    n_pass = sum(1 for result in finite if result["score"] >= threshold)
    n_na = len(results) - len(finite)
    threshold_text = _fmt_threshold(threshold)
    consensus = f"{n_pass}/{len(finite)} metrics ≥ {threshold_text}"
    if n_na:
        consensus += f" · {n_na} N/A"

    primary, reason = _primary_result(results)
    if primary is None:
        theme.note(f"No primary score ({PRIMARY_METHOD} / {FALLBACK_METHOD}): {reason}. "
                   f"{consensus}.")
        return primary

    score = primary["score"]
    passed = score >= threshold
    status = f"{'≥' if passed else '<'} {threshold_text}"
    if reason is None:
        lead = f"Primary · {status}"
    else:
        lead = f"Fallback; {PRIMARY_METHOD} unavailable: {reason} · {status}"
    theme.verdict(
        f"{primary['method']} · {score:.4f}",
        f"{lead} · {consensus}",
        theme.BLUE if passed else theme.MUTED,
    )
    return primary


def _render_results(state, threshold):
    results = state["results"]
    primary = _render_verdict(results, threshold)
    peaks_a, peaks_b = state["peaks_a"], state["peaks_b"]
    precursor_a, precursor_b = state["precursor_a"], state["precursor_b"]
    tolerance = state["tolerance"]

    theme.section("Summary")
    if match_summary is not None:
        summary = match_summary(peaks_a, peaks_b, tolerance)
        theme.metric_cards([
            ("Matched pairs", str(summary["matched_pairs"])),
            (f"Query matched · {summary['query_matched_pct']:.0f}%",
             f"{summary['query_matched']}/{summary['query_total']}"),
            (f"Reference matched · {summary['ref_matched_pct']:.0f}%",
             f"{summary['ref_matched']}/{summary['ref_total']}"),
            ("Tolerance", f"±{tolerance:g} Da"),
        ])
    else:
        theme.metric_cards([
            ("Query peaks", str(len(peaks_a))),
            ("Reference peaks", str(len(peaks_b))),
            ("Tolerance", f"±{tolerance:g} Da"),
        ])

    if precursor_a is not None or precursor_b is not None:
        parts = [f"query {precursor_a:.4f}" if precursor_a is not None else "query —",
                 f"reference {precursor_b:.4f}" if precursor_b is not None else "reference —"]
        if precursor_a is not None and precursor_b is not None:
            parts.append(f"Δ {precursor_b - precursor_a:+.4f} Da")
        theme.note("Precursor m/z: " + " · ".join(parts))

    theme.section("Mirror plot")
    if render_mirror is not None:
        theme.chips([("Matched", theme.PLOT_DARK), ("Unmatched", theme.PLOT_LIGHT)])
        hint = "Query up; reference down. Hover for details; drag to zoom."
        if precursor_a is not None or precursor_b is not None:
            hint += " Dotted lines: precursor m/z."
        theme.note(hint)
        render_mirror(
            peaks_a, peaks_b,
            precursor_mz1=precursor_a,
            precursor_mz2=precursor_b,
            tolerance=tolerance,
            height=520,
        )
    else:
        st.warning("Mirror plot unavailable.")

    theme.section("Scores")
    groups = []
    for category in _ordered_categories(results):
        rows = []
        for result in results:
            if result.get("category") == category:
                rows.append(dict(result, primary=True) if result is primary else result)
        groups.append((category, rows))
    theme.similarity_score_cards(groups, threshold=threshold)


def render():
    with theme.panel("ms2sim-setup"):
        col_query, col_ref = st.columns(2)
        with col_query:
            spec_a_text = st.text_area(
                "Query spectrum",
                value=DEFAULT_SPEC_A,
                height=139,
                help="Lines: `m/z intensity`; or pairs: `m/z:intensity`.",
                key="ms2sim_spec_a",
            )
            _peak_list_note(spec_a_text)
            precursor_a_text = st.text_input(
                "Query precursor m/z (optional)",
                value="",
                placeholder="e.g. 516.2800",
                help="Enables neutral-loss and modified-cosine scores.",
                key="ms2sim_precursor_a",
            )

        with col_ref:
            spec_b_text = st.text_area(
                "Reference spectrum",
                value=DEFAULT_SPEC_B,
                height=139,
                help="Lines: `m/z intensity`; or pairs: `m/z:intensity`.",
                key="ms2sim_spec_b",
            )
            _peak_list_note(spec_b_text)
            precursor_b_text = st.text_input(
                "Reference precursor m/z (optional)",
                value="",
                placeholder="e.g. 516.2800",
                help="Enables neutral-loss and modified-cosine scores.",
                key="ms2sim_precursor_b",
            )

        c1, c2, c3, c4 = st.columns(4, vertical_alignment="bottom")
        tolerance = c1.number_input(
            "Tolerance (Da)", value=0.02, min_value=0.001, max_value=1.0,
            step=0.005, format="%.3f",
            help="Peak matching tolerance and bin width (± Da).",
            key="ms2sim_tolerance",
        )
        threshold = c2.number_input(
            "Score threshold", value=0.70, min_value=-1.0, max_value=1.0,
            step=0.05, format="%.2f",
            help="Pass threshold; updates results live.",
            key="ms2sim_threshold",
        )
        ionmode = c3.selectbox(
            "Ion mode", ["positive", "negative"], format_func=str.capitalize,
            help="Used by MS2DeepScore; Spec2Vec is positive-only.",
            key="ms2sim_ionmode",
        )
        use_deep = c4.checkbox(
            "Deep scores", value=False,
            help="Needs both precursors; downloads ~700 MB on first use.",
            key="ms2sim_deep",
        )
        run = st.button("Compare", type="primary", key="ms2sim_run")

    signature = (spec_a_text, spec_b_text, precursor_a_text.strip(), precursor_b_text.strip(),
                 float(tolerance), ionmode, bool(use_deep))

    if run:
        st.session_state.pop("ms2sim", None)
        try:
            peaks_a = parse_peaks(spec_a_text)
            peaks_b = parse_peaks(spec_b_text)
            if not peaks_a or not peaks_b:
                st.error("Both spectra need a valid peak.")
                return

            precursor_a = _parse_optional_float(precursor_a_text, "Query precursor m/z")
            precursor_b = _parse_optional_float(precursor_b_text, "Reference precursor m/z")
            spinner = "Loading models and computing…" if use_deep else "Computing…"
            with st.spinner(spinner):
                results = engine.compute_all_similarities(
                    peaks_a, peaks_b,
                    precursor_mz1=precursor_a,
                    precursor_mz2=precursor_b,
                    tolerance=tolerance,
                    use_deep=use_deep,
                    ionmode=ionmode,
                )
            st.session_state["ms2sim"] = {
                "results": results,
                "peaks_a": peaks_a,
                "peaks_b": peaks_b,
                "precursor_a": precursor_a,
                "precursor_b": precursor_b,
                "tolerance": tolerance,
                "ionmode": ionmode,
                "use_deep": use_deep,
                "signature": signature,
            }
        except Exception as exc:
            st.error(f"Comparison failed: {exc}")
            return

    state = st.session_state.get("ms2sim")
    if state is None:
        st.info("Enter two spectra, then click Compare.")
        return

    theme.divider()
    if state.get("signature") != signature:
        st.warning("Inputs changed. Compare again to update; threshold changes apply live.")
    _render_results(state, threshold)
