"""MS² Similarity tab.

Compare two MS² spectra across a battery of similarity metrics powered by
`core/ms2_similarity.py`, then visualise them as an interactive mirror plot.

The peak-list input uses the same parser as the MS² Interpretation tab
(``metfrag_engine.parse_spectrum``):

    m/z intensity          (one pair per line, whitespace/comma/tab separated)
    m/z:intensity m/z:intensity ...   (space-separated tokens)

Precursor m/z values are optional; supplying them enables modified-cosine /
neutral-loss metrics and MS2DeepScore.

The result banner always headlines the weighted spectral-entropy similarity
(falling back to CosineGreedy, clearly labelled, if entropy is unavailable) plus
a consensus count of how many scored metrics reach the score threshold.
"""

from __future__ import annotations

import math

import streamlit as st

import theme
import ms2_similarity as engine
from metfrag_engine import parse_spectrum

try:
    from ._ms2_mirror import match_summary, render_mirror
except Exception:  # pragma: no cover - mirror module is authored separately
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
    """Parse a peak list into a sorted ``[(mz, intensity), ...]`` list.

    Uses the app-wide parser (``metfrag_engine.parse_spectrum``): format detected per
    line, header/malformed lines and tokens skipped, non-finite values, m/z <= 0 and
    intensity <= 0 dropped, duplicate m/z collapsed keeping the most intense.
    """
    return parse_spectrum(text)


def _parse_optional_float(value, label):
    """Return a positive float from a non-empty string, or None when blank."""
    value = (value or "").strip()
    if not value:
        return None
    number = _to_float(value)
    if number is None or not math.isfinite(number) or number <= 0:
        raise ValueError(f"{label} must be a positive number or left blank.")
    return number


def _peak_list_note(text):
    """One-line live parse summary shown under a peak-list box."""
    peaks = parse_peaks(text)
    if peaks:
        theme.note(f"{len(peaks)} valid peak{'s' if len(peaks) != 1 else ''} · "
                   f"m/z {peaks[0][0]:.4f}–{peaks[-1][0]:.4f}")
    else:
        theme.note("No valid peaks yet (one “mz intensity” pair per line, or mz:intensity pairs).")


def _is_finite_score(r):
    score = r.get("score")
    return isinstance(score, (int, float)) and math.isfinite(score)


def _fmt_threshold(value):
    """Two decimals, or more when needed (0.999 must not read as 1.00)."""
    return f"{value:.2f}" if round(value, 2) == value else f"{value:g}"


def _ordered_categories(results):
    ordered = []
    for r in results:
        cat = r.get("category", "Other")
        if cat not in ordered:
            ordered.append(cat)
    return ordered


def _primary_result(results):
    """Return ``(row, reason)``: weighted entropy, else CosineGreedy (reason says why)."""
    by_name = {r.get("method"): r for r in results}
    primary = by_name.get(PRIMARY_METHOD)
    if primary is not None and _is_finite_score(primary):
        return primary, None
    reason = (primary.get("note") if primary else "") or "not available"
    fallback = by_name.get(FALLBACK_METHOD)
    if fallback is not None and _is_finite_score(fallback):
        return fallback, reason
    return None, reason


def _render_verdict(results, threshold):
    """Result banner: primary metric + live consensus count at the current threshold."""
    finite = [r for r in results if _is_finite_score(r)]
    n_pass = sum(1 for r in finite if r["score"] >= threshold)
    n_na = len(results) - len(finite)
    thr = _fmt_threshold(threshold)
    consensus = f"{n_pass}/{len(finite)} scored metrics ≥ {thr}"
    if n_na:
        consensus += f" · {n_na} N/A (not computed)"

    primary, reason = _primary_result(results)
    if primary is None:
        theme.note(f"No finite primary score ({PRIMARY_METHOD} / {FALLBACK_METHOD}) — "
                   f"{reason}. Consensus: {consensus}.")
        return primary

    score = primary["score"]
    passed = score >= threshold
    status = f"{'≥' if passed else '<'} threshold {thr}"
    if reason is None:
        lead = f"Primary metric · {status}"
    else:
        lead = f"Fallback primary metric ({PRIMARY_METHOD} unavailable: {reason}) · {status}"
    theme.verdict(
        f"{primary['method']} · {score:.4f}",
        f"{lead} · Consensus {consensus}",
        theme.BLUE if passed else theme.MUTED,
    )
    return primary


def _render_results(state, threshold):
    results = state["results"]
    primary = _render_verdict(results, threshold)
    if not state.get("use_deep"):
        theme.note("Deep-learning scores (MS2DeepScore, Spec2Vec) are off — tick "
                   "Deep scores and compare again to include them.")

    peaks_a, peaks_b = state["peaks_a"], state["peaks_b"]
    precursor_a, precursor_b = state["precursor_a"], state["precursor_b"]
    tolerance = state["tolerance"]

    theme.section("Overview")
    if match_summary is not None:
        ms = match_summary(peaks_a, peaks_b, tolerance)
        theme.metric_cards([
            ("Matched pairs", str(ms["matched_pairs"])),
            (f"Query matched · {ms['query_matched_pct']:.0f}%",
             f"{ms['query_matched']}/{ms['query_total']}"),
            (f"Ref matched · {ms['ref_matched_pct']:.0f}%",
             f"{ms['ref_matched']}/{ms['ref_total']}"),
            ("Tolerance", f"±{tolerance:g} Da"),
        ])
    else:
        theme.metric_cards([
            ("Query peaks", str(len(peaks_a))),
            ("Ref peaks", str(len(peaks_b))),
            ("Tolerance", f"±{tolerance:g} Da"),
        ])
    if precursor_a is not None or precursor_b is not None:
        parts = [f"query {precursor_a:.4f}" if precursor_a is not None else "query —",
                 f"ref {precursor_b:.4f}" if precursor_b is not None else "ref —"]
        if precursor_a is not None and precursor_b is not None:
            parts.append(f"Δ {precursor_b - precursor_a:+.4f} Da")
        theme.note("Precursor m/z — " + " · ".join(parts))

    theme.section("Mirror MS² plot")
    if render_mirror is not None:
        theme.chips([
            ("Matched peak", theme.PLOT_DARK),
            ("Unmatched peak", theme.PLOT_LIGHT),
        ])
        precursor_hint = (" Dotted lines mark the precursor m/z."
                          if precursor_a is not None or precursor_b is not None else "")
        theme.note(
            f"Query points up, reference down. Matched = CosineGreedy peak assignment within "
            f"±{tolerance:g} Da (same count as its score card). Hover near a peak for its "
            "m/z, relative intensity, neutral loss and matched partner; drag to zoom, "
            "double-click to reset." + precursor_hint
        )
        render_mirror(
            peaks_a,
            peaks_b,
            precursor_mz1=precursor_a,
            precursor_mz2=precursor_b,
            tolerance=tolerance,
            height=520,
        )
    else:
        theme.note("Mirror plot module unavailable.")

    theme.section("Similarity scores")
    theme.note("Meter length = |score|; the notch marks the score threshold "
               f"({_fmt_threshold(threshold)}), so a row reaching past it passes. "
               "A dashed meter means the metric was not computed — the reason is under the row.")
    groups = []
    for cat in _ordered_categories(results):
        rows = []
        for r in results:
            if r.get("category") == cat:
                rows.append(dict(r, primary=True) if r is primary else r)
        groups.append((cat, rows))
    theme.similarity_score_cards(groups, threshold=threshold)


def render():
    theme.tab_header(
        "MS² Similarity",
        subtitle=(
            "Compare two MS² spectra across classic peak-matching, spectral-entropy, "
            "custom and (optionally) deep-learning similarity metrics."
        ),
        kicker="similarity · mirror plot · spectral comparison",
    )

    with theme.panel("ms2sim-setup"):
        theme.section("Comparison setup")

        col_query, col_ref = st.columns(2)
        with col_query:
            spec_a_text = st.text_area(
                "Query Spec (peak list)",
                value=DEFAULT_SPEC_A,
                height=139,  # 5 rows: 5*20 + 39 — see "Peak-list heights" in the README
                help="`mz intensity` per line, or space-separated `mz:intensity` pairs.",
                key="ms2sim_spec_a",
            )
            _peak_list_note(spec_a_text)
            precursor_a_text = st.text_input(
                "Query precursor m/z (optional)",
                value="",
                placeholder="e.g. 516.2800",
                key="ms2sim_precursor_a",
            )

        with col_ref:
            spec_b_text = st.text_area(
                "Ref Spec (peak list)",
                value=DEFAULT_SPEC_B,
                height=139,  # 5 rows: 5*20 + 39 — see "Peak-list heights" in the README
                help="`mz intensity` per line, or space-separated `mz:intensity` pairs.",
                key="ms2sim_spec_b",
            )
            _peak_list_note(spec_b_text)
            precursor_b_text = st.text_input(
                "Ref precursor m/z (optional)",
                value="",
                placeholder="e.g. 516.2800",
                key="ms2sim_precursor_b",
            )

        theme.divider()

        c1, c2, c3, c4 = st.columns(4, vertical_alignment="bottom")
        tolerance = c1.number_input(
            "Tolerance (Da)",
            value=0.02,
            min_value=0.001,
            max_value=1.0,
            step=0.005,
            format="%.3f",
            help="Fragment m/z matching tolerance (± Da); also the bin width of the "
                 "binned metrics. Allowed range 0.001–1 Da.",
            key="ms2sim_tolerance",
        )
        threshold = c2.number_input(
            "Score threshold",
            value=0.70,
            min_value=-1.0,
            max_value=1.0,
            step=0.05,
            format="%.2f",
            help="Scores at or above this value are emphasised and counted in the "
                 "consensus. Updates the results live.",
            key="ms2sim_threshold",
        )
        ionmode = c3.selectbox(
            "Ion mode",
            ["positive", "negative"],
            format_func=str.capitalize,
            help="Passed to MS2DeepScore. The Spec2Vec model is positive-mode only.",
            key="ms2sim_ionmode",
        )
        use_deep = c4.checkbox(
            "Deep scores",
            value=False,
            help="MS2DeepScore (needs both precursor m/z) and Spec2Vec. Models (~700 MB) "
                 "are downloaded to ~/.cache/ms2_similarity on first use; the first run "
                 "also takes a while to load them.",
            key="ms2sim_deep",
        )

        # The note and the CTA share the last row: alone under the controls the button
        # left two thirds of the panel's bottom line empty.
        c_hint, c_run = st.columns([2.4, 1], vertical_alignment="center")
        with c_hint:
            theme.note(
                "Supplying precursor m/z unlocks modified-cosine and neutral-loss metrics, "
                "and is required for MS2DeepScore."
            )
        run = c_run.button("Compare spectra", type="primary", key="ms2sim_run",
                           width="stretch")

    signature = (spec_a_text, spec_b_text, precursor_a_text.strip(), precursor_b_text.strip(),
                 float(tolerance), ionmode, bool(use_deep))

    if run:
        st.session_state.pop("ms2sim", None)
        try:
            peaks_a = parse_peaks(spec_a_text)
            peaks_b = parse_peaks(spec_b_text)
            if not peaks_a or not peaks_b:
                st.error("Both spectra need at least one valid peak.")
                return

            precursor_a = _parse_optional_float(precursor_a_text, "Query precursor m/z")
            precursor_b = _parse_optional_float(precursor_b_text, "Ref precursor m/z")

            spinner = ("Computing similarities (the first deep-learning run loads the "
                       "models and can take a while)…" if use_deep else "Computing similarities…")
            with st.spinner(spinner):
                results = engine.compute_all_similarities(
                    peaks_a,
                    peaks_b,
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
            st.error(f"Could not compare the spectra: {exc}")
            return

    state = st.session_state.get("ms2sim")
    if state is None:
        st.info("Enter two spectra in the comparison panel above, then click Compare spectra.")
        return

    theme.divider()
    if state.get("signature") != signature:
        st.warning("The inputs have changed since these results were computed — click "
                   "Compare spectra to update them. (The score threshold applies live.)")
    _render_results(state, threshold)
