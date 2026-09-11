"""MS2 Similarity tab.

Compare two MS/MS spectra across a battery of similarity metrics powered by
`core/ms2_similarity.py`, then visualise them as an interactive mirror plot.

The peak-list input accepts the same formats as the MS² Interpretation tab:

    m/z intensity          (one pair per line, whitespace/comma separated)
    m/z:intensity m/z:intensity ...   (space-separated tokens)

Precursor m/z values are optional; supplying them enables modified-cosine /
neutral-loss metrics and the deep-learning scores.
"""

from __future__ import annotations

import re

import streamlit as st

import theme
import ms2_similarity as engine

try:
    from ._ms2_mirror import match_summary, render_mirror
except Exception:  # pragma: no cover - mirror module is authored separately
    match_summary = None
    render_mirror = None


DEFAULT_SPEC_A = "197.04643:10 255.22775:92 256.23184:31 277.60822:10 283.24814:20"

DEFAULT_SPEC_B = "197.04643:8 255.22775:100 256.23184:28 277.60822:7 283.24814:22"


def parse_peaks(text):
    """Parse a peak list into ``[(mz, intensity), ...]`` pairs."""
    text = (text or "").strip()
    if not text:
        return []

    peaks = []
    # Colon form: `mz:intensity mz:intensity ...`
    if ":" in text:
        for tok in text.split():
            bits = tok.split(":")
            try:
                mz = float(bits[0])
            except (ValueError, IndexError):
                continue
            intensity = 1.0
            if len(bits) > 1 and bits[1]:
                try:
                    intensity = float(bits[1])
                except ValueError:
                    intensity = 1.0
            peaks.append((mz, intensity))
    else:
        for line in text.splitlines():
            parts = [p for p in re.split(r"[\s,]+", line.strip()) if p]
            if not parts:
                continue
            try:
                mz = float(parts[0])
            except ValueError:
                continue
            intensity = float(parts[1]) if len(parts) > 1 else 1.0
            peaks.append((mz, intensity))

    # Collapse duplicate m/z (rounded to 5 dp) keeping the strongest intensity.
    aggregated = {}
    for mz, intensity in peaks:
        key = round(mz, 5)
        aggregated[key] = max(aggregated.get(key, 0.0), intensity)
    return sorted(aggregated.items())


def _parse_optional_float(value, label):
    """Return a float from a non-empty string, or None when blank."""
    value = (value or "").strip()
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        raise ValueError(f"{label} must be a number or left blank.")


def _ordered_categories(results):
    ordered = []
    for r in results:
        cat = r.get("category", "Other")
        if cat not in ordered:
            ordered.append(cat)
    return ordered


def _render_results(results, threshold):
    theme.divider()

    scored = [r for r in results if isinstance(r.get("score"), (int, float))]
    if scored:
        best = max(scored, key=lambda r: r["score"])
        theme.verdict(
            best["method"],
            f"Best similarity · {best['score']:.4f}",
            theme.BLUE,
        )
    else:
        theme.note("No finite similarity scores were produced for these spectra.")

    n_na = sum(1 for r in results if r.get("score") is None)
    if n_na:
        theme.note(
            f"{n_na} method(s) were skipped or failed — see the note in each score card. "
            "Deep-learning rows are N/A unless enabled and the required models are installed."
        )

    groups = []
    for cat in _ordered_categories(results):
        groups.append((cat, [r for r in results if r.get("category") == cat]))

    theme.section("Similarity scores")
    theme.similarity_score_cards(groups, threshold=threshold)


def render():
    theme.tab_header(
        "MS2 Similarity",
        subtitle=(
            "Compare two MS/MS spectra across classic peak-matching, spectral-entropy, "
            "custom and (optionally) deep-learning similarity metrics."
        ),
        kicker="similarity · mirror plot · spectral comparison",
    )

    with st.container(border=True):
        theme.section("Comparison setup")

        col_query, col_ref = st.columns(2)
        with col_query:
            theme.section("Query Spec")
            spec_a_text = st.text_area(
                "Query Spec (peak list)",
                value=DEFAULT_SPEC_A,
                height=170,
                help="`mz intensity` per line, or space-separated `mz:intensity` pairs.",
                key="ms2sim_spec_a",
            )
            precursor_a_text = st.text_input(
                "Query precursor m/z (optional)",
                value="",
                placeholder="e.g. 516.2800",
                key="ms2sim_precursor_a",
            )

        with col_ref:
            theme.section("Ref Spec")
            spec_b_text = st.text_area(
                "Ref Spec (peak list)",
                value=DEFAULT_SPEC_B,
                height=170,
                help="`mz intensity` per line, or space-separated `mz:intensity` pairs.",
                key="ms2sim_spec_b",
            )
            precursor_b_text = st.text_input(
                "Ref precursor m/z (optional)",
                value="",
                placeholder="e.g. 516.2800",
                key="ms2sim_precursor_b",
            )

        theme.divider()

        c1, c2, c3, c4, c5 = st.columns([1.05, 1.05, 1, 1, 1.05])
        tolerance = c1.number_input(
            "Tolerance (Da)",
            value=0.02,
            min_value=0.0,
            step=0.01,
            format="%.2f",
            key="ms2sim_tolerance",
        )
        threshold = c2.number_input(
            "Score threshold",
            value=0.70,
            min_value=-1.0,
            max_value=1.0,
            step=0.05,
            format="%.2f",
            key="ms2sim_threshold",
        )
        ionmode = c3.selectbox(
            "Ion mode",
            ["positive", "negative"],
            key="ms2sim_ionmode",
        )
        use_deep = c4.checkbox(
            "Deep scores",
            value=False,
            help="MS2DeepScore / Spec2Vec — may download models on first run.",
            key="ms2sim_deep",
        )
        with c5:
            st.write("")
            run = st.button(
                "Compare spectra",
                type="primary",
                width="stretch",
                key="ms2sim_run",
            )

    theme.note(
        "Supplying precursor m/z unlocks modified-cosine and neutral-loss metrics, "
        "and is required for MS2DeepScore."
    )

    if run:
        try:
            peaks_a = parse_peaks(spec_a_text)
            peaks_b = parse_peaks(spec_b_text)
            if not peaks_a or not peaks_b:
                st.error("Both spectra need at least one valid peak.")
                return

            precursor_a = _parse_optional_float(precursor_a_text, "Query precursor m/z")
            precursor_b = _parse_optional_float(precursor_b_text, "Ref precursor m/z")

            with st.spinner("Computing similarities…"):
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
                    "threshold": threshold,
                }
        except Exception as exc:
            st.error(str(exc))
            return

    if "ms2sim" in st.session_state:
        state = st.session_state["ms2sim"]
        peaks_a = state["peaks_a"]
        peaks_b = state["peaks_b"]
        precursor_a = state["precursor_a"]
        precursor_b = state["precursor_b"]
        tolerance = state["tolerance"]
        results = state["results"]

        theme.divider()
        theme.section("Overview")
        theme.metric_cards([
            ("Query peaks", str(len(peaks_a))),
            ("Ref peaks", str(len(peaks_b))),
            ("Query precursor m/z", f"{precursor_a:.4f}" if precursor_a is not None else "—"),
            ("Ref precursor m/z", f"{precursor_b:.4f}" if precursor_b is not None else "—"),
            ("Tolerance", f"{tolerance:.3f} Da"),
        ])

        if match_summary is not None:
            ms = match_summary(peaks_a, peaks_b, tolerance)
            theme.section("Matched peaks")
            theme.metric_cards([
                ("Matched pairs", str(ms["matched_pairs"])),
                ("Query matched", f"{ms['query_matched']}/{ms['query_total']}"),
                ("Query match %", f"{ms['query_matched_pct']:.1f}%"),
                ("Ref matched", f"{ms['ref_matched']}/{ms['ref_total']}"),
                ("Ref match %", f"{ms['ref_matched_pct']:.1f}%"),
            ])

        theme.section("Mirror MS² plot")
        if render_mirror is not None:
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

        _render_results(results, threshold)
    else:
        st.info("Enter two spectra in the comparison panel above, then click Compare spectra.")
