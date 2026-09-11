"""BUDDY molecular formula annotation engine (Structure Annotation tab).

Wraps `msbuddy` (Xing et al., Nat. Methods 2023) to provide:
  - top-5 molecular formula candidates for a precursor m/z (+ optional MS1
    isotopic pattern, MS² spectrum, retention time, adduct)
  - per-peak MS² fragment-formula explanation for ANY chosen formula via
    `msbuddy.assign_subformula` (same bottom-up subformula assignment used
    internally, applied here to the raw user peak list).

All public functions are Streamlit-cached; heavy msbuddy calls only re-run
when their inputs actually change.
"""
import numpy as np
import plotly.graph_objects as go
import streamlit as st
from metfrag_engine import parse_spectrum  # shared peak-list parser (core/ is on sys.path)
from msbuddy import Msbuddy, assign_subformula
from msbuddy.base import MetaFeature, Spectrum

ADDUCTS_POSITIVE = ["[M+H]+", "[M+Na]+", "[M+NH4]+"]
ADDUCTS_NEGATIVE = ["[M-H]-", "[M+Cl]-", "[M+HCOO]-", "[M+CH3COO]-"]

# instrument -> (ms1 tol ppm, ms2 tol ppm), per msbuddy defaults
INSTRUMENT_TOLS = {"orbitrap": (5, 10), "qtof": (10, 20), "fticr": (2, 5)}

MAX_CANDIDATES = 5


def parse_peak_list(text):
    """Parse 'm/z intensity' text into [(mz, intensity), ...] (delegates to the MS² tab parser)."""
    return parse_spectrum(text)


def formula_is_valid(formula_str):
    """True if the string parses as a CHNOPSFClBrINaK formula understood by msbuddy."""
    from msbuddy import read_formula
    return bool(formula_str.strip()) and read_formula(formula_str.strip()) is not None


@st.cache_data(show_spinner=False)
def annotate_formula(precursor_mz, positive, adduct, rt, formula_input,
                     ms2_peaks, ms1_peaks,
                     ms_instr="orbitrap", halogen=False, timeout_secs=300):
    """Run BUDDY annotation for one query.

    Parameters map 1:1 to the tab inputs; ``ms2_peaks`` / ``ms1_peaks`` are
    tuples of (mz, intensity) pairs (or empty tuples). If ``formula_input`` is
    a non-empty formula string it is used as the sole candidate and BUDDY's
    ML ranking is skipped entirely.

    Returns a dict:
      candidates : [{rank, formula, prob, normed_prob, fdr, isotope_sim}]
      adduct_used, from_input(bool), warning(optional str)
    """
    if formula_input:
        return {
            "candidates": [{"rank": 1, "formula": formula_input, "prob": None,
                            "normed_prob": None, "fdr": None, "isotope_sim": None}],
            "adduct_used": adduct,
            "from_input": True,
        }

    charge = 1 if positive else -1
    engine = Msbuddy()
    engine.update_config(
        ms_instr=ms_instr,
        halogen=bool(halogen),
        timeout_secs=float(timeout_secs),
    )
    mf = MetaFeature(
        identifier=0,
        mz=float(precursor_mz),
        charge=charge,
        rt=(float(rt) if rt is not None else None),
        adduct=adduct,
        ms1=Spectrum(mz_array=np.array([p[0] for p in ms1_peaks]),
                     int_array=np.array([p[1] for p in ms1_peaks]))
        if ms1_peaks else None,
        ms2=Spectrum(mz_array=np.array([p[0] for p in ms2_peaks]),
                     int_array=np.array([p[1] for p in ms2_peaks]))
        if ms2_peaks else None,
    )
    engine.add_data([mf])
    engine.annotate_formula()

    meta = engine.data[0]
    candidates = []
    for i, cand in enumerate(meta.candidate_formula_list[:MAX_CANDIDATES]):
        candidates.append({
            "rank": i + 1,
            "formula": str(cand.formula),
            "prob": cand.estimated_prob,
            "normed_prob": cand.normed_estimated_prob,
            "fdr": cand.estimated_fdr,
            "isotope_sim": cand.ms1_isotope_similarity,
        })
    return {
        "candidates": candidates,
        "adduct_used": str(meta.adduct.string) if meta.adduct is not None else adduct,
        "from_input": False,
    }


@st.cache_data(show_spinner=False)
def explain_fragments(ms2_mz_tuple, precursor_formula, adduct, ms2_tol_ppm=10.0,
                      dbe_cutoff=0.0):
    """Assign fragment formulas to raw MS² peaks for one precursor formula.

    Returns {peak_index: {"formula": str, "ppm": float}} for explained peaks.
    """
    results = assign_subformula(
        list(ms2_mz_tuple), precursor_formula, adduct,
        ms2_tol=float(ms2_tol_ppm), ppm=True, dbe_cutoff=float(dbe_cutoff))
    explained = {}
    if not results:
        return explained
    for res in results:
        subs = getattr(res, "subform_list", None) or []
        if not subs:
            continue
        best = min(subs, key=lambda f: abs(f.mass_error_ppm))
        explained[int(res.idx)] = {"formula": str(best.formula),
                                   "ppm": float(best.mass_error_ppm)}
    return explained


def build_annotation_rows(peaks, explained, precursor_mz=None):
    """Merge raw peaks with explain_fragments output into display rows."""
    max_inten = max((i for _, i in peaks), default=1.0)
    rows = []
    for idx, (mz, inten) in enumerate(peaks):
        hit = explained.get(idx)
        rows.append({
            "mz": mz,
            "intensity": inten,
            "rel": inten / max_inten * 100.0,
            "annotated": hit is not None,
            "formula": hit["formula"] if hit else "",
            "ppm": hit["ppm"] if hit else None,
            "nl": (precursor_mz - mz) if precursor_mz else None,
        })
    return rows


def build_figure(rows, dark="#1B3A8F", light="#A9BEEA"):
    """MS² bar figure; dark bars are formula-annotated and hover shows the fragment formula."""
    annotated = [r for r in rows if r["annotated"]]
    unannotated = [r for r in rows if not r["annotated"]]
    if not rows:
        return go.Figure()
    mz_min = min(r["mz"] for r in rows)
    mz_max = max(r["mz"] for r in rows)
    width = max((mz_max - mz_min) * 0.005, 0.4)

    fig = go.Figure()
    if unannotated:
        fig.add_trace(go.Bar(
            name="unannotated",
            x=[r["mz"] for r in unannotated], y=[r["rel"] for r in unannotated],
            marker_color=light, width=width,
            hovertemplate="m/z %{x:.4f}<br>rel. int %{y:.1f}%<extra></extra>"))
    if annotated:
        fig.add_trace(go.Bar(
            name="annotated",
            x=[r["mz"] for r in annotated], y=[r["rel"] for r in annotated],
            marker_color=dark, width=width,
            customdata=[[r["formula"], r["nl"], r["ppm"]] for r in annotated],
            hovertemplate=(
                "<b>%{customdata[0]}</b><br>"
                "m/z %{x:.4f}<br>"
                "NL %{customdata[1]:.4f} · %{customdata[2]:.1f} ppm"
                "<extra>rel. int %{y:.1f}%</extra>")))
    fig.update_layout(
        template="simple_white", bargap=0, showlegend=False, hovermode="closest",
        hoverdistance=20, margin=dict(l=55, r=20, t=15, b=45), height=470,
        xaxis_title="m/z", yaxis_title="Relative intensity (%)",
        xaxis_range=[mz_min - 15, mz_max + 15], yaxis_range=[0, 108])
    return fig
