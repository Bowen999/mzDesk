"""BUDDY molecular formula annotation engine (Structure Annotation tab).

Wraps `msbuddy` (Xing et al., Nat. Methods 2023) to provide:
  - top-5 molecular formula candidates for a precursor m/z (+ optional MS¹
    isotope pattern, MS² spectrum, retention time, adduct);
  - per-peak MS² fragment-formula explanation for ANY chosen formula via
    `msbuddy.assign_subformula` (the subformula enumeration BUDDY uses
    internally, applied here to the raw user peak list);
  - the theoretical ion m/z of a formula + adduct (consistency checks, mass errors);
  - the interactive MS² fragment-formula figure.

msbuddy 0.3.9 behaviour relied upon (checked against the installed source):
  - `annotate_formula` drops precursors outside 1 < m/z < 1500 ("No data loaded.");
  - `timeout_secs` is only enforced in parallel (batch) mode, single queries run untimed;
  - retention time is metadata only (seconds); it is not an ML feature;
  - `FormulaResult.mass_error_ppm` is (theoretical - observed) / observed; this module
    reports the conventional (observed - theoretical) / theoretical, like the MS² tab;
  - `assign_subformula` never explains the intact precursor ion (its neutral loss is
    empty), so peaks at the precursor m/z are matched here directly;
  - `Msbuddy` is a process-wide singleton that holds the query data, so runs are serialised.

The heavy msbuddy calls are cached (functools.lru_cache, so their arguments must be hashable:
peak lists as tuples of (mz, intensity) tuples); callers get a copy they may change.
"""
import copy
import math
import threading
from functools import lru_cache

import numpy as np
import plotly.graph_objects as go
from mzdesk.shared.spectrum import parse_spectrum  # shared peak-list parser
from msbuddy import Msbuddy, MsbuddyConfig, assign_subformula, form_arr_to_str, read_formula
from msbuddy.base import Adduct, MetaFeature, Spectrum, check_adduct

ADDUCTS_POSITIVE = ["[M+H]+", "[M+Na]+", "[M+NH4]+", "[M+K]+", "[M+H-H2O]+"]
ADDUCTS_NEGATIVE = ["[M-H]-", "[M+Cl]-", "[M+HCOO]-", "[M+CH3COO]-", "[M+Br]-"]

# instrument -> (ms1 tol ppm, ms2 tol ppm), per msbuddy defaults
INSTRUMENT_TOLS = {"orbitrap": (5, 10), "qtof": (10, 20), "fticr": (2, 5)}

MAX_CANDIDATES = 5

# msbuddy only annotates precursors with 1 < m/z < 1500
PRECURSOR_MZ_RANGE = (1.0, 1500.0)

# msbuddy's element order with monoisotopic masses (float64 -> sub-ppm theoretical m/z)
_ALPHABET = ("C", "H", "Br", "Cl", "F", "I", "K", "N", "Na", "O", "P", "S")
_MONO = np.array([12.0, 1.00782503207, 78.9183371, 34.96885268, 18.99840322, 126.904473,
                  38.96370668, 14.0030740048, 22.9897692809, 15.99491461956, 30.97376163,
                  31.97207100])
_ELECTRON = 0.00054857990946

# Msbuddy is a singleton holding the loaded query: one annotation at a time
_BUDDY_LOCK = threading.Lock()

# plot styling, shared with the MS² tab's spectrum (the tab passes the theme.py colours)
from mzdesk.shared.colors import PLOT_FONT as _FONT, PLOT_MONO as _MONO_FONT   # app fonts
_Y_FLOOR = 1.6    # sticks start slightly below 0 % so even tiny peaks stay visible
_Y_TOP = 116.0    # headroom for the peak labels above the base peak


def default_adduct(positive):
    """Adduct BUDDY assumes when the Adduct input is left on Auto."""
    return "[M+H]+" if positive else "[M-H]-"


def parse_peak_list(text):
    """Parse a peak list into a sorted [(mz, intensity), ...] list.

    Delegates to the shared MS² parser, then enforces the cross-tab contract
    (finite values, m/z > 0, intensity > 0) before anything reaches msbuddy.
    """
    return [(float(mz), float(inten)) for mz, inten in parse_spectrum(text)
            if math.isfinite(mz) and math.isfinite(inten) and mz > 0 and inten > 0]


def _safe_read_formula(formula_str):
    """msbuddy's formula parser, hardened for live validation of arbitrary user text.

    read_formula() raises (OverflowError) on huge element counts, silently wraps counts
    that overflow its int16 array, and returns an all-zero array for strings with no
    recognisable element ("**bold**"), so guard all three here.
    """
    s = (formula_str or "").strip()
    if not s or len(s) > 200:
        return None
    try:
        arr = read_formula(s)
    except Exception:
        return None
    if arr is None:
        return None
    arr = np.asarray(arr)
    return None if (arr < 0).any() or not arr.any() else arr


def formula_is_valid(formula_str):
    """True if the string parses as a CHNOPSFClBrINaK formula understood by msbuddy."""
    return _safe_read_formula(formula_str) is not None


def _parse_adduct(adduct):
    """msbuddy Adduct for a parseable adduct string (e.g. '[M+H]+'), else None."""
    if not adduct:
        return None
    try:
        ok, pos_mode = check_adduct(adduct)
        return Adduct(adduct, pos_mode, report_invalid=True) if ok else None
    except Exception:  # msbuddy raises ValueError / TypeError on malformed strings
        return None


def ion_composition(formula, adduct):
    """(element counts in msbuddy order, charge) of neutral `formula` observed as `adduct`.

    Returns (None, None) if the formula or adduct cannot be parsed, or if the adduct
    removes atoms the formula does not have (e.g. [M-H]- of a formula without H).
    """
    arr = _safe_read_formula(formula)
    ion = _parse_adduct(adduct)
    if arr is None or ion is None:
        return None, None
    counts = (np.asarray(arr, dtype=np.int64) * int(ion.m)
              + np.asarray(ion.net_formula.array, dtype=np.int64))
    if (counts < 0).any():
        return None, None
    return counts, int(ion.charge)


def _mz_from_counts(counts, charge):
    return (float(np.dot(np.asarray(counts, dtype=float), _MONO)) - charge * _ELECTRON) / abs(charge)


def ion_mz(formula, adduct):
    """Theoretical m/z of neutral `formula` observed as `adduct` (None if either is invalid)."""
    counts, charge = ion_composition(formula, adduct)
    return None if counts is None else _mz_from_counts(counts, charge)


def ion_mz_and_formula(formula, adduct):
    """(theoretical m/z, ion formula with the charge sign) of neutral `formula` as `adduct`; (None, None) if the adduct does not apply.

    Unlike ion_mz this also handles [M]+ / [M]-: msbuddy reads [M]+ as [M+H]+, so the radical / permanent ion is done here.
    The ion formula carries only the sign (a doubly charged adduct is still written "+" / "-").
    """
    if adduct in ("[M]+", "[M]-"):
        arr = _safe_read_formula(formula)
        if arr is None:
            return None, None
        charge = 1 if adduct == "[M]+" else -1
        return _mz_from_counts(arr, charge), form_arr_to_str(np.asarray(arr, dtype=np.int16)) + adduct[-1]
    counts, charge = ion_composition(formula, adduct)
    if counts is None:
        return None, None
    return _mz_from_counts(counts, charge), form_arr_to_str(counts.astype(np.int16)) + ("+" if charge > 0 else "-")


def ppm_error(observed, theoretical):
    """Conventional mass error (observed - theoretical) / theoretical, in ppm."""
    return (observed - theoretical) / theoretical * 1e6


def _float(value):
    return None if value is None else float(value)


def _with_mass_error(candidates, precursor_mz, adduct):
    for cand in candidates:
        theo = ion_mz(cand["formula"], adduct)
        cand["mz_theo"] = theo
        cand["ppm"] = None if theo is None else ppm_error(precursor_mz, theo)
    return candidates


def annotate_formula(precursor_mz, positive, adduct, rt, formula_input,
                     ms2_peaks, ms1_peaks, ms_instr="orbitrap", halogen=False):
    """Run BUDDY annotation for one query.

    Parameters map 1:1 to the tab inputs; ``ms2_peaks`` / ``ms1_peaks`` are
    tuples of (mz, intensity) pairs (or empty tuples); ``adduct=None`` means Auto
    ([M+H]+ / [M-H]-); ``rt`` is in seconds and stored as metadata only. If
    ``formula_input`` is a non-empty formula string it is used as the sole
    candidate and BUDDY's ML ranking is skipped entirely.

    Returns a dict:
      candidates  : [{rank, formula, prob, normed_prob, fdr, isotope_sim, mz_theo, ppm}]
      adduct_used : adduct the candidates refer to
      from_input  : bool
      ms1_used    : MS¹ isotope peaks msbuddy matched to the precursor (0 = pattern unused,
                    None = not applicable)
    Raises ValueError for precursors outside msbuddy's supported m/z range.
    """
    return copy.deepcopy(_annotate_formula(precursor_mz, positive, adduct, rt, formula_input,
                                           tuple(map(tuple, ms2_peaks)), tuple(map(tuple, ms1_peaks)),
                                           ms_instr, halogen))


@lru_cache(maxsize=32)
def _annotate_formula(precursor_mz, positive, adduct, rt, formula_input,
                      ms2_peaks, ms1_peaks, ms_instr, halogen):
    adduct_used = adduct or default_adduct(positive)
    if formula_input:
        cands = [{"rank": 1, "formula": formula_input, "prob": None,
                  "normed_prob": None, "fdr": None, "isotope_sim": None}]
        return {"candidates": _with_mass_error(cands, precursor_mz, adduct_used),
                "adduct_used": adduct_used, "from_input": True, "ms1_used": None}

    lo, hi = PRECURSOR_MZ_RANGE
    if not lo < float(precursor_mz) < hi:
        raise ValueError(f"BUDDY supports precursor m/z between {lo:.0f} and {hi:.0f}.")

    def _spectrum(peaks):
        if not peaks:
            return None
        return Spectrum(mz_array=np.array([p[0] for p in peaks], dtype=float),
                        int_array=np.array([p[1] for p in peaks], dtype=float))

    # An MS¹ pattern whose peaks all miss the precursor cannot contribute isotope similarity,
    # but msbuddy would still switch to its "with MS1" model and score that similarity as 0.
    ms1_tol_ppm = INSTRUMENT_TOLS.get(ms_instr, (5, 10))[0]
    if ms1_peaks and not any(abs(ppm_error(mz, float(precursor_mz))) <= ms1_tol_ppm
                             for mz, _ in ms1_peaks):
        ms1_peaks = ()

    mf = MetaFeature(
        identifier=0,
        mz=float(precursor_mz),
        charge=1 if positive else -1,
        rt=(float(rt) if rt is not None else None),
        adduct=adduct,
        ms1=_spectrum(ms1_peaks),
        ms2=_spectrum(ms2_peaks),
    )
    with _BUDDY_LOCK:
        engine = Msbuddy(MsbuddyConfig(ms_instr=ms_instr, halogen=bool(halogen)))
        engine.add_data([mf])
        engine.annotate_formula()
        meta = engine.data[0]
        engine.clear_data()

    cands = []
    for i, cand in enumerate(list(meta.candidate_formula_list or [])[:MAX_CANDIDATES]):
        cands.append({
            "rank": i + 1,
            "formula": str(cand.formula),
            "prob": _float(cand.estimated_prob),
            "normed_prob": _float(cand.normed_estimated_prob),
            "fdr": _float(cand.estimated_fdr),
            "isotope_sim": _float(cand.ms1_isotope_similarity),
        })
    adduct_used = str(meta.adduct.string)
    return {
        "candidates": _with_mass_error(cands, precursor_mz, adduct_used),
        "adduct_used": adduct_used,
        "from_input": False,
        "ms1_used": len(meta.ms1_processed) if meta.ms1_processed is not None else 0,
    }


@lru_cache(maxsize=1024)
def _counts_from_formula(formula):
    """{'C': 8, 'H': 10, ...} from a Hill-style formula string (charge signs ignored)."""
    counts, el, num = {}, "", ""
    for ch in (formula or "") + "|":
        if ch.isupper() or ch == "|":
            if el:
                counts[el] = counts.get(el, 0) + (int(num) if num else 1)
            el, num = ("" if ch == "|" else ch), ""
        elif ch.islower():
            el += ch
        elif ch.isdigit():
            num += ch
    return counts


def _is_even_electron(ion_formula):
    """True when a singly charged ion of this composition is even-electron.

    2·RDBE = 2C - (H + halogens) + (N + P) + 2 is odd for even-electron ions of either
    polarity and even for odd-electron (radical) ions. CID fragments are almost always
    even-electron, so those are preferred when several compositions fit one peak.
    """
    c = _counts_from_formula(ion_formula)
    halogens = c.get("F", 0) + c.get("Cl", 0) + c.get("Br", 0) + c.get("I", 0)
    return (2 * c.get("C", 0) - (c.get("H", 0) + halogens)
            + c.get("N", 0) + c.get("P", 0)) % 2 == 1


def explain_fragments(ms2_mz_tuple, precursor_formula, adduct, ms2_tol_ppm=10.0,
                      dbe_cutoff=0.0):
    """Assign fragment-ion formulas to raw MS² peaks for one precursor formula.

    Returns {peak_index: {"formula": str, "ppm": float, "precursor": bool,
                          "alternatives": [(formula, ppm), ...]}} for explained peaks.
    ppm = (observed - theoretical) / theoretical. When several subformulas fit within
    the tolerance, even-electron compositions are preferred and the smallest |ppm| wins
    among them; the rest are listed as alternatives. Peaks at the precursor m/z are
    matched to the intact precursor ion.
    """
    return copy.deepcopy(_explain_fragments(tuple(ms2_mz_tuple), precursor_formula, adduct,
                                            float(ms2_tol_ppm), float(dbe_cutoff)))


@lru_cache(maxsize=32)
def _explain_fragments(ms2_mz_tuple, precursor_formula, adduct, ms2_tol_ppm, dbe_cutoff):
    mzs = [float(m) for m in ms2_mz_tuple]
    counts, charge = ion_composition(precursor_formula, adduct)
    if not mzs or counts is None:
        return {}
    tol = float(ms2_tol_ppm)

    explained = {}
    precursor_ion = form_arr_to_str(counts.astype(np.int16))
    precursor_theo = _mz_from_counts(counts, charge)
    for k, mz in enumerate(mzs):
        err = ppm_error(mz, precursor_theo)
        if abs(err) <= tol:
            explained[k] = {"formula": precursor_ion, "ppm": err,
                            "precursor": True, "alternatives": []}

    results = assign_subformula(mzs, precursor_formula.strip(), adduct,
                                ms2_tol=tol, ppm=True, dbe_cutoff=float(dbe_cutoff))
    for res in results or []:
        k = int(res.idx)
        if k in explained or not res.subform_list:
            continue
        fits = []
        for sub in res.subform_list:
            sub_counts = read_formula(str(sub.formula))
            if sub_counts is not None:
                fits.append((str(sub.formula),
                             ppm_error(mzs[k], _mz_from_counts(sub_counts, charge))))
        if fits:
            fits.sort(key=lambda fit: (not _is_even_electron(fit[0]), abs(fit[1])))
            explained[k] = {"formula": fits[0][0], "ppm": fits[0][1],
                            "precursor": False, "alternatives": fits[1:]}
    return explained


def build_annotation_rows(peaks, explained, precursor_mz=None):
    """Merge raw peaks with explain_fragments output into display rows."""
    max_inten = max((i for _, i in peaks), default=0.0) or 1.0
    rows = []
    for idx, (mz, inten) in enumerate(peaks):
        hit = explained.get(idx)
        nl = None
        if precursor_mz is not None and mz <= precursor_mz + 0.02:
            nl = max(precursor_mz - mz, 0.0)   # peaks above the precursor have no neutral loss
        rows.append({
            "mz": mz,
            "intensity": inten,
            "rel": inten / max_inten * 100.0,
            "annotated": hit is not None,
            "formula": hit["formula"] if hit else "",
            "ppm": hit["ppm"] if hit else None,
            "precursor": bool(hit and hit["precursor"]),
            "alternatives": [f for f, _ in hit["alternatives"]] if hit else [],
            "nl": nl,
        })
    return rows


def _label_rows(rows, x_range, max_labels=10, min_rel=3.0, min_sep_frac=0.09):
    """Strongest peaks whose labels stay >= 9 % of the axis apart (labels never overlap)."""
    min_sep = (x_range[1] - x_range[0]) * min_sep_frac
    chosen = []
    for row in sorted(rows, key=lambda r: -r["rel"]):
        if row["rel"] < min_rel or len(chosen) >= max_labels:
            break
        if all(abs(row["mz"] - c["mz"]) >= min_sep for c in chosen):
            chosen.append(row)
    return chosen


def _hover_text(r):
    nl = "—" if r["nl"] is None else f"{r['nl']:.4f}"
    if not r["annotated"]:
        return (f"<b>m/z {r['mz']:.4f}</b><br>rel. int {r['rel']:.1f}% · NL {nl}"
                "<br>no fragment formula within tolerance")
    text = (f"<b>{r['formula']}</b>{' · precursor ion' if r['precursor'] else ''}"
            f"<br>m/z {r['mz']:.4f} · {r['ppm']:+.1f} ppm<br>NL {nl} · rel. int {r['rel']:.1f}%")
    alts = r["alternatives"]
    if alts:
        text += "<br>also fits: " + ", ".join(alts[:3]) + (
            f" +{len(alts) - 3} more" if len(alts) > 3 else "")
    return text


def build_figure(rows, precursor_mz=None, adduct="", dark="#1B3A8F", light="#A9BEEA",
                 ink="#0A1633", muted="#5C6B8A", grid="#E9EDF7", highlight="#3B5BDB",
                 height=460):
    """MS² stick spectrum in the MS² tab's visual language.

    dark = explained (fragment formula assigned), light = unexplained. Sticks are lines
    of constant pixel width (from slightly below 0 %), so neighbouring peaks never merge
    when zooming and every peak stays visible zoomed out. One invisible marker per peak
    drives the hover label; hovermode "x" snaps to the stem nearest the cursor (within
    16 px) at any height, and a spike line highlights it. The x-range spans every peak
    and the precursor (dotted marker) with 5 % padding (at least 5 Da).
    """
    fig = go.Figure()
    if not rows:
        return fig
    xs = [r["mz"] for r in rows]
    lo, hi = min(xs), max(xs)
    if precursor_mz:
        lo, hi = min(lo, precursor_mz), max(hi, precursor_mz)
    pad = max((hi - lo) * 0.05, 5.0)
    x_range = [max(0.0, lo - pad), hi + pad]

    for name, color, group in (("unexplained", light, [r for r in rows if not r["annotated"]]),
                               ("explained", dark, [r for r in rows if r["annotated"]])):
        if not group:
            continue
        sx, sy = [], []
        for r in group:
            sx += [r["mz"], r["mz"], None]
            sy += [-_Y_FLOOR, r["rel"], None]
        fig.add_trace(go.Scatter(x=sx, y=sy, mode="lines", name=name, hoverinfo="skip",
                                 line=dict(color=color, width=2)))
    fig.add_trace(go.Scatter(x=xs, y=[r["rel"] for r in rows], mode="markers", name="peaks",
                             marker=dict(size=10, color="rgba(0,0,0,0)"),
                             hovertext=[_hover_text(r) for r in rows],
                             hovertemplate="%{hovertext}<extra></extra>"))

    labelled = _label_rows(rows, x_range)
    annotations, shapes = [], []
    for r in labelled:
        text = f"{r['mz']:.4f}"
        if r["annotated"]:
            text += f"<br><span style='color:{dark}'>{r['formula']}</span>"
        annotations.append(dict(x=r["mz"], y=r["rel"], xref="x", yref="y", yanchor="bottom",
                                yshift=3, showarrow=False, text=text,
                                font=dict(family=_MONO_FONT, size=10,
                                          color=ink if r["annotated"] else muted)))
    if precursor_mz and x_range[0] <= precursor_mz <= x_range[1]:
        # dotted marker above the precursor peak (and its label), or full height if absent
        prec = [r for r in rows if r["precursor"]]
        if prec:
            top = max(prec, key=lambda r: r["rel"])
            y_start = top["rel"] + (13.0 if top in labelled else 4.0)
        else:
            y_start = -_Y_FLOOR
        if y_start < _Y_TOP:
            shapes.append(dict(type="line", xref="x", yref="y", x0=precursor_mz,
                               x1=precursor_mz, y0=y_start, y1=_Y_TOP, layer="below",
                               line=dict(color=muted, width=1, dash="dot")))
        frac = (precursor_mz - x_range[0]) / (x_range[1] - x_range[0])
        annotations.append(dict(
            x=precursor_mz, y=1, xref="x", yref="paper", yanchor="bottom", yshift=4,
            showarrow=False, text=f"precursor {adduct} {precursor_mz:.4f}".replace("  ", " "),
            xanchor="right" if frac > 0.8 else ("left" if frac < 0.2 else "center"),
            font=dict(family=_MONO_FONT, size=10, color=muted)))

    axis = dict(showline=True, linecolor=ink, linewidth=1.5, ticks="outside", tickcolor=ink,
                ticklen=5, tickfont=dict(color=muted, size=11), zeroline=False)
    fig.update_layout(
        height=height, margin=dict(l=62, r=22, t=30, b=50),
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family=_FONT, color=ink, size=12),
        showlegend=False, hovermode="x", hoverdistance=16, dragmode="zoom",
        hoverlabel=dict(bgcolor="#FFFFFF", bordercolor=ink, align="left",
                        font=dict(family=_FONT, size=12, color=ink)),
        shapes=shapes, annotations=annotations,
        xaxis=dict(axis, title=dict(text="m/z", font=dict(color=muted, size=12)),
                   range=x_range, showgrid=False, hoverformat=".4f",
                   showspikes=True, spikemode="toaxis", spikesnap="data",
                   spikecolor=highlight, spikethickness=3, spikedash="solid"),
        yaxis=dict(axis, title=dict(text="Relative intensity (%)", font=dict(color=muted, size=12)),
                   range=[-_Y_FLOOR, _Y_TOP], tickvals=[0, 20, 40, 60, 80, 100],
                   showgrid=True, gridcolor=grid),
    )
    return fig
