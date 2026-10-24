"""MS2 小类：相似度 + 镜像图（matchms / ms_entropy / 自定义），碎片解释（MetFrag）。"""
import copy
import math
from functools import lru_cache

from mzdesk.api import common
from mzdesk.api.common import EngineError, InputError, fmt, signed
from mzdesk.shared import colors

PRIMARY_METHOD = "SpectralEntropy (weighted)"
FALLBACK_METHOD = "CosineGreedy"


def _finite(score):
    return isinstance(score, (int, float)) and math.isfinite(score)


def run_similarity(inputs):
    from mzdesk.ms2 import mirror_plot
    from mzdesk.ms2 import ms2_similarity as sim

    warnings = []
    query, query_file_mz = common.spectrum(inputs, "query", "Query spectrum", warnings)
    ref, ref_file_mz = common.spectrum(inputs, "reference", "Reference spectrum", warnings)
    query_mz = common.number(inputs, "query_precursor", "Query precursor m/z", positive=True) or query_file_mz
    ref_mz = common.number(inputs, "reference_precursor", "Reference precursor m/z", positive=True) or ref_file_mz
    tol = common.number(inputs, "tol_da", "Fragment tolerance", default=0.02, positive=True)
    threshold = common.number(inputs, "threshold", "Score threshold", default=0.7, minimum=-1, maximum=1)
    ionmode = "negative" if inputs.get("ion_mode") == "Negative" else "positive"
    if query_mz is None or ref_mz is None:
        warnings.append("Without both precursor m/z values, modified-cosine, neutral-loss and deep-learning scores are N/A.")

    try:
        results = sim.compute_all_similarities(query, ref, query_mz, ref_mz, tolerance=tol,
                                               use_deep=bool(inputs.get("deep")), ionmode=ionmode)
    except ValueError as exc:
        raise InputError(str(exc)) from None

    by_name = {r["method"]: r for r in results}
    primary = by_name.get(PRIMARY_METHOD)
    if primary is None or not _finite(primary.get("score")):
        primary = by_name.get(FALLBACK_METHOD)
        if primary is not None and _finite(primary.get("score")):
            warnings.append(f"{PRIMARY_METHOD} is unavailable; the primary score is {FALLBACK_METHOD}.")
        else:
            primary = None
    finite = [r for r in results if _finite(r.get("score"))]
    n_pass = sum(1 for r in finite if r["score"] >= threshold)
    summary = mirror_plot.match_summary(query, ref, tol)
    fig = mirror_plot.build_mirror_figure(query, ref, query_mz, ref_mz, tolerance=tol, labels=("Query", "Reference"), height=480)
    return {"outputs": {
        "summary": {"values": {
            "Primary score": {"value": fmt(primary["score"], 4), "note": primary["method"]} if primary else common.MISSING,
            "Metrics ≥ threshold": {"value": f"{n_pass} / {len(finite)}", "note": f"threshold {threshold:g}"},
            "Matched peaks": str(summary["matched_pairs"]),
            "Query matched": {"value": f"{summary['query_matched_pct']:.0f}", "note": f"{summary['query_matched']} / {summary['query_total']} peaks"},
            "Reference matched": {"value": f"{summary['ref_matched_pct']:.0f}", "note": f"{summary['ref_matched']} / {summary['ref_total']} peaks"}}},
        "mirror": common.plotly(fig),
        "scores": {"rows": [[r["category"], r["method"], fmt(r.get("score"), 4),
                             common.MISSING if r.get("matches") is None else str(r["matches"]),
                             ("yes" if r["score"] >= threshold else "no") if _finite(r.get("score")) else common.MISSING,
                             r.get("note") or ""] for r in results]},
    }, "warnings": warnings}


@lru_cache(maxsize=16)
def _metfrag(smiles, peaks, adduct, precursor_mz, tree_depth, tol_abs, tol_ppm):
    from mzdesk.ms2 import metfrag_engine as mf

    try:
        return mf.run_metfrag(smiles, list(peaks), adduct, precursor_mz, tree_depth, tol_abs, tol_ppm)
    except mf.MetFragError as exc:
        raise EngineError(str(exc), exc.log) from None
    except ValueError as exc:
        raise InputError(str(exc)) from None


_KIND = {"precursor": "Precursor ion", "fragment": "MetFrag fragment", "": "Not explained"}


def run_fragments(inputs):
    from mzdesk.ms2 import metfrag_engine as mf

    warnings = []
    smiles = common.structure(inputs, label="Candidate structure")["smiles"]
    precursor_mz, positive, adduct = common.precursor(inputs)
    peaks, file_mz = common.spectrum(inputs, "ms2", "MS2 spectrum", warnings)
    if precursor_mz is None and file_mz:
        precursor_mz = file_mz
        warnings.append(f"Precursor m/z {file_mz:.4f} is taken from the file.")
    adduct = adduct or ("[M+H]+" if positive else "[M-H]-")
    if adduct not in mf.ADDUCTS:
        raise InputError(f"MetFrag does not support {adduct}. Supported: {', '.join(mf.ADDUCTS)}.")
    if mf.ADDUCTS[adduct]["positive"] != positive:
        raise InputError(f"{adduct} does not match the ion mode.")
    depth = int(common.number(inputs, "tree_depth", "Tree depth", default=2, minimum=1, maximum=4))
    tol_abs = common.number(inputs, "tol_abs", "Absolute tolerance", default=0.002, minimum=0)
    tol_ppm = common.number(inputs, "tol_ppm", "Relative tolerance", default=10.0, minimum=0)

    result = copy.deepcopy(_metfrag(smiles, tuple(peaks), adduct, precursor_mz, depth, tol_abs, tol_ppm))
    warnings += result["warnings"]
    ann = result["annotation"]
    total = sum(a["intensity"] for a in ann) or 1.0
    explained = [a for a in ann if a["annotated"]]
    sign = "+" if result["positive"] else "-"

    def card(a):
        rows = []
        if a["ion_formula"]:
            rows.append(["Ion formula", a["ion_formula"]])
        if a["calc_mz"] is not None:
            rows.append(["Theor. m/z", fmt(a["calc_mz"], 4)])
        if a["ppm"] is not None:
            rows.append(["Mass error", f"{signed(a['ppm'], 1)} ppm"])
        if a["kind"] == "fragment" and a["nl"] is not None:
            rows.append(["Neutral loss", fmt(a["nl"], 4)])
        rows.append(["Rel. intensity", f"{fmt(a['rel'], 1)} %"])
        kind = f"Precursor ion · {adduct}" if a["kind"] == "precursor" else _KIND[a["kind"]]
        return {"mz": a["mz"], "rel": a["rel"], "on": a["annotated"], "kind": kind, "rows": rows,
                "img": f"data:image/png;base64,{a['b64']}" if a["b64"] else ""}

    fig = mf.build_figure(result, dark=colors.PLOT_DARK, light=colors.PLOT_LIGHT, ink=colors.INK,
                          muted=colors.MUTED, grid=colors.SURFACE_2)
    return {"outputs": {
        "summary": {"values": {
            "Explained peaks": f"{len(explained)} / {len(ann)}",
            "Explained intensity": f"{sum(a['intensity'] for a in explained) / total * 100:.1f}",
            "Theoretical precursor m/z": {"value": fmt(result["theo_precursor_mz"], 4), "note": adduct},
            "Ion formula": result["ion_formula"] + sign}},
        "spectrum": common.plotly(fig, peaks=[card(a) for a in ann]),
        "peaks": {"rows": [[fmt(a["mz"], 4), fmt(a["rel"], 1), _KIND[a["kind"]], a["ion_formula"] or common.MISSING,
                            fmt(a["calc_mz"], 4), signed(a["ppm"], 1) if a["ppm"] is not None else common.MISSING,
                            fmt(a["nl"], 4) if a["kind"] == "fragment" else common.MISSING, common.png(a["b64"])] for a in ann]},
    }, "warnings": warnings}


RUNNERS = {"ms2_similarity": run_similarity, "ms2_fragments": run_fragments}
