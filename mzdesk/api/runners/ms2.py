"""MS2 小类：清洗与可视化，库检索（MassBank / GNPS 本地索引），相似度 + 镜像图（matchms / ms_entropy / 自定义），碎片解释（MetFrag）。"""
import copy
import math
from functools import lru_cache

from mzdesk.api import common
from mzdesk.api.common import EngineError, InputError, fmt, signed
from mzdesk.shared import ms2_clean, ms2_formats, ms2_plot
from mzdesk.shared.spectrum import detect_format, first_record

PRIMARY_METHOD = "SpectralEntropy (weighted)"
FALLBACK_METHOD = "CosineGreedy"


def _finite(score):
    return isinstance(score, (int, float)) and math.isfinite(score)


_Y_TITLE = {"None": "Intensity", "Base peak = 100": "Relative intensity (%)"}


def run_clean(inputs):
    warnings = []
    peaks, file_mz = common.spectrum(inputs, "spectrum", "MS2 spectrum", warnings)
    raw = str((inputs.get("spectrum") or {}).get("value") or "")
    precursor_mz = common.number(inputs, "precursor_mz", "Precursor m/z", positive=True)
    if precursor_mz is None and file_mz:
        precursor_mz = file_mz
        warnings.append(f"Precursor m/z {file_mz:.4f} is read from the spectrum.")
    ion_mode = "Negative" if inputs.get("ion_mode") == "Negative" else "Positive"
    noise = common.number(inputs, "noise_pct", "Noise threshold", default=0.0, minimum=0, maximum=100)
    tol = common.number(inputs, "precursor_tol", "Precursor tolerance", default=0.02, positive=True)
    normalize = inputs.get("normalize") or "None"
    sort_by = inputs.get("sort_by") or "m/z"
    if normalize not in ms2_clean.NORMALIZE or sort_by not in ms2_clean.SORT_BY:
        raise InputError("Unknown normalization or sort order.")
    remove = bool(inputs.get("remove_precursor"))
    if remove and precursor_mz is None:
        warnings.append("Remove precursor is on but there is no precursor m/z; no peak was removed.")

    res = ms2_clean.clean_spectrum(peaks, noise_pct=noise, normalize=normalize, precursor_mz=precursor_mz,
                                   remove_precursor=remove, precursor_tol=tol, sort_by=sort_by)
    rows = res["rows"]
    n_noise = sum(r["reason"] == "noise" for r in rows)
    n_prec = sum(r["reason"] == "precursor" for r in rows)
    if not res["n_out"]:
        warnings.append("Every peak was removed.")
    fig, cards = ms2_plot.spectrum_figure(
        [(r["mz"], r["intensity"]) for r in rows], precursor_mz=precursor_mz, y_title=_Y_TITLE.get(normalize, "Normalized intensity"),
        dim={i: r["reason"] for i, r in enumerate(rows) if not r["kept"]})
    return {"outputs": {
        "summary": {"values": {
            "Detected format": detect_format(raw), "Peaks in": str(res["n_in"]), "Peaks out": str(res["n_out"]),
            "Removed": {"value": str(res["n_in"] - res["n_out"]), "note": f"{n_noise} noise · {n_prec} precursor"}}},
        "spectrum": common.plotly(fig, peaks=cards),
        "peaks": {"rows": [[fmt(r["mz"], 4), ms2_plot.fmt_intensity(r["intensity"]), fmt(r["rel"], 1),
                            "yes" if r["kept"] else "no", r["reason"] or common.MISSING] for r in rows]},
        "formats": {"formats": ms2_formats.all_formats(res["peaks"], precursor_mz=precursor_mz, ion_mode=ion_mode,
                                                       name=first_record(raw)["name"])},
    }, "warnings": warnings}


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
        warnings.append(f"Precursor m/z {file_mz:.4f} is read from the spectrum.")
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

    fig = mf.build_figure(result)
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


# 库谱的网页：MassBank 官方记录页；GNPS 用经典库页面（外部 gnps2 的 /gnpslibraryspectrum 地址现在返回 404）
RECORD_URL = {"MassBank": "https://massbank.eu/MassBank/RecordDisplay?id={}",
              "GNPS": "https://gnps.ucsd.edu/ProteoSAFe/gnpslibraryspectrum.jsp?SpectrumID={}"}
SHORT_METRIC = {"SpectralEntropy (weighted)": "entropy", "SpectralEntropy (unweighted)": "unweighted entropy"}
LIBRARY_PAGE_SIZE = 100


def _names(inputs, key, label, allowed, default):
    value = inputs.get(key)
    if value is None:
        return list(default)
    if isinstance(value, str):
        value = [value]
    if not isinstance(value, (list, tuple)) or not value:
        raise InputError(f"Select at least one {label}.")
    bad = [v for v in value if v not in allowed]
    if bad:
        raise InputError(f"Unknown {label}: {', '.join(map(str, bad))}.")
    return list(value)


def run_library(inputs):
    from urllib.parse import quote

    from mzdesk.ms2 import library_index as lx
    from mzdesk.ms2 import library_search as ls
    from mzdesk.ms2 import metric_names as mn
    from mzdesk.ms2 import mirror_plot
    from mzdesk.spec import LIBRARIES, LIBRARY_DEFAULT_METRICS

    warnings = []
    spectrum_type = inputs.get("spectrum_type") or "MS/MS"
    if spectrum_type not in ls.SPECTRUM_TYPES:
        raise InputError("Spectrum type must be MS/MS or EI.")
    ei = spectrum_type == "EI"
    peaks, file_mz = common.spectrum(inputs, "ms2", "Query spectrum", warnings)
    field_mz, positive, adduct = common.precursor(inputs)
    libraries = _names(inputs, "libraries", "library", LIBRARIES, LIBRARIES)
    metrics = _names(inputs, "metrics", "metric", mn.LIBRARY_ORDER, LIBRARY_DEFAULT_METRICS)
    threshold = common.number(inputs, "threshold", "Score threshold", default=0.7, minimum=-1, maximum=1)
    tol_ppm = common.number(inputs, "tol_ppm", "Precursor tolerance", default=10.0, positive=True, maximum=1000)
    tol_da = common.number(inputs, "tol_da", "Fragment tolerance", default=0.02, positive=True, maximum=5)

    precursor_mz = field_mz or file_mz
    if ei:
        if precursor_mz is not None:
            warnings.append("The precursor m/z is not used for EI spectra; it was ignored.")
        precursor_mz = None
    elif precursor_mz is None:
        raise InputError("Precursor m/z is required for an MS/MS library search. Enter it, or use a spectrum whose MGF / MSP header has it "
                         "(for an EI spectrum, set the spectrum type to EI).")
    elif field_mz is None:
        warnings.append(f"Precursor m/z {file_mz:.4f} is read from the spectrum.")

    try:
        res = ls.search(peaks, spectrum_type, libraries, metrics, threshold, precursor_mz=precursor_mz, positive=positive, adduct=adduct,
                        precursor_tol_ppm=tol_ppm, tolerance=tol_da, top_n=LIBRARY_PAGE_SIZE)
    except lx.LibraryMissing as exc:
        raise EngineError(str(exc)) from None
    except ValueError as exc:
        raise InputError(str(exc)) from None
    except ImportError as exc:
        raise EngineError(f"{exc.name or 'A required package'} is not installed on this server, so the library search cannot run.") from None
    warnings += res["warnings"]
    warnings += [f"{name}: {text}" for name, text in res["notes"].items()]

    used, tol = res["metrics"], res["tolerance"]
    hits = res["hits"]
    rank = res["rank_metric"] or used[0]               # 第一个对命中有分数的所选指标（第一个若全是 N/A 就跳过）
    if hits and rank != used[0]:
        skipped = used[:used.index(rank)]
        warnings.append(f"Ranked by {rank}, the first selected metric that applies ({', '.join(skipped)} "
                        f"{'is' if len(skipped) == 1 else 'are'} N/A for every hit).")
    ms_cols = ["Rank", "Library", "Accession", "Name"] + ([] if ei else ["Adduct", "Precursor m/z", "Δ ppm"]) + ["Matched peaks"]
    columns = ms_cols + used + ["Instrument / CE"]

    def score_cell(v):
        if v is None:
            return common.MISSING
        return {"badge": f"{v:.4f}", "tone": "good"} if v >= threshold else f"{v:.4f}"

    rows = []
    for h in hits:
        acc = h["accession"]
        cells = [str(h["rank"]), h["library"], {"text": acc, "url": RECORD_URL[h["library"]].format(quote(acc, safe=""))} if acc else common.MISSING,
                 h["name"] or common.MISSING]
        if not ei:
            cells += [h["adduct"] or common.MISSING, fmt(h["precursor_mz"], 4), signed(h["ppm"], 1)]
        cells += [str(h["matches"])] + [score_cell(h["scores"][m]) for m in used]
        cells.append(" · ".join(x for x in (h["instrument"], h["collision_energy"]) if x) or common.MISSING)
        rows.append(cells)

    libs = res["libraries"]
    pre = res["prefilter"]
    n_scored = res["n_candidates"]
    if ei:
        scored_note = (f"top {pre['kept']:,} of {pre['total']:,} EI spectra by entropy" if pre and pre["dropped"]
                       else (f"all {pre['total']:,} EI spectra" if pre else "no EI spectra in the libraries"))
    else:
        scored_note = (f"within {tol_ppm:g} ppm" + (f", {adduct}" if adduct else "")) if n_scored else "no spectrum within the precursor window"
        if res["hidden_copies"]:
            scored_note += f" · {res['hidden_copies']:,} GNPS copies of MassBank spectra hidden"
    summary = {
        "Best hit": ({"value": hits[0]["name"] or hits[0]["accession"], "note": f"{hits[0]['library']} {hits[0]['accession']}"} if hits else common.MISSING),
        "Rank score": ({"value": fmt(hits[0]["scores"][rank], 4), "note": rank} if hits and hits[0]["scores"][rank] is not None else common.MISSING),
        "Hits ≥ threshold": {"value": str(res["n_kept"]), "note": f"any of {len(used)} metric{'s' if len(used) != 1 else ''} ≥ {threshold:g}"},
        "Candidates scored": {"value": f"{n_scored:,}", "note": scored_note},
        "Libraries": {"value": " · ".join(l["name"] for l in libs), "note": " · ".join(l["version"] or "unknown version" for l in libs)},
    }
    outputs = {"summary": {"values": summary}, "hits": {"rows": rows, "columns": columns, "count": res["n_kept"]}}
    if hits and res["top_peaks"] is not None:
        top_mz, top_rel = res["top_peaks"]
        label = f"{hits[0]['library']} {hits[0]['accession']}"
        fig = mirror_plot.build_mirror_figure(peaks, list(zip(top_mz.tolist(), top_rel.tolist())), precursor_mz, res["top_precursor"],
                                              tolerance=tol, labels=("Query", label), height=480)
        outputs["mirror"] = common.plotly(fig)
    return {"outputs": outputs, "warnings": warnings}


RUNNERS = {"ms2_clean": run_clean, "ms2_library": run_library, "ms2_similarity": run_similarity, "ms2_fragments": run_fragments}
