"""CCS 页：观测 CCS 和本地索引（CCSbase 源数据）里同一化合物、同一加合物的 N2 参考值比较；没有参考值时用本地训练的模型预测。
引擎和比较规则见 mzdesk/orthogonal/ccs.py，数据准备见 ccs_setup.py。

Layout: five metric cards (reference CCS, source, observed, Δ, verdict), the reference values (one row per measured value, with its tolerance and verdict),
the other adducts of the compound (information only), the prediction (only when the index has no value), and the assumptions.
Nothing goes online at run time. The index is not installed: the page still runs and says how to install it.
"""
from mzdesk.api import blocks as B
from mzdesk.api import common
from mzdesk.api.common import EngineError, InputError

MISSING = common.MISSING
TONES = {"Consistent": "good", "Inconsistent": "warn", "Not comparable": "neutral", "No observed value": "neutral", "No reference": "neutral"}


def _num(v, nd=1):
    return f"{v:.{nd}f}"


def _tol(v):
    return f"±{v:.1f}"


def _obs(v):
    """观测值按输入的样子显示（152.37 不被四舍五入成 152.4）。"""
    return f"{v:g}"


# ---------------------------------------------------------------- tables
def _source_cell(res, row):
    url = res.index.source(row).get("url") or ""
    return B.link(row.src_tag, url) if url.startswith(("http://", "https://")) else row.src_tag


def _row_help(res, row, comp):
    from mzdesk.orthogonal import ccs
    lines = [f"Reference: {res.index.source(row).get('reference') or MISSING}",
             f"Name in the source: {row.name}. Adduct {row.adduct}, m/z {_num(row.mz, 4)} as listed."]
    mismatch = ccs.mz_mismatch(row, res.cfg)
    if mismatch:
        lines.append(f"Check this row: the listed m/z ({_num(mismatch[0], 4)}) differs from the calculated m/z ({_num(mismatch[1], 4)}).")
    if comp.delta_pct is not None:
        lines.append(f"Closest observed value: {_obs(comp.observed)} Å². Δ = (observed − reference) / reference.")
    return lines


def _reference_table(res):
    from mzdesk.orthogonal import ccs
    rows = [B.row([_source_cell(res, r), ccs.TYPE_LABEL[r.ccs_type], r.ccs_method or MISSING, _num(r.ccs),
                   common.signed(c.delta_pct, 2) if c.delta_pct is not None else MISSING,
                   _tol(c.tolerance) if c.tolerance is not None else MISSING,
                   B.badge(c.verdict, TONES[c.verdict])], help=_row_help(res, r, c)) for r, c in zip(res.rows, res.comparisons)]
    return B.table_data(rows, count=len(rows))


def _other_ions_table(res):
    from mzdesk.orthogonal import ccs
    by = {}
    for r in res.others:
        by.setdefault(r.adduct, []).append(r)
    rows = []
    for adduct, rs in sorted(by.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        values = sorted(r.ccs for r in rs)
        shown = " / ".join(_num(v) for v in values[:4]) + (" / …" if len(values) > 4 else "")
        types = ", ".join(ccs.TYPE_LABEL[t] for t in sorted({r.ccs_type for r in rs}, key=ccs.TYPE_ORDER.get))
        rows.append([adduct, shown, types, ", ".join(sorted({r.src_tag for r in rs}))])
    return B.table_data(rows, count=len(rows))


def _prediction_table(res):
    p = res.prediction
    if p is None:
        return B.table_data([])
    info = res.meta["model"]
    split = info["split"]
    thr = res.cfg.domain_min_tanimoto
    help_ = [f"{p.model} trained on {split['train_rows'] + split['val_rows']} measured values of {split['trainval_compounds']} compounds (m/z, adduct, 42 MQN descriptors).",
             f"Test error on {split['test_rows']} values of {split['test_compounds']} unseen compounds: median {p.median_pct:.2f} %, 95th percentile {p.q95_pct:.2f} %.",
             f"Nearest training compound: Tanimoto {p.max_tanimoto:.2f}; below {thr:g} the compound is outside the training domain."]
    domain = B.badge("Inside" if p.in_domain else "Outside", "good" if p.in_domain else "warn",
                     help=f"Nearest training compound: Tanimoto {p.max_tanimoto:.2f}; threshold {thr:g}.")
    return B.table_data([B.row([_num(p.ccs), _tol(p.tolerance_pct), f"{p.median_pct:.2f}", f"{p.q95_pct:.2f}", domain], help=help_)], count=1)


# ---------------------------------------------------------------- notes
def _notes(res):
    from mzdesk.orthogonal import ccs
    cfg = res.cfg
    pct = lambda v: f"{v:.1f} %"      # noqa: E731
    items = []
    if not res.installed:
        items.append(f"The CCS reference data is not installed. Run `{ccs.SETUP_COMMAND}`.")
    items += [
        "Compared: the same compound (InChIKey first block) and the same adduct; other adducts are for information.",
        "N₂ values only: He, TWIMS calibrated with He, and FAIMS / DMS are not comparable. For TIMS, enter CCS in Å², not 1/K0.",
        f"Tolerance: same type {pct(cfg.tol_same)}; DT and TIMS {pct(cfg.tol_dt_tims)}; DT and TW {pct(cfg.tol_dt_tw)}; TIMS and TW {pct(cfg.tol_tims_tw)}.",
    ]
    if not res.installed:
        return {"items": items}
    meta, model = res.meta, res.meta.get("model")
    if model:
        test = model["test"]
        items.append(f"Prediction (no value in the index): median test error {test['median_pct']:.2f} %, 95th percentile {test['q95_pct']:.2f} %, "
                     f"tolerance {_tol(model['tolerance_pct'])} %.")
    if res.prediction_note:
        items.append(res.prediction_note)
    items.append(f"Data: CCSbase source data (Ross et al., Anal. Chem. 2020; c3sdb commit {str(meta.get('commit', ''))[:10]}, {meta.get('licence')} licence).")
    return {"items": items}


# ---------------------------------------------------------------- summary cards
def _summary(res):
    from mzdesk.orthogonal import ccs
    used, p = res.used, res.prediction
    observed = {"value": ", ".join(_obs(v) for v in res.observed), "note": res.platform + (", He" if res.gas == "He" and res.obs_type else "")} if res.observed \
        else {"value": MISSING}
    if not res.installed:
        return {"Reference CCS": {"value": MISSING}, "Source": {"value": "Not installed", "note": f"run {ccs.SETUP_COMMAND}"}, "Observed": observed,
                "Δ": {"value": MISSING}, "Verdict": {"value": ccs.NO_REFERENCE, "note": "CCS data not installed"}}
    # reference CCS and where it comes from
    if res.source == "measured":
        n, values = len(res.rows), [r.ccs for r in res.rows]
        if used is not None:                      # compared: the reference value closest to the observed one (relative to its tolerance)
            ref = {"value": _num(used.ccs), "note": f"closest of {n} measured values" if n > 1 else "1 measured value"}
            source = {"value": "Measured", "note": f"{ccs.TYPE_LABEL[used.ref_type]} · {used.row.src_tag}"}
        else:                                     # nothing to compare: the median of the measured values
            ref = {"value": _num(res.ref_ccs), "note": f"median of {n} values ({_num(min(values))}–{_num(max(values))})" if n > 1 else "1 measured value"}
            source = {"value": "Measured", "note": ", ".join(ccs.TYPE_LABEL[t] for t in sorted({r.ccs_type for r in res.rows}, key=ccs.TYPE_ORDER.get))}
    elif res.source == "predicted":
        ref = {"value": _num(p.ccs), "note": f"predicted, {_tol(p.tolerance_pct)} %"}
        source = {"value": "Predicted", "note": f"{p.model}" + ("" if p.in_domain else " · outside its training domain")}
    else:
        ref = {"value": MISSING, "note": res.prediction_note or "none in the index"}
        source = {"value": MISSING}
    # Δ
    delta = {"value": MISSING}
    if res.verdict in (ccs.CONSISTENT, ccs.INCONSISTENT):
        note = f"tolerance {_tol(used.tolerance)} %" + (f"; closest {_obs(used.observed)} Å²" if len(res.observed) > 1 else "")
        delta = {"value": common.signed(used.delta_pct, 2), "note": note}
    # verdict
    verdict = {"value": res.verdict}
    if res.verdict in (ccs.CONSISTENT, ccs.INCONSISTENT):
        verdict["note"] = "vs measured reference" if res.source == "measured" else "vs prediction"
    if res.verdict == ccs.INCONSISTENT:
        verdict["help"] = "Outside the tolerance."
    elif res.verdict == ccs.NOT_COMPARABLE:
        verdict.update(note=res.reason_short, help=res.reason)
    elif res.verdict == ccs.NO_REFERENCE:
        verdict["note"] = res.prediction_note or "none in the index"
    elif res.verdict == ccs.NO_OBSERVED:
        verdict["note"] = "enter an observed CCS to compare"
    return {"Reference CCS": ref, "Source": source, "Observed": observed, "Δ": delta, "Verdict": verdict}


def run_ccs(inputs):
    from mzdesk.orthogonal import ccs

    smiles = common.structure(inputs)["smiles"]
    adduct = common.text(inputs, "adduct") or "[M+H]+"
    platform = common.text(inputs, "platform") or next(iter(ccs.PLATFORMS))
    gas = common.text(inputs, "gas") or ccs.GASES[0]
    tw = common.text(inputs, "tw_calibration") or ccs.TW_CALIBRATIONS[0]
    for value, options, label in ((platform, tuple(ccs.PLATFORMS), "Platform"), (gas, ccs.GASES, "Drift gas"), (tw, ccs.TW_CALIBRATIONS, "TWIMS calibration")):
        if value not in options:
            raise InputError(f"{label} must be one of: {', '.join(options)}.")
    try:
        observed = ccs.parse_observed(inputs.get("observed"))
        res = ccs.analyse(smiles, adduct, observed, platform, gas, tw)
    except ccs.IndexProblem as exc:
        raise EngineError(str(exc)) from None
    except ValueError as exc:                 # QueryError is a ValueError
        raise InputError(str(exc)) from None

    warnings = [f"Structure: {n}." for n in res.notes]
    if platform != "TWIMS" and tw == ccs.TW_CALIBRATIONS[1]:
        warnings.append("The TWIMS calibration choice only applies to TWIMS and was ignored.")
    if not res.installed:
        warnings.append(f"The CCS reference data is not installed. Run {ccs.SETUP_COMMAND}.")
    elif res.prediction is not None and not res.prediction.in_domain:
        warnings.append(f"Outside the model's training domain: the most similar training compound has Tanimoto {res.prediction.max_tanimoto:.2f} "
                        f"(threshold {res.cfg.domain_min_tanimoto:g}).")
    empty = B.table_data([])
    outputs = {
        "summary": {"values": _summary(res)},
        "reference": _reference_table(res) if res.installed else empty,
        "other_ions": _other_ions_table(res) if res.installed else empty,
        "prediction": _prediction_table(res) if res.installed else empty,
        "notes": _notes(res),
    }
    return {"outputs": outputs, "warnings": warnings}


RUNNERS = {"ccs": run_ccs}
