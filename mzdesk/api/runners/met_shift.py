"""Mass-shift search: parent ion + measured features + fixed biotransformation shifts."""
import plotly.graph_objects as go

from mzdesk.api import blocks, common
from mzdesk.api.common import InputError
from mzdesk.shared import colors as theme
from mzdesk.shared.ms2_plot import base_layout


def _selected_shifts(inputs):
    raw = inputs.get("shifts")
    if isinstance(raw, str):
        selected = [raw] if raw.strip() else []
    elif isinstance(raw, (list, tuple)):
        selected = [str(item) for item in raw if str(item).strip()]
    else:
        selected = []
    if not selected:
        raise InputError("Choose at least one mass shift.")
    return selected


def _intensity(value):
    if value is None:
        return common.MISSING
    return f"{value:,.0f}" if abs(value) >= 1e4 else f"{value:.4g}"


def _features(inputs, warnings):
    value = inputs.get("features_file") or {}
    file_text = str(value.get("text") or "")
    grid = common.text(inputs, "features")
    if file_text.strip() and grid:
        warnings.append("Both a feature file and pasted features were given; the file is used.")
    text = file_text if file_text.strip() else grid
    if not text.strip():
        return []
    from mzdesk.calculator import mass_shift

    try:
        features, parse_warnings = mass_shift.parse_features(text)
    except ValueError as exc:
        raise InputError(str(exc)) from None
    warnings.extend(parse_warnings)
    return features


def _figure(features, rows):
    if not features:
        fig = go.Figure()
        fig.add_annotation(text="No feature table loaded", showarrow=False,
                           font=dict(family=theme.PLOT_FONT, color=theme.MUTED, size=13))
        fig.update_layout(**base_layout(height=360, x_range=[0, 1], y_range=[0, 1], y_title="Intensity"))
        return common.plotly(fig)
    has_rt = all(feature["rt"] is not None for feature in features)
    index_of = {feature["line"]: index for index, feature in enumerate(features)}
    xs = [feature["rt"] if has_rt else index for index, feature in enumerate(features)]
    ys = [feature["intensity"] if feature["intensity"] is not None else 1.0 for feature in features]
    matched_lines = {row["feature"]["line"] for row in rows if row["feature"] is not None}
    fig = go.Figure()
    unmatched = [(xs[index], ys[index], feature["mz"]) for index, feature in enumerate(features)
                 if feature["line"] not in matched_lines]
    if unmatched:
        fig.add_trace(go.Scatter(
            x=[item[0] for item in unmatched], y=[item[1] for item in unmatched],
            mode="markers", name="Unmatched", showlegend=False,
            marker=dict(size=7, color=theme.DATA_GREY),
            customdata=[item[2] for item in unmatched],
            hovertemplate="Unmatched<br>m/z %{customdata:.4f}<extra></extra>",
        ))
    for order, label in enumerate(dict.fromkeys(row["label"] for row in rows)):
        matched = [row for row in rows if row["label"] == label and row["feature"] is not None]
        if not matched:
            continue
        fig.add_trace(go.Scatter(
            x=[xs[index_of[row["feature"]["line"]]] for row in matched],
            y=[row["feature"]["intensity"] if row["feature"]["intensity"] is not None else 1.0 for row in matched],
            mode="markers", name=label, showlegend=False,
            marker=dict(size=11, color=theme.SERIES[order % len(theme.SERIES)], line=dict(color=theme.PAPER, width=1)),
            customdata=[[row["feature"]["mz"], row["ppm"], row["expected"]] for row in matched],
            hovertemplate=f"{label}<br>found %{{customdata[0]:.4f}}<br>expected %{{customdata[2]:.4f}}<br>%{{customdata[1]:+.2f}} ppm<extra></extra>",
        ))
    fig.update_layout(**base_layout(height=420, x_range=None, y_range=None, y_title="Intensity"))
    fig.update_xaxes(title_text="RT (min)" if has_rt else "Feature index")
    return common.plotly(fig)


def run_met_shift(inputs):
    from mzdesk.calculator import calculations, mass_shift

    warnings = []
    struct = common.structure(inputs, "structure", "Parent compound or formula", allow_formula=True)
    identity = mass_shift.parent_identity(struct)
    precursor_mz, positive, adduct = common.precursor(inputs, label="Measured precursor m/z")
    adduct = adduct or ("[M+H]+" if positive else "[M-H]-")
    try:
        charge = calculations.adduct_charge(adduct)
        theoretical_mz = calculations.adduct_mz(identity["formula"], adduct)
    except ValueError as exc:
        raise InputError(str(exc)) from None
    base_mz = precursor_mz or theoretical_mz
    if precursor_mz is not None:
        error = calculations.ppm_error(precursor_mz, theoretical_mz)
        if abs(error) > 10:
            warnings.append(f"The measured precursor is {common.signed(error, 2)} ppm from the theoretical {adduct} ion ({theoretical_mz:.5f}).")
    selected = _selected_shifts(inputs)
    tol_ppm = common.number(inputs, "tol_ppm", "Match tolerance", default=10.0, positive=True)
    if tol_ppm > 1000:
        warnings.append("The match tolerance is very wide (> 1000 ppm); review the matches manually.")
    features = _features(inputs, warnings)
    try:
        result = mass_shift.search(features, base_mz, charge, selected, tol_ppm)
    except ValueError as exc:
        raise InputError(str(exc)) from None

    rows = []
    for row in result["rows"]:
        feature = row["feature"]
        rows.append([
            row["label"], common.signed(row["delta"], 4), common.fmt(row["expected"], 5),
            common.fmt(feature["mz"], 5) if feature else common.MISSING,
            common.signed(row["ppm"], 2) if feature else common.MISSING,
            common.fmt(feature["rt"], 2) if feature and feature["rt"] is not None else common.MISSING,
            _intensity(feature["intensity"]) if feature else common.MISSING,
            row["reaction"],
        ])
    outputs = {
        "hits": blocks.table_data(rows, count=len(rows)),
        "features_plot": _figure(features, result["rows"]),
        "notes": {"items": [
            f"Parent: {identity['formula']}, monoisotopic mass {identity['mass']:.5f} Da; {adduct} ion m/z {theoretical_mz:.5f}.",
            "Expected m/z = base m/z + Δmass / |z|. The measured precursor is used as the base when supplied; otherwise the theoretical ion m/z is used.",
            "Feature table: one row per line as m/z [RT] [intensity], with an optional header. Comma, tab, semicolon and whitespace are accepted.",
            "A shift may match more than one feature; unmatched shifts remain in the table as expected-only rows.",
            "Mass shifts are fixed biotransformation presets; they are hypotheses, not compound identifications.",
        ]},
    }
    if features:
        n_features = len(features)
        verdict = (f"Parent <b>{adduct} m/z {base_mz:.4f}</b> ({identity['formula']}). "
                   f"<b>{result['matched_shifts']}/{len(selected)}</b> shifts matched; "
                   f"<b>{len(result['matched_features'])}/{n_features}</b> features assigned.")
    else:
        verdict = f"Parent <b>{adduct} m/z {base_mz:.4f}</b> ({identity['formula']}); load a feature table to search."
    return {"outputs": outputs, "warnings": warnings, "verdict": verdict}


RUNNERS = {"met_shift": run_met_shift}
