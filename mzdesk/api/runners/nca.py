"""NCA page: a concentration–time profile, dose and route -> non-compartmental pharmacokinetic parameters."""
import math

import plotly.graph_objects as go

from mzdesk.api import blocks, common
from mzdesk.api.common import InputError
from mzdesk.shared import colors as theme
from mzdesk.shared.ms2_plot import base_layout


def _sig(value, digits=3):
    """3–4 significant digits, thousands separator, MISSING for a missing / non-finite value."""
    if value is None or not math.isfinite(value):
        return common.MISSING
    if value == 0:
        return "0"
    magnitude = math.floor(math.log10(abs(value)))
    factor = 10 ** (magnitude - digits + 1)
    rounded = round(value / factor) * factor
    decimals = max(0, digits - 1 - magnitude)
    if abs(rounded) >= 1e6:
        return f"{rounded:.{digits - 1}e}"
    text = f"{rounded:,.{decimals}f}"
    if "." in text:                     # drop places that carry no information ("4.00" -> "4")
        text = text.rstrip("0").rstrip(".")
    return text


def _card(value, note=None):
    return {"value": value, "note": note} if note else value


def _figure(points, result):
    times = [point[0] for point in points]
    concentrations = [point[1] for point in points]
    span = times[-1] - times[0] or 1.0
    x_range = [min(0.0, times[0]), times[-1] + 0.04 * span]

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=times, y=concentrations, mode="lines+markers", name="Observed",
        line=dict(color=theme.DATA_BLUE, width=2),
        marker=dict(color=theme.DATA_BLUE, size=8, line=dict(color=theme.PAPER, width=1)),
        hovertemplate="%{x:.4g} h<br>%{y:.4g} ng/mL<extra>Observed</extra>",
    ))
    fig.add_trace(go.Scatter(
        x=result["lambda_times"], y=result["lambda_predicted"], mode="lines", name="λz fit",
        line=dict(color=theme.DATA_AMBER, width=2, dash="dash"),
        hovertemplate="%{x:.4g} h<br>%{y:.4g} ng/mL<extra>λz fit</extra>",
    ))
    annotation = dict(
        x=0.985, y=0.96, xref="paper", yref="paper", xanchor="right", yanchor="top", align="right", showarrow=False,
        text=f"λz fit · t½ = {_sig(result['t_half'], 3)} h<br>R² = {result['r2']:.3f} · {result['n_lambda']} points",
        font=dict(family=theme.PLOT_FONT, size=11, color=theme.MUTED),
    )
    fig.update_layout(**base_layout(height=430, x_range=x_range, y_range=None,
                                    y_title="Concentration (ng/mL, log scale)", annotations=[annotation]))
    fig.update_xaxes(title_text="Time (h)", hoverformat=".4g")
    fig.update_yaxes(type="log")
    return common.plotly(fig)


def run_nca(inputs):
    from mzdesk.dmpk import nca

    dose = common.number(inputs, "dose", "Dose", required=True, positive=True)
    route = common.text(inputs, "route") or "IV bolus"
    raw_points = common.number(inputs, "lambda_points", "Points for λz", default=3, minimum=3)
    if raw_points != int(raw_points):
        raise InputError("Points for λz must be an integer.")
    lambda_points = int(raw_points)

    profile = common.text(inputs, "profile")
    try:
        points = nca.parse_profile(profile)
        result = nca.analyse(points, dose, route, lambda_points)
    except ValueError as exc:
        raise InputError(str(exc)) from None

    iv = route in ("IV bolus", "IV infusion")
    if route == "IV infusion":
        mrt_note, mrt_method = "h · includes T_inf/2, not corrected", "AUMC0-∞ / AUC0-∞; includes T_inf/2, not corrected"
    elif route == "Oral":
        mrt_note, mrt_method = "h · MAT + MRT", "AUMC0-∞ / AUC0-∞; MAT + MRT"
    else:
        mrt_note, mrt_method = "h · AUMC0-∞ / AUC0-∞", "AUMC0-∞ / AUC0-∞"
    if iv:
        cl_note, cl_method = "L/h · IV: this is CL", "Dose / AUC0-∞ (IV: CL)"
        vz_note, vz_method = "L · IV: this is Vz", "CL / λz (IV: Vz)"
    else:
        cl_note, cl_method = "L/h · oral: apparent CL/F", "Dose / AUC0-∞ (oral: apparent CL/F)"
        vz_note, vz_method = "L · oral: apparent Vz/F", "CL/F ÷ λz (oral: apparent Vz/F)"

    outputs = {
        "summary": {"values": {
            "Cmax": _sig(result["cmax"], 3),
            "Tmax": _sig(result["tmax"], 3),
            "AUC0-t": _sig(result["auc_t"], 3),
            "AUC0-∞": _sig(result["auc_inf"], 3),
            "t½": _sig(result["t_half"], 3),
            "CL/F": _card(_sig(result["cl"], 3), cl_note),
            "Vz/F": _card(_sig(result["vz"], 3), vz_note),
            "MRT": _card(_sig(result["mrt"], 3), mrt_note),
        }},
        "parameters": blocks.table_data([
            ["Cmax", _sig(result["cmax"], 4), "ng/mL", "Largest observed concentration (first on ties)"],
            ["Tmax", _sig(result["tmax"], 4), "h", "Time of Cmax"],
            ["AUC0-t", _sig(result["auc_t"], 4), "ng·h/mL", "Linear trapezoid between adjacent observations"],
            ["AUC0-∞", _sig(result["auc_inf"], 4), "ng·h/mL", "AUC0-t + C_last / λz"],
            ["AUC extrapolated (%)", _sig(result["auc_extrapolated_pct"], 3), "%", "(C_last / λz) / AUC0-∞ × 100"],
            ["t½", _sig(result["t_half"], 4), "h", "ln2 / λz"],
            ["λz", _sig(result["lambda_z"], 4), "1/h", f"Unweighted OLS of ln C vs t (last {result['n_lambda']} points)"],
            ["CL/F", _sig(result["cl"], 4), "L/h", cl_method],
            ["Vz/F", _sig(result["vz"], 4), "L", vz_method],
            ["MRT", _sig(result["mrt"], 4), "h", mrt_method],
            ["Points used", str(result["n_lambda"]), "points", "Observations in the λz window"],
            ["R²", f"{result['r2']:.4f}", common.MISSING, "Goodness of fit of the λz regression (ln C vs t)"],
        ], count=12),
        "profile_plot": _figure(points, result),
        "notes": {"items": [
            "AUC0-t and AUMC0-t use the linear trapezoidal rule between adjacent observations, so AUC0-t spans the observed interval only; no logarithmic interpolation and no C0 back-extrapolation.",
            f"λz is the unweighted OLS slope of ln C vs t over the last {result['n_lambda']} observations; t½ = ln2/λz. The window is chosen by the user, not by an algorithm.",
            "AUC0-∞ = AUC0-t + C_last/λz and AUMC0-∞ = AUMC0-t + t_last·C_last/λz + C_last/λz², assuming a mono-exponential terminal phase.",
            "AUC extrapolated = (C_last/λz) / AUC0-∞; values above 20% are flagged.",
            "Concentrations must be > 0. BLQ values are not modelled: replace or exclude them before pasting, and note that a pre-dose 0 ng/mL sample cannot be entered.",
            "No steady-state handling: one dose is assumed, with time measured from administration; times are used as given.",
            "NCA assumes linear kinetics. CL/F and Vz/F are apparent values; they equal CL and Vz only for IV dosing.",
            "MRT = AUMC0-∞/AUC0-∞. For IV infusion it includes half the unknown infusion time and is not corrected; for oral dosing it is MAT + MRT.",
            "Units: time h, concentration ng/mL, dose mg. CL/F = dose (mg) × 10³ / AUC0-∞ (ng·h/mL) in L/h; Vz/F = CL/F ÷ λz in L.",
            "Cmax is the largest observed concentration (first on ties); Tmax is its sampling time. No smoothing or interpolation is applied.",
        ]},
    }
    warnings = []
    if result["auc_extrapolated_pct"] > 20:
        warnings.append(f"AUC0-∞ extrapolates {_sig(result['auc_extrapolated_pct'], 3)}% of the total exposure; "
                        "the terminal phase may be too short for a reliable AUC0-∞.")
    if result["r2"] < 0.9:
        warnings.append(f"The λz regression R² is {result['r2']:.3f} (< 0.90); the terminal phase is not well described "
                        "by a single exponential.")
    if result["n_lambda"] <= 3:
        warnings.append(f"λz uses only {result['n_lambda']} points; a {result['n_lambda']}-point terminal window is sensitive "
                        "to a single observation.")
    return {"outputs": outputs, "warnings": warnings}


RUNNERS = {"nca": run_nca}
