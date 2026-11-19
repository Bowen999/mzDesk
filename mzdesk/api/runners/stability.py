"""Metabolic stability page: ln-linear fit of an in vitro time course, CLint and well-stirred hepatic clearance."""
import math

import plotly.graph_objects as go

from mzdesk.api import blocks, common
from mzdesk.api.common import InputError
from mzdesk.shared import colors as theme
from mzdesk.shared.ms2_plot import base_layout


def _sig(value, digits=3):
    if value is None or not isinstance(value, (int, float)) or not math.isfinite(value):
        return common.MISSING
    return f"{value:.{digits}g}"


def _card(value, note=None):
    return {"value": value, "note": note} if note else value


def _figure(points, result):
    """ln(% remaining) vs time: observed points plus the OLS fit line."""
    times = [point[0] for point in points]
    ln_percent = [math.log(point[1]) for point in points]
    line_x = [min(times), max(times)]
    line_y = [math.log(100.0) + result["intercept"] + result["slope"] * x for x in line_x]
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=line_x, y=line_y, mode="lines", name="OLS fit", hoverinfo="skip",
        line=dict(color=theme.DATA_AMBER, width=2, dash="dash"),
    ))
    fig.add_trace(go.Scatter(
        x=times, y=ln_percent, mode="markers", name="Observed",
        marker=dict(size=9, color=theme.DATA_BLUE, line=dict(color=theme.PAPER, width=1)),
        customdata=[point[1] for point in points],
        hovertemplate="%{x:.4g} min<br>ln(% remaining) %{y:.4f}<br>% remaining %{customdata:.4g}%<extra></extra>",
    ))
    fig.update_layout(**base_layout(height=400, x_range=None, y_range=None, y_title="ln(% remaining)"))
    fig.update_xaxes(title_text="Time (min)")
    fig.add_annotation(
        x=0.02, y=0.98, xref="paper", yref="paper", xanchor="left", yanchor="top", showarrow=False,
        text=f"t½ {_sig(result['t_half'])} min · R² {result['r2']:.4f}",
        font=dict(family=theme.PLOT_FONT, color=theme.MUTED, size=12),
    )
    return common.plotly(fig)


def _parameter_rows(result):
    """Rows of the scaling-factor table: parameter, value, unit and the assumption / source behind it."""
    system = result["system"]
    plasma = system == "Plasma"
    rows = [
        ["k", _sig(result["k"]), "min⁻¹",
         f"Negative slope of the unweighted least-squares fit of ln(% remaining / 100) vs time (n = {result['n']} points)."],
        ["t½", _sig(result["t_half"]), "min", "ln 2 / k."],
        ["R²", f"{result['r2']:.4f}", "—", "Variance explained by the same fit; below 0.9 the estimate is flagged."],
    ]
    if plasma:
        reason = "Not applicable: plasma stability yields t½ only."
        for label in ("CLint (in vitro)", "MPPGL / HPGL", "Liver weight fraction", "CLint scaled", "Qh", "fu",
                      "Predicted hepatic CL", "Extraction ratio"):
            rows.append([label, common.MISSING, "—", reason])
        return rows

    microsomal = system in ("Liver microsomes", "S9")
    per_gram_label = "MPPGL" if microsomal else "HPGL"
    per_gram = result["mppgl"] if microsomal else result["hpgl"]
    per_gram_unit = "mg protein/g liver" if microsomal else "10⁶ cells/g liver"
    if microsomal:
        per_gram_note = (f"Barter et al. 2007, human consensus ({per_gram:g} mg/g liver), applied to all species.")
        if system == "S9":
            per_gram_note += " S9 reuses the microsomal value: a simplification."
    else:
        per_gram_note = (f"Barter et al. 2007, human consensus ({per_gram:g} × 10⁶ cells/g liver; hepatocyte data "
                         "from Wilson et al. 2003), applied to all species.")

    if result["clint"] is None:
        clint_value, clint_note = common.MISSING, "Protein concentration or cell count must be positive."
        clint_unit = "—"
    else:
        clint_value, clint_unit = _sig(result["clint"]), result["clint_unit"]
        clint_note = ("k × 1000 / protein (mg/mL)." if microsomal
                      else "k × 1000 / cell count (10⁶ cells/mL).")

    if result["clint_scaled"] is None:
        scaled_value, scaled_note = common.MISSING, "Needs a positive protein concentration or cell count."
    else:
        scaled_value = _sig(result["clint_scaled"])
        scaled_note = (f"CLint × {per_gram_label} × liver weight ({result['liver_g_per_kg']:g} g/kg) / 1000."
                       + (" S9 reuses microsomal scaling: a simplification." if system == "S9" else ""))

    rows += [
        ["CLint (in vitro)", clint_value, clint_unit, clint_note],
        [per_gram_label, _sig(per_gram), per_gram_unit, per_gram_note],
        ["Liver weight fraction", _sig(result["liver_fraction"]), "g liver/g body weight",
         f"Davies & Morris 1993 ({result['liver_g_per_kg']:g} g liver/kg body weight)."],
        ["CLint scaled", scaled_value, "mL/min/kg", scaled_note],
        ["Qh", _sig(result["qh"]), "mL/min/kg", "Davies & Morris 1993, hepatic blood flow."],
        ["fu", _sig(result["fu"]), "—",
         "No fu input: assumed 1 (no plasma protein binding). CLh is overestimated for bound compounds."],
        ["Predicted hepatic CL", _sig(result["clh"]) if result["clh"] is not None else common.MISSING, "mL/min/kg",
         "Well-stirred: Qh · fu · CLint_scaled / (Qh + fu · CLint_scaled)." if result["clh"] is not None
         else "Needs a scaled CLint."],
        ["Extraction ratio", _sig(result["er"]) if result["er"] is not None else common.MISSING, "—",
         "CLh / Qh; 1 = flow-limited." if result["er"] is not None else "Needs a scaled CLint."],
    ]
    return rows


def _summary(result):
    system = result["system"]
    plasma = system == "Plasma"
    protein = result["protein_or_cells"]
    if plasma:
        clint_note = scaled_note = clh_note = er_note = "plasma stability: t½ only"
    elif result["clint"] is None:
        clint_note = "protein / cell count not positive"
        scaled_note = clh_note = er_note = "needs a positive protein / cell count"
    else:
        amount = (f"{protein:g} × 10⁶ cells/mL" if system == "Hepatocytes" else f"{protein:g} mg/mL")
        clint_note = f"{result['clint_unit']} · {amount}"
        hepatocyte = system == "Hepatocytes"
        per_gram = f"{'HPGL' if hepatocyte else 'MPPGL'} {result['hpgl' if hepatocyte else 'mppgl']:g}"
        scaled_note = f"{result['species']} · {per_gram} × liver {result['liver_g_per_kg']:g} g/kg"
        clh_note, er_note = "well-stirred · fu = 1", "CLh / Qh"
    values = {
        "t½": _card(_sig(result["t_half"]), f"R² {result['r2']:.4f} · {result['n']} points"),
        "CLint": _card(_sig(result["clint"]), clint_note),
        "Scaled CLint": _card(_sig(result["clint_scaled"]), scaled_note),
        "Predicted hepatic CL": _card(_sig(result["clh"]), clh_note),
        "Extraction ratio": _card(_sig(result["er"]), er_note),
    }
    return {"values": values}


def run_stability(inputs):
    from mzdesk.dmpk import stability

    warnings = []
    raw = common.text(inputs, "timecourse")
    if not raw.strip():
        raise InputError("Time course: paste at least three time points (time and % remaining).")
    system = common.text(inputs, "system") or "Liver microsomes"
    species = common.text(inputs, "species") or "Human"
    protein = common.number(inputs, "protein", "Protein (mg/mL) or cells (10⁶/mL)", default=0.5)
    try:
        points = stability.parse_timecourse(raw)
        result = stability.compute(points, system, protein, species)
    except ValueError as exc:
        raise InputError(str(exc)) from None

    if result["r2"] < 0.9:
        warnings.append(f"Poor ln-linear fit (R² = {result['r2']:.4f}, below 0.90); t½ and CLint are uncertain.")
    if result["last_percent"] > 80:
        warnings.append(f"The last point is {result['last_percent']:g}% remaining (above 80%); t½ is poorly "
                        "determined. Extend the incubation.")
    if result["t_half"] > 2 * result["last_time"]:
        warnings.append(f"t½ ({_sig(result['t_half'])} min) is longer than twice the incubation time "
                        f"({2 * result['last_time']:g} min); extend the incubation to determine it.")
    if system != "Plasma" and result["clint"] is None:
        warnings.append("Protein concentration or cell count must be positive to convert k to CLint; "
                        "clearance values are not reported.")

    outputs = {
        "summary": _summary(result),
        "fit_plot": _figure(points, result),
        "parameters": blocks.table_data(_parameter_rows(result)),
        "notes": {"items": [
            "k is the negative slope of an unweighted ordinary-least-squares fit of ln(% remaining / 100) against "
            "time, and t½ = ln 2 / k. The figure plots ln(% remaining); the fit line is the same regression "
            "shifted by ln 100.",
            "CLint = k × 1000 / protein for microsomes and S9 (µL/min/mg protein), and k × 1000 / cell count for "
            "hepatocytes (µL/min/10⁶ cells). Plasma stability has no protein or cell basis, so it reports t½ only.",
            "Scaling: CLint_scaled = CLint × MPPGL (or HPGL) × liver weight fraction. MPPGL = "
            f"{result['mppgl']:g} mg/g and HPGL = {result['hpgl']:g} × 10⁶ cells/g are the human consensus values "
            "of Barter et al. 2007, applied to every species; liver weight fraction and hepatic blood flow are "
            "species-specific (Davies & Morris 1993).",
            f"S9 reuses the microsomal MPPGL ({result['mppgl']:g} mg/g) and scaling. S9 is not microsomal protein "
            "and its recovery differs, so treat an S9 result as an order-of-magnitude estimate.",
            "The predicted hepatic clearance uses the well-stirred model with fu = 1: no plasma protein binding, "
            "no blood-to-plasma partitioning and no active transport. fu = 1 overestimates CLh for bound "
            "compounds; this is the largest assumption on the page.",
            "Linear kinetics, a well-mixed liver and hepatic metabolism only: no renal or extrahepatic clearance, "
            "no microsomal binding and no transporter effects.",
            "A low R², a last point above 80% remaining, or a t½ longer than twice the incubation time make the "
            "estimate unreliable; the warnings say which applies.",
        ]},
    }
    return {"outputs": outputs, "warnings": warnings}


RUNNERS = {"stability": run_stability}
