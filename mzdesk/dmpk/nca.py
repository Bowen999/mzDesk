"""Non-compartmental analysis (NCA) of one concentration–time profile.

Pure functions only: input parsing, linear trapezoids, the terminal log-linear fit and the derived
pharmacokinetic parameters. Every convention that is not forced by the data is documented on the page;
unsupported input raises ValueError with a message that can be shown to the user unchanged.
"""
import math

ROUTES = ("IV bolus", "IV infusion", "Oral")
MIN_POINTS = 3


# ---------------------------------------------------------------- parsing
def _split(line):
    """Fields of one pasted line: comma-separated if a comma is present, otherwise whitespace (tabs included)."""
    parts = line.split(",") if "," in line else line.split()
    return [part.strip() for part in parts if part.strip()]


def _float(token):
    try:
        return float(token)
    except (TypeError, ValueError):
        return None


def parse_profile(text):
    """Paste table -> [(time, concentration), ...]; every error names the source line.

    Accepts an optional header line, comma / tab / whitespace separators, blank lines and ``#`` comments
    (whole-line or trailing). Requires exactly two numeric columns, at least three rows, strictly increasing
    times and positive concentrations.
    """
    lines = []
    for line_number, raw in enumerate((text or "").lstrip("\ufeff").splitlines(), start=1):
        raw = raw.split("#", 1)[0].strip()
        if raw:
            lines.append((line_number, raw))
    if not lines:
        raise ValueError("The concentration–time table is empty. Paste one row per line: time and concentration.")

    # An optional header: a first line with no numeric field at all (e.g. "Time (h) Concentration (ng/mL)").
    first_tokens = _split(lines[0][1])
    has_header = bool(first_tokens) and all(_float(token) is None for token in first_tokens)
    rows = lines[1:] if has_header else lines
    if not rows:
        raise ValueError(f"At least {MIN_POINTS} concentration–time rows are required; got 0.")

    points = []
    for line_number, line in rows:
        tokens = _split(line)
        if len(tokens) != 2:
            raise ValueError(f"Row {line_number}: expected 2 values (time and concentration), found {len(tokens)}.")
        values = []
        for token in tokens:
            value = _float(token)
            if value is None:
                raise ValueError(f"Row {line_number}: '{token}' is not a number.")
            if not math.isfinite(value):
                raise ValueError(f"Row {line_number}: time and concentration must be finite numbers.")
            values.append(value)
        time, concentration = values
        if points and time <= points[-1][0]:
            raise ValueError(f"Row {line_number}: time must strictly increase; {time:g} h is not after {points[-1][0]:g} h.")
        if concentration <= 0:
            raise ValueError(f"Row {line_number}: concentration must be > 0 ng/mL (got {concentration:g}).")
        points.append((time, concentration))
    if len(points) < MIN_POINTS:
        raise ValueError(f"At least {MIN_POINTS} concentration–time rows are required; got {len(points)}.")
    return points


# ---------------------------------------------------------------- formulas
def linear_trapezoid(times, values):
    """Sum of (t[i+1] - t[i]) * (v[i] + v[i+1]) / 2 over adjacent points."""
    times = [float(value) for value in times]
    values = [float(value) for value in values]
    if len(times) != len(values):
        raise ValueError("Times and values must have the same length.")
    if len(times) < 2:
        raise ValueError("At least 2 points are required for a trapezoid.")
    return sum((times[i + 1] - times[i]) * (values[i] + values[i + 1]) / 2 for i in range(len(times) - 1))


def fit_lambda_z(times, concentrations, points):
    """Unweighted OLS of ln C vs t over the last ``points`` observations -> λz, intercept, R², fit values."""
    times = [float(value) for value in times]
    concentrations = [float(value) for value in concentrations]
    if len(times) != len(concentrations):
        raise ValueError("Times and concentrations must have the same length.")
    if isinstance(points, bool) or not isinstance(points, int):
        raise ValueError("Points for λz must be an integer.")
    if points < MIN_POINTS:
        raise ValueError(f"Points for λz must be at least {MIN_POINTS}.")
    if points > len(times):
        raise ValueError(f"Points for λz cannot exceed the number of observations ({len(times)}).")
    window_t = times[-points:]
    window_c = concentrations[-points:]
    if any(value <= 0 for value in window_c):
        raise ValueError("All concentrations must be > 0 for the log-linear λz fit.")
    if any(window_t[i + 1] <= window_t[i] for i in range(points - 1)):
        raise ValueError("Times must strictly increase for the λz fit.")

    log_c = [math.log(value) for value in window_c]
    mean_t = sum(window_t) / points
    mean_log_c = sum(log_c) / points
    sxx = sum((value - mean_t) ** 2 for value in window_t)
    sxy = sum((window_t[i] - mean_t) * (log_c[i] - mean_log_c) for i in range(points))
    slope = sxy / sxx
    if slope >= 0:
        raise ValueError("The concentration does not decline over the chosen points.")
    intercept = mean_log_c - slope * mean_t
    predicted_log = [intercept + slope * value for value in window_t]
    ss_res = sum((log_c[i] - predicted_log[i]) ** 2 for i in range(points))
    ss_tot = sum((value - mean_log_c) ** 2 for value in log_c)
    r2 = 1 - ss_res / ss_tot
    return {
        "lambda_z": -slope,
        "slope": slope,
        "intercept": intercept,
        "r2": r2,
        "n_points": points,
        "times": window_t,
        "predicted": [math.exp(value) for value in predicted_log],
    }


def analyse(points, dose_mg, route, lambda_points=3):
    """Full NCA -> dict of numbers (times in h, concentrations in ng/mL, AUC in ng·h/mL, CL in L/h)."""
    points = [(float(time), float(concentration)) for time, concentration in points]
    if len(points) < MIN_POINTS:
        raise ValueError(f"At least {MIN_POINTS} concentration–time points are required; got {len(points)}.")
    previous = None
    for index, (time, concentration) in enumerate(points, start=1):
        if not math.isfinite(time) or not math.isfinite(concentration):
            raise ValueError(f"Point {index}: time and concentration must be finite numbers.")
        if previous is not None and time <= previous:
            raise ValueError(f"Point {index}: time must strictly increase.")
        if concentration <= 0:
            raise ValueError(f"Point {index}: concentration must be > 0 ng/mL.")
        previous = time
    dose_mg = float(dose_mg)
    if not math.isfinite(dose_mg) or dose_mg <= 0:
        raise ValueError("Dose must be a positive number (mg).")
    if route not in ROUTES:
        raise ValueError(f"Unknown route '{route}'; choose IV bolus, IV infusion or Oral.")

    times = [point[0] for point in points]
    concentrations = [point[1] for point in points]
    fit = fit_lambda_z(times, concentrations, lambda_points)

    cmax_index = 0
    for index in range(1, len(concentrations)):
        if concentrations[index] > concentrations[cmax_index]:
            cmax_index = index
    cmax, tmax = concentrations[cmax_index], times[cmax_index]

    auc_t = linear_trapezoid(times, concentrations)
    aumc_t = linear_trapezoid(times, [time * concentration for time, concentration in points])
    lambda_z = fit["lambda_z"]
    c_last, t_last = concentrations[-1], times[-1]
    terminal_auc = c_last / lambda_z
    auc_inf = auc_t + terminal_auc
    aumc_inf = aumc_t + t_last * c_last / lambda_z + c_last / lambda_z ** 2
    t_half = math.log(2) / lambda_z
    # dose mg -> ng: x 1e6; AUC ng·h/mL -> CL mL/h; /1000 -> L/h (so dose x 1e3 / AUC).
    cl = dose_mg * 1e6 / auc_inf / 1000
    vz = cl / lambda_z
    return {
        "cmax": cmax,
        "tmax": tmax,
        "auc_t": auc_t,
        "auc_inf": auc_inf,
        "auc_extrapolated_pct": 100 * terminal_auc / auc_inf,
        "aumc_t": aumc_t,
        "aumc_inf": aumc_inf,
        "mrt": aumc_inf / auc_inf,
        "lambda_z": lambda_z,
        "lambda_intercept": fit["intercept"],
        "lambda_times": fit["times"],
        "lambda_predicted": fit["predicted"],
        "t_half": t_half,
        "r2": fit["r2"],
        "n_points": len(points),
        "n_lambda": fit["n_points"],
        "cl": cl,
        "vz": vz,
        "dose_mg": dose_mg,
        "route": route,
        "t_last": t_last,
        "c_last": c_last,
    }
