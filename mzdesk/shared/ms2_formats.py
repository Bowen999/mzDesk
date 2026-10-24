"""MS2 峰表写出器：同一张谱写成 MS-DIAL、GNPS、MassBank、HMDB、MSP、MGF 等常见格式。每种输出都能被 parse_spectrum 读回。"""
import json

FORMATS = ["Peak list", "MS-DIAL / MoNA", "GNPS peaks_json", "Two arrays", "MSP", "MGF", "MassBank",
           "SIRIUS .ms", "HMDB XML", "CSV", "JSON"]


def _mz(x):
    """m/z to 5 dp (the parser's duplicate resolution), trailing zeros stripped."""
    return f"{x:.5f}".rstrip("0").rstrip(".")


def _int(y):
    if float(y).is_integer():
        return f"{y:.0f}"
    s = f"{y:.6g}"
    return f"{y:.2f}" if "e" in s else s


def _line(name):
    return " ".join(str(name).split()) or "Unknown"


def to_format(peaks, fmt, *, precursor_mz=None, name=None, ion_mode=None):
    """Write ``peaks`` [(mz, intensity)] as text in format ``fmt`` (one of FORMATS), in the given order.

    ``precursor_mz``, ``name`` and ``ion_mode`` ("Positive" / "Negative") go into the headers of the
    record formats (MSP, MGF, MassBank, SIRIUS) when given.
    """
    p = [(_mz(a), _int(b)) for a, b in peaks]
    title = _line(name or "Unknown")
    pre = None if precursor_mz is None else _mz(precursor_mz)
    pos = None if not ion_mode else str(ion_mode).lower().startswith("pos")
    lines = [f"{a}\t{b}" for a, b in p]
    if fmt == "Peak list":
        return "\n".join(lines)
    if fmt == "MS-DIAL / MoNA":
        return " ".join(f"{a}:{b}" for a, b in p)
    if fmt == "GNPS peaks_json":
        return "[" + ", ".join(f"[{a}, {b}]" for a, b in p) + "]"
    if fmt == "Two arrays":
        return f"mz = [{', '.join(a for a, _ in p)}]\nintensities = [{', '.join(b for _, b in p)}]"
    if fmt == "MSP":
        head = [f"NAME: {title}"] + ([f"PRECURSORMZ: {pre}"] if pre else []) \
            + ([f"IONMODE: {'Positive' if pos else 'Negative'}"] if pos is not None else []) + [f"Num Peaks: {len(p)}"]
        return "\n".join(head + lines)
    if fmt == "MGF":
        head = ["BEGIN IONS", f"TITLE={title}"] + ([f"PEPMASS={pre}"] if pre else []) \
            + ([f"CHARGE=1{'+' if pos else '-'}"] if pos is not None else [])
        return "\n".join(head + lines + ["END IONS"])
    if fmt == "MassBank":
        top = max((b for _, b in peaks), default=1.0)
        head = [f"CH$NAME: {title}"] + ([f"MS$FOCUSED_ION: PRECURSOR_M/Z {pre}"] if pre else []) \
            + ([f"AC$MASS_SPECTROMETER: ION_MODE {'POSITIVE' if pos else 'NEGATIVE'}"] if pos is not None else []) \
            + [f"PK$NUM_PEAK: {len(p)}", "PK$PEAK: m/z int. rel.int."]
        rows = [f"  {a} {b} {max(1, round(y / top * 999))}" for (a, b), (_, y) in zip(p, peaks)]
        return "\n".join(head + rows + ["//"])
    if fmt == "SIRIUS .ms":
        head = [f">compound {title}"] + ([f">parentmass {pre}"] if pre else []) + [">ms2peaks"]
        return "\n".join(head + lines)
    if fmt == "HMDB XML":
        body = [f"  <ms-ms-peak>\n    <mass-charge>{a}</mass-charge>\n    <intensity>{b}</intensity>\n  </ms-ms-peak>" for a, b in p]
        return "\n".join(["<ms-ms-peaks>"] + body + ["</ms-ms-peaks>"])
    if fmt == "CSV":
        return "\n".join(["mz,intensity"] + [f"{a},{b}" for a, b in p])
    if fmt == "JSON":
        rows = [f'  {{"mz": {a}, "intensity": {b}}}' for a, b in p]
        return "[\n" + ",\n".join(rows) + "\n]" if rows else "[]"
    raise ValueError(f"unknown format {fmt!r}; choose one of {', '.join(FORMATS)}")


def all_formats(peaks, **kw):
    """{label: text} for every format in FORMATS (keyword arguments as in ``to_format``)."""
    return {f: to_format(peaks, f, **kw) for f in FORMATS}
