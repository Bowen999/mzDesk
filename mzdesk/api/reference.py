"""One reference spectrum, for a clicked row of the Reference Spectra table (Spectrum · Reference Spectra): the Selected spectrum block
under the table shows it.

POST /api/reference/spectrum answers at once (no job): the record's metadata, its scores against my spectrum when the page had one, the plot
(the spectrum, or the mirror plot against my spectrum), the peak list and the spectrum in every text format. The row says which record
({"library", "row", "accession"}); the block's `query` data (mzdesk/api/runners/ms2.py, run_reference) says what every row shares.

Request: {"library", "row", "accession", "inchikey"?, "compound"?, "peaks"?: [[m/z, intensity]], "precursor_mz"?, "tol_da"?}.
Answer: {"record": {... as reference_lookup.lookup, "url"}, "spectrum": figure data, "formats": {"formats": {...}}, "peak_list": [[m/z, rel. %]]}.
"""
import math

from fastapi import APIRouter, HTTPException

from mzdesk.api.common import clean

router = APIRouter()

MAX_PEAKS = 20000
MAX_TEXT = 200


def _text(body, key, limit=MAX_TEXT):
    value = body.get(key) or ""
    if not isinstance(value, str) or len(value) > limit:
        raise HTTPException(400, f"{key} must be text of at most {limit} characters.")
    return value.strip()


def _number(body, key, default, low, high):
    value = body.get(key)
    if value is None:
        return default
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not low <= value <= high:
        raise HTTPException(400, f"{key} must be a number from {low:g} to {high:g}.")
    return float(value)


def parse_request(body):
    """Request JSON -> checked values, or HTTP 400 with a message for the user."""
    from mzdesk.spec import LIBRARIES

    library = body.get("library")
    if library not in LIBRARIES:
        raise HTTPException(400, f"library must be one of {', '.join(LIBRARIES)}.")
    row = body.get("row")
    if isinstance(row, bool) or not isinstance(row, int) or row < 0:
        raise HTTPException(400, "row must be a record number of 0 or more.")
    raw = body.get("peaks")
    peaks = None
    if raw is not None:
        if not isinstance(raw, list) or len(raw) > MAX_PEAKS:
            raise HTTPException(400, f"peaks must be a list of at most {MAX_PEAKS} [m/z, intensity] pairs.")
        peaks = []
        for p in raw:
            if not (isinstance(p, (list, tuple)) and len(p) == 2
                    and all(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) for v in p) and p[0] > 0 and p[1] >= 0):
                raise HTTPException(400, "Every peak must be [m/z, intensity] with m/z above 0 and intensity 0 or more.")
            peaks.append((float(p[0]), float(p[1])))
        peaks = peaks or None
    return {"library": library, "row": row, "accession": _text(body, "accession"), "inchikey": _text(body, "inchikey", 40),
            "compound": _text(body, "compound"), "peaks": peaks, "precursor_mz": _number(body, "precursor_mz", None, 1, 20000),
            "tol_da": _number(body, "tol_da", 0.02, 0.0001, 5)}


@router.post("/api/reference/spectrum")
def reference_spectrum(body: dict):
    from mzdesk.api.runners.ms2 import reference_view
    from mzdesk.ms2 import library_index as lx
    from mzdesk.ms2 import reference_lookup as rl

    req = parse_request(body)
    try:
        rec = rl.record(req["library"], req["row"], req["accession"], inchikey=req["inchikey"], query=req["peaks"], tolerance=req["tol_da"])
    except lx.LibraryMissing as exc:
        raise HTTPException(503, str(exc)) from None
    except LookupError as exc:
        raise HTTPException(409, str(exc)) from None
    except ValueError as exc:
        raise HTTPException(400, str(exc)) from None
    except ImportError as exc:
        raise HTTPException(503, f"{exc.name or 'A required package'} is not installed on this server, so the spectrum cannot be scored.") from None
    return clean(reference_view(rec, req["peaks"], req["precursor_mz"], req["tol_da"], req["compound"]))
