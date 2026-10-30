"""Spectrum figures for the workflow; these never participate in annotation."""
import logging

from mzdesk.api.common import plotly

log = logging.getLogger("mzdesk.annotation")


def spectrum_figures(peaks, precursor_mz, tolerance, reference=None):
    """Return the experimental spectrum and, if available, the best judged library hit."""
    from mzdesk.shared.ms2_plot import spectrum_figure

    result = {"query": None, "library": None}
    try:
        figure, cards = spectrum_figure(peaks, precursor_mz=precursor_mz, y_title="Intensity")
        result["query"] = plotly(figure, peaks=cards)
    except Exception:  # noqa: BLE001 - display failures must not stop a workflow
        log.exception("query spectrum figure failed")
    if reference and reference.get("peaks"):
        try:
            from mzdesk.ms2.mirror_plot import build_mirror_figure

            figure = build_mirror_figure(peaks, reference["peaks"], precursor_mz, reference.get("precursor_mz"),
                                         tolerance=tolerance, labels=("Experimental", "Library"))
            result["library"] = {"figure": plotly(figure), "library": reference["library"],
                                 "accession": reference["accession"], "name": reference["name"]}
        except Exception:  # noqa: BLE001
            log.exception("library mirror figure failed")
    return result
