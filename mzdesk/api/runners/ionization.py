"""Ionization page: ESI polarity and response (LC-MS), EI usability (GC-MS). Rules and thresholds: mzdesk/platform_fit/ionization_rules.py.

Layout: five coloured chips (Polarity, ESI+, ESI-, EI, M+•; the deciding reason is each chip's "?"), the structure with the ionizable
sites highlighted, and the rule hits in a collapsed note. Nothing explanatory is visible.
"""
from mzdesk.api import blocks as B
from mzdesk.api import common
from mzdesk.api.common import InputError

ESI_LEVEL = {"strong": "good", "moderate": "possible", "weak": "poor"}
EI_LEVEL = {"yes": "good", "after derivatization": "possible", "no": "poor"}
EI_VALUE = {"yes": "usable", "after derivatization": "after derivatization", "no": "not usable"}


def _lines(*parts):
    """Help lines: drop the empty ones."""
    return [p for p in parts if p]


def run_ionization(inputs):
    from mzdesk.platform_fit import ionization_rules as ir
    from mzdesk.platform_fit import ionization_viz as viz

    smiles = common.structure(inputs)["smiles"]
    ph = common.number(inputs, "ph", "Mobile phase pH", default=2.7, minimum=ir.PH_MIN, maximum=ir.PH_MAX)
    try:
        a = ir.assess(smiles, ph)
    except ValueError as exc:
        raise InputError(f"Structure: {exc}") from None

    ion_levels = {ir.VISIBLE: "good", ir.WEAK_ION: "possible", ir.ABSENT: "poor"}      # M+• visible / weak / not expected
    ion_value = a.ion.level or a.ion.value
    if a.ion.value in (ir.NEIMS_MISSING, "not available"):
        ion_value += " (rule of thumb)"                     # not a NEIMS prediction: say so on the chip
    if a.gc["gc"] != "yes":
        ion_value += " (underivatized)" if a.gc["gc"] == "after derivatization" else " (hypothetical)"      # not the molecule the detector sees
    why = {ir.POS: "ESI+ responds better than ESI−", ir.NEG: "ESI− responds better than ESI+",
           ir.BOTH: "both modes respond at least moderately", ir.POOR: "both modes respond weakly"}[a.polarity]
    solution = a.polarity_note.partition("; ")[2]          # for "both": which mode the chosen pH favours, or "try both modes"
    chips = B.chips_data([
        B.chip("Polarity", "poor" if a.polarity == ir.POOR else "good", a.polarity, _lines(why, solution)),
        B.chip("ESI+", ESI_LEVEL[a.pos.level], a.pos.level, a.pos.reason),
        B.chip(f"ESI{ir.MINUS}", ESI_LEVEL[a.neg.level], a.neg.level, a.neg.reason),
        B.chip("EI", EI_LEVEL[a.gc["gc"]], EI_VALUE[a.gc["gc"]], [r.removeprefix("Rule: ") for r in a.gc["reasons"]]),
        B.chip("M+•", ion_levels.get(a.ion.level, "poor"), ion_value, _lines(a.ion.note, ir.ion_caveat(a.gc))),
    ])

    outputs = {
        "chips": chips,
        "sites_structure": {"svg": viz.sites_svg(a.profile)},
        "rule_hits": {"items": a.hits},
    }
    return {"outputs": outputs, "warnings": a.warnings}


RUNNERS = {"ionization": run_ionization}
