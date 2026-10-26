"""Builders for the runtime data of the generic output blocks (spec helpers and contracts: mzdesk/spec.py).

Every builder returns plain JSON-able dicts / lists and leaves out keys that are None or empty, so the frontend sees exactly the
documented shape. The closed vocabularies (tile status, chip level) are checked here, so a typo fails in the runner or its test
instead of silently showing as a grey tile in the browser. `help` is a string or a list of strings everywhere.

    portal_data([tile("Solid", "yes", "GC-ready", chips=["TMS"], help="...")])
    chips_data([chip("ESI+", "good", "pKa 9.1")])
    table_data([row([link("JP18", url), badge("JP", "jp")], help="why")], count=1)
    grouptable_data([group("Caffeine", [[...]], meta="0.66 · 3 methods")])
    reaction_data([reaction_option("CIL · dansyl", [step("R–NH₂", "Dns-Cl", "R–NH–Dns", 2)], formula="C₁₂H₁₂N₂O₂S", mass=248.0619)], svg=svg)
"""

STATUSES = ("yes", "partial", "no")
LEVELS = ("good", "possible", "poor", "strong", "moderate", "weak", "yes", "no")


def _keep(d):
    """Drop None / empty-string / empty-list entries (0 and False are kept). The main list of each block (rows, tiles, ...) is always written out."""
    return {k: v for k, v in d.items() if v is not None and v != "" and v != []}


# ---- table / grouptable cells and rows ----
def link(text, url, help=None):
    return _keep({"text": text, "url": url, "help": help})


def badge(text, tone="neutral", url=None, help=None):
    return _keep({"badge": text, "tone": tone, "url": url, "help": help})


def hinted(text, help):
    """Text with an inline "?"."""
    return _keep({"text": text, "help": help})


def row(cells, help=None):
    """A row with its own "?" (shown with spec row_help=True). Without help it is just the list of cells."""
    return {"cells": list(cells), "help": help} if help else list(cells)


def table_data(rows, count=None, help=None):
    return {"rows": list(rows), **_keep({"count": count, "help": help})}


def group(title, rows, meta=None, help=None):
    return {"title": title, "rows": list(rows), **_keep({"meta": meta, "help": help})}


def grouptable_data(groups, count=None, help=None):
    return {"groups": list(groups), **_keep({"count": count, "help": help})}


# ---- portal ----
def tile(label, status, value=None, chips=None, help=None):
    if status not in STATUSES:
        raise ValueError(f"tile status must be one of {STATUSES}, got {status!r}")
    return _keep({"label": label, "status": status, "value": value, "chips": list(chips or []), "help": help})


def portal_data(tiles, help=None):
    return {"tiles": list(tiles), **_keep({"help": help})}


# ---- chips ----
def chip(label, level, value=None, help=None):
    if level not in LEVELS:
        raise ValueError(f"chip level must be one of {LEVELS}, got {level!r}")
    return _keep({"label": label, "level": level, "value": value, "help": help})


def chips_data(chips, help=None):
    return {"chips": list(chips), **_keep({"help": help})}


# ---- reaction ----
def step(from_, reagent, to, count=1):
    """One functional-group step: from —reagent→ to, ×count (count 1 is not drawn)."""
    return _keep({"from": from_, "reagent": reagent, "to": to, "count": count if count and count > 1 else None})


def reaction_option(label, steps, formula=None, mass=None, heavy=None, help=None, svg=None):
    steps = list(steps)
    if not steps:
        raise ValueError("a reaction option needs at least one step")
    return _keep({"label": label, "steps": steps, "formula": formula, "mass": mass, "heavy": heavy, "help": help, "svg": svg})


def reaction_data(options, svg=None, help=None):
    options = list(options)
    if not options:
        raise ValueError("a reaction block needs at least one option")
    return _keep({"svg": svg, "options": options, "help": help})
