"""Non-MS detection 页：UV / PDA、原生荧光、CAD / ELSD / RI（规则与阈值见 mzdesk/platform_fit/non_ms_rules.py）。

Layout: one tile per detector (yes / partial / no; the reason is the tile's "?"; the UV tile also carries the chromophore strength and λmax,
the FLD tile the native fluorescence), the structure with the chromophores highlighted (its "?" lists them with their atoms), and the rule
hits in a collapsed note.
"""
from rdkit import Chem

from mzdesk.api import blocks as B
from mzdesk.api import common
from mzdesk.api.common import InputError

STATUS = {"good": "yes", "possible": "partial", "poor": "no"}


def _lines(*parts):
    """Help lines: drop the empty ones."""
    return [p for p in parts if p]


def _atoms(idx):
    """0-based atom indices -> "2–5, 7" (numbered from 1, as in the SMILES and in the figure)."""
    nums, out = sorted(i + 1 for i in idx), []
    for n in nums:
        if out and n == out[-1][1] + 1:
            out[-1][1] = n
        else:
            out.append([n, n])
    return ", ".join(str(a) if a == b else f"{a}–{b}" if b > a + 1 else f"{a}, {b}" for a, b in out)


def run_non_ms(inputs):
    from mzdesk.platform_fit import non_ms_rules as nm

    smiles = common.structure(inputs)["smiles"]
    phase = inputs.get("mobile_phase") or nm.DEFAULT_MOBILE_PHASE
    try:
        a = nm.analyze(smiles, phase)
    except nm.NonMsError as exc:
        raise InputError(str(exc)) from None

    usable = a.strength != nm.NONE

    # detector tiles: UV carries chromophore strength + λmax, FLD the native fluorescence level
    detail = {name: (suit, reason) for name, suit, reason in a.detectors}
    if not usable:
        uv_value = "none"
    elif a.strength == nm.WEAK and a.lam < nm.TH.end_absorption_nm:
        uv_value = "weak · end absorption"
    else:
        uv_value = f"{a.strength} · ~{a.lam:.0f} nm"
    cad_value = "volatile" if a.gc == "yes" else "non-volatile" if a.mw >= nm.TH.mw_cad_min else "small, non-volatile"
    values = {"UV/PDA": uv_value, "FLD": a.fluor, "CAD": cad_value, "ELSD": cad_value, "RI": "isocratic only"}
    helps = {name: _lines(reason) for name, (_, reason) in detail.items()}      # the UV reason names the chromophore class
    if usable:
        helps["UV/PDA"] += [f"λmax by {a.lam_method}: approximate ({nm.LAMBDA_ERROR}); solvent and pH are not modelled"]
    if a.strength != nm.STRONG:
        helps["UV/PDA"] += _lines(nm.deriv_hint(a.reagents))          # none / weak chromophore: OPA / FMOC / DNPH (the FLD reason already carries it)
    helps["FLD"] += ["a structural alert only: quantum yield depends on solvent, pH and quenchers"]
    tiles = [B.tile(name, STATUS[detail[name][0]], values[name], help=helps[name]) for name in ("UV/PDA", "FLD", "CAD", "ELSD", "RI")]

    # figure: the chromophore list explains the atom numbers
    found = [f"{c.name}{f' ×{c.count}' if c.count > 1 else ''}: {c.text}, {c.strength}; atoms {_atoms(c.atoms)}" for c in a.chromophores]
    svg = a.svg or common.depiction_svg(Chem.MolToSmiles(a.mol))
    outputs = {
        "detectors": B.portal_data(tiles),
        "rule_hits": {"items": a.notes},
    }
    if svg:
        outputs["chromophore_structure"] = {"svg": svg, "help": found or "no chromophore found"}
    return {"outputs": outputs, "warnings": a.warnings}


RUNNERS = {"non_ms": run_non_ms}
