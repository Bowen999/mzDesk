"""Adducts 页：结构 + ESI 极性 + 流动相添加剂 → 预测的离子（加合物、源内丢失、二聚体），按预期强度排序；规则见 mzdesk/platform_fit/adduct_rules.py。

Layout: two metric cards (dominant ion, its m/z; the reason is the "?"), a ranked ion table whose Likelihood column is a coloured badge
(each row's "?" gives the reason), and the rule hits / assumptions (assumed pH, lipid class, fallbacks) in a collapsed note.
"""
from mzdesk.api import blocks as B
from mzdesk.api import common
from mzdesk.api.common import InputError, fmt

SIGN = {"+": "positive", "-": "negative"}
TONE = {"high": "good", "medium": "warn", "low": "neutral"}


def run_adducts(inputs):
    from mzdesk.platform_fit import adduct_rules as ar

    smiles = common.structure(inputs)["smiles"]
    mode = inputs.get("ion_mode") or "Positive"         # the API also accepts "Both"
    if mode not in ar.POLARITIES:
        raise InputError(f"Ion mode must be one of {', '.join(ar.POLARITIES)}.")
    add = inputs.get("additive") or "Formic acid"
    if add not in ar.ADDITIVES:
        raise InputError(f"Mobile phase additive must be one of {', '.join(ar.ADDITIVES)}.")
    alkali = inputs.get("alkali")
    alkali = True if alkali is None else bool(alkali)
    try:
        res = ar.predict(smiles, mode, add, alkali)
    except ValueError as exc:
        raise InputError(f"Structure: {exc}") from None

    pols = ("+", "-") if mode == "Both" else ("+",) if mode == "Positive" else ("-",)
    tops = [res.best[p] for p in pols]
    if all(t is None for t in tops):
        base = {"Dominant ion": {"value": "none expected", "note": "", "help": f"nothing forms in {mode.lower()} ESI (see Rule hits)"},
                "m/z": {"value": common.MISSING, "note": ""}}
    elif len(pols) == 1:
        t = tops[0]
        base = {"Dominant ion": {"value": t.adduct, "note": f"{t.likelihood} likelihood", "help": t.reason},
                "m/z": {"value": fmt(t.mz, 4), "note": t.formula or ""}}
    else:
        def part(t, what):
            return "none" if t is None else (t.adduct if what == "ion" else fmt(t.mz, 4))
        base = {"Dominant ion": {"value": " / ".join(part(t, "ion") for t in tops),
                                 "note": " / ".join("none" if t is None else t.likelihood for t in tops) + " likelihood",
                                 "help": [f"{SIGN[p]}: " + ("no ion expected" if t is None else t.reason) for t, p in zip(tops, pols)]},
                "m/z": {"value": " / ".join(part(t, "mz") for t in tops),
                        "note": " / ".join((t.formula or "") if t else "—" for t in tops)}}
    rows = [B.row([i.adduct, i.formula or common.MISSING, fmt(i.mz, 4), B.badge(i.likelihood, TONE[i.likelihood])], help=i.reason) for i in res.ions]
    warnings = list(res.warnings)
    if rows:
        common.note_isotope_labels(warnings, smiles, "The ion formulas and m/z values are")
    return {"outputs": {"base": {"values": base}, "ions": B.table_data(rows, count=len(rows)), "rule_hits": {"items": res.notes}},
            "warnings": warnings}


RUNNERS = {"adducts": run_adducts}
