"""English display layer for the (Chinese-internal) classifier.

The core `lipid_metabolite_classifier` keeps its Chinese category strings and
reasons untouched. This module maps categories to English and re-derives the
explanation in English from `Result.features`.

The FEATURE lines are formatted from `Features`, but the closing verdict is
rendered from `Result.rule` -- the key of the branch that actually fired. This
module used to re-implement the decision tree to produce that line, which meant
every change to `classify()` had to be mirrored here by hand or the English
"Why" would quietly contradict the verdict shown above it.
"""
from __future__ import annotations

from typing import List, Optional

from lipid_metabolite_classifier import Result, Features, Thresholds


# Canonical Chinese category -> short English label.
CATEGORY_EN = {
    "脂质组": "Lipidome",
    "代谢组": "Metabolome",
    "都有可能": "Both",
    "都不可能": "Neither",
    "无效SMILES": "Invalid",
}

# English label -> one-line subtitle.
CATEGORY_SUBTITLE_EN = {
    "Lipidome": "Likely detected by lipidomics",
    "Metabolome": "Likely detected by metabolomics",
    "Both": "Detectable by both platforms",
    "Neither": "Outside typical LC-MS detection",
    "Invalid": "SMILES could not be parsed",
}

# Core isoprenoid type -> display value.
ISOPRENOID_EN = {
    "sterol": "Sterol",
    "triterpenoid": "Triterpenoid",
    "terpenoid/prenol": "Terpenoid / prenol",
}

# Core lipid class -> display value.
LIPID_CLASS_EN = {
    "glycerophospholipid": "Glycerophospholipid",
    "glycerolipid": "Glycerolipid",
    "sphingolipid": "Sphingolipid",
    "sterol lipid": "Sterol lipid",
    "fatty acyl": "Fatty acyl",
    "prenol lipid": "Prenol lipid",
    "hydrocarbon/fatty alcohol": "Hydrocarbon / fatty alcohol",
}


def category_en(category: str) -> str:
    """Chinese core category -> English label."""
    return CATEGORY_EN.get(category, category)


def category_subtitle_en(english_category: str) -> str:
    """English label -> subtitle string."""
    return CATEGORY_SUBTITLE_EN.get(english_category, "")


def chain_label(thr: Optional[Thresholds] = None) -> str:
    """Label of the lipid-tail feature; follows the tunable chain threshold."""
    return f"≥C{(thr or Thresholds()).acyl_min_carbons} tail"


def lipid_class_en(features: Features) -> str:
    """English name of the matched lipid backbone, or 'none'."""
    if not features.lipid_class:
        return "none"
    return LIPID_CLASS_EN.get(features.lipid_class, features.lipid_class)


def _pretty_heads(features: Features) -> str:
    if not features.polar_head_types:
        return "none"
    return ", ".join(t.replace("_", " ") for t in features.polar_head_types)


def feature_rows_en(features: Features, thr: Optional[Thresholds] = None) -> List[tuple]:
    """Return [(label, value_str), ...] with English labels for the stat cards."""
    f = features
    if not f.valid:
        return [("SMILES", "invalid")]
    if f.has_long_acyl:
        linked = "Yes" if f.is_acyl_linked else "No"
    else:
        linked = "—"  # no long tail to test
    return [
        ("logP", f"{f.logp:.2f}"),
        ("MW (Da)", f"{f.mw:.1f}"),
        ("Lipid class", lipid_class_en(f)),
        ("Longest chain", f"C{f.longest_chain_carbons}"),
        ("Lipid tail", f"C{f.lipid_tail_carbons}"),
        (chain_label(thr), "Yes" if f.has_long_acyl else "No"),
        ("Acyl-linked", linked),
        ("Isoprenoid / sterol", ISOPRENOID_EN.get(f.isoprenoid_type, f.isoprenoid_type or "No")),
        ("Polar head groups", str(f.n_polar_head)),
        ("Glycosyl rings", str(f.n_glycosyl)),
    ]


def _verdict_line(result: Result, thr: Thresholds) -> str:
    """The closing '→ ... → Category.' line for the rule that actually fired."""
    f = result.features
    rule = result.rule

    if rule == "no_carbon":
        return "→ No carbon atom: not an organic small molecule → Neither."
    if rule == "outside_window":
        return "→ Outside the LC-MS detection window → Neither."
    if rule == "steroid_both":
        return (
            f"→ Steroid nucleus with no lipophilic tail and logP ≤ "
            f"{thr.polar_lipid_logp_max:g}: a hormone, bile acid or oxysterol, "
            "which steroid panels and lipidomics both report → Both."
        )
    if rule == "boundary_both":
        return (
            f"→ Medium-chain free fatty acyl (≤C{thr.medium_chain_max} in total), "
            f"logP within [{thr.boundary_logp_low:g}, {thr.boundary_logp_high:g}], "
            f"MW ≤ {thr.boundary_mw_max:g} Da and no polar head: "
            "detectable by both → Both."
        )
    if rule == "lipid_class":
        return f"→ Canonical {lipid_class_en(f).lower()} backbone → Lipidome."
    if rule == "polar_metabolite":
        return (
            f"→ No lipid-class backbone and logP < {thr.logp_polar_max:g}: "
            "polar metabolite → Metabolome."
        )
    if rule == "nonpolar_both":
        return (
            f"→ No lipid-class backbone but logP ≥ {thr.logp_polar_max:g}: "
            "could be either → Both."
        )
    return f"→ {category_en(result.category)}."


def english_reasons(result: Result, thr: Optional[Thresholds] = None) -> List[str]:
    """Re-derive the explanation in English.

    Feature lines are formatted from `Result.features`; the last entry is the
    decision, taken from `Result.rule`, and ends with the English category
    ("→ ... → Lipidome.").
    """
    thr = thr or Thresholds()
    f = result.features
    if not f.valid:
        return ["The SMILES string could not be parsed by RDKit."]

    reasons: List[str] = [
        f"logP = {f.logp:.2f}; molecular weight = {f.mw:.1f} Da."
    ]

    if f.n_fragments > 1:
        reasons.append(
            f"The input has {f.n_fragments} fragments (a salt or mixture); logP and "
            "MW are measured on the largest carbon-containing fragment."
        )

    window = f"[{thr.mw_min:g}, {thr.mw_max:g}] Da"
    if f.mw < thr.mw_min or f.mw > thr.mw_max:
        reasons.append(f"MW is OUTSIDE the typical LC-MS window {window}.")
    else:
        reasons.append(f"MW is within the typical LC-MS window {window}.")

    if f.has_long_acyl:
        link = ("ending in an ester/amide/acid carbonyl"
                if f.is_acyl_linked else "with no carbonyl at either end")
        reasons.append(
            f"Has a ≥C{thr.acyl_min_carbons} mostly-hydrocarbon tail "
            f"(tail C{f.lipid_tail_carbons} of a C{f.longest_chain_carbons} chain), {link}."
        )
    else:
        reasons.append(
            f"No ≥C{thr.acyl_min_carbons} mostly-hydrocarbon tail "
            f"(tail C{f.lipid_tail_carbons} of a C{f.longest_chain_carbons} chain)."
        )

    if f.has_isoprenoid:
        reasons.append(
            "Isoprenoid / sterol scaffold detected "
            f"({ISOPRENOID_EN.get(f.isoprenoid_type, f.isoprenoid_type).lower()}; "
            f"{f.n_isoprene_units} isoprene unit(s))."
        )
    else:
        reasons.append(
            f"No isoprenoid / sterol scaffold ({f.n_isoprene_units} isoprene unit(s))."
        )

    reasons.append(
        f"Polar head groups: {f.n_polar_head} ({_pretty_heads(f)}); "
        f"glycosyl rings: {f.n_glycosyl}."
    )

    if f.lipid_class:
        reasons.append(f"Matches a canonical lipid backbone: {lipid_class_en(f)}.")
    else:
        reasons.append("Matches no canonical lipid backbone.")

    reasons.append(_verdict_line(result, thr))
    return reasons
