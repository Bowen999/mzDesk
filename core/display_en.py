"""English display layer for the (Chinese-internal) classifier.

The core `lipid_metabolite_classifier` keeps its Chinese category strings and
reasons untouched. This module maps categories to English and RE-DERIVES the
explanation in English from `Result.features` so the English text always mirrors
the core decision logic exactly. Used by both the Streamlit app and the CLI.
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


def category_en(category: str) -> str:
    """Chinese core category -> English label."""
    return CATEGORY_EN.get(category, category)


def category_subtitle_en(english_category: str) -> str:
    """English label -> subtitle string."""
    return CATEGORY_SUBTITLE_EN.get(english_category, "")


def _pretty_heads(features: Features) -> str:
    if not features.polar_head_types:
        return "none"
    return ", ".join(t.replace("_", " ") for t in features.polar_head_types)


def feature_rows_en(features: Features) -> List[tuple]:
    """Return [(label, value_str), ...] with English labels for the stat cards."""
    f = features
    if not f.valid:
        return [("SMILES", "invalid")]
    return [
        ("logP", f"{f.logp:.2f}"),
        ("MW (Da)", f"{f.mw:.1f}"),
        ("Longest chain", f"C{f.longest_chain_carbons}"),
        ("≥C12 acyl chain", "Yes" if f.has_long_acyl else "No"),
        ("Isoprenoid / sterol", f.isoprenoid_type or "No"),
        ("Polar head groups", str(f.n_polar_head)),
        ("Glycosyl rings", str(f.n_glycosyl)),
    ]


def english_reasons(result: Result, thr: Optional[Thresholds] = None) -> List[str]:
    """Re-derive the explanation in English, mirroring classify() exactly."""
    thr = thr or Thresholds()
    f = result.features
    if not f.valid:
        return ["The SMILES string could not be parsed by RDKit."]

    reasons: List[str] = [
        f"logP = {f.logp:.2f}; molecular weight = {f.mw:.1f} Da."
    ]

    if f.mw < thr.mw_min or f.mw > thr.mw_max:
        reasons.append(
            f"MW is OUTSIDE the typical LC-MS window "
            f"[{thr.mw_min:.0f}, {thr.mw_max:.0f}] Da."
        )
    else:
        reasons.append(
            f"MW is within the typical LC-MS window "
            f"[{thr.mw_min:.0f}, {thr.mw_max:.0f}] Da."
        )

    if f.has_long_acyl:
        link = ("linked to an ester/amide/acid carbonyl"
                if f.is_acyl_linked else "not on a carbonyl")
        reasons.append(
            f"Has a >=C{thr.acyl_min_carbons} aliphatic chain "
            f"(longest C{f.longest_chain_carbons}), {link}."
        )
    else:
        reasons.append(
            f"No >=C{thr.acyl_min_carbons} aliphatic chain "
            f"(longest C{f.longest_chain_carbons})."
        )

    if f.has_isoprenoid:
        reasons.append(f"Isoprenoid / sterol scaffold detected ({f.isoprenoid_type}).")
    else:
        reasons.append("No isoprenoid / sterol scaffold.")

    reasons.append(
        f"Polar head groups: {f.n_polar_head} ({_pretty_heads(f)}); "
        f"glycosyl rings: {f.n_glycosyl}."
    )

    # Final rationale mirrors classify()'s decision tree.
    if not f.has_carbon:
        reasons.append("=> No carbon atom: not an organic small molecule => Neither.")
    elif f.mw < thr.mw_min or f.mw > thr.mw_max:
        reasons.append("=> Outside the LC-MS detection window => Neither.")
    else:
        lipid_scaffold = f.has_long_acyl or f.has_isoprenoid
        amphipathic = f.n_polar_head + f.n_glycosyl
        if lipid_scaffold:
            boundary = (
                f.has_long_acyl
                and not f.has_isoprenoid
                and f.longest_chain_carbons <= thr.medium_chain_max
                and thr.boundary_logp_low <= f.logp <= thr.boundary_logp_high
                and f.mw <= thr.boundary_mw_max
                and amphipathic == 0
            )
            if boundary:
                reasons.append(
                    "=> Medium-length acyl chain, low polarity, no head group: "
                    "detectable by both => Both."
                )
            else:
                reasons.append("=> Clear lipid scaffold => Lipidome.")
        else:
            if f.logp < thr.logp_polar_max:
                reasons.append(
                    f"=> No lipid scaffold and logP < {thr.logp_polar_max}: "
                    "polar metabolite => Metabolome."
                )
            else:
                reasons.append(
                    f"=> No lipid scaffold but logP >= {thr.logp_polar_max}: "
                    "could be either => Both."
                )

    return reasons
