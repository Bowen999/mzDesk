from __future__ import annotations

from typing import List, Optional

from mzdesk.platform_fit.lipid_metabolite_classifier import Result, Features, Thresholds

CATEGORY_EN = {
    "脂质组": "Lipidome",
    "代谢组": "Metabolome",
    "都有可能": "Both",
    "都不可能": "Neither",
    "无效SMILES": "Invalid",
}

CATEGORY_SUBTITLE_EN = {
    "Lipidome": "Likely detected by lipidomics",
    "Metabolome": "Likely detected by metabolomics",
    "Both": "Detectable by both platforms",
    "Neither": "Outside typical LC-MS detection",
    "Invalid": "SMILES could not be parsed",
}

ISOPRENOID_EN = {
    "sterol": "Sterol",
    "triterpenoid": "Triterpenoid",
    "terpenoid/prenol": "Terpenoid / prenol",
}

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
    return CATEGORY_EN.get(category, category)


def category_subtitle_en(english_category: str) -> str:
    return CATEGORY_SUBTITLE_EN.get(english_category, "")


def chain_label(thr: Optional[Thresholds] = None) -> str:
    return f"≥C{(thr or Thresholds()).acyl_min_carbons} tail"


def lipid_class_en(features: Features) -> str:
    if not features.lipid_class:
        return "none"
    return LIPID_CLASS_EN.get(features.lipid_class, features.lipid_class)


def _pretty_heads(features: Features) -> str:
    if not features.polar_head_types:
        return "none"
    return ", ".join(head.replace("_", " ") for head in features.polar_head_types)


def feature_rows_en(features: Features, thr: Optional[Thresholds] = None) -> List[tuple]:
    if not features.valid:
        return [("SMILES", "invalid")]
    if features.has_long_acyl:
        linked = "Yes" if features.is_acyl_linked else "No"
    else:
        linked = "—"
    return [
        ("logP", f"{features.logp:.2f}"),
        ("MW (Da)", f"{features.mw:.1f}"),
        ("Lipid class", lipid_class_en(features)),
        ("Longest chain", f"C{features.longest_chain_carbons}"),
        ("Lipid tail", f"C{features.lipid_tail_carbons}"),
        (chain_label(thr), "Yes" if features.has_long_acyl else "No"),
        ("Acyl-linked", linked),
        ("Isoprenoid / sterol", ISOPRENOID_EN.get(
            features.isoprenoid_type, features.isoprenoid_type or "No")),
        ("Polar head groups", str(features.n_polar_head)),
        ("Glycosyl rings", str(features.n_glycosyl)),
    ]


def _verdict_line(result: Result, thr: Thresholds) -> str:
    features = result.features
    rule = result.rule
    if rule == "no_carbon":
        return "No carbon → Neither."
    if rule == "outside_window":
        return "MW outside the LC-MS window → Neither."
    if rule == "steroid_both":
        return f"Tail-free steroid with logP ≤ {thr.polar_lipid_logp_max:g} → Both."
    if rule == "boundary_both":
        return (f"Free fatty acyl ≤C{thr.medium_chain_max}; "
                f"logP [{thr.boundary_logp_low:g}, {thr.boundary_logp_high:g}]; "
                f"MW ≤ {thr.boundary_mw_max:g} Da; no polar head → Both.")
    if rule == "lipid_class":
        return f"{lipid_class_en(features)} backbone → Lipidome."
    if rule == "polar_metabolite":
        return f"No lipid backbone; logP < {thr.logp_polar_max:g} → Metabolome."
    if rule == "nonpolar_both":
        return f"No lipid backbone; logP ≥ {thr.logp_polar_max:g} → Both."
    return f"{category_en(result.category)}."


def english_reasons(result: Result, thr: Optional[Thresholds] = None) -> List[str]:
    thr = thr or Thresholds()
    features = result.features
    if not features.valid:
        return ["Invalid SMILES."]

    reasons = [f"logP {features.logp:.2f}; MW {features.mw:.1f} Da."]
    if features.n_fragments > 1:
        reasons.append(
            f"{features.n_fragments} fragments; logP and MW use the largest carbon-containing fragment."
        )

    window = f"[{thr.mw_min:g}, {thr.mw_max:g}] Da"
    if features.mw < thr.mw_min or features.mw > thr.mw_max:
        reasons.append(f"MW outside {window}.")
    else:
        reasons.append(f"MW within {window}.")

    if features.has_long_acyl:
        link = "acyl-linked" if features.is_acyl_linked else "no terminal carbonyl"
        reasons.append(
            f"≥C{thr.acyl_min_carbons} tail: C{features.lipid_tail_carbons} of "
            f"C{features.longest_chain_carbons}; {link}."
        )
    else:
        reasons.append(
            f"No ≥C{thr.acyl_min_carbons} tail "
            f"(C{features.lipid_tail_carbons} of C{features.longest_chain_carbons})."
        )

    if features.has_isoprenoid:
        isoprenoid = ISOPRENOID_EN.get(features.isoprenoid_type, features.isoprenoid_type)
        reasons.append(
            f"Isoprenoid/sterol: {isoprenoid.lower()}; {features.n_isoprene_units} isoprene units."
        )
    else:
        reasons.append(f"No isoprenoid/sterol scaffold; {features.n_isoprene_units} isoprene units.")

    reasons.append(
        f"Polar heads {features.n_polar_head} ({_pretty_heads(features)}); "
        f"glycosyl rings {features.n_glycosyl}."
    )
    if features.lipid_class:
        reasons.append(f"Lipid backbone: {lipid_class_en(features)}.")
    else:
        reasons.append("No canonical lipid backbone.")

    reasons.append(_verdict_line(result, thr))
    return reasons
