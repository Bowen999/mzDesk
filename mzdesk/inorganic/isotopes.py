"""ICP 同位素选择规则与质量 / 第一电离能匹配的内标排序。"""
from mzdesk.inorganic.core import all_isotopes, average_masses, element_isotopes, parse_isotope
from mzdesk.inorganic.data import INTERNAL_STANDARDS, IONIZATION_EV, MATRIX_ELEMENTS
from mzdesk.inorganic.interferences import find_interferences


def element_report(element, matrix=None):
    matrix = MATRIX_ELEMENTS if matrix is None else matrix
    overlaps = {}
    for isotope in all_isotopes():
        overlaps.setdefault(isotope.mass_no, []).append(isotope)
    rows = []
    for isotope in element_isotopes(element):
        isobars = [other for other in overlaps[isotope.mass_no] if other.element != element and other.abundance > 0]
        poly = [row for row in find_interferences(isotope.label, matrix, 3, (1,))["interferences"]
                if row["type"] == "polyatomic"]
        strongest = poly[0] if poly else None
        burden = sum(row["score"] for row in poly)
        score = isotope.abundance / (1.0 + burden)
        rows.append({
            "isotope": isotope,
            "isobars": isobars,
            "poly_count": len(poly),
            "strongest_poly": strongest,
            "poly_burden": burden,
            "score": score,
        })
    ranked = sorted(rows, key=lambda row: (bool(row["isobars"]), -row["score"], -row["isotope"].abundance, row["isotope"].mass_no))
    recommendations = []
    for rank, row in enumerate(ranked, 1):
        isotope = row["isotope"]
        reason = ("No atomic isobar" if not row["isobars"] else
                  "Atomic isobar: " + ", ".join(i.display for i in row["isobars"]))
        reason += f"; natural abundance {isotope.abundance * 100:.4f}%; polyatomic likelihood sum {row['poly_burden']:.4g}"
        if row["strongest_poly"]:
            top = row["strongest_poly"]
            reason += f"; strongest {top['formula']} (score {top['score']:.4g})"
        else:
            reason += "; no listed polyatomic interference"
        recommendations.append({**row, "rank": rank, "reason": reason})
    return {"isotopes": rows, "recommendations": recommendations}


def internal_standards(element):
    if element not in IONIZATION_EV:
        raise ValueError(f"No first ionization energy is available for {element}.")
    averages = average_masses()
    analyte_mass = averages[element]
    analyte_ie = IONIZATION_EV[element]
    rows = []
    for label in INTERNAL_STANDARDS:
        isotope = parse_isotope(label)
        if isotope.element == element:
            continue
        mass_difference = abs(isotope.mass - analyte_mass)
        ie_difference = abs(IONIZATION_EV[isotope.element] - analyte_ie)
        score = 0.5 * (mass_difference / analyte_mass + ie_difference / analyte_ie)
        rows.append({
            "isotope": isotope,
            "mass_difference": mass_difference,
            "ie_difference": ie_difference,
            "ionization_energy": IONIZATION_EV[isotope.element],
            "score": score,
        })
    rows.sort(key=lambda row: (row["score"], row["isotope"].mass))
    for rank, row in enumerate(rows, 1):
        row["rank"] = rank
    return {"analyte_mass": analyte_mass, "analyte_ie": analyte_ie, "standards": rows}
