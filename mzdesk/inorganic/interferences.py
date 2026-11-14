"""名义质量筛选 + 同位素精确质量：ICP-MS 原子、双电荷与多原子干扰。"""
from functools import lru_cache
from itertools import combinations_with_replacement
import math

from mzdesk.inorganic.core import _ELECTRON, all_isotopes, element_isotopes, isotope_label, parse_isotope, polyatomic_label
from mzdesk.inorganic.data import ELEMENT_SYMBOLS, MATRIX_ELEMENTS


def _charges(values):
    if isinstance(values, (str, int)):
        values = [values]
    charges = set()
    for value in values:
        text = str(value).strip()
        charge = int(text[:-1] if text.endswith("+") else text)
        if charge not in (1, 2):
            raise ValueError("Charge states must be 1+ and/or 2+.")
        charges.add(charge)
    if not charges:
        raise ValueError("Choose at least one charge state.")
    return tuple(sorted(charges))


def _record(target_mz, mass, charge, formula, kind, score):
    mz = (mass - charge * _ELECTRON) / charge
    delta = mz - target_mz
    return {
        "formula": formula,
        "type": kind,
        "mz": mz,
        "delta_mda": delta * 1000.0,
        "required_r": None if delta == 0 else target_mz / abs(delta),
        "score": score,
    }


def _polyatomics(target, matrix, max_atoms, target_mz):
    by_element = {element: element_isotopes(element) for element in matrix}
    limits = {element: (min(i.mass_no for i in rows), max(i.mass_no for i in rows)) for element, rows in by_element.items()}
    found = []
    for n_atoms in range(2, max_atoms + 1):
        for elements in combinations_with_replacement(matrix, n_atoms):
            def visit(position, remaining_mass, selected):
                if position == len(elements):
                    if remaining_mass == 0:
                        score = math.prod(i.abundance for i in selected)
                        found.append(_record(target_mz, sum(i.mass for i in selected), 1,
                                             polyatomic_label(selected), "polyatomic", score))
                    return
                rest = elements[position + 1:]
                min_rest = sum(limits[e][0] for e in rest)
                max_rest = sum(limits[e][1] for e in rest)
                for isotope in by_element[elements[position]]:
                    left = remaining_mass - isotope.mass_no
                    if min_rest <= left <= max_rest:
                        visit(position + 1, left, selected + (isotope,))
            visit(0, target.mass_no, ())
    return found


@lru_cache(maxsize=512)
def _find(target_label, matrix, max_atoms, charges):
    target = parse_isotope(target_label)
    target_mz = target.mass - _ELECTRON
    found = []
    if 1 in charges:
        for isotope in all_isotopes():
            if isotope.element != target.element and isotope.mass_no == target.mass_no:
                found.append(_record(target_mz, isotope.mass, 1, isotope_label(isotope) + "⁺", "isobar", isotope.abundance))
        found.extend(_polyatomics(target, matrix, max_atoms, target_mz))
    if 2 in charges:
        for isotope in all_isotopes():
            if isotope.mass_no == 2 * target.mass_no:
                found.append(_record(target_mz, isotope.mass, 2, isotope_label(isotope) + "²⁺", "doubly charged", isotope.abundance))
    found.sort(key=lambda row: (-row["score"], abs(row["delta_mda"]), row["formula"]))
    return target, tuple(found)


def find_interferences(target, matrix=None, max_atoms=3, charge_states=(1, 2)):
    target = parse_isotope(target)
    matrix = MATRIX_ELEMENTS if matrix is None else matrix
    if isinstance(matrix, str):
        matrix = [matrix]
    unknown = [element for element in matrix if element not in ELEMENT_SYMBOLS]
    if unknown:
        raise ValueError(f"Unknown matrix element(s): {', '.join(unknown)}.")
    matrix = tuple(dict.fromkeys(matrix))
    max_atoms = int(max_atoms)
    if max_atoms not in (2, 3):
        raise ValueError("Maximum atoms per polyatomic ion must be 2 or 3.")
    charges = _charges(charge_states)
    target, rows = _find(target.label, matrix, max_atoms, charges)
    return {"target": target, "interferences": [dict(row) for row in rows]}
