"""IsoSpecPy 同位素表的唯一读取层，以及标签、组成和摩尔质量计算。"""
from dataclasses import dataclass
from functools import lru_cache
import math
import re

from mzdesk.inorganic.data import ELEMENT_SYMBOLS
from mzdesk.ms1.isotopes import _ELECTRON

_SUP = str.maketrans("0123456789+-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻")
_FROM_SUP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻", "0123456789+-")


@dataclass(frozen=True)
class Isotope:
    element: str
    mass_no: int
    mass: float
    abundance: float

    @property
    def label(self):
        return f"{self.mass_no}{self.element}"

    @property
    def display(self):
        return isotope_label(self)


@lru_cache(maxsize=1)
def isotope_table():
    from IsoSpecPy import PeriodicTbl
    masses = PeriodicTbl.symbol_to_masses
    if set(masses) - set(ELEMENT_SYMBOLS) - {"D", "E", "Me", "Pn", "T"}:
        raise ValueError("IsoSpecPy contains elements missing from the Inorganic MS element list.")
    table = {}
    for element in ELEMENT_SYMBOLS:
        if element not in masses:
            raise ValueError(f"IsoSpecPy does not provide {element} isotopes.")
        rows = []
        for mass_no, mass, abundance in zip(PeriodicTbl.symbol_to_massNo[element], masses[element], PeriodicTbl.symbol_to_probs[element]):
            if abundance > 0:
                rows.append(Isotope(element, int(mass_no), float(mass), float(abundance)))
        if not rows or not math.isclose(sum(i.abundance for i in rows), 1.0, rel_tol=0, abs_tol=1e-6):
            raise ValueError(f"Invalid natural isotope abundances for {element}.")
        table[element] = tuple(rows)
    return table


def element_isotopes(element):
    if element not in ELEMENT_SYMBOLS:
        raise ValueError(f"Unknown or unsupported element {element!r}.")
    return isotope_table()[element]


def all_isotopes():
    return tuple(isotope for element in ELEMENT_SYMBOLS for isotope in element_isotopes(element))


def isotope_label(value):
    if isinstance(value, Isotope):
        value = value.label
    return str(value).translate(_SUP)


def parse_isotope(value, element=None):
    text = str(value or "").strip().translate(_FROM_SUP)
    match = re.fullmatch(r"(?:(\d+)\s*([A-Z][a-z]?)|([A-Z][a-z]?)\s*(?:-\s*|\s+)?(\d+))", text)
    if not match:
        raise ValueError(f"Invalid isotope {value!r}; use a label such as 75As or As-75.")
    symbol = match.group(2) or match.group(3)
    mass_no = int(match.group(1) or match.group(4))
    if element is not None and symbol != element:
        raise ValueError(f"Isotope {value!r} does not belong to {element}.")
    for isotope in element_isotopes(symbol):
        if isotope.mass_no == mass_no:
            return isotope
    raise ValueError(f"{symbol} has no natural isotope with mass number {mass_no} in IsoSpecPy.")


def natural_fraction(element, mass_no):
    return next(i.abundance for i in element_isotopes(element) if i.mass_no == int(mass_no))


@lru_cache(maxsize=1)
def average_masses():
    from IsoSpecPy import PeriodicTbl
    return {element: float(PeriodicTbl.symbol_to_avg_mass[element]) for element in ELEMENT_SYMBOLS}


def formula_counts(formula):
    from mzdesk.ms1.isotopes import parse_formula
    counts, charge = parse_formula(formula)
    if charge:
        raise ValueError("Particle composition must be a neutral formula without a charge suffix.")
    if any(count <= 0 for count in counts.values()):
        raise ValueError("Particle composition must contain positive atom counts.")
    if any(element not in ELEMENT_SYMBOLS for element in counts):
        raise ValueError("The composition contains an unsupported element.")
    return counts


def formula_average_mass(counts):
    averages = average_masses()
    return sum(averages[element] * count for element, count in counts.items())


def isotopic_molar_mass(element, a, b, fraction_a, fraction_b):
    """未列出的同位素保持天然相对比例，并缩放到剩余原子分数。"""
    rows = element_isotopes(element)
    remaining = 1.0 - fraction_a - fraction_b
    if remaining < -1e-10:
        raise ValueError("The two atom fractions cannot sum to more than 1.")
    remaining = max(0.0, remaining)
    others = [i for i in rows if i.mass_no not in (a.mass_no, b.mass_no)]
    other_natural = sum(i.abundance for i in others)
    if remaining > 1e-10 and not others:
        raise ValueError("The two atom fractions must sum to 1 for this two-isotope element.")
    fractions = {a.mass_no: fraction_a, b.mass_no: fraction_b}
    for isotope in others:
        fractions[isotope.mass_no] = remaining * isotope.abundance / other_natural
    return sum(i.mass * fractions[i.mass_no] for i in rows)


def polyatomic_label(isotopes):
    """Ar 在最前，金属 / 其他原子其次，O、H 最后，便于读出氧化物和羟基离子。"""
    def key(item):
        index, isotope = item
        group = 0 if isotope.element == "Ar" else 2 if isotope.element == "O" else 3 if isotope.element == "H" else 1
        return group, index
    return "".join(isotope_label(i) for _, i in sorted(enumerate(isotopes), key=key)) + "⁺"


def require_finite(value, label):
    if not math.isfinite(value):
        raise ValueError(f"{label} must be finite.")
    return value
