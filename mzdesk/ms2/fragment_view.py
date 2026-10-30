"""Fragment view: which part of the structure each explained peak comes from (the pop-up of Alternative explanations and ESI MS/MS Annotation).

One structure, one spectrum, one adduct -> MetFrag explains the peaks (metfrag_engine.run_metfrag), each explained fragment is placed on the
structure, and the answer is drawn three ways: the structure with every shown fragment in its own colour, the spectrum with the explained peaks
in the same colours (m/z and mass error on each), and one structure per fragment with its atoms highlighted.

Placing a fragment is a guess made after MetFrag: MetFrag reports the fragment as a SMILES (with the hydrogens it moved), not as atom numbers.
The skeleton of that SMILES (elements and connections only) is searched in the structure; if several places fit, the one whose hydrogen counts
agree best is drawn and the number of equally good alternatives is reported. A fragment whose atoms do not add up to its ion formula, or
that fits nowhere, is not placed: its own structure is drawn instead. Nothing here changes which peaks are explained.

The peaks are the ones the isomer check considers (query()): at least 1 % of the base peak, without the precursor and its isotope peaks.
"""
import json
import math
import re
import time
from collections import Counter

from rdkit import Chem
from rdkit.Chem import rdqueries
from rdkit.Chem.Draw import rdMolDraw2D

from mzdesk import adducts as adduct_names
from mzdesk.exclusivity import isomers
from mzdesk.ms2 import metfrag_engine as mf
from mzdesk.shared import colors as theme
from mzdesk.shared import depict, ms2_plot

MAX_QUERY_PEAKS = 2000          # same cap as the PubChem isomer comparison
MAX_FRAGMENTS = 8               # explained peaks drawn in colour (the most intense); the others stay dark in the spectrum
MAX_MATCHES = 2000              # substructure matches examined per fragment
OVERVIEW_SIZE = (440, 330)
CARD_SIZE = (250, 190)
HIGHLIGHT_ALPHA = 0.6
PALETTE = list(theme.FRAGMENT_COLORS)            # equal-weight colours (mzdesk/shared/colors.py); the labels use the darker FRAGMENT_TEXT
Y_HEADROOM = 1.45               # y axis top = tallest peak x this: room for the two-line labels


def _finite(x):
    return isinstance(x, (int, float)) and math.isfinite(x)


def considered_peaks(peaks, precursor_mz, adduct, tol_da, tol_ppm, max_peaks=MAX_QUERY_PEAKS):
    """[(m/z, % of the base peak)] the isomer check considers: >= 1 % of the base peak, the precursor and its isotope peaks left out
    (the same rule as pubchem_isomers.discriminate, so counts agree with the Isomers table)."""
    clean = [(float(mz), float(i)) for mz, i in peaks or () if _finite(mz) and _finite(i) and mz > 0 and i > 0]
    if not clean or not precursor_mz:
        return []
    base = max(i for _, i in clean)
    kept = sorted(sorted((p for p in clean if p[1] >= isomers.MIN_PEAK * base), key=lambda p: -p[1])[:max_peaks])
    z = isomers._charge(adduct_names.canonical(adduct) or adduct)
    window = max(tol_da, precursor_mz * tol_ppm * 1e-6)
    return [(mz, 100.0 * i / base) for mz, i in kept
            if not any(abs(mz - (precursor_mz + k * isomers.ISOTOPE_STEP / z)) <= window for k in range(isomers.ISOTOPE_PEAKS + 1))
            and (z > 1 or mz < precursor_mz)]


def query(peaks, precursor_mz, adduct, *, tol_da, tol_ppm):
    """The spectrum half of a fragment-view request, JSON-able: {"peaks": [[m/z, % of the base peak]], "precursor_mz", "tol_da"};
    None when there is no precursor or no fragment peak left. The structure half is the table cell (api/blocks.fragment_cell)."""
    kept = considered_peaks(peaks, precursor_mz, adduct, tol_da, tol_ppm)
    if not kept:
        return None
    return {"peaks": [[round(mz, 5), round(rel, 4)] for mz, rel in kept], "precursor_mz": float(precursor_mz), "tol_da": float(tol_da)}


def metfrag_adduct(adduct):
    """The adduct in MetFrag's spelling, or None when MetFrag cannot use it."""
    name = adduct_names.engine_name(adduct, "metfrag") if adduct else None
    return name if name in mf.ADDUCTS else None


# ---------------------------------------------------------------- placing a fragment on the structure
def _skeleton_query(fragment):
    """RDKit query of a fragment's heavy-atom skeleton: atoms match by element, bonds by connection only (MetFrag's SMILES moves hydrogens
    and may change bond orders). Returns (query, atoms of the fragment in query order)."""
    atoms = [a for a in fragment.GetAtoms() if a.GetAtomicNum() > 1]
    index = {a.GetIdx(): n for n, a in enumerate(atoms)}
    q = Chem.RWMol()
    for a in atoms:
        q.AddAtom(rdqueries.AtomNumEqualsQueryAtom(a.GetAtomicNum()))
    any_bond = Chem.MolFromSmarts("*~*").GetBondWithIdx(0)
    for b in fragment.GetBonds():
        i, j = index.get(b.GetBeginAtomIdx()), index.get(b.GetEndAtomIdx())
        if i is not None and j is not None:
            q.ReplaceBond(q.AddBond(i, j, Chem.BondType.SINGLE) - 1, any_bond)
    return q, atoms


def fragment_formula(metfrag_ion):
    """The formula of the fragment itself in MetFrag's ion string, without the hydrogens it moved and without the adduct it carries:
    '[C5H10O4]+HCOO-' -> 'C5H10O4', '[C3H7O4P+H]-' -> 'C3H7O4P'. '' if the string is not one."""
    m = re.match(r"\s*\[([^\]]+)\]", metfrag_ion or "")
    if not m:
        return ""
    shift = re.match(r"^(.*?)[+-]\d*H$", m.group(1).strip())
    return (shift.group(1) if shift and shift.group(1) else m.group(1)).strip()


def fragment_atoms(parent, fragment_smiles, fragment_formula=""):
    """Where a MetFrag fragment sits on the structure -> {"atoms": sorted atom indices, "alternatives": other equally good places}
    or None (no fragment structure, no fit, or the atoms do not add up to the fragment's formula, see fragment_formula()).
    Among the places that fit, the one whose hydrogen counts differ least from the fragment's is taken (ties: lowest atom numbers)."""
    if not fragment_smiles:
        return None
    fragment = Chem.MolFromSmiles(fragment_smiles) or Chem.MolFromSmiles(fragment_smiles, sanitize=False)
    if fragment is None:
        return None
    fragment.UpdatePropertyCache(strict=False)
    query, atoms = _skeleton_query(fragment)
    if not atoms:
        return None
    ion = {el: n for el, n in mf._formula_counts(fragment_formula).items() if el != "H"}
    scored = []
    for match in parent.GetSubstructMatches(query, uniquify=True, maxMatches=MAX_MATCHES):
        if fragment_formula:
            have = Counter(parent.GetAtomWithIdx(i).GetSymbol() for i in match)
            if have["C"] != ion.get("C", 0) or any(n > ion.get(el, 0) for el, n in have.items()):
                continue
        h_gap = sum(abs(atoms[n].GetTotalNumHs() - parent.GetAtomWithIdx(i).GetTotalNumHs()) for n, i in enumerate(match))
        scored.append((h_gap, tuple(sorted(match))))
    if not scored:
        return None
    best = min(scored)
    return {"atoms": list(best[1]), "alternatives": sum(1 for s in scored if s[0] == best[0]) - 1}


def shown_rows(rows):
    """The explained fragment peaks that get a colour: the MAX_FRAGMENTS most intense, in m/z order. rows: MetFrag annotation rows."""
    return sorted(sorted(rows, key=lambda r: -r.get("intensity", 0))[:MAX_FRAGMENTS], key=lambda r: r["mz"])


def atom_coverage(smiles, rows, *, allow_charge=False):
    """How much of the structure the coloured fragments cover -> {"atoms": covered, "heavy": all, "fraction"} (non-hydrogen atoms), or None when
    the structure cannot be read. rows: the explained fragment rows of a MetFrag run; the same rule as the fragment view (shown_rows), so the
    Isomers table and the pop-up give the same number. A fragment that cannot be placed covers nothing."""
    parent, _dropped, _error = mf._prepare_molecule(smiles, allow_charge=allow_charge)
    if parent is None or not parent.GetNumHeavyAtoms():
        return None
    covered = set()
    for row in shown_rows(rows):
        place = fragment_atoms(parent, row.get("smiles", ""), fragment_formula(row.get("metfrag_ion", "")))
        if place:
            covered.update(place["atoms"])
    heavy = parent.GetNumHeavyAtoms()
    return {"atoms": len(covered), "heavy": heavy, "fraction": round(len(covered) / heavy, 4)}


# ---------------------------------------------------------------- drawing
def _rgba(color):
    return theme._rgb(color) + (HIGHLIGHT_ALPHA,)


def _draw(mol, groups, size):
    """SVG of the (prepared) structure in pure black with each group (hex colour, atom indices) highlighted in its colour; bonds inside a group too."""
    atom_colors, bond_colors = {}, {}
    for color, atoms in groups:
        inside = set(atoms)
        for a in atoms:
            atom_colors.setdefault(a, []).append(_rgba(color))
        for b in mol.GetBonds():
            if b.GetBeginAtomIdx() in inside and b.GetEndAtomIdx() in inside:
                bond_colors.setdefault(b.GetIdx(), []).append(_rgba(color))
    drawer = depict.svg_drawer(*size)
    opts = drawer.drawOptions()
    opts.useBWAtomPalette()                         # the whole structure in one ink: oxygen red and nitrogen blue would fight the colour blocks
    opts.fillHighlights = True
    opts.highlightBondWidthMultiplier = 14
    drawer.DrawMoleculeWithHighlights(mol, "", atom_colors, bond_colors, {a: 0.38 for a in atom_colors}, {}, -1)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def _figure(ann, colors, theo_mz, adduct, y_top):
    """Spectrum: shown fragments in their colour, other explained peaks dark, the rest light; dotted precursor line."""
    import plotly.graph_objects as go

    fig = go.Figure()
    plain = [(a, "dark" if a["annotated"] else "light") for i, a in enumerate(ann) if i not in colors]
    for key, opacity, name in (("light", 0.28, "unexplained"), ("dark", 0.8, "explained")):          # the whole spectrum in pure ink, paler when unexplained
        rows = [a for a, k in plain if k == key]
        if rows:
            fig.add_trace(ms2_plot.stick_trace([a["mz"] for a in rows], [a["rel"] for a in rows], theme.PURE_INK, name).update(opacity=opacity))
    for i, color in colors.items():
        a = ann[i]
        fig.add_trace(ms2_plot.stick_trace([a["mz"]], [a["rel"]], color, f"{a['mz']:.4f}"))
    xs = [a["mz"] for a in ann]
    lo, hi = min(xs + [theo_mz]), max(xs + [theo_mz])
    x_range = [max(0.0, lo - max(0.1 * (hi - lo), 25.0)), hi + max(0.05 * (hi - lo), 5.0)]       # room on the left for a label beside the first peak
    fig.add_trace(ms2_plot.hover_markers(xs, [a["rel"] for a in ann]))
    fig.update_layout(**ms2_plot.base_layout(
        430, x_range=x_range, y_range=[0, y_top], y_title="Relative intensity (%)",
        shapes=[ms2_plot.precursor_line(theo_mz, 0, y_top)],
        annotations=[ms2_plot.precursor_label(theo_mz, f"{adduct} {theo_mz:.4f}", x_range)]))
    return json.loads(fig.to_json())


def _ppm(value):
    """'+0.4' / '-0.7'; a value that rounds to zero is '0.0', not '-0.0'."""
    text = f"{value:+.1f}"
    return "0.0" if float(text) == 0 else text


def _card(a, color=None):
    """Hover card of one peak (Figure.tsx PeakCard)."""
    rows = []
    if a["ion_formula"]:
        rows.append(["Ion formula", a["ion_formula"]])
    if a["calc_mz"] is not None:
        rows.append(["Theor. m/z", f"{a['calc_mz']:.4f}"])
    if a["ppm"] is not None:
        rows.append(["Mass error", f"{_ppm(a['ppm'])} ppm"])
    if a["annotated"] and a["nl"] is not None:
        rows.append(["Neutral loss", f"{a['nl']:.4f}"])
    rows.append(["Rel. intensity", f"{a['rel']:.1f} %"])
    card = {"mz": a["mz"], "rel": a["rel"], "on": a["annotated"], "kind": "Explained fragment" if a["annotated"] else "Not explained", "rows": rows}
    if a["b64"]:
        card["img"] = f"data:image/png;base64,{a['b64']}"
    if color:
        card["color"] = color
        card["ink"] = theme.FRAGMENT_TEXT[PALETTE.index(color)]
        card["label"] = f"{a['mz']:.4f}" + (f"<br>Δ {_ppm(a['ppm'])} ppm" if a["ppm"] is not None else "")
    return card


# ---------------------------------------------------------------- the view
def build(smiles, peaks, adduct, precursor_mz, *, tol_da=0.01, tol_ppm=0.0, tree_depth=2, name="", progress=None):
    """Run MetFrag on one structure and build the fragment view (JSON-able).

    peaks: [(m/z, intensity)] as query() returns them. adduct: MetFrag's spelling. progress(fraction 0..1, text) is called as the work goes on.
    Raises ValueError (bad input) or metfrag_engine.MetFragError (engine failure). Returns
      {"name", "smiles", "formula", "adduct", "ion_formula", "precursor_mz", "theo_mz", "peaks", "explained", "explained_fraction",
       "overview": SVG, "fragments": [{"color", "mz", "rel", "ion_formula", "calc_mz", "ppm", "nl", "svg", "placed", "alternatives", "img"}],
       "figure": {"plotly", "peaks"}, "notes": [str], "seconds"}
    fragments are the shown ones (at most MAX_FRAGMENTS, the most intense explained peaks) in m/z order; "img" is the fragment's own
    structure (PNG data URL), given only when it could not be placed.
    """
    t0 = time.monotonic()
    tick = progress or (lambda fraction, text: None)
    tick(0.02, "Checking the structure")
    parent, _dropped, error = mf._prepare_molecule(smiles, allow_charge=(adduct == "[M]+"))
    if parent is None:
        raise ValueError(error)

    res = mf.run_metfrag(smiles, peaks, adduct, precursor_mz, tree_depth, tol_da, tol_ppm, progress=tick)
    ann = res["annotation"]
    notes = list(res["warnings"])
    if not res["scored"]:
        notes.append("MetFrag returned no result for this structure, so no peak could be explained.")
    explained = [i for i, a in enumerate(ann) if a["annotated"] and a["kind"] == "fragment"]
    total = sum(a["intensity"] for a in ann) or 1.0

    # the shown fragments: the most intense explained peaks, colours in m/z order
    shown = [ann.index(r) for r in shown_rows([ann[i] for i in explained])]
    if len(explained) > len(shown):
        notes.append(f"Colour is given to the {len(shown)} most intense of {len(explained)} explained peaks; the others are dark in the spectrum.")
    colors = {i: PALETTE[n] for n, i in enumerate(shown)}

    tick(0.82, "Placing fragments on the structure")
    drawn = Chem.Mol(parent)
    rdMolDraw2D.PrepareMolForDrawing(drawn)        # one set of coordinates for every drawing, so all views have the same orientation
    fragments, placed_groups, covered_atoms = [], [], set()
    for n, i in enumerate(shown):
        a = ann[i]
        place = fragment_atoms(parent, a["smiles"], fragment_formula(a["metfrag_ion"]))
        item = {"color": colors[i], "mz": a["mz"], "rel": round(a["rel"], 2), "ion_formula": a["ion_formula"], "calc_mz": a["calc_mz"], "ppm": a["ppm"],
                "nl": a["nl"], "placed": place is not None, "alternatives": place["alternatives"] if place else 0}
        if place:
            item["svg"] = _draw(drawn, [(colors[i], place["atoms"])], CARD_SIZE)
            placed_groups.append((colors[i], place["atoms"]))
            covered_atoms.update(place["atoms"])
        else:
            item["img"] = f"data:image/png;base64,{a['b64']}" if a["b64"] else ""
        fragments.append(item)
        tick(0.82 + 0.15 * (n + 1) / len(shown), "Placing fragments on the structure")
    unplaced = [f for f in fragments if not f["placed"]]
    if unplaced:
        notes.append(f"{len(unplaced)} fragment{'s' if len(unplaced) != 1 else ''} could not be placed on the structure (no fit, or the atoms do not "
                     "match the ion formula); the fragment's own structure is drawn instead.")
    ambiguous = [f for f in fragments if f["alternatives"]]
    if ambiguous:
        notes.append("Where the same fragment fits several places (symmetric or repeated units), the drawing marks the one whose hydrogen counts "
                     "agree best; the card says how many other places fit equally well.")

    tick(0.98, "Drawing the spectrum")
    y_top = max([a["rel"] for a in ann] + [1.0]) * Y_HEADROOM
    cards = [_card(a, colors.get(i)) for i, a in enumerate(ann)]
    figure = {"plotly": _figure(ann, colors, res["theo_precursor_mz"], adduct, y_top), "peaks": cards}
    tick(1.0, "Done")
    return {
        "name": name, "smiles": Chem.MolToSmiles(parent), "formula": res["mol_formula"], "adduct": adduct, "ion_formula": res["ion_formula"],
        "precursor_mz": res["precursor_mz"], "theo_mz": res["theo_precursor_mz"], "positive": res["positive"],
        "peaks": len(ann), "explained": len(explained),
        "explained_fraction": round(sum(ann[i]["intensity"] for i in explained) / total, 4),
        "coverage": {"atoms": len(covered_atoms), "heavy": parent.GetNumHeavyAtoms(), "fraction": round(len(covered_atoms) / max(1, parent.GetNumHeavyAtoms()), 4)},
        "overview": _draw(drawn, placed_groups, OVERVIEW_SIZE), "fragments": fragments, "figure": figure, "notes": notes,
        "settings": {"tol_da": tol_da, "tol_ppm": tol_ppm, "tree_depth": tree_depth, "engine": "MetFrag 2.5.2"},
        "seconds": round(time.monotonic() - t0, 2),
    }
