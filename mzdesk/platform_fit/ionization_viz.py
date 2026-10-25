"""Ionization page figure: the neutral molecule with its ESI-relevant ionizable sites highlighted (RDKit SVG).

Same colours as the Separation page (red: acid, teal: base, violet: both). A permanent charge is marked on its atom (teal).
Acid sites above ``acid_moderate_pka`` (alcohol OH, amide NH ...) are left out: too weak to deprotonate in ESI-.
"""
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from mzdesk.platform_fit.ionization_rules import THRESHOLDS
from mzdesk.platform_fit.separation_viz import ACID, BASE, BOTH


def sites_svg(prof, th=THRESHOLDS, size=(560, 420)):
    """SVG of ``prof.mol`` with the ionizable sites highlighted; plain depiction when there is nothing to highlight."""
    mol = Chem.Mol(prof.mol)
    by_atom = {}
    for s in prof.sites:
        if s.kind == "base" or s.pka <= th.acid_moderate_pka:
            by_atom.setdefault(s.atom, []).append(s)
    colors = {}
    for idx, group in by_atom.items():
        kinds = {s.kind for s in group}
        colors[idx] = BOTH if len(kinds) == 2 else (ACID if kinds == {"acid"} else BASE)
        mol.GetAtomWithIdx(idx).SetProp("atomNote", f"{group[0].label} " + "/".join(f"{s.pka:.1f}" for s in group))
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() > 0 and atom.GetIdx() not in colors and not any(n.GetFormalCharge() < 0 for n in atom.GetNeighbors()):
            colors[atom.GetIdx()] = BASE
            atom.SetProp("atomNote", f"{atom.GetSymbol()}{atom.GetIdx() + 1} permanent +")
    drawer = rdMolDraw2D.MolDraw2DSVG(*size)
    opts = drawer.drawOptions()
    opts.clearBackground = False
    opts.addStereoAnnotation = False
    opts.annotationFontScale = 0.75
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol, highlightAtoms=list(colors), highlightAtomColors=colors)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()
