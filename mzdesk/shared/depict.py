"""House style for RDKit structure drawings: thin bonds in the UI ink colour, muted heteroatom colours from the app palette,
no background (the page shows through). Use `svg_drawer(w, h)` instead of MolDraw2DSVG, or `style(drawer)` on any drawer."""
from rdkit.Chem.Draw import rdMolDraw2D

from mzdesk.shared import colors as theme

# element -> colour; everything not listed (C, H, …) is ink
ATOM_COLORS = {7: "#2D5BD0", 8: "#C2452D", 9: "#2F8F5B", 15: "#C27C1A", 16: "#A88412", 17: "#2F8F5B", 35: "#8E3B2A", 53: "#6B3FA0"}


def style(drawer):
    opts = drawer.drawOptions()
    ink = theme._rgb(theme.INK)
    palette = {-1: ink, 0: ink, 1: ink, 6: ink}
    palette.update({z: theme._rgb(c) for z, c in ATOM_COLORS.items()})
    opts.setAtomPalette(palette)
    opts.bondLineWidth = 1.6
    opts.clearBackground = False
    opts.addStereoAnnotation = False
    opts.padding = 0.08
    opts.annotationColour = theme._rgb(theme.MUTED)
    opts.highlightRadius = 0.32
    opts.fillHighlights = True
    return drawer


def svg_drawer(width, height):
    return style(rdMolDraw2D.MolDraw2DSVG(width, height))
