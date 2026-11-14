"""示意结构的画法（RDKit SVG）：猜出来的部分高亮，图下一行小色标。

颜色用主题色系（蓝 / 青 / 紫的浅色调），不用琥珀、珊瑚红：强调色留给警告和报错。浅色调在深色主题下由前端 darkColors 自动换成深色调。
色标和文字都由 RDKit 画成路径，不依赖浏览器字体。
"""
from rdkit import Chem
from rdkit.Geometry import Point2D

from mzdesk.shared import colors as theme, depict

HL_DB = theme._rgb(theme.DATA_BLUE, .2)          # 默认 C=C 位置（亮度低于 0.62，深色主题下 darkColors 会提亮而不是压暗，两种主题都看得清）
HL_FG = theme._rgb(theme.DATA_TEAL, .25)         # 默认官能团位置
HL_SPLIT = theme._rgb(theme.DATA_VIOLET, .55)    # 按默认规则拆出来的链（只标原子，不描键：整条链描键太重）
KEY = (("db", HL_DB, "C=C position assumed"), ("fg", HL_FG, "Group position assumed"), ("split", HL_SPLIT, "Chains split by default rule"))

UNIT_PX = 15                  # 每个坐标单位（1.5 = 一个键长）的像素数
MIN_W, MAX_W, MIN_H = 420, 980, 150
PAD = 20
LEGEND_FONT, LEGEND_ROW, SWATCH_R, ENTRY_GAP = 13, 24, 5.5, 26


def _highlights(ill):
    """原子和键的高亮颜色；优先级：官能团 > C=C > 拆分的链。C=C 和官能团连键一起描（双键、C=O 看得清），拆分的链只标原子。"""
    mol = ill.mol
    atoms, bonds = {}, {}
    for atom_ids, color, with_bonds in ((ill.split_atoms, HL_SPLIT, False), (ill.db_atoms, HL_DB, True), (ill.fg_atoms, HL_FG, True)):
        ids = set(atom_ids)
        for i in ids:
            atoms[i] = color
        if with_bonds:
            for b in mol.GetBonds():
                if b.GetBeginAtomIdx() in ids and b.GetEndAtomIdx() in ids:
                    bonds[b.GetIdx()] = color
    return atoms, bonds


def _entries(ill):
    used = {"db": ill.db_atoms, "fg": ill.fg_atoms, "split": ill.split_atoms}
    return [(color, label) for name, color, label in KEY if used[name]]


def _legend_rows(entries, width):
    """色标按宽度折行；每项宽度按字数估（路径文字，没有现成的量宽接口）。"""
    rows, row, x = [], [], PAD
    for color, label in entries:
        w = 2 * SWATCH_R + 8 + 0.56 * LEGEND_FONT * len(label)
        if row and x + w > width - PAD:
            rows.append(row)
            row, x = [], PAD
        row.append((x, color, label))
        x += w + ENTRY_GAP
    if row:
        rows.append(row)
    return rows


def illustration_svg(ill):
    """Illustration -> SVG 字符串。画布按分子的长宽比定，图下留一行（或几行）色标。"""
    mol = Chem.Mol(ill.mol)
    Chem.rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    xs = [conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())]
    ys = [conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())]
    ext_w, ext_h = max(xs) - min(xs), max(ys) - min(ys)
    width = int(min(MAX_W, max(MIN_W, ext_w * UNIT_PX + 2 * PAD)))
    scale = (width - 2 * PAD) / max(ext_w, 1e-6)
    height = int(max(MIN_H, ext_h * scale + 2 * PAD))
    rows = _legend_rows(_entries(ill), width)
    legend_h = LEGEND_ROW * len(rows) + (10 if rows else 0)

    drawer = depict.svg_drawer(width, height + legend_h)
    opts = drawer.drawOptions()
    opts.highlightRadius = 0.3
    opts.padding = 0.03
    opts.baseFontSize = 0.72
    opts.highlightBondWidthMultiplier = 7
    if rows:
        opts.legendFraction = legend_h / (height + legend_h)
    atoms, bonds = _highlights(ill)
    drawer.DrawMolecule(mol, legend=" " if rows else "", highlightAtoms=list(atoms), highlightAtomColors=atoms,
                        highlightBonds=list(bonds), highlightBondColors=bonds)
    drawer.SetFontSize(LEGEND_FONT)
    for r, row in enumerate(rows):
        y = height + 10 + LEGEND_ROW * r + LEGEND_ROW / 2
        for x, color, label in row:
            drawer.SetFillPolys(True)
            drawer.SetColour(color)
            drawer.DrawEllipse(Point2D(x, y - SWATCH_R), Point2D(x + 2 * SWATCH_R, y + SWATCH_R), True)
            drawer.DrawString(label, Point2D(x + 2 * SWATCH_R + 8, y), 1, True)
    drawer.FinishDrawing()
    # DrawString 的文字固定画成黑色（SetColour 和各种 colour 选项都不管它），而分子的原子标签用的是调色板颜色，所以整张图里 #000000 只有色标文字：换成次要文字色
    return drawer.GetDrawingText().replace("fill='#000000'", f"fill='{theme.MUTED}'")
