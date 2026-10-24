"""lipid_metabolite_viz.py

Render a 3-panel report for a classification Result:
  Panel 1 - molecular structure with substructures highlighted by category
  Panel 2 - feature numeric panel (logP, MW + detection-window bar, counts)
  Panel 3 - decision panel (category + rule-hit reasons)
"""
from __future__ import annotations

import io
import textwrap
from typing import Optional, Dict, List

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import font_manager

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from lipid_metabolite_classifier import Result, Features, Thresholds


# RGB (0..1) colors for the four highlight categories.
# Harmonised with theme.FRAG_HL (cool-white + indigo blue design system).
HIGHLIGHT_COLORS: Dict[str, tuple] = {
    "acyl":       (0.106, 0.227, 0.561),  # indigo    - >=C12 acyl / long chain
    "isoprenoid": (0.055, 0.486, 0.608),  # teal      - sterol / isoprenoid
    "head":       (0.333, 0.400, 0.800),  # periwinkle- polar head group
    "glycosyl":   (0.902, 0.624, 0.208),  # amber     - sugar ring (warm anchor)
}

# Display color per output category.
CATEGORY_COLORS: Dict[str, str] = {
    "脂质组":   "#2f6fb0",
    "代谢组":   "#2f9e5f",
    "都有可能": "#c08a1e",
    "都不可能": "#8a8a8a",
    "无效SMILES": "#b03030",
}

# CJK font candidates (macOS first), tried in order.
_CJK_FONTS: List[str] = [
    "PingFang SC", "Heiti SC", "STHeiti", "Arial Unicode MS",
    "Songti SC", "Hiragino Sans GB", "Noto Sans CJK SC",
    "Microsoft YaHei", "SimHei",
]


def _setup_cjk_font() -> Optional[str]:
    """Pick the first installed CJK font and set matplotlib rcParams.

    ALGORITHM: collect installed font names via font_manager; return the first
    name in _CJK_FONTS that is installed and set plt.rcParams['font.sans-serif']
    to put it first, and rcParams['axes.unicode_minus']=False. Return the chosen
    name (or None if none found)."""
    installed = set()
    try:
        installed = {f.name for f in font_manager.fontManager.ttflist}
    except Exception:
        pass

    # Fallback scan in case the font manager cache is not yet populated.
    if not any(name in installed for name in _CJK_FONTS):
        try:
            for path in font_manager.findSystemFonts():
                try:
                    installed.add(font_manager.FontProperties(fname=path).get_name())
                except Exception:
                    continue
        except Exception:
            pass

    for name in _CJK_FONTS:
        if name in installed:
            current = list(matplotlib.rcParams.get("font.sans-serif", []))
            matplotlib.rcParams["font.sans-serif"] = [name] + [f for f in current if f != name]
            matplotlib.rcParams["axes.unicode_minus"] = False
            return name
    return None


def category_color(category: str) -> str:
    """Return the hex color for a category (default gray)."""
    return CATEGORY_COLORS.get(category, "#8a8a8a")


def draw_structure_highlight(result: Result, size=(560, 420)):
    """Render the molecule with substructures highlighted by category color.

    Returns a numpy RGB image array (H, W, 3) uint8 suitable for imshow.
    ALGORITHM:
      - if not result.features.valid: return a small white image array.
      - mol = Chem.MolFromSmiles(result.features.smiles).
      - Build highlightAtoms list and highlightAtomColors dict from
        result.features.highlight. Priority when an atom is in several categories:
        head > glycosyl > isoprenoid > acyl (so head/sugar colors win).
      - Also highlight bonds whose both endpoints share the same category
        (compute per category; put into highlightBonds + highlightBondColors).
      - Use rdMolDraw2D.MolDraw2DCairo(size[0], size[1]); drawOptions:
        addStereoAnnotation False; set useBWAtomPalette off (default colors ok).
        DrawMolecule(mol, highlightAtoms=..., highlightAtomColors=...,
        highlightBonds=..., highlightBondColors=...); FinishDrawing().
      - Convert PNG bytes (GetDrawingText()) via matplotlib.image or PIL to an
        RGB numpy array and return it. (Use io.BytesIO + plt.imread; take [:,:,:3];
        if float in 0..1, scale to 0..255 uint8.)
    """
    features = getattr(result, "features", None)
    if features is None or not getattr(features, "valid", False):
        return np.full((size[1] // 2, size[0] // 2, 3), 255, dtype=np.uint8)

    smi = getattr(features, "smiles", None)
    mol = Chem.MolFromSmiles(smi) if smi else None
    if mol is None:
        return np.full((size[1] // 2, size[0] // 2, 3), 255, dtype=np.uint8)

    highlight = getattr(features, "highlight", None) or {}
    if not isinstance(highlight, dict):
        try:
            highlight = dict(highlight)
        except Exception:
            highlight = {}

    atom_sets: Dict[str, set] = {}
    for cat, vals in highlight.items():
        if cat not in HIGHLIGHT_COLORS:
            continue
        try:
            if isinstance(vals, (int, np.integer)):
                atom_sets[cat] = {int(vals)}
            else:
                atom_sets[cat] = {int(v) for v in vals}
        except Exception:
            atom_sets[cat] = set()

    priority = ["acyl", "isoprenoid", "glycosyl", "head"]

    highlight_atoms: List[int] = []
    highlight_atom_colors: Dict[int, tuple] = {}
    for cat in priority:
        color = HIGHLIGHT_COLORS[cat]
        for aidx in atom_sets.get(cat, set()):
            if 0 <= aidx < mol.GetNumAtoms():
                if aidx not in highlight_atoms:
                    highlight_atoms.append(aidx)
                highlight_atom_colors[aidx] = color

    highlight_bonds: List[int] = []
    highlight_bond_colors: Dict[int, tuple] = {}
    for cat in priority:
        color = HIGHLIGHT_COLORS[cat]
        atoms = atom_sets.get(cat, set())
        if not atoms:
            continue
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in atoms and a2 in atoms:
                bidx = bond.GetIdx()
                if bidx not in highlight_bonds:
                    highlight_bonds.append(bidx)
                highlight_bond_colors[bidx] = color

    return _render(mol, size, highlight_atoms, highlight_atom_colors,
                   highlight_bonds, highlight_bond_colors)


def draw_structure(smiles: str, size=(560, 420)):
    """Plain 2D depiction without substructure highlights (RGB uint8 array)."""
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return np.full((size[1] // 2, size[0] // 2, 3), 255, dtype=np.uint8)
    return _render(mol, size, [], {}, [], {})


def _render(mol, size, highlight_atoms, highlight_atom_colors,
            highlight_bonds, highlight_bond_colors):
    """Draw `mol` with RDKit's Cairo backend and return an RGB uint8 array."""
    w, h = int(size[0]), int(size[1])
    drawer = rdMolDraw2D.MolDraw2DCairo(w, h)
    opts = drawer.drawOptions()
    opts.addStereoAnnotation = False
    try:
        opts.useBWAtomPalette = False
    except Exception:
        pass

    drawer.DrawMolecule(
        mol,
        highlightAtoms=highlight_atoms,
        highlightAtomColors=highlight_atom_colors,
        highlightBonds=highlight_bonds,
        highlightBondColors=highlight_bond_colors,
    )
    drawer.FinishDrawing()

    img = plt.imread(io.BytesIO(drawer.GetDrawingText()))
    if img.ndim == 3:
        arr = img[:, :, :3]
    else:
        arr = img
    if arr.dtype.kind == "f":
        arr = (np.clip(arr, 0.0, 1.0) * 255).astype(np.uint8)
    else:
        arr = np.asarray(arr, dtype=np.uint8)
    return arr


def _draw_feature_panel(ax, result: Result, thr: Optional[Thresholds] = None):
    """Draw the numeric feature panel onto a matplotlib Axes `ax`.

    Shows, as a clean text table (ax.axis('off') + ax.text lines):
      logP, 分子量(MW), 最长脂肪链碳数, 是否≥C12酰基链, 异戊二烯骨架(类型),
      极性头基数(类型), 糖基数.
    Then draw a horizontal 'detection window' bar spanning [thr.mw_min, thr.mw_max]
    with a marker at MW (green if inside window, red if outside). Use a small
    inset or the lower part of the axes. Title: '特征指标'."""
    ax.set_title("特征指标")
    ax.axis("off")

    features = getattr(result, "features", None)

    def _get(*names):
        if features is None:
            return None
        for n in names:
            if hasattr(features, n):
                v = getattr(features, n)
                if v is not None:
                    return v
        return None

    def _fmt(v):
        if v is None:
            return "N/A"
        if isinstance(v, (bool, np.bool_)):
            return "是" if bool(v) else "否"
        if isinstance(v, (float, np.floating)):
            return f"{float(v):.2f}"
        return str(v)

    def _is_true(v):
        if isinstance(v, str):
            return v.strip().lower() in ("y", "yes", "true", "1", "是")
        return bool(v)

    logp = _get("logp")
    mw = _get("mw")
    max_chain = _get("longest_chain_carbons")
    is_long = _get("has_long_acyl")
    has_iso = _get("has_isoprenoid")
    iso_type = _get("isoprenoid_type")
    n_head = _get("n_polar_head")
    head_types = _get("polar_head_types")
    n_gly = _get("n_glycosyl")

    is_long_text = "N/A" if is_long is None else ("是" if _is_true(is_long) else "否")
    if not has_iso:
        iso_text = "否"
    else:
        iso_text = f"是（{iso_type}）" if iso_type else "是"
    if n_head is None:
        head_text = "N/A"
    elif head_types:
        head_text = f"{n_head}（{'、'.join(map(str, head_types))}）"
    else:
        head_text = str(n_head)

    rows = [
        ("logP", logp),
        ("分子量(MW)", mw),
        ("最长脂肪链碳数", max_chain),
        ("≥C12酰基链", is_long_text),
        ("异戊二烯骨架", iso_text),
        ("极性头基", head_text),
        ("糖基数", n_gly),
    ]

    y = 0.96
    for label, value in rows:
        ax.text(
            0.05, y, f"{label}：{_fmt(value)}",
            transform=ax.transAxes, va="top", fontsize=9.5,
        )
        y -= 0.10

    mw_min = getattr(thr, "mw_min", None) if thr is not None else None
    mw_max = getattr(thr, "mw_max", None) if thr is not None else None

    if mw_min is not None and mw_max is not None:
        lo, hi = sorted([mw_min, mw_max])
        inset = ax.inset_axes([0.12, 0.06, 0.76, 0.13])
        inset.axvspan(lo, hi, color="#4c72b0", alpha=0.15)
        inset.axhline(0, color="#888888", lw=3)
        inset.set_xlim(lo, hi)
        inset.set_ylim(-0.6, 0.6)
        inset.set_yticks([])
        inset.set_xticks([lo, hi])
        inset.set_xticklabels([f"{lo:.1f}", f"{hi:.1f}"], fontsize=8)
        inset.spines["top"].set_visible(False)
        inset.spines["right"].set_visible(False)
        inset.spines["left"].set_visible(False)
        inset.set_title("检测窗口", fontsize=9)

        if mw is not None:
            color = "#2f9e5f" if lo <= mw <= hi else "#d43f3f"
            inset.plot(mw, 0, marker="o", color=color, markersize=9, zorder=5)
    else:
        ax.text(
            0.5, 0.10, "检测窗口：阈值未提供",
            transform=ax.transAxes, ha="center", fontsize=9, color="#888888",
        )


def _draw_decision_panel(ax, result: Result):
    """Draw the decision panel onto a matplotlib Axes `ax`.

    ax.axis('off'); put a large colored category label (color=category_color)
    at top ('判定：<category>'), then list result.reasons as wrapped bullet
    lines below. Title: '判定理由'."""
    ax.set_title("判定理由")
    ax.axis("off")

    category = getattr(result, "category", "") or ""
    ax.text(
        0.5, 0.88, f"判定：{category}",
        transform=ax.transAxes, ha="center", va="top",
        fontsize=16, fontweight="bold", color=category_color(category),
    )

    reasons = getattr(result, "reasons", None) or []
    y = 0.80
    step = 0.052
    max_reasons = 10

    for idx, reason in enumerate(reasons[:max_reasons]):
        if y < 0.02:
            break
        lines = textwrap.wrap(str(reason), width=30) or [""]
        for i, line in enumerate(lines):
            if y < 0.02:
                break
            prefix = "• " if i == 0 else "   "
            ax.text(
                0.04, y, prefix + line,
                transform=ax.transAxes, va="top", fontsize=8.5, color="#222222",
            )
            y -= step

    if len(reasons) > max_reasons:
        ax.text(0.05, y, "…", transform=ax.transAxes, fontsize=12)


def render_report(result: Result, thr: Optional[Thresholds] = None, figsize=(12, 6)):
    """Compose the full 3-panel figure and return the matplotlib Figure.

    Layout: left column = structure image (imshow of draw_structure_highlight,
    axis off, title '结构 + 高亮'); right column split top/bottom into
    _draw_feature_panel and _draw_decision_panel. Add a legend mapping the four
    highlight colors to 酰基链/异戊二烯/极性头基/糖基. Call _setup_cjk_font()
    first. Use fig.suptitle with the SMILES. tight_layout. Return fig."""
    _setup_cjk_font()
    if thr is None:
        thr = Thresholds()

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(
        2, 2,
        width_ratios=[1, 1.15],
        height_ratios=[1, 1],
        hspace=0.35,
        wspace=0.2,
    )

    ax_struct = fig.add_subplot(gs[:, 0])
    ax_feature = fig.add_subplot(gs[0, 1])
    ax_decision = fig.add_subplot(gs[1, 1])

    img = draw_structure_highlight(result)
    ax_struct.imshow(img)
    ax_struct.axis("off")
    ax_struct.set_title("结构 + 高亮")

    from matplotlib.patches import Patch

    legend_handles = [
        Patch(facecolor=HIGHLIGHT_COLORS["acyl"], label="酰基链"),
        Patch(facecolor=HIGHLIGHT_COLORS["isoprenoid"], label="异戊二烯"),
        Patch(facecolor=HIGHLIGHT_COLORS["head"], label="极性头基"),
        Patch(facecolor=HIGHLIGHT_COLORS["glycosyl"], label="糖基"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=4,
        fontsize=9,
        frameon=False,
    )

    features = getattr(result, "features", None)
    smiles = getattr(features, "smiles", "") if features is not None else ""
    fig.suptitle(f"SMILES: {smiles}", fontsize=13)

    _draw_feature_panel(ax_feature, result, thr)
    _draw_decision_panel(ax_decision, result)

    # Manual margins (tight_layout is incompatible with the inset axes).
    fig.subplots_adjust(left=0.03, right=0.97, top=0.90, bottom=0.12,
                        wspace=0.15, hspace=0.40)
    return fig
