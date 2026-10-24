"""Basic tab: chemical identity, adduct m/z, descriptors and the lipidome /
metabolome classification of a single SMILES."""
import math
from typing import List, Optional, Tuple

import streamlit as st

import theme
import molecular_properties
from lipid_metabolite_classifier import analyze, Thresholds
from lipid_metabolite_viz import draw_structure, draw_structure_highlight, HIGHLIGHT_COLORS
from display_en import (category_en, category_subtitle_en, chain_label,
                        english_reasons, feature_rows_en)


PLACEHOLDER = "— pick an example —"
DEFAULT_EXAMPLE = "Cholesterol (lipid)"

EXAMPLES = {
    PLACEHOLDER: "",
    "Cholesterol (lipid)": "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
    "Palmitic acid (lipid)": "CCCCCCCCCCCCCCCC(=O)O",
    "POPC phospholipid (lipid)": "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
    "Retinol / Vitamin A (lipid)": "OC/C=C(C)/C=C/C=C(C)/C=C/C1=C(C)CCCC1(C)C",
    "Glucose (metabolite)": "OCC1OC(O)C(O)C(O)C1O",
    "Caffeine (metabolite)": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Citric acid (metabolite)": "OC(=O)CC(O)(C(=O)O)CC(=O)O",
    "Ibuprofen (both)": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "Anthracene (both)": "c1ccc2cc3ccccc3cc2c1",
    "Water (neither)": "O",
}

LEGEND_KEYS = ["acyl", "isoprenoid", "head", "glycosyl"]
# The "acyl" legend label follows the chain-length threshold (display_en.chain_label).
LEGEND_LABELS = {
    "isoprenoid": "Isoprenoid / sterol",
    "head": "Polar head group",
    "glycosyl": "Glycosyl ring",
}

OUTSIDE_RED = theme.CATEGORY_ACCENT["Invalid"]

# Classifier thresholds as (session key, Thresholds field, label, min, max, step, help),
# one row per rule they tune. Default values come from Thresholds().
THRESHOLD_ROWS = [
    [
        ("clf_mw_min", "mw_min", "MW min (Da)", 0.0, 2000.0, 10.0,
         "Average MW below this is outside the typical LC-MS window → Neither."),
        ("clf_mw_max", "mw_max", "MW max (Da)", 100.0, 10000.0, 100.0,
         "Average MW above this is outside the typical LC-MS window → Neither."),
        ("clf_logp_polar_max", "logp_polar_max", "LogP polar max", -10.0, 20.0, 0.1,
         "Without a lipid scaffold: logP below this → Metabolome, otherwise → Both."),
    ],
    [
        ("clf_acyl_min_carbons", "acyl_min_carbons", "Acyl min carbons", 4, 30, 1,
         "Shortest mostly-hydrocarbon run in an acyclic chain that counts as a lipid tail."),
        ("clf_medium_chain_max", "medium_chain_max", "Medium chain max", 4, 30, 1,
         "Largest total carbon count still treated as a medium-chain fatty acid by the "
         "boundary (Both) rule."),
        ("clf_isoprene_min_units", "isoprene_min_units", "Isoprene min units", 1, 20, 1,
         "Methyl-bearing C=C bonds needed for an isoprenoid chain, one per isoprene unit; "
         "one saturated two-unit motif always counts, and sterol/triterpene nuclei are "
         "detected separately."),
    ],
    [
        ("clf_boundary_logp_low", "boundary_logp_low", "Boundary logP low", -10.0, 20.0, 0.1,
         "Boundary rule: a medium-chain fatty acyl without isoprenoid, polar head or sugar "
         "and with logP in [low, high] → Both."),
        ("clf_boundary_logp_high", "boundary_logp_high", "Boundary logP high", -10.0, 20.0, 0.1,
         "Upper logP bound of the boundary (Both) rule."),
        ("clf_boundary_mw_max", "boundary_mw_max", "Boundary MW max (Da)", 50.0, 5000.0, 10.0,
         "Size ceiling of the boundary (Both) rule."),
    ],
    [
        ("clf_chain_max_polar_fraction", "chain_max_polar_fraction",
         "Chain max polar fraction", 0.0, 1.0, 0.05,
         "Largest share of a chain's carbons that may carry a heteroatom and still count "
         "as a lipid tail. Keeps oxygenated fatty acids (HETEs, leukotriene B4) while "
         "rejecting polyols and sugar alcohols of the same carbon count."),
        ("clf_lipid_min_mw", "lipid_min_mw", "Lipid min MW (Da)", 0.0, 2000.0, 10.0,
         "Size floor for lipid evidence that is only a chain or a methylation count "
         "(fatty acyl, prenol). A recognised GP/GL/SP/ST backbone ignores it."),
        ("clf_polar_lipid_logp_max", "polar_lipid_logp_max", "Polar lipid logP max",
         -10.0, 20.0, 0.1,
         "At or below this logP a steroid nucleus with no lipophilic tail (hormone, "
         "bile acid, oxysterol) is reported as Both rather than Lipidome."),
    ],
]


def _rgb(t):
    r, g, b = t
    return f"rgb({int(r*255)},{int(g*255)},{int(b*255)})"


# --- input callbacks ---------------------------------------------------------
def _on_example_change():
    sel = st.session_state.get("clf_example")
    if sel and sel != PLACEHOLDER:
        st.session_state["clf_smiles"] = EXAMPLES[sel]


def _on_smiles_change():
    """A hand-edited SMILES no longer is the example shown in the picker: reset
    the picker so choosing that example again re-applies it."""
    sel = st.session_state.get("clf_example", PLACEHOLDER)
    if sel != PLACEHOLDER and st.session_state.get("clf_smiles", "").strip() != EXAMPLES.get(sel):
        st.session_state["clf_example"] = PLACEHOLDER


def _reset_thresholds():
    for row in THRESHOLD_ROWS:
        for key, *_ in row:
            st.session_state.pop(key, None)


# --- helpers -------------------------------------------------------------------
def _threshold_inputs(columns: int = 3) -> Thresholds:
    """Threshold editor. ``columns`` is the grid width — 1 for a narrow rail."""
    defaults = Thresholds()
    values = {}
    with st.expander("Classifier thresholds"):
        for row in THRESHOLD_ROWS:
            if columns > 1:
                for col, (key, field, label, lo, hi, step, help_text) in zip(
                        st.columns(columns), row):
                    values[field] = col.number_input(
                        label, min_value=lo, max_value=hi, value=getattr(defaults, field),
                        step=step, key=key, help=help_text)
            else:
                for key, field, label, lo, hi, step, help_text in row:
                    values[field] = st.number_input(
                        label, min_value=lo, max_value=hi, value=getattr(defaults, field),
                        step=step, key=key, help=help_text)
        st.button("Reset to defaults", key="clf_thr_reset", on_click=_reset_thresholds,
                  width="stretch")
    return Thresholds(**values)


def _thresholds_from_state() -> Thresholds:
    """The thresholds as they stand *before* this run draws the editor.

    The query panel shows a live structure preview whose legend follows the acyl
    chain threshold, but the editor is drawn further down the same panel — so read
    the widget values straight out of the session instead of reordering the panel.
    """
    defaults = Thresholds()
    values = {}
    for row in THRESHOLD_ROWS:
        for key, field, *_rest in row:
            values[field] = st.session_state.get(key, getattr(defaults, field))
    return Thresholds(**values)


def _tuned_thresholds(thr: Thresholds) -> List[str]:
    """Fields the user has moved off the default, as 'label value' strings."""
    defaults = Thresholds()
    moved = []
    for row in THRESHOLD_ROWS:
        for _key, field, label, *_rest in row:
            if getattr(thr, field) != getattr(defaults, field):
                moved.append(f"{label} → {getattr(thr, field):g}")
    return moved


def _remember(smiles: str, name: str, category: str) -> None:
    """Push a finished analysis onto the session history (most recent first)."""
    history = st.session_state.setdefault("clf_history", [])
    history[:] = [h for h in history if h["smiles"] != smiles]
    history.insert(0, {"smiles": smiles, "name": name, "category": category})
    del history[8:]


def _recall(smiles: str) -> None:
    """History click: put that SMILES back in the input and drop the example picker."""
    st.session_state["clf_smiles"] = smiles
    st.session_state["clf_example"] = PLACEHOLDER


def _history_panel(current: str) -> None:
    """Session history as full-width button rows, below the results.

    It used to be a vertical list in the query rail, where it was the only thing
    under a panel that ended a third of the way down the page.
    """
    history = st.session_state.get("clf_history", [])
    if len(history) < 2:
        return
    theme.section("This session", aux=f"{len(history)} analysed")
    per_row = 4
    for start in range(0, len(history), per_row):
        chunk = history[start:start + per_row]
        for col, item in zip(st.columns(per_row, gap="small"), chunk):
            with col:
                st.button(
                    f"{item['name']} · {item['category']}",
                    key=f"clf_hist_{start}_{item['smiles']}",
                    on_click=_recall,
                    args=(item["smiles"],),
                    width="stretch",
                    disabled=item["smiles"] == current,
                    help=item["smiles"],
                )


def _threshold_problems(thr: Thresholds) -> List[str]:
    """Threshold combinations that make a rule impossible to satisfy."""
    problems = []
    if thr.mw_min >= thr.mw_max:
        problems.append("MW min ≥ MW max: the detection window is empty, so every molecule is Neither.")
    if thr.medium_chain_max < thr.acyl_min_carbons:
        problems.append("Medium chain max < Acyl min carbons: the boundary (Both) rule can never apply.")
    if thr.boundary_logp_low > thr.boundary_logp_high:
        problems.append("Boundary logP low > Boundary logP high: the boundary (Both) rule can never apply.")
    if thr.boundary_mw_max < thr.mw_min:
        problems.append("Boundary MW max < MW min: the boundary (Both) rule can never apply.")
    if thr.lipid_min_mw > thr.mw_max:
        problems.append("Lipid min MW > MW max: no molecule can ever be a fatty acyl or a "
                        "prenol lipid.")
    if thr.chain_max_polar_fraction >= 1.0:
        problems.append("Chain max polar fraction ≥ 1: every carbon chain counts as a lipid "
                        "tail, so polyols and sugar alcohols read as lipids.")
    return problems


def _adduct_fields(mass: Optional[float], positive: bool) -> List[tuple]:
    """Identity-card rows (ion, m/z, copyable, Δ note) for one polarity."""
    rows = []
    for label, (delta, charge) in molecular_properties.ADDUCTS.items():
        if (charge > 0) != positive:
            continue
        mz = molecular_properties.adduct_exact_mz(mass, label)
        rows.append((label, "—" if mz is None else f"{mz:.4f}", mz is not None,
                     f"Δ {delta:+.4f}"))
    return rows


def _ion_label(charge: int) -> str:
    """'[M]+', '[M]2-', ... for an input that carries a net charge."""
    n = abs(charge)
    return f"[M]{n if n > 1 else ''}{'+' if charge > 0 else '-'}"


def _log_ticks(lo: float, hi: float) -> List[float]:
    """1-2-5 ticks strictly inside (lo, hi); only decades if that gives too many."""
    ticks = [(m, m * 10.0 ** e)
             for e in range(math.floor(math.log10(lo)), math.ceil(math.log10(hi)) + 1)
             for m in (1, 2, 5)]
    ticks = [(m, v) for m, v in ticks if lo < v < hi]
    if len(ticks) > 6:
        ticks = [(m, v) for m, v in ticks if m == 1]
    return [v for _, v in ticks]


def _window_bar_html(mw: float, lo: float, hi: float) -> Tuple[str, str]:
    """Detection-window bar for [lo, hi] Da with a marker at `mw`.

    Log scale when 0 < lo < hi (on a linear 50-2000 Da axis most metabolites
    crowd the left edge), linear otherwise. Returns (html, scale name).
    """
    use_log = 0 < lo < hi

    def pos(v: float) -> float:
        if hi <= lo:
            return 0.0
        if use_log:
            t = (math.log10(max(v, 1e-9)) - math.log10(lo)) / (math.log10(hi) - math.log10(lo))
        else:
            t = (v - lo) / (hi - lo)
        return 100.0 * min(1.0, max(0.0, t))

    inside = lo <= mw <= hi
    color = theme.BLUE if inside else OUTSIDE_RED
    status = "" if inside else (" · below window" if mw < lo else " · above window")
    p = pos(mw)
    ticks = _log_ticks(lo, hi) if use_log else []
    tick_marks = "".join(
        f'<div style="position:absolute;top:0;bottom:0;left:{pos(t):.2f}%;width:1px;'
        f'background:{theme.MUTED};opacity:.4;"></div>' for t in ticks)
    tick_labels = "".join(
        f'<span style="position:absolute;left:{pos(t):.2f}%;transform:translateX(-50%);">{t:g}</span>'
        for t in ticks if 8.0 <= pos(t) <= 92.0)
    mono = "ui-monospace,'SF Mono',Menlo,monospace"
    html = (
        '<div class="nb-mwbar" style="position:relative;padding-top:32px;margin:0 6px 4px 0;">'
        # value flag: follows the marker; translateX(-p%) keeps it inside the bar
        f'<div style="position:absolute;top:0;left:{p:.2f}%;transform:translateX(-{p:.2f}%);'
        f'white-space:nowrap;background:{color};color:#fff;border:2px solid {theme.INK};'
        f'border-radius:2px;padding:2px 7px;font:800 .72rem/1.4 {mono};">'
        f'MW {mw:.1f} Da{status}</div>'
        # track, ticks and the marker (clamped so it never pokes out of the ends)
        f'<div style="position:relative;height:22px;border:3px solid {theme.INK};'
        f'background:{theme.SURFACE_2};border-radius:2px;">{tick_marks}'
        f'<div style="position:absolute;top:-7px;width:13px;height:30px;'
        f'left:clamp(-3px, calc({p:.2f}% - 6.5px), calc(100% - 10px));'
        f'border:3px solid {theme.INK};background:{color};border-radius:2px;"></div></div>'
        # axis labels
        f'<div style="position:relative;height:1.2em;margin-top:7px;font-size:.7rem;'
        f'font-weight:700;color:{theme.MUTED};">'
        f'<span style="position:absolute;left:0;">{lo:g}</span>{tick_labels}'
        f'<span style="position:absolute;right:0;">{hi:g}</span></div></div>'
    )
    return html, ("log scale" if use_log else "linear scale")




# --- sections ------------------------------------------------------------------
STRUCTURE_VIEWS = ["2D · highlights", "3D conformer"]


def _keep_view():
    """Segmented controls allow deselection; re-select the previous view instead."""
    if st.session_state.get("clf_view") is None:
        st.session_state["clf_view"] = st.session_state.get("_clf_view_last", STRUCTURE_VIEWS[0])
    else:
        st.session_state["_clf_view_last"] = st.session_state["clf_view"]


def _display_name(smiles: str, props) -> str:
    """Short label for the session history: the example's name, else the formula."""
    for name, smi in EXAMPLES.items():
        if smi and smi == smiles:
            return name.split(" (")[0]
    if props.formula:
        return props.formula
    return smiles if len(smiles) <= 18 else smiles[:17] + "…"


def _structure_figure(smiles: str, result) -> None:
    """2D depiction with substructure highlights, or the 3D conformer on demand.

    The conformer is only embedded when that view is selected — RDKit's embedding
    is the slowest step on this tab, and most lookups never need it.
    """
    theme.subsection("Structure", aux="RDKit")
    st.segmented_control(
        "View", STRUCTURE_VIEWS, key="clf_view", default=STRUCTURE_VIEWS[0],
        on_change=_keep_view, label_visibility="collapsed",
    )
    view = st.session_state.get("clf_view") or STRUCTURE_VIEWS[0]

    if view == STRUCTURE_VIEWS[1]:
        with st.spinner("Embedding 3D conformer…"):
            molblock = molecular_properties.structure_3d(smiles)
        if molblock:
            theme.molecule_3d(molblock, height=430, caption="computed conformer · RDKit")
            return
        theme.note("3D embedding failed for this structure — showing the 2D depiction.")

    st.image(draw_structure_highlight(result, size=(640, 500)), width="stretch")


def _legend(thr: Thresholds) -> None:
    labels = dict(LEGEND_LABELS, acyl=chain_label(thr))
    theme.chips([(labels[k], _rgb(HIGHLIGHT_COLORS[k])) for k in LEGEND_KEYS])


def _features(result, thr: Thresholds) -> None:
    rows = feature_rows_en(result.features, thr)
    theme.section("Features", aux=f"{len(rows)} rows")
    theme.stat_grid(rows, wide=True, max_cols=5)

    bar_html, scale = _window_bar_html(result.features.mw, float(thr.mw_min), float(thr.mw_max))
    theme.subsection("Detection window", aux=f"MW · {scale}")
    st.markdown(bar_html, unsafe_allow_html=True)


def _why(result, thr: Thresholds) -> None:
    reasons = english_reasons(result, thr)
    theme.section("Why", aux=f"{len(reasons)} rules")
    theme.reasons(reasons)


def _identity_and_adducts(props: molecular_properties.MolecularProperties) -> None:
    """Three cards across: the long identifiers, then one card per polarity.

    The identities used to span the full width with the two adduct tables half of it
    each; on the full-width canvas that left the ion / m/z / Δ triplet stretched over
    ~600 px of empty row. Identity, ESI+ and ESI− side by side keeps every row the
    width it was designed for and makes the section a single even band.
    """
    theme.section("Identity & adducts")

    mass = props.exact_mass
    col_id, col_pos, col_neg = st.columns([1.15, 1, 1], gap="medium")

    with col_id:
        theme.copy_fields(
            [
                ("InChI Key", props.inchikey or "—", bool(props.inchikey)),
                ("Formula", props.formula or "—", bool(props.formula)),
                ("Canonical SMILES", props.canonical_smiles or "—", bool(props.canonical_smiles)),
                ("InChI", props.inchi or "—", bool(props.inchi)),
            ],
            heading="Molecular identity",
        )

    if props.formal_charge:
        z = props.formal_charge
        mz = None if mass is None else mass / abs(z)
        with col_pos:
            theme.copy_fields(
                [(_ion_label(z), "—" if mz is None else f"{mz:.4f}", mz is not None,
                  f"z = {z:+d}")],
                heading="Charged input", variant="adducts",
            )
        with col_neg:
            theme.note(f"The input carries a net charge of {z:+d}, so it is observed as is: "
                       "m/z = exact mass / |z|. Adducts assume a neutral molecule — enter "
                       "the neutral form to list them.")
    else:
        with col_pos:
            theme.copy_fields(_adduct_fields(mass, positive=True),
                              heading="ESI+ · positive mode", variant="adducts")
        with col_neg:
            theme.copy_fields(_adduct_fields(mass, positive=False),
                              heading="ESI− · negative mode", variant="adducts")

    if not props.formal_charge:
        m_txt = "n/a" if mass is None else f"{mass:.4f} Da"
        theme.note(f"m/z = (M + Δ) / |z|, where M is the neutral monoisotopic mass ({m_txt}) "
                   "and Δ includes the electron mass.")


def _descriptors(props: molecular_properties.MolecularProperties) -> None:
    rows = molecular_properties.property_rows(props)
    theme.section("Molecular properties", aux=f"{len(rows)} descriptors")
    theme.stat_grid(rows, wide=True, max_cols=5)


def _query_panel(thr: Thresholds, smiles: str, result) -> None:
    """Full-width query panel: input beside its live structure preview, then the
    threshold editor.

    The tab used to put the query in a one-third rail, which left a column of empty
    page for the ~1 500 px the results run below it, and squeezed the results into a
    ~750 px canvas. Input and preview now sit side by side across the top, which is
    also how the other three tabs group their inputs.

    Everything here renders on every run — including for an empty or invalid SMILES —
    so tuned thresholds and the example picker survive a bad input (Streamlit drops
    the state of widgets that are not drawn).
    """
    with theme.panel("clf-query"):
        col_in, col_fig = st.columns([1.4, 1], gap="large")

        with col_in:
            theme.subsection("Input")
            st.selectbox("Examples", list(EXAMPLES), key="clf_example",
                         on_change=_on_example_change)
            st.text_area("SMILES", key="clf_smiles", on_change=_on_smiles_change,
                         height=159,  # 6 rows: 6*20 + 39 — see "Peak-list heights" in the README
                         placeholder="e.g. CCCCCCCCCCCCCCCC(=O)O")
            theme.note("Updates as you edit — press ⌘↵ or click away to commit.")

        with col_fig:
            if result is not None:
                _structure_figure(smiles, result)
                _legend(thr)
            else:
                theme.subsection("Structure", aux="RDKit")
                theme.note("The 2D depiction with its substructure highlights appears here "
                           "as soon as the SMILES parses.")

        _threshold_inputs(columns=3)
        tuned = _tuned_thresholds(thr)
        if tuned:
            theme.note("Tuned: " + " · ".join(tuned))


def render():
    theme.tab_header(
        "Basic",
        subtitle="Lipidome/metabolome classification, chemical identity, descriptors and "
                 "adduct m/z from a SMILES string.",
        kicker="SMILES classifier",
    )

    if "clf_smiles" not in st.session_state:
        st.session_state["clf_smiles"] = EXAMPLES[DEFAULT_EXAMPLE]
        st.session_state.setdefault("clf_example", DEFAULT_EXAMPLE)

    # Read the thresholds before drawing the editor: the structure preview in the
    # panel is labelled with the current acyl-chain threshold.
    thr = _thresholds_from_state()
    smiles = st.session_state["clf_smiles"].strip()
    props = molecular_properties.from_smiles(smiles) if smiles else None
    valid = bool(props and props.valid)
    result = analyze(smiles, thr) if valid else None

    _query_panel(thr, smiles, result)

    problems = _threshold_problems(thr)
    if problems:
        st.warning("Inconsistent thresholds:\n\n" + "\n".join(f"- {p}" for p in problems))

    if not smiles:
        st.info("Enter a SMILES string or choose an example in the query panel.")
        _history_panel(smiles)
        return

    if not valid:
        st.error("RDKit could not parse this SMILES. "
                 + (molecular_properties.parse_error(smiles) or ""))
        _history_panel(smiles)
        return

    if props.n_fragments > 1:
        st.info(f"The input has {props.n_fragments} disconnected components (salt or "
                "mixture): identifiers, masses, descriptors and the classification "
                "refer to the whole input.")

    # The headline first: the classification is what this tab is asked for.
    en = category_en(result.category)
    theme.verdict(en, category_subtitle_en(en), theme.CATEGORY_ACCENT.get(en, theme.BLUE))
    _remember(smiles, _display_name(smiles, props), en)

    _features(result, thr)
    _why(result, thr)
    _identity_and_adducts(props)
    _descriptors(props)
    _history_panel(smiles)

    theme.divider()
    theme.note("Heuristic tool — lipid vs metabolite has no hard boundary. "
               "Tune the classifier thresholds in the query panel.")
