import math
from typing import List, Optional, Tuple

import streamlit as st

from mzdesk.shared import theme
from mzdesk.platform_fit import molecular_properties
from mzdesk.platform_fit.lipid_metabolite_classifier import analyze, Thresholds
from mzdesk.platform_fit.lipid_metabolite_viz import draw_structure_highlight, HIGHLIGHT_COLORS
from mzdesk.platform_fit.display_en import (category_en, category_subtitle_en, chain_label,
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
LEGEND_LABELS = {
    "isoprenoid": "Isoprenoid / sterol",
    "head": "Polar head group",
    "glycosyl": "Glycosyl ring",
}
OUTSIDE_RED = theme.CATEGORY_ACCENT["Invalid"]

THRESHOLD_ROWS = [
    [
        ("clf_mw_min", "mw_min", "MW min (Da)", 0.0, 2000.0, 10.0,
         "Below this MW, classify as Neither."),
        ("clf_mw_max", "mw_max", "MW max (Da)", 100.0, 10000.0, 100.0,
         "Above this MW, classify as Neither."),
        ("clf_logp_polar_max", "logp_polar_max", "LogP polar max", -10.0, 20.0, 0.1,
         "Without a lipid backbone: below → Metabolome; otherwise → Both."),
    ],
    [
        ("clf_acyl_min_carbons", "acyl_min_carbons", "Min acyl carbons", 4, 30, 1,
         "Minimum mostly hydrocarbon chain length for a lipid tail."),
        ("clf_medium_chain_max", "medium_chain_max", "Max medium-chain C", 4, 30, 1,
         "Maximum total carbons for the medium-chain Both rule."),
        ("clf_isoprene_min_units", "isoprene_min_units", "Min isoprene units", 1, 20, 1,
         "Minimum isoprene units for an isoprenoid chain."),
    ],
    [
        ("clf_boundary_logp_low", "boundary_logp_low", "Boundary logP min", -10.0, 20.0, 0.1,
         "Minimum logP for the medium-chain Both rule."),
        ("clf_boundary_logp_high", "boundary_logp_high", "Boundary logP max", -10.0, 20.0, 0.1,
         "Maximum logP for the medium-chain Both rule."),
        ("clf_boundary_mw_max", "boundary_mw_max", "Boundary MW max (Da)", 50.0, 5000.0, 10.0,
         "Maximum MW for the medium-chain Both rule."),
    ],
    [
        ("clf_chain_max_polar_fraction", "chain_max_polar_fraction",
         "Max chain polar fraction", 0.0, 1.0, 0.05,
         "Maximum heteroatom-bearing fraction of chain carbons."),
        ("clf_lipid_min_mw", "lipid_min_mw", "Min lipid MW (Da)", 0.0, 2000.0, 10.0,
         "Minimum MW for chain-only lipid evidence; backbone matches ignore it."),
        ("clf_polar_lipid_logp_max", "polar_lipid_logp_max", "Polar lipid logP max",
         -10.0, 20.0, 0.1,
         "Tail-free steroids at or below this logP are Both."),
    ],
]

STRUCTURE_VIEWS = ["2D", "3D"]


def _rgb(color):
    r, g, b = color
    return f"rgb({int(r * 255)},{int(g * 255)},{int(b * 255)})"


def _on_example_change():
    selected = st.session_state.get("clf_example")
    if selected and selected != PLACEHOLDER:
        st.session_state["clf_smiles"] = EXAMPLES[selected]


def _on_smiles_change():
    selected = st.session_state.get("clf_example", PLACEHOLDER)
    if (selected != PLACEHOLDER and
            st.session_state.get("clf_smiles", "").strip() != EXAMPLES.get(selected)):
        st.session_state["clf_example"] = PLACEHOLDER


def _reset_thresholds():
    for row in THRESHOLD_ROWS:
        for key, *_ in row:
            st.session_state.pop(key, None)


def _threshold_inputs(columns: int = 3) -> Thresholds:
    defaults = Thresholds()
    values = {}
    with st.expander("Thresholds"):
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
        st.button("Reset", key="clf_thr_reset", on_click=_reset_thresholds, width="stretch")
    return Thresholds(**values)


def _thresholds_from_state() -> Thresholds:
    defaults = Thresholds()
    values = {}
    for row in THRESHOLD_ROWS:
        for key, field, *_ in row:
            values[field] = st.session_state.get(key, getattr(defaults, field))
    return Thresholds(**values)


def _remember(smiles: str, name: str, category: str) -> None:
    history = st.session_state.setdefault("clf_history", [])
    history[:] = [item for item in history if item["smiles"] != smiles]
    history.insert(0, {"smiles": smiles, "name": name, "category": category})
    del history[8:]


def _recall(smiles: str) -> None:
    st.session_state["clf_smiles"] = smiles
    st.session_state["clf_example"] = PLACEHOLDER


def _history_panel(current: str) -> None:
    history = st.session_state.get("clf_history", [])
    if len(history) < 2:
        return
    theme.section("History")
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
    problems = []
    if thr.mw_min >= thr.mw_max:
        problems.append("MW min must be below MW max.")
    if thr.medium_chain_max < thr.acyl_min_carbons:
        problems.append("Medium-chain max must be ≥ min acyl carbons.")
    if thr.boundary_logp_low > thr.boundary_logp_high:
        problems.append("Boundary logP min must be ≤ max.")
    if thr.boundary_mw_max < thr.mw_min:
        problems.append("Boundary MW max must be ≥ MW min.")
    if thr.lipid_min_mw > thr.mw_max:
        problems.append("Lipid min MW must be ≤ MW max.")
    if thr.chain_max_polar_fraction >= 1.0:
        problems.append("Chain polar fraction must be below 1.")
    return problems


def _adduct_fields(mass: Optional[float], positive: bool) -> List[tuple]:
    rows = []
    for label, (delta, charge) in molecular_properties.ADDUCTS.items():
        if (charge > 0) != positive:
            continue
        mz = molecular_properties.adduct_exact_mz(mass, label)
        rows.append((label, "—" if mz is None else f"{mz:.4f}", mz is not None,
                     f"Δ {delta:+.4f}"))
    return rows


def _ion_label(charge: int) -> str:
    n = abs(charge)
    return f"[M]{n if n > 1 else ''}{'+' if charge > 0 else '-'}"


def _log_ticks(lo: float, hi: float) -> List[float]:
    ticks = [(m, m * 10.0 ** e)
             for e in range(math.floor(math.log10(lo)), math.ceil(math.log10(hi)) + 1)
             for m in (1, 2, 5)]
    ticks = [(m, value) for m, value in ticks if lo < value < hi]
    if len(ticks) > 6:
        ticks = [(m, value) for m, value in ticks if m == 1]
    return [value for _, value in ticks]


def _window_bar_html(mw: float, lo: float, hi: float) -> Tuple[str, str]:
    use_log = 0 < lo < hi

    def pos(value: float) -> float:
        if hi <= lo:
            return 0.0
        if use_log:
            t = ((math.log10(max(value, 1e-9)) - math.log10(lo)) /
                 (math.log10(hi) - math.log10(lo)))
        else:
            t = (value - lo) / (hi - lo)
        return 100.0 * min(1.0, max(0.0, t))

    inside = lo <= mw <= hi
    color = theme.BLUE if inside else OUTSIDE_RED
    status = "" if inside else (" · below" if mw < lo else " · above")
    p = pos(mw)
    ticks = _log_ticks(lo, hi) if use_log else []
    tick_marks = "".join(
        f'<div style="position:absolute;top:0;bottom:0;left:{pos(t):.2f}%;width:1px;'
        f'background:{theme.MUTED};opacity:.35;"></div>' for t in ticks)
    tick_labels = "".join(
        f'<span style="position:absolute;left:{pos(t):.2f}%;transform:translateX(-50%);">{t:g}</span>'
        for t in ticks if 8.0 <= pos(t) <= 92.0)
    mono = "ui-monospace,'SF Mono',Menlo,monospace"
    html = (
        '<div style="position:relative;padding-top:30px;margin:0 6px 4px 0;">'
        f'<div style="position:absolute;top:0;left:{p:.2f}%;transform:translateX(-{p:.2f}%);'
        f'white-space:nowrap;background:{color};color:#fff;border-radius:4px;'
        f'padding:2px 6px;font:600 .72rem/1.4 {mono};">MW {mw:.1f} Da{status}</div>'
        f'<div style="position:relative;height:16px;border:1px solid {theme.LINE};'
        f'background:{theme.SURFACE_2};border-radius:4px;">{tick_marks}'
        f'<div style="position:absolute;top:-4px;width:9px;height:22px;'
        f'left:clamp(-2px,calc({p:.2f}% - 4.5px),calc(100% - 7px));'
        f'border:1px solid {color};background:{color};border-radius:3px;"></div></div>'
        f'<div style="position:relative;height:1.2em;margin-top:6px;font-size:.7rem;'
        f'color:{theme.MUTED};"><span style="position:absolute;left:0;">{lo:g}</span>'
        f'{tick_labels}<span style="position:absolute;right:0;">{hi:g}</span></div></div>'
    )
    return html, ("log" if use_log else "linear")


def _keep_view():
    if st.session_state.get("clf_view") is None:
        st.session_state["clf_view"] = st.session_state.get("_clf_view_last", STRUCTURE_VIEWS[0])
    else:
        st.session_state["_clf_view_last"] = st.session_state["clf_view"]


def _display_name(smiles: str, props) -> str:
    for name, example_smiles in EXAMPLES.items():
        if example_smiles and example_smiles == smiles:
            return name.split(" (")[0]
    if props.formula:
        return props.formula
    return smiles if len(smiles) <= 18 else smiles[:17] + "…"


def _structure_figure(smiles: str, result) -> None:
    theme.subsection("Structure")
    st.segmented_control(
        "View", STRUCTURE_VIEWS, key="clf_view", default=STRUCTURE_VIEWS[0],
        on_change=_keep_view, label_visibility="collapsed",
    )
    view = st.session_state.get("clf_view") or STRUCTURE_VIEWS[0]

    if view == STRUCTURE_VIEWS[1]:
        with st.spinner("Building 3D…"):
            molblock = molecular_properties.structure_3d(smiles)
        if molblock:
            theme.molecule_3d(molblock, height=430)
            return
        st.warning("3D unavailable; showing 2D.")

    st.image(draw_structure_highlight(result, size=(640, 500)), width="stretch")


def _legend(thr: Thresholds) -> None:
    labels = dict(LEGEND_LABELS, acyl=chain_label(thr))
    theme.chips([(labels[key], _rgb(HIGHLIGHT_COLORS[key])) for key in LEGEND_KEYS])


def _features(result, thr: Thresholds) -> None:
    rows = feature_rows_en(result.features, thr)
    theme.section("Features")
    theme.stat_grid(rows, wide=True, max_cols=5)

    bar_html, scale = _window_bar_html(result.features.mw, float(thr.mw_min), float(thr.mw_max))
    theme.subsection("MW window", aux=scale)
    st.markdown(bar_html, unsafe_allow_html=True)


def _why(result, thr: Thresholds) -> None:
    theme.section("Why")
    theme.reasons(english_reasons(result, thr))


def _identity_and_adducts(props: molecular_properties.MolecularProperties) -> None:
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
            heading="Identity",
        )

    if props.formal_charge:
        charge = props.formal_charge
        mz = None if mass is None else mass / abs(charge)
        with col_pos:
            theme.copy_fields(
                [(_ion_label(charge), "—" if mz is None else f"{mz:.4f}", mz is not None,
                  f"z = {charge:+d}")],
                heading="Charged input", variant="adducts",
            )
        with col_neg:
            theme.note(f"Net charge {charge:+d}: m/z = exact mass / |z|. "
                       "Enter the neutral form for adducts.")
    else:
        with col_pos:
            theme.copy_fields(_adduct_fields(mass, positive=True),
                              heading="ESI+", variant="adducts")
        with col_neg:
            theme.copy_fields(_adduct_fields(mass, positive=False),
                              heading="ESI−", variant="adducts")

    if not props.formal_charge:
        mass_text = "n/a" if mass is None else f"{mass:.4f} Da"
        theme.note(f"m/z = (M + Δ) / |z|; M = {mass_text}; Δ includes electron mass.")


def _descriptors(props: molecular_properties.MolecularProperties) -> None:
    theme.section("Properties")
    theme.stat_grid(molecular_properties.property_rows(props), wide=True, max_cols=5)


def _query_panel(smiles: str, result) -> None:
    with theme.panel("clf-query"):
        col_in, col_fig = st.columns([1.4, 1], gap="large")
        with col_in:
            st.selectbox("Example", list(EXAMPLES), key="clf_example",
                         on_change=_on_example_change)
            st.text_area("SMILES", key="clf_smiles", on_change=_on_smiles_change,
                         height=119, placeholder="e.g. CCCCCCCCCCCCCCCC(=O)O")
        with col_fig:
            if result is not None:
                _structure_figure(smiles, result)
                _legend(_thresholds_from_state())
            else:
                theme.subsection("Structure")
        _threshold_inputs(columns=3)


def render():
    if "clf_smiles" not in st.session_state:
        st.session_state["clf_smiles"] = EXAMPLES[DEFAULT_EXAMPLE]
        st.session_state.setdefault("clf_example", DEFAULT_EXAMPLE)

    thr = _thresholds_from_state()
    smiles = st.session_state["clf_smiles"].strip()
    props = molecular_properties.from_smiles(smiles) if smiles else None
    valid = bool(props and props.valid)
    result = analyze(smiles, thr) if valid else None

    _query_panel(smiles, result)

    problems = _threshold_problems(thr)
    if problems:
        st.warning("Check thresholds:\n\n" + "\n".join(f"- {p}" for p in problems))

    if not smiles:
        st.info("Enter SMILES or choose an example.")
        _history_panel(smiles)
        return

    if not valid:
        st.error("Invalid SMILES. " + (molecular_properties.parse_error(smiles) or ""))
        _history_panel(smiles)
        return

    if props.n_fragments > 1:
        st.info(f"{props.n_fragments} disconnected components; results use the whole input.")

    english_category = category_en(result.category)
    theme.verdict(english_category, category_subtitle_en(english_category),
                  theme.CATEGORY_ACCENT.get(english_category, theme.BLUE))
    _remember(smiles, _display_name(smiles, props), english_category)

    _features(result, thr)
    _why(result, thr)
    _identity_and_adducts(props)
    _descriptors(props)
    _history_panel(smiles)
    theme.note("Heuristic classification; adjust thresholds as needed.")
