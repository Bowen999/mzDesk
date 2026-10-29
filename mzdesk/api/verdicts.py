"""The one-line answer shown above a module's results ("verdict"). Built from the finished outputs, so runners stay as they are.

build(module_id, inputs, outputs) -> str | None. The text may use <b> (the page renders <sub> <sup> <i> <b> <code>).
A runner may also return its own "verdict"; that one wins. None = the page shows no verdict line.
"""
import html
import re

MINUS = "−"


def _card(outputs, block, label):
    """Value of one card (a plain string or the "value" of {"value", "note"}); None when missing or blank."""
    v = ((outputs.get(block) or {}).get("values") or {}).get(label)
    if isinstance(v, dict):
        v = v.get("value")
    return None if v in (None, "", "—") else str(v)


def _note(outputs, block, label):
    v = ((outputs.get(block) or {}).get("values") or {}).get(label)
    return v.get("note") if isinstance(v, dict) else None


def _join(items, word="and"):
    items = [i for i in items if i]
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + f" {word} " + items[-1] if items else ""


def _b(text):
    return f"<b>{text}</b>"


def _rows(outputs, block):
    rows = (outputs.get(block) or {}).get("rows") or []
    return [r["cells"] if isinstance(r, dict) else r for r in rows]


def _text(cell):
    if isinstance(cell, dict):
        return cell.get("text") or cell.get("badge") or ""
    return "" if cell is None else str(cell)


# ---------------------------------------------------------------- one function per module
def identity(i, o):
    name, formula, cid = _card(o, "resolved", "Name"), _card(o, "resolved", "Formula"), _card(o, "resolved", "PubChem CID")
    if not name and not formula:
        return None
    parts = [_b(name) if name else None, formula, f"PubChem CID {cid}" if cid else None]
    return " · ".join(p for p in parts if p)


def properties(i, o):
    formula, mass = _card(o, "summary", "Formula"), _card(o, "summary", "Monoisotopic mass")
    if not formula:
        return None
    desc = {r[1]: r[2] for r in _rows(o, "descriptors") if len(r) >= 3}
    desc = {k: str(x).replace("-", MINUS) for k, x in desc.items()}
    extra = [f"logP {desc['logP (Crippen)']}" if desc.get("logP (Crippen)") not in (None, "—") else None,
             f"TPSA {desc['TPSA (Å²)']} Å²" if desc.get("TPSA (Å²)") not in (None, "—") else None]
    text = f"{_b(formula)}, monoisotopic mass {_b(mass + ' Da')}" if mass else _b(formula)
    extra = [e for e in extra if e]
    return text + ("; " + ", ".join(extra) if extra else "") + "."


def platform_fit(i, o):
    tiles = (o.get("platform") or {}).get("tiles") or []
    yes = [t["label"] for t in tiles if t.get("status") == "yes"]
    partial = [t for t in tiles if t.get("status") == "partial"]
    if not yes and not partial:
        return "None of the four workflows is expected to measure it."
    yes = [x.lower() if x in ("Lipidomics", "Metabolomics") else x for x in yes]
    text = ("Measurable by " + _join([_b(x) for x in yes]) + ".") if yes else "No workflow works as is."
    for t in partial:
        via = f" ({', '.join(t['chips'])})" if t.get("chips") else ""
        text += f" {t['label']} only {_b('after derivatization')}{via}."
    return text


def separation(i, o):
    chips = (o.get("modes") or {}).get("chips") or []
    good = [c["label"] for c in chips if c.get("level") == "good"]
    possible = [c["label"] for c in chips if c.get("level") == "possible"]
    ph = i.get("ph", 2.7)
    neutral, logd, charge = _card(o, "ph_summary", "Fraction neutral"), _card(o, "ph_summary", "logD at pH"), _card(o, "ph_summary", "Dominant charge state")
    if good:
        head = f"{_join([_b(g) for g in good], 'or')} {'fits' if len(good) == 1 else 'fit'} best at pH {ph:g}"
    elif possible:
        head = f"No clear best mode at pH {ph:g}; possible: {_join([_b(p) for p in possible], 'or')}"
    else:
        head = f"No separation mode stands out at pH {ph:g}"
    state = []
    try:
        frac = float(re.sub(r"[^\d.]", "", neutral or "")) if neutral and "<" not in neutral else 0.0
    except ValueError:
        frac = 0.0
    if frac >= 50:
        state.append(f"mostly neutral ({neutral})")
    elif charge is not None:
        state.append(f"net charge {charge}")
    if logd:
        state.append(f"logD {logd}")
    return head + (": " + ", ".join(state) if state else "") + "."


def ionization(i, o):
    chips = {c["label"]: c for c in (o.get("chips") or {}).get("chips") or []}
    pol = chips.get("Polarity")
    if not pol:
        return None
    v = (pol.get("value") or "").strip()
    if pol.get("level") == "poor":
        text = "Weak ESI response in both polarities."
    elif v in ("±", "+/−", "+/-", "both"):
        text = f"Works in {_b('both ESI+ and ESI−')}."
    elif v.startswith("+"):
        text = f"Best in {_b('ESI+')}" + (f"; ESI− {chips['ESI−'].get('value')}" if chips.get("ESI−", {}).get("value") else "") + "."
    else:
        text = f"Best in {_b('ESI−')}" + (f"; ESI+ {chips['ESI+'].get('value')}" if chips.get("ESI+", {}).get("value") else "") + "."
    ei, mi = chips.get("EI"), chips.get("M+•")
    if ei and ei.get("value"):
        text += f" EI {ei['value']}" + (f", molecular ion {mi['value']}" if mi and mi.get("value") else "") + "."
    return text


def adducts(i, o):
    ion, mz = _card(o, "base", "Dominant ion"), _card(o, "base", "m/z")
    if not ion:
        return None
    n = len(_rows(o, "ions"))
    likely = _note(o, "base", "Dominant ion")
    return (f"Expect {_b(ion)} at {_b('m/z ' + mz) if mz else ''} as the strongest ion" + (f" ({likely})" if likely else "")
            + (f"; {n} ions in total." if n > 1 else "."))


def non_ms(i, o):
    tiles = (o.get("detectors") or {}).get("tiles") or []
    yes = [t for t in tiles if t.get("status") == "yes"]
    partial = [t for t in tiles if t.get("status") == "partial"]
    if not yes and not partial:
        return "Not expected to be detectable without MS."
    def name(t):
        return _b(t["label"]) + (f" ({t['value']})" if t.get("value") else "")
    text = ("Detectable by " + _join([name(t) for t in yes]) + ".") if yes else "No detector works well."
    if partial:
        text += " Possible: " + _join([t["label"] + (f" ({t['value']})" if t.get("value") else "") for t in partial]) + "."
    return text


def derivatization(i, o):
    rows = _rows(o, "deriv_rows")
    if not rows:
        return "None of the selected reagents reacts with this structure."
    names = [_b(_text(r[0])) for r in rows]
    most = max(rows, key=lambda r: int(_text(r[1]) or 0))
    return f"{len(rows)} reagent{'s' if len(rows) > 1 else ''} react: {_join(names)}. Most labels: {_text(most[1])} ({_text(most[0])})."


def known_methods(i, o):
    match = (_card(o, "summary", "Match") or "").lower()
    n_exact = len(_rows(o, "exact"))
    groups = (o.get("similar") or {}).get("groups") or []
    n_sim = sum(len(g.get("rows") or []) for g in groups)
    sim = f"{n_sim} more from {len(groups)} similar compound{'s' if len(groups) != 1 else ''}" if groups else ""
    if match.startswith("exact") and n_exact:
        return f"{_b(f'{n_exact} published method' + ('s' if n_exact != 1 else ''))} for this exact compound" + (f", plus {sim}." if sim else ".")
    if groups:
        return f"No published method for the compound itself; {_b(sim)}."
    return "No published method found for it or for similar compounds."


def ms1_theory(i, o):
    mz, formula = _card(o, "theory", "Precursor m/z"), _card(o, "theory", "Ion formula")
    if not mz:
        return None
    adduct = _note(o, "theory", "Precursor m/z")
    top, top_mz = _card(o, "theory", "Most abundant peak"), _note(o, "theory", "Most abundant peak")
    n = len(_rows(o, "pattern"))
    text = f"{_b(adduct) + ' at ' if adduct else ''}{_b('m/z ' + mz)}" + (f" ({formula})" if formula else "")
    if top:
        text += f"; {n} isotope peak{'s' if n != 1 else ''}, most abundant {_b(top)}" + (f" ({top_mz})" if top_mz else "")
    return text + "."


def ms1_formula(i, o):
    top = _card(o, "summary", "Top formula")
    if not top:
        return "No formula fits the precursor within the tolerance."
    adduct = _card(o, "summary", "Adduct used")
    rows = _rows(o, "candidates")
    prob = _text(rows[0][4]) if rows and len(rows[0]) > 4 else None
    return (f"Top formula {_b(top)}" + (f" as {adduct}" if adduct else "") + (f", probability {prob}" if prob and prob != "—" else "")
            + (f"; {len(rows)} candidates." if len(rows) > 1 else "."))


def ms1_search(i, o):
    n = _card(o, "summary", "Compounds")
    if n is None:
        return None
    try:
        tol = f"{float(i.get('tol_ppm') or 5):g}"
    except (TypeError, ValueError):
        tol = "5"
    note_all = _note(o, "summary", "Compounds") or ""
    hid = re.search(r"(\d+) compounds? in (\d+) formula groups? hidden: isotope score below ([\d.]+)", note_all)
    if n == "0" and hid:
        return f"All {hid.group(1)} compounds within {tol} ppm fail the isotope score (≥ {hid.group(3)})."
    if n == "0":
        dbs = _card(o, "summary", "Databases searched")
        return f"Nothing matched within {tol} ppm" + (f" in {_join(dbs.split(', '))}." if dbs else ".")
    name, note = _card(o, "summary", "Closest hit"), _note(o, "summary", "Closest hit")
    formulas = _card(o, "summary", "Formulas")
    capped = "capped" in (_note(o, "summary", "Compounds") or "")
    text = _b(f"{n} compound{'s' if n != '1' else ''}") + (f" ({formulas} formula{'s' if formulas != '1' else ''})" if formulas else "") + f" within {tol} ppm"
    text += ", list capped" if capped else ""
    if hid:
        text += f", {hid.group(2)} formula group{'s' if hid.group(2) != '1' else ''} hidden by isotope score"
    if name and note:
        ion, _, ppm = note.partition(" · ")
        text += f"; closest: {_b(name)} {ion}" + (f" ({ppm})" if ppm else "")
    return text + "."


def ms2_clean(i, o):
    n_in, n_out = _card(o, "summary", "Peaks in"), _card(o, "summary", "Peaks out")
    if not n_in:
        return None
    removed = _note(o, "summary", "Removed")
    fmt = _card(o, "summary", "Detected format")
    return (f"Kept {_b(f'{n_out} of {n_in} peaks')}" + (f" ({removed} removed)" if removed and n_out != n_in else "")
            + (f"; read as {fmt}." if fmt else "."))


def ms2_ei(i, o):
    mi, base = _card(o, "summary", "Molecular ion m/z"), _card(o, "summary", "Base peak m/z")
    if not base:
        return None
    if mi and mi == base:
        return f"The molecular ion {_b('m/z ' + mi)} is the base peak."
    return f"Base peak {_b('m/z ' + base)}" + (f"; molecular ion at m/z {mi}." if mi else ".")


def ms2_esi(i, o):
    peaks, mz = _card(o, "summary", "Peaks"), _card(o, "summary", "Precursor m/z")
    if not peaks:
        return None
    engine, formula = _card(o, "summary", "Engine"), _card(o, "summary", "Ion formula")
    return (f"{_b(peaks + ' predicted peaks')}" + (f" for {formula}" if formula else "") + (f" at m/z {mz}" if mz else "")
            + (f" ({engine})." if engine else "."))


_SHORT_METRIC = {"SpectralEntropy (weighted)": "entropy", "SpectralEntropy (unweighted)": "unweighted entropy"}


def ms2_library(i, o):
    n, scored = _card(o, "summary", "Hits ≥ threshold"), _card(o, "summary", "Candidates scored")
    if n is None or scored is None:
        return None
    if scored == "0":
        if i.get("spectrum_type") == "EI":
            return "The selected libraries have no EI spectra."
        return "No library spectrum has a matching precursor m/z, so nothing was scored."
    m = re.search(r"≥ ([\d.\-]+)$", _note(o, "summary", "Hits ≥ threshold") or "")
    thr = m.group(1) if m else (f"{float(i.get('threshold') or 0.7):g}")
    if n == "0":
        return f"No hit: none of the {scored} candidates reaches {thr} on any selected metric."
    name, note = _card(o, "summary", "Best hit"), _note(o, "summary", "Best hit")
    score, metric = _card(o, "summary", "Rank score"), _note(o, "summary", "Rank score")
    parts = [f"Best hit {_b(name)}" + (f" ({note})" if note else "") if name else None,
             f"{_SHORT_METRIC.get(metric, metric)} {score}" if score and metric else None,
             f"{n} hit{'s' if n != '1' else ''} ≥ {thr}"]
    return " · ".join(p for p in parts if p) + "."


def ms2_similarity(i, o):
    score = _card(o, "summary", "Primary score")
    try:
        s = float(score)
    except (TypeError, ValueError):
        return None
    word = "Very similar" if s >= 0.8 else "Similar" if s >= 0.6 else "Weakly similar" if s >= 0.4 else "Not similar"
    passed, matched = _card(o, "summary", "Metrics ≥ threshold"), _card(o, "summary", "Matched peaks")
    extra = [f"{passed} metrics above threshold" if passed else None, f"{matched} matched peaks" if matched else None]
    return f"{_b(word)}: entropy similarity {_b(score)}" + ("; " + ", ".join(e for e in extra if e) if any(extra) else "") + "."


def ms2_fragments(i, o):
    peaks, inten = _card(o, "summary", "Explained peaks"), _card(o, "summary", "Explained intensity")
    if not peaks:
        return None
    return f"The structure explains {_b(peaks.replace(' / ', ' of ') + ' peaks')}" + (f", {_b(inten + ' %')} of the intensity." if inten else ".")


def ms2_structure(i, o):
    rows = _rows(o, "structures")
    formula = _card(o, "summary", "Top formula")
    if not rows:
        if formula:
            return f"No structure in MS-FINDER's local database for {_b(formula)} and the other top formulas; the formula table has the candidates."
        return "MS-FINDER found no molecular formula for this precursor."
    name, form, score = _text(rows[0][2]), _text(rows[0][4]), _text(rows[0][5])
    n = (_card(o, "summary", "Candidates") or "").split(" / ")[0]
    try:
        score = f"{float(score):.2f}"
    except ValueError:
        pass
    return (f"Top candidate: {_b(name)} ({form}), score {score}" + (f" · {n} candidates" if n.isdigit() and n != "1" else "") + ".")


def ms2_denovo(i, o):
    rows = _rows(o, "candidates")
    formula = _card(o, "summary", "Formula used")
    if not rows:
        return f"FRIGID generated no formula-matched structure for {_b(formula)} (see the warnings)." if formula else None
    smiles, score = _text(rows[0][2]), _text(rows[0][4])
    short = html.escape(smiles if len(smiles) <= 60 else smiles[:59] + "…")
    n = len(rows)
    return (f"Top generated candidate: {_b(short)}" + (f", similarity {score}" if score not in ("", "—") else "")
            + f" · {n} formula-matched candidate{'s' if n != 1 else ''} (hypotheses).")


def stability(i, o):
    t_half = _card(o, "summary", "t½")
    if not t_half:
        return None
    system = i.get("system") or "Liver microsomes"
    unit = "µL/min/10⁶ cells" if system == "Hepatocytes" else "µL/min/mg"
    clint, clh = _card(o, "summary", "CLint"), _card(o, "summary", "Predicted hepatic CL")
    parts = [f"t½ {_b(t_half + ' min')}"]
    if clint:
        parts.append(f"CLint {_b(clint + ' ' + unit)}")
    if clh:
        parts.append(f"predicted hepatic CL {_b(clh + ' mL/min/kg')}")
    text = " · ".join(parts)
    if not clint and system == "Plasma":
        text += " · plasma stability: t½ only"
    return text + "."


BUILDERS = {f.__name__: f for f in (identity, properties, platform_fit, separation, ionization, adducts, non_ms, derivatization, known_methods,
                                    ms1_theory, ms1_formula, ms1_search, ms2_clean, ms2_ei, ms2_esi, ms2_library, ms2_similarity, ms2_fragments,
                                    ms2_structure, ms2_denovo, stability)}
NO_VERDICT = {"structure_formats"}    # the page itself is the answer (drawings + formats)


def build(module_id, inputs, outputs):
    f = BUILDERS.get(module_id)
    if not f:
        return None
    try:
        return f(inputs or {}, outputs or {})
    except Exception:            # a verdict is a convenience: never fail a run because of it
        return None
