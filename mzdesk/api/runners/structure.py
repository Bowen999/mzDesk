"""结构类模块：结构与格式、分子属性（RDKit）、平台适配里的可检测性分类（脂质组 / 代谢组启发式）。"""
import re

from rdkit import Chem

from mzdesk.api import blocks as bl
from mzdesk.api import common
from mzdesk.api.common import EngineError, InputError, fmt


def _check_parts(p, warnings, noun):
    if p.n_fragments > 1:
        warnings.append(f"{p.n_fragments} disconnected parts (salt or mixture): the {noun} the whole input.")
    if p.formal_charge:
        warnings.append(f"Net formal charge {p.formal_charge:+d}: the input is already an ion.")


def run_structure_formats(inputs):
    """Up to three compounds: one square 2D depiction + one square 3D viewer each, plus every format per compound."""
    from mzdesk.platform_fit import molecular_properties as mp

    warnings, compounds, sets = [], [], []
    for n, key in ((1, "structure"), (2, "structure_2"), (3, "structure_3")):
        value = inputs.get(key) if isinstance(inputs.get(key), dict) else {}
        # the optional label input of the structure field; without it the searched name, then "Compound n"
        given = " ".join(str(value.get("label") or "").split())
        label = given or f"Compound {n}"
        parsed = common.structure(inputs, key, label, required=(n == 1))
        if parsed is None:
            continue
        smiles = parsed["smiles"]
        p = mp.from_smiles(smiles)
        if not p.valid:
            raise InputError(f"{label}: RDKit could not parse it.")
        local = []
        _check_parts(p, local, "formats cover")
        mol3d = mp.structure_3d(smiles)
        if mol3d is None:
            local.append("3D embedding failed for this structure; only the 2D depiction is shown.")
        typed = str(value.get("value") or "").strip() if value.get("mode") == "Name / ID / Formula" else ""
        name = given or typed or label
        mol = Chem.MolFromSmiles(smiles)
        texts = {"SMILES": p.canonical_smiles, "InChI": p.inchi, "InChIKey": p.inchikey,
                 "MOL block": Chem.MolToMolBlock(mol) if mol is not None else ""}
        compounds.append({"label": name, "formula": p.formula or "",
                          "svg": common.depiction_svg(smiles, size=(420, 420)) or "", "mol3d": mol3d or ""})
        sets.append({"label": name, "formats": {k: v for k, v in texts.items() if v}})
        warnings.extend(f"{name}: {w}" for w in local)
    if not compounds:
        raise InputError("Enter at least one compound.")
    return {"outputs": {"compare": {"compounds": compounds}, "formats": {"sets": sets}}, "warnings": warnings}


def run_properties(inputs):
    from mzdesk.platform_fit import functional_groups as fg
    from mzdesk.platform_fit import molecular_properties as mp

    warnings = []
    smiles = common.structure(inputs)["smiles"]
    p = mp.from_smiles(smiles)
    if not p.valid:
        raise InputError("Structure: RDKit could not parse it.")
    _check_parts(p, warnings, "properties are for")

    def count(v):
        return common.MISSING if v is None else str(v)

    groups = [("Polarity", [("logP (Crippen)", fmt(p.logp, 2)), ("TPSA (Å²)", fmt(p.tpsa, 1)),
                            ("H-bond donors", count(p.hbd)), ("H-bond acceptors", count(p.hba))]),
              ("Size & topology", [("Heavy atoms", count(p.heavy_atoms)), ("Rotatable bonds", count(p.rotatable_bonds)),
                                   ("Rings", count(p.rings)), ("Aromatic rings", count(p.aromatic_rings))])]
    rows = [[group if i == 0 else "", name, value] for group, items in groups for i, (name, value) in enumerate(items)]
    views = fg.group_views(smiles)
    outputs = {
        "summary": {"values": {"Formula": p.formula or common.MISSING, "Monoisotopic mass": fmt(p.exact_mass, 4),
                               "Average MW": fmt(p.mw, 4), "Net charge": f"{p.formal_charge:+d}" if p.formal_charge else "0"}},
        "descriptors": {"rows": rows},
        "functional_groups": bl.groupmap_data(views["svg"], views["groups"]),
    }
    return {"outputs": outputs, "warnings": warnings}


def _cap(text):
    return text[:1].upper() + text[1:]


_SUBSCRIPT = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


def _formula_text(formula):
    """C21H52O6Si5 -> C₂₁H₅₂O₆Si₅ (unicode subscripts: no markup, so a table cell stays short)."""
    return re.sub(r"\d+", lambda m: m.group().translate(_SUBSCRIPT), formula)


def _reaction_options(results, *, detail, mol=None, groups=None):
    """Reaction-block options for the derivatization results that apply (n > 0). detail=True (Platform page): formula / mass / heavy mass and the EI ions,
    caveats and label note in the option's "?"; detail=False (Derivatization page: the table has them): the scheme only, the "?" explains the abbreviations.
    When mol and groups are given, each option carries its own reactive-site SVG; all four CIL channel options also carry a complete product SVG."""
    from mzdesk.platform_fit import derivatization as deriv

    options = []
    for r in results:
        if not r["n"]:
            continue
        steps = [bl.step(frm, reagent, to, n) for frm, reagent, to, n in r["steps"]]
        abbr = deriv.SCHEMES[r["name"]]["abbr"]
        label = deriv.reagent_label(r["reagent"])
        svg = product_svg = None
        if mol is not None and groups is not None:
            svg = deriv.highlight_svg(mol, deriv.reacting_groups([r], groups))
            product_svg = deriv.reaction_product_svg(mol, r, groups)
        if not detail:
            options.append(bl.reaction_option(label, steps, help=abbr, svg=svg, product_svg=product_svg))
            continue
        help_ = []
        if r["ions"]:
            help_.append("EI key ions: " + "; ".join(r["ions"]) + ".")
        help_ += [_cap(w) + "." for w in r["warns"]]
        if r["heavy"] is not None:
            help_.append(f"Heavy mass: ¹³C₂-labelled reagent, +{2 * deriv.C13_SHIFT:.4f} Da per label ({r['n']} label{'s' if r['n'] > 1 else ''}).")
        help_.append(abbr + ".")
        options.append(bl.reaction_option(label, steps, formula=_formula_text(r["formula"]), mass=round(r["mass"], 4),
                                          heavy=round(r["heavy"], 4) if r["heavy"] is not None else None, help=help_, svg=svg, product_svg=product_svg))
    return options


def run_platform_fit(inputs):
    """Platform 页：脂质组 / 代谢组分类、CIL LC-MS、GC-MS 四个 tile；CIL 或 GC 适用时再给一个官能团反应块（衍生化引擎）。"""
    from mzdesk.platform_fit import derivatization as deriv
    from mzdesk.platform_fit import display_en as en
    from mzdesk.platform_fit import lipid_metabolite_classifier as clf

    smiles = common.structure(inputs)["smiles"]
    res = clf.analyze(smiles)
    f = res.features
    if not f.valid:
        raise InputError("Structure: the classifier could not parse it.")
    category = en.category_en(res.category)
    lipid_help, metab_help = en.tile_reasons_en(res)
    lipid_on, metab_on = category in ("Lipidome", "Both"), category in ("Metabolome", "Both")
    tiles = [bl.tile("Lipidomics", "yes" if lipid_on else "no", value=en.lipid_class_en(f) if lipid_on and f.lipid_class else None, help=lipid_help),
             bl.tile("Metabolomics", "yes" if metab_on else "no", help=metab_help)]

    mol, neutralized = deriv.neutralize(Chem.MolFromSmiles(smiles))
    warnings = ["Input neutralized (charges removed) for the CIL / GC-MS rules."] if neutralized else []
    groups = deriv.find_groups(mol)
    cil = deriv.cil_assess(groups)
    gc = deriv.gc_assess(mol, groups, f.n_glycosyl)
    cil_help = [("CIL rule: labelable via " + "; ".join(cil["channels"]) + ".") if cil["cil"]
                else "CIL rule: no amine, phenol, hydroxyl, carboxyl or carbonyl group."]
    if category == "Lipidome":
        cil_help.append("Lipid: usually measured by lipidomics, not CIL.")
    tiles.append(bl.tile("CIL LC-MS", "yes" if cil["cil"] else "no", chips=cil["short"], help=cil_help))

    gc_option = deriv.gc_option(mol, groups)
    status = {"yes": "yes", "after derivatization": "partial", "no": "no"}[gc["gc"]]
    value = {"yes": "Volatile as is", "partial": "After derivatization"}.get(status) or _cap(gc["volatility"])
    gc_help = [x.replace(" -> ", " → ") for x in gc["reasons"]]
    if status != "no" and gc["stability"] != "OK":
        gc_help.append(f"Thermal stability: {gc['stability']}.")
    tiles.append(bl.tile("GC-MS", status, value=value, chips=[deriv.REAGENTS[gc_option][0]] if status == "partial" else None, help=gc_help))
    outputs = {"platform": bl.portal_data(tiles)}

    # why and how: one option per CIL channel that applies, plus the GC route (TMS, or MeOx + TMS when there is a carbonyl or a reducing-sugar ring)
    wanted = list(cil["options"]) + ([gc_option] if gc["gc"] != "no" else [])
    options = []
    if wanted:
        derived = deriv.derive(smiles, wanted)
        options = _reaction_options(derived["results"], detail=True, mol=derived["mol"], groups=derived["groups"])
        if gc["gc"] == "yes" and options and options[-1]["label"].startswith("GC"):
            options[-1]["help"] = ["Optional: the molecule is volatile enough to run directly."] + options[-1].get("help", [])
    if options:
        svg = deriv.highlight_svg(derived["mol"], deriv.reacting_groups(derived["results"], derived["groups"]))
        outputs["reactions"] = bl.reaction_data(options, svg=svg)
    else:             # nothing reacts: the block is still returned, with no options (the page hides an empty reaction block)
        outputs["reactions"] = {"options": []}
    return {"outputs": outputs, "warnings": warnings}


def run_derivatization(inputs):
    from mzdesk.platform_fit import derivatization as deriv

    smiles = common.structure(inputs)["smiles"]
    if Chem.MolFromSmiles(smiles) is None:
        raise InputError("Structure: RDKit could not parse it.")
    chosen = inputs.get("reagents")                    # missing = the default reagents; an empty list = none selected
    res = deriv.derive(smiles, chosen)
    warnings = ["Input neutralized (charges removed) before derivatization."] if res["neutralized"] else []
    shown = [r for r in res["results"] if r["n"]]
    if not res["results"]:
        warnings.append("No reagent selected.")
    elif not res["groups"]:
        warnings.append("No reactive group found in this structure.")
    elif not shown:
        warnings.append("None of the selected reagents reacts with the groups found in this structure.")

    tones = {"GC": "jp", "CIL": "report", "UV/FLD": "neutral"}         # badge colours by platform (blue / violet / grey)
    rows = []
    for r in shown:
        help_ = [f"Groups labelled: {r['text']}."]
        if r["ions"]:
            help_.append("EI key ions: " + "; ".join(r["ions"]) + ".")
        help_ += [_cap(w) + "." for w in r["warns"]]
        rows.append(bl.row([bl.badge(deriv.reagent_label(r["reagent"]), tones[r["platform"]]), str(r["n"]), _formula_text(r["formula"]), f"{r['mass']:.4f}",
                            f"{r['heavy']:.4f}" if r["heavy"] is not None else common.MISSING], help=help_))
    skipped = [deriv.reagent_label(r["reagent"]) for r in res["results"] if not r["n"]]
    outputs = {"deriv_rows": bl.table_data(rows, count=len(rows), help=("No reactive group for: " + ", ".join(skipped) + ".") if skipped else None)}
    options = _reaction_options(shown, detail=False, mol=res["mol"], groups=res["groups"])
    if options:
        svg = deriv.highlight_svg(res["mol"], deriv.reacting_groups(shown, res["groups"]))
        outputs["deriv_reactions"] = bl.reaction_data(options, svg=svg, help="Reactive groups found: " + ", ".join(deriv.group_text(res["groups"])) + ".")
    else:             # nothing reacts: the block is still returned, with no options (the page hides an empty reaction block; the warnings say why)
        outputs["deriv_reactions"] = {"options": []}
    return {"outputs": outputs, "warnings": warnings}


def run_separation(inputs):
    """Separation 页：MolGpKa 的位点 pKa → 各 pH 的电荷与 logD → 分离模式规则（规则见 separation_rules.py）。"""
    from mzdesk.platform_fit import pka_engine as pe
    from mzdesk.platform_fit import separation_rules as sr
    from mzdesk.platform_fit import separation_viz as viz

    smiles = common.structure(inputs)["smiles"]
    ph = common.number(inputs, "ph", "Mobile phase pH", default=2.7, minimum=sr.PH_MIN, maximum=sr.PH_MAX)
    try:
        a = pe.analyze(smiles)
    except pe.PkaError as exc:
        raise EngineError(str(exc)) from None
    except ValueError as exc:
        raise InputError(f"Structure: {exc}") from None
    st = a.at(ph)
    q, f = st["charge"], st["f_neutral"]

    if st["logd"] is None:
        logd_card = {"value": common.MISSING, "note": "permanent charge: no neutral form"}
    else:
        shift = st["logd"] - a.logp
        logd_card = {"value": sr.logd_text(st["logd"]),
                     "note": f"pH {ph:g}, logP {sr.logd_text(a.logp)}" + (f" ({sr.logd_text(shift)} from ionization)" if abs(shift) >= 0.05 else ", no ionization")}
    zwitterion = q == 0 and f < 0.5 * st["p_charge"]
    charge_note = f"{100 * st['p_charge']:.0f} % of molecules" + ("; mostly zwitterion" if zwitterion else "")
    if f >= 0.001:
        f_text = f"{100 * f:.1f} %"
    elif f > 0:
        f_text = "< 0.1 %"
    else:
        f_text = "0 %"

    sites = sorted(a.sites, key=lambda s: (s.pka, s.atom))
    rows = sr.mode_rows(ph=ph, logd=st["logd"], tpsa=a.tpsa, mw=a.mw, charge=q, f_neutral=f, permanent=a.permanent,
                        fixed_charge=a.fixed_charge, n_stereo=a.n_stereo, n_unassigned=a.n_unassigned, meso=a.meso)
    # details behind the "?" of the Notes block (the logD formula itself is the "?" of the logD card)
    hits = [f"{s.label} {why}: not basic, base site ignored (MolGpKa predicted pKa {s.pka:.1f})" for s, why in a.ignored]
    hits += sr.ph_notes(sites, ph, a.ph, a.f_neutral, a.logp, a.permanent)
    if zwitterion and st["logd"] is not None:
        hits.append(sr.ZWITTERION_NOTE)
    hits.append("Assumptions: sites are independent (Henderson–Hasselbalch per site, so successive ionizations of one group, such as phosphate pKa1 / pKa2, "
                "are not modelled); base sites on amide, sulfonamide and pyrrole-type aromatic N are dropped as not basic, but other weakly basic sites "
                "(aniline N, for example) still count; predicted pKa values carry errors of about one unit, sometimes more (validation/molgpka).")

    outputs = {
        "ph_summary": {"values": {"logD at pH": logd_card, "Dominant charge state": {"value": sr.charge_text(q), "note": charge_note},
                                  "Fraction neutral": {"value": f_text}}},
        "ph_curves": common.plotly(viz.ph_figure(a, ph, st)),
        "modes": bl.chips_data([bl.chip(label, level, help=reason) for label, level, reason in sr.mode_chips(rows)]),
        "sites": bl.table_data([[s.label, f"{s.pka:.1f}", s.kind] for s in sites], count=len(sites)),
        "rule_hits": {"items": hits},
    }
    svg = viz.sites_svg(a)
    if svg:
        outputs["sites_structure"] = {"svg": svg}
    return {"outputs": outputs, "warnings": a.notes}


def run_known_methods(inputs):
    """Known methods 页：同一化合物的已发表色谱方法（JP18、RepoRT），再加结构 / 性质最相近化合物的方法（引擎见 mzdesk/methods/known_methods.py）。"""
    from mzdesk.methods import known_methods as km

    smiles = common.structure(inputs)["smiles"]
    chosen = inputs.get("sources")
    chosen = list(km.SOURCES) if chosen is None else [chosen] if isinstance(chosen, str) else list(chosen)
    unknown = [s for s in chosen if s not in km.SOURCES]
    if unknown:
        raise InputError(f"Unknown source: {', '.join(map(str, unknown))}.")
    if not chosen:
        raise InputError("Choose at least one source.")
    try:
        res = km.search(smiles, chosen)
    except km.QueryError as exc:
        raise InputError(f"Structure: {exc}") from None
    except km.IndexProblem as exc:
        raise EngineError(str(exc)) from None
    return {"outputs": km.page_outputs(res), "warnings": res.warnings}


RUNNERS = {"structure_formats": run_structure_formats, "properties": run_properties, "platform_fit": run_platform_fit,
           "derivatization": run_derivatization, "separation": run_separation, "known_methods": run_known_methods}
