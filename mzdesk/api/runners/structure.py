"""结构类模块：结构与格式、分子属性（RDKit）、平台适配里的可检测性分类（脂质组 / 代谢组启发式）。"""
from rdkit import Chem

from mzdesk.api import common
from mzdesk.api.common import EngineError, InputError, fmt


def _check_parts(p, warnings, noun):
    if p.n_fragments > 1:
        warnings.append(f"{p.n_fragments} disconnected parts (salt or mixture): the {noun} the whole input.")
    if p.formal_charge:
        warnings.append(f"Net formal charge {p.formal_charge:+d}: the input is already an ion.")


def run_structure_formats(inputs):
    from mzdesk.platform_fit import molecular_properties as mp

    warnings = []
    smiles = common.structure(inputs)["smiles"]
    p = mp.from_smiles(smiles)
    if not p.valid:
        raise InputError("Structure: RDKit could not parse it.")
    _check_parts(p, warnings, "formats cover")
    mol3d = mp.structure_3d(smiles)
    if mol3d is None:
        warnings.append("3D embedding failed for this structure; only the 2D depiction is shown.")
    texts = {"SMILES": p.canonical_smiles, "InChI": p.inchi, "InChIKey": p.inchikey,
             "MOL block": Chem.MolToMolBlock(Chem.MolFromSmiles(smiles))}
    outputs = {"depiction": {"svg": common.depiction_svg(smiles)},
               "formats": {"formats": {k: v for k, v in texts.items() if v}}}
    if mol3d:
        outputs["viewer3d"] = {"mol3d": mol3d}
    return {"outputs": outputs, "warnings": warnings}


def run_properties(inputs):
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
    outputs = {
        "summary": {"values": {"Formula": p.formula or common.MISSING, "Monoisotopic mass": fmt(p.exact_mass, 4),
                               "Average MW": fmt(p.mw, 4), "Net charge": f"{p.formal_charge:+d}" if p.formal_charge else "0"}},
        "descriptors": {"rows": rows},
    }
    return {"outputs": outputs, "warnings": warnings}


def run_platform_fit(inputs):
    """Platform 页：脂质组 / 代谢组分类，加上 CIL LC-MS 与 GC-MS 两节（衍生化引擎）。"""
    from mzdesk.platform_fit import derivatization as deriv
    from mzdesk.platform_fit import display_en as en
    from mzdesk.platform_fit import lipid_metabolite_classifier as clf
    from mzdesk.platform_fit import lipid_metabolite_viz as viz

    smiles = common.structure(inputs)["smiles"]
    res = clf.analyze(smiles)
    if not res.features.valid:
        raise InputError("Structure: the classifier could not parse it.")
    category = en.category_en(res.category)
    outputs = {
        "detect": {"values": {"Best platform": {"value": category, "note": en.category_subtitle_en(category)}}},
        "detect_features": {"rows": [[k, v] for k, v in en.feature_rows_en(res.features)]},
        "detect_reasons": {"items": en.english_reasons(res)},
    }
    mol, neutralized = deriv.neutralize(Chem.MolFromSmiles(smiles))
    warnings = ["Input neutralized (charges removed) for the CIL / GC-MS rules."] if neutralized else []
    groups = deriv.find_groups(mol)
    cil = deriv.cil_assess(groups)
    gc = deriv.gc_assess(mol, groups, res.features.n_glycosyl)
    lipid = category == "Lipidome"
    cil_note = "lipid: usually measured by lipidomics, not CIL" if lipid else ""
    outputs["cil"] = {"values": {"CIL LC-MS": {"value": "yes" if cil["cil"] else "no", "note": cil_note},
                                 "Labelable groups": {"value": ", ".join(cil["channels"]) or common.MISSING, "note": ""}}}
    outputs["gc"] = {"values": {"GC-MS": {"value": gc["gc"], "note": ""}, "Volatility": {"value": gc["volatility"], "note": ""},
                                "Thermal stability": {"value": gc["stability"], "note": ""}}}
    cil_rule = ("CIL rule: labelable via " + "; ".join(cil["channels"])) if cil["cil"] else "CIL rule: no amine, phenol, hydroxyl, carboxyl or carbonyl group"
    outputs["gc_cil_reasons"] = {"items": [cil_rule] + ([f"CIL rule: {cil_note}"] if lipid else []) + gc["reasons"]}
    svg = viz.highlight_svg(res)
    if svg:
        outputs["detect_structure"] = {"svg": svg}
    return {"outputs": outputs, "warnings": warnings}


def run_derivatization(inputs):
    from mzdesk.platform_fit import derivatization as deriv

    smiles = common.structure(inputs)["smiles"]
    if Chem.MolFromSmiles(smiles) is None:
        raise InputError("Structure: RDKit could not parse it.")
    res = deriv.derive(smiles, inputs.get("reagents") or [])
    warnings = ["Input neutralized (charges removed) before derivatization."] if res["neutralized"] else []
    if not res["groups"]:
        warnings.append("No reactive group found in this structure.")
    outputs = {"deriv_groups": {"rows": res["group_rows"]}, "deriv_rows": {"rows": res["rows"]},
               "deriv_structure": {"svg": res["svg"]}}
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
                     "note": f"logP {sr.logd_text(a.logp)}" + (f", {sr.logd_text(shift)} from ionization" if abs(shift) >= 0.05 else ", no ionization at this pH")}
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
    hits = [f"{s.label} {why}: not basic, base site ignored (MolGpKa predicted pKa {s.pka:.1f})" for s, why in a.ignored]
    hits += sr.ph_notes(sites, ph, a.ph, a.f_neutral, a.logp, a.permanent)
    if zwitterion and st["logd"] is not None:
        hits.append(sr.ZWITTERION_NOTE)
    hits.append("Assumptions: sites are independent (Henderson–Hasselbalch per site, so successive ionizations of one group, such as phosphate pKa1 / pKa2, "
                "are not modelled); base sites on amide, sulfonamide and pyrrole-type aromatic N are dropped as not basic, but other weakly basic sites "
                "(aniline N, for example) still count; predicted pKa values carry errors of about one unit, sometimes more (validation/molgpka).")
    hits.append("logD(pH) = logP (Crippen) + log10(fraction fully uncharged): ions are assumed not to partition and zwitterions count as charged, "
                "so logD is a lower bound at low or high pH.")

    outputs = {
        "ph_summary": {"values": {"logD at pH": logd_card, "Dominant charge state": {"value": sr.charge_text(q), "note": charge_note},
                                  "Fraction neutral": {"value": f_text, "note": "fully uncharged"}}},
        "sites": {"rows": [[s.label, f"{s.pka:.1f}", s.kind] for s in sites]},
        "ph_curves": common.plotly(viz.ph_figure(a, ph, st)),
        "modes": {"rows": rows},
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
