"""结构类模块：结构与格式、分子属性（RDKit）、平台适配里的可检测性分类（脂质组 / 代谢组启发式）。"""
from rdkit import Chem

from mzdesk.api import common
from mzdesk.api.common import InputError, fmt


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
    """只接了可检测性分类；pH、色谱柱、电离、加合物预测还是占位。"""
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
    svg = viz.highlight_svg(res)
    if svg:
        outputs["detect_structure"] = {"svg": svg}
    return {"outputs": outputs, "warnings": []}


RUNNERS = {"structure_formats": run_structure_formats, "properties": run_properties, "platform_fit": run_platform_fit}
