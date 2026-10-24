"""结构类模块：分子属性（RDKit）、平台适配里的可检测性分类（脂质组 / 代谢组启发式）。"""
from rdkit import Chem

from mzdesk.api import common
from mzdesk.api.common import InputError, fmt


def run_properties(inputs):
    from mzdesk.platform_fit import molecular_properties as mp

    warnings = []
    smiles = common.structure(inputs)["smiles"]
    p = mp.from_smiles(smiles)
    if not p.valid:
        raise InputError("Structure: RDKit could not parse it.")
    if p.n_fragments > 1:
        warnings.append(f"{p.n_fragments} disconnected parts (salt or mixture): the properties are for the whole input.")
    if p.formal_charge:
        warnings.append(f"Net formal charge {p.formal_charge:+d}: the input is already an ion.")
    mol3d = mp.structure_3d(smiles)
    if mol3d is None:
        warnings.append("3D embedding failed for this structure; only the 2D depiction is shown.")
    mol = Chem.MolFromSmiles(smiles)
    outputs = {
        "depiction": {"svg": common.depiction_svg(smiles)},
        "summary": {"values": {"Formula": p.formula or common.MISSING, "Exact mass": fmt(p.exact_mass, 4),
                               "logP": fmt(p.logp, 2), "Charge": f"{p.formal_charge:+d}" if p.formal_charge else "0"}},
        "descriptors": {"rows": [[k, v] for k, v in mp.property_rows(p)]},
        "formats": {"rows": [["SMILES (canonical)", p.canonical_smiles or common.MISSING], ["InChI", p.inchi or common.MISSING],
                             ["InChIKey", p.inchikey or common.MISSING]]},
        "molblock": {"text": Chem.MolToMolBlock(mol)},
    }
    if mol3d:
        outputs["viewer3d"] = {"mol3d": mol3d}
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
        "detect": {"values": {"Detectability": {"value": category, "note": en.category_subtitle_en(category)}}},
        "detect_features": {"rows": [[k, v] for k, v in en.feature_rows_en(res.features)]},
        "detect_reasons": {"items": en.english_reasons(res)},
    }
    svg = viz.highlight_svg(res)
    if svg:
        outputs["detect_structure"] = {"svg": svg}
    return {"outputs": outputs, "warnings": []}


RUNNERS = {"properties": run_properties, "platform_fit": run_platform_fit}
