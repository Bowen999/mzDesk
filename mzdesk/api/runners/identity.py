"""身份模块：名称 / 数据库 ID / InChIKey 等 -> 结构与各库标识符（mzdesk.identity.resolve，联网）。"""
from mzdesk.api import common
from mzdesk.api.common import MISSING, InputError

ID_ROWS = ["HMDB", "KEGG", "ChEBI", "LIPID MAPS", "DrugBank", "CompTox", "ChEMBL"]


def run_identity(inputs):
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    from mzdesk.identity import resolve

    raw = inputs.get("structure") or {}
    query = str(raw.get("value") or "").strip() if isinstance(raw, dict) else str(raw).strip()
    if not query:
        raise InputError("Compound is required.")
    try:
        r = resolve.resolve(query, inputs.get("kind") or "Auto")
    except resolve.ResolveError as exc:
        raise InputError(str(exc)) from None

    formula = MISSING
    mol = Chem.MolFromSmiles(r.smiles) if r.smiles else None
    if mol is not None:
        formula = rdMolDescriptors.CalcMolFormula(mol)

    def row(label, values):
        """[label, value(s) joined by "; ", link cell or MISSING]"""
        if values is None or isinstance(values, str):
            values = [values] if values else []
        text = "; ".join(str(v) for v in values)
        href = r.links.get("PubChem" if label == "PubChem CID" else label)  # resolver keys the link "PubChem"
        return [label, text or MISSING, common.link(text, href) if href else MISSING]

    cid = str(r.cid) if r.cid is not None else None
    rows = [row("Name", r.name), row("IUPAC", r.iupac), row("CAS", r.ids.get("CAS")), row("SMILES", r.smiles),
            row("InChI", r.inchi), row("InChIKey", r.inchikey), row("PubChem CID", cid)]
    rows += [row(label, r.ids.get(label)) for label in ID_ROWS]

    warnings = list(r.warnings)
    if r.stereo_fallback:
        warnings.append("No exact stereoisomer match; the result is the flat (stereo-agnostic) structure.")
    outputs = {
        "resolved": {"values": {"Input read as": r.kind, "Name": r.name or MISSING,
                                "PubChem CID": cid or MISSING, "Formula": formula}},
        "identifiers": {"rows": rows},
        "candidates": {"rows": [[str(cid), title] for cid, title in r.candidates[1:]]},
        "synonyms": {"rows": [[s] for s in r.synonyms]},
    }
    svg = common.depiction_svg(r.smiles) if r.smiles else None
    if svg:
        outputs["depiction"] = {"svg": svg}
    return {"outputs": outputs, "warnings": warnings}


RUNNERS = {"identity": run_identity}
