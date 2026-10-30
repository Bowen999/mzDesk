"""Alternative explanations 页（exclusivity）：三张表。

1. Same-formula candidates：PubChem 里分子式完全相同的结构（在线），按 InChIKey 前 14 位合并连接方式，可看库归属。
2. Isobaric interferences：同一 m/z 窗口里分子式不同的化合物（PubChem 精确质量检索，在线），每个分子式一行。
3. Isomers：MS2 能不能把化合物和它的同分子式异构体分开（MS-FINDER 本地库，离线；规则见 mzdesk/exclusivity/isomers.py）。

输入可以是结构、分子式，或只有 precursor m/z：只有 m/z 时跳过前两张表的分子式部分和异构体表。两张在线表失败只给 warning，
不影响异构体表；异构体表需要结构和 MS2。
"""
import html
import re

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from mzdesk.api import common
from mzdesk.api.blocks import badge, depicted, hinted, table_data
from mzdesk.api.common import EngineError, InputError, fmt, signed
from mzdesk.spec import ADDUCTS_NEG, ADDUCTS_POS, DATABASES

PEAKS_SHOWN = 4         # fragments written in a cell; all of them are in the cell's "?"
NAMES_IN_VERDICT = 3
NAME_CHARS = 40         # longer database names are cut in the verdict line (the table has them in full)
SMILES_DRAWING = (260, 182)   # 2D structure drawn on hover over a SMILES cell


def _smiles_cell(smiles):
    """The SMILES text with its 2D structure drawn on hover; plain text when there is nothing to draw."""
    return depicted(smiles, common.depiction_svg(smiles, size=SMILES_DRAWING)) if smiles else common.MISSING


def _peaks_cell(peaks):
    """'195.0451 (6.7 %) · …' with m/z, intensity, ion formula and fragment SMILES of every peak in the "?"."""
    if not peaks:
        return common.MISSING
    text = " · ".join(f"{p['mz']:.4f} ({p['rel']:.1f} %)" for p in peaks[:PEAKS_SHOWN])
    if len(peaks) > PEAKS_SHOWN:
        text += f" · +{len(peaks) - PEAKS_SHOWN}"
    return hinted(text, [f"m/z {p['mz']:.4f} · {p['rel']:.1f} % · {p['formula'] or '?'} · {p['smiles'] or '?'}" for p in peaks])


def _rows(res):
    by_id = {i["id"]: i for i in res["isomers"]}
    rows = []
    for c in res["candidates"]:
        explained = (f"{len(c['explained'])} / {res['peaks_considered']} · {c['explained_fraction'] * 100:.0f} %"
                     if c["evaluated"] else "not scored")
        name = c["name"] or common.MISSING
        if c["id"] == res["top_id"]:
            rows.append([name, _smiles_cell(c["smiles"]), explained, badge("your compound"), common.MISSING, common.MISSING])
            continue
        iso = by_id[c["id"]]
        verdict = (badge("not scored", help="MS-FINDER did not score this structure, so it cannot be ruled out.") if iso["distinguished"] is None
                   else badge("yes", "good") if iso["distinguished"] else badge("no", "warn"))
        rows.append([name, _smiles_cell(c["smiles"]), explained, verdict, _peaks_cell(iso["distinguishing_peaks"]), _peaks_cell(iso["isomer_only_peaks"])])
    return rows


def _same_rows(res):
    rows = []
    for g in res["rows"]:
        name = hinted(g["name"] or common.MISSING, g["names"]) if g["names"] else (g["name"] or common.MISSING)
        rows.append([g["key"], name, ", ".join(g["databases"]) or common.MISSING,
                     str(g["stereoisomers"]) if g["stereoisomers"] else "0"])
    return rows


def _isobaric_rows(res):
    rows = []
    for g in res["rows"]:
        notes = list(g["names"])
        if g["compounds"] > 1:
            notes.append(f"{g['compounds']} compounds with this formula; the closest is shown.")
        name = hinted(g["name"] or g["formula"], notes) if notes else (g["name"] or g["formula"])
        rows.append([g["formula"], name, ", ".join(g["databases"]) or common.MISSING, fmt(g["mz"], 4), signed(g["ppm"], 1)])
    return rows


def _short(name):
    return name if len(name) <= NAME_CHARS else name[:NAME_CHARS - 1] + "…"


def _verdict(res, compound):
    names = {c["id"]: html.escape(_short(c["name"] or c["smiles"])) for c in res["candidates"]}
    me = f"<b>{html.escape(compound)}</b>"
    n = len(res["isomers"])
    if res["verdict"] == "single":
        who = "the one isomer considered fails" if n == 1 else f"each of the {n} isomers considered fails"
        return f"{me} is the single answer: {who} to explain at least one of its fragments."
    if res["verdict"] == "indistinguishable":
        others = [names[i] for i in res["group"][1:]]
        shown = ", ".join(others[:NAMES_IN_VERDICT]) + (f" and {len(others) - NAMES_IN_VERDICT} more" if len(others) > NAMES_IN_VERDICT else "")
        return f"<b>These cannot be separated</b> by this spectrum: {me} and {shown} (of {n} isomer{'s' if n != 1 else ''} considered)."
    return f"<b>No isomer verdict</b>: {html.escape(res['notes'][0].rstrip('.'))}."


def _neutral_formula(smiles, warnings):
    """结构 -> 中性分子式（去掉电荷后缀）；RDKit 读不了时返回 None。"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    if len(Chem.GetMolFrags(mol)) > 1:
        warnings.append("The structure has several disconnected parts; the formula covers all of them.")
    charge = Chem.GetFormalCharge(mol)
    if charge:
        warnings.append(f"The structure carries a net charge ({charge:+d}); its neutral formula is used for the same-formula search.")
    return re.sub(r"[+-]\d*$", "", rdMolDescriptors.CalcMolFormula(mol))


def _databases(inputs):
    dbs = inputs.get("databases")
    dbs = ["HMDB", "PubChem"] if dbs is None else [dbs] if isinstance(dbs, str) else list(dbs)
    if not dbs:
        raise InputError("Choose at least one database.")
    unknown = [d for d in dbs if d not in DATABASES]
    if unknown:
        raise InputError(f"Unknown database: {', '.join(unknown)}.")
    return dbs


def run_exclusivity(inputs):
    from mzdesk.exclusivity import isomers
    from mzdesk.ms1 import exact_mass as em

    warnings = []
    struct = common.structure(inputs, label="Compound or formula", allow_formula=True, required=False, warnings=warnings)
    field_mz, positive, adduct = common.precursor(inputs)
    peaks, file_mz = common.spectrum(inputs, "ms2", "MS2 spectrum", warnings, required=False)
    precursor_mz = field_mz or file_mz
    if field_mz is None and file_mz:
        warnings.append(f"Precursor m/z {file_mz:.4f} is read from the spectrum.")
    elif field_mz and file_mz and abs(file_mz - field_mz) > 0.01:
        warnings.append(f"The MS2 file says precursor m/z {file_mz:.4f}; the form value {field_mz:.4f} is used.")
    if struct is None and precursor_mz is None:
        raise InputError("Enter a compound or formula, or a precursor m/z.")

    tol_ppm = common.number(inputs, "tol_ppm", "Tolerance", default=5.0, positive=True, maximum=100)
    tol_da = common.number(inputs, "tol_da", "Fragment tolerance", default=isomers.DEFAULT_TOL_DA, positive=True, maximum=0.1)
    dbs = _databases(inputs)
    pool = (ADDUCTS_POS if positive else ADDUCTS_NEG)[1:]
    if adduct is not None and adduct not in pool:
        raise InputError(f"{adduct} is not a {'positive' if positive else 'negative'}-mode adduct.")

    formula, name = None, None
    if struct and struct["kind"] == "formula":
        formula = re.sub(r"[+-]\d*$", "", struct["formula"])       # 电荷由加合物表示，搜索用中性式
        name = (inputs.get("structure") or {}).get("value") or formula
    elif struct:
        formula = _neutral_formula(struct["smiles"], warnings)
        name = common.record_name(inputs, struct["smiles"])
        if formula is None:
            warnings.append("The structure could not be converted to a molecular formula; the same-formula table is skipped.")
    else:
        name = f"m/z {precursor_mz:.4f}"

    if precursor_mz is None and formula:
        from mzdesk.ms1 import msbuddy_engine as bd
        used = adduct or bd.default_adduct(positive)
        theoretical, _ = bd.ion_mz_and_formula(formula, used)
        if theoretical:
            precursor_mz = theoretical
            warnings.append(f"No precursor m/z: the theoretical m/z {theoretical:.4f} ({used}) is used.")

    # ---- same formula (PubChem, online)
    same_rows, same_notes = [], []
    if formula:
        try:
            same = em.formula_groups(formula, dbs)
        except em.ExactMassError as exc:
            warnings.append(f"Same-formula candidates: {exc}")
            same_notes.append(f"PubChem could not be reached: {exc}")
        except ValueError as exc:
            raise InputError(str(exc)) from None
        else:
            same_rows = _same_rows(same)
            warnings.extend(same["warnings"])
            same_notes.extend(same["warnings"])
            if not same["examined"]:
                same_notes.append(f"PubChem lists no structure with the formula {formula}.")
            elif not same_rows:
                same_notes.append(f"{same['examined']} structures have the formula {formula}, but none is listed in the selected databases.")
    else:
        same_notes.append("Enter a structure or a molecular formula to search same-formula candidates.")

    # ---- isobaric interferences (PubChem, online)
    iso_rows, iso_notes = [], []
    if precursor_mz is not None:
        adducts = pool if adduct is None else [adduct]
        if not em.mass_windows(precursor_mz, adducts, tol_ppm):
            warnings.append("The precursor m/z is too small for the chosen adducts; isobaric interferences are skipped.")
            iso_notes.append("The precursor m/z is too small for the chosen adducts.")
        else:
            try:
                iso = em.isobaric(precursor_mz, adducts, tol_ppm, dbs, exclude_formula=formula)
            except em.ExactMassError as exc:
                warnings.append(f"Isobaric interferences: {exc}")
                iso_notes.append(f"PubChem could not be reached: {exc}")
            else:
                iso_rows = _isobaric_rows(iso)
                warnings.extend(iso["warnings"])
                iso_notes.extend(iso["warnings"])
                if not iso_rows:
                    iso_notes.append(f"No compound with a different formula is within {tol_ppm:g} ppm.")
    else:
        iso_notes.append("Enter a precursor m/z to search isobaric interferences.")

    # ---- isomers (MS-FINDER, local)
    iso_res, iso_note = None, None
    if struct is None:
        iso_note = "No structure: the isomer check needs a structure and an MS2 spectrum."
    elif struct["kind"] == "formula":
        iso_note = "Formula input: the isomer check needs a structure."
    elif not peaks:
        iso_note = "No MS2 peaks: the isomer check needs an MS2 spectrum."
    else:
        iso_res = isomers.discriminate([{"id": "compound", "smiles": struct["smiles"], "name": name}], peaks, precursor_mz,
                                       adduct or ("[M+H]+" if positive else "[M-H]-"), tol_da=tol_da, tol_ppm=tol_ppm, local_isomers=True)
        if iso_res["error"]:
            raise EngineError(iso_res["error"], iso_res["log"])

    if iso_res:
        verdict = _verdict(iso_res, name)
    else:
        counts = []
        if formula:
            what = "other structure" if struct and struct["kind"] == "smiles" else "structure"
            counts.append(f"{len(same_rows)} {what}{'s' if len(same_rows) != 1 else ''} share {formula}")
        if precursor_mz is not None:
            counts.append(f"{len(iso_rows)} isobaric formula{'s' if len(iso_rows) != 1 else ''} within {tol_ppm:g} ppm")
        verdict = ("; ".join(counts) + ". " if counts else "") + f"<b>No isomer verdict</b>: {html.escape(iso_note or 'the isomer check did not run.')}"

    return {"outputs": {
        "same_formula": table_data(same_rows, count=len(same_rows), help=same_notes),
        "isobaric": table_data(iso_rows, count=len(iso_rows), help=iso_notes),
        "isomers": table_data(_rows(iso_res) if iso_res else [], count=len(iso_res["isomers"]) if iso_res else 0,
                              help=(iso_res["notes"] if iso_res else [iso_note])),
    }, "warnings": warnings, "verdict": verdict}


RUNNERS = {"exclusivity": run_exclusivity}
