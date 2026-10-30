"""Alternative explanations 页（exclusivity）：三张表。

1. Same-formula candidates：PubChem 里分子式完全相同的结构（在线），按 InChIKey 前 14 位合并连接方式，可看库归属。
2. Isobaric interferences：同一 m/z 窗口里分子式不同的化合物（PubChem 精确质量检索，在线），每个分子式一行。
3. Isomers：PubChem 同分子式候选（在线），用本地 MetFrag 碎片预测比较；脂质不做此检查。每行有碎片视图按钮，并给出结构中被彩色碎片覆盖的原子比例。

输入可以是结构、分子式，或只有 precursor m/z：只有 m/z 时跳过前两张表的分子式部分和异构体表。两张在线表失败只给 warning，
异构体表需要结构、MS2 和可用的 PubChem 候选。
"""
import html
import re

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from mzdesk.api import common
from mzdesk.api.blocks import badge, depicted, fragment_cell, hinted, table_data
from mzdesk.api.common import InputError, fmt, signed
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


def _coverage_cell(c):
    """'62 % · 13/21': non-hydrogen atoms inside the coloured fragments of the structure's fragment view, as a share of all of them."""
    cov = c.get("atom_coverage")
    return f"{cov['fraction'] * 100:.0f} % · {cov['atoms']}/{cov['heavy']}" if cov else common.MISSING


def _view_cell(c, adduct):
    """The row's "View" button (fragment view); a dash when there is no adduct MetFrag can use or no spectrum to show."""
    return fragment_cell(c["smiles"], adduct, c["name"]) if adduct else common.MISSING


def _rows(res, adduct):
    """adduct: MetFrag's spelling, the one the comparison used (None: no "View" buttons). The Fragments column is the second one: the page freezes it next to the first."""
    by_id = {i["id"]: i for i in res["isomers"]}
    rows = []
    for c in res["candidates"]:
        explained = (f"{len(c['explained'])} / {res['peaks_considered']} · {c['explained_fraction'] * 100:.0f} %"
                     if c["evaluated"] else "not scored")
        name = c["name"] or common.MISSING
        if c["id"] == res["top_id"]:
            rows.append([name, _view_cell(c, adduct), c["source"], _smiles_cell(c["smiles"]), explained, badge("your compound"), common.MISSING,
                         _coverage_cell(c)])
            continue
        iso = by_id.get(c["id"], {"distinguished": None, "distinguishing_peaks": [], "isomer_only_peaks": []})
        verdict = (badge("not scored", help="MetFrag did not score this structure, so it cannot be ruled out.") if iso["distinguished"] is None
                   else badge("yes", "good") if iso["distinguished"] else badge("no", "warn"))
        rows.append([name, _view_cell(c, adduct), c["source"], _smiles_cell(c["smiles"]), explained, verdict, _peaks_cell(iso["distinguishing_peaks"]),
                     _coverage_cell(c)])
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
        text = f"MetFrag predicts {me} can be distinguished from {n} tested isomer{'s' if n != 1 else ''}; this is not a unique identification."
    elif res["verdict"] == "indistinguishable":
        others = [names[i] for i in res["group"][1:]]
        shown = ", ".join(others[:NAMES_IN_VERDICT]) + (f" and {len(others) - NAMES_IN_VERDICT} more" if len(others) > NAMES_IN_VERDICT else "")
        text = f"<b>Not distinguished or not scored</b> by MetFrag: {me} and {shown} (of {n} isomer{'s' if n != 1 else ''} considered)."
    else:
        text = f"<b>No isomer verdict</b>: {html.escape(res['notes'][0].rstrip('.'))}."
    return text + (" " + html.escape(res["source_summary"]) if res.get("source_summary") else "")


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
    from mzdesk.exclusivity import isomers, pubchem_isomers
    from mzdesk.ms1 import exact_mass as em
    from mzdesk.ms2 import fragment_view

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
    same_rows, same_notes, same_ok = [], [], False
    if formula:
        try:
            same = em.formula_groups(formula, dbs)
        except em.ExactMassError as exc:
            warnings.append(f"Same-formula candidates: {exc}")
            same_notes.append(f"PubChem could not be reached: {exc}")
        except ValueError as exc:
            raise InputError(str(exc)) from None
        else:
            same_ok = True
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
    iso_rows, iso_notes, isobaric_ok = [], [], False
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
                isobaric_ok = True
                iso_rows = _isobaric_rows(iso)
                warnings.extend(iso["warnings"])
                iso_notes.extend(iso["warnings"])
                if not iso_rows:
                    iso_notes.append(f"No compound with a different formula is within {tol_ppm:g} ppm.")
    else:
        iso_notes.append("Enter a precursor m/z to search isobaric interferences.")

    # ---- isomers (PubChem online candidates, local MetFrag)
    iso_res, iso_note, used_adduct = None, None, None
    lipid = False
    if struct and struct["kind"] == "smiles":
        from mzdesk.platform_fit import lipid_metabolite_classifier as clf
        lipid = bool(clf.extract_features(struct["smiles"]).lipid_class)
    if struct is None:
        iso_note = "No structure: the isomer check needs a structure and an MS2 spectrum."
    elif struct["kind"] == "formula":
        iso_note = "Formula input: the isomer check needs a structure."
    elif lipid:
        iso_note = ("Isomer check skipped for this structure recognised as a lipid. Conventional MS/MS usually cannot reliably determine "
                    "C=C double-bond or sn-positions; specific methods or additional evidence may distinguish them.")
    elif not peaks:
        iso_note = "No MS2 peaks: the isomer check needs an MS2 spectrum."
    else:
        try:
            candidates = em.formula_structures(formula)
        except em.ExactMassError as exc:
            iso_note = f"PubChem isomer candidate search failed: {exc} Other isomers are not ruled out."
            warnings.append(iso_note)
        except ValueError as exc:
            iso_note = f"PubChem isomer candidate search could not run: {exc} Other isomers are not ruled out."
        else:
            used_adduct = adduct or ("[M+H]+" if positive else "[M-H]-")
            iso_res = pubchem_isomers.discriminate({"id": "compound", "smiles": struct["smiles"], "name": name}, candidates, peaks,
                                                   precursor_mz, used_adduct, tol_da=tol_da, tol_ppm=tol_ppm)
            warnings.extend(n for n in iso_res["notes"] if n.startswith(("Unavailable or unusable", "Not scored")))

    if iso_res:
        verdict = _verdict(iso_res, name)
    else:
        counts = []
        if formula and same_ok:
            counts.append(f"{len(same_rows)} connectivity group{'s' if len(same_rows) != 1 else ''} listed for {formula}")
        if precursor_mz is not None and isobaric_ok:
            counts.append(f"{len(iso_rows)} isobaric formula{'s' if len(iso_rows) != 1 else ''} within {tol_ppm:g} ppm")
        verdict = ("; ".join(counts) + ". " if counts else "") + f"<b>No isomer verdict</b>: {html.escape(iso_note or 'the isomer check did not run.')}"

    # the "View" buttons of the Isomers table: MetFrag's spelling of the adduct, and the peaks the comparison considered
    view_adduct = fragment_view.metfrag_adduct(used_adduct) if iso_res else None
    view_query = (fragment_view.query(peaks, precursor_mz, used_adduct, tol_da=tol_da, tol_ppm=tol_ppm) if view_adduct else None)
    return {"outputs": {
        "same_formula": table_data(same_rows, count=len(same_rows), help=same_notes),
        "isobaric": table_data(iso_rows, count=len(iso_rows), help=iso_notes),
        "isomers": table_data(_rows(iso_res, view_adduct if view_query else None) if iso_res else [], count=max(0, len(iso_res["candidates"]) - 1) if iso_res else 0,
                              help=(iso_res["notes"] if iso_res else [iso_note]), fragments=view_query),
    }, "warnings": warnings, "verdict": verdict}


RUNNERS = {"exclusivity": run_exclusivity}
