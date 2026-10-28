#!/usr/bin/env python3
"""Build ground-truth table for CASMI 2022 scoring -> results/truth.tsv

Source of truth: metadata/source/casmi2022_key_final.xlsx (official key, which already
contains the corrections). The challenge CSV and *_corrections.xlsx are cross-checked.
Checks: CSV vs key; RDKit InChIKey/formula from SMILES; expected precursor m/z vs given.
Run: python3 03_build_truth.py   (system python3 has rdkit, pandas, openpyxl)
"""
import os, re, sys
from collections import Counter
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors, Descriptors
from rdkit.Chem.inchi import MolToInchiKey
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.abspath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(BASE, "metadata", "source")
OUT = os.path.join(HERE, "..", "results", "truth.tsv")

E_MASS = 0.00054858
H, NA, CL = 1.00782503, 22.98976928, 34.96885268
NH4 = 14.0030740 + 4 * H
FA = 46.00547930           # HCOOH
ADD = {  # neutral M -> precursor m/z
    "[M+H]+": H - E_MASS, "[M-H]-": -(H - E_MASS)*1 - 0 + 0, "[M+Na]+": NA - E_MASS,
    "[M+NH4]+": NH4 - E_MASS, "[M+Cl]-": CL + E_MASS, "[M+FA-H]-": FA - (H - E_MASS),
}
ADD["[M-H]-"] = -(H - E_MASS)
# SIRIUS notation. SIRIUS's own profiles accept [M+NH4]+ and [M+FA-H]-; its canonical
# formula-style names are [M+H3N+H]+ / [M+CH2O2-H]- (see adduct_sirius_formula).
SIR = {"[M+NH4]+": "[M+NH4]+", "[M+FA-H]-": "[M+FA-H]-"}
SIR_FORMULA = {"[M+NH4]+": "[M+H3N+H]+", "[M+FA-H]-": "[M+CH2O2-H]-"}

def load_csv_key():
    k = pd.read_excel(os.path.join(SRC, "casmi2022_key_final.xlsx"), sheet_name=0, header=1)
    k = k.rename(columns={"Compound Number": "challenge_id", "Priority/Bonus": "priority",
                          "Precursor m/z (Da)": "precursor_mz", "RT [min]": "rt_min", "File": "file",
                          "Adduct": "adduct", "InChIKey": "inchikey", "SMILES": "smiles",
                          "Formula": "formula", "Monoisotopic Mass": "monoisotopic_mass"})
    for c in ["file", "adduct", "inchikey", "smiles", "formula", "priority"]:
        k[c] = k[c].astype(str).str.strip()
    k["challenge_id"] = k["challenge_id"].astype(int)
    return k

def hill(formula):
    return formula  # RDKit CalcMolFormula already returns Hill order

def main():
    key = load_csv_key()
    csv = pd.read_csv(os.path.join(BASE, "metadata", "casmi2022_challenges.csv"))
    for c in ["file", "adduct", "inchikey", "smiles", "formula", "priority"]:
        csv[c] = csv[c].astype(str)
    log = []
    # ---- check 1: CSV vs key
    m = key.merge(csv, on="challenge_id", suffixes=("_k", "_c"), how="outer", indicator=True)
    log.append(f"[1] rows key={len(key)} csv={len(csv)} unmatched ids={(m['_merge']!='both').sum()}")
    for col in ["file", "inchikey", "smiles", "formula", "adduct", "priority"]:
        a, b = m[col + "_k"], m[col + "_c"]
        raw = (a != b)
        strip = (a.str.strip() != b.str.strip())
        log.append(f"[1] {col}: diff_raw={raw.sum()} diff_after_strip={strip.sum()}")
        for _, r in m[strip].iterrows():
            log.append(f"     id={r.challenge_id} key={r[col+'_k']!r} csv={r[col+'_c']!r}")
    for col in ["precursor_mz", "monoisotopic_mass", "rt_min"]:
        d = (m[col + "_k"] - m[col + "_c"]).abs()
        log.append(f"[1] {col}: max abs diff={d.max():.3g}, n>1e-6={(d>1e-6).sum()}")
    # corrections files vs key
    for f in ["priority_corrections", "bonus_corrections"]:
        c = pd.read_excel(os.path.join(SRC, f + ".xlsx"), header=1)
        c = c.rename(columns={"Compound Number": "challenge_id", "Precursor m/z (Da)": "precursor_mz",
                              "Priority/Bonus": "priority", "File": "file", "RT [min]": "rt_min"})
        for _, r in c.iterrows():
            kr = key[key.challenge_id == r.challenge_id].iloc[0]
            cr = csv[csv.challenge_id == r.challenge_id].iloc[0]
            log.append(f"[1] {f} id={r.challenge_id}: file {r.file}|key {kr.file}|csv {cr.file}; "
                       f"mz {r.precursor_mz}|key {kr.precursor_mz}|csv {cr.precursor_mz}; "
                       f"rt {r.rt_min}|key {kr.rt_min}|csv {cr.rt_min}")
    # ---- build
    rows = []
    for _, r in key.sort_values("challenge_id").iterrows():
        notes = []
        if r.adduct != str(r.adduct).strip(): pass
        smi = r.smiles
        mol = Chem.MolFromSmiles(smi)
        ok = mol is not None
        ik_calc = f_calc = None
        mass_calc = None
        if ok:
            ik_calc = MolToInchiKey(mol)
            f_calc = rdMolDescriptors.CalcMolFormula(mol)
            mass_calc = Descriptors.ExactMolWt(mol)
            if "." in smi: notes.append("multi-component")
            net = sum(a.GetFormalCharge() for a in mol.GetAtoms())
            if net != 0: notes.append(f"net_charge_{net}")
            elif any(a.GetFormalCharge() != 0 for a in mol.GetAtoms()): notes.append("zwitterion_or_nitro_noxide(net0)")
            if re.search(r"[@/\\]", smi): notes.append("smiles_has_stereo")
            if ik_calc != r.inchikey:
                if ik_calc[:14] != r.inchikey[:14]: notes.append("inchikey14_MISMATCH")
                elif ik_calc[15:25] == "UHFFFAOYSA": notes.append("stereo_missing_in_smiles")
                else: notes.append("inchikey_block2_mismatch")
            if f_calc != r.formula: notes.append(f"formula_mismatch(rdkit={f_calc})")
            if abs(mass_calc - r.monoisotopic_mass) > 0.001: notes.append(f"mass_mismatch(rdkit={mass_calc:.4f})")
        else:
            notes.append("rdkit_parse_failed")
        ad = r.adduct.replace(" ", "")
        M = r.monoisotopic_mass
        exp = M + ADD[ad]
        ppm = (r.precursor_mz - exp) / exp * 1e6
        if abs(ppm) > 10: notes.append(f"precursor_off_{ppm:.1f}ppm")
        rows.append(dict(
            challenge_id=r.challenge_id, priority=r.priority,
            polarity="negative" if ad.endswith("-") else "positive",
            precursor_mz=r.precursor_mz, adduct=ad, adduct_sirius=SIR.get(ad, ad),
            adduct_sirius_formula=SIR_FORMULA.get(ad, ad),
            inchikey=r.inchikey, inchikey14=r.inchikey.split("-")[0], smiles=smi,
            formula=r.formula, monoisotopic_mass=M, rdkit_ok=ok,
            expected_mz=round(exp, 5), ppm_err=round(ppm, 2), notes=";".join(notes)))
    t = pd.DataFrame(rows)
    t.to_csv(OUT, sep="\t", index=False)
    # polarity cross-check with csv
    pc = t.merge(csv[["challenge_id", "polarity"]], on="challenge_id", suffixes=("", "_csv"))
    log.append(f"[1] polarity adduct-vs-csv diffs: {(pc.polarity != pc.polarity_csv).sum()}")
    # ---- reports
    log.append(f"[2] rdkit parse fail={(~t.rdkit_ok).sum()}; notes flagged: " +
               str(Counter(n.split('(')[0] for x in t.notes for n in x.split(';') if n and not n.startswith('precursor'))))
    log.append(f"[2] smiles w/ stereo marks: {t.notes.str.contains('smiles_has_stereo').sum()}; official keys with stereo block: {(~t.inchikey.str.contains('UHFFFAOYSA')).sum()}")
    for _, r in t[t.notes.str.contains("MISMATCH|block2|multi|failed|charge|zwitter|formula_|mass_")].iterrows():
        log.append(f"     id={r.challenge_id} {r.notes}")
    log.append(f"[3] ppm err: max abs={t.ppm_err.abs().max():.2f}, median abs={t.ppm_err.abs().median():.2f}, >10ppm: {(t.ppm_err.abs()>10).sum()}")
    for _, r in t[t.ppm_err.abs() > 10].iterrows():
        log.append(f"     id={r.challenge_id} {r.adduct} given={r.precursor_mz} exp={r.expected_mz} ppm={r.ppm_err}")
    log.append("[4] adducts: " + str(dict(t.adduct.value_counts())))
    els = []
    for f in t.formula:
        els.append(set(re.findall(r"[A-Z][a-z]?", f)))
    cnt = Counter(e for s in els for e in s)
    log.append("[4] element presence (n formulas): " + str(dict(sorted(cnt.items(), key=lambda x: -x[1]))))
    nonCHNO = [(i, f) for i, f, s in zip(t.challenge_id, t.formula, els) if s - set("CHNO")]
    log.append(f"[4] formulas with non-CHNO: {len(nonCHNO)}; halogen(F/Cl/Br/I): {sum(1 for s in els if s & set(['F','Cl','Br','I']))}; "
               f"S:{cnt['S']} P:{cnt['P']} Si:{cnt['Si']}")
    other = {e: n for e, n in cnt.items() if e not in set("C H N O S P F Cl Br I Si".split())}
    log.append(f"[4] other elements: {other}")
    dup = t[t.duplicated("inchikey14", keep=False)].sort_values("inchikey14")
    log.append(f"[4] duplicate inchikey14: {dup.inchikey14.nunique()} keys over {len(dup)} challenges; "
               f"duplicate full inchikey: {t.duplicated('inchikey', keep=False).sum()}")
    for k14, g in dup.groupby("inchikey14"):
        log.append(f"     {k14}: ids={list(g.challenge_id)} adducts={list(g.adduct)} full_same={g.inchikey.nunique()==1}")
    print("\n".join(log))

if __name__ == "__main__":
    main()
