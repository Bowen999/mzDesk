"""Bounded PubChem candidate comparison for Alternative explanations.

The workflow's MS-FINDER isomer check remains separate. This comparison is a
prediction for the selected structures, not proof of a unique identification.
"""
import math
import time

from mzdesk import adducts
from mzdesk.exclusivity import isomers
from mzdesk.ms2 import metfrag_engine as mf

MAX_CANDIDATES = 10              # query plus at most nine different connectivities
MAX_PEAKS = 2000
TREE_DEPTH = 2
TIME_BUDGET_S = 120


def discriminate(top, search, peaks, precursor_mz, adduct, *, tol_da=0.01, tol_ppm=5):
    """Compare a user structure to the bounded PubChem sample using local MetFrag."""
    summary = (f"PubChem same-formula search: {search['retrieved']} CIDs retrieved (limit {search['cap']}; "
               f"{'truncated' if search['truncated'] else 'not truncated at this limit'}). "
               f"At most {MAX_CANDIDATES - 1} other connectivities can be compared. Not exhaustive; untested structures are not ruled out.")
    result = {"verdict": "insufficient", "top_id": top["id"], "group": [], "isomers": [], "candidates": [],
              "peaks_considered": 0, "notes": [summary], "source_summary": summary}

    def stop(reason):
        result["notes"].insert(0, reason)
        return result

    canonical = adducts.canonical(adduct)
    mf_adduct = adducts.engine_name(adduct, "metfrag")
    if not mf_adduct:
        return stop(f"MetFrag cannot use the adduct {adduct or '(none given)'}.")
    clean = [(float(mz), float(i)) for mz, i in peaks or ()
             if isinstance(mz, (int, float)) and isinstance(i, (int, float)) and math.isfinite(mz) and math.isfinite(i) and mz > 0 and i > 0]
    if not clean:
        return stop("No MS2 peaks.")
    prep = isomers._prepare(top.get("smiles"), canonical == "[M]+")
    if isinstance(prep, str):
        return stop(f"The compound cannot be used: {prep}.")
    if prep["formula"].rstrip("+-") != search["formula"]:
        return stop("The structure's formula changes on preparation; a same-formula comparison cannot be made safely.")
    from rdkit import Chem
    if len(Chem.GetMolFrags(Chem.MolFromSmiles(top["smiles"]))) != 1:
        return stop("The structure has disconnected parts; a same-formula comparison cannot be made safely.")

    cands = [{**top, **prep, "source": "Your compound"}]
    seen = {prep["inchikey"][:14]}
    invalid = duplicates = over_cap = 0
    for candidate in search["candidates"]:
        candidate_mol = Chem.MolFromSmiles(candidate.get("smiles") or "")
        if candidate_mol is None or len(Chem.GetMolFrags(candidate_mol)) != 1:
            invalid += 1
            continue
        p = isomers._prepare(candidate.get("smiles"), canonical == "[M]+")
        if isinstance(p, str) or p["formula"].rstrip("+-") != prep["formula"].rstrip("+-"):
            invalid += 1
            continue
        block = p["inchikey"][:14]
        if block in seen:
            duplicates += 1
            continue
        seen.add(block)
        if len(cands) >= MAX_CANDIDATES:
            over_cap += 1
            continue
        cands.append({**candidate, **p})

    result["source_summary"] = (
        f"PubChem same-formula search: {search['retrieved']} CIDs retrieved (limit {search['cap']}; "
        f"{'truncated' if search['truncated'] else 'not truncated at this limit'}); "
        f"{len(cands) - 1} other connectivities selected (limit {MAX_CANDIDATES - 1}, "
        f"{'truncated' if over_cap else 'not truncated'}). Not exhaustive; untested structures are not ruled out.")
    result["notes"][0] = result["source_summary"]
    if search["missing"] or search["wrong_formula"] or invalid:
        result["notes"].append(f"Unavailable or unusable PubChem structures: {search['missing']} missing properties, "
                               f"{search['wrong_formula']} different formula, {invalid} invalid or changed on preparation; not ruled out.")
    if duplicates:
        result["notes"].append(f"{duplicates} records with an already-seen connectivity were not separately compared; "
                               "stereoisomers are not assessed.")
    if len(cands) < 2:
        return stop("No other usable same-formula connectivity was found to compare; other isomers are not ruled out.")

    from mzdesk.ms1 import msbuddy_engine
    theo, _ = msbuddy_engine.ion_mz_and_formula(prep["formula"], canonical)
    if not theo:
        return stop(f"{prep['formula']} cannot form {canonical}.")
    if precursor_mz and abs(precursor_mz - theo) / theo * 1e6 > tol_ppm:
        return stop(f"The compound ({prep['formula']}, {canonical} m/z {theo:.4f}) does not match the precursor m/z "
                    f"{precursor_mz:.4f} within {tol_ppm:g} ppm.")
    prec = precursor_mz or theo
    base = max(i for _, i in clean)
    eligible = [p for p in clean if p[1] >= isomers.MIN_PEAK * base]
    kept = sorted(sorted(eligible, key=lambda p: -p[1])[:MAX_PEAKS])
    z = isomers._charge(canonical)
    window = max(tol_da, prec * tol_ppm * 1e-6)
    considered = [(mz, i) for mz, i in kept
                  if not any(abs(mz - (prec + k * isomers.ISOTOPE_STEP / z)) <= window for k in range(isomers.ISOTOPE_PEAKS + 1))
                  and (z > 1 or mz < prec)]
    result["peaks_considered"] = len(considered)
    if not considered:
        return stop("No fragment peak of at least 1 % of the base peak remains after excluding the precursor and its isotopes.")
    if len(eligible) > MAX_PEAKS:
        result["notes"].append(f"Only the {MAX_PEAKS} most intense peaks were used.")

    explanations = {}
    failures = []
    deadline = time.monotonic() + TIME_BUDGET_S
    for c in cands:
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            failures.append(f"{c['name'] or c['id']}: comparison time limit reached")
            continue
        try:
            run = mf.run_metfrag(c["smiles"], considered, mf_adduct, prec, TREE_DEPTH, tol_da, 0.0, timeout_s=remaining)
        except (mf.MetFragError, ValueError) as exc:
            failures.append(f"{c['name'] or c['id']}: {exc}")
            continue
        if not run.get("scored", False):
            failures.append(f"{c['name'] or c['id']}: MetFrag returned no result")
            continue
        explanations[c["id"]] = {i: row for i, row in enumerate(run["annotation"])
                                  if row.get("annotated") and row.get("kind") == "fragment"}
    if failures:
        result["notes"].append("Not scored, not ruled out: " + "; ".join(failures) + ".")
    if top["id"] not in explanations:
        return stop("MetFrag did not score your compound; no isomer verdict is possible.")

    def peak(i, frag):
        mz, intensity = considered[i]
        return {"mz": mz, "intensity": intensity, "rel": round(100 * intensity / base, 2),
                "formula": frag.get("ion_formula", ""), "smiles": frag.get("smiles", "")}

    total = sum(i for _, i in considered)
    for c in cands:
        e = explanations.get(c["id"])
        result["candidates"].append({"id": c["id"], "name": c.get("name", ""), "smiles": c["smiles"], "source": c["source"],
                                     "evaluated": e is not None, "explained": [peak(i, f) for i, f in sorted((e or {}).items())],
                                     "explained_fraction": round(sum(considered[i][1] for i in e) / total, 4) if e is not None else None})
    if len(explanations) < 2:
        return stop("MetFrag did not score any of the other candidates; no isomer verdict is possible.")
    mine = explanations[top["id"]]
    strong = {i for i, (_, intensity) in enumerate(considered) if intensity >= isomers.MIN_DISTINGUISHING * base}
    group = [top["id"]]
    for c in cands[1:]:
        e = explanations.get(c["id"])
        only_top = [peak(i, mine[i]) for i in sorted(strong & (mine.keys() - e.keys()))] if e is not None else []
        only_other = [peak(i, e[i]) for i in sorted(strong & (e.keys() - mine.keys()))] if e is not None else []
        distinguished = bool(only_top) if e is not None else None
        result["isomers"].append({"id": c["id"], "distinguished": distinguished, "distinguishing_peaks": only_top,
                                  "isomer_only_peaks": only_other})
        if not distinguished:
            group.append(c["id"])
    result["group"] = group
    result["verdict"] = "single" if len(group) == 1 else "indistinguishable"
    result["notes"].append(f"MetFrag depth {TREE_DEPTH}, fragment tolerance {tol_da:g} Da, relative tolerance 0 ppm. "
                           "Predicted fragments are model-dependent, not proof of identity. "
                           "Peaks ≥1% of the base peak are considered; distinguishing peaks must be ≥5%.")
    return result
