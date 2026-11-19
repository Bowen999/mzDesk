"""Isomer check: can the MS/MS spectrum tell the top candidate apart from the other candidates with the same formula?

Step 7 of the ESI MS/MS Annotation workflow, and the Isomers table of Alternative explanations (module `exclusivity`).
Owner's rule, read pairwise: the top-1 candidate is "distinguished" from isomer j when some observed peak is explained by the top-1 and
not by j. Distinguished from every isomer -> verdict "single" (one answer). Otherwise "indistinguishable": the group is the top-1 plus
every isomer it is not distinguished from ("these cannot be separated"). Isomers that explain observed peaks the top-1 cannot are
"challengers": information only, they never change the verdict. "insufficient": no usable MS2 peak, no isomer to compare with, or the
engine could not run; the check never guesses. Every verdict only holds for the candidates that were considered; the notes say so.

Peaks: a peak counts when it is >= 1 % of the base peak; the precursor, its 13C isotope peaks and (singly charged) anything above them are
left out. A distinguishing or challenger peak must be >= 5 % of the base peak, so a fragment in the noise cannot decide the answer.
Which peaks a structure explains comes from MS-FINDER's in-silico fragmenter (mzdesk/ms2/msfinder_engine.explain_candidates), never from
rules of our own. An isomer MS-FINDER did not score cannot be ruled out, so it stays in the group.

Why MS-FINDER and not MetFrag (test run 2026-10-28): real library spectra of kaempferol (5) and luteolin (2), and of leucine and isoleucine
(1 each), each scored with the true and the wrong isomer first. Both engines gave the wrong isomer a "single" answer twice
(MS-FINDER 2 of 18 runs at 0.01 and 0.02 Da; MetFrag 2 of 27 at depth 2 / 3 and 0.01 / 0.02 Da), each time with the true compound listed
as a challenger. MS-FINDER kept leucine and isoleucine inseparable in every run, as expected, and its answers hardly moved with the
tolerance; MetFrag separated them at its usual depth 2 (0.01 Da) and changed several answers with depth and tolerance (plain bond
cleavage, no hydrogen-rearrangement rules). MetFrag is cheaper (one 1.4 s launch for all candidates); MS-FINDER needs one 5-12 s, 2.3 GB
run under the heavy-job lock, which is acceptable for at most 10 candidates.
"""
import math
import re

from mzdesk import adducts

MIN_PEAK = 0.01              # a peak counts when it is >= 1 % of the base peak
MIN_DISTINGUISHING = 0.05    # a distinguishing / challenger peak must be >= 5 % of the base peak
ISOTOPE_STEP = 1.003355      # 13C - 12C, Da
ISOTOPE_PEAKS = 2            # the precursor's M+1 and M+2 are left out
ECHO_TOL = 1e-4              # MS-FINDER echoes the m/z of the query file (5 decimals)
MAX_CANDIDATES = 10
DEFAULT_TOL_DA = 0.01        # MS-FINDER's own fragment tolerance (shipped parameter file)
SINGLE, INDISTINGUISHABLE, INSUFFICIENT = "single", "indistinguishable", "insufficient"
DATABASE_SOURCE = "MS-FINDER database"


def _finite(x):
    return isinstance(x, (int, float)) and math.isfinite(x)


def _charge(adduct):
    """|charge| of an adduct name: '[M+2H]2+' -> 2, '[M+H]+' -> 1."""
    m = re.search(r"\](\d*)[+-]$", adduct or "")
    return int(m.group(1)) if m and m.group(1) else 1


def _prepare(smiles, charged):
    """SMILES -> {"smiles", "inchikey", "formula", "mass"} of the largest fragment, neutralised unless the adduct is [M]+ / [M]-; or an error string."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, rdMolDescriptors
    from rdkit.Chem.MolStandardize import rdMolStandardize

    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles or "")
    if mol is None:
        return "the SMILES cannot be read"
    mol = rdMolStandardize.LargestFragmentChooser().choose(mol)
    if not charged:
        mol = rdMolStandardize.Uncharger().uncharge(mol)
        if Chem.GetFormalCharge(mol) != 0:
            return "the structure keeps a net charge, so it cannot form this adduct"
    inchikey = Chem.MolToInchiKey(mol)
    if not inchikey:
        return "no InChIKey could be made for it"
    return {"smiles": Chem.MolToSmiles(mol), "inchikey": inchikey, "formula": rdMolDescriptors.CalcMolFormula(mol).rstrip("+-"),
            "mass": Descriptors.ExactMolWt(mol)}


def _label(c):
    return c.get("name") or c.get("smiles") or str(c.get("id"))


def _first_name(name):
    """MS-FINDER joins synonyms with ';': keep the first."""
    return (name or "").split(";")[0].strip()


def discriminate(candidates, peaks, precursor_mz, adduct, *, tol_da=DEFAULT_TOL_DA, tol_ppm=10.0, engine="msfinder",
                 max_candidates=MAX_CANDIDATES, local_isomers=False):
    """Isomer verdict for the top-1 candidate (owner's pairwise rule, see the module docstring).

    candidates: [{"id", "smiles", "name"}], rank-1 first; only those with the top-1's formula are compared (at most max_candidates in all,
    by rank). peaks: [(m/z, intensity)]. precursor_mz: measured, or None (the theoretical m/z is then used). adduct: any spelling
    mzdesk/adducts.py knows. tol_da: fragment tolerance (Da); tol_ppm: precursor tolerance (ppm). engine: "msfinder", the only engine.
    local_isomers: also compare with the best-scoring same-formula structures of MS-FINDER's local databases (stereoisomers of a candidate
    are skipped: MS/MS cannot separate them), up to max_candidates in all; the Alternative explanations page uses this.
    Never raises for bad data or a missing engine: the verdict is then "insufficient" and notes[0] says why. Returns (JSON-able)
      {"verdict": single | indistinguishable | insufficient, "top_id",
       "group": ids that cannot be separated, top-1 first ([top_id] when single, [] when insufficient),
       "isomers": [{"id", "name", "smiles", "evaluated", "distinguished": True | False | None (not evaluated),
                    "distinguishing_peaks": [peak]: explained by the top-1, not by this isomer,
                    "isomer_only_peaks": [peak]: explained by this isomer, not by the top-1}],
       "challengers": [{"id", "name", "peaks": [peak]}],
       "candidates": [{"id", "name", "smiles", "inchikey", "formula", "source": input | MS-FINDER database, "evaluated", "score",
                       "explained": [peak], "explained_fraction": share of the considered intensity, or None}],
       "peaks_considered", "engine", "settings", "notes": [str], "error": engine failure message or None, "log", "seconds"}
    peak = {"mz", "intensity", "rel": % of the base peak, "formula": fragment ion formula, "smiles": fragment structure}.
    """
    if engine != "msfinder":
        raise ValueError(f"Unknown engine {engine!r}: the isomer check runs on MS-FINDER (see mzdesk/exclusivity/isomers.py for why).")
    from mzdesk.ms2 import msfinder_engine as mf

    res = {"verdict": INSUFFICIENT, "top_id": None, "group": [], "isomers": [], "challengers": [], "candidates": [], "peaks_considered": 0,
           "engine": f"MS-FINDER {mf.VERSION}", "notes": [], "error": None, "log": "", "seconds": 0.0,
           "settings": {"tol_da": tol_da, "tol_ppm": tol_ppm, "min_peak_pct": MIN_PEAK * 100, "min_distinguishing_pct": MIN_DISTINGUISHING * 100}}

    def stop(why):
        res["notes"].insert(0, why)
        return res

    clean = [(float(mz), float(i)) for mz, i in peaks or () if _finite(mz) and _finite(i) and mz > 0 and i > 0]
    if not clean:
        return stop("No MS2 peaks.")
    candidates = list(candidates or [])
    if not candidates:
        return stop("No candidate structure.")
    res["top_id"] = candidates[0].get("id")
    canonical, msf_adduct, sign = adducts.canonical(adduct), adducts.engine_name(adduct, "msfinder"), adducts.polarity(adduct)
    if msf_adduct is None:
        return stop(f"MS-FINDER cannot use the adduct {adduct or '(none given)'}.")
    charged = canonical in ("[M]+", "[M]-")

    # ---- candidates: the top-1, then the isomers by rank
    cands, other_formula, over_cap = [], 0, 0
    for rank, c in enumerate(candidates):
        prep = _prepare(c.get("smiles"), charged)
        if isinstance(prep, str):
            if rank == 0:
                return stop(f"The top candidate cannot be used: {prep}.")
            res["notes"].append(f"{_label(c)} is not compared: {prep}.")
        elif cands and prep["formula"] != cands[0]["formula"]:
            other_formula += 1                                  # not an isomer of the top-1
        elif any(prep["inchikey"] == d["inchikey"] for d in cands):
            res["notes"].append(f"{_label(c)} is not compared: the same structure as a candidate ranked higher.")
        elif len(cands) >= max_candidates:
            over_cap += 1
        else:
            cands.append({"id": c.get("id"), "name": c.get("name") or "", "source": "input", **prep})
    top = cands[0]
    if other_formula:
        res["notes"].append(f"{other_formula} candidate(s) with another formula than {top['formula']} are not isomers and are not compared.")
    if over_cap:
        res["notes"].append(f"{over_cap} lower-ranked isomer(s) are not compared: at most {max_candidates} candidates are.")

    from mzdesk.ms1 import msbuddy_engine
    theo, _ = msbuddy_engine.ion_mz_and_formula(top["formula"], canonical)
    if theo is None:
        return stop(f"{top['formula']} cannot form {canonical}.")
    if precursor_mz and abs(precursor_mz - theo) / theo * 1e6 > tol_ppm:
        return stop(f"The top candidate ({top['formula']}, {canonical} m/z {theo:.4f}) does not match the precursor m/z {precursor_mz:.4f} "
                    f"within {tol_ppm:g} ppm.")
    prec = precursor_mz or theo

    # ---- peaks: what MS-FINDER sees (>= 1 % of the base peak) and what can decide (without the precursor and its isotopes)
    base = max(i for _, i in clean)
    kept = sorted(sorted((p for p in clean if p[1] >= MIN_PEAK * base), key=lambda p: -p[1])[:mf.MAX_PEAKS])
    z = _charge(canonical)
    window = max(tol_da, prec * tol_ppm * 1e-6)
    considered = [(mz, i) for mz, i in kept
                  if not any(abs(mz - (prec + k * ISOTOPE_STEP / z)) <= window for k in range(ISOTOPE_PEAKS + 1)) and (z > 1 or mz < prec)]
    res["peaks_considered"] = len(considered)
    if not considered:
        return stop("No fragment peak of at least 1 % of the base peak (the precursor and its isotope peaks are left out).")
    if len(cands) < 2 and not local_isomers:
        return stop(f"No isomer to compare with: no other candidate has the formula {top['formula']}.")

    # ---- one MS-FINDER run for all candidates
    try:
        run = mf.explain_candidates(theo, sign == "Positive", msf_adduct, kept, cands, tol_da=tol_da, local_db=local_isomers)
    except mf.MsFinderError as exc:
        res["error"], res["log"] = str(exc), exc.log
        return stop(f"MS-FINDER could not run: {exc}")
    except ValueError as exc:
        return stop(f"MS-FINDER cannot use this input: {exc}")
    res["seconds"] = run["seconds"]
    scored = {s["id"]: s for s in run["structures"] if s["id"] is not None}
    if top["id"] not in scored:
        reached = top["formula"] in {f["formula"] for f in run["formulas"][:mf.CANDIDATE_TOP_FORMULAS]}
        return stop("MS-FINDER did not score the top candidate: " + ("it could not fragment the structure." if reached else
                    f"{top['formula']} is not among the {mf.CANDIDATE_TOP_FORMULAS} formulas it ranked best for m/z {theo:.4f}."))

    if local_isomers:
        given = {c["inchikey"][:14] for c in cands}
        seen, stereo, n_local = set(given), 0, 0
        for s in sorted((s for s in run["structures"] if s["id"] is None and s["formula"] == top["formula"]), key=lambda s: -(s["score"] or 0)):
            block = s["inchikey"][:14]
            if block in given and s["inchikey"] not in {c["inchikey"] for c in cands}:
                stereo += 1
            if not block or block in seen or len(cands) >= max_candidates:
                continue
            seen.add(block)
            cid = f"msfinder:{s['inchikey']}"
            cands.append({"id": cid, "name": _first_name(s["name"]), "smiles": s["smiles"], "inchikey": s["inchikey"], "formula": s["formula"],
                          "mass": None, "source": DATABASE_SOURCE})
            scored[cid] = s
            n_local += 1
        res["notes"].append(f"Isomers from MS-FINDER's local structure databases: the {n_local} best-scoring other structures with the formula "
                            f"{top['formula']}" + (f"; {stereo} stereoisomer(s) of a candidate are not compared (MS/MS does not tell them apart)."
                                                   if stereo else "."))
        if len(cands) < 2:
            return stop(f"No isomer to compare with: MS-FINDER's local databases have no other structure with the formula {top['formula']}.")

    # ---- which considered peak each structure explains
    def explained(s):
        out = {}
        for f in s["fragments"]:
            k = min(range(len(considered)), key=lambda n: abs(considered[n][0] - f["mz"])) if _finite(f["mz"]) else None
            if k is not None and abs(considered[k][0] - f["mz"]) <= ECHO_TOL:
                out.setdefault(k, f)
        return out

    def peak(k, frag):
        mz, i = considered[k]
        return {"mz": mz, "intensity": i, "rel": round(100 * i / base, 2), "formula": frag["formula"], "smiles": frag["smiles"]}

    expl = {c["id"]: explained(scored[c["id"]]) for c in cands if c["id"] in scored}
    if len(expl) < 2:
        return stop("MS-FINDER did not score any of the isomers, so none could be compared.")
    total = sum(i for _, i in considered)
    strong = {k for k, (_, i) in enumerate(considered) if i >= MIN_DISTINGUISHING * base}
    for c in cands:
        e = expl.get(c["id"])
        res["candidates"].append({"id": c["id"], "name": c["name"], "smiles": c["smiles"], "inchikey": c["inchikey"], "formula": c["formula"],
                                  "source": c["source"], "evaluated": e is not None, "score": scored[c["id"]]["score"] if e is not None else None,
                                  "explained": [peak(k, f) for k, f in sorted((e or {}).items())],
                                  "explained_fraction": round(sum(considered[k][1] for k in e) / total, 4) if e is not None else None})

    # ---- the pairwise rule
    mine = expl[top["id"]]
    group, unscored = [top["id"]], []
    for c in cands[1:]:
        other = expl.get(c["id"])
        if other is None:
            unscored.append(_label(c))
            only_top, only_other = [], []
        else:
            only_top = [peak(k, mine[k]) for k in sorted(strong & mine.keys() - other.keys())]
            only_other = [peak(k, other[k]) for k in sorted(strong & other.keys() - mine.keys())]
        distinguished = bool(only_top) if other is not None else None
        res["isomers"].append({"id": c["id"], "name": c["name"], "smiles": c["smiles"], "evaluated": other is not None,
                               "distinguished": distinguished, "distinguishing_peaks": only_top, "isomer_only_peaks": only_other})
        if not distinguished:
            group.append(c["id"])
        if only_other:
            res["challengers"].append({"id": c["id"], "name": c["name"], "peaks": only_other})
    res["verdict"], res["group"] = (SINGLE, [top["id"]]) if len(group) == 1 else (INDISTINGUISHABLE, group)

    n = len(cands) - 1
    res["notes"][:0] = [f"Relative to the {n + 1} candidates considered (the top candidate and {n} isomer{'s' if n != 1 else ''}).",
                        f"{len(considered)} of {len(clean)} peaks considered (≥ 1 % of the base peak, no precursor or isotope peaks); "
                        "distinguishing peaks need ≥ 5 %."]
    if unscored:
        res["notes"].append(f"MS-FINDER did not score {', '.join(unscored)}; kept in the group.")
    res["notes"].append(f"MS-FINDER {mf.VERSION}: fragments within {tol_da:g} Da, tree depth 2; formula search from the theoretical m/z {theo:.4f}.")
    return res
