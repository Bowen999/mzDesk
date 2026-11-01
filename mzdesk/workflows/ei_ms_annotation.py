"""EI MS Annotation workflow: one EI spectrum in, candidate metabolites (named as the original, underivatized compounds) out, with evidence.

Steps (spec.py stages), each reusing an existing engine; the result contract and the step runner are in mzdesk/workflows/common.py:
  1 Spectrum cleanup        ms2_clean.clean_spectrum: merged to unit mass if it has decimals, peaks under NOISE_PCT of the base peak dropped
                            (the molecular ion never), base peak = 100
  2 EI library search       library_search.search(..., "EI", ...) on the installed libraries, the default metrics; hits from LIST_THRESHOLD (the ESI
                            workflow's list threshold) on form the library tier
  3 Derivatization markers  ms2/ei_markers.py: m/z 73 (backed by 75 or 147), 147, 115, M-15, M-31, M-57 (M-90 and M-180 are only shown) -> TMS evidence,
                            suggested TMS count, inferred derivatization (Auto)
  4 Candidates              library hits (best 10 compounds, derivative entries mapped back to the original) + the user's list (SMILES, or a name looked up
                            online), merged by the first block of the original's InChIKey (ei_candidates.py)
  5 Derivatives             platform_fit/derivatization.derivatize: the derivative structure of each candidate for the chosen reagent
                            (Auto with TMS but MeOx undecided: TMS and MeOx + TMS for compounds with a C=O); a library entry that is itself a derivative
                            of the same kind keeps its own structure
  6 NEIMS prediction        ms2/neims_engine.predict on every derivative (skipped when the weights are not installed)
  7 Similarity and rerank   the query against each NEIMS spectrum with the library search's metric family (unit mass); ei_rerank.rerank:
                            library tier above NEIMS-only, then M consistency, then score; ei_rerank.ei_level: Level 2a / 3 / none
  8 Kovats RI               calculator.calculations.kovats_ri (computed and reported only; no RI prediction and no RI matching)
  9 Answer                  the top candidate named as the original compound, its level, and the evidence

Result: the generic one (common.make_result). answer.statement names the top candidate; answer.facts: Derivatization, TMS groups, Molecular ion,
Kovats RI, Top library match. blocks: the ranked "candidates" table and a mirror plot of the query against the top candidate's reference.
No statement here makes an accuracy claim: the workflow reports what matched, not how often such a match is right.

Form field ids: spectrum (the spectrum field: {"source", "value"}), derivatization, molecular_ion, candidates (textarea), alkanes (grid text,
"C8 8.00" per line), sample_rt.
"""
from __future__ import annotations

import logging
import math
import os
from dataclasses import dataclass, field
from urllib.parse import quote

from mzdesk.api import blocks as blk
from mzdesk.api import common
from mzdesk.api.common import InputError
from mzdesk.spec import EI_DERIVATIZATIONS, PIPELINES
from mzdesk.workflows.common import (STATUS_DONE, STATUS_SKIPPED, StepRunner, block, make_answer, make_result, notes_block, table_block)

log = logging.getLogger("mzdesk.workflows")

PIPELINE_ID = "ei_ms_annotation"
STEP_TITLES = {s["id"]: s["title"] for s in next(p for p in PIPELINES if p["id"] == PIPELINE_ID)["stages"]}
AUTO_DERIVATIZATION = EI_DERIVATIZATIONS[0]
MAX_INPUT_PEAKS = 5000          # a longer spectrum is rejected (400): it is not a cleaned EI spectrum
MAX_CANDIDATES = 50             # lines of the candidate list
NOISE_PCT = 0.5                 # step 1: peaks under this % of the base peak are dropped (the molecular ion is kept whatever its height)
UNIT_MASS_TOL = 0.02            # step 1: a peak further than this from an integer m/z makes the spectrum "high resolution": every peak is merged to its nominal mass
SEARCH_TOP_N = 100              # step 2: library spectra taken from the search before they are reduced to one per compound
LIBRARY_CANDIDATES = 10         # step 4: compounds taken from the library hits
TABLE_ROWS = 10                 # rows of the headline candidate table
SMILES_DRAWING = (260, 182)     # 2D structure drawn on hover over a name / derivative cell (as in the ESI workflow)
RECORD_URL = {"MassBank": "https://massbank.eu/MassBank/RecordDisplay?id={}",
              "GNPS": "https://gnps.ucsd.edu/ProteoSAFe/gnpslibraryspectrum.jsp?SpectrumID={}"}
NO_NEIMS = "NEIMS n/a"


@dataclass
class Query:
    """The form, parsed. ei: [(m/z, intensity)] in the order given. derivatization None = Auto (read the marker ions), else one of
    spec.EI_DERIVATIZATIONS[1:] ("None", "TMS", "MeOx + TMS", "TBDMS"). molecular_ion: m/z of the derivative's M+• or None. candidates: the
    non-blank lines of the list, unique, still SMILES or names (a name is looked up online by the step that uses it). alkanes: [(carbon number,
    RT)] sorted by carbon number, [] when there is no retention index, then sample_rt is None as well."""
    ei: list
    derivatization: str | None = None
    molecular_ion: float | None = None
    candidates: list = field(default_factory=list)
    alkanes: list = field(default_factory=list)
    sample_rt: float | None = None
    warnings: list = field(default_factory=list)


def parse_inputs(inputs):
    """Form JSON -> Query. Raises InputError (HTTP 400) with a message for the user."""
    from mzdesk.calculator import calculations

    warnings = []
    ei, _ = common.spectrum(inputs, "spectrum", "EI spectrum", warnings)
    if len(ei) > MAX_INPUT_PEAKS:
        raise InputError(f"The EI spectrum has {len(ei)} peaks; at most {MAX_INPUT_PEAKS} are accepted. Keep the strongest peaks.")
    choice = common.text(inputs, "derivatization") or AUTO_DERIVATIZATION
    if choice not in EI_DERIVATIZATIONS:
        raise InputError(f"Derivatization must be one of: {', '.join(EI_DERIVATIZATIONS)}.")
    molecular_ion = common.number(inputs, "molecular_ion", "Molecular ion m/z", positive=True)
    lines = (common.text(inputs, "candidates") or "").splitlines()
    candidates = list(dict.fromkeys(s for s in (line.strip() for line in lines) if s and not s.startswith("#")))
    if len(candidates) > MAX_CANDIDATES:
        raise InputError(f"At most {MAX_CANDIDATES} candidate structures are accepted; there are {len(candidates)}.")
    alkanes_text = common.text(inputs, "alkanes")
    sample_rt = common.number(inputs, "sample_rt", "Sample RT", positive=True)
    alkanes = []
    if alkanes_text and sample_rt is None:
        raise InputError("Sample RT is required when an alkane series is given.")
    if sample_rt is not None and not alkanes_text:
        raise InputError("The alkane series is required when a sample RT is given.")
    if alkanes_text:
        try:
            alkanes = sorted(calculations.parse_alkanes(alkanes_text))
            calculations.kovats_ri(alkanes, sample_rt)                    # only to reject a bad series now (repeated carbon number, RT not increasing)
        except ValueError as exc:
            raise InputError(f"Alkane series: {exc}") from None
    return Query(ei=ei, derivatization=None if choice == AUTO_DERIVATIZATION else choice, molecular_ion=molecular_ion, candidates=candidates,
                 alkanes=alkanes, sample_rt=sample_rt, warnings=warnings)


# ================================================================ run
def run_query(q, progress=None):
    """Run the nine steps. Every step is guarded: a step that fails or cannot run says so in its note and the later steps use what is there."""
    steps = StepRunner(STEP_TITLES, progress)
    ctx = {"q": q, "warnings": list(q.warnings), "peaks": [], "hits": [], "primary": None, "analysis": None, "mode": q.derivatization or "None",
           "basis": "selected" if q.derivatization else "assumed", "try_both": False, "candidates": [], "neims_ready": False, "ranked": [], "ri": None, "blocks": {}, "failed": []}
    for n, step in enumerate((_cleanup, _library, _markers, _candidates, _derivatives, _neims, _rerank, _kovats), 1):
        steps.guarded(n, lambda step=step, steps=steps: step(q, ctx, steps))
    ctx["failed"] = [r["title"] for r in steps.failed()]
    answer = steps.guarded(9, lambda: _answer_step(q, ctx))
    if answer is None:                              # the answer step itself failed: say so instead of showing nothing
        answer = make_answer("No answer could be built; see the failed step.")
    warnings = ctx["warnings"] + [r["note"] for r in steps.failed()]
    headline = [b for b in (ctx["blocks"].get("candidates"), ctx["blocks"].get("mirror")) if b]
    return make_result(PIPELINE_ID, steps, answer, headline, warnings)


# ---------------------------------------------------------------- helpers
def _score_text(value):
    return common.fmt(value, 3)


def _depicted(text, smiles):
    """Table cell with the structure drawn on hover."""
    svg = common.depiction_svg(smiles, SMILES_DRAWING) if smiles else None
    return blk.depicted(text, svg) if svg else text


def derivative_text(variant):
    """"3TMS", "1MeOx + 5TMS", "2TBDMS", "none" for an underivatized form."""
    if not variant:
        return common.MISSING
    n_silyl, n_meox, mode = variant["n_silyl"], variant["n_meox"], variant["mode"]
    if not n_silyl and not n_meox:
        return "none"
    silyl = "TBDMS" if mode == "TBDMS" else "TMS"
    return " + ".join(x for x in (f"{n_meox}MeOx" if n_meox else "", f"{n_silyl}{silyl}" if n_silyl else "") if x)


def _round_half_up(x):
    return int(math.floor(x + 0.5))


# ---------------------------------------------------------------- step 1
def clean_query(ei, molecular_ion=None):
    """Step 1 on a peak list [(m/z, intensity)] -> {"peaks" (m/z ascending, base peak = 100), "n_in", "n_merged" (peaks left after the merge to unit
    mass; equal to n_in for a unit-mass spectrum), "unit_mass" (False: the input had decimals and was merged), "n_noise" (peaks under NOISE_PCT),
    "m_restored" (the molecular ion was below the noise floor and was put back), "n_out"}."""
    from mzdesk.shared import ms2_clean

    valid = [(float(m), float(i)) for m, i in ei if math.isfinite(m) and math.isfinite(i) and m > 0 and i > 0]
    unit = all(abs(m - round(m)) <= UNIT_MASS_TOL for m, _ in valid)
    work = valid
    if not unit:
        merged = {}
        for m, i in valid:
            merged[_round_half_up(m)] = merged.get(_round_half_up(m), 0.0) + i
        work = [(float(k), v) for k, v in sorted(merged.items())]
    out = ms2_clean.clean_spectrum(work, noise_pct=NOISE_PCT, normalize="Base peak = 100", sort_by="m/z")
    peaks = list(out["peaks"])
    restored = False
    if molecular_ion is not None and work:
        base = max(i for _, i in work)
        near = [(m, i) for m, i in work if abs(m - molecular_ion) <= 0.5]
        if near and not any(abs(m - molecular_ion) <= 0.5 for m, _ in peaks):
            m, i = max(near, key=lambda p: p[1])
            peaks = sorted(peaks + [(m, i / base * 100.0)])
            restored = True
    return {"peaks": peaks, "n_in": len(ei), "n_merged": len(work), "unit_mass": unit, "n_noise": sum(1 for r in out["rows"] if r["reason"] == "noise"),
            "m_restored": restored, "n_out": len(peaks)}


def _cleanup(q, ctx, steps):
    from mzdesk.shared import ms2_plot

    res = clean_query(q.ei, q.molecular_ion)
    ctx["peaks"] = res["peaks"]
    if not res["peaks"]:
        raise ValueError("the spectrum has no usable peaks")
    notes = [f"{res['n_in']} peaks in, {res['n_out']} kept."]
    if not res["unit_mass"]:
        notes.append(f"The spectrum has decimal m/z values; the peaks were merged to unit mass ({res['n_merged']} nominal masses).")
    if res["n_noise"]:
        notes.append(f"{res['n_noise']} peaks under {NOISE_PCT:g} % of the base peak were dropped.")
    if res["m_restored"]:
        notes.append(f"The molecular ion m/z {q.molecular_ion:g} is below that floor and was kept.")
    notes.append("Intensities are scaled to a base peak of 100.")
    fig, cards = ms2_plot.spectrum_figure(res["peaks"], y_title="Relative intensity (%)")
    blocks = [block("figure", "query_spectrum", None, common.plotly(fig, peaks=cards), label="Query spectrum"), notes_block("cleanup_notes", "Cleanup", notes, collapsed=True)]
    note = f"{res['n_in']} peaks in, {res['n_out']} out"
    return STATUS_DONE, note, blocks, res


# ---------------------------------------------------------------- step 2
def _installed_libraries():
    """Names of the installed libraries that hold EI spectra (spec.LIBRARIES order) and the messages for the ones that are not installed."""
    from mzdesk import spec
    from mzdesk.ms2 import library_index as lx

    found, missing = [], []
    for name in spec.LIBRARIES:
        try:
            lib = lx.load(name)
        except lx.LibraryMissing as exc:
            missing.append(str(exc))
            continue
        if (lib.meta.get("counts") or {}).get("ei", 0) > 0:
            found.append(name)
    return found, missing


def _library(q, ctx, steps):
    from mzdesk import spec
    from mzdesk.ms2 import library_index as lx
    from mzdesk.ms2 import library_search as ls
    from mzdesk.ms2 import metric_names as mn
    from mzdesk.workflows import ei_rerank

    if not ctx["peaks"]:
        return STATUS_SKIPPED, "The cleaned spectrum has no peaks.", [], None
    libraries, missing = _installed_libraries()
    if not libraries:
        text = missing[0] if missing else "No installed library holds EI spectra."
        ctx["warnings"].append("No EI library spectra were searched: " + text)
        return STATUS_SKIPPED, "No EI library is installed, so no library spectrum was searched.", [], None
    metrics = list(spec.LIBRARY_DEFAULT_METRICS)
    list_threshold = ei_rerank.gate()[2]
    try:
        res = ls.search(ctx["peaks"], "EI", libraries, metrics, list_threshold, tolerance=ls.EI_TOLERANCE, top_n=SEARCH_TOP_N)
    except ImportError:                                           # ms_entropy missing: the cosine alone
        metrics = [mn.COSINE_GREEDY]
        res = ls.search(ctx["peaks"], "EI", libraries, metrics, list_threshold, tolerance=ls.EI_TOLERANCE, top_n=SEARCH_TOP_N)
    except lx.LibraryMissing as exc:
        ctx["warnings"].append(f"No EI library spectra were searched: {exc}")
        return STATUS_SKIPPED, "No EI library is installed, so no library spectrum was searched.", [], None
    primary = res["rank_metric"] or res["metrics"][0]
    pool = [h for h in res["hits"] if (h["scores"].get(primary) or 0.0) >= list_threshold]
    ctx.update(hits=pool, primary=primary, metrics=res["metrics"], library_names=libraries)
    cols = ["Rank", "Library", "Accession", "Name", "Matched peaks"] + res["metrics"]
    rows = []
    for h in pool[:TABLE_ROWS]:
        acc = h["accession"]
        rows.append([str(h["rank"]), h["library"], {"text": acc, "url": RECORD_URL[h["library"]].format(quote(acc, safe=""))} if acc else common.MISSING,
                     h["name"] or common.MISSING, str(h["matches"])] + [common.fmt(h["scores"].get(m), 3) for m in res["metrics"]])
    pre = res["prefilter"]
    info = [f"{', '.join(libraries)}: {pre['total']:,} EI spectra" + (f", the top {pre['kept']:,} by entropy scored with the selected metrics" if pre and pre["dropped"] else "") + "."
            if pre else "No EI spectra in the libraries.",
            f"Unit-mass match (±{ls.EI_TOLERANCE:g} Da). A hit joins the library tier from {primary} {list_threshold:.2f}."]
    info += [w for w in res["warnings"] if "EI: all" not in w] + [f"{name}: {text}" for name, text in res["notes"].items()]
    blocks = [table_block("library_hits", "EI library hits", cols, rows, count=len(pool)), notes_block("library_notes", "Search", info, collapsed=True)]
    note = f"{len(pool)} library spectra from {primary} {list_threshold:.2f}" if pool else f"no library spectrum reaches {primary} {list_threshold:.2f}"
    return STATUS_DONE, note, blocks, res


# ---------------------------------------------------------------- step 3
def _markers(q, ctx, steps):
    from mzdesk.ms2 import ei_markers as em

    a = em.analyze(ctx["peaks"], q.molecular_ion)
    ctx["analysis"] = a
    if q.derivatization:
        ctx.update(mode=q.derivatization, basis="selected")
        ctx["warnings"] += em.check_selection(a, q.derivatization)
    else:
        ctx.update(mode=a["derivatization"], basis="inferred")
        ctx["try_both"] = a["derivatization"] == em.AUTO_TMS and a["meox"] == em.MEOX_UNDECIDED
    ctx["warnings"] += a["warnings"]
    rows = []
    for m in a["markers"]:
        if m["rel"] is None:
            found, rel = blk.badge("needs M", "neutral"), common.MISSING
        else:
            found, rel = (blk.badge("yes", "good") if m["present"] else "no"), f"{m['rel']:.1f}"
        mz = common.MISSING if m["mz"] is None else f"{m['mz']:g}"
        rows.append([m["label"], mz, rel, f"≥ {m['min']:g}", found, m["meaning"]])
    detail = list(a["reasons"])
    if q.derivatization:
        detail.append(f"Derivatization {q.derivatization} was chosen by you and is used.")
    elif ctx["try_both"]:
        detail.append("Compounds with a C=O are tried both as TMS and as MeOx + TMS; the better match stays.")
    blocks = [table_block("markers", "Marker ions", ["Marker", "m/z", "% of base", "Needs (%)", "Found", "Means"], rows),
              notes_block("marker_notes", "Reading", detail, collapsed=True)]
    shown = f"{ctx['mode']} ({ctx['basis']})"
    note = shown + (f", suggests {a['tms_min']} or more TMS" if a["tms_min"] >= 2 and ctx["mode"] != "None" else "")
    return STATUS_DONE, note, blocks, a


# ---------------------------------------------------------------- step 4
def _candidates(q, ctx, steps):
    from mzdesk.identity import resolve
    from mzdesk.ms2 import ei_markers as em
    from mzdesk.ms2 import library_index as lx
    from mzdesk.workflows import ei_candidates as cands

    from_library = cands.from_hits(ctx["hits"], ctx["primary"], limit=LIBRARY_CANDIDATES) if ctx["hits"] else []
    for c in from_library:                           # the derivatization state the reference spectrum itself shows
        lib = c.library
        if lib["entry_mode"]:
            lib["state"] = em.family(lib["entry_mode"])
        elif lib["row"] is not None:
            try:
                mz, rel = lx.load(lib["library"]).peaks(lib["row"])
                lib["state"] = em.state_of_spectrum(list(zip(mz.tolist(), (rel * 100.0).tolist())))
            except Exception:                        # noqa: BLE001 - the state is a flag, not a requirement
                lib["state"] = None
    from_list, list_warnings = cands.from_lines(q.candidates, resolve.to_smiles) if q.candidates else ([], [])
    ctx["warnings"] += list_warnings
    labelled = [c.user_text for c in from_list if common.has_isotope_labels(c.parent_smiles)]
    if labelled:        # the derivative formula and mass of step 5 come from the unlabelled formula (the NEIMS prediction of step 6 uses the labelled structure)
        who = ("Candidate " if len(labelled) == 1 else "Candidates ") + ", ".join(f"“{line}”" for line in labelled) + (" has" if len(labelled) == 1 else " have")
        ctx["warnings"].append(common.isotope_note("The derivative formulas and masses are", who=who))
    candidates = cands.merge(from_library, from_list)
    lookup = cands.library_names([c.key for c in candidates if len(c.key) == 14 and (c.name_note or c.name == c.parent_smiles)])
    for c in candidates:                             # a name the derivative entry did not give, or a SMILES the user typed, gets the library's plain name
        found = lookup.get(c.key)
        if found and (c.name_note or c.name == c.parent_smiles):
            c.name, c.name_note = found, ""
    ctx["candidates"] = candidates
    if not candidates:
        return STATUS_DONE, "no candidate compound", [], candidates
    rows = []
    for c in candidates:
        origin = " + ".join({"library": "EI library", "list": "your list"}[s] for s in c.sources)
        hit = c.library
        rows.append([_depicted(c.name, c.parent_smiles), origin,
                     {"text": hit["accession"], "url": RECORD_URL[hit["library"]].format(quote(hit["accession"], safe=""))} if hit and hit["accession"] else common.MISSING,
                     _score_text(hit["score"]) if hit else common.MISSING])
    blocks = [table_block("candidate_pool", "Candidate compounds", ["Compound", "From", "Library spectrum", ctx["primary"] or "Score"], rows)]
    n_lib, n_list = sum("library" in c.sources for c in candidates), sum("list" in c.sources for c in candidates)
    return STATUS_DONE, f"{len(candidates)} compounds ({n_lib} from the library, {n_list} from your list)", blocks, candidates


# ---------------------------------------------------------------- step 5
def _variant(smiles, mode, own):
    """A derivative structure as a dict from a SMILES that already is the derivative (a library entry) or from the derivatization engine."""
    from rdkit import Chem

    from mzdesk.platform_fit import derivatization as d

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("RDKit could not parse the derivative")
    formula = d.parent_formula(mol)
    info = d.underivatize_info(smiles)
    return {"mode": mode, "smiles": smiles, "formula": d.formula_string(formula), "mass": d.formula_mass(formula), "nominal": d.nominal_mass(formula),
            "n_silyl": info["n_tms"] + info["n_tbdms"], "n_meox": info["n_meox"], "own": own, "reagent": d.EI_MODES.get(mode)}


def build_variants(c, modes, molecular_ion=None):
    """The derivative structures of Candidate ``c`` for each mode of ``modes`` ("None" | "TMS" | "MeOx + TMS" | "TBDMS"): a list of variant dicts
    {mode, smiles, formula, mass, nominal, n_silyl, n_meox, own, reagent, m_status, m_k}. A library entry whose stored structure is the derivative of
    that very kind (same reagent, or TMS for a compound without a C=O under MeOx + TMS) keeps its own structure."""
    from mzdesk.platform_fit import derivatization as d
    from mzdesk.workflows import ei_rerank

    out = []
    if not c.parent_smiles:
        return out
    entry = c.library["smiles"] if c.library and c.library["entry_mode"] else None
    entry_mode = c.library["entry_mode"] if entry else None
    for mode in modes:
        if mode == "None":
            v = _variant(c.parent_smiles, "None", False)
        else:
            reagent = d.EI_MODES[mode]
            if entry and (entry_mode == mode or (mode == "MeOx + TMS" and entry_mode == "TMS" and not d.carbonyl_count(c.parent_smiles))):
                v = _variant(entry, mode, True)
            else:
                dv = d.derivatize(c.parent_smiles, reagent)
                v = _variant(dv["smiles"], mode, False)
        v["m_status"], v["m_k"] = ei_rerank.m_status(v, molecular_ion)
        out.append(v)
    return out


def _derivatives(q, ctx, steps):
    from mzdesk.platform_fit import derivatization as d

    candidates, mode = ctx["candidates"], ctx["mode"]
    if not candidates:
        return STATUS_DONE, "no candidates to derivatize", [], None
    rows, failed = [], 0
    for c in candidates:
        if not c.parent_smiles:
            c.flags.append("no structure")
            continue
        modes = [mode]
        if ctx["try_both"] and d.carbonyl_count(c.parent_smiles):
            modes = ["TMS", "MeOx + TMS"]
        try:
            c.variants = build_variants(c, modes, q.molecular_ion)
        except Exception as exc:                     # noqa: BLE001 - one structure the engine cannot derivatize must not stop the others
            log.warning("derivative of %s failed: %s", c.name, exc)
            c.flags.append("derivative failed")
            failed += 1
            continue
        for v in c.variants:
            fit = {"ok": "fits M", "partial": f"fits M with {v['m_k']} of {v['n_silyl']} silyl groups", "mismatch": "does not fit M"}.get(v["m_status"], common.MISSING)
            rows.append([_depicted(c.name, c.parent_smiles), v["mode"] if v["n_silyl"] or v["n_meox"] else "none", _depicted(derivative_text(v), v["smiles"]),
                         v["formula"], str(v["nominal"]), fit])
    cols = ["Compound", "Derivatization", "Groups", "Formula", "Nominal mass", "Typed M"]
    blocks = [table_block("derivatives", "Derivatives", cols, rows)] if rows else []
    note = f"{len(rows)} derivative structures" + (f", {failed} failed" if failed else "")
    return STATUS_DONE, note, blocks, None


# ---------------------------------------------------------------- step 6
def _neims(q, ctx, steps):
    from mzdesk.ms2 import neims_engine as ne
    from mzdesk.ms2 import neims_setup

    variants = [v for c in ctx["candidates"] for v in c.variants]
    if not variants:
        return STATUS_SKIPPED, "No derivative structure to predict.", [], None
    if not os.path.isfile(ne.weights_path()):
        ctx["warnings"].append("NEIMS weights are not installed, so candidates are ranked by the library evidence only. Run: " + neims_setup.SETUP_COMMAND)
        return STATUS_SKIPPED, "NEIMS weights are not installed.", [], None
    cache, rows, errors = {}, [], 0
    for i, v in enumerate(variants, 1):
        smiles = v["smiles"]
        if smiles not in cache:
            try:
                cache[smiles] = ne.predict(smiles)
            except ne.NEIMSError as exc:
                ctx["warnings"].append(f"NEIMS could not run: {exc}")
                return STATUS_SKIPPED, "NEIMS could not run.", [], None
            except ValueError as exc:
                cache[smiles] = {"error": str(exc)}
            if i % 5 == 0:
                steps.emit(6, "running", f"predicted {i}/{len(variants)}")
        res = cache[smiles]
        if "error" in res:
            v["neims_error"], v["neims_peaks"] = res["error"], None
            errors += 1
        else:
            v["neims_error"], v["neims_peaks"] = "", res["peaks"]
    for c in ctx["candidates"]:
        for v in c.variants:
            peaks = v.get("neims_peaks")
            base = max(peaks, key=lambda p: p[1])[0] if peaks else None
            rows.append([c.name, derivative_text(v), str(len(peaks)) if peaks else common.MISSING, str(base) if base is not None else common.MISSING,
                         v["neims_error"] or "predicted"])
    ctx["neims_ready"] = True
    blocks = [table_block("neims", "NEIMS predictions", ["Compound", "Derivative", "Peaks", "Base peak m/z", "Status"], rows)]
    return STATUS_DONE, f"{len(variants) - errors} of {len(variants)} spectra predicted", blocks, None


# ---------------------------------------------------------------- step 7
def _scorer(ctx):
    """The query's scorer (unit mass, the library search's metrics) and the primary metric name."""
    from mzdesk import spec
    from mzdesk.ms2 import library_search as ls
    from mzdesk.ms2 import metric_names as mn
    from mzdesk.ms2 import ms2_similarity as sim

    metrics = ls.order_metrics(ctx.get("metrics") or spec.LIBRARY_DEFAULT_METRICS)
    q_mz, q_int = sim.clean_peaks(*sim.normalize_peaks_input(ctx["peaks"]))
    try:
        scorer = ls._Scorer(q_mz, q_int, None, ls.EI_TOLERANCE, metrics, True)
    except ImportError:
        metrics = [mn.COSINE_GREEDY]
        scorer = ls._Scorer(q_mz, q_int, None, ls.EI_TOLERANCE, metrics, True)
    primary = ctx.get("primary") if ctx.get("primary") in metrics else metrics[0]
    return scorer, primary


def score_peaks(scorer, primary, peaks):
    """Similarity (primary metric) and matched peak count of a reference spectrum [(mz, intensity)] against the query."""
    from mzdesk.ms2 import ms2_similarity as sim

    mz, inten = sim.clean_peaks(*sim.normalize_peaks_input(peaks))
    got, matches = scorer.score(mz, inten, None)
    return got.get(primary), matches


def _flags(c, ctx):
    q = ctx["q"]
    out = []
    v = c.chosen
    if c.m_status == "mismatch" and v:
        out.append(f"M mismatch ({v['nominal']} vs {q.molecular_ion:g})")
    elif c.m_status == "partial":
        out.append(f"partial silylation ({c.m_partial} of {v['n_silyl']})")
    if c.library and not c.state_ok:
        out.append("reference " + {"None": "underivatized"}.get(c.library["state"], str(c.library["state"])))
    if c.parent_smiles is None:
        out.append("no structure")
    if c.variants and v and v.get("neims_error"):
        out.append(NO_NEIMS)
    elif c.variants and not ctx["neims_ready"]:
        out.append(NO_NEIMS)
    out += [f for f in c.flags if f not in out]
    return out


def _rerank(q, ctx, steps):
    from mzdesk.workflows import ei_rerank as rr

    candidates = ctx["candidates"]
    if not candidates:
        return STATUS_DONE, "no candidates to rank", [], []
    scorer, primary = _scorer(ctx)
    ctx["rank_metric"] = primary
    rank_of = {None: 0, "ok": 0, "partial": 1, "mismatch": 2}
    for c in candidates:
        for v in c.variants:
            v["neims_score"], v["neims_matches"] = (score_peaks(scorer, primary, v["neims_peaks"]) if v.get("neims_peaks") else (None, 0))
        if c.variants:                                          # the form that stays: one that fits M, then the better NEIMS match
            c.chosen = min(c.variants, key=lambda v: (rank_of[v["m_status"]], -(v["neims_score"] if v["neims_score"] is not None else -1.0)))
            c.neims_score, c.m_status, c.m_partial = c.chosen["neims_score"], c.chosen["m_status"], c.chosen["m_k"]
        elif q.molecular_ion is not None:
            c.m_status = rr.M_UNKNOWN                            # no derivative structure to compare with the typed M
        c.tier = rr.TIER_LIBRARY if c.library else rr.TIER_NEIMS
        c.lib_score = c.library["score"] if c.library else None
        c.state_ok = bool(c.library is None or rr.state_matches(ctx["mode"], c.library["state"]))
        c.flags = _flags(c, ctx)
        c.level, c.level_reason = rr.ei_level(c, len(ctx["peaks"]))
    ranked = rr.rerank(candidates)
    ctx["ranked"] = ranked
    table = _candidate_table(ctx, ranked)
    ctx["blocks"]["candidates"] = table
    note = f"{len(ranked)} candidates ranked: {sum(c.tier == rr.TIER_LIBRARY for c in ranked)} with a library spectrum, {sum(c.tier == rr.TIER_NEIMS for c in ranked)} predicted only"
    return STATUS_DONE, note, [table], ranked


def _candidate_table(ctx, ranked):
    from mzdesk.annotation.levels import level_label
    from mzdesk.workflows import ei_rerank as rr

    gate = rr.gate()[0]
    rows = []
    for c in ranked[:TABLE_ROWS]:
        v = c.chosen
        derivative = (_depicted(f"{derivative_text(v)} · M {v['nominal']}" if v["n_silyl"] or v["n_meox"] else "none", v["smiles"])) if v else common.MISSING
        evidence = "EI library" if c.tier == rr.TIER_LIBRARY else "NEIMS only"
        lib = _score_text(c.lib_score) if c.lib_score is not None else common.MISSING
        if c.lib_score is not None and c.lib_score >= gate:
            lib = blk.badge(lib, "good")
        level = blk.badge(c.level, "good" if c.level == "2a" else "neutral") if c.level else common.MISSING
        help_lines = [c.level_reason]
        if c.level:
            help_lines.append(level_label(c.level))
        help_lines.append("From: " + " + ".join({"library": f"EI library ({c.n_spectra} spectra)", "list": "your list"}[s] for s in c.sources) + ".")
        if c.library:
            help_lines.append(f"Best library spectrum: {c.library['library']} {c.library['accession']} {c.library['name']}.")
        if v and v.get("neims_error"):
            help_lines.append(f"NEIMS: {v['neims_error']}")
        rows.append(blk.row([str(c.rank), _depicted(c.name, c.parent_smiles), derivative, evidence, lib,
                             _score_text(c.neims_score) if c.neims_score is not None else common.MISSING, level, "; ".join(c.flags) or common.MISSING],
                            help=help_lines))
    cols = ["Rank", "Compound", "Derivative", "Evidence", "Library score", "NEIMS score", "Level", "Flags"]
    return table_block("candidates", "Candidates", cols, rows, count=len(ranked), row_help=True,
                       help=["Compounds are named as the original, underivatized metabolites.",
                             "A compound with an experimental library spectrum ranks above one with only a NEIMS-predicted spectrum; inside a group the score decides.",
                             f"Scores are {ctx.get('rank_metric') or 'similarity'} at unit mass."])


# ---------------------------------------------------------------- step 8
def kovats(q):
    """Kovats RI of the sample (van den Dool-Kratz, linear) -> the calculator's dict, or None without an alkane series."""
    from mzdesk.calculator import calculations

    return calculations.kovats_ri(q.alkanes, q.sample_rt) if q.alkanes and q.sample_rt is not None else None


def _kovats(q, ctx, steps):
    ri = kovats(q)
    if ri is None:
        return STATUS_SKIPPED, "No alkane series was given.", [], None
    ctx["ri"] = ri
    (n1, t1), (n2, t2) = ri["lower"], ri["upper"]
    rows = [["Kovats RI", f"{ri['ri']:.1f}"], ["Lower alkane", f"C{n1} at {t1:g}"], ["Upper alkane", f"C{n2} at {t2:g}"], ["Sample RT", f"{q.sample_rt:g}"]]
    if ri["extrapolated"]:
        rows.append(["Extrapolated", "yes: the sample lies outside the alkane series"])
        ctx["warnings"].append("The Kovats RI is extrapolated: the sample RT lies outside the alkane series.")
    blocks = [table_block("kovats", "Kovats RI", ["Item", "Value"], rows,
                          help="Linear (van den Dool-Kratz) index from the two alkanes around the sample. It is reported only: "
                               "the libraries have no retention index, so it does not change the ranking.")]
    return STATUS_DONE, f"RI {ri['ri']:.1f}" + (" (extrapolated)" if ri["extrapolated"] else ""), blocks, ri


# ---------------------------------------------------------------- step 9
def statement_for(c, n_query_peaks):
    """The answer sentence for the top candidate ``c``: "Glycine (as 3TMS derivative), Level 2a: MassBank EI match MSBNK-..., score 0.91."."""
    from mzdesk.annotation.levels import level_name
    from mzdesk.workflows import ei_rerank as rr

    v = c.chosen
    form = ("library entry without a structure" if c.parent_smiles is None else
            "derivative structure not built" if v is None else
            f"as {derivative_text(v)} derivative" if v["n_silyl"] or v["n_meox"] else "underivatized")
    head = f"{c.name} ({form})"
    if c.level is None:
        if c.lib_score is None and c.neims_score is None:
            return "No identification: no spectrum could be compared (no library hit, no NEIMS prediction for the candidates)."
        return f"No identification: the best candidate is {head}. {c.level_reason}"
    if c.tier == rr.TIER_LIBRARY:
        lib = c.library
        evidence = f"{lib['library']} EI match {lib['accession']}, score {lib['score']:.2f}"
        problems = rr.library_problems(c, n_query_peaks) if c.level == "3" else []
        if problems:
            evidence += f" ({'; '.join(problems)})"
    else:
        evidence = f"NEIMS-predicted spectrum match only, score {c.neims_score:.2f}"
    return f"{head}, {level_name(c.level)}: {evidence}."


def _answer_step(q, ctx):
    answer = _answer(q, ctx)
    notes = [c.level_reason and f"{c.rank}. {c.name}: {c.level_reason}" for c in ctx["ranked"][:TABLE_ROWS]]
    blocks = [notes_block("answer_notes", "Evidence", [n for n in notes if n], collapsed=True)] if ctx["ranked"] else []
    top = ctx["ranked"][0] if ctx["ranked"] else None
    from mzdesk.annotation.levels import level_name
    note = f"{top.name}, {level_name(top.level)}" if top and top.level else "no identification"
    return STATUS_DONE, note, blocks, answer


def _answer(q, ctx):
    """The result's answer (common.make_answer) from the ranked candidates and the markers; also builds the mirror plot block."""
    from mzdesk.annotation.levels import level_label, level_name

    ranked, a = ctx["ranked"], ctx["analysis"]
    facts = []
    if a is not None:
        top = ranked[0] if ranked else None
        mode = ctx["mode"]
        shown_mode = top.chosen["mode"] if top and top.chosen and ctx["try_both"] else mode
        extra = "; MeOx undecided" if ctx["try_both"] and not (top and top.chosen and top.chosen["mode"] != "TMS") else ""
        facts.append(["Derivatization", f"{shown_mode} ({ctx['basis']}{extra})"])
        if top and top.chosen and (top.chosen["n_silyl"] or top.chosen["n_meox"]):
            v = top.chosen
            label = "TBDMS groups" if v["mode"] == "TBDMS" else "TMS groups"
            count = str(v["n_silyl"]) + (f" + {v['n_meox']} MeOx" if v["n_meox"] else "")
            if a["tms_min"] >= 2 and v["mode"] != "TBDMS":
                count += f" (m/z 147 suggests {a['tms_min']} or more)"
            facts.append([label, count])
        if q.molecular_ion is not None:
            seen = {m["id"]: m for m in a["markers"]}
            hint = "M-15 seen" if seen["M-15"]["present"] else "M seen" if seen["M"]["present"] else "no M ion seen"
            facts.append(["Molecular ion", f"{q.molecular_ion:g} ({hint})"])
    if ctx["ri"] is not None:
        facts.append(["Kovats RI", f"{ctx['ri']['ri']:.1f}" + (" (extrapolated)" if ctx["ri"]["extrapolated"] else "")])
    best = max((c.library for c in ctx["candidates"] if c.library), key=lambda lib: lib["score"], default=None)
    facts.append(["Top library match", f"{best['name']} ({best['library']} {best['accession']}), score {best['score']:.2f}" if best else "none"])

    mirror = _mirror(ctx, ranked[0]) if ranked else None
    if mirror:
        ctx["blocks"]["mirror"] = mirror
    if not ranked:
        if ctx.get("failed"):
            why = "A step failed (" + ", ".join(ctx["failed"]) + "); see the warnings."
        else:
            why = ("No library spectrum reaches the list threshold and no candidate structure was given."
                   if ctx.get("primary") else "No EI library was searched and no candidate structure was given.")
        return make_answer("No candidate. " + why, None, None, facts)
    top = ranked[0]
    return make_answer(statement_for(top, len(ctx['peaks'])), level_name(top.level) if top.level else None, level_label(top.level) if top.level else None, facts)


def _mirror(ctx, top):
    """Mirror plot of the query against the top candidate's reference: its library spectrum, else its NEIMS prediction."""
    from mzdesk.ms2 import library_index as lx
    from mzdesk.ms2 import library_search as ls
    from mzdesk.ms2 import mirror_plot

    try:
        if top.library and top.library.get("row") is not None:
            mz, rel = lx.load(top.library["library"]).peaks(top.library["row"])
            reference, axis, label = list(zip(mz.tolist(), (rel * 100.0).tolist())), "Library", f"{top.library['library']} {top.library['accession']}"
        elif top.chosen and top.chosen.get("neims_peaks"):
            reference, axis, label = top.chosen["neims_peaks"], "NEIMS", "NEIMS prediction"
        else:
            return None
        fig, cards = mirror_plot.build_mirror_figure(ctx["peaks"], reference, tolerance=ls.EI_TOLERANCE, labels=("Query", axis), height=420)
        return block("figure", "mirror", None, common.plotly(fig, peaks=cards), label=f"Query vs {top.name} ({label})")
    except Exception:                                           # noqa: BLE001 - a plot must not stop the answer
        log.exception("mirror plot failed")
        return None
