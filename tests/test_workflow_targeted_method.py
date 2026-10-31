"""Targeted method workflow: the pure rules (mzdesk/workflows/targeted_rules.py) and the nine steps end to end (mzdesk/workflows/targeted_method.py).

The end-to-end tests run the real engines with the Identity lookup replaced (no internet). MolGpKa weights are optional: the tests that need them
skip cleanly (python -m mzdesk.platform_fit.pka_setup); the NEIMS weights and RepoRT are optional too and never change what is asserted.
"""
import json
import re
import time
from pathlib import Path
from types import SimpleNamespace

import pytest
from fastapi.testclient import TestClient
from rdkit import Chem
from rdkit.Chem import inchi

from mzdesk import spec
from mzdesk.api import app as app_module
from mzdesk.api.runners import adducts as adducts_runner
from mzdesk.api.runners import ionization as ionization_runner
from mzdesk.identity import resolve
from mzdesk.methods import report
from mzdesk.platform_fit import pka_engine
from mzdesk.workflows import targeted_method as tm
from mzdesk.workflows import targeted_rules as tr

client = TestClient(app_module.app)
TARGETED = next(p for p in spec.PIPELINES if p["id"] == "targeted_method")
NO_MOLGPKA = not Path(pka_engine.weights_path()).is_file()
needs_molgpka = pytest.mark.skipif(NO_MOLGPKA, reason="MolGpKa weights not installed (python -m mzdesk.platform_fit.pka_setup)")
needs_report = pytest.mark.skipif(report.load() is None, reason=f"RepoRT not installed ({report.SETUP_COMMAND})")
CAFFEINE = spec.EXAMPLE_STRUCTURES["caffeine"]
GLUCOSE = spec.EXAMPLE_STRUCTURES["D-glucose"]
TOLUENE = "Cc1ccccc1"
SHEET_ITEMS = ["Can it be measured", "Platform", "Column", "Mobile phase / pH", "Polarity", "Main ion", "Theoretical m/z", "Derivatization",
               "Non-MS alternative", "Published methods"]
FORBIDDEN = re.compile(r"\bMRM\b|\bSRM\b|collision|\bCE\b|transition", re.I)         # the owner's decision: no transitions, no collision energies


@pytest.fixture(autouse=True)
def clean_jobs():
    app_module._jobs.clear()
    yield
    app_module._jobs.clear()


class FakeResolved:
    """What resolve.resolve returns for caffeine (the fields the Identity step reads)."""

    def __init__(self, **kw):
        self.name, self.inchikey, self.cid = "Caffeine", "RYYVLZVUVIJVGH-UHFFFAOYSA-N", 2519
        self.ids = {"HMDB": ["HMDB0001847"], "KEGG": ["C07481"], "ChEBI": ["CHEBI:27732"], "CAS": ["58-08-2"]}
        self.links = {"PubChem": "https://pubchem.ncbi.nlm.nih.gov/compound/2519", "HMDB": "https://hmdb.ca/metabolites/HMDB0001847"}
        self.warnings = []
        self.__dict__.update(kw)


def _offline_resolve(query, kind="Auto", **kw):
    raise resolve.ResolveError("Network lookup failed (PubChem/UniChem unreachable): test")


NAMES = {CAFFEINE: "Caffeine", GLUCOSE: "D-Glucose", TOLUENE: "Toluene"}
KEYS = {inchi.MolToInchiKey(Chem.MolFromSmiles(smiles)): name for smiles, name in NAMES.items()}


def _known_resolve(query, kind="Auto", **kw):
    """The Identity lookup answered from a table (no internet), by SMILES or InChIKey: a name for the compounds above, 'unreachable' for any other."""
    name = KEYS.get(query) if kind == "InChIKey" else NAMES.get(query)
    if name is None:
        return _offline_resolve(query, kind)
    return FakeResolved(name=name, ids={}, links={})


@pytest.fixture
def offline(monkeypatch):
    monkeypatch.setattr(resolve, "resolve", _offline_resolve)


@pytest.fixture(scope="module")
def offline_module():
    mp = pytest.MonkeyPatch()
    mp.setattr(resolve, "resolve", _known_resolve)
    yield
    mp.undo()


def form(smiles, **over):
    base = {"structure": {"mode": "SMILES / InChI", "value": smiles}, "ph": "", "additive": "Formic acid", "instrument": "Any"}
    return {**base, **over}


def run(smiles, progress=None, **over):
    return tm.run_query(tm.parse_inputs(form(smiles, **over)), progress)


_cache = {}


def cached(smiles, **over):
    """A real run, kept for the other tests of the module (the engines take a few seconds)."""
    key = (smiles, tuple(sorted(over.items())))
    if key not in _cache:
        _cache[key] = run(smiles, **over)
    return _cache[key]


def text(cell):
    return cell if isinstance(cell, str) else cell.get("badge") or cell.get("text") or ""


def sheet(result):
    """{item: (recommendation text, basis)} of the Method sheet."""
    block = result["blocks"][0]
    return {row[0]: (text(row[1]), row[2]) for row in block["data"]["rows"]}


def statuses(result):
    return {s["id"]: s["status"] for s in result["steps"]}


# ================================================================ rules: can it be measured
def test_helpers_read_the_engine_words():
    assert tr.esi_best("weak", "strong") == "strong" and tr.esi_best("moderate", "weak") == "moderate" and tr.esi_best(None, None) is None
    assert [tr.polarity_char(p) for p in ("+", "−", "-", "both", "poor", "x")] == ["+", "-", "-", "both", "poor", None]
    assert (tr.esi_label("+"), tr.esi_label("-")) == ("ESI+", "ESI−")
    assert [tr.omics_label(a, b) for a, b in ((True, True), (True, False), (False, True), (False, False))] == [
        "lipidomics and metabolomics", "lipidomics", "metabolomics", ""]
    assert tr.range_text([(0, 0), (4.1, 14)]) == "0 or 4.1–14"


@pytest.mark.parametrize("esi, gc, non_ms, verdict, route", [
    ("strong", "yes", None, "Yes", "LC-MS"), ("moderate", "after derivatization", None, "Yes", "LC-MS"), ("moderate", "no", None, "Yes", "LC-MS"),
    ("weak", "yes", None, "Yes", "GC-MS"),
    ("weak", "after derivatization", None, "Possible", "GC-MS"),
    ("weak", "no", ("CAD", "good"), "Possible", "non-MS"),
    ("weak", "no", ("RI", "possible"), "Unlikely", "none"), ("weak", "no", ("UV/PDA", "possible"), "Unlikely", "none"), ("weak", "no", None, "Unlikely", "none"),
])
def test_decide_any_takes_the_best_platform(esi, gc, non_ms, verdict, route):
    d = tr.decide("Any", esi, gc, non_ms)
    assert (d.verdict, d.route) == (verdict, route)
    assert d.ms_route == (tr.ROUTE_GC if route == "GC-MS" else tr.ROUTE_LC)


def test_decide_never_says_unlikely_for_a_missing_input():
    assert tr.decide("Any", None, None).verdict is None
    assert tr.decide("Any", "weak", None).verdict is None and tr.decide("Any", None, "no").verdict is None
    assert tr.decide("Any", None, "yes").verdict == "Yes" and tr.decide("Any", None, "yes").route == "GC-MS"
    assert tr.decide("LC-MS", None, "yes").verdict is None and tr.decide("GC-MS", "strong", None).verdict is None


def test_decide_names_the_lipid_or_metabolite_workflow_and_the_platform_label():
    assert tr.decide("Any", "strong", "yes", omics="lipidomics").platform == "LC-MS (lipidomics)"
    assert tr.decide("Any", "strong", "yes").platform == "LC-MS"
    assert tr.decide("Any", "weak", "after derivatization").platform == "GC-MS (derivatized)"
    assert tr.decide("Any", "weak", "no", ("CAD", "good")).platform == "HPLC-CAD (no MS)"
    assert tr.decide("Any", "weak", "no").platform == "None suitable by ESI or GC"


def test_decide_judges_only_the_chosen_instrument_and_says_so():
    lc_weak = tr.decide("LC-MS", "weak", "yes")
    assert (lc_weak.verdict, lc_weak.route, lc_weak.ms_route) == ("Unlikely", "LC-MS", "LC-MS")
    assert "APCI" in lc_weak.reason and "GC-MS would be volatile as is" in lc_weak.reason          # unsuitable, and what would suit
    assert tr.decide("LC-MS", "strong", "no").verdict == "Yes" and "judged on LC-MS only" in tr.decide("LC-MS", "strong", "no").reason
    # the same compound is "Yes" on Any (LC-MS) and only "Possible" when the instrument is GC-MS
    assert tr.decide("Any", "strong", "after derivatization").verdict == "Yes"
    gc = tr.decide("GC-MS", "strong", "after derivatization")
    assert (gc.verdict, gc.route, gc.ms_route) == ("Possible", "GC-MS", "GC-MS") and "judged on GC-MS only" in gc.reason
    assert tr.decide("GC-MS", "weak", "yes").verdict == "Yes"
    no = tr.decide("GC-MS", "strong", "no")
    assert no.verdict == "Unlikely" and "LC-MS would suit" in no.reason
    assert tr.decide("GC-MS", "weak", "no", ("CAD", "good")).verdict == "Unlikely"           # a filter ignores the non-MS fallback


# ================================================================ rules: separation, column, pH
ROWS = [["RP", "poor", "logD low"], ["HILIC", "good", "polar"], ["Ion-pair RP", "poor", "x"], ["Mixed-mode / ion exchange", "poor", "x"],
        ["Chiral", "not needed", "no stereocentre"], ["SEC", "poor", "x"]]


def test_best_mode_takes_the_highest_rating_and_breaks_ties_in_order():
    levels = tr.mode_levels(ROWS)
    assert levels == {"RP": "poor", "HILIC": "good", "Ion-pair RP": "poor", "Mixed-mode / ion exchange": "poor"}          # chiral and SEC are not LC modes here
    assert tr.best_mode(levels) == ("HILIC", "good")
    assert tr.best_mode({"RP": "good", "HILIC": "good", "Ion-pair RP": "possible", "Mixed-mode / ion exchange": "poor"}) == ("RP", "good")
    assert tr.best_mode({"RP": "poor", "HILIC": "possible", "Ion-pair RP": "possible", "Mixed-mode / ion exchange": "possible"}) == ("HILIC", "possible")
    assert tr.best_mode({"RP": "poor", "HILIC": "poor", "Ion-pair RP": "possible", "Mixed-mode / ion exchange": "poor"}) == ("Ion-pair RP", "possible")
    assert tr.best_mode({m: "poor" for m in tr.MODE_ORDER}) == ("RP", "poor")                          # nothing rated well: reversed phase is the default
    assert tr.best_mode(None) == (None, None)


@pytest.mark.parametrize("column, cls", [
    ("Waters ACQUITY UPLC HSS T3 (L1), 2.1 × 100 mm, 1.8 µm", "RP"), ("Waters ACQUITY UPLC BEH C18 (L1), 2.1 × 100 mm", "RP"),
    ("Restek Raptor Biphenyl (L11), 2.1 × 100 mm", "RP"), ("C18, 4.6 × 150 mm, 5 µm, 40 °C", "RP"), ("phenyl, 4.6 × 250 mm", "RP"),
    ("Merck SeQuant ZIC-HILIC (L114), 2.1 × 100 mm", "HILIC"), ("Waters XBridge BEH Amide (L68), 2.1 × 100 mm", "HILIC"),
    ("Thermo Scientific Accucore HILIC (L3), 2.1 × 100 mm", "HILIC"), ("cation-exchange, 7.8 × 300 mm, 9 µm, 85 °C", "IEX"),
    ("GC-PEG, 3.2 mm × 3 m", "GC"), ("GC-other, 3 mm × 1.5 m, 105 °C", "GC"), ("chiral, 4.6 × 250 mm", "chiral"),
    ("2.1 × 50 mm, 1.7 µm", None), ("silica, 4.6 × 250 mm", None), ("", None), (None, None)])
def test_column_class_reads_published_column_texts(column, cls):
    assert tr.column_class(column) == cls


def m(column, detector="MS", ref="RepoRT 0001", rt=None):
    return SimpleNamespace(column=column, detector=detector, reference=(ref, None), rt=rt)


C18, HILIC_COL, IEX_COL, GC_COL = "Waters BEH C18 (L1), 2.1 × 100 mm", "Merck ZIC-HILIC (L114), 2.1 × 100 mm", "cation-exchange, 7.8 × 300 mm", "GC-PEG, 3.2 mm × 3 m"


def test_lc_column_is_the_first_published_method_the_rules_do_not_reject():
    levels = {"RP": "good", "HILIC": "possible", "Ion-pair RP": "poor", "Mixed-mode / ion exchange": "poor"}
    pick = tr.choose_column("LC-MS", levels, [m(HILIC_COL), m(C18)], n_stereo=0)
    assert (pick.text, pick.source, pick.kind) == (C18, "published", "C18")                       # the better-rated class wins over the list order
    pick = tr.choose_column("LC-MS", levels, [m(C18, "MS"), m(C18.replace("BEH", "HSS"), "MS")])
    assert pick.text == C18                                                                       # same class: list order (best documented first)
    pick = tr.choose_column("LC-MS", levels, [m(C18.replace("BEH", "UV one"), "UV"), m(C18, "MS")])
    assert pick.text == C18 and pick.method.detector == "MS"                                      # an MS method is preferred over a UV one of the same class


def test_lc_column_never_uses_a_published_class_the_rules_rate_poor():
    levels = {"RP": "poor", "HILIC": "good", "Ion-pair RP": "poor", "Mixed-mode / ion exchange": "poor"}
    pick = tr.choose_column("LC-MS", levels, [m(C18), m(HILIC_COL)])
    assert (pick.text, pick.source, pick.kind) == (HILIC_COL, "published", "HILIC")
    pick = tr.choose_column("LC-MS", levels, [m(C18), m(IEX_COL), m("2.1 × 50 mm, 1.7 µm"), m(GC_COL)])
    assert pick.source == "rule" and pick.text == tr.TYPICAL_COLUMN["HILIC"] and pick.method is None   # nothing acceptable published: the typical column of the best mode


def test_published_retention_beats_the_predicted_rating_and_an_unretained_peak_is_not_used():
    # caffeine: the rules rate RP poor (logD -1.1), but a published C18 method retains it (RT 20 min) and the HILIC dataset does not (RT 1.1 min)
    levels = {"RP": "poor", "HILIC": "good", "Ion-pair RP": "poor", "Mixed-mode / ion exchange": "poor"}
    pick = tr.choose_column("LC-MS", levels, [m(HILIC_COL, rt="1.06 / 1.33"), m(C18, rt="20.37")])
    assert (pick.text, pick.source, pick.kind) == (C18, "published", "C18")
    # glucose: a C18 gradient elutes it at the dead time (1 min), the HILIC method retains it
    levels = {"RP": "poor", "HILIC": "good", "Ion-pair RP": "poor", "Mixed-mode / ion exchange": "poor"}
    assert tr.choose_column("LC-MS", levels, [m(C18, rt="0.98 / 0.99 / 1.01"), m(HILIC_COL, rt="6.86")]).text == HILIC_COL
    # an unretained peak is rejected even when the rules rate its class good; nothing retained -> the typical column of the best mode
    levels = {"RP": "good", "HILIC": "possible", "Ion-pair RP": "poor", "Mixed-mode / ion exchange": "poor"}
    pick = tr.choose_column("LC-MS", levels, [m(C18, rt="1.2")])
    assert pick.source == "rule" and pick.text == tr.TYPICAL_COLUMN["RP"]
    # both retained: the class the rules rate higher first; no retention time (JP18) needs a good or possible rating
    levels = {"RP": "possible", "HILIC": "good", "Ion-pair RP": "poor", "Mixed-mode / ion exchange": "poor"}
    assert tr.choose_column("LC-MS", levels, [m(C18, rt="5"), m(HILIC_COL, rt="3")]).text == HILIC_COL
    assert tr.choose_column("LC-MS", levels, [m("C18, 4.6 × 150 mm", "UV", rt=None)]).source == "published"
    assert tr.choose_column("LC-MS", {**levels, "RP": "poor"}, [m("C18, 4.6 × 150 mm", "UV", rt=None)]).source == "rule"
    assert tr.retention_min(m(C18, rt="1.06 / 1.33")) == 1.33 and tr.retention_min(m(C18, rt="20.37")) == 20.37 and tr.retention_min(m(C18)) is None
    assert tr.MIN_RETAINED_RT == 2.0


@pytest.mark.parametrize("mode, text, kind", [
    ("RP", "C18 (reversed phase)", "C18"), ("HILIC", "Amide or zwitterionic HILIC", "HILIC"),
    ("Ion-pair RP", "C18 with an ion-pair reagent (not MS-friendly)", "C18"), ("Mixed-mode / ion exchange", "Mixed-mode reversed-phase / ion-exchange", "mixed-mode")])
def test_lc_column_map_from_the_best_mode(mode, text, kind):
    levels = {k: ("good" if k == mode else "poor") for k in tr.MODE_ORDER}
    pick = tr.choose_column("LC-MS", levels, [])
    assert (pick.text, pick.kind, pick.source) == (text, kind, "rule")


def test_lc_column_notes_chiral_poor_modes_and_a_missing_separation_step():
    pick = tr.choose_column("LC-MS", {m_: "good" if m_ == "RP" else "poor" for m_ in tr.MODE_ORDER}, [], n_stereo=2)
    assert "chiral column only for enantiomer separation" in pick.note
    pick = tr.choose_column("LC-MS", {m_: "poor" for m_ in tr.MODE_ORDER}, [])
    assert pick.text == "C18 (reversed phase)" and "no mode is rated well" in pick.note
    # no Separation step: any recognised published LC column, MS first; none at all -> not available
    assert tr.choose_column("LC-MS", None, [m(GC_COL), m(HILIC_COL, "UV"), m(C18, "MS")]).text == C18
    none = tr.choose_column("LC-MS", None, [])
    assert none.text is None and none.source is None


def test_gc_column_is_a_published_gc_column_only_without_derivatization():
    pub = tr.choose_column("GC-MS", None, [m(C18), m(GC_COL, "FID"), m("GC-5%phenyl, 0.25 mm × 30 m", "MS")], gc_direct=True)
    assert (pub.source, pub.kind, pub.text) == ("published", "GC", "GC-5%phenyl, 0.25 mm × 30 m")      # MS detector first
    assert tr.choose_column("GC-MS", None, [m(GC_COL, "FID")], gc_direct=True).text == GC_COL
    derivatized = tr.choose_column("GC-MS", None, [m(GC_COL, "FID")], gc_direct=False)
    assert (derivatized.text, derivatized.source) == (tr.GC_COLUMN, "rule") and "derivatives" in derivatized.note
    assert tr.choose_column("GC-MS", None, [m(C18)], gc_direct=True).text == tr.GC_COLUMN
    assert "5 % phenyl" in tr.GC_COLUMN and "DB-5ms" in tr.GC_COLUMN


RANGES = [(0.0, 0.0), (4.1, 14.0)]


def test_ph_choice_says_whether_the_ph_is_inside_the_robust_range():
    inside = tr.ph_choice(6.8, True, RANGES)
    assert (inside.ph, inside.inside) == (6.8, True) and inside.text == "pH 6.8 is inside the robust range pH 0 or 4.1–14"
    outside = tr.ph_choice(2.7, False, RANGES)
    assert outside.inside is False and outside.text.startswith("default pH 2.7 is outside the robust range pH 0 or 4.1–14")
    assert "pH 4.1–10 is more robust" in outside.text                    # the robust range inside the practical 2-10 window
    given = tr.ph_choice(2.7, True, RANGES)
    assert given.text.startswith("pH 2.7 is outside") and not given.text.startswith("default")


def test_ph_choice_edge_cases():
    assert tr.ph_choice(2.7, False, None).inside is None and "not available" in tr.ph_choice(2.7, False, None).text
    assert tr.ph_choice(2.7, False, [], n_sites=3).inside is False and "buffer well" in tr.ph_choice(2.7, False, [], n_sites=3).text
    none = tr.ph_choice(2.7, False, [(0.0, 14.0)], n_sites=0)
    assert none.inside is True and "no ionizable site" in none.text
    assert "more robust" not in tr.ph_choice(2.7, False, [(10.8, 14.0), (4.3, 4.4)]).text               # windows narrower than 0.5 or outside 2-10 are not suggested
    assert tr.ph_choice(2.0, False, [(0.0, 2.0)]).inside is True                                    # a range edge counts as inside
    assert tr.ADDITIVE_TEXT["Formic acid"] == "0.1 % formic acid" and set(tr.ADDITIVE_TEXT) == {"Formic acid", "Acetic acid", "Ammonium formate", "Ammonium acetate", "None"}


# ================================================================ rules: main ion, derivatization, detectors, methods
def ion(adduct, score, likelihood="high"):
    return SimpleNamespace(adduct=adduct, score=score, likelihood=likelihood)


def test_main_ion_follows_the_recommended_polarity():
    pos, neg = ion("[M+H]+", 65.0), ion("[M+HCOO]-", 15.0)
    assert tr.pick_main_ion("+", pos, neg)[:2] == ("+", pos)
    assert tr.pick_main_ion("−", pos, neg)[:2] == ("-", neg) and tr.pick_main_ion("-", pos, neg)[:2] == ("-", neg)
    assert tr.pick_main_ion("+", None, neg) == ("+", None, "no ion expected in ESI+")


def test_main_ion_for_both_and_poor_is_the_higher_score_and_says_so():
    sign, best, note = tr.pick_main_ion("both", ion("[M+Na]+", 63.8), ion("[M+HCOO]-", 70.0))
    assert (sign, best.adduct) == ("-", "[M+HCOO]-") and "both polarities respond" in note and "ESI−" in note and "70 vs 63.8" in note
    assert tr.pick_main_ion("both", ion("[M+H]+", 70.0), ion("[M-H]-", 70.0))[0] == "+"            # a tie goes to ESI+
    sign, best, note = tr.pick_main_ion("poor", ion("[M]+", 5.0), None)
    assert (sign, best.adduct) == ("+", "[M]+") and "weak in both" in note
    assert tr.pick_main_ion("both", None, None)[:2] == (None, None)
    assert tr.pick_main_ion(None, ion("[M+H]+", 1.0), ion("[M-H]-", 2.0))[0] == "-"                  # an unknown polarity is treated as both


def test_lc_derivatization_is_needed_only_when_the_esi_response_is_weak():
    channels = [("amine / phenol", "Dansyl chloride"), ("carboxyl", "DmPA bromide")]
    for level in ("strong", "moderate"):
        assert tr.lc_derivatization(level, channels).kind == "not needed" and tr.lc_derivatization(level, channels).text == "Not needed"
    weak = tr.lc_derivatization("weak", channels)
    assert weak.kind == "suggested" and weak.text == "CIL: Dansyl chloride (amine / phenol)" and "also carboxyl" in weak.basis
    none = tr.lc_derivatization("weak", [])
    assert none.kind == "none" and none.text == "None available"
    assert tr.lc_derivatization(None, channels) is None


def test_gc_derivatization_from_the_gc_verdict():
    assert tr.gc_derivatization("yes", "TMS").text == "Not needed"
    needed = tr.gc_derivatization("after derivatization", "MeOx + TMS", "needs derivatization (MeOx + TMS, 5 TMS → 569.29 Da)")
    assert (needed.kind, needed.text) == ("needed", "Needed: MeOx + TMS") and "569.29" in needed.basis
    assert tr.gc_derivatization("no", "TMS").kind == "not feasible" and tr.gc_derivatization(None, "TMS") is None


def test_best_non_ms_detector():
    rows = [["UV/PDA", "poor", "none"], ["FLD", "possible", "x"], ["CAD", "good", "y"], ["ELSD", "good", "y"], ["RI", "possible", "z"]]
    assert tr.best_non_ms(rows) == ("CAD", "good", "y")                                              # ties go to the earlier detector
    assert tr.best_non_ms([["UV/PDA", "good", "a"], ["CAD", "good", "b"]])[0] == "UV/PDA"
    assert tr.best_non_ms([["UV/PDA", "poor", "a"], ["RI", "possible", "z"]])[0] == "RI"
    assert tr.best_non_ms([]) is None


def test_methods_count_text():
    assert tr.methods_count_text(1, 150) == "151 exact (JP18 1, RepoRT 150)"
    assert tr.methods_count_text(0, 6) == "6 exact (RepoRT 6)" and tr.methods_count_text(3, 0) == "3 exact (JP18 3)"
    assert tr.methods_count_text(0, 0, n_similar=10) == "None for this compound; methods of 10 similar compounds"
    assert tr.methods_count_text(0, 0, n_similar=1).endswith("1 similar compound") and tr.methods_count_text(0, 0) == "None for this compound"


def test_silyl_fragment_is_the_usually_stronger_ion_of_the_derivative():
    # glycine 3 TMS: M = 291.1506, M+• 291.1501 (MassBank MSBNK-RIKEN-PR010014 has the base-peak family at 73 / 174 and a 276 peak)
    name, mz = tr.silyl_fragment("TMS", 291.15061)
    assert name == "M−15" and abs(mz - 276.1266) < 0.0002
    assert tr.silyl_fragment("MeOx + TMS", 569.28750)[0] == "M−15" and abs(tr.silyl_fragment("MeOx + TMS", 569.28750)[1] - 554.2635) < 0.0002
    name, mz = tr.silyl_fragment("TBDMS", 400.0)                                                   # tert-butyl is lost: C4H9 = 57.0704
    assert name == "M−57" and abs(mz - (400.0 - 57.0704253 - tr.ELECTRON)) < 1e-9
    assert tr.silyl_fragment("TFAA", 300.0) is None and tr.silyl_fragment("TMS", None) is None


def test_gc_methods_text_counts_gc_methods_and_names_the_lc_ones_apart():
    assert tr.gc_methods_text(1, 150) == "1 GC exact (150 LC methods also exist)" and tr.gc_methods_text(2, 1) == "2 GC exact (1 LC method also exists)"
    assert tr.gc_methods_text(0, 90) == "None for GC (90 LC methods exist)" and tr.gc_methods_text(0, 1) == "None for GC (1 LC method exists)"
    assert tr.gc_methods_text(0, 0) == "None for GC" and tr.gc_methods_text(3, 0) == "3 GC exact"


def test_method_counts_are_the_true_totals_of_the_search_result_not_the_listed_rows():
    cut = SimpleNamespace(exact=[SimpleNamespace(source="JP18")] + [SimpleNamespace(source="RepoRT")] * 25, exact_totals={"JP18": 1, "RepoRT": 150})
    assert tm.method_counts(cut) == (1, 150)                                                       # 26 rows listed, 151 methods exist
    assert tm.method_counts(SimpleNamespace(exact=[], exact_totals={"JP18": 0, "RepoRT": 0})) == (0, 0)
    assert tm.method_counts(SimpleNamespace(exact=[], exact_totals={})) == (0, 0)                 # a source that was not searched is absent


def test_the_statement_says_what_the_sheet_says():
    parts = {"sign": "+", "mz": 100.1234, "ion_label": "[M+H]+", "column_kind": "C18", "column_source": "rule", "ph": 2.7, "deriv_kind": "none",
             "deriv_text": "None available", "methods": 0, "methods_gc": False}
    # nothing suitable: no 'by None suitable by ESI or GC', no ion and no column for a compound that has no route
    none = tm.statement("X", tr.decide("Any", "weak", "no"), parts)
    assert none.startswith("X: unlikely to be measured by LC-MS or GC-MS; ") and "None suitable" not in none and "m/z" not in none and "column" not in none
    # HPLC with a non-MS detector: no ESI ion, the LC column and pH stay
    nonms = tm.statement("X", tr.decide("Any", "weak", "no", ("CAD", "good")), parts)
    assert nonms.startswith("X: possible by HPLC-CAD (no MS) on a C18 column at pH 2.7;") and "m/z" not in nonms and "derivatization" not in nonms
    # an LC-MS or GC-MS route keeps the ion, and the GC column is the one the sheet gives
    assert "ESI+ [M+H]+ m/z 100.1234 on a C18 column at pH 2.7" in tm.statement("X", tr.decide("Any", "strong", "no"), parts)
    gc = {**parts, "ion_label": "M+•", "column_kind": "GC", "methods": 2, "methods_gc": True}
    assert "EI M+• m/z 100.1234 on a DB-5ms-type column" in tm.statement("X", tr.decide("Any", "weak", "yes"), gc)
    published = tm.statement("X", tr.decide("Any", "weak", "yes"), {**gc, "column_source": "published"})
    assert "on the published GC column" in published and "DB-5ms" not in published and published.endswith("; 2 published GC methods.")
    assert tm.statement("X", tr.decide("Any", "weak", "yes"), {**gc, "methods": 0}).endswith("; no published GC method.")
    assert tm.statement("X", tr.decide("Any", "strong", "yes"), {**parts, "methods": 1}).endswith("; 1 published method.")


# ================================================================ the sheet from pieces (no engines)
def test_an_empty_sheet_says_not_available_in_every_row():
    q = tm.Query(smiles="CCO")
    decision = tm._decide(q, {})
    assert decision.verdict is None
    rows, parts = tm.build_sheet(q, {}, decision)
    assert [r[0] for r in rows] == SHEET_ITEMS and all(len(r) == 3 for r in rows)
    na = {r[0] for r in rows if text(r[1]) == "not available"}
    assert {"Can it be measured", "Column", "Polarity", "Main ion", "Theoretical m/z", "Derivatization", "Non-MS alternative", "Published methods"} <= na
    assert "incomplete" in tm.statement("Ethanol", decision, parts)


# ================================================================ end to end
def assert_sheet_shape(result):
    assert len(result["blocks"]) == 1
    block = result["blocks"][0]
    assert (block["type"], block["id"], block["title"]) == ("table", "method_sheet", "Method sheet") and block["columns"] == ["Item", "Recommendation", "Basis"]
    rows = block["data"]["rows"]
    assert [r[0] for r in rows] == SHEET_ITEMS and all(len(r) == 3 and r[1] and r[2] for r in rows)
    assert result["answer"]["level"] is None and [f[0] for f in result["answer"]["facts"]] == ["Name", "Formula", "Monoisotopic mass", "InChIKey"]


def assert_no_transitions(*parts):
    blob = json.dumps(parts, ensure_ascii=False)
    assert not FORBIDDEN.findall(blob), FORBIDDEN.findall(blob)


@needs_molgpka
def test_caffeine_is_measurable_by_lc_ms_in_esi_positive(offline_module):
    result = cached(CAFFEINE)
    assert_sheet_shape(result)
    assert statuses(result) == {n: "done" for n in range(1, 10)} and result["warnings"] == []
    s = sheet(result)
    assert s["Can it be measured"][0] == "Yes" and s["Platform"][0].startswith("LC-MS")
    assert s["Polarity"][0] == "ESI+" and s["Main ion"][0] == "[M+H]+"
    assert abs(float(s["Theoretical m/z"][0]) - 195.0877) <= 0.0005 and "C8H11N4O2" in s["Theoretical m/z"][1] and "agrees" in s["Theoretical m/z"][1]
    assert s["Derivatization"][0] == "Not needed" and s["Mobile phase / pH"][0] == "0.1 % formic acid, pH 2.7"
    assert s["Published methods"][0].split(" ")[0].rstrip("+").isdigit() and int(s["Published methods"][0].split(" ")[0].rstrip("+")) >= 1
    assert s["Non-MS alternative"][0].startswith("UV/PDA at ~")
    assert s["Column"][0] not in ("", "not available") and all("(step" in basis for _rec, basis in s.values())
    answer = result["answer"]
    assert answer["statement"].startswith("Caffeine: measurable by LC-MS") and "ESI+ [M+H]+ m/z 195.0877" in answer["statement"]
    assert "pH 2.7" in answer["statement"] and "no derivatization" in answer["statement"] and "published methods." in answer["statement"]
    assert dict(answer["facts"]) == {"Name": "Caffeine", "Formula": "C8H10N4O2", "Monoisotopic mass": "194.0804", "InChIKey": "RYYVLZVUVIJVGH-UHFFFAOYSA-N"}
    assert_no_transitions(result)


@needs_molgpka
def test_caffeine_main_ion_and_polarity_agree_with_the_standalone_pages(offline_module):
    s = sheet(cached(CAFFEINE))
    inputs = {"structure": {"mode": "SMILES / InChI", "value": CAFFEINE}}
    base = adducts_runner.run_adducts({**inputs, "ion_mode": "Positive", "additive": "Formic acid"})["outputs"]["base"]["values"]
    assert base["Dominant ion"]["value"] == s["Main ion"][0] and base["m/z"]["value"] == f"{float(s['Theoretical m/z'][0]):.4f}"
    chips = {c["label"]: c for c in ionization_runner.run_ionization({**inputs, "ph": 2.7})["outputs"]["chips"]["chips"]}
    assert f"ESI{chips['Polarity']['value']}" == s["Polarity"][0]
    assert f"ESI+ {chips['ESI+']['value']}, ESI− {chips['ESI−']['value']}" in s["Polarity"][1]


@needs_molgpka
@needs_report
def test_published_retention_gives_caffeine_a_reversed_phase_column_and_glucose_a_hilic_column(offline_module):
    caffeine, glucose = cached(CAFFEINE), cached(GLUCOSE)
    assert tr.column_class(sheet(caffeine)["Column"][0]) == "RP" and "RT" in sheet(caffeine)["Column"][1] and "C18 column" in caffeine["answer"]["statement"]
    assert tr.column_class(sheet(glucose)["Column"][0]) == "HILIC" and "HILIC column" in glucose["answer"]["statement"]
    # the published-method count is the real number behind the cut lists, not the 25 listed per source
    assert int(sheet(caffeine)["Published methods"][0].split(" ")[0]) > 26 and "RepoRT" in sheet(caffeine)["Published methods"][0]


@needs_molgpka
def test_glucose_main_ion_agrees_with_the_standalone_adducts_page_in_negative_mode(offline_module):
    s = sheet(cached(GLUCOSE))
    assert s["Platform"][0].startswith("LC-MS") and s["Polarity"][0] == "ESI−" and "both respond" in s["Polarity"][1]
    base = adducts_runner.run_adducts({"structure": {"mode": "SMILES / InChI", "value": GLUCOSE}, "ion_mode": "Negative", "additive": "Formic acid"})["outputs"]["base"]["values"]
    assert base["Dominant ion"]["value"] == s["Main ion"][0] == "[M+HCOO]-"
    assert base["m/z"]["value"] == f"{float(s['Theoretical m/z'][0]):.4f}" == "225.0616"
    # LC-MS needs no derivatization; the GC route would
    assert s["Derivatization"][0] == "Not needed" and "GC-MS would need MeOx + TMS" in s["Derivatization"][1]


def test_a_volatile_small_compound_goes_to_gc_ms(offline_module):
    result = cached(TOLUENE)
    assert_sheet_shape(result)
    s = sheet(result)
    st = statuses(result)
    assert s["Can it be measured"][0] == "Yes" and s["Platform"][0] == "GC-MS"
    assert st[5] == "skipped" and st[9] == "done" and st[7] == "done" and st[2] == "done"
    assert s["Derivatization"][0] == "Not needed" and s["Mobile phase / pH"][0] == "Not applicable (GC)"
    assert s["Polarity"][0].startswith("EI") and s["Main ion"][0] == "M+•"
    assert abs(float(s["Theoretical m/z"][0]) - 92.0621) <= 0.0005                                      # C7H8 monoisotopic 92.0626 minus an electron
    assert "DB-5ms" in s["Column"][0] or s["Column"][0].startswith("GC")
    assert result["answer"]["statement"].startswith("Toluene: measurable by GC-MS") and "EI M+• m/z 92.0621" in result["answer"]["statement"]
    assert_no_transitions(result)


def test_glucose_on_gc_ms_needs_meox_tms_derivatization(offline_module):
    result = cached(GLUCOSE, instrument="GC-MS")
    assert_sheet_shape(result)
    s = sheet(result)
    st = statuses(result)
    assert (st[3], st[5]) == ("skipped", "skipped") and st[9] == "done"
    assert s["Can it be measured"][0] == "Possible" and s["Platform"][0] == "GC-MS (derivatized)"
    assert s["Derivatization"][0].startswith("Needed: MeOx + TMS") and "C22H55NO6Si5" in s["Derivatization"][0]
    assert s["Column"][0] == tr.GC_COLUMN and "same column" in s["Column"][1]                          # a derivative runs on the default phase, not the published FID column
    assert s["Main ion"][0] == "M+• (MeOx + TMS derivative)" and abs(float(s["Theoretical m/z"][0]) - 569.2870) <= 0.0005
    assert "C22H55NO6Si5+•" in s["Theoretical m/z"][1]
    assert result["answer"]["statement"].startswith("D-Glucose: possible by GC-MS (derivatized), EI M+• m/z 569.2870")
    assert "derivatization: MeOx + TMS" in result["answer"]["statement"]
    assert_no_transitions(result)


def test_the_gc_main_ion_names_the_stronger_fragment_of_the_derivative(offline_module):
    s = sheet(cached(GLUCOSE, instrument="GC-MS"))
    assert "[M−15]+ 554.2635 is usually stronger" in s["Main ion"][1] and s["Main ion"][1].endswith("(step 7)")      # 569.2870 - CH3 (15.0235)
    assert abs(float(s["Theoretical m/z"][0]) - 15.0235 - 554.2635) <= 0.0005
    glycine = sheet(run("NCC(=O)O", instrument="GC-MS"))                                          # the MassBank example of the EI MS Annotation workflow: 3 TMS, 276 present
    assert "[M−15]+ 276.1266" in glycine["Main ion"][1] and abs(float(glycine["Theoretical m/z"][0]) - 291.1501) <= 0.0005
    assert "usually stronger" not in sheet(cached(TOLUENE))["Main ion"][1]                       # no derivative, no fragment


def test_on_the_gc_route_the_published_methods_are_the_gc_ones(offline_module):
    glucose = cached(GLUCOSE, instrument="GC-MS")
    cell, basis = sheet(glucose)["Published methods"]
    assert re.fullmatch(r"None for GC \(\d+ LC methods exist\)", cell) and "exact matches are LC methods" in basis      # glucose is in JP18 / RepoRT, but only with LC columns
    assert glucose["answer"]["statement"].endswith("; no published GC method.")
    assert sheet(cached(TOLUENE))["Published methods"] == ("None for GC", "no exact match (step 8)")
    # caffeine has a GC method in JP18: it is counted, named and its column is the one the sheet and the statement give
    caffeine = cached(CAFFEINE, instrument="GC-MS")
    s = sheet(caffeine)
    assert re.match(r"[1-9]\d* GC exact", s["Published methods"][0]) and s["Published methods"][1].startswith("top: JP18")
    assert tr.column_class(s["Column"][0]) == "GC" and s["Column"][0] != tr.GC_COLUMN
    assert "on the published GC column" in caffeine["answer"]["statement"] and "DB-5ms" not in caffeine["answer"]["statement"]
    assert re.search(r"; [1-9]\d* published GC methods?\.$", caffeine["answer"]["statement"])
    # the LC route still counts every method
    assert re.match(r"[1-9]\d* exact \(JP18 1", sheet(cached(CAFFEINE))["Published methods"][0])


def test_the_gc_column_basis_cites_the_step_it_rests_on(offline_module):
    assert sheet(cached(TOLUENE))["Column"][1] == "no published GC column; default phase (steps 2, 8)"           # a volatile compound: the GC verdict of step 2
    assert sheet(cached(GLUCOSE, instrument="GC-MS"))["Column"][1].endswith("(steps 7, 8)")                       # a derivative: the derivatization of step 7


def test_the_instrument_filter_judges_only_the_chosen_platform(offline_module):
    lc = cached(TOLUENE, instrument="LC-MS")
    s = sheet(lc)
    assert s["Can it be measured"][0] == "Unlikely" and "APCI" in s["Can it be measured"][1] and "GC-MS would be volatile as is" in s["Can it be measured"][1]
    assert s["Platform"][0] == "LC-MS (ESI weak)" and statuses(lc)[5] == "done" and "weak" in s["Polarity"][0] + s["Polarity"][1]
    assert lc["answer"]["statement"].startswith("Toluene: unlikely to be measured by LC-MS")
    # the derivatization row still answers for LC-MS: weak ESI and no labelable group
    assert s["Derivatization"][0] == "None available"
    gc = cached(CAFFEINE, instrument="GC-MS")
    assert sheet(gc)["Can it be measured"][0] == "Yes" and sheet(gc)["Platform"][0] == "GC-MS"


@needs_molgpka
def test_a_given_ph_and_additive_are_used_and_the_inside_or_outside_statement_follows(offline):
    result = run(CAFFEINE, ph="6.8", additive="Ammonium acetate")
    s = sheet(result)
    assert s["Mobile phase / pH"][0] == "5-10 mM ammonium acetate, pH 6.8"
    assert s["Mobile phase / pH"][1].startswith("pH 6.8 is inside the robust range") and not s["Mobile phase / pH"][1].startswith("default")
    assert result["warnings"] == []                                                          # a pH was given: nothing to warn about
    blank = run(CAFFEINE, additive="Ammonium acetate")
    assert sheet(blank)["Mobile phase / pH"][0].endswith("pH 2.7") and sheet(blank)["Mobile phase / pH"][1].startswith("default pH 2.7")
    assert any("No pH was given" in w and "about pH 6.8" in w for w in blank["warnings"])


# ---------------------------------------------------------------- identity
def test_identity_online_lists_the_names_and_database_ids(monkeypatch):
    seen = {}

    def fake(query, kind="Auto", **kw):
        seen.update(query=query, kind=kind, **kw)
        return FakeResolved()
    monkeypatch.setattr(resolve, "resolve", fake)
    status, note, blocks, ident = tm.step_identity(tm.Query(smiles=CAFFEINE, name="caffeine"))
    assert status == "done" and (seen["kind"], seen["query"]) == ("InChIKey", "RYYVLZVUVIJVGH-UHFFFAOYSA-N") and seen["timeout"] == tm.ID_TIMEOUT
    assert note == "Caffeine · C8H10N4O2" and ident["online"] and ident["name"] == "Caffeine"
    rows = {r[0]: r[1] for r in blocks[0]["data"]["rows"]}
    assert (blocks[0]["id"], blocks[0]["columns"]) == ("identity", ["Item", "Value"])
    assert rows["Formula"] == "C8H10N4O2" and rows["Monoisotopic mass"] == "194.0804" and rows["InChIKey"] == "RYYVLZVUVIJVGH-UHFFFAOYSA-N"
    assert rows["PubChem CID"] == {"text": "2519", "url": "https://pubchem.ncbi.nlm.nih.gov/compound/2519"}
    assert rows["HMDB"] == {"text": "HMDB0001847", "url": "https://hmdb.ca/metabolites/HMDB0001847"} and rows["KEGG"] == "C07481" and rows["CAS"] == "58-08-2"


def test_identity_offline_falls_back_to_rdkit_values_and_the_workflow_goes_on(offline):
    status, note, blocks, ident = tm.step_identity(tm.Query(smiles=CAFFEINE, name="caffeine"))
    assert status == "done" and "RDKit values only (PubChem unreachable)" in note and not ident["online"]
    assert (ident["name"], ident["formula"], ident["inchikey"]) == ("caffeine", "C8H10N4O2", "RYYVLZVUVIJVGH-UHFFFAOYSA-N") and abs(ident["mass"] - 194.0804) < 1e-4
    assert {r[0] for r in blocks[0]["data"]["rows"]} == {"Name", "Formula", "Monoisotopic mass", "InChIKey"}          # no database rows without the lookup
    status, note, _blocks, nameless = tm.step_identity(tm.Query(smiles=CAFFEINE))
    assert nameless["name"] == "C8H10N4O2" and note.startswith("C8H10N4O2 · RDKit values only")      # no name typed and none found: the formula names it, once


def test_identity_looks_the_compound_up_by_inchikey_so_metformin_is_not_an_isotopologue(monkeypatch):
    """A SMILES search answered '[14C]metformin' for plain metformin (the first of several PubChem records); the InChIKey finds the parent compound."""
    seen = []

    def fake(query, kind="Auto", **kw):
        seen.append(kind)
        return FakeResolved(name="Metformin" if kind == "InChIKey" else "[14C]metformin", cid=4091 if kind == "InChIKey" else 152743144)
    monkeypatch.setattr(resolve, "resolve", fake)
    _status, note, _blocks, ident = tm.step_identity(tm.Query(smiles="CN(C)C(=N)NC(=N)N"))
    assert ident["name"] == "Metformin" and ident["cid"] == 4091 and note.startswith("Metformin") and seen == ["InChIKey"]


def test_identity_falls_back_to_the_smiles_search_only_when_pubchem_does_not_know_the_key(monkeypatch):
    seen = []

    def not_found(query, kind="Auto", **kw):
        seen.append(kind)
        if kind == "InChIKey":
            raise resolve.ResolveError(f"'{query}' not found in PubChem.")
        return FakeResolved(name="Found by SMILES")
    monkeypatch.setattr(resolve, "resolve", not_found)
    assert tm.step_identity(tm.Query(smiles=CAFFEINE))[3]["name"] == "Found by SMILES" and seen == ["InChIKey", "SMILES"]
    seen.clear()

    def unreachable(query, kind="Auto", **kw):
        seen.append(kind)
        raise resolve.ResolveError("Network lookup failed (PubChem/UniChem unreachable): timed out")
    monkeypatch.setattr(resolve, "resolve", unreachable)
    status, note, _blocks, ident = tm.step_identity(tm.Query(smiles=CAFFEINE))
    assert status == "done" and not ident["online"] and "PubChem unreachable" in note and seen == ["InChIKey"]               # one failed request, not two


def test_identity_for_a_structure_pubchem_does_not_know_keeps_the_rdkit_values(monkeypatch):
    monkeypatch.setattr(resolve, "resolve", lambda q, kind="Auto", **kw: FakeResolved(name=None, cid=None, ids={}, links={}))
    status, note, _blocks, ident = tm.step_identity(tm.Query(smiles="CCO"))
    assert status == "done" and "not in PubChem" in note and ident["formula"] == "C2H6O" and ident["name"] == "C2H6O"


def test_the_platform_step_keeps_an_odd_structure_that_does_not_survive_a_smiles_round_trip():
    status, note, blocks, value = tm.step_platform(tm.Query(smiles="[Fe+2].c1cc[cH-]c1.c1cc[cH-]c1"))             # ferrocene written as ions
    assert status == "done" and blocks[0]["id"] == "platform" and value["gc"]["gc"] in ("yes", "after derivatization", "no")


def test_warnings_drop_preparation_notes_and_repeated_remarks():
    steps = tm.StepRunner(tm.STEP_TITLES)
    got = {3: {"warnings": ["Permanent charge (net +1, e.g. quaternary N+): counted as fixed.", "2 disconnected parts (salt or mixture): only the largest is analysed."]},
           4: {"warnings": ["Input neutralized (charges that can be removed by adding or losing H).", "2 disconnected parts (salt or mixture): only the largest part (9 atoms) is analysed."]},
           5: {"warnings": ["m/z not computed for some ions of C4H12ClN5."]}}
    q = tm.Query(smiles="CCO", warnings=["'CO' is both a valid molecular formula and a valid SMILES"])
    out = tm._warnings(q, steps, got, tr.decide("Any", "strong", "yes"))
    assert out == ["'CO' is both a valid molecular formula and a valid SMILES", "2 disconnected parts (salt or mixture): only the largest is analysed.",
                   "m/z not computed for some ions of C4H12ClN5."]


def _warn(q, got, decision):
    return tm._warnings(q, tm.StepRunner(tm.STEP_TITLES), got, decision)


def test_neims_remarks_belong_to_the_gc_route_only():
    got = {4: {"warnings": ["The structure has a net charge of +1. NEIMS predicts EI of neutral molecules, so this prediction is unreliable.", "Another remark."]}}
    q = tm.Query(smiles="CCO")
    assert _warn(q, got, tr.decide("Any", "strong", "no")) == ["Another remark."]                         # LC-MS: EI is not used
    assert len(_warn(q, got, tr.decide("Any", "weak", "yes"))) == 2                                          # GC-MS: the EI prediction is the main ion


def test_choline_on_lc_ms_has_no_ei_remark(offline_module):
    result = run("C[N+](C)(C)CCO")
    assert sheet(result)["Main ion"][0] == "[M]+" and abs(float(sheet(result)["Theoretical m/z"][0]) - 104.1070) <= 0.0005       # C5H14NO+ minus an electron
    assert not any("NEIMS" in w for w in result["warnings"])


def test_a_structure_without_carbon_is_flagged(offline_module):
    result = run("[Na+].[Cl-]")                                                                              # NaCl is analysed as HCl: not a sensible GC-MS result
    assert any("contains no carbon" in w and "(HCl)" in w and "ion chromatography" in w for w in result["warnings"])
    assert not any("contains no carbon" in w for w in run("CCO")["warnings"])
    assert not any("contains no carbon" in w for w in run("CC(=O)[O-].[Na+]")["warnings"])                     # sodium acetate is analysed as acetic acid


def test_isotope_labels_are_flagged_on_the_lc_route_only(offline_module):
    d3 = "[2H]C([2H])([2H])n1cnc2c1c(=O)n(C)c(=O)n2C"
    result = run(d3)
    assert sheet(result)["Main ion"][0] == "[M+H]+" and abs(float(sheet(result)["Theoretical m/z"][0]) - 195.0877) <= 0.0005       # the labels do not reach the ESI m/z
    assert any("isotope labels" in w and "unlabelled" in w for w in result["warnings"])
    assert not any("isotope labels" in w for w in run(CAFFEINE)["warnings"])
    assert not any("isotope labels" in w for w in run(d3, instrument="GC-MS")["warnings"])                  # the EI M+• is calculated from the labelled mass


@pytest.mark.parametrize("ph, additive, instrument, expect", [
    ("7", "Formic acid", "LC-MS", "pH 7 does not match formic acid, which gives about pH 2.7."),
    ("2.7", "Ammonium acetate", "LC-MS", "pH 2.7 does not match ammonium acetate, which gives about pH 6.8."),
    ("", "None", "LC-MS", "No pH was given, so pH 2.7 was used for the separation and ionization rules, but no additive (unbuffered) gives about pH 7."),
    ("", "Ammonium formate", "LC-MS", "but ammonium formate gives about pH 6.5. Enter the pH to match.")])
def test_a_ph_that_does_not_fit_the_additive_is_a_warning_whether_or_not_it_was_typed(ph, additive, instrument, expect):
    q = tm.parse_inputs(form("CCO", ph=ph, additive=additive, instrument=instrument))
    out = _warn(q, {}, tr.decide("Any", "strong", "yes"))
    assert len(out) == 1 and expect in out[0] and "none" not in out[0].lower().split("but")[-1]


def test_a_ph_that_fits_the_additive_or_a_gc_route_gives_no_ph_warning():
    ok = tm.parse_inputs(form("CCO", ph="6.8", additive="Ammonium acetate"))
    assert _warn(ok, {}, tr.decide("Any", "strong", "yes")) == []
    assert _warn(tm.parse_inputs(form("CCO")), {}, tr.decide("Any", "strong", "yes")) == []                    # blank pH with the default additive
    assert _warn(tm.parse_inputs(form("CCO", ph="7")), {}, tr.decide("Any", "weak", "yes")) == []              # GC-MS: no mobile phase


# ---------------------------------------------------------------- failing and missing engines
def boom(*args, **kwargs):
    raise RuntimeError("engine exploded")


@pytest.mark.parametrize("module_path, target, step, rows_na", [
    ("mzdesk.platform_fit.non_ms_rules", "analyze", 6, ["Non-MS alternative"]),
    ("mzdesk.methods.known_methods", "search", 8, ["Published methods"]),
    ("mzdesk.platform_fit.adduct_rules", "predict", 5, ["Main ion", "Theoretical m/z"]),
    ("mzdesk.platform_fit.ionization_rules", "assess", 4, [])])
@needs_molgpka
def test_a_failing_engine_fails_its_step_and_the_sheet_still_comes_back(offline, monkeypatch, module_path, target, step, rows_na):
    import importlib
    monkeypatch.setattr(importlib.import_module(module_path), target, boom)
    result = run(CAFFEINE)
    assert_sheet_shape(result)
    st = statuses(result)
    assert st[step] == "failed"
    assert "engine exploded" in next(s["note"] for s in result["steps"] if s["id"] == step)
    assert any(w.startswith(f"Step {step} (") and "engine exploded" in w for w in result["warnings"])
    s = sheet(result)
    for item in rows_na:
        assert s[item][0] == "not available", item
    assert all(status in ("done", "skipped") for n, status in st.items() if n != step)          # every other step still ran
    assert result["answer"]["statement"] and not result["answer"]["statement"].startswith("The method sheet")
    if step == 8:                                                                              # no published methods: the column falls back to the rule map
        assert s["Column"][1].startswith("rule:")
    assert_no_transitions(result)


def test_missing_molgpka_weights_skip_separation_with_a_note_and_the_sheet_goes_on(offline, monkeypatch, tmp_path):
    monkeypatch.setenv("MZDESK_MOLGPKA_WEIGHTS", str(tmp_path / "missing.npz"))
    result = run(CAFFEINE)
    assert_sheet_shape(result)
    step3 = result["steps"][2]
    assert step3["status"] == "skipped" and "MolGpKa weights not installed" in step3["note"]
    s = sheet(result)
    assert "robust range is not available" in s["Mobile phase / pH"][1]
    assert s["Column"][0] != "not available" and "published:" in s["Column"][1]               # no separation rules: the first recognised published LC column
    assert s["Polarity"][0] == "ESI+" and s["Main ion"][0] == "[M+H]+"                         # ionization falls back to textbook group pKa values
    assert not result["warnings"] or all("failed" not in w for w in result["warnings"])


def test_an_isotope_check_that_disagrees_is_a_warning(offline, monkeypatch):
    from mzdesk.ms1 import isotopes
    real = isotopes.ion_mz
    monkeypatch.setattr(isotopes, "ion_mz", lambda formula, adduct: (real(formula, adduct)[0] + 0.01, real(formula, adduct)[1]))
    result = run("CCO", instrument="GC-MS")
    s = sheet(result)
    assert "differs by +0.0100" in s["Theoretical m/z"][1]
    assert any("disagree on the m/z of M+•" in w for w in result["warnings"])


# ---------------------------------------------------------------- the job API and the page
def wait_done(job_id, timeout=120):
    end = time.time() + timeout
    while time.time() < end:
        body = client.get(f"/api/pipelines/jobs/{job_id}").json()
        if body["status"] not in ("running", "queued"):
            return body
        time.sleep(0.05)
    raise AssertionError("job did not finish")


def test_the_spec_has_two_examples_that_resolve_offline():
    assert [e["label"] for e in TARGETED["examples"]] == ["Caffeine, any platform", "D-Glucose, GC-MS"]
    ids = {f["id"] for f in TARGETED["inputs"] if f.get("id")}
    for example in TARGETED["examples"]:
        assert set(example["values"]) <= ids
        value = example["values"]["structure"]
        assert value["mode"] == "Name / ID / Formula" and value["value"] in spec.EXAMPLE_STRUCTURES
        q = tm.parse_inputs(example["values"])                                                 # tests/conftest.py resolves the names offline
        assert q.smiles == spec.EXAMPLE_STRUCTURES[value["value"]] and q.name == value["value"]


@needs_molgpka
@pytest.mark.parametrize("index, platform, ion", [(0, "LC-MS", "[M+H]+"), (1, "GC-MS", "M+•")])
def test_both_examples_run_through_the_job_api_with_live_step_blocks(offline, index, platform, ion):
    example = TARGETED["examples"][index]["values"]
    r = client.post("/api/pipelines/targeted_method/run", json={"ph": "", **example})
    assert r.status_code == 200
    body = wait_done(r.json()["job_id"])
    assert body["status"] == "done" and body["error"] is None
    result = body["result"]
    assert result["pipeline"] == "targeted_method" and set(result) == {"pipeline", "steps", "answer", "blocks", "warnings"}
    assert_sheet_shape(result)
    s = sheet(result)
    assert s["Platform"][0].startswith(platform) and s["Main ion"][0].startswith(ion)
    # every step reports its blocks while the run goes on; every block has an id, a type and data, and the ids are unique across the steps
    blocks = [b for step in body["steps"] for b in step["blocks"]]
    assert all(b["id"] and b["type"] and isinstance(b["data"], dict) for b in blocks)
    assert len({b["id"] for b in blocks}) == len(blocks)
    by_step = {step["id"]: [b["id"] for b in step["blocks"]] for step in body["steps"]}
    assert by_step[1] == ["identity"] and by_step[2] == ["platform"] and by_step[6] == ["non_ms_detectors"] and by_step[9] == ["main_ion"]
    assert by_step[4] == ["ionization_chips", "ionization_notes"] and "derivatization" in by_step[7] and by_step[8][:2] == ["methods_summary", "methods_exact"]
    platform_block = next(b for b in blocks if b["id"] == "platform")
    assert (platform_block["type"], platform_block["title"]) == ("portal", "Platform") and [t["label"] for t in platform_block["data"]["tiles"]] == [
        "Lipidomics", "Metabolomics", "CIL LC-MS", "GC-MS"]
    assert_no_transitions(body)
    json.dumps(body)                                                                           # the whole job is plain JSON


def test_progress_reports_every_step_once_per_state_in_order(offline_module):
    events = []
    run(TOLUENE, progress=lambda *a: events.append(a))
    started = [e[0] for e in events if e[1] == "running"]
    finished = [e[0] for e in events if e[1] in ("done", "skipped", "failed")]
    assert started == [1, 2, 3, 4, 6, 7, 8, 9] and finished == [1, 2, 3, 4, 5, 6, 7, 8, 9]       # step 5 is skipped without running
    assert all(e[3] is not None for e in events if e[1] in ("done", "skipped"))
