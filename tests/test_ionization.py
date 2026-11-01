"""Ionization page tests: rule branches, thresholds, functional-group fallback, M+• levels, page outputs, and a validation set of
compounds with chemistry-based expectations.

The validation expectations come from general LC-MS / GC-MS knowledge (a one-line rationale each), not from what the rules return.
Genuine misses are kept as strict xfail with the reason, not tuned away: if a later change fixes one, the test turns into an
unexpected pass and has to be moved to the hit list.

Tests that need the MolGpKa weights (python -m mzdesk.platform_fit.pka_setup) or the NEIMS weights (python -m mzdesk.ms2.neims_setup)
skip when they are absent. The rule, fallback and API tests do not need them.
"""
from functools import lru_cache
from pathlib import Path

import pytest
from fastapi.testclient import TestClient
from rdkit import Chem
from rdkit.Chem import Crippen

from mzdesk.api.app import app
from mzdesk.ms2 import neims_engine as ne
from mzdesk.platform_fit import ionization_rules as R
from mzdesk.platform_fit import pka_engine as pe
from mzdesk.spec import MODULES

client = TestClient(app)
needs_pka = pytest.mark.skipif(not Path(pe.weights_path()).is_file(), reason="MolGpKa weights not installed (python -m mzdesk.platform_fit.pka_setup)")
needs_neims = pytest.mark.skipif(not Path(ne.weights_path()).is_file(), reason="NEIMS weights not installed (python -m mzdesk.ms2.neims_setup)")
MINUS = R.MINUS
LIDOCAINE = "CCN(CC)CC(=O)Nc1c(C)cccc1C"
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
GLUCOSE = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"
PC = "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCCCC"


@lru_cache(maxsize=None)
def A(smiles, ph=2.7):
    return R.assess(smiles, ph)


def _run(smiles, ph=2.7):
    return client.post("/api/modules/ionization/run", json={"structure": {"mode": "SMILES / InChI", "value": smiles}, "ph": ph})


def _chips(smiles, ph=2.7):
    """label -> chip of the page outputs."""
    return {c["label"]: c for c in _run(smiles, ph).json()["outputs"]["chips"]["chips"]}


def _help(item):
    """Help of a chip / card as a list of lines (the page sends a string for one line)."""
    h = item.get("help") or []
    return [h] if isinstance(h, str) else h


def site(kind, pka, atom=0, element="O"):
    return pe.Site(atom=atom, element=element, kind=kind, pka=pka)


def prof(smiles, sites=(), fixed=0, logp=None):
    """A Profile with hand-made sites: the rules can be tested without the MolGpKa weights."""
    mol = Chem.MolFromSmiles(smiles)
    return R.Profile(mol, [], list(sites), [], fixed, fixed != 0, Crippen.MolLogP(mol) if logp is None else logp, R.METHOD_MOLGPKA)


def resp(level):
    return R.Response(level, "")


# ================================================================ rules, no weights needed
# ---------------------------------------------------------------- ESI+
def test_esi_plus_permanent_charge_is_strong():
    r = R.esi_positive(prof("C[N+](C)(C)CCO", fixed=1), 2.7)
    assert r.level == R.STRONG and "permanent" in r.reason


def test_esi_plus_basic_site_protonated_in_solution_is_strong_from_ph_plus_one():
    p = lambda pka: prof("CCN(CC)CC", [site("base", pka, 2, "N")])                   # noqa: E731
    assert R.esi_positive(p(4.0), 3.0).level == R.STRONG                              # pKa = pH + 1
    assert R.esi_positive(p(3.99), 3.0).level == R.MODERATE                           # just below
    assert "protonated in solution" in R.esi_positive(p(7.6), 2.7).reason


def test_esi_plus_basic_site_near_ph_is_moderate_and_far_below_ph_still_moderate():
    p = prof("CCN(CC)CC", [site("base", 7.6, 2, "N")])
    near = R.esi_positive(p, 7.0)
    assert near.level == R.MODERATE and "partly protonated" in near.reason
    far = R.esi_positive(p, 11.0)                                                     # pH 11: 7.6 < pH - 1, not protonated in solution
    assert far.level == R.MODERATE and "not protonated in solution" in far.reason


def test_esi_plus_no_basic_site_uses_heteroatoms_and_weakly_basic_groups():
    assert R.esi_positive(prof(GLUCOSE), 2.7).level == R.MODERATE                    # 6 O: Na+ / NH4+ adducts
    assert R.esi_positive(prof("CCOC(=O)C"), 2.7).level == R.MODERATE                # ester: carbonyl and ether
    assert R.esi_positive(prof("CC(=O)C"), 2.7).level == R.MODERATE                  # one carbonyl in a small molecule
    assert R.esi_positive(prof("CCCCCCCCCC"), 2.7).level == R.WEAK                   # C and H only
    assert R.esi_positive(prof("ClC(Cl)(Cl)Cl"), 2.7).level == R.WEAK               # halogens only
    assert R.esi_positive(prof("[O-][N+](=O)c1ccccc1"), 2.7).level == R.WEAK         # nitro is not a basic group
    assert R.esi_positive(prof("OC(=O)c1ccccc1"), 2.7).level == R.WEAK               # carboxyl atoms belong to the acid group
    r = R.esi_positive(prof("CCCO"), 2.7)
    assert r.level == R.WEAK and "no carbonyl, amide or ether" in r.reason           # a lone alcohol is not in the weakly basic list


def test_esi_plus_single_heteroatom_in_a_large_nonpolar_molecule_is_weak():
    r = R.esi_positive(prof("CCCCCCCCCCCCCCCCO"), 2.7)                               # hexadecanol: 17 heavy atoms, logP 5
    assert r.level == R.WEAK and "large non-polar" in r.reason
    assert R.esi_positive(prof("CCCCCCCCCCCCCCCCC=O"), 2.7).level == R.WEAK          # even a carbonyl
    assert R.esi_positive(prof("CCCCCCCCCCCCCCCCO", logp=1.0), 2.7).level == R.WEAK  # no weak group and not 3 atoms: still weak
    assert R.esi_positive(prof("CCCCCCCCCCCCC=O"), 2.7).level == R.MODERATE          # 14 heavy atoms: not "large"


def test_esi_plus_adduct_route_is_capped_by_an_acidic_site():
    aspirin = "CC(=O)Oc1ccccc1C(=O)O"
    assert R.esi_positive(prof(aspirin), 2.7).level == R.MODERATE                    # ester carbonyl and ether, nothing acidic known
    capped = R.esi_positive(prof(aspirin, [site("acid", 3.4)]), 2.7)
    assert capped.level == R.WEAK and "anion former" in capped.reason
    assert R.esi_positive(prof(aspirin, [site("acid", 12.0)]), 2.7).level == R.MODERATE   # an acid above pKa 10 does not cap


def test_esi_plus_acid_that_is_a_counterion_does_not_cap_and_a_basic_site_beats_the_cap():
    assert R.esi_positive(prof("CCN(CC)CC(=O)O", [site("acid", 3.4), site("base", 9.0, 3, "N")]), 2.7).level == R.STRONG


# ---------------------------------------------------------------- ESI-
def test_esi_minus_acid_pka_cutoffs():
    p = lambda pka: prof("CC(=O)O", [site("acid", pka)])                              # noqa: E731
    assert R.esi_negative(p(6.0), 2.7).level == R.STRONG                              # pKa <= 6, even at an acidic pH
    assert R.esi_negative(p(6.01), 2.7).level == R.MODERATE
    assert R.esi_negative(p(10.0), 2.7).level == R.MODERATE
    r = R.esi_negative(p(10.01), 2.7)
    assert r.level == R.WEAK and "not deprotonated" in r.reason
    assert R.esi_negative(prof("CCC"), 2.7).level == R.WEAK                           # no acidic site at all


def test_esi_minus_acid_deprotonated_in_solution_is_strong():
    phenol = prof("Oc1ccccc1", [site("acid", 10.0)])
    assert R.esi_negative(phenol, 7.0).level == R.MODERATE
    r = R.esi_negative(phenol, 12.0)
    assert r.level == R.STRONG and "deprotonated in solution" in r.reason


def test_esi_minus_polyol_is_moderate_unless_the_molecule_has_a_base():
    assert R.esi_negative(prof(GLUCOSE), 2.7).level == R.MODERATE
    assert R.esi_negative(prof("OCC(O)CO"), 2.7).level == R.MODERATE                  # glycerol: 3 OH
    assert R.esi_negative(prof("OCCO"), 2.7).level == R.WEAK                          # 2 OH
    amino_sugar = prof("N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O", [site("base", 8.3, 0, "N")])
    r = R.esi_negative(amino_sugar, 2.7)
    assert r.level == R.WEAK and "favours ESI+" in r.reason


def test_esi_minus_acid_next_to_a_permanent_cation_is_its_counterion():
    carnitine = prof("C[N+](C)(C)CC(O)CC(=O)O", [site("acid", 3.8)], fixed=1)
    r = R.esi_negative(carnitine, 2.7)
    assert r.level == R.WEAK and "counterion" in r.reason
    two = prof("C[N+](C)(C)CC(CC(=O)O)C(=O)O", [site("acid", 3.0), site("acid", 4.0, 1)], fixed=1)
    assert R.esi_negative(two, 2.7).level == R.STRONG                                  # the second acid is free
    assert R.esi_negative(prof("OS(=O)(=O)CC", fixed=-1), 2.7).level == R.STRONG       # permanent anion


def test_thresholds_are_one_dataclass_and_change_the_verdicts():
    p = prof("CC(=O)O", [site("acid", 5.0)])
    assert R.esi_negative(p, 2.7).level == R.STRONG
    assert R.esi_negative(p, 2.7, R.Thresholds(acid_strong_pka=4.0)).level == R.MODERATE
    b = prof("CCN(CC)CC", [site("base", 5.0, 2, "N")])
    assert R.esi_positive(b, 2.7).level == R.STRONG
    assert R.esi_positive(b, 2.7, R.Thresholds(base_protonated=3.0)).level == R.MODERATE


# ---------------------------------------------------------------- polarity
@pytest.mark.parametrize("pos,neg,expected", [
    (R.STRONG, R.STRONG, R.BOTH), (R.STRONG, R.MODERATE, R.BOTH), (R.MODERATE, R.STRONG, R.BOTH), (R.MODERATE, R.MODERATE, R.BOTH),
    (R.STRONG, R.WEAK, R.POS), (R.MODERATE, R.WEAK, R.POS), (R.WEAK, R.STRONG, R.NEG), (R.WEAK, R.MODERATE, R.NEG), (R.WEAK, R.WEAK, R.POOR),
])
def test_polarity_table(pos, neg, expected):
    assert R.polarity(resp(pos), resp(neg))[0] == expected


def test_polarity_both_names_the_mode_the_ph_favours():
    glycine = prof("NCC(=O)O", [site("acid", 2.1), site("base", 9.6, 0, "N")])
    assert R.solution_favour(glycine, 2.7) == R.POS                                    # base protonated, acid not ionized
    assert R.solution_favour(glycine, 7.0) is None                                     # both formed in solution
    assert R.solution_favour(glycine, 12.0) == R.NEG
    note = R.polarity(resp(R.STRONG), resp(R.STRONG), glycine, 2.7)[1]
    assert "ESI+ is usually the more sensitive" in note
    assert "try both modes" in R.polarity(resp(R.STRONG), resp(R.STRONG), glycine, 7.0)[1]


# ---------------------------------------------------------------- EI: M+• levels
@pytest.mark.parametrize("rel,level", [(100.0, R.VISIBLE), (5.0, R.VISIBLE), (4.99, R.WEAK_ION), (1.0, R.WEAK_ION), (0.99, R.ABSENT), (0.0, R.ABSENT)])
def test_mol_ion_level_thresholds(rel, level):
    assert R.ion_level(rel) == level


def _fake_neims(monkeypatch, peaks, mz=100):
    monkeypatch.setattr(ne, "predict", lambda smiles: {"smiles": smiles, "formula": "", "molecular_ion_mz": mz, "model_mass": mz, "peaks": peaks, "warnings": []})


@pytest.mark.parametrize("rel,level", [(5.0, "visible"), (4.9, "weak"), (1.0, "weak"), (0.9, "not expected")])
def test_molecular_ion_reads_the_neims_intensity_at_the_nominal_mass(monkeypatch, rel, level):
    _fake_neims(monkeypatch, [(57, 100.0), (99, 20.0), (100, rel)])
    ion = R.molecular_ion("CCCC", Chem.MolFromSmiles("CCCC"))
    assert ion.value == level and ion.rel == rel and ion.mz == 100 and "NEIMS" in ion.note


def test_molecular_ion_missing_from_the_predicted_peaks_is_not_expected(monkeypatch):
    _fake_neims(monkeypatch, [(57, 100.0), (99, 20.0)])                                # NEIMS drops peaks below 0.1 %
    ion = R.molecular_ion("CCCC", Chem.MolFromSmiles("CCCC"))
    assert ion.value == R.ABSENT and ion.rel == 0.0


def test_molecular_ion_without_neims_weights_says_so_and_gives_the_rule_of_thumb(monkeypatch):
    monkeypatch.setenv("MZDESK_NEIMS_WEIGHTS", "/nonexistent/neims_weights.npz")
    ion = R.molecular_ion("c1ccccc1", Chem.MolFromSmiles("c1ccccc1"))
    assert ion.value == "NEIMS not installed" and ion.rel is None
    assert ion.note.startswith("rule of thumb: M+• visible")
    row = R.ei_reason({"gc": "yes", "reasons": ["Rule: no active H, no carbonyl, MW 78"]}, ion)
    assert "NEIMS not installed; rule of thumb" in row


def test_molecular_ion_structure_outside_the_neims_range_falls_back_to_the_rule_of_thumb(monkeypatch):
    def too_heavy(smiles):
        raise ValueError("The molecular ion m/z 1200 is above the model's maximum (m/z 999).")
    monkeypatch.setattr(ne, "predict", too_heavy)
    ion = R.molecular_ion("C" * 80, Chem.MolFromSmiles("C" * 80))
    assert ion.value == "not available" and "above the model's maximum" in ion.note


@pytest.mark.parametrize("smiles,level", [
    ("c1ccccc1", R.VISIBLE),                       # aromatic
    ("CC=CC=CC=O", R.VISIBLE),                     # conjugated
    ("C1CCCCC1", R.VISIBLE),                       # ring
    ("CCCCCCCCCCCCCCCCCCO", R.ABSENT),             # long-chain alcohol
    ("CC(C)CC(C)(C)CC(C)(C)CC(C)(C)C", R.ABSENT),  # highly branched alkane
    ("CCCCCCCC", R.WEAK_ION),                      # acyclic alkane
])
def test_rule_of_thumb_for_the_molecular_ion(smiles, level):
    assert R.rule_of_thumb(Chem.MolFromSmiles(smiles))[0] == level


# ---------------------------------------------------------------- GC verdict = the Platform page's
@pytest.mark.parametrize("smiles", [LIDOCAINE, IBUPROFEN, GLUCOSE, "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "C[N+](C)(C)CCO", PC, "CC(N)Cc1ccccc1",
                                    "c1cc2ccc3cccc4ccc(c1)c2c34", "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
                                    "CCOP(=O)(OCC)OCC", "COP(C)(=O)OC", "OC(=O)CNCP(=O)(O)O", "C[N+](C)(C)CC(=O)[O-]", "[Na+].[O-]C(=O)c1ccccc1"])
def test_gc_verdict_is_the_one_the_platform_page_shows(smiles):
    platform = client.post("/api/modules/platform_fit/run", json={"structure": {"mode": "SMILES / InChI", "value": smiles}}).json()
    gc, _mol = R.gc_verdict(smiles)
    tile = next(t for t in platform["outputs"]["platform"]["tiles"] if t["label"] == "GC-MS")        # yes / partial / no tile
    assert {"yes": "yes", "after derivatization": "partial", "no": "no"}[gc["gc"]] == tile["status"]


@pytest.mark.parametrize("smiles,expected,why", [
    (LIDOCAINE, "yes", "small amide-amine, one active H: routine GC drug"),
    (IBUPROFEN, "yes", "small carboxylic acid, MW 206: runs directly (often derivatized for peak shape)"),
    (GLUCOSE, "after derivatization", "five OH: MeOx-TMS is the standard GC route"),
    ("CC(N)Cc1ccccc1", "after derivatization", "primary amine: amphetamine is acylated or silylated for GC"),
    ("C[N+](C)(C)CCO", "no", "quaternary ammonium: ionic, not volatile"),
    (PC, "no", "phospholipid: ionic, MW 734"),
    ("c1cc2ccc3cccc4ccc(c1)c2c34", "yes", "PAH: volatile and stable"),
    ("COP(=O)(OC)OC=C(Cl)Cl", "yes", "dichlorvos: a phosphate triester has no P-OH and runs in GC"),
    ("CCOP(=O)(OCC)OCC", "yes", "triethyl phosphate"),
    ("COP(C)(=O)OC", "yes", "DMMP: a phosphonate diester"),
    ("OC(=O)CNCP(=O)(O)O", "no", "glyphosate: free P-OH, ionic"),
])
def test_gc_amenability_expectations(smiles, expected, why):
    assert R.gc_verdict(smiles)[0]["gc"] == expected, why


def test_gc_verdict_uses_the_largest_component_of_a_salt():
    free, _ = R.gc_verdict(LIDOCAINE)
    salt, _ = R.gc_verdict(LIDOCAINE + ".Cl")
    assert salt["gc"] == free["gc"]


# ---------------------------------------------------------------- functional-group fallback (MolGpKa weights missing)
def _kinds(smiles):
    return sorted((s.kind, s.pka) for s in R.group_sites(Chem.MolFromSmiles(smiles)))


def test_group_sites_cover_the_main_groups():
    assert _kinds(LIDOCAINE) == [("base", 10.0)]                                       # the anilide N is not a site
    assert _kinds(IBUPROFEN) == [("acid", 4.5)]
    assert _kinds("NCC(=O)O") == [("acid", 4.5), ("base", 10.0)]
    assert _kinds("Nc1ccccc1") == [("base", 4.6)]
    assert _kinds("c1ccncc1") == [("base", 5.0)]
    assert _kinds("c1ncccn1") == [("base", 1.5), ("base", 1.5)]                        # pyrimidine: two weakly basic ring N
    assert _kinds("Oc1ccccc1") == [("acid", 10.0)]
    assert _kinds("C[N+](C)(C)[O-]") == [("base", 4.7)] and _kinds("[O-][n+]1ccccc1") == [("base", 0.8)]      # N-oxides (Refinement E)
    assert _kinds("[O-][N+](=O)c1ccccc1") == []                                         # nitro is not an N-oxide
    assert _kinds("Cc1ccc(cc1)S(O)(=O)=O") == [("acid", 0.0)]
    assert _kinds("NC(N)=N") == [("base", 12.0)]
    assert _kinds("CCCCCC") == []


@needs_pka
def test_n_oxide_gets_a_textbook_base_site_next_to_the_predicted_ones():
    a = R.assess("C[N+](C)(C)[O-]", 2.7)
    assert [(s.kind, s.pka) for s in a.profile.sites] == [("base", 4.7)] and len(a.profile.textbook) == 1
    assert a.pos.level == R.STRONG and a.polarity == R.POS
    assert any("N-oxide site O5" in h and "not predicted" in h for h in a.hits)
    assert R.assess("[O-][N+](=O)c1ccccc1", 2.7).profile.textbook == []


def test_assess_without_molgpka_weights_falls_back_to_functional_groups(monkeypatch):
    monkeypatch.setenv("MZDESK_MOLGPKA_WEIGHTS", "/nonexistent/molgpka_weights.npz")
    a = R.assess(LIDOCAINE, 2.7)
    assert a.profile.method == R.METHOD_GROUPS and a.polarity == R.POS
    assert "MolGpKa weights not available" in a.hits[0] and "textbook values" in a.hits[0]
    b = R.assess(IBUPROFEN, 2.7)
    assert b.polarity == R.NEG and b.pos.level == R.WEAK and b.neg.level == R.STRONG
    c = R.assess("C[N+](C)(C)CCO", 2.7)                                                # permanent charge needs no pKa
    assert c.polarity == R.POS and c.pos.level == R.STRONG


def test_page_runs_without_molgpka_and_without_neims(monkeypatch):
    monkeypatch.setenv("MZDESK_MOLGPKA_WEIGHTS", "/nonexistent/molgpka_weights.npz")
    monkeypatch.setenv("MZDESK_NEIMS_WEIGHTS", "/nonexistent/neims_weights.npz")
    r = _run(LIDOCAINE)
    assert r.status_code == 200, r.text
    out = r.json()["outputs"]
    chips = {c["label"]: c for c in out["chips"]["chips"]}
    assert chips["Polarity"]["value"] == "+"
    assert chips["EI"]["value"] == "usable"
    ion = chips["M+•"]
    assert ion["level"] == "good" and ion["value"] == "visible (rule of thumb)" and _help(ion)[0].startswith("rule of thumb")
    assert any("MolGpKa weights not available" in h for h in out["rule_hits"]["items"])
    assert any("NEIMS weights not installed" in h for h in out["rule_hits"]["items"])


# ================================================================ page
def test_spec_page_is_live_with_the_requested_inputs_and_outputs():
    m = next(m for m in MODULES if m["id"] == "ionization")
    assert m["live"] and m["subgroup"] == "detectability"
    assert m["one_liner"] == "Compound + pH → ESI polarity and EI usability"
    ids = [f["id"] for f in m["inputs"]]
    assert ids == ["structure", "ph"] and "platform" not in ids
    ph = m["inputs"][1]
    assert (ph["type"], ph["min"], ph["max"], ph["value"]) == ("slider", 1, 13, 2.7)
    blocks = {b["id"]: b for b in m["outputs"]}
    assert [b["id"] for b in m["outputs"]] == ["chips", "sites_structure", "rule_hits"]               # polarity is the first chip, not a card
    assert blocks["chips"]["type"] == "chips" and blocks["chips"]["title"] == "Ionization" and blocks["chips"]["help"]
    assert blocks["sites_structure"]["type"] == "figure" and blocks["sites_structure"]["help"]
    assert blocks["rule_hits"]["type"] == "notes" and blocks["rule_hits"]["title"] == "Rule hits" and blocks["rule_hits"]["collapsed"] is True
    assert m["example"]["ph"] == 2.7


@needs_pka
@needs_neims
def test_page_outputs_for_lidocaine():
    r = _run(LIDOCAINE)
    assert r.status_code == 200, r.text
    body = r.json()
    out = body["outputs"]
    assert set(out) == {"chips", "sites_structure", "rule_hits"}
    chips = out["chips"]["chips"]
    assert [c["label"] for c in chips] == ["Polarity", "ESI+", f"ESI{MINUS}", "EI", "M+•"]
    pol, esi_pos, esi_neg, ei, ion = chips
    assert (esi_pos["level"], esi_pos["value"]) == ("good", "strong") and (esi_neg["level"], esi_neg["value"]) == ("poor", "weak")
    assert "N3" in esi_pos["help"] and "pKa" in esi_pos["help"]                          # the deciding reason is the chip's help
    assert (ei["level"], ei["value"]) == ("good", "usable") and "small molecule" in _help(ei)[0]
    assert ion["level"] in ("good", "possible", "poor") and ion["value"] in ("visible", "weak", "not expected") and "NEIMS" in _help(ion)[0]
    assert (pol["level"], pol["value"]) == ("good", "+") and _help(pol) == ["ESI+ responds better than ESI−"]
    assert out["sites_structure"]["svg"].startswith("<?xml") or "<svg" in out["sites_structure"]["svg"]
    items = out["rule_hits"]["items"]
    assert any("MolGpKa" in t for t in items) and any("N3 base" in t for t in items) and any("amide N" in t for t in items)
    assert any("basic pKa ≥" in t and "acidic pKa ≤" in t for t in items)                  # the pKa thresholds of the rule; no general disclaimers
    assert body["warnings"] == []


@needs_pka
@needs_neims
def test_derivatization_note_and_not_gc_amenable_note():
    chips = _chips(GLUCOSE)
    assert (chips["EI"]["level"], chips["EI"]["value"]) == ("possible", "after derivatization")
    assert chips["M+•"]["value"].endswith("(underivatized)") and any("underivatized molecule" in h for h in _help(chips["M+•"]))
    hits = _run(GLUCOSE).json()["outputs"]["rule_hits"]["items"]
    assert any("NEIMS predicts the underivatized molecule" in t for t in hits)
    choline = _chips("C[N+](C)(C)CCO")
    assert (choline["EI"]["level"], choline["EI"]["value"]) == ("poor", "not usable")
    assert choline["M+•"]["value"].endswith("(hypothetical)") and any("not GC-amenable" in h for h in _help(choline["M+•"]))


@needs_pka
@needs_neims
def test_ph_changes_the_esi_plus_verdict_of_a_base():
    low, high = _chips(LIDOCAINE, 2.7)["ESI+"], _chips(LIDOCAINE, 11.0)["ESI+"]
    assert (low["level"], low["value"]) == ("good", "strong")
    assert (high["level"], high["value"]) == ("possible", "moderate") and "not protonated in solution" in high["help"]


@needs_pka
@needs_neims
def test_zwitterion_note_names_the_favoured_mode_and_salt_is_reduced_to_its_largest_part():
    out = _chips("NCC(=O)O")["Polarity"]
    assert (out["level"], out["value"]) == ("good", "both") and any("ESI+ is usually the more sensitive" in h for h in _help(out))
    r = _run(LIDOCAINE + ".Cl")
    assert r.status_code == 200 and any("disconnected parts" in w for w in r.json()["warnings"])
    assert r.json()["outputs"]["chips"]["chips"][0]["value"] == "+"


@needs_pka
@needs_neims
def test_polarity_chip_levels_and_help():
    chips = _chips(GLUCOSE)
    assert chips["Polarity"]["value"] == "both" and chips["Polarity"]["level"] == "good"
    assert _help(chips["Polarity"]) == ["both modes respond at least moderately", "try both modes at pH 2.7"]            # no favoured mode: try both
    assert _chips(IBUPROFEN)["Polarity"]["value"] == R.NEG and _chips(IBUPROFEN)["Polarity"]["level"] == "good"
    poor = _chips("CCCCCC")["Polarity"]                                                                                  # hydrocarbon: weak in both modes
    assert (poor["level"], poor["value"]) == ("poor", "poor") and _help(poor) == ["both modes respond weakly"]


def _hex(rgb):
    return "".join(f"{int(c * 255):02X}" for c in rgb)   # RDKit truncates


from mzdesk.shared import colors as _theme  # noqa: E402
ACID_HEX, BASE_HEX = _hex(_theme.HL_ACID), _hex(_theme.HL_BASE)   # highlight colours of the figure (acid coral, base teal), as RDKit writes them


def test_figure_highlights_acid_and_base_sites_and_a_permanent_charge():
    from mzdesk.platform_fit import ionization_viz as viz
    glycine = viz.sites_svg(prof("NCC(=O)O", [site("acid", 2.1, 3), site("base", 9.6, 0, "N")]))
    assert ACID_HEX in glycine and BASE_HEX in glycine
    assert ACID_HEX in viz.sites_svg(prof("CC(=O)O", [site("acid", 4.5, 3)])) and BASE_HEX not in viz.sites_svg(prof("CC(=O)O", [site("acid", 4.5, 3)]))
    assert BASE_HEX in viz.sites_svg(prof("C[N+](C)(C)CCO", fixed=1))                  # permanent charge marked like a base
    glucose = viz.sites_svg(prof(GLUCOSE, [site("acid", 12.0, 1)]))                    # acid above pKa 10: not highlighted
    assert ACID_HEX not in glucose and "<svg" in glucose
    assert "<svg" in viz.sites_svg(prof("CCCCCC"))                                     # nothing to highlight: plain depiction


@needs_pka
def test_figure_of_a_quaternary_ammonium_marks_the_permanent_charge():
    assert BASE_HEX in _run("C[N+](C)(C)CCO").json()["outputs"]["sites_structure"]["svg"]


def test_input_errors_are_400():
    assert _run("C1CC").status_code == 400                                              # unparsable SMILES
    assert _run("").status_code == 400
    assert _run(LIDOCAINE, 15).status_code == 400
    assert _run(LIDOCAINE, -1).status_code == 400
    assert client.post("/api/modules/ionization/run", json={"structure": {"mode": "SMILES / InChI", "value": LIDOCAINE}, "ph": "abc"}).status_code == 400


def test_ph_defaults_to_2_7_when_missing():
    r = client.post("/api/modules/ionization/run", json={"structure": {"mode": "SMILES / InChI", "value": IBUPROFEN}})
    assert r.status_code == 200
    assert any("pH 2.7" in t for t in r.json()["outputs"]["rule_hits"]["items"])


# ================================================================ validation set (needs MolGpKa)
# (name, SMILES, allowed polarity, allowed ESI+ levels, allowed ESI- levels, rationale). Expectations are chemistry, not rule output.
P, N, B, O = R.POS, R.NEG, R.BOTH, R.POOR
S, M, W = R.STRONG, R.MODERATE, R.WEAK
LC_HITS = [
    # bases
    ("lidocaine", LIDOCAINE, P, (S,), (W,), "tertiary amine, pKa about 7.9: protonated at pH 2.7; the anilide NH is not acidic"),
    ("amphetamine", "CC(N)Cc1ccccc1", P, (S,), (W,), "primary aliphatic amine, pKa about 9.9"),
    ("propranolol", "CC(C)NCC(O)COc1cccc2ccccc12", P, (S,), (W,), "secondary amine, pKa about 9.5; the aryl ether has no acidic group"),
    ("imipramine", "CN(C)CCCN1c2ccccc2CCc2ccccc12", P, (S,), (W,), "tricyclic with a tertiary alkylamine, pKa about 9.5"),
    ("nicotine", "CN1CCCC1c1cccnc1", P, (S,), (W,), "pyrrolidine N pKa about 8 (pyridine N about 3)"),
    ("diphenhydramine", "CN(C)CCOC(c1ccccc1)c1ccccc1", P, (S,), (W,), "tertiary amine, pKa about 9"),
    ("quinine", "COc1ccc2nccc([C@H](O)[C@H]3C[C@@H]4CCN3C[C@@H]4C=C)c2c1", P, (S,), (W,), "quinuclidine N pKa about 8.5 and quinoline N about 4.2"),
    ("trimethylamine N-oxide", "C[N+](C)(C)[O-]", P, (S,), (W,), "amine oxide, the O protonates (pKa about 4.7): TMAO is a positive-mode analyte (MolGpKa has no N-oxide pattern: textbook pKa)"),
    ("cocaine", "COC(=O)[C@H]1[C@@H]2CC[C@H](C[C@@H]1OC(=O)c1ccccc1)N2C", P, (S,), (W,), "tropane tertiary amine, pKa 8.6"),
    ("codeine", "COc1ccc2C[C@H]3N(C)CC[C@@]45[C@@H](Oc1c24)[C@@H](O)C=C[C@@H]35", P, (S,), (W,), "tertiary amine, pKa 8.2; opioids are positive-mode analytes"),
    ("dextromethorphan", "COc1ccc2C[C@H]3N(C)CC[C@@]4(CCCC[C@@H]34)c2c1", P, (S,), (W,), "tertiary amine, pKa 9.2"),
    ("diazepam", "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc12", P, (M, S), (W,), "benzodiazepine imine N, pKa about 3: [M+H]+ is routine"),
    ("atrazine", "CCNc1nc(Cl)nc(NC(C)C)n1", P, (M, S), (W,), "triazine N, pKa 1.7, ethylamino / isopropylamino: [M+H]+ is routine"),
    ("melamine", "Nc1nc(N)nc(N)n1", P, (M, S), (W,), "triazine N, pKa 5.0: [M+H]+ is the routine mode"),
    ("creatinine", "CN1CC(=O)NC1=N", P, (S,), (W,), "amidine-type N, pKa 4.8: [M+H]+"),
    ("carbamazepine", "NC(=O)N1c2ccccc2C=Cc2ccccc12", P, (M, S), (W,), "urea-type carbonyl protonates: [M+H]+ is routine; no acidic site"),
    ("aflatoxin B1", "COC1=C2C3=C(C(=O)OC3)C(=O)OC2=C2[C@@H]3C=CO[C@@H]3OC2=C1", P, (M, S), (W,), "lactones, ethers: [M+H]+ / Na+ adducts in positive mode, no acid"),
    ("tripalmitin", "CCCCCCCCCCCCCCCC(=O)OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCC", P, (M, S), (W,), "triacylglycerol: NH4+ / Na+ adducts in positive mode"),
    ("ceramide d18:1/16:0", "CCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)/C=C/CCCCCCCCCCCCC", P, (M, S), (W, M), "amide and two OH: [M+H]+ / [M+H−H2O]+ in positive mode; acetate adduct only in negative"),
    ("ascorbic acid", "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O", N, (W,), (S,), "enol acid pKa 4.2: [M−H]−"),
    ("cholic acid", "C[C@H](CCC(O)=O)[C@H]1CC[C@H]2[C@@H]3[C@H](O)C[C@@H]4C[C@H](O)CC[C@]4(C)[C@H]3C[C@H](O)[C@]12C", N, (W,), (S,), "bile acid: carboxylic acid, [M−H]−"),
    ("taurocholic acid", "C[C@H](CCC(=O)NCCS(O)(=O)=O)[C@H]1CC[C@H]2[C@@H]3[C@H](O)C[C@@H]4C[C@H](O)CC[C@]4(C)[C@H]3C[C@H](O)[C@]12C", N, (W,), (S,), "taurine conjugate: sulfonic acid"),
    ("pentachlorophenol", "Oc1c(Cl)c(Cl)c(Cl)c(Cl)c1Cl", N, (W,), (S,), "phenol with five Cl, pKa 4.7: [M−H]− is the standard"),
    ("PFOA", "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", N, (W,), (S,), "perfluorinated acid, pKa below 1: [M−H]−"),
    ("naproxen", "COc1ccc2cc(ccc2c1)[C@H](C)C(O)=O", (N, B), (W, M), (S,), "carboxylic acid pKa 4.2: [M−H]− is the standard (positive mode also works)"),
    ("indomethacin", "COc1ccc2n(C(=O)c3ccc(Cl)cc3)c(C)c(CC(O)=O)c2c1", (N, B), (W, M), (S,), "carboxylic acid pKa 4.5: [M−H]− is the standard"),
    ("warfarin", "CC(=O)CC(c1ccccc1)C1=C(O)c2ccccc2OC1=O", (N, B), (W, M), (S,), "4-hydroxycoumarin, pKa 5: [M−H]− is the standard"),
    ("glyphosate", "OC(=O)CNCP(O)(=O)O", (N, B), (S,), (S,), "phosphonic acid, carboxylic acid and a secondary amine: mostly negative mode"),
    ("tyrosine", "N[C@@H](Cc1ccc(O)cc1)C(O)=O", B, (S,), (S, M), "amino acid; phenol is irrelevant"),
    ("phosphatidylethanolamine 34:1", "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OCCN)OC(=O)CCCCCCC/C=C\\CCCCCCCC", B, (S,), (S,), "phosphate and primary amine: lipidomics runs PE in both modes"),
    ("sucrose", "OC[C@H]1O[C@@](CO)(O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@@H]1O", B, (M,), (M,), "disaccharide, no ionizable site: Na+ / NH4+ adducts (+), Cl− / formate adducts (−)"),
    ("vitamin D3", "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]\\2[C@@]1(CCC/C2=C\\C=C3/C[C@H](CCC3=C)O)C", O, (W,), (W,), "secosteroid with one OH: poor in ESI, APCI is used"),
    ("benzo[a]pyrene", "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34", O, (W,), (W,), "polycyclic aromatic hydrocarbon, no heteroatom"),
    ("DDT", "ClC(Cl)(Cl)C(c1ccc(Cl)cc1)c1ccc(Cl)cc1", O, (W,), (W,), "organochlorine, no polar group: GC-MS or APCI"),
    ("adenosine", "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O", P, (M, S), (W,), "adenine N1 pKa 3.5: [M+H]+ is the routine mode for nucleosides; ribose OH are not acidic enough for ESI−"),
    ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", P, (M, S), (W,), "purine N: pKa below 1 but high proton affinity, [M+H]+ is routine; no acidic H"),
    # acids
    ("ibuprofen", IBUPROFEN, N, (W,), (S,), "carboxylic acid pKa 4.4: [M−H]− is the standard mode; nothing to protonate"),
    ("benzoic acid", "OC(=O)c1ccccc1", N, (W,), (S,), "carboxylic acid pKa 4.2, no basic site"),
    ("p-toluenesulfonic acid", "Cc1ccc(cc1)S(O)(=O)=O", N, (W,), (S,), "sulfonic acid, pKa below 0"),
    ("salicylic acid", "OC(=O)c1ccccc1O", N, (W,), (S,), "carboxylic acid pKa 3.0 (phenol pKa 13 is irrelevant)"),
    ("palmitic acid", "CCCCCCCCCCCCCCCC(O)=O", N, (W,), (S,), "fatty acid: [M−H]− standard, positive mode poor"),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O", N, (W,), (S,), "carboxylic acid pKa 3.5; the ester carbonyl alone does not make a positive-mode analyte"),
    ("phenobarbital", "CCC1(c2ccccc2)C(=O)NC(=O)NC1=O", N, (W,), (M, S), "barbiturate imide NH pKa 7.4: barbiturates are run in negative mode"),
    ("bisphenol A", "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1", N, (W,), (M, S), "two phenols, pKa about 10: negative mode is the standard"),
    ("citric acid", "OC(=O)CC(O)(CC(O)=O)C(O)=O", N, (W,), (S,), "tricarboxylic acid"),
    ("glucose-6-phosphate", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1OP(O)(O)=O", N, (W,), (S,), "phosphate monoester, pKa about 1.5: phosphorylated sugars are measured in negative mode"),
    ("sodium dodecyl sulfate", "CCCCCCCCCCCCOS(O)(=O)=O", N, (W,), (S,), "alkyl sulfate, pKa below 0"),
    ("uracil", "O=c1cc[nH]c(=O)[nH]1", N, (W,), (M, S), "imide NH pKa 9.5: [M−H]− is the usual mode, [M+H]+ is poor"),
    # amphoteric
    ("phenylalanine", "N[C@@H](Cc1ccccc1)C(=O)O", B, (S,), (S,), "ammonium (pKa about 9) and carboxylic acid (pKa about 2.2): both modes work"),
    ("levofloxacin", "C[C@H]1COc2c(N3CCN(C)CC3)c(F)cc3c(=O)c(C(=O)O)cn1c23", B, (S,), (S,), "quinolone zwitterion: amine and carboxylic acid, both modes used"),
    ("cetirizine", "OC(=O)COCCN1CCN(CC1)C(c1ccccc1)c1ccc(Cl)cc1", B, (S,), (S,), "zwitterionic antihistamine (piperazine N and carboxylic acid), both modes used"),
    ("sulfamethoxazole", "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1", (P, B), (M, S), (W, M), "aniline N for ESI+ (standard), sulfonamide NH pKa about 6 for ESI−"),
    ("acetaminophen", "CC(=O)Nc1ccc(O)cc1", (P, B), (M, S), (W, M), "amide and phenol, no basic site: [M+H]+ via the amide / Na+ adducts, phenol pKa 9.5 for ESI−"),
    # neutral polar
    ("glucose", GLUCOSE, B, (M,), (M,), "no ionizable site; sugars appear only as Na+ / NH4+ adducts (+) and Cl− / formate adducts (−)"),
    # quaternary and zwitterionic quaternary
    ("choline", "C[N+](C)(C)CCO", P, (S,), (W,), "quaternary ammonium: permanent cation, no acid"),
    ("carnitine", "C[N+](C)(C)CC(O)CC([O-])=O", P, (S,), (W,), "quaternary N+ with a carboxylate counterion: no [M−H]−"),
    ("tetramethylammonium", "C[N+](C)(C)C", P, (S,), (W,), "quaternary ammonium"),
    ("phosphatidylcholine 32:0", PC, P, (S,), (W, M), "permanent cation; the phosphate is the counterion, ESI− only as anion adducts"),
    # hydrocarbons / nonpolar
    ("pyrene", "c1cc2ccc3cccc4ccc(c1)c2c34", O, (W,), (W,), "PAH, no heteroatom: nothing to protonate or deprotonate (APPI / APCI territory)"),
    ("naphthalene", "c1ccc2ccccc2c1", O, (W,), (W,), "aromatic hydrocarbon"),
    ("hexadecane", "CCCCCCCCCCCCCCCC", O, (W,), (W,), "saturated hydrocarbon"),
    ("nitrobenzene", "[O-][N+](=O)c1ccccc1", O, (W,), (W,), "nitroarenes are APCI-negative analytes; ESI is poor in both modes"),
    ("phenol", "Oc1ccccc1", O, (W,), (W,), "small phenol, pKa 10 and one O: poor in both ESI modes"),
    # steroids
    ("cholesterol", "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C", O, (W,), (W,), "sterol with one OH: ESI is poor, APCI is used"),
    ("testosterone", "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C", P, (M, S), (W,), "enone protonates: [M+H]+ is the routine mode; no acid"),
    ("cortisol", "C[C@]12CCC(=O)C=C1CC[C@@H]1[C@@H]2[C@@H](O)C[C@@]2(C)[C@H]1CC[C@]2(O)C(=O)CO", (P, B), (M, S), (W, M), "polyoxygenated steroid: [M+H]+ / Na+ in positive mode, formate adduct in negative"),
]
# Genuine misses: (name, SMILES, what the chemistry says, field, why the rules disagree). Strict xfail: a fix shows up as an unexpected pass.
LC_MISSES = [
    ("triethylamine", "CCN(CC)CC", P, (S,), (W,), "aliphatic tertiary amine, pKa 10.7: ESI+ strong",
     "MolGpKa predicts pKa 1.7 for Me3N / Et3N / Et2NH (real 9.8 / 10.7 / 11): the base looks unprotonated, so ESI+ is only moderate (polarity still +)"),
    ("4-nitrophenol", "Oc1ccc(cc1)[N+](=O)[O-]", N, (W,), (S,), "pKa 7.2: [M−H]− is a strong signal even from slightly acidic eluents",
     "the spec ladder puts pKa 6-10 at moderate; nitrophenols are among the acids that behave as strong (polarity - is right)"),
    ("morphine", "CN1CC[C@]23[C@@H]4Oc5c3c(C[C@@H]1[C@@H]2C=C[C@@H]4O)ccc5O", P, (S,), (W,), "tertiary amine pKa 8.2: positive mode only, ESI− is not used",
     "the phenol (pKa 9.4) is rated moderate for ESI− and strong + moderate gives both"),
    ("dopamine", "NCCc1ccc(O)c(O)c1", P, (S,), (W,), "primary amine pKa 10: positive mode only",
     "the catechol (pKa 9.7) is rated moderate for ESI−, so the polarity comes out both"),
    ("triclosan", "Oc1cc(Cl)ccc1Oc1ccc(Cl)cc1Cl", N, (W,), (S,), "chlorophenol, pKa 7.9: [M−H]− is a strong signal",
     "same ladder issue as 4-nitrophenol: acids with pKa 6-10 are rated moderate, but phenols at pKa 7-8 respond as strong ones (polarity - is right)"),
    ("chloramphenicol", "OC[C@@H](NC(=O)C(Cl)Cl)[C@H](O)c1ccc(cc1)[N+]([O-])=O", N, (W,), (M, S), "[M−H]− (and Cl− adducts) is the standard mode, positive mode is poor",
     "no site with pKa <= 10 is found (dichloroacetamide NH, nitroarene, two aliphatic OH), so ESI− is weak while ESI+ is moderate (amide / carbonyl count): polarity +"),
    ("penicillin G", "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(O)=O", (P, B), (M, S), (S,), "beta-lactams are routinely run as [M+H]+ as well as [M−H]−",
     "Refinement A: the acid (pKa 2.6) caps the weakly basic amide / lactam route at weak, so the polarity is -"),
    ("hydrochlorothiazide", "NS(=O)(=O)c1cc2c(NCNS2(=O)=O)cc1Cl", N, (W,), (M, S), "[M−H]− is the standard mode (sulfonamide NH pKa 7.9 and 9.2), positive mode is poor",
     "MolGpKa gives the ring aniline-type N a base pKa of 2.1 (partly protonated at pH 2.7), so ESI+ is moderate and the polarity both"),
    ("metformin", "CN(C)C(=N)NC(N)=N", P, (S,), (W,), "biguanide, pKa 12.4: positive mode only",
     "MolGpKa predicts an acid site pKa 9.9 on a guanidine NH (real pKa above 13), so ESI− is moderate and the polarity both"),
    ("estradiol", "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCc1cc(O)ccc13", N, (W,), (M,), "phenol: negative mode is the usual one for underivatized estrogens, sensitivity modest",
     "MolGpKa phenol pKa 10.1 sits just above the pKa 10 cut-off, so ESI− is weak and the polarity poor (borderline)"),
]


def _check_lc(smiles, polarity, pos, neg, ph=2.7):
    a = A(smiles, ph)
    polarity = (polarity,) if isinstance(polarity, str) else polarity
    assert a.polarity in polarity, f"polarity {a.polarity!r}, expected {polarity}: {a.polarity_note}"
    assert a.pos.level in pos, f"ESI+ {a.pos.level!r}, expected {pos}: {a.pos.reason}"
    assert a.neg.level in neg, f"ESI− {a.neg.level!r}, expected {neg}: {a.neg.reason}"


@needs_pka
@pytest.mark.parametrize("name,smiles,polarity,pos,neg,why", [pytest.param(*c, id=c[0]) for c in LC_HITS])
def test_lc_esi_validation(name, smiles, polarity, pos, neg, why):
    _check_lc(smiles, polarity, pos, neg)


@needs_pka
@pytest.mark.parametrize("name,smiles,polarity,pos,neg,truth,reason", [
    pytest.param(*c, id=c[0], marks=pytest.mark.xfail(strict=True, reason=c[6])) for c in LC_MISSES])
def test_lc_esi_known_misses(name, smiles, polarity, pos, neg, truth, reason):
    _check_lc(smiles, polarity, pos, neg)


def test_validation_set_has_at_least_20_compounds_with_a_rationale():
    assert len(LC_HITS) + len(LC_MISSES) >= 20
    assert all(len(c[5]) > 10 for c in LC_HITS) and all(len(c[5]) > 10 and len(c[6]) > 10 for c in LC_MISSES)
    assert len({c[1] for c in LC_HITS + LC_MISSES}) == len(LC_HITS) + len(LC_MISSES)


# ---------------------------------------------------------------- EI validation (needs NEIMS)
EI_HITS = [
    ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", ("visible",), "purine ring system: M+• is the base peak"),
    ("pyrene", "c1cc2ccc3cccc4ccc(c1)c2c34", ("visible",), "four fused aromatic rings: M+• is the base peak"),
    ("naphthalene", "c1ccc2ccccc2c1", ("visible",), "aromatic: M+• is the base peak"),
    ("benzoic acid", "OC(=O)c1ccccc1", ("visible",), "aromatic acid: strong M+•"),
    ("imipramine", "CN(C)CCCN1c2ccccc2CCc2ccccc12", ("visible",), "tricyclic: M+• at m/z 280 clearly present"),
    ("1-hexadecanol", "CCCCCCCCCCCCCCCCO", ("weak", "not expected"), "long-chain alcohol: loses H2O, M+• usually absent"),
    ("1-octadecanol", "CCCCCCCCCCCCCCCCCCO", ("weak", "not expected"), "long-chain alcohol: M+• usually absent"),
    ("heptamethylnonane", "CC(C)CC(C)(C)CC(C)(C)CC(C)(C)C", ("weak", "not expected"), "highly branched alkane: cleaves at the branch points"),
    ("tert-butanol", "CC(C)(C)O", ("weak", "not expected"), "tertiary alcohol: M+• at m/z 74 is almost absent"),
    ("cocaine", "COC(=O)[C@H]1[C@@H]2CC[C@H](C[C@@H]1OC(=O)c1ccccc1)N2C", ("visible",), "tropane ester: M+• at m/z 303 clearly present"),
    ("amphetamine", "CC(N)Cc1ccccc1", ("weak", "not expected"), "alpha-cleavage gives m/z 44 almost exclusively: M+• at about 1 %"),
    ("methamphetamine", "CNC(C)Cc1ccccc1", ("weak", "not expected"), "alpha-cleavage gives m/z 58: M+• at about 1 %"),
]
EI_MISSES = [
    ("glucose", GLUCOSE, ("not expected",), "sugars fragment completely under EI: no M+•", "NEIMS predicts 7 % at m/z 180 (a cyclic polyol outside what the training spectra of sugars show)"),
    ("pristane", "CC(C)CCCC(C)CCCC(C)CCCC(C)C", ("weak", "not expected"), "branched isoprenoid alkane: M+• below 1 %", "NEIMS predicts about 17 % at m/z 268"),
    ("neopentane", "CC(C)(C)C", ("weak", "not expected"), "quaternary carbon: M+• at m/z 72 is practically absent", "NEIMS predicts 42 % at m/z 72"),
    ("isooctane", "CC(C)CC(C)(C)C", ("weak", "not expected"), "2,2,4-trimethylpentane: M+• at m/z 114 below 1 %", "NEIMS predicts 11 % at m/z 114"),
]


@needs_neims
@pytest.mark.parametrize("name,smiles,levels,why", [pytest.param(*c, id=c[0]) for c in EI_HITS])
def test_ei_molecular_ion_validation(name, smiles, levels, why):
    ion = R.molecular_ion(smiles, Chem.MolFromSmiles(smiles))
    assert ion.value in levels, f"{ion.value} ({ion.rel} %), expected {levels}: {why}"


@needs_neims
@pytest.mark.parametrize("name,smiles,levels,truth,reason", [pytest.param(*c, id=c[0], marks=pytest.mark.xfail(strict=True, reason=c[4])) for c in EI_MISSES])
def test_ei_molecular_ion_known_misses(name, smiles, levels, truth, reason):
    ion = R.molecular_ion(smiles, Chem.MolFromSmiles(smiles))
    assert ion.value in levels, f"{ion.value} ({ion.rel} %), expected {levels}: {truth}"


@needs_neims
def test_nominal_mass_is_the_peak_that_is_read():
    res = ne.predict("Cn1cnc2c1c(=O)n(C)c(=O)n2C")
    ion = R.molecular_ion("Cn1cnc2c1c(=O)n(C)c(=O)n2C", Chem.MolFromSmiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C"))
    assert ion.mz == res["molecular_ion_mz"] == 194 and ion.rel == pytest.approx(dict(res["peaks"])[194])
