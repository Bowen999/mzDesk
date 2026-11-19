"""Separation 页的测试：位点匹配与特征、numpy 前向对照原代码的参考值、Henderson–Hasselbalch 与 logD、分离模式规则、页面输出、错误输入。

参考值 tests/fixtures/molgpka_reference.json 来自原代码 + 原权重（PyTorch + torch_geometric，validation/molgpka/）。
权重没安装（python -m mzdesk.platform_fit.pka_setup）时，需要权重的测试跳过；其余测试（规则、数学、中和、位点匹配）不需要权重。
"""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
from fastapi.testclient import TestClient
from rdkit import Chem

from mzdesk.api.app import app
from mzdesk.platform_fit import molgpka_smarts, pka_engine as E, pka_setup
from mzdesk.platform_fit import separation_rules as R

client = TestClient(app)
ROOT = Path(__file__).parent.parent
REF = json.loads((Path(__file__).parent / "fixtures" / "molgpka_reference.json").read_text())["molecules"]
needs_weights = pytest.mark.skipif(not Path(E.weights_path()).is_file(), reason="MolGpKa weights not installed (python -m mzdesk.platform_fit.pka_setup)")
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
GLYCINE = "NCC(=O)O"
PC = "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCCCC"
GLUCOSE = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"


def _run(smiles, ph=2.7):
    return client.post("/api/modules/separation/run", json={"structure": {"mode": "SMILES / InChI", "value": smiles}, "ph": ph})


def _modes(out):
    """{chip label: chip} of the Separation mode chips."""
    return {c["label"]: c for c in out["modes"]["chips"]}


def _levels(out):
    return {label: c["level"] for label, c in _modes(out).items()}


def _site(kind, pka, atom=0, element="O"):
    return E.Site(atom=atom, element=element, kind=kind, pka=pka)


# ---------------------------------------------------------------- the SMARTS table and site matching (no weights)
def test_smarts_table_is_the_original_one():
    assert len(molgpka_smarts.PATTERNS) == 143
    assert sum(p[3] == "A" for p in molgpka_smarts.PATTERNS) == 54 and sum(p[3] == "B" for p in molgpka_smarts.PATTERNS) == 89
    assert all(Chem.MolFromSmarts(p[1]) is not None for p in molgpka_smarts.PATTERNS)


def test_acid_site_is_the_lost_hydrogen_and_base_site_the_heavy_atom():
    mol_h = Chem.AddHs(Chem.MolFromSmiles(IBUPROFEN))
    acid, base = E.find_site_atoms(mol_h)
    assert base == [] and len(acid) == 1
    h = mol_h.GetAtomWithIdx(acid[0])
    assert h.GetSymbol() == "H" and h.GetNeighbors()[0].GetSymbol() == "O"
    mol_h = Chem.AddHs(Chem.MolFromSmiles("CCN(CC)CC"))
    assert E.find_site_atoms(mol_h) == ([], [2])                        # triethylamine N


def test_molecule_without_ionizable_site_has_no_sites():
    assert E.find_site_atoms(Chem.AddHs(Chem.MolFromSmiles("CCCCCC"))) == ([], [])


def test_static_features_follow_the_original_layout():
    mol_h = Chem.AddHs(Chem.MolFromSmiles("CC(=O)O"))
    x = E.static_features(mol_h)
    assert x.shape == (mol_h.GetNumAtoms(), 27) and x.dtype == np.float32
    assert x[:, 0].sum() == 2 and x[:, 1].sum() == 4 and x[:, 2].sum() == 2     # C, H, O one-hot
    assert not x[:, 21:23].any()                                                # donor / acceptor columns: always 0, as in the original code
    assert x[0, 10] == 4 and x[1, 10] == 3                                      # degree of CH3 and C(=O)
    a_hat = E.normalized_adjacency(mol_h)
    assert np.allclose(a_hat, a_hat.T) and np.all(np.diag(a_hat) > 0)


# ---------------------------------------------------------------- neutralization and fixed charges (no weights)
@pytest.mark.parametrize("smiles,neutral,fixed,permanent", [
    ("[NH3+]CC([O-])=O", "NCC(=O)O", 0, False),                                  # glycine zwitterion
    ("CC(=O)[O-].[Na+]", "CC(=O)O", 0, False),                                   # salt: the largest fragment
    ("C[N+](C)(C)CC(O)CC([O-])=O", "C[N+](C)(C)CC(O)CC(=O)O", 1, True),          # carnitine: the carboxylate is protonated, N+ stays
    ("C[N+](=O)[O-]", "C[N+](=O)[O-]", 0, False),                                # nitro: charge-separated but neutral
    ("C[N+](C)(C)[O-]", "C[N+](C)(C)[O-]", 0, False),                            # N-oxide
    ("C[n+]1ccccc1", "C[n+]1ccccc1", 1, True),
    ("[O-]c1ccc([N+](=O)[O-])cc1", "O=[N+]([O-])c1ccc(O)cc1", 0, False),          # phenolate -> phenol, nitro stays
])
def test_prepare_neutralizes_what_can_be_neutralized(smiles, neutral, fixed, permanent):
    p = E.prepare(smiles)
    assert Chem.MolToSmiles(p.mol) == Chem.MolToSmiles(Chem.MolFromSmiles(neutral))
    assert p.fixed_charge == fixed and p.permanent is permanent


def test_prepare_phosphatidylcholine_keeps_only_the_quaternary_nitrogen_charged():
    p = E.prepare(PC)
    assert p.permanent and p.fixed_charge == 1
    assert any("Permanent charge" in n for n in p.notes) and any("neutralized" in n for n in p.notes)
    assert "P(=O)(O)" in Chem.MolToSmiles(p.mol)


def test_prepare_keeps_the_atom_order_of_a_single_fragment_input():
    p = E.prepare("CC(=O)Nc1ccc(O)cc1")
    assert [a.GetSymbol() for a in p.mol.GetAtoms()] == [a.GetSymbol() for a in Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1").GetAtoms()]


def test_prepare_warns_about_salts_unusual_elements_and_size():
    assert any("2 disconnected parts" in n for n in E.prepare("CC(=O)O.CCN").notes)
    assert any("Element(s) Si" in n for n in E.prepare("C[Si](C)(C)CC(=O)O").notes)
    assert any("heavy atoms" in n for n in E.prepare("C" * 120).notes)
    with pytest.raises(ValueError, match="limited to"):
        E.prepare("C" * (E.MAX_HEAVY_ATOMS + 1))
    with pytest.raises(ValueError, match="could not parse"):
        E.prepare("C1CC")


# ---------------------------------------------------------------- chemically non-basic N (no weights: synthetic sites on real molecules)
LIDOCAINE = "CCN(CC)CC(=O)Nc1c(C)cccc1C"                 # N3 amine (index 2), N9 anilide N (index 8)
CAFFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"                 # index 1, 8, 12 are pyrrole-type n; index 3 is the pyridine-type n
SULFAMETHOXAZOLE = "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1"    # index 4 sulfonamide N, index 12 aniline N


def _filter(smiles, *sites):
    return E.drop_nonbasic_sites(Chem.MolFromSmiles(smiles), list(sites))


def test_amide_n_is_not_a_base_but_stays_an_acid():
    amide_base, amine_base, amide_acid = _site("base", 5.2, 8, "N"), _site("base", 7.6, 2, "N"), _site("acid", 13.9, 8, "N")
    kept, dropped = _filter(LIDOCAINE, amide_base, amine_base, amide_acid)
    assert kept == [amine_base, amide_acid] and dropped == [(amide_base, "amide N")]


@pytest.mark.parametrize("smiles,atom,reason", [
    ("CC(=O)NC", 3, "amide N"), ("CC(=O)NC(C)=O", 3, "amide N"), ("COC(=O)NC", 4, "amide N"), ("CNC(N)=O", 1, "amide N"),   # amide, imide, carbamate, urea
    ("CNC(N)=O", 3, "amide N"), ("CC(=S)NC", 3, "amide N"), ("O=C1CCCN1", 5, "amide N"),                               # urea NH2, thioamide, lactam
    ("CS(=O)(=O)NC", 4, "sulfonamide N"), ("CS(=O)(=O)N(C)C", 4, "sulfonamide N"),
    ("c1ccc2[nH]ccc2c1", 4, "pyrrole-type aromatic N"), ("Cn1ccnc1", 1, "pyrrole-type aromatic N"), ("c1c[nH]cn1", 2, "pyrrole-type aromatic N"),
])
def test_nonbasic_n_rules(smiles, atom, reason):
    site = _site("base", 5.0, atom, "N")
    assert Chem.MolFromSmiles(smiles).GetAtomWithIdx(atom).GetSymbol() == "N"
    kept, dropped = _filter(smiles, site)
    assert kept == [] and dropped == [(site, reason)]


@pytest.mark.parametrize("smiles,atom", [
    ("c1ccncc1", 3), ("CCN(CC)CC", 2), ("Nc1ccccc1", 0), ("CN(C)C(=N)N", 1), ("CC(=N)N", 3), ("CC(=N)N", 2),     # pyridine, amine, aniline, guanidine / amidine N
    ("c1c[nH]cn1", 4),                                                                           # the pyridine-type n of imidazole
])
def test_basic_n_is_kept(smiles, atom):
    site = _site("base", 6.0, atom, "N")
    assert Chem.MolFromSmiles(smiles).GetAtomWithIdx(atom).GetSymbol() in ("N", "n")
    assert _filter(smiles, site) == ([site], [])


def test_only_base_sites_on_n_are_filtered():
    acid_o = _site("acid", 4.0, 2, "O")
    kept, dropped = _filter("CC(=O)O", acid_o)
    assert kept == [acid_o] and dropped == []
    sulfonamide_acid = _site("acid", 7.0, 4, "N")
    assert _filter(SULFAMETHOXAZOLE, sulfonamide_acid) == ([sulfonamide_acid], [])         # sulfonamide NH acid stays


def test_caffeine_and_sulfamethoxazole_filter():
    ring = [_site("base", 2.0 + i / 10, a, "N") for i, a in enumerate((1, 3, 8, 12))]
    kept, dropped = _filter(CAFFEINE, *ring)
    assert [s.atom for s in kept] == [3] and sorted(s.atom for s, _ in dropped) == [1, 8, 12]
    sulfonamide_base, aniline_base = _site("base", 4.3, 4, "N"), _site("base", 4.1, 12, "N")
    kept, dropped = _filter(SULFAMETHOXAZOLE, sulfonamide_base, aniline_base)
    assert kept == [aniline_base] and dropped == [(sulfonamide_base, "sulfonamide N")]


# ---------------------------------------------------------------- Henderson-Hasselbalch, independent sites (no weights)
def test_one_acid_is_half_ionized_at_its_pka():
    acid = _site("acid", 4.4)
    assert E.charged_fraction(acid, 4.4) == pytest.approx(0.5)
    assert E.charged_fraction(acid, 6.4) == pytest.approx(100 / 101)               # 2 units above: 99 % ionized
    base = _site("base", 9.0)
    assert E.charged_fraction(base, 9.0) == pytest.approx(0.5)
    assert E.charged_fraction(base, 7.0) == pytest.approx(100 / 101)


def test_charge_distribution_sums_to_one_and_matches_the_formula():
    sites = [_site("acid", 4.0), _site("acid", 8.0, atom=1), _site("base", 9.0, atom=2, element="N")]
    ph = np.array([0.0, 2.0, 6.0, 9.0, 14.0])
    charges, dist, uncharged = E.charge_states(sites, ph)
    assert list(charges) == [-2, -1, 0, 1]
    assert np.allclose(dist.sum(axis=1), 1.0)
    pa, pb, pc = (E.charged_fraction(s, ph) for s in sites)
    assert np.allclose(uncharged, (1 - pa) * (1 - pb) * (1 - pc))
    assert np.allclose(dist[:, 3], (1 - pa) * (1 - pb) * pc)                        # only the base charged -> +1
    assert np.allclose(dist @ charges, pc - pa - pb)                                # expected charge = sum of the site charges
    assert dist[0].argmax() == 3 and dist[-1].argmax() == 0                         # pH 0: +1; pH 14: -2 and the base neutral -> index 0 is -2


def test_zwitterion_counts_as_net_zero_but_not_as_uncharged():
    sites = [_site("acid", 2.3), _site("base", 9.6, atom=1, element="N")]            # glycine
    charges, dist, uncharged = E.charge_states(sites, [6.0])
    assert list(charges) == [-1, 0, 1]
    assert dist[0, 1] > 0.999                                                        # net 0
    pa, pb = (E.charged_fraction(s, 6.0) for s in sites)
    assert uncharged[0] == pytest.approx((1 - pa) * (1 - pb), rel=1e-6) and uncharged[0] < 1e-6      # the uncharged microstate is very rare


def test_fixed_charge_shifts_the_distribution():
    charges, dist, _ = E.charge_states([_site("acid", 2.0)], [7.0], fixed_charge=1)
    assert list(charges) == [0, 1] and dist[0, 0] > 0.99                             # N+ with an ionized acid: net 0
    charges, dist, _ = E.charge_states([], [7.0], fixed_charge=1)
    assert list(charges) == [1] and dist[0, 0] == 1.0


def test_no_sites_means_neutral_at_every_ph():
    charges, dist, uncharged = E.charge_states([], E.PH_GRID)
    assert list(charges) == [0] and np.all(dist == 1) and np.all(uncharged == 1)


# ---------------------------------------------------------------- logD from the formula (sites mocked: no weights)
@pytest.fixture
def ibuprofen_sites(monkeypatch):
    monkeypatch.setattr(E, "predict_sites", lambda mol: [_site("acid", 4.4, atom=14)])


def test_logd_is_logp_plus_log_fraction_neutral(ibuprofen_sites):
    a = E.analyze(IBUPROFEN)
    st = a.at(2.7)
    f = 1 / (1 + 10 ** (2.7 - 4.4))
    assert st["f_neutral"] == pytest.approx(f) and st["logd"] == pytest.approx(a.logp + np.log10(f))
    assert st["charge"] == 0 and st["p_charge"] == pytest.approx(f)
    i = int(np.argmin(abs(a.ph - 7.4)))
    assert a.logd[i] == pytest.approx(a.logp + np.log10(1 / (1 + 10 ** (7.4 - 4.4))))
    assert len(a.ph) == 141 and a.ph[0] == 0 and a.ph[-1] == 14 and a.ph[1] == 0.1
    assert a.logd[0] == pytest.approx(a.logp, abs=1e-3) and a.logd[-1] < a.logp - 9


def test_permanent_charge_has_no_logd(monkeypatch):
    monkeypatch.setattr(E, "predict_sites", lambda mol: [_site("acid", 1.2, atom=24)])
    a = E.analyze(PC)
    assert a.permanent and np.all(np.isnan(a.logd)) and np.all(a.f_neutral == 0)
    st = a.at(2.7)
    assert st["logd"] is None and st["charge"] == 0                                  # N+ and the ionized phosphate
    assert a.at(0.0)["charge"] == 1


# ---------------------------------------------------------------- rules (no weights)
def _rules(**kw):
    base = dict(ph=2.7, logd=0.5, tpsa=50.0, mw=300.0, charge=0, f_neutral=1.0, permanent=False, fixed_charge=0, n_stereo=0)
    return {r[0]: r for r in R.mode_rows(**{**base, **kw})}


def test_every_mode_is_listed_in_order():
    assert [r[0] for r in R.mode_rows(ph=2.7, logd=1.0, tpsa=50, mw=300, charge=0, f_neutral=1.0, permanent=False, fixed_charge=0, n_stereo=0)] == \
        ["RP", "HILIC", "Ion-pair RP", "Mixed-mode / ion exchange", "Chiral", "SEC"]


@pytest.mark.parametrize("logd,expected", [(3.0, "good"), (1.0, "good"), (0.99, "possible"), (-1.0, "possible"), (-1.01, "poor"), (-6.0, "poor")])
def test_rp_thresholds(logd, expected):
    assert _rules(logd=logd)["RP"][1] == expected


@pytest.mark.parametrize("logd,tpsa,expected", [(-0.01, 20, "good"), (0.0, 20, "possible"), (0.5, 99.9, "possible"), (0.5, 100.0, "good"),
                                                (1.0, 150, "possible"), (2.0, 20, "possible"), (2.01, 20, "poor"), (2.5, 150, "poor")])
def test_hilic_thresholds(logd, tpsa, expected):
    assert _rules(logd=logd, tpsa=tpsa)["HILIC"][1] == expected


def test_ion_pair_needs_charge_and_low_logd():
    r = _rules(logd=-2.0, charge=1, f_neutral=0.01)["Ion-pair RP"]
    assert r[1] == "possible" and "not MS-friendly" in r[2]
    assert _rules(logd=-2.0, charge=0, f_neutral=0.9)["Ion-pair RP"][1] == "poor"          # neutral: nothing to pair
    assert _rules(logd=0.3, charge=-1, f_neutral=0.01)["Ion-pair RP"][1] == "poor"         # charged but plain RP retains it
    assert _rules(logd=-2.0, charge=0, f_neutral=0.01)["Ion-pair RP"][1] == "possible"     # zwitterion: net 0 but charged groups


def test_mixed_mode_for_permanent_charge_or_large_net_charge():
    assert _rules(charge=2, f_neutral=0.0)["Mixed-mode / ion exchange"][1] == "possible"
    assert _rules(charge=-2, f_neutral=0.0)["Mixed-mode / ion exchange"][2].startswith("net charge −2")
    assert _rules(charge=1, f_neutral=0.01)["Mixed-mode / ion exchange"][1] == "poor"
    r = _rules(logd=None, charge=1, f_neutral=0.0, permanent=True, fixed_charge=1)["Mixed-mode / ion exchange"]
    assert r[1] == "possible" and "permanent charge (+1)" in r[2]


def test_permanent_charge_rules():
    rules = _rules(logd=None, charge=0, f_neutral=0.0, permanent=True, fixed_charge=1)
    assert rules["RP"][1] == "possible" and "logD is not defined" in rules["RP"][2]
    assert rules["HILIC"][1] == "good" and rules["Ion-pair RP"][1] == "possible"


def test_chiral_counts_stereocentres_and_knows_meso():
    assert _rules(n_stereo=0)["Chiral"][1:] == ["not needed", "no stereocentre"]
    r = _rules(n_stereo=2, n_unassigned=2)["Chiral"]
    assert r[1] == "needed only for enantiomer separation" and "2 stereocentres (2 not specified" in r[2]
    assert _rules(n_stereo=2, meso=True)["Chiral"][1] == "not needed"


@pytest.mark.parametrize("mw,expected", [(300, "poor"), (5000, "poor"), (5000.5, "good"), (66000, "good")])
def test_sec_only_above_5000(mw, expected):
    assert _rules(mw=mw)["SEC"][1] == expected


def test_mode_chips_map_recommendations_to_levels_and_short_labels():
    rows = R.mode_rows(ph=2.7, logd=1.5, tpsa=50, mw=300, charge=0, f_neutral=1.0, permanent=False, fixed_charge=0, n_stereo=2, n_unassigned=1)
    chips = R.mode_chips(rows)
    assert [(c[0], c[1]) for c in chips] == [("RP", "good"), ("HILIC", "possible"), ("Ion-pair RP", "poor"), ("Mixed-mode / IEX", "poor"),
                                             ("Chiral", "possible"), ("SEC", "poor")]
    assert chips[0][2] == rows[0][2]                                              # the help is the reason
    assert chips[4][2] == "Needed only for enantiomer separation: " + rows[4][2]
    flat = R.mode_chips(R.mode_rows(ph=2.7, logd=1.5, tpsa=50, mw=300, charge=0, f_neutral=1.0, permanent=False, fixed_charge=0, n_stereo=0))
    assert flat[4] == ("Chiral", "poor", "no stereocentre")


def test_all_thresholds_live_in_one_dataclass():
    th = R.Thresholds(rp_good=2.0)
    assert _rules(logd=1.5)["RP"][1] == "good" and R.mode_rows(ph=2.7, logd=1.5, tpsa=50, mw=300, charge=0, f_neutral=1.0, permanent=False,
                                                                 fixed_charge=0, n_stereo=0, th=th)[0][1] == "possible"


def test_ph_windows_are_pka_plus_minus_two():
    sites = [_site("acid", 4.4), _site("base", 7.6, atom=1, element="N")]
    assert R.site_window_text(sites[0]) == "O1 (acid, pKa 4.4): >99 % neutral at pH ≤ 2.4 / anion at pH ≥ 6.4"
    assert R.site_window_text(sites[1]) == "N2 (base, pKa 7.6): >99 % cation at pH ≤ 5.6 / neutral at pH ≥ 9.6"
    assert R.robust_ph_ranges(sites, E.PH_GRID) == [(0.0, 2.4), (9.6, 14.0)]
    assert R.robust_ph_ranges([_site("acid", 7.0)], E.PH_GRID) == [(0.0, 5.0), (9.0, 14.0)]
    assert R.robust_ph_ranges([_site("acid", 4.0), _site("base", 6.0, atom=1, element="N")], E.PH_GRID) == [(0.0, 2.0), (8.0, 14.0)]
    assert R.robust_ph_ranges([_site("acid", 5.0), _site("base", 7.0, atom=1, element="N"), _site("acid", 9.0, atom=2)], E.PH_GRID) == [(0.0, 3.0), (11.0, 14.0)]


def test_ph_notes_say_when_the_chosen_ph_is_on_a_slope():
    sites = [_site("acid", 4.4)]
    f = E.charge_states(sites, E.PH_GRID)[2]
    notes = R.ph_notes(sites, 3.5, E.PH_GRID, f, 3.0, False)
    assert any("pH 3.5 is within ± 2 of the pKa of O1 (4.4)" in n for n in notes)
    assert any(n.startswith("Robust mobile-phase pH") and "pH 0–2.4 or 6.4–14" in n for n in notes)
    assert not any("is within" in n for n in R.ph_notes(sites, 2.0, E.PH_GRID, f, 3.0, False))
    assert "No ionizable site" in R.ph_notes([], 2.7, E.PH_GRID, np.ones(141), 1.0, False)[0]


def test_ph_notes_list_only_the_closest_sites_but_use_all_for_the_ranges():
    sites = [_site("acid", 1.0 + i, atom=i) for i in range(14)]
    f = E.charge_states(sites, E.PH_GRID)[2]
    notes = R.ph_notes(sites, 7.0, E.PH_GRID, f, 1.0, False)
    assert sum("(acid, pKa" in n for n in notes) == 10 and any(n.startswith("4 more sites") for n in notes)
    assert any(n.startswith("No pH between 0 and 14 keeps every site") for n in notes)


# ---------------------------------------------------------------- the network against the original code
@needs_weights
@pytest.mark.parametrize("name", list(REF))
def test_numpy_network_matches_original_code(name):
    mol_h = Chem.AddHs(Chem.MolFromSmiles(REF[name]["smiles"]))
    acid, base = E.find_site_atoms(mol_h)
    assert acid == sorted(int(a) for a in REF[name]["same_graph"]["acid"]) and base == sorted(int(a) for a in REF[name]["same_graph"]["base"])
    for kind, aids in (("acid", acid), ("base", base)):
        got = E.predict_aids(mol_h, kind, aids)
        for aid in aids:
            assert got[aid] == pytest.approx(REF[name]["same_graph"][kind][str(aid)], abs=1e-3)


@needs_weights
def test_sites_are_one_prediction_per_heavy_atom_sorted_by_pka():
    sites = E.predict_sites(Chem.MolFromSmiles("Nc1ccccc1"))                        # NH2: two equivalent H atoms, one acid site
    assert [(s.kind, s.element) for s in sites].count(("acid", "N")) == 1
    pkas = [s.pka for s in sites]
    assert pkas == sorted(pkas)


@needs_weights
def test_site_sets_after_the_nonbasic_n_filter():
    def sites(smiles):
        a = E.analyze(smiles)
        return sorted((x.label, x.kind) for x in a.sites), sorted((x.label, x.kind, why) for x, why in a.ignored)
    kept, ignored = sites(LIDOCAINE)                                       # only the amine is basic; the anilide NH stays an (very weak) acid
    assert kept == [("N3", "base"), ("N9", "acid")] and ignored == [("N9", "base", "amide N")]
    kept, ignored = sites(CAFFEINE)                                        # only the pyridine-type N (index 3)
    assert kept == [("N4", "base")] and [i[0] for i in ignored] == ["N13", "N2", "N9"]
    assert {i[2] for i in ignored} == {"pyrrole-type aromatic N"}
    kept, ignored = sites(SULFAMETHOXAZOLE)                                # sulfonamide N: no base, still an acid; the aniline base is a model error that stays
    assert ("N5", "acid") in kept and ("N5", "base") not in kept and ("N13", "base") in kept
    assert ignored == [("N5", "base", "sulfonamide N")]
    a = E.analyze(SULFAMETHOXAZOLE)
    assert next(s.pka for s in a.sites if (s.label, s.kind) == ("N5", "acid")) == pytest.approx(7.02, abs=0.05)
    assert next(s.pka for s in a.sites if (s.label, s.kind) == ("N13", "base")) == pytest.approx(4.13, abs=0.05)


@needs_weights
def test_raw_predict_sites_still_returns_every_pattern_hit():
    raw = E.predict_sites(Chem.MolFromSmiles(LIDOCAINE))
    assert ("N9", "base") in [(s.label, s.kind) for s in raw] and len(raw) == 3


# ---------------------------------------------------------------- literature pKa (known deviations are in validation/molgpka, not asserted)
def _literature_module():
    path = ROOT / "validation" / "molgpka" / "literature.py"
    if not path.is_file():
        pytest.skip("validation/molgpka/literature.py not present")
    spec = importlib.util.spec_from_file_location("molgpka_literature", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@needs_weights
def test_literature_pka_for_the_groups_the_model_gets_right():
    lit = _literature_module()
    within_one = {"ibuprofen", "aspirin", "lidocaine", "morphine (amine)", "morphine (phenol)", "glycine (acid)", "glycine (amine)",
                  "salicylic acid", "procainamide"}
    seen = set()
    for name, smiles, smarts, kind, lo, hi, _ in lit.CASES:
        site = lit.predict(smiles, smarts, kind)
        assert site is not None, name                                                # the group was found as an ionizable site of the right type
        if name in within_one:
            assert abs(lit.error(site.pka, lo, hi)) <= 1.0, (name, site.pka)
            seen.add(name)
    assert seen == within_one


# ---------------------------------------------------------------- the page
@needs_weights
def test_ibuprofen_at_low_and_neutral_ph():
    body = _run(IBUPROFEN, 2.7).json()
    out = body["outputs"]
    assert body["status"] == "ok" and body["warnings"] == []
    v = out["ph_summary"]["values"]
    assert v["Dominant charge state"]["value"] == "0" and v["Fraction neutral"]["value"].endswith("%")
    assert 2.8 < float(v["logD at pH"]["value"]) < 3.3
    assert [r[0] for r in out["sites"]["rows"]] == ["O15"] and out["sites"]["rows"][0][2] == "acid" and out["sites"]["count"] == 1
    assert 3.9 < float(out["sites"]["rows"][0][1]) < 5.4
    assert _levels(out) == {"RP": "good", "HILIC": "poor", "Ion-pair RP": "poor", "Mixed-mode / IEX": "poor", "Chiral": "possible", "SEC": "poor"}
    assert list(_modes(out)) == ["RP", "HILIC", "Ion-pair RP", "Mixed-mode / IEX", "Chiral", "SEC"]            # the order of the rules
    assert "well retained" in _modes(out)["RP"]["help"] and _modes(out)["Chiral"]["help"].startswith("Needed only for enantiomer separation: 1 stereocentre")
    assert out["sites_structure"]["svg"].startswith("<?xml") or "<svg" in out["sites_structure"]["svg"]
    assert any(n.startswith("Robust mobile-phase pH") for n in out["rule_hits"]["items"])

    v = _run(IBUPROFEN, 7.4).json()["outputs"]
    assert v["ph_summary"]["values"]["Dominant charge state"]["value"] == "−1"
    assert float(v["ph_summary"]["values"]["logD at pH"]["value"]) < 1.0
    assert _levels(v)["RP"] == "possible"


@needs_weights
def test_lidocaine_page_ignores_the_anilide_n_and_is_singly_charged():
    body = _run(LIDOCAINE, 2.7).json()
    out = body["outputs"]
    assert [(r[0], r[2]) for r in out["sites"]["rows"]] == [("N3", "base"), ("N9", "acid")]
    assert out["ph_summary"]["values"]["Dominant charge state"]["value"] == "+1"
    assert -2.6 < float(out["ph_summary"]["values"]["logD at pH"]["value"].replace("\u2212", "-")) < -2.0
    assert out["ph_summary"]["values"]["logD at pH"]["note"].startswith("pH 2.7, logP ")
    assert out["rule_hits"]["items"][0] == "N9 amide N: not basic, base site ignored (MolGpKa predicted pKa 5.2)"
    assert out["sites_structure"]["svg"]
    assert R.ZWITTERION_NOTE not in out["rule_hits"]["items"]


@needs_weights
def test_curves_figure_has_both_curves_and_the_marker():
    fig = _run(IBUPROFEN, 4.0).json()["outputs"]["ph_curves"]["plotly"]
    names = [t.get("name") for t in fig["data"]]
    assert names[:2] == ["logD", "Net charge"]
    logd, charge = fig["data"][0], fig["data"][1]
    assert len(logd["x"]) == len(logd["y"]) == len(charge["y"]) == 141 and logd["x"][0] == 0 and logd["x"][-1] == 14
    assert charge["y"][0] == pytest.approx(0, abs=1e-3) and charge["y"][-1] == pytest.approx(-1, abs=1e-3)
    assert [t["x"] for t in fig["data"][2:]] == [[4.0], [4.0]]                      # one marker per curve at the chosen pH
    assert fig["layout"]["shapes"][0]["x0"] == 4.0


@needs_weights
def test_zwitterion_is_net_zero_but_charged():
    out = _run(GLYCINE, 6.0).json()["outputs"]
    v = out["ph_summary"]["values"]
    assert v["Dominant charge state"]["value"] == "0" and "zwitterion" in v["Dominant charge state"]["note"]
    assert v["Fraction neutral"]["value"] == "< 0.1 %"
    assert _levels(out)["Ion-pair RP"] == "possible"
    assert R.ZWITTERION_NOTE in out["rule_hits"]["items"] and R.ZWITTERION_NOTE.startswith("zwitterion: logD likely strongly underestimated")
    assert R.ZWITTERION_NOTE not in _run(GLYCINE, 2.0).json()["outputs"]["rule_hits"]["items"]          # net +1 there, not a zwitterion


@needs_weights
def test_permanently_charged_lipid_has_no_logd():
    body = _run(PC, 2.7).json()
    v = body["outputs"]["ph_summary"]["values"]
    assert v["logD at pH"]["value"] == "—" and "permanent charge" in v["logD at pH"]["note"] and v["Fraction neutral"]["value"] == "0 %"
    modes = _levels(body["outputs"])
    assert modes["RP"] == "possible" and modes["HILIC"] == "good" and modes["Mixed-mode / IEX"] == "possible"
    assert any("Permanent charge" in w for w in body["warnings"])
    fig = body["outputs"]["ph_curves"]["plotly"]
    assert [t.get("name") for t in fig["data"]][:1] == ["Net charge"]               # logD curve omitted


@needs_weights
def test_molecule_without_sites_has_flat_curves_and_no_structure_block():
    body = _run("CCCCCC", 2.7).json()
    out = body["outputs"]
    assert out["sites"]["rows"] == [] and "sites_structure" not in out
    v = out["ph_summary"]["values"]
    assert v["Dominant charge state"]["value"] == "0" and v["Fraction neutral"]["value"] == "100.0 %"
    assert any("No ionizable site" in n for n in out["rule_hits"]["items"])
    assert _levels(out)["Chiral"] == "poor" and _modes(out)["Chiral"]["help"] == "no stereocentre"
    assert out["ph_summary"]["values"]["logD at pH"]["value"] == f"{E.crippen_logp(Chem.MolFromSmiles('CCCCCC')):.1f}"       # logD = logP when nothing ionizes


@needs_weights
def test_alcohols_count_as_very_weak_acids():
    rows = _run(GLUCOSE, 2.7).json()["outputs"]["sites"]["rows"]                       # MolGpKa has a sugar-OH pattern (pKa about 12)
    assert rows and all(r[2] == "acid" and 11 < float(r[1]) < 14.5 for r in rows)


@needs_weights
def test_salt_input_is_reduced_to_the_largest_neutral_fragment():
    body = _run("CC(=O)[O-].[Na+]").json()
    assert body["status"] == "ok" and any("disconnected parts" in w for w in body["warnings"])
    assert [r[2] for r in body["outputs"]["sites"]["rows"]] == ["acid"]


def test_page_lists_the_ignored_sites_with_mocked_predictions(monkeypatch):
    monkeypatch.setattr(E, "predict_sites", lambda mol: [_site("base", 5.2, 8, "N"), _site("base", 7.6, 2, "N"), _site("acid", 13.9, 8, "N")])
    out = _run(LIDOCAINE, 2.7).json()["outputs"]
    assert out["rule_hits"]["items"][0] == "N9 amide N: not basic, base site ignored (MolGpKa predicted pKa 5.2)"
    assert [r[0] + r[2] for r in out["sites"]["rows"]] == ["N3base", "N9acid"]
    assert out["ph_summary"]["values"]["Dominant charge state"]["value"] == "+1"
    assert not any(n.startswith("Assumptions") for n in out["rule_hits"]["items"])             # the general assumptions live in the Docs page


def test_zwitterion_note_with_mocked_predictions(monkeypatch):
    monkeypatch.setattr(E, "predict_sites", lambda mol: [_site("acid", 2.3, 4, "O"), _site("base", 9.6, 0, "N")])
    items = _run(GLYCINE, 6.0).json()["outputs"]["rule_hits"]["items"]
    assert R.ZWITTERION_NOTE in items
    assert R.ZWITTERION_NOTE not in _run(GLYCINE, 1.0).json()["outputs"]["rule_hits"]["items"]


def test_page_with_mocked_sites_needs_no_weights(monkeypatch):
    monkeypatch.setattr(E, "predict_sites", lambda mol: [_site("acid", 4.4, atom=14)])
    body = _run(IBUPROFEN, 2.7).json()
    assert body["status"] == "ok"
    assert set(body["outputs"]) == {"ph_summary", "sites", "ph_curves", "modes", "rule_hits", "sites_structure"}
    assert body["outputs"]["sites"] == {"rows": [["O15", "4.4", "acid"]], "count": 1}
    assert body["outputs"]["ph_summary"]["values"]["Dominant charge state"]["value"] == "0"


def test_spec_block_ids_match_the_runner():
    from mzdesk.spec import MODULES
    mod = next(m for m in MODULES if m["id"] == "separation")
    assert mod["live"] and mod["subgroup"] == "detectability"
    assert [f["id"] for f in mod["inputs"]] == ["structure", "ph"] and mod["inputs"][1]["value"] == 2.7
    assert [(b["type"], b["id"]) for b in mod["outputs"]] == [("cards", "ph_summary"), ("figure", "ph_curves"), ("figure", "sites_structure"),
                                                              ("chips", "modes"), ("table", "sites"), ("notes", "rule_hits")]
    blocks = {b["id"]: b for b in mod["outputs"]}
    assert blocks["modes"]["title"] == "Separation mode" and blocks["sites"]["columns"] == ["Site", "pKa", "Type"]
    assert blocks["ph_summary"]["style"] == "metric" and blocks["rule_hits"]["collapsed"]
    items = blocks["ph_summary"]["items"]
    assert [i["label"] for i in items] == ["logD at pH", "Dominant charge state", "Fraction neutral"] and all(i["help"] for i in items)
    assert "logP" in items[0]["help"] and "neutral fraction" in items[0]["help"]                                     # what logD is lives in the card's ?


# ---------------------------------------------------------------- errors
def test_missing_weights_is_a_502_that_names_the_setup_command(monkeypatch, tmp_path):
    monkeypatch.setenv("MZDESK_MOLGPKA_WEIGHTS", str(tmp_path / "nope.npz"))
    r = _run(IBUPROFEN)
    assert r.status_code == 502 and pka_setup.SETUP_COMMAND in r.json()["detail"]


def test_corrupt_weights_is_a_502(monkeypatch, tmp_path):
    bad = tmp_path / "bad.npz"
    np.savez(bad, nothing=np.zeros(3))
    monkeypatch.setenv("MZDESK_MOLGPKA_WEIGHTS", str(bad))
    r = _run(IBUPROFEN)
    assert r.status_code == 502 and "unreadable" in r.json()["detail"]


@pytest.mark.parametrize("inputs", [
    {"structure": {"mode": "SMILES / InChI", "value": "C1CC"}, "ph": 2.7},
    {"structure": {"mode": "SMILES / InChI", "value": ""}, "ph": 2.7},
    {"structure": {"mode": "SMILES / InChI", "value": IBUPROFEN}, "ph": 15},
    {"structure": {"mode": "SMILES / InChI", "value": IBUPROFEN}, "ph": "acidic"},
    {"structure": {"mode": "SMILES / InChI", "value": "C" * 400}, "ph": 2.7},
])
def test_bad_input_is_400(inputs, monkeypatch):
    monkeypatch.setattr(E, "predict_sites", lambda mol: [])
    r = client.post("/api/modules/separation/run", json=inputs)
    assert r.status_code == 400 and r.json()["detail"]


def test_blank_ph_uses_the_default(monkeypatch):
    monkeypatch.setattr(E, "predict_sites", lambda mol: [_site("acid", 4.4, atom=14)])
    r = client.post("/api/modules/separation/run", json={"structure": {"mode": "SMILES / InChI", "value": IBUPROFEN}, "ph": ""})
    assert r.status_code == 200
    assert "pH 2.7" in _modes(r.json()["outputs"])["RP"]["help"]


# ---------------------------------------------------------------- setup
def test_setup_paths_and_pinned_files(monkeypatch, tmp_path):
    assert pka_setup.default_weights().name == "molgpka_weights.npz"
    assert pka_setup.default_weights().parent == Path.home() / ".mzdesk" / "models" / "molgpka"
    monkeypatch.setenv("MZDESK_MOLGPKA_DIR", str(tmp_path))
    assert pka_setup.default_weights() == tmp_path / "molgpka_weights.npz"
    assert len(pka_setup.COMMIT) == 40 and all(len(s) == 64 for s in pka_setup.FILES.values())
    assert list(pka_setup.FILES) == ["models/weight_acid.pth", "models/weight_base.pth"] and pka_setup.KINDS == ("acid", "base")


def test_setup_rejects_a_file_with_the_wrong_checksum(monkeypatch, tmp_path):
    monkeypatch.setenv("MZDESK_MOLGPKA_DIR", str(tmp_path))
    (tmp_path / "weight_acid.pth").write_bytes(b"not the weights")
    with pytest.raises(SystemExit, match="expected"):
        pka_setup.get_file("models/weight_acid.pth", pka_setup.FILES["models/weight_acid.pth"])
