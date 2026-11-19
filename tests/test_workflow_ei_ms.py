"""EI MS Annotation workflow (mzdesk/workflows/ei_ms_annotation.py, ei_candidates.py, ei_rerank.py).

The pure parts (parent names, candidate merging, the ranking, the levels, the cleanup, the Kovats RI) are tested directly. The whole workflow runs on a
tiny synthetic EI library (tests/fixtures/library_ei_ms.msp) indexed into a temporary MZDESK_LIBRARY_DIR, with NEIMS replaced by a small fake; one
test uses the real NEIMS weights when they are installed. Nothing here reads the real MassBank index or the internet. The tests check what the
workflow does with spectra and structures, not how often it is right.
"""
import time
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from mzdesk import spec
from mzdesk.api import app as app_module
from mzdesk.ms2 import library_index as lx
from mzdesk.ms2 import library_setup as setup_mod
from mzdesk.ms2 import neims_engine as ne
from mzdesk.platform_fit import derivatization as d
from mzdesk.workflows import ei_ms_annotation as ems
from mzdesk.workflows import ei_candidates as cd
from mzdesk.workflows import ei_rerank as rr

FIXTURE = Path(__file__).parent / "fixtures" / "library_ei_ms.msp"
client = TestClient(app_module.app)
ENTROPY = "SpectralEntropy (weighted)"

LACTIC_TMS = "45 30\n73 999\n75 150\n117 600\n118 60\n147 90\n191 20\n219 40"                    # = MSBNK-FIXGC-0001 (LACTIC ACID-DI-TMS, M 234)
GLYCINE_TMS = "73 999\n86 380\n100 210\n102 40\n147 300\n174 700\n175 140\n248 140\n276 45"          # = MSBNK-FIXGC-0002 (Glycine, 3 TMS spectrum, M 291)
UREA = "17 25\n28 60\n30 80\n43 120\n44 999\n60 700"                                              # = MSBNK-FIXGC-0003
HEXANOATE_TMS = "43 120\n73 999\n75 400\n117 500\n129 200\n145 80\n173 100"                      # = MSBNK-FIXGC-0004 (TRIMETHYLSILYL HEXANOATE, M 188)
LACTIC_SMILES, ALANINE, SUCCINIC, GLYCINE, PYRUVIC = "CC(O)C(=O)O", "C[C@H](N)C(=O)O", "OC(=O)CCC(O)=O", "NCC(=O)O", "CC(=O)C(=O)O"


def form(spectrum_text, **over):
    base = {"spectrum": {"source": "Paste", "value": spectrum_text}, "derivatization": "Auto (from marker ions)", "molecular_ion": "", "candidates": "",
            "alkanes": "", "sample_rt": ""}
    return {**base, **over}


def run(spectrum_text, progress=None, **over):
    return ems.run_query(ems.parse_inputs(form(spectrum_text, **over)), progress)


def steps_of(result):
    return {s["id"]: s for s in result["steps"]}


def table(result, block_id="candidates"):
    block = next(b for b in result["blocks"] if b["id"] == block_id)
    return [r["cells"] if isinstance(r, dict) else r for r in block["data"]["rows"]]


def text_of(cell):
    return cell if isinstance(cell, str) else cell.get("text") or cell.get("badge") or ""


def fake_peaks(formula, nominal):
    table_ = {"C9H22O3Si2": [(45, 25), (73, 999), (75, 120), (117, 650), (147, 100), (191, 15), (219, 50)],         # lactic acid 2TMS: close to LACTIC_TMS
              "C12H31NO2Si3": [(45, 30), (73, 999), (75, 150), (117, 600), (118, 60), (147, 90), (191, 20), (219, 40)],   # alanine 3TMS: the query's twin
              "C11H29NO2Si3": [(73, 999), (86, 350), (100, 200), (147, 320), (174, 720), (248, 150), (276, 40)],
              "C9H20O2Si": [(43, 100), (73, 999), (75, 380), (117, 520), (129, 180), (173, 90)],
              "CH4N2O": [(17, 20), (30, 70), (43, 100), (44, 999), (60, 650)]}
    return table_.get(formula, [(41, 100), (nominal - 15, 300), (nominal, 999)])


@pytest.fixture
def fake_neims(monkeypatch, tmp_path):
    """NEIMS replaced by a table of spectra keyed by the derivative's formula (deterministic, no weights needed)."""
    weights = tmp_path / "neims_weights.npz"
    weights.write_bytes(b"fake")
    calls = []

    def predict(smiles):
        mol = d.Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid or empty SMILES.")
        formula = d.parent_formula(mol)
        nominal = d.nominal_mass(formula)
        if nominal > ne.MAX_MZ:
            raise ValueError(f"The molecular ion m/z {nominal} is above the model's maximum (m/z {ne.MAX_MZ}).")
        calls.append(smiles)
        return {"smiles": smiles, "formula": d.formula_string(formula), "molecular_ion_mz": nominal, "model_mass": nominal,
                "peaks": [(float(m), float(i) / 9.99) for m, i in fake_peaks(d.formula_string(formula), nominal)], "warnings": []}
    monkeypatch.setattr(ne, "weights_path", lambda: str(weights))
    monkeypatch.setattr(ne, "predict", predict)
    return calls


@pytest.fixture(scope="module")
def ei_libs(tmp_path_factory):
    """The synthetic EI library indexed as MassBank into a temporary library directory (module-wide)."""
    mp = pytest.MonkeyPatch()
    mp.setenv("MZDESK_LIBRARY_DIR", str(tmp_path_factory.mktemp("ei_libs")))
    lx.clear_cache()
    setup_mod.build_from("MassBank", FIXTURE, {"url": "file://fixture", "version": "fixture-ei", "date": "2026-01-01"}, progress=lambda *_: None)
    yield
    mp.undo()
    lx.clear_cache()


@pytest.fixture
def no_libs(monkeypatch, tmp_path):
    monkeypatch.setenv("MZDESK_LIBRARY_DIR", str(tmp_path / "empty"))
    lx.clear_cache()
    yield
    lx.clear_cache()


@pytest.fixture
def names(monkeypatch):
    """Name lookup answered from a small table instead of PubChem."""
    from mzdesk.identity import resolve
    known = {"alanine": ALANINE, "succinic acid": SUCCINIC, "glycine": GLYCINE, "lactic acid": LACTIC_SMILES}

    def to_smiles(query, kind="Auto"):
        if query.strip().lower() not in known:
            raise resolve.ResolveError(f"'{query}' not found in PubChem.")
        return known[query.strip().lower()]
    monkeypatch.setattr(resolve, "to_smiles", to_smiles)


# ================================================================ pure parts
@pytest.mark.parametrize("name, parent", [
    ("TARTARIC ACID-TETRA-TMS", "TARTARIC ACID"), ("BENZOIC ACID MONO-TMS", "BENZOIC ACID"), ("ETHINYLESTRADIOL-3-TMS", "ETHINYLESTRADIOL"),
    ("LACTIC ACID-DI-TMS", "LACTIC ACID"), ("P-HYDROXYBENZOIC ACID TMS ETHER TMS ESTER", "P-HYDROXYBENZOIC ACID"),
    ("ERUCIC ACID TRIMETHYLSILYL ESTER", "ERUCIC ACID"), ("METHYLSUCCINIC ACID BIS(TRIMETHYLSILYL) ESTER", "METHYLSUCCINIC ACID"),
    ("ERGOSTEROL-TRIMETHYLSILYL-ETHER", "ERGOSTEROL"), ("Glycine, N,N-bis(trimethylsilyl)-, trimethylsilyl ester", "Glycine"),
    ("Lactic acid, 2TMS", "Lactic acid"), ("Alanine (3TMS)", "Alanine"), ("Citric acid, TMS derivative", "Citric acid"), ("Glucose 1MEOX 5TMS", "Glucose"),
    ("D-Glucose, 2,3,4,5,6-pentakis-O-(trimethylsilyl)-, methyloxime", "D-Glucose"), ("Valine, TBDMS derivative", "Valine"),
    ("TRIMETHYLSILYL HEXANOATE", "HEXANOIC ACID"), ("TRIMETHYLSILYL OCTADECANOATE", "OCTADECANOIC ACID"),
    ("METHYL 2-TRIMETHYLSILOXY-OCTADECANOATE", "METHYL 2-HYDROXY-OCTADECANOATE"), ("2-(TRIMETHYLSILYL)OXYDODECANOIC ACID TRIMETHYLSILYL ESTER", "2-HYDROXYDODECANOIC ACID"),
    ("1,3-DI-O-TRIMETHYLSILYL DECASPHINGANINE", "DECASPHINGANINE"), ("N-TRIMETHYLSILYLPHENYLACETYLGLYCINE TRIMETHYLSILYL ESTER", "PHENYLACETYLGLYCINE"),
    ("L-Alanine", "L-Alanine"), ("Caffeine", "Caffeine"), ("TARTARIC ACID", "TARTARIC ACID")])
def test_parent_name_removes_the_derivative_words(name, parent):
    assert cd.parent_name(name) == (parent, True)


@pytest.mark.parametrize("name, parent", [
    ("(2S,3S)-4-BENZYLOXY-1-(TERT-BUTYLDIMETHYLSILYLOXY)-2-METHOXYETHANE", "(2S,3S)-4-BENZYLOXY-1-HYDROXY-2-METHOXYETHANE"),      # the closing bracket goes with the group
    ("1-(TRIMETHYLSILYLOXY)BUTANE", "1-HYDROXYBUTANE"), ("(TRIMETHYLSILYL)OXYACETIC ACID", "HYDROXYACETIC ACID"),
    ("3,4-BIS(TRIMETHYLSILYL)OXY-BENZOIC ACID TRIMETHYLSILYL ESTER", "3,4-DIHYDROXY-BENZOIC ACID"),
    ("METHYL 9,10-BIS(TERT-BUTYLDIMETHYLSILOXY)OCTADEC-12(Z)-ENOATE", "METHYL 9,10-DIHYDROXYOCTADEC-12(Z)-ENOATE"),
    ("3,4,5-TRI(TRIMETHYLSILYL)OXYBENZOIC ACID", "3,4,5-TRIHYDROXYBENZOIC ACID"), ("2-(2,3-DITRIMETHYLSILYLOXYBUTOXY)ETHANOL", "2-(2,3-DIHYDROXYBUTOXY)ETHANOL")])
def test_parent_name_turns_a_silyloxy_group_into_hydroxy_with_balanced_brackets(name, parent):
    assert cd.parent_name(name) == (parent, True)


def test_tidy_name_keeps_stereo_descriptors_in_capitals():
    assert cd.tidy_name("(R)-(-)-2-OCTANOL") == "(R)-(-)-2-octanol" and cd.tidy_name("(2S,3S)-BUTANEDIOL") == "(2S,3S)-butanediol"
    assert cd.tidy_name("(4R*,5R*)-4-ETHYLHEXANOL") == "(4R*,5R*)-4-ethylhexanol" and cd.tidy_name("(E)-2-BUTENE") == "(E)-2-butene"
    assert cd.tidy_name("(SE)-FOO") == "(se)-foo" and cd.tidy_name("(12)-FOO") == "(12)-foo"                                     # not a descriptor


@pytest.mark.parametrize("name", ["TRIMETHYLSILYLSERINE", "ETHOXYTRIMETHYLSILANE", "N-OCTYLBIS(TRIMETHYLSILYL)AMINE", "L4", "TMS"])
def test_a_name_that_still_holds_derivative_words_is_not_trusted(name):
    assert cd.parent_name(name)[1] is False


def test_tidy_name_lowercases_capital_names_but_keeps_stereo_prefixes():
    assert cd.tidy_name("TARTARIC ACID") == "Tartaric acid" and cd.tidy_name("DL-ALANINE") == "DL-alanine" and cd.tidy_name("L-ALANINE") == "L-alanine"
    assert cd.tidy_name("Glycine") == "Glycine" and cd.tidy_name("beta-Alanine") == "beta-Alanine" and cd.tidy_name("") == ""
    assert cd.tidy_name("DHP") == "DHP" and cd.tidy_name("UREA") == "Urea"                                            # short acronyms stay


def hit(library, accession, name, smiles, score, matches=8, **extra):
    return {"library": library, "accession": accession, "name": name, "smiles": smiles, "inchikey": extra.pop("inchikey", ""), "matches": matches, "row": extra.pop("row", 0),
            "scores": {ENTROPY: score}}


def test_hits_become_one_candidate_per_original_compound():
    hits = [hit("MassBank", "A1", "LACTIC ACID-DI-TMS", "CC(O[Si](C)(C)C)C(=O)O[Si](C)(C)C", 0.95),
            hit("MassBank", "A2", "L-Lactic acid", "C[C@H](O)C(O)=O", 0.90, matches=9),
            hit("MassBank", "A3", "GLYCINE", "NCC(O)=O", 0.85),
            hit("MassBank", "A4", "Mystery", "", 0.80, inchikey="ABCDEFGHIJKLMN-UHFFFAOYSA-N"),
            hit("MassBank", "A5", "Some other", "", 0.75),
            hit("MassBank", "A6", "No score", "CCO", None)]
    found = cd.from_hits(hits, ENTROPY, limit=10)
    assert [c.name for c in found] == ["Lactic acid", "Glycine", "Mystery", "Some other"]          # the A6 hit has no score and is ignored
    lactic = found[0]
    assert lactic.n_spectra == 2 and lactic.library["accession"] == "A1" and lactic.library["entry_mode"] == "TMS" and lactic.library["score"] == 0.95
    assert d.Chem.MolToInchiKey(d.Chem.MolFromSmiles(lactic.parent_smiles))[:14] == "JVTAAEKCZFNVCJ" and lactic.sources == ["library"]
    assert found[1].library["entry_mode"] is None and found[1].parent_smiles == "NCC(O)=O"
    assert found[2].parent_smiles is None and found[2].key == "ABCDEFGHIJKLMN" and found[3].key == "some other"           # no structure: the InChIKey, else the name
    assert [c.name for c in cd.from_hits(hits, ENTROPY, limit=2)] == ["Lactic acid", "Glycine"]
    assert cd.from_hits(hits, ENTROPY, limit=1)[0].n_spectra == 2                                  # a limit never drops further spectra of a compound that is in


def test_a_derivative_name_without_a_structure_loses_its_derivative_words_and_joins_the_plain_entry():
    hits = [hit("MassBank", "B1", "TRANS-11-OCTADECENOIC ACID TRIMETHYLSILYLESTER", "", 0.90), hit("MassBank", "B2", "Trans-11-octadecenoic acid", "", 0.85),
            hit("MassBank", "B3", "TRIMETHYLSILYLSERINE", "", 0.80)]
    found = cd.from_hits(hits, ENTROPY)
    assert [c.name for c in found] == ["Trans-11-octadecenoic acid", "Trimethylsilylserine"] and found[0].n_spectra == 2 and found[0].library["accession"] == "B1"
    assert found[0].name_note == "" and found[1].name_note == "derivative name"                                    # an untrusted name stays as it is, with the note


def test_a_siloxane_or_silazane_entry_is_a_compound_of_its_own_not_water_or_ammonia():
    hits = [hit("MassBank", "S1", "HEXAMETHYLDISILOXANE", "C[Si](C)(C)O[Si](C)(C)C", 0.95), hit("MassBank", "S2", "HEXAMETHYLDISILAZANE", "C[Si](C)(C)N[Si](C)(C)C", 0.90),
            hit("MassBank", "S3", "TRIMETHYLSILANOL", "C[Si](C)(C)O", 0.85)]
    found = cd.from_hits(hits, ENTROPY)
    assert [c.parent_smiles for c in found] == ["C[Si](C)(C)O[Si](C)(C)C", "C[Si](C)(C)N[Si](C)(C)C", "C[Si](C)(C)O"] and len({c.key for c in found}) == 3
    assert all(c.library["entry_mode"] is None for c in found)


def test_the_best_spectrum_of_a_compound_is_kept():
    hits = [hit("MassBank", "A1", "Alanine", "CC(N)C(O)=O", 0.80, matches=5), hit("MassBank", "A2", "L-Alanine", "C[C@H](N)C(O)=O", 0.90, matches=3),
            hit("GNPS", "A3", "L-Alanine", "C[C@H](N)C(O)=O", 0.90, matches=7)]
    (c,) = cd.from_hits(hits, ENTROPY)
    assert (c.library["accession"], c.library["library"], c.n_spectra) == ("A3", "GNPS", 3)       # equal scores: more matched peaks win


def test_user_lines_are_smiles_or_names_and_a_derivative_is_taken_back_to_its_compound():
    lookup = {"alanine": ALANINE}

    def to_smiles(text):
        if text not in lookup:
            raise ValueError(f"'{text}' not found")
        return lookup[text]
    lines = ["CCO", "alanine", "nonsense name", "[NH3+]CC(=O)[O-]", "O=C(O[Si](C)(C)C)C(C)O[Si](C)(C)C", "CC(=O)[O-].[Na+]"]
    found, warnings = cd.from_lines(lines, to_smiles)
    assert [c.name for c in found] == ["CCO", "alanine", "[NH3+]CC(=O)[O-]", "CC(O)C(=O)O", "CC(=O)[O-]"]
    assert [c.user_text for c in found][:3] == ["CCO", "alanine", "[NH3+]CC(=O)[O-]"] and all(c.sources == ["list"] for c in found)
    assert found[1].parent_smiles == ALANINE and found[3].key == "JVTAAEKCZFNVCJ"
    assert len(warnings) == 2 and "“nonsense name” was skipped" in warnings[0] and "TMS derivative" in warnings[1]


def test_a_line_with_a_space_is_a_name_not_a_smiles_plus_a_comment():
    """RDKit reads "PC 33:1" as the SMILES "PC" (methylphosphine) with a name; the workflow must look it up as a name instead."""
    looked_up = []

    def to_smiles(text):
        looked_up.append(text)
        raise ValueError("not found")
    found, warnings = cd.from_lines(["PC 33:1", "CC(=O)O acetic acid", "CCO"], to_smiles)
    assert looked_up == ["PC 33:1", "CC(=O)O acetic acid"] and [c.name for c in found] == ["CCO"] and len(warnings) == 2


def test_merging_joins_the_sources_and_prefers_a_typed_name():
    from_library = cd.from_hits([hit("MassBank", "A1", "L-ALANINE", "C[C@H](N)C(O)=O", 0.9)], ENTROPY)
    typed, _ = cd.from_lines(["alanine", "CCO"], lambda text: ALANINE)
    merged = cd.merge(from_library, typed)
    assert [c.name for c in merged] == ["alanine", "CCO"]
    assert merged[0].sources == ["library", "list"] and merged[0].name == "alanine" and merged[0].library["accession"] == "A1" and merged[0].user_text == "alanine"
    only_smiles, _ = cd.from_lines([ALANINE], lambda text: text)
    kept = cd.merge(from_library, only_smiles)[0]
    assert kept.name == "L-alanine" and kept.sources == ["library", "list"]                         # a SMILES line does not rename the library compound


def candidate(name, tier, lib=None, neims=None, m=None, **extra):
    c = cd.Candidate(name=name, key=name, parent_smiles="CCO", tier=tier, lib_score=lib, neims_score=neims, m_status=m)
    c.library = {"score": lib, "matches": 8, "library": "MassBank", "accession": name, "smiles": "CCO", "entry_mode": None, "state": None, "row": 0, "name": name} if lib is not None else None
    for k, v in extra.items():
        setattr(c, k, v)
    return c


def test_ranking_puts_library_evidence_first_then_m_consistency_then_the_score():
    cs = [candidate("n_high", "neims", neims=0.99), candidate("lib_mismatch", "library", lib=0.95, neims=0.9, m="mismatch"), candidate("lib_low", "library", lib=0.82, neims=0.3),
          candidate("n_mid", "neims", neims=0.5), candidate("lib_mid", "library", lib=0.90, neims=0.2), candidate("n_none", "neims"),
          candidate("n_mismatch", "neims", neims=0.9, m="mismatch"), candidate("lib_partial", "library", lib=0.85, neims=0.4, m="partial"),
          candidate("lib_unknown", "library", lib=0.99, m="unknown")]
    before = [c.name for c in cs]
    ranked = rr.rerank(cs)
    assert [c.name for c in ranked] == ["lib_mid", "lib_partial", "lib_low", "lib_unknown", "lib_mismatch", "n_high", "n_mid", "n_none", "n_mismatch"]
    assert [c.rank for c in ranked] == list(range(1, 10)) and [c.name for c in cs] == before                       # the input list keeps its order
    assert rr.rerank([]) == []


def test_ranking_breaks_a_library_tie_with_the_neims_score_and_then_the_name():
    cs = [candidate("b", "library", lib=0.9, neims=0.4), candidate("a", "library", lib=0.9, neims=0.7), candidate("c", "library", lib=0.9, neims=0.7)]
    assert [c.name for c in rr.rerank(cs)] == ["a", "c", "b"]
    assert [c.name for c in rr.rerank([candidate("z", "library", lib=0.9), candidate("y", "library", lib=0.9, neims=0.1)])] == ["y", "z"]


def variant(nominal, n_silyl, mode="TMS", mass=None):
    return {"mode": mode, "nominal": nominal, "n_silyl": n_silyl, "mass": nominal + 0.1 if mass is None else mass}


def test_m_status_compares_the_derivative_with_the_typed_molecular_ion():
    assert rr.m_status(variant(291, 3), None) == (None, None) and rr.m_status(None, 291) == (None, None) and rr.m_status({"nominal": None}, 291) == (None, None)
    assert rr.m_status(variant(291, 3), 291) == ("ok", 3) and rr.m_status(variant(291, 3), 291.1508) == ("ok", 3)
    assert rr.m_status(variant(291, 3), 219) == ("partial", 2) and rr.m_status(variant(291, 3), 147) == ("partial", 1)         # 2 TMS, 1 TMS
    assert rr.m_status(variant(291, 3), 75) == ("mismatch", None)                                                              # the bare compound is not a derivative of it
    assert rr.m_status(variant(291, 3), 305) == ("mismatch", None) and rr.m_status(variant(291, 3), 290)[0] == "mismatch"
    assert rr.m_status(variant(399, 3, "TBDMS"), 285) == ("partial", 2) and rr.m_status(variant(75, 0, "None"), 291)[0] == "mismatch"
    assert rr.m_status(variant(569, 6, "MeOx + TMS"), 497) == ("partial", 5)
    assert rr.m_status(variant(300, 0, "None", mass=300.4), 300.2) == ("ok", 0)                                                  # a typed exact mass within 0.5 Da


def test_the_level_gate_is_the_esi_workflows():
    from mzdesk.annotation import library_step
    assert rr.gate() == (library_step.GATE_SIMILARITY, library_step.GATE_MIN_MATCHED, library_step.LIST_THRESHOLD) == (0.80, 4, 0.70)
    assert 0.0 < rr.NEIMS_MATCH_MIN < rr.gate()[0]


def level(c, n=30):
    return rr.ei_level(c, n)


def test_levels_follow_the_esi_scheme():
    from mzdesk.annotation.levels import level_label, level_name
    assert level(candidate("a", "library", lib=0.91))[0] == "2a" and level(candidate("a", "library", lib=0.80))[0] == "2a"       # the gate is inclusive
    lvl, reason = level(candidate("a", "library", lib=0.79))
    assert lvl == "3" and "score below 0.80" in reason
    assert level(candidate("a", "library", lib=0.95, library_matches=0))[0] == "2a"
    few = candidate("a", "library", lib=0.95)
    few.library["matches"] = 3
    assert level(few, 30)[0] == "3" and "only 3 matched peaks (needs 4)" in level(few, 30)[1] and level(few, 2)[0] == "2a" and level(few, 3)[0] == "2a"
    nostruct = candidate("a", "library", lib=0.95, parent_smiles=None)
    assert level(nostruct)[0] == "3" and "no structure" in level(nostruct)[1]
    assert level(candidate("a", "library", lib=0.95, state_ok=False))[0] == "3"
    assert level(candidate("a", "library", lib=0.95, m_status="mismatch"))[0] == "3" and level(candidate("a", "library", lib=0.95, m_status="partial"))[0] == "2a"
    assert level(candidate("a", "neims", neims=0.45))[0] == "3" and level(candidate("a", "neims", neims=rr.NEIMS_MATCH_MIN))[0] == "3"
    lvl, reason = level(candidate("a", "neims", neims=0.39))
    assert lvl is None and "below 0.40" in reason
    assert level(candidate("a", "neims", neims=0.9, m_status="mismatch"))[0] is None and level(candidate("a", "neims"))[0] is None
    assert level(candidate("a", "neims", neims=0.9, m_status="partial"))[0] == "3"
    assert (level_name("2a"), level_name("3")) == ("Level 2a", "Level 3") and level_label("2a").startswith("Level 2a: probable structure")


def test_there_is_no_level_one_or_h_in_the_ei_workflow():
    for c in (candidate("a", "library", lib=1.0), candidate("b", "neims", neims=1.0), candidate("c", "library", lib=0.7), candidate("d", "neims")):
        assert level(c)[0] in ("2a", "3", None)


# ---------------------------------------------------------------- step 1
def test_cleanup_drops_noise_but_never_the_molecular_ion():
    ei = [(73.0, 5000.0), (147.0, 1500.0), (200.0, 10.0), (291.0, 20.0)]
    res = ems.clean_query(ei, 291)
    assert res["peaks"] == [(73.0, 100.0), (147.0, 30.0), (291.0, pytest.approx(0.4))] and res["m_restored"] and res["unit_mass"]
    assert (res["n_in"], res["n_noise"], res["n_out"]) == (4, 2, 3)
    no_m = ems.clean_query(ei)
    assert [m for m, _ in no_m["peaks"]] == [73.0, 147.0] and not no_m["m_restored"]
    assert ems.clean_query([(73.0, 5000.0), (291.0, 200.0)], 291)["m_restored"] is False                 # the ion is already above the floor
    assert ems.NOISE_PCT == 0.5


def test_cleanup_merges_a_high_resolution_spectrum_to_unit_mass():
    res = ems.clean_query([(73.0468, 1000.0), (73.0470, 500.0), (147.0441, 300.0), (147.3, 10.0), (-5.0, 100.0), (88.0, 0.0)])
    assert not res["unit_mass"] and res["n_merged"] == 2 and res["peaks"][0] == (73.0, 100.0) and res["peaks"][1][0] == 147.0
    assert res["peaks"][1][1] == pytest.approx(310.0 / 1500.0 * 100.0)
    assert ems.clean_query([(73.0, 10.0), (147.0, 5.0)])["unit_mass"]                                    # integers (and 73.01) stay as they are
    assert ems.clean_query([(73.01, 10.0), (147.0, 5.0)])["unit_mass"]


# ---------------------------------------------------------------- step 5
def test_variants_keep_the_structure_of_a_library_derivative_of_the_same_kind():
    entry = "CC(O[Si](C)(C)C)C(=O)O[Si](C)(C)C"
    lactic = cd.from_hits([hit("MassBank", "A1", "LACTIC ACID-DI-TMS", entry, 0.95)], ENTROPY)[0]
    (own,) = ems.build_variants(lactic, ["TMS"], 234)
    assert own["own"] is True and own["smiles"] == entry and own["formula"] == "C9H22O3Si2" and (own["nominal"], own["n_silyl"], own["n_meox"]) == (234, 2, 0)
    assert (own["m_status"], own["m_k"]) == ("ok", 2)
    (tbdms,) = ems.build_variants(lactic, ["TBDMS"], 234)
    assert tbdms["own"] is False and tbdms["formula"] == "C15H34O3Si2" and tbdms["m_status"] == "mismatch"
    (none,) = ems.build_variants(lactic, ["None"], 234)
    assert none["smiles"] == lactic.parent_smiles and (none["n_silyl"], none["nominal"]) == (0, 90) and none["m_status"] == "mismatch"
    meox_tms = ems.build_variants(lactic, ["TMS", "MeOx + TMS"], None)
    assert [v["own"] for v in meox_tms] == [True, True] and meox_tms[1]["mode"] == "MeOx + TMS" and meox_tms[1]["m_status"] is None     # no C=O: the TMS entry is the MeOx + TMS form


def test_both_forms_are_built_for_a_compound_with_a_carbonyl():
    pyruvic, _ = cd.from_lines([PYRUVIC], lambda t: t)
    tms, meox = ems.build_variants(pyruvic[0], ["TMS", "MeOx + TMS"], 161)
    assert (tms["formula"], tms["nominal"], tms["m_status"]) == ("C6H12O3Si", 160, "mismatch")
    assert (meox["formula"], meox["nominal"], meox["n_silyl"], meox["n_meox"], meox["m_status"]) == ("C7H15NO3Si", 189, 1, 1, "mismatch")
    assert ems.build_variants(pyruvic[0], ["MeOx + TMS"], 189)[0]["m_status"] == "ok"
    glycine, _ = cd.from_lines([GLYCINE], lambda t: t)
    assert ems.build_variants(glycine[0], ["TMS"], 219)[0]["m_status"] == "ok" or ems.build_variants(glycine[0], ["TMS"], 219)[0]["m_status"] == "partial"
    assert ems.build_variants(glycine[0], ["TMS"], 219)[0]["m_k"] == 2
    assert ems.build_variants(cd.Candidate(name="x", key="x"), ["TMS"], 100) == []                       # no structure, no derivative


def test_derivative_text():
    assert ems.derivative_text({"n_silyl": 3, "n_meox": 0, "mode": "TMS"}) == "3TMS" and ems.derivative_text({"n_silyl": 5, "n_meox": 1, "mode": "MeOx + TMS"}) == "1MeOx + 5TMS"
    assert ems.derivative_text({"n_silyl": 2, "n_meox": 0, "mode": "TBDMS"}) == "2TBDMS" and ems.derivative_text({"n_silyl": 0, "n_meox": 0, "mode": "None"}) == "none"
    assert ems.derivative_text(None) == "—"


# ---------------------------------------------------------------- step 8
def test_kovats_ri_by_hand():
    q = ems.parse_inputs(form(LACTIC_TMS, alkanes="C10 10.0\nC8 8.0\nC9 9.1", sample_rt="9.5"))
    ri = ems.kovats(q)
    assert ri["ri"] == pytest.approx(900 + 100 * (9.5 - 9.1) / (10.0 - 9.1)) and ri["ri"] == pytest.approx(944.44, abs=0.01)           # van den Dool-Kratz between C9 and C10
    assert ri["lower"] == (9, 9.1) and ri["upper"] == (10, 10.0) and ri["extrapolated"] is False
    below = ems.kovats(ems.parse_inputs(form(LACTIC_TMS, alkanes="C8 8.0\nC9 9.1\nC10 10.0", sample_rt="7.0")))
    assert below["ri"] == pytest.approx(800 + 100 * (7.0 - 8.0) / (9.1 - 8.0)) and below["extrapolated"] is True
    assert ems.kovats(ems.parse_inputs(form(LACTIC_TMS))) is None
    on_alkane = ems.kovats(ems.parse_inputs(form(LACTIC_TMS, alkanes="C8 8.0\nC9 9.1", sample_rt="9.1")))
    assert on_alkane["ri"] == pytest.approx(900.0)


# ================================================================ the whole workflow
def test_a_derivative_entry_of_the_library_is_named_as_its_original_compound(ei_libs, fake_neims):
    r = run(LACTIC_TMS, molecular_ion="234")
    a = r["answer"]
    assert a["statement"].startswith("Lactic acid (as 2TMS derivative), Level 2a: MassBank EI match MSBNK-FIXGC-0001, score 1.00")
    assert a["level"] == "Level 2a" and a["level_label"] == "Level 2a: probable structure (library spectrum match)"
    facts = dict(a["facts"])
    assert facts["Derivatization"] == "TMS (inferred; MeOx undecided)" and facts["TMS groups"] == "2 (m/z 147 suggests 2 or more)"
    assert facts["Molecular ion"] == "234 (M-15 seen)" and facts["Top library match"].startswith("LACTIC ACID-DI-TMS (MassBank MSBNK-FIXGC-0001), score 1.00")
    assert "Kovats RI" not in facts
    rows = table(r)
    assert text_of(rows[0][1]) == "Lactic acid" and text_of(rows[0][2]) == "2TMS · M 234" and rows[0][3] == "EI library" and rows[0][6]["badge"] == "2a"
    assert [s["status"] for s in r["steps"]] == ["done"] * 7 + ["skipped", "done"]
    assert [b["id"] for b in r["blocks"]] == ["candidates", "mirror"] and r["warnings"] == []


def test_a_hit_without_a_structure_is_listed_but_cannot_reach_level_2a(ei_libs, fake_neims):
    r = run(LACTIC_TMS, molecular_ion="234")
    rows = {text_of(row[1]): row for row in table(r)}
    mystery = rows["Unknown silyl compound"]
    assert mystery[3] == "EI library" and mystery[6]["badge"] == "3" and "no structure" in mystery[7] and mystery[2] == "—"          # no derivative without a structure
    assert [text_of(row[1]) for row in table(r)].index("Lactic acid") < [text_of(row[1]) for row in table(r)].index("Unknown silyl compound")


def test_a_parent_structure_entry_with_a_tms_spectrum_gets_its_derivative_from_the_rule(ei_libs, fake_neims):
    r = run(GLYCINE_TMS, molecular_ion="291")
    assert r["answer"]["statement"].startswith("Glycine (as 3TMS derivative), Level 2a: MassBank EI match MSBNK-FIXGC-0002")
    top = table(r)[0]
    assert text_of(top[2]) == "3TMS · M 291" and top[7] == "—" and dict(r["answer"]["facts"])["Molecular ion"] == "291 (M-15 seen)"
    partial = run(GLYCINE_TMS, molecular_ion="219")
    top = table(partial)[0]
    assert text_of(top[1]) == "Glycine" and "partial silylation (2 of 3)" in top[7] and partial["answer"]["level"] == "Level 2a"


def test_an_underivatized_spectrum_is_named_without_a_derivative(ei_libs, fake_neims):
    r = run(UREA)
    assert r["answer"]["statement"].startswith("Urea (underivatized), Level 2a: MassBank EI match MSBNK-FIXGC-0003")
    facts = dict(r["answer"]["facts"])
    assert facts["Derivatization"] == "None (inferred)" and "TMS groups" not in facts and "Molecular ion" not in facts
    assert text_of(table(r)[0][2]) == "none"


def test_a_prefix_style_derivative_name_is_turned_into_the_acid(ei_libs, fake_neims):
    r = run(HEXANOATE_TMS, molecular_ion="188")
    assert r["answer"]["statement"].startswith("Hexanoic acid (as 1TMS derivative), Level 2a: MassBank EI match MSBNK-FIXGC-0004")
    assert {text_of(row[1]) for row in table(r)} >= {"Hexanoic acid"} and len([row for row in table(r) if text_of(row[1]) == "Hexanoic acid"]) == 1     # the plain entry merges in


def by_name(result, name, block_id="candidates"):
    """The row of the candidate table named ``name`` (case does not matter; a typed name keeps the case it was typed in)."""
    return next(row for row in table(result, block_id) if text_of(row[1]).lower() == name.lower())


def test_a_typed_derivatization_that_clashes_with_the_markers_is_kept_with_a_warning(ei_libs, fake_neims):
    r = run(LACTIC_TMS, derivatization="None")
    assert any("You chose None" in w and "Your choice is kept" in w for w in r["warnings"])
    top = by_name(r, "Lactic acid")
    assert text_of(top[2]) == "none" and "reference TMS" in top[7]
    assert top[6]["badge"] == "3" and table(r)[0][0] == "1" and text_of(table(r)[0][1]) == "Lactic acid" and r["answer"]["level"] == "Level 3"
    assert "another derivatization state" in r["answer"]["statement"] and dict(r["answer"]["facts"])["Derivatization"] == "None (selected)"
    ok = run(LACTIC_TMS, derivatization="TMS", molecular_ion="234")
    assert ok["answer"]["level"] == "Level 2a" and not any("You chose" in w for w in ok["warnings"])


def test_a_molecular_ion_that_contradicts_the_derivative_demotes_the_hit(ei_libs, fake_neims):
    r = run(LACTIC_TMS, molecular_ion="250")
    lactic = by_name(r, "Lactic acid")
    assert "M mismatch (234 vs 250)" in lactic[7] and lactic[6]["badge"] == "3"
    order = [text_of(row[1]) for row in table(r)]
    assert order.index("Unknown silyl compound") < order.index("Lactic acid")                                    # nothing to compare beats a contradiction
    assert not any("M mismatch" in row[7] for row in table(r) if text_of(row[1]) == "Unknown silyl compound")
    ok = run(LACTIC_TMS, molecular_ion="234")
    assert text_of(table(ok)[0][1]) == "Lactic acid"


def test_user_candidates_are_ranked_under_the_library_hits(ei_libs, fake_neims, names):
    r = run(LACTIC_TMS, candidates=f"{ALANINE}\nsuccinic acid\nlactic acid\nnot in pubchem")
    rows = table(r)
    tiers = [row[3] for row in rows]                     # the library hits first (how many clear the threshold depends on the ms_entropy version:
    assert tiers[:2] == ["EI library", "EI library"] and tiers[-2:] == ["NEIMS only", "NEIMS only"]          # 1.5 scores the hexanoate 0.71, 1.3 0.69)
    assert tiers == sorted(tiers, key=lambda t: t != "EI library")
    assert [text_of(row[1]) for row in rows if row[3] == "NEIMS only"] == [ALANINE, "succinic acid"]               # "not in pubchem" gives no row
    assert text_of(rows[0][1]).lower() == "lactic acid" and text_of(rows[0][1]) == "lactic acid"               # the typed name is kept
    assert sum(text_of(row[1]).lower() == "lactic acid" for row in rows) == 1 and rows[0][6]["badge"] == "2a"       # typing lactic acid again makes no second row
    alanine, succinic = rows[-2], rows[-1]
    assert text_of(alanine[1]) == ALANINE and alanine[4] == "—" and float(alanine[5]) == pytest.approx(1.0, abs=0.01) and alanine[6]["badge"] == "3"
    assert text_of(succinic[1]) == "succinic acid" and succinic[6] == "—"
    assert any("“not in pubchem” was skipped" in w for w in r["warnings"])
    assert r["answer"]["statement"].startswith("lactic acid (as 2TMS derivative), Level 2a")


def test_a_candidate_that_only_neims_supports_reaches_level_3_and_a_wrong_one_gets_none(ei_libs, fake_neims, names):
    r = run(LACTIC_TMS, candidates=f"{ALANINE}\n{SUCCINIC}")
    alanine, succinic = by_name(r, ALANINE), by_name(r, SUCCINIC)
    assert alanine[6]["badge"] == "3" and float(alanine[5]) > 0.9 and succinic[6] == "—" and float(succinic[5]) < rr.NEIMS_MATCH_MIN
    typed_m = run(LACTIC_TMS, candidates=ALANINE, molecular_ion="234")                                          # alanine 3TMS is 305, not 234
    row = by_name(typed_m, ALANINE)
    assert "M mismatch (305 vs 234)" in row[7] and row[6] == "—"


def test_without_the_library_the_user_candidates_are_still_ranked(no_libs, fake_neims, names):
    r = run(LACTIC_TMS, candidates=f"{ALANINE}\nlactic acid\n{SUCCINIC}")
    s = steps_of(r)
    assert s[2]["status"] == "skipped" and "No EI library is installed" in s[2]["note"] and s[4]["status"] == "done" and s[7]["status"] == "done"
    assert any("No EI library spectra were searched" in w for w in r["warnings"])
    rows = table(r)
    assert [row[3] for row in rows] == ["NEIMS only"] * 3
    assert [text_of(row[1]) for row in rows][:2] == [ALANINE, "lactic acid"] or [text_of(row[1]) for row in rows][:2] == ["lactic acid", ALANINE]
    assert float(rows[0][5]) >= float(rows[1][5]) >= float(rows[2][5])
    assert r["answer"]["level"] == "Level 3" and "NEIMS-predicted spectrum match only" in r["answer"]["statement"]
    assert dict(r["answer"]["facts"])["Top library match"] == "none"


def test_a_labelled_candidate_gets_one_note_about_its_derivative_mass(no_libs, fake_neims, names):
    """The derivative formula and mass of step 5 come from the unlabelled formula (the NEIMS spectrum of step 6 uses the labelled structure)."""
    glycine_13c, glycine_d2 = "N[13CH2]C(=O)O", "[2H]C([2H])(N)C(=O)O"

    def notes(candidates):
        return [w for w in run(GLYCINE_TMS, candidates=candidates)["warnings"] if "isotope labels" in w]

    assert notes(glycine_13c) == [f"Candidate “{glycine_13c}” has isotope labels. The derivative formulas and masses are calculated for the unlabelled compound; "
                                  "use Calculator › Isotope labelling for labelled masses."]
    two = notes(f"{glycine_13c}\n{ALANINE}\n{glycine_d2}")                       # one note for the list, naming the labelled lines only
    assert len(two) == 1 and f"Candidates “{glycine_13c}”, “{glycine_d2}” have isotope labels." in two[0] and ALANINE not in two[0]
    assert notes(f"{GLYCINE}\n{ALANINE}") == [] and notes("") == []


def test_without_the_neims_weights_the_library_evidence_alone_ranks(ei_libs, monkeypatch, tmp_path, names):
    monkeypatch.setattr(ne, "weights_path", lambda: str(tmp_path / "missing.npz"))
    r = run(LACTIC_TMS, molecular_ion="234", candidates=ALANINE)
    s = steps_of(r)
    assert s[6]["status"] == "skipped" and "weights are not installed" in s[6]["note"] and s[7]["status"] == "done"
    assert any("NEIMS weights are not installed" in w and "neims_setup" in w for w in r["warnings"])
    lactic = by_name(r, "Lactic acid")
    assert table(r)[0] == lactic and lactic[5] == "—" and "NEIMS n/a" in lactic[7] and r["answer"]["level"] == "Level 2a"
    alanine = by_name(r, ALANINE)
    assert alanine[6] == "—" and "NEIMS n/a" in alanine[7]


def test_a_neims_failure_for_one_candidate_does_not_stop_the_others(ei_libs, fake_neims, monkeypatch, names):
    real_predict = ne.predict

    def predict(smiles):
        if "Si" not in smiles:
            raise ValueError("The molecular ion m/z 1200 is above the model's maximum (m/z 999).")
        return real_predict(smiles)
    monkeypatch.setattr(ne, "predict", predict)
    r = run(UREA, candidates="succinic acid", derivatization="None")
    assert steps_of(r)[6]["status"] == "done" and "of" in steps_of(r)[6]["note"]
    assert all("NEIMS n/a" in row[7] for row in table(r)) and r["answer"]["level"] == "Level 2a"
    assert text_of(table(r)[0][1]) == "Urea" and any("succinic acid" in text_of(row[1]) for row in table(r))


@pytest.mark.parametrize("step", ["_candidates", "_rerank"])
def test_a_failed_step_is_named_in_the_answer_instead_of_blaming_the_library(ei_libs, fake_neims, monkeypatch, step):
    def boom(*args, **kwargs):
        raise RuntimeError("injected")
    monkeypatch.setattr(ems, step, boom)
    r = run(LACTIC_TMS, molecular_ion="234")
    title = {"_candidates": "Candidates", "_rerank": "Similarity and rerank"}[step]
    assert r["answer"]["statement"] == f"No candidate. A step failed ({title}); see the warnings." and r["answer"]["level"] is None
    assert any("failed: RuntimeError: injected" in w for w in r["warnings"]) and [s["status"] for s in r["steps"]].count("failed") == 1


def test_a_candidate_whose_derivative_could_not_be_built_is_not_called_underivatized(ei_libs, fake_neims, monkeypatch):
    def boom(smiles, reagent):
        raise ValueError("cannot derivatize")
    monkeypatch.setattr(d, "derivatize", boom)
    r = run(GLYCINE_TMS, molecular_ion="291")                       # glycine is a parent-structure entry: its derivative must come from derivatize()
    top = table(r)[0]
    assert text_of(top[1]) == "Glycine" and "derivative failed" in top[7] and top[2] == "—"
    assert r["answer"]["statement"].startswith("Glycine (derivative structure not built)") and "underivatized" not in r["answer"]["statement"]


def test_a_spectrum_alone_still_gives_a_valid_result(no_libs, monkeypatch, tmp_path):
    monkeypatch.setattr(ne, "weights_path", lambda: str(tmp_path / "missing.npz"))
    r = run(LACTIC_TMS)
    assert set(r) == {"pipeline", "steps", "answer", "blocks", "warnings"} and r["pipeline"] == "ei_ms_annotation"
    assert r["answer"]["statement"].startswith("No candidate.") and r["answer"]["level"] is None and r["blocks"] == []
    s = steps_of(r)
    assert [s[n]["status"] for n in range(1, 10)] == ["done", "skipped", "done", "done", "done", "skipped", "done", "skipped", "done"]
    assert dict(r["answer"]["facts"])["Derivatization"] == "TMS (inferred; MeOx undecided)" and dict(r["answer"]["facts"])["Top library match"] == "none"
    assert len(r["warnings"]) == 1 and "No EI library spectra were searched" in r["warnings"][0]


def test_a_search_with_no_hit_over_the_list_threshold_is_not_an_error(ei_libs, fake_neims):
    r = run("20 100\n21 50\n22 30\n23 20")
    assert steps_of(r)[2]["status"] == "done" and "no library spectrum reaches" in steps_of(r)[2]["note"] and r["answer"]["statement"].startswith("No candidate.")


def test_the_kovats_ri_is_reported_in_the_answer_and_the_step(ei_libs, fake_neims):
    r = run(LACTIC_TMS, molecular_ion="234", alkanes="C8 8.0\nC9 9.1\nC10 10.0", sample_rt="9.5")
    s = steps_of(r)[8]
    assert s["status"] == "done" and s["note"] == "RI 944.4"
    assert dict(r["answer"]["facts"])["Kovats RI"] == "944.4"
    assert r["answer"]["statement"].startswith("Lactic acid")                                          # the RI is reported only: it does not enter the ranking
    extrapolated = run(LACTIC_TMS, alkanes="C8 8.0\nC9 9.1", sample_rt="12.0")
    assert dict(extrapolated["answer"]["facts"])["Kovats RI"].endswith("(extrapolated)") and any("extrapolated" in w for w in extrapolated["warnings"])


def test_the_mirror_plot_shows_the_library_spectrum_or_else_the_neims_prediction(ei_libs, fake_neims):
    r = run(LACTIC_TMS, molecular_ion="234")
    mirror = next(b for b in r["blocks"] if b["id"] == "mirror")
    assert mirror["type"] == "figure" and mirror["label"] == "Query vs Lactic acid (MassBank MSBNK-FIXGC-0001)" and "plotly" in mirror["data"] and mirror["data"]["peaks"]


def test_a_neims_only_top_candidate_is_drawn_against_its_prediction(no_libs, fake_neims, names):
    r = run(LACTIC_TMS, candidates=ALANINE)
    mirror = next(b for b in r["blocks"] if b["id"] == "mirror")
    assert "NEIMS prediction" in mirror["label"]


def run_collecting(spectrum_text, **over):
    """-> (result, {block id: block}) with every block the steps handed to the progress listener."""
    got = {}
    result = run(spectrum_text, progress=lambda n, status, note, blocks: got.update({b["id"]: b for b in (blocks or [])}), **over)
    return result, got


def test_both_reagents_are_tried_for_a_compound_with_a_carbonyl_when_meox_is_undecided(ei_libs, fake_neims, names):
    r, blocks = run_collecting(LACTIC_TMS, candidates=f"{PYRUVIC}\n{GLYCINE}")                              # Auto, TMS inferred, M not given: MeOx undecided
    rows = [row[:3] for row in (x["cells"] if isinstance(x, dict) else x for x in blocks["derivatives"]["data"]["rows"])]
    pyruvic = [row for row in rows if text_of(row[0]) == PYRUVIC]
    assert [row[1] for row in pyruvic] == ["TMS", "MeOx + TMS"] and [text_of(row[2]) for row in pyruvic] == ["1TMS", "1MeOx + 1TMS"]
    assert len([row for row in rows if text_of(row[0]) == "Glycine"]) == 1                                # no carbonyl: one form (named from the library)
    fixed, blocks = run_collecting(LACTIC_TMS, candidates=PYRUVIC, derivatization="TMS")
    assert [x["cells"][1] if isinstance(x, dict) else x[1] for x in blocks["derivatives"]["data"]["rows"] if text_of((x["cells"] if isinstance(x, dict) else x)[0]) == PYRUVIC] == ["TMS"]
    meox, blocks = run_collecting(LACTIC_TMS, candidates=PYRUVIC, derivatization="MeOx + TMS")
    assert [text_of((x["cells"] if isinstance(x, dict) else x)[2]) for x in blocks["derivatives"]["data"]["rows"] if text_of((x["cells"] if isinstance(x, dict) else x)[0]) == PYRUVIC] == ["1MeOx + 1TMS"]


def test_the_variant_that_fits_the_typed_molecular_ion_is_the_one_kept(ei_libs, fake_neims):
    r = run(LACTIC_TMS, candidates=PYRUVIC, molecular_ion="189")                                         # pyruvic acid MeOx + TMS is 189; TMS alone is 160
    row = by_name(r, PYRUVIC)
    assert text_of(row[2]).startswith("1MeOx + 1TMS") and "M mismatch" not in row[7]

def test_the_step_blocks_reach_the_progress_listener(ei_libs, fake_neims):
    events = []
    run(LACTIC_TMS, progress=lambda *a: events.append(a), molecular_ion="234")
    done = {e[0]: e for e in events if e[1] == "done"}
    ids = {n: [b["id"] for b in (e[3] or [])] for n, e in done.items()}
    assert ids[1] == ["query_spectrum", "cleanup_notes"] and ids[2] == ["library_hits", "library_notes"] and ids[3] == ["markers", "marker_notes"]
    assert ids[4] == ["candidate_pool"] and ids[5] == ["derivatives"] and ids[6] == ["neims"] and ids[7] == ["candidates"] and ids[9] == ["answer_notes"]
    assert [e[0] for e in events if e[1] == "skipped"] == [8]
    every = [b for e in events for b in (e[3] or [])]
    assert all(b["id"] and b["type"] and isinstance(b["data"], dict) for b in every) and len({b["id"] for b in every}) == len(every)       # every block has a unique id


def test_the_candidate_table_shape(ei_libs, fake_neims):
    r = run(LACTIC_TMS, molecular_ion="234")
    block = next(b for b in r["blocks"] if b["id"] == "candidates")
    assert block["columns"] == ["Rank", "Compound", "Derivative", "Evidence", "Library score", "NEIMS score", "Level", "Flags"] and block["row_help"] is True
    assert block["data"]["count"] == len(block["data"]["rows"]) and block["title"] == "Candidates"
    first = block["data"]["rows"][0]
    assert first["cells"][0] == "1" and first["cells"][4]["tone"] == "good" and any("EI library spectrum match" in line for line in first["help"])
    assert "Level 2a" in " ".join(first["help"])
    no_accuracy = " ".join([r["answer"]["statement"], *[str(h) for h in block.get("help", [])]]).lower()
    assert "accura" not in no_accuracy and "%" not in no_accuracy                                       # the results make no accuracy claim


# ---------------------------------------------------------------- through the job API
def wait_done(job_id, timeout=120):
    end = time.time() + timeout
    while time.time() < end:
        body = client.get(f"/api/pipelines/jobs/{job_id}").json()
        if body["status"] not in ("running", "queued"):
            return body
        time.sleep(0.05)
    raise AssertionError("job did not finish")


def test_the_job_api_runs_the_workflow_and_streams_the_step_blocks(ei_libs, fake_neims):
    app_module._jobs.clear()
    r = client.post("/api/pipelines/ei_ms_annotation/run", json=form(GLYCINE_TMS, molecular_ion="291", alkanes="C8 8.0\nC10 10.0", sample_rt="9.0"))
    assert r.status_code == 200
    body = wait_done(r.json()["job_id"])
    assert body["status"] == "done" and body["error"] is None
    result = body["result"]
    assert result["answer"]["statement"].startswith("Glycine (as 3TMS derivative), Level 2a") and dict(result["answer"]["facts"])["Kovats RI"] == "900.0"
    assert [s["status"] for s in body["steps"]] == ["done"] * 9 and all(s["blocks"] for s in body["steps"])
    assert body["steps"][1]["blocks"][0]["id"] == "library_hits" and body["steps"][6]["blocks"][0]["id"] == "candidates"
    assert [b["id"] for b in result["blocks"]] == ["candidates", "mirror"]
    app_module._jobs.clear()


def test_the_job_api_rejects_a_bad_form_before_starting_a_job():
    app_module._jobs.clear()
    assert client.post("/api/pipelines/ei_ms_annotation/run", json={}).status_code == 400
    assert not app_module._jobs


# ---------------------------------------------------------------- the spec
def test_the_examples_are_real_spectra_the_workflow_accepts():
    pipeline = next(p for p in spec.PIPELINES if p["id"] == "ei_ms_annotation")
    assert [e["label"] for e in pipeline["examples"]] == ["Glycine, 3 TMS (MassBank MSBNK-RIKEN-PR010014)",
                                                          "Dipentyl phthalate, high resolution (MassBank MSBNK-NILU-NL0043)"]
    for example in pipeline["examples"]:
        q = ems.parse_inputs({**{f["id"]: f.get("value") for f in pipeline["inputs"] if f.get("id")}, **example["values"]})
        assert len(q.ei) >= 12 and q.derivatization is None
    first = ems.parse_inputs({**{f["id"]: f.get("value") for f in pipeline["inputs"] if f.get("id")}, **pipeline["examples"][0]["values"]})
    assert first.molecular_ion == 291.0 and (73.0, 999.0) in first.ei and (174.0, 750.0) in first.ei


def test_the_candidate_field_declares_its_online_lookup():
    pipeline = next(p for p in spec.PIPELINES if p["id"] == "ei_ms_annotation")
    field = next(f for f in pipeline["inputs"] if f.get("id") == "candidates")
    assert field["network_names"] is True and field["type"] == "textarea"


# ---------------------------------------------------------------- real NEIMS (skipped when the weights are not installed)
@pytest.mark.skipif(not Path(ne.weights_path()).is_file(), reason="NEIMS weights are not installed")
def test_real_neims_ranks_the_candidate_whose_derivative_made_the_spectrum(no_libs):
    der = d.derivatize(ALANINE, d.EI_MODES["TMS"])
    predicted = ne.predict(der["smiles"])["peaks"]
    text = "\n".join(f"{mz} {rel * 9.99 * (0.85 + 0.3 * ((mz * 7) % 10) / 10):.1f}" for mz, rel in predicted)            # the prediction, each peak jittered by up to 15 %
    r = run(text, candidates=f"{GLYCINE}\nC[C@H](O)C(=O)O\n{ALANINE}\nOC(=O)CCC(O)=O\nN[C@@H](CO)C(=O)O", molecular_ion=str(der["mass"] and d.nominal_mass(der["formula"])))
    rows = table(r)
    assert text_of(rows[0][1]) == ALANINE and rows[0][3] == "NEIMS only" and text_of(rows[0][6]) == "3" and float(rows[0][5]) > 0.8
    assert r["answer"]["statement"].startswith(f"{ALANINE} (as 3TMS derivative), Level 3: NEIMS-predicted spectrum match only")
    assert steps_of(r)[6]["status"] == "done" and steps_of(r)[2]["status"] == "skipped"
