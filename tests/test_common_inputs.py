"""api.common input handling: the warnings a pasted peak table gets (decimal / thousands commas, rows with more than two numbers) and the note
for a structure with isotope labels, which only the pages whose outputs ignore the labels show (once each)."""
import copy
import shutil
from pathlib import Path

import pytest
from fastapi.testclient import TestClient
from rdkit import Chem

from mzdesk import spec
from mzdesk.api import app as app_module
from mzdesk.api import common
from mzdesk.ms2 import cfmid_engine as cf
from mzdesk.shared.spectrum import NOTE_DECIMAL_COMMA, NOTE_EXTRA_NUMBERS, NOTE_THOUSANDS
from mzdesk.workflows import targeted_method as tm

client = TestClient(app_module.app)
GLUCOSE = spec.EXAMPLE_STRUCTURES["D-glucose"]
GLUCOSE_13C6 = "O[13CH2][13C@H]1O[13CH](O)[13C@H](O)[13C@@H](O)[13C@@H]1O"
CAFFEINE_D3 = "[2H]C([2H])([2H])n1cnc2c1c(=O)n(C)c(=O)n2C"


# ================================================================ peak tables
def read(text, label="MS2 spectrum"):
    warnings = []
    peaks, _ = common.spectrum({"ms2": {"source": "Paste", "value": text}}, "ms2", label, warnings)
    return peaks, warnings


def test_a_peak_table_with_decimal_commas_is_read_and_flagged():
    two = [(138.0662, 100.0), (195.0877, 32.0)]
    for text in ("138,0662\t100\n195,0877\t32", "138,0662 100\n195,0877 32", "138,0662;100\n195,0877;32", "m/z;intensity\n138,0662;100\n195,0877;32"):
        assert read(text) == (two, [f"MS2 spectrum: {NOTE_DECIMAL_COMMA}"]), text


def test_a_semicolon_csv_is_read_as_mz_and_intensity_without_a_warning():
    assert read("m/z;intensity\n195.0877;100\n138.0662;60") == ([(138.0662, 60.0), (195.0877, 100.0)], [])
    assert read("195.0877;100\n138.0662;60") == ([(138.0662, 60.0), (195.0877, 100.0)], [])


def test_thousands_commas_and_extra_numbers_are_flagged_not_silent():
    assert read("100.0 1,200") == ([(100.0, 1200.0)], [f"MS2 spectrum: {NOTE_THOUSANDS}"])
    assert read("195.0877 100 5\n138.0662 60 3", "Query spectrum")[1] == [f"Query spectrum: {NOTE_EXTRA_NUMBERS}"]


def test_integer_mz_tables_for_ei_need_no_warning():
    three = [(41.0, 100.0), (55.0, 85.0), (69.0, 60.0)]
    for text in ("41,100\n55,85\n69,60", "41\t100\n55\t85\n69\t60", "41;100\n55;85\n69;60", "41 100\n55 85\n69 60", "41,100;55,85;69,60"):
        assert read(text) == (three, []), text


def test_other_formats_give_no_warning():
    mgf = "BEGIN IONS\nPEPMASS=195.0877 1000\nCHARGE=1+\n138.0662 100\n110.0713 18\nEND IONS"
    assert read(mgf) == ([(110.0713, 18.0), (138.0662, 100.0)], [])
    assert read("[[138.0662, 100], [110.0713, 18]]")[1] == [] and read("138.0662:100 110.0713:18")[1] == []


def test_a_page_reads_a_decimal_comma_table_like_the_dot_version():
    """Spectrum cleaning of a table pasted from Excel in a German or French locale (decimal commas, tab separated) against the same table with dots."""
    dots = "63.9626 319\n92.0507 63\n156.0129 999\n157.0136 15\n157.0166 54\n158.0089 20\n213.8606 11\n253.0486 37"
    commas = "\n".join(line.replace(".", ",", 1).replace(" ", "\t") for line in dots.splitlines())
    assert commas.startswith("63,9626\t319")

    def run(text):
        body = {"spectrum": {"source": "Paste", "value": text}, "precursor_mz": "252.0448", "ion_mode": "Negative", "noise_pct": "1",
                "normalize": "Base peak = 100", "remove_precursor": True}
        r = client.post("/api/modules/ms2_clean/run", json=body)
        assert r.status_code == 200, r.text
        return r.json()

    plain, decimal = run(dots), run(commas)
    assert decimal["outputs"]["peaks"] == plain["outputs"]["peaks"] and decimal["outputs"]["summary"] == plain["outputs"]["summary"]      # the same 8 peaks
    assert [w for w in decimal["warnings"] if "Decimal commas" in w] == [f"MS2 spectrum: {NOTE_DECIMAL_COMMA}"]
    assert plain["warnings"] == []


# ================================================================ isotope labels
# A page says so only when it shows outputs that ignore the labels, and says it in its own words. Facts behind the table (13C6-glucose, d3-caffeine,
# each against the unlabelled compound):
#   follow the labels   structure_formats (SMILES / InChI / InChIKey / MOL block / 3D), properties (monoisotopic mass, MW), separation and non_ms (MW in the
#                       rule texts), ionization (MW, NEIMS), ms2_ei (the whole spectrum), identity (looks up the isotopologue record)
#   unlabelled          adducts (ion formulas, m/z), ms1_theory (ion formula, m/z, isotope pattern), derivatization and platform_fit (derivative formula and
#                       mass; platform_fit's MW follows the labels), exclusivity (formula, theoretical m/z), ms2_fragments (fragment ions; the theoretical
#                       precursor follows the labels), ms2_esi (precursor m/z, ion formula)
#   not worth a note    ccs and known_methods show no mass or m/z: they find the reference CCS / published methods of the unlabelled compound (the first
#                       InChIKey block), which is what a chemist wants
LABELLED = {      # module -> (labelled structure, notes expected)
    "identity": (GLUCOSE_13C6, 0), "structure_formats": (GLUCOSE_13C6, 0), "properties": (GLUCOSE_13C6, 0), "platform_fit": (GLUCOSE_13C6, 1),
    "separation": (CAFFEINE_D3, 0), "ionization": (CAFFEINE_D3, 0), "adducts": (CAFFEINE_D3, 1), "non_ms": (CAFFEINE_D3, 0),
    "derivatization": (GLUCOSE_13C6, 1), "known_methods": (GLUCOSE_13C6, 0), "ms1_theory": (CAFFEINE_D3, 1), "ms2_ei": (CAFFEINE_D3, 0),
    "ms2_esi": (CAFFEINE_D3, 1), "ms2_fragments": (CAFFEINE_D3, 1), "exclusivity": (CAFFEINE_D3, 1), "ccs": (GLUCOSE_13C6, 0),
}
STRUCTURE_PAGES = [m["id"] for m in spec.MODULES if not spec.is_deferred(m) and any(f["type"] == "structure" for f in m["inputs"])]
HAS_JAVA = shutil.which("java") is not None


def unlabelled(smiles):
    mol = Chem.MolFromSmiles(smiles)
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)
    return Chem.MolToSmiles(Chem.RemoveHs(mol))


def run_page(module, smiles, monkeypatch):
    """The page's example with the structure replaced; the parts that need the internet, Docker or a local engine are stubbed or dropped."""
    from mzdesk.identity import enrich, resolve
    if module == "identity":                 # the resolver is online: a canned answer (the page never calls common.structure)
        from test_api import _fake_resolved
        monkeypatch.setattr(resolve, "resolve", lambda q, kind="Auto", **kw: _fake_resolved(query=q))
        monkeypatch.setattr(enrich, "enrich", lambda r: {"outputs": enrich.empty_outputs(), "warnings": []})
    if module == "ms2_esi":                  # CFM-ID runs in Docker: a recorded cfm-predict log; the precursor m/z and ion formula are still calculated
        monkeypatch.setattr(cf, "check_docker", lambda: None)
        monkeypatch.setattr(cf, "_msrb", lambda tmp, smi, adduct: ((Path(__file__).parent / "fixtures" / "cfmid_cfm_caffeine.log").read_text(encoding="utf-8"), ""))
        cf._predict.cache_clear()
    if module == "ms2_fragments" and not HAS_JAVA:
        pytest.skip("MetFrag needs Java")
    inputs = copy.deepcopy(next(m for m in spec.MODULES if m["id"] == module)["example"])
    inputs["structure"] = {"mode": "SMILES / InChI", "value": smiles}
    if module == "exclusivity":              # the isomer check needs MS-FINDER; the formula and theoretical m/z do not
        inputs.pop("ms2")
    response = client.post(f"/api/modules/{module}/run", json=inputs)
    if response.status_code == 502:          # an engine or data file that is not installed here
        pytest.skip(response.json()["detail"][:120])
    assert response.status_code == 200, response.text
    return response.json()


def isotope_notes(body):
    return [w for w in body["warnings"] if "isotope labels" in w]


def test_every_open_structure_page_has_a_decision():
    assert sorted(STRUCTURE_PAGES) == sorted(LABELLED)


@pytest.mark.parametrize("module", STRUCTURE_PAGES)
def test_a_page_says_so_when_its_outputs_ignore_isotope_labels(module, monkeypatch):
    smiles, expected = LABELLED[module]
    assert "[" in smiles and common.has_isotope_labels(smiles) and not common.has_isotope_labels(unlabelled(smiles))
    notes = isotope_notes(run_page(module, smiles, monkeypatch))
    assert len(notes) == expected, notes
    assert len(set(notes)) == len(notes)                                            # never the same note twice
    assert not isotope_notes(run_page(module, unlabelled(smiles), monkeypatch))


def test_the_notes_name_what_the_page_ignores_the_labels_in(monkeypatch):
    said = {m: isotope_notes(run_page(m, LABELLED[m][0], monkeypatch)) for m in ("adducts", "ms1_theory", "derivatization", "platform_fit", "exclusivity", "ms2_fragments", "ms2_esi")
            if m != "ms2_fragments" or HAS_JAVA}
    assert "ion formulas and m/z values are calculated for the unlabelled compound" in said["adducts"][0]
    assert "ion formula, m/z values and isotope pattern are" in said["ms1_theory"][0]
    assert "derivative formulas and masses are" in said["derivatization"][0] and said["derivatization"] == said["platform_fit"]       # platform_fit's MW is labelled
    assert "same-formula search and the theoretical precursor m/z" in said["exclusivity"][0]
    assert "precursor m/z and ion formula are" in said["ms2_esi"][0]
    assert all("Calculator › Isotope labelling" in n for notes in said.values() for n in notes)


@pytest.mark.parametrize("module", ["platform_fit", "derivatization"])
def test_no_derivative_no_note(module, monkeypatch):
    """d3-caffeine has no group any reagent reacts with, so no derivative mass is shown and nothing ignores the labels."""
    assert not isotope_notes(run_page(module, CAFFEINE_D3, monkeypatch))


def test_the_pages_that_say_so_really_ignore_the_labels(monkeypatch):
    def shown(body, module):
        outputs = body["outputs"]
        return {"adducts": lambda: outputs["ions"]["rows"], "ms1_theory": lambda: outputs["pattern"]["rows"],
                "derivatization": lambda: outputs["deriv_rows"]["rows"], "ms2_esi": lambda: outputs["summary"]["values"]["Precursor m/z"]}[module]()

    for module in ("adducts", "ms1_theory", "derivatization", "ms2_esi"):
        smiles = LABELLED[module][0]
        assert shown(run_page(module, smiles, monkeypatch), module) == shown(run_page(module, unlabelled(smiles), monkeypatch), module), module


def test_the_pages_that_stay_silent_follow_the_labels(monkeypatch):
    mass = lambda body: body["outputs"]["summary"]["values"]["Monoisotopic mass"]
    assert mass(run_page("properties", GLUCOSE_13C6, monkeypatch)) == "186.0835" and mass(run_page("properties", GLUCOSE, monkeypatch)) == "180.0634"


def test_the_ei_spectrum_follows_the_labels(monkeypatch):
    ion = lambda smiles: run_page("ms2_ei", smiles, monkeypatch)["outputs"]["summary"]["values"]["Molecular ion m/z"]["value"]
    assert ion(CAFFEINE_D3) == "197" and ion(unlabelled(CAFFEINE_D3)) == "194"


# ---------------------------------------------------------------- the helper
def test_note_isotope_labels_adds_one_note_for_any_label():
    for smiles in (GLUCOSE_13C6, CAFFEINE_D3, "[13CH4]", "C[15NH2]"):
        warnings = []
        common.note_isotope_labels(warnings, smiles, "The ion formulas are")
        common.note_isotope_labels(warnings, smiles, "The ion formulas are")         # a page that checks twice: still one note
        assert warnings == [common.isotope_note("The ion formulas are")], smiles
    assert common.isotope_note("The ion formulas are") == ("The structure has isotope labels. The ion formulas are calculated for the unlabelled compound; "
                                                           "use Calculator › Isotope labelling for labelled masses.")


def test_note_isotope_labels_is_silent_without_labels():
    for smiles in (GLUCOSE, "CCO", "[Na+].[Cl-]", "[H][H]", "c1ccccc1", "", "not a smiles"):
        warnings = []
        common.note_isotope_labels(warnings, smiles, "The ion formulas are")
        assert warnings == [], smiles


def test_structure_keeps_the_labels_and_adds_no_note_of_its_own():
    """The page decides (note_isotope_labels); structure() only reads the field. A label survives an InChI."""
    inchi = Chem.MolToInchi(Chem.MolFromSmiles(GLUCOSE_13C6))
    assert "/i" in inchi
    for value in (GLUCOSE_13C6, CAFFEINE_D3, inchi):
        warnings = []
        got = common.structure({"structure": {"mode": "SMILES / InChI", "value": value}}, warnings=warnings)
        assert warnings == [] and common.has_isotope_labels(got["smiles"]), value


def test_targeted_method_says_it_in_its_own_words_and_only_once():
    # The route-specific note of Targeted Method (LC route only: the EI M+. of the GC route is calculated from the labelled mass) comes from the workflow,
    # not from the form reader. The full run with the count is in tests/test_workflow_targeted_method.py.
    q = tm.parse_inputs({"structure": {"mode": "SMILES / InChI", "value": CAFFEINE_D3}, "ph": "", "additive": "Formic acid", "instrument": "Any"})
    assert q.warnings == []
