"""Shared structure preparation (mzdesk/platform_fit/structure_prep.py): the pKa / Ionization / Adducts / Separation pages and the Platform /
Derivatization pages choose the same fragment, neutralize the same way and give the same warning; Targeted Method does not repeat it."""
import pytest
from rdkit import Chem

from mzdesk.identity import resolve
from mzdesk.platform_fit import derivatization as d
from mzdesk.platform_fit import pka_engine as pe
from mzdesk.platform_fit import structure_prep as sp
from mzdesk.workflows import targeted_method as tm

PARTS_NOTE = "2 disconnected parts (salt or mixture): only the largest is analysed."


def test_the_pka_engine_uses_the_shared_helpers():
    assert pe._neutralize_anions is sp.neutralize_anions and pe._largest_fragment is sp.largest_fragment
    assert sp.parts_note(2) == PARTS_NOTE and sp.parts_note(3).startswith("3 disconnected parts")


def test_largest_fragment_prefers_the_organic_part():
    mol, n = sp.largest_fragment(Chem.MolFromSmiles("[Na+].[O-]C(=O)c1ccccc1.[Cl-]"))
    assert n == 3 and Chem.MolToSmiles(mol) == "O=C([O-])c1ccccc1"
    mol, n = sp.largest_fragment(Chem.MolFromSmiles("CCO"))
    assert n == 1 and Chem.MolToSmiles(mol) == "CCO"


def test_neutralize_anions_leaves_charge_separated_groups_alone():
    assert Chem.MolToSmiles(sp.neutralize_anions(Chem.MolFromSmiles("C[N+](C)(C)CC(=O)[O-]"))[0]) == "C[N+](C)(C)CC(=O)O"
    assert sp.neutralize_anions(Chem.MolFromSmiles("c1ccccc1[N+](=O)[O-]"))[1] is False
    assert sp.neutralize_anions(Chem.MolFromSmiles("[O-][n+]1ccccc1"))[1] is False           # pyridine N-oxide: the anion sits on the cation


@pytest.mark.parametrize("smiles", ["NCCc1ccc(O)c(O)c1.Cl", "[Na+].[O-]C(=O)c1ccccc1", "C[N+](C)(C)CC(=O)[O-]"])
def test_prepare_is_what_the_other_pages_use(smiles):
    """Same molecule as pka_engine.prepare (atom order and all), except that pKa's own heavy-atom-only form drops the explicit hydrogens."""
    mol, n_parts, _ = d.prepare(Chem.MolFromSmiles(smiles))
    prep = pe.prepare(smiles)
    assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(prep.mol) and n_parts == len(Chem.GetMolFrags(Chem.MolFromSmiles(smiles)))


def test_targeted_method_gives_a_salt_warning_once(monkeypatch):
    def offline(query, kind="Auto", **kw):
        raise resolve.ResolveError("Network lookup failed (PubChem/UniChem unreachable): test")
    monkeypatch.setattr(resolve, "resolve", offline)
    form = {"structure": {"mode": "SMILES / InChI", "value": "NCCc1ccc(O)c(O)c1.Cl"}, "ph": "", "additive": "Formic acid", "instrument": "Any"}
    result = tm.run_query(tm.parse_inputs(form))
    assert [w for w in result["warnings"] if "disconnected parts" in w] == [PARTS_NOTE]
    platform = next(s for s in result["steps"] if s["id"] == 2)
    assert platform["status"] == "done" and "GC-MS after derivatization" in platform["note"]


def test_targeted_method_platform_step_for_a_quaternary_acid_and_a_phosphate_ester():
    from mzdesk.workflows import targeted_method as t

    class Q:
        pass
    for smiles, gc, cil in (("C[N+](C)(C)CC(=O)[O-]", "no", True), ("CCOP(=O)(OCC)OCC", "yes", False), ("COP(=O)(OC)OC=C(Cl)Cl", "yes", False)):
        q = Q()
        q.smiles = smiles
        _, _, _, value = t.step_platform(q)
        assert value["gc"]["gc"] == gc and value["cil"]["cil"] is cil, smiles
