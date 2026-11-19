"""Lipid Short Name 页：拷进来的 lipid_name_map（版本、和源包逐项一致、数据文件、没有绝对导入）、runner 输出、示意结构（分子式、手性、高亮、SVG）、错误输入。

和源包的对照用固定的 golden 文件（tests/fixtures/lipid/name_map_golden.json，由源包 0.4.0 生成，421 个名字），不依赖仓库外的路径。
示意结构的 golden（illustration_golden.json）是当前代码的输出，防止画法被无意改动；改画法时重新生成。
分子式有独立的参照：reference_formulas.json 是 LIPID+ 参考表里已知结构的分子式（按 L1 名），画出的结构必须同分子式；画出来之前还有一道按"类常数 + 总量"的核对。
"""
import json
import re
import time
import xml.dom.minidom
from pathlib import Path

import pytest
from fastapi.testclient import TestClient
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from mzdesk import spec
from mzdesk.api import verdicts
from mzdesk.api.app import app
from mzdesk.api.common import InputError
from mzdesk.api.runners import RUNNERS
from mzdesk.api.runners.lipid import run_lipid_short_name
from mzdesk.lipid import draw, structure as ls, texts
from mzdesk.lipid_name_map import __version__ as NAME_MAP_VERSION, contract, lm, normalize, reasons

HERE = Path(__file__).parent
GOLDEN = json.loads((HERE / "fixtures/lipid/name_map_golden.json").read_text(encoding="utf-8"))
ILLUSTRATIONS = json.loads((HERE / "fixtures/lipid/illustration_golden.json").read_text(encoding="utf-8"))
ILLUSTRATIONS_ASSUMPTIONS = [a for v in ILLUSTRATIONS.values() for a in v.get("assumptions", [])]
REFERENCE_FORMULAS = json.loads((HERE / "fixtures/lipid/reference_formulas.json").read_text(encoding="utf-8"))
PACKAGE = Path(ls.__file__).resolve().parents[1] / "lipid_name_map"
client = TestClient(app)


def canon(smiles):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def run(name):
    return run_lipid_short_name({"name": name})["outputs"]


def cip(mol):
    m = Chem.Mol(mol)
    rdCIPLabeler.AssignCIPLabels(m)
    return [a.GetProp("_CIPCode") for a in m.GetAtoms() if a.HasProp("_CIPCode")], [b.GetProp("_CIPCode") for b in m.GetBonds() if b.HasProp("_CIPCode")]


# ---------------------------------------------------------------- the copied package
def test_vendored_package_version_and_location():
    assert NAME_MAP_VERSION == "0.4.0"
    assert (PACKAGE / "__init__.py").is_file() and (PACKAGE / "VENDORED.md").is_file()
    assert "mzdesk" in str(PACKAGE) and PACKAGE.name == "lipid_name_map"


def test_vendored_package_matches_the_source_package_on_golden_names():
    assert GOLDEN["version"] == "0.4.0" and len(GOLDEN["results"]) >= 200
    wrong = []
    for want in GOLDEN["results"]:
        got = normalize(want["input"])
        for key, value in want.items():
            if got[key] != value:
                wrong.append((want["input"], key, got[key], value))
    assert not wrong, wrong[:10]


def test_golden_covers_every_status_and_level():
    statuses = {r["status"] for r in GOLDEN["results"]}
    levels = {r["level"] for r in GOLDEN["results"]}
    assert statuses == {"ok", "rejected", "unsupported"} and levels == {0, 1, 2, 3}


def test_vendored_package_has_no_absolute_self_import_and_finds_its_data():
    """Relative imports only (a vendored copy must not need a top-level lipid_name_map), and every data file is inside the package."""
    absolute = re.compile(r"^\s*(?:from|import)\s+lipid_name_map\b", re.M)
    offenders = [str(p.relative_to(PACKAGE)) for p in PACKAGE.rglob("*.py") if absolute.search(p.read_text(encoding="utf-8"))]
    assert not offenders, offenders
    data = [p for p in (PACKAGE / "data").rglob("*") if p.is_file()]
    assert len(data) >= 19 and {p.suffix for p in data} == {".csv", ".json"}
    for module_path in (lm.CSV_PATH, contract.VOCAB_PATH):
        assert PACKAGE in module_path.parents and module_path.is_file()


def test_package_data_pattern_covers_the_data_files():
    text = (HERE.parent / "pyproject.toml").read_text(encoding="utf-8")
    assert '"lipid_name_map/data/**/*"' in text


# ---------------------------------------------------------------- names and texts
def test_every_reason_code_and_flag_has_a_sentence():
    for name in ("PC 16:0/18:1", "SM 34:1;O2", "Cer 18:1;O2/24:0", "TG 52:3", "FA 20:4;OH", "CE 18:1", "CoQ9", "ST 24:1;O4;Phe"):
        normalize(name)                                  # loads every family, which registers its flags
    assert set(reasons.REASONS) <= set(texts.REASON_TEXT), set(reasons.REASONS) - set(texts.REASON_TEXT)
    assert set(contract.all_flags()) <= set(texts.FLAG_TEXT), set(contract.all_flags()) - set(texts.FLAG_TEXT)
    assert set(lm.LM_REASONS) <= set(texts.LM_REASON_TEXT), set(lm.LM_REASONS) - set(texts.LM_REASON_TEXT)
    assert all(t.endswith(".") and len(t) < 110 for t in [*texts.REASON_TEXT.values(), *texts.FLAG_TEXT.values(), *texts.LM_REASON_TEXT.values()])


@pytest.mark.parametrize("name, plus, maps", [
    ("PC 16:0/18:1", {"L1": "PC 34:1", "L2": "PC 18:1_16:0", "L3": "PC 18:1_16:0"}, {"Species": "PC 34:1", "Molecular species": "PC 16:0_18:1"}),
    ("PC 34:1", {"L1": "PC 34:1"}, {"Species": "PC 34:1"}),
    ("LPC 18:1", {"L1": "PC-lyso 18:1", "L2": "PC-lyso 18:1", "L3": "PC-lyso 18:1"}, {"Species": "LPC 18:1", "Molecular species": "LPC 18:1"}),
    ("PE P-18:0/20:4", {"L1": "PE-ether 38:5", "L2": "PE-ether 20:4_18:1", "L3": "PE-ether 20:4_O-18:1;P"}, {"Species": "PE O-38:5", "Molecular species": "PE P-18:0_20:4"}),
    ("Cer(d18:1/24:0)", {"L1": "Cer 42:1;O2", "L2": "Cer 24:0_18:1;(O2)", "L3": "Cer 24:0_18:1;OH2"}, {"Species": "Cer 42:1;O2", "Molecular species": "Cer 18:1;O2/24:0"}),
    ("Cer 24:0_18:1;OH2", {"L1": "Cer 42:1;O2", "L2": "Cer 24:0_18:1;(O2)", "L3": "Cer 24:0_18:1;OH2"}, {"Species": "Cer 42:1;O2", "Molecular species": "Cer 18:1;O2/24:0"}),
    ("SM 34:1;O2", {"L1": "SM 34:1;O2"}, {"Species": "SM 34:1;O2"}),
    ("TG 52:3", {"L1": "TG 52:3"}, {"Species": "TG 52:3"}),
    ("CE 18:1", {"L1": "ST-oacyl 45:2;O", "L2": "ST-oacyl 27:1_18:1;(O)", "L3": "ST-oacyl 27:1;OH_18:1"}, {"Species": "CE 18:1", "Molecular species": "CE 18:1"}),
    ("SE 27:1/18:1", {"L1": "ST-oacyl 45:2;O", "L2": "ST-oacyl 27:1_18:1;(O)", "L3": "ST-oacyl 27:1;OH_18:1"}, {"Species": "SE 27:1/18:1", "Molecular species": "SE 27:1/18:1"}),
    ("FA 20:4;OH", {"L1": "FA 20:4;O", "L2": "FA 20:4;(O)", "L3": "FA 20:4;OH"}, {"Species": "FA 20:4;O", "Molecular species": "FA 20:4;O"}),
    ("ST 27:1;O", {"L1": "ST 27:1;O", "L2": "ST 27:1;(O)"}, {"Species": "ST 27:1;O"}),
])
def test_runner_names(name, plus, maps):
    o = run(name)
    assert o["lipid_plus"]["values"] == plus
    got = {k: (v["value"] if isinstance(v, dict) else v) for k, v in o["lipid_maps"]["values"].items()}
    assert {k: v for k, v in got.items() if v != "—"} == maps
    for label in ("Species", "Molecular species"):          # a blank LIPID MAPS name says why in its "?"
        if label not in maps:
            assert o["lipid_maps"]["values"][label]["value"] == "—" and o["lipid_maps"]["values"][label]["help"]


def test_runner_blocks_match_the_spec():
    m = next(m for m in spec.MODULES if m["id"] == "lipid_short_name")
    ids = [b["id"] for b in m["outputs"]]
    assert ids == ["lipid_plus", "lipid_maps", "structure", "no_structure", "assumptions", "flags"]
    assert m["live"] and m["subgroup"] == "molecular" and m["example"] == {"name": "PC 16:0/18:1"} and not m.get("network")
    assert [f["id"] for f in m["inputs"]] == ["name"] and m["inputs"][0]["type"] == "text" and m["inputs"][0]["required"]
    assert all(b.get("style") == "metric" for b in m["outputs"] if b["type"] == "cards")
    assert "lipid_short_name" in RUNNERS
    for name in ("PC 16:0/18:1", "ST 27:1;O"):
        assert set(run(name)) >= {b for b in ids if b != "structure"}


def test_names_without_a_drawing_say_so_in_one_line():
    for name in ("ST 27:1;O", "ST 24:1;O4;Phe", "BMP 18:1_18:1", "TG-oacyl 16:0;OH_16:0_16:0_16:0"):
        o = run(name)
        assert "structure" not in o, name
        assert o["no_structure"]["items"] == ["No illustrative structure for this lipid class."], name
        assert o["assumptions"]["items"] == []


def test_flags_are_listed_with_a_short_meaning():
    """English sentences only (no internal code), and nothing the "?" of a blank LIPID MAPS card already says."""
    assert run("PC 16:0/18:1")["flags"]["items"] == []
    assert run("PC 34:1")["flags"]["items"] == []                     # sum_only: the "—" of Molecular species explains it
    assert run("PC 34:1")["lipid_maps"]["values"]["Molecular species"]["help"] == "The name gives only the sum composition."
    assert run("PC 16:0/18:1 [M+H]+")["flags"]["items"] == ["Adduct removed from the name."]
    ce = run("CE 18:1")["flags"]["items"]
    assert "Class read from an alias." in ce and len(ce) == 4
    assert run("FAHFA 36:1;O2")["flags"]["items"] == ["Only the sum composition is given."]       # its "—" says something else (O-acyl class), so the flag stays


def test_flags_never_show_an_internal_code():
    sentences = set(texts.FLAG_TEXT.values())
    for want in GOLDEN["results"]:
        if want["status"] != "ok":
            continue
        for line in run_lipid_short_name({"name": want["input"]})["outputs"]["flags"]["items"]:
            assert line in sentences and not re.match(r"^[a-zA-Z_]+: ", line), (want["input"], line)


def test_a_flag_said_by_a_blank_card_is_not_repeated():
    r = normalize("PE P-36:3")
    assert r["meta"]["lm"]["molecular"] == "sum_only" and "sum_only" in r["flags"]
    assert "Only the sum composition is given." not in run("PE P-36:3")["flags"]["items"]
    assert run("PE P-36:3")["flags"]["items"] == ["P- chain: the vinyl C=C counts in the double bonds."]


def test_the_example_adduct_and_spacing_are_accepted():
    assert run("  PC 16:0/18:1  ")["lipid_plus"]["values"]["L1"] == "PC 34:1"
    assert run("PC(16:0/18:1)+H")["lipid_maps"]["values"]["Molecular species"] == "PC 16:0_18:1"


COMMON_NAMES = "Common names are not supported. Use a shorthand such as PC 16:0/18:1."


@pytest.mark.parametrize("name, fragment", [
    ("Phosphatidylcholine", COMMON_NAMES),
    ("palmitic acid", COMMON_NAMES),
    ("POPC", COMMON_NAMES),
    ("PGE2", COMMON_NAMES),
    ("oleic acid", COMMON_NAMES),
    ("cholesterol", COMMON_NAMES),
    ("PC", COMMON_NAMES),
    ("α-linolenic acid", COMMON_NAMES),
    ("palmitic acid (16:0)", "Common and long class names are not accepted"),      # has a chain composition: the package's own sentence
    ("oleic acid, 18:1", "Common and long class names are not accepted"),
    ("", "Enter a lipid name."),
    ("   ", "Enter a lipid name."),
    ("12345", "database ID or a number"),
    ("16:0_16:0", "Chains only"),
    ("lysoPC 18:1", "Unknown class name."),
    ("PC-d7 16:0/18:1", "Isotope-labelled"),
    ("PC 16:0/18:1/20:4/22:6", "number of chains"),
])
def test_rejected_names_give_a_readable_reason(name, fragment):
    with pytest.raises(InputError) as exc:
        run_lipid_short_name({"name": name})
    assert fragment in str(exc.value)
    r = client.post("/api/modules/lipid_short_name/run", json={"name": name})
    assert r.status_code == 400 and fragment in r.json()["detail"]


def test_missing_or_odd_input_is_a_400_not_a_500():
    for payload in ({}, {"name": None}, {"name": 12345}, {"name": ["PC 34:1"]}, {"name": "x" * 5000}, {"name": "PC 34:1\u0000"}):
        r = client.post("/api/modules/lipid_short_name/run", json=payload)
        assert r.status_code in (200, 400), (payload, r.status_code, r.text[:200])


def test_every_golden_name_runs_without_a_server_error():
    """The golden names include rejected, unsupported and every level: the page answers each with 200 or a 400 reason, never a 500."""
    drawn = 0
    for want in GOLDEN["results"]:
        if want["status"] != "ok":
            with pytest.raises(InputError):
                run_lipid_short_name({"name": want["input"]})
            continue
        try:                                               # the drawing code itself must only ever refuse with Unsupported
            ls.illustrate(want["input"])
            drawn += 1
        except ls.Unsupported:
            pass
    assert drawn > 60


def test_the_page_answers_a_spread_of_golden_names():
    ok = [w for w in GOLDEN["results"] if w["status"] == "ok"][::7]
    assert len(ok) >= 40
    for want in ok:
        out = run_lipid_short_name({"name": want["input"]})["outputs"]
        assert out["lipid_plus"]["values"]["L1"] == want["L1"] and ("structure" in out) != bool(out["no_structure"]["items"]), want["input"]


def test_example_through_the_api_gives_a_verdict_and_a_drawing():
    r = client.post("/api/modules/lipid_short_name/run", json=next(m for m in spec.MODULES if m["id"] == "lipid_short_name")["example"])
    body = r.json()
    assert r.status_code == 200 and body["status"] == "ok"
    assert body["verdict"] == "<b>PC 18:1_16:0</b> in LIPID+, <b>PC 16:0_18:1</b> in LIPID MAPS."
    assert body["outputs"]["structure"]["svg"].lstrip().startswith("<?xml") and "<svg" in body["outputs"]["structure"]["svg"]


def test_verdict_wording():
    assert verdicts.build("lipid_short_name", {}, run("PC 34:1")) == "<b>PC 34:1</b> in LIPID+ and LIPID MAPS."
    assert verdicts.build("lipid_short_name", {}, run("PC 16:0/18:1")) == "<b>PC 18:1_16:0</b> in LIPID+, <b>PC 16:0_18:1</b> in LIPID MAPS."
    assert verdicts.build("lipid_short_name", {}, {}) is None


# ---------------------------------------------------------------- illustrative structures: formulas
KNOWN_FORMULAS = [
    ("PC 16:0/18:1", "C42H82NO8P"),                  # POPC
    ("PE 16:0/18:1", "C39H76NO8P"),
    ("PS 16:0/18:1", "C40H76NO10P"),
    ("PG 16:0/18:1", "C40H77O10P"),
    ("PI 16:0/18:1", "C43H81O13P"),
    ("PA 16:0/18:1", "C37H71O8P"),
    ("LPC 18:1", "C26H52NO7P"),
    ("LPE 18:1", "C23H46NO7P"),
    ("PE P-18:0/20:4", "C43H78NO7P"),
    ("PC O-16:0/18:1", "C42H84NO7P"),
    ("DG 16:0/18:1", "C37H70O5"),
    ("MG 18:1", "C21H40O4"),
    ("TG 52:3", "C55H100O6"),
    ("TG 16:0_18:1_18:2", "C55H100O6"),
    ("CL 18:2_18:2_18:2_18:2", "C81H142O17P2"),
    ("Cer(d18:1/24:0)", "C42H83NO3"),
    ("Cer 42:1;O2", "C42H83NO3"),
    ("SM d18:1/16:0", "C39H79N2O6P"),
    ("SM 34:1;O2", "C39H79N2O6P"),
    ("HexCer d18:1/16:0", "C40H77NO8"),
    ("SPB 18:1;O2", "C18H37NO2"),
    ("SPB 18:0;O3", "C18H39NO3"),
    ("CE 18:1", "C45H78O2"),
    ("FA 20:4", "C20H32O2"),
    ("FA 20:4;OH", "C20H32O3"),
    ("FOH 18:1", "C18H36O"),
    ("CAR 16:0", "C23H45NO4"),
    ("FAHFA 16:0_18:0;OH", "C34H66O4"),
]


@pytest.mark.parametrize("name, formula", KNOWN_FORMULAS)
def test_illustration_formula_matches_the_known_molecule(name, formula):
    ill = ls.illustrate(name)
    assert ill.formula == formula
    assert CalcMolFormula(Chem.MolFromSmiles(ill.smiles)) == formula


@pytest.mark.parametrize("sum_name, chain_name", [
    ("PC 34:1", "PC 16:0/18:1"), ("PE 38:4", "PE 18:0/20:4"), ("TG 52:3", "TG 16:0/18:1/18:2"), ("Cer 42:1;O2", "Cer d18:1/24:0"), ("SM 34:1;O2", "SM d18:1/16:0"),
    ("PC-lyso 18:1", "PC-lyso 18:1"), ("DG 34:1", "DG 16:0/18:1"), ("PE-ether 38:5", "PE P-18:0/20:4"), ("CE 18:1", "SE 27:1/18:1"), ("CL 72:8", "CL 18:2_18:2_18:2_18:2"),
])
def test_a_sum_name_and_its_chain_name_draw_the_same_formula(sum_name, chain_name):
    assert ls.illustrate(sum_name).formula == ls.illustrate(chain_name).formula


def test_drawn_formula_equals_the_known_structure_of_the_same_l1_name():
    """Independent of structure.py: the formula of a LIPID+ reference structure with this L1 name (446 names, all drawable classes, with and without ;O)."""
    drawn, wrong = 0, []
    for name, formula in REFERENCE_FORMULAS.items():
        try:
            ill = ls.illustrate(name)
        except ls.Unsupported:
            continue
        drawn += 1
        if ill.formula != formula:
            wrong.append((name, ill.formula, formula))
    assert drawn > 400 and not wrong, wrong[:10]


def test_expected_elements_are_class_constant_plus_the_sum():
    r = normalize("PC 34:1")
    assert ls.expected_elements(r) == dict(C=42, H=82, N=1, O=8, P=1, S=0)
    assert ls.expected_elements(normalize("PE P-18:0/20:4")) == dict(C=43, H=78, N=1, O=7, P=1, S=0)          # ether: H +2, O -1; the P- vinyl C=C is in the DB
    assert ls.expected_elements(normalize("CL 72:8;O4")) == dict(C=81, H=142, N=0, O=21, P=2, S=0)           # each ;O is one oxygen
    assert ls.expected_elements(normalize("ST 27:1;O")) is None                                                # a class without a constant
    for name in ("Cer 42:1;O2", "SM 34:1;O2", "TG 52:3", "FAHFA 34:0;O", "CE 18:1", "CAR 16:0", "HexCer 34:1;O2", "FA 20:4;OH"):
        r = normalize(name)
        assert ls.expected_elements(r) == ls.drawn_elements(ls.illustrate(r).mol), name


def test_a_drawing_that_disagrees_with_the_name_is_refused(monkeypatch):
    """The gate itself: whatever build() returns, the atoms have to add up to the name."""
    real = ls.build
    monkeypatch.setattr(ls, "build", lambda r, raw=None: ("CCCCCCCC(=O)OC", [], False))
    with pytest.raises(ls.Unsupported, match="does not match"):
        ls.illustrate("PC 16:0/18:1")
    monkeypatch.setattr(ls, "build", real)
    assert ls.illustrate("PC 16:0/18:1").formula == "C42H82NO8P"


@pytest.mark.parametrize("name", [
    "Cer 16:0_18:1;O",        # one untyped O: the sphingoid base was drawn as d type (two OH), one oxygen too many
    "DG O-16:1", "PC 4:1", "FA 2:1", "PA 2:1", "TG 3:2",      # the chain is too short for its C=C: the C=C used to be dropped silently
])
def test_a_name_that_cannot_be_drawn_to_its_formula_gets_no_drawing(name):
    with pytest.raises(ls.Unsupported):
        ls.illustrate(name)
    o = run(name)
    assert "structure" not in o and o["no_structure"]["items"] == ["No illustrative structure for this name."] and o["assumptions"]["items"] == []
    assert o["lipid_plus"]["values"]["L1"]                                 # the names are still the answer


@pytest.mark.parametrize("name, l1", [
    ("FA 17:0(15Me)", "FA 17:0"),            # was drawn as the straight C17 acid (C17H34O2); the real molecule has one more carbon
    ("FA iso-17:0", "FA 17:0"),
    ("FA anteiso-17:0", "FA 17:0"),
    ("FA 16:0(2Me)", "FA 16:0"),
    ("PC 16:0/17:0(15Me)", "PC 33:0"),
])
def test_a_methyl_branched_chain_gets_no_drawing_and_says_why(name, l1):
    """A straight chain would be a different molecule, so the page names the lipid but does not draw it. The names are unchanged."""
    with pytest.raises(ls.BranchedChain):
        ls.illustrate(name)
    r = client.post("/api/modules/lipid_short_name/run", json={"name": name})
    assert r.status_code == 200
    o = r.json()["outputs"]
    assert "structure" not in o and o["no_structure"]["items"] == ["Branched chain: structure not drawn."]
    assert o["assumptions"]["items"] == []
    assert o["lipid_plus"]["values"]["L1"] == l1 and run(name) == o


METHYL_FLAG = "Methyl branch carbon not counted at the sum level."
ISO_FLAG = "Branch type (iso/anteiso) not kept at the sum level."


@pytest.mark.parametrize("name, lines", [
    ("FA 17:0(15Me)", [METHYL_FLAG]),                  # positional methyl: the CH3 is not in L1
    ("FA 16:0(2Me)", [METHYL_FLAG]),
    ("FA 17:0-15me", [METHYL_FLAG]),
    ("FA iso-17:0", [ISO_FLAG]),                       # iso / anteiso: the carbon count already holds the branch carbon, the branch type is what is lost
    ("FA anteiso-17:0", [ISO_FLAG]),
    ("PC 16:0/iso-17:0", [ISO_FLAG]),
    ("PC iso-16:0/17:0(15Me)", [METHYL_FLAG, ISO_FLAG]),     # one of each, on two chains
    ("PC 16:0/17:0(15Me) [M+H]+", ["Adduct removed from the name.", METHYL_FLAG]),
])
def test_the_branch_flag_says_what_the_name_lost(name, lines):
    assert run(name)["flags"]["items"] == lines


@pytest.mark.parametrize("name", ["PC 34:1", "PC 16:0/18:1", "TG 52:2", "PE P-18:0/20:4", "Cer 18:1;O2/16:0"])
def test_common_unbranched_lipids_are_still_drawn(name):
    o = run(name)
    assert "structure" in o and o["no_structure"]["items"] == []
    assert ls.illustrate(name).formula


@pytest.mark.parametrize("name, formula", [
    ("CL 72:8;O4", "C81H142O21P2"),         # used to lose all 4 oxygens
    ("CL 18:2_18:2_18:2_18:2;O", "C81H142O18P2"),
    ("CE 18:2;O", "C45H76O3"),
    ("FAHFA 36:1;O2", "C36H68O5"),
    ("FAHFA 16:0_18:0;O2", "C34H66O5"),
    ("TG 52:3;O", "C55H100O7"),
    ("FA 20:4;O2", "C20H32O4"),
])
def test_untyped_oxygen_is_drawn_as_oh_and_the_formula_still_matches(name, formula):
    ill = ls.illustrate(name)
    assert ill.formula == formula
    assert "Untyped O drawn as OH" in ill.assumptions
    assert ill.mol.HasSubstructMatch(Chem.MolFromSmarts("[CX4][OX2H1]"))
    assert "Untyped O drawn as OH" in run(name)["assumptions"]["items"]


def test_sterol_ester_with_two_sterol_like_chains_takes_the_one_with_the_oh():
    """CE 28:2 = 28:2 acyl + cholesterol (the class implies cholesterol), not a campesterol-type 28:2 sterol with a 27:1 acyl chain."""
    ill = ls.illustrate("CE 28:2")
    assert ill.formula == "C55H96O2"
    cholesteryl = ls.illustrate("SE 27:1/28:2")
    assert ill.smiles == cholesteryl.smiles
    assert ill.mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)O[C@H]1CC[C@@]2(C)C(=CC[C@H]3[C@@H]4CC[C@H]([C@H](C)CCCC(C)C)[C@@]4(C)CC[C@@H]32)C1"))
    r = normalize("CE 28:2")
    assert [(c["C"], c["DB"], c["toks"]) for c in r["meta"]["chains"]] == [(28, 2, {}), (27, 1, {"OH": 1})]
    # both chains look like sterols: the OH marks the sterol; with no OH, 27:1 (cholesterol) wins over the others
    assert ls.illustrate("SE 28:2/28:2").formula == "C56H96O2"
    assert ls.illustrate("CE 29:1").formula == "C56H100O2" and ls.illustrate("CE 29:1").smiles == ls.illustrate("SE 27:1/29:1").smiles
    assert ls.illustrate("CE 27:1").formula == "C54H96O2"
    assert ls.illustrate("SE 18:1/27:1").smiles == ls.illustrate("SE 27:1/18:1").smiles
    assert ls.illustrate("SE 29:1/28:1").formula == "C57H102O2"          # the 29:1 carries the OH: it is the sterol


def test_illustrations_match_the_stored_golden():
    for name, want in ILLUSTRATIONS.items():
        try:
            got = ls.illustrate(name)
        except ls.Unsupported as exc:
            assert want == {"error": type(exc).__name__}, name
            continue
        assert (got.smiles, got.formula, got.assumptions) == (want["smiles"], want["formula"], want["assumptions"]), name
        assert (len(got.db_atoms), len(got.fg_atoms), len(got.split_atoms)) == (want["db"], want["fg"], want["split"]), name


# ---------------------------------------------------------------- illustrative structures: stereochemistry
def test_popc_and_dopg_equal_the_pubchem_structures_used_elsewhere_in_mzdesk():
    assert ls.illustrate("PC 16:0/18:1").smiles == canon(spec.EXAMPLE_STRUCTURES["POPC"])
    assert ls.illustrate("PG 18:1_18:1").smiles == canon(spec.EXAMPLE_STRUCTURES["DOPG"])


def test_cholesteryl_oleate_equals_cholesterol_esterified():
    """The sterol skeleton is checked centre by centre against the cholesterol SMILES the app already uses (PubChem)."""
    oleate = spec.CHOLESTEROL.replace("(C4)O)", "(C4)OC(=O)CCCCCCC/C=C\\CCCCCCCC)")
    assert oleate != spec.CHOLESTEROL
    assert ls.illustrate("CE 18:1").smiles == canon(oleate) == ls.illustrate("SE 27:1/18:1").smiles
    # the other sterol side chains share the skeleton and the side-chain start; only the carbon and C=C counts differ
    for name, formula in (("SE 28:1/18:1", "C46H80O2"), ("SE 29:1/16:0", "C45H80O2"), ("SE 29:2/18:2", "C47H78O2"), ("SE 27:2/18:1", "C45H76O2")):
        assert ls.illustrate(name).formula == formula


def test_tg_has_the_natural_sn_configuration():
    """Derived from the PubChem POPC SMILES: the same glycerol atom order with an acyl in place of the phosphocholine."""
    popc = spec.EXAMPLE_STRUCTURES["POPC"]
    head = "P(=O)([O-])OCC[N+](C)(C)C"
    assert head in popc
    tg = popc.replace("OP(=O)([O-])OCC[N+](C)(C)C", "OC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC")
    assert tg.count("C(=O)") == 3
    assert ls.illustrate("TG 16:0_18:1_18:2").smiles == canon(tg)


def test_natural_configurations_by_cip():
    assert cip(ls.illustrate("PC 16:0/18:1").mol)[0] == ["R"]                       # glycerol C2 of sn-glycerophospholipids
    assert cip(ls.illustrate("PS 16:0_18:1").mol)[0] == ["R", "S"]                   # + L-serine
    centres, bonds = cip(ls.illustrate("Cer(d18:1/24:0)").mol)
    assert sorted(centres) == ["R", "S"] and bonds == ["E"]                          # (2S,3R,4E)
    sm_centres, sm_bonds = cip(ls.illustrate("SM 34:1;O2").mol)
    assert sorted(sm_centres) == ["R", "S"] and sm_bonds == ["E"]
    assert sorted(cip(ls.illustrate("SPB 18:0;OH3").mol)[0]) == ["R", "S", "S"]       # phytosphingosine (2S,3S,4R)


def test_cardiolipin_outer_glycerols_have_the_same_configuration():
    centres, _ = cip(ls.illustrate("CL 18:2_18:2_18:2_18:2").mol)
    assert centres == ["R", "R"]


def test_c_c_double_bonds_are_cis_except_the_sphingoid_base():
    for name, z, e in (("PC 16:0/18:1", 1, 0), ("FA 20:4", 4, 0), ("Cer(d18:1/24:0)", 0, 1), ("PE P-18:0/20:4", 5, 0)):
        _, bonds = cip(ls.illustrate(name).mol)
        assert (bonds.count("Z"), bonds.count("E")) == (z, e), name


# ---------------------------------------------------------------- illustrative structures: what is flagged as guessed
def atoms_symbols(ill, atoms):
    return sorted({ill.mol.GetAtomWithIdx(i).GetSymbol() for i in atoms})


def test_guessed_parts_are_flagged():
    exact = ls.illustrate("PC 16:0/18:1")                      # chains known, the position of the C=C is not
    assert len(exact.db_atoms) == 2 and exact.fg_atoms == [] and exact.split_atoms == []
    assert exact.assumptions == []                             # the C=C guess is in the key under the drawing; the sn order is the name's
    total = ls.illustrate("PC 34:1")                           # only the sum: the chains are a default split
    assert len(total.split_atoms) == 34 and atoms_symbols(total, total.split_atoms) == ["C"]      # every chain carbon, nothing of the head group
    assert total.assumptions == ["sn positions assumed"]       # the split is in the key too
    hete = ls.illustrate("FA 20:4;OH")
    assert atoms_symbols(hete, hete.fg_atoms) == ["C", "O"] and len(hete.db_atoms) == 8
    assert hete.assumptions == ["Group positions assumed"]
    assert ls.illustrate("PC 16:0_16:0").db_atoms == [] and ls.illustrate("PC 16:0_16:0").assumptions == []


def test_assumptions_do_not_repeat_the_colour_key():
    """"C=C positions assumed" and "Chains split by default rule" are the key labels under the drawing; they are not listed again."""
    key_labels = {label for _, _, label in draw.KEY}
    assert "C=C positions assumed" not in key_labels and not any(a in ("C=C positions assumed", "Chains split by default rule") for a in ILLUSTRATIONS_ASSUMPTIONS)
    for name in ("PC 34:1", "TG 52:3", "CL 72:8", "FA 20:4", "Cer 42:1;O2", "CE 18:1"):
        assert not set(ls.illustrate(name).assumptions) & {"C=C positions assumed", "Chains split by default rule"}, name
    assert run("PC 16:0/18:1")["assumptions"]["items"] == []                     # nothing to say: the block is not shown (web Outputs.tsx hides an empty notes block)


def test_plasmalogen_vinyl_ether_is_not_a_guess():
    ill = ls.illustrate("PE P-18:0/20:4")
    assert len(ill.db_atoms) == 8                              # the four C=C of 20:4 only; the vinyl C=C of P-18:0 is given by the name
    assert ill.mol.HasSubstructMatch(Chem.MolFromSmarts("[CH2]O[CH]=[CH][CH2]"))
    assert "Ether chain drawn as O-alkyl" not in ill.assumptions and "sn positions assumed" not in ill.assumptions


def test_a_sum_level_plasmalogen_is_drawn_as_a_vinyl_ether():
    """PE P-36:3 (ether_linkage P): the O-CH=CH- is given by the name, so it is neither "O-alkyl" nor a guessed C=C."""
    r = normalize("PE P-36:3")
    assert r["level"] == 1 and r["meta"]["ether_linkage"] == "P"
    ill = ls.illustrate("PE P-36:3")
    assert ill.mol.HasSubstructMatch(Chem.MolFromSmarts("[CH2]O[CH]=[CH][CH2]"))
    assert ill.formula == "C41H76NO7P" == ls.illustrate("PE P-16:0/20:3").formula          # C41, 4 C=C in the label: the vinyl one + 3
    assert ill.assumptions == [] and len(ill.db_atoms) == 6                                  # the 3 C=C of the acyl chain only
    assert len(ill.split_atoms) == 36
    for name, formula in (("PC P-34:1", "C42H82NO7P"), ("LPE P-16:0", "C21H44NO6P"), ("PE P-34:2", "C39H74NO7P")):
        d = ls.illustrate(name)
        assert d.formula == formula and d.mol.HasSubstructMatch(Chem.MolFromSmarts("[CH2]O[CH]=[CH]")), name
    assert "Ether chain drawn as O-alkyl" in ls.illustrate("PE-ether 36:3").assumptions       # linkage unknown: still O-alkyl


def test_ether_class_without_linkage_says_it_is_drawn_as_o_alkyl():
    assert "Ether chain drawn as O-alkyl" in ls.illustrate("PE-ether 34:2").assumptions
    assert "Ether chain drawn as O-alkyl" not in ls.illustrate("PC O-16:0/18:1").assumptions


def test_sn_note_only_when_the_chains_differ():
    assert "sn positions assumed" not in ls.illustrate("PC 16:0_16:0").assumptions
    assert "sn positions assumed" in ls.illustrate("PC 16:0_18:1").assumptions
    assert "sn positions assumed" in ls.illustrate("LPC 18:1").assumptions


@pytest.mark.parametrize("name, note", [
    ("PC 16:0/18:1", None),                      # the name gives sn-1 / sn-2 and the drawing has them that way: nothing to say
    ("PC 18:1/16:0", "sn order in the name is not used"),
    ("PC(18:1/16:0)", "sn order in the name is not used"),
    ("PC 16:0_18:1", "sn positions assumed"),    # no order in the name
    ("PC 34:1", "sn positions assumed"),
    ("TG 16:0/18:1/18:2", None),
    ("TG 18:1/16:0/18:1", "sn order in the name is not used"),
    ("DG 16:0/18:1/0:0", None),
    ("DG 18:1/16:0/0:0", "sn order in the name is not used"),
    ("LPC 18:1/0:0", None),
    ("LPC 0:0/18:1", "sn order in the name is not used"),
    ("MG 0:0/18:1/0:0", "sn order in the name is not used"),
    ("PC 16:0/16:0", None),
    ("PC O-16:0/18:1", None),
    ("PC 18:1/O-16:0", "sn order in the name is not used"),      # the ether chain is always drawn at sn-1
    ("CL 18:2_18:2_18:2_18:2", None),
    ("CL 16:0_18:1_18:1_18:1", "sn positions assumed"),
])
def test_sn_order_given_by_a_slash_is_reported_as_used_or_not(name, note):
    notes = [a for a in ls.illustrate(name).assumptions if a.startswith("sn ")]
    assert notes == ([note] if note else []), name
    assert [a for a in run(name)["assumptions"]["items"] if a.startswith("sn ")] == notes          # the page passes the typed text on


def test_pg_head_group_glycerol_has_no_stereochemistry():
    for name in ("PG 16:0_18:1", "PG 34:1", "LPG 18:1", "PG O-16:0/18:1"):
        assert "Head-group glycerol drawn without stereochemistry" in ls.illustrate(name).assumptions, name
    assert "Inositol drawn without stereochemistry" in ls.illustrate("PI 38:4").assumptions
    assert "Head-group glycerol drawn without stereochemistry" not in ls.illustrate("PC 16:0_18:1").assumptions


def test_positions_given_in_the_name_are_reported_as_not_used():
    assert ls.positions_in_name("FA 18:1(11Z)") and ls.positions_in_name("FA 18:1(9Z);12OH") and ls.positions_in_name("FA 20:4;5Ep")
    assert not ls.positions_in_name("PC 16:0/18:1") and not ls.positions_in_name("FA 18:2;(2OH)") and not ls.positions_in_name("TG 52:3")
    items = run("FA 18:1(11Z)")["assumptions"]["items"]
    assert items[0] == "Positions given in the name are not used"
    assert "Positions given in the name are not used" not in run("FA 18:1")["assumptions"]["items"]


def test_unsupported_names_are_named_by_kind():
    with pytest.raises(ls.NoStructureForClass):
        ls.illustrate("ST 27:1;O")
    with pytest.raises(ls.Unsupported):
        ls.illustrate("Phosphatidylcholine")
    huge = run("TG 200:50")                                    # too large to split by the default rule: names still answer
    assert "structure" not in huge and huge["no_structure"]["items"] == ["No illustrative structure for this name."]


def test_a_drawing_never_takes_long():
    t0 = time.time()
    for name in ("CL 80:12", "TG 66:9", "CL 72:8", "TG 18:1_18:1_18:1"):
        run(name)
    assert time.time() - t0 < 15


# ---------------------------------------------------------------- the drawing (SVG)
@pytest.mark.parametrize("name", ["PC 16:0/18:1", "PC 34:1", "Cer 24:0_18:1;OH2", "FA 20:4;OH", "CE 18:1", "TG 52:3"])
def test_svg_is_valid_and_themeable(name):
    svg = run(name)["structure"]["svg"]
    xml.dom.minidom.parseString(svg.encode())
    assert len(svg) > 5000
    # the dark theme recolours "#RRGGBB" values only (web/src/darkColors.ts): no rgb() colours, no named colours other than none
    assert "rgb(" not in svg and "rgba(" not in svg
    fills = set(re.findall(r"(?:fill|stroke):\s*([^;'\"\s]+)", svg)) | set(re.findall(r"(?:fill|stroke)='([^']+)'", svg))
    assert all(f == "none" or re.fullmatch(r"#[0-9A-Fa-f]{6}", f) for f in fills), fills


def test_the_key_lists_only_what_is_highlighted():
    def colours(name):
        svg = run(name)["structure"]["svg"].upper()
        return {key: ("#%02X%02X%02X" % tuple(int(255 * c) for c in colour)) in svg for key, colour, _ in draw.KEY}      # RDKit truncates to 8 bit
    assert colours("PC 16:0/18:1") == {"db": True, "fg": False, "split": False}
    assert colours("FA 20:4;OH") == {"db": True, "fg": True, "split": False}
    assert colours("PC 34:1") == {"db": True, "fg": False, "split": True}
    assert colours("PC 16:0_16:0") == {"db": False, "fg": False, "split": False}


def test_highlight_colours_are_theme_colours_not_warning_colours():
    """Tints of blue / teal / violet; amber and coral stay for warnings and errors."""
    import colorsys
    from mzdesk.shared import colors

    def hue(rgb):
        return colorsys.rgb_to_hsv(*rgb)[0] * 360

    def near(a, b):
        return min(abs(a - b), 360 - abs(a - b))

    theme = [hue(colors._rgb(c)) for c in (colors.DATA_BLUE, colors.DATA_TEAL, colors.DATA_VIOLET)]
    warn = [hue(colors._rgb(c)) for c in (colors.DATA_AMBER, colors.DATA_CORAL)]
    for _, colour, _ in draw.KEY:
        assert min(near(hue(colour), h) for h in theme) < 6, colour
        assert min(near(hue(colour), h) for h in warn) > 25, colour


def test_the_drawing_is_deterministic():
    assert run("TG 52:3")["structure"]["svg"] == run("TG 52:3")["structure"]["svg"]


def test_the_page_icon_exists():
    icons = (HERE.parent / "web/src/assets/icons.svg").read_text(encoding="utf-8")
    assert 'id="lipid_short_name"' in icons
    assert '"lipid_short_name"' in (HERE.parent / "web/src/components/Icon.tsx").read_text(encoding="utf-8")
