"""Known methods 页的测试：规范化与精确匹配、相似度分数、JP18（随包数据，始终运行）、RepoRT（用 2-3 个假数据集建的小索引，不联网）、页面输出、错误输入。

真实 RepoRT 的测试在没安装时跳过（python -m mzdesk.methods.report_setup）。
"""
import time
from dataclasses import fields
from pathlib import Path

import numpy as np
import pytest
from fastapi.testclient import TestClient
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

from mzdesk.api.app import app
from mzdesk.methods import compound_index as ci
from mzdesk.methods import jp18, known_methods as km, report, report_setup
from mzdesk.spec import MODULES

client = TestClient(app)
JP_ONLY = [km.SOURCE_JP]
CAPTOPRIL = "C[C@H](CS)C(=O)N1CCC[C@H]1C(=O)O"
CAPTOPRIL_FLAT = "CC(CS)C(=O)N1CCCC1C(=O)O"
CAFFEINE = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
NOT_IN_JP = "CC(C)(C)c1ccc(cc1)C(=O)NCCN1CCOCC1"               # made-up amide, in neither source
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
S_IBUPROFEN = "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O"
XANTHINES = ["Cn1c(=O)c2[nH]cnc2n(C)c1=O", "Cn1cnc2c1c(=O)[nH]c(=O)n2C", "Cn1c(=O)[nH]c2ncn(C)c2c1=O", "O=c1[nH]c(=O)c2[nH]cnc2[nH]1",
             "Cn1c(=O)[nH]c(=O)c2[nH]cnc21"]                   # theophylline, theobromine, paraxanthine, xanthine, 3-methylxanthine
needs_real_report = pytest.mark.skipif(not report.installed(), reason="RepoRT not installed (python -m mzdesk.methods.report_setup)")


@pytest.fixture(scope="session")
def jp_dir(tmp_path_factory):
    """JP18 index cache of this test session (built once, in a temp folder, so the user's cache is not touched)."""
    return tmp_path_factory.mktemp("jp_cache")


def _tsv(path, header, rows):
    path.write_text("\t".join(header) + "\n" + "".join("\t".join(str(c) for c in row) + "\n" for row in rows), encoding="utf-8")


def _dataset(root, ds_id, name, kind, url, meta, gradient, entries):
    d = root / "processed_data" / ds_id
    d.mkdir(parents=True)
    _tsv(d / f"{ds_id}_info.tsv", ["id", "name", "method.type", "source", "url", "pmid", "authors", "missing information", "comment"],
         [[ds_id, name, kind, "Publication - Fake", url, "", "A. Author", "", ""]])
    _tsv(d / f"{ds_id}_metadata.tsv", list(meta), [list(meta.values())])
    _tsv(d / f"{ds_id}_gradient.tsv", ["t [min]", "A [%]", "B [%]", "C [%]", "D [%]", "flow rate [ml/min]"], gradient)
    _tsv(d / f"{ds_id}_rtdata_canonical_success.tsv", ["id", "name", "formula", "rt", "smiles.std", "inchi.std", "inchikey.std"],
         [[f"{ds_id}_{i:05d}", n, "", rt, smi, "", ""] for i, (n, rt, smi) in enumerate(entries, 1)])


@pytest.fixture(scope="session")
def fake_report(tmp_path_factory):
    """A tiny RepoRT index built from three fake datasets (functional-group charge classes: no MolGpKa, no network)."""
    root = tmp_path_factory.mktemp("fake_report_raw")
    rp = {"column.name": "Test C18 RP", "column.usp.code": "L1", "column.length": 100, "column.id": 2.1, "column.particle.size": 1.7,
          "column.temperature": 40, "column.flowrate": 0.4, "eluent.A.h2o": 100, "eluent.A.formic": 0.1, "eluent.A.formic.unit": "%", "eluent.A.pH": 3,
          "eluent.B.acn": 100, "eluent.B.formic": 0.1, "eluent.B.formic.unit": "%", "eluent.B.pH": 3}
    hilic = {"column.name": "Test HILIC", "column.usp.code": "L3", "column.length": 150, "column.id": 2.1, "column.particle.size": 1.7,
             "column.temperature": 35, "column.flowrate": 0.3, "eluent.A.h2o": 100, "eluent.A.nh4form": 10, "eluent.A.nh4form.unit": "mM",
             "eluent.B.acn": 95, "eluent.B.h2o": 5, "eluent.B.nh4form": 10, "eluent.B.nh4form.unit": "mM"}
    bare = {"column.name": "Test C8"}
    grad = [[0, 95, 5, 0, 0, 0.4], [20, 5, 95, 0, 0, 0.4], [25, 5, 95, 0, 0, 0.4]]
    _dataset(root, "0901", "FAKE_RP", "RP", "https://doi.org/10.1000/fake.rp", rp, grad,
             [("caffeine", 3.2, CAFFEINE), ("theophylline", 2.9, XANTHINES[0]), ("S-ibuprofen", 9.1, S_IBUPROFEN), ("captopril", 7.7, CAPTOPRIL)])
    _dataset(root, "0902", "FAKE_HILIC", "HILIC", "10.1000/fake.hilic", hilic, [[0, 0, 100, 0, 0, 0.3], [10, 50, 50, 0, 0, 0.3]],
             [("caffeine", 1.1, CAFFEINE), ("caffeine", 1.4, CAFFEINE), ("glycine", 6.5, "NCC(=O)O")])
    _dataset(root, "0903", "FAKE_MANY", "RP", "", bare, [],
             [(n, 5 + i, s) for i, (n, s) in enumerate(zip(["theophylline", "theobromine", "paraxanthine", "xanthine", "3-methylxanthine"], XANTHINES))])
    out = tmp_path_factory.mktemp("fake_report_index")
    stats = report_setup.build_index(root, out, commit="f" * 40, method=ci.METHOD_FG, workers=1, log=lambda *_: None)
    assert stats["datasets"] == 3
    return out


# ---------------------------------------------------------------- config, normalisation, math
def test_every_weight_and_cutoff_is_in_the_config():
    c = ci.CONFIG
    assert (c.w_structure, c.w_properties) == (0.6, 0.4)
    assert (c.scale_logp, c.scale_tpsa, c.scale_mw, c.ph_acidic, c.ph_neutral) == (1.0, 30.0, 100.0, 2.7, 7.4)
    assert (c.fp_radius, c.fp_bits, c.min_score, c.max_similar, c.max_methods_per_similar, c.max_similar_per_dataset) == (2, 2048, 0.5, 10, 3, 2)
    assert c.min_structure == 0.25
    assert {f.name for f in fields(c)} >= {"prop_weights", "min_structure", "max_exact_per_source"}


def test_normalize_takes_the_largest_fragment_and_neutralises():
    base = ci.normalize(CAPTOPRIL)
    salt = ci.normalize("C[C@H](CS)C(=O)N1CCC[C@H]1C(=O)[O-].[Na+]")
    assert salt.inchikey == base.inchikey and salt.notes
    assert ci.normalize("CC(CS)C(=O)N1CCCC1C(=O)O").key14 == base.key14 != ""
    assert ci.normalize(CAPTOPRIL_FLAT).stereo_block != base.stereo_block
    with pytest.raises(ValueError):
        ci.normalize("not a smiles")
    with pytest.raises(ValueError):
        ci.normalize("")


def test_popcount_and_tanimoto_match_python_and_rdkit():
    rng = np.random.default_rng(0)
    words = rng.integers(0, 2**63, size=(40, 32), dtype=np.uint64)
    assert ci.popcount(words).sum(axis=1).tolist() == [sum(bin(int(w)).count("1") for w in row) for row in words]
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    smiles = [CAPTOPRIL, CAFFEINE, IBUPROFEN, NOT_IN_JP, "CCO"]
    fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in smiles]
    packed = np.stack([np.frombuffer(ci.fingerprint(Chem.MolFromSmiles(s)), dtype=np.uint64) for s in smiles])
    pops = ci.popcount(packed).sum(axis=1).astype(np.uint16)
    got = ci.tanimoto(packed, pops, packed[3], int(pops[3]))
    assert got == pytest.approx([DataStructs.TanimotoSimilarity(f, fps[3]) for f in fps], abs=1e-6)


def test_functional_group_charge_classes():
    def cls(smiles):
        return ci.functional_group_classes(ci.normalize(smiles).mol)
    assert cls(IBUPROFEN) == (0, -1)                          # carboxylic acid: neutral at 2.7, negative at 7.4
    assert cls("CCCCS(=O)(=O)O") == (-1, -1)                  # sulfonic acid
    assert cls("OP(O)(=O)OCC(O)CO") == (-1, -1)               # phosphoric acid
    assert cls("CCN(CC)CC") == (1, 1)                         # aliphatic amine
    assert cls("CCCCN=C(N)N") == (1, 1)                       # guanidine
    assert cls("C[N+](C)(C)CCO") == (1, 1)                    # quaternary N
    assert cls("NCC(=O)O") == (1, 0)                          # amino acid: + at 2.7, zwitterion (net 0) at 7.4
    assert cls(CAFFEINE) == (0, 0) and cls("CC(=O)NC") == (0, 0) and cls("c1ccccc1N") == (0, 0)


needs_weights = pytest.mark.skipif(not ci.molgpka_available()[0], reason="MolGpKa weights not installed (python -m mzdesk.platform_fit.pka_setup)")


@needs_weights
def test_molgpka_charge_classes_agree_with_the_engine():
    from mzdesk.platform_fit import pka_engine as pe
    for smi in (IBUPROFEN, "NCC(=O)O", "CCN(CC)CC", CAPTOPRIL):
        a = pe.analyze(smi)
        want = (int(np.sign(a.at(2.7)["charge"])), int(np.sign(a.at(7.4)["charge"])))
        assert ci.molgpka_classes(ci.normalize(smi).smiles, (2.7, 7.4)) == want, smi


# ---------------------------------------------------------------- JP18 only (bundled data)
def _jp(smiles, jp_dir):
    return km.search(smiles, JP_ONLY, jp_dir=jp_dir)


def test_captopril_exact_in_jp_through_the_monograph_substance(jp_dir):
    res = _jp(CAPTOPRIL, jp_dir)
    assert res.match == "exact" and res.searched == [km.SOURCE_JP]
    m = res.exact[0]
    # the only JP18 LC block of the captopril monograph measures the disulfide impurity: captopril is matched as the monograph substance
    assert m.source == "JP18" and m.compound == "Captopril" and m.match == f"exact · {SUBSTANCE_ONLY} · same stereo"
    assert m.column.startswith("C18, 3.9 × 300 mm, 10 µm, 25 °C") and m.detector == "UV 220 nm" and m.mobile.startswith("isocratic; methanol")
    assert m.reference[0] == "JP18 p. 623 (Official Monographs A–L)" and m.reference[1].endswith("000904440.pdf#page=225")
    assert m.block_id == "JP18-ML-0211" and "Detector: An ultraviolet absorption photometer" in m.condition and m.rt is None


def test_exact_ignores_stereo_and_salt_form_but_says_so(jp_dir):
    flat = _jp(CAPTOPRIL_FLAT, jp_dir)
    assert flat.match == "exact" and flat.exact[0].match == f"exact · {SUBSTANCE_ONLY} · stereo differs / not specified"
    salt = _jp("C[C@H](CS)C(=O)N1CCC[C@H]1C(=O)[O-].[Na+]", jp_dir)
    assert salt.exact and salt.exact[0].match.endswith("· same stereo") and any("largest fragment" in w for w in salt.warnings)


def test_crude_drug_marker_paeoniflorin_merges_identical_conditions(jp_dir):
    jp = jp18.load(cache_dir=jp_dir)
    rec = next(r for r in jp.records if r["block_id"] == "JP18-CD-0126")        # Peony Root, assay of paeoniflorin
    smiles = rec["analyte_structures"][0]["smiles"]
    res = _jp(smiles, jp_dir)
    assert res.match == "exact" and all(m.compound.lower() == "paeoniflorin" for m in res.exact)
    assert {m.match for m in res.exact} == {"exact · same stereo"}                  # paeoniflorin is what these tests measure: plain text
    n_records = sum(1 for r in jp.records if any(a.get("smiles") == smiles for a in r["analyte_structures"]))
    assert n_records >= 10 and 1 <= len(res.exact) < n_records                  # 17 monographs, a handful of distinct conditions
    assert any("more with identical conditions" in m.test for m in res.exact)
    assert {m.mobile.split(";")[0] for m in res.exact} == {"isocratic"} or len(res.exact) >= 1


def test_substance_only_methods_are_marked_listed_after_measuring_ones_and_not_merged(jp_dir):
    jp = jp18.load(cache_dir=jp_dir)
    idx = jp.index
    both = [c for c in range(len(idx)) if {0, 1} <= set(int(x) for x in idx.roles(c))]
    assert both, "expected a JP compound that some tests measure and other monograph tests only use as the substance"
    comp = both[0]
    res = _jp(idx.smiles[comp], jp_dir)
    flags = ["monograph substance" in m.match for m in res.exact]
    assert True in flags and False in flags and flags == sorted(flags)                # measuring tests first, substance-only after
    assert all("(test measures: " in m.match for m in res.exact if "monograph substance" in m.match)
    assert all(m.match in ("exact · same stereo", "exact · stereo differs / not specified") for m in res.exact if "monograph substance" not in m.match)
    assert all(r["analytes"] for r in jp.records if r["block_id"] in {m.block_id.split(" ")[0] for m in res.exact if "monograph substance" in m.match})


def test_jp_rows_take_analytes_and_substances_once_per_record_and_compound():
    rec = {"analyte_structures": [{"smiles": "CCO", "name": "ethanol"}, {"smiles": None, "name": "x"}],
           "substances": [{"smiles": "CCO", "name": "ethanol"}, {"smiles": "CCN", "name": "ethylamine"}]}
    rows = jp18._rows([rec])
    assert [(r[0], r[4]) for r in rows] == [("CCO", 0), ("CCO", 1), ("CCN", 1)]
    index, skipped, _ = ci.CompoundIndex.build(rows, {"source": "test"}, ci.METHOD_FG, dedupe_refs=True)
    assert skipped == 0 and len(index) == 2 and index.meta["n_entries"] == 2
    ethanol = index.exact(ci.normalize("CCO").key14)[0]
    amine = index.exact(ci.normalize("CCN").key14)[0]
    assert index.roles(ethanol).tolist() == [0] and index.roles(amine).tolist() == [1]          # measured wins; the substance-only one keeps its flag


def test_similar_hits_for_a_compound_not_in_jp(jp_dir):
    res = _jp(NOT_IN_JP, jp_dir)
    assert res.exact == [] and res.match in ("similar", "none")
    scores = [s.score for s in res.similar]
    assert scores == sorted(scores, reverse=True) and all(0.5 <= s <= 1.0 for s in scores)
    assert len(res.similar) <= 10 and all(1 <= len(s.methods) <= 3 for s in res.similar)
    assert all(m.match.startswith("similar") and m.similarity == s.score for s in res.similar for m in s.methods)


def test_similar_list_leaves_out_the_compound_itself_and_scores_stay_in_range(jp_dir):
    res = _jp(CAPTOPRIL, jp_dir)
    q = res.query.norm.key14
    assert res.similar and all(s.key14 != q for s in res.similar)
    assert all(0.0 <= s.tanimoto <= 1.0 and 0.5 <= s.score <= 1.0 for s in res.similar)
    jp = jp18.load(cache_dir=jp_dir)
    scores = jp.index.score(res.query)
    assert scores.total.min() >= 0 and scores.total.max() <= 1 + 1e-6
    own = jp.index.exact(q)
    assert len(own) == 1 and scores.total[own[0]] == pytest.approx(1.0, abs=1e-5)       # identical structure and properties


def test_min_score_and_list_size_follow_the_config(jp_dir):
    strict = ci.Config(min_score=0.9, max_similar=2)
    res = km.search(NOT_IN_JP, JP_ONLY, cfg=strict, jp_dir=jp_dir)
    assert len(res.similar) <= 2 and all(s.score >= 0.9 for s in res.similar)
    loose = km.search(CAFFEINE, JP_ONLY, cfg=ci.Config(max_similar=3, max_methods_per_similar=1), jp_dir=jp_dir)
    assert len(loose.similar) <= 3 and all(len(s.methods) == 1 for s in loose.similar)


HEPTYL_CAPTOPRIL = "CCCCCCCC(CS)C(=O)N1CCCC1C(=O)O"
DISULFIDE = "1,1′-[3,3′-dithiobis(2-methyl-1-oxopropyl)]-l-diproline"       # what JP18's captopril LC test measures
SUBSTANCE_ONLY = f"monograph substance (test measures: {DISULFIDE})"


def test_structure_floor_keeps_out_hits_that_only_match_in_properties(jp_dir):
    default = km.search(HEPTYL_CAPTOPRIL, JP_ONLY, jp_dir=jp_dir)
    assert default.similar and all(s.tanimoto >= 0.25 and s.score >= 0.5 for s in default.similar)
    assert not any("mitiglinide" in s.name.lower() for s in default.similar)
    captopril = next(s for s in default.similar if s.name == "Captopril")
    assert [m.match for m in captopril.methods if m.source == "JP18"] == [f"similar · {SUBSTANCE_ONLY}"]          # same marking for similar hits
    no_floor = km.search(HEPTYL_CAPTOPRIL, JP_ONLY, cfg=ci.Config(min_structure=0.0), jp_dir=jp_dir)
    weak = [s for s in no_floor.similar if s.tanimoto < 0.25]
    assert any("mitiglinide" in s.name.lower() for s in weak)                              # property agreement alone carried it over the cut-off
    notes = km.page_outputs(default)["similar"]["help"]
    assert "Tanimoto ≥ 0.25" in notes[0]


def test_jp_texts_for_gc_and_fluorescence_records(jp_dir):
    jp = jp18.load(cache_dir=jp_dir)
    gc = next(r for r in jp.records if r["technique"] == "GC" and r["carrier_gas"] and r["split_ratio"])
    assert jp18.mobile_phase_text(gc).startswith("GC; carrier gas") and "split" in jp18.mobile_phase_text(gc)
    fl = next(r for r in jp.records if r["detector_type"] == "fluorescence")
    assert jp18.detector_text(fl).startswith("fluorescence ex ")
    assert jp18.column_text({"packing_class": "C18", "column_id_mm": 4.6, "column_length_mm": 150, "particle_or_film_um": 5, "column_temp_c": 40}) == "C18, 4.6 × 150 mm, 5 µm, 40 °C"
    assert jp18.mobile_phase_text({"technique": "LC", "elution_mode": "isocratic", "organic_solvents": ["acetonitrile"], "mobile_phase_ph": 3}) == "isocratic; acetonitrile; pH 3.0"


def test_jp_index_cache_is_reused_and_rebuilt_when_stale(jp_dir):
    jp18._cache.clear()
    first = jp18.load(cache_dir=jp_dir)
    assert (Path(jp_dir) / jp18.CACHE_FILE).is_file()
    jp18._cache.clear()
    again = jp18.load(cache_dir=jp_dir)
    assert again.index.meta == first.index.meta and np.array_equal(again.index.ik14, first.index.ik14)
    ci.CompoundIndex.load(Path(jp_dir) / jp18.CACHE_FILE).meta["data_sha"]


# ---------------------------------------------------------------- RepoRT (fake datasets)
def test_report_index_roundtrip_and_merging(fake_report):
    rep = report.load(fake_report)
    assert rep is not None and len(rep.datasets) == 3 and rep.commit == "f" * 40
    idx = rep.index
    assert idx.charge_method == ci.METHOD_FG and idx.meta["n_compounds"] == len(idx) == 9      # caffeine (3 entries) and theophylline (2) count once
    again = ci.CompoundIndex.load(Path(fake_report) / report.INDEX_FILE)
    assert np.array_equal(again.fp, idx.fp) and again.smiles[0] == idx.smiles[0] and again.names[0] == idx.names[0]
    idx.check()
    with pytest.raises(ValueError):
        idx.check(ci.Config(fp_bits=1024))
    assert report.load(Path(fake_report) / "nothing here") is None


def test_report_texts():
    ds = {"id": "0001", "name": "X", "type": "RP", "url": "https://doi.org/10.1000/x", "source": "S",
          "column": {"name": "Waters ACQUITY UPLC HSS T3", "usp": "L1", "length": 150, "id": 2.1, "particle": 1.8, "temp": 40, "flow": 0.28},
          "eluents": {"A": {"solvents": {"h2o": 100}, "additives": [["formic", 0.1, "%"]], "pH": 3},
                      "B": {"solvents": {"acn": 100}, "additives": [["formic", 0.1, "%"], ["nh4form", 10, "mM"]], "pH": None}},
          "gradient": {"t": [0, 1, 20, 25], "A": [95, 95, 5, 5], "B": [5, 5, 95, 95]}}
    assert report.column_text(ds) == "Waters ACQUITY UPLC HSS T3 (L1), 2.1 × 150 mm, 1.8 µm, 40 °C"
    assert report.eluents_text(ds) == "A: water + 0.1% formic (pH 3); B: ACN + 0.1% formic + 10 mM ammonium formate"
    assert report.gradient_text(ds) == "5→95% B in 20 min"
    assert report.mobile_phase_text(ds).endswith("5→95% B in 20 min; 0.28 mL/min")
    assert report.reference(ds) == ("RepoRT 0001 · 10.1000/x", "https://doi.org/10.1000/x")
    assert report.reference({"id": "0002", "url": "", "source": "PredRet"}) == ("RepoRT 0002 · PredRet", None)
    assert report.reference({"id": "0003", "url": "10.1016/j.x.1", "source": "S"})[1] == "https://doi.org/10.1016/j.x.1"
    ds["gradient"] = {"t": [0, 10], "A": [40, 40], "B": [60, 60]}
    assert report.gradient_text(ds) == "isocratic 40% A, 60% B"
    ds["gradient"] = {"t": [0, 5, 20], "A": [5, 5, 43], "B": [95, 95, 57]}                  # HILIC: B falls
    assert report.gradient_text(ds) == "95→57% B in 20 min"
    ds["gradient"] = {}
    ds["start"], ds["end"] = {"A": 95, "B": 5}, {"A": 5, "B": 95}
    assert report.gradient_text(ds) == "5→95% B"


def test_report_exact_hits_list_every_dataset_with_retention_times(fake_report):
    res = km.search(CAFFEINE, [km.SOURCE_REPORT], report_path=fake_report)
    assert res.match == "exact" and res.report_installed and res.n_datasets == 3
    rows = {m.test.split()[0]: m for m in res.exact}
    assert set(rows) == {"0901", "0902"}
    rp, hilic = rows["0901"], rows["0902"]
    assert rp.source == "RepoRT" and rp.detector == "MS" and rp.rt == "3.2" and hilic.rt == "1.1 / 1.4"
    assert rp.column == "Test C18 RP (L1), 2.1 × 100 mm, 1.7 µm, 40 °C"
    assert rp.mobile == "A: water + 0.1% formic (pH 3); B: ACN + 0.1% formic (pH 3); 5→95% B in 20 min; 0.4 mL/min"
    assert rp.reference == ("RepoRT 0901 · 10.1000/fake.rp", "https://doi.org/10.1000/fake.rp")
    assert hilic.reference[1] == "https://doi.org/10.1000/fake.hilic" and hilic.test == "0902 FAKE_HILIC (HILIC)"
    assert "10 mM ammonium formate" in hilic.mobile and hilic.match == "exact · same stereo"


def test_report_exact_stereo_note(fake_report):
    flat = km.search(IBUPROFEN, [km.SOURCE_REPORT], report_path=fake_report)
    assert flat.match == "exact" and flat.exact[0].match == "exact · stereo differs / not specified"
    same = km.search(S_IBUPROFEN, [km.SOURCE_REPORT], report_path=fake_report)
    assert same.exact[0].match == "exact · same stereo"
    other = km.search("CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O", [km.SOURCE_REPORT], report_path=fake_report)
    assert other.exact[0].match == "exact · stereo differs / not specified"


def test_report_similar_cap_per_dataset_and_per_compound(fake_report):
    res = km.search(CAFFEINE, [km.SOURCE_REPORT], report_path=fake_report)
    assert res.similar and all(s.key14 != res.query.norm.key14 for s in res.similar)
    per_dataset = {}
    for m in (m for s in res.similar for m in s.methods):
        per_dataset[m.test] = per_dataset.get(m.test, 0) + 1
    assert per_dataset and max(per_dataset.values()) <= 2
    assert all(1 <= len(s.methods) <= 3 for s in res.similar)
    assert "0903 FAKE_MANY (RP)" in per_dataset                                    # the bare dataset is used, but only twice
    scores = [s.score for s in res.similar]
    assert scores == sorted(scores, reverse=True) and all(0.5 <= x <= 1 for x in scores)
    assert all(m.rt for s in res.similar for m in s.methods)                      # RT of the similar compound in that dataset


def test_report_missing_metadata_shows_dashes(fake_report):
    res = km.search("CC(C)(C)c1ccc(cc1)C(=O)NCCN1CCOCC1", [km.SOURCE_REPORT], report_path=fake_report)
    assert res.exact == []
    xanth = km.search("Cn1c(=O)c2[nH]c(Cl)nc2n(C)c1=O", [km.SOURCE_REPORT], report_path=fake_report)           # 8-chlorotheophylline
    bare = [m for s in xanth.similar for m in s.methods if m.test.startswith("0903")]
    assert bare and all(m.column == "Test C8" and m.mobile == "—" for m in bare)


def test_both_sources_together(fake_report, jp_dir):
    res = km.search(CAPTOPRIL, km.SOURCES, report_path=fake_report, jp_dir=jp_dir)           # in JP18 and in the fake dataset 0901
    assert {m.source for m in res.exact} == {"JP18", "RepoRT"}
    assert [m.source for m in res.exact] == sorted((m.source for m in res.exact), key=("JP18", "RepoRT").index)
    out = km.page_outputs(res)
    assert out["summary"]["values"]["Sources"]["value"] == "JP18, RepoRT"
    assert {_cells(r)[0]["badge"] for r in out["exact"]["rows"]} == {"JP18", "RepoRT"}


# ---------------------------------------------------------------- the page
BLOCKS = {"summary", "exact", "similar"}
COLUMNS = ["Source", "Test / dataset", "Column", "Mobile phase", "Detector", "RT (min)"]


def _run(inputs):
    return client.post("/api/modules/known_methods/run", json=inputs)


def _inputs(smiles, sources):
    return {"structure": {"mode": "SMILES / InChI", "value": smiles}, "sources": sources}


def _lines(help_):
    return [help_] if isinstance(help_, str) else list(help_ or [])


def _cells(row):
    return row["cells"] if isinstance(row, dict) else row


def test_spec_page_is_live_with_the_agreed_blocks():
    mod = next(m for m in MODULES if m["id"] == "known_methods")
    assert mod["live"] and {b["id"] for b in mod["outputs"] if b.get("id")} == BLOCKS
    sources = next(f for f in mod["inputs"] if f["id"] == "sources")
    assert sources["options"] == ["Japanese Pharmacopoeia", "RepoRT"]
    summary, exact, similar = mod["outputs"]
    assert (summary["type"], summary["style"], summary["title"]) == ("cards", "metric", "Summary")
    assert [i["label"] for i in summary["items"]] == ["Match", "Methods found", "Sources"] and all(i["help"] for i in summary["items"])
    assert (exact["type"], exact["title"], exact["row_help"], exact["columns"]) == ("table", "Exact matches", True, COLUMNS)
    assert (similar["type"], similar["title"], similar["row_help"], similar["columns"]) == ("grouptable", "Similar compounds", True, COLUMNS)
    assert mod["example"]["structure"]["value"] == "captopril"


def test_page_for_captopril_jp_only():
    body = _run(_inputs(CAPTOPRIL, ["Japanese Pharmacopoeia"])).json()
    out = body["outputs"]
    assert set(out) == BLOCKS
    cards = out["summary"]["values"]
    n_methods = len(out["exact"]["rows"]) + sum(len(g["rows"]) for g in out["similar"]["groups"])
    assert cards["Match"] == {"value": "exact", "note": "monograph substance only; same stereo"} and cards["Methods found"]["value"] == str(n_methods)
    assert cards["Sources"]["value"] == "JP18" and "note" not in cards["Sources"]
    assert out["exact"]["count"] == len(out["exact"]["rows"]) == 1
    row = out["exact"]["rows"][0]
    assert row["cells"][:2] == [{"badge": "JP18", "tone": "jp"}, {"text": "Purity (4) · p. 623", "url": "https://www.mhlw.go.jp/content/11120000/000904440.pdf#page=225"}]
    assert row["cells"][2:] == ["C18, 3.9 × 300 mm, 10 µm, 25 °C", "isocratic · methanol", "UV 220 nm", "—"]
    help_ = _lines(row["help"])                                                                           # match details, test, reference, printed conditions
    assert help_[0] == f"Monograph substance: the test measures {DISULFIDE}."
    assert "Reference: JP18 p. 623 (Official Monographs A–L), block JP18-ML-0211" in help_
    assert any(h.startswith("Printed conditions: Detector: An ultraviolet absorption photometer (wavelength: 220 nm). Column:") for h in help_)
    assert out["similar"]["count"] == len(out["similar"]["groups"]) and out["similar"]["groups"]
    notes = _lines(out["similar"]["help"])                                                               # the score formula and the query, behind the block's ?
    assert notes[0].startswith("Score = 0.6 × Tanimoto") and notes[1].startswith("Query: logP")
    scores = []
    for g in out["similar"]["groups"]:
        assert set(g) == {"title", "rows", "meta", "help"} and g["rows"]
        score, _, count = g["meta"].partition(" · ")
        assert count == f"{len(g['rows'])} method{'' if len(g['rows']) == 1 else 's'}" and 0.5 <= float(score) <= 1
        scores.append(float(score))
        gh = " ".join(_lines(g["help"]))
        assert "Tanimoto" in gh and "ΔlogP" in gh and "charge class" in gh.lower() and "In JP18" in gh and "SMILES: " in gh
        assert all(len(_cells(r)) == 6 and _cells(r)[0] == {"badge": "JP18", "tone": "jp"} for r in g["rows"])
    assert scores == sorted(scores, reverse=True)
    licence = " ".join(_lines(cards["Sources"]["help"]))
    assert "1924 LC / GC methods" in licence and "Ministerial Notification" in licence


def test_jp_test_cell_names_the_monograph_only_when_it_is_not_the_compound(jp_dir):
    jp = jp18.load(cache_dir=jp_dir)
    paeoniflorin = next(r for r in jp.records if r["block_id"] == "JP18-CD-0126")["analyte_structures"][0]["smiles"]       # the marker of Peony Root and 16 more
    res = km.search(paeoniflorin, JP_ONLY, jp_dir=jp_dir)
    texts = [_cells(r)[1]["text"] for r in km.page_outputs(res)["exact"]["rows"]]
    assert any(t.startswith("Peony Root: Assay · p. ") for t in texts)                                               # monograph differs from the compound: it is named
    assert any("(+" in t for t in texts)                                                                             # merged identical conditions: (+n)
    own = km.page_outputs(_jp(CAPTOPRIL, jp_dir))["exact"]["rows"][0]
    assert _cells(own)[1]["text"] == "Purity (4) · p. 623"                                                           # captopril's own monograph: just the test and the page


def test_page_reports_a_missing_report_install(tmp_path, monkeypatch):
    monkeypatch.setenv("MZDESK_REPORT_DIR", str(tmp_path / "empty"))
    body = _run(_inputs(CAPTOPRIL, ["Japanese Pharmacopoeia", "RepoRT"])).json()
    card = body["outputs"]["summary"]["values"]["Sources"]
    assert card["value"] == "JP18" and "python -m mzdesk.methods.report_setup" in card["note"]
    assert any("python -m mzdesk.methods.report_setup" in h for h in card["help"])
    assert body["outputs"]["summary"]["values"]["Match"]["value"] == "exact"


def test_page_with_the_fake_report(fake_report, monkeypatch):
    monkeypatch.setenv("MZDESK_REPORT_DIR", str(fake_report))
    out = _run(_inputs(CAFFEINE, ["RepoRT"])).json()["outputs"]
    assert out["summary"]["values"]["Match"]["value"] == "exact" and out["summary"]["values"]["Sources"]["value"] == "RepoRT"
    exact = out["exact"]["rows"]
    assert out["exact"]["count"] == 2 and [_cells(r)[1]["text"].split()[0] for r in exact] == ["0901", "0902"]
    c0, c1 = _cells(exact[0]), _cells(exact[1])
    assert c0[0] == {"badge": "RepoRT", "tone": "report"} and c0[1] == {"text": "0901 FAKE_RP (RP)", "url": "https://doi.org/10.1000/fake.rp"}
    assert c0[4] == "MS" and c0[5] == "3.2" and c1[5] == "1.1 / 1.4"
    help0 = _lines(exact[0]["help"])
    assert "Reference: RepoRT 0901 · 10.1000/fake.rp" in help0 and "Gradient, % B: 0 min 5, 20 min 95, 25 min 95" in help0       # the gradient behind the row's ?
    assert any("commit ffffffffff" in h for h in out["summary"]["values"]["Sources"]["help"])
    assert not any("Monograph" in h for h in help0)                       # same-stereo hit: no match note


def test_exact_row_help_says_when_the_stereo_differs(fake_report, monkeypatch):
    monkeypatch.setenv("MZDESK_REPORT_DIR", str(fake_report))
    out = _run(_inputs(IBUPROFEN, ["RepoRT"])).json()["outputs"]
    assert "Stereo differs / not specified." in _lines(out["exact"]["rows"][0]["help"])
    same = _run(_inputs(S_IBUPROFEN, ["RepoRT"])).json()["outputs"]
    assert "Stereo differs / not specified." not in _lines(same["exact"]["rows"][0]["help"])


def test_similar_groups_with_the_fake_report(fake_report, monkeypatch):
    monkeypatch.setenv("MZDESK_REPORT_DIR", str(fake_report))
    out = _run(_inputs(CAFFEINE, ["RepoRT"])).json()["outputs"]
    groups = out["similar"]["groups"]
    assert out["similar"]["count"] == len(groups) and groups
    titles = [g["title"] for g in groups]
    assert "theophylline" in titles and "caffeine" not in titles                                                # the query itself is not listed
    for g in groups:
        assert all(_cells(r)[0] == {"badge": "RepoRT", "tone": "report"} and _cells(r)[5] for r in g["rows"])
        assert "In RepoRT" in " ".join(_lines(g["help"]))


def test_page_with_no_match_says_none(monkeypatch, fake_report):
    monkeypatch.setenv("MZDESK_REPORT_DIR", str(fake_report))
    body = _run(_inputs("C" * 40, ["RepoRT"])).json()
    out = body["outputs"]
    assert out["summary"]["values"]["Match"]["value"] == "none" and out["exact"] == {"rows": [], "count": 0}
    assert out["similar"]["groups"] == [] and out["similar"]["count"] == 0
    assert _lines(out["similar"]["help"])[-1].startswith("No compound reached the score cut-off")


def test_group_title_prefers_a_readable_name_and_cuts_a_smiles_name():
    assert ci.is_id_like("HMDB0034131") and ci.is_id_like("CHEBI:1234") and ci.is_id_like("62571-86-2") and not ci.is_id_like("caffeine")
    assert km._title("caffeine") == "caffeine" and len(km._title("C" * 100)) == 60 and km._title("C" * 100).endswith("…")


def test_page_errors():
    assert _run(_inputs("not a smiles", ["Japanese Pharmacopoeia"])).status_code == 400
    assert _run(_inputs(CAPTOPRIL, [])).status_code == 400
    r = _run(_inputs(CAPTOPRIL, ["International Pharmacopoeia"]))
    assert r.status_code == 400 and "Unknown source" in r.json()["detail"]
    assert _run({"structure": {"mode": "SMILES / InChI", "value": ""}, "sources": ["RepoRT"]}).status_code == 400


# ---------------------------------------------------------------- the real RepoRT (skipped when not installed)
@needs_real_report
def test_real_report_caffeine_and_speed():
    rep = report.load()
    assert len(rep.datasets) > 300 and len(rep.index) > 50000
    res = km.search(CAFFEINE, [km.SOURCE_REPORT])
    assert res.match == "exact" and all(m.detector == "MS" and m.reference[0].startswith("RepoRT ") for m in res.exact)
    assert any(m.rt for m in res.exact) and len(res.exact) <= ci.CONFIG.max_exact_per_source
    q = ci.make_query(NOT_IN_JP)
    started = time.perf_counter()
    scores = rep.index.score(q)
    rows = rep.index.candidates(scores, q.norm.key14)
    assert time.perf_counter() - started < 0.5, "Tanimoto + property score over the whole index should take well under 0.5 s"
    assert len(rows) >= 0 and scores.total.shape == (len(rep.index),)


@needs_real_report
def test_real_report_similar_list_obeys_the_caps():
    res = km.search("CC(C)(C)c1ccc(cc1)C(=O)NCCN1CCOCC1", [km.SOURCE_REPORT])
    assert res.exact == [] and res.similar
    assert len(res.similar) <= 10 and all(1 <= len(s.methods) <= 3 for s in res.similar)
    per_dataset = {}
    for s in res.similar:
        for m in s.methods:
            per_dataset[m.test] = per_dataset.get(m.test, 0) + 1
    assert max(per_dataset.values()) <= 2


def test_mobile_phase_cell_has_one_line_per_eluent():
    from mzdesk.methods.known_methods import _phase_cell
    assert _phase_cell("A: water + 0.1% formic (pH 3); B: ACN + 0.1% formic (pH 3); 5→95% B in 20 min; 0.4 mL/min") == \
        "A: water + 0.1% formic (pH 3)\nB: ACN + 0.1% formic (pH 3)\n5→95% B in 20 min · 0.4 mL/min"
    assert _phase_cell("isocratic; acetonitrile; pH 3.0") == "isocratic · acetonitrile · pH 3.0"
    assert _phase_cell("—") == "—"
