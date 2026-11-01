"""CCS 页的测试：容差按平台对选择、只比较同加合物、多个观测值取最近、He / TWIMS N2→He / FAIMS 不可比、预测后备（假模型）、
索引没装的提示、runner 输出的块形状、加合物写法映射、建索引（假数据集）、训练（合成小数据）。全部离线，不下载。

真实索引的测试在没安装时跳过（python -m mzdesk.orthogonal.ccs_setup）。
"""
import dataclasses
import gzip
import json
import shutil
import subprocess

import numpy as np
import pytest
from fastapi.testclient import TestClient

from mzdesk import spec
from mzdesk.api import verdicts
from mzdesk.api.app import app
from mzdesk.methods import compound_index as ci
from mzdesk.orthogonal import ccs, ccs_setup
from mzdesk.spec import MODULES

client = TestClient(app)
CAFFEINE = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
ASPIRIN = "CC(=O)OC1=CC=CC=C1C(=O)O"
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
DECANOIC = "CCCCCCCCCC(=O)O"
needs_real_ccs = pytest.mark.skipif(not ccs.installed(), reason="CCS index not installed (python -m mzdesk.orthogonal.ccs_setup)")
BLOCKS = {"summary", "reference", "other_ions", "prediction", "notes"}


def _meta(path):
    return json.loads((path / ccs.META_FILE).read_text(encoding="utf-8"))


# ---------------------------------------------------------------- 假索引：3 个假数据集，经过 ccs_setup.build_rows
def _datasets():
    def meta(tag, kind, method):
        return {"src_tag": tag, "ccs_type": kind, "ccs_method": method, "reference": f"Fake, A. et al. {tag} (2020).", "url": f"https://doi.org/10.1000/{tag}"}
    dt = [{"name": "Caffeine", "adduct": "[M+H]+", "mz": "195.0877", "ccs": "140.5", "smi": CAFFEINE},
          {"name": "Caffeine", "adduct": "[M+Na]+", "mz": 217.0696, "ccs": 142.9, "smi": CAFFEINE},
          {"name": "Caffeine", "adduct": "[M+H]+", "mz": 195.0877, "ccs": 140.5, "smi": CAFFEINE},          # 重复：跳过
          {"name": "Nameless", "adduct": "[M+H]+", "mz": 100.0, "ccs": 120.0},                               # 没有 SMILES：跳过
          {"name": "Bad", "adduct": "[M+H]+", "mz": 100.0, "ccs": 120.0, "smi": "not a smiles"},          # 结构用不了：跳过
          {"name": "Broken", "adduct": "[M+H]+", "mz": 100.0, "ccs": "n/a", "smi": CAFFEINE}]             # 数值不对：跳过
    tw = [{"name": "Caffeine", "adduct": "[M+H]+", "mz": 195.0877, "ccs": 138.3, "smi": CAFFEINE},
          {"name": "Aspirin", "adduct": "[M+H]+", "mz": 181.0495, "ccs": 134.7, "smi": ASPIRIN},
          {"name": "Aspirin", "adduct": "[M+H20-H]-", "mz": 197.0, "ccs": 130.0},                          # 加合物写错，名字在 SMILES 缓存里
          {"name": "Cached acid", "adduct": "[M+]+", "mz": 171.1385, "ccs": 150.0}]
    tims = [{"name": "Caffeine", "adduct": "[M+H]+", "mz": 195.0877, "ccs": 141.0, "smi": CAFFEINE}]
    return [(meta("fakeDT", "DT", "single field, calibrated"), dt), (meta("fakeTW", "TW", "calibrated with polyalanine"), tw),
            (meta("fakeTIMS", "TIMS", "single field, calibrated"), tims)], {"Cached acid": DECANOIC, "Aspirin": ASPIRIN}


MODEL_INFO = {"name": "SVR", "params": {"C": 10.0}, "adducts": {"[M+H]+": 500, "[M+Na]+": 300, "[M+K]+": 120}, "min_adduct_rows": 100, "sklearn": "0",
              "split": {"by": "connectivity key", "train_rows": 700, "val_rows": 150, "test_rows": 150, "trainval_compounds": 400, "test_compounds": 80},
              "test": {"n": 150, "median_pct": 1.2, "q95_pct": 5.51, "by_type": {}}, "tolerance_pct": 5.6,
              "domain": {"threshold": 0.5, "train_compounds": 400, "test_inside": {"n": 130, "median_pct": 1.0, "q95_pct": 4.5},
                         "test_outside": {"n": 20, "median_pct": 2.5, "q95_pct": 9.0}}}


def make_index(path, with_model=False):
    datasets, cache = _datasets()
    rows, sources, stats = ccs_setup.build_rows(datasets, cache, log=lambda *_: None)
    ccs.save_index(path / ccs.INDEX_FILE, rows, sources)
    meta = {"commit": "f" * 40, "licence": "MIT", "built": "2026-01-01", "schema": ccs.SCHEMA, **stats}
    if with_model:
        meta["model"] = MODEL_INFO
    ccs_setup.write_meta(path, meta)
    return stats


@pytest.fixture
def fake_ccs(tmp_path, monkeypatch):
    """假索引（没有模型）。"""
    monkeypatch.setenv("MZDESK_CCS_DIR", str(tmp_path))
    make_index(tmp_path)
    return tmp_path


class _Estimator:
    """假模型：永远预测 150.0，并记下收到的特征矩阵。"""
    def __init__(self):
        self.seen = []

    def predict(self, X):
        self.seen.append(np.asarray(X))
        return np.full(len(X), 150.0)


@pytest.fixture
def fake_model(tmp_path, monkeypatch):
    """假索引 + 假预测模型（训练集只有咖啡因和阿司匹林）。"""
    monkeypatch.setenv("MZDESK_CCS_DIR", str(tmp_path))
    make_index(tmp_path, with_model=True)
    fps = np.stack([np.frombuffer(ci.fingerprint(ci.normalize(s).mol), dtype=np.uint8) for s in (CAFFEINE, ASPIRIN)]).view(np.uint64)
    est = _Estimator()
    monkeypatch.setattr(ccs, "load_predictor", lambda path=None, index=None: (ccs.Predictor(est, MODEL_INFO, fps), ""))
    return est


def _inputs(smiles=CAFFEINE, adduct="[M+H]+", observed="", platform="DTIMS (single field)", gas="N₂", tw="N₂ reference values"):
    return {"structure": {"mode": "SMILES / InChI", "value": smiles}, "adduct": adduct, "observed": observed, "platform": platform, "gas": gas, "tw_calibration": tw}


def _run(**kw):
    return client.post("/api/modules/ccs/run", json=_inputs(**kw))


def _card(body, label):
    v = body["outputs"]["summary"]["values"][label]
    return v if isinstance(v, dict) else {"value": v}


def _cells(row):
    return row["cells"] if isinstance(row, dict) else row


# ---------------------------------------------------------------- spec 和 Config
def test_the_page_is_live_and_orthogonal_is_ccs_only():
    mod = next(m for m in MODULES if m["id"] == "ccs")
    assert mod["live"] and mod["subgroup"] == "orthogonal" and mod["example"]
    assert [m["id"] for m in MODULES if m["subgroup"] == "orthogonal"] == ["ccs"]
    assert "rt" not in {m["id"] for m in MODULES}
    assert [f["id"] for f in mod["inputs"]] == ["structure", "adduct", "observed", "platform", "gas", "tw_calibration"]
    assert {b["id"] for b in mod["outputs"] if b.get("id")} == BLOCKS
    assert all(b.get("help") for b in mod["outputs"]), "every block has a short help"
    by_id = {f["id"]: f for f in mod["inputs"]}
    assert by_id["platform"]["options"] == list(ccs.PLATFORMS)
    assert tuple(by_id["gas"]["options"]) == ccs.GASES and tuple(by_id["tw_calibration"]["options"]) == ccs.TW_CALIBRATIONS
    assert by_id["adduct"]["options"] == spec.ADDUCTS_POS[1:] + spec.ADDUCTS_NEG[1:]
    assert by_id["observed"]["type"] == "text" and "comma-separated" in by_id["observed"]["label"]
    summary = next(b for b in mod["outputs"] if b.get("id") == "summary")
    assert (summary["style"], [i["label"] for i in summary["items"]]) == ("metric", ["Reference CCS", "Source", "Observed", "Δ", "Verdict"])
    assert all(i["help"] for i in summary["items"])
    assert client.post("/api/modules/rt/run", json={}).status_code == 404


def test_every_number_is_in_the_config():
    c = ccs.CONFIG
    assert (c.tol_same, c.tol_dt_tims, c.tol_dt_tw, c.tol_tims_tw) == (1.0, 1.0, 2.0, 3.0)
    assert (c.domain_min_tanimoto, c.min_adduct_rows, c.tol_round) == (0.5, 100, 0.1)
    assert {f.name for f in dataclasses.fields(c)} >= {"seed", "test_fraction", "val_fraction", "observed_min", "observed_max", "mz_check_da"}
    assert dataclasses.is_dataclass(c) and c.__dataclass_params__.frozen


@pytest.mark.parametrize("obs,ref,want", [("DT", "DT", 1.0), ("TW", "TW", 1.0), ("TIMS", "TIMS", 1.0), ("DT", "TIMS", 1.0), ("TIMS", "DT", 1.0),
                                          ("DT", "TW", 2.0), ("TW", "DT", 2.0), ("TIMS", "TW", 3.0), ("TW", "TIMS", 3.0)])
def test_tolerance_follows_the_platform_pair(obs, ref, want):
    assert ccs.tolerance(obs, ref) == want
    custom = dataclasses.replace(ccs.CONFIG, tol_same=0.5, tol_dt_tims=0.7, tol_dt_tw=1.5, tol_tims_tw=2.5)
    assert ccs.tolerance(obs, ref, custom) == {1.0: (0.7 if obs != ref else 0.5), 2.0: 1.5, 3.0: 2.5}[want]


def test_platform_labels_map_to_types():
    assert ccs.PLATFORMS == {"DTIMS (stepped field)": "DT", "DTIMS (single field)": "DT", "TIMS": "TIMS", "TWIMS": "TW", "FAIMS / DMS": None}


# ---------------------------------------------------------------- 加合物
@pytest.mark.parametrize("raw,want", [("[M+]+", "[M]+"), ("M+NH4]+", "[M+NH4]+"), ("[M+H]+*", "[M+H]+"), ("[M+Na]+*", "[M+Na]+"), ("[M+H20-H]-", "[M+H2O-H]-"),
                                      ("[M-H2O-H]-", "[M-H-H2O]-"), ("[M-H2O+H]+", "[M+H-H2O]+"), ("[M+H3C2O2]-", "[M+CH3COO]-"), ("[M+FA-H]-", "[M+HCOO]-"),
                                      (" [M+H]+ ", "[M+H]+"), ("[M+Na-H+Cl]-", "[M+Cl+Na-H]-")])
def test_adduct_names_of_the_source_map_to_mzdesk_names(raw, want):
    assert ccs.map_adduct(raw) == want


def test_spec_adducts_are_their_own_canonical_names_and_unknown_ones_stay():
    for a in spec.ADDUCTS_POS[1:] + spec.ADDUCTS_NEG[1:]:
        assert ccs.map_adduct(a) == a, a
        assert ccs.adduct_mz(100.0, a) is not None, a
    assert ccs.map_adduct("[M+OAcO]-") == "[M+OAcO]-"          # 认不出的保持原样
    assert ccs.adduct_mz(100.0, "[M+OAcO]-") is None and ccs.adduct_mz(100.0, "nonsense") is None


def test_adduct_mz_and_charge():
    m = 194.080376                                              # caffeine, monoisotopic
    for adduct, want in (("[M+H]+", 195.0877), ("[M+Na]+", 217.0696), ("[M-H]-", 193.0731), ("[M+HCOO]-", 239.0786), ("[M+NH4]+", 212.1142),
                         ("[2M+H]+", 389.1680), ("[M+H-H2O]+", 177.0771)):
        assert ccs.adduct_mz(m, adduct) == pytest.approx(want, abs=2e-4), adduct
    assert ccs.adduct_mz(m, "[M-2H]2-") == pytest.approx((m - 2 * 1.007276) / 2, abs=1e-4)
    assert [ccs.adduct_charge(a) for a in ("[M+H]+", "[M-H]-", "[M-2H]2-", "[M+2H]2+", "[M]+", "junk")] == [1, -1, -2, 2, 1, 0]


# ---------------------------------------------------------------- 输入
def test_parse_observed():
    assert ccs.parse_observed("") == [] and ccs.parse_observed(None) == [] and ccs.parse_observed("  ") == []
    assert ccs.parse_observed("152.3") == [152.3] and ccs.parse_observed(152.3) == [152.3]
    assert ccs.parse_observed("152.3, 160.1;170 171.5") == [152.3, 160.1, 170.0, 171.5]
    for bad, text in (("abc", "not a number"), ("152.3, x", "not a number"), ("nan", "finite"), ("0.85", "1/K0"), ("1500", "too large"), (",".join(["150"] * 11), "at most")):
        with pytest.raises(ValueError, match=text):
            ccs.parse_observed(bad)


def test_comparability_priority():
    assert ccs.comparability("DTIMS (single field)", "N₂", ccs.TW_CALIBRATIONS[0])[0]
    assert ccs.comparability("TWIMS", "N₂", ccs.TW_CALIBRATIONS[0])[0]
    ok, short, long = ccs.comparability("FAIMS / DMS", "He", ccs.TW_CALIBRATIONS[1])         # FAIMS 先于 He
    assert not ok and "no CCS" in short and "compensation voltage" in long
    ok, short, long = ccs.comparability("TIMS", "He", ccs.TW_CALIBRATIONS[1])
    assert not ok and short == "measured in He" and "N₂ only" in long
    ok, short, long = ccs.comparability("TWIMS", "N₂", ccs.TW_CALIBRATIONS[1])
    assert not ok and "N₂→He" in short and "He reference values" in long
    assert ccs.comparability("TIMS", "N₂", ccs.TW_CALIBRATIONS[1])[0]                       # 校准选项只对 TWIMS 有效
    with pytest.raises(ValueError):
        ccs.comparability("Magic", "N₂", "")


# ---------------------------------------------------------------- 建索引（假数据集）
def test_build_rows_counts_and_skips(tmp_path):
    datasets, cache = _datasets()
    rows, sources, stats = ccs_setup.build_rows(datasets, cache, log=lambda *_: None)
    assert stats["rows_in_source"] == 11 and stats["rows_kept"] == len(rows) == 7
    assert stats["rows_skipped"] == {"bad_value": 1, "no_structure": 1, "duplicate": 1, "bad_structure": 1}
    assert stats["rows_by_type"] == {"DT": 2, "TW": 4, "TIMS": 1} and stats["compounds"] == 3 and stats["datasets"] == 3
    assert stats["compounds_by_type"] == {"DT": 1, "TW": 3, "TIMS": 1}
    assert {r.adduct for r in rows} == {"[M+H]+", "[M+Na]+", "[M+H2O-H]-", "[M]+"}                   # c3sdb 的修正在建索引时就做了
    key = ci.normalize(CAFFEINE).key14
    assert {r.key for r in rows if r.name == "Caffeine"} == {key}
    cached = next(r for r in rows if r.name == "Cached acid")
    assert cached.key == ci.normalize(DECANOIC).key14 and cached.smiles == ci.normalize(DECANOIC).smiles
    assert sources["fakeDT"] == {"reference": "Fake, A. et al. fakeDT (2020).", "url": "https://doi.org/10.1000/fakeDT", "ccs_type": "DT",
                                 "ccs_method": "single field, calibrated", "n": 2}
    assert {r.z for r in rows} == {1}


def test_index_round_trip_and_old_format(tmp_path):
    make_index(tmp_path)
    index = ccs.load(tmp_path)
    assert len(index) == 7 and index.meta["commit"] == "f" * 40 and ccs.installed(tmp_path)
    assert ccs.load(tmp_path) is index                                                  # 缓存
    key = ci.normalize(CAFFEINE).key14
    assert index.compound_name(key) == "Caffeine"
    same, others = index.lookup(key, "[M+H]+")
    assert [(r.ccs_type, r.src_tag, r.ccs) for r in same] == [("DT", "fakeDT", 140.5), ("TIMS", "fakeTIMS", 141.0), ("TW", "fakeTW", 138.3)]     # DT, TIMS, TW
    assert [(r.adduct, r.ccs) for r in others] == [("[M+Na]+", 142.9)]
    assert ccs.load(tmp_path / "nowhere") is None and not ccs.installed(tmp_path / "nowhere")
    with gzip.open(tmp_path / ccs.INDEX_FILE, "rt", encoding="utf-8") as f:
        payload = json.load(f)
    payload["schema"] = 99
    with gzip.open(tmp_path / ccs.INDEX_FILE, "wt", encoding="utf-8") as f:
        json.dump(payload, f)
    ccs._cache.clear()
    with pytest.raises(ccs.IndexProblem, match="old format"):
        ccs.load(tmp_path)


# ---------------------------------------------------------------- 比较
def test_only_the_same_adduct_is_compared(fake_ccs):
    res = ccs.analyse(CAFFEINE, "[M+H]+", [142.9], "DTIMS (single field)")
    assert {r.adduct for r in res.rows} == {"[M+H]+"} and len(res.rows) == 3
    assert [(r.adduct, r.ccs) for r in res.others] == [("[M+Na]+", 142.9)]                # 另一个加合物只作参考，哪怕数值刚好吻合
    assert res.verdict == ccs.INCONSISTENT and res.used.row.src_tag == "fakeTIMS"           # +1.71 % / +1.34 % / +3.3 % -> 比容差 1.71 / 1.34 / 1.65 倍：TIMS 行最近
    na = ccs.analyse(CAFFEINE, "[M+Na]+", [142.9], "DTIMS (single field)")
    assert len(na.rows) == 1 and na.verdict == ccs.CONSISTENT and na.used.delta_pct == pytest.approx(0)
    assert [r.adduct for r in na.others] == ["[M+H]+"] * 3
    assert res.name == "Caffeine" and res.key == ci.normalize(CAFFEINE).key14


def test_a_salt_or_charged_form_is_the_same_compound(fake_ccs):
    res = ccs.analyse(CAFFEINE + ".[Na+].[Cl-]", "[M+H]+", [140.5], "DTIMS (single field)")
    assert len(res.rows) == 3 and res.verdict == ccs.CONSISTENT and any("largest fragment" in n for n in res.notes)


def test_multiple_observed_values_closest_wins(fake_ccs):
    res = ccs.analyse(CAFFEINE, "[M+H]+", [160.0, 140.2], "DTIMS (single field)")
    dt = next(c for c in res.comparisons if c.ref_type == "DT")
    assert dt.observed == 140.2 and dt.delta_pct == pytest.approx((140.2 - 140.5) / 140.5 * 100) and dt.verdict == ccs.CONSISTENT
    assert res.verdict == ccs.CONSISTENT and res.used.observed == 140.2
    far = ccs.analyse(CAFFEINE, "[M+H]+", [160.0, 170.0], "DTIMS (single field)")
    assert all(c.observed == 160.0 for c in far.comparisons) and far.verdict == ccs.INCONSISTENT
    # 一个参考值只和离它最近的观测值比：这里 140.2 对 TW 行 138.3 是 +1.37 %，对 DT 容差 1 % 一致、对 TW 容差 2 % 也一致
    assert [c.verdict for c in res.comparisons] == [ccs.CONSISTENT] * 3


def test_tolerance_depends_on_the_platform_pair(fake_ccs):
    dt_row, tims_row, tw_row = [(c, c.row) for c in ccs.analyse(CAFFEINE, "[M+H]+", [], "TIMS").comparisons]
    assert [c.row.ccs_type for c, _ in (dt_row, tims_row, tw_row)] == ["DT", "TIMS", "TW"]
    obs = 143.0                                                                           # DT 140.5: +1.78 %, TIMS 141.0: +1.42 %, TW 138.3: +3.40 %
    by_platform = {}
    for platform in ("DTIMS (stepped field)", "TIMS", "TWIMS"):
        res = ccs.analyse(CAFFEINE, "[M+H]+", [obs], platform)
        by_platform[platform] = {c.ref_type: (c.tolerance, c.verdict) for c in res.comparisons}
    assert by_platform["DTIMS (stepped field)"] == {"DT": (1.0, ccs.INCONSISTENT), "TIMS": (1.0, ccs.INCONSISTENT), "TW": (2.0, ccs.INCONSISTENT)}
    assert by_platform["TIMS"] == {"DT": (1.0, ccs.INCONSISTENT), "TIMS": (1.0, ccs.INCONSISTENT), "TW": (3.0, ccs.INCONSISTENT)}
    assert by_platform["TWIMS"] == {"DT": (2.0, ccs.CONSISTENT), "TIMS": (3.0, ccs.CONSISTENT), "TW": (1.0, ccs.INCONSISTENT)}
    assert ccs.analyse(CAFFEINE, "[M+H]+", [139.0], "TIMS").comparisons[2].verdict == ccs.CONSISTENT       # TIMS vs TW: +0.5 % (3 %)
    assert ccs.analyse(CAFFEINE, "[M+H]+", [138.3 * 1.029], "TIMS").comparisons[2].verdict == ccs.CONSISTENT
    assert ccs.analyse(CAFFEINE, "[M+H]+", [138.3 * 1.031], "TIMS").comparisons[2].verdict == ccs.INCONSISTENT


def test_the_boundary_is_inclusive(fake_ccs):
    res = ccs.analyse(CAFFEINE, "[M+H]+", [140.5 * 1.01], "DTIMS (single field)")
    assert res.comparisons[0].delta_pct == pytest.approx(1.0) and res.comparisons[0].verdict == ccs.CONSISTENT


@pytest.mark.parametrize("platform,gas,tw", [("DTIMS (single field)", "He", "N₂ reference values"), ("TIMS", "He", "N₂ reference values"),
                                              ("TWIMS", "N₂", "He reference values (N₂→He)"), ("TWIMS", "He", "N₂ reference values"),
                                              ("FAIMS / DMS", "N₂", "N₂ reference values")])
def test_not_comparable_cases_compare_nothing_but_still_list_references(fake_ccs, platform, gas, tw):
    res = ccs.analyse(CAFFEINE, "[M+H]+", [140.5], platform, gas, tw)
    assert res.verdict == ccs.NOT_COMPARABLE and not res.comparable and res.reason
    assert len(res.rows) == 3 and all(c.verdict == ccs.NOT_COMPARABLE and c.delta_pct is None and c.tolerance is None for c in res.comparisons)
    assert res.used is None and res.ref_ccs == pytest.approx(140.5)                       # 参考值仍然给出（中位数），只是不比较
    body = _run(observed="140.5", platform=platform, gas=gas, tw=tw).json()
    assert _card(body, "Verdict")["value"] == "Not comparable" and _card(body, "Δ")["value"] == "—"
    rows = body["outputs"]["reference"]["rows"]
    assert len(rows) == 3 and all(_cells(r)[4:] == ["—", "—", {"badge": "Not comparable", "tone": "neutral"}] for r in rows)
    assert body["verdict"].startswith("<b>Not comparable</b>. ") and _card(body, "Verdict")["help"] in body["verdict"]


def test_faims_says_so_plainly(fake_ccs):
    body = _run(observed="140.5", platform="FAIMS / DMS").json()
    verdict = _card(body, "Verdict")
    assert verdict["note"] == "FAIMS / DMS gives no CCS" and "compensation voltage, not a CCS" in verdict["help"]
    assert body["verdict"] == "<b>Not comparable</b>. FAIMS / DMS gives a compensation voltage, not a CCS, so nothing is compared."


def test_no_observed_value_gives_the_median_reference(fake_ccs):
    res = ccs.analyse(CAFFEINE, "[M+H]+", [], "TIMS")
    assert res.verdict == ccs.NO_OBSERVED and res.ref_ccs == 140.5 and res.used is None            # 138.3, 140.5, 141.0 的中位数
    assert [c.verdict for c in res.comparisons] == [ccs.NO_OBSERVED] * 3 and [c.tolerance for c in res.comparisons] == [1.0, 1.0, 3.0]
    body = _run().json()
    assert _card(body, "Verdict")["value"] == "No observed value" and _card(body, "Reference CCS") == {"value": "140.5", "note": "median of 3 values (138.3–141.0)"}
    assert body["verdict"] == "No observed CCS entered. N₂ reference CCS for [M+H]+: <b>140.5 Å²</b> (median of 3 values (138.3–141.0))."


def test_no_reference_and_no_model(fake_ccs):
    res = ccs.analyse(IBUPROFEN, "[M+H]+", [145.0], "TIMS")
    assert res.verdict == ccs.NO_REFERENCE and res.rows == [] and res.prediction is None and "No prediction model is installed" in res.prediction_note
    body = _run(smiles=IBUPROFEN, observed="145").json()
    assert _card(body, "Verdict")["value"] == "No reference" and body["outputs"]["reference"]["rows"] == [] and body["outputs"]["prediction"]["rows"] == []
    assert body["verdict"].startswith("No reference CCS for [M+H]+ of this compound in the data. No prediction model is installed.")
    assert any("No prediction model" in n for n in body["outputs"]["notes"]["items"])


# ---------------------------------------------------------------- 预测后备（假模型）
def test_prediction_is_used_only_without_a_reference_value(fake_model):
    ref = ccs.analyse(CAFFEINE, "[M+H]+", [140.5], "DTIMS (single field)")
    assert ref.source == "measured" and ref.prediction is None and not fake_model.seen              # 有参考值：模型根本不被调用
    res = ccs.analyse(CAFFEINE, "[M+K]+", [151.0], "DTIMS (single field)")                           # 咖啡因有数据，但 [M+K]+ 没有
    assert res.rows == [] and res.source == "predicted" and res.ref_ccs == 150.0
    p = res.prediction
    assert (p.ccs, p.tolerance_pct, p.median_pct, p.q95_pct, p.n_test, p.model) == (150.0, 5.6, 1.2, 5.51, 150, "SVR")
    assert p.in_domain and p.max_tanimoto == pytest.approx(1.0)                                     # 咖啡因在训练集里
    c = res.prediction_comparison
    assert c.tolerance == 5.6 and c.delta_pct == pytest.approx((151.0 - 150.0) / 150.0 * 100) and res.verdict == ccs.CONSISTENT
    x = fake_model.seen[-1]
    assert x.shape == (1, 1 + len(MODEL_INFO["adducts"]) + 42)                                      # m/z、加合物 one-hot、42 个 MQN
    assert x[0, 0] == pytest.approx(ccs.adduct_mz(194.080376, "[M+K]+"), abs=1e-3) and list(x[0, 1:4]) == [0.0, 0.0, 1.0]


def test_prediction_tolerance_is_the_test_q95(fake_model):
    inside = ccs.analyse(CAFFEINE, "[M+K]+", [150.0 * 1.055], "TIMS")
    outside = ccs.analyse(CAFFEINE, "[M+K]+", [150.0 * 1.057], "TIMS")
    assert inside.verdict == ccs.CONSISTENT and outside.verdict == ccs.INCONSISTENT


def test_prediction_flags_the_training_domain(fake_model):
    res = ccs.analyse(DECANOIC, "[M+H]+", [], "TIMS")                                               # 和咖啡因、阿司匹林都不像
    assert res.prediction.max_tanimoto < ccs.CONFIG.domain_min_tanimoto and not res.prediction.in_domain
    body = _run(smiles=DECANOIC).json()
    assert any("Outside the model's training domain" in w for w in body["warnings"])
    row = _cells(body["outputs"]["prediction"]["rows"][0])
    assert row[0] == "150.0" and row[1] == "±5.6" and row[2:4] == ["1.20", "5.51"] and row[4]["badge"] == "Outside" and row[4]["tone"] == "warn"
    assert "outside its training domain" in _card(body, "Source")["note"]


def test_prediction_is_labelled_as_a_prediction(fake_model):
    body = _run(smiles=CAFFEINE, adduct="[M+K]+", observed="149").json()
    assert _card(body, "Source") == {"value": "Predicted", "note": "SVR"}
    assert _card(body, "Reference CCS") == {"value": "150.0", "note": "predicted, ±5.6 %"}
    verdict = _card(body, "Verdict")
    assert verdict == {"value": "Consistent", "note": "vs prediction (weaker)"}
    assert "predicted N₂ value" in body["verdict"] and "A prediction is weaker than a measured value." in body["verdict"]
    assert body["outputs"]["reference"]["rows"] == []                                              # 没有实测行
    row = body["outputs"]["prediction"]["rows"][0]
    assert _cells(row)[4]["badge"] == "Inside" and any("Test error: 150 values of 80 compounds" in h for h in row["help"])
    assert not body["warnings"]


def test_no_prediction_for_a_rare_adduct(fake_model):
    res = ccs.analyse(CAFFEINE, "[M+HCOO]-", [150.0], "TIMS")                                       # 模型只有 [M+H]+ / [M+Na]+ / [M+K]+
    assert res.prediction is None and res.verdict == ccs.NO_REFERENCE and "No prediction for this adduct" in res.prediction_note
    body = _run(adduct="[M+HCOO]-", observed="150").json()
    assert _card(body, "Verdict")["note"].startswith("No prediction for this adduct")
    assert not fake_model.seen


def test_not_comparable_still_shows_the_prediction(fake_model):
    res = ccs.analyse(CAFFEINE, "[M+K]+", [150.0], "TIMS", gas="He")
    assert res.verdict == ccs.NOT_COMPARABLE and res.prediction is not None and res.prediction_comparison.verdict == ccs.NOT_COMPARABLE


# ---------------------------------------------------------------- 页面（runner）
def test_page_outputs_for_a_compound_with_references(fake_ccs):
    body = _run(observed="140.2").json()
    assert body["status"] == "ok" and set(body["outputs"]) == BLOCKS and body["warnings"] == []
    out = body["outputs"]
    assert out["summary"]["values"] == {
        "Reference CCS": {"value": "140.5", "note": "closest of 3 measured values"},
        "Source": {"value": "Measured", "note": "DTCCS N2 · fakeDT"},
        "Observed": {"value": "140.2", "note": "DTIMS (single field)"},
        "Δ": {"value": "−0.21", "note": "tolerance ±1.0 %"},
        "Verdict": {"value": "Consistent", "note": "vs measured reference"}}
    ref = out["reference"]
    assert ref["count"] == 3 and len(ref["rows"]) == 3
    first = ref["rows"][0]
    assert _cells(first) == [{"text": "fakeDT", "url": "https://doi.org/10.1000/fakeDT"}, "DTCCS N2", "single field, calibrated", "140.5", "−0.21", "±1.0",
                             {"badge": "Consistent", "tone": "good"}]
    assert first["help"][0] == "Reference: Fake, A. et al. fakeDT (2020)." and any("Closest observed value: 140.2" in h for h in first["help"])
    assert [_cells(r)[1] for r in ref["rows"]] == ["DTCCS N2", "TIMSCCS N2", "TWCCS N2"]
    assert [_cells(r)[6]["tone"] for r in ref["rows"]] == ["good", "good", "good"]
    assert [_cells(r)[5] for r in ref["rows"]] == ["±1.0", "±1.0", "±2.0"]
    assert out["other_ions"] == {"rows": [["[M+Na]+", "142.9", "DTCCS N2", "fakeDT"]], "count": 1}
    assert out["prediction"] == {"rows": [], "count": None} or out["prediction"]["rows"] == []
    notes = " ".join(out["notes"]["items"])
    for text in ("same adduct", "N₂ scale only", "TWCCS N₂→He", "compensation voltage", "not 1/K0", "same type 1.0 %", "DT and TW 2.0 %", "TIMS and TW 3.0 %",
                 "Stow et al.", "Feuerstein et al.", "never used to exclude", "commit ffffffffff", "MIT licence", "doi.org/10.1021/acs.analchem.9b05772"):
        assert text in notes, text
    assert body["verdict"] == "CCS <b>140.2 Å²</b> is <b>consistent</b> with the measured N₂ reference 140.5 Å² (Δ −0.21 %, tolerance ±1.0 %)."


def test_inconsistent_never_excludes(fake_ccs):
    body = _run(observed="150").json()
    assert _card(body, "Verdict")["value"] == "Inconsistent" and "does not exclude the compound" in _card(body, "Verdict")["help"]
    assert [_cells(r)[6] for r in body["outputs"]["reference"]["rows"]] == [{"badge": "Inconsistent", "tone": "warn"}] * 3
    assert body["verdict"].endswith("This does not exclude the compound.") and "<b>inconsistent</b>" in body["verdict"]


def test_page_with_several_observed_values_shows_the_closest(fake_ccs):
    body = _run(observed="160.1, 140.2").json()
    assert _card(body, "Observed") == {"value": "160.1, 140.2", "note": "DTIMS (single field)"}
    assert _card(body, "Δ")["note"] == "tolerance ±1.0 %; closest 140.2 Å²"
    assert "CCS <b>140.2 Å²</b> is <b>consistent</b>" in body["verdict"]


def test_page_notes_the_mz_that_does_not_match(fake_ccs):
    # 假数据里 Caffeine 的 [M+H]+ m/z 都对；把 TW 行的 m/z 改错，行的 ? 要提示
    index = ccs.load()
    bad = dataclasses.replace(next(r for r in index.rows if r.src_tag == "fakeTW" and r.name == "Caffeine"), mz=196.2)
    assert ccs.mz_mismatch(bad) == (196.2, pytest.approx(195.0877, abs=1e-3))
    assert ccs.mz_mismatch(next(r for r in index.rows if r.src_tag == "fakeDT")) is None
    ccs.save_index(fake_ccs / ccs.INDEX_FILE, [bad if r.src_tag == "fakeTW" and r.name == "Caffeine" else r for r in index.rows], index.sources)
    ccs._cache.clear()
    rows = _run(observed="140.2").json()["outputs"]["reference"]["rows"]
    tw = next(r for r in rows if _cells(r)[1] == "TWCCS N2")
    assert any(h.startswith("Check this row: the listed m/z (196.2000) differs from the m/z calculated") for h in tw["help"])
    assert not any(h.startswith("Check this row") for r in rows if r is not tw for h in r["help"])


def test_page_runs_without_the_index_and_says_how_to_install_it(tmp_path, monkeypatch):
    monkeypatch.setenv("MZDESK_CCS_DIR", str(tmp_path / "empty"))
    body = _run(observed="140.2").json()
    assert body["status"] == "ok" and set(body["outputs"]) == BLOCKS
    assert _card(body, "Source") == {"value": "Not installed", "note": "run python -m mzdesk.orthogonal.ccs_setup"}
    assert _card(body, "Verdict") == {"value": "No reference", "note": "CCS data not installed"}
    assert _card(body, "Observed")["value"] == "140.2"
    assert body["verdict"] == "CCS reference data is not installed: run <code>python -m mzdesk.orthogonal.ccs_setup</code>."
    assert any("python -m mzdesk.orthogonal.ccs_setup" in w for w in body["warnings"])
    assert "ccs_setup" in body["outputs"]["notes"]["items"][0] and "Comparison rule" in " ".join(body["outputs"]["notes"]["items"])
    assert all(body["outputs"][b]["rows"] == [] for b in ("reference", "other_ions", "prediction"))


def test_page_errors(fake_ccs):
    assert _run(observed="abc").status_code == 400
    r = _run(observed="0.85")
    assert r.status_code == 400 and "1/K0" in r.json()["detail"]
    assert _run(smiles="not a smiles").status_code == 400
    r = _run(platform="Magic")
    assert r.status_code == 400 and "Platform must be one of" in r.json()["detail"]
    assert _run(gas="Ar").status_code == 400 and _run(tw="x").status_code == 400
    assert _run(adduct="[M+Zz]+").status_code == 400
    assert client.post("/api/modules/ccs/run", json={"structure": {"mode": "SMILES / InChI", "value": ""}}).status_code == 400


def test_old_index_format_is_a_server_side_problem(fake_ccs):
    with gzip.open(fake_ccs / ccs.INDEX_FILE, "wt", encoding="utf-8") as f:
        json.dump({"schema": 0, "columns": [], "rows": [], "sources": {}}, f)
    ccs._cache.clear()
    r = _run()
    assert r.status_code == 502 and "old format" in r.json()["detail"]


def test_twims_calibration_is_ignored_for_other_platforms(fake_ccs):
    body = _run(observed="140.2", platform="TIMS", tw="He reference values (N₂→He)").json()
    assert _card(body, "Verdict")["value"] == "Consistent" and any("only applies to TWIMS" in w for w in body["warnings"])


def test_defaults_when_the_optional_fields_are_missing(fake_ccs):
    body = client.post("/api/modules/ccs/run", json={"structure": {"mode": "SMILES / InChI", "value": CAFFEINE}}).json()
    assert _card(body, "Verdict")["value"] == "No observed value" and body["outputs"]["reference"]["count"] == 3


# ---------------------------------------------------------------- 一句话结论
def test_verdict_function_reads_the_cards():
    def cards(**v):
        return {"summary": {"values": v}}
    assert verdicts.build("ccs", {}, {"summary": "not a dict"}) is None
    assert verdicts.build("ccs", {}, cards(Verdict={"value": "—"})) is None
    base = {"Reference CCS": {"value": "151.4", "note": "closest of 3 measured values"}, "Source": {"value": "Measured", "note": "DTCCS N2 · x"}, "Observed": "152.3",
            "Δ": {"value": "+0.59", "note": "tolerance ±1.0 %"}}
    assert verdicts.build("ccs", {}, cards(**base, Verdict={"value": "Consistent"})) == \
        "CCS <b>152.3 Å²</b> is <b>consistent</b> with the measured N₂ reference 151.4 Å² (Δ +0.59 %, tolerance ±1.0 %)."
    out = verdicts.build("ccs", {}, cards(**{**base, "Δ": {"value": "-3.2", "note": "tolerance ±1.0 %"}}, Verdict={"value": "Inconsistent"}))
    assert "<b>inconsistent</b>" in out and "Δ −3.2 %" in out and out.endswith("This does not exclude the compound.")
    out = verdicts.build("ccs", {"adduct": "[M-H]-"}, cards(**{**base, "Verdict": {"value": "No reference", "note": "No prediction for this adduct (fewer than 100 training rows)."},
                                                              "Reference CCS": "—", "Source": "—"}))
    assert out == "No reference CCS for [M-H]- of this compound in the data. No prediction for this adduct (fewer than 100 training rows)."
    assert verdicts.build("ccs", {}, cards(**base, Verdict={"value": "Not comparable", "note": "measured in He"})) == "<b>Not comparable</b>. measured in He"
    assert verdicts.build("ccs", {}, cards(**base, Verdict={"value": "Not comparable", "note": "measured in He", "help": "A CCS measured in He is not comparable."})) == \
        "<b>Not comparable</b>. A CCS measured in He is not comparable."
    out = verdicts.build("ccs", {"adduct": "[M+Na]+"}, cards(**{**base, "Source": {"value": "Predicted", "note": "SVR"}}, Verdict={"value": "No observed value"}))
    assert out.startswith("No observed CCS entered. No reference value in the data; predicted CCS for [M+Na]+: <b>151.4 Å²</b>")


# ---------------------------------------------------------------- 训练（合成小数据，需要 scikit-learn）
def test_training_end_to_end_on_synthetic_rows(tmp_path, monkeypatch):
    pytest.importorskip("sklearn")
    import joblib   # noqa: F401
    rows = []
    for n in range(2, 62):
        for smiles in ("C" * n + "O", "C" * n + "N", "C" * (n // 2) + "C(C)C" + "C" * (n - n // 2) + "O"):
            norm = ci.normalize(smiles)
            mass = 12.0 * n + 20
            for adduct, shift in (("[M+H]+", 1.0073), ("[M+Na]+", 22.9892)):
                value = 60 + 1.9 * mass + (5 if adduct == "[M+Na]+" else 0) + (n % 3)
                rows.append(ccs.Row(norm.key14, f"c{n}", norm.smiles, adduct, mass + shift, 1, value, "DT" if n % 2 else "TW", "m", "t"))
    keys = {r.key for r in rows}
    cfg = dataclasses.replace(ccs.CONFIG, min_adduct_rows=40)
    lines = []
    grid = [("SVR", {"C": 10.0, "gamma": 0.01, "epsilon": 0.02}), ("HistGradientBoosting", {"learning_rate": 0.1, "max_iter": 50, "max_leaf_nodes": 8})]
    info = ccs_setup.train(rows, tmp_path, cfg, log=lines.append, grid=grid)
    assert [v["name"] for v in info["validation"]] == ["SVR", "HistGradientBoosting"]
    assert (tmp_path / ccs.MODEL_FILE).is_file() and (tmp_path / ccs.DOMAIN_FILE).is_file()
    assert info["name"] in ("SVR", "HistGradientBoosting") and list(info["adducts"]) == ["[M+H]+", "[M+Na]+"] or sorted(info["adducts"]) == ["[M+H]+", "[M+Na]+"]
    split = info["split"]
    assert split["by"] == "connectivity key" and split["test_compounds"] > 0 and split["trainval_compounds"] + split["test_compounds"] <= len(keys)
    test = info["test"]
    assert test["n"] == split["test_rows"] > 0 and 0 < test["median_pct"] <= test["q95_pct"] and set(test["by_type"]) == {"DT", "TW", "TIMS"}
    assert info["tolerance_pct"] >= test["q95_pct"] and info["tolerance_pct"] - test["q95_pct"] < 0.1 + 1e-9      # 向上取整到 0.1 %
    assert round(info["tolerance_pct"] * 10) == pytest.approx(info["tolerance_pct"] * 10)
    assert info["seed"] == cfg.seed and info["domain"]["threshold"] == cfg.domain_min_tanimoto and len(info["domain"]["test_bins"]) == 6
    assert any("Split by compound" in line for line in lines)
    # 训练集和测试集没有同一个化合物：重新按同样的种子划分
    from sklearn.model_selection import GroupShuffleSplit
    kk = np.array([r.key for r in rows])
    tv, te = next(GroupShuffleSplit(n_splits=1, test_size=cfg.test_fraction, random_state=cfg.seed).split(np.zeros(len(rows)), groups=kk))
    assert not set(kk[tv]) & set(kk[te])
    # 装好以后能预测
    index_dir = tmp_path
    ccs.save_index(index_dir / ccs.INDEX_FILE, rows, {})
    ccs_setup.write_meta(index_dir, {"commit": "x", "model": info})
    monkeypatch.setenv("MZDESK_CCS_DIR", str(index_dir))
    predictor, why = ccs.load_predictor(index_dir)
    assert predictor is not None and why == ""
    p = predictor.predict(ci.normalize("C" * 20 + "O"), "[M+H]+", 12.0 * 20 + 20 + 1.0073)
    assert 60 + 1.9 * 260 - 40 < p.ccs < 60 + 1.9 * 260 + 40 and p.in_domain and p.tolerance_pct == info["tolerance_pct"]
    assert predictor.predict(ci.normalize("C" * 20 + "O"), "[M-H]-", 259.0) is None             # 训练集里没有的加合物
    info2 = dict(info, sklearn="0.0.1")
    ccs_setup.write_meta(index_dir, {"commit": "x", "model": info2})
    ccs._models.clear()
    none, why = ccs.load_predictor(index_dir)
    assert none is None and "scikit-learn 0.0.1" in why


def test_rounding_up_to_a_tenth():
    assert ccs_setup.round_up(6.55, 0.1) == 6.6 and ccs_setup.round_up(6.5, 0.1) == 6.5 and ccs_setup.round_up(6.51, 0.1) == 6.6
    assert ccs_setup.error_summary([]) == {"n": 0, "median_pct": None, "q95_pct": None}
    assert ccs_setup.error_summary([0.01, 0.02, 0.03])["median_pct"] == 2.0


@pytest.mark.skipif(shutil.which("git") is None, reason="git not installed")
def test_fetch_checks_out_only_the_data_files_of_the_pinned_commit(tmp_path, monkeypatch, capsys):
    """用本地假上游仓库试稀疏克隆：只检出 LICENSE / README / json，克隆里的 Python 文件不会出现；读 json 不执行任何东西。"""
    up = tmp_path / "upstream"
    (up / "c3sdb" / "_include" / "src_data").mkdir(parents=True)
    (up / "c3sdb" / "_include" / "smiles_search_cache.json").write_text(json.dumps({"Cached acid": DECANOIC}))
    datasets, _ = _datasets()
    (up / "c3sdb" / "_include" / "src_data" / "fakeDT.json").write_text(json.dumps({"metadata": datasets[0][0], "data": datasets[0][1]}))
    (up / "c3sdb" / "__init__.py").write_text("raise SystemExit('planted code must never run')\n")
    (up / "LICENSE.txt").write_text("MIT License\nCopyright 2020 Dylan H. Ross\n")
    (up / "README.md").write_text("readme\n")
    (up / "setup.cfg").write_text("[metadata]\n")
    git = lambda *a: subprocess.run(["git", "-C", str(up), "-c", "user.name=t", "-c", "user.email=t@t", *a], check=True, capture_output=True, text=True).stdout.strip()   # noqa: E731
    git("init", "-q")
    git("config", "uploadpack.allowFilter", "true")
    git("config", "uploadpack.allowAnySHA1InWant", "true")
    git("add", "-A")
    git("commit", "-q", "-m", "init")
    sha = git("rev-parse", "HEAD")
    monkeypatch.setattr(ccs_setup, "REPO_GIT", str(up))
    repo = tmp_path / "clone"
    ccs_setup.fetch(repo, sha)
    files = sorted(p.relative_to(repo).as_posix() for p in repo.rglob("*") if p.is_file() and ".git" not in p.parts)
    assert files == ["LICENSE.txt", "README.md", "c3sdb/_include/smiles_search_cache.json", "c3sdb/_include/src_data/fakeDT.json"]
    ccs_setup.fetch(repo, sha)
    assert "already downloaded" in capsys.readouterr().out
    loaded, cache = ccs_setup.read_source(repo)
    assert [m["src_tag"] for m, _ in loaded] == ["fakeDT"] and cache == {"Cached acid": DECANOIC}
    with pytest.raises(SystemExit):
        ccs_setup.fetch(tmp_path / "clone2", "0" * 40)                                          # 不存在的提交


def test_setup_pins_the_commit_and_the_files():
    assert ccs_setup.COMMIT == "5f1393eec6f0acea9af9c592b911d5602a9236ea" and ccs_setup.REPO_GIT == "https://github.com/dylanhross/c3sdb.git"
    assert "/c3sdb/_include/src_data/*.json" in ccs_setup.SPARSE and "/LICENSE.txt" in ccs_setup.SPARSE
    assert ccs.SETUP_COMMAND == "python -m mzdesk.orthogonal.ccs_setup"


def test_ccs_dir_follows_the_environment(monkeypatch, tmp_path):
    monkeypatch.setenv("MZDESK_CCS_DIR", str(tmp_path / "x"))
    assert ccs.ccs_dir() == tmp_path / "x"
    monkeypatch.delenv("MZDESK_CCS_DIR")
    assert ccs.ccs_dir().parts[-3:] == (".mzdesk", "data", "ccs")


# ---------------------------------------------------------------- 真实索引（没安装时跳过）
@needs_real_ccs
def test_real_index_counts_and_licence():
    index = ccs.load()
    meta = index.meta
    assert meta["commit"] == ccs_setup.COMMIT and meta["licence"] == "MIT" and meta["datasets"] == 26
    assert meta["rows_kept"] == len(index) == sum(meta["rows_by_type"].values()) > 14000
    assert (ccs.ccs_dir() / "LICENSE.txt").is_file() and "Dylan" in (ccs.ccs_dir() / "LICENSE.txt").read_text(encoding="utf-8")
    assert {r.ccs_type for r in index.rows} == {"DT", "TW", "TIMS"} and all(r.ccs > 0 for r in index.rows)


@needs_real_ccs
def test_real_example_runs():
    mod = next(m for m in MODULES if m["id"] == "ccs")
    body = client.post("/api/modules/ccs/run", json=mod["example"]).json()
    assert body["status"] == "ok" and set(body["outputs"]) == BLOCKS
    assert _card(body, "Verdict")["value"] == "Consistent" and _card(body, "Source")["value"] == "Measured"
    assert len(body["outputs"]["reference"]["rows"]) >= 5 and body["verdict"].startswith("CCS <b>140.4 Å²</b> is <b>consistent</b>")


@needs_real_ccs
def test_real_model_predicts_a_compound_that_is_not_in_the_data():
    model, why = ccs.load_predictor()
    if model is None:
        pytest.skip(why)
    res = ccs.analyse("CCCCCCCCCCCCCC(=O)NCCO", "[M+H]+", [], "TIMS")               # N-myristoylethanolamine: a lipid, not in the data
    if res.rows:
        pytest.skip("the compound is in the index")
    assert res.source == "predicted" and 150 < res.prediction.ccs < 200
    assert res.prediction.tolerance_pct == res.meta["model"]["tolerance_pct"] > res.meta["model"]["test"]["q95_pct"] - 1e-9
