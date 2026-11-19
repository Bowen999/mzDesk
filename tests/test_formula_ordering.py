"""Regression checks for formula order on the pages that share the annotation engines."""
from types import SimpleNamespace

from mzdesk.api.runners import ms1, ms2_denovo, ms2_structure
from mzdesk.ms1 import msbuddy_engine as bd
from mzdesk.ms2 import msfinder_engine as mf


def test_buddy_formula_page_and_denovo_fallback_keep_engine_order(monkeypatch):
    # The scores round to the same displayed value; neither page may re-sort on that text.
    cands = [{"rank": i, "formula": formula, "prob": prob, "mz_theo": 195.0, "ppm": 0.0,
              "fdr": None, "isotope_sim": None}
             for i, (formula, prob) in enumerate([("C8H10N4O2", 0.90009), ("C9H14N2O3", 0.90008)], 1)]
    monkeypatch.setattr(bd, "annotate_formula", lambda *a, **k: {"candidates": cands, "adduct_used": "[M+H]+",
                                                                "from_input": False, "ms1_used": 0})
    inputs = {"precursor": {"precursor_mz": "195.0877", "ion_mode": "Positive", "adduct": "[M+H]+"}}
    outputs = ms1.run_formula(inputs)["outputs"]
    rows = outputs["candidates"]["rows"]
    assert [(row[0], row[1], row[4]) for row in rows] == [("1", "C8H10N4O2", "0.9001"), ("2", "C9H14N2O3", "0.9001")]
    assert outputs["summary"]["values"]["Top formula"]["value"] == rows[0][1]
    warnings = []
    formula, note = ms2_denovo._formula_from_buddy(SimpleNamespace(check_formula=lambda f: f), 195.0877, "[M+H]+", [], [], False, warnings)
    assert formula == rows[0][1] and "top-1" in note and formula in warnings[0]


def test_msfinder_formula_page_uses_numeric_full_precision_scores_and_missing_last():
    # Scores 10.00001 and 10.00000 display alike; 10 must not sort below 9 as a string.
    text = "Formula\tTheoretical mass\tMass error\tFormula score\n" + "\n".join([
        "C1\t100\t0\t1.9", "C2\t100\t0\t10.00000", "C3\t100\t0\t9.99996",
        "C4\t100\t0\t10.00001", "C5\t100\t0\t", "C6\t100\t0\t10.00000"])
    parsed = mf.parse_formula_results(text)["rows"]
    rows = ms2_structure._formula_rows(parsed)
    assert [r["formula"] for r in parsed] == ["C4", "C2", "C6", "C3", "C1", "C5"]
    assert [(r[0], r[1]) for r in rows] == [(str(i), f["formula"]) for i, f in enumerate(parsed, 1)]
    assert [r[4] for r in rows[:4]] == ["10.000"] * 4
    assert rows[-1][4] == "—"
