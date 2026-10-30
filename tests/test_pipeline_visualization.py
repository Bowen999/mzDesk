"""Workflow display data are tied to the input and the best judged library hit."""
import json

import numpy as np

from mzdesk.annotation import library_step, visualization


def test_reference_is_the_best_judged_hit_not_the_first_similarity_hit(monkeypatch):
    def hit(row, adduct, similarity):
        return {"row": row, "library": "MassBank", "accession": f"hit-{row}", "name": f"hit {row}",
                "adduct": adduct, "precursor_mz": 195.0877, "ppm": 0.0, "matches": 4,
                "scores": {library_step.PRIMARY_METRIC: similarity}, "instrument": "QTOF", "collision_energy": "30",
                "inchikey": "", "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"}

    hits = [hit(3, "unknown", 0.99), hit(7, "[M+H]+", 0.93)]
    search = {"hits": hits, "n_candidates": 2, "n_kept": 2, "warnings": [], "notes": {}}
    monkeypatch.setattr(library_step, "_search", lambda *a: (search, library_step.PRIMARY_METRIC))

    class Index:
        def peaks(self, row):
            return np.array([float(row), 69.0]), np.array([1.0, 0.5])

    monkeypatch.setattr(library_step.lx, "load", lambda name: Index())
    peaks = [(69, 1), (85, 1), (110, 1), (138, 1)]
    layer, summary = library_step.run_library(peaks, 195.0877, True, ["MassBank"], adducts=["[M+H]+"], formulas=[])
    assert summary.passed and layer.candidates[0].evidence["library"]["accession"] == "hit-7"
    assert summary.reference["accession"] == "hit-7" and summary.reference["peaks"][0][0] == 7
    json.dumps(summary.to_dict(), allow_nan=False)


def test_spectrum_figures_are_json_safe_and_optional():
    peaks = [(69.0, 20.0), (85.0, 100.0)]
    reference = {"peaks": [(69.01, 1.0), (100.0, 0.5)], "name": "example", "library": "MassBank", "accession": "A1", "precursor_mz": 195.0}
    result = visualization.spectrum_figures(peaks, 195.0, 0.02, reference)
    assert len(result["query"]["peaks"]) == 2 and result["library"]["accession"] == "A1"
    assert result["library"]["figure"]["plotly"]["data"]
    json.dumps(result, allow_nan=False)
    assert visualization.spectrum_figures(peaks, 195.0, 0.02)["library"] is None


def test_missing_library_reference_does_not_raise(monkeypatch):
    monkeypatch.setattr(library_step.lx, "load", lambda name: (_ for _ in ()).throw(OSError("unavailable")))
    assert library_step._reference_spectrum({"library": "MassBank", "row": 1}) is None
