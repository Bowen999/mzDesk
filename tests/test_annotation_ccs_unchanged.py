"""Without an observed CCS, ESI MS/MS Annotation gives exactly the result it gave before the CCS check existed.

Three hand-made runs (the Rig of test_pipeline.py replaces the engines; one run uses the real re-rank) are snapshotted: answer, every candidate field,
the steps (status and note), the warnings, the context, and what step 6 showed while the job ran (its note and blocks). The snapshot was taken with
the code of v0.3.0, before the CCS check, and is stored in tests/fixtures/esi_without_ccs.json. This file imports nothing from the CCS check, so the
same file can write the snapshot from an older checkout:  python tests/test_annotation_ccs_unchanged.py --write

If a later change moves a rank, a score, a level, a note or a block on purpose, regenerate the file in the same commit and say so.
"""
import json
import sys
from pathlib import Path

import pytest

from mzdesk.annotation import pipeline
from test_pipeline import CAFFEINE, ETHANOL, ETHER, REAL_RERANK, Rig, form, isomer_result, library_candidate, structure_candidate

GOLDEN = Path(__file__).parent / "fixtures" / "esi_without_ccs.json"
SCENARIOS = ("library gate passed", "MS-FINDER and FRIGID", "real re-rank")


def _run(name, mp):
    """One hand-made run -> (result, what step 6 emitted when it finished)."""
    rig = Rig(mp)
    if name == "library gate passed":
        rig.library([library_candidate("lib1", CAFFEINE, "Caffeine", 0.93, True), library_candidate("lib2", ETHER, "Dimethyl ether", 0.72, False, formula="C2H6O")],
                    passed=True)
        rig.isomers(isomer_result("single", "c1"))
        inputs = form(libraries=["MassBank"])
    elif name == "MS-FINDER and FRIGID":
        rig.msfinder([structure_candidate("msfinder", "m1", CAFFEINE, 0.8), structure_candidate("msfinder", "m2", ETHANOL, 0.76, formula="C2H6O"),
                      structure_candidate("msfinder", "m3", ETHER, 0.3, formula="C2H6O")], good_enough=False)
        rig.frigid([structure_candidate("frigid", "f1", ETHER, 0.2, formula="C8H10N4O2")])
        inputs = form()
    else:
        mp.setattr(pipeline, "rerank", REAL_RERANK)
        rig.msfinder([structure_candidate("msfinder", "m1", CAFFEINE, 0.8), structure_candidate("msfinder", "m2", ETHANOL, 0.78, formula="C2H6O"),
                      structure_candidate("msfinder", "m3", ETHER, 0.3, formula="C2H6O")], good_enough=True)
        inputs = form(omics="Metabolomics")
    emitted = {}
    result = pipeline.run_query(pipeline.parse_inputs(inputs), progress=lambda n, status, note, blocks: emitted.update({n: (status, note, blocks)}) if n == 6 else None)
    return result, emitted[6]


def snapshot(result, step6):
    """The result without timings and without the display-only figures, as plain JSON."""
    def strip(x):
        if isinstance(x, dict):
            return {k: (0 if k == "seconds" else strip(v)) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [strip(v) for v in x]
        return x
    body = {k: v for k, v in result.items() if k != "visualization"}
    return json.loads(json.dumps({"result": strip(body), "step6": strip(list(step6))}, default=str))


@pytest.mark.parametrize("name", SCENARIOS)
def test_a_run_without_an_observed_ccs_is_what_it_was_before(name, monkeypatch):
    golden = json.loads(GOLDEN.read_text(encoding="utf-8"))
    result, step6 = _run(name, monkeypatch)
    assert snapshot(result, step6) == golden[name]


@pytest.mark.parametrize("name", SCENARIOS)
def test_no_candidate_row_or_step_mentions_ccs_when_none_was_given(name, monkeypatch):
    result, step6 = _run(name, monkeypatch)
    assert all("ccs" not in c for c in result["candidates"])
    assert "CCS" not in json.dumps([result["steps"], result["warnings"], step6], ensure_ascii=False)


if __name__ == "__main__":
    if "--write" not in sys.argv:
        sys.exit("usage: python tests/test_annotation_ccs_unchanged.py --write   (writes tests/fixtures/esi_without_ccs.json from the code that is importable now)")
    out = {}
    for scenario in SCENARIOS:
        mp = pytest.MonkeyPatch()
        try:
            out[scenario] = snapshot(*_run(scenario, mp))
        finally:
            mp.undo()
    GOLDEN.write_text(json.dumps(out, indent=1, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {GOLDEN} ({GOLDEN.stat().st_size} bytes)")
