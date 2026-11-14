"""Smoke test of the engines on a real machine: Windows (.github/workflows/windows.yml), macOS and Linux (.github/workflows/posix.yml).
Also fine to run by hand after `mzdesk setup`.

    python validation/windows/smoke.py msfinder|frigid|neims|repo|cfmid|ccs|redu|buddy|deep

Caffeine, [M+H]+ 195.0877. Each check prints what it got and exits 1 when the answer is wrong."""
import os
import sys

# `python validation/windows/smoke.py` puts only this folder on sys.path; the checkout root holds the mzdesk package (CI does not pip-install it)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

CAFFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
PEAKS = [(138.0662, 100.0), (110.0713, 18.0), (123.0553, 9.0), (83.0604, 7.0), (195.0877, 40.0), (69.0447, 5.0)]


def msfinder():
    from mzdesk.ms2 import msfinder_engine as mf
    out = mf.predict(195.0877, True, "[M+H]+", PEAKS)
    formulas = [f.get("formula") for f in out["formulas"]]
    print("formulas:", formulas[:5], "structures:", len(out["structures"]), "install copy:", out["install_copy"], f"{out['seconds']:.1f} s")
    assert "C8H10N4O2" in formulas, "caffeine's formula is not among the formulas"
    again = mf.predict(195.0877, True, "[M+H]+", PEAKS)
    scores = lambda r: [round(s["score"], 4) for s in r["structures"][:10]]          # noqa: E731 - the thread fix: same scores on every run
    assert scores(out) == scores(again), "structure scores differ between two runs (thread fix missing?)"


def frigid():
    from mzdesk.ms2 import frigid_engine as fe
    # a CI runner has 4 slow cores and no GPU: give it 15 minutes (the page allows 300 s) and print how long it really took
    out = fe.predict("C8H10N4O2", 195.0877, "[M+H]+", PEAKS, "quick", timeout=900)
    print("device:", out["device"], "candidates:", len(out["candidates"]), f"{out['seconds']:.0f} s", "peak RSS MB:", out["peak_rss_mb"])
    assert out["candidates"], "FRIGID returned no candidate"


def neims():
    from mzdesk.ms2 import neims_engine as ne
    out = ne.predict(CAFFEINE)
    print("neims keys:", sorted(out)[:6])
    assert out


def repo():
    from mzdesk.methods import report
    assert report.installed(), "the RepoRT index is not installed"
    print("RepoRT index:", report.report_dir())


def buddy():
    from mzdesk.ms1 import msbuddy_engine as bd
    import msbuddy
    out = bd.annotate_formula(195.0877, True, "[M+H]+", None, "", tuple(PEAKS), ())
    formulas = [c["formula"] for c in out["candidates"]]
    print("msbuddy", getattr(msbuddy, "__version__", "?"), "formulas:", formulas[:5])
    assert "C8H10N4O2" in formulas, "BUDDY did not rank caffeine's formula"


def cfmid():
    from mzdesk.ms2 import cfmid_engine as cf
    out = cf.run_cfmid(CAFFEINE, "[M+H]+")
    sizes = [len(s["peaks"]) for s in out["spectra"]]
    print("engine:", out["engine"], "precursor:", round(out["precursor_mz"], 4), "peaks per energy:", sizes)
    assert abs(out["precursor_mz"] - 195.0877) < 0.001 and out["spectra"] and all(sizes), "CFM-ID returned no spectrum"


def ccs():
    from mzdesk.orthogonal import ccs as c
    assert c.installed(), "the CCS index is not installed"
    predictor = c.load_predictor()
    print("CCS index:", c.ccs_dir(), "predictor loaded:", predictor is not None)
    assert predictor is not None, "the CCS predictor did not load"


def redu():
    from mzdesk.gnps import compound_db_check as g
    assert g.redu_installed(), "the Pan-ReDU metadata is not downloaded"
    import gzip
    with gzip.open(g.REDU_SLIM, "rt", encoding="utf-8") as f:
        rows = sum(1 for _ in f) - 1
    print("Pan-ReDU:", g.REDU_SLIM, "rows:", rows, "snapshot:", g.redu_cache_date(g.REDU_SLIM))
    assert rows > 100000, "the Pan-ReDU table is far too small"


def deep():
    """MS2DeepScore and Spec2Vec (torch / gensim; the model downloads on first use). Errors become a note, so check for a score."""
    from mzdesk.ms2 import ms2_similarity as s
    other = [(138.0662, 90.0), (110.0713, 25.0), (195.0877, 30.0), (83.0604, 5.0)]
    rows = [r for r in s.compute_all_similarities(PEAKS, other, 195.0877, 195.0877) if r["category"] == s.mn.CAT_DEEP]
    for r in rows:
        print(r["method"], r["score"], r["note"][:80])
    if s.ms2deepscore_unsupported():                # Intel Mac: MS2DeepScore says why it is left out, Spec2Vec still scores
        rows = [r for r in rows if not (r["method"] == s.mn.MS2DEEPSCORE and r["note"].startswith("skipped: not available on Intel Macs"))]
        assert len(rows) == 1, "Spec2Vec is missing, or MS2DeepScore did not say it is left out on an Intel Mac"
    assert rows and all(isinstance(r["score"], float) for r in rows), "a deep-learning score is missing"


if __name__ == "__main__":
    name = sys.argv[1] if len(sys.argv) > 1 else ""
    checks = {"msfinder": msfinder, "frigid": frigid, "neims": neims, "repo": repo, "cfmid": cfmid, "ccs": ccs, "redu": redu, "buddy": buddy, "deep": deep}
    if name not in checks:
        sys.exit("usage: python validation/windows/smoke.py " + "|".join(checks))
    try:
        checks[name]()
    except Exception as exc:
        if getattr(exc, "log", ""):                 # MsFinderError / FrigidError / CFMIDError carry the engine's own output: show it, CI has nothing else
            print("---- engine output ----\n" + exc.log + "\n---- end ----", flush=True)
        raise
    print("ok")
