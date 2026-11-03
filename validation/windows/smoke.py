"""Smoke test of the engines on a real Windows machine (used by .github/workflows/windows.yml, also fine to run by hand after `mzdesk setup`).

    python validation/windows/smoke.py msfinder|frigid|neims|repo

Caffeine, [M+H]+ 195.0877. Each check prints what it got and exits 1 when the answer is wrong."""
import sys

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
    out = fe.predict("C8H10N4O2", 195.0877, "[M+H]+", PEAKS, "quick")
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


if __name__ == "__main__":
    name = sys.argv[1] if len(sys.argv) > 1 else ""
    checks = {"msfinder": msfinder, "frigid": frigid, "neims": neims, "repo": repo}
    if name not in checks:
        sys.exit("usage: python validation/windows/smoke.py " + "|".join(checks))
    checks[name]()
    print("ok")
