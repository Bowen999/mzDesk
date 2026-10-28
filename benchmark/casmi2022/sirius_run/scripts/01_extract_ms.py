#!/usr/bin/env python3
"""Blind extraction of MS1 + MS2 for CASMI 2022 challenges -> SIRIUS .ms files.

BLIND: only the columns challenge_id, file, rt_min, precursor_mz, polarity,
priority, mzml_path are read from the challenge table (see USECOLS). No
annotation columns are touched.

Run:  python3 01_extract_ms.py [--parentmass ms1|given] [--outdir DIR]
      (rerunnable; overwrites ms/*.ms)
--parentmass ms1   (default) >parentmass = measured MS1 monoisotopic m/z (falls back to the
                   given CASMI value if the monoisotopic peak is not found in MS1)
--parentmass given >parentmass = rounded CASMI precursor_mz (can be up to 8 ppm off)
The given value is always kept in extraction_report.tsv (precursor_given); the value written
to the .ms file is in parentmass_used / parentmass_source.
Needs: pyteomics, numpy, pandas  (system python3 has them)
"""
import os, sys, time, math, argparse
from collections import defaultdict
import numpy as np
import pandas as pd
from pyteomics import mzml

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.abspath(os.path.join(HERE, "..", ".."))          # benchmark/casmi2022
OUT = os.path.join(HERE, "..", "ms")
USECOLS = ["challenge_id", "file", "rt_min", "precursor_mz", "polarity", "priority", "mzml_path"]

PPM_TOL = 10.0          # precursor match (MS2) and isotope match (MS1)
ABS_TOL = 0.01          # Da fallback for MS2 precursor match
RT_TOLS = (0.3, 0.5)    # min; widen if nothing found
APEX_WIN = 0.1          # min; scans of other NCEs must lie within this of the apex scan
MS2_UPPER = 0.5         # keep fragments up to precursor + this (Da)
C13 = 1.00336
N_ISO = 5               # k = 0..4


def ppm(obs, ref):
    return (obs - ref) / ref * 1e6


def read_file(path):
    """One pass over the mzML. Returns (ms1 dict id->(rt,mz,int), ms2 list of dicts)."""
    ms1, ms2 = {}, []
    with mzml.read(path) as r:
        for s in r:
            lvl = s["ms level"]
            rt = float(s["scanList"]["scan"][0]["scan start time"])   # minutes
            mz = np.asarray(s["m/z array"], dtype=float)
            it = np.asarray(s["intensity array"], dtype=float)
            if lvl == 1:
                ms1[s["id"]] = (rt, mz, it)
            elif lvl == 2:
                p = s["precursorList"]["precursor"][0]
                si = p["selectedIonList"]["selectedIon"][0]
                iw = p.get("isolationWindow", {})
                ce = p.get("activation", {}).get("collision energy")
                ms2.append(dict(id=s["id"], rt=rt, mz=mz, it=it,
                                pmz=float(si["selected ion m/z"]),
                                tgt=float(iw.get("isolation window target m/z", si["selected ion m/z"])),
                                pint=float(si.get("peak intensity", 0.0)),
                                ce=float(ce) if ce is not None else None,
                                parent=p.get("spectrumRef")))
    return ms1, ms2


def pick_ms2(ms2, rt, pmz):
    """Return (selected scans, note, tolerance-label)."""
    notes = []
    for mode in ("ppm", "abs"):
        for rtol in RT_TOLS:
            c = []
            for s in ms2:
                d = s["pmz"] - pmz
                ok = abs(ppm(s["pmz"], pmz)) <= PPM_TOL if mode == "ppm" else abs(d) <= ABS_TOL
                if ok and abs(s["rt"] - rt) <= rtol:
                    c.append(s)
            if c:
                if mode == "abs":
                    notes.append("precursor matched by 0.01 Da (not 10 ppm)")
                if rtol != RT_TOLS[0]:
                    notes.append("RT window widened to +-%.1f min" % rtol)
                return c, notes
    return [], notes


def select_apex_set(c):
    """apex = most intense precursor; then per NCE the scan closest in time to the apex."""
    apex = max(c, key=lambda s: s["pint"])
    best = {}
    for s in c:
        if abs(s["rt"] - apex["rt"]) > APEX_WIN:
            continue
        k = s["ce"]
        if k not in best or abs(s["rt"] - apex["rt"]) < abs(best[k]["rt"] - apex["rt"]):
            best[k] = s
    sel = sorted(best.values(), key=lambda s: (s["ce"] is None, s["ce"]))
    return apex, sel


def iso_pattern(ms1scan, pmz):
    rt, mz, it = ms1scan
    peaks = []
    for k in range(N_ISO):
        t = pmz + C13 * k
        tol = t * PPM_TOL * 1e-6
        i = np.searchsorted(mz, t - tol), np.searchsorted(mz, t + tol)
        if i[1] > i[0]:
            j = i[0] + int(np.argmax(it[i[0]:i[1]]))
            if it[j] > 0:
                peaks.append((float(mz[j]), float(it[j]), k))
    return peaks


def nearest_ms1(ms1, rt):
    return min(ms1.items(), key=lambda kv: abs(kv[1][0] - rt))


def main():
    global OUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--parentmass", choices=["ms1", "given"], default="ms1")
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    if args.outdir:
        OUT = os.path.abspath(args.outdir)
    t0 = time.time()
    df = pd.read_csv(os.path.join(BASE, "metadata", "casmi2022_challenges.csv"), usecols=USECOLS)
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for fpath, grp in df.groupby("mzml_path", sort=False):
        ms1, ms2 = read_file(os.path.join(BASE, fpath))
        for r in grp.itertuples(index=False):
            cid = int(r.challenge_id); rt = float(r.rt_min); pmz = float(r.precursor_mz)
            neg = str(r.polarity).lower().startswith("neg")
            notes = []
            c, n = pick_ms2(ms2, rt, pmz)
            notes += n
            row = dict(id=cid, priority=r.priority, file=r.file, rt_given=rt, precursor_given=pmz,
                       polarity=r.polarity, n_ms2_scans_used=0, ms2_rt_used="", ms2_precursor_ppm_error="",
                       ms2_rt_offset_min="", ms1_rt_used="", n_ms1_isotope_peaks=0, n_ms2_peaks=0)
            sel, apex, iso, ms1rt = [], None, [], None
            if not c:
                notes.append("NO_MS2")
            else:
                apex, sel = select_apex_set(c)
                # MS1: survey scan that triggered the apex MS2 (spectrumRef), else nearest in RT
                key = apex["parent"] if apex["parent"] in ms1 else None
                if key is None:
                    key = nearest_ms1(ms1, apex["rt"])[0]; notes.append("MS1 by nearest RT")
                scan = ms1[key]
                iso = iso_pattern(scan, pmz); ms1rt = scan[0]
                if not iso or iso[0][2] != 0:
                    # try neighbouring survey scans (+-3) for one that contains the precursor
                    near = sorted(ms1.items(), key=lambda kv: abs(kv[1][0] - apex["rt"]))[:7]
                    cand = [(iso_pattern(v, pmz), v[0]) for _, v in near]
                    cand = [x for x in cand if x[0] and x[0][0][2] == 0]
                    if cand:
                        iso, ms1rt = max(cand, key=lambda x: x[0][0][1]); notes.append("MS1 from neighbouring scan")
                    else:
                        notes.append("PRECURSOR_NOT_IN_MS1" if not iso else "MONOISOTOPIC_NOT_IN_MS1")
            # ---- precursor m/z written to the .ms file
            pm_out, pm_src = pmz, "given"
            if args.parentmass == "ms1" and iso and iso[0][2] == 0:
                pm_out, pm_src = iso[0][0], "ms1_mono"
            row["parentmass_used"] = "%.5f" % pm_out
            row["parentmass_source"] = pm_src
            # ---- write .ms
            lines = [">compound challenge_%03d" % cid,
                     ">parentmass %.5f" % pm_out,
                     ">ionization [M+?]%s" % ("-" if neg else "+"),
                     ">rt %.1fs" % (rt * 60.0), ""]
            if iso:
                lines.append(">ms1")
                lines += ["%.6f %.2f" % (m, i) for m, i, _ in iso]
                lines.append("")
            tot = 0
            for s in sel:
                keep = (s["it"] > 0) & (s["mz"] <= pmz + MS2_UPPER)
                m, i = s["mz"][keep], s["it"][keep]
                if s["ce"] is not None:
                    lines.append(">collision %g" % s["ce"])
                else:
                    lines.append(">ms2")
                lines += ["%.6f %.2f" % (a, b) for a, b in zip(m, i)]
                lines.append("")
                tot += len(m)
            if sel and tot == 0:
                notes.append("MS2_EMPTY")
            with open(os.path.join(OUT, "challenge_%03d.ms" % cid), "w") as f:
                f.write("\n".join(lines))
            if sel:
                row.update(n_ms2_scans_used=len(sel),
                           ms2_rt_used=";".join("%.4f" % s["rt"] for s in sel),
                           ms2_precursor_ppm_error="%.2f" % ppm(apex["pmz"], pmz),
                           ms2_rt_offset_min="%.4f" % (apex["rt"] - rt),
                           n_ms2_peaks=tot)
                if ms1rt is not None:
                    row["ms1_rt_used"] = "%.4f" % ms1rt
                row["n_ms1_isotope_peaks"] = len(iso)
            row["status"] = "OK" if (sel and iso and iso[0][2] == 0 and tot > 0) else "FAIL"
            if notes and row["status"] == "OK":
                row["status"] = "OK_WITH_NOTES"
            row["notes"] = "; ".join(notes)
            rows.append(row)
        print("%-45s %3d challenges  %.0fs" % (os.path.basename(fpath), len(grp), time.time() - t0), flush=True)
    rep = pd.DataFrame(rows).sort_values("id")
    cols = ["id", "priority", "file", "rt_given", "precursor_given", "polarity", "n_ms2_scans_used",
            "ms2_rt_used", "ms2_precursor_ppm_error", "ms2_rt_offset_min", "ms1_rt_used",
            "n_ms1_isotope_peaks", "n_ms2_peaks", "parentmass_used", "parentmass_source", "status", "notes"]
    rep[cols].to_csv(os.path.join(OUT, "extraction_report.tsv"), sep="\t", index=False)
    print("done: %d challenges, %.0f s" % (len(rep), time.time() - t0))


if __name__ == "__main__":
    main()
