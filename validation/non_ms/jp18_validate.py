"""Non-MS detection rules against the JP18 LC methods (mzdesk/data/known_methods/jp18_methods.jsonl.gz).

  python validation/non_ms/jp18_validate.py            # prints the three checks
  python validation/non_ms/jp18_validate.py --misses   # also lists the largest wavelength differences and the analytes predicted "none"

(1) analytes measured by UV / PDA: how often a chromophore (strong or weak) is predicted
(2) analytes measured by RI, post-column derivatization (reaction coil, o-phthalaldehyde, ninhydrin) or other non-optical detection:
    how often "none" / "weak" is predicted.  JP18 has no ELSD / CAD method; conductivity is listed apart (it detects ions, not chromophores)
(3) |JP wavelength - estimated lambda_max| for the UV analytes. JP often picks a lower wavelength than lambda_max for sensitivity
    (254 nm is a common default, 200-230 nm for end absorption), so this is context, not an error rate.
Counts are per unique analyte (SMILES); (3) is per unique analyte x wavelength. Records whose text mentions derivatization are not counted as UV.
"""
import collections
import gzip
import json
import re
import statistics
import sys
from pathlib import Path

from rdkit import RDLogger

from mzdesk.platform_fit import non_ms_rules as nm

RDLogger.DisableLog("rdApp.*")
DATA = Path(__file__).resolve().parents[2] / "mzdesk" / "data" / "known_methods" / "jp18_methods.jsonl.gz"
DERIVATIZED = re.compile(r"reaction coil|post-column|derivati[sz]|ninhydrin|o-phthal", re.I)
INORGANIC = re.compile(r"Cl|Ca|Cu|Co|Cr|Cd|Cs|Ce")
_cache = {}


def analyze(smiles):
    if smiles not in _cache:
        _cache[smiles] = nm.analyze(smiles)
    return _cache[smiles]


def has_carbon(smiles):
    return bool(re.search(r"[Cc]", INORGANIC.sub("", smiles)))


def load():
    """-> (UV items, universal items): dicts with name, smiles, wavelength, detector, group."""
    uv, universal = [], []
    for r in (json.loads(line) for line in gzip.open(DATA, "rt")):
        if r["technique"] != "LC":
            continue
        text = r["condition_text"] or ""
        derivatized = bool(DERIVATIZED.search(text))
        detector, wl = r["detector_type"], r["wavelength_nm"]
        for s in r["analyte_structures"] or []:
            smi = s.get("smiles")
            if not smi or not has_carbon(smi):
                continue
            item = {"name": s["name"], "smiles": smi, "wl": wl, "detector": detector, "monograph": r["monograph"]}
            if detector in ("UV", "PDA") and not derivatized:
                uv.append(item)
            elif detector == "RI":
                universal.append({**item, "group": "RI"})
            elif derivatized or (detector == "visible" and wl == 570):
                universal.append({**item, "group": "post-column / derivatization"})
            elif detector == "conductivity":
                universal.append({**item, "group": "conductivity (ions)"})
    return uv, universal


def unique(items, key):
    seen = {}
    for i in items:
        seen.setdefault(key(i), i)
    return list(seen.values())


def tally(items):
    c = collections.Counter(analyze(i["smiles"]).strength for i in items)
    n = sum(c.values())
    return n, c


def distribution(rows, label):
    d = [abs(x) for x, *_ in rows]
    if not d:
        return
    bins = [(0, 10), (10, 20), (20, 30), (30, 50), (50, 1e9)]
    parts = [f"{lo}-{hi if hi < 1e8 else 'inf'}: {100 * sum(lo <= v < hi for v in d) / len(d):.0f} %" for lo, hi in bins]
    signed = statistics.median(x for x, *_ in rows)
    print(f"  {label}: n={len(d)}, median |d| {statistics.median(d):.0f} nm, mean {statistics.mean(d):.0f} nm, median JP - ours {signed:+.0f} nm | " + ", ".join(parts))


def main(show_misses=False):
    uv, universal = load()
    print("(1) UV / PDA analytes: predicted chromophore")
    u = unique(uv, lambda i: i["smiles"])
    n, c = tally(u)
    print(f"  unique analytes n={n}: strong {c['strong']} ({100 * c['strong'] / n:.1f} %), weak {c['weak']} ({100 * c['weak'] / n:.1f} %), "
          f"none {c['none']} ({100 * c['none'] / n:.1f} %); chromophore predicted for {100 * (c['strong'] + c['weak']) / n:.1f} %")
    n, c = tally(uv)
    print(f"  records x analytes n={n}: chromophore predicted for {100 * (c['strong'] + c['weak']) / n:.1f} %")

    print("(2) analytes measured without a UV chromophore detector: predicted none / weak")
    groups = collections.defaultdict(list)
    for i in universal:
        groups[i["group"]].append(i)
    for g, items in groups.items():
        n, c = tally(unique(items, lambda i: i["smiles"]))
        print(f"  {g}: n={n}, none {c['none']}, weak {c['weak']}, strong {c['strong']}")
    optical = [i for g in ("RI", "post-column / derivatization") for i in groups[g]]
    n, c = tally(unique(optical, lambda i: i["smiles"]))
    print(f"  RI + derivatization together: n={n}, none / weak {100 * (c['none'] + c['weak']) / n:.0f} %")

    print("(3) |JP wavelength - lambda_max| (UV / PDA analytes with a predicted chromophore)")
    rows = [(i["wl"] - analyze(i["smiles"]).lam, i, analyze(i["smiles"])) for i in uv if i["wl"] and analyze(i["smiles"]).lam]
    uq = list({(i["smiles"], i["wl"]): (d, i, a) for d, i, a in rows}.values())
    distribution(uq, "unique analyte x wavelength, all")
    distribution([r for r in uq if r[1]["wl"] >= 230], "JP wavelength >= 230 nm")
    distribution([r for r in uq if r[1]["wl"] < 230], "JP wavelength < 230 nm (end-absorption choice)")
    for m in ("Woodward–Fieser", "Fieser–Kuhn", "Scott", "class table", "generic"):
        distribution([r for r in uq if r[1]["wl"] >= 230 and r[2].lam_method == m], f"  method {m}, JP >= 230 nm")

    if show_misses:
        print("UV analytes predicted none:")
        for i in unique([i for i in uv if analyze(i["smiles"]).strength == "none"], lambda i: i["smiles"]):
            print(f"  {i['name']} | JP {i['wl']} nm | {i['monograph'][:30]}")
        print("largest wavelength differences (JP >= 230 nm):")
        for d, i, a in sorted((r for r in uq if r[1]["wl"] >= 230), key=lambda r: -abs(r[0]))[:25]:
            print(f"  {d:+.0f} nm  JP {i['wl']}, ours {a.lam:.0f} ({a.lam_method}) {i['name']}: {[c.name for c in a.chromophores if c.strength == a.strength][:2]}")


if __name__ == "__main__":
    main("--misses" in sys.argv)
