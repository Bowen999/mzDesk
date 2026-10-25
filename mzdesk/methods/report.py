"""RepoRT (retention-time repository, CC BY-SA 4.0) as a source of chromatography methods.

``python -m mzdesk.methods.report_setup`` downloads the needed files of a pinned RepoRT commit and builds the local index
(``report_index.npz`` + ``report_datasets.json``) in ``report_dir()``. This module loads it and turns a dataset's metadata into the
short Column / Mobile phase / Reference texts of the Known methods page.

RepoRT: Kretschmer, Harrieder, Hoffmann, Boecker, Witting. Nat. Methods 21, 153-155 (2024), https://doi.org/10.1038/s41592-023-02143-z
Data: https://github.com/michaelwitting/RepoRT, licence CC BY-SA 4.0.
"""
import json
import os
import re
from dataclasses import dataclass
from pathlib import Path

from mzdesk.methods import compound_index as ci

SETUP_COMMAND = "python -m mzdesk.methods.report_setup"
CITATION = ("F. Kretschmer, E.-M. Harrieder, M. A. Hoffmann, S. Böcker, M. Witting. RepoRT: a comprehensive repository for small "
            "molecule retention times. Nat. Methods 21, 153-155 (2024). https://doi.org/10.1038/s41592-023-02143-z")
LICENCE = "CC BY-SA 4.0"
LICENCE_URL = "https://creativecommons.org/licenses/by-sa/4.0/"
REPO_URL = "https://github.com/michaelwitting/RepoRT"
INDEX_FILE = "report_index.npz"
DATASETS_FILE = "report_datasets.json"

SOLVENTS = (("h2o", "water"), ("meoh", "MeOH"), ("acn", "ACN"), ("iproh", "IPA"), ("acetone", "acetone"), ("hex", "hexane"),
            ("chcl3", "CHCl3"), ("ch2cl2", "CH2Cl2"), ("hept", "heptane"))
ADDITIVES = (("formic", "formic"), ("acetic", "acetic"), ("trifluoroacetic", "TFA"), ("phosphor", "phosphoric"),
             ("nh4ac", "ammonium acetate"), ("nh4form", "ammonium formate"), ("nh4carb", "ammonium carbonate"),
             ("nh4bicarb", "ammonium bicarbonate"), ("nh4f", "ammonium fluoride"), ("nh4oh", "NH4OH"), ("trieth", "triethylamine"),
             ("triprop", "tripropylamine"), ("tribut", "tributylamine"), ("nndimethylhex", "N,N-dimethylhexylamine"),
             ("medronic", "medronic acid"), ("heptafluorobutyric", "HFBA"))


def report_dir():
    """Where the RepoRT index lives. MZDESK_REPORT_DIR overrides it (default: <methods dir>/report)."""
    return Path(os.environ.get("MZDESK_REPORT_DIR") or ci.methods_dir() / "report")


def installed(path=None):
    d = Path(path) if path else report_dir()
    return (d / INDEX_FILE).is_file() and (d / DATASETS_FILE).is_file()


@dataclass
class Report:
    index: ci.CompoundIndex
    datasets: list                 # dataset dicts (see report_setup.parse_dataset), position = the index's entry reference
    meta: dict

    @property
    def commit(self):
        return self.meta.get("commit", "")


_cache = {}


def load(path=None):
    """The installed RepoRT index (cached; reloaded when the files change), or None when it is not installed."""
    d = Path(path) if path else report_dir()
    if not installed(d):
        return None
    stamp = ((d / INDEX_FILE).stat().st_mtime_ns, (d / DATASETS_FILE).stat().st_mtime_ns)
    hit = _cache.get(str(d))
    if hit and hit[0] == stamp:
        return hit[1]
    index = ci.CompoundIndex.load(d / INDEX_FILE)
    with open(d / DATASETS_FILE, encoding="utf-8") as f:
        payload = json.load(f)
    report = Report(index=index, datasets=payload["datasets"], meta=index.meta)
    _cache[str(d)] = (stamp, report)
    return report


# ---------------------------------------------------------------- texts
def _g(v):
    return f"{v:g}"


def column_text(ds):
    """e.g. "Waters ACQUITY UPLC HSS T3 (L1), 2.1 × 150 mm, 1.8 µm, 40 °C"."""
    c = ds.get("column") or {}
    name = c.get("name") or ""
    if c.get("usp"):
        name = f"{name} ({c['usp']})" if name else c["usp"]
    parts = [name] if name else []
    if c.get("id") and c.get("length"):
        parts.append(f"{_g(c['id'])} × {_g(c['length'])} mm")
    elif c.get("length"):
        parts.append(f"{_g(c['length'])} mm long")
    if c.get("particle"):
        parts.append(f"{_g(c['particle'])} µm")
    if c.get("temp") is not None:
        parts.append(f"{_g(c['temp'])} °C")
    return ", ".join(parts) or "—"


def _eluent_text(e):
    solvents = [(abbr, e["solvents"][key]) for key, abbr in SOLVENTS if e["solvents"].get(key)]
    if len(solvents) == 1:
        mix = solvents[0][0]
    else:
        mix = "/".join(a for a, _ in solvents) + " " + ":".join(_g(round(v, 1)) for _, v in solvents)
    adds = {key: (v, unit) for key, v, unit in e.get("additives", [])}
    extra = []
    for key, label in ADDITIVES:
        if key in adds:
            v, unit = adds[key]
            extra.append(f"{_g(v)}{'' if unit == '%' else ' '}{unit} {label}")
    text = mix + "".join(" + " + x for x in extra)
    if e.get("pH"):
        text += f" (pH {_g(e['pH'])})"
    return text


def eluents_text(ds):
    """"A: water + 0.1% formic; B: ACN + 0.1% formic" (only the eluents that have a composition)."""
    out = []
    for letter in "ABCD":
        e = (ds.get("eluents") or {}).get(letter)
        if e and e.get("solvents"):
            out.append(f"{letter}: {_eluent_text(e)}")
    return "; ".join(out)


def _furthest(v):
    """Index of the first point whose value is furthest from the start value."""
    return max(range(len(v)), key=lambda i: (abs(v[i] - v[0]), -i))


def gradient_text(ds):
    """"5→95% B in 20 min", "isocratic 40% B"; from the gradient table, else from the start / end proportions of the metadata."""
    g = ds.get("gradient") or {}
    t = g.get("t") or []
    letters = [L for L in "ABCD" if g.get(L)]
    if len(t) >= 2 and letters:
        varying = [L for L in letters if max(g[L]) - min(g[L]) > 1e-9]
        if not varying:
            return "isocratic " + ", ".join(f"{_g(round(g[L][0], 1))}% {L}" for L in letters if g[L][0])
        show = ["B"] if varying in (["A", "B"], ["B"]) else varying              # with two eluents A is the complement of B
        t_end = t[0]
        parts = []
        for L in show:
            i = _furthest(g[L])
            parts.append(f"{_g(round(g[L][0], 1))}→{_g(round(g[L][i], 1))}% {L}")
            t_end = max(t_end, t[i])
        return f"{', '.join(parts)} in {_g(round(t_end - t[0], 1))} min"
    start, end = ds.get("start") or {}, ds.get("end") or {}
    moved = [L for L in "ABCD" if start.get(L, 0) != end.get(L, 0)]
    if moved:
        L = "B" if "B" in moved else moved[0]
        return f"{_g(start.get(L, 0))}→{_g(end.get(L, 0))}% {L}"
    return ""


def mobile_phase_text(ds):
    """"A: water + 0.1% formic; B: ACN + 0.1% formic; 5→95% B in 20 min; 0.4 mL/min"."""
    parts = [p for p in (eluents_text(ds), gradient_text(ds)) if p]
    flow = (ds.get("column") or {}).get("flow")
    if flow:
        parts.append(f"{_g(flow)} mL/min")
    return "; ".join(parts) or "—"


_DOI = re.compile(r"10\.\d{4,9}/[^\s,;]+")


def reference(ds):
    """(text, href or None): "RepoRT 0001 · 10.1007/s11306-014-0638-x" linking to the DOI, else the source name (linked when the URL is a web address)."""
    url = (ds.get("url") or "").strip()
    doi = _DOI.search(url)
    if doi:
        label = doi.group(0).rstrip(".")
        href = f"https://doi.org/{label}"
    else:
        label = ds.get("source") or url
        href = url if url.startswith(("http://", "https://")) else None
        if len(label) > 70:
            label = label[:67] + "…"
    return f"RepoRT {ds['id']}" + (f" · {label}" if label else ""), href


def dataset_title(ds):
    """"0001 FEM_short (RP)"."""
    return f"{ds['id']} {ds.get('name') or ''}".strip() + (f" ({ds['type']})" if ds.get("type") else "")


def metadata_gaps(ds):
    """How much a dataset leaves unsaid (used to list the better documented datasets first)."""
    c = ds.get("column") or {}
    gaps = sum(1 for k in ("name", "length", "id", "particle", "temp", "flow") if not c.get(k) and c.get(k) != 0)
    return gaps + (0 if (ds.get("gradient") or {}).get("t") else 1) + (0 if eluents_text(ds) else 1)
