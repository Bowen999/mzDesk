"""Japanese Pharmacopoeia 18th edition (bundled, mzdesk/data/known_methods/jp18_methods.jsonl.gz) as a source of chromatography methods.

The JP is a Ministerial Notification (not protected by copyright, Japanese Copyright Act Art. 13). One record per "Operating conditions" block
(1924 LC / GC methods); each record names the compounds it measures (`analyte_structures`). This module indexes those structures (same normalisation
and descriptors as RepoRT; the index is cached in the methods folder because MolGpKa takes ~10 s for the ~1000 structures) and writes the short
Column / Mobile phase / Detector / Reference texts of the Known methods page.
"""
import gzip
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path

from mzdesk.methods import compound_index as ci

DATA_FILE = Path(__file__).resolve().parents[1] / "data" / "known_methods" / "jp18_methods.jsonl.gz"
CACHE_FILE = "jp18_index.npz"
PDF_URL = "https://www.mhlw.go.jp/content/11120000/{file}"
VOLUMES = {"000904440.pdf": "Official Monographs A–L", "000904449.pdf": "Official Monographs M–Z", "000904450.pdf": "Crude Drugs and Related Drugs",
           "000945683.pdf": "General Notices and General Tests", "000904454.pdf": "General Information"}
LICENCE_NOTE = ("Japanese Pharmacopoeia, 18th edition (English version) with Supplements I and II, Ministry of Health, Labour and Welfare "
                "(Ministerial Notification; not protected by copyright under Article 13 of the Japanese Copyright Act). Parameters were extracted from the "
                "PDF text layer of the 'Operating conditions' blocks, with the printed page; the English version is a convenience translation, so check "
                "the printed page before relying on a value (errata published after the supplements are not applied).")


@dataclass
class JP:
    index: ci.CompoundIndex
    records: list
    by_block: dict


_cache = {}


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_records(path=DATA_FILE):
    with gzip.open(path, "rt", encoding="utf-8") as f:
        return [json.loads(line) for line in f if line.strip()]


def _rows(records):
    """(smiles, record, nan, name, role) for every structure a record gives: role 0 = what the test measures (`analyte_structures`),
    role 1 = the monograph's own substance (`substances`; an impurity or purity test is still a method for it). One entry per record and compound
    (the lowest role wins; the index merges the structures by connectivity key)."""
    rows = []
    for i, r in enumerate(records):
        for role, key in ((0, "analyte_structures"), (1, "substances")):
            seen = set()
            for a in r.get(key) or []:
                smiles = a.get("smiles") or a.get("parent_smiles")
                if smiles and smiles not in seen:
                    seen.add(smiles)
                    rows.append((smiles, i, float("nan"), a.get("pubchem_title") or a.get("name") or "", role))
    return rows


def load(path=DATA_FILE, cache_dir=None):
    """The JP18 records and their compound index. The index is cached in the methods folder (rebuilt when the data file or the charge method changes)."""
    cache_dir = Path(cache_dir) if cache_dir else ci.methods_dir()
    method = ci.best_charge_method()
    key = (str(path), Path(path).stat().st_mtime_ns, method, str(cache_dir))
    if key in _cache:
        return _cache[key]
    records = load_records(path)
    sha = _sha(path)
    cache = cache_dir / CACHE_FILE
    index = None
    if cache.is_file():
        try:
            cached = ci.CompoundIndex.load(cache)
            m = cached.meta
            if m.get("data_sha") == sha and m.get("charge_method") == method and m.get("schema") == ci.SCHEMA:
                cached.check()
                index = cached
        except Exception:  # noqa: BLE001 - unreadable or stale cache: rebuild
            index = None
    if index is None:
        index, _skipped, _fallbacks = ci.CompoundIndex.build(_rows(records), {"source": "JP18", "data_sha": sha}, method, dedupe_refs=True)
        try:
            index.save(cache)
        except OSError:
            pass                                    # read-only home: keep it in memory only
    jp = JP(index=index, records=records, by_block={r["block_id"]: r for r in records})
    _cache.clear()
    _cache[key] = jp
    return jp


# ---------------------------------------------------------------- texts
def _g(v):
    return f"{v:g}"


def _g1(v):
    s = f"{v:g}"
    return s if "." in s else s + ".0"


def column_text(r):
    """e.g. "C18, 4.6 × 150 mm, 5 µm, 40 °C"; GC: "GC-PEG, 0.32 mm × 30 m, 0.25 µm, temperature programme"."""
    parts = []
    if r.get("packing_class"):
        parts.append(r["packing_class"])
    cid, length = r.get("column_id_mm"), r.get("column_length_mm")
    if cid and length:
        if r.get("technique") == "GC" or length >= 1000:
            parts.append(f"{_g(cid)} mm × {_g(length / 1000)} m")
        else:
            parts.append(f"{_g(cid)} × {_g(length)} mm")
    if r.get("particle_or_film_um"):
        parts.append(f"{_g(r['particle_or_film_um'])} µm")
    temp = r.get("column_temp_c")
    if isinstance(temp, (int, float)):
        parts.append(f"{_g(temp)} °C")
    elif temp == "program":
        parts.append("temperature programme")
    elif temp == "room":
        parts.append("room temperature")
    return ", ".join(parts) or "—"


def mobile_phase_text(r):
    """LC: "isocratic; acetonitrile; pH 3.0; 1.0 mL/min". GC: "carrier gas helium; split 1:20"."""
    parts = []
    if r.get("technique") == "GC":
        parts.append("GC")
        if r.get("carrier_gas"):
            parts.append(f"carrier gas {r['carrier_gas']}")
        if r.get("split_ratio"):
            parts.append(f"split {r['split_ratio']}")
    else:
        if r.get("elution_mode"):
            parts.append(r["elution_mode"])
        if r.get("organic_solvents"):
            parts.append(" + ".join(r["organic_solvents"]))
        elif r.get("organic_solvents") == [] and r.get("elution_mode"):
            parts.append("no organic solvent")
        if r.get("mobile_phase_ph") is not None:
            parts.append(f"pH {_g1(r['mobile_phase_ph'])}")
    if r.get("flow_rate_value"):
        parts.append(f"{_g(r['flow_rate_value'])} {r.get('flow_rate_unit') or ''}".strip())
    return "; ".join(parts) or "—"


_FLUOR = re.compile(r"excitation wavelength:?\s*([\d.]+)\s*nm.*?emission wavelength:?\s*([\d.]+)\s*nm", re.S | re.I)


def detector_text(r):
    """"UV 220 nm", "fluorescence ex 380 / em 520 nm", "FID"."""
    d = r.get("detector_type")
    if not d:
        return "—"
    if d == "fluorescence":
        m = _FLUOR.search(r.get("condition_text") or "")
        return f"fluorescence ex {m.group(1)} / em {m.group(2)} nm" if m else "fluorescence"
    return f"{d} {_g(r['wavelength_nm'])} nm" if r.get("wavelength_nm") else d


def test_text(r):
    return f"{r.get('monograph') or ''}: {r.get('test_label') or ''}".strip(": ")


def reference(r):
    """(text, PDF link): "JP18 p. 623 (Official Monographs A–L)"; the link opens the MHLW PDF at that page."""
    volume = VOLUMES.get(r.get("file"))
    text = f"{r.get('source') or 'JP18'} p. {r.get('page_printed') or r.get('page')}" + (f" ({volume})" if volume else "")
    href = PDF_URL.format(file=r["file"]) + (f"#page={r['page']}" if r.get("page") else "") if r.get("file") else None
    return text, href


def _flat(text):
    return re.sub(r"\s+", " ", text or "").strip()


def condition_text(r, by_block):
    """The printed conditions on one line; for a cross-referenced block the referenced block's text follows."""
    text = _flat(r.get("condition_text"))
    ref = by_block.get(r.get("resolved_from") or "")
    if ref and ref is not r:
        text += f" ▸ Referenced conditions ({ref['block_id']}, {test_text(ref)}): {_flat(ref.get('condition_text'))}"
    return text


def _priority(r):
    """Assay first, then purity tests, then the rest; a plain "Assay" (the substance's own monograph) before numbered items of multi-component extracts."""
    label = (r.get("test_label") or "").lower()
    return (0 if label.startswith("assay") else 1 if label.startswith("purity") else 2, len(label))


def measures_text(r):
    """What the test measures, from the record's analytes (shown when a method is only for the monograph substance)."""
    names = [n for n in (r.get("analytes") or []) if n]
    text = ", ".join(names[:3]) + (" …" if len(names) > 3 else "")
    return text if len(text) <= 80 else text[:79] + "…"


def groups(records, rec_ids, roles=None):
    """Records (indices) merged by identical displayed conditions, best test first -> [[index, ...], ...].

    Several monographs often print the same conditions (17 extracts assay paeoniflorin with 3 different conditions); the page lists each distinct
    condition once and says how many monographs share it. ``roles`` {record: 0 or 1}: tests that measure the compound come before tests of its
    monograph that measure something else, and the two are never merged (nor are substance-only tests that measure different things).
    """
    roles = roles or {}
    out = {}
    for i in sorted(set(rec_ids), key=lambda i: (roles.get(i, 0), _priority(records[i]), i)):
        r = records[i]
        role = roles.get(i, 0)
        key = (role, measures_text(r) if role else None, column_text(r), mobile_phase_text(r), detector_text(r), _flat(r.get("condition_text")), r.get("resolved_from"))
        out.setdefault(key, []).append(i)
    return sorted(out.values(), key=lambda g: (roles.get(g[0], 0), _priority(records[g[0]]), g[0]))
