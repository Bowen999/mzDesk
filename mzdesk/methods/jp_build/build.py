"""Orchestration: PDFs -> text -> blocks -> records -> merged 'current JP18' record set (no network)."""
from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from . import analytes as an
from . import header as hd
from . import merge as mg
from . import parse as ps
from . import segment as sg

SOURCES = [  # (stem, source label, default part, file description, id prefix, order)
    ("000945683", "JP18", "gen", "General Notices + General Tests", "JP18"),
    ("000904440", "JP18", "mono", "Official Monographs A-L", "JP18"),
    ("000904449", "JP18", "mono", "Official Monographs M-Z", "JP18"),
    ("000904450", "JP18", "crude", "Crude Drugs and Related Drugs", "JP18"),
    ("000904454", "JP18", "geninfo", "General Information", "JP18"),
    ("001109611", "JP18 Supplement I", "gen", "Supplement I", "JP18S1"),
    ("001354067", "JP18 Supplement II", "gen", "Supplement II", "JP18S2"),
]
ALT_NAMES = {"001354067": ["jp18_supp2_en"]}          # local copies of the Supplement II PDF may use this name
SOURCE_URL = {"000945683": "https://www.mhlw.go.jp/content/11120000/000945683.pdf",
              "000904440": "https://www.mhlw.go.jp/content/11120000/000904440.pdf",
              "000904449": "https://www.mhlw.go.jp/content/11120000/000904449.pdf",
              "000904450": "https://www.mhlw.go.jp/content/11120000/000904450.pdf",
              "000904454": "https://www.mhlw.go.jp/content/11120000/000904454.pdf",
              "001109611": "https://www.mhlw.go.jp/content/11120000/001109611.pdf",
              "001354067": "https://www.mhlw.go.jp/content/11120000/001354067.pdf"}
PART_CODE = {"mono": "ML", "crude": "CD", "gen": "GT", "geninfo": "GI"}
PART_NAME = {"mono": "monograph", "crude": "crude drug / extract", "gen": "general test", "geninfo": "general information"}
SOURCE_ORDER = {"JP18": 0, "JP18 Supplement I": 1, "JP18 Supplement II": 2}
SECTION_CANON = {"Uniformity of dosage unit": "Uniformity of dosage units", "Residue on Ignition": "Residue on ignition",
                 "Method of Preparation": "Method of preparation"}


def find_pdf(folder: Path, stem: str) -> Path | None:
    for name in [stem] + ALT_NAMES.get(stem, []):
        for cand in (folder / f"{name}.pdf", folder / f"{name}_en.pdf"):
            if cand.exists():
                return cand
    return None


def extract_all(folder: Path, work: Path, force: bool = False, log=print) -> dict[str, Path]:
    from .extract import extract_pdf
    text_dir = work / "text"
    text_dir.mkdir(parents=True, exist_ok=True)
    out = {}
    for stem, source, _part, desc, _pfx in SOURCES:
        pdf = find_pdf(folder, stem)
        if pdf is None:
            raise FileNotFoundError(f"{stem}.pdf not found in {folder}")
        txt = text_dir / f"{stem}.txt"
        code_mtime = Path(__file__).with_name("extract.py").stat().st_mtime       # re-extract when the extractor changed
        if force or not txt.exists() or txt.stat().st_mtime < max(pdf.stat().st_mtime, code_mtime):
            log(f"extracting {pdf.name} ...")
            extract_pdf(pdf, txt)
        out[stem] = txt
    return out


def condense(text: str) -> str:
    """Block text -> display text: one logical line per field, single spaces, table cells joined by ' | '."""
    lines = text.split("\n")
    first = lines[0]
    heading = re.sub(r"[\s—\d]", "", first) == "Operatingconditions"
    body = lines[1:] if heading else lines
    out: list[str] = []
    for ln in body:
        if not out or ps.LINE_LABEL_RE.match(ln) or "  " in ln or re.match(r"^(Proceed|Perform)\b", ln) and not out[-1].rstrip().endswith(("and", ",")) \
                or (out and "  " in out[-1] and re.match(r"^[A-Z]", ln) and not re.match(r"^of sample|^sample", ln)):
            out.append(ln.strip())
        else:
            out[-1] += " " + ln.strip()
    return "\n".join(re.sub(r"\s{2,}", " | ", o) if "  " in o else re.sub(r"\s+", " ", o) for o in out)


def build_records(text_files: dict[str, Path], log=print):
    """Segment every file and return (records, events_by_source, header_by_title, stats)."""
    records_by_source: dict[str, list[dict]] = defaultdict(list)
    events_by_source: dict[str, list[dict]] = defaultdict(list)
    headers_by_title: dict[str, dict] = {}
    stats = {"blocks_found": Counter(), "events": Counter()}
    counters: Counter = Counter()
    all_blocks: list[dict] = []
    for stem, source, part, desc, pfx in SOURCES:
        rows = sg.load_rows(text_files[stem].read_text(encoding="utf-8"))
        blocks, events, monos = sg.segment(rows, source=source, file=stem, default_part=part)
        log(f"  {stem}: {len(blocks)} 'Operating conditions' blocks, {len(events)} revision events")
        stats["blocks_found"][stem] = len(blocks)
        # headers of the monographs in this file
        hdrs = []
        for m in monos:
            h = hd.parse_header(m["header_rows"], m["monograph"]) if m["part"] in ("mono", "crude") else {}
            hdrs.append(h)
            if h and (source == "JP18" or h.get("formula") or h.get("cas")):      # main volume first; supplements add new monographs
                headers_by_title.setdefault(mg.norm(m["monograph"]), {**h, "title": m["monograph"]})
        for b in blocks:
            code = PART_CODE[b["part"]]
            counters[(pfx, code)] += 1
            b["block_id"] = f"{pfx}-{code}-{counters[(pfx, code)]:04d}"
            b["header"] = hdrs[b["mono_idx"]] if b["mono_idx"] >= 0 and b["mono_idx"] < len(hdrs) else {}
            all_blocks.append(b)
            records_by_source[source].append(_record(b))
        for e in events:
            events_by_source[source].append(e.__dict__)
    return records_by_source, events_by_source, headers_by_title, all_blocks, stats


def _record(b: dict) -> dict:
    p = ps.parse_block(b)
    sec = SECTION_CANON.get(b["section"], b["section"])
    rec = {
        "block_id": b["block_id"], "source": b["source"], "file": b["file"] + ".pdf", "page": b["page"], "page_printed": b["page_printed"] or None,
        "part": b["part"], "monograph": b["monograph"], "entry": b["entry"], "section": sec,
        "item_no": b["item_no"], "item_title": b["item_title"], "cond_label": b["cond_label"], "sub_label": b.get("sub_label"), "sub_disp": b.get("sub_disp"),
        "event_id": b["event_id"], "technique_hint": b["technique_hint"], "ctx_text": b["ctx_text"],
        "_header": b["header"], "condition_text": condense(b["text"]), "raw_text": b["text"], **p,
    }
    if rec["technique"] is None and b["technique_hint"] in ("LC", "GC"):
        rec["technique"] = b["technique_hint"]
    return rec


def item_label(rec: dict) -> tuple[str | None, str | None]:
    """(item, test_label) e.g. ('(3) Related substances', 'Purity (3) Related substances')."""
    no, title = rec.get("item_no"), rec.get("item_title")
    item = None
    if no is not None and title:
        item = f"({no}) {title}"
    elif no is not None:
        item = f"({no})"
    elif title:
        item = title
    if rec.get("sub_disp"):
        item = f"{item} {rec['sub_disp']}" if item else rec["sub_disp"]
    sec = rec.get("section")
    if rec["part"] in ("gen", "geninfo") and rec.get("entry") and not re.match(r"^9\.41", rec["monograph"]):
        base = f"{rec['entry']}: {sec}" if sec else rec["entry"]
    else:
        base = sec
    label = " ".join(x for x in (base, item) if x) or (rec.get("entry") if rec["part"] in ("gen", "geninfo") else None)
    if label and sec and item and sec in item:
        label = " ".join(x for x in (base.replace(sec, "").strip() or sec, item) if x)
    return item, label


def lookup_bulk(bulk_index: dict, name: str):
    """Bulk monograph for a substance name; preparations are titled with the anhydrous / salt name, bulk monographs often 'X Hydrate'."""
    n = mg.norm(name)
    for key in (n, n + "hydrate", n + "monohydrate", "anhydrous" + n, mg.norm(re.sub(r"(?i)\s+hydrate$", "", name))):
        if key in bulk_index:
            return bulk_index[key]
    return None


def header_fields(rec: dict, headers_by_title: dict, bulk_index: dict) -> dict:
    """Monograph-level header info for a record (JP name, formula, MW, chemical name, CAS) + how its analytes are named."""
    h = rec.pop("_header", {}) or {}
    mono_key = mg.norm(rec["monograph"])
    base = headers_by_title.get(mono_key) or {}          # supplement partial revisions repeat only the title: use the main-volume header
    own_formula = bool(h.get("formula") or base.get("formula"))
    out = {"jp_name": h.get("jp_name") or base.get("jp_name"), "formula": h.get("formula") or base.get("formula"),
           "mw": h.get("mw") if h.get("mw") is not None else base.get("mw"),
           "chemical_name": h.get("chemical_name") or base.get("chemical_name"),
           "cas": list(h.get("cas") or base.get("cas") or []), "substance": h.get("substance") or base.get("substance"),
           "formula_source": "header" if own_formula else None, "header_inherited_from": None,
           "substances": [x["name"] for x in (h.get("substances") or base.get("substances") or [])],
           "alt_names": h.get("alt_names") or base.get("alt_names") or [], "prep_substances": None, "analyte_cas": {}, "bulk": False}
    if not out["formula"] and (h.get("formula_stmt") or base.get("formula_stmt")):       # "(C21H18ClNO6: 415.82)" in the content statement
        out["formula"] = h.get("formula_stmt") or base.get("formula_stmt")
        out["mw"] = h.get("mw_stmt") if h.get("mw_stmt") is not None else base.get("mw_stmt")
        out["formula_source"] = "content statement"
    is_prep = bool(re.search(r"\s+" + hd.FORMS + r"$", rec["monograph"]))
    if rec["part"] == "mono" and not is_prep:                    # bulk drug monograph: analyte = the title, CAS from its header
        out["bulk_name"] = hd.bulk_name(rec["monograph"])
        out["analyte_cas"] = {out["bulk_name"].lower(): list(out["cas"])}
        out["bulk"] = True
    elif rec["part"] == "mono":                                  # preparation: find the bulk monograph(s) behind the title
        base_title = hd.bulk_name(re.sub(r"\s+" + hd.FORMS + r"$", "", rec["monograph"]).strip())
        parts = [base_title] if lookup_bulk(bulk_index, base_title) else re.split(r"\s+and\s+|,\s+", base_title)
        hits = [lookup_bulk(bulk_index, p) for p in parts]
        if base_title != rec["monograph"] and hits and all(hits):
            out["prep_substances"] = [(b["title"].lower(), list(b["header"].get("cas", []))) for b in hits]
            out["analyte_cas"] = dict(out["prep_substances"])
            if len(hits) == 1:
                b = hits[0]
                out.update(cas=list(b["header"].get("cas", [])), chemical_name=b["header"].get("chemical_name"),
                           formula=b["header"].get("formula"), mw=b["header"].get("mw"),
                           formula_source=f"bulk monograph {b['title']}", header_inherited_from=b["title"])
            else:
                out.update(cas=[], chemical_name=None, formula=None, mw=None, formula_source="combination preparation (see analytes)")
        elif len(out["substances"]) > 1:
            out.update(cas=[], chemical_name=None, formula=None, mw=None, formula_source="combination preparation (see analytes)")
        elif out["substance"]:
            bulk = lookup_bulk(bulk_index, out["substance"])
            if bulk and bulk["title"] != rec["monograph"]:
                out["cas"] = list(bulk["header"].get("cas", []))
                out["chemical_name"] = out["chemical_name"] or bulk["header"].get("chemical_name")
                out["analyte_cas"] = {out["substance"].lower(): list(out["cas"])}
                if not out["formula"]:
                    out["formula"], out["mw"], out["formula_source"] = bulk["header"].get("formula"), bulk["header"].get("mw"), f"bulk monograph {bulk['title']}"
                out["header_inherited_from"] = bulk["title"]
    out["monograph_is_bulk"] = bool(out.get("bulk"))
    return out


def bulk_monograph_index(headers_by_title: dict, rec_titles: dict) -> dict:
    """norm(title) -> {'title', 'header'} for bulk-substance monographs (own formula line or own CAS RN)."""
    idx = {}
    for k, h in headers_by_title.items():
        if h.get("formula") or h.get("cas"):
            idx[k] = {"title": h.get("title") or rec_titles.get(k, k), "header": h}
    return idx


def finalize(records_by_source, events_by_source, headers_by_title, log=print):
    """Merge editions, resolve cross-references, add analytes/labels. Returns (final, dropped_nonchrom, superseded, merge_log)."""
    # titles for bulk lookup
    rec_titles = {}
    for recs in records_by_source.values():
        for r in recs:
            rec_titles.setdefault(mg.norm(r["monograph"]), r["monograph"])
    bulk_index = bulk_monograph_index(headers_by_title, rec_titles)
    for recs in records_by_source.values():
        for r in recs:
            r["header_info"] = header_fields(r, headers_by_title, bulk_index)
    # keep only chromatographic blocks (LC / GC / unknown technique of prose blocks), drop NMR / other
    dropped = []
    chrom: dict[str, list[dict]] = {}
    for src, recs in records_by_source.items():
        keep = []
        for r in recs:
            if r["technique"] in ("NMR", "other"):
                dropped.append(r)
            else:
                keep.append(r)
        chrom[src] = keep
    mlog: list[dict] = []
    main = chrom.get("JP18", [])
    supps = [(chrom.get(s, []), events_by_source.get(s, [])) for s in ("JP18 Supplement I", "JP18 Supplement II")]
    current, superseded = mg.apply_supplements(main, supps, mlog)
    mg.resolve_xrefs(current, mlog)
    final = []
    for r in current:
        hi = r.pop("header_info")
        r["item"], r["test_label"] = item_label(r)
        ans, how = an.analytes_for(r, hi)
        r["analytes"], r["analyte_source"] = ans, how
        r["substance_names"] = an.substance_for(r, hi)
        if r.get("_wl_options"):                     # "231 nm for A and B; 254 nm for C": the wavelength stated for this block's analyte
            for a_ in ans:
                tok = an._first_tok(a_)
                hit = next((w for w, label in r["_wl_options"] if tok in label), None)
                if a_ != "related substances" and hit is not None:
                    r["wavelength_nm"] = hit
                    break
        for k in ("jp_name", "formula", "mw", "chemical_name", "cas", "substance", "formula_source", "header_inherited_from", "alt_names", "analyte_cas"):
            r[k] = hi.get(k)
        final.append(r)
    for r in dropped + superseded:
        r.pop("header_info", None)
    return final, dropped, superseded, mlog
