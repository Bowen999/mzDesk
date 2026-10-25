"""python -m mzdesk.methods.jp_build <pdf folder> [--work DIR] [--out DIR] [--offline] [--skip-mapping] [--force]

Rebuilds mzdesk/data/known_methods/jp18_methods.jsonl.gz from the official MHLW PDFs (see README.md in that folder).
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import Counter
from datetime import date
from pathlib import Path

from . import build as bd
from . import qa as qa_mod
from . import resolve as rs

DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "known_methods"

OUT_FIELDS = [
    "block_id", "source", "file", "page", "page_printed",
    "monograph", "entry", "section", "item", "item_no", "test_label", "kind",
    "jp_name", "formula", "mw", "cas", "chemical_name", "analytes", "analyte_source",
    "technique", "block_type", "xref_target", "xref_fields", "xref_time_span", "xref_except", "resolved", "resolved_from",
    "inherited_fields", "resolve_note",
    "detector_type", "wavelength_nm", "column_count", "column_material", "column_id_mm", "column_length_mm", "packing_class",
    "particle_or_film_um", "column_temp_c", "injection_temp_c", "detector_temp_c", "elution_mode", "n_mobile_phases",
    "mobile_phase_ph", "organic_solvents", "gradient_table", "carrier_gas", "split_ratio", "flow_rate_value", "flow_rate_unit",
    "flow_target_compound", "flow_target_rt_min", "run_time_min", "run_time_rt_multiple",
    "cid", "smiles", "inchikey", "parent_inchikey", "map_method", "map_reason", "analyte_structures", "substances",
    "revised_in", "replaces", "cond_label", "notes", "qa", "condition_text",
]


def to_output(r: dict) -> dict:
    out = {k: r.get(k) for k in OUT_FIELDS if k != "kind"}
    out["kind"] = bd.PART_NAME.get(r["part"], r["part"])
    if r["part"] == "gen" and r["monograph"].startswith("9.41"):
        out["kind"] = "reagent / test solution"
    out["resolved"] = r.get("resolved") if r["block_type"] in ("full_xref", "partial_xref") else None
    if r["block_type"] == "prose":
        out["resolved"] = r.get("resolved")
    for k in ("analytes", "qa", "cas", "inherited_fields", "replaces", "substances"):
        out[k] = out[k] or []
    return {k: out[k] for k in OUT_FIELDS}


def write_jsonl_gz(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9) as gz:
        for row in rows:
            gz.write((json.dumps(row, ensure_ascii=False, separators=(",", ":")) + "\n").encode("utf-8"))


def write_jsonl(path: Path, rows: list[dict]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for row in rows:
            f.write(json.dumps(row, ensure_ascii=False, default=str) + "\n")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m mzdesk.methods.jp_build", description=__doc__)
    ap.add_argument("pdf_folder", help="folder with the MHLW JP18 PDFs (000945683, 000904440, 000904449, 000904450, 000904454, 001109611, 001354067)")
    ap.add_argument("--work", help="scratch folder for text cache, debug files, PubChem cache, QA report (default: <pdf folder>/_work)")
    ap.add_argument("--out", help=f"output folder (default: {DATA_DIR})")
    ap.add_argument("--offline", action="store_true", help="use only the cached PubChem responses")
    ap.add_argument("--skip-mapping", action="store_true", help="do not map analytes to structures")
    ap.add_argument("--force", action="store_true", help="re-extract the text from the PDFs")
    a = ap.parse_args(argv)

    folder = Path(a.pdf_folder).expanduser().resolve()
    work = Path(a.work).expanduser().resolve() if a.work else folder / "_work"
    out_dir = Path(a.out).expanduser().resolve() if a.out else DATA_DIR
    work.mkdir(parents=True, exist_ok=True)

    print("1/5 extracting text ...")
    text_files = bd.extract_all(folder, work, force=a.force)
    print("2/5 segmenting blocks and parsing conditions ...")
    rbs, ebs, hdrs, all_blocks, stats = bd.build_records(text_files)
    print("3/5 merging supplements, resolving cross-references ...")
    final, dropped, superseded, mlog = bd.finalize(rbs, ebs, hdrs)
    qa_mod.annotate(final)
    map_stats, unmapped = None, None
    if not a.skip_mapping:
        print("4/5 mapping analytes to structures (PubChem, cached) ...")
        cache = rs.Cache(work / "pubchem_cache.json", offline=a.offline)
        map_stats = rs.map_records(final, cache)
        unmapped = unmapped_list(final)
        with open(work / "unmapped_analytes.tsv", "w", encoding="utf-8") as f:
            f.write("analyte\tmonograph\tn_records\treason\n")
            for u in unmapped:
                f.write(f"{u['analyte']}\t{u['monograph']}\t{u['n']}\t{u['reason']}\n")
    else:
        for r in final:
            r.update(cid=None, smiles=None, inchikey=None, parent_inchikey=None, map_method="none", map_reason="mapping skipped",
                     analyte_structures=[], substances=[])
    print("5/5 writing outputs ...")
    final.sort(key=lambda r: (bd.SOURCE_ORDER[r["source"]], r["block_id"]))
    rows = [to_output(r) for r in final]
    write_jsonl_gz(out_dir / "jp18_methods.jsonl.gz", rows)
    (out_dir / "README.md").write_text(readme(final, stats, map_stats), encoding="utf-8")
    write_jsonl(work / "blocks_all.jsonl", [{k: v for k, v in b.items() if k != "header"} | {"header": b["header"].get("substance")} for b in all_blocks])
    write_jsonl(work / "superseded.jsonl", [{k: v for k, v in r.items() if k in ("block_id", "source", "page", "monograph", "entry", "section",
                                                                            "item_no", "superseded_by", "condition_text")} for r in superseded])
    (work / "merge_log.json").write_text(json.dumps(mlog, ensure_ascii=False, indent=1), encoding="utf-8")
    (work / "QA_REPORT.md").write_text(qa_mod.report(final, dropped, superseded, mlog, stats, map_stats, unmapped, all_blocks), encoding="utf-8")
    size = (out_dir / "jp18_methods.jsonl.gz").stat().st_size
    print(f"done: {len(rows)} records -> {out_dir / 'jp18_methods.jsonl.gz'} ({size / 1024:.0f} KiB); QA report: {work / 'QA_REPORT.md'}")
    return 0


def unmapped_list(final: list[dict]) -> list[dict]:
    seen: dict[str, dict] = {}
    for r in final:
        for s in r.get("analyte_structures", []):
            if s["cid"] is None:
                u = seen.setdefault(s["name"], {"analyte": s["name"], "monograph": r["monograph"], "n": 0, "reason": s["reason"] or ""})
                u["n"] += 1
    return sorted(seen.values(), key=lambda u: (u["reason"], u["analyte"]))


def readme(final: list[dict], stats: dict, map_stats: dict | None) -> str:
    from .build import SOURCES, SOURCE_URL
    n_src = Counter(r["source"] for r in final)
    tech = Counter(r["technique"] for r in final)
    mapped = sum(1 for r in final if r["cid"])
    fields = "\n".join(f"- `{f}`" for f in OUT_FIELDS)
    src_rows = "\n".join(f"| {s} | {desc} | {SOURCE_URL[stem]} |" for stem, s, _p, desc, _x in SOURCES)
    return f"""# Known methods: JP18 chromatography conditions

`jp18_methods.jsonl.gz` has one JSON object per "Operating conditions" block of the Japanese Pharmacopoeia, 18th edition (English
version) that describes a liquid-chromatography (LC) or gas-chromatography (GC) method. It is the dataset behind the
"Known methods" page (compound -> published methods). Rebuild: `python -m mzdesk.methods.jp_build <folder with the PDFs>`.

* Edition: JP18 (Ministerial Notification No. 220, 2021), English version, with Supplement I (Notification No. 355, 2022) and Supplement II
  (Notification No. 238, 2024) merged into a "current" set. Errata (正誤表) published later by MHLW are NOT applied.
* Build date: {date.today().isoformat()}.  Records: {len(final)} ({n_src.get("JP18", 0)} from the main volume, {n_src.get("JP18 Supplement I", 0)} from
  Supplement I, {n_src.get("JP18 Supplement II", 0)} from Supplement II); LC {tech.get("LC", 0)}, GC {tech.get("GC", 0)}.
  NMR (qNMR) blocks and non-chromatographic blocks (viscometer, capillary electrophoresis, thermogravimetry, polarimeter) are dropped.
* Structure mapping: {mapped} records have a PubChem structure for their primary analyte (see `map_method`).

## Sources (official MHLW PDFs)

| edition | content | URL |
|---|---|---|
{src_rows}

Index page: https://www.mhlw.go.jp/stf/seisakunitsuite/bunya/0000066530.html

## Licence

The Japanese Pharmacopoeia is a Ministerial Notification (厚生労働省告示). Under Article 13 of the Japanese Copyright Act, laws, notifications
and the like issued by national or local authorities are not protected by copyright. This dataset therefore stores the extracted
parameters plus the text of each condition block (`condition_text`) and page references (`file`, `page`, `page_printed`); the PDFs are not redistributed.
The English version is a convenience translation: where it differs from the Japanese original, the Japanese is authoritative. Check the
printed page before relying on a value.

## Records

One object per block. Cross-referenced blocks ("Proceed as directed in the operating conditions in the Assay") are resolved: the missing fields are
copied from the referenced block (`resolved`, `resolved_from`, `inherited_fields`); the original reference stays in `xref_target` / `xref_fields`.
Supplement II overrides Supplement I overrides the main volume (matched by monograph, section and item); superseded blocks are not in this file
(`revised_in`, `replaces` mark the winning records). `null` = the block does not state the value.

Fields (in order):

{fields}

Notes on selected fields

* `file`, `page`, `page_printed`: PDF file name, page index in that PDF (1-based) and the page number printed in the running header (differs by a constant per file).
* `kind`: monograph / crude drug / extract / general test / reagent / general information.
* `item`, `item_no`, `test_label`: e.g. `(3) Related substances`, 3, `Purity (3) Related substances`. Sub-procedure headings are kept when they name the block:
  `(3) Related substances (iii)`, `(ii) Determination of molecular mass`; long titles that wrap over two printed lines are joined.
* `cas`: CAS RN(s) in brackets in the monograph header; preparations (tablets ...) inherit the CAS of the bulk-drug monograph.
* `analytes`: lower-case names of what the test MEASURES. Assay / dissolution / content uniformity: the monograph substance (bulk monograph = its title, a
  preparation = the bulk monograph it is made from, a combination preparation = the active(s) the test text talks about). A Purity item whose title names
  compounds ("Methanol", "Ethylene glycol, diethylene glycol and related substances", "1-Vinyl-2-pyrrolidone") measures those compounds and not the substance
  (the substance is kept only when the test compares the impurity peak with the substance's own peak). "Residual solvents" take the solvents named in the
  quantity formulae. `related substances` is added for related-substance tests (monographs and reagent entries). Crude drugs / extracts: marker compounds from
  the item title, else the quantity formula or the peak statement. Reagent entries: the reagent compound. Unnamed impurity classes ("Diastereomer",
  "Multimers") keep the substance.
* `detector_type`: UV, PDA, visible, fluorescence, RI, conductivity, MS, FID, TCD, ECD, other. `packing_class`: C18, C8, C4, phenyl, CN, NH2, silica, diol, PFP, anion-exchange, cation-exchange, SEC, chiral, polymer-RP, GC-PEG, GC-6%cyanopropylphenyl, GC-5%phenyl, GC-dimethylpolysiloxane, GC-other, other. `silica` is only the plain silica gel for liquid chromatography; any bonded or modified gel (carbamoyl, trimethylsilanized,
  sulfonamide / hexadecylsilanized, ...) is `other`, albumin / glycoprotein / cellulose / cyclodextrin phases are `chiral`, "sulfo" alone is not a cation exchanger.
* `column_temp_c`: number, `"program"` (temperature programme) or `"room"`. `elution_mode`: isocratic / gradient (a gradient table, or several mobile phases switched).
* `gradient_table`: rows `[t_start, t_end, A_start, A_end, B_start, B_end]` in min / vol%; further mobile phases append `[C_start, C_end]` ...
* `flow_rate_value` + `flow_rate_unit` (mL/min or cm/s) when a flow is stated; `flow_target_compound` + `flow_target_rt_min` when the text says "adjust so that the retention time of X is about N minutes" (midpoint of a range).
* `xref_time_span`: the cross-reference also covers the time span of measurement. `xref_except`: text after "except".
* `cid`, `smiles` (isomeric), `inchikey`, `parent_inchikey` (largest fragment, neutralised; use it for matching), `map_method` (`cas`, `name`, `name-stripped`, `none`) describe the primary analyte; `analyte_structures` has one entry per analyte. `substances` (drug monographs only) lists the structure of the monograph substance the
  sample consists of, which differs from `analytes` for impurity tests (the Captopril Purity (4) test measures captopril disulfide; `substances` still names captopril). Biologics, polymers and mixtures have no structure (`map_reason`).
* `qa`: automatic QA flags (missing fields, out-of-range values, unresolved references ...); see QA_REPORT.md in the build work folder.

## Known limitations

* Only blocks introduced by an "Operating conditions" heading are covered. A few monographs print the conditions under another heading (the Assay of
  Propiverine Hydrochloride has "System suitability" there), so references to them stay unresolved (`resolved: false`, listed in the QA report).
* The mobile-phase composition, buffer concentrations, injection volume and system-suitability limits are not parsed; read `condition_text`.
* Alternatives take the first value ("0.32 mm (or 0.53 mm)"); ranges of retention time or flow take the midpoint; for a guard / pre-column plus a
  separation column the analytical column is reported (`notes` says so). Stepwise switching of several mobile phases counts as `gradient`.
* Cross-references that name a sub-procedure only by its position (e.g. "prescription 1)", "Procedure (ii)") are resolved when the position matches a
  block of the same test, otherwise left unresolved rather than guessed.
* Text comes from the PDF text layer; glyph fixes (percent, degree, prime, Greek letters, fi / fl ligatures) are in `extract.py`. Printed page numbers are
  those of the running header.
"""


if __name__ == "__main__":
    sys.exit(main())
