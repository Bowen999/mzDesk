"""Automatic QA: flag suspicious records and write QA_REPORT.md (counts, missing fields, out-of-range values,
unresolved cross-references, supplement merge decisions, structure-mapping statistics)."""
from __future__ import annotations

from collections import Counter, defaultdict
from datetime import date

RANGES_LC = {  # field -> (low, high, label)
    "wavelength_nm": (190, 800, "wavelength 190-800 nm"),
    "column_id_mm": (1, 10, "column ID 1-10 mm"),
    "column_length_mm": (30, 300, "column length 30-300 mm"),
    "particle_or_film_um": (1.5, 10, "particle 1.5-10 um"),
    "column_temp_c": (5, 80, "column temperature 5-80 C"),
}
RANGES_GC = {"column_length_mm": (1000, 100000, "GC column length 1-100 m")}
FILE_LABEL = {"000945683": "General Notices + General Tests", "000904440": "Official Monographs A-L",
              "000904449": "Official Monographs M-Z", "000904450": "Crude Drugs and Related Drugs",
              "000904454": "General Information", "001109611": "Supplement I", "001354067": "Supplement II"}


def annotate(records: list[dict]) -> None:
    """Add `qa` (list of reasons) to every record."""
    for r in records:
        q = r.setdefault("qa", [])
        tech = r["technique"]
        if tech == "LC" and r["block_type"] in ("full", "partial_xref", "full_xref") and (r["block_type"] == "full" or r.get("resolved")):
            miss = []
            if r["detector_type"] is None:
                miss.append("detector")
            if r["column_count"] is None:
                miss.append("column")
            if r["n_mobile_phases"] is None:
                miss.append("mobile phase")
            if r["flow_rate_value"] is None and r["flow_target_rt_min"] is None and "programmed flow" not in (r.get("notes") or ""):
                miss.append("flow (value or target)")
            if miss:
                q.append("missing " + ", ".join(miss))
        if tech == "GC" and r["block_type"] in ("full", "partial_xref", "full_xref") and (r["block_type"] == "full" or r.get("resolved")):
            miss = [n for n, k in (("detector", "detector_type"), ("column", "column_count"), ("carrier gas", "carrier_gas")) if r[k] is None]
            if r["flow_rate_value"] is None and r["flow_target_rt_min"] is None:
                miss.append("flow (value or target)")
            if miss:
                q.append("missing " + ", ".join(miss))
        ranges = RANGES_LC if tech == "LC" else RANGES_GC if tech == "GC" else {}
        for f, (lo, hi, label) in ranges.items():
            v = r.get(f)
            if isinstance(v, (int, float)) and not isinstance(v, bool) and not (lo <= v <= hi):
                q.append(f"out of range: {label} (value {v})")
        if tech == "LC" and r.get("flow_rate_unit") == "mL/min" and r.get("flow_rate_value") is not None \
                and not (0.1 <= r["flow_rate_value"] <= 3):
            q.append(f"out of range: LC flow 0.1-3 mL/min (value {r['flow_rate_value']})")
        if r.get("mobile_phase_ph") is not None and not (0 <= r["mobile_phase_ph"] <= 14):
            q.append(f"out of range: pH {r['mobile_phase_ph']}")
        if r["block_type"] in ("full_xref", "partial_xref") and r.get("resolved") is False:
            q.append("unresolved cross-reference: " + (r.get("resolve_note") or "?"))
        if r["block_type"] == "prose" and r.get("resolved") is False:
            q.append("unresolved prose reference: " + (r.get("resolve_note") or "?"))
        if not r.get("analytes"):
            q.append("no analyte identified")
        if tech is None:
            q.append("technique unknown")
        if r.get("resolved") and r.get("resolve_note") and "multiple identical" in r["resolve_note"]:
            q.append("xref resolved among identical candidates")


def _table(rows, header):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join(["---"] * len(header)) + "|"]
    for row in rows:
        out.append("| " + " | ".join(str(c) for c in row) + " |")
    return "\n".join(out)


def report(final: list[dict], dropped: list[dict], superseded: list[dict], mlog: list[dict], stats: dict, map_stats: dict | None,
           unmapped: list[dict] | None, all_blocks: list[dict]) -> str:
    L: list[str] = [f"# JP18 chromatography methods: QA report", "", f"Build date: {date.today().isoformat()}", ""]
    # ---- counts
    L += ["## 1. Blocks found, kept and dropped", ""]
    found = stats["blocks_found"]
    dropped_by = defaultdict(Counter)
    for d in dropped:
        dropped_by[d["file"][:-4]][d["technique"]] += 1
    final_by = Counter(r["file"][:-4] for r in final)
    sup_by = Counter(r["file"][:-4] for r in superseded)
    rows = []
    tot = Counter()
    for stem, label in FILE_LABEL.items():
        nm, no = dropped_by[stem].get("NMR", 0), dropped_by[stem].get("other", 0)
        kept = final_by[stem] + sup_by[stem]
        rows.append((stem, label, found.get(stem, 0), nm, no, kept, sup_by[stem], final_by[stem]))
        tot.update({"f": found.get(stem, 0), "n": nm, "o": no, "k": kept, "s": sup_by[stem], "t": final_by[stem]})
    rows.append(("**total**", "", tot["f"], tot["n"], tot["o"], tot["k"], tot["s"], tot["t"]))
    L += [_table(rows, ["file", "content", "'Operating conditions' blocks found", "dropped: NMR", "dropped: other (viscometer, CE, TG, polarimeter ...)",
                        "LC/GC kept (before merge)", "superseded by a later edition", "in the final set"]), ""]
    L += [f"Superseded or deleted by a supplement (not in the final set; listed in superseded.jsonl): {len(superseded)}", "",
          f"Final records (current JP18 = main volume + Supplement I + II): **{len(final)}**", ""]
    L += ["### By source and technique", "", _table(sorted(Counter((r["source"], r["technique"]) for r in final).items()) and
                                                   [(s, t, n) for (s, t), n in sorted(Counter((r["source"], r["technique"]) for r in final).items())],
                                                   ["source", "technique", "records"]), ""]
    ft = Counter((r["file"][:-4], r["technique"]) for r in final)
    L += ["### By file and technique (final set)", "", _table([(stem, label, ft[(stem, "LC")], ft[(stem, "GC")]) for stem, label in FILE_LABEL.items()],
                                                            ["file", "content", "LC", "GC"]), ""]
    L += ["### By section", "", _table(Counter(r["section"] or "(none)" for r in final).most_common(40), ["section", "records"]), ""]
    L += ["### By part and block type", "", _table([(p, b, n) for (p, b), n in sorted(Counter((r["part"], r["block_type"]) for r in final).items())],
                                                  ["part", "block type", "records"]), ""]
    # ---- flags
    flagged = [r for r in final if r.get("qa")]
    reasons = Counter()
    for r in flagged:
        for q in r["qa"]:
            reasons[q.split(":")[0].split(" (")[0] if q.startswith("out of range") else q.split(":")[0]] += 1
    L += ["## 2. Flag summary", "", f"Records with at least one flag: {len(flagged)} of {len(final)}", "",
          _table(reasons.most_common(), ["reason", "records"]), ""]
    groups = defaultdict(list)
    for r in flagged:
        for q in r["qa"]:
            key = q.split(":")[0] if not q.startswith("out of range") else "out of range"
            groups[key].append((r, q))
    for key in sorted(groups):
        L += [f"### {key} ({len(groups[key])})", ""]
        L += [_table([(r["block_id"], r["page_printed"] or r["page"], f"{r['monograph'][:40]}{' / ' + r['entry'] if r.get('entry') else ''}",
                       r["test_label"] or "", q) for r, q in groups[key]], ["block_id", "page", "monograph", "test", "reason"]), ""]
    # ---- xrefs
    xr = [r for r in final if r["block_type"] in ("full_xref", "partial_xref")]
    L += ["## 3. Cross-references", "", f"Blocks with a cross-reference: {len(xr)}; resolved: {sum(1 for r in xr if r.get('resolved'))}; "
          f"unresolved: {sum(1 for r in xr if not r.get('resolved'))}", ""]
    bad = [r for r in xr if not r.get("resolved")]
    if bad:
        L += [_table([(r["block_id"], r["page_printed"] or r["page"], r["monograph"][:40], r["xref_target"], r.get("resolve_note")) for r in bad],
                     ["block_id", "page", "monograph", "xref target", "reason"]), ""]
    prose = [r for r in final if r["block_type"] == "prose"]
    if prose:
        L += [f"Prose references ('operating conditions fundamentally follow the procedure A'): {len(prose)}; resolved: "
              f"{sum(1 for r in prose if r.get('resolved'))}", ""]
    # ---- merge
    L += ["## 4. Supplement merge decisions", "", _table(Counter(m["decision"] for m in mlog).most_common(), ["decision", "events"]), ""]
    L += ["Each decision: supplement, page, monograph (entry), instruction, replaced blocks -> new blocks.", ""]
    L += [_table([(m["source"].replace("JP18 ", ""), m["page"], (m["monograph"][:36] + (f" / {m['entry']}" if m.get("entry") else "")),
                   m["instruction"][:60], m["decision"], ", ".join(m["replaced_blocks"]) or "-", ", ".join(m["new_blocks"]) or "-",
                   m.get("note") or "") for m in mlog],
                 ["source", "page", "monograph", "instruction", "decision", "replaced", "new", "note"]), ""]
    kept = [m for m in mlog if m["decision"] == "kept_both"]
    L += [f"Ambiguous instructions (kept both, flagged): {len(kept)}", ""]
    # ---- mapping
    L += ["## 5. Structure mapping", ""]
    if map_stats:
        by_method = Counter()
        for r in final:
            by_method[r["map_method"]] += 1
        L += [f"PubChem requests this run: {map_stats.get('requests')}; network errors: {map_stats.get('errors')}; "
              f"distinct analyte queries: {map_stats.get('unique')}", "",
              "Records by mapping method of the primary analyte:", "", _table(by_method.most_common(), ["map_method", "records"]), ""]
        slot = Counter()
        for r in final:
            for s in r.get("analyte_structures", []):
                slot[s["map_method"]] += 1
        L += ["Analyte slots by method:", "", _table(slot.most_common(), ["map_method", "analyte slots"]), ""]
        if unmapped is not None:
            reasons = Counter(u["reason"] for u in unmapped)
            L += [f"Distinct unmapped analytes: {len(unmapped)} (list: unmapped_analytes.tsv)", "",
                  _table(reasons.most_common(15), ["reason", "analytes"]), ""]
    else:
        L += ["Structure mapping was not run.", ""]
    return "\n".join(L) + "\n"
