"""Merge the three editions into a 'current JP18' record set and resolve cross-references.

1. `apply_supplements`: Supplement II overrides Supplement I overrides the main volume, matched by
   (monograph, section, item). New monographs are added, deleted monographs removed. Revision instructions
   ("Change the Purity (3) as follows:") were parsed in segment.py; their scope decides which earlier blocks are replaced.
2. `resolve_xrefs`: 'Proceed as directed in the operating conditions in the Assay (3) under X' is resolved to the
   referenced block of the current set and the missing fields are filled in.
"""
from __future__ import annotations

import difflib
import re
from collections import defaultdict

from .parse import XREF_GROUPS, RUN_TIME

KEY_FIELDS = ["detector_type", "wavelength_nm", "column_id_mm", "column_length_mm", "packing_class", "particle_or_film_um",
              "column_temp_c", "elution_mode", "n_mobile_phases", "carrier_gas"]


def norm(s) -> str:
    """Comparison key for names: lower case, no spaces / punctuation ('L- α -Phenethicillin' == 'L-α-phenethicillin')."""
    return re.sub(r"[^0-9a-zα-ωα-ω]+", "", (s or "").lower())


def norm_section(s) -> str:
    n = norm(s)
    return re.sub(r"s$", "", n) if n.startswith("uniformityofdosageunit") else n


# ----------------------------------------------------------------------------------------------------------
# supplements
# ----------------------------------------------------------------------------------------------------------
def _item_match(scope_item, rec_item) -> bool:
    if scope_item is None or rec_item is None:
        return True
    return scope_item == rec_item


def _in_scope(rec: dict, scope: list) -> bool:
    for sec, item in scope:
        if norm_section(sec) == norm_section(rec.get("section")) and _item_match(item, rec.get("item_no")):
            return True
    return False


def apply_supplements(main: list[dict], supps: list[tuple[list[dict], list[dict]]], log: list[dict]):
    """main: records of the main volume. supps: [(records, events)] for Supplement I then II (events as dicts).
    Returns (current_records, superseded_records). Every decision is appended to `log`."""
    current = list(main)
    superseded: list[dict] = []
    for recs, events in supps:
        by_event = defaultdict(list)
        for r in recs:
            by_event[r.get("event_id")].append(r)
        added: list[dict] = []
        handled_blocks = set()
        for ev in events:
            new = by_event.get(ev["event_id"], [])
            for r in new:
                handled_blocks.add(id(r))
            mono_n = norm(ev["monograph"])
            ent_n = norm(ev["entry"]) if ev.get("entry") else None
            act, kind = ev["action"], ev["kind"]
            note = None
            victims: list[dict] = []
            if kind == "container" or (act in ("other", "main") and kind in ("whole", "other")) or ev.get("ambiguous"):
                if new:
                    note = "no replacement scope (instruction ambiguous or absent): kept both"
                    log.append(_log(ev, new, [], "kept_both", note))
                added.extend(new)
                continue
            if kind == "whole":
                victims = [c for c in current if norm(c["monograph"]) == mono_n and (ent_n is None or norm(c.get("entry")) == ent_n)]
                if not victims and act in ("change", "delete"):
                    fz = _fuzzy_victims(current, mono_n, ent_n)
                    if fz:
                        victims, note = fz, "fuzzy name match (spelling differs between editions)"
                if act == "add" and victims:
                    note = "'Add' instruction collides with existing record: replaced"
                decision = {"add": "added", "change": "replaced" if victims else "added", "delete": "deleted"}.get(act, act)
            elif kind == "partial":
                if ev["scope"] and act != "add":               # "Add the following next to X" inserts text; nothing is replaced
                    victims = [c for c in current if norm(c["monograph"]) == mono_n and _in_scope(c, ev["scope"])]
                for sec, item in ev.get("delete_items", []):
                    victims += [c for c in current if norm(c["monograph"]) == mono_n and norm_section(c.get("section")) == norm_section(sec)
                                and c.get("item_no") == item and c not in victims]
                decision = {"add": "added_section", "change": "replaced" if victims else "added", "delete": "deleted"}.get(act, act)
                if not ev["scope"] and new:
                    note = f"instruction names no known section ({ev['raw'][:60]!r}): kept both"
                    log.append(_log(ev, new, [], "kept_both", note))
                    added.extend(new)
                    continue
                for r in new:                                      # block outside the stated scope -> flag
                    if ev["scope"] and not _in_scope(r, ev["scope"]):
                        r.setdefault("qa", []).append("block outside instruction scope")
            else:
                added.extend(new)
                continue
            if act == "delete":
                for v in victims:
                    v["superseded_by"] = f"DELETED by {ev['source']} p{ev['page']}"
                    superseded.append(v)
                current = [c for c in current if c not in victims]
                if victims or new:
                    log.append(_log(ev, new, victims, "deleted", note))
                continue
            if victims:
                for v in victims:
                    v["superseded_by"] = (new[0]["block_id"] if new else f"removed by {ev['source']} p{ev['page']} (no chromatographic block in replacement)")
                    superseded.append(v)
                current = [c for c in current if c not in victims]
                for r in new:
                    r["revised_in"] = ev["source"]
                    r["replaces"] = [v["block_id"] for v in victims]
            if victims or new:
                log.append(_log(ev, new, victims, decision, note))
            added.extend(new)
        # blocks whose event is unknown (should not happen) are added with a flag
        for r in recs:
            if id(r) not in handled_blocks:
                r.setdefault("qa", []).append("block without revision event")
                added.append(r)
        current.extend(added)
    return current, superseded


def _fuzzy_victims(current: list[dict], mono_n: str, ent_n: str | None) -> list[dict]:
    """Spelling differs between editions ('Rhyncophylline' vs 'Rhynchophylline'): closest entry (or monograph) name, unique and clear."""
    if ent_n:
        pool = {norm(c.get("entry")) for c in current if norm(c["monograph"]) == mono_n and c.get("entry")}
        target, field = ent_n, "entry"
    else:
        pool = {norm(c["monograph"]) for c in current}
        target, field = mono_n, "monograph"
    scored = sorted(((difflib.SequenceMatcher(None, target, p).ratio(), p) for p in pool if p != target), reverse=True)
    if scored and scored[0][0] >= 0.93 and (len(scored) == 1 or scored[0][0] - scored[1][0] >= 0.03):
        best = scored[0][1]
        return [c for c in current if norm(c["monograph"] if field == "monograph" else c.get("entry")) == best
                and (field == "monograph" or norm(c["monograph"]) == mono_n)]
    return []


def _log(ev, new, victims, decision, note):
    return {"source": ev["source"], "page": ev["page"], "monograph": ev["monograph"], "entry": ev.get("entry"),
            "instruction": ev["raw"], "scope": ev["scope"], "decision": decision, "note": note,
            "new_blocks": [r["block_id"] for r in new], "replaced_blocks": [v["block_id"] for v in victims]}


# ----------------------------------------------------------------------------------------------------------
# cross-references
# ----------------------------------------------------------------------------------------------------------
SECTION_RE = r"(?P<sec>Assay|Purity|Identification|Dissolution|Uniformity of dosage units?|Content ratio of the active principle|" \
             r"Isomer ratio|Unity of peak|Molecular mass|Viscosity|Optical purity|Composition)"
TARGET_RE = re.compile(r"^(?:the\s+)?" + SECTION_RE + r"(?:\s*\((?P<item>\d+)\))?(?P<sub>(?:\s*\(?(?:[ivx]+|\d)\)\s*)*)"
                       r"(?:\s*\((?P<paren>[^)]*)\))?(?:\s+under\s+(?:the\s+)?(?P<mono>.+))?$", re.I)


def parse_target(t: str | None) -> dict | None:
    if not t:
        return None
    t = re.sub(r"\s+", " ", t).strip(" .")
    m = re.match(r"^operating conditions (\d+)$", t, re.I)
    if m:
        return {"cond": m.group(1)}
    m = TARGET_RE.match(t)
    if m:
        sec = m.group("sec")
        return {"sec": re.sub(r"(?i)unit$", "units", sec), "item": int(m.group("item")) if m.group("item") else None,
                "sub": (m.group("sub") or "").strip() or None, "paren": m.group("paren"), "mono": (m.group("mono") or "").strip() or None}
    m = re.match(r"^(?:(?:\d+\.\s*)?Purity\s+)?(?P<num>\d+(?:\.\d+)*)\.?\s*(?:under\s+(?:the\s+)?(?P<mono>.+))?$", t)
    if m and (m.group("mono") or "." in m.group("num")):
        return {"chapter_item": m.group("num"), "mono": (m.group("mono") or "").strip() or None}
    m = re.match(r"^Procedure\s*(?P<sub>\(?[ivx]+\)|\(\d+\)|\d\))\s*under\s+(?:the\s+)?" + SECTION_RE + r"(?:\s*\((?P<item>\d+)\))?$", t, re.I)
    if m:
        return {"sec": m.group("sec"), "item": int(m.group("item")) if m.group("item") else None, "sub": m.group("sub"), "paren": None, "mono": None}
    m = re.match(r"^(?P<sub>\(?[ivx]+\)|\(\d+\)|\d\))(?:\s.*)?$", t, re.I)
    if m:
        return {"sub_only": m.group("sub")}
    m = re.match(r"^(?P<sec>[A-Za-zα-ω\- ]+?)(?:\s*\((?P<item>\d+)\))?$", t)
    if m:
        return {"sec": m.group("sec"), "item": int(m.group("item")) if m.group("item") else None, "sub": None, "paren": None, "mono": None}
    return None


def _same_scope(a: dict, b: dict) -> bool:
    return norm(a["monograph"]) == norm(b["monograph"]) and norm(a.get("entry")) == norm(b.get("entry"))


def resolve_xrefs(records: list[dict], log: list[dict]) -> None:
    by_mono = defaultdict(list)
    for r in records:
        by_mono[norm(r["monograph"])].append(r)
    mono_titles = {norm(r["monograph"]): r["monograph"] for r in records}

    def find_mono(name: str) -> str | None:
        n = norm(re.sub(r"<[\d.]+>", "", name))
        if n in by_mono:
            return n
        cands = [k for k in by_mono if k.endswith(n)]            # "Crude Drugs Test" -> "5.01 Crude Drugs Test"
        if len(cands) == 1:
            return cands[0]
        cands = [k for k in by_mono if k.startswith(n) and len(k) - len(n) < 4]
        if len(cands) == 1:
            return cands[0]
        close = difflib.get_close_matches(n, list(by_mono), n=2, cutoff=0.9)       # typos in the source: "Ephedrin", "Processed Aconitine Root"
        return close[0] if len(close) == 1 or (close and difflib.SequenceMatcher(None, n, close[0]).ratio() - difflib.SequenceMatcher(None, n, close[1]).ratio() > 0.03) else None

    state: dict[str, str] = {}

    def _sub(x: str | None) -> str | None:
        return re.sub(r"[()\s]", "", x or "").lower() or None

    def candidates(rec: dict, tgt: dict) -> tuple[list[dict], str | None]:
        same = [c for c in by_mono[norm(rec["monograph"])] if _same_scope(c, rec) and c is not rec]
        if "cond" in tgt:
            return [c for c in same if c.get("cond_label") == tgt["cond"]
                    and norm_section(c.get("section")) == norm_section(rec.get("section")) and c.get("item_no") == rec.get("item_no")], None
        if "chapter_item" in tgt:
            mkey = find_mono(tgt["mono"]) if tgt.get("mono") else None
            if not mkey:
                return [], "target chapter not found"
            num = tgt["chapter_item"].rstrip(".")
            return [c for c in by_mono[mkey] if re.match(r"^" + re.escape(num) + r"\.?(\s|$)", c.get("entry") or "")], None
        if "sub_only" in tgt:
            sub = _sub(tgt["sub_only"])
            if sub and sub.isdigit():                              # "(1)" = item (1) of the same test section
                return [c for c in same if norm_section(c.get("section")) == norm_section(rec.get("section"))
                        and c.get("item_no") == int(sub)], None
            pool = [c for c in same if norm_section(c.get("section")) == norm_section(rec.get("section")) and c.get("item_no") == rec.get("item_no")]
            return [c for c in pool if c.get("sub_label") == sub] or pool, None
        sec = tgt["sec"]
        if tgt.get("mono"):
            mkey = find_mono(tgt["mono"])
            if not mkey:
                return [], f"target monograph '{tgt['mono']}' not found"
            pool = by_mono[mkey]
            if mkey != norm(rec["monograph"]) and not any(True for _ in pool):
                return [], "target monograph has no chromatographic block"
        else:
            pool = same
        pool = [c for c in pool if c is not rec and norm_section(c.get("section")) == norm_section(sec)]
        if tgt.get("item") is not None:
            pool = [c for c in pool if c.get("item_no") in (tgt["item"], None)] if len(pool) == 1 else [c for c in pool if c.get("item_no") == tgt["item"]]
        sub = _sub(tgt.get("sub"))
        if sub:
            sub_pool = [c for c in pool if c.get("sub_label") == sub]
            pool = sub_pool or pool
        return pool, None

    def effective(rec: dict, depth=0) -> dict | None:
        """The record whose own fields are complete (follows xref chains)."""
        if rec["block_type"] == "full" or depth > 5:
            return rec
        if rec.get("resolved") and rec.get("resolved_from"):
            return rec
        return rec

    def pick(rec: dict, pool: list[dict]) -> tuple[dict | None, str | None]:
        pool = [c for c in pool if c["technique"] in ("LC", "GC")]
        if not pool:
            return None, "no chromatographic block at target"
        full = [c for c in pool if c["block_type"] == "full"]
        use = full or pool
        if len(use) == 1:
            return use[0], None
        sig = {tuple(str(c.get(k)) for k in KEY_FIELDS) for c in use}
        if len(sig) == 1:
            return use[0], "multiple identical candidates"
        return None, f"ambiguous: {len(use)} candidate blocks ({', '.join(c['block_id'] for c in use[:4])})"

    order = [r for r in records if r["block_type"] in ("full_xref", "partial_xref", "prose")]
    # resolve in passes so that chains (xref -> xref -> full) settle
    for _pass in range(4):
        progress = False
        for rec in order:
            if rec.get("resolved") or rec.get("_unres_final"):
                continue
            if rec["block_type"] == "prose":
                _resolve_prose(rec, by_mono, log)
                continue
            tgt = parse_target(rec.get("xref_target"))
            if tgt is None:
                rec["resolved"], rec["resolve_note"] = False, "cross-reference target not understood"
                rec["_unres_final"] = True
                continue
            pool, err = candidates(rec, tgt)
            if err:
                rec["resolved"], rec["resolve_note"] = False, err
                rec["_unres_final"] = True
                continue
            if not pool and not tgt.get("mono") and "sec" in tgt and rec.get("entry") is None:
                rec["resolved"], rec["resolve_note"] = False, f"target '{rec['xref_target']}' not found in monograph"
                rec["_unres_final"] = True
                continue
            tr, err = pick(rec, pool)
            if tr is None:
                rec["resolved"], rec["resolve_note"] = False, err
                rec["_unres_final"] = True
                continue
            if tr["block_type"] != "full" and not tr.get("resolved"):
                if _pass < 3:
                    continue                                        # wait for the chain's end
                rec["resolved"], rec["resolve_note"] = False, "target is itself an unresolved cross-reference"
                rec["_unres_final"] = True
                continue
            _fill(rec, tr)
            progress = True
        if not progress:
            break
    for rec in order:
        if rec["block_type"] != "prose" and not rec.get("resolved"):
            rec["resolved"] = False
            rec.setdefault("resolve_note", "target not resolved")
        rec.pop("_unres_final", None)


def _fill(rec: dict, tr: dict) -> None:
    groups = sorted(set(rec.get("xref_fields") or []))
    inherited = []
    for g in groups:
        for f in XREF_GROUPS.get(g, []):
            if rec.get(f) is None and tr.get(f) is not None:
                rec[f] = tr[f]
                if not f.startswith("_"):
                    inherited.append(f)
    if rec.get("xref_time_span"):
        for f in RUN_TIME:
            if rec.get(f) is None and tr.get(f) is not None:
                rec[f] = tr[f]
                inherited.append(f)
    if rec["technique"] is None or rec["technique"] == "LC" and tr["technique"] == "GC":
        rec["technique"] = tr["technique"]
    rec["resolved"] = True
    rec["resolved_from"] = tr["block_id"]
    rec["inherited_fields"] = inherited


def _resolve_prose(rec: dict, by_mono, log) -> None:
    """'Operating conditions and system suitability fundamentally follow the procedure A.'"""
    t = rec.get("xref_target")
    if not t:
        rec["resolved"], rec["resolve_note"] = None, None
        return
    m = re.match(r"^Procedure ([A-Z])$", t, re.I)
    pool = by_mono[norm(rec["monograph"])]
    prefix = re.match(r"^(\d+\.\d+)\.", rec.get("entry") or "")
    cands = [c for c in pool if c is not rec and re.search(r"Procedure " + (m.group(1) if m else "") + r"$", c.get("entry") or "")
             and c["block_type"] == "full" and (not prefix or (c.get("entry") or "").startswith(prefix.group(1) + "."))]
    if len(cands) == 1:
        tr = cands[0]
        rec["xref_fields"] = ["column", "column_temperature", "detector", "flow_rate", "mobile_phase"]
        _fill(rec, tr)
        rec["technique"] = tr["technique"]
    else:
        rec["resolved"], rec["resolve_note"] = False, f"procedure '{t}': {len(cands)} candidate blocks"
        rec["_unres_final"] = True
