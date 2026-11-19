"""Cut the cleaned JP text into one record per 'Operating conditions' occurrence, with monograph / entry /
test-section context, page numbers, and (for supplements) the revision instruction that governs the block.

Text conventions come from extract.py: '### ' titles, '**bold**', '@@I ' supplement instruction lines,
'@@PAGE n [printed]' page markers.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

# ----------------------------------------------------------------------------------------------------------
# vocabulary
# ----------------------------------------------------------------------------------------------------------
TESTS = ["Identification", "Purity", "Assay", "Dissolution", "Uniformity of dosage units", "Content uniformity",
         "Isomer ratio", "Optical purity", "Residual solvents", "Unity of peak", "Content", "Particle size",
         "Composition", "Molecular mass", "Peptide map", "Water", "Related substances", "Enantiomer",
         "Ratio", "Specific activity", "Sialic acid", "Amino acid", "Description", "Others", "Containers and storage",
         "Loss on drying", "Residue on ignition", "Optical rotation", "Melting point", "Absorbance", "pH",
         "Total ash", "Acid-insoluble ash", "Extract content", "Essential oil content", "Disintegration",
         "Foreign insoluble matter", "Insoluble particulate matter", "Sterility", "Bacterial endotoxins",
         "Method of preparation", "Microbial limit", "Heavy metals", "Arsenic", "Content ratio of the active principle",
         "Assay of", "Purity test", "Impurities", "Related substance", "Chiral purity", "Stability"]
_TEST_ALT = "|".join(sorted(map(re.escape, TESTS), key=len, reverse=True))
TEST_RE = re.compile(r"^(?:\*\*)?(" + _TEST_ALT + r")(?:\*\*)?(?=[\s—<(*]|$)")
INLINE_TEST_RE = re.compile(r"\.\s+(" + _TEST_ALT + r")\s*(?:<[\d.]+>)?\s*(?:\(\d+\))?\s*[^—.]{0,80}—")
ITEM_RE = re.compile(r"(?:^|\s)(?:\*\*)?\((\d+)\)(?:\*\*)?\s+([A-Za-z][^—]{1,80}?)—")
FIELD_RE = re.compile(r"^(Detector|Column|Mobile phase|Flowing of mobile phase|Flow rate|Time span|Injection|Carrier gas|"
                      r"Split ratio|Apparatus|Measuring|Reaction|Reagent|Headspace|Equilibration|Temperature|"
                      r"Spectral|Nucleus|Gradient|Post-column|Excitation|Fluorescence|Rate of|Detection sensitivity|"
                      r"Sample|Make|Selection of column|Ion|Mass|Gas|Sensitivity|Reaction coil)\b[^:]{0,110}:")
TABLE_RE = re.compile(r"^(Time after injection|of sample \(min\)|\(min\)|Mobile phase [A-D]\b(?!:)|\(vol%\)|"
                      r"Time \(min\)|Time\b.*\(min\)|Mobile phase|[\d.]+\s*[–-]\s*[\d.]+\s|[\d.]+\s*→|\d+\s{2,}|"
                      r"[\d.]+\s*[–-]\s*[\d.]+\s*$)")
STOP_RE = re.compile(r"^(System [Ss]uitability|System repeatability|System performance|Test for required detectability|### |\*\*|@@I)")
PARA_RE = re.compile(r"^(\d{1,2}\)\s)")
# a wrapped label line such as "Detector, column, column temperature, mobile phase, and" (no colon yet)
LABELISH_RE = re.compile(r"^(?:Detector|Column|Mobile phase|Flow rate|Flowing|Time span)[A-Za-z ,()\-/']*(?:,|\band)$")
GENERIC_FIELD_RE = re.compile(r"^[A-Z][a-z]+(?: [a-z0-9()]+){0,5}:\s")
LEAD_BOLD_RE = re.compile(r"^((?:\*\*[^*]+\*\*\s*(?:[^\s*]{1,3}(?=\s*\*\*)\s*)?)+)")
SUB_RE = re.compile(r"^(?:\(([ivx]{1,4})\)|([ivx]{1,4})\)|(\d{1,2})\))\s")
OPCOND_RE = re.compile(r"^(?:\*\*)?Operating conditions\b")
CHAPTER_RE = re.compile(r"^\d+\.\d+\s")
PART_TITLES = {"Official Monographs": "mono", "Crude Drugs and Related Drugs": "crude", "General Information": "geninfo"}
KEEP_HYPHEN = {"n", "sec", "tert", "cis", "trans", "o", "m", "p", "dl", "d", "l", "meso", "endo", "exo", "syn", "anti",
               "threo", "erythro", "iso", "non", "pre", "post", "non-", "co", "re", "x", "ex"}

# ----------------------------------------------------------------------------------------------------------
# row helpers
# ----------------------------------------------------------------------------------------------------------
MARKS_RE = re.compile(r"[◆◇]")
# "N **-Acetylneuraminic acid** C11H19NO9": the italic prefix of a reagent name sits outside the bold run
ITALIC_PREFIX_RE = re.compile(r"^((?:[A-Za-z\u0370-\u03ff]{1,3}|\([A-Za-z]\))(?:,[A-Za-z])?(?: \*\*[^*]{0,3}\*\* [A-Za-z\u0370-\u03ff\u2032]{1,3})*)"
                              r" \*\*([-(\[,][^*]*)\*\*")


@dataclass
class Row:
    page: int
    printed: str
    text: str
    kind: str = "T"      # T text, H title, I instruction


def load_rows(text: str) -> list[Row]:
    rows: list[Row] = []
    page, printed = 0, ""
    for ln in text.split("\n"):
        if ln.startswith("@@PAGE "):
            parts = ln.split()
            page = int(parts[1])
            printed = parts[2] if len(parts) > 2 else ""
            continue
        if not ln.strip():
            continue
        if ln.startswith("@@I "):
            rows.append(Row(page, printed, ln[4:].strip(), "I"))
        elif ln.startswith("### "):
            rows.append(Row(page, printed, ln[4:].strip(), "H"))
        else:
            if "◆" in ln or "◇" in ln:
                ln = MARKS_RE.sub("", ln).replace("  ", " ").strip()
            ln = ITALIC_PREFIX_RE.sub(lambda m: "**" + re.sub(r" ?\*\* ?", "", m.group(1)) + m.group(2) + "**", ln, count=1)
            rows.append(Row(page, printed, ln))
    return merge_instructions(dehyphen(rows))


def merge_instructions(rows: list[Row]) -> list[Row]:
    out: list[Row] = []
    for r in rows:
        if r.kind == "I" and out and out[-1].kind == "I" and out[-1].page == r.page \
                and not (re.match(r"^(Add|Change|Delete|Remove)\b", r.text) and out[-1].text.rstrip().endswith((":", "."))):
            prev = out[-1].text
            out[-1].text = (prev[:-1] + r.text) if (prev.endswith("-") and r.text[:1].islower()) else prev + " " + r.text
        else:
            out.append(r)
    return out


def _strip_marks(s: str) -> str:
    return s.replace("**", "")


def dehyphen(rows: list[Row]) -> list[Row]:
    """Join words split across lines ('diame-' + 'ter'); keep real hyphens (2-aminopyridine, (E)-ferulic, tert-butyl)."""
    out: list[Row] = []
    for r in rows:
        if out and r.kind == "T" and out[-1].kind == "T":
            prev = out[-1].text
            p0 = prev[:-2] if prev.endswith("**") else prev
            n0 = r.text[2:] if r.text.startswith("**") else r.text
            if p0.endswith("-") and n0[:1].islower() and not prev.rstrip().endswith("--"):
                prev_word = re.split(r"[\s(\[{,/]", p0[:-1])[-1].lower() if p0[:-1] else ""
                if re.search(r"[a-z]$", p0[:-1]) and prev_word not in KEEP_HYPHEN:
                    joined = p0[:-1] + n0
                    if prev.endswith("**") and r.text.startswith("**"):
                        pass                                  # bold run continues: marks cancel out
                    elif prev.endswith("**"):
                        joined = p0[:-1] + "**" + n0
                    elif r.text.startswith("**"):
                        joined = p0[:-1] + "**" + n0
                    out[-1] = Row(out[-1].page, out[-1].printed, joined)
                    continue
                joined = p0 + n0
                if prev.endswith("**") != r.text.startswith("**"):
                    joined = p0 + ("**" if prev.endswith("**") else "") + n0
                out[-1] = Row(out[-1].page, out[-1].printed, joined)
                continue
        out.append(r)
    return out


# ----------------------------------------------------------------------------------------------------------
# supplement instructions
# ----------------------------------------------------------------------------------------------------------
SECTION_NAMES = ["Description", "Identification", "Purity", "Assay", "Dissolution", "Uniformity of dosage units",
                 "Containers and storage", "Residue on ignition", "Loss on drying", "Content ratio of the active principle",
                 "Shelf life", "Water", "pH", "Optical rotation", "Method of preparation", "Isomer ratio",
                 "Molecular mass", "Viscosity", "Total ash", "Acid-insoluble ash", "Extract content", "Essential oil content",
                 "Congealing point", "Melting point", "Specific gravity", "Absorbance", "Osmotic pressure ratio",
                 "Disintegration", "Microbial limit", "Chemical name", "Acetic acid", "Antigenicity", "Pyrogen",
                 "Sterility", "Bacterial endotoxins", "Foreign insoluble matter", "Insoluble particulate matter",
                 "Extractable volume", "Others", "Specific activity", "Constituent amino acids", "Acid value"]
_SEC_ALT = "|".join(sorted(map(re.escape, SECTION_NAMES), key=len, reverse=True))
_SEC_ITEM_RE = re.compile(r"(" + _SEC_ALT + r")((?:\s*,?\s*(?:and\s+)?\(\d+\))*)", re.I)


def parse_instruction(text: str) -> dict:
    """Classify a supplement instruction line.

    Returns {raw, action: add|change|delete|other, kind: whole|partial|other, scope: [(section, item_no|None)],
             delete_items: [(section, item_no)], ambiguous: bool, note}
    'whole' means the text that follows (titles / entries) is added / replaced / deleted as a unit; 'partial' means
    only the named sections of the current monograph are replaced.
    """
    raw = re.sub(r"\s+", " ", text).strip()
    low = raw.lower()
    out = {"raw": raw, "action": "other", "kind": "other", "scope": [], "delete_items": [], "ambiguous": False, "note": ""}
    if re.match(r"add the following( spectra)?:?$", low):
        out.update(action="add", kind="whole")
    elif re.match(r"change the following( paragraph)? (as follows|read as follows|as follows)?:?$", low) or \
            re.match(r"change the following read as follows:?$", low) or re.match(r"change the following as follows:?$", low):
        out.update(action="change", kind="whole")
    elif re.match(r"delete the following( monographs| spectra)?:?$", low):
        out.update(action="delete", kind="whole")
    elif re.match(r"(add|delete) the following (to|from) section", low):
        out.update(action="other", kind="other", note="reference standards list")
    elif low.startswith("add the following next to") or low.startswith("add the following to"):
        out.update(action="add", kind="partial")
        m = re.search(r"next to the (.+?):?$", raw, re.I)
        out["scope"] = _scope_from(m.group(1) if m else "")[0]
    elif low.startswith("change") or low.startswith("delete") or low.startswith("remove"):
        act = "delete" if (low.startswith("delete") or low.startswith("remove")) and "change" not in low else "change"
        out.update(action=act, kind="partial")
        scope, dels = _scope_from(raw)
        out["scope"], out["delete_items"] = scope, dels
        if re.search(r"move up|sequentially", low):
            out["note"] = "renumbering: whole section replaced"
            out["scope"] = [(s, None) for s, _ in scope]
        if "arsenic item" in low or re.search(r"\bitem\b", low) and act == "delete":
            out["note"] = "item deletion without chromatographic effect"
            out["scope"] = []
        if act == "delete" and not scope:
            out["note"] = out["note"] or "deletion of unnamed text"
    else:
        out["note"] = "unrecognised instruction"
        out["ambiguous"] = True
    return out


def _scope_from(text: str):
    """Section names (+ item numbers) mentioned in an instruction -> ([(section, item_no|None)], delete_items)."""
    scope, dels = [], []
    # "delete (8)" after a section mention
    low = text
    dm = re.search(r"delete\s*\(?(\d+)\)?", low, re.I)
    for m in _SEC_ITEM_RE.finditer(text):
        name = next(s for s in SECTION_NAMES if s.lower() == m.group(1).lower())
        nums = [int(x) for x in re.findall(r"\((\d+)\)", m.group(2))]
        if nums:
            scope.extend((name, n) for n in nums)
        else:
            scope.append((name, None))
    if dm and scope and re.search(r"delete\s*\(\d+\)", text, re.I):
        sec = scope[-1][0]
        dels.append((sec, int(dm.group(1))))
    return scope, dels


# ----------------------------------------------------------------------------------------------------------
# segmentation
# ----------------------------------------------------------------------------------------------------------
def block_end(rows: list[Row], start: int) -> list[str]:
    first = rows[start].text
    k = start + 1
    if re.sub(r"[\s—*\d]", "", first) != "Operatingconditions":
        # prose mention, e.g. "Operating conditions and system suitability fundamentally follow ..."
        body = [first]
        while k < len(rows) and rows[k].kind == "T" and not (rows[k - 1].text.endswith(".") and rows[k].text[:1].isupper()):
            body.append(rows[k].text)
            k += 1
        return body
    body = [first]
    table_mode = False
    table_lines = 0
    while k < len(rows):
        if rows[k].kind != "T":
            break
        ln = rows[k].text
        prev = body[-1]
        sentence_end = prev.endswith((".", ".)", ".]", ")."))
        if STOP_RE.match(ln):
            break
        is_label = FIELD_RE.match(ln) or GENERIC_FIELD_RE.match(ln) or re.match(r"^\d+[A-Z][a-z]* [a-z]+:", ln) or LABELISH_RE.match(ln)
        if table_mode:
            if is_label and not TABLE_RE.match(ln):
                table_mode = False
            elif table_lines < 70 and not (TEST_RE.match(ln) and sentence_end and not re.match(r"^\S+\s{2}", ln)):
                body.append(ln)
                table_lines += 1
                k += 1
                continue
            else:
                break
        if re.match(r"^(?:Use a column|Proceed with \d)", ln):
            body.append(ln)
            k += 1
            continue
        if sentence_end and (TEST_RE.match(ln) or PARA_RE.match(ln)):
            break
        if sentence_end and not is_label and not TABLE_RE.match(ln) and not prev.endswith("table."):
            ahead = [r.text for r in rows[k + 1:k + 16] if r.kind == "T"]
            more = False
            for a_ in ahead:
                if OPCOND_RE.match(a_):                      # the next Operating-conditions block starts: this one is over
                    break
                if STOP_RE.match(a_) or TEST_RE.match(a_) or PARA_RE.match(a_):
                    more = a_.startswith(("System suitability", "System Suitability")); break
                if FIELD_RE.match(a_) or GENERIC_FIELD_RE.match(a_):
                    more = True; break
            if not more:
                break
        body.append(ln)
        if re.search(r"following tables?|table below|tables? shown below|shown in the (?:following )?table|in the table", ln):
            table_mode, table_lines = True, 0
        k += 1
    return body


@dataclass
class Event:
    """One revision action in a supplement (or a plain presence record in a main volume)."""
    event_id: str
    source: str
    action: str          # add / change / delete / other / main
    kind: str            # whole / partial / other
    monograph: str | None
    entry: str | None
    scope: list
    delete_items: list
    page: int
    raw: str
    ambiguous: bool = False
    note: str = ""
    part: str = ""


def lead_info(ln: str):
    """Leading bold run of a line -> (lead_text, remainder)."""
    bm = LEAD_BOLD_RE.match(ln)
    if not bm:
        return "", ln
    lead = re.sub(r"\s+", " ", bm.group(1).replace("**", "")).strip()
    lead = re.sub(r"\(\s*([A-Z])\s*\)", r"(\1)", lead)
    return lead, ln[bm.end():].strip()


def clean_title(s: str) -> str:
    s = re.sub(r"\s*<[\d.]+>.*$", "", s)
    s = re.sub(r"\s+", " ", s).strip(" —-.:")
    return s


def norm_name(s: str) -> str:
    """Tidy a reagent / entry name rebuilt from differently-fonted fragments: 'Benz[ a ]anthracene', 'n -caprylate'."""
    s = re.sub(r"\s+", " ", s).strip()
    s = re.sub(r"([\[({])\s+", r"\1", s)
    s = re.sub(r"\s+([\])},])", r"\1", s)
    s = re.sub(r"\s*′\s*", "′", s)
    s = re.sub(r"(?<![A-Za-z0-9])([A-Za-z\u0370-\u03ff]|dl|sec|tert|cis|trans)\s+-", r"\1-", s)
    s = re.sub(r"-\s+([\u0370-\u03ff])", r"-\1", s)
    s = re.sub(r"([\u0370-\u03ff])\s+-", r"\1-", s)
    s = re.sub(r"\(\s*\+\s*\)", "(+)", s)
    s = re.sub(r"\(\s*-\s*\)", "(-)", s)
    return s


def _head_sub(text: str) -> str | None:
    """'Related substances (i)—Dissolve ...' -> 'i'; '(i) Apparatus—' -> 'i'"""
    m = re.match(r"^\s*\(([ivx]{1,4})\)\s+[^—]{1,60}—", text) or re.match(r"^\s*[^—]{1,90}?\(([ivx]{1,4})\)\s*—", text)
    return m.group(1).lower() if m else None


GENERIC_TITLES = {"apparatus", "procedure", "preparation", "reagents", "reagent", "method", "equipment", "conditions", "operation", "test"}
# sub-procedure headings that only describe preparation steps; they do not name what the Operating-conditions block belongs to
GENERIC_SUB_RE = re.compile(r"^(?:apparatus|procedure|preparation(?: of .*)?|reagents?|method|equipment|conditions|operation|test|calculation(?: method)?|"
                            r"(?:standard|sample|test|control|reference|internal standard)\s+solutions?|system suitability|detection|measurement)$", re.I)
_VERB_RE = re.compile(r"\b(?:weigh|dissolve|perform|pipet|transfer|add|take|measure|proceed|use|mix|prepare|shake)\b", re.I)


def _head_title(rest_text: str) -> str | None:
    rest_text = re.sub(r"\s*\([ivx]{1,4}\)\s*(?=—)", "", rest_text, count=1)
    rest_text = re.sub(r"^\s*\(?[ivx]{1,4}\)\s+", "", rest_text, count=1)       # "(i) Apparatus—" : sub-procedure label, not a title
    mt = re.match(r"^\s*([^—]{1,170}?)—", rest_text)
    if mt and mt.group(1).strip().lower() in GENERIC_TITLES:
        return None
    if mt and not _VERB_RE.search(mt.group(1)) and not re.search(r"\.\s", mt.group(1)) \
            and not re.search(r"\b(?:the|of the|and the)\s*$", mt.group(1)):
        return clean_title(mt.group(1)) or None
    return None


def _sub_title(ln: str) -> str | None:
    """'(ii) Determination of molecular mass—Dissolve ...' -> '(ii) Determination of molecular mass' (preparation-step headings excluded)."""
    m = re.match(r"^\(([ivx]{1,4})\)\s+([^—]{1,100}?)—", ln)
    if m and not GENERIC_SUB_RE.match(m.group(2).strip()) and not re.search(r"\.\s", m.group(2)):
        return f"({m.group(1)}) {clean_title(m.group(2))}"
    return None


def segment(rows: list[Row], *, source: str, file: str, default_part: str) -> tuple[list[dict], list[Event], list[dict]]:
    """Return (blocks, events, monograph_records)."""
    blocks: list[dict] = []
    events: list[Event] = []
    monos: list[dict] = []
    part = default_part
    monograph = entry = section = None
    item_no: int | None = None
    item_title: str | None = None
    sub_label: str | None = None
    sub_disp: str | None = None        # shown in `item`: "(iii)" taken from "Related substances (iii)—", or "(ii) Determination of molecular mass"
    hint = None
    sec_start = 0
    whole = None            # active whole-scope instruction (applies to titles / entries that follow)
    partial = None          # active partial instruction (applies to the current monograph)
    cur_event: Event | None = None
    mono_rec: dict | None = None
    in_hdr = False
    i = 0
    n = len(rows)
    is_supp = source.startswith("JP18 S")

    def new_event(action, kind, mono, ent, instr, page):
        ev = Event(f"{source}#E{len(events) + 1}", source, action, kind, mono, ent,
                   instr["scope"] if instr else [], instr["delete_items"] if instr else [], page,
                   instr["raw"] if instr else "", bool(instr and instr["ambiguous"]), instr["note"] if instr else "", part)
        events.append(ev)
        return ev

    while i < n:
        r = rows[i]
        ln = r.text
        if r.kind == "I":
            ins = parse_instruction(ln)
            if ins["kind"] == "whole":
                whole, partial, cur_event = ins, None, None
            elif ins["kind"] == "partial":
                partial, whole = ins, None
                cur_event = new_event(ins["action"], "partial", monograph, entry, ins, r.page)
            else:
                whole, partial = None, None
                cur_event = new_event(ins["action"], "other", monograph, entry, ins, r.page)
            i += 1
            continue
        if r.kind == "H":
            parts = [ln]
            j = i + 1
            while j < n and rows[j].kind == "H" and ln not in PART_TITLES and rows[j].text not in PART_TITLES \
                    and not CHAPTER_RE.match(rows[j].text):
                parts.append(rows[j].text)
                j += 1
            title = norm_name(re.sub(r"(\d)\s+%", r"\1%", re.sub(r"%(?=[A-Za-z])", "% ", " ".join(parts))))
            i = j
            if title in PART_TITLES:
                part = PART_TITLES[title]
                monograph = entry = section = item_title = None
                item_no = None
                in_hdr = False
                continue
            monograph, entry, section, item_no, item_title, sub_label, sub_disp = title, None, None, None, None, None, None
            sec_start, hint = i, None
            partial = None
            partial_follows = False          # "### X / <Japanese name> / Change the Purity (3) as follows:" = partial revision of X
            for r2 in rows[i:i + 4]:
                if r2.kind == "H":
                    break
                if r2.kind == "I":
                    partial_follows = parse_instruction(r2.text)["kind"] == "partial"
                    break
            if partial_follows:
                cur_event = None
            else:
                cur_event = new_event(whole["action"] if whole else ("other" if is_supp else "main"), "whole", title, None, whole, r.page)
                if re.match(r"^9\.\d+\b", title):        # reference standards / reagents / ... : entries are added or replaced one by one
                    cur_event.kind = "container"
            mono_rec = {"source": source, "file": file, "part": part, "monograph": title, "page": r.page,
                        "event_id": cur_event.event_id if cur_event else None, "header_rows": []}
            monos.append(mono_rec)
            in_hdr = part in ("mono", "crude")
            continue

        general_chapter = part in ("gen", "geninfo") or bool(monograph and re.match(r"^\d+\.\d+\s", monograph))
        if in_hdr and mono_rec is not None:
            if re.match(r"^\*\*(Description|Identification|Purity|Assay|Method of preparation|Containers|Origin|"
                        r"Preparation|Manufacture|Shelf|Loss on|Total ash|Content)", ln) or \
                    len(mono_rec["header_rows"]) >= 40:
                in_hdr = False
            else:
                mono_rec["header_rows"].append(ln)

        tm = TEST_RE.match(ln)
        if tm and not ln.startswith("**") and not re.match(
                re.escape(tm.group(1)) + r"\s*(?:<[\d.]+>)?\s*(?:\(\d+\))?\s*[^—.]{0,80}—", ln):
            tm = None                                              # "Content (%) = ..." formula line, "Assay" mid-sentence
        run_in = None if tm else INLINE_TEST_RE.search(ln)
        lead, rest = lead_info(ln)
        if lead and "—" not in rest:                  # a long title can wrap: "(3) Total alkaloids (A and B, or C" / "and D)—Weigh ..."
            rest = " ".join([rest] + [rows[k].text for k in range(i + 1, min(i + 3, n)) if rows[k].kind == "T" and not rows[k].text.startswith("**")])

        if general_chapter:
            if tm:
                section, item_no, sec_start = tm.group(1), None, i
                m_it = re.search(r"\((\d+)\)", ln[tm.end():tm.end() + 12])
                if m_it:
                    item_no = int(m_it.group(1))
                item_title = _head_title(ln[tm.end():])
            elif run_in:
                section, item_no, sec_start = run_in.group(1), None, i
                m_it = re.search(r"\((\d+)\)", ln[run_in.end(1):run_in.end(1) + 12])
                if m_it:
                    item_no = int(m_it.group(1))
                item_title = _head_title(ln[run_in.end(1):])
            elif lead and re.match(r"\(\d+\)$", lead):
                item_no = int(re.match(r"\((\d+)\)", lead).group(1))
                item_title = _head_title(rest)
            elif lead and not re.match(r"^(Fig|Table|Figure)\b", lead):
                prev_ln = rows[i - 1] if i else None
                if prev_ln is not None and prev_ln.kind == "T" and prev_ln.text.startswith("**") and prev_ln.text.endswith("**") \
                        and entry and norm_name(_strip_marks(prev_ln.text)) == entry \
                        and not re.match(r"^\d+(\.\d+)*\.?\s", lead):
                    entry = norm_name(f"{entry} {lead}")
                else:
                    entry = norm_name(lead)
                    if whole and part == "gen":
                        cur_event = new_event(whole["action"], "whole", monograph, entry, whole, r.page)
                section, item_no, item_title, sec_start, hint = None, None, None, i, None
            elif not lead and section:
                m_it = re.match(r"^\((\d+)\)\s", ln)
                if m_it:
                    item_no = int(m_it.group(1))
                    t = _head_title(ln[m_it.end():])
                    if t:
                        item_title = t
        else:
            if lead:
                m = re.match(r"^(.*?)\s*\((\d+)\)$", lead)
                if re.match(r"^\(\d+\)$", lead):
                    item_no, item_title, sec_start = int(lead.strip("()")), _head_title(rest), i
                elif m and m.group(1):
                    section, item_no, item_title, sec_start = m.group(1).strip(), int(m.group(2)), _head_title(rest), i
                else:
                    section, item_no, item_title, sec_start = norm_name(lead), None, _head_title(rest), i
            elif tm or run_in:
                mm = tm or run_in
                section, item_no, sec_start = mm.group(1), None, i
                tail = ln[tm.end():] if tm else ln[run_in.end(1):]
                m_it = re.search(r"\((\d+)\)", tail[:12])
                if m_it:
                    item_no = int(m_it.group(1))
                item_title = _head_title(tail)
            else:
                m_it = re.match(r"^\((\d+)\)\s+([^—]{1,90}?)—", ln)
                if m_it and section:
                    item_no, item_title = int(m_it.group(1)), clean_title(m_it.group(2))

        if lead or tm or run_in:
            sub_label = sub_disp = None                        # new section / item / entry
            tail = rest if lead else (ln[tm.end():] if tm else ln[run_in.end(1):])
            sub_label = _head_sub(tail)
            if sub_label and re.match(r"^\s*\(", tail) is None:
                sub_disp = f"({sub_label})"
            elif sub_label:
                sub_disp = _sub_title(tail.strip())
        ms = SUB_RE.match(ln) if not lead else None
        if ms:
            sub_label = (ms.group(1) or ms.group(2) or ms.group(3)).lower()
            sub_disp = _sub_title(ln) if ms.group(1) or ms.group(2) else None
        if re.search(r"Gas Chromatography|by gas chromatography", ln, re.I):
            hint = "GC"
        if re.search(r"Liquid Chromatography|liquid chromatography", ln, re.I):
            hint = "LC"

        if OPCOND_RE.match(ln):
            body = block_end(rows, i)
            m_lab = re.match(r"^(?:\*\*)?Operating conditions\s*(\d+)?", ln)
            ctx_lines = [x.text for x in rows[max(sec_start, i - 120):i] if x.kind == "T"]
            ins = partial or whole
            blocks.append({
                "source": source, "file": file, "page": r.page, "page_printed": r.printed, "part": part,
                "monograph": monograph, "entry": entry, "section": section, "item_no": item_no, "item_title": item_title,
                "cond_label": m_lab.group(1) if m_lab and m_lab.group(1) else None, "mono_idx": len(monos) - 1,
                "sub_label": sub_label, "sub_disp": sub_disp,
                "technique_hint": hint, "event_id": cur_event.event_id if cur_event else None,
                "instr": {"action": ins["action"], "kind": ins["kind"], "raw": ins["raw"]} if ins else None,
                "ctx_text": "\n".join(ctx_lines)[-3000:], "text": "\n".join(body).replace("**", ""),
                "after": [x.text for x in rows[i + len(body):i + len(body) + 3]],
            })
            i += len(body)
            continue
        i += 1
    # several Operating-conditions blocks of one test that belong to different roman-numeral sub-procedures: show which one in `item`
    groups: dict[tuple, list[dict]] = {}
    for b in blocks:
        groups.setdefault((b["mono_idx"], b["section"], b["item_no"], b["entry"]), []).append(b)
    for grp in groups.values():
        subs = [b["sub_label"] for b in grp if b["sub_label"] and re.fullmatch(r"[ivx]{1,4}", b["sub_label"])]
        if len(grp) > 1 and len(set(subs)) > 1:
            for b in grp:
                if not b["sub_disp"] and b["sub_label"] and re.fullmatch(r"[ivx]{1,4}", b["sub_label"]):
                    b["sub_disp"] = f"({b['sub_label']})"
    return blocks, events, monos


def _entry_name(lead: str) -> str:
    return re.sub(r"\s+", " ", lead).strip()
