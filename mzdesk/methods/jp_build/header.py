"""Monograph header parser: the lines between a monograph title and its first bold section heading.

    Acebutolol Hydrochloride            <- title (### line)
    アセブトロール塩酸塩                   <- Japanese name
    C18H28N2O4.HCl: 372.89              <- molecular formula: molecular weight
    N-{3-Acetyl-4-[(2RS)-2-hydroxy-     <- chemical name (wrapped)
    3-(1-methylethyl)aminopropyloxy]phenyl}butanamide monohydrochloride
    [34381-68-5]                        <- CAS RN(s) in brackets
    Acebutolol Hydrochloride, when dried, contains not less than 98.0% ... of acebutolol hydrochloride (C18H28N2O4.HCl).
"""
from __future__ import annotations

import re

CJK_RE = re.compile(r"[぀-ヿ㐀-鿿]")
CAS_RE = re.compile(r"\b(\d{2,7}-\d{2}-\d)\b")
CAS_BRACKET_RE = re.compile(r"\[(\d{2,7}-\d{2}-\d(?:\s*,\s*[^\]]*)?)\]")
FORMULA_RE = re.compile(r"^(?P<f>[A-Z(][A-Za-z0-9().·,\[\]{}+\-ξ]*\d[A-Za-z0-9().·,\[\]{}+\-ξ]*|\([A-Za-z0-9]+\)n[A-Za-z0-9.·]*)"
                        r"(?::\s*(?P<mw>\d[\d,]*(?:\.\d+)?))?$")
FORMULA_ANY_RE = re.compile(r"[(\[](?P<f>[A-Z(][A-Za-z0-9().·]*?)(?::\s*(?P<mw>\d[\d,]*(?:\.\d+)?))?[)\]]")
FORMS = r"(?:Tablets?|Capsules?|Injection|Granules|Fine Granules|Syrup|Ointment|Cream|Gel|Lotion|Suppositor(?:y|ies)|Powder|" \
        r"Solution|Oral Solution|Ophthalmic Solution|Ophthalmic Ointment|Enteric-coated Tablets|Orally Disintegrating Tablets|" \
        r"Delayed-release (?:Fine )?Granules|Extended-release Tablets|Tablets for [A-Za-z ]+|for Injection|for Syrup|" \
        r"Dry Syrup|Elixir|Spirit|Plaster|Patch|Tape|Inhalation|Eye Drops|Nasal Drops|Aerosol|Liquid|Suspension|Spray|Cataplasm|" \
        r"Combination Granules|Chewable Tablets|Liniment|Injectable Aqueous Suspension)(?: [A-Z])?"
STATEMENT_RES = [
    re.compile(r"(?:labeled|labelled)\s+(?:amounts?|potency)\s+of\s+(?:the\s+)?(?P<n>[^()\[]+?)\s*[(\[]", re.I),
    re.compile(r"\(potency\)\s+of\s+(?:the\s+)?(?P<n>[^()\[]+?)\s*[(\[]", re.I),
    re.compile(r"(?:contains?|containing)\s+not\s+less\s+than\s+[\d.,]+\s*(?:%|w/v%|vol%|w/w%)[^()]*?\bof\s+(?:the\s+)?"
               r"(?P<n>[A-Za-z0-9][^()\[]*?)\s*[(\[]\s*(?=[A-Z(])", re.I),
]
NEXT_SUBST_RE = re.compile(r"^\s*[)\]]?\s*(?:,\s*)?(?:and\s+)(?P<n>[A-Za-z][^()\[]*?)\s*[(\[]\s*(?P<f>[A-Z(][^)\]]*)[)\]]")


def _clean_name(n: str) -> str:
    n = re.sub(r"\s+", " ", n).strip(" ,;:")
    n = re.sub(r"^(?:the\s+)?(?:labeled|labelled)\s+(?:amount|potency)\s+of\s+", "", n, flags=re.I)
    n = re.sub(r"^the\s+", "", n, flags=re.I)
    n = re.sub(r"^an amount of\s+", "", n, flags=re.I)
    n = re.split(r"\s+(?:corresponding|expressed|calculated)\b", n)[0]
    return n.strip(" ,;:")


def _paren_end(s: str) -> int:
    """index just after the bracket group that starts at s[0] ('(' or '[')"""
    if not s or s[0] not in "([":
        return 0
    depth = 0
    for i, ch in enumerate(s):
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
            if depth == 0:
                return i + 1
    return 0


def bulk_name(title: str) -> str:
    """Substance named by a bulk monograph title: 'Isosorbide Mononitrate 70%/Lactose 30%' -> 'Isosorbide Mononitrate'."""
    t = re.sub(r"\s+\d+(?:\.\d+)?\s*%\s*/.*$", "", title)
    t = re.sub(r"^(?:Biphasic\s+)?Isophane\s+", "", t)
    return t.strip()


def clean_join(lines: list[str]) -> str:
    out = ""
    for ln in lines:
        ln = ln.replace("**", "").strip()
        if not out:
            out = ln
        elif out.endswith("-"):
            out += ln                                   # chemical names break after real hyphens only
        else:
            out += " " + ln
    return re.sub(r"\s+", " ", out).strip()


def _chem_join(lines: list[str]) -> str:
    return clean_join(lines)


def parse_header(rows: list[str], title: str) -> dict:
    out = {"jp_name": None, "alt_names": [], "formula": None, "mw": None, "chemical_name": None, "cas": [],
           "cas_note": None, "substance": None, "statement": None}
    if not rows:
        return out
    lines = [r.replace("**", "").strip() for r in rows if r.strip()]
    # Japanese name = first line with kana/kanji; whatever precedes it (Latin / alternative names) are alt names
    j = next((k for k, l in enumerate(lines) if CJK_RE.search(l)), None)
    pos = 0
    if j is not None:
        out["jp_name"] = re.sub(r"\s+", "", lines[j])
        out["alt_names"] = [l for l in lines[:j] if not re.match(r"^\d", l)]
        pos = j + 1
    # formula line
    k = pos
    chem_start = None
    while k < len(lines):
        l = lines[k]
        m = FORMULA_RE.match(l.replace(" ", ""))
        if m and not re.search(r"[a-z]{4,}", l) and (m.group("mw") or re.search(r"\d", l)) and len(l) < 60:
            out["formula"] = m.group("f")
            out["mw"] = float(m.group("mw").replace(",", "")) if m.group("mw") else None
            chem_start = k + 1
            break
        if k - pos >= 3:
            break
        k += 1
    start = chem_start if chem_start is not None else pos
    # chemical name: lines up to and including the line that carries the CAS bracket
    cas_line = next((q for q in range(start, len(lines)) if CAS_BRACKET_RE.search(lines[q])), None)
    if out["formula"] is not None and cas_line is not None:
        name_lines = lines[start:cas_line + 1]
        name = _chem_join(name_lines)
        brackets = CAS_BRACKET_RE.findall(name)
        name = CAS_BRACKET_RE.sub("", name).strip()
        out["chemical_name"] = name or None
        for b in brackets:
            out["cas"] += [c for c in CAS_RE.findall(b) if c not in out["cas"]]
        notes = [re.sub(CAS_RE, "", b).strip(" ,") for b in brackets]
        out["cas_note"] = "; ".join(n for n in notes if n) or None
    elif cas_line is not None:
        # CAS without formula (e.g. peptide drugs): take it, keep the name
        for b in CAS_BRACKET_RE.findall(lines[cas_line]):
            out["cas"] += [c for c in CAS_RE.findall(b) if c not in out["cas"]]
    # statement paragraph (first sentence stating the content)
    body = clean_join(lines[(cas_line + 1) if cas_line is not None else (chem_start if chem_start is not None else pos):])
    out["statement"] = body[:600] or None
    out["substances"] = []
    for rx in STATEMENT_RES:
        for m in rx.finditer(body):
            name = _clean_name(m.group("n"))
            if not name or len(name) > 90 or re.search(r"\b(is|are|per|expressed|labeled|potency)\b", name, re.I):
                continue
            out["substance"] = name
            bp = m.end("n")
            while bp < len(body) and body[bp] not in "([":
                bp += 1
            fm = FORMULA_ANY_RE.match(body[bp:bp + 90])
            subs = [{"name": name, "formula": fm.group("f") if fm else None,
                     "mw": float(fm.group("mw").replace(",", "")) if fm and fm.group("mw") else None}]
            if fm and out["formula"] is None:
                out["formula_stmt"] = fm.group("f")
                if fm.group("mw"):
                    out["mw_stmt"] = float(fm.group("mw").replace(",", ""))
            # combination preparations: "... of the labeled amount of A (C..: ..) and B (C..: ..)"
            rest = body[bp:]
            close = _paren_end(rest)
            while close:
                rest = rest[close:]
                nm = NEXT_SUBST_RE.match(rest)
                if not nm:
                    break
                n2 = _clean_name(nm.group("n"))
                if not n2 or len(n2) > 60 or re.search(r"\b(is|are|per|expressed|calculated|labeled)\b|%|not (?:less|more)", n2, re.I):
                    break
                f2 = re.match(r"([A-Z(][A-Za-z0-9().·]*)(?::\s*([\d,.]+))?", nm.group("f"))
                subs.append({"name": n2, "formula": f2.group(1) if f2 else None,
                             "mw": float(f2.group(2).replace(",", "")) if f2 and f2.group(2) else None})
                rest = rest[nm.end() - 1:]
                close = _paren_end(rest)
            for m2 in STATEMENT_RES[0].finditer(body):             # "... of the labeled amount of A (..) and ... of the labeled amount of B (..)"
                n3 = _clean_name(m2.group("n"))
                if n3 and len(n3) < 70 and n3.lower() not in {x["name"].lower() for x in subs} and not re.search(r"\b(is|are|per|expressed)\b|%", n3, re.I):
                    bp2 = m2.end("n")
                    while bp2 < len(body) and body[bp2] not in "([":
                        bp2 += 1
                    f3 = FORMULA_ANY_RE.match(body[bp2:bp2 + 90])
                    subs.append({"name": n3, "formula": f3.group("f") if f3 else None,
                                 "mw": float(f3.group("mw").replace(",", "")) if f3 and f3.group("mw") else None})
            out["substances"] = subs
            break
        if out["substances"]:
            break
    if out["substance"] is None:
        base = re.sub(r"\s+" + FORMS + r"$", "", title).strip()
        if base != title:
            out["substance"] = base
            out["substance_from_title"] = True
    return out
