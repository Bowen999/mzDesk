"""What does a method measure? -> `analytes` (lower-case names).

* drug monographs: the monograph substance ("related substances" is added for Purity / related-substance tests)
* crude drugs / Kampo extracts: marker compound(s) named in the item title ("Assay (1) Paeoniflorin"), else the
  compound named in the test's quantity formula ("Amount (mg) of barbaloin = ...")
* reagent entries: the reagent compound ("Acemetacin for assay" -> acemetacin)
* a few general chapters: fixed by hand (ethanol, vinyl chloride, ...)
"""
from __future__ import annotations

import re

REAGENT_SUFFIX = re.compile(
    r"\s+(?:for|to be used for)\s+(?:the\s+)?(?:assay|purity|resolution check|thin-layer chromatography|liquid chromatography|"
    r"gas chromatography|identification|quantitative analysis|system suitability|sulbactam penicillamine|"
    r"crude drugs?\s+purity\s+test|component determination|component analysis|content determination|"
    r"liquid chromatography for crude drugs|nuclear magnetic resonance|qNMR)\b.*$", re.I)
REAGENT_TS = re.compile(r"\b(?:TS|VS|RS)\b")
FORM_WORDS = re.compile(r"\s+(?:for injection|tablets?|capsules?|injection|granules|fine granules|syrup|for syrup|ointment|"
                        r"ophthalmic ointment|powder|oral solution|solution|suppositor(?:y|ies)|cream|gel|lotion)$", re.I)
GENERIC = {"related substances": "related substances", "related substance": "related substances"}
MANUAL = {  # (monograph, entry-prefix) -> analytes
    ("Test for Trace Amounts of Aluminum in Total Parenteral Nutrition (TPN) Solutions", None): ["aluminum"],
    ("1.01 Alcohol Number Determination", None): ["ethanol"],
    ("5.01 Crude Drugs Test", "4.3."): ["bhc isomers", "ddt isomers"],
    ("2.46 Residual Solvents", None): ["residual solvents"],
}
GENERAL_ITEM_STOP = re.compile(r"^(?:related substances?|purity|assay|content|identification|total|components?)$", re.I)


def split_markers(title: str) -> list[str]:
    """'Ephedrine and pseudoephedrine' / 'Baicalin, glycyrrhizic acid' -> ['ephedrine', 'pseudoephedrine'];
    'Total alkaloids (ephedrine and pseudoephedrine)' -> ['ephedrine', 'pseudoephedrine']; 'Berberine [as berberine chloride]' -> ['berberine']."""
    t = re.sub(r"\s+", " ", title).strip()
    m = re.match(r"^(?:total\s+[a-z ]+?|[a-z ]*(?:alkaloid|saponin|flavonoid|curcuminoid|capsaicin|ginsenoside|glycoside|sennoside)s?)\s*[(\[]([^)\]]+)[)\]]$", t, re.I)
    if m and not re.match(r"^(?:as|calculated)\b", m.group(1), re.I):
        t = m.group(1)
    t = re.sub(r"\s*[(\[](?:as|calculated)\b[^)\]]*[)\]]", "", t, flags=re.I)
    parts = re.split(r",\s+(?:and\s+)?|\s+and\s+|\s*&\s*|\s+\+\s+", t)
    out = [re.sub(r"^or\s+", "", p.strip(" ,;")) for p in parts if p and p.strip(" ,;")]
    return list(dict.fromkeys(p.lower() for p in out))


_STOP = re.compile(r"\b(?:is|are|not|by|the|equation|equations|obtained|designate|following|sum|total as|each|using|taken|solution|"
                   r"sample|standard|internal|calculated|dried|basis|anhydrous)\b", re.I)


def _split_names(txt: str) -> list[str]:
    parts = re.split(r",\s+(?:and\s+)?|\s+and\s+|;\s*", txt)
    out = []
    for p in parts:
        p = re.sub(r"\s+", " ", p).strip(" .;:,")
        p = re.sub(r"(?<=\S)\s+[\[(].*$", "", p)              # "berberine [as berberine chloride" -> berberine
        p = re.sub(r"\s+for\s+.*$", "", p, flags=re.I)          # "emetine hydrochloride for assay of crude drugs"
        p = re.sub(r"\s+and$", "", p.strip(" .;:,")).strip(" .;:,")
        if p and len(p) < 50 and not _STOP.search(p):
            out.append(p.lower())
    return out


def _formula_names(tail: str) -> list[str]:
    """'Amount (mg) of barbaloin = ...', 'The amount (mg) of hyoscyamine (C17H23NO3) = ...',
    'Amount (mg) of total curcuminoids (curcumin, demethoxycurcumin and bisdemethoxycurcumin)'."""
    names: list[str] = []
    for m in re.finditer(r"(?m)^(?:[Tt]he )?[Aa]mounts?\s*(?:\[[^\]]*\]|\([^)]*\))?\s*of\s+(?P<n>[A-Za-z0-9\[\]\-, αβγ′()]+?)\s*"
                         r"(?P<p>\((?:C\d|[A-Z][a-z]?\d)[^)]*\)|\([^)]{3,80}\))?\s*(?:=|＝|\n|$)", tail):
        n = m.group("n").strip()
        par = m.group("p")
        if re.match(r"total\b", n, re.I) and par and not re.match(r"\((?:C\d|[A-Z][a-z]?\d)", par):
            names += _split_names(par.strip("()"))
            if not names:
                names.append(n.lower())
        else:
            names += _split_names(n)
    return list(dict.fromkeys(x for x in names if x))


def _peak_names(tail: str) -> list[str]:
    m = re.search(r"peaks? corresponding to (?P<n>[A-Za-z0-9\[\]\-, αβγ′]+?)(?: in the| from| obtained| shows?| of the)", tail)
    if m:
        return _split_names(m.group(1))
    m = re.search(r"peak areas?(?:,\s*[A-Z][A-Za-z]*(?:,? and [A-Z][A-Za-z]*)?,)? of (?P<n>[A-Za-z0-9\[\]\-, αβγ′]+?) in the", tail)
    if m:
        return _split_names(m.group(1))
    m = re.findall(r"dissolve\s+[\d.]+\s*(?:mg|g)\s+(?:each\s+)?of\s+(?P<n>[A-Za-z0-9\[\]\-, αβγ′]+?)\s+(?:for [a-z\- ]+ )?(?:in|and)\b", tail)
    if m:
        return _split_names(m[-1])
    return []


def _first_tok(n: str) -> str:
    w = re.findall(r"[a-z]{4,}", n.lower())
    return w[0] if w else n.lower()


def _pick_combo(rec: dict, subs: list[str]) -> list[str]:
    """Combination preparation: the active substance(s) this test is about (item title > what the test text mentions > flow target)."""
    title = (rec.get("item_title") or "").lower()
    by_title = [n for n in subs if _first_tok(n) in title]
    if by_title:
        return by_title
    if title and not GENERAL_ITEM_STOP.match(title) and not re.search(r"related|purity|assay|dissolution|uniformity", title):
        return [title]                                   # e.g. "Assay (2) Calcium pantothenate" when the content statement lacks it
    tail = rec.get("ctx_text", "")[-1200:].lower()
    mentioned = [n for n in subs if _first_tok(n) in tail]          # substances the test text talks about (e.g. "retention time of cefoperazone")
    if mentioned:
        return mentioned
    ft = (rec.get("flow_target_compound") or "").lower()
    by_ft = [n for n in subs if _first_tok(n) in ft]
    if by_ft:
        return by_ft
    return subs


# words that mark a title as a class / property / material rather than a compound name
_NOT_COMPOUND = re.compile(
    r"related|substance|impurit|contamin|multimer|dimer|polymer|isomer|enantiomer|diastereomer|degradation|assay|content|total|residual|"
    r"solvent|heavy|ratio|composition|purity|unity|clarity|color|colour|acidity|alkalin|molecular|peak|identification|protein|dissolved|"
    r"desamido|free |amount|specific|limit|organic|volatile|component|derivative|pyrogen|endotoxin|particle|foreign|arsenic|"
    r"\bash\b|loss|residue|viscosity|^ph$|conductivity|sulfated|halogen|oxidizing|reducing|peroxide|other|\bany\b|\beach\b|"
    r"\bforms?\b|constituent|amino acids?|isoelectric|\bpoint\b|distribution|activity|potency|\boil\b|\bfat\b|\bwax\b|\bextract\b|\bresin\b|"
    r"\bpowder\b|^the |preparation|solution|sample|mixture|^[ivx]+$|capacity|\bgroups?\b", re.I)
STANDARD_SECTIONS = {"identification", "purity", "assay", "dissolution", "uniformity of dosage units", "uniformity of dosage unit", "description",
                     "content ratio of the active principle", "isomer ratio", "unity of peak", "content", "water", "loss on drying"}
_REFERENCE_MENTION = (r"relative retention time(?:[^.]|\.\d){{0,40}}\bto\s+{tok}", r"peak areas? of {tok}(?:[^.]|\.\d){{0,30}}from the standard")


def _looks_like_compound(title: str) -> bool:
    t = title.strip()
    return bool(t) and len(t) <= 90 and len(t.split()) <= 8 and not _NOT_COMPOUND.search(t)


def title_compounds(title: str) -> tuple[list[str], bool]:
    """(compound names, mentions related substances) of an item title: 'Ethylene glycol, diethylene glycol and related substances'."""
    title = re.sub(r"^(?:relative\s+)?(?:content|amount|ratio)s?\s+of\s+", "", re.sub(r"\s*\((?:[ivx]{1,4}|\d)\)\s*$", "", title), flags=re.I)
    pieces = [re.sub(r"^(?:or|free)\s+", "", p).strip() for p in split_markers(title)]
    rs = any(re.search(r"related substance|impurit|degradation", p) for p in pieces)
    comps = [p for p in dict.fromkeys(pieces) if _looks_like_compound(p)]
    return comps, rs


_MEASURED_PATTERNS = (
    r"calculate the amounts? of (?P<n>[^:.]+?)(?: by the following|, and calculate|[:.])",
    r"peak heights?,?(?: [A-Z][A-Za-z0-9]+ and [A-Z][A-Za-z0-9]+,)? of (?P<n>[^.:;]+?) (?:of|in|from|obtained)",
    r"determine the peak areas?,?(?: [A-Z][A-Za-z0-9]+ and [A-Z][A-Za-z0-9]+,)? of (?P<n>[^.:;]+?) in (?:each|the|both)",
    r"calculate the ratios[^.]{0,80}peak areas? of (?P<n>[^.:;]+?) (?:and )?to the peak area of the internal standard",
    r"peak areas? of (?P<n>[^.:;]+?) obtained from the (?:sample|test) solution",
)


def measured_names(tail: str) -> list[str]:
    """Compounds a test quantifies, from its own text: 'Amount (%) of X' lines, 'calculate the amounts of X, Y and Z', 'peak heights of X' ..."""
    names = [n for n in _formula_names(tail) if _looks_like_compound(n)]
    if names:
        return names
    for pat in _MEASURED_PATTERNS:
        m = re.search(pat, tail)
        if m:
            got = [n for n in _split_names(m.group("n")) if _looks_like_compound(n)]
            if got:
                return got
    return []


def _reference_mention(tail: str, substance: str) -> bool:
    """True when the test compares the impurity with the substance's own peak ('relative retention time of about 0.8 to paroxetine',
    'the peak area of paroxetine from the standard solution')."""
    tok = _first_tok(substance)
    return any(re.search(p.format(tok=re.escape(tok)), tail, re.I) for p in _REFERENCE_MENTION)


def _base_substances(rec: dict, header: dict | None) -> tuple[list[str], list[str]]:
    """The substance(s) a drug monograph is about, as seen by this test: (base, all substances of a combination preparation)."""
    mono = rec["monograph"]
    subst = (header or {}).get("substance")
    subs = [x.lower() for x in (header or {}).get("substances") or []]
    ps = (header or {}).get("prep_substances")
    if ps:                                          # preparation of bulk monograph(s): "Cefteram Pivoxil Fine Granules" -> cefteram pivoxil
        subs, subst = [n for n, _ in ps], None
        if len(subs) == 1:
            subs, subst = [], subs[0]
    elif (header or {}).get("monograph_is_bulk"):
        subs, subst = [], ((header or {}).get("bulk_name") or mono)   # bulk drug monograph: the monograph substance is the title
    base: list[str] = []
    if len(subs) > 1 and not (header or {}).get("formula_source") == "header":
        base = _pick_combo(rec, subs)               # combination preparation: the active substance(s) this test is about
    elif subst:
        base = [subst.lower()]
    return base, subs


def substance_for(rec: dict, header: dict | None) -> list[str]:
    """Monograph substance(s) of the test's sample (not necessarily what is measured: a Purity test of Captopril measures its disulfide)."""
    return _base_substances(rec, header)[0] if rec["part"] == "mono" else []


def analytes_for(rec: dict, header: dict | None) -> tuple[list[str], str]:
    """Return (analytes, source-note). The analytes are what the test MEASURES."""
    mono, entry, sec = rec["monograph"], rec.get("entry"), (rec.get("section") or "")
    part = rec["part"]
    title = rec.get("item_title") or ""
    rs_title = bool(re.search(r"related substance", f"{title} {sec}", re.I))
    # --- general chapters / reagents
    if part in ("gen", "geninfo"):
        for (m, e), val in MANUAL.items():
            if mono == m and (e is None or (entry or "").startswith(e)):
                return list(val), "manual"
        if re.match(r"^9\.41", mono) and entry:
            e = REAGENT_SUFFIX.sub("", entry).strip()
            if REAGENT_TS.search(entry) or re.search(r"\bTS\b", e):
                return [e.lower()], "reagent solution"
            return [e.lower()] + (["related substances"] if rs_title else []), "reagent name"
        if title and not GENERAL_ITEM_STOP.match(title):
            return split_markers(title), "item title"
        if entry:
            return [re.sub(r"^[\d.]+\s*", "", entry).lower()], "entry"
        return [], "none"
    tail = rec.get("ctx_text", "")[-1800:]
    # --- crude drugs / extracts: marker compound(s)
    if part == "crude":
        if title and not GENERAL_ITEM_STOP.match(title) and not re.match(r"^(?:foreign|heavy|total ash|loss|extract|residual)", title, re.I) \
                and not re.search(r"\b(?:oil|fat|wax|resin)\b", title, re.I):
            return split_markers(title), "item title"
        names = _formula_names(tail)
        if names:
            return names, "quantity formula"
        names = _peak_names(tail) or measured_names(tail)
        if names:
            return names, "peak statement"
        if title:                                    # material named in the title (e.g. an adulterant) and nothing better
            return split_markers(title), "item title"
        return [], "none"
    # --- drug monographs
    out: list[str] = []
    base, subs = _base_substances(rec, header)
    # compounds named by the item title (or by a non-standard section heading such as "Acetone", "Ethanol") win over the substance
    heading = title or (sec if sec.lower() not in STANDARD_SECTIONS else "")
    comps, rs = title_compounds(heading) if heading else ([], False)
    how = "monograph substance"
    if comps:
        if subs:                                    # combination preparation: keep the preparation's own name for the active
            comps = [next((n for n in subs if _first_tok(n) == _first_tok(c)), c) for c in comps]
        out = list(comps)
        how = "item title"
        if sec == "Purity" and base and not any(_first_tok(b) == _first_tok(c) for b in base for c in comps) \
                and _reference_mention(tail, base[0]):
            out += [b for b in base if b not in out]          # impurity measured against the substance's own peak
    elif heading and not rs_title and not rs and re.search(r"residual solvent|volatile|organic", heading, re.I):
        got = measured_names(tail)                  # "Residual solvents": the solvents the test quantifies
        out = got or base
        how = "test text" if got else "monograph substance"
    else:
        out = list(base)
    if not out:
        out = [mono.lower()]
    if (rs_title or rs or (sec == "Purity" and re.search(r"related|impurit|degradation", title, re.I))) and "related substances" not in out:
        out.append("related substances")
    return out, how
