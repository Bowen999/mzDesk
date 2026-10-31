"""Candidates of the EI MS Annotation workflow (steps 4-5): the original compound behind every library hit and every user line.

A library hit is the compound the spectrum belongs to. MassBank stores a GC derivative in one of two ways:
  * as the DERIVATIVE structure ("LACTIC ACID-DI-TMS", SMILES with O-Si(CH3)3) - the way back to the original is derivatization.underivatize_info();
  * as the ORIGINAL structure with the derivative's spectrum (the GC-TOF metabolite entries: "Glycine", SMILES NCC(=O)O, base peak m/z 174).
Both end as one candidate per original compound, merged by the first block of the original's InChIKey (annotation.merge.structure_keys).

Parent names (parent_name): the name of a derivative entry minus its derivative words. The words come in the forms found in the local
MassBank EI index and in NIST-style names; PARENT_TAILS lists them (a regular expression each, applied to the END of the name until nothing
changes):  "-DI-TMS", " MONO-TMS", " 3TMS", " (2TMS)", "-3-TMS", " TMS ETHER", " TMS ESTER", ", TMS derivative"  |  " TRIMETHYLSILYL ESTER",
" BIS(TRIMETHYLSILYL) ESTER", ", N,O-bis(trimethylsilyl)-"  |  " MEOX", " methoxime". A name that still holds silyl / siloxy / TMS / TBDMS / MeOx after
that (a prefix form "TRIMETHYLSILYL HEXANOATE", an inner "(TRIMETHYLSILYL)OXY") is not trusted: the name of a non-derivative entry with the same
original structure is taken from the library instead, and only if there is none the stripped name stays.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field, replace

# ---------------------------------------------------------------- parent names
_COUNT = r"(?:\d+|mono|di|bis|tri|tris|tetra|tetrakis|penta|pentakis|hexa|hexakis|hepta|octa)"
_SILYL = r"(?:tert-butyldimethylsilyl|trimethylsilyl)"
PARENT_TAILS = [
    re.compile(rf"[\s,;:-]*\(?\s*(?:{_COUNT}[\s-]*)?(?:TMS|TBDMS|TBS|TMSi)\s*\)?(?:[\s-]*(?:derivative|ether|ester))?\s*$", re.I),
    re.compile(rf"[\s,;:-]*(?:[\d,]+-)?(?:{_COUNT}[\s-]*)?(?:[NOS](?:,[NOS])*-)?(?:{_COUNT}[\s-]*)?\(?{_SILYL}\)?(?:[\s-]*(?:derivative|ether|ester))?\s*-?\s*$", re.I),
    re.compile(r"[\s,;:-]*(?:\d+\s*)?(?:MEOX|MOX|methoxime|methyloxime|methoxyamine)(?:\s+derivative)?\s*$", re.I),
]
# inside the name: "2-(TRIMETHYLSILYL)OXY..." / "(TERT-BUTYLDIMETHYLSILYLOXY)" / "BIS(TRIMETHYLSILYL)OXY" / "...SILOXY..." is a hydroxy group (a multiplier before it
# makes dihydroxy / trihydroxy), "1,3-DI-O-TRIMETHYLSILYL " / "N-TRIMETHYLSILYL" a substituent
_SILOXY = r"(?:trimethyl|tert-butyldimethyl)siloxy"
_HYDROXY = re.compile(rf"(?:\b(?P<mult>bis|tris|di|tri)[\s-]*)?(?:\(\s*{_SILYL}\s*OXY\s*\)|\(\s*{_SILOXY}\s*\)|\(\s*{_SILYL}\s*\)\s*OXY|{_SILYL}\s*OXY|{_SILOXY})", re.I)
_HYDROXY_WORD = {None: "hydroxy", "bis": "dihydroxy", "di": "dihydroxy", "tris": "trihydroxy", "tri": "trihydroxy"}
_SUBSTITUENT = re.compile(rf"(?:[\w,]+-)?(?:(?:di|tri|tetra|penta|bis|tris)-?)?[NOS]-{_SILYL}\s*", re.I)
_PREFIX_ESTER = re.compile(rf"^{_SILYL}\s+(.+?)(oate|ate)$", re.I)
_LEFTOVER = re.compile(r"silyl|siloxy|silane|\bTMS|\bTBDMS|\bTBS\b|\bMEOX|\bMOX\b|methoxime|methyloxime", re.I)
_STEREO_PREFIX = re.compile(r"\b(dl|d|l|n|o|s|r|e|z)(?=-)", re.I)
_STEREO_GROUP = re.compile(r"\(\s*\d*[rsez]\*?(?:\s*,\s*\d*[rsez]\*?)*\s*\)", re.I)          # "(R)", "(2S,3S)", "(4R*,5R*)", "(E)"


def tidy_name(name):
    """MassBank writes many names in capitals ("TARTARIC ACID"): lower-case them with a leading capital and keep d/l/n/o/s/r/e/z prefixes in capitals.
    A single capital word of up to three letters (DHP, TNT) is an acronym and stays."""
    text = (name or "").strip()
    if text and text.upper() == text and not (len(text) <= 3 and text.isalpha()):
        text = _STEREO_GROUP.sub(lambda m: m.group(0).upper(), _STEREO_PREFIX.sub(lambda m: m.group(1).upper(), text.capitalize()))
    return text


def parent_name(name):
    """The name of a derivative entry without its derivative words -> (name, clean). clean is False when derivative words are left (see the module
    docstring): the caller should look the compound up instead. An ordinary name comes back unchanged (clean True)."""
    text = (name or "").strip()
    ester = _PREFIX_ESTER.match(text)
    if ester:                                                   # "TRIMETHYLSILYL HEXANOATE" -> "HEXANOIC ACID"
        text = ester.group(1) + ("OIC ACID" if ester.group(2).lower() == "oate" else "IC ACID")
        text = text if text.isupper() or not name.isupper() else text.upper()
    for _ in range(4):
        before = text
        for tail in PARENT_TAILS:
            text = tail.sub("", text)
        upper = text.isupper()
        text = _HYDROXY.sub(lambda m: _HYDROXY_WORD[(m.group("mult") or "").lower() or None].upper() if upper else _HYDROXY_WORD[(m.group("mult") or "").lower() or None], text)
        text = _SUBSTITUENT.sub("", text)
        if text == before:
            break
    text = text.strip(" ,;:-")
    return text, len(text) >= 3 and not _LEFTOVER.search(text)


def library_names(ik14s, libraries=("MassBank",)):
    """{first InChIKey block: name} for compounds the local libraries know under a plain (non-derivative) record; {} when no library is installed.
    The most frequent tidy name wins, the shorter one on a tie."""
    wanted = {k for k in ik14s if k}
    if not wanted:
        return {}
    import sqlite3

    from mzdesk.ms2 import library_index as lx

    counts = {}
    for library in libraries:
        try:
            lib = lx.load(library)
            con = sqlite3.connect(f"{(lib.path / 'meta.sqlite').resolve().as_uri()}?mode=ro", uri=True)
        except Exception:                                   # noqa: BLE001 - names are a convenience
            continue
        try:
            marks = ",".join("?" * len(wanted))
            for key, name in con.execute(f"SELECT substr(inchikey, 1, 14), name FROM meta WHERE substr(inchikey, 1, 14) IN ({marks}) "
                                         "AND smiles NOT LIKE '%Si%'", sorted(wanted)):
                if name and not _LEFTOVER.search(name):
                    counts.setdefault(key, {}).setdefault(tidy_name(name), 0)
                    counts[key][tidy_name(name)] += 1
        finally:
            con.close()
    return {key: min(names, key=lambda n: (-names[n], len(n), n)) for key, names in counts.items()}


# ---------------------------------------------------------------- the candidate
@dataclass
class Candidate:
    """One original compound. name: the parent name; key: first block of the original's InChIKey (else the lower-case name); parent_smiles: the
    original structure, None for a library hit that has none. sources: "library" and / or "list". library: the best EI library spectrum
    {library, accession, name, score, matches, metric, smiles (as stored), entry_mode (derivative the stored structure is, or None), state
    (derivatization the spectrum itself shows: TMS | TBDMS | None), row} or None; n_spectra: how many hits of this compound. variants: the
    derivative structures of step 5, chosen: the one that stays (step 7). The rest is filled by the ranking."""
    name: str
    key: str
    parent_smiles: str | None = None
    sources: list = field(default_factory=list)
    user_text: str = ""
    library: dict | None = None
    n_spectra: int = 0
    variants: list = field(default_factory=list)
    chosen: dict | None = None
    tier: str = "neims"
    lib_score: float | None = None
    neims_score: float | None = None
    m_status: str | None = None
    m_partial: int | None = None
    state_ok: bool = True
    flags: list = field(default_factory=list)
    level: str | None = None
    level_reason: str = ""
    rank: int = 0
    name_note: str = ""


def largest_fragment(smiles):
    """The largest component of a SMILES (a salt or a mixture), as given when it has one component; None when RDKit cannot read it."""
    from rdkit import Chem
    from rdkit.rdBase import BlockLogs

    with BlockLogs():
        mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None or not mol.GetNumAtoms():
        return None
    parts = Chem.GetMolFrags(mol, asMols=True)
    if len(parts) == 1:
        return smiles
    return Chem.MolToSmiles(max(parts, key=lambda p: p.GetNumHeavyAtoms()))


def _original(smiles):
    """-> (original SMILES, derivative mode of the given structure or None, n_changes)."""
    from mzdesk.platform_fit import derivatization as d

    info = d.underivatize_info(smiles)
    return (info["parent"], info["mode"], info["n_tms"] + info["n_tbdms"] + info["n_meox"]) if info["changed"] else (smiles, None, 0)


def from_hits(hits, primary, *, limit=10):
    """Library hits (library_search.search()["hits"], best first) -> up to ``limit`` Candidates, one per original compound, best spectrum kept.
    primary: the metric the hits are ranked by. The hit score is hit["scores"][primary]."""
    from mzdesk.annotation.merge import structure_keys

    found = {}
    for hit in hits:
        score = (hit.get("scores") or {}).get(primary)
        if score is None:
            continue
        smiles = largest_fragment(hit.get("smiles") or "")
        entry_mode, original = None, smiles
        if smiles and "Si" in smiles:
            original, entry_mode, _ = _original(smiles)
        ik14 = structure_keys(original)[1] if original else None
        if not ik14:
            ik14 = (hit.get("inchikey") or "")[:14] or None
        name, clean = hit.get("name") or "", True
        if entry_mode or not smiles:                            # a derivative entry, or one without a structure whose name may still carry the derivative words
            stripped, clean = parent_name(name)
            if entry_mode or clean:
                name = stripped
        key = ik14 or name.strip().lower() or f"{hit['library']}:{hit['accession']}"
        record = {"library": hit["library"], "accession": hit["accession"], "name": hit.get("name") or "", "score": float(score), "matches": int(hit.get("matches") or 0),
                  "metric": primary, "smiles": smiles, "entry_mode": entry_mode, "state": None, "row": hit.get("row"), "scores": dict(hit.get("scores") or {})}
        cand = found.get(key)
        if cand is None:
            if len(found) >= limit:
                continue
            found[key] = Candidate(name=tidy_name(name) or original or hit["accession"], key=key, parent_smiles=original, sources=["library"], library=record, n_spectra=1,
                                   name_note="" if clean else "derivative name")
        else:
            cand.n_spectra += 1
            if (record["score"], record["matches"]) > (cand.library["score"], cand.library["matches"]):
                cand.library = record
    return list(found.values())


def from_lines(lines, to_smiles):
    """The candidate list of the form -> (Candidates, warnings). A line that RDKit reads is a SMILES; any other line is a name looked up with
    ``to_smiles`` (online); a line that cannot be resolved is skipped with a warning. A structure that is itself a GC derivative becomes its
    original compound."""
    from rdkit import Chem
    from rdkit.rdBase import BlockLogs

    from mzdesk.annotation.merge import structure_keys

    out, warnings = [], []
    for line in lines:
        with BlockLogs():
            is_smiles = not any(ch.isspace() for ch in line) and Chem.MolFromSmiles(line) is not None       # RDKit reads "PC 33:1" as the SMILES "PC" plus a name
        name = ""
        if is_smiles:
            smiles = line
        else:
            try:
                smiles = to_smiles(line)
            except Exception as exc:                        # noqa: BLE001 - a lookup that fails skips this line only
                warnings.append(f"Candidate “{line}” was skipped: {exc}")
                continue
            name = line
        smiles = largest_fragment(smiles)
        if smiles is None:
            warnings.append(f"Candidate “{line}” was skipped: no usable structure.")
            continue
        original, mode, _ = _original(smiles)
        if mode:
            warnings.append(f"Candidate “{line}” is a {mode} derivative; its original compound is used.")
        ik14 = structure_keys(original)[1]
        out.append(Candidate(name=name or original, key=ik14 or original, parent_smiles=original, sources=["list"], user_text=line))
    return out, warnings


def merge(*groups):
    """Candidates of several sources -> one per key, in order of first appearance (library hits first). A typed name beats a SMILES or a
    library name as the display name; sources are joined."""
    merged = {}
    for cand in (c for group in groups for c in group):
        have = merged.get(cand.key)
        if have is None:
            merged[cand.key] = replace(cand, sources=list(cand.sources), flags=list(cand.flags))         # the inputs are not changed
            continue
        for s in cand.sources:
            if s not in have.sources:
                have.sources.append(s)
        if cand.library and not have.library:
            have.library, have.n_spectra = cand.library, cand.n_spectra
        if cand.user_text and not have.user_text:
            have.user_text = cand.user_text
        if have.parent_smiles is None:
            have.parent_smiles = cand.parent_smiles
        if cand.user_text and cand.name != cand.parent_smiles:          # the user typed a name (a SMILES line is named by its SMILES)
            have.name, have.name_note = cand.name, ""
    return list(merged.values())
