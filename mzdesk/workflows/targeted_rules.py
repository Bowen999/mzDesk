"""Decision rules of the Targeted method workflow (mzdesk/workflows/targeted_method.py). Pure functions on plain data: no engine, no network.

The workflow chains the existing engines (ESI response, GC amenability, separation modes, adducts, detectors, published methods); each engine
answers one question. The rules below only combine those answers into the lines of the method sheet. Every constant is in this file.

Can it be measured (decide)        ESI response = the better of ESI+ / ESI-; GC = the Platform page's GC-MS verdict.
  instrument Any    ESI strong or moderate           Yes       LC-MS
                    else GC "yes" (volatile as is)   Yes       GC-MS
                    else GC "after derivatization"   Possible  GC-MS after derivatization
                    else a non-MS detector rated good (RI never counts: it is always "possible")   Possible   HPLC with that detector, no MS
                    else                             Unlikely  (ESI weak and not GC-amenable: APCI or GC is the next thing to try)
  instrument LC-MS  strong / moderate -> Yes, weak -> Unlikely (judged on LC-MS only; the reason says what the other platform would give)
  instrument GC-MS  yes -> Yes, after derivatization -> Possible, no -> Unlikely (judged on GC-MS only)
  "Possible" = a route exists but needs an extra step (derivatization) or gives up the MS; "Yes" = a direct route. A missing input never
  gives "Unlikely": the verdict is then None (the sheet says "not available").
  The lipid / metabolite tiles of the Platform page only name the LC-MS workflow (lipidomics, metabolomics); they do not change the verdict.

Column (choose_column)  best separation mode = the highest rating among RP, HILIC, mixed-mode, ion-pair (ties in that order; ion-pair last
  because ion-pair reagents are not MS-friendly). LC: a published exact-match method (best documented first, MS detector preferred, the class
  the rules rate higher first) is used when it shows the compound is retained: a retention time of at least MIN_RETAINED_RT minutes (an earlier
  peak is within about twice a typical UHPLC dead time, i.e. unretained: glucose on a C18 gradient). A method without a retention time (JP18)
  counts when the rules rate its column class good or possible. Published retention beats the predicted logD (caffeine runs on C18 although
  the rules rate RP poor). Without such a method: the typical column of the best mode (TYPICAL_COLUMN).
  GC: a published exact-match GC column when the compound runs without derivatization; otherwise the GC_COLUMN default.

pH (ph_choice)  the pH used by all steps is the user's pH or the default 2.7; the row says whether it is inside the robust range of the Separation page
  (every ionizable site > 99 % one form) and, when it is not, which part of the robust range is practical for silica columns (PRACTICAL_PH).

Main ion (pick_main_ion)  the top-scoring ion of the recommended polarity; for "both" (and "poor") the higher-scoring ion of the two polarities.
Derivatization  LC-MS: not needed when the ESI response in the chosen polarity is strong / moderate, else the first open CIL channel.
                GC-MS: from the GC verdict (yes / after derivatization / no) and the reagent the Platform page shows.
"""
from __future__ import annotations

import re
from collections import namedtuple
from dataclasses import dataclass

MINUS = "−"
VERDICT_YES, VERDICT_POSSIBLE, VERDICT_UNLIKELY = "Yes", "Possible", "Unlikely"
ROUTE_LC, ROUTE_GC, ROUTE_NONMS, ROUTE_NONE = "LC-MS", "GC-MS", "non-MS", "none"
ESI_WORKABLE = ("strong", "moderate")            # an ESI response at least this good counts as measurable by LC-MS
ESI_RANK = {"weak": 0, "moderate": 1, "strong": 2}
LEVEL_RANK = {"good": 2, "possible": 1, "poor": 0}
MODE_ORDER = ("RP", "HILIC", "Mixed-mode / ion exchange", "Ion-pair RP")
TYPICAL_COLUMN = {
    "RP": "C18 (reversed phase)",
    "HILIC": "Amide or zwitterionic HILIC",
    "Ion-pair RP": "C18 with an ion-pair reagent (not MS-friendly)",
    "Mixed-mode / ion exchange": "Mixed-mode reversed-phase / ion-exchange",
}
CHIRAL_COLUMN = "Chiral column"
GC_COLUMN = "5 % phenyl-methylpolysiloxane (DB-5ms type)"
MIN_RETAINED_RT = 2.0                            # min: a published retention time below this = not retained (UHPLC dead time is about 0.6-1 min)
PRACTICAL_PH = (2.0, 10.0)                       # pH range silica-based LC-MS columns tolerate
MIN_HINT_WIDTH = 0.5                             # a robust pH window narrower than this (pH units) is not worth suggesting
PH_TOLERANCE = 0.5                               # an additive within this of the used pH is "consistent" with it
ADDITIVE_TEXT = {"Formic acid": "0.1 % formic acid", "Acetic acid": "0.1 % acetic acid", "Ammonium formate": "5-10 mM ammonium formate",
                 "Ammonium acetate": "5-10 mM ammonium acetate", "None": "no additive"}
NON_MS_ORDER = ("UV/PDA", "FLD", "CAD", "ELSD", "RI")
NOT_AVAILABLE = "not available"


# ---------------------------------------------------------------- small helpers
def esi_best(pos, neg):
    """The better of the two ESI response levels (strong / moderate / weak), or None when neither is known."""
    known = [x for x in (pos, neg) if x in ESI_RANK]
    return max(known, key=ESI_RANK.get) if known else None


def polarity_char(polarity):
    """The engine's polarity word ('+', '−' or '-', 'both', 'poor') -> '+', '-', 'both', 'poor' (None when unknown)."""
    if polarity in ("+", "-", MINUS):
        return "-" if polarity in ("-", MINUS) else "+"
    return polarity if polarity in ("both", "poor") else None


def esi_label(sign):
    """'+' -> 'ESI+', '-' -> 'ESI−'."""
    return "ESI+" if sign == "+" else f"ESI{MINUS}"


def level_in(sign, pos, neg):
    """ESI response level in the chosen polarity ('+' or '-')."""
    return pos if sign == "+" else neg


def omics_label(lipid, metabolite):
    """Name of the LC-MS workflow the Platform tiles point to ('' when the compound is in neither range)."""
    if lipid and metabolite:
        return "lipidomics and metabolomics"
    return "lipidomics" if lipid else "metabolomics" if metabolite else ""


def range_text(ranges):
    """[(lo, hi), ...] -> '0 or 4.1-14'."""
    return " or ".join(f"{lo:g}–{hi:g}" if hi > lo else f"{lo:g}" for lo, hi in ranges)


# ---------------------------------------------------------------- can it be measured
@dataclass(frozen=True)
class Decision:
    """verdict: Yes / Possible / Unlikely / None (not available). route: the platform (ROUTE_*). platform: the sheet's label.
    reason: one short sentence for the Basis column. ms_route: the ion route the MS steps follow (ROUTE_LC or ROUTE_GC)."""
    verdict: str | None
    route: str
    platform: str
    reason: str
    ms_route: str


def _lc_label(omics):
    return f"LC-MS ({omics})" if omics else "LC-MS"


def decide(instrument, esi, gc, non_ms=None, omics=""):
    """Can it be measured, and on which platform. esi: best ESI level or None; gc: 'yes' / 'after derivatization' / 'no' or None;
    non_ms: (detector, suitability) of the best non-MS detector or None; omics: omics_label()."""
    lc = None if esi is None else esi in ESI_WORKABLE
    gc_text = {"yes": "volatile as is", "after derivatization": "needs derivatization", "no": "not GC-amenable"}.get(gc)
    other_gc = f"GC-MS would be {gc_text}" if gc_text else ""
    other_lc = f"LC-MS would suit (ESI {esi})" if lc else ""

    if instrument == "LC-MS":
        if lc is None:
            return Decision(None, ROUTE_LC, NOT_AVAILABLE, "ESI response not available", ROUTE_LC)
        if lc:
            return Decision(VERDICT_YES, ROUTE_LC, _lc_label(omics), f"ESI response {esi}; judged on LC-MS only", ROUTE_LC)
        return Decision(VERDICT_UNLIKELY, ROUTE_LC, "LC-MS (ESI weak)",
                        "; ".join(x for x in ("ESI response weak on LC-MS; consider APCI", other_gc) if x), ROUTE_LC)

    if instrument == "GC-MS":
        if gc is None:
            return Decision(None, ROUTE_GC, NOT_AVAILABLE, "GC verdict not available", ROUTE_GC)
        if gc == "yes":
            return Decision(VERDICT_YES, ROUTE_GC, "GC-MS", "volatile as is; judged on GC-MS only", ROUTE_GC)
        if gc == "after derivatization":
            return Decision(VERDICT_POSSIBLE, ROUTE_GC, "GC-MS (derivatized)", "GC-MS works after derivatization; judged on GC-MS only", ROUTE_GC)
        return Decision(VERDICT_UNLIKELY, ROUTE_GC, "GC-MS (not GC-amenable)",
                        "; ".join(x for x in ("not volatile or thermally labile", other_lc) if x), ROUTE_GC)

    # Any: the best platform overall
    if lc:
        also = {"yes": "; GC-MS also works", "after derivatization": "; GC-MS needs derivatization"}.get(gc, "")
        return Decision(VERDICT_YES, ROUTE_LC, _lc_label(omics), f"ESI response {esi}{also}", ROUTE_LC)
    if gc == "yes":
        return Decision(VERDICT_YES, ROUTE_GC, "GC-MS", "ESI response weak; volatile as is" if lc is False else "volatile as is; ESI response not available", ROUTE_GC)
    if gc == "after derivatization":
        return Decision(VERDICT_POSSIBLE, ROUTE_GC, "GC-MS (derivatized)", "ESI response weak; GC-MS works after derivatization" if lc is False
                        else "GC-MS works after derivatization; ESI response not available", ROUTE_GC)
    if non_ms and non_ms[1] == "good" and non_ms[0] != "RI":
        return Decision(VERDICT_POSSIBLE, ROUTE_NONMS, f"HPLC-{non_ms[0]} (no MS)", f"no MS route (ESI weak, not GC-amenable); {non_ms[0]} is suitable", ROUTE_LC)
    if lc is False and gc == "no":
        return Decision(VERDICT_UNLIKELY, ROUTE_NONE, "None suitable by ESI or GC", "ESI response weak and not GC-amenable; consider APCI", ROUTE_LC)
    return Decision(None, ROUTE_NONE, NOT_AVAILABLE, "ESI response or GC verdict not available", ROUTE_LC)


# ---------------------------------------------------------------- separation: best mode and column
def mode_levels(rows):
    """separation_rules.mode_rows -> {mode: level} for the four LC modes (RP, HILIC, Ion-pair RP, Mixed-mode / ion exchange)."""
    return {mode: rec for mode, rec, _ in rows if mode in MODE_ORDER}


def best_mode(levels):
    """-> (mode, level): the highest-rated mode, ties in MODE_ORDER. When nothing is rated above poor, RP is still returned (level 'poor')."""
    if not levels:
        return None, None
    mode = max(MODE_ORDER, key=lambda m: (LEVEL_RANK.get(levels.get(m), -1), -MODE_ORDER.index(m)))
    return mode, levels.get(mode)


_GC = re.compile(r"^GC\b|^GC-", re.I)
_CHIRAL = re.compile(r"chiral|crownpak|cyclobond", re.I)
_IEX = re.compile(r"exchange|mixed[- ]mode|\b(?:WAX|WCX|SAX|SCX)\b", re.I)
_HILIC = re.compile(r"HILIC|\bamide\b|\bZIC\b|\(L(?:3|68|109|114|122)\)", re.I)
_RP = re.compile(r"\bC(?:18|8|4)\b|phenyl|polymer-RP|\bT3\b|\bPFP\b|\(L(?:1|7|10|11|43)\)", re.I)
CLASS_MODE = {"RP": "RP", "HILIC": "HILIC", "IEX": "Mixed-mode / ion exchange"}


def column_class(text):
    """Class of a published column text: 'GC', 'chiral', 'IEX', 'HILIC', 'RP' or None (unknown, or normal phase / SEC)."""
    t = text or ""
    if _GC.search(t):
        return "GC"
    for name, pattern in (("chiral", _CHIRAL), ("IEX", _IEX), ("HILIC", _HILIC), ("RP", _RP)):
        if pattern.search(t):
            return name
    return None


def column_kind(cls, text=""):
    """Short name of the column for the statement: 'C18', 'reversed-phase', 'HILIC', 'mixed-mode', 'chiral' or 'GC'."""
    if cls == "RP":
        return "C18" if re.search(r"\bC18\b|\(L1\)|\bT3\b", text or "") else "reversed-phase"
    return {"HILIC": "HILIC", "IEX": "mixed-mode", "chiral": "chiral", "GC": "GC"}.get(cls, "")


ColumnChoice = namedtuple("ColumnChoice", "text kind source method note")      # source: 'published' | 'rule' | None


def _is_ms(m):
    return (getattr(m, "detector", "") or "").strip().upper().startswith("MS")


def retention_min(m):
    """Longest retention time (min) of a published method ('20.37', '1.06 / 1.33'), None when it has none (JP18)."""
    times = [float(x) for x in re.findall(r"\d+(?:\.\d+)?", str(getattr(m, "rt", None) or ""))]
    return max(times) if times else None


def choose_column(route, levels, exact_methods, gc_direct=True, n_stereo=0):
    """Column for the platform. route: ROUTE_LC or ROUTE_GC; levels: {mode: level} (None when the Separation step is unavailable);
    exact_methods: known_methods Method objects of the exact match (attributes column, detector, rt); gc_direct: the compound runs on GC without
    derivatization. -> ColumnChoice(text, kind, source, method, note); text None = not available."""
    if route == ROUTE_GC:
        cands = [m for m in exact_methods if column_class(m.column) == "GC"] if gc_direct else []
        cands.sort(key=lambda m: not _is_ms(m))
        if cands:
            return ColumnChoice(cands[0].column, "GC", "published", cands[0], "")
        return ColumnChoice(GC_COLUMN, "GC", "rule", None, "" if gc_direct else "derivatives run on the same phase")

    mode, level = best_mode(levels)
    chiral = f"; {CHIRAL_COLUMN.lower()} only for enantiomer separation" if n_stereo else ""
    ranked = []
    for i, m in enumerate(exact_methods):
        cls = column_class(m.column)
        if cls not in CLASS_MODE:
            continue
        rank = 1 if levels is None else LEVEL_RANK.get(levels.get(CLASS_MODE[cls]), 0)      # no Separation step: any recognised LC class
        rt = retention_min(m)
        if (rt >= MIN_RETAINED_RT) if rt is not None else rank >= 1:
            ranked.append((-rank, not _is_ms(m), i, m, cls))
    if ranked:
        _, _, _, m, cls = min(ranked, key=lambda r: r[:3])
        return ColumnChoice(m.column, column_kind(cls, m.column), "published", m, chiral.lstrip("; "))
    if mode is None:
        return ColumnChoice(None, "", None, None, "")
    cls = {"RP": "RP", "Ion-pair RP": "RP", "HILIC": "HILIC", "Mixed-mode / ion exchange": "IEX"}[mode]
    note = ("no mode is rated well; reversed phase is the default" if level == "poor" else "") + chiral
    return ColumnChoice(TYPICAL_COLUMN[mode], "C18" if cls == "RP" else column_kind(cls), "rule", None, note.lstrip("; "))


# ---------------------------------------------------------------- pH
PhChoice = namedtuple("PhChoice", "ph inside text")


def ph_choice(ph, given, ranges, n_sites=1):
    """pH row. ph: the pH used by every step; given: the user typed it; ranges: robust pH ranges [(lo, hi)] of the Separation page (None when the
    step is unavailable); n_sites: number of ionizable sites. -> PhChoice(ph, inside True / False / None, basis text without the step number)."""
    who = "pH" if given else "default pH"
    if ranges is None:
        return PhChoice(ph, None, f"{who} {ph:g}; the robust range is not available")
    if not n_sites:
        return PhChoice(ph, True, f"{who} {ph:g}; no ionizable site, retention does not depend on pH")
    if not ranges:
        return PhChoice(ph, False, f"{who} {ph:g}; no pH keeps every site one form, so buffer well")
    inside = any(lo - 1e-9 <= ph <= hi + 1e-9 for lo, hi in ranges)
    if inside:
        return PhChoice(ph, True, f"{who} {ph:g} is inside the robust range pH {range_text(ranges)}")
    lo_p, hi_p = PRACTICAL_PH
    practical = [(max(lo, lo_p), min(hi, hi_p)) for lo, hi in ranges if min(hi, hi_p) - max(lo, lo_p) >= MIN_HINT_WIDTH]
    hint = f"; pH {range_text(practical)} is more robust" if practical else ""
    return PhChoice(ph, False, f"{who} {ph:g} is outside the robust range pH {range_text(ranges)}: retention is pH-sensitive, buffer well{hint}")


# ---------------------------------------------------------------- main ion
def pick_main_ion(polarity, pos, neg):
    """Main ion for the recommended polarity. pos / neg: the top predicted ion of each polarity (attributes adduct, score, likelihood) or None.
    polarity: '+', '-' ('−'), 'both' or 'poor'. -> (sign '+' / '-' / None, ion or None, note)."""
    pol = polarity_char(polarity) or "both"           # an unknown polarity is treated as both
    if pol in ("+", "-"):
        ion = pos if pol == "+" else neg
        return (pol, ion, "") if ion is not None else (pol, None, f"no ion expected in {esi_label(pol)}")
    if pos is None and neg is None:
        return None, None, "no ion expected in either polarity"
    sp, sn = (pos.score if pos else -1.0), (neg.score if neg else -1.0)
    sign = "+" if sp >= sn else "-"
    ion = pos if sign == "+" else neg
    which = "both polarities respond" if pol == "both" else "ESI response is weak in both polarities"
    other = neg if sign == "+" else pos
    versus = f" vs {other.score:g} in {esi_label('-' if sign == '+' else '+')}" if other is not None else ""
    return sign, ion, f"{which}; {esi_label(sign)} has the higher ion score ({ion.score:g}{versus})"


# ---------------------------------------------------------------- derivatization
Deriv = namedtuple("Deriv", "kind text basis")          # kind: 'not needed' | 'suggested' | 'needed' | 'none' | 'not feasible'


def lc_derivatization(level, channels):
    """LC-MS. level: ESI response in the chosen polarity; channels: [(channel name, reagent name)] of the open CIL channels, best first."""
    if level is None:
        return None
    if level in ESI_WORKABLE:
        return Deriv("not needed", "Not needed", f"ESI response {level}")
    if channels:
        channel, reagent = channels[0]
        more = f"; also {', '.join(c for c, _ in channels[1:])}" if len(channels) > 1 else ""
        return Deriv("suggested", f"CIL: {reagent} ({channel})", f"ESI response weak; labelling the {channel} group adds an ionizable tag{more}")
    return Deriv("none", "None available", "ESI response weak and no amine, phenol, hydroxyl, carboxyl or carbonyl group to label")


def gc_derivatization(gc, reagent, volatility=""):
    """GC-MS. gc: 'yes' / 'after derivatization' / 'no' (None unknown); reagent: the Platform page's reagent name (MeOx + TMS or TMS)."""
    if gc == "yes":
        return Deriv("not needed", "Not needed", volatility or "volatile as is")
    if gc == "after derivatization":
        return Deriv("needed", f"Needed: {reagent}", volatility or "active H or carbonyl groups need silylation")
    if gc == "no":
        return Deriv("not feasible", "Not feasible by GC", volatility or "not volatile or thermally labile")
    return None


# ---------------------------------------------------------------- non-MS detectors
def best_non_ms(rows):
    """non_ms_rules detectors rows [[name, suitability good / possible / poor, reason], ...] -> (name, suitability, reason) of the best one
    (ties in NON_MS_ORDER), or None."""
    rows = [r for r in rows if r and r[0] in NON_MS_ORDER]
    if not rows:
        return None
    best = max(rows, key=lambda r: (LEVEL_RANK.get(r[1], -1), -NON_MS_ORDER.index(r[0])))
    return best[0], best[1], best[2] if len(best) > 2 else ""


# ---------------------------------------------------------------- published methods
def methods_count_text(n_jp, n_report, n_similar=0):
    """'151 exact (JP18 1, RepoRT 150)'. n_jp / n_report: the true totals of the search (known_methods Result.exact_totals), not only the listed rows.
    No exact match -> 'None for this compound' (+ the number of similar compounds that have methods)."""
    n = n_jp + n_report
    if not n:
        return "None for this compound" + (f"; methods of {n_similar} similar compound{'s' if n_similar != 1 else ''}" if n_similar else "")
    parts = ([f"JP18 {n_jp}"] if n_jp else []) + ([f"RepoRT {n_report}"] if n_report else [])
    return f"{n} exact ({', '.join(parts)})"


def gc_methods_text(n_gc, n_lc):
    """Published methods on the GC route. n_gc: exact methods that use a GC column; n_lc: the other exact methods (LC, mostly RepoRT), which do not
    help a GC-MS user: '1 GC exact (150 LC methods also exist)', 'None for GC (93 LC methods exist)', 'None for GC'."""
    head = f"{n_gc} GC exact" if n_gc else "None for GC"
    if not n_lc:
        return head
    return f"{head} ({n_lc} LC method{'s' if n_lc != 1 else ''} {'also ' if n_gc else ''}exist{'s' if n_lc == 1 else ''})"


# ---------------------------------------------------------------- GC: the ion that usually dominates the derivative's EI spectrum
ELECTRON = 0.00054857990946     # u
SILYL_LOSS = {"TMS": ("M−15", 15.0234750), "MeOx + TMS": ("M−15", 15.0234750),         # CH3 from a TMS group: C + 3 H
              "TBDMS": ("M−57", 57.0704253)}                                           # C4H9 (tert-butyl) from a TBDMS group: 4 C + 9 H


def silyl_fragment(reagent, mass):
    """The usually stronger EI ion of a silyl derivative: ('M−15', m/z) for TMS and MeOx + TMS (CH3 lost), ('M−57', m/z) for TBDMS (tert-butyl lost);
    mass: the neutral derivative's monoisotopic mass. None for any other reagent. The ion is M+• minus the radical, so its m/z is mass − electron − loss."""
    loss = SILYL_LOSS.get(reagent)
    return None if loss is None or mass is None else (loss[0], mass - ELECTRON - loss[1])
