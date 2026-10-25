"""Rule-based parser: JP 'Operating conditions' block text -> structured fields (annotation_guide semantics).

Fields (31): technique, block_type, xref_target, xref_fields, detector_type, wavelength_nm, column_count, column_material,
column_id_mm, column_length_mm, packing_class, particle_or_film_um, column_temp_c, injection_temp_c, detector_temp_c,
elution_mode, n_mobile_phases, mobile_phase_ph, organic_solvents, gradient_table, carrier_gas, split_ratio,
flow_rate_value, flow_rate_unit, flow_target_compound, flow_target_rt_min, run_time_min, run_time_rt_multiple,
plus block_id / monograph / entry / section (carried by the caller) and `notes`.
Numbers are JSON numbers; null = not stated in the block (cross-referenced values are filled in later by merge.py).
"""
from __future__ import annotations

import re

NUM = r"(\d+(?:\.\d+)?)"
FIELDS = ["technique", "block_type", "xref_target", "xref_fields", "detector_type", "wavelength_nm", "column_count",
          "column_material", "column_id_mm", "column_length_mm", "packing_class", "particle_or_film_um", "column_temp_c",
          "injection_temp_c", "detector_temp_c", "elution_mode", "n_mobile_phases", "mobile_phase_ph", "organic_solvents",
          "gradient_table", "carrier_gas", "split_ratio", "flow_rate_value", "flow_rate_unit", "flow_target_compound",
          "flow_target_rt_min", "run_time_min", "run_time_rt_multiple", "notes"]
# fields that each cross-referenced group fills in (merge.py)
XREF_GROUPS = {
    "detector": ["detector_type", "wavelength_nm", "detector_temp_c", "_wl_options"],
    "column": ["column_count", "column_material", "column_id_mm", "column_length_mm", "packing_class", "particle_or_film_um"],
    "column_temperature": ["column_temp_c", "injection_temp_c"],
    "mobile_phase": ["elution_mode", "n_mobile_phases", "mobile_phase_ph", "organic_solvents", "gradient_table",
                     "carrier_gas", "split_ratio"],
    "flow_rate": ["flow_rate_value", "flow_rate_unit", "flow_target_compound", "flow_target_rt_min"],
}
RUN_TIME = ["run_time_min", "run_time_rt_multiple"]

PACKING = [  # (pattern, class), checked in order. Bonded / modified silica gels never fall back to "silica" (see PLAIN_SILICA_RE).
    (r"triacontylsil", "other"), (r"glycol etherified", "diol"),
    (r"octade[a-z]{0,8}sil[a-z]*|octadecyl[- ]?silyl|octadecyl", "C18"), (r"octylsil|octyl[- ]?silyl", "C8"), (r"butylsil", "C4"),
    (r"pentafluorophenyl", "PFP"), (r"phenyl(?:ated|sil|hexyl)|phenylsilanized", "phenyl"), (r"cyanopropylsil|cyanosil", "CN"),
    (r"aminopropylsil|aminosil", "NH2"), (r"\bdiol", "diol"),
    (r"diethylaminoethyl|quaternary ammonium|anion[- ]exchange|strong(?:ly)? basic|weakly basic|aminated", "anion-exchange"),
    (r"sulfonic acid|sulfopropyl|sulfoethyl|sulfonated|cation[- ]exchange|strongly acidic|weakly acidic|slightly acidic|acidic ion|strong acid|carboxymethyl", "cation-exchange"),
    (r"gel filtration|size[- ]exclusion|for gel|hydrophilic silica gel for|hydrophilic polymer", "SEC"),
    (r"chiral|optical resolution|optically active|albumin|glycoprotein|ovomucoid|cellulose|amylose|cyclodextrin|teicoplanin", "chiral"),
    (r"st[yi]rene[- ]divinylbenzene|stylene|polystyrene", "polymer-RP"),
]
# only the plain "silica gel for liquid chromatography" is class "silica"; anything bonded / modified ("carbamoyl group bound silica gel",
# "trimethylsilanized silica gel", "sulfonamide group bound to hexadecylsilanized silica gel") is "other"
PLAIN_SILICA_RE = re.compile(r"(?:packed with|filled with)\s+(?:an?\s+)?(?:(?:porous|spherical|totally porous|fully porous|superficially porous)\s+)*silica gel\b")
GC_PACKING = [
    (r"polyethylene glycol|peg\b|carbowax", "GC-PEG"),
    (r"(?<!\d)6\s*%\s*cyanopropyl\s*phenyl", "GC-6%cyanopropylphenyl"),
    (r"cyanopropyl", "GC-other"),                                # 14 % cyanopropylphenyl-86 % dimethyl silicone polymer etc.
    (r"(?<!\d)5\s*% (?:di)?phenyl", "GC-5%phenyl"),
    (r"\d{2}\s*% (?:di)?phenyl", "GC-other"),                  # 35 % / 50 % phenyl-methyl polysiloxane
    (r"dimethylpolysiloxane|dimethyl silicone|methyl silicone|methylsilicone|methylpolysiloxane", "GC-dimethylpolysiloxane"),
]
DETECTORS = [(r"photodiode array|diode array", "PDA"), (r"ultraviolet|uv\b", "UV"), (r"visible", "visible"),
             (r"flame[- ]ionization", "FID"), (r"thermal[- ]conductivity", "TCD"), (r"fluorescen|fluorophotometer|fluorometr|spectrofluorometer", "fluorescence"),
             (r"refract", "RI"), (r"mass spectro", "MS"), (r"(?<!thermal )(?<!thermal-)conductivity", "conductivity"),
             (r"electron capture", "ECD"), (r"flame photometric|nitrogen[- ]phosph", "other")]
SOLVENTS = ["acetonitrile", "methanol", "ethanol", "2-propanol", "1-propanol", "tetrahydrofuran", "acetone", "hexane",
            "heptane", "dichloromethane", "chloroform", "ethyl acetate", "diethyl ether", "1-butanol", "2-butanol",
            "dimethylsulfoxide", "N,N-dimethylformamide", "cyclohexane", "1,4-dioxane", "isooctane"]
XREF_FIELD_WORDS = [("detector", "detector"), ("column temperature", "column_temperature"),
                    ("mobile phase", "mobile_phase"), ("flow rate", "flow_rate")]
LINE_LABEL_RE = re.compile(r"^([A-Z][A-Za-z()\-/' ,]{0,220}?):(?:\s|$)")
XREF_VAL_RE = re.compile(r"^(?:\(.*?\)\s*)?(?:Proceed|Perform)(?:\s+as directed)?\s+(?:in|the|under|as)\b", re.I)
PREAMBLE_XREF_RE = re.compile(r"^\s*(?:Proceed|Perform)\b", re.I)
INLINE_LABEL_XREF_RE = re.compile(r"(?<=\.)\s+([A-Z][A-Za-z()\-/' ,]{0,200}?):\s+(?=(?:Proceed|Perform)\b)")


def first_temp(val: str):
    """First temperature in '85 ± 2 °C' / 'about 40 °C' -> 85 / 40 (the tolerance is not the value)."""
    m = re.search(NUM + r"\s*(?:±\s*" + NUM + r"\s*)?°C", val)
    return fnum(m.group(1)) if m else None


def fnum(s):
    v = float(s)
    return int(v) if v.is_integer() else v


def join_lines(lines: list[str]) -> str:
    return re.sub(r"\s+", " ", " ".join(lines)).strip()


def split_fields(text: str) -> tuple[str, list[tuple[str, str]]]:
    """(preamble, [(label, value)]) of an Operating-conditions block (header line excluded). The preamble is any text
    before the first 'Label:' line, e.g. 'Proceed as directed in the operating conditions in the Assay.'"""
    lines = text.split("\n")[1:]
    merged: list[str] = []
    for ln in lines:      # wrapped label: "Detector, column, column temperature, mobile phase, and" + "flow rate: Proceed ..."
        if merged and ":" not in merged[-1] and re.match(r"^[A-Z][A-Za-z ,()\-/']*$", merged[-1]) \
                and not merged[-1].endswith(".") and re.match(r"^[a-z][A-Za-z ,()\-/']*(?::(?:\s|$)|$)", ln):
            merged[-1] += " " + ln
        else:
            merged.append(ln)
    fields: list[list[str]] = []
    pre: list[str] = []
    for ln in merged:
        m = LINE_LABEL_RE.match(ln)
        if m and not ln.startswith(("Proceed", "Perform", "Time after")) and len(m.group(1).split()) <= 40:
            fields.append([m.group(1).strip(), ln[m.end():].strip()])
        elif fields:
            fields[-1][1] += " " + ln.strip()
        else:
            pre.append(ln)
    out: list[tuple[str, str]] = []
    for l, v in fields:                  # "Detector: A differential refractometer. Column, ...: Proceed as directed ..."
        m = INLINE_LABEL_XREF_RE.search(v)
        if m:
            out.append((l, v[:m.start()].strip()))
            out.append((m.group(1).strip(), v[m.end():].strip()))
        else:
            out.append((l, v.strip()))
    return " ".join(pre).strip(), out


def label_groups(label: str) -> set[str]:
    """Which condition groups a (possibly composite) label like 'Detector, column, and column temperature' names."""
    ll = re.sub(r"^the\s+", "", label.lower())
    groups = set()
    for tok in re.split(r",|\band\b(?! b\b)", ll):
        t = tok.strip()
        if not t:
            continue
        if t.startswith(("detector", "apparatus", "reaction coil", "reaction temperature", "reagent", "reaction reagent",
                         "detector temperature", "suppressor", "chemical reaction")):
            groups.add("detector")
        elif t.startswith(("column temperature", "temperature of column", "injection port", "temperature of injection")):
            groups.add("column_temperature")
        elif t.startswith(("column", "columns", "selection of column", "pretreatment")):
            groups.add("column")
        elif t.startswith(("mobile phase", "mobile phases", "phase", "flowing of", "changing of", "switching of", "carrier gas", "split ratio",
                           "gradient")):
            groups.add("mobile_phase")
        elif t.startswith(("flow rate", "flow of", "flow rate of")):
            groups.add("flow_rate")
        elif t.startswith(("time span", "time of measurement")):
            groups.add("time_span")
    # "flow rate of the mobile phase" / "mobile phases A and B" must not leak extra groups
    if re.search(r"flow rate of (?:the )?mobile phase", ll):
        groups.add("flow_rate")
    return groups


def gradient_table(text: str):
    """Rows of a solvent-gradient table: [t_start, t_end, A_start, A_end, B_start, B_end, (C_start, C_end ...)]."""
    rows = []
    lines = text.split("\n")
    in_flow_table = False
    for ln in lines:
        if re.match(r"^Time after", ln):
            in_flow_table = False
        if re.search(r"\(mL/min(?:ute)?\)", ln) and not re.search(r"vol%", ln):
            in_flow_table = True
        if in_flow_table:
            continue
        m = re.match(r"^\s*" + NUM + r"\s*[–-]\s*" + NUM + r"\*?\s{2,}(.+)$", ln)
        if not m:
            continue
        cells = re.split(r"\s{2,}", m.group(3).strip())
        vals = []
        ok = True
        for c in cells:
            c = c.strip()
            a = re.match(NUM + r"\s*→\s*" + NUM + r"$", c)
            if a:
                vals += [fnum(a.group(1)), fnum(a.group(2))]
            elif re.fullmatch(NUM, c):
                vals += [fnum(c), fnum(c)]
            else:
                ok = False
                break
        if ok and len(vals) in (4, 6, 8, 10):
            rows.append([fnum(m.group(1)), fnum(m.group(2))] + vals)
    return rows or None


def detector_info(val: str, lv: str):
    found = [(m.start(), t) for p, t in DETECTORS for m in [re.search(p, lv)] if m]
    dtype = min(found)[1] if found else "other"
    wl = None
    if dtype == "fluorescence":
        m = re.search(r"excitation\s+wave\s?length[^:\d]*[:\s]\s*" + NUM, lv)
        wl = fnum(m.group(1)) if m else None
    if wl is None:
        m = re.search(r"wave\s?length[s]?\s*(?:of\s*)?[:\s]?\s*(?:a constant wavelength between\s*)?\[?\s*" + NUM, lv)
        if m:
            wl = fnum(m.group(1))
    return dtype, wl


def column_info(val: str, lv: str, gc: bool) -> dict:
    out = {}
    m = re.search(r"\b(two|three|four|2|3|4)\b(?![.\d])[^.]*?columns?", lv)
    out["column_count"] = {"two": 2, "three": 3, "four": 4, "2": 2, "3": 3, "4": 4}[m.group(1)] if m else 1
    m = re.search(r"(stainless steel|fused silica|glass|quartz)", lv)
    out["column_material"] = m.group(1) if m else None
    dims = []
    for md in re.finditer(NUM + r"\s*mm\b(?:\s*\(or [^)]*\))?\s*(?:in\s*)?(?:inside\s*)?diam[a-z.]*\s*and\s*(?:about\s*)?" + NUM + r"\s*(mm|cm|m)\b\s*in\s*length", lv):
        pre, post = lv[max(0, md.start() - 45):md.start()], lv[md.end():md.end() + 130]
        kind = "tube" if re.search(r"\btubes?\b|\btubing\b", pre) and not re.search(r"column", pre[-25:]) else \
            "guard" if re.search(r"as (?:a|the) (?:guard column|pre-?column)", post[:90]) else "column"
        dims.append((fnum(md.group(1)), fnum(float(md.group(2)) * {"mm": 1, "cm": 10, "m": 1000}[md.group(3)]), kind))
    cols = [d for d in dims if d[2] != "tube"]
    if cols and (len(cols) > 1 or len(dims) > len(cols)):
        paired = re.search(r"(?:pre-?column|guard column)[^.]{0,30}(?:separation|analytical) column[^.]{0,15}respectively|"
                           r"(?:pre-?column|guard column) and the separation column", lv)
        main = [d for d in cols if d[2] == "column"] or cols
        best = cols[-1] if paired else max(main, key=lambda d: d[1])
        if len(cols) > 1 or len({(d[0], d[1]) for d in dims}) > 1:
            out["column_id_mm"], out["column_length_mm"] = best[0], best[1]
            out["column_count"] = len(cols)
            skipped = [d for d in dims if d[2] == "tube"]
            out["_column_note"] = ("columns of different sizes (" + "; ".join(f"{a} x {b} mm" for a, b, _ in cols) + "): " +
                                   ("separation column (last) reported" if paired else "analytical column (longest) reported") +
                                   ("; connecting tube ignored" if skipped else ""))
        else:
            out["column_id_mm"], out["column_length_mm"] = best[0], best[1]
    else:
        m = re.search(NUM + r"\s*(?:(?:to|–|-)\s*" + NUM + r"\s*)?mm\b(?:\s*\(or [^)]*\))?\s*(?:in\s*)?(?:inside\s*)?diam", lv) or \
            re.search(r"(?:inside\s*)?diameters?\s*(?:of|:)?\s*(?:about\s*)?" + NUM + r"\s*mm", lv)
        out["column_id_mm"] = fnum(m.group(1)) if m else None
        m = re.search(NUM + r"\s*(?:(?:to|–|-)\s*" + NUM + r"\s*)?(mm|cm|m)\b\s*(?:in length|long)", lv)
        if m:
            out["column_length_mm"] = fnum(float(m.group(1)) * {"mm": 1, "cm": 10, "m": 1000}[m.group(3)])
        elif re.search(r"lengths?\s*(?:of|:)?\s*(?:about\s*)?" + NUM + r"\s*(mm|cm|m)\b", lv):
            m = re.search(r"lengths?\s*(?:of|:)?\s*(?:about\s*)?" + NUM + r"\s*(mm|cm|m)\b", lv)
            out["column_length_mm"] = fnum(float(m.group(1)) * {"mm": 1, "cm": 10, "m": 1000}[m.group(2)])
        else:
            m = re.search(NUM + r"\s*mm\s*[×x]\s*" + NUM + r"\s*(mm|cm|m)\b", lv)
            if m:
                out["column_id_mm"] = out["column_id_mm"] or fnum(m.group(1))
                out["column_length_mm"] = fnum(float(m.group(2)) * {"mm": 1, "cm": 10, "m": 1000}[m.group(3)])
            else:
                out["column_length_mm"] = None
    pack = None
    for p, c in (GC_PACKING if gc else []) + PACKING:
        if re.search(p, lv):
            pack = c
            break
    if pack is None:
        pack = "GC-other" if gc else ("silica" if PLAIN_SILICA_RE.search(lv) else "other")
    elif gc and pack in ("C18", "C8", "C4", "phenyl", "CN", "NH2", "diol", "PFP", "polymer-RP"):
        pack = "GC-other"
    out["packing_class"] = pack
    m = re.search(NUM + r"\s*-?\s*(?:to|–|-)\s*" + NUM + r"\s*-?\s*μm", lv) or \
        re.search(r"(?:layer|film)\s*(?:of\s*)?(?:about\s*)?" + NUM + r"\s*-?\s*μm", lv) or \
        re.search(NUM + r"\s*-?\s*μm\s*(?:thick|layer|in film)", lv) or \
        re.search(r"\(\s*" + NUM + r"\s*μm\s*(?:in )?particle", lv) or \
        re.search(NUM + r"\s*μm\s*(?:\(or [^)]*\)\s*)?(?:in )?(?:film )?thickness", lv) or \
        re.search(r"thickness\s*(?:of)?\s*" + NUM + r"\s*μm", lv) or \
        re.search(r"(?:particle (?:diameter|size)\s*(?:of)?[:\s]*)" + NUM + r"\s*μm", lv) or \
        re.search(NUM + r"\s*μm\s*(?:in )?particle", lv) or \
        re.search(r"in\s*" + NUM + r"\s*μm\s*(?:film )?thick", lv) or \
        re.search(r"\(\s*" + NUM + r"\s*μm\s*\)", lv) or \
        re.search(r"in\s+" + NUM + r"\s*μm\b", lv)
    out["particle_or_film_um"] = fnum(m.group(1)) if m else None
    if not m:
        mt = re.search(NUM + r"\s*(?:to|–|-)\s*" + NUM + r"\s*mm\s+in\s+particle", lv)       # typo in the source: "5 to 10 mm in particle diameter"
        if mt:
            out["particle_or_film_um"] = fnum(mt.group(1))
            out["_column_note"] = "particle size printed in mm in the source, read as um"
    return out


_UNIT_AFTER = r"(?!\d|\.\d|\s*(?:mol/l|mmol/l|mol|%|g\b|mg|ml|μ|nm|mm|cm|°|volumes?|times|in\s+\d))"


def ph_value(mt: str):
    """First pH value stated in lower-cased mobile-phase text: 'adjust to pH 3.0', 'adjust the pH of 1000 mL of water to 3.0',
    'phosphate buffer (pH 7.0)'. Concentrations such as 'pH of 0.05 mol/L ...' are not pH values."""
    cands = []
    for m in re.finditer(r"\bph\s*(?:of\s*)?(?:about\s*|between\s*|to\s*)?" + NUM + _UNIT_AFTER, mt):
        cands.append((m.start(), m.group(1)))
    for m in re.finditer(r"\bph\s+(?:of\s+)?(?:[^.;]|\.\d){0,220}?\bto\s+(?:a\s+ph\s+of\s+|ph\s+)?(?:about\s+|between\s+)?" + NUM + _UNIT_AFTER, mt):
        cands.append((m.start(), m.group(1)))
    for m in re.finditer(r"\bph\b(?:[^.;]|\.\d){0,40}?\bis\s+(?:about\s+)?" + NUM + _UNIT_AFTER, mt):
        cands.append((m.start(), m.group(1)))
    return fnum(min(cands)[1]) if cands else None


_MP_END_RE = re.compile(r"^(?:Flowing of|Changing|Switching|Reaction|Reagent|Flow rate|Flow of|Time span|Detection|Selection|Column|Detector|System)\b")


def table_text(text: str) -> str:
    """Lowercased lines of the mobile-phase composition table that follows the 'Mobile phase:' field ("Ethanol (95)  60.0 mL ..."),
    up to the next condition label; empty when there is no such table."""
    lines = text.split("\n")
    out, on = [], False
    for ln in lines[1:]:
        if re.match(r"^Mobile phase\b", ln):
            on = True
            continue
        if on:
            if _MP_END_RE.match(ln):
                break
            if re.search(r"\d\s*(?:mL|g)\b", ln) and "  " in ln:
                out.append(ln.lower())
    return " ".join(out)


def wavelength_options(val: str) -> list[tuple]:
    """'231 nm for A and B; 254 nm for C' -> [(231, 'a and b'), (254, 'c')]; nothing unless every wavelength carries a label."""
    ms = list(re.finditer(NUM + r"\s*nm", val))
    if len(ms) < 2:
        return []
    opts = []
    for i, m in enumerate(ms):
        end = ms[i + 1].start() if i + 1 < len(ms) else len(val)
        label = re.sub(r"^\s*(?:for|\(|\[)\s*|[\s;,)\]]+$|\s+and\s*$|\s*;\s*$", "", val[m.end():end]).strip(" ;,)(]")
        label = re.sub(r"\s+and$", "", label)
        if not label or re.fullmatch(r"(?:and|or|,)", label):
            return []
        opts.append((fnum(m.group(1)), label.lower()))
    return opts if len({w for w, _ in opts}) > 1 else []


def parse_block(b: dict) -> dict:
    text = b["text"]
    out: dict = {k: None for k in FIELDS}
    out["xref_time_span"] = False
    first = text.split("\n")[0]
    low = text.lower()
    hint = b.get("technique_hint")
    if re.sub(r"[\s—\d]", "", first) != "Operatingconditions":
        out["block_type"] = "prose"
        m = re.search(r"follow the (procedure [A-Z])", text, re.I)
        if m:
            out["xref_target"] = "Procedure " + m.group(1).split()[-1].upper()
        out["technique"] = "GC" if ("gas chromatography" in low or "headspace" in low) else hint
        return out
    preamble, fields = split_fields(text)
    labels = " ".join(l.lower() for l, _ in fields)

    # ---- classification: NMR / other (non-chromatographic) / LC / GC
    if "nuclear magnetic" in low or "resonance frequency" in low or re.search(r"target nucle|pulse angle|spinning|dummy scan", labels):
        out["technique"], out["block_type"] = "NMR", "full"
        return out
    if re.search(r"electrolyte solution|running temperature|migration|capillary tube", labels + " " + low) and not re.search(r"mobile phase|carrier gas", labels):
        out["technique"], out["block_type"] = "other", "full"
        return out
    if re.search(r"viscometer|rotor no|brookfield", low[:300]):
        out["technique"], out["block_type"] = "other", "full"
        return out
    if re.search(r"heating rate|atmospheric gas|temperature range|gel support|polyacrylamide|electrophoresis|spotting|electrolyte|cathode|anode",
                 labels + " " + low[:400]) and not re.search(r"\bcolumn\b|mobile phase|carrier gas", labels):
        out["technique"], out["block_type"] = "other", "full"
        return out
    gc_markers = bool(re.search(r"carrier gas|injection port|split ratio|detector temperature|temperature of injection|sample vaporization", labels + " " + preamble.lower())
                      or re.search(r"flame-ionization|thermal.conductivity|electron capture|for gas chromatography", low))
    if gc_markers:
        out["technique"] = "GC"
    elif re.search(r"mobile phase|column|detector|flow rate|reaction coil", labels + " " + preamble.lower()):
        out["technique"] = "LC"
    elif PREAMBLE_XREF_RE.match(preamble):
        out["technique"] = hint or "LC"                                # whole block is "Proceed as directed in ..."
    else:
        out["technique"], out["block_type"] = "other", "full"          # e.g. polarimeter "Wavelength: 340.458 nm"
        return out
    if out["technique"] == "LC" and hint == "GC" and not re.search(r"mobile phase|mL per minute", labels + low):
        out["technique"] = "GC"
    gc = out["technique"] == "GC"

    groups: set[str] = set()
    xref_text = None
    stated = 0
    mobile: list[str] = []
    stock_text: list[str] = []
    stated_time = False
    if preamble and PREAMBLE_XREF_RE.match(preamble):        # whole-block cross-reference
        groups |= {"detector", "column", "column_temperature", "mobile_phase", "flow_rate", "time_span"}
        xref_text = preamble
    for label, val in fields:
        ll, lv = label.lower(), val.lower()
        if XREF_VAL_RE.match(val):
            groups |= label_groups(label)
            xref_text = xref_text or val
            continue
        stated += 1
        if ll.startswith("detector temperature") or ll.startswith("temperature of detector"):
            out["detector_temp_c"] = first_temp(val)
        elif ll.startswith("interface temperature"):
            pass
        elif ll.startswith("injection port and detector temperature"):
            m = re.findall(NUM + r"\s*(?:±\s*" + NUM + r"\s*)?°C", val)
            if m:
                out["injection_temp_c"] = fnum(m[0][0])
                out["detector_temp_c"] = fnum(m[-1][0])
        elif ll.startswith(("detector", "ditector")):
            out["detector_type"], out["wavelength_nm"] = detector_info(val, lv)
            if out["detector_type"] in ("UV", "visible", "PDA"):
                wo = wavelength_options(val)
                if wo:
                    out["_wl_options"] = wo
            kinds = {t for p, t in DETECTORS if re.search(p, lv)}
            if len(kinds) > 1:
                out["notes"] = "multiple detectors named: " + ", ".join(sorted(kinds))
            m = re.search(r"detector temperature[^:]*?:?[^\d]*" + NUM + r"\s*°C", lv)
            if m and "refract" in lv:
                out["detector_temp_c"] = fnum(m.group(1))
        elif ll.startswith("column temperature"):
            if re.search(r"raise|rate of|program|at beginning|change to|then maintain|increase|\bvari(?:ed|es)\b", lv):
                out["column_temp_c"] = "program"
            elif "room temperature" in lv or "ordinary temperature" in lv:
                out["column_temp_c"] = "room"
            else:
                out["column_temp_c"] = first_temp(val)
        elif ll.startswith(("column", "columns")):
            ci = column_info(val, lv, gc)
            note = ci.pop("_column_note", None)
            out.update(ci)
            if note:
                out["notes"] = (out["notes"] + "; " if out["notes"] else "") + note
        elif ll.startswith(("reaction reagent flow", "reaction regent flow")):
            pass
        elif ll.startswith(("mobile phase", "mobil phase")) and "flow" not in ll:
            if ll.startswith("mobile phase stock"):
                stock_text.append(val)
            else:
                mobile.append(val)
                if re.match(r"mobile phases? [a-d]\b", ll) and not re.match(r"mobile phase [a-d]$", ll):
                    out["_phase_letters"] = re.findall(r"\b([a-d])\b", ll.split("phase", 1)[1])
        elif ll.startswith(("flowing of", "changing of mobile", "switching of mobile", "changing mobile", "switching mobile",
                            "switching mobile phases")):
            if "gradient" in lv:
                out["elution_mode"] = "gradient"
        elif ll.startswith(("flow rate", "flow of mobile", "mobile phase flow")):
            if ll.startswith(("flow rate of reaction", "flow rate of reagent", "flow rate of refresh", "flow rate of wash",
                              "flow rate of atmos", "flow rate of reaction regent")):
                continue
            m2 = re.search(r"between\s*" + NUM + r"\s*mL\s*and\s*" + NUM + r"\s*mL\s*per min", val)
            m = re.search(NUM + r"\s*(?:(?:–|-|to)\s*" + NUM + r")?\s*mL\s*(?:per min|/min)", val)
            if m2:
                out["flow_rate_value"], out["flow_rate_unit"] = fnum((float(m2.group(1)) + float(m2.group(2))) / 2), "mL/min"
            elif m:
                out["flow_rate_value"] = fnum((float(m.group(1)) + float(m.group(2))) / 2) if m.group(2) else fnum(m.group(1))
                out["flow_rate_unit"] = "mL/min"
            m = re.search(NUM + r"\s*cm\s*(?:per second|/s(?:ec(?:ond)?)?\b)", val)
            if m:
                out["flow_rate_value"], out["flow_rate_unit"] = fnum(m.group(1)), "cm/s"
            m = re.search(r"(?:retention|reaction) times? (?:of|for) (?:the (?:peak|peaks) (?:of|showing [a-z ]*of the two peaks of) )?(.+?) "
                          r"(?:is|are|becomes|become|will be|should be) (?:about|between|approximately)?\s*"
                          r"(" + NUM + r"(?:\s*(?:,|and|to|–|-)\s*" + NUM + r")*)\s*minutes?", val)
            if m:
                comp = re.sub(r"^the ", "", m.group(1).strip().lower())
                nums = [float(x) for x in re.findall(NUM, m.group(2))]
                several = bool(re.search(r"times|respectively", val[m.start():m.end() + 20])) or len(nums) > 2 or \
                    (len(nums) == 2 and re.search(r"\d\s*,\s*\d", m.group(2)))
                if several:
                    comp = re.split(r",| and ", comp)[0].strip()      # several compounds / times: keep the first pair
                    comp = re.sub(r"^the ", "", comp)
                out["flow_target_compound"] = comp
                out["flow_target_rt_min"] = fnum((nums[0] + nums[1]) / 2) if len(nums) == 2 and not several else fnum(nums[0])
            if re.search(r"directed in the following table|as shown in the following table|following table", lv) and \
                    out["flow_rate_value"] is None and out["flow_target_rt_min"] is None:
                out["notes"] = (out["notes"] + "; " if out["notes"] else "") + "programmed flow (table)"
        elif ll.startswith(("time span", "time of measurement", "time for measurement", "time span for measurement")):
            stated_time = True
            words = {"twice": 2, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6, "seven": 7, "eight": 8, "ten": 10}
            mult = None
            if "retention time" in lv:
                m = re.search(NUM + r"\s*times?\b", lv) or re.search(r"\b(twice|two|three|four|five|six|seven|eight|ten)\s+times\b", lv)
                if m:
                    mult = fnum(m.group(1)) if re.fullmatch(NUM, m.group(1)) else words[m.group(1)]
                elif re.search(r"\btwice\b", lv):
                    mult = 2
            out["run_time_rt_multiple"] = mult
            if mult is None:
                m = re.search(r"from\s*" + NUM + r"\s*minutes? to\s*" + NUM + r"\s*minutes", lv)
                if m:
                    out["run_time_min"] = fnum(m.group(2))
                else:
                    m = re.search(NUM + r"\s*minutes?\b(?!\s*after (?:the )?elution)", lv)
                    if m and not re.search(r"retention time", lv[:m.start()]):
                        out["run_time_min"] = fnum(m.group(1))
        elif ll.startswith(("injection port temperature", "temperature of injection port", "temperature of sample vaporization")):
            out["injection_temp_c"] = first_temp(val)
        elif ll.startswith("carrier gas"):
            out["carrier_gas"] = val.rstrip(".").lower()
        elif ll.startswith("split ratio"):
            m = re.search(r"(\d+\s*:\s*\d+)", val)
            out["split_ratio"] = m.group(1).replace(" ", "") if m else ("splitless" if "splitless" in lv else None)

    if mobile:
        mt = " ".join(mobile).lower()
        n_m = len(out.get("_phase_letters") or []) or len(mobile)
        ml = re.search(r"(?:mobile phases?|buffer solutions?|eluents?)\s+([A-E](?:\s*,\s*[A-E])*(?:\s*,?\s*and\s+[A-E])?)\b", " ".join(mobile), re.I)
        if ml and n_m == 1:
            n_m = len(set(re.findall(r"\b([A-E])\b", ml.group(1))))
        out["n_mobile_phases"] = n_m
        if out["elution_mode"] is None:
            stepwise = re.search(r"step[- ]?wise|step gradient|linear(?: concentration)? gradient|concentration gradient", mt)
            out["elution_mode"] = "gradient" if (n_m > 1 or stepwise) else "isocratic"
        out["mobile_phase_ph"] = ph_value(mt) if ph_value(mt) is not None else ph_value(" ".join(stock_text).lower())
        out["organic_solvents"] = sorted(s for s in SOLVENTS if re.search(r"(?<![\w-])" + re.escape(s.lower()), mt + " " + table_text(text)))
    out.pop("_phase_letters", None)
    out["gradient_table"] = gradient_table(text)
    if out["gradient_table"]:
        out["elution_mode"] = "gradient"
        if out["n_mobile_phases"] is None:
            out["n_mobile_phases"] = (len(out["gradient_table"][0]) - 2) // 2
    if groups:
        if xref_text and re.search(r"\bexcept\b", xref_text, re.I):
            exc = re.search(r"\bexcept\b(.*)$", xref_text, re.I).group(1).lower()
            out["xref_except"] = exc.strip(" .")
            if "time span" in exc:
                groups.discard("time_span")
        if stated_time:
            groups.discard("time_span")
        out["xref_time_span"] = "time_span" in groups
        names = sorted(g for g in groups if g != "time_span")
        out["xref_fields"] = names or None
        out["block_type"] = "partial_xref" if stated else "full_xref"
        out["xref_target"] = _xref_target(xref_text or "")
        if not names and out["xref_time_span"]:
            out["xref_fields"] = None
    else:
        out["block_type"] = "full"
    return out


def _xref_target(s: str) -> str | None:
    """Normalise the target of 'Proceed as directed in the operating conditions in the Assay (3) under X.'"""
    s = re.sub(r"\s+", " ", s)
    s = re.sub(r"\s+except\b.*$", "", s, flags=re.I)
    m = re.search(r"operating conditions?\s*(\d+)?\s*(?:in|of|under)?\s*(?:the\s+)?(.*)$", s, re.I)
    if m:
        tgt = m.group(2)
        if not tgt.strip(" .") and m.group(1):
            return "operating conditions " + m.group(1)
    else:
        m = re.search(r"(?:Proceed|Perform)\s*(?:as directed\s*)?(?:in|under)?\s*(?:the\s+)?(.*)$", s, re.I)
        if not m:
            return None
        tgt = m.group(1)
    tgt = re.sub(r"<[\d.]+>", "", tgt)
    tgt = re.sub(r"\s+", " ", tgt).strip(" .")
    return tgt or None
