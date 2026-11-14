"""共享峰表解析器。MS1、MS2 各页都从这里读谱，不各自解析。"""
import json
import math
import os
import re
import urllib.error
import urllib.parse
import urllib.request

_SEP = re.compile(r"[\s,]+")
# '184.07 : 100' -> '184.07:100' (only between two numbers; never merges '184.07: 200.1:5')
_COLON_GAP = re.compile(r"(?<=[\d.])\s*:\s*(?=[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?(?:[\s,]|$))")

# Universal Spectrum Identifier resolver. GNPS2 covers GNPS / MassIVE USIs; a local mirror can
# be selected with MZDESK_USI_ENDPOINT (the USI is appended URL-encoded).
USI_ENDPOINT = "https://metabolomics-usi.gnps2.org/json/?usi1="
USI_TIMEOUT_S = 25
USI_MAX_BYTES = 20 << 20


class USIError(ValueError):
    """The USI could not be resolved (bad syntax, network error, no spectrum, malformed response)."""


def fetch_usi(usi):
    """Resolve a USI -> ([(mz, intensity), ...], precursor_mz or None).

    Only the public GNPS/MassIVE resolver is used. The peaks are cleaned with the same rules
    as parse_spectrum (positive finite m/z and intensity, duplicate m/z collapsed).
    """
    usi = (usi or "").strip()
    parts = usi.split(":")
    if not usi.startswith("mzspec:") or len(parts) < 5:
        raise USIError("A USI must look like mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB... "
                       "(mzspec:collection:accession:spectrum).")
    endpoint = os.environ.get("MZDESK_USI_ENDPOINT") or USI_ENDPOINT
    url = endpoint + urllib.parse.quote(usi, safe="")
    request = urllib.request.Request(url, headers={"User-Agent": "mzDesk USI resolver"})
    try:
        with urllib.request.urlopen(request, timeout=USI_TIMEOUT_S) as response:
            data = response.read(USI_MAX_BYTES + 1)
    except urllib.error.HTTPError as exc:
        raise USIError(f"The USI service answered HTTP {exc.code} for {usi}. "
                       "Check the USI, or paste the peak list instead.") from None
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        raise USIError(f"Could not reach the USI service ({exc}). Check the network, or paste the peak list instead.") from None
    if len(data) > USI_MAX_BYTES:
        raise USIError("The USI response is too large (over 20 MB).")
    try:
        obj = json.loads(data.decode("utf-8", "replace"))
    except (ValueError, UnicodeDecodeError):
        raise USIError("The USI service returned a response that is not JSON.") from None
    pairs = obj.get("peaks") if isinstance(obj, dict) else None
    if not isinstance(pairs, list) or not pairs:
        raise USIError(f"The USI service returned no peaks for {usi}.")
    agg = {}
    for pair in pairs:
        if not isinstance(pair, (list, tuple)) or len(pair) < 2:
            continue
        mz, inten = _finite(pair[0]), _finite(pair[1])
        if mz is None or inten is None or mz <= 0 or inten <= 0:
            continue
        key = round(mz, 5)
        if inten > agg.get(key, 0.0):
            agg[key] = inten
    if not agg:
        raise USIError(f"The USI service returned no usable peaks for {usi}.")
    precursor = _finite(obj.get("precursor_mz")) if isinstance(obj, dict) else None
    return sorted(agg.items()), (precursor if precursor and precursor > 0 else None)



def _finite(tok):
    try:
        v = float(tok)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


_ARR = r"""["']?\b{}\w*["']?\s*[=:]\s*(?:(?:np\.)?array\()?\s*\[([^\]]*)\]"""
_MZ_ARR = re.compile(_ARR.format(r"(?:m/z|mz)"), re.I)
_INT_ARR = re.compile(_ARR.format(r"(?:int|abundance)"), re.I)
_XML_PEAK = re.compile(r"<ms-ms-peak\b[^>]*>(.*?)</ms-ms-peak>", re.I | re.S)
_XML_MZ = re.compile(r"<mass-charge>\s*([^<\s]+)\s*</mass-charge>", re.I)
_XML_INT = re.compile(r"<intensity>\s*([^<\s]+)\s*</intensity>", re.I)
_PK_PEAK = re.compile(r"^\s*PK\$PEAK\s*:", re.I)
_HAS_PK = re.compile(r"^\s*PK\$PEAK\s*:", re.I | re.M)
_HAS_MGF = re.compile(r"^\s*BEGIN\s+IONS\s*$", re.I | re.M)
_HAS_MSP = re.compile(r"^\s*(?:NAME|NUM\s*PEAKS)\s*:", re.I | re.M)
_HAS_SIRIUS = re.compile(r"^\s*>(?:compound|parentmass|ms2)", re.I | re.M)
_SIRIUS_MS2 = re.compile(r"^\s*>ms2peaks\b", re.I | re.M)
_SIRIUS_NEXT = re.compile(r"^\s*>", re.M)
_CSV_HEAD = re.compile(r"\A\s*\"?m/?z\"?\s*[,;\t]", re.I)
_MZ_KEYS, _INT_KEYS = ("mz", "m/z"), ("intensity", "int", "abundance")


def _json_pairs(s):
    """(label, pairs) for a JSON list of [mz, int] pairs, {mz, intensity} objects, or a record
    object such as {"peaks": [[mz, int], ...], "precursor_mz": ...} (outer [] optional)."""
    if s[0] not in "[{":
        return None
    for cand in (s, "[" + s + "]"):
        try:
            obj = json.loads(cand)
        except ValueError:
            continue
        record = isinstance(obj, dict) and isinstance(obj.get("peaks"), list)
        if record:
            obj = obj["peaks"]
        if not isinstance(obj, list):
            continue
        pairs, label = [], "JSON record" if record else None
        for e in obj:
            if isinstance(e, (list, tuple)) and len(e) >= 2:
                mz, inten = _finite(e[0]), _finite(e[1])
                label = label or "GNPS peaks_json"
            elif isinstance(e, dict):
                low = {str(k).lower(): v for k, v in e.items()}
                mz = _finite(next((low[k] for k in _MZ_KEYS if k in low), None))
                inten = _finite(next((low[k] for k in _INT_KEYS if k in low), None))
                label = label or "JSON objects"
            else:
                continue
            if mz is not None and inten is not None:
                pairs.append((mz, inten))
        if pairs:
            return label, pairs
    return None


def _array_pairs(s):
    """Pairs from ``mz = [..]`` + ``intensities = [..]`` (matchms / Python style or a JSON object)."""
    m, n = _MZ_ARR.search(s), _INT_ARR.search(s)
    if not (m and n):
        return None
    mz = [_finite(t) for t in _SEP.split(m.group(1).strip()) if t]
    inten = [_finite(t) for t in _SEP.split(n.group(1).strip()) if t]
    if not mz or len(mz) != len(inten):
        return None
    return [(a, b) for a, b in zip(mz, inten) if a is not None and b is not None] or None


def _xml_pairs(s):
    out = []
    for body in _XML_PEAK.findall(s):
        m, i = _XML_MZ.search(body), _XML_INT.search(body)
        if m and i and _finite(m.group(1)) is not None and _finite(i.group(1)) is not None:
            out.append((_finite(m.group(1)), _finite(i.group(1))))
    return out


def _massbank_pairs(s):
    """Only the lines under ``PK$PEAK:`` (m/z, intensity, rel.int.); header numbers never count."""
    out, inside = [], False
    for ln in s.splitlines():
        if _PK_PEAK.match(ln):
            inside = True
        elif inside and re.match(r"\s*[\d.+-]", ln):
            toks = ln.split()
            mz, inten = _finite(toks[0]), _finite(toks[1]) if len(toks) > 1 else None
            if mz is not None and inten is not None:
                out.append((mz, inten))
        elif inside:
            inside = False
    return out


_NUM_START = re.compile(r"\s*[\d.+-]")
_PLAIN = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")      # 138.0662  1e5
_COMMA_NUM = re.compile(r"[-+]?\d+,\d+")                                # 138,0662 (decimal comma), 1,200 (thousands), or "mz,int" when nothing else separates the columns
_THOUSANDS = re.compile(r"[-+]?\d{1,3}(?:,\d{3})+(?:\.\d+)?")           # 1,200  1,234,567  1,234.5
_GROUP3 = re.compile(r"[-+]?\d{1,3},\d{3}(?:\.\d+)?")                  # 1,200  1,234.5: one thousands group; could also be 1.2 (three decimals) or the EI csv row 41,100.5
_ZERO_FRAC = re.compile(r"[-+]?\d+,0\d+")                               # 195,0877: nobody writes an intensity 0877, so this comma is a decimal one
NOTE_DECIMAL_COMMA = "Decimal commas (as in 138,0662) were read as decimal points."
NOTE_THOUSANDS = "Commas in numbers like 1,200 were read as thousands separators (1200)."
NOTE_EXTRA_NUMBERS = "Some lines have more than two numbers; only the first two (m/z, intensity) were used."


def _note(notes, message):
    if notes is not None and message not in notes:
        notes.append(message)


def _certain_thousands(cell):
    """1,234,567 can only be thousands groups; 1,200 could also be 1.2 written with three decimals, and 41,100.5 an EI csv row."""
    return bool(_THOUSANDS.fullmatch(cell)) and cell.count(",") > 1


def _is_plain(cell):
    """A number whose commas (if any) cannot be decimal commas."""
    return bool(_PLAIN.fullmatch(cell)) or _certain_thousands(cell)


def _cells(text):
    """The tab / space separated cells of a line, without the cells that are only commas ("195.0877 , 100")."""
    return [c for c in text.split() if c.strip(",")]


def _starts_number(cells):
    return bool(cells) and bool(_NUM_START.match(cells[0].lstrip(",")))


def _proves_decimal_comma(cells):
    """True when this row shows that the text writes decimals with a comma: 195,0877 (a 0 right after the comma), or a comma number
    standing among plain ones (138,0662<TAB>100; 100 32,5) that is not a bare thousands group like 1,200."""
    commas = [c for c in cells if _COMMA_NUM.fullmatch(c)]
    if not commas:
        return False
    if any(_ZERO_FRAC.fullmatch(c) and not _GROUP3.fullmatch(c) for c in commas):
        return True
    if not any(_is_plain(c) for c in cells):
        return False
    return bool(_COMMA_NUM.fullmatch(cells[0])) or any(not _GROUP3.fullmatch(c) for c in commas)


def _row_numbers(cells, decimal, notes):
    """The cells of one peak row (split by tab / space / ';') -> number strings, once every comma has been given its meaning.

    138,0662 is a decimal comma when another cell stands next to it (138,0662<TAB>100) or the text proves decimal commas elsewhere
    (``decimal``); 1,200 after a dot-decimal m/z is a thousands group. With nothing else separating the columns a comma still
    separates them: "41,100" is m/z 41, intensity 100 (the usual EI csv)."""
    plain = any(_is_plain(c) for c in cells)
    comma_mz = bool(_COMMA_NUM.fullmatch(cells[0]))
    out = []
    for i, cell in enumerate(cells):
        if "," not in cell:
            out.append(cell)
        elif _certain_thousands(cell) or (i > 0 and plain and not (decimal or comma_mz) and _GROUP3.fullmatch(cell)):
            out.append(cell.replace(",", ""))
            _note(notes, NOTE_THOUSANDS)
        elif not _COMMA_NUM.fullmatch(cell):
            out += [t for t in cell.split(",") if t]                    # "195.0877,100": a csv pair
        elif plain or decimal:
            out.append(cell.replace(",", "."))
            _note(notes, NOTE_DECIMAL_COMMA)
        else:
            out += cell.split(",")
    return out


def _line_pairs(s, notes=None):
    """(pairs, any ``mz:int`` token) for line-based text; headers (lines not starting with a number) are skipped.

    ``notes`` (a list, optional) receives what had to be guessed: decimal / thousands commas, rows with more than two numbers."""
    lines = s.splitlines()
    n_numeric = sum(1 for raw in lines if _NUM_START.match(raw))
    rows = []                                            # ("colon", "mz:int ..." line) | ("cells", cells of one peak) | ("cols", cells of a ';' column row)
    for raw in lines:                                    # ';' separates peaks only on lines that start with a number
        parts = [q.strip() for q in raw.split(";")] if _NUM_START.match(raw) else [raw]
        while parts and not parts[-1]:                   # "195.0877;100;" (a trailing ';')
            parts.pop()
        # "195.0877;100" and "195.0877;100;[M+H]+" are csv rows with ';' columns (m/z; intensity; …).
        # One line of numbers only ("100.1;200.2;300.3") is a ';' separated list of m/z.
        numeric = [bool(_PLAIN.fullmatch(p) or _COMMA_NUM.fullmatch(p) or _THOUSANDS.fullmatch(p)) for p in parts]
        if len(parts) > 1 and numeric[0] and numeric[1] and (not any(numeric[2:]) or (all(numeric) and n_numeric > 1)):
            rows.append(("cols", parts))
            continue
        for part in parts:
            line = _COLON_GAP.sub(":", part.strip())
            rows.append(("colon", line) if ":" in line else ("cells", _cells(line)))
    decimal = any(_proves_decimal_comma(r[1]) for r in rows if r[0] != "colon" and _starts_number(r[1]))
    peaks, colon = [], False
    for kind, data in rows:
        if kind == "colon":
            toks = [t for t in _SEP.split(data) if t]
            if not toks or not re.match(r"[\d.+-]", toks[0]):
                continue
            pairs = []
            for tok in toks:
                if ":" not in tok:
                    continue
                mz_s, _, int_s = tok.partition(":")
                mz = _finite(mz_s)
                inten = 1.0 if int_s == "" else _finite(int_s)
                if mz is not None and inten is not None:
                    pairs.append((mz, inten))
            colon = colon or bool(pairs)
            if not pairs:                                # plain "mz [intensity] [extra…]" line
                mz = _finite(toks[0])
                inten = _finite(toks[1]) if len(toks) > 1 else 1.0
                if mz is not None and inten is not None:
                    pairs.append((mz, inten))
            peaks.extend(pairs)
            continue
        # "41,100;55,85" is two "mz,int" peaks, not m/z 41.1 with intensity 55.85
        split = kind == "cols" and not decimal and not any(_is_plain(c) for c in data)
        for cells in [_cells(c) for c in data] if split else [data]:
            if not _starts_number(cells):
                continue
            nums = _row_numbers(cells, decimal, notes)
            mz = _finite(nums[0])
            inten = _finite(nums[1]) if len(nums) > 1 else 1.0
            if mz is not None and inten is not None:
                peaks.append((mz, inten))
                if sum(_finite(t) is not None for t in nums) > 2:
                    _note(notes, NOTE_EXTRA_NUMBERS)
    return peaks, colon


def _extract(text, notes=None):
    """(format label, raw ``(mz, intensity)`` pairs) of a peak list in any supported format; ``notes`` as in _line_pairs."""
    s = (text or "").strip()
    if not s:
        return "Unknown", []
    got = _json_pairs(s)
    if got:
        return got
    pairs = _array_pairs(s)
    if pairs:
        return "Two arrays (matchms)", pairs
    pairs = _xml_pairs(s)
    if pairs:
        return "HMDB XML", pairs
    if _HAS_PK.search(s):
        return "MassBank record", _massbank_pairs(s)
    # SIRIUS .ms can hold >ms1peaks and >ms2peaks: only the MS2 section belongs to an MS2 input.
    sirius = _SIRIUS_MS2.search(s)
    if sirius:
        body = s[sirius.end():]
        nxt = _SIRIUS_NEXT.search(body)
        if nxt:
            body = body[:nxt.start()]
        pairs, _ = _line_pairs(body, notes)
        if pairs:
            return "SIRIUS .ms", pairs
    pairs, colon = _line_pairs(s, notes)
    if not pairs:
        return "Unknown", []
    if _HAS_MGF.search(s):
        label = "MGF"
    elif _HAS_SIRIUS.search(s):
        label = "SIRIUS .ms"
    elif _HAS_MSP.search(s):
        label = "MSP"
    elif _CSV_HEAD.match(s):
        label = "CSV"
    else:
        label = "MS-DIAL / MoNA (mz:int)" if colon else "Peak list (one per line)"
    return label, pairs


def detect_format(text):
    """Human label of the input format (``"Unknown"`` when no peaks are found)."""
    return _extract(text)[0]


def parse_spectrum_notes(text):
    """``(peaks, notes)``: parse_spectrum plus the list of guesses made while reading the text (decimal or thousands commas, rows with more
    than two numbers). Shown to the user as page warnings by api.common.spectrum."""
    notes = []
    agg = {}
    for mz, inten in _extract(text, notes)[1]:
        if mz <= 0 or inten <= 0:
            continue
        key = round(mz, 5)
        if inten > agg.get(key, 0.0):
            agg[key] = inten
    return sorted(agg.items()), notes


def parse_spectrum(text):
    """Parse a peak list into a sorted list of ``(mz, intensity)`` tuples.

    Auto-detects: one ``mz intensity`` pair per line (whitespace-, comma- or tab-separated; extra
    columns ignored; a bare ``mz`` gets intensity 1.0, so HMDB peak lists, CSV and SIRIUS .ms
    work), ``mz;intensity`` csv rows, ``mz:intensity`` tokens (MS-DIAL / MoNA; several per line, ``mz:`` = 1.0), MSP, MGF,
    GNPS ``[[mz, int], …]``, JSON objects, ``mz = [...]`` / ``intensities = [...]`` arrays,
    HMDB ``<ms-ms-peak>`` XML and MassBank records (only the ``PK$PEAK`` lines). In line-based
    text, lines that do not start with a number are headers (``Name: PG 36:2`` must not yield a
    peak at m/z 36; ``PEPMASS=…``, ``Num Peaks: 13``). A comma is a column separator when nothing else separates the columns
    (``195.0877,100``); next to a tab, space or ``;`` it is a decimal comma (``138,0662<TAB>100``) or, after a dot-decimal m/z, a
    thousands separator (``100.0 1,200``). Malformed entries are skipped,
    non-finite values, m/z <= 0 and intensity <= 0 are dropped, and duplicate m/z (rounded to
    5 dp) are collapsed keeping the most intense.
    """
    return parse_spectrum_notes(text)[0]


_MGF_BEGIN = re.compile(r"^\s*BEGIN\s+IONS\s*$", re.I)
_MGF_END = re.compile(r"^\s*END\s+IONS\s*$", re.I)
_MSP_NAME = re.compile(r"^\s*(?:NAME\s*:|>compound\b)", re.I)       # an MSP record starts at NAME:, a SIRIUS one at >compound
_MB_END = re.compile(r"^\s*//\s*$")                                 # a MassBank record ends at //
_NAME = re.compile(r"^\s*(?:(?:TITLE|NAME|CH\$NAME)\s*[=:]|>compound\s)\s*(.+?)\s*$", re.I)
_PRECURSOR = re.compile(r"^\s*(?:(?:PEPMASS|PRECURSORMZ|PRECURSOR_MZ|PRECURSOR\s+M/Z)\s*[=:]|MS\$FOCUSED_ION:\s*PRECURSOR_M/Z|>parentmass)\s*"
                        r"([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)", re.I)


def first_record(text):
    """First spectrum of an MGF (BEGIN IONS … END IONS), MSP ('Name:'), SIRIUS ('>compound') or MassBank (ends with '//') text.

    The precursor is also read from MassBank ``MS$FOCUSED_ION: PRECURSOR_M/Z`` and SIRIUS ``>parentmass``.

    Returns {"peaks_text", "precursor_mz" (float or None), "name" (str or None), "n_records"}.
    Text without record markers is one record. ``peaks_text`` still holds the header lines;
    parse_spectrum skips them.
    """
    lines = (text or "").splitlines()
    records = []
    if any(_MGF_BEGIN.match(ln) for ln in lines):
        cur = None
        for ln in lines:
            if _MGF_BEGIN.match(ln):
                cur = []
            elif _MGF_END.match(ln):
                if cur is not None:
                    records.append(cur)
                cur = None
            elif cur is not None:
                cur.append(ln)
        if cur:                                  # unterminated last record
            records.append(cur)
    else:
        cur = []
        for ln in lines:
            # Only a block that already holds a peak line is a record; a leading comment /
            # license block before "Name:" must not become the first (empty) record.
            if _MSP_NAME.match(ln) and any(re.match(r"\s*[\d.+-]", x) for x in cur):
                records.append(cur)
                cur = []
            cur.append(ln)
            if _MB_END.match(ln):
                records.append(cur)
                cur = []
        if any(x.strip() for x in cur):
            records.append(cur)
    rec = records[0] if records else []
    name = precursor = None
    for ln in rec:
        m = _NAME.match(ln)
        if m and name is None:
            name = m.group(1)
        m = _PRECURSOR.match(ln)
        if m and precursor is None:
            value = _finite(m.group(1))
            precursor = value if value and value > 0 else None
    if (precursor is None or name is None) and (text or "").strip().startswith("{"):
        try:
            obj = json.loads((text or "").strip())
        except ValueError:
            obj = None
        if isinstance(obj, dict):
            value = _finite(obj.get("precursor_mz"))
            precursor = precursor if precursor is not None else (value if value and value > 0 else None)
            if name is None and isinstance(obj.get("name"), str):
                name = obj["name"].strip() or None
    return {"peaks_text": "\n".join(rec), "precursor_mz": precursor, "name": name, "n_records": max(len(records), 1)}
