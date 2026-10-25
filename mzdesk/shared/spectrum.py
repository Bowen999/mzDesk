"""共享峰表解析器。MS1、MS2 各页都从这里读谱，不各自解析。"""
import json
import math
import re

_SEP = re.compile(r"[\s,]+")
# '184.07 : 100' -> '184.07:100' (only between two numbers; never merges '184.07: 200.1:5')
_COLON_GAP = re.compile(r"(?<=[\d.])\s*:\s*(?=[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?(?:[\s,]|$))")


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
_CSV_HEAD = re.compile(r"\A\s*\"?m/?z\"?\s*[,;\t]", re.I)
_MZ_KEYS, _INT_KEYS = ("mz", "m/z"), ("intensity", "int", "abundance")


def _json_pairs(s):
    """(label, pairs) for a JSON list of [mz, int] pairs or {mz, intensity} objects (outer [] optional)."""
    if s[0] not in "[{":
        return None
    for cand in (s, "[" + s + "]"):
        try:
            obj = json.loads(cand)
        except ValueError:
            continue
        if not isinstance(obj, list):
            continue
        pairs, label = [], None
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


def _line_pairs(s):
    """(pairs, any ``mz:int`` token) for line-based text; headers (lines not starting with a number) are skipped."""
    peaks, colon = [], False
    segments = []
    for raw in s.splitlines():                           # ';' separates peaks only on lines that start with a number
        segments += raw.split(";") if re.match(r"\s*[\d.+-]", raw) else [raw]
    for raw in segments:
        line = _COLON_GAP.sub(":", raw.strip())
        toks = [t for t in _SEP.split(line) if t]
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
        if not pairs:                                    # plain "mz [intensity] [extra…]" line
            mz = _finite(toks[0])
            inten = _finite(toks[1]) if len(toks) > 1 else 1.0
            if mz is not None and inten is not None:
                pairs.append((mz, inten))
        peaks.extend(pairs)
    return peaks, colon


def _extract(text):
    """(format label, raw ``(mz, intensity)`` pairs) of a peak list in any supported format."""
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
    pairs, colon = _line_pairs(s)
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


def parse_spectrum(text):
    """Parse a peak list into a sorted list of ``(mz, intensity)`` tuples.

    Auto-detects: one ``mz intensity`` pair per line (whitespace-, comma- or tab-separated; extra
    columns ignored; a bare ``mz`` gets intensity 1.0, so HMDB peak lists, CSV and SIRIUS .ms
    work), ``mz:intensity`` tokens (MS-DIAL / MoNA; several per line, ``mz:`` = 1.0), MSP, MGF,
    GNPS ``[[mz, int], …]``, JSON objects, ``mz = [...]`` / ``intensities = [...]`` arrays,
    HMDB ``<ms-ms-peak>`` XML and MassBank records (only the ``PK$PEAK`` lines). In line-based
    text, lines that do not start with a number are headers (``Name: PG 36:2`` must not yield a
    peak at m/z 36; ``PEPMASS=…``, ``Num Peaks: 13``). Malformed entries are skipped,
    non-finite values, m/z <= 0 and intensity <= 0 are dropped, and duplicate m/z (rounded to
    5 dp) are collapsed keeping the most intense.
    """
    agg = {}
    for mz, inten in _extract(text)[1]:
        if mz <= 0 or inten <= 0:
            continue
        key = round(mz, 5)
        if inten > agg.get(key, 0.0):
            agg[key] = inten
    return sorted(agg.items())


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
            if _MSP_NAME.match(ln) and any(x.strip() for x in cur):
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
    return {"peaks_text": "\n".join(rec), "precursor_mz": precursor, "name": name, "n_records": max(len(records), 1)}
