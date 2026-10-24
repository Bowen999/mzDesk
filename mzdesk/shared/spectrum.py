"""共享峰表解析器。MS1、MS2 各页都从这里读谱，不各自解析。"""
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


def parse_spectrum(text):
    """Parse a peak list into a sorted list of ``(mz, intensity)`` tuples.

    Accepts one ``mz intensity`` pair per line (whitespace-, comma- or tab-separated;
    extra columns are ignored; a bare ``mz`` gets intensity 1.0) and/or space-separated
    ``mz:intensity`` tokens (several per line allowed; ``mz:`` = intensity 1.0). The format
    is detected per line, and lines that do not start with a number are skipped as headers
    (``Name: PG 36:2`` must not yield a peak at m/z 36; ``PEPMASS=…``, ``Num Peaks: 13``).
    Malformed lines/tokens are skipped, non-finite values, m/z <= 0 and intensity <= 0 are
    dropped, and duplicate m/z (rounded to 5 dp) are collapsed keeping the most intense.
    """
    peaks = []
    for raw in (text or "").splitlines():
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
        if not pairs:                                    # plain "mz [intensity] [extra…]" line
            mz = _finite(toks[0])
            inten = _finite(toks[1]) if len(toks) > 1 else 1.0
            if mz is not None and inten is not None:
                pairs.append((mz, inten))
        peaks.extend(pairs)
    agg = {}
    for mz, inten in peaks:
        if mz <= 0 or inten <= 0:
            continue
        key = round(mz, 5)
        if inten > agg.get(key, 0.0):
            agg[key] = inten
    return sorted(agg.items())


_MGF_BEGIN = re.compile(r"^\s*BEGIN\s+IONS\s*$", re.I)
_MGF_END = re.compile(r"^\s*END\s+IONS\s*$", re.I)
_MSP_NAME = re.compile(r"^\s*NAME\s*:", re.I)
_NAME = re.compile(r"^\s*(?:TITLE|NAME)\s*[=:]\s*(.+?)\s*$", re.I)
_PRECURSOR = re.compile(r"^\s*(?:PEPMASS|PRECURSORMZ|PRECURSOR_MZ|PRECURSOR\s+M/Z)\s*[=:]\s*([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)", re.I)


def first_record(text):
    """First spectrum of an MGF (BEGIN IONS … END IONS) or MSP (records start with 'Name:') text.

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
