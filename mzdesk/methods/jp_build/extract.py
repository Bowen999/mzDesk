"""PDF -> cleaned text (JP18 official English PDFs).

Rebuilds visual lines from PyMuPDF word fragments (justified text is split into many pieces), maps the
special-font glyphs back to real characters (% degree prime Greek letters ...), marks bold runs as **...**,
monograph / chapter titles (size >= 12) as '### ', supplement instruction lines (font Go075Pif, e.g.
"Change the Purity (3) as follows:") as '@@I ', and drops running headers / footers by position.

Output text format: one visual line per text line; '@@PAGE n' marks the start of PDF page n (1-based).
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

# (font, char) -> real char. Gj202Pif is the "%" / degree font, Oh107Pi carries primes.
GLYPH = {("Gj202Pif", "z"): "%", ("Gj202Pif", "9"): "°", ("Oh107Pi", "?"): "′", ("Oh107Pi", "!"): "′",
         ("Gj200Pif", "s"): "[", ("Gj200Pif", "t"): "]"}
# Symbol-encoded Greek (Go068Pif upright, Go069Pif italic, Oh102Pi capitals)
GREEK = {"a": "α", "b": "β", "c": "χ", "d": "δ", "e": "ε", "f": "φ", "g": "γ", "h": "η", "i": "ι", "j": "φ",
         "k": "κ", "l": "λ", "m": "μ", "n": "ν", "p": "π", "q": "θ", "r": "ρ", "s": "σ", "t": "τ", "u": "υ",
         "v": "ω", "w": "ω", "x": "ξ", "y": "ψ", "z": "ζ", "D": "Δ", "P": "Π", "S": "Σ", "Q": "Ω", "F": "Φ"}
GREEK_FONTS = {"Go068Pif", "Go069Pif", "Oh102Pi"}
# ligature / punctuation glyphs that the PDF fonts encode as odd characters (e.g. "Identiˆcation" = "Identification")
LIGATURE = {"ˆ": "fi", "‰": "fl", "Š": "ff", "‹": "ffi", "Œ": "ffl", "‚": ".", "„": "–",
            "¿": "–"}
FULLWIDTH = {"＝": "=", "＋": "+", "－": "-", "＜": "<", "＞": ">", "≦": "≤", "≧": "≥", "･": "·", "・": "·",
             "\ue2ea": "%", "（": "(", "）": ")", "：": ":", "，": ",", "．": ".", "／": "/", "％": "%"}
BOLD = {"Go096Pif"}
TITLE_SIZE = 13.5   # bold text at least this large = monograph / chapter title
PART_SIZE = 25.0    # volume part markers ("Official Monographs") are set very large
INSTR_FONT = "Go075Pif"


def page_bounds(page) -> tuple[float, float]:
    """(header_y, footer_y): body text lies between them. Main volumes are 603x845 pt, supplements 595x842."""
    return (55.0, 800.0) if page.rect.height > 843 else (70.0, 775.0)


def span_text(s) -> str:
    out = []
    font = s["font"]
    for ch in s["text"]:
        if (font, ch) in GLYPH:
            out.append(GLYPH[(font, ch)])
        elif font in GREEK_FONTS and ch in GREEK:
            out.append(GREEK[ch])
        elif ch in LIGATURE and (font not in GREEK_FONTS or ch not in GREEK):
            out.append(LIGATURE[ch])
        else:
            out.append(FULLWIDTH.get(ch, ch))
    return "".join(out)


def visual_lines(page) -> list[str]:
    w = page.rect.width
    header_y, footer_y = page_bounds(page)
    frags = []
    for b in page.get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            x0, y0, x1, y1 = l["bbox"]
            if y0 < header_y or y0 > footer_y:
                continue
            spans = [s for s in l["spans"] if s["text"]]
            if not spans:
                continue
            col = 0 if x0 < w / 2 - 10 else 1
            if x1 - x0 > w * 0.55:
                col = -1  # full-width line
            base = max(s["origin"][1] for s in spans)
            frags.append((col, base, x0, x1, spans))
    frags.sort(key=lambda f: (f[0], f[1], f[2]))
    lines = []
    for col, base, x0, x1, spans in frags:
        if lines and lines[-1]["col"] == col and abs(lines[-1]["base"] - base) <= 3.0:
            lines[-1]["parts"].append((x0, x1, spans))
        else:
            lines.append({"col": col, "base": base, "parts": [(x0, x1, spans)]})
    out = []
    for ln in sorted(lines, key=lambda l: (l["col"] if l["col"] >= 0 else -1, l["base"])):
        pieces, size = [], 0
        prev_x1 = None
        instr_chars = other_chars = title_chars = big_chars = 0
        for x0, x1, spans in sorted(ln["parts"], key=lambda p: p[0]):
            if prev_x1 is not None:
                gap = x0 - prev_x1
                pieces.append(("  " if gap > 15 else (" " if gap > 0.8 else ""), None))
            prev_x1 = x1
            for s in spans:
                t = span_text(s)
                if t.strip():
                    size = max(size, s["size"])
                    if s["font"] in BOLD and s["size"] >= TITLE_SIZE:
                        title_chars += len(t.strip())
                    if s["size"] >= PART_SIZE:
                        big_chars += len(t.strip())
                    if s["font"] == INSTR_FONT:
                        instr_chars += len(t.strip())
                    else:
                        other_chars += len(t.strip())
                pieces.append((t, s["font"] in BOLD))
        txt, in_bold = "", False
        for t, bold in pieces:
            if bold is None:
                txt += t
                continue
            if bold and not in_bold:
                txt += (" " if txt and not txt.endswith(" ") else "") + "**"
            if not bold and in_bold:
                txt = txt.rstrip() + "** "
            txt += t
            in_bold = bold
        if in_bold:
            txt = txt.rstrip() + "**"
        txt = re.sub(r" {3,}", "  ", txt).strip()
        if not txt:
            continue
        n_chars = instr_chars + other_chars
        if title_chars >= 0.5 * n_chars or big_chars >= 0.5 * n_chars:
            txt = "### " + txt.replace("**", "")
        elif instr_chars and instr_chars >= 0.8 * (instr_chars + other_chars):
            txt = "@@I " + txt.replace("**", "")
        out.append(txt)
    return out


def printed_page(page) -> str:
    """Printed page number from the running header (digits only; '' when absent, e.g. front matter)."""
    header_y, _ = page_bounds(page)
    for b in page.get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            if l["bbox"][1] < header_y:
                t = "".join(s["text"] for s in l["spans"]).strip()
                if re.fullmatch(r"\d{1,5}", t):
                    return t
    return ""


def extract_pdf(pdf: str | Path, out: str | Path, *, first: int = 1, last: int | None = None) -> int:
    import fitz
    doc = fitz.open(str(pdf))
    res = []
    n = len(doc)
    for pno in range(first, (last or n) + 1):
        res.append(f"@@PAGE {pno} {printed_page(doc[pno - 1])}".rstrip())
        res.extend(visual_lines(doc[pno - 1]))
    text = "\n".join(res)
    text = re.sub(r"(\d)\s*°\s*C\b", r"\1 °C", text)
    Path(out).write_text(text, encoding="utf-8")
    return n


if __name__ == "__main__":
    extract_pdf(sys.argv[1], sys.argv[2])
