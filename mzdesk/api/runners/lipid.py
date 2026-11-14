"""Lipid Short Name 页：一个脂质 shorthand 名 -> LIPID+ 的 L1 / L2 / L3、LIPID MAPS 的 species / molecular species、一张示意结构。
名字规范化是 mzdesk.lipid_name_map（直接拷进来的 LIPID+ 名字包，纯本地）；示意结构按模板拼，见 mzdesk/lipid/structure.py。
名字被拒绝（俗名、不认识的类……）是输入错误：400 加一句原因，不是 500。画不了结构的类只给名字，加一句说明。
"""
import logging
import re

from mzdesk.api.common import InputError

log = logging.getLogger("mzdesk.api")

# 被拒绝的原因里，这几种是输入的形式错误（不是"写了个俗名"），保留各自的原因句
_FORM_ERRORS = {"empty", "not_string", "too_long", "id_or_number"}
COMMON_NAME_TEXT = "Common names are not supported. Use a shorthand such as PC 16:0/18:1."


def _present(**kv):
    return {k: v for k, v in kv.items() if v}


def _blank_code(result, key):
    """LIPID MAPS 名写不出来时的原因码（meta.lm 里）。"""
    return (result["meta"].get("lm") or {}).get(key)


def _maps_value(result, field, key):
    from mzdesk.lipid import texts

    text = result[field]
    if text:
        return text
    code = _blank_code(result, key)
    return {"value": "—", "help": texts.LM_REASON_TEXT.get(code, code)} if code else None


# branch_dropped 只有一个码，两种写法的意思不同：位置甲基 (15Me) / -15me 的碳不计入 L1；iso- / anteiso- 的碳数已含支链碳。按原文分开说
_POS_METHYL = re.compile(r"\d{1,3}-?[Mm]e(?![A-Za-z])")
_ISO = re.compile(r"(?<![A-Za-z])(?:anteiso|iso|ai)-?(?![a-z])")


def _flag_codes(result, raw):
    """flag 码列表；branch_dropped 按原文换成 branch_dropped（位置甲基）和 / 或 branch_type_dropped（iso / anteiso）。"""
    codes = []
    for f in result["flags"]:
        if f != "branch_dropped":
            codes.append(f)
            continue
        methyl, iso = _POS_METHYL.search(raw), _ISO.search(raw)
        codes += (["branch_dropped"] if methyl or not iso else []) + (["branch_type_dropped"] if iso else [])
    return codes


def flag_lines(result, raw=""):
    """Flags 只给英文短句（不显示内部码）；卡片 "—" 的 "?" 已经说过的（同一个码）不再重复。"""
    from mzdesk.lipid import texts

    said = {c for c in (_blank_code(result, "species") if not result["LM_species"] else None,
                        _blank_code(result, "molecular") if not result["LM_molecular"] else None) if c}
    said_text = {texts.LM_REASON_TEXT.get(c, c) for c in said}
    lines = [(f, texts.FLAG_TEXT.get(f, f)) for f in _flag_codes(result, raw)]
    return list(dict.fromkeys(t for f, t in lines if f not in said and t not in said_text))


def _reject_text(result, raw):
    """被拒绝时的原因句。没有 "数字:数字" 的链组成时，几乎都是俗名 / 类缩写 / 全名：统一给一句话。"""
    from mzdesk.lipid import texts

    if result["status"] == "rejected" and result["reason"] not in _FORM_ERRORS and not re.search(r"\d+:\d+", raw):
        return COMMON_NAME_TEXT
    return texts.REASON_TEXT.get(result["reason"], result["reason"])


def illustrate(result, raw):
    """(Illustration or None, svg or None, assumptions, note). When nothing can be drawn, the note is one sentence and the rest is empty.
    The Lipid Fragment Explanation workflow (mzdesk/workflows/lipid_fragments.py) explains the spectrum with this same structure."""
    from mzdesk.lipid import draw, structure as ls

    try:
        ill = ls.illustrate(result, raw)
    except ls.NoStructureForClass:
        return None, None, [], "No illustrative structure for this lipid class."
    except ls.BranchedChain:
        return None, None, [], "Branched chain: structure not drawn."
    except ls.Unsupported:
        return None, None, [], "No illustrative structure for this name."
    except Exception:                       # noqa: BLE001 - the names are still the answer; a drawing bug must not hide them
        log.exception("illustrative structure failed for %r", raw)
        return None, None, [], "No illustrative structure for this name."
    assumptions = list(ill.assumptions)
    if ls.positions_in_name(raw):
        assumptions.insert(0, "Positions given in the name are not used")
    return ill, draw.illustration_svg(ill), assumptions, None


def read_name(raw):
    """Lipid shorthand -> the normalize() result. A rejected or unsupported name is an input error (400) with one sentence saying why."""
    from mzdesk.lipid_name_map import normalize

    raw = str(raw or "").strip()
    if not raw:
        raise InputError("Enter a lipid name.")
    r = normalize(raw)
    if r["status"] != "ok":
        why = _reject_text(r, raw)
        shown = raw if len(raw) <= 60 else raw[:57] + "..."
        raise InputError(f"Cannot read “{shown}”. {why}" if r["status"] == "rejected" else f"“{shown}” is not supported. {why}")
    return r


def name_cards(r):
    """The LIPID+ and LIPID MAPS card data of a normalize() result: {"lipid_plus": {...}, "lipid_maps": {...}}."""
    return {
        "lipid_plus": {"values": _present(L1=r["L1"], L2=r["L2"], L3=r["L3"])},
        "lipid_maps": {"values": {k: v for k, v in (("Species", _maps_value(r, "LM_species", "species")),
                                                    ("Molecular species", _maps_value(r, "LM_molecular", "molecular"))) if v}},
    }


def run_lipid_short_name(inputs):
    raw = str(inputs.get("name") or "").strip()
    r = read_name(raw)
    _, svg, assumptions, no_structure = illustrate(r, raw)
    outputs = {
        **name_cards(r),
        "no_structure": {"items": [no_structure] if no_structure else []},
        "assumptions": {"items": assumptions},
        "flags": {"items": flag_lines(r, raw)},
    }
    if svg:
        outputs["structure"] = {"svg": svg}
    return {"outputs": outputs, "warnings": []}


RUNNERS = {"lipid_short_name": run_lipid_short_name}
