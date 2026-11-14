"""The page text that needs a translation: every display string of the page list (mzdesk/spec.py, served at /api/spec).

Which strings count is decided once, in web/src/i18n/spec-rules.json. The browser reads the same file (localizeSpec in
web/src/i18n/index.ts), so this list and what the page translates cannot drift apart. A word that means something else in one block
("Species" among the LIPID MAPS cards) is listed there as "Species|lipid": see "context" in the rules.

The titles of the blocks that workflow steps send back (mzdesk/annotation, mzdesk/workflows) are not in the page list; run_titles() reads them out of
the Python source, so a new title shows up in the list (and in the coverage test) as soon as it is written.

    python -m mzdesk.spec_i18n      writes web/src/i18n/source/spec.json and run.json (the lists the translators work from)
"""
import ast
import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
I18N_DIR = REPO / "web" / "src" / "i18n"
RULES_PATH = I18N_DIR / "spec-rules.json"
SOURCE_PATH = I18N_DIR / "source" / "spec.json"
RUN_PATH = I18N_DIR / "source" / "run.json"
LOCALES_DIR = I18N_DIR / "locales"
UI_PATH = I18N_DIR / "source" / "ui.json"
LANGUAGES = ("zh", "ja", "de", "fr", "es")


def _rules():
    return json.loads(RULES_PATH.read_text(encoding="utf-8"))


def _strings(value):
    """The strings in a text value: one string, or a list of them (help, columns)."""
    if isinstance(value, str):
        return [value]
    if isinstance(value, list):
        return [x for x in value if isinstance(x, str)]
    return []


def _keyed(strings, context):
    """The dictionary keys of strings inside a block that has contexts: 'Species' there is 'Species|lipid'."""
    return [f"{s}|{context[s]}" if s in context else s for s in strings]


def _visit(node, kind, rules, text, lookup, context=None):
    spec = rules["kinds"][kind]
    if not isinstance(node, dict):
        return
    if kind == "block":
        context = rules.get("context", {}).get(node.get("id"), context)       # the items of a block share its context
    context = context or {}
    for key in spec.get("text", []):
        text.extend(_keyed(_strings(node.get(key)), context))
    for key in spec.get("keyed", []):
        keyed = node.get(key)
        if isinstance(keyed, dict):
            for name, value in keyed.items():
                text.extend(_keyed([name], context))
                text.extend(_keyed(_strings(value), context))
    for key in spec.get("lookup", []):
        lookup.extend(_keyed(_strings(node.get(key)), context))
    for key, child in spec.get("children", {}).items():
        for item in node.get(key) or []:
            _visit(item, child, rules, text, lookup, context)


def _split(spec_dict):
    rules = _rules()
    text, lookup = [], []
    for key, kind in rules["root"].items():
        for item in spec_dict.get(key) or []:
            _visit(item, kind, rules, text, lookup)
    return text, lookup


def _clean(strings):
    """Unique and sorted; a string without a letter ("", "2.7", "%") has nothing to translate."""
    return sorted({s for s in strings if any(c.isalpha() for c in s)})


def displayed_strings(spec_dict):
    """The strings localizeSpec translates in place."""
    return _clean(_split(spec_dict)[0])


def lookup_strings(spec_dict):
    """The strings the page code also uses as keys (section titles, cards item labels): translated by the component that draws them."""
    return _clean(_split(spec_dict)[1])


def translatable_strings(spec_dict):
    """Every English string of the page list that needs a translation: sorted, unique."""
    text, lookup = _split(spec_dict)
    return _clean(text + lookup)


# ---------------------------------------------------------------- titles of the blocks a workflow run sends back
RUN_SOURCES = ("mzdesk/annotation", "mzdesk/workflows")
# the functions that build a block, and the position of the title among their positional arguments
TITLED = {"_table": 1, "_notes": 1, "_engine_blocks": 1, "table_block": 1, "notes_block": 1, "block": 2}


def _call_name(call):
    return call.func.id if isinstance(call.func, ast.Name) else call.func.attr if isinstance(call.func, ast.Attribute) else None


def _literal(node):
    return node.value if isinstance(node, ast.Constant) and isinstance(node.value, str) else None


def run_titles(root=REPO):
    """(titles, unexplained): the block titles written as string literals in the workflow code (also the label of a figure), sorted and unique,
    and the places where a title is built some other way (a variable that is not the `title` parameter of a helper, an f-string): there the
    list cannot be complete, so tests/test_i18n.py wants this to be empty."""
    found, unexplained = [], []
    for folder in RUN_SOURCES:
        for path in sorted((Path(root) / folder).glob("*.py")):
            tree = ast.parse(path.read_text(encoding="utf-8"))
            where = path.relative_to(root).as_posix()
            for node in ast.walk(tree):
                if isinstance(node, ast.Call) and _call_name(node) in TITLED:
                    at = TITLED[_call_name(node)]
                    arg = node.args[at] if len(node.args) > at else None
                    if arg is not None and not (isinstance(arg, ast.Constant) and arg.value is None):
                        if _literal(arg) is not None:
                            found.append(_literal(arg))
                        elif not (isinstance(arg, ast.Name) and arg.id == "title"):       # a helper passing its own parameter on is fine
                            unexplained.append(f"{where}:{node.lineno}  {ast.unparse(arg)}")
                    for kw in node.keywords:
                        if kw.arg == "label" and _call_name(node) == "block" and _literal(kw.value) is not None:
                            found.append(_literal(kw.value))
                elif isinstance(node, ast.Dict):           # {"type": ..., "title": ...}
                    keys = [_literal(k) for k in node.keys]
                    if "type" in keys and "title" in keys:
                        value = node.values[keys.index("title")]
                        if _literal(value) is not None:
                            found.append(_literal(value))
                        elif not (isinstance(value, ast.Name) and value.id == "title"):
                            unexplained.append(f"{where}:{node.lineno}  {ast.unparse(value)}")
    return _clean(found), unexplained


def _write(path, strings):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(strings, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return strings


def write_source(path=SOURCE_PATH):
    from mzdesk.spec import spec

    return _write(path, translatable_strings(spec()))


def write_run_source(path=RUN_PATH):
    return _write(path, run_titles()[0])


if __name__ == "__main__":
    print(f"{len(write_source())} strings -> {SOURCE_PATH}")
    print(f"{len(write_run_source())} strings -> {RUN_PATH}")
