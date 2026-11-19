"""The language switch: the lists the translators work from, and the five dictionaries (web/src/i18n/locales/<lang>.json).

source/spec.json = the page text of mzdesk/spec.py that needs a translation (python -m mzdesk.spec_i18n writes it);
source/ui.json   = the interface strings, the first argument of t("...") in web/src (node web/scripts/extract-i18n.mjs writes it);
source/run.json  = the titles of the blocks a workflow run sends back, read out of mzdesk/annotation and mzdesk/workflows (python -m mzdesk.spec_i18n writes it).
A word with two meanings has one entry per meaning: t("System", { context: "theme" }) and the page text of a block named under "context" in
spec-rules.json are listed as "System|theme"; a language without such an entry falls back to the plain "System".
"""
import json
import re

import pytest

from mzdesk import spec as spec_module
from mzdesk import spec_i18n
from mzdesk.spec import spec

PLACEHOLDER = re.compile(r"\{(\w+)\}")
RELEASE = (spec_module.DEFERRED_SUBGROUPS, spec_module.DEFERRED_MODULES)     # as the module defines them, before conftest.py switches the held-back pages on


@pytest.fixture(autouse=True)
def release_pages(monkeypatch):
    """The lists are made from the pages the release offers; conftest.py switches DMPK, ICP-MS and MRM on for the other test files."""
    monkeypatch.setattr(spec_module, "DEFERRED_SUBGROUPS", RELEASE[0])
    monkeypatch.setattr(spec_module, "DEFERRED_MODULES", RELEASE[1])


def _json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def _locale(lang):
    return _json(spec_i18n.LOCALES_DIR / f"{lang}.json")


def _holes(text):
    return sorted(PLACEHOLDER.findall(text))


# ---------------------------------------------------------------- the dictionaries
@pytest.mark.parametrize("lang", spec_i18n.LANGUAGES)
def test_every_string_is_translated(lang):
    needed = set(spec_i18n.translatable_strings(spec())) | set(_json(spec_i18n.UI_PATH)) | set(_json(spec_i18n.RUN_PATH))
    words = _locale(lang)
    missing = sorted(k for k in needed if not str(words.get(k, "")).strip())
    assert not missing, f"{lang}.json lacks {len(missing)} of {len(needed)} strings, e.g. {missing[:5]}"


@pytest.mark.parametrize("lang", spec_i18n.LANGUAGES)
def test_a_translation_keeps_the_placeholders_of_its_english_text(lang):
    wrong = {k: v for k, v in _locale(lang).items() if _holes(k) != _holes(v)}
    assert not wrong, f"{lang}.json changes the {{placeholders}} of: {list(wrong)[:5]}"


@pytest.mark.parametrize("lang", spec_i18n.LANGUAGES)
def test_a_dictionary_is_a_flat_object_of_strings(lang):
    words = _locale(lang)
    assert isinstance(words, dict) and all(isinstance(k, str) and isinstance(v, str) for k, v in words.items())


# ---------------------------------------------------------------- the lists
def test_source_list_of_page_text_is_up_to_date():
    assert _json(spec_i18n.SOURCE_PATH) == spec_i18n.translatable_strings(spec()), "run: python -m mzdesk.spec_i18n"


def test_source_list_of_run_titles_is_up_to_date_and_complete():
    titles, unexplained = spec_i18n.run_titles()
    assert _json(spec_i18n.RUN_PATH) == titles, "run: python -m mzdesk.spec_i18n"
    assert unexplained == [], "a block title that is not a string literal cannot be listed: write it as one, or teach spec_i18n.run_titles"
    assert titles == sorted(set(titles)) and len(titles) > 20
    assert {"Ranked candidates", "Library hits", "Method sheet", "Marker ions", "Query spectrum"} <= set(titles)


def test_run_titles_are_read_from_the_calls_that_build_a_block(tmp_path):
    folder = tmp_path / "mzdesk" / "workflows"
    folder.mkdir(parents=True)
    (folder / "demo.py").write_text(
        'def table_block(block_id, title, columns, rows):\n    return block("table", block_id, title, {})\n'
        'def one(ctx, name):\n'
        '    a = table_block("t1", "First table", ["A"], [])\n'
        '    b = notes_block("n1", "A note title", [])\n'
        '    c = block("figure", "f1", None, {}, label="A figure")\n'
        '    d = table_block("t2", f"Hits for {name}", ["A"], [])\n'
        '    e = {"type": "notes", "id": "x", "title": "Dict title"}\n'
        '    f = table_block("t3", "no letters 123".replace("a", "b"), ["A"], [])\n', encoding="utf-8")
    (tmp_path / "mzdesk" / "annotation").mkdir()
    titles, unexplained = spec_i18n.run_titles(tmp_path)
    assert titles == ["A figure", "A note title", "Dict title", "First table"]
    assert [u.split("  ")[1] for u in unexplained] == ["f'Hits for {name}'", "'no letters 123'.replace('a', 'b')"]
    assert all(u.startswith("mzdesk/workflows/demo.py:") for u in unexplained)


def test_source_list_of_interface_text_is_sorted_and_unique():
    ui = _json(spec_i18n.UI_PATH)
    assert ui == sorted(set(ui)) and all(isinstance(s, str) and s.strip() for s in ui)


def test_page_text_list_is_sorted_unique_and_has_something_to_translate():
    strings = spec_i18n.translatable_strings(spec())
    assert strings == sorted(set(strings)) and len(strings) > 300
    assert all(any(c.isalpha() for c in s) for s in strings), "a string without a letter has nothing to translate"


def test_lookup_keys_are_listed_too_because_the_components_translate_them():
    s = spec()
    lookup = spec_i18n.lookup_strings(s)
    assert {"Advanced", "Optional verification", "Best hit"} <= set(lookup)
    assert set(lookup) <= set(spec_i18n.translatable_strings(s))


def test_what_is_a_value_or_an_id_is_not_page_text():
    """Engine names, ids, the adduct lists of a precursor field and example inputs name something or are sent to the backend: they are never in the list.
    (The options of a select / radio / multiselect are in it: they are drawn through t() and sent unchanged; see the next test.)"""
    strings = set(spec_i18n.translatable_strings(spec()))
    for value in ("RDKit", "ms2_clean", "esi_ms2_annotation", "[M+2H]2+", "[M+ACN+H]+", "[M+FA-H]-", "L-leucine", "Mus musculus"):
        assert value not in strings, value


def test_the_options_of_a_select_are_listed_as_lookup_text():
    """Options are also the values sent to the backend, so localizeSpec leaves them alone and the field that draws one calls t() on it."""
    s = spec()
    lookup = set(spec_i18n.lookup_strings(s))
    shown = set(spec_i18n.displayed_strings(s))
    for option in ("Positive", "Negative", "None", "Any", "Either / unknown", "Formic acid", "MassBank", "[M+H]+", "Orbitrap 60k", "DTIMS (stepped field)", "Water/acetonitrile + phosphate"):
        assert option in lookup | shown, option
    assert "" not in set(spec_i18n.translatable_strings(s)), "the empty choice of a select has nothing to translate"


def test_every_option_with_a_letter_is_in_the_list():
    s = spec()
    listed = set(spec_i18n.translatable_strings(s))
    options = {o for page in s["modules"] + s["pipelines"] for f in page.get("inputs", []) for o in f.get("options", []) if isinstance(o, str)}
    assert len(options) > 100
    assert sorted(o for o in options if any(c.isalpha() for c in o) and o not in listed) == []


def test_rules_are_read_the_same_way_as_in_the_browser():
    """tests/i18n.test.mjs runs localizeSpec on the real page list and compares with displayed_strings(); this pins the reading of the rules."""
    fake = {
        "groups": [{"id": "g", "title": "Modules"}], "subgroups": [],
        "pipelines": [{"id": "p", "title": "Flow", "one_liner": "A → B", "based_on": ["MetFrag"], "stages": [{"id": 1, "title": "Adduct"}],
                       "inputs": [{"type": "select", "id": "s", "label": "Mode", "options": ["", "Quick", "100"], "value": "Quick", "section": "Advanced", "help": ["Tip one", "Tip two"]},
                              {"type": "precursor", "id": "p", "ion_modes": ["Negative"], "adducts": {"Negative": ["[M-H]-"]}, "visible_when": {"field": "s", "equals": "Full"}}],
                       "examples": [{"label": "Caffeine", "values": {"s": "Quick"}}]}],
        "modules": [{"id": "m", "title": "Peaks", "one_liner": "x → y", "network": "Sent online.", "example": {"a": "Positive"},
                     "inputs": [{"type": "table", "id": "t", "label": "Alkanes", "columns": ["Carbon", "RT"], "caption": "Per line", "placeholder": "2.7", "auto_label": "Auto"}],
                     "outputs": [{"type": "table", "id": "o", "title": "Hits", "columns": ["Name", "Score"], "help": "Tip", "column_help": {"Score": "Higher is better"}},
                                 {"type": "cards", "id": "c", "title": "Totals", "items": [{"label": "Best", "unit": "ppm", "help": "Hint"}, {"label": "Count", "unit": "%"}]},
                                 {"type": "figure", "label": "Plot"}]}],
    }
    assert spec_i18n.displayed_strings(fake) == sorted([
        "Modules", "Flow", "A → B", "Mode", "Tip one", "Tip two", "Caffeine", "Peaks", "x → y", "Sent online.", "Alkanes", "Carbon", "RT", "Per line", "Auto",
        "Hits", "Name", "Score", "Tip", "Higher is better", "Totals", "ppm", "Hint", "Plot"])
    # lookup: section titles, cards item labels and the options of a select (without the empty choice and what has no letter); not adducts, ion modes or visible_when
    assert spec_i18n.lookup_strings(fake) == ["Advanced", "Best", "Count", "Quick"]


def test_a_context_is_read_per_block_and_listed_as_text_and_context(monkeypatch):
    """'context' in the rules: a word inside that block is listed as 'word|context' (columns, column help keys, cards labels); other blocks keep the plain word."""
    fake = {"groups": [], "subgroups": [], "pipelines": [], "modules": [{"id": "m", "title": "Page", "outputs": [
        {"type": "table", "id": "a", "title": "Hits", "columns": ["Column", "Value"], "column_help": {"Column": "Which one"}},
        {"type": "table", "id": "b", "title": "Methods", "columns": ["Column"]},
        {"type": "cards", "id": "c", "title": "Totals", "items": [{"label": "Species", "unit": "Chains"}]}]}]}
    rules = spec_i18n._rules()
    rules["context"] = {"a": {"Column": "metadata"}, "c": {"Species": "lipid"}}
    monkeypatch.setattr(spec_i18n, "_rules", lambda: rules)
    assert spec_i18n.displayed_strings(fake) == ["Chains", "Column", "Column|metadata", "Hits", "Methods", "Page", "Totals", "Value", "Which one"]
    assert spec_i18n.lookup_strings(fake) == ["Species|lipid"]


def test_a_context_names_one_block_of_one_page_and_a_word_the_page_shows():
    s = spec()
    contexts = _json(spec_i18n.RULES_PATH)["context"]
    blocks = [b.get("id") for page in s["modules"] + s["pipelines"] for b in page.get("outputs", [])]
    listed = set(spec_i18n.translatable_strings(s)) | set(_json(spec_i18n.UI_PATH))
    for block, words in contexts.items():
        assert blocks.count(block) == 1, f"{block}: a block id in 'context' must be the id of exactly one block"
        for word, context in words.items():
            assert f"{word}|{context}" in listed, f"{block}: '{word}' is not shown there, or is not text"
            assert word in listed, f"'{word}' has no plain entry for the languages without '{word}|{context}'"


@pytest.mark.parametrize("lang", spec_i18n.LANGUAGES)
def test_a_context_entry_has_its_plain_text_too(lang):
    """'System|theme' falls back to 'System' where it is missing, so the plain entry must exist."""
    words = _locale(lang)
    pairs = {k: k.split("|")[0] for k in words if "|" in k}
    assert pairs, "the contexts of the interface and of the page text"
    assert [k for k, plain in pairs.items() if plain not in words] == []


def test_every_context_has_a_different_translation_in_some_language():
    for key in (k for k in _locale("zh") if "|" in k):
        plain = key.split("|")[0]
        assert any(_locale(lang)[key] != _locale(lang)[plain] for lang in spec_i18n.LANGUAGES), f"{key}: the same as '{plain}' everywhere, so it needs no context"


def test_every_kind_the_rules_name_exists():
    rules = _json(spec_i18n.RULES_PATH)
    kinds = rules["kinds"]
    assert set(rules["root"].values()) <= set(kinds)
    assert all(set(k.get("children", {}).values()) <= set(kinds) for k in kinds.values())
