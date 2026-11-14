"""Shared test setup: the page examples are compound names; resolve them from spec.EXAMPLE_STRUCTURES so tests stay offline."""
import os

import pytest

from mzdesk import spec
from mzdesk.identity import resolve


@pytest.fixture(autouse=True)
def stats_in_tmp(tmp_path, monkeypatch):
    """Usage stats go to a throwaway file in tests, never to ~/.mzdesk."""
    from mzdesk.api import stats
    monkeypatch.setenv("MZDESK_STATS_DB", str(tmp_path / "stats.sqlite"))
    monkeypatch.setattr(stats, "_con", None)


@pytest.fixture(autouse=True)
def offline_example_names(monkeypatch):
    real = resolve.to_smiles
    monkeypatch.setattr(resolve, "to_smiles", lambda q, kind="Auto": spec.EXAMPLE_STRUCTURES.get((q or "").strip()) or real(q, kind))


@pytest.fixture(autouse=True)
def offline_pubchem(request, monkeypatch):
    """The exact-mass search answers from tests/fake_pubchem.py, never from PubChem (tests marked `network` use the real service)."""
    from mzdesk.ms1 import exact_mass
    exact_mass._get_json.cache_clear()
    if not request.node.get_closest_marker("network"):
        from fake_pubchem import FakePubChem
        monkeypatch.setattr(exact_mass, "_http", FakePubChem())
    yield
    exact_mass._get_json.cache_clear()


@pytest.fixture(autouse=True)
def held_back_pages_on(request, monkeypatch):
    """DMPK, ICP-MS and MRM Transitions wait for 2.0 and are not offered (spec.DEFERRED_*), but their engines stay tested: switch them on here.
    tests/test_release_pages.py runs without this, to check what the release really offers."""
    if request.module.__name__.endswith("test_release_pages"):
        yield
        return
    from mzdesk.api import app as app_module
    monkeypatch.setattr(spec, "DEFERRED_SUBGROUPS", frozenset())
    monkeypatch.setattr(spec, "DEFERRED_MODULES", frozenset())
    monkeypatch.setattr(app_module, "_modules", {m["id"]: m for m in spec.MODULES})
    yield


@pytest.fixture(autouse=True)
def reports_in_tmp(tmp_path, monkeypatch):
    """Feedback and crash reports go to a throwaway folder in tests, never to ~/.mzdesk; nothing is forwarded; the admin lock starts fresh."""
    from mzdesk.api import admin
    monkeypatch.setenv("MZDESK_FEEDBACK_DIR", str(tmp_path / "feedback"))
    for name in ("MZDESK_RESEND_API_KEY", "MZDESK_RESEND_URL", "MZDESK_MAIL_FROM", "MZDESK_NTFY_URL", "MZDESK_FEEDBACK_TO", "MZDESK_ADMIN_SALT", "MZDESK_ADMIN_HASH",
                 "MZDESK_SMTP_HOST", "MZDESK_SMTP_PORT", "MZDESK_SMTP_USER", "MZDESK_SMTP_PASSWORD", "MZDESK_SMTP_FROM", "MZDESK_SMTP_SECURITY"):
        monkeypatch.delenv(name, raising=False)
    monkeypatch.delenv("MZDESK_ANNOTATION_FRIGID", raising=False)          # the local default: MS/MS Annotation may run FRIGID
    admin.reset_throttle()
    yield


@pytest.fixture(autouse=True)
def a_machine_that_runs_every_component(monkeypatch):
    """The component tests describe a computer where every component is supported. On an Intel Mac FRIGID is not (components.py), which
    would change their answers; the tests of that rule set the machine themselves."""
    from mzdesk.local import components
    monkeypatch.setattr(components, "current_machine", lambda: "arm64")


# ---------------------------------------------------------------- tests that need `python dev_setup.py`
# They use what a fresh git checkout does not have: the MetFrag jar, the built web interface, or a component (FRIGID, RepoRT, MolGpKa,
# Pan-ReDU, the spectral libraries). Found by running the suite on clean Linux and macOS machines (CI, 2026-11-13). With
# MZDESK_TESTS_WITHOUT_SETUP=1 (the CI run on every push) they are skipped; everywhere else they run, and the CI "full" run installs
# everything with dev_setup.py first and runs them all. Prefixes of pytest node ids.
NEEDS_SETUP = (
    'tests/test_adducts.py::test_likelihood_badge_tones',
    'tests/test_annotation_ccs_unchanged.py::test_a_run_without_an_observed_ccs_is_what_it_was_before[real re-rank]',   # MolGpKa pKa in the reasons
    'tests/test_api.py::test_example_runs_and_fills_every_block[gnps_redu]',
    'tests/test_api.py::test_example_runs_and_fills_every_block[ms2_fragments]',
    'tests/test_api.py::test_metfrag_matches_engine',
    'tests/test_fragment_view.py::test_a_structure_with_nothing_explained_still_gives_a_view',
    'tests/test_fragment_view.py::test_charged_and_adduct_structures_are_explained_and_placed[C[N+](C)(C)CCO-[M]+-104.107-peaks0]',
    'tests/test_fragment_view.py::test_charged_and_adduct_structures_are_explained_and_placed[Cn1c(=O)c2c(ncn2C)n(C)c1=O-[M+Na]+-217.0696-peaks2]',
    'tests/test_fragment_view.py::test_charged_and_adduct_structures_are_explained_and_placed[OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O-[M+HCOO]--225.0616-peaks1]',
    'tests/test_fragment_view.py::test_real_job_end_to_end',
    'tests/test_fragment_view.py::test_real_view_of_the_flavan_dimer',
    'tests/test_fragment_view.py::test_the_coverage_in_the_table_is_the_one_in_the_view',
    'tests/test_known_methods.py::test_page_for_captopril_jp_only',
    'tests/test_known_methods.py::test_similar_list_leaves_out_the_compound_itself_and_scores_stay_in_range',
    'tests/test_security.py::test_every_answer_carries_the_security_headers',
    'tests/test_spec_copy.py::test_the_example_of_every_shown_page_runs_with_short_plain_text[gnps_redu]',
    'tests/test_spec_copy.py::test_the_example_of_every_shown_page_runs_with_short_plain_text[ionization]',
    'tests/test_spec_copy.py::test_the_example_of_every_shown_page_runs_with_short_plain_text[ms2_fragments]',
    'tests/test_versions.py::test_nothing_installed_gives_null_versions',
    'tests/test_workflow_lipid_fragments.py::test_a_far_precursor_falls_back_with_one_warning',
    'tests/test_workflow_lipid_fragments.py::test_a_name_runs_end_to_end',
    'tests/test_workflow_lipid_fragments.py::test_a_smiles_runs_with_dash_names_and_the_same_explanation',
    'tests/test_workflow_lipid_fragments.py::test_the_job_api_runs_the_workflow',
    'tests/test_workflow_targeted_method.py::test_missing_molgpka_weights_skip_separation_with_a_note_and_the_sheet_goes_on',
    'tests/test_workflow_targeted_method.py::test_on_the_gc_route_the_published_methods_are_the_gc_ones',
)


def pytest_collection_modifyitems(config, items):
    if os.environ.get("MZDESK_TESTS_WITHOUT_SETUP") != "1":
        return
    skip = pytest.mark.skip(reason="needs `python dev_setup.py` (skipped because MZDESK_TESTS_WITHOUT_SETUP=1)")
    for item in items:
        if item.nodeid.startswith(NEEDS_SETUP):
            item.add_marker(skip)
