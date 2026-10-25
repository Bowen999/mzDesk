"""Shared test setup: the page examples are compound names; resolve them from spec.EXAMPLE_STRUCTURES so tests stay offline."""
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
