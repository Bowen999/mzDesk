"""Library download: the MassBank release is found even when the GitHub API refuses (rate limit on a shared address)."""
import io
import urllib.error

import pytest

from mzdesk.ms2 import library_setup as ls


class _Page(io.BytesIO):
    def __init__(self, url):
        super().__init__(b"")
        self._url = url

    def geturl(self):
        return self._url

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _fake(code):
    def urlopen(req, timeout=None):
        if req.full_url == ls.MASSBANK_API:
            raise urllib.error.HTTPError(req.full_url, code, "rate limit exceeded", {}, None)
        assert req.full_url == "https://github.com/MassBank/MassBank-data/releases/latest"
        return _Page("https://github.com/MassBank/MassBank-data/releases/tag/2026.11")
    return urlopen


@pytest.mark.parametrize("code", [403, 429])
def test_a_refused_api_falls_back_to_the_release_page(monkeypatch, code):
    monkeypatch.setattr(ls.urllib.request, "urlopen", _fake(code))
    info = ls.source_info("MassBank")
    assert info == {"url": "https://github.com/MassBank/MassBank-data/releases/download/2026.11/MassBank_NISTformat.msp",
                    "file": "MassBank_NISTformat.msp", "version": "2026.11", "date": "", "size": None}


def test_other_api_errors_are_not_hidden(monkeypatch):
    monkeypatch.setattr(ls.urllib.request, "urlopen", _fake(500))
    with pytest.raises(urllib.error.HTTPError):
        ls.source_info("MassBank")
