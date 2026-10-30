"""Feedback endpoint: validation, the .eml fallback, SMTP delivery and the honeypot.

The endpoint is network-free here: SMTP is either unset (message goes to the spool folder) or a fake smtplib.SMTP.
"""
import base64
import email
import email.policy
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from mzdesk.api import feedback
from mzdesk.api.app import app

client = TestClient(app)

# a real 1x1 PNG: the server sniffs the bytes instead of trusting the data-URL type
PNG = base64.b64decode(
    b"iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+A8AAQUBAScY42YAAAAASUVORK5CYII="
)
PNG_URL = "data:image/png;base64," + base64.b64encode(PNG).decode()


@pytest.fixture(autouse=True)
def isolated(tmp_path, monkeypatch):
    """No real SMTP, a private spool folder, and no rate-limit history between tests."""
    feedback._HITS.clear()
    monkeypatch.setenv("MZDESK_FEEDBACK_DIR", str(tmp_path / "feedback"))
    for name in ("MZDESK_SMTP_HOST", "MZDESK_SMTP_PORT", "MZDESK_SMTP_USER", "MZDESK_SMTP_PASSWORD",
                 "MZDESK_SMTP_FROM", "MZDESK_SMTP_SECURITY", "MZDESK_FEEDBACK_TO"):
        monkeypatch.delenv(name, raising=False)
    yield
    feedback._HITS.clear()


def _payload(**over):
    body = {"message": "The Run button does nothing on the MS1 search page.",
            "page": "/modules/ms1_search", "website": "", "images": []}
    body.update(over)
    return body


def _spooled():
    return sorted(feedback.spool_dir().glob("*.eml"))


def test_feedback_is_saved_as_eml_without_smtp():
    r = client.post("/api/feedback", json=_payload(images=[{"name": "shot.png", "data": PNG_URL}]))
    assert r.status_code == 200, r.text
    assert r.json() == {"ok": True, "delivered": False,
                        "note": "Your message was saved on the server; email delivery is not configured yet."}

    files = _spooled()
    assert len(files) == 1
    msg = email.message_from_bytes(files[0].read_bytes(), policy=email.policy.default)
    assert msg["To"] == "feedback@example.org"
    assert "The Run button does nothing" in msg.get_body(preferencelist=("plain",)).get_content()
    attachments = list(msg.iter_attachments())
    assert len(attachments) == 1
    assert attachments[0].get_filename() == "shot.png"
    assert attachments[0].get_content_type() == "image/png"
    assert attachments[0].get_content() == PNG


def test_feedback_is_emailed_when_smtp_is_configured(monkeypatch):
    sent = []

    class FakeSMTP:
        def __init__(self, host, port, timeout=None):
            self.host, self.port = host, port

        def starttls(self, context=None):
            pass

        def login(self, user, password):
            assert (user, password) == ("mailer@example.org", "secret")

        def send_message(self, msg):
            sent.append(msg)

        def quit(self):
            pass

    monkeypatch.setattr(feedback.smtplib, "SMTP", FakeSMTP)
    monkeypatch.setenv("MZDESK_SMTP_HOST", "smtp.example.org")
    monkeypatch.setenv("MZDESK_SMTP_USER", "mailer@example.org")
    monkeypatch.setenv("MZDESK_SMTP_PASSWORD", "secret")

    r = client.post("/api/feedback", json=_payload())
    assert r.status_code == 200, r.text
    assert r.json()["delivered"] is True
    assert len(sent) == 1
    assert sent[0]["To"] == "feedback@example.org"
    assert sent[0]["From"] == "mailer@example.org"
    assert not _spooled()


def test_smtp_failure_keeps_the_message(monkeypatch):
    class BrokenSMTP:
        def __init__(self, *a, **kw):
            raise OSError("connection refused")

    monkeypatch.setattr(feedback.smtplib, "SMTP", BrokenSMTP)
    monkeypatch.setenv("MZDESK_SMTP_HOST", "smtp.example.org")
    monkeypatch.setenv("MZDESK_SMTP_USER", "mailer@example.org")
    monkeypatch.setenv("MZDESK_SMTP_PASSWORD", "secret")

    r = client.post("/api/feedback", json=_payload())
    assert r.status_code == 200
    assert r.json()["delivered"] is False
    assert "saved" in r.json()["note"]
    assert len(_spooled()) == 1


def test_honeypot_pretends_success_and_sends_nothing():
    r = client.post("/api/feedback", json=_payload(website="https://spam.example"))
    assert r.status_code == 200
    assert r.json()["delivered"] is True
    assert not _spooled()


@pytest.mark.parametrize("body", [
    _payload(message=""),
    _payload(message="  x "),
    _payload(message="x" * 4001),
    _payload(images=[{"name": f"{i}.png", "data": PNG_URL} for i in range(4)]),
    _payload(images=[{"name": "x.svg", "data": "data:image/svg+xml;base64," + base64.b64encode(b"<svg/>").decode()}]),
    _payload(images=[{"name": "x.png", "data": "data:image/png;base64,not base64!"}]),
])
def test_bad_feedback_is_400(body):
    r = client.post("/api/feedback", json=body)
    assert r.status_code == 400, r.text
    assert not _spooled()


def test_rate_limit_after_the_first_message():
    assert client.post("/api/feedback", json=_payload()).status_code == 200
    r = client.post("/api/feedback", json=_payload(message="Another bug report."))
    assert r.status_code == 429
    assert "wait" in r.json()["detail"].lower()
    assert len(_spooled()) == 1


def test_feedback_dir_can_be_configured(monkeypatch, tmp_path):
    target = tmp_path / "elsewhere"
    monkeypatch.setenv("MZDESK_FEEDBACK_DIR", str(target))
    assert feedback.spool_dir() == Path(target)
    r = client.post("/api/feedback", json=_payload())
    assert r.status_code == 200
    assert len(list(target.glob("*.eml"))) == 1
