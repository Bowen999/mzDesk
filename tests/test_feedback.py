"""Feedback endpoint: validation, the store, the forwarding channels (email through Resend / SMTP, ntfy) and the honeypot.

Nothing here touches the network: the Resend and ntfy channels talk to a small local HTTP server, SMTP is a fake smtplib.SMTP.
"""
import base64
import json
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer

import pytest
from fastapi.testclient import TestClient

from mzdesk.api import feedback, reports
from mzdesk.api.app import app

client = TestClient(app)

# a real 1x1 PNG: the server sniffs the bytes instead of trusting the data-URL type
PNG = base64.b64decode(
    b"iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+A8AAQUBAScY42YAAAAASUVORK5CYII="
)
PNG_URL = "data:image/png;base64," + base64.b64encode(PNG).decode()


REAL_FORWARD_ASYNC = reports.forward_async


@pytest.fixture(autouse=True)
def isolated(monkeypatch):
    """No rate-limit history between tests (the folder and the mail settings are reset by conftest.reports_in_tmp).
    The background forward is recorded instead of started, so each test sends exactly what it asks for."""
    feedback._HITS.clear()
    feedback._CRASH_HITS.clear()
    queued = []
    monkeypatch.setattr(reports, "forward_async", queued.append)
    yield queued
    feedback._HITS.clear()
    feedback._CRASH_HITS.clear()


def _payload(**over):
    body = {"message": "The Run button does nothing on the MS1 search page.",
            "page": "/modules/ms1_search", "website": "", "images": []}
    body.update(over)
    return body


class Sink:
    """A local HTTP server that records what it is sent (stands in for api.resend.com and ntfy.sh)."""

    def __init__(self, status=200):
        self.requests = []
        sink = self

        class Handler(BaseHTTPRequestHandler):
            def do_POST(self):
                body = self.rfile.read(int(self.headers.get("Content-Length", 0)))
                sink.requests.append({"path": self.path, "headers": dict(self.headers), "body": body})
                self.send_response(status)
                self.end_headers()
                self.wfile.write(b"{}")

            def log_message(self, *args):
                pass

        self.server = HTTPServer(("127.0.0.1", 0), Handler)
        self.url = f"http://127.0.0.1:{self.server.server_port}"
        threading.Thread(target=self.server.serve_forever, daemon=True).start()

    def close(self):
        self.server.shutdown()
        self.server.server_close()


@pytest.fixture
def sink():
    s = Sink()
    yield s
    s.close()


def test_feedback_is_stored_with_its_screenshot():
    r = client.post("/api/feedback", json=_payload(images=[{"name": "shot.png", "data": PNG_URL}]))
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["ok"] and body["delivered"] and body["note"] == "" and body["ref"].startswith("F-")

    [row] = reports.listing()
    assert row["id"] == body["ref"] and row["kind"] == "feedback" and row["page"] == "/modules/ms1_search" and row["images"] == 1
    full = reports.detail(body["ref"])
    assert "The Run button does nothing" in full["body"]
    assert full["images"] == [{"name": "shot.png", "src": PNG_URL}]


def test_a_report_is_stored_before_anything_is_forwarded_and_survives_a_dead_mail_service(monkeypatch):
    monkeypatch.setenv("MZDESK_RESEND_API_KEY", "re_test")
    monkeypatch.setenv("MZDESK_RESEND_URL", "http://127.0.0.1:9/emails")           # nothing listens there
    ref = client.post("/api/feedback", json=_payload()).json()["ref"]
    assert reports.forward(ref) == []                                              # the channel failed...
    assert reports.detail(ref)["body"].startswith("The Run button")               # ...and the report is still there


def test_email_through_resend(monkeypatch, sink):
    monkeypatch.setenv("MZDESK_RESEND_API_KEY", "re_test_key")
    monkeypatch.setenv("MZDESK_RESEND_URL", sink.url + "/emails")
    monkeypatch.setenv("MZDESK_FEEDBACK_TO", "me@example.org")
    ref = client.post("/api/feedback", json=_payload(images=[{"name": "shot.png", "data": PNG_URL}])).json()["ref"]
    assert reports.forward(ref) == ["resend"]
    [req] = sink.requests
    assert req["headers"]["Authorization"] == "Bearer re_test_key" and req["headers"]["User-Agent"].startswith("mzDesk/")
    mail = json.loads(req["body"])
    assert mail["to"] == ["me@example.org"] and mail["from"] == "mzDesk <onboarding@resend.dev>"
    assert mail["subject"] == f"[mzDesk] Feedback {ref}: The Run button does nothing on the MS1 search page."
    assert "Page: /modules/ms1_search" in mail["text"]
    assert base64.b64decode(mail["attachments"][0]["content"]) == PNG and mail["attachments"][0]["filename"] == "shot.png"


def test_a_new_report_is_queued_for_forwarding_and_the_background_thread_sends_it(monkeypatch, sink, isolated):
    monkeypatch.setenv("MZDESK_RESEND_API_KEY", "re_test_key")
    monkeypatch.setenv("MZDESK_RESEND_URL", sink.url + "/emails")
    ref = client.post("/api/feedback", json=_payload()).json()["ref"]
    assert isolated == [ref]
    threads = []
    real = threading.Thread
    monkeypatch.setattr(reports.threading, "Thread", lambda **kw: threads.append(real(**kw)) or threads[-1])
    REAL_FORWARD_ASYNC(ref)                                    # starts its daemon thread itself
    assert len(threads) == 1 and threads[0].daemon
    threads[0].join(5)
    assert len(sink.requests) == 1


def test_without_a_channel_nothing_is_started(monkeypatch):
    started = []
    monkeypatch.setattr(reports.threading, "Thread", lambda **kw: started.append(kw))
    REAL_FORWARD_ASYNC("F-000000")
    assert started == [] and reports.channels() == []


def test_push_through_ntfy_carries_no_message_text(monkeypatch, sink):
    monkeypatch.setenv("MZDESK_NTFY_URL", sink.url + "/secret-topic")
    ref = client.post("/api/feedback", json=_payload(message="Private result: the unknown at 301.141 is my patient's drug.")).json()["ref"]
    assert reports.forward(ref) == ["ntfy"]
    [req] = sink.requests
    assert req["path"] == "/secret-topic" and req["headers"]["Title"] == "mzDesk feedback"
    assert ref in req["body"].decode() and b"patient" not in req["body"]


def test_email_through_smtp(monkeypatch):
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

    monkeypatch.setattr(reports.smtplib, "SMTP", FakeSMTP)
    monkeypatch.setenv("MZDESK_SMTP_HOST", "smtp.example.org")
    monkeypatch.setenv("MZDESK_SMTP_USER", "mailer@example.org")
    monkeypatch.setenv("MZDESK_SMTP_PASSWORD", "secret")
    ref = client.post("/api/feedback", json=_payload()).json()["ref"]
    assert reports.channels() == ["smtp"] and reports.forward(ref) == ["smtp"]
    assert sent[0]["To"] == "feedback@example.org" and sent[0]["From"] == "mailer@example.org"


def test_a_failing_channel_does_not_stop_the_next_one(monkeypatch, sink):
    monkeypatch.setenv("MZDESK_RESEND_API_KEY", "re_test_key")
    monkeypatch.setenv("MZDESK_RESEND_URL", "http://127.0.0.1:9/emails")
    monkeypatch.setenv("MZDESK_NTFY_URL", sink.url + "/t")
    ref = client.post("/api/feedback", json=_payload()).json()["ref"]
    assert reports.forward(ref) == ["ntfy"]


def test_honeypot_pretends_success_and_stores_nothing():
    r = client.post("/api/feedback", json=_payload(website="https://spam.example"))
    assert r.status_code == 200
    assert r.json()["delivered"] is True
    assert reports.listing() == []


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
    assert reports.listing() == []


def test_rate_limit_after_the_first_message():
    assert client.post("/api/feedback", json=_payload()).status_code == 200
    r = client.post("/api/feedback", json=_payload(message="Another bug report."))
    assert r.status_code == 429
    assert "wait" in r.json()["detail"].lower()
    assert len(reports.listing()) == 1


def test_the_store_folder_can_be_configured(monkeypatch, tmp_path):
    target = tmp_path / "elsewhere"
    monkeypatch.setenv("MZDESK_FEEDBACK_DIR", str(target))
    assert reports.spool_dir() == target
    assert client.post("/api/feedback", json=_payload()).status_code == 200
    assert (target / "reports.sqlite").is_file()
