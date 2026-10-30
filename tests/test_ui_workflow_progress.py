"""Opt-in browser tests: MZDESK_UI_TEST_URL=http://localhost:5173 pytest tests/test_ui_workflow_progress.py.

The workflow API is mocked; no engines, network lookups or server jobs are started.
"""
import os

import pytest

from test_ui_navigation import browser, button


URL = os.environ.get("MZDESK_UI_TEST_URL")
pytestmark = pytest.mark.skipif(not URL, reason="set MZDESK_UI_TEST_URL to a running frontend")
TITLES = ["Adduct", "Formula", "Library search", "MS-FINDER structures", "FRIGID structures", "Merge and re-rank", "Isomer check", "Answer"]


def snapshot(status, statuses, note="", error=None, result=None):
    return {"status": status, "note": note, "error": error, "result": result,
            "steps": [{"id": i + 1, "title": TITLES[i], "status": s, "seconds": 0, "note": "", "blocks": []} for i, s in enumerate(statuses)]}


def result(statuses):
    return {"answer": {"type": "none", "ids": [], "statement": "Workflow result retained", "top_id": None, "isomer_verdict": "not run"},
            "candidates": [], "warnings": [], "context": {"precursor_mz": 195.0877, "polarity": "Positive", "gate_passed": False}, "isomer": None,
            "steps": snapshot("done", statuses)["steps"]}


def open_workflow(browser, first):
    browser.get(URL.rstrip("/") + "/pipelines/esi_ms2_annotation")
    button(browser, "Load example: Caffeine, library match").click()
    browser.execute_script("""
        window.__progressJob = arguments[0];
        window.__holdProgressStart = true;
        window.__releaseProgressStart = null;
        window.__failProgressPoll = false;
        window.__transientProgressFailures = 0;
        window.__holdProgressPoll = false;
        window.__releaseProgressPoll = null;
        window.__progressPollReleased = false;
        const originalFetch = window.fetch.bind(window);
        const response = (data, status = 200) => new Response(JSON.stringify(data), {status, headers: {'Content-Type': 'application/json'}});
        window.fetch = async (input, init) => {
          const url = new URL(typeof input === 'string' ? input : input.url, location.href);
          if (url.pathname.endsWith('/api/pipelines/esi_ms2_annotation/run')) {
            if (window.__holdProgressStart) await new Promise(resolve => { window.__releaseProgressStart = resolve; });
            return response({job_id: 'progress-job', status: 'running'});
          }
          if (url.pathname.endsWith('/api/pipelines/jobs/progress-job')) {
            if (window.__holdProgressPoll) {
              await new Promise(resolve => { window.__releaseProgressPoll = resolve; });
              window.__progressPollReleased = true;
            }
            if (window.__transientProgressFailures > 0) {
              window.__transientProgressFailures--;
              throw new TypeError('temporary network failure');
            }
            if (window.__failProgressPoll) return response({detail: 'job unavailable'}, 404);
            return response(window.__progressJob);
          }
          return originalFetch(input, init);
        };
    """, first)
    button(browser, "Run workflow").click()


def wait(browser, predicate):
    from selenium.webdriver.support.ui import WebDriverWait
    return WebDriverWait(browser, 12).until(predicate)


def progress(browser):
    from selenium.webdriver.common.by import By
    return browser.find_element(By.CSS_SELECTOR, '[data-pipeline="esi_ms2_annotation"] .wf-progress')


def wait_state(browser, text):
    return wait(browser, lambda d: text in progress(d).text)


def release_start(browser):
    wait(browser, lambda d: d.execute_script("return !!window.__releaseProgressStart"))
    browser.execute_script("window.__holdProgressStart = false; window.__releaseProgressStart()")


def set_job(browser, job):
    browser.execute_script("window.__progressJob = arguments[0]", job)


def test_progress_follows_the_job_and_preserves_step_details(browser):
    from selenium.webdriver.common.by import By

    browser.set_window_size(390, 844)
    open_workflow(browser, snapshot("queued", ["queued"] * 8, "waiting for the previous run (position 2 in the queue)"))
    wait_state(browser, "Submitting workflow")
    bar = progress(browser).find_element(By.CSS_SELECTOR, '[role="progressbar"]')
    assert bar.get_attribute("aria-valuenow") is None
    release_start(browser)
    wait_state(browser, "position 2 in the queue")
    assert progress(browser).find_element(By.CSS_SELECTOR, '[role="progressbar"]').get_attribute("aria-valuenow") == "0"

    set_job(browser, snapshot("running", ["running", "running"] + ["queued"] * 6))
    wait_state(browser, "Running: Adduct + Formula")
    stage6 = snapshot("running", ["done", "done", "skipped", "failed", "done", "running", "queued", "queued"])
    stage6["steps"][5]["note"] = "scored 3/9"
    stage6["steps"][1]["blocks"] = [{"id": "formula", "type": "table", "title": "Formula detail", "columns": ["Formula"],
                                     "data": {"rows": [["C8H10N4O2"]]}}]
    set_job(browser, stage6)
    wait_state(browser, "62.5%")
    area = progress(browser)
    assert "Merge and re-rank: scored 3/9" in area.text
    assert "not elapsed-time progress" in area.text
    assert area.find_element(By.CSS_SELECTOR, '[role="progressbar"]').get_attribute("aria-valuenow") == "62.5"
    assert area.find_element(By.CSS_SELECTOR, 'a[href="#workflow-step-6"]')
    assert "scored 3/9" in browser.find_element(By.ID, "workflow-step-6").text
    assert "C8H10N4O2" in browser.find_element(By.ID, "workflow-step-2").text
    area.find_element(By.CSS_SELECTOR, 'a[href="#workflow-step-2"]').click()
    assert browser.execute_script("return location.hash") == "#workflow-step-2"
    assert browser.execute_script("return document.documentElement.scrollWidth <= window.innerWidth")

    ended = ["done", "done", "skipped", "failed", "done", "done", "done", "done"]
    set_job(browser, snapshot("running", ended))
    wait_state(browser, "Finalizing result")
    assert "100%" in progress(browser).text
    set_job(browser, snapshot("done", ended, result=result(ended)))
    wait_state(browser, "Workflow completed with 1 failed stage")
    assert "Workflow result retained" in browser.find_element(By.CSS_SELECTOR, '[data-pipeline="esi_ms2_annotation"]').text
    assert len(browser.find_elements(By.CSS_SELECTOR, '[data-pipeline="esi_ms2_annotation"] [role="progressbar"]')) == 1
    clear = button(browser, "Clear")
    browser.execute_script("arguments[0].scrollIntoView({block: 'center'})", clear)
    clear.click()
    assert not browser.find_elements(By.CSS_SELECTOR, '[data-pipeline="esi_ms2_annotation"] .wf-progress')


def test_stopped_polling_marks_the_snapshot_as_stale(browser):
    open_workflow(browser, snapshot("running", ["done", "running"] + ["queued"] * 6))
    release_start(browser)
    wait_state(browser, "Running: Formula")
    browser.execute_script("window.__failProgressPoll = true")
    wait_state(browser, "Progress tracking stopped")
    assert "last update" in progress(browser).text
    assert "12.5%" in progress(browser).text
    assert "Running: Formula" not in progress(browser).text


def test_poll_failure_before_first_snapshot_does_not_invent_progress(browser):
    from selenium.webdriver.common.by import By

    open_workflow(browser, snapshot("running", ["running"] + ["queued"] * 7))
    browser.execute_script("window.__failProgressPoll = true")
    release_start(browser)
    wait_state(browser, "Progress tracking stopped")
    assert "Stage progress unavailable" in progress(browser).text
    assert progress(browser).find_element(By.CSS_SELECTOR, '[role="progressbar"]').get_attribute("aria-valuenow") is None


def test_transient_poll_failure_is_retried(browser):
    open_workflow(browser, snapshot("running", ["running"] + ["queued"] * 7))
    browser.execute_script("window.__transientProgressFailures = 1")
    release_start(browser)
    wait_state(browser, "Waiting for workflow progress")
    wait_state(browser, "Running: Adduct")
    assert "Progress tracking stopped" not in progress(browser).text


def test_job_error_is_distinct_from_completed_stages(browser):
    open_workflow(browser, snapshot("error", ["failed"] + ["queued"] * 7, error="engine exploded"))
    release_start(browser)
    wait_state(browser, "Workflow failed")
    assert "12.5%" in progress(browser).text
    assert "Workflow completed" not in progress(browser).text


def test_unknown_or_missing_data_does_not_fabricate_progress(browser):
    from selenium.webdriver.common.by import By

    first = snapshot("running", ["running", "queued"])
    first["steps"][0]["seconds"] = None
    first["steps"][1]["status"] = None
    open_workflow(browser, first)
    release_start(browser)
    wait_state(browser, "Running: Adduct")
    assert "unknown" in progress(browser).text
    assert progress(browser).find_element(By.CSS_SELECTOR, '[role="progressbar"]').get_attribute("aria-valuenow") is None
    assert " s" not in browser.find_element(By.ID, "workflow-step-1").text

    set_job(browser, {"status": "running", "steps": None, "note": "", "error": None, "result": None})
    wait_state(browser, "stage details unavailable")
    assert progress(browser).find_element(By.CSS_SELECTOR, '[role="progressbar"]').get_attribute("aria-valuenow") is None

    set_job(browser, snapshot("unexpected", ["running", "queued"]))
    wait_state(browser, "Workflow status unavailable; progress tracking stopped")
    assert "Running: Adduct" not in progress(browser).text
    button(browser, "Run workflow")


def test_clear_ignores_a_late_poll_response(browser):
    from selenium.webdriver.common.by import By

    open_workflow(browser, snapshot("running", ["running"] + ["queued"] * 7))
    release_start(browser)
    wait_state(browser, "Running: Adduct")
    browser.execute_script("window.__holdProgressPoll = true")
    wait(browser, lambda d: d.execute_script("return !!window.__releaseProgressPoll"))
    button(browser, "Clear").click()
    browser.execute_script("window.__releaseProgressPoll()")
    wait(browser, lambda d: d.execute_script("return window.__progressPollReleased"))
    assert not browser.find_elements(By.CSS_SELECTOR, '[data-pipeline="esi_ms2_annotation"] .wf-progress')
    assert "previous workflow may still finish" in browser.find_element(By.CSS_SELECTOR, '[data-pipeline="esi_ms2_annotation"]').text
