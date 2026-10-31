"""Opt-in browser tests: MZDESK_UI_TEST_URL=http://localhost:5173 pytest tests/test_ui_workflow_progress.py.

The workflow API is mocked; no engines, network lookups or server jobs are started.
"""
import os
import re
import time

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


def percent(browser):
    """aria-valuenow of the single progress bar, as a number (None while the bar is indeterminate)."""
    from selenium.webdriver.common.by import By
    value = progress(browser).find_element(By.CSS_SELECTOR, '[role="progressbar"]').get_attribute("aria-valuenow")
    return None if value is None else float(value)


def release_start(browser):
    wait(browser, lambda d: d.execute_script("return !!window.__releaseProgressStart"))
    browser.execute_script("window.__holdProgressStart = false; window.__releaseProgressStart()")


def set_job(browser, job):
    browser.execute_script("window.__progressJob = arguments[0]", job)


def test_progress_is_one_block_that_follows_the_job_and_never_goes_back(browser):
    from selenium.webdriver.common.by import By

    browser.set_window_size(390, 844)
    open_workflow(browser, snapshot("queued", ["queued"] * 8, "waiting for the previous run (position 2 in the queue)"))
    wait_state(browser, "Submitting workflow")
    assert 2 <= percent(browser) < 6, "the bar has a head start while the run is submitted"
    release_start(browser)
    wait_state(browser, "Queued")
    assert "position 2 in the queue" in progress(browser).text
    assert 2 <= percent(browser) < 6
    seen = [percent(browser)]

    set_job(browser, snapshot("running", ["running", "running"] + ["queued"] * 6))
    wait_state(browser, "Adduct + Formula")
    assert "Steps 1, 2 of 8" in progress(browser).text
    seen.append(percent(browser))
    stage6 = snapshot("running", ["done", "done", "skipped", "failed", "done", "running", "queued", "queued"])
    stage6["steps"][5]["note"] = "scored 3/9"
    stage6["steps"][1]["blocks"] = [{"id": "formula", "type": "table", "title": "Formula detail", "columns": ["Formula"],
                                     "data": {"rows": [["C8H10N4O2"]]}}]
    set_job(browser, stage6)
    wait_state(browser, "Step 6 of 8 \u00b7 Merge and re-rank")
    area = progress(browser)
    assert "Merge and re-rank" in area.text and "scored 3/9" in area.text
    assert "elapsed-time" not in area.text and "Stage-count" not in area.text
    assert not area.find_elements(By.CSS_SELECTOR, "ol, li, a"), "no stage card grid"
    wait(browser, lambda d: percent(d) >= 59)            # the bar glides to its target at no more than 15 points per second
    first = percent(browser)
    assert 59 <= first < 73.2, "four of the seven stages that were not skipped have ended: inside the fifth of them"
    assert re.fullmatch(r"\d+%", area.find_element(By.CSS_SELECTOR, ".wf-progress-percent").text)
    assert "scored 3/9" in browser.find_element(By.ID, "workflow-step-6").text
    assert "C8H10N4O2" not in browser.find_element(By.ID, "workflow-step-2").text, "steps are folded: the table is not shown yet"
    browser.find_element(By.CSS_SELECTOR, "#workflow-step-2 > summary").click()
    assert "C8H10N4O2" in browser.find_element(By.ID, "workflow-step-2").text
    assert browser.execute_script("return document.documentElement.scrollWidth <= window.innerWidth")
    seen.append(first)
    wait(browser, lambda d: percent(d) > first)          # the running stage keeps filling its slice between polls
    seen.append(percent(browser))
    set_job(browser, snapshot("running", ["done", "done", "skipped", "failed", "done", "done", "queued", "queued"]))
    wait(browser, lambda d: percent(d) > seen[-1])
    seen.append(percent(browser))
    assert seen == sorted(seen), f"progress went backwards: {seen}"

    ended = ["done", "done", "skipped", "failed", "done", "done", "done", "done"]
    set_job(browser, snapshot("running", ended))
    wait_state(browser, "Finalizing result")
    wait(browser, lambda d: percent(d) == 99)
    set_job(browser, snapshot("done", ended, result=result(ended)))
    wait_state(browser, "Completed with 1 failed stage")
    assert percent(browser) == 100 and "100%" in progress(browser).text
    assert "Workflow result retained" in browser.find_element(By.CSS_SELECTOR, '[data-pipeline="esi_ms2_annotation"]').text
    assert len(browser.find_elements(By.CSS_SELECTOR, '[data-pipeline="esi_ms2_annotation"] [role="progressbar"]')) == 1
    assert len(browser.find_elements(By.CSS_SELECTOR, '[data-pipeline="esi_ms2_annotation"] .wf-progress')) == 1
    clear = button(browser, "Clear")
    browser.execute_script("arguments[0].scrollIntoView({block: 'center'})", clear)
    clear.click()
    assert not browser.find_elements(By.CSS_SELECTOR, '[data-pipeline="esi_ms2_annotation"] .wf-progress')


def test_stopped_polling_marks_the_snapshot_as_stale(browser):
    open_workflow(browser, snapshot("running", ["done", "running"] + ["queued"] * 6))
    release_start(browser)
    wait_state(browser, "Step 2 of 8 \u00b7 Formula")
    wait(browser, lambda d: percent(d) >= 18)             # one stage ended and the second one is filling
    browser.execute_script("window.__failProgressPoll = true")
    wait_state(browser, "Tracking stopped")
    assert "last update" in progress(browser).text
    frozen = percent(browser)
    assert frozen >= 18, "one stage ended: the bar keeps what it had"
    time.sleep(1.2)
    assert percent(browser) == frozen, "a stopped bar does not move"
    assert "Step 2 of 8" not in progress(browser).text


def test_poll_failure_before_first_snapshot_does_not_invent_stages(browser):
    open_workflow(browser, snapshot("running", ["running"] + ["queued"] * 7))
    browser.execute_script("window.__failProgressPoll = true")
    release_start(browser)
    wait_state(browser, "Tracking stopped")
    assert "Step 1 of" not in progress(browser).text
    assert 2 <= percent(browser) < 6, "only the head start is shown"


def test_transient_poll_failure_is_retried(browser):
    open_workflow(browser, snapshot("running", ["running"] + ["queued"] * 7))
    browser.execute_script("window.__transientProgressFailures = 1")
    release_start(browser)
    wait_state(browser, "Waiting for workflow progress")
    wait_state(browser, "Step 1 of 8 \u00b7 Adduct")
    assert "Tracking stopped" not in progress(browser).text


def test_job_error_is_distinct_from_completed_stages(browser):
    open_workflow(browser, snapshot("error", ["failed"] + ["queued"] * 7, error="engine exploded"))
    release_start(browser)
    wait_state(browser, "Failed at Adduct")
    assert 2 <= percent(browser) < 99, "an error leaves the bar where it was"
    assert "Completed" not in progress(browser).text


def test_unknown_or_missing_data_does_not_fabricate_progress(browser):
    from selenium.webdriver.common.by import By

    first = snapshot("running", ["running", "queued"])
    first["steps"][0]["seconds"] = None
    first["steps"][1]["status"] = None
    open_workflow(browser, first)
    release_start(browser)
    wait_state(browser, "stage details unavailable")
    assert percent(browser) is None
    assert " s" not in browser.find_element(By.ID, "workflow-step-1").text

    set_job(browser, {"status": "running", "steps": None, "note": "", "error": None, "result": None})
    wait_state(browser, "stage details unavailable")
    assert percent(browser) is None

    set_job(browser, snapshot("unexpected", ["running", "queued"]))
    wait_state(browser, "Status unavailable")
    assert "Step 1 of" not in progress(browser).text
    button(browser, "Run workflow")


def test_clear_ignores_a_late_poll_response(browser):
    from selenium.webdriver.common.by import By

    open_workflow(browser, snapshot("running", ["running"] + ["queued"] * 7))
    release_start(browser)
    wait_state(browser, "Step 1 of 8 \u00b7 Adduct")
    browser.execute_script("window.__holdProgressPoll = true")
    wait(browser, lambda d: d.execute_script("return !!window.__releaseProgressPoll"))
    button(browser, "Clear").click()
    browser.execute_script("window.__releaseProgressPoll()")
    wait(browser, lambda d: d.execute_script("return window.__progressPollReleased"))
    assert not browser.find_elements(By.CSS_SELECTOR, '[data-pipeline="esi_ms2_annotation"] .wf-progress')
    assert "previous workflow may still finish" in browser.find_element(By.CSS_SELECTOR, '[data-pipeline="esi_ms2_annotation"]').text
