"""Opt-in browser regression: MZDESK_UI_TEST_URL=http://localhost:5173 pytest tests/test_ui_navigation.py.

The API is mocked in the browser; no real online lookup or workflow is started.
Requires Selenium, Chrome and a compatible chromedriver.
"""
import json
import os

import pytest

from mzdesk import __version__
from mzdesk.spec import spec


URL = os.environ.get("MZDESK_UI_TEST_URL")
pytestmark = pytest.mark.skipif(not URL, reason="set MZDESK_UI_TEST_URL to a running frontend")


@pytest.fixture
def browser():
    pytest.importorskip("selenium")
    from selenium import webdriver
    from selenium.webdriver.chrome.options import Options
    from selenium.webdriver.chrome.service import Service

    options = Options()
    options.add_argument("--headless=new")
    options.add_argument("--window-size=1440,1000")
    driver_path = os.environ.get("MZDESK_CHROMEDRIVER")
    driver = webdriver.Chrome(service=Service(executable_path=driver_path) if driver_path else Service(), options=options)
    data = {**spec(), "version": __version__}
    script = """
    const spec = SPEC_DATA;
    const originalFetch = window.fetch.bind(window);
    const response = data => new Response(JSON.stringify(data), {status: 200, headers: {'Content-Type': 'application/json'}});
    window.__delayModule = false;
    window.__releaseModule = null;
    window.__moduleCalls = 0;
    window.__delayWorkflow = false;
    window.__releaseWorkflow = null;
    window.__workflowPolls = 0;
    window.fetch = async (input, init) => {
      const url = new URL(typeof input === 'string' ? input : input.url, location.href);
      if (url.pathname.endsWith('/api/spec')) return response(spec);
      if (url.pathname.endsWith('/api/track')) return response({});
      if (url.pathname.endsWith('/api/modules/structure_formats/run')) {
        window.__moduleCalls++;
        const delayed = window.__delayModule;
        if (delayed) await new Promise(resolve => { window.__releaseModule = resolve; });
        window.__oldModuleDone = delayed;
        return response({status:'ok', verdict:delayed ? 'Old structure result' : 'Structure result retained', warnings:[], outputs:{
          compare:{compounds:[{label:'Test compound', formula:'C2H6'}]},
          formats:{formats:{SMILES:'CC'}}
        }});
      }
      if (url.pathname.endsWith('/api/pipelines/esi_ms2_annotation/run')) {
        if (window.__delayWorkflow) await new Promise(resolve => { window.__releaseWorkflow = resolve; });
        return response({job_id:'test-job', status:'running'});
      }
      if (url.pathname.endsWith('/api/pipelines/jobs/test-job')) {
        window.__workflowPolls++;
        return response({status:'done', note:'', steps:[], error:null, result:{
          answer:{type:'none', ids:[], statement:'Workflow result retained', top_id:null, isomer_verdict:'not run'},
          candidates:[], warnings:[], context:{precursor_mz:195.0877, polarity:'Positive', libraries:[], gate_passed:false},
          isomer:null
        }});
      }
      return originalFetch(input, init);
    };
    """.replace("SPEC_DATA", json.dumps(data))
    driver.execute_cdp_cmd("Page.addScriptToEvaluateOnNewDocument", {"source": script})
    yield driver
    driver.quit()


def button(driver, text):
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support import expected_conditions as EC
    from selenium.webdriver.support.ui import WebDriverWait

    return WebDriverWait(driver, 15).until(EC.element_to_be_clickable((By.XPATH, f"//main//button[not(ancestor::*[@hidden]) and normalize-space(.)='{text}']")))


def visible_text(driver, text):
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support import expected_conditions as EC
    from selenium.webdriver.support.ui import WebDriverWait

    return WebDriverWait(driver, 15).until(EC.visibility_of_element_located((By.XPATH, f"//main//*[not(ancestor::*[@hidden]) and contains(text(), '{text}')]")))


def test_results_survive_navigation_and_clear_only_current_page(browser):
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.ui import WebDriverWait

    browser.get(URL.rstrip("/") + "/modules/structure_formats")
    button(browser, "Load example").click()
    button(browser, "Run Structure & formats").click()
    visible_text(browser, "Structure result retained")
    browser.find_element(By.CSS_SELECTOR, 'a[href$="/pipelines/esi_ms2_annotation"]').click()
    button(browser, "Load example: Caffeine, library match").click()
    button(browser, "Run workflow").click()
    visible_text(browser, "Workflow result retained")

    browser.back()
    visible_text(browser, "Structure result retained")
    browser.forward()
    visible_text(browser, "Workflow result retained")
    button(browser, "Clear").click()
    WebDriverWait(browser, 5).until(lambda d: not d.find_elements(By.XPATH, "//main//*[not(ancestor::*[@hidden]) and contains(text(), 'Workflow result retained')]"))

    browser.back()
    visible_text(browser, "Structure result retained")
    button(browser, "Clear").click()
    WebDriverWait(browser, 5).until(lambda d: not d.find_elements(By.XPATH, "//main//*[not(ancestor::*[@hidden]) and contains(text(), 'Structure result retained')]"))
    field = browser.find_element(By.CSS_SELECTOR, '[data-module="structure_formats"] [data-field="structure"] input[type="text"]')
    assert field.get_attribute("value") == ""


def test_run_finishes_while_another_page_is_open(browser):
    browser.get(URL.rstrip("/") + "/modules/structure_formats")
    browser.execute_script("window.__delayModule = true")
    button(browser, "Load example").click()
    button(browser, "Run Structure & formats").click()
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.ui import WebDriverWait
    WebDriverWait(browser, 5).until(lambda d: d.execute_script("return !!window.__releaseModule"))
    browser.find_element(By.CSS_SELECTOR, 'a[href$="/pipelines/esi_ms2_annotation"]').click()
    button(browser, "Run workflow")
    browser.execute_script("window.__releaseModule()")
    browser.back()
    visible_text(browser, "Old structure result")


def test_clear_ignores_a_late_result(browser):
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.ui import WebDriverWait

    browser.get(URL.rstrip("/") + "/modules/structure_formats")
    browser.execute_script("window.__delayModule = true")
    button(browser, "Load example").click()
    button(browser, "Run Structure & formats").click()
    WebDriverWait(browser, 5).until(lambda d: d.execute_script("return !!window.__releaseModule"))
    button(browser, "Clear").click()
    visible_text(browser, "previous run may still finish")
    browser.execute_script("window.__delayModule = false; window.__releaseModule()")
    WebDriverWait(browser, 5).until(lambda d: d.execute_script("return window.__oldModuleDone"))
    assert not browser.find_elements(By.XPATH, "//main//*[not(ancestor::*[@hidden]) and contains(text(), 'Old structure result')]")
    # Run a fresh request: if the old one can overwrite it, the page's state is not isolated.
    button(browser, "Load example").click()
    button(browser, "Run Structure & formats").click()
    visible_text(browser, "Structure result retained")
    assert len(browser.find_elements(By.CSS_SELECTOR, '[data-module="structure_formats"] .verdict')) == 1


def test_clear_resets_uploaded_file(browser, tmp_path):
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.ui import WebDriverWait

    sample = tmp_path / "sample.mgf"
    sample.write_text("BEGIN IONS\n100 50\nEND IONS\n")
    browser.get(URL.rstrip("/") + "/modules/gnps_massql")
    file_input = WebDriverWait(browser, 15).until(lambda d: d.find_element(By.CSS_SELECTOR, '[data-module="gnps_massql"] input[type="file"]'))
    file_input.send_keys(str(sample))
    visible_text(browser, "sample.mgf")
    button(browser, "Clear").click()
    assert browser.find_element(By.CSS_SELECTOR, '[data-module="gnps_massql"] input[type="file"]').get_attribute("value") == ""
    assert not browser.find_elements(By.XPATH, "//main//*[not(ancestor::*[@hidden]) and text()='sample.mgf']")


def test_radio_selections_are_independent_between_pages(browser):
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.ui import WebDriverWait

    browser.get(URL.rstrip("/") + "/modules/structure_formats")
    first = WebDriverWait(browser, 15).until(lambda d: d.find_element(By.CSS_SELECTOR, '[data-module="structure_formats"] [data-field="structure"]'))
    first.find_element(By.XPATH, ".//label[contains(., 'SMILES / InChI')]").click()
    first.find_element(By.CSS_SELECTOR, 'input[type="text"]').send_keys("CC")
    browser.find_element(By.CSS_SELECTOR, 'a[href$="/modules/properties"]').click()
    second = WebDriverWait(browser, 15).until(lambda d: d.find_element(By.CSS_SELECTOR, '[data-module="properties"] [data-field="structure"]'))
    second.find_element(By.XPATH, ".//label[contains(., 'SMILES / InChI')]").click()
    browser.back()
    WebDriverWait(browser, 5).until(lambda d: d.find_element(By.CSS_SELECTOR, '[data-module="structure_formats"]').is_displayed())
    assert first.find_element(By.XPATH, ".//label[contains(., 'SMILES / InChI')]/input").is_selected()
    assert first.find_element(By.CSS_SELECTOR, 'input[type="text"]').get_attribute("value") == "CC"


def test_clear_workflow_ignores_a_late_start(browser):
    from selenium.webdriver.support.ui import WebDriverWait

    browser.get(URL.rstrip("/") + "/pipelines/esi_ms2_annotation")
    browser.execute_script("window.__delayWorkflow = true")
    button(browser, "Load example: Caffeine, library match").click()
    button(browser, "Run workflow").click()
    WebDriverWait(browser, 5).until(lambda d: d.execute_script("return !!window.__releaseWorkflow"))
    button(browser, "Clear").click()
    visible_text(browser, "previous workflow may still finish")
    browser.execute_script("window.__delayWorkflow = false; window.__releaseWorkflow()")
    # A new run is allowed, and only that run should be polled and displayed.
    button(browser, "Load example: Caffeine, library match").click()
    button(browser, "Run workflow").click()
    visible_text(browser, "Workflow result retained")
    assert browser.execute_script("return window.__workflowPolls") == 1


def test_switching_pages_dismisses_pending_network_confirmation(browser):
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.ui import WebDriverWait

    browser.get(URL.rstrip("/") + "/modules/structure_formats")
    browser.find_element(By.CSS_SELECTOR, 'a[href$="/pipelines/esi_ms2_annotation"]').click()
    button(browser, "Run workflow")
    browser.back()
    button(browser, "Load example").click()
    browser.find_element(By.CSS_SELECTOR, 'button[aria-label="Internet search"]').click()
    button(browser, "Run Structure & formats").click()
    WebDriverWait(browser, 5).until(lambda d: d.find_elements(By.CSS_SELECTOR, '[role="dialog"][aria-modal="true"]'))
    browser.forward()
    button(browser, "Run workflow")
    WebDriverWait(browser, 5).until(lambda d: not d.find_elements(By.CSS_SELECTOR, '[role="dialog"][aria-modal="true"]'))
    assert browser.execute_script("return window.__moduleCalls") == 0


def test_each_requested_page_has_its_own_header_icon(browser):
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.ui import WebDriverWait

    wanted = {"molecular", "ms1", "ms2", "gnps", "tools"}
    modules = [m for m in spec()["modules"] if m["subgroup"] in wanted]
    assert len(modules) == 21
    for module in modules:
        browser.get(URL.rstrip("/") + "/modules/" + module["id"])
        page = WebDriverWait(browser, 15).until(lambda d: d.find_element(By.CSS_SELECTOR, f'[data-module="{module["id"]}"]'))
        assert page.find_element(By.TAG_NAME, "h1").text == module["title"]
        assert page.find_element(By.CSS_SELECTOR, ".module-kicker").text
        assert page.find_element(By.CSS_SELECTOR, ".module-emblem use").get_attribute("href").endswith("#" + module["id"])

    browser.get(URL.rstrip("/") + "/pipelines/esi_ms2_annotation")
    page = WebDriverWait(browser, 15).until(lambda d: d.find_element(By.CSS_SELECTOR, '[data-pipeline="esi_ms2_annotation"]'))
    assert page.find_element(By.CSS_SELECTOR, ".module-kicker").text == "WORKFLOWS"
    assert page.find_element(By.CSS_SELECTOR, ".module-emblem use").get_attribute("href").endswith("#esi_ms2_annotation")


@pytest.mark.parametrize("path", ["modules/structure_formats", "modules/ms2_similarity", "modules/gnps_massql", "pipelines/esi_ms2_annotation"])
def test_pages_fit_a_narrow_viewport(browser, path):
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.ui import WebDriverWait

    browser.set_window_size(390, 844)
    browser.get(URL.rstrip("/") + "/" + path)
    WebDriverWait(browser, 15).until(lambda d: d.find_element(By.CSS_SELECTOR, "main .module-heading"))
    assert browser.execute_script("return document.documentElement.scrollWidth <= window.innerWidth")
