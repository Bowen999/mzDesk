"""Opt-in browser regression: MZDESK_UI_TEST_URL=http://localhost:5173 pytest tests/test_ui_ms2_hover.py.

The API is mocked, but Plotly and the shared React figure run in a real browser.
"""
import json
import os

import pytest

from mzdesk.api.common import plotly
from mzdesk.ms2.mirror_plot import build_mirror_figure
from mzdesk.shared.ms2_plot import spectrum_figure
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
    single, cards = spectrum_figure([(100, 100), (150, 1), (200, 0)], dim={1: "noise"})
    cards[0]["color"] = "#00aabc"
    cards[0]["img"] = "data:image/svg+xml," + "<svg xmlns='http://www.w3.org/2000/svg' width='220' height='190'/>"
    mirror, mirror_cards = build_mirror_figure([(100, 100), (150, 10)], [(100, 100), (160, 30)])
    fragment_table = {
        "columns": ["Candidate", "Fragments"],
        "rows": [["Test", {"fragment": {"smiles": "CC", "adduct": "[M+H]+", "name": "Test"}}]],
        "fragments": {"peaks": [[100, 100]], "precursor_mz": 200, "tol_da": 0.01},
    }
    figures = {"ms2_clean": {"spectrum": plotly(single, peaks=cards)},
               "ms2_similarity": {"mirror": plotly(mirror, peaks=mirror_cards), "scores": fragment_table}}
    svg = '<svg xmlns="http://www.w3.org/2000/svg" width="100" height="100"></svg>'
    view = {
        "name": "Test", "smiles": "CC", "formula": "C2H6", "adduct": "[M+H]+", "ion_formula": "C2H7",
        "precursor_mz": 200, "theo_mz": 200, "positive": True, "peaks": 3, "explained": 1, "explained_fraction": 0.9,
        "coverage": {"atoms": 1, "heavy": 2, "fraction": 0.5}, "overview": svg,
        "fragments": [{"color": "#00aabc", "mz": 100, "rel": 100, "ion_formula": "C2H7", "ppm": 0, "placed": True,
                       "alternatives": 0, "svg": svg}],
        "figure": plotly(single, peaks=cards), "notes": [],
        "settings": {"tol_da": 0.01, "tol_ppm": 0, "tree_depth": 2, "engine": "MetFrag"}, "seconds": 1,
    }
    script = """
    const spec = SPEC_DATA, figures = FIGURE_DATA, view = VIEW_DATA;
    const originalFetch = window.fetch.bind(window);
    const response = data => new Response(JSON.stringify(data), {status: 200, headers: {'Content-Type': 'application/json'}});
    window.fetch = async (input, init) => {
      const path = new URL(typeof input === 'string' ? input : input.url, location.href).pathname;
      if (path.endsWith('/api/spec')) return response(spec);
      if (path.endsWith('/api/track')) return response({});
      if (path.endsWith('/api/fragments/run')) return response({job_id:'test-fragments', status:'queued'});
      if (path.endsWith('/api/fragments/jobs/test-fragments')) return response({status:'done', result:view});
      for (const [name, outputs] of Object.entries(figures)) {
        if (path.endsWith(`/api/modules/${name}/run`)) return response({status:'ok', warnings:[], outputs});
      }
      return originalFetch(input, init);
    };
    """.replace("SPEC_DATA", json.dumps(spec())).replace("FIGURE_DATA", json.dumps(figures)).replace("VIEW_DATA", json.dumps(view))
    driver.execute_cdp_cmd("Page.addScriptToEvaluateOnNewDocument", {"source": script})
    yield driver
    driver.quit()


def open_figure(browser, module, title):
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.ui import WebDriverWait

    wait = WebDriverWait(browser, 20)
    browser.get(URL.rstrip("/") + "/modules/" + module)
    wait.until(lambda d: d.find_element(By.XPATH, "//main//button[normalize-space(.)='Load example']")).click()
    wait.until(lambda d: d.find_element(By.XPATH, f"//main//button[normalize-space(.)='Run {title}']")).click()
    gd = wait.until(lambda d: d.find_element(By.CSS_SELECTOR, ".js-plotly-plot"))
    wait.until(lambda d: d.execute_script("return !!arguments[0]._fullLayout?.xaxis", gd))
    wait.until(lambda d: not d.find_elements(By.CSS_SELECTOR, ".fig-loading"))
    center_plot(browser, gd)
    return gd


def center_plot(browser, gd):
    # scroll events are delivered on a later frame; let them finish before moving the pointer.
    browser.execute_async_script("""
      const done = arguments[arguments.length - 1];
      arguments[0].scrollIntoView({block:'center', behavior:'instant'});
      requestAnimationFrame(() => requestAnimationFrame(done));
    """, gd)


def move_to_peak(browser, gd, mz, y):
    from selenium.webdriver.common.action_chains import ActionChains

    x, y = browser.execute_script("""
      const gd = arguments[0], rect = gd.getBoundingClientRect(), {xaxis, yaxis} = gd._fullLayout;
      return [Math.round(rect.left + xaxis._offset + xaxis.l2p(arguments[1])),
              Math.round(rect.top + yaxis._offset + yaxis.l2p(arguments[2]))];
    """, gd, mz, y)
    actions = ActionChains(browser)
    actions.w3c_actions.pointer_action.move_to_location(x, y)
    actions.perform()


def tooltip(browser):
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.ui import WebDriverWait

    return WebDriverWait(browser, 5).until(lambda d: d.find_element(By.CSS_SELECTOR, ".hover-card"))


def test_single_stick_highlight_details_and_clipping(browser):
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.ui import WebDriverWait

    gd = open_figure(browser, "ms2_clean", "Spectrum Cleanup")
    # A transformed, clipping ancestor must not trap the card (as in the fragment dialog).
    browser.execute_script("arguments[0].parentElement.style.transform='translateZ(0)'", gd)
    move_to_peak(browser, gd, 100, 50)  # middle of the stick, not its tip
    card = tooltip(browser)
    assert "m/z 100" in card.text and "Intensity" in card.text
    assert browser.execute_script("return arguments[0].parentElement === document.body", card)
    assert browser.find_element(By.CSS_SELECTOR, ".fig-peak-line").get_attribute("stroke") == "#00aabc"
    assert browser.find_element(By.CSS_SELECTOR, ".fig-peak-highlight circle")
    assert browser.execute_script("""
      const r = arguments[0].getBoundingClientRect();
      return r.left >= 0 && r.top >= 0 && r.right <= innerWidth && r.bottom <= innerHeight;
    """, card)
    move_to_peak(browser, gd, 150, 0.5)
    assert "Removed: noise" in tooltip(browser).text
    move_to_peak(browser, gd, 170, 50)
    WebDriverWait(browser, 5).until(lambda d: not d.find_elements(By.CSS_SELECTOR, ".hover-card, .fig-peak-highlight"))


def test_mirror_half_switch_zoom_and_dark_theme(browser):
    from selenium.webdriver.common.by import By
    from selenium.webdriver.common.action_chains import ActionChains
    from selenium.webdriver.support.ui import WebDriverWait

    gd = open_figure(browser, "ms2_similarity", "Similarity")
    move_to_peak(browser, gd, 100, 50)
    assert "Query · matched" in tooltip(browser).text
    move_to_peak(browser, gd, 100, -50)
    assert "Reference · matched" in tooltip(browser).text
    line = browser.find_element(By.CSS_SELECTOR, ".fig-peak-line")
    assert float(line.get_attribute("y2")) > float(line.get_attribute("y1"))

    # Real Plotly drag-to-zoom clears the hover and keeps hit testing aligned with the new axes.
    start, end = browser.execute_script("""
      const gd = arguments[0], r = gd.getBoundingClientRect(), {xaxis:x, yaxis:y} = gd._fullLayout;
      const pt = (mz, rel) => [Math.round(r.left + x._offset + x.l2p(mz)), Math.round(r.top + y._offset + y.l2p(rel))];
      return [pt(95, 110), pt(130, -110)];
    """, gd)
    actions = ActionChains(browser)
    actions.w3c_actions.pointer_action.move_to_location(*start)
    actions.w3c_actions.pointer_action.pointer_down()
    actions.w3c_actions.pointer_action.move_to_location(*end)
    actions.w3c_actions.pointer_action.pointer_up()
    actions.perform()
    WebDriverWait(browser, 5).until(lambda d: d.execute_script("return arguments[0]._fullLayout.xaxis.range[1] < 150", gd))
    move_to_peak(browser, gd, 100, -50)
    assert "Reference · matched" in tooltip(browser).text

    browser.execute_script("localStorage.setItem('mzdesk.theme', 'dark')")
    gd = open_figure(browser, "ms2_similarity", "Similarity")
    move_to_peak(browser, gd, 100, 50)
    assert "Query · matched" in tooltip(browser).text
    assert browser.execute_script("return document.documentElement.getAttribute('data-theme')") == "dark"
    assert browser.find_element(By.CSS_SELECTOR, ".fig-peak-line").get_attribute("stroke") == "var(--accent)"


def test_fragment_dialog_shows_colored_peak_and_card_above_overlay(browser):
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.ui import WebDriverWait

    open_figure(browser, "ms2_similarity", "Similarity")
    wait = WebDriverWait(browser, 20)
    wait.until(lambda d: d.find_element(By.CSS_SELECTOR, ".fv-open")).click()
    gd = wait.until(lambda d: d.find_element(By.CSS_SELECTOR, ".fv-card .js-plotly-plot"))
    wait.until(lambda d: d.execute_script("return !!arguments[0]._fullLayout?.xaxis", gd))
    wait.until(lambda d: not d.find_elements(By.CSS_SELECTOR, ".fv-card .fig-loading"))
    center_plot(browser, gd)
    move_to_peak(browser, gd, 100, 50)
    card = tooltip(browser)
    assert "m/z 100" in card.text and "Intensity" in card.text
    assert browser.find_element(By.CSS_SELECTOR, ".fv-card .fig-peak-line").get_attribute("stroke") == "#00aabc"
    assert browser.execute_script("""
      const r = arguments[0].getBoundingClientRect(), overlay = document.querySelector('.fv-overlay');
      return arguments[0].parentElement === document.body &&
        Number(getComputedStyle(arguments[0]).zIndex) > Number(getComputedStyle(overlay).zIndex) &&
        r.left >= 0 && r.top >= 0 && r.right <= innerWidth && r.bottom <= innerHeight;
    """, card)
    browser.find_element(By.CSS_SELECTOR, ".fv-close").click()
    wait.until(lambda d: not d.find_elements(By.CSS_SELECTOR, ".hover-card, .fv-card"))


def test_tooltip_fits_short_narrow_viewport_and_clears_on_scroll(browser):
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.ui import WebDriverWait

    browser.set_window_size(400, 500)
    gd = open_figure(browser, "ms2_clean", "Spectrum Cleanup")
    move_to_peak(browser, gd, 100, 50)
    card = tooltip(browser)
    assert browser.execute_script("""
      const r = arguments[0].getBoundingClientRect();
      return r.left >= 0 && r.top >= 0 && r.right <= innerWidth && r.bottom <= innerHeight &&
        arguments[0].scrollHeight <= arguments[0].clientHeight;
    """, card)
    browser.execute_script("window.scrollBy(0, 40)")
    WebDriverWait(browser, 5).until(lambda d: not d.find_elements(By.CSS_SELECTOR, ".hover-card, .fig-peak-highlight"))
